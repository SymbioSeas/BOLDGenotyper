"""
Unit tests for genotype assignment module.

Tests cover:
- Edit distance calculation (edlib and pure Python)
- ProcessID extraction from FASTA headers
- Best match finding with ties and low-confidence assignments
- Parallel processing (single and multi-threaded)
- Missing sample handling
- Diagnostics generation
- Full integration workflow
"""

import unittest
from unittest.mock import patch, Mock
import tempfile
import shutil
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from boldgenotyper import genotype_assignment


class TestEditDistance(unittest.TestCase):
    """Test edit distance calculations."""

    def test_levenshtein_identical_sequences(self):
        """Test Levenshtein distance for identical sequences."""
        dist = genotype_assignment.levenshtein_distance("ACGT", "ACGT")
        self.assertEqual(dist, 0)

    def test_levenshtein_single_substitution(self):
        """Test Levenshtein distance with one substitution."""
        dist = genotype_assignment.levenshtein_distance("ACGT", "AGGT")
        self.assertEqual(dist, 1)

    def test_levenshtein_multiple_changes(self):
        """Test Levenshtein distance with multiple changes."""
        dist = genotype_assignment.levenshtein_distance("ACGT", "TGCA")
        self.assertEqual(dist, 4)

    def test_levenshtein_insertion(self):
        """Test Levenshtein distance with insertion."""
        dist = genotype_assignment.levenshtein_distance("ACGT", "ACGGT")
        self.assertEqual(dist, 1)

    def test_levenshtein_deletion(self):
        """Test Levenshtein distance with deletion."""
        dist = genotype_assignment.levenshtein_distance("ACGGT", "ACGT")
        self.assertEqual(dist, 1)

    def test_levenshtein_empty_strings(self):
        """Test Levenshtein distance with empty strings."""
        dist = genotype_assignment.levenshtein_distance("", "ACGT")
        self.assertEqual(dist, 4)
        dist = genotype_assignment.levenshtein_distance("ACGT", "")
        self.assertEqual(dist, 4)

    def test_levenshtein_symmetry(self):
        """Test that Levenshtein distance is symmetric."""
        dist1 = genotype_assignment.levenshtein_distance("ACGT", "TGCA")
        dist2 = genotype_assignment.levenshtein_distance("TGCA", "ACGT")
        self.assertEqual(dist1, dist2)

    def test_calculate_edit_distance_pure_python(self):
        """Test edit distance with pure Python (no edlib)."""
        dist = genotype_assignment.calculate_edit_distance(
            "ACGT", "AGGT", use_edlib=False
        )
        self.assertEqual(dist, 1)

    def test_calculate_edit_distance_with_edlib(self):
        """Test edit distance with edlib if available."""
        if genotype_assignment.EDLIB_AVAILABLE:
            dist = genotype_assignment.calculate_edit_distance(
                "ACGT", "AGGT", use_edlib=True
            )
            self.assertEqual(dist, 1)
        else:
            self.skipTest("edlib not available")

    def test_edlib_and_python_agree(self):
        """Test that edlib and Python implementations agree."""
        if not genotype_assignment.EDLIB_AVAILABLE:
            self.skipTest("edlib not available")

        test_cases = [
            ("ACGT", "ACGT"),
            ("ACGT", "AGGT"),
            ("ACGTACGT", "ACGTACGA"),
            ("AAAA", "TTTT"),
        ]

        for seq1, seq2 in test_cases:
            dist_edlib = genotype_assignment.calculate_edit_distance(
                seq1, seq2, use_edlib=True
            )
            dist_python = genotype_assignment.calculate_edit_distance(
                seq1, seq2, use_edlib=False
            )
            self.assertEqual(
                dist_edlib, dist_python,
                f"edlib and Python disagree for {seq1} vs {seq2}"
            )


class TestIdentityCalculation(unittest.TestCase):
    """Test sequence identity calculation."""

    def test_calculate_identity_perfect(self):
        """Test identity for perfect match."""
        identity = genotype_assignment.calculate_identity(0, 100, 100)
        self.assertEqual(identity, 1.0)

    def test_calculate_identity_ninety_percent(self):
        """Test identity for 90% match."""
        identity = genotype_assignment.calculate_identity(10, 100, 100)
        self.assertEqual(identity, 0.9)

    def test_calculate_identity_different_lengths(self):
        """Test identity with different sequence lengths."""
        # Distance 5, lengths 100 and 95
        # Identity = 1 - (5 / 100) = 0.95
        identity = genotype_assignment.calculate_identity(5, 100, 95)
        self.assertEqual(identity, 0.95)

    def test_calculate_identity_zero_length(self):
        """Test identity with zero-length sequences."""
        identity = genotype_assignment.calculate_identity(0, 0, 0)
        self.assertEqual(identity, 1.0)

    def test_calculate_identity_complete_mismatch(self):
        """Test identity for complete mismatch."""
        identity = genotype_assignment.calculate_identity(100, 100, 100)
        self.assertEqual(identity, 0.0)


class TestProcessIDExtraction(unittest.TestCase):
    """Test processid extraction from FASTA headers."""

    def test_extract_processid_standard_format(self):
        """Test processid extraction from standard BOLD format."""
        header = "Sphyrna_lewini_ANGBF11456-15.COI-5P"
        pid = genotype_assignment.extract_processid_from_header(header)
        self.assertEqual(pid, "ANGBF11456-15")

    def test_extract_processid_with_leading_angle_bracket(self):
        """Test processid extraction with leading >."""
        header = ">Genus_species_BOLD12345.marker"
        pid = genotype_assignment.extract_processid_from_header(header)
        self.assertEqual(pid, "BOLD12345")

    def test_extract_processid_simple_format(self):
        """Test processid extraction from simple format."""
        header = "species_SAMPLE001"
        pid = genotype_assignment.extract_processid_from_header(header)
        self.assertEqual(pid, "SAMPLE001")

    def test_extract_processid_with_spaces(self):
        """Test processid extraction with spaces after processid."""
        header = "genus_species_ABC123 some description"
        pid = genotype_assignment.extract_processid_from_header(header)
        self.assertEqual(pid, "ABC123")

    def test_extract_processid_no_match(self):
        """Test processid extraction when pattern doesn't match."""
        # Header with no underscores - regex won't match
        header = "nounderscore"
        pid = genotype_assignment.extract_processid_from_header(header)
        self.assertIsNone(pid)

    def test_extract_processid_empty_string(self):
        """Test processid extraction from empty string."""
        pid = genotype_assignment.extract_processid_from_header("")
        self.assertIsNone(pid)


class TestBestMatchFinding(unittest.TestCase):
    """Test best consensus match finding."""

    def test_find_best_match_perfect(self):
        """Test finding best match with perfect identity."""
        consensus = [
            ("c1_n10", "ACGTACGT"),
            ("c2_n5", "TGCATGCA"),
        ]
        result = genotype_assignment.find_best_consensus_match(
            "ACGTACGT", consensus, min_identity=0.90
        )

        self.assertEqual(result['best_group'], "c1_n10")
        self.assertEqual(result['best_identity'], 1.0)
        self.assertEqual(result['runner_up_group'], "c2_n5")
        self.assertFalse(result['is_tie'])

    def test_find_best_match_below_threshold(self):
        """Test that sequences below threshold are unassigned."""
        consensus = [
            ("c1_n10", "AAAAAAAAAA"),
            ("c2_n5", "TTTTTTTTTT"),
        ]
        result = genotype_assignment.find_best_consensus_match(
            "CCCCCCCCCC", consensus, min_identity=0.90
        )

        self.assertIsNone(result['best_group'])
        self.assertLess(result['best_identity'], 0.90)

    def test_find_best_match_with_tie(self):
        """Test detection of ties between nearly equal matches."""
        consensus = [
            ("c1_n10", "ACGTACGTACGTACGT"),
            ("c2_n5", "ACGTACGTACGTACGA"),  # 1 bp different
        ]
        # Query is 1 bp different from each (equidistant)
        query = "ACGTACGTACGTACGG"

        result = genotype_assignment.find_best_consensus_match(
            query, consensus, min_identity=0.80
        )

        # Should detect tie
        self.assertTrue(result['is_tie'])
        self.assertLess(
            abs(result['best_identity'] - result['runner_up_identity']),
            0.01
        )

    def test_find_best_match_low_confidence(self):
        """Test detection of low-confidence assignments."""
        consensus = [
            ("c1_n10", "ACGTACGTACGTACGT"),  # 16 bp
        ]
        # Create query with ~93.75% identity (just above 90% threshold)
        # 1 bp different -> 93.75% identity, which is in low-confidence range (90-95%)
        query = "ACGTACGTACGTACGA"  # Only last bp different

        result = genotype_assignment.find_best_consensus_match(
            query, consensus, min_identity=0.90
        )

        self.assertEqual(result['best_group'], "c1_n10")
        self.assertTrue(result['is_low_confidence'])

    def test_find_best_match_multiple_consensus(self):
        """Test with multiple consensus sequences."""
        consensus = [
            ("c1_n10", "ACGTACGT"),
            ("c2_n8", "TGCATGCA"),
            ("c3_n5", "GGCCGGCC"),
            ("c4_n3", "AATTAATT"),
        ]

        result = genotype_assignment.find_best_consensus_match(
            "GGCCGGCC", consensus, min_identity=0.90
        )

        self.assertEqual(result['best_group'], "c3_n5")
        self.assertEqual(result['best_identity'], 1.0)


class TestFASTAReading(unittest.TestCase):
    """Test FASTA file reading."""

    def setUp(self):
        """Create temporary directory for test files."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_read_fasta_simple(self):
        """Test reading simple FASTA file."""
        fasta_path = self.temp_dir / "test.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq1\nACGT\n>seq2\nTGCA\n")

        records = genotype_assignment.read_fasta_simple(str(fasta_path))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0], ("seq1", "ACGT"))
        self.assertEqual(records[1], ("seq2", "TGCA"))

    def test_read_fasta_multiline_sequences(self):
        """Test reading FASTA with multiline sequences."""
        fasta_path = self.temp_dir / "test.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq1\nACGT\nACGT\n>seq2\nTGCA\nTGCA\n")

        records = genotype_assignment.read_fasta_simple(str(fasta_path))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0], ("seq1", "ACGTACGT"))
        self.assertEqual(records[1], ("seq2", "TGCATGCA"))

    def test_read_fasta_with_empty_lines(self):
        """Test reading FASTA with empty lines."""
        fasta_path = self.temp_dir / "test.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq1\n\nACGT\n\n>seq2\nTGCA\n")

        records = genotype_assignment.read_fasta_simple(str(fasta_path))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0], ("seq1", "ACGT"))

    def test_read_fasta_uppercase_conversion(self):
        """Test that sequences are converted to uppercase."""
        fasta_path = self.temp_dir / "test.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq1\nacgt\n")

        records = genotype_assignment.read_fasta_simple(str(fasta_path))
        self.assertEqual(records[0][1], "ACGT")

    def test_read_fasta_nonexistent_file(self):
        """Test that reading nonexistent file raises error."""
        with self.assertRaises(FileNotFoundError):
            genotype_assignment.read_fasta_simple("nonexistent.fasta")

    def test_read_fasta_empty_file(self):
        """Test that empty FASTA file raises error."""
        fasta_path = self.temp_dir / "empty.fasta"
        fasta_path.touch()

        with self.assertRaises(genotype_assignment.GenotypeAssignmentError):
            genotype_assignment.read_fasta_simple(str(fasta_path))


class TestParallelProcessing(unittest.TestCase):
    """Test parallel genotype assignment."""

    def setUp(self):
        """Create temporary directory and test data."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self._create_test_data()

    def tearDown(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def _create_test_data(self):
        """Create test files for parallel processing tests."""
        # Create metadata
        metadata = pd.DataFrame({
            'processid': [f'SAMPLE{i:03d}' for i in range(1, 11)],
            'species': ['Species_A'] * 10
        })
        self.metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(self.metadata_path, sep='\t', index=False)

        # Create raw sequences
        raw_records = []
        for i in range(1, 11):
            # Create slight variations of ACGT repeated pattern
            seq = "ACGTACGT" * 5
            if i % 2 == 0:
                seq = seq[:-1] + "A"  # Slight variation
            record = SeqRecord(
                Seq(seq),
                id=f"Species_A_SAMPLE{i:03d}.COI",
                description=""
            )
            raw_records.append(record)

        self.fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, self.fasta_path, 'fasta')

        # Create consensus sequences
        consensus_records = [
            SeqRecord(Seq("ACGTACGT" * 5), id="consensus_c1", description=""),
        ]
        self.consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, self.consensus_path, 'fasta')

    def test_single_process_assignment(self):
        """Test genotype assignment with single process."""
        output_path = self.temp_dir / 'output_single.tsv'

        results = genotype_assignment.assign_genotypes(
            metadata_path=str(self.metadata_path),
            fasta_path=str(self.fasta_path),
            consensus_path=str(self.consensus_path),
            output_path=str(output_path),
            n_processes=1
        )

        self.assertEqual(results['total'], 10)
        self.assertEqual(results['assigned'], 10)
        self.assertTrue(output_path.exists())

    def test_multi_process_assignment(self):
        """Test genotype assignment with multiple processes."""
        output_path = self.temp_dir / 'output_multi.tsv'

        results = genotype_assignment.assign_genotypes(
            metadata_path=str(self.metadata_path),
            fasta_path=str(self.fasta_path),
            consensus_path=str(self.consensus_path),
            output_path=str(output_path),
            n_processes=2
        )

        self.assertEqual(results['total'], 10)
        self.assertEqual(results['assigned'], 10)
        self.assertTrue(output_path.exists())

    def test_parallel_results_consistency(self):
        """Test that parallel and serial produce same results."""
        output_single = self.temp_dir / 'output_single.tsv'
        output_multi = self.temp_dir / 'output_multi.tsv'

        # Run with single process
        genotype_assignment.assign_genotypes(
            metadata_path=str(self.metadata_path),
            fasta_path=str(self.fasta_path),
            consensus_path=str(self.consensus_path),
            output_path=str(output_single),
            n_processes=1
        )

        # Run with multiple processes
        genotype_assignment.assign_genotypes(
            metadata_path=str(self.metadata_path),
            fasta_path=str(self.fasta_path),
            consensus_path=str(self.consensus_path),
            output_path=str(output_multi),
            n_processes=2
        )

        # Compare results
        df_single = pd.read_csv(output_single, sep='\t')
        df_multi = pd.read_csv(output_multi, sep='\t')

        pd.testing.assert_frame_equal(df_single, df_multi)


class TestMissingSamples(unittest.TestCase):
    """Test handling of missing samples."""

    def setUp(self):
        """Create temporary directory and test data."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_missing_sequences_in_fasta(self):
        """Test handling when some processids missing from FASTA."""
        # Metadata has 5 samples
        metadata = pd.DataFrame({
            'processid': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003', 'SAMPLE004', 'SAMPLE005'],
            'species': ['Species_A'] * 5
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        # FASTA only has 3 samples
        raw_records = [
            SeqRecord(Seq('ACGTACGT'), id='Species_A_SAMPLE001.COI', description=''),
            SeqRecord(Seq('ACGTACGT'), id='Species_A_SAMPLE002.COI', description=''),
            SeqRecord(Seq('ACGTACGT'), id='Species_A_SAMPLE003.COI', description=''),
        ]
        fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, fasta_path, 'fasta')

        # Consensus
        consensus_records = [
            SeqRecord(Seq('ACGTACGT'), id='consensus_c1', description=''),
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        # Run assignment
        output_path = self.temp_dir / 'output.tsv'
        results = genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path)
        )

        # Check results
        self.assertEqual(results['total'], 5)
        self.assertEqual(results['assigned'], 3)
        self.assertEqual(results['no_sequence'], 2)

        # Check output file
        df = pd.read_csv(output_path, sep='\t')
        self.assertTrue(df.loc[df['processid'] == 'SAMPLE004', 'consensus_group'].isna().all())
        self.assertTrue(df.loc[df['processid'] == 'SAMPLE005', 'consensus_group'].isna().all())

    def test_unparseable_fasta_headers(self):
        """Test handling when FASTA headers can't be parsed for processid."""
        metadata = pd.DataFrame({
            'processid': ['SAMPLE001'],
            'species': ['Species_A']
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        # FASTA with unparseable header
        fasta_path = self.temp_dir / 'sequences.fasta'
        with open(fasta_path, 'w') as f:
            f.write(">no_processid_here\nACGT\n")

        consensus_records = [
            SeqRecord(Seq('ACGT'), id='consensus_c1', description=''),
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        # Run assignment
        output_path = self.temp_dir / 'output.tsv'
        results = genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path)
        )

        self.assertEqual(results['no_sequence'], 1)

    def test_below_threshold_assignments(self):
        """Test that sequences below threshold are marked as unassigned."""
        metadata = pd.DataFrame({
            'processid': ['SAMPLE001'],
            'species': ['Species_A']
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        # Sequence very different from consensus
        raw_records = [
            SeqRecord(Seq('TTTTTTTT'), id='Species_A_SAMPLE001.COI', description=''),
        ]
        fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, fasta_path, 'fasta')

        consensus_records = [
            SeqRecord(Seq('AAAAAAAA'), id='consensus_c1', description=''),
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        # Run with high threshold
        output_path = self.temp_dir / 'output.tsv'
        results = genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path),
            min_identity=0.90
        )

        self.assertEqual(results['below_threshold'], 1)


class TestDiagnostics(unittest.TestCase):
    """Test diagnostics generation."""

    def setUp(self):
        """Create temporary directory and test data."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_diagnostics_file_created(self):
        """Test that diagnostics file is created when requested."""
        metadata = pd.DataFrame({
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'species': ['Species_A', 'Species_A']
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        raw_records = [
            SeqRecord(Seq('ACGTACGT'), id='Species_A_SAMPLE001.COI', description=''),
            SeqRecord(Seq('ACGTACGA'), id='Species_A_SAMPLE002.COI', description=''),
        ]
        fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, fasta_path, 'fasta')

        consensus_records = [
            SeqRecord(Seq('ACGTACGT'), id='consensus_c1', description=''),
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        output_path = self.temp_dir / 'output.tsv'
        diagnostics_path = self.temp_dir / 'diagnostics.csv'

        genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path),
            diagnostics_path=str(diagnostics_path)
        )

        # Check diagnostics file exists
        self.assertTrue(diagnostics_path.exists())

        # Read and validate diagnostics
        diag_df = pd.read_csv(diagnostics_path)
        self.assertEqual(len(diag_df), 2)
        self.assertIn('processid', diag_df.columns)
        self.assertIn('consensus_group', diag_df.columns)
        self.assertIn('identity', diag_df.columns)
        self.assertIn('runner_up_group', diag_df.columns)
        self.assertIn('is_tie', diag_df.columns)
        self.assertIn('is_low_confidence', diag_df.columns)
        self.assertIn('status', diag_df.columns)

    def test_diagnostics_tie_detection(self):
        """Test that ties are correctly flagged in diagnostics."""
        metadata = pd.DataFrame({
            'processid': ['SAMPLE001'],
            'species': ['Species_A']
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        # Create query that's equidistant from two consensus
        raw_records = [
            SeqRecord(Seq('ACGTACGTACGTACGG'), id='Species_A_SAMPLE001.COI', description=''),
        ]
        fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, fasta_path, 'fasta')

        # Two consensus sequences equidistant from query
        consensus_records = [
            SeqRecord(Seq('ACGTACGTACGTACGT'), id='consensus_c1', description=''),
            SeqRecord(Seq('ACGTACGTACGTACGA'), id='consensus_c2', description=''),
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        output_path = self.temp_dir / 'output.tsv'
        diagnostics_path = self.temp_dir / 'diagnostics.csv'

        genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path),
            diagnostics_path=str(diagnostics_path),
            min_identity=0.80
        )

        # Check for tie flag
        diag_df = pd.read_csv(diagnostics_path)
        self.assertTrue(diag_df.loc[0, 'is_tie'])

    def test_diagnostics_low_confidence_detection(self):
        """Test that low-confidence assignments are flagged."""
        metadata = pd.DataFrame({
            'processid': ['SAMPLE001'],
            'species': ['Species_A']
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        # Sequence with ~93.75% identity (barely above 90% threshold)
        # 1 bp different from consensus -> 93.75% identity (in low-confidence range 90-95%)
        raw_records = [
            SeqRecord(Seq('ACGTACGTACGTACGA'), id='Species_A_SAMPLE001.COI', description=''),
        ]
        fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, fasta_path, 'fasta')

        consensus_records = [
            SeqRecord(Seq('ACGTACGTACGTACGT'), id='consensus_c1', description=''),
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        output_path = self.temp_dir / 'output.tsv'
        diagnostics_path = self.temp_dir / 'diagnostics.csv'

        genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path),
            diagnostics_path=str(diagnostics_path),
            min_identity=0.90
        )

        # Check for low-confidence flag
        diag_df = pd.read_csv(diagnostics_path)
        self.assertTrue(diag_df.loc[0, 'is_low_confidence'])


class TestIntegration(unittest.TestCase):
    """Integration tests using real test data."""

    @classmethod
    def setUpClass(cls):
        """Set up test data paths."""
        cls.test_data_dir = Path(__file__).parent / "data"
        cls.test_fasta = cls.test_data_dir / "test_sequences.fasta"

    def setUp(self):
        """Create temporary directory for test outputs."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_integration_with_test_fasta(self):
        """Test full workflow with test FASTA from dereplication tests."""
        if not self.test_fasta.exists():
            self.skipTest(f"Test FASTA not found: {self.test_fasta}")

        # Read test sequences
        records = list(SeqIO.parse(self.test_fasta, 'fasta'))
        self.assertGreater(len(records), 0)

        # Create metadata from test sequences
        metadata = pd.DataFrame({
            'processid': [rec.id.replace('_cluster1', '').replace('_cluster2', '').replace('_cluster3', '')
                         for rec in records],
            'species': ['Test_species'] * len(records)
        })
        metadata_path = self.temp_dir / 'metadata.tsv'
        metadata.to_csv(metadata_path, sep='\t', index=False)

        # Create raw FASTA with processids in headers
        raw_records = []
        for rec in records:
            processid = rec.id.replace('_cluster1', '').replace('_cluster2', '').replace('_cluster3', '')
            new_id = f"Test_species_{processid}.COI"
            raw_records.append(SeqRecord(rec.seq, id=new_id, description=''))

        fasta_path = self.temp_dir / 'sequences.fasta'
        SeqIO.write(raw_records, fasta_path, 'fasta')

        # Create consensus sequences (one per cluster)
        consensus_records = [
            SeqRecord(records[0].seq, id='consensus_c1', description=''),  # Cluster 1
            SeqRecord(records[4].seq, id='consensus_c2', description=''),  # Cluster 2
            SeqRecord(records[7].seq, id='consensus_c3', description=''),  # Cluster 3
        ]
        consensus_path = self.temp_dir / 'consensus.fasta'
        SeqIO.write(consensus_records, consensus_path, 'fasta')

        # Run assignment
        output_path = self.temp_dir / 'output.tsv'
        diagnostics_path = self.temp_dir / 'diagnostics.csv'

        results = genotype_assignment.assign_genotypes(
            metadata_path=str(metadata_path),
            fasta_path=str(fasta_path),
            consensus_path=str(consensus_path),
            output_path=str(output_path),
            diagnostics_path=str(diagnostics_path),
            min_identity=0.90
        )

        # Verify results
        self.assertEqual(results['total'], len(records))
        self.assertGreater(results['assigned'], 0)

        # Verify output files
        self.assertTrue(output_path.exists())
        self.assertTrue(diagnostics_path.exists())

        # Verify output format
        output_df = pd.read_csv(output_path, sep='\t')
        self.assertIn('consensus_group', output_df.columns)
        self.assertEqual(len(output_df), len(records))


class TestErrorHandling(unittest.TestCase):
    """Test error handling."""

    def setUp(self):
        """Create temporary directory."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up temporary directory."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def test_missing_metadata_file(self):
        """Test error when metadata file doesn't exist."""
        with self.assertRaises(FileNotFoundError):
            genotype_assignment.assign_genotypes(
                metadata_path="nonexistent.tsv",
                fasta_path="fake.fasta",
                consensus_path="fake.fasta",
                output_path="output.tsv"
            )

    def test_missing_fasta_file(self):
        """Test error when FASTA file doesn't exist."""
        metadata_path = self.temp_dir / 'metadata.tsv'
        pd.DataFrame({'processid': ['A']}).to_csv(metadata_path, sep='\t', index=False)

        with self.assertRaises(FileNotFoundError):
            genotype_assignment.assign_genotypes(
                metadata_path=str(metadata_path),
                fasta_path="nonexistent.fasta",
                consensus_path="fake.fasta",
                output_path="output.tsv"
            )

    def test_missing_processid_column(self):
        """Test error when metadata missing processid column."""
        metadata_path = self.temp_dir / 'metadata.tsv'
        pd.DataFrame({'sample_id': ['A']}).to_csv(metadata_path, sep='\t', index=False)

        fasta_path = self.temp_dir / 'seq.fasta'
        fasta_path.write_text(">seq1\nACGT\n")

        consensus_path = self.temp_dir / 'cons.fasta'
        consensus_path.write_text(">c1\nACGT\n")

        with self.assertRaises(genotype_assignment.GenotypeAssignmentError):
            genotype_assignment.assign_genotypes(
                metadata_path=str(metadata_path),
                fasta_path=str(fasta_path),
                consensus_path=str(consensus_path),
                output_path=str(self.temp_dir / 'output.tsv')
            )

    def test_invalid_min_identity(self):
        """Test error with invalid min_identity parameter."""
        metadata_path = self.temp_dir / 'metadata.tsv'
        pd.DataFrame({'processid': ['A']}).to_csv(metadata_path, sep='\t', index=False)

        fasta_path = self.temp_dir / 'seq.fasta'
        fasta_path.write_text(">seq1\nACGT\n")

        consensus_path = self.temp_dir / 'cons.fasta'
        consensus_path.write_text(">c1\nACGT\n")

        with self.assertRaises(ValueError):
            genotype_assignment.assign_genotypes(
                metadata_path=str(metadata_path),
                fasta_path=str(fasta_path),
                consensus_path=str(consensus_path),
                output_path=str(self.temp_dir / 'output.tsv'),
                min_identity=1.5  # Invalid
            )

    def test_invalid_n_processes(self):
        """Test error with invalid n_processes parameter."""
        metadata_path = self.temp_dir / 'metadata.tsv'
        pd.DataFrame({'processid': ['A']}).to_csv(metadata_path, sep='\t', index=False)

        fasta_path = self.temp_dir / 'seq.fasta'
        fasta_path.write_text(">seq1\nACGT\n")

        consensus_path = self.temp_dir / 'cons.fasta'
        consensus_path.write_text(">c1\nACGT\n")

        with self.assertRaises(ValueError):
            genotype_assignment.assign_genotypes(
                metadata_path=str(metadata_path),
                fasta_path=str(fasta_path),
                consensus_path=str(consensus_path),
                output_path=str(self.temp_dir / 'output.tsv'),
                n_processes=0  # Invalid
            )


if __name__ == '__main__':
    unittest.main()
