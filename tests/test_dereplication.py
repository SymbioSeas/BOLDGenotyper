"""
Unit tests for sequence dereplication module.

Tests cover:
- External tool availability checking
- Distance calculation with known sequence pairs
- Consensus generation with various scenarios
- Sequence clustering with different thresholds
- MAFFT/trimAl wrappers with mocking
- Full integration test with real FASTA data
"""

import unittest
from unittest.mock import Mock, patch, mock_open, MagicMock
import tempfile
import shutil
from pathlib import Path
import subprocess

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from boldgenotyper import dereplication


class TestExternalTools(unittest.TestCase):
    """Test external tool availability checking."""

    def test_check_external_tools_mafft_available(self):
        """Test that MAFFT availability is correctly detected."""
        with patch('shutil.which') as mock_which:
            mock_which.side_effect = lambda x: '/usr/bin/mafft' if x == 'mafft' else None
            tools = dereplication.check_external_tools()
            self.assertTrue(tools['mafft'])
            self.assertFalse(tools['trimal'])

    def test_check_external_tools_trimal_available(self):
        """Test that trimAl availability is correctly detected."""
        with patch('shutil.which') as mock_which:
            mock_which.side_effect = lambda x: '/usr/bin/trimal' if x == 'trimal' else None
            tools = dereplication.check_external_tools()
            self.assertFalse(tools['mafft'])
            self.assertTrue(tools['trimal'])

    def test_check_external_tools_both_available(self):
        """Test when both tools are available."""
        with patch('shutil.which') as mock_which:
            mock_which.return_value = '/usr/bin/tool'
            tools = dereplication.check_external_tools()
            self.assertTrue(tools['mafft'])
            self.assertTrue(tools['trimal'])

    def test_check_external_tools_none_available(self):
        """Test when neither tool is available."""
        with patch('shutil.which') as mock_which:
            mock_which.return_value = None
            tools = dereplication.check_external_tools()
            self.assertFalse(tools['mafft'])
            self.assertFalse(tools['trimal'])


class TestDistanceCalculation(unittest.TestCase):
    """Test pairwise distance calculation."""

    def test_compute_distance_identical_sequences(self):
        """Test distance between identical sequences is 0."""
        seq1 = "ACGTACGT"
        seq2 = "ACGTACGT"
        distance = dereplication._compute_distance(seq1, seq2)
        self.assertEqual(distance, 0.0)

    def test_compute_distance_completely_different(self):
        """Test distance between completely different sequences is 1."""
        seq1 = "AAAAAAAA"
        seq2 = "TTTTTTTT"
        distance = dereplication._compute_distance(seq1, seq2)
        self.assertEqual(distance, 1.0)

    def test_compute_distance_half_different(self):
        """Test distance for sequences with partial differences."""
        seq1 = "ACGTACGT"
        seq2 = "ACGTTTTT"
        distance = dereplication._compute_distance(seq1, seq2)
        # ACGT matches, then A≠T, C≠T, G≠T, T=T
        # 5 matches out of 8 = 62.5% identity, distance = 37.5%
        self.assertEqual(distance, 0.375)

    def test_compute_distance_ignores_gaps(self):
        """Test that gaps are ignored in distance calculation."""
        seq1 = "ACGTACGT"
        seq2 = "ACGT----"
        distance = dereplication._compute_distance(seq1, seq2)
        # Only first 4 positions compared, all match
        self.assertEqual(distance, 0.0)

    def test_compute_distance_ignores_ns(self):
        """Test that N bases are ignored in distance calculation."""
        seq1 = "ACGTACGT"
        seq2 = "ACGTNNNN"
        distance = dereplication._compute_distance(seq1, seq2)
        # Only first 4 positions compared, all match
        self.assertEqual(distance, 0.0)

    def test_compute_distance_uppercase_required(self):
        """Test that _compute_distance requires uppercase input."""
        # _compute_distance is an internal function that expects uppercase
        # The public calculate_pairwise_distances() does the uppercasing
        seq1 = "ACGTACGT"
        seq2 = "ACGTACGT"
        distance = dereplication._compute_distance(seq1, seq2)
        self.assertEqual(distance, 0.0)

        # Verify that lowercase is ignored (not in "ACGT")
        seq1_lower = "acgtacgt"
        seq2_lower = "acgtacgt"
        distance_lower = dereplication._compute_distance(seq1_lower, seq2_lower)
        # No valid sites since lowercase not in "ACGT"
        self.assertEqual(distance_lower, 1.0)

    def test_compute_distance_no_valid_sites(self):
        """Test that distance is 1.0 when no valid comparison sites exist."""
        seq1 = "--------"
        seq2 = "NNNNNNNN"
        distance = dereplication._compute_distance(seq1, seq2)
        self.assertEqual(distance, 1.0)

    def test_calculate_pairwise_distances_shape(self):
        """Test that distance matrix has correct shape."""
        alignment = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
            SeqRecord(Seq("TTTT"), id="seq3"),
        ]
        distances = dereplication.calculate_pairwise_distances(alignment)
        # For 3 sequences: n*(n-1)/2 = 3*2/2 = 3 distances
        self.assertEqual(len(distances), 3)
        self.assertIsInstance(distances, np.ndarray)

    def test_calculate_pairwise_distances_empty_alignment(self):
        """Test that empty alignment raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            dereplication.calculate_pairwise_distances([])
        self.assertIn("empty", str(cm.exception).lower())

    def test_calculate_pairwise_distances_single_sequence(self):
        """Test that single sequence raises ValueError."""
        alignment = [SeqRecord(Seq("ACGT"), id="seq1")]
        with self.assertRaises(ValueError) as cm:
            dereplication.calculate_pairwise_distances(alignment)
        self.assertIn("at least 2", str(cm.exception).lower())

    def test_calculate_pairwise_distances_unequal_lengths(self):
        """Test that sequences of different lengths raise ValueError."""
        alignment = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGTACGT"), id="seq2"),
        ]
        with self.assertRaises(ValueError) as cm:
            dereplication.calculate_pairwise_distances(alignment)
        self.assertIn("same length", str(cm.exception).lower())

    def test_calculate_pairwise_distances_case_insensitive(self):
        """Test that calculate_pairwise_distances handles mixed case."""
        # The public API should handle case conversion
        alignment = [
            SeqRecord(Seq("acgt"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
        ]
        distances = dereplication.calculate_pairwise_distances(alignment)
        # Both sequences are identical after uppercasing
        self.assertEqual(len(distances), 1)
        self.assertEqual(distances[0], 0.0)


class TestConsensusGeneration(unittest.TestCase):
    """Test consensus sequence generation."""

    def test_generate_consensus_identical_sequences(self):
        """Test consensus from identical sequences."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
            SeqRecord(Seq("ACGT"), id="seq3"),
        ]
        consensus = dereplication.generate_consensus(sequences, cluster_id=1)
        self.assertEqual(str(consensus.seq), "ACGT")
        self.assertEqual(consensus.id, "consensus_c1_n3")

    def test_generate_consensus_majority_rule(self):
        """Test majority rule consensus calling."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
            SeqRecord(Seq("ACGT"), id="seq3"),
            SeqRecord(Seq("TCGT"), id="seq4"),  # 75% A, 25% T at position 0
        ]
        consensus = dereplication.generate_consensus(
            sequences, cluster_id=1, frequency_cutoff=0.7
        )
        # A appears in 75% of sequences, meets 70% cutoff
        self.assertEqual(str(consensus.seq), "ACGT")

    def test_generate_consensus_below_cutoff(self):
        """Test that positions below frequency cutoff get N."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
            SeqRecord(Seq("TCGT"), id="seq3"),  # 66.7% A at position 0
        ]
        consensus = dereplication.generate_consensus(
            sequences, cluster_id=1, frequency_cutoff=0.7
        )
        # A appears in 66.7%, below 70% cutoff, so N
        self.assertEqual(str(consensus.seq), "NCGT")

    def test_generate_consensus_ignores_gaps(self):
        """Test that gaps are ignored in consensus generation."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
            SeqRecord(Seq("A---"), id="seq3"),
        ]
        consensus = dereplication.generate_consensus(sequences, cluster_id=1)
        # Position 1-3: only 2 sequences have valid bases, both are CGT
        self.assertEqual(str(consensus.seq), "ACGT")

    def test_generate_consensus_ignores_ns(self):
        """Test that N bases are ignored in consensus generation."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2"),
            SeqRecord(Seq("ANNN"), id="seq3"),
        ]
        consensus = dereplication.generate_consensus(sequences, cluster_id=1)
        self.assertEqual(str(consensus.seq), "ACGT")

    def test_generate_consensus_all_gaps_gives_n(self):
        """Test that position with all gaps gives N."""
        sequences = [
            SeqRecord(Seq("A---"), id="seq1"),
            SeqRecord(Seq("T---"), id="seq2"),
        ]
        consensus = dereplication.generate_consensus(sequences, cluster_id=1)
        # First position: A=50%, T=50%, below 70% cutoff
        # Positions 2-4: all gaps, no valid bases
        self.assertEqual(str(consensus.seq), "NNNN")

    def test_generate_consensus_empty_sequences(self):
        """Test that empty sequence list raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            dereplication.generate_consensus([], cluster_id=1)
        self.assertIn("empty", str(cm.exception).lower())

    def test_generate_consensus_invalid_cutoff(self):
        """Test that invalid frequency cutoff raises ValueError."""
        sequences = [SeqRecord(Seq("ACGT"), id="seq1")]
        with self.assertRaises(ValueError):
            dereplication.generate_consensus(sequences, cluster_id=1, frequency_cutoff=1.5)
        with self.assertRaises(ValueError):
            dereplication.generate_consensus(sequences, cluster_id=1, frequency_cutoff=-0.1)

    def test_generate_consensus_unequal_lengths(self):
        """Test that sequences of different lengths raise ValueError."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGTACGT"), id="seq2"),
        ]
        with self.assertRaises(ValueError) as cm:
            dereplication.generate_consensus(sequences, cluster_id=1)
        self.assertIn("same length", str(cm.exception).lower())

    def test_generate_consensus_cluster_naming(self):
        """Test consensus sequence naming convention."""
        sequences = [SeqRecord(Seq("ACGT"), id=f"seq{i}") for i in range(5)]
        consensus = dereplication.generate_consensus(sequences, cluster_id=42)
        self.assertEqual(consensus.id, "consensus_c42_n5")
        self.assertIn("cluster 42", consensus.description)
        self.assertIn("5 sequences", consensus.description)


class TestClustering(unittest.TestCase):
    """Test sequence clustering."""

    def test_cluster_sequences_basic(self):
        """Test basic clustering functionality."""
        # Create distance matrix for 3 sequences: 2 close, 1 far
        # seq1-seq2: distance 0.0 (identical)
        # seq1-seq3: distance 1.0 (completely different)
        # seq2-seq3: distance 1.0 (completely different)
        distances = np.array([0.0, 1.0, 1.0])
        labels = dereplication.cluster_sequences(distances, threshold=0.01)

        # Should form 2 clusters: {seq1, seq2} and {seq3}
        self.assertEqual(len(labels), 3)
        self.assertEqual(labels[0], labels[1])  # seq1 and seq2 in same cluster
        self.assertNotEqual(labels[0], labels[2])  # seq3 in different cluster

    def test_cluster_sequences_threshold_effect(self):
        """Test that threshold affects clustering."""
        # Sequences with intermediate distances
        # Assume 3 sequences with distances: 0.005, 0.015, 0.02
        distances = np.array([0.005, 0.015, 0.02])

        # Strict threshold: all separate
        labels_strict = dereplication.cluster_sequences(distances, threshold=0.001)
        self.assertEqual(len(set(labels_strict)), 3)

        # Moderate threshold: some cluster
        labels_moderate = dereplication.cluster_sequences(distances, threshold=0.01)
        # First two should cluster (distance 0.005)
        self.assertEqual(labels_moderate[0], labels_moderate[1])

        # Relaxed threshold: all cluster
        labels_relaxed = dereplication.cluster_sequences(distances, threshold=0.1)
        # With average linkage, might still separate depending on structure

    def test_cluster_sequences_linkage_methods(self):
        """Test different linkage methods."""
        distances = np.array([0.01, 0.5, 0.5])

        for method in ['single', 'complete', 'average', 'weighted']:
            labels = dereplication.cluster_sequences(distances, threshold=0.05, method=method)
            self.assertEqual(len(labels), 3)

    def test_cluster_sequences_empty_matrix(self):
        """Test that empty distance matrix raises error."""
        with self.assertRaises(dereplication.ClusteringError):
            dereplication.cluster_sequences(np.array([]))

    def test_cluster_sequences_invalid_threshold(self):
        """Test that invalid threshold raises error."""
        distances = np.array([0.5])
        with self.assertRaises(dereplication.ClusteringError):
            dereplication.cluster_sequences(distances, threshold=-0.1)
        with self.assertRaises(dereplication.ClusteringError):
            dereplication.cluster_sequences(distances, threshold=1.5)

    def test_cluster_sequences_invalid_method(self):
        """Test that invalid linkage method raises error."""
        distances = np.array([0.5])
        with self.assertRaises(dereplication.ClusteringError):
            dereplication.cluster_sequences(distances, threshold=0.01, method='invalid')


class TestMAFFTAlignment(unittest.TestCase):
    """Test MAFFT alignment wrapper with mocking."""

    @patch('shutil.which')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_run_mafft_alignment_success(self, mock_file, mock_run, mock_which):
        """Test successful MAFFT alignment."""
        mock_which.return_value = '/usr/bin/mafft'
        mock_run.return_value = Mock(returncode=0)

        dereplication.run_mafft_alignment(
            input_fasta="input.fasta",
            output_fasta="output.fasta"
        )

        # Check that subprocess.run was called with correct arguments
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        self.assertEqual(args[0], 'mafft')
        self.assertIn('--auto', args)
        self.assertIn('input.fasta', args)

    @patch('shutil.which')
    def test_run_mafft_alignment_not_found(self, mock_which):
        """Test that missing MAFFT raises AlignmentError."""
        mock_which.return_value = None

        with self.assertRaises(dereplication.AlignmentError) as cm:
            dereplication.run_mafft_alignment("input.fasta", "output.fasta")
        self.assertIn("not found", str(cm.exception).lower())

    @patch('shutil.which')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_run_mafft_alignment_failure(self, mock_file, mock_run, mock_which):
        """Test MAFFT alignment failure."""
        mock_which.return_value = '/usr/bin/mafft'
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd='mafft', stderr='Error message'
        )

        with self.assertRaises(dereplication.AlignmentError):
            dereplication.run_mafft_alignment("input.fasta", "output.fasta")

    @patch('shutil.which')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_run_mafft_alignment_custom_options(self, mock_file, mock_run, mock_which):
        """Test MAFFT with custom options."""
        mock_which.return_value = '/usr/bin/mafft'
        mock_run.return_value = Mock(returncode=0)

        custom_options = ['--maxiterate', '1000', '--localpair']
        dereplication.run_mafft_alignment(
            input_fasta="input.fasta",
            output_fasta="output.fasta",
            mafft_options=custom_options
        )

        args = mock_run.call_args[0][0]
        self.assertIn('--maxiterate', args)
        self.assertIn('1000', args)
        self.assertIn('--localpair', args)


class TestTrimAlTrimming(unittest.TestCase):
    """Test trimAl trimming wrapper with mocking."""

    @patch('shutil.which')
    @patch('subprocess.run')
    def test_run_trimal_trimming_success(self, mock_run, mock_which):
        """Test successful trimAl trimming."""
        mock_which.return_value = '/usr/bin/trimal'
        mock_run.return_value = Mock(returncode=0)

        dereplication.run_trimal_trimming(
            input_fasta="input.fasta",
            output_fasta="output.fasta"
        )

        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        self.assertEqual(args[0], 'trimal')
        self.assertIn('-automated1', args)
        self.assertIn('-in', args)
        self.assertIn('-out', args)

    @patch('shutil.which')
    def test_run_trimal_trimming_not_found(self, mock_which):
        """Test that missing trimAl raises TrimmingError."""
        mock_which.return_value = None

        with self.assertRaises(dereplication.TrimmingError) as cm:
            dereplication.run_trimal_trimming("input.fasta", "output.fasta")
        self.assertIn("not found", str(cm.exception).lower())

    @patch('shutil.which')
    @patch('subprocess.run')
    def test_run_trimal_trimming_failure(self, mock_run, mock_which):
        """Test trimAl trimming failure."""
        mock_which.return_value = '/usr/bin/trimal'
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd='trimal', stderr='Error message'
        )

        with self.assertRaises(dereplication.TrimmingError):
            dereplication.run_trimal_trimming("input.fasta", "output.fasta")

    @patch('shutil.which')
    @patch('subprocess.run')
    def test_run_trimal_trimming_custom_options(self, mock_run, mock_which):
        """Test trimAl with custom options."""
        mock_which.return_value = '/usr/bin/trimal'
        mock_run.return_value = Mock(returncode=0)

        custom_options = ['-gt', '0.8', '-st', '0.001']
        dereplication.run_trimal_trimming(
            input_fasta="input.fasta",
            output_fasta="output.fasta",
            trimal_options=custom_options
        )

        args = mock_run.call_args[0][0]
        self.assertIn('-gt', args)
        self.assertIn('0.8', args)
        self.assertIn('-st', args)


class TestDereplicationIntegration(unittest.TestCase):
    """Integration tests for full dereplication workflow."""

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

    def test_integration_with_real_fasta(self):
        """Test full workflow with real test FASTA file."""
        # Skip if external tools not available
        tools = dereplication.check_external_tools()
        if not (tools['mafft'] and tools['trimal']):
            self.skipTest("MAFFT and/or trimAl not available")

        # Skip if test FASTA doesn't exist
        if not self.test_fasta.exists():
            self.skipTest(f"Test FASTA not found: {self.test_fasta}")

        # Read test FASTA and create TSV
        test_tsv = self.temp_dir / "test_data.tsv"
        records = list(SeqIO.parse(self.test_fasta, "fasta"))

        import pandas as pd
        data = {
            'processid': [rec.id for rec in records],
            'nucleotides': [str(rec.seq) for rec in records]
        }
        df = pd.DataFrame(data)
        df.to_csv(test_tsv, sep='\t', index=False)

        # Run dereplication
        output_dir = self.temp_dir / "output"
        consensus_seqs = dereplication.dereplicate_sequences(
            tsv_path=str(test_tsv),
            output_dir=str(output_dir),
            threshold=0.01,
            frequency_cutoff=0.7
        )

        # Verify outputs
        self.assertIsInstance(consensus_seqs, dict)
        self.assertGreater(len(consensus_seqs), 0)

        # Expected: 3 clusters based on test data
        # Cluster 1: seq1-seq4 (4 sequences, ACGT pattern)
        # Cluster 2: seq5-seq7 (3 sequences, TGCA pattern)
        # Cluster 3: seq8-seq10 (3 sequences, GGCC pattern)
        self.assertEqual(len(consensus_seqs), 3,
                        f"Expected 3 clusters, got {len(consensus_seqs)}")

        # Check consensus sequence IDs
        for cons_id in consensus_seqs:
            self.assertTrue(cons_id.startswith('consensus_c'))
            self.assertIn('_n', cons_id)

        # Check that consensus FASTA was created
        consensus_fasta = output_dir / "test_data_consensus.fasta"
        self.assertTrue(consensus_fasta.exists())

        # Verify intermediate files exist (since cleanup not enabled)
        aligned_fasta = output_dir / "test_data_aligned.fasta"
        trimmed_fasta = output_dir / "test_data_trimmed.fasta"
        self.assertTrue(aligned_fasta.exists())
        self.assertTrue(trimmed_fasta.exists())

    def test_integration_with_cleanup(self):
        """Test that cleanup_intermediates removes intermediate files."""
        tools = dereplication.check_external_tools()
        if not (tools['mafft'] and tools['trimal']):
            self.skipTest("MAFFT and/or trimAl not available")

        if not self.test_fasta.exists():
            self.skipTest(f"Test FASTA not found: {self.test_fasta}")

        # Create test TSV
        test_tsv = self.temp_dir / "test_data.tsv"
        records = list(SeqIO.parse(self.test_fasta, "fasta"))

        import pandas as pd
        data = {
            'processid': [rec.id for rec in records],
            'nucleotides': [str(rec.seq) for rec in records]
        }
        df = pd.DataFrame(data)
        df.to_csv(test_tsv, sep='\t', index=False)

        # Run with cleanup
        output_dir = self.temp_dir / "output"
        consensus_seqs = dereplication.dereplicate_sequences(
            tsv_path=str(test_tsv),
            output_dir=str(output_dir),
            threshold=0.01,
            cleanup_intermediates=True
        )

        # Verify consensus file exists
        consensus_fasta = output_dir / "test_data_consensus.fasta"
        self.assertTrue(consensus_fasta.exists())

        # Verify intermediate files were removed
        aligned_fasta = output_dir / "test_data_aligned.fasta"
        trimmed_fasta = output_dir / "test_data_trimmed.fasta"
        sequences_fasta = output_dir / "test_data_sequences.fasta"
        self.assertFalse(aligned_fasta.exists())
        self.assertFalse(trimmed_fasta.exists())
        self.assertFalse(sequences_fasta.exists())

    def test_dereplicate_missing_tsv(self):
        """Test that missing TSV file raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            dereplication.dereplicate_sequences(
                tsv_path="nonexistent.tsv",
                output_dir=str(self.temp_dir)
            )

    def test_dereplicate_missing_external_tools(self):
        """Test that missing external tools raises DereplicationError."""
        # Create minimal TSV
        test_tsv = self.temp_dir / "test.tsv"
        with open(test_tsv, 'w') as f:
            f.write("processid\tnucleotides\nseq1\tACGT\n")

        with patch('boldgenotyper.dereplication.check_external_tools') as mock_check:
            mock_check.return_value = {'mafft': False, 'trimal': False}
            with self.assertRaises(dereplication.DereplicationError):
                dereplication.dereplicate_sequences(
                    tsv_path=str(test_tsv),
                    output_dir=str(self.temp_dir)
                )

    def test_dereplicate_empty_tsv(self):
        """Test that TSV with no sequences raises error."""
        test_tsv = self.temp_dir / "empty.tsv"
        with open(test_tsv, 'w') as f:
            f.write("processid\tnucleotides\n")

        with self.assertRaises(dereplication.DereplicationError):
            dereplication.dereplicate_sequences(
                tsv_path=str(test_tsv),
                output_dir=str(self.temp_dir)
            )

    def test_dereplicate_missing_nucleotides_column(self):
        """Test that TSV without nucleotides column raises error."""
        test_tsv = self.temp_dir / "bad.tsv"
        with open(test_tsv, 'w') as f:
            f.write("processid\tother\nseq1\tvalue\n")

        with self.assertRaises(dereplication.DereplicationError):
            dereplication.dereplicate_sequences(
                tsv_path=str(test_tsv),
                output_dir=str(self.temp_dir)
            )


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_consensus_with_single_sequence(self):
        """Test consensus generation with single sequence."""
        sequences = [SeqRecord(Seq("ACGT"), id="seq1")]
        consensus = dereplication.generate_consensus(sequences, cluster_id=1)
        self.assertEqual(str(consensus.seq), "ACGT")
        self.assertEqual(consensus.id, "consensus_c1_n1")

    def test_distance_calculation_two_sequences(self):
        """Test distance calculation with minimum number of sequences."""
        alignment = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("ACGT"), id="seq2")
        ]
        distances = dereplication.calculate_pairwise_distances(alignment)
        self.assertEqual(len(distances), 1)
        self.assertEqual(distances[0], 0.0)

    def test_clustering_single_distance(self):
        """Test clustering with single pairwise distance."""
        distances = np.array([0.005])  # Distance between 2 sequences
        labels = dereplication.cluster_sequences(distances, threshold=0.01)
        # Both sequences should be in same cluster
        self.assertEqual(len(labels), 2)
        self.assertEqual(labels[0], labels[1])

    def test_consensus_all_ambiguous(self):
        """Test consensus when all positions are ambiguous."""
        sequences = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("TGCA"), id="seq2"),
        ]
        consensus = dereplication.generate_consensus(
            sequences, cluster_id=1, frequency_cutoff=0.9
        )
        # Each position is 50-50, below 90% cutoff
        self.assertEqual(str(consensus.seq), "NNNN")

    def test_very_long_sequences(self):
        """Test with longer sequences similar to real COI barcodes."""
        # Create 500bp sequences
        seq_length = 500
        sequences = [
            SeqRecord(Seq("A" * seq_length), id="seq1"),
            SeqRecord(Seq("A" * seq_length), id="seq2"),
        ]
        distances = dereplication.calculate_pairwise_distances(sequences)
        self.assertEqual(len(distances), 1)
        self.assertEqual(distances[0], 0.0)


if __name__ == '__main__':
    unittest.main()
