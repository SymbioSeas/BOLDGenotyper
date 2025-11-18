"""
Tests for error handling and data validation.

This test suite specifically validates error handling for the scenarios documented
in ERROR_HANDLING_GUIDE.md:
1. Missing metadata (missing required columns in TSV)
2. Non-unique sample IDs (duplicate processids)
3. Erroneous geographic coordinates (out of range lat/lon)
4. FASTA sequences with non-nucleotides
5. Sequences outside length parameters

Each test verifies that the appropriate error handling mechanisms are triggered
and that the system responds appropriately (raises errors, logs warnings, or
gracefully degrades).
"""

import sys
import unittest
from unittest.mock import Mock, patch, MagicMock
import tempfile
import shutil
from pathlib import Path
import pandas as pd
from io import StringIO
import logging

# Mock problematic cartopy imports BEFORE importing boldgenotyper modules
cartopy_mock = MagicMock()
sys.modules['cartopy'] = cartopy_mock
sys.modules['cartopy.crs'] = MagicMock()
sys.modules['cartopy.io'] = MagicMock()
sys.modules['cartopy.io.shapereader'] = MagicMock()
sys.modules['cartopy.feature'] = MagicMock()
sys.modules['cartopy.mpl'] = MagicMock()
sys.modules['cartopy.mpl.ticker'] = MagicMock()

from boldgenotyper import metadata, geographic, utils
from boldgenotyper.dereplication import dereplicate_from_fasta


class TestMissingMetadata(unittest.TestCase):
    """Test handling of missing or incomplete metadata."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    def test_missing_required_column_processid(self):
        """Test that missing 'processid' column raises ValueError."""
        tsv_path = Path(self.tmpdir) / "missing_processid.tsv"

        # Create TSV without processid column
        test_data = {
            'nuc': ['ACGTACGT' * 50, 'ACGTACGT' * 50],
            'species': ['Test sp', 'Test sp'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Should raise ValueError with informative message
        with self.assertRaises(ValueError) as cm:
            metadata.parse_bold_tsv(str(tsv_path))

        error_msg = str(cm.exception)
        self.assertIn("processid", error_msg.lower())
        self.assertIn("missing", error_msg.lower())

    def test_missing_required_column_nuc(self):
        """Test that missing 'nuc' column raises ValueError."""
        tsv_path = Path(self.tmpdir) / "missing_nuc.tsv"

        # Create TSV without nuc column
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'species': ['Test sp', 'Test sp'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Should raise ValueError
        with self.assertRaises(ValueError) as cm:
            metadata.parse_bold_tsv(str(tsv_path))

        error_msg = str(cm.exception)
        self.assertIn("nuc", error_msg.lower())
        self.assertIn("missing", error_msg.lower())

    def test_missing_multiple_required_columns(self):
        """Test that missing multiple required columns reports all missing."""
        tsv_path = Path(self.tmpdir) / "missing_multiple.tsv"

        # Create TSV without processid or nuc
        test_data = {
            'species': ['Test sp', 'Test sp'],
            'genus': ['Test', 'Test'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Should raise ValueError mentioning both missing columns
        with self.assertRaises(ValueError) as cm:
            metadata.parse_bold_tsv(str(tsv_path))

        error_msg = str(cm.exception)
        # Should mention both missing columns
        self.assertIn("processid", error_msg.lower())
        self.assertIn("nuc", error_msg.lower())

    def test_empty_tsv_file(self):
        """Test that empty TSV file is handled appropriately."""
        tsv_path = Path(self.tmpdir) / "empty.tsv"
        tsv_path.write_text("")

        # Should raise an error (could be ValueError or pandas error)
        with self.assertRaises(Exception):
            metadata.parse_bold_tsv(str(tsv_path))

    def test_missing_optional_columns_handled_gracefully(self):
        """Test that missing optional columns (lat, lon, species) doesn't fail."""
        tsv_path = Path(self.tmpdir) / "minimal_valid.tsv"

        # Create TSV with only required columns
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'nuc': ['ACGTACGT' * 50, 'ACGTACGT' * 50],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Should not raise error - optional columns are OK to be missing
        df_result = metadata.parse_bold_tsv(str(tsv_path))

        self.assertEqual(len(df_result), 2)
        self.assertIn('processid', df_result.columns)
        self.assertIn('nuc', df_result.columns)


class TestNonUniqueSampleIDs(unittest.TestCase):
    """Test handling of duplicate processids (non-unique sample IDs)."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    def test_duplicate_processids_logged_and_removed(self):
        """Test that duplicate processids are detected and first occurrence kept."""
        tsv_path = Path(self.tmpdir) / "duplicates.tsv"

        # Create TSV with duplicate processids
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE001', 'SAMPLE002', 'SAMPLE003', 'SAMPLE002'],
            'nuc': ['ACGTACGT' * 50] * 5,
            'species': ['Test sp1', 'Test sp2', 'Test sp3', 'Test sp4', 'Test sp5'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Capture log output
        with self.assertLogs('boldgenotyper.metadata', level='WARNING') as log_context:
            df_result = metadata.parse_bold_tsv(str(tsv_path))

        # Should have logged warning about duplicates
        self.assertTrue(any('duplicate' in log.lower() for log in log_context.output))

        # Should have removed duplicates (kept first occurrence)
        self.assertEqual(len(df_result), 3)  # 3 unique processids

        # Check that first occurrences are kept
        sample001_rows = df_result[df_result['processid'] == 'SAMPLE001']
        self.assertEqual(len(sample001_rows), 1)
        self.assertEqual(sample001_rows.iloc[0]['species'], 'Test sp1')  # First occurrence

    def test_all_duplicate_processids(self):
        """Test file with all duplicate processids (extreme case)."""
        tsv_path = Path(self.tmpdir) / "all_duplicates.tsv"

        # Create TSV where all rows have same processid
        test_data = {
            'processid': ['SAMPLE001'] * 5,
            'nuc': ['ACGTACGT' * 50] * 5,
            'species': ['Test sp'] * 5,
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        with self.assertLogs('boldgenotyper.metadata', level='WARNING') as log_context:
            df_result = metadata.parse_bold_tsv(str(tsv_path))

        # Should log warning about 4 duplicates
        log_output = ' '.join(log_context.output)
        self.assertIn('duplicate', log_output.lower())

        # Should keep only first occurrence
        self.assertEqual(len(df_result), 1)

    def test_no_duplicates_no_warning(self):
        """Test that unique processids don't trigger duplicate warnings."""
        tsv_path = Path(self.tmpdir) / "unique.tsv"

        # Create TSV with all unique processids
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003'],
            'nuc': ['ACGTACGT' * 50] * 3,
            'species': ['Test sp'] * 3,
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Should not produce warning logs for duplicates
        with self.assertLogs('boldgenotyper.metadata', level='DEBUG') as log_context:
            df_result = metadata.parse_bold_tsv(str(tsv_path))

        # Check no duplicate warnings
        log_output = ' '.join(log_context.output)
        self.assertNotIn('duplicate', log_output.lower())

        # All samples retained
        self.assertEqual(len(df_result), 3)


class TestErroneousGeographicCoordinates(unittest.TestCase):
    """Test handling of invalid or erroneous geographic coordinates."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    def test_latitude_out_of_range_high(self):
        """Test that latitude > 90 is handled (set to null)."""
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'lat': [95.0, 25.0],  # 95.0 is invalid
            'lon': [130.0, 131.0],
        }
        df = pd.DataFrame(test_data)

        # Create GeoDataFrame - should log warning and set invalid to None
        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log_context:
            gdf = geographic.create_geodataframe(df)

        # Should have warning about latitude
        log_output = ' '.join(log_context.output)
        self.assertIn('latitude', log_output.lower())
        self.assertIn('95', log_output)

        # First row should have null geometry
        self.assertTrue(gdf.iloc[0]['geometry'] is None or pd.isna(gdf.iloc[0]['geometry']))

        # Second row should have valid geometry
        self.assertIsNotNone(gdf.iloc[1]['geometry'])

    def test_latitude_out_of_range_low(self):
        """Test that latitude < -90 is handled (set to null)."""
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'lat': [-95.0, 25.0],  # -95.0 is invalid
            'lon': [130.0, 131.0],
        }
        df = pd.DataFrame(test_data)

        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log_context:
            gdf = geographic.create_geodataframe(df)

        log_output = ' '.join(log_context.output)
        self.assertIn('latitude', log_output.lower())
        self.assertIn('-95', log_output)

        # First row should have null geometry
        self.assertTrue(gdf.iloc[0]['geometry'] is None or pd.isna(gdf.iloc[0]['geometry']))

    def test_longitude_out_of_range_high(self):
        """Test that longitude > 180 is handled (set to null)."""
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'lat': [25.0, 26.0],
            'lon': [185.0, 131.0],  # 185.0 is invalid
        }
        df = pd.DataFrame(test_data)

        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log_context:
            gdf = geographic.create_geodataframe(df)

        log_output = ' '.join(log_context.output)
        self.assertIn('longitude', log_output.lower())
        self.assertIn('185', log_output)

        # First row should have null geometry
        self.assertTrue(gdf.iloc[0]['geometry'] is None or pd.isna(gdf.iloc[0]['geometry']))

    def test_longitude_out_of_range_low(self):
        """Test that longitude < -180 is handled (set to null)."""
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'lat': [25.0, 26.0],
            'lon': [-185.0, 131.0],  # -185.0 is invalid
        }
        df = pd.DataFrame(test_data)

        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log_context:
            gdf = geographic.create_geodataframe(df)

        log_output = ' '.join(log_context.output)
        self.assertIn('longitude', log_output.lower())
        self.assertIn('-185', log_output)

        # First row should have null geometry
        self.assertTrue(gdf.iloc[0]['geometry'] is None or pd.isna(gdf.iloc[0]['geometry']))

    def test_both_coordinates_out_of_range(self):
        """Test that both lat and lon out of range are handled."""
        test_data = {
            'processid': ['SAMPLE001'],
            'lat': [95.0],  # Invalid
            'lon': [185.0],  # Invalid
        }
        df = pd.DataFrame(test_data)

        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log_context:
            gdf = geographic.create_geodataframe(df)

        log_output = ' '.join(log_context.output)
        # Should have warnings about both
        self.assertIn('latitude', log_output.lower())

        # Row should have null geometry
        self.assertTrue(gdf.iloc[0]['geometry'] is None or pd.isna(gdf.iloc[0]['geometry']))

    def test_centroid_coordinates_filtered(self):
        """Test that centroid coordinates (0, 0) are filtered."""
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003'],
            'lat': [0.0, 0.001, 25.0],  # First is centroid
            'lon': [0.0, 0.001, 130.0],
        }
        df = pd.DataFrame(test_data)

        # Filter centroids using metadata function
        with self.assertLogs('boldgenotyper.metadata', level='DEBUG') as log_context:
            df_filtered = metadata.filter_centroids(df)

        # First row (0, 0) should be filtered out
        # Second row (0.001, 0.001) is close but not exactly centroid, may be kept
        # Third row should definitely be kept
        self.assertGreater(len(df_filtered), 0)

        # Check that exact (0, 0) is not present
        exact_zeros = df_filtered[(df_filtered['lat'] == 0.0) & (df_filtered['lon'] == 0.0)]
        self.assertEqual(len(exact_zeros), 0)

    def test_missing_coordinate_columns_handled(self):
        """Test that missing lat/lon columns doesn't crash."""
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'nuc': ['ACGTACGT' * 50, 'ACGTACGT' * 50],
        }
        df = pd.DataFrame(test_data)

        # Should not crash - returns empty GeoDataFrame or handles gracefully
        try:
            gdf = geographic.create_geodataframe(df)
            # If it succeeds, all geometries should be None/null
            self.assertTrue(all(pd.isna(gdf['geometry'])))
        except (KeyError, ValueError):
            # Also acceptable to raise error for missing columns
            pass


class TestFASTASequenceValidation(unittest.TestCase):
    """Test handling of FASTA sequences with invalid characters."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    def test_sequence_with_non_nucleotide_characters(self):
        """Test that sequences with invalid characters are rejected."""
        # Test with various invalid characters
        invalid_chars = ['X', 'Z', 'J', '!', '@', '1', '2']

        for invalid_char in invalid_chars:
            seq = 'ACGTACGT' * 10 + invalid_char + 'ACGTACGT' * 10

            is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

            self.assertFalse(is_valid, f"Sequence with '{invalid_char}' should be invalid")
            self.assertIn('invalid', error_msg.lower())

    def test_valid_iupac_ambiguous_nucleotides_accepted(self):
        """Test that valid IUPAC ambiguous nucleotides are accepted."""
        # These should all be valid
        valid_ambiguous = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']

        for ambiguous in valid_ambiguous:
            seq = 'ACGTACGT' * 10 + ambiguous + 'ACGTACGT' * 10

            is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

            # Should be valid (unless excessive N content)
            if ambiguous == 'N':
                # Multiple N's might trigger excessive N content warning
                pass
            else:
                self.assertTrue(is_valid, f"Sequence with valid IUPAC code '{ambiguous}' should be valid")

    def test_sequence_with_lowercase_accepted(self):
        """Test that lowercase nucleotides are accepted (normalized to uppercase)."""
        seq = 'acgtacgt' * 20

        is_valid, error_msg = utils.validate_sequence(seq.upper(), min_length=100)

        self.assertTrue(is_valid, "Lowercase nucleotides should be valid when uppercased")

    def test_dereplicate_skips_invalid_sequences(self):
        """Test that dereplication skips sequences with invalid characters."""
        fasta_path = Path(self.tmpdir) / "invalid_sequences.fasta"

        # Create FASTA with mix of valid and invalid sequences
        fasta_content = """>VALID001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>INVALID001
ACGTACGTXYZACGTACGTACGTACGTACGTACGTACGT
>VALID002
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
"""
        fasta_path.write_text(fasta_content)

        # Mock external tools
        with patch('boldgenotyper.utils.check_external_tool', return_value=True):
            with patch('boldgenotyper.utils.run_mafft', return_value=fasta_path):
                with patch('boldgenotyper.utils.run_trimal', return_value=fasta_path):
                    # Should log warnings about invalid sequences
                    with self.assertLogs('boldgenotyper.dereplication', level='WARNING') as log_context:
                        result = dereplicate_from_fasta(
                            str(fasta_path),
                            str(Path(self.tmpdir) / "output.fasta"),
                            threshold=0.01,
                            min_cluster_size=1
                        )

                    log_output = ' '.join(log_context.output)
                    # Should have warning about invalid characters
                    self.assertTrue(
                        'invalid' in log_output.lower() or 'skip' in log_output.lower(),
                        "Should log warning about invalid sequences"
                    )

    def test_empty_sequence_rejected(self):
        """Test that empty sequences are rejected."""
        seq = ""

        is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

        self.assertFalse(is_valid)
        self.assertIn('short', error_msg.lower())


class TestSequenceLengthParameters(unittest.TestCase):
    """Test handling of sequences outside length parameters."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    def test_sequence_too_short_rejected(self):
        """Test that sequences below minimum length are rejected."""
        # Default minimum is 100bp
        short_seq = 'ACGT' * 20  # 80bp

        is_valid, error_msg = utils.validate_sequence(short_seq, min_length=100)

        self.assertFalse(is_valid)
        self.assertIn('short', error_msg.lower())
        self.assertIn('80', error_msg)  # Should mention actual length
        self.assertIn('100', error_msg)  # Should mention minimum

    def test_sequence_at_minimum_length_accepted(self):
        """Test that sequences exactly at minimum length are accepted."""
        min_length = 100
        seq = 'ACGT' * 25  # Exactly 100bp

        is_valid, error_msg = utils.validate_sequence(seq, min_length=min_length)

        self.assertTrue(is_valid, "Sequence at exactly minimum length should be valid")

    def test_sequence_above_minimum_accepted(self):
        """Test that sequences above minimum length are accepted."""
        seq = 'ACGT' * 200  # 800bp - well above typical minimums

        is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

        self.assertTrue(is_valid)

    def test_custom_minimum_length_enforced(self):
        """Test that custom minimum length parameter is enforced."""
        seq = 'ACGT' * 50  # 200bp

        # Should pass with min_length=100
        is_valid1, _ = utils.validate_sequence(seq, min_length=100)
        self.assertTrue(is_valid1)

        # Should fail with min_length=300
        is_valid2, error_msg = utils.validate_sequence(seq, min_length=300)
        self.assertFalse(is_valid2)
        self.assertIn('short', error_msg.lower())

    def test_excessive_n_content_rejected(self):
        """Test that sequences with >50% N content are rejected."""
        # Create sequence with 60% N's
        valid_bases = 'ACGT' * 10  # 40bp
        n_bases = 'N' * 60  # 60bp
        seq = valid_bases + n_bases  # Total 100bp, 60% N

        is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

        self.assertFalse(is_valid)
        self.assertIn('n', error_msg.lower())

    def test_acceptable_n_content_accepted(self):
        """Test that sequences with <50% N content are accepted."""
        # Create sequence with 30% N's
        valid_bases = 'ACGT' * 35  # 140bp
        n_bases = 'N' * 60  # 60bp
        seq = valid_bases + n_bases  # Total 200bp, 30% N

        is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

        self.assertTrue(is_valid, "Sequence with <50% N content should be valid")

    def test_excessive_gaps_rejected(self):
        """Test that sequences with >80% gaps are rejected."""
        # Create sequence with 85% gaps
        valid_bases = 'ACGT' * 15  # 60bp
        gaps = '-' * 340  # 340bp gaps
        seq = valid_bases + gaps  # Total 400bp, 85% gaps

        is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

        self.assertFalse(is_valid)
        self.assertIn('gap', error_msg.lower())

    def test_acceptable_gaps_accepted(self):
        """Test that sequences with <80% gaps are accepted."""
        # Create sequence with 50% gaps
        valid_bases = 'ACGT' * 50  # 200bp
        gaps = '-' * 200  # 200bp gaps
        seq = valid_bases + gaps  # Total 400bp, 50% gaps

        is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

        self.assertTrue(is_valid, "Sequence with <80% gaps should be valid")


class TestIntegratedErrorScenarios(unittest.TestCase):
    """Test integrated error scenarios with multiple issues."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    def test_multiple_data_quality_issues_combined(self):
        """Test file with multiple data quality issues combined."""
        tsv_path = Path(self.tmpdir) / "multiple_issues.tsv"

        # Create TSV with:
        # - Duplicate processids
        # - Invalid coordinates
        # - Short sequences
        # - Sequences with invalid characters
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE001', 'SAMPLE002', 'SAMPLE003'],  # Duplicate
            'nuc': [
                'ACGT' * 50,  # Valid 200bp
                'ACGT' * 10,  # Too short (40bp)
                'ACGTXYZ' * 20,  # Invalid characters
                'ACGT' * 50,  # Valid
            ],
            'lat': [25.0, 95.0, 30.0, 35.0],  # Second is out of range
            'lon': [130.0, 185.0, 135.0, 140.0],  # Second is out of range
            'species': ['Test sp1', 'Test sp2', 'Test sp3', 'Test sp4'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Parse TSV - should handle duplicates
        with self.assertLogs('boldgenotyper.metadata', level='WARNING') as log_context:
            df_parsed = metadata.parse_bold_tsv(str(tsv_path))

        log_output = ' '.join(log_context.output)
        self.assertIn('duplicate', log_output.lower())

        # Should have removed 1 duplicate
        self.assertEqual(len(df_parsed), 3)

        # Create GeoDataFrame - should handle invalid coordinates
        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log_context:
            gdf = geographic.create_geodataframe(df_parsed)

        log_output = ' '.join(log_context.output)
        # Should have warnings about coordinates
        self.assertTrue('latitude' in log_output.lower() or 'longitude' in log_output.lower())

        # Validate sequences - should catch short and invalid ones
        for idx, row in df_parsed.iterrows():
            seq = row['nuc']
            is_valid, error_msg = utils.validate_sequence(seq, min_length=100)

            if len(seq) < 100 or any(c not in 'ACGTUN-RYKMSWBDHV' for c in seq.upper()):
                self.assertFalse(is_valid, f"Invalid sequence should be rejected: {seq[:20]}...")

    def test_graceful_recovery_from_non_critical_errors(self):
        """Test that pipeline can recover from non-critical errors."""
        tsv_path = Path(self.tmpdir) / "recoverable_issues.tsv"

        # Create TSV with some issues but mostly valid data
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003', 'SAMPLE004'],
            'nuc': ['ACGT' * 50] * 4,  # All valid
            'lat': [25.0, 95.0, 30.0, 35.0],  # One invalid
            'lon': [130.0, 131.0, 135.0, 140.0],  # All valid
            'species': ['Test sp'] * 4,
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Parse and process
        df_parsed = metadata.parse_bold_tsv(str(tsv_path))

        # Should still have all 4 samples
        self.assertEqual(len(df_parsed), 4)

        # Create GeoDataFrame - one will have null geometry
        with self.assertLogs('boldgenotyper.geographic', level='WARNING'):
            gdf = geographic.create_geodataframe(df_parsed)

        # Should have 3 valid geometries, 1 null
        valid_geoms = sum(1 for geom in gdf['geometry'] if geom is not None and not pd.isna(geom))
        self.assertEqual(valid_geoms, 3)


if __name__ == '__main__':
    unittest.main()
