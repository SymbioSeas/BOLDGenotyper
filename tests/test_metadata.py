"""
Unit tests for boldgenotyper.metadata module

Tests cover:
1. BOLD TSV parsing with various edge cases
2. Coordinate extraction and validation
3. Coordinate quality filtering
4. Edit distance calculations
5. Error handling for malformed data
6. Genotype assignment workflow

Author: Steph Smith (steph.smith@unc.edu)
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os

# Import functions to test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from boldgenotyper.metadata import (
    parse_bold_tsv,
    validate_required_columns,
    extract_coordinates,
    parse_coordinates_column,
    filter_by_coordinate_quality,
    get_coordinate_quality_stats,
    levenshtein_distance,
    calculate_edit_distance,
    validate_metadata_for_analysis,
)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def test_bold_tsv(test_data_dir):
    """Return path to test BOLD TSV file."""
    return test_data_dir / "test_bold.tsv"


@pytest.fixture
def sample_dataframe():
    """Create a small sample DataFrame for testing."""
    data = {
        'processid': ['TEST001', 'TEST002', 'TEST003', 'TEST004', 'TEST005'],
        'nuc': [
            'ACTGACTGACTG',
            'ACTGACTGACTG',
            'ACTGGGGGACTG',
            'ACTGACTGACTG',
            ''
        ],
        'coord': [
            '[23.5863, 81.173004]',
            '[0, 0]',
            '[-15.7835, -47.8912]',
            '',
            None
        ],
        'coord_source': [
            'GPS',
            None,
            'GPS',
            'Country Centroid',
            None
        ],
        'country/ocean': ['Mexico', 'Unknown', 'Brazil', 'Mexico', 'Australia']
    }
    return pd.DataFrame(data)


@pytest.fixture
def temp_tsv_file():
    """Create a temporary TSV file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        # Write header
        f.write('processid\tnuc\tcoord\n')
        # Write data
        f.write('TEST001\tACTG\t[10.5, -20.3]\n')
        f.write('TEST002\tGCTA\t[0, 0]\n')
        temp_path = f.name

    yield temp_path

    # Cleanup
    if os.path.exists(temp_path):
        os.unlink(temp_path)


# ============================================================================
# Test TSV Parsing
# ============================================================================

class TestBoldTSVParsing:
    """Tests for BOLD TSV file parsing."""

    def test_parse_valid_tsv(self, test_bold_tsv):
        """Test parsing a valid BOLD TSV file."""
        df = parse_bold_tsv(test_bold_tsv)

        # Check basic properties
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 20  # 20 test samples
        assert 'processid' in df.columns
        assert 'nuc' in df.columns

        # Check processids are unique
        assert df['processid'].nunique() == 20

    def test_parse_tsv_required_columns(self, test_bold_tsv):
        """Test that required columns are validated."""
        df = parse_bold_tsv(test_bold_tsv, required_columns=['processid', 'nuc'])
        assert 'processid' in df.columns
        assert 'nuc' in df.columns

    def test_parse_tsv_missing_required_columns(self, temp_tsv_file):
        """Test error when required column is missing."""
        with pytest.raises(ValueError, match="missing required columns"):
            parse_bold_tsv(temp_tsv_file, required_columns=['processid', 'missing_col'])

    def test_parse_nonexistent_file(self):
        """Test error when file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            parse_bold_tsv("nonexistent_file.tsv")

    def test_duplicate_processids(self):
        """Test error when processids are not unique."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('processid\tnuc\n')
            f.write('TEST001\tACTG\n')
            f.write('TEST001\tGCTA\n')  # Duplicate
            temp_path = f.name

        try:
            with pytest.raises(ValueError, match="duplicates"):
                parse_bold_tsv(temp_path)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_encoding_fallback(self):
        """Test that latin-1 encoding fallback works."""
        # This test would require a file with latin-1 encoding
        # For now, just verify UTF-8 works
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv',
                                         delete=False, encoding='utf-8') as f:
            f.write('processid\tnuc\n')
            f.write('TEST001\tACTG\n')
            temp_path = f.name

        try:
            df = parse_bold_tsv(temp_path, encoding='utf-8')
            assert len(df) == 1
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_whitespace_cleanup(self):
        """Test that whitespace is cleaned up."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('processid\tnuc\n')
            f.write('  TEST001  \t  ACTG  \n')  # Extra whitespace
            temp_path = f.name

        try:
            df = parse_bold_tsv(temp_path)
            assert df.iloc[0]['processid'] == 'TEST001'
            assert df.iloc[0]['nuc'] == 'ACTG'
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)


class TestColumnValidation:
    """Tests for column validation."""

    def test_validate_required_columns_present(self, sample_dataframe):
        """Test validation when all required columns are present."""
        result = validate_required_columns(sample_dataframe, ['processid', 'nuc'])
        assert result is True

    def test_validate_required_columns_missing(self, sample_dataframe):
        """Test validation when required columns are missing."""
        with pytest.raises(ValueError, match="missing required columns"):
            validate_required_columns(sample_dataframe, ['processid', 'missing_col'])


# ============================================================================
# Test Coordinate Extraction
# ============================================================================

class TestCoordinateExtraction:
    """Tests for coordinate parsing and extraction."""

    def test_extract_standard_format(self):
        """Test extraction of standard [lat, lon] format."""
        coords = extract_coordinates('[23.5863, 81.173004]')
        assert coords == (23.5863, 81.173004)

    def test_extract_without_brackets(self):
        """Test extraction without brackets."""
        coords = extract_coordinates('23.5863, 81.173004')
        assert coords == (23.5863, 81.173004)

    def test_extract_space_separated(self):
        """Test extraction of space-separated coordinates."""
        coords = extract_coordinates('23.5863 81.173004')
        assert coords == (23.5863, 81.173004)

    def test_extract_negative_coordinates(self):
        """Test extraction of negative coordinates."""
        coords = extract_coordinates('[-15.7835, -47.8912]')
        assert coords == (-15.7835, -47.8912)

    def test_extract_zero_coordinates(self):
        """Test extraction of [0, 0] coordinates."""
        coords = extract_coordinates('[0, 0]')
        assert coords == (0.0, 0.0)

    def test_extract_invalid_format(self):
        """Test that invalid formats return None."""
        assert extract_coordinates('invalid') is None
        assert extract_coordinates('23.5') is None  # Only one coordinate
        assert extract_coordinates('[23.5, 81.173, 45]') is None  # Too many

    def test_extract_empty_string(self):
        """Test that empty strings return None."""
        assert extract_coordinates('') is None
        assert extract_coordinates('   ') is None

    def test_extract_none(self):
        """Test that None returns None."""
        assert extract_coordinates(None) is None

    def test_extract_nan(self):
        """Test that NaN returns None."""
        assert extract_coordinates(float('nan')) is None
        assert extract_coordinates(np.nan) is None

    def test_extract_out_of_range_latitude(self):
        """Test that out-of-range latitudes return None."""
        assert extract_coordinates('[100, 0]') is None  # Lat > 90
        assert extract_coordinates('[-100, 0]') is None  # Lat < -90

    def test_extract_out_of_range_longitude(self):
        """Test that out-of-range longitudes return None."""
        assert extract_coordinates('[0, 200]') is None  # Lon > 180
        assert extract_coordinates('[0, -200]') is None  # Lon < -180

    def test_extract_edge_values(self):
        """Test extraction at edge of valid ranges."""
        assert extract_coordinates('[90, 180]') == (90.0, 180.0)
        assert extract_coordinates('[-90, -180]') == (-90.0, -180.0)


class TestCoordinateParsing:
    """Tests for parsing coordinates in DataFrame."""

    def test_parse_coordinates_column(self, sample_dataframe):
        """Test parsing coordinates from DataFrame column."""
        df = parse_coordinates_column(sample_dataframe)

        assert 'latitude' in df.columns
        assert 'longitude' in df.columns

        # Check first row (valid coordinates)
        assert df.iloc[0]['latitude'] == 23.5863
        assert df.iloc[0]['longitude'] == 81.173004

        # Check second row ([0, 0])
        assert df.iloc[1]['latitude'] == 0.0
        assert df.iloc[1]['longitude'] == 0.0

        # Check rows with missing coordinates
        assert pd.isna(df.iloc[3]['latitude'])
        assert pd.isna(df.iloc[4]['latitude'])

    def test_parse_missing_coord_column(self, sample_dataframe):
        """Test parsing when coord column is missing."""
        df = sample_dataframe.drop(columns=['coord'])
        df = parse_coordinates_column(df)

        # Should add columns filled with NaN
        assert 'latitude' in df.columns
        assert 'longitude' in df.columns
        assert df['latitude'].isna().all()


# ============================================================================
# Test Coordinate Filtering
# ============================================================================

class TestCoordinateFiltering:
    """Tests for coordinate quality filtering."""

    def test_filter_missing_coordinates(self):
        """Test filtering of missing coordinates."""
        df = pd.DataFrame({
            'processid': ['A', 'B', 'C'],
            'latitude': [10.0, np.nan, 20.0],
            'longitude': [20.0, 30.0, np.nan]
        })

        filtered = filter_by_coordinate_quality(df, exclude_missing=True)
        assert len(filtered) == 1  # Only row 0 has both coords

    def test_filter_zero_coordinates(self):
        """Test filtering of [0, 0] coordinates."""
        df = pd.DataFrame({
            'processid': ['A', 'B', 'C'],
            'latitude': [10.0, 0.0, 20.0],
            'longitude': [20.0, 0.0, 30.0],
            'coord_source': ['GPS', 'GPS', 'GPS']
        })

        filtered = filter_by_coordinate_quality(df, exclude_zero_coords=True)
        assert len(filtered) == 2  # Exclude row with [0, 0]
        assert filtered.iloc[0]['processid'] == 'A'
        assert filtered.iloc[1]['processid'] == 'C'

    def test_filter_centroid_coordinates(self):
        """Test filtering of centroid coordinates."""
        df = pd.DataFrame({
            'processid': ['A', 'B', 'C', 'D'],
            'latitude': [10.0, 20.0, 30.0, 40.0],
            'longitude': [20.0, 30.0, 40.0, 50.0],
            'coord_source': ['GPS', 'centroid', 'Country Centroid', 'GPS']
        })

        filtered = filter_by_coordinate_quality(
            df,
            exclude_centroids=True,
            exclude_missing=False,
            exclude_zero_coords=False
        )
        assert len(filtered) == 2  # Keep A and D
        assert 'B' not in filtered['processid'].values
        assert 'C' not in filtered['processid'].values

    def test_filter_combined_criteria(self, sample_dataframe):
        """Test filtering with all criteria combined."""
        df = parse_coordinates_column(sample_dataframe)
        filtered = filter_by_coordinate_quality(
            df,
            exclude_missing=True,
            exclude_zero_coords=True,
            exclude_centroids=True
        )

        # Should exclude:
        # - TEST002: [0, 0]
        # - TEST004: empty coord
        # - TEST005: None coord
        # Should keep:
        # - TEST001: valid GPS
        # - TEST003: valid GPS

        assert len(filtered) == 2
        assert 'TEST001' in filtered['processid'].values
        assert 'TEST003' in filtered['processid'].values

    def test_filter_no_filtering(self):
        """Test that disabling all filters keeps all rows."""
        df = pd.DataFrame({
            'processid': ['A', 'B'],
            'latitude': [10.0, 0.0],
            'longitude': [20.0, 0.0]
        })

        filtered = filter_by_coordinate_quality(
            df,
            exclude_missing=False,
            exclude_zero_coords=False,
            exclude_centroids=False
        )
        assert len(filtered) == 2


class TestCoordinateStats:
    """Tests for coordinate quality statistics."""

    def test_get_coordinate_quality_stats(self):
        """Test calculation of coordinate quality statistics."""
        df = pd.DataFrame({
            'latitude': [10.0, 0.0, np.nan, 20.0],
            'longitude': [20.0, 0.0, 30.0, 40.0],
            'coord_source': ['GPS', 'GPS', 'centroid', 'GPS']
        })

        stats = get_coordinate_quality_stats(df)

        assert stats['total_samples'] == 4
        assert stats['with_coordinates'] == 3
        assert stats['missing_coordinates'] == 1
        assert stats['zero_coordinates'] == 1
        assert stats['centroid_coordinates'] == 1


# ============================================================================
# Test Edit Distance Calculations
# ============================================================================

class TestEditDistance:
    """Tests for edit distance calculations."""

    def test_levenshtein_identical(self):
        """Test Levenshtein distance for identical sequences."""
        dist = levenshtein_distance('ACTG', 'ACTG')
        assert dist == 0

    def test_levenshtein_one_substitution(self):
        """Test Levenshtein distance with one substitution."""
        dist = levenshtein_distance('ACTG', 'ACTC')
        assert dist == 1

    def test_levenshtein_one_insertion(self):
        """Test Levenshtein distance with one insertion."""
        dist = levenshtein_distance('ACTG', 'ACCTG')
        assert dist == 1

    def test_levenshtein_one_deletion(self):
        """Test Levenshtein distance with one deletion."""
        dist = levenshtein_distance('ACTG', 'ACG')
        assert dist == 1

    def test_levenshtein_empty_strings(self):
        """Test Levenshtein distance with empty strings."""
        assert levenshtein_distance('', '') == 0
        assert levenshtein_distance('ACTG', '') == 4
        assert levenshtein_distance('', 'ACTG') == 4

    def test_levenshtein_completely_different(self):
        """Test Levenshtein distance for completely different sequences."""
        dist = levenshtein_distance('AAAA', 'TTTT')
        assert dist == 4

    def test_calculate_edit_distance_with_identity(self):
        """Test edit distance calculation with identity score."""
        dist, identity = calculate_edit_distance('ACTG', 'ACTG', use_edlib=False)
        assert dist == 0
        assert identity == 1.0

    def test_calculate_edit_distance_one_diff(self):
        """Test edit distance with one difference."""
        dist, identity = calculate_edit_distance('ACTG', 'ACTC', use_edlib=False)
        assert dist == 1
        assert identity == 0.75  # 1 - (1/4)

    def test_calculate_edit_distance_different_lengths(self):
        """Test edit distance with different length sequences."""
        dist, identity = calculate_edit_distance('ACTG', 'ACTGG', use_edlib=False)
        assert dist == 1
        # Identity based on max length: 1 - (1/5) = 0.8
        assert abs(identity - 0.8) < 0.01

    def test_calculate_edit_distance_case_insensitive(self):
        """Test that edit distance is case-insensitive."""
        dist, identity = calculate_edit_distance('actg', 'ACTG', use_edlib=False)
        assert dist == 0
        assert identity == 1.0


# ============================================================================
# Test Metadata Validation
# ============================================================================

class TestMetadataValidation:
    """Tests for metadata validation."""

    def test_validate_valid_metadata(self, sample_dataframe):
        """Test validation of valid metadata."""
        is_valid, issues = validate_metadata_for_analysis(sample_dataframe)
        # Will have issues due to empty/missing sequences
        assert isinstance(is_valid, bool)
        assert isinstance(issues, list)

    def test_validate_empty_dataframe(self):
        """Test validation of empty DataFrame."""
        df = pd.DataFrame()
        is_valid, issues = validate_metadata_for_analysis(df)
        assert is_valid is False
        assert 'empty' in ' '.join(issues).lower()

    def test_validate_missing_required_columns(self):
        """Test validation with missing required columns."""
        df = pd.DataFrame({'id': ['A', 'B']})
        is_valid, issues = validate_metadata_for_analysis(df)
        assert is_valid is False
        assert any('missing' in issue.lower() for issue in issues)

    def test_validate_duplicate_processids(self):
        """Test validation with duplicate processids."""
        df = pd.DataFrame({
            'processid': ['A', 'A', 'B'],
            'nuc': ['ACTG', 'GCTA', 'TTAA']
        })
        is_valid, issues = validate_metadata_for_analysis(df)
        assert is_valid is False
        assert any('duplicate' in issue.lower() for issue in issues)

    def test_validate_all_empty_sequences(self):
        """Test validation when all sequences are empty."""
        df = pd.DataFrame({
            'processid': ['A', 'B'],
            'nuc': ['', '']
        })
        is_valid, issues = validate_metadata_for_analysis(df)
        assert is_valid is False
        assert any('empty' in issue.lower() for issue in issues)


# ============================================================================
# Integration Tests
# ============================================================================

class TestIntegration:
    """Integration tests for complete workflows."""

    def test_full_parsing_and_filtering_workflow(self, test_bold_tsv):
        """Test complete workflow from parsing to filtering."""
        # Parse TSV
        df = parse_bold_tsv(test_bold_tsv)
        assert len(df) == 20

        # Parse coordinates
        df = parse_coordinates_column(df)
        assert 'latitude' in df.columns
        assert 'longitude' in df.columns

        # Get stats before filtering
        stats_before = get_coordinate_quality_stats(df)

        # Apply filtering
        df_filtered = filter_by_coordinate_quality(
            df,
            exclude_missing=True,
            exclude_zero_coords=True,
            exclude_centroids=True
        )

        # Check that filtering reduced the dataset
        assert len(df_filtered) < len(df)

        # Verify no [0, 0] coordinates remain
        zero_coords = (df_filtered['latitude'] == 0) & (df_filtered['longitude'] == 0)
        assert not zero_coords.any()

        # Verify no centroid coordinates remain
        if 'coord_source' in df_filtered.columns:
            has_centroid = df_filtered['coord_source'].str.contains(
                'centroid', case=False, na=False
            )
            assert not has_centroid.any()

    def test_coordinate_edge_cases_from_file(self, test_bold_tsv):
        """Test specific edge cases from test file."""
        df = parse_bold_tsv(test_bold_tsv)
        df = parse_coordinates_column(df)

        # TESTID001: Standard format [23.5863, 81.173004]
        row1 = df[df['processid'] == 'TESTID001'].iloc[0]
        assert abs(row1['latitude'] - 23.5863) < 0.0001

        # TESTID004: [0, 0] coordinates
        row4 = df[df['processid'] == 'TESTID004'].iloc[0]
        assert row4['latitude'] == 0.0
        assert row4['longitude'] == 0.0

        # TESTID006: Missing coordinates
        row6 = df[df['processid'] == 'TESTID006'].iloc[0]
        assert pd.isna(row6['latitude'])

        # TESTID008: Space-separated format
        row8 = df[df['processid'] == 'TESTID008'].iloc[0]
        assert abs(row8['latitude'] - 23.5863) < 0.0001

        # TESTID013: Invalid format
        row13 = df[df['processid'] == 'TESTID013'].iloc[0]
        assert pd.isna(row13['latitude'])

        # TESTID016: Out of range (200, 300)
        row16 = df[df['processid'] == 'TESTID016'].iloc[0]
        assert pd.isna(row16['latitude'])


# ============================================================================
# Edge Case Tests
# ============================================================================

class TestEdgeCases:
    """Tests for specific edge cases and error conditions."""

    def test_coordinate_with_extra_whitespace(self):
        """Test coordinate extraction with extra whitespace."""
        coords = extract_coordinates('[  23.5863  ,  81.173004  ]')
        assert coords == (23.5863, 81.173004)

    def test_coordinate_with_parentheses(self):
        """Test coordinate extraction with parentheses instead of brackets."""
        coords = extract_coordinates('(23.5863, 81.173004)')
        assert coords == (23.5863, 81.173004)

    def test_high_precision_coordinates(self):
        """Test extraction of high-precision coordinates."""
        coords = extract_coordinates('[23.58631234, 81.17300456]')
        assert coords is not None
        assert abs(coords[0] - 23.58631234) < 0.00000001

    def test_centroid_case_insensitive(self):
        """Test that centroid detection is case-insensitive."""
        df = pd.DataFrame({
            'processid': ['A', 'B', 'C'],
            'latitude': [10.0, 20.0, 30.0],
            'longitude': [20.0, 30.0, 40.0],
            'coord_source': ['GPS', 'CENTROID', 'Centroid']
        })

        filtered = filter_by_coordinate_quality(df, exclude_centroids=True)
        assert len(filtered) == 1  # Only GPS remains

    def test_very_long_sequence(self):
        """Test edit distance with very long sequences."""
        seq1 = 'ACTG' * 100  # 400bp
        seq2 = 'ACTG' * 100

        dist, identity = calculate_edit_distance(seq1, seq2, use_edlib=False)
        assert dist == 0
        assert identity == 1.0

    def test_sequence_with_gaps(self):
        """Test edit distance with gap characters."""
        # Gaps should be treated as regular characters for edit distance
        dist, identity = calculate_edit_distance('ACT-G', 'ACTAG', use_edlib=False)
        assert dist == 1  # - replaced with A


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
