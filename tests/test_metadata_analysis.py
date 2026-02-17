"""
Unit tests for metadata_analysis.py module.

Tests metadata coverage analysis, categorical field analysis,
temporal analysis, and statistical association tests.

Author: Steph Smith (symbioseas@outlook.com)
"""

import unittest
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
import tempfile
import pandas as pd
import numpy as np
import json
import shutil

from boldgenotyper import metadata_analysis


class TestMetadataCoverage(unittest.TestCase):
    """Tests for metadata coverage analysis functions."""

    def setUp(self):
        """Create sample test data."""
        self.sample_df = pd.DataFrame({
            'processid': ['A', 'B', 'C', 'D', 'E'],
            'haplotype_sp': ['h1', 'h1', 'h2', 'h2', 'h3'],
            'sex': ['M', 'F', 'Male', None, 'F'],
            'life_stage': ['Adult', 'Juvenile', 'Adult', 'Adult', None],
            'country/ocean': ['USA', 'USA', 'Mexico', 'Mexico', 'USA'],
            'realm': ['Nearctic', 'Nearctic', None, 'Neotropic', None],
            'collection_date_start': ['2020-01-15', '2019-06-20', '2018-03-10', None, '2021-08-05'],
        })

    def test_analyze_metadata_coverage_basic(self):
        """Test basic coverage analysis."""
        fields = ['sex', 'life_stage', 'country/ocean', 'realm']
        result = metadata_analysis.analyze_metadata_coverage(self.sample_df, fields)

        self.assertEqual(result['total_samples'], 5)
        self.assertIn('fields', result)
        self.assertIn('sex', result['fields'])

    def test_coverage_calculation(self):
        """Test coverage percentage calculation."""
        fields = ['sex', 'realm']
        result = metadata_analysis.analyze_metadata_coverage(self.sample_df, fields)

        # sex has 4 non-null values out of 5
        self.assertEqual(result['fields']['sex']['n_with_value'], 4)
        self.assertAlmostEqual(result['fields']['sex']['pct_coverage'], 80.0, places=1)

        # realm has 3 non-null values out of 5 (Nearctic, Nearctic, Neotropic)
        self.assertEqual(result['fields']['realm']['n_with_value'], 3)
        self.assertAlmostEqual(result['fields']['realm']['pct_coverage'], 60.0, places=1)

    def test_missing_field_handling(self):
        """Test handling of missing fields."""
        fields = ['sex', 'nonexistent_field']
        result = metadata_analysis.analyze_metadata_coverage(self.sample_df, fields)

        self.assertIn('warning', result['fields']['nonexistent_field'])
        self.assertIn('not found', result['fields']['nonexistent_field']['warning'])


class TestValueVariantDetection(unittest.TestCase):
    """Tests for value variant detection."""

    def test_detect_sex_variants(self):
        """Test detection of sex value variants."""
        values = pd.Series(['M', 'F', 'Male', 'Female', 'M', 'F'])
        variants = metadata_analysis.detect_value_variants(values)

        # Should detect M/Male and F/Female as variants
        self.assertTrue(len(variants) > 0)

    def test_detect_case_variants(self):
        """Test detection of case variants."""
        values = pd.Series(['Adult', 'ADULT', 'adult', 'Juvenile'])
        variants = metadata_analysis.detect_value_variants(values)

        # Should detect Adult/ADULT/adult as variants
        variant_found = any('Adult' in v or 'adult' in v for v in variants)
        self.assertTrue(variant_found)

    def test_no_variants(self):
        """Test when no variants exist."""
        values = pd.Series(['A', 'B', 'C'])
        variants = metadata_analysis.detect_value_variants(values)

        self.assertEqual(len(variants), 0)


class TestSexNormalization(unittest.TestCase):
    """Tests for sex value normalization."""

    def setUp(self):
        """Create test data with various sex values."""
        self.df = pd.DataFrame({
            'processid': ['A', 'B', 'C', 'D', 'E', 'F', 'G'],
            'sex': ['M', 'F', 'Male', 'Female', 'm', 'f', None]
        })

    def test_normalize_sex_male(self):
        """Test normalization of male variants."""
        result = metadata_analysis.normalize_sex_values(self.df.copy())

        male_mask = result['sex'] == 'Male'
        # 'M', 'Male', and 'm' all normalize to 'Male'
        self.assertEqual(male_mask.sum(), 3)

    def test_normalize_sex_female(self):
        """Test normalization of female variants."""
        result = metadata_analysis.normalize_sex_values(self.df.copy())

        female_mask = result['sex'] == 'Female'
        # 'F', 'Female', and 'f' all normalize to 'Female'
        self.assertEqual(female_mask.sum(), 3)

    def test_preserve_null_values(self):
        """Test that null values are preserved."""
        result = metadata_analysis.normalize_sex_values(self.df.copy())

        self.assertTrue(result['sex'].iloc[-1] is None or pd.isna(result['sex'].iloc[-1]))

    def test_missing_column(self):
        """Test handling when sex column is missing."""
        df_no_sex = self.df.drop(columns=['sex'])
        result = metadata_analysis.normalize_sex_values(df_no_sex)

        self.assertNotIn('sex', result.columns)


class TestCategoricalAnalysis(unittest.TestCase):
    """Tests for categorical field analysis."""

    def setUp(self):
        """Create test data."""
        self.df = pd.DataFrame({
            'processid': list('ABCDEFGHIJ'),
            'haplotype_sp': ['h1', 'h1', 'h1', 'h2', 'h2', 'h2', 'h2', 'h3', 'h3', 'h3'],
            'sex': ['M', 'M', 'F', 'M', 'F', 'F', 'F', 'M', 'M', 'F'],
            'realm': ['A', 'A', 'B', 'A', 'A', 'B', 'B', 'A', 'B', 'B']
        })

    def test_analyze_categorical_basic(self):
        """Test basic categorical analysis."""
        result = metadata_analysis.analyze_categorical_field(
            self.df, 'sex', 'haplotype_sp'
        )

        self.assertFalse(result.empty)
        self.assertIn('haplotype_sp', result.columns)
        self.assertIn('value', result.columns)
        self.assertIn('n_samples', result.columns)

    def test_categorical_percentages(self):
        """Test percentage calculations."""
        result = metadata_analysis.analyze_categorical_field(
            self.df, 'sex', 'haplotype_sp'
        )

        # Check that percentages sum correctly per haplotype
        for haplotype in result['haplotype_sp'].unique():
            hap_data = result[result['haplotype_sp'] == haplotype]
            # pct_of_haplotype should sum to ~100%
            self.assertAlmostEqual(hap_data['pct_of_haplotype'].sum(), 100.0, places=0)

    def test_missing_field(self):
        """Test handling of missing field."""
        result = metadata_analysis.analyze_categorical_field(
            self.df, 'nonexistent', 'haplotype_sp'
        )

        self.assertTrue(result.empty)


class TestAssociationTests(unittest.TestCase):
    """Tests for statistical association tests."""

    def setUp(self):
        """Create test data with clear associations."""
        # Create data where sex is strongly associated with haplotype
        self.associated_df = pd.DataFrame({
            'processid': list(range(100)),
            'haplotype_sp': ['h1'] * 50 + ['h2'] * 50,
            'sex': ['M'] * 40 + ['F'] * 10 + ['F'] * 40 + ['M'] * 10,  # Strong association
        })

        # Create data with no association
        self.independent_df = pd.DataFrame({
            'processid': list(range(100)),
            'haplotype_sp': ['h1'] * 50 + ['h2'] * 50,
            'sex': ['M', 'F'] * 50,  # No association
        })

    def test_association_test_returns_expected_fields(self):
        """Test that association test returns expected fields."""
        result = metadata_analysis.test_haplotype_association(
            self.associated_df, 'sex', 'haplotype_sp'
        )

        self.assertIn('field', result)
        self.assertIn('test_type', result)
        self.assertIn('statistic', result)
        self.assertIn('p_value', result)
        self.assertIn('n_samples', result)

    def test_significant_association(self):
        """Test detection of significant association."""
        result = metadata_analysis.test_haplotype_association(
            self.associated_df, 'sex', 'haplotype_sp'
        )

        # Strong association should yield low p-value
        self.assertLess(result['p_value'], 0.05)

    def test_no_association(self):
        """Test when there's no association."""
        result = metadata_analysis.test_haplotype_association(
            self.independent_df, 'sex', 'haplotype_sp'
        )

        # No association should yield higher p-value
        # (though this is probabilistic, with 100 samples should be clear)
        self.assertGreater(result['p_value'], 0.01)

    def test_insufficient_categories(self):
        """Test handling of insufficient categories."""
        df_single_value = pd.DataFrame({
            'processid': ['A', 'B', 'C'],
            'haplotype_sp': ['h1', 'h1', 'h1'],
            'sex': ['M', 'M', 'M']
        })

        result = metadata_analysis.test_haplotype_association(
            df_single_value, 'sex', 'haplotype_sp'
        )

        self.assertIn('warning', result)
        self.assertTrue(pd.isna(result['p_value']))


class TestTemporalAnalysis(unittest.TestCase):
    """Tests for temporal analysis functions."""

    def setUp(self):
        """Create test data with various date formats."""
        self.df = pd.DataFrame({
            'processid': list('ABCDEFGH'),
            'haplotype_sp': ['h1', 'h1', 'h2', 'h2', 'h3', 'h3', 'h3', 'h1'],
            'collection_date_start': [
                '2020-01-15',
                '2019-06-20',
                '2018-03-10',
                '2021-08-05',
                '2020-12-01',
                '12/14/98',  # Different format
                None,
                'invalid_date'
            ],
            'country/ocean': ['USA', 'USA', 'Mexico', 'Mexico', 'USA', 'USA', 'USA', 'Canada']
        })

    def test_parse_dates_basic(self):
        """Test basic date parsing."""
        result_df, stats = metadata_analysis.parse_collection_dates(
            self.df.copy(), 'collection_date_start'
        )

        self.assertIn('parsed_date', result_df.columns)
        self.assertIn('collection_year', result_df.columns)
        self.assertIn('collection_month', result_df.columns)
        self.assertIn('collection_decade', result_df.columns)

    def test_parse_stats_accuracy(self):
        """Test accuracy of parsing statistics."""
        result_df, stats = metadata_analysis.parse_collection_dates(
            self.df.copy(), 'collection_date_start'
        )

        self.assertEqual(stats['n_total'], 8)
        # 7 rows have date strings, but only 6 are valid (1 is 'invalid_date')
        self.assertGreater(stats['n_parsed'], 0)

    def test_haplotype_emergence(self):
        """Test haplotype emergence calculation."""
        result_df, _ = metadata_analysis.parse_collection_dates(
            self.df.copy(), 'collection_date_start'
        )

        emergence_df = metadata_analysis.calculate_haplotype_emergence(
            result_df, 'haplotype_sp'
        )

        self.assertFalse(emergence_df.empty)
        self.assertIn('haplotype_sp', emergence_df.columns)
        self.assertIn('first_collection_date', emergence_df.columns)
        self.assertIn('first_collection_year', emergence_df.columns)

    def test_temporal_distribution(self):
        """Test temporal distribution analysis."""
        result_df, _ = metadata_analysis.parse_collection_dates(
            self.df.copy(), 'collection_date_start'
        )

        temporal_df = metadata_analysis.analyze_temporal_distribution(
            result_df, 'haplotype_sp', 'year'
        )

        self.assertFalse(temporal_df.empty)
        self.assertIn('time_period', temporal_df.columns)
        self.assertIn('haplotype_sp', temporal_df.columns)
        self.assertIn('n_samples', temporal_df.columns)
        self.assertIn('pct_of_period', temporal_df.columns)


class TestVisualization(unittest.TestCase):
    """Tests for visualization functions."""

    def setUp(self):
        """Create test data and temp directory."""
        self.temp_dir = tempfile.mkdtemp()
        self.df = pd.DataFrame({
            'processid': list('ABCDEFGHIJ'),
            'haplotype_sp': ['h1', 'h1', 'h1', 'h2', 'h2', 'h2', 'h2', 'h3', 'h3', 'h3'],
            'sex': ['M', 'M', 'F', 'M', 'F', 'F', 'F', 'M', 'M', 'F'],
        })

    def tearDown(self):
        """Clean up temp directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_get_haplotype_colors(self):
        """Test haplotype color generation."""
        colors = metadata_analysis.get_haplotype_colors(self.df, 'haplotype_sp')

        self.assertIsInstance(colors, dict)
        self.assertEqual(len(colors), 3)  # 3 haplotypes
        for hap, color in colors.items():
            self.assertTrue(color.startswith('#'))

    def test_plot_metadata_coverage(self):
        """Test coverage plot generation."""
        coverage_stats = metadata_analysis.analyze_metadata_coverage(
            self.df, ['sex']
        )

        output_path = Path(self.temp_dir) / 'coverage'
        metadata_analysis.plot_metadata_coverage(coverage_stats, output_path)

        # Check that files were created
        self.assertTrue((Path(self.temp_dir) / 'coverage.png').exists())

    def test_plot_categorical_by_haplotype(self):
        """Test categorical plot generation."""
        output_path = Path(self.temp_dir) / 'categorical'
        colors = metadata_analysis.get_haplotype_colors(self.df, 'haplotype_sp')

        metadata_analysis.plot_categorical_by_haplotype(
            self.df, 'sex', output_path, colors, 'haplotype_sp'
        )

        # Check that files were created
        self.assertTrue((Path(self.temp_dir) / 'categorical.png').exists())


class TestExport(unittest.TestCase):
    """Tests for export functions."""

    def setUp(self):
        """Create test data and temp directory."""
        self.temp_dir = tempfile.mkdtemp()
        self.df = pd.DataFrame({
            'processid': list('ABCDEFGHIJ'),
            'haplotype_sp': ['h1', 'h1', 'h1', 'h2', 'h2', 'h2', 'h2', 'h3', 'h3', 'h3'],
            'sex': ['M', 'M', 'F', 'M', 'F', 'F', 'F', 'M', 'M', 'F'],
            'life_stage': ['Adult'] * 5 + ['Juvenile'] * 5,
            'collection_date_start': ['2020-01-15'] * 5 + ['2019-06-20'] * 5,
            'country/ocean': ['USA'] * 6 + ['Mexico'] * 4,
        })

    def tearDown(self):
        """Clean up temp directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_export_metadata_analysis(self):
        """Test full metadata analysis export."""
        output_dir = Path(self.temp_dir)

        outputs = metadata_analysis.export_metadata_analysis(
            df=self.df,
            output_dir=output_dir,
            organism='TestOrganism',
            fields=['sex', 'life_stage'],
            normalize_sex=False,
            temporal_analysis=True
        )

        # Check that output directory was created
        metadata_dir = output_dir / 'metadata_analysis'
        self.assertTrue(metadata_dir.exists())

        # Check that key files were created
        self.assertIn('coverage_csv', outputs)
        self.assertIn('coverage_json', outputs)

    def test_run_metadata_analysis(self):
        """Test main entry point function."""
        output_dir = Path(self.temp_dir)

        results = metadata_analysis.run_metadata_analysis(
            annotated_df=self.df,
            output_dir=output_dir,
            organism='TestOrganism',
            fields=['sex'],
            normalize_sex=False,
            temporal_analysis=False
        )

        self.assertIn('output_files', results)
        self.assertIn('n_fields_analyzed', results)


class TestLastDetection(unittest.TestCase):
    """Tests for first and last detection tracking."""

    def setUp(self):
        """Create test data with date range."""
        self.df = pd.DataFrame({
            'processid': list('ABCDEFGH'),
            'haplotype_sp': ['h1', 'h1', 'h1', 'h2', 'h2', 'h3', 'h3', 'h3'],
            'collection_date_start': [
                '2015-01-15',  # h1 first
                '2018-06-20',
                '2020-12-01',  # h1 last
                '2016-03-10',  # h2 first
                '2019-08-05',  # h2 last
                '2017-04-01',  # h3 first
                '2021-02-15',
                '2022-11-20',  # h3 last
            ],
            'country/ocean': ['USA', 'USA', 'Mexico', 'Mexico', 'USA', 'USA', 'Canada', 'Canada']
        })

    def test_emergence_includes_last_detection(self):
        """Test that emergence data includes last detection date."""
        result_df, _ = metadata_analysis.parse_collection_dates(
            self.df.copy(), 'collection_date_start'
        )

        emergence_df = metadata_analysis.calculate_haplotype_emergence(
            result_df, 'haplotype_sp'
        )

        self.assertIn('last_collection_date', emergence_df.columns)
        self.assertIn('last_collection_year', emergence_df.columns)
        self.assertIn('years_observed', emergence_df.columns)

    def test_years_observed_calculation(self):
        """Test years_observed is calculated correctly."""
        result_df, _ = metadata_analysis.parse_collection_dates(
            self.df.copy(), 'collection_date_start'
        )

        emergence_df = metadata_analysis.calculate_haplotype_emergence(
            result_df, 'haplotype_sp'
        )

        # h1: 2015 to 2020 = 5 years
        h1_row = emergence_df[emergence_df['haplotype_sp'] == 'h1'].iloc[0]
        self.assertEqual(h1_row['years_observed'], 5)

        # h3: 2017 to 2022 = 5 years
        h3_row = emergence_df[emergence_df['haplotype_sp'] == 'h3'].iloc[0]
        self.assertEqual(h3_row['years_observed'], 5)


class TestMissingFieldsTracking(unittest.TestCase):
    """Tests for tracking missing metadata fields."""

    def setUp(self):
        """Create test data with limited fields."""
        self.df = pd.DataFrame({
            'processid': list('ABCDE'),
            'haplotype_sp': ['h1', 'h1', 'h2', 'h2', 'h3'],
            'sex': ['M', 'F', 'M', 'F', 'M'],
            # Note: missing life_stage, realm, etc.
        })
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temp directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_run_metadata_analysis_returns_missing_fields(self):
        """Test that run_metadata_analysis tracks missing fields."""
        output_dir = Path(self.temp_dir)

        results = metadata_analysis.run_metadata_analysis(
            annotated_df=self.df,
            output_dir=output_dir,
            organism='TestOrganism',
            fields=['sex', 'life_stage', 'realm'],  # life_stage and realm are missing
            normalize_sex=False,
            temporal_analysis=False
        )

        self.assertIn('missing_fields', results)
        self.assertIn('available_fields', results)
        self.assertIn('life_stage', results['missing_fields'])
        self.assertIn('realm', results['missing_fields'])
        self.assertIn('sex', results['available_fields'])


class TestEdgeCases(unittest.TestCase):
    """Tests for edge cases and error handling."""

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        empty_df = pd.DataFrame()

        result = metadata_analysis.analyze_metadata_coverage(empty_df, ['sex'])
        self.assertEqual(result['total_samples'], 0)

    def test_all_null_values(self):
        """Test handling of all-null values."""
        df = pd.DataFrame({
            'haplotype_sp': ['h1', 'h2', 'h3'],
            'sex': [None, None, None]
        })

        coverage = metadata_analysis.analyze_metadata_coverage(df, ['sex'])
        self.assertEqual(coverage['fields']['sex']['n_with_value'], 0)
        self.assertEqual(coverage['fields']['sex']['pct_coverage'], 0.0)

    def test_single_haplotype(self):
        """Test handling of single haplotype."""
        df = pd.DataFrame({
            'haplotype_sp': ['h1', 'h1', 'h1'],
            'sex': ['M', 'F', 'M']
        })

        result = metadata_analysis.test_haplotype_association(df, 'sex', 'haplotype_sp')
        self.assertIn('warning', result)

    def test_special_characters_in_field_names(self):
        """Test handling of fields with special characters."""
        df = pd.DataFrame({
            'haplotype_sp': ['h1', 'h2'],
            'country/ocean': ['USA', 'Mexico']
        })

        coverage = metadata_analysis.analyze_metadata_coverage(df, ['country/ocean'])
        self.assertIn('country/ocean', coverage['fields'])


if __name__ == '__main__':
    unittest.main()
