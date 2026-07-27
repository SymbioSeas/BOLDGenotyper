"""
Tests for early CLI validation of custom shapefile + field arguments.
"""
import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil
import sys


class TestCustomShapefileEarlyValidation(unittest.TestCase):
    """Tests that bad --shp-field fails at startup with a clear message."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.shp_path = Path(self.test_dir) / 'test.shp'
        self.shp_path.touch()

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    @patch('geopandas.read_file')
    def test_bad_shp_field_exits_nonzero(self, mock_read_file):
        """When --shp-field is not in shapefile, sys.exit with nonzero code."""
        import geopandas as gpd
        from shapely.geometry import Point

        mock_gdf = gpd.GeoDataFrame({
            'HYBAS_ID': [1, 2],
            'geometry': [Point(0, 0), Point(1, 1)]
        }, crs='EPSG:4326')
        mock_read_file.return_value = mock_gdf

        from boldgenotyper.cli import validate_shapefile_field

        with self.assertRaises(SystemExit) as cm:
            validate_shapefile_field(self.shp_path, 'wrong_field')

        self.assertNotEqual(cm.exception.code, 0)

    @patch('geopandas.read_file')
    def test_valid_shp_field_passes_validation(self, mock_read_file):
        """When --shp-field matches a shapefile column, validation passes silently."""
        import geopandas as gpd
        from shapely.geometry import Point

        mock_gdf = gpd.GeoDataFrame({
            'HYBAS_ID': [1, 2],
            'geometry': [Point(0, 0), Point(1, 1)]
        }, crs='EPSG:4326')
        mock_read_file.return_value = mock_gdf

        from boldgenotyper.cli import validate_shapefile_field

        # Should not raise or exit
        validate_shapefile_field(self.shp_path, 'HYBAS_ID')

    @patch('geopandas.read_file')
    def test_none_shp_field_passes_validation(self, mock_read_file):
        """When --shp-field is None (auto-detect), validation passes silently."""
        import geopandas as gpd
        from shapely.geometry import Point

        mock_gdf = gpd.GeoDataFrame({
            'HYBAS_ID': [1, 2],
            'geometry': [Point(0, 0), Point(1, 1)]
        }, crs='EPSG:4326')
        mock_read_file.return_value = mock_gdf

        from boldgenotyper.cli import validate_shapefile_field

        # Should not raise or exit
        validate_shapefile_field(self.shp_path, None)


class TestGeoAnalysisFailureSurfacing(unittest.TestCase):
    """Tests that geographic analysis failure produces stderr output, not just log."""

    def test_print_geo_failure_outputs_to_stderr(self):
        """_print_geo_failure_to_stderr writes a visible warning to stderr."""
        import io
        from unittest.mock import patch as mock_patch
        import boldgenotyper.cli as cli_module

        captured_stderr = io.StringIO()

        with mock_patch('sys.stderr', captured_stderr):
            cli_module._print_geo_failure_to_stderr(
                shapefile_path='/path/to/test.shp',
                error=ValueError("Field 'name' not found. Available fields: HYBAS_ID")
            )

        output = captured_stderr.getvalue()
        self.assertIn('geographic analysis', output.lower())
        self.assertIn('HYBAS_ID', output)


if __name__ == '__main__':
    unittest.main()
