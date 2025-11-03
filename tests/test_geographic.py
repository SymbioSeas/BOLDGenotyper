"""
Unit tests for geographic.py module.

Tests ocean basin assignment, coordinate validation, and spatial operations
using mocked GeoPandas functionality to avoid dependency on actual shapefiles.

Author: Steph Smith (steph.smith@unc.edu)
"""

import unittest
from unittest.mock import Mock, MagicMock, patch, PropertyMock
from pathlib import Path
import tempfile
import pytest
import pandas as pd
import numpy as np

# Import the module
from boldgenotyper import geographic


# Mock shapely Point for tests without geopandas
class MockPoint:
    """Mock shapely Point geometry."""
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def buffer(self, distance):
        """Mock buffer method."""
        return self

    def __repr__(self):
        return f"POINT ({self.x} {self.y})"


class MockGeometry:
    """Mock geometry series for GeoDataFrame."""
    def __init__(self, geometries):
        self.geometries = geometries

    def contains(self, point):
        """Mock contains method."""
        # Simple mock: check if point is in our predefined basins
        results = []
        for geom in self.geometries:
            # Mock logic: North Atlantic contains points with lon < 0 and lat > 0
            if hasattr(geom, 'basin_name'):
                if geom.basin_name == "North Atlantic Ocean":
                    results.append(point.x < 0 and point.y > 0)
                elif geom.basin_name == "South Pacific Ocean":
                    results.append(point.x < 0 and point.y < 0)
                elif geom.basin_name == "Indian Ocean":
                    results.append(point.x > 0 and point.y < 0)
                else:
                    results.append(False)
            else:
                results.append(False)
        return pd.Series(results)

    def intersects(self, point):
        """Mock intersects method."""
        return self.contains(point)

    def any(self):
        """Mock any method."""
        return any(self.geometries)

    def notna(self):
        """Mock notna method."""
        return pd.Series([g is not None for g in self.geometries])

    @property
    def empty(self):
        """Mock empty property."""
        return len(self.geometries) == 0


def create_mock_goas_data():
    """
    Create a mock GeoDataFrame representing ocean basins.

    Returns
    -------
    Mock GeoDataFrame with three ocean basins:
    - North Atlantic Ocean: lon < 0, lat > 0
    - South Pacific Ocean: lon < 0, lat < 0
    - Indian Ocean: lon > 0, lat < 0
    """
    # Create mock geometries with basin names
    mock_geom1 = Mock()
    mock_geom1.basin_name = "North Atlantic Ocean"

    mock_geom2 = Mock()
    mock_geom2.basin_name = "South Pacific Ocean"

    mock_geom3 = Mock()
    mock_geom3.basin_name = "Indian Ocean"

    # Create mock GeoDataFrame
    mock_gdf = Mock()
    mock_gdf.__len__ = Mock(return_value=3)
    mock_gdf.columns = ['basin_name', 'geometry']
    mock_gdf.crs = Mock()
    mock_gdf.crs.to_epsg = Mock(return_value=4326)

    # Mock geometry series - create a proper Mock that is not None and not empty
    mock_geometry_series = MockGeometry([mock_geom1, mock_geom2, mock_geom3])
    type(mock_gdf).geometry = PropertyMock(return_value=mock_geometry_series)

    # Mock getitem for column selection
    def getitem(self, key):
        if isinstance(key, list):
            # When selecting multiple columns, return a mock with those columns
            subset_mock = Mock()
            subset_mock.__len__ = Mock(return_value=3)
            subset_mock.columns = key
            subset_mock.crs = self.crs
            subset_mock.geometry = mock_geometry_series
            subset_mock.__getitem__ = lambda k: subset_mock
            return subset_mock
        elif key == 'basin_name':
            return pd.Series(["North Atlantic Ocean", "South Pacific Ocean", "Indian Ocean"])
        return self

    mock_gdf.__getitem__ = getitem

    # Mock unique method for basin names
    basin_series = pd.Series(["North Atlantic Ocean", "South Pacific Ocean", "Indian Ocean"])
    basin_series.unique = Mock(return_value=["North Atlantic Ocean", "South Pacific Ocean", "Indian Ocean"])

    return mock_gdf


class TestGeoPandasAvailability(unittest.TestCase):
    """Test GeoPandas availability checking."""

    def test_check_geopandas_available_when_installed(self):
        """Test check passes when GeoPandas is installed."""
        # GeoPandas should be available in test environment
        if geographic.GEOPANDAS_AVAILABLE:
            # Should not raise
            geographic.check_geopandas_available()
        else:
            # Should raise if not available
            with self.assertRaises(geographic.GeospatialLibraryError):
                geographic.check_geopandas_available()

    @patch('boldgenotyper.geographic.GEOPANDAS_AVAILABLE', False)
    def test_check_geopandas_raises_when_not_available(self):
        """Test check raises error when GeoPandas not available."""
        with self.assertRaises(geographic.GeospatialLibraryError) as cm:
            geographic.check_geopandas_available()

        self.assertIn("GeoPandas", str(cm.exception))
        self.assertIn("pip install geopandas", str(cm.exception))


class TestLoadGoasData(unittest.TestCase):
    """Test loading GOaS shapefile data."""

    @patch('boldgenotyper.geographic.GEOPANDAS_AVAILABLE', False)
    def test_load_goas_raises_without_geopandas(self):
        """Test that loading GOaS data raises error without GeoPandas."""
        with self.assertRaises(geographic.GeospatialLibraryError):
            geographic.load_goas_data("fake_path.shp")

    def test_load_goas_raises_for_nonexistent_file(self):
        """Test that loading raises FileNotFoundError for missing file."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        with self.assertRaises(FileNotFoundError) as cm:
            geographic.load_goas_data("nonexistent_shapefile.shp")

        self.assertIn("GOaS shapefile not found", str(cm.exception))
        self.assertIn("setup_goas.py", str(cm.exception))

    @patch('boldgenotyper.geographic.gpd.read_file')
    def test_load_goas_validates_empty_geodataframe(self, mock_read_file):
        """Test that empty GeoDataFrame raises error."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Mock empty GeoDataFrame
        mock_gdf = Mock()
        mock_gdf.__len__ = Mock(return_value=0)
        mock_read_file.return_value = mock_gdf

        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix='.shp', delete=False) as f:
            temp_path = f.name

        try:
            with self.assertRaises(geographic.GOaSDataError) as cm:
                geographic.load_goas_data(temp_path)

            self.assertIn("empty", str(cm.exception))
        finally:
            Path(temp_path).unlink()

    @patch('boldgenotyper.geographic.gpd.read_file')
    def test_load_goas_handles_missing_crs(self, mock_read_file):
        """Test that missing CRS is handled with warning."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Mock GeoDataFrame without CRS
        mock_gdf = create_mock_goas_data()

        # Ensure geometry is not None or empty
        mock_gdf.geometry.__bool__ = Mock(return_value=True)  # geometry is truthy

        # Set CRS to None initially
        mock_gdf.crs = None

        # Mock set_crs to actually set the CRS
        def set_crs_side_effect(*args, **kwargs):
            mock_crs = Mock()
            mock_crs.to_epsg = Mock(return_value=4326)
            mock_gdf.crs = mock_crs

        mock_gdf.set_crs = Mock(side_effect=set_crs_side_effect)
        mock_read_file.return_value = mock_gdf

        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix='.shp', delete=False) as f:
            temp_path = f.name

        try:
            with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log:
                result = geographic.load_goas_data(temp_path)

            # Check that warning was logged
            self.assertTrue(any("no CRS" in message for message in log.output))
            self.assertTrue(mock_gdf.set_crs.called)
        finally:
            Path(temp_path).unlink()

    @patch('boldgenotyper.geographic.gpd.read_file')
    def test_load_goas_detects_basin_name_column(self, mock_read_file):
        """Test that basin name column is detected from common variations."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Mock GeoDataFrame with 'NAME' column instead of 'basin_name'
        mock_gdf = create_mock_goas_data()

        # Ensure geometry is not None or empty
        mock_gdf.geometry.__bool__ = Mock(return_value=True)

        # Set columns to use 'NAME' instead of 'basin_name'
        mock_gdf.columns = ['NAME', 'geometry']

        # Mock rename to return the same mock_gdf with updated columns
        def mock_rename(columns):
            mock_gdf.columns = ['basin_name', 'geometry']
            return mock_gdf

        mock_gdf.rename = Mock(side_effect=mock_rename)
        mock_read_file.return_value = mock_gdf

        # Create temporary file
        with tempfile.NamedTemporaryFile(suffix='.shp', delete=False) as f:
            temp_path = f.name

        try:
            result = geographic.load_goas_data(temp_path)

            # Check that rename was called
            self.assertTrue(mock_gdf.rename.called)
        finally:
            Path(temp_path).unlink()


class TestCreatePointsGeoDataFrame(unittest.TestCase):
    """Test creating point geometries from coordinates."""

    def test_create_points_raises_without_geopandas(self):
        """Test that point creation raises error without GeoPandas."""
        if geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas is available")

        df = pd.DataFrame({'lat': [23.5], 'lon': [81.2]})

        with self.assertRaises(geographic.GeospatialLibraryError):
            geographic.create_points_geodataframe(df)

    def test_create_points_raises_for_missing_columns(self):
        """Test that missing columns raise ValueError."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({'lat': [23.5]})

        with self.assertRaises(ValueError) as cm:
            geographic.create_points_geodataframe(df)

        self.assertIn("lon", str(cm.exception))

    def test_create_points_with_valid_coordinates(self):
        """Test creating points from valid coordinates."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({
            'sample_id': ['S1', 'S2', 'S3'],
            'lat': [23.5, -15.8, 35.7],
            'lon': [81.2, -47.9, -120.1]
        })

        gdf = geographic.create_points_geodataframe(df)

        # Check that GeoDataFrame was created
        self.assertEqual(len(gdf), 3)
        self.assertIn('geometry', gdf.columns)

        # Check that all geometries are valid
        self.assertEqual(gdf.geometry.notna().sum(), 3)

    def test_create_points_handles_missing_coordinates(self):
        """Test that missing coordinates result in null geometries."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({
            'sample_id': ['S1', 'S2', 'S3'],
            'lat': [23.5, np.nan, 35.7],
            'lon': [81.2, -47.9, np.nan]
        })

        gdf = geographic.create_points_geodataframe(df)

        # S1 should have valid geometry
        self.assertIsNotNone(gdf.loc[0, 'geometry'])

        # S2 and S3 should have null geometry
        self.assertIsNone(gdf.loc[1, 'geometry'])
        self.assertIsNone(gdf.loc[2, 'geometry'])

    def test_create_points_validates_coordinate_ranges(self):
        """Test that out-of-range coordinates are flagged."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({
            'sample_id': ['S1', 'S2', 'S3'],
            'lat': [23.5, 100.0, 35.7],  # S2 has invalid latitude
            'lon': [81.2, -47.9, -200.0]  # S3 has invalid longitude
        })

        with self.assertLogs('boldgenotyper.geographic', level='WARNING') as log:
            gdf = geographic.create_points_geodataframe(df)

        # Check that warnings were logged
        self.assertTrue(any("out of range" in message for message in log.output))

        # S1 should have valid geometry
        self.assertIsNotNone(gdf.loc[0, 'geometry'])

        # S2 and S3 should have null geometry due to invalid ranges
        self.assertIsNone(gdf.loc[1, 'geometry'])
        self.assertIsNone(gdf.loc[2, 'geometry'])


class TestAssignOceanBasins(unittest.TestCase):
    """Test ocean basin assignment."""

    def test_assign_basins_raises_without_geopandas(self):
        """Test that assignment raises error without GeoPandas."""
        if geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas is available")

        df = pd.DataFrame({'lat': [23.5], 'lon': [81.2]})

        with self.assertRaises(geographic.GeospatialLibraryError):
            geographic.assign_ocean_basins(df)

    def test_assign_basins_raises_without_data_or_path(self):
        """Test that assignment requires either goas_data or shapefile_path."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({'lat': [23.5], 'lon': [81.2]})

        with self.assertRaises(ValueError) as cm:
            geographic.assign_ocean_basins(df)

        self.assertIn("goas_data or shapefile_path", str(cm.exception))

    @patch('boldgenotyper.geographic.gpd.sjoin')
    @patch('boldgenotyper.geographic.create_points_geodataframe')
    def test_assign_basins_with_valid_coordinates(self, mock_create_points, mock_sjoin):
        """Test basin assignment for samples with valid coordinates."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Input data
        df = pd.DataFrame({
            'sample_id': ['S1', 'S2', 'S3'],
            'lat': [-30.0, -20.0, 40.0],
            'lon': [-60.0, 70.0, -30.0]
        })

        # Mock GOaS data
        mock_goas = create_mock_goas_data()

        # Mock points GeoDataFrame
        import geopandas as gpd
        from shapely.geometry import Point

        mock_points = gpd.GeoDataFrame(
            df.copy(),
            geometry=[Point(-60, -30), Point(70, -20), Point(-30, 40)],
            crs="EPSG:4326"
        )
        mock_create_points.return_value = mock_points

        # Mock spatial join result
        joined_df = mock_points.copy()
        joined_df['basin_name'] = ['South Pacific Ocean', 'Indian Ocean', 'North Atlantic Ocean']
        joined_df['index_right'] = [0, 1, 2]
        mock_sjoin.return_value = joined_df

        # Perform assignment
        result = geographic.assign_ocean_basins(df, goas_data=mock_goas)

        # Check that ocean_basin column was added
        self.assertIn('ocean_basin', result.columns)

        # Check assignments
        self.assertEqual(result.loc[0, 'ocean_basin'], 'South Pacific Ocean')
        self.assertEqual(result.loc[1, 'ocean_basin'], 'Indian Ocean')
        self.assertEqual(result.loc[2, 'ocean_basin'], 'North Atlantic Ocean')

    @patch('boldgenotyper.geographic.gpd.sjoin')
    @patch('boldgenotyper.geographic.create_points_geodataframe')
    def test_assign_basins_handles_samples_outside_basins(self, mock_create_points, mock_sjoin):
        """Test that samples outside all basins are assigned 'Unknown'."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Input data
        df = pd.DataFrame({
            'sample_id': ['S1'],
            'lat': [0.0],
            'lon': [0.0]
        })

        # Mock GOaS data
        mock_goas = create_mock_goas_data()

        # Mock points GeoDataFrame
        import geopandas as gpd
        from shapely.geometry import Point

        mock_points = gpd.GeoDataFrame(
            df.copy(),
            geometry=[Point(0, 0)],
            crs="EPSG:4326"
        )
        mock_create_points.return_value = mock_points

        # Mock spatial join result (no match)
        joined_df = mock_points.copy()
        joined_df['basin_name'] = [np.nan]
        joined_df['index_right'] = [np.nan]
        mock_sjoin.return_value = joined_df

        # Perform assignment
        result = geographic.assign_ocean_basins(df, goas_data=mock_goas)

        # Check that sample is assigned to 'Unknown'
        self.assertEqual(result.loc[0, 'ocean_basin'], 'Unknown')

    @patch('boldgenotyper.geographic.gpd.sjoin')
    @patch('boldgenotyper.geographic.create_points_geodataframe')
    def test_assign_basins_fallback_to_country_ocean(self, mock_create_points, mock_sjoin):
        """Test fallback to country/ocean column for samples without coordinates."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Input data with missing coordinates
        df = pd.DataFrame({
            'sample_id': ['S1', 'S2'],
            'lat': [np.nan, 40.0],
            'lon': [np.nan, -30.0],
            'country/ocean': ['Brazil', 'USA']
        })

        # Mock GOaS data
        mock_goas = create_mock_goas_data()

        # Mock points GeoDataFrame
        import geopandas as gpd
        from shapely.geometry import Point

        mock_points = gpd.GeoDataFrame(
            df.copy(),
            geometry=[None, Point(-30, 40)],
            crs="EPSG:4326"
        )
        mock_create_points.return_value = mock_points

        # Mock spatial join result (only S2 has geometry)
        joined_df = mock_points[mock_points.geometry.notna()].copy()
        joined_df['basin_name'] = ['North Atlantic Ocean']
        joined_df['index_right'] = [0]
        mock_sjoin.return_value = joined_df

        # Perform assignment with fallback
        result = geographic.assign_ocean_basins(
            df,
            goas_data=mock_goas,
            country_ocean_col='country/ocean',
            fallback_to_country_ocean=True
        )

        # S1 should get basin from country/ocean column
        self.assertEqual(result.loc[0, 'ocean_basin'], 'Brazil')

        # S2 should get basin from spatial join
        self.assertEqual(result.loc[1, 'ocean_basin'], 'North Atlantic Ocean')


class TestGetBasinCounts(unittest.TestCase):
    """Test basin counting functionality."""

    def test_get_basin_counts_basic(self):
        """Test basic basin counting."""
        df = pd.DataFrame({
            'ocean_basin': ['North Atlantic Ocean', 'North Atlantic Ocean',
                          'South Pacific Ocean', 'Indian Ocean']
        })

        counts = geographic.get_basin_counts(df)

        self.assertEqual(counts['North Atlantic Ocean'], 2)
        self.assertEqual(counts['South Pacific Ocean'], 1)
        self.assertEqual(counts['Indian Ocean'], 1)

    def test_get_basin_counts_sorting(self):
        """Test that counts are sorted correctly."""
        df = pd.DataFrame({
            'ocean_basin': ['A', 'B', 'B', 'B', 'C', 'C']
        })

        # Sort by count (descending)
        counts = geographic.get_basin_counts(df, sort_by='count')
        keys = list(counts.keys())
        self.assertEqual(keys[0], 'B')  # 3 samples
        self.assertEqual(keys[1], 'C')  # 2 samples
        self.assertEqual(keys[2], 'A')  # 1 sample

        # Sort by name
        counts = geographic.get_basin_counts(df, sort_by='name')
        keys = list(counts.keys())
        self.assertEqual(keys, ['A', 'B', 'C'])

    def test_get_basin_counts_raises_for_missing_column(self):
        """Test that missing basin column raises error."""
        df = pd.DataFrame({'lat': [23.5], 'lon': [81.2]})

        with self.assertRaises(ValueError) as cm:
            geographic.get_basin_counts(df)

        self.assertIn("ocean_basin", str(cm.exception))


class TestGetBasinSummaryStats(unittest.TestCase):
    """Test basin summary statistics."""

    def test_get_basin_summary_stats(self):
        """Test summary statistics generation."""
        df = pd.DataFrame({
            'ocean_basin': ['Basin A', 'Basin A', 'Basin B'],
            'lat': [10.0, 20.0, -30.0],
            'lon': [-40.0, -50.0, 60.0]
        })

        stats = geographic.get_basin_summary_stats(df)

        # Check structure
        self.assertEqual(len(stats), 2)
        self.assertIn('basin', stats.columns)
        self.assertIn('n_samples', stats.columns)
        self.assertIn('lat_min', stats.columns)
        self.assertIn('lat_max', stats.columns)
        self.assertIn('lon_min', stats.columns)

        # Check Basin A stats
        basin_a = stats[stats['basin'] == 'Basin A'].iloc[0]
        self.assertEqual(basin_a['n_samples'], 2)
        self.assertEqual(basin_a['lat_min'], 10.0)
        self.assertEqual(basin_a['lat_max'], 20.0)
        self.assertEqual(basin_a['lat_mean'], 15.0)

    def test_get_basin_summary_stats_handles_missing_coords(self):
        """Test that missing coordinates are handled."""
        df = pd.DataFrame({
            'ocean_basin': ['Basin A', 'Basin B'],
            'lat': [10.0, np.nan],
            'lon': [-40.0, np.nan]
        })

        stats = geographic.get_basin_summary_stats(df)

        # Basin A should have stats
        basin_a = stats[stats['basin'] == 'Basin A'].iloc[0]
        self.assertEqual(basin_a['n_samples'], 1)
        self.assertIsNotNone(basin_a['lat_min'])

        # Basin B should have None for coordinate stats
        basin_b = stats[stats['basin'] == 'Basin B'].iloc[0]
        self.assertEqual(basin_b['n_samples'], 1)
        self.assertTrue(pd.isna(basin_b['lat_min']))


class TestValidateMarineCoordinates(unittest.TestCase):
    """Test single coordinate validation."""

    def test_validate_marine_raises_without_geopandas(self):
        """Test that validation raises error without GeoPandas."""
        if geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas is available")

        with self.assertRaises(geographic.GeospatialLibraryError):
            geographic.validate_marine_coordinates(23.5, 81.2, None)

    def test_validate_marine_raises_for_invalid_ranges(self):
        """Test that invalid coordinate ranges raise ValueError."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        mock_goas = create_mock_goas_data()

        # Invalid latitude
        with self.assertRaises(ValueError) as cm:
            geographic.validate_marine_coordinates(100.0, 0.0, mock_goas)
        self.assertIn("Latitude", str(cm.exception))

        # Invalid longitude
        with self.assertRaises(ValueError) as cm:
            geographic.validate_marine_coordinates(0.0, 200.0, mock_goas)
        self.assertIn("Longitude", str(cm.exception))

    @patch('boldgenotyper.geographic.Point')
    def test_validate_marine_coordinates_marine(self, mock_point_class):
        """Test validation of marine coordinates."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Mock Point
        mock_point = MockPoint(-30.0, 40.0)
        mock_point_class.return_value = mock_point

        # Mock GOaS data that contains this point
        mock_goas = Mock()
        mock_goas.geometry = Mock()
        mock_goas.geometry.contains = Mock(return_value=pd.Series([True]))
        mock_goas.geometry.contains.return_value.any = Mock(return_value=True)

        result = geographic.validate_marine_coordinates(40.0, -30.0, mock_goas)

        self.assertTrue(result)

    @patch('boldgenotyper.geographic.Point')
    def test_validate_marine_coordinates_terrestrial(self, mock_point_class):
        """Test validation of terrestrial coordinates."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Mock Point
        mock_point = MockPoint(0.0, 0.0)
        mock_point_class.return_value = mock_point

        # Mock GOaS data that doesn't contain this point
        mock_goas = Mock()
        mock_goas.geometry = Mock()
        mock_goas.geometry.contains = Mock(return_value=pd.Series([False]))
        mock_goas.geometry.contains.return_value.any = Mock(return_value=False)

        result = geographic.validate_marine_coordinates(0.0, 0.0, mock_goas)

        self.assertFalse(result)


class TestFilterMarineSamples(unittest.TestCase):
    """Test filtering marine vs terrestrial samples."""

    def test_filter_marine_raises_without_geopandas(self):
        """Test that filtering raises error without GeoPandas."""
        if geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas is available")

        df = pd.DataFrame({'lat': [23.5], 'lon': [81.2]})

        with self.assertRaises(geographic.GeospatialLibraryError):
            geographic.filter_marine_samples(df, None)

    @patch('boldgenotyper.geographic.gpd.sjoin')
    @patch('boldgenotyper.geographic.create_points_geodataframe')
    def test_filter_marine_removes_terrestrial(self, mock_create_points, mock_sjoin):
        """Test that terrestrial samples are removed."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Input data
        df = pd.DataFrame({
            'sample_id': ['S1', 'S2', 'S3'],
            'lat': [40.0, 0.0, -30.0],
            'lon': [-30.0, 0.0, -60.0]
        })

        # Mock GOaS data
        mock_goas = create_mock_goas_data()

        # Mock points GeoDataFrame
        import geopandas as gpd
        from shapely.geometry import Point

        mock_points = gpd.GeoDataFrame(
            df.copy(),
            geometry=[Point(-30, 40), Point(0, 0), Point(-60, -30)],
            crs="EPSG:4326"
        )
        mock_create_points.return_value = mock_points

        # Mock spatial join result (S1 and S3 are marine, S2 is terrestrial)
        joined_df = mock_points.copy()
        joined_df['index_right'] = [0, np.nan, 1]  # S2 has NaN
        mock_sjoin.return_value = joined_df

        # Filter with removal
        result = geographic.filter_marine_samples(df, mock_goas, remove_terrestrial=True)

        # Should only have 2 samples (S1 and S3)
        self.assertEqual(len(result), 2)
        self.assertIn('S1', result['sample_id'].values)
        self.assertIn('S3', result['sample_id'].values)
        self.assertNotIn('S2', result['sample_id'].values)

    @patch('boldgenotyper.geographic.gpd.sjoin')
    @patch('boldgenotyper.geographic.create_points_geodataframe')
    def test_filter_marine_flags_terrestrial(self, mock_create_points, mock_sjoin):
        """Test that terrestrial samples are flagged when not removing."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        # Input data
        df = pd.DataFrame({
            'sample_id': ['S1', 'S2'],
            'lat': [40.0, 0.0],
            'lon': [-30.0, 0.0]
        })

        # Mock GOaS data
        mock_goas = create_mock_goas_data()

        # Mock points GeoDataFrame
        import geopandas as gpd
        from shapely.geometry import Point

        mock_points = gpd.GeoDataFrame(
            df.copy(),
            geometry=[Point(-30, 40), Point(0, 0)],
            crs="EPSG:4326"
        )
        mock_create_points.return_value = mock_points

        # Mock spatial join result
        joined_df = mock_points.copy()
        joined_df['index_right'] = [0, np.nan]
        mock_sjoin.return_value = joined_df

        # Filter without removal
        result = geographic.filter_marine_samples(df, mock_goas, remove_terrestrial=False)

        # Should have all samples with is_marine column
        self.assertEqual(len(result), 2)
        self.assertIn('is_marine', result.columns)
        self.assertTrue(result.loc[0, 'is_marine'])
        self.assertFalse(result.loc[1, 'is_marine'])


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""

    @patch('boldgenotyper.geographic.gpd.sjoin')
    @patch('boldgenotyper.geographic.create_points_geodataframe')
    def test_empty_dataframe(self, mock_create_points, mock_sjoin):
        """Test handling of empty DataFrame."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame(columns=['lat', 'lon'])

        # Mock GOaS data
        mock_goas = create_mock_goas_data()

        # Mock empty points GeoDataFrame
        import geopandas as gpd
        mock_points = gpd.GeoDataFrame(df.copy(), geometry=[], crs="EPSG:4326")
        mock_create_points.return_value = mock_points

        # Empty join result
        mock_sjoin.return_value = mock_points

        result = geographic.assign_ocean_basins(df, goas_data=mock_goas)

        # Should return empty DataFrame with ocean_basin column
        self.assertEqual(len(result), 0)
        self.assertIn('ocean_basin', result.columns)

    def test_all_missing_coordinates(self):
        """Test DataFrame with all missing coordinates."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({
            'sample_id': ['S1', 'S2'],
            'lat': [np.nan, np.nan],
            'lon': [np.nan, np.nan],
            'country/ocean': ['Mexico', 'Brazil']
        })

        mock_goas = create_mock_goas_data()

        result = geographic.assign_ocean_basins(
            df,
            goas_data=mock_goas,
            country_ocean_col='country/ocean',
            fallback_to_country_ocean=True
        )

        # Should fall back to country/ocean column
        self.assertEqual(result.loc[0, 'ocean_basin'], 'Mexico')
        self.assertEqual(result.loc[1, 'ocean_basin'], 'Brazil')

    def test_boundary_coordinates(self):
        """Test coordinates at valid boundaries."""
        if not geographic.GEOPANDAS_AVAILABLE:
            self.skipTest("GeoPandas not available")

        df = pd.DataFrame({
            'sample_id': ['S1', 'S2', 'S3', 'S4'],
            'lat': [90.0, -90.0, 0.0, 0.0],
            'lon': [180.0, -180.0, 0.0, 0.0]
        })

        # Should not raise for valid boundary coordinates
        gdf = geographic.create_points_geodataframe(df)

        # All should have valid geometries
        self.assertEqual(gdf.geometry.notna().sum(), 4)


if __name__ == '__main__':
    unittest.main()
