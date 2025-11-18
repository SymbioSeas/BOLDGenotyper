"""
Unit and integration tests for core pipeline orchestration.

Tests cover:
- Pipeline orchestration through all 7 phases
- Input validation and error handling
- Directory structure creation
- Organism name extraction
- Configuration handling
- Result dictionary structure
- Integration with all pipeline modules
- Error recovery and graceful degradation
- Geographic analysis with and without GOaS data
- Phylogenetic tree building (optional)
- File output verification

Test Categories:
- Unit tests: Test helper functions and validation
- Integration tests: Test full pipeline with mocked external tools
- Error handling tests: Test various failure scenarios

Dependencies:
- External tools: MAFFT, trimAl (mocked in tests)
- Test data: Small TSV and FASTA files in test_data/
"""

import sys
import unittest
from unittest.mock import Mock, patch, MagicMock, mock_open, call
import tempfile
import shutil
from pathlib import Path
import pandas as pd
from io import StringIO

# Mock problematic cartopy imports BEFORE importing boldgenotyper modules
# This prevents cartopy/pyproj compatibility issues during testing
cartopy_mock = MagicMock()
sys.modules['cartopy'] = cartopy_mock
sys.modules['cartopy.crs'] = MagicMock()
sys.modules['cartopy.io'] = MagicMock()
sys.modules['cartopy.io.shapereader'] = MagicMock()
sys.modules['cartopy.feature'] = MagicMock()
sys.modules['cartopy.mpl'] = MagicMock()
sys.modules['cartopy.mpl.ticker'] = MagicMock()

from boldgenotyper import core
from boldgenotyper.config import get_default_config


def create_mock_assign_genotypes(tsv_path, n_assigned=3, n_total=3):
    """
    Create a mock assign_genotypes function that creates output files.

    This is a helper function used across multiple test classes to mock
    the assign_genotypes function and create expected output files.

    Parameters
    ----------
    tsv_path : Path
        Path to the input TSV file
    n_assigned : int
        Number of samples successfully assigned
    n_total : int
        Total number of samples

    Returns
    -------
    callable
        Mock function that can be used as side_effect
    """
    def mock_assign_side_effect(*args, **kwargs):
        output_path = kwargs.get('output_path')
        diagnostics_path = kwargs.get('diagnostics_path')

        # Create annotated TSV
        df_assigned = pd.read_csv(tsv_path, sep='\t')
        df_assigned['consensus_group'] = ['consensus_1', 'consensus_1', 'consensus_2'][:n_total]
        df_assigned['identity'] = [0.99, 0.98, 0.97][:n_total]
        df_assigned['status'] = ['assigned'] * n_assigned + ['below_threshold'] * (n_total - n_assigned)
        df_assigned.to_csv(output_path, sep='\t', index=False)

        # Create diagnostics CSV
        diagnostics = pd.DataFrame({
            'processid': df_assigned['processid'],
            'consensus_group': df_assigned['consensus_group'],
            'identity': df_assigned['identity'],
            'status': df_assigned['status']
        })
        diagnostics.to_csv(diagnostics_path, index=False)

        return {
            'total': n_total,
            'assigned': n_assigned,
            'below_threshold': n_total - n_assigned,
            'ties': 0,
            'no_sequence': 0
        }

    return mock_assign_side_effect


class TestHelperFunctions(unittest.TestCase):
    """Test helper functions for directory setup and organism extraction."""

    def test_extract_organism_name_genus_species(self):
        """Test extraction from Genus_species.tsv format."""
        path = Path("Euprymna_scolopes.tsv")
        result = core._extract_organism_name(path)
        self.assertEqual(result, "Euprymna_scolopes")

    def test_extract_organism_name_single_genus(self):
        """Test extraction from Genus.tsv format."""
        path = Path("Euprymna.tsv")
        result = core._extract_organism_name(path)
        self.assertEqual(result, "Euprymna")

    def test_extract_organism_name_with_path(self):
        """Test extraction from file with full path."""
        path = Path("data/specimens/Sphyrna_lewini.tsv")
        result = core._extract_organism_name(path)
        self.assertEqual(result, "Sphyrna_lewini")

    def test_extract_organism_name_ignores_data_suffix(self):
        """Test that 'data' suffix is ignored."""
        path = Path("Euprymna_data.tsv")
        result = core._extract_organism_name(path)
        self.assertEqual(result, "Euprymna")

    def test_extract_organism_name_ignores_bold_suffix(self):
        """Test that 'bold' suffix is ignored."""
        path = Path("Sphyrna_bold.tsv")
        result = core._extract_organism_name(path)
        self.assertEqual(result, "Sphyrna")

    def test_extract_organism_name_capitalizes(self):
        """Test that single names are capitalized."""
        path = Path("euprymna.tsv")
        result = core._extract_organism_name(path)
        self.assertEqual(result, "Euprymna")

    def test_setup_directories_creates_structure(self):
        """Test that directory structure is created correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "test_output"
            dirs = core._setup_directories(base_path)

            # Check all expected directories exist
            expected_dirs = [
                'base', 'intermediate', 'dereplication', 'consensus',
                'assignments', 'taxonomy', 'geographic', 'phylogenetic', 'reports'
            ]
            for dir_name in expected_dirs:
                self.assertIn(dir_name, dirs)
                self.assertTrue(dirs[dir_name].exists())
                self.assertTrue(dirs[dir_name].is_dir())

            # Check hierarchy
            self.assertEqual(dirs['base'], base_path)
            self.assertTrue(dirs['dereplication'].parent == dirs['intermediate'])

    def test_setup_directories_existing_dirs(self):
        """Test that setup handles existing directories gracefully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "existing"
            base_path.mkdir(parents=True)

            # Create some directories manually
            (base_path / "intermediate").mkdir()

            # Should not raise error
            dirs = core._setup_directories(base_path)
            self.assertTrue(dirs['intermediate'].exists())


class TestInputValidation(unittest.TestCase):
    """Test input validation and error handling."""

    def test_run_pipeline_missing_input_file(self):
        """Test that missing input file raises FileNotFoundError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(FileNotFoundError):
                core.run_pipeline(
                    input_tsv="nonexistent.tsv",
                    output_dir=tmpdir
                )

    def test_run_pipeline_invalid_cluster_threshold_too_low(self):
        """Test that cluster_threshold <= 0 raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "test.tsv"
            tsv_path.write_text("processid\tnuc\n")

            with self.assertRaises(ValueError) as cm:
                core.run_pipeline(
                    input_tsv=str(tsv_path),
                    output_dir=tmpdir,
                    cluster_threshold=0.0
                )
            self.assertIn("cluster_threshold", str(cm.exception))

    def test_run_pipeline_invalid_cluster_threshold_too_high(self):
        """Test that cluster_threshold >= 1 raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "test.tsv"
            tsv_path.write_text("processid\tnuc\n")

            with self.assertRaises(ValueError) as cm:
                core.run_pipeline(
                    input_tsv=str(tsv_path),
                    output_dir=tmpdir,
                    cluster_threshold=1.5
                )
            self.assertIn("cluster_threshold", str(cm.exception))

    def test_run_pipeline_invalid_min_identity_too_low(self):
        """Test that min_identity <= 0 raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "test.tsv"
            tsv_path.write_text("processid\tnuc\n")

            with self.assertRaises(ValueError) as cm:
                core.run_pipeline(
                    input_tsv=str(tsv_path),
                    output_dir=tmpdir,
                    min_identity=0.0
                )
            self.assertIn("min_identity", str(cm.exception))

    def test_run_pipeline_invalid_threads(self):
        """Test that threads < 1 raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "test.tsv"
            tsv_path.write_text("processid\tnuc\n")

            with self.assertRaises(ValueError) as cm:
                core.run_pipeline(
                    input_tsv=str(tsv_path),
                    output_dir=tmpdir,
                    threads=0
                )
            self.assertIn("threads", str(cm.exception))


class TestPipelineIntegration(unittest.TestCase):
    """Test full pipeline integration with mocked external tools."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()
        self.tsv_path = Path(self.tmpdir) / "Euprymna_test.tsv"
        self.output_dir = Path(self.tmpdir) / "output"

        # Create minimal test TSV
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003'],
            'nuc': ['ACGTACGTACGT' * 50, 'ACGTACGTACGT' * 50, 'ACGTACGTACGT' * 50],
            'species': ['Euprymna scolopes', 'Euprymna scolopes', 'Euprymna morsei'],
            'genus': ['Euprymna', 'Euprymna', 'Euprymna'],
            'lat': [25.0, 26.0, 35.0],
            'lon': [130.0, 131.0, 140.0],
            'coord': ['25.0, 130.0', '26.0, 131.0', '35.0, 140.0'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(self.tsv_path, sep='\t', index=False)

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.visualization.plot_ocean_basin_abundance')
    @patch('boldgenotyper.visualization.plot_distribution_map')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.geographic.assign_ocean_basins')
    @patch('boldgenotyper.geographic.load_goas_data')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_minimal_success(
        self, mock_check_tool, mock_load_goas, mock_assign_basins,
        mock_dereplicate, mock_assign_genotypes, mock_plot_map,
        mock_plot_bar, mock_plot_identity, mock_generate_summary
    ):
        """Test successful pipeline run with minimal configuration."""
        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock GOaS data
        mock_load_goas.return_value = Mock()

        # Mock basin assignment
        mock_df = pd.read_csv(self.tsv_path, sep='\t')
        mock_df['ocean_basin'] = 'North Pacific'
        mock_assign_basins.return_value = mock_df

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
            'consensus_2': SeqRecord(Seq('ACGT' * 150), id='consensus_2'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment - create expected output files
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(self.tsv_path, n_assigned=3, n_total=3)

        # Run pipeline
        results = core.run_pipeline(
            input_tsv=str(self.tsv_path),
            output_dir=str(self.output_dir),
            cluster_threshold=0.01,
            min_identity=0.90,
            threads=1,
            skip_geographic=False
        )

        # Verify results structure
        self.assertTrue(results['success'])
        self.assertEqual(results['organism'], 'Euprymna_test')
        self.assertEqual(results['n_samples_input'], 3)
        self.assertGreater(results['n_genotypes'], 0)
        self.assertIsInstance(results['files'], dict)
        self.assertIsInstance(results['errors'], list)

        # Verify key files were created
        self.assertIn('annotated_data', results['files'])
        self.assertTrue(results['files']['annotated_data'].exists())

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_skip_geographic(
        self, mock_check_tool, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test pipeline with geographic analysis skipped."""
        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(self.tsv_path, n_assigned=2, n_total=3)

        # Run pipeline with geographic analysis skipped
        results = core.run_pipeline(
            input_tsv=str(self.tsv_path),
            output_dir=str(self.output_dir),
            skip_geographic=True
        )

        # Verify geographic analysis was not performed
        self.assertFalse(results['geographic_analysis'])
        self.assertTrue(results['success'])

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.geographic.load_goas_data')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_missing_goas_shapefile(
        self, mock_check_tool, mock_load_goas, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test pipeline gracefully handles missing GOaS shapefile."""
        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock GOaS data loading to fail
        mock_load_goas.side_effect = FileNotFoundError("GOaS shapefile not found")

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(self.tsv_path, n_assigned=3, n_total=3)

        # Create custom config with non-existent GOaS path
        cfg = get_default_config()
        cfg = cfg.update(
            geographic__goas_shapefile_path=Path("/nonexistent/goas.shp")
        )

        # Run pipeline - should not fail, just skip geographic analysis
        results = core.run_pipeline(
            input_tsv=str(self.tsv_path),
            output_dir=str(self.output_dir),
            config_obj=cfg
        )

        # Should still succeed but without geographic analysis
        self.assertTrue(results['success'])
        self.assertFalse(results['geographic_analysis'])
        self.assertIn("GOaS shapefile not found", results['errors'][0])

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_phylogenetic_tree')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.phylogenetics.relabel_tree_and_alignment')
    @patch('boldgenotyper.phylogenetics.build_phylogeny')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_with_phylogeny(
        self, mock_check_tool, mock_dereplicate, mock_assign_genotypes,
        mock_build_tree, mock_relabel, mock_plot_identity,
        mock_plot_tree, mock_generate_summary
    ):
        """Test pipeline with phylogenetic tree building enabled."""
        # Mock external tool checks (MAFFT and FastTree available)
        def mock_tool_check(tool_name):
            return tool_name in ['mafft', 'fasttree', 'trimal']
        mock_check_tool.side_effect = mock_tool_check

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
            'consensus_2': SeqRecord(Seq('ACGT' * 150), id='consensus_2'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(self.tsv_path, n_assigned=3, n_total=3)

        # Mock phylogeny building
        mock_tree = Mock()
        mock_build_tree.return_value = mock_tree

        # Create dummy tree file
        tree_file = self.output_dir / "phylogenetic" / "Euprymna_test_tree.nwk"
        tree_file.parent.mkdir(parents=True, exist_ok=True)
        tree_file.write_text("(consensus_1:0.1,consensus_2:0.1);")

        # Create dummy alignment file for relabeling
        aln_file = self.output_dir / "phylogenetic" / "Euprymna_test_aligned.fasta"
        aln_file.write_text(">consensus_1\nACGT\n>consensus_2\nACGT\n")

        # Run pipeline with phylogeny enabled
        results = core.run_pipeline(
            input_tsv=str(self.tsv_path),
            output_dir=str(self.output_dir),
            build_phylogeny=True,
            skip_geographic=True
        )

        # Verify phylogeny was built
        self.assertTrue(results['success'])
        self.assertIsNotNone(results['phylogenetic_tree'])
        self.assertIn('phylogenetic_tree', results['files'])

        # Verify tree building was called
        mock_build_tree.assert_called_once()

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_custom_organism_name(
        self, mock_check_tool, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test pipeline with custom organism name."""
        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(self.tsv_path, n_assigned=3, n_total=3)

        # Run pipeline with custom organism name
        custom_name = "MyCustomOrganism"
        results = core.run_pipeline(
            input_tsv=str(self.tsv_path),
            output_dir=str(self.output_dir),
            organism_name=custom_name,
            skip_geographic=True
        )

        # Verify custom organism name was used
        self.assertEqual(results['organism'], custom_name)
        self.assertTrue(results['success'])

        # Verify files use custom name
        annotated_file = results['files']['annotated_data']
        self.assertIn(custom_name, annotated_file.name)

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_custom_config(
        self, mock_check_tool, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test pipeline with custom configuration object."""
        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(self.tsv_path, n_assigned=3, n_total=3)

        # Create custom config
        custom_cfg = get_default_config()
        custom_cfg = custom_cfg.update(
            dereplication__clustering_threshold=0.005,
            genotype_assignment__min_identity=0.95
        )

        # Run pipeline with custom config
        results = core.run_pipeline(
            input_tsv=str(self.tsv_path),
            output_dir=str(self.output_dir),
            config_obj=custom_cfg,
            skip_geographic=True
        )

        # Verify success
        self.assertTrue(results['success'])

        # Verify dereplicate was called with custom threshold
        call_args = mock_dereplicate.call_args
        self.assertEqual(call_args[1]['threshold'], 0.005)

        # Verify assign_genotypes was called with custom threshold
        call_args = mock_assign_genotypes.call_args
        self.assertEqual(call_args[1]['min_identity'], 0.95)


class TestErrorHandling(unittest.TestCase):
    """Test error handling and recovery."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    @patch('boldgenotyper.metadata.parse_bold_tsv')
    def test_run_pipeline_phase1_failure(self, mock_parse_tsv):
        """Test pipeline handles Phase 1 failure gracefully."""
        # Create test TSV
        tsv_path = Path(self.tmpdir) / "test.tsv"
        tsv_path.write_text("processid\tnuc\nSAMPLE001\tACGT\n")

        # Mock parse_bold_tsv to raise error
        mock_parse_tsv.side_effect = ValueError("Invalid TSV format")

        # Should raise the exception
        with self.assertRaises(ValueError):
            core.run_pipeline(
                input_tsv=str(tsv_path),
                output_dir=str(Path(self.tmpdir) / "output")
            )

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_visualization_failure_non_critical(
        self, mock_check_tool, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test that visualization errors don't fail the pipeline."""
        # Create test TSV
        tsv_path = Path(self.tmpdir) / "test.tsv"
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002'],
            'nuc': ['ACGTACGT' * 50, 'ACGTACGT' * 50],
            'species': ['Test sp', 'Test sp'],
            'genus': ['Test', 'Test'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(tsv_path, n_assigned=2, n_total=2)

        # Mock visualization to fail
        mock_plot_identity.side_effect = Exception("Plotting error")

        # Run pipeline - should succeed despite visualization error
        results = core.run_pipeline(
            input_tsv=str(tsv_path),
            output_dir=str(Path(self.tmpdir) / "output"),
            skip_geographic=True
        )

        # Pipeline should still succeed
        self.assertTrue(results['success'])
        # But should have error recorded
        self.assertTrue(len(results['errors']) > 0)

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_run_pipeline_report_generation_failure_non_critical(
        self, mock_check_tool, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test that report generation errors don't fail the pipeline."""
        # Create test TSV
        tsv_path = Path(self.tmpdir) / "test.tsv"
        test_data = {
            'processid': ['SAMPLE001'],
            'nuc': ['ACGTACGT' * 50],
            'species': ['Test sp'],
            'genus': ['Test'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(tsv_path, n_assigned=1, n_total=1)

        # Mock report generation to fail
        mock_generate_summary.side_effect = Exception("Report error")

        # Run pipeline - should succeed despite report error
        results = core.run_pipeline(
            input_tsv=str(tsv_path),
            output_dir=str(Path(self.tmpdir) / "output"),
            skip_geographic=True
        )

        # Pipeline should still succeed
        self.assertTrue(results['success'])
        # But should have error recorded
        self.assertTrue(len(results['errors']) > 0)
        self.assertTrue(any('Report generation failed' in err for err in results['errors']))


class TestResultsDictionary(unittest.TestCase):
    """Test the structure and contents of the results dictionary."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.tmpdir)

    @patch('boldgenotyper.reports.generate_assignment_summary')
    @patch('boldgenotyper.visualization.plot_identity_distribution')
    @patch('boldgenotyper.genotype_assignment.assign_genotypes')
    @patch('boldgenotyper.dereplication.dereplicate_from_fasta')
    @patch('boldgenotyper.utils.check_external_tool')
    def test_results_dict_structure(
        self, mock_check_tool, mock_dereplicate,
        mock_assign_genotypes, mock_plot_identity, mock_generate_summary
    ):
        """Test that results dictionary has expected structure."""
        # Create test TSV
        tsv_path = Path(self.tmpdir) / "test.tsv"
        test_data = {
            'processid': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003'],
            'nuc': ['ACGT' * 150, 'ACGT' * 150, 'ACGT' * 150],
            'species': ['Test sp', 'Test sp', 'Test sp'],
            'genus': ['Test', 'Test', 'Test'],
        }
        df = pd.DataFrame(test_data)
        df.to_csv(tsv_path, sep='\t', index=False)

        # Mock external tool checks
        mock_check_tool.return_value = True

        # Mock dereplication
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        mock_consensus = {
            'consensus_1': SeqRecord(Seq('ACGT' * 150), id='consensus_1'),
            'consensus_2': SeqRecord(Seq('ACGT' * 150), id='consensus_2'),
        }
        mock_dereplicate.return_value = mock_consensus

        # Mock genotype assignment
        mock_assign_genotypes.side_effect = create_mock_assign_genotypes(tsv_path, n_assigned=2, n_total=3)

        # Run pipeline
        results = core.run_pipeline(
            input_tsv=str(tsv_path),
            output_dir=str(Path(self.tmpdir) / "output"),
            skip_geographic=True
        )

        # Check required keys
        required_keys = [
            'success', 'output_dir', 'organism', 'n_samples_input',
            'n_samples_filtered', 'n_genotypes', 'n_assigned',
            'assignment_rate', 'files', 'geographic_analysis',
            'phylogenetic_tree', 'errors'
        ]

        for key in required_keys:
            self.assertIn(key, results, f"Missing key: {key}")

        # Check data types
        self.assertIsInstance(results['success'], bool)
        self.assertIsInstance(results['output_dir'], Path)
        self.assertIsInstance(results['organism'], str)
        self.assertIsInstance(results['n_samples_input'], int)
        self.assertIsInstance(results['n_samples_filtered'], int)
        self.assertIsInstance(results['n_genotypes'], int)
        self.assertIsInstance(results['n_assigned'], int)
        self.assertIsInstance(results['assignment_rate'], float)
        self.assertIsInstance(results['files'], dict)
        self.assertIsInstance(results['geographic_analysis'], bool)
        self.assertIsInstance(results['errors'], list)

        # Check value ranges
        self.assertTrue(results['success'])
        self.assertGreaterEqual(results['n_samples_input'], 0)
        self.assertGreaterEqual(results['n_samples_filtered'], 0)
        self.assertGreaterEqual(results['n_genotypes'], 0)
        self.assertGreaterEqual(results['n_assigned'], 0)
        self.assertGreaterEqual(results['assignment_rate'], 0.0)
        self.assertLessEqual(results['assignment_rate'], 1.0)


if __name__ == '__main__':
    unittest.main()
