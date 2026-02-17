"""
Unit tests for the dereplication module (haplotype-first workflow).

Tests cover:
- MAFFT alignment wrapper (run_mafft_alignment)
"""

import unittest
from unittest.mock import patch, mock_open
import tempfile
from pathlib import Path
import subprocess

from boldgenotyper import dereplication


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


if __name__ == '__main__':
    unittest.main()
