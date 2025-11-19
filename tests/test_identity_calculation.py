"""
Unit tests for improved identity calculation using CIGAR alignment paths.

Tests cover:
- CIGAR string parsing
- Identity calculation with edlib CIGAR paths
- Target-based vs classic identity metrics
- Configuration validation for identity_method
- Edge cases (empty sequences, identical sequences, length differences)
- Method selection in find_best_consensus_match
"""

import unittest
from unittest.mock import patch, Mock
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from boldgenotyper.genotype_assignment import (
    parse_cigar,
    calculate_identity_with_cigar,
    find_best_consensus_match,
    EDLIB_AVAILABLE
)
from boldgenotyper.config import GenotypeAssignmentConfig


class TestCigarParsing(unittest.TestCase):
    """Test CIGAR string parsing."""

    def test_parse_cigar_all_matches(self):
        """Test parsing CIGAR with only matches."""
        result = parse_cigar("10=")
        self.assertEqual(result['='], 10)
        self.assertEqual(result['X'], 0)
        self.assertEqual(result['I'], 0)
        self.assertEqual(result['D'], 0)

    def test_parse_cigar_matches_and_mismatches(self):
        """Test parsing CIGAR with matches and mismatches."""
        result = parse_cigar("5=1X4=")
        self.assertEqual(result['='], 9)
        self.assertEqual(result['X'], 1)
        self.assertEqual(result['I'], 0)
        self.assertEqual(result['D'], 0)

    def test_parse_cigar_with_insertions(self):
        """Test parsing CIGAR with insertions."""
        result = parse_cigar("10=5I")
        self.assertEqual(result['='], 10)
        self.assertEqual(result['X'], 0)
        self.assertEqual(result['I'], 5)
        self.assertEqual(result['D'], 0)

    def test_parse_cigar_with_deletions(self):
        """Test parsing CIGAR with deletions."""
        result = parse_cigar("5=2D3=")
        self.assertEqual(result['='], 8)
        self.assertEqual(result['X'], 0)
        self.assertEqual(result['I'], 0)
        self.assertEqual(result['D'], 2)

    def test_parse_cigar_complex(self):
        """Test parsing complex CIGAR string."""
        result = parse_cigar("5=1X2=3I4=2D5=")
        self.assertEqual(result['='], 16)
        self.assertEqual(result['X'], 1)
        self.assertEqual(result['I'], 3)
        self.assertEqual(result['D'], 2)

    def test_parse_cigar_empty(self):
        """Test parsing empty CIGAR string."""
        result = parse_cigar("")
        self.assertEqual(result['='], 0)
        self.assertEqual(result['X'], 0)
        self.assertEqual(result['I'], 0)
        self.assertEqual(result['D'], 0)


class TestIdentityWithCigar(unittest.TestCase):
    """Test identity calculation using CIGAR alignment paths."""

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_identical_sequences(self):
        """Test identity calculation for identical sequences."""
        result = calculate_identity_with_cigar("ACTGACTG", "ACTGACTG")

        self.assertEqual(result['edit_distance'], 0)
        self.assertEqual(result['matches'], 8)
        self.assertEqual(result['mismatches'], 0)
        self.assertEqual(result['insertions'], 0)
        self.assertEqual(result['deletions'], 0)
        self.assertEqual(result['target_identity'], 1.0)
        self.assertEqual(result['classic_identity'], 1.0)
        self.assertIn(result['method'], ['identical', 'edlib_cigar'])

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_single_mismatch(self):
        """Test identity calculation with one mismatch."""
        result = calculate_identity_with_cigar("ACTGACTG", "ACTGAGTG")

        self.assertEqual(result['edit_distance'], 1)
        self.assertEqual(result['matches'], 7)
        self.assertEqual(result['mismatches'], 1)
        self.assertEqual(result['insertions'], 0)
        self.assertEqual(result['deletions'], 0)
        self.assertAlmostEqual(result['target_identity'], 7/8)
        self.assertAlmostEqual(result['classic_identity'], 7/8)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_sample_with_noisy_end(self):
        """Test identity calculation for sample with noisy 3' end."""
        # Sample has 4 extra noisy bases
        result = calculate_identity_with_cigar("ACTGACTGNNNN", "ACTGACTG")

        self.assertEqual(result['edit_distance'], 4)
        self.assertEqual(result['matches'], 8)
        self.assertEqual(result['mismatches'], 0)
        self.assertEqual(result['insertions'], 4)
        self.assertEqual(result['deletions'], 0)

        # Target-based: 8 matches / 8 target length = 100%
        self.assertEqual(result['target_identity'], 1.0)

        # Classic: 1 - (4 / 12) = 66.7%
        self.assertAlmostEqual(result['classic_identity'], 2/3, places=5)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_partial_sequence(self):
        """Test identity calculation for partial sequence."""
        # Sample is only first half
        result = calculate_identity_with_cigar("ACTG", "ACTGACTG")

        self.assertEqual(result['edit_distance'], 4)
        self.assertEqual(result['matches'], 4)
        self.assertEqual(result['mismatches'], 0)
        self.assertEqual(result['insertions'], 0)
        self.assertEqual(result['deletions'], 4)

        # Target-based: 4 matches / 8 target length = 50%
        self.assertEqual(result['target_identity'], 0.5)

        # Classic: 1 - (4 / 8) = 50%
        self.assertEqual(result['classic_identity'], 0.5)

    def test_empty_sequences(self):
        """Test identity calculation with empty sequences."""
        with self.assertRaises(ValueError):
            calculate_identity_with_cigar("", "")

    def test_one_empty_sequence(self):
        """Test identity calculation when one sequence is empty."""
        result = calculate_identity_with_cigar("ACTG", "")

        self.assertEqual(result['target_identity'], 0.0)
        self.assertEqual(result['classic_identity'], 0.0)
        self.assertEqual(result['method'], 'empty_sequence')

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_realistic_coi_scenario(self):
        """Test realistic COI barcode scenario with noisy ends."""
        # 580bp consensus, sample has 70bp noisy 3' end
        consensus = "A" * 580
        sample = "A" * 580 + "N" * 70

        result = calculate_identity_with_cigar(sample, consensus)

        # Should have 580 matches, 70 insertions
        self.assertEqual(result['matches'], 580)
        self.assertEqual(result['insertions'], 70)

        # Target-based: 580/580 = 100%
        self.assertEqual(result['target_identity'], 1.0)

        # Classic: 1 - (70/650) ≈ 89.2%
        self.assertAlmostEqual(result['classic_identity'], 580/650, places=5)


class TestFindBestConsensusMatch(unittest.TestCase):
    """Test find_best_consensus_match with different identity methods."""

    def setUp(self):
        """Set up test consensus groups."""
        self.consensus_groups = [
            ("consensus_c1", "ACTGACTGACTG"),
            ("consensus_c2", "TGCATGCATGCA"),
        ]

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_perfect_match_target_based(self):
        """Test perfect match with target-based method."""
        result = find_best_consensus_match(
            "ACTGACTGACTG",
            self.consensus_groups,
            min_identity=0.90,
            identity_method="target_based"
        )

        self.assertEqual(result['best_group'], 'consensus_c1')
        self.assertEqual(result['best_identity'], 1.0)
        self.assertEqual(result['identity_method'], 'target_based')
        self.assertIsNotNone(result['best_group'])

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_perfect_match_classic(self):
        """Test perfect match with classic method."""
        result = find_best_consensus_match(
            "ACTGACTGACTG",
            self.consensus_groups,
            min_identity=0.90,
            identity_method="classic"
        )

        self.assertEqual(result['best_group'], 'consensus_c1')
        self.assertEqual(result['best_identity'], 1.0)
        self.assertEqual(result['identity_method'], 'classic')

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_noisy_end_target_based_passes(self):
        """Test that sample with noisy end passes with target-based method."""
        # Sample has perfect match + 4bp noisy end
        result = find_best_consensus_match(
            "ACTGACTGACTGNNNN",
            self.consensus_groups,
            min_identity=0.90,
            identity_method="target_based"
        )

        # Should pass 90% threshold with target-based
        self.assertIsNotNone(result['best_group'])
        self.assertEqual(result['best_group'], 'consensus_c1')
        self.assertEqual(result['best_identity'], 1.0)
        self.assertGreaterEqual(result['best_identity'], 0.90)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_noisy_end_classic_fails(self):
        """Test that sample with noisy end fails with classic method."""
        # Sample has perfect match + 4bp noisy end
        result = find_best_consensus_match(
            "ACTGACTGACTGNNNN",
            self.consensus_groups,
            min_identity=0.90,
            identity_method="classic"
        )

        # Should fail 90% threshold with classic (75% identity)
        self.assertIsNone(result['best_group'])
        self.assertLess(result['best_identity'], 0.90)
        self.assertAlmostEqual(result['classic_identity'], 0.75, places=2)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_both_metrics_always_available(self):
        """Test that both identity metrics are always returned."""
        result = find_best_consensus_match(
            "ACTGACTGACTGNNNN",
            self.consensus_groups,
            min_identity=0.90,
            identity_method="target_based"
        )

        # Both metrics should be present
        self.assertIn('target_identity', result)
        self.assertIn('classic_identity', result)
        self.assertIsNotNone(result['target_identity'])
        self.assertIsNotNone(result['classic_identity'])

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_length_discrepancy_reported(self):
        """Test that length discrepancy is reported."""
        result = find_best_consensus_match(
            "ACTGACTGACTGNNNN",  # 16bp
            self.consensus_groups,  # consensus is 12bp
            min_identity=0.90,
            identity_method="target_based"
        )

        self.assertEqual(result['length_discrepancy'], 4)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_alignment_details_present(self):
        """Test that alignment details are present in result."""
        result = find_best_consensus_match(
            "ACTGACTGACTG",
            self.consensus_groups,
            min_identity=0.90,
            identity_method="target_based"
        )

        self.assertIn('matches', result)
        self.assertIn('mismatches', result)
        self.assertIn('insertions', result)
        self.assertIn('deletions', result)
        self.assertIn('edit_distance', result)
        self.assertIn('cigar', result)

    def test_invalid_identity_method(self):
        """Test that invalid identity method raises error."""
        with self.assertRaises(ValueError) as cm:
            find_best_consensus_match(
                "ACTGACTGACTG",
                self.consensus_groups,
                min_identity=0.90,
                identity_method="invalid"
            )

        self.assertIn("identity_method must be", str(cm.exception))


class TestConfigValidation(unittest.TestCase):
    """Test configuration validation for identity_method."""

    def test_config_target_based(self):
        """Test config with target_based method."""
        config = GenotypeAssignmentConfig(identity_method="target_based")
        self.assertEqual(config.identity_method, "target_based")

    def test_config_classic(self):
        """Test config with classic method."""
        config = GenotypeAssignmentConfig(identity_method="classic")
        self.assertEqual(config.identity_method, "classic")

    def test_config_invalid_method(self):
        """Test that invalid method raises error."""
        with self.assertRaises(ValueError) as cm:
            GenotypeAssignmentConfig(identity_method="invalid")

        self.assertIn("identity_method must be", str(cm.exception))

    def test_config_default_is_target_based(self):
        """Test that default identity_method is target_based."""
        config = GenotypeAssignmentConfig()
        self.assertEqual(config.identity_method, "target_based")


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_very_long_sequences(self):
        """Test identity calculation with very long sequences."""
        # 10kb sequences
        seq1 = "ACTG" * 2500
        seq2 = "ACTG" * 2500

        result = calculate_identity_with_cigar(seq1, seq2)

        self.assertEqual(result['target_identity'], 1.0)
        self.assertEqual(result['classic_identity'], 1.0)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_very_short_sequences(self):
        """Test identity calculation with very short sequences."""
        result = calculate_identity_with_cigar("AC", "AT")

        self.assertEqual(result['matches'], 1)
        self.assertEqual(result['mismatches'], 1)
        self.assertEqual(result['target_identity'], 0.5)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_all_mismatches(self):
        """Test identity calculation when all bases mismatch."""
        result = calculate_identity_with_cigar("AAAA", "TTTT")

        self.assertEqual(result['matches'], 0)
        self.assertEqual(result['mismatches'], 4)
        self.assertEqual(result['target_identity'], 0.0)
        self.assertEqual(result['classic_identity'], 0.0)

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_ambiguous_bases(self):
        """Test identity calculation with ambiguous bases (N, R, Y)."""
        # edlib treats N/R/Y as mismatches
        result = calculate_identity_with_cigar("ACNR", "ACTG")

        # Should have 2 matches (A, C), 2 mismatches (N, R)
        self.assertEqual(result['matches'], 2)
        self.assertEqual(result['mismatches'], 2)


class TestMethodComparison(unittest.TestCase):
    """Test comparison between target-based and classic methods."""

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_methods_equal_for_same_length(self):
        """Test that both methods give same result for same-length sequences."""
        consensus_groups = [("c1", "ACTGACTGACTG")]
        sample = "ACTGAGTGACTG"  # One mismatch

        result_target = find_best_consensus_match(
            sample, consensus_groups, identity_method="target_based"
        )
        result_classic = find_best_consensus_match(
            sample, consensus_groups, identity_method="classic"
        )

        # Should be equal for same-length sequences
        self.assertAlmostEqual(
            result_target['best_identity'],
            result_classic['best_identity'],
            places=5
        )

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_target_more_lenient_for_length_diff(self):
        """Test that target-based is more lenient for length differences."""
        consensus_groups = [("c1", "ACTGACTGACTG")]
        sample = "ACTGACTGACTGNNNN"  # 4bp longer

        result_target = find_best_consensus_match(
            sample, consensus_groups, identity_method="target_based"
        )
        result_classic = find_best_consensus_match(
            sample, consensus_groups, identity_method="classic"
        )

        # Target-based should be more lenient (higher identity)
        self.assertGreater(
            result_target['best_identity'],
            result_classic['best_identity']
        )

    @unittest.skipIf(not EDLIB_AVAILABLE, "edlib not available")
    def test_realistic_impact_70bp_noisy_end(self):
        """Test realistic scenario: 580bp consensus with 70bp noisy 3' end."""
        consensus_groups = [("c1", "A" * 580)]
        sample = "A" * 580 + "N" * 70

        result_target = find_best_consensus_match(
            sample, consensus_groups, min_identity=0.90,
            identity_method="target_based"
        )
        result_classic = find_best_consensus_match(
            sample, consensus_groups, min_identity=0.90,
            identity_method="classic"
        )

        # Target-based should pass (100%)
        self.assertIsNotNone(result_target['best_group'])
        self.assertEqual(result_target['best_identity'], 1.0)

        # Classic should fail (~89.2%)
        self.assertIsNone(result_classic['best_group'])
        self.assertLess(result_classic['best_identity'], 0.90)


if __name__ == '__main__':
    unittest.main()
