"""
Tests for Haplotype Query Module

Tests the haplotype_query module that allows querying new COI sequences
against previously identified haplotypes from a completed BOLD analysis.

Author: Steph Smith (symbioseas@outlook.com)
"""

import pytest
from pathlib import Path
import tempfile
import os

from boldgenotyper import haplotype_query


class TestAlignmentFunctions:
    """Test core alignment and classification functions."""

    def test_align_query_to_haplotype_perfect_match(self):
        """Test alignment with perfect match."""
        query_seq = "ATCGATCGATCG"
        haplotype_seq = "ATCGATCGATCG"

        result = haplotype_query.align_query_to_haplotype(
            query_seq, "query1",
            haplotype_seq, "hap1"
        )

        assert result.query_id == "query1"
        assert result.haplotype_id == "hap1"
        assert result.identity_pct == 100.0
        assert result.differences == 0
        assert result.matches == 12
        assert result.divergence == 0.0

    def test_align_query_to_haplotype_partial_match(self):
        """A partial match: identity below 100% but well above 50%.

        Uses a non-repetitive sequence: on a periodic toy sequence such as
        ATCGATCG... the overlap aligner can legitimately shift the overlap by one
        period (lower coverage, perfect identity) instead of comparing the tail.
        """
        from tests.synthetic_data import A
        nxt = {"A": "C", "C": "G", "G": "T", "T": "A"}
        query_seq = A[:40]
        haplotype_seq = "".join(nxt[c] if i in (5, 15, 25, 39) else c   # 4 substitutions,
                                for i, c in enumerate(A[:40]))           # one at the last base

        result = haplotype_query.align_query_to_haplotype(
            query_seq, "query1",
            haplotype_seq, "hap1"
        )

        assert result.identity_pct < 100.0
        assert result.identity_pct > 50.0
        assert result.n_mismatches == 4           # the overlap alignment counts end mismatches
        assert result.identity_pct == 90.0
        assert result.overlap_coverage == 1.0

    def test_align_query_to_haplotype_length_difference(self):
        """Test alignment with different length sequences."""
        # Query longer than haplotype (common case)
        query_seq = "ATCGATCGATCGATCGATCGATCG"  # 24 bp
        haplotype_seq = "ATCGATCG"  # 8 bp

        result = haplotype_query.align_query_to_haplotype(
            query_seq, "query1",
            haplotype_seq, "hap1"
        )

        # Should find perfect match in the overlapping region
        assert result.identity_pct == 100.0
        assert result.aligned_length == 8


class TestRealData:
    """Test with real Sphyrnidae data."""

    @pytest.mark.skipif(
        not Path("data/Sphyrnidae_test_2").exists(),
        reason="Sphyrnidae test data not available"
    )
    def test_sphyrnidae_query(self):
        """
        Test querying against Sphyrnidae haplotypes.

        This integration test uses real data from the Sphyrnidae_test_2
        analysis to verify the query functionality works end-to-end.
        """
        # Create a test query sequence (use a known haplotype sequence)
        haplotype_file = Path("data/Sphyrnidae_test_2/haplotypes/Sphyrnidae_haplotypes.fasta")

        if not haplotype_file.exists():
            pytest.skip("Haplotype file not found")

        # Read first haplotype to use as query
        from Bio import SeqIO
        haplotypes = list(SeqIO.parse(haplotype_file, "fasta"))

        if not haplotypes:
            pytest.skip("No haplotypes in file")

        # Use first haplotype as query
        query_seq = str(haplotypes[0].seq)
        query_id = haplotypes[0].id

        # Create temporary query file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">{query_id}_query\n")
            f.write(f"{query_seq}\n")
            temp_query = f.name

        try:
            # Run query
            results, metadata, _verdicts = haplotype_query.query_against_haplotypes(
                query_fasta=Path(temp_query),
                haplotype_fasta=haplotype_file,
                analysis_dir=Path("data/Sphyrnidae_test_2"),
                top_n=5
            )

            # Verify results
            assert len(results) > 0, "No results returned"
            best = results[0]

            # Should match itself perfectly
            assert best.identity_pct == 100.0, f"Expected 100% identity, got {best.identity_pct}%"
            assert best.match_quality == "perfect"
            assert best.haplotype_id == query_id

        finally:
            # Cleanup
            if os.path.exists(temp_query):
                os.remove(temp_query)


class TestOutputFormatting:
    """Test output formatting functions."""

    def test_results_to_dataframe(self):
        """Test converting results to DataFrame."""
        # Create mock results
        result1 = haplotype_query.AlignmentResult(
            query_id="query1",
            haplotype_id="hap1",
            score=100.0,
            identity_pct=98.5,
            matches=197,
            aligned_length=200,
            query_length=200,
            haplotype_length=200,
            divergence=0.015,
            aligned_query="A" * 200,
            aligned_haplotype="A" * 200
        )

        result2 = haplotype_query.AlignmentResult(
            query_id="query1",
            haplotype_id="hap2",
            score=95.0,
            identity_pct=95.0,
            matches=190,
            aligned_length=200,
            query_length=200,
            haplotype_length=200,
            divergence=0.05,
            aligned_query="A" * 200,
            aligned_haplotype="A" * 200
        )

        results = [result1, result2]

        # Convert to DataFrame
        df = haplotype_query.results_to_dataframe(results)

        # Verify structure
        assert len(df) == 2
        assert 'query_id' in df.columns
        assert 'haplotype_id' in df.columns
        assert 'rank' in df.columns
        assert 'identity_pct' in df.columns
        assert 'differences' in df.columns
        assert 'match_quality' not in df.columns

        # Verify ranking
        assert df.iloc[0]['rank'] == 1
        assert df.iloc[1]['rank'] == 2

    def test_format_alignment_display(self):
        """Test alignment display formatting."""
        result = haplotype_query.AlignmentResult(
            query_id="query1",
            haplotype_id="hap1",
            score=20.0,
            identity_pct=100.0,
            matches=10,
            aligned_length=10,
            query_length=10,
            haplotype_length=10,
            divergence=0.0,
            aligned_query="ATCGATCGAT",
            aligned_haplotype="ATCGATCGAT"
        )

        formatted = haplotype_query.format_alignment_display(result)

        # Should contain query, haplotype, and match lines
        assert "Query:" in formatted
        assert "Haplotype:" in formatted
        assert "Match:" in formatted
        assert "ATCGATCGAT" in formatted
        assert "||||||||||" in formatted  # All matches


def test_query_validation():
    """Test input validation."""
    with pytest.raises(FileNotFoundError):
        haplotype_query.query_against_haplotypes(
            query_fasta=Path("/nonexistent/file.fasta"),
            haplotype_fasta=Path("/another/nonexistent.fasta")
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
