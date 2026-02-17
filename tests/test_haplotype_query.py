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

    def test_classify_match_quality(self):
        """Test match quality classification."""
        assert haplotype_query.classify_match_quality(100.0) == "perfect"
        assert haplotype_query.classify_match_quality(99.8) == "high"
        assert haplotype_query.classify_match_quality(97.5) == "good"
        assert haplotype_query.classify_match_quality(95.5) == "moderate"
        assert haplotype_query.classify_match_quality(90.0) == "low"

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
        assert result.match_quality == "perfect"
        assert result.matches == 12
        assert result.divergence == 0.0

    def test_align_query_to_haplotype_partial_match(self):
        """Test alignment with partial match."""
        query_seq = "ATCGATCGATCG"
        haplotype_seq = "ATCGATCGXXXX"  # 8/12 match

        result = haplotype_query.align_query_to_haplotype(
            query_seq, "query1",
            haplotype_seq, "hap1"
        )

        assert result.identity_pct < 100.0
        assert result.identity_pct > 50.0
        assert result.match_quality in ["high", "good", "moderate", "low"]

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
            results, metadata = haplotype_query.query_against_haplotypes(
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
            match_quality="high",
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
            match_quality="moderate",
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
        assert 'match_quality' in df.columns

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
            match_quality="perfect",
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
