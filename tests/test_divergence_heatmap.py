"""2.1.1: the divergence heatmap is drawn only when it is kept (--extended) and small enough.

It is an extended output; drawn on every run and sweep threshold, it took about 25 minutes for
a family-level Pieridae download (~1,860 haplotypes) and risked exhausting memory.
"""

import logging

import numpy as np
import pandas as pd

from boldgenotyper import divergence_analysis as dv


def _matrix(n):
    names = [f"H{i:04d}" for i in range(n)]
    m = np.random.default_rng(1).uniform(0.001, 0.1, size=(n, n))
    m = (m + m.T) / 2
    np.fill_diagonal(m, 0)
    return pd.DataFrame(m, index=names, columns=names)


def _run(tmp_path, monkeypatch, extended):
    calls = []
    monkeypatch.setattr(dv, "align_consensus_sequences", lambda *a, **k: tmp_path / "aligned.fasta")
    monkeypatch.setattr(dv, "calculate_pairwise_divergence_matrix", lambda *a, **k: _matrix(3))
    monkeypatch.setattr(dv, "generate_divergence_summary", lambda *a, **k: pd.DataFrame())
    monkeypatch.setattr(dv, "create_divergence_heatmap", lambda m, p: calls.append(p))
    taxonomy = tmp_path / "taxonomy.csv"
    pd.DataFrame({"haplotype_id": ["H0000"], "haplotype_sp": ["H0000"]}).to_csv(taxonomy, index=False)
    dv.generate_divergence_analysis(tmp_path / "c.fasta", taxonomy, tmp_path / "out", "Org", extended=extended)
    return calls


def test_divergence_heatmap_only_with_extended(tmp_path, monkeypatch):
    assert _run(tmp_path, monkeypatch, extended=False) == []
    assert len(_run(tmp_path, monkeypatch, extended=True)) == 1


def test_divergence_heatmap_skipped_when_too_large(tmp_path, caplog):
    out = tmp_path / "heatmap.pdf"
    with caplog.at_level(logging.INFO):
        dv.create_divergence_heatmap(_matrix(200), out)       # 40,000 cells
    assert not out.exists() and "skipped" in caplog.text
    dv.create_divergence_heatmap(_matrix(5), out)
    assert out.exists()
