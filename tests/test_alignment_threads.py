"""--threads reaches the MAFFT alignments of the records and of the haplotypes (2.3.1).

Up to 2.3.0 the pipeline's --threads reached only the tree; the record alignment (half of a
large run) and the divergence alignment ran on one thread. MAFFT gives the same alignment for
any number of threads, so the results must not change.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from boldgenotyper import dereplication, divergence_analysis

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")


class _Stop(Exception):
    pass


def _record_commands(monkeypatch):
    """subprocess.run that records each command and writes an empty alignment."""
    commands = []

    def fake_run(cmd, stdout=None, **kwargs):
        commands.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
    return commands


def test_record_alignment_passes_threads_to_mafft(tmp_path, monkeypatch):
    commands = _record_commands(monkeypatch)
    dereplication.run_mafft_alignment("in.fasta", str(tmp_path / "out.fasta"), threads=8)
    dereplication.run_mafft_alignment("in.fasta", str(tmp_path / "out.fasta"))
    assert commands == [["mafft", "--auto", "--thread", "8", "in.fasta"], ["mafft", "--auto", "in.fasta"]]


def test_record_alignment_refuses_threads_below_one(tmp_path):
    with pytest.raises(ValueError, match="threads"):
        dereplication.run_mafft_alignment("in.fasta", str(tmp_path / "out.fasta"), threads=0)


def test_divergence_alignment_passes_threads_to_mafft(tmp_path, monkeypatch):
    commands = _record_commands(monkeypatch)
    divergence_analysis.align_consensus_sequences(Path("haps.fasta"), tmp_path, "Synth", threads=8)
    divergence_analysis.align_consensus_sequences(Path("haps.fasta"), tmp_path, "Synth")
    assert commands == [["mafft", "--auto", "--quiet", "--thread", "8", "haps.fasta"],
                        ["mafft", "--auto", "--quiet", "haps.fasta"]]


def test_identify_haplotypes_hands_its_threads_to_the_alignment(tmp_path, monkeypatch):
    seen = []

    def fake_alignment(input_fasta, output_fasta, mafft_options=None, threads=1):
        seen.append(threads)
        raise _Stop

    monkeypatch.setattr(dereplication, "run_mafft_alignment", fake_alignment)
    with pytest.raises(_Stop):
        dereplication.identify_haplotypes("x.tsv", "x.fasta", tmp_path, threads=8)
    assert seen == [8]


def test_divergence_analysis_hands_its_threads_to_the_alignment(tmp_path, monkeypatch):
    seen = []

    def fake_alignment(consensus_fasta, output_dir, organism, threads=1):
        seen.append(threads)
        raise _Stop

    monkeypatch.setattr(divergence_analysis, "align_consensus_sequences", fake_alignment)
    (tmp_path / "tax.csv").write_text("haplotype_id\nH0001\n")      # read before the alignment
    with pytest.raises(_Stop):
        divergence_analysis.generate_divergence_analysis(tmp_path / "haps.fasta", tmp_path / "tax.csv",
                                                         tmp_path, "Synth", threads=8)
    assert seen == [8]


@needs_mafft
def test_pipeline_threads_reach_both_alignments_and_change_no_result(tmp_path, monkeypatch):
    from tests.test_pipeline_synthetic import run_synthetic

    threads_seen = {"records": [], "divergence": []}
    real_records, real_divergence = dereplication.run_mafft_alignment, divergence_analysis.align_consensus_sequences

    def records(*args, threads=1, **kwargs):
        threads_seen["records"].append(threads)
        return real_records(*args, threads=threads, **kwargs)

    def divergence(*args, threads=1, **kwargs):
        threads_seen["divergence"].append(threads)
        return real_divergence(*args, threads=threads, **kwargs)

    ok_one, one = run_synthetic(tmp_path / "one", n_threads=1)
    monkeypatch.setattr(dereplication, "run_mafft_alignment", records)
    monkeypatch.setattr(divergence_analysis, "align_consensus_sequences", divergence)
    ok_two, two = run_synthetic(tmp_path / "two", n_threads=2)
    assert ok_one and ok_two
    assert threads_seen == {"records": [2], "divergence": [2]}
    results = sorted(p.relative_to(one) for p in one.rglob("*")
                     if p.is_file() and p.suffix in {".csv", ".fasta", ".tsv"})
    assert results, "no result files to compare"
    for relative in results:
        assert (two / relative).read_bytes() == (one / relative).read_bytes(), str(relative)
