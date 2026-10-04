"""The run report (report.html), built from the final run folder.

generate_html_report reads the run folder's tables, figure manifests and
provenance (report.data), and writes one self-contained page whose figures
are linked from the folder: Overview (with key findings), Records and
haplotypes, Taxonomy and labels, Phylogeny, Geography, Sampling and
metadata, Divergence (--extended), Methods, Output files.
"""

from pathlib import Path
from typing import Optional

from . import (data, divergence, files, geography, methods, overview, phylogeny, records, sampling,
               taxonomy)
from .page import render_page


def generate_html_report(run_dir, organism: str, provenance: Optional[dict] = None,
                         version: Optional[str] = None) -> Path:
    """Write run_dir/report.html and return its path."""
    from .. import __version__

    run = data.load_run(run_dir, organism, provenance)
    version = version or run.provenance.get("boldgenotyper_version") or __version__
    sections = [overview.build(run), records.build(run), taxonomy.build(run), phylogeny.build(run),
                geography.build(run), sampling.build(run)]
    if divergence.applies(run):
        sections.append(divergence.build(run))
    sections.append(methods.build(run, version))
    sections.append(files.build(run))
    started = run.provenance.get("started_at") or ""
    subtitle = (f"BOLDGenotyper v{version}" + (f" | run started {started}" if started else "")
                + " | open this report from the run folder: figures are linked, not embedded")
    path = Path(run_dir) / "report.html"
    path.write_text(render_page(f"{organism}: BOLDGenotyper report", subtitle, sections), encoding="utf-8")
    return path
