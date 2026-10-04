"""Everything the report reads, loaded once from the final run folder."""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

FIGURE_FOLDERS = ["phylogeny", "geography", "extended/geography", "extended/metadata/figures"]


def _csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, low_memory=False, **kwargs)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError):
        return pd.DataFrame()


@dataclass
class RunData:
    run_dir: Path
    organism: str
    provenance: dict
    summary: dict
    records: pd.DataFrame
    haplotypes: pd.DataFrame
    label_conflicts: pd.DataFrame
    coverage: pd.DataFrame
    excluded_tips: pd.DataFrame
    manifests: Dict[str, dict]
    region_col: Optional[str]
    _non_monophyletic: Optional[List[str]] = field(default=None, repr=False)

    @property
    def applied(self) -> dict:
        return self.provenance.get("applied", {}) or {}

    @property
    def config(self) -> dict:
        return self.provenance.get("config", {}) or {}

    @property
    def extended(self) -> bool:
        return bool(self.config.get("extended"))

    def file(self, rel: str) -> Path:
        return self.run_dir / rel

    def non_monophyletic(self) -> List[str]:
        """Taxa whose haplotypes do not form a clade in phylogeny/tree.nwk."""
        if self._non_monophyletic is None:
            self._non_monophyletic = []
            tree_path = self.file("phylogeny/tree.nwk")
            if tree_path.exists() and not self.haplotypes.empty:
                from Bio import Phylo
                from ..tree_figures import _haplotype_info, collapse_groups
                tree = Phylo.read(str(tree_path), "newick")
                _, taxon_of, records_of = _haplotype_info(self.haplotypes)
                self._non_monophyletic = collapse_groups(tree, taxon_of, records_of)[1]
        return self._non_monophyletic


def load_run(run_dir, organism: str, provenance: Optional[dict] = None) -> RunData:
    """Load a run folder for the report; provenance defaults to run_provenance.json."""
    run_dir = Path(run_dir)
    if provenance is None:
        path = run_dir / "run_provenance.json"
        provenance = json.loads(path.read_text()) if path.exists() else {}
    summary_df = _csv(run_dir / "summary.csv")
    summary = (summary_df.iloc[0].to_dict() if len(summary_df)
               else dict((provenance.get("applied") or {}).get("assignment_summary") or {}))
    records = _csv(run_dir / "records.csv")
    manifests = {}
    for folder in FIGURE_FOLDERS:
        path = run_dir / folder / "manifest.json"
        if path.exists():
            try:
                manifests[folder] = json.loads(path.read_text())
            except ValueError:
                pass
    applied = provenance.get("applied") or {}
    col = applied.get("geo_category")
    region_col = col if applied.get("geographic_analysis") and col in records.columns else None
    return RunData(
        run_dir=run_dir, organism=organism, provenance=provenance, summary=summary,
        records=records, haplotypes=_csv(run_dir / "haplotypes.csv"),
        label_conflicts=_csv(run_dir / "label_conflicts.csv"),
        coverage=_csv(run_dir / "metadata_coverage.csv", keep_default_na=False),
        excluded_tips=_csv(run_dir / "phylogeny" / "excluded_tips.txt", sep="\t"),
        manifests=manifests, region_col=region_col,
    )
