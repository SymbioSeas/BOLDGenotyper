#!/usr/bin/env python3
"""
Generate Figure 6: Multi-environment comparison summary.

Reads outputs from output/[taxon]_full/ directories produced by
scripts/run_full_datasets.sh and produces a 3-panel comparison figure:

  Panel A — Haplotype yield (haplotypes per input sequence)
  Panel B — Assignment coverage (% samples assigned to a haplotype)
  Panel C — Geographic region count (number of distinct regions with ≥1 sample)

Run from the BOLDGenotyper project root AFTER run_full_datasets.sh completes:
    python scripts/generate_figure6.py

Output:
    benchmarking/Figure_6_MultiEnvironment_Comparison.png
    benchmarking/Figure_6_MultiEnvironment_Comparison.pdf
"""

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Dataset metadata — must match run_full_datasets.sh output directory names
# ---------------------------------------------------------------------------
DATASETS = [
    # (output_subdir, organism_name_in_files, display_label, environment, taxonomic_level)
    ("Panulirus_full",         "Panulirus",         "Panulirus\n(genus)",          "marine",       755),
    ("Sphyrnidae_full",        "Sphyrnidae",        "Sphyrnidae\n(family)",        "marine",      1087),
    ("Carcharhiniformes_full", "Carcharhiniformes", "Carcharhiniformes\n(order)",  "marine",      9244),
    ("Salmonidae_full",        "Salmonidae",        "Salmonidae\n(family)",        "freshwater",  5546),
    ("Pieridae_full",          "Pieridae",          "Pieridae\n(family)",          "terrestrial", 3484),
]

ENV_COLORS = {
    "marine":      "#3182bd",
    "freshwater":  "#31a354",
    "terrestrial": "#a63603",
}

OUTPUT_ROOT = Path("output")
BENCHMARK_DIR = Path("benchmarking")

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def find_file(directory: Path, pattern: str) -> Path | None:
    """Return the first file matching glob pattern, or None."""
    matches = list(directory.glob(pattern))
    return matches[0] if matches else None


def parse_n_input(log_path: Path) -> int | None:
    """Parse number of input sequences from pipeline log."""
    text = log_path.read_text(errors="replace")
    m = re.search(r"Read (\d+) rows", text)
    return int(m.group(1)) if m else None


def parse_n_qc(log_path: Path) -> int | None:
    """Parse number of sequences retained after QC from pipeline log."""
    text = log_path.read_text(errors="replace")
    m = re.search(r"Retained (\d+) samples after QC", text)
    return int(m.group(1)) if m else None


def load_haplotype_stats(out_dir: Path, organism: str) -> pd.DataFrame | None:
    """Load haplotype_stats.csv; try haplotypes/ subdir first then root."""
    for subdir in ["haplotypes", "genotype_assignments", "."]:
        path = find_file(out_dir / subdir, f"*haplotype_stats*")
        if path is None:
            path = find_file(out_dir / subdir, f"*genotype_stats*")
        if path is not None:
            return pd.read_csv(path)
    return None


def load_haplotype_mapping(out_dir: Path, organism: str) -> pd.DataFrame | None:
    """Load haplotype_mapping.csv."""
    for subdir in ["haplotypes", "genotype_assignments", "."]:
        path = find_file(out_dir / subdir, f"*haplotype_mapping*")
        if path is None:
            path = find_file(out_dir / subdir, f"*genotype_mapping*")
        if path is not None:
            return pd.read_csv(path)
    return None


def get_n_regions(out_dir: Path, organism: str) -> int | None:
    """
    Count distinct geographic regions with ≥1 assigned sample.
    Tries (in order):
      1. visualization/json/*distribution_bar_data.json — 'basins' key is list of region names
      2. geographic_analysis/*_region*.csv              — column 'region' or 'basin_name'
    """
    # Try JSON bar data — actual structure: {"basins": [...], "counts": {...}, ...}
    for subdir in ["visualization/json", "visualization", "."]:
        path = find_file(out_dir / subdir, "*distribution_bar_data*")
        if path is not None:
            try:
                with open(path) as f:
                    data = json.load(f)
                if isinstance(data, dict) and "basins" in data:
                    regions = [b for b in data["basins"] if b and b != "Unknown"]
                    if regions:
                        return len(set(regions))
            except Exception:
                pass

    # Try region summary CSV
    for subdir in ["geographic_analysis", "visualization", "."]:
        for pattern in ["*region*summary*", "*region*", "*basin*"]:
            path = find_file(out_dir / subdir, pattern)
            if path is not None:
                try:
                    df = pd.read_csv(path)
                    for col in ("region", "basin_name", "ocean_basin"):
                        if col in df.columns:
                            return df[col].nunique()
                except Exception:
                    pass

    return None


# ---------------------------------------------------------------------------
# Collect metrics for each dataset
# ---------------------------------------------------------------------------

records = []
missing = []

for out_subdir, organism, label, env, n_input_expected in DATASETS:
    out_dir = OUTPUT_ROOT / out_subdir
    rec = {
        "label":     label,
        "env":       env,
        "n_input":   n_input_expected,   # fall back to expected if log unavailable
        "n_haplotypes": None,
        "haplotype_yield": None,
        "coverage_pct":   None,
        "n_regions":      None,
    }

    if not out_dir.exists():
        missing.append(out_subdir)
        records.append(rec)
        continue

    # 1. Pipeline log
    log_path = find_file(out_dir, "*.log")
    if log_path:
        n_input = parse_n_input(log_path)
        if n_input:
            rec["n_input"] = n_input

    # 2. Haplotype stats — count total filtered haplotypes
    stats_df = load_haplotype_stats(out_dir, organism)
    if stats_df is not None:
        # Exclude singletons marked as suspect if column exists
        if "is_suspect" in stats_df.columns:
            stats_df = stats_df[~stats_df["is_suspect"].fillna(False)]
        rec["n_haplotypes"] = len(stats_df)

    # 3. Assignment coverage — haplotype_mapping rows with assigned haplotype
    mapping_df = load_haplotype_mapping(out_dir, organism)
    if mapping_df is not None:
        hap_col = next(
            (c for c in mapping_df.columns if "haplotype_id" in c or "genotype_id" in c),
            None
        )
        if hap_col:
            total = len(mapping_df)
            assigned = mapping_df[hap_col].notna() & (mapping_df[hap_col] != "")
            rec["coverage_pct"] = 100.0 * assigned.sum() / total if total > 0 else None

    # 4. Haplotype yield
    if rec["n_haplotypes"] is not None and rec["n_input"] > 0:
        rec["haplotype_yield"] = rec["n_haplotypes"] / rec["n_input"]

    # 5. Geographic region count
    rec["n_regions"] = get_n_regions(out_dir, organism)

    records.append(rec)

# Report what's missing
if missing:
    print("WARNING: The following output directories were not found.")
    print("Run scripts/run_full_datasets.sh first, then re-run this script.")
    for d in missing:
        print(f"  {d}")
    print()

# Check if we have enough data to plot
available = [r for r in records if any(
    r[k] is not None for k in ["haplotype_yield", "coverage_pct", "n_regions"]
)]
if not available:
    print("ERROR: No data found in output/. Run run_full_datasets.sh first.")
    raise SystemExit(1)

# ---------------------------------------------------------------------------
# Plot — 3 panels
# ---------------------------------------------------------------------------
plt.style.use("seaborn-v0_8-paper")
fig, axes = plt.subplots(1, 3, figsize=(13, 5))

x = np.arange(len(records))
bar_width = 0.55
x_labels = [r["label"] for r in records]

def _plot_panel(ax, values, ylabel, title, letter, fmt="{:.3f}"):
    """Generic bar panel with missing-data hatching."""
    for i, (rec, val) in enumerate(zip(records, values)):
        color = ENV_COLORS[rec["env"]]
        if val is not None:
            bar = ax.bar(i, val, width=bar_width, color=color,
                         edgecolor="white", linewidth=0.8, zorder=3)
            ax.text(i, val + max(v for v in values if v) * 0.02,
                    fmt.format(val), ha="center", va="bottom",
                    fontsize=8.5, fontweight="bold", color="#222222")
        else:
            # Hatch bar at 0 to indicate missing data
            ax.bar(i, 0.001, width=bar_width, color="#cccccc",
                   edgecolor="#999999", linewidth=0.8, hatch="///", zorder=3)
            ax.text(i, 0.005, "pending", ha="center", va="bottom",
                    fontsize=8, color="#888888", style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(f"{letter}. {title}", fontsize=11, fontweight="bold", loc="left")
    ax.grid(axis="y", alpha=0.35, zorder=1)
    ax.set_axisbelow(True)
    ax.set_xlim(-0.6, len(records) - 0.4)

    # N-samples annotation
    for i, rec in enumerate(records):
        ax.text(i, -0.07, f"n={rec['n_input']:,}", ha="center", va="top",
                fontsize=7.5, color="#666666",
                transform=ax.get_xaxis_transform())

# Panel A: Haplotype yield
yields = [r["haplotype_yield"] for r in records]
_plot_panel(axes[0], yields,
            "Haplotype yield\n(haplotypes / input sequences)",
            "Haplotype yield", "A", fmt="{:.4f}")

# Panel B: Assignment coverage
coverages = [r["coverage_pct"] for r in records]
_plot_panel(axes[1], coverages,
            "Assignment coverage (%)",
            "Assignment coverage", "B", fmt="{:.1f}%")
# Add 100% reference line
axes[1].axhline(100, color="#888888", linestyle="--", linewidth=0.8, zorder=2)
axes[1].set_ylim(0, 115)

# Panel C: Geographic region count
regions = [r["n_regions"] for r in records]
_plot_panel(axes[2], regions,
            "Geographic regions\n(regions with ≥1 assigned sample)",
            "Geographic complexity", "C", fmt="{:.0f}")

# Shared environment legend
legend_patches = [
    mpatches.Patch(facecolor=ENV_COLORS[e], edgecolor="white",
                   label=e.capitalize())
    for e in ["marine", "freshwater", "terrestrial"]
]
fig.legend(handles=legend_patches, loc="lower center", ncol=3,
           fontsize=10, framealpha=0.9, edgecolor="#cccccc",
           bbox_to_anchor=(0.5, -0.04))

plt.tight_layout(rect=[0, 0.05, 1, 1])

# Save
BENCHMARK_DIR.mkdir(exist_ok=True)
out_base = BENCHMARK_DIR / "Figure_6_MultiEnvironment_Comparison"
fig.savefig(str(out_base) + ".png", dpi=300, bbox_inches="tight")
fig.savefig(str(out_base) + ".pdf", bbox_inches="tight")
print(f"Saved: {out_base}.png")
print(f"Saved: {out_base}.pdf")

# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
print("\nSummary table:")
print(f"{'Dataset':<26} {'Env':<12} {'N input':>8}  {'Haplotypes':>10}  "
      f"{'Yield':>8}  {'Coverage':>9}  {'Regions':>8}")
print("-" * 92)
for rec in records:
    name = rec["label"].replace("\n", " ")
    yield_str    = f"{rec['haplotype_yield']:.4f}" if rec["haplotype_yield"] is not None else "pending"
    coverage_str = f"{rec['coverage_pct']:.1f}%"   if rec["coverage_pct"]   is not None else "pending"
    haps_str     = str(rec["n_haplotypes"])          if rec["n_haplotypes"]   is not None else "pending"
    regions_str  = str(rec["n_regions"])             if rec["n_regions"]      is not None else "pending"
    print(f"{name:<26} {rec['env']:<12} {rec['n_input']:>8}  {haps_str:>10}  "
          f"{yield_str:>8}  {coverage_str:>9}  {regions_str:>8}")
