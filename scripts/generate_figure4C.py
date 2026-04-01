#!/usr/bin/env python3
"""
Generate Figure 4C: Cross-dataset parameter sweep summary.

Reads recommendations.txt files from benchmarking/parameter_sweep/[taxon]/
and produces a bar chart showing the optimal singleton distance threshold
per dataset, colour-coded by ecological environment.

Run from the BOLDGenotyper project root:
    python scripts/generate_figure4C.py

Output:
    benchmarking/parameter_sweep/Figure_4C_Threshold_Summary.png
    benchmarking/parameter_sweep/Figure_4C_Threshold_Summary.pdf
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Configuration — dataset metadata (order determines x-axis layout)
# ---------------------------------------------------------------------------
DATASETS = [
    # (directory_name, display_label, environment, taxonomic_level)
    ("Panulirus",        "Panulirus\n(genus)",         "marine",      755),
    ("Sphyrnidae",       "Sphyrnidae\n(family)",       "marine",    1_087),
    ("Carcharhiniformes","Carcharhiniformes\n(order)",  "marine",    9_244),
    ("Salmonidae",       "Salmonidae\n(family)",        "freshwater", 5_546),
    ("Pieridae",         "Pieridae\n(family)",          "terrestrial",3_484),
]

ENV_COLORS = {
    "marine":      "#3182bd",   # blue
    "freshwater":  "#31a354",   # green
    "terrestrial": "#a63603",   # brown-orange
}

SWEEP_DIR = Path("benchmarking/parameter_sweep")

# ---------------------------------------------------------------------------
# Parse a recommendations.txt file for threshold and acceptable range
# ---------------------------------------------------------------------------
def parse_recommendations(taxon: str) -> dict:
    """Extract elbow threshold and acceptable range from recommendations.txt."""
    path = SWEEP_DIR / taxon / "recommendations.txt"
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")
    text = path.read_text()

    # Primary recommendation:  "PRIMARY:  min_singleton_distance = 0.015"
    m = re.search(r"PRIMARY:\s+min_singleton_distance\s*=\s*([\d.]+)", text)
    if not m:
        raise ValueError(f"Could not parse PRIMARY threshold from {path}")
    threshold = float(m.group(1))

    # Acceptable range:  "ACCEPTABLE RANGE:  0.010 – 0.080"
    m = re.search(r"ACCEPTABLE RANGE:\s+([\d.]+)\s*[–-]\s*([\d.]+)", text)
    range_lo = float(m.group(1)) if m else threshold
    range_hi = float(m.group(2)) if m else threshold

    # Assignment coverage at elbow
    m = re.search(r"Assignment coverage at elbow:\s+([\d.]+)%", text)
    coverage = float(m.group(1)) if m else None

    # Haplotype count at elbow
    m = re.search(r"At [\d.]+ \(elbow\):\s+([\d.]+) haplotypes", text)
    n_haplotypes = int(float(m.group(1))) if m else None

    return {
        "threshold": threshold,
        "range_lo": range_lo,
        "range_hi": range_hi,
        "coverage": coverage,
        "n_haplotypes": n_haplotypes,
    }

# ---------------------------------------------------------------------------
# Collect data
# ---------------------------------------------------------------------------
records = []
for dirname, label, env, n_samples in DATASETS:
    data = parse_recommendations(dirname)
    records.append({
        "label":       label,
        "env":         env,
        "n_samples":   n_samples,
        "threshold":   data["threshold"],
        "range_lo":    data["range_lo"],
        "range_hi":    data["range_hi"],
        "coverage":    data["coverage"],
        "n_haplotypes":data["n_haplotypes"],
    })

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))

x = np.arange(len(records))
bar_width = 0.55

bars = []
for i, rec in enumerate(records):
    color = ENV_COLORS[rec["env"]]
    bar = ax.bar(
        i, rec["threshold"],
        width=bar_width,
        color=color,
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )
    bars.append(bar)

    # Acceptable range shown as a thin error-bar style bracket
    err_lo = rec["threshold"] - rec["range_lo"]
    err_hi = rec["range_hi"] - rec["threshold"]
    ax.errorbar(
        i, rec["threshold"],
        yerr=[[err_lo], [err_hi]],
        fmt="none",
        color="#333333",
        capsize=5,
        capthick=1.2,
        elinewidth=1.2,
        zorder=4,
    )

    # Annotate threshold value above each bar
    ax.text(
        i, rec["threshold"] + rec["range_hi"] - rec["threshold"] + 0.004,
        f"{rec['threshold']:.3f}",
        ha="center", va="bottom",
        fontsize=9, fontweight="bold",
        color="#222222",
    )

# Horizontal reference line at 0.015 (the consensus marine/freshwater threshold)
ax.axhline(0.015, color="#888888", linestyle="--", linewidth=1.0,
           zorder=2, label="0.015 (marine/freshwater consensus)")

# Axis formatting
ax.set_xticks(x)
ax.set_xticklabels([r["label"] for r in records], fontsize=10)
ax.set_ylabel("Optimal singleton distance threshold\n(min_singleton_distance)", fontsize=11)
ax.set_ylim(0, 0.115)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.3f}"))
ax.set_yticks([0, 0.015, 0.030, 0.045, 0.060, 0.075, 0.085, 0.100])
ax.grid(axis="y", alpha=0.35, zorder=1)
ax.set_axisbelow(True)

# N-samples annotation below x-axis labels
for i, rec in enumerate(records):
    ax.text(
        i, -0.011,
        f"n = {rec['n_samples']:,}",
        ha="center", va="top",
        fontsize=8, color="#555555",
        transform=ax.get_xaxis_transform(),
    )

# Legend — environments
legend_patches = [
    mpatches.Patch(facecolor=ENV_COLORS[e], edgecolor="white",
                   label=e.capitalize())
    for e in ["marine", "freshwater", "terrestrial"]
]
legend_patches.append(
    plt.Line2D([0], [0], color="#888888", linestyle="--", linewidth=1.2,
               label="0.015 consensus threshold")
)
ax.legend(
    handles=legend_patches,
    loc="upper left",
    fontsize=9,
    framealpha=0.9,
    edgecolor="#cccccc",
)

# Error-bar explanation note
ax.text(
    0.99, 0.97,
    "Error bars = acceptable threshold range",
    ha="right", va="top",
    transform=ax.transAxes,
    fontsize=8, color="#666666",
    style="italic",
)

plt.tight_layout(rect=[0, 0.04, 1, 1])

out_base = SWEEP_DIR / "Figure_4C_Threshold_Summary"
fig.savefig(str(out_base) + ".png", dpi=300, bbox_inches="tight")
fig.savefig(str(out_base) + ".pdf", bbox_inches="tight")
print(f"Saved: {out_base}.png")
print(f"Saved: {out_base}.pdf")

# ---------------------------------------------------------------------------
# Print data table for verification
# ---------------------------------------------------------------------------
print("\nSummary table:")
print(f"{'Dataset':<20} {'Env':<12} {'N':>6}  {'Threshold':>9}  {'Range':>12}  {'Coverage':>9}  {'Haplotypes':>10}")
print("-" * 88)
for rec in records:
    name = rec["label"].replace("\n", " ")
    print(
        f"{name:<20} {rec['env']:<12} {rec['n_samples']:>6}  "
        f"{rec['threshold']:>9.3f}  "
        f"{rec['range_lo']:.3f}–{rec['range_hi']:.3f}  "
        f"{rec['coverage']:>8.1f}%  "
        f"{rec['n_haplotypes']:>10}"
    )
