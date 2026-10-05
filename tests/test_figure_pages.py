"""Figure paging helpers (Phase 11)."""

import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from boldgenotyper import config, figure_pages as fp


def _n_pdf_pages(path):
    return len(re.findall(rb"/Type\s*/Page[^s]", path.read_bytes()))


def test_cap_items():
    counts = pd.Series({f"H{i:02d}": 20 - i for i in range(14)})
    labels = fp.cap_items(counts, 12, noun="haplotypes")
    assert labels["H00"] == "H00" and labels["H10"] == "H10"
    assert set(labels[["H11", "H12", "H13"]]) == {"Other (3 haplotypes)"}
    assert (fp.cap_items(counts.head(12), 12) == counts.head(12).index.to_series()).all()


def test_colors_and_groups():
    haps = pd.DataFrame({"haplotype_id": ["H1", "H2", "H3", "H4"],
                         "assigned_sp": ["Genus beta", "Genus alpha", "Genus alpha", None],
                         "n_members": [5, 2, 9, 1]})
    groups = fp.species_groups(haps)
    assert list(groups) == ["Genus alpha", "Genus beta", "Unassigned"]
    assert groups["Genus alpha"] == ["H3", "H2"]
    colors = fp.item_colors(["a", "b", "Other (2)"])
    assert colors["Other (2)"] == fp.OTHER_COLOR and colors["a"] != colors["b"]
    assert set(fp.haplotype_colors(haps)) == {"H1", "H2", "H3", "H4"}
    assert fp.overview_level(40) == "haplotype" and fp.overview_level(41) == "species"
    assert fp.slug("Genus alpha (genus level)") == "Genus_alpha_genus_level"


def test_write_paged_figure_and_manifest(tmp_path):
    def fig(title):
        f, ax = plt.subplots()
        ax.set_title(title)
        return f

    groups = {"Genus alpha": ["H1"], "Genus beta": ["H2"], "Genus gamma": ["H3"]}
    entry = fp.write_paged_figure("distribution_map", tmp_path, lambda: fig("overview"),
                                  lambda g, ids: None if g == "Genus gamma" else fig(g), groups)
    assert entry["overview"] == ["distribution_map.pdf", "distribution_map.svg"]
    assert entry["pages_pdf"] == "distribution_map_by_species.pdf"
    assert entry["pages"] == {"Genus alpha": "by_species/distribution_map/Genus_alpha.svg",
                              "Genus beta": "by_species/distribution_map/Genus_beta.svg"}
    assert _n_pdf_pages(tmp_path / "distribution_map_by_species.pdf") == 2
    assert "<text" in (tmp_path / "by_species" / "distribution_map" / "Genus_alpha.svg").read_text()
    fp.write_manifest(tmp_path, {"distribution_map": entry})
    fp.write_manifest(tmp_path, {"other": {"overview": ["other.pdf"]}})
    assert set(json.loads((tmp_path / "manifest.json").read_text())) == {"distribution_map", "other"}


def test_visualization_config_thresholds():
    v = config.get_default_config().visualization
    assert (v.max_overview_haplotypes, v.max_colored_items, v.max_overview_species,
            v.max_heatmap_species) == (40, 12, 15, 30)


def test_count_agrees_in_number():
    assert fp.count(1, "record") == "1 record"
    assert fp.count(1234, "record") == "1,234 records"
    assert fp.count(2, "taxon", "taxa") == "2 taxa"
