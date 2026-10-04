"""A config file with the removed genotype_assignment section still loads."""

import logging

from boldgenotyper import config


def test_removed_section_is_ignored_with_a_warning(caplog):
    with caplog.at_level(logging.WARNING):
        cfg = config._dict_to_config({"genotype_assignment": {"min_identity": 0.9}})
    assert not hasattr(cfg, "genotype_assignment")
    assert "genotype_assignment" in caplog.text


def test_removed_visualization_options_are_ignored_with_a_warning(caplog):
    with caplog.at_level(logging.WARNING):
        cfg = config._dict_to_config({"visualization": {"figure_dpi": 600, "facet_by": "genotype",
                                                        "map_background_detail": "simple"}})
    v = cfg.visualization
    assert v.map_background_detail == "simple"
    assert not hasattr(v, "figure_dpi") and not hasattr(v, "facet_by")
    assert "figure_dpi" in caplog.text and "facet_by" in caplog.text
    assert set(vars(v)) == {"map_background_detail", "max_overview_haplotypes", "max_colored_items",
                            "max_overview_species", "max_heatmap_species"}
