"""A config file with the removed genotype_assignment section still loads."""

import logging

from boldgenotyper import config


def test_removed_section_is_ignored_with_a_warning(caplog):
    with caplog.at_level(logging.WARNING):
        cfg = config._dict_to_config({"genotype_assignment": {"min_identity": 0.9}})
    assert not hasattr(cfg, "genotype_assignment")
    assert "genotype_assignment" in caplog.text
