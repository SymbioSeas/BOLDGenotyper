"""Task 3.1: --extended gates interpretive and descriptive extras."""

from boldgenotyper import cli, config


def test_extended_defaults_off_and_parses():
    assert config.get_default_config().extended is False
    args = cli._build_main_parser().parse_args(["in.tsv", "--extended"])
    assert args.extended is True
    assert cli._build_main_parser().parse_args(["in.tsv"]).extended is False


def test_extended_is_in_the_serialized_config():
    cfg = config.get_default_config().update(extended=True)
    assert cfg.to_dict()["extended"] is True
