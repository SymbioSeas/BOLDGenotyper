"""Pytest configuration for the BOLDGenotyper test suite.

Some files under tests/ are standalone scripts rather than test modules: they
define no test functions and run code at import time (call sys.exit(), write
logs, load map data, or read local run outputs that are not in the
repository), or they import modules that no longer exist. Excluding them lets
a bare `pytest` collect the suite without side effects. The files are kept for
manual use (`python tests/<file>.py`).
"""

collect_ignore = [
    # Scripts with no test functions
    "test_coord_inference_local.py",
    "test_full_pipeline.py",
    "test_label_fix.py",            # exits unless tests/Euprymna_output/ exists
    "test_my_data.py",              # exits at import
    "test_phylo_viz_pipeline.py",
    "test_pipeline_script.py",      # writes tests/test_run.log at import
    "test_snap_to_coast_local.py",  # loads Natural Earth land geometries at import
    "test_tree_scaling.py",
    "test_viz_warnings.py",         # needs a local run in tests/Carcharhinus_output/
]
