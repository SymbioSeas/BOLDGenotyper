#!/usr/bin/env python3
"""
Verify GOaS (Global Oceans and Seas) reference data setup.

This script checks whether the GOaS v1 shapefile dataset is correctly installed
for use with BOLDGenotyper. Because Marine Regions requires a registration form
before download, the shapefile must be obtained manually.

Usage:
    python -m boldgenotyper.goas_downloader [--data-dir PATH]

The GOaS dataset provides standardized ocean basin boundaries for
biogeographic analysis. Citation:

Costello MJ, Tsai P, Wong PS, Cheung AKL, Basher Z, Chaudhary C (2017)
Marine biogeographic realms and species endemicity. Nature Communications 8: 1057.
https://doi.org/10.1038/s41467-017-01121-2

Author: Steph Smith (symbioseas@outlook.com)
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# Expected files in the GOaS dataset
EXPECTED_FILES = [
    "World_Seas_IHO_v3.shp",
    "World_Seas_IHO_v3.shx",
    "World_Seas_IHO_v3.dbf",
    "World_Seas_IHO_v3.prj",
    "World_Seas_IHO_v3.cpg",
]

# Marine Regions download page (requires form submission)
GOAS_DOWNLOAD_PAGE = "https://www.marineregions.org/downloads.php"
GOAS_FILE_NAME = "GOaS_v1_20211214.zip"


def verify_files(data_dir: Path) -> bool:
    """
    Verify that all expected GOaS files are present.

    Returns True if all files exist, False otherwise.
    """
    logger.info("Verifying GOaS files...")
    missing = []

    for filename in EXPECTED_FILES:
        file_path = data_dir / filename
        if not file_path.exists():
            missing.append(filename)

    if missing:
        logger.error(f"Missing files: {missing}")
        return False

    logger.info("All GOaS files present")
    return True


def print_manual_instructions(data_dir: Path) -> None:
    """Print manual download instructions."""
    logger.info("")
    logger.info("=" * 60)
    logger.info("GOaS shapefile must be downloaded manually.")
    logger.info("Marine Regions requires a registration form before download.")
    logger.info("=" * 60)
    logger.info("")
    logger.info("Steps:")
    logger.info(f"  1. Visit: {GOAS_DOWNLOAD_PAGE}")
    logger.info(f"  2. Search for or locate '{GOAS_FILE_NAME}'")
    logger.info("  3. Fill in the required registration form")
    logger.info(f"  4. Download '{GOAS_FILE_NAME}'")
    logger.info(f"  5. Extract the ZIP and place the files into:")
    logger.info(f"       {data_dir}")
    logger.info("")
    logger.info("  After extraction, the directory should contain:")
    for f in EXPECTED_FILES:
        logger.info(f"    - {f}")
    logger.info("")
    logger.info("  Then re-run this script to verify the installation:")
    logger.info("    python -m boldgenotyper.goas_downloader")
    logger.info("")
    logger.info("  Or skip geographic analysis entirely with --no-geo:")
    logger.info("    boldgenotyper data/my_species.tsv --no-geo")
    logger.info("=" * 60)


def setup_goas(data_dir: Path) -> bool:
    """
    Verify GOaS data is present and print setup instructions if not.

    Parameters
    ----------
    data_dir : Path
        Directory where GOaS data should be installed.

    Returns
    -------
    bool
        True if all required files are present, False otherwise.
    """
    if verify_files(data_dir):
        logger.info("GOaS data is installed and ready.")
        logger.info(f"Location: {data_dir}")
        return True

    print_manual_instructions(data_dir)
    return False


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Verify GOaS reference data installation for BOLDGenotyper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
The GOaS shapefile must be downloaded manually from Marine Regions
because their site requires a registration form before download.

Steps:
  1. Visit: {GOAS_DOWNLOAD_PAGE}
  2. Download: {GOAS_FILE_NAME}
  3. Extract to: shapefiles/GOaS_v1_20211214/

For more information about the GOaS dataset, visit:
https://www.marineregions.org/
        """
    )

    parser.add_argument(
        '--data-dir',
        type=Path,
        default=Path(__file__).parent.parent / "shapefiles" / "GOaS_v1_20211214",
        help='Directory where GOaS data is installed (default: shapefiles/GOaS_v1_20211214)'
    )

    args = parser.parse_args()

    success = setup_goas(args.data_dir)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
