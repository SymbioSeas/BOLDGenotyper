#!/usr/bin/env python3
"""
Download and setup GOaS (Global Oceans and Seas) reference data.

This script downloads the GOaS v1 shapefile dataset from Marine Regions
and sets it up for use with BOLDGenotyper.

Usage:
    python setup_goas.py [--data-dir PATH]

The GOaS dataset provides standardized ocean basin boundaries for
biogeographic analysis. Citation:

Costello MJ, Tsai P, Wong PS, Cheung AKL, Basher Z, Chaudhary C (2017)
Marine biogeographic realms and species endemicity. Nature Communications 8: 1057.
https://doi.org/10.1038/s41467-017-01121-2

Author: Steph Smith (steph.smith@unc.edu)
"""

import argparse
import hashlib
import logging
import shutil
import sys
import zipfile
from pathlib import Path
from urllib.request import urlretrieve
from urllib.error import URLError

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# GOaS v1 download information
GOAS_URL = "https://www.marineregions.org/download_file.php?name=World_Seas_IHO_v3.zip"
GOAS_ALTERNATIVE_URL = "https://github.com/iobis/mregions-static/raw/master/shapefiles/World_Seas_IHO_v3.zip"

# Expected files in the GOaS dataset
EXPECTED_FILES = [
    "World_Seas_IHO_v3.shp",
    "World_Seas_IHO_v3.shx",
    "World_Seas_IHO_v3.dbf",
    "World_Seas_IHO_v3.prj",
    "World_Seas_IHO_v3.cpg",
]

# MD5 checksum for verification (optional but recommended)
# Update this if the file changes
EXPECTED_MD5 = None  # Set to actual MD5 if known


def calculate_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Calculate MD5 checksum of a file."""
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)
    return md5.hexdigest()


def download_file(url: str, output_path: Path) -> bool:
    """
    Download file from URL with progress reporting.
    
    Returns True if successful, False otherwise.
    """
    try:
        logger.info(f"Downloading from {url}")
        logger.info("This may take a few minutes...")
        
        def progress_hook(block_num, block_size, total_size):
            downloaded = block_num * block_size
            if total_size > 0:
                percent = min(100, (downloaded / total_size) * 100)
                if block_num % 50 == 0:  # Report every ~4MB
                    logger.info(f"Downloaded: {percent:.1f}% ({downloaded / 1e6:.1f} MB)")
        
        urlretrieve(url, output_path, reporthook=progress_hook)
        logger.info(f"Download complete: {output_path}")
        return True
        
    except URLError as e:
        logger.error(f"Download failed: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error during download: {e}")
        return False


def extract_zip(zip_path: Path, extract_dir: Path) -> bool:
    """
    Extract ZIP file to directory.
    
    Returns True if successful, False otherwise.
    """
    try:
        logger.info(f"Extracting {zip_path.name}...")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        logger.info(f"Extraction complete: {extract_dir}")
        return True
        
    except zipfile.BadZipFile:
        logger.error(f"Invalid ZIP file: {zip_path}")
        return False
    except Exception as e:
        logger.error(f"Extraction failed: {e}")
        return False


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
    
    logger.info("✓ All GOaS files present")
    return True


def create_citation_file(data_dir: Path) -> None:
    """Create a citation file for the GOaS dataset."""
    citation = """GOaS (Global Oceans and Seas) Dataset v1
=====================================================

Citation:
Costello MJ, Tsai P, Wong PS, Cheung AKL, Basher Z, Chaudhary C (2017)
Marine biogeographic realms and species endemicity. 
Nature Communications 8: 1057.
https://doi.org/10.1038/s41467-017-01121-2

Data Source:
Marine Regions - https://www.marineregions.org/

License:
Creative Commons Attribution 4.0 International (CC BY 4.0)
https://creativecommons.org/licenses/by/4.0/

Downloaded: {date}

When using this dataset, please cite the original publication and
acknowledge Marine Regions.
""".format(date=Path(__file__).stat().st_mtime)

    citation_file = data_dir / "CITATION.txt"
    with open(citation_file, 'w') as f:
        f.write(citation)
    
    logger.info(f"Created citation file: {citation_file}")


def setup_goas(data_dir: Path, force: bool = False) -> bool:
    """
    Main setup function for GOaS data.
    
    Parameters
    ----------
    data_dir : Path
        Directory to install GOaS data
    force : bool
        Force re-download even if files exist
    
    Returns
    -------
    bool
        True if setup successful, False otherwise
    """
    # Check if already installed
    if not force and verify_files(data_dir):
        logger.info("GOaS data already installed and verified")
        logger.info(f"Location: {data_dir}")
        return True
    
    # Create data directory
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Download ZIP file
    zip_path = data_dir / "goas_download.zip"
    
    # Try primary URL
    success = download_file(GOAS_URL, zip_path)
    
    # Try alternative URL if primary fails
    if not success:
        logger.info("Trying alternative download location...")
        success = download_file(GOAS_ALTERNATIVE_URL, zip_path)
    
    if not success:
        logger.error("Failed to download GOaS data from all sources")
        logger.error("Please download manually from:")
        logger.error("  https://www.marineregions.org/downloads.php")
        return False
    
    # Verify download with checksum if available
    if EXPECTED_MD5:
        logger.info("Verifying download...")
        actual_md5 = calculate_md5(zip_path)
        if actual_md5 != EXPECTED_MD5:
            logger.error(f"Checksum mismatch! Expected {EXPECTED_MD5}, got {actual_md5}")
            logger.error("The downloaded file may be corrupted")
            return False
        logger.info("✓ Checksum verified")
    
    # Extract files
    if not extract_zip(zip_path, data_dir):
        return False
    
    # Verify extracted files
    if not verify_files(data_dir):
        return False
    
    # Create citation file
    create_citation_file(data_dir)
    
    # Clean up ZIP file
    try:
        zip_path.unlink()
        logger.info("Cleaned up download file")
    except Exception as e:
        logger.warning(f"Could not delete ZIP file: {e}")
    
    logger.info("=" * 60)
    logger.info("GOaS setup complete!")
    logger.info(f"Data location: {data_dir}")
    logger.info("=" * 60)
    
    return True


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Download and setup GOaS reference data for BOLDGenotyper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Install to default location
  python setup_goas.py
  
  # Install to custom location
  python setup_goas.py --data-dir /path/to/data
  
  # Force re-download
  python setup_goas.py --force

For more information about the GOaS dataset, visit:
https://www.marineregions.org/
        """
    )
    
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=Path(__file__).parent / "GOaS_v1_20211214",
        help='Directory to install GOaS data (default: ./GOaS_v1_20211214)'
    )
    
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-download even if files exist'
    )
    
    args = parser.parse_args()
    
    # Run setup
    success = setup_goas(args.data_dir, force=args.force)
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
