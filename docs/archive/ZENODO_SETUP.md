# Zenodo DOI Setup Guide

This guide explains how to set up automatic DOI generation for BOLDGenotyper releases using Zenodo.

## What is Zenodo?

[Zenodo](https://zenodo.org/) is a research data repository that:
- Provides **permanent Digital Object Identifiers (DOIs)** for software releases
- Ensures **long-term preservation** of your code
- Integrates seamlessly with GitHub
- Is **free** for open-source projects
- Is recommended by many journals for software citation

## Why Use Zenodo for BOLDGenotyper?

1. **Citability**: Provides a permanent DOI that can be cited in publications
2. **Archival**: Ensures BOLDGenotyper is preserved even if GitHub changes
3. **Versioning**: Creates a separate DOI for each release
4. **Discoverability**: Makes BOLDGenotyper findable by researchers worldwide
5. **Publication Requirements**: Many journals require a DOI for software citations

---

## Setup Instructions

### Step 1: Create a Zenodo Account

1. Go to [https://zenodo.org/](https://zenodo.org/)
2. Click **"Sign up"** or **"Log in"**
3. Choose **"Log in with GitHub"** for easiest integration
4. Authorize Zenodo to access your GitHub account

### Step 2: Link BOLDGenotyper Repository to Zenodo

1. After logging in, click on your **username** in the top right
2. Select **"GitHub"** from the dropdown menu
3. You'll see a list of your GitHub repositories
4. Find **"BOLDGenotyper"** in the list
5. Toggle the switch to **"ON"** for BOLDGenotyper
   - This enables automatic DOI creation for new releases

**Important**: You must have admin access to the repository to enable Zenodo integration.

### Step 3: Create a GitHub Release

Zenodo creates a new DOI every time you create a GitHub release.

#### First Release (v0.1.0)

1. Go to your GitHub repository: `https://github.com/SymbioSeas/BOLDGenotyper`
2. Click on **"Releases"** (right sidebar)
3. Click **"Create a new release"**
4. Fill in the release information:

**Tag version**: `v0.1.0`
- Format: `vX.Y.Z` (semantic versioning)
- X = major version, Y = minor version, Z = patch version

**Release title**: `BOLDGenotyper v0.1.0 - Initial Release`

**Description**:
```markdown
# BOLDGenotyper v0.1.0 - Initial Release

This is the first stable release of BOLDGenotyper, a Python package for automated COI sequence genotyping and biogeographic analysis.

## Key Features
- Complete pipeline from BOLD TSV to annotated results
- Sequence dereplication and consensus generation
- Genotype assignment with CIGAR-based identity calculation
- Geographic analysis with GOaS integration
- Phylogenetic tree building (optional)
- Publication-ready visualizations

## Installation
```bash
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper
conda env create -f boldgenotyper_env.yml
conda activate depredation
```

## Quick Start
See [README.md](https://github.com/SymbioSeas/BOLDGenotyper/blob/main/README.md) for detailed instructions.

## Citation
If you use this software, please cite:
[Citation information from CITATION.cff]

## Full Changelog
See [CHANGELOG.md](https://github.com/SymbioSeas/BOLDGenotyper/blob/main/CHANGELOG.md)
```

5. **Attach files** (optional):
   - Pre-built packages
   - Example datasets
   - Documentation PDFs

6. Check **"This is a pre-release"** if appropriate (for v0.x.x versions)

7. Click **"Publish release"**

### Step 4: Retrieve Your DOI from Zenodo

1. After creating the GitHub release, go to [Zenodo](https://zenodo.org/)
2. Click on **"Upload"** → **"GitHub"** in the top menu
3. You should see BOLDGenotyper listed with your new release
4. Click on the release to view details
5. Copy the **DOI** (format: `10.5281/zenodo.XXXXXXX`)

**The DOI badge will look like**:
```
DOI: 10.5281/zenodo.XXXXXXX
```

### Step 5: Add DOI to Your Repository

Once you have the DOI, update the following files:

#### 1. Update README.md Badge

Replace the DOI badge at the top of README.md:

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

#### 2. Update CITATION.cff

Update the `doi` field in CITATION.cff:

```yaml
doi: 10.5281/zenodo.XXXXXXX
```

#### 3. Update README Citation Section

In the Software Citation BibTeX entry:

```bibtex
@software{boldgenotyper2025,
  author = {Smith, Steph},
  title = {BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis},
  year = {2025},
  version = {0.1.0},
  url = {https://github.com/SymbioSeas/BOLDGenotyper},
  doi = {10.5281/zenodo.XXXXXXX}  ← Update this
}
```

#### 4. Commit and Push Changes

```bash
git add README.md CITATION.cff
git commit -m "Add Zenodo DOI badge and citation info"
git push origin main
```

---

## Zenodo Metadata

When you create a Zenodo entry, you can customize the metadata:

### Recommended Metadata for BOLDGenotyper

1. **Title**: BOLDGenotyper: Automated COI Sequence Genotyping and Biogeographic Analysis
2. **Authors**:
   - Steph Smith (University of North Carolina, Institute of Marine Sciences)
3. **Description**: Use the abstract from CITATION.cff
4. **Version**: 0.1.0
5. **License**: MIT License
6. **Keywords**:
   - bioinformatics
   - DNA barcoding
   - COI
   - genotyping
   - phylogeography
   - biogeography
   - BOLD database
7. **Related identifiers**:
   - GitHub repository: https://github.com/SymbioSeas/BOLDGenotyper
   - Related publication DOI (when available)

### Editing Metadata

1. Go to your Zenodo record
2. Click **"Edit"**
3. Update metadata fields
4. Click **"Save"** → **"Publish"**

**Note**: Some fields cannot be changed after publication. Edit carefully before publishing.

---

## Best Practices

### Version Numbering

Follow [Semantic Versioning](https://semver.org/):
- **MAJOR** (X.0.0): Incompatible API changes
- **MINOR** (0.X.0): New features, backwards compatible
- **PATCH** (0.0.X): Bug fixes, backwards compatible

Examples:
- `v0.1.0`: Initial release
- `v0.1.1`: Bug fix
- `v0.2.0`: New features (terrestrial support)
- `v1.0.0`: Stable API, ready for publication

### Release Timing

Create a new release when:
- Initial stable version is ready (v0.1.0) ✓
- Significant new features are added (v0.2.0)
- Critical bug fixes are made (v0.1.1)
- Preparing for publication (v1.0.0)

### DOI Usage

**Concept DOI vs Version DOI**:
- **Concept DOI**: Always points to latest version (e.g., `10.5281/zenodo.1234567`)
- **Version DOI**: Specific to one version (e.g., `10.5281/zenodo.1234568`)

**For citations**:
- Use **Version DOI** for reproducibility
- Use **Concept DOI** for general reference to the software

### GitHub Release Notes

Good release notes should include:
- **What's New**: Key features and improvements
- **Breaking Changes**: Incompatible changes (if any)
- **Bug Fixes**: Fixed issues
- **Installation**: How to install this version
- **Citation**: How to cite this version
- **Full Changelog**: Link to CHANGELOG.md

---

## Troubleshooting

### "Repository not showing in Zenodo"

**Possible causes**:
1. Repository is private (Zenodo only works with public repos)
2. Zenodo-GitHub sync is delayed (wait 5-10 minutes)
3. Need to refresh the page

**Solution**: Make sure repository is public, wait a few minutes, refresh the Zenodo GitHub page.

### "DOI not created after release"

**Possible causes**:
1. Zenodo integration wasn't enabled before creating the release
2. Release was created as a draft
3. Zenodo webhook failed

**Solution**:
- Check that toggle is ON in Zenodo → GitHub
- Publish the release (not just save as draft)
- Delete and recreate the release if necessary

### "Want to update DOI metadata"

**Solution**:
1. Go to Zenodo record
2. Click "Edit"
3. Update fields
4. Click "Publish"

**Note**: Some fields (like DOI, publication date) cannot be changed.

### "Multiple DOIs for same version"

**Cause**: Created multiple releases with same tag

**Solution**:
- Delete duplicate releases from GitHub
- Keep only the correct one
- The DOI is permanent, but you can mark old versions as "superseded"

---

## Integration with Publication

### For Journal Submission

When submitting your manuscript that uses BOLDGenotyper:

1. **Methods Section**: Cite the software with version and DOI
   ```
   "Genotyping was performed using BOLDGenotyper v0.1.0 (Smith, 2025;
   DOI: 10.5281/zenodo.XXXXXXX), an automated pipeline for..."
   ```

2. **Code Availability**: Include DOI in code availability statement
   ```
   "BOLDGenotyper is freely available at https://github.com/SymbioSeas/BOLDGenotyper
   and archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX)."
   ```

3. **Supplementary Materials**: Can include Zenodo link for permanent access

### For Software Paper (JOSS, etc.)

If submitting to Journal of Open Source Software:

1. JOSS requires a DOI (Zenodo is recommended)
2. Include DOI in `paper.md`
3. JOSS review will check that DOI resolves correctly

---

## Automation (Advanced)

### Automated Releases with GitHub Actions

You can automate the release process:

Create `.github/workflows/release.yml`:

```yaml
name: Create Release

on:
  push:
    tags:
      - 'v*'

jobs:
  release:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout code
        uses: actions/checkout@v3

      - name: Create Release
        uses: actions/create-release@v1
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        with:
          tag_name: ${{ github.ref }}
          release_name: Release ${{ github.ref }}
          body_path: CHANGELOG.md
          draft: false
          prerelease: false
```

Then to create a release:

```bash
git tag -a v0.1.0 -m "Release version 0.1.0"
git push origin v0.1.0
```

GitHub Actions will automatically create the release, and Zenodo will automatically create a DOI.

---

## Checklist for First Release

- [ ] Zenodo account created and linked to GitHub
- [ ] BOLDGenotyper repository toggle enabled in Zenodo
- [ ] All code committed and pushed to main branch
- [ ] CHANGELOG.md updated with v0.1.0 changes
- [ ] README.md complete and up-to-date
- [ ] LICENSE file present
- [ ] CITATION.cff present
- [ ] GitHub release created (v0.1.0)
- [ ] DOI generated by Zenodo
- [ ] DOI badge added to README.md
- [ ] DOI added to CITATION.cff
- [ ] Citation section updated in README.md
- [ ] Changes committed and pushed

---

## Resources

- **Zenodo Help**: https://help.zenodo.org/
- **GitHub Releases**: https://docs.github.com/en/repositories/releasing-projects-on-github
- **Semantic Versioning**: https://semver.org/
- **Software Citation Principles**: https://www.force11.org/software-citation-principles
- **Making Your Code Citable**: https://guides.github.com/activities/citable-code/

---

## Contact

If you have questions about Zenodo setup for BOLDGenotyper:
- **Email**: steph.smith@unc.edu
- **GitHub Issues**: https://github.com/SymbioSeas/BOLDGenotyper/issues

---

**Last Updated**: 2025-11-17
**Version**: 1.0
