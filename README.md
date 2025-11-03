## Installation

### 1. Clone the repository
```bash
git clone https://github.com/SymbioSeas/BOLDGenotyper.git
cd BOLDGenotyper
```

### 2. Create conda environment
```bash
conda env create -f environment.yml
conda activate boldgenotyper
```

### 3. Install the package
```bash
pip install -e .
```

### 4. Download GOaS reference data
```bash
python setup_goas.py
```

This downloads the Global Oceans and Seas (GOaS) dataset (~120 MB) required for ocean basin assignment.

**Alternative**: If you already have GOaS data, you can specify a custom path in your configuration file.
