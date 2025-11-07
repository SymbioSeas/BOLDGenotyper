from setuptools import setup, find_packages

setup(
    name="boldgenotyper",
    version="0.1.0",
    description="Automated COI genotyping from BOLD database data",
    author="Steph Smith",
    author_email="steph.smith@unc.edu",
    url="https://github.com/SymbioSeas/BOLDGenotyper",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
    ],
    entry_points={
        'console_scripts': [
            'boldgenotyper=boldgenotyper.cli:main',
        ],
    },
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
