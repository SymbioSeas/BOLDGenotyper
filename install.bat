@echo off
REM BOLDGenotyper Installation Script for Windows
REM This script sets up the conda environment and installs the package

echo =========================================
echo BOLDGenotyper Installation
echo =========================================
echo.

REM Check for conda
where conda >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Error: conda not found. Please install Miniconda or Anaconda first.
    echo Download from: https://docs.conda.io/en/latest/miniconda.html
    exit /b 1
)

echo Step 1: Creating conda environment...
conda env create -f environment.yml
if %ERRORLEVEL% NEQ 0 (
    echo Failed to create conda environment
    exit /b 1
)

echo.
echo Step 2: Activating environment...
call conda activate boldgenotyper

echo.
echo Step 3: Installing boldgenotyper package...
pip install -e .
if %ERRORLEVEL% NEQ 0 (
    echo Failed to install package
    exit /b 1
)

echo.
echo Step 4: Verifying installation...
boldgenotyper --version
if %ERRORLEVEL% EQU 0 (
    echo.
    echo =========================================
    echo Installation successful!
    echo =========================================
    echo.
    echo To use boldgenotyper, activate the environment:
    echo   conda activate boldgenotyper
    echo.
    echo Then run:
    echo   boldgenotyper --help
    echo.
) else (
    echo.
    echo Installation completed, but verification failed.
    echo Please activate the environment and try running boldgenotyper manually:
    echo   conda activate boldgenotyper
    echo   boldgenotyper --version
    exit /b 1
)
