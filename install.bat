@echo off
REM BOLDGenotyper Installation Helper for Windows
REM
REM BOLDGenotyper requires bioconda packages (MAFFT, trimAl, FastTree)
REM which are only available on Linux and macOS. On Windows, you must
REM use WSL2 (Windows Subsystem for Linux).
REM
REM This script checks for WSL2 and launches the Linux installer inside it.

echo =========================================
echo BOLDGenotyper Installation (Windows)
echo =========================================
echo.

REM Check for WSL
where wsl >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: WSL2 is required to run BOLDGenotyper on Windows.
    echo.
    echo BOLDGenotyper depends on bioconda packages (MAFFT, trimAl, FastTree^)
    echo which are not available on native Windows.
    echo.
    echo To install WSL2:
    echo   1. Open PowerShell as Administrator
    echo   2. Run: wsl --install
    echo   3. Restart your computer
    echo   4. Install Miniconda inside WSL:
    echo      https://docs.conda.io/en/latest/miniconda.html#linux-installers
    echo   5. Re-run this script
    echo.
    echo For details: https://learn.microsoft.com/en-us/windows/wsl/install
    exit /b 1
)

echo WSL2 detected. Launching Linux installer inside WSL...
echo.

REM Run the Linux install script inside WSL
wsl bash -c "cd '%cd%' && bash install.sh"
if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Installation failed. Please check the output above for errors.
    echo.
    echo You can also install manually inside WSL:
    echo   1. Open WSL terminal
    echo   2. cd to the BOLDGenotyper directory
    echo   3. Run: bash install.sh
    exit /b 1
)

echo.
echo =========================================
echo Installation complete!
echo =========================================
echo.
echo To use BOLDGenotyper, open WSL and run:
echo   conda activate boldgenotyper
echo   boldgenotyper --help
echo.
