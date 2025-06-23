@echo off
setlocal enabledelayedexpansion

:: Settings
set PYTHON_VERSION=3.9.13
:: Install Python in the same directory as the BAT file
set SCRIPT_DIR=%~dp0
set INSTALL_DIR=%SCRIPT_DIR%Python_BiomiX
:: Environment name and path
set ENV_NAME=BiomiX-env
set ENV_DIR=%SCRIPT_DIR%%ENV_NAME%

:: Detect system architecture
set ARCHITECTURE=x86
if "%PROCESSOR_ARCHITECTURE%"=="AMD64" set ARCHITECTURE=amd64

:: Choose installer based on architecture
set INSTALLER_32=python-%PYTHON_VERSION%.exe
set INSTALLER_64=python-%PYTHON_VERSION%-amd64.exe

if "%ARCHITECTURE%"=="amd64" (
    set INSTALLER=%INSTALLER_64%
) else (
    set INSTALLER=%INSTALLER_32%
)

echo %ARCHITECTURE%

:: Verify installer exists
if not exist "%~dp0%INSTALLER%" (
    echo Python installer "%INSTALLER%" not found in current directory.
    pause
    exit /b 1
)

:: Install Python silently
echo Installing Python %PYTHON_VERSION% to %INSTALL_DIR%
"%~dp0%INSTALLER%" /quiet InstallAllUsers=0 PrependPath=1 TargetDir="%INSTALL_DIR%" Include_launcher=1

if not exist "%INSTALL_DIR%\python.exe" (
    echo Installation failed. Exiting.
    pause
    exit /b 1
)

:: Create virtual environment
echo Creating virtual environment: %ENV_NAME%
"%INSTALL_DIR%\python.exe" -m venv "%ENV_DIR%"

:: Activate the environment
call "%ENV_DIR%\Scripts\activate.bat"

:: Upgrade pip and install packages
python -m pip install --upgrade pip setuptools wheel

:: Install packages
python -m pip install mamba
python -m pip install PyQt5==5.12.3 pandas==2.2.3 scikit-learn==1.6.1 ^
    xlrd==2.0.1 openpyxl==3.1.5 mofapy2==0.7.1 numpy==1.26.4

echo.
echo Python Environment setup complete!
echo To activate it in the future, run:
echo call "%ENV_DIR%\Scripts\activate"
pause
