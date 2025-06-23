@echo off
setlocal

:: Get the folder of this .bat script
set "SCRIPT_DIR=%~dp0"
set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"

:: Path to BiomiX Python environment
set "PYTHON_ENV_DIR=%SCRIPT_DIR%\_INSTALL\BiomiX-env"
set "PYTHON_EXE=%PYTHON_ENV_DIR%\Scripts\python.exe"
set "PYTHON_SCRIPTS=%PYTHON_ENV_DIR%\Scripts"

:: Path to launcher script
set "LAUNCHER_SCRIPT=%SCRIPT_DIR%\MODULE_WINDOWS.py"

:: Check if python.exe exists
if not exist "%PYTHON_EXE%" (
    echo Python executable not found in: %PYTHON_EXE%
    pause
    exit /b 1
)

:: Confirm Python version and path
echo 🔍 Verifying Python environment...
"%PYTHON_EXE%" -c "import sys; print('Using:', sys.executable); print('site-packages:', sys.path)" || (
    echo Python failed to run
    pause
    exit /b 1
)

:: Set PATH to ensure the correct environment is loaded
set "PATH=%PYTHON_ENV_DIR%;%PYTHON_SCRIPTS%;%PATH%"
set "PYTHONPATH="

:: Launch BiomiX interface
echo Launching BiomiX...
"%PYTHON_EXE%" "%LAUNCHER_SCRIPT%"

pause
