@echo off
setlocal

:: Get full path to this script's folder
set SCRIPT_DIR=%~dp0

:: Set Rscript path relative to this script
set R_EXE=%SCRIPT_DIR%R_BiomiX\bin\Rscript.exe

:: Check R exists
if not exist "%R_EXE%" (
    echo Rscript not found at %R_EXE%
    pause
    exit /b 1
)

:: Test R is working
"%R_EXE%" -e "cat('R_BiomiX Rscript detected')" || echo Failed to run R_BiomiX Rscript

:: Run the script
echo Running install_from_github.R...
"%R_EXE%" "%SCRIPT_DIR%MODULE_WINDOWS_R_packages.r"

pause
