@echo off
setlocal

:: Set installer name and install path
set R_VERSION=4.4.1
set R_INSTALLER=R-%R_VERSION%-win.exe
set SCRIPT_DIR=%~dp0
set R_INSTALL_PATH=%SCRIPT_DIR%R_BiomiX

:: Check installer exists
if not exist "%~dp0%R_INSTALLER%" (
    echo R installer "%R_INSTALLER%" not found in current directory.
    pause
    exit /b 1
)

:: Install R silently
echo Installing R %R_VERSION% to %R_INSTALL_PATH%...
"%~dp0%R_INSTALLER%" /VERYSILENT /DIR="%R_INSTALL_PATH%"

if not exist "%R_INSTALL_PATH%\bin\Rscript.exe" (
    echo R installation failed.
    pause
    exit /b 1
)

:: Install renv using Rscript
echo Installing renv package...
"%R_INSTALL_PATH%\bin\Rscript.exe" -e "install.packages('renv', repos='https://cloud.r-project.org')"

:: Optional: initialize renv in project folder
echo Initializing renv in project folder (optional)...
"%R_INSTALL_PATH%\bin\Rscript.exe" -e "renv::init(bare = TRUE)"

echo.
echo R and renv successfully installed.
echo R location: %R_INSTALL_PATH%
pause
