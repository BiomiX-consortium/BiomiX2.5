@echo off
setlocal enabledelayedexpansion

echo -----------------------------------------
echo Uninstalling BiomiX components...
echo -----------------------------------------

:: Set base directory where this script is located
set SCRIPT_DIR=%~dp0

:: Define paths
set R_DIR=%SCRIPT_DIR%R_BiomiX
set PYTHON_DIR=%SCRIPT_DIR%Python_BiomiX
set PYTHON_ENV_DIR=%SCRIPT_DIR%BiomiX-env
set RENV_DIR=%SCRIPT_DIR%renv
set RPROFILE_FILE=%SCRIPT_DIR%.Rprofile

set PYTHON_VERSION=3.9.13
set SCRIPT_DIR=%~dp0
set INSTALL_DIR=%SCRIPT_DIR%Python_BiomiX

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

"%~dp0%INSTALLER%" 
:: Function to delete folder if exists
call :delete_if_exists "%R_DIR%" "R_BiomiX"
call :delete_if_exists "%PYTHON_ENV_DIR%" "Python BiomiX-env"
call :delete_if_exists "%RENV_DIR%" "renv environment"
call :delete_file_if_exists "%RPROFILE_FILE%" ".Rprofile file"

echo.
echo Uninstallation complete.
pause
exit /b 0

:delete_if_exists
set TARGET=%~1
set NAME=%~2
if exist "!TARGET!" (
    echo Deleting %NAME%...
    rmdir /s /q "!TARGET!"
    echo %NAME% removed.
) else (
    echo %NAME% not found — skipping.
)
exit /b

:delete_file_if_exists
set FILE=%~1
set NAME=%~2
if exist "!FILE!" (
    echo Deleting %NAME%...
    del /q "!FILE!"
    echo %NAME% removed.
) else (
    echo %NAME% not found — skipping.
)
exit /b
