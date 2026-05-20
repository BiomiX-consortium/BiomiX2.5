@echo off
setlocal

:: Get the folder of this .bat script
set "SCRIPT_DIR=%~dp0"
set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"

:: Path to R environment
set "R_ENV_DIR=%SCRIPT_DIR%\_INSTALL\R_BiomiX"
set "R_EXE=%R_ENV_DIR%\bin\Rscript.exe"

:: Path to the R app script
set "APP_SCRIPT=%SCRIPT_DIR%\app.R"

:: Path to renv library
set "R_LIB=%SCRIPT_DIR%\_INSTALL\renv\library\windows\R-4.4\x86_64-w64-mingw32"

:: Check if Rscript.exe exists
if not exist "%R_EXE%" (
    echo R executable not found in: %R_EXE%
    pause
    exit /b 1
)

:: Check if app.R exists
if not exist "%APP_SCRIPT%" (
    echo app.R not found in: %APP_SCRIPT%
    pause
    exit /b 1
)

:: Confirm R version and path
echo Verifying R environment...
"%R_EXE%" -e "cat('Using R:', R.version$major, R.version$minor, '\n'); cat('R home:', R.home(), '\n')" || (
    echo R failed to run
    pause
    exit /b 1
)

:: Set PATH and R_HOME
set "R_HOME=%R_ENV_DIR%"
set "PATH=%R_ENV_DIR%\bin;%PATH%"

:: Disable renv auto-activate (we set the lib manually)
set "RENV_CONFIG_AUTOLOADER_ENABLED=FALSE"

:: Convert backslashes to forward slashes for R
set "R_LIB_R=%R_LIB:\=/%"
set "APP_SCRIPT_R=%APP_SCRIPT:\=/%"
set "SCRIPT_DIR_R=%SCRIPT_DIR:\=/%"

:: Launch the Shiny app opening the browser automatically
echo Launching R app...
"%R_EXE%" --vanilla -e ".libPaths('%R_LIB_R%'); shiny::runApp('%SCRIPT_DIR_R%', launch.browser=TRUE)"

pause