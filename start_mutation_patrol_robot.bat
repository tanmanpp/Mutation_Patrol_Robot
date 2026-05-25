@echo off
setlocal

set "PROJECT_WIN=%~dp0"
if "%PROJECT_WIN:~-1%"=="\" set "PROJECT_WIN=%PROJECT_WIN:~0,-1%"
set "APP_URL=http://localhost:8000"

echo Starting Mutation Patrol Robot through WSL + conda...
echo Project: %PROJECT_WIN%
echo.

for /f "usebackq delims=" %%i in (`wsl wslpath -a "%PROJECT_WIN%"`) do set "PROJECT_WSL=%%i"

if "%PROJECT_WSL%"=="" (
  echo Failed to translate the project path into a WSL path.
  echo Make sure WSL is installed and available from Windows.
  pause
  exit /b 1
)

echo WSL project path: %PROJECT_WSL%
echo.
echo The first run may take a while because conda/npm dependencies are installed.
echo When the server is ready, open:
echo   %APP_URL%
echo.

start "" "%APP_URL%"

wsl bash -lc "cd '%PROJECT_WSL%' && bash scripts/run_app_wsl_conda.sh"

echo.
echo Mutation Patrol Robot stopped.
pause
