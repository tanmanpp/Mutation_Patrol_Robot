@echo off
setlocal

set "PROJECT_WIN=%~dp0"
if "%PROJECT_WIN:~-1%"=="\" set "PROJECT_WIN=%PROJECT_WIN:~0,-1%"
set "APP_URL=http://localhost:8000"

echo Starting Mutation Patrol Robot through WSL + conda...
echo Project: %PROJECT_WIN%
echo.

wsl --status >nul 2>&1
if errorlevel 1 (
  echo [ERROR] WSL is not installed or is not ready.
  echo Run "wsl --install" from an Administrator PowerShell, restart Windows,
  echo then open this launcher again.
  pause
  exit /b 1
)

for /f "usebackq delims=" %%i in (`wsl wslpath -a "%PROJECT_WIN%"`) do set "PROJECT_WSL=%%i"

if "%PROJECT_WSL%"=="" (
  echo Failed to translate the project path into a WSL path.
  echo Make sure WSL is installed and available from Windows.
  pause
  exit /b 1
)

echo WSL project path: %PROJECT_WSL%
echo.
echo The launcher will verify the conda environment and required tools first.
echo The first run may take a while while dependencies are installed.
echo When the server is ready, open:
echo   %APP_URL%
echo.

wsl bash -lc "cd '%PROJECT_WSL%' && bash scripts/run_app_wsl_conda.sh"

echo.
echo Mutation Patrol Robot stopped.
pause
