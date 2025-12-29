@echo off
setlocal

REM Check that WSL exists
wsl --status >nul 2>&1 || (
  echo WSL is not installed.
  echo Install it with: wsl --install
  pause
  exit /b 1
)

REM Directory containing this BAT file
set SCRIPT_DIR=%~dp0

REM Convert Windows path to WSL path using wslpath
for /f "delims=" %%i in ('wsl wslpath "%SCRIPT_DIR%"') do set WSL_DIR=%%i

echo Running installer in WSL...
echo.

wsl bash -lc "cd '%WSL_DIR%' && chmod +x install-sibeliaz-wsl.sh && ./install-sibeliaz-wsl.sh"

echo.
echo If this failed:
echo 1) Open WSL manually
echo 2) cd %WSL_DIR%
echo 3) ./install-sibeliaz-wsl.sh
pause
