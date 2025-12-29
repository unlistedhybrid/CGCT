@echo off
setlocal

REM Get directory of this BAT file
set SCRIPT_DIR=%~dp0

REM Convert "C:\path\to\CGCT\" -> "/mnt/c/path/to/CGCT"
set DRIVE=%SCRIPT_DIR:~0,1%
set REST=%SCRIPT_DIR:~3%
set REST=%REST:\=/%
set WSL_DIR=/mnt/%DRIVE%/%REST%

echo Running installer in WSL...
echo.

wsl bash -lc "cd '%WSL_DIR%' && chmod +x install-sibeliaz-wsl.sh && ./install-sibeliaz-wsl.sh"

echo.
echo If this failed:
echo 1) Open WSL manually
echo 2) cd %WSL_DIR%
echo 3) ./install-sibeliaz-wsl.sh
pause
