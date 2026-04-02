@echo off
setlocal EnableDelayedExpansion

cd /d "%~dp0"
title RNA Hypothesis Workbench

set "ROOT_DIR=%CD%"
set "WSL_DRIVE=%ROOT_DIR:~0,1%"
if /I "%WSL_DRIVE%"=="A" set "WSL_DRIVE=a"
if /I "%WSL_DRIVE%"=="B" set "WSL_DRIVE=b"
if /I "%WSL_DRIVE%"=="C" set "WSL_DRIVE=c"
if /I "%WSL_DRIVE%"=="D" set "WSL_DRIVE=d"
if /I "%WSL_DRIVE%"=="E" set "WSL_DRIVE=e"
if /I "%WSL_DRIVE%"=="F" set "WSL_DRIVE=f"
if /I "%WSL_DRIVE%"=="G" set "WSL_DRIVE=g"
if /I "%WSL_DRIVE%"=="H" set "WSL_DRIVE=h"
if /I "%WSL_DRIVE%"=="I" set "WSL_DRIVE=i"
if /I "%WSL_DRIVE%"=="J" set "WSL_DRIVE=j"
if /I "%WSL_DRIVE%"=="K" set "WSL_DRIVE=k"
if /I "%WSL_DRIVE%"=="L" set "WSL_DRIVE=l"
if /I "%WSL_DRIVE%"=="M" set "WSL_DRIVE=m"
if /I "%WSL_DRIVE%"=="N" set "WSL_DRIVE=n"
if /I "%WSL_DRIVE%"=="O" set "WSL_DRIVE=o"
if /I "%WSL_DRIVE%"=="P" set "WSL_DRIVE=p"
if /I "%WSL_DRIVE%"=="Q" set "WSL_DRIVE=q"
if /I "%WSL_DRIVE%"=="R" set "WSL_DRIVE=r"
if /I "%WSL_DRIVE%"=="S" set "WSL_DRIVE=s"
if /I "%WSL_DRIVE%"=="T" set "WSL_DRIVE=t"
if /I "%WSL_DRIVE%"=="U" set "WSL_DRIVE=u"
if /I "%WSL_DRIVE%"=="V" set "WSL_DRIVE=v"
if /I "%WSL_DRIVE%"=="W" set "WSL_DRIVE=w"
if /I "%WSL_DRIVE%"=="X" set "WSL_DRIVE=x"
if /I "%WSL_DRIVE%"=="Y" set "WSL_DRIVE=y"
if /I "%WSL_DRIVE%"=="Z" set "WSL_DRIVE=z"
set "WSL_REST=%ROOT_DIR:~2%"
set "WSL_REST=%WSL_REST:\=/%"
set "WSL_DIR=/mnt/%WSL_DRIVE%%WSL_REST%"

echo ============================================================
echo RNA Hypothesis Workbench
echo ============================================================
echo.
echo Trying the best available runtime automatically:
echo - Windows R first
echo - WSL R second
echo - setup guide if neither is available
echo.

netstat -ano | findstr /R /C:":3838 .*LISTENING" >nul 2>nul
if not errorlevel 1 goto already_running

where Rscript >nul 2>nul
if not errorlevel 1 goto run_windows

where wsl.exe >nul 2>nul
if errorlevel 1 goto no_runtime

wsl.exe bash -lc "command -v Rscript >/dev/null 2>&1"
if not errorlevel 1 goto run_wsl

goto no_runtime

:already_running
echo The workbench already appears to be running on http://127.0.0.1:3838
echo.
echo The launcher will open that existing app instead of starting a second copy.
echo.
start "" "http://127.0.0.1:3838"
pause
goto end_ok

:run_windows
echo Starting with Windows R...
echo.
echo This terminal will stay open while the app runs.
echo If the browser does not open, visit http://127.0.0.1:3838
echo.
Rscript deseq2_workbench.R
goto stopped

:run_wsl
echo Windows R was not found, so the launcher is falling back to WSL R.
echo.
echo WSL workspace:
echo   !WSL_DIR!
echo.
echo This terminal will stay open while the app runs.
echo If the browser does not open, visit http://127.0.0.1:3838
echo.
wsl.exe bash -lc "cd '!WSL_DIR!' && bash launchers/start_workbench.sh"
goto stopped

:no_runtime
echo No usable R runtime was found.
echo.
echo Please open:
echo   START_HERE.txt
echo.
if exist "START_HERE.txt" start "" "START_HERE.txt"
echo.
pause
goto end_ok

:stopped
echo.
echo The app has stopped.
pause

:end_ok
endlocal
exit /b 0
