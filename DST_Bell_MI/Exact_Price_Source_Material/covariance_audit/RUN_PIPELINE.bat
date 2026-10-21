@echo off
REM ============================================================================
REM  RUN_PIPELINE.bat -- self-driving terrain search (pipeline_v1.py), Oct 18 2026
REM  Put this file in covariance_audit\ next to: pipeline_v1.py, breed_product_v1.py,
REM  flight_time.py, dither_sweep_rho.py, product_census_v3.py, greedy_results.json
REM  Double-click it (or run it from a cmd window). It:
REM    0. smoke-tests all three scripts at toy scale (~2-5 min) and STOPS if anything is broken
REM    A. breeds a terrain for a station-local census from the greedy start (ROUNDS_A rounds, ~2 h each)
REM    V. if promising: validates the champion at full scale (~4 h banks + ~2 h LPs)
REM    B. else breeds from a one-term start (ROUNDS_B rounds), then validates or concludes
REM  Watch it:  pipeline_status.txt (every transition)   logs\<stage>.txt (live output)
REM  Stop it after the current stage:  create an empty file named STOP.txt in this folder
REM  Result:    pipeline_FINAL.txt
REM ============================================================================
cd /d "%~dp0"
set CORES=56
set ROUNDS_A=8
set ROUNDS_B=8
set PROMISING_RATIO=1.3
set PIPELINE_TINY=0
echo [%date% %time%] pipeline starting, console copy in pipeline_console.txt
python pipeline_v1.py
echo.
echo [%date% %time%] pipeline finished (exit %errorlevel%). See pipeline_FINAL.txt
pause
