@echo off
REM RUN_PARALLEL3.bat -- batch 3: eps=0.10 row repair + certificates by value + miss curve (census_bandlimit_v4.py)
REM Run from the covariance_audit folder AFTER the swap job has finished (not alongside it: page-file exhaustion on 22 Oct).
REM Needs: census_bandlimit_v4.py, run_parallel_v3.py, the four rho_nsh32_*_eps0.01_nang161.json files and their banks_*.npz.
cd /d "%~dp0"
echo batch 3 started %date% %time%  -- per-job logs bl_^<tag^>.txt, status parallel3_status.txt, console parallel3_console.txt
python run_parallel_v3.py > parallel3_console.txt 2>&1
echo batch 3 finished %date% %time%
type parallel3_console.txt
pause
