@echo off
REM RUN_PARALLEL.bat -- launches the bandwidth replication / eps-dependence jobs in parallel (run_parallel_v1.py)
REM Put next to census_bandlimit_v3.py, run_parallel_v1.py and the four rho files + banks in covariance_audit\
cd /d "%~dp0"
echo [%date% %time%] launching parallel jobs; watch parallel_status.txt and bl_<tag>.txt
python run_parallel_v1.py > parallel_console.txt 2>&1
echo [%date% %time%] done; see parallel_console.txt (SUMMARY of every job at the end)
pause
