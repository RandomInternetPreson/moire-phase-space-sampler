@echo off
REM RUN_PARALLEL2.bat -- second parallel batch: K=5 certificates on strat0/strat1 and the eps=0.02 / 0.10 rows of the (eps,K) surface
REM Put next to census_bandlimit_v3.py, run_parallel_v1.py and the four rho files + banks in covariance_audit\
cd /d "%~dp0"
echo [%date% %time%] launching parallel jobs; watch parallel2_status.txt and the bl_ logs
python run_parallel_v2.py > parallel2_console.txt 2>&1
echo [%date% %time%] done; see parallel2_console.txt (SUMMARY of every job at the end)
pause
