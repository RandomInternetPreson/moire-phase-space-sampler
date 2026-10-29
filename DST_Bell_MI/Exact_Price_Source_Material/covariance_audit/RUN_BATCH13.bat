@echo off
REM RUN_BATCH13.bat -- batch 13: the crest test (census_bumps_v1.py on the lattice census), 6 jobs, light; the nc32 ones take hours.
REM Runs directly (under the 8 GiB rule); fine alongside job 399. Run from covariance_audit.
REM Needs: census_bumps_v1.py, census_bandlimit_v4.py, run_parallel_v13.py, rho_nsh32_grid0_eps0.01_nang161.json + banks_nsh32_grid0_nang161_nw2048.npz.
cd /d "%~dp0"
echo batch 13 started %date% %time%  -- logs bumps_cr_^<tag^>.txt, status parallel13_status.txt, console parallel13_console.txt
python run_parallel_v13.py > parallel13_console.txt 2>&1
echo batch 13 finished %date% %time%
type parallel13_console.txt
pause
