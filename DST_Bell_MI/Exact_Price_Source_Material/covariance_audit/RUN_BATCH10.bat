@echo off
REM RUN_BATCH10.bat -- batch 10: birth_model_v3.py on the lattice banks, 9 jobs x 3 h, light (< 2 GB each). Overnight.
REM Runs directly (under the 8 GiB rule); fine alongside job 399 and batch 9. Run from covariance_audit.
REM Needs: birth_model_v3.py, run_parallel_v10.py, rho_nsh32_grid0_eps0.01_nang161.json + banks_nsh32_grid0_nang161_nw2048.npz,
REM   census_bandlimit_v4_n32_k6_rho_nsh32_grid0_eps0.01_nang161.json.
cd /d "%~dp0"
echo batch 10 started %date% %time%  -- logs bm3_^<tag^>.txt, status parallel10_status.txt, console parallel10_console.txt
python run_parallel_v10.py > parallel10_console.txt 2>&1
echo batch 10 finished %date% %time%
type parallel10_console.txt
pause
