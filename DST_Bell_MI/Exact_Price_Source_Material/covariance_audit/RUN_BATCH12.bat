@echo off
REM RUN_BATCH12.bat -- batch 12: the crease test, birth_model_v4.py (sharp lobes allowed) on the lattice banks, 6 jobs x 3 h, light.
REM Runs directly (under the 8 GiB rule); fine alongside job 399 and batch 9. Run from covariance_audit.
REM Needs: birth_model_v4.py, run_parallel_v12.py, rho_nsh32_grid0_eps0.01_nang161.json + banks_nsh32_grid0_nang161_nw2048.npz,
REM   census_bandlimit_v4_n32_k6_rho_nsh32_grid0_eps0.01_nang161.json.
cd /d "%~dp0"
echo batch 12 started %date% %time%  -- logs bm4_^<tag^>.txt, status parallel12_status.txt, console parallel12_console.txt
python run_parallel_v12.py > parallel12_console.txt 2>&1
echo batch 12 finished %date% %time%
type parallel12_console.txt
pause
