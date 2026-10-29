@echo off
REM RUN_BATCH9.bat -- batch 9: the (eps, K) surface on the lattice census: K=5,6,8,10 x eps=0.02,0.05,0.10 (12 light jobs, ~1 h).
REM Runs directly (under the 8 GiB rule); fine alongside job 399. Run from covariance_audit.
REM Needs: census_bandlimit_v4.py, run_parallel_v9.py, rho_nsh32_grid0_eps0.01_nang161.json + banks_nsh32_grid0_nang161_nw2048.npz.
cd /d "%~dp0"
echo batch 9 started %date% %time%  -- logs bl_g_e^<eps^>_k^<K^>.txt, status parallel9_status.txt, console parallel9_console.txt
python run_parallel_v9.py > parallel9_console.txt 2>&1
echo batch 9 finished %date% %time%
type parallel9_console.txt
pause
