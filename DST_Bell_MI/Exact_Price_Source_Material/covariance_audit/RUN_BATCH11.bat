@echo off
REM RUN_BATCH11.bat -- batch 11: the pleat scan on the lattice census (census_pleat_v1.py): 14 "lines" jobs (~10 min) + 3 "blocks" jobs (~4 h).
REM 17 light processes (< 3 GB each), runs directly; fine alongside job 399. Run from covariance_audit.
REM Needs: census_pleat_v1.py, run_parallel_v11.py, rho_nsh32_grid0_eps0.01_nang161.json + banks_nsh32_grid0_nang161_nw2048.npz.
cd /d "%~dp0"
echo batch 11 started %date% %time%  -- logs bl_pl_N^<N^>_^<mode^>.txt, status parallel11_status.txt, console parallel11_console.txt
python run_parallel_v11.py > parallel11_console.txt 2>&1
echo batch 11 finished %date% %time%
type parallel11_console.txt
pause
