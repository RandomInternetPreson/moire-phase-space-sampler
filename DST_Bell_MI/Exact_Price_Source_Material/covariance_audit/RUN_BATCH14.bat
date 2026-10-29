@echo off
REM RUN_BATCH14.bat -- batch 14: the partition control (census_bumps_v2.py, six perturbed vertices), ~20 min each, light.
REM Runs directly; fine alongside job 399. Run from covariance_audit. Afterwards run the beads_control_v1.py line in my message.
REM Needs: census_bumps_v2.py, census_bandlimit_v4.py, run_parallel_v14.py, rho_nsh32_grid0_eps0.01_nang161.json + banks.
cd /d "%~dp0"
echo batch 14 started %date% %time%  -- logs bumps_ctl_p^<seed^>.txt, status parallel14_status.txt, console parallel14_console.txt
python run_parallel_v14.py > parallel14_console.txt 2>&1
echo batch 14 finished %date% %time%
type parallel14_console.txt
pause
