@echo off
REM RUN_BATCH15.bat -- batch 15: the sorting rule on random0 and the 4096 census (census_bumps_v2.py, three vertices each), ~1 h.
REM Runs directly; fine alongside job 399. Run from covariance_audit. Afterwards run the beads_control_v1.py line in my message.
REM Needs: census_bumps_v2.py, census_bandlimit_v4.py, run_parallel_v15.py, rho_nsh32_random0 and rho_nsh64_random0 jsons + their banks.
cd /d "%~dp0"
echo batch 15 started %date% %time%  -- logs bumps_r0ctl_p and bumps_n64ctl_p^<seed^>.txt, status parallel15_status.txt, console parallel15_console.txt
python run_parallel_v15.py > parallel15_console.txt 2>&1
echo batch 15 finished %date% %time%
type parallel15_console.txt
pause
