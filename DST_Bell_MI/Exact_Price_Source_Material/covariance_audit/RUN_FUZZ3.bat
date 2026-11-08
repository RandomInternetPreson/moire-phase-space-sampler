@echo off
REM RUN_FUZZ3.bat -- the disc's edge, tractable pass (5 Nov): the three atom-scope LPs that have been crawling in jobs 405/407, re-solved
REM   with census_fuzz_v3.py -- rim harmonics band-limited to K = 3 (--Krim 3), positivity at 16 points per harmonic order (the clipped
REM   census gives the matching upper bound), HiGHS tolerances 1e-7, dual-simplex retry.  Uses the particle store and the Lmax=3
REM   contraction already in the folder; no trajectories.  Three launcher jobs, 30 GiB each; expected hours, not days.
REM Run from covariance_audit.  Needs: census_fuzz_v3.py, census_bandlimit_v4.py, part_nsh32_grid0_nang161_nw2048.npz,
REM   fuzz_contracted_Lmax3_part_nsh32_grid0_nang161_nw2048.npz.
REM Logs (written by the script): fuzz3_L2_K6_r3_atom.txt, fuzz3_L1_K8_r3_atom.txt, fuzz3_L2_K8_r3_atom.txt.
REM Read the SUMMARY line of each: "MD/floor x.xxxx ... clipped upper bound y.yyyy" -- if x and y agree, the price is bracketed and done.
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 30G --name fuzz3_L2K6 -- python "%~dp0census_fuzz_v3.py" part_nsh32_grid0_nang161_nw2048.npz --L 2 --K 6 --scope atom --eps 0.01 --Lmax 3 --Krim 3
python P:\SharedComputerResources\launch.py submit --owner claude --mem 30G --name fuzz3_L1K8 -- python "%~dp0census_fuzz_v3.py" part_nsh32_grid0_nang161_nw2048.npz --L 1 --K 8 --scope atom --eps 0.01 --Lmax 3 --Krim 3
python P:\SharedComputerResources\launch.py submit --owner claude --mem 30G --name fuzz3_L2K8 -- python "%~dp0census_fuzz_v3.py" part_nsh32_grid0_nang161_nw2048.npz --L 2 --K 8 --scope atom --eps 0.01 --Lmax 3 --Krim 3
echo.
python P:\SharedComputerResources\launch.py status
pause
