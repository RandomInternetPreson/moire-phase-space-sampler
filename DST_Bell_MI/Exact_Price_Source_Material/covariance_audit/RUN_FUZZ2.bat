@echo off
REM RUN_FUZZ2.bat -- the disc's edge, diagnostics pass (2 Nov): re-solve the two atom-scope LPs that matter with census_fuzz_v2.py,
REM   which saves the solution, prints the per-atom rim (modulation depth; which approach phase the rim's lean follows) and a clipped
REM   upper bound on the price.  Uses the particle store and the Lmax=3 contraction already in the folder (job 405): no trajectories.
REM   Two launcher jobs (30 + 40 GiB); the second waits on the RAM ledger if 405 is still running.  ~11 h and ~15 h.
REM Run from covariance_audit.  Needs: census_fuzz_v2.py, census_bandlimit_v4.py, part_nsh32_grid0_nang161_nw2048.npz,
REM   fuzz_contracted_Lmax3_part_nsh32_grid0_nang161_nw2048.npz.   Logs (written by the script itself): fuzz2_L1_K6_atom.txt, fuzz2_L2_K6_atom.txt.
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 30G --name fuzz2_L1K6 -- python "%~dp0census_fuzz_v2.py" part_nsh32_grid0_nang161_nw2048.npz --L 1 --K 6 --scope atom --eps 0.01 --Lmax 3
python P:\SharedComputerResources\launch.py submit --owner claude --mem 40G --name fuzz2_L2K6 -- python "%~dp0census_fuzz_v2.py" part_nsh32_grid0_nang161_nw2048.npz --L 2 --K 6 --scope atom --eps 0.01 --Lmax 3
echo.
python P:\SharedComputerResources\launch.py status
pause
