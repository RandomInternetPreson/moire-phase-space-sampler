@echo off
REM RUN_FUZZ4.bat -- the disc's edge, the clean tests (6 Nov): the rim at the CENSUS'S OWN degree (Krim = K), which is what pays
REM   (degree-3 rims bought almost nothing: 1.647 / 1.412 against 1.373 for the degree-6 rim at K = 6).  census_fuzz_v3.py with the
REM   lighter positivity grid and tolerances, and a 48 h hard time limit per LP (--tlim 172800): if the solve caps, the script reports
REM   what it has and the cap is recorded as a negative.  Two launcher jobs, 40 GiB each.  Job 409 (L1 K8, Krim 3) keeps running.
REM Run from covariance_audit.  Needs: census_fuzz_v3.py, census_bandlimit_v4.py, part_nsh32_grid0_nang161_nw2048.npz,
REM   fuzz_contracted_Lmax3_part_nsh32_grid0_nang161_nw2048.npz.
REM Logs (written by the script): fuzz3_L2_K6_r6_atom.txt, fuzz3_L1_K8_r8_atom.txt.   Read:  findstr /C:"MD/floor" fuzz3_*.txt
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 40G --name fuzz4_L2K6r6 -- python "%~dp0census_fuzz_v3.py" part_nsh32_grid0_nang161_nw2048.npz --L 2 --K 6 --scope atom --eps 0.01 --Lmax 3 --Krim 6 --tlim 172800
python P:\SharedComputerResources\launch.py submit --owner claude --mem 40G --name fuzz4_L1K8r8 -- python "%~dp0census_fuzz_v3.py" part_nsh32_grid0_nang161_nw2048.npz --L 1 --K 8 --scope atom --eps 0.01 --Lmax 3 --Krim 8 --tlim 172800
echo.
python P:\SharedComputerResources\launch.py status
pause
