@echo off
REM RUN_FUZZ.bat -- the disc's edge (1 Nov): rebuild the 32x32 lattice trajectories ONCE keeping every particle's final bin
REM   (part_nsh32_grid0_nang161_nw2048.npz, 5.2 GB in memory, ~1-3 GB on disk; ~7.5 h on 56 cores, same trajectories as job 403 --
REM   checked bit for bit against hist_nsh32_grid0_nang161_nw2048.npz), contract once (Lmax = 3), then the joint-birth LPs:
REM   L = 0 (the paper), 1, 2 lambda_0 harmonics  x  K = 6, 8, 10  x  scope cell / atom  -> 15 LPs, 6 at a time, fuzz_L<L>_K<K>_<scope>.txt.
REM ONE launcher job, 100 GiB.  Run from covariance_audit.  Needs there: pipeline_fuzz_v1.py, census_fuzz_v1.py, census_bandlimit_v4.py,
REM   flight_time.py, dither_sweep_rho.py, greedy_results.json  (+ hist_nsh32_grid0_nang161_nw2048.npz for the check).
REM Progress:  type pipeline_fuzz_console.txt      the line "check against hist_...: 0 (cell, setting) blocks differ" is the proof the
REM   particle store is the deposited dynamics.  L=0 rows of the final table must read 1.700 / 1.429 / 1.258.
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 100G --name fuzz_grid32 -- python "%~dp0pipeline_fuzz_v1.py" --nsh 32 --mode grid --seed 0 --nang 161 --nw 2048 --eps 0.01 --cores 56 --Llist 0,1,2 --Lmax 3 --klist 6,8,10 --scopes cell,atom --maxpar 6
echo.
python P:\SharedComputerResources\launch.py status
echo.
echo Submitted. Progress:  type pipeline_fuzz_console.txt
pause
