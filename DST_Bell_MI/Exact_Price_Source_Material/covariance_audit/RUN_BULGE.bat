@echo off
REM RUN_BULGE.bat -- the bulge test (31 Oct): rebuild the 32x32 lattice trajectories ONCE storing every lambda-histogram
REM   (hist_nsh32_grid0_nang161_nw2048.npz, ~0.5-1.5 GB), then price the census under 14 detection profiles (the paper's step
REM   as the control, the smooth bulge, four blunted steps, notched and anti-notched steps), each: S_max, relaxed miss,
REM   full price, K = 6,8,10,12,14, marginals.  ONE launcher job, 80 GiB: ~4 h build + ~4-8 h sweep (14 profiles in parallel).
REM Run from covariance_audit.  Needs there: pipeline_bulge_v1.py, census_profile_v1.py, census_bandlimit_v4.py,
REM   flight_time.py, dither_sweep_rho.py, greedy_results.json  (+ banks_nsh32_grid0_nang161_nw2048.npz for the bit-for-bit check).
REM Progress:  type pipeline_bulge_console.txt      per profile:  type bulge_sign.txt  (bulge_smooth, bulge_tanh1 ... bulge_anotch5_e10)
REM The line "against banks_...: max|Svec - banks| ... max|Ecoef - banks| ..." in bulge_sign.txt must read ~1e-15: that is the
REM   proof the rebuilt histograms are the deposited banks.  If the build is interrupted, rerun this file: the sweep skips finished profiles.
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 80G --name bulge_grid32 -- python "%~dp0pipeline_bulge_v1.py" --nsh 32 --mode grid --seed 0 --nang 161 --nw 2048 --eps 0.01 --cores 56 --klist 6,8,10,12,14 --maxpar 14
echo.
python P:\SharedComputerResources\launch.py status
echo.
echo Submitted. Progress:  type pipeline_bulge_console.txt
pause
