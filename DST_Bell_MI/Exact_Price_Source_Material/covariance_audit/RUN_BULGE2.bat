@echo off
REM RUN_BULGE2.bat -- bulge sweep 2 (1 Nov): the mesa's second-order geometry -- tilted cap, front/back walls (+ twin upside down),
REM   kinked walls, wall position following the approach phase (square vs circle).  12 readers, each: S_max, relaxed miss, full price,
REM   K = 6, 8, 10, marginals, joint MD.  Runs on the stored histograms from job 403 (no trajectories), 12 in parallel, ONE launcher job 60 GiB,
REM   ~2-5 h.  Run from covariance_audit.  Needs there: run_bulge2_v1.py, census_profile_v3.py, census_bandlimit_v4.py,
REM   hist_nsh32_grid0_nang161_nw2048.npz (+ banks_nsh32_grid0_nang161_nw2048.npz for the Q check).
REM Progress:  type bulge2_console.txt     per reader:  type bulge2_tilt0p05.txt  (bulge2_awall3p6_7p2, bulge2_mawall3p6_14p3, bulge2_bwall1p8_7p2, bulge2_shift10 ...)
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 60G --name bulge_sweep2 -- python "%~dp0run_bulge2_v1.py"
echo.
python P:\SharedComputerResources\launch.py status
echo.
echo Submitted. Progress:  type bulge2_console.txt
pause
