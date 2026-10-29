@echo off
REM RUN_TERRAIN.bat -- the grain test: two cut terrains (max(m,n) <= 3 and <= 2), lattice banks, LP, K-curve K=3..14.
REM Two launcher jobs of 60 GiB each (~4 h build + LP + ~4 h K-curve apiece); they run side by side. Run from covariance_audit.
REM Needs there: pipeline_terrain_v1.py, flight_time.py, greedy_results.json, dither_sweep_rho.py, census_bandlimit_v4.py,
REM   census_structure_v2.py, census_sumdiff_v1.py.  Each job works in its own subfolder terrain_ord3\ or terrain_ord2\.
REM Progress: type terrain_ord3\pipeline_terrain_ord3_console.txt   (and ord2).  The worker-check line near the top must say
REM   "workers see N terms (parent N)" with the same N -- that is the proof the cut terrain is the one being evolved.
cd /d "%~dp0"
python P:\SharedComputerResources\launch.py submit --owner claude --mem 60G --name terrain_ord3 -- python "%~dp0pipeline_terrain_v1.py" --mmax 3 --tag ord3 --nsh 32 --mode grid --seed 0 --nang 161 --nw 2048 --eps 0.01 --cores 28 --klist 3,4,5,6,8,10,12,14
python P:\SharedComputerResources\launch.py submit --owner claude --mem 60G --name terrain_ord2 -- python "%~dp0pipeline_terrain_v1.py" --mmax 2 --tag ord2 --nsh 32 --mode grid --seed 0 --nang 161 --nw 2048 --eps 0.01 --cores 28 --klist 3,4,5,6,8,10,12,14
echo.
python P:\SharedComputerResources\launch.py status
pause
