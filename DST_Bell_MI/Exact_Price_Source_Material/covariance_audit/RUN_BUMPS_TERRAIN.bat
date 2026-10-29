@echo off
REM RUN_BUMPS_TERRAIN.bat -- the crest LP (K=6 + 6-degree bumps, 16x16) on the two CUT terrains, for comparison with batch 13's
REM   full-terrain run: does the bump count / position set track the mesh's beats?  Two light jobs (~30 min each), run directly.
REM Run from covariance_audit. Copies census_bumps_v1.py into terrain_ord3\ and terrain_ord2\ (census_bandlimit_v4.py and the
REM   banks are already there) and runs the LP inside each folder. Outputs: terrain_ord3\bumps_cr_ord3_k6_w6_nc16.txt (and ord2).
cd /d "%~dp0"
copy /Y census_bumps_v1.py terrain_ord3\ >nul
copy /Y census_bumps_v1.py terrain_ord2\ >nul
set K=6
set BUMP_W=6
set BUMP_NC=16
set EPS=0.01
set METHOD=highs-ipm
cd terrain_ord3
set OUT_TAG=cr_ord3_k6_w6_nc16
start "bumps ord3" /B python census_bumps_v1.py rho_nsh32_grid0_eps0.01_nang161.json > bumps_cr_ord3_k6_w6_nc16.txt 2>&1
cd ..\terrain_ord2
set OUT_TAG=cr_ord2_k6_w6_nc16
start "bumps ord2" /B python census_bumps_v1.py rho_nsh32_grid0_eps0.01_nang161.json > bumps_cr_ord2_k6_w6_nc16.txt 2>&1
cd ..
echo started both; check with:  findstr /C:"SUMMARY" terrain_ord3\bumps_cr_*.txt terrain_ord2\bumps_cr_*.txt
pause
