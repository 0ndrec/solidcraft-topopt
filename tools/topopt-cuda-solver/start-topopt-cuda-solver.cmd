@echo off
rem TopOpt CUDA Solver: optimizes on an NVIDIA GPU (and the CPU) for the TopOpt plugin in SolidCraft.
rem The folder must be the one granted to the plugin as "results" (Preferences > Plugins > TopOpt).
rem Default: the "results" folder next to this script's folder (build\results for build\solver).
rem Pass another folder as the first argument, or edit RESULTS below.
title TopOpt CUDA Solver
set "RESULTS=%~1"
if "%RESULTS%"=="" set "RESULTS=%~dp0..\results"
if not exist "%RESULTS%" mkdir "%RESULTS%"
"%~dp0topopt-cuda-solver.exe" watch "%RESULTS%"
pause
