#!/bin/sh

###To run application locally on a Mac

arg1=2
arg2=5000
arg3=500
arg4=0.75

/Applications/Mathematica2.app/Contents/MacOS/MathKernel math -noprompt -run "Repl=$arg1;EvolSteps=$arg2;ParamSettingSamplingRate=$arg3;thresholdFitness=$arg4;<<runMCMClikeEvolGRCMathScriptFile.m"