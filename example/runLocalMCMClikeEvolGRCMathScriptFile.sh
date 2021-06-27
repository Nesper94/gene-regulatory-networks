#!/bin/sh

# To run application locally

arg1=2
arg2=5000
arg3=500
arg4=0.95

echo "Running with values:"
echo "  EvolSteps = $arg2"
echo "  ParamSettingSamplingRate = $arg3"
echo "  thresholdFitness = $arg4"

for NUM_NUCLEI in 10 20 40 50; do

echo "Running algorithm with NumNuclei = $NUM_NUCLEI"

#/Applications/Mathematica2.app/Contents/MacOS/MathKernel math -noprompt -run "Repl=$arg1;EvolSteps=$arg2;ParamSettingSamplingRate=$arg3;thresholdFitness=$arg4;<<runMCMClikeEvolGRCMathScriptFile.m"
math -noprompt -run "uNumNuclei=$NUM_NUCLEI;Repl=$arg1;EvolSteps=$arg2;ParamSettingSamplingRate=$arg3;thresholdFitness=$arg4;<<runMCMClikeEvolGRCMathScriptFile.m"

done
