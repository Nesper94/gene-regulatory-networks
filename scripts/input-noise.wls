#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
<<"../modules/EvolAlgorithm.m"
<<"../modules/importData.m"

grns=importGRNs[];

A0s:=RandomVariate[NormalDistribution[1,.3],100];

header={{"Mean fitness","GRN","Topology"}};

fitnessTable=Table[
	{Mean[AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[SSInputMorphogen[#], grns[[i,1]]]&/@A0s],
	i, grns[[i,2]]}
	,{i,Length[grns]}];
	
table=Join[header,fitnessTable];
Export["../data/inputNoiseFitness.tsv", table]
