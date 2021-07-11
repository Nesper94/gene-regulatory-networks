#!/usr/bin/env wolframscript
(* ::Package:: *)

<<"../modules/EvolAlgorithm.m"
<<"../modules/importData.m"

grns=importGRNs[];

A0s:=RandomVariate[NormalDistribution[1,.3],100];

header={"A0","Fitness","GRN","Topology"};

table=Table[initialConc=A0s;
	{#, AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[SSInputMorphogen[#], grns[[i,1]]], i, grns[[i,2]]}&/@initialConc
	,{i, Length[grns]}];

table=Flatten[table, 1];

PrependTo[table, header];
Export["../data/inputNoiseFitness.tsv", table]
