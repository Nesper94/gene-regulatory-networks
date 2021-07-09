#!/usr/bin/env wolframscript
(* ::Package:: *)

<<"../modules/importData.m"
<<"../modules/EvolAlgorithm.m"


grns=importGRNs[];


(*Set diffusion parameters to 0 in all GRNs*)
Do[grns[[i,1,2]]={0,0,0}, {i, Length[grns]}];


zeroDiffusionFitness = Table[{AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[grns[[i,1]]] , i, grns[[i,2]]}, {i, Length[grns]}];


(*Add header to the table*)
PrependTo[zeroDiffusionFitness, {"Fitness","GRN","Topology"}];


Export["../data/zero-diffusion-fitness.tsv", zeroDiffusionFitness]
