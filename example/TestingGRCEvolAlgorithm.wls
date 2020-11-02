#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
SetDirectory["/home/nesper/docs/biologia-de-sistemas/source"];
<<EvolAlgorithm4PatternFormingGRCModel2MathPackageV9.m
<<DesignMorphogeneResponsiveGRCs.m
(*<<graphGRN.m*)
(*<<ExpressionDynamics.m*)


(* ::Subsubsection:: *)
(*Run an instance of the evolutionary-inspired (MCMC) optimization algorithm used to engineered morephogen interpretation GRCs. *)
(*NOTE: This can take several hours before the algorithm converges upon a high-fitness solution with pattern forming capabilities. Check shell script: runLocalMCMClikeEvolGRCMathScriptFile.sh to run in batch mode on a Mac!*)


(* ::Input:: *)
Repl=2;
EvolSteps=20000;
ParamSettingSamplingRate=500;
thresholdFitness=0.7;
Quiet[EngineerMorphInducibleGRCGenotypes4SFGRM[ToExpression[Repl],ToExpression[EvolSteps],ToExpression[ParamSettingSamplingRate],ToExpression[thresholdFitness]]];


(* ::Subsubsection:: *)
(*Graph the resulting Gene Regulatory Networks*)


(* ::Input:: *)
(*(*Import resulting GRNs*)*)
(*grns=ReadList["GRCModelWithSumAndFilterCRIF/SampledMorphInducibleGRCGenotypes/ParamGenotypeGRC_Repl_1.dat"][[-3;;]];*)


(* ::Input:: *)
(*tb=Table[{GraphGRN[grns[[i,1]]],GraphGRNInteractionStrength[grns[[i,1]]]},{i, Length[grns]}];*)
(*tb=Prepend[tb,{"GRN Topology","GRN with interaction strength"}];*)
(*Grid[tb, Frame->All]*)


(* ::Subsubsection:: *)
(*Graph the GRN expression dynamics*)


(* ::Input:: *)
(*CreatSpaceTimePlots4GRC[*)
(*spatioTemporalExpressionProfile[grns[[2]]]*)
(*]*)


(* ::Input:: *)
(*SolveGRCDynsBeforeAfterMorphInput[#][[{1,-1}]]&/@grns[[All]]*)
