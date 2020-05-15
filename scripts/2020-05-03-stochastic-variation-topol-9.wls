#!/usr/bin/env wolframscript -script
(* ::Package:: *)

(*Set directory and import modules*)
SetDirectory[FileNameJoin[{ParentDirectory[DirectoryName @ $InputFileName], "/modules"}]];
<<GRCPerturbutations`
<<importData`
<<dibujarGRN`
<<fitness2`


(*Import GRNs and select only those that present topology number 9*)
grns = importGRNs[];
Print["Topologies imported"]
topol9grns = Select[grns, #[[2]] == 9 &];
Print["Topology 9 selected"]


(* Stochastic variation of morphogen max concentration (A0) *)
fitnessStochasticA0 = Table[a0=RandomVariate[ NormalDistribution[1,.3] ];
 Max[ fitness[topol9grns[[#,1]],a0]]
, 1] & /@ Range[ Length[topol9grns] ];

(*Export list of fitness to fitness-perturb-a0.txt*)
Export[ FileNameJoin[{ ParentDirectory[ DirectoryName @ $InputFileName ], "results/perturbations-topol-9/fitness-perturb-a0.txt" }], fitnessStochasticA0]
Print["Exported fitnessStochasticA0 to results/perturbations-topol-9/fitness-perturb-a0.txt"]

(*Flatten list and calculate percentage of times fitness was \[GreaterEqual] than 0.95, this will be our measure of robustness*)
robustnessA0 = Count[ Flatten[ fitnessStochasticA0 ], u_/; u>=0.95]/Length[Flatten[fitnessStochasticA0]];
Print["Robustness to A0 stochastic fluctations: ",robustnessA0]


(* Stochastic variation of morphogen gradient decay parameter (h) *)
fitnessStochastich = Table[h = RandomVariate[ NormalDistribution[.4, .3*.4] ];
  Max[ fitness[topol9grns[[#, 1]], 1, h]]
 , 1] & /@ Range[ Length[topol9grns] ];

(*Export list of fitness to fitness-perturb-h.txt*)
Export[ FileNameJoin[{ ParentDirectory[ DirectoryName @ $InputFileName ], "results/perturbations-topol-9/fitness-perturb-h.txt" }], fitnessStochastich]
Print["Exported fitnessStochastich to results/perturbations-topol-9/fitness-perturb-h.txt"]

(*Flatten list and calculate percentage of times fitness was \[GreaterEqual] than 0.95, this will be our measure of robustness*)
robustnessh = Count[ Flatten[ fitnessStochastich ], u_/; u>=0.95]/Length[Flatten[fitnessStochastich]];
Print["Robustness to h stochastic fluctations: ",robustnessh]


(* Stochastic variation of Initial Concentrations of Gene Products (ICs)
Remember that default values are A0=1, h=0.4, ICs=ConstantArray[1,90] *)
fitnessStochasticICs = Table[ICs = Table[RandomVariate[NormalDistribution[1, .3]], 90];
   Max[ fitness[ topol9grns[[#, 1]], 1, .4, ICs ] ]
 , 1] & /@ Range[ Length[topol9grns] ];

(*Export list of fitness to fitness-perturb-h.txt*)
Export[ FileNameJoin[{ ParentDirectory[ DirectoryName @ $InputFileName ], "results/perturbations-topol-9/fitness-perturb-ics.txt" }], fitnessStochasticICs]
Print["Exported fitnessStochasticICs to results/perturbations-topol-9/fitness-perturb-ics.txt"]

(*Flatten list and calculate percentage of times fitness was \[GreaterEqual] than 0.95, this will be our measure of robustness*)
robustnessICs = Count[ Flatten[ fitnessStochasticICs ], u_/; u>=0.95]/Length[Flatten[fitnessStochasticICs]];
Print["Robustness to ICs stochastic fluctations: ",robustnessICs]
