(* ::Package:: *)

Needs["parameters`"]

(*USE THE FOLLOWING METHODS TO PROBE THE DESIGN SPACE OF 3-NODE MORPHOGENE RESPONSIVE GRCs*)
(************************************************************************************************************************)
(************************************************************************************************************************)

PatternFilterChecking[ExpProfiles_List]:=(1/Length[ExpProfiles])(EuclideanDistance[ExpProfiles,ConstantArray[Mean[ExpProfiles],Length[ExpProfiles]]])

Q1[x_]:=(x^5)/( (x^5) + (0.1^5) )

(************************************************************************************************************************)

AssessExpPattern4SingleStripeFeat[ExpProfiles_List]:=Block[{},

(*Optimal expression pattern represented in a scale between 1-10 in expression level, being 10 an expression level which is > 90% of the maximal level observed along the 1D field of cells for the output genes*)
OptimalPattern=Flatten[{ConstantArray[1,10],ConstantArray[10,10],ConstantArray[1,10]}];
(*Maximal discrepancy achievable for a given pattern w.r.t to the optimal pattern above*)
Dmax=9*30;

ThresholdedExpValues={{0.`,0.1`},{0.1`,0.2`},{0.2`,0.3`},{0.3`,0.4`},{0.4`,0.5`},{0.5`,0.6`},{0.6`,0.7`},{0.7`,0.8`},{0.8`,0.9`},{0.9`,1.001`}};

testInt[Rng_,Value_]:= Rng[[1]]<=Value< Rng[[2]];

NormEquilbPattern=ExpProfiles/Max[ExpProfiles];

DiscretizedPattern=Table[Position[testInt[#,NormEquilbPattern[[x]]]&/@ThresholdedExpValues,True][[1,1]],{x,Length[NormEquilbPattern]}];

PFeff = 1-(ManhattanDistance[OptimalPattern,DiscretizedPattern]/Dmax)
]

(************************************************************************************************************************)

FitnessScore4StripePattern[TimeSeriesSpatialExpOutput_,EndPointExpPatterns_]:=Block[{MaxGeneExpVarianceAcrossLastTimeWindow,SSThrValue,SSCondition,PF,FS},

(*Testing GRC dynamics for stationarity. Threshold value <= 10^-6*)
MaxGeneExpVarianceAcrossLastTimeWindow=Max[AssessTargetGeneExpVarianceAcrossLastTimeWindow/@TimeSeriesSpatialExpOutput];

SSThrValue=10^-6;
SSCondition=If[MaxGeneExpVarianceAcrossLastTimeWindow<= SSThrValue,1,0];

(*Testing for spatial heterogeneity in the steady state expression pattern, which is also an indicative of whether the expression profiles are sufficiently high so that the cross regulatory interactions can be effective in controlling the expression putput.*)
PF=PatternFilterChecking/@EndPointExpPatterns;

FS=Max[(AssessExpPattern4SingleStripeFeat[#]*Q1[Min[PF]]*SSCondition)&/@EndPointExpPatterns]
]

(************************************************************************************************************************)
(*2019.03.07: En la siguiente función modifiqué los argumentos ICs y MorphInput de tal manera que tuvieran valores por defecto*)
AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs_List:ConstantArray[1, 3*NumNuclei],MorphInput_List:SSInputMorphogen[1],{Wmatrix_List,DiffParams_List,DegParams_List}]:=Block[
{TimeSeriesSpatialExpOutput,EndPointExpPatterns,SSExpValues,MaxFS},

TimeSeriesSpatialExpOutput = StripeFormingGRCs4SSMorpGradientSFGRM[ICs,MorphInput,{Wmatrix,DiffParams,DegParams}];



EndPointExpPatterns=Last/@TimeSeriesSpatialExpOutput;

SSExpValues=Flatten[Thread[EndPointExpPatterns]];

(*Check if there is any negative gene expression value, which is not biologically realistic. If there is, then Fitness must be 0. *)

MaxFS=If[
	Count[Flatten[TimeSeriesSpatialExpOutput],x_/;x<=0]>0,0,
	FitnessScore4StripePattern[TimeSeriesSpatialExpOutput,EndPointExpPatterns]
	];

{MaxFS,SSExpValues}
		]

(************************************************************************************************************************)

AssessExpDistance[ExpProfile_,EndPointExpProfile_]:=(1/Length[ExpProfile])(EuclideanDistance[ExpProfile,EndPointExpProfile])

AssessTargetGeneExpVarianceAcrossLastTimeWindow[TimeSeriesTargetGeneSpatialExp_]:=Block[{NormalizedTimeSeriesTargetGeneSpatialExp,NucleiWiseAvgGeneExpAcrossInterval},
NormalizedTimeSeriesTargetGeneSpatialExp=TimeSeriesTargetGeneSpatialExp/Max[Flatten[TimeSeriesTargetGeneSpatialExp]];
(*Computing average expression across interval*)
NucleiWiseAvgGeneExpAcrossInterval=Mean/@Transpose[NormalizedTimeSeriesTargetGeneSpatialExp];
(*Computing mean deviation across time interval with respect to average expression: this is analogous to a variance over the last t time steps simulated. This metric is borrowed from Pujato etal, 2013. Values <= 0.0001 can be taken as indicative of stationarity*)
Mean[AssessExpDistance[#,NucleiWiseAvgGeneExpAcrossInterval]&/@NormalizedTimeSeriesTargetGeneSpatialExp]
]

(************************************************************************************************************************)

StripeFormingGRCs4SSMorpGradientSFGRM[ICs_List,MorphInput_List,{Wmatrix_List,DiffParams_List,DegParams_List}]:=Module[
{Prot,ProtStateVariab,ProtInitExpStat,DynSyst,Sol,GRCPhenotReadout,IntTime,sigmoidThreshold,alpha,NumNuclei,GRNSize},

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*Transfer function for mapping integrated regulatory inputs to transcriptional outputs. In this model the SigmoidSteepness = alpha = 5, while the parameter b, which gives the location of the threshold value, is set to 1 *)
SumAndFilterF[SigmoidSteepness_,Threshold_,IntegratedRegInput_]:=1./(1.+Exp[(SigmoidSteepness- (Threshold*IntegratedRegInput))]);
(*A function for summing the weigthed regulatory inputs*)
AdditiveRegContributionF[a_,n_,t_,W_List,Morph_]:=Total[W[[a]]*Prepend[Table[Prot[b,n,t],{b,1,GRNSize}],Morph]];
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)

GRNSize = Length[Wmatrix];

(*Initial conditions for all variables, including boundary conditions set = 0*)
ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];
ProtInitExpStat=Flatten[{Thread[Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,NumNuclei}]]==ICs],Flatten[Table[Table[Prot[a,n,0]==0.0,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

DynSyst=Flatten[
{Table[
Table[D[Prot[a,n,t],t]==
( ( SumAndFilterF[alpha,sigmoidThreshold,AdditiveRegContributionF[a,n,t,Wmatrix,MorphInput[[n]]]] )+ ( DiffParams[[a]]( Prot[a,n-1,t]+Prot[a,n+1,t]-2*Prot[a,n,t] ) ) - ( DegParams[[a]]*Prot[a,n,t] )),{a,1,GRNSize}],
{n,NumNuclei}],Table[Table[D[Prot[a,n,t],t]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]}
];

Sol=NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

GRCPhenotReadout= Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,IntTime-100,IntTime,2}],{i,3}]
]

(************************************************************************************************************************)

AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[EvalGRCParamGenotype_List]:=Block[{ICs,NullMorpInput,MorpInputProfile,PreMorpInputFS,SSExpValuesPreMorpInput,FS4GRCInResponse2MorpInput},
(*Set ICs for all genes in all nuclei and run system without Morphogene input*)
ICs=ConstantArray[1, 3*NumNuclei];
NullMorpInput=ConstantArray[0, NumNuclei];
MorpInputProfile=SSInputMorphogen[1.0];

{PreMorpInputFS,SSExpValuesPreMorpInput}=AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,NullMorpInput,EvalGRCParamGenotype];
(*Set the pre-Morphogen SS expression levels as ICs for simulating the dynamics of the system in the presence of the Morphogene*)

(*Check that fitness before applying the Morphogene input is >0 (reaches stationary behavior) and is <= 0.5 (is unable to generate any striped pattern at all), if this condition is met then assess the ability of the system to interpret the input and generate a striped pattern for any of the genes*)
FS4GRCInResponse2MorpInput=If[
0<PreMorpInputFS<=0.5,
AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[SSExpValuesPreMorpInput,MorpInputProfile,EvalGRCParamGenotype][[1]],0
			]
]

(************************************************************************************************************************)

Sol4StripeFormingGRCs4SSMorpGradientSFGRM[EvalGRCParamGenotype_List]:=Block[{ICs,NullMorpInput,MorpInputProfile,EndPointExpPatterns,SSExpValues,TimeSeriesSpatialExpOutput},

ICs=ConstantArray[1, 3*NumNuclei];
NullMorpInput=ConstantArray[0, NumNuclei];
MorpInputProfile=SSInputMorphogen[1.0];

EndPointExpPatterns=Last/@StripeFormingGRCs4SSMorpGradientSFGRM[ICs,NullMorpInput,EvalGRCParamGenotype];

SSExpValues=Flatten[Thread[EndPointExpPatterns]];

TimeSeriesSpatialExpOutput=StripeFormingGRCs4SSMorpGradientSFGRM[SSExpValues,MorpInputProfile,EvalGRCParamGenotype]
]


(************************************************************************************************************************)

(*NOTE: to assess the degenerate structure of a GRC we have to slightly modify the previous methods used to assess a GRC
 for its ability to	generate a stripe expression pattern*)
AssessDegenerateStruc4EngineeredMorpResponsiveGRC[GRCGenotype_List,
  FitnessThrs_, MaxNumDeletedEdges_, ReportData_String] :=
 Block[{ParamGRCWiring, MorphogeneInput, DiffParams, DegParams,
   GeneGeneRegEdgesIDs, EdgeDeletionCombsList, EdgeDeletedGRCWsList,
   FocalGRCFitness, MutantGRCFitnessList, ViableWiringDesigns,
   BestPerformingWiringDesignsList},

  {ParamGRCWiring, DiffParams, DegParams} = GRCGenotype;

  GeneGeneRegEdgesIDs =
   Flatten[Position[Flatten[ParamGRCWiring], x_ /; x != 0]];

  (*EdgeDeletionCombsList=Cases[Subsets[GeneGeneRegEdgesIDs],x_List/;
  0<Length[x]\[LessEqual](Length[GeneGeneRegEdgesIDs]-4)];*)

  EdgeDeletionCombsList =
   Cases[Subsets[GeneGeneRegEdgesIDs],
    x_List /; 0 < Length[x] <= MaxNumDeletedEdges];

  EdgeDeletedGRCWsList =
   Partition[
      ReplacePart[Flatten[ParamGRCWiring],
       Thread[# -> ConstantArray[0, Length[#]]]], 4] & /@
    EdgeDeletionCombsList;

  FocalGRCFitness =
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[{\
ParamGRCWiring, DiffParams, DegParams}];

  MutantGRCFitnessList =
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM2[{\
#, DiffParams, DegParams}] & /@ EdgeDeletedGRCWsList;

  ViableWiringDesigns =
   Flatten[Position[MutantGRCFitnessList, x_ /; x >= FitnessThrs]];

  BestPerformingWiringDesignsList =
   Reverse[SortBy[
      Thread[{EdgeDeletedGRCWsList[[ViableWiringDesigns]],
        MutantGRCFitnessList[[ViableWiringDesigns]]}], Last]][[All,
     1]];

  If[Length[BestPerformingWiringDesignsList] > 0,
   If[ReportData ==
     "AllFeats", {Count[(Flatten[
         BestPerformingWiringDesignsList[[1]]] -
        Flatten[ParamGRCWiring]),
      x_ /; x != 0], {BestPerformingWiringDesignsList[[1]],
      DiffParams,
      DegParams}, {DisplayNonWeightedHaploidGraphWithMorphInputs@
       ParamGRCWiring,
      DisplayNonWeightedHaploidGraphWithMorphInputs@
       BestPerformingWiringDesignsList[[1]]}}, {Count[(Flatten[
         BestPerformingWiringDesignsList[[1]]] -
        Flatten[ParamGRCWiring]),
      x_ /; x != 0], {BestPerformingWiringDesignsList[[1]],
      DiffParams, DegParams}}], {}]

  ]


AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM2[
  EvalGRCParamGenotype_List] :=
 Block[{ICs, NullMorpInput, MorpInputProfile, PreMorpInputFS,
   SSExpValuesPreMorpInput, FS4GRCInResponse2MorpInput},
  (*Set ICs for all genes in all nuclei and run system without \
Morphogene input*)
  ICs = ConstantArray[1, 3*NumNuclei];
  NullMorpInput = ConstantArray[0, NumNuclei];
  MorpInputProfile = SSInputMorphogen[1.0];

  {PreMorpInputFS, SSExpValuesPreMorpInput} =
   AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,
    NullMorpInput, EvalGRCParamGenotype];

  (*Set the pre-
  Morphogen SS expression levels as ICs for simulating the dynamics \
of the system in the presence of the Morphogene*)(*Check that fitness \
before applying the Morphogene input is>
  0 (reaches stationary behavior) and is\[LessEqual]0.5 (is unable to \
generate any striped pattern at all),
  if this condition is met then assess the ability of the system to \
interpret the input and generate a striped pattern for any of the \
genes*)FS4GRCInResponse2MorpInput =
   If[0 < PreMorpInputFS <= 0.5,
    AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM2[
      SSExpValuesPreMorpInput, MorpInputProfile,
      EvalGRCParamGenotype][[1]], 0]]

AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM2[ICs_List,
  MorphInput_List, {Wmatrix_List, DiffParams_List, DegParams_List}] :=
  Block[{TimeSeriesSpatialExpOutput, EndPointExpPatterns, SSExpValues,
    MaxFS},

  TimeSeriesSpatialExpOutput =
   StripeFormingGRCs4SSMorpGradientSFGRM[ICs,
    MorphInput, {Wmatrix, DiffParams, DegParams}];

  EndPointExpPatterns = Last /@ TimeSeriesSpatialExpOutput;

  SSExpValues = Flatten[Thread[EndPointExpPatterns]];

  (*Check if there is any negative gene expression value,
  which is not biologically realistic.If there is,
  then Fitness must be 0.*)
  MaxFS = If[
    Count[Flatten[TimeSeriesSpatialExpOutput], x_ /; x <= 0] > 0, 0,
    FitnessScore4StripePattern2[TimeSeriesSpatialExpOutput,
     EndPointExpPatterns]];
  {MaxFS, SSExpValues}]


  FitnessScore4StripePattern2[TimeSeriesSpatialExpOutput_,
  EndPointExpPatterns_] :=
 Block[{MaxGeneExpVarianceAcrossLastTimeWindow, SSThrValue,
   SSCondition, PF, FS},

  (*Testing GRC dynamics for stationarity.Threshold \
value\[LessEqual]10^-6*)
  MaxGeneExpVarianceAcrossLastTimeWindow =
   Max[AssessTargetGeneExpVarianceAcrossLastTimeWindow /@
     TimeSeriesSpatialExpOutput];

  SSThrValue = 10^-6;

  SSCondition =
   If[MaxGeneExpVarianceAcrossLastTimeWindow <= SSThrValue, 1, 0];

  (*Testing for spatial heterogeneity in the steady state expression \
pattern,which is also an indicative of whether the expression \
profiles are sufficiently high so that the cross regulatory \
interactions can be effective in controlling the expression putput.*)
  PF = PatternFilterChecking /@ EndPointExpPatterns;
  FS = Max[
    MapIndexed[(AssessExpPattern4SingleStripeFeat[#]*
        Q1[PF[[#2[[1]]]]]*SSCondition) &, EndPointExpPatterns]]]


(************************************************************************************************************************)

(*Engineer morphogene responsive GRCs using a MCMC-like algorithm*)
EngineerMorphInducibleGRCGenotypes4SFGRM[Repl_,PSR_,ThrFS_]:=Block[{EvolSteps,thresholdFitness,ParamSettingSamplingRate,
MutaProbList,GRCWiring,MorphogenInput,DiffRates,ProtDegRates,GetInputParamSet,GRCDyns,RefInputParamSetting,Fb,
Fa,DeltaF,ThisMutId,InputParamSetting,InputParamSettingAllin,OutputFileName,PatternFormingSeeds,DestinationFolder,
ParamsSamplingCond,SamplingPointsList,BoolWMatrix,WContainer},

EvolSteps=5000000;
ParamSettingSamplingRate=PSR;(*e.g. PSR = 500*)
thresholdFitness=ThrFS;(*e.g. ThrFS = 0.925*)
SamplingPointsList=Drop[Range[0,EvolSteps,ParamSettingSamplingRate],{1}];

(*Container to append high-fitness-scoring, topologically distinct wirings as they are discovered*)
WContainer = {};

(*Only mutate the wiring (with the highest probab), as well as diffusion and protein degradation rates*)
MutId:= RandomSample[{0.7,0.15,0.15}-> {1,2,3},1][[1]];

(*Select number of mutations to perform between 1-10, acording to a power-law function*)
MutaProbList = (0.5^#)&/@Range[10];
NumMutants:= RandomSample[MutaProbList-> Range[10],1][[1]];

(***********************************************************************)
(***********************************************************************)
DestinationFolder = "GRCModelWithSumAndFilterCRIF/SampledMorphInducibleGRCGenotypes";
(*To ensure that the destination folders exist, if don't then create them!*)

If[
	FileExistsQ[DestinationFolder]==False,
	CreateDirectory[DestinationFolder]
	];

(*To read in within mathematica these files do: ReadList[FitnessFeaturesFile]*)
OutputFileName = DestinationFolder<>"/ParamGenotypeGRC_Repl_"<>ToString[Repl]<>".dat";
(***********************************************************************)
(***********************************************************************)


(*PatternFormingSeeds=DeleteDuplicates[ReadList["GRCModelWithSumAndFilterCRIF/PatternFormingSeeds4SFGRModelExploration.dat"]];*)

(*Start exploration of the design space of GRCs from randomly chosen parameters/genotypes*)
{GRCWiring,DiffRates,ProtDegRates} = {GetRndComposedWiringMatrix2,GetDiffRates[3],GetProteinDegrad[3]};

RefInputParamSetting={GRCWiring,DiffRates,ProtDegRates};

Fb=AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[RefInputParamSetting];

(*Optimization cycle to explore design space of GRCs*)
Do[

Clear[InputParamSetting,InputParamSettingAllin,GRCDyns,Fa,DeltaF];

ThisMutId = MutId;

InputParamSetting=Nest[MutatorF2[ThisMutId,#]&,RefInputParamSetting,NumMutants];

Fa=AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[InputParamSetting];

DeltaF = Fa-Fb;

(*
If[DeltaF>= 0.001,
Print[{Fa,iter}]
];
*)

(*Implementing a Metropolis algorithm*)
If[
	DeltaF>=0,
	RefInputParamSetting=InputParamSetting;
	Fb=Fa,
	If[DeltaF<0,
               If[Exp[DeltaF/0.005]>= Random[],
	           RefInputParamSetting=InputParamSetting;
	            Fb=Fa,
	           RefInputParamSetting=RefInputParamSetting;
	           Fb=Fb]
	              ]
			];

(*assessing if current Boolean wiring was already tested, if not then save to file*)
BoolWMatrix=Sign/@RefInputParamSetting[[1]];
(*Checking sampling times*)
ParamsSamplingCond = MemberQ[SamplingPointsList,iter];

(*Save current parameter setting if the following 3 conditions are met:
1)not previously seen BoolW, 2) Fitness >= thresholdFitness, 3) Sampling time*)
If[
    (MemberQ[WContainer,BoolWMatrix]==False),
         If[Fb>= thresholdFitness,
	  If[ParamsSamplingCond,
		(PutAppend[RefInputParamSetting,OutputFileName];
                   AppendTo[WContainer,BoolWMatrix])]
			]
				]

,{iter,0,EvolSteps}]
]

(************************************************************************************************************************)

(*Use this function to inspect the dynamics of an engineered GRC before and after applying the morphogene input*)
SolveGRCDynsBeforeAfterMorphInput[TestGRCParamGenotype_]:=Block[{ICs,NullMorpInput,PreMorpInputFS,SSExpValuesPreMorpInput,PatternExpBeforeMorpInput,GRCW,P1,P2,MorpInducedExpPattern},
(*Set ICs to 1 for all genes in all nuclei and run system without Morphogene input*)
ICs=ConstantArray[1, 3*NumNuclei];
NullMorpInput=ConstantArray[0, NumNuclei];

{PreMorpInputFS,SSExpValuesPreMorpInput} = AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,NullMorpInput,TestGRCParamGenotype];
(*Set the pre-Morphogen SS expression levels as ICs for simulating the dynamics of the system in the presence of the Morphogene*)

(*Before applying the morphogene input profile*)
PatternExpBeforeMorpInput =StripeFormingGRCs4SSMorpGradientSFGRM[ICs,NullMorpInput,TestGRCParamGenotype];

P1=MapIndexed[ListPlot[Last[#],Joined->True,ImageSize->200,PlotStyle->{Thickness[0.015],{Red,Blue,Green}[[#2[[1]]]]},Frame->True,FrameLabel->{Style["Nuclei",9],Style["Expression level",9]}]&,PatternExpBeforeMorpInput];

(*Upon applying the morphogene input profile*)
MorpInducedExpPattern =StripeFormingGRCs4SSMorpGradientSFGRM[SSExpValuesPreMorpInput,SSInputMorphogen[1.0],TestGRCParamGenotype];

P2=MapIndexed[ListPlot[Last[#],Joined->True,ImageSize->200,PlotStyle->{Thickness[0.015],{Red,Blue,Green}[[#2[[1]]]]},Frame->True,FrameLabel->{Style["Nuclei",9],Style["Expression level",9]}]&,MorpInducedExpPattern];

GRCW=Show[DisplayNonWeightedHaploidGraphWithMorphInputs@First[TestGRCParamGenotype],ImageSize->200];

{GRCW,P1,P2}
]

(*return a long series spatial expression profiles for each gene at different time points*)
ReturnTimeSeriesStripeFormingGRCs4SSMorpGradientSFGRM[ICs_List,MorphInput_List,{Wmatrix_List,DiffParams_List,DegParams_List}]:=Module[
{Prot,ProtStateVariab,ProtInitExpStat,DynSyst,Sol,GRCPhenotReadout,IntTime,sigmoidThreshold,alpha,NumNuclei,GRNSize},

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*Transfer function for mapping integrated regulatory inputs to transcriptional outputs. In this model the SigmoidSteepness = alpha = 5, while the parameter b, which gives the location of the threshold value, is set to 1 *)
SumAndFilterF[SigmoidSteepness_,Threshold_,IntegratedRegInput_]:=1./(1.+Exp[(SigmoidSteepness- (Threshold*IntegratedRegInput))]);
(*A function for summing the weigthed regulatory inputs*)
AdditiveRegContributionF[a_,n_,t_,W_List,Morph_]:=Total[W[[a]]*Prepend[Table[Prot[b,n,t],{b,1,GRNSize}],Morph]];
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)

GRNSize = Length[Wmatrix];

(*Initial conditions for all variables, including boundary conditions set = 0*)
ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];
ProtInitExpStat=Flatten[{Thread[Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,NumNuclei}]]==ICs],Flatten[Table[Table[Prot[a,n,0]==0.0,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

DynSyst=Flatten[
{Table[
Table[D[Prot[a,n,t],t]==
( ( SumAndFilterF[alpha,sigmoidThreshold,AdditiveRegContributionF[a,n,t,Wmatrix,MorphInput[[n]]]] )+ ( DiffParams[[a]]( Prot[a,n-1,t]+Prot[a,n+1,t]-2*Prot[a,n,t] ) ) - ( DegParams[[a]]*Prot[a,n,t] )),{a,1,GRNSize}],
{n,NumNuclei}],Table[Table[D[Prot[a,n,t],t]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]}
];

Sol=NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

GRCPhenotReadout= Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,0,IntTime,2}],{i,3}]
]

(*Use this function to inspect the dynamics of a GRC genotype when perturbed*)
InspectGRCExpressionDynamics4SSMorpGradientSFGRM[EvalGRCParamGenotype_List]:=Block[{ICs,NullMorpInput,MorpInputProfile,PreMorpInputFS,SSExpValuesPreMorpInput,FS4GRCInResponse2MorpInput},
(*Set ICs for all genes in all nuclei and run system without Morphogene input*)
ICs=ConstantArray[1, 3*NumNuclei];
NullMorpInput=ConstantArray[0, NumNuclei];
MorpInputProfile=SSInputMorphogen[1.0];

{PreMorpInputFS,SSExpValuesPreMorpInput}=AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,NullMorpInput,EvalGRCParamGenotype];
(*Set the pre-Morphogen SS expression levels as ICs for simulating the dynamics of the system in the presence of the Morphogene*)

(*Check that fitness before applying the Morphogene input is >0 (reaches stationary behavior) and is <= 0.5 (is unable to generate any striped pattern at all), if this condition is met then assess the ability of the system to interpret the input and generate a striped pattern for any of the genes*)
FS4GRCInResponse2MorpInput=If[
0<PreMorpInputFS<=0.5,
ReturnTimeSeriesStripeFormingGRCs4SSMorpGradientSFGRM[SSExpValuesPreMorpInput,MorpInputProfile,EvalGRCParamGenotype],
0
	]
]
