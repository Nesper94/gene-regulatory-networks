(* ::Package:: *)

BeginPackage["EvolAlgorithm4PatternFormingGRCModel2MathPackageV9`",{"ErrorBarPlots`"}]

(*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Package developed by Jayson Guti\[EAcute]rrez on 15/9/2015%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*)


SSInputMorphogen::usage = "Generate input morphogen profile";

GetRndComposedWiringMatrix::usage = "generate wiring";

GetRndComposedWiringMatrix2::usage = "generate wiring";

GetProteinDegrad::usage = "generate protein degradation";

GetDiffRates::usage = "generate diffusion";

DisplayNonWeightedHaploidGraphWithMorphInputs::usage = "Function for displaying non-weigthed regulatory wirings (GRC design templates), including the input Morphogen";

DisplayWeightedHaploidGraphWithMorphInputs::usage = "Function for displaying weigthed regulatory wirings, including the input Morphogen";

SpaceTimeFeats4StripeFormingGRCs4SSMorpGradientWithBiophyGroundedCRIF0::usage = "Function to assess spatio-temporal expression portraits of the GRCs implementing the thermodynamically-grounded CRIF";

SpaceTimeFeats4StripeFormingGRCs4SSMorpGradientWithSumAndFilterModel0::usage = "Function to assess spatio-temporal expression portraits of the GRCs implementing the Sum&Filter CRIF";

CreatSpaceTimePlots4GRC::usage="Function to create the spatio-temporal expression portraits of the GRCs";

EquilibFilter::usage = "This function is to assess whether the expression profile of the output node in the GRC models
						reaches a quasi-steady state. GRC dynamics are simulated over 500 time steps, and the spatial
						profile over the total number of nuclei considered attained at time step = 500 is compared to
						that attained previously at time step = 250. The difference between expression profiles is
						assessed based of the normalized Euclidean Distance. We focus on GRCs capable of achieving a
						relatively fast developmental patterning task (single stripe formation) because:
						1) computing time restrictions; and 2) the speed of developmental patterning is likely to be
						critical during embryogenesis (potential fitness correlate).
						Ideally, the value of this score should be \[LessEqual] 0.001, which is indicative that the expression
						profile at time tfinal has reached steady state";

PatternFilter::usage = "A patterning score describes how much \[OpenCurlyQuote]pattern\[CloseCurlyQuote] or spatial heterogeneity there is. Ideally, a value
						for this score should be \[GreaterEqual] 6, so that one can assure that there is sufficient spatial heterogeneity,
						and that the expression profiles are appreciably high across the field of nuclei considered so that
						genes can exert effective regulatory control on their targets";

Q1::usage = "Function to assess the quality of a given phenotypic function based on the pattern filter score ";

Q2::usage = "Function to assess the quality of a given phenotypic function based on the stability filter score";

AssessExpPattern4Stripe::usage = "In this function: the end point expression profile for the output gene is normalized
								  and then each expression value within a cell i is discretized so that it is assigned
								  a value on [1,10]. Based on this discretized list of expression values we compare
								  this with the ideal target pattern where the expression level of the output gene
								  within each cell lying in the middle part of the morphogen gradient (i=11 to i=20)
								  must attain > 90% of the maximal expression level observed along the 1D field of
								  cells";

FitnessF2::usage = "Same as FitnessF, but in this case the pattern filter is scored just for the output node";

FitnessF3::usage = "Similar to FitnessF2, but in this case we consider only the Q1[Min[PF]]";

AssessFitnessScoreStripeFormingGRCs4SSMorpGradient::usage = "This function implements a function that scores the fitness of a
															GRC based on the expression pattern of any gene in the system";

ComputeSigmoidCRIFWithLeakTerm::usage = "thermodynamically grounded cis-regulatory input function (CRIF), which is built upon principles described by Sherman and Cohen (2012) PloS Comput Biol";

StripeFormingGRCs4SSMorpGradientWithBiophyGroundedCRIF0::usage = "GRC model implementing the new CRIF function of gene regulation";

StripeFormingGRCs4SSMorpGradientWithSumAndFilterGeneRegModel::usage = "GRC model implementing the Sum & Filter sigmoidal transfer function for gene regulation";

Mutator4GRCW::usage = "Use this function to mutate the W matrix for a thermodynamically-grounded CRIF model. Note that a mutation can be either a change in parameter values (with the highest probab) or removal of one existing interaction";

Mutator4GRCW2::usage = "Use this function to mutate the W matrix for the Sum & filter model. Note that a mutation can be either a change in parameter values (with the highest probab) or removal of one existing interaction"

MutateParamVector::usage = "Use this function to mutate a parameter setting for Diffusion or protein degradation rates";

MutatorF::usage = "Use this function to mutate a thermodynamically-grounded CRIF GRC model in any of its component parameter settings";

MutatorF2::usage = "Use this function to mutate a Sum & Filter GRC model in any one of its component parameter settings";

RunRndWalk4GRCModelWithBiophyGroundedCRIF::usage = "Use the function below to perform an evolutionary exploration of parameter space of the GRC model implementing the thermodynamic state ensemble
													modeling approach. This algorithm performs a constrained mutational walk across parameter space (for one single configuration/GRC genotype) in
													search for a viable path along which the fitness landscape of GRC models will be climbed up stochastically (with neutral mutations being allowed
													over the course of evolution). The walk is performed via most frequently single changes in the parameter setting of the GRC, with double and triple
													mutations per step taken being allowed with a small probability. In this algorithm we implement the Metropolis rule to sporadically accept solutions
													where DeltaF < 0, with a given prob which is proportional to Exp[DeltaF/0.0001]";

(**************************************************************************************************************)

Begin["`Private`"]

(**************************************************************************************************************)

(*Load Viridis ColorMap*)
ClearAll[MPLColorMap] << "http://pastebin.com/raw/pFsb4ZBS";

ColorF01[Node_] := Switch[Node,"A",Darker[Gray,0.2],"B",Darker[Gray,0.2],"C",Darker[Gray,0.2],"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"A\",\nFontSize->0]\)\)]\)",Black,"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"B\",\nFontSize->0]\)\)]\)",Black,"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"C\",\nFontSize->0]\)\)]\)",Black];
EdgeColoringF[RegEffect_]:=Switch[RegEffect,1,Blue,-1,Red]


DisplayNonWeightedHaploidGraphWithMorphInputs[mKVsMatrix_List]:=Block[{},
HapRegMotifs=Prepend[Partition[(#[[1]]-> #[[2]])&/@Tuples[{"A","B","C"},2],3],{"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"A\",\nFontSize->0]\)\)]\)"->"A","\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"B\",\nFontSize->0]\)\)]\)"->"B","\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"C\",\nFontSize->0]\)\)]\)"->"C"}];

NodeA="A"-> {-1,1};
NodeB="B"-> {0,-0.25};
NodeC="C"-> {1,1};
MonA = "\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"A\",\nFontSize->0]\)\)]\)"-> {-1,1.65};
MonB = "\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"B\",\nFontSize->0]\)\)]\)"-> {-0.65,-0.25};
MonC = "\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"C\",\nFontSize->0]\)\)]\)"-> {1,1.65};

QuantitativeMatrix=Transpose[mKVsMatrix];
MaskMatrix=(If[#!=0,1,0]&/@#)&/@QuantitativeMatrix;
GraphRules=Flatten@Pick[HapRegMotifs,MaskMatrix,1];
Tags=Flatten@Pick[QuantitativeMatrix,MaskMatrix,1];
RegEffectList=Sign[Tags];
{EdgeComponents,Edge2ColorMap}={({#[[1]],#[[2]]}&/@GraphRules),(EdgeColoringF/@RegEffectList)};

CombinedGraphRules=Thread[{Flatten[GraphRules],Abs[Sign[Tags]]*0.0095}];

GraphPlot[CombinedGraphRules,VertexLabeling->True,DirectedEdges->{True,"ArrowheadsSize"->0.0001},
SelfLoopStyle->0.25,VertexCoordinateRules->{MonA,NodeA,MonB,NodeB,MonC,NodeC},
EdgeRenderingFunction->({Edge2ColorMap[[Flatten[Position[EdgeComponents,#2]][[1]]]],Thickness[#3],Arrow[#1,0.15]}&),
Method->{"Automatic","Rotation"->2 Pi},VertexRenderingFunction->({EdgeForm[ColorF01[#2]],ColorF01[#2],Disk[#1,0.15],Darker[Gray,0.009],Text[Style[#2,13,Bold,White],#1]}&),ImageSize->300]
]

(**************************************************************************************************************)

DisplayWeightedHaploidGraphWithMorphInputs[mKVsMatrix_List]:=Block[{},
HapRegMotifs=Prepend[Partition[(#[[1]]-> #[[2]])&/@Tuples[{"A","B","C"},2],3],{"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"A\",\nFontSize->0]\)\)]\)"->"A","\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"B\",\nFontSize->0]\)\)]\)"->"B","\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"C\",\nFontSize->0]\)\)]\)"->"C"}];

NodeA="A"-> {-1,1};
NodeB="B"-> {0,-0.25};
NodeC="C"-> {1,1};
MonA = "\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"A\",\nFontSize->0]\)\)]\)"-> {-1,1.65};
MonB = "\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"B\",\nFontSize->0]\)\)]\)"-> {-0.65,-0.25};
MonC = "\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"C\",\nFontSize->0]\)\)]\)"-> {1,1.65};

QuantitativeMatrix=Transpose[mKVsMatrix];
MaskMatrix=(If[#!=0,1,0]&/@#)&/@QuantitativeMatrix;
GraphRules=Flatten@Pick[HapRegMotifs,MaskMatrix,1];
Tags=Flatten@Pick[QuantitativeMatrix,MaskMatrix,1];
RegEffectList=Sign[Tags];
{EdgeComponents,Edge2ColorMap}={({#[[1]],#[[2]]}&/@GraphRules),(EdgeColoringF/@RegEffectList)};
MaxThickness=0.0095;
minKV=Min[Abs[Tags]];maxKV=Max[Abs[Tags]];
NormalizedKVs=Abs[((#-minKV)/(maxKV-minKV))]&/@Tags;
ScaledKVs=If[#<=0.01,(MaxThickness*0.95),MaxThickness*#]&/@NormalizedKVs;
CombinedGraphRules=Thread[{Flatten[GraphRules],ScaledKVs}];

GraphPlot[CombinedGraphRules,VertexLabeling->True,DirectedEdges->{True,"ArrowheadsSize"->0.0001},
SelfLoopStyle->0.25,VertexCoordinateRules->{MonA,NodeA,MonB,NodeB,MonC,NodeC},
EdgeRenderingFunction->({Edge2ColorMap[[Flatten[Position[EdgeComponents,#2]][[1]]]],Thickness[#3],Arrow[#1,0.15]}&),
Method->{"Automatic","Rotation"->2 Pi},VertexRenderingFunction->({EdgeForm[ColorF01[#2]],ColorF01[#2],Disk[#1,0.15],Darker[Gray,0.009],Text[Style[#2,13,Bold,White],#1]}&),ImageSize->300]
]



(**************************************************************************************************************)

SpaceTimeFeats4StripeFormingGRCs4SSMorpGradientWithBiophyGroundedCRIF0[LeakTerm_,HillCoeff_,{Params1_,Params2_,Params3_,Params4_}]:=Module[
{Prot,ProtStateVariab,ProtInitExpStat,Sol,DynSyst,GRCPhenotReadout},

GRNSize = Length[Params1];
(*RNApConc = 5.5; KRNAp=7.5;*)
NumNuclei=30;
IntTime = 500.;

(*Initial conditions for all variables, including boundary conditions set = 0*)
ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

ProtInitExpStat=Flatten[{Table[Table[Prot[a,n,0]==0.01,{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

DynSyst=Flatten[
{Table[
Table[D[Prot[a,n,t],t]==
( ( ComputeSigmoidCRIFWithLeakTerm[LeakTerm,HillCoeff,{Params2[[n]],Prot[1,n,t],Prot[2,n,t],Prot[3,n,t]},Params1[[a]]] ) + Params3[[a]]( Prot[a,n-1,t]+Prot[a,n+1,t]-2*Prot[a,n,t] ) - ( Params4[[a]]*Prot[a,n,t] )),
{a,1,GRNSize}],
{n,NumNuclei}],
Table[Table[D[Prot[a,n,t],t]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]}
];

Sol=NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

GRCPhenotReadout= Transpose[Thread[Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,Range[9,270,9]}],{i,3}]]]

]


SpaceTimeFeats4StripeFormingGRCs4SSMorpGradientWithSumAndFilterModel0[{Wmatrix_List,MorphInput_List,DiffParams_List,DegParams_List}]:=Module[
{Prot,ProtStateVariab,ProtInitExpStat,DynSyst,Sol,GRCPhenotReadout},

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*Transfer function for mapping integrated regulatory inputs to transcriptional outputs. In this model the SigmoidSteepness = alpha = 5, while the parameter b, which gives the location of the threshold value, is set to 1 *)
SumAndFilterF[SigmoidSteepness_,Threshold_,IntegratedRegInput_]:=1./(1.+Exp[(SigmoidSteepness- (Threshold*IntegratedRegInput))]);
(*A function for summing the weigthed regulatory inputs*)
AdditiveRegContributionF[a_,n_,t_,W_List,Morph_]:=Total[W[[a]]*Prepend[Table[Prot[b,n,t],{b,1,GRNSize}],Morph]];
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)

GRNSize = Length[Wmatrix];
(*As set in the paper referenced above*)
NumNuclei=30;
(*As set in the paper referenced above*)
alpha=5.;
sigmoidThreshold=1;
IntTime = 500.;

(*Initial conditions for all variables, including boundary conditions set = 0*)
ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];
ProtInitExpStat=Flatten[{Table[Table[Prot[a,n,0]==0.1,{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

DynSyst=Flatten[
{Table[
Table[D[Prot[a,n,t],t]==
( ( SumAndFilterF[alpha,sigmoidThreshold,AdditiveRegContributionF[a,n,t,Wmatrix,MorphInput[[n]]]] )+ ( DiffParams[[a]]( Prot[a,n-1,t]+Prot[a,n+1,t]-2*Prot[a,n,t] ) ) - ( DegParams[[a]]*Prot[a,n,t] )),{a,1,GRNSize}],
{n,NumNuclei}],Table[Table[D[Prot[a,n,t],t]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]}
];

Sol=NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

GRCPhenotReadout= Transpose[Thread[Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,Range[9,270,9]}],{i,3}]]]

]

(**************************************************************************************************************)

colorbar[{min_,max_},divs_: 150]:=DensityPlot[y,{x,0,0.1},{y,min,max},AspectRatio->12,PlotRangePadding->0,PlotPoints->{2,divs},MaxRecursion->0,FrameTicks->{None,Automatic,None,None},ColorFunction->MPLColorMap["Viridis"]]

CreatSpaceTimePlots4GRC[TimeSeries_List]:=GraphicsRow[MapIndexed[With[{opts={ImageSize->{Automatic,130},ImagePadding->5}},ArrayPlot[#,ColorFunction->MPLColorMap["Viridis"],FrameTicks->None,opts]]&,((#/Max[Flatten[#]])&/@TimeSeries)],Spacings->{-10,0}]

(**************************************************************************************************************)
(**************************************************************************************************************)
(**************************************************************************************************************)


AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[{Wmatrix_List,MorphInput_List,DiffParams_List,DegParams_List}]:=Block[
{GRCPhenotReadout},

GRCPhenotReadout = StripeFormingGRCs4SSMorpGradientWithSumAndFilterGeneRegModel[{Wmatrix,MorphInput,DiffParams,DegParams}];

(*This condition here is to check if there is any negative gene expression value, which is not biologically realistic. If there is, then Fitness must be 0. *)

If[
	Count[Flatten[GRCPhenotReadout],x_/;x<=0]>0,0,
	FitnessF4SingleStripe[GRCPhenotReadout]
	]
		]

(**************************************************************************************************************)
(**************************************************************************************************************)
AssessExpPattern4SingleStripe[ExpProfiles_List]:=Block[{},

(*Optimal expression pattern represented in a scale between 1-10 in expression level, being 10 an expression level which is > 90% of the maximal level observed along the 1D field of cells for the output genes*)
OptimalPattern=Flatten[{ConstantArray[1,10],ConstantArray[10,10],ConstantArray[1,10]}];
(*Maximal discrepancy achievable for a given pattern w.r.t to the optimal pattern above*)
Dmax=9*30;

ThresholdedExpValues={{0.`,0.1`},{0.1`,0.2`},{0.2`,0.3`},{0.3`,0.4`},{0.4`,0.5`},{0.5`,0.6`},{0.6`,0.7`},{0.7`,0.8`},{0.8`,0.9`},{0.9`,1.001`}};

testInt[Rng_,Value_]:= Rng[[1]]<=Value< Rng[[2]];

NormEquilbPattern=Last[ExpProfiles]/Max[Last[ExpProfiles]];

DiscretizedPattern=Table[Position[testInt[#,NormEquilbPattern[[x]]]&/@ThresholdedExpValues,True][[1,1]],{x,Length[NormEquilbPattern]}];

PFeff = 1-(ManhattanDistance[OptimalPattern,DiscretizedPattern]/Dmax)
]

(**************************************************************************************************************)
Q1[x_]:=(x^5)/( (x^5) + (0.1^5) ) (*(x^10)/( (x^10) + (0.15^10) );*)
Q2[x_]:=(0.1^2)/( (x^2) + (0.1^2) )
(**************************************************************************************************************)
PatternFilter[ExpProfiles_List]:=(1/Length[Last[ExpProfiles]])(EuclideanDistance[Last[ExpProfiles],ConstantArray[Mean[Last[ExpProfiles]],Length[Last[ExpProfiles]]]])

EquilibFilter[ExpProfiles_List]:=(1/Length[Last[ExpProfiles]])(EuclideanDistance[Last[ExpProfiles],First[ExpProfiles]])


FitnessF4SingleStripe[GRCPhenotReadout_]:=Block[{},

(*Testing for stability of the expression pattern after 500 time steps*)
EF=((1/3)*Total[EquilibFilter/@GRCPhenotReadout]);
(*Testing for spatial heterogeneity in the steady state expression pattern, which is also an indicative of whether the expression profiles are sufficiently high so that the cross regulatory interactions can be effective in controlling the expression putput. Note that unlike the FitnessF implemented before, in this FitnessF2 we only checked for the heterogeineity in the output node, and not in the whole set of nodes considered*)
PF=PatternFilter/@GRCPhenotReadout;

FS=Max[(AssessExpPattern4SingleStripe[#]*Q1[Min[PF]]*Q2[EF])&/@GRCPhenotReadout]
]

(**************************************************************************************************************)
(**************************************************************************************************************)

(*Modeling a steady state (exponentially decaying from the anterior-to-posterior axis) morphogen gradient. This morphogen distribution profile is similar to that exhibited by Bicoid in the Drosophila developing embryo*)
MGradientF[A0_,CellIndex_,DecayRate_]:=A0*Exp[-CellIndex/DecayRate]
NumNuclei=30;
SSInputMorphogen[A0_]:=MGradientF[A0,#,0.4]&/@Drop[Range[0,1,1/NumNuclei],-1];

GenerateComposedWiringMatrixTemplate:=Block[{},
RndWiring=Partition[RandomChoice[{-1,0,1},12],{4}];
RndWiring[[All,1]]={0,0,0};
RndWiring[[1,1]]=1;
RndWiring
]


GetProteinHalfLife[S_]:=RandomReal[{5,50},S];
GetProteinDegrad[S_]:=Log[2]/GetProteinHalfLife[S];
GetDiffRates[S_]:=RandomReal[{0,0.1},S]


GetRndComposedWiringMatrix:=Block[{},
RndW=GenerateComposedWiringMatrixTemplate*Partition[RandomReal[{0,1000},12],4];
RndW[[All,1]]={0,0,0};
RndW[[1,1]]=RandomReal[{0,1000},1][[1]];
RndW
]

GetRndComposedWiringMatrix2:=Block[{},
RndW=GenerateComposedWiringMatrixTemplate*Partition[RandomReal[{-10,10},12],4];
RndW[[All,1]]={0,0,0};
RndW[[1,1]]=RandomReal[{-10,10},1][[1]];
RndW
]

(**************************************************************************************************************)
(**************************************************************************************************************)

ComputeSigmoidCRIFWithLeakTerm[LeakTerm_,HillCoeff_,TFConc_List,AffFactors_List]:=Block[
{RegEffecList,Summands,RNApSummand,IndexPosReg,NumeratorPolynBase,NonContributingFactors,NumeratorPolyn,PartitionFunc,FracOccupancy},

RegEffecList=Sign[AffFactors];

Summands=((#)^HillCoeff)&/@(TFConc*Abs[AffFactors]);

IndexPosReg=Flatten[Position[RegEffecList,1]];

NumeratorPolynBase =Apply[Times,1+Summands[[IndexPosReg]]];

NonContributingFactors=Apply[Times,(1+Summands[[IndexPosReg]])];

NumeratorPolyn=NumeratorPolynBase-1+LeakTerm;

PartitionFunc=Apply[Times,1+Summands];

FracOccupancy=(NumeratorPolyn/PartitionFunc)

]

(**************************************************************************************************************)
(**************************************************************************************************************)

StripeFormingGRCs4SSMorpGradientWithSumAndFilterGeneRegModel[{Wmatrix_List,MorphInput_List,DiffParams_List,DegParams_List}]:=Module[
{Prot,ProtStateVariab,ProtInitExpStat,DynSyst,Sol,GRCPhenotReadout},

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*Transfer function for mapping integrated regulatory inputs to transcriptional outputs. In this model the SigmoidSteepness = alpha = 5, while the parameter b, which gives the location of the threshold value, is set to 1 *)
SumAndFilterF[SigmoidSteepness_,Threshold_,IntegratedRegInput_]:=1./(1.+Exp[(SigmoidSteepness- (Threshold*IntegratedRegInput))]);
(*A function for summing the weigthed regulatory inputs*)
AdditiveRegContributionF[a_,n_,t_,W_List,Morph_]:=Total[W[[a]]*Prepend[Table[Prot[b,n,t],{b,1,GRNSize}],Morph]];
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)

GRNSize = Length[Wmatrix];
(*As set in the paper referenced above*)
NumNuclei=30;
(*As set in the paper referenced above*)
alpha=5.;
sigmoidThreshold=1;
IntTime = 500.;

(*Initial conditions for all variables, including boundary conditions set = 0*)
ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];
ProtInitExpStat=Flatten[{Table[Table[Prot[a,n,0]==0.1,{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

DynSyst=Flatten[
{Table[
Table[D[Prot[a,n,t],t]==
( ( SumAndFilterF[alpha,sigmoidThreshold,AdditiveRegContributionF[a,n,t,Wmatrix,MorphInput[[n]]]] )+ ( DiffParams[[a]]( Prot[a,n-1,t]+Prot[a,n+1,t]-2*Prot[a,n,t] ) ) - ( DegParams[[a]]*Prot[a,n,t] )),{a,1,GRNSize}],
{n,NumNuclei}],Table[Table[D[Prot[a,n,t],t]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]}
];

Sol=NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

GRCPhenotReadout= Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,IntTime-250,IntTime,250}],{i,3}]
]

(**************************************************************************************************************)
(**************************************************************************************************************)

StripeFormingGRCs4SSMorpGradientWithBiophyGroundedCRIF0[LeakTerm_,HillCoeff_,{Params1_,Params2_,Params3_,Params4_}]:=Module[
{Prot,ProtStateVariab,ProtInitExpStat,Sol,DynSyst,GRCPhenotReadout},

GRNSize = Length[Params1];
(*RNApConc = 5.5; KRNAp=7.5;*)
NumNuclei=30;
IntTime = 500.;

(*Initial conditions for all variables, including boundary conditions set = 0*)
ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

ProtInitExpStat=Flatten[{Table[Table[Prot[a,n,0]==0.01,{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

DynSyst=Flatten[
{Table[
Table[D[Prot[a,n,t],t]==
( ( ComputeSigmoidCRIFWithLeakTerm[LeakTerm,HillCoeff,{Params2[[n]],Prot[1,n,t],Prot[2,n,t],Prot[3,n,t]},Params1[[a]]] ) + Params3[[a]]( Prot[a,n-1,t]+Prot[a,n+1,t]-2*Prot[a,n,t] ) - ( Params4[[a]]*Prot[a,n,t] )),
{a,1,GRNSize}],
{n,NumNuclei}],
Table[Table[D[Prot[a,n,t],t]==0.00,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]}
];

Sol=NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

GRCPhenotReadout= Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,IntTime-250,IntTime,250}],{i,3}]

]

(**************************************************************************************************************)
(**************************************************************************************************************)

Mutator4GRCW[GRCWG_]:=Block[{FlattenGRCWG,Positions2Mutate,MutID},
FlattenGRCWG=Flatten[GRCWG];
Positions2Mutate=Range[12];
MutID=RandomSample[Positions2Mutate,1][[1]];
MutStrength=RandomSample[{0.2,0.8}-> {0,RandomReal[{-1000,1000}]},1][[1]];
Partition[ReplacePart[FlattenGRCWG,MutID-> MutStrength],4]
]

Mutator4GRCW2[GRCWG_]:=Block[{FlattenGRCWG,Positions2Mutate,MutID},
FlattenGRCWG=Flatten[GRCWG];
Positions2Mutate=Range[12];
MutID=RandomSample[Positions2Mutate,1][[1]];
MutStrength=RandomSample[{0.2,0.8}-> {0,RandomReal[{-10,10}]},1][[1]];
Partition[ReplacePart[FlattenGRCWG,MutID-> MutStrength],4]
]

(**************************************************************************************************************)
(**************************************************************************************************************)

MutateParamVector[InputParams_List,MID_]:=Block[{InputParamsVect,MutID},
NewParam=Switch[MID,2,GetDiffRates[1][[1]],3,GetProteinDegrad[1][[1]]];
InputParamsVect=InputParams;
MutID=RandomSample[{1,2,3},1][[1]];
InputParamsVect[[MutID]]=NewParam;
InputParamsVect
]

(**************************************************************************************************************)
(**************************************************************************************************************)


(**************************************************************************************************************)
(**************************************************************************************************************)

(*Use this mutator F for the thermodynamically-grounded CRIF model*)
MutatorF[MutID_,InputParamSet_List]:=Block[{},

ParamSet=InputParamSet;

MutatedParam=If[MutID==1,Mutator4GRCW[ParamSet[[MutID]]],
			      MutateParamVector[ParamSet[[MutID]],MutID]
			       ];

ParamSet[[MutID]]=MutatedParam;

ParamSet
]

(*Use this mutator F for the Sum & Filter model*)
MutatorF2[MutID_,InputParamSet_List]:=Block[{},

ParamSet=InputParamSet;

MutatedParam=If[MutID==1,Mutator4GRCW2[ParamSet[[MutID]]],
			      MutateParamVector[ParamSet[[MutID]],MutID]
			       ];

ParamSet[[MutID]]=MutatedParam;

ParamSet
]

(**************************************************************************************************************)
(**************************************************************************************************************)

RunRndWalk4SumAndFilterGRCModel[Repl_]:=Block[{EvolSteps,thresholdFitness,ParamSettingSamplingRate,
MutaProbList,GRCWiring,MorphogenInput,DiffRates,ProtDegRates,GetInputParamSet,GRCDyns,RefInputParamSetting,Fb,
Fa,DeltaF,ThisMutId,InputParamSetting,InputParamSettingAllin,OutputFileName,PatternFormingSeeds,DestinationFolder,
ParamsSamplingCond,SamplingPointsList,BoolWMatrix,WContainer},

EvolSteps=5000000;
ParamSettingSamplingRate=500;
thresholdFitness=0.95;
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
DestinationFolder = "GRCModelWithSumAndFilterCRIF/SampledGRCTopologies";
(*To ensure that the destination folders exist, if don't then create them!*)

If[
	FileExistsQ[DestinationFolder]==False,
	CreateDirectory[DestinationFolder]
	];

(*To read in within mathematica these files do: ReadList[FitnessFeaturesFile]*)
OutputFileName = DestinationFolder<>"/ParamGenotypeGRC_Repl_"<>ToString[Repl]<>".dat";
(***********************************************************************)
(***********************************************************************)

(*Start exploration of the design space of GRCs using previously, highly-optimized (F>=0.95), parameters/genotypes*)
PatternFormingSeeds=DeleteDuplicates[ReadList["GRCModelWithSumAndFilterCRIF/PatternFormingSeeds4SFGRModelExploration.dat"]];

{GRCWiring,MorphogenInput,DiffRates,ProtDegRates} = PatternFormingSeeds[[Repl]];(*{GetRndComposedWiringMatrix2,SSInputMorphogen[1],GetDiffRates[3],GetProteinDegrad[3]};*)

BoolWMatrix=Sign/@GRCWiring;
AppendTo[WContainer,BoolWMatrix];

(*GRCDyns = StripeFormingGRCs4SSMorpGradientSumAndFilterModel0[{GRCWiring,MorphogenInput,DiffRates,ProtDegRates}];*)

Fb=AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[{GRCWiring,MorphogenInput,DiffRates,ProtDegRates}];

RefInputParamSetting={GRCWiring,DiffRates,ProtDegRates};

(*Optimization cycle to explore design space of GRCs*)
Do[

Clear[InputParamSetting,InputParamSettingAllin,GRCDyns,Fa,DeltaF];

ThisMutId = MutId;

InputParamSetting=Nest[MutatorF2[ThisMutId,#]&,RefInputParamSetting,NumMutants];

InputParamSettingAllin=Insert[InputParamSetting,MorphogenInput,2];

Fa=AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[InputParamSettingAllin];

DeltaF = Fa-Fb;

(*If[DeltaF> 0.001,Print[{Fa,iter}]];*)

(*Implementing a Metropolis algorithm*)
If[
	DeltaF>=0,
	RefInputParamSetting=InputParamSetting;
	Fb=Fa,
	If[DeltaF<0,
      If[Exp[DeltaF/0.01]>= Random[],
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


(**************************************************************************************************************)
(**************************************************************************************************************)

RunRndWalk4SumAndFilterGRCModel2[Repl_, EvolSteps_:5000000, ParamSettingSamplingRate_:500, thresholdFitness_:0.95]:=Block[{MutaProbList,
	GRCWiring,MorphogenInput,DiffRates,ProtDegRates,GetInputParamSet,GRCDyns,RefInputParamSetting,Fb,Fa,DeltaF,ThisMutId,InputParamSetting,
	InputParamSettingAllin,OutputFileName,PatternFormingSeeds,DestinationFolder,ParamsSamplingCond,SamplingPointsList,BoolWMatrix,WContainer},

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
DestinationFolder = "GRCModelWithSumAndFilterCRIF/SampledGRCTopologies1";
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
{GRCWiring,MorphogenInput,DiffRates,ProtDegRates} = {GetRndComposedWiringMatrix2,SSInputMorphogen[1],GetDiffRates[3],GetProteinDegrad[3]};

(*GRCDyns = StripeFormingGRCs4SSMorpGradientSumAndFilterModel0[{GRCWiring,MorphogenInput,DiffRates,ProtDegRates}];*)

Fb=AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[{GRCWiring,MorphogenInput,DiffRates,ProtDegRates}];

RefInputParamSetting={GRCWiring,DiffRates,ProtDegRates};

(*Optimization cycle to explore design space of GRCs*)
Do[

Clear[InputParamSetting,InputParamSettingAllin,GRCDyns,Fa,DeltaF];

ThisMutId = MutId;

InputParamSetting=Nest[MutatorF2[ThisMutId,#]&,RefInputParamSetting,NumMutants];

InputParamSettingAllin=Insert[InputParamSetting,MorphogenInput,2];

Fa=AssessFitnessScoreStripeFormingGRCs4SSMorpGradient[InputParamSettingAllin];

DeltaF = Fa-Fb;

If[DeltaF> 0.001,Print[{Fa,iter}]];

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


(**************************************************************************************************************)
(**************************************************************************************************************)

RunRndWalk4GRCModelWithBiophyGroundedCRIF[Repl_,LeakTerm_,HillCoeff_]:=Block[{GetInputParamSet,GRCDyns,RefInputParamSetting},

EvolSteps=600000;
thresholdFitness=0.925;


(*Only mutate the wiring (with the highest probab), as well as diffusion and protein degradation rates*)
MutId:=RandomSample[{0.7,0.15,0.15}-> {1,2,3},1][[1]];

(*Select number of mutations to perform between 1-3*)
NumMutants:=RandomSample[{0.9,0.1}-> {1,2},1][[1]];

(***********************************************************************)
(***********************************************************************)
DestinationFolder="GRCModelWithBiophyGroundedCRIF/SimTests/LeakTerm_"<>ToString[LeakTerm]<>"_HillCoeff_"<>ToString[HillCoeff];
(*To ensure that the destination folders exist, if don't then create them!*)

If[
FileExistsQ[DestinationFolder]==False,
CreateDirectory[DestinationFolder]
];

(*To read in within mathematica these files do: ReadList[FitnessFeaturesFile]*)
OutputFileName = DestinationFolder<>"/OutputParams4Repl_"<>Repl<>".dat";
(***********************************************************************)
(***********************************************************************)

{GRCWiring,MorphogenInput,DiffRates,ProtDegRates}={GetRndComposedWiringMatrix,SSInputMorphogen[1],GetDiffRates[3],GetProteinDegrad[3]};

GRCDyns=StripeFormingGRCs4SSMorpGradientWithBiophyGroundedCRIF0[LeakTerm,HillCoeff,{GRCWiring,MorphogenInput,DiffRates,ProtDegRates}];

Fb=If[
		Count[Flatten[GRCDyns],x_/;x<=0]>0,0,
	    FitnessF3[GRCDyns]
	];

RefInputParamSetting={GRCWiring,DiffRates,ProtDegRates};

Do[

Clear[InputParamSetting,InputParamSettingAllin,GRCDyns,Fa,DeltaF];

ThisMutId = MutId;

InputParamSetting=Nest[MutatorF[ThisMutId,#]&,RefInputParamSetting,NumMutants];

InputParamSettingAllin=Insert[InputParamSetting,MorphogenInput,2];

GRCDyns=StripeFormingGRCs4SSMorpGradientWithBiophyGroundedCRIF0[LeakTerm,HillCoeff,InputParamSettingAllin];

Fa=If[Count[Flatten[GRCDyns],x_/;x<=0]>0,0,
	  FitnessF3[GRCDyns]
	];

DeltaF = Fa-Fb;

If[Fb>=thresholdFitness && DeltaF>= 0.001,PutAppend[InputParamSettingAllin,OutputFileName]];

If[
DeltaF>=0,
RefInputParamSetting=InputParamSetting;
Fb=Fa,
If[DeltaF<0,
     If[Exp[DeltaF/0.0005]>= Random[],
	 RefInputParamSetting=InputParamSetting;
	 Fb=Fa,
	 RefInputParamSetting=RefInputParamSetting;
	 Fb=Fb]
	                     ]
				];
(*
If[iter==EvolSteps && DeltaF>= 0.005 && Fb>  thresholdFitness,
PutAppend[InputParamSettingAllin,OutputFileName]
]
*)

If[MemberQ[Drop[Range[0,EvolSteps,6000],{1}],iter] && Fb>= thresholdFitness,
PutAppend[InputParamSetting,OutputFileName]
]

,{iter,0,EvolSteps}]
]


(**************************************************************************************************************)
(**************************************************************************************************************)





End[]

EndPackage[]
(**************************************************************************************************************)
