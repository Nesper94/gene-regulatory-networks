(* ::Package:: *)
BeginPackage["VideosExpressionDynamics`"]

(*Package to analyze the expression dynamics of each genotype.
This package depends on EvolAlgorithm4PatternFormingGRCModel2MathPackageV9*)

expEspacioTemporal::usage = "Create expression profile.";
videoExp::usage = "videoExp[ genotipo ] create a video of the expression \
dynamics of a GRN, it takes as argument 'gen', with options 'a','b' and 'c', \
image size and 'fps' (frames per second).";

Begin[ "`Private`"]
Needs["parameters`"]
Needs["EvolAlgorithm4PatternFormingGRCModel2MathPackageV9`"]
Needs["DesignMorphogeneResponsiveGRCs`"]

(*Import list of classified, ordered, non isomorphic motifs*)
motivosClasificadosOrdenadosNoIsomorfos = Flatten[ToExpression[#],1] &/@
Import["../data/ordered-non-isomorph-grns.txt","TSV"];

(*This function shows the positions of GRNs with some topology*)
posici\[OAcute]nMotivos[topology_Integer]:=Flatten[{First[#],Last[#]}&/@{Position[motivosClasificadosOrdenadosNoIsomorfos[[All,2]],topology]}]

(*Here the function "SpaceTimeFeats4StripeFormingGRCs4SSMorpGradientWithSumAndFilterModel0" is redefined so that it takes GRNs as input without need to supply the morphogen input*)
expEspacioTemporal[{Wmatrix_List, DiffParams_List, DegParams_List},ti_:0,tf_:300,step_:6] :=
  Module[
  {Prot, ProtStateVariab, ProtInitExpStat, DynSyst, Sol, GRCPhenotReadout},

  (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
  (*Transfer function for mapping integrated regulatory inputs to transcriptional
  outputs. In this model the SigmoidSteepness = alpha = 5, while the parameter b,
  which gives the location of the threshold value, is set to 1 *)
  ICs=ConstantArray[1, 3*NumNuclei];
  NullMorpInput=ConstantArray[0, NumNuclei];
  {PreMorpInputFS,SSExpValuesPreMorpInput} = AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,NullMorpInput,{Wmatrix, DiffParams, DegParams}];

  SumAndFilterF[SigmoidSteepness_, Threshold_, IntegratedRegInput_] := 1./(1. + Exp[(SigmoidSteepness - (Threshold*IntegratedRegInput))]);
  (*A function for summing the weigthed regulatory inputs*)

  AdditiveRegContributionF[a_,n_,t_,W_List,Morph_] := Total[W[[a]]*Prepend[Table[Prot[b, n, t], {b,1,GRNSize}], Morph]];
  (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
  MorphInput = SSInputMorphogen[1];  (*Line added by J.C.A.R*)

  GRNSize = Length[Wmatrix];

  (*Initial conditions for all variables, including boundary conditions set = 0*)

  ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

  ProtInitExpStat=Flatten[{Thread[Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,NumNuclei}]]==SSExpValuesPreMorpInput],Flatten[Table[Table[Prot[a,n,0]==0.0,{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

  DynSyst=Flatten[
    {Table[
      Table[D[Prot[a,n,t],t]==
        ( ( SumAndFilterF[alpha,sigmoidThreshold,AdditiveRegContributionF[a,n,t,Wmatrix,MorphInput[[n]]]] ) + (
           DiffParams[[a]] ( Prot[a,n-1,t] + Prot[a,n+1,t] - 2*Prot[a,n,t] ) ) - (
           DegParams[[a]]*Prot[a, n, t] )), {a, 1, GRNSize}],
      {n, NumNuclei}],
     Table[Table[ D[Prot[a, n, t], t] == 0.00, {a, 1, GRNSize}], {n, {0, NumNuclei + 1}}]}
    ];

  Sol = NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,
  {t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

  GRCPhenotReadout= Transpose[Thread[Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,Range[ti,tf,step]}], {i, 3}]]]
  ]

(******************************************************************************)

(*Function that creates videos of expression dynamics*)
videoExp[ genotipo_List, OptionsPattern[]] := Module[{patron},
  Options[videoExp] = {gen -> 0, tama\[NTilde]o -> UpTo[300], PlotRange -> {0,25}, fps -> Automatic};
    patron = expEspacioTemporal[ genotipo[[1]] ];
  If[ SubsetQ[ {a, b, c}, Flatten[{ OptionValue[gen] }] ],
   ListAnimate[
    Table[ ListPlot[ Table[
       Prepend[{patron[[ # /. {a -> 1, b -> 2, c -> 3} ]] [[n]] [[i]]}, i]
       , {i, NumNuclei}], Joined -> True, PlotRange -> OptionValue[PlotRange],
      ImageSize -> OptionValue[tama\[NTilde]o] ], {n, NumNuclei}], OptionValue[fps] ] &/@
      Flatten[{OptionValue[gen]} ],
   Table[
    ListAnimate[
     Table[ ListPlot[ Table[ Prepend[ {patron[[k]][[n]][[i]]}, i]
        , {i, NumNuclei}], Joined -> True, PlotRange -> OptionValue[PlotRange],
        ImageSize -> OptionValue[tama\[NTilde]o] ], {n, NumNuclei}], OptionValue[fps] ]
    , {k, 3} ] ] ]

(********************************************************************************************************************)

(*Redefine function "AssessFitnessScoreStripeFormingGRCs4SSMorpGradient"
to use it without specifying the morphogen input*)
AssessFitnessScoreStripeFormingGRCs4SSMorpGradientJCAR[{Wmatrix_List,DiffParams_List,DegParams_List}]:=Block[
{GRCPhenotReadout},
MorphInput = SSInputMorphogen[1];  (*Line added by J.C.A.R*)
GRCPhenotReadout = StripeFormingGRCs4SSMorpGradientWithSumAndFilterGeneRegModel[{Wmatrix,MorphInput,DiffParams,DegParams}];

(*This condition here is to check if there is any negative gene expression
value, which is not biologically realistic. If there is, then Fitness must be 0. *)

If[
	Count[Flatten[GRCPhenotReadout],x_/;x<=0]>0,0,
	FitnessF4SingleStripe[GRCPhenotReadout]
	]
		]

End[]
EndPackage[]
