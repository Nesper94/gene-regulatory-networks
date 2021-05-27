(* ::Package:: *)
(* Author: Juan Camilo Arboleda Rivera based on code by Jayson Gutiérrez*)

Needs["parameters`"]

(*This module depends on DesignMorphogeneResponsiveGRCs.m and
EvolAlgorithm4PatternFormingGRCModel2MathPackageV9.m*)
<< EvolAlgorithm4PatternFormingGRCModel2MathPackageV9`
<< DesignMorphogeneResponsiveGRCs`

(*The following is basically the function
AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM2*)
fitness::usage="fitness[ GRN_List, A0_:1, h_:calculateMorphDecay[NumNuclei, A0], ICs_:ConstantArray[1, 3*NumNuclei] ] calculates the fitness of the Gene Regulatory Network (GRN) using A0 as the max morphogen concentration, h as the morphogen decay parameter and ICs as the initial concentrations of the gene products in each one of the cells of the morphogenetic field."

fitness[GRN_List, A0_:1, h_:calculateMorphDecay[NumNuclei, A0], ICs_:ConstantArray[1, 3*NumNuclei] ] :=
 Block[{NullMorpInput, MorpInputProfile, PreMorpInputFS,
   SSExpValuesPreMorpInput, FS4GRCInResponse2MorpInput, NumNuclei},
  (*Set ICs for all genes in all nuclei and run system without Morphogene input*)
  NullMorpInput = ConstantArray[0, NumNuclei];
  MorpInputProfile = SSInputMorphogen[A0, h];

  {PreMorpInputFS, SSExpValuesPreMorpInput} =
   AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,
    NullMorpInput, GRN];

  (*Set the pre-Morphogen SS expression levels as ICs for simulating the
  dynamics of the system in the presence of the Morphogene*)(*Check that
  fitness before applying the Morphogene input is>0 (reaches stationary
  behavior) and is\[LessEqual]0.5 (is unable to generate any striped pattern at
  all), if this condition is met then assess the ability of the system to
  interpret the input and generate a striped pattern for any of the genes*)
  FS4GRCInResponse2MorpInput =
   If[0 < PreMorpInputFS <= 0.5,
    AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM2[
      SSExpValuesPreMorpInput, MorpInputProfile,
      GRN][[1]], 0]
]


modifiedAssessFitnessScore4StripePattern4SSMorpGradient4SFGRM2[ICs_List,
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

        (*FS = Max[
          MapIndexed[(AssessExpPattern4SingleStripeFeat[#]*
              Q1[PF[[#2[[1]]]]]*SSCondition) &, EndPointExpPatterns]]]*)

          FS = MapIndexed[(AssessExpPattern4SingleStripeFeat[#]*
                  Q1[PF[[#2[[1]]]]]*SSCondition) &, EndPointExpPatterns]]
