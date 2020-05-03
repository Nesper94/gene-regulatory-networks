(* ::Package:: *)

(*USE THE FOLLOWING METHODS TO PERFORM SCALAR AND PAIRWISE GENE-GENE INTERACTION PERTURBATIONS FOR 3-NODE MORPHOGENE RESPONSIVE GRCs*)
(************************************************************************************************************************)
(************************************************************************************************************************)

(*With this function we generate a full description of which entries \
in the W interaction matrix are perturbed (scalar perturbation), \
either up or down the nominal value, and to what extent this value is \
perturbed*)

GeneratePerturbSetDescription[FocalParamSet_, PerturbStrength_] := 
 Block[{GeneGeneInteractionMatrixWSymbs, 
   GeneGeneInteractionMatrixWIDs, ExistingWIDs, 
   ExistingGeneGeneInteractionMatrixWSymbs, 
   ExistingGeneGeneInteractionMatrixWIDs, UpPerturbIDs, 
   DownPerturbIDs},
  
  GeneGeneInteractionMatrixWSymbs = 
   Flatten[Transpose[Array[Subscript[w, # -> #2] &, {4, 3}]]];
  
  GeneGeneInteractionMatrixWIDs = Flatten[Array[{#, #2} &, {3, 4}], 1];
  
  ExistingWIDs = 
   Flatten[Position[Flatten[First[FocalParamSet]], _?(# != 0 &)]];
  
  ExistingGeneGeneInteractionMatrixWSymbs = 
   GeneGeneInteractionMatrixWSymbs[[ExistingWIDs]];
  
  ExistingGeneGeneInteractionMatrixWIDs = 
   GeneGeneInteractionMatrixWIDs[[ExistingWIDs]];
  
  UpPerturbIDs = 
   Append[#, PerturbStrength] & /@ (Append[#, "up"] & /@ 
      Thread[{ExistingGeneGeneInteractionMatrixWIDs, 
        ExistingGeneGeneInteractionMatrixWSymbs}]);
  
  DownPerturbIDs = 
   Append[#, PerturbStrength] & /@ (Append[#, "down"] & /@ 
      Thread[{ExistingGeneGeneInteractionMatrixWIDs, 
        ExistingGeneGeneInteractionMatrixWSymbs}]);
  
  Flatten[Thread[{UpPerturbIDs, DownPerturbIDs}], 1]
  ]


(*Use this function to create sort of an allelic series from a given \
reference parameterized wiring (genotype) of a GRC. The set of \
perturbations involved increased and decreased (by a given predefined \
Perturbation Strength) single gene-gene interactions. Note: each \
existing gene-gene interaction in the reference wiring is perturbed \
up/down and the list of perturbed wirings returned*)

GenerateSingleGeneGeneIntPertSet[FocalParamSet_, PerturbStrength_] := 
 Block[{GeneGeneIntParams, PositNonNullGeneGeneInt, PertUp, PertDown, 
   SingleGeneGeneIntPertSet},
  
  GeneGeneIntParams = Flatten[First[FocalParamSet]];
  
  PositNonNullGeneGeneInt = 
   Flatten[Position[GeneGeneIntParams, x_ /; x != 0]];
  
  {PertUp, PertDown} = 
   Transpose[{(# (1 + PerturbStrength)), (# (1 - 
           PerturbStrength))} & /@ 
     GeneGeneIntParams[[PositNonNullGeneGeneInt]]];
  
  SingleGeneGeneIntPertSet = 
   Flatten[Thread[{Partition[#, 4] & /@ 
       MapIndexed[
        ReplacePart[GeneGeneIntParams, # -> PertUp[[#2[[1]]]] ] &, 
        PositNonNullGeneGeneInt], 
      Partition[#, 4] & /@ 
       MapIndexed[
        ReplacePart[GeneGeneIntParams, # -> PertDown[[#2[[1]]]] ] &, 
        PositNonNullGeneGeneInt]}], 1];
  
  SingleGeneGeneIntPertSet
  ]



(*Generating up and down scalar perturbations for each single wiring \
parameter, where nominal parameters are increased/decreased in 1 \
percent w.r.t. the nominal parameters. In the following we increase \
and decrease (in 1% every time) every single parameter from 1% to \
100%. This creates two lists, which are sort of mutational variants \
of a GRC*)
RunScalarPerturbationExps[FocalParamSet_, 
  PerturbationStrenghts_] := 
 Block[{NumOfGeneGeneInteracts, UpPerturbs, DownPerturbs, 
   FitnessEffects4UpPerturbs, FitnessEffects4DownPerturbs, RefFitness},
  
  NumOfGeneGeneInteracts = 
   Count[Flatten[First[FocalParamSet]], x_ /; x != 0];
  
  {UpPerturbs, DownPerturbs} = 
   Catenate /@ 
    Transpose[
     Table[Partition[
       ReplacePart[FocalParamSet, 1 -> #] & /@ 
        GenerateSingleGeneGeneIntPertSet[FocalParamSet, x], 
       NumOfGeneGeneInteracts], {x, PerturbationStrenghts}]];
  
  RefFitness = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[
    FocalParamSet];
  
  FitnessEffects4UpPerturbs = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM /@
     UpPerturbs;
  
  FitnessEffects4DownPerturbs = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM /@
     DownPerturbs;
  
  (*return relative fitness effects of both types (up/
  down) of perturbations*)
  {(FitnessEffects4UpPerturbs/
     RefFitness), (FitnessEffects4DownPerturbs/RefFitness)}
  ]



(*Run a set of scalar perturbations*)

RunScalarPerturbationExps[FocalParamSet_, PerturbationStrenghts_] := 
 Block[{NumOfGeneGeneInteracts, ScalarPerturbedWMatrices, 
   FitnessEffects4UpDownPerturbs, RefFitness},
  
  NumOfGeneGeneInteracts = 
   Count[Flatten[First[FocalParamSet]], x_ /; x != 0];
  
  RefFitness = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[
    FocalParamSet];
  
  ScalarPerturbedWMatrices = 
   ReplacePart[FocalParamSet, 1 -> #] & /@ 
    Catenate[
     GenerateSingleGeneGeneIntPertSet[FocalParamSet, #] & /@ 
      PerturbationStrenghts];
  
  FitnessEffects4UpDownPerturbs = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM /@
     ScalarPerturbedWMatrices;
  
  (*return relative fitness effects of both types (up/
  down) of perturbations*)
  {FitnessEffects4UpDownPerturbs, 
   Abs[FitnessEffects4UpDownPerturbs - RefFitness]/RefFitness}
  ]



(*Generate a database of perturbations upon simulating a large set of \
(scalar) single-point mutations on a target GRC genotype*)

GenerateScalarPertubationDataSet4Genotype[FocalParamSet_, 
  PerturbationStrenghts_] := 
 Block[{CreatePerturbSetDescription, dataset0, dataset1, 
   FitnessEffects4UpDownPerturbs, RelFitnessChanges4UpDownPerturbs},
  
  CreatePerturbSetDescription = 
   Prepend[Catenate[
     GeneratePerturbSetDescription[FocalParamSet, #] & /@ 
      PerturbationStrenghts], {"MatrixID", "WInteraction", 
     "Perturbation", "Magnitude"}];
  
  dataset0 = 
   With[{header = First@CreatePerturbSetDescription}, 
     AssociationThread[header -> #] & /@ 
      Rest@CreatePerturbSetDescription] // Dataset;
  
  {FitnessEffects4UpDownPerturbs, RelFitnessChanges4UpDownPerturbs} = 
   RunScalarPerturbationExps[FocalParamSet, PerturbationStrenghts];
  
  dataset1 = 
   dataset0[
    MapThread[
      Append[#1, "RelFitnessChange" -> #2] &, {#, 
       RelFitnessChanges4UpDownPerturbs}] &];
  
  dataset1
  ]


(*Use this function to generate a list of pairwise gene-gene \
interaction perturbations for a fixed up and down magnitude around \
nominal values*)

GeneratePairwiseGeneGeneIntPertSet[FocalParamSet_, PerturbStrength_] :=
  Block[{GeneGeneIntParams, PositNonNullGeneGeneInt, PertUp, PertDown,
    GeneGeneInteractionMatrixWSymbs, 
   ExistingGeneGeneInteractionMatrixWSymbs, 
   ExistingEntriesInGeneGeneInteractionMatrixW, Feats4Up, Feats4Down, 
   PairwisePerturbationIDsSet, RefPositions, WIDs, WArrayIDs, 
   PerturbedWs, PerturbTypes, PerturbMagnitude, 
   PairwisePerturbedGRCWiring},
  
  GeneGeneIntParams = Flatten[First[FocalParamSet]];
  
  PositNonNullGeneGeneInt = 
   Flatten[Position[GeneGeneIntParams, x_ /; x != 0]];
  
  {PertUp, PertDown} = 
   Transpose[{(# (1 + PerturbStrength)), (# (1 - 
           PerturbStrength))} & /@ 
     GeneGeneIntParams[[PositNonNullGeneGeneInt]]];
  
  GeneGeneInteractionMatrixWSymbs = 
   Flatten[Transpose[Array[Subscript[w, # -> #2] &, {4, 3}]]];
  
  ExistingGeneGeneInteractionMatrixWSymbs = 
   GeneGeneInteractionMatrixWSymbs[[PositNonNullGeneGeneInt]];
  
  ExistingEntriesInGeneGeneInteractionMatrixW = 
   First[Position[
       Partition[GeneGeneInteractionMatrixWSymbs, 4], #]] & /@ 
    ExistingGeneGeneInteractionMatrixWSymbs;
  
  Feats4Up = 
   Thread[{PositNonNullGeneGeneInt, 
     ExistingEntriesInGeneGeneInteractionMatrixW, 
     ExistingGeneGeneInteractionMatrixWSymbs, PertUp, 
     ConstantArray["up", Length[PositNonNullGeneGeneInt]]}];
  
  Feats4Down = 
   Thread[{PositNonNullGeneGeneInt, 
     ExistingEntriesInGeneGeneInteractionMatrixW, 
     ExistingGeneGeneInteractionMatrixWSymbs, PertDown, 
     ConstantArray["down", Length[PositNonNullGeneGeneInt]]}];
  
  PairwisePerturbationIDsSet = 
   DeleteCases[Subsets[Join[Feats4Up, Feats4Down], {2}], 
    x_ /; (First[x[[1]]] == First[x[[2]]])];
  
  {RefPositions, WIDs, WArrayIDs, PerturbedWs, PerturbTypes} = 
   Transpose[Thread /@ PairwisePerturbationIDsSet];
  PerturbMagnitude = 
   ConstantArray[PerturbStrength, 
    Length[PairwisePerturbationIDsSet]];
  PairwisePerturbedGRCWiring = 
   Partition[ReplacePart[GeneGeneIntParams, #], 4] & /@ 
    Table[Thread[RefPositions[[i]] -> PerturbedWs[[i]]], {i, 
      Length[PairwisePerturbationIDsSet]}];
  
  Thread[{WIDs, WArrayIDs, PerturbTypes, PerturbMagnitude, 
    PairwisePerturbedGRCWiring}]
  ]



(*Use this function to generate a list of pairwise gene-gene \
interaction perturbations for a fixed up and down magnitude around \
nominal values*)
GeneratePairwiseGeneGeneIntPertSet2[FocalParamSet_, MutationalStep_] :=
  Block[{GeneGeneIntParams, PositNonNullGeneGeneInt, PertUp, PertDown,
    GeneGeneInteractionMatrixWSymbs, 
   ExistingGeneGeneInteractionMatrixWSymbs, 
   ExistingEntriesInGeneGeneInteractionMatrixW, Feats4Up, Feats4Down, 
   PairwisePerturbationIDsSet, RefPositions, WIDs, WArrayIDs, 
   PerturbedWs, PerturbTypes, PerturbMagnitude, 
   PairwisePerturbedGRCWiring, AllPairwisePerturbStrengths, 
   AppliedAllPairwisePerturbStrengths, PerturbDirectionIDs, 
   AllPairwiseGeneIDs, AllPairwisePerturbedParameterizedWirings, 
   ExistingCombinationsEntriesInGeneGeneInteractionMatrixW, 
   ExistingCombinationsGeneGeneInteractionMatrixWSymbs},
  
  GeneGeneIntParams = Flatten[First[FocalParamSet]];
  
  PositNonNullGeneGeneInt = 
   Flatten[Position[GeneGeneIntParams, x_ /; x != 0]];
  
  (*Generate all possible pairwise perturbations ranging from -0.95 \
to 0.95 at intervals of 0.05*)
  
  AllPairwisePerturbStrengths = 
   Tuples[DeleteCases[Range[-0.95, 0.95, MutationalStep], 
     x_ /; (x > -MutationalStep) && (x < MutationalStep)], 2];
  
  (*And the resulting relative change is applied here, 
  which can be above or below the nominal values*)
  
  AppliedAllPairwisePerturbStrengths = (1 + 
     AllPairwisePerturbStrengths);
  
  PerturbDirectionIDs = 
   Partition[
    If[# < 1, "down", "up"] & /@ 
     Flatten[AppliedAllPairwisePerturbStrengths], 2];
  
  (*Generating all possible pairs of gene-
  gene interactions to be perturbed simultaneously*)
  
  AllPairwiseGeneIDs = Subsets[PositNonNullGeneGeneInt, {2}];
  
  (*Generating all conceivable pairwise perturbed wirings, 
  within a given perturbation range*)
  
  AllPairwisePerturbedParameterizedWirings = 
   Catenate[
    Table[Partition[
        ReplacePart[GeneGeneIntParams, 
         Thread[AllPairwiseGeneIDs[[x]] -> \
(GeneGeneIntParams[[AllPairwiseGeneIDs[[x]]]]*#)]], 4] & /@ 
      AppliedAllPairwisePerturbStrengths, {x, 
      Length[AllPairwiseGeneIDs]}]];
  
  GeneGeneInteractionMatrixWSymbs = 
   Flatten[Transpose[Array[Subscript[w, # -> #2] &, {4, 3}]]];
  
  ExistingGeneGeneInteractionMatrixWSymbs = 
   GeneGeneInteractionMatrixWSymbs[[PositNonNullGeneGeneInt]];
  
  ExistingCombinationsGeneGeneInteractionMatrixWSymbs = 
   Subsets[ExistingGeneGeneInteractionMatrixWSymbs, {2}];
  
  ExistingEntriesInGeneGeneInteractionMatrixW = 
   First[Position[
       Partition[GeneGeneInteractionMatrixWSymbs, 4], #]] & /@ 
    ExistingGeneGeneInteractionMatrixWSymbs;
  
  ExistingCombinationsEntriesInGeneGeneInteractionMatrixW = 
   Subsets[ExistingEntriesInGeneGeneInteractionMatrixW, {2}];
  
  {WIDs, WArrayIDs, PerturbTypes, PerturbMagnitude} = 
   Transpose[
    Catenate[
     Table[
      Thread[{ConstantArray[
         ExistingCombinationsEntriesInGeneGeneInteractionMatrixW[[x]],
          Length[AllPairwisePerturbStrengths]], 
        ConstantArray[
         ExistingCombinationsGeneGeneInteractionMatrixWSymbs[[x]], 
         Length[AllPairwisePerturbStrengths]], PerturbDirectionIDs, 
        Abs[AllPairwisePerturbStrengths]}], {x, 
       Length[AllPairwiseGeneIDs]}]]];
  
  Thread[{WIDs, WArrayIDs, PerturbTypes, PerturbMagnitude, 
    AllPairwisePerturbedParameterizedWirings}]
  ]


(*Use function to run a full pairwise perturbation experiment,at the \
end of which a well structured Dataset is returned*)
GeneratePairwisePertubationDataSet4Genotype[FocalParamSet_, 
  MutationalStep_] := 
 Block[{PairwiseGeneGeneIntPertSet, PairwiseGeneGeneIntPertDataSet, 
   dataset0, PairwiseMutantGRCGenotypes, RefFitness, 
   FitnessEffects4PairwiseUpDownPerturbs, 
   RelFitnessEffects4PairwiseUpDownPerturbs, dataset1},
  
  PairwiseGeneGeneIntPertSet = 
   GeneratePairwiseGeneGeneIntPertSet2[FocalParamSet, 
    MutationalStep];
  
  PairwiseGeneGeneIntPertDataSet = 
   Prepend[PairwiseGeneGeneIntPertSet[[All, 1 ;; 4]], {"MatrixID", 
     "WInteraction", "Perturbation", "Magnitude"}];
  
  dataset0 = 
   With[{header = First@PairwiseGeneGeneIntPertDataSet}, 
     AssociationThread[header -> #] & /@ 
      Rest@PairwiseGeneGeneIntPertDataSet] // Dataset;
  
  PairwiseMutantGRCGenotypes = 
   ReplacePart[FocalParamSet, 1 -> #] & /@ 
    PairwiseGeneGeneIntPertSet[[All, -1]];
  
  RefFitness = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM[
    FocalParamSet];
  
  FitnessEffects4PairwiseUpDownPerturbs = 
   AssessFS4GRC2GenerateStripedPatternInducedBySSMorpGradient4SFGRM /@
     PairwiseMutantGRCGenotypes;
  
  RelFitnessEffects4PairwiseUpDownPerturbs = 
   Abs[FitnessEffects4PairwiseUpDownPerturbs - RefFitness]/
    RefFitness;
  
  dataset1 = 
   dataset0[
    MapThread[
      Append[#1, "RelFitnessChange" -> #2] &, {#, 
       RelFitnessEffects4PairwiseUpDownPerturbs}] &];
  dataset1]



