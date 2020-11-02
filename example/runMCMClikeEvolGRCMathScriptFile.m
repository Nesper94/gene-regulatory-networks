SetDirectory[DirectoryName[$InputFileName]];

SeedRandom[RandomInteger[10^7, 1][[1]]]

<<EvolAlgorithm4PatternFormingGRCModel2MathPackageV9`
<< DesignMorphogeneResponsiveGRCs`

EngineerMorphInducibleGRCGenotypes4SFGRM[ToExpression[Repl],ToExpression[EvolSteps],ToExpression[ParamSettingSamplingRate],ToExpression[thresholdFitness]]

Exit[];
