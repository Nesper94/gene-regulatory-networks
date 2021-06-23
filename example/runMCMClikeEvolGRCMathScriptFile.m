SetDirectory[DirectoryName[$InputFileName]];

SeedRandom[RandomInteger[10^7, 1][[1]]]

<< "../modules/EvolAlgorithm4PatternFormingGRCModel2MathPackageV9.m"
<< "../modules/DesignMorphogeneResponsiveGRCs.m"
<< "../modules/parameters.m"

NumNuclei = uNumNuclei

EngineerMorphInducibleGRCGenotypes4SFGRM[ToExpression[Repl],ToExpression[EvolSteps],ToExpression[ParamSettingSamplingRate],ToExpression[thresholdFitness]]

Exit[];
