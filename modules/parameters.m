BeginPackage["parameters`"]

NumNuclei::usage="";
alpha::usage="";
sigmoidThreshold::usage="";
IntTime::usage="";
OptimalPattern::usage="";
Dmax::usage="";

Begin[ "`Private`"]

(*Global parameters to be used in the algorithm*)
NumNuclei = 30;
alpha = 5.;
sigmoidThreshold = 1;
IntTime = 500.;

(*Optimal expression pattern represented in a scale between 1-10 in
expression level, being 10 an expression level which is > 90% of the
maximal level observed along the 1D field of cells for the output genes*)
OptimalPattern =
 Flatten[{ConstantArray[1, Quotient[NumNuclei, 3]],
   ConstantArray[10, Quotient[NumNuclei, 3] + Mod[NumNuclei, 3]],
   ConstantArray[1, Quotient[NumNuclei, 3]]}];

(*Maximal discrepancy achievable for a given pattern w.r.t to the optimal pattern above*)
Dmax=9*30;

End[]
EndPackage[]
