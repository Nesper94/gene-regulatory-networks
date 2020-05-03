(*WARNING: This module is deprecated. For the new version see fitness2.m*)

fitnessGenes[ {Wmatrix_List, DiffParams_List, DegParams_List} ]:=Block[{GRCPhenotReadout, MorphInput},
	(*La siguiente línea me genera el GRCPhenotReadout*)
	MorphInput = SSInputMorphogen[1];
	GRCPhenotReadout = StripeFormingGRCs4SSMorpGradientWithSumAndFilterGeneRegModel[{Wmatrix,MorphInput,DiffParams,DegParams}];
	(*GRCPhenotReadout*)
	(************************************************************************************************)

	(*Testing for stability of the expression pattern after 500 time steps*)
	EF=((1/3)*Total[EquilibFilter/@GRCPhenotReadout]);
	(*EF*)

	(*Testing for spatial heterogeneity in the steady state expression pattern,
	which is also an indicative of whether the expression profiles are
	sufficiently high so that the cross regulatory interactions can be
	effective in controlling the expression putput. Note that unlike the
	FitnessF implemented before, in this FitnessF2 we only checked for the
	heterogeineity in the output node, and not in the whole set of nodes
	considered*)
	PF = PatternFilter /@ GRCPhenotReadout;
	(*PF*)

	FS=(AssessExpPattern4SingleStripe[#]*Q1[Min[PF]]*Q2[EF])&/@GRCPhenotReadout
	(*{Q1[Min[PF]], Q2[EF], AssessExpPattern4SingleStripe[ GRCPhenotReadout[[1]] ]}*)

	]


(*La siguiente función fue copiada de EvolAlgorithm4PatternFormingGRCModel2MathPackageV9. No se por qué, pero tiene que estar en el mismo paquete que la función de arriba para que esta funcione.*)

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
