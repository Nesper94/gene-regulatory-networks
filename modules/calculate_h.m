calculateMorphDecay::usage="calculateMorphDecay[NumNuclei_:30, A0_:1] calculates the value of the morphogen decay parameter (h) for morphogenetic fields with different numbers of nuclei given by NumNuclei so that the concentration in the last cell is equal to A0*Exp[-29/(30*0.4)]"
calculateMorphDecay[NumNuclei_:30, A0_:1] := Block[
    {rightEndConcentration, h},

  rightEndConcentration = A0*Exp[-29/(30*0.4)];

  h = -(NumNuclei - 1)/(NumNuclei*Log[(rightEndConcentration/A0)]);

  Return[h]
]
