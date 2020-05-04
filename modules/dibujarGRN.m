(*2019.03.07: Voy a cambiar la presentación de los gráficos de manera que la
inhibición se represente como una flecha con cabeza plana, según la convención
general. Este módulo es básicamente la función
DisplayNonWeightedHaploidGraphWithMorphInputs del módulo
EvolAlgorithm4PatternFormingGRCModel2MathPackageV9 pero con las modificaciones
necesarias para la nueva representación.*)
ColorF01[Node_] := Switch[Node,"A",Darker[Gray,0.2],"B",Darker[Gray,0.2],
    "C",Darker[Gray,0.2],"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"A\",\nFontSize->0]\)\)]\)",Black,"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"B\",\nFontSize->0]\)\)]\)",Black,"\!\(\*SubscriptBox[\(M\), \(\!\(\*
StyleBox[\"C\",\nFontSize->0]\)\)]\)",Black];
EdgeColoringF::usage="Función que hace que una interacción de represión (-1)
sea roja y una de activación (1) sea azul."
EdgeColoringF[RegEffect_]:=Switch[RegEffect,1,Blue,-1,Red]

(*La siguiente es la función creada para generar la flecha de inhibición*)
ArrowheadStyleFunction[RegEffect_] :=
 Switch[RegEffect, 1, Arrowheads[Automatic], -1,
  Arrowheads[{{0.25, Automatic, Graphics[Line[{{0, 1/12}, {0, -1/12}}]]  }}]  ]

dibujarGRN[mKVsMatrix_List]:=Block[{EdgeComponents,Edge2ColorMap,ArrowHeadMap},
HapRegMotifs=Prepend[Partition[(#[[1]]-> #[[2]])&/@Tuples[{"A","B","C"},2],3],
{"\!\(\*SubscriptBox[\(M\), \(\!\(\*
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
{EdgeComponents,Edge2ColorMap}={({#[[1]],#[[2]]}&/@GraphRules),
(EdgeColoringF/@RegEffectList)};

(*Aplicamos la función a la lista de efectos de regulación*)
ArrowHeadMap = ArrowheadStyleFunction /@ RegEffectList;

CombinedGraphRules=Thread[{Flatten[GraphRules],Abs[Sign[Tags]]*0.0095}];

(*GraphPlot[CombinedGraphRules,VertexLabeling->True,DirectedEdges->{True,
"ArrowheadsSize"->0.0001},*)
GraphPlot[CombinedGraphRules,VertexLabeling->True,DirectedEdges->{True,
    "ArrowheadsSize"->0.0001},
SelfLoopStyle->0.25,VertexCoordinateRules->{MonA,NodeA,MonB,NodeB,MonC,NodeC},
(*2019.03.07: Voy a cambiar la presentación de los gráficos de manera que la
inhibición se represente como una flecha con cabeza plana, según la convención
general. Bajo esta línea se encuentra la línea original en forma de comentario,
la línea que le sigue es la línea modificada*)
(*EdgeRenderingFunction->({Edge2ColorMap[[Flatten[Position[
EdgeComponents,#2]][[1]]]],Thickness[#3],Arrow[#1,0.15]}&),*)
EdgeRenderingFunction->({
  ArrowHeadMap[[  Flatten[ Position[EdgeComponents, #2]   ][[1]]   ]],
  Edge2ColorMap[[Flatten[Position[EdgeComponents,#2]][[1]]]],Thickness[#3],
  Arrow[#1,0.2]}&),
Method->{"Automatic","Rotation"->2 Pi},VertexRenderingFunction->
({EdgeForm[ColorF01[#2]],ColorF01[#2],Disk[#1,0.15],Darker[Gray,0.009],
    Text[Style[#2,13,Bold,White],#1]}&),ImageSize->300]
]

(******************************************************************************)
(*Esta es la función para graficar las topologías teniendo en cuenta la fuerza
de las interacciones*)

dibujarGRNfuerza[mKVsMatrix_List]:=Block[{},
HapRegMotifs=Prepend[Partition[(#[[1]]-> #[[2]])&/@Tuples[{"A","B","C"},2],3],
{"\!\(\*SubscriptBox[\(M\), \(\!\(\*
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
{EdgeComponents,Edge2ColorMap}={({#[[1]],#[[2]]}&/@GraphRules),
(EdgeColoringF/@RegEffectList)};
(*La siguiente línea fue insertada por JCAR el 2019-03-15*)
ArrowHeadMap = ArrowheadStyleFunction /@ RegEffectList;
MaxThickness=0.0095;
minKV=Min[Abs[Tags]];maxKV=Max[Abs[Tags]];
NormalizedKVs=Abs[((#-minKV)/(maxKV-minKV))]&/@Tags;
ScaledKVs=If[#<=0.01,(MaxThickness*0.95),MaxThickness*#]&/@NormalizedKVs;
CombinedGraphRules=Thread[{Flatten[GraphRules],ScaledKVs}];

GraphPlot[CombinedGraphRules,VertexLabeling->True,DirectedEdges->
{True,"ArrowheadsSize"->0.0001},
SelfLoopStyle->0.25,VertexCoordinateRules->{MonA,NodeA,MonB,NodeB,MonC,NodeC},
EdgeRenderingFunction->({
ArrowHeadMap[[  Flatten[ Position[EdgeComponents, #2]   ][[1]]   ]],
(*En la siguiente línea se cambia Arrow[#1,0.15] por Arrow[#1,0.2] para que la
punta de las flechas no sea tapada por los nodos.*)
  Edge2ColorMap[[Flatten[Position[EdgeComponents,#2]][[1]]]],Thickness[#3],
  Arrow[#1,0.2]}&),
Method->{"Automatic","Rotation"->2 Pi},VertexRenderingFunction->
({EdgeForm[ColorF01[#2]],ColorF01[#2],Disk[#1,0.15],Darker[Gray,0.009],
    Text[Style[#2,13,Bold,White],#1]}&),ImageSize->300]
]
