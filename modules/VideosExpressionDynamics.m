(*Paquete para analizar la dinámica de expresión de cada genotipo.
This package depends on EvolAlgorithm4PatternFormingGRCModel2MathPackageV9*)
<<EvolAlgorithm4PatternFormingGRCModel2MathPackageV9`

videoExp::usage = "videoExp[ genotipo ] crea un video de la dinámica de expresión de un genotipo, acepta como argumentos 'gen', con las opciones 'a','b' y 'c', tamaño del gráfico y 'fps' (cuadros por segundo)."

(*Esta línea importa la lista de motivos clasificados, ordenados no isomorfos*)
motivosClasificadosOrdenadosNoIsomorfos = Flatten[ToExpression[#],1]&/@Import["/home/nesper94/Documentos/biologia-de-sistemas/RobustnessModularityProject/RobustnessModularityProject/motivos_clasificados_ordenados_no_isomorfos.txt","TSV"];

(*Esta función me permite saber en qué posiciones están los genotipos que pertenecen a determinado motivo de red*)
posiciónMotivos[motivo_Integer]:=Flatten[{First[#],Last[#]}&/@{Position[motivosClasificadosOrdenadosNoIsomorfos[[All,2]],motivo]}]

(*Aquí redefino la función "SpaceTimeFeats4StripeFormingGRCs4SSMorpGradientWithSumAndFilterModel0" de tal manera que acepte los genotipos sin necesidad de usar como argumento el input de morfógeno*)
(*Nota [2018.11.08]: Esta función no concuerda totalmente con los resultados de "SolveGRCDynsBeforeAfterMorphInput", por lo tanto debe ser revisada y modificada*)
(*Nota [2018.11.16]: La función fue modificada y ahora concuerda con los resultados de "SolveGRCDynsBeforeAfterMorphInput"*)
expEspacioTemporal[{Wmatrix_List, DiffParams_List, DegParams_List},ti_:0,tf_:300,step_:6] :=
  Module[
  {Prot, ProtStateVariab, ProtInitExpStat, DynSyst, Sol, GRCPhenotReadout},

  (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
  (*Transfer function for mapping integrated regulatory inputs to transcriptional outputs. In this model the SigmoidSteepness = alpha = 5, while the parameter b, which gives the location of the threshold value, is set to 1 *)
  ICs=ConstantArray[1,90];
  NullMorpInput=ConstantArray[0,30];
  {PreMorpInputFS,SSExpValuesPreMorpInput} = AssessFitnessScore4StripePattern4SSMorpGradient4SFGRM[ICs,NullMorpInput,{Wmatrix, DiffParams, DegParams}];

  SumAndFilterF[SigmoidSteepness_, Threshold_, IntegratedRegInput_] := 1./(1. + Exp[(SigmoidSteepness - (Threshold*IntegratedRegInput))]);
  (*A function for summing the weigthed regulatory inputs*)

  AdditiveRegContributionF[a_,n_,t_,W_List,Morph_] := Total[W[[a]]*Prepend[Table[Prot[b, n, t], {b,1,GRNSize}], Morph]];
  (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
  MorphInput = SSInputMorphogen[1];  (*Esta es la línea de código que agregué yo (J.C.A.R)*)

  GRNSize = Length[Wmatrix];
  (*As set in the paper referenced above*)
  NumNuclei = 30;
  (*As set in the paper referenced above*)
  alpha = 5.;
  sigmoidThreshold = 1;
  IntTime = 500.;
  (*Initial conditions for all variables, including boundary conditions set = 0*)

  (*Línea original de este código*)
  (*ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];*)

  (*Línea proveniente de 'StripeFormingGRCs4SSMorpGradientSFGRM'*)
  ProtStateVariab=Flatten[{Table[Table[Prot[a,n,t],{a,1,GRNSize}],{n,NumNuclei}],Flatten[Table[Table[Prot[a,n,0],{a,1,GRNSize}],{n,{0,NumNuclei+1}}]]}];

  (*Línea original de este código*)
  (*ProtInitExpStat=Flatten[{Table[Table[Prot[a,n,0]== 0.1, {a, 1, GRNSize}], {n, NumNuclei}],Flatten[Table[Table[Prot[a, n, 0] == 0.00, {a, 1, GRNSize}], {n, {0,   NumNuclei + 1}}]]}];*)

  (*Línea proveniente de 'StripeFormingGRCs4SSMorpGradientSFGRM'*)
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

  Sol = NDSolve[Flatten[{DynSyst,ProtInitExpStat}],ProtStateVariab,{t,0.,IntTime},Method->{"EquationSimplification"->"Residual"}];

  (*La siguiente es la línea sacada de 'StripeFormingGRCs4SSMorpGradientSFGRM'*)
  (*GRCPhenotReadout= Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,IntTime-100,IntTime,2}],{i,3}]*)

  (*Esta es la línea original de este código*)
  GRCPhenotReadout= Transpose[Thread[Table[Flatten[Table[Flatten[Evaluate[Prot[i,n,t]/.Sol/.t->#]],{n,NumNuclei}]]&/@Table[x,{x,Range[ti,tf,step]}], {i, 3}]]]
  ]

(**************************************************************************************************************)

(*Aquí defino la función que me crea los videos de la dinámica de expresión*)
videoExp[ genotipo_List, OptionsPattern[]] := Module[{patron},
  Options[videoExp] = {gen -> 0, tamaño -> UpTo[300], PlotRange -> {0,25}, fps -> Automatic};
    patron = expEspacioTemporal[ genotipo[[1]] ];
  If[ SubsetQ[ {a, b, c}, Flatten[{ OptionValue[gen] }] ],
   ListAnimate[
    Table[ ListPlot[ Table[
       Prepend[{patron[[ # /. {a -> 1, b -> 2, c -> 3} ]] [[n]] [[i]]}, i]
       , {i, 30}], Joined -> True, PlotRange -> OptionValue[PlotRange],
      ImageSize -> OptionValue[tamaño] ], {n, 30}], OptionValue[fps] ] &/@ Flatten[{OptionValue[gen]} ],
   Table[
    ListAnimate[
     Table[ ListPlot[ Table[ Prepend[ {patron[[k]][[n]][[i]]}, i]
        , {i, 30}], Joined -> True, PlotRange -> OptionValue[PlotRange], ImageSize -> OptionValue[tamaño] ], {n, 30}], OptionValue[fps] ]
    , {k, 3} ] ] ]

generarIsomorfos::usage = "generarIsomorfos[ genotipo(matriz) ] Esta función genera los diferentes isomorfos de una matriz en forma de matrices signo que indican las interacciones entre los nodos en las matrices de adyacencia."

generarIsomorfos[genotipo_List]:= Module[{ordfilas,ordcolumnas,isomorfos},
  (*Definimos las matrices que nos indican el orden de las filas y las columnas de las matrices*)
  ordfilas = Permutations[{1, 2, 3}];
  (*ordcolumnas debe ser igual a ordfilas pero con la columna del morfógeno siempre como primera columna*)
  ordcolumnas = Table[ Prepend[ordfilas[[i]] + 1, 1] , {i, Length[ordfilas]}];
  (*Generamos los isomorfos topológicos de la matriz*)
  isomorfos = Table[
    Sign[genotipo][[ ordfilas[[i]], ordcolumnas[[i]] ]]
    , {i, Length[ordfilas]}] ]

isSubgraphQ::usage = "isSubgraphQ[ lista_1, lista_2 ] Da Verdadero si lista_1 es un subgrafo de lista_2 y Falso si lo contrario, i.e. dice si un patrón de interconexiones entre nodos determinado se encuentra en una topología de red. Las listas deben ser matrices de adyacencia."

isSubgraphQ[subgraph_List,graph_List]:= Module[{isomorfos},
  isomorfos = generarIsomorfos[subgraph];
  Do[ If[ i > Length[isomorfos], Break[False] ];
      If[
    SubsetQ[ Position[ Sign[graph],1], Position[isomorfos[[i]],1] ] && SubsetQ[ Position[ Sign[graph],-1], Position[isomorfos[[i]],-1] ], Break[True] ];
    ,{i, Length[isomorfos]+1}
    ]]

(********************************************************************************************************************)

(*Aquí quiero redefinir la función "AssessFitnessScoreStripeFormingGRCs4SSMorpGradient" para que sea posible usarla sin necesidad de especificar el input de morfógeno*)
AssessFitnessScoreStripeFormingGRCs4SSMorpGradientJCAR[{Wmatrix_List,DiffParams_List,DegParams_List}]:=Block[
{GRCPhenotReadout},
MorphInput = SSInputMorphogen[1];  (*Esta es la línea de código que agregué yo (J.C.A.R)*)
GRCPhenotReadout = StripeFormingGRCs4SSMorpGradientWithSumAndFilterGeneRegModel[{Wmatrix,MorphInput,DiffParams,DegParams}];

(*This condition here is to check if there is any negative gene expression value, which is not biologically realistic. If there is, then Fitness must be 0. *)

If[
	Count[Flatten[GRCPhenotReadout],x_/;x<=0]>0,0,
	FitnessF4SingleStripe[GRCPhenotReadout]
	]
		]
