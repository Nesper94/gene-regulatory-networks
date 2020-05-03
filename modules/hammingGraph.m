(*La idea con este módulo es generar el grafo de distancias de Hamming pero uniendo topologías separadas hasta por dos pasos mutacionales*)
importClassifiedTopologies::usage = "importClassifiedTopologies[] importa la lista de topologías de red con su abundancia."
importClassifiedTopologies[]:= {
  Flatten[ToExpression[#], 1] & /@
   Import["/home/nesper94/Documentos/Biologia_de_sistemas/RobustnessModularityProject/RobustnessModularityProject/motivos_de_red.txt", "TSV"];
}

isomorfos::usage = "isomorfos[matriz] genera los isomorfos de una matriz dada."
(*Definimos la función que presenta los isomorfos de cada topología*)
isomorfos[matriz_] := Module[{ordfilas, ordcolumnas},
  ordfilas = Permutations[{1, 2, 3}];
  ordcolumnas = Table[Prepend[ordfilas[[i]] + 1, 1], {i, Length[ordfilas]}];
  Table[ matriz[[ordfilas[[i]], ordcolumnas[[i]] ]], {i, Length[ordfilas]}] ]

(*La siguiente es la misma función de antes que calcula la menor distancia de Hamming entre dos matrices, pero con la nueva parte del código*)
jaminnn::usage = "jaminnn[list1, list2] compara los isomorfos de las dos listas y arroja la menor distancia de Hamming entre ellas."
jaminnn[list1_, list2_] := Module[{comp},
  comp = Flatten[#, 1] & /@ Table[ {{ (Flatten[#]&/@ isomorfos[list1])[[i]] }, {Flatten[list2]} }, {i, Length[ (Flatten[#]&/@isomorfos[list1]) ] } ] ; (*Esta línea me crea tuplas compuestas por la lista 2 con cada uno de los isomorfos de la lista 1*)
  Min[Apply[HammingDistance, #]&/@ comp (*Esta línea me calcula la distancia de Hamming a cada una de las tuplas anteriores y luego me elige la menor distancia*)
   ]]

calcularSimilitudesNoClasificadas::usage = "calcularSimilitudes[lista_de_topologías] da como resultado una lista con todas las posibles comparaciones entre pares de topologías y las distancias de Hamming entre ellas. Usa como input una lista de topologías NO clasificadas."
calcularSimilitudesNoClasificadas[topologies_List]:= Module[{},
  Table[{a, b, jaminnn[topologies[[a]], topologies[[b]]]}, {a, 1, Length[topologies]}, {b, a + 1, Length[topologies]}] ]

calcularSimilitudesClasificadas::usage = "calcularSimilitudes[lista_de_topologías] da como resultado una lista con todas las posibles comparaciones entre pares de topologías y las distancias de Hamming entre ellas. Usa como input una lista de topologías clasificadas."
calcularSimilitudesClasificadas[topologies_List]:= Module[{topologies},
  Table[{a, b, jaminnn[topologies[[a,1]], topologies[[b,1]]]}, {a, 1, Length[topologies]}, {b, a + 1, Length[topologies]}] ]

vecinos::usage = "vecinos[similitudes, distancia] de la lista de distancias de Hamming entre topologías, filtra aquellas que estén a la distancia establecida o menor."
vecinos[similitudes_List, distancia_]:= Module[{},
  Cases[DeleteCases[SortBy[Flatten[similitudes, 1], Last], {_, _, c_} /; c > distancia], {a_, b_, c_} :> {a, b}] ]
