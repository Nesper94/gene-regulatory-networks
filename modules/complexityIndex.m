(*Módulo para calcular el índice de complejidad*)

vertexDegrees::usage="vertexDegrees[ adjacencyMatrix ] calcula la cantidad de enlaces salientes a partir de una matriz de adyacencia en un grafo dirigido."
vertexDegrees[adjacencyMatrix_] := Total[ Abs[Sign[ adjacencyMatrix ]] ]

complexityIndex::usage="complexityIndex[ adjacencyMatrix ] calcula el índice de complejidad de la red representada por una matriz de adyacencia."
complexityIndex[adjacencyMatrix_] := Module[{transp, new},
  transp = Prepend[vertexDegrees[Transpose[adjacencyMatrix]], 0];
  new = vertexDegrees[adjacencyMatrix] + transp;
  N[Sum[ Subscript[a, i]*Log[2, Subscript[a, i]],
      {Subscript[a, i], Cases[ new, Except[0] ]} ]] /; MatrixQ[adjacencyMatrix] ]
