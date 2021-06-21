(* ::Package:: *)

BeginPackage["GraphUtils`"]

getIsomorphs::usage =
    "Generates the different isomorphs of a graph as adjacency matrices."

isSubgraphQ::usage =
    "Gives True if graph1 is "<>
    "subgraph of graph2 and False in other case. Both lists must be "<>
    "adjacency matrices."

isGraphQ::usage =
    "isGraphQ[ graph1_List, graph2_List ] Gives True if graph1 is the "<>
    "same graph as graph2 (or an isomorph) and False in other case. "<>
    "Both lists must be adjacency matrices."

getSupergraphs::usage =
    "getSupergraphs[topology_List, topologies_List] Returns a list of the "<>
    "positions of topologies that are supergraphs of topology."

findGraph::usage =
    "findGraph[topology_List, topologies_List] Returns a list of the "<>
    "positions of topologies that are isomorphic with topology."

Begin["Private`"]

getIsomorphs[genotipo_List] :=
 Module[{rowOrder, columnOrder, isomorphs},
  (*Define matrices indicating the order of rows and columns in the
    matrices*)
  rowOrder = Permutations[{1, 2, 3}];

  (*columnOrder must be equal to rowOrder but with the morphogen
    always as the first column*)
  columnOrder =
   Table[Prepend[rowOrder[[i]] + 1, 1], {i, Length[rowOrder]}];

  (*Generate isomorphs for the matrix*)
  isomorphs =
   Table[Sign[genotipo][[rowOrder[[i]], columnOrder[[i]]]], {i,
     Length[rowOrder]}]
  ]

isSubgraphQ[graph1_List, graph2_List] :=
 Module[{isomorphs},
  isomorphs = getIsomorphs[graph1];
  Do[If[i > Length[isomorphs], Break[False]];
   If[SubsetQ[Position[Sign[graph2], 1], Position[isomorphs[[i]], 1]] &&
      SubsetQ[Position[Sign[graph2], -1],
      Position[isomorphs[[i]], -1]], Break[True]];, {i,
    Length[isomorphs] + 1}]
  ]

isGraphQ[graph1_List, graph2_List] :=
 Module[{isomorphs},
  isomorphs = getIsomorphs[graph1];
  Do[If[i > Length[isomorphs], Break[False]];
   If[graph2 == isomorphs[[i]], Break[True]];
   , {i, Length[isomorphs] + 1}]
  ]

getSupergraphs[topology_List, topologies_List] :=
 Module[{supergraphs},
  supergraphs = {};
  Do[ If[ isSubgraphQ[topology, topologies[[i, 1]] ],
    AppendTo[supergraphs, i] ]
   , {i, Length[topologies]}];
  Return[supergraphs]
  ]

findGraph[topology_List, topologies_List] :=
 Module[{positions},
  positions = {};
  Do[
   If[isGraphQ[topology, topologies[[i, 1]]], AppendTo[positions, i]]
   , {i, Length[topologies]}];
  Return[positions]
  ]

End[]

EndPackage[]
