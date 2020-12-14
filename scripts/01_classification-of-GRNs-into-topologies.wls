#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*SetDirectory[ NotebookDirectory[]];*)
(*<<"../modules/importData.m"*)


(* ::Input:: *)
(*(*Import genotypes*)*)
(*matrices=Flatten[ToExpression[#],1]&/@Import["../data/raw-grns.txt","TSV"];*)


(* ::Input:: *)
(*(*Define matrices indicating row and column order*)*)
(*rowOrder=Permutations[{1,2,3}];*)
(*(*columnOrder must be equal to rowOrder but with morphogen's column always as first column*)*)
(*columnOrder=Table[*)
(*Prepend[rowOrder[[i]]+1,1]*)
(*,{i,Length[rowOrder]}];*)
(**)
(*(*Get Sign for each one of the matrices generating a new list*)*)
(*elimination=Sign[First[#]]&/@matrices;*)
(**)
(*(*Create empty list that will contain the topologies and its frequency*)*)
(*motifs={} ;*)


(* ::Input:: *)
(*While[Length[elimination]>0,{*)
(*(*Choose one of the matrices*)*)
(*matrizAleatoria=elimination[[   RandomInteger[{1,Length[elimination]}]   ]] ,*)
(**)
(*(*Generate matrix's isomorphs*)*)
(*isomorphs=Table[*)
(*matrizAleatoria[[ rowOrder[[i]],columnOrder[[i]] ]]*)
(*,{i,Length[rowOrder]}] ,*)
(**)
(*(*Count the number of time a matrix or its isomorphs appears in the list of matrices*)*)
(*frecMatriz=Count[elimination,  Alternatives@@isomorphs],*)
(**)
(*(*Eliminamos la matriz y sus isomorphs de la lista de matrices*)*)
(*elimination=DeleteCases[elimination,Alternatives @@isomorphs] ,*)
(**)
(*(*Add the matrix and its frequency to the list "motifs"*)*)
(*AppendTo[motifs,{isomorphs[[1]],frecMatriz} ]*)
(*}   ];*)


(* ::Input:: *)
(*(*Sort list by frequency*)*)
(*motifs=Sort[motifs,#1[[2]]>#2[[2]]&];*)


(* ::Input:: *)
(*(*Export results in "motifs" to a text file:*)*)
(*Export["../data/topologies-with-abundances.txt",motifs]*)
