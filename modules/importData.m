(*Package to import data
Author: Juan Camilo Arboleda Rivera
Last modified: 2020-05-04*)

importGRNs::usage="importGRNs[  ] importa la lista de GRNs clasificadas no
isomórficas -> {GRN, número de topología}."

importGRNs[]:= Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[NotebookDirectory[]],"data","ordered-non-isomorph-grns.txt"}],"TSV"]

importTopologies::usage="importTopologies[  ] importa la lista de topologías
ordenadas según su respectiva abundancia -> {topología, ranking de abundancia}."

importTopologies[]:= Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[NotebookDirectory[]],"data","topologies.txt"}],"TSV"]

importAbundanceTopol::usage="importAbundanceTopol[  ] importa la lista de
topologías y su abundancia -> {topología, abundancia}."

importAbundanceTopol[]:= Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[NotebookDirectory[]],"data","topologies-with-abundances.txt"}],"TSV"]

importProcrustedGRNs::usage="importProcrustedGRNs[  ]"

importProcrustedGRNs[]:=Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[NotebookDirectory[]],"data","procrusted-grns.txt"}],"TSV"]
