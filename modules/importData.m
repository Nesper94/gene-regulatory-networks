(*Package to import data
Author: Juan Camilo Arboleda Rivera
Last modified: 2020-05-14*)

(* Get module path *)
path = $InputFileName /. "" :> NotebookFileName[];
(* parentDir = DirectoryName @ parentPath; *)

importGRNs::usage="importGRNs[  ] importa la lista de GRNs clasificadas no
isomórficas -> {GRN, número de topología}."

importGRNs[]:= Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[ DirectoryName @ path ],"data","ordered-non-isomorph-grns.txt"}],"TSV"]

importTopologies::usage="importTopologies[  ] importa la lista de topologías
ordenadas según su respectiva abundancia -> {topología, ranking de abundancia}."

importTopologies[]:= Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[ DirectoryName @ path ],"data","topologies.txt"}],"TSV"]

importAbundanceTopol::usage="importAbundanceTopol[  ] importa la lista de
topologías y su abundancia -> {topología, abundancia}."

importAbundanceTopol[]:= Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[ DirectoryName @ path ],"data","topologies-with-abundances.txt"}],"TSV"]

importProcrustedGRNs::usage="importProcrustedGRNs[  ]"

importProcrustedGRNs[]:=Flatten[ToExpression[#],1] &/@ Import[ FileNameJoin[
{ParentDirectory[ DirectoryName @ path ],"data","procrusted-grns.txt"}],"TSV"]
