(*Paquete para importar datos.*)

importGRNs::usage="importGRNs[  ] importa la lista de GRNs clasificadas no isomórficas -> {GRN, número de topología}."

importGRNs[]:= Flatten[ToExpression[#],1]&/@Import["/home/nesper94/Documentos/biologia-de-sistemas/RobustnessModularityProject/RobustnessModularityProject/motivos_clasificados_ordenados_no_isomorfos.txt","TSV"]

importTopologies::usage="importTopologies[  ] importa la lista de topologías ordenadas según su respectiva abundancia -> {topología, ranking de abundancia}."

importTopologies[]:= Flatten[ToExpression[#],1]&/@Import["/home/nesper94/Documentos/biologia-de-sistemas/RobustnessModularityProject/RobustnessModularityProject/lista_de_topologías.txt","TSV"]

importAbundanceTopol::usage="importAbundanceTopol[  ] importa la lista de topologías y su abundancia -> {topología, abundancia}."

importAbundanceTopol[]:= Flatten[ToExpression[#],1]&/@Import["/home/nesper94/Documentos/biologia-de-sistemas/RobustnessModularityProject/RobustnessModularityProject/motivos_de_red.txt","TSV"]

importProcrustedGRNs::usage="importProcrustedGRNs[  ]"

importProcrustedGRNs[]:=Flatten[ToExpression[#],1]&/@Import["/home/nesper94/Documentos/biologia-de-sistemas/RobustnessModularityProject/RobustnessModularityProject/procrusted_GRNs.txt","TSV"]
