const variable = [
  "Empty"
  "E"
  "H"
  "E1"
  "E2"
  "E3"
  "E4"
  "E_E1"
  "H_E3"
  "H_E4"
  "E11"
  "E12"
  "E23"
  "E34"
  "E24"
  "E44"
  "E11P"
  "E12P"
  "E23P"
  "E34P"
  "E24P"
  "E44P"
  "G"
  "S"
  "I"
  "R"
  "O"
  "A"
  "E11G"
  "E11S"
  "E11R"
  "E12G"
  "E12S"
  "E12R"
  "E23G"
  "E23S"
  "E23I"
  "E23R"
  "E34G"
  "E34S"
  "E34I"
  "E34R"
  "E24G"
  "E24S"
  "E24I"
  "E24R"
  "E44G"
  "E44S"
  "E44I"
  "E44R"
  "sigmaG"
  "sigmaS"
  "sigmaI"
  "sigmaR"
  "sigmaA"
  "sigmaSP"
  "sigmaAP"
  "sigmaG_O"
  "sigmaG_A"
  "sigmaSP_G"
  "sigmaAP_S"
  "sigmaAP_I"
  "sigmaAP_R"
  "P3_A"
  "P2"
  "P3"
  "Akt"
  "Aktstar"
  "RsD"
  "RsT"
  "sigmaRP"
  "Raf"
  "Rafstar"
  "MEK"
  "MEKstar"
  "ERK"
  "ERKstar"
  "OP"
  "AP"
  "A_sigmaG_O"
  "sigmaA_G"
  "sigmaA_G_O"
  "sigmaO"
  "E13"
  "E14"
  "E13P"
  "E14P"
  "E13G"
  "E13S"
  "E13I"
  "E13R"
  "E14G"
  "E14S"
  "E14I"
  "E14R"
  "fint"
  "T"
  "E11T"
  "E12T"
  "E23T"
  "E34T"
  "E24T"
  "E44T"
  "E13T"
  "E14T"
  "sigmaAP_T"
  "sigmaT"
  "E1_PT"
  "E2_PT"
  "E4_PT"
  "E_E1_PT"
  "H_E4_PT"
  "pERK"
  "ERK_MEKstar"
  "pERK_MEKstar"
  "ERKPpase"
  "ERKstar_ERKPpase"
  "pERK_ERKPpase"
];

for (index,value) in enumerate(variable)
  eval(Meta.parse("const $value=$index"));
end