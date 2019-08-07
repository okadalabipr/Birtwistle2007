function diffeq(du,u,p,t)

    # extracellular volume to cytoplasmic volume ratio
    VeVc = 33.3;

    # fraction definitions
    if u[V.sigmaS] + u[V.sigmaSP] + u[V.sigmaSP_G] > 0.0
        fsigmaS = u[V.sigmaS]/(u[V.sigmaS] + u[V.sigmaSP] + u[V.sigmaSP_G]);
    else
        fsigmaS = 0.0;
    end

    if u[V.sigmaG] + u[V.sigmaG_A] + u[V.sigmaG_O] + u[V.A_sigmaG_O] > 0.0
        fsigmaG = u[V.sigmaG]/(u[V.sigmaG] + u[V.sigmaG_A] + u[V.sigmaG_O] + u[V.A_sigmaG_O]);
    else
        fsigmaG = 0.0;
    end

    if u[V.sigmaA] + u[V.sigmaAP] + u[V.sigmaAP_S] + u[V.sigmaAP_R] + u[V.sigmaAP_I] + u[V.sigmaAP_T] > 0.0
        fsigmaA = u[V.sigmaA]/(u[V.sigmaA] + u[V.sigmaAP] + u[V.sigmaAP_S] + u[V.sigmaAP_R] + u[V.sigmaAP_I] + u[V.sigmaAP_T]);
    else
        fsigmaA = 0.0;
    end

    if u[V.sigmaR] + u[V.sigmaRP] > 0.0
        fsigmaR = u[V.sigmaR]/(u[V.sigmaR] + u[V.sigmaRP]);
    else
        fsigmaR = 0.0;
    end

    sigmaEP = u[V.E11P] + u[V.E12P] + u[V.E23P] + u[V.E24P] + u[V.E34P] + u[V.E44P] + u[V.E13P] + u[V.E14P];
    if sigmaEP > 0.0
        f11 = u[V.E11P]/sigmaEP;
    else
        f11 = 0.0;
    end

    # E11 degradation
    v10d = p[C.kdeg]*u[V.E11P];
    v16d = p[C.kdeg]*u[V.E11G];
    v17d = p[C.kdeg]*u[V.E11S];
    v18d = p[C.kdeg]*u[V.E11R];
    v73d = p[C.kdeg]*u[V.E11T];

    v::Vector{Float64} = zeros(96);

    v[1] = (p[C.kon1]*u[V.E]*u[V.E1] - p[C.EGF_off]*u[V.E_E1]);
    v[2] = (p[C.kon2]*u[V.H]*u[V.E3] - p[C.HRGoff_3]*u[V.H_E3]);
    v[3] = (p[C.kon3]*u[V.H]*u[V.E4] - p[C.HRGoff_4]*u[V.H_E4]);
    v[4] = (p[C.kon4]*u[V.E_E1]*u[V.E_E1] - p[C.koff4]*u[V.E11]);
    v[5] = (p[C.kon5]*u[V.E_E1]*u[V.E2] - p[C.koff5]*u[V.E12]);
    v[6] = (p[C.kon6]*u[V.H_E3]*u[V.E2] - p[C.koff6]*u[V.E23]);
    v[7] = (p[C.kon7]*u[V.H_E3]*u[V.H_E4] - p[C.koff7]*u[V.E34]);
    v[8] = (p[C.kon8]*u[V.H_E4]*u[V.E2] - p[C.koff8]*u[V.E24]);
    v[9] = (p[C.kon9]*u[V.H_E4]*u[V.H_E4] - p[C.koff9]*u[V.E44]);
    v[10] = (p[C.kf10]*u[V.E11] - p[C.VmaxPY]*u[V.E11P]/(p[C.KmPY] + u[V.E11P]) - p[C.kPTP10]*u[V.sigmaT]*u[V.E11P]);
    v[11] = (p[C.kf11]*u[V.E12] - p[C.VmaxPY]*u[V.E12P]/(p[C.KmPY] + u[V.E12P]) - p[C.kPTP11]*u[V.sigmaT]*u[V.E12P]);
    v[12] = (p[C.kf12]*u[V.E23] - p[C.VmaxPY]*u[V.E23P]/(p[C.KmPY] + u[V.E23P]) - p[C.kPTP12]*u[V.sigmaT]*u[V.E23P]);
    v[13] = (p[C.kf13]*u[V.E34] - p[C.VmaxPY]*u[V.E34P]/(p[C.KmPY] + u[V.E34P]) - p[C.kPTP13]*u[V.sigmaT]*u[V.E34P]);
    v[14] = (p[C.kf14]*u[V.E24] - p[C.VmaxPY]*u[V.E24P]/(p[C.KmPY] + u[V.E24P]) - p[C.kPTP14]*u[V.sigmaT]*u[V.E24P]);
    v[15] = (p[C.kf15]*u[V.E44] - p[C.VmaxPY]*u[V.E44P]/(p[C.KmPY] + u[V.E44P]) - p[C.kPTP15]*u[V.sigmaT]*u[V.E44P]);
    v[16] = (4*p[C.kon16]*u[V.E11P]*u[V.G] - p[C.koff16]*fsigmaG*u[V.E11G]);
    v[17] = (8*p[C.kon17]*u[V.E11P]*u[V.S] - p[C.koff17]*fsigmaS*u[V.E11S]);
    v[18] = (2*p[C.kon18]*u[V.E11P]*u[V.R] - p[C.koff18]*fsigmaR*u[V.E11R]);
    v[19] = (3*p[C.kon19]*u[V.E12P]*u[V.G] - p[C.koff19]*fsigmaG*u[V.E12G]);
    v[20] = (6*p[C.kon20]*u[V.E12P]*u[V.S] - p[C.koff20]*fsigmaS*u[V.E12S]);
    v[21] = (2*p[C.kon21]*u[V.E12P]*u[V.R] - p[C.koff21]*fsigmaR*u[V.E12R]);
    v[22] = (3*p[C.kon22]*u[V.E23P]*u[V.G] - p[C.koff22]*fsigmaG*u[V.E23G]);
    v[23] = (3*p[C.kon23]*u[V.E23P]*u[V.S] - p[C.koff23]*fsigmaS*u[V.E23S]);
    v[24] = (3*p[C.kon24]*u[V.E23P]*u[V.I] - p[C.koff24]*u[V.E23I]);
    v[25] = (2*p[C.kon25]*u[V.E23P]*u[V.R] - p[C.koff25]*fsigmaR*u[V.E23R]);
    v[26] = (4*p[C.kon26]*u[V.E34P]*u[V.G] - p[C.koff26]*fsigmaG*u[V.E34G]);
    v[27] = (3*p[C.kon27]*u[V.E34P]*u[V.S] - p[C.koff27]*fsigmaS*u[V.E34S]);
    v[28] = (4*p[C.kon28]*u[V.E34P]*u[V.I] - p[C.koff28]*u[V.E34I]);
    v[29] = (2*p[C.kon29]*u[V.E34P]*u[V.R] - p[C.koff29]*fsigmaR*u[V.E34R]);
    v[30] = (3*p[C.kon30]*u[V.E24P]*u[V.G] - p[C.koff30]*fsigmaG*u[V.E24G]);
    v[31] = (4*p[C.kon31]*u[V.E24P]*u[V.S] - p[C.koff31]*fsigmaS*u[V.E24S]);
    v[32] = (1*p[C.kon32]*u[V.E24P]*u[V.I] - p[C.koff32]*u[V.E24I]);
    v[33] = (2*p[C.kon33]*u[V.E24P]*u[V.R] - p[C.koff33]*fsigmaR*u[V.E24R]);
    v[34] = (4*p[C.kon34]*u[V.E44P]*u[V.G] - p[C.koff34]*fsigmaG*u[V.E44G]);
    v[35] = (4*p[C.kon35]*u[V.E44P]*u[V.S] - p[C.koff35]*fsigmaS*u[V.E44S]);
    v[36] = (2*p[C.kon36]*u[V.E44P]*u[V.I] - p[C.koff36]*u[V.E44I]);
    v[37] = (2*p[C.kon37]*u[V.E44P]*u[V.R] - p[C.koff37]*fsigmaR*u[V.E44R]);
    v[38] = (p[C.kf38]*u[V.sigmaS]*sigmaEP - p[C.VmaxPY]*u[V.sigmaSP]/(p[C.KmPY] + u[V.sigmaSP]) - p[C.kPTP38]*u[V.sigmaT]*u[V.sigmaSP]);
    v[39] = (p[C.kf39]*u[V.sigmaA]*sigmaEP - p[C.VmaxPY]*u[V.sigmaAP]/(p[C.KmPY] + u[V.sigmaAP]) - p[C.kPTP39]*u[V.sigmaT]*u[V.sigmaAP]);
    v[40] = (p[C.kon40]*u[V.sigmaG]*u[V.O] - p[C.koff40]*u[V.sigmaG_O]);
    v[41] = (p[C.kon41]*u[V.sigmaG]*u[V.A] - p[C.koff41]*u[V.sigmaG_A]*fsigmaA);
    v[42] = (p[C.kon42]*u[V.sigmaSP]*u[V.G] - p[C.koff42]*u[V.sigmaSP_G]*fsigmaG);
    v[43] = (3*p[C.kon43]*u[V.sigmaAP]*u[V.S] - p[C.koff43]*u[V.sigmaAP_S]*fsigmaS);
    v[44] = (3*p[C.kon44]*u[V.sigmaAP]*u[V.I] - p[C.koff44]*u[V.sigmaAP_I]);
    v[45] = (2*p[C.kon45]*u[V.sigmaAP]*u[V.R] - p[C.koff45]*u[V.sigmaAP_R]*fsigmaR);
    v[46] = (p[C.kon46]*u[V.P3]*u[V.A] - p[C.koff46]*u[V.P3_A]*fsigmaA);
    v[47] = (p[C.kf47]*u[V.P3]*u[V.Akt]/(p[C.Kmf47] + u[V.Akt]) - p[C.Vmaxr47]*u[V.Aktstar]/(p[C.Kmr47] + u[V.Aktstar]));
    v[48] = (p[C.kf48]*(1 - u[V.fint]*f11)*u[V.sigmaI]*u[V.P2]/(p[C.Kmf48] + u[V.P2]) - 3*p[C.PTEN]*u[V.P3]/(p[C.Kmr48] + u[V.P3]));
    v[49] = (p[C.kf49]*u[V.sigmaO]*u[V.RsD]/(p[C.Kmf49] + u[V.RsD]) - p[C.kr49]*u[V.sigmaR]*u[V.RsT]/(p[C.Kmr49] + u[V.RsT]) - p[C.kr49b]*u[V.sigmaRP]*u[V.RsT]/(p[C.Kmr49b] + u[V.RsT]) - p[C.kcon49]*u[V.RsT]);
    v[50] = (p[C.kf50]*u[V.sigmaR]*sigmaEP - p[C.VmaxPY]*u[V.sigmaRP]/(p[C.KmPY] + u[V.sigmaRP]) - p[C.kPTP50]*u[V.sigmaT]*u[V.sigmaRP]);
    v[51] = (p[C.kf51]*u[V.RsT]*u[V.Raf]/(p[C.Kmf51] + u[V.Raf]) - p[C.Vmaxr51]*u[V.Rafstar]/(p[C.Kmrb51] + u[V.Rafstar]));
    v[52] = (p[C.kf52]*u[V.Rafstar]*u[V.MEK]/(p[C.Kmf52] + u[V.MEK]) - p[C.Vmaxr52]*u[V.MEKstar]/(p[C.Kmr52] + u[V.MEKstar]));
    # v[53]: No reaction
    v[54] = (p[C.kf54]*u[V.O]*u[V.ERKstar]/(p[C.Kmf54] + u[V.O]) - p[C.Vmaxr54]*u[V.OP]/(p[C.Kmr54] + u[V.OP]));
    v[55] = (p[C.kf55]*u[V.A]*u[V.ERKstar]/(p[C.Kmf55] + u[V.A]) - p[C.Vmaxr55]*u[V.AP]/(p[C.Kmr55] + u[V.AP]));
    # v[56]: No reaction
    v[57] = (p[C.kon57]*u[V.P3_A]*u[V.G] - p[C.koff57]*u[V.sigmaA_G]);
    v[58] = (p[C.kon58]*u[V.sigmaA_G]*u[V.O] - p[C.koff58]*u[V.sigmaA_G_O]);
    v[59] = (p[C.kon59]*u[V.sigmaG_O]*u[V.A] - p[C.koff59]*u[V.A_sigmaG_O]*fsigmaA);
    v[60] = (p[C.kon60]*u[V.sigmaG_A]*u[V.O] - p[C.koff60]*u[V.A_sigmaG_O]);
    v[61] = (p[C.kon61]*u[V.H_E3]*u[V.E_E1] - p[C.koff61]*u[V.E13]);
    v[62] = (p[C.kon62]*u[V.H_E4]*u[V.E_E1] - p[C.koff62]*u[V.E14]);
    v[63] = (p[C.kf63]*u[V.E13] - p[C.VmaxPY]*u[V.E13P]/(p[C.KmPY] + u[V.E13P]) - p[C.kPTP63]*u[V.sigmaT]*u[V.E13P]);
    v[64] = (p[C.kf64]*u[V.E14] - p[C.VmaxPY]*u[V.E14P]/(p[C.KmPY] + u[V.E14P]) - p[C.kPTP64]*u[V.sigmaT]*u[V.E14P]);
    v[65] = (4*p[C.kon65]*u[V.E13P]*u[V.G] - p[C.koff65]*fsigmaG*u[V.E13G]);
    v[66] = (5*p[C.kon66]*u[V.E13P]*u[V.S] - p[C.koff66]*fsigmaS*u[V.E13S]);
    v[67] = (3*p[C.kon67]*u[V.E13P]*u[V.I] - p[C.koff67]*u[V.E13I]);
    v[68] = (2*p[C.kon68]*u[V.E13P]*u[V.R] - p[C.koff68]*fsigmaR*u[V.E13R]);
    v[69] = (4*p[C.kon69]*u[V.E14P]*u[V.G] - p[C.koff69]*fsigmaG*u[V.E14G]);
    v[70] = (6*p[C.kon70]*u[V.E14P]*u[V.S] - p[C.koff70]*fsigmaS*u[V.E14S]);
    v[71] = (1*p[C.kon71]*u[V.E14P]*u[V.I] - p[C.koff71]*u[V.E14I]);
    v[72] = (2*p[C.kon72]*u[V.E14P]*u[V.R] - p[C.koff72]*fsigmaR*u[V.E14R]);
    v[73] = (4*p[C.kon73]*u[V.E11P]*u[V.T] - p[C.koff73]*u[V.E11T]);
    v[74] = (3*p[C.kon74]*u[V.E12P]*u[V.T] - p[C.koff74]*u[V.E12T]);
    v[75] = (2*p[C.kon75]*u[V.E23P]*u[V.T] - p[C.koff75]*u[V.E23T]);
    v[76] = (2*p[C.kon76]*u[V.E34P]*u[V.T] - p[C.koff76]*u[V.E34T]);
    v[77] = (2*p[C.kon77]*u[V.E24P]*u[V.T] - p[C.koff77]*u[V.E24T]);
    v[78] = (2*p[C.kon78]*u[V.E44P]*u[V.T] - p[C.koff78]*u[V.E44T]);
    v[79] = (3*p[C.kon79]*u[V.E13P]*u[V.T] - p[C.koff79]*u[V.E13T]);
    v[80] = (3*p[C.kon80]*u[V.E14P]*u[V.T] - p[C.koff80]*u[V.E14T]);
    v[81] = (p[C.kf81]*u[V.E1]*u[V.ERKstar]/(p[C.Kmf81] + u[V.E1]) - p[C.Vmaxr81]*u[V.E1_PT]/(p[C.Kmr81] + u[V.E1_PT]));
    v[82] = (p[C.kf82]*u[V.E2]*u[V.ERKstar]/(p[C.Kmf82] + u[V.E2]) - p[C.Vmaxr82]*u[V.E2_PT]/(p[C.Kmr82] + u[V.E2_PT]));
    v[83] = (p[C.kf83]*u[V.E4]*u[V.ERKstar]/(p[C.Kmf83] + u[V.E4]) - p[C.Vmaxr83]*u[V.E4_PT]/(p[C.Kmr83] + u[V.E4_PT]));
    v[84] = (p[C.kf84]*u[V.E_E1]*u[V.ERKstar]/(p[C.Kmf84] + u[V.E_E1]) - p[C.Vmaxr84]*u[V.E_E1_PT]/(p[C.Kmr84] + u[V.E_E1_PT]));
    v[85] = (p[C.kf85]*u[V.H_E4]*u[V.ERKstar]/(p[C.Kmf85] + u[V.H_E4]) - p[C.Vmaxr85]*u[V.H_E4_PT]/(p[C.Kmr85] + u[V.H_E4_PT]));
    v[86] = (p[C.kon86]*u[V.E]*u[V.E1_PT] - p[C.EGF_off]*u[V.E_E1_PT]);
    v[87] = (p[C.kon87]*u[V.H]*u[V.E4_PT] - p[C.HRGoff_4]*u[V.H_E4_PT]);
    v[88] = (2*p[C.kon88]*u[V.sigmaAP]*u[V.T] - p[C.koff88]*u[V.sigmaAP_T]);
    v[89] = (p[C.kon89]*u[V.ERK]*u[V.MEKstar] - p[C.koff89]*u[V.ERK_MEKstar]);
    v[90] = p[C.kcat90]*u[V.ERK_MEKstar];
    v[91] = (p[C.kon91]*u[V.pERK]*u[V.MEKstar] - p[C.koff91]*u[V.pERK_MEKstar]);
    v[92] = p[C.kcat92]*u[V.pERK_MEKstar];
    v[93] = (p[C.kon93]*u[V.ERKstar]*u[V.ERKPpase] - p[C.koff93]*u[V.ERKstar_ERKPpase]);
    v[94] = p[C.kcat94]*u[V.ERKstar_ERKPpase];
    v[95] = (p[C.kon95]*u[V.pERK]*u[V.ERKPpase] - p[C.koff95]*u[V.pERK_ERKPpase]);
    v[96] = p[C.kcat96]*u[V.pERK_ERKPpase];

    du[V.E1] = -v[1] - v[81];
    du[V.E2] = -v[5] - v[6] - v[8] - v[82];
    du[V.E3] = -v[2];
    du[V.E4] = -v[3] - v[83];
    du[V.E_E1] = v[1] - v[4] - v[4] - v[5] - v[61] - v[62] - v[84];
    du[V.H_E3] = v[2] - v[6] - v[7] - v[61];
    du[V.H_E4] = v[3] - v[7] - v[8] - v[9] - v[9] - v[62] - v[85];
    du[V.E11] = v[4] - v[10];
    du[V.E12] = v[5] - v[11];
    du[V.E23] = v[6] - v[12];
    du[V.E34] = v[7] - v[13];
    du[V.E24] = v[8] - v[14];
    du[V.E44] = v[9] - v[15];
    du[V.E11P] = v[10] - v[16] - v[17] - v[18] - v[73] - v10d;
    du[V.E12P] = v[11] - v[19] - v[20] - v[21] - v[74];
    du[V.E23P] = v[12] - v[22] - v[23] - v[24] - v[25] - v[75];
    du[V.E34P] = v[13] - v[26] - v[27] - v[28] - v[29] - v[76];
    du[V.E24P] = v[14] - v[30] - v[31] - v[32] - v[33] - v[77];
    du[V.E44P] = v[15] - v[34] - v[35] - v[36] - v[37] - v[78];
    du[V.G] = -v[16] - v[19] - v[22] - v[26] - v[30] - v[34] - v[42] - v[57] - v[65] - v[69] + v16d;
    du[V.S] = -v[17] - v[20] - v[23] - v[27] - v[31] - v[35] - v[43] - v[66] - v[70] + v17d;
    du[V.I] = -v[24] - v[28] - v[32] - v[36] - v[44] - v[67] - v[71];
    du[V.R] = -v[18] - v[21] - v[25] - v[29] - v[33] - v[37] - v[45] - v[68] - v[72] + v18d;
    du[V.O] = -v[40] - v[54] - v[58] - v[60];
    du[V.A] = -v[41] - v[46] - v[55] - v[59];
    du[V.E11G] = v[16] - v16d;
    du[V.E11S] = v[17] - v17d;
    du[V.E11R] = v[18] - v18d;
    du[V.E12G] = v[19];
    du[V.E12S] = v[20];
    du[V.E12R] = v[21];
    du[V.E23G] = v[22];
    du[V.E23S] = v[23];
    du[V.E23I] = v[24];
    du[V.E23R] = v[25];
    du[V.E34G] = v[26];
    du[V.E34S] = v[27];
    du[V.E34I] = v[28];
    du[V.E34R] = v[29];
    du[V.E24G] = v[30];
    du[V.E24S] = v[31];
    du[V.E24I] = v[32];
    du[V.E24R] = v[33];
    du[V.E44G] = v[34];
    du[V.E44S] = v[35];
    du[V.E44I] = v[36];
    du[V.E44R] = v[37];
    du[V.sigmaG] = v[16] + v[19] + v[22] + v[26] + v[30] + v[34] - v[40] - v[41] + v[42] + v[65] + v[69] - v16d;
    du[V.sigmaS] = v[17] + v[20] + v[23] + v[27] + v[31] + v[35] - v[38] + v[43] + v[66] + v[70] - v17d;
    du[V.sigmaI] = v[24] + v[28] + v[32] + v[36] + v[44] + v[67] + v[71];
    du[V.sigmaR] = v[18] + v[21] + v[25] + v[29] + v[33] + v[37] + v[45] - v[50] + v[68] + v[72] - v18d;
    du[V.sigmaA] = -v[39] + v[41] + v[46] + v[59];
    du[V.sigmaSP] = v[38] - v[42];
    du[V.sigmaAP] = v[39] - v[43] - v[44] - v[45] - v[88];
    du[V.sigmaG_O] = v[40] - v[59];
    du[V.sigmaG_A] = v[41] - v[60];
    du[V.sigmaSP_G] = v[42];
    du[V.sigmaAP_S] = v[43];
    du[V.sigmaAP_I] = v[44];
    du[V.sigmaAP_R] = v[45];
    du[V.P3_A] = v[46] - v[57];
    du[V.P2] = -v[48];
    du[V.P3] = -v[46] + v[48];
    du[V.Akt] = -v[47];
    du[V.RsD] = -v[49];
    du[V.RsT] = v[49];
    du[V.sigmaRP] = v[50];
    du[V.Raf] = -v[51];
    du[V.Rafstar] = v[51];
    du[V.MEK] = -v[52];
    du[V.MEKstar] = v[52] - v[89] + v[90] - v[91] + v[92];
    du[V.ERK] = -v[89] + v[96];
    du[V.ERKstar] = v[92] - v[93];
    du[V.OP] = v[54];
    du[V.AP] = v[55];
    du[V.A_sigmaG_O] = v[59] + v[60];
    du[V.sigmaA_G] = v[57] - v[58];
    du[V.sigmaA_G_O] = v[58];
    du[V.sigmaO] = v[40] + v[58] + v[60];
    du[V.E13] = v[61] - v[63];
    du[V.E14] = v[62] - v[64];
    du[V.E13P] = v[63] - v[65] - v[66] - v[67] - v[68] - v[79];
    du[V.E14P] = v[64] - v[69] - v[70] - v[71] - v[72] - v[80];
    du[V.E13G] = v[65];
    du[V.E13S] = v[66];
    du[V.E13I] = v[67];
    du[V.E13R] = v[68];
    du[V.E14G] = v[69];
    du[V.E14S] = v[70];
    du[V.E14I] = v[71];
    du[V.E14R] = v[72];
    du[V.T] = -v[73] - v[74] - v[75] - v[76] - v[77] - v[78] - v[79] - v[80] - v[88] + v73d;
    du[V.E11T] = v[73] - v73d;
    du[V.E12T] = v[74];
    du[V.E23T] = v[75];
    du[V.E34T] = v[76];
    du[V.E24T] = v[77];
    du[V.E44T] = v[78];
    du[V.E13T] = v[79];
    du[V.E14T] = v[80];
    du[V.sigmaT] = v[73] + v[74] + v[75] + v[76] + v[77] + v[78] + v[79] + v[80] + v[88] - v73d;
    du[V.E1_PT] = v[81] - v[86];
    du[V.E2_PT] = v[82];
    du[V.E4_PT] = v[83] - v[87];
    du[V.E_E1_PT] = v[84] + v[86];
    du[V.H_E4_PT] = v[85] + v[87];
    du[V.Aktstar] = v[47];
    du[V.sigmaAP_T] = v[88];
    du[V.E] = (-v[1] - v[86])/VeVc;
    du[V.H] = (-v[2] - v[3] - v[87])/VeVc;
    du[V.fint] = p[C.a98]*(-u[V.fint] + p[C.b98]);
    du[V.pERK] = v[90] - v[91] + v[94] - v[95];
    du[V.ERK_MEKstar] = v[89] - v[90];
    du[V.pERK_MEKstar] = v[91] - v[92];
    du[V.pERK_ERKPpase] = v[95] - v[96];
    du[V.ERKPpase] = -v[93] + v[94] - v[95] + v[96];
    du[V.ERKstar_ERKPpase] = v[93] - v[94];

end