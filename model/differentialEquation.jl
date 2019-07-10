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

    EGF_Binding_to_ErbB1 = (p[C.kon1]*u[V.E]*u[V.E1] - p[C.EGF_off]*u[V.E_E1]);
    HRG_Binding_to_ErbB3 = (p[C.kon2]*u[V.H]*u[V.E3] - p[C.HRGoff_3]*u[V.H_E3]);
    HRG_Binding_to_ErbB4 = (p[C.kon3]*u[V.H]*u[V.E4] - p[C.HRGoff_4]*u[V.H_E4]);
    ErbB1_Dimerization = (p[C.kon4]*u[V.E_E1]*u[V.E_E1] - p[C.koff4]*u[V.E11]);
    ErbB1ErbB2_Heterodimerization = (p[C.kon5]*u[V.E_E1]*u[V.E2] - p[C.koff5]*u[V.E12]);
    ErbB2ErbB3_Heterodimerization = (p[C.kon6]*u[V.H_E3]*u[V.E2] - p[C.koff6]*u[V.E23]);
    ErbB3ErbB4_Heterodimerization = (p[C.kon7]*u[V.H_E3]*u[V.H_E4] - p[C.koff7]*u[V.E34]);
    ErbB2ErbB4_Heterodimerization = (p[C.kon8]*u[V.H_E4]*u[V.E2] - p[C.koff8]*u[V.E24]);
    ErbB4_Homodimerization = (p[C.kon9]*u[V.H_E4]*u[V.H_E4] - p[C.koff9]*u[V.E44]);
    E11_Phosphorylation = (p[C.kf10]*u[V.E11] - p[C.VmaxPY]*u[V.E11P]/(p[C.KmPY] + u[V.E11P]) - p[C.kPTP10]*u[V.sigmaT]*u[V.E11P]);
    E12_Phosphorylation = (p[C.kf11]*u[V.E12] - p[C.VmaxPY]*u[V.E12P]/(p[C.KmPY] + u[V.E12P]) - p[C.kPTP11]*u[V.sigmaT]*u[V.E12P]);
    E23_Phosphorylation = (p[C.kf12]*u[V.E23] - p[C.VmaxPY]*u[V.E23P]/(p[C.KmPY] + u[V.E23P]) - p[C.kPTP12]*u[V.sigmaT]*u[V.E23P]);
    E34_Phosphorylation = (p[C.kf13]*u[V.E34] - p[C.VmaxPY]*u[V.E34P]/(p[C.KmPY] + u[V.E34P]) - p[C.kPTP13]*u[V.sigmaT]*u[V.E34P]);
    E24_Phosphorylation = (p[C.kf14]*u[V.E24] - p[C.VmaxPY]*u[V.E24P]/(p[C.KmPY] + u[V.E24P]) - p[C.kPTP14]*u[V.sigmaT]*u[V.E24P]);
    E44_Phosphorylation = (p[C.kf15]*u[V.E44] - p[C.VmaxPY]*u[V.E44P]/(p[C.KmPY] + u[V.E44P]) - p[C.kPTP15]*u[V.sigmaT]*u[V.E44P]);
    E11PGrb2_binding = (4*p[C.kon16]*u[V.E11P]*u[V.G] - p[C.koff16]*fsigmaG*u[V.E11G]);
    E11PShc_binding = (8*p[C.kon17]*u[V.E11P]*u[V.S] - p[C.koff17]*fsigmaS*u[V.E11S]);
    E11PRasGAP_binding = (2*p[C.kon18]*u[V.E11P]*u[V.R] - p[C.koff18]*fsigmaR*u[V.E11R]);
    E12PGrb2_binding = (3*p[C.kon19]*u[V.E12P]*u[V.G] - p[C.koff19]*fsigmaG*u[V.E12G]);
    E12PShc_binding = (6*p[C.kon20]*u[V.E12P]*u[V.S] - p[C.koff20]*fsigmaS*u[V.E12S]);
    E12PRasGAP_binding = (2*p[C.kon21]*u[V.E12P]*u[V.R] - p[C.koff21]*fsigmaR*u[V.E12R]);
    E23PGrb2_binding = (3*p[C.kon22]*u[V.E23P]*u[V.G] - p[C.koff22]*fsigmaG*u[V.E23G]);
    E23PShc_binding = (3*p[C.kon23]*u[V.E23P]*u[V.S] - p[C.koff23]*fsigmaS*u[V.E23S]);
    E23PPI3K_binding = (3*p[C.kon24]*u[V.E23P]*u[V.I] - p[C.koff24]*u[V.E23I]);
    E23PRasGAP_binding = (2*p[C.kon25]*u[V.E23P]*u[V.R] - p[C.koff25]*fsigmaR*u[V.E23R]);
    E34PGrb2_binding = (4*p[C.kon26]*u[V.E34P]*u[V.G] - p[C.koff26]*fsigmaG*u[V.E34G]);
    E34PShc_binding = (3*p[C.kon27]*u[V.E34P]*u[V.S] - p[C.koff27]*fsigmaS*u[V.E34S]);
    E34PPI3K_binding = (4*p[C.kon28]*u[V.E34P]*u[V.I] - p[C.koff28]*u[V.E34I]);
    E34PRasGAP_binding = (2*p[C.kon29]*u[V.E34P]*u[V.R] - p[C.koff29]*fsigmaR*u[V.E34R]);
    E24PGrb2_binding = (3*p[C.kon30]*u[V.E24P]*u[V.G] - p[C.koff30]*fsigmaG*u[V.E24G]);
    E24PShc_binding = (4*p[C.kon31]*u[V.E24P]*u[V.S] - p[C.koff31]*fsigmaS*u[V.E24S]);
    E24PPI3K_binding = (1*p[C.kon32]*u[V.E24P]*u[V.I] - p[C.koff32]*u[V.E24I]);
    E24PRasGAP_binding = (2*p[C.kon33]*u[V.E24P]*u[V.R] - p[C.koff33]*fsigmaR*u[V.E24R]);
    E44PRasGAP_binding_1 = (4*p[C.kon34]*u[V.E44P]*u[V.G] - p[C.koff34]*fsigmaG*u[V.E44G]);
    E44PShc_binding = (4*p[C.kon35]*u[V.E44P]*u[V.S] - p[C.koff35]*fsigmaS*u[V.E44S]);
    E44PPI3K_binding = (2*p[C.kon36]*u[V.E44P]*u[V.I] - p[C.koff36]*u[V.E44I]);
    E44PRasGAP_binding_2 = (2*p[C.kon37]*u[V.E44P]*u[V.R] - p[C.koff37]*fsigmaR*u[V.E44R]);
    Shc_Phosphorylation = (p[C.kf38]*u[V.sigmaS]*sigmaEP - p[C.VmaxPY]*u[V.sigmaSP]/(p[C.KmPY] + u[V.sigmaSP]) - p[C.kPTP38]*u[V.sigmaT]*u[V.sigmaSP]);
    Gab1_Phosphorylation = (p[C.kf39]*u[V.sigmaA]*sigmaEP - p[C.VmaxPY]*u[V.sigmaAP]/(p[C.KmPY] + u[V.sigmaAP]) - p[C.kPTP39]*u[V.sigmaT]*u[V.sigmaAP]);
    Grb2SOS_binding = (p[C.kon40]*u[V.sigmaG]*u[V.O] - p[C.koff40]*u[V.sigmaG_O]);
    Grb2Gab1_binding = (p[C.kon41]*u[V.sigmaG]*u[V.A] - p[C.koff41]*u[V.sigmaG_A]*fsigmaA);
    ShcGrb2_binding = (p[C.kon42]*u[V.sigmaSP]*u[V.G] - p[C.koff42]*u[V.sigmaSP_G]*fsigmaG);
    Gab1Shc_binding = (3*p[C.kon43]*u[V.sigmaAP]*u[V.S] - p[C.koff43]*u[V.sigmaAP_S]*fsigmaS);
    Gab1PI3K_binding = (3*p[C.kon44]*u[V.sigmaAP]*u[V.I] - p[C.koff44]*u[V.sigmaAP_I]);
    Gab1RasGAP_binding = (2*p[C.kon45]*u[V.sigmaAP]*u[V.R] - p[C.koff45]*u[V.sigmaAP_R]*fsigmaR);
    Gab1PIP3_binding = (p[C.kon46]*u[V.P3]*u[V.A] - p[C.koff46]*u[V.P3_A]*fsigmaA);
    Akt_Activation = (p[C.kf47]*u[V.P3]*u[V.Akt]/(p[C.Kmf47] + u[V.Akt]) - p[C.Vmaxr47]*u[V.Aktstar]/(p[C.Kmr47] + u[V.Aktstar]));
    PIP3_Production = (p[C.kf48]*(1 - u[V.fint]*f11)*u[V.sigmaI]*u[V.P2]/(p[C.Kmf48] + u[V.P2]) - 3*p[C.PTEN]*u[V.P3]/(p[C.Kmr48] + u[V.P3]));
    RasGTP_Production = (p[C.kf49]*u[V.sigmaO]*u[V.RsD]/(p[C.Kmf49] + u[V.RsD]) - p[C.kr49]*u[V.sigmaR]*u[V.RsT]/(p[C.Kmr49] + u[V.RsT]) - p[C.kr49b]*u[V.sigmaRP]*u[V.RsT]/(p[C.Kmr49b] + u[V.RsT]) - p[C.kcon49]*u[V.RsT]);
    RasGAP_Phosphorylation = (p[C.kf50]*u[V.sigmaR]*sigmaEP - p[C.VmaxPY]*u[V.sigmaRP]/(p[C.KmPY] + u[V.sigmaRP]) - p[C.kPTP50]*u[V.sigmaT]*u[V.sigmaRP]);
    Raf_Activation = (p[C.kf51]*u[V.RsT]*u[V.Raf]/(p[C.Kmf51] + u[V.Raf]) - p[C.Vmaxr51]*u[V.Rafstar]/(p[C.Kmrb51] + u[V.Rafstar]));
    MEK_Activation = (p[C.kf52]*u[V.Rafstar]*u[V.MEK]/(p[C.Kmf52] + u[V.MEK]) - p[C.Vmaxr52]*u[V.MEKstar]/(p[C.Kmr52] + u[V.MEKstar]));

    SOS_ST_Phosphorylation = (p[C.kf54]*u[V.O]*u[V.ERKstar]/(p[C.Kmf54] + u[V.O]) - p[C.Vmaxr54]*u[V.OP]/(p[C.Kmr54] + u[V.OP]));
    Gab1_ST_Phosphorylation = (p[C.kf55]*u[V.A]*u[V.ERKstar]/(p[C.Kmf55] + u[V.A]) - p[C.Vmaxr55]*u[V.AP]/(p[C.Kmr55] + u[V.AP]));

    Grb2_binding_to_P3_A = (p[C.kon57]*u[V.P3_A]*u[V.G] - p[C.koff57]*u[V.sigmaA_G]);
    SOS_binding_to_SigAG = (p[C.kon58]*u[V.sigmaA_G]*u[V.O] - p[C.koff58]*u[V.sigmaA_G_O]);
    Gab1_binding_to_SigGO = (p[C.kon59]*u[V.sigmaG_O]*u[V.A] - p[C.koff59]*u[V.A_sigmaG_O]*fsigmaA);
    SOS_binding_to_SigGA = (p[C.kon60]*u[V.sigmaG_A]*u[V.O] - p[C.koff60]*u[V.A_sigmaG_O]);
    ErbB1ErbB3_Heterodimerization = (p[C.kon61]*u[V.H_E3]*u[V.E_E1] - p[C.koff61]*u[V.E13]);
    ErbB1ErbB4_Heterodimerization = (p[C.kon62]*u[V.H_E4]*u[V.E_E1] - p[C.koff62]*u[V.E14]);
    E13_Phosphorylation = (p[C.kf63]*u[V.E13] - p[C.VmaxPY]*u[V.E13P]/(p[C.KmPY] + u[V.E13P]) - p[C.kPTP63]*u[V.sigmaT]*u[V.E13P]);
    E14_Phosphorylation = (p[C.kf64]*u[V.E14] - p[C.VmaxPY]*u[V.E14P]/(p[C.KmPY] + u[V.E14P]) - p[C.kPTP64]*u[V.sigmaT]*u[V.E14P]);
    E13Grb2_Binding = (4*p[C.kon65]*u[V.E13P]*u[V.G] - p[C.koff65]*fsigmaG*u[V.E13G]);
    E13Shc_Binding = (5*p[C.kon66]*u[V.E13P]*u[V.S] - p[C.koff66]*fsigmaS*u[V.E13S]);
    E13PI3K_Binding = (3*p[C.kon67]*u[V.E13P]*u[V.I] - p[C.koff67]*u[V.E13I]);
    E13RasGAP_Binding = (2*p[C.kon68]*u[V.E13P]*u[V.R] - p[C.koff68]*fsigmaR*u[V.E13R]);
    E14Grb2_Binding = (4*p[C.kon69]*u[V.E14P]*u[V.G] - p[C.koff69]*fsigmaG*u[V.E14G]);
    E14Shc_Binding = (6*p[C.kon70]*u[V.E14P]*u[V.S] - p[C.koff70]*fsigmaS*u[V.E14S]);
    E14PI3K_Binding = (1*p[C.kon71]*u[V.E14P]*u[V.I] - p[C.koff71]*u[V.E14I]);
    E14RasGAP_Binding = (2*p[C.kon72]*u[V.E14P]*u[V.R] - p[C.koff72]*fsigmaR*u[V.E14R]);
    E11PTP1B_Binding = (4*p[C.kon73]*u[V.E11P]*u[V.T] - p[C.koff73]*u[V.E11T]);
    E12PTP1B_Binding = (3*p[C.kon74]*u[V.E12P]*u[V.T] - p[C.koff74]*u[V.E12T]);
    E23PTP1B_Binding = (2*p[C.kon75]*u[V.E23P]*u[V.T] - p[C.koff75]*u[V.E23T]);
    E34PTP1B_Binding = (2*p[C.kon76]*u[V.E34P]*u[V.T] - p[C.koff76]*u[V.E34T]);
    E24PTP1B_Binding = (2*p[C.kon77]*u[V.E24P]*u[V.T] - p[C.koff77]*u[V.E24T]);
    E44PTP1B_Binding = (2*p[C.kon78]*u[V.E44P]*u[V.T] - p[C.koff78]*u[V.E44T]);
    E13PTP1B_Binding = (3*p[C.kon79]*u[V.E13P]*u[V.T] - p[C.koff79]*u[V.E13T]);
    E14PTP1B_Binding = (3*p[C.kon80]*u[V.E14P]*u[V.T] - p[C.koff80]*u[V.E14T]);
    E1_ST_Phosphorylation = (p[C.kf81]*u[V.E1]*u[V.ERKstar]/(p[C.Kmf81] + u[V.E1]) - p[C.Vmaxr81]*u[V.E1_PT]/(p[C.Kmr81] + u[V.E1_PT]));
    E2_ST_Phosphorylation = (p[C.kf82]*u[V.E2]*u[V.ERKstar]/(p[C.Kmf82] + u[V.E2]) - p[C.Vmaxr82]*u[V.E2_PT]/(p[C.Kmr82] + u[V.E2_PT]));
    E4_ST_Phosphorylation = (p[C.kf83]*u[V.E4]*u[V.ERKstar]/(p[C.Kmf83] + u[V.E4]) - p[C.Vmaxr83]*u[V.E4_PT]/(p[C.Kmr83] + u[V.E4_PT]));
    E_E1_ST_Phosphorylation = (p[C.kf84]*u[V.E_E1]*u[V.ERKstar]/(p[C.Kmf84] + u[V.E_E1]) - p[C.Vmaxr84]*u[V.E_E1_PT]/(p[C.Kmr84] + u[V.E_E1_PT]));
    H_E4_ST_Phosphorylation = (p[C.kf85]*u[V.H_E4]*u[V.ERKstar]/(p[C.Kmf85] + u[V.H_E4]) - p[C.Vmaxr85]*u[V.H_E4_PT]/(p[C.Kmr85] + u[V.H_E4_PT]));
    EGF_binding_to_E1PT = (p[C.kon86]*u[V.E]*u[V.E1_PT] - p[C.EGF_off]*u[V.E_E1_PT]);
    HRG_binding_to_E4PT = (p[C.kon87]*u[V.H]*u[V.E4_PT] - p[C.HRGoff_4]*u[V.H_E4_PT]);
    PTP1B_binding_to_SigAP = (2*p[C.kon88]*u[V.sigmaAP]*u[V.T] - p[C.koff88]*u[V.sigmaAP_T]);
    E11P_Degradation = p[C.kdeg]*u[V.E11P];
    E11G_Degradation = p[C.kdeg]*u[V.E11G];
    E11S_Degradation = p[C.kdeg]*u[V.E11S];
    E11R_Degradation = p[C.kdeg]*u[V.E11R];
    E11T_Degradation = p[C.kdeg]*u[V.E11T];
    PI3K_Fractional_Multiplier = p[C.a98]*(-u[V.fint] + p[C.b98]);
    EGF_in_EC_compartment = (p[C.kon1]*u[V.E]*u[V.E1] - p[C.EGF_off]*u[V.E_E1] + p[C.kon86]*u[V.E]*u[V.E1_PT] - p[C.EGF_off]*u[V.E_E1_PT])/VeVc;
    HRG_in_EC_compartment = (p[C.kon87]*u[V.H]*u[V.E4_PT] - p[C.HRGoff_4]*u[V.H_E4_PT] + p[C.kon2]*u[V.H]*u[V.E3] - p[C.HRGoff_3]*u[V.H_E3] + p[C.kon3]*u[V.H]*u[V.E4] - p[C.HRGoff_4]*u[V.H_E4])/VeVc;
    ERK_binding_to_MEKstar_1 = (p[C.kon89]*u[V.ERK]*u[V.MEKstar] - p[C.koff89]*u[V.ERK_MEKstar]);
    pERK_production = p[C.kcat90]*u[V.ERK_MEKstar];
    ERK_binding_to_MEKstar_2 = (p[C.kon91]*u[V.pERK]*u[V.MEKstar] - p[C.koff91]*u[V.pERK_MEKstar]);
    ERKstar_production = p[C.kcat92]*u[V.pERK_MEKstar];
    ERKstar_binding_to_Phosphatase = (p[C.kon93]*u[V.ERKstar]*u[V.ERKPpase] - p[C.koff93]*u[V.ERKstar_ERKPpase]);
    ERKstar_dephosphorylation = p[C.kcat94]*u[V.ERKstar_ERKPpase];
    pERK_binding_to_Phosphatase = (p[C.kon95]*u[V.pERK]*u[V.ERKPpase] - p[C.koff95]*u[V.pERK_ERKPpase]);
    pERK_dephosphorylation = p[C.kcat96]*u[V.pERK_ERKPpase];

    du[V.E1] = ((-1.0 * EGF_Binding_to_ErbB1) + (-1.0 * E1_ST_Phosphorylation));
    du[V.E2] = ((-1.0 * ErbB1ErbB2_Heterodimerization) + (-1.0 * ErbB2ErbB3_Heterodimerization) + (-1.0 * ErbB2ErbB4_Heterodimerization) + (-1.0 * E2_ST_Phosphorylation));
    du[V.E3] = ((-1.0 * HRG_Binding_to_ErbB3));
    du[V.E4] = ((-1.0 * HRG_Binding_to_ErbB4) + (-1.0 * E4_ST_Phosphorylation));
    du[V.E_E1] = (( 1.0 * EGF_Binding_to_ErbB1) + (-1.0 * ErbB1_Dimerization) + (-1.0 * ErbB1_Dimerization) + (-1.0 * ErbB1ErbB2_Heterodimerization) + (-1.0 * ErbB1ErbB3_Heterodimerization) + (-1.0 * ErbB1ErbB4_Heterodimerization) + (-1.0 * E_E1_ST_Phosphorylation));
    du[V.H_E3] = (( 1.0 * HRG_Binding_to_ErbB3) + (-1.0 * ErbB2ErbB3_Heterodimerization) + (-1.0 * ErbB3ErbB4_Heterodimerization) + (-1.0 * ErbB1ErbB3_Heterodimerization));
    du[V.H_E4] = (( 1.0 * HRG_Binding_to_ErbB4) + (-1.0 * ErbB3ErbB4_Heterodimerization) + (-1.0 * ErbB2ErbB4_Heterodimerization) + (-1.0 * ErbB4_Homodimerization) + (-1.0 * ErbB4_Homodimerization) + (-1.0 * ErbB1ErbB4_Heterodimerization) + (-1.0 * H_E4_ST_Phosphorylation));
    du[V.E11] = (( 1.0 * ErbB1_Dimerization) + (-1.0 * E11_Phosphorylation));
    du[V.E12] = (( 1.0 * ErbB1ErbB2_Heterodimerization) + (-1.0 * E12_Phosphorylation));
    du[V.E23] = (( 1.0 * ErbB2ErbB3_Heterodimerization) + (-1.0 * E23_Phosphorylation));
    du[V.E34] = (( 1.0 * ErbB3ErbB4_Heterodimerization) + (-1.0 * E34_Phosphorylation));
    du[V.E24] = (( 1.0 * ErbB2ErbB4_Heterodimerization) + (-1.0 * E24_Phosphorylation));
    du[V.E44] = (( 1.0 * ErbB4_Homodimerization) + (-1.0 * E44_Phosphorylation));
    du[V.E11P] = (( 1.0 * E11_Phosphorylation) + (-1.0 * E11PGrb2_binding) + (-1.0 * E11PShc_binding) + (-1.0 * E11PRasGAP_binding) + (-1.0 * E11PTP1B_Binding) + (-1.0 * E11P_Degradation));
    du[V.E12P] = (( 1.0 * E12_Phosphorylation) + (-1.0 * E12PGrb2_binding) + (-1.0 * E12PShc_binding) + (-1.0 * E12PRasGAP_binding) + (-1.0 * E12PTP1B_Binding));
    du[V.E23P] = (( 1.0 * E23_Phosphorylation) + (-1.0 * E23PGrb2_binding) + (-1.0 * E23PShc_binding) + (-1.0 * E23PPI3K_binding) + (-1.0 * E23PRasGAP_binding) + (-1.0 * E23PTP1B_Binding));
    du[V.E34P] = (( 1.0 * E34_Phosphorylation) + (-1.0 * E34PGrb2_binding) + (-1.0 * E34PShc_binding) + (-1.0 * E34PPI3K_binding) + (-1.0 * E34PRasGAP_binding) + (-1.0 * E34PTP1B_Binding));
    du[V.E24P] = (( 1.0 * E24_Phosphorylation) + (-1.0 * E24PGrb2_binding) + (-1.0 * E24PShc_binding) + (-1.0 * E24PPI3K_binding) + (-1.0 * E24PRasGAP_binding) + (-1.0 * E24PTP1B_Binding));
    du[V.E44P] = (( 1.0 * E44_Phosphorylation) + (-1.0 * E44PRasGAP_binding_1) + (-1.0 * E44PShc_binding) + (-1.0 * E44PPI3K_binding) + (-1.0 * E44PRasGAP_binding_2) + (-1.0 * E44PTP1B_Binding));
    du[V.G] = ((-1.0 * E11PGrb2_binding) + (-1.0 * E12PGrb2_binding) + (-1.0 * E23PGrb2_binding) + (-1.0 * E34PGrb2_binding) + (-1.0 * E24PGrb2_binding) + (-1.0 * E44PRasGAP_binding_1) + (-1.0 * ShcGrb2_binding) + (-1.0 * Grb2_binding_to_P3_A) + (-1.0 * E13Grb2_Binding) + (-1.0 * E14Grb2_Binding) + ( 1.0 * E11G_Degradation));
    du[V.S] = ((-1.0 * E11PShc_binding) + (-1.0 * E12PShc_binding) + (-1.0 * E23PShc_binding) + (-1.0 * E34PShc_binding) + (-1.0 * E24PShc_binding) + (-1.0 * E44PShc_binding) + (-1.0 * Gab1Shc_binding) + (-1.0 * E13Shc_Binding) + (-1.0 * E14Shc_Binding) + ( 1.0 * E11S_Degradation));
    du[V.I] = ((-1.0 * E23PPI3K_binding) + (-1.0 * E34PPI3K_binding) + (-1.0 * E24PPI3K_binding) + (-1.0 * E44PPI3K_binding) + (-1.0 * Gab1PI3K_binding) + (-1.0 * E13PI3K_Binding) + (-1.0 * E14PI3K_Binding));
    du[V.R] = ((-1.0 * E11PRasGAP_binding) + (-1.0 * E12PRasGAP_binding) + (-1.0 * E23PRasGAP_binding) + (-1.0 * E34PRasGAP_binding) + (-1.0 * E24PRasGAP_binding) + (-1.0 * E44PRasGAP_binding_2) + (-1.0 * Gab1RasGAP_binding) + (-1.0 * E13RasGAP_Binding) + (-1.0 * E14RasGAP_Binding) + ( 1.0 * E11R_Degradation));
    du[V.O] = ((-1.0 * Grb2SOS_binding) + (-1.0 * SOS_ST_Phosphorylation) + (-1.0 * SOS_binding_to_SigAG) + (-1.0 * SOS_binding_to_SigGA));
    du[V.A] = ((-1.0 * Grb2Gab1_binding) + (-1.0 * Gab1PIP3_binding) + (-1.0 * Gab1_ST_Phosphorylation) + (-1.0 * Gab1_binding_to_SigGO));
    du[V.E11G] = (( 1.0 * E11PGrb2_binding) + (-1.0 * E11G_Degradation));
    du[V.E11S] = (( 1.0 * E11PShc_binding) + (-1.0 * E11S_Degradation));
    du[V.E11R] = (( 1.0 * E11PRasGAP_binding) + (-1.0 * E11R_Degradation));
    du[V.E12G] = (( 1.0 * E12PGrb2_binding));
    du[V.E12S] = (( 1.0 * E12PShc_binding));
    du[V.E12R] = (( 1.0 * E12PRasGAP_binding));
    du[V.E23G] = (( 1.0 * E23PGrb2_binding));
    du[V.E23S] = (( 1.0 * E23PShc_binding));
    du[V.E23I] = (( 1.0 * E23PPI3K_binding));
    du[V.E23R] = (( 1.0 * E23PRasGAP_binding));
    du[V.E34G] = (( 1.0 * E34PGrb2_binding));
    du[V.E34S] = (( 1.0 * E34PShc_binding));
    du[V.E34I] = (( 1.0 * E34PPI3K_binding));
    du[V.E34R] = (( 1.0 * E34PRasGAP_binding));
    du[V.E24G] = (( 1.0 * E24PGrb2_binding));
    du[V.E24S] = (( 1.0 * E24PShc_binding));
    du[V.E24I] = (( 1.0 * E24PPI3K_binding));
    du[V.E24R] = (( 1.0 * E24PRasGAP_binding));
    du[V.E44G] = (( 1.0 * E44PRasGAP_binding_1));
    du[V.E44S] = (( 1.0 * E44PShc_binding));
    du[V.E44I] = (( 1.0 * E44PPI3K_binding));
    du[V.E44R] = (( 1.0 * E44PRasGAP_binding_2));
    du[V.sigmaG] = (( 1.0 * E11PGrb2_binding) + ( 1.0 * E12PGrb2_binding) + ( 1.0 * E23PGrb2_binding) + ( 1.0 * E34PGrb2_binding) + ( 1.0 * E24PGrb2_binding) + ( 1.0 * E44PRasGAP_binding_1) + (-1.0 * Grb2SOS_binding) + (-1.0 * Grb2Gab1_binding) + ( 1.0 * ShcGrb2_binding) + ( 1.0 * E13Grb2_Binding) + ( 1.0 * E14Grb2_Binding) + (-1.0 * E11G_Degradation));
    du[V.sigmaS] = (( 1.0 * E11PShc_binding) + ( 1.0 * E12PShc_binding) + ( 1.0 * E23PShc_binding) + ( 1.0 * E34PShc_binding) + ( 1.0 * E24PShc_binding) + ( 1.0 * E44PShc_binding) + (-1.0 * Shc_Phosphorylation) + ( 1.0 * Gab1Shc_binding) + ( 1.0 * E13Shc_Binding) + ( 1.0 * E14Shc_Binding) + (-1.0 * E11S_Degradation));
    du[V.sigmaI] = (( 1.0 * E23PPI3K_binding) + ( 1.0 * E34PPI3K_binding) + ( 1.0 * E24PPI3K_binding) + ( 1.0 * E44PPI3K_binding) + ( 1.0 * Gab1PI3K_binding) + ( 1.0 * E13PI3K_Binding) + ( 1.0 * E14PI3K_Binding));
    du[V.sigmaR] = (( 1.0 * E11PRasGAP_binding) + ( 1.0 * E12PRasGAP_binding) + ( 1.0 * E23PRasGAP_binding) + ( 1.0 * E34PRasGAP_binding) + ( 1.0 * E24PRasGAP_binding) + ( 1.0 * E44PRasGAP_binding_2) + ( 1.0 * Gab1RasGAP_binding) + (-1.0 * RasGAP_Phosphorylation) + ( 1.0 * E13RasGAP_Binding) + ( 1.0 * E14RasGAP_Binding) + (-1.0 * E11R_Degradation));
    du[V.sigmaA] = ((-1.0 * Gab1_Phosphorylation) + ( 1.0 * Grb2Gab1_binding) + ( 1.0 * Gab1PIP3_binding) + ( 1.0 * Gab1_binding_to_SigGO));
    du[V.sigmaSP] = (( 1.0 * Shc_Phosphorylation) + (-1.0 * ShcGrb2_binding));
    du[V.sigmaAP] = (( 1.0 * Gab1_Phosphorylation) + (-1.0 * Gab1Shc_binding) + (-1.0 * Gab1PI3K_binding) + (-1.0 * Gab1RasGAP_binding) + (-1.0 * PTP1B_binding_to_SigAP));
    du[V.sigmaG_O] = (( 1.0 * Grb2SOS_binding) + (-1.0 * Gab1_binding_to_SigGO));
    du[V.sigmaG_A] = (( 1.0 * Grb2Gab1_binding) + (-1.0 * SOS_binding_to_SigGA));
    du[V.sigmaSP_G] = (( 1.0 * ShcGrb2_binding));
    du[V.sigmaAP_S] = (( 1.0 * Gab1Shc_binding));
    du[V.sigmaAP_I] = (( 1.0 * Gab1PI3K_binding));
    du[V.sigmaAP_R] = (( 1.0 * Gab1RasGAP_binding));
    du[V.Empty] = (( 1.0 * E11P_Degradation) + (-1.0 * PI3K_Fractional_Multiplier) + ( 1.0 * EGF_in_EC_compartment) + ( 1.0 * HRG_in_EC_compartment));
    du[V.P3_A] = (( 1.0 * Gab1PIP3_binding) + (-1.0 * Grb2_binding_to_P3_A));
    du[V.P2] = ((-1.0 * PIP3_Production));
    du[V.P3] = ((-1.0 * Gab1PIP3_binding) + ( 1.0 * PIP3_Production));
    du[V.Akt] = ((-1.0 * Akt_Activation));
    du[V.RsD] = ((-1.0 * RasGTP_Production));
    du[V.RsT] = (( 1.0 * RasGTP_Production));
    du[V.sigmaRP] = (( 1.0 * RasGAP_Phosphorylation));
    du[V.Raf] = ((-1.0 * Raf_Activation));
    du[V.Rafstar] = (( 1.0 * Raf_Activation));
    du[V.MEK] = ((-1.0 * MEK_Activation));
    du[V.MEKstar] = (( 1.0 * MEK_Activation) + (-1.0 * ERK_binding_to_MEKstar_1) + ( 1.0 * pERK_production) + (-1.0 * ERK_binding_to_MEKstar_2) + ( 1.0 * ERKstar_production));
    du[V.ERK] = ((-1.0 * ERK_binding_to_MEKstar_1) + ( 1.0 * pERK_dephosphorylation));
    du[V.ERKstar] = (( 1.0 * ERKstar_production) + (-1.0 * ERKstar_binding_to_Phosphatase));
    du[V.OP] = (( 1.0 * SOS_ST_Phosphorylation));
    du[V.AP] = (( 1.0 * Gab1_ST_Phosphorylation));
    du[V.A_sigmaG_O] = (( 1.0 * Gab1_binding_to_SigGO) + ( 1.0 * SOS_binding_to_SigGA));
    du[V.sigmaA_G] = (( 1.0 * Grb2_binding_to_P3_A) + (-1.0 * SOS_binding_to_SigAG));
    du[V.sigmaA_G_O] = (( 1.0 * SOS_binding_to_SigAG));
    du[V.sigmaO] = (( 1.0 * Grb2SOS_binding) + ( 1.0 * SOS_binding_to_SigAG) + ( 1.0 * SOS_binding_to_SigGA));
    du[V.E13] = (( 1.0 * ErbB1ErbB3_Heterodimerization) + (-1.0 * E13_Phosphorylation));
    du[V.E14] = (( 1.0 * ErbB1ErbB4_Heterodimerization) + (-1.0 * E14_Phosphorylation));
    du[V.E13P] = (( 1.0 * E13_Phosphorylation) + (-1.0 * E13Grb2_Binding) + (-1.0 * E13Shc_Binding) + (-1.0 * E13PI3K_Binding) + (-1.0 * E13RasGAP_Binding) + (-1.0 * E13PTP1B_Binding));
    du[V.E14P] = (( 1.0 * E14_Phosphorylation) + (-1.0 * E14Grb2_Binding) + (-1.0 * E14Shc_Binding) + (-1.0 * E14PI3K_Binding) + (-1.0 * E14RasGAP_Binding) + (-1.0 * E14PTP1B_Binding));
    du[V.E13G] = (( 1.0 * E13Grb2_Binding));
    du[V.E13S] = (( 1.0 * E13Shc_Binding));
    du[V.E13I] = (( 1.0 * E13PI3K_Binding));
    du[V.E13R] = (( 1.0 * E13RasGAP_Binding));
    du[V.E14G] = (( 1.0 * E14Grb2_Binding));
    du[V.E14S] = (( 1.0 * E14Shc_Binding));
    du[V.E14I] = (( 1.0 * E14PI3K_Binding));
    du[V.E14R] = (( 1.0 * E14RasGAP_Binding));
    du[V.T] = ((-1.0 * E11PTP1B_Binding) + (-1.0 * E12PTP1B_Binding) + (-1.0 * E23PTP1B_Binding) + (-1.0 * E34PTP1B_Binding) + (-1.0 * E24PTP1B_Binding) + (-1.0 * E44PTP1B_Binding) + (-1.0 * E13PTP1B_Binding) + (-1.0 * E14PTP1B_Binding) + (-1.0 * PTP1B_binding_to_SigAP) + ( 1.0 * E11T_Degradation));
    du[V.E11T] = (( 1.0 * E11PTP1B_Binding) + (-1.0 * E11T_Degradation));
    du[V.E12T] = (( 1.0 * E12PTP1B_Binding));
    du[V.E23T] = (( 1.0 * E23PTP1B_Binding));
    du[V.E34T] = (( 1.0 * E34PTP1B_Binding));
    du[V.E24T] = (( 1.0 * E24PTP1B_Binding));
    du[V.E44T] = (( 1.0 * E44PTP1B_Binding));
    du[V.E13T] = (( 1.0 * E13PTP1B_Binding));
    du[V.E14T] = (( 1.0 * E14PTP1B_Binding));
    du[V.sigmaT] = (( 1.0 * E11PTP1B_Binding) + ( 1.0 * E12PTP1B_Binding) + ( 1.0 * E23PTP1B_Binding) + ( 1.0 * E34PTP1B_Binding) + ( 1.0 * E24PTP1B_Binding) + ( 1.0 * E44PTP1B_Binding) + ( 1.0 * E13PTP1B_Binding) + ( 1.0 * E14PTP1B_Binding) + ( 1.0 * PTP1B_binding_to_SigAP) + (-1.0 * E11T_Degradation));
    du[V.E1_PT] = (( 1.0 * E1_ST_Phosphorylation) + (-1.0 * EGF_binding_to_E1PT));
    du[V.E2_PT] = (( 1.0 * E2_ST_Phosphorylation));
    du[V.E4_PT] = (( 1.0 * E4_ST_Phosphorylation) + (-1.0 * HRG_binding_to_E4PT));
    du[V.E_E1_PT] = (( 1.0 * E_E1_ST_Phosphorylation) + ( 1.0 * EGF_binding_to_E1PT));
    du[V.H_E4_PT] = (( 1.0 * H_E4_ST_Phosphorylation) + ( 1.0 * HRG_binding_to_E4PT));
    du[V.Aktstar] = (( 1.0 * Akt_Activation));
    du[V.sigmaAP_T] = (( 1.0 * PTP1B_binding_to_SigAP));
    du[V.E] = ((-1.0 * EGF_in_EC_compartment));
    du[V.H] = ((-1.0 * HRG_in_EC_compartment));
    du[V.fint] = (( 1.0 * PI3K_Fractional_Multiplier));
    du[V.pERK] = (( 1.0 * pERK_production) + (-1.0 * ERK_binding_to_MEKstar_2) + ( 1.0 * ERKstar_dephosphorylation) + (-1.0 * pERK_binding_to_Phosphatase));
    du[V.ERK_MEKstar] = (( 1.0 * ERK_binding_to_MEKstar_1) + (-1.0 * pERK_production));
    du[V.pERK_MEKstar] = (( 1.0 * ERK_binding_to_MEKstar_2) + (-1.0 * ERKstar_production));
    du[V.pERK_ERKPpase] = (( 1.0 * pERK_binding_to_Phosphatase) + (-1.0 * pERK_dephosphorylation));
    du[V.ERKPpase] = ((-1.0 * ERKstar_binding_to_Phosphatase) + ( 1.0 * ERKstar_dephosphorylation) + (-1.0 * pERK_binding_to_Phosphatase) + ( 1.0 * pERK_dephosphorylation));
    du[V.ERKstar_ERKPpase] = (( 1.0 * ERKstar_binding_to_Phosphatase) + (-1.0 * ERKstar_dephosphorylation));

end