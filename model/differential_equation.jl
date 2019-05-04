function diffeq(t,u)

    # extracellular volume to cytoplasmic volume ratio
    VeVc = 33.3;

    # fraction definitions
    if u[sigmaS] + u[sigmaSP] + u[sigmaSP_G] > 0.0
        fsigmaS = u[sigmaS]/(u[sigmaS] + u[sigmaSP] + u[sigmaSP_G]);
    else
        fsigmaS = 0.0;
    end

    if u[sigmaG] + u[sigmaG_A] + u[sigmaG_O] + u[A_sigmaG_O] > 0.0
        fsigmaG = u[sigmaG]/(u[sigmaG] + u[sigmaG_A] + u[sigmaG_O] + u[A_sigmaG_O]);
    else
        fsigmaG = 0.0;
    end

    if u[sigmaA] + u[sigmaAP] + u[sigmaAP_S] + u[sigmaAP_R] + u[sigmaAP_I] + u[sigmaAP_T] > 0.0
        fsigmaA = u[sigmaA]/(u[sigmaA] + u[sigmaAP] + u[sigmaAP_S] + u[sigmaAP_R] + u[sigmaAP_I] + u[sigmaAP_T]);
    else
        fsigmaA = 0.0;
    end

    if u[sigmaR] + u[sigmaRP] > 0.0
        fsigmaR = u[sigmaR]/(u[sigmaR] + u[sigmaRP]);
    else
        fsigmaR = 0.0;
    end

    sigmaEP = u[E11P] + u[E12P] + u[E23P] + u[E24P] + u[E34P] + u[E44P] + u[E13P] + u[E14P];
    if sigmaEP > 0.0
        f11 = u[E11P]/sigmaEP;
    else
        f11 = 0.0;
    end

    EGF_Binding_to_ErbB1 = (p[kon1]*u[E]*u[E1] - p[EGF_off]*u[E_E1]);
    HRG_Binding_to_ErbB3 = (p[kon2]*u[H]*u[E3] - p[HRGoff_3]*u[H_E3]);
    HRG_Binding_to_ErbB4 = (p[kon3]*u[H]*u[E4] - p[HRGoff_4]*u[H_E4]);
    ErbB1_Dimerization = (p[kon4]*u[E_E1]*u[E_E1] - p[koff4]*u[E11]);
    ErbB1ErbB2_Heterodimerization = (p[kon5]*u[E_E1]*u[E2] - p[koff5]*u[E12]);
    ErbB2ErbB3_Heterodimerization = (p[kon6]*u[H_E3]*u[E2] - p[koff6]*u[E23]);
    ErbB3ErbB4_Heterodimerization = (p[kon7]*u[H_E3]*u[H_E4] - p[koff7]*u[E34]);
    ErbB2ErbB4_Heterodimerization = (p[kon8]*u[H_E4]*u[E2] - p[koff8]*u[E24]);
    ErbB4_Homodimerization = (p[kon9]*u[H_E4]*u[H_E4] - p[koff9]*u[E44]);
    E11_Phosphorylation = (p[kf10]*u[E11] - p[VmaxPY]*u[E11P]/(p[KmPY] + u[E11P]) - p[kPTP10]*u[sigmaT]*u[E11P]);
    E12_Phosphorylation = (p[kf11]*u[E12] - p[VmaxPY]*u[E12P]/(p[KmPY] + u[E12P]) - p[kPTP11]*u[sigmaT]*u[E12P]);
    E23_Phosphorylation = (p[kf12]*u[E23] - p[VmaxPY]*u[E23P]/(p[KmPY] + u[E23P]) - p[kPTP12]*u[sigmaT]*u[E23P]);
    E34_Phosphorylation = (p[kf13]*u[E34] - p[VmaxPY]*u[E34P]/(p[KmPY] + u[E34P]) - p[kPTP13]*u[sigmaT]*u[E34P]);
    E24_Phosphorylation = (p[kf14]*u[E24] - p[VmaxPY]*u[E24P]/(p[KmPY] + u[E24P]) - p[kPTP14]*u[sigmaT]*u[E24P]);
    E44_Phosphorylation = (p[kf15]*u[E44] - p[VmaxPY]*u[E44P]/(p[KmPY] + u[E44P]) - p[kPTP15]*u[sigmaT]*u[E44P]);
    E11PGrb2_binding = (4*p[kon16]*u[E11P]*u[G] - p[koff16]*fsigmaG*u[E11G]);
    E11PShc_binding = (8*p[kon17]*u[E11P]*u[S] - p[koff17]*fsigmaS*u[E11S]);
    E11PRasGAP_binding = (2*p[kon18]*u[E11P]*u[R] - p[koff18]*fsigmaR*u[E11R]);
    E12PGrb2_binding = (3*p[kon19]*u[E12P]*u[G] - p[koff19]*fsigmaG*u[E12G]);
    E12PShc_binding = (6*p[kon20]*u[E12P]*u[S] - p[koff20]*fsigmaS*u[E12S]);
    E12PRasGAP_binding = (2*p[kon21]*u[E12P]*u[R] - p[koff21]*fsigmaR*u[E12R]);
    E23PGrb2_binding = (3*p[kon22]*u[E23P]*u[G] - p[koff22]*fsigmaG*u[E23G]);
    E23PShc_binding = (3*p[kon23]*u[E23P]*u[S] - p[koff23]*fsigmaS*u[E23S]);
    E23PPI3K_binding = (3*p[kon24]*u[E23P]*u[I] - p[koff24]*u[E23I]);
    E23PRasGAP_binding = (2*p[kon25]*u[E23P]*u[R] - p[koff25]*fsigmaR*u[E23R]);
    E34PGrb2_binding = (4*p[kon26]*u[E34P]*u[G] - p[koff26]*fsigmaG*u[E34G]);
    E34PShc_binding = (3*p[kon27]*u[E34P]*u[S] - p[koff27]*fsigmaS*u[E34S]);
    E34PPI3K_binding = (4*p[kon28]*u[E34P]*u[I] - p[koff28]*u[E34I]);
    E34PRasGAP_binding = (2*p[kon29]*u[E34P]*u[R] - p[koff29]*fsigmaR*u[E34R]);
    E24PGrb2_binding = (3*p[kon30]*u[E24P]*u[G] - p[koff30]*fsigmaG*u[E24G]);
    E24PShc_binding = (4*p[kon31]*u[E24P]*u[S] - p[koff31]*fsigmaS*u[E24S]);
    E24PPI3K_binding = (1*p[kon32]*u[E24P]*u[I] - p[koff32]*u[E24I]);
    E24PRasGAP_binding = (2*p[kon33]*u[E24P]*u[R] - p[koff33]*fsigmaR*u[E24R]);
    E44PRasGAP_binding_1 = (4*p[kon34]*u[E44P]*u[G] - p[koff34]*fsigmaG*u[E44G]);
    E44PShc_binding = (4*p[kon35]*u[E44P]*u[S] - p[koff35]*fsigmaS*u[E44S]);
    E44PPI3K_binding = (2*p[kon36]*u[E44P]*u[I] - p[koff36]*u[E44I]);
    E44PRasGAP_binding_2 = (2*p[kon37]*u[E44P]*u[R] - p[koff37]*fsigmaR*u[E44R]);
    Shc_Phosphorylation = (p[kf38]*u[sigmaS]*sigmaEP - p[VmaxPY]*u[sigmaSP]/(p[KmPY] + u[sigmaSP]) - p[kPTP38]*u[sigmaT]*u[sigmaSP]);
    Gab1_Phosphorylation = (p[kf39]*u[sigmaA]*sigmaEP - p[VmaxPY]*u[sigmaAP]/(p[KmPY] + u[sigmaAP]) - p[kPTP39]*u[sigmaT]*u[sigmaAP]);
    Grb2SOS_binding = (p[kon40]*u[sigmaG]*u[O] - p[koff40]*u[sigmaG_O]);
    Grb2Gab1_binding = (p[kon41]*u[sigmaG]*u[A] - p[koff41]*u[sigmaG_A]*fsigmaA);
    ShcGrb2_binding = (p[kon42]*u[sigmaSP]*u[G] - p[koff42]*u[sigmaSP_G]*fsigmaG);
    Gab1Shc_binding = (3*p[kon43]*u[sigmaAP]*u[S] - p[koff43]*u[sigmaAP_S]*fsigmaS);
    Gab1PI3K_binding = (3*p[kon44]*u[sigmaAP]*u[I] - p[koff44]*u[sigmaAP_I]);
    Gab1RasGAP_binding = (2*p[kon45]*u[sigmaAP]*u[R] - p[koff45]*u[sigmaAP_R]*fsigmaR);
    Gab1PIP3_binding = (p[kon46]*u[P3]*u[A] - p[koff46]*u[P3_A]*fsigmaA);
    Akt_Activation = (p[kf47]*u[P3]*u[Akt]/(p[Kmf47] + u[Akt]) - p[Vmaxr47]*u[Aktstar]/(p[Kmr47] + u[Aktstar]));
    PIP3_Production = (p[kf48]*(1 - u[fint]*f11)*u[sigmaI]*u[P2]/(p[Kmf48] + u[P2]) - 3*p[PTEN]*u[P3]/(p[Kmr48] + u[P3]));
    RasGTP_Production = (p[kf49]*u[sigmaO]*u[RsD]/(p[Kmf49] + u[RsD]) - p[kr49]*u[sigmaR]*u[RsT]/(p[Kmr49] + u[RsT]) - p[kr49b]*u[sigmaRP]*u[RsT]/(p[Kmr49b] + u[RsT]) - p[kcon49]*u[RsT]);
    RasGAP_Phosphorylation = (p[kf50]*u[sigmaR]*sigmaEP - p[VmaxPY]*u[sigmaRP]/(p[KmPY] + u[sigmaRP]) - p[kPTP50]*u[sigmaT]*u[sigmaRP]);
    Raf_Activation = (p[kf51]*u[RsT]*u[Raf]/(p[Kmf51] + u[Raf]) - p[Vmaxr51]*u[Rafstar]/(p[Kmrb51] + u[Rafstar]));
    MEK_Activation = (p[kf52]*u[Rafstar]*u[MEK]/(p[Kmf52] + u[MEK]) - p[Vmaxr52]*u[MEKstar]/(p[Kmr52] + u[MEKstar]));

    SOS_ST_Phosphorylation = (p[kf54]*u[O]*u[ERKstar]/(p[Kmf54] + u[O]) - p[Vmaxr54]*u[OP]/(p[Kmr54] + u[OP]));
    Gab1_ST_Phosphorylation = (p[kf55]*u[A]*u[ERKstar]/(p[Kmf55] + u[A]) - p[Vmaxr55]*u[AP]/(p[Kmr55] + u[AP]));

    Grb2_binding_to_P3_A = (p[kon57]*u[P3_A]*u[G] - p[koff57]*u[sigmaA_G]);
    SOS_binding_to_SigAG = (p[kon58]*u[sigmaA_G]*u[O] - p[koff58]*u[sigmaA_G_O]);
    Gab1_binding_to_SigGO = (p[kon59]*u[sigmaG_O]*u[A] - p[koff59]*u[A_sigmaG_O]*fsigmaA);
    SOS_binding_to_SigGA = (p[kon60]*u[sigmaG_A]*u[O] - p[koff60]*u[A_sigmaG_O]);
    ErbB1ErbB3_Heterodimerization = (p[kon61]*u[H_E3]*u[E_E1] - p[koff61]*u[E13]);
    ErbB1ErbB4_Heterodimerization = (p[kon62]*u[H_E4]*u[E_E1] - p[koff62]*u[E14]);
    E13_Phosphorylation = (p[kf63]*u[E13] - p[VmaxPY]*u[E13P]/(p[KmPY] + u[E13P]) - p[kPTP63]*u[sigmaT]*u[E13P]);
    E14_Phosphorylation = (p[kf64]*u[E14] - p[VmaxPY]*u[E14P]/(p[KmPY] + u[E14P]) - p[kPTP64]*u[sigmaT]*u[E14P]);
    E13Grb2_Binding = (4*p[kon65]*u[E13P]*u[G] - p[koff65]*fsigmaG*u[E13G]);
    E13Shc_Binding = (5*p[kon66]*u[E13P]*u[S] - p[koff66]*fsigmaS*u[E13S]);
    E13PI3K_Binding = (3*p[kon67]*u[E13P]*u[I] - p[koff67]*u[E13I]);
    E13RasGAP_Binding = (2*p[kon68]*u[E13P]*u[R] - p[koff68]*fsigmaR*u[E13R]);
    E14Grb2_Binding = (4*p[kon69]*u[E14P]*u[G] - p[koff69]*fsigmaG*u[E14G]);
    E14Shc_Binding = (6*p[kon70]*u[E14P]*u[S] - p[koff70]*fsigmaS*u[E14S]);
    E14PI3K_Binding = (1*p[kon71]*u[E14P]*u[I] - p[koff71]*u[E14I]);
    E14RasGAP_Binding = (2*p[kon72]*u[E14P]*u[R] - p[koff72]*fsigmaR*u[E14R]);
    E11PTP1B_Binding = (4*p[kon73]*u[E11P]*u[T] - p[koff73]*u[E11T]);
    E12PTP1B_Binding = (3*p[kon74]*u[E12P]*u[T] - p[koff74]*u[E12T]);
    E23PTP1B_Binding = (2*p[kon75]*u[E23P]*u[T] - p[koff75]*u[E23T]);
    E34PTP1B_Binding = (2*p[kon76]*u[E34P]*u[T] - p[koff76]*u[E34T]);
    E24PTP1B_Binding = (2*p[kon77]*u[E24P]*u[T] - p[koff77]*u[E24T]);
    E44PTP1B_Binding = (2*p[kon78]*u[E44P]*u[T] - p[koff78]*u[E44T]);
    E13PTP1B_Binding = (3*p[kon79]*u[E13P]*u[T] - p[koff79]*u[E13T]);
    E14PTP1B_Binding = (3*p[kon80]*u[E14P]*u[T] - p[koff80]*u[E14T]);
    E1_ST_Phosphorylation = (p[kf81]*u[E1]*u[ERKstar]/(p[Kmf81] + u[E1]) - p[Vmaxr81]*u[E1_PT]/(p[Kmr81] + u[E1_PT]));
    E2_ST_Phosphorylation = (p[kf82]*u[E2]*u[ERKstar]/(p[Kmf82] + u[E2]) - p[Vmaxr82]*u[E2_PT]/(p[Kmr82] + u[E2_PT]));
    E4_ST_Phosphorylation = (p[kf83]*u[E4]*u[ERKstar]/(p[Kmf83] + u[E4]) - p[Vmaxr83]*u[E4_PT]/(p[Kmr83] + u[E4_PT]));
    E_E1_ST_Phosphorylation = (p[kf84]*u[E_E1]*u[ERKstar]/(p[Kmf84] + u[E_E1]) - p[Vmaxr84]*u[E_E1_PT]/(p[Kmr84] + u[E_E1_PT]));
    H_E4_ST_Phosphorylation = (p[kf85]*u[H_E4]*u[ERKstar]/(p[Kmf85] + u[H_E4]) - p[Vmaxr85]*u[H_E4_PT]/(p[Kmr85] + u[H_E4_PT]));
    EGF_binding_to_E1PT = (p[kon86]*u[E]*u[E1_PT] - p[EGF_off]*u[E_E1_PT]);
    HRG_binding_to_E4PT = (p[kon87]*u[H]*u[E4_PT] - p[HRGoff_4]*u[H_E4_PT]);
    PTP1B_binding_to_SigAP = (2*p[kon88]*u[sigmaAP]*u[T] - p[koff88]*u[sigmaAP_T]);
    E11P_Degradation = p[kdeg]*u[E11P];
    E11G_Degradation = p[kdeg]*u[E11G];
    E11S_Degradation = p[kdeg]*u[E11S];
    E11R_Degradation = p[kdeg]*u[E11R];
    E11T_Degradation = p[kdeg]*u[E11T];
    PI3K_Fractional_Multiplier = p[a98]*(-u[fint] + p[b98]);
    EGF_in_EC_compartment = (p[kon1]*u[E]*u[E1] - p[EGF_off]*u[E_E1] + p[kon86]*u[E]*u[E1_PT] - p[EGF_off]*u[E_E1_PT])/VeVc;
    HRG_in_EC_compartment = (p[kon87]*u[H]*u[E4_PT] - p[HRGoff_4]*u[H_E4_PT] + p[kon2]*u[H]*u[E3] - p[HRGoff_3]*u[H_E3] + p[kon3]*u[H]*u[E4] - p[HRGoff_4]*u[H_E4])/VeVc;
    ERK_binding_to_MEKstar_1 = (p[kon89]*u[ERK]*u[MEKstar] - p[koff89]*u[ERK_MEKstar]);
    pERK_production = p[kcat90]*u[ERK_MEKstar];
    ERK_binding_to_MEKstar_2 = (p[kon91]*u[pERK]*u[MEKstar] - p[koff91]*u[pERK_MEKstar]);
    ERKstar_production = p[kcat92]*u[pERK_MEKstar];
    ERKstar_binding_to_Phosphatase = (p[kon93]*u[ERKstar]*u[ERKPpase] - p[koff93]*u[ERKstar_ERKPpase]);
    ERKstar_dephosphorylation = p[kcat94]*u[ERKstar_ERKPpase];
    pERK_binding_to_Phosphatase = (p[kon95]*u[pERK]*u[ERKPpase] - p[koff95]*u[pERK_ERKPpase]);
    pERK_dephosphorylation = p[kcat96]*u[pERK_ERKPpase];

    du = similar(u);

    du[E1] = ((-1.0 * EGF_Binding_to_ErbB1) + (-1.0 * E1_ST_Phosphorylation));
    du[E2] = ((-1.0 * ErbB1ErbB2_Heterodimerization) + (-1.0 * ErbB2ErbB3_Heterodimerization) + (-1.0 * ErbB2ErbB4_Heterodimerization) + (-1.0 * E2_ST_Phosphorylation));
    du[E3] = ((-1.0 * HRG_Binding_to_ErbB3));
    du[E4] = ((-1.0 * HRG_Binding_to_ErbB4) + (-1.0 * E4_ST_Phosphorylation));
    du[E_E1] = (( 1.0 * EGF_Binding_to_ErbB1) + (-1.0 * ErbB1_Dimerization) + (-1.0 * ErbB1_Dimerization) + (-1.0 * ErbB1ErbB2_Heterodimerization) + (-1.0 * ErbB1ErbB3_Heterodimerization) + (-1.0 * ErbB1ErbB4_Heterodimerization) + (-1.0 * E_E1_ST_Phosphorylation));
    du[H_E3] = (( 1.0 * HRG_Binding_to_ErbB3) + (-1.0 * ErbB2ErbB3_Heterodimerization) + (-1.0 * ErbB3ErbB4_Heterodimerization) + (-1.0 * ErbB1ErbB3_Heterodimerization));
    du[H_E4] = (( 1.0 * HRG_Binding_to_ErbB4) + (-1.0 * ErbB3ErbB4_Heterodimerization) + (-1.0 * ErbB2ErbB4_Heterodimerization) + (-1.0 * ErbB4_Homodimerization) + (-1.0 * ErbB4_Homodimerization) + (-1.0 * ErbB1ErbB4_Heterodimerization) + (-1.0 * H_E4_ST_Phosphorylation));
    du[E11] = (( 1.0 * ErbB1_Dimerization) + (-1.0 * E11_Phosphorylation));
    du[E12] = (( 1.0 * ErbB1ErbB2_Heterodimerization) + (-1.0 * E12_Phosphorylation));
    du[E23] = (( 1.0 * ErbB2ErbB3_Heterodimerization) + (-1.0 * E23_Phosphorylation));
    du[E34] = (( 1.0 * ErbB3ErbB4_Heterodimerization) + (-1.0 * E34_Phosphorylation));
    du[E24] = (( 1.0 * ErbB2ErbB4_Heterodimerization) + (-1.0 * E24_Phosphorylation));
    du[E44] = (( 1.0 * ErbB4_Homodimerization) + (-1.0 * E44_Phosphorylation));
    du[E11P] = (( 1.0 * E11_Phosphorylation) + (-1.0 * E11PGrb2_binding) + (-1.0 * E11PShc_binding) + (-1.0 * E11PRasGAP_binding) + (-1.0 * E11PTP1B_Binding) + (-1.0 * E11P_Degradation));
    du[E12P] = (( 1.0 * E12_Phosphorylation) + (-1.0 * E12PGrb2_binding) + (-1.0 * E12PShc_binding) + (-1.0 * E12PRasGAP_binding) + (-1.0 * E12PTP1B_Binding));
    du[E23P] = (( 1.0 * E23_Phosphorylation) + (-1.0 * E23PGrb2_binding) + (-1.0 * E23PShc_binding) + (-1.0 * E23PPI3K_binding) + (-1.0 * E23PRasGAP_binding) + (-1.0 * E23PTP1B_Binding));
    du[E34P] = (( 1.0 * E34_Phosphorylation) + (-1.0 * E34PGrb2_binding) + (-1.0 * E34PShc_binding) + (-1.0 * E34PPI3K_binding) + (-1.0 * E34PRasGAP_binding) + (-1.0 * E34PTP1B_Binding));
    du[E24P] = (( 1.0 * E24_Phosphorylation) + (-1.0 * E24PGrb2_binding) + (-1.0 * E24PShc_binding) + (-1.0 * E24PPI3K_binding) + (-1.0 * E24PRasGAP_binding) + (-1.0 * E24PTP1B_Binding));
    du[E44P] = (( 1.0 * E44_Phosphorylation) + (-1.0 * E44PRasGAP_binding_1) + (-1.0 * E44PShc_binding) + (-1.0 * E44PPI3K_binding) + (-1.0 * E44PRasGAP_binding_2) + (-1.0 * E44PTP1B_Binding));
    du[G] = ((-1.0 * E11PGrb2_binding) + (-1.0 * E12PGrb2_binding) + (-1.0 * E23PGrb2_binding) + (-1.0 * E34PGrb2_binding) + (-1.0 * E24PGrb2_binding) + (-1.0 * E44PRasGAP_binding_1) + (-1.0 * ShcGrb2_binding) + (-1.0 * Grb2_binding_to_P3_A) + (-1.0 * E13Grb2_Binding) + (-1.0 * E14Grb2_Binding) + ( 1.0 * E11G_Degradation));
    du[S] = ((-1.0 * E11PShc_binding) + (-1.0 * E12PShc_binding) + (-1.0 * E23PShc_binding) + (-1.0 * E34PShc_binding) + (-1.0 * E24PShc_binding) + (-1.0 * E44PShc_binding) + (-1.0 * Gab1Shc_binding) + (-1.0 * E13Shc_Binding) + (-1.0 * E14Shc_Binding) + ( 1.0 * E11S_Degradation));
    du[I] = ((-1.0 * E23PPI3K_binding) + (-1.0 * E34PPI3K_binding) + (-1.0 * E24PPI3K_binding) + (-1.0 * E44PPI3K_binding) + (-1.0 * Gab1PI3K_binding) + (-1.0 * E13PI3K_Binding) + (-1.0 * E14PI3K_Binding));
    du[R] = ((-1.0 * E11PRasGAP_binding) + (-1.0 * E12PRasGAP_binding) + (-1.0 * E23PRasGAP_binding) + (-1.0 * E34PRasGAP_binding) + (-1.0 * E24PRasGAP_binding) + (-1.0 * E44PRasGAP_binding_2) + (-1.0 * Gab1RasGAP_binding) + (-1.0 * E13RasGAP_Binding) + (-1.0 * E14RasGAP_Binding) + ( 1.0 * E11R_Degradation));
    du[O] = ((-1.0 * Grb2SOS_binding) + (-1.0 * SOS_ST_Phosphorylation) + (-1.0 * SOS_binding_to_SigAG) + (-1.0 * SOS_binding_to_SigGA));
    du[A] = ((-1.0 * Grb2Gab1_binding) + (-1.0 * Gab1PIP3_binding) + (-1.0 * Gab1_ST_Phosphorylation) + (-1.0 * Gab1_binding_to_SigGO));
    du[E11G] = (( 1.0 * E11PGrb2_binding) + (-1.0 * E11G_Degradation));
    du[E11S] = (( 1.0 * E11PShc_binding) + (-1.0 * E11S_Degradation));
    du[E11R] = (( 1.0 * E11PRasGAP_binding) + (-1.0 * E11R_Degradation));
    du[E12G] = (( 1.0 * E12PGrb2_binding));
    du[E12S] = (( 1.0 * E12PShc_binding));
    du[E12R] = (( 1.0 * E12PRasGAP_binding));
    du[E23G] = (( 1.0 * E23PGrb2_binding));
    du[E23S] = (( 1.0 * E23PShc_binding));
    du[E23I] = (( 1.0 * E23PPI3K_binding));
    du[E23R] = (( 1.0 * E23PRasGAP_binding));
    du[E34G] = (( 1.0 * E34PGrb2_binding));
    du[E34S] = (( 1.0 * E34PShc_binding));
    du[E34I] = (( 1.0 * E34PPI3K_binding));
    du[E34R] = (( 1.0 * E34PRasGAP_binding));
    du[E24G] = (( 1.0 * E24PGrb2_binding));
    du[E24S] = (( 1.0 * E24PShc_binding));
    du[E24I] = (( 1.0 * E24PPI3K_binding));
    du[E24R] = (( 1.0 * E24PRasGAP_binding));
    du[E44G] = (( 1.0 * E44PRasGAP_binding_1));
    du[E44S] = (( 1.0 * E44PShc_binding));
    du[E44I] = (( 1.0 * E44PPI3K_binding));
    du[E44R] = (( 1.0 * E44PRasGAP_binding_2));
    du[sigmaG] = (( 1.0 * E11PGrb2_binding) + ( 1.0 * E12PGrb2_binding) + ( 1.0 * E23PGrb2_binding) + ( 1.0 * E34PGrb2_binding) + ( 1.0 * E24PGrb2_binding) + ( 1.0 * E44PRasGAP_binding_1) + (-1.0 * Grb2SOS_binding) + (-1.0 * Grb2Gab1_binding) + ( 1.0 * ShcGrb2_binding) + ( 1.0 * E13Grb2_Binding) + ( 1.0 * E14Grb2_Binding) + (-1.0 * E11G_Degradation));
    du[sigmaS] = (( 1.0 * E11PShc_binding) + ( 1.0 * E12PShc_binding) + ( 1.0 * E23PShc_binding) + ( 1.0 * E34PShc_binding) + ( 1.0 * E24PShc_binding) + ( 1.0 * E44PShc_binding) + (-1.0 * Shc_Phosphorylation) + ( 1.0 * Gab1Shc_binding) + ( 1.0 * E13Shc_Binding) + ( 1.0 * E14Shc_Binding) + (-1.0 * E11S_Degradation));
    du[sigmaI] = (( 1.0 * E23PPI3K_binding) + ( 1.0 * E34PPI3K_binding) + ( 1.0 * E24PPI3K_binding) + ( 1.0 * E44PPI3K_binding) + ( 1.0 * Gab1PI3K_binding) + ( 1.0 * E13PI3K_Binding) + ( 1.0 * E14PI3K_Binding));
    du[sigmaR] = (( 1.0 * E11PRasGAP_binding) + ( 1.0 * E12PRasGAP_binding) + ( 1.0 * E23PRasGAP_binding) + ( 1.0 * E34PRasGAP_binding) + ( 1.0 * E24PRasGAP_binding) + ( 1.0 * E44PRasGAP_binding_2) + ( 1.0 * Gab1RasGAP_binding) + (-1.0 * RasGAP_Phosphorylation) + ( 1.0 * E13RasGAP_Binding) + ( 1.0 * E14RasGAP_Binding) + (-1.0 * E11R_Degradation));
    du[sigmaA] = ((-1.0 * Gab1_Phosphorylation) + ( 1.0 * Grb2Gab1_binding) + ( 1.0 * Gab1PIP3_binding) + ( 1.0 * Gab1_binding_to_SigGO));
    du[sigmaSP] = (( 1.0 * Shc_Phosphorylation) + (-1.0 * ShcGrb2_binding));
    du[sigmaAP] = (( 1.0 * Gab1_Phosphorylation) + (-1.0 * Gab1Shc_binding) + (-1.0 * Gab1PI3K_binding) + (-1.0 * Gab1RasGAP_binding) + (-1.0 * PTP1B_binding_to_SigAP));
    du[sigmaG_O] = (( 1.0 * Grb2SOS_binding) + (-1.0 * Gab1_binding_to_SigGO));
    du[sigmaG_A] = (( 1.0 * Grb2Gab1_binding) + (-1.0 * SOS_binding_to_SigGA));
    du[sigmaSP_G] = (( 1.0 * ShcGrb2_binding));
    du[sigmaAP_S] = (( 1.0 * Gab1Shc_binding));
    du[sigmaAP_I] = (( 1.0 * Gab1PI3K_binding));
    du[sigmaAP_R] = (( 1.0 * Gab1RasGAP_binding));
    du[Empty] = (( 1.0 * E11P_Degradation) + (-1.0 * PI3K_Fractional_Multiplier) + ( 1.0 * EGF_in_EC_compartment) + ( 1.0 * HRG_in_EC_compartment));
    du[P3_A] = (( 1.0 * Gab1PIP3_binding) + (-1.0 * Grb2_binding_to_P3_A));
    du[P2] = ((-1.0 * PIP3_Production));
    du[P3] = ((-1.0 * Gab1PIP3_binding) + ( 1.0 * PIP3_Production));
    du[Akt] = ((-1.0 * Akt_Activation));
    du[RsD] = ((-1.0 * RasGTP_Production));
    du[RsT] = (( 1.0 * RasGTP_Production));
    du[sigmaRP] = (( 1.0 * RasGAP_Phosphorylation));
    du[Raf] = ((-1.0 * Raf_Activation));
    du[Rafstar] = (( 1.0 * Raf_Activation));
    du[MEK] = ((-1.0 * MEK_Activation));
    du[MEKstar] = (( 1.0 * MEK_Activation) + (-1.0 * ERK_binding_to_MEKstar_1) + ( 1.0 * pERK_production) + (-1.0 * ERK_binding_to_MEKstar_2) + ( 1.0 * ERKstar_production));
    du[ERK] = ((-1.0 * ERK_binding_to_MEKstar_1) + ( 1.0 * pERK_dephosphorylation));
    du[ERKstar] = (( 1.0 * ERKstar_production) + (-1.0 * ERKstar_binding_to_Phosphatase));
    du[OP] = (( 1.0 * SOS_ST_Phosphorylation));
    du[AP] = (( 1.0 * Gab1_ST_Phosphorylation));
    du[A_sigmaG_O] = (( 1.0 * Gab1_binding_to_SigGO) + ( 1.0 * SOS_binding_to_SigGA));
    du[sigmaA_G] = (( 1.0 * Grb2_binding_to_P3_A) + (-1.0 * SOS_binding_to_SigAG));
    du[sigmaA_G_O] = (( 1.0 * SOS_binding_to_SigAG));
    du[sigmaO] = (( 1.0 * Grb2SOS_binding) + ( 1.0 * SOS_binding_to_SigAG) + ( 1.0 * SOS_binding_to_SigGA));
    du[E13] = (( 1.0 * ErbB1ErbB3_Heterodimerization) + (-1.0 * E13_Phosphorylation));
    du[E14] = (( 1.0 * ErbB1ErbB4_Heterodimerization) + (-1.0 * E14_Phosphorylation));
    du[E13P] = (( 1.0 * E13_Phosphorylation) + (-1.0 * E13Grb2_Binding) + (-1.0 * E13Shc_Binding) + (-1.0 * E13PI3K_Binding) + (-1.0 * E13RasGAP_Binding) + (-1.0 * E13PTP1B_Binding));
    du[E14P] = (( 1.0 * E14_Phosphorylation) + (-1.0 * E14Grb2_Binding) + (-1.0 * E14Shc_Binding) + (-1.0 * E14PI3K_Binding) + (-1.0 * E14RasGAP_Binding) + (-1.0 * E14PTP1B_Binding));
    du[E13G] = (( 1.0 * E13Grb2_Binding));
    du[E13S] = (( 1.0 * E13Shc_Binding));
    du[E13I] = (( 1.0 * E13PI3K_Binding));
    du[E13R] = (( 1.0 * E13RasGAP_Binding));
    du[E14G] = (( 1.0 * E14Grb2_Binding));
    du[E14S] = (( 1.0 * E14Shc_Binding));
    du[E14I] = (( 1.0 * E14PI3K_Binding));
    du[E14R] = (( 1.0 * E14RasGAP_Binding));
    du[T] = ((-1.0 * E11PTP1B_Binding) + (-1.0 * E12PTP1B_Binding) + (-1.0 * E23PTP1B_Binding) + (-1.0 * E34PTP1B_Binding) + (-1.0 * E24PTP1B_Binding) + (-1.0 * E44PTP1B_Binding) + (-1.0 * E13PTP1B_Binding) + (-1.0 * E14PTP1B_Binding) + (-1.0 * PTP1B_binding_to_SigAP) + ( 1.0 * E11T_Degradation));
    du[E11T] = (( 1.0 * E11PTP1B_Binding) + (-1.0 * E11T_Degradation));
    du[E12T] = (( 1.0 * E12PTP1B_Binding));
    du[E23T] = (( 1.0 * E23PTP1B_Binding));
    du[E34T] = (( 1.0 * E34PTP1B_Binding));
    du[E24T] = (( 1.0 * E24PTP1B_Binding));
    du[E44T] = (( 1.0 * E44PTP1B_Binding));
    du[E13T] = (( 1.0 * E13PTP1B_Binding));
    du[E14T] = (( 1.0 * E14PTP1B_Binding));
    du[sigmaT] = (( 1.0 * E11PTP1B_Binding) + ( 1.0 * E12PTP1B_Binding) + ( 1.0 * E23PTP1B_Binding) + ( 1.0 * E34PTP1B_Binding) + ( 1.0 * E24PTP1B_Binding) + ( 1.0 * E44PTP1B_Binding) + ( 1.0 * E13PTP1B_Binding) + ( 1.0 * E14PTP1B_Binding) + ( 1.0 * PTP1B_binding_to_SigAP) + (-1.0 * E11T_Degradation));
    du[E1_PT] = (( 1.0 * E1_ST_Phosphorylation) + (-1.0 * EGF_binding_to_E1PT));
    du[E2_PT] = (( 1.0 * E2_ST_Phosphorylation));
    du[E4_PT] = (( 1.0 * E4_ST_Phosphorylation) + (-1.0 * HRG_binding_to_E4PT));
    du[E_E1_PT] = (( 1.0 * E_E1_ST_Phosphorylation) + ( 1.0 * EGF_binding_to_E1PT));
    du[H_E4_PT] = (( 1.0 * H_E4_ST_Phosphorylation) + ( 1.0 * HRG_binding_to_E4PT));
    du[Aktstar] = (( 1.0 * Akt_Activation));
    du[sigmaAP_T] = (( 1.0 * PTP1B_binding_to_SigAP));
    du[E] = ((-1.0 * EGF_in_EC_compartment));
    du[H] = ((-1.0 * HRG_in_EC_compartment));
    du[fint] = (( 1.0 * PI3K_Fractional_Multiplier));
    du[pERK] = (( 1.0 * pERK_production) + (-1.0 * ERK_binding_to_MEKstar_2) + ( 1.0 * ERKstar_dephosphorylation) + (-1.0 * pERK_binding_to_Phosphatase));
    du[ERK_MEKstar] = (( 1.0 * ERK_binding_to_MEKstar_1) + (-1.0 * pERK_production));
    du[pERK_MEKstar] = (( 1.0 * ERK_binding_to_MEKstar_2) + (-1.0 * ERKstar_production));
    du[pERK_ERKPpase] = (( 1.0 * pERK_binding_to_Phosphatase) + (-1.0 * pERK_dephosphorylation));
    du[ERKPpase] = ((-1.0 * ERKstar_binding_to_Phosphatase) + ( 1.0 * ERKstar_dephosphorylation) + (-1.0 * pERK_binding_to_Phosphatase) + ( 1.0 * pERK_dephosphorylation));
    du[ERKstar_ERKPpase] = (( 1.0 * ERKstar_binding_to_Phosphatase) + (-1.0 * ERKstar_dephosphorylation));

    return du
end