from .name2idx import parameters as C
from .name2idx import variables as V

def diffeq(t,y,x):

    # extracellular volume to cytoplasmic volume ratio
    VeVc = 33.3

    # fraction definitions
    if y[V.sigmaS] + y[V.sigmaSP] + y[V.sigmaSP_G] > 0.0:
        fsigmaS = y[V.sigmaS]/(y[V.sigmaS] + y[V.sigmaSP] + y[V.sigmaSP_G])
    else:
        fsigmaS = 0.0
    

    if y[V.sigmaG] + y[V.sigmaG_A] + y[V.sigmaG_O] + y[V.A_sigmaG_O] > 0.0:
        fsigmaG = y[V.sigmaG]/(y[V.sigmaG] + y[V.sigmaG_A] + y[V.sigmaG_O] + y[V.A_sigmaG_O])
    else:
        fsigmaG = 0.0
    

    if y[V.sigmaA] + y[V.sigmaAP] + y[V.sigmaAP_S] + y[V.sigmaAP_R] + y[V.sigmaAP_I] + y[V.sigmaAP_T] > 0.0:
        fsigmaA = y[V.sigmaA]/(y[V.sigmaA] + y[V.sigmaAP] + y[V.sigmaAP_S] + y[V.sigmaAP_R] + y[V.sigmaAP_I] + y[V.sigmaAP_T])
    else:
        fsigmaA = 0.0
    

    if y[V.sigmaR] + y[V.sigmaRP] > 0.0:
        fsigmaR = y[V.sigmaR]/(y[V.sigmaR] + y[V.sigmaRP])
    else:
        fsigmaR = 0.0
    

    sigmaEP = y[V.E11P] + y[V.E12P] + y[V.E23P] + y[V.E24P] + y[V.E34P] + y[V.E44P] + y[V.E13P] + y[V.E14P]
    if sigmaEP > 0.0:
        f11 = y[V.E11P]/sigmaEP
    else:
        f11 = 0.0
    

    # E11 degradation
    v10d = x[C.kdeg]*y[V.E11P]
    v16d = x[C.kdeg]*y[V.E11G]
    v17d = x[C.kdeg]*y[V.E11S]
    v18d = x[C.kdeg]*y[V.E11R]
    v73d = x[C.kdeg]*y[V.E11T]

    v = [0]*97

    v[1] = (x[C.kon1]*y[V.E]*y[V.E1] - x[C.EGF_off]*y[V.E_E1])
    v[2] = (x[C.kon2]*y[V.H]*y[V.E3] - x[C.HRGoff_3]*y[V.H_E3])
    v[3] = (x[C.kon3]*y[V.H]*y[V.E4] - x[C.HRGoff_4]*y[V.H_E4])
    v[4] = (x[C.kon4]*y[V.E_E1]*y[V.E_E1] - x[C.koff4]*y[V.E11])
    v[5] = (x[C.kon5]*y[V.E_E1]*y[V.E2] - x[C.koff5]*y[V.E12])
    v[6] = (x[C.kon6]*y[V.H_E3]*y[V.E2] - x[C.koff6]*y[V.E23])
    v[7] = (x[C.kon7]*y[V.H_E3]*y[V.H_E4] - x[C.koff7]*y[V.E34])
    v[8] = (x[C.kon8]*y[V.H_E4]*y[V.E2] - x[C.koff8]*y[V.E24])
    v[9] = (x[C.kon9]*y[V.H_E4]*y[V.H_E4] - x[C.koff9]*y[V.E44])
    v[10] = (x[C.kf10]*y[V.E11] - x[C.VmaxPY]*y[V.E11P]/(x[C.KmPY] + y[V.E11P]) - x[C.kPTP10]*y[V.sigmaT]*y[V.E11P])
    v[11] = (x[C.kf11]*y[V.E12] - x[C.VmaxPY]*y[V.E12P]/(x[C.KmPY] + y[V.E12P]) - x[C.kPTP11]*y[V.sigmaT]*y[V.E12P])
    v[12] = (x[C.kf12]*y[V.E23] - x[C.VmaxPY]*y[V.E23P]/(x[C.KmPY] + y[V.E23P]) - x[C.kPTP12]*y[V.sigmaT]*y[V.E23P])
    v[13] = (x[C.kf13]*y[V.E34] - x[C.VmaxPY]*y[V.E34P]/(x[C.KmPY] + y[V.E34P]) - x[C.kPTP13]*y[V.sigmaT]*y[V.E34P])
    v[14] = (x[C.kf14]*y[V.E24] - x[C.VmaxPY]*y[V.E24P]/(x[C.KmPY] + y[V.E24P]) - x[C.kPTP14]*y[V.sigmaT]*y[V.E24P])
    v[15] = (x[C.kf15]*y[V.E44] - x[C.VmaxPY]*y[V.E44P]/(x[C.KmPY] + y[V.E44P]) - x[C.kPTP15]*y[V.sigmaT]*y[V.E44P])
    v[16] = (4*x[C.kon16]*y[V.E11P]*y[V.G] - x[C.koff16]*fsigmaG*y[V.E11G])
    v[17] = (8*x[C.kon17]*y[V.E11P]*y[V.S] - x[C.koff17]*fsigmaS*y[V.E11S])
    v[18] = (2*x[C.kon18]*y[V.E11P]*y[V.R] - x[C.koff18]*fsigmaR*y[V.E11R])
    v[19] = (3*x[C.kon19]*y[V.E12P]*y[V.G] - x[C.koff19]*fsigmaG*y[V.E12G])
    v[20] = (6*x[C.kon20]*y[V.E12P]*y[V.S] - x[C.koff20]*fsigmaS*y[V.E12S])
    v[21] = (2*x[C.kon21]*y[V.E12P]*y[V.R] - x[C.koff21]*fsigmaR*y[V.E12R])
    v[22] = (3*x[C.kon22]*y[V.E23P]*y[V.G] - x[C.koff22]*fsigmaG*y[V.E23G])
    v[23] = (3*x[C.kon23]*y[V.E23P]*y[V.S] - x[C.koff23]*fsigmaS*y[V.E23S])
    v[24] = (3*x[C.kon24]*y[V.E23P]*y[V.I] - x[C.koff24]*y[V.E23I])
    v[25] = (2*x[C.kon25]*y[V.E23P]*y[V.R] - x[C.koff25]*fsigmaR*y[V.E23R])
    v[26] = (4*x[C.kon26]*y[V.E34P]*y[V.G] - x[C.koff26]*fsigmaG*y[V.E34G])
    v[27] = (3*x[C.kon27]*y[V.E34P]*y[V.S] - x[C.koff27]*fsigmaS*y[V.E34S])
    v[28] = (4*x[C.kon28]*y[V.E34P]*y[V.I] - x[C.koff28]*y[V.E34I])
    v[29] = (2*x[C.kon29]*y[V.E34P]*y[V.R] - x[C.koff29]*fsigmaR*y[V.E34R])
    v[30] = (3*x[C.kon30]*y[V.E24P]*y[V.G] - x[C.koff30]*fsigmaG*y[V.E24G])
    v[31] = (4*x[C.kon31]*y[V.E24P]*y[V.S] - x[C.koff31]*fsigmaS*y[V.E24S])
    v[32] = (1*x[C.kon32]*y[V.E24P]*y[V.I] - x[C.koff32]*y[V.E24I])
    v[33] = (2*x[C.kon33]*y[V.E24P]*y[V.R] - x[C.koff33]*fsigmaR*y[V.E24R])
    v[34] = (4*x[C.kon34]*y[V.E44P]*y[V.G] - x[C.koff34]*fsigmaG*y[V.E44G])
    v[35] = (4*x[C.kon35]*y[V.E44P]*y[V.S] - x[C.koff35]*fsigmaS*y[V.E44S])
    v[36] = (2*x[C.kon36]*y[V.E44P]*y[V.I] - x[C.koff36]*y[V.E44I])
    v[37] = (2*x[C.kon37]*y[V.E44P]*y[V.R] - x[C.koff37]*fsigmaR*y[V.E44R])
    v[38] = (x[C.kf38]*y[V.sigmaS]*sigmaEP - x[C.VmaxPY]*y[V.sigmaSP]/(x[C.KmPY] + y[V.sigmaSP]) - x[C.kPTP38]*y[V.sigmaT]*y[V.sigmaSP])
    v[39] = (x[C.kf39]*y[V.sigmaA]*sigmaEP - x[C.VmaxPY]*y[V.sigmaAP]/(x[C.KmPY] + y[V.sigmaAP]) - x[C.kPTP39]*y[V.sigmaT]*y[V.sigmaAP])
    v[40] = (x[C.kon40]*y[V.sigmaG]*y[V.O] - x[C.koff40]*y[V.sigmaG_O])
    v[41] = (x[C.kon41]*y[V.sigmaG]*y[V.A] - x[C.koff41]*y[V.sigmaG_A]*fsigmaA)
    v[42] = (x[C.kon42]*y[V.sigmaSP]*y[V.G] - x[C.koff42]*y[V.sigmaSP_G]*fsigmaG)
    v[43] = (3*x[C.kon43]*y[V.sigmaAP]*y[V.S] - x[C.koff43]*y[V.sigmaAP_S]*fsigmaS)
    v[44] = (3*x[C.kon44]*y[V.sigmaAP]*y[V.I] - x[C.koff44]*y[V.sigmaAP_I])
    v[45] = (2*x[C.kon45]*y[V.sigmaAP]*y[V.R] - x[C.koff45]*y[V.sigmaAP_R]*fsigmaR)
    v[46] = (x[C.kon46]*y[V.P3]*y[V.A] - x[C.koff46]*y[V.P3_A]*fsigmaA)
    v[47] = (x[C.kf47]*y[V.P3]*y[V.Akt]/(x[C.Kmf47] + y[V.Akt]) - x[C.Vmaxr47]*y[V.Aktstar]/(x[C.Kmr47] + y[V.Aktstar]))
    v[48] = (x[C.kf48]*(1 - y[V.fint]*f11)*y[V.sigmaI]*y[V.P2]/(x[C.Kmf48] + y[V.P2]) - 3*x[C.PTEN]*y[V.P3]/(x[C.Kmr48] + y[V.P3]))
    v[49] = (x[C.kf49]*y[V.sigmaO]*y[V.RsD]/(x[C.Kmf49] + y[V.RsD]) - x[C.kr49]*y[V.sigmaR]*y[V.RsT]/(x[C.Kmr49] + y[V.RsT]) - x[C.kr49b]*y[V.sigmaRP]*y[V.RsT]/(x[C.Kmr49b] + y[V.RsT]) - x[C.kcon49]*y[V.RsT])
    v[50] = (x[C.kf50]*y[V.sigmaR]*sigmaEP - x[C.VmaxPY]*y[V.sigmaRP]/(x[C.KmPY] + y[V.sigmaRP]) - x[C.kPTP50]*y[V.sigmaT]*y[V.sigmaRP])
    v[51] = (x[C.kf51]*y[V.RsT]*y[V.Raf]/(x[C.Kmf51] + y[V.Raf]) - x[C.Vmaxr51]*y[V.Rafstar]/(x[C.Kmrb51] + y[V.Rafstar]))
    v[52] = (x[C.kf52]*y[V.Rafstar]*y[V.MEK]/(x[C.Kmf52] + y[V.MEK]) - x[C.Vmaxr52]*y[V.MEKstar]/(x[C.Kmr52] + y[V.MEKstar]))
    # v[53]: No reaction
    v[54] = (x[C.kf54]*y[V.O]*y[V.ERKstar]/(x[C.Kmf54] + y[V.O]) - x[C.Vmaxr54]*y[V.OP]/(x[C.Kmr54] + y[V.OP]))
    v[55] = (x[C.kf55]*y[V.A]*y[V.ERKstar]/(x[C.Kmf55] + y[V.A]) - x[C.Vmaxr55]*y[V.AP]/(x[C.Kmr55] + y[V.AP]))
    # v[56]: No reaction
    v[57] = (x[C.kon57]*y[V.P3_A]*y[V.G] - x[C.koff57]*y[V.sigmaA_G])
    v[58] = (x[C.kon58]*y[V.sigmaA_G]*y[V.O] - x[C.koff58]*y[V.sigmaA_G_O])
    v[59] = (x[C.kon59]*y[V.sigmaG_O]*y[V.A] - x[C.koff59]*y[V.A_sigmaG_O]*fsigmaA)
    v[60] = (x[C.kon60]*y[V.sigmaG_A]*y[V.O] - x[C.koff60]*y[V.A_sigmaG_O])
    v[61] = (x[C.kon61]*y[V.H_E3]*y[V.E_E1] - x[C.koff61]*y[V.E13])
    v[62] = (x[C.kon62]*y[V.H_E4]*y[V.E_E1] - x[C.koff62]*y[V.E14])
    v[63] = (x[C.kf63]*y[V.E13] - x[C.VmaxPY]*y[V.E13P]/(x[C.KmPY] + y[V.E13P]) - x[C.kPTP63]*y[V.sigmaT]*y[V.E13P])
    v[64] = (x[C.kf64]*y[V.E14] - x[C.VmaxPY]*y[V.E14P]/(x[C.KmPY] + y[V.E14P]) - x[C.kPTP64]*y[V.sigmaT]*y[V.E14P])
    v[65] = (4*x[C.kon65]*y[V.E13P]*y[V.G] - x[C.koff65]*fsigmaG*y[V.E13G])
    v[66] = (5*x[C.kon66]*y[V.E13P]*y[V.S] - x[C.koff66]*fsigmaS*y[V.E13S])
    v[67] = (3*x[C.kon67]*y[V.E13P]*y[V.I] - x[C.koff67]*y[V.E13I])
    v[68] = (2*x[C.kon68]*y[V.E13P]*y[V.R] - x[C.koff68]*fsigmaR*y[V.E13R])
    v[69] = (4*x[C.kon69]*y[V.E14P]*y[V.G] - x[C.koff69]*fsigmaG*y[V.E14G])
    v[70] = (6*x[C.kon70]*y[V.E14P]*y[V.S] - x[C.koff70]*fsigmaS*y[V.E14S])
    v[71] = (1*x[C.kon71]*y[V.E14P]*y[V.I] - x[C.koff71]*y[V.E14I])
    v[72] = (2*x[C.kon72]*y[V.E14P]*y[V.R] - x[C.koff72]*fsigmaR*y[V.E14R])
    v[73] = (4*x[C.kon73]*y[V.E11P]*y[V.T] - x[C.koff73]*y[V.E11T])
    v[74] = (3*x[C.kon74]*y[V.E12P]*y[V.T] - x[C.koff74]*y[V.E12T])
    v[75] = (2*x[C.kon75]*y[V.E23P]*y[V.T] - x[C.koff75]*y[V.E23T])
    v[76] = (2*x[C.kon76]*y[V.E34P]*y[V.T] - x[C.koff76]*y[V.E34T])
    v[77] = (2*x[C.kon77]*y[V.E24P]*y[V.T] - x[C.koff77]*y[V.E24T])
    v[78] = (2*x[C.kon78]*y[V.E44P]*y[V.T] - x[C.koff78]*y[V.E44T])
    v[79] = (3*x[C.kon79]*y[V.E13P]*y[V.T] - x[C.koff79]*y[V.E13T])
    v[80] = (3*x[C.kon80]*y[V.E14P]*y[V.T] - x[C.koff80]*y[V.E14T])
    v[81] = (x[C.kf81]*y[V.E1]*y[V.ERKstar]/(x[C.Kmf81] + y[V.E1]) - x[C.Vmaxr81]*y[V.E1_PT]/(x[C.Kmr81] + y[V.E1_PT]))
    v[82] = (x[C.kf82]*y[V.E2]*y[V.ERKstar]/(x[C.Kmf82] + y[V.E2]) - x[C.Vmaxr82]*y[V.E2_PT]/(x[C.Kmr82] + y[V.E2_PT]))
    v[83] = (x[C.kf83]*y[V.E4]*y[V.ERKstar]/(x[C.Kmf83] + y[V.E4]) - x[C.Vmaxr83]*y[V.E4_PT]/(x[C.Kmr83] + y[V.E4_PT]))
    v[84] = (x[C.kf84]*y[V.E_E1]*y[V.ERKstar]/(x[C.Kmf84] + y[V.E_E1]) - x[C.Vmaxr84]*y[V.E_E1_PT]/(x[C.Kmr84] + y[V.E_E1_PT]))
    v[85] = (x[C.kf85]*y[V.H_E4]*y[V.ERKstar]/(x[C.Kmf85] + y[V.H_E4]) - x[C.Vmaxr85]*y[V.H_E4_PT]/(x[C.Kmr85] + y[V.H_E4_PT]))
    v[86] = (x[C.kon86]*y[V.E]*y[V.E1_PT] - x[C.EGF_off]*y[V.E_E1_PT])
    v[87] = (x[C.kon87]*y[V.H]*y[V.E4_PT] - x[C.HRGoff_4]*y[V.H_E4_PT])
    v[88] = (2*x[C.kon88]*y[V.sigmaAP]*y[V.T] - x[C.koff88]*y[V.sigmaAP_T])
    v[89] = (x[C.kon89]*y[V.ERK]*y[V.MEKstar] - x[C.koff89]*y[V.ERK_MEKstar])
    v[90] = x[C.kcat90]*y[V.ERK_MEKstar]
    v[91] = (x[C.kon91]*y[V.pERK]*y[V.MEKstar] - x[C.koff91]*y[V.pERK_MEKstar])
    v[92] = x[C.kcat92]*y[V.pERK_MEKstar]
    v[93] = (x[C.kon93]*y[V.ERKstar]*y[V.ERKPpase] - x[C.koff93]*y[V.ERKstar_ERKPpase])
    v[94] = x[C.kcat94]*y[V.ERKstar_ERKPpase]
    v[95] = (x[C.kon95]*y[V.pERK]*y[V.ERKPpase] - x[C.koff95]*y[V.pERK_ERKPpase])
    v[96] = x[C.kcat96]*y[V.pERK_ERKPpase]

    dydt = [0]*V.len_f_vars
    
    dydt[V.E1] = -v[1] - v[81]
    dydt[V.E2] = -v[5] - v[6] - v[8] - v[82]
    dydt[V.E3] = -v[2]
    dydt[V.E4] = -v[3] - v[83]
    dydt[V.E_E1] = v[1] - v[4] - v[4] - v[5] - v[61] - v[62] - v[84]
    dydt[V.H_E3] = v[2] - v[6] - v[7] - v[61]
    dydt[V.H_E4] = v[3] - v[7] - v[8] - v[9] - v[9] - v[62] - v[85]
    dydt[V.E11] = v[4] - v[10]
    dydt[V.E12] = v[5] - v[11]
    dydt[V.E23] = v[6] - v[12]
    dydt[V.E34] = v[7] - v[13]
    dydt[V.E24] = v[8] - v[14]
    dydt[V.E44] = v[9] - v[15]
    dydt[V.E11P] = v[10] - v[16] - v[17] - v[18] - v[73] - v10d
    dydt[V.E12P] = v[11] - v[19] - v[20] - v[21] - v[74]
    dydt[V.E23P] = v[12] - v[22] - v[23] - v[24] - v[25] - v[75]
    dydt[V.E34P] = v[13] - v[26] - v[27] - v[28] - v[29] - v[76]
    dydt[V.E24P] = v[14] - v[30] - v[31] - v[32] - v[33] - v[77]
    dydt[V.E44P] = v[15] - v[34] - v[35] - v[36] - v[37] - v[78]
    dydt[V.G] = -v[16] - v[19] - v[22] - v[26] - v[30] - v[34] - v[42] - v[57] - v[65] - v[69] + v16d
    dydt[V.S] = -v[17] - v[20] - v[23] - v[27] - v[31] - v[35] - v[43] - v[66] - v[70] + v17d
    dydt[V.I] = -v[24] - v[28] - v[32] - v[36] - v[44] - v[67] - v[71]
    dydt[V.R] = -v[18] - v[21] - v[25] - v[29] - v[33] - v[37] - v[45] - v[68] - v[72] + v18d
    dydt[V.O] = -v[40] - v[54] - v[58] - v[60]
    dydt[V.A] = -v[41] - v[46] - v[55] - v[59]
    dydt[V.E11G] = v[16] - v16d
    dydt[V.E11S] = v[17] - v17d
    dydt[V.E11R] = v[18] - v18d
    dydt[V.E12G] = v[19]
    dydt[V.E12S] = v[20]
    dydt[V.E12R] = v[21]
    dydt[V.E23G] = v[22]
    dydt[V.E23S] = v[23]
    dydt[V.E23I] = v[24]
    dydt[V.E23R] = v[25]
    dydt[V.E34G] = v[26]
    dydt[V.E34S] = v[27]
    dydt[V.E34I] = v[28]
    dydt[V.E34R] = v[29]
    dydt[V.E24G] = v[30]
    dydt[V.E24S] = v[31]
    dydt[V.E24I] = v[32]
    dydt[V.E24R] = v[33]
    dydt[V.E44G] = v[34]
    dydt[V.E44S] = v[35]
    dydt[V.E44I] = v[36]
    dydt[V.E44R] = v[37]
    dydt[V.sigmaG] = v[16] + v[19] + v[22] + v[26] + v[30] + v[34] - v[40] - v[41] + v[42] + v[65] + v[69] - v16d
    dydt[V.sigmaS] = v[17] + v[20] + v[23] + v[27] + v[31] + v[35] - v[38] + v[43] + v[66] + v[70] - v17d
    dydt[V.sigmaI] = v[24] + v[28] + v[32] + v[36] + v[44] + v[67] + v[71]
    dydt[V.sigmaR] = v[18] + v[21] + v[25] + v[29] + v[33] + v[37] + v[45] - v[50] + v[68] + v[72] - v18d
    dydt[V.sigmaA] = -v[39] + v[41] + v[46] + v[59]
    dydt[V.sigmaSP] = v[38] - v[42]
    dydt[V.sigmaAP] = v[39] - v[43] - v[44] - v[45] - v[88]
    dydt[V.sigmaG_O] = v[40] - v[59]
    dydt[V.sigmaG_A] = v[41] - v[60]
    dydt[V.sigmaSP_G] = v[42]
    dydt[V.sigmaAP_S] = v[43]
    dydt[V.sigmaAP_I] = v[44]
    dydt[V.sigmaAP_R] = v[45]
    dydt[V.P3_A] = v[46] - v[57]
    dydt[V.P2] = -v[48]
    dydt[V.P3] = -v[46] + v[48]
    dydt[V.Akt] = -v[47]
    dydt[V.RsD] = -v[49]
    dydt[V.RsT] = v[49]
    dydt[V.sigmaRP] = v[50]
    dydt[V.Raf] = -v[51]
    dydt[V.Rafstar] = v[51]
    dydt[V.MEK] = -v[52]
    dydt[V.MEKstar] = v[52] - v[89] + v[90] - v[91] + v[92]
    dydt[V.ERK] = -v[89] + v[96]
    dydt[V.ERKstar] = v[92] - v[93]
    dydt[V.OP] = v[54]
    dydt[V.AP] = v[55]
    dydt[V.A_sigmaG_O] = v[59] + v[60]
    dydt[V.sigmaA_G] = v[57] - v[58]
    dydt[V.sigmaA_G_O] = v[58]
    dydt[V.sigmaO] = v[40] + v[58] + v[60]
    dydt[V.E13] = v[61] - v[63]
    dydt[V.E14] = v[62] - v[64]
    dydt[V.E13P] = v[63] - v[65] - v[66] - v[67] - v[68] - v[79]
    dydt[V.E14P] = v[64] - v[69] - v[70] - v[71] - v[72] - v[80]
    dydt[V.E13G] = v[65]
    dydt[V.E13S] = v[66]
    dydt[V.E13I] = v[67]
    dydt[V.E13R] = v[68]
    dydt[V.E14G] = v[69]
    dydt[V.E14S] = v[70]
    dydt[V.E14I] = v[71]
    dydt[V.E14R] = v[72]
    dydt[V.T] = -v[73] - v[74] - v[75] - v[76] - v[77] - v[78] - v[79] - v[80] - v[88] + v73d
    dydt[V.E11T] = v[73] - v73d
    dydt[V.E12T] = v[74]
    dydt[V.E23T] = v[75]
    dydt[V.E34T] = v[76]
    dydt[V.E24T] = v[77]
    dydt[V.E44T] = v[78]
    dydt[V.E13T] = v[79]
    dydt[V.E14T] = v[80]
    dydt[V.sigmaT] = v[73] + v[74] + v[75] + v[76] + v[77] + v[78] + v[79] + v[80] + v[88] - v73d
    dydt[V.E1_PT] = v[81] - v[86]
    dydt[V.E2_PT] = v[82]
    dydt[V.E4_PT] = v[83] - v[87]
    dydt[V.E_E1_PT] = v[84] + v[86]
    dydt[V.H_E4_PT] = v[85] + v[87]
    dydt[V.Aktstar] = v[47]
    dydt[V.sigmaAP_T] = v[88]
    dydt[V.E] = (-v[1] - v[86])/VeVc
    dydt[V.H] = (-v[2] - v[3] - v[87])/VeVc
    dydt[V.fint] = x[C.a98]*(-y[V.fint] + x[C.b98])
    dydt[V.pERK] = v[90] - v[91] + v[94] - v[95]
    dydt[V.ERK_MEKstar] = v[89] - v[90]
    dydt[V.pERK_MEKstar] = v[91] - v[92]
    dydt[V.pERK_ERKPpase] = v[95] - v[96]
    dydt[V.ERKPpase] = -v[93] + v[94] - v[95] + v[96]
    dydt[V.ERKstar_ERKPpase] = v[93] - v[94]
    
    return dydt