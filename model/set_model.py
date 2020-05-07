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


def f_params():
    
    x = [0]*C.len_f_params

    x[C.VmaxPY] = 223.8776
    x[C.KmPY] = 486.1398
    x[C.kdeg] = 0.0259
    x[C.kf47] = 24.6048
    x[C.Vmaxr47] = 590.5058
    x[C.Kmf47] = 698.6022
    x[C.Kmr47] = 483.8622
    x[C.kf48] = 16.833
    x[C.Kmf48] = 715.5688
    x[C.Kmr48] = 324.9294
    x[C.PTEN] = 693.5786
    x[C.kf49] = 44.3501
    x[C.kr49] = 552.6746
    x[C.Kmf49] = 343.2483
    x[C.Kmr49] = 753.1667
    x[C.Kmr49b] = 381.2208
    x[C.kr49b] = 640.8212
    x[C.kf51] = 3.6515
    x[C.Vmaxr51] = 16.737
    x[C.Kmf51] = 599.7076
    x[C.Kmr51] = 346.4779
    x[C.Kmrb51] = 988.4496
    x[C.kf52] = 0.7745
    x[C.Vmaxr52] = 199.2773
    x[C.Kmf52] = 545.4408
    x[C.Kmr52] = 675.2994
    x[C.kf54] = 0.0538
    x[C.Vmaxr54] = 588.2671
    x[C.Kmf54] = 457.9645
    x[C.Kmr54] = 336.183
    x[C.kf55] = 0.2256
    x[C.Vmaxr55] = 646.9003
    x[C.Kmf55] = 460.9446
    x[C.Kmr55] = 643.9247
    x[C.kf38] = 279.9929
    x[C.kf39] = 385.7428
    x[C.kf50] = 389.1061
    x[C.a98] = 0.0849
    x[C.b98] = 0.1833
    x[C.koff46] = 0.5194
    x[C.EGF_off] = 0.0175
    x[C.HRGoff_3] = 9.0E-4
    x[C.HRGoff_4] = 0.0973
    x[C.koff4] = 0.1717
    x[C.koff5] = 4.3985
    x[C.koff6] = 2.6619
    x[C.koff7] = 8.0557
    x[C.koff8] = 9.1034
    x[C.koff9] = 5.5425
    x[C.koff57] = 0.4526
    x[C.koff16] = 0.5737
    x[C.koff17] = 4.6874
    x[C.koff18] = 2.2768
    x[C.koff19] = 2.3361
    x[C.koff20] = 0.6761
    x[C.koff21] = 4.7291
    x[C.koff22] = 3.6962
    x[C.koff23] = 2.3619
    x[C.koff24] = 4.4226
    x[C.koff25] = 2.225
    x[C.koff26] = 0.0103
    x[C.koff27] = 1.8922
    x[C.koff28] = 4.6432
    x[C.koff29] = 2.0148
    x[C.koff30] = 4.9936
    x[C.koff31] = 1.2204
    x[C.koff32] = 3.8752
    x[C.koff33] = 1.2817
    x[C.koff34] = 3.2036
    x[C.koff35] = 1.8696
    x[C.koff36] = 1.2567
    x[C.koff37] = 0.4059
    x[C.koff65] = 0.1185
    x[C.koff66] = 2.2988
    x[C.koff67] = 1.6142
    x[C.koff40] = 3.1051
    x[C.koff41] = 7.0487
    x[C.koff42] = 3.5195
    x[C.koff43] = 0.5441
    x[C.koff44] = 0.4265
    x[C.koff45] = 3.9967
    x[C.koff58] = 6.3059
    x[C.koff59] = 9.172
    x[C.koff68] = 2.8871
    x[C.kPTP10] = 29.8531
    x[C.kPTP11] = 78.204
    x[C.kPTP12] = 11.4211
    x[C.kPTP13] = 55.2104
    x[C.kPTP14] = 57.7506
    x[C.kPTP15] = 60.2628
    x[C.kPTP63] = 7.4766
    x[C.kPTP64] = 48.6335
    x[C.koff73] = 3.0048
    x[C.koff74] = 1.2496
    x[C.koff75] = 1.4323
    x[C.koff76] = 2.1542
    x[C.koff77] = 1.2237
    x[C.koff78] = 0.2007
    x[C.koff79] = 1.1852
    x[C.koff80] = 2.9373
    x[C.kPTP38] = 83.4465
    x[C.kPTP39] = 79.6132
    x[C.koff88] = 3.9255
    x[C.kPTP50] = 96.5716
    x[C.kf81] = 1.361
    x[C.Vmaxr81] = 242.6034
    x[C.Kmf81] = 485.2626
    x[C.Kmr81] = 323.4012
    x[C.kf82] = 6.9987
    x[C.Vmaxr82] = 398.1931
    x[C.Kmf82] = 781.4374
    x[C.Kmr82] = 595.8395
    x[C.kf83] = 1.7628
    x[C.Vmaxr83] = 534.0531
    x[C.Kmf83] = 609.4766
    x[C.Kmr83] = 653.5184
    x[C.kf84] = 4.6894
    x[C.Vmaxr84] = 634.1626
    x[C.Kmf84] = 622.3847
    x[C.Kmr84] = 258.4637
    x[C.kf85] = 6.7591
    x[C.Vmaxr85] = 369.2261
    x[C.Kmf85] = 179.6486
    x[C.Kmr85] = 290.7667
    x[C.kcon49] = 9.9783
    x[C.kon1] = 1.0086E-4
    x[C.kon86] = 0.0038
    x[C.kon2] = 0.0059
    x[C.kon3] = 0.0334
    x[C.kon87] = 8.0E-4
    x[C.kon4] = 0.5005
    x[C.kon5] = 2.5427
    x[C.kon6] = 0.2283
    x[C.kon7] = 1.0606
    x[C.kon8] = 1.0259
    x[C.kon9] = 2.2868
    x[C.kon57] = 0.0039
    x[C.kf10] = 0.6496
    x[C.kf11] = 0.3721
    x[C.kf12] = 1.8012
    x[C.kf13] = 0.8875
    x[C.kf14] = 6.1726
    x[C.kf15] = 1.3565
    x[C.kf63] = 0.9297
    x[C.kf64] = 1.2083
    x[C.kon16] = 0.0097
    x[C.kon17] = 0.0049
    x[C.kon18] = 0.0117
    x[C.kon73] = 0.0116
    x[C.kon19] = 0.0896
    x[C.kon20] = 0.0478
    x[C.kon21] = 0.0114
    x[C.kon74] = 0.0133
    x[C.kon22] = 7.0E-4
    x[C.kon23] = 0.0138
    x[C.kon24] = 0.005
    x[C.kon25] = 0.0995
    x[C.kon75] = 0.0137
    x[C.kon26] = 0.0355
    x[C.kon27] = 0.0201
    x[C.kon28] = 0.0074
    x[C.kon29] = 0.0346
    x[C.kon76] = 0.0053
    x[C.kon30] = 0.002
    x[C.kon31] = 0.0032
    x[C.kon32] = 9.0E-4
    x[C.kon33] = 0.0335
    x[C.kon77] = 0.0101
    x[C.kon34] = 1.0E-4
    x[C.kon35] = 0.0602
    x[C.kon36] = 0.0043
    x[C.kon37] = 0.0791
    x[C.kon78] = 0.0076
    x[C.kon79] = 0.0078
    x[C.kon65] = 0.0123
    x[C.kon66] = 1.9264E-4
    x[C.kon67] = 6.6667E-5
    x[C.kon80] = 2.0E-4
    x[C.kon40] = 0.0191
    x[C.kon41] = 0.0051
    x[C.kon42] = 0.0023
    x[C.kon43] = 0.0127
    x[C.kon44] = 0.0122
    x[C.kon45] = 0.0028
    x[C.kon88] = 0.0108
    x[C.kon46] = 0.0148
    x[C.kon58] = 0.0215
    x[C.kon59] = 0.0077
    x[C.kon60] = 1.1994E-4
    x[C.koff60] = 4.9981
    x[C.koff61] = 5.229
    x[C.kon61] = 0.8048
    x[C.kon62] = 1.782
    x[C.koff62] = 5.5142
    x[C.kon68] = 0.0045
    x[C.kon69] = 0.0084
    x[C.koff69] = 3.97
    x[C.kon70] = 0.0116
    x[C.koff70] = 2.6069
    x[C.kon71] = 0.0078
    x[C.koff71] = 2.2988
    x[C.kon72] = 0.0355
    x[C.koff72] = 0.907
    x[C.kon89] = 0.1997
    x[C.koff89] = 99.9637
    x[C.kcat90] = 20.0037
    x[C.kon91] = 0.1966
    x[C.koff91] = 99.9983
    x[C.kcat92] = 0.2004
    x[C.kon93] = 0.2003
    x[C.koff93] = 100.0037
    x[C.kcat94] = 0.9966
    x[C.kon95] = 0.1993
    x[C.koff95] = 100.0023
    x[C.kcat96] = 19.9851

    return x


def initial_values():
    
    y0 = [0]*V.len_f_vars

    y0[V.E1] = 274.0
    y0[V.E2] = 158.0
    y0[V.E3] = 294.0
    y0[V.E4] = 399.0
    y0[V.G] = 82.4
    y0[V.S] = 11.5
    y0[V.I] = 46.4
    y0[V.R] = 93.6
    y0[V.T] = 500.0
    y0[V.O] = 82.3
    y0[V.A] = 43.1
    y0[V.P2] = 197
    y0[V.Akt] = 444.0
    y0[V.RsD] = 95.7
    y0[V.Raf] = 743.0
    y0[V.MEK] = 772.0
    y0[V.ERK] = 750.0
    y0[V.ERKPpase] = 35.0

    return y0