const F_P = [
    "VmaxPY"
    "KmPY"
    "kdeg"
    "kf47"
    "Vmaxr47"
    "Kmf47"
    "Kmr47"
    "kf48"
    "Kmf48"
    "Kmr48"
    "PTEN"
    "kf49"
    "kr49"
    "Kmf49"
    "Kmr49"
    "Kmr49b"
    "kr49b"
    "kf51"
    "Vmaxr51"
    "Kmf51"
    "Kmr51"
    "Kmrb51"
    "kf52"
    "Vmaxr52"
    "Kmf52"
    "Kmr52"
    "kf54"
    "Vmaxr54"
    "Kmf54"
    "Kmr54"
    "kf55"
    "Vmaxr55"
    "Kmf55"
    "Kmr55"
    "kf38"
    "kf39"
    "kf50"
    "a98"
    "b98"
    "koff46"
    "EGF_off"
    "HRGoff_3"
    "HRGoff_4"
    "koff4"
    "koff5"
    "koff6"
    "koff7"
    "koff8"
    "koff9"
    "koff61"
    "koff62"
    "koff16"
    "koff17"
    "koff18"
    "koff19"
    "koff20"
    "koff21"
    "koff22"
    "koff23"
    "koff24"
    "koff25"
    "koff26"
    "koff27"
    "koff28"
    "koff29"
    "koff30"
    "koff31"
    "koff32"
    "koff33"
    "koff34"
    "koff35"
    "koff36"
    "koff37"
    "koff65"
    "koff66"
    "koff67"
    "koff68"
    "koff69"
    "koff70"
    "koff71"
    "koff72"
    "koff40"
    "koff41"
    "koff42"
    "koff43"
    "koff44"
    "koff45"
    "koff57"
    "koff58"
    "koff59"
    "koff60"
    "kPTP10"
    "kPTP11"
    "kPTP12"
    "kPTP13"
    "kPTP14"
    "kPTP15"
    "kPTP63"
    "kPTP64"
    "koff73"
    "koff74"
    "koff75"
    "koff76"
    "koff77"
    "koff78"
    "koff79"
    "koff80"
    "kPTP38"
    "kPTP39"
    "koff88"
    "kPTP50"
    "kf81"
    "Vmaxr81"
    "Kmf81"
    "Kmr81"
    "kf82"
    "Vmaxr82"
    "Kmf82"
    "Kmr82"
    "kf83"
    "Vmaxr83"
    "Kmf83"
    "Kmr83"
    "kf84"
    "Vmaxr84"
    "Kmf84"
    "Kmr84"
    "kf85"
    "Vmaxr85"
    "Kmf85"
    "Kmr85"
    "kcon49"
    "kon1"
    "kon86"
    "kon2"
    "kon3"
    "kon87"
    "kon4"
    "kon5"
    "kon6"
    "kon7"
    "kon8"
    "kon9"
    "kon61"
    "kon62"
    "kf10"
    "kf11"
    "kf12"
    "kf13"
    "kf14"
    "kf15"
    "kf63"
    "kf64"
    "kon16"
    "kon17"
    "kon18"
    "kon73"
    "kon19"
    "kon20"
    "kon21"
    "kon74"
    "kon22"
    "kon23"
    "kon24"
    "kon25"
    "kon75"
    "kon26"
    "kon27"
    "kon28"
    "kon29"
    "kon76"
    "kon30"
    "kon31"
    "kon32"
    "kon33"
    "kon77"
    "kon34"
    "kon35"
    "kon36"
    "kon37"
    "kon78"
    "kon65"
    "kon66"
    "kon67"
    "kon68"
    "kon79"
    "kon69"
    "kon70"
    "kon71"
    "kon72"
    "kon80"
    "kon40"
    "kon41"
    "kon42"
    "kon43"
    "kon44"
    "kon45"
    "kon88"
    "kon46"
    "kon57"
    "kon58"
    "kon59"
    "kon60"
    "kon89"
    "koff89"
    "kcat90"
    "kon91"
    "koff91"
    "kcat92"
    "kon93"
    "koff93"
    "kcat94"
    "kon95"
    "koff95"
    "kcat96"
];

for (index,value) in enumerate(F_P)
    eval(Meta.parse("const $value=$index"));
end

function f_params()
    p::Vector{Float64} = zeros(length(F_P));

    p[VmaxPY] = 223.8776;
    p[KmPY] = 486.1398;
    p[kdeg] = 0.0259;
    p[kf47] = 24.6048;
    p[Vmaxr47] = 590.5058;
    p[Kmf47] = 698.6022;
    p[Kmr47] = 483.8622;
    p[kf48] = 16.833;
    p[Kmf48] = 715.5688;
    p[Kmr48] = 324.9294;
    p[PTEN] = 693.5786;
    p[kf49] = 44.3501;
    p[kr49] = 552.6746;
    p[Kmf49] = 343.2483;
    p[Kmr49] = 753.1667;
    p[Kmr49b] = 381.2208;
    p[kr49b] = 640.8212;
    p[kf51] = 3.6515;
    p[Vmaxr51] = 16.737;
    p[Kmf51] = 599.7076;
    p[Kmr51] = 346.4779;
    p[Kmrb51] = 988.4496;
    p[kf52] = 0.7745;
    p[Vmaxr52] = 199.2773;
    p[Kmf52] = 545.4408;
    p[Kmr52] = 675.2994;
    p[kf54] = 0.0538;
    p[Vmaxr54] = 588.2671;
    p[Kmf54] = 457.9645;
    p[Kmr54] = 336.183;
    p[kf55] = 0.2256;
    p[Vmaxr55] = 646.9003;
    p[Kmf55] = 460.9446;
    p[Kmr55] = 643.9247;
    p[kf38] = 279.9929;
    p[kf39] = 385.7428;
    p[kf50] = 389.1061;
    p[a98] = 0.0849;
    p[b98] = 0.1833;
    p[koff46] = 0.5194;
    p[EGF_off] = 0.0175;
    p[HRGoff_3] = 9.0E-4;
    p[HRGoff_4] = 0.0973;
    p[koff4] = 0.1717;
    p[koff5] = 4.3985;
    p[koff6] = 2.6619;
    p[koff7] = 8.0557;
    p[koff8] = 9.1034;
    p[koff9] = 5.5425;
    p[koff57] = 0.4526;
    p[koff16] = 0.5737;
    p[koff17] = 4.6874;
    p[koff18] = 2.2768;
    p[koff19] = 2.3361;
    p[koff20] = 0.6761;
    p[koff21] = 4.7291;
    p[koff22] = 3.6962;
    p[koff23] = 2.3619;
    p[koff24] = 4.4226;
    p[koff25] = 2.225;
    p[koff26] = 0.0103;
    p[koff27] = 1.8922;
    p[koff28] = 4.6432;
    p[koff29] = 2.0148;
    p[koff30] = 4.9936;
    p[koff31] = 1.2204;
    p[koff32] = 3.8752;
    p[koff33] = 1.2817;
    p[koff34] = 3.2036;
    p[koff35] = 1.8696;
    p[koff36] = 1.2567;
    p[koff37] = 0.4059;
    p[koff65] = 0.1185;
    p[koff66] = 2.2988;
    p[koff67] = 1.6142;
    p[koff40] = 3.1051;
    p[koff41] = 7.0487;
    p[koff42] = 3.5195;
    p[koff43] = 0.5441;
    p[koff44] = 0.4265;
    p[koff45] = 3.9967;
    p[koff58] = 6.3059;
    p[koff59] = 9.172;
    p[koff68] = 2.8871;
    p[kPTP10] = 29.8531;
    p[kPTP11] = 78.204;
    p[kPTP12] = 11.4211;
    p[kPTP13] = 55.2104;
    p[kPTP14] = 57.7506;
    p[kPTP15] = 60.2628;
    p[kPTP63] = 7.4766;
    p[kPTP64] = 48.6335;
    p[koff73] = 3.0048;
    p[koff74] = 1.2496;
    p[koff75] = 1.4323;
    p[koff76] = 2.1542;
    p[koff77] = 1.2237;
    p[koff78] = 0.2007;
    p[koff79] = 1.1852;
    p[koff80] = 2.9373;
    p[kPTP38] = 83.4465;
    p[kPTP39] = 79.6132;
    p[koff88] = 3.9255;
    p[kPTP50] = 96.5716;
    p[kf81] = 1.361;
    p[Vmaxr81] = 242.6034;
    p[Kmf81] = 485.2626;
    p[Kmr81] = 323.4012;
    p[kf82] = 6.9987;
    p[Vmaxr82] = 398.1931;
    p[Kmf82] = 781.4374;
    p[Kmr82] = 595.8395;
    p[kf83] = 1.7628;
    p[Vmaxr83] = 534.0531;
    p[Kmf83] = 609.4766;
    p[Kmr83] = 653.5184;
    p[kf84] = 4.6894;
    p[Vmaxr84] = 634.1626;
    p[Kmf84] = 622.3847;
    p[Kmr84] = 258.4637;
    p[kf85] = 6.7591;
    p[Vmaxr85] = 369.2261;
    p[Kmf85] = 179.6486;
    p[Kmr85] = 290.7667;
    p[kcon49] = 9.9783;
    p[kon1] = 1.0086E-4;
    p[kon86] = 0.0038;
    p[kon2] = 0.0059;
    p[kon3] = 0.0334;
    p[kon87] = 8.0E-4;
    p[kon4] = 0.5005;
    p[kon5] = 2.5427;
    p[kon6] = 0.2283;
    p[kon7] = 1.0606;
    p[kon8] = 1.0259;
    p[kon9] = 2.2868;
    p[kon57] = 0.0039;
    p[kf10] = 0.6496;
    p[kf11] = 0.3721;
    p[kf12] = 1.8012;
    p[kf13] = 0.8875;
    p[kf14] = 6.1726;
    p[kf15] = 1.3565;
    p[kf63] = 0.9297;
    p[kf64] = 1.2083;
    p[kon16] = 0.0097;
    p[kon17] = 0.0049;
    p[kon18] = 0.0117;
    p[kon73] = 0.0116;
    p[kon19] = 0.0896;
    p[kon20] = 0.0478;
    p[kon21] = 0.0114;
    p[kon74] = 0.0133;
    p[kon22] = 7.0E-4;
    p[kon23] = 0.0138;
    p[kon24] = 0.005;
    p[kon25] = 0.0995;
    p[kon75] = 0.0137;
    p[kon26] = 0.0355;
    p[kon27] = 0.0201;
    p[kon28] = 0.0074;
    p[kon29] = 0.0346;
    p[kon76] = 0.0053;
    p[kon30] = 0.002;
    p[kon31] = 0.0032;
    p[kon32] = 9.0E-4;
    p[kon33] = 0.0335;
    p[kon77] = 0.0101;
    p[kon34] = 1.0E-4;
    p[kon35] = 0.0602;
    p[kon36] = 0.0043;
    p[kon37] = 0.0791;
    p[kon78] = 0.0076;
    p[kon79] = 0.0078;
    p[kon65] = 0.0123;
    p[kon66] = 1.9264E-4;
    p[kon67] = 6.6667E-5;
    p[kon80] = 2.0E-4;
    p[kon40] = 0.0191;
    p[kon41] = 0.0051;
    p[kon42] = 0.0023;
    p[kon43] = 0.0127;
    p[kon44] = 0.0122;
    p[kon45] = 0.0028;
    p[kon88] = 0.0108;
    p[kon46] = 0.0148;
    p[kon58] = 0.0215;
    p[kon59] = 0.0077;
    p[kon60] = 1.1994E-4;
    p[koff60] = 4.9981;
    p[koff61] = 5.229;
    p[kon61] = 0.8048;
    p[kon62] = 1.782;
    p[koff62] = 5.5142;
    p[kon68] = 0.0045;
    p[kon69] = 0.0084;
    p[koff69] = 3.97;
    p[kon70] = 0.0116;
    p[koff70] = 2.6069;
    p[kon71] = 0.0078;
    p[koff71] = 2.2988;
    p[kon72] = 0.0355;
    p[koff72] = 0.907;
    p[kon89] = 0.1997;
    p[koff89] = 99.9637;
    p[kcat90] = 20.0037;
    p[kon91] = 0.1966;
    p[koff91] = 99.9983;
    p[kcat92] = 0.2004;
    p[kon93] = 0.2003;
    p[koff93] = 100.0037;
    p[kcat94] = 0.9966;
    p[kon95] = 0.1993;
    p[koff95] = 100.0023;
    p[kcat96] = 19.9851;

    return p
end