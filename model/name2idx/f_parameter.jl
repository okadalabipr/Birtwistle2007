module C

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

const len_f_params = length(F_P);

end  # module