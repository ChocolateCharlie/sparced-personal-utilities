#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "SPARCED_x_rdata.h"
#include "SPARCED_p.h"
#include "SPARCED_tcl.h"
#include "SPARCED_dtotal_cldx_rdata.h"

namespace amici {
namespace model_SPARCED {

void dtotal_cldx_rdata_SPARCED(realtype *dtotal_cldx_rdata, const realtype *x_rdata, const realtype *p, const realtype *k, const realtype *tcl){
    dtotal_cl140_dm_TP53 = 1.0;  // dtotal_cldx_rdata[0]
    dtotal_cl139_dm_MDM2 = 1.0;  // dtotal_cldx_rdata[1]
    dtotal_cl138_dm_PPM1D = 1.0;  // dtotal_cldx_rdata[2]
    dtotal_cl137_dm_ATM = 1.0;  // dtotal_cldx_rdata[3]
    dtotal_cl136_dm_ATR = 1.0;  // dtotal_cldx_rdata[4]
    dtotal_cl135_dm_RB1 = 1.0;  // dtotal_cldx_rdata[5]
    dtotal_cl134_dm_E2F1 = 1.0;  // dtotal_cldx_rdata[6]
    dtotal_cl133_dm_E2F2 = 1.0;  // dtotal_cldx_rdata[7]
    dtotal_cl132_dm_E2F3 = 1.0;  // dtotal_cldx_rdata[8]
    dtotal_cl131_dm_CCND1 = 1.0;  // dtotal_cldx_rdata[9]
    dtotal_cl130_dm_CCND2 = 1.0;  // dtotal_cldx_rdata[10]
    dtotal_cl129_dm_CCND3 = 1.0;  // dtotal_cldx_rdata[11]
    dtotal_cl128_dm_CCNE1 = 1.0;  // dtotal_cldx_rdata[12]
    dtotal_cl127_dm_CCNE2 = 1.0;  // dtotal_cldx_rdata[13]
    dtotal_cl126_dm_SKP2 = 1.0;  // dtotal_cldx_rdata[14]
    dtotal_cl125_dm_CDC25A = 1.0;  // dtotal_cldx_rdata[15]
    dtotal_cl124_dm_CDC25B = 1.0;  // dtotal_cldx_rdata[16]
    dtotal_cl123_dm_CDC25C = 1.0;  // dtotal_cldx_rdata[17]
    dtotal_cl122_dm_CCNA2 = 1.0;  // dtotal_cldx_rdata[18]
    dtotal_cl121_dm_CDKN1B = 1.0;  // dtotal_cldx_rdata[19]
    dtotal_cl120_dm_CDH1 = 1.0;  // dtotal_cldx_rdata[20]
    dtotal_cl119_dm_CCNB1 = 1.0;  // dtotal_cldx_rdata[21]
    dtotal_cl118_dm_CDC20 = 1.0;  // dtotal_cldx_rdata[22]
    dtotal_cl117_dm_WEE1 = 1.0;  // dtotal_cldx_rdata[23]
    dtotal_cl116_dm_CHEK1 = 1.0;  // dtotal_cldx_rdata[24]
    dtotal_cl115_dm_CDKN1A = 1.0;  // dtotal_cldx_rdata[25]
    dtotal_cl114_dm_CDK1 = 1.0;  // dtotal_cldx_rdata[26]
    dtotal_cl113_dm_CDK2 = 1.0;  // dtotal_cldx_rdata[27]
    dtotal_cl112_dm_CDK4 = 1.0;  // dtotal_cldx_rdata[28]
    dtotal_cl111_dm_CDK6 = 1.0;  // dtotal_cldx_rdata[29]
    dtotal_cl110_dm_TNFSF10 = 1.0;  // dtotal_cldx_rdata[30]
    dtotal_cl109_dm_TNFRSF10A = 1.0;  // dtotal_cldx_rdata[31]
    dtotal_cl108_dm_TNFRSF10B = 1.0;  // dtotal_cldx_rdata[32]
    dtotal_cl107_dm_CFLAR = 1.0;  // dtotal_cldx_rdata[33]
    dtotal_cl106_dm_CASP8 = 1.0;  // dtotal_cldx_rdata[34]
    dtotal_cl105_dm_CASP10 = 1.0;  // dtotal_cldx_rdata[35]
    dtotal_cl104_dm_BFAR = 1.0;  // dtotal_cldx_rdata[36]
    dtotal_cl103_dm_CASP3 = 1.0;  // dtotal_cldx_rdata[37]
    dtotal_cl102_dm_CASP7 = 1.0;  // dtotal_cldx_rdata[38]
    dtotal_cl101_dm_CASP6 = 1.0;  // dtotal_cldx_rdata[39]
    dtotal_cl100_dm_XIAP = 1.0;  // dtotal_cldx_rdata[40]
    dtotal_cl99_dm_PARP1 = 1.0;  // dtotal_cldx_rdata[41]
    dtotal_cl98_dm_BID = 1.0;  // dtotal_cldx_rdata[42]
    dtotal_cl97_dm_BCL2 = 1.0;  // dtotal_cldx_rdata[43]
    dtotal_cl96_dm_BCL2L1 = 1.0;  // dtotal_cldx_rdata[44]
    dtotal_cl95_dm_MCL1 = 1.0;  // dtotal_cldx_rdata[45]
    dtotal_cl94_dm_BAX = 1.0;  // dtotal_cldx_rdata[46]
    dtotal_cl93_dm_CYCS = 1.0;  // dtotal_cldx_rdata[47]
    dtotal_cl92_dm_DIABLO = 1.0;  // dtotal_cldx_rdata[48]
    dtotal_cl91_dm_APAF1 = 1.0;  // dtotal_cldx_rdata[49]
    dtotal_cl90_dm_CASP9 = 1.0;  // dtotal_cldx_rdata[50]
    dtotal_cl89_dm_BAD = 1.0;  // dtotal_cldx_rdata[51]
    dtotal_cl88_dm_BBC3 = 1.0;  // dtotal_cldx_rdata[52]
    dtotal_cl87_dm_PMAIP1 = 1.0;  // dtotal_cldx_rdata[53]
    dtotal_cl86_dm_BCL2L11 = 1.0;  // dtotal_cldx_rdata[54]
    dtotal_cl85_dm_EGF = 1.0;  // dtotal_cldx_rdata[55]
    dtotal_cl84_dm_NRG1 = 1.0;  // dtotal_cldx_rdata[56]
    dtotal_cl83_dm_EGFR = 1.0;  // dtotal_cldx_rdata[57]
    dtotal_cl82_dm_ERBB2 = 1.0;  // dtotal_cldx_rdata[58]
    dtotal_cl81_dm_ERBB3 = 1.0;  // dtotal_cldx_rdata[59]
    dtotal_cl80_dm_ERBB4 = 1.0;  // dtotal_cldx_rdata[60]
    dtotal_cl79_dm_EGFRvIII = 1.0;  // dtotal_cldx_rdata[61]
    dtotal_cl78_dm_MET = 1.0;  // dtotal_cldx_rdata[62]
    dtotal_cl77_dm_HGF = 1.0;  // dtotal_cldx_rdata[63]
    dtotal_cl76_dm_PDGFRA = 1.0;  // dtotal_cldx_rdata[64]
    dtotal_cl75_dm_PDGFRB = 1.0;  // dtotal_cldx_rdata[65]
    dtotal_cl74_dm_PDGFB = 1.0;  // dtotal_cldx_rdata[66]
    dtotal_cl73_dm_SPRY2 = 1.0;  // dtotal_cldx_rdata[67]
    dtotal_cl72_dm_CBL = 1.0;  // dtotal_cldx_rdata[68]
    dtotal_cl71_dm_GRB2 = 1.0;  // dtotal_cldx_rdata[69]
    dtotal_cl70_dm_PLCG1 = 1.0;  // dtotal_cldx_rdata[70]
    dtotal_cl69_dm_PLCG2 = 1.0;  // dtotal_cldx_rdata[71]
    dtotal_cl68_dm_PIK3CA = 1.0;  // dtotal_cldx_rdata[72]
    dtotal_cl67_dm_PIK3CB = 1.0;  // dtotal_cldx_rdata[73]
    dtotal_cl66_dm_PIK3CG = 1.0;  // dtotal_cldx_rdata[74]
    dtotal_cl65_dm_PIK3CD = 1.0;  // dtotal_cldx_rdata[75]
    dtotal_cl64_dm_PIK3R1 = 1.0;  // dtotal_cldx_rdata[76]
    dtotal_cl63_dm_PIK3R2 = 1.0;  // dtotal_cldx_rdata[77]
    dtotal_cl62_dm_PIK3R3 = 1.0;  // dtotal_cldx_rdata[78]
    dtotal_cl61_dm_PIK3R4 = 1.0;  // dtotal_cldx_rdata[79]
    dtotal_cl60_dm_PIK3C2A = 1.0;  // dtotal_cldx_rdata[80]
    dtotal_cl59_dm_RASGRP1 = 1.0;  // dtotal_cldx_rdata[81]
    dtotal_cl58_dm_RASGRP3 = 1.0;  // dtotal_cldx_rdata[82]
    dtotal_cl57_dm_NRAS = 1.0;  // dtotal_cldx_rdata[83]
    dtotal_cl56_dm_KRAS = 1.0;  // dtotal_cldx_rdata[84]
    dtotal_cl55_dm_HRAS = 1.0;  // dtotal_cldx_rdata[85]
    dtotal_cl54_dm_NF1 = 1.0;  // dtotal_cldx_rdata[86]
    dtotal_cl53_dm_RAF1 = 1.0;  // dtotal_cldx_rdata[87]
    dtotal_cl52_dm_BRAF = 1.0;  // dtotal_cldx_rdata[88]
    dtotal_cl51_dm_MAP2K1 = 1.0;  // dtotal_cldx_rdata[89]
    dtotal_cl50_dm_MAP2K2 = 1.0;  // dtotal_cldx_rdata[90]
    dtotal_cl49_dm_DUSP6 = 1.0;  // dtotal_cldx_rdata[91]
    dtotal_cl48_dm_RPS6KA1 = 1.0;  // dtotal_cldx_rdata[92]
    dtotal_cl47_dm_RPS6KA2 = 1.0;  // dtotal_cldx_rdata[93]
    dtotal_cl46_dm_RPS6KA3 = 1.0;  // dtotal_cldx_rdata[94]
    dtotal_cl45_dm_RPS6KA4 = 1.0;  // dtotal_cldx_rdata[95]
    dtotal_cl44_dm_DUSP1 = 1.0;  // dtotal_cldx_rdata[96]
    dtotal_cl43_dm_FOS = 1.0;  // dtotal_cldx_rdata[97]
    dtotal_cl42_dm_JUN = 1.0;  // dtotal_cldx_rdata[98]
    dtotal_cl41_dm_MYC = 1.0;  // dtotal_cldx_rdata[99]
    dtotal_cl40_dm_CTNNB1 = 1.0;  // dtotal_cldx_rdata[100]
    dtotal_cl39_dm_PTEN = 1.0;  // dtotal_cldx_rdata[101]
    dtotal_cl38_dm_AKT1 = 1.0;  // dtotal_cldx_rdata[102]
    dtotal_cl37_dm_AKT2 = 1.0;  // dtotal_cldx_rdata[103]
    dtotal_cl36_dm_PDPK1 = 1.0;  // dtotal_cldx_rdata[104]
    dtotal_cl35_dm_RICTOR = 1.0;  // dtotal_cldx_rdata[105]
    dtotal_cl34_dm_MTOR = 1.0;  // dtotal_cldx_rdata[106]
    dtotal_cl33_dm_GSK3B = 1.0;  // dtotal_cldx_rdata[107]
    dtotal_cl32_dm_TSC1 = 1.0;  // dtotal_cldx_rdata[108]
    dtotal_cl31_dm_TSC2 = 1.0;  // dtotal_cldx_rdata[109]
    dtotal_cl30_dm_PRKCA = 1.0;  // dtotal_cldx_rdata[110]
    dtotal_cl29_dm_PRKCB = 1.0;  // dtotal_cldx_rdata[111]
    dtotal_cl28_dm_PRKCG = 1.0;  // dtotal_cldx_rdata[112]
    dtotal_cl27_dm_PRKCD = 1.0;  // dtotal_cldx_rdata[113]
    dtotal_cl26_dm_PEBP1 = 1.0;  // dtotal_cldx_rdata[114]
    dtotal_cl25_dm_MAPK1 = 1.0;  // dtotal_cldx_rdata[115]
    dtotal_cl24_dm_MAPK3 = 1.0;  // dtotal_cldx_rdata[116]
    dtotal_cl23_dm_FOXO3 = 1.0;  // dtotal_cldx_rdata[117]
    dtotal_cl22_dm_RHEB = 1.0;  // dtotal_cldx_rdata[118]
    dtotal_cl21_dm_RPTOR = 1.0;  // dtotal_cldx_rdata[119]
    dtotal_cl20_dm_RPS6KB1 = 1.0;  // dtotal_cldx_rdata[120]
    dtotal_cl19_dm_RPS6KB2 = 1.0;  // dtotal_cldx_rdata[121]
    dtotal_cl18_dm_EIF4EBP1 = 1.0;  // dtotal_cldx_rdata[122]
    dtotal_cl17_dm_SOS1 = 1.0;  // dtotal_cldx_rdata[123]
    dtotal_cl16_dm_CDKN2A = 1.0;  // dtotal_cldx_rdata[124]
    dtotal_cl15_dm_MDM4 = 1.0;  // dtotal_cldx_rdata[125]
    dtotal_cl14_dm_FGFR1 = 1.0;  // dtotal_cldx_rdata[126]
    dtotal_cl13_dm_FGFR2 = 1.0;  // dtotal_cldx_rdata[127]
    dtotal_cl12_dm_FGF1 = 1.0;  // dtotal_cldx_rdata[128]
    dtotal_cl11_dm_FGF2 = 1.0;  // dtotal_cldx_rdata[129]
    dtotal_cl10_dm_EIF4E = 1.0;  // dtotal_cldx_rdata[130]
    dtotal_cl9_dm_IRS1 = 1.0;  // dtotal_cldx_rdata[131]
    dtotal_cl8_dm_IRS2 = 1.0;  // dtotal_cldx_rdata[132]
    dtotal_cl7_dm_IGF1 = 1.0;  // dtotal_cldx_rdata[133]
    dtotal_cl6_dm_IGF2 = 1.0;  // dtotal_cldx_rdata[134]
    dtotal_cl5_dm_IGF1R = 1.0;  // dtotal_cldx_rdata[135]
    dtotal_cl4_dm_MSH6 = 1.0;  // dtotal_cldx_rdata[136]
    dtotal_cl3_dm_BRCA2 = 1.0;  // dtotal_cldx_rdata[137]
    dtotal_cl2_dm_MGMT = 1.0;  // dtotal_cldx_rdata[138]
    dtotal_cl1_dm_INSR = 1.0;  // dtotal_cldx_rdata[139]
    dtotal_cl0_dm_INS = 1.0;  // dtotal_cldx_rdata[140]
}

} // namespace model_SPARCED
} // namespace amici
