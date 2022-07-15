#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>
#include <algorithm>

#include "SPARCED_p.h"
#include "SPARCED_y.h"
#include "SPARCED_sigmay.h"
#include "SPARCED_my.h"

namespace amici {
namespace model_SPARCED {

void dJydsigma_SPARCED(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_yRibosome - 1.0*std::pow(-myRibosome + yRibosome, 2)/std::pow(sigma_yRibosome, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_yp53inac - 1.0*std::pow(-myp53inac + yp53inac, 2)/std::pow(sigma_yp53inac, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_yp53ac - 1.0*std::pow(-myp53ac + yp53ac, 2)/std::pow(sigma_yp53ac, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_yMDM2 - 1.0*std::pow(-myMDM2 + yMDM2, 2)/std::pow(sigma_yMDM2, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_yWip1 - 1.0*std::pow(-myWip1 + yWip1, 2)/std::pow(sigma_yWip1, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_yATMP - 1.0*std::pow(-myATMP + yATMP, 2)/std::pow(sigma_yATMP, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_yATRac - 1.0*std::pow(-myATRac + yATRac, 2)/std::pow(sigma_yATRac, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_yMDM2product1 - 1.0*std::pow(-myMDM2product1 + yMDM2product1, 2)/std::pow(sigma_yMDM2product1, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_yMDM2product2 - 1.0*std::pow(-myMDM2product2 + yMDM2product2, 2)/std::pow(sigma_yMDM2product2, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_yMDM2product3 - 1.0*std::pow(-myMDM2product3 + yMDM2product3, 2)/std::pow(sigma_yMDM2product3, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_yMDM2product4 - 1.0*std::pow(-myMDM2product4 + yMDM2product4, 2)/std::pow(sigma_yMDM2product4, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_yMDM2product5 - 1.0*std::pow(-myMDM2product5 + yMDM2product5, 2)/std::pow(sigma_yMDM2product5, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_yMDM2product6 - 1.0*std::pow(-myMDM2product6 + yMDM2product6, 2)/std::pow(sigma_yMDM2product6, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_yMDM2product7 - 1.0*std::pow(-myMDM2product7 + yMDM2product7, 2)/std::pow(sigma_yMDM2product7, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_yMDM2product8 - 1.0*std::pow(-myMDM2product8 + yMDM2product8, 2)/std::pow(sigma_yMDM2product8, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_yMDM2product9 - 1.0*std::pow(-myMDM2product9 + yMDM2product9, 2)/std::pow(sigma_yMDM2product9, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_yMDM2pro - 1.0*std::pow(-myMDM2pro + yMDM2pro, 2)/std::pow(sigma_yMDM2pro, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_yWip1product1 - 1.0*std::pow(-myWip1product1 + yWip1product1, 2)/std::pow(sigma_yWip1product1, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_yWip1product2 - 1.0*std::pow(-myWip1product2 + yWip1product2, 2)/std::pow(sigma_yWip1product2, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_yWip1product3 - 1.0*std::pow(-myWip1product3 + yWip1product3, 2)/std::pow(sigma_yWip1product3, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_yWip1product4 - 1.0*std::pow(-myWip1product4 + yWip1product4, 2)/std::pow(sigma_yWip1product4, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_yWip1product5 - 1.0*std::pow(-myWip1product5 + yWip1product5, 2)/std::pow(sigma_yWip1product5, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_yWip1product6 - 1.0*std::pow(-myWip1product6 + yWip1product6, 2)/std::pow(sigma_yWip1product6, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_yWip1product7 - 1.0*std::pow(-myWip1product7 + yWip1product7, 2)/std::pow(sigma_yWip1product7, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_yWip1product8 - 1.0*std::pow(-myWip1product8 + yWip1product8, 2)/std::pow(sigma_yWip1product8, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_yWip1product9 - 1.0*std::pow(-myWip1product9 + yWip1product9, 2)/std::pow(sigma_yWip1product9, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_yWip1pro - 1.0*std::pow(-myWip1pro + yWip1pro, 2)/std::pow(sigma_yWip1pro, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_yBRCA2 - 1.0*std::pow(-myBRCA2 + yBRCA2, 2)/std::pow(sigma_yBRCA2, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_yMSH6 - 1.0*std::pow(-myMSH6 + yMSH6, 2)/std::pow(sigma_yMSH6, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_yMGMT - 1.0*std::pow(-myMGMT + yMGMT, 2)/std::pow(sigma_yMGMT, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_ydamageDSB - 1.0*std::pow(-mydamageDSB + ydamageDSB, 2)/std::pow(sigma_ydamageDSB, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_ydamageSSB - 1.0*std::pow(-mydamageSSB + ydamageSSB, 2)/std::pow(sigma_ydamageSSB, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_yppAKT_MDM2 - 1.0*std::pow(-myppAKT_MDM2 + yppAKT_MDM2, 2)/std::pow(sigma_yppAKT_MDM2, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_ypMDM2 - 1.0*std::pow(-mypMDM2 + ypMDM2, 2)/std::pow(sigma_ypMDM2, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_yARF - 1.0*std::pow(-myARF + yARF, 2)/std::pow(sigma_yARF, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_yMDM4 - 1.0*std::pow(-myMDM4 + yMDM4, 2)/std::pow(sigma_yMDM4, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_yp53ac_MDM4 - 1.0*std::pow(-myp53ac_MDM4 + yp53ac_MDM4, 2)/std::pow(sigma_yp53ac_MDM4, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_yATMinac - 1.0*std::pow(-myATMinac + yATMinac, 2)/std::pow(sigma_yATMinac, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_yATRinac - 1.0*std::pow(-myATRinac + yATRinac, 2)/std::pow(sigma_yATRinac, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_ypRB - 1.0*std::pow(-mypRB + ypRB, 2)/std::pow(sigma_ypRB, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_ypRBp - 1.0*std::pow(-mypRBp + ypRBp, 2)/std::pow(sigma_ypRBp, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_ypRBpp - 1.0*std::pow(-mypRBpp + ypRBpp, 2)/std::pow(sigma_ypRBpp, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_yE2F - 1.0*std::pow(-myE2F + yE2F, 2)/std::pow(sigma_yE2F, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_yCd - 1.0*std::pow(-myCd + yCd, 2)/std::pow(sigma_yCd, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_yMdi - 1.0*std::pow(-myMdi + yMdi, 2)/std::pow(sigma_yMdi, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_yMd - 1.0*std::pow(-myMd + yMd, 2)/std::pow(sigma_yMd, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_yMdp27 - 1.0*std::pow(-myMdp27 + yMdp27, 2)/std::pow(sigma_yMdp27, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_yCe - 1.0*std::pow(-myCe + yCe, 2)/std::pow(sigma_yCe, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_yMei - 1.0*std::pow(-myMei + yMei, 2)/std::pow(sigma_yMei, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_yMe - 1.0*std::pow(-myMe + yMe, 2)/std::pow(sigma_yMe, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_ySkp2 - 1.0*std::pow(-mySkp2 + ySkp2, 2)/std::pow(sigma_ySkp2, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_yMep27 - 1.0*std::pow(-myMep27 + yMep27, 2)/std::pow(sigma_yMep27, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_yPe - 1.0*std::pow(-myPe + yPe, 2)/std::pow(sigma_yPe, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_yPai - 1.0*std::pow(-myPai + yPai, 2)/std::pow(sigma_yPai, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_yPei - 1.0*std::pow(-myPei + yPei, 2)/std::pow(sigma_yPei, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_yPbi - 1.0*std::pow(-myPbi + yPbi, 2)/std::pow(sigma_yPbi, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_yCa - 1.0*std::pow(-myCa + yCa, 2)/std::pow(sigma_yCa, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_yMai - 1.0*std::pow(-myMai + yMai, 2)/std::pow(sigma_yMai, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_yMa - 1.0*std::pow(-myMa + yMa, 2)/std::pow(sigma_yMa, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_yMap27 - 1.0*std::pow(-myMap27 + yMap27, 2)/std::pow(sigma_yMap27, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_yp27 - 1.0*std::pow(-myp27 + yp27, 2)/std::pow(sigma_yp27, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_yCdh1i - 1.0*std::pow(-myCdh1i + yCdh1i, 2)/std::pow(sigma_yCdh1i, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_yCdh1a - 1.0*std::pow(-myCdh1a + yCdh1a, 2)/std::pow(sigma_yCdh1a, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_yE2Fp - 1.0*std::pow(-myE2Fp + yE2Fp, 2)/std::pow(sigma_yE2Fp, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_yp27p - 1.0*std::pow(-myp27p + yp27p, 2)/std::pow(sigma_yp27p, 3);
            break;
        case 65:
            dJydsigma[65] = 1.0/sigma_yPa - 1.0*std::pow(-myPa + yPa, 2)/std::pow(sigma_yPa, 3);
            break;
        case 66:
            dJydsigma[66] = 1.0/sigma_yCb - 1.0*std::pow(-myCb + yCb, 2)/std::pow(sigma_yCb, 3);
            break;
        case 67:
            dJydsigma[67] = 1.0/sigma_yMbi - 1.0*std::pow(-myMbi + yMbi, 2)/std::pow(sigma_yMbi, 3);
            break;
        case 68:
            dJydsigma[68] = 1.0/sigma_yMb - 1.0*std::pow(-myMb + yMb, 2)/std::pow(sigma_yMb, 3);
            break;
        case 69:
            dJydsigma[69] = 1.0/sigma_yCdc20i - 1.0*std::pow(-myCdc20i + yCdc20i, 2)/std::pow(sigma_yCdc20i, 3);
            break;
        case 70:
            dJydsigma[70] = 1.0/sigma_yCdc20a - 1.0*std::pow(-myCdc20a + yCdc20a, 2)/std::pow(sigma_yCdc20a, 3);
            break;
        case 71:
            dJydsigma[71] = 1.0/sigma_yPb - 1.0*std::pow(-myPb + yPb, 2)/std::pow(sigma_yPb, 3);
            break;
        case 72:
            dJydsigma[72] = 1.0/sigma_yWee1 - 1.0*std::pow(-myWee1 + yWee1, 2)/std::pow(sigma_yWee1, 3);
            break;
        case 73:
            dJydsigma[73] = 1.0/sigma_yWee1p - 1.0*std::pow(-myWee1p + yWee1p, 2)/std::pow(sigma_yWee1p, 3);
            break;
        case 74:
            dJydsigma[74] = 1.0/sigma_yMbp27 - 1.0*std::pow(-myMbp27 + yMbp27, 2)/std::pow(sigma_yMbp27, 3);
            break;
        case 75:
            dJydsigma[75] = 1.0/sigma_yChk1 - 1.0*std::pow(-myChk1 + yChk1, 2)/std::pow(sigma_yChk1, 3);
            break;
        case 76:
            dJydsigma[76] = 1.0/sigma_ypRBc1 - 1.0*std::pow(-mypRBc1 + ypRBc1, 2)/std::pow(sigma_ypRBc1, 3);
            break;
        case 77:
            dJydsigma[77] = 1.0/sigma_ypRBc2 - 1.0*std::pow(-mypRBc2 + ypRBc2, 2)/std::pow(sigma_ypRBc2, 3);
            break;
        case 78:
            dJydsigma[78] = 1.0/sigma_yp21 - 1.0*std::pow(-myp21 + yp21, 2)/std::pow(sigma_yp21, 3);
            break;
        case 79:
            dJydsigma[79] = 1.0/sigma_yMdp21 - 1.0*std::pow(-myMdp21 + yMdp21, 2)/std::pow(sigma_yMdp21, 3);
            break;
        case 80:
            dJydsigma[80] = 1.0/sigma_yMep21 - 1.0*std::pow(-myMep21 + yMep21, 2)/std::pow(sigma_yMep21, 3);
            break;
        case 81:
            dJydsigma[81] = 1.0/sigma_yMap21 - 1.0*std::pow(-myMap21 + yMap21, 2)/std::pow(sigma_yMap21, 3);
            break;
        case 82:
            dJydsigma[82] = 1.0/sigma_yMbp21 - 1.0*std::pow(-myMbp21 + yMbp21, 2)/std::pow(sigma_yMbp21, 3);
            break;
        case 83:
            dJydsigma[83] = 1.0/sigma_yL - 1.0*std::pow(-myL + yL, 2)/std::pow(sigma_yL, 3);
            break;
        case 84:
            dJydsigma[84] = 1.0/sigma_yR - 1.0*std::pow(-myR + yR, 2)/std::pow(sigma_yR, 3);
            break;
        case 85:
            dJydsigma[85] = 1.0/sigma_yL_R - 1.0*std::pow(-myL_R + yL_R, 2)/std::pow(sigma_yL_R, 3);
            break;
        case 86:
            dJydsigma[86] = 1.0/sigma_yRactive - 1.0*std::pow(-myRactive + yRactive, 2)/std::pow(sigma_yRactive, 3);
            break;
        case 87:
            dJydsigma[87] = 1.0/sigma_yflip - 1.0*std::pow(-myflip + yflip, 2)/std::pow(sigma_yflip, 3);
            break;
        case 88:
            dJydsigma[88] = 1.0/sigma_yRactive_flip - 1.0*std::pow(-myRactive_flip + yRactive_flip, 2)/std::pow(sigma_yRactive_flip, 3);
            break;
        case 89:
            dJydsigma[89] = 1.0/sigma_ypC8 - 1.0*std::pow(-mypC8 + ypC8, 2)/std::pow(sigma_ypC8, 3);
            break;
        case 90:
            dJydsigma[90] = 1.0/sigma_yRactive_pC8 - 1.0*std::pow(-myRactive_pC8 + yRactive_pC8, 2)/std::pow(sigma_yRactive_pC8, 3);
            break;
        case 91:
            dJydsigma[91] = 1.0/sigma_yC8 - 1.0*std::pow(-myC8 + yC8, 2)/std::pow(sigma_yC8, 3);
            break;
        case 92:
            dJydsigma[92] = 1.0/sigma_yBar - 1.0*std::pow(-myBar + yBar, 2)/std::pow(sigma_yBar, 3);
            break;
        case 93:
            dJydsigma[93] = 1.0/sigma_yC8_Bar - 1.0*std::pow(-myC8_Bar + yC8_Bar, 2)/std::pow(sigma_yC8_Bar, 3);
            break;
        case 94:
            dJydsigma[94] = 1.0/sigma_ypC3 - 1.0*std::pow(-mypC3 + ypC3, 2)/std::pow(sigma_ypC3, 3);
            break;
        case 95:
            dJydsigma[95] = 1.0/sigma_yC8_pC3 - 1.0*std::pow(-myC8_pC3 + yC8_pC3, 2)/std::pow(sigma_yC8_pC3, 3);
            break;
        case 96:
            dJydsigma[96] = 1.0/sigma_yC3 - 1.0*std::pow(-myC3 + yC3, 2)/std::pow(sigma_yC3, 3);
            break;
        case 97:
            dJydsigma[97] = 1.0/sigma_ypC6 - 1.0*std::pow(-mypC6 + ypC6, 2)/std::pow(sigma_ypC6, 3);
            break;
        case 98:
            dJydsigma[98] = 1.0/sigma_yC3_pC6 - 1.0*std::pow(-myC3_pC6 + yC3_pC6, 2)/std::pow(sigma_yC3_pC6, 3);
            break;
        case 99:
            dJydsigma[99] = 1.0/sigma_yC6 - 1.0*std::pow(-myC6 + yC6, 2)/std::pow(sigma_yC6, 3);
            break;
        case 100:
            dJydsigma[100] = 1.0/sigma_yC6_C8 - 1.0*std::pow(-myC6_C8 + yC6_C8, 2)/std::pow(sigma_yC6_C8, 3);
            break;
        case 101:
            dJydsigma[101] = 1.0/sigma_yXIAP - 1.0*std::pow(-myXIAP + yXIAP, 2)/std::pow(sigma_yXIAP, 3);
            break;
        case 102:
            dJydsigma[102] = 1.0/sigma_yC3_XIAP - 1.0*std::pow(-myC3_XIAP + yC3_XIAP, 2)/std::pow(sigma_yC3_XIAP, 3);
            break;
        case 103:
            dJydsigma[103] = 1.0/sigma_yPARP - 1.0*std::pow(-myPARP + yPARP, 2)/std::pow(sigma_yPARP, 3);
            break;
        case 104:
            dJydsigma[104] = 1.0/sigma_yC3_PARP - 1.0*std::pow(-myC3_PARP + yC3_PARP, 2)/std::pow(sigma_yC3_PARP, 3);
            break;
        case 105:
            dJydsigma[105] = 1.0/sigma_ycPARP - 1.0*std::pow(-mycPARP + ycPARP, 2)/std::pow(sigma_ycPARP, 3);
            break;
        case 106:
            dJydsigma[106] = 1.0/sigma_yBid - 1.0*std::pow(-myBid + yBid, 2)/std::pow(sigma_yBid, 3);
            break;
        case 107:
            dJydsigma[107] = 1.0/sigma_yC8_Bid - 1.0*std::pow(-myC8_Bid + yC8_Bid, 2)/std::pow(sigma_yC8_Bid, 3);
            break;
        case 108:
            dJydsigma[108] = 1.0/sigma_ytBid - 1.0*std::pow(-mytBid + ytBid, 2)/std::pow(sigma_ytBid, 3);
            break;
        case 109:
            dJydsigma[109] = 1.0/sigma_yBcl2c - 1.0*std::pow(-myBcl2c + yBcl2c, 2)/std::pow(sigma_yBcl2c, 3);
            break;
        case 110:
            dJydsigma[110] = 1.0/sigma_ytBid_Bcl2c - 1.0*std::pow(-mytBid_Bcl2c + ytBid_Bcl2c, 2)/std::pow(sigma_ytBid_Bcl2c, 3);
            break;
        case 111:
            dJydsigma[111] = 1.0/sigma_yBax - 1.0*std::pow(-myBax + yBax, 2)/std::pow(sigma_yBax, 3);
            break;
        case 112:
            dJydsigma[112] = 1.0/sigma_ytBid_Bax - 1.0*std::pow(-mytBid_Bax + ytBid_Bax, 2)/std::pow(sigma_ytBid_Bax, 3);
            break;
        case 113:
            dJydsigma[113] = 1.0/sigma_yBaxactive - 1.0*std::pow(-myBaxactive + yBaxactive, 2)/std::pow(sigma_yBaxactive, 3);
            break;
        case 114:
            dJydsigma[114] = 1.0/sigma_yBaxm - 1.0*std::pow(-myBaxm + yBaxm, 2)/std::pow(sigma_yBaxm, 3);
            break;
        case 115:
            dJydsigma[115] = 1.0/sigma_yBcl2 - 1.0*std::pow(-myBcl2 + yBcl2, 2)/std::pow(sigma_yBcl2, 3);
            break;
        case 116:
            dJydsigma[116] = 1.0/sigma_yBaxm_Bcl2 - 1.0*std::pow(-myBaxm_Bcl2 + yBaxm_Bcl2, 2)/std::pow(sigma_yBaxm_Bcl2, 3);
            break;
        case 117:
            dJydsigma[117] = 1.0/sigma_yBax2 - 1.0*std::pow(-myBax2 + yBax2, 2)/std::pow(sigma_yBax2, 3);
            break;
        case 118:
            dJydsigma[118] = 1.0/sigma_yBax2_Bcl2 - 1.0*std::pow(-myBax2_Bcl2 + yBax2_Bcl2, 2)/std::pow(sigma_yBax2_Bcl2, 3);
            break;
        case 119:
            dJydsigma[119] = 1.0/sigma_yBax4 - 1.0*std::pow(-myBax4 + yBax4, 2)/std::pow(sigma_yBax4, 3);
            break;
        case 120:
            dJydsigma[120] = 1.0/sigma_yBax4_Bcl2 - 1.0*std::pow(-myBax4_Bcl2 + yBax4_Bcl2, 2)/std::pow(sigma_yBax4_Bcl2, 3);
            break;
        case 121:
            dJydsigma[121] = 1.0/sigma_yM - 1.0*std::pow(-myM + yM, 2)/std::pow(sigma_yM, 3);
            break;
        case 122:
            dJydsigma[122] = 1.0/sigma_yBax4_M - 1.0*std::pow(-myBax4_M + yBax4_M, 2)/std::pow(sigma_yBax4_M, 3);
            break;
        case 123:
            dJydsigma[123] = 1.0/sigma_yMactive - 1.0*std::pow(-myMactive + yMactive, 2)/std::pow(sigma_yMactive, 3);
            break;
        case 124:
            dJydsigma[124] = 1.0/sigma_yCytoCm - 1.0*std::pow(-myCytoCm + yCytoCm, 2)/std::pow(sigma_yCytoCm, 3);
            break;
        case 125:
            dJydsigma[125] = 1.0/sigma_yMactive_CytoCm - 1.0*std::pow(-myMactive_CytoCm + yMactive_CytoCm, 2)/std::pow(sigma_yMactive_CytoCm, 3);
            break;
        case 126:
            dJydsigma[126] = 1.0/sigma_yCytoCr - 1.0*std::pow(-myCytoCr + yCytoCr, 2)/std::pow(sigma_yCytoCr, 3);
            break;
        case 127:
            dJydsigma[127] = 1.0/sigma_ySmacm - 1.0*std::pow(-mySmacm + ySmacm, 2)/std::pow(sigma_ySmacm, 3);
            break;
        case 128:
            dJydsigma[128] = 1.0/sigma_yMactive_Smacm - 1.0*std::pow(-myMactive_Smacm + yMactive_Smacm, 2)/std::pow(sigma_yMactive_Smacm, 3);
            break;
        case 129:
            dJydsigma[129] = 1.0/sigma_ySmacr - 1.0*std::pow(-mySmacr + ySmacr, 2)/std::pow(sigma_ySmacr, 3);
            break;
        case 130:
            dJydsigma[130] = 1.0/sigma_yCytoC - 1.0*std::pow(-myCytoC + yCytoC, 2)/std::pow(sigma_yCytoC, 3);
            break;
        case 131:
            dJydsigma[131] = 1.0/sigma_yApaf - 1.0*std::pow(-myApaf + yApaf, 2)/std::pow(sigma_yApaf, 3);
            break;
        case 132:
            dJydsigma[132] = 1.0/sigma_yCytoC_Apaf - 1.0*std::pow(-myCytoC_Apaf + yCytoC_Apaf, 2)/std::pow(sigma_yCytoC_Apaf, 3);
            break;
        case 133:
            dJydsigma[133] = 1.0/sigma_yApafactive - 1.0*std::pow(-myApafactive + yApafactive, 2)/std::pow(sigma_yApafactive, 3);
            break;
        case 134:
            dJydsigma[134] = 1.0/sigma_ypC9 - 1.0*std::pow(-mypC9 + ypC9, 2)/std::pow(sigma_ypC9, 3);
            break;
        case 135:
            dJydsigma[135] = 1.0/sigma_yApop - 1.0*std::pow(-myApop + yApop, 2)/std::pow(sigma_yApop, 3);
            break;
        case 136:
            dJydsigma[136] = 1.0/sigma_yApop_C3 - 1.0*std::pow(-myApop_C3 + yApop_C3, 2)/std::pow(sigma_yApop_C3, 3);
            break;
        case 137:
            dJydsigma[137] = 1.0/sigma_ySmac - 1.0*std::pow(-mySmac + ySmac, 2)/std::pow(sigma_ySmac, 3);
            break;
        case 138:
            dJydsigma[138] = 1.0/sigma_yApop_XIAP - 1.0*std::pow(-myApop_XIAP + yApop_XIAP, 2)/std::pow(sigma_yApop_XIAP, 3);
            break;
        case 139:
            dJydsigma[139] = 1.0/sigma_ySmac_XIAP - 1.0*std::pow(-mySmac_XIAP + ySmac_XIAP, 2)/std::pow(sigma_ySmac_XIAP, 3);
            break;
        case 140:
            dJydsigma[140] = 1.0/sigma_yC3_Ub - 1.0*std::pow(-myC3_Ub + yC3_Ub, 2)/std::pow(sigma_yC3_Ub, 3);
            break;
        case 141:
            dJydsigma[141] = 1.0/sigma_yBAD - 1.0*std::pow(-myBAD + yBAD, 2)/std::pow(sigma_yBAD, 3);
            break;
        case 142:
            dJydsigma[142] = 1.0/sigma_yPUMA - 1.0*std::pow(-myPUMA + yPUMA, 2)/std::pow(sigma_yPUMA, 3);
            break;
        case 143:
            dJydsigma[143] = 1.0/sigma_yNOXA - 1.0*std::pow(-myNOXA + yNOXA, 2)/std::pow(sigma_yNOXA, 3);
            break;
        case 144:
            dJydsigma[144] = 1.0/sigma_yBcl2c_BAD - 1.0*std::pow(-myBcl2c_BAD + yBcl2c_BAD, 2)/std::pow(sigma_yBcl2c_BAD, 3);
            break;
        case 145:
            dJydsigma[145] = 1.0/sigma_yBcl2c_PUMA - 1.0*std::pow(-myBcl2c_PUMA + yBcl2c_PUMA, 2)/std::pow(sigma_yBcl2c_PUMA, 3);
            break;
        case 146:
            dJydsigma[146] = 1.0/sigma_yBcl2c_NOXA - 1.0*std::pow(-myBcl2c_NOXA + yBcl2c_NOXA, 2)/std::pow(sigma_yBcl2c_NOXA, 3);
            break;
        case 147:
            dJydsigma[147] = 1.0/sigma_yBIM - 1.0*std::pow(-myBIM + yBIM, 2)/std::pow(sigma_yBIM, 3);
            break;
        case 148:
            dJydsigma[148] = 1.0/sigma_yBIM_Bax - 1.0*std::pow(-myBIM_Bax + yBIM_Bax, 2)/std::pow(sigma_yBIM_Bax, 3);
            break;
        case 149:
            dJydsigma[149] = 1.0/sigma_yBcl2c_BIM - 1.0*std::pow(-myBcl2c_BIM + yBcl2c_BIM, 2)/std::pow(sigma_yBcl2c_BIM, 3);
            break;
        case 150:
            dJydsigma[150] = 1.0/sigma_yppERK_BIM - 1.0*std::pow(-myppERK_BIM + yppERK_BIM, 2)/std::pow(sigma_yppERK_BIM, 3);
            break;
        case 151:
            dJydsigma[151] = 1.0/sigma_ypBIM - 1.0*std::pow(-mypBIM + ypBIM, 2)/std::pow(sigma_ypBIM, 3);
            break;
        case 152:
            dJydsigma[152] = 1.0/sigma_yppAKT_BAD - 1.0*std::pow(-myppAKT_BAD + yppAKT_BAD, 2)/std::pow(sigma_yppAKT_BAD, 3);
            break;
        case 153:
            dJydsigma[153] = 1.0/sigma_ypBAD - 1.0*std::pow(-mypBAD + ypBAD, 2)/std::pow(sigma_ypBAD, 3);
            break;
        case 154:
            dJydsigma[154] = 1.0/sigma_yppERK_BAD - 1.0*std::pow(-myppERK_BAD + yppERK_BAD, 2)/std::pow(sigma_yppERK_BAD, 3);
            break;
        case 155:
            dJydsigma[155] = 1.0/sigma_yE - 1.0*std::pow(-myE + yE, 2)/std::pow(sigma_yE, 3);
            break;
        case 156:
            dJydsigma[156] = 1.0/sigma_yH - 1.0*std::pow(-myH + yH, 2)/std::pow(sigma_yH, 3);
            break;
        case 157:
            dJydsigma[157] = 1.0/sigma_yHGF - 1.0*std::pow(-myHGF + yHGF, 2)/std::pow(sigma_yHGF, 3);
            break;
        case 158:
            dJydsigma[158] = 1.0/sigma_yP - 1.0*std::pow(-myP + yP, 2)/std::pow(sigma_yP, 3);
            break;
        case 159:
            dJydsigma[159] = 1.0/sigma_yF - 1.0*std::pow(-myF + yF, 2)/std::pow(sigma_yF, 3);
            break;
        case 160:
            dJydsigma[160] = 1.0/sigma_yI - 1.0*std::pow(-myI + yI, 2)/std::pow(sigma_yI, 3);
            break;
        case 161:
            dJydsigma[161] = 1.0/sigma_yINS - 1.0*std::pow(-myINS + yINS, 2)/std::pow(sigma_yINS, 3);
            break;
        case 162:
            dJydsigma[162] = 1.0/sigma_yE1 - 1.0*std::pow(-myE1 + yE1, 2)/std::pow(sigma_yE1, 3);
            break;
        case 163:
            dJydsigma[163] = 1.0/sigma_ypE1 - 1.0*std::pow(-mypE1 + ypE1, 2)/std::pow(sigma_ypE1, 3);
            break;
        case 164:
            dJydsigma[164] = 1.0/sigma_yE2 - 1.0*std::pow(-myE2 + yE2, 2)/std::pow(sigma_yE2, 3);
            break;
        case 165:
            dJydsigma[165] = 1.0/sigma_ypE2 - 1.0*std::pow(-mypE2 + ypE2, 2)/std::pow(sigma_ypE2, 3);
            break;
        case 166:
            dJydsigma[166] = 1.0/sigma_yE3 - 1.0*std::pow(-myE3 + yE3, 2)/std::pow(sigma_yE3, 3);
            break;
        case 167:
            dJydsigma[167] = 1.0/sigma_yE4 - 1.0*std::pow(-myE4 + yE4, 2)/std::pow(sigma_yE4, 3);
            break;
        case 168:
            dJydsigma[168] = 1.0/sigma_ypE4 - 1.0*std::pow(-mypE4 + ypE4, 2)/std::pow(sigma_ypE4, 3);
            break;
        case 169:
            dJydsigma[169] = 1.0/sigma_yEv3 - 1.0*std::pow(-myEv3 + yEv3, 2)/std::pow(sigma_yEv3, 3);
            break;
        case 170:
            dJydsigma[170] = 1.0/sigma_yMet - 1.0*std::pow(-myMet + yMet, 2)/std::pow(sigma_yMet, 3);
            break;
        case 171:
            dJydsigma[171] = 1.0/sigma_yPr - 1.0*std::pow(-myPr + yPr, 2)/std::pow(sigma_yPr, 3);
            break;
        case 172:
            dJydsigma[172] = 1.0/sigma_yFr - 1.0*std::pow(-myFr + yFr, 2)/std::pow(sigma_yFr, 3);
            break;
        case 173:
            dJydsigma[173] = 1.0/sigma_yIr - 1.0*std::pow(-myIr + yIr, 2)/std::pow(sigma_yIr, 3);
            break;
        case 174:
            dJydsigma[174] = 1.0/sigma_yIsr - 1.0*std::pow(-myIsr + yIsr, 2)/std::pow(sigma_yIsr, 3);
            break;
        case 175:
            dJydsigma[175] = 1.0/sigma_yE1E1 - 1.0*std::pow(-myE1E1 + yE1E1, 2)/std::pow(sigma_yE1E1, 3);
            break;
        case 176:
            dJydsigma[176] = 1.0/sigma_yE1E2 - 1.0*std::pow(-myE1E2 + yE1E2, 2)/std::pow(sigma_yE1E2, 3);
            break;
        case 177:
            dJydsigma[177] = 1.0/sigma_yE1E3 - 1.0*std::pow(-myE1E3 + yE1E3, 2)/std::pow(sigma_yE1E3, 3);
            break;
        case 178:
            dJydsigma[178] = 1.0/sigma_yE1E4 - 1.0*std::pow(-myE1E4 + yE1E4, 2)/std::pow(sigma_yE1E4, 3);
            break;
        case 179:
            dJydsigma[179] = 1.0/sigma_yE2E2 - 1.0*std::pow(-myE2E2 + yE2E2, 2)/std::pow(sigma_yE2E2, 3);
            break;
        case 180:
            dJydsigma[180] = 1.0/sigma_yE2E3 - 1.0*std::pow(-myE2E3 + yE2E3, 2)/std::pow(sigma_yE2E3, 3);
            break;
        case 181:
            dJydsigma[181] = 1.0/sigma_yE2E4 - 1.0*std::pow(-myE2E4 + yE2E4, 2)/std::pow(sigma_yE2E4, 3);
            break;
        case 182:
            dJydsigma[182] = 1.0/sigma_yE3E4 - 1.0*std::pow(-myE3E4 + yE3E4, 2)/std::pow(sigma_yE3E4, 3);
            break;
        case 183:
            dJydsigma[183] = 1.0/sigma_yE4E4 - 1.0*std::pow(-myE4E4 + yE4E4, 2)/std::pow(sigma_yE4E4, 3);
            break;
        case 184:
            dJydsigma[184] = 1.0/sigma_yMet_Met - 1.0*std::pow(-myMet_Met + yMet_Met, 2)/std::pow(sigma_yMet_Met, 3);
            break;
        case 185:
            dJydsigma[185] = 1.0/sigma_yFrFr - 1.0*std::pow(-myFrFr + yFrFr, 2)/std::pow(sigma_yFrFr, 3);
            break;
        case 186:
            dJydsigma[186] = 1.0/sigma_yIrIr - 1.0*std::pow(-myIrIr + yIrIr, 2)/std::pow(sigma_yIrIr, 3);
            break;
        case 187:
            dJydsigma[187] = 1.0/sigma_yIsr_Isr - 1.0*std::pow(-myIsr_Isr + yIsr_Isr, 2)/std::pow(sigma_yIsr_Isr, 3);
            break;
        case 188:
            dJydsigma[188] = 1.0/sigma_yEE1 - 1.0*std::pow(-myEE1 + yEE1, 2)/std::pow(sigma_yEE1, 3);
            break;
        case 189:
            dJydsigma[189] = 1.0/sigma_yHE3 - 1.0*std::pow(-myHE3 + yHE3, 2)/std::pow(sigma_yHE3, 3);
            break;
        case 190:
            dJydsigma[190] = 1.0/sigma_yHE4 - 1.0*std::pow(-myHE4 + yHE4, 2)/std::pow(sigma_yHE4, 3);
            break;
        case 191:
            dJydsigma[191] = 1.0/sigma_yHGF_Met - 1.0*std::pow(-myHGF_Met + yHGF_Met, 2)/std::pow(sigma_yHGF_Met, 3);
            break;
        case 192:
            dJydsigma[192] = 1.0/sigma_yPPr - 1.0*std::pow(-myPPr + yPPr, 2)/std::pow(sigma_yPPr, 3);
            break;
        case 193:
            dJydsigma[193] = 1.0/sigma_yFFr - 1.0*std::pow(-myFFr + yFFr, 2)/std::pow(sigma_yFFr, 3);
            break;
        case 194:
            dJydsigma[194] = 1.0/sigma_yEE1E2 - 1.0*std::pow(-myEE1E2 + yEE1E2, 2)/std::pow(sigma_yEE1E2, 3);
            break;
        case 195:
            dJydsigma[195] = 1.0/sigma_yEE1Ev3 - 1.0*std::pow(-myEE1Ev3 + yEE1Ev3, 2)/std::pow(sigma_yEE1Ev3, 3);
            break;
        case 196:
            dJydsigma[196] = 1.0/sigma_yEE1E1 - 1.0*std::pow(-myEE1E1 + yEE1E1, 2)/std::pow(sigma_yEE1E1, 3);
            break;
        case 197:
            dJydsigma[197] = 1.0/sigma_yEE1E3 - 1.0*std::pow(-myEE1E3 + yEE1E3, 2)/std::pow(sigma_yEE1E3, 3);
            break;
        case 198:
            dJydsigma[198] = 1.0/sigma_yEE1E4 - 1.0*std::pow(-myEE1E4 + yEE1E4, 2)/std::pow(sigma_yEE1E4, 3);
            break;
        case 199:
            dJydsigma[199] = 1.0/sigma_yE2HE3 - 1.0*std::pow(-myE2HE3 + yE2HE3, 2)/std::pow(sigma_yE2HE3, 3);
            break;
        case 200:
            dJydsigma[200] = 1.0/sigma_yE1HE3 - 1.0*std::pow(-myE1HE3 + yE1HE3, 2)/std::pow(sigma_yE1HE3, 3);
            break;
        case 201:
            dJydsigma[201] = 1.0/sigma_yHE3E3 - 1.0*std::pow(-myHE3E3 + yHE3E3, 2)/std::pow(sigma_yHE3E3, 3);
            break;
        case 202:
            dJydsigma[202] = 1.0/sigma_yHE3Ev3 - 1.0*std::pow(-myHE3Ev3 + yHE3Ev3, 2)/std::pow(sigma_yHE3Ev3, 3);
            break;
        case 203:
            dJydsigma[203] = 1.0/sigma_yHE3E4 - 1.0*std::pow(-myHE3E4 + yHE3E4, 2)/std::pow(sigma_yHE3E4, 3);
            break;
        case 204:
            dJydsigma[204] = 1.0/sigma_yE2HE4 - 1.0*std::pow(-myE2HE4 + yE2HE4, 2)/std::pow(sigma_yE2HE4, 3);
            break;
        case 205:
            dJydsigma[205] = 1.0/sigma_yHE4Ev3 - 1.0*std::pow(-myHE4Ev3 + yHE4Ev3, 2)/std::pow(sigma_yHE4Ev3, 3);
            break;
        case 206:
            dJydsigma[206] = 1.0/sigma_yE1HE4 - 1.0*std::pow(-myE1HE4 + yE1HE4, 2)/std::pow(sigma_yE1HE4, 3);
            break;
        case 207:
            dJydsigma[207] = 1.0/sigma_yE3HE4 - 1.0*std::pow(-myE3HE4 + yE3HE4, 2)/std::pow(sigma_yE3HE4, 3);
            break;
        case 208:
            dJydsigma[208] = 1.0/sigma_yHE4E4 - 1.0*std::pow(-myHE4E4 + yHE4E4, 2)/std::pow(sigma_yHE4E4, 3);
            break;
        case 209:
            dJydsigma[209] = 1.0/sigma_yHGF_Met_Met - 1.0*std::pow(-myHGF_Met_Met + yHGF_Met_Met, 2)/std::pow(sigma_yHGF_Met_Met, 3);
            break;
        case 210:
            dJydsigma[210] = 1.0/sigma_yPPrPr - 1.0*std::pow(-myPPrPr + yPPrPr, 2)/std::pow(sigma_yPPrPr, 3);
            break;
        case 211:
            dJydsigma[211] = 1.0/sigma_yFFrFr - 1.0*std::pow(-myFFrFr + yFFrFr, 2)/std::pow(sigma_yFFrFr, 3);
            break;
        case 212:
            dJydsigma[212] = 1.0/sigma_yIIrIr - 1.0*std::pow(-myIIrIr + yIIrIr, 2)/std::pow(sigma_yIIrIr, 3);
            break;
        case 213:
            dJydsigma[213] = 1.0/sigma_yINS_Isr_Isr - 1.0*std::pow(-myINS_Isr_Isr + yINS_Isr_Isr, 2)/std::pow(sigma_yINS_Isr_Isr, 3);
            break;
        case 214:
            dJydsigma[214] = 1.0/sigma_yEE1EE1 - 1.0*std::pow(-myEE1EE1 + yEE1EE1, 2)/std::pow(sigma_yEE1EE1, 3);
            break;
        case 215:
            dJydsigma[215] = 1.0/sigma_yEE1HE3 - 1.0*std::pow(-myEE1HE3 + yEE1HE3, 2)/std::pow(sigma_yEE1HE3, 3);
            break;
        case 216:
            dJydsigma[216] = 1.0/sigma_yEE1HE4 - 1.0*std::pow(-myEE1HE4 + yEE1HE4, 2)/std::pow(sigma_yEE1HE4, 3);
            break;
        case 217:
            dJydsigma[217] = 1.0/sigma_yHE3HE3 - 1.0*std::pow(-myHE3HE3 + yHE3HE3, 2)/std::pow(sigma_yHE3HE3, 3);
            break;
        case 218:
            dJydsigma[218] = 1.0/sigma_yHE3HE4 - 1.0*std::pow(-myHE3HE4 + yHE3HE4, 2)/std::pow(sigma_yHE3HE4, 3);
            break;
        case 219:
            dJydsigma[219] = 1.0/sigma_yHE4HE4 - 1.0*std::pow(-myHE4HE4 + yHE4HE4, 2)/std::pow(sigma_yHE4HE4, 3);
            break;
        case 220:
            dJydsigma[220] = 1.0/sigma_yHGF_Met_HGF_Met - 1.0*std::pow(-myHGF_Met_HGF_Met + yHGF_Met_HGF_Met, 2)/std::pow(sigma_yHGF_Met_HGF_Met, 3);
            break;
        case 221:
            dJydsigma[221] = 1.0/sigma_yPPrPPr - 1.0*std::pow(-myPPrPPr + yPPrPPr, 2)/std::pow(sigma_yPPrPPr, 3);
            break;
        case 222:
            dJydsigma[222] = 1.0/sigma_yFFrFFr - 1.0*std::pow(-myFFrFFr + yFFrFFr, 2)/std::pow(sigma_yFFrFFr, 3);
            break;
        case 223:
            dJydsigma[223] = 1.0/sigma_yIIrIrI - 1.0*std::pow(-myIIrIrI + yIIrIrI, 2)/std::pow(sigma_yIIrIrI, 3);
            break;
        case 224:
            dJydsigma[224] = 1.0/sigma_yINS_Isr_Isr_INS - 1.0*std::pow(-myINS_Isr_Isr_INS + yINS_Isr_Isr_INS, 2)/std::pow(sigma_yINS_Isr_Isr_INS, 3);
            break;
        case 225:
            dJydsigma[225] = 1.0/sigma_yE1_ppERK - 1.0*std::pow(-myE1_ppERK + yE1_ppERK, 2)/std::pow(sigma_yE1_ppERK, 3);
            break;
        case 226:
            dJydsigma[226] = 1.0/sigma_yE2_ppERK - 1.0*std::pow(-myE2_ppERK + yE2_ppERK, 2)/std::pow(sigma_yE2_ppERK, 3);
            break;
        case 227:
            dJydsigma[227] = 1.0/sigma_yE4_ppERK - 1.0*std::pow(-myE4_ppERK + yE4_ppERK, 2)/std::pow(sigma_yE4_ppERK, 3);
            break;
        case 228:
            dJydsigma[228] = 1.0/sigma_ypEE1E2 - 1.0*std::pow(-mypEE1E2 + ypEE1E2, 2)/std::pow(sigma_ypEE1E2, 3);
            break;
        case 229:
            dJydsigma[229] = 1.0/sigma_ypEE1Ev3 - 1.0*std::pow(-mypEE1Ev3 + ypEE1Ev3, 2)/std::pow(sigma_ypEE1Ev3, 3);
            break;
        case 230:
            dJydsigma[230] = 1.0/sigma_ypEE1E1 - 1.0*std::pow(-mypEE1E1 + ypEE1E1, 2)/std::pow(sigma_ypEE1E1, 3);
            break;
        case 231:
            dJydsigma[231] = 1.0/sigma_ypEE1EE1 - 1.0*std::pow(-mypEE1EE1 + ypEE1EE1, 2)/std::pow(sigma_ypEE1EE1, 3);
            break;
        case 232:
            dJydsigma[232] = 1.0/sigma_ypEE1E3 - 1.0*std::pow(-mypEE1E3 + ypEE1E3, 2)/std::pow(sigma_ypEE1E3, 3);
            break;
        case 233:
            dJydsigma[233] = 1.0/sigma_ypEE1HE3 - 1.0*std::pow(-mypEE1HE3 + ypEE1HE3, 2)/std::pow(sigma_ypEE1HE3, 3);
            break;
        case 234:
            dJydsigma[234] = 1.0/sigma_ypEE1E4 - 1.0*std::pow(-mypEE1E4 + ypEE1E4, 2)/std::pow(sigma_ypEE1E4, 3);
            break;
        case 235:
            dJydsigma[235] = 1.0/sigma_ypEE1HE4 - 1.0*std::pow(-mypEE1HE4 + ypEE1HE4, 2)/std::pow(sigma_ypEE1HE4, 3);
            break;
        case 236:
            dJydsigma[236] = 1.0/sigma_ypE2HE3 - 1.0*std::pow(-mypE2HE3 + ypE2HE3, 2)/std::pow(sigma_ypE2HE3, 3);
            break;
        case 237:
            dJydsigma[237] = 1.0/sigma_ypHE3Ev3 - 1.0*std::pow(-mypHE3Ev3 + ypHE3Ev3, 2)/std::pow(sigma_ypHE3Ev3, 3);
            break;
        case 238:
            dJydsigma[238] = 1.0/sigma_ypE1HE3 - 1.0*std::pow(-mypE1HE3 + ypE1HE3, 2)/std::pow(sigma_ypE1HE3, 3);
            break;
        case 239:
            dJydsigma[239] = 1.0/sigma_ypHE3E4 - 1.0*std::pow(-mypHE3E4 + ypHE3E4, 2)/std::pow(sigma_ypHE3E4, 3);
            break;
        case 240:
            dJydsigma[240] = 1.0/sigma_ypHE3HE4 - 1.0*std::pow(-mypHE3HE4 + ypHE3HE4, 2)/std::pow(sigma_ypHE3HE4, 3);
            break;
        case 241:
            dJydsigma[241] = 1.0/sigma_ypE2HE4 - 1.0*std::pow(-mypE2HE4 + ypE2HE4, 2)/std::pow(sigma_ypE2HE4, 3);
            break;
        case 242:
            dJydsigma[242] = 1.0/sigma_ypHE4Ev3 - 1.0*std::pow(-mypHE4Ev3 + ypHE4Ev3, 2)/std::pow(sigma_ypHE4Ev3, 3);
            break;
        case 243:
            dJydsigma[243] = 1.0/sigma_ypE1HE4 - 1.0*std::pow(-mypE1HE4 + ypE1HE4, 2)/std::pow(sigma_ypE1HE4, 3);
            break;
        case 244:
            dJydsigma[244] = 1.0/sigma_ypE3HE4 - 1.0*std::pow(-mypE3HE4 + ypE3HE4, 2)/std::pow(sigma_ypE3HE4, 3);
            break;
        case 245:
            dJydsigma[245] = 1.0/sigma_ypHE4E4 - 1.0*std::pow(-mypHE4E4 + ypHE4E4, 2)/std::pow(sigma_ypHE4E4, 3);
            break;
        case 246:
            dJydsigma[246] = 1.0/sigma_ypHE4HE4 - 1.0*std::pow(-mypHE4HE4 + ypHE4HE4, 2)/std::pow(sigma_ypHE4HE4, 3);
            break;
        case 247:
            dJydsigma[247] = 1.0/sigma_ypHGF_Met_Met - 1.0*std::pow(-mypHGF_Met_Met + ypHGF_Met_Met, 2)/std::pow(sigma_ypHGF_Met_Met, 3);
            break;
        case 248:
            dJydsigma[248] = 1.0/sigma_ypHGF_Met_HGF_Met - 1.0*std::pow(-mypHGF_Met_HGF_Met + ypHGF_Met_HGF_Met, 2)/std::pow(sigma_ypHGF_Met_HGF_Met, 3);
            break;
        case 249:
            dJydsigma[249] = 1.0/sigma_ypPPrPPr - 1.0*std::pow(-mypPPrPPr + ypPPrPPr, 2)/std::pow(sigma_ypPPrPPr, 3);
            break;
        case 250:
            dJydsigma[250] = 1.0/sigma_ypPPrPr - 1.0*std::pow(-mypPPrPr + ypPPrPr, 2)/std::pow(sigma_ypPPrPr, 3);
            break;
        case 251:
            dJydsigma[251] = 1.0/sigma_ypFFrFFr - 1.0*std::pow(-mypFFrFFr + ypFFrFFr, 2)/std::pow(sigma_ypFFrFFr, 3);
            break;
        case 252:
            dJydsigma[252] = 1.0/sigma_ypFFrFr - 1.0*std::pow(-mypFFrFr + ypFFrFr, 2)/std::pow(sigma_ypFFrFr, 3);
            break;
        case 253:
            dJydsigma[253] = 1.0/sigma_ypIIrIr - 1.0*std::pow(-mypIIrIr + ypIIrIr, 2)/std::pow(sigma_ypIIrIr, 3);
            break;
        case 254:
            dJydsigma[254] = 1.0/sigma_ypINS_Isr_Isr - 1.0*std::pow(-mypINS_Isr_Isr + ypINS_Isr_Isr, 2)/std::pow(sigma_ypINS_Isr_Isr, 3);
            break;
        case 255:
            dJydsigma[255] = 1.0/sigma_ypIIrIrI - 1.0*std::pow(-mypIIrIrI + ypIIrIrI, 2)/std::pow(sigma_ypIIrIrI, 3);
            break;
        case 256:
            dJydsigma[256] = 1.0/sigma_ypINS_Isr_Isr_INS - 1.0*std::pow(-mypINS_Isr_Isr_INS + ypINS_Isr_Isr_INS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS, 3);
            break;
        case 257:
            dJydsigma[257] = 1.0/sigma_ypIIrIr_IRS - 1.0*std::pow(-mypIIrIr_IRS + ypIIrIr_IRS, 2)/std::pow(sigma_ypIIrIr_IRS, 3);
            break;
        case 258:
            dJydsigma[258] = 1.0/sigma_ypINS_Isr_Isr_IRS - 1.0*std::pow(-mypINS_Isr_Isr_IRS + ypINS_Isr_Isr_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS, 3);
            break;
        case 259:
            dJydsigma[259] = 1.0/sigma_ypIIrIrI_IRS - 1.0*std::pow(-mypIIrIrI_IRS + ypIIrIrI_IRS, 2)/std::pow(sigma_ypIIrIrI_IRS, 3);
            break;
        case 260:
            dJydsigma[260] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS + ypINS_Isr_Isr_INS_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS, 3);
            break;
        case 261:
            dJydsigma[261] = 1.0/sigma_ySp_EE1E2 - 1.0*std::pow(-mySp_EE1E2 + ySp_EE1E2, 2)/std::pow(sigma_ySp_EE1E2, 3);
            break;
        case 262:
            dJydsigma[262] = 1.0/sigma_ySp_EE1Ev3 - 1.0*std::pow(-mySp_EE1Ev3 + ySp_EE1Ev3, 2)/std::pow(sigma_ySp_EE1Ev3, 3);
            break;
        case 263:
            dJydsigma[263] = 1.0/sigma_ySp_EE1E1 - 1.0*std::pow(-mySp_EE1E1 + ySp_EE1E1, 2)/std::pow(sigma_ySp_EE1E1, 3);
            break;
        case 264:
            dJydsigma[264] = 1.0/sigma_ySp_EE1EE1 - 1.0*std::pow(-mySp_EE1EE1 + ySp_EE1EE1, 2)/std::pow(sigma_ySp_EE1EE1, 3);
            break;
        case 265:
            dJydsigma[265] = 1.0/sigma_ySp_EE1E3 - 1.0*std::pow(-mySp_EE1E3 + ySp_EE1E3, 2)/std::pow(sigma_ySp_EE1E3, 3);
            break;
        case 266:
            dJydsigma[266] = 1.0/sigma_ySp_EE1HE3 - 1.0*std::pow(-mySp_EE1HE3 + ySp_EE1HE3, 2)/std::pow(sigma_ySp_EE1HE3, 3);
            break;
        case 267:
            dJydsigma[267] = 1.0/sigma_ySp_EE1E4 - 1.0*std::pow(-mySp_EE1E4 + ySp_EE1E4, 2)/std::pow(sigma_ySp_EE1E4, 3);
            break;
        case 268:
            dJydsigma[268] = 1.0/sigma_ySp_EE1HE4 - 1.0*std::pow(-mySp_EE1HE4 + ySp_EE1HE4, 2)/std::pow(sigma_ySp_EE1HE4, 3);
            break;
        case 269:
            dJydsigma[269] = 1.0/sigma_ySp_E2HE3 - 1.0*std::pow(-mySp_E2HE3 + ySp_E2HE3, 2)/std::pow(sigma_ySp_E2HE3, 3);
            break;
        case 270:
            dJydsigma[270] = 1.0/sigma_ySp_HE3Ev3 - 1.0*std::pow(-mySp_HE3Ev3 + ySp_HE3Ev3, 2)/std::pow(sigma_ySp_HE3Ev3, 3);
            break;
        case 271:
            dJydsigma[271] = 1.0/sigma_ySp_E1HE3 - 1.0*std::pow(-mySp_E1HE3 + ySp_E1HE3, 2)/std::pow(sigma_ySp_E1HE3, 3);
            break;
        case 272:
            dJydsigma[272] = 1.0/sigma_ySp_HE3E4 - 1.0*std::pow(-mySp_HE3E4 + ySp_HE3E4, 2)/std::pow(sigma_ySp_HE3E4, 3);
            break;
        case 273:
            dJydsigma[273] = 1.0/sigma_ySp_HE3HE4 - 1.0*std::pow(-mySp_HE3HE4 + ySp_HE3HE4, 2)/std::pow(sigma_ySp_HE3HE4, 3);
            break;
        case 274:
            dJydsigma[274] = 1.0/sigma_ySp_E2HE4 - 1.0*std::pow(-mySp_E2HE4 + ySp_E2HE4, 2)/std::pow(sigma_ySp_E2HE4, 3);
            break;
        case 275:
            dJydsigma[275] = 1.0/sigma_ySp_HE4Ev3 - 1.0*std::pow(-mySp_HE4Ev3 + ySp_HE4Ev3, 2)/std::pow(sigma_ySp_HE4Ev3, 3);
            break;
        case 276:
            dJydsigma[276] = 1.0/sigma_ySp_E1HE4 - 1.0*std::pow(-mySp_E1HE4 + ySp_E1HE4, 2)/std::pow(sigma_ySp_E1HE4, 3);
            break;
        case 277:
            dJydsigma[277] = 1.0/sigma_ySp_E3HE4 - 1.0*std::pow(-mySp_E3HE4 + ySp_E3HE4, 2)/std::pow(sigma_ySp_E3HE4, 3);
            break;
        case 278:
            dJydsigma[278] = 1.0/sigma_ySp_HE4E4 - 1.0*std::pow(-mySp_HE4E4 + ySp_HE4E4, 2)/std::pow(sigma_ySp_HE4E4, 3);
            break;
        case 279:
            dJydsigma[279] = 1.0/sigma_ySp_HE4HE4 - 1.0*std::pow(-mySp_HE4HE4 + ySp_HE4HE4, 2)/std::pow(sigma_ySp_HE4HE4, 3);
            break;
        case 280:
            dJydsigma[280] = 1.0/sigma_ySp_HGF_Met_Met - 1.0*std::pow(-mySp_HGF_Met_Met + ySp_HGF_Met_Met, 2)/std::pow(sigma_ySp_HGF_Met_Met, 3);
            break;
        case 281:
            dJydsigma[281] = 1.0/sigma_ySp_HGF_Met_HGF_Met - 1.0*std::pow(-mySp_HGF_Met_HGF_Met + ySp_HGF_Met_HGF_Met, 2)/std::pow(sigma_ySp_HGF_Met_HGF_Met, 3);
            break;
        case 282:
            dJydsigma[282] = 1.0/sigma_ySp_PPrPPr - 1.0*std::pow(-mySp_PPrPPr + ySp_PPrPPr, 2)/std::pow(sigma_ySp_PPrPPr, 3);
            break;
        case 283:
            dJydsigma[283] = 1.0/sigma_ySp_PPrPr - 1.0*std::pow(-mySp_PPrPr + ySp_PPrPr, 2)/std::pow(sigma_ySp_PPrPr, 3);
            break;
        case 284:
            dJydsigma[284] = 1.0/sigma_ySp_FFrFFr - 1.0*std::pow(-mySp_FFrFFr + ySp_FFrFFr, 2)/std::pow(sigma_ySp_FFrFFr, 3);
            break;
        case 285:
            dJydsigma[285] = 1.0/sigma_ySp_FFrFr - 1.0*std::pow(-mySp_FFrFr + ySp_FFrFr, 2)/std::pow(sigma_ySp_FFrFr, 3);
            break;
        case 286:
            dJydsigma[286] = 1.0/sigma_ySp_IIrIr - 1.0*std::pow(-mySp_IIrIr + ySp_IIrIr, 2)/std::pow(sigma_ySp_IIrIr, 3);
            break;
        case 287:
            dJydsigma[287] = 1.0/sigma_ySp_INS_Isr_Isr - 1.0*std::pow(-mySp_INS_Isr_Isr + ySp_INS_Isr_Isr, 2)/std::pow(sigma_ySp_INS_Isr_Isr, 3);
            break;
        case 288:
            dJydsigma[288] = 1.0/sigma_ySp_IIrIrI - 1.0*std::pow(-mySp_IIrIrI + ySp_IIrIrI, 2)/std::pow(sigma_ySp_IIrIrI, 3);
            break;
        case 289:
            dJydsigma[289] = 1.0/sigma_ySp_INS_Isr_Isr_INS - 1.0*std::pow(-mySp_INS_Isr_Isr_INS + ySp_INS_Isr_Isr_INS, 2)/std::pow(sigma_ySp_INS_Isr_Isr_INS, 3);
            break;
        case 290:
            dJydsigma[290] = 1.0/sigma_yEE1E2int - 1.0*std::pow(-myEE1E2int + yEE1E2int, 2)/std::pow(sigma_yEE1E2int, 3);
            break;
        case 291:
            dJydsigma[291] = 1.0/sigma_yEE1Ev3int - 1.0*std::pow(-myEE1Ev3int + yEE1Ev3int, 2)/std::pow(sigma_yEE1Ev3int, 3);
            break;
        case 292:
            dJydsigma[292] = 1.0/sigma_yEE1E1int - 1.0*std::pow(-myEE1E1int + yEE1E1int, 2)/std::pow(sigma_yEE1E1int, 3);
            break;
        case 293:
            dJydsigma[293] = 1.0/sigma_yEE1EE1int - 1.0*std::pow(-myEE1EE1int + yEE1EE1int, 2)/std::pow(sigma_yEE1EE1int, 3);
            break;
        case 294:
            dJydsigma[294] = 1.0/sigma_yEE1E3int - 1.0*std::pow(-myEE1E3int + yEE1E3int, 2)/std::pow(sigma_yEE1E3int, 3);
            break;
        case 295:
            dJydsigma[295] = 1.0/sigma_yEE1HE3int - 1.0*std::pow(-myEE1HE3int + yEE1HE3int, 2)/std::pow(sigma_yEE1HE3int, 3);
            break;
        case 296:
            dJydsigma[296] = 1.0/sigma_yEE1E4int - 1.0*std::pow(-myEE1E4int + yEE1E4int, 2)/std::pow(sigma_yEE1E4int, 3);
            break;
        case 297:
            dJydsigma[297] = 1.0/sigma_yEE1HE4int - 1.0*std::pow(-myEE1HE4int + yEE1HE4int, 2)/std::pow(sigma_yEE1HE4int, 3);
            break;
        case 298:
            dJydsigma[298] = 1.0/sigma_yE2HE3int - 1.0*std::pow(-myE2HE3int + yE2HE3int, 2)/std::pow(sigma_yE2HE3int, 3);
            break;
        case 299:
            dJydsigma[299] = 1.0/sigma_yHE3Ev3int - 1.0*std::pow(-myHE3Ev3int + yHE3Ev3int, 2)/std::pow(sigma_yHE3Ev3int, 3);
            break;
        case 300:
            dJydsigma[300] = 1.0/sigma_yE1HE3int - 1.0*std::pow(-myE1HE3int + yE1HE3int, 2)/std::pow(sigma_yE1HE3int, 3);
            break;
        case 301:
            dJydsigma[301] = 1.0/sigma_yHE3E4int - 1.0*std::pow(-myHE3E4int + yHE3E4int, 2)/std::pow(sigma_yHE3E4int, 3);
            break;
        case 302:
            dJydsigma[302] = 1.0/sigma_yHE3HE4int - 1.0*std::pow(-myHE3HE4int + yHE3HE4int, 2)/std::pow(sigma_yHE3HE4int, 3);
            break;
        case 303:
            dJydsigma[303] = 1.0/sigma_yE2HE4int - 1.0*std::pow(-myE2HE4int + yE2HE4int, 2)/std::pow(sigma_yE2HE4int, 3);
            break;
        case 304:
            dJydsigma[304] = 1.0/sigma_yHE4Ev3int - 1.0*std::pow(-myHE4Ev3int + yHE4Ev3int, 2)/std::pow(sigma_yHE4Ev3int, 3);
            break;
        case 305:
            dJydsigma[305] = 1.0/sigma_yE1HE4int - 1.0*std::pow(-myE1HE4int + yE1HE4int, 2)/std::pow(sigma_yE1HE4int, 3);
            break;
        case 306:
            dJydsigma[306] = 1.0/sigma_yE3HE4int - 1.0*std::pow(-myE3HE4int + yE3HE4int, 2)/std::pow(sigma_yE3HE4int, 3);
            break;
        case 307:
            dJydsigma[307] = 1.0/sigma_yHE4E4int - 1.0*std::pow(-myHE4E4int + yHE4E4int, 2)/std::pow(sigma_yHE4E4int, 3);
            break;
        case 308:
            dJydsigma[308] = 1.0/sigma_yHE4HE4int - 1.0*std::pow(-myHE4HE4int + yHE4HE4int, 2)/std::pow(sigma_yHE4HE4int, 3);
            break;
        case 309:
            dJydsigma[309] = 1.0/sigma_yHGF_Met_Metint - 1.0*std::pow(-myHGF_Met_Metint + yHGF_Met_Metint, 2)/std::pow(sigma_yHGF_Met_Metint, 3);
            break;
        case 310:
            dJydsigma[310] = 1.0/sigma_yHGF_Met_HGF_Metint - 1.0*std::pow(-myHGF_Met_HGF_Metint + yHGF_Met_HGF_Metint, 2)/std::pow(sigma_yHGF_Met_HGF_Metint, 3);
            break;
        case 311:
            dJydsigma[311] = 1.0/sigma_yPPrPPrint - 1.0*std::pow(-myPPrPPrint + yPPrPPrint, 2)/std::pow(sigma_yPPrPPrint, 3);
            break;
        case 312:
            dJydsigma[312] = 1.0/sigma_yPPrPrint - 1.0*std::pow(-myPPrPrint + yPPrPrint, 2)/std::pow(sigma_yPPrPrint, 3);
            break;
        case 313:
            dJydsigma[313] = 1.0/sigma_yFFrFFrint - 1.0*std::pow(-myFFrFFrint + yFFrFFrint, 2)/std::pow(sigma_yFFrFFrint, 3);
            break;
        case 314:
            dJydsigma[314] = 1.0/sigma_yFFrFrint - 1.0*std::pow(-myFFrFrint + yFFrFrint, 2)/std::pow(sigma_yFFrFrint, 3);
            break;
        case 315:
            dJydsigma[315] = 1.0/sigma_yIIrIr_int - 1.0*std::pow(-myIIrIr_int + yIIrIr_int, 2)/std::pow(sigma_yIIrIr_int, 3);
            break;
        case 316:
            dJydsigma[316] = 1.0/sigma_yINS_Isr_Isr_int - 1.0*std::pow(-myINS_Isr_Isr_int + yINS_Isr_Isr_int, 2)/std::pow(sigma_yINS_Isr_Isr_int, 3);
            break;
        case 317:
            dJydsigma[317] = 1.0/sigma_yIIrIrI_int - 1.0*std::pow(-myIIrIrI_int + yIIrIrI_int, 2)/std::pow(sigma_yIIrIrI_int, 3);
            break;
        case 318:
            dJydsigma[318] = 1.0/sigma_yINS_Isr_Isr_INS_int - 1.0*std::pow(-myINS_Isr_Isr_INS_int + yINS_Isr_Isr_INS_int, 2)/std::pow(sigma_yINS_Isr_Isr_INS_int, 3);
            break;
        case 319:
            dJydsigma[319] = 1.0/sigma_ypEE1E2int - 1.0*std::pow(-mypEE1E2int + ypEE1E2int, 2)/std::pow(sigma_ypEE1E2int, 3);
            break;
        case 320:
            dJydsigma[320] = 1.0/sigma_ypEE1Ev3int - 1.0*std::pow(-mypEE1Ev3int + ypEE1Ev3int, 2)/std::pow(sigma_ypEE1Ev3int, 3);
            break;
        case 321:
            dJydsigma[321] = 1.0/sigma_ypEE1E1int - 1.0*std::pow(-mypEE1E1int + ypEE1E1int, 2)/std::pow(sigma_ypEE1E1int, 3);
            break;
        case 322:
            dJydsigma[322] = 1.0/sigma_ypEE1EE1int - 1.0*std::pow(-mypEE1EE1int + ypEE1EE1int, 2)/std::pow(sigma_ypEE1EE1int, 3);
            break;
        case 323:
            dJydsigma[323] = 1.0/sigma_ypEE1E3int - 1.0*std::pow(-mypEE1E3int + ypEE1E3int, 2)/std::pow(sigma_ypEE1E3int, 3);
            break;
        case 324:
            dJydsigma[324] = 1.0/sigma_ypEE1HE3int - 1.0*std::pow(-mypEE1HE3int + ypEE1HE3int, 2)/std::pow(sigma_ypEE1HE3int, 3);
            break;
        case 325:
            dJydsigma[325] = 1.0/sigma_ypEE1E4int - 1.0*std::pow(-mypEE1E4int + ypEE1E4int, 2)/std::pow(sigma_ypEE1E4int, 3);
            break;
        case 326:
            dJydsigma[326] = 1.0/sigma_ypEE1HE4int - 1.0*std::pow(-mypEE1HE4int + ypEE1HE4int, 2)/std::pow(sigma_ypEE1HE4int, 3);
            break;
        case 327:
            dJydsigma[327] = 1.0/sigma_ypE2HE3int - 1.0*std::pow(-mypE2HE3int + ypE2HE3int, 2)/std::pow(sigma_ypE2HE3int, 3);
            break;
        case 328:
            dJydsigma[328] = 1.0/sigma_ypHE3Ev3int - 1.0*std::pow(-mypHE3Ev3int + ypHE3Ev3int, 2)/std::pow(sigma_ypHE3Ev3int, 3);
            break;
        case 329:
            dJydsigma[329] = 1.0/sigma_ypE1HE3int - 1.0*std::pow(-mypE1HE3int + ypE1HE3int, 2)/std::pow(sigma_ypE1HE3int, 3);
            break;
        case 330:
            dJydsigma[330] = 1.0/sigma_ypHE3E4int - 1.0*std::pow(-mypHE3E4int + ypHE3E4int, 2)/std::pow(sigma_ypHE3E4int, 3);
            break;
        case 331:
            dJydsigma[331] = 1.0/sigma_ypHE3HE4int - 1.0*std::pow(-mypHE3HE4int + ypHE3HE4int, 2)/std::pow(sigma_ypHE3HE4int, 3);
            break;
        case 332:
            dJydsigma[332] = 1.0/sigma_ypE2HE4int - 1.0*std::pow(-mypE2HE4int + ypE2HE4int, 2)/std::pow(sigma_ypE2HE4int, 3);
            break;
        case 333:
            dJydsigma[333] = 1.0/sigma_ypHE4Ev3int - 1.0*std::pow(-mypHE4Ev3int + ypHE4Ev3int, 2)/std::pow(sigma_ypHE4Ev3int, 3);
            break;
        case 334:
            dJydsigma[334] = 1.0/sigma_ypE1HE4int - 1.0*std::pow(-mypE1HE4int + ypE1HE4int, 2)/std::pow(sigma_ypE1HE4int, 3);
            break;
        case 335:
            dJydsigma[335] = 1.0/sigma_ypE3HE4int - 1.0*std::pow(-mypE3HE4int + ypE3HE4int, 2)/std::pow(sigma_ypE3HE4int, 3);
            break;
        case 336:
            dJydsigma[336] = 1.0/sigma_ypHE4E4int - 1.0*std::pow(-mypHE4E4int + ypHE4E4int, 2)/std::pow(sigma_ypHE4E4int, 3);
            break;
        case 337:
            dJydsigma[337] = 1.0/sigma_ypHE4HE4int - 1.0*std::pow(-mypHE4HE4int + ypHE4HE4int, 2)/std::pow(sigma_ypHE4HE4int, 3);
            break;
        case 338:
            dJydsigma[338] = 1.0/sigma_ypHGF_Met_Metint - 1.0*std::pow(-mypHGF_Met_Metint + ypHGF_Met_Metint, 2)/std::pow(sigma_ypHGF_Met_Metint, 3);
            break;
        case 339:
            dJydsigma[339] = 1.0/sigma_ypHGF_Met_HGF_Metint - 1.0*std::pow(-mypHGF_Met_HGF_Metint + ypHGF_Met_HGF_Metint, 2)/std::pow(sigma_ypHGF_Met_HGF_Metint, 3);
            break;
        case 340:
            dJydsigma[340] = 1.0/sigma_ypPPrPPrint - 1.0*std::pow(-mypPPrPPrint + ypPPrPPrint, 2)/std::pow(sigma_ypPPrPPrint, 3);
            break;
        case 341:
            dJydsigma[341] = 1.0/sigma_ypPPrPrint - 1.0*std::pow(-mypPPrPrint + ypPPrPrint, 2)/std::pow(sigma_ypPPrPrint, 3);
            break;
        case 342:
            dJydsigma[342] = 1.0/sigma_ypFFrFFrint - 1.0*std::pow(-mypFFrFFrint + ypFFrFFrint, 2)/std::pow(sigma_ypFFrFFrint, 3);
            break;
        case 343:
            dJydsigma[343] = 1.0/sigma_ypFFrFrint - 1.0*std::pow(-mypFFrFrint + ypFFrFrint, 2)/std::pow(sigma_ypFFrFrint, 3);
            break;
        case 344:
            dJydsigma[344] = 1.0/sigma_ypIIrIr_int - 1.0*std::pow(-mypIIrIr_int + ypIIrIr_int, 2)/std::pow(sigma_ypIIrIr_int, 3);
            break;
        case 345:
            dJydsigma[345] = 1.0/sigma_ypINS_Isr_Isr_int - 1.0*std::pow(-mypINS_Isr_Isr_int + ypINS_Isr_Isr_int, 2)/std::pow(sigma_ypINS_Isr_Isr_int, 3);
            break;
        case 346:
            dJydsigma[346] = 1.0/sigma_ypIIrIrI_int - 1.0*std::pow(-mypIIrIrI_int + ypIIrIrI_int, 2)/std::pow(sigma_ypIIrIrI_int, 3);
            break;
        case 347:
            dJydsigma[347] = 1.0/sigma_ypINS_Isr_Isr_INS_int - 1.0*std::pow(-mypINS_Isr_Isr_INS_int + ypINS_Isr_Isr_INS_int, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int, 3);
            break;
        case 348:
            dJydsigma[348] = 1.0/sigma_ypIIrIr_int_IRS - 1.0*std::pow(-mypIIrIr_int_IRS + ypIIrIr_int_IRS, 2)/std::pow(sigma_ypIIrIr_int_IRS, 3);
            break;
        case 349:
            dJydsigma[349] = 1.0/sigma_ypINS_Isr_Isr_int_IRS - 1.0*std::pow(-mypINS_Isr_Isr_int_IRS + ypINS_Isr_Isr_int_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_int_IRS, 3);
            break;
        case 350:
            dJydsigma[350] = 1.0/sigma_ypIIrIrI_int_IRS - 1.0*std::pow(-mypIIrIrI_int_IRS + ypIIrIrI_int_IRS, 2)/std::pow(sigma_ypIIrIrI_int_IRS, 3);
            break;
        case 351:
            dJydsigma[351] = 1.0/sigma_ypINS_Isr_Isr_INS_int_IRS - 1.0*std::pow(-mypINS_Isr_Isr_INS_int_IRS + ypINS_Isr_Isr_INS_int_IRS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS, 3);
            break;
        case 352:
            dJydsigma[352] = 1.0/sigma_ypEE1E2_G2_SOS - 1.0*std::pow(-mypEE1E2_G2_SOS + ypEE1E2_G2_SOS, 2)/std::pow(sigma_ypEE1E2_G2_SOS, 3);
            break;
        case 353:
            dJydsigma[353] = 1.0/sigma_ypEE1Ev3_G2_SOS - 1.0*std::pow(-mypEE1Ev3_G2_SOS + ypEE1Ev3_G2_SOS, 2)/std::pow(sigma_ypEE1Ev3_G2_SOS, 3);
            break;
        case 354:
            dJydsigma[354] = 1.0/sigma_ypEE1E1_G2_SOS - 1.0*std::pow(-mypEE1E1_G2_SOS + ypEE1E1_G2_SOS, 2)/std::pow(sigma_ypEE1E1_G2_SOS, 3);
            break;
        case 355:
            dJydsigma[355] = 1.0/sigma_ypEE1EE1_G2_SOS - 1.0*std::pow(-mypEE1EE1_G2_SOS + ypEE1EE1_G2_SOS, 2)/std::pow(sigma_ypEE1EE1_G2_SOS, 3);
            break;
        case 356:
            dJydsigma[356] = 1.0/sigma_ypEE1E3_G2_SOS - 1.0*std::pow(-mypEE1E3_G2_SOS + ypEE1E3_G2_SOS, 2)/std::pow(sigma_ypEE1E3_G2_SOS, 3);
            break;
        case 357:
            dJydsigma[357] = 1.0/sigma_ypEE1HE3_G2_SOS - 1.0*std::pow(-mypEE1HE3_G2_SOS + ypEE1HE3_G2_SOS, 2)/std::pow(sigma_ypEE1HE3_G2_SOS, 3);
            break;
        case 358:
            dJydsigma[358] = 1.0/sigma_ypEE1E4_G2_SOS - 1.0*std::pow(-mypEE1E4_G2_SOS + ypEE1E4_G2_SOS, 2)/std::pow(sigma_ypEE1E4_G2_SOS, 3);
            break;
        case 359:
            dJydsigma[359] = 1.0/sigma_ypEE1HE4_G2_SOS - 1.0*std::pow(-mypEE1HE4_G2_SOS + ypEE1HE4_G2_SOS, 2)/std::pow(sigma_ypEE1HE4_G2_SOS, 3);
            break;
        case 360:
            dJydsigma[360] = 1.0/sigma_ypE2HE3_G2_SOS - 1.0*std::pow(-mypE2HE3_G2_SOS + ypE2HE3_G2_SOS, 2)/std::pow(sigma_ypE2HE3_G2_SOS, 3);
            break;
        case 361:
            dJydsigma[361] = 1.0/sigma_ypHE3Ev3_G2_SOS - 1.0*std::pow(-mypHE3Ev3_G2_SOS + ypHE3Ev3_G2_SOS, 2)/std::pow(sigma_ypHE3Ev3_G2_SOS, 3);
            break;
        case 362:
            dJydsigma[362] = 1.0/sigma_ypE1HE3_G2_SOS - 1.0*std::pow(-mypE1HE3_G2_SOS + ypE1HE3_G2_SOS, 2)/std::pow(sigma_ypE1HE3_G2_SOS, 3);
            break;
        case 363:
            dJydsigma[363] = 1.0/sigma_ypHE3E4_G2_SOS - 1.0*std::pow(-mypHE3E4_G2_SOS + ypHE3E4_G2_SOS, 2)/std::pow(sigma_ypHE3E4_G2_SOS, 3);
            break;
        case 364:
            dJydsigma[364] = 1.0/sigma_ypHE3HE4_G2_SOS - 1.0*std::pow(-mypHE3HE4_G2_SOS + ypHE3HE4_G2_SOS, 2)/std::pow(sigma_ypHE3HE4_G2_SOS, 3);
            break;
        case 365:
            dJydsigma[365] = 1.0/sigma_ypE2HE4_G2_SOS - 1.0*std::pow(-mypE2HE4_G2_SOS + ypE2HE4_G2_SOS, 2)/std::pow(sigma_ypE2HE4_G2_SOS, 3);
            break;
        case 366:
            dJydsigma[366] = 1.0/sigma_ypHE4Ev3_G2_SOS - 1.0*std::pow(-mypHE4Ev3_G2_SOS + ypHE4Ev3_G2_SOS, 2)/std::pow(sigma_ypHE4Ev3_G2_SOS, 3);
            break;
        case 367:
            dJydsigma[367] = 1.0/sigma_ypE1HE4_G2_SOS - 1.0*std::pow(-mypE1HE4_G2_SOS + ypE1HE4_G2_SOS, 2)/std::pow(sigma_ypE1HE4_G2_SOS, 3);
            break;
        case 368:
            dJydsigma[368] = 1.0/sigma_ypE3HE4_G2_SOS - 1.0*std::pow(-mypE3HE4_G2_SOS + ypE3HE4_G2_SOS, 2)/std::pow(sigma_ypE3HE4_G2_SOS, 3);
            break;
        case 369:
            dJydsigma[369] = 1.0/sigma_ypHE4E4_G2_SOS - 1.0*std::pow(-mypHE4E4_G2_SOS + ypHE4E4_G2_SOS, 2)/std::pow(sigma_ypHE4E4_G2_SOS, 3);
            break;
        case 370:
            dJydsigma[370] = 1.0/sigma_ypHE4HE4_G2_SOS - 1.0*std::pow(-mypHE4HE4_G2_SOS + ypHE4HE4_G2_SOS, 2)/std::pow(sigma_ypHE4HE4_G2_SOS, 3);
            break;
        case 371:
            dJydsigma[371] = 1.0/sigma_ypHGF_Met_Met_G2_SOS - 1.0*std::pow(-mypHGF_Met_Met_G2_SOS + ypHGF_Met_Met_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_Met_G2_SOS, 3);
            break;
        case 372:
            dJydsigma[372] = 1.0/sigma_ypHGF_Met_HGF_Met_G2_SOS - 1.0*std::pow(-mypHGF_Met_HGF_Met_G2_SOS + ypHGF_Met_HGF_Met_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS, 3);
            break;
        case 373:
            dJydsigma[373] = 1.0/sigma_ypPPrPPr_G2_SOS - 1.0*std::pow(-mypPPrPPr_G2_SOS + ypPPrPPr_G2_SOS, 2)/std::pow(sigma_ypPPrPPr_G2_SOS, 3);
            break;
        case 374:
            dJydsigma[374] = 1.0/sigma_ypPPrPr_G2_SOS - 1.0*std::pow(-mypPPrPr_G2_SOS + ypPPrPr_G2_SOS, 2)/std::pow(sigma_ypPPrPr_G2_SOS, 3);
            break;
        case 375:
            dJydsigma[375] = 1.0/sigma_ypFFrFFr_G2_SOS - 1.0*std::pow(-mypFFrFFr_G2_SOS + ypFFrFFr_G2_SOS, 2)/std::pow(sigma_ypFFrFFr_G2_SOS, 3);
            break;
        case 376:
            dJydsigma[376] = 1.0/sigma_ypFFrFr_G2_SOS - 1.0*std::pow(-mypFFrFr_G2_SOS + ypFFrFr_G2_SOS, 2)/std::pow(sigma_ypFFrFr_G2_SOS, 3);
            break;
        case 377:
            dJydsigma[377] = 1.0/sigma_ypIIrIr_IRS_G2_SOS - 1.0*std::pow(-mypIIrIr_IRS_G2_SOS + ypIIrIr_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIr_IRS_G2_SOS, 3);
            break;
        case 378:
            dJydsigma[378] = 1.0/sigma_ypINS_Isr_Isr_IRS_G2_SOS - 1.0*std::pow(-mypINS_Isr_Isr_IRS_G2_SOS + ypINS_Isr_Isr_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS, 3);
            break;
        case 379:
            dJydsigma[379] = 1.0/sigma_ypIIrIrI_IRS_G2_SOS - 1.0*std::pow(-mypIIrIrI_IRS_G2_SOS + ypIIrIrI_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIrI_IRS_G2_SOS, 3);
            break;
        case 380:
            dJydsigma[380] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_G2_SOS + ypINS_Isr_Isr_INS_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS, 3);
            break;
        case 381:
            dJydsigma[381] = 1.0/sigma_ypEE1E2int_G2_SOS - 1.0*std::pow(-mypEE1E2int_G2_SOS + ypEE1E2int_G2_SOS, 2)/std::pow(sigma_ypEE1E2int_G2_SOS, 3);
            break;
        case 382:
            dJydsigma[382] = 1.0/sigma_ypEE1Ev3int_G2_SOS - 1.0*std::pow(-mypEE1Ev3int_G2_SOS + ypEE1Ev3int_G2_SOS, 2)/std::pow(sigma_ypEE1Ev3int_G2_SOS, 3);
            break;
        case 383:
            dJydsigma[383] = 1.0/sigma_ypEE1E1int_G2_SOS - 1.0*std::pow(-mypEE1E1int_G2_SOS + ypEE1E1int_G2_SOS, 2)/std::pow(sigma_ypEE1E1int_G2_SOS, 3);
            break;
        case 384:
            dJydsigma[384] = 1.0/sigma_ypEE1EE1int_G2_SOS - 1.0*std::pow(-mypEE1EE1int_G2_SOS + ypEE1EE1int_G2_SOS, 2)/std::pow(sigma_ypEE1EE1int_G2_SOS, 3);
            break;
        case 385:
            dJydsigma[385] = 1.0/sigma_ypEE1E3int_G2_SOS - 1.0*std::pow(-mypEE1E3int_G2_SOS + ypEE1E3int_G2_SOS, 2)/std::pow(sigma_ypEE1E3int_G2_SOS, 3);
            break;
        case 386:
            dJydsigma[386] = 1.0/sigma_ypEE1HE3int_G2_SOS - 1.0*std::pow(-mypEE1HE3int_G2_SOS + ypEE1HE3int_G2_SOS, 2)/std::pow(sigma_ypEE1HE3int_G2_SOS, 3);
            break;
        case 387:
            dJydsigma[387] = 1.0/sigma_ypEE1E4int_G2_SOS - 1.0*std::pow(-mypEE1E4int_G2_SOS + ypEE1E4int_G2_SOS, 2)/std::pow(sigma_ypEE1E4int_G2_SOS, 3);
            break;
        case 388:
            dJydsigma[388] = 1.0/sigma_ypEE1HE4int_G2_SOS - 1.0*std::pow(-mypEE1HE4int_G2_SOS + ypEE1HE4int_G2_SOS, 2)/std::pow(sigma_ypEE1HE4int_G2_SOS, 3);
            break;
        case 389:
            dJydsigma[389] = 1.0/sigma_ypE2HE3int_G2_SOS - 1.0*std::pow(-mypE2HE3int_G2_SOS + ypE2HE3int_G2_SOS, 2)/std::pow(sigma_ypE2HE3int_G2_SOS, 3);
            break;
        case 390:
            dJydsigma[390] = 1.0/sigma_ypHE3Ev3int_G2_SOS - 1.0*std::pow(-mypHE3Ev3int_G2_SOS + ypHE3Ev3int_G2_SOS, 2)/std::pow(sigma_ypHE3Ev3int_G2_SOS, 3);
            break;
        case 391:
            dJydsigma[391] = 1.0/sigma_ypE1HE3int_G2_SOS - 1.0*std::pow(-mypE1HE3int_G2_SOS + ypE1HE3int_G2_SOS, 2)/std::pow(sigma_ypE1HE3int_G2_SOS, 3);
            break;
        case 392:
            dJydsigma[392] = 1.0/sigma_ypHE3E4int_G2_SOS - 1.0*std::pow(-mypHE3E4int_G2_SOS + ypHE3E4int_G2_SOS, 2)/std::pow(sigma_ypHE3E4int_G2_SOS, 3);
            break;
        case 393:
            dJydsigma[393] = 1.0/sigma_ypHE3HE4int_G2_SOS - 1.0*std::pow(-mypHE3HE4int_G2_SOS + ypHE3HE4int_G2_SOS, 2)/std::pow(sigma_ypHE3HE4int_G2_SOS, 3);
            break;
        case 394:
            dJydsigma[394] = 1.0/sigma_ypE2HE4int_G2_SOS - 1.0*std::pow(-mypE2HE4int_G2_SOS + ypE2HE4int_G2_SOS, 2)/std::pow(sigma_ypE2HE4int_G2_SOS, 3);
            break;
        case 395:
            dJydsigma[395] = 1.0/sigma_ypHE4Ev3int_G2_SOS - 1.0*std::pow(-mypHE4Ev3int_G2_SOS + ypHE4Ev3int_G2_SOS, 2)/std::pow(sigma_ypHE4Ev3int_G2_SOS, 3);
            break;
        case 396:
            dJydsigma[396] = 1.0/sigma_ypE1HE4int_G2_SOS - 1.0*std::pow(-mypE1HE4int_G2_SOS + ypE1HE4int_G2_SOS, 2)/std::pow(sigma_ypE1HE4int_G2_SOS, 3);
            break;
        case 397:
            dJydsigma[397] = 1.0/sigma_ypE3HE4int_G2_SOS - 1.0*std::pow(-mypE3HE4int_G2_SOS + ypE3HE4int_G2_SOS, 2)/std::pow(sigma_ypE3HE4int_G2_SOS, 3);
            break;
        case 398:
            dJydsigma[398] = 1.0/sigma_ypHE4E4int_G2_SOS - 1.0*std::pow(-mypHE4E4int_G2_SOS + ypHE4E4int_G2_SOS, 2)/std::pow(sigma_ypHE4E4int_G2_SOS, 3);
            break;
        case 399:
            dJydsigma[399] = 1.0/sigma_ypHE4HE4int_G2_SOS - 1.0*std::pow(-mypHE4HE4int_G2_SOS + ypHE4HE4int_G2_SOS, 2)/std::pow(sigma_ypHE4HE4int_G2_SOS, 3);
            break;
        case 400:
            dJydsigma[400] = 1.0/sigma_ypHGF_Met_Metint_G2_SOS - 1.0*std::pow(-mypHGF_Met_Metint_G2_SOS + ypHGF_Met_Metint_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_Metint_G2_SOS, 3);
            break;
        case 401:
            dJydsigma[401] = 1.0/sigma_ypHGF_Met_HGF_Metint_G2_SOS - 1.0*std::pow(-mypHGF_Met_HGF_Metint_G2_SOS + ypHGF_Met_HGF_Metint_G2_SOS, 2)/std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS, 3);
            break;
        case 402:
            dJydsigma[402] = 1.0/sigma_ypPPrPPrint_G2_SOS - 1.0*std::pow(-mypPPrPPrint_G2_SOS + ypPPrPPrint_G2_SOS, 2)/std::pow(sigma_ypPPrPPrint_G2_SOS, 3);
            break;
        case 403:
            dJydsigma[403] = 1.0/sigma_ypPPrPrint_G2_SOS - 1.0*std::pow(-mypPPrPrint_G2_SOS + ypPPrPrint_G2_SOS, 2)/std::pow(sigma_ypPPrPrint_G2_SOS, 3);
            break;
        case 404:
            dJydsigma[404] = 1.0/sigma_ypFFrFFrint_G2_SOS - 1.0*std::pow(-mypFFrFFrint_G2_SOS + ypFFrFFrint_G2_SOS, 2)/std::pow(sigma_ypFFrFFrint_G2_SOS, 3);
            break;
        case 405:
            dJydsigma[405] = 1.0/sigma_ypFFrFrint_G2_SOS - 1.0*std::pow(-mypFFrFrint_G2_SOS + ypFFrFrint_G2_SOS, 2)/std::pow(sigma_ypFFrFrint_G2_SOS, 3);
            break;
        case 406:
            dJydsigma[406] = 1.0/sigma_ypIIrIr_int_IRS_G2_SOS - 1.0*std::pow(-mypIIrIr_int_IRS_G2_SOS + ypIIrIr_int_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIr_int_IRS_G2_SOS, 3);
            break;
        case 407:
            dJydsigma[407] = 1.0/sigma_ypINS_Isr_Isr_int_IRS_G2_SOS - 1.0*std::pow(-mypINS_Isr_Isr_int_IRS_G2_SOS + ypINS_Isr_Isr_int_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS, 3);
            break;
        case 408:
            dJydsigma[408] = 1.0/sigma_ypIIrIrI_int_IRS_G2_SOS - 1.0*std::pow(-mypIIrIrI_int_IRS_G2_SOS + ypIIrIrI_int_IRS_G2_SOS, 2)/std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS, 3);
            break;
        case 409:
            dJydsigma[409] = 1.0/sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS - 1.0*std::pow(-mypINS_Isr_Isr_INS_int_IRS_G2_SOS + ypINS_Isr_Isr_INS_int_IRS_G2_SOS, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS, 3);
            break;
        case 410:
            dJydsigma[410] = 1.0/sigma_ypEE1E2_PLCg - 1.0*std::pow(-mypEE1E2_PLCg + ypEE1E2_PLCg, 2)/std::pow(sigma_ypEE1E2_PLCg, 3);
            break;
        case 411:
            dJydsigma[411] = 1.0/sigma_ypEE1Ev3_PLCg - 1.0*std::pow(-mypEE1Ev3_PLCg + ypEE1Ev3_PLCg, 2)/std::pow(sigma_ypEE1Ev3_PLCg, 3);
            break;
        case 412:
            dJydsigma[412] = 1.0/sigma_ypEE1E1_PLCg - 1.0*std::pow(-mypEE1E1_PLCg + ypEE1E1_PLCg, 2)/std::pow(sigma_ypEE1E1_PLCg, 3);
            break;
        case 413:
            dJydsigma[413] = 1.0/sigma_ypEE1EE1_PLCg - 1.0*std::pow(-mypEE1EE1_PLCg + ypEE1EE1_PLCg, 2)/std::pow(sigma_ypEE1EE1_PLCg, 3);
            break;
        case 414:
            dJydsigma[414] = 1.0/sigma_ypEE1E3_PLCg - 1.0*std::pow(-mypEE1E3_PLCg + ypEE1E3_PLCg, 2)/std::pow(sigma_ypEE1E3_PLCg, 3);
            break;
        case 415:
            dJydsigma[415] = 1.0/sigma_ypEE1HE3_PLCg - 1.0*std::pow(-mypEE1HE3_PLCg + ypEE1HE3_PLCg, 2)/std::pow(sigma_ypEE1HE3_PLCg, 3);
            break;
        case 416:
            dJydsigma[416] = 1.0/sigma_ypEE1E4_PLCg - 1.0*std::pow(-mypEE1E4_PLCg + ypEE1E4_PLCg, 2)/std::pow(sigma_ypEE1E4_PLCg, 3);
            break;
        case 417:
            dJydsigma[417] = 1.0/sigma_ypEE1HE4_PLCg - 1.0*std::pow(-mypEE1HE4_PLCg + ypEE1HE4_PLCg, 2)/std::pow(sigma_ypEE1HE4_PLCg, 3);
            break;
        case 418:
            dJydsigma[418] = 1.0/sigma_ypE2HE3_PLCg - 1.0*std::pow(-mypE2HE3_PLCg + ypE2HE3_PLCg, 2)/std::pow(sigma_ypE2HE3_PLCg, 3);
            break;
        case 419:
            dJydsigma[419] = 1.0/sigma_ypHE3Ev3_PLCg - 1.0*std::pow(-mypHE3Ev3_PLCg + ypHE3Ev3_PLCg, 2)/std::pow(sigma_ypHE3Ev3_PLCg, 3);
            break;
        case 420:
            dJydsigma[420] = 1.0/sigma_ypE1HE3_PLCg - 1.0*std::pow(-mypE1HE3_PLCg + ypE1HE3_PLCg, 2)/std::pow(sigma_ypE1HE3_PLCg, 3);
            break;
        case 421:
            dJydsigma[421] = 1.0/sigma_ypHE3E4_PLCg - 1.0*std::pow(-mypHE3E4_PLCg + ypHE3E4_PLCg, 2)/std::pow(sigma_ypHE3E4_PLCg, 3);
            break;
        case 422:
            dJydsigma[422] = 1.0/sigma_ypHE3HE4_PLCg - 1.0*std::pow(-mypHE3HE4_PLCg + ypHE3HE4_PLCg, 2)/std::pow(sigma_ypHE3HE4_PLCg, 3);
            break;
        case 423:
            dJydsigma[423] = 1.0/sigma_ypE2HE4_PLCg - 1.0*std::pow(-mypE2HE4_PLCg + ypE2HE4_PLCg, 2)/std::pow(sigma_ypE2HE4_PLCg, 3);
            break;
        case 424:
            dJydsigma[424] = 1.0/sigma_ypHE4Ev3_PLCg - 1.0*std::pow(-mypHE4Ev3_PLCg + ypHE4Ev3_PLCg, 2)/std::pow(sigma_ypHE4Ev3_PLCg, 3);
            break;
        case 425:
            dJydsigma[425] = 1.0/sigma_ypE1HE4_PLCg - 1.0*std::pow(-mypE1HE4_PLCg + ypE1HE4_PLCg, 2)/std::pow(sigma_ypE1HE4_PLCg, 3);
            break;
        case 426:
            dJydsigma[426] = 1.0/sigma_ypE3HE4_PLCg - 1.0*std::pow(-mypE3HE4_PLCg + ypE3HE4_PLCg, 2)/std::pow(sigma_ypE3HE4_PLCg, 3);
            break;
        case 427:
            dJydsigma[427] = 1.0/sigma_ypHE4E4_PLCg - 1.0*std::pow(-mypHE4E4_PLCg + ypHE4E4_PLCg, 2)/std::pow(sigma_ypHE4E4_PLCg, 3);
            break;
        case 428:
            dJydsigma[428] = 1.0/sigma_ypHE4HE4_PLCg - 1.0*std::pow(-mypHE4HE4_PLCg + ypHE4HE4_PLCg, 2)/std::pow(sigma_ypHE4HE4_PLCg, 3);
            break;
        case 429:
            dJydsigma[429] = 1.0/sigma_ypHGF_Met_Met_PLCg - 1.0*std::pow(-mypHGF_Met_Met_PLCg + ypHGF_Met_Met_PLCg, 2)/std::pow(sigma_ypHGF_Met_Met_PLCg, 3);
            break;
        case 430:
            dJydsigma[430] = 1.0/sigma_ypHGF_Met_HGF_Met_PLCg - 1.0*std::pow(-mypHGF_Met_HGF_Met_PLCg + ypHGF_Met_HGF_Met_PLCg, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PLCg, 3);
            break;
        case 431:
            dJydsigma[431] = 1.0/sigma_ypPPrPPr_PLCg - 1.0*std::pow(-mypPPrPPr_PLCg + ypPPrPPr_PLCg, 2)/std::pow(sigma_ypPPrPPr_PLCg, 3);
            break;
        case 432:
            dJydsigma[432] = 1.0/sigma_ypPPrPr_PLCg - 1.0*std::pow(-mypPPrPr_PLCg + ypPPrPr_PLCg, 2)/std::pow(sigma_ypPPrPr_PLCg, 3);
            break;
        case 433:
            dJydsigma[433] = 1.0/sigma_ypFFrFFr_PLCg - 1.0*std::pow(-mypFFrFFr_PLCg + ypFFrFFr_PLCg, 2)/std::pow(sigma_ypFFrFFr_PLCg, 3);
            break;
        case 434:
            dJydsigma[434] = 1.0/sigma_ypFFrFr_PLCg - 1.0*std::pow(-mypFFrFr_PLCg + ypFFrFr_PLCg, 2)/std::pow(sigma_ypFFrFr_PLCg, 3);
            break;
        case 435:
            dJydsigma[435] = 1.0/sigma_ypIIrIr_IRS_PLCg - 1.0*std::pow(-mypIIrIr_IRS_PLCg + ypIIrIr_IRS_PLCg, 2)/std::pow(sigma_ypIIrIr_IRS_PLCg, 3);
            break;
        case 436:
            dJydsigma[436] = 1.0/sigma_ypINS_Isr_Isr_IRS_PLCg - 1.0*std::pow(-mypINS_Isr_Isr_IRS_PLCg + ypINS_Isr_Isr_IRS_PLCg, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg, 3);
            break;
        case 437:
            dJydsigma[437] = 1.0/sigma_ypIIrIrI_IRS_PLCg - 1.0*std::pow(-mypIIrIrI_IRS_PLCg + ypIIrIrI_IRS_PLCg, 2)/std::pow(sigma_ypIIrIrI_IRS_PLCg, 3);
            break;
        case 438:
            dJydsigma[438] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_PLCg - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_PLCg + ypINS_Isr_Isr_INS_IRS_PLCg, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg, 3);
            break;
        case 439:
            dJydsigma[439] = 1.0/sigma_ypEE1E2_PI3K1 - 1.0*std::pow(-mypEE1E2_PI3K1 + ypEE1E2_PI3K1, 2)/std::pow(sigma_ypEE1E2_PI3K1, 3);
            break;
        case 440:
            dJydsigma[440] = 1.0/sigma_ypEE1Ev3_PI3K1 - 1.0*std::pow(-mypEE1Ev3_PI3K1 + ypEE1Ev3_PI3K1, 2)/std::pow(sigma_ypEE1Ev3_PI3K1, 3);
            break;
        case 441:
            dJydsigma[441] = 1.0/sigma_ypEE1E1_PI3K1 - 1.0*std::pow(-mypEE1E1_PI3K1 + ypEE1E1_PI3K1, 2)/std::pow(sigma_ypEE1E1_PI3K1, 3);
            break;
        case 442:
            dJydsigma[442] = 1.0/sigma_ypEE1EE1_PI3K1 - 1.0*std::pow(-mypEE1EE1_PI3K1 + ypEE1EE1_PI3K1, 2)/std::pow(sigma_ypEE1EE1_PI3K1, 3);
            break;
        case 443:
            dJydsigma[443] = 1.0/sigma_ypEE1E3_PI3K1 - 1.0*std::pow(-mypEE1E3_PI3K1 + ypEE1E3_PI3K1, 2)/std::pow(sigma_ypEE1E3_PI3K1, 3);
            break;
        case 444:
            dJydsigma[444] = 1.0/sigma_ypEE1HE3_PI3K1 - 1.0*std::pow(-mypEE1HE3_PI3K1 + ypEE1HE3_PI3K1, 2)/std::pow(sigma_ypEE1HE3_PI3K1, 3);
            break;
        case 445:
            dJydsigma[445] = 1.0/sigma_ypEE1E4_PI3K1 - 1.0*std::pow(-mypEE1E4_PI3K1 + ypEE1E4_PI3K1, 2)/std::pow(sigma_ypEE1E4_PI3K1, 3);
            break;
        case 446:
            dJydsigma[446] = 1.0/sigma_ypEE1HE4_PI3K1 - 1.0*std::pow(-mypEE1HE4_PI3K1 + ypEE1HE4_PI3K1, 2)/std::pow(sigma_ypEE1HE4_PI3K1, 3);
            break;
        case 447:
            dJydsigma[447] = 1.0/sigma_ypE2HE3_PI3K1 - 1.0*std::pow(-mypE2HE3_PI3K1 + ypE2HE3_PI3K1, 2)/std::pow(sigma_ypE2HE3_PI3K1, 3);
            break;
        case 448:
            dJydsigma[448] = 1.0/sigma_ypHE3Ev3_PI3K1 - 1.0*std::pow(-mypHE3Ev3_PI3K1 + ypHE3Ev3_PI3K1, 2)/std::pow(sigma_ypHE3Ev3_PI3K1, 3);
            break;
        case 449:
            dJydsigma[449] = 1.0/sigma_ypE1HE3_PI3K1 - 1.0*std::pow(-mypE1HE3_PI3K1 + ypE1HE3_PI3K1, 2)/std::pow(sigma_ypE1HE3_PI3K1, 3);
            break;
        case 450:
            dJydsigma[450] = 1.0/sigma_ypHE3E4_PI3K1 - 1.0*std::pow(-mypHE3E4_PI3K1 + ypHE3E4_PI3K1, 2)/std::pow(sigma_ypHE3E4_PI3K1, 3);
            break;
        case 451:
            dJydsigma[451] = 1.0/sigma_ypHE3HE4_PI3K1 - 1.0*std::pow(-mypHE3HE4_PI3K1 + ypHE3HE4_PI3K1, 2)/std::pow(sigma_ypHE3HE4_PI3K1, 3);
            break;
        case 452:
            dJydsigma[452] = 1.0/sigma_ypE2HE4_PI3K1 - 1.0*std::pow(-mypE2HE4_PI3K1 + ypE2HE4_PI3K1, 2)/std::pow(sigma_ypE2HE4_PI3K1, 3);
            break;
        case 453:
            dJydsigma[453] = 1.0/sigma_ypHE4Ev3_PI3K1 - 1.0*std::pow(-mypHE4Ev3_PI3K1 + ypHE4Ev3_PI3K1, 2)/std::pow(sigma_ypHE4Ev3_PI3K1, 3);
            break;
        case 454:
            dJydsigma[454] = 1.0/sigma_ypE1HE4_PI3K1 - 1.0*std::pow(-mypE1HE4_PI3K1 + ypE1HE4_PI3K1, 2)/std::pow(sigma_ypE1HE4_PI3K1, 3);
            break;
        case 455:
            dJydsigma[455] = 1.0/sigma_ypE3HE4_PI3K1 - 1.0*std::pow(-mypE3HE4_PI3K1 + ypE3HE4_PI3K1, 2)/std::pow(sigma_ypE3HE4_PI3K1, 3);
            break;
        case 456:
            dJydsigma[456] = 1.0/sigma_ypHE4E4_PI3K1 - 1.0*std::pow(-mypHE4E4_PI3K1 + ypHE4E4_PI3K1, 2)/std::pow(sigma_ypHE4E4_PI3K1, 3);
            break;
        case 457:
            dJydsigma[457] = 1.0/sigma_ypHE4HE4_PI3K1 - 1.0*std::pow(-mypHE4HE4_PI3K1 + ypHE4HE4_PI3K1, 2)/std::pow(sigma_ypHE4HE4_PI3K1, 3);
            break;
        case 458:
            dJydsigma[458] = 1.0/sigma_ypHGF_Met_Met_PI3K1 - 1.0*std::pow(-mypHGF_Met_Met_PI3K1 + ypHGF_Met_Met_PI3K1, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K1, 3);
            break;
        case 459:
            dJydsigma[459] = 1.0/sigma_ypHGF_Met_HGF_Met_PI3K1 - 1.0*std::pow(-mypHGF_Met_HGF_Met_PI3K1 + ypHGF_Met_HGF_Met_PI3K1, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1, 3);
            break;
        case 460:
            dJydsigma[460] = 1.0/sigma_ypPPrPPr_PI3K1 - 1.0*std::pow(-mypPPrPPr_PI3K1 + ypPPrPPr_PI3K1, 2)/std::pow(sigma_ypPPrPPr_PI3K1, 3);
            break;
        case 461:
            dJydsigma[461] = 1.0/sigma_ypPPrPr_PI3K1 - 1.0*std::pow(-mypPPrPr_PI3K1 + ypPPrPr_PI3K1, 2)/std::pow(sigma_ypPPrPr_PI3K1, 3);
            break;
        case 462:
            dJydsigma[462] = 1.0/sigma_ypFFrFFr_PI3K1 - 1.0*std::pow(-mypFFrFFr_PI3K1 + ypFFrFFr_PI3K1, 2)/std::pow(sigma_ypFFrFFr_PI3K1, 3);
            break;
        case 463:
            dJydsigma[463] = 1.0/sigma_ypFFrFr_PI3K1 - 1.0*std::pow(-mypFFrFr_PI3K1 + ypFFrFr_PI3K1, 2)/std::pow(sigma_ypFFrFr_PI3K1, 3);
            break;
        case 464:
            dJydsigma[464] = 1.0/sigma_ypIIrIr_IRS_PI3K1 - 1.0*std::pow(-mypIIrIr_IRS_PI3K1 + ypIIrIr_IRS_PI3K1, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K1, 3);
            break;
        case 465:
            dJydsigma[465] = 1.0/sigma_ypINS_Isr_Isr_IRS_PI3K1 - 1.0*std::pow(-mypINS_Isr_Isr_IRS_PI3K1 + ypINS_Isr_Isr_IRS_PI3K1, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1, 3);
            break;
        case 466:
            dJydsigma[466] = 1.0/sigma_ypIIrIrI_IRS_PI3K1 - 1.0*std::pow(-mypIIrIrI_IRS_PI3K1 + ypIIrIrI_IRS_PI3K1, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K1, 3);
            break;
        case 467:
            dJydsigma[467] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_PI3K1 - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K1 + ypINS_Isr_Isr_INS_IRS_PI3K1, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1, 3);
            break;
        case 468:
            dJydsigma[468] = 1.0/sigma_ypEE1E2_PI3K2 - 1.0*std::pow(-mypEE1E2_PI3K2 + ypEE1E2_PI3K2, 2)/std::pow(sigma_ypEE1E2_PI3K2, 3);
            break;
        case 469:
            dJydsigma[469] = 1.0/sigma_ypEE1Ev3_PI3K2 - 1.0*std::pow(-mypEE1Ev3_PI3K2 + ypEE1Ev3_PI3K2, 2)/std::pow(sigma_ypEE1Ev3_PI3K2, 3);
            break;
        case 470:
            dJydsigma[470] = 1.0/sigma_ypEE1E1_PI3K2 - 1.0*std::pow(-mypEE1E1_PI3K2 + ypEE1E1_PI3K2, 2)/std::pow(sigma_ypEE1E1_PI3K2, 3);
            break;
        case 471:
            dJydsigma[471] = 1.0/sigma_ypEE1EE1_PI3K2 - 1.0*std::pow(-mypEE1EE1_PI3K2 + ypEE1EE1_PI3K2, 2)/std::pow(sigma_ypEE1EE1_PI3K2, 3);
            break;
        case 472:
            dJydsigma[472] = 1.0/sigma_ypEE1E3_PI3K2 - 1.0*std::pow(-mypEE1E3_PI3K2 + ypEE1E3_PI3K2, 2)/std::pow(sigma_ypEE1E3_PI3K2, 3);
            break;
        case 473:
            dJydsigma[473] = 1.0/sigma_ypEE1HE3_PI3K2 - 1.0*std::pow(-mypEE1HE3_PI3K2 + ypEE1HE3_PI3K2, 2)/std::pow(sigma_ypEE1HE3_PI3K2, 3);
            break;
        case 474:
            dJydsigma[474] = 1.0/sigma_ypEE1E4_PI3K2 - 1.0*std::pow(-mypEE1E4_PI3K2 + ypEE1E4_PI3K2, 2)/std::pow(sigma_ypEE1E4_PI3K2, 3);
            break;
        case 475:
            dJydsigma[475] = 1.0/sigma_ypEE1HE4_PI3K2 - 1.0*std::pow(-mypEE1HE4_PI3K2 + ypEE1HE4_PI3K2, 2)/std::pow(sigma_ypEE1HE4_PI3K2, 3);
            break;
        case 476:
            dJydsigma[476] = 1.0/sigma_ypE2HE3_PI3K2 - 1.0*std::pow(-mypE2HE3_PI3K2 + ypE2HE3_PI3K2, 2)/std::pow(sigma_ypE2HE3_PI3K2, 3);
            break;
        case 477:
            dJydsigma[477] = 1.0/sigma_ypHE3Ev3_PI3K2 - 1.0*std::pow(-mypHE3Ev3_PI3K2 + ypHE3Ev3_PI3K2, 2)/std::pow(sigma_ypHE3Ev3_PI3K2, 3);
            break;
        case 478:
            dJydsigma[478] = 1.0/sigma_ypE1HE3_PI3K2 - 1.0*std::pow(-mypE1HE3_PI3K2 + ypE1HE3_PI3K2, 2)/std::pow(sigma_ypE1HE3_PI3K2, 3);
            break;
        case 479:
            dJydsigma[479] = 1.0/sigma_ypHE3E4_PI3K2 - 1.0*std::pow(-mypHE3E4_PI3K2 + ypHE3E4_PI3K2, 2)/std::pow(sigma_ypHE3E4_PI3K2, 3);
            break;
        case 480:
            dJydsigma[480] = 1.0/sigma_ypHE3HE4_PI3K2 - 1.0*std::pow(-mypHE3HE4_PI3K2 + ypHE3HE4_PI3K2, 2)/std::pow(sigma_ypHE3HE4_PI3K2, 3);
            break;
        case 481:
            dJydsigma[481] = 1.0/sigma_ypE2HE4_PI3K2 - 1.0*std::pow(-mypE2HE4_PI3K2 + ypE2HE4_PI3K2, 2)/std::pow(sigma_ypE2HE4_PI3K2, 3);
            break;
        case 482:
            dJydsigma[482] = 1.0/sigma_ypHE4Ev3_PI3K2 - 1.0*std::pow(-mypHE4Ev3_PI3K2 + ypHE4Ev3_PI3K2, 2)/std::pow(sigma_ypHE4Ev3_PI3K2, 3);
            break;
        case 483:
            dJydsigma[483] = 1.0/sigma_ypE1HE4_PI3K2 - 1.0*std::pow(-mypE1HE4_PI3K2 + ypE1HE4_PI3K2, 2)/std::pow(sigma_ypE1HE4_PI3K2, 3);
            break;
        case 484:
            dJydsigma[484] = 1.0/sigma_ypE3HE4_PI3K2 - 1.0*std::pow(-mypE3HE4_PI3K2 + ypE3HE4_PI3K2, 2)/std::pow(sigma_ypE3HE4_PI3K2, 3);
            break;
        case 485:
            dJydsigma[485] = 1.0/sigma_ypHE4E4_PI3K2 - 1.0*std::pow(-mypHE4E4_PI3K2 + ypHE4E4_PI3K2, 2)/std::pow(sigma_ypHE4E4_PI3K2, 3);
            break;
        case 486:
            dJydsigma[486] = 1.0/sigma_ypHE4HE4_PI3K2 - 1.0*std::pow(-mypHE4HE4_PI3K2 + ypHE4HE4_PI3K2, 2)/std::pow(sigma_ypHE4HE4_PI3K2, 3);
            break;
        case 487:
            dJydsigma[487] = 1.0/sigma_ypHGF_Met_Met_PI3K2 - 1.0*std::pow(-mypHGF_Met_Met_PI3K2 + ypHGF_Met_Met_PI3K2, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K2, 3);
            break;
        case 488:
            dJydsigma[488] = 1.0/sigma_ypHGF_Met_HGF_Met_PI3K2 - 1.0*std::pow(-mypHGF_Met_HGF_Met_PI3K2 + ypHGF_Met_HGF_Met_PI3K2, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2, 3);
            break;
        case 489:
            dJydsigma[489] = 1.0/sigma_ypPPrPPr_PI3K2 - 1.0*std::pow(-mypPPrPPr_PI3K2 + ypPPrPPr_PI3K2, 2)/std::pow(sigma_ypPPrPPr_PI3K2, 3);
            break;
        case 490:
            dJydsigma[490] = 1.0/sigma_ypPPrPr_PI3K2 - 1.0*std::pow(-mypPPrPr_PI3K2 + ypPPrPr_PI3K2, 2)/std::pow(sigma_ypPPrPr_PI3K2, 3);
            break;
        case 491:
            dJydsigma[491] = 1.0/sigma_ypFFrFFr_PI3K2 - 1.0*std::pow(-mypFFrFFr_PI3K2 + ypFFrFFr_PI3K2, 2)/std::pow(sigma_ypFFrFFr_PI3K2, 3);
            break;
        case 492:
            dJydsigma[492] = 1.0/sigma_ypFFrFr_PI3K2 - 1.0*std::pow(-mypFFrFr_PI3K2 + ypFFrFr_PI3K2, 2)/std::pow(sigma_ypFFrFr_PI3K2, 3);
            break;
        case 493:
            dJydsigma[493] = 1.0/sigma_ypIIrIr_IRS_PI3K2 - 1.0*std::pow(-mypIIrIr_IRS_PI3K2 + ypIIrIr_IRS_PI3K2, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K2, 3);
            break;
        case 494:
            dJydsigma[494] = 1.0/sigma_ypINS_Isr_Isr_IRS_PI3K2 - 1.0*std::pow(-mypINS_Isr_Isr_IRS_PI3K2 + ypINS_Isr_Isr_IRS_PI3K2, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2, 3);
            break;
        case 495:
            dJydsigma[495] = 1.0/sigma_ypIIrIrI_IRS_PI3K2 - 1.0*std::pow(-mypIIrIrI_IRS_PI3K2 + ypIIrIrI_IRS_PI3K2, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K2, 3);
            break;
        case 496:
            dJydsigma[496] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_PI3K2 - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K2 + ypINS_Isr_Isr_INS_IRS_PI3K2, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2, 3);
            break;
        case 497:
            dJydsigma[497] = 1.0/sigma_ypEE1E2int_G2_SOS_RasD - 1.0*std::pow(-mypEE1E2int_G2_SOS_RasD + ypEE1E2int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E2int_G2_SOS_RasD, 3);
            break;
        case 498:
            dJydsigma[498] = 1.0/sigma_ypEE1Ev3int_G2_SOS_RasD - 1.0*std::pow(-mypEE1Ev3int_G2_SOS_RasD + ypEE1Ev3int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1Ev3int_G2_SOS_RasD, 3);
            break;
        case 499:
            dJydsigma[499] = 1.0/sigma_ypEE1E1int_G2_SOS_RasD - 1.0*std::pow(-mypEE1E1int_G2_SOS_RasD + ypEE1E1int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E1int_G2_SOS_RasD, 3);
            break;
        case 500:
            dJydsigma[500] = 1.0/sigma_ypEE1EE1int_G2_SOS_RasD - 1.0*std::pow(-mypEE1EE1int_G2_SOS_RasD + ypEE1EE1int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1EE1int_G2_SOS_RasD, 3);
            break;
        case 501:
            dJydsigma[501] = 1.0/sigma_ypEE1E3int_G2_SOS_RasD - 1.0*std::pow(-mypEE1E3int_G2_SOS_RasD + ypEE1E3int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E3int_G2_SOS_RasD, 3);
            break;
        case 502:
            dJydsigma[502] = 1.0/sigma_ypEE1HE3int_G2_SOS_RasD - 1.0*std::pow(-mypEE1HE3int_G2_SOS_RasD + ypEE1HE3int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE3int_G2_SOS_RasD, 3);
            break;
        case 503:
            dJydsigma[503] = 1.0/sigma_ypEE1E4int_G2_SOS_RasD - 1.0*std::pow(-mypEE1E4int_G2_SOS_RasD + ypEE1E4int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E4int_G2_SOS_RasD, 3);
            break;
        case 504:
            dJydsigma[504] = 1.0/sigma_ypEE1HE4int_G2_SOS_RasD - 1.0*std::pow(-mypEE1HE4int_G2_SOS_RasD + ypEE1HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE4int_G2_SOS_RasD, 3);
            break;
        case 505:
            dJydsigma[505] = 1.0/sigma_ypE2HE3int_G2_SOS_RasD - 1.0*std::pow(-mypE2HE3int_G2_SOS_RasD + ypE2HE3int_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE3int_G2_SOS_RasD, 3);
            break;
        case 506:
            dJydsigma[506] = 1.0/sigma_ypHE3Ev3int_G2_SOS_RasD - 1.0*std::pow(-mypHE3Ev3int_G2_SOS_RasD + ypHE3Ev3int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3Ev3int_G2_SOS_RasD, 3);
            break;
        case 507:
            dJydsigma[507] = 1.0/sigma_ypE1HE3int_G2_SOS_RasD - 1.0*std::pow(-mypE1HE3int_G2_SOS_RasD + ypE1HE3int_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE3int_G2_SOS_RasD, 3);
            break;
        case 508:
            dJydsigma[508] = 1.0/sigma_ypHE3E4int_G2_SOS_RasD - 1.0*std::pow(-mypHE3E4int_G2_SOS_RasD + ypHE3E4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3E4int_G2_SOS_RasD, 3);
            break;
        case 509:
            dJydsigma[509] = 1.0/sigma_ypHE3HE4int_G2_SOS_RasD - 1.0*std::pow(-mypHE3HE4int_G2_SOS_RasD + ypHE3HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3HE4int_G2_SOS_RasD, 3);
            break;
        case 510:
            dJydsigma[510] = 1.0/sigma_ypE2HE4int_G2_SOS_RasD - 1.0*std::pow(-mypE2HE4int_G2_SOS_RasD + ypE2HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE4int_G2_SOS_RasD, 3);
            break;
        case 511:
            dJydsigma[511] = 1.0/sigma_ypHE4Ev3int_G2_SOS_RasD - 1.0*std::pow(-mypHE4Ev3int_G2_SOS_RasD + ypHE4Ev3int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4Ev3int_G2_SOS_RasD, 3);
            break;
        case 512:
            dJydsigma[512] = 1.0/sigma_ypE1HE4int_G2_SOS_RasD - 1.0*std::pow(-mypE1HE4int_G2_SOS_RasD + ypE1HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE4int_G2_SOS_RasD, 3);
            break;
        case 513:
            dJydsigma[513] = 1.0/sigma_ypE3HE4int_G2_SOS_RasD - 1.0*std::pow(-mypE3HE4int_G2_SOS_RasD + ypE3HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypE3HE4int_G2_SOS_RasD, 3);
            break;
        case 514:
            dJydsigma[514] = 1.0/sigma_ypHE4E4int_G2_SOS_RasD - 1.0*std::pow(-mypHE4E4int_G2_SOS_RasD + ypHE4E4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4E4int_G2_SOS_RasD, 3);
            break;
        case 515:
            dJydsigma[515] = 1.0/sigma_ypHE4HE4int_G2_SOS_RasD - 1.0*std::pow(-mypHE4HE4int_G2_SOS_RasD + ypHE4HE4int_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4HE4int_G2_SOS_RasD, 3);
            break;
        case 516:
            dJydsigma[516] = 1.0/sigma_ypHGF_Met_Metint_G2_SOS_RasD - 1.0*std::pow(-mypHGF_Met_Metint_G2_SOS_RasD + ypHGF_Met_Metint_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_Metint_G2_SOS_RasD, 3);
            break;
        case 517:
            dJydsigma[517] = 1.0/sigma_ypHGF_Met_HGF_Metint_G2_SOS_RasD - 1.0*std::pow(-mypHGF_Met_HGF_Metint_G2_SOS_RasD + ypHGF_Met_HGF_Metint_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_HGF_Metint_G2_SOS_RasD, 3);
            break;
        case 518:
            dJydsigma[518] = 1.0/sigma_ypPPrPPrint_G2_SOS_RasD - 1.0*std::pow(-mypPPrPPrint_G2_SOS_RasD + ypPPrPPrint_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPPrint_G2_SOS_RasD, 3);
            break;
        case 519:
            dJydsigma[519] = 1.0/sigma_ypPPrPrint_G2_SOS_RasD - 1.0*std::pow(-mypPPrPrint_G2_SOS_RasD + ypPPrPrint_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPrint_G2_SOS_RasD, 3);
            break;
        case 520:
            dJydsigma[520] = 1.0/sigma_ypFFrFFrint_G2_SOS_RasD - 1.0*std::pow(-mypFFrFFrint_G2_SOS_RasD + ypFFrFFrint_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFFrint_G2_SOS_RasD, 3);
            break;
        case 521:
            dJydsigma[521] = 1.0/sigma_ypFFrFrint_G2_SOS_RasD - 1.0*std::pow(-mypFFrFrint_G2_SOS_RasD + ypFFrFrint_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFrint_G2_SOS_RasD, 3);
            break;
        case 522:
            dJydsigma[522] = 1.0/sigma_ypIIrIr_int_IRS_G2_SOS_RasD - 1.0*std::pow(-mypIIrIr_int_IRS_G2_SOS_RasD + ypIIrIr_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIr_int_IRS_G2_SOS_RasD, 3);
            break;
        case 523:
            dJydsigma[523] = 1.0/sigma_ypINS_Isr_Isr_int_IRS_G2_SOS_RasD - 1.0*std::pow(-mypINS_Isr_Isr_int_IRS_G2_SOS_RasD + ypINS_Isr_Isr_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_int_IRS_G2_SOS_RasD, 3);
            break;
        case 524:
            dJydsigma[524] = 1.0/sigma_ypIIrIrI_int_IRS_G2_SOS_RasD - 1.0*std::pow(-mypIIrIrI_int_IRS_G2_SOS_RasD + ypIIrIrI_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIrI_int_IRS_G2_SOS_RasD, 3);
            break;
        case 525:
            dJydsigma[525] = 1.0/sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD - 1.0*std::pow(-mypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD + ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD, 3);
            break;
        case 526:
            dJydsigma[526] = 1.0/sigma_ypEE1E2_G2_SOS_RasD - 1.0*std::pow(-mypEE1E2_G2_SOS_RasD + ypEE1E2_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E2_G2_SOS_RasD, 3);
            break;
        case 527:
            dJydsigma[527] = 1.0/sigma_ypEE1Ev3_G2_SOS_RasD - 1.0*std::pow(-mypEE1Ev3_G2_SOS_RasD + ypEE1Ev3_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1Ev3_G2_SOS_RasD, 3);
            break;
        case 528:
            dJydsigma[528] = 1.0/sigma_ypEE1E1_G2_SOS_RasD - 1.0*std::pow(-mypEE1E1_G2_SOS_RasD + ypEE1E1_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E1_G2_SOS_RasD, 3);
            break;
        case 529:
            dJydsigma[529] = 1.0/sigma_ypEE1EE1_G2_SOS_RasD - 1.0*std::pow(-mypEE1EE1_G2_SOS_RasD + ypEE1EE1_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1EE1_G2_SOS_RasD, 3);
            break;
        case 530:
            dJydsigma[530] = 1.0/sigma_ypEE1E3_G2_SOS_RasD - 1.0*std::pow(-mypEE1E3_G2_SOS_RasD + ypEE1E3_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E3_G2_SOS_RasD, 3);
            break;
        case 531:
            dJydsigma[531] = 1.0/sigma_ypEE1HE3_G2_SOS_RasD - 1.0*std::pow(-mypEE1HE3_G2_SOS_RasD + ypEE1HE3_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE3_G2_SOS_RasD, 3);
            break;
        case 532:
            dJydsigma[532] = 1.0/sigma_ypEE1E4_G2_SOS_RasD - 1.0*std::pow(-mypEE1E4_G2_SOS_RasD + ypEE1E4_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1E4_G2_SOS_RasD, 3);
            break;
        case 533:
            dJydsigma[533] = 1.0/sigma_ypEE1HE4_G2_SOS_RasD - 1.0*std::pow(-mypEE1HE4_G2_SOS_RasD + ypEE1HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypEE1HE4_G2_SOS_RasD, 3);
            break;
        case 534:
            dJydsigma[534] = 1.0/sigma_ypE2HE3_G2_SOS_RasD - 1.0*std::pow(-mypE2HE3_G2_SOS_RasD + ypE2HE3_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE3_G2_SOS_RasD, 3);
            break;
        case 535:
            dJydsigma[535] = 1.0/sigma_ypHE3Ev3_G2_SOS_RasD - 1.0*std::pow(-mypHE3Ev3_G2_SOS_RasD + ypHE3Ev3_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3Ev3_G2_SOS_RasD, 3);
            break;
        case 536:
            dJydsigma[536] = 1.0/sigma_ypE1HE3_G2_SOS_RasD - 1.0*std::pow(-mypE1HE3_G2_SOS_RasD + ypE1HE3_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE3_G2_SOS_RasD, 3);
            break;
        case 537:
            dJydsigma[537] = 1.0/sigma_ypHE3E4_G2_SOS_RasD - 1.0*std::pow(-mypHE3E4_G2_SOS_RasD + ypHE3E4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3E4_G2_SOS_RasD, 3);
            break;
        case 538:
            dJydsigma[538] = 1.0/sigma_ypHE3HE4_G2_SOS_RasD - 1.0*std::pow(-mypHE3HE4_G2_SOS_RasD + ypHE3HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE3HE4_G2_SOS_RasD, 3);
            break;
        case 539:
            dJydsigma[539] = 1.0/sigma_ypE2HE4_G2_SOS_RasD - 1.0*std::pow(-mypE2HE4_G2_SOS_RasD + ypE2HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypE2HE4_G2_SOS_RasD, 3);
            break;
        case 540:
            dJydsigma[540] = 1.0/sigma_ypHE4Ev3_G2_SOS_RasD - 1.0*std::pow(-mypHE4Ev3_G2_SOS_RasD + ypHE4Ev3_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4Ev3_G2_SOS_RasD, 3);
            break;
        case 541:
            dJydsigma[541] = 1.0/sigma_ypE1HE4_G2_SOS_RasD - 1.0*std::pow(-mypE1HE4_G2_SOS_RasD + ypE1HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypE1HE4_G2_SOS_RasD, 3);
            break;
        case 542:
            dJydsigma[542] = 1.0/sigma_ypE3HE4_G2_SOS_RasD - 1.0*std::pow(-mypE3HE4_G2_SOS_RasD + ypE3HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypE3HE4_G2_SOS_RasD, 3);
            break;
        case 543:
            dJydsigma[543] = 1.0/sigma_ypHE4E4_G2_SOS_RasD - 1.0*std::pow(-mypHE4E4_G2_SOS_RasD + ypHE4E4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4E4_G2_SOS_RasD, 3);
            break;
        case 544:
            dJydsigma[544] = 1.0/sigma_ypHE4HE4_G2_SOS_RasD - 1.0*std::pow(-mypHE4HE4_G2_SOS_RasD + ypHE4HE4_G2_SOS_RasD, 2)/std::pow(sigma_ypHE4HE4_G2_SOS_RasD, 3);
            break;
        case 545:
            dJydsigma[545] = 1.0/sigma_ypHGF_Met_Met_G2_SOS_RasD - 1.0*std::pow(-mypHGF_Met_Met_G2_SOS_RasD + ypHGF_Met_Met_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_Met_G2_SOS_RasD, 3);
            break;
        case 546:
            dJydsigma[546] = 1.0/sigma_ypHGF_Met_HGF_Met_G2_SOS_RasD - 1.0*std::pow(-mypHGF_Met_HGF_Met_G2_SOS_RasD + ypHGF_Met_HGF_Met_G2_SOS_RasD, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_G2_SOS_RasD, 3);
            break;
        case 547:
            dJydsigma[547] = 1.0/sigma_ypPPrPPr_G2_SOS_RasD - 1.0*std::pow(-mypPPrPPr_G2_SOS_RasD + ypPPrPPr_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPPr_G2_SOS_RasD, 3);
            break;
        case 548:
            dJydsigma[548] = 1.0/sigma_ypPPrPr_G2_SOS_RasD - 1.0*std::pow(-mypPPrPr_G2_SOS_RasD + ypPPrPr_G2_SOS_RasD, 2)/std::pow(sigma_ypPPrPr_G2_SOS_RasD, 3);
            break;
        case 549:
            dJydsigma[549] = 1.0/sigma_ypFFrFFr_G2_SOS_RasD - 1.0*std::pow(-mypFFrFFr_G2_SOS_RasD + ypFFrFFr_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFFr_G2_SOS_RasD, 3);
            break;
        case 550:
            dJydsigma[550] = 1.0/sigma_ypFFrFr_G2_SOS_RasD - 1.0*std::pow(-mypFFrFr_G2_SOS_RasD + ypFFrFr_G2_SOS_RasD, 2)/std::pow(sigma_ypFFrFr_G2_SOS_RasD, 3);
            break;
        case 551:
            dJydsigma[551] = 1.0/sigma_ypIIrIr_IRS_G2_SOS_RasD - 1.0*std::pow(-mypIIrIr_IRS_G2_SOS_RasD + ypIIrIr_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIr_IRS_G2_SOS_RasD, 3);
            break;
        case 552:
            dJydsigma[552] = 1.0/sigma_ypINS_Isr_Isr_IRS_G2_SOS_RasD - 1.0*std::pow(-mypINS_Isr_Isr_IRS_G2_SOS_RasD + ypINS_Isr_Isr_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_G2_SOS_RasD, 3);
            break;
        case 553:
            dJydsigma[553] = 1.0/sigma_ypIIrIrI_IRS_G2_SOS_RasD - 1.0*std::pow(-mypIIrIrI_IRS_G2_SOS_RasD + ypIIrIrI_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypIIrIrI_IRS_G2_SOS_RasD, 3);
            break;
        case 554:
            dJydsigma[554] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_G2_SOS_RasD + ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_G2_SOS_RasD, 3);
            break;
        case 555:
            dJydsigma[555] = 1.0/sigma_ypEE1E2_PLCg_PIP2 - 1.0*std::pow(-mypEE1E2_PLCg_PIP2 + ypEE1E2_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E2_PLCg_PIP2, 3);
            break;
        case 556:
            dJydsigma[556] = 1.0/sigma_ypEE1Ev3_PLCg_PIP2 - 1.0*std::pow(-mypEE1Ev3_PLCg_PIP2 + ypEE1Ev3_PLCg_PIP2, 2)/std::pow(sigma_ypEE1Ev3_PLCg_PIP2, 3);
            break;
        case 557:
            dJydsigma[557] = 1.0/sigma_ypEE1E1_PLCg_PIP2 - 1.0*std::pow(-mypEE1E1_PLCg_PIP2 + ypEE1E1_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E1_PLCg_PIP2, 3);
            break;
        case 558:
            dJydsigma[558] = 1.0/sigma_ypEE1EE1_PLCg_PIP2 - 1.0*std::pow(-mypEE1EE1_PLCg_PIP2 + ypEE1EE1_PLCg_PIP2, 2)/std::pow(sigma_ypEE1EE1_PLCg_PIP2, 3);
            break;
        case 559:
            dJydsigma[559] = 1.0/sigma_ypEE1E3_PLCg_PIP2 - 1.0*std::pow(-mypEE1E3_PLCg_PIP2 + ypEE1E3_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E3_PLCg_PIP2, 3);
            break;
        case 560:
            dJydsigma[560] = 1.0/sigma_ypEE1HE3_PLCg_PIP2 - 1.0*std::pow(-mypEE1HE3_PLCg_PIP2 + ypEE1HE3_PLCg_PIP2, 2)/std::pow(sigma_ypEE1HE3_PLCg_PIP2, 3);
            break;
        case 561:
            dJydsigma[561] = 1.0/sigma_ypEE1E4_PLCg_PIP2 - 1.0*std::pow(-mypEE1E4_PLCg_PIP2 + ypEE1E4_PLCg_PIP2, 2)/std::pow(sigma_ypEE1E4_PLCg_PIP2, 3);
            break;
        case 562:
            dJydsigma[562] = 1.0/sigma_ypEE1HE4_PLCg_PIP2 - 1.0*std::pow(-mypEE1HE4_PLCg_PIP2 + ypEE1HE4_PLCg_PIP2, 2)/std::pow(sigma_ypEE1HE4_PLCg_PIP2, 3);
            break;
        case 563:
            dJydsigma[563] = 1.0/sigma_ypE2HE3_PLCg_PIP2 - 1.0*std::pow(-mypE2HE3_PLCg_PIP2 + ypE2HE3_PLCg_PIP2, 2)/std::pow(sigma_ypE2HE3_PLCg_PIP2, 3);
            break;
        case 564:
            dJydsigma[564] = 1.0/sigma_ypHE3Ev3_PLCg_PIP2 - 1.0*std::pow(-mypHE3Ev3_PLCg_PIP2 + ypHE3Ev3_PLCg_PIP2, 2)/std::pow(sigma_ypHE3Ev3_PLCg_PIP2, 3);
            break;
        case 565:
            dJydsigma[565] = 1.0/sigma_ypE1HE3_PLCg_PIP2 - 1.0*std::pow(-mypE1HE3_PLCg_PIP2 + ypE1HE3_PLCg_PIP2, 2)/std::pow(sigma_ypE1HE3_PLCg_PIP2, 3);
            break;
        case 566:
            dJydsigma[566] = 1.0/sigma_ypHE3E4_PLCg_PIP2 - 1.0*std::pow(-mypHE3E4_PLCg_PIP2 + ypHE3E4_PLCg_PIP2, 2)/std::pow(sigma_ypHE3E4_PLCg_PIP2, 3);
            break;
        case 567:
            dJydsigma[567] = 1.0/sigma_ypHE3HE4_PLCg_PIP2 - 1.0*std::pow(-mypHE3HE4_PLCg_PIP2 + ypHE3HE4_PLCg_PIP2, 2)/std::pow(sigma_ypHE3HE4_PLCg_PIP2, 3);
            break;
        case 568:
            dJydsigma[568] = 1.0/sigma_ypE2HE4_PLCg_PIP2 - 1.0*std::pow(-mypE2HE4_PLCg_PIP2 + ypE2HE4_PLCg_PIP2, 2)/std::pow(sigma_ypE2HE4_PLCg_PIP2, 3);
            break;
        case 569:
            dJydsigma[569] = 1.0/sigma_ypHE4Ev3_PLCg_PIP2 - 1.0*std::pow(-mypHE4Ev3_PLCg_PIP2 + ypHE4Ev3_PLCg_PIP2, 2)/std::pow(sigma_ypHE4Ev3_PLCg_PIP2, 3);
            break;
        case 570:
            dJydsigma[570] = 1.0/sigma_ypE1HE4_PLCg_PIP2 - 1.0*std::pow(-mypE1HE4_PLCg_PIP2 + ypE1HE4_PLCg_PIP2, 2)/std::pow(sigma_ypE1HE4_PLCg_PIP2, 3);
            break;
        case 571:
            dJydsigma[571] = 1.0/sigma_ypE3HE4_PLCg_PIP2 - 1.0*std::pow(-mypE3HE4_PLCg_PIP2 + ypE3HE4_PLCg_PIP2, 2)/std::pow(sigma_ypE3HE4_PLCg_PIP2, 3);
            break;
        case 572:
            dJydsigma[572] = 1.0/sigma_ypHE4E4_PLCg_PIP2 - 1.0*std::pow(-mypHE4E4_PLCg_PIP2 + ypHE4E4_PLCg_PIP2, 2)/std::pow(sigma_ypHE4E4_PLCg_PIP2, 3);
            break;
        case 573:
            dJydsigma[573] = 1.0/sigma_ypHE4HE4_PLCg_PIP2 - 1.0*std::pow(-mypHE4HE4_PLCg_PIP2 + ypHE4HE4_PLCg_PIP2, 2)/std::pow(sigma_ypHE4HE4_PLCg_PIP2, 3);
            break;
        case 574:
            dJydsigma[574] = 1.0/sigma_ypHGF_Met_Met_PLCg_PIP2 - 1.0*std::pow(-mypHGF_Met_Met_PLCg_PIP2 + ypHGF_Met_Met_PLCg_PIP2, 2)/std::pow(sigma_ypHGF_Met_Met_PLCg_PIP2, 3);
            break;
        case 575:
            dJydsigma[575] = 1.0/sigma_ypHGF_Met_HGF_Met_PLCg_PIP2 - 1.0*std::pow(-mypHGF_Met_HGF_Met_PLCg_PIP2 + ypHGF_Met_HGF_Met_PLCg_PIP2, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PLCg_PIP2, 3);
            break;
        case 576:
            dJydsigma[576] = 1.0/sigma_ypPPrPPr_PLCg_PIP2 - 1.0*std::pow(-mypPPrPPr_PLCg_PIP2 + ypPPrPPr_PLCg_PIP2, 2)/std::pow(sigma_ypPPrPPr_PLCg_PIP2, 3);
            break;
        case 577:
            dJydsigma[577] = 1.0/sigma_ypPPrPr_PLCg_PIP2 - 1.0*std::pow(-mypPPrPr_PLCg_PIP2 + ypPPrPr_PLCg_PIP2, 2)/std::pow(sigma_ypPPrPr_PLCg_PIP2, 3);
            break;
        case 578:
            dJydsigma[578] = 1.0/sigma_ypFFrFFr_PLCg_PIP2 - 1.0*std::pow(-mypFFrFFr_PLCg_PIP2 + ypFFrFFr_PLCg_PIP2, 2)/std::pow(sigma_ypFFrFFr_PLCg_PIP2, 3);
            break;
        case 579:
            dJydsigma[579] = 1.0/sigma_ypFFrFr_PLCg_PIP2 - 1.0*std::pow(-mypFFrFr_PLCg_PIP2 + ypFFrFr_PLCg_PIP2, 2)/std::pow(sigma_ypFFrFr_PLCg_PIP2, 3);
            break;
        case 580:
            dJydsigma[580] = 1.0/sigma_ypIIrIr_IRS_PLCg_PIP2 - 1.0*std::pow(-mypIIrIr_IRS_PLCg_PIP2 + ypIIrIr_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypIIrIr_IRS_PLCg_PIP2, 3);
            break;
        case 581:
            dJydsigma[581] = 1.0/sigma_ypINS_Isr_Isr_IRS_PLCg_PIP2 - 1.0*std::pow(-mypINS_Isr_Isr_IRS_PLCg_PIP2 + ypINS_Isr_Isr_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PLCg_PIP2, 3);
            break;
        case 582:
            dJydsigma[582] = 1.0/sigma_ypIIrIrI_IRS_PLCg_PIP2 - 1.0*std::pow(-mypIIrIrI_IRS_PLCg_PIP2 + ypIIrIrI_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypIIrIrI_IRS_PLCg_PIP2, 3);
            break;
        case 583:
            dJydsigma[583] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_PLCg_PIP2 - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_PLCg_PIP2 + ypINS_Isr_Isr_INS_IRS_PLCg_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PLCg_PIP2, 3);
            break;
        case 584:
            dJydsigma[584] = 1.0/sigma_ypEE1E2_PI3K1_PIP2 - 1.0*std::pow(-mypEE1E2_PI3K1_PIP2 + ypEE1E2_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E2_PI3K1_PIP2, 3);
            break;
        case 585:
            dJydsigma[585] = 1.0/sigma_ypEE1Ev3_PI3K1_PIP2 - 1.0*std::pow(-mypEE1Ev3_PI3K1_PIP2 + ypEE1Ev3_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1Ev3_PI3K1_PIP2, 3);
            break;
        case 586:
            dJydsigma[586] = 1.0/sigma_ypEE1E1_PI3K1_PIP2 - 1.0*std::pow(-mypEE1E1_PI3K1_PIP2 + ypEE1E1_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E1_PI3K1_PIP2, 3);
            break;
        case 587:
            dJydsigma[587] = 1.0/sigma_ypEE1EE1_PI3K1_PIP2 - 1.0*std::pow(-mypEE1EE1_PI3K1_PIP2 + ypEE1EE1_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1EE1_PI3K1_PIP2, 3);
            break;
        case 588:
            dJydsigma[588] = 1.0/sigma_ypEE1E3_PI3K1_PIP2 - 1.0*std::pow(-mypEE1E3_PI3K1_PIP2 + ypEE1E3_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E3_PI3K1_PIP2, 3);
            break;
        case 589:
            dJydsigma[589] = 1.0/sigma_ypEE1HE3_PI3K1_PIP2 - 1.0*std::pow(-mypEE1HE3_PI3K1_PIP2 + ypEE1HE3_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1HE3_PI3K1_PIP2, 3);
            break;
        case 590:
            dJydsigma[590] = 1.0/sigma_ypEE1E4_PI3K1_PIP2 - 1.0*std::pow(-mypEE1E4_PI3K1_PIP2 + ypEE1E4_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1E4_PI3K1_PIP2, 3);
            break;
        case 591:
            dJydsigma[591] = 1.0/sigma_ypEE1HE4_PI3K1_PIP2 - 1.0*std::pow(-mypEE1HE4_PI3K1_PIP2 + ypEE1HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypEE1HE4_PI3K1_PIP2, 3);
            break;
        case 592:
            dJydsigma[592] = 1.0/sigma_ypE2HE3_PI3K1_PIP2 - 1.0*std::pow(-mypE2HE3_PI3K1_PIP2 + ypE2HE3_PI3K1_PIP2, 2)/std::pow(sigma_ypE2HE3_PI3K1_PIP2, 3);
            break;
        case 593:
            dJydsigma[593] = 1.0/sigma_ypHE3Ev3_PI3K1_PIP2 - 1.0*std::pow(-mypHE3Ev3_PI3K1_PIP2 + ypHE3Ev3_PI3K1_PIP2, 2)/std::pow(sigma_ypHE3Ev3_PI3K1_PIP2, 3);
            break;
        case 594:
            dJydsigma[594] = 1.0/sigma_ypE1HE3_PI3K1_PIP2 - 1.0*std::pow(-mypE1HE3_PI3K1_PIP2 + ypE1HE3_PI3K1_PIP2, 2)/std::pow(sigma_ypE1HE3_PI3K1_PIP2, 3);
            break;
        case 595:
            dJydsigma[595] = 1.0/sigma_ypHE3E4_PI3K1_PIP2 - 1.0*std::pow(-mypHE3E4_PI3K1_PIP2 + ypHE3E4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE3E4_PI3K1_PIP2, 3);
            break;
        case 596:
            dJydsigma[596] = 1.0/sigma_ypHE3HE4_PI3K1_PIP2 - 1.0*std::pow(-mypHE3HE4_PI3K1_PIP2 + ypHE3HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE3HE4_PI3K1_PIP2, 3);
            break;
        case 597:
            dJydsigma[597] = 1.0/sigma_ypE2HE4_PI3K1_PIP2 - 1.0*std::pow(-mypE2HE4_PI3K1_PIP2 + ypE2HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypE2HE4_PI3K1_PIP2, 3);
            break;
        case 598:
            dJydsigma[598] = 1.0/sigma_ypHE4Ev3_PI3K1_PIP2 - 1.0*std::pow(-mypHE4Ev3_PI3K1_PIP2 + ypHE4Ev3_PI3K1_PIP2, 2)/std::pow(sigma_ypHE4Ev3_PI3K1_PIP2, 3);
            break;
        case 599:
            dJydsigma[599] = 1.0/sigma_ypE1HE4_PI3K1_PIP2 - 1.0*std::pow(-mypE1HE4_PI3K1_PIP2 + ypE1HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypE1HE4_PI3K1_PIP2, 3);
            break;
        case 600:
            dJydsigma[600] = 1.0/sigma_ypE3HE4_PI3K1_PIP2 - 1.0*std::pow(-mypE3HE4_PI3K1_PIP2 + ypE3HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypE3HE4_PI3K1_PIP2, 3);
            break;
        case 601:
            dJydsigma[601] = 1.0/sigma_ypHE4E4_PI3K1_PIP2 - 1.0*std::pow(-mypHE4E4_PI3K1_PIP2 + ypHE4E4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE4E4_PI3K1_PIP2, 3);
            break;
        case 602:
            dJydsigma[602] = 1.0/sigma_ypHE4HE4_PI3K1_PIP2 - 1.0*std::pow(-mypHE4HE4_PI3K1_PIP2 + ypHE4HE4_PI3K1_PIP2, 2)/std::pow(sigma_ypHE4HE4_PI3K1_PIP2, 3);
            break;
        case 603:
            dJydsigma[603] = 1.0/sigma_ypHGF_Met_Met_PI3K1_PIP2 - 1.0*std::pow(-mypHGF_Met_Met_PI3K1_PIP2 + ypHGF_Met_Met_PI3K1_PIP2, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K1_PIP2, 3);
            break;
        case 604:
            dJydsigma[604] = 1.0/sigma_ypHGF_Met_HGF_Met_PI3K1_PIP2 - 1.0*std::pow(-mypHGF_Met_HGF_Met_PI3K1_PIP2 + ypHGF_Met_HGF_Met_PI3K1_PIP2, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K1_PIP2, 3);
            break;
        case 605:
            dJydsigma[605] = 1.0/sigma_ypPPrPPr_PI3K1_PIP2 - 1.0*std::pow(-mypPPrPPr_PI3K1_PIP2 + ypPPrPPr_PI3K1_PIP2, 2)/std::pow(sigma_ypPPrPPr_PI3K1_PIP2, 3);
            break;
        case 606:
            dJydsigma[606] = 1.0/sigma_ypPPrPr_PI3K1_PIP2 - 1.0*std::pow(-mypPPrPr_PI3K1_PIP2 + ypPPrPr_PI3K1_PIP2, 2)/std::pow(sigma_ypPPrPr_PI3K1_PIP2, 3);
            break;
        case 607:
            dJydsigma[607] = 1.0/sigma_ypFFrFFr_PI3K1_PIP2 - 1.0*std::pow(-mypFFrFFr_PI3K1_PIP2 + ypFFrFFr_PI3K1_PIP2, 2)/std::pow(sigma_ypFFrFFr_PI3K1_PIP2, 3);
            break;
        case 608:
            dJydsigma[608] = 1.0/sigma_ypFFrFr_PI3K1_PIP2 - 1.0*std::pow(-mypFFrFr_PI3K1_PIP2 + ypFFrFr_PI3K1_PIP2, 2)/std::pow(sigma_ypFFrFr_PI3K1_PIP2, 3);
            break;
        case 609:
            dJydsigma[609] = 1.0/sigma_ypIIrIr_IRS_PI3K1_PIP2 - 1.0*std::pow(-mypIIrIr_IRS_PI3K1_PIP2 + ypIIrIr_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K1_PIP2, 3);
            break;
        case 610:
            dJydsigma[610] = 1.0/sigma_ypINS_Isr_Isr_IRS_PI3K1_PIP2 - 1.0*std::pow(-mypINS_Isr_Isr_IRS_PI3K1_PIP2 + ypINS_Isr_Isr_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K1_PIP2, 3);
            break;
        case 611:
            dJydsigma[611] = 1.0/sigma_ypIIrIrI_IRS_PI3K1_PIP2 - 1.0*std::pow(-mypIIrIrI_IRS_PI3K1_PIP2 + ypIIrIrI_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K1_PIP2, 3);
            break;
        case 612:
            dJydsigma[612] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2 - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K1_PIP2 + ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K1_PIP2, 3);
            break;
        case 613:
            dJydsigma[613] = 1.0/sigma_ypEE1E2_PI3K2_PIP - 1.0*std::pow(-mypEE1E2_PI3K2_PIP + ypEE1E2_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E2_PI3K2_PIP, 3);
            break;
        case 614:
            dJydsigma[614] = 1.0/sigma_ypEE1Ev3_PI3K2_PIP - 1.0*std::pow(-mypEE1Ev3_PI3K2_PIP + ypEE1Ev3_PI3K2_PIP, 2)/std::pow(sigma_ypEE1Ev3_PI3K2_PIP, 3);
            break;
        case 615:
            dJydsigma[615] = 1.0/sigma_ypEE1E1_PI3K2_PIP - 1.0*std::pow(-mypEE1E1_PI3K2_PIP + ypEE1E1_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E1_PI3K2_PIP, 3);
            break;
        case 616:
            dJydsigma[616] = 1.0/sigma_ypEE1EE1_PI3K2_PIP - 1.0*std::pow(-mypEE1EE1_PI3K2_PIP + ypEE1EE1_PI3K2_PIP, 2)/std::pow(sigma_ypEE1EE1_PI3K2_PIP, 3);
            break;
        case 617:
            dJydsigma[617] = 1.0/sigma_ypEE1E3_PI3K2_PIP - 1.0*std::pow(-mypEE1E3_PI3K2_PIP + ypEE1E3_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E3_PI3K2_PIP, 3);
            break;
        case 618:
            dJydsigma[618] = 1.0/sigma_ypEE1HE3_PI3K2_PIP - 1.0*std::pow(-mypEE1HE3_PI3K2_PIP + ypEE1HE3_PI3K2_PIP, 2)/std::pow(sigma_ypEE1HE3_PI3K2_PIP, 3);
            break;
        case 619:
            dJydsigma[619] = 1.0/sigma_ypEE1E4_PI3K2_PIP - 1.0*std::pow(-mypEE1E4_PI3K2_PIP + ypEE1E4_PI3K2_PIP, 2)/std::pow(sigma_ypEE1E4_PI3K2_PIP, 3);
            break;
        case 620:
            dJydsigma[620] = 1.0/sigma_ypEE1HE4_PI3K2_PIP - 1.0*std::pow(-mypEE1HE4_PI3K2_PIP + ypEE1HE4_PI3K2_PIP, 2)/std::pow(sigma_ypEE1HE4_PI3K2_PIP, 3);
            break;
        case 621:
            dJydsigma[621] = 1.0/sigma_ypE2HE3_PI3K2_PIP - 1.0*std::pow(-mypE2HE3_PI3K2_PIP + ypE2HE3_PI3K2_PIP, 2)/std::pow(sigma_ypE2HE3_PI3K2_PIP, 3);
            break;
        case 622:
            dJydsigma[622] = 1.0/sigma_ypHE3Ev3_PI3K2_PIP - 1.0*std::pow(-mypHE3Ev3_PI3K2_PIP + ypHE3Ev3_PI3K2_PIP, 2)/std::pow(sigma_ypHE3Ev3_PI3K2_PIP, 3);
            break;
        case 623:
            dJydsigma[623] = 1.0/sigma_ypE1HE3_PI3K2_PIP - 1.0*std::pow(-mypE1HE3_PI3K2_PIP + ypE1HE3_PI3K2_PIP, 2)/std::pow(sigma_ypE1HE3_PI3K2_PIP, 3);
            break;
        case 624:
            dJydsigma[624] = 1.0/sigma_ypHE3E4_PI3K2_PIP - 1.0*std::pow(-mypHE3E4_PI3K2_PIP + ypHE3E4_PI3K2_PIP, 2)/std::pow(sigma_ypHE3E4_PI3K2_PIP, 3);
            break;
        case 625:
            dJydsigma[625] = 1.0/sigma_ypHE3HE4_PI3K2_PIP - 1.0*std::pow(-mypHE3HE4_PI3K2_PIP + ypHE3HE4_PI3K2_PIP, 2)/std::pow(sigma_ypHE3HE4_PI3K2_PIP, 3);
            break;
        case 626:
            dJydsigma[626] = 1.0/sigma_ypE2HE4_PI3K2_PIP - 1.0*std::pow(-mypE2HE4_PI3K2_PIP + ypE2HE4_PI3K2_PIP, 2)/std::pow(sigma_ypE2HE4_PI3K2_PIP, 3);
            break;
        case 627:
            dJydsigma[627] = 1.0/sigma_ypHE4Ev3_PI3K2_PIP - 1.0*std::pow(-mypHE4Ev3_PI3K2_PIP + ypHE4Ev3_PI3K2_PIP, 2)/std::pow(sigma_ypHE4Ev3_PI3K2_PIP, 3);
            break;
        case 628:
            dJydsigma[628] = 1.0/sigma_ypE1HE4_PI3K2_PIP - 1.0*std::pow(-mypE1HE4_PI3K2_PIP + ypE1HE4_PI3K2_PIP, 2)/std::pow(sigma_ypE1HE4_PI3K2_PIP, 3);
            break;
        case 629:
            dJydsigma[629] = 1.0/sigma_ypE3HE4_PI3K2_PIP - 1.0*std::pow(-mypE3HE4_PI3K2_PIP + ypE3HE4_PI3K2_PIP, 2)/std::pow(sigma_ypE3HE4_PI3K2_PIP, 3);
            break;
        case 630:
            dJydsigma[630] = 1.0/sigma_ypHE4E4_PI3K2_PIP - 1.0*std::pow(-mypHE4E4_PI3K2_PIP + ypHE4E4_PI3K2_PIP, 2)/std::pow(sigma_ypHE4E4_PI3K2_PIP, 3);
            break;
        case 631:
            dJydsigma[631] = 1.0/sigma_ypHE4HE4_PI3K2_PIP - 1.0*std::pow(-mypHE4HE4_PI3K2_PIP + ypHE4HE4_PI3K2_PIP, 2)/std::pow(sigma_ypHE4HE4_PI3K2_PIP, 3);
            break;
        case 632:
            dJydsigma[632] = 1.0/sigma_ypHGF_Met_Met_PI3K2_PIP - 1.0*std::pow(-mypHGF_Met_Met_PI3K2_PIP + ypHGF_Met_Met_PI3K2_PIP, 2)/std::pow(sigma_ypHGF_Met_Met_PI3K2_PIP, 3);
            break;
        case 633:
            dJydsigma[633] = 1.0/sigma_ypHGF_Met_HGF_Met_PI3K2_PIP - 1.0*std::pow(-mypHGF_Met_HGF_Met_PI3K2_PIP + ypHGF_Met_HGF_Met_PI3K2_PIP, 2)/std::pow(sigma_ypHGF_Met_HGF_Met_PI3K2_PIP, 3);
            break;
        case 634:
            dJydsigma[634] = 1.0/sigma_ypPPrPPr_PI3K2_PIP - 1.0*std::pow(-mypPPrPPr_PI3K2_PIP + ypPPrPPr_PI3K2_PIP, 2)/std::pow(sigma_ypPPrPPr_PI3K2_PIP, 3);
            break;
        case 635:
            dJydsigma[635] = 1.0/sigma_ypPPrPr_PI3K2_PIP - 1.0*std::pow(-mypPPrPr_PI3K2_PIP + ypPPrPr_PI3K2_PIP, 2)/std::pow(sigma_ypPPrPr_PI3K2_PIP, 3);
            break;
        case 636:
            dJydsigma[636] = 1.0/sigma_ypFFrFFr_PI3K2_PIP - 1.0*std::pow(-mypFFrFFr_PI3K2_PIP + ypFFrFFr_PI3K2_PIP, 2)/std::pow(sigma_ypFFrFFr_PI3K2_PIP, 3);
            break;
        case 637:
            dJydsigma[637] = 1.0/sigma_ypFFrFr_PI3K2_PIP - 1.0*std::pow(-mypFFrFr_PI3K2_PIP + ypFFrFr_PI3K2_PIP, 2)/std::pow(sigma_ypFFrFr_PI3K2_PIP, 3);
            break;
        case 638:
            dJydsigma[638] = 1.0/sigma_ypIIrIr_IRS_PI3K2_PIP - 1.0*std::pow(-mypIIrIr_IRS_PI3K2_PIP + ypIIrIr_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypIIrIr_IRS_PI3K2_PIP, 3);
            break;
        case 639:
            dJydsigma[639] = 1.0/sigma_ypINS_Isr_Isr_IRS_PI3K2_PIP - 1.0*std::pow(-mypINS_Isr_Isr_IRS_PI3K2_PIP + ypINS_Isr_Isr_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypINS_Isr_Isr_IRS_PI3K2_PIP, 3);
            break;
        case 640:
            dJydsigma[640] = 1.0/sigma_ypIIrIrI_IRS_PI3K2_PIP - 1.0*std::pow(-mypIIrIrI_IRS_PI3K2_PIP + ypIIrIrI_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypIIrIrI_IRS_PI3K2_PIP, 3);
            break;
        case 641:
            dJydsigma[641] = 1.0/sigma_ypINS_Isr_Isr_INS_IRS_PI3K2_PIP - 1.0*std::pow(-mypINS_Isr_Isr_INS_IRS_PI3K2_PIP + ypINS_Isr_Isr_INS_IRS_PI3K2_PIP, 2)/std::pow(sigma_ypINS_Isr_Isr_INS_IRS_PI3K2_PIP, 3);
            break;
        case 642:
            dJydsigma[642] = 1.0/sigma_yIRS - 1.0*std::pow(-myIRS + yIRS, 2)/std::pow(sigma_yIRS, 3);
            break;
        case 643:
            dJydsigma[643] = 1.0/sigma_ySp - 1.0*std::pow(-mySp + ySp, 2)/std::pow(sigma_ySp, 3);
            break;
        case 644:
            dJydsigma[644] = 1.0/sigma_yCbl - 1.0*std::pow(-myCbl + yCbl, 2)/std::pow(sigma_yCbl, 3);
            break;
        case 645:
            dJydsigma[645] = 1.0/sigma_yG2 - 1.0*std::pow(-myG2 + yG2, 2)/std::pow(sigma_yG2, 3);
            break;
        case 646:
            dJydsigma[646] = 1.0/sigma_yG2_SOS - 1.0*std::pow(-myG2_SOS + yG2_SOS, 2)/std::pow(sigma_yG2_SOS, 3);
            break;
        case 647:
            dJydsigma[647] = 1.0/sigma_yG2_pSOS - 1.0*std::pow(-myG2_pSOS + yG2_pSOS, 2)/std::pow(sigma_yG2_pSOS, 3);
            break;
        case 648:
            dJydsigma[648] = 1.0/sigma_yPLCg - 1.0*std::pow(-myPLCg + yPLCg, 2)/std::pow(sigma_yPLCg, 3);
            break;
        case 649:
            dJydsigma[649] = 1.0/sigma_yPI3KC1 - 1.0*std::pow(-myPI3KC1 + yPI3KC1, 2)/std::pow(sigma_yPI3KC1, 3);
            break;
        case 650:
            dJydsigma[650] = 1.0/sigma_yPI3KR1 - 1.0*std::pow(-myPI3KR1 + yPI3KR1, 2)/std::pow(sigma_yPI3KR1, 3);
            break;
        case 651:
            dJydsigma[651] = 1.0/sigma_yPI3K1 - 1.0*std::pow(-myPI3K1 + yPI3K1, 2)/std::pow(sigma_yPI3K1, 3);
            break;
        case 652:
            dJydsigma[652] = 1.0/sigma_ypPI3K1 - 1.0*std::pow(-mypPI3K1 + ypPI3K1, 2)/std::pow(sigma_ypPI3K1, 3);
            break;
        case 653:
            dJydsigma[653] = 1.0/sigma_yPI3K2 - 1.0*std::pow(-myPI3K2 + yPI3K2, 2)/std::pow(sigma_yPI3K2, 3);
            break;
        case 654:
            dJydsigma[654] = 1.0/sigma_ymTORC1 - 1.0*std::pow(-mymTORC1 + ymTORC1, 2)/std::pow(sigma_ymTORC1, 3);
            break;
        case 655:
            dJydsigma[655] = 1.0/sigma_ymTORC1active - 1.0*std::pow(-mymTORC1active + ymTORC1active, 2)/std::pow(sigma_ymTORC1active, 3);
            break;
        case 656:
            dJydsigma[656] = 1.0/sigma_yPIP - 1.0*std::pow(-myPIP + yPIP, 2)/std::pow(sigma_yPIP, 3);
            break;
        case 657:
            dJydsigma[657] = 1.0/sigma_yPI3P - 1.0*std::pow(-myPI3P + yPI3P, 2)/std::pow(sigma_yPI3P, 3);
            break;
        case 658:
            dJydsigma[658] = 1.0/sigma_yDAG - 1.0*std::pow(-myDAG + yDAG, 2)/std::pow(sigma_yDAG, 3);
            break;
        case 659:
            dJydsigma[659] = 1.0/sigma_yGRP - 1.0*std::pow(-myGRP + yGRP, 2)/std::pow(sigma_yGRP, 3);
            break;
        case 660:
            dJydsigma[660] = 1.0/sigma_yDAG_GRP - 1.0*std::pow(-myDAG_GRP + yDAG_GRP, 2)/std::pow(sigma_yDAG_GRP, 3);
            break;
        case 661:
            dJydsigma[661] = 1.0/sigma_yRasT - 1.0*std::pow(-myRasT + yRasT, 2)/std::pow(sigma_yRasT, 3);
            break;
        case 662:
            dJydsigma[662] = 1.0/sigma_yRasD - 1.0*std::pow(-myRasD + yRasD, 2)/std::pow(sigma_yRasD, 3);
            break;
        case 663:
            dJydsigma[663] = 1.0/sigma_yNF1 - 1.0*std::pow(-myNF1 + yNF1, 2)/std::pow(sigma_yNF1, 3);
            break;
        case 664:
            dJydsigma[664] = 1.0/sigma_ypNF1 - 1.0*std::pow(-mypNF1 + ypNF1, 2)/std::pow(sigma_ypNF1, 3);
            break;
        case 665:
            dJydsigma[665] = 1.0/sigma_ypCRaf - 1.0*std::pow(-mypCRaf + ypCRaf, 2)/std::pow(sigma_ypCRaf, 3);
            break;
        case 666:
            dJydsigma[666] = 1.0/sigma_yCRaf - 1.0*std::pow(-myCRaf + yCRaf, 2)/std::pow(sigma_yCRaf, 3);
            break;
        case 667:
            dJydsigma[667] = 1.0/sigma_yRasT_CRaf - 1.0*std::pow(-myRasT_CRaf + yRasT_CRaf, 2)/std::pow(sigma_yRasT_CRaf, 3);
            break;
        case 668:
            dJydsigma[668] = 1.0/sigma_yBRaf - 1.0*std::pow(-myBRaf + yBRaf, 2)/std::pow(sigma_yBRaf, 3);
            break;
        case 669:
            dJydsigma[669] = 1.0/sigma_yRasT_CRaf_BRaf - 1.0*std::pow(-myRasT_CRaf_BRaf + yRasT_CRaf_BRaf, 2)/std::pow(sigma_yRasT_CRaf_BRaf, 3);
            break;
        case 670:
            dJydsigma[670] = 1.0/sigma_yMEK - 1.0*std::pow(-myMEK + yMEK, 2)/std::pow(sigma_yMEK, 3);
            break;
        case 671:
            dJydsigma[671] = 1.0/sigma_ypMEK - 1.0*std::pow(-mypMEK + ypMEK, 2)/std::pow(sigma_ypMEK, 3);
            break;
        case 672:
            dJydsigma[672] = 1.0/sigma_yppMEK - 1.0*std::pow(-myppMEK + yppMEK, 2)/std::pow(sigma_yppMEK, 3);
            break;
        case 673:
            dJydsigma[673] = 1.0/sigma_yMKP3 - 1.0*std::pow(-myMKP3 + yMKP3, 2)/std::pow(sigma_yMKP3, 3);
            break;
        case 674:
            dJydsigma[674] = 1.0/sigma_yERKnuc - 1.0*std::pow(-myERKnuc + yERKnuc, 2)/std::pow(sigma_yERKnuc, 3);
            break;
        case 675:
            dJydsigma[675] = 1.0/sigma_yppERKnuc - 1.0*std::pow(-myppERKnuc + yppERKnuc, 2)/std::pow(sigma_yppERKnuc, 3);
            break;
        case 676:
            dJydsigma[676] = 1.0/sigma_yRSK - 1.0*std::pow(-myRSK + yRSK, 2)/std::pow(sigma_yRSK, 3);
            break;
        case 677:
            dJydsigma[677] = 1.0/sigma_ypRSK - 1.0*std::pow(-mypRSK + ypRSK, 2)/std::pow(sigma_ypRSK, 3);
            break;
        case 678:
            dJydsigma[678] = 1.0/sigma_ypRSKnuc - 1.0*std::pow(-mypRSKnuc + ypRSKnuc, 2)/std::pow(sigma_ypRSKnuc, 3);
            break;
        case 679:
            dJydsigma[679] = 1.0/sigma_yMKP1 - 1.0*std::pow(-myMKP1 + yMKP1, 2)/std::pow(sigma_yMKP1, 3);
            break;
        case 680:
            dJydsigma[680] = 1.0/sigma_ypMKP1 - 1.0*std::pow(-mypMKP1 + ypMKP1, 2)/std::pow(sigma_ypMKP1, 3);
            break;
        case 681:
            dJydsigma[681] = 1.0/sigma_ycFos - 1.0*std::pow(-mycFos + ycFos, 2)/std::pow(sigma_ycFos, 3);
            break;
        case 682:
            dJydsigma[682] = 1.0/sigma_ypcFos - 1.0*std::pow(-mypcFos + ypcFos, 2)/std::pow(sigma_ypcFos, 3);
            break;
        case 683:
            dJydsigma[683] = 1.0/sigma_ycJun - 1.0*std::pow(-mycJun + ycJun, 2)/std::pow(sigma_ycJun, 3);
            break;
        case 684:
            dJydsigma[684] = 1.0/sigma_ypcFos_cJun - 1.0*std::pow(-mypcFos_cJun + ypcFos_cJun, 2)/std::pow(sigma_ypcFos_cJun, 3);
            break;
        case 685:
            dJydsigma[685] = 1.0/sigma_ycMyc - 1.0*std::pow(-mycMyc + ycMyc, 2)/std::pow(sigma_ycMyc, 3);
            break;
        case 686:
            dJydsigma[686] = 1.0/sigma_ybCATENINnuc - 1.0*std::pow(-mybCATENINnuc + ybCATENINnuc, 2)/std::pow(sigma_ybCATENINnuc, 3);
            break;
        case 687:
            dJydsigma[687] = 1.0/sigma_ybCATENIN - 1.0*std::pow(-mybCATENIN + ybCATENIN, 2)/std::pow(sigma_ybCATENIN, 3);
            break;
        case 688:
            dJydsigma[688] = 1.0/sigma_ypbCATENIN - 1.0*std::pow(-mypbCATENIN + ypbCATENIN, 2)/std::pow(sigma_ypbCATENIN, 3);
            break;
        case 689:
            dJydsigma[689] = 1.0/sigma_yIP3 - 1.0*std::pow(-myIP3 + yIP3, 2)/std::pow(sigma_yIP3, 3);
            break;
        case 690:
            dJydsigma[690] = 1.0/sigma_yPIP2 - 1.0*std::pow(-myPIP2 + yPIP2, 2)/std::pow(sigma_yPIP2, 3);
            break;
        case 691:
            dJydsigma[691] = 1.0/sigma_yPIP3 - 1.0*std::pow(-myPIP3 + yPIP3, 2)/std::pow(sigma_yPIP3, 3);
            break;
        case 692:
            dJydsigma[692] = 1.0/sigma_yPTEN - 1.0*std::pow(-myPTEN + yPTEN, 2)/std::pow(sigma_yPTEN, 3);
            break;
        case 693:
            dJydsigma[693] = 1.0/sigma_yPIP3_AKT - 1.0*std::pow(-myPIP3_AKT + yPIP3_AKT, 2)/std::pow(sigma_yPIP3_AKT, 3);
            break;
        case 694:
            dJydsigma[694] = 1.0/sigma_yAKT - 1.0*std::pow(-myAKT + yAKT, 2)/std::pow(sigma_yAKT, 3);
            break;
        case 695:
            dJydsigma[695] = 1.0/sigma_ypAKT - 1.0*std::pow(-mypAKT + ypAKT, 2)/std::pow(sigma_ypAKT, 3);
            break;
        case 696:
            dJydsigma[696] = 1.0/sigma_yppAKT - 1.0*std::pow(-myppAKT + yppAKT, 2)/std::pow(sigma_yppAKT, 3);
            break;
        case 697:
            dJydsigma[697] = 1.0/sigma_yPDK1 - 1.0*std::pow(-myPDK1 + yPDK1, 2)/std::pow(sigma_yPDK1, 3);
            break;
        case 698:
            dJydsigma[698] = 1.0/sigma_yPIP3_PDK1 - 1.0*std::pow(-myPIP3_PDK1 + yPIP3_PDK1, 2)/std::pow(sigma_yPIP3_PDK1, 3);
            break;
        case 699:
            dJydsigma[699] = 1.0/sigma_yPIP3_pAKT - 1.0*std::pow(-myPIP3_pAKT + yPIP3_pAKT, 2)/std::pow(sigma_yPIP3_pAKT, 3);
            break;
        case 700:
            dJydsigma[700] = 1.0/sigma_yRictor - 1.0*std::pow(-myRictor + yRictor, 2)/std::pow(sigma_yRictor, 3);
            break;
        case 701:
            dJydsigma[701] = 1.0/sigma_ymTOR - 1.0*std::pow(-mymTOR + ymTOR, 2)/std::pow(sigma_ymTOR, 3);
            break;
        case 702:
            dJydsigma[702] = 1.0/sigma_ymTORC2 - 1.0*std::pow(-mymTORC2 + ymTORC2, 2)/std::pow(sigma_ymTORC2, 3);
            break;
        case 703:
            dJydsigma[703] = 1.0/sigma_yPIP3_ppAKT - 1.0*std::pow(-myPIP3_ppAKT + yPIP3_ppAKT, 2)/std::pow(sigma_yPIP3_ppAKT, 3);
            break;
        case 704:
            dJydsigma[704] = 1.0/sigma_yGSK3b - 1.0*std::pow(-myGSK3b + yGSK3b, 2)/std::pow(sigma_yGSK3b, 3);
            break;
        case 705:
            dJydsigma[705] = 1.0/sigma_ypGSK3b - 1.0*std::pow(-mypGSK3b + ypGSK3b, 2)/std::pow(sigma_ypGSK3b, 3);
            break;
        case 706:
            dJydsigma[706] = 1.0/sigma_yTSC1 - 1.0*std::pow(-myTSC1 + yTSC1, 2)/std::pow(sigma_yTSC1, 3);
            break;
        case 707:
            dJydsigma[707] = 1.0/sigma_yTSC2 - 1.0*std::pow(-myTSC2 + yTSC2, 2)/std::pow(sigma_yTSC2, 3);
            break;
        case 708:
            dJydsigma[708] = 1.0/sigma_ypTSC2 - 1.0*std::pow(-mypTSC2 + ypTSC2, 2)/std::pow(sigma_ypTSC2, 3);
            break;
        case 709:
            dJydsigma[709] = 1.0/sigma_yTSC - 1.0*std::pow(-myTSC + yTSC, 2)/std::pow(sigma_yTSC, 3);
            break;
        case 710:
            dJydsigma[710] = 1.0/sigma_yPKC - 1.0*std::pow(-myPKC + yPKC, 2)/std::pow(sigma_yPKC, 3);
            break;
        case 711:
            dJydsigma[711] = 1.0/sigma_yDAG_PKC - 1.0*std::pow(-myDAG_PKC + yDAG_PKC, 2)/std::pow(sigma_yDAG_PKC, 3);
            break;
        case 712:
            dJydsigma[712] = 1.0/sigma_ypRKIP - 1.0*std::pow(-mypRKIP + ypRKIP, 2)/std::pow(sigma_ypRKIP, 3);
            break;
        case 713:
            dJydsigma[713] = 1.0/sigma_yRKIP - 1.0*std::pow(-myRKIP + yRKIP, 2)/std::pow(sigma_yRKIP, 3);
            break;
        case 714:
            dJydsigma[714] = 1.0/sigma_yRKIP_CRaf - 1.0*std::pow(-myRKIP_CRaf + yRKIP_CRaf, 2)/std::pow(sigma_yRKIP_CRaf, 3);
            break;
        case 715:
            dJydsigma[715] = 1.0/sigma_yERK - 1.0*std::pow(-myERK + yERK, 2)/std::pow(sigma_yERK, 3);
            break;
        case 716:
            dJydsigma[716] = 1.0/sigma_ypERK - 1.0*std::pow(-mypERK + ypERK, 2)/std::pow(sigma_ypERK, 3);
            break;
        case 717:
            dJydsigma[717] = 1.0/sigma_yppERK - 1.0*std::pow(-myppERK + yppERK, 2)/std::pow(sigma_yppERK, 3);
            break;
        case 718:
            dJydsigma[718] = 1.0/sigma_yFOXO - 1.0*std::pow(-myFOXO + yFOXO, 2)/std::pow(sigma_yFOXO, 3);
            break;
        case 719:
            dJydsigma[719] = 1.0/sigma_ypFOXO - 1.0*std::pow(-mypFOXO + ypFOXO, 2)/std::pow(sigma_ypFOXO, 3);
            break;
        case 720:
            dJydsigma[720] = 1.0/sigma_yRhebD - 1.0*std::pow(-myRhebD + yRhebD, 2)/std::pow(sigma_yRhebD, 3);
            break;
        case 721:
            dJydsigma[721] = 1.0/sigma_yRhebT - 1.0*std::pow(-myRhebT + yRhebT, 2)/std::pow(sigma_yRhebT, 3);
            break;
        case 722:
            dJydsigma[722] = 1.0/sigma_yRaptor - 1.0*std::pow(-myRaptor + yRaptor, 2)/std::pow(sigma_yRaptor, 3);
            break;
        case 723:
            dJydsigma[723] = 1.0/sigma_yS6K - 1.0*std::pow(-myS6K + yS6K, 2)/std::pow(sigma_yS6K, 3);
            break;
        case 724:
            dJydsigma[724] = 1.0/sigma_ypS6K - 1.0*std::pow(-mypS6K + ypS6K, 2)/std::pow(sigma_ypS6K, 3);
            break;
        case 725:
            dJydsigma[725] = 1.0/sigma_yEIF4EBP1 - 1.0*std::pow(-myEIF4EBP1 + yEIF4EBP1, 2)/std::pow(sigma_yEIF4EBP1, 3);
            break;
        case 726:
            dJydsigma[726] = 1.0/sigma_ypEIF4EBP1 - 1.0*std::pow(-mypEIF4EBP1 + ypEIF4EBP1, 2)/std::pow(sigma_ypEIF4EBP1, 3);
            break;
        case 727:
            dJydsigma[727] = 1.0/sigma_ySOS - 1.0*std::pow(-mySOS + ySOS, 2)/std::pow(sigma_ySOS, 3);
            break;
        case 728:
            dJydsigma[728] = 1.0/sigma_yG2_SOS_ppERK - 1.0*std::pow(-myG2_SOS_ppERK + yG2_SOS_ppERK, 2)/std::pow(sigma_yG2_SOS_ppERK, 3);
            break;
        case 729:
            dJydsigma[729] = 1.0/sigma_yCRaf_ppERK - 1.0*std::pow(-myCRaf_ppERK + yCRaf_ppERK, 2)/std::pow(sigma_yCRaf_ppERK, 3);
            break;
        case 730:
            dJydsigma[730] = 1.0/sigma_yRasD_DAG_GRP - 1.0*std::pow(-myRasD_DAG_GRP + yRasD_DAG_GRP, 2)/std::pow(sigma_yRasD_DAG_GRP, 3);
            break;
        case 731:
            dJydsigma[731] = 1.0/sigma_yRasT_NF1 - 1.0*std::pow(-myRasT_NF1 + yRasT_NF1, 2)/std::pow(sigma_yRasT_NF1, 3);
            break;
        case 732:
            dJydsigma[732] = 1.0/sigma_yNF1_ppERK - 1.0*std::pow(-myNF1_ppERK + yNF1_ppERK, 2)/std::pow(sigma_yNF1_ppERK, 3);
            break;
        case 733:
            dJydsigma[733] = 1.0/sigma_yMEK_RasT_CRaf_BRaf - 1.0*std::pow(-myMEK_RasT_CRaf_BRaf + yMEK_RasT_CRaf_BRaf, 2)/std::pow(sigma_yMEK_RasT_CRaf_BRaf, 3);
            break;
        case 734:
            dJydsigma[734] = 1.0/sigma_ypMEK_RasT_CRaf_BRaf - 1.0*std::pow(-mypMEK_RasT_CRaf_BRaf + ypMEK_RasT_CRaf_BRaf, 2)/std::pow(sigma_ypMEK_RasT_CRaf_BRaf, 3);
            break;
        case 735:
            dJydsigma[735] = 1.0/sigma_yERK_ppMEK - 1.0*std::pow(-myERK_ppMEK + yERK_ppMEK, 2)/std::pow(sigma_yERK_ppMEK, 3);
            break;
        case 736:
            dJydsigma[736] = 1.0/sigma_ypERK_ppMEK - 1.0*std::pow(-mypERK_ppMEK + ypERK_ppMEK, 2)/std::pow(sigma_ypERK_ppMEK, 3);
            break;
        case 737:
            dJydsigma[737] = 1.0/sigma_yRSK_ppERK - 1.0*std::pow(-myRSK_ppERK + yRSK_ppERK, 2)/std::pow(sigma_yRSK_ppERK, 3);
            break;
        case 738:
            dJydsigma[738] = 1.0/sigma_ypRSKnuc_MKP1 - 1.0*std::pow(-mypRSKnuc_MKP1 + ypRSKnuc_MKP1, 2)/std::pow(sigma_ypRSKnuc_MKP1, 3);
            break;
        case 739:
            dJydsigma[739] = 1.0/sigma_yppERKnuc_MKP1 - 1.0*std::pow(-myppERKnuc_MKP1 + yppERKnuc_MKP1, 2)/std::pow(sigma_yppERKnuc_MKP1, 3);
            break;
        case 740:
            dJydsigma[740] = 1.0/sigma_ycFos_pRSKnuc - 1.0*std::pow(-mycFos_pRSKnuc + ycFos_pRSKnuc, 2)/std::pow(sigma_ycFos_pRSKnuc, 3);
            break;
        case 741:
            dJydsigma[741] = 1.0/sigma_ycFos_ppERKnuc - 1.0*std::pow(-mycFos_ppERKnuc + ycFos_ppERKnuc, 2)/std::pow(sigma_ycFos_ppERKnuc, 3);
            break;
        case 742:
            dJydsigma[742] = 1.0/sigma_yRKIP_DAG_PKC - 1.0*std::pow(-myRKIP_DAG_PKC + yRKIP_DAG_PKC, 2)/std::pow(sigma_yRKIP_DAG_PKC, 3);
            break;
        case 743:
            dJydsigma[743] = 1.0/sigma_yPIP3_PTEN - 1.0*std::pow(-myPIP3_PTEN + yPIP3_PTEN, 2)/std::pow(sigma_yPIP3_PTEN, 3);
            break;
        case 744:
            dJydsigma[744] = 1.0/sigma_yPIP3_AKT_PIP3_PDK1 - 1.0*std::pow(-myPIP3_AKT_PIP3_PDK1 + yPIP3_AKT_PIP3_PDK1, 2)/std::pow(sigma_yPIP3_AKT_PIP3_PDK1, 3);
            break;
        case 745:
            dJydsigma[745] = 1.0/sigma_yPIP3_pAKT_mTORC2 - 1.0*std::pow(-myPIP3_pAKT_mTORC2 + yPIP3_pAKT_mTORC2, 2)/std::pow(sigma_yPIP3_pAKT_mTORC2, 3);
            break;
        case 746:
            dJydsigma[746] = 1.0/sigma_yGSK3b_ppAKT - 1.0*std::pow(-myGSK3b_ppAKT + yGSK3b_ppAKT, 2)/std::pow(sigma_yGSK3b_ppAKT, 3);
            break;
        case 747:
            dJydsigma[747] = 1.0/sigma_yTSC2_ppAKT - 1.0*std::pow(-myTSC2_ppAKT + yTSC2_ppAKT, 2)/std::pow(sigma_yTSC2_ppAKT, 3);
            break;
        case 748:
            dJydsigma[748] = 1.0/sigma_yTSC2_ppERK - 1.0*std::pow(-myTSC2_ppERK + yTSC2_ppERK, 2)/std::pow(sigma_yTSC2_ppERK, 3);
            break;
        case 749:
            dJydsigma[749] = 1.0/sigma_yRhebT_TSC - 1.0*std::pow(-myRhebT_TSC + yRhebT_TSC, 2)/std::pow(sigma_yRhebT_TSC, 3);
            break;
        case 750:
            dJydsigma[750] = 1.0/sigma_yEIF4EBP1_mTORC1active - 1.0*std::pow(-myEIF4EBP1_mTORC1active + yEIF4EBP1_mTORC1active, 2)/std::pow(sigma_yEIF4EBP1_mTORC1active, 3);
            break;
        case 751:
            dJydsigma[751] = 1.0/sigma_yS6K_mTORC1active - 1.0*std::pow(-myS6K_mTORC1active + yS6K_mTORC1active, 2)/std::pow(sigma_yS6K_mTORC1active, 3);
            break;
        case 752:
            dJydsigma[752] = 1.0/sigma_yFOXO_ppAKT - 1.0*std::pow(-myFOXO_ppAKT + yFOXO_ppAKT, 2)/std::pow(sigma_yFOXO_ppAKT, 3);
            break;
        case 753:
            dJydsigma[753] = 1.0/sigma_yPI3K1_mTORC1active - 1.0*std::pow(-myPI3K1_mTORC1active + yPI3K1_mTORC1active, 2)/std::pow(sigma_yPI3K1_mTORC1active, 3);
            break;
        case 754:
            dJydsigma[754] = 1.0/sigma_ypERK_MKP3 - 1.0*std::pow(-mypERK_MKP3 + ypERK_MKP3, 2)/std::pow(sigma_ypERK_MKP3, 3);
            break;
        case 755:
            dJydsigma[755] = 1.0/sigma_yppERK_MKP3 - 1.0*std::pow(-myppERK_MKP3 + yppERK_MKP3, 2)/std::pow(sigma_yppERK_MKP3, 3);
            break;
        case 756:
            dJydsigma[756] = 1.0/sigma_yppERKnuc_pMKP1 - 1.0*std::pow(-myppERKnuc_pMKP1 + yppERKnuc_pMKP1, 2)/std::pow(sigma_yppERKnuc_pMKP1, 3);
            break;
        case 757:
            dJydsigma[757] = 1.0/sigma_yRasT_BRaf - 1.0*std::pow(-myRasT_BRaf + yRasT_BRaf, 2)/std::pow(sigma_yRasT_BRaf, 3);
            break;
        case 758:
            dJydsigma[758] = 1.0/sigma_yRasT_BRaf_BRaf - 1.0*std::pow(-myRasT_BRaf_BRaf + yRasT_BRaf_BRaf, 2)/std::pow(sigma_yRasT_BRaf_BRaf, 3);
            break;
        case 759:
            dJydsigma[759] = 1.0/sigma_yMEK_RasT_BRaf_BRaf - 1.0*std::pow(-myMEK_RasT_BRaf_BRaf + yMEK_RasT_BRaf_BRaf, 2)/std::pow(sigma_yMEK_RasT_BRaf_BRaf, 3);
            break;
        case 760:
            dJydsigma[760] = 1.0/sigma_ypMEK_RasT_BRaf_BRaf - 1.0*std::pow(-mypMEK_RasT_BRaf_BRaf + ypMEK_RasT_BRaf_BRaf, 2)/std::pow(sigma_ypMEK_RasT_BRaf_BRaf, 3);
            break;
        case 761:
            dJydsigma[761] = 1.0/sigma_yEIF4E - 1.0*std::pow(-myEIF4E + yEIF4E, 2)/std::pow(sigma_yEIF4E, 3);
            break;
        case 762:
            dJydsigma[762] = 1.0/sigma_yEIF4EBP1_EIF4E - 1.0*std::pow(-myEIF4EBP1_EIF4E + yEIF4EBP1_EIF4E, 2)/std::pow(sigma_yEIF4EBP1_EIF4E, 3);
            break;
        case 763:
            dJydsigma[763] = 1.0/sigma_yRasT_CRaf_CRaf - 1.0*std::pow(-myRasT_CRaf_CRaf + yRasT_CRaf_CRaf, 2)/std::pow(sigma_yRasT_CRaf_CRaf, 3);
            break;
        case 764:
            dJydsigma[764] = 1.0/sigma_yMEK_RasT_CRaf_CRaf - 1.0*std::pow(-myMEK_RasT_CRaf_CRaf + yMEK_RasT_CRaf_CRaf, 2)/std::pow(sigma_yMEK_RasT_CRaf_CRaf, 3);
            break;
        case 765:
            dJydsigma[765] = 1.0/sigma_ypMEK_RasT_CRaf_CRaf - 1.0*std::pow(-mypMEK_RasT_CRaf_CRaf + ypMEK_RasT_CRaf_CRaf, 2)/std::pow(sigma_ypMEK_RasT_CRaf_CRaf, 3);
            break;
        case 766:
            dJydsigma[766] = 1.0/sigma_yFOXOnuc - 1.0*std::pow(-myFOXOnuc + yFOXOnuc, 2)/std::pow(sigma_yFOXOnuc, 3);
            break;
        case 767:
            dJydsigma[767] = 1.0/sigma_yMEKi - 1.0*std::pow(-myMEKi + yMEKi, 2)/std::pow(sigma_yMEKi, 3);
            break;
        case 768:
            dJydsigma[768] = 1.0/sigma_yMEKi_ppMEK - 1.0*std::pow(-myMEKi_ppMEK + yMEKi_ppMEK, 2)/std::pow(sigma_yMEKi_ppMEK, 3);
            break;
        case 769:
            dJydsigma[769] = 1.0/sigma_yAKTi - 1.0*std::pow(-myAKTi + yAKTi, 2)/std::pow(sigma_yAKTi, 3);
            break;
        case 770:
            dJydsigma[770] = 1.0/sigma_yAKTi_AKT - 1.0*std::pow(-myAKTi_AKT + yAKTi_AKT, 2)/std::pow(sigma_yAKTi_AKT, 3);
            break;
        case 771:
            dJydsigma[771] = 1.0/sigma_ymT - 1.0*std::pow(-mymT + ymT, 2)/std::pow(sigma_ymT, 3);
            break;
        case 772:
            dJydsigma[772] = 1.0/sigma_yEIF4E_mT - 1.0*std::pow(-myEIF4E_mT + yEIF4E_mT, 2)/std::pow(sigma_yEIF4E_mT, 3);
            break;
        case 773:
            dJydsigma[773] = 1.0/sigma_ym_TP53 - 1.0*std::pow(-mym_TP53 + ym_TP53, 2)/std::pow(sigma_ym_TP53, 3);
            break;
        case 774:
            dJydsigma[774] = 1.0/sigma_ym_MDM2 - 1.0*std::pow(-mym_MDM2 + ym_MDM2, 2)/std::pow(sigma_ym_MDM2, 3);
            break;
        case 775:
            dJydsigma[775] = 1.0/sigma_ym_PPM1D - 1.0*std::pow(-mym_PPM1D + ym_PPM1D, 2)/std::pow(sigma_ym_PPM1D, 3);
            break;
        case 776:
            dJydsigma[776] = 1.0/sigma_ym_ATM - 1.0*std::pow(-mym_ATM + ym_ATM, 2)/std::pow(sigma_ym_ATM, 3);
            break;
        case 777:
            dJydsigma[777] = 1.0/sigma_ym_ATR - 1.0*std::pow(-mym_ATR + ym_ATR, 2)/std::pow(sigma_ym_ATR, 3);
            break;
        case 778:
            dJydsigma[778] = 1.0/sigma_ym_RB1 - 1.0*std::pow(-mym_RB1 + ym_RB1, 2)/std::pow(sigma_ym_RB1, 3);
            break;
        case 779:
            dJydsigma[779] = 1.0/sigma_ym_E2F1 - 1.0*std::pow(-mym_E2F1 + ym_E2F1, 2)/std::pow(sigma_ym_E2F1, 3);
            break;
        case 780:
            dJydsigma[780] = 1.0/sigma_ym_E2F2 - 1.0*std::pow(-mym_E2F2 + ym_E2F2, 2)/std::pow(sigma_ym_E2F2, 3);
            break;
        case 781:
            dJydsigma[781] = 1.0/sigma_ym_E2F3 - 1.0*std::pow(-mym_E2F3 + ym_E2F3, 2)/std::pow(sigma_ym_E2F3, 3);
            break;
        case 782:
            dJydsigma[782] = 1.0/sigma_ym_CCND1 - 1.0*std::pow(-mym_CCND1 + ym_CCND1, 2)/std::pow(sigma_ym_CCND1, 3);
            break;
        case 783:
            dJydsigma[783] = 1.0/sigma_ym_CCND2 - 1.0*std::pow(-mym_CCND2 + ym_CCND2, 2)/std::pow(sigma_ym_CCND2, 3);
            break;
        case 784:
            dJydsigma[784] = 1.0/sigma_ym_CCND3 - 1.0*std::pow(-mym_CCND3 + ym_CCND3, 2)/std::pow(sigma_ym_CCND3, 3);
            break;
        case 785:
            dJydsigma[785] = 1.0/sigma_ym_CCNE1 - 1.0*std::pow(-mym_CCNE1 + ym_CCNE1, 2)/std::pow(sigma_ym_CCNE1, 3);
            break;
        case 786:
            dJydsigma[786] = 1.0/sigma_ym_CCNE2 - 1.0*std::pow(-mym_CCNE2 + ym_CCNE2, 2)/std::pow(sigma_ym_CCNE2, 3);
            break;
        case 787:
            dJydsigma[787] = 1.0/sigma_ym_SKP2 - 1.0*std::pow(-mym_SKP2 + ym_SKP2, 2)/std::pow(sigma_ym_SKP2, 3);
            break;
        case 788:
            dJydsigma[788] = 1.0/sigma_ym_CDC25A - 1.0*std::pow(-mym_CDC25A + ym_CDC25A, 2)/std::pow(sigma_ym_CDC25A, 3);
            break;
        case 789:
            dJydsigma[789] = 1.0/sigma_ym_CDC25B - 1.0*std::pow(-mym_CDC25B + ym_CDC25B, 2)/std::pow(sigma_ym_CDC25B, 3);
            break;
        case 790:
            dJydsigma[790] = 1.0/sigma_ym_CDC25C - 1.0*std::pow(-mym_CDC25C + ym_CDC25C, 2)/std::pow(sigma_ym_CDC25C, 3);
            break;
        case 791:
            dJydsigma[791] = 1.0/sigma_ym_CCNA2 - 1.0*std::pow(-mym_CCNA2 + ym_CCNA2, 2)/std::pow(sigma_ym_CCNA2, 3);
            break;
        case 792:
            dJydsigma[792] = 1.0/sigma_ym_CDKN1B - 1.0*std::pow(-mym_CDKN1B + ym_CDKN1B, 2)/std::pow(sigma_ym_CDKN1B, 3);
            break;
        case 793:
            dJydsigma[793] = 1.0/sigma_ym_CDH1 - 1.0*std::pow(-mym_CDH1 + ym_CDH1, 2)/std::pow(sigma_ym_CDH1, 3);
            break;
        case 794:
            dJydsigma[794] = 1.0/sigma_ym_CCNB1 - 1.0*std::pow(-mym_CCNB1 + ym_CCNB1, 2)/std::pow(sigma_ym_CCNB1, 3);
            break;
        case 795:
            dJydsigma[795] = 1.0/sigma_ym_CDC20 - 1.0*std::pow(-mym_CDC20 + ym_CDC20, 2)/std::pow(sigma_ym_CDC20, 3);
            break;
        case 796:
            dJydsigma[796] = 1.0/sigma_ym_WEE1 - 1.0*std::pow(-mym_WEE1 + ym_WEE1, 2)/std::pow(sigma_ym_WEE1, 3);
            break;
        case 797:
            dJydsigma[797] = 1.0/sigma_ym_CHEK1 - 1.0*std::pow(-mym_CHEK1 + ym_CHEK1, 2)/std::pow(sigma_ym_CHEK1, 3);
            break;
        case 798:
            dJydsigma[798] = 1.0/sigma_ym_CDKN1A - 1.0*std::pow(-mym_CDKN1A + ym_CDKN1A, 2)/std::pow(sigma_ym_CDKN1A, 3);
            break;
        case 799:
            dJydsigma[799] = 1.0/sigma_ym_CDK1 - 1.0*std::pow(-mym_CDK1 + ym_CDK1, 2)/std::pow(sigma_ym_CDK1, 3);
            break;
        case 800:
            dJydsigma[800] = 1.0/sigma_ym_CDK2 - 1.0*std::pow(-mym_CDK2 + ym_CDK2, 2)/std::pow(sigma_ym_CDK2, 3);
            break;
        case 801:
            dJydsigma[801] = 1.0/sigma_ym_CDK4 - 1.0*std::pow(-mym_CDK4 + ym_CDK4, 2)/std::pow(sigma_ym_CDK4, 3);
            break;
        case 802:
            dJydsigma[802] = 1.0/sigma_ym_CDK6 - 1.0*std::pow(-mym_CDK6 + ym_CDK6, 2)/std::pow(sigma_ym_CDK6, 3);
            break;
        case 803:
            dJydsigma[803] = 1.0/sigma_ym_TNFSF10 - 1.0*std::pow(-mym_TNFSF10 + ym_TNFSF10, 2)/std::pow(sigma_ym_TNFSF10, 3);
            break;
        case 804:
            dJydsigma[804] = 1.0/sigma_ym_TNFRSF10A - 1.0*std::pow(-mym_TNFRSF10A + ym_TNFRSF10A, 2)/std::pow(sigma_ym_TNFRSF10A, 3);
            break;
        case 805:
            dJydsigma[805] = 1.0/sigma_ym_TNFRSF10B - 1.0*std::pow(-mym_TNFRSF10B + ym_TNFRSF10B, 2)/std::pow(sigma_ym_TNFRSF10B, 3);
            break;
        case 806:
            dJydsigma[806] = 1.0/sigma_ym_CFLAR - 1.0*std::pow(-mym_CFLAR + ym_CFLAR, 2)/std::pow(sigma_ym_CFLAR, 3);
            break;
        case 807:
            dJydsigma[807] = 1.0/sigma_ym_CASP8 - 1.0*std::pow(-mym_CASP8 + ym_CASP8, 2)/std::pow(sigma_ym_CASP8, 3);
            break;
        case 808:
            dJydsigma[808] = 1.0/sigma_ym_CASP10 - 1.0*std::pow(-mym_CASP10 + ym_CASP10, 2)/std::pow(sigma_ym_CASP10, 3);
            break;
        case 809:
            dJydsigma[809] = 1.0/sigma_ym_BFAR - 1.0*std::pow(-mym_BFAR + ym_BFAR, 2)/std::pow(sigma_ym_BFAR, 3);
            break;
        case 810:
            dJydsigma[810] = 1.0/sigma_ym_CASP3 - 1.0*std::pow(-mym_CASP3 + ym_CASP3, 2)/std::pow(sigma_ym_CASP3, 3);
            break;
        case 811:
            dJydsigma[811] = 1.0/sigma_ym_CASP7 - 1.0*std::pow(-mym_CASP7 + ym_CASP7, 2)/std::pow(sigma_ym_CASP7, 3);
            break;
        case 812:
            dJydsigma[812] = 1.0/sigma_ym_CASP6 - 1.0*std::pow(-mym_CASP6 + ym_CASP6, 2)/std::pow(sigma_ym_CASP6, 3);
            break;
        case 813:
            dJydsigma[813] = 1.0/sigma_ym_XIAP - 1.0*std::pow(-mym_XIAP + ym_XIAP, 2)/std::pow(sigma_ym_XIAP, 3);
            break;
        case 814:
            dJydsigma[814] = 1.0/sigma_ym_PARP1 - 1.0*std::pow(-mym_PARP1 + ym_PARP1, 2)/std::pow(sigma_ym_PARP1, 3);
            break;
        case 815:
            dJydsigma[815] = 1.0/sigma_ym_BID - 1.0*std::pow(-mym_BID + ym_BID, 2)/std::pow(sigma_ym_BID, 3);
            break;
        case 816:
            dJydsigma[816] = 1.0/sigma_ym_BCL2 - 1.0*std::pow(-mym_BCL2 + ym_BCL2, 2)/std::pow(sigma_ym_BCL2, 3);
            break;
        case 817:
            dJydsigma[817] = 1.0/sigma_ym_BCL2L1 - 1.0*std::pow(-mym_BCL2L1 + ym_BCL2L1, 2)/std::pow(sigma_ym_BCL2L1, 3);
            break;
        case 818:
            dJydsigma[818] = 1.0/sigma_ym_MCL1 - 1.0*std::pow(-mym_MCL1 + ym_MCL1, 2)/std::pow(sigma_ym_MCL1, 3);
            break;
        case 819:
            dJydsigma[819] = 1.0/sigma_ym_BAX - 1.0*std::pow(-mym_BAX + ym_BAX, 2)/std::pow(sigma_ym_BAX, 3);
            break;
        case 820:
            dJydsigma[820] = 1.0/sigma_ym_CYCS - 1.0*std::pow(-mym_CYCS + ym_CYCS, 2)/std::pow(sigma_ym_CYCS, 3);
            break;
        case 821:
            dJydsigma[821] = 1.0/sigma_ym_DIABLO - 1.0*std::pow(-mym_DIABLO + ym_DIABLO, 2)/std::pow(sigma_ym_DIABLO, 3);
            break;
        case 822:
            dJydsigma[822] = 1.0/sigma_ym_APAF1 - 1.0*std::pow(-mym_APAF1 + ym_APAF1, 2)/std::pow(sigma_ym_APAF1, 3);
            break;
        case 823:
            dJydsigma[823] = 1.0/sigma_ym_CASP9 - 1.0*std::pow(-mym_CASP9 + ym_CASP9, 2)/std::pow(sigma_ym_CASP9, 3);
            break;
        case 824:
            dJydsigma[824] = 1.0/sigma_ym_BAD - 1.0*std::pow(-mym_BAD + ym_BAD, 2)/std::pow(sigma_ym_BAD, 3);
            break;
        case 825:
            dJydsigma[825] = 1.0/sigma_ym_BBC3 - 1.0*std::pow(-mym_BBC3 + ym_BBC3, 2)/std::pow(sigma_ym_BBC3, 3);
            break;
        case 826:
            dJydsigma[826] = 1.0/sigma_ym_PMAIP1 - 1.0*std::pow(-mym_PMAIP1 + ym_PMAIP1, 2)/std::pow(sigma_ym_PMAIP1, 3);
            break;
        case 827:
            dJydsigma[827] = 1.0/sigma_ym_BCL2L11 - 1.0*std::pow(-mym_BCL2L11 + ym_BCL2L11, 2)/std::pow(sigma_ym_BCL2L11, 3);
            break;
        case 828:
            dJydsigma[828] = 1.0/sigma_ym_EGF - 1.0*std::pow(-mym_EGF + ym_EGF, 2)/std::pow(sigma_ym_EGF, 3);
            break;
        case 829:
            dJydsigma[829] = 1.0/sigma_ym_NRG1 - 1.0*std::pow(-mym_NRG1 + ym_NRG1, 2)/std::pow(sigma_ym_NRG1, 3);
            break;
        case 830:
            dJydsigma[830] = 1.0/sigma_ym_EGFR - 1.0*std::pow(-mym_EGFR + ym_EGFR, 2)/std::pow(sigma_ym_EGFR, 3);
            break;
        case 831:
            dJydsigma[831] = 1.0/sigma_ym_ERBB2 - 1.0*std::pow(-mym_ERBB2 + ym_ERBB2, 2)/std::pow(sigma_ym_ERBB2, 3);
            break;
        case 832:
            dJydsigma[832] = 1.0/sigma_ym_ERBB3 - 1.0*std::pow(-mym_ERBB3 + ym_ERBB3, 2)/std::pow(sigma_ym_ERBB3, 3);
            break;
        case 833:
            dJydsigma[833] = 1.0/sigma_ym_ERBB4 - 1.0*std::pow(-mym_ERBB4 + ym_ERBB4, 2)/std::pow(sigma_ym_ERBB4, 3);
            break;
        case 834:
            dJydsigma[834] = 1.0/sigma_ym_EGFRvIII - 1.0*std::pow(-mym_EGFRvIII + ym_EGFRvIII, 2)/std::pow(sigma_ym_EGFRvIII, 3);
            break;
        case 835:
            dJydsigma[835] = 1.0/sigma_ym_MET - 1.0*std::pow(-mym_MET + ym_MET, 2)/std::pow(sigma_ym_MET, 3);
            break;
        case 836:
            dJydsigma[836] = 1.0/sigma_ym_HGF - 1.0*std::pow(-mym_HGF + ym_HGF, 2)/std::pow(sigma_ym_HGF, 3);
            break;
        case 837:
            dJydsigma[837] = 1.0/sigma_ym_PDGFRA - 1.0*std::pow(-mym_PDGFRA + ym_PDGFRA, 2)/std::pow(sigma_ym_PDGFRA, 3);
            break;
        case 838:
            dJydsigma[838] = 1.0/sigma_ym_PDGFRB - 1.0*std::pow(-mym_PDGFRB + ym_PDGFRB, 2)/std::pow(sigma_ym_PDGFRB, 3);
            break;
        case 839:
            dJydsigma[839] = 1.0/sigma_ym_PDGFB - 1.0*std::pow(-mym_PDGFB + ym_PDGFB, 2)/std::pow(sigma_ym_PDGFB, 3);
            break;
        case 840:
            dJydsigma[840] = 1.0/sigma_ym_SPRY2 - 1.0*std::pow(-mym_SPRY2 + ym_SPRY2, 2)/std::pow(sigma_ym_SPRY2, 3);
            break;
        case 841:
            dJydsigma[841] = 1.0/sigma_ym_CBL - 1.0*std::pow(-mym_CBL + ym_CBL, 2)/std::pow(sigma_ym_CBL, 3);
            break;
        case 842:
            dJydsigma[842] = 1.0/sigma_ym_GRB2 - 1.0*std::pow(-mym_GRB2 + ym_GRB2, 2)/std::pow(sigma_ym_GRB2, 3);
            break;
        case 843:
            dJydsigma[843] = 1.0/sigma_ym_PLCG1 - 1.0*std::pow(-mym_PLCG1 + ym_PLCG1, 2)/std::pow(sigma_ym_PLCG1, 3);
            break;
        case 844:
            dJydsigma[844] = 1.0/sigma_ym_PLCG2 - 1.0*std::pow(-mym_PLCG2 + ym_PLCG2, 2)/std::pow(sigma_ym_PLCG2, 3);
            break;
        case 845:
            dJydsigma[845] = 1.0/sigma_ym_PIK3CA - 1.0*std::pow(-mym_PIK3CA + ym_PIK3CA, 2)/std::pow(sigma_ym_PIK3CA, 3);
            break;
        case 846:
            dJydsigma[846] = 1.0/sigma_ym_PIK3CB - 1.0*std::pow(-mym_PIK3CB + ym_PIK3CB, 2)/std::pow(sigma_ym_PIK3CB, 3);
            break;
        case 847:
            dJydsigma[847] = 1.0/sigma_ym_PIK3CG - 1.0*std::pow(-mym_PIK3CG + ym_PIK3CG, 2)/std::pow(sigma_ym_PIK3CG, 3);
            break;
        case 848:
            dJydsigma[848] = 1.0/sigma_ym_PIK3CD - 1.0*std::pow(-mym_PIK3CD + ym_PIK3CD, 2)/std::pow(sigma_ym_PIK3CD, 3);
            break;
        case 849:
            dJydsigma[849] = 1.0/sigma_ym_PIK3R1 - 1.0*std::pow(-mym_PIK3R1 + ym_PIK3R1, 2)/std::pow(sigma_ym_PIK3R1, 3);
            break;
        case 850:
            dJydsigma[850] = 1.0/sigma_ym_PIK3R2 - 1.0*std::pow(-mym_PIK3R2 + ym_PIK3R2, 2)/std::pow(sigma_ym_PIK3R2, 3);
            break;
        case 851:
            dJydsigma[851] = 1.0/sigma_ym_PIK3R3 - 1.0*std::pow(-mym_PIK3R3 + ym_PIK3R3, 2)/std::pow(sigma_ym_PIK3R3, 3);
            break;
        case 852:
            dJydsigma[852] = 1.0/sigma_ym_PIK3R4 - 1.0*std::pow(-mym_PIK3R4 + ym_PIK3R4, 2)/std::pow(sigma_ym_PIK3R4, 3);
            break;
        case 853:
            dJydsigma[853] = 1.0/sigma_ym_PIK3C2A - 1.0*std::pow(-mym_PIK3C2A + ym_PIK3C2A, 2)/std::pow(sigma_ym_PIK3C2A, 3);
            break;
        case 854:
            dJydsigma[854] = 1.0/sigma_ym_RASGRP1 - 1.0*std::pow(-mym_RASGRP1 + ym_RASGRP1, 2)/std::pow(sigma_ym_RASGRP1, 3);
            break;
        case 855:
            dJydsigma[855] = 1.0/sigma_ym_RASGRP3 - 1.0*std::pow(-mym_RASGRP3 + ym_RASGRP3, 2)/std::pow(sigma_ym_RASGRP3, 3);
            break;
        case 856:
            dJydsigma[856] = 1.0/sigma_ym_NRAS - 1.0*std::pow(-mym_NRAS + ym_NRAS, 2)/std::pow(sigma_ym_NRAS, 3);
            break;
        case 857:
            dJydsigma[857] = 1.0/sigma_ym_KRAS - 1.0*std::pow(-mym_KRAS + ym_KRAS, 2)/std::pow(sigma_ym_KRAS, 3);
            break;
        case 858:
            dJydsigma[858] = 1.0/sigma_ym_HRAS - 1.0*std::pow(-mym_HRAS + ym_HRAS, 2)/std::pow(sigma_ym_HRAS, 3);
            break;
        case 859:
            dJydsigma[859] = 1.0/sigma_ym_NF1 - 1.0*std::pow(-mym_NF1 + ym_NF1, 2)/std::pow(sigma_ym_NF1, 3);
            break;
        case 860:
            dJydsigma[860] = 1.0/sigma_ym_RAF1 - 1.0*std::pow(-mym_RAF1 + ym_RAF1, 2)/std::pow(sigma_ym_RAF1, 3);
            break;
        case 861:
            dJydsigma[861] = 1.0/sigma_ym_BRAF - 1.0*std::pow(-mym_BRAF + ym_BRAF, 2)/std::pow(sigma_ym_BRAF, 3);
            break;
        case 862:
            dJydsigma[862] = 1.0/sigma_ym_MAP2K1 - 1.0*std::pow(-mym_MAP2K1 + ym_MAP2K1, 2)/std::pow(sigma_ym_MAP2K1, 3);
            break;
        case 863:
            dJydsigma[863] = 1.0/sigma_ym_MAP2K2 - 1.0*std::pow(-mym_MAP2K2 + ym_MAP2K2, 2)/std::pow(sigma_ym_MAP2K2, 3);
            break;
        case 864:
            dJydsigma[864] = 1.0/sigma_ym_DUSP6 - 1.0*std::pow(-mym_DUSP6 + ym_DUSP6, 2)/std::pow(sigma_ym_DUSP6, 3);
            break;
        case 865:
            dJydsigma[865] = 1.0/sigma_ym_RPS6KA1 - 1.0*std::pow(-mym_RPS6KA1 + ym_RPS6KA1, 2)/std::pow(sigma_ym_RPS6KA1, 3);
            break;
        case 866:
            dJydsigma[866] = 1.0/sigma_ym_RPS6KA2 - 1.0*std::pow(-mym_RPS6KA2 + ym_RPS6KA2, 2)/std::pow(sigma_ym_RPS6KA2, 3);
            break;
        case 867:
            dJydsigma[867] = 1.0/sigma_ym_RPS6KA3 - 1.0*std::pow(-mym_RPS6KA3 + ym_RPS6KA3, 2)/std::pow(sigma_ym_RPS6KA3, 3);
            break;
        case 868:
            dJydsigma[868] = 1.0/sigma_ym_RPS6KA4 - 1.0*std::pow(-mym_RPS6KA4 + ym_RPS6KA4, 2)/std::pow(sigma_ym_RPS6KA4, 3);
            break;
        case 869:
            dJydsigma[869] = 1.0/sigma_ym_DUSP1 - 1.0*std::pow(-mym_DUSP1 + ym_DUSP1, 2)/std::pow(sigma_ym_DUSP1, 3);
            break;
        case 870:
            dJydsigma[870] = 1.0/sigma_ym_FOS - 1.0*std::pow(-mym_FOS + ym_FOS, 2)/std::pow(sigma_ym_FOS, 3);
            break;
        case 871:
            dJydsigma[871] = 1.0/sigma_ym_JUN - 1.0*std::pow(-mym_JUN + ym_JUN, 2)/std::pow(sigma_ym_JUN, 3);
            break;
        case 872:
            dJydsigma[872] = 1.0/sigma_ym_MYC - 1.0*std::pow(-mym_MYC + ym_MYC, 2)/std::pow(sigma_ym_MYC, 3);
            break;
        case 873:
            dJydsigma[873] = 1.0/sigma_ym_CTNNB1 - 1.0*std::pow(-mym_CTNNB1 + ym_CTNNB1, 2)/std::pow(sigma_ym_CTNNB1, 3);
            break;
        case 874:
            dJydsigma[874] = 1.0/sigma_ym_PTEN - 1.0*std::pow(-mym_PTEN + ym_PTEN, 2)/std::pow(sigma_ym_PTEN, 3);
            break;
        case 875:
            dJydsigma[875] = 1.0/sigma_ym_AKT1 - 1.0*std::pow(-mym_AKT1 + ym_AKT1, 2)/std::pow(sigma_ym_AKT1, 3);
            break;
        case 876:
            dJydsigma[876] = 1.0/sigma_ym_AKT2 - 1.0*std::pow(-mym_AKT2 + ym_AKT2, 2)/std::pow(sigma_ym_AKT2, 3);
            break;
        case 877:
            dJydsigma[877] = 1.0/sigma_ym_PDPK1 - 1.0*std::pow(-mym_PDPK1 + ym_PDPK1, 2)/std::pow(sigma_ym_PDPK1, 3);
            break;
        case 878:
            dJydsigma[878] = 1.0/sigma_ym_RICTOR - 1.0*std::pow(-mym_RICTOR + ym_RICTOR, 2)/std::pow(sigma_ym_RICTOR, 3);
            break;
        case 879:
            dJydsigma[879] = 1.0/sigma_ym_MTOR - 1.0*std::pow(-mym_MTOR + ym_MTOR, 2)/std::pow(sigma_ym_MTOR, 3);
            break;
        case 880:
            dJydsigma[880] = 1.0/sigma_ym_GSK3B - 1.0*std::pow(-mym_GSK3B + ym_GSK3B, 2)/std::pow(sigma_ym_GSK3B, 3);
            break;
        case 881:
            dJydsigma[881] = 1.0/sigma_ym_TSC1 - 1.0*std::pow(-mym_TSC1 + ym_TSC1, 2)/std::pow(sigma_ym_TSC1, 3);
            break;
        case 882:
            dJydsigma[882] = 1.0/sigma_ym_TSC2 - 1.0*std::pow(-mym_TSC2 + ym_TSC2, 2)/std::pow(sigma_ym_TSC2, 3);
            break;
        case 883:
            dJydsigma[883] = 1.0/sigma_ym_PRKCA - 1.0*std::pow(-mym_PRKCA + ym_PRKCA, 2)/std::pow(sigma_ym_PRKCA, 3);
            break;
        case 884:
            dJydsigma[884] = 1.0/sigma_ym_PRKCB - 1.0*std::pow(-mym_PRKCB + ym_PRKCB, 2)/std::pow(sigma_ym_PRKCB, 3);
            break;
        case 885:
            dJydsigma[885] = 1.0/sigma_ym_PRKCG - 1.0*std::pow(-mym_PRKCG + ym_PRKCG, 2)/std::pow(sigma_ym_PRKCG, 3);
            break;
        case 886:
            dJydsigma[886] = 1.0/sigma_ym_PRKCD - 1.0*std::pow(-mym_PRKCD + ym_PRKCD, 2)/std::pow(sigma_ym_PRKCD, 3);
            break;
        case 887:
            dJydsigma[887] = 1.0/sigma_ym_PEBP1 - 1.0*std::pow(-mym_PEBP1 + ym_PEBP1, 2)/std::pow(sigma_ym_PEBP1, 3);
            break;
        case 888:
            dJydsigma[888] = 1.0/sigma_ym_MAPK1 - 1.0*std::pow(-mym_MAPK1 + ym_MAPK1, 2)/std::pow(sigma_ym_MAPK1, 3);
            break;
        case 889:
            dJydsigma[889] = 1.0/sigma_ym_MAPK3 - 1.0*std::pow(-mym_MAPK3 + ym_MAPK3, 2)/std::pow(sigma_ym_MAPK3, 3);
            break;
        case 890:
            dJydsigma[890] = 1.0/sigma_ym_FOXO3 - 1.0*std::pow(-mym_FOXO3 + ym_FOXO3, 2)/std::pow(sigma_ym_FOXO3, 3);
            break;
        case 891:
            dJydsigma[891] = 1.0/sigma_ym_RHEB - 1.0*std::pow(-mym_RHEB + ym_RHEB, 2)/std::pow(sigma_ym_RHEB, 3);
            break;
        case 892:
            dJydsigma[892] = 1.0/sigma_ym_RPTOR - 1.0*std::pow(-mym_RPTOR + ym_RPTOR, 2)/std::pow(sigma_ym_RPTOR, 3);
            break;
        case 893:
            dJydsigma[893] = 1.0/sigma_ym_RPS6KB1 - 1.0*std::pow(-mym_RPS6KB1 + ym_RPS6KB1, 2)/std::pow(sigma_ym_RPS6KB1, 3);
            break;
        case 894:
            dJydsigma[894] = 1.0/sigma_ym_RPS6KB2 - 1.0*std::pow(-mym_RPS6KB2 + ym_RPS6KB2, 2)/std::pow(sigma_ym_RPS6KB2, 3);
            break;
        case 895:
            dJydsigma[895] = 1.0/sigma_ym_EIF4EBP1 - 1.0*std::pow(-mym_EIF4EBP1 + ym_EIF4EBP1, 2)/std::pow(sigma_ym_EIF4EBP1, 3);
            break;
        case 896:
            dJydsigma[896] = 1.0/sigma_ym_SOS1 - 1.0*std::pow(-mym_SOS1 + ym_SOS1, 2)/std::pow(sigma_ym_SOS1, 3);
            break;
        case 897:
            dJydsigma[897] = 1.0/sigma_ym_CDKN2A - 1.0*std::pow(-mym_CDKN2A + ym_CDKN2A, 2)/std::pow(sigma_ym_CDKN2A, 3);
            break;
        case 898:
            dJydsigma[898] = 1.0/sigma_ym_MDM4 - 1.0*std::pow(-mym_MDM4 + ym_MDM4, 2)/std::pow(sigma_ym_MDM4, 3);
            break;
        case 899:
            dJydsigma[899] = 1.0/sigma_ym_FGFR1 - 1.0*std::pow(-mym_FGFR1 + ym_FGFR1, 2)/std::pow(sigma_ym_FGFR1, 3);
            break;
        case 900:
            dJydsigma[900] = 1.0/sigma_ym_FGFR2 - 1.0*std::pow(-mym_FGFR2 + ym_FGFR2, 2)/std::pow(sigma_ym_FGFR2, 3);
            break;
        case 901:
            dJydsigma[901] = 1.0/sigma_ym_FGF1 - 1.0*std::pow(-mym_FGF1 + ym_FGF1, 2)/std::pow(sigma_ym_FGF1, 3);
            break;
        case 902:
            dJydsigma[902] = 1.0/sigma_ym_FGF2 - 1.0*std::pow(-mym_FGF2 + ym_FGF2, 2)/std::pow(sigma_ym_FGF2, 3);
            break;
        case 903:
            dJydsigma[903] = 1.0/sigma_ym_EIF4E - 1.0*std::pow(-mym_EIF4E + ym_EIF4E, 2)/std::pow(sigma_ym_EIF4E, 3);
            break;
        case 904:
            dJydsigma[904] = 1.0/sigma_ym_IRS1 - 1.0*std::pow(-mym_IRS1 + ym_IRS1, 2)/std::pow(sigma_ym_IRS1, 3);
            break;
        case 905:
            dJydsigma[905] = 1.0/sigma_ym_IRS2 - 1.0*std::pow(-mym_IRS2 + ym_IRS2, 2)/std::pow(sigma_ym_IRS2, 3);
            break;
        case 906:
            dJydsigma[906] = 1.0/sigma_ym_IGF1 - 1.0*std::pow(-mym_IGF1 + ym_IGF1, 2)/std::pow(sigma_ym_IGF1, 3);
            break;
        case 907:
            dJydsigma[907] = 1.0/sigma_ym_IGF2 - 1.0*std::pow(-mym_IGF2 + ym_IGF2, 2)/std::pow(sigma_ym_IGF2, 3);
            break;
        case 908:
            dJydsigma[908] = 1.0/sigma_ym_IGF1R - 1.0*std::pow(-mym_IGF1R + ym_IGF1R, 2)/std::pow(sigma_ym_IGF1R, 3);
            break;
        case 909:
            dJydsigma[909] = 1.0/sigma_ym_MSH6 - 1.0*std::pow(-mym_MSH6 + ym_MSH6, 2)/std::pow(sigma_ym_MSH6, 3);
            break;
        case 910:
            dJydsigma[910] = 1.0/sigma_ym_BRCA2 - 1.0*std::pow(-mym_BRCA2 + ym_BRCA2, 2)/std::pow(sigma_ym_BRCA2, 3);
            break;
        case 911:
            dJydsigma[911] = 1.0/sigma_ym_MGMT - 1.0*std::pow(-mym_MGMT + ym_MGMT, 2)/std::pow(sigma_ym_MGMT, 3);
            break;
        case 912:
            dJydsigma[912] = 1.0/sigma_ym_INSR - 1.0*std::pow(-mym_INSR + ym_INSR, 2)/std::pow(sigma_ym_INSR, 3);
            break;
        case 913:
            dJydsigma[913] = 1.0/sigma_ym_INS - 1.0*std::pow(-mym_INS + ym_INS, 2)/std::pow(sigma_ym_INS, 3);
            break;
        case 914:
            dJydsigma[914] = 1.0/sigma_yCytoplasm - 1.0*std::pow(-myCytoplasm + yCytoplasm, 2)/std::pow(sigma_yCytoplasm, 3);
            break;
        case 915:
            dJydsigma[915] = 1.0/sigma_yExtracellular - 1.0*std::pow(-myExtracellular + yExtracellular, 2)/std::pow(sigma_yExtracellular, 3);
            break;
        case 916:
            dJydsigma[916] = 1.0/sigma_yNucleus - 1.0*std::pow(-myNucleus + yNucleus, 2)/std::pow(sigma_yNucleus, 3);
            break;
        case 917:
            dJydsigma[917] = 1.0/sigma_yMitochondrion - 1.0*std::pow(-myMitochondrion + yMitochondrion, 2)/std::pow(sigma_yMitochondrion, 3);
            break;
    }
}

} // namespace model_SPARCED
} // namespace amici
