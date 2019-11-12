/*! \class VZDCFData
    \brief  data class for ZDCF plotting


*/

#include "VZDCFData.h"

VZDCFData::VZDCFData()
{
    tau = 0.;
    sigma_tau_neg = 0.;
    sigma_tau_pos = 0.;
    dcf = 0.;
    dcf_error_low = 0.;
    dcf_error_up = 0.;
    nbins = 0;
}

void VZDCFData::print()
{
    cout << "Tau: " << tau;
    cout << "\t -sig(tau) " << sigma_tau_neg;
    cout << "\t +sig(tau) " << sigma_tau_pos;
    cout << "\t dcf " << dcf;
    cout << "\t -err(dcf) " << dcf_error_low;
    cout << "\t +err(dcf) " << dcf_error_up;
    cout << "\t nbins " << nbins;
    cout << endl;
}
