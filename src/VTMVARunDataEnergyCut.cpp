/*! \class VTMVARunDataEnergyCut
    \brief  VTMVARunDataEnergyCut data class for TMVA energy cuts


*/

#include "VTMVARunDataEnergyCut.h"

VTMVARunDataEnergyCut::VTMVARunDataEnergyCut()
{
    fEnergyCutBin = 0;
    fEnergyCut_Log10TeV_min = 0.;
    fEnergyCut_Log10TeV_max = 0;
    fEnergyCut = "";
    fEnergyReconstructionMethod = 1;
}

void VTMVARunDataEnergyCut::print()
{
    cout << "energy bin " << fEnergyCutBin;
    cout << ": log10(TeV) [" << fEnergyCut_Log10TeV_min;
    cout << ", " << fEnergyCut_Log10TeV_max << "]";
    cout << " method " << fEnergyReconstructionMethod << endl;
    cout << "\t cuts: " << fEnergyCut.GetTitle() << endl;
}


