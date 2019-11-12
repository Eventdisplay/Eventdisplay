/*! \class VTMVARunDataZenithCut
    \brief  VTMVARunDataZenithCut data class for TMVA energy cuts


*/

#include "VTMVARunDataZenithCut.h"

VTMVARunDataZenithCut::VTMVARunDataZenithCut()
{
    fZenithCutBin = 0;
    fZenithCut_min = 0.;
    fZenithCut_max = 0;
    fZenithCut = "";
}

double VTMVARunDataZenithCut::getSecant( double iZ )
{
    if( iZ > 89. )
    {
        return 1. / cos( 89. * TMath::DegToRad() );
    }
    return 1. / cos( iZ * TMath::DegToRad() );
}

void VTMVARunDataZenithCut::print()
{
    cout << "zenith bin " << fZenithCutBin;
    cout << ":  [" << fZenithCut_min;
    cout << ", " << fZenithCut_max << "] deg";
    cout << " (secant " << getSecant( fZenithCut_min );
    cout << ", " << getSecant( fZenithCut_max ) << ")";
    cout << "\t cuts: " << fZenithCut.GetTitle() << endl;
}


