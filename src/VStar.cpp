/* \class VStar
   \brief data class for star catalogues

*/

#include "VStar.h"

VStar::VStar()
{
    fStarID = 0;
    fStarName = "";
    fDec2000 = 0.;
    fRA2000 = 0.;
    fDecCurrentEpoch = 0.;
    fRACurrentEpoch = 0.;
    fRunGalLong1958 = 0.;
    fRunGalLat1958 = 0.;
    fBrightness_V = 0.;
    fBrightness_B = 0.;
    fMajorDiameter = 0.;                        // this is either the source diameter or the positional error
    fMinorDiameter = 0.;
    fPositionAngle = 0.;
    fMajorDiameter_68 = 0.;
    fMinorDiameter_68 = 0.;
    fPositionAngle_68 = 0.;
    fSignificance = 0.;
    fSpectralIndex = 0.;
    fSpectralIndexError = 0.;
    fSpectrumType = "PowerLaw";
    fCutOff_MeV = -9999.;
    fCutOffError_MeV = -9999.;
    fVariability = false;
    fType = "";
    fQualityFlag = 0;
}

/*

   return angular separation of two stars in [deg]

*/
double VStar::getDistance( VStar* iStar1, VStar* iStar2 )
{
    if( iStar1 == 0 )
    {
        return -99.;
    }
    
    if( iStar2 == 0 )
    {
        return TMath::RadToDeg() * VAstronometry::vlaDsep( iStar1->fRA2000 * TMath::DegToRad(), iStar1->fDec2000 * TMath::DegToRad(),
                fRA2000 * TMath::DegToRad(), fDec2000 * TMath::DegToRad() );
    }
    
    return TMath::RadToDeg() * VAstronometry::vlaDsep( iStar1->fRA2000 * TMath::DegToRad(), iStar1->fDec2000 * TMath::DegToRad(),
            iStar2->fRA2000 * TMath::DegToRad(), iStar2->fDec2000 * TMath::DegToRad() );
}

void VStar::printStar()
{
    cout << fStarName << "\t" << fRA2000 << "\t" << fDec2000 << " B: " << fBrightness_B << " V: " << fBrightness_V << endl;
}
