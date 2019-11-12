/*! \class VOrbitalPhaseData
    \brief data element / calculation of orbital phases

*/

#include "VOrbitalPhaseData.h"

VOrbitalPhaseData::VOrbitalPhaseData()
{
    fName = "";
    fZeroPhase_MJD = 0.;
    fOrbit_days = 0.;
    fOrbit_days_error_low = 0.;
    fOrbit_days_error_high = 0.;
}

bool VOrbitalPhaseData::isSet()
{
    if( fOrbit_days > 0 && fZeroPhase_MJD >= 0. )
    {
        return true;
    }
    return false;
}

void VOrbitalPhaseData::print()
{
    if( fName.size() > 0 )
    {
        cout << fName << ": ";
    }
    cout << "MJD_0 = " << fZeroPhase_MJD << ", binary orbit [days]: ";
    cout << fOrbit_days << "+-[" << fOrbit_days_error_low << ", " << fOrbit_days_error_high << "]";
    cout << endl;
}

double VOrbitalPhaseData::getOrbitalPhase( double iMJD )
{
    return VFluxAndLightCurveUtilities::getOrbitalPhase( iMJD, fZeroPhase_MJD, fOrbit_days );
}
