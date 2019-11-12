//! helper functions for light curve analysis

#ifndef VFluxAndLightCurveUtilities_H
#define VFluxAndLightCurveUtilities_H

#include "TMath.h"

using namespace std;

namespace VFluxAndLightCurveUtilities
{

    double convertPhotonFlux_to_CrabUnits( double iEnergy_TeV, double iFlux, double iSpectralIndex = 2.49, bool bLinerarEnergyScale = true );
    double convertPhotonFlux_to_Ergs( double iEnergy_TeV, double iFlux, bool bLin = true );
    double convertEnergy_keV_to_Hz( double iEnergy_keV );
    double convertEnergy_TeV_to_Hz( double iEnergy_TeV );
    
    double getOrbitalPhase( double iMJD, double iZeroPhase_MJD, double iPhase_Days );
    double getOrbitalPhaseError( double iMJD, double iZeroPhase_MJD, double iOrbitalPeriod_Days, double iOrbitalPeriod_loE_days = 0., double iOrbitalPeriod_hiE_days = 0. );
}

#endif
