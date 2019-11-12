/* \file  VFluxAndLightCurveUtilities.cpp
   \brief helper functions for flux and light curve analysis

*/

#include "VFluxAndLightCurveUtilities.h"



/*

   get flux in Crab Nebula units; using Whipple 1998 Crab spectrum

   (note: spectral index should be adjusted)

*/
double VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( double iEnergy_TeV, double iFlux, double iSpectralIndex, bool bLin )
{
    if( !bLin )
    {
        iEnergy_TeV = TMath::Power( 10., iEnergy_TeV );
    }
    
    // Whipple (Hillas 1998)
    if( iSpectralIndex != -1. )
    {
        double iCrab = 3.2e-7 / ( iSpectralIndex - 1. ) * TMath::Power( iEnergy_TeV, -iSpectralIndex + 1. ) / 1.e4;
        return iFlux / iCrab;
    }
    
    return 0.;
}

/*
    convert fluxes in 1/cm2/s to ergs/cm2/s

    e    :    energy in TeV
    f    :    flux per energy bin in 1/cm2/s
*/
double VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( double iEnergy_TeV, double iFlux, bool bLin )
{
    if( bLin )
    {
        iEnergy_TeV *= 1.e12;
    }
    else
    {
        iEnergy_TeV = TMath::Power( 10., iEnergy_TeV ) * 1.e12;
    }
    // eV / cm2 / s
    iFlux *= iEnergy_TeV * iEnergy_TeV / 1.e12;
    
    // eV -> J
    iFlux *= TMath::Qe();
    // J -> ergs
    iFlux /= 1.e-7;
    
    return iFlux;
}


/*
   convert energies from TeV to Hz

   linear energy
*/
double VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( double iEnergy_TeV )
{
    double constant_TeVtoHz = TMath::C() * 1.e12 /  1.239841875e-6;
    return iEnergy_TeV * constant_TeVtoHz;
}

/*
   convert energies from keV to Hz

   linear energy
*/
double VFluxAndLightCurveUtilities::convertEnergy_keV_to_Hz( double iEnergy_keV )
{
    return convertEnergy_TeV_to_Hz( iEnergy_keV * 1.e-9 );
}


/*

   calculate orbital phase for binaries

   all input values in MJDs and days

*/
double VFluxAndLightCurveUtilities::getOrbitalPhase( double iMJD, double iZeroPhase_MJD, double iOrbitalPeriod_Days )
{
    if( iOrbitalPeriod_Days > 0. )
    {
        iMJD = ( iMJD - iZeroPhase_MJD ) / iOrbitalPeriod_Days;
        iMJD =   iMJD - TMath::Floor( iMJD );
        return iMJD;
    }
    return -99.;
}


/*

   propagate error on orbital period into phase error

*/
double VFluxAndLightCurveUtilities::getOrbitalPhaseError( double iMJD, double iZeroPhase_MJD, double iOrbitalPeriod_Days, double iOrbitalPeriod_loE_days, double iOrbitalPeriod_hiE_days )
{
    // calculate average error from upper and lower value
    double iError = sqrt( iOrbitalPeriod_loE_days * iOrbitalPeriod_loE_days + iOrbitalPeriod_hiE_days * iOrbitalPeriod_hiE_days );
    
    double iP = 0.;
    if( iOrbitalPeriod_Days > 0. )
    {
        iP  = ( iMJD - iZeroPhase_MJD ) * ( iMJD - iZeroPhase_MJD ) / iOrbitalPeriod_Days / iOrbitalPeriod_Days / iOrbitalPeriod_Days / iOrbitalPeriod_Days;
        iP *= iError * iError;
        return sqrt( iP );
    }
    
    return -99.;
}
