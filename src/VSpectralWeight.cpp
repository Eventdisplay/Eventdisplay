/*! \file   getEnergyWeighting
    \brief calculate energy dependent spectral weight

    general assumption: MC is simulated according to a power law between min and max energy

    use power laws only (but can easily be expanded for arbitrary spectral shapes)

*/

#include "VSpectralWeight.h"

VSpectralWeight::VSpectralWeight()
{
    fDebug = false;
    
    fIndex = 2.;
    fSpectralWeightAlpha = 1.;
    
    fMCSpectralIndex = 2.;
    fMCMinEnergy_TeV_Lin = 0.03;
    fMCMaxEnergy_TeV_Lin = 200.;
    fMCMaxConeAngle_deg = 0.;
    fMCSimulatedEvents = 0.;
    fMCFluxConstant = 1.;
    fEnergySpectrum = 0;
}

/*

    spectral index iMCSpectralIndex > 0 ( E^{-iMCSpectralIndex} )

*/
void VSpectralWeight::setMCParameter( double iMCSpectralIndex, double iMCEnergy_min_TeV_Lin, double iMCEnergy_max_TeV_Lin,
                                      double iMCMaxConeAngle_deg, double iMCSimulatedEvents )
{
    fMCSpectralIndex = iMCSpectralIndex;
    fMCMinEnergy_TeV_Lin = iMCEnergy_min_TeV_Lin;
    fMCMaxEnergy_TeV_Lin = iMCEnergy_max_TeV_Lin;
    fMCMaxConeAngle_deg = iMCMaxConeAngle_deg;
    fMCSimulatedEvents = iMCSimulatedEvents;
    
    calculateMCFluxConstant();
    
    setSpectralIndex( fIndex );
}

void VSpectralWeight::calculateMCFluxConstant()
{
    if( fMCSimulatedEvents < 0. || fMCMaxConeAngle_deg < 0. )
    {
        fMCFluxConstant = 1.;
        return;
    }
    // solid angle
    double iS = 2. * TMath::Pi() * ( 1. - cos( fMCMaxConeAngle_deg * TMath::DegToRad() ) );
    
    double iP = TMath::Power( fMCMaxEnergy_TeV_Lin, -1.*fMCSpectralIndex + 1 ) - TMath::Power( fMCMinEnergy_TeV_Lin, -1.*fMCSpectralIndex + 1 );
    if( iS > 0. && iP > 0. )
    {
        fMCFluxConstant = fMCSimulatedEvents * ( -1.*fMCSpectralIndex + 1 ) / iS / iP;
    }
    else
    {
        fMCFluxConstant = 1.;
    }
}

/*
    weights will be always > 1

    spectral index iG > 0 ( E^{-iG} )
*/
void VSpectralWeight::setSpectralIndex( double iG, bool iPrint )
{
    fIndex = iG;
    
    if( iPrint )
    {
        cout << "weighting events to spectral index of " << fIndex << endl;
    }
    
    if( fabs( fIndex - fMCSpectralIndex ) < 0.02 )
    {
        fSpectralWeightAlpha = 1.;
    }
    else if( fIndex > fMCSpectralIndex )
    {
        fSpectralWeightAlpha = TMath::Power( fMCMaxEnergy_TeV_Lin, -1.*fMCSpectralIndex ) / TMath::Power( fMCMaxEnergy_TeV_Lin, -1.*fIndex );
    }
    else if( fIndex < fMCSpectralIndex )
    {
        fSpectralWeightAlpha = TMath::Power( fMCMinEnergy_TeV_Lin, -1.*fMCSpectralIndex ) / TMath::Power( fMCMinEnergy_TeV_Lin, -1.*fIndex );
    }
}

/*

    energy in TeV (linear scale)

*/
double VSpectralWeight::getSpectralWeight( double iE_TeV_lin )
{
    // weighting for IRFs
    if( fEnergySpectrum == 0 )
    {
        if( fabs( fIndex - fMCSpectralIndex ) < 0.01 )
        {
            return 1.;
        }
        else
        {
            return fSpectralWeightAlpha * TMath::Power( iE_TeV_lin, -1.*fIndex + fMCSpectralIndex );
        }
    }
    // exact weighting for CR spectra (arbitrary function)
    else
    {
        double d = fMCFluxConstant * TMath::Power( iE_TeV_lin, -1.*fMCSpectralIndex );
        double w = 1.;
        if( d > 0. )
        {
            w = fEnergySpectrum->Eval( iE_TeV_lin );  // is energy input log or lin??
            return w / d;
        }
        return w;
    }
    return 1.;
}

void VSpectralWeight::print()
{
    cout << "VSpectralWeight: ";
    cout << "expect input (MC) energy spectra with index " << fMCSpectralIndex;
    cout << " in energy range [" << fMCMinEnergy_TeV_Lin << ", " << fMCMaxEnergy_TeV_Lin << "] TeV";
    cout << endl;
    if( fMCMaxConeAngle_deg > 0. && fMCSimulatedEvents > 0. )
    {
        cout << "\t MC cone angle: " << fMCMaxConeAngle_deg << "[deg]";
        cout << ", # of MC events: " << fMCSimulatedEvents;
        cout << "\t MC flux constant " << fMCFluxConstant << endl;
    }
}
