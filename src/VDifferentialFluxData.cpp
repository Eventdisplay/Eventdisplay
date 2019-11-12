/*! \class VDifferentialFluxData
    \brief data and converter class for differential flux values


*/

#include "VDifferentialFluxData.h"

VDifferentialFluxData::VDifferentialFluxData()
{
    MJD_min = 0.;
    MJD_max = 0.;
    Energy = 0.;
    Energy_lowEdge = 0.;
    Energy_upEdge = 0.;
    Energy_lowEdge_bin = 0;
    Energy_upEdge_bin = 0;
    Energy_Hz = 0.;
    EnergyWeightedMean = 0.;
    dE = 0.;
    
    resetCountsAndFluxes();
}


void VDifferentialFluxData::resetCountsAndFluxes()
{
    EffectiveArea = 0.;
    DifferentialFlux = 0.;
    DifferentialFluxError = 0.;
    DifferentialFluxError_low = 0.;
    DifferentialFluxError_up = 0.;
    DifferentialFlux_vFv = 0.;
    DifferentialFluxError_vFv = 0.;
    DifferentialFluxError_low_vFv = 0.;
    DifferentialFluxError_up_vFv = 0.;
    DifferentialFluxWeighting = 0.;
    ObsTime = 0.;
    ExposureTime = 0.;
    NOn = 0.;
    NOn_error = 0.;
    NOff = 0.;
    NOff_error = 0.;
    NOff_alpha = 1.;
    Significance = 0.;
}


/*
    print all data

    bSED = true: print vF_v fluxes to be used in VSpectralEnergyDistribution
*/
void VDifferentialFluxData::print( bool bSED )
{
    ios::fmtflags f( cout.flags() );
    
    if( !bSED )
    {
        cout << "E: " << setprecision( 2 ) << setw( 4 ) << Energy << " [TeV]";
        cout << " (dE = (" << Energy_lowEdge << "-" << Energy_upEdge << ")TeV";
        cout << " = " << dE << " TeV)";
        if( DifferentialFluxError > 0. )
        {
            cout << scientific << setprecision( 2 ) <<  "\tdiff flux: " << DifferentialFlux << " +- " << DifferentialFluxError;
        }
        else
        {
            cout << scientific << setprecision( 2 ) << "\tUL: " << DifferentialFlux;
        }
        cout << " [1/cm^2/s/TeV]" << endl;
        cout << setw( 7 ) << fixed << setprecision( 1 ) << "NOn: " << NOn << "+-" << NOn_error;
        cout << setw( 7 ) << fixed << setprecision( 1 ) << "\tNOff: " << NOff << "+-" << NOff_error;
        cout << setprecision( 2 ) << " (alpha=" << NOff_alpha << ")";
        cout << setw( 7 ) << fixed << setprecision( 1 ) << "\tSign.: " << Significance << " sigma";
        cout << setw( 7 ) << fixed << setprecision( 1 ) << "\tObs.Time: " << ObsTime << "[s]";
        cout << endl;
    }
    else
    {
        if( Energy_Hz > 0. )
        {
            cout << scientific;
            cout << setprecision( 8 ) << MJD_min << "\t" << MJD_max << "\t";
            cout << setprecision( 3 ) << Energy_Hz << "\t";
            cout << DifferentialFlux_vFv << "\t";
            cout << DifferentialFluxError_vFv;
            cout << endl;
        }
    }
    
    cout.flags( f );
}


/*
    print all data but without texts, so that the output can be copied
    for printing externally. e.g. in gnuplot

    bSED = true: print vF_v fluxes to be used in VSpectralEnergyDistribution
*/
void VDifferentialFluxData::printClean( bool bSED, bool b_PrintError_low_up )
{
    if( !bSED )
    {
        cout <<  setprecision( 3 ) << fixed << setw( 9 ) << Energy;
        // 		cout << "    " << Energy_lowEdge << "    " << Energy_upEdge ;
        // 		cout << "    " << dE ;
        cout << setw( 9 )  << fixed << Energy_lowEdge << setw( 8 )  << fixed << Energy_upEdge ;
        cout << setw( 9 )  << fixed << dE ;
        if( DifferentialFluxError > 0. && b_PrintError_low_up )
        {
            cout << scientific << setprecision( 2 ) <<  "      " << DifferentialFlux;
            cout << "  [" << DifferentialFluxError_low << ", " << DifferentialFluxError_up << "] ";
        }
        else if( DifferentialFluxError > 0. )
        {
            cout << scientific << setprecision( 2 ) <<  "      " << DifferentialFlux;
            cout << " [" << DifferentialFluxError << ", " << DifferentialFluxError << "] ";
        }
        else
        {
            cout << scientific << setprecision( 2 ) << "      " << DifferentialFlux;
            cout << "                     ";
        }
        cout << setw( 6 ) << right << fixed << setprecision( 0 ) << NOn;
        cout << setw( 6 ) << right << fixed << setprecision( 1 ) << NOn_error;
        cout << setw( 8 ) << right << fixed << setprecision( 1 ) << NOff << " ";
        cout << setw( 6 ) << right << fixed << setprecision( 1 ) << NOff_error;
        cout << setw( 6 ) << right << fixed << setprecision( 2 ) << NOff_alpha;
        cout << setw( 6 ) << right << fixed << setprecision( 1 ) << Significance;
        cout << setw( 12 ) << right << fixed << setprecision( 1 ) << ObsTime;
        if( EffectiveArea > 0. )
        {
            cout << setw( 9 ) << right << fixed << EffectiveArea;
        }
        cout << endl;
    }
    else
    {
        if( Energy_Hz > 0. )
        {
            cout << scientific;
            cout << setprecision( 8 ) << MJD_min << "\t" << MJD_max << "\t";
            cout << setprecision( 3 ) << Energy_Hz << "\t";
            cout << DifferentialFlux_vFv << "\t";
            cout << DifferentialFluxError_vFv;
            cout << endl;
        }
    }
}


/*
    calculate energy in Hz, vF_v, etc...
*/
void VDifferentialFluxData::fillEvent( double iMJD_min, double iMJD_max )
{
    MJD_min = iMJD_min;
    MJD_max = iMJD_max;
    
    Energy_Hz = VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( Energy );
    
    DifferentialFlux_vFv      = VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( Energy, DifferentialFlux );
    DifferentialFluxError_vFv = VFluxAndLightCurveUtilities::convertPhotonFlux_to_Ergs( Energy, DifferentialFluxError );
}


/*

   calculate nuFnu from an integral flux point

   energy in eV

   assume power laws dN/dE = c * E^{gamma}

   return value in [erg/cm2/s]

*/
double VDifferentialFluxData::nuFnu( double F, double gamma, double e1, double e2, double e3 )
{
    if( e1 >= e2 )
    {
        cout << "invalid energy interval " << e1 << ", " << e2 << endl;
        return -99.;
    }
    
    // calculate constant
    double c = 0.;
    if( gamma != -1. )
    {
        c = F * ( gamma + 1. ) * TMath::Power( e3, gamma ) / ( TMath::Power( e2, gamma + 1. ) - TMath::Power( e1, gamma + 1. ) );
    }
    // (not correct)
    else
    {
        c = F * ( log( e2 ) - log( e1 ) );
    }
    
    // calculate nuFu
    double nF = 0.;
    if( e3 > 0. )
    {
        nF = c * TMath::Power( e3 / e3, gamma ) * e3 * e3;
    }
    // Following A.Tramacere (Fermi Saas Fee analysis session; 2010)
    // http://www.isdc.unige.ch/sf2010/fermi
    else
    {
        nF = c * TMath::Power( sqrt( e1 * e2 ) / e3, gamma + 2. );
    }
    
    // from eV to ergs
    nF *= TMath::Qe();
    nF /= 1.e-7;
    
    return nF;
}
