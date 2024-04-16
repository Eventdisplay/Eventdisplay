/*! \file VMonteCarloRunHeader.cpp
    \brief MC parameter class

    use same parameters and parameter definition as sim_telarray

*/

#include "VMonteCarloRunHeader.h"


VMonteCarloRunHeader::VMonteCarloRunHeader()
{
    reset();
}

void VMonteCarloRunHeader::reset()
{
    runnumber = -1;
    shower_prog_id = 1;      ///< CORSIKA=1, ALTAI=2, KASCADE=3, MOCCA=4.
    shower_prog_vers = 0;    ///< version * 1000
    detector_prog_id = 0;    ///< sim_telarray=1, ...
    detector_prog_vers = 0;  ///< version * 1000
    converter_prog_vers = 0;  ///< version * 1000

    shower_date = 0;
    detector_date = 0;

    obsheight = 9.;        ///< Height of simulated observation level.
    num_showers = 0;         ///< Number of showers simulated.
    num_use = 0;             ///< Number of uses of each shower.
    primary_id = 0;
    core_pos_mode = 0;       ///< Core position fixed/circular/rectangular/...
    core_range[0] = 0.;    ///< rmin+rmax or dx+dy.
    core_range[1] = 0.;    ///< rmin+rmax or dx+dy.
    az_range[0] = 0.;      ///< Range of shower azimuth [rad, N->E].
    az_range[1] = 0.;      ///< Range of shower azimuth [rad, N->E].
    alt_range[0] = 0.;     ///< Range of shower altitude [rad].
    alt_range[1] = 0.;     ///< Range of shower altitude [rad].
    diffuse = 0;             ///< Diffuse mode off/on
    viewcone[0] = 0.;      ///< Min.+max. opening angle for diffuse mode [rad].
    viewcone[1] = 0.;      ///< Min.+max. opening angle for diffuse mode [rad].
    E_range[0] = 0.;       ///< Energy range [TeV] of simulated showers.
    E_range[1] = 0.;       ///< Energy range [TeV] of simulated showers.
    spectral_index = 0.;   ///< Power-law spectral index of spectrum (<0).
    B_total = 0.;          ///< Total geomagnetic field assumed [microT].
    B_inclination = 0.;    ///< Inclination of geomagnetic field [rad].
    B_declination = 0.;    ///< Declination of geomagnetic field [rad].
    injection_height = 0.; ///< Height of particle injection [m].
    fixed_int_depth = 0.;  ///< Fixed depth of first interaction or 0 [g/cm^2].
    atmosphere = 0;          ///< Atmospheric model number.
    corsika_iact_options = 0;
    corsika_low_E_model = 0;
    corsika_high_E_model = 0;
    corsika_low_high_E = 0.;
    corsika_bunchsize = 0.;
    corsika_wlen_min = 0.;
    corsika_wlen_max = 0.;
    corsika_low_E_detail = 0;
    corsika_high_E_detail = 0;
    detector_Simulator = "NOSET";

    fFADC_hilo_multipler = -999.;

    combined_runHeader = false;
}

bool VMonteCarloRunHeader::VOLUMEDET_set()
{
    bitset<32> EVTH76 = corsika_iact_options;
    for( unsigned int i = 0; i < 10; i++ )
    {
        EVTH76.set( i, 0 );
    }
    unsigned int iATM_tab = EVTH76.to_ulong() / 1024;
    EVTH76 = corsika_iact_options - 1024 * iATM_tab;
    return EVTH76.test( 5 );
}

void VMonteCarloRunHeader::printRunNumber()
{
    cout << runnumber << endl;
}

void VMonteCarloRunHeader::printMCAz( bool iLowerLimit )
{
    if( iLowerLimit )
    {
        cout <<  TMath::Nint( az_range[0] * 45. / atan( 1. ) ) << endl;
        return;
    }

    cout << TMath::Nint( az_range[1] * 45. / atan( 1. ) ) << endl;
}


void VMonteCarloRunHeader::print()
{
    cout << endl;
    cout << "Monte Carlo run header" << endl;
    cout << "======================" << endl;
    if( runnumber >= 0 )
    {
        cout << "run number: " << runnumber << endl;
    }
    cout << "code version: shower prog " << shower_prog_id << " (" << shower_prog_vers << "), ";
    cout << "detector prog " << detector_prog_id << " (" << detector_prog_vers << "), ";
    cout << "convert " << converter_prog_vers << endl;
    cout << "date: " << shower_date << "\t" << detector_date << endl;
    cout << "number of showers: " << num_showers << " (each shower used " << num_use << " times)" << endl;
    cout << "Primary " << primary_id << endl;
    cout << "Energy range: [" << E_range[0] << ", " << E_range[1] << "] TeV, powerlaw index " << spectral_index << endl;
    cout << "Core scattering: " << core_range[0] << "\t" << core_range[1] << " [m]";
    if( core_pos_mode == 0 )
    {
        cout << " (fixed)";
    }
    else if( core_pos_mode == 1 )
    {
        cout << " (circular)";
    }
    else if( core_pos_mode == 2 )
    {
        cout << " (rectangular)";
    }
    cout << endl;
    cout << "Azimuth range: [" << az_range[0] * 45. / atan( 1. ) << ", " << az_range[1] * 45. / atan( 1. ) << "]" << endl;
    cout << "Zenith range: [" << 90. - alt_range[0] * 45. / atan( 1. ) << ", " << 90. - alt_range[1] * 45. / atan( 1. ) << "]" << endl;
    cout << "Viewcone: [" << viewcone[0] << ", " << viewcone[1] << "] (" << diffuse << ")" << endl;
    cout << "Observatory height " << obsheight << " [m]" << endl;
    cout << "B-Field: " << B_total << " microT (" << B_inclination * 45. / atan( 1. );
    cout << "," << B_declination * 45. / atan( 1. ) << ")" << endl;
    cout << "Atmospheric model: " << atmosphere << endl;
    cout << "Cherenkov photon wavelength range: [" << corsika_wlen_min << ", " << corsika_wlen_max << "]" << endl;
    cout << "CORSIKA interaction detail: lowE " << corsika_low_E_detail << ", highE " << corsika_high_E_detail;
    cout << ", interaction models: lowE " << corsika_low_E_model << ", highE " << corsika_high_E_model;
    cout << ", transition energy " << corsika_low_high_E << " GeV" << endl;
    // print CHERENKOV FLAG
    cout << "CORSIKA iact options: " << corsika_iact_options << endl;
    bitset<32> EVTH76 = corsika_iact_options;
    for( unsigned int i = 0; i < 10; i++ )
    {
        EVTH76.set( i, 0 );
    }
    unsigned int iATM_tab = EVTH76.to_ulong() / 1024;
    EVTH76 = corsika_iact_options - 1024 * iATM_tab;
    cout << "CERENKOV " << EVTH76.test( 0 ) << endl;
    cout << "IACT " << EVTH76.test( 1 ) << endl;
    cout << "CEFFIC " << EVTH76.test( 2 ) << endl;
    cout << "ATMEXT " << EVTH76.test( 3 ) << endl;
    cout << "ATMEXT with refraction " << EVTH76.test( 4 ) << endl;
    cout << "VOLUMEDET " << EVTH76.test( 5 ) << endl;
    cout << "CURVED " << EVTH76.test( 6 ) << endl;
    cout << "SLANT " << EVTH76.test( 8 ) << endl;
    if( combined_runHeader )
    {
        cout << "(note that this is a MC run header filled from several simulation input files.";
        cout << " Some parameters might not be correct for a heterogen mix of simulation files)" << endl;
    }
    cout << endl << endl;
}

double VMonteCarloRunHeader::getMeanZenithAngle_Deg()
{
    return 0.5 * ( 90. - alt_range[0] * 45. / atan( 1. ) + 90. - alt_range[1] * 45. / atan( 1. ) );
}
