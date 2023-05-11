/*! \class VSpectralEnergyDistribution
    \brief calculate spectral energy distribution

    (example code)

    // create new instance with name "WComae"
    VSpectralEnergyDistribution *a = new VSpectralEnergyDistribution( "WComae" );

    // read photometric data
    a->readPhotoMetricBands("multiwavelengthdata/photometricBands-JMC.dat" );

    // read some optical data from file 2008_JuneFlare/Tuorla_R.dat
    a->readOpticalData("Tuorla (R)", "2008_JuneFlare/Tuorla_R.dat", "R", true, 20, 2 );

    // read some Swift data from file 2008_JuneFlare/SwiftUVOT_SED_2008-06-08.txt for MJD 54625.25
    a->readSwiftData( "Swift XRT", "2008_JuneFlare/SwiftXRT_SED_2008-06-08.txt", 54625.25, 54625.25, true, 21, 3 );

    // read some VERITAS data from file 2008_JuneFlare/VERITAS_Flare2_d20080607.txt
    a->readTeVEvndispData( "VERITAS", "2008_JuneFlare/VERITAS_Flare2_d20080607.txt", true, 22, 4 );

    // plot SED
    a->plot();

    (end example code)

*/

#include "VSpectralEnergyDistribution.h"

/*
     constructor

     iname   :   name used for naming of histograms, canvases, etc.

*/
VSpectralEnergyDistribution::VSpectralEnergyDistribution( string iname )
{
    fName = iname;
    
    fDebug = false;
    
    setTimeRange();
    setPlottingEnergyRange_Hz();
    setPlottingFluxRange();
}


/*
    set time range for plotting

    (values outside of this range are still read into the data vector (fSpectralFlux))
*/
void VSpectralEnergyDistribution::setTimeRange( double iMJDmin, double iMJDmax )
{
    fMJDMin = iMJDmin;
    fMJDMax = iMJDmax;
}


/*
  read TeV data from an ascii file

  use VEnergySpectrum::printDifferentialFluxes( true ); to get this data from your anasum output

*/
bool VSpectralEnergyDistribution::readTeVEvndispData( string name, string ifile, bool bPrint, int imarker, int icolor )
{
    return readDataFile( name, ifile, -1., -1., bPrint, imarker, icolor );
}


/*
    read FermiLAT data from an ascii file

    iname   :   data descriptor
    ifile   :   file name (ascii format)
    MJD_min :   MJD of begin of measurements
    MJD_max:    MJD of end of measurements
    bPrint  :   print values to screen
    imarker :   marker style for plotting
    icolor  :   marker and line color for plotting

all results are added to the vector of photon fluxes (sPhotonFlux)

expected file format (0):
<frequency (GeV)> <flux (erg cm^-2 s^-1)> <error in Flux (erg cm^-2 s^-1)>

alternative file format (1):
<frequency (GeV)> <flux (erg cm^-2 s^-1)> <exl> <exh> <eyl> <eyh>

An upper limit is drawn for fluxes with errors <= 0.

*/
bool VSpectralEnergyDistribution::readFermiData( string name, string ifile, double MJD_min, double MJD_max, bool bPrint, int imarker, int icolor, int iFormat )
{
    sPhotonFlux i_pF_temp;
    vector< sPhotonFlux > i_photonFlux;
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading data from " << ifile << endl;
        return false;
    }
    string is_line;
    string is_temp;
    int z = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        if( is_line.substr( 0, 1 ) == "*" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        i_pF_temp.name = name;
        if( MJD_min > 0. && MJD_max > 0. )
        {
            i_pF_temp.MJD_min = MJD_min;
            i_pF_temp.MJD_max = MJD_max;
        }
        else
        {
            is_stream >> is_temp;
            double iT_min = atof( is_temp.c_str() );
            is_stream >> is_temp;
            double iT_max = atof( is_temp.c_str() );
            
            // check if following line is at a different time period (if yes, clear the vector and start again)
            if( z > 0 && ( TMath::Abs( i_photonFlux.back().MJD_min - iT_min ) > 1.e-2 || ( TMath::Abs( i_photonFlux.back().MJD_max - iT_max ) > 1.e-2 ) ) )
            {
                i_pF_temp.energy_Hz.clear();
                i_pF_temp.energy_Hz_min.clear();
                i_pF_temp.energy_Hz_max.clear();
                i_pF_temp.energy_eV.clear();
                i_pF_temp.energy_eV_min.clear();
                i_pF_temp.energy_eV_max.clear();
                i_pF_temp.flux_ergscms.clear();
                i_pF_temp.flux_error_up_ergscms.clear();
                i_pF_temp.flux_error_down_ergscms.clear();
            }
            i_pF_temp.MJD_min = iT_min;
            i_pF_temp.MJD_max = iT_max;
        }
        i_pF_temp.Marker = imarker;
        i_pF_temp.Color = icolor;
        
        // energy [GeV]
        double i_eV = 0.;
        is_stream >> is_temp;
        i_eV = atof( is_temp.c_str() ) * 1.e9;
        if( i_eV  <= 0. )
        {
            continue;
        }
        i_pF_temp.energy_eV.push_back( i_eV );
        double iHz = i_eV / 1.239841875e-6 * TMath::C();
        i_pF_temp.energy_Hz.push_back( iHz );
        
        if( iFormat == 0 )
        {
            i_pF_temp.energy_eV_min.push_back( i_pF_temp.energy_eV.back() );
            i_pF_temp.energy_eV_max.push_back( i_pF_temp.energy_eV.back() );
            i_pF_temp.energy_Hz_min.push_back( i_pF_temp.energy_Hz.back() );
            i_pF_temp.energy_Hz_max.push_back( i_pF_temp.energy_Hz.back() );
        }
        is_stream >> is_temp;
        i_pF_temp.flux_ergscms.push_back( atof( is_temp.c_str() ) );
        
        if( iFormat == 0 )
        {
            is_stream >> is_temp;
            i_pF_temp.flux_error_up_ergscms.push_back( atof( is_temp.c_str() ) );
            i_pF_temp.flux_error_down_ergscms.push_back( atof( is_temp.c_str() ) );
        }
        else if( iFormat == 1 )
        {
            double i_eV_low = 0.;
            is_stream >> is_temp;
            i_eV_low = atof( is_temp.c_str() ) * 1.e9;
            
            i_pF_temp.energy_eV_min.push_back( i_eV - i_eV_low );
            i_pF_temp.energy_Hz_min.push_back( ( i_eV - i_eV_low ) / 1.239841875e-6 * TMath::C() );
            
            double i_eV_high = 0.;
            is_stream >> is_temp;
            i_eV_high = atof( is_temp.c_str() ) * 1.e9;
            
            i_pF_temp.energy_eV_max.push_back( i_eV + i_eV_high );
            i_pF_temp.energy_Hz_max.push_back( ( i_eV + i_eV_high ) / 1.239841875e-6 * TMath::C() );
            
            is_stream >> is_temp;
            i_pF_temp.flux_error_up_ergscms.push_back( atof( is_temp.c_str() ) );
            is_stream >> is_temp;
            i_pF_temp.flux_error_down_ergscms.push_back( atof( is_temp.c_str() ) );
        }
        i_photonFlux.push_back( i_pF_temp );
        
        z++;
    }
    fSpectralFlux.push_back( i_photonFlux );
    
    if( bPrint )
    {
        cout << endl;
        cout << "reading " << name << " data for MJD = (" << MJD_min << ", " << MJD_max << ")" << " from " << ifile << endl;
        cout << "------------------------------------------------------------------------------------------------------" << endl;
        for( unsigned int i = 0; i < i_photonFlux.size(); i++ )
        {
            cout << "\t";
            cout << setprecision( 8 ) << i_photonFlux[i].MJD_min << " " << i_photonFlux[i].MJD_max;
            cout << setprecision( 6 ) << " " << i_photonFlux[i].energy_Hz[0] << " " << i_photonFlux[i].energy_eV[0] << " ";
            cout << i_photonFlux[i].flux_ergscms[0] << " ";
            cout << i_photonFlux[i].flux_error_up_ergscms[0] << " ";
            cout << i_photonFlux[i].flux_error_down_ergscms[0] << " ";
            cout << endl;
        }
        cout << endl;
    }
    
    return true;
}


/*
   read complete SED for different time periods from an ascii file

   marker styles and colors are according to the appearance in the file

   Expected format:

   MJD 54624 to 54626
       <energy [Hz]>  <vFv [erg/s/cm2]>  <lower/upper error on vFv>

   values with zero errors are treated as upper flux limits
*/
bool VSpectralEnergyDistribution::readSED( string ifile )
{
    // set up the data structure
    sPhotonFlux i_pF_temp;
    vector< sPhotonFlux > i_photonFlux;
    
    // read ascii file
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading SED from " << ifile << endl;
        return false;
    }
    string is_line;
    string is_temp;
    int z = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        if( is_line.substr( 0, 1 ) == "*" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> is_temp;
        if( is_temp == "MJD" )
        {
            if( z != 0 )
            {
                i_photonFlux.push_back( i_pF_temp );
                fSpectralFlux.push_back( i_photonFlux );
                i_pF_temp.energy_Hz.clear();
                i_pF_temp.energy_Hz_min.clear();
                i_pF_temp.energy_Hz_max.clear();
                i_pF_temp.energy_eV.clear();
                i_pF_temp.energy_eV_min.clear();
                i_pF_temp.energy_eV_max.clear();
                i_pF_temp.flux_ergscms.clear();
                i_pF_temp.flux_error_up_ergscms.clear();
                i_pF_temp.flux_error_down_ergscms.clear();
                i_photonFlux.clear();
            }
            
            is_stream >> is_temp;
            i_pF_temp.MJD_min = atof( is_temp.c_str() );
            is_stream >> is_temp;
            is_stream >> is_temp;
            i_pF_temp.MJD_max = atof( is_temp.c_str() );
            
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> is_temp;
                i_pF_temp.Color = atoi( is_temp.c_str() );
            }
            else
            {
                i_pF_temp.Color = 1 + z;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> is_temp;
                i_pF_temp.Marker = atoi( is_temp.c_str() );
            }
            else
            {
                i_pF_temp.Marker = 20 + z;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                i_pF_temp.name = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
            }
            
            z++;
        }
        else
        {
            bool bXMM = false;
            if( is_temp == "keV" )
            {
                is_stream >> is_temp;
                i_pF_temp.energy_eV_min.push_back( atof( is_temp.c_str() ) * 1.e3 );
                is_stream >> is_temp;
                i_pF_temp.energy_eV_max.push_back( atof( is_temp.c_str() ) * 1.e3 );
                i_pF_temp.energy_eV.push_back( 0.5 * ( i_pF_temp.energy_eV_min.back() + i_pF_temp.energy_eV_max.back() ) );
                i_pF_temp.energy_Hz.push_back( i_pF_temp.energy_eV.back() * TMath::C() / 1.239841875e-6 );
                i_pF_temp.energy_Hz_min.push_back( i_pF_temp.energy_eV_min.back() * TMath::C() / 1.239841875e-6 );
                i_pF_temp.energy_Hz_max.push_back( i_pF_temp.energy_eV_max.back() * TMath::C() / 1.239841875e-6 );
                bXMM = true;
            }
            else
            {
                i_pF_temp.energy_Hz.push_back( atof( is_temp.c_str() ) );
                i_pF_temp.energy_Hz_min.push_back( i_pF_temp.energy_Hz.back() );
                i_pF_temp.energy_Hz_max.push_back( i_pF_temp.energy_Hz.back() );
                i_pF_temp.energy_eV.push_back( i_pF_temp.energy_Hz.back() * 1.239841875e-6 / TMath::C() );
                i_pF_temp.energy_eV_min.push_back( i_pF_temp.energy_eV.back() );
                i_pF_temp.energy_eV_max.push_back( i_pF_temp.energy_eV.back() );
            }
            
            is_stream >> is_temp;
            i_pF_temp.flux_ergscms.push_back( atof( is_temp.c_str() ) );
            is_stream >> is_temp;
            i_pF_temp.flux_error_down_ergscms.push_back( atof( is_temp.c_str() ) );
            if( TMath::Abs( i_pF_temp.flux_error_down_ergscms.back() ) < 1.e-16 )
            {
                i_pF_temp.flux_error_down_ergscms.back() = 0.;
            }
            if( ( is_stream >> std::ws ).eof() )
            {
                i_pF_temp.flux_error_up_ergscms.push_back( atof( is_temp.c_str() ) );
                if( TMath::Abs( i_pF_temp.flux_error_up_ergscms.back() ) < 1.e-16 )
                {
                    i_pF_temp.flux_error_up_ergscms.back() = 0.;
                }
            }
            else
            {
                is_stream >> is_temp;
                i_pF_temp.flux_error_up_ergscms.push_back( atof( is_temp.c_str() ) );
                if( TMath::Abs( i_pF_temp.flux_error_up_ergscms.back() ) < 1.e-16 )
                {
                    i_pF_temp.flux_error_up_ergscms.back() = 0.;
                }
            }
            // no upper limits for < X-ray frequencies (work around)
            if( i_pF_temp.energy_Hz.back() < 5.e18 && TMath::Abs( i_pF_temp.flux_error_down_ergscms.back() ) < 1.e-20 && TMath::Abs( i_pF_temp.flux_error_up_ergscms.back() ) < 1.e-20 )
            {
                i_pF_temp.flux_error_down_ergscms.back() = 1.e-19;
                i_pF_temp.flux_error_up_ergscms.back() = 1.e-19;
            }
            if( bXMM )
            {
                double iMult = 1.;
                if( i_pF_temp.energy_eV_max.back() - i_pF_temp.energy_eV_min.back() > 0. )
                {
                    iMult = i_pF_temp.energy_eV.back() / ( i_pF_temp.energy_eV_max.back() - i_pF_temp.energy_eV_min.back() );
                }
                i_pF_temp.flux_ergscms.back() *= iMult;
                i_pF_temp.flux_error_down_ergscms.back() *= iMult;
                i_pF_temp.flux_error_up_ergscms.back() *= iMult;
            }
            
        }
    }
    i_photonFlux.push_back( i_pF_temp );
    fSpectralFlux.push_back( i_photonFlux );
    
    return true;
}


/*
    read swift data from an ascii file
*/
bool VSpectralEnergyDistribution::readSwiftData( string name, string ifile, double MJD_min, double MJD_max, bool bPrint, int imarker, int icolor )
{
    return readDataFile( name, ifile, MJD_min, MJD_max, bPrint, imarker, icolor );
}


/*
    read data from an ascii file

    iname   :   data descriptor
    ifile   :   file name (ascii format)
    MJD_min :   MJD of begin of measurements
    MJD_max:    MJD of end of measurements
    bPrint  :   print values to screen
    imarker :   marker style for plotting
    icolor  :   marker and line color for plotting

all results are added to the vector of photon fluxes (sPhotonFlux)

expected file format:

<frequency (Hz)> <flux (erg cm^-2 s^-1)> <error in Flux (erg cm^-2 s^-1)>

An upper limit is drawn for fluxes with errors <= 0.

*/
bool VSpectralEnergyDistribution::readDataFile( string name, string ifile, double MJD_min, double MJD_max, bool bPrint, int imarker, int icolor )
{
    // set up the data structure
    sPhotonFlux i_pF_temp;
    vector< sPhotonFlux > i_photonFlux;
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading data from " << ifile << endl;
        return false;
    }
    string is_line;
    string is_temp;
    int z = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        if( is_line.substr( 0, 1 ) == "*" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        i_pF_temp.name = name;
        if( MJD_min > 0. && MJD_max > 0. )
        {
            i_pF_temp.MJD_min = MJD_min;
            i_pF_temp.MJD_max = MJD_max;
        }
        else
        {
            is_stream >> is_temp;
            double iT_min = atof( is_temp.c_str() );
            is_stream >> is_temp;
            double iT_max = atof( is_temp.c_str() );
            
            // check if following line is at a different time period (if yes, clear the vector and start again)
            if( z > 0 && ( TMath::Abs( i_photonFlux.back().MJD_min - iT_min ) > 1.e-2 || ( TMath::Abs( i_photonFlux.back().MJD_max - iT_max ) > 1.e-2 ) ) )
            {
                i_pF_temp.energy_Hz.clear();
                i_pF_temp.energy_Hz_min.clear();
                i_pF_temp.energy_Hz_max.clear();
                i_pF_temp.energy_eV.clear();
                i_pF_temp.energy_eV_min.clear();
                i_pF_temp.energy_eV_max.clear();
                i_pF_temp.flux_ergscms.clear();
                i_pF_temp.flux_error_up_ergscms.clear();
                i_pF_temp.flux_error_down_ergscms.clear();
            }
            i_pF_temp.MJD_min = iT_min;
            i_pF_temp.MJD_max = iT_max;
        }
        i_pF_temp.Marker = imarker;
        i_pF_temp.Color = icolor;
        
        // energy [Hz]
        is_stream >> is_temp;
        double iHz = atof( is_temp.c_str() );
        if( iHz <= 0. )
        {
            continue;
        }
        i_pF_temp.energy_Hz.push_back( iHz );
        i_pF_temp.energy_Hz_min.push_back( i_pF_temp.energy_Hz.back() );
        i_pF_temp.energy_Hz_max.push_back( i_pF_temp.energy_Hz.back() );
        double i_eV = 1.239841875e-6 / TMath::C() * iHz;
        i_pF_temp.energy_eV.push_back( i_eV );
        i_pF_temp.energy_eV_min.push_back( i_pF_temp.energy_eV.back() );
        i_pF_temp.energy_eV_max.push_back( i_pF_temp.energy_eV.back() );
        
        is_stream >> is_temp;
        i_pF_temp.flux_ergscms.push_back( atof( is_temp.c_str() ) );
        is_stream >> is_temp;
        i_pF_temp.flux_error_up_ergscms.push_back( atof( is_temp.c_str() ) );
        i_pF_temp.flux_error_down_ergscms.push_back( atof( is_temp.c_str() ) );
        
        i_photonFlux.push_back( i_pF_temp );
        
        z++;
    }
    fSpectralFlux.push_back( i_photonFlux );
    
    if( bPrint )
    {
        cout << endl;
        cout << "reading " << name << " data for MJD = (" << MJD_min << ", " << MJD_max << ")" << " from " << ifile << endl;
        cout << "------------------------------------------------------------------------------------------------------" << endl;
        for( unsigned int i = 0; i < i_photonFlux.size(); i++ )
        {
            cout << "\t";
            cout << setprecision( 8 ) << i_photonFlux[i].MJD_min << " " << i_photonFlux[i].MJD_max;
            cout << setprecision( 6 ) << " " << i_photonFlux[i].energy_Hz[i] << " " << i_photonFlux[i].energy_eV[i] << " ";
            cout << i_photonFlux[i].flux_ergscms[i] << " ";
            cout << i_photonFlux[i].flux_error_up_ergscms[i] << " ";
            cout << i_photonFlux[i].flux_error_down_ergscms[i] << " ";
            cout << endl;
        }
        cout << endl;
    }
    
    return true;
}


/*
    read XMM data

    iname   :   data descriptor
    ifile   :   file name (ascii format)
    MJD_min :   MJD of begin of measurements
    MJD_max:    MJD of end of measurements
    bPrint  :   print values to screen
    imarker :   marker style for plotting
    icolor  :   marker and line color for plotting

all results are added to the vector of photon fluxes (sPhotonFlux)

expected file format:

<frequency (keV)> <flux (cm^-2 s^-1)> <error (low) in Flux (cm^-2 s^-1)> <error (upper) in Flux (cm^-2 s^-1)>

An upper limit is drawn for fluxes with errors <= 0.
*/
bool VSpectralEnergyDistribution::readXMMData( string name, string ifile, double MJD_min, double MJD_max, bool bPrint, int imarker, int icolor, bool iModel, int iFormat )
{
    // set up the data structure
    sPhotonFlux i_pF_temp;
    vector< sPhotonFlux > i_photonFlux;
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading XMM data from " << ifile << endl;
        return false;
    }
    string is_line;
    string is_temp;
    int z = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        // use data from one instrument only
        if( is_line.substr( 0, 2 ) == "NO" )
        {
            break;
        }
        
        if( z == 0 || is_line.substr( 0, 2 ) == "NO" )
        {
            if( z != 0 )
            {
                i_photonFlux.push_back( i_pF_temp );
            }
            i_pF_temp.name = name;
            i_pF_temp.MJD_min = MJD_min;
            i_pF_temp.MJD_max = MJD_max;
            i_pF_temp.Marker = imarker;
            i_pF_temp.Color = icolor;
            
            i_pF_temp.energy_eV.clear();
            i_pF_temp.energy_Hz_min.clear();
            i_pF_temp.energy_Hz_max.clear();
            i_pF_temp.energy_eV.clear();
            i_pF_temp.energy_eV_min.clear();
            i_pF_temp.energy_eV_max.clear();
            i_pF_temp.flux_ergscms.clear();
            i_pF_temp.flux_error_up_ergscms.clear();
            i_pF_temp.flux_error_down_ergscms.clear();
        }
        
        // energy [keV]
        double i_eV = 0.;
        if( iFormat == 0 )
        {
            is_stream >> is_temp;
            i_eV = atof( is_temp.c_str() ) * 1.e3;
        }
        else if( iFormat == 1 )
        {
            is_stream >> is_temp;
            i_eV = atof( is_temp.c_str() ) * 1.e3;
            is_stream >> is_temp;
            i_eV += atof( is_temp.c_str() ) * 1.e3;
            i_eV /= 2.;
        }
        if( i_eV  <= 0. )
        {
            continue;
        }
        i_pF_temp.energy_eV.push_back( i_eV );
        i_pF_temp.energy_eV_min.push_back( i_pF_temp.energy_eV.back() );
        i_pF_temp.energy_eV_max.push_back( i_pF_temp.energy_eV.back() );
        double iHz = i_eV / 1.239841875e-6 * TMath::C();
        i_pF_temp.energy_Hz.push_back( iHz );
        i_pF_temp.energy_Hz_min.push_back( i_pF_temp.energy_Hz.back() );
        i_pF_temp.energy_Hz_max.push_back( i_pF_temp.energy_Hz.back() );
        // energy bin width (not used)
        if( iFormat == 0 )
        {
            is_stream >> is_temp;
        }
        // fluxes
        if( iFormat == 0 )
        {
            is_stream >> is_temp;
            double i_flMeasured = atof( is_temp.c_str() );
            // values are keV^2 / cm2 / keV / s -> transform to eV / cm2 / s
            // convert eV to ergs
            double i_fl = i_flMeasured * TMath::Qe() * 1.e3 / 1.e-7;
            i_pF_temp.flux_ergscms.push_back( i_fl );
            is_stream >> is_temp;
            double i_fle = atof( is_temp.c_str() )  * TMath::Qe() * 1.e3 / 1.e-7;
            i_pF_temp.flux_error_down_ergscms.push_back( i_fle );
            i_pF_temp.flux_error_up_ergscms.push_back( i_fle );
            if( iModel )
            {
                is_stream >> is_temp;
                i_fl = atof( is_temp.c_str() ) * TMath::Qe() * 1.e3 / 1.e-7;
                i_pF_temp.flux_ergscms.back() = i_fl;
            }
        }
        else if( iFormat == 1 )
        {
            is_stream >> is_temp;
            double i_flMeasured = atof( is_temp.c_str() ) * 1.e-12;
            i_pF_temp.flux_ergscms.push_back( i_flMeasured );
            is_stream >> is_temp;
            i_pF_temp.flux_error_down_ergscms.push_back( i_flMeasured - atof( is_temp.c_str() ) * 1.e-12 );
            is_stream >> is_temp;
            i_pF_temp.flux_error_up_ergscms.push_back( atof( is_temp.c_str() ) * 1.e-12 - i_flMeasured );
        }
        
        z++;
    }
    i_photonFlux.push_back( i_pF_temp );
    fSpectralFlux.push_back( i_photonFlux );
    
    if( bPrint )
    {
        cout << endl;
        cout << "reading " << name << " data for MJD = (" << MJD_min << ", " << MJD_max << ")" << " from " << ifile << endl;
        cout << "------------------------------------------------------------------------------------------------------" << endl;
        for( unsigned int i = 0; i < i_photonFlux.size(); i++ )
        {
            cout << "\t";
            cout << setprecision( 8 ) << i_photonFlux[i].MJD_min << " " << i_photonFlux[i].MJD_max;
            cout << setprecision( 6 ) << " " << i_photonFlux[i].energy_Hz[i] << " " << i_photonFlux[i].energy_eV[i] << " ";
            cout << i_photonFlux[i].flux_ergscms[i] << " ";
            cout << i_photonFlux[i].flux_error_up_ergscms[i] << " ";
            cout << i_photonFlux[i].flux_error_down_ergscms[i] << " ";
            cout << endl;
        }
        cout << endl;
    }
    
    return true;
}


/*
    read data in optical wavelengths from an ascii file

    iname   :   data descriptor
    ifile   :   file name (ascii format)
    iband   :   photometric band (e.g. R, I, J)
    bPrint  :   print values to screen
    imarker :   marker style for plotting
    icolor  :   marker and line color for plotting

    all results are added to the vector of photon fluxes (sPhotonFlux)

expected file format:

<MJD>  <magnitude>  <error in magnitude>

additional data needed:

i.  file with wavelengths and zero-points for photometric bands (e.g. $OBS_EVNDISP_AUX_DIR/AstroData/Multiwavelengthdata/photometricBands.dat)

ii. file with galactic extinction correction (for dereddening)

*/
TGraphErrors* VSpectralEnergyDistribution::readOpticalData( string iname, string ifile, string iband,
        bool bPrint, int imarker, int icolor, bool bAverage,
        double iPlotMagnitudeMultiplier, bool bCorrection, string icorfile )
{
    // set up the data structure
    sPhotonFlux i_pF_temp;
    vector< sPhotonFlux > i_photonFlux;
    
    i_pF_temp.name = iname;
    i_pF_temp.energy_Hz.push_back( getEffectiveWavelength( iband, "Hz" ) );
    i_pF_temp.energy_eV.push_back( getEffectiveWavelength( iband, "eV" ) );
    i_pF_temp.energy_Hz_min.push_back( i_pF_temp.energy_Hz.back() );
    i_pF_temp.energy_Hz_max.push_back( i_pF_temp.energy_Hz.back() );
    i_pF_temp.energy_eV_min.push_back( i_pF_temp.energy_eV.back() );
    i_pF_temp.energy_eV_max.push_back( i_pF_temp.energy_eV.back() );
    i_pF_temp.Marker = imarker;
    i_pF_temp.Color = icolor;
    i_pF_temp.MJD_min = 0.;
    i_pF_temp.MJD_max = 0.;
    double flux_ergscms = 0.;
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading optical data from " << ifile << endl;
        return 0;
    }
    string is_line;
    string is_temp;
    
    TGraphErrors* g = new TGraphErrors( 1 );
    g->SetMarkerStyle( imarker );
    g->SetMarkerColor( icolor );
    g->SetLineColor( icolor );
    
    int z = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        if( is_line.substr( 0, 1 ) == "*" )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> is_temp;
        if( i_pF_temp.MJD_min > 2400000. )
        {
            i_pF_temp.MJD_min -= 2400000.5;
        }
        double MJD_min = atof( is_temp.c_str() );
        if( MJD_min > 2400000. )
        {
            MJD_min -= 2400000.5;
        }
        is_stream >> is_temp;
        double imagnitude = atof( is_temp.c_str() );
        is_stream >> is_temp;
        double imagnitudeError =  atof( is_temp.c_str() );
        
        if( bCorrection )
        {
            readGalacticExtinction( icorfile, bPrint );
            double icor = getGalacticExtinctionCorrection( iband );
            imagnitude -= icor;
        }
        
        
        if( !bAverage )
        {
            i_pF_temp.MJD_min = MJD_min;
            i_pF_temp.MJD_max = MJD_min;
            i_pF_temp.flux_ergscms.push_back( getFluxfromMagnitude( imagnitude, iband ) * i_pF_temp.energy_Hz.back() );
            double imagnitude_H = imagnitude + imagnitudeError;
            i_pF_temp.flux_error_up_ergscms.push_back( i_pF_temp.flux_ergscms.back() - getFluxfromMagnitude( imagnitude_H, iband ) * i_pF_temp.energy_Hz.back() );
            imagnitude_H = imagnitude - atof( is_temp.c_str() );
            i_pF_temp.flux_error_down_ergscms.push_back( getFluxfromMagnitude( imagnitude_H, iband ) * i_pF_temp.energy_Hz.back() - i_pF_temp.flux_ergscms.back() );
            
            i_photonFlux.push_back( i_pF_temp );
        }
        else
        {
            i_pF_temp.MJD_min += MJD_min;
            i_pF_temp.MJD_max += MJD_min;
            flux_ergscms      += getFluxfromMagnitude( imagnitude, iband ) * i_pF_temp.energy_Hz.back();
        }
        
        g->SetPoint( z, i_pF_temp.MJD_min, imagnitude * iPlotMagnitudeMultiplier );
        g->SetPointError( z, 0., imagnitudeError );
        
        z++;
    }
    
    if( bAverage )
    {
        if( z > 0 )
        {
            i_pF_temp.MJD_min /= z;
            i_pF_temp.MJD_max /= z;
            i_pF_temp.flux_ergscms.push_back( flux_ergscms / z );
            // no error on fluxes for average values (keep it small and non-zero)
            i_pF_temp.flux_error_up_ergscms.push_back( 1.e-20 );
            i_pF_temp.flux_error_down_ergscms.push_back( 1.e-20 );
            i_photonFlux.push_back( i_pF_temp );
        }
    }
    
    fSpectralFlux.push_back( i_photonFlux );
    
    if( bPrint )
    {
        cout << endl;
        cout << "reading optical data for " << iname << " from " << ifile;
        cout << "  (" << iband << " band )" << endl;
        cout << "------------------------------------------------------------------------------------------------------" << endl;
        for( unsigned int i = 0; i < i_photonFlux.size(); i++ )
        {
            cout << "\t";
            cout << setprecision( 8 ) << i_photonFlux[i].MJD_min;
            cout << setprecision( 6 ) << " " << i_photonFlux[i].energy_Hz[i] << " " << i_photonFlux[i].energy_eV[i] << " ";
            cout << i_photonFlux[i].flux_ergscms[i] << " ";
            cout << i_photonFlux[i].flux_error_up_ergscms[i] << " ";
            cout << i_photonFlux[i].flux_error_down_ergscms[i] << " ";
            cout << endl;
        }
        cout << endl;
    }
    
    return g;
}


/*!
    read wavelengths and zero-points for photometric bands from files

    ifile   :   file name (ascii format)

    expected file format:

    <band>  <effective wavelenght in micron> <F_0 [Jy] (CIT)> <<F_0 [Jy] (UKIRT)>

*/
bool VSpectralEnergyDistribution::readPhotoMetricBands( string ifile, bool iPrint )
{
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VSpectralEnergyDistribution: error reading photometric bands from " << ifile << endl;
        return false;
    }
    string is_line;
    string is_temp;
    
    fPhotoMetricBand.clear();
    sPhotoMetricBand iTemp;
    
    if( iPrint )
    {
        cout << "VSpectralEnergyDistribution: reading photometric bands from " << ifile << endl;
    }
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> iTemp.fBand;
        if( iTemp.fBand == "*" )
        {
            continue;
        }
        is_stream >> is_temp;
        iTemp.fEffectiveWavelength_micron = atof( is_temp.c_str() );
        is_stream >> is_temp;
        iTemp.fF0_CIT_Jy = atof( is_temp.c_str() );
        is_stream >> is_temp;
        iTemp.fF0_UKIRT_Jy = atof( is_temp.c_str() );
        if( TMath::Abs( iTemp.fF0_UKIRT_Jy ) < 1.e-8 )
        {
            iTemp.fF0_UKIRT_Jy = 0.;
        }
        
        fPhotoMetricBand.push_back( iTemp );
    }
    
    if( iPrint )
    {
        cout << endl;
        cout << "\teffective wavelengths and zero-points for photometric bands" << endl;
        cout << "\tband \t lambda [micron] \t F_0 [Jy] (CIT) \t F_0 [Jy] (UKIRT)" << endl;
        cout << "\t---------------------------------------------------------------------------------------------------------------" << endl;
        for( unsigned int i = 0; i < fPhotoMetricBand.size(); i++ )
        {
            cout << "\t";
            cout << fPhotoMetricBand[i].fBand << "\t\t" << fPhotoMetricBand[i].fEffectiveWavelength_micron << "\t\t";
            cout << fPhotoMetricBand[i].fF0_CIT_Jy << "\t\t\t" << fPhotoMetricBand[i].fF0_UKIRT_Jy << endl;
        }
    }
    
    return true;
}


/*
       get effective wavelength or energy for a given photometric band

       iBand    :  photometric band (e.g. "V")
       iUnit    :  frequency ("Hz") or energy ("eV")

*/
double VSpectralEnergyDistribution::getEffectiveWavelength( string iband, string iUnit )
{
    double iWL = 0.;
    for( unsigned int i = 0; i < fPhotoMetricBand.size(); i++ )
    {
        if( fPhotoMetricBand[i].fBand == iband )
        {
            iWL = fPhotoMetricBand[i].fEffectiveWavelength_micron;
            break;
        }
    }
    // micrometer to meter
    iWL /= 1.e6;
    
    if( fDebug )
    {
        cout << "getEffectiveWavelength: Wavelength: [m] " << iWL << " " << iband << " " << iUnit << endl;
    }
    
    if( iWL <= 0. )
    {
        return -99.;
    }
    
    if( iUnit == "eV" )
    {
        iWL = 1.239841875e-6 / iWL;               // m
        if( fDebug )
        {
            cout << "\t eV " << iWL << endl;
        }
        return iWL;
    }
    else if( iUnit == "Hz" )
    {
        iWL = TMath::C() / iWL;
        if( fDebug )
        {
            cout << "\t Hz " << iWL << endl;
        }
        return iWL;
    }
    else
    {
        cout << "error: unknown unit in getEffectiveWavelength: " << iUnit << endl;
        return -99.;
    }
    
    return 0.;
}



/*
  Copy the galactic extinction numbers from NED
  The file should then look like this:

U     B     V     R     I     J     H     K     L'
0.130 0.103 0.079 0.064 0.046 0.022 0.014 0.009 0.004

*/
bool VSpectralEnergyDistribution::readGalacticExtinction( string ifile, bool iPrint )
{
    ifstream is( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading file" << ifile << endl;
        return false;
    }
    if( iPrint )
    {
        cout << "successfully opened file " << ifile << endl;
    }
    
    //    vector< sGalacticExtinction > fGalacticExtinction;
    sGalacticExtinction iGalExt;
    
    string iband[9];
    double icorrection[9];
    
    string is_line;
    string is_temp;
    
    int i = 0;
    while( getline( is, is_line ) )
    {
        //	if( iPrint ) cout << is_line << endl;
        stringstream is_stream( is_line );
        
        int k = 0;
        while( !( is_stream >> std::ws ).eof() )
        {
            is_stream >> is_temp;
            if( i == 0 )
            {
                iband[k] = is_temp.c_str();
            }
            if( i == 1 )
            {
                icorrection[k] = atof( is_temp.c_str() );
            }
            k++;
        }
        i++;
    }
    
    int length = sizeof( iband ) / sizeof( string );
    
    for( int l = 0; l < length; l++ )
    {
        iGalExt.fBand = iband[l];
        iGalExt.fCorrection = icorrection[l];
        
        fGalacticExtinction.push_back( iGalExt );
    }
    
    if( iPrint )
    {
        cout << endl;
        cout << "\tgalactic extinction corrections " << endl;
        cout << "\tband \t correction [mag]" << endl;
        cout << "\t--------------------------------" << endl;
        for( unsigned int i = 0; i < fGalacticExtinction.size(); i++ )
        {
            cout << "\t";
            cout << fGalacticExtinction[i].fBand << "\t\t"
                 << fGalacticExtinction[i].fCorrection << endl;
        }
    }
    
    return true;
}



double VSpectralEnergyDistribution::getGalacticExtinctionCorrection( string iband )
{
    double cor = 0;
    
    for( unsigned int i = 0; i < fGalacticExtinction.size(); i++ )
    {
        if( fGalacticExtinction[i].fBand == iband )
        {
            cor = fGalacticExtinction[i].fCorrection;
            break;
        }
    }
    
    return cor;
}


/*
    get photometric flux in erg/s/m2/Hz

    magnitude   :   magnitude
    iband       :   photometric band
    system      :   system for flux conversion ("CIT" or "UKIRT" )
*/
double VSpectralEnergyDistribution::getFluxfromMagnitude( double magnitude, string iband, string system )
{
    double f_0 = 0.;
    for( unsigned int i = 0; i < fPhotoMetricBand.size(); i++ )
    {
        if( fPhotoMetricBand[i].fBand == iband )
        {
            if( system == "CIT" )
            {
                f_0 = fPhotoMetricBand[i].fF0_CIT_Jy;
            }
            if( system == "UKIRT" )
            {
                f_0 = fPhotoMetricBand[i].fF0_UKIRT_Jy;
            }
            break;
        }
    }
    if( f_0 <= 0. )
    {
        return 0.;
    }
    
    double f_vu = 0;
    
    // calculate spectral flux density (see Skinner 1996 page 3)
    f_vu = TMath::Power( 10., -1.* magnitude / 2.5 ) * f_0;
    // convert Jansky to erg/s/m2/Hz
    f_vu *= 1.e-19;
    // convert from 1/m2 to 1/cm2
    f_vu *= 1.e-4;
    
    return f_vu;
}

TCanvas* VSpectralEnergyDistribution::plotCanvas( int canvas_x, int canvas_y )
{
    char hname[800];
    sprintf( hname, "c_%s", fName.c_str() );
    
    TCanvas* cSP = new TCanvas( hname, fName.c_str(), 10, 10, canvas_x, canvas_y );
    cSP->SetGridx( 0 );
    cSP->SetGridy( 0 );
    cSP->SetLogy( 1 );
    if( canvas_x <= 600 )
    {
        cSP->SetLeftMargin( 0.15 );
    }
    cSP->Draw();
    
    sprintf( hname, "h_%s", fName.c_str() );
    TH1D* hnull = new TH1D( hname, "", 100, log10( fPlotting_EnergyRange_min_Hz ), log10( fPlotting_EnergyRange_max_Hz ) );
    hnull->SetStats( 0 );
    hnull->SetXTitle( "log_{10} #nu [Hz]" );
    hnull->SetYTitle( "#nu F_{#nu} [erg s^{-1} cm^{-2}]" );
    hnull->SetMinimum( fPlotting_FluxRange_min );
    hnull->SetMaximum( fPlotting_FluxRange_max );
    if( canvas_x <= 600 )
    {
        hnull->GetYaxis()->SetTitleOffset( 1.6 );
    }
    hnull->Draw();
    
    return cSP;
}

/*
    plot all data into one SED

    c    :   plot into this canvas (default = 0: create new canvas)

*/
TCanvas* VSpectralEnergyDistribution::plot( TCanvas* c, int bLegend, int canvas_x, int canvas_y, bool bPlotErrorX, bool bPlotName )
{

    TCanvas* cSP = c;
    if( cSP == 0 )
    {
        cSP = plotCanvas( canvas_x, canvas_y );
    }
    TLegend* iL = 0;
    if( bLegend == 1 )
    {
        iL = new TLegend( 0.75, 0.75, 0.95, 0.95 );
    }
    else if( bLegend == 2 )
    {
        iL = new TLegend( 0.31, 0.77, 0.5, 0.85 );
    }
    else if( bLegend == 3 )
    {
        iL = new TLegend( 0.5, 0.73, 0.95, 0.95 );
    }
    
    // set up graphs and draw them
    for( unsigned int i = 0; i < fSpectralFlux.size(); i++ )
    {
        if( fSpectralFlux[i].size() == 0 )
        {
            continue;
        }
        
        for( unsigned int j = 0; j < fSpectralFlux[i].size(); j++ )
        {
            // check if observations are in allowed MJD range
            if( fSpectralFlux[i][j].MJD_min < fMJDMin && fSpectralFlux[i][j].MJD_max < fMJDMin )
            {
                continue;
            }
            if( fSpectralFlux[i][j].MJD_min > fMJDMax && fSpectralFlux[i][j].MJD_max > fMJDMax )
            {
                continue;
            }
            
            TGraphAsymmErrors* g = new TGraphAsymmErrors( 1 );
            g->SetLineColor( fSpectralFlux[i][j].Color );
            g->SetMarkerColor( fSpectralFlux[i][j].Color );
            g->SetMarkerSize( 1 );
            g->SetMarkerStyle( fSpectralFlux[i][j].Marker );
            
            // loop over all data points
            int z = 0;
            for( unsigned int k = 0; k < fSpectralFlux[i][j].energy_Hz.size(); k++ )
            {
                if( fSpectralFlux[i][j].energy_Hz[k] > 0. && fSpectralFlux[i][j].flux_error_up_ergscms[k] > 0. )
                {
                    g->SetPoint( z, log10( fSpectralFlux[i][j].energy_Hz[k] ), fSpectralFlux[i][j].flux_ergscms[k] );
                    g->SetPointEYhigh( z, fSpectralFlux[i][j].flux_error_up_ergscms[k] );
                    g->SetPointEYlow( z, fSpectralFlux[i][j].flux_error_down_ergscms[k] );
                    if( bPlotErrorX )
                    {
                        g->SetPointEXlow( z, log10( fSpectralFlux[i][j].energy_Hz[k] ) - log10( fSpectralFlux[i][j].energy_Hz_min[k] ) );
                        g->SetPointEXhigh( z, log10( fSpectralFlux[i][j].energy_Hz_max[k] ) - log10( fSpectralFlux[i][j].energy_Hz[k] ) );
                    }
                    z++;
                }
                else if( fSpectralFlux[i][j].flux_error_up_ergscms[k] <= 0. )
                {
                    TArrow* i_aw = new TArrow( log10( fSpectralFlux[i][j].energy_Hz[k] ), fSpectralFlux[i][j].flux_ergscms[k], log10( fSpectralFlux[i][j].energy_Hz[k] ), 0.5 * fSpectralFlux[i][j].flux_ergscms[k],  0.01, "|-|>" );
                    i_aw->SetLineColor( fSpectralFlux[i][j].Color );
                    i_aw->SetFillColor( fSpectralFlux[i][j].Color );
                    i_aw->Draw();
                }
            }
            if( z > 0 )
            {
                g->Draw( "p" );
            }
            
            if( bLegend )
            {
                char hname[600];
                if( bPlotName )
                {
                    iL->AddEntry( g, fSpectralFlux[i][j].name.c_str(), "pl" );
                }
                else
                {
                    sprintf( hname, "MJD %.1f-%.1f", fSpectralFlux[i][j].MJD_min, fSpectralFlux[i][j].MJD_max );
                    iL->AddEntry( g, hname, "pl" );
                }
            }
        }
    }
    if( bLegend )
    {
        iL->Draw();
    }
    
    return cSP;
}


void VSpectralEnergyDistribution::printASCII()
{
    cout << "MJD " << fMJDMin << " to " << fMJDMax << endl;
    
    multimap< double, double > iFlux;
    multimap< double, double > iFlux_U;
    multimap< double, double > iFlux_D;
    
    ////////////////////////////////////
    // loop over all data sets
    for( unsigned int i = 0; i < fSpectralFlux.size(); i++ )
    {
        if( fSpectralFlux[i].size() == 0 )
        {
            continue;
        }
        
        for( unsigned int j = 0; j < fSpectralFlux[i].size(); j++ )
        {
            // check if observations are in allowed MJD range
            if( fSpectralFlux[i][j].MJD_min < fMJDMin && fSpectralFlux[i][j].MJD_max < fMJDMin )
            {
                continue;
            }
            if( fSpectralFlux[i][j].MJD_min > fMJDMax && fSpectralFlux[i][j].MJD_max > fMJDMax )
            {
                continue;
            }
            
            iFlux.insert( make_pair( fSpectralFlux[i][j].energy_Hz[j], fSpectralFlux[i][j].flux_ergscms[j] ) );
            
            if( fSpectralFlux[i][j].flux_error_up_ergscms[j] < 1.e-19 )
            {
                iFlux_U.insert( make_pair( fSpectralFlux[i][j].energy_Hz[j], 0. ) );
            }
            else
            {
                iFlux_U.insert( make_pair( fSpectralFlux[i][j].energy_Hz[j], fSpectralFlux[i][j].flux_error_up_ergscms[j] ) );
            }
            
            if( fSpectralFlux[i][j].flux_error_down_ergscms[j] < 1.e-19 )
            {
                iFlux_D.insert( make_pair( fSpectralFlux[i][j].energy_Hz[j], 0. ) );
            }
            else
            {
                iFlux_D.insert( make_pair( fSpectralFlux[i][j].energy_Hz[j], fSpectralFlux[i][j].flux_error_down_ergscms[j] ) );
            }
        }
    }
    
    // print everthing sorted by energy
    // (terrible coding...)
    vector< string > itemp;
    char hname[600];
    
    typedef multimap< double, double >::const_iterator I;
    for( I p = iFlux.begin(); p != iFlux.end(); ++p )
    {
        sprintf( hname, "%.5e \t %.5e \t", p->first, p->second );
        itemp.push_back( hname );
    }
    unsigned int z = 0;
    for( I p = iFlux_U.begin(); p != iFlux_U.end(); ++p )
    {
        sprintf( hname, "%.5e \t", p->second );
        if( z < itemp.size() )
        {
            itemp[z] += hname;
        }
        z++;
    }
    z = 0;
    for( I p = iFlux_D.begin(); p != iFlux_D.end(); ++p )
    {
        sprintf( hname, "%.5e \t", p->second );
        if( z < itemp.size() )
        {
            itemp[z] += hname;
        }
        z++;
    }
    
    for( unsigned int i = 0; i < itemp.size(); i++ )
    {
        cout << itemp[i] << endl;
    }
    
}


TGraph* VSpectralEnergyDistribution::plotModel( TCanvas* c, string ifile, int icolor, int ilinestyle, int ilinewidth, bool isJyHz )
{
    if( !c || ifile.size() == 0 )
    {
        return 0;
    }
    
    TGraph* g = new TGraph( 1 );
    g->SetLineColor( icolor );
    g->SetLineStyle( ilinestyle );
    g->SetLineWidth( ilinewidth );
    
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading model file : " << ifile << endl;
        return 0;
    }
    string is_line;
    string is_temp;
    string is_temp2;
    int z = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() <= 0 )
        {
            continue;
        }
        
        istringstream is_stream( is_line );
        
        is_stream >> is_temp;
        is_stream >> is_temp2;
        
        if( atof( is_temp2.c_str() ) > 5.e-50 && atof( is_temp.c_str() ) > 0. )
        {
            double flux;
            if( isJyHz ) // convert from JyHz to erg s^-1 cm^-2
            {
                flux = atof( is_temp2.c_str() ) * 1e-23 ;
            }
            else
            {
                flux = atof( is_temp2.c_str() );
            }
            g->SetPoint( z, log10( atof( is_temp.c_str() ) ), flux );
            z++;
        }
    }
    is.close();
    
    c->cd();
    
    g->Draw( "c" );
    
    return g;
}

/*!

   plot a power law function on top of the current canvas

   (differential flux in photons / m^2 / s / TeV)

*/
TCanvas* VSpectralEnergyDistribution::plotPowerLaw( TCanvas* c, string iName, double iEMin_TeV, double iEMax_TeV,
        double iNorm, double iGamma, double iNormEnergy_TeV,
        bool bPlotButterfly, double iNormError, double iGammaError,
        int iLineColor, int iLineStyle )
{
    TCanvas* cSP = c;
    if( cSP == 0 )
    {
        cSP = plotCanvas();
    }
    cSP->cd();
    
    TF1* f = new TF1( iName.c_str(), "[0] * TMath::H() * 1.e7 * TMath::Power( 4.14e-27 * TMath::Power( 10., x )/[2], [1] ) * TMath::Power( 10., x )",
                      log10( VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( iEMin_TeV ) ),
                      log10( VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( iEMax_TeV ) ) );
    f->SetParameter( 0, iNorm );
    f->SetParameter( 1, -1.*iGamma + 1. );
    f->SetParameter( 2, iNormEnergy_TeV );
    f->SetLineColor( iLineColor );
    f->SetLineStyle( iLineStyle );
    
    // plot power law only
    if( !bPlotButterfly )
    {
        f->Draw( "same" );
    }
    // plot a butterfly taking the error on flux normalization and index into account
    // NOTE: assume that both are uncorrelated;
    //       if correlated: do not use, error will be underestimated
    //
    else
    {
        TGraphErrors* g = new TGraphErrors( 1 );
        g->SetLineColor( iLineColor );
        g->SetFillColor( iLineColor );
        for( int i = 0; i < 10; i++ )
        {
            double e = TMath::Power( 10., log10( iEMin_TeV ) + i * ( log10( iEMax_TeV ) - log10( iEMin_TeV ) ) / 10. );
            e = VFluxAndLightCurveUtilities::convertEnergy_TeV_to_Hz( e );
            
            double fl = f->Eval( log10( e ) );
            
            double c1 = TMath::H() * 1.e7 * TMath::Power( 4.14e-27 * e / iNormEnergy_TeV, -1.*iGamma + 1. ) * e;
            double c2 = 4.14e-27 * e / iNormEnergy_TeV;
            
            double er = 0.;
            
            er += c1 * c1 * iNormError * iNormError;
            er += ( iNorm * TMath::H() * 1.e7 * e * log( c2 ) * TMath::Power( c2, -1.*iGamma + 1. ) )
                  * ( iNorm * TMath::H() * 1.e7 * e * log( c2 ) * TMath::Power( c2, -1.*iGamma + 1. ) )
                  * iGammaError * iGammaError;
                  
            er = sqrt( er );
            
            g->SetPoint( i, log10( e ), fl );
            g->SetPointError( i, 0., er );
        }
        
        g->Draw( "3" );
        g->Print();
    }
    
    return cSP;
}
