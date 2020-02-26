/*! \class VGlobalRunParameter
    \brief global definitions for eventdisplay package

*/

#include "VGlobalRunParameter.h"

VGlobalRunParameter::VGlobalRunParameter( bool bSetGlobalParameter )
{
    // read global parameters
    if( bSetGlobalParameter )
    {
        if( !bReadRunParameter )
        {
            // set directories
            setDirectories();
            if( !readRunparameterFile( getDirectory_EVNDISPParameterFiles() + "EVNDISP.global.runparameter" ) )
            {
                cout << "VGlobalRunParameter: error while reading parameter file with global run parameters" << endl;
                cout << endl;
                cout << "Note: several parameter files are needed to run the eventdisplay analysis" << endl;
                cout << "      they are provided in form of a tar ball and are observatory dependent (CTA or VERITAS)" << endl;
                cout << "The file not found is EVNDISP.global.runparameter" << endl;
                cout << endl;
                cout << "Parameter files are expected to be in the following directory: " << endl;
                cout << getDirectory_EVNDISPParameterFiles() << endl;
                exit( EXIT_FAILURE );
            }
            else
            {
                bReadRunParameter = true;
            }
        }
    }
}

VGlobalRunParameter::~VGlobalRunParameter()
{

}

bool VGlobalRunParameter::readRunparameterFile( string i_filename )
{
    ifstream is;
    is.open( i_filename.c_str(), ifstream::in );
    if( !is )
    {
        const char* evn_dir = gSystem->Getenv( "EVNDISPDATA" );
        if( evn_dir )
        {
            string itemp = evn_dir;
            itemp += "/" + i_filename;
            is.open( itemp.c_str(), ifstream::in );
            if( !is )
            {
                return false;
            }
            i_filename = itemp;
        }
    }
    string is_line;
    string temp;
    string temp2;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            // print runparameter to stdout
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> temp;
                if( temp == "OBSERVATORY" )
                {
                    fObservatory = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
                }
                else if( temp == "OBSERVATORY_COORDINATES" )
                {
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> fObservatory_Latitude_deg;
                    }
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> fObservatory_Longitude_deg;
                    }
                    if( !(is_stream>>std::ws).eof() )
                    {
                        is_stream >> fObservatory_Height_m;
                    }
                }
                else if( temp == "DBSERVER" )
                {
                    fDBServer = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
                    // remove all white spaces
                    string iTT;
                    for( unsigned int i = 0; i < fDBServer.length(); i++ )
                    {
                        if( fDBServer[i] != ' ' )
                        {
                            iTT += fDBServer[i];
                        }
                    }
                    fDBServer = iTT;
                }
                else if( temp == "VTSRAWDATA" )
                {
                    fRawDataServer  = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
                    // remove all white spaces
                    string iTT;
                    for( unsigned int i = 0; i < fRawDataServer.length(); i++ )
                    {
                        if( fRawDataServer[i] != ' ' )
                        {
                            iTT += fRawDataServer[i];
                        }
                    }
                    fRawDataServer = iTT;
                }
            }
            else
            {
                cout << "VGlobalRunParameter::readRunparameterFile error while reading parameters" << endl;
                return false;
            }
        }
    }
    
    is.close();
    
    return true;
}

/*!

   set calibration directory
   (require full path here)

*/
bool VGlobalRunParameter::setDirectory_EVNDISPCalibrationData( string iDir )
{
    fEVNDISPCalibrationDataDirectory = iDir;
    if( gSystem->AccessPathName( fEVNDISPCalibrationDataDirectory.c_str() ) )
    {
        cout << "VGlobalRunParameter::setDirectory_EVNDISPCalibrationData(): creating calibration directory in " << endl;
        cout << "   " << fEVNDISPCalibrationDataDirectory << endl;
        gSystem->MakeDirectory( fEVNDISPCalibrationDataDirectory.c_str() );
        if( gSystem->AccessPathName( fEVNDISPCalibrationDataDirectory.c_str() ) )
        {
            cout << "VGlobalRunParameter::setDirectory_EVNDISPCalibrationData(): error creating calibration directory in " << endl;
            cout << "   " << fEVNDISPCalibrationDataDirectory << endl;
            return false;
        }
    }
    cout << "VGlobalRunParameter::setDirectory_EVNDISPCalibrationData(): calibration directory is " << endl;
    cout << "   " << fEVNDISPCalibrationDataDirectory << endl;
    
    return true;
}

bool VGlobalRunParameter::setDirectories()
{
    //////////////////////////////////////////////////////////////////////
    // get directories with all analysis data
    string data_dir = "";
    if( gSystem->Getenv( "OBS_EVNDISP_AUX_DIR" ) )
    {
        data_dir = gSystem->Getenv( "OBS_EVNDISP_AUX_DIR" );
    }
    // assume that by default a VERITAS analysis is requested
    else if( gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" ) )
    {
        data_dir = gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" );
    }
    if( data_dir.size() != 0 )
    {
        fEVNDISPAnaDataDirectory = data_dir;
        fEVNDISPAnaDataDirectory += "/";
    }
    // test if directory exists
    if( gSystem->AccessPathName( fEVNDISPAnaDataDirectory.c_str() ) )
    {
        cout << "VGlobalRunParameter::setDirectories(): cannot find directory with EVNDISP aux data" << endl;
        cout << "\t looking for " << fEVNDISPAnaDataDirectory << endl;
        cout << "\t is environmental variable $OBS_EVNDISP_AUX_DIR (or $VERITAS_EVNDISP_AUX_DIR or $CTA_EVNDISP_AUX_DIR) set?" << endl;
        cout << "\t (see README/INSTALL)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    // by default: calibration directory = fEVNDISPAnaDataDirectory
    if( fEVNDISPCalibrationDataDirectory.size() == 0 )
    {
        fEVNDISPCalibrationDataDirectory = fEVNDISPAnaDataDirectory;
    }
    //////////////////////////////////////////////////////////////////////
    // raw data is expected to be here
    const char* raw_dir = gSystem->Getenv( "VERITAS_DATA_DIR" );
    if( raw_dir )
    {
        fVBFRawDataDirectory = raw_dir;
    }
    else
    {
        fVBFRawDataDirectory = "./data/";
    }
    
    //////////////////////////////////////////////////////////////////////
    // output is written to this directory (unless stated otherwise on command line)
    const char* cal_out = gSystem->Getenv( "OBS_USER_DATA_DIR" );
    if( cal_out )
    {
        fEVNDISPOutputDirectory = cal_out;
        fEVNDISPOutputDirectory += "/";
    }
    else
    {
        const char* vcal_out = gSystem->Getenv( "VERITAS_USER_DATA_DIR" );
        if( vcal_out )
        {
            fEVNDISPOutputDirectory = vcal_out;
            fEVNDISPOutputDirectory += "/";
        }
    }
    if( gSystem->AccessPathName( fEVNDISPOutputDirectory.c_str() ) )
    {
        cout << "VGlobalRunParameter::setDirectories(): cannot find directory for EVNDISP output data" << endl;
        cout << "\t looking for " << fEVNDISPOutputDirectory << endl;
        cout << "\t is environmental variable $OBS_USER_DATA_DIR set?" << endl;
        cout << "\t (see README/INSTALL)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    return true;
}


/*
   get eventdisplay tree version
*/
unsigned int VGlobalRunParameter::getEVNDISP_TREE_VERSION( TTree* t )
{
    if( !t )
    {
        return 0;
    }
    
    string itemp = t->GetTitle();
    if( itemp.find( "VERSION" ) < itemp.size() )
    {
        return ( unsigned int )( atoi( itemp.substr( itemp.find( "VERSION" ) + 7, itemp.size() ).c_str() ) );
    }
    return 0;
}

bool VGlobalRunParameter::getEVNDISP_TREE_isShort( TTree* t )
{
    if( !t )
    {
        return false;
    }
    
    string itemp = t->GetTitle();
    if( itemp.find( "short tree" ) < itemp.size() )
    {
        return true;
    }
    
    return false;
}

/*!

*/
bool VGlobalRunParameter::update( TChain* c )
{
    if( !c )
    {
        return false;
    }
    
    fEVNDISP_TREE_VERSION = getEVNDISP_TREE_VERSION( c->GetTree() );
    
    return true;
}

void VGlobalRunParameter::printGlobalRunParameter()
{
    if( !bReadRunParameter )
    {
        cout << "VGlobalRunParameter::printGlobalRunParameter(): no global run parameters read" << endl;
        return;
    }
    
    cout << "reading global run parameters from " << getDirectory_EVNDISPParameterFiles();
    cout << "EVNDISP.global.runparameter" << endl;
    cout << endl;
    cout << "VERSION " << fEVNDISP_VERSION << " (tree version " << fEVNDISP_TREE_VERSION << ")";
    cout << endl;
    cout << "Observatory: " << fObservatory;
    if( TMath::Abs( fObservatory_Longitude_deg ) > 1.e-5
            && TMath::Abs( fObservatory_Latitude_deg ) > 1.e-5
            && TMath::Abs( fObservatory_Height_m ) > 1.e-5 )
    {
        cout << " (long " << fObservatory_Longitude_deg << " deg, lat " << fObservatory_Latitude_deg << " deg,";
        cout << " altitude " << fObservatory_Height_m << " m)";
    }
    cout << endl;
    if( fDBServer.size() > 0 )
    {
        cout << "DB server " << fDBServer << endl;
    }
    if( fRawDataServer.size() > 0 )
    {
        cout << "Raw data server " << fRawDataServer << endl;
    }
    cout << "Directories: " << endl;
    if( fEVNDISPAnaDataDirectory.size() > 0 )
    {
        cout << "    for EVNDISP aux data: \t\t" << fEVNDISPAnaDataDirectory << endl;
    }
    if( fVBFRawDataDirectory.size() > 0 )
    {
        cout << "    for VBF raw: \t\t\t" << fVBFRawDataDirectory << endl;
    }
    if( getDirectory_EVNDISPCalibrationData().size() > 0 )
    {
        cout << "    for Calibration data: \t\t" << getDirectory_EVNDISPCalibrationData() << endl;
    }
    if( getDirectory_EVNDISPCalibrationData_perRun().size() > 0 )
    {
        cout << "    for Calibration data (per run): \t" << getDirectory_EVNDISPCalibrationData_perRun() << endl;
    }
    if( fEVNDISPOutputDirectory.size() > 0 )
    {
        cout << "    for EVNDISP output: \t\t" << fEVNDISPOutputDirectory << endl;
    }
    cout << endl;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////

string VGlobalRunParameter::fObservatory = "Whipple";
bool VGlobalRunParameter::bReadRunParameter = false;
unsigned int VGlobalRunParameter::fEVNDISP_TREE_VERSION = 9;
string VGlobalRunParameter::fEVNDISP_VERSION = "g500";
string VGlobalRunParameter::fDBServer = "";
string VGlobalRunParameter::fRawDataServer = "";
string VGlobalRunParameter::fEVNDISPAnaDataDirectory = "";
string VGlobalRunParameter::fEVNDISPCalibrationDataDirectory = "";
string VGlobalRunParameter::fVBFRawDataDirectory = "";
string VGlobalRunParameter::fEVNDISPOutputDirectory = "";
double VGlobalRunParameter::fObservatory_Longitude_deg = 0.;
double VGlobalRunParameter::fObservatory_Latitude_deg = 0.;
double VGlobalRunParameter::fObservatory_Height_m = 0.;


