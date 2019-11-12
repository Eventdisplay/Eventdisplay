/*! \file  makeDISPTables.cpp
    \brief generate tables for angular reconstruction (disp method)


    stereo reconstruction with disp method, see Hofmann et al (1999) paper

    creates tables (disp vs log size vs width/length) for each ze, az, woff, pedvar

    tables are then used in VArrayAnalyzer for direction reconstruction

*/

#include "VGlobalRunParameter.h"
#include "VDispTable.h"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

#include <TStopwatch.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////
// parameters to be read in from runparameter file
// (here defaults are defined; are overwritting after reading run parameter file)
////////////////////////////////////////////////////////////////////////////
// number of telescopes
unsigned int fNTel = 4;
// total number of events to use for table generation (per simulation file; -1 = use all events)
int fNTotEvents = -1;
// quality cuts
int f_ntubes_min = 4;
double f_size_min = 0.;
double f_length_min = 0.;
double f_loss_max = 0.5;
// scaling parameters
double fLengthScaleParameter = 0.21;
double fWidthScaleParameter  = 0.08;
// input directory pattern
string fInputDir;
// zenith angle bins
vector< string > f_ze;
// wobble offsets
vector< string > f_woff;
// noise bins (grisudet units)
vector< string > f_noise;
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
/*!
      see makeDISPTables.input.dat for an example of a runparameter file
*/
////////////////////////////////////////////////////////////////////////////
bool readInputParameter( string i_filename )
{
    ifstream is;
    is.open( i_filename.c_str(), ifstream::in );
    if( !is )
    {
        cout << "no file found to read run parameters: " << i_filename << endl;
        exit( 0 );
    }
    cout << "reading run parameters from " << i_filename << " :" << endl;
    cout << endl;
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
            cout << is_line << endl;
            if( is_stream.eof() )
            {
                cout << "error reading runparameter file" << endl;
                return false;
            }
            is_stream >> temp;
            if( is_stream.eof() )
            {
                cout << "error reading runparameter file" << endl;
                return false;
            }
            if( temp == "SIMUDIRECTORY" )
            {
                is_stream >> fInputDir;
            }
            else if( temp == "WIDTHSCALING" )
            {
                is_stream >> temp;
                fWidthScaleParameter = atoi( temp.c_str() );
            }
            else if( temp == "LENGTHSCALING" )
            {
                is_stream >> temp;
                fLengthScaleParameter = atof( temp.c_str() );
            }
            else if( temp == "USE_N_EVENTS" )
            {
                is_stream >> temp;
                fNTotEvents = atoi( temp.c_str() );
            }
            else if( temp == "NUMBEROFTELESCOPES" )
            {
                is_stream >> temp;
                fNTel = ( unsigned int )atoi( temp.c_str() );
            }
            else if( temp == "QUALITYCUTS" )
            {
                is_stream >> temp;
                f_ntubes_min = atoi( temp.c_str() );
                is_stream >> temp;
                f_size_min = atof( temp.c_str() );
                is_stream >> temp;
                f_length_min = atof( temp.c_str() );
                is_stream >> temp;
                f_loss_max = atof( temp.c_str() );
            }
            else if( temp == "ZENITHANGLE" )
            {
                while( !is_stream.eof() )
                {
                    is_stream >> temp;
                    f_ze.push_back( temp );
                }
            }
            else if( temp == "WOBBLEOFFSET" )
            {
                while( !is_stream.eof() )
                {
                    is_stream >> temp;
                    f_woff.push_back( temp );
                }
            }
            else if( temp == "NOISELEVEL" )
            {
                while( !is_stream.eof() )
                {
                    is_stream >> temp;
                    f_noise.push_back( temp );
                }
            }
        }
    }
    cout << endl;
    cout << "number of telescopes: " << fNTel << ", scaling parameters: length:" <<  fLengthScaleParameter << ", width: " << fWidthScaleParameter << endl;
    cout << "quality cut: ntubes > " << f_ntubes_min << ", size > " << f_size_min << ", length > " << f_length_min << ", loss < " << f_loss_max << endl;
    cout << "read in " << f_ze.size() << " zenith bins, " << f_woff.size() << " wobble offset bins and " << f_noise.size() << " noise level bins" << endl;
    cout << "(total number of events to read from each MC data file: " << fNTotEvents << ")" << endl;
    cout << endl;
    
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{
    cout << endl;
    cout << "makeDISPTables " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "=======================================================" << endl;
    
    if( argc != 3 && argc != 4 )
    {
        cout << "makeDISPTables " << "<runparameter file> <output root file> [options]" << endl;
        cout << endl;
        cout << "look for makeDISPTables.input.dat as an example for a runparameter file" << endl;
        cout << endl;
        cout << "[options]:" << endl << endl;
        cout << "    -printpaths     print input simulation data paths and quit" << endl;
        cout << endl;
        cout << "(ARGC " << argc << ")" << endl;
        exit( 0 );
    }
    
    // run timing
    TStopwatch fStopWatch;
    fStopWatch.Start();
    
    if( !readInputParameter( argv[1] ) )
    {
        exit( 0 );
    }
    
    // output root file
    string iTableFile = argv[2];
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // disp table
    VDispTable* fDisp = new VDispTable( fNTel, iTableFile );
    fDisp->setQualityCuts( f_ntubes_min, f_size_min, f_length_min, f_loss_max );
    fDisp->setWidthLengthScalingParameters( fWidthScaleParameter, fLengthScaleParameter );
    
    // add four azimuth bins
    fDisp->addAzBin( 135., -135. );
    fDisp->addAzBin( -135., -45. );
    fDisp->addAzBin( -45., 45. );
    fDisp->addAzBin( 45., 135. );
    
    fDisp->setDataVectors( f_ze, f_woff, f_noise );
    
    // loop over all files
    for( unsigned int i = 0; i < f_ze.size(); i++ )
    {
        for( unsigned int j = 0; j < f_woff.size(); j++ )
        {
            for( unsigned int k = 0; k < f_noise.size(); k++ )
            {
                string iName = fInputDir;
                if( iName.find( "__ZENITHANGLE__" ) != string::npos )
                {
                    iName.replace( iName.find( "__ZENITHANGLE__" ), 15, f_ze[i] );
                }
                if( iName.find( "__WOBBLEOFFSET__" ) != string::npos )
                {
                    iName.replace( iName.find( "__WOBBLEOFFSET__" ), 16, f_woff[j] );
                }
                if( iName.find( "__NOISELEVEL__" ) != string::npos )
                {
                    iName.replace( iName.find( "__NOISELEVEL__" ), 14, f_noise[k] );
                }
                // hardwired file name (!)
                iName += "100" + f_ze[i] + ".root";
                
                float iNoise = fDisp->fillTable( iName, atof( f_ze[i].c_str() ), atof( f_woff[j].c_str() ), fNTotEvents );
                fDisp->setNoise( k, iNoise );
            }
        }
    }
    fDisp->terminate();
    
    fStopWatch.Stop();
    fStopWatch.Print();
    
    cout << "exiting...(success)" << endl;
    
    return 0;
}
