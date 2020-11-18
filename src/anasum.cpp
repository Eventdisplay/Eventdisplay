/*! \file anasum.cpp
    \brief main program to create an analysis summary (VERITAS data analysis chain)

*/

#include "VAnaSum.h"
#include "VGlobalRunParameter.h"

#include <getopt.h>
#include <iostream>
#include <string>

using namespace std;

int parseOptions( int argc, char* argv[] );

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters read in from command line
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// input run list
string listfilename = "";
// output file with all anasum results
string outfile = "output.ansum.root";
// run types:
//   0:  sequentiell analysis of a list file
//   1:  combine a list of runs and do a combined analysis
unsigned int runType = 0;
// location of data files (might be mscw file (run type 0) or anasum result file (run type 1 )
string datadir = "";
// run parameter file
string fRunParameterfile = "ANASUM.runparameter";
// for usage of random generators: see VStereoMaps.cpp
int fRandomSeed = 17;
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

   count number of command line arguments

*/
bool testCommandlineArguments()
{
    // require a runlist file
    if( listfilename.size() < 1 )
    {
        cout << "error: missing required command line argument --runlist (-l)" << endl;
        return false;
    }
    // require data directory
    if( datadir.size() < 1 )
    {
        cout << "error: missing required command line argument --datadir (-d)" << endl;
        return false;
    }
    return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{

    // print version only
    if( argc == 2 )
    {
        string fCommandLine = argv[1];
        if( fCommandLine == "-v" || fCommandLine == "--version" )
        {
            VGlobalRunParameter fRunPara;
            cout << fRunPara.getEVNDISP_VERSION() << endl;
            exit( EXIT_SUCCESS );
        }
    }
    
    cout << endl << "VERITAS Analysis Summary (University of Delaware & DESY) ";
    cout << " (version " << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout <<         "==========================================================================" << endl;
    cout << endl;
    
    parseOptions( argc, argv );
    
    if( !testCommandlineArguments() )
    {
        exit( EXIT_FAILURE );
    }
    
    // initialize analysis
    VAnaSum* anasum = new VAnaSum( datadir );
    anasum->initialize( listfilename, runType, outfile, fRandomSeed, fRunParameterfile );
    cout << endl;
    
    // stereo analysis (default)
    anasum->doStereoAnalysis();
    
    // clean up and write results to disk
    anasum->terminate();
    
    cout << endl << "analysis results written to " << outfile << endl;
    
    return 0;
}

/*
 * read command line options
 */
int parseOptions( int argc, char* argv[] )
{
    while( 1 )
    {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"runlist", required_argument, 0, 'l'},
            {"outfile", required_argument, 0, 'o'},
            {"datadir", required_argument, 0, 'd'},
            {"randomseed", required_argument, 0, 'r'},
            {"runType", required_argument, 0, 'i'},
            {"parameterfile",  required_argument, 0, 'f'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        int c = getopt_long( argc, argv, "h:l:k:m:o:d:s:r:i:u:f:g", long_options, &option_index );
        if( optopt != 0 )
        {
            cout << "error: unknown option" << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        if( argc == 1 )
        {
            c = 'h';
        }
        if( c == -1 )
        {
            break;
        }
        switch( c )
        {
            case 0:
                if( long_options[option_index].flag != 0 )
                {
                    break;
                }
                printf( "option %s", long_options[option_index].name );
                if( optarg )
                {
                    printf( " with arg %s", optarg );
                }
                printf( "\n" );
                break;
            case 'h':
                if( gSystem->Getenv( "EVNDISPSYS" ) )
                {
                    int i_s = system( "cat $EVNDISPSYS/README/README.ANASUM" );
                    if( i_s == -1 )
                    {
                         cout << "error, README.ANASUM not found" << endl;
                    }
                }
                else
                {
                    cout << " no help find (environmental variable EVNDISPSYS not set)" << endl;
                }
                exit( EXIT_FAILURE );
                break;
            case 'd':
                datadir = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'i':
                runType = ( unsigned int )atoi( optarg );
                break;
            case 'l':
                listfilename = optarg;
                break;
            case 'r':
                fRandomSeed = atoi( optarg );
                break;
            case 'f':
                fRunParameterfile = optarg;
                break;
            case '?':
                break;
            default:
                abort();
        }
    }
    return optind;
}
