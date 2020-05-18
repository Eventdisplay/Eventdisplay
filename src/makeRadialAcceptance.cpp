/*! \file  makeRadialAcceptance
 *  \brief calculate radial acceptance from data
 *
 *   use off events which pass all gamma/hadron separation cuts
 *
 */

#include "CData.h"
#include "TFile.h"

#include "VEvndispRunParameter.h"
#include "VRadialAcceptance.h"
#include "VGammaHadronCuts.h"
#include "VAnaSumRunParameter.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>

using namespace std;

int parseOptions( int argc, char* argv[] );
string listfilename = "";
string cutfilename = "";
string simpleListFileName = "";
string outfile = "acceptance.root";
unsigned int ntel = 4;                   // this shouldn't be changed unless you really unterstand why
string datadir = "../eventdisplay/output";
int entries = -1;
string histdir = "" ;
struct stat sb ;
string fInstrumentEpoch = "NOT_SET";
double fMaxDistanceAllowed = 2.0;
string teltoanastring = "1234";
vector<unsigned int> teltoana;
// IRF production IO (shorter)
bool   production_shortIO = false;
// exclude regions with stars from radial acceptance
// calculations (must be given in an anasum-style
// runparameter file
string exclusionregionfile = "";
// target region excluded
bool fRemoveTargetRegionFromAcceptanceFilling = true;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
    // print program version only
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
    
    ///////////////////////////////////////////
    // use VAnaSumRunParameter to read run lists and list of
    // exclusion regions
    VAnaSumRunParameter* fRunPara = new VAnaSumRunParameter();
    
    ///////////////////////////////////////////
    // print general stuff and version numbers
    cout << endl;
    cout << "makeRadialAcceptance (" << fRunPara->getEVNDISP_VERSION() << ")" << endl << endl;
    cout << "determine radial acceptance from off events (after cuts)" << endl;
    cout << endl;
    
    ///////////////////////////////////////////
    // read options from command line
    parseOptions( argc, argv );
    
    ///////////////////////////////////////////
    // read exclusion regions
    // (note that exclusion regions per run are added later)
    if( exclusionregionfile.size() > 0 )
    {
        cout << "Reading exclusion regions from " << exclusionregionfile << endl;
        fRunPara->readRunParameter( exclusionregionfile, true );
    }
    else
    {
       cout << "No exclusion regions defined" << endl;
    }
    if( fRemoveTargetRegionFromAcceptanceFilling )
    {
       cout << "Removing target region from radial acceptance filling" << endl;
    }
    
    ///////////////////////////////////////////
    // find telescopes to analyse
    for( unsigned int i = 0; i < teltoanastring.length(); i++ )
    {
        char tel = teltoanastring.at( i );
        teltoana.push_back( atoi( &tel ) - 1 );
    }
    
    ///////////////////////////////////////////
    // read file list from run list file
    if( listfilename.size() > 0 )
    {
        fRunPara->loadLongFileList( listfilename, true );
        if( fRunPara->getRunListVersion() < 4 )
        {
            cout << "require run list >= 4" << endl;
            cout << "...exiting" << endl;
            exit( EXIT_FAILURE );
        }
        fRunPara->getEventdisplayRunParameter( datadir );
    }
    else if( simpleListFileName.size() > 0 )
    {
        fRunPara->loadSimpleFileList( simpleListFileName );
    }
    else
    {
        cout << "error reading run list" << endl;
        cout << "exiting..";
        exit( EXIT_FAILURE );
    }
    
    ///////////////////////////////////////////
    // initialize gamma/hadron cuts
    VGammaHadronCuts* fCuts = 0;
    if( cutfilename.size() > 0 && cutfilename.find( ".root" ) != string::npos )
    {
        // read cuts from effective area file
        fCuts = fRunPara->getGammaHadronCuts( cutfilename );
        if( !fCuts )
        {
            cout << "error reading gamma/hadron ctus from " << cutfilename << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        // set reconstruction type to default if not set
        if( fCuts->fReconstructionType == NOT_SET )
        {
            fRunPara->fReconstructionType = GEO;
            fCuts->setReconstructionType( fRunPara->fReconstructionType );
        }
    }
    else if( cutfilename.size() > 0 )
    {
        // read gamma/hadron cuts from cut file
        // (ascii file)
        fCuts = new VGammaHadronCuts();
        fCuts->setInstrumentEpoch( fInstrumentEpoch );
        fCuts->setTelToAnalyze( teltoana );
        fCuts->setNTel( ntel );
        // set reconstruction type (e.g. GEO, DISP, FROGS, ...
        fCuts->setReconstructionType( fRunPara->fReconstructionType );
        fCuts->readCuts( cutfilename );
    }
    else
    {
        cout << "error: no gamma/hadron cut file given" << endl;
        cout << "(command line option -c)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( !fCuts )
    {
        cout << "error reading cut file: " << cutfilename << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "Max distance of events to camera centre: " << fMaxDistanceAllowed << " [deg]" << endl;
    cout << "Instrument epoch is " << fInstrumentEpoch << endl;
    cout << "Telescopes to analyse: " << teltoanastring << endl;
    if( production_shortIO )
    {
        cout << "Write histograms required for production only" << endl;
    }
    else
    {
        cout << "Write long list of histograms" << endl;
    }
    cout << "total number of files to read: " << fRunPara->fRunList.size() << endl;
    
    char ifile[1800];
    
    //////////////////////////////////////////////////////////////
    // create output file
    TFile* fo = new TFile( outfile.c_str(), "RECREATE" );
    if( fo->IsZombie() )
    {
        cout << "makeRadialAcceptances: error opening output file " << outfile << endl;
        exit( EXIT_FAILURE );
    }
    cout << endl << "open output file for acceptance curves: " << fo->GetName() << endl;
    TDirectory* facc_dir = ( TDirectory* )fo;
    
    // create acceptance object
    VRadialAcceptance* facc = new VRadialAcceptance( fCuts, fRunPara, fMaxDistanceAllowed );
    facc->setProductionIO( production_shortIO );
    
    // set facc to write extra histograms if necessary
    if( histdir.size() > 0 )
    {
        if( stat( histdir.c_str(), &sb ) == 0 && S_ISDIR( sb.st_mode ) ) // then directory 'histdir' exists
        {
            //facc->SetExtraHistogramMode( 1 ) ;
            facc->SetExtraHistogramDirectory( histdir ) ;
        }
        else
        {
            cout << "Error, directory specified by makeRadialAcceptance -w option '";
            cout << histdir << "' does not exist, exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    // az dependent radial acceptance curves
    vector< VRadialAcceptance* > facc_az;
    vector< string > fDirName;
    vector< string > fDirTitle;
    vector< TDirectory* > facc_az_dir;
    vector< double > iAz_min;
    vector< double > iAz_max;
    
    if( !production_shortIO )
    {
        iAz_min.push_back( 337.5 );
        iAz_max.push_back( 22.5 );
        for( unsigned int i = 1; i < 8; i++ )
        {
            iAz_min.push_back( -22.5 + 45. * ( double )i );
            iAz_max.push_back( 22.5 + 45. * ( double )i );
        }
    }
    for( unsigned int i = 0; i < iAz_min.size(); i++ )
    {
        if( facc_dir )
        {
            facc_dir->cd();
        }
        sprintf( ifile, "az_%u", i );
        fDirName.push_back( ifile );
        sprintf( ifile, "AZ dependend radial acceptance, %.2f < az < %.2f", iAz_min[i], iAz_max[i] );
        fDirTitle.push_back( ifile );
        facc_az_dir.push_back( fo->mkdir( fDirName.back().c_str(), fDirTitle.back().c_str() ) );
        if( facc_az_dir.back() )
        {
            facc_az_dir.back()->cd();
            facc_az.push_back( new VRadialAcceptance( fCuts, fRunPara ) );
            facc_az.back()->setAzCut( iAz_min[i], iAz_max[i] );
            cout << "initializing azimuth dependend radial acceptance class for ";
            cout << iAz_min[i] << " < az <= " <<   iAz_max[i] << endl;
        }
    }
    
    //////////////////////////////////////////////////////////////
    // main loop over all runs
    //   - fill radial acceptances per run
    //   - fill average radial acceptance
    for( unsigned int i = 0; i < fRunPara->fRunList.size(); i++ )
    {
        sprintf( ifile, "%s/%d.mscw.root", datadir.c_str(), fRunPara->fRunList[i].fRunOff );
        cout << "now chaining " << ifile;
        cout << " (wobble offset " << -1.*fRunPara->fRunList[i].fWobbleNorth;
        cout << ", " << fRunPara->fRunList[i].fWobbleWest << ")" << endl;
        
        // open data (mscw_energy) file
        TFile fDataFile( ifile );
        if( fDataFile.IsZombie() )
        {
            cout << "error: file not found, " << ifile << endl;
            exit( EXIT_FAILURE );
        }
        // get data tree
        TTree* c = ( TTree* )fDataFile.Get( "data" );
        CData* d = new CData( c, false, true );
        if( !d )
        {
            cout << "makeRadialAcceptance: no data tree defined: run " << fRunPara->fRunList[i].fRunOff << endl;
            exit( EXIT_FAILURE );
        }
        // set reconstruction type (e.g. GEO, DISP, ...)
        d->setReconstructionType( fCuts->fReconstructionType );
        
        /////////////////////////////////////
        // read run parameters from mscw file
        VEvndispRunParameter* iParV2 = ( VEvndispRunParameter* )fDataFile.Get( "runparameterV2" );
        double iMJD = 0.;
        if( iParV2 )
        {
            ostringstream iTel_temp;
            for( unsigned int i = 0; i < iParV2->fTelToAnalyze.size(); i++ )
            {
                iTel_temp << iParV2->fTelToAnalyze[i] + 1;
            }
            if( teltoanastring != iTel_temp.str() )
            {
                cout << endl;
                cout << "error: Requested telescopes " << teltoanastring << "(" << teltoanastring.size() << ")";
                cout << " do not equal telescopes in run " << ifile << endl;
                cout << "(found: " << iTel_temp.str() << ")" << endl;
                exit( EXIT_FAILURE );
            }
            
            fRunPara->fRunList[i].fTargetRAJ2000       = iParV2->fTargetRA;
            fRunPara->fRunList[i].fTargetDecJ2000      = iParV2->fTargetDec;
            fRunPara->fRunList[i].fWobbleNorth         = iParV2->fWobbleNorth;
            fRunPara->fRunList[i].fWobbleWest          = -1.*iParV2->fWobbleEast;
            iMJD = 0.5 * ( iParV2->fDBDataStoppTimeMJD + iParV2->fDBDataStartTimeMJD );
            if( fRunPara->fRunList[i].fTargetDecJ2000 < -89.99 )
            {
                cout << "ERROR in makeRadialAcceptance: invalid target coordinates " << endl;
                cout << "\t run " << fRunPara->fRunList[i].fRunOn << "\t";
                cout << fRunPara->fRunList[i].fTargetDecJ2000 << "\t" << fRunPara->fRunList[i].fTargetShiftDecJ2000 << endl;
                exit( EXIT_FAILURE );
            }
        }
        else
        {
            cout << "error reading run parameters from mscw file" << endl;
            exit( EXIT_FAILURE );
        }
        
        // initialize gamma/hadron cuts
        fCuts->initializeCuts( fRunPara->fRunList[i].fRunOff, datadir );
        fCuts->setDataTree( d );
        // data trees and cuts
        int nentries = d->fChain->GetEntries();
        if( entries > 0 )
        {
            nentries = entries;
        }
        cout << "filling acceptance curves with " << nentries << " events (before gamma/hadron cuts)" << endl;
        
        ////////////////////////////////////////////
        // calculate coordinates of camera centre:
        
        double raJ2000_deg = 0.;
        double decJ2000_deg = 0.;
        VSkyCoordinatesUtilities::getCameraCentreCoordinates_J2000( iMJD,
                fRunPara->fRunList[i].fTargetRAJ2000,
                fRunPara->fRunList[i].fTargetDecJ2000,
                fRunPara->fRunList[i].fWobbleNorth,
                -1.*fRunPara->fRunList[i].fWobbleWest,    // require wobble EAST
                raJ2000_deg, decJ2000_deg );
                
        cout << "Camera centre at (ra,dec) J2000 : ( " << raJ2000_deg << ", " << decJ2000_deg << ")" << endl;
        ////////////////////////////////////////////
        // initialize exclusion regions

        // exclude source direction
        if( fRunPara->getExclusionRegions() )
        {
              fRunPara->getExclusionRegions()->addExclusionRegion( -1., -1.,
                                       fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                                       0.25, 0.25, 0., "on source" );
        }
        if( exclusionregionfile.size() > 0 )
        {
            cout << "Exclusion regions: " << endl;
            double iMax = fMaxDistanceAllowed * 1.3;
            VStarCatalogue iStarCatalogue;
            iStarCatalogue.init( iMJD, fRunPara->getStarCatalogue() );
            iStarCatalogue.setFOV( raJ2000_deg, decJ2000_deg, iMax, iMax, true );
            // (note that sky map and camera centre is the same in this context)
            fRunPara->initializeExclusionRegions( i, &iStarCatalogue,
                                                  raJ2000_deg, decJ2000_deg,
                                                  raJ2000_deg, decJ2000_deg );
            facc->setRegionToExcludeAcceptance( fRunPara->getExclusionRegions( i ) );
            for( unsigned int a = 0; a < facc_az.size(); a++ )
            {
                 if( facc_az[a] )
                 {
                      facc_az[a]->setRegionToExcludeAcceptance( fRunPara->getExclusionRegions( i ) );
                 }
            }
        }
        
        if( fCuts )
        {
            fCuts->printCutSummary();
        }
        
        int neventStats = 0;
        int i_entries_after_cuts = 0;
        
        double x_rotJ2000 = 0.;
        double y_rotJ2000 = 0.;
        
        ////////////////////////////////////////////////////////////////////
        // loop over all entries in data trees and fill acceptance curves
        for( int n = 0; n < nentries; n++ )
        {
            d->GetEntry( n );
            
            // printout for MC
            if( n == 0 and d->isMC() )
            {
                cout << "\t (analysing MC data)" << endl;
            }
            
            // convert de-rotated coordinates to J2000
            // (exclusion regions are in J2000)
            x_rotJ2000 = d->getXoff_derot();
            y_rotJ2000 = d->getYoff_derot();
            VSkyCoordinatesUtilities::convert_derotatedCoordinates_to_J2000(
                d->MJD, d->Time,
                d->ArrayPointing_Azimuth, d->ArrayPointing_Elevation,
                x_rotJ2000, y_rotJ2000 );
                
            // check if event is inside an exclusion region
            if( exclusionregionfile.size() > 0 && facc->isExcludedfromBackground( x_rotJ2000, y_rotJ2000 ) )
            {
                continue;
            }
            
            ////////////////////////////////
            // fill acceptances
            
            neventStats = facc->fillAcceptanceFromData( d, n, x_rotJ2000, y_rotJ2000 );
            
            for( unsigned int a = 0; a < facc_az.size(); a++ )
            {
                if( facc_az[a] )
                {
                    facc_az[a]->fillAcceptanceFromData( d, n, x_rotJ2000, y_rotJ2000 );
                }
            }
            
            if( neventStats < 0 )
            {
                break;
            }
            i_entries_after_cuts += neventStats;
        }
        cout << "total number of entries after cuts: " << i_entries_after_cuts << endl;
        cout << endl << endl;
        
        fDataFile.Close();
        facc->correctRadialAcceptancesForExclusionRegions( facc_dir, fRunPara->fRunList[i].fRunOff );
        for( unsigned int a = 0; a < facc_az.size(); a++ )
        {
             if( facc_az[a] )
             {
                  if( facc_az[a] )
                  {
                      facc_az[a]->correctRadialAcceptancesForExclusionRegions( facc_az_dir[a],
                                                                               fRunPara->fRunList[i].fRunOff );
                  }
             }
        }
    }
    
    /////////////////////////////////////////
    // write acceptance files to disk
    facc->calculateAverageRadialAcceptanceCurveFromRuns( facc_dir );
    facc->calculate2DBinNormalizationConstant() ;
    facc->terminate( facc_dir );
    for( unsigned int a = 0; a < facc_az_dir.size(); a++ )
    {
        if( facc_az[a] )
        {
            facc_az[a]->calculateAverageRadialAcceptanceCurveFromRuns( facc_az_dir[a] );
            facc_az[a]->terminate( facc_az_dir[a] );
        }
    }
    
    fo->Close();
    cout << "closing radial acceptance file: " << fo->GetName() << endl;
    
    cout << "exiting.." << endl;
}


/*!

    read in command line parameters

*/
int parseOptions( int argc, char* argv[] )
{
    while( 1 )
    {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"runlist", required_argument, 0, 'l'},
            {"srunlist", required_argument, 0, 's'},
            {"cutfile", required_argument, 0, 'c'},
            {"instrumentepoch", required_argument, 0, 'i'},
            {"maxdistance", required_argument, 0, 'm'},
            {"outfile", required_argument, 0, 'o'},
            {"entries", required_argument, 0, 'n'},
            {"datadir", required_argument, 0, 'd'},
            {"writehists", optional_argument, 0, 'w'},
            {"teltoana", required_argument, 0, 't'},
            {"productionIO", required_argument, 0, 'p'},
            {"--remove_target", required_argument, 0, 'r'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        int c = getopt_long( argc, argv, "ht:s:f:p:l:e:m:r:o:i:d:n:c:w:t:", long_options, &option_index );
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
                cout << endl;
                cout << "Options are:" << endl;
                cout << "-l --runlist [anasum-style run list file name, runlist on/off like]" << endl;
                cout << "-s --srunlist [simple run list file name]" << endl;
                cout << "-c --cutfile [cut file name]" << endl;
                cout << "-i --instrumentepoch [instrument epoch (e.g. V6)" << endl;
                cout << "-d --datadir [directory for input mscw root files]" << endl;
                cout << "-o --outfile [output ROOT file name]" << endl;
                cout << "-e --entries [number of entries]" << endl;
                cout << "-m --maxdist [max distance from camera centre (deg)]" << endl;
                cout << "-w --writehists [directory]" << endl ;
                cout << "-r --remove_target region from acceptance filling [default=true]" << endl;
                cout << "-t --teltoana <telescopes>" << endl;
                cout << "-p --productionIO [0/1] " << endl;
                cout << "-f --exclusionregionfile [file with exclusion regions]" << endl;
                cout << endl;
                exit( EXIT_SUCCESS );
                break;
            case 'd':
                cout << "Directory for input Files is " << optarg << endl;
                datadir = optarg;
                break;
            case 'o':
                cout << "Output File Name is " << optarg << endl;
                outfile = optarg;
                break;
            case 'l':
                cout << "Run List File Name is " << optarg << endl;
                listfilename = optarg;
                break;
            case 's':
                cout << "Simple List File Name is " << optarg << endl;
                simpleListFileName = optarg;
                break;
            case 'c':
                cutfilename = optarg;
                cout << "Cuts are taken from " << cutfilename << endl;
                if( cutfilename == "IGNOREEFFECTIVEAREA" )
                {
                     cout << "error: cannot read cuts (IGNOREEFFECTIVEAREA given)" << endl;
                     cout << "exiting..." << endl;
                     exit( EXIT_FAILURE );
                }  
                break;
            case 'm':
                fMaxDistanceAllowed = atof( optarg );
                break;
            case 'i':
                fInstrumentEpoch = optarg;
                break;
            case 'e':
                entries = ( int )atoi( optarg );
                break;
            case 'w':
                histdir = optarg;
                cout << "Extra histograms will be written to " << histdir << endl;
                break;
            case 't':
                teltoanastring = optarg;
                break;
            case 'p':
                production_shortIO = ( int )atoi( optarg );
            case 'f':
                exclusionregionfile = optarg;
                break;
            case 'r':
                fRemoveTargetRegionFromAcceptanceFilling = (bool)atoi(optarg);
                break;
            case '?':
                break;
            default:
                exit( EXIT_SUCCESS );
        }
    }
    return optind;
}
