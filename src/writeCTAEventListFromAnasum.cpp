#include "writeCTAEventListFromAnasum.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////
//
//  :: NAVIGATING THIS FILE ::
//
// search for the following phrases to jump to the relevant sections:
//
//  -evlist  : gamma ray event list table
//  -psf     : point spread function table
//  -effarea : effective area table
//  -eres    : energy resolution table
//  -bkgnd   : background table
//
//////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
    // root option to ignore "no dictionary for class X" warning
    gErrorIgnoreLevel = kError ;
    
    // setup vars for input args
    bool helptext              = false ;
    bool overwriteFlag         = false ;
    bool maxEventsFlag         = false ;
    bool printEventsFlag       = false ;
    bool writedebugroothists   = false ;
    bool printeffareafilenames = false ;
    bool skipeffarea           = false ;
    bool skippsf               = false ;
    unsigned int maxEvents     = -1    ;
    double chunkLength         = 0.0   ; // units=seconds
    char anasumInputFileChar[200] = "" ;
    char outputDirectoryChar[200] = "" ;
    string anasumInputFile ;
    string outputDirectory ;
    int psfTable = -1 ;
    
    char argflag[200] ;
    for( int iarg = 1 ; iarg < argc ; iarg++ )
    {
        sprintf( argflag, "%s", argv[iarg] ) ;
        if( strcmp( "-i", argflag ) == 0 || strcmp( "--inputfile",   argflag ) == 0 )
        {
            //cout << "-i: read in flag '" << argflag << "' with 2nd argument '" << argv[iarg+1] << "'" << endl;
            checkArgBounds( argflag, argc, iarg ) ;
            sprintf( anasumInputFileChar, "%s", argv[iarg + 1] ) ;
            anasumInputFile = anasumInputFileChar ;
            iarg++ ;
        }
        else if( strcmp( "-o", argflag ) == 0 || strcmp( "--outputdir",   argflag ) == 0 )
        {
            //cout << "-o: read in flag '" << argflag << "' with 2nd argument '" << argv[iarg+1]  << "'"<< endl;
            checkArgBounds( argflag, argc, iarg ) ;
            sprintf( outputDirectoryChar, "%s", argv[iarg + 1] ) ;
            outputDirectory = outputDirectoryChar ;
            iarg++ ;
        }
        else if( strcmp( "-n", argflag ) == 0 || strcmp( "--nevents",     argflag ) == 0 )
        {
            //cout << "-n: read in flag '" << argflag << "' with 2nd argument '" << argv[iarg+1]  << "'"<< endl;
            checkArgBounds( argflag, argc, iarg ) ;
            maxEventsFlag = true ;
            maxEvents     = atoi( argv[iarg + 1] ) ;
            iarg++ ;
        }
        else if( strcmp( "-p", argflag ) == 0 || strcmp( "--psf", argflag ) == 0 )
        {
            checkArgBounds( argflag, argc, iarg ) ;
            if( strcmp( "king" , argv[iarg + 1] ) == 0 )
            {
                psfTable = PSFTABLE_BASIC ;
                iarg++ ;
            }
            else if( strcmp( "gauss", argv[iarg + 1] ) == 0 )
            {
                psfTable = PSFTABLE_GAUSS ;
                iarg++ ;
            }
            else
            {
                printf( "Unrecognized argument '%s' for flag '%s', must be either 'gauss' or 'king'\n", argv[iarg + 1], argv[iarg] ) ;
                exit( 1 );
            }
        }
        else if( strcmp( "-c", argflag ) == 0 || strcmp( "--chunklength", argflag ) == 0 )
        {
            //cout << "-c: read in flag '" << argflag << "' with 2nd argument '" << argv[iarg+1]  << "'"<< endl;
            checkArgBounds( argflag, argc, iarg ) ;
            chunkLength = parseTimeArg( argv[iarg + 1] ) ;
            if( chunkLength < 3.0 )
            {
                printf( "Error, -c option must be greater than 3 seconds, exiting...\n" );
                exit( 1 ) ;
            }
            iarg++;
        }
        else if( strcmp( "-h", argflag ) == 0 || strcmp( "--help",        argflag ) == 0 )
        {
            //cout << "-h: read in flag '" << argflag << "'" << endl;
            helptext = true ;
        }
        else if( strcmp( "-v", argflag ) == 0 || strcmp( "--verbose", argflag ) == 0 )
        {
            //cout << "-p: read in flag '" << argflag << "'" << endl;
            printEventsFlag = true ;
        }
        else if( strcmp( "-f", argflag ) == 0 || strcmp( "--forceoverwrite", argflag ) == 0 )
        {
            //cout << "-f: read in flag '" << argflag << "'" << endl;
            overwriteFlag = true ;
        }
        else if( strcmp( "-r", argflag ) == 0 || strcmp( "--responserootfile", argflag ) == 0 )
        {
            writedebugroothists = true ;
        }
        else if( strcmp( "-e", argflag ) == 0 || strcmp( "--effectiveareafiles", argflag ) == 0 )
        {
            printeffareafilenames = true ;
        }
        else if( strcmp( "--skipeffarea", argflag ) == 0 )
        {
            skipeffarea = true ;
        }
        else if( strcmp( "--skippsf", argflag ) == 0 )
        {
            skippsf = true ;
        }
        else
        {
            printf( "Error, unrecognized flag '%s', exiting...\n", argflag ) ;
            exit( 1 ) ;
        }
    }
    
    // figure out example command string
    char examplecmd[150] = "" ;
    sprintf( examplecmd, "$ %s -i 45538.anasum.root -o myoutputdirectory/", argv[0] ) ;
    
    // print help text if needed
    if( argc < 7 || helptext )
    {
        cout << endl << argv[0] << endl;
        cout << "   For converting eventdisplay anasum files to the 'ctools fits' format." << endl;
        cout << "   Must specify at least the input filename, output directory, and effective area filename." << endl << endl ;
        printf( "   example: \n" );
        printf( "      %s\n\n", examplecmd ) ;
        cout << "   Input Arguments are:"                                  << endl ;
        cout << "      -i, --inputfile   <fname> : input anasum root file to convert ***" << endl ;
        cout << "      -o, --outputdir   <fname> : output directory to save fits files to"                      << endl ;
        cout << "   Optional:"                                                                                  << endl ;
        cout << "      -h, --help                : prints help text"                                            << endl ;
        cout << "      -c, --chunklength <time>  : break up each gti into <tim> long chunks."                   << endl ;
        cout << "                                  <time> should have the format ####[smh], like: 14s, 22m, or" << endl ;
        cout << "                                  3h, etc.  default: only break up the run by their gti's"     << endl ;
        cout << "      -p, --psf <type>          : write psf table with specific format, default to gauss"      << endl ;
        cout << "                                  <type>=gauss : fit psf with gauss function"                  << endl ;
        cout << "                                  <type>=king  : fit psf with king function (warning, experimental, won't have proper king function in most of parameter space due to poor fitting)" << endl;
        cout << "      -n, --nevents <int>       : only read at most <int> events total, for debugging."        << endl ;
        cout << "      -f, --forceoverwrite      : will overwrite existing output fits file."                   << endl ;
        cout << "      -r, --responserootfile    : will save all used IRF histograms to <outputdir>/<runnumber>.debughists.root" << endl ;
        cout << "                                  otherwise error out if file already exists."                 << endl ;
        cout << "      --skipeffarea             : will skip writing effective area table"                      << endl ;
        cout << "      --skippsf                 : will skip writing psf table"                                 << endl ;
        cout << "      -e, --effectiveareafiles  : will print list of needed effective area files, then exit."  << endl ;
        cout << "                                  will not write any files at all."                            << endl ;
        cout << endl;
        cout << "   *** Note: input anasum file must have been created with the '* WRITEEVENTTREEFORCTOOLS 1'" << endl;
        cout << "        option in its ANASUM.runparameter file, which is not on by default.  The output anasum" << endl;
        cout << "        file will then have the extra needed information in it." << endl;
        cout << endl;
        return 1 ;
    }
    
    // print info about which arguments were read
    cout << "input anasum file     : " << anasumInputFile.substr( anasumInputFile.find_last_of( "\\/" ) + 1, -1 ) << endl;
    cout << "output directory      : " << outputDirectory << endl;
    cout << "overwrite output files: " ;
    if( overwriteFlag )
    {
        cout << "yes" << endl ;
    }
    else
    {
        cout << "no" << endl ;
    }
    if( maxEventsFlag )
    {
        cout << "events to read in     : " << maxEvents  << endl;
    }
    if( chunkLength != 0.0 )
    {
        cout << "chunk length          : " << chunkLength << endl;
    }
    if( psfTable == -1 )
    {
        psfTable = PSFTABLE_GAUSS ;
    }
    if( psfTable == PSFTABLE_GAUSS )
    {
        cout << "psf table             : gauss" << endl;
    }
    else if( psfTable == PSFTABLE_BASIC )
    {
        cout << "psf table             : king " << endl;
    }
    
    // check input anasum file
    ifstream acheck( anasumInputFile.c_str() ) ;
    if( strcmp( anasumInputFile.c_str(), "" ) == 0 )
    {
        cout << "Need to specify input anasum file with -i" << endl;
        printf( "   ex: %s\n", examplecmd ) ;
        return 1 ;
    }
    else if( ! acheck.good() )
    {
        cout << "Error, could not read input anasum file '" << anasumInputFile << "', exiting." << endl;
        return 1 ;
    }
    acheck.close();
    
    // check anasum file is a root file
    TFile* aroot = new TFile( anasumInputFile.c_str() ) ;
    if( aroot->IsZombie() || aroot->TestBit( TFile::kRecovered ) )
    {
        cout << "Error, anasum file '" << anasumInputFile << "' wasn't closed correctly or is corrupt.  Exiting..." << endl;
        return 1 ;
    }
    
    // check output directory name is valid
    if( strcmp( outputDirectory.c_str(),    "" ) == 0 )
    {
        cout << "Need to specify output directory to store fits files with -o" << endl;
        printf( "   ex: %s\n", examplecmd ) ;
        return 1 ;
    }
    
    // load the input anasum root file
    TFile* anasumfile = new TFile( anasumInputFile.c_str() );
    // see root.cern.ch/drupal/content/how-read-objects-file
    // and ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch11s02.html
    // for more explaination of the below TIter and TKey algorithm
    TIter next( anasumfile->GetListOfKeys() );
    TKey* key ;
    string runRootDir = "" ;
    
    // loop over objects in the root file,
    // looking for the first run_##### directory
    while( ( key = ( TKey* )next() ) )
    {
    
        // skip all objects that aren't TDirectoryFile()'s
        if( strcmp( key->GetClassName(), "TDirectoryFile" ) != 0 )
        {
            continue ;
        }
        
        // look for TDirectoryFile() with a name like 'run_#####',
        // skip if it doesn't match
        string nam( key->GetName() ) ;
        // God help me, but there is no regex libs in the C++ standard.
        // So now this only works if the directory name contains 'run_' ,
        // rather than matches the robust '^run_\d{5,6}$' regex string.
        // I'm going to cry.
        if( nam.find( "run_" ) == std::string::npos )
        {
            continue ;
        }
        
        // get the first run_##### listed, assuming
        runRootDir = nam ;
        break ;
    }
    anasumfile->Close() ;
    
    // get crude runid from the TDirectory name
    int runid = 0 ;
    sscanf( runRootDir.c_str(), "run_%d", &runid ) ;
    
    ////////////////////////////////////////////////////////////////////////////////////
    // setup the TFile
    anasumfile = new TFile( anasumInputFile.c_str(), "READ" ) ;
    if( !anasumfile )
    {
        cout << "  Error, Unable to read file " << anasumInputFile << " , exiting..." << endl;
        return 1 ;
    }
    else
    {
        cout << "  Loaded TFile from " << anasumInputFile << "..." << endl;
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    // get VAnaSumRunParameter object from file
    char objName[200]   = "" ;
    sprintf( objName, "run_%d/stereo/VAnaSumRunParameter", runid ) ;
    VAnaSumRunParameter* anasumRunPar = ( VAnaSumRunParameter* )anasumfile->Get( objName ) ;
    if( !anasumRunPar )
    {
        cout << "  Error, Unable to load VAnaSumRunParameter object " << objName << " from file " << anasumInputFile << " , exiting..." << endl;
        return 1 ;
    }
    else
    {
        cout << "  Loaded VAnaSumRunParameter object '" << objName << "' ..." << endl;
    }
    string anasumEffFile = anasumRunPar->fRunList[0].fEffectiveAreaFile ;
    printf( "    anasum file was analyzed using effective area file '%s'\n", anasumEffFile.c_str() ) ;
    
    if( printeffareafilenames )
    {
        cout << "requires " << anasumEffFile << endl;
        exit( 0 );
    }
    
    // check which directory our eff file is stored in
    string EVNAUXDIR = get_env_var( "VERITAS_EVNDISP_AUX_DIR" ) ;
    string PWD       = get_env_var( "PWD" ) ;
    char fullEffFileName[1000] = "" ;
    // first look in $PWD (like a copied file in a batch job temp dir)
    sprintf( fullEffFileName, "%s", anasumEffFile.c_str() ) ;
    if( ! file_exists( fullEffFileName ) )
    {
        // then try looking in $VERITAS_EVNDISP_AUX_DIR/EffectiveAreas/
        sprintf( fullEffFileName, "%s/EffectiveAreas/%s", EVNAUXDIR.c_str(), anasumEffFile.c_str() ) ;
        if( ! file_exists( fullEffFileName ) )
        {
            printf( "Error, could not find needed effective area file in either $PWD='%s' or $VERITAS_EVNDISP_AUX_DIR/EffectiveAreas/='%s/EffectiveAreas/', exiting...\n", PWD.c_str(), EVNAUXDIR.c_str() ) ;
            exit( 1 );
        }
    }
    string effFile = fullEffFileName ;
    cout << "    using effective area file '" << effFile << "'" << endl;
    
    // Storage class for effective area data
    VArchiveEffArea* archiveEffArea = new VArchiveEffArea( effFile ) ;
    
    // load TTrees from anasum file
    char treeName[100] ;
    TChain* chainEventList = new TChain() ;
    TChain* chainTelConfig = new TChain() ;
    TChain* chainPointData = new TChain() ;
    
    // load event list TTree
    printf( "\nLoading needed objects from root file...\n" );
    sprintf( treeName, "%s/%s/stereo/TreeWithEventsForCtools", anasumInputFile.c_str(), runRootDir.c_str() ) ;
    chainEventList->Add( treeName ) ;
    cout << "  Loaded TTree " << treeName << "..." << endl;
    
    // load telescope config TTree 'telconfig'
    sprintf( treeName, "%s/%s/stereo/telconfig", anasumInputFile.c_str(), runRootDir.c_str() ) ;
    cout << "  Loaded TTree " << treeName << "..." << endl;
    chainTelConfig->Add( treeName ) ;
    
    // load TTree 'pointingDataReduced' containing pointing info for the run
    sprintf( treeName, "%s/%s/stereo/pointingDataReduced", anasumInputFile.c_str(), runRootDir.c_str() ) ;
    cout << "  Loaded TTree " << treeName << "..." << endl;
    chainPointData->Add( treeName ) ;
    
    // FITS records' template file
    string EVNDISPSYS = get_env_var( "EVNDISPSYS" ) ;
    char eventListTemplateChar[200] = "" ;
    sprintf( eventListTemplateChar, "%s/templates/EventList.tpl", EVNDISPSYS.c_str() ) ;
    string eventListTemplate = eventListTemplateChar ;
    
    // Test if template exists?
    if( ! file_exists( eventListTemplate ) )
    {
        cout << "  Error, template file '" << eventListTemplate << "' doesn't exist, exiting..." << endl;
        return 1 ;
    }
    
    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // get VGammaHadronCuts object from file
    sprintf( objName, "run_%d/stereo/GammaHadronCuts", runid ) ;
    VGammaHadronCuts * gammaHadronCuts = ( VGammaHadronCuts* )anasumfile->Get( objName ) ;
    if ( !gammaHadronCuts )
    {
      cout << "  Error, Unable to load VGammaHadronCuts object '" << objName << "' from file '" << anasumInputFile << "' , exiting..." << endl;
      return 1 ;
    }
    else
    {
      cout << "  Loaded VGammaHadronCuts object '" << objName << "'..." << endl;
    }
    cout << endl;
    gammaHadronCuts->printCutSummary();
    cout << endl;
    */
    
    // figure out ed version used to analyze the anasum root file
    string edverstr = command_line_regex( anasumEffFile, "grep -oP \"\\-v\\d+\\-\" | grep -oP \"\\d+\" " ) ;
    int  edver  = atoi( edverstr.c_str() ) ; // 470 ;
    printf( "    analyzed using ed version v%d\n", edver ) ;
    
    // figure out cut type that was used
    string cutss = command_line_regex( anasumEffFile, "grep -oP \"Cut-NTel\\d-\\w+-\\w+\"" ) ;
    char cuts[1000] = "" ; // aka "Cut-NTel2-PointSource-Soft" ;
    sprintf( cuts, "%s", cutss.c_str() ) ;
    printf( "    found cuts string '%s'\n", cuts ) ;
    
    // figure out sim type that was used
    //string simtypes = get_simtype_string( anasumEffFile ) ;
    string simtypes = command_line_regex( anasumEffFile, "grep -oP \"CARE_June1425|GRISU-SW6\"" ) ;
    char sims[1000] = "" ; // aka "GRISU-SW6" ;
    sprintf( sims, "%s" , simtypes.c_str() ) ;
    printf( "    found sims string '%s'\n", sims ) ;
    
    string telcombo = command_line_regex( anasumEffFile, "grep -oP \"T[1234]+\"" ) ;
    printf( "    found telcombo string '%s'\n", telcombo.c_str() ) ;
    
    // load dead time object
    sprintf( objName, "run_%d/stereo/vdeadtime", runid ) ;
    VDeadTime* deadTime = ( VDeadTime* )anasumfile->Get( objName ) ;
    if( !deadTime )
    {
        cout << "  Error, Unable to load VDeadTime object " << objName << " from file " << anasumInputFile << " , exiting..." << endl;
        return 1 ;
    }
    else
    {
        cout << "  Loaded VDeadTime object '" << objName << "' ..." << endl;
    }
    sprintf( objName, "run_%d/stereo/deadTimeHistograms", runid ) ;
    TDirectoryFile* deadTimeDir = ( TDirectoryFile* ) anasumfile->Get( objName ) ;
    deadTime->readHistograms( deadTimeDir ) ;
    deadTime->calculateDeadTime();
    deadTime->printDeadTime();
    
    // store event display version number
    char progname[100] = "" ;
    int  svnrev        = 0  ;
    sscanf( anasumRunPar->getSVN_VERSION().c_str(), "$Revision: %d $", &svnrev ) ;
    sprintf( progname, "Event Display %s, svn rev:%d", anasumRunPar->getEVNDISP_VERSION().c_str(), svnrev ) ;
    
    /////////////////////////////////////////////////////////////////////////////////////
    // get VEvndispRunParameter object from file
    sprintf( objName, "run_%d/stereo/runparameterV2", runid ) ;
    VEvndispRunParameter* evndispRunPar = ( VEvndispRunParameter* )anasumfile->Get( objName ) ;
    if( !evndispRunPar )
    {
        cout << "  Error, Unable to load VEvndispRunParameter object " << objName << " from file " << anasumInputFile << " , exiting..." << endl;
        return 1 ;
    }
    else
    {
        cout << "  Loaded VEvndispRunParameter object '" << objName << "' ..." << endl;
    }
    int  epoch  =   6 ;
    sscanf( evndispRunPar->fInstrumentEpoch.c_str(), "V%d", &epoch ) ;
    printf( "    found instrument epoch V%d\n", epoch ) ;
    int  atmo   =  evndispRunPar->fAtmosphereID ;
    printf( "    found atmo ATM%d\n", atmo ) ;
    
    // get run's pedvar from the run summary
    sprintf( objName, "total_1/stereo/tRunSummary" ) ;
    TTree* trunsummary = ( TTree* ) anasumfile->Get( objName ) ;
    if( !trunsummary )
    {
        cout << "  Error, Unable to load Run Summary TTree '" << objName << "'..." << endl;
        return 1 ;
    }
    else
    {
        cout << "  Loaded Run Summary TTree '" << objName << "'..." << endl;
    }
    
    // print basic not-obvious info about this run from the run summary
    double pedvar = getPedVarFromRunSummary( trunsummary, runid ) ;
    printf( "    Run %d has a pedestal variation of %f (units=?? )\n", runid, pedvar ) ;
    double wobble_radius = getWobbleRadiusFromRunSummary( trunsummary, runid ) ;
    printf( "    Run %d has a wobble radius of %.2f degrees\n", runid,  wobble_radius ) ;
    
    cout << endl;
    
    // run start/stop date strings, for human reading
    char dateobs[30] = "null" ;
    char timeobs[30] = "null" ;
    char dateend[30] = "null" ;
    char timeend[30] = "null" ;
    getHumanDateTimes( evndispRunPar, dateobs, timeobs, dateend, timeend ) ;
    
    // Establish Reference MJD (0.0 in MET time)
    // This is defined by:
    //   "The beginning of the UTC day of the earliest run in the database"
    //
    // Earliest Run: 30576 , 2006-05-01, UTC 08:28
    int    mjdrefi = 0   ; // day in MJD when MET is 0.0  (integer part only)
    double mjdreff = 0.0 ; // decimal part of mjdrefi
    int MJDREF_YEAR = 2006 ;
    int MJDREF_MONT =    5 ;
    int MJDREF_DAYY =    1 ;
    printf( "  Calculating Mission Elapsed Time(MET) starting at the beginning of UTC day %04d-%02d-%02d...\n", MJDREF_YEAR, MJDREF_MONT, MJDREF_DAYY ) ;
    mjdrefi = VSkyCoordinatesUtilities::getMJD( MJDREF_YEAR, MJDREF_MONT, MJDREF_DAYY ) ;
    
    // MET start, stop, and duration of run (seconds)
    double tstart   = 0.0 ; // seconds since MET reference
    double tstop    = 0.0 ; // seconds since MET reference
    double telapse  = 0.0 ;
    double met_zero = mjdrefi + mjdreff ; // MJD(days) - met_zero = MET (days)
    double startMJD = 0.0 ;
    double stopMJD  = 0.0 ;
    getRunTimes( evndispRunPar, met_zero, tstart, tstop, telapse, startMJD, stopMJD ) ;
    //printf("run:%5d tstart:%.6f tstop:%.6f dur:%.1f\n", runid, tstart, tstop, tstop-tstart);
    
    // calculate average RA, Dec, Alt, Az
    // just by a straight average
    // veritas always tracks one RA/Dec coordinate
    double ra_pnt  ;
    double dec_pnt ;
    double alt_pnt ;
    double az_pnt  ;
    
    // setup VTimeMask
    VTimeMask* vtm = 0 ;
    
    // set vtmmode
    // vtmmode: 1 = use hardcoded debugging timemask file (only for debugging!)
    //          2 = load VTimeMask object straight from the anasum file (preferred)
    int vtmmode = 2 ;
    
    // load time mask, get good time intervals
    if( vtmmode == 1 )
    {
        string debugTMfile = "/afs/ifh.de/user/n/nkelhos/scratch/tmp3V6Crab/anasum.timemask.dat" ;
        cout << endl << "  Warning, using hardcoded timemask file '" << debugTMfile << "' , and ignoring the the input anasum root file's 'vtimemask' object.  This should only be used for debugging, ***NOT FOR SCIENCE!***" << endl << endl;
        vtm  = new VTimeMask( runid, startMJD, stopMJD, debugTMfile ) ;
        sprintf( objName, "/run_%d/stereo/timeMask", runid ) ;
        anasumfile->cd( objName ) ;
        TDirectory* tdir = anasumfile->CurrentDirectory();
        vtm->readObjects( tdir ) ;
        vtm->setMask( runid, startMJD, stopMJD, debugTMfile ) ;
    }
    else if( vtmmode == 2 )
    {
        cout << "    ...from streamed object" << endl;
        sprintf( objName, "run_%d/stereo/vtimemask", runid ) ;
        vtm = ( VTimeMask* )anasumfile->Get( objName ) ;
        vtm->setMask() ;
    }
    else
    {
        cout << "  Error, vtmmode=" << vtmmode << endl;
        return 1 ;
    }
    
    // check that time mask was loaded fine
    if( ! vtm )
    {
        cout << "  Error, Unable to load VTimeMask directory " << objName << " from file '" << anasumInputFile << "' , exiting..." << endl;
        return 1 ;
    }
    else
    {
        cout << "  Loaded VTimeMask directory " << objName << "..." << endl;
    }
    
    // print our timemask
    //cout << endl;
    //cout << "displayMask()" ;
    //vtm->displayMask( cout ) ;
    //cout << endl;
    
    // sample the time mask every 'tInterval' seconds
    // to get the Good Time Intervals (hereafter: GTIs)
    double tInterval = 0.0001 ;
    
    // break out of search loop automatically after 'maxSeconds' seconds
    double maxSeconds = 3000.0 ;
    
    // vectors for storing the times of the
    // beginning and end of each GTI
    vector<double> gti_Beginning ;
    vector<double> gti_Ending ;
    
    // find GTI edges
    calcGTIsFromTimeMask( vtm, tInterval, maxSeconds, startMJD, stopMJD, gti_Beginning, gti_Ending ) ;
    
    // if no gti's were found in VTimeMask, add the whole run manually
    if( gti_Beginning.size() == 0 && gti_Ending.size() == 0 )
    {
        gti_Beginning.push_back( startMJD ) ;
        gti_Ending.push_back( stopMJD ) ;
    }
    
    printf( "    Loaded gti list from MJD %.4f to %.4f\n", gti_Beginning[0], gti_Ending[gti_Ending.size() - 1] )  ;
    //cout.precision( 4 ) ;
    //cout << "    startMJD:" << startMJD << endl;
    //cout << "    stopMJD :" << stopMJD  << endl;
    if( gti_Beginning.size() > 0 )
    {
        //printf( "    Edges:  Begin    - Ending\n" );
        for( unsigned int i = 0 ; i < gti_Beginning.size() ; i++ )
        {
            printf( "    gti[%2d]: %6.1f - %6.1f\n", i + 1,
                    ( gti_Beginning[i] - met_zero ) * 86400.0 - tstart,
                    ( gti_Ending[i]   - met_zero ) * 86400.0 - tstart )  ;
            //printf( "           %9.4f - %9.4f \n", ( gti_Beginning[i] - startMJD ) * 24 * 60 * 60, ( gti_Ending[i] - startMJD ) * 24 * 60 * 60 ) ;
        }
    }
    
    // break each gti up into chunks
    // we must be careful to avoid creating empty chunks outside the gti
    vector <Ntimespan> chunklist ;
    //cout << endl;
    int nchunks ;
    for( unsigned int i_gti = 0 ; i_gti < gti_Beginning.size() ; i_gti++ )
    {
        // begining and ending of this GTI
        double gti_beg = ( gti_Beginning[i_gti] - met_zero ) * 86400.0 ; // MET
        double gti_end = ( gti_Ending[   i_gti] - met_zero ) * 86400.0 ; // MET
        
        if( chunkLength > 0 )
        {
            // loop over each chunk
            nchunks = ceil( ( gti_end - gti_beg ) / chunkLength ) ;
            //printf( "gti:   %6.1f - %6.1f ( %d chunks )\n", gti_beg-tstart, gti_end-tstart, nchunks ) ;
            double ch_beg, ch_end ;
            for( int i_chunk = 0 ; i_chunk < nchunks ; i_chunk++ )
            {
                // figure out chunk boundaries
                ch_beg = gti_beg + ( i_chunk    * chunkLength ) ;
                ch_end = ch_beg  + chunkLength ;
                
                // don't go past the end of the gti with this chunk
                if( ch_end > gti_end )
                {
                    ch_end = gti_end ;
                }
                
                // dump the tiny runs we get from some
                // poorly chosen chunk lengths and rounding stuff
                if( ch_end - ch_beg < 0.001 )
                {
                    continue ;
                }
                
                // create timespan and add to chunklist
                Ntimespan tmpspan ;
                tmpspan.beg = ch_beg ;
                tmpspan.end = ch_end ;
                chunklist.push_back( tmpspan ) ;
                //printf( " ch%2d: %6.1f - %6.1f (%fs long)\n", (int)chunklist.size(), tmpspan.beg-tstart, tmpspan.end-tstart, tmpspan.end-tmpspan.beg ) ;
            }
        }
        else
        {
            // create timespan and add to chunklist
            Ntimespan tmpspan ;
            tmpspan.beg = gti_beg ;
            tmpspan.end = gti_end ;
            chunklist.push_back( tmpspan ) ;
            //printf( " g%2d: %6.1f - %6.1f (%fs long)\n", (int)chunklist.size(), tmpspan.beg-tstart, tmpspan.end-tstart, tmpspan.end-tmpspan.beg ) ;
            
        }
        
    }
    
    // get $VERITAS_DATA_DIR environment variable
    char veritas_data_dir[200] = "" ; //"/lustre/fs5/group/cta/VERITAS" ;
    sprintf( veritas_data_dir, "%s", get_env_var( "VERITAS_DATA_DIR" ).c_str() ) ;
    
    // extra settings
    double specind = 1.5 ;
    
    // other vars needed later
    double offset = 0.0 ;
    
    // Load list of usable azimuths, zeniths, and offsets
    // within effective area files
    VMagnetZenith*   magZe     = new VMagnetZenith( effFile ) ;
    VMagnetOffset*   magOffset = new VMagnetOffset( effFile )  ;
    VMagnetAzimuth* magAz     = new VMagnetAzimuth( effFile ) ;
    //magOffset->printAll();
    cout << endl;
    cout << "Loading Run-wise Magnets:" << endl;
    cout << "  offsets  : " << magOffset->print() << endl;
    cout << "  azimuths : " << magAz->print() << endl;
    cout << "  zeniths  : " << magZe->print() << endl;
    
    // root file for debug histograms
    char fnameRoot[200] ;
    TFile* debugRootFile = 0 ;
    FILE* debugAsciiFile = 0 ;
    if( writedebugroothists )
    {
        sprintf( fnameRoot, "%s/VR%d.debughists.root", outputDirectory.c_str(), runid ) ;
        debugRootFile = new TFile( fnameRoot, "RECREATE" ) ;
        
        sprintf( fnameRoot, "%s/VR%d.debughists.dat", outputDirectory.c_str(), runid ) ;
        debugAsciiFile = fopen( fnameRoot, "w" ) ;
    }
    
    // list of chunk files produced, saved to a file
    char chunkListName[1000] = "" ;
    sprintf( chunkListName, "%s/VR%d.chunklist", outputDirectory.c_str(), runid ) ;
    FILE* chunkList = 0 ;
    chunkList = fopen( chunkListName, "w" ) ;
    
    // print some info about the chunks
    //cout << endl;
    //printChunkList( tstart, tstop, chunklist ) ;
    
    /////////////////////////////////////////////////////////////////////////
    // LOOP OVER EACH CHUNK
    //printf("\nrun:%5d tstart:%.1f tstop:%.1f dur:%.1f\n", runid, tstart, tstop, tstop-tstart);
    //printf("chunklength:%.1f nchunks:%d\n", chunkLength, (int)chunklist.size() ) ;
    double ch_beg, ch_end, ch_dur, ch_MJDbeg, ch_MJDend ;
    for( unsigned int i_chunk = 0 ; i_chunk < chunklist.size() ; i_chunk++ )
    {
        // beginning and ending of this chunk
        ch_beg    = chunklist[i_chunk].beg ;
        ch_end    = chunklist[i_chunk].end ;
        ch_dur    = ch_end - ch_beg ;
        ch_MJDbeg = ( ch_beg / 86400.0 ) + met_zero ;
        ch_MJDend = ( ch_end / 86400.0 ) + met_zero ;
        telapse   = ch_dur ;
        
        // calculate the average az/alt/ra/dec for this chunk
        ra_pnt  = 0.0 ;
        dec_pnt = 0.0 ;
        alt_pnt = 0.0 ;
        az_pnt  = 0.0 ;
        pointinginput  pntinp ;
        pointingoutput pntout ;
        pntinp.ch       = chainPointData ;
        pntinp.beg      = ch_MJDbeg ;
        pntinp.end      = ch_MJDend ;
        pntout.ra  = -999.0 ;
        pntout.dec = -999.0 ;
        pntout.alt = -999.0 ;
        pntout.az  = -999.0 ;
        calcAvgPointing( pntinp, pntout ) ;
        ra_pnt  = pntout.ra  ;
        dec_pnt = pntout.dec ;
        alt_pnt = pntout.alt ;
        az_pnt  = pntout.az  ;
        //printf("CH%2d Time(%.1f,%.1f)  RADec(%f,%f)  AltAz(%f,%f)\n", i_chunk+1, ch_beg-tstart, ch_end-tstart, ra_pnt, dec_pnt, alt_pnt, az_pnt ) ;
        
        // figure out chunk's fits file name
        char fnameChar[200] ;
        sprintf( fnameChar, "VR%d.chunk%d.fits", runid, i_chunk + 1 ) ;
        if( chunkList )
        {
            fprintf( chunkList, "%s\n", fnameChar ) ;
        }
        char fnameFullChar[200] ;
        sprintf( fnameFullChar, "%s/%s", outputDirectory.c_str(), fnameChar ) ;
        string fitsFileName = fnameFullChar ;
        
        // check if output file exists
        if( file_exists( fitsFileName.c_str() ) )
        {
            // if clobber, then delete
            if( overwriteFlag )
            {
                if( remove( fitsFileName.c_str() ) != 0 )
                {
                    printf( "Error deleting file '%s', exiting...\n", fitsFileName.c_str() ) ;
                    return 1 ;
                }
            }
            else
            {
                cout << "Error, file '" << fitsFileName << "' already exists.  Use the '-f' option to overwrite this existing file." << endl;
                return 1 ;
            }
        }
        
        
        // get the deadtime for this chunk
        vector<bool> deadTimeMask ;
        for( double isec = 0.0 ; isec < tstop - tstart ; isec++ )
        {
            if( isec > ch_beg - tstart && isec < ch_end - tstart )
            {
                deadTimeMask.push_back( true ) ;
            }
            else
            {
                deadTimeMask.push_back( false ) ;
            }
        }
        //printDeadTimeMask( i_chunk+1, 0.0, tstop-tstart, deadTimeMask ) ;
        
        // how much of the run's time was spent reading out the electronics
        bool useTimeIndependentDeadTime = false ;
        double deadtimefrac = deadTime->getDeadTimeFraction( deadTimeMask, useTimeIndependentDeadTime ) ;
        
        // how much time was actually spent observing for gamma rays
        double livetime = ch_dur * ( 1.0 - deadtimefrac ) ;
        
        // figure out the needed effective area files
        printf( "\n\n#############\nChunk %d #####\n#############\n", i_chunk + 1 ) ;
        //printf("  Figure out needed effective area files...(alt_pnt=%.3f)\n", alt_pnt );
        double avgzen = 90.0 - alt_pnt ; // ze:0-90, alt:90-0
        int zenith = magZe->magnetize( avgzen ) ;
        int azbin  = magAz->magnetize( az_pnt ) ;
        //printf("    alt_pnt(%.1f) -> avgzen(%.1f) -> zenith(%2d)\n", alt_pnt, avgzen, zenith ) ;
        
        // figure out noises
        
        VMagnetNoise* magNo = new VMagnetNoise( effFile, azbin, zenith ) ;
        int noise = magNo->magnetize( pedvar ) ;
        //cout << "    noise list:" << magNo->printNoises() << endl;
        printf( "  selected azbin(%2d) zenith(%2d) pedvar(%.1f) -> noise(%d)\n", azbin, zenith, pedvar, noise ) ;
        
        // place to keep track of this chunk's effArea files
        vector<archiveStruct> effAreaChunk ;
        
        // loop over offsets
        //printf("\n    Looping over offsets...\n");
        printf( "assembling effarea archives for this chunk...\n" ) ;
        for( unsigned int i_offset = 0 ; i_offset < magOffset->offsets.size() ; i_offset++ )
        {
        
            // get this loop's offset
            offset = magOffset->offsets[i_offset] ;
            
            // store this effArea set in the archive,
            archiveStruct archline ;
            archline.zenith  = zenith  ;
            archline.azbin   = azbin   ;
            archline.noise   = noise   ;
            archline.specind = specind ;
            archline.offset  = offset  ;
            archiveEffArea->addArchive( archline ) ;
            
            // save this effArea set for later use in this chunk's loop
            effAreaChunk.push_back( archline ) ;
            
        }
        
        // load them if we haven't already
        printf( "  Loading needed irfs...\n" );
        archiveEffArea->loadAllArchives() ;
        
        // setup FITSrecords for events, telescope, and gti tables
        // -evlist
        FITSRecord recEVENTS( fitsFileName, eventListTemplate, "EVENTS" );
        recEVENTS.setVerbose( 2 );   // Dont show us adding every single column mapping.
        
        // setup EVENTS table headers
        recEVENTS.writeHeader( "CONV_RA" , 9999.9 ) ;
        recEVENTS.writeHeader( "CONV_DEC", 9999.9 ) ;
        recEVENTS.writeHeader( "CREATOR" , string( progname ) ) ;
        recEVENTS.writeHeader( "DATE_END", string( dateend ) ) ;          // human readable run-end   date string
        recEVENTS.writeHeader( "DATE_OBS", string( dateobs ) ) ;          // human readable run-start date string
        recEVENTS.writeHeader( "EQUINOX" , 2000.0 ) ;                     // which JXXXX RA/Dec epoch to use for all RA/Dec in here
        recEVENTS.writeHeader( "EXTNAME" , string( "EVENTS" ) ) ;
        recEVENTS.writeHeader( "OBJECT"  , evndispRunPar->fTargetName ) ;
        recEVENTS.writeHeader( "OBS_ID"  , runid ) ;
        recEVENTS.writeHeader( "RA_OBJ"  , evndispRunPar->fTargetRA ) ;
        recEVENTS.writeHeader( "DEC_OBJ" , evndispRunPar->fTargetDec ) ;  // Dec of the target object in view
        // fk5 = something related to a precession-rotation standard
        recEVENTS.writeHeader( "RADECSYS", string( "fk5" ) ) ;
        recEVENTS.writeHeader( "TELESCOP", string( "VERITAS" ) ) ;
        recEVENTS.writeHeader( "TIME_END", string( timeend ) ) ;
        recEVENTS.writeHeader( "TIME_OBS", string( timeobs ) ) ;
        recEVENTS.writeHeader( "TIMEREF" , string( "local" ) ) ;
        recEVENTS.writeHeader( "TIMESYS" , string( "TT" ) ) ;
        recEVENTS.writeHeader( "TIMEUNIT", string( "s" ) ) ;
        recEVENTS.writeHeader( "MJDREFI" , mjdrefi ) ;
        recEVENTS.writeHeader( "MJDREFF" , mjdreff ) ;
        recEVENTS.writeHeader( "TSTART"  , tstart ) ;
        recEVENTS.writeHeader( "TSTOP"   , tstop ) ;
        recEVENTS.writeHeader( "TELAPSE" , telapse ) ;
        recEVENTS.writeHeader( "ONTIME"  , telapse ) ;
        recEVENTS.writeHeader( "RA_PNT"  , ra_pnt ) ;                     // average pointing RA
        recEVENTS.writeHeader( "DEC_PNT" , dec_pnt ) ;                    // averate pointing Dec
        recEVENTS.writeHeader( "ALT_PNT" , alt_pnt ) ;
        recEVENTS.writeHeader( "AZ_PNT"  , az_pnt ) ;
        recEVENTS.writeHeader( "DEADC"   , deadtimefrac ) ;
        recEVENTS.writeHeader( "LIVETIME", livetime ) ;
        
        // energy/skydir/time boundaries for this run
        char dsval2[100] = "" ;
        char dsval3[100] = "" ;
        sprintf( dsval2, "CIRCLE(%f,%f,%f)" , ra_pnt, dec_pnt, 2.25 ) ;
        sprintf( dsval3, "0.02:500.0" ) ; // 20GeV - 500TeV
        recEVENTS.writeHeader( "DSVAL2", string( dsval2 ) ) ;
        recEVENTS.writeHeader( "DSVAL3", string( dsval3 ) ) ;
        
        // VERITAS Event Parameters (start with caps?)
        int    RunNumber        =     0   ;
        int    EventNumber      =     0   ;
        int    NImages          =     0   ;
        int    MJD              =     0   ;
        double Time             =     0.0 ;
        double RA               =     0.0 ;
        double DEC              =     0.0 ;
        double XGroundCore      =     0.0 ;
        double YGroundCore      =     0.0 ;
        double Energy_ErecS     =     0.0 ;
        double Energy_ErecS_Err =     0.0 ;
        double Az               =     0.0 ;
        double El               =     0.0 ;
        double MSCW             = -9999.0 ;
        double MSCL             = -9999.0 ;
        double EmissionHeight   = -9999.0 ;
        double Xoff             =     0.0 ;
        double Yoff             =     0.0 ;
        //int    ImgSel           =     0   ;
        ULong64_t    ImgSel           =     0   ;
        
        // assign c++ variables to root branches
        chainEventList->SetBranchAddress( "runNumber"     , &RunNumber ) ;
        chainEventList->SetBranchAddress( "eventNumber"   , &EventNumber ) ;
        chainEventList->SetBranchAddress( "dayMJD"        , &MJD ) ;
        chainEventList->SetBranchAddress( "timeOfDay"     , &Time ) ;
        chainEventList->SetBranchAddress( "NImages"       , &NImages ) ;
        chainEventList->SetBranchAddress( "RA"            , &RA ) ;
        chainEventList->SetBranchAddress( "DEC"           , &DEC ) ;
        chainEventList->SetBranchAddress( "XGroundCore"   , &XGroundCore ) ;
        chainEventList->SetBranchAddress( "YGroundCore"   , &YGroundCore ) ;
        chainEventList->SetBranchAddress( "EnergyS"       , &Energy_ErecS ) ;
        chainEventList->SetBranchAddress( "EnergyS_Err"   , &Energy_ErecS_Err ) ;
        chainEventList->SetBranchAddress( "Az"            , &Az ) ;
        chainEventList->SetBranchAddress( "El"            , &El ) ;
        chainEventList->SetBranchAddress( "MSCW"          , &MSCW ) ;
        chainEventList->SetBranchAddress( "MSCL"          , &MSCL ) ;
        chainEventList->SetBranchAddress( "EmissionHeight", &EmissionHeight ) ;
        chainEventList->SetBranchAddress( "Xoff"          , &Xoff ) ;
        chainEventList->SetBranchAddress( "Yoff"          , &Yoff ) ;
        chainEventList->SetBranchAddress( "ImgSel"        , &ImgSel ) ;
        
        // ctools event list variables
        long  int obs_id      = 0   ;
        long  int event_id    = 0   ;
        double    event_time  = 0.0 ;
        short int multip      = 0   ;
        float     ra          = 0.0 ;
        float     dec         = 0.0 ;
        float     dir_err     = 0.0 ;
        float     corex       = 0.0 ;
        float     corey       = 0.0 ;
        //float     cor_err     = 0.0 ;
        float     energy      = 0.0 ;
        float     energyerr   = 0.0 ;
        float     az          = 0.0 ;
        float     alt         = 0.0 ;
        float     hil_msw     = 0.0 ;
        float     hil_msl     = 0.0 ;
        float     hil_msw_err = 0.0 ;
        float     hil_msl_err = 0.0 ;
        float     xmax        = 0.0 ;
        float     xmax_err    = 0.0 ;
        bool*     telmask     = new bool[NTEL];
        float     detx        = 0.0 ;
        float     dety        = 0.0 ;
        //float     shwidth     = 0.0 ;
        //float     shlength    = 0.0 ;
        
        // align c++ variables to event list variables
        recEVENTS.mapColumnToVar( "OBS_ID"     , obs_id ) ;
        recEVENTS.mapColumnToVar( "EVENT_ID"   , event_id ) ;
        recEVENTS.mapColumnToVar( "TIME"       , event_time ) ;
        recEVENTS.mapColumnToVar( "MULTIP"     , multip ) ;
        recEVENTS.mapColumnToVar( "RA"         , ra ) ;
        recEVENTS.mapColumnToVar( "DEC"        , dec ) ;
        recEVENTS.mapColumnToVar( "DIR_ERR"    , dir_err ) ;
        recEVENTS.mapColumnToVar( "COREX"      , corex ) ;
        recEVENTS.mapColumnToVar( "COREY"      , corey ) ;
        recEVENTS.mapColumnToVar( "ENERGY"     , energy ) ;
        recEVENTS.mapColumnToVar( "ENERGY_ERR" , energyerr ) ;
        recEVENTS.mapColumnToVar( "AZ"         , az ) ;
        recEVENTS.mapColumnToVar( "ALT"        , alt ) ;
        recEVENTS.mapColumnToVar( "HIL_MSW"    , hil_msw ) ;
        recEVENTS.mapColumnToVar( "HIL_MSW_ERR", hil_msw_err ) ;
        recEVENTS.mapColumnToVar( "HIL_MSL"    , hil_msl ) ;
        recEVENTS.mapColumnToVar( "HIL_MSL_ERR", hil_msl_err ) ;
        recEVENTS.mapColumnToVar( "XMAX"       , xmax ) ;
        recEVENTS.mapColumnToVar( "XMAX_ERR"   , xmax_err ) ;
        recEVENTS.mapColumnToVar( "TELMASK"    , telmask ) ;
        recEVENTS.mapColumnToVar( "DETX"       , detx ) ;
        recEVENTS.mapColumnToVar( "DETY"       , dety ) ;
        
        // variables for figuring out alt/az bounds
        bool altAzEmpty = true ;
        double az_min  = 0.0 ;
        double az_max  = 0.0 ;
        double alt_min = 0.0 ;
        double alt_max = 0.0 ;
        
        // fill events to FITSRecord EVENTS Table
        unsigned int totalEvents = chainEventList->GetEntries() ;
        printf( "looping over %d photons...\n", totalEvents ) ;
        
        int eventcounter = 0 ;
        if( totalEvents > 0 )
        {
            for( unsigned int i = 0; i <= totalEvents - 1 ; i++ )
            {
                chainEventList->GetEntry( i );
                
                // convert MJD and Time to UTC
                event_time  = ( ( MJD - met_zero ) * 86400.0 ) + Time ; // MET
                
                // check that the event is in our chunk's time window
                if( event_time < ch_beg || event_time > ch_end )
                {
                    continue ;
                }
                
                // assign other event parameters
                // variable names in all lowercase are CTOOLs
                // variable names with any capital letters are VERITAS variables
                obs_id      = RunNumber      ;
                event_id    = EventNumber    ;
                multip      = NImages        ;
                ra          = RA             ;
                dec         = DEC            ;
                dir_err     = -9999.9        ;
                corex       = XGroundCore    ;
                corey       = YGroundCore    ;
                //cor_err     = -9999.9        ;
                energy      = Energy_ErecS   ;
                energyerr   = -9999.9        ; // Energy_ErecS_Err ?
                az          = Az             ;
                alt         = El             ;
                hil_msw     = MSCW           ;
                hil_msw_err = -9999.9        ;
                hil_msl     = MSCL           ;
                hil_msl_err = -9999.9        ;
                xmax        = -9999.9        ;
                xmax_err    = -9999.9        ;
                detx        = Xoff           ;
                dety        = Yoff           ;
                //shwidth     = -9999.9        ;
                //shlength    = -9999.9        ;
                //mjd         = ( double )MJD + ( Time / 86400.0 ) ; // decimal MJD
                
                // figure out alt/az event range
                if( altAzEmpty )
                {
                    // set initial values, if not already set
                    az_min  = az  ;
                    az_max  = az  ;
                    alt_min = alt ;
                    alt_max = alt ;
                    altAzEmpty = false ;
                }
                else
                {
                    // otherwise, check to see if we've broken any records so far
                    if( az  < az_min )
                    {
                        az_min  = az  ;
                    }
                    if( az  > az_max )
                    {
                        az_max  = az  ;
                    }
                    if( alt < alt_min )
                    {
                        alt_min = alt ;
                    }
                    if( alt > alt_max )
                    {
                        alt_max = alt ;
                    }
                }
                
                // print raw event info, if asked
                if( printEventsFlag )
                {
                    printf( "EVENT %lu %lu %f %f %f\n", obs_id, event_id, ra, dec, event_time ) ;
                }
                
                // count the number of events in this chunk
                eventcounter++ ;
                //if( eventcounter < 10 ) printf( "    id:%6lu beg:%6.1f ev:%7.2f end:%6.1f\n", event_id, ch_beg-tstart, event_time-tstart, ch_end-tstart ) ;
                
                // set the telmask
                telmask_clear( telmask ) ;
                for( unsigned int j_tel = 1 ; j_tel <= NTEL ; j_tel++ )
                {
                    telmask_set( telmask, j_tel, readImgSel( ImgSel, j_tel ) ) ;
                }
                //printf("%5d T%d%d%d%d\n", event_id, telmask[0], telmask[1], telmask[2], telmask[3] ) ;
                
                
                // write the event
                recEVENTS.write() ;
                
            } // endfor: loop over events EVENTS table
            
        } // endif: if totalEvents > 0
        else
        {
            printf( "Warning, no events available to be written to the EVENTS table for chunk %d...\n", i_chunk + 1 ) ;
        }
        
        recEVENTS.finishWriting();
        printf( "recEVENTS written...\n" );
        
        ///////////////////////////////////////////
        // add in extra keywords to EVENTS header
        ///////////////////////////////////////////
        int status ;
        fitsfile* newfptr ;
        char ds_pos[100] = "" ;
        sprintf( ds_pos, "CIRCLE(%f,%f,2.0)", ra_pnt, dec_pnt ) ; // for DSVAL2, 2.0 deg radius ~ 3.5 deg fov
        
        // open our fits file
        fits_open_file( &newfptr, fitsFileName.c_str(), READWRITE, &status ) ;
        printf( " fits_open_file() status=%d\n", status ) ;
        
        // move to the EVENTS table
        fits_movnam_hdu( newfptr, BINARY_TBL, "EVENTS", 0, &status ) ;
        printf( " fits_movnam_hdu() status=%d\n", status ) ;
        
        // create keys
        fits_write_key( newfptr, TSTRING, "DSTYP1", ( void* )"TIME"       , "data time selection"          , &status ) ; //printf(" fits_write_key(DSTYP1) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSUNI1", ( void* )"s"          , "data time selection unit"     , &status ) ; //printf(" fits_write_key(DSUNI1) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSVAL1", ( void* )"TABLE"      , "data time selection value"    , &status ) ; //printf(" fits_write_key(DSVAL1) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSREF1", ( void* )":GTI"       , "data time selection reference", &status ) ; //printf(" fits_write_key(DSREF1) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSTYP2", ( void* )"POS(RA,DEC)", "data skydir selection"        , &status ) ; //printf(" fits_write_key(DSTYP2) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSUNI2", ( void* )"deg"        , "data skydir selection unit"   , &status ) ; //printf(" fits_write_key(DSUNI2) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSVAL2", ( void* )ds_pos       , "data skydir selection value"  , &status ) ; //printf(" fits_write_key(DSVAL2) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSTYP3", ( void* )"ENERGY"     , "data energy selection"        , &status ) ; //printf(" fits_write_key(DSTYP3) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSUNI3", ( void* )"TeV"        , "data energy selection unit"   , &status ) ; //printf(" fits_write_key(DSUNI3) status=%d\n", status ) ;
        fits_write_key( newfptr, TSTRING, "DSVAL3", ( void* )"0.03:300.0" , "data energy selection value"  , &status ) ; //printf(" fits_write_key(DSVAL3) status=%d\n", status ) ;
        
        // close file
        //fits_close_file( newfptr, &status ) ; printf(" fits_close_file() status=%d\n", status ) ;
        
        printf( "finished writing chunk %d , wrote %d events , spanning %.1f seconds\n", i_chunk + 1, eventcounter, ch_dur ) ;
        
        ///////////////////////////////////
        // setup TELARRAY Table
        // write header keywords
        FITSRecord recTELARRAY( fitsFileName, eventListTemplate, "TELARRAY" );
        recTELARRAY.setVerbose( 0 );
        recTELARRAY.writeHeader( "TELESCOP", string( "VERITAS" ) ) ;
        recTELARRAY.writeHeader( "ARRAY"   , string( "null" ) ) ;
        
        // setup FITSrecords columns
        short int telid    = 0   ;
        long long int subclass = 0   ;
        double posx  = 0.0 ;
        double posy  = 0.0 ;
        double posz  = 0.0 ;
        double mirror_area = 0.0 ;
        double cam_area    = 0.0 ;
        double focal_length = 0.0 ;
        double fieldofview  = 0.0 ;
        recTELARRAY.mapColumnToVar( "TELID"    , telid );
        recTELARRAY.mapColumnToVar( "SUBCLASS" , subclass );
        recTELARRAY.mapColumnToVar( "POSX"     , posx );
        recTELARRAY.mapColumnToVar( "POSY"     , posy );
        recTELARRAY.mapColumnToVar( "POSZ"     , posz );
        recTELARRAY.mapColumnToVar( "MIRAREA"  , mirror_area );
        recTELARRAY.mapColumnToVar( "CAMAREA"  , cam_area );
        recTELARRAY.mapColumnToVar( "FOCLEN"   , focal_length );
        recTELARRAY.mapColumnToVar( "FOV"      , fieldofview );
        
        // loop over telescopes
        totalEvents = chainTelConfig->GetEntries() ;
        for( unsigned int i = 0 ; i <= totalEvents - 1 ; i++ )
        {
            if( i > 10 )
            {
                break ;
            }
            chainTelConfig->GetEntry( i );
            telid        = 0 ;
            subclass     = 0 ;
            posx         = 0.0 ;
            posy         = 0.0 ;
            posz         = 0.0 ;
            mirror_area  = 0.0 ;
            cam_area     = 0.0 ;
            focal_length = 0.0 ;
            fieldofview  = 0.0 ;
            recTELARRAY.write() ;
        }
        
        recTELARRAY.finishWriting();
        
        //printf( "  chunk %2d - beg:%6.1f end:%6.1f nevents:%5d fits:%s - min/max event: { alt(%.1f-%.1f) az(%.1f-%.1f) }\n", i_chunk+1, ch_beg-tstart, ch_end-tstart, eventcounter, fitsFileName.c_str(), alt_min, alt_max, az_min, az_max ) ;
        //printf("ch%2d dead time fraction %f\n", i_chunk+1, deadtimefrac ) ;
        printf( "  chunk %2d - beg:%6.1f end:%6.1f nevents:%5d fits:%-38s - RADec(%.5f,%.5f) AltAz(%.1f,%.1f) - dead time frac:%6.4f\n", i_chunk + 1, ch_beg - tstart, ch_end - tstart, eventcounter, fitsFileName.c_str(), ra_pnt, dec_pnt, alt_pnt, az_pnt, deadtimefrac ) ;
        
        ///////////////////////////////////
        // setup GTI Table
        // Good Time Interval
        // in VERITAS, this is the INVERSE table of Time Cuts
        //   VERITAS Time Cuts: list of times that are *UNUSABLE*
        //   GTI: list of times that are *USABLE*
        //   start = when to start the cut time period, from db
        //   stop  = when to stop the cut time period, from db
        //recGTI.writeHeader( "" , ) ;
        FITSRecord recGTI( fitsFileName, eventListTemplate, "GTI" );
        recGTI.setVerbose( 0 );
        recGTI.writeHeader( "MJDREFI" , mjdrefi ) ;
        recGTI.writeHeader( "MJDREFF" , mjdreff ) ;
        recGTI.writeHeader( "TSTART"  , tstart ) ;
        recGTI.writeHeader( "TSTOP"   , tstop ) ;
        recGTI.writeHeader( "TIMESYS" , string( "TT" ) ) ;
        recGTI.writeHeader( "TIMEREF" , string( "local" ) ) ;
        recGTI.writeHeader( "TIMEUNIT", string( "s" ) ) ;
        double gti_start  = 0.0 ;
        double gti_stop   = 0.0 ;
        recGTI.mapColumnToVar( "START" , gti_start );
        recGTI.mapColumnToVar( "STOP"  , gti_stop );
        
        // add our chunk's time to the gti
        // since we first looped over each Good Time Interval,
        // we know that this entire chunk is valid.
        gti_start = ch_beg ;
        gti_stop  = ch_end ;
        recGTI.write() ;
        
        // cleanup FITSrecords
        recGTI.finishWriting();
        
        // declare some variables for adding irf tables
        fitsfile* chunkptr = 0 ;
        status = 0    ; // 0=function ok, >0 is an error
        char tblname[100] ;
        int  nrows        ;
        int  tfields = 0;
        
        ///////////////////////////////////////////////////
        // POINT SPREAD FUNCTION table
        // -psf table
        // move to a single while-loop, so variables can be redeclared when
        // adding each IRF, and code can just be copied between IRF sections
        bool writePOINTSPREADFUNCTIONtable = true;
        if( skippsf )
        {
            writePOINTSPREADFUNCTIONtable = false ;
        }
        int psfmode = PSFTABLE_BASIC ;
        int psffunc = PSF_KING ; //PSF_GAUSS ;
        int kingmode = KING_NATIVE ; // KING_GUESS
        if( psffunc != PSF_GAUSS && psffunc != PSF_KING )
        {
            printf( "Error, unrecognized psffunc(%d) value, exiting...\n", psffunc ) ;
            exit( 1 );
        }
        while( writePOINTSPREADFUNCTIONtable )
        {
            if( psfmode == PSFTABLE_BASIC )
            {
                // open fits file for editing
                if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup table columns
                sprintf( tblname, "%s", "POINT SPREAD FUNCTION" ) ;
                nrows     = 1  ; // psf only has one row
                
                // figure out how many columns to work with
                if( psffunc == PSF_GAUSS )
                {
                    tfields = 10 ;
                }
                else if( psffunc == PSF_KING )
                {
                    tfields =  6 ;
                }
                
                // how many elements should our PSF table have?
                int archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[0] ) ; // look at first offset's psf vector for the energy bins
                int nEnergyRows  = archiveEffArea->fPsf[archiveIndex].energy.size() ;
                int nOffsetRows  = magOffset->offsets.size() ;
                int nelem_xy     = nEnergyRows * nOffsetRows ; // number of elements in our 2D arrays
                
                // min and max energy for the histograms
                // will be broken up into <nEnergyRows> bins of equal log10-width
                double lowenergybin     = archiveEffArea->fPsf[archiveIndex].energy[ 0]    ; // center of first  energy bin (log10 TeV)
                double nextenergybin    = archiveEffArea->fPsf[archiveIndex].energy[ 1]    ; // center of second energy bin (log10 TeV)
                double highenergybin    = archiveEffArea->fPsf[archiveIndex].energy.back() ; // center of last   energy bin (log10 TeV)
                double logenergyhalfwid = ( nextenergybin - lowenergybin ) / 2.0           ; // half the width of one energy bin (log10 TeV)
                float  energyMin        = pow( 10.0, lowenergybin  - logenergyhalfwid )    ; // lower edge of the lowest  energy bin (TeV)
                float  energyMax        = pow( 10.0, highenergybin + logenergyhalfwid )    ; // upper edge of the highest energy bin (TeV)
                printf( "\npsf table will have %d bins from %.3f to %.3f TeV...\n", nEnergyRows, energyMin, energyMax ) ;
                
                // figure out column names, formats, and units
                char form_x [11] ;
                char form_y [11] ;
                char form_xy[11] ;
                sprintf( form_x,  "%dE", nEnergyRows ) ;
                sprintf( form_y,  "%dE", nOffsetRows ) ;
                sprintf( form_xy, "%dE", nelem_xy ) ;
                // 21E = each row element in ENERG_LO and ENERG_HI will be a list of 21 floats('E')
                // 9E  = each row element in THETA_LO and THETA_HI will be a list if  9 floats('E')
                
                // setup variables for column numbers
                int C_ENERG_LO =  1 ;
                int C_ENERG_HI =  2 ;
                int C_THETA_LO =  3 ;
                int C_THETA_HI =  4 ;
                
                // psf_gauss columns
                int C_SCALE    =  5 ;
                int C_SIGMA_1  =  6 ;
                int C_AMPL_2   =  7 ;
                int C_SIGMA_2  =  8 ;
                int C_AMPL_3   =  9 ;
                int C_SIGMA_3  = 10 ;
                
                // psf king
                int C_GAMMA    =  5 ;
                int C_SIGMA    =  6 ;
                
                // define our table column names, types, and units
                char* ttype[] = { "ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "SCALE", "SIGMA_1", "AMPL_2", "SIGMA_2", "AMPL_3", "SIGMA_3" } ;
                char* tform[] = { form_x    , form_x    , form_y    , form_y    , form_xy, form_xy  , form_xy , form_xy  , form_xy , form_xy   } ;
                char* tunit[] = { "TeV"     , "TeV"     , "deg"     , "deg"     , "\0"   , "deg"    , "\0"   , "deg"     , "\0"    , "deg"     } ;
                // use "\0" for no unit
                
                // change some of the columns if we're using a king function
                if( psffunc == PSF_KING )
                {
                    ttype[C_GAMMA - 1] = "GAMMA" ;
                    ttype[C_SIGMA - 1] = "SIGMA" ;
                }
                
                // create table
                if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // set some columns to have a dimensionality <nEnergyRows>x<nOffsetRows>
                // what this means: each row element in the GAMMA, SCALE, SIGMA* and AMPL* columns
                // was set to be a list of 189 floats (189E) (in the ttform variable above), but below we
                // reorganize these 189 floats into a 2D array of 21x9 floats (21*9=189)
                int  naxis  = 2   ; // (XX,XX), I think naxis=3 would be (XX,XX,XX), etc.
                long naxes[naxis] ;
                //naxes[0] = nOffsetRows ;
                //naxes[1] = nEnergyRows ;
                naxes[0] = nEnergyRows ;
                naxes[1] = nOffsetRows ;
                if( psffunc == PSF_GAUSS )
                {
                    if( fits_write_tdim( chunkptr, C_SCALE  , naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_tdim( chunkptr, C_SIGMA_1, naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_tdim( chunkptr, C_AMPL_2 , naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_tdim( chunkptr, C_SIGMA_2, naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_tdim( chunkptr, C_AMPL_3 , naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_tdim( chunkptr, C_SIGMA_3, naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                }
                else
                {
                    if( fits_write_tdim( chunkptr, C_GAMMA  , naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_tdim( chunkptr, C_SIGMA  , naxis, naxes, &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                }
                
                // create data maps for each column
                float data_ENERG_LO[nEnergyRows] ;
                float data_ENERG_HI[nEnergyRows] ;
                float data_THETA_LO[nOffsetRows] ;
                float data_THETA_HI[nOffsetRows] ;
                float data_SCALE   [nOffsetRows][nEnergyRows] ; // gauss
                float data_SIGMA_1 [nOffsetRows][nEnergyRows] ; // gauss
                float data_AMPL_2  [nOffsetRows][nEnergyRows] ; // gauss
                float data_SIGMA_2 [nOffsetRows][nEnergyRows] ; // gauss
                float data_AMPL_3  [nOffsetRows][nEnergyRows] ; // gauss
                float data_SIGMA_3 [nOffsetRows][nEnergyRows] ; // gauss
                float data_GAMMA   [nOffsetRows][nEnergyRows] ; // king
                float data_SIGMA   [nOffsetRows][nEnergyRows] ; // king
                
                // loop variables
                double rad68  = 0.0 ;
                double rad80  = 0.0 ;
                double sigma1 = 0.0 ; // gauss
                //double scale  = 0.0 ; // gauss
                double sigma  = 0.0 ; // king
                double gamma  = 0.0 ; // king
                double contain68check = 0.0 ;
                double contain80check = 0.0 ;
                Double_t pars[2] ;
                
                // double arrays for storing and checking the psf parameters
                // before they are loaded into the fits table
                VPsfBlock* psfBlock = new VPsfBlock( nOffsetRows, nEnergyRows ) ;
                int goodpsffitpoints  = 0 ;
                int totalpsffitpoints = 0 ;
                
                // loop over psf parameter space, fitting for a king function
                // first do fit of king function, apply interpolated values later
                for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
                {
                
                    // get this offset's archive
                    archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                    
                    // loop over energies
                    for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                    {
                        // set default block values
                        psfBlock->fPsfBlock[i_offset][i_en].gamma = 0.0 ;
                        psfBlock->fPsfBlock[i_offset][i_en].sigma = 0.0 ;
                        psfBlock->fPsfBlock[i_offset][i_en].rad68 = 0.0 ;
                        psfBlock->fPsfBlock[i_offset][i_en].rad80 = 0.0 ;
                        psfBlock->fPsfBlock[i_offset][i_en].i_offset  = i_offset ;
                        psfBlock->fPsfBlock[i_offset][i_en].i_en      = i_en     ;
                        psfBlock->fPsfBlock[i_offset][i_en].archindex = archiveIndex ;
                        psfBlock->fPsfBlock[i_offset][i_en].goodfit         = false ;
                        psfBlock->fPsfBlock[i_offset][i_en].wasinterpolated = false ;
                        psfBlock->fPsfBlock[i_offset][i_en].nativeking_sigma = 0.0 ;
                        psfBlock->fPsfBlock[i_offset][i_en].nativeking_gamma = 0.0 ;
                        
                        // get native fit king sigma/gamma parameters
                        psfBlock->fPsfBlock[i_offset][i_en].nativeking_sigma = archiveEffArea->fPsf[archiveIndex].nativeking_sigma[i_en] ;
                        psfBlock->fPsfBlock[i_offset][i_en].nativeking_gamma = archiveEffArea->fPsf[archiveIndex].nativeking_gamma[i_en] ;
                        
                        // get containment radii and do TMinuit fit to get king parameters
                        psfBlock->fPsfBlock[i_offset][i_en].rad68 = archiveEffArea->fPsf[archiveIndex].rad68[i_en] ;
                        psfBlock->fPsfBlock[i_offset][i_en].rad80 = archiveEffArea->fPsf[archiveIndex].rad80[i_en] ;
                        sigma = 0.0 ;
                        gamma = 0.0 ;
                        calcKingParametersRootFit( psfBlock->fPsfBlock[i_offset][i_en].rad68,
                                                   psfBlock->fPsfBlock[i_offset][i_en].rad80,
                                                   sigma, gamma ) ;
                        psfBlock->fPsfBlock[i_offset][i_en].sigma = sigma ;
                        psfBlock->fPsfBlock[i_offset][i_en].gamma = gamma ;
                        
                        // plug the fitted sigma/gamma back in to verify the fit converged
                        contain68check = 0.0 ;
                        contain80check = 0.0 ;
                        pars[0] = psfBlock->fPsfBlock[i_offset][i_en].sigma ;
                        pars[1] = psfBlock->fPsfBlock[i_offset][i_en].gamma ;
                        contain68check = kingfunc( psfBlock->fPsfBlock[i_offset][i_en].rad68, pars ) ;
                        contain80check = kingfunc( psfBlock->fPsfBlock[i_offset][i_en].rad80, pars ) ;
                        totalpsffitpoints += 1 ;
                        
                        // check if we got a good fit
                        if( fabs( contain68check - 0.68 ) < 0.015 &&
                                fabs( contain80check - 0.80 ) < 0.015 &&
                                psfBlock->fPsfBlock[i_offset][i_en].rad68 > 0.0 &&
                                psfBlock->fPsfBlock[i_offset][i_en].rad80 > 0.0 )
                        {
                            psfBlock->fPsfBlock[i_offset][i_en].goodfit = true  ;
                            goodpsffitpoints += 1 ;
                        }
                        
                        // save this point in the psf space
                        string str = formatPsfBlockStructString( psfBlock->fPsfBlock[i_offset][i_en] ) ;
                        printf( "saving %s | check68=%6.4f check80=%6.4f\n", str.c_str(), contain68check, contain80check ) ;
                    }
                }
                
                //if ( goodpsffitpoints != totalpsffitpoints )
                //{
                //  printf("run %d chunk %d : warning, only %d of %d points in psf parameter space were fit, the rest were extrapolated...\n", runid, i_chunk+1,  goodpsffitpoints, totalpsffitpoints );
                //}
                //else if ( goodpsffitpoints == 0 )
                //{
                //  printf("run %d chunk %d : error, could not read any points in the psf parameter space from effective area file %s , exiting...\n", runid, i_chunk+1, effFile.c_str() );
                //  exit(1);
                //}
                
                // check each parameter space's fit,
                // if bad fit found, scan nearby parameter space and just
                // copy that value found
                for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
                {
                    // loop over energies
                    for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                    {
                        // if we didn't get a good psf fit at this point in the parameter space,
                        // scan neighboring points, and just use the first one found
                        if( ! psfBlock->fPsfBlock[i_offset][i_en].goodfit )
                        {
                            psfBlock->scanBlockForAlternateFit( i_offset, i_en ) ;
                        }
                    }
                }
                
                //printf("\nfit fixing:\n");
                //psfBlock->formatBlock_goodfit() ;
                //psfBlock->formatBlock_wasinterpolated() ;
                //psfBlock->formatBlock_zerorad() ;
                //psfBlock->formatBlock_sanefit() ;
                
                for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
                {
                    // offset bin edges
                    data_THETA_LO[i_offset] = magOffset->offsetloweredges[i_offset] ;
                    data_THETA_HI[i_offset] = magOffset->offsetupperedges[i_offset] ;
                    
                    // create empty histograms
                    char hname1[200], hname2[200], hdesc1[200], hdesc2[200] ;
                    sprintf( hname1, "angres_chunk%d_offset%d_containment68", i_chunk, i_offset ) ;
                    sprintf( hname2, "angres_chunk%d_offset%d_containment80", i_chunk, i_offset ) ;
                    sprintf( hdesc1, "angular resolution at 68%% containment, for offset between %.2frad and %.2frad", data_THETA_LO[i_offset], data_THETA_HI[i_offset] );
                    sprintf( hdesc2, "angular resolution at 80%% containment, for offset between %.2frad and %.2frad", data_THETA_LO[i_offset], data_THETA_HI[i_offset] );
                    TH1F* angres68 = new TH1F( hname1, hdesc1, nEnergyRows, log10( energyMin ), log10( energyMax ) ) ;
                    TH1F* angres80 = 0 ;
                    if( psffunc == PSF_KING )
                    {
                        angres80 = new TH1F( hname2, hdesc2, nEnergyRows, log10( energyMin ), log10( energyMax ) ) ;
                        angres68->SetXTitle( "log_{10} (E/TeV)" );
                        angres80->SetXTitle( "log_{10} (E/TeV)" );
                        angres68->SetYTitle( "containment radius (68%) [deg]" );
                        angres80->SetYTitle( "containment radius (80%) [deg]" );
                    }
                    
                    // get this offset's archive
                    //if ( nOffsetRows != (int) effAreaChunk.size() ) { printf("Error, size mismatch between nOffsetRows and effAreaChunk.size(), exiting...\n"); exit(1) ; }
                    //archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                    
                    // fill the angular resolution histograms
                    for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                    {
                        angres68->SetBinContent( i_en, psfBlock->fPsfBlock[i_offset][i_en].rad68 ) ;
                        angres80->SetBinContent( i_en, psfBlock->fPsfBlock[i_offset][i_en].rad80 ) ;
                        //printf("filling angres68/80: %2d %.3f %.3f\n", i_en, archiveEffArea->fPsf[archiveIndex].rad68[i_en], archiveEffArea->fPsf[archiveIndex].rad80[i_en] );
                    }
                    //archiveEffArea->irf[archiveIndex]->fillResolutionHistogram( angres68, "68", "t_angular_resolution" ) ;
                    //if ( psffunc == PSF_KING ) archiveEffArea->irf[archiveIndex]->fillResolutionHistogram( angres80, "80", "t_angular_resolution" ) ;
                    //double minthresh68 = angres68->GetBinLowEdge( angres68->FindFirstBinAbove(0.0) ) ;
                    //double maxthresh68 = angres68->GetBinLowEdge( angres68->FindLastBinAbove( 0.0) ) + angres68->GetBinWidth( angres68->FindLastBinAbove(0.0) ) ;
                    //printf("    angres68 - min/max energies: %7.3f/%7.2f\n", pow(10,minthresh68), pow(10,maxthresh68) ) ;
                    
                    // write these filled histograms to a debug root file, if needed
                    if( debugRootFile )
                    {
                        debugRootFile->cd();
                        angres68->Write();
                        if( psffunc == PSF_KING )
                        {
                            angres80->Write();
                        }
                    }
                    
                    // loop over energy bins (log space)
                    printf( "\n" );
                    for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                    {
                    
                        // energy bin edges
                        data_ENERG_LO[i_en] = pow( 10, angres68->GetXaxis()->GetBinLowEdge( i_en + 1 ) ) ;
                        data_ENERG_HI[i_en] = pow( 10, angres68->GetXaxis()->GetBinUpEdge( i_en + 1 ) ) ;
                        
                        // 68% containment radius
                        rad68 = psfBlock->fPsfBlock[i_offset][i_en].rad68 ;
                        
                        // calculate gauss parameters
                        if( psffunc == PSF_GAUSS )
                        {
                            // calculate gauss scale and radius from the 68% containment radius
                            // taken from the $CTOOLS/scripts/cta_root2caldb.py: root2psf_gauss() function
                            // PDF(r) = norm * Exp( - r^2 / SIGMA_1^2 )
                            //   norm is automatically calculated by GCTAPsf2D
                            sigma1 = rad68 * R68_TO_1SIGMA ;
                            
                            // save parameters to our data arrays
                            data_SCALE  [i_en][i_offset] = 0.0 ; // no SCALE parameter needed, automtically handled by GCTAPsf2D Class as per email with Jurgen on 30/12/2014
                            data_SIGMA_1[i_offset][i_en] = sigma1 ;
                            data_AMPL_2 [i_offset][i_en] = 0.0 ;
                            data_SIGMA_2[i_offset][i_en] = 0.0 ;
                            data_AMPL_3 [i_offset][i_en] = 0.0 ;
                            data_SIGMA_3[i_offset][i_en] = 0.0 ;
                            data_SIGMA  [i_offset][i_en] = 0.0 ;
                            data_GAMMA  [i_offset][i_en] = 0.0 ;
                            printf( "  chunk(%2d)  offset(%4.2f:%4.2f)  energy(%5.1f:%5.1f) : psf gauss : scale=%.3f sigma1=%.3f\n", i_chunk + 1,
                                    data_THETA_LO      [i_offset],
                                    data_THETA_HI      [i_offset],
                                    data_ENERG_LO[i_en]          ,
                                    data_ENERG_HI[i_en]          ,
                                    data_SCALE   [i_en][i_offset],
                                    data_SIGMA_1 [i_en][i_offset] ) ;
                        }
                        else if( psffunc == PSF_KING )
                        {
                            // As per email from Jurgen on 30/12/2014
                            // PDF(r) = norm * ( ( 1.0 + 1/(2*GAMMA) ) * r^2 / SIGMA^2 )^-GAMMA
                            //   norm is calculated by Integrate[ 2*Pi * PDF(r) * Sin(r) , {r,0,inf} ]
                            // collect needed variables
                            rad80 = psfBlock->fPsfBlock[i_offset][i_en].rad80 ;
                            sigma = psfBlock->fPsfBlock[i_offset][i_en].sigma ;
                            gamma = psfBlock->fPsfBlock[i_offset][i_en].gamma ;
                            
                            if( kingmode == KING_NATIVE )
                            {
                                sigma = psfBlock->fPsfBlock[i_offset][i_en].nativeking_sigma ;
                                gamma = psfBlock->fPsfBlock[i_offset][i_en].nativeking_gamma ;
                                printf( "using native king values i_off:%2d  i_en:%2d  sigma:%5.3f  gamma:%7.2f\n", i_offset, i_en, sigma, gamma ) ;
                            }
                            
                            // add spike for debugging
                            //   ien
                            // i X000
                            // o 0X0X
                            // f 00X0
                            //char staticflag[100]="" ;
                            /*
                            printf("dingbling i_en=%d i_offset=%d!\n", i_en, i_offset );
                            if ( ( i_offset==i_en && i_en<3 ) || ( i_offset==2 && i_en==3 ) )
                            {
                              gamma = 155.55    ;
                              sigma =   0.09999 ;
                              printf("i_en=%d i_offset=%d forced: gamma=%f sigma=%f\n", i_en, i_offset, gamma, sigma ) ;
                              sprintf( staticflag, "%sstatic%s", KYEL, KNRM ) ;
                            }
                            printf("dingding!\n");
                            */
                            
                            // save our parameters to the data array
                            data_SCALE  [i_offset][i_en] = 0.0 ;
                            data_SIGMA_1[i_offset][i_en] = 0.0 ;
                            data_AMPL_2 [i_offset][i_en] = 0.0 ;
                            data_SIGMA_2[i_offset][i_en] = 0.0 ;
                            data_AMPL_3 [i_offset][i_en] = 0.0 ;
                            data_SIGMA_3[i_offset][i_en] = 0.0 ;
                            data_SIGMA  [i_offset][i_en] = sigma ;
                            data_GAMMA  [i_offset][i_en] = gamma ;
                            char kingwarn[100] = "" ;
                            if( sigma == 0 || gamma == 0 || rad68 == 0 || rad80 == 0 ) // || fabs(contain68check-0.68)>0.01 || fabs(contain80check-0.8)>0.01 )
                            {
                                sprintf( kingwarn, " %shey! listen!%s", KRED, KNRM ) ;
                            }
                            char rad68warn    [100] = "" ;
                            char rad80warn    [100] = "" ;
                            char sigmawarn    [100] = "" ;
                            char gammawarn    [100] = "" ;
                            if( rad68 == 0 )
                            {
                                sprintf( rad68warn, "%srad68=%.4f%s" , KRED, rad68, KNRM ) ;
                            }
                            else
                            {
                                sprintf( rad68warn,   "rad68=%.4f"   ,       rad68 ) ;
                            }
                            if( rad80 == 0 )
                            {
                                sprintf( rad80warn, "%srad80=%.4f%s" , KRED, rad80, KNRM ) ;
                            }
                            else
                            {
                                sprintf( rad80warn,   "rad80=%.4f"   ,       rad80 ) ;
                            }
                            if( sigma == 0 )
                            {
                                sprintf( sigmawarn, "%ssigma=%.4f%s" , KRED, sigma, KNRM ) ;
                            }
                            else
                            {
                                sprintf( sigmawarn,   "sigma=%.4f"   ,       sigma ) ;
                            }
                            if( gamma == 0 )
                            {
                                sprintf( gammawarn, "%sgamma=%8.4f%s", KRED, gamma, KNRM ) ;
                            }
                            else
                            {
                                sprintf( gammawarn,   "gamma=%8.4f"  ,       gamma ) ;
                            }
                            
                            // flag if it was fit properly
                            char fittag[100] = "" ;
                            if( psfBlock->fPsfBlock[i_offset][i_en].goodfit )
                            {
                                sprintf( fittag, "goodfit" ) ;
                            }
                            else if( psfBlock->fPsfBlock[i_offset][i_en].wasinterpolated )
                            {
                                sprintf( fittag, "stolen" ) ;
                            }
                            
                            //printf("  chunk(%2d)  offset(%d,%4.2f:%4.2f)  energy(%2d,%8.4f:%8.4f) : psf king tminuit : %s %s %s %s : fit=%-7s : %s : %s\n", i_chunk+1,
                            //  i_offset,
                            //  data_THETA_LO[i_offset],
                            //  data_THETA_HI[i_offset],
                            //  i_en,
                            //  data_ENERG_LO[i_en] ,
                            //  data_ENERG_HI[i_en] ,
                            //  rad68warn, rad80warn,
                            //  sigmawarn, gammawarn,
                            //  fittag,
                            //  kingwarn,
                            //  staticflag ) ;
                            
                        }
                    }
                }
                
                // write column arrays to their table columns
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_LO, 1, 1, nEnergyRows, &data_ENERG_LO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_HI, 1, 1, nEnergyRows, &data_ENERG_HI, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_THETA_LO, 1, 1, nOffsetRows, &data_THETA_LO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_THETA_HI, 1, 1, nOffsetRows, &data_THETA_HI, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( psffunc == PSF_GAUSS )
                {
                    if( fits_write_col( chunkptr, TFLOAT, C_SCALE   , 1, 1, nelem_xy   , &data_SCALE   , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_col( chunkptr, TFLOAT, C_SIGMA_1 , 1, 1, nelem_xy   , &data_SIGMA_1 , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_col( chunkptr, TFLOAT, C_AMPL_2  , 1, 1, nelem_xy   , &data_AMPL_2  , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_col( chunkptr, TFLOAT, C_SIGMA_2 , 1, 1, nelem_xy   , &data_SIGMA_2 , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_col( chunkptr, TFLOAT, C_AMPL_3  , 1, 1, nelem_xy   , &data_AMPL_3  , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_col( chunkptr, TFLOAT, C_SIGMA_3 , 1, 1, nelem_xy   , &data_SIGMA_3 , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                }
                else if( psffunc == PSF_KING )
                {
                    if( fits_write_col( chunkptr, TFLOAT, C_GAMMA   , 1, 1, nelem_xy   , &data_GAMMA   , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    if( fits_write_col( chunkptr, TFLOAT, C_SIGMA   , 1, 1, nelem_xy   , &data_SIGMA   , &status ) )
                    {
                        print_fits_error( status ) ;
                    }
                    
                }
                
                // close fits file
                if( fits_close_file( chunkptr, &status ) )
                {
                    print_fits_error( status ) ;
                }
            }
            else
            {
                printf( "Error, unrecognized pdfmode=%d, exiting...\n", psfmode ) ;
                exit( 1 );
            }
            
            break ;
        }
        
        
        
        printf( "\n\nstarting psf table for chunk %d...\n", i_chunk ) ;
        bool writePSFtable = false ;
        ///////////////////////////////////////////////////
        // POINT SPREAD FUNCTION table
        // -psf table
        // saves psf info to a fits table for gammalib.GCTAPsfKing
        // move to a single while-loop, so variables can be redeclared when
        // adding each IRF, and code can just be copied between IRF sections
        while( writePSFtable )
        {
        
            // open fits file for editing
            if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // setup table configuration
            sprintf( tblname, "%s", "POINT SPREAD FUNCTION" ) ;
            nrows   = 1 ;
            tfields = 6 ;
            
            int archiveIndex = 0 ;
            int nOffsetRows  = magOffset->offsets.size() ;
            
            // figure out energy binning (different offsets have different energy ranges)
            // since all offsets need to go into one table
            vector<double> e_bins ;
            for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
            {
                archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                
                // loop over each energy in this offset
                for( unsigned int i = 0 ; i < archiveEffArea->fPsf[archiveIndex].energy.size() ; i++ )
                {
                    e_bins.push_back( archiveEffArea->fPsf[archiveIndex].energy[ i ] ) ;
                }
            }
            
            // get sorted list of uniq energies
            sort( e_bins.begin(), e_bins.end() ) ;
            e_bins.erase( unique( e_bins.begin(), e_bins.end() ), e_bins.end() ) ;
            for( unsigned int i = 0 ; i < e_bins.size() ; i++ )
            {
                printf( "e_bins[%2d] = %.2f\n", i, e_bins[i] ) ;
            }
            
            // figure out smallest bin difference
            float en_logdelta = 1.0 ;
            float sld = 1.0 ;
            for( unsigned int i = 0 ; i < e_bins.size() - 1 ; i++ )
            {
                sld = e_bins[i + 1] - e_bins[i] ;
                if( sld < en_logdelta )
                {
                    en_logdelta = sld ;
                }
            }
            
            // figure out energy binning for the psf table
            double energyMin   = pow( 10, e_bins[0]               - ( en_logdelta / 2.0 ) ) ; // TeV
            double energyMax   = pow( 10, e_bins[e_bins.size() - 1] + ( en_logdelta / 2.0 ) ) ; // TeV
            double en_logmin   = log10( energyMin ) ;
            double en_logmax   = log10( energyMax ) ;
            int    nEnergyRows = e_bins.size() ;
            int    floored_nEn = floor( ( en_logmax - en_logmin ) / en_logdelta ) ;
            printf( "\n" );
            printf( "log energyMin : %.3f\n", en_logmin ) ;
            printf( "log energyMax : %.3f\n", en_logmax ) ;
            printf( "en_logdelta   : %.3f\n", en_logdelta ) ;
            printf( "nEnergyRows   : %d\n"  , nEnergyRows ) ;
            printf( "floored nEn   : %d\n"  , floored_nEn ) ;
            printf( "\n" );
            
            // check for weird error case
            if( nEnergyRows != floored_nEn )
            {
                printf( "Warning, psf energy bins don't quite match up, suspect a bin got skipped or a float error happened...\n" );
            }
            
            // print out our energy bins
            for( int i = 0 ; i < nEnergyRows ; i++ )
            {
                printf( "en bin %d : %.3f - %.3f (%7.3f - %7.3f TeV)\n", i, en_logmin + ( i * en_logdelta ), en_logmin + ( ( i + 1 )*en_logdelta ), pow( 10, en_logmin + ( i * en_logdelta ) ), pow( 10, en_logmin + ( ( i + 1 )*en_logdelta ) ) ) ;
            }
            printf( "\npsf table will have %d bins from %.3f to %.3f TeV...\n", nEnergyRows, energyMin, energyMax ) ;
            
            // how many elements should our PSF array columns have?
            int nelem_xy = nEnergyRows * nOffsetRows ; // number of elements in our 2D arrays
            
            // figure out column names, formats, and units
            char form_x [11] ;
            char form_y [11] ;
            char form_xy[11] ;
            sprintf( form_x,  "%dE", nEnergyRows ) ;
            sprintf( form_y,  "%dE", nOffsetRows ) ;
            sprintf( form_xy, "%dE", nelem_xy ) ;
            
            // define our table column names, types, and units
            char* ttype[] = { "ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "GAMMA", "SIGMA"  } ;
            char* tform[] = { form_x    , form_x    , form_y    , form_y    , form_xy, form_xy  } ;
            char* tunit[] = { "TeV"     , "TeV"     , "deg"     , "deg"     , "\0"   , "deg"    } ;
            // use "\0" for no unit
            
            // create table
            if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // setup variables for column numbers
            int C_ENERG_LO =  1 ;
            int C_ENERG_HI =  2 ;
            int C_THETA_LO =  3 ;
            int C_THETA_HI =  4 ;
            int C_GAMMA    =  5 ;
            int C_SIGMA    =  6 ;
            
            // set some columns are arrays, so we have to specify their dimensionality <nEnergyRows>x<nOffsetRows>
            int  naxis  = 2   ; // (XX,XX), I think naxis=3 would be (XX,XX,XX), etc.
            long naxes[naxis] ;
            naxes[0] = nEnergyRows ;
            naxes[1] = nOffsetRows ;
            if( fits_write_tdim( chunkptr, C_GAMMA  , naxis, naxes, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_tdim( chunkptr, C_SIGMA  , naxis, naxes, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // create data maps for each column
            float data_ENERG_LO[nEnergyRows] ;
            float data_ENERG_HI[nEnergyRows] ;
            float data_THETA_LO[nOffsetRows] ;
            float data_THETA_HI[nOffsetRows] ;
            float data_GAMMA   [nOffsetRows][nEnergyRows] ;
            float data_SIGMA   [nOffsetRows][nEnergyRows] ;
            bool  data_INTERP  [nOffsetRows][nEnergyRows] ;
            
            // set THETA column edges
            for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
            {
                data_THETA_LO[i_offset] = magOffset->offsetloweredges[i_offset] ;
                data_THETA_HI[i_offset] = magOffset->offsetupperedges[i_offset] ;
            }
            
            // set ENERG column edges
            for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
            {
                data_ENERG_LO[i_en] = pow( 10, en_logmin + ( i_en   * en_logdelta ) ) ;
                data_ENERG_HI[i_en] = pow( 10, en_logmin + ( ( i_en + 1 ) * en_logdelta ) ) ;
            }
            
            // fill GAMMA and SIGMA arrays
            int    i = 0 ;
            double min_i_diff = 100.0 ;
            double en_logcent = 0.0 ;
            double en_logtmp  = 0.0 ;
            double sigma = 0.0 ;
            double gamma = 0.0 ;
            double targ_offset = 0.0 ;
            double targ_en     = 0.0 ;
            for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
            {
                archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                printf( "\ni_offset=%d\n", i_offset ) ;
                
                for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                {
                    en_logcent = en_logmin + ( ( i_en + 0.5 ) * en_logdelta ) ;
                    //printf("  i_en=%2d  en_cent=%4.2f\n", i_en, en_logcent ) ;
                    i = 0 ;
                    min_i_diff = 100.0 ;
                    
                    // find energy bin closest to this en_logcent
                    for( unsigned int j_en = 0 ; j_en < archiveEffArea->fPsf[archiveIndex].energy.size() ; j_en++ )
                    {
                        en_logtmp = archiveEffArea->fPsf[archiveIndex].energy[j_en] ;
                        if( fabs( en_logcent - en_logtmp ) < min_i_diff )
                        {
                            min_i_diff = fabs( en_logcent - en_logtmp ) ;
                            i = j_en ;
                            //printf("    new min i=%2d  diff=%6.4f\n", i, min_i_diff ) ;
                        }
                    }
                    if( min_i_diff >= en_logdelta / 2.0 )
                    {
                        printf( "warning, not enough statistics in psf table at offset=%.2fdeg , log10(en)=%.2f\n", magOffset->offsets[i_offset], en_logcent );
                        gamma = -1.0 ;
                        sigma = -1.0 ;
                        data_INTERP[i_offset][i_en] = true ;
                    }
                    else
                    {
                        printf( "  i_en=%2d  using min i=%2d  diff=%6.4f  psftarg=%.2f  loading=%.2f\n", i_en, i, min_i_diff, en_logcent, archiveEffArea->fPsf[archiveIndex].energy[i] ) ;
                        
                        // load the sigma and gamma values from that closest energy bin
                        sigma = archiveEffArea->fPsf[archiveIndex].nativeking_sigma[i] ;
                        gamma = archiveEffArea->fPsf[archiveIndex].nativeking_gamma[i] ;
                        data_INTERP[i_offset][i_en] = false ;
                    }
                    
                    targ_offset =  0.5 ; // deg
                    targ_en     = 10.0 ; // TeV
                    if( ( offset > magOffset->offsetloweredges[i_offset] ) && ( offset < magOffset->offsetupperedges[i_offset] ) )
                    {
                        if( fabs( en_logcent - log10( targ_en ) ) < en_logdelta / 2.0 )
                        {
                            printf( "TARG EN=%.1f  OF=%.1f  SIGMA=%f  GAMMA=%f\n", pow( 10, en_logcent ), targ_offset, sigma, gamma ) ;
                            printf( "EN RANGE %f %f\n", pow( 10, en_logmin + ( i_en * en_logdelta ) ), pow( 10, en_logmax + ( ( i_en + 1 )*en_logdelta ) ) ) ;
                            printf( "OF RANGE %f %f\n", magOffset->offsetloweredges[i_offset], magOffset->offsetupperedges[i_offset] ) ;
                        }
                        
                    }
                    
                    // fill the psf table withour gamma and sigma
                    data_SIGMA[i_offset][i_en] = sigma ;
                    data_GAMMA[i_offset][i_en] = gamma ;
                }
            }
            
            // for all data_SIGMA/_GAMMA bins with data_INTERP set to 'true' ,
            // use the average value from un-interpolated neighboring bins
            int ndelta = 2 ; // scan all bins this far away for averaging
            for( int i = -ndelta ; i <= ndelta ; i++ )
            {
                printf( "scanning neighbors: %d\n", i ) ;
            }
            
            // do several passes, since once may not reach all bins
            for( int i_pass = 0 ; i_pass < 4 ; i_pass++ )
            {
            
                // for each table array index
                for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
                {
                    for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                    {
                    
                        // if it needs interpolation
                        if( data_INTERP[i_offset][i_en] )
                        {
                            float interp_gamma_tot = 0.0 ;
                            float interp_sigma_tot = 0.0 ;
                            int ninterps = 0 ;
                            int jof = 0 ;
                            int jen = 0 ;
                            
                            // loop over all neighboring bins
                            for( int dx = -ndelta ; dx <= ndelta ; dx++ )
                            {
                                for( int dy = -ndelta ; dy <= ndelta ; dy++ )
                                {
                                    jen = i_en + dx ;
                                    jof = i_offset + dy ;
                                    
                                    // make sure our indexes are still in the array, and are un-interpolated
                                    if( jof >= 0 && jof < nOffsetRows && jen >= 0 && jen < nEnergyRows )
                                    {
                                        // only use uninterpolated bins so far
                                        // this also excludes the target bin itself (i_offset,i_en) as well
                                        if( ! data_INTERP[jof][jen] )
                                        {
                                            interp_gamma_tot += data_GAMMA[jof][jen] ;
                                            interp_sigma_tot += data_SIGMA[jof][jen] ;
                                            ninterps += 1 ;
                                        }
                                    }
                                }
                            }
                            
                            // average all neighboring bins that met our criteria
                            if( ninterps > 0 && interp_gamma_tot > 0.0 && interp_sigma_tot )
                            {
                                sigma = interp_sigma_tot / ninterps ;
                                gamma = interp_gamma_tot / ninterps ;
                                printf( "(pass %d) successfully averaged values for iof=%2d ien=%2d (%2d neighbors) : sigma=%5.3f  gamma=%5.3f\n", i_pass, i_offset, i_en, ninterps, sigma, gamma ) ;
                                data_SIGMA[ i_offset][i_en] = sigma ;
                                data_GAMMA[ i_offset][i_en] = gamma ;
                                data_INTERP[i_offset][i_en] = false ;
                            }
                            
                        }
                    }
                }
            }
            
            // check that all table points have valid values
            int n_missing = 0 ;
            for( int iof = 0 ; iof < nOffsetRows ; iof++ )
            {
                for( int ien = 0 ; ien < nEnergyRows ; ien++ )
                {
                    if( data_INTERP[iof][ien] )
                    {
                        printf( "still missing a value: (iof,ien) : (%2d,%2d)\n", iof, ien ) ;
                        n_missing += 1 ;
                    }
                }
            }
            printf( "in total data_GAMMA/SIGMA are still missing %d table bin values...\n", n_missing ) ;
            
            // write columns and arrays to fits table
            if( fits_write_col( chunkptr, TFLOAT, C_ENERG_LO, 1, 1, nEnergyRows, &data_ENERG_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_ENERG_HI, 1, 1, nEnergyRows, &data_ENERG_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_THETA_LO, 1, 1, nOffsetRows, &data_THETA_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_THETA_HI, 1, 1, nOffsetRows, &data_THETA_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_GAMMA   , 1, 1, nelem_xy   , &data_GAMMA   , &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_SIGMA   , 1, 1, nelem_xy   , &data_SIGMA   , &status ) )
            {
                print_fits_error( status ) ;
            }
            
            if( fits_close_file( chunkptr, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // exit out of this psf-table-writing code block
            break ;
            
        } // endwhile
        
        writePSFtable = true ;
        while( writePSFtable )
        {
        
            printf( "staring buffered psf table\n" ) ;
            // open fits file for editing
            if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // setup table configuration
            sprintf( tblname, "%s", "POINT SPREAD FUNCTION" ) ;
            nrows   = 1 ;
            tfields = 6 ;
            
            // number of buffer rows to add to each end of the parameter space
            int nbuffer = 3 ;
            
            int archiveIndex = 0 ;
            int nOffsetRows  = magOffset->offsets.size() ;
            nOffsetRows += nbuffer ;
            
            // figure out energy binning (different offsets have different energy ranges)
            // since all offsets need to go into one table
            printf( "figure out energy binning\n" ) ;
            vector<double> e_bins ;
            for( unsigned int i_offset = 0 ; i_offset < magOffset->offsets.size() ; i_offset++ )
            {
                archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                
                // loop over each energy in this offset
                for( unsigned int i = 0 ; i < archiveEffArea->fPsf[archiveIndex].energy.size() ; i++ )
                {
                    e_bins.push_back( archiveEffArea->fPsf[archiveIndex].energy[ i ] ) ;
                }
            }
            
            // get sorted list of uniq energies
            sort( e_bins.begin(), e_bins.end() ) ;
            e_bins.erase( unique( e_bins.begin(), e_bins.end() ), e_bins.end() ) ;
            for( unsigned int i = 0 ; i < e_bins.size() ; i++ )
            {
                printf( "e_bins[%2d] = %.2f\n", i, e_bins[i] ) ;
            }
            
            // figure out smallest bin difference
            //printf("figure out energy bin widths\n") ;
            float en_logdelta = 1.0 ;
            float sld = 1.0 ;
            for( unsigned int i = 0 ; i < e_bins.size() - 1 ; i++ )
            {
                sld = e_bins[i + 1] - e_bins[i] ;
                if( sld < en_logdelta )
                {
                    en_logdelta = sld ;
                }
            }
            printf( "energy bins are %.2f log10(TeV) wide...\n", en_logdelta ) ;
            
            // add energy bins to the beginning until we at least cover -0.5
            double tmp_en = 5 ;
            int n_ebins = 0 ;
            double tmp_en_min = e_bins[0] ;
            //printf("\n tmp_en_min = %f\n", tmp_en_min ) ;
            //printf(" preq: %f\n", ( tmp_en_min - (-0.5) ) / en_logdelta ) ;
            n_ebins = ( int )( ( ( tmp_en_min - ( -0.5 ) ) / en_logdelta ) - 0.5 + nbuffer ) ;
            printf( "adding %d extra low-energy bins...", n_ebins );
            for( int i = 0 ; i < n_ebins ; i++ )
            {
                tmp_en = tmp_en_min - ( ( i + 1 ) * en_logdelta ) ;
                e_bins.insert( e_bins.begin(), tmp_en ) ;
                //printf("adding energy bin %.2f...\n", tmp_en ) ;
            }
            
            // add energy bins to the end until we at least cover 2.5
            double tmp_en_max = e_bins[ e_bins.size() - 1 ] ;
            //printf("\n preq: %f\n", ( 2.5 - tmp_en_max ) / en_logdelta ) ;
            n_ebins = ( int )( ( ( 2.5 - tmp_en_max ) / en_logdelta ) + 0.5 + nbuffer ) ;
            printf( "adding %d extra high-energy bins...", n_ebins );
            for( int i = 0 ; i < n_ebins ; i++ )
            {
                tmp_en = tmp_en_max + ( ( i + 1 ) * en_logdelta ) ;
                e_bins.insert( e_bins.end(), tmp_en ) ;
                //printf("adding energy bin %.2f...\n", tmp_en ) ;
                
            }
            
            //printf("\n");
            //for ( unsigned int i = 0 ; i < e_bins.size() ; i++ )
            //{
            //  printf("e_bins[%2d] = %5.2f\n", i, e_bins[i] ) ;
            //}
            
            double energyMin   = pow( 10, e_bins[0]               - ( en_logdelta / 2.0 ) ) ; // TeV
            double energyMax   = pow( 10, e_bins[e_bins.size() - 1] + ( en_logdelta / 2.0 ) ) ; // TeV
            double en_logmin   = log10( energyMin ) ;
            double en_logmax   = log10( energyMax ) ;
            int    nEnergyRows = e_bins.size() ;
            printf( "\n" );
            printf( "log energyMin : %.3f\n", en_logmin ) ;
            printf( "log energyMax : %.3f\n", en_logmax ) ;
            printf( "en_logdelta   : %.3f\n", en_logdelta ) ;
            printf( "nEnergyRows   : %d\n"  , nEnergyRows ) ;
            printf( "\n" );
            
            // print out our energy bins
            for( int i = 0 ; i < nEnergyRows ; i++ )
            {
                printf( "en bin %d : %.3f - %.3f (%7.3f - %7.3f TeV)\n", i, en_logmin + ( i * en_logdelta ), en_logmin + ( ( i + 1 )*en_logdelta ), pow( 10, en_logmin + ( i * en_logdelta ) ), pow( 10, en_logmin + ( ( i + 1 )*en_logdelta ) ) ) ;
            }
            printf( "\npsf table will have %d bins from %.3f to %.3f TeV...\n", nEnergyRows, energyMin, energyMax ) ;
            
            // how many elements should our PSF array columns have?
            int nelem_xy = nEnergyRows * nOffsetRows ; // number of elements in our 2D arrays
            
            // figure out column names, formats, and units
            char form_x [11] ;
            char form_y [11] ;
            char form_xy[11] ;
            sprintf( form_x,  "%dE", nEnergyRows ) ;
            sprintf( form_y,  "%dE", nOffsetRows ) ;
            sprintf( form_xy, "%dE", nelem_xy ) ;
            
            // define our table column names, types, and units
            char* ttype[] = { "ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "GAMMA", "SIGMA"  } ;
            char* tform[] = { form_x    , form_x    , form_y    , form_y    , form_xy, form_xy  } ;
            char* tunit[] = { "TeV"     , "TeV"     , "deg"     , "deg"     , "\0"   , "deg"    } ;
            // use "\0" for no unit
            
            // create table
            if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // setup variables for column numbers
            int C_ENERG_LO =  1 ;
            int C_ENERG_HI =  2 ;
            int C_THETA_LO =  3 ;
            int C_THETA_HI =  4 ;
            int C_GAMMA    =  5 ;
            int C_SIGMA    =  6 ;
            
            // set some columns are arrays, so we have to specify their dimensionality <nEnergyRows>x<nOffsetRows>
            int  naxis  = 2   ; // (XX,XX), I think naxis=3 would be (XX,XX,XX), etc.
            long naxes[naxis] ;
            naxes[0] = nEnergyRows ;
            naxes[1] = nOffsetRows ;
            if( fits_write_tdim( chunkptr, C_GAMMA  , naxis, naxes, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_tdim( chunkptr, C_SIGMA  , naxis, naxes, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // create data maps for each column
            float data_ENERG_LO[nEnergyRows] ;
            float data_ENERG_HI[nEnergyRows] ;
            float data_THETA_LO[nOffsetRows] ;
            float data_THETA_HI[nOffsetRows] ;
            float data_GAMMA   [nOffsetRows][nEnergyRows] ;
            float data_SIGMA   [nOffsetRows][nEnergyRows] ;
            bool  data_INTERP  [nOffsetRows][nEnergyRows] ;
            
            // initialize all data_INTERP elements to true
            for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
            {
                for( int j_energy = 0 ; j_energy < nEnergyRows ; j_energy++ )
                {
                    data_GAMMA [i_offset][j_energy] = 0.0  ;
                    data_SIGMA [i_offset][j_energy] = 0.0  ;
                    data_INTERP[i_offset][j_energy] = true ;
                }
            }
            
            // set THETA column edges
            for( unsigned int i_offset = 0 ; i_offset < magOffset->offsets.size() ; i_offset++ )
            {
                data_THETA_LO[i_offset] = magOffset->offsetloweredges[i_offset] ;
                data_THETA_HI[i_offset] = magOffset->offsetupperedges[i_offset] ;
            }
            
            // set ENERG column edges
            for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
            {
                data_ENERG_LO[i_en] = pow( 10, en_logmin + ( i_en   * en_logdelta ) ) ;
                data_ENERG_HI[i_en] = pow( 10, en_logmin + ( ( i_en + 1 ) * en_logdelta ) ) ;
            }
            
            // fill GAMMA and SIGMA arrays
            bool   found_match = false ;
            int    i = 0 ;
            double min_i_diff = 100.0 ;
            double en_logcent = 0.0 ;
            double en_logtmp  = 0.0 ;
            double sigma = 0.0 ;
            double gamma = 0.0 ;
            for( unsigned int i_offset = 0 ; i_offset < magOffset->offsets.size() ; i_offset++ )
            {
                archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                printf( "\ni_offset=%d\n", i_offset ) ;
                
                for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                {
                    en_logcent = en_logmin + ( ( i_en + 0.5 ) * en_logdelta ) ;
                    //printf("  i_en=%2d  en_cent=%4.2f\n", i_en, en_logcent ) ;
                    i = 0 ;
                    min_i_diff = 100.0 ;
                    found_match = false ;
                    sigma = 0.0 ;
                    gamma = 0.0 ;
                    
                    // find energy bin closest to this en_logcent
                    for( unsigned int j_en = 0 ; j_en < archiveEffArea->fPsf[archiveIndex].energy.size() ; j_en++ )
                    {
                        en_logtmp = archiveEffArea->fPsf[archiveIndex].energy[j_en] ;
                        if( fabs( en_logcent - en_logtmp ) < min_i_diff )
                        {
                            min_i_diff = fabs( en_logcent - en_logtmp ) ;
                            i = j_en ;
                            //printf("    new min i=%2d  diff=%6.4f\n", i, min_i_diff ) ;
                        }
                    }
                    if( min_i_diff < en_logdelta / 2.0 )
                    {
                        printf( "  i_en=%2d  using min i=%2d  diff=%6.4f  psftarg=%.2f  loading=%.2f\n", i_en, i, min_i_diff, en_logcent, archiveEffArea->fPsf[archiveIndex].energy[i] ) ;
                        
                        // load the sigma and gamma values from that closest energy bin
                        sigma = archiveEffArea->fPsf[archiveIndex].nativeking_sigma[i] ;
                        gamma = archiveEffArea->fPsf[archiveIndex].nativeking_gamma[i] ;
                        found_match = true ;
                    }
                    
                    if( found_match )
                    {
                        // fill the psf table withour gamma and sigma
                        data_SIGMA[ i_offset][i_en] = sigma ;
                        data_GAMMA[ i_offset][i_en] = gamma ;
                        data_INTERP[i_offset][i_en] = false ;
                    }
                }
            }
            
            // print block of which data_ elements have native values
            // 0 = has value
            // - = needs interpolation
            printf( "\ndata_INTERP     i_en:0123456789012345678901234567890123456789\n" );
            for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
            {
                printf( "data_INTERP i_off=%2d ", i_offset ) ;
                for( int j_energy = 0 ; j_energy < nEnergyRows ; j_energy++ )
                {
                    if( data_INTERP[i_offset][j_energy] )
                    {
                        printf( "-" );
                    }
                    else
                    {
                        printf( "O" );
                    }
                }
                printf( "\n" );
            }
            printf( "\n" );
            
            
            // write columns and arrays to fits table
            if( fits_write_col( chunkptr, TFLOAT, C_ENERG_LO, 1, 1, nEnergyRows, &data_ENERG_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_ENERG_HI, 1, 1, nEnergyRows, &data_ENERG_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_THETA_LO, 1, 1, nOffsetRows, &data_THETA_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_THETA_HI, 1, 1, nOffsetRows, &data_THETA_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_GAMMA   , 1, 1, nelem_xy   , &data_GAMMA   , &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_SIGMA   , 1, 1, nelem_xy   , &data_SIGMA   , &status ) )
            {
                print_fits_error( status ) ;
            }
            
            if( fits_close_file( chunkptr, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            break ;
            
        } // end while
        
        printf( "\nending psf table for chunk %d...\n\n\n", i_chunk ) ;
        
        
        ///////////////////////////////////////////////////
        // add EFFECTIVE AREA table
        // -effarea
        // move to while loop
        bool writeEFFECTIVEAREAtable = true ;
        if( skipeffarea )
        {
            writeEFFECTIVEAREAtable = false ;    // user may have flagged to skip the effective area part
        }
        while( writeEFFECTIVEAREAtable )
        {
            // What format table should be written?
            int  tablemode = EFFAREATABLE_BASIC ;
            if( tablemode == EFFAREATABLE_BASIC )
            {
                // open fits file
                if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup table info
                sprintf( tblname, "%s", "EFFECTIVE AREA" ) ;
                nrows   = 1 ;
                tfields = 6 ;
                
                // how many elements will we have in each axis?
                int nEnergyRows = 20 ; // manually chosen
                int nOffsetRows = magOffset->offsets.size() ;
                int nelem_xy    = nEnergyRows * nOffsetRows ; // number of elements in our 2D arrays
                
                float energyMin =   0.0125 ; // 12.5 GeV
                float energyMax = 500.0    ; // 200  TeV
                
                // figure out column names, formats, and units
                char form_x [11] ;
                char form_y [11] ;
                char form_xy[11] ;
                sprintf( form_x,  "%dE", nEnergyRows ) ;
                sprintf( form_y,  "%dE", nOffsetRows ) ;
                sprintf( form_xy, "%dE", nelem_xy ) ;
                // 21E = each row element in ENERG_LO and ENERG_HI will be a list of 21 floats('E')
                // 9E  = each row element in THETA_LO and THETA_HI will be a list if  9 floats('E')
                char* ttype[]   = { "ENERG_LO", "ENERG_HI", "THETA_LO", "THETA_HI", "EFFAREA", "EFFAREA_RECO" } ;
                char* tform[]   = { form_x    , form_x    , form_y    , form_y    , form_xy  , form_xy        } ;
                char* tunit[]   = { "TeV"     , "TeV"     , "deg"     , "deg"     , "m2"     , "m2"           } ;
                // use "\0" for no unit
                
                // create new table in fits file
                if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup variables for column numbers
                int  C_ENERG_LO     =  1 ;
                int  C_ENERG_HI     =  2 ;
                int  C_THETA_LO     =  3 ;
                int  C_THETA_HI     =  4 ;
                int  C_EFFAREA      =  5 ;
                int  C_EFFAREA_RECO =  6 ;
                
                // set some columns to have dimensionality (21,9)
                // what this means: each row element in the SCALE, SIGMA* and AMPL* columns
                // was set to be a list of 189 floats (189E) (in the ttform variable above), but below we
                // reorganize these 189 floats into a 2D array of 21x9 floats (21*9=189)
                int  naxis  = 2   ; // (XX,XX), I think naxis=3 would be (XX,XX,XX), etc.
                long naxes[naxis] ;
                naxes[0] = nEnergyRows ; // effareatag
                naxes[1] = nOffsetRows ; // effareatag
                //naxes[0] = nOffsetRows ; // effareatag
                //naxes[1] = nEnergyRows ; // effareatag
                if( fits_write_tdim( chunkptr, C_EFFAREA     , naxis, naxes, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_tdim( chunkptr, C_EFFAREA_RECO, naxis, naxes, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // create data maps for each column
                float data_ENERG_LO    [nEnergyRows] ;
                float data_ENERG_HI    [nEnergyRows] ;
                float data_THETA_LO    [nOffsetRows] ;
                float data_THETA_HI    [nOffsetRows] ;
                //float data_EFFAREA     [nEnergyRows][nOffsetRows] ; // effareatag
                //float data_EFFAREA_RECO[nEnergyRows][nOffsetRows] ; // effareatag
                float data_EFFAREA     [nOffsetRows][nEnergyRows] ; // effareatag
                float data_EFFAREA_RECO[nOffsetRows][nEnergyRows] ; // effareatag
                
                // loop over each offset
                for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
                {
                    // offset bin edges
                    data_THETA_LO[i_offset] = magOffset->offsetloweredges[i_offset] ;
                    data_THETA_HI[i_offset] = magOffset->offsetupperedges[i_offset] ;
                    
                    // create empty histograms
                    char hname1[200], hname2[200], hdesc1[200], hdesc2[200] ;
                    sprintf( hname1, "effarea_chunk%d_offset%d_vsMC"  , i_chunk, i_offset ) ;
                    sprintf( hname2, "effarea_chunk%d_offset%d_vsRECO", i_chunk, i_offset ) ;
                    sprintf( hdesc1, "effective areas vs mc energy, for offset between %.2fdeg and %.2fdeg"                  , data_THETA_LO[i_offset], data_THETA_HI[i_offset] );
                    sprintf( hdesc2, "effective areas vs true (reconstructed) energy, for offset between %.2fdeg and %.2fdeg", data_THETA_LO[i_offset], data_THETA_HI[i_offset] );
                    TH1F* hist_mc = new TH1F( hname1, hdesc1, nEnergyRows, log10( energyMin ), log10( energyMax ) ) ;
                    TH1F* hist_tr = new TH1F( hname2, hdesc2, nEnergyRows, log10( energyMin ), log10( energyMax ) ) ;
                    hist_tr->SetXTitle( "log_{10} (E_{True}/TeV)" ) ;
                    hist_mc->SetXTitle( "log_{10} (E_{MC}/TeV)" ) ;
                    hist_tr->SetYTitle( "effective area (m^{2})" ) ;
                    hist_mc->SetYTitle( "effective area (m^{2})" ) ;
                    hist_mc->SetLineColor( 2 ) ;
                    
                    // get this offset's archive
                    if( nOffsetRows != ( int ) effAreaChunk.size() )
                    {
                        printf( "Error, size mismatch between nOffsetRows and effAreaChunk.size(), exiting...\n" );
                        exit( 1 ) ;
                    }
                    int archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                    
                    // fill the histograms
                    archiveEffArea->irf[archiveIndex]->fillEffectiveAreasHistograms( hist_tr, "", hist_mc ) ;
                    if( debugRootFile )
                    {
                        debugRootFile->cd();
                        hist_mc->Write();
                        hist_tr->Write();
                    }
                    
                    // loop over energy bins (log space)
                    for( int i_en = 0 ; i_en < nEnergyRows ; i_en++ )
                    {
                        // store this bin's energy range
                        data_ENERG_LO[i_en] = pow( 10, hist_tr->GetXaxis()->GetBinLowEdge( i_en + 1 ) ) ;
                        data_ENERG_HI[i_en] = pow( 10, hist_tr->GetXaxis()->GetBinUpEdge( i_en + 1 ) ) ;
                        
                        // store effective area data to tables
                        //data_EFFAREA     [i_en][i_offset] = hist_mc->GetBinContent( i_en+1 ) ; // effareatag
                        //data_EFFAREA_RECO[i_en][i_offset] = hist_tr->GetBinContent( i_en+1 ) ; // effareatag
                        data_EFFAREA     [i_offset][i_en] = hist_mc->GetBinContent( i_en + 1 ) ; // effareatag
                        data_EFFAREA_RECO[i_offset][i_en] = hist_tr->GetBinContent( i_en + 1 ) ; // effareatag
                        //printf("  chunk(%2d)  offset(%4.2f:%4.2f)  energy(%5.1f:%5.1f) : effarea :  mc=%9.2f  reco=%9.2f\n", i_chunk+1,
                        //  data_THETA_LO           [i_offset],
                        //  data_THETA_HI           [i_offset],
                        //  log10(data_ENERG_LO[i_en])        ,
                        //  log10(data_ENERG_HI[i_en])        ,
                        //  data_EFFAREA     [i_offset] [i_en], // effareatag
                        //  data_EFFAREA_RECO[i_offset] [i_en] ) ; // effareatag
                        //data_EFFAREA     [i_en] [i_offset], // effareatag
                        //data_EFFAREA_RECO[i_en] [i_offset] ) ; // effareatag
                    }
                }
                
                // write data to columns
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_LO     , 1, 1, nEnergyRows, &data_ENERG_LO    , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_HI     , 1, 1, nEnergyRows, &data_ENERG_HI    , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_THETA_LO     , 1, 1, nOffsetRows, &data_THETA_LO    , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_THETA_HI     , 1, 1, nOffsetRows, &data_THETA_HI    , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_EFFAREA      , 1, 1, nelem_xy   , &data_EFFAREA     , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_EFFAREA_RECO , 1, 1, nelem_xy   , &data_EFFAREA_RECO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // close fits file
                if( fits_close_file( chunkptr, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
            }
            else
            {
                printf( "Error, unrecognized table mode, exiting...\n" );
                exit( 1 );
            }
            
            // break out of while loop
            break ;
        } // endwhile
        
        
        /////////////////
        // add proper ENERGY RESOLUTION table
        // -eres
        // to be read in by gammalib.GCTAEdisp2D()
        bool writeEDISPtable = true ;
        while( writeEDISPtable )
        {
        
            // open fits file for editing
            if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
            {
                print_fits_error( status ) ;
            }
            sprintf( tblname, "%s", "ENERGY DISPERSION" ) ;
            
            int nOffsetRows = magOffset->offsets.size() ;
            
            // figure out energy migration matrix energy rows
            //char hname1[200] ;
            //int archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[ 0] ) ; // look at first offset's psf vector for the energy bins
            int archiveIndex = 0 ;
            
            //TH2D * migmatrix = (TH2D*) archiveEffArea->fenmigmatrix[archiveIndex]->Clone("migmatrix_ioffset_0") ;
            //TH2D * migmatrix = archiveEffArea->fenmigmatrix[archiveIndex] ;
            TH2D* migmatrix = 0 ;
            
            //printf("migra_\n");
            float migra_min = 10000.0 ;
            float migra_max = 0.00001 ;
            
            // the MIGRA column in GCTAEdisp2D is a ratio of E_reconstructed / E_true,
            // so we have to find the min and max ratios to be used with this table
            int   ibin  = 0 ;
            float binc  = 0.0 ;
            float elo   = 0.0 ;
            float ehi   = 0.0 ;
            float etrue = 0.0 ;
            float ereco = 0.0 ;
            float ratio = 0.0 ;
            for( int i_offset = 0 ; i_offset < nOffsetRows ; i_offset++ )
            {
                archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[ i_offset] ) ; // look at first offset's psf vector for the energy bins
                migmatrix = archiveEffArea->fenmigmatrix[archiveIndex] ;
                for( int i_etrue = 0 ; i_etrue < migmatrix->GetNbinsY() ; i_etrue++ )
                {
                    for( int i_ereco = 0 ; i_ereco < migmatrix->GetNbinsX() ; i_ereco++ )
                    {
                        ibin = migmatrix->GetBin( i_ereco, i_etrue ) ;
                        binc = migmatrix->GetBinContent( ibin ) ;
                        if( binc > 0 )
                        {
                            elo = migmatrix->GetYaxis()->GetBinLowEdge( i_etrue ) ;
                            ehi = migmatrix->GetYaxis()->GetBinUpEdge( i_etrue ) ;
                            etrue = pow( 10, ( ehi + elo ) / 2.0 ) ;
                            elo = migmatrix->GetXaxis()->GetBinLowEdge( i_ereco ) ;
                            ehi = migmatrix->GetXaxis()->GetBinUpEdge( i_ereco ) ;
                            ereco = pow( 10, ( ehi + elo ) / 2.0 ) ;
                            ratio = ereco / etrue ;
                            //if ( i_chunk == 0 )
                            //{
                            //printf("migra value %f\n", ratio ) ;
                            //}
                            if( ratio < migra_min )
                            {
                                //printf("new migra_min at (%3d,%3d): ereco/etrue = %f / %f = %f \n", i_ereco, i_etrue, ereco, etrue, ratio ) ;
                                migra_min = ratio ;
                            }
                            if( ratio > migra_max )
                            {
                                //printf("new migra_max at (%3d,%3d): ereco/etrue = %f / %f = %f \n", i_ereco, i_etrue, ereco, etrue, ratio ) ;
                                migra_max = ratio ;
                            }
                        }
                    }
                }
            }
            
            // check if something went wrong
            if( migra_min >= migra_max )
            {
                printf( "warning, migra_min (%.2f) is >= migra_max (%.2f), which means something went wrong!\n", migra_min, migra_max ) ;
                exit( 1 );
            }
            
            // expand our bin range just a little to get the edge events
            migra_min *= 0.98 ;
            migra_max *= 1.02 ;
            //printf( "migra_min: %f\n", migra_min ) ;
            //printf( "migra_max: %f\n", migra_max ) ;
            
            archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[ 0] ) ; // look at first offset's psf vector for the energy bins
            migmatrix    = archiveEffArea->fenmigmatrix[archiveIndex] ;
            
            int nMigraRows  = 200 ;
            int nEtrueRows  = migmatrix->GetYaxis()->GetNbins() ; // energy_{MC}
            
            float migra_delta = ( log10( migra_max ) - log10( migra_min ) ) / nMigraRows ;
            //printf("migra_delta : %f\n", migra_delta ) ;
            
            //printf("nEtrueRows : %d\n", nEtrueRows  ) ;
            //printf("nEmigrRows : %d\n", nMigraRows  ) ;
            //printf("nOffsetRows: %d\n", nOffsetRows ) ;
            
            // set the fits table column lengths
            // if you swap these, you also have to change the
            // column loops below!
            int nelem_x   = nEtrueRows ;
            int nelem_y   = nMigraRows ;
            int nelem_z   = nOffsetRows ;
            int nelem_xyz = nelem_x * nelem_y * nelem_z ;
            
            tfields = 7 ; // fits table has this many columns
            nrows   = 1 ;
            
            // figure out column names, formats, and units
            char form_x  [11] ;
            char form_y  [11] ;
            char form_z  [11] ;
            char form_xyz[11] ;
            sprintf( form_x,   "%dE", nelem_x ) ;
            sprintf( form_y,   "%dE", nelem_y ) ;
            sprintf( form_z,   "%dE", nelem_z ) ;
            sprintf( form_xyz, "%dE", nelem_xyz ) ;
            char* ttype[]   = { "ETRUE_LO", "ETRUE_HI", "MIGRA_LO", "MIGRA_HI", "THETA_LO", "THETA_HI", "MATRIX" } ;
            char* tform[]   = { form_x    , form_x    , form_y    , form_y    , form_z    , form_z    , form_xyz } ;
            char* tunit[]   = { "TeV"     , "TeV"     , "\0"     , "\0"     , "deg"     , "deg"     , "counts" } ;
            // use "\0" for no unit
            
            // create table
            if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // setup variables for column numbers
            int  C_ETRUE_LO =  1 ;
            int  C_ETRUE_HI =  2 ;
            int  C_MIGRA_LO =  3 ;
            int  C_MIGRA_HI =  4 ;
            int  C_THETA_LO =  5 ;
            int  C_THETA_HI =  6 ;
            int  C_MATRIX   =  7 ;
            
            // set the format of the MATRIX table
            int  naxis  = 3   ; // (XX,XX), I think naxis=3 would be (XX,XX,XX), etc.
            long naxes[naxis] ;
            naxes[0] = nelem_x ; // ETRUE axis
            naxes[1] = nelem_y ; // MIGRA axis
            naxes[2] = nelem_z ; // THETA axis
            if( fits_write_tdim( chunkptr, C_MATRIX , naxis, naxes, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            float data_ETRUE_LO[ nelem_x ] ;
            float data_ETRUE_HI[ nelem_x ] ;
            float data_MIGRA_LO[ nelem_y ] ;
            float data_MIGRA_HI[ nelem_y ] ;
            float data_THETA_LO[ nelem_z ] ;
            float data_THETA_HI[ nelem_z ] ;
            float data_MATRIX  [ nelem_z ][ nelem_y ][ nelem_x ] ;
            
            // since we add ints to data_MATRIX, make sure its all zeros first
            for( int iz = 0 ; iz < nelem_z ; iz++ )
            {
                for( int iy = 0 ; iy < nelem_y ; iy++ )
                {
                    for( int ix = 0 ; ix < nelem_x ; ix++ )
                    {
                        data_MATRIX[iz][iy][ix] = 0.0 ;
                    }
                }
            }
            
            ///////////////////////////
            // SET ETRUE COLUMN BINS //
            ///////////////////////////
            //printf("\nFITS Columns:\n");
            // check for axis consistancy
            if( nelem_x > migmatrix->GetNbinsY() )
            {
                printf( "error, fits table ETRUE column set to use %d bins, but migmatrix only has %d bins, exiting...\n", nelem_x, migmatrix->GetNbinsY() ) ;
                exit( 1 );
            }
            
            // set the ETRUE column bin edge values
            for( int ix = 0 ; ix < nelem_x ; ix++ )
            {
                // LO = lower edge of the energy bin for this node
                // HI = upper edge of the energy bin for this node
                // migmatrix Y axis = true energy axis
                data_ETRUE_LO[ix] = pow( 10, migmatrix->GetYaxis()->GetBinLowEdge( ix ) ) ;
                data_ETRUE_HI[ix] = pow( 10, migmatrix->GetYaxis()->GetBinUpEdge( ix ) ) ;
                //printf("ETRUE : %.3f - %.3f TeV\n", data_ETRUE_LO[ix], data_ETRUE_HI[ix] ) ;
            }
            
            ///////////////////////////
            // SET MIGRA COLUMN BINS //
            ///////////////////////////
            for( int ix = 0 ; ix < nelem_y ; ix++ )
            {
                // migmatrix X axis = reconstructed energy axis
                data_MIGRA_LO[ix] = pow( 10, log10( migra_min ) + ( ix   * migra_delta ) ) ;
                data_MIGRA_HI[ix] = pow( 10, log10( migra_min ) + ( ( ix + 1 ) * migra_delta ) ) ;
                //printf("MIGRA : %.3f - %.3f (E_reco/E_true)\n", data_MIGRA_LO[ix], data_MIGRA_HI[ix] ) ;
            }
            
            ///////////////////////////
            // SET THETA COLUMN BINS //
            ///////////////////////////
            for( int i_offset = 0 ; i_offset < nelem_z ; i_offset++ )
            {
                data_THETA_LO[i_offset] = magOffset->offsetloweredges[i_offset] ;
                data_THETA_HI[i_offset] = magOffset->offsetupperedges[i_offset] ;
                //printf("THETA: %.1f - %.1f deg\n", data_THETA_LO[i_offset], data_THETA_HI[i_offset] ) ;
            }
            
            ////////////////////////////
            // SET MATRIX CUBE VALUES //
            ////////////////////////////
            
            // loop over each theta bin
            for( int i_offset = 0 ; i_offset < nelem_z ; i_offset++ )
            {
            
                // load migration matrix from effarea file
                int archiveIndex = archiveEffArea->checkIfLoaded( effAreaChunk[i_offset] ) ;
                
                int ibin = 0 ;
                char hname4[200] ;
                sprintf( hname4, "migmatrix_chunk%d_offset%d_archive%d", i_chunk, i_offset, archiveIndex ) ;
                //printf("loading migration matrix %-35s %-s\n", hname4, archiveEffArea->fenmigmatrix[archiveIndex]->GetTitle() ) ;
                if( debugRootFile )
                {
                    debugRootFile->cd();
                    TH2D* migmatrix = ( TH2D* ) archiveEffArea->fenmigmatrix[archiveIndex]->Clone( hname4 ) ;
                    //ibin = migmatrix->GetBin(100,50) ;
                    //printf("get bin content: %d\n", migmatrix->GetBinContent( ibin ) ) ;
                    migmatrix->Write();
                }
                //printf("after saving migmatrix to debug root file...\n");
                
                // loop over each E_True and E_Reco bin
                int ix = 0 ;
                int iy = 0 ;
                for( int i_etrue = 0 ; i_etrue < migmatrix->GetNbinsY() ; i_etrue++ )
                {
                    //printf("i_etrue = %d\n", i_etrue );
                    for( int i_ereco = 0 ; i_ereco < migmatrix->GetNbinsX() ; i_ereco++ )
                    {
                        //printf("i_ereco = %d\n", i_ereco );
                        
                        // only try to add this bin to MATRIX if it at least one event in it
                        ibin = archiveEffArea->fenmigmatrix[archiveIndex]->GetBin( i_ereco, i_etrue ) ;
                        binc = archiveEffArea->fenmigmatrix[archiveIndex]->GetBinContent( ibin ) ;
                        if( binc > 0 )
                        {
                            // figure out the E_Reco/E_True ratio of this energy bin
                            elo   = archiveEffArea->fenmigmatrix[archiveIndex]->GetYaxis()->GetBinLowEdge( i_etrue ) ;
                            ehi   = archiveEffArea->fenmigmatrix[archiveIndex]->GetYaxis()->GetBinUpEdge( i_etrue ) ;
                            etrue = pow( 10, ( ehi + elo ) / 2.0 ) ;
                            elo   = archiveEffArea->fenmigmatrix[archiveIndex]->GetXaxis()->GetBinLowEdge( i_ereco ) ;
                            ehi   = archiveEffArea->fenmigmatrix[archiveIndex]->GetXaxis()->GetBinUpEdge( i_ereco ) ;
                            ereco = pow( 10, ( ehi + elo ) / 2.0 ) ;
                            ratio = ereco / etrue ;
                            
                            // add our counts to the ETRUE/MIGRA bin
                            ix = i_etrue ;
                            iy = floor( ( log10( ratio ) - log10( migra_min ) ) / migra_delta ) ;
                            //printf("chunk %d  offset %d : etrue=%.3f  i_etrue=%d  ereco=%.3f  ratio=%.3f  iy=%d (nelem_y=%d)\n", i_chunk, i_offset, etrue, i_etrue, ereco, ratio, iy, nelem_y );
                            data_MATRIX[i_offset][iy][ix] += binc  ;
                            
                        } // endif : bincontent is > 0
                    } // endfor : i_reco loop
                } // endfor : i_etrue loop
                //printf("after writing this offset's data_MATRIX...\n");
            } // endfor : offset loop
            //printf("done filling data_MATRIX...\n");
            
            // write data to fits columns
            if( fits_write_col( chunkptr, TFLOAT, C_ETRUE_LO, 1, 1, nelem_x  , &data_ETRUE_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_ETRUE_HI, 1, 1, nelem_x  , &data_ETRUE_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_MIGRA_LO, 1, 1, nelem_y  , &data_MIGRA_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_MIGRA_HI, 1, 1, nelem_y  , &data_MIGRA_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_THETA_LO, 1, 1, nelem_z  , &data_THETA_LO, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_THETA_HI, 1, 1, nelem_z  , &data_THETA_HI, &status ) )
            {
                print_fits_error( status ) ;
            }
            if( fits_write_col( chunkptr, TFLOAT, C_MATRIX  , 1, 1, nelem_xyz, &data_MATRIX  , &status ) )
            {
                print_fits_error( status ) ;
            }
            
            // close fits file
            if( fits_close_file( chunkptr, &status ) )
            {
                print_fits_error( status ) ;
            }
            
            //printf("done writing fits table...\n");
            
            // break out of while loop
            break ;
        }
        
        
        ///////////////////////////////////////////////////
        // add BACKGROUND table
        // -bkgnd
        // move to a single while-loop, so variables can be redeclared when
        // adding each IRF, and code can just be copied between IRF sections
        bool writeBACKGROUNDtable = true ;
        while( writeBACKGROUNDtable )
        {
            int tablemode = BACKGROUNDTABLE_3D ;
            int tableflavor = BCK_FLAT ;
            
            // BCK_FLAT settings, only shape is used, not raw counts
            double flat_value = 1e1 ;
            // BCK_OFFSET settings, only shape is used, not actual counts
            double offset_peak  = 1e-7   ; // maximum counts at offset=0 (1/(s*MeV*sr))
            double offset_zero  = 10e-12 ; // minimum counts at offset=infinity (1/(s*MeV*sr))
            double offset_sigma = 0.75   ; // gaussian sigma radius (deg)
            
            if( tablemode == BACKGROUNDTABLE_BASIC )
            {
            
                // open fits file for editing
                if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup table columns
                sprintf( tblname, "%s", "BACKGROUND" ) ;
                nrows     = 1 ; // 1 row
                tfields   = 8 ; // 8 columns
                
                // how many elements should our background table have?
                int nelem_x = 18 ;  // x axis: DETX, or camera coordinate x axis
                int nelem_y = nelem_x ;  // y axis: DETY, or camera coordinate y axis
                int nelem_z = 21 ;  // z axis: energy
                int nelem_xyz = nelem_x * nelem_y * nelem_z ;
                if( nelem_x != nelem_y )
                {
                    cout << "Error, for BACKGROUND table, nelem_x=" << nelem_x << " and nelem_y=" << nelem_y << " must be equal, exiting..." << endl;
                    exit( 1 ) ;
                }
                
                // camera coordinate min and max
                float cam_max = 2.5 ; // distance from camera center to side of square, deg
                float cam_wid = 2 * cam_max / nelem_x ;
                
                // energy min and max
                float en_min =   0.0125 ; // 12.5 GeV
                float en_max = 200.0    ; // 200  TeV
                float en_logwid = ( log10( en_max ) - log10( en_min ) ) / nelem_x ;
                
                // figure out column names, formats, and units
                char form_x  [11] ;
                char form_y  [11] ;
                char form_z  [11] ;
                char form_xyz[11] ;
                sprintf( form_x,   "%dE", nelem_x ) ;
                sprintf( form_y,   "%dE", nelem_y ) ;
                sprintf( form_z,   "%dE", nelem_z ) ;
                sprintf( form_xyz, "%dE", nelem_xyz ) ;
                char* ttype[]   = { "DETX_LO", "DETX_HI", "DETY_LO", "DETY_HI", "ENERG_LO", "ENERG_HI", "BGD"       , "BGD_RECO"   } ;
                char* tform[]   = { form_x   , form_x   , form_y   , form_y   , form_z    , form_z    , form_xyz    , form_xyz     } ;
                char* tunit[]   = { "deg"    , "deg"    , "deg"    , "deg"    , "TeV"     , "TeV"     , "1/s/MeV/sr", "1/s/MeV/sr" } ;
                // use "\0" for no unit
                
                // create table
                if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup variables for column numbers
                int C_DETX_LO  = 1 ;
                int C_DETX_HI  = 2 ;
                int C_DETY_LO  = 3 ;
                int C_DETY_HI  = 4 ;
                int C_ENERG_LO = 5 ;
                int C_ENERG_HI = 6 ;
                int C_BGD      = 7 ;
                int C_BGD_RECO = 8 ;
                
                // set some columns to have dimensionality (18,18,21)
                // what this means: each row element in the BGD and BGD_RECO columns
                // was set to be a list of 6804 floats (6804E in the ttform variable above), but below we
                // reorganize these 6804 floats into a 3D array of 18x18x21 float (18*18*21=6804)
                int  naxis  = 3   ; // (XX,XX,XX)
                long naxes[naxis] ;
                naxes[0] = nelem_x ;
                naxes[1] = nelem_y ;
                naxes[2] = nelem_z ;
                if( fits_write_tdim( chunkptr, C_BGD     , naxis, naxes, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_tdim( chunkptr, C_BGD_RECO, naxis, naxes, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // create data maps for each column
                float data_DETX_LO [nelem_x] ;
                float data_DETX_HI [nelem_x] ;
                float data_DETY_LO [nelem_y] ;
                float data_DETY_HI [nelem_y] ;
                float data_ENERG_LO[nelem_z] ;
                float data_ENERG_HI[nelem_z] ;
                float data_BGD     [nelem_x][nelem_y][nelem_z] ;
                float data_BGD_RECO[nelem_x][nelem_y][nelem_z] ;
                
                // loop over each bin in energy and camera coord axis
                // and calculate the background params for that bin
                for( int ix = 0 ; ix < nelem_x ; ix++ )
                {
                    // calculate camera x bin edges
                    data_DETX_LO[ix] = -1 * cam_max + ( ix    * cam_wid ) ;
                    data_DETX_HI[ix] = -1 * cam_max + ( ( ix + 1 ) * cam_wid ) ;
                    
                    // loop over camera y bins
                    for( int iy = 0 ; iy < nelem_y ; iy++ )
                    {
                        // calculate camera y bin edges
                        data_DETY_LO[iy] = -1 * cam_max + ( iy    * cam_wid ) ;
                        data_DETY_HI[iy] = -1 * cam_max + ( ( iy + 1 ) * cam_wid ) ;
                        
                        // loop over energy bins
                        for( int iz = 0 ; iz < nelem_z ; iz++ )
                        {
                            // calculate energy bin edges
                            data_ENERG_LO[iz] = pow( 10, log10( en_min ) + ( iz    * en_logwid ) ) ;
                            data_ENERG_HI[iz] = pow( 10, log10( en_min ) + ( ( iz + 1 ) * en_logwid ) ) ;
                            
                            // setup structs to calculate background
                            inputbackground  back_inp ;
                            outputbackground back_out ;
                            back_inp.camx_min = data_DETX_LO [ix] ;
                            back_inp.camx_max = data_DETX_HI [ix] ;
                            back_inp.camy_min = data_DETY_LO [iy] ;
                            back_inp.camy_max = data_DETY_HI [iy] ;
                            back_inp.emin     = data_ENERG_LO[iz] ;
                            back_inp.emax     = data_ENERG_HI[iz] ;
                            
                            // calclate background values
                            calcBACKGROUNDvalues( back_inp, back_out ) ;
                            
                            // save valuse to data map
                            data_BGD     [ix][iy][iz] = back_out.bgd      ;
                            data_BGD_RECO[ix][iy][iz] = back_out.bgd_reco ;
                            
                        } // endfor: loop over energy bins
                    } // endfor: loop over camera y bins
                } // endfor: loop over camera x bins
                
                // write data to columns
                if( fits_write_col( chunkptr, TFLOAT, C_DETX_LO , 1, 1, nelem_x  , &data_DETX_LO , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_DETX_HI , 1, 1, nelem_x  , &data_DETX_HI , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_DETY_LO , 1, 1, nelem_y  , &data_DETY_LO , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_DETY_HI , 1, 1, nelem_y  , &data_DETY_HI , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_LO, 1, 1, nelem_z  , &data_ENERG_LO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_HI, 1, 1, nelem_z  , &data_ENERG_HI, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_BGD     , 1, 1, nelem_xyz, &data_BGD     , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_BGD_RECO, 1, 1, nelem_xyz, &data_BGD_RECO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // close fits file
                if( fits_close_file( chunkptr, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
            }
            else if( tablemode == BACKGROUNDTABLE_3D )
            {
            
                // open fits file for editing
                if( fits_open_file( &chunkptr, fitsFileName.c_str(), READWRITE, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup table columns
                sprintf( tblname, "%s", "BACKGROUND" ) ;
                nrows     = 1 ; // 1 row
                tfields   = 8 ; // 8 columns
                
                // how many bins in each axis in camera coordinates should we have?
                int nDETXYrows  = 18 ;
                
                // how many log10(energy) bins should we have?
                int nENERGYrows = 20 ;
                
                // how many elements should our background table have?
                int nelem_x = nDETXYrows ; // x axis: DETX, or camera coordinate x axis
                int nelem_y = nDETXYrows ; // y axis: DETY, or camera coordinate y axis
                int nelem_z = nENERGYrows; // z axis: energy
                int nelem_xyz = nelem_x * nelem_y * nelem_z ;
                if( nelem_x != nelem_y )
                {
                    cout << "Error, for BACKGROUND table, nelem_x=" << nelem_x << " and nelem_y=" << nelem_y << " must be equal, exiting..." << endl;
                    exit( 1 ) ;
                }
                
                // camera coordinate min and max
                float cam_max = 2.5 ; // distance from camera center to side of square, deg
                float cam_wid = 2 * cam_max / nelem_x ; // camera coordinate bin width
                
                // energy min and max
                float en_min =   0.0125 ; // 12.5 GeV
                float en_max = 200.0    ; // 200  TeV
                float en_logwid = ( log10( en_max ) - log10( en_min ) ) / nelem_x ;
                
                // figure out column names, formats, and units
                char form_x  [11] ;
                char form_y  [11] ;
                char form_z  [11] ;
                char form_xyz[11] ;
                sprintf( form_x,   "%dE", nelem_x ) ;
                sprintf( form_y,   "%dE", nelem_y ) ;
                sprintf( form_z,   "%dE", nelem_z ) ;
                sprintf( form_xyz, "%dE", nelem_xyz ) ;
                char* ttype[]   = { "DETX_LO", "DETX_HI", "DETY_LO", "DETY_HI", "ENERG_LO", "ENERG_HI", "BGD"       , "BGD_RECO"   } ;
                char* tform[]   = { form_x   , form_x   , form_y   , form_y   , form_z    , form_z    , form_xyz    , form_xyz     } ;
                char* tunit[]   = { "deg"    , "deg"    , "deg"    , "deg"    , "TeV"     , "TeV"     , "1/s/MeV/sr", "1/s/MeV/sr" } ;
                // use "\0" for no unit
                
                // create table
                if( fits_create_tbl( chunkptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, tblname, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // setup variables for column numbers
                int C_DETX_LO  = 1 ;
                int C_DETX_HI  = 2 ;
                int C_DETY_LO  = 3 ;
                int C_DETY_HI  = 4 ;
                int C_ENERG_LO = 5 ;
                int C_ENERG_HI = 6 ;
                int C_BGD      = 7 ;
                int C_BGD_RECO = 8 ;
                
                // set some columns to have dimensionality (18,18,21)
                // what this means: each row element in the BGD and BGD_RECO columns
                // was set to be a list of 6804 floats (6804E in the ttform variable above), but below we
                // reorganize these 6804 floats into a 3D array of 18x18x21 float (18*18*21=6804)
                int  naxis  = 3   ; // (XX,XX,XX)
                long naxes[naxis] ;
                naxes[0] = nelem_x ;
                naxes[1] = nelem_y ;
                naxes[2] = nelem_z ;
                if( fits_write_tdim( chunkptr, C_BGD     , naxis, naxes, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_tdim( chunkptr, C_BGD_RECO, naxis, naxes, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // create data maps for each column
                float data_DETX_LO [nelem_x] ;
                float data_DETX_HI [nelem_x] ;
                float data_DETY_LO [nelem_y] ;
                float data_DETY_HI [nelem_y] ;
                float data_ENERG_LO[nelem_z] ;
                float data_ENERG_HI[nelem_z] ;
                float data_BGD     [nelem_x][nelem_y][nelem_z] ;
                float data_BGD_RECO[nelem_x][nelem_y][nelem_z] ;
                
                // extra in-loop variables
                double detx = 0.0 ;
                double dety = 0.0 ;
                double offset = 0.0 ;
                
                // loop over each bin in energy and camera coord axis
                // and calculate the background params for that bin
                for( int ix = 0 ; ix < nelem_x ; ix++ )
                {
                    // calculate camera x bin edges
                    data_DETX_LO[ix] = -1 * cam_max + ( ix    * cam_wid ) ;
                    data_DETX_HI[ix] = -1 * cam_max + ( ( ix + 1 ) * cam_wid ) ;
                    
                    // loop over camera y bins
                    for( int iy = 0 ; iy < nelem_y ; iy++ )
                    {
                        // calculate camera y bin edges
                        data_DETY_LO[iy] = -1 * cam_max + ( iy    * cam_wid ) ;
                        data_DETY_HI[iy] = -1 * cam_max + ( ( iy + 1 ) * cam_wid ) ;
                        
                        // loop over energy bins
                        for( int iz = 0 ; iz < nelem_z ; iz++ )
                        {
                        
                            // calculate energy bin edges
                            data_ENERG_LO[iz] = pow( 10, log10( en_min ) + ( iz    * en_logwid ) ) ;
                            data_ENERG_HI[iz] = pow( 10, log10( en_min ) + ( ( iz + 1 ) * en_logwid ) ) ;
                            
                            // pick background flavor (shape of the background)
                            if( tableflavor == BCK_FLAT )
                            {
                                // save valuse to data map
                                data_BGD     [ix][iy][iz] = flat_value ;
                                data_BGD_RECO[ix][iy][iz] = flat_value ;
                            }
                            else if( tableflavor == BCK_OFFSET )
                            {
                                detx   = ( data_DETX_HI[ix] + data_DETX_LO[ix] ) / 2.0 ;
                                dety   = ( data_DETY_HI[ix] + data_DETY_LO[ix] ) / 2.0 ;
                                offset = sqrt( detx * detx + dety * dety ) ;
                                data_BGD     [ix][iy][iz] = ( ( offset_peak - offset_zero ) * exp( -1.0 * ( offset * offset ) / ( offset_sigma * offset_sigma ) ) ) + offset_zero ;
                                data_BGD_RECO[ix][iy][iz] = ( ( offset_peak - offset_zero ) * exp( -1.0 * ( offset * offset ) / ( offset_sigma * offset_sigma ) ) ) + offset_zero ;
                            }
                            else
                            {
                                printf( "Error, unrecognized tableflavor '%d', exiting...", tableflavor ) ;
                                exit( 1 );
                            }
                        } // endfor: loop over energy bins
                    } // endfor: loop over camera y bins
                } // endfor: loop over camera x bins
                // write data to columns
                if( fits_write_col( chunkptr, TFLOAT, C_DETX_LO , 1, 1, nelem_x  , &data_DETX_LO , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_DETX_HI , 1, 1, nelem_x  , &data_DETX_HI , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_DETY_LO , 1, 1, nelem_y  , &data_DETY_LO , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_DETY_HI , 1, 1, nelem_y  , &data_DETY_HI , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_LO, 1, 1, nelem_z  , &data_ENERG_LO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_ENERG_HI, 1, 1, nelem_z  , &data_ENERG_HI, &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_BGD     , 1, 1, nelem_xyz, &data_BGD     , &status ) )
                {
                    print_fits_error( status ) ;
                }
                if( fits_write_col( chunkptr, TFLOAT, C_BGD_RECO, 1, 1, nelem_xyz, &data_BGD_RECO, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
                // close fits file
                if( fits_close_file( chunkptr, &status ) )
                {
                    print_fits_error( status ) ;
                }
                
            }
            
            // break out of while loop
            break ;
            
        } // endwhile: scope for writing the background table to the chunk's fits file
        
        // check output file
        ifstream output_file( fitsFileName.c_str() ) ;
        if( ! output_file.good() )
        {
            cout << "Warning, output file '" << fitsFileName ;
            cout << "' may not have been closed successfully." << endl;
        }
        
    }
    
    //printf("\nArchive Summary:\n");
    //archiveEffArea->printAllArchiveFilenames() ;
    
    // close debug root file if needed
    if( debugRootFile )
    {
        debugRootFile->Close();
        fclose( debugAsciiFile ) ;
    }
    
    // cleanup
    anasumfile->Close();
    if( chunkList )
    {
        fclose( chunkList ) ;
    }
    cout << "writeCTAEventListFromAnasum has finished converting." << endl;
    return 0;
    
}

// met_zero : MJD time (in days) when MET is 0.0
void getRunTimes( VEvndispRunParameter* evrunpar, double met_zero, double& tstart, double& tstop, double& telapse, double& startMJD, double& stopMJD )
{
    bool debug = false ;
    if( debug )
    {
        cout << "getRunTimes()" << endl;
    }
    int sYe = 0 ;
    int sMo = 0 ;
    int sDa = 0 ;
    int sHo = 0 ;
    int sMi = 0 ;
    int sSe = 0 ;
    
    // run start
    if( debug )
    {
        cout << "fDBStart    :" << evrunpar->fDBRunStartTimeSQL.c_str() << endl;
    }
    sscanf( evrunpar->fDBRunStartTimeSQL.c_str(), "%4d-%2d-%2d %2d:%2d:%2d", &sYe, &sMo, &sDa, &sHo, &sMi, &sSe ) ;
    if( debug )
    {
        printf( "start string: %4d %2d %2d - %2d %2d %2d\n", sYe, sMo, sDa, sHo, sMi, sSe ) ;
    }
    int sMJDint = VSkyCoordinatesUtilities::getMJD( sYe, sMo, sDa ) ;
    startMJD = ( ( double )sMJDint ) + ( sHo / 24.0 ) + ( sMi / ( 24.*60. ) ) + ( sSe / ( 24.*60.*60. ) ) ;
    if( debug )
    {
        printf( "startMJD    : %10.4f\n", startMJD ) ;
    }
    
    // run end
    if( debug )
    {
        cout << "fDBStop     :" << evrunpar->fDBRunStoppTimeSQL.c_str() << endl;
    }
    sscanf( evrunpar->fDBRunStoppTimeSQL.c_str(), "%4d-%2d-%2d %2d:%2d:%2d", &sYe, &sMo, &sDa, &sHo, &sMi, &sSe ) ;
    if( debug )
    {
        printf( "stop string : %4d %2d %2d - %2d %2d %2d\n", sYe, sMo, sDa, sHo, sMi, sSe ) ;
    }
    sMJDint = VSkyCoordinatesUtilities::getMJD( sYe, sMo, sDa ) ;
    stopMJD = ( ( double )sMJDint ) + ( sHo / 24.0 ) + ( sMi / ( 24.*60. ) ) + ( sSe / ( 24.*60.*60. ) ) ;
    if( debug )
    {
        printf( "stopMJD     : %10.4f\n", stopMJD ) ;
    }
    
    // final times
    double met_corrected_start = ( startMJD - met_zero ) * 24 * 60 * 60 ;
    double met_corrected_stop  = ( stopMJD  - met_zero ) * 24 * 60 * 60 ;
    double met_duration        = met_corrected_stop - met_corrected_start ;
    if( debug )
    {
        cout << endl;
    }
    if( debug )
    {
        printf( "  met_corrected_start:%14.2f (s)\n", met_corrected_start ) ;
    }
    if( debug )
    {
        printf( "  met_corrected_stop :%14.2f (s)\n", met_corrected_stop ) ;
    }
    if( debug )
    {
        printf( "  met_duration       :%14.2f (s)\n", met_duration ) ;
    }
    tstart  = met_corrected_start ;
    tstop   = met_corrected_stop  ;
    telapse = met_duration       ;
    
}


string exec_shell( char* cmd )
{
    FILE* pipe = popen( cmd, "r" );
    if( !pipe )
    {
        printf( "Error: exec_shell(\"%s\"), unable to open pipe\n", cmd );
        return "ERROR" ;
    }
    char buffer[128];
    string result = "";
    while( !feof( pipe ) )
    {
        if( fgets( buffer, 128, pipe ) != NULL )
        {
            result += buffer;
        }
    }
    pclose( pipe );
    return result;
}

// get environment variable 'key'
// http://stackoverflow.com/questions/631664/accessing-environment-variables-in-c
std::string get_env_var( std::string const& key )
{
    char* val;
    val = getenv( key.c_str() );
    std::string retval = "";
    if( val != NULL )
    {
        retval = val;
    }
    return retval;
}

// quickly check if file exists
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool file_exists( const std::string& name )
{
    struct stat buffer;
    return ( stat( name.c_str(), &buffer ) == 0 );
}

void getHumanDateTimes( VEvndispRunParameter* evrunpar, char dateobs[], char timeobs[], char dateend[], char timeend[] )
{
    bool debug = false ;
    if( debug )
    {
        cout << "getHumanDateTimes()" << endl;
    }
    if( debug )
    {
        cout << "   startSQL:" << evrunpar->fDBRunStartTimeSQL << endl ;
    }
    if( debug )
    {
        cout << "   stopSQL :" << evrunpar->fDBRunStoppTimeSQL << endl ;
    }
    
    //cout << "fDBRunStartTimeSQL:'" << evrunpar->fDBRunStartTimeSQL << "'" << endl;
    sscanf( evrunpar->fDBRunStartTimeSQL.c_str(), "%s %s", dateobs, timeobs ) ;
    sscanf( evrunpar->fDBRunStoppTimeSQL.c_str(), "%s %s", dateend, timeend ) ;
    if( debug )
    {
        cout << "   dateobs:|" << dateobs << "|" << endl ;
    }
    if( debug )
    {
        cout << "   timeobs:|" << timeobs << "|" << endl ;
    }
    if( debug )
    {
        cout << "   dateend:|" << dateend << "|" << endl ;
    }
    if( debug )
    {
        cout << "   timeend:|" << timeend << "|" << endl ;
    }
    
}

// convert '##U' to time in seconds, where
//   ## is an integer, and
//   U  is either 's', 'm', or 'h',
//      for seconds, minutes, and hours, respectivly
double parseTimeArg( char timechar[] )
{
    string timestr = timechar ;
    //printf(" parseTimeArg(%s)\n", timestr.c_str() ) ;
    double multi = 0   ;
    double time  = 0.0 ;
    std::size_t found ;
    
    found = timestr.find( "s" ) ;
    if( found != std::string::npos )
    {
        multi = 1.0 ;
    }
    
    found = timestr.find( "m" ) ;
    if( found != std::string::npos )
    {
        multi = 60.0 ;
    }
    
    found = timestr.find( "h" ) ;
    if( found != std::string::npos )
    {
        multi = 3600.0 ;
    }
    
    if( multi == 0.0 )
    {
        printf( "error, parseTimeArg('%s'): need to specify either 's', 'm', or 'h' as the time unit for the -c option, exiting...", timestr.c_str() );
        exit( 1 );
    }
    
    time = atoi( timestr.c_str() ) ;
    time *= multi ;
    //printf(" parseTimeArg: multi=%f, time=%f\n", multi,time ) ;
    return time ;
}

// find the simple mean ra/dec/alt/az for a given span of time
void calcAvgPointing( const pointinginput inp, pointingoutput& out )
{
    // setup branches for reading from pointing data TChain
    UInt_t readMJD  = 0   ;
    double readTime = 0.0 ;
    double readRA   = 0.0 ;
    double readDE   = 0.0 ;
    double readAlt  = 0.0 ;
    double readAz   = 0.0 ;
    inp.ch->SetBranchAddress( "MJD"         , &readMJD ) ;
    inp.ch->SetBranchAddress( "Time"        , &readTime ) ;
    inp.ch->SetBranchAddress( "TelAzimuth"  , &readAz ) ;
    inp.ch->SetBranchAddress( "TelElevation", &readAlt ) ;
    inp.ch->SetBranchAddress( "TelRAJ2000"  , &readRA ) ;
    inp.ch->SetBranchAddress( "TelDecJ2000" , &readDE ) ;
    
    // setup loop variables
    double row_time     = 0.0 ;
    unsigned int npts   = inp.ch->GetEntries() ;
    unsigned int totpts = 0 ;
    double ramu  = 0.0 ;
    double decmu = 0.0 ;
    double altmu = 0.0 ;
    double azmu  = 0.0 ;
    
    int iearly = -1 ;
    
    // loop over rows of pointingDataReduced TChain
    for( unsigned int i = 0 ; i < npts ; i++ )
    {
        // load row
        inp.ch->GetEntry( i ) ;
        row_time = readMJD + ( readTime / 86400.0 ) ;
        
        // save last index before our time segment, in case this chunk's span is too small
        if( row_time < inp.beg )
        {
            iearly = i ;
        }
        
        // check if we're in the time chunk
        if( row_time > inp.beg && row_time < inp.end )
        {
            //printf( "i:%d sec:%.1f  rtime:%.6f  chbeg:%.6f  chend:%.6f !!\n", i, readTime, row_time, inp.beg, inp.end ) ;
            ramu   += readRA  ;
            decmu  += readDE  ;
            altmu  += readAlt ;
            azmu   += readAz  ;
            totpts += 1       ;
        }
        
    }
    
    // check if we got points to average over,
    // if not, try just using the last point before the span
    if( totpts <= 0 )
    {
        //printf("Warning, no pointing data found in tree 'pointingDataReduced' in input anasum file, for time chunk between %f and %f (spanning %.1f seconds). Attempting to use last good pointing datum (iearly=%d) before this span...", inp.beg, inp.end, (inp.end-inp.beg)*24*60*60.0, iearly);
        
        // check that our last point before the span is still good
        if( iearly >= 0 && ( unsigned int )iearly < npts )
        {
            inp.ch->GetEntry( iearly ) ;
            row_time = readMJD + ( readTime / 86400.0 ) ;
            printf( "calcAvgPointing(): warning, using one pointing datum (i=%d) for entire chunk (between %f and %f, spanning %.1f seconds) : RA=%.3f Dec=%.3f Alt=%.4f Az=%.4f ...\n", iearly, inp.beg, inp.end, ( inp.end - inp.beg ) * 24 * 60 * 60.0, readRA, readDE, readAlt, readAz );
            out.ra   = readRA  ;
            out.dec  = readDE  ;
            out.alt  = readAlt ;
            out.az   = readAz  ;
        }
        else
        {
            cout << "calcAvgPointing() Error: pointing row iearly=%d was not within bounds (0-%d), can't use last good pointing datum, exiting..." << endl;
            exit( 1 );
        }
    }
    else
    {
        // save output
        out.ra  = ( ramu  / totpts ) * TMath::RadToDeg() ;
        out.dec = ( decmu / totpts ) * TMath::RadToDeg() ;
        out.alt = ( altmu / totpts ) ;
        out.az  = ( azmu  / totpts ) ;
    }
    
    // ResetBranchAddresses() of the chain we just used.
    // because if we dont, it will silently fail,
    // and will rear its ugly head later when you try to do new
    // SetBranchAddress.  God Damn You ROOT.
    inp.ch->ResetBranchAddresses() ;
    
}


// for a given time mask 'vtm', find beginning and end of all time-regions which pass cuts (GTI's).
// will search for at most 'maxSecond' seconds from the beginning of the run
// will check for the edge of a cut-passing time region every 'tInterval' seconds
// will return a vector of the beginning and ends of the GTIs
// the returned GTI times will be correct to within +-(0.5*tInterval)
void calcGTIsFromTimeMask( VTimeMask* vtm, const double tInterval, const double maxSeconds, const double startMJD, const double stopMJD, vector<double>& gti_beg, vector<double>& gti_end )
{

    // time we are currently checking
    double searchTime = 0.0 ;
    double IsearchTime = 0.0 ;
    
    // last iteration's
    bool searchStatus  = false ; // false= time does not pass cuts
    bool lastStatus    = false ; // true = time does pass cuts
    bool IsearchStatus = false ;
    bool IlastStatus   = false ;
    
    // maximum number of loops to try
    int maxIter = ( int )( maxSeconds / tInterval ) ;
    
    // flag to exit both loops
    bool breakFlag = false ;
    
    // Loop over entire run duration, searching for GTIs
    for( int i = 0 ; i < maxIter ; i++ )
    {
    
        // instant in time we're checking
        searchTime = startMJD + ( ( i * tInterval ) / ( 24.*60.*60. ) ) ;
        
        // if we're at the end of the run, then break out of the loop
        if( searchTime > stopMJD )
        {
            break ;
        }
        
        // see if this point in time passes the cuts
        searchStatus = vtm->checkAgainstMask( searchTime ) ;
        
        // if we go from a does-not-pass-cuts instant to a does-pass-cuts instant
        // then its the start of a GTI
        // and we should start searching forwards for the end, aka
        // when we go from a does-pass-cuts instant to a does-not-pass-cuts instant
        if( searchStatus && !lastStatus )
        {
        
            // reset interior loop variables
            IsearchTime   = 0 ;
            IsearchStatus = false ;
            IlastStatus   = false ;
            for( int j = 0 ; j < maxIter ; j++ )
            {
                IsearchTime = searchTime + ( ( j * tInterval ) / ( 24.*60.*60. ) ) ;
                
                // if we're at the end of the run, then break out of the loop
                if( IsearchTime > stopMJD )
                {
                    breakFlag = true ;
                    break ;
                }
                
                // see if this point in time passes the cut mask
                IsearchStatus = vtm->checkAgainstMask( IsearchTime ) ;
                
                // if we go from a does-pass-cuts instant to a does-not-pass-cuts instant
                // OR, if we get to the end of the run
                // then we have reached the end of this GTI
                if( ( !IsearchStatus && IlastStatus ) || IsearchTime > stopMJD )
                {
                
                    // congrats, you've found a complete GTI
                    gti_beg.push_back( searchTime ) ;  // -((0.5*tInterval)/(24.*60.*60.)) ) ;
                    gti_end.push_back( IsearchTime ) ; // -((0.5*tInterval)/(24.*60.*60.)) ) ;
                    
                    // continue searching for a next GTI where the old one left off
                    i = i + j ;
                    
                    // if we're at the end of the run, jump out of both loops
                    // break out of this loop
                    break ;
                }
                IlastStatus = IsearchStatus ;
            }
            
            // if the inner loop set 'breakFlag', we should exit the outer loop too
            if( breakFlag )
            {
                break ;
            }
        }
        
        lastStatus = searchStatus ;
    }
    
    if( gti_beg.size() != gti_end.size() )
    {
        cout << "  Warning, gti_beg and gti_end aren't the same size." << endl;
        cout << "  Problem finding/saving GTI edge?" << endl;
    }
    if( gti_beg.size() == 0 || gti_end.size() == 0 )
    {
        //cout << "Warning, gti_beg/gti_end is size 0." << endl;
        cout << "    No GTI's found, assuming entire run is usable..." << endl;
    }
}

// sets telmask to all off
void telmask_clear( bool telmask[] )
{
    //cout << "   telmask_clear(): Clearing telmask to 0000" << endl;
    for( unsigned int i = 0 ; i < NTEL ; i++ )
    {
        //printf( "  testing index %3d : %d\n", i, telmask[i] ) ;
        telmask[i] = false ;
        //printf( "                now : %d\n", telmask[i] ) ;
    }
}

// tel: veritas telescope number, T1234 = 1, 2, 3, or 4
void telmask_set( bool telmask[], int tel, bool state )
{
    //cout << "telmask_set(): Setting T" << tel << " to state " << state << "..." << endl;
    if( state )
    {
        telmask[tel - 1] = true  ;
    }
    else
    {
        telmask[tel - 1] = false ;
    }
}

// tel: veritas telescope number, T1234 = 1, 2, 3, or 4
bool readImgSel( int ImgSel, int tel )
{
    bitset<8 * sizeof( ULong64_t )> imgsel = ImgSel;
    bool out = imgsel.test( tel - 1 ) ;
    //printf( "   readImgSel( %2d , %d ): %d\n", ImgSel, tel, out ) ;
    return out ;
}


// compare time span t1 and t2, and reduce them to to their smallest
// overlapping timespan (essentially AND the two timespans)
//
void reduceTimeSpans( const double t1_beg,  const double t1_end,
                      const double t2_beg,  const double t2_end,
                      double&       red_beg, double&       red_end )
{

    // find latest beginning time
    if( t1_beg > t2_beg )
    {
        red_beg = t1_beg ;
    }
    else if( t2_beg > t1_beg )
    {
        red_beg = t2_beg ;
    }
    else
    {
        red_beg = t1_beg ;
    }
    
    // find earliest ending time
    if( t1_end < t2_end )
    {
        red_end = t1_end ;
    }
    else if( t2_end < t1_end )
    {
        red_end = t2_end ;
    }
    else
    {
        red_end = t2_end ;
    }
    
}

void reduceTimeSpans( const Ntimespan t1, const Ntimespan t2, Ntimespan& reduced )
{

    // find latest beginning time
    if( t1.beg > t2.beg )
    {
        reduced.beg = t1.beg ;
    }
    else if( t2.beg > t1.beg )
    {
        reduced.beg = t2.beg ;
    }
    else
    {
        reduced.beg = t1.beg ;
    }
    
    // find earliest ending time
    if( t1.end < t2.end )
    {
        reduced.end = t1.end ;
    }
    else if( t2.end < t1.end )
    {
        reduced.end = t2.end ;
    }
    else
    {
        reduced.end = t2.end ;
    }
    
    printf( "\nreduceTimeSpans():\n" ) ;
    printf( "   t1:  beg:%6.1f end:%6.1f\n"  , t1.beg     , t1.end );
    printf( "   t2:  beg:%6.1f end:%6.1f\n"  , t2.beg     , t2.end );
    printf( "   red: beg:%6.1f end:%6.1f\n\n", reduced.beg, reduced.end );
    
    // sanity check
    if( reduced.beg > reduced.end )
    {
        printf( "Warning, reduceTimeSpans( Ntimespan, Ntimespan, Ntimespan ): beg(%f) > end(%f)\n", reduced.beg, reduced.end );
    }
}

void printChunkList( const double tmin, const double tmax, vector<Ntimespan>& chunklist )
{
    /*
    printf("\nchunklist:\n");
    for ( unsigned int i_chunk = 0 ; i_chunk < chunklist.size() ; i_chunk++ )
    {
      printf( " ch%2d: %6.1f - %6.1f (%6.1fs long)\n", i_chunk+1,
          chunklist[i_chunk].beg-tmin,
          chunklist[i_chunk].end-tmin,
          chunklist[i_chunk].end-chunklist[i_chunk].beg ) ;
    }
    */
    
    // setup printing variables
    printf( "printChunkList()\n" );
    double pc_res  = 60.0 ; // print run in 'pc_res'-second chunks
    int npc = ceil( ( tmax - tmin ) / pc_res ) ;
    double pc_beg, pc_end ;
    
    // character names to give each chunk
    char chunkNames[200] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789" ;
    
    // loop over each 60-second line
    for( int i_pc = 0 ; i_pc < npc ; i_pc++ )
    {
        // beginning and ending times of the line
        pc_beg = tmin + ( i_pc * pc_res ) ;
        pc_end = pc_beg + pc_res ;
        if( pc_end > tmax )
        {
            pc_end = tmax ;
        }
        
        // loop over each 1-second digit of the line
        double digtime = 0.0 ;
        string digits( 60, '~' ) ;
        int chunkiter ;
        for( unsigned int i_dig = 0 ; i_dig < 60 ; i_dig++ )
        {
            // find the time of the particular digit,
            // and figure out which chunk it belongs to
            digtime = pc_beg + i_dig ;
            chunkiter = checkChunkList( chunklist, digtime ) ;
            
            // assign it a chunk if it belongs to one,
            // or '~' if it is an unused part of the run,
            // or ' ' if its past the end of the run
            if( chunkiter >= 0 )
            {
                // loop back to the beginning of 'chunkNames' if
                // we go past the end ( 26 letters + 9 numbers = 35 )
                digits[i_dig] = chunkNames[chunkiter % 35] ;
            }
            else
            {
                if( digtime >= pc_end - 0.001 )
                {
                    digits[i_dig] = ' ' ;
                }
                else
                {
                    digits[i_dig] = '~' ;
                }
            }
        }
        
        // print our line
        printf( "min %2d - sec %6.1f - %-60s\n", i_pc, pc_beg - tmin, digits.c_str() ) ;
        
    }
    
}

// check to see which chunk time 'time' belongs to
// returns 0 to (nchunks-1) if it found the containing chunk
// returns -1 if 'time' does not belong to any chunk
// time: MET time
int checkChunkList( vector<Ntimespan>& chunklist, double time )
{
    for( unsigned int i = 0 ; i < chunklist.size() ; i++ )
    {
        if( chunklist[i].beg <= time  && time < chunklist[i].end )
        {
            return i ;
        }
    }
    return -1 ;
}

// print out cfitsio error messages and exit program
void print_fits_error( int status )
{
    if( status )
    {
        fits_report_error( stderr, status );
        exit( status ) ;
    }
    return ;
}

// calculate background values
void calcBACKGROUNDvalues( const inputbackground inp, outputbackground& out )
{
    // TODO: calculate proper background!
    // TODO: incorporate dead pixel information!
    out.bgd      = 1.2345678E-99 ; // background should not have zeros if no counts,
    out.bgd_reco = 2.3456789E-99 ; // should instead have really tiny numbers
}

// check that the 'iarg'th bound isnt >= 'narg',
// else throw error about unpaired 'flag' flag
// otherwise we crash if the user gives -n or -c as the last argument
void checkArgBounds( char flag[], int narg, int iarg )
{
    if( iarg + 1 >= narg )
    {
        printf( "Error, %d'th flag '%s' needs an argument, exiting...\n", iarg, flag )  ;
        exit( 1 ) ;
    }
}

// print the deadtimemask that is used to calculate the deadtime for a specific chunk of a run
// ichunk: the 'ith' chunk, only displayed
// beg   : beginning time of chunk, only displayed
// end   : ending time of chunk, only displayed
// timemask: vector bools, true for each second (since beginning of the run) where mask is valid,
//           false for each second is invalid
void printDeadTimeMask( const int ichunk, const double beg, const double end, const vector<bool> timemask )
{
    for( unsigned int i_sec = 0 ; i_sec < timemask.size() ; i_sec++ )
    {
        if( i_sec % 60 == 0 )
        {
            printf( "\nch:%2d (%.1f-%.1f) min:%2d - sec:%4d - ", ichunk, beg, end, i_sec / 60, i_sec ) ;
        }
        if( timemask[i_sec] )
        {
            printf( "T" );
        }
        else
        {
            printf( "-" );
        }
    }
    printf( "\n" );
    
}

string getEffAreaFileName( char* veritas_data_dir, unsigned int run_edver, char* run_epoch, int run_atmo, char* run_cut, float run_wobble, char* sim_type, char* eff_particle, int eff_zenith, int eff_noise )
{

    char efffilechar[500] = "" ;
    snprintf( efffilechar, sizeof efffilechar,
              "%s/IRFPRODUCTION/v%d/%s/%s_ATM%d_%s/EffectiveAreas_%s/EffArea-%s-%s-ID0-Ze%ddeg-%.1fwob-%d-%s.root",
              veritas_data_dir, run_edver, sim_type, run_epoch, run_atmo, eff_particle, run_cut, sim_type, run_epoch, eff_zenith,
              run_wobble, eff_noise, run_cut ) ;
              
    string efffile = efffilechar ;
    printf( "  efffile: '%s'\n", efffile.c_str() ) ;
    if( ! file_exists( efffile ) )
    {
        printf( "\nWarning, can't find requested file '%s'...\n\n", efffile.c_str() ) ;
    }
    else
    {
        printf( "found effective area file '%s'\n", efffile.c_str() ) ;
    }
    return efffile ;
    
}

double getPedVarFromRunSummary( TTree* runsum, int runid )
{
    int    tree_runid  = 0   ;
    double tree_pedvar = -100.0 ;
    runsum->SetBranchAddress( "runOn"    , &tree_runid ) ;
    runsum->SetBranchAddress( "pedvarsOn", &tree_pedvar ) ;
    double pedvar = 0.0 ;
    for( int i = 0 ; i < runsum->GetEntries() ; i++ )
    {
        runsum->GetEntry( i ) ;
        if( tree_runid == runid )
        {
            pedvar = tree_pedvar ;
            break ;
        }
    }
    if( pedvar < 0.1 || pedvar > 20.0 )
    {
        cout << "Warning: getPedVarFromRunSummary(): pedvar=" << pedvar << " out of bounds 0.1 < pedvar < 20.0, may not be trustable..." << endl;
    }
    runsum->ResetBranchAddresses() ;
    return pedvar ;
}

double getWobbleRadiusFromRunSummary( TTree* runsum, int runid )
{
    int    tree_runid   =    0   ;
    double wobrad       =    0.0 ;
    double tree_wobnort = -100.0 ;
    double tree_wobwest = -100.0 ;
    runsum->SetBranchAddress( "runOn"      , &tree_runid ) ;
    runsum->SetBranchAddress( "WobbleNorth", &tree_wobnort ) ;
    runsum->SetBranchAddress( "WobbleWest" , &tree_wobwest ) ;
    for( int i = 0 ; i < runsum->GetEntries() ; i++ )
    {
        runsum->GetEntry( i ) ;
        if( tree_runid == runid )
        {
            wobrad = sqrt( ( tree_wobnort * tree_wobnort ) + ( tree_wobwest * tree_wobwest ) ) ;
            break ;
        }
    }
    if( wobrad < 0.0 || wobrad > 4.0 )
    {
        cout << "Warning: getWobbleRadiusFromRunSummary(): wobrad=" << wobrad << " out of bounds 0.0 < wobrad < 4.0, may not be trustable..." << endl;
    }
    runsum->ResetBranchAddresses() ;
    return wobrad ;
}

VMagnetAzimuth::VMagnetAzimuth( string efile )
{
    fEffFile = efile ;
    
    // load up the fEffArea TTree
    TFile* tfile = new TFile( fEffFile.c_str(), "READ" ) ;
    TTree* ttree = ( TTree* ) tfile->Get( "fEffArea" ) ;
    ttree->SetEstimate( -1 ) ; // get all entries
    ttree->Draw( "az:azMin:azMax", "", "goff" ) ;
    int sentries = ttree->GetSelectedRows() ;
    
    // find all unique azimuth bin numbers
    double* values  = ttree->GetV1() ;
    double* valmins = ttree->GetV2() ;
    double* valmaxs = ttree->GetV3() ;
    vector<int>    uniqs_bins ; // azbin
    vector<double> uniqs_mins ; // azMin
    vector<double> uniqs_maxs ; // azMax
    int val ;
    double valmin, valmax ;
    for( int i = 0 ; i < sentries; i++ )
    {
        val    = ( int ) values[ i] ;
        valmin =       valmins[i] ;
        valmax =       valmaxs[i] ;
        
        // most azbins (0-7) are 90deg windows, but azbin 8 is the entire 360 deg range
        if( val == 8 )
        {
            continue ;    // skip azbin 8, since its the entire range
        }
        if( uniqs_bins.size() == 0 )
        {
            uniqs_bins.push_back( val ) ;
            uniqs_mins.push_back( valmin ) ;
            uniqs_maxs.push_back( valmax ) ;
        }
        else
        {
            // if 'val' is not in 'uniqs', add it to 'uniqs'
            if( std::find( uniqs_bins.begin(), uniqs_bins.end(), val ) == uniqs_bins.end() )
            {
                uniqs_bins.push_back( val ) ;
                uniqs_mins.push_back( valmin ) ;
                uniqs_maxs.push_back( valmax ) ;
            }
        }
    }
    
    for( unsigned int i = 0 ; i < uniqs_bins.size() ; i++ )
    {
        //printf("VMagnetAzimuth: azbin: %2d: %6.1f - %6.1f\n", uniqs_bins[i], uniqs_mins[i], uniqs_maxs[i] ) ;
        
        // calculate angle halfway between begaz and endaz (average the X and Y components, then arctan)
        double avgx = ( cos( uniqs_mins[i] * toRad ) + cos( uniqs_maxs[i] * toRad ) ) / 2.0 ;
        double avgy = ( sin( uniqs_mins[i] * toRad ) + sin( uniqs_maxs[i] * toRad ) ) / 2.0 ;
        double cent = atan2( avgy, avgx ) * toDeg ;
        
        //printf("VMagnetAzimuth: azbin: %2d: %7.1f to %7.1f (cent=%7.1f)\n", uniqs_bins[i], uniqs_mins[i], uniqs_maxs[i], cent) ;
        azBins.push_back( uniqs_bins[i] ) ;
        azMins.push_back( uniqs_mins[i] ) ;
        azMaxs.push_back( uniqs_maxs[i] ) ;
        azCens.push_back( cent ) ;
        
        //printf("azbins[%2d]: %.1f %.1f %.1f\n", bin, begaz, endaz, cent ) ;
    }
    
}

string VMagnetAzimuth::print()
{
    string result = "(" ;
    double cent = -1.0 ;
    char   res[20] = "" ;
    for( unsigned int i = 0 ; i < azCens.size() ; i++ )
    {
        if( i != 0 )
        {
            result += ",";
        }
        cent = azCens[i] ;
        if( cent < 0.0 )
        {
            cent = cent + 360.0 ;
        }
        sprintf( res, "%.1f", cent ) ;
        result += string( res ) ;
    }
    result += ")" ;
    return result ;
}

int VMagnetAzimuth::magnetize( double az )
{
    //double az = this->getAvgAzimuth() ;
    double angdist    =   0.0 ;
    double closestang = 180.0 ;
    int    closestbin =  -1   ;
    //vector<azbinspec> azbins = this->getAzBinRanges() ;
    for( unsigned int i = 0 ; i < azBins.size() ; i++ )
    {
        angdist = get_ang_dist( az, azCens[i] ) ;
        if( i == 0 )
        {
            closestang = angdist ;
            closestbin = i       ;
        }
        if( angdist < closestang )
        {
            closestang = angdist         ;
            closestbin = i               ;
        }
    }
    if( closestbin == -1 )
    {
        cout << "VMagnetAzimuth::magnitize(" << az << "): warning, closest bin not found..." << endl;
    }
    return closestbin ;
}

VMagnetZenith::VMagnetZenith( string efile )
{
    fEffFile = efile ;
    vector<double> zenithsd = get_unique_values_in_branch( fEffFile, "fEffArea", "ze" ) ;
    for( unsigned int i = 0 ; i < zenithsd.size() ; i++ )
    {
        zeniths.push_back( ( int ) zenithsd[i] ) ; // convert to int
    }
    
    sort( zeniths.begin(), zeniths.end() ) ;
    
}

string VMagnetZenith::print()
{
    string result  = "(" ;
    char   res[20] = ""  ;
    int    ang     = 99  ;
    for( unsigned int i = 0 ; i < zeniths.size() ; i++ )
    {
        if( i != 0 )
        {
            result += ",";
        }
        ang = zeniths[i] ;
        sprintf( res, "%d", ang ) ;
        result += string( res ) ;
    }
    result += ")" ;
    return result ;
}

int VMagnetZenith::magnetize( double zen )
{
    double angdist     =    0.0 ;
    double closestang  =  180.0 ;
    double closestdist = 1000.0 ;
    for( unsigned int i = 0 ; i < zeniths.size() ; i++ )
    {
        angdist = get_ang_dist( zen, ( double )( zeniths[i] ) ) ;
        if( angdist < closestdist )
        {
            closestdist = angdist ;
            closestang  = zeniths[i]  ;
        }
    }
    return closestang ;
}

VMagnetNoise::VMagnetNoise( string efile, int azbin, int zenith )
{

    fEffFile = efile  ;
    fAzbin   = azbin  ;
    fZenith  = zenith ;
    
    // load up the fEffArea TTree
    TFile* tfile = new TFile( fEffFile.c_str(), "READ" ) ;
    TTree* ttree = ( TTree* ) tfile->Get( "fEffArea" ) ;
    
    unsigned int nentries = ttree->GetEntries() ;
    printf( "nentries: %d\n", nentries ) ;
    
    // build list of unique noises and pedvars
    TChain* chainEffArea = new TChain() ;
    char treeName[1000] ;
    sprintf( treeName, "%s/fEffArea", fEffFile.c_str() ) ;
    chainEffArea->Add( treeName ) ;
    int    readNoise  = 0   ;
    float  readPedvar = 0.0 ;
    int    readAzbin  = 0   ;
    float  readZe     = 0   ;
    ttree->SetBranchAddress( "noise" , &readNoise );
    ttree->SetBranchAddress( "pedvar", &readPedvar );
    ttree->SetBranchAddress( "az"    , &readAzbin );
    ttree->SetBranchAddress( "ze"    , &readZe );
    ttree->SetBranchStatus( "hResponseMatrixFine", 0 ) ;
    
    for( unsigned int i = 0 ; i < nentries ; i++ )
    {
        ttree->GetEntry( i ) ;
        
        // skip all entries that aren't in this azimuth/zenith range
        if( !( readAzbin == azbin && readZe == zenith ) )
        {
            continue ;
        }
        
        if( noises.size() == 0 )
        {
            printf( "noises and pedvars vectors empty, adding (%d,%f)\n", readNoise, readPedvar ) ;
            noises.push_back( readNoise ) ;
            pedvars.push_back( readPedvar ) ;
        }
        // search for element readNoise in vector noises
        else if( find( noises.begin(), noises.end(), readNoise ) == noises.end() )
        {
            // if we didnt' find this noise, add it and its pedvar
            printf( "didnt find noise %d, adding it (and pedvar %f)\n", readNoise, readPedvar ) ;
            noises.push_back( readNoise ) ;
            pedvars.push_back( readPedvar ) ;
        }
    }
    
    ttree->ResetBranchAddresses() ;
    
    /*
    ttree->SetEstimate( -1 ) ; // get all entries
    char cutstring[500] = "" ;
    sprintf( cutstring, "az==%d && ze==%d", fAzbin, fZenith ); // only get the noise/pedvar pairs for a given telescope config
    ttree->Draw( "noise:pedvar", cutstring, "goff") ;
    int sentries = ttree->GetSelectedRows() ;
    //printf("cuts with '%s' went from %d events to %d events...\n", cutstring, (int) ttree->GetEntries(), sentries ) ;
    
    // find all unique noise/pedvar combinations
    double * val_noises  = ttree->GetV1() ;
    double * val_pedvars = ttree->GetV2() ;
    vector<noiseAndPedvarStruct> uncleanNAP ;
    for ( int i=1 ; i<sentries ; i++ )
    {
    
      // if ( uncleanNAP.size() == 0 || ( val_noises[i] != val_noises[i-1] && val_pedvars[i] != val_pedvars[i-1] ) )
      if ( val_noises[i] != val_noises[i-1] && val_pedvars[i] != val_pedvars[i-1] )
      {
        noiseAndPedvarStruct tmp ;
        tmp.noise  = (int) val_noises[ i] ;
        tmp.pedvar =       val_pedvars[i] ;
        uncleanNAP.push_back( tmp ) ;
      }
    }
    
    // add first entry
    printf("uncleanNAP.size() : %d\n", uncleanNAP.size() ) ;
    printf( "uncleanNAP[0].noise: %d\n", uncleanNAP[0].noise );
    noises.push_back(  uncleanNAP[0].noise  ) ;
    pedvars.push_back( uncleanNAP[0].pedvar ) ;
    //printf("new n/p: %4d/%7.4f\n", uncleanNAP[0].noise, uncleanNAP[0].pedvar ) ;
    
    for ( unsigned int i=1 ; i<uncleanNAP.size() ; i++ )
    {
      if ( std::find( noises.begin(), noises.end(), uncleanNAP[i].noise ) == noises.end() )
      {
        noises.push_back(  uncleanNAP[i].noise  ) ;
        pedvars.push_back( uncleanNAP[i].pedvar ) ;
      }
    }
    */
    
    /*
    for ( unsigned int i=0 ; i<noises.size() ; i++ )
    {
      printf(" uniq p/n: %4d, %7.4f\n", noises[i], pedvars[i] ) ;
    }
    */
    
    
}

void VMagnetNoise::printAll()
{
    printf( "%5s %6s\n", "noise", "pedvar" );
    
    for( unsigned int i = 0 ; i < pedvars.size() ; i++ )
    {
        printf( "%5d %6.2f\n", noises[i], pedvars[i] ) ;
    }
}

int VMagnetNoise::magnetize( double pedvar )
{
    double pedvardist   = 100.0 ;
    double closestdist  = 100.0 ;
    int    closestindex = -1    ;
    
    for( unsigned int i = 0 ; i < pedvars.size() ; i++ )
    {
        pedvardist = abs( pedvars[i] - pedvar ) ;
        //printf( " pedvardist i:%d abs( %.1f - %.1f )=%.1f\n", i, pedvars[i], pedvar, pedvardist ) ;
        if( pedvardist < closestdist )
        {
            closestdist = pedvardist ;
            closestindex = i ;
            //printf( "  update: closestindex=%d , closestdist=%.1f\n", closestindex, closestdist ) ;
        }
    }
    
    //printf(" closestindex %.2f %d %d\n", pedvar, closestindex, noises[closestindex] ) ;
    return noises[closestindex] ;
}

string VMagnetNoise::printNoises()
{
    char   res[20] = ""  ;
    string result  = "(" ;
    for( unsigned int i = 0 ; i < noises.size() ; i++ )
    {
        if( i != 0 )
        {
            result += "," ;
        }
        sprintf( res, "%d", noises[i] ) ;
        result += string( res ) ;
    }
    result += ")" ;
    return result ;
}

VMagnetOffset::VMagnetOffset( string efile )
{
    //printf("\nVMagnetOffset()\n");
    fEffArea = efile ;
    
    // this list has ugly floats (like 1.24999999999 and 2.00000000005), we want clean(er) 0.25, 1.75, etc
    vector<double> woffdirty = get_unique_values_in_branch( fEffArea, "fEffArea", "Woff" )  ;
    
    double twoff    ;
    int    twoffint ;
    double twoffclean ;
    vector<double> woff ;
    for( unsigned int i = 0 ; i < woffdirty.size() ; i++ )
    {
    
        // clean the float by multiplying by 4, since the offset is always a multiple of 0.25
        twoff = woffdirty[i] * 4.0 ;  // Since offset is always a multiple of 0.25, scale it up to ints
        twoff += 0.1 ;                // Make sure all dirty woff's are *above* their actual value
        twoffint = ( int ) twoff ;    // Round them all down the same way
        twoffclean = twoffint / 4.0 ; // Return them back to their 0.25 scale
        
        //printf(" twoffint[%d]:%f -> %d -> %f\n", i, woffdirty[i], twoffint, twoffclean ) ;
        if( woff.size() == 0 )
        {
            //printf("  cleaner: first element in woff:%f\n", twoffclean ) ;
            woff.push_back( twoffclean ) ;
        }
        else
        {
            if( std::find( woff.begin(), woff.end(), twoffclean ) == woff.end() )
            {
                //printf("  cleaner: adding %f to woff\n", twoffclean );
                woff.push_back( twoffclean ) ;
            }
        }
    }
    
    // sort woff so they end up in the right order for later fits table writing
    sort( woff.begin(), woff.end() ) ;
    
    binwidth = 0.125 ;
    double bincent = 0.0 ;
    double binupp  = 0.0 ;
    double binlow  = 0.0 ;
    for( unsigned int i = 0 ; i < woff.size() ; i++ )
    {
        bincent = woff[i] ;
        binupp  = bincent + binwidth ;
        binlow  = bincent - binwidth ;
        if( binlow < 0.0 )
        {
            binlow = 0.0 ;
        }
        //printf("  woff[%d]: %5.3f -> %4.2f -> %5.3f\n", i, binlow, bincent, binupp ) ;
        offsets.push_back( bincent ) ;
        offsetupperedges.push_back( binupp ) ;
        offsetloweredges.push_back( binlow ) ;
    }
    
}

string VMagnetOffset::print()
{
    string result = "(" ;
    char   res[20] = "" ;
    for( unsigned int i = 0 ; i < offsets.size() ; i++ )
    {
        if( i != 0 )
        {
            result += "," ;
        }
        sprintf( res, "%.2f", offsets[i] ) ;
        result += string( res ) ;
    }
    result += ")" ;
    return result ;
}

void VMagnetOffset::printAll()
{
    printf( "%6s %6s %6s\n", "binlow", "offset", "binupp" ) ;
    for( unsigned int i = 0 ; i < offsets.size() ; i++ )
    {
        printf( "%-6.3f %-6.2f %-6.3f\n", offsetloweredges[i], offsets[i], offsetupperedges[i] ) ;
    }
}

string VMagnetOffset::gernotFormat( double offset )
{
    char offsetstr[100] = "" ;
    if( offset == 0.5 )
    {
        snprintf( offsetstr, sizeof offsetstr, "%.1f", offset ) ;
    }
    else
    {
        snprintf( offsetstr, sizeof offsetstr, "%.2f", offset ) ;
    }
    return string( offsetstr ) ;
}

VArchiveEffArea::VArchiveEffArea( string efile )
{
    fEffFile = efile ;
    if( ! file_exists( fEffFile ) )
    {
        printf( "Error: VArchiveEffArea::VArchiveEffArea('%s'): unable to locate file '%s', exiting...\n", fEffFile.c_str(), fEffFile.c_str() ) ;
    }
}

int VArchiveEffArea::checkIfLoaded( archiveStruct archline )
{
    for( unsigned int i_arch = 0 ; i_arch < fArchive.size() ; i_arch++ )
    {
        // check if we have an entry in the archive
        if( archline.zenith   == fArchive[i_arch].zenith   &&
                archline.azbin    == fArchive[i_arch].azbin    &&
                archline.noise    == fArchive[i_arch].noise    &&
                archline.specind  == fArchive[i_arch].specind  &&
                archline.offset   == fArchive[i_arch].offset )
        {
            if( ! fLoaded[i_arch] )
            {
                // requested data in archive but not loaded to memory
                return DATA_IN_ARCHIVE_BUT_NOT_LOADED ;
            }
            else
            {
                // requested data is in archive and loaded, here's the index
                return i_arch ;
            }
        }
    }
    
    // requested data not in archive
    return DATA_NOT_IN_ARCHIVE ;
}

bool VArchiveEffArea::checkIfInArchive( archiveStruct archline )
{
    for( unsigned int i_arch = 0 ; i_arch < fArchive.size() ; i_arch++ )
    {
        if( archline.zenith   == fArchive[i_arch].zenith   &&
                archline.azbin    == fArchive[i_arch].azbin    &&
                archline.noise    == fArchive[i_arch].noise    &&
                archline.specind  == fArchive[i_arch].specind  &&
                archline.offset   == fArchive[i_arch].offset )
        {
            return true ;
        }
    }
    return false ;
}

void VArchiveEffArea::addArchive( archiveStruct archline )
{
    if( ! this->checkIfInArchive( archline ) )
    {
        string archstr = getArchiveROOTString( archline ) ;
        fArchive.push_back( archline ) ;
        fLoaded. push_back( false ) ;
        irf.     push_back( new VInstrumentResponseFunctionReader() ) ;
        psfDataStruct tmppsf ;
        fPsf.    push_back( tmppsf ) ;
        fenmigmatrix.push_back( 0 ) ;
        fEffEntry.push_back( 0 ) ;
    }
}

void VArchiveEffArea::loadAllArchives()
{

    for( unsigned int i_arch = 0 ; i_arch < fArchive.size() ; i_arch++ )
    {
        if( ! fLoaded[i_arch] )
        {
        
            // load one archive
            printf( "    now loading %s\n", getArchiveString( fArchive[i_arch] ).c_str() ) ; //zen(%2d) azbin(%2d) noise(%3d) wobble(%.2f) offset(%.2f) fname(...)\n",
            
            // load the VInstrumentResponseFunctionReader object
            irf[  i_arch]->fillData(
                fEffFile.c_str()         ,
                fArchive[i_arch].zenith  ,
                fArchive[i_arch].offset  ,
                fArchive[i_arch].azbin   ,
                fArchive[i_arch].specind ,
                fArchive[i_arch].noise   ,
                "A_MC" ) ;
                
            // read the effective area file, looping through the parameter space
            printf( "loading psf data... '%s'\n", getArchiveString( fArchive[i_arch] ).c_str() ) ;
            TFile* tfile = new TFile( this->fEffFile.c_str(), "READ" ) ;
            TTree* ttree = ( TTree* ) tfile->Get( "fEffArea" ) ;
            CEffArea* ceffarea = new CEffArea( ttree ) ;
            
            // loop through all points in the parameter space, saving all points that match our conditions
            vector<int   > matchingEntries      ;
            vector<double> phaseSpaceSimilarity ;
            
            // don't load the heavy hResponseMatrixFine histogram branch, since we're just
            // scanning for similar parameter space points first
            ceffarea->fChain->SetBranchStatus( "hResponseMatrixFine", 0 ) ;
            
            
            // each increment of these amounts counts for 1 unit of similarity
            // roughly chosen to be the bin size in each dimension
            double similarity  = 0.0 ;
            printf( "scanning through all points in the parameter space...\n" );
            printf( "ceffarea->fChain->GetEntries() : %lld\n", ceffarea->fChain->GetEntries() ) ;
            for( int i = 0 ; i < ceffarea->fChain->GetEntries() ; i++ )
            {
                // get a point in the parameter space
                ceffarea->GetEntry( i ) ;
                //printf("i:%d az:%d\n", i, ceffarea->az ) ;
                
                // calculate this point in the phaseSpace's similarity to the one we want.
                // generally, every bin-width further from the target counts as a similarity of 1,
                // lower similarity = better = closer to target
                similarity = 0 ;
                similarity += similarity_azbin( fArchive[i_arch].azbin  , ceffarea->az ) ;
                similarity += similarity_wobble( fArchive[i_arch].offset , ceffarea->Woff ) ;
                similarity += similarity_specind( fArchive[i_arch].specind, ceffarea->index ) ;
                similarity += similarity_noise( fArchive[i_arch].noise  , ceffarea->noise ) ;
                similarity += similarity_zenith( fArchive[i_arch].zenith , ceffarea->ze ) ;
                
                // check if we've reached the right point in the parameter space
                //if ( ceffarea->az    !=     fArchive[i_arch].azbin            ) continue ;
                //if ( fabs(ceffarea->index - fArchive[i_arch].specind ) > 0.05 ) continue ;
                //if ( fabs(ceffarea->Woff  - fArchive[i_arch].offset  ) > 0.05 ) continue ;
                //if ( ceffarea->noise !=     fArchive[i_arch].noise            ) continue ;
                //if ( fabs(ceffarea->ze    - fArchive[i_arch].zenith  ) > 3.0  ) continue ;
                
                matchingEntries.push_back( i ) ;
                phaseSpaceSimilarity.push_back( similarity ) ;
                
            }
            printf( "scanned %d entries in the parameter space...\n", ( int ) matchingEntries.size() ) ;
            
            // among all entries, find the most similar one
            int    mostSimilarEntry =  -1   ;
            double lowestSimilarity = 100.0 ;
            for( unsigned int ime = 0 ; ime < matchingEntries.size() ; ime ++ )
            {
                //printf( "matchingEntries[%2d]:%3d    phaseSpaceSimilarity[%2d]:%f \n", ime, matchingEntries[ime], ime, phaseSpaceSimilarity[ime] ) ;
                if( phaseSpaceSimilarity[ime] < lowestSimilarity )
                {
                    lowestSimilarity = phaseSpaceSimilarity[ime] ;
                    mostSimilarEntry = matchingEntries[ime]      ;
                }
            }
            
            //printf("\n") ;
            //printf("most similar entry is %5d with a similarity of %f ...\n", mostSimilarEntry, lowestSimilarity ) ;
            
            // now that we've found the most similar parameter space point,
            // reneable the heavy energy migration matrix branch
            ceffarea->fChain->SetBranchStatus( "hResponseMatrixFine", 1 ) ;
            
            // ceffarea->GetEntry( mostSimilarEntry ) ;
            
            // at this point, we know there is only one matching point in the parameter space
            //printf("  matching entry %d\n", mostSimilarEntry ) ;
            fEffEntry[i_arch] = mostSimilarEntry ;
            printf( "fEffEntry[%d] = %d\n", i_arch, mostSimilarEntry ) ;
            ceffarea->GetEntry( mostSimilarEntry ) ;
            
            // loop over energy bins to load our psf containment radii
            for( int i_en = 0 ; i_en < ceffarea->Rec_nbins ; i_en++ )
            {
                printf( "loading irf %5d into archive %2d : ", mostSimilarEntry, i_arch ) ;
                printf( "i_en:%2d  "    , i_en ) ;
                printf( "e0:%5.2f  "    , ceffarea->Rec_e0[i_en] ) ;
                printf( "eff:%6.0f  "   , ceffarea->Rec_eff[i_en] ) ;
                //printf("rad68:%5.3f   ", ceffarea->Rec_angRes_p68[i_en] ) ;
                printf( "sigma:%.4f  ", ceffarea->Rec_angRes_kingSigma[i_en] ) ;
                printf( "gamma:%.4f  ", ceffarea->Rec_angRes_kingGamma[i_en] ) ;
                printf( "\n" );
                fPsf[i_arch].energy.push_back( ceffarea->Rec_e0[        i_en] ) ;
                fPsf[i_arch].rad68.push_back( ceffarea->Rec_angRes_p68[i_en] ) ;
                fPsf[i_arch].rad80.push_back( ceffarea->Rec_angRes_p80[i_en] ) ;
                fPsf[i_arch].nativeking_sigma.push_back( ceffarea->Rec_angRes_kingSigma[i_en] ) ;
                fPsf[i_arch].nativeking_gamma.push_back( ceffarea->Rec_angRes_kingGamma[i_en] ) ;
            }
            printf( "  read in %d energy bins\n", ( int ) fPsf[i_arch].energy.size() );
            printf( "\n" );
            
            int    el = 90 - fArchive[i_arch].zenith ;
            int    az = fArchive[i_arch].azbin ;
            int    no = fArchive[i_arch].noise ;
            //double si = fArchive[i_arch].specind ;
            double of = fArchive[i_arch].offset ;
            double en ;
            for( int i_en = 0 ; i_en < ceffarea->Rec_nbins ; i_en++ )
            {
                en = ceffarea->Rec_e0[i_en] ;
                printf( "psf energies entry=%5d el%2d az%d offset%4.2f noise%4d i_arch=%2d  i_en=%2d  %5.3f = %6.2fTeV\n", mostSimilarEntry, el, az, of, no, i_arch, i_en, en, pow( 10, en ) ) ;
            }
            
            /////////////////////////////
            // ENERGY MIGRATION MATRIX //
            /////////////////////////////
            printf( "cloning energy migration matrix...\n" ) ;
            char hname[200] ;
            sprintf( hname, "migmatrix_arch%d", i_arch );
            if( ceffarea->hResponseMatrixFine )
            {
                fenmigmatrix[i_arch] = ( TH2D* ) ceffarea->hResponseMatrixFine->Clone( hname ) ;
                sprintf( hname, "migration matrix for elev=%2d azbin=%2d specind=%4.2f offset=%4.2f noise=%4d", 90 - fArchive[i_arch].zenith, fArchive[i_arch].azbin, fArchive[i_arch].specind, fArchive[i_arch].offset, fArchive[i_arch].noise ) ;
                fenmigmatrix[i_arch]->SetTitle( hname ) ;
                //int ibin = 0 ;
                //ibin = fenmigmatrix[i_arch]->GetBin(  200, 205 ) ;
                //fenmigmatrix[i_arch]->SetBinContent( ibin, 300 ) ;
                //ibin = fenmigmatrix[i_arch]->GetBin(  200, 206 ) ;
                //fenmigmatrix[i_arch]->SetBinContent( ibin, 300 ) ;
                //ibin = fenmigmatrix[i_arch]->GetBin(  201, 205 ) ;
                //fenmigmatrix[i_arch]->SetBinContent( ibin, 300 ) ;
                //ibin = fenmigmatrix[i_arch]->GetBin(  201, 206 ) ;
                //fenmigmatrix[i_arch]->SetBinContent( ibin, 300 ) ;
            }
            else
            {
                printf( "error, mig matrix missing...\n" );
                printf( " -> effarea file : %s\n", fEffFile.c_str() ) ;
                printf( " -> entry : %d\n", mostSimilarEntry ) ;
                printf( " -> hResponseMatrixFine\n" ) ;
                exit( 1 );
            }
            
            // record that this archive was loaded to memory
            fLoaded[i_arch] = true ;
            
        }
    }
}

// make sure we only have one parameter point in 'matent',
// all other quantities of parameter points are errors
void VArchiveEffArea::onlyOneParameterPointAllowed( archiveStruct archline, vector<int> matent, CEffArea* c )
{

    if( matent.size() != 1 )
    {
        cout << "VArchiveEffArea::getPsfData( " << getArchiveString( archline ) << " ) : Error, was expecting 1 point in parameter space, found " << matent.size() << ", exiting..." << endl;
        if( matent.size() > 1 )
        {
            int sentries = matent.size() ;
            for( int i = 0 ; i < sentries ; i++ )
            {
                c->GetEntry( matent[i] ) ;
                printf( "    many matching effarea points:  i:%4d e:%6d  =  ze%2.0f  azbin%2d  wob%4.2f  noise%04d  specind%3.1f\n", i, matent[i], c->ze, c->az, c->Woff, c->noise, c->index ) ;
            }
        }
        exit( 1 );
    }
    
}

void VArchiveEffArea::printAllArchiveFilenames()
{
    for( unsigned int i_arch = 0 ; i_arch < fArchive.size() ; i_arch++ )
    {
        printf( "requires %2d '%s'\n", fArchive[i_arch].azbin, fEffFile.c_str() ) ;
    }
}

double get_pedvar_from_eff( char* fname )
{
    if( ! file_exists( fname ) )
    {
        printf( "Error: get_pedvar_from_eff(fname='%s'): file '%s' doesn't seem to exist, exiting...\n", fname, fname ) ;
        exit( 1 ) ;
    }
    double pedvar = -1.0 ;
    TFile* searchfile = new TFile( fname, "READ" ) ;
    TTree* searchtree = ( TTree* ) searchfile->Get( "fEffArea" ) ;
    searchtree->SetBranchAddress( "pedvar", &pedvar );
    searchtree->GetEntry( 2 ) ;
    //printf("get_pedvar_from_eff(): pedvar:%5.2f  file:'%s'\n", pedvar, fname ) ;
    searchtree->ResetBranchAddresses() ;
    return pedvar ;
}

// angular distance between two angles (deg,deg, return deg)
// arccos( dotproduct(ang1,ang2) )
double      get_ang_dist( double ang1, double ang2 )
{
    return acos( ( cos( ang1 * toRad ) * cos( ang2 * toRad ) ) + ( sin( ang1 * toRad ) * sin( ang2 * toRad ) ) ) * toDeg ;
}


int calcKingParametersRootFit( double rad68, double rad80, double& sigma, double& gamma )
{
    KingFitter* kf = new KingFitter( rad68, rad80 ) ;
    kf->fit() ;
    sigma = kf->getSigma() ;
    gamma = kf->getGamma() ;
    return 0 ;
}


KingFitter::KingFitter( double r68, double r80 )
{
    fRad68 = r68 ;
    fRad80 = r80 ;
    fSigma = 0.0 ;
    fGamma = 0.0 ;
}

// integrate the psf from 0->r,
// return the containment fraction, from 0.0 to 1.0
// par is a double[2] list of sigma, gamma
// NOTE: THIS IS NOT THE PSF, its the integral of the psf
double kingfunc( Float_t r, Double_t* par )
{
    // from mathematica snippet in KingPsf.nb
    return 1 - pow( 2, -1 + par[1] ) * pow( par[0], -2 + 2 * par[1] ) * pow( par[1] / ( pow( r, 2 ) + 2 * par[1] * pow( par[0], 2 ) ), -1 + par[1] ) ;
}

void kingfcn( Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag )
{
    //calculate chisquare
    double chisq = 0;
    double delta;
    for( int i = 0 ; i < KINGFITTER_NPAR ; i++ )
    {
        delta  = ( KingFitYArray[i] - kingfunc( KingFitXArray[i], par ) ) / KingFitYArrayErr[i] ;
        chisq += delta * delta;
    }
    f = chisq;
}

void KingFitter::fit()
{
    KingFitXArray[0] = fRad68 ;
    KingFitXArray[1] = fRad80 ;
    KingFitYArray[0] = 0.68   ;
    KingFitYArray[1] = 0.80   ;
    KingFitYArrayErr[0] = 0.0015 ; // TODO: 0.15% error? why choose this?
    KingFitYArrayErr[1] = 0.0015 ;
    
    // init TMinuit object
    TMinuit* gMinuit = new TMinuit( KINGFITTER_NPAR ) ;
    gMinuit->SetPrintLevel( -1 ) ;
    gMinuit->SetFCN( kingfcn ) ;
    
    Double_t arglist[10] ; // TODO: why is this 10?
    Int_t ierflg = 0 ;
    arglist[0] = 1;
    gMinuit->mnexcm( "SET ERR", arglist, 1, ierflg ) ;
    
    // set starting values
    static Double_t vstart[KINGFITTER_NPAR] = { 0.03 , 2.0 } ;
    static Double_t step[  KINGFITTER_NPAR] = { 0.001, 0.1 } ;
    gMinuit->mnparm( 0, "a1", vstart[0], step[0], 0, 0, ierflg ) ;
    gMinuit->mnparm( 1, "a2", vstart[1], step[1], 0, 0, ierflg ) ;
    
    // minimization step
    // in terminal, do
    arglist[0] = 2000 ; // maximum number of minimization loops
    arglist[1] = 0.1  ; // minimization tolerance
    gMinuit->mnexcm( "MIGRAD", arglist, 2, ierflg ) ;
    
    // Print results
    Double_t amin, edm, errdef ;
    Int_t nvpar, nparx, icstat ;
    gMinuit->mnstat( amin, edm, errdef, nvpar, nparx, icstat ) ;
    //gMinuit->mnprin( 3, amin ) ;
    
    
    double sigerr, gamerr ;
    gMinuit->GetParameter( 0, fSigma, sigerr );
    gMinuit->GetParameter( 1, fGamma, gamerr );
    //Double_t pars[2] = { fSigma, fGamma } ;
    //double contain68 = kingfunc( fRad68, pars ) ;
    //double contain80 = kingfunc( fRad80, pars ) ;
    //printf("rootfit: r68=%5.3f  r80=%5.3f  :  ", fRad68, fRad80 ) ;
    //printf("status=%2d amin=%13.6f edm=%21.8f errdef=%7f nvpar=%3d nparx=%3d icstat=%3d  :  ", ierflg, amin, edm, errdef, nvpar, nparx, icstat ) ;
    //printf("sigma=%6.4f sigmaerr=%9.4f gamma=%8.4f gammaerr=%11.4f  :  ", fSigma, sigerr, fGamma, gamerr );
    //printf("cont68=%6.4f cont80=%6.4f", contain68, contain80 ) ;
    //printf("\n");
}


// from 68 and 80 % containment radii (in degrees), calculate the sigma and gamma
// parameters for a normalized king function (WITH THE SECANT METHOD, fails if r80 < r68*1.1)
int calcKingParametersViaSecant( double rad68, double rad80, double& sigma, double& gamma )
{
    // copied almost verbatim from
    // $CTOOLS/scripts/cta_root2caldb.py's root2psf_king() function
    
    // no point in continuing if only 0's
    if( rad68 <= 0.0 || rad80 <= 0.0 )
    {
        printf( "                                           calcKingParameters(%.2f %.2f %.2f %.2f): Warning, a containment radius is 0, setting g/s to 0s...\n", rad68, rad80, sigma, gamma ) ;
        gamma = 0.0 ;
        sigma = 0.0 ;
        return 1 ;
    }
    
    // initial ratios
    double a = 1.0 - 0.68 ;
    double b = 1.0 - 0.80 ;
    double c = ( rad68 * rad68 ) / ( rad80 * rad80 ) ;
    
    // solve equation (a^x-1)/(b^x-1)=c for x using secant
    // method.  Stop when we are better than 1e-6;
    double x, f ;
    double x1 = -0.5 ;
    double x2 = -10000.0 ; // if too close to -0.5, will fail to fit some combinations of rad68 and rad80
    double f1 = ( ( pow( a, x1 ) - 1.0 ) / ( pow( b, x1 ) - 1.0 ) ) - c ;
    double f2 = ( ( pow( a, x2 ) - 1.0 ) / ( pow( b, x2 ) - 1.0 ) ) - c ;
    int iter  = 0 ;
    while( true )
    {
        x = x1 - f1 * ( x1 - x2 ) / ( f1 - f2 ) ;
        f = ( pow( a, x ) - 1.0 ) / ( pow( b, x ) - 1.0 ) - c ;
        iter += 1 ;
        if( iter == 1 )
        {
            printf( "\n" );
        }
        
        char xwarn[100] = "" ;
        if( x > 0.0 )
        {
            sprintf( xwarn, "%sx:%.8f%s", KRED, x, KNRM ) ;
        }
        else
        {
            sprintf( xwarn, "x:%.8f"    , x ) ;
        }
        printf( "              calcKingParameters(%7.3f %7.3f %7.3f %7.3f): iter:%d  diff:%.8f  %s  f1=%f  f2=%f...\n", rad68, rad80, sigma, gamma, iter, fabs( f ), xwarn, f1, f2 ) ;
        
        if( fabs( f ) < 1.0e-6 )
        {
            //printf("              calcKingParameters(%7.3f %7.3f %7.3f %7.3f): Warning, breaking!...\n", rad68, rad80, sigma, gamma ) ;
            break ;
        }
        else
        {
            f2 = f1 ;
            x2 = x1 ;
            f1 = f  ;
            x1 = x  ;
        }
        
        // compute gamma
        if( x < 0.0 )
        {
            gamma = 1.0 - ( 1.0 / x ) ;
        }
        else
        {
            gamma = 1.0 ;
            //printf("              calcKingParameters(%7.3f %7.3f %7.3f %7.3f): Warning, x(%7.4f) >= 0, setting gamma to 1...\n", rad68, rad80, sigma, gamma, x ) ;
        }
        
        // compute sigma
        double denom = 2.0 * gamma * ( pow( b, x ) - 1.0 ) ;
        if( denom > 0.0 )
        {
            sigma = rad68 * sqrt( 1.0 / denom ) ;
        }
        else
        {
            //printf("              calcKingParameters(%7.3f %7.3f %7.3f %7.3f): Warning, first denom is 0, attempting 2nd method...\n", rad68, rad80, sigma, gamma ) ;
            denom = 2.0 * gamma * ( pow( b, x ) - 1.0 ) ;
            if( denom > 0.0 )
            {
                sigma = rad80 * sqrt( 1.0 / denom ) ;
            }
            else
            {
                //printf("              calcKingParameters(%7.3f %7.3f %7.3f %7.3f): Warning, denom is still 0.0, setting g/s to 0s...\n", rad68, rad80, sigma, gamma ) ;
                gamma = 0.0 ;
                sigma = 0.0 ;
                return 0 ;
            }
        }
    }
    
    return 0 ;
}

// from 68 and 80 % containment radii (in degrees), calculate the sigma and gamma
// parameters for a normalized king function (WITH THE NEWTON METHOD)
int calcKingParameters( double rad68, double rad80, double& sigma, double& gamma )
{

    // copied almost verbatim from
    // $CTOOLS/scripts/cta_root2caldb.py's root2psf_king() function
    
    // no point in continuing if only 0's
    if( rad68 <= 0.0 || rad80 <= 0.0 )
    {
        gamma = 0.0 ;
        sigma = 0.0 ;
        return 1 ;
    }
    
    // initial ratios
    double a = 1.0 - 0.68 ;
    double b = 1.0 - 0.80 ;
    double c = ( rad68 * rad68 ) / ( rad80 * rad80 ) ;
    
    // f1 = f(x1)
    // g(x1) = f'(x1)
    double x, f1, g1 ;
    double x1 = 0.1 ;
    
    int iter  = 0 ;
    while( true )
    {
        f1 = ( pow( a, x1 ) - 1.0 ) / ( pow( b, x1 ) - 1.0 ) - c ;
        g1 = ( ( pow( a, x1 ) * log( a ) ) / ( pow( b, x1 ) - 1.0 ) ) - ( ( ( pow( a, x1 ) - 1 ) * pow( b, x1 ) * log( b ) ) / pow( pow( b, x1 ) - 1, 2 ) ) ;
        x = x1 - ( f1 / g1 ) ;
        iter += 1 ;
        
        if( iter == 0 )
        {
            printf( "\n" );
        }
        printf( "  calcKingParameters(%7.3f %7.3f %7.3f %7.3f): iter:%d  diff:%.8f  x:%.8f  x1:%f  f1:%f  g1:%f ...\n", rad68, rad80, sigma, gamma, iter, fabs( f1 ), x, x1, f1, g1 ) ;
        
        if( fabs( f1 ) < 1.0e-6 )
        {
            break ;
        }
        else
        {
            x1 = x ;
        }
        
        // compute gamma
        if( x < 0.0 )
        {
            gamma = 1.0 - ( 1.0 / x ) ;
        }
        else
        {
            gamma = 1.0 ;
        }
        
        // compute sigma
        double denom = 2.0 * gamma * ( pow( b, x ) - 1.0 ) ;
        if( denom > 0.0 )
        {
            sigma = rad68 * sqrt( 1.0 / denom ) ;
        }
        else
        {
            denom = 2.0 * gamma * ( pow( b, x ) - 1.0 ) ;
            if( denom > 0.0 )
            {
                sigma = rad80 * sqrt( 1.0 / denom ) ;
            }
            else
            {
                gamma = 0.0 ;
                sigma = 0.0 ;
                return 0 ;
            }
        }
    }
    
    return 0 ;
}


string command_line_regex( string input, string regcmd )
{
    //printf("command_line_regex( input='%s', regcmd='%s'\n", input.c_str(), regcmd.c_str() ) ;
    char cmd[1000] = "" ;
    
    // create a command-line command which uses a grep-regex to extract the needed part of the string
    sprintf( cmd, "echo \"%s\" | %s", input.c_str(), regcmd.c_str() ) ;
    //printf("cmd:'%s'\n", cmd ) ;
    
    string result = exec_shell( cmd ) ;
    result.erase( std::remove( result.begin(), result.end(), '\n' ), result.end() ) ;
    //printf("result:'%s'\n", result.c_str() ) ;
    if( result.size() == 0 )
    {
        printf( "Error get_command_line_regex( input='%s', regcmd='%s' ): no resulting string, exiting...\n", input.c_str(), regcmd.c_str() ) ;
        exit( 1 );
    }
    return result ;
}

vector<double> get_unique_values_in_branch( string rootfile, string treename, string branchname )
{
    //printf("\n");
    //printf("get_uniq_values_in_branch( rootfile='...', treename='%s', branchname='%s' )\n", treename.c_str(), branchname.c_str() ) ;
    
    TFile* tfile = new TFile( rootfile.c_str(), "READ" ) ;
    TTree* ttree = ( TTree* ) tfile->Get( treename.c_str() ) ;
    //int nentries  = ttree->GetEntries() ;
    ttree->SetEstimate( -1 ) ; // get all entries
    //printf("read fEffArea: TTree with %d entries...\n", nentries ) ;
    
    ttree->Draw( branchname.c_str(), "", "goff" ) ;
    int sentries = ttree->GetSelectedRows() ;
    //printf("selected %d entries...\n", nentries ) ;
    
    double* values = ttree->GetV1() ;
    vector<double> uniqs ;
    for( int i = 0 ; i < sentries; i++ )
    {
        if( std::find( uniqs.begin(), uniqs.end(), values[i] ) == uniqs.end() )
        {
            uniqs.push_back( values[i] );
            //printf("get_unique_values_in_branch(...,%s,%s): adding: %f\n", treename.c_str(), branchname.c_str(), values[i]);
        }
    }
    
    /*
    for ( unsigned int i=0 ; i<uniqs.size() ; i++ )
    {
      printf("get_unique_values_in_branch(...,tree='%s',branch='%s'):%f\n", treename.c_str(), branchname.c_str(), uniqs[i]) ;
    }
    */
    
    tfile->Close() ;
    return uniqs ;
}


int VPsfBlock::scanBlockForAlternateFit( int i_offset, int i_en )
{

    string str = formatPsfBlockStructString( fPsfBlock[i_offset][i_en] ) ;
    //printf("scan %s\n", str.c_str() ) ;
    
    // number of elements in each dimension
    int nOffsetRows = fPsfBlock.size() ;
    int nEnergyRows = fPsfBlock[0].size() ;
    
    // scan direction/radius, relative to i_offset and i_en
    int scanradius ;
    int scandirect ;
    
    // calculated param space to scan
    int targ_offset = 0 ;
    int targ_en     = 0 ;
    
    // scan through nearby nodes
    int maxradius = nOffsetRows ;
    if( nEnergyRows > maxradius )
    {
        maxradius = nEnergyRows ;
    }
    for( scanradius = 1 ; scanradius < maxradius ; scanradius++ )
    {
        // loop through each parameter space direction
        for( scandirect = 0 ; scandirect < 4 ; scandirect++ )
        {
            //printf("scan  rad=%d  dir=%d\n", scanradius, scandirect ) ;
            // calculate the target node to scan
            targ_offset = i_offset ;
            targ_en     = i_en     ;
            if( scandirect == 0 )
            {
                targ_offset += scanradius ;
            }
            else if( scandirect == 1 )
            {
                targ_en     += scanradius ;
            }
            else if( scandirect == 2 )
            {
                targ_offset -= scanradius ;
            }
            else
            {
                targ_en     -= scanradius ;
            }
            //printf("scan    initial: targ_offset=%d  targ_en=%d\n", targ_offset, targ_en ) ;
            
            // check that targ_offset and targ_en are within sane bounds,
            // skip them if they aren't
            if( targ_offset <  0 )
            {
                continue ;
            }
            if( targ_offset >= nOffsetRows )
            {
                continue ;
            }
            if( targ_en     <  0 )
            {
                continue ;
            }
            if( targ_en     >= nEnergyRows )
            {
                continue ;
            }
            //printf("scan    bounded: targ_offset=%d  targ_en=%d\n", targ_offset, targ_en ) ;
            //printf("scan    rad=%2d dir=%d : targ_off=%2d targ_en=%2d : targ goodfit:%s\n", scanradius, scandirect, targ_offset, targ_en, fPsfBlock[targ_offset][targ_en].goodfit ? "true" : "false" ) ;
            
            // scan the target node,
            // skip to the next target node if we're still at the same node
            if( targ_offset == i_offset && targ_en == i_en )
            {
                continue ;
            }
            // skip to the next target node if this one didn't have a converged fit
            if( ! fPsfBlock[targ_offset][targ_en].goodfit )
            {
                continue ;
            }
            
            // if it did have a good fit, then copy the sigma/gamma from the target node
            fPsfBlock[i_offset][i_en].sigma           = fPsfBlock[targ_offset][targ_en].sigma ;
            fPsfBlock[i_offset][i_en].gamma           = fPsfBlock[targ_offset][targ_en].gamma ;
            fPsfBlock[i_offset][i_en].wasinterpolated = true ;
            printf( "scan: alternate found: %s (will use king fit form fPsfBlock[%d][%d])...\n", str.c_str(), targ_offset, targ_en ) ;
            return 0 ;
            
        }
    }
    //printf(" interp! %s failed to find interpolation!\n", str.c_str() ) ;
    return 1 ;
    
}

// class for storing and checking the psf parameter space
VPsfBlock::VPsfBlock( int nOffRows, int nEnRows )
{
    // fill our vector
    fPsfBlock.resize( nOffRows ) ;
    for( int i = 0 ; i < nOffRows ; i++ )
    {
        fPsfBlock[i].resize( nEnRows ) ;
    }
}

// take all the information in a psfBlockPartStruct,
// and format it to a string for easy and unified printing
string formatPsfBlockStructString( psfBlockPartStruct psf )
{
    char chaar[1000] = "" ;
    sprintf( chaar, "psfblock i_off=%2d i_en=%2d | r68=%5.3f r80=%5.3f | gfit=%d wasint=%d | sigma=%6.4f gamma=%8.4f",
             psf.i_offset, psf.i_en, psf.rad68, psf.rad80, psf.goodfit, psf.wasinterpolated, psf.sigma, psf.gamma ) ;
    string str = chaar ;
    return str ;
}

// print a 2d array of the parameter space
// at each (offset,energy) node,
//   print 1 if the initial TMinuit King Function fit converged,
//   print 0 otherwise
void VPsfBlock::formatBlock_goodfit()
{
    printf( "formatBlock_goodfit     i_en:01234567890123456789012345678901\n" );
    for( unsigned int i_off = 0 ; i_off < fPsfBlock.size() ; i_off++ )
    {
        printf( "formatBlock_goodfit i_off=%2d ", i_off ) ;
        for( unsigned int i_en = 0 ; i_en < fPsfBlock[i_off].size() ; i_en++ )
        {
            printf( "%d", fPsfBlock[i_off][i_en].goodfit ) ;
        }
        printf( "\n" );
    }
    printf( "\n" );
}

// print a 2d array of the parameter space
// at each (offset,energy) node,
//   print 1 if this node successfully found an alternate
//     value from a nearby node
//   print 0 otherwise
void VPsfBlock::formatBlock_wasinterpolated()
{
    printf( "formatBlock_wasinte     i_en:01234567890123456789012345678901\n" );
    for( unsigned int i_off = 0 ; i_off < fPsfBlock.size() ; i_off++ )
    {
        printf( "formatBlock_wasinte i_off=%2d ", i_off ) ;
        for( unsigned int i_en = 0 ; i_en < fPsfBlock[i_off].size() ; i_en++ )
        {
            printf( "%d", fPsfBlock[i_off][i_en].wasinterpolated ) ;
        }
        printf( "\n" );
    }
    printf( "\n" );
}

// print a 2d array of the parameter space
// at each (offset,energy) node,
//   print 1 if this node's 68 and 80% containment radii
//     are both above 0.0
//   print 0 otherwise
void VPsfBlock::formatBlock_zerorad()
{
    printf( "formatBlock_zerorad     i_en:01234567890123456789012345678901\n" );
    for( unsigned int i_off = 0 ; i_off < fPsfBlock.size() ; i_off++ )
    {
        printf( "formatBlock_zerorad i_off=%2d ", i_off ) ;
        for( unsigned int i_en = 0 ; i_en < fPsfBlock[i_off].size() ; i_en++ )
        {
            if( fPsfBlock[i_off][i_en].rad68 > 0.0 && fPsfBlock[i_off][i_en].rad80 > 0.0 )
            {
                printf( "1" );
            }
            else
            {
                printf( "0" );
            }
        }
        printf( "\n" );
    }
    printf( "\n" );
}

// print a 2d array of the parameter space
// at each (offset,energy) node,
//   print 1 if this node's sigma and gamma (its king function parameters)
//     are both above 0.0
//   print 0 otherwise
void VPsfBlock::formatBlock_sanefit()
{
    printf( "formatBlock_sanefit     i_en:01234567890123456789012345678901\n" );
    for( unsigned int i_off = 0 ; i_off < fPsfBlock.size() ; i_off++ )
    {
        printf( "formatBlock_sanefit i_off=%2d ", i_off ) ;
        for( unsigned int i_en = 0 ; i_en < fPsfBlock[i_off].size() ; i_en++ )
        {
            if( fPsfBlock[i_off][i_en].sigma > 0.0 && fPsfBlock[i_off][i_en].gamma > 0.0 )
            {
                printf( "1" );
            }
            else
            {
                printf( "0" );
            }
        }
        printf( "\n" );
    }
    printf( "\n" );
}


// find similarity rating for a target azimuth bin (0-8) and a comparison azimuth bin
// 0 [135,-135] 1 [-180,-90] 2 [-135,-45] 3 [-90,0] 4 [-45,45] 5 [0,90] 6 [45,135] 7 [90,180] 8 [-1000,1000]
float similarity_azbin( int target, int compare )
{
    float sim     = 100.0 ;
    
    // handle cases involving either azimuth bin being zero
    if( target == 8 && compare == 8 )
    {
        sim = 0.0 ;
    }
    if( target == 8 && compare != 8 )
    {
        sim = 1.0 ;
    }
    if( target != 8 && compare == 8 )
    {
        sim = 1.0 ;
    }
    if( target != 8 && compare != 8 )
    {
        // find angle between target and compare
        // if angle range is 360, equation is
        // sim = 180 - abs( abs(target-compare) - 180 ) ;
        // but since these are azimuth bins, they go from 1-8
        sim = 4.0 - fabs( fabs( ( float )target - ( float )compare ) - 4.0 ) ;
    }
    
    return sim ;
}

// find similarity between two wobbles
float similarity_wobble( float target, float compare )
{
    return fabs( target - compare ) / 0.25 ; // 0.25 = wobble bin width
}

// find similarity between two spectral indexes
float similarity_specind( float target, float compare )
{
    return fabs( target - compare ) / 0.1 ; // 0.1 = spectral index bin width
}

// find similarity between two noises (50-450)
float similarity_noise( int target, int compare )
{
    return fabs( ( float )target - ( float )compare ) / 60.0 ; // 60.0 = average noise bin width
}

// find similarity between two zenith angles 0-90
float similarity_zenith( int target, int compare )
{
    return fabs( ( float )target - ( float )compare ) / 5.0 ; // 5.0 = average zenith bin width
}


