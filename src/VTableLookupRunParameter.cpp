/*! \class VTableLookupRunParameter
    \brief parameter storage class

*/

#include "VTableLookupRunParameter.h"


VTableLookupRunParameter::VTableLookupRunParameter()
{
    fDebug = 0;

    outputfile = "";
    tablefile = "";
    ze = -1.;
    isMC = false;
    fUpdateInstrumentEpoch = true;
    fUseMedianEnergy = 1;
    fPE = false;
    fWriteTables = false;
    writeoption = "recreate";
    fMinRequiredShowerPerBin = 10.;
    bNoNoTrigger = true;
    fUseSelectedImagesOnly = true;
    bWriteReconstructedEventsOnly = 1;
    bShortTree = false;
    fWritePixelLists = false;
    bWriteMCPars = false;
    rec_method = 0;
    fWrite1DHistograms = false;
    fSpectralIndex = 2.0;
    fWobbleOffset = 500;     // integer of wobble offset * 100
    fNoiseLevel = 250;
    fTableFilling_useStereoMCParameter = false;
    fTableFillingCut_NImages_min = 2;
    fTableFillingCut_WobbleCut_max = 15.;
    fminsize = 0.;
    fmaxdist = 50000.;
    fmaxdistfraction = -1.;
    fmaxloss = 1.;
    fminfui = 0.;
    fSelectRandom = -1.;
    fSelectRandomSeed = 17;
    fRerunStereoReconstruction = false;
    fRerunStereoReconstruction_minAngle = -1.;
    fRerunStereoReconstruction_BDTNImages_max = 4;
    fRerunStereoReconstruction_BDTFileName = "";
    fEnergyReconstruction_BDTFileName = "";
    fCoreReconstruction_BDTFileName = "";
    fDispError_BDTFileName = "";
    fDispError_BDTWeight = 5.;
    fTelescopeList_sim_telarray_Counting = "";
    fTelescopeType_weightFile = "";
    fRunParameterFile = "";
    fQualityCutLevel = 0;

    fUsetimeGradientLookupTables = false;

    fMC_distance_to_cameracenter_min =  0.;
    fMC_distance_to_cameracenter_max =  1.e10;

    fEventSelectionCut_lossCutMax = 1.e9;
    fEventSelectionCut_distanceCutMax = 1.e9;

    fNentries = TChain::kBigNumber;
    fMaxRunTime = 1.e9;

    printpara = "";

    meanpedvars = 0.;
}


bool VTableLookupRunParameter::fillParameters( int argc, char* argv[] )
{
    // check number of command line parameters
    if( argc < 2 )
    {
        printHelp();
        return false;
    }
    // =============================================
    // reading command line parameters
    // =============================================
    // read command line parameters
    int i = 1;
    while( i++ < argc )
    {
        string iTemp = argv[i - 1];
        string iTemp2 = "";
        if( i < argc )
        {
            iTemp2 = argv[i];
        }
        if( iTemp.find( "-help" ) < iTemp.size() )
        {
            printHelp();
            return false;
        }
        if( ( iTemp.find( "-input" ) < iTemp.size() || iTemp.find( "-sourcefile" ) < iTemp.size() )
                && !( iTemp.find( "-inputfilelist" ) < iTemp.size() ) )
        {
            if( iTemp2.size() > 0 )
            {
                inputfile.push_back( iTemp2 );
                i++;
            }
        }
        else if( iTemp.find( "-inputfilelist" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fillInputFile_fromList( iTemp2 );
                i++;
            }
        }
        else if( iTemp.find( "-o" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                outputfile = iTemp2;
                i++;
            }
        }
        else if( iTemp.find( "printrunparameters" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                printpara = iTemp2;
                i++;
            }
            return true;
        }
        else if( iTemp.find( "useMedian" ) < iTemp.size() )
        {
            fUseMedianEnergy = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "updateEpoch" ) < iTemp.size() )
        {
            fUpdateInstrumentEpoch = ( bool )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "noise" ) < iTemp.size() )
        {
            fNoiseLevel = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "minshowerperbin" ) < iTemp.size() )
        {
            fMinRequiredShowerPerBin = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "woff" ) < iTemp.size() )
        {
            fWobbleOffset = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            // wobble offset is given in float
            if( fWobbleOffset < 10 )
            {
                fWobbleOffset = ( int )( atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) * 1000 + 0.5 );
            }
        }
        else if( iTemp.find( "-table" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                tablefile = iTemp2;
                i++;
            }
        }
        else if( iTemp.find( "-fill" ) < iTemp.size() )
        {
            int iT = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( iT == 1 )
            {
                fWriteTables = true;
            }
            else if( iT == 0 )
            {
                fWriteTables = false;
            }
            else
            {
                cout << "unknown parameter, choose 1=fill or 2=read lookup tables" << endl;
                return false;
            }
        }
        // rerun the stero reconstruction
        else if( iTemp.find( "-redo_stereo_reconstruction" ) < iTemp.size() )
        {
            fRerunStereoReconstruction = true;
        }
        // new minimum angle between image axes for simple stereo reconstruction
        else if( iTemp.find( "-minangle_stereo_reconstruction" ) < iTemp.size() )
        {
            fRerunStereoReconstruction_minAngle = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        // BDT directory and file name for disp stereo reconstruction (direction)
        else if( iTemp.find( "-tmva_filename_stereo_reconstruction" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRerunStereoReconstruction_BDTFileName = iTemp2;
                i++;
            }
        }
        // DISP BDT reconstruction is applied for images with up to this multiplicity
        else if( iTemp.find( "-tmva_nimages_max_stereo_reconstruction" ) < iTemp.size() )
        {
            fRerunStereoReconstruction_BDTNImages_max = ( unsigned int )( atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
            if( fRerunStereoReconstruction_BDTNImages_max > 40000 )
            {
                cout << "VTableLookupRunParameter::fillParameters() error:";
                cout << " maximum number of images for TMVA disp reconstruction is 4";
                cout << " (selection was " << fRerunStereoReconstruction_BDTNImages_max << ")" << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
        }
        // BDT directory and file name for disp stereo reconstruction (energy)
        else if( iTemp.find( "-tmva_filename_energy_reconstruction" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fEnergyReconstruction_BDTFileName = iTemp2;
                i++;
            }
        }
        // BDT directory and file name for disp stereo reconstruction (core)
        else if( iTemp.find( "-tmva_filename_core_reconstruction" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fCoreReconstruction_BDTFileName = iTemp2;
                i++;
            }
        }
        // BDT directory and file name for disp stereo reconstruction (disperror)
        else if( iTemp.find( "-tmva_filename_disperror_reconstruction" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fDispError_BDTFileName = iTemp2;
                i++;
            }
        }
        else if( iTemp.find( "-tmva_disperror_weight" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fDispError_BDTWeight = atof( iTemp2.c_str() );
                i++;
            }
        }
        // read in a list of telescopes (sim_telarray counting)
        else if( iTemp.find( "-sub_array_sim_telarray_counting" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fTelescopeList_sim_telarray_Counting = iTemp2;
                i++;
            }
        }
        // run parameters from file
        else if( iTemp.find( "-runparameter" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fRunParameterFile = iTemp2;
                i++;
            }
        }
        // telescope type dependent weights
        else if( iTemp.find( "-teltypeweightfile" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fTelescopeType_weightFile = iTemp2;
                i++;
            }
        }
        else if( iTemp.find( "-qualitycutlevel" ) < iTemp.size() )
        {
            if( iTemp2.size() > 0 )
            {
                fQualityCutLevel = ( unsigned int )( atoi( iTemp2.c_str() ) );
                i++;
            }
        }
        else if( iTemp.find( "-ze" ) < iTemp.size() )
        {
            ze = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-selectRandom" ) < iTemp.size() && !( iTemp.find( "-selectRandomSeed" ) < iTemp.size() ) )
        {
            fSelectRandom = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( fSelectRandom > 1. || fSelectRandom < 0. )
            {
                cout << "Error: probability has to be in [0,1]: " << fSelectRandom << endl;
                return false;
            }
        }
        else if( iTemp.find( "-selectRandomSeed" ) < iTemp.size() )
        {
            fSelectRandomSeed = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-loss_energyCuts" ) < iTemp.size() )
        {
            fEventSelectionCut_lossCutMax = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-distance_energyCuts" ) < iTemp.size() )
        {
            fEventSelectionCut_distanceCutMax = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-minImages" ) < iTemp.size() )
        {
            fTableFillingCut_NImages_min = ( unsigned int )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-spectralIndex" ) < iTemp.size() )
        {
            fSpectralIndex = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-maxdist" ) < iTemp.size()
                 && !( iTemp.find( "-maxdistancetocameracenter" ) < iTemp.size() )
                 && !( iTemp.find( "-maxdistfraction" ) < iTemp.size() ) )
        {
            fmaxdist = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-maxdistfraction" ) < iTemp.size() )
        {
            fmaxdistfraction = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-maxloss" ) < iTemp.size() )
        {
            fmaxloss = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-minfui" ) < iTemp.size() )
        {
            fminfui = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-maxdistancetocameracenter" ) < iTemp.size() )
        {
            fMC_distance_to_cameracenter_max  = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-mindistancetocameracenter" ) < iTemp.size() )
        {
            fMC_distance_to_cameracenter_min  = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            // looking at squared differences!!
            if( fMC_distance_to_cameracenter_min < 0. )
            {
                fMC_distance_to_cameracenter_min = 0.;
            }
        }
        else if( iTemp.find( "-CTAoffAxisBins" ) < iTemp.size() )
        {
            setCTA_MC_offaxisBins();
        }
        else if( iTemp.find( "-add_mc_spectral_index" ) < iTemp.size() )
        {
            fAddMC_spectral_index.push_back( atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() ) );
        }
        else if( iTemp.find( "-minsize" ) < iTemp.size() )
        {
            fminsize = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-debug" ) < iTemp.size() )
        {
            fDebug = ( unsigned int )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-arrayrecid" ) < iTemp.size() || iTemp.find( "-recid" ) < iTemp.size() )
        {
            rec_method = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-update" ) < iTemp.size() )
        {
            bool iT = ( bool )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            if( iT )
            {
                writeoption = "update";
            }
            else
            {
                writeoption = "recreate";
            }
        }
        else if( iTemp.find( "-noNo" ) < iTemp.size() )
        {
            bNoNoTrigger = false;
        }
        else if( iTemp.find( "-writeReconstructedEventsOnly" ) < iTemp.size() )
        {
            if( iTemp.rfind( "=" ) != string::npos )
            {
                bWriteReconstructedEventsOnly = atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
            }
            else
            {
                bWriteReconstructedEventsOnly = 0;
            }
        }
        else if( iTemp.find( "use_mc_parameters" ) < iTemp.size() )
        {
            fTableFilling_useStereoMCParameter = true;
        }
        else if( iTemp.find( "use_tgrad_tables" ) < iTemp.size() )
        {
            fUsetimeGradientLookupTables = true;
        }
        else if( iTemp.find( "-short" ) < iTemp.size() )
        {
            bShortTree = true;
        }
        else if( iTemp.find( "-pixellist" ) < iTemp.size() )
        {
            fWritePixelLists = true;
        }
        else if( iTemp.find( "-pe" ) < iTemp.size() )
        {
            fPE = true;
        }
        else if( iTemp.find( "-nomctree" ) < iTemp.size() )
        {
            bWriteMCPars = false;
        }
        else if( iTemp.find( "-mctree" ) < iTemp.size() )
        {
            bWriteMCPars = true;
        }
        else if( iTemp.find( "-write1DHistograms" ) < iTemp.size() )
        {
            fWrite1DHistograms = true;
        }
        else if( iTemp.find( "maxnevents" ) < iTemp.size() )
        {
            fNentries = ( Long64_t )atoi( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "maxruntime" ) < iTemp.size() )
        {
            fMaxRunTime = atof( iTemp.substr( iTemp.rfind( "=" ) + 1, iTemp.size() ).c_str() );
        }
        else if( iTemp.find( "-limitEnergyReconstruction" ) < iTemp.size() )
        {
            cout << "obsolete run parameter -limitEnergyReconstruction; ignored" << endl;
        }
        else
        {
            cout << "Error: unknown run parameter: " << iTemp << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    // filling of tables requires Monte Carlo
    if( fWriteTables )
    {
        isMC = true;
    }
    // =============================================
    // end of reading command line parameters
    // =============================================

    // require inputfile name
    if( inputfile.size() == 0 )
    {
        cout << "error: no input file" << endl;
        cout << "...exiting" << endl;
        return false;
    }
    // require table file
    if( tablefile.size() == 0 )
    {
        cout << "error: no lookup table file" << endl;
        cout << "...exiting" << endl;
        return false;
    }

    // set output file name (mainly for VTS analysis with a single inputfile)
    if( outputfile.size() == 0 && inputfile.size() == 1 )
    {
        // wildcards for input file
        if( inputfile[0].find( "*" ) < inputfile[0].size() )
        {
            outputfile = "mscw.root";
        }
        // no wildcards for input file
        else
        {
            outputfile = inputfile[0].substr( 0, inputfile[0].rfind( "." ) );
            outputfile += ".mscw.root";
        }
    }
    // run parameters from file (optional)
    if( fRunParameterFile.size() > 0 )
    {
        if( !readRunParameters( fRunParameterFile ) )
        {
            exit( EXIT_FAILURE );
        }
    }
    // read telescope type dependent weights (optional)
    if( fTelescopeType_weightFile.size() > 0 )
    {
        readTelTypeDepdendentWeights( fTelescopeType_weightFile );
    }
    // for VTS analysis with a single inputfile: get telescope combinations
    if( fTelescopeList_sim_telarray_Counting.size() == 0 )
    {
        if( !readTelescopeToAnalyze( inputfile[0] ) )
        {
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    // CTA analysis might provide additional a list of telescope to analyse
    else if( fTelescopeList_sim_telarray_Counting.size() > 0 && inputfile[0].find( "*" ) == string::npos )
    {
        if( !readTelescopeToAnalyze( fTelescopeList_sim_telarray_Counting, inputfile[0] ) )
        {
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    else if( inputfile[0].find( "*" ) == string::npos )
    {
        cout << "Error: unable to read list of telescopes." << endl;
        cout << "Provide a list with command line parameter -sub_array_sim_telarray_counting <telescope list>" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }

    fillTelescopeTypeDependentWeights();


    return true;
}

/*
 *  file telescope type dependent weights
 *
 *  weights from external file overwrites the weights from the telescope lists
 *
*/
void VTableLookupRunParameter::fillTelescopeTypeDependentWeights()
{
    if( fTelescopeType_weight.size() > 0 )
    {
        for( unsigned int i = 0; i < fTelToAnalyzeData.size(); i++ )
        {
            if( fTelToAnalyzeData[i]
                    && fTelescopeType_weight.find( fTelToAnalyzeData[i]->fTelType ) != fTelescopeType_weight.end() )
            {
                fTelToAnalyzeData[i]->fWeight = fTelescopeType_weight[fTelToAnalyzeData[i]->fTelType];
            }
        }
    }
}


/*
 * read runparameters from a parameter file
 *
*/
bool VTableLookupRunParameter::readRunParameters( string iFile )
{

    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VTableLookupRunParameter::readRunParameters() error:";
        cout << "file with runparameters not found, ";
        cout << iFile << endl;
        return false;
    }
    string iLine;
    string iS1;
    string iS2;

    cout << "Reading run parameters from " << iFile << endl;
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            istringstream is_stream( iLine );
            is_stream >> iS1;
            // a valid line starts with a '*'
            if( iS1 != "*" )
            {
                continue;
            }
            is_stream >> iS2;
            if( VUtilities::lowerCase( iS2 ) == "maxloss" )
            {
                is_stream >> fEventSelectionCut_lossCutMax;
            }
            else if( VUtilities::lowerCase( iS2 ) == "maxdist" )
            {
                is_stream >> fmaxdist;
            }
            else if( VUtilities::lowerCase( iS2 ) == "minfui" )
            {
                is_stream >> fminfui;
            }
            else if( VUtilities::lowerCase( iS2 ) == "maxdist_energy" )
            {
                is_stream >> fEventSelectionCut_distanceCutMax;;
            }
            else if( VUtilities::lowerCase( iS2 ) == "median_energy" )
            {
                is_stream >> fUseMedianEnergy;
            }
            else if( VUtilities::lowerCase( iS2 ) == "minangle_stereoreconstruction" )
            {
                is_stream >> fRerunStereoReconstruction_minAngle;
            }
        }
    }

    return true;
}


/*
 * read telescope type dependent weights (optional)
 *
*/
bool VTableLookupRunParameter::readTelTypeDepdendentWeights( string iFile )
{
    fTelescopeType_weight.clear();

    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VTableLookupRunParameter::readTelTypeDepdendentWeights error:";
        cout << "file with telescope type dependent weights not found: ";
        cout << iFile << endl;
        return false;
    }
    string iLine;
    string iS1;
    ULong64_t iT1 = 0;
    double iT2 = 0.;

    cout << "Reading telescope type dependent weights from " << iFile << endl;
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            istringstream is_stream( iLine );
            is_stream >> iS1;
            // a valid line starts with a '*'
            if( iS1 != "*" )
            {
                continue;
            }
            is_stream >> iT1;
            is_stream >> iT2;
            fTelescopeType_weight[iT1] = iT2;
        }
    }

    // telescope type dependent weights
    if( fTelescopeType_weight.size() > 0 )
    {
        cout << "Telescope type dependent weights: " << endl;
        for( map< ULong64_t, double >::iterator iT  = fTelescopeType_weight.begin();
                iT != fTelescopeType_weight.end();
                ++iT )
        {
            cout << "\t weight for telescope type " << iT->first << ": ";
            cout << iT->second << endl;
        }
    }


    return true;
}

/*
 * read telescope combinations for analysis
 *
 * sim_telarray counting
 *
 */
bool VTableLookupRunParameter::readTelescopeToAnalyze( string iTelescopeList_sim_telarray_Counting,
        string iEvndispRootFile )
{
    fTelToAnalyzeData.clear();
    // (this vector will have later the length of the number
    //  of telescopes of the full array)

    ////////////////////////////////////////////////////////
    // read list of all telescopes from first evndisp root file

    TFile iF( iEvndispRootFile.c_str() );
    if( iF.IsZombie() )
    {
        cout << "VTableLookupRunParameter::readTelescopeToAnalyze warning: could not open input file to read run parameters" << endl;
        cout << "\t " << iEvndispRootFile << endl;
        return false;
    }
    TTree* iT = ( TTree* )iF.Get( "telconfig" );
    if( !iT )
    {
        cout << "VTableLookupRunParameter::readTelescopeToAnalyze error reading telconfig tree" << endl;
        cout << "\t " << iEvndispRootFile << endl;
        return false;
    }
    cout << "reading telescope configuration from " << iEvndispRootFile << endl;
    Ctelconfig* itelconfig = new Ctelconfig( iT );
    itelconfig->GetEntry( 0 );
    for( unsigned int i = 0; i < itelconfig->fChain->GetEntries(); i++ )
    {
        itelconfig->GetEntry( i );

        fTelToAnalyzeData.push_back( new VTableLookupTelToAnalyze() );
        // note: by default are all telescopes off!!
        fTelToAnalyzeData.back()->fTelToAnalyze     = 0;
        fTelToAnalyzeData.back()->fTelID            = itelconfig->TelID;
        fTelToAnalyzeData.back()->fTelID_hyperArray = itelconfig->TelID_hyperArray;
        fTelToAnalyzeData.back()->fTelType          = itelconfig->TelType;
        fTelToAnalyzeData.back()->fWeight           = 1.;
    }

    ////////////////////////////////////////////////////////
    // read list of telescope to analyse
    ifstream is;
    is.open( iTelescopeList_sim_telarray_Counting.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VTableLookupRunParameter::readTelescopeToAnalyze error: file with telescope list not found: ";
        cout << iTelescopeList_sim_telarray_Counting << endl;
        return false;
    }
    string iLine;
    unsigned int iT1 = 0;
    double iW = 1;
    string iT2;
    cout << "reading sub-array configuration from " << iTelescopeList_sim_telarray_Counting << endl;

    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            iW = 1.;
            unsigned int z = 0;
            istringstream is_stream( iLine );

            while( is_stream >> iT2 )
            {
                if( z == 0 )
                {
                    iT1 = atoi( iT2.c_str() );
                }
                // expected weighting for stereo reconstruction
                // in column 5
                else if( z == 5 )
                {
                    iW = atof( iT2.c_str() );
                }
                z++;
            }
            // find this telescope and switch it on
            for( unsigned int t = 0; t < fTelToAnalyzeData.size(); t++ )
            {
                if( fTelToAnalyzeData[t]
                        && fTelToAnalyzeData[t]->fTelID_hyperArray == iT1 )
                {
                    fTelToAnalyzeData[t]->fTelToAnalyze = 1;
                    fTelToAnalyzeData[t]->fWeight = iW;
                }
            }
        }
    }
    is.close();

    return true;
}

/*
 * read telescope combination for analysis
 *
 * this works only for a small number of telescopes (<10)
 */
bool VTableLookupRunParameter::readTelescopeToAnalyze( string iEvndispRootFile )
{
    fTelToAnalyse.clear();

    // use chain to get list of files
    TChain iTel( "telconfig" );
    int iNFil = iTel.Add( iEvndispRootFile.c_str() );
    if( iNFil < 1 )
    {
        return false;
    }
    TFile* iF = iTel.GetFile();
    if( !iF || iF->IsZombie() )
    {
        cout << "VTableLookupRunParameter::readTelescopeToAnalyze warning: could not open input file to read run parameters" << endl;
        cout << "\t " << iEvndispRootFile.c_str() << endl;
        cout << "\t file is possibly empty?" << endl;
        return false;
    }

    // read telescopes to analyse from eventdisplay run parameter list
    vector< unsigned int > iRunParT;
    VEvndispRunParameter* iPar = ( VEvndispRunParameter* )iF->Get( "runparameterV2" );
    if( iPar )
    {
        iRunParT = iPar->fTelToAnalyze;
        if( iPar->getObservatory().find( "VERITAS" ) == string::npos )
        {
            cout << "VTableLookupRunParameter::readTelescopeToAnalyze warning: ";
            cout << "reading without telescope lists not enabled for non-VERITAS observatories";
            cout << endl;
            return false;
        }
    }
    else
    {
        cout << "VTableLookupRunParameter::readTelescopeToAnalyze warning: could not find evndisp run parameters (runparameterV2)" << endl;
        return false;
    }
    // cross check if one of the has been switched off in the analysis
    VEvndispReconstructionParameter* iRecPar = ( VEvndispReconstructionParameter* )iF->Get( "EvndispReconstructionParameter" );
    if( iRecPar )
    {
        // this works only if number of telescopes = number of telescope types
        if( iRecPar->getReconstructionParameterData( rec_method )
                && iPar->fNTelescopes == iRecPar->getReconstructionParameterData( rec_method )->fLocalUseImage.size() )
        {
            for( unsigned int i = 0; i < iRunParT.size(); i++ )
            {

                if( iRunParT[i] < iRecPar->getReconstructionParameterData( rec_method )->fLocalUseImage.size()
                        && iRecPar->getReconstructionParameterData( rec_method )->fLocalUseImage[iRunParT[i]] )
                {
                    fTelToAnalyse.push_back( iRunParT[i] );
                }
            }
            cout << "Following telescopes are included in the analysis: ";
            for( unsigned int i = 0; i < fTelToAnalyse.size(); i++ )
            {
                cout << " T" << fTelToAnalyse[i] + 1;
            }
            cout << endl;
        }
    }
    else
    {
        cout << "VTableLookupRunParameter::readTelescopeToAnalyze warning: ";
        cout << "could not find evndisp reconstruction parameters (EvndispReconstructionParameter)" << endl;
        return false;
    }

    // fill secondary list used in many analysis part
    fTelToAnalyzeData.clear();
    for( unsigned int i = 0; i < iPar->fNTelescopes; i++ )
    {
        fTelToAnalyzeData.push_back( new VTableLookupTelToAnalyze() );
        fTelToAnalyzeData.back()->fTelToAnalyze     = 0;
        fTelToAnalyzeData.back()->fTelID            = i + 1;
        fTelToAnalyzeData.back()->fTelID_hyperArray = i + 1;
        fTelToAnalyzeData.back()->fTelType          = i + 1;
        fTelToAnalyzeData.back()->fWeight           = 1.;
    }
    for( unsigned int j = 0; j < fTelToAnalyse.size(); j++ )
    {
        if( fTelToAnalyse[j] < fTelToAnalyzeData.size() )
        {
            fTelToAnalyzeData[fTelToAnalyse[j]]->fTelToAnalyze = 1;
        }
    }

    return true;
}


void VTableLookupRunParameter::printHelp()
{
    if( gSystem->Getenv( "EVNDISPSYS" ) )
    {
        int syst_ret = system( "cat $EVNDISPSYS/README/README.MSCW_ENERGY" );
        if( syst_ret == -1 )
        {
            cout << "VTableLookupRunParameter::printHelp() error: could not find helper file in README directory" << endl;
        }
    }
    else
    {
        cout << "VTableLookupRunParameter::printHelp() no help available (environmental variable EVNDISPSYS not set)" << endl;
    }
    return;
}


void VTableLookupRunParameter::print( int iP )
{
    cout << "mscw_energy VERSION " << getEVNDISP_VERSION() << endl;
    cout << endl;
    cout << "debug level " << fDebug << endl;
    cout << "lookuptable: " << tablefile << endl;
    cout << endl;
    cout << "evndisp reconstruction parameter ID: " << rec_method << endl;
    cout << endl;
    printCTA_MC_offaxisBins();
    cout << endl;
    // (do not print more than 10 input files)
    if( inputfile.size() < 10 )
    {
        cout << "input file(s): " << endl;
        for( unsigned int i = 0; i < inputfile.size(); i++ )
        {
            cout << "\t" << inputfile[i] << endl;
        }
    }
    else
    {
        cout << "input files (print only first 10 of " << inputfile.size() << " files): " << endl;
        for( unsigned int i = 0; i < 10; i++ )
        {
            cout << "\t" << inputfile[i] << endl;
        }
        cout << "\t ..." << endl;
    }

    if( isMC )
    {
        cout << " (input data is MC)";
    }
    if( fPE )
    {
        cout << " (input data is PE)";
    }
    cout << endl;
    if( !fWriteTables )
    {
        cout << "output file: " << outputfile << endl;
        if( bWriteReconstructedEventsOnly >= 0 )
        {
            cout << "writing reconstructed events only (" << bWriteReconstructedEventsOnly << ")" << endl;
        }
    }
    else
    {
        cout << "filling lookup tables for: ";
        cout << " zenith " << ze << ", direction offset " << fWobbleOffset << "(x0.01) [deg], ";
        cout << "noise level " << fNoiseLevel << ", spectral index " << fSpectralIndex << endl;
        if( fWrite1DHistograms )
        {
            cout << "write 1D histograms to disk" << endl;
        }
        cout << "\t minimum telescope multiplicity: " << fTableFillingCut_NImages_min << endl;
        cout << "\t distance to camera: > " << fMC_distance_to_cameracenter_min << " [deg], <";
        cout << fMC_distance_to_cameracenter_max << " [deg]" << endl;
        if( fTableFilling_useStereoMCParameter )
        {
            cout << "\t using MC core and direction for table filling " << endl;
        }
    }
    if( fRerunStereoReconstruction )
    {
        cout << "\t rerunning stereo reconstruction" << endl;
        if( fRerunStereoReconstruction_BDTFileName.size() > 0 )
        {
            cout << "\t reading BDT TMVA files from " << fRerunStereoReconstruction_BDTFileName << endl;
            cout << "\t BDT TMVA stereo reconstruction is applied for events with <= ";
            cout << fRerunStereoReconstruction_BDTNImages_max << " images" << endl;
            if( fmaxdist < 1.e3 )
            {
                cout << "\t BDT TMVA stereo reconstruction distance cut < " << fmaxdist << endl;
            }
            if( fmaxdistfraction > 0. )
            {
                cout << "\t BDT TMVA stereo reconstruction distance cut (fraction of FOV) < ";
                cout << fmaxdistfraction << endl;
            }
            if( fmaxloss < 1. )
            {
                cout << "\t BDT TMVA stereo reconstruction loss cut < " << fmaxloss << endl;
            }
            if( fminfui > 0. )
            {
                cout << "\t BDT TMVA stereo reconstruction fui cut < " << fminfui << endl;
            }
        }
    }
    cout << "\t maximum loss for energy calculation: " << fEventSelectionCut_lossCutMax << endl;
    cout << "\t maximum distance for energy calculation: " << fEventSelectionCut_distanceCutMax << endl;
    if( iP == 2 && isMC )
    {
        cout << "zenith angle " << ze << " [deg], wobble offset " << fWobbleOffset / 100. << " [deg], noise level " << fNoiseLevel << endl;
    }
    if( fSelectRandom > 0. )
    {
        cout << "random event selection: " << fSelectRandom << ", seed:" << fSelectRandomSeed << endl;
    }
    if( fUseSelectedImagesOnly )
    {
        cout << "\t use evndisp image selection" << endl;
    }
    else
    {
        cout << "\t use all images" << endl;
    }
    if( fWriteTables )
    {
        cout << "minimum number of showers required per lookup table bin: " << fMinRequiredShowerPerBin << endl;
    }
    if( fUseMedianEnergy == 1 )
    {
        cout << "use median of energy distributions" << endl;
    }
    else if( fUseMedianEnergy == 2 )
    {
        cout << "use median+mpv of energy distributions" << endl;
    }
    else
    {
        cout << "use mean of energy distributions" << endl;
    }
    if( fUpdateInstrumentEpoch )
    {
        cout << "updating instrument epoch from default epoch file" << endl;
    }

    if( iP >= 1 )
    {
        cout << endl;
        if( meanpedvars > 0. )
        {
            cout << "mean pedvars: " << meanpedvars << endl;
            cout << "mean pedvars per telescope: ";
            for( unsigned int i = 0; i < pedvars.size(); i++ )
            {
                cout << pedvars[i] << "/";
            }
            cout << endl;
        }
        else
        {
            cout << "no pedvar information available" << endl;
        }
    }
}

bool VTableLookupRunParameter::fillInputFile_fromList( string iList )
{
    ifstream is;
    is.open( iList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VTableLookupRunParameter::fillInputFile_fromList() error reading list of input files: " << endl;
        cout << iList << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "VTableLookupRunParameter::fillInputFile_fromList() reading input file list: " << endl;
    cout << iList << endl;
    string iLine;
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            inputfile.push_back( iLine );
        }
    }
    is.close();

    cout << "total number of input files " << inputfile.size() << endl;

    return true;
}

void VTableLookupRunParameter::setCTA_MC_offaxisBins()
{
    fCTA_MC_offaxisBin_min.clear();
    fCTA_MC_offaxisBin_max.clear();

    /* New binning definition 2015-11-10
     *
     * choose smaller number of bins to increase statistics
     * per off-axis bins
     *
     */

    fCTA_MC_offaxisBin_min.push_back( 0.0 );
    fCTA_MC_offaxisBin_max.push_back( 1.0 );

    fCTA_MC_offaxisBin_min.push_back( 1.0 );
    fCTA_MC_offaxisBin_max.push_back( 2.0 );

    fCTA_MC_offaxisBin_min.push_back( 2.0 );
    fCTA_MC_offaxisBin_max.push_back( 3.0 );

    fCTA_MC_offaxisBin_min.push_back( 3.0 );
    fCTA_MC_offaxisBin_max.push_back( 4.0 );

    fCTA_MC_offaxisBin_min.push_back( 4.0 );
    fCTA_MC_offaxisBin_max.push_back( 5.0 );

    fCTA_MC_offaxisBin_min.push_back( 5.0 );
    fCTA_MC_offaxisBin_max.push_back( 6.0 );

    /* Binning definition until 2015-11-10

    	fCTA_MC_offaxisBin_min.push_back( 0.0 );
    	fCTA_MC_offaxisBin_max.push_back( 1.0 );

    	fCTA_MC_offaxisBin_min.push_back( 1.0 );
    	fCTA_MC_offaxisBin_max.push_back( 2.0 );

    	fCTA_MC_offaxisBin_min.push_back( 2.0 );
    	fCTA_MC_offaxisBin_max.push_back( 3.0 );

    	fCTA_MC_offaxisBin_min.push_back( 3.0 );
    	fCTA_MC_offaxisBin_max.push_back( 3.5 );

    	fCTA_MC_offaxisBin_min.push_back( 3.5 );
    	fCTA_MC_offaxisBin_max.push_back( 4.0 );

    	fCTA_MC_offaxisBin_min.push_back( 4.0 );
    	fCTA_MC_offaxisBin_max.push_back( 4.5 );

    	fCTA_MC_offaxisBin_min.push_back( 4.5 );
    	fCTA_MC_offaxisBin_max.push_back( 5.0 );

    	fCTA_MC_offaxisBin_min.push_back( 5.0 );
    	fCTA_MC_offaxisBin_max.push_back( 5.5 );

    	fCTA_MC_offaxisBin_min.push_back( 5.5 );
    	fCTA_MC_offaxisBin_max.push_back( 6.0 );

    */
}

void VTableLookupRunParameter::printCTA_MC_offaxisBins()
{
    if( fCTA_MC_offaxisBin_min.size() == 0 )
    {
        return;
    }

    cout << "setting the following off-axis bins for CTA analysis: " << endl;
    for( unsigned int i = 0; i < fCTA_MC_offaxisBin_min.size(); i++ )
    {
        cout << "   bin " << i << "\t min " << fCTA_MC_offaxisBin_min[i] << " deg, max " << fCTA_MC_offaxisBin_max[i] << " deg" << endl;
    }
}


/////////////////////////////////////////////////////////////////////////////
// data class for list of telescopes to be used in the analysis
VTableLookupTelToAnalyze::VTableLookupTelToAnalyze()
{
    fTelToAnalyze     = 1;
    fTelID            = 0;
    fTelID_hyperArray = 0;
    fTelType          = 0;
    fWeight           = 1.;
}

void VTableLookupTelToAnalyze::print()
{
    cout << "\t telescope " << fTelID << " (hyperarray ID: " << fTelID_hyperArray << ", type " << fTelType << ")";
    cout << ", weight: " << fWeight;
    cout << endl;
}
