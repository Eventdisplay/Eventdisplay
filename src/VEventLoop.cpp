/*!  \class VEventLoop
  \brief  main event loop, steering of analysis and event display

*/

#include "VEventLoop.h"

//! standard constructor
/*!
    \param irunparameter Pointer to run parameters (from command line or config file)
*/
VEventLoop::VEventLoop( VEvndispRunParameter* irunparameter )
{
    fDebug = irunparameter->fDebug;
    if( fDebug )
    {
        cout << "VEventLoop::VEventLoop()" << endl;
    }
    
    fRunPar = irunparameter;
    
    // total number of telescopes
    fNTel = fRunPar->fNTelescopes;
    if( fNTel < 1 )
    {
        cout << "VEventLoop::VEventLoop error: no telescopes defined" << endl;
        exit( EXIT_FAILURE );
    }
    if( fNTel >= 10 )
    {
        fRunPar->fPrintSmallArray = false;
    }
    
    // data readers
    fReader = 0;                                  // this pointer is used in the program for accesing any data source (raw or MC)
#ifndef NOVBF
    fRawDataReader = 0;
#endif
    fGrIsuReader = 0;
    fDSTReader = 0;
    fPEReader = 0;
    
    bMCSetAtmosphericID = false;
    fBoolPrintSample.assign( fNTel, true );
    fGPSClockWarnings.assign( fNTel, 0 );
    fTimeCut_RunStartSeconds = 0;
    
    fAnalyzeMode = true;
    fRunMode = ( E_runmode )fRunPar->frunmode;
    fEventNumber = 0;
    fNextEventStatus = false;
    fTimeCutsfNextEventStatus = true;
    fEndCalibrationRunNow = false;
    fBoolSumWindowChangeWarning = 0;
    fLowGainMultiplierWarning = 0;
    
    setRunNumber( fRunPar->frunnumber );
    
    // get the detector settings from the configuration files and set the cameras
    bool iMakeNeighbourList = setDetectorGeometry( fNTel, fRunPar->fcamera, fRunPar->getDirectory_EVNDISPDetectorGeometry() );
    
    // set which telescopes should be analyzed
    setTeltoAna( fRunPar->fTelToAnalyze );
    
    // set dead channel text
    setDeadChannelText();
    
    // read reconstruction parameters
    if( !get_reconstruction_parameters( fRunPar->freconstructionparameterfile, iMakeNeighbourList ) )
    {
        cout << "VEventLoop error while reading file with reconstruction parameters:" << endl;
        cout << fRunPar->freconstructionparameterfile << endl;
        exit( EXIT_FAILURE );
    }
    
    // set tracehandler
    fTraceHandler = new VTraceHandler();
    if( getRunParameter()->fTraceIntegrationMethod.size() > 0 )
    {
        fTraceHandler->setTraceIntegrationmethod( getRunParameter()->fTraceIntegrationMethod[0] );
    }
    fTraceHandler->setMC_FADCTraceStart( getRunParameter()->fMC_FADCTraceStart );
    fTraceHandler->setPulseTimingLevels( getRunParameter()->fpulsetiminglevels );
    
    // initialize calibrator (one for all telescopes)
    if( fCalibrated.size() == 0 ) for( unsigned int i = 0; i < fNTel; i++ )
        {
            fCalibrated.push_back( false );
        }
    fCalibrator = new VCalibrator();
    
    // create data summarizer
    fDST = 0;
    if( fRunMode == R_DST )
    {
        fDST = new VDST( ( fRunMode == R_DST ), ( fRunPar->fsourcetype == 1 || fRunPar->fsourcetype == 2 || fRunPar->fsourcetype == 6 ) );
    }
    
    // create analyzer (one for all telescopes)
    fAnalyzer = new VImageAnalyzer();
    // create new pedestal calculator
    fPedestalCalculator = new VPedestalCalculator();
    
    // create array analyzer
    fArrayAnalyzer = new VArrayAnalyzer();
    
    // create dead time calculator
    fDeadTime = new VDeadTime();
    fDeadTime->defineHistograms( fRunPar->fRunDuration );
    
    // create dead pixel organizer
    if( fRunPar->fSaveDeadPixelRegistry )
    {
        fDeadPixelOrganizer = new VDeadPixelOrganizer(
            fRunPar->fNTelescopes, fDetectorGeo->getNumChannels(), fDetectorGeo,
            fRunPar->fDBDataStartTimeMJD, fRunPar->fDBDataStartTimeSecOfDay,
            fRunPar->fDBDataStoppTimeMJD, fRunPar->fDBDataStoppTimeSecOfDay,
            "deadPixelRegistry", getRunNumber() ) ;
    }
    else
    {
        fDeadPixelOrganizer = 0 ;
    }
    
#ifndef NOGSL
    // FROGS
    if( fRunPar->ffrogsmode )
    {
        fFrogs = new VFrogs();
    }
#endif
    // Model3D
    if( fRunPar->fUseModel3D )
    {
        fModel3D = new VModel3D();
    }
    // reset cut strings and variables
    resetRunOptions();
}


VEventLoop::~VEventLoop()
{
}


/*!
    print basic run infos (file name, runnumber, etc.) to standard output
*/
void VEventLoop::printRunInfos()
{
    if( fDebug )
    {
        cout << "VEventLoop::printRunInfos()" << endl;
    }
    // astronometry type
    cout << "Using " << VAstronometry::getAstronometryLibrary() << " for astronometry" << endl << endl;
    // telescope parameters
    cout << endl << "Analysis parameters: " << endl;
    for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
    {
        setTelID( fRunPar->fTelToAnalyze[i] );
        
        cout << "Telescope " << fRunPar->fTelToAnalyze[i] + 1;
        if( i < getDetectorGeometry()->getTelType().size() )
        {
            cout << " (type " << getDetectorGeometry()->getTelType()[i] << ")";
        }
        cout << endl;
        if( fRunPar->fTraceIntegrationMethod[fRunPar->fTelToAnalyze[i]] )
        {
            cout << "\t trace integration method: \t" << fRunPar->fTraceIntegrationMethod[fRunPar->fTelToAnalyze[i]];
            if( fRunPar->fDoublePass[fRunPar->fTelToAnalyze[i]] )
            {
                cout << "  (doublepass, integration method pass 1: " << fRunPar->fTraceIntegrationMethod_pass1[fRunPar->fTelToAnalyze[i]] << ")";
            }
            cout << endl;
            if( fRunPar->fDF_DigitalFilter[fRunPar->fTelToAnalyze[i]] != 0 )
            {
                cout << "\t digital filter:\t\t upsample: " << fRunPar->fDF_UpSample[fRunPar->fTelToAnalyze[i]];
                cout << ", pole-zero cancellation parameter: " << fRunPar->fDF_PoleZero[fRunPar->fTelToAnalyze[i]];
                cout << endl;
            }
            cout << "\t start of summation window: \t" << fRunPar->fsumfirst[fRunPar->fTelToAnalyze[i]];
            cout << "\t(shifted by " << fRunPar->fTraceWindowShift[fRunPar->fTelToAnalyze[i]] << " samples)" << endl;
            cout << "\t length of summation window: \t" << fRunPar->fsumwindow_1[fRunPar->fTelToAnalyze[i]];
            cout << "/" << fRunPar->fsumwindow_2[fRunPar->fTelToAnalyze[i]];
            if( fRunPar->fDoublePass[fRunPar->fTelToAnalyze[i]] )
            {
                cout << "\t length of first pass summation window (double pass): \t" << fRunPar->fsumwindow_pass1[fRunPar->fTelToAnalyze[i]];
            }
            cout << endl;
        }
        else
        {
            cout << "\t no trace integration ";
            cout << "[pulse timing method " << fRunPar->fsumfirst_startingMethod[fRunPar->fTelToAnalyze[i]] << "]" << endl;
        }
        // print image cleaning parameters
        if( getImageCleaningParameter( false ) )
        {
            getImageCleaningParameter( false )->print();
        }
        if( isDoublePass() && getImageCleaningParameter( true ) )
        {
            getImageCleaningParameter( true )->print();
        }
        if( getCalData()->getLowGainMultiplierDistribution() && getCalData()->getLowGainMultiplierDistribution()->GetEntries() > 0 )
        {
            cout << "\t low gain multiplier: \t" << setprecision( 3 ) << getCalData()->getLowGainMultiplierDistribution()->GetMean();
            if( getCalData()->getLowGainMultiplierDistribution()->GetRMS() > 1.e-3 )
            {
                cout << "+-" << getCalData()->getLowGainMultiplierDistribution()->GetRMS();
            }
        }
        else
        {
            cout << "\t (no low gain multiplier distributions)";
        }
        cout << endl;
        if( TMath::Abs( fRunPar->fGainCorrection[fRunPar->fTelToAnalyze[i]] ) - 1. > 1.e-2 )
        {
            cout << "\t additional gain correction: " << fRunPar->fGainCorrection[fRunPar->fTelToAnalyze[i]];
        }
        cout << "\t LL edge fit: \t\t " << fRunPar->fLogLikelihoodLoss_min[i] << " < loss < " << fRunPar->fLogLikelihoodLoss_max[i];
        cout << "\t ntubes > " << fRunPar->fLogLikelihood_Ntubes_min[i] << endl;
    }
}


/*!
   reset all run parameters, open files, set eventnumbers, etc.

  handle with file names, get run numbers, set analyzer default values
*/
bool VEventLoop::initEventLoop()
{
    return initEventLoop( fRunPar->fsourcefile );
}


/*!
   reset all run parameters, open files, set eventnumbers, etc.

   handle with file names, get run numbers, set analyzer default values, read calibration files

   \param iFileName data source file
*/
bool VEventLoop::initEventLoop( string iFileName )
{
    if( fDebug )
    {
        cout << "VEventLoop::initEventLoop()" << endl;
    }
    fRunPar->fsourcefile = iFileName;
    fEventNumber = 0;
    
    // check if file exists (bizzare return value)
    if( gSystem->AccessPathName( iFileName.c_str() ) && fRunPar->fsourcetype != 5 )
    {
        cout << endl;
        cout << "VEventLoop::initEventLoop error; sourcefile not found: |" << iFileName << "|" << endl;
        exit( EXIT_FAILURE );
    }
    
    // set the data readers and open data files
    //     (different file formats)
#ifndef NOVBF
    try
    {
        // ============================
        // sourcefile has raw data format (prototype or vbf)
        if( fRunPar->fsourcetype == 0 || fRunPar->fsourcetype == 2 || fRunPar->fsourcetype == 3 )
        {
            if( fRawDataReader != 0 )
            {
                delete fRawDataReader;
            }
            if( fRunPar->fsourcetype == 0 )
            {
                fRawDataReader = new VRawDataReader( fRunPar->fsourcefile, fRunPar->fsourcetype, fRunPar->fNTelescopes, fDebug );
            }
            else
            {
                fRawDataReader = new VBFDataReader( fRunPar->fsourcefile, fRunPar->fsourcetype, fRunPar->fNTelescopes, fDebug, fRunPar->fPrintGrisuHeader );
                /////////////////////////////////////////////////////////////////////
                // open temporary file (do make sure that event numbering is correct)
                // get number of samples
                VBFDataReader i_tempReader( fRunPar->fsourcefile, fRunPar->fsourcetype, fRunPar->fNTelescopes, fDebug );
                // loop over several events
                // (note: MC might have no telescope event data until the first triggered event)
                // (note: zero suppressed data)
                vector< bool > i_nSampleSet( fRunPar->fTelToAnalyze.size(), false );
                unsigned int i_counter = 0;
                for( ;; )
                {
                    i_tempReader.getNextEvent();
                    for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
                    {
                        i_tempReader.setTelescopeID( fRunPar->fTelToAnalyze[i] );
                        // set number of samples
                        if( i_tempReader.getNumSamples() > 4 )
                        {
                            setNSamples( fRunPar->fTelToAnalyze[i], i_tempReader.getNumSamples() );
                            i_nSampleSet[i] = true;
                        }
                    }
                    unsigned int z = 0;
                    for( unsigned int i = 0; i < i_nSampleSet.size(); i++ )
                    {
                        if( i_nSampleSet[i] )
                        {
                            z++;
                        }
                    }
                    // found samples for all telescopes
                    if( z == i_nSampleSet.size() )
                    {
                        // require that sample length is the same in all telescopes
                        // (note: this is a strong requirement, but should be ok for
                        //  VTS (VBF is assumed to be VTS)
                        unsigned int i_tempSampeLength = 0;
                        bool iAll_are_the_same = true;
                        if( fRunPar->fTelToAnalyze.size() > 0 )
                        {
                            i_tempSampeLength = getNSamples( fRunPar->fTelToAnalyze[0] );
                        }
                        for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
                        {
                            if( i_tempSampeLength != getNSamples( fRunPar->fTelToAnalyze[i] ) )
                            {
                                iAll_are_the_same = false;
                            }
                        }
                        if( iAll_are_the_same )
                        {
                            break;
                        }
                    }
                    i_counter++;
                    if( i_counter == 1000 )
                    {
                        cout << "VEventLoop warning: could not find number of samples in the first 1000 events";
                        cout << " (this is normal for runs with event losses at the beginning)" << endl;
                    }
                    if( i_counter > 999999 )
                    {
                        cout << "VEventLoop warning: could not find number of samples in the first 999999 events" << endl;
                        break;
                    }
                }
                if( getNSamples() == 0 || i_counter > 99999 )
                {
                    cout << "VEventLoop::initEventLoop error: could not find any telescope events to determine sample length" << endl;
                    cout << "exiting..." << endl;
                    exit( EXIT_FAILURE );
                }
                cout << "Found consistent number of samples in vbf file after " << i_counter + 1 << " event(s): ";
                for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
                {
                    cout << "T" << fRunPar->fTelToAnalyze[i] + 1 << ": " << getNSamples( fRunPar->fTelToAnalyze[i] ) << " ";
                }
                cout << endl;
                ///////////////////////////////////////////////////////////////
            }
            // sourcefile is MC vbf file; noise is read from separate file
            if( fRawDataReader && fRunPar->fsourcetype == 2 && fRunPar->fsimu_pedestalfile.size() > 0 )
            {
                fRawDataReader->initTraceNoiseGenerator( 0, fRunPar->fsimu_pedestalfile, getDetectorGeo(), fRunPar->fsumwindow_1,
                        fDebug, fRunPar->fgrisuseed, fRunPar->fsimu_pedestalfile_DefaultPed, fRunPar->fGainCorrection );
            }
            if( fRawDataReader && fRunPar->finjectGaussianNoise > 0. )
            {
                fRawDataReader->injectGaussianNoise( fRunPar->finjectGaussianNoise, fRunPar->finjectGaussianNoiseSeed );
            }
        }
    }
    // something went wrong, probably wrong filename
    catch( VFileException ex )
    {
        cout << ex.what() << endl;
        //      cout << "data file not found, exiting: " << fRunPar->fsourcefile << endl;
        // !!! no solution, should be something else
        if( !fRunPar->fdisplaymode )
        {
            exit( EXIT_FAILURE );
        }
        else
        {
            return false;
        }
    }
#endif
    // ============================
    // sourcefile has MC grisu format
    if( fRunPar->fsourcetype == 1 )
    {
        if( fGrIsuReader != 0 )
        {
            delete fGrIsuReader;
        }
        fGrIsuReader = new VGrIsuReader( getDetectorGeo(), getDetectorGeo()->getNumTelescopes(), fRunPar->fsourcefile, fRunPar->fsumwindow_1,
                                         fRunPar->ftelescopeNOffset, fRunPar->fsampleoffset, fRunPar->fMCScale,
                                         fDebug, fRunPar->fgrisuseed,
                                         fRunPar->fsimu_pedestalfile, fRunPar->fIgnoreCFGversions );
        fGrIsuReader->setTraceFile( fRunPar->ftracefile );
    }
    
    // ============================
    // source has MC grisu format (multiple files)
    else if( fRunPar->fsourcetype == 5 )
    {
        if( fMultipleGrIsuReader != 0 )
        {
            delete fMultipleGrIsuReader;
        }
        fMultipleGrIsuReader = new VMultipleGrIsuReader( fNTel, fRunPar->fTelToAnalyze, fDebug );
        fMultipleGrIsuReader->init( getDetectorGeo(), fRunPar->fsourcefile, fRunPar->fsumwindow_1,
                                    fRunPar->ftelescopeNOffset, fRunPar->fsampleoffset, fRunPar->fMCScale,
                                    fRunPar->fgrisuseed, fRunPar->fsimu_pedestalfile, true, fRunPar->fsimu_pedestalfile_DefaultPed );
        fMultipleGrIsuReader->setTraceFile( fRunPar->ftracefile );
    }
    // ============================
    // sourcefile has DST format
    else if( fRunPar->fsourcetype == 4 || fRunPar->fsourcetype == 7 )
    {
        if( fDSTReader != 0 )
        {
            delete fDSTReader;
        }
        fDSTReader = new VDSTReader( fRunPar->fsourcefile, fRunPar->fIsMC, fRunPar->fNTelescopes, fDebug );
        if( fDSTReader->isMC() && fRunPar->fIsMC == 0 )
        {
            fRunPar->fIsMC = 1;
        }
        for( unsigned int i = 0; i <  fRunPar->fTelToAnalyze.size(); i++ )
        {
            fDSTReader->setNumSamples( fRunPar->fTelToAnalyze[i], getNSamples( fRunPar->fTelToAnalyze[i] ) );
        }
    }
    // sourcefile has PE format
    else if( fRunPar->fsourcetype == 6 )
    {
        if( fPEReader != 0 )
        {
            delete fPEReader;
        }
        fPEReader = new VPEReader( fRunPar->fsourcefile, fRunPar->fTelToAnalyze, getDetectorGeo(), fDebug );
    }
    // ============================
    // set the data readers for all inherent classes
    initializeDataReader();
    
    // ============================
    // read pixel values from DB
    fDB_PixelDataReader = 0;
    if( fRunPar->useDB() && getDetectorGeo() )
    {
        fDB_PixelDataReader = new VDB_PixelDataReader( getDetectorGeo()->getNumChannelVector() );
        fDB_PixelDataReader->setDebug( fRunPar->fDebug );
        fDB_PixelDataReader->readFromDB( fRunPar->getDBServer(), fRunPar->frunnumber,
                                         fRunPar->fDBRunStartTimeSQL, fRunPar->fDBRunStoppTimeSQL );
    }
    
    // set event number vector
    fTelescopeEventNumber.assign( fNTel, 0 );
    // set event times
    fEventMJD.assign( fNTel, 0 );
    fEventTime.assign( fNTel, 0. );
    
    // set number of channels
    for( unsigned int i = 0; i <  fRunPar->fTelToAnalyze.size(); i++ )
    {
        setTelID( fRunPar->fTelToAnalyze[i] );
        // this is to ignore the photodiode
        if( fReader->getMaxChannels() - getNChannels() != 1 )
        {
            setNChannels( fReader->getMaxChannels() );
        }
        // check telescope configuration
        if( fReader->getMaxChannels() == 0 )
        {
            cout << "VEventLoop::initEventLoop error: telescope " << fRunPar->fTelToAnalyze[i] << " with 0 channels" << endl;
            exit( EXIT_FAILURE );
        }
    }
    // initialize analyzers (output files are created as well here)
    initializeAnalyzers();
    
    
    // create calibrators, analyzers, etc. at first event
    if( fCalibrator )
    {
        fCalibrator->initialize();
    }
    
    // initialize pedestal calculator
    if( fPedestalCalculator && fRunPar->fPedestalsInTimeSlices )
    {
        fPedestalCalculator->initialize( ( fRunMode == R_PED ),  getNChannels(), fRunPar->fPedestalsLengthOfTimeSlice,
                                         fRunPar->fCalibrationSumFirst, fRunPar->fCalibrationSumWindow,
                                         fRunPar->fDBDataStartTimeSecOfDay, fRunPar->fDBDataStoppTimeSecOfDay );
    }
    // print run informations
    printRunInfos();
    
    ////////////////////////////////////////////////////////////////////////////////
    // set array pointing (values valid for all telescope)
    ////////////////////////////////////////////////////////////////////////////////
    cout << endl;
    cout << "----------------------" << endl;
    cout << "Initialize pointing..." << endl;
    fArrayPointing = new VArrayPointing();
    fArrayPointing->setObservatory( fRunPar->getObservatory_Longitude_deg(), fRunPar->getObservatory_Latitude_deg() );
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Monte Carlo file
    if( fRunPar->fIsMC != 0 )
    {
        fArrayPointing->setMC();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    // data
    else
    {
        ///////////////////////////////////////////////////////////////////////////////////////////
        // set target coordinates from command line or from DB
        ///////////////////////////////////////////////////////////////////////////////////////////
        if( fRunPar->fTargetDec > -99. && fRunPar->fTargetRA > -99. )
        {
            fArrayPointing->setTargetName( fRunPar->fTargetName );
            fArrayPointing->setTargetJ2000( fRunPar->fTargetDec, fRunPar->fTargetRA );
        }
    }
    // add any offsets to the pointing [J2000]
    fArrayPointing->setPointingOffset( fRunPar->fTargetRAOffset, fRunPar->fTargetDecOffset );
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // set pointing for all telescopes
    ///////////////////////////////////////////////////////////////////////////////////////////
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( isTeltoAna( i ) )
        {
            fPointing.push_back( new VPointing( i ) );
            fPointing.back()->setObservatory( fRunPar->getObservatory_Longitude_deg(), fRunPar->getObservatory_Latitude_deg() );
            ///////////////////////////////////////////////////////////////////////////////////////////
            // Monte Carlo file
            if( fRunPar->fIsMC != 0 )
            {
                fPointing.back()->setMC();
            }
            ///////////////////////////////////////////////////////////////////////////////////////////
            // data
            else
            {
                ///////////////////////////////////////////////////////////////////////////////////////////
                // set target coordinates from command line or from DB
                ///////////////////////////////////////////////////////////////////////////////////////////
                if( fRunPar->fTargetDec > -99. && fRunPar->fTargetRA > -99. )
                {
                    fPointing.back()->setTargetName( fRunPar->fTargetName );
                    fPointing.back()->setTargetJ2000( fRunPar->fTargetDec, fRunPar->fTargetRA );
                }
            }
            // add any offsets to the pointing [J2000]
            fPointing.back()->setPointingOffset( fRunPar->fTargetRAOffset, fRunPar->fTargetDecOffset );
            // set pointing error
            if( fRunPar->fDBTracking )
            {
                fPointing.back()->getPointingFromDB( fRunPar->frunnumber, fRunPar->fDBTrackingCorrections, fRunPar->fPMTextFileDirectory,
                                                     fRunPar->fDBVPM, fRunPar->fDBUncalibratedVPM );
            }
            else
            {
                fPointing.back()->setPointingError( fRunPar->fPointingErrorX[i], fRunPar->fPointingErrorY[i] );
            }
        }
        else
        {
            fPointing.push_back( 0 );
        }
    }
    // set coordinates in run parameter
    // (J2000)
    if( getArrayPointing() )
    {
        getRunParameter()->fTargetDec = getArrayPointing()->getTargetDecJ2000();
        getRunParameter()->fTargetRA  = getArrayPointing()->getTargetRAJ2000();
    }
    
    return true;
}


void VEventLoop::initializeAnalyzers()
{
    if( fDebug )
    {
        cout << "VEventLoop::initializeAnalyzers()" << endl;
    }
    
    // initialize the image analyzers
    if( fAnalyzer
            && fRunMode != R_PED && fRunMode != R_PEDLOW             // no pedestal analysis
            && fRunMode != R_GTO && fRunMode != R_GTOLOW             // no timing analysis
            && fRunMode != R_TZERO && fRunMode != R_TZEROLOW )       // no toffset analysis
    {
        fAnalyzer->initializeDataReader();
        fAnalyzer->initOutput();
    }
    // initialize the array analyzers
    if( fArrayAnalyzer
            && fRunMode != R_PED && fRunMode != R_PEDLOW             // no pedestal analysis
            && fRunMode != R_GTO && fRunMode != R_GTOLOW             // no timing analysis
            && fRunMode != R_TZERO && fRunMode != R_TZEROLOW )       // no toffset analysis
    {
        fArrayAnalyzer->initializeDataReader();
        fArrayAnalyzer->initOutput();
        fArrayAnalyzer->initTree();
    }
    // initialize the DST class
    if( fDST )
    {
        fDST->initialize();
    }
    
    // initialize the model3D analysis
    if( fRunPar->fUseModel3D && fModel3D )
    {
        if( getOutputFile() )
        {
            getOutputFile()->cd();
            fModel3D->initialize();
            if( fRunPar->fCreateLnLTable )
            {
                fModel3D->createLnLTable();
                fRunPar->fUseModel3D = false;
                exit( EXIT_FAILURE );
            }
        }
        else
        {
            cout << "Error initialzing Model3D trees; no output file available" << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    // set analysis data storage classes
    // (slight inconsistency, produce VImageAnalyzerData for all telescopes,
    //  not only for the requested ones (in teltoana))
    if( fAnaData.size() == 0 )
    {
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            if( i < fAnaDir.size() && fAnaDir[i] )
            {
                fAnaDir[i]->cd();
            }
            setTelID( i );
            fAnaData.push_back( new VImageAnalyzerData( i, fRunPar->fShortTree, ( fRunMode == R_PED || fRunMode == R_PEDLOW ||
                                fRunMode == R_GTO || fRunMode == R_GTOLOW ||
                                fRunMode == R_TZERO || fRunMode == R_TZEROLOW ) ) );
            int iseed = fRunPar->fMCNdeadSeed;
            if( iseed != 0 )
            {
                iseed += i;
            }
            fAnaData.back()->initialize( getNChannels(), getReader()->getMaxChannels(),
                                         getDebugFlag(), iseed, getNSamples(),
                                         getRunParameter()->fpulsetiminglevels.size(), getRunParameter()->fpulsetiming_tzero_index,
                                         getRunParameter()->fpulsetiming_width_index, getRunParameter()->fpulsetiming_triggertime_index );
            if( fRunMode == R_DST )
            {
                fAnaData.back()->initializeMeanPulseHistograms();
                fAnaData.back()->initializeIntegratedChargeHistograms();
            }
            fAnaData.back()->setTraceIntegrationMethod( getRunParameter()->fTraceIntegrationMethod[i] );
        }
        // reading special channels for all requested telescopes
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            if( getTeltoAna()[i] < fAnaData.size() && fAnaData[getTeltoAna()[i]] )
            {
                fAnaData[getTeltoAna()[i]]->readSpecialChannels( getRunNumber(),
                        fRunPar->fsetSpecialChannels,
                        getRunParameter()->getDirectory_EVNDISPParameterFiles() );
            }
        }
        // initialize cleaning
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            setTelID( i );
            if( !getImageCleaningParameter() || !getImageCleaningParameter( true ) )
            {
                cout << "VEventLoop::initializeAnalyzers() error initializing image cleaning for telescope " << getTelID() + 1 << endl;
                exit( EXIT_FAILURE );
            }
        }
    }
}


/*!

     clean up and write all data to disk

*/
void VEventLoop::shutdown()
{
    // additional output for writing to disk (MC only)
    bool fDebug_writing = false;
    if( isMC() )
    {
        fDebug_writing = true;
    }
    if( fDebug )
    {
        cout << "VEventLoop::shutdown()" << endl;
        fDebug_writing = fDebug;
    }
    endOfRunInfo();
    cout << endl << "-----------------------------------------------" << endl;
    
    // if we have the proper settings,
    // print the dead pixel information
    if( ( fRunPar->frunmode == R_ANA || fRunPar->frunmode == R_GTO ) && fRunPar->fprintdeadpixelinfo )  // DEADCHAN
    {
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            cout << endl;
            setTelID( getTeltoAna()[i] );
            printDeadChannels( false , true );
            printDeadChannels( true , true );
        }
    }
    
    // write detector parameter tree to disk
    if( fOutputfile != 0 && fRunPar->foutputfileName != "-1" )
    {
        fOutputfile->cd();
    }
    else if( fRunPar->frunmode == R_DST && fDST && fDST->getDSTFile() )
    {
        fDST->getDSTFile()->cd();
    }
    else if( fRunPar->frunmode != R_PED && fRunPar->frunmode != R_PEDLOW
             && fRunMode != R_GTO && fRunMode != R_GTOLOW
             && fRunMode != R_TZERO && fRunMode != R_TZEROLOW )
    {
        cout << "VEventLoop::shutdown: Error accessing output file" << endl;
    }
    // write run parameter to disk
    if( fRunPar->frunmode != R_PED && fRunPar->frunmode != R_GTO && fRunPar->frunmode != R_GTOLOW
            && fRunPar->frunmode != R_PEDLOW && fRunPar->frunmode != R_TZERO && fRunPar->frunmode != R_TZEROLOW )
    {
        int i_nbytes = fRunPar->Write();
        if( fDebug_writing )
        {
            cout << "WRITEDEBUG: runparameters (nbytes " << i_nbytes << "):";
            if( fOutputfile )
            {
                cout << fOutputfile->Get( "runparameterV2" );
            }
            cout << endl;
        }
    }
    // analysis or trace library mode
    if( fRunPar->frunmode == R_ANA )
    {
        // write information about detector to disk
        if( getDetectorTree() )
        {
            if( fDebug )
            {
                cout << "\t writing detector tree: " << getDetectorTree()->GetName() << endl;
            }
            int i_nbytes = getDetectorTree()->Write();
            if( fDebug_writing )
            {
                cout << "WRITEDEBUG: detector tree (nbytes " << i_nbytes << "):";
                if( fOutputfile )
                {
                    cout << fOutputfile->Get( "telconfig" );
                }
                cout << endl;
            }
        }
        // write pedestal variation calculations to output file
        if( ( fRunPar->frunmode == R_ANA ) && fRunPar->fPedestalsInTimeSlices && fPedestalCalculator )
        {
            fPedestalCalculator->terminate( true, fDebug_writing );
        }
        // calculate and write deadtime calculation to disk
        if( fDeadTime && !isMC() )
        {
            fDeadTime->calculateDeadTime();
            fDeadTime->printDeadTime();
            fDeadTime->writeHistograms( fDebug_writing );
        }
        // write MC run header to output file
        if( getReader()->getMonteCarloHeader() )
        {
            fOutputfile->cd();
            getReader()->getMonteCarloHeader()->print();
            int i_nbytes = getReader()->getMonteCarloHeader()->Write();
            if( fDebug_writing )
            {
                cout << "WRITEDEBUG: MC run header(nbytes " << i_nbytes << "):";
                if( fOutputfile )
                {
                    cout << fOutputfile->Get( "MC_runheader" );
                }
                cout << endl;
            }
        }
        
        // orgainze and write out tree
        if( fRunPar->frunmode == R_ANA && fDeadPixelOrganizer )
        {
            fDeadPixelOrganizer->printSummary() ;
            
            // copy tree to output root file
            fDeadPixelOrganizer->finalize() ;
            
        }
        
        // write array analysis results to output file
        if( fArrayAnalyzer )
        {
            fArrayAnalyzer->terminate( fDebug_writing );
        }
        if( fRunPar->fUseModel3D && fModel3D )
        {
            fModel3D->terminate();
        }
#ifndef NOGSL
        if( fRunPar->ffrogsmode )
        {
            fFrogs->terminate();
        }
#endif
        // write analysis results for each telescope to output file
        if( fAnalyzer )
        {
            for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
            {
                fAnalyzer->setTelID( fRunPar->fTelToAnalyze[i] );
                fAnalyzer->terminate( fDebug_writing );
            }
        }
        // close output file here (!! CLOSE OUTPUT FILE FOREVER !!)
        fAnalyzer->shutdown();
    }
    // write calibration/analysis results for each telescope
    else if( fRunPar->frunmode == R_PED || fRunPar->frunmode == R_GTO || fRunPar->frunmode == R_GTOLOW || fRunPar->frunmode == R_PEDLOW
             || fRunPar->frunmode == R_TZERO || fRunPar->frunmode == R_TZEROLOW )
    {
        VPedestalCalculator* iP = 0;
        if( fRunPar->frunmode == R_PED && fRunPar->fPedestalsInTimeSlices && fPedestalCalculator )
        {
            iP = fPedestalCalculator;
            fPedestalCalculator->terminate( false );
        }
        if( fCalibrator )
        {
            fCalibrator->terminate( iP );
        }
    }
    // write data summary
    else if( fDST && fRunPar->frunmode == R_DST )
    {
        fDST->terminate();
    }
    // delete readers
    if( fRunPar->fsourcetype != 0 && fGrIsuReader )
    {
        delete fGrIsuReader;
    }
    if( fDebug )
    {
        cout << "VEventLoop::shutdown() ... finished" << endl;
    }
    // final check of output file; just open and close it again
    if( fRunMode == R_ANA )
    {
        if( fDebug )
        {
            cout << "VEventLoop::shutdown: final check of output file" << endl;
        }
        TFile f( fRunPar->foutputfileName.c_str() );
        //       f.Recover();
        if( f.TestBit( TFile::kRecovered ) )
        {
            cout << "Warning: output file has been recovered" << endl;
        }
        if( f.IsZombie() )
        {
            cout << "Error: problem with eventdisplay output file: " << fRunPar->foutputfileName << endl;
        }
        else
        {
            cout << endl << "Final checks on result file (seems to be OK): " << fRunPar->foutputfileName << endl;
        }
        // FROGS finishing here
        // (GM) not clear why this has to happen at this point in the program
        // (logically wrong)
#ifndef NOGSL
        if( fRunPar->ffrogsmode )
        {
            fFrogs->finishFrogs( &f );
        }
#endif
        f.Close();
    }
    // end of analysis
    cout << endl;
    cout << "END OF ANALYSIS, exiting..." << endl;
    cout << "===========================" << endl;
    cout << endl;
}


/*!
  \param gEv goto this event number
  (gEv==0 means reset file and goto event number 1 )
*/
void VEventLoop::gotoEvent( int gEv )
{
    if( fDebug )
    {
        cout << "VEventLoop::gotoEvent() " << gEv << endl;
    }
    bool i_res = false;
    // goto event number 0, which means, reset file and look for first event
    if( gEv == 0 )
    {
        // reset file, intialize calibrator and analyzer
        initEventLoop( fRunPar->fsourcefile );
        return;
    }
    // goto event number gEv (backward in sourcefile)
    // event number is smaller than current eventnumber
    else if( gEv - int( fEventNumber ) <= 0 )
    {
        // reset file, start at the beginning and search for this event
        initEventLoop( fRunPar->fsourcefile );
        gotoEvent( gEv );
        return;
    }
    // goto eventnumber gEv (forward in sourcefile)
    else
    {
        fAnalyzeMode = false;
        // go forward in file and search for event gEv
        while( ( int )fEventNumber != gEv )
        {
            i_res = nextEvent();
            if( fReader->getEventStatus() > 998 || !fTimeCutsfNextEventStatus )
            {
                i_res = -1;
                break;
            }
        }
        // event number larger than number of events in file
        if( i_res <= 0 )
        {
            cout <<  "VEventLoop::gotoEvent( int gEv ): event not found: " << gEv << endl;
            return;
        }
        // event found, analyze it
        if( i_res )
        {
            analyzeEvent();
        }
    }
}


/*!
   \param iEvents number of events to be analyzed (iEvents = -1, analyze all events in current file)

   \return is always true
*/
bool VEventLoop::loop( int iEvents )
{
    if( fDebug )
    {
        cout << "VEventLoop::loop()" << endl;
    }
    // print a statement every 5000 events
    int i = 0;
    bool iEventStatus = true;
    fNumberofIncompleteEvents = 0;
    fNumberofGoodEvents = 0;
    
    // Skip to start eventnumber
    if( fRunPar->fFirstEvent > 0 )
    {
        gotoEvent( fRunPar->fFirstEvent );
    }
    
    while( ( i < iEvents || iEvents < 0 ) && iEventStatus )
    {
        iEventStatus = nextEvent();
        if( !iEventStatus )
        {
            if( fReader->getEventStatus() > 998 || fEndCalibrationRunNow )
            {
                break;
            }
            else
            {
                fNumberofIncompleteEvents++;
            }
            iEventStatus = true;
        }
        else
        {
            fNumberofGoodEvents++;
        }
        if( i == 0 )
        {
            cout << endl;
            cout << "##########################################" << endl;
            cout << "########  starting analysis  ############# " << endl;
            cout << "##########################################" << endl;
            cout << endl;
        }
        else if( fRunPar->fPrintAnalysisProgress > 0 && i % fRunPar->fPrintAnalysisProgress == 0 )
        {
            cout << "\t now at event " << i << endl;
        }
        i++;
    }
    terminate( i );
    return true;
}


/*!
  checking event cuts only in analysis mode

  \return
    true if getting the next event was succesful
*/
bool VEventLoop::nextEvent()
{
    if( fDebug )
    {
        cout << "VEventLoop::nextEvent()" << endl;
    }
    int i_Analysis_cut = 1;
    int i_Time_cut = 1;
    do                                            // => while( i_Analysis_cut == 0 );
    {
        // get next event from data reader and check
        // if there is a next event (or EOF) ??
        if( !fReader->getNextEvent() )
        {
            // check if this getNextEvent() failed due to an invalid event
            if( fReader->getEventStatus() < 999 )
            {
                return false;
            }
            else
            {
                cout << "!!! void VEventLoop::nextEvent(): no next event (end of file)" << endl;
                // if the display is run in the loop mode, goto event 0 and start again
                if( fRunPar->floopmode )
                {
                    gotoEvent( 0 );
                    nextEvent();
                    continue;
                }
                else
                {
                    return false;
                }
            }
            return false;
        }
        fReader->setTelescopeID( getTeltoAna()[0] );
        if( !fReader->hasFADCTrace() )
        {
            for( unsigned int i = 0; i < getRunParameter()->fTraceIntegrationMethod.size(); i++ )
            {
                getRunParameter()->fTraceIntegrationMethod[i] = 0;
            }
        }
        // grisu sims only (currently)
        // set FADC hilo mulitplier as read from simulation run header in the vbf file
        // do this only for the first event
        if( getRunParameter()->fsimu_HILO_from_simFile && fReader->getMonteCarloHeader() )
        {
            if( fReader->getMonteCarloHeader()->fFADC_hilo_multipler > 0 )
            {
                getRunParameter()->fsimu_HILO_from_simFile = false;
                for( unsigned int i = 0; i < getNTel(); i++ )
                {
                    setLowGainMultiplier_Trace( i, fReader->getMonteCarloHeader()->fFADC_hilo_multipler );
                    getDetectorGeometry()->setLowGainMultiplier_Trace( i, fReader->getMonteCarloHeader()->fFADC_hilo_multipler );
                }
                cout << "Lowgain multiplier (trace) read from MC run header: " << fReader->getMonteCarloHeader()->fFADC_hilo_multipler << endl;
            }
        }
        fillTriggerVectors();
        ///////////////////////////////////////////////////////////
        // set eventnumbers
        ///////////////////////////////////////////////////////////
        // event numbers for array event
        if( fReader->getArrayTrigger() )
        {
            fEventNumber = int( fReader->getArrayTrigger()->getEventNumber() );
        }
        else if( fReader->isMC() || fReader->isDST() )
        {
            fEventNumber = int( fReader->getEventNumber() );
        }
        else
        {
            fEventNumber = 99999999;
        }
        // event numbers for telescope events
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            fReader->setTelescopeID( getTeltoAna()[i] );
            getTelescopeEventNumber()[getTeltoAna()[i]] = fReader->getEventNumber();
        }
        // set event time from data reader
        setEventTimeFromReader();
        // check time into the run
        i_Time_cut = checkTimeCuts();
        if( i_Time_cut == 2 )
        {
            fTimeCutsfNextEventStatus = 0;
            return false;
        }
        else if( i_Time_cut == 1 )
        {
            continue;
        }
        ///////////////////////////////////////////////////////////
        // in displaymode, look for user interaction
        if( fRunPar->fdisplaymode )
        {
            gSystem->ProcessEvents();
        }
        // analyze event ( do this always except if searching for a specifing event number in the file)
        if( fAnalyzeMode )
        {
            i_Analysis_cut = analyzeEvent();
        }
    }
    while( i_Analysis_cut == 0 && fNextEventStatus );
    fTimeCutsfNextEventStatus = true;
    // user cut failed
    if( i_Analysis_cut < 0 )
    {
        return false;
    }
    return true;
}


/*!
     check run mode and call the analyzers
*/
int VEventLoop::analyzeEvent()
{
    if( fDebug )
    {
        cout << "VEventLoop::analyzeEvent()" << endl;
        cout << "\t now at event " << getEventNumber() << endl;
        cout << "----------------------------------------" << endl;
    }
    // analysis is running
    fAnalyzeMode = true;
    int i_cut = 0;
    int i_cutTemp = 0;
    
    // short cut for dst writing
    if( fRunMode == R_DST && fDST )
    {
#ifndef NOVBF
        if( fReader->getATEventType() != VEventType::PED_TRIGGER )
#endif
        {
            fDST->fill();
            return 1;
        }
    }
    
    ////////////////////////////////////
    // analyze all requested telescopes
    ////////////////////////////////////
    for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
    {
        setTelID( fRunPar->fTelToAnalyze[i] );
        
        if( isMC() && !bMCSetAtmosphericID )
        {
        
            if( fRunPar->fAtmosphereID == 0 && getReader()->getMonteCarloHeader() )
            {
                fRunPar->fAtmosphereID = getReader()->getMonteCarloHeader()->atmosphere;
            }
            bMCSetAtmosphericID = true;
        }
        
        // check number of samples
        if( getTelID() < fBoolPrintSample.size() && fBoolPrintSample[getTelID()] && !isDST_MC() )
        {
            cout << "setting sample length for telescope " << getTelID() + 1 << " to " << getNSamples() << endl;
            fBoolPrintSample[getTelID()] = false;
            // make sure that calibration sum window is not too long
            if( fRunMode == R_PED || fRunMode == R_TZERO )
            {
                if( fRunPar->fCalibrationSumFirst + fRunPar->fCalibrationSumWindow > ( int )getNSamples() )
                {
                    cout << "VEventLoop::analyzeEvent: resetting calibration sum window from ";
                    cout << fRunPar->fCalibrationSumWindow;
                    fRunPar->fCalibrationSumWindow = getNSamples() - fRunPar->fCalibrationSumFirst;
                    cout << " to " << fRunPar->fCalibrationSumWindow;
                    cout << " (sum first at " << fRunPar->fCalibrationSumFirst;
                    cout << ", min sum per channel " << fRunPar->fCalibrationIntSumMin;
                    cout << ")" << endl;
                }
            }
        }
        // quit when number of samples if set to '0'
        if( getNSamples() == 0 && fRunPar->fsourcetype != 7 )
        {
            cout << "VEventLoop::analyzeEvent() error: retrieved sample length of zero" << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        
        // check the requested sumwindow is not larger than the number of samples.
        // Also check that correct low gain multipliers were read in for all 'reset' sum windows.
        if( ( int )getNSamples() < ( int )fRunPar->fsumwindow_1[fRunPar->fTelToAnalyze[i]] )
        {
            if( fBoolSumWindowChangeWarning < 1 && fRunPar->fsourcetype != 7 && fRunPar->fsourcetype != 6 && fRunPar->fsourcetype != 4 )
            {
                cout << "VEventLoop::analyzeEvent: resetting summation window 1 from ";
                cout << fRunPar->fsumwindow_1[fRunPar->fTelToAnalyze[i]] << " to " << getNSamples() << endl;
            }
            fRunPar->fsumwindow_1[fRunPar->fTelToAnalyze[i]] = getNSamples();
            fBoolSumWindowChangeWarning = 1;
        }
        // look for low gain multiplier with nominal window = sumwindow_1. If not found, exit for data analysis.
        //jJust warn for display, calib runs, or if nocalibnoproblem is enabled.
        checkLowGainMultipliers( fRunPar->fTraceIntegrationMethod[fRunPar->fTelToAnalyze[i]], fRunPar->fsumwindow_1[fRunPar->fTelToAnalyze[i]], "sumwindow 1" );
        
        if( ( int )getNSamples() < ( int )fRunPar->fsumwindow_2[fRunPar->fTelToAnalyze[i]] )
        {
            if( fBoolSumWindowChangeWarning < 1 && fRunPar->fsourcetype != 7 && fRunPar->fsourcetype != 6 && fRunPar->fsourcetype != 4 )
            {
                cout << "VEventLoop::analyzeEvent: resetting summation window 2 from ";
                cout << fRunPar->fsumwindow_2[fRunPar->fTelToAnalyze[i]] << " to " << getNSamples() << endl;
            }
            fRunPar->fsumwindow_2[fRunPar->fTelToAnalyze[i]] = getNSamples();
            fBoolSumWindowChangeWarning = 1;
        }
        // look for low gain multiplier with nominal window = sumwindow_2. If not found, exit for data analysis.
        // Just warn for display, calib runs, or if nocalibnoproblem is enabled.
        checkLowGainMultipliers( fRunPar->fTraceIntegrationMethod[fRunPar->fTelToAnalyze[i]], fRunPar->fsumwindow_2[fRunPar->fTelToAnalyze[i]], "sumwindow 2" );
        
        
        if( ( int )getNSamples() < fRunPar->fsumwindow_pass1[fRunPar->fTelToAnalyze[i]] && fRunPar->fDoublePass[fRunPar->fTelToAnalyze[i]] )
        {
            if( fBoolSumWindowChangeWarning < 2 && fRunPar->fsourcetype != 7 && fRunPar->fsourcetype != 6 && fRunPar->fsourcetype != 4 )
            {
                cout << "VEventLoop::analyzeEvent: resetting double pass summation window from ";
                cout << fRunPar->fsumwindow_pass1[fRunPar->fTelToAnalyze[i]] << " to " << getNSamples() << endl;
            }
            fRunPar->fsumwindow_pass1[fRunPar->fTelToAnalyze[i]] = getNSamples();
            fBoolSumWindowChangeWarning = 2;
        }
        
        // look for low gain multiplier with nominal window = sumwindow_pass1. If not found, exit for data analysis.
        // Just warn for display, calib runs, or if nocalibnoproblem is enabled.
        if( fRunPar->fDoublePass[fRunPar->fTelToAnalyze[i]] )
        {
            checkLowGainMultipliers( fRunPar->fTraceIntegrationMethod_pass1[fRunPar->fTelToAnalyze[i]], fRunPar->fsumwindow_pass1[fRunPar->fTelToAnalyze[i]], "sumwindow pass1" );
        }
        
        switch( fRunMode )
        {
            /////////////////
            // analysis
            case R_ANA:                           // analysis mode
                // ignore pedestal events (important for VBF only)
#ifndef NOVBF
                if( fReader->getATEventType() != VEventType::PED_TRIGGER )
#endif
                {
                    if( !fRunPar->fWriteTriggerOnly || fReader->hasArrayTrigger() )
                    {
                        fAnalyzer->doAnalysis();
                    }
                    // check user cuts
                    // (single telescope cuts don't work
                    // what to do:
                    // determine cut result for each telescope
                    // AND then with vector of cut telescopes )
                    i_cutTemp = checkCuts();
                    if( i_cut > 0 && !fCutTelescope )
                    {
                        i_cut = 1;
                    }
                    else
                    {
                        i_cut = i_cutTemp;
                    }
                }
                break;
            /////////////////
            // pedestal calculation
            case R_PED:
                i_cut = 1;
                if( fRunPar->fUsePedEvents )
                {
                    // look for pedestal events
                    // (for VBF only: identify pedestal events by event type)
#ifndef NOVBF
                    if( fReader->getATEventType() == VEventType::PED_TRIGGER )
#endif
                    {
                        fCalibrator->calculatePedestals();
                    }
                }
                else
                {
                    if( !fReader->wasLossyCompressed() )
                    {
                        // for DSTs, require a local trigger
                        if( fReader->hasLocalTrigger( getTelID() ) )
                        {
                            fCalibrator->calculatePedestals();
                        }
                    }
                }
                break;
            /////////////////
            // low gain pedestal calculation
            case R_PEDLOW:
                i_cut = 1;
                if( !fReader->wasLossyCompressed() )
                {
                    fCalibrator->calculatePedestals( true );
                }
                break;
            /////////////////
            // gain and toffset calculation
            case R_GTO:
                i_cut = 1;
                // don't use pedestal events for gain calculation (vbf file only)
#ifndef NOVBF
                if( fReader->getATEventType() != VEventType::PED_TRIGGER )
#endif
                {
                    fCalibrator->calculateGainsAndTOffsets();
                }
                break;
            /////////////////
            // gains/toffset calculation for low gain channels
            case R_GTOLOW:
                i_cut = 1;
                // don't use pedestal events for gain calculation (vbf file only)
#ifndef NOVBF
                if( fReader->getATEventType() != VEventType::PED_TRIGGER )
#endif
                {
                    fCalibrator->calculateGainsAndTOffsets( true );
                }
                break;
            /////////////////
            // mean tzero calculation
            case R_TZERO:
                i_cut = 1;
                // don't use pedestal events for tzero calculation (vbf file only)
#ifndef NOVBF
                if( fReader->getATEventType() != VEventType::PED_TRIGGER &&  fReader->hasArrayTrigger() )
#endif
                {
                    fCalibrator->calculateAverageTZero();
                }
                break;
            /////////////////
            // mean tzero calculation (low gain)
            case R_TZEROLOW:
                i_cut = 1;
                // don't use pedestal events for tzero calculation (vbf file only)
#ifndef NOVBF
                if( fReader->getATEventType() != VEventType::PED_TRIGGER )
#endif
                {
                    fCalibrator->calculateAverageTZero( true );
                }
                break;
            /////////////////
            // this should not happen
            default:
                break;
                
        }
    }
    /////////////////////////////////////////////////////////////////////////
    // ARRAY ANALYSIS
    if( fRunMode != R_PED && fRunMode != R_GTO && fRunMode != R_GTOLOW && fRunMode != R_PEDLOW && fRunMode != R_TZERO && fRunMode != R_TZEROLOW )
    {
#ifndef NOVBF
        if( fReader->getATEventType() != VEventType::PED_TRIGGER )
#endif
        {
            fArrayAnalyzer->doAnalysis();
            // Model3D analysis
            if( fRunPar->fUseModel3D && fReader->hasArrayTrigger() )
            {
                fModel3D->doModel3D();
            }
            // Frogs Analysis
#ifndef NOGSL
            if( fRunPar->ffrogsmode )
            {
                string fArrayEpoch = getRunParameter()->fInstrumentEpoch;
                fFrogs->doFrogsStuff( fEventNumber, fArrayEpoch );
            }
#endif
        }
    }
    
    /////////////////////////////////////////////////////////////////////////
    // dead time calculation
    if( !isMC() && fDeadTime )
    {
        if( fReader->getArrayTrigger() && fReader->getArrayTrigger()->hasTenMHzClockArray() )
        {
            fDeadTime->fillDeadTime( getEventTime(), fReader->getArrayTrigger()->getTenMHzClockArray() );
        }
        else
        {
            fDeadTime->fillDeadTime( getEventTime(), 0 );
        }
    }
    
    /////////////////////////////////////////////////////////////////////////
    // calculate pedestals in time slices
    // (this should be done at this point, as MJD time from array timing is used)
    if( fRunMode == R_ANA || fRunMode == R_PED )
    {
        for( unsigned int i = 0; i < fRunPar->fTelToAnalyze.size(); i++ )
        {
            // for pedestal calculation: require pedestal event (vbf file only)
#ifndef NOVBF
            if( fRunMode == R_PED && fRunPar->fUsePedEvents && fReader->getATEventType() != VEventType::PED_TRIGGER )
            {
                continue;
            }
#endif
            setTelID( fRunPar->fTelToAnalyze[i] );
            if( fRunPar->fPedestalsInTimeSlices && !fReader->wasLossyCompressed() )
            {
                fPedestalCalculator->doAnalysis( fRunMode == R_PEDLOW );
            }
        }
    }
    // these cuts are important for display mode only
    if( fRunMode != R_PED && fRunMode != R_GTO && fRunMode != R_GTOLOW && fRunMode != R_PEDLOW && fRunMode != R_TZERO && fRunMode != R_TZEROLOW )
    {
        i_cut = int( checkArrayCuts() == 1 && i_cut > 0 );
    }
    
    // number of event in calibration runs
    if( fRunMode == R_PED || fRunMode == R_PEDLOW || fRunMode == R_TZERO || fRunMode == R_TZEROLOW )
    {
        if( fRunPar->fNCalibrationEvents > 0 )
        {
            if( ( int )fCalibrator->getNumberOfEventsUsedInCalibration( -1, fRunMode ) > fRunPar->fNCalibrationEvents )
            {
                i_cut = -1;
                fEndCalibrationRunNow = true;
            }
        }
    }
    
    // add dead pixel states to VDeadPixelOrganizer
    if( fDeadPixelOrganizer )
    {
        // get this event's info
        int    eventMJD    = fArrayEventMJD  ;
        double eventTime   = fArrayEventTime ;
        int    eventNumber = fEventNumber    ;
        
        // set up some initial variables
        bool   higGain   = false ;
        PixelStateInt lowGainState = 121212 ;
        
        
        // loop over all telescopes
        // itel is from 0-3
        for( unsigned int itel = 0 ; itel < getTeltoAna().size() ; itel++ )
        {
            setTelID( getTeltoAna()[itel] ) ;
            
            // loop over all? lowgain pixels
            // ipix is from 0-498
            //cout << coutprefix << "getDead(highGain).size:" << getDead(higGain).size() << endl;
            for( unsigned int ipix = 0 ; ipix < getDead( higGain ).size() ; ipix++ )
            {
                lowGainState  = ( PixelStateInt ) getDead( higGain )[ipix] ;
                
                fDeadPixelOrganizer->UpdatePixelState( itel + 1, ipix + 1, higGain,  eventMJD, eventTime,  lowGainState ) ;
            }
            
        } // endfor: no more telescopes
        
        fDeadPixelOrganizer->updatePreviousEventInfo( eventNumber, eventMJD, eventTime ) ;
        
    } // endif: fDeadPixelOrganizer doesn't exist
    
    
    //}
    
    return i_cut;
}


int VEventLoop::checkArrayCuts()
{
    if( fDebug )
    {
        cout << "VEventLoop::checkArrayCuts()" << endl;
    }
    // donnot apply array cuts when analysing one telescopes only
    if( getTeltoAna().size() < 2 )
    {
        return 1;
    }
    // array cuts
    if( ( int )getShowerParameters()->fNTrig < fNCutNArrayTrigger )
    {
        return 0;
    }
    if( ( int )getShowerParameters()->fShowerNumImages[0] < fNCutArrayImages )
    {
        return 0;
    }
    
    return 1;
}


/*!
   using mechanism of ROOT trees in method TTree::Draw( .., char *selection, ..)
   (see http://root.cern.ch/root/html400/TTree.html#TTree:Draw)

   example: plot only events with alpha<10 and more than 10 tubes
   write into the cut box: alpha<10&&ntubes>10.
*/
int VEventLoop::checkCuts()
{
    if( fDebug )
    {
        cout << "VEventLoop::checkCuts()" << endl;
    }
    if( !fChangedCut )
    {
        return 1;    // no cuts are applied
    }
    
    int i_cut;
    // number of triggered channels
    int i_numtrig = 0;
    i_numtrig = getImageParameters()->ntubes;
    if( i_numtrig < 0 || i_numtrig < fNCutNTrigger[getTelID()] )
    {
        return 0;
    }
    
    // very ugly, but no better idea... (as well no better idea from ROOT people)
    // define temporarly a root tree and use its selection mechanism
    if( fStringCut[getTelID()].length() > 0 )
    {
        TTree i_tree( "i_tree", "" );
        float cen_x = 0., cen_y = 0., length = 0., width = 0., size = 0., azwidth = 0.;
        float dist = 0., asymmetry = 0.;
        float muonSize = 0., muonRadius = 0., muonRSigma = 0., muonIPCorrectedSize = 0.;   // muon (Iterative fit muon analysis)
        double houghAP = 0, houghTD = 0., houghCN = 0., houghContained = 0.; // muon (Hough)
        int houghMuonValid = 0, houghNpix = 0; //Hough
        short int fLocalTrigger = 0;
        float MCenergy = 0.;
        unsigned short int ntubes = 0, bad = 0, nlowgain = 0, nsat = 0;
        int muonValid = 0;
        unsigned short int eventType = 0;
        int fitStat = 0;
        float loss = 0.;
        
        i_tree.Branch( "eventType", &eventType, "eventType/s" );
        i_tree.Branch( "cen_x", &cen_x, "cen_x/F" );
        i_tree.Branch( "cen_y", &cen_y, "cen_y/F" );
        i_tree.Branch( "length", &length, "length/F" );
        i_tree.Branch( "size", &size, "size/F" );
        i_tree.Branch( "width", &width, "width/F" );
        i_tree.Branch( "dist", &dist, "dist/F" );
        i_tree.Branch( "asymmetry", &asymmetry, "asymmetry/F" );
        i_tree.Branch( "azwidth", &azwidth, "azwidth/F" );
        i_tree.Branch( "loss", &loss, "loss/F" );
        i_tree.Branch( "fitStat", &fitStat, "fitStat/I" );
        i_tree.Branch( "nlowgain", &nlowgain, "nlowgain/s" );
        i_tree.Branch( "nsat", &nsat, "nsat/s" );
        i_tree.Branch( "ntubes", &ntubes, "ntubes/s" );
        i_tree.Branch( "bad", &bad, "bad/s" );
        i_tree.Branch( "MCenergy", &MCenergy, "MCenergy/F" );
        i_tree.Branch( "fLocalTrigger", &fLocalTrigger, "fLocalTrigger/S" );
        if( i_tree.GetBranchStatus( "muonRadius" ) )
        {
            i_tree.Branch( "muonRadius", &muonRadius, "muonRadius/F" );
        }
        if( i_tree.GetBranchStatus( "muonRSigma" ) )
        {
            i_tree.Branch( "muonRSigma", &muonRSigma, "muonRSigma/F" );
        }
        if( i_tree.GetBranchStatus( "muonSize" ) )
        {
            i_tree.Branch( "muonSize", &muonSize, "muonSize/F" );
        }
        if( i_tree.GetBranchStatus( "muonIPCorrectedSize" ) )
        {
            i_tree.Branch( "muonIPCorrectedSize", &muonIPCorrectedSize, "muonIPCorrectedSize/F" );
        }
        if( i_tree.GetBranchStatus( "muonValid" ) )
        {
            i_tree.Branch( "muonValid", &muonValid, "muonValid/I" );
        }
        if( i_tree.GetBranchStatus( "houghMuonValid" ) )
        {
            i_tree.Branch( "houghMuonValid", &houghMuonValid, "houghMuonValid/I" );
        }
        if( i_tree.GetBranchStatus( "houghAP" ) )
        {
            i_tree.Branch( "houghAP", &houghAP, "houghAP/D" );
        }
        if( i_tree.GetBranchStatus( "houghTD" ) )
        {
            i_tree.Branch( "houghTD", &houghTD, "houghTD/D" );
        }
        if( i_tree.GetBranchStatus( "houghNpix" ) )
        {
            i_tree.Branch( "houghNpix", &houghNpix, "houghNpix/I" );
        }
        if( i_tree.GetBranchStatus( "houghCN" ) )
        {
            i_tree.Branch( "houghCN", &houghCN, "houghCN/D" );
        }
        if( i_tree.GetBranchStatus( "houghContained" ) )
        {
            i_tree.Branch( "houghContained", &houghContained, "houghContained/D" );
        }
        
        eventType = fAnalyzer->getImageParameters()->eventType;
        cen_x = fAnalyzer->getImageParameters()->cen_x;
        cen_y = fAnalyzer->getImageParameters()->cen_y;
        length = fAnalyzer->getImageParameters()->length;
        width = fAnalyzer->getImageParameters()->width;
        size = fAnalyzer->getImageParameters()->size;
        dist = fAnalyzer->getImageParameters()->dist;
        asymmetry = fAnalyzer->getImageParameters()->asymmetry;
        azwidth = fAnalyzer->getImageParameters()->azwidth;
        loss = fAnalyzer->getImageParameters()->loss;
        fitStat = fAnalyzer->getImageParameters()->Fitstat;
        ntubes = fAnalyzer->getImageParameters()->ntubes;
        nsat = fAnalyzer->getImageParameters()->nsat;
        nlowgain = fAnalyzer->getImageParameters()->nlowgain;
        bad = fAnalyzer->getImageParameters()->bad;
        MCenergy = fAnalyzer->getImageParameters()->MCenergy;
        fLocalTrigger = fAnalyzer->getImageParameters()->fLocalTrigger;
        muonSize = fAnalyzer->getImageParameters()->muonSize;
        muonIPCorrectedSize = fAnalyzer->getImageParameters()->muonIPCorrectedSize;
        muonRadius = fAnalyzer->getImageParameters()->muonRadius;
        muonRSigma = fAnalyzer->getImageParameters()->muonRSigma;
        muonValid = fAnalyzer->getImageParameters()->muonValid;
        houghMuonValid = fAnalyzer->getImageParameters()->houghMuonValid;
        houghAP = fAnalyzer->getImageParameters()->houghAP;
        houghTD = fAnalyzer->getImageParameters()->houghTD;
        houghNpix = fAnalyzer->getImageParameters()->houghNpix;
        houghCN = fAnalyzer->getImageParameters()->houghCN;
        houghContained = fAnalyzer->getImageParameters()->houghContained;
        
        i_tree.Fill();
        i_cut = int( i_tree.Draw( "cen_x", fStringCut[getTelID()].c_str(), "goff" ) );
        return i_cut;
    }
    // end off ugly stuff
    return 1;
}


/*!
     reset cut strings and variables
 */
void VEventLoop::resetRunOptions()
{
    if( fDebug )
    {
        cout << "VEventLoop::resetRunOptions()" << endl;
    }
    fCutTelescope = false;
    fChangedCut = false;
    fNCutNArrayTrigger = 0;
    fNCutArrayImages = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        fNCutNTrigger.push_back( 0 );
        fStringCut.push_back( "" );
    }
}


/*!
      will be implemented as soon as all features of the vbf reader is used
*/
void VEventLoop::previousEvent()
{
    cout << "not implemented, fFileRead->getPrevEvent crashes always" << endl;
}


/*!
    analyze only events with more than iNtrigger triggered PMTs

   \param iNtrigger more than iNtrigger tubes triggered in this event
*/
void VEventLoop::setCutNTrigger( int iNtrigger )
{
    if( fDebug )
    {
        cout << "VEventLoop::setCutNTrigger()" << endl;
    }
    fNCutNTrigger[fTelID] = iNtrigger;
    fChangedCut = true;
}


void VEventLoop::setCutNArrayTrigger( int iNtrigger )
{
    if( fDebug )
    {
        cout << "VEventLoop::setCutNArrayTrigger()" << endl;
    }
    fNCutNArrayTrigger = iNtrigger;
    fChangedCut = true;
}


void VEventLoop::setCutNArrayImages( int iNImages )
{
    if( fDebug )
    {
        cout << "VEventLoop::setCutNArrayImages()" << endl;
    }
    fNCutArrayImages = iNImages;
    fChangedCut = true;
}


/*!
    set the user cut string

   \param iCut string with ROOT tree cut (like in TTree::draw()), used in checkCuts()
*/
void VEventLoop::setCutString( string iCut )
{
    if( fDebug )
    {
        cout << "VEventLoop::setCutString()" << endl;
    }
    fStringCut[getTelID()] = iCut;
    fChangedCut = true;
}


void VEventLoop::fillTriggerVectors()
{
    // first check how many telescope participate
    unsigned int iNTrig = 0;
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( getReader()->hasLocalTrigger( getTeltoAna()[i] ) )
        {
            iNTrig++;
        }
    }
    // now fill the corresponding fields per telescope
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( getReader()->hasLocalTrigger( getTeltoAna()[i] ) && iNTrig < getTriggeredTel()[getTeltoAna()[i]].size() )
        {
            getTriggeredTel()[getTeltoAna()[i]][iNTrig]++;
        }
    }
    // fill total number of multiplicities
    if( iNTrig < getTriggeredTelN().size() )
    {
        getTriggeredTelN()[iNTrig]++;
    }
}


void VEventLoop::terminate( int iAna )
{
    cout << endl;
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
    cout << "End of event loop, finishing up..." << endl;
    cout << endl;
    if( fReader->getNIncompleteEvents().size() > 0 && fReader->getNIncompleteEvents().size() >= getTeltoAna().size() )
    {
        cout << "Number of incomplete events:" << endl;
        for( unsigned int t = 0; t < getTeltoAna().size(); t++ )
        {
            cout << "\t Telescope " << getTeltoAna()[t] + 1 << ": " << fReader->getNIncompleteEvents()[getTeltoAna()[t]] << endl;
        }
    }
    if( !isMC() )
    {
        cout << "Number of events with GPS faults (status bit set): " << endl;
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            cout << "\t Telescope " << getTeltoAna()[i] + 1 << ": ";
            if( i < fGPSClockWarnings.size() )
            {
                cout << fGPSClockWarnings[i];
            }
            cout << endl;
        }
    }
    
    cout << endl << "Analyzed " << iAna << " events" << endl;
}

/*

    setting event times from reader

    take care of non-working GPS clocks (apply majority rules)

*/
void VEventLoop::setEventTimeFromReader()
{
    /////////////////////////////////////////////////////////////////////////
    // ignore event time calculation for DSTs
    if( fReader->getDataFormatNum() == 4 || fReader->getDataFormatNum() == 6 )
    {
        return;
    }
    
    /////////////////////////////////////////////////////////////////////////
    // event times setting for VBF sources
#ifndef NOVBF
    /////////////////////////////////////////////////////////////////////////
    unsigned int iCurrentTelID = getTelID();
    
    VGPSDecoder fGPS;
    ///////////////// /////////////////////////////////// /////////////////
    // test if all times are the same, apply majority rule otherwise
    ///////////////// /////////////////////////////////// /////////////////
    // allow differences between event times of this value
    double i_max_time_diff = 1.e-4;
    // check number of telescopes
    if( getNTel() > VDST_MAXTELESCOPES )
    {
        cout << " VEventLoop::setEventTimeFromReader: warning, cannot apply majority rule to times, too many telescopes: ";
        cout << VDST_MAXTELESCOPES << " " << getNTel() << endl;
        return;
    }
    map< unsigned int, double > i_telescope_time;
    map< unsigned int, int > i_telescope_timeN;
    map< unsigned int, double > i_MJD;
    
    // get times for each telescope
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        getReader()->setTelescopeID( getTeltoAna()[i] );
        
        fGPS.decode( fReader->getGPS0(), fReader->getGPS1(), fReader->getGPS2(), fReader->getGPS3(), fReader->getGPS4() );
        // check status of GPS clock
        if( fGPS.getStatus() != 0 )
        {
            if( i < fGPSClockWarnings.size() )
            {
                fGPSClockWarnings[i]++;
            }
            if( fGPSClockWarnings[i] < 30 && fDebug )
            {
            
                cout << " VEventLoop::setEventTimeFromReader: info, event with GPS error status in telescope " << getTeltoAna()[i] + 1 << endl;
                cout << "\t error status " << fGPS.getStatus();
                if( fGPS.getStatus() & 0x0001 )
                {
                    cout << " (unlocked from GPS time source)" << endl;
                }
                else if( fGPS.getStatus() & 0x0002 )
                {
                    cout << " (offset in time with respect to UTC)" << endl;
                }
                else if( fGPS.getStatus() & 0x0004 )
                {
                    cout << " (large frequency offset)" << endl;
                }
            }
            else if( fGPSClockWarnings[i] == 30 )
            {
                cout << " VEventLoop::setEventTimeFromReader: info, more than 30 events with GPS status set";
                cout << " (Telescope " << getTeltoAna()[i] + 1 << ")" << endl;
            }
        }
        
        // time is given in seconds per day
        if( getTelID() < fEventTime.size() )
        {
            fEventTime[getTelID()] = fGPS.getHrs() * 60.*60. + fGPS.getMins() * 60. + fGPS.getSecs();
            fArrayEventTime = fEventTime[getTelID()];
        }
        i_telescope_time[getTeltoAna()[i]] = fGPS.getHrs() * 60.*60. + fGPS.getMins() * 60. + fGPS.getSecs();
        i_telescope_timeN[getTeltoAna()[i]] = 0;
        //! fGPS.getDays() gives day in year, so I calculate the MJD for
        //! 1st January of this year, then add fGPS.getDays()-1.
        int  j = 0;
        double dMJD = 0.;
        VAstronometry::vlaCldj( fReader->getATGPSYear() + 2000, 1, 1, &dMJD, &j );
        dMJD += fGPS.getDays() - 1.;
        if( fReader->isGrisuMC() && dMJD == 51543 )
        {
            dMJD = 54383;
        }
        // horrible fudge to deal with broken T3 clock (for MJD=54101 only)
        if( getTelID() == 2 && dMJD == 54101 )
        {
            dMJD += 100;
        }
        i_MJD[getTeltoAna()[i]] = dMJD;
        // set MJD
        if( getTelID() < fEventMJD.size() )
        {
            fEventMJD[getTelID()] = ( int )dMJD;
            fArrayEventMJD = ( int )dMJD;
        }
    }
    
    
    ///// Time in [s] of the day /////
    // count equal times
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        for( unsigned int j = 0; j < getTeltoAna().size(); j++ )
        {
            if( i == j )
            {
                continue;
            }
            
            if( fabs( i_telescope_time[getTeltoAna()[i]] - i_telescope_time[getTeltoAna()[j]] ) < i_max_time_diff )
            {
                i_telescope_timeN[getTeltoAna()[i]]++;
            }
        }
    }
    // get time with most majority votes
    int z_max = getTeltoAna()[0];
    int i_nold = 0;
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( i_telescope_timeN[getTeltoAna()[i]] > i_nold )
        {
            z_max = getTeltoAna()[i];
            i_nold = i_telescope_timeN[getTeltoAna()[i]];
        }
    }
    if( fabs( i_telescope_time[z_max] - fEventTime[getTelID()] ) > i_max_time_diff )
    {
        fEventTime[getTelID()] = i_telescope_time[z_max];
        fArrayEventTime = i_telescope_time[z_max];
    }
    //// MJD ////
    // check if all MJDs are the same (use same routines as for time)
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        i_telescope_timeN[getTeltoAna()[i]] = 0;
    }
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        for( unsigned int j = 0; j < getTeltoAna().size(); j++ )
        {
            if( i == j )
            {
                continue;
            }
            
            if( fabs( i_MJD[getTeltoAna()[i]] - i_MJD[getTeltoAna()[j]] ) < i_max_time_diff )
            {
                i_telescope_timeN[getTeltoAna()[i]]++;
            }
        }
    }
    // get MJD with majority rule
    z_max = getTeltoAna()[0];
    i_nold = 0;
    for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
    {
        if( i_telescope_timeN[getTeltoAna()[i]] > i_nold )
        {
            z_max = getTeltoAna()[i];
            i_nold = i_telescope_timeN[getTeltoAna()[i]];
        }
    }
    if( fabs( i_MJD[z_max] - fEventMJD[getTelID()] ) > i_max_time_diff )
    {
        fEventMJD[getTelID()] = ( int )i_MJD[z_max];
        fArrayEventMJD = fEventMJD[getTelID()];
    }
    // check if MJD of current event is different from the value of the previous event
    // this is only ok if we are close to midnight
    if( !isMC() && fArrayPreviousEventMJD > 0 && fArrayPreviousEventMJD != fArrayEventMJD && fArrayEventTime < 86400. - 30. )
    {
        cout << "VEventLoop::setEventTimeFromReader: warning,";
        cout << " sudden jump in MJD between previous and current event";
        cout << " (Telescope " << getTelID() + 1 << ", eventnumber " << getEventNumber() << "): " << endl;
        cout << "\t current event: MJD " << fArrayEventMJD << ", Time " << fArrayEventTime << endl;
        cout << "\t previous event: MJD " << fArrayPreviousEventMJD << endl;
        cout << "\t using MJD of previous event" << endl;
        fEventMJD[getTelID()] = fArrayPreviousEventMJD;
        fArrayEventMJD = fArrayPreviousEventMJD;
        cout << "\t GPS clock status: " << fGPS.getStatus() << endl;
    }
    fArrayPreviousEventMJD = fArrayEventMJD;
    /////////////////////////////////////////////////////////////////////
    // end of time fixes
    /////////////////////////////////////////////////////////////////////
    setTelID( iCurrentTelID );
#endif
}

int VEventLoop::checkTimeCuts()
{
    if( fTimeCut_RunStartSeconds == 0 )
    {
        fTimeCut_RunStartSeconds = fArrayEventTime;
    }
    
    if( getRunParameter()->fTimeCutsMin_min > 0 && ( fArrayEventTime - fTimeCut_RunStartSeconds ) < getRunParameter()->fTimeCutsMin_min * 60 )
    {
        return 1;
    }
    if( getRunParameter()->fTimeCutsMin_max > 0 && ( fArrayEventTime - fTimeCut_RunStartSeconds ) > getRunParameter()->fTimeCutsMin_max * 60 )
    {
        return 2;
    }
    
    return 0;
}


// check that low gain multipliers for a given integration method/integration window have been read in.
// for data analysis, exit if low gain multipliers are not found.
// for display, -nocalibnoproblem, or MC files, just print a warning.

void VEventLoop::checkLowGainMultipliers( unsigned int iMethod, int iSumWindow, TString iName )
{
    if( fLowGainMultiplierWarning < 3 && ( fRunPar->frunmode == 0 || fRunPar->frunmode == 3 || fRunPar->frunmode == 4 ) )
    {
        bool found = false;
        for( unsigned int j = 0; j < getLowGainDefaultSumWindows().size() ; j++ )
        {
            if( getLowGainDefaultSumWindows()[j].first  == iMethod && getLowGainDefaultSumWindows()[j].second == iSumWindow )
            {
                found = true;
                break;
            }
        }
        if( !found )
        {
            if( fRunPar->frunmode == 0 && !fRunPar->fdisplaymode && !fRunPar->fNoCalibNoPb && !fRunPar->fIsMC )
            {
                cout << "VEventLoop::analyzeEvent error: No low gain multipliers available for " << iName << " (";
                cout << iSumWindow << " samples); trace integration method " << iMethod << ". Exiting. " << endl;
                cout << "Use option -nocalibnoproblem to force analysis." << endl;
                exit( EXIT_FAILURE );
            }
            fLowGainMultiplierWarning++;
            cout << "VEventLoop::analyzeEvent warning: No low gain multipliers available for " << iName << " (";
            cout << iSumWindow << " samples); trace integration method " << iMethod << endl;
            if( fRunPar->getObservatory() == "VERITAS" )
            {
                cout << "\tThis is expected for Grisu MC, but bad if you are analysing data or Care MC." << endl;
            }
        }
    }
}


