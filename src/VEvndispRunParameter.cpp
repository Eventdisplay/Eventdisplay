/*! \class VEvndispRunParameter
  \brief input parameter storage class

*/

#include "VEvndispRunParameter.h"


VEvndispRunParameter::VEvndispRunParameter( bool bSetGlobalParameter ) : VGlobalRunParameter( bSetGlobalParameter )
{
    SetName( "runparameterV2" );
    SetTitle( getEVNDISP_VERSION().c_str() );
    fEventDisplayUser = "";
    fEventDisplayHost = "";
    fEventDisplayDate = "";
    // system parameters
    fEventDisplayBuildCompiler = "";
    fEventDisplayBuildCompilerVersion = "";
    fEventDisplayBuildArch = "";
    fEventDisplayBuildNode = "";
    fEventDisplayBuildDir = "";
    fEventDisplayBuildROOTVersion = "";
    fEventDisplayBuildROOTVersionInt = 0;
    fEventDisplaySystemInfo = 0;

    fSGE_TASK_ID = 0;

    // debug parameters
    fDebug = false;
    fPrintSmallArray = true;

    // run parameters
#ifdef RUNWITHDB
    fuseDB = true;
#else
    fuseDB = false;
#endif
    frunmode = 0;
    fRunIsZeroSuppressed = false;
    frunnumber = -1;
    fRunTitle  = "";
    fsourcetype = 3;           // 0 = rawdata, 1 = GrIsu simulation, 2 = MC in VBF format,
    // 3 = rawdata in VBF, 4 = DST (data), 5 = multiple GrIsu file,
    // 6 = PE file, 7 = DST (MC)
    fsourcefile = "";

    fDBRunType = "";
    fDBDataStartTimeMJD = 0.;
    fDBDataStoppTimeMJD = 0.;
    fDBDataStartTimeSecOfDay = 0.;
    fDBDataStoppTimeSecOfDay = 0.;
    fDBRunStartTimeSQL = "";
    fDBRunStoppTimeSQL = "";
    fsimu_pedestalfile = "";
    fsimu_HILO_from_simFile = false;
    fsimu_noiselevel   = 250;
    fsimu_pedestalfile_DefaultPed = 20.;
    fsimu_lowgain_pedestal_DefaultPed = -999.;
    fCombineChannelsForPedestalCalculation = 0;
    fPedestalSingleRootFile = false;
    fnevents = -10000;
    fFirstEvent = -10000;
    fTimeCutsMin_min = -99;
    fTimeCutsMin_max = -99;
    fIsMC = 0;
    fIgnoreCFGversions = false;
    fPrintAnalysisProgress = 25000;
    fRunDuration = 1. * 3600.;        // default run duration is 1 h (reset by DBRunInfo)
    fPrintGrisuHeader = 0;
    finjectGaussianNoise = -1.;
    finjectGaussianNoiseSeed = 0;

    fprintdeadpixelinfo = false ; // DEADCHAN if true, print list of dead pixels to evndisp.log

    fSaveDeadPixelRegistry = false;

    // geometry/calibration parameters
    fNTelescopes = 4;                             // there is always at least one telescope
    fcamera.push_back( "EVN_V4_Autumn2007_20130110.txt" );
    fCalibrationDataType = 1;  // should be 0 for e.g. CTA DSTs
    fcalibrationfile = "";
    fLowGainCalibrationFile = "calibrationlist.LowGain.dat";
    fcalibrationrun = false;
    fNCalibrationEvents = -1;
    fLaserSumMin = 50000.;
    fGainFileNumber.push_back( 0 );
    fTOffFileNumber.push_back( 0 );
    fPedFileNumber.push_back( 0 );
    fTZeroFileNumber.push_back( 0 );
    fLowGainMultiplierFileNumber.push_back( 0 );
    fPedLowGainFileNumber.push_back( 0 );
    fPedLowGainFile = "";
    fGainLowGainFileNumber.push_back( 0 );
    fTOffLowGainFileNumber.push_back( 0 );
    fTZeroLowGainFileNumber.push_back( 0 );
    fPixFileNumber.push_back( 0 );
    fIgnoreDSTGains = false;
    faverageTZeroFiducialRadius = 1.5;

    fTelToAnalyze.push_back( 0 );

    fDeadChannelFile = "EVNDISP.validchannels.dat";

    fEpochFile = "VERITAS.Epochs.runparameter";
    fInstrumentEpoch = "noepoch";
    fAtmosphereID = 0;

    fCameraCoordinateTransformX = 1.;
    fCameraCoordinateTransformY = 1.;
    //calibration read from DB
    freadCalibfromDB = false;
    freadCalibfromDB_versionquery = -1000;
    freadCalibfromDB_save_file = false;
    fNoCalibNoPb = false;
    fNextDayGainHack = false;

    // pointing parameters
    fTargetName = "NONAME";
    fTargetDec = -999.;
    fTargetRA = -999.;
    fTargetDecOffset = 0;
    fTargetRAOffset = 0;
    felevation = -999.;
    fazimuth = -999.;
    fWobbleNorth = 0.;
    fWobbleEast = 0.;
    fCheckPointing = 999.;
#ifdef RUNWITHDB
    fDBTracking = true;
    fDBCameraRotationMeasurements = true;
    fDBVPM = true;
    fDBUncalibratedVPM = false;
#else
    fDBTracking = false;
    fDBCameraRotationMeasurements = false;
    fDBVPM = false;
    fDBUncalibratedVPM = false;
#endif
    fPointingErrorX.push_back( 0. );
    fPointingErrorY.push_back( 0. );
    // star catalogue
    fStarCatalogueName = "";
    fMinStarBrightness_B = 7.;
    fMinStarPixelDistance_deg = -1.;
    fMinStarNTubes = 100000;

    fGainCorrection.push_back( 1. );

    fLowGainPeds = true;

    // analyzer parameters
    fImageCleaningParameters.push_back( new VImageCleaningRunParameter( "cleaningParameters_mainPass" ) );
    fImageCleaningParameters_DB_Pass1.push_back( new VImageCleaningRunParameter( "cleaningParameters_DP_Pass1" ) );

    fsumfirst.push_back( 2 );
    fSearchWindowLast.push_back( 9999 );
    fsumwindow_1.push_back( 12 );
    fsumwindow_2.push_back( 12 );
    fsumwindow_pass1.push_back( 18 );
    fNeighbourDistanceFactor.push_back( 1. );
    fSquarePixels.push_back( false );
    fNeighbourMultiplicity.push_back( 0 );
    fImageLL = 0;
    fLogLikelihoodLoss_min.push_back( 1.e3 );
    fLogLikelihoodLoss_max.push_back( 1.e3 );
    fLogLikelihood_Ntubes_min.push_back( 0 );
    fForceLLImageFit = false;
    fMinimizeTimeGradient = false;
    fMinimizeTimeGradient_minGradforFit = 2.5;
    fMinimizeTimeGradient_minLoss = 0.1;
    fMinimizeTimeGradient_minNtubes = 15;
    fSquaredImageCalculation = false;
    fImageAnalysisFUIFactor = 2.;
    fFixWindowStart_sumwindow2 = false;
    fDoublePass.push_back( true );
    fDoublePassErrorWeighting2005.push_back( true );
    frecoverImagePixelNearDeadPixel = true;
    fFillImageBorderNeighbours = true;
    fTraceWindowShift.push_back( -1. );
    fsumfirst_startingMethod.push_back( 1 );
    fsumfirst_maxT0startDiff.push_back( 0 );
    fSumWindow_searchmaxreverse.push_back( true );
    fTraceIntegrationMethod.push_back( 1 );
    fTraceIntegrationMethod_pass1.push_back( 1 );
    fDF_DigitalFilter.push_back( 0 );
    fDF_UpSample.push_back( 4 );
    fDF_PoleZero.push_back( 0.75 );
    fSumWindowMaxTimedifferenceToDoublePassPosition.push_back( -4. );
    fSumWindowMaxTimeDifferenceLGtoHG.push_back( -5. );
    fSmoothDead = false;
    fUsePedEvents = true;
    fFADCChargeUnit = "DC";
    // pedestal calculation in time slices
    fLowGainUsePedestalsInTimeSlices = false;
    fUsePedestalsInTimeSlices = true;
    fPedestalsInTimeSlices = true;
    fPedestalsLengthOfTimeSlice = 180.;            //!< [s]
    fCalibrationSumWindow = 16;
    fCalibrationSumFirst = 0;
    fCalibrationSumWindowAverageTime = 6;
    fCalibrationIntSumMin = 20.;
    fL2TimeCorrect = true;
    fsetSpecialChannels = "EVNDISP.specialchannels.dat";
    fthroughputCorrectionFile = "";
    ftraceamplitudecorrectionFile = "";
    freconstructionparameterfile = "EVNDISP.reconstruction.runparameter";

    ////////////////////////////////////////////////////////////////////////////////
    // pulse timing (fraction of maximum where times are determined)
    // OBSERVE: there should be a timing level with value 0.5 (for tzero calculations)
    // OBSERVE: this vector is always symmetric around the entry at 1
    //    fpulsetiminglevels.push_back( 0.2 );
    fpulsetiminglevels.push_back( 0.5 );
    //    fpulsetiminglevels.push_back( 0.8 );
    fpulsetiminglevels.push_back( 1.0 );
    unsigned int i_fps = fpulsetiminglevels.size();
    if( i_fps > 1 )
    {
        for( unsigned int i = 0; i < i_fps - 1; i++ )
        {
            fpulsetiminglevels.push_back( fpulsetiminglevels[i_fps - i - 2] );
        }
    }
    // get index for tzero and width
    setPulseZeroIndex();

    fWriteTriggerOnly = true;
    fShortTree = 1;
    fwriteMCtree = 0;
    fFillMCHistos = true;

    // muon parameters
    fmuonmode = false;
    fhoughmuonmode = false;

    // output parameters
    ffillhistos = false;                          // obsolete
    foutputfileName = "";
    fWriteExtraCalibTree = false;
    fWriteImagePixelList = false;
    // MC parameters
    // offset in telescope numbering (0 for old grisudet version (<3.0.0))
    ftelescopeNOffset = 1;
    fsampleoffset = 0;
    fUseVBFSampleLength = false;
    fMC_FADCTraceStart = 0;
    fgrisuseed = 0;
    ftracefile = "";
    fMCnoDead = false;
    fMCNdead = 0;
    fMCNdeadSeed = 0;
    fMCNdeadboard = 0;
    fMCScale = 1.;

    // display parameters
    fdisplaymode = false;
    floopmode = false;
    fh = 680;                                     // height of main window
    fw = 1000;                                    // width of main window
    fh = 650;                                     // height of main window
    fw = 975;                                     // width of main window
    fPlotRaw = false;
    fPlotPaper = false;
    fPlotAllInOneMethod = 0;

    // writing of laser pulses
    fwriteLaserPulseN = 0;
    fwriteAverageLaserPulse = false;

    // dst parameters
    fdstfile = "";
    fdstminntubes = -1;
    fdstwriteallpixel = true;
    fdstcalibration   = false;

    // NN cleaning parameters
    ifWriteGraphsToFile = false;
    ifReadIPRfromDatabase = false;
    ifCreateIPRdatabase = false;
    ifReadIPRfromDSTFile = false;
    for( unsigned int i = 0; i < VDST_MAXTELTYPES; i++ )
    {
        fFADCtoPhe[i] = -1.;
    }

    fNNGraphsFile = "";
    fIPRdatabase = "";
    fIPRdatabaseFile = "";

    // revisit default parameters depending on the observatory
    setDefaultParameters( getObservatory() );
}

VEvndispRunParameter::~VEvndispRunParameter()
{
    for( unsigned int i = 0; i < fImageCleaningParameters.size(); i++ )
    {
        if( fImageCleaningParameters[i] )
        {
            delete fImageCleaningParameters[i];
        }
    }


}

/*
 * this is the only function where observatory (e.g. VTS or CTA) dependent parameter settings should appear
*/
void VEvndispRunParameter::setDefaultParameters( string iObservatory )
{
    if( iObservatory.find( "CTA" ) != string::npos )
    {
        fsetSpecialChannels = "nofile";
        faverageTZeroFiducialRadius = 0.;
        //      fIgnoreDSTGains = true;
    }
}



void VEvndispRunParameter::print()
{
    print( 1 );
}


/*!
   \param iEv = 0: print important parameters only
              = 1: print in eventdisplay run
              = 2: print all parameters
*/
void VEvndispRunParameter::print( int iEv )
{
    if( fDebug )
    {
        cout << "VEvndispRunParameter::printParams()" << endl;
    }

    // print less to screen for some variables for large arrays
    if( fTelToAnalyze.size() >= 10 )
    {
        fPrintSmallArray = false;
    }

    cout << endl;
    if( iEv == 1 )
    {
        cout << "\t ---------------------------------------------" << endl;
        if( frunmode == 0 && !fdisplaymode )
        {
            cout << "\t       ANALYZING DATA" << endl;
        }
        else if( frunmode == 0 && fdisplaymode )
        {
            cout << "\t       DISPLAYING DATA" << endl;
        }
        else if( frunmode == 1 )
        {
            cout << "\t       CALCULATING PEDESTALS" << endl;
        }
        else if( frunmode == 2 )
        {
            cout << "\t       CALCULATING GAINS AND TIME OFFSETS (high gain channels)" << endl;
        }
        else if( frunmode == 3 )
        {
            cout << "\t       WRITING TRACE LIBRARY" << endl;
        }
        else if( frunmode == 4 )
        {
            cout << "\t       WRITING DATA SUMMARY FILES" << endl;
        }
        else if( frunmode == 5 )
        {
            cout << "\t       CALCULATING GAINS AND TIME OFFSETS (low gain channels)" << endl;
        }
        else if( frunmode == 6 )
        {
            cout << "\t       CALCULATING PEDESTALS (low gain channels)" << endl;
        }
        else if( frunmode == 7 )
        {
            cout << "\t       CALCULATING TZEROS (high gain channels)" << endl;
        }
        else if( frunmode == 8 )
        {
            cout << "\t       CALCULATING TZEROS (low gain channels)" << endl;
        }
        cout << "\t ---------------------------------------------" << endl << endl;
    }
    if( iEv == 2 )
    {
        cout << "Eventdisplay version: " << getEVNDISP_VERSION() << endl;
        cout << "============================" << endl << endl;
    }

    cout << "RUN " << frunnumber;
    if( fRunTitle.size() > 0 )
    {
        cout << " (" << fRunTitle << ")";
    }
    cout << endl;
    cout << "Observatory: " << getObservatory();
    cout << " (lat " << getObservatory_Latitude_deg() << ", long " << getObservatory_Longitude_deg();
    cout << ", height " << getObservatory_Height_m() << "m)";
    cout << endl;
    cout << "File: " << fsourcefile << " (sourcetype " << fsourcetype;
    cout << ")" << endl;
    cout << "===========" << endl;
    cout << fEventDisplayDate;
    if( fIsMC )
    {
        cout << ", MC";
    }
    if( fEventDisplayUser.size() > 0 )
    {
        cout << ", User: " << fEventDisplayUser << " ";
    }
    cout << "Host: " << fEventDisplayHost << endl;
    if( fEventDisplayBuildNode.size() > 0 )
    {
        cout << fEventDisplayBuildNode << endl;
    }
    if( fEventDisplayBuildArch.size() > 0 )
    {
        cout << fEventDisplayBuildArch << "\t";
    }
    if( fEventDisplayBuildCompiler.size() > 0 )
    {
        cout << fEventDisplayBuildCompiler << "\t";
    }
    if( fEventDisplayBuildCompilerVersion.size() > 0 )
    {
        cout << fEventDisplayBuildCompilerVersion << endl;
    }
    if( fEventDisplayBuildDir.size() > 0 )
    {
        cout << fEventDisplayBuildDir << endl;
    }
    if( fEventDisplaySystemInfo )
    {
        cout << fEventDisplaySystemInfo->fModel << ", " << fEventDisplaySystemInfo->fOS;
        cout << ", " << fEventDisplaySystemInfo->fCpuType << ", " << fEventDisplaySystemInfo->fCpuSpeed << " MHz";
        cout << ", " << fEventDisplaySystemInfo->fPhysRam << " MB" << endl;
    }
    if( fEventDisplayBuildROOTVersion.size() > 0 )
    {
        cout << "ROOT version " << fEventDisplayBuildROOTVersion << endl;
    }
    if( fSGE_TASK_ID > 0 )
    {
        cout << "SGE TASK ID " << fSGE_TASK_ID << endl;
    }

    cout << endl;
    if( fTargetName.size() > 0 )
    {
        cout << "Target: " << fTargetName;
    }
    if( fTargetDec > -99 )
    {
        cout << "\t Target: ( J2000 dec=" << fTargetDec << ", ra=" << fTargetRA << ")" << endl;
    }
    cout << "\t offsets (ra,dec): " <<  fTargetRAOffset << ", " << fTargetDecOffset << endl;
    cout << "\t wobble (north,east): " << fWobbleNorth << ", " << fWobbleEast << endl;
    if( fTelToAnalyze.size() < 20 )
    {
        cout << "\t pointing corrections (x,y): ";
        if( !fDBTracking )
        {
            for( unsigned int i = 0; i < fTelToAnalyze.size(); i++ )
            {
                cout << "\t T" << fTelToAnalyze[i] + 1 << ": " << fPointingErrorX[fTelToAnalyze[i]] << ", " << fPointingErrorY[fTelToAnalyze[i]];
            }
            cout << endl;
        }
        else
        {
            cout << " use database" << endl;
        }
    }
    if( fDBRunType.size() > 0 )
    {
        cout << "Run type: " << fDBRunType << endl;
    }
    cout << "Instrument epoch: " << fInstrumentEpoch << "  Atmosphere (corsika ID): " << fAtmosphereID << endl;
    if( fEpochFile.size() > 0 )
    {
        cout << "(epochs read from " << fEpochFile << ")" << endl;
    }
    if( fDBCameraRotationMeasurements )
    {
        cout << "using camera rotation values from DB" << endl;
    }
    cout << endl;
    cout << "analyzing following telescope: ";
    for( unsigned int i = 0; i < fTelToAnalyze.size(); i++ )
    {
        cout << fTelToAnalyze[i] + 1;
        if( i != fTelToAnalyze.size() - 1 )
        {
            cout << ", ";
        }
    }
    cout << endl;
    cout << "detector configuration file: " << fcamera[0];
    if( fIgnoreCFGversions )
    {
        cout << " (ignoring cfg version numbering)";
    }
    cout << endl;
    if( fUseVBFSampleLength )
    {
        cout << "\t using number of FADC samples from cfg file" << endl;
    }
    cout << endl;
    cout << "runmode: " << frunmode << endl;
    if( fnevents > 0 )
    {
        cout << "number of events to analyse: " << fnevents << endl;
    }
    if( fFirstEvent > 0 )
    {
        cout << "starting analysis at event:  " << fFirstEvent << endl;
    }
    if( fTimeCutsMin_min > 0 )
    {
        cout << "start analysing at minute " << fTimeCutsMin_min << endl;
    }
    if( fTimeCutsMin_max > 0 )
    {
        cout << "stop analysing at minute " << fTimeCutsMin_max << endl;
    }
    if( fNCalibrationEvents > 0 )
    {
        cout << "number of events in calibration analysis: " << fNCalibrationEvents << endl;
    }
    if( frunmode == 4 )
    {
        cout << "dstfile: " << fdstfile << " (mintubes: " << fdstminntubes << ")" << endl;
    }
    cout << endl;
    if( fcalibrationfile.size() > 0 )
    {
        cout << "calibration file: " << fcalibrationfile << endl;
    }
    if( fLowGainCalibrationFile.size() > 0 )
    {
        cout << "calibration file (low gain): " << fLowGainCalibrationFile << endl;
    }
    else if( frunmode != 2 && frunmode != 5 && !fIsMC )
    {
        cout << "reading laser/flasher run numbers from database" << endl;
    }
    if( frunmode == 2 )
    {
        cout << "Minimum size required for laser events (lasermin): " << fLaserSumMin << " [dc]" << endl;
    }
    if( fIgnoreDSTGains )
    {
        cout << "Ignoring gains from DST file" << endl;
    }
    if( frunmode == 1 )
    {
        if( ( fsourcetype == 1 || fsourcetype == 2 || fsourcetype == 5 ) && fsimu_pedestalfile.size() > 0 )
        {
            cout << "calculate pedestals from " << fsimu_pedestalfile;
            cout << " with noise level " << fsimu_noiselevel;
            cout << " (default ped: " << fsimu_pedestalfile_DefaultPed;
            cout << ")" << endl;
        }
        else if( fsourcetype == 1 )
        {
            cout << "calculate pedestals from " << fsourcefile << endl;
        }
        if( fPedestalsInTimeSlices )
        {
            cout << "calculating time dependent pedestals" << endl;
        }
        else
        {
            cout << "no time dependent pedestals" << endl;
        }
        if( fUsePedEvents )
        {
            cout << "using pedestal events for pedestal calculation" << endl;
        }
    }
    if( finjectGaussianNoise > 0. )
    {
        cout << "Injecting Gaussian noise with standard deviation " << finjectGaussianNoise;
        cout << " (seed " << finjectGaussianNoiseSeed << ")";
        cout << endl;
    }
    if( fsimu_HILO_from_simFile )
    {
        cout << "reading hilo multiplier from MC run header" << endl;
    }
    if( frunmode == 6 )
    {
        cout << "using low gain events only for pedestal calculation" << endl;
    }
    cout << endl;

    if( fCalibrationDataType == 0 )
    {
        cout << "no calibration data available or calibration data is read from DST file" << endl;
    }
    cout << "signal charge unit is " << fFADCChargeUnit << endl;

    if( frunmode == 0 || frunmode == 4 )
    {
        if( fUsePedestalsInTimeSlices || fLowGainUsePedestalsInTimeSlices )
        {
            cout << "pedestals in time slices: ";
        }
        else
        {
            cout << "pedestals constant over run" << endl;
        }
        if( fUsePedestalsInTimeSlices )
        {
            cout << "high gain";
        }
        if( fLowGainUsePedestalsInTimeSlices )
        {
            cout << " low gain";
        }
        cout << endl;
        cout << "setting special channels (e.g. with L2 signal): " << fsetSpecialChannels << endl;
        if( fthroughputCorrectionFile.size() > 0 )
        {
            cout << "setting throughput correction from file: " << fthroughputCorrectionFile << endl;
        }
        if( ftraceamplitudecorrectionFile.size() )
        {
            cout << "setting throughput correction from file (FADC): ";
            cout << ftraceamplitudecorrectionFile << endl;
        }
        cout << "pulse timing levels: ";
        for( unsigned int i = 0; i < fpulsetiminglevels.size(); i++ )
        {
            cout << fpulsetiminglevels[i] << ", ";
        }
        cout << " (tzero index: " << fpulsetiming_tzero_index;
        cout << ", max index: " << fpulsetiming_max_index;
        if( fpulsetiming_triggertime_index != 9999 )
        {
            cout << ", trigger time index: " << fpulsetiming_triggertime_index;
        }
        cout << endl;
        if( fL2TimeCorrect )
        {
            cout << "Correcting FADC times for create jitter with L2 signals" << endl;
        }
        else
        {
            cout << "No correcting FADC times for create jitter with L2 signals" << endl;
        }
        if( fFixWindowStart_sumwindow2 )
        {
            cout << "using fixed window start for summation window 2" << endl;
        }
        if( fMC_FADCTraceStart > 0 )
        {
            cout << "MC trace start: " << fMC_FADCTraceStart  << endl;
        }
        if( fSmoothDead )
        {
            cout << "smoothing dead pixels" << endl;
        }
        if( fmuonmode )
        {
            cout << "muon ring analysis: " << fmuonmode << endl;
        }
        if( fhoughmuonmode )
        {
            cout << "Hough transform muon ring analysis: " << fhoughmuonmode << endl;
        }

        if( fImageLL != 0 )
        {
            cout << "loglikelihood fitting of images: " << fImageLL;
            cout << " (using these images for the array reconstruction)";
            cout << endl;
        }
        cout << "Fraction of image/border pixel under image ellipse fact (FUI-factor): " << fImageAnalysisFUIFactor << endl;
        if( fSquaredImageCalculation )
        {
            cout << "Use squared weighting for image calculation" << endl;
        }
        // time gradient fitting
        if( fMinimizeTimeGradient )
        {
            cout << "LL with time gradient fitting (min time grad: ";
            cout << fMinimizeTimeGradient_minGradforFit << " deg/s, min loss: ";
            cout << fMinimizeTimeGradient_minLoss << " min ntubes: ";
            cout << fMinimizeTimeGradient_minNtubes << ")" << endl;
        }
    }
    if( ftracefile.size() > 0 )
    {
        cout << "\t tracelib file: " << ftracefile << endl;
    }
    if( fsimu_lowgain_pedestal_DefaultPed > 0. )
    {
        cout << "Low gain ped assumed: " << fsimu_lowgain_pedestal_DefaultPed;
    }
    if( fsourcetype == 1 || fsourcetype == 5 )
    {
        cout << "telescope numbering offset: " << ftelescopeNOffset << endl;
    }
    if( fMCNdead > 0 )
    {
        cout << "Random dead channels: " << fMCNdead << " (seed " <<  fMCNdeadSeed << "), " << fMCNdeadboard << endl;
    }
    if( fPlotPaper )
    {
        cout << " (paper plotting mode)";
    }
    cout << endl;

    cout << "directories:" << endl;
    cout << "\t analysis data: " << getDirectory_EVNDISPAnaData() << endl;
    if( fsourcetype == 0 || fsourcetype == 2 || fsourcetype == 3 )
    {
        cout << "\t raw data: " << getDirectory_VBFRawData() << endl;
    }
    cout << "\t output: " << getDirectory_EVNDISPOutput() << endl;
    if( fShortTree )
    {
        cout << endl << "shortened tree output " << endl;
    }
    if( fwriteMCtree )
    {
        cout << "writing full MC tree with all MC events " << endl;
    }

    // print analysis parameters
    if( iEv == 2 )
    {
        cout << endl;
        for( unsigned int i = 0; i < fTelToAnalyze.size(); i++ )
        {
            cout << "Telescope " << fTelToAnalyze[i] + 1 << endl;

            // trace integration method
            if( fTraceIntegrationMethod[fTelToAnalyze[i]] != 0 )
            {
                cout << "\t trace integration method: \t" << fTraceIntegrationMethod[fTelToAnalyze[i]];
                if( fTelToAnalyze[i] < fDoublePass.size() && fDoublePass[fTelToAnalyze[i]] )
                {
                    cout << "\t (doublepass, integration method pass 1: " << fTraceIntegrationMethod_pass1[fTelToAnalyze[i]] << ")";
                }
                cout << endl;
                if( fDF_DigitalFilter[fTelToAnalyze[i]] != 0 )
                {
                    cout << "\t digital filter: upsample: " << fDF_UpSample[fTelToAnalyze[i]];
                    cout << ", pole-zero cancellation parameter: " << fDF_PoleZero[fTelToAnalyze[i]];
                    cout << endl;
                }
                cout << "\t start of summation window: \t" << fsumfirst[fTelToAnalyze[i]];
                cout << "\t (shifted by " << fTraceWindowShift[fTelToAnalyze[i]] << " samples";
                cout << " [method-" << fsumfirst_startingMethod[fTelToAnalyze[i]] << "], ";
                cout << "max diff to T0 in int start: " << fsumfirst_maxT0startDiff[fTelToAnalyze[i]] << ")";
                if( fSearchWindowLast[fTelToAnalyze[i]] < 9999 )
                {
                    cout << ", last sample for sliding window search: " << fSearchWindowLast[fTelToAnalyze[i]];
                }
                cout << endl;
                cout << "\t length of summation window: \t" << fsumwindow_1[fTelToAnalyze[i]];
                cout << "/" << fsumwindow_2[fTelToAnalyze[i]];
                if( fTelToAnalyze[i] < fDoublePass.size() && fDoublePass[fTelToAnalyze[i]] )
                {
                    cout << "\t length of first pass summation window (double pass): \t" << fsumwindow_pass1[fTelToAnalyze[i]];
                }
                if( fTelToAnalyze[i] < fSumWindow_searchmaxreverse.size() && !fSumWindow_searchmaxreverse[fTelToAnalyze[i]] )
                {
                    cout << " (no reverse tmax search for low gain)";
                }
                cout << endl;
            }
            else
            {
                cout << "\t no trace integration ";
                cout << "[pulse timing method " << fsumfirst_startingMethod[fTelToAnalyze[i]] << "]" << endl;
            }
            // image cleaning method and values
            if( i < fDoublePass.size() && fDoublePass[i] && i < fImageCleaningParameters_DB_Pass1.size() )
            {
                fImageCleaningParameters_DB_Pass1[i]->print();
            }
            if( i < fImageCleaningParameters.size() )
            {
                fImageCleaningParameters[i]->print();
            }
            cout << "\t LL edge fit: \t\t\tloss > " << fLogLikelihoodLoss_min[i];
            cout << "\t loss < " << fLogLikelihoodLoss_max[i];
            cout << "\t ntubes > " << fLogLikelihood_Ntubes_min[i] << endl;

            // calibration
            if( fTelToAnalyze[i] < fGainCorrection.size() && TMath::Abs( fGainCorrection[fTelToAnalyze[i]] ) - 1. > 1.e-2 )
            {
                cout << "\t additional gain correction: " << fGainCorrection[fTelToAnalyze[i]];
            }
            cout << "\t pedestal file: " << fPedFileNumber[i];
            if( i < fPedLowGainFileNumber.size() && fPedLowGainFileNumber[i] > 0 )
            {
                cout << ", low gain pedestal file: " << fPedLowGainFileNumber[i];
            }
            if( i < fGainFileNumber.size() )
            {
                cout << ", gain file: " << fGainFileNumber[i];
            }
            if( i < fGainLowGainFileNumber.size() && fGainLowGainFileNumber[i] > 0 )
            {
                cout << ", low gain gain file: " << fGainLowGainFileNumber[i];
            }
            if( i < fGainCorrection.size() && TMath::Abs( fGainCorrection[i] - 1. ) > 0.001 )
            {
                cout << " (gain correction: " << fGainCorrection[i]  << ")";
            }
            if( i < fTOffFileNumber.size() )
            {
                cout << ", toff file: " << fTOffFileNumber[i];
            }
            if( i < fTOffLowGainFileNumber.size() && fTOffLowGainFileNumber[i] > 0 )
            {
                cout << ", low gain toff file: " << fTOffLowGainFileNumber[i];
            }
            if( i < fPixFileNumber.size() )
            {
                cout << ", pixel file: " << fPixFileNumber[i];
            }
            if( i < fTZeroFileNumber.size() )
            {
                cout << ", tzero file: " << fTZeroFileNumber[i];
            }
            if( i < fTZeroLowGainFileNumber.size() && fTZeroLowGainFileNumber[i] > 0 )
            {
                cout << ", low gain tzero file: " << fTZeroLowGainFileNumber[i];
            }
            cout << endl;
        }
    }

}

void VEvndispRunParameter::setPulseZeroIndex()
{
    fpulsetiming_tzero_index = 9999;
    fpulsetiming_width_index = 9999;
    fpulsetiming_max_index   = 9999;
    // trigger time index is fixed
    fpulsetiming_triggertime_index = fpulsetiminglevels.size();
    for( unsigned int i = 0; i < fpulsetiminglevels.size(); i++ )
    {
        if( TMath::Abs( fpulsetiminglevels[i] - 1. ) < 1.e-4 && fpulsetiming_max_index == 9999 )
        {
            fpulsetiming_max_index = i;
            break;
        }
    }
    for( unsigned int i = 0; i < fpulsetiminglevels.size(); i++ )
    {
        if( TMath::Abs( fpulsetiminglevels[i] - 0.5 ) < 1.e-4 && fpulsetiming_tzero_index == 9999 )
        {
            fpulsetiming_tzero_index = i;
        }
        else if( TMath::Abs( fpulsetiminglevels[i] - 0.5 ) < 1.e-4 && fpulsetiming_tzero_index != 9999 )
        {
            fpulsetiming_width_index = i;
            break;
        }
    }
}

void VEvndispRunParameter::printCTA_DST()
{
    cout << "Eventdisplay version: " << getEVNDISP_VERSION() << endl;
    cout << "============================" << endl << endl;
    cout << fEventDisplayDate << endl;

    cout << "Observatory: " << getObservatory() << endl;

    cout << "RUN " << frunnumber << endl;
    cout << endl;
    cout << "source file " << fsourcefile << endl;
    cout << "number of telescope: " << fNTelescopes << endl;
    cout << "pulse timing levels: ";
    for( unsigned int i = 0; i < fpulsetiminglevels.size(); i++ )
    {
        cout << fpulsetiminglevels[i] << ", ";
    }
    cout << " (tzero index: " << fpulsetiming_tzero_index;
    cout << ", max index: " << fpulsetiming_max_index;
    if( fpulsetiming_triggertime_index != 9999 )
    {
        cout << ", trigger time index: " << fpulsetiming_max_index;
    }
    cout << ")" << endl;
    cout << "L2 timing correct: " << fL2TimeCorrect << endl;
    cout << endl;
}

void VEvndispRunParameter::setSystemParameters()
{
    // get date
    TDatime t_time;
    fEventDisplayDate = t_time.AsSQLString();

    // get host name
    fEventDisplayHost = gSystem->HostName();;
    // get user name
    UserGroup_t* i_userGroup = gSystem->GetUserInfo();
    if( i_userGroup )
    {
        fEventDisplayUser = i_userGroup->fUser;
        delete i_userGroup;
    }
    // get system parameters
    fEventDisplayBuildCompiler = gSystem->GetBuildCompiler();
    fEventDisplayBuildCompilerVersion = gSystem->GetBuildCompilerVersion();
    fEventDisplayBuildArch = gSystem->GetBuildArch();
    fEventDisplayBuildNode = gSystem->GetBuildNode();
    fEventDisplayBuildDir = gSystem->GetBuildDir();
    gSystem->GetSysInfo( fEventDisplaySystemInfo );
    // get root info
    fEventDisplayBuildROOTVersion = gROOT->GetVersion();
    fEventDisplayBuildROOTVersionInt = gROOT->GetVersionInt();

    const char* i_sge = gSystem->Getenv( "SGE_TASK_ID" );
    if( i_sge )
    {
        fSGE_TASK_ID = atoi( i_sge );
    }
}

/*
 * return instrument epoch
 *
 * for VTS, expect this to be in the format
 * MAJOR_MINOR epoch, e.g. V6_2016
 *
 */
string VEvndispRunParameter::getInstrumentEpoch( bool iMajor, bool iUpdateInstrumentEpoch )
{
    // re-read instrument epoch
    if( iUpdateInstrumentEpoch )
    {
        updateInstrumentEpochFromFile();
    }
    if( iMajor )
    {
        return fInstrumentEpoch.substr( 0, fInstrumentEpoch.find( "_" ) );
    }
    return fInstrumentEpoch;
}


/*
   read instrument epoch from file
   (typically VERITAS.Epochs.runparameter)
*/
bool VEvndispRunParameter::updateInstrumentEpochFromFile( string iEpocheFile,
        bool iReadInstrumentEpoch )
{
    if( iEpocheFile != "usedefault" )
    {
        fEpochFile = iEpocheFile;
    }
    if( fEpochFile.size() == 0 )
    {
        return true;
    }

    ifstream is;
    is.open( fEpochFile.c_str(), ifstream::in );
    if( !is )
    {
        string iTemp = getDirectory_EVNDISPParameterFiles() + fEpochFile;
        is.open( iTemp.c_str(), ifstream::in );
        if( !is )
        {
            cout << "error opening epoch parameter file " << fEpochFile << endl;
            cout << iTemp << endl;
            exit( EXIT_FAILURE );
        }
    }
    string is_line;
    string temp;
    string itemp_epoch = "not_found";
    int itemp_atmo = 0;
    int run_min = 0;
    int run_max = 0;
    while( getline( is, is_line ) )
    {
        if( is_line.size() == 0 )
        {
            continue;
        }
        istringstream is_stream( is_line );
        is_stream >> temp >> temp;
        if( iReadInstrumentEpoch && temp == "EPOCH" )
        {
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> itemp_epoch;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> run_min;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> run_max;
            }
            if( frunnumber >= run_min && frunnumber <= run_max )
            {
                break;
            }
        }
        else if( temp == "ATMOSPHERE" )
        {
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> itemp_atmo;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> temp;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> temp;
            }
            int mjd_min = 0;
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> mjd_min;
            }
            int mjd_max = 0;
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> mjd_max;
            }
            if( fDBDataStartTimeMJD >= mjd_min && fDBDataStoppTimeMJD <= mjd_max )
            {
                break;
            }
        }
    }
    if( iReadInstrumentEpoch )
    {
        fInstrumentEpoch = itemp_epoch;
    }
    else
    {
        fAtmosphereID = itemp_atmo;
    }
    is.close();
    return true;
}

unsigned int VEvndispRunParameter::getAtmosphereID( bool iUpdateInstrumentEpoch )
{
    if( iUpdateInstrumentEpoch )
    {
        updateInstrumentEpochFromFile( "usedefault", false );
    }

    return fAtmosphereID;
}
