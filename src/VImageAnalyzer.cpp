/*! \class VImageAnalyzer
    \brief class for analyzing VERITAS data

    image analysis (per telescope)

*/

#include "VImageAnalyzer.h"

VImageAnalyzer::VImageAnalyzer()
{
    fDebug = getDebugFlag();
    if( fDebug )
    {
        cout << "VImageAnalyzer::VImageAnalyzer()" << endl;
    }
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fCalibrated.push_back( false );
    }
    fRaw = false;
    fOutputfile = 0;
    fInit = false;
    
    // image cleaning
    fVImageCleaning = new VImageCleaning( getData() );
    
    // image parameterisation (Hillas parameters)
    fVImageParameterCalculation = new VImageParameterCalculation( getRunParameter()->fShortTree, getData() );
    fVImageParameterCalculation->setDebug( fDebug );
    fVImageParameterCalculation->setDetectorGeometry( getDetectorGeometry() );
    fVImageParameterCalculation->setDetectorGeometry( getDetectorGeometry() );
    
    // initalize minuit for log likelihood method
    if( fRunPar->fImageLL )
    {
        fVImageParameterCalculation->initMinuit( fRunPar->fImageLL );
    }
    else
    {
        fVImageParameterCalculation->initMinuit( 1 );
    }
    
    // initialize root directories
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fAnaDir.push_back( 0 );
    }
    // initialize tgraphs (used for double pass method)
    // (save also graph from first pass)
    char hname[200];
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        sprintf( hname, "xgraphDP1_%u", i + 1 );
        fXGraph.push_back( new TGraphErrors( 1 ) );
        fXGraph.back()->SetName( hname );
        sprintf( hname, "xgraphDP2_%u", i + 1 );
        fXGraphDP2.push_back( new TGraphErrors( 1 ) );
        fXGraphDP2.back()->SetName( hname );
    }
    
    //If -hough is specified on the command line, run the hough transform initialization method in
    //VImageParameterCalculation
    if( fRunPar->fhoughmuonmode || fRunPar->fmuonmode )
    {
    
#ifndef NOGSL
        cout << "Using GSL libraries for muon analysis." << endl;
#else
        cout << "Warning! No GSL libraries found. Muon impact parameter corrected Size will not be calculated." << endl;
#endif
        
        cout << "" << endl;
        if( fRunPar->fhoughmuonmode )
        {
            fVImageParameterCalculation->houghInitialization();
        }
        
    }
    
}


VImageAnalyzer::~VImageAnalyzer()
{
    delete fVImageCleaning;
    delete fVImageParameterCalculation;
}


/*
 *  run through different steps of image cleaning and image parameterization
 *
 *  this is the main loop
 *
 */
void VImageAnalyzer::doAnalysis()
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::doAnalysis() for telescope " << getTelID() + 1 << endl;
    }
    setDebugLevel( 0 );
    
    if( getTelID() < getAnalysisTelescopeEventStatus().size() )
    {
        getAnalysisTelescopeEventStatus()[getTelID()] = 0;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // first call -> initialize and read calibration data
    //            -> find dead channels
    //            -> book trees
    // reset debug level
    if( fReader->getMaxChannels() > 0 && !getImageParameters()->getTree() )
    {
        initAnalysis();
        // set special channel (L2 signal into each of the FADC crates)
        setSpecialChannels();
        // find dead channels, fill tree for this set only (but: dead channels are time dependent)
        findDeadChans( false, true );
        findDeadChans( true, true );
        fInit = true;
    }
    // find dead channels (time dependence)
    else
    {
        findDeadChans( false, false );
        findDeadChans( true, false );
    }
    initEvent();
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // don't do analysis if init event failed or if there was no array trigger
    bool bTrigger = getReader()->hasArrayTrigger();
    // for DSTs, require additionally a local trigger
    if( getRunParameter()->fsourcetype == 7 )
    {
        bTrigger = getReader()->hasArrayTrigger() && getReader()->hasLocalTrigger( getTelID() );
    }
    if( !bTrigger )
    {
        fillOutputTree();
        setSums( 0. );
        setPulseTiming( 0., true );
        setPulseTiming( 0., false );
        setImage( false );
        setBorder( false );
        setImageBorderNeighbour( false );
        setHiLo( false );
        setZeroSuppressed( 0 );
        return;
    }
    if( getDebugFlag() )
    {
        cout << "VImageAnalyzer:doAnalysis array trigger: " << getReader()->hasArrayTrigger() << endl;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // integrate pulses and calculate timing parameters
    if( isDoublePass() )
    {
        calcTZerosSums( getSumFirst(), getSumFirst() + getSumWindow_Pass1(), getTraceIntegrationMethod_pass1() );
    }
    // no double pass (e.g. NN cleaning)
    else
    {
        calcTZerosSums( getSumFirst(), getSumFirst() + getSumWindow(), getTraceIntegrationMethod() );
    }
    
    // number of saturated channels
    getImageParameters()->nsat = fillSaturatedChannels();
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // correct for FADC crate jitter timing
    if( fRunPar->fL2TimeCorrect )
    {
        FADCStopCorrect();
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // apply timing correction from laser calibration (flatfielding in time)
    timingCorrect();
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // image cleaning & gain correction
    imageCleaning( isDoublePass() );
    // print image and border pixels from double pass 1
    if( isDoublePass() && getDebugFlag() )
    {
        for( unsigned int i = 0; i < getImage().size(); i++ )
        {
            if( getImage()[i] || getBorder()[i] )
            {
                cout << "DOUBLE PASS1 (tel " << getTelID() + 1 << "): pixel " << i;
                cout << " is a ";
                if( getImage()[i] )
                {
                    cout << "is a image pixel" << endl;
                }
                else
                {
                    cout << "is a border pixel" << endl;
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // MC trigger parameter
    if( fReader->isGrisuMC() )
    {
        setNTrigger();
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // parallax width cleaning (not default)
    if( fRunPar->fPWmethod > -1 && getImageCleaningParameter() )
    {
        fVImageCleaning->cleanTriggerFixed( getImageCleaningParameter() );
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // set parameters required for image parameter calculation
    fVImageParameterCalculation->setDetectorGeometry( getDetectorGeometry() );
    fVImageParameterCalculation->setParameters( getImageParameters() );
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // parallax width trigger parameter calculation (not default)
    if( fRunPar->fPWmethod > -1 )
    {
        fVImageParameterCalculation->calcTriggerParameters( getTrigger() );
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // image parameter calculation
    fVImageParameterCalculation->calcParameters();
    fVImageParameterCalculation->calcTimingParameters( false );
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // here: no double pass trace integration
    // do a log likelihood image fitting on events on the camera edge only
    if( !isDoublePass() || !hasFADCData() )
    {
        if( getImageParameters()->ntubes > fRunPar->fLogLikelihood_Ntubes_min[getTelID()]
                && getImageParameters()->loss > fRunPar->fLogLikelihoodLoss_min[getTelID()]
                && getImageParameters()->loss < fRunPar->fLogLikelihoodLoss_max[getTelID()] )
        {
            if( !isEqualSummationWindows() )
            {
                fVImageParameterCalculation->setParametersLogL( getImageParameters() );
                setLLEst( fVImageParameterCalculation->calcLL( true, true ) );    // sum2
            }
            fVImageParameterCalculation->setParametersLogL( getImageParameters() );
            setLLEst( fVImageParameterCalculation->calcLL( false, true, isEqualSummationWindows() ) );   // sum
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // muon ring analysis
    if( fRunPar->fmuonmode && !isDoublePass() )
    {
        muonRingAnalysis();
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Hough transform muon ring analysis
    if( fRunPar->fhoughmuonmode && !isDoublePass() )
    {
        houghMuonRingAnalysis();
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // second pass of image calculation using time gradient to adjust sum window (default)
    if( isDoublePass() && hasFADCData() )
    {
        ///////////////////////////////////////////////////////////////////////////////////////////
        // integrate pulses and calculate timing parameters taking time gradients over images into account
        if( fVImageParameterCalculation->getboolCalcGeo() && fVImageParameterCalculation->getboolCalcTiming() )
        {
            calcSecondTZerosSums();
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // smoothing of dead or disabled pixels (first part)
        if( fRunPar->fSmoothDead )
        {
            smoothDeadTubes();
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // image cleaning & gain correction (second pass)
        imageCleaning( false );
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // muon ring analysis (second pass)
        if( fRunPar->fmuonmode )
        {
            muonRingAnalysis();
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        //Hough transform muon ring analysis (second pass)
        
        if( fRunPar->fhoughmuonmode )
        {
            houghMuonRingAnalysis();
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // smoothing of dead or disabled pixels (second part)
        if( fRunPar->fSmoothDead )
        {
            for( unsigned int i = 0; i < getDead().size(); i++ )
            {
                getDead()[i] = savedDead[i];
                getGains()[i] = savedGains[i];
                getDead( true )[i] = savedDeadLow[i];
                getGains( true )[i] = savedGainsLow[i];
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////
        // parallax width trigger parameter calculation
        if( fRunPar->fPWmethod > -1 && getImageCleaningParameter() )
        {
            fVImageCleaning->cleanTriggerFixed( getImageCleaningParameter() );
            fVImageParameterCalculation->calcTriggerParameters( getTrigger() );
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // image parameter calculation
        fVImageParameterCalculation->calcParameters();
        fVImageParameterCalculation->calcTimingParameters( true );
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        // do a log likelihood image fitting on events on the camera edge only
        if( getImageParameters()->ntubes > fRunPar->fLogLikelihood_Ntubes_min[getTelID()]
                && ( fRunPar->fForceLLImageFit ||
                     ( getImageParameters()->loss > fRunPar->fLogLikelihoodLoss_min[getTelID()] &&  getImageParameters()->loss < fRunPar->fLogLikelihoodLoss_max[getTelID()] ) ) ) // FORCELL
        {
            if( !isEqualSummationWindows() )
            {
                fVImageParameterCalculation->setParametersLogL( getImageParameters() );
                setLLEst( fVImageParameterCalculation->calcLL( true, true ) );
            }
            fVImageParameterCalculation->setParametersLogL( getImageParameters() );
            setLLEst( fVImageParameterCalculation->calcLL( false, true, isEqualSummationWindows() ) );
        }
        
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // log likelihood image analysis (not default)
    // calculate image parameters (mainly for edge images)
    if( fRunPar->fImageLL )
    {
        // do this only if geometrical calculation found border/image channels
        if( getImageParameters()->ntubes > 0 )
        {
            fVImageParameterCalculation->setParametersLogL( getImageParametersLogL() );
            setLLEst( fVImageParameterCalculation->calcLL( true ) );
            fVImageParameterCalculation->setParametersLogL( getImageParametersLogL() );
            setLLEst( fVImageParameterCalculation->calcLL( false ) );
        }
        else
        {
            getImageParameters( fRunPar->fImageLL )->reset();
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // fill results into output tree
    fillOutputTree();
}


void VImageAnalyzer::fillOutputTree()
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::fillOutputTree()" << endl;
    }
    
    // fill some run quality histograms
    if( !fReader->isMC() )
    {
        getAnaHistos()->fillL2DiagnosticTree( getRunNumber(), getTelescopeEventNumber( getTelID() ),
                                              0, 0, getFADCstopTZero(), getFADCstopSums() );
    }
    
    // fill some basic run parameters
    getImageParameters()->fsumfirst = getSumFirst();
    getImageParameters()->fsumwindow = getSumWindow();
    getImageParameters()->fsumwindow_2 = getSumWindow_2();
    if( fReader->hasLocalTrigger( getTelID() ) )
    {
        getImageParameters()->fLocalTrigger = 1;
    }
    else
    {
        getImageParameters()->fLocalTrigger = 0;
    }
    //  there is always a trigger in the real data, but not in MC
    if( !fReader->isMC() )
    {
        getImageParameters()->fLocalTrigger = 1;
    }
    getImageParameters()->bad = ( unsigned short int )getNDead();
    getImageParameters()->badLow = ( unsigned short int )getNDead( true );
    // get MC information
    if( fReader->isMC() )
    {
        getImageParameters()->setMC();
        getImageParameters()->MCprimary = fReader->getMC_primary();
        getImageParameters()->MCenergy  = fReader->getMC_energy();
        getImageParameters()->MCxcore   = fReader->getMC_X();
        getImageParameters()->MCycore   = fReader->getMC_Y();
        getImageParameters()->MCxcos    = fReader->getMC_Xcos();
        getImageParameters()->MCycos    = fReader->getMC_Ycos();
        getImageParameters()->MCLocalTriggerTime = fReader->getLocalTriggerTime( getTelID() );
        getImageParameters()->MCLocalDelayedTriggerTime = fReader->getLocalDelayedTriggerTime( getTelID() );
        getImageParameters()->MCTel_Xoff = fReader->getMC_Xoffset();
        getImageParameters()->MCTel_Yoff = fReader->getMC_Yoffset();
    }
    getImageParameters()->eventStatus = getAnalysisTelescopeEventStatus()[getTelID()];
    // fill the trees with the results
    getImageParameters()->fill();
    if( fRunPar->fImageLL )
    {
        getImageParametersLogL()->eventStatus = getAnalysisTelescopeEventStatus()[getTelID()];
        getImageParametersLogL()->fill();
    }
}


void VImageAnalyzer::initAnalysis()
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::initAnalysis()" << endl;
    }
    
    fVImageParameterCalculation->setParameters( getImageParameters() );
    if( fRunPar->fImageLL )
    {
        fVImageParameterCalculation->setParametersLogL( getImageParametersLogL() );
    }
    
    initOutput();
    initTrees();
    
    // temporary vectors for dead pixel smoothing
    savedDead.assign( getDead().size(), false );
    savedDeadLow.assign( getDead( true ).size(), false );
    savedGains.resize( getGains().size(), 1. );
    savedGainsLow.resize( getGains( true ).size(), 1. );
}


void VImageAnalyzer::terminate( bool iDebug_IO )
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::terminate()" << endl;
    }
    if( fOutputfile != 0 && fRunPar->foutputfileName != "-1" )
    {
        fOutputfile->cd();
        fAnaDir[fTelID]->cd();
        // write calibration summaries
        if( getRunParameter()->fsourcetype != 7 )
        {
            getCalibrationData()->terminate( getDead( false ), getDead( true ), getRunParameter()->fTraceIntegrationMethod[fTelID] );
        }
        // write dead channel tree
        // note: this writes the dead channel list of the last event to this tree
        //       (dead channel list is time dependent!)
        TTree* iT = makeDeadChannelTree();
        if( iT )
        {
            iT->Write();
        }
        // write DB pixel data (if available)
        if( getDBPixelDataReader() && getDBPixelDataReader()->getDBStatus() )
        {
            getDBPixelDataReader()->writeDataTree( getTelID() );
        }
        // write histograms
        if( !isMC() )
        {
            getAnaHistos()->terminate( fOutputfile );
        }
        // write pointing data from db to disk (if available)
        if( getTelID() < getPointing().size() && getPointing()[getTelID()] )
        {
            getPointing()[getTelID()]->terminate( isMC() );
        }
        // write main output trees
        if( getImageParameters()->getTree() )
        {
            int i_nbytes = getImageParameters()->getTree()->Write();
            if( iDebug_IO )
            {
                cout << "WRITEDEBUG: tpars trees (nbytes " << i_nbytes << ")" << endl;
            }
        }
        if( fRunPar->fImageLL )
        {
            getImageParametersLogL()->getTree()->Write();
        }
    }
    
    fOutputfile->cd();
}


/*!
     create ROOT output file

     create output directories (for each telescope one, Tel_1, Tel_2, ...)

*/
void VImageAnalyzer::initOutput()
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::initOutput()" << endl;
    }
    // check if root outputfile exist, otherwise create it
    if( fOutputfile != 0 )
    {
        return;
    }
    if( fRunPar->foutputfileName != "-1" )
    {
        // tree versioning numbers used in mscw_energy
        stringstream i_textTitle;
        i_textTitle << "VERSION " << fRunPar->getEVNDISP_TREE_VERSION();
        if( getRunParameter()->fShortTree )
        {
            i_textTitle << "(short tree)";
        }
        fOutputfile = new TFile( fRunPar->foutputfileName.c_str(), "RECREATE", i_textTitle.str().c_str() );
        if( fOutputfile->IsZombie() )
        {
            cout << endl;
            cout << "ERROR: unable to create eventdisplay output file: " << fRunPar->foutputfileName.c_str() << endl;
            cout << "exiting..." << endl;
            cout << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    // creating the directories in the root output file
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        setAnaDir( i );
    }
}


void VImageAnalyzer::setAnaDir( unsigned int iTel )
{
    if( !fOutputfile )
    {
        return;
    }
    // create output directories in root files
    TDirectory* i_Directory = gDirectory;
    char i_rootDir[200];
    sprintf( i_rootDir, "Tel_%d", ( int )iTel + 1 );
    char i_rootDirTit[200];
    sprintf( i_rootDirTit, "analysis results for telescope %d (type %lld)", ( int )iTel + 1, getDetectorGeometry()->getTelType()[iTel] );
    fOutputfile->cd();
    if( !fOutputfile->FindObject( i_rootDir ) )
    {
        setRootDir( iTel, fOutputfile->mkdir( i_rootDir, i_rootDirTit ) );
    }
    i_Directory->cd();
}


/*!
   initialize the output trees

*/
void VImageAnalyzer::initTrees()
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::initTrees() " <<  endl;
    }
    
    if( !fOutputfile->cd() )
    {
        cout << "VImageAnalyzer::initTrees() unable to enter output file " << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fAnaDir[getTelID()]->cd();
    if( getTelID() < fAnaDir.size() && fAnaDir[getTelID()] && !fAnaDir[getTelID()]->cd() )
    {
        cout << "VImageAnalyzer::initTrees() unable to enter directory: " << fAnaDir[getTelID()]->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else if( getTelID() >= fAnaDir.size() )
    {
        cout << "VImageAnalyzer::initTrees() error: analysis directory vector not large enough" << endl;
        cout << "\t should be " << getNTel() << ", is " << fAnaDir.size() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else if( !fAnaDir[getTelID()] )
    {
        cout << "VImageAnalyzer::initTrees() error: analysis directory not set" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    // now book the trees (for MC with additional MC block)
    // tree versioning numbers used in mscw_energy
    char i_text[300];
    sprintf( i_text, "tpars" );
    ostringstream iSTRText;
    iSTRText << "Event Parameters (Telescope " << getTelID() + 1;
    iSTRText << ", VERSION " << fRunPar->getEVNDISP_TREE_VERSION() << ")";
    if( getRunParameter()->fShortTree )
    {
        iSTRText << " (short tree)";
    }
    fVImageParameterCalculation->getParameters()->initTree( i_text, iSTRText.str().c_str(), fReader->isMC(), false, fRunPar->fmuonmode, fRunPar->fhoughmuonmode );
    
    // for log likelihood method, book a second image parameter tree
    if( fRunPar->fImageLL )
    {
        sprintf( i_text, "lpars" );
        char i_textTitle[300];
        sprintf( i_textTitle, "Event Parameters, loglikelihood (Telescope %d)", getTelID() + 1 );
        fVImageParameterCalculation->getLLParameters()->initTree( i_text, i_textTitle, fReader->isMC(), true, fRunPar->fmuonmode, fRunPar->fhoughmuonmode );
    }
}


bool VImageAnalyzer::initEvent()
{
    if( fDebug )
    {
        cout << "VImageAnalyzer::initEvent" << endl;
    }
    setSums( 0. );
    setPulseTiming( 0., true );
    setPulseTiming( 0., false );
    setTCorrectedSumFirst( getSumFirst() );
    setTCorrectedSumLast( getSumFirst() + getSumWindow() );
    setCurrentSummationWindow( getSumWindow(), false );
    setCurrentSummationWindow( getSumWindow_2(), true );
    setImage( false );
    setBorder( false );
    setImageBorderNeighbour( false );
    
    // set all variables in getImageParameters() to zero
    getImageParameters()->reset();
    if( fRunPar->fImageLL )
    {
        getImageParameters( fRunPar->fImageLL )->reset();
    }
    
    //!Initialize this event
    getImageParameters()->fTelID = getTelID();
    
    getImageParameters()->MJD  = getEventMJD();
    getImageParameters()->time = getEventTime();
    
    // get time since run start
    if( getAnaData()->fTimeSinceRunStart < 0. )
    {
        getAnaData()->fTimeSinceRunStart = 0.;
        getAnaData()->fTimeRunStart = getImageParameters()->time;
    }
    else
    {
        double idiff = getImageParameters()->time - getAnaData()->fTimeRunStart;
        if( idiff >= 0. )
        {
            getAnaData()->fTimeSinceRunStart = idiff;
        }
        else
        {
            getAnaData()->fTimeSinceRunStart = 86400. - idiff;
        }
    }
    
    getImageParameters()->eventNumber = getTelescopeEventNumber( getTelID() );
    getImageParameters()->nsamples = getNSamples();
    getImageParameters()->runNumber = getRunNumber();
    getImageParameters()->eventType = ( unsigned short int )getReader()->getNewEventType();
    
    getImageParameters()->fTrig_type = fReader->getLocalTriggerType( getTelID() );
    if( fReader->isGrisuMC() )
    {
        setNTrigger();
    }
    
    // fill high/low gain vector
    getImageParameters()->nlowgain = fillHiLo();
    // fill vector with zero suppressed channels
    // (return number of zero suppressed channels)
    getImageParameters()->nzerosuppressed = fillZeroSuppressed();
    
    if( fRunPar->fImageLL )
    {
        getImageParametersLogL()->time = getImageParameters()->time;
        getImageParametersLogL()->MJD = getImageParameters()->MJD;
        getImageParametersLogL()->eventNumber = getImageParameters()->eventNumber;
        getImageParametersLogL()->nsamples = getImageParameters()->nsamples;
        getImageParametersLogL()->runNumber = getImageParameters()->runNumber;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // get the telescope pointing
    //////////////////////////////////////////////////////////////////////////////////////////////
    // MC readers
    // source file is Grisu MC, MC DST, or MC Pe
    if( fReader->isMC() )
    {
        getPointing()[getTelID()]->setMC();
        if( getTelID() < fReader->getTelElevation().size() )
        {
            getPointing()[getTelID()]->setTelElevation( fReader->getTelElevation()[getTelID()] );
        }
        else
        {
            getPointing()[getTelID()]->setTelElevation( 0. );
        }
        if( getTelID() < fReader->getTelAzimuth().size() )
        {
            getPointing()[getTelID()]->setTelAzimuth( fReader->getTelAzimuth()[getTelID()] );
        }
        else
        {
            getPointing()[getTelID()]->setTelAzimuth( 0. );
        }
    }
    // set pointing direction from command line
    else if( getRunParameter()->felevation > 0. && getRunParameter()->fazimuth > 0. )
    {
        getPointing()[getTelID()]->setTelElevation( getRunParameter()->felevation );
        getPointing()[getTelID()]->setTelAzimuth( getRunParameter()->fazimuth );
    }
    // set pointing direction with target
    // target is set via command line, getPointing() is initiated in run/VEventLoop::initEventLoop()
    else
    {
        bool iSet = true;
        if( getPointing()[getTelID()]->isSet() && getPointing()[getTelID()]->getTargetName() != "laser" )
        {
            // this calculates telescope elevation and azimuth from telescope pointing in ra and dec
            getPointing()[getTelID()]->setTelPointing( getImageParameters()->MJD, getImageParameters()->time, getRunParameter()->fDBTracking );
            if( !getPointing()[getTelID()]->isPrecessed() )
            {
                getPointing()[getTelID()]->precessTarget( getImageParameters()->MJD, getTelID() );
                // set wobble offsets
                getPointing()[getTelID()]->setWobbleOffset( getRunParameter()->fWobbleNorth, getRunParameter()->fWobbleEast, getTelID(), getImageParameters()->MJD );
                getPointing()[getTelID()]->setTelPointing( getImageParameters()->MJD, getImageParameters()->time );
            }
        }
        else
        {
            iSet = false;
        }
        if( !iSet )
        {
            if( fNTel != 1 && getTelescopeEventNumber( getTelID() ) < 1 )
            {
                cout << "VImageAnalyzer::initEvent(): no telescope pointing" << endl;
                setNoPointing( true );
            }
            getPointing()[getTelID()]->setTelElevation( 0. );
            getPointing()[getTelID()]->setTelAzimuth( 0. );
        }
    }
    
    return true;
}

/*!
This reduces the effect of dead tubes by assigning them charge and timing values
derived by taking the average of the neighbouring tubes
*/
void VImageAnalyzer::smoothDeadTubes()
{
    unsigned int i_nchannel = getNChannels();
    
    // reset vectors
    for( unsigned int i = 0; i < getDead().size(); i++ )
    {
        getDeadRecovered()[i] = false;
        savedDead[i] = getDead()[i];
        savedGains[i] = getGains()[i];
        savedDeadLow[i] = getDead( true )[i];
        savedGainsLow = getGains( true )[i];
    }
    
    // loop over all channels
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        // select dead channels
        if( getDetectorGeo()->getAnaPixel()[i] < 1 || !getDead()[i] )
        {
            continue;
        }
        // exclude low gain channel from this estimation
        // (GM) not clear why low-gain channels are excluded
        // if( getHiLo()[i] )
        // {
        //	continue;
        // }
        
        double ave_sum = 0.;
        double ave_pedvar = 0.;
        double ave_gain = 0.;
        double ave_tzero = 0.;
        double ave_toffvars = 0.;
        double count = 0.;
        // loop over all neighbours of this dead pixel
        for( unsigned int j = 0; j < getDetectorGeo()->getNeighbours()[i].size(); j++ )
        {
            unsigned int k = getDetectorGeo()->getNeighbours()[i][j];
            
            if( getDead()[k] )
            {
                continue;
            }
            
            ave_sum += getSums()[k];
            ave_pedvar += getPedvars()[k];
            ave_gain += getGains()[k];
            ave_tzero += getTZeros()[k];
            ave_toffvars += getTOffsetvars()[k];
            count += 1.;
        }
        if( count > 0. )
        {
            ave_sum /= count;
            ave_pedvar /= count;
            ave_gain /= count;
            ave_tzero /= count;
            ave_toffvars /= count;
            
            // this does not work when dead channel list is time dependent
            getSums()[i] = ave_sum;
            getPedvars()[i] = ave_pedvar;
            getGains()[i] = ave_gain;
            getTZeros()[i] = ave_tzero;
            getTOffsetvars()[i] = ave_toffvars;
            getDead()[i] = 0;
            getDeadRecovered()[i] = true;
        }
    }
    
}


/*!
   this is very preliminary and ugly

   useful debugging routine (and for some Leeds project students)
*/
void VImageAnalyzer::printTrace( int iChannel )
{
    // calculate mean sum over all pixels
    double i_total = 0;
    for( unsigned int i = 0; i < getNChannels(); i++ )
    {
        i_total += getSums()[i];
    }
    //if( i_total < 30000 ) return;
    double i_totSum = 0;
    double i_totSumN = 0;
    unsigned int nhits = fReader->getNumChannelsHit();
    for( unsigned int i = 0; i < nhits; i++ )
    {
        unsigned int chanID = 0;
        try
        {
            chanID = fReader->getHitID( i );
            // exclude channels with fadc stop signal
            if( getTOffsets()[chanID] > 6. )
            {
                getDead()[chanID] = 8;
            }
            if( !getDead()[chanID] )
            {
                i_totSum += getSums()[chanID];
                i_totSumN++;
            }
        }
        catch( ... )
        {
            cout << "VImageAnalyzer::printTrace(), index out of range (fReader->getHitID) ";
            cout << i << "(Telescope " << getTelID() + 1 << ", event " << getEventNumber() << ")" << endl;
            continue;
        }
    }
    if( i_totSumN > 0. )
    {
        i_totSum /= i_totSumN;
    }
    
    // now print everything out
    cout << getTelescopeEventNumber( getTelID() ) << " " << getSums()[fReader->getHitID( iChannel )] << " " << i_totSum << "    ";
    unsigned int chanID = fReader->getHitID( iChannel );
    fReader->selectHitChan( chanID );
    if( fReader->getNumSamples() > 0 )
    {
        for( unsigned int j = 0; j < fReader->getNumSamples(); j++ )
        {
            if( fReader->has16Bit() )
            {
                cout << fReader->getSamplesVec16Bit()[j] << " ";
            }
            else
            {
                cout << fReader->getSamplesVec()[j] << " ";
            }
        }
    }
    cout << endl;
    
}


void VImageAnalyzer::shutdown()
{
    if( fDebug )
    {
        cout << "void VImageAnalyzer::shutdown()" << endl;
    }
    if( fOutputfile && fOutputfile->IsOpen() )
    {
        fOutputfile->Flush();
        cout << "closing evndisp output file, final contents: ";
        if( fReader->isMC() )
        {
            fOutputfile->ls();
        }
        cout << "\t file size: " << fOutputfile->GetSize() << endl;
        fOutputfile->Close();
    }
}


/*!
   calculate trigger vector

   (needed for grisu MC only)
*/
void VImageAnalyzer::setNTrigger()
{
    if( getDebugFlag() )
    {
        cout << "VImageAnalyzer::setNTrigger" << endl;
    }
    getImageParameters()->ntrig = getReader()->getNumberofFullTrigger();
    
    if( getRunParameter()->fsourcetype == 6 )
    {
        return;
    }
    
    std::vector<bool> triggered = getReader()->getFullTrigVec();
    unsigned int triggered_size = triggered.size();
    unsigned short max_num_in_patch = 0;
    
    for( unsigned int i = 0; i < triggered_size; i++ )
    {
        if( i < getDetectorGeo()->getNumChannels() && triggered[i] )
        {
            //! find position of triggered tube
            float xi = getDetectorGeo()->getX()[i];
            float yi = getDetectorGeo()->getY()[i];
            unsigned short num_in_patch = 1;
            //! see how many other triggered tubes are within 0.3 degrees
            for( unsigned int j = 0; j < triggered_size; j++ )
            {
                if( triggered[j] && j != i )
                {
                    bool i_bBreak = false;
                    for( unsigned int t = 0; t < getFADCstopTrig().size(); t++ )
                    {
                        if( j == getFADCstopTrig()[t] )
                        {
                            i_bBreak = true;
                        }
                    }
                    if( i_bBreak )
                    {
                        continue;
                    }
                    float xj = getDetectorGeo()->getX()[j];
                    float yj = getDetectorGeo()->getY()[j];
                    if( ( xj - xi ) * ( xj - xi ) + ( yj - yi ) * ( yj - yi ) < 0.09 )
                    {
                        num_in_patch += 1;
                    }
                }
            }
            if( num_in_patch > max_num_in_patch )
            {
                max_num_in_patch = num_in_patch;
            }
        }
    }
    getImageParameters()->ntrig_per_patch = max_num_in_patch;
}

/*
 *   run the image cleaning
 *
 *   - select the correct image cleaning method and pass parameters
 *
 *   - TODO: simplify the interface to the image cleaning class by passing the
 *           pointer to the image cleaning run parameters
 *
 */
void VImageAnalyzer::imageCleaning( bool iDoublePassParameters )
{
    if( !getImageCleaningParameter( iDoublePassParameters ) )
    {
        return;
    }
    
    /////////////////////////////
    // fixed threshold cleaning
    if( getImageCleaningParameter( iDoublePassParameters )->fUseFixedThresholds )
    {
        gainCorrect();
        if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "TIMECLUSTERCLEANING" )
        {
            fVImageCleaning->cleanImageFixedWithTiming( getImageCleaningParameter( iDoublePassParameters ) );
        }
        // time-next-neighbour cleaning
        else if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "TIMENEXTNEIGHBOUR" )
        {
            fVImageCleaning->cleanNNImageFixed( getImageCleaningParameter( iDoublePassParameters ) );
        }
        //cluster cleaning
        else if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "CLUSTERCLEANING" )
        {
            fVImageCleaning->cleanImageWithClusters( getImageCleaningParameter( iDoublePassParameters ), true );
        }
        // fixed cleaning levels (classic image/border)
        else
        {
            fVImageCleaning->cleanImageFixed( getImageCleaningParameter( iDoublePassParameters ) );
        }
    }
    /////////////////////////////
    // signal/noise cleaning thresholds
    else
    {
        if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "TIMECLUSTERCLEANING" )
        {
            fVImageCleaning->cleanImagePedvarsWithTiming( getImageCleaningParameter( iDoublePassParameters ) );
        }
        else if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "TWOLEVELANDCORRELATION" )
        {
            fVImageCleaning->cleanImageTraceCorrelate( getImageCleaningParameter( iDoublePassParameters ) );
        }
        // simple time two-level cleaning
        else if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "TIMETWOLEVEL" )
        {
            fVImageCleaning->cleanImagePedvarsTimeDiff( getImageCleaningParameter( iDoublePassParameters ) );
        }
        //cluster cleaning
        else if( getImageCleaningParameter( iDoublePassParameters )->getImageCleaningMethod() == "CLUSTERCLEANING" )
        {
            fVImageCleaning->cleanImageWithClusters( getImageCleaningParameter( iDoublePassParameters ), false );
        }
        else
        {
            fVImageCleaning->cleanImagePedvars( getImageCleaningParameter( iDoublePassParameters ) );
        }
        gainCorrect();
    }
}

void VImageAnalyzer::muonRingAnalysis()
{
    fVImageParameterCalculation->muonRingFinder();
    fVImageParameterCalculation->muonPixelDistribution();
    fVImageParameterCalculation->sizeInMuonRing();
}

void VImageAnalyzer::houghMuonRingAnalysis()
{

    // Iterative fit muon analysis
    fVImageParameterCalculation->muonRingFinder();
    
    //Hough transform based muon parametrization algorithm invoked here
    //fVImageParameterCalculation->houghMuonRingFinder();
    
    // Iterative fit muon analysis
    fVImageParameterCalculation->muonPixelDistribution();
    
    //Hough transform muon ID technique
    fVImageParameterCalculation->houghMuonPixelDistribution();
    
    //Hough transform based size calculation algorithm
    //fVImageParameterCalculation->houghSizeInMuonRing();
    
    // Iterative fit muon analysis
    fVImageParameterCalculation->sizeInMuonRing();
    
}


