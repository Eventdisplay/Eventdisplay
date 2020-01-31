/*! \class VStereoAnalysis
    \brief class for producing histograms from parameterized stereo VERITAS data

*/

#include "VStereoAnalysis.h"

VStereoAnalysis::VStereoAnalysis( bool ion, string i_hsuffix, VAnaSumRunParameter* irunpara, vector< TDirectory* > iDirRun,
                                  TDirectory* iDirTot, string iDataDir, int iRandomSeed, bool iTotalAnalysisOnly )
{
    fDebug = false;
    
    fDataFile = 0;
    fInstrumentEpoch = "NOT_SET";
    fDirTot = iDirTot;
    fDirTotRun = iDirRun;
    bTotalAnalysisOnly = iTotalAnalysisOnly;
    
    // set default values
    fIsOn = ion;
    // do full sky plots
    fNoSkyPlots = false;
    
    gMeanEffectiveArea = 0;
    gTimeBinnedMeanEffectiveArea = 0;
    gMeanEsys_MC = 0;
    
    gMeanEffectiveAreaMC = 0;

    fHisCounter = 0;
    fTotCount = 0.;
    
    fMeanAzimuth = 0.;
    fMeanElevation = 0.;
    fNMeanElevation = 0.;
    
    fTreeSelectedEvents = 0;
    fCDataTreeClone = 0;
    
    fRunPara = irunpara;
    fTreeWithEventsForCtools = 0 ; // WRITEEVENTTREEFORCTOOLS
    fDeadTimeStorage = 0.0 ;
    
    fVsky = new VSkyCoordinates() ;
    fVsky->supressStdoutText( true ) ;
    fVsky->setObservatory( VGlobalRunParameter::getObservatory_Longitude_deg(), VGlobalRunParameter::getObservatory_Latitude_deg() );
    
    // calculating run start, end and duration (verifies data trees)
    if( !bTotalAnalysisOnly )
    {
        setRunTimes();
    }
    
    // targets and exclusion regions
    if( !bTotalAnalysisOnly )
    {
        defineAstroSource();
    }
    
    ///////////////////////////////
    // define histograms
    
    // combined results
    iDirTot->cd();
    fHistoTot = new VStereoHistograms( i_hsuffix, fRunPara->fSkyMapBinSize, fRunPara->fSkyMapBinSizeUC,
                                       fRunPara->fEnergySpectrumBinSize, fRunPara->fTimeIntervall, -1, -1, fIsOn );
    fHistoTot->setSkyMapSize( fRunPara->fSkyMapSizeXmin, fRunPara->fSkyMapSizeXmax, fRunPara->fSkyMapSizeYmin, fRunPara->fSkyMapSizeYmax );
    
    // one set of histograms for each run
    if( iDirRun.size() != fRunPara->fRunList.size() )
    {
        cout << "VStereoAnalysis::VStereoAnalysis fatal error, directory and run list different ";
        cout << iDirRun.size() << "\t" << fRunPara->fRunList.size() << endl;
        exit( EXIT_FAILURE );
    }
    // define histograms and rate counters
    vector< double > i_v;
    for( unsigned int i = 0; i < fRunPara->fRunList.size(); i++ )
    {
        iDirRun[i]->cd();
        fHisto.push_back( new VStereoHistograms( i_hsuffix, fRunPara->fSkyMapBinSize, fRunPara->fSkyMapBinSizeUC,
                          fRunPara->fEnergySpectrumBinSize, fRunPara->fTimeIntervall,
                          f_t_in_s_min[fIsOn ? fRunPara->fRunList[i].fRunOn : fRunPara->fRunList[i].fRunOff],
                          f_t_in_s_max[fIsOn ? fRunPara->fRunList[i].fRunOn : fRunPara->fRunList[i].fRunOff], fIsOn ) );
        fHisto.back()->setSkyMapSize( fRunPara->fSkyMapSizeXmin, fRunPara->fSkyMapSizeXmax,
                                      fRunPara->fSkyMapSizeYmin, fRunPara->fSkyMapSizeYmax );
        if( fIsOn )
        {
            fHisto.back()->setRunNumber( fRunPara->fRunList[i].fRunOn );
        }
        else
        {
            fHisto.back()->setRunNumber( fRunPara->fRunList[i].fRunOff );
        }
        
        // define dead time calculators
        fDeadTime.push_back( new VDeadTime( fIsOn ) );
        
        // rate plots
        fRateCounts.push_back( i_v );
        fRateTime.push_back( i_v );
        fRateTimeIntervall.push_back( i_v );
        
    }
    
    // define the time mask
    fTimeMask = new VTimeMask();
    
    // define the cuts
    fCuts = new VGammaHadronCuts();
    char hname[200];
    if( fIsOn )
    {
        sprintf( hname, "GammaHadronCuts" );
    }
    else
    {
        sprintf( hname, "GammaHadronCuts_off" );
    }
    fCuts->SetName( hname );
    fCuts->resetCutValues();
    fCuts->setDataTree( 0 );
    fCuts->setDataDirectory( iDataDir );
    fCuts->setReconstructionType( fRunPara->fReconstructionType );
    // define the background model
    fMap   = new VStereoMaps( false, iRandomSeed, fRunPara->fTMPL_RE_RemoveOffRegionsRandomly );
    fMapUC = new VStereoMaps( true,  iRandomSeed, fRunPara->fTMPL_RE_RemoveOffRegionsRandomly );
}

VStereoAnalysis::~VStereoAnalysis()
{
    if( fVsky )
    {
        delete fVsky;
    }
    if( fHistoTot )
    {
        delete fHistoTot;
    }
    for( unsigned int i = 0; i < fHisto.size(); i++ )
    {
        if( fHisto[i] )
        {
            delete fHisto[i];
        }
    }
    
    for( unsigned int i = 0; i < fDeadTime.size(); i++ )
    {
        if( fDeadTime[i] )
        {
            delete fDeadTime[i];
        }
    }
    if( fCuts )
    {
        delete fCuts;
    }
    if( fMap )
    {
        delete fMap;
    }
    if( fMapUC )
    {
        delete fMapUC;
    }
}


/*!
 *
 * establish run times (mean, start, end and duration) for a list of runs
 *
 */
void VStereoAnalysis::setRunTimes()
{
    cout << endl << "-----------------------------------------------------------------------" << endl;
    cout << "Checking data trees " << ( fIsOn ? "(ON runs)" : "(OFF runs)" ) << endl;
    cout << "\t Run \t| Start (MJD : secs)\t| End (MJD : secs)\t| Duration (secs [mins])" << endl;
    
    for( unsigned int i = 0; i < fRunPara->fRunList.size(); i++ )
    {
        int i_run = fIsOn ? fRunPara->fRunList[i].fRunOn : fRunPara->fRunList[i].fRunOff;
        
        CData* c = getDataFromFile( i_run );
        
        cout << setRunTimes( c ) << endl;
        if( fIsOn )
        {
            fRunPara->fRunList[i].fMJDOn = getMJD( i_run );
        }
        else
        {
            fRunPara->fRunList[i].fMJDOff = getMJD( i_run );
        }
        
        closeDataFile();
    }
}


/*!
 *
 * establish run times (mean, start, end and duration) of a single run
 *
 */
string  VStereoAnalysis::setRunTimes( CData* iData )
{
    ostringstream ires( "" );
    ires.setf( ios_base::fixed, ios_base::floatfield );
    
    double i_min, i_minMJD, i_minUTC = 0.;
    double i_max, i_maxMJD, i_maxUTC = 0.;
    double i_dur = 0.;
    
    fDataRun = iData;
    int i_run = getDataRunNumber();
    
    fDataRun->GetEntry( 1 );
    i_min = fDataRun->Time;
    f_t_in_s_min[i_run] = i_min;
    i_minMJD = fDataRun->MJD;
    i_minUTC = VSkyCoordinatesUtilities::getUTC( ( int )i_minMJD, i_min );
    
    int i_nentries = ( int )fDataRun->fChain->GetEntries() - 2;
    fDataRun->GetEntry( i_nentries );
    i_max = fDataRun->Time;
    f_t_in_s_max[i_run] = i_max;
    i_maxMJD = fDataRun->MJD;
    i_maxUTC = VSkyCoordinatesUtilities::getUTC( ( int )i_maxMJD, i_max );
    
    i_dur = ( i_maxUTC - i_minUTC ) * 24 * 60 * 60;
    
    fRunMJDStart[i_run] = i_minUTC;
    fRunMJDStopp[i_run] = i_maxUTC;
    fRunMJD[i_run] = ( i_maxUTC + i_minUTC ) / 2.;
    fRunDuration[i_run] = i_dur;
    
    ires.precision( 0 );
    ires << "\t " << i_run << "\t| ";
    
    ires <<  i_minMJD;
    ires.precision( 2 );
    ires << " : " << i_min;
    
    ires.precision( 0 );
    ires << "\t| " << i_maxMJD;
    ires.precision( 2 );
    ires << " : " << i_max;
    
    ires << "\t| " << i_dur  << " [" << i_dur / 60. << "]";
    
    return ires.str();
}


/*
 *
 * Get run number from current data tree
 *
 */
int VStereoAnalysis::getDataRunNumber() const
{
    if( !fDataRun )
    {
        cout << "VStereoAnalysis::getDataRunNumber error: no data tree." << endl;
    }
    else
    {
        if( fDataRun->GetEntry( 0 ) == -1 )
        {
            cout << "VStereoAnalysis::getDataRunNumber error: can't seem to access tree." << endl;
        }
        else if( fDataRun->GetEntry( 0 ) == 0 )
        {
            cout << "VStereoAnalysis::getDataRunNumber error: tree is empty." << endl;
            fDataRun->fChain->Print();
        }
        else
        {
            return fDataRun->runNumber;
        }
    }
    
    exit( EXIT_FAILURE );
    return 0;
}


/*!
 *
 *  fill all histograms and maps
 *  (main event loop)
 *
 *  int irun  runnumber to be analyzed, set irun = -1 to analyze all runs in data chain
 *
 *  return number of events passing all cuts
 *
 */
double VStereoAnalysis::fillHistograms( int icounter, int irun, double iAzMin, double iAzMax, double iPedVar )
{
    if( fDebug )
    {
        cout << "DEBUG double VStereoAnalysis::fillHistograms() "  << icounter << "\t" << irun << endl;
    }
    // this flag increases significantly the debug output
    bool fDebugCuts = false;
    
    fHisCounter = icounter;
    
    ////////////////////////////////////////////////
    // combine all histograms from all runs
    if( irun < 0 )
    {
        fHisCounter = -1;
        return combineHistograms();
    }
    ////////////////////////////////////////////////
    
    ////////////////////////////////////////////////
    // analyze individual run
    if( fIsOn )
    {
        cout << endl << "------------------------------------------------------------------------" << endl;
        cout << "Filling ON histograms for run " << irun << " -----------------------------" << endl;
    }
    else
    {
        cout << endl << "------------------------------------------------------------------------" << endl;
        cout << "Filling OFF histograms for run " << irun << " -----------------------------" << endl;
    }
    
    // set pointer to data tree (run wise)
    fDataRun = getDataFromFile( irun );
    if( fDataRun == 0 || fDataRunTree == 0 )
    {
        cout << "VStereoAnalysis::fillHistograms error, no data tree " << endl;
        cout << "\t" << fDataRun << "\t" << fDataRunTree << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( fHisCounter > ( int )fHisto.size() )
    {
        cout << "VStereoAnalysis::fillHistograms invalid run number " << irun << "\t" << fHisCounter << "\t" << fHisto.size() << endl;
        exit( EXIT_FAILURE );
    }
    
    double iMJDStart = 0.;
    double iMJDStopp = 0.;
    if( getDataRunNumber() != irun )
    {
        cout << "VStereoAnalysis::fillHistograms warning: given run (" << irun << ") does not match run of given tree (" << getDataRunNumber() << ")" << endl;
    }
    else
    {
        if( fRunMJDStart.find( irun ) == fRunMJDStart.end() )
        {
            setRunTimes( fDataRun );
        }
        iMJDStart = fRunMJDStart[irun];
        if( fRunMJDStopp.find( irun ) == fRunMJDStopp.end() )
        {
            setRunTimes( fDataRun );
        }
        iMJDStopp = fRunMJDStopp[irun];
    }
    //////////////////////////////////////////
    // boolean for gamma/hadron cuts
    
    // event is gamma-ray like according to VGammaHadronCuts
    bIsGamma = false;
    // event direction is inside search region (e.g. reflected region)
    bool bDirectionCuts = false;
    
    // rate vectors
    vector< double > iRateCounts;
    vector< double > iRateTime;
    vector< double > iRateTimeIntervall;
    
    // initialize time mask
    fTimeMask->setMask( irun, iMJDStart, iMJDStopp, fRunPara->fTimeMaskFile );
    
    // initialize cuts
    setCuts( fRunPara->fRunList[fHisCounter], irun );
    fDataRun->setReconstructionType( fCuts->fReconstructionType );
    
    // define histograms
    fDirTotRun[fHisCounter]->cd();
    fHisto[fHisCounter]->setRunNumber( irun );
    fHisto[fHisCounter]->defineHistograms();
    
    // adjust dead time calculator
    TDirectoryFile* iDeadtimeDirectory = 0;
    if( fDataFile )
    {
        iDeadtimeDirectory = ( TDirectoryFile* )fDataFile->Get( "deadTimeHistograms" );
    }
    if( iDeadtimeDirectory )
    {
        fDeadTime[fHisCounter]->readHistograms( iDeadtimeDirectory );
    }
    else
    {
        fDeadTime[fHisCounter]->defineHistograms();
    }
    
    // adjust axis in rate histograms
    fHisto[fHisCounter]->makeRateHistograms( iMJDStart, iMJDStopp );
    
    // set map properties and exclusion regions
    fMap->setData( fDataRun );
    fMap->setRunList( fRunPara->fRunList[fHisCounter] );
    fMap->setTargetShift( fRunPara->fRunList[fHisCounter].fTargetShiftWest, fRunPara->fRunList[fHisCounter].fTargetShiftNorth );
    fMap->setNoSkyPlots( fNoSkyPlots );
    fMap->setRegionToExclude( fRunPara->getExclusionRegions( fHisCounter ) );
    fMap->setHistograms( fHisto[fHisCounter]->hmap_stereo, fHisto[fHisCounter]->hmap_alpha, fHisto[fHisCounter]->hmap_MeanSignalBackgroundAreaRatio );
    
    fMapUC->setData( fDataRun );
    fMapUC->setRunList( fRunPara->fRunList[fHisCounter] );
    fMapUC->setTargetShift( fRunPara->fRunList[fHisCounter].fTargetShiftWest, fRunPara->fRunList[fHisCounter].fTargetShiftNorth );
    fMapUC->setNoSkyPlots( fNoSkyPlots );
    fMapUC->setRegionToExclude( fRunPara->getExclusionRegions( fHisCounter ) );
    fMapUC->setHistograms( fHisto[fHisCounter]->hmap_stereoUC, fHisto[fHisCounter]->hmap_alphaUC, 0 );
    
    // initialize gamma/hadron cuts
    fCuts->setDataTree( fDataRun );
    
    // tree with selected events
    init_TreeWithSelectedEvents( irun, fIsOn );
    
    if( fIsOn && fRunPara->fWriteEventTreeForCtools )  // WRITEEVENTTREEFORCTOOLS
    {
        init_TreeWithEventsForCtools( irun );
    }
    
    // spectral energy reconstruction (effective areas, etc.)
    // effective area class
    VEffectiveAreaCalculator fEnergy( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, iAzMin, iAzMax, iPedVar,
                                      fRunPara->fEnergyReconstructionSpectralIndex, fRunPara->fMCZe,
                                      fRunPara->fEnergyEffectiveAreaSmoothingIterations,
                                      fRunPara->fEnergyEffectiveAreaSmoothingThreshold, fRunPara->fEffectiveAreaVsEnergyMC,
                                      fRunPara->fLikelihoodAnalysis,
                                      fIsOn);
    double iEnergyWeighting = 1.;
    double iErec = 0.;
    double iErecChi2 = 0.;
    double iPedVar_temp = 0.;
    double iXoff = 0.;
    double iYoff = 0.;
    // for time-check save old-time and new-time
    // variable to set the real duration of each time bin
    double time_of_EVENT = 0;
    int index_time_bin_NOW  = 1;
    
    double i_UTC = 0.;
    double i_xderot = -99.;
    double i_yderot = -99.;
    double i_theta2 = -99.;
    
    // mean direction values
    fMeanAzimuth = 0.;
    fMeanElevation = 0.;
    fNMeanElevation = 0.;
    double iDirectionOffset = 0.;
    
    // get number of entries from data tree
    Int_t nentries = Int_t( fDataRun->fChain->GetEntries() );
    if( fDebug )
    {
        cout << "DEBUG double VStereoAnalysis::fillHistograms() reading chain " << fDataRun->fChain->GetName() << "\t" << nentries << endl;
    }
    cout << "\t number of entries for this run: " << nentries << endl;
    
    double i_count = 0.;
    int nentries_run = 0;
    
    /////////////////////////////////////////////////////////////////////
    // loop over all entries/events in the data tree
    for( int i = 0; i < nentries; i++ )
    {
        fDataRun->GetEntry( i );
        
        if( fDebugCuts )
        {
            cout << "DEBUGCUTS: ENTRY " << i << " / " << fDataRun->eventNumber << "\t" << fDataRun->runNumber << "\t" << irun << endl;
        }
        
        if( fDataRun->runNumber == irun )
        {
            // count how many entries are in this run
            nentries_run++;
            
            // UTC time
            i_UTC = VSkyCoordinatesUtilities::getUTC( fDataRun->MJD, fDataRun->Time );
            
            // make time cut
            if( !fTimeMask->checkAgainstMask( i_UTC ) )
            {
                if( fDebugCuts )
                {
                    cout << "DEBUGCUTS: time mask cut " << i << " / " << fDataRun->eventNumber << "\t UTC " << i_UTC << endl;
                }
                continue;
            }
            
            // fill rate histograms
            fHisto[fHisCounter]->hrate_1sec->Fill( i_UTC );
            fHisto[fHisCounter]->hrate_10sec->Fill( i_UTC );
            fHisto[fHisCounter]->hrate_1min->Fill( i_UTC );
            
            // dead time calculation
            if( !iDeadtimeDirectory )
            {
                fDeadTime[fHisCounter]->fillDeadTime( fDataRun->Time );
            }
            
            // get energy (depending on energy reconstruction method)
            iErec     = fDataRun->getEnergy_TeV( );
            iErecChi2 = fDataRun->getEnergyChi2( );
            // get shower direction (depending on shower reconstruction method)
            iXoff = fDataRun->getXoff();
            iYoff = fDataRun->getYoff();
            // direction offset
            iDirectionOffset = sqrt( iXoff * iXoff + iYoff * iYoff );
            
            ////////////////////////////////////////////////
            // apply all quality cuts
            //
            // check if event is outside fiducial area
            if( !fCuts->applyInsideFiducialAreaCut() )
            {
                if( fDebugCuts )
                {
                    cout << "DEBUGCUTS: applyInsideFiducialAreaCut " << i << " / " << fDataRun->eventNumber << "\t";
                    cout << "offset " << iDirectionOffset << endl;
                }
                continue;
            }
            
            // stereo quality cuts (e.g. successful direction, mscw, mscl reconstruction)
            if( !fCuts->applyStereoQualityCuts( fRunPara->fEnergyReconstructionMethod, false, i , fIsOn ) )
            {
                if( fDebugCuts )
                {
                    cout << "DEBUGCUTS: applyStereoQualityCuts " << i << " / " << fDataRun->eventNumber << endl;
                }
                continue;
            }
            
            
            // require valid effective area valid for this event
            // effective areas depend on ze, wobble offset, pedestal variation, etc
            if( fDataRun->meanPedvar_Image > 0. )
            {
                iPedVar_temp = fDataRun->meanPedvar_Image;
            }
            else
            {
                iPedVar_temp = iPedVar;
            }
            
            // fill image and trigger pattern histograms
            fHisto[fHisCounter]->hTriggerPatternBeforeCuts->Fill( fDataRun->LTrig );
            fHisto[fHisCounter]->hImagePatternBeforeCuts->Fill( fDataRun->getImgSel() );
            
            // derotate coordinates
            getDerotatedCoordinates( icounter, fDataRun, i_xderot, i_yderot );
            
            // gamma/hadron cuts && successful energy reconstruction
            bIsGamma = fCuts->isGamma( i, false, fIsOn ) && fCuts->applyEnergyReconstructionQualityCuts( fRunPara->fEnergyReconstructionMethod );
            
            // fill on/offstereo maps and direction cut
            i_theta2 = -99;
            bDirectionCuts = fMap->fill( fIsOn, i_xderot, i_yderot, fCuts->getTheta2Cut_max( iErec ),
                                         fDataRun->runNumber, bIsGamma, i_theta2 );
            bDirectionCuts = fMapUC->fill( fIsOn, i_xderot, i_yderot, fCuts->getTheta2Cut_max( iErec ),
                                           fDataRun->runNumber, bIsGamma, i_theta2 );
                                           
                                           
            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // following histograms (theta2, mscw, mscl, core position, etc.)  assume source at given target position
            
            // theta2 ---
            fMap->calculateTheta2( fIsOn, i_xderot, i_yderot );
            // theta2 values for debugging
            for( unsigned int dex = 0; dex < 25; dex++ )
            {
                fDataRun->theta2_All[dex] = fMap->getTheta2_All()[dex];
            }
            
            for( unsigned int t = 0; t < fMap->getTheta2_length(); t++ )
            {
                fHisto[fHisCounter]->htheta2->Fill( fMap->getTheta2()[t], fMap->getTheta2_weigth()[t] );
            }
            // fill theta for tree with selected events
            if( fMap->getTheta2_length() > 0 )
            {
                //	fDataRun->theta2 = fMap->getTheta2()[0];
            }
            else
            {
                //	fDataRun->theta2 = -1;
            }
            
            if( fDebugCuts )
            {
                cout << "DEBUGCUTS " << i << " / " << fDataRun->eventNumber << ": Direction " << bDirectionCuts;
                cout << ", IsGamma " << bIsGamma << endl;
            }
            
            /////////////////////////////////////////////////////////
            // histograms after direction cuts only
            //
            if( bDirectionCuts )
            {
            
                /////////////////////////////////////////////
                // filling effective areas for each time bin
                
                double time_of_previous_EVENT = time_of_EVENT;
                time_of_EVENT = ( ( double )fDataRun->Time - f_t_in_s_min[irun] );
                double index_time_bin_PREVIOUS_EVENT = index_time_bin_NOW;
                index_time_bin_NOW = fHisto[fHisCounter]->hRealDuration1DtimeBinned->FindFixBin( time_of_EVENT );
                if( time_of_previous_EVENT > time_of_EVENT )
                {
                    cout << "VStereoAnalysis::fillHistograms: filling of time binned effective areas";
                    cout << ", error, events are not ordered chronologically" << endl;
                }
                
                // get 1 / effective area
                iEnergyWeighting = fEnergy.getEffectiveArea( iErec, fDataRun->getZe(),
                                   iDirectionOffset, iPedVar_temp,
                                   fRunPara->fEnergyReconstructionSpectralIndex, true );
                                   
                if( index_time_bin_PREVIOUS_EVENT != index_time_bin_NOW )
                {
                    //--- we just got into a new time bin
                    // getting the effective area for the time bin we just left
                    fEnergy.setTimeBinnedMeanEffectiveArea( fHisto[fHisCounter]->hRealDuration1DtimeBinned->GetBinCenter( index_time_bin_PREVIOUS_EVENT ) );
                }
                
                
                // mean width/length/distance histograms
                fHisto[fHisCounter]->hmean_width->Fill( fCuts->getMeanImageWidth() );
                fHisto[fHisCounter]->hmean_length->Fill( fCuts->getMeanImageLength() );
                fHisto[fHisCounter]->hmean_dist->Fill( fCuts->getMeanImageDistance() );
                if( fDataRun->MSCW > -50. )
                {
                    fHisto[fHisCounter]->hmscw->Fill( fDataRun->MSCW );
                }
                if( fDataRun->MSCL > -50. )
                {
                    fHisto[fHisCounter]->hmscl->Fill( fDataRun->MSCL );
                }
                if( fDataRun->MSCW > -50. && fDataRun->MSCL > -50. )
                {
                    fHisto[fHisCounter]->hmsc->Fill( fDataRun->MSCW, fDataRun->MSCL );
                }
                // probability threshold cuts
                if( fCuts->getProbabilityCut_Selector() > 0. )
                {
                    fHisto[fHisCounter]->hrf->Fill( fCuts->getProbabilityCut_Selector() );
                }
                // mean emission height histograms
                if( fDataRun->EmissionHeight > 0. )
                {
                    fHisto[fHisCounter]->hemiss->Fill( fDataRun->EmissionHeight );
                    fHisto[fHisCounter]->hemissC2->Fill( fDataRun->EmissionHeightChi2 );
                }
                // chi2 of energy reconstruction
                if( iErecChi2 > 0. )
                {
                    fHisto[fHisCounter]->herecChi2->Fill( iErecChi2 );
                }
            }
            // fill a tree with the selected events (possibly after direction cut only)
            if( fRunPara->fWriteEventTree > 0 )
            {
                if( fRunPara->fWriteEventTree == 1 || ( fRunPara->fWriteEventTree == 2 && bDirectionCuts ) )
                {
                    fill_TreeWithSelectedEvents( fDataRun, i_xderot, i_yderot, i_theta2, bDirectionCuts );
                }
            }
            
            // fill a tree with current event for ctools converter
            if( fIsOn && bIsGamma && fRunPara->fWriteEventTreeForCtools )
            {
                fill_TreeWithEventsForCtools( fDataRun, i_xderot, i_yderot );
            }
            /////////////////////////////////////////////////////////
            // histograms after gamma and energy reconstruction cuts
            if( bIsGamma )
            {
                // solid angle of this bin
                double i_ymax = fHisto[fHisCounter]->herecCounts2D_vs_distance->GetYaxis()->GetBinUpEdge(
                                    fHisto[fHisCounter]->herecCounts2D_vs_distance->GetYaxis()->FindBin( iDirectionOffset ) );
                double i_ymin = fHisto[fHisCounter]->herecCounts2D_vs_distance->GetYaxis()->GetBinLowEdge(
                                    fHisto[fHisCounter]->herecCounts2D_vs_distance->GetYaxis()->FindBin( iDirectionOffset ) );
                double iSoli = 2. * TMath::Pi() * ( 1. - cos( i_ymax * TMath::Pi() / 180. ) );
                iSoli       -= 2. * TMath::Pi() * ( 1. - cos( i_ymin * TMath::Pi() / 180. ) );
                // number of events as expected in a theta2 circle at the given offset
                double iWeight = fCuts->getTheta2Cut_max( iErec );
                if( iWeight > 0. )
                {
                    iWeight = 2. * TMath::Pi() * ( 1. - cos( sqrt( iWeight )  * TMath::Pi() / 180. ) );
                    iWeight /= iSoli;
                    fHisto[fHisCounter]->herecCounts2D_vs_distance->Fill( log10( iErec ), iDirectionOffset, iWeight );
                }
            }
            
            /////////////////////////////////////////////////////////
            // histograms after all cuts ( shape and direction cuts )
            //
            if( bIsGamma && bDirectionCuts )
            {
                // image and trigger pattern
                fHisto[fHisCounter]->hTriggerPatternAfterCuts->Fill( fDataRun->LTrig );
                fHisto[fHisCounter]->hImagePatternAfterCuts->Fill( fDataRun->getImgSel() );
                // make core plots
                fHisto[fHisCounter]->hcore->Fill( fDataRun->getXcore_M(), fDataRun->getYcore_M() );
                // ##################################
                // spectral energy reconstruction
                
                // fill energy histograms: require a valid effective area value
                if( iEnergyWeighting > 0. )
                {
                    // energy histogram (counts per bin)
                    fHisto[fHisCounter]->herecCounts->Fill( log10( iErec ) );
                    fHisto[fHisCounter]->herecCounts2DtimeBinned->Fill( log10( iErec ), ( ( double )fDataRun->Time - f_t_in_s_min[irun] ) );
                    fHisto[fHisCounter]->hLinerecCounts->Fill( iErec );
                    fHisto[fHisCounter]->hLinerecCounts2DtimeBinned->Fill( iErec , ( ( double )fDataRun->Time - f_t_in_s_min[irun] ) );
                    fHisto[fHisCounter]->herecWeights->Fill( log10( iErec ), log10( 1. / iEnergyWeighting ) );
                    fHisto[fHisCounter]->hLinerecWeights->Fill( iErec, log10( 1. / iEnergyWeighting ) );
                }
                
                // mean azimuth and elevation
                // (get first running telescope)
                fMeanElevation += fDataRun->ArrayPointing_Elevation;
                fMeanAzimuth    = VSkyCoordinatesUtilities::addToMeanAzimuth( fMeanAzimuth, fDataRun->ArrayPointing_Azimuth );
                fNMeanElevation++;
            }
            // event counter
            if( bIsGamma && bDirectionCuts )
            {
                i_count++;
                fTimeMask->countOn( i_UTC );      // keep track of gamma ON counts for rate plots
            }
        }
    }
    // END: loop over all entries/events in the data tree
    /////////////////////////////////////////////////////////////////////
    
    // filling the effective area for last time bin
    // fill energy histograms: require a valid effective area value
    if( iEnergyWeighting > 0. )
    {
        fEnergy.setTimeBinnedMeanEffectiveArea( fHisto[fHisCounter]->hRealDuration1DtimeBinned->GetBinCenter( index_time_bin_NOW ) );
    }
    // get dead time
    if( fHisCounter == 0 )
    {
        fDeadTimeStorage = fDeadTime[fHisCounter]->calculateDeadTime();
    }
    else
    {
        fDeadTime[fHisCounter]->calculateDeadTime();
    }
    fDeadTime[fHisCounter]->checkStatus();
    fDeadTime[fHisCounter]->printDeadTime();
    
    // filling the histo with the duration of the time bin
    // looping over the mask seconds
    for( unsigned int i_s = 0 ; i_s < fTimeMask->getMaskSize() ; i_s++ )
    {
        if( fTimeMask->getMask()[i_s] )
        {
            // dead time is taken into account for each second
            double dead_time_fraction = fDeadTime[fHisCounter]->getDeadTimeFraction( ( double )i_s + 0.5, fRunPara->fDeadTimeCalculationMethod );
            fHisto[fHisCounter]->hRealDuration1DtimeBinned->Fill( i_s, 1 - dead_time_fraction );
            fHisto[fHisCounter]->hDuration1DtimeBinned->Fill( i_s );
        }
    }
    
    // fill rate vectors
    fTimeMask->getIntervalRates( iRateCounts, iRateTime, iRateTimeIntervall, fRunPara->fTimeIntervall );
    fRateCounts[fHisCounter] = iRateCounts;
    fRateTime[fHisCounter] = iRateTime;
    fRateTimeIntervall[fHisCounter] = iRateTimeIntervall;
    
    // finalize sky maps
    fMap->finalize( fIsOn, fCuts->getProbabilityCutAlpha( fIsOn ) );
    fMapUC->finalize( fIsOn, fCuts->getProbabilityCutAlpha( fIsOn ) );
    
    fTotCount += i_count;
    
    // calculate mean elevation
    if( fNMeanElevation > 0. )
    {
        fMeanAzimuth   /= fNMeanElevation;
        fMeanElevation /= fNMeanElevation;
    }
    
    // get mean effective area
    gMeanEffectiveArea = ( TGraphAsymmErrors* )fEnergy.getMeanEffectiveArea();
    if( gMeanEffectiveArea )
    {
        gMeanEffectiveArea = ( TGraphAsymmErrors* )gMeanEffectiveArea->Clone();
    }
    // get mean energy systematic histogram (needed possibly for energy threshold determination)
    gMeanEsys_MC = ( TGraphErrors* )fEnergy.getMeanSystematicErrorHistogram();
    if( gMeanEsys_MC )
    {
        gMeanEsys_MC = ( TGraphErrors* )gMeanEsys_MC->Clone();
    }
    
    if ( fRunPara->fLikelihoodAnalysis && fIsOn )
    {
      cout << "\t\tVStereoAnalysis::fillHistograms Getting Effective Area MC" << endl;
      gMeanEffectiveAreaMC = ( TGraphAsymmErrors* )fEnergy.getMeanEffectiveAreaMC();
      if( gMeanEffectiveAreaMC )
      {
          gMeanEffectiveAreaMC = ( TGraphAsymmErrors* )gMeanEffectiveAreaMC->Clone();
      }

      cout << "\t\tVStereoAnalysis::fillHistograms Getting Response Matrix" << endl;
      hMeanResponseMatrix = (TH2D*)fEnergy.getMeanResponseMatrix();
      cout << "\t\tVStereoAnalysis::fillHistograms Got Response Matrix" << endl;
      cout << "\t\tVStereoAnalysis::fillHistograms Cloning Response Matrix " << hMeanResponseMatrix << endl;

      if( hMeanResponseMatrix )
      {
          hMeanResponseMatrix = ( TH2D* )hMeanResponseMatrix->Clone();
      }
      cout << "\t\tVStereoAnalysis::fillHistograms Cloned Response Matrix" << endl;

    }
    // get mean effective area for TIME BINs
    gTimeBinnedMeanEffectiveArea = ( TGraph2DErrors* )fEnergy.getTimeBinnedMeanEffectiveArea()->Clone();
    // get mean run times after time cuts
    fRunExposure[irun] = fTimeMask->getEffectiveDuration();
    fRunMJD[irun] = fTimeMask->getMeanUTC_Mask();
    fTimeMask->printMask( 100, kTRUE );
    fTimeMask->printMeanTime( kTRUE );
    
    return i_count;
}


/*
 *
 * write created histograms to the appropriate directories and tidy up
 *
 */
void VStereoAnalysis::writeHistograms( bool bOn )
{
    if( fDebug )
    {
        cout << "DEBUG void VStereoAnalysis::writeHistograms()" << endl;
    }
    if( fHisCounter < 0 )
    {
        fHistoTot->writeHistograms();
    }
    else
    {
        if( fCuts )
        {
            fCuts->Write();
        }
        fTimeMask->writeObjects();
        if( bOn )
        {
        
            // write the VTimeMask object to the root file
            VTimeMask* iTimeMask = ( VTimeMask* )fTimeMask->Clone();
            iTimeMask->Write( "vtimemask" );
            delete iTimeMask;
            
            // write the VDeadTime object to the root file
            VDeadTime* iDeadTime = ( VDeadTime* )( fDeadTime[fHisCounter]->Clone() ) ;
            iDeadTime->Write( "vdeadtime" ) ;
            delete iDeadTime;
        }
        fHisto[fHisCounter]->writeHistograms();
        
        // need to grab fScalarDeadTimeFrac while fDeadTime histograms are intact,
        // fScalarDeadTimeFrac is needed in  save_TreeWithEventsForCtools()
        fRunPara->fScalarDeadTimeFrac = fDeadTime[fHisCounter]->getDeadTimeFraction( fTimeMask->getMask(), fRunPara->fDeadTimeCalculationMethod );
        
        fDeadTime[fHisCounter]->writeHistograms();
        // copy effective areas and radial acceptance to anasum output file
        if( bOn )
        {
            if( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile != "IGNOREEFFECTIVEAREA" )
            {
                fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", gMeanEffectiveArea );
                fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", gTimeBinnedMeanEffectiveArea );
                fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", gMeanEsys_MC );
            }
            if (fRunPara->fLikelihoodAnalysis)
            {
              cout << "\t\tVStereoAnalysis::writeHistograms Writing histograms" << endl;
              fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", gMeanEffectiveAreaMC );
              fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", hMeanResponseMatrix );

            }
            if( fRunPara->fRunList[fHisCounter].fAcceptanceFile.size() > 0
                    && fRunPara->fRunList[fHisCounter].fAcceptanceFile != "IGNOREACCEPTANCE" )
            {
                fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fAcceptanceFile, "RadialAcceptances", 0 );
            }
        }
        else
        {
            char hname[1000];
            if( gMeanEffectiveArea )
            {
                sprintf( hname, "%s_off", gMeanEffectiveArea->GetName() );
                gMeanEffectiveArea->SetName( hname );
            }
            if( gTimeBinnedMeanEffectiveArea )
            {
                sprintf( hname, "%s_off", gTimeBinnedMeanEffectiveArea->GetName() );
                gTimeBinnedMeanEffectiveArea->SetName( hname );
            }
            if( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile != "IGNOREEFFECTIVEAREA" )
            {
                fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", gMeanEffectiveArea );
                fHisto[fHisCounter]->writeObjects( fRunPara->fRunList[fHisCounter].fEffectiveAreaFile, "EffectiveAreas", gTimeBinnedMeanEffectiveArea );
            }
        }
        if( fTreeSelectedEvents )
        {
            fTreeSelectedEvents->AutoSave();
        }
        
        if( fTreeWithEventsForCtools && fIsOn && fRunPara->fWriteEventTreeForCtools )  // WRITEEVENTTREEFORCTOOLS
        {
            save_TreeWithEventsForCtools() ;
        }
    }
}


void VStereoAnalysis::writeDebugHistograms()
{
    if( fDebug )
    {
        cout << "DEBUG void VStereoAnalysis::writeDebugHistograms()" << endl;
    }
    
    TDirectory* iDir = gDirectory;
    
    if( iDir->mkdir( "debug" )->cd() )
    {
        if( fMap && fMap->getAux_hisList() )
        {
            fMap->getAux_hisList()->Write();
            fMap->getAux_hisList()->Delete();
        }
    }
    
    iDir->cd();
}

/*

    divide on by off alpha histograms

    this should only be called for OFF stereo analysis

*/
void VStereoAnalysis::scaleAlpha( TH2D* halpha_on, bool bUC )
{
    if( fIsOn )
    {
        cout << "VStereoAnalysis::scaleAlpha() error: this function should only be called for OFF stereo analysis" << endl;
        cout << "(this must be a coding error, please report)" << endl;
        exit( EXIT_FAILURE );
    }
    TH2D* halpha_off = 0;
    TH2D* hmap_alphaNorm = 0;
    // uncorrelated maps
    if( bUC )
    {
        if( fHisCounter < 0 )
        {
            halpha_off = fHistoTot->hmap_alphaUC;
            hmap_alphaNorm = fHistoTot->hmap_alphaNormUC;
        }
        else
        {
            halpha_off = fHisto[fHisCounter]->hmap_alphaUC;
            hmap_alphaNorm = fHisto[fHisCounter]->hmap_alphaNormUC;
        }
    }
    // correlated maps
    else
    {
        if( fHisCounter < 0 )
        {
            halpha_off = fHistoTot->hmap_alpha;
            hmap_alphaNorm = fHistoTot->hmap_alphaNorm;
        }
        else
        {
            halpha_off = fHisto[fHisCounter]->hmap_alpha;
            hmap_alphaNorm = fHisto[fHisCounter]->hmap_alphaNorm;
        }
    }
    if( !halpha_off || !hmap_alphaNorm )
    {
        cout << "VStereoAnalysis::scaleAlpha: fatal error, cannot find histograms ";
        cout << halpha_off << "\t" << hmap_alphaNorm << endl;
        exit( EXIT_FAILURE );
    }
    
    // halpha_on: on alpha histogram
    // halpha_off: off alpha histogram
    // hmap_alphaNorm: alpha histogram used in significance calculations (alphaNorm)
    // (this is slightly different for average alpha calculation)
    for( int i = 1; i <= halpha_off->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= halpha_off->GetNbinsY(); j++ )
        {
            if( halpha_off->GetBinContent( i, j ) > 0. )
            {
                // this one is used for the sky maps
                hmap_alphaNorm->SetBinContent( i, j, halpha_on->GetBinContent( i, j ) / halpha_off->GetBinContent( i, j ) );
            }
            else
            {
                hmap_alphaNorm->SetBinContent( i, j, 0. );
            }
        }
    }
}


/*!
 *   combine histograms from all runs
 *
 *   this function shall be called at the end of the analysis
 */
double VStereoAnalysis::combineHistograms()
{
    unsigned int n_histo = fHisto.size();
    
    TDirectory* iDir = gDirectory;
    fDirTot->cd();
    fHistoTot->defineHistograms();
    
    // list of trees with selected events
    iDir->cd();
    
    ///////////////////////////////////////////////////
    // loop over all runs (= all available histograms = n_histo)
    for( unsigned h = 0; h < n_histo; h++ )
    {
        fDirTotRun[h]->cd();
        // read in sky plots from disk
        fHisto[h]->readSkyPlots();
        
        /////////////////////////////
        // UNCORRELATED PLOTS
        int nxbin = fHistoTot->hmap_stereoUC->GetNbinsX();
        int nybin = fHistoTot->hmap_stereoUC->GetNbinsY();
        for( int i = 1; i <= nxbin; i++ )
        {
            for( int j = 1; j <= nybin; j++ )
            {
                // calculate average normalization (alpha) factor
                if( fHisto[h]->hmap_alphaUC && fHisto[h]->hmap_alphaUC->GetBinContent( i, j ) > 0.
                        && fHisto[h]->h_combine_map_alpha_offUC && fHisto[h]->h_combine_map_alpha_offUC->GetBinContent( i, j ) > 0. )
                {
                    fHistoTot->hmap_stereoUC->SetBinContent( i, j, fHisto[h]->hmap_stereoUC->GetBinContent( i, j )
                            + fHistoTot->hmap_stereoUC->GetBinContent( i, j ) );
                    // calculate average alpha
                    if( fHisto[h]->h_combine_map_stereo_onUC && fHisto[h]->h_combine_map_stereo_offUC )
                    {
                        float alphaUC = 0.;
                        if( fIsOn && fHisto[h]->h_combine_map_alpha_offUC->GetBinContent( i, j ) != -1. )
                        {
                            alphaUC = fHisto[h]->h_combine_map_alpha_offUC->GetBinContent( i, j )
                                      / ( 1. + fHisto[h]->h_combine_map_alpha_offUC->GetBinContent( i, j ) );
                        }
                        else if( !fIsOn )
                        {
                            alphaUC = 1. / ( 1. + fHisto[h]->h_combine_map_alpha_offUC->GetBinContent( i, j ) );
                        }
                        alphaUC *= ( fHisto[h]->h_combine_map_stereo_onUC->GetBinContent( i, j )
                                     + fHisto[h]->h_combine_map_stereo_offUC->GetBinContent( i, j ) );
                        fHistoTot->hmap_alphaUC->SetBinContent( i, j, fHistoTot->hmap_alphaUC->GetBinContent( i, j ) + alphaUC );
                    }
                }
            }
        }
        //////////////////////////////
        // CORRELATED PLOTS
        nxbin = fHistoTot->hmap_stereo->GetNbinsX();
        nybin = fHistoTot->hmap_stereo->GetNbinsY();
        for( int i = 1; i <= nxbin; i++ )
        {
            for( int j = 1; j <= nybin; j++ )
            {
                // calculate average normalization (alpha) factor
                if( fHisto[h]->hmap_alpha && fHisto[h]->hmap_alpha->GetBinContent( i, j ) > 0.
                        && fHisto[h]->h_combine_map_alpha_off && fHisto[h]->h_combine_map_alpha_off->GetBinContent( i, j ) > 0. )
                {
                    fHistoTot->hmap_stereo->SetBinContent( i, j, fHisto[h]->hmap_stereo->GetBinContent( i, j )
                                                           + fHistoTot->hmap_stereo->GetBinContent( i, j ) );
                    // calculate average alpha
                    if( fHisto[h]->h_combine_map_stereo_on && fHisto[h]->h_combine_map_stereo_off )
                    {
                        float alpha = 0.;
                        if( fIsOn && fHisto[h]->h_combine_map_alpha_off->GetBinContent( i, j ) != -1. )
                        {
                            alpha = fHisto[h]->h_combine_map_alpha_off->GetBinContent( i, j )
                                    / ( 1. + fHisto[h]->h_combine_map_alpha_off->GetBinContent( i, j ) );
                        }
                        else if( !fIsOn )
                        {
                            alpha = 1. / ( 1. + fHisto[h]->h_combine_map_alpha_off->GetBinContent( i, j ) );
                        }
                        alpha *= ( fHisto[h]->h_combine_map_stereo_on->GetBinContent( i, j )
                                   + fHisto[h]->h_combine_map_stereo_off->GetBinContent( i, j ) );
                        fHistoTot->hmap_alpha->SetBinContent( i, j, fHistoTot->hmap_alpha->GetBinContent( i, j ) + alpha );
                    }
                }
            }
        }
        fHisto[h]->deleteSkyPlots();
        iDir->cd();
    }  // (end loop over all histograms)
    
    //////////////////////////////////////
    // errors in sky maps (counting error)
    for( int i = 1; i <= fHistoTot->hmap_stereoUC->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= fHistoTot->hmap_stereoUC->GetNbinsY(); j++ )
        {
            if( fHistoTot->hmap_stereoUC->GetBinContent( i, j ) > 0 )
            {
                fHistoTot->hmap_stereoUC->SetBinError( i, j, sqrt( fHistoTot->hmap_stereoUC->GetBinContent( i, j ) ) );
            }
        }
    }
    for( int i = 1; i <= fHistoTot->hmap_stereo->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= fHistoTot->hmap_stereo->GetNbinsY(); j++ )
        {
            if( fHistoTot->hmap_stereo->GetBinContent( i, j ) > 0 )
            {
                fHistoTot->hmap_stereo->SetBinError( i, j, sqrt( fHistoTot->hmap_stereo->GetBinContent( i, j ) ) );
            }
        }
    }
    
    //////////////////////////////////////
    // combine parameter (1D) histograms
    for( unsigned int h = 0; h < n_histo; h++ )
    {
        iDir->cd();
        fDirTotRun[h]->cd();
        fHisto[h]->readParameterHistograms();
        TIter next( fHistoTot->hListParameterHistograms );
        TIter nexth( fHisto[h]->hListParameterHistograms );
        while( TH1* h1 = ( TH1* )next() )
        {
            TH1* h2 = ( TH1* )nexth();
            if( !h1 || !h2 )
            {
                continue;
            }
            
            string iTemp = h1->GetName();
            if( iTemp.find( "2D" ) != string::npos )
            {
                continue;
            }
            if( iTemp.find( "Duration1DtimeBinned" ) != string::npos )
            {
                continue;
            }
            h1->Add( h2 );
            
        }
        fHisto[h]->deleteParameterHistograms();
    }
    iDir->cd();
    
    // combine rate vectors (in time intervals)
    for( unsigned int h = 0; h < n_histo; h++ )
    {
        for( unsigned int i = 0; i <  fRateCounts[h].size(); i++ )
        {
            fRateCountsTot.push_back( fRateCounts[h][i] );
        }
        for( unsigned int i = 0; i <  fRateTime[h].size(); i++ )
        {
            fRateTimeTot.push_back( fRateTime[h][i] );
        }
        for( unsigned int i = 0; i <  fRateTimeIntervall[h].size(); i++ )
        {
            fRateTimeIntervallTot.push_back( fRateTimeIntervall[h][i] );
        }
    }
    
    iDir->cd();
    return fTotCount;
}

TH1D* VStereoAnalysis::getMeanSignalBackgroundAreaRatio()
{
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_MeanSignalBackgroundAreaRatio;
    }
    
    return fHisto[fHisCounter]->hmap_MeanSignalBackgroundAreaRatio;
}

TH1D* VStereoAnalysis::getMeanSignalBackgroundAreaRatioUC()
{
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_MeanSignalBackgroundAreaRatioUC;
    }
    
    return fHisto[fHisCounter]->hmap_MeanSignalBackgroundAreaRatioUC;
}


TH2D* VStereoAnalysis::getAlpha()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getAlpha() " << fHisCounter << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_alpha;
    }
    
    return fHisto[fHisCounter]->hmap_alpha;
}


TH2D* VStereoAnalysis::getAlphaUC()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getAlphaUC() " << fHisCounter << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_alphaUC;
    }
    
    return fHisto[fHisCounter]->hmap_alphaUC;
}


TH2D* VStereoAnalysis::getAlphaNorm()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getAlphaNorm() " << fHisCounter << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_alphaNorm;
    }
    
    return fHisto[fHisCounter]->hmap_alphaNorm;
}


TH2D* VStereoAnalysis::getAlphaNormUC()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getAlphaNormUC() " << fHisCounter << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_alphaNormUC;
    }
    
    return fHisto[fHisCounter]->hmap_alphaNormUC;
}


TList* VStereoAnalysis::getHisList()
{
    if( fHisCounter < 0 )
    {
        return fHistoTot->hisList;
    }
    
    return fHisto[fHisCounter]->hisList;
}


TList* VStereoAnalysis::getSkyHistograms( bool bUC )
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getSkyHistograms() " << fHisCounter << "\t" << bUC << endl;
    }
    
    // uncorrelated plot
    if( bUC )
    {
        if( fHisCounter < 0 )
        {
            return fHistoTot->hListSkyMapsUC;
        }
        return fHisto[fHisCounter]->hListSkyMapsUC;
    }
    // correlated plots
    if( fHisCounter < 0 )
    {
        return fHistoTot->hListSkyMaps;
    }
    return fHisto[fHisCounter]->hListSkyMaps;
}


TList* VStereoAnalysis::getParameterHistograms()
{
    if( fHisCounter < 0 )
    {
        return fHistoTot->hListParameterHistograms;
    }
    
    return fHisto[fHisCounter]->hListParameterHistograms;
}


TH2D* VStereoAnalysis::getStereoSkyMapUC()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getStereoSkyMapUC()" << "\t" << fHisCounter << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_stereoUC;
    }
    
    return fHisto[fHisCounter]->hmap_stereoUC;
}


TH2D* VStereoAnalysis::getStereoSkyMap()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getStereoSkyMap()" << "\t" << fHisCounter << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return fHistoTot->hmap_stereo;
    }
    
    return fHisto[fHisCounter]->hmap_stereo;
}


double VStereoAnalysis::getDeadTimeFraction()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::getDeadTimeFraction()" << endl;
    }
    
    if( fHisCounter < 0 )
    {
        return 0.;
    }
    
    if( fHisCounter < ( int )fDeadTime.size() )
    {
        // dead time depending on time mask
        if( fTimeMask && fTimeMask->getMask().size() > 0 )
        {
            return fDeadTime[fHisCounter]->getDeadTimeFraction( fTimeMask->getMask(), fRunPara->fDeadTimeCalculationMethod );
        }
        return fDeadTime[fHisCounter]->getDeadTimeFraction( -1, fRunPara->fDeadTimeCalculationMethod );
    }
    
    return 0.;
}


/*
  this function is called for each run

  - set targets
  - set sky map centres
  - set area to calculate 1D histograms and energy spectra
  - set exclusion regions for background calculation

  Note that throughout the analysis and mapfilling, the coordinate system is J2000

*/
void VStereoAnalysis::defineAstroSource()
{
    if( fDebug )
    {
        cout << "VStereoAnalysis::defineAstroSource()" << endl;
    }
    
    //////////////////////////////////////
    // loop over all runs in runlist
    for( unsigned int i = 0; i < fRunPara->fRunList.size(); i++ )
    {
        if( fIsOn )
        {
            cout << endl << "-----------------------------------------------------------------------" << endl;
            cout << "Defining targets for on run " << fRunPara->fRunList[i].fRunOn << endl;
        }
        else
        {
            cout << endl << "-----------------------------------------------------------------------" << endl;
            cout << "Defining targets and exclusion regions for off run " << fRunPara->fRunList[i].fRunOff << endl;
        }
        
        /////////////////////////////////////////////////////////
        // check source coordinates
        // (this is the target of observation)
        /////////////////////////////////////////////////////////
        if( fRunPara->fRunList[i].fTargetDecJ2000 < -89.99 )
        {
            cout << "ERROR in VStereoAnalysis::defineAstroSource: invalid target " << fRunPara->fRunList[i].fTarget << endl;
            cout << "\t run " << fRunPara->fRunList[i].fRunOn << "\t" << fRunPara->fRunList[i].fTarget;
            cout << fRunPara->fRunList[i].fTargetDecJ2000 << "\t" << fRunPara->fRunList[i].fTargetShiftDecJ2000 << endl;
            exit( EXIT_FAILURE );
        }
        
        bool iUserSetSkyMapCentre = false;
        /////////////////////////////////////////////////////////
        // from runparameter file: set the sky map centre as xy offset [deg]
        if( TMath::Abs( fRunPara->fSkyMapCentreNorth ) > 1.e-8 || TMath::Abs( fRunPara->fSkyMapCentreWest ) > 1.e-8 )
        {
            fRunPara->fRunList[i].fSkyMapCentreWest  = fRunPara->fSkyMapCentreWest;
            fRunPara->fRunList[i].fSkyMapCentreNorth = fRunPara->fSkyMapCentreNorth;
            double i_decDiff =  0.;
            double i_raDiff = 0.;
            VSkyCoordinatesUtilities::getWobbleOffset_in_RADec( fRunPara->fRunList[i].fSkyMapCentreNorth,
                    fRunPara->fRunList[i].fSkyMapCentreWest,
                    fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                    i_decDiff, i_raDiff );
            fRunPara->fRunList[i].fSkyMapCentreRAJ2000  = fRunPara->fRunList[i].fTargetRAJ2000 + i_raDiff;
            fRunPara->fRunList[i].fSkyMapCentreDecJ2000 = fRunPara->fRunList[i].fTargetDecJ2000 + i_decDiff;
            if( TMath::Abs( i_raDiff ) > 1.e-2 || TMath::Abs( i_decDiff ) > 1.e-2 )
            {
                iUserSetSkyMapCentre = true;
            }
        }
        // from runparameter file: set the sky map centre in J2000
        // (this is in almost all analysis the usual/default case)
        else if( TMath::Abs( fRunPara->fSkyMapCentreRAJ2000 ) > 1.e-8 )
        {
            fRunPara->fRunList[i].fSkyMapCentreRAJ2000  = fRunPara->fSkyMapCentreRAJ2000;
            fRunPara->fRunList[i].fSkyMapCentreDecJ2000 = fRunPara->fSkyMapCentreDecJ2000;
            fRunPara->fRunList[i].fSkyMapCentreWest =
                VSkyCoordinatesUtilities::getTargetShiftWest( fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                        fRunPara->fSkyMapCentreRAJ2000, fRunPara->fSkyMapCentreDecJ2000 ) * -1.;
            fRunPara->fRunList[i].fSkyMapCentreNorth =
                VSkyCoordinatesUtilities::getTargetShiftNorth( fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                        fRunPara->fSkyMapCentreRAJ2000, fRunPara->fSkyMapCentreDecJ2000 );
            if( TMath::Abs( fRunPara->fRunList[i].fSkyMapCentreWest ) < 1.e-4 )
            {
                fRunPara->fRunList[i].fSkyMapCentreWest = 0.;
            }
            if( TMath::Abs( fRunPara->fRunList[i].fSkyMapCentreNorth ) < 1.e-4 )
            {
                fRunPara->fRunList[i].fSkyMapCentreNorth = 0.;
            }
            if( TMath::Abs( fRunPara->fRunList[i].fSkyMapCentreWest ) > 1.e-2
                    ||  TMath::Abs( fRunPara->fRunList[i].fSkyMapCentreNorth ) > 1.e-2 )
            {
                iUserSetSkyMapCentre = true;
            }
        }
        // if not set in runparameter file: set to target direction
        else
        {
            fRunPara->fRunList[i].fSkyMapCentreRAJ2000  = fRunPara->fRunList[i].fTargetRAJ2000;
            fRunPara->fRunList[i].fSkyMapCentreDecJ2000 = fRunPara->fRunList[i].fTargetDecJ2000;
            fRunPara->fSkyMapCentreRAJ2000              = fRunPara->fRunList[i].fSkyMapCentreRAJ2000;
            fRunPara->fSkyMapCentreDecJ2000             = fRunPara->fRunList[i].fSkyMapCentreDecJ2000;
        }
        
        /////////////////////////////////////////////////////////
        // from runparameter file: set and get target shifts
        // (calculated relative to sky map centre)
        // (this is the position where all 1D histograms (theta2, energy spectra, etc) are calculated)
        if( fIsOn )
        {
            if( TMath::Abs( fRunPara->fTargetShiftDecJ2000 ) > 1.e-8 || TMath::Abs( fRunPara->fTargetShiftRAJ2000 ) > 1.e-8 )
            {
                fRunPara->fRunList[i].fTargetShiftWest = VSkyCoordinatesUtilities::getTargetShiftWest( fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                        fRunPara->fTargetShiftRAJ2000, fRunPara->fTargetShiftDecJ2000 );
                fRunPara->fRunList[i].fTargetShiftNorth = -1.*VSkyCoordinatesUtilities::getTargetShiftNorth( fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                        fRunPara->fTargetShiftRAJ2000, fRunPara->fTargetShiftDecJ2000 );
                        
                fRunPara->fRunList[i].fTargetShiftWest  += fRunPara->fRunList[i].fSkyMapCentreWest;
                fRunPara->fRunList[i].fTargetShiftNorth += fRunPara->fRunList[i].fSkyMapCentreNorth;
                if( TMath::Abs( fRunPara->fRunList[i].fTargetShiftWest ) < 1.e-4 )
                {
                    fRunPara->fRunList[i].fTargetShiftWest = 0.;
                }
                if( TMath::Abs( fRunPara->fRunList[i].fTargetShiftNorth ) < 1.e-4 )
                {
                    fRunPara->fRunList[i].fTargetShiftNorth = 0.;
                }
            }
            else
            {
                fRunPara->fRunList[i].fTargetShiftWest  = fRunPara->fTargetShiftWest;
                fRunPara->fRunList[i].fTargetShiftNorth = fRunPara->fTargetShiftNorth;
            }
            fRunPara->fRunList[i].fTargetShiftWest *= -1.;
            fRunPara->fTargetShiftWest = fRunPara->fRunList[i].fTargetShiftWest;
            fRunPara->fTargetShiftNorth = fRunPara->fRunList[i].fTargetShiftNorth;
            fRunPara->setTargetShifts( i, fRunPara->fRunList[i].fTargetShiftWest, fRunPara->fRunList[i].fTargetShiftNorth,
                                       fRunPara->fTargetShiftRAJ2000, fRunPara->fTargetShiftDecJ2000 );
        }
        /////////////////////////////////////////////////////////
        // precess target coordinates from J2000 to current epoch
        // (direction of telescope pointing)
        double i_dec = fRunPara->fRunList[i].fTargetDecJ2000 * TMath::DegToRad();
        double i_ra  = fRunPara->fRunList[i].fTargetRAJ2000 * TMath::DegToRad();
        double iMJD = ( double )fRunPara->fRunList[i].fMJDOn;
        if( !fIsOn )
        {
            iMJD = ( double )fRunPara->fRunList[i].fMJDOff;
        }
        // (i_dec and i_ra are in current epoch coordinates in the following, not J2000)
        VSkyCoordinatesUtilities::precessTarget( iMJD, i_ra, i_dec );
        
        // print some information on targeting/pointing to screen
        if( fIsOn )
        {
            cout << "Run " << fRunPara->fRunList[i].fRunOn << " ---------------------------" << endl;
            // set target coordinates into run parameter list
            fRunPara->setTargetRADec_currentEpoch( i, i_ra * TMath::RadToDeg(), i_dec * TMath::RadToDeg() );
            // print target info to screen
            cout << "\tTarget: " << fRunPara->fRunList[i].fTarget << " (ra,dec)=(";
            cout << fRunPara->fRunList[i].fTargetRA << ", " << fRunPara->fRunList[i].fTargetDec << ")";
            cout << " (precessed, MJD=" << iMJD << "), ";
            cout << "(ra,dec (J2000)) = (" << fRunPara->fRunList[i].fTargetRAJ2000 << ", " << fRunPara->fRunList[i].fTargetDecJ2000 << ")";
            if( TMath::Abs( fRunPara->fRunList[i].fPairOffset ) > 1.e-2 )
            {
                cout << ", pair offset [min]: " << fRunPara->fRunList[i].fPairOffset;
            }
            cout << endl;
        }
        /////////////////////////////////////
        // calculate wobble offsets in J2000
        // (this might be overcomplicated)
        /////////////////////////////////////
        // calculate wobble offset in ra/dec for current epoch
        double i_decDiff = 0.;
        double i_raDiff = 0.;
        VSkyCoordinatesUtilities::getWobbleOffset_in_RADec( fRunPara->fRunList[i].fWobbleNorth, -1.*fRunPara->fRunList[i].fWobbleWest,
                i_dec * TMath::RadToDeg(), i_ra * TMath::RadToDeg(), i_decDiff, i_raDiff );
        if( i_raDiff < -180. )
        {
            i_raDiff += 360.;
        }
        // ra/dec of pointing direction in current epoch
        double i_decWobble = i_dec * TMath::RadToDeg() + i_decDiff;
        double i_raWobble  = i_ra * TMath::RadToDeg()  + i_raDiff;
        // correct for precession (from current epoch to J2000=MJD51544)
        VSkyCoordinatesUtilities::precessTarget( 51544., i_raWobble, i_decWobble, iMJD, true );
        double i_WobbleJ2000_West = VSkyCoordinatesUtilities::getTargetShiftWest( fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                                    i_raWobble, i_decWobble ) * -1.;
        if( TMath::Abs( i_WobbleJ2000_West ) < 1.e-4 )
        {
            i_WobbleJ2000_West = 0.;
        }
        double i_WobbleJ2000_North = VSkyCoordinatesUtilities::getTargetShiftNorth( fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000,
                                     i_raWobble, i_decWobble );
        if( TMath::Abs( i_WobbleJ2000_North ) < 1.e-4 )
        {
            i_WobbleJ2000_North = 0.;
        }
        // modify wobble offsets for centering of sky maps
        fRunPara->fRunList[i].fWobbleNorthMod = i_WobbleJ2000_North - fRunPara->fRunList[i].fSkyMapCentreNorth;
        fRunPara->fRunList[i].fWobbleWestMod  = i_WobbleJ2000_West  - fRunPara->fRunList[i].fSkyMapCentreWest;
        
        // fill run parameter values
        fRunPara->setTargetRADecJ2000( i, fRunPara->fRunList[i].fTargetRAJ2000, fRunPara->fRunList[i].fTargetDecJ2000 );
        fRunPara->setTargetShifts( i, fRunPara->fRunList[i].fTargetShiftWest, fRunPara->fRunList[i].fTargetShiftNorth,
                                   fRunPara->fTargetShiftRAJ2000, fRunPara->fTargetShiftDecJ2000 );
        fRunPara->setSkyMapCentreJ2000( i, fRunPara->fRunList[i].fSkyMapCentreRAJ2000, fRunPara->fRunList[i].fSkyMapCentreDecJ2000 );
        
        ///////////////////////////////////////////////////////////////////
        // some printout
        if( fIsOn )
        {
            cout << "\tWobble offsets (currE): N: " << fRunPara->fRunList[i].fWobbleNorth << " W: " << fRunPara->fRunList[i].fWobbleWest;
            cout << ",  RA " << i_raDiff << ", Dec " << i_decDiff << endl;
            cout << "\tWobble offsets (J2000): N: " << i_WobbleJ2000_North << " W: " << i_WobbleJ2000_West << endl;
            cout << "\tSky maps centred at (ra,dec (J2000)) (";
            cout << fRunPara->fRunList[i].fSkyMapCentreRAJ2000 << ", " << fRunPara->fRunList[i].fSkyMapCentreDecJ2000 << ")" << endl;
            cout << "\tTelescopes pointing to: (ra,dec (J2000)) (" << i_raWobble << ", " << i_decWobble << ")";
            cout << ", N: " << fRunPara->fRunList[i].fWobbleNorthMod << " W: " << fRunPara->fRunList[i].fWobbleWestMod << endl;
            cout << "\t1D-histograms calculated at (x,y): " << fRunPara->fRunList[i].fTargetShiftNorth << ", " << fRunPara->fRunList[i].fTargetShiftWest;
            if( TMath::Abs( fRunPara->fTargetShiftDecJ2000 ) > 1.e-8 &&  TMath::Abs( fRunPara->fTargetShiftRAJ2000 ) > 1.e-8 )
            {
                cout << " (ra,dec (J2000)) " << fRunPara->fTargetShiftRAJ2000 << ", " << fRunPara->fTargetShiftDecJ2000;
            }
            cout << endl;
        }
        
        //////////////////////////////////
        
        // =============================================================
        // define source and tracking class
        fAstro.push_back( new VSkyCoordinates() );
        // get wobble offsets in ra,dec
        if( fIsOn )
        {
            i_dec = fRunPara->fRunList[i].fTargetDecJ2000 * TMath::DegToRad();
            i_ra  = fRunPara->fRunList[i].fTargetRAJ2000 * TMath::DegToRad();
        }
        else
        {
            i_dec = fRunPara->fRunList[i].fOff_TargetDecJ2000 * TMath::DegToRad();
            i_ra  = fRunPara->fRunList[i].fOff_TargetRAJ2000 * TMath::DegToRad();
        }
        // precess target to current epoch
        VSkyCoordinatesUtilities::precessTarget( iMJD, i_ra, i_dec );
        double idec_T = 0.;
        double ira_T = 0.;
        if( fIsOn )
        {
            VSkyCoordinatesUtilities::getWobbledDirection( fRunPara->fRunList[i].fWobbleNorth, fRunPara->fRunList[i].fWobbleWest,
                    i_dec * TMath::RadToDeg(), i_ra * TMath::RadToDeg(), idec_T, ira_T );
        }
        else
        {
            VSkyCoordinatesUtilities::getWobbledDirection( fRunPara->fRunList[i].fOff_WobbleNorth, fRunPara->fRunList[i].fOff_WobbleWest,
                    i_dec * TMath::RadToDeg(), i_ra * TMath::RadToDeg(), idec_T, ira_T );
        }
        // setting telescope coordinates (in current epoch)
        // (ignore pointing errors here (very small impact))
        fAstro.back()->setTelDec_deg( idec_T );
        fAstro.back()->setTelRA_deg( ira_T );
        // set observatory position
        fAstro.back()->setObservatory( fRunPara->getObservatory_Longitude_deg(), fRunPara->getObservatory_Latitude_deg() );
        
        // =============================================================
        // set up star catalogue and exclusion regions
        // (all coordinates in J2000)
        if( !fIsOn )
        {
            double i_decC = fRunPara->fRunList[i].fSkyMapCentreDecJ2000;
            double i_raC = fRunPara->fRunList[i].fSkyMapCentreRAJ2000;
            // on/off or matched run analysis
            if( fRunPara->fRunList[i].fRunOn != fRunPara->fRunList[i].fRunOff )
            {
                if( !iUserSetSkyMapCentre )
                {
                    i_decC = fRunPara->fRunList[i].fOff_TargetDecJ2000;
                    i_raC  = fRunPara->fRunList[i].fOff_TargetRAJ2000;
                }
                else
                {
                    cout << "VStereoAnalysis::defineAstroSource warning: matched run analysis and sky map centre shift are not implemented yet" << endl;
                    cout << "\t the location of the exclusion regions due to bright stars will be wrong" << endl;
                }
            }
            // star catalogue is initialized around the pointing direction of the background run
            // (as stars are excluded from the background estimation in the off maps)
            // assume a much larger search region then the actual sky map
            fAstro.back()->initStarCatalogue( fRunPara->getStarCatalogue(), iMJD,
                                              fRunPara->fSkyMapSizeXmin - fRunPara->getLargestStarExlusionRadius(),
                                              fRunPara->fSkyMapSizeXmax + fRunPara->getLargestStarExlusionRadius(),
                                              fRunPara->fSkyMapSizeYmin - fRunPara->getLargestStarExlusionRadius(),
                                              fRunPara->fSkyMapSizeYmax + fRunPara->getLargestStarExlusionRadius(),
                                              i_raC, i_decC );
            // initialize exclusion regions
            // (even if there is no star catalogue given,
            //  there might be user given exclusion regions)
            fRunPara->initializeExclusionRegions( i, fAstro.back()->getStarCatalogue(),
                                                  i_raC, i_decC,
                                                  i_raWobble, i_decWobble );
        }
    } // end loop over all runs
}
/////////////////////////////////////////////////////////


void VStereoAnalysis::setCuts( VAnaSumRunParameterDataClass iL, int irun )
{
    if( iL.fCutFile != "" )
    {
        // read cuts from effective area root file
        if( iL.fCutFile.find( ".root" ) != string::npos
                && iL.fCutFile.find( "IGNOREEFFECTIVEAREA" ) == string::npos )
        {
            string iEffFile = VUtilities::testFileLocation( iL.fCutFile, "EffectiveAreas", true );
            
            TFile* iF  = new TFile( iEffFile.c_str() );
            if( iF->IsZombie() )
            {
                cout << "VStereoAnalysis::setCuts error opening file to read cuts: " << endl;
                cout << "\t" << iEffFile << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            VGammaHadronCuts* iC = ( VGammaHadronCuts* )iF->Get( "GammaHadronCuts" );
            if( !iC )
            {
                cout << "VStereoAnalysis::setCuts error reading cuts from file: " << endl;
                cout << "\t" << iEffFile << endl;
                cout << "exciting..." << endl;
                exit( EXIT_FAILURE );
            }
            fCuts = iC;
            cout << "Reading gamma/hadron cuts from effective area file " << iEffFile << endl;
            iF->Close();
        }
        // read cuts from text file
        else
        {
            fCuts->setNTel( iL.fMaxTelID );
            fCuts->setInstrumentEpoch( fInstrumentEpoch );
            fCuts->setTelToAnalyze( fTelToAnalyze );
            fCuts->readCuts( iL.fCutFile );
            fCuts->setTheta2Cut( iL.fSourceRadius );
        }
    }
    else
    {
        fCuts->resetCutValues();
    }
    fCuts->initializeCuts( irun );
    fCuts->printCutSummary();
}


vector< double > VStereoAnalysis::getRateCounts()
{
    if( fHisCounter < 0 )
    {
        return fRateCountsTot;
    }
    
    if( fHisCounter < ( int )fRateCounts.size() )
    {
        return fRateCounts[fHisCounter];
    }
    
    // this shouldn't happen
    vector< double > f;
    return f;
}


vector< double > VStereoAnalysis::getRateTime()
{
    if( fHisCounter < 0 )
    {
        return fRateTimeTot;
    }
    
    if( fHisCounter < ( int )fRateTime.size() )
    {
        return fRateTime[fHisCounter];
    }
    
    // this shouldn't happen
    vector< double > f;
    return f;
}


vector< double > VStereoAnalysis::getRateTimeIntervall()
{
    if( fHisCounter < 0 )
    {
        return fRateTimeIntervallTot;
    }
    
    if( fHisCounter < ( int )fRateTimeIntervall.size() )
    {
        return fRateTimeIntervall[fHisCounter];
    }
    
    // this shouldn't happen
    vector< double > f;
    return f;
}


TList* VStereoAnalysis::getEnergyHistograms()
{
    if( fHisCounter < 0 )
    {
        return fHistoTot->hListEnergyHistograms;
    }
    else if( fHisCounter < ( int )fHisto.size() )
    {
        return fHisto[fHisCounter]->hListEnergyHistograms;
    }
    
    return 0;
}


TH1D* VStereoAnalysis::getTheta2()
{
    if( fHisCounter < 0 )
    {
        return fHistoTot->htheta2;
    }
    else if( fHisCounter < ( int )fHisto.size() )
    {
        return fHisto[fHisCounter]->htheta2;
    }
    
    return 0;
}


double VStereoAnalysis::getRawRate()
{
    if( fHisCounter < 0 )
    {
        return 0.;
    }
    else if( fHisCounter < ( int )fHisto.size() )
    {
        return fHisto[fHisCounter]->hrate_1sec->GetEntries();
    }
    
    return 0.;
}


CData* VStereoAnalysis::getDataFromFile( int i_runNumber )
{
    cout << "VStereoAnalysis::getDataFromFile Getting Data from file!" << endl;
    CData* c = 0;
    for( unsigned int i = 0; i < fRunPara->fRunList.size(); i++ )
    {
        int i_run = fIsOn ? fRunPara->fRunList[i].fRunOn : fRunPara->fRunList[i].fRunOff;
        
        if( i_runNumber > 0 && i_runNumber != i_run )
        {
            continue;
        }
        
        string iFileName = fIsOn ? fRunPara->fRunList[i].fRunOnFileName : fRunPara->fRunList[i].fRunOffFileName;
        
        fDataFile = new TFile( iFileName.c_str() );
        if( fDataFile->IsZombie() )
        {
            cout << "VStereoAnalysis::getDataFromFile() error opening file " << iFileName << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        fDataRunTree = ( TTree* )fDataFile->Get( "data" );
        if( !fDataRunTree )
        {
            cout << "VStereoAnalysis::getDataFromFile() error: cannot find data tree in " << iFileName << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        fDataFrogsTree = ( TTree* )fDataFile->Get( "frogspars" );
        if( fRunPara->fRunList[i].fIsFrogs )
        {
            if( !fDataFrogsTree )
            {
                cout << "VStereoAnalysis::getDataFromFile() info: cannot find frogspars tree in " << iFileName << endl;
                cout << "(this will lead to a failure for the frogs analysis, but is otherwise not a problem)" << endl;
            }
            else
            {
                cout << "VStereoAnalysis::getDataFromFile() info: found frogspars tree" << endl;
                //check that frogs tree and data tree have same number of events.
                if( fDataRunTree->GetEntries() == fDataFrogsTree->GetEntries() )
                {
                    fDataRunTree->AddFriend( fDataFrogsTree );
                }
                else
                {
                    cout << "VStereoAnalysis::getDataFromFile() Error: data tree has " << fDataRunTree->GetEntries() ;
                    cout << " entries; frogs tree has " << fDataFrogsTree->GetEntries() << " entries; that's not good. Exiting now." << endl;
                    exit( EXIT_FAILURE );
                }
                
            }
        }
        c = new CData( fDataRunTree );
        // read current epoch from data file
        VEvndispRunParameter* i_runPara = ( VEvndispRunParameter* )fDataFile->Get( "runparameterV2" );
        if( i_runPara )
        {
            fInstrumentEpoch = i_runPara->fInstrumentEpoch;
            fTelToAnalyze = i_runPara->fTelToAnalyze;
        }
        else
        {
            cout << "VStereoAnalysis::getDataFromFile() warning: epoch of current file " << endl;
            cout << "and active telescope combination cannot be determined; " << endl;
            cout << "this might lead to a wrong choice in the gamma/hadron cuts - please check" << endl;
            fInstrumentEpoch = "NOT_FOUND";
        }
    }
    return c;
}


bool VStereoAnalysis::closeDataFile()
{
    if( fDataFile )
    {
        fDataFile->Close();
    }
    
    return true;
}


bool VStereoAnalysis::terminate()
{
    closeDataFile();
    
    return true;
}

bool VStereoAnalysis::init_TreeWithSelectedEvents( int irun, bool isOn )
{
    //NOTE: This tree is currently allways filled with the eventdisplay reconstruction results. If you want to have frogs results, you need to fix that.
    if( fTreeSelectedEvents )
    {
        delete fTreeSelectedEvents;
    }
    if( !fRunPara )
    {
        return false;
    }
    
    char hname[200];
    char htitle[200];
    if( isOn )
    {
        sprintf( hname, "data_on" );
        sprintf( htitle, "selected events (on) for run %d", irun );
    }
    else
    {
        sprintf( hname, "data_off" );
        sprintf( htitle, "selected events (off) for run %d", irun );
    }
    fTreeSelectedEvents = new TTree( hname, htitle );
    if( !isOn )
    {
        fTreeSelectedEvents->SetLineColor( 2 );
    }
    
    fTreeSelectedEvents->Branch( "runNumber", &fTreeSelected_runNumber, "runNumber/I" );
    fTreeSelectedEvents->Branch( "eventNumber", &fTreeSelected_eventNumber, "eventNumber/I" );
    fTreeSelectedEvents->Branch( "MJD", &fTreeSelected_MJD, "MJD/I" );
    fTreeSelectedEvents->Branch( "Time", &fTreeSelected_Time, "Time/D" );
    fTreeSelectedEvents->Branch( "NImages", &fTreeSelected_NImages, "NImages/I" );
    fTreeSelectedEvents->Branch( "ImgSel", &fTreeSelected_ImgSel, "ImgSel/I" );
    fTreeSelectedEvents->Branch( "theta2", &fTreeSelected_theta2, "theta2/D" );
    fTreeSelectedEvents->Branch( "Xoff", &fTreeSelected_Xoff, "Xoff/D" );
    fTreeSelectedEvents->Branch( "Yoff", &fTreeSelected_Yoff, "Yoff/D" );
    fTreeSelectedEvents->Branch( "Xoff_derot", &fTreeSelected_Xoff_derot, "Xoff_derot/D" );
    fTreeSelectedEvents->Branch( "Yoff_derot", &fTreeSelected_Yoff_derot, "Yoff_derot/D" );
    fTreeSelectedEvents->Branch( "Xcore", &fTreeSelected_Xcore, "Xcore/D" );
    fTreeSelectedEvents->Branch( "Ycore", &fTreeSelected_Ycore, "Ycore/D" );
    fTreeSelectedEvents->Branch( "MSCW", &fTreeSelected_MSCW, "MSCW/D" );
    fTreeSelectedEvents->Branch( "MSCL", &fTreeSelected_MSCL, "MSCL/D" );
    fTreeSelectedEvents->Branch( "MWR", &fTreeSelected_MWR, "MWR/D" );
    fTreeSelectedEvents->Branch( "MLR", &fTreeSelected_MLR, "MLR/D" );
    fTreeSelectedEvents->Branch( "ErecS", &fTreeSelected_ErecS, "ErecS/D" );
    fTreeSelectedEvents->Branch( "EChi2S", &fTreeSelected_EChi2S, "EChi2S/D" );
    fTreeSelectedEvents->Branch( "dES", &fTreeSelected_dES, "dES/D" );
    fTreeSelectedEvents->Branch( "EmissionHeight", &fTreeSelected_EmissionHeight, "EmissionHeight/F" );
    fTreeSelectedEvents->Branch( "EmissionHeightChi2", &fTreeSelected_EmissionHeightChi2, "EmissionHeightChi2/F" );
    fTreeSelectedEvents->Branch( "SizeSecondMax", &fTreeSelected_SizeSecondMax, "SizeSecondMax/D" );
    fTreeSelectedEvents->Branch( "Az", &fTreeSelected_Az, "Az/D" );
    fTreeSelectedEvents->Branch( "El", &fTreeSelected_El, "El/D" );
    fTreeSelectedEvents->Branch( "MVA", &fTreeSelected_MVA, "MVA/D" );
    fTreeSelectedEvents->Branch( "IsGamma", &fTreeSelected_IsGamma, "IsGamma/i" );
    fTreeSelectedEvents->Branch( "passedDirectionCut", &fTreeSelected_DirectionCut, "passedDirectionCut/i" );
    
    if( fCuts && fCuts->useFrogsCuts() )
    {
        fTreeSelectedEvents->Branch( "frogsEventID", &fTreeSelescted_frogsEventID, "frogsEventID/I" );
        fTreeSelectedEvents->Branch( "frogsGSLConStat", &fTreeSelescted_frogsGSLConStat, "frogsGSLConStat/I" );
        fTreeSelectedEvents->Branch( "frogsNB_iter", &fTreeSelescted_frogsNB_iter, "frogsNB_iter/I" );
        fTreeSelectedEvents->Branch( "frogsNImages", &fTreeSelescted_frogsNImages, "frogsNImages/I" );
        fTreeSelectedEvents->Branch( "frogsXS", &fTreeSelescted_frogsXS, "frogsXS/D" );
        fTreeSelectedEvents->Branch( "frogsXSerr", &fTreeSelescted_frogsXSerr, "frogsXSerr/D" );
        fTreeSelectedEvents->Branch( "frogsYS", &fTreeSelescted_frogsYS, "frogsYS/D" );
        fTreeSelectedEvents->Branch( "frogsYSerr", &fTreeSelescted_frogsYSerr, "frogsYSerr/D" );
        fTreeSelectedEvents->Branch( "frogsXP", &fTreeSelescted_frogsXP, "frogsXP/D" );
        fTreeSelectedEvents->Branch( "frogsXPerr", &fTreeSelescted_frogsXPerr, "frogsXPerr/D" );
        fTreeSelectedEvents->Branch( "frogsYP", &fTreeSelescted_frogsYP, "frogsYP/D" );
        fTreeSelectedEvents->Branch( "frogsYPerr", &fTreeSelescted_frogsYPerr, "frogsYPerr/D" );
        fTreeSelectedEvents->Branch( "frogsXPGC", &fTreeSelescted_frogsXPGC, "frogsXPGC/D" );
        fTreeSelectedEvents->Branch( "frogsYPYC", &fTreeSelescted_frogsYPGC, "frogsYPGC/D" );
        fTreeSelectedEvents->Branch( "frogsEnergy", &fTreeSelescted_frogsEnergy, "frogsEnergy/D" );
        fTreeSelectedEvents->Branch( "frogsEnergyerr", &fTreeSelescted_frogsEnergyerr, "frogsEnergyerr/D" );
        fTreeSelectedEvents->Branch( "frogsLambda", &fTreeSelescted_frogsLambda, "frogsLambda/D" );
        fTreeSelectedEvents->Branch( "frogsLambdaerr", &fTreeSelescted_frogsLambdaerr, "frogsLambdaerr/D" );
        fTreeSelectedEvents->Branch( "frogsGoodnessImg", &fTreeSelescted_frogsGoodnessImg, "frogsGoodnessImg/D" );
        fTreeSelectedEvents->Branch( "frogsNpixImg", &fTreeSelescted_frogsNpixImg, "frogsNpixImg/I" );
        fTreeSelectedEvents->Branch( "frogsGoodnessBkg", &fTreeSelescted_frogsGoodnessBkg, "frogsGoodnessBkg/D" );
        fTreeSelectedEvents->Branch( "frogsNpixBkg", &fTreeSelescted_frogsNpixBkg, "frogsNpixBkg/I" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessImg0", &fTreeSelescted_frogsTelGoodnessImg0, "frogsTelGoodnessImg0/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessImg1", &fTreeSelescted_frogsTelGoodnessImg1, "frogsTelGoodnessImg1/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessImg2", &fTreeSelescted_frogsTelGoodnessImg2, "frogsTelGoodnessImg2/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessImg3", &fTreeSelescted_frogsTelGoodnessImg3, "frogsTelGoodnessImg3/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessBkg0", &fTreeSelescted_frogsTelGoodnessBkg0, "frogsTelGoodnessBkg0/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessBkg1", &fTreeSelescted_frogsTelGoodnessBkg1, "frogsTelGoodnessBkg1/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessBkg2", &fTreeSelescted_frogsTelGoodnessBkg2, "frogsTelGoodnessBkg2/D" );
        fTreeSelectedEvents->Branch( "frogsTelGoodnessBkg3", &fTreeSelescted_frogsTelGoodnessBkg3, "frogsTelGoodnessBkg3/D" );
        fTreeSelectedEvents->Branch( "frogsXS_derot", &fTreeSelescted_frogsXS_derot, "frogsXS_derot/D" );
        fTreeSelectedEvents->Branch( "frogsYS_derot", &fTreeSelescted_frogsYS_derot, "frogsYS_derot/D" );
        fTreeSelectedEvents->Branch( "frogs_theta2", &fTreeSelescted_frogs_theta2, "frogs_theta2/D" );
    }
    
    reset_TreeWithSelectedEvents();
    
    return true;
}

void VStereoAnalysis::reset_TreeWithSelectedEvents()
{
    fTreeSelected_runNumber = 0;
    fTreeSelected_eventNumber = 0;
    fTreeSelected_MJD = 0;
    fTreeSelected_Time = 0.;
    fTreeSelected_NImages = 0;
    fTreeSelected_ImgSel = 0;
    fTreeSelected_theta2 = 0.;
    fTreeSelected_Xoff = 0.;
    fTreeSelected_Yoff = 0.;
    fTreeSelected_Xoff_derot = 0.;
    fTreeSelected_Yoff_derot = 0.;
    fTreeSelected_Xcore = 0.;
    fTreeSelected_Ycore = 0.;
    fTreeSelected_MSCW = 0.;
    fTreeSelected_MSCL = 0.;
    fTreeSelected_MWR = 0.;
    fTreeSelected_MLR = 0.;
    fTreeSelected_ErecS = 0.;
    fTreeSelected_EChi2S = 0.;
    fTreeSelected_dES = 0.;
    fTreeSelected_EmissionHeight = 0.;
    fTreeSelected_EmissionHeightChi2 = 0.;
    fTreeSelected_SizeSecondMax = 0.;
    fTreeSelected_Az = 0.;
    fTreeSelected_El = 0.;
    fTreeSelected_MVA = -99.;
    fTreeSelected_IsGamma = 0;
    fTreeSelected_DirectionCut = 0;
    
    /// frogs ///
    fTreeSelescted_frogsEventID = 0;
    fTreeSelescted_frogsGSLConStat = 0;
    fTreeSelescted_frogsNB_iter = 0;
    fTreeSelescted_frogsNImages = 0;
    fTreeSelescted_frogsXS = 0.;
    fTreeSelescted_frogsXSerr = 0.;
    fTreeSelescted_frogsYS = 0.;
    fTreeSelescted_frogsYSerr = 0.;
    fTreeSelescted_frogsXP = 0.;
    fTreeSelescted_frogsXPerr = 0.;
    fTreeSelescted_frogsYP = 0.;
    fTreeSelescted_frogsYPerr = 0.;
    fTreeSelescted_frogsXPGC = 0.;
    fTreeSelescted_frogsYPGC = 0.;
    fTreeSelescted_frogsEnergy = 0.;
    fTreeSelescted_frogsEnergyerr = 0.;
    fTreeSelescted_frogsLambda = 0.;
    fTreeSelescted_frogsLambdaerr = 0.;
    fTreeSelescted_frogsGoodnessImg = 0.;
    fTreeSelescted_frogsNpixImg = 0;
    fTreeSelescted_frogsGoodnessBkg = 0.;
    fTreeSelescted_frogsNpixBkg = 0;
    fTreeSelescted_frogsTelGoodnessImg0 = 0.;
    fTreeSelescted_frogsTelGoodnessImg1 = 0.;
    fTreeSelescted_frogsTelGoodnessImg2 = 0.;
    fTreeSelescted_frogsTelGoodnessImg3 = 0.;
    fTreeSelescted_frogsTelGoodnessBkg0 = 0.;
    fTreeSelescted_frogsTelGoodnessBkg1 = 0.;
    fTreeSelescted_frogsTelGoodnessBkg2 = 0.;
    fTreeSelescted_frogsTelGoodnessBkg3 = 0.;
    fTreeSelescted_frogsXS_derot = 0.;
    fTreeSelescted_frogsYS_derot = 0.;
    fTreeSelescted_frogs_theta2  = 0.;
}

void VStereoAnalysis::fill_TreeWithSelectedEvents( CData* c, double i_xderot, double i_yderot, double i_theta2, bool i_bDirectionCut )
{
    if( !c )
    {
        return;
    }
    
    fTreeSelected_runNumber = c->runNumber;
    fTreeSelected_eventNumber = c->eventNumber;
    fTreeSelected_MJD = c->MJD;
    fTreeSelected_Time = c->Time;
    fTreeSelected_NImages = c->NImages;
    fTreeSelected_ImgSel = ( int )c->ImgSel;
    fTreeSelected_theta2 = c->theta2;
    fTreeSelected_Xoff = c->Xoff;
    fTreeSelected_Yoff = c->Yoff;
    fTreeSelected_Xoff_derot = c->Xoff_derot;
    fTreeSelected_Yoff_derot = c->Yoff_derot;
    fTreeSelected_Xcore = c->Xcore;
    fTreeSelected_Ycore = c->Ycore;
    fTreeSelected_MSCW = c->MSCW;
    fTreeSelected_MSCL = c->MSCL;
    fTreeSelected_MWR = c->MWR;
    fTreeSelected_MLR = c->MLR;
    fTreeSelected_ErecS = c->ErecS;
    fTreeSelected_EChi2S = c->EChi2S;
    fTreeSelected_dES = c->dES;
    fTreeSelected_EmissionHeight = c->EmissionHeight;
    fTreeSelected_EmissionHeightChi2 = c->EmissionHeightChi2;
    fTreeSelected_SizeSecondMax = c->SizeSecondMax;
    fTreeSelected_Az = c->Az;
    fTreeSelected_El = 90. - c->Ze;
    if( fCuts )
    {
        fTreeSelected_MVA = fCuts->getTMVA_EvaluationResult();
    }
    else
    {
        fTreeSelected_MVA = -99.;
    }
    
    if( bIsGamma )
    {
        fTreeSelected_IsGamma = 1;
    }
    else
    {
        fTreeSelected_IsGamma = 0;
    }
    if( i_bDirectionCut )
    {
        fTreeSelected_DirectionCut = 1;
    }
    else
    {
        fTreeSelected_DirectionCut = 0;
    }
    
    /// frogs ///
    fTreeSelescted_frogsEventID     = c->frogsEventID;
    fTreeSelescted_frogsGSLConStat  = c->frogsGSLConStat;
    fTreeSelescted_frogsNB_iter     = c->frogsNB_iter;
    fTreeSelescted_frogsNImages     = c->frogsNImages;
    fTreeSelescted_frogsXS          = c->frogsXS;
    fTreeSelescted_frogsXSerr       = c->frogsXSerr;
    fTreeSelescted_frogsYS          = c->frogsYS;
    fTreeSelescted_frogsYSerr       = c->frogsYSerr;
    fTreeSelescted_frogsXP          = c->frogsXP;
    fTreeSelescted_frogsXPerr       = c->frogsXPerr;
    fTreeSelescted_frogsYP          = c->frogsYP;
    fTreeSelescted_frogsYPerr       = c->frogsYPerr;
    fTreeSelescted_frogsXPGC        = c->frogsXPGC;
    fTreeSelescted_frogsYPGC        = c->frogsYPGC;
    fTreeSelescted_frogsEnergy      = c->frogsEnergy;
    fTreeSelescted_frogsEnergyerr   = c->frogsEnergyerr;
    fTreeSelescted_frogsLambda      = c->frogsLambda;
    fTreeSelescted_frogsLambdaerr   = c->frogsLambdaerr;
    fTreeSelescted_frogsGoodnessImg = c->frogsGoodnessImg;
    fTreeSelescted_frogsNpixImg     = c->frogsNpixImg;
    fTreeSelescted_frogsGoodnessBkg = c->frogsGoodnessBkg;
    fTreeSelescted_frogsNpixBkg     = c->frogsNpixBkg;
    fTreeSelescted_frogsTelGoodnessImg0 = c->frogsTelGoodnessImg0;
    fTreeSelescted_frogsTelGoodnessImg1 = c->frogsTelGoodnessImg1;
    fTreeSelescted_frogsTelGoodnessImg2 = c->frogsTelGoodnessImg2;
    fTreeSelescted_frogsTelGoodnessImg3 = c->frogsTelGoodnessImg3;
    fTreeSelescted_frogsTelGoodnessBkg0 = c->frogsTelGoodnessBkg0;
    fTreeSelescted_frogsTelGoodnessBkg1 = c->frogsTelGoodnessBkg1;
    fTreeSelescted_frogsTelGoodnessBkg2 = c->frogsTelGoodnessBkg2;
    fTreeSelescted_frogsTelGoodnessBkg3 = c->frogsTelGoodnessBkg3;
    fTreeSelescted_frogsXS_derot = i_xderot;
    fTreeSelescted_frogsYS_derot = i_yderot;
    fTreeSelescted_frogs_theta2  = i_theta2;
    
    if( fTreeSelectedEvents )
    {
        fTreeSelectedEvents->Fill();
    }
    
}

/*
 * get derotated camera coordinates in J2000
 *
 */
void VStereoAnalysis::getDerotatedCoordinates( unsigned int icounter, CData* iDataRun, double& x_derot, double& y_derot )
{
    if( !iDataRun )
    {
        return;
    }
    if( icounter >= fRunPara->fRunList.size() )
    {
        return;
    }
    x_derot = iDataRun->getXoff_derot();
    y_derot = iDataRun->getYoff_derot();
    
    // convert de-rotated camera coordinates to J2000
    VSkyCoordinatesUtilities::convert_derotatedCoordinates_to_J2000(
        fDataRun->MJD, fDataRun->Time,
        iDataRun->ArrayPointing_Azimuth, iDataRun->ArrayPointing_Elevation,
        x_derot, y_derot );
}

double VStereoAnalysis::getWobbleNorth()
{
    if( fRunPara && fHisCounter >= 0 && fHisCounter < ( int )fRunPara->fRunList.size() )
    {
        return fRunPara->fRunList[fHisCounter].fWobbleNorthMod;
    }
    
    return 0.;
}

double VStereoAnalysis::getWobbleWest()
{
    if( fRunPara && fHisCounter >= 0 && fHisCounter < ( int )fRunPara->fRunList.size() )
    {
        return fRunPara->fRunList[fHisCounter].fWobbleWestMod;
    }
    
    return 0.;
}

bool VStereoAnalysis::init_TreeWithEventsForCtools( int irun ) // WRITEEVENTTREEFORCTOOLS
{
    //NOTE: This tree is currently allways filled with the eventdisplay reconstruction results. If you want to have frogs results, you need to fix that.
    cout << endl;
    cout << " :: init_TreeWithEventsForCtools( " << irun << " )" << endl;
    cout << endl;
    
    
    if( fTreeWithEventsForCtools )
    {
        delete fTreeWithEventsForCtools;
    }
    if( !fRunPara )
    {
        return false;
    }
    
    char hname[200];
    char htitle[200];
    sprintf( hname, "TreeWithEventsForCtools" );
    sprintf( htitle, "all gamma events with X,Y and Time for run %d, in a format for CTOOL's Event List format", irun );
    
    fTreeWithEventsForCtools = new TTree( hname, htitle );
    
    fTreeWithEventsForCtools->Branch( "runNumber",      &fTreeSelected_runNumber,      "runNumber/I" );      // runNumber
    fTreeWithEventsForCtools->Branch( "eventNumber",    &fTreeSelected_eventNumber,    "eventNumber/I" );    // eventNumber
    fTreeWithEventsForCtools->Branch( "timeOfDay",      &fTreeSelected_Time,           "timeOfDay/D" );      // Time
    fTreeWithEventsForCtools->Branch( "dayMJD",         &fTreeSelected_MJD,            "dayMJD/I" );         // MJD
    fTreeWithEventsForCtools->Branch( "EnergyS",        &fTreeSelected_ErecS,          "ErecS/D" );          // ErecS
    fTreeWithEventsForCtools->Branch( "EnergyS_Err",    &fTreeSelected_dES,      "ErecS_Err/D" );      // dES
    fTreeWithEventsForCtools->Branch( "XGroundCore",    &fTreeSelected_Xcore,    "XGroundCore/D" );    // Xcore
    fTreeWithEventsForCtools->Branch( "YGroundCore",    &fTreeSelected_Ycore,    "YGroundCore/D" );    // Ycore
    fTreeWithEventsForCtools->Branch( "Xderot",         &fTreeSelected_Xoff_derot,         "Xderot/D" );         // Xoff_derot
    fTreeWithEventsForCtools->Branch( "Yderot",         &fTreeSelected_Yoff_derot,         "Yderot/D" );         // Yoff_derot
    fTreeWithEventsForCtools->Branch( "NImages",        &fTreeSelected_NImages,        "NImages/I" );        // NImages
    fTreeWithEventsForCtools->Branch( "ImgSel",         &fTreeSelected_ImgSel,         "ImgSel/l" );         // ImgSel
    fTreeWithEventsForCtools->Branch( "MSCW",           &fTreeSelected_MSCW,           "MSCW/D" );           // MSCW
    fTreeWithEventsForCtools->Branch( "MSCL",           &fTreeSelected_MSCL,           "MSCL/D" );           // MSCL
    fTreeWithEventsForCtools->Branch( "MWR"           , &fTreeSelected_MWR,            "MWR/D" );            // MWR
    fTreeWithEventsForCtools->Branch( "MLR"           , &fTreeSelected_MLR,            "MLR/D" );            // MLR
    //	fTreeWithEventsForCtools->Branch( "TargetRA"      , &fTreeCTOOLS_TargetRA,       "TargetRA/D" );
    //	fTreeWithEventsForCtools->Branch( "TargetDEC"     , &fTreeCTOOLS_TargetDEC,      "TargetDEC/D" );
    fTreeWithEventsForCtools->Branch( "RA"            , &fTreeCTOOLS_RA,             "RA/D" );
    fTreeWithEventsForCtools->Branch( "DEC"           , &fTreeCTOOLS_DEC,            "DEC/D" );
    fTreeWithEventsForCtools->Branch( "Az"            , &fTreeCTOOLS_Az,             "Az/D" );
    fTreeWithEventsForCtools->Branch( "El"            , &fTreeCTOOLS_El,             "El/D" );
    fTreeWithEventsForCtools->Branch( "EmissionHeight", &fTreeSelected_EmissionHeight, "EmissionHeight/D" );   // EmissionHeight
    fTreeWithEventsForCtools->Branch( "Xoff"          , &fTreeSelected_Xoff, "Xoff/D" );             // Xoff
    fTreeWithEventsForCtools->Branch( "Yoff"          , &fTreeSelected_Yoff, "Yoff/D" );             // Yoff
    //	fTreeWithEventsForCtools->Branch( "GregYear"      , &fTreeCTOOLS_GregYear      , "GregYear/D" );
    //	fTreeWithEventsForCtools->Branch( "GregMonth"     , &fTreeCTOOLS_GregMonth     , "GregMonth/D" );
    //	fTreeWithEventsForCtools->Branch( "GregDay"       , &fTreeCTOOLS_GregDay       , "GregDay/D" );
    //	fTreeWithEventsForCtools->Branch( "Acceptance"    , &fTreeCTOOLS_Acceptance    , "Acceptance/D" );
    cout << endl;
    
    // init acceptance critter
    fCTOOLSAcceptance = new VRadialAcceptance( fRunPara->fRunList[0].fAcceptanceFile ) ;
    fCTOOLSAcceptance->Set2DAcceptanceMode( fRunPara->fRunList[0].f2DAcceptanceMode ) ;
    
    cout << " :: init_TreeWithEventsForCtools()" << endl;
    cout << endl;
    
    return true;
}

void VStereoAnalysis::fill_TreeWithEventsForCtools( CData* c , double i_xderot, double i_yderot ) // WRITEEVENTTREEFORCTOOLS
{
    if( !c )
    {
        return;
    }
    
    fTreeSelected_runNumber    = c->runNumber;    // Run Number
    fTreeSelected_eventNumber  = c->eventNumber;  // Event Number
    printf( "filling event %d to TreeWithEventsForCtools...\n", fTreeSelected_eventNumber ) ;
    fTreeSelected_Time         = c->Time;         // Time of day (seconds) of gamma ray event
    fTreeSelected_MJD          = c->MJD;          // Day of epoch (days)
    fTreeSelected_Xoff         = c->getXoff();         // Gamma Point-Of-Origin, in camera coodinates (deg)
    fTreeSelected_Yoff         = c->getYoff();         // Gamma Point-Of-Origin, in camera coodinates (deg)
    fTreeSelected_Xoff_derot   = i_xderot;        // Derotated Gamma Point-Of-Origin (deg, RA)
    fTreeSelected_Yoff_derot   = i_yderot;        // Derotated Gamma Point-Of-Origin (deg, DEC)
    fTreeSelected_ErecS        = c->getEnergy_TeV();        // Reconstructed Gamma Energy (TeV)
    fTreeSelected_dES          = c->getEnergyDelta();          // Reconstructed Gamma Energy (TeV) Error
    fTreeSelected_Xcore        = c->getXcore_M();        // Gamma Ray Core-Ground intersection location (north?)
    fTreeSelected_Ycore        = c->getYcore_M();        // Gamma Ray Core-Ground intersection location (east?)
    fTreeSelected_NImages      = c->getNImages();      // Number of images used in reconstruction?
    fTreeSelected_ImgSel       = c->getImgSel();       // 4-bit binary code describing which telescopes had images
    fTreeSelected_MSCW         = c->MSCW ;        // mean scaled width
    fTreeSelected_MSCL         = c->MSCL ;        // mean scaled length
    fTreeSelected_MWR          = c->MWR ;
    fTreeSelected_MLR          = c->MLR ;
    fTreeSelected_EmissionHeight = c->EmissionHeight ; // height of shower maximum (in km) above telescope z-plane
    fTreeCTOOLS_Acceptance     = fCTOOLSAcceptance->getAcceptance( i_xderot, i_yderot ) ;
    
    // pointing target
    fTreeCTOOLS_TargetRA  = fRunPara->fRunList[0].fTargetRAJ2000  ;
    fTreeCTOOLS_TargetDEC = fRunPara->fRunList[0].fTargetDecJ2000 ;
    
    // get telescope pointing in the current epoch
    double CenterPoint_RA  = fAstro.back()->getTelRA() ;  // radians
    double CenterPoint_DEC = fAstro.back()->getTelDec() ; // radians
    
    // precess telescope pointing to J2000 epoch
    VSkyCoordinatesUtilities::precessTarget( 51544., CenterPoint_RA, CenterPoint_DEC, c->MJD, false ) ;
    
    // Convert derotated detx/dety/telpointing to spherical RA/Dec
    double Spherical_RA  = 0.0 ;
    double Spherical_DEC = 0.0 ;
    VAstronometry::vlaDtp2s( fTreeSelected_Xoff_derot * TMath::DegToRad(),
              fTreeSelected_Yoff_derot * TMath::DegToRad(),
              CenterPoint_RA,   CenterPoint_DEC,
              &Spherical_RA,    &Spherical_DEC ) ;
              
    // at this point, the ra/dec is the unprecessed 'J2000' coordinates,
    // ctools needs the 'fk5' radec coordinates
    
    // Convert from spherical RA and DEC to Azimuth and Zenith
    double Az_deg = 0.0 ;
    double El_deg = 0.0 ;
    
    // convert to degrees and do calculation
    double Spherical_RA_deg  = Spherical_RA  * TMath::RadToDeg() ;
    double Spherical_DEC_deg = Spherical_DEC * TMath::RadToDeg() ;
    fVsky->setTargetJ2000( Spherical_DEC_deg , Spherical_RA_deg ) ;
    fVsky->precessTarget( fTreeSelected_MJD, 0 ) ;
    fVsky->updatePointing( fTreeSelected_MJD, fTreeSelected_Time ) ;
    
    // this should be the event's rad/dec coordinates in the precessed 'fk5' system
    fTreeCTOOLS_RA  = fVsky->getTargetRA()  ;
    fTreeCTOOLS_DEC = fVsky->getTargetDec() ;
    
    // calculate new param
    Az_deg = fVsky->getTargetAzimuth() ;
    El_deg = fVsky->getTargetElevation() ;
    fTreeCTOOLS_Az = Az_deg ;
    fTreeCTOOLS_El = El_deg ;
    
    // Convert MJD to Year, Month, and Day
    fTreeCTOOLS_GregYear  = 0 ;
    fTreeCTOOLS_GregMonth = 0 ;
    fTreeCTOOLS_GregDay   = 0 ;
    double junk1 = 0.0 ;
    int junk2 = 0 ;
    VAstronometry::vlaDjcl( c->MJD, &fTreeCTOOLS_GregYear, &fTreeCTOOLS_GregMonth, &fTreeCTOOLS_GregDay, &junk1, &junk2 ) ;
    
    
    if( fTreeWithEventsForCtools )
    {
        fTreeWithEventsForCtools->Fill();
        printf( "CTOOLS EVENT %d %d %f %f %f %f %f\n",
                fTreeSelected_runNumber,
                fTreeSelected_eventNumber,
                fTreeSelected_Xoff,
                fTreeSelected_Yoff,
                fTreeSelected_ErecS,
                fTreeCTOOLS_RA,
                fTreeCTOOLS_DEC ) ;
    }
    
}

void VStereoAnalysis::save_TreeWithEventsForCtools() // WRITEEVENTTREEFORCTOOLS
{
    // save our ctools tree
    printf( "Preparing to write TreeWithEventsForCtools with %lld events...\n", fTreeWithEventsForCtools->GetEntries() ) ;
    fTreeWithEventsForCtools->Write() ; // or maybe ->AutoSave() ?
    
    fRunPara->SetName( "VAnaSumRunParameter" );
    fRunPara->Write() ;
    
    fCDataTreeClone = fDataRunTree->CloneTree() ;
    fCDataTreeClone->SetName( "cdatatree" );
    fCDataTreeClone->Write();
}


