/*! \class VAnaSum
    brief class for producing an analysis summary from parameterized veritas data

*/

#include "VAnaSum.h"

VAnaSum::VAnaSum( string i_datadir )
{
    fAnalysisRunMode = 0;
    
    fDatadir = i_datadir + "/";
    fPrefix = "";
    fSuffix = ".mscw.root";
    
    fRunPara = 0;
    fStereoOn = 0;
    fStereoOff = 0;
    fRunSummary = 0;
    
    fOPfile = 0;
    fTotalDir = 0;
    fTotalDirName = "total_1";
    fStereoTotalDir = 0;
}

/*
   initialize data analysis

   check run mode

*/
void VAnaSum::initialize( string i_LongListFilename, unsigned int iRunType, string i_outfile, int iRandomSeed, string iRunParameterfile )
{
    char i_temp[2000];
    char i_title[200];
    
    ///////////////////////////////////////////////////////////////////////////////
    // check analysis run mode
    fAnalysisRunMode = iRunType;
    // merging analysis
    if( fAnalysisRunMode == 1 )
    {
        cout << "merging analysis results" << endl << endl;
    }
    // sequentiell analysis
    else
    {
        cout << "sequentiell analysis" << endl << endl;
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // define run parameter, load list of runs
    fRunPara = new VAnaSumRunParameter();
    if( !fRunPara->readRunParameter( iRunParameterfile ) )
    {
        cout << "error while reading run parameters" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    int i_npair = fRunPara->loadLongFileList( i_LongListFilename, false, ( fAnalysisRunMode == 1 ) );
    if( i_npair == 0 && iRunType == 1 )
    {
        cout << "\t run list not of long type, try simple format" << endl;
        i_npair = fRunPara->loadSimpleFileList( i_LongListFilename );
        cout << "\t found " << i_npair << " runs in simple run list" << endl;
    }
    if( i_npair == 0 )
    {
        cout << "VAnaSum error: no files found in runlist" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    if( fAnalysisRunMode != 1 )
    {
        cout << "Random seed for stereo maps: " << iRandomSeed << endl;
    }
    cout << endl;
    cout << "File with list of runs: " << i_LongListFilename;
    cout << " (file format version " << fRunPara->getInputFileVersionNumber() << ", " << i_npair << " runs found in list)" << endl;
    
    if( fAnalysisRunMode != 1 )
    {
        if( fRunPara->getInputFileVersionNumber() > 3 )
        {
            fRunPara->getEventdisplayRunParameter( fDatadir );
        }
        else
        {
            fRunPara->getWobbleOffsets( fDatadir );
        }
    }
    fNumTels = fRunPara->getMaxNumberofTelescopes();
    
    ///////////////////////////////////////////////////////////////////////////////
    // open output root file for all the results and create directory structure
    fOPfile = new TFile( i_outfile.c_str(), "recreate" );
    if( fOPfile->IsZombie() )
    {
        cout << "Error: cannot open output file " << i_outfile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    // make directories with combined result
    fOPfile->cd();
    fTotalDir = fOPfile->mkdir( fTotalDirName.c_str(), "combined results from all runs" );
    if( !fTotalDir )
    {
        cout << "VAnaSum::initialize error creating directory " << fTotalDirName << " in output file " << fOPfile->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    fTotalDir->cd();
    fStereoTotalDir = fTotalDir->mkdir( "stereo", "combined stereo results" );
    if( !fStereoTotalDir )
    {
        cout << "VAnaSum::initialize error creating directory stereo in output file " << fOPfile->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    
    // directories with run wise results
    for( unsigned int j = 0; j < fRunPara->fRunList.size(); j++ )
    {
        sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), fRunPara->fRunList[j].fRunOn, fSuffix.c_str() );
        fRunPara->fRunList[j].fRunOnFileName = i_temp;
        sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), fRunPara->fRunList[j].fRunOff, fSuffix.c_str() );
        fRunPara->fRunList[j].fRunOffFileName = i_temp;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // merging analysis
    // (combine several anasum files into one analysis file)
    ///////////////////////////////////////////////////////////////////////////////////////////
    if( fAnalysisRunMode == 1 )
    {
        // loop over all files in run list and copy histograms
        for( unsigned int j = 0; j < fRunPara->fRunList.size(); j++ )
        {
            // open input file
            sprintf( i_temp, "%s/%d.anasum.root", fDatadir.c_str(), fRunPara->fRunList[j].fRunOn );
            TFile iAnasumInputFile( i_temp );
            if( iAnasumInputFile.IsZombie() )
            {
                cout << "VAnasum::initialize error opening anasum file for run  " << fRunPara->fRunList[j].fRunOn << ": " << endl;
                cout << iAnasumInputFile.GetName() << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            // copying histograms
            cout << "reading in histograms from " << iAnasumInputFile.GetName() << endl;
            fOPfile->cd();
            sprintf( i_temp, "run_%d", fRunPara->fRunList[j].fRunOn );
            TDirectory* iDir = iAnasumInputFile.GetDirectory( i_temp );
            if( iDir )
            {
                cout << "\t copy histograms for run " << fRunPara->fRunList[j].fRunOn << endl;
                copyDirectory( iDir );
                fOldRunList.insert( fRunPara->fRunList[j].fRunOn );
            }
            else
            {
                cout << endl;
                cout << "Fatal error: directory not found for run " << fRunPara->fRunList[j].fRunOn << endl;
                cout << "...exiting" << endl;
                exit( EXIT_FAILURE );
            }
            iDir = fOPfile->GetDirectory( i_temp );
            fStereoRunDir.push_back( iDir );
            
            // get list of excluded regions
            // (note: list is read only from first file)
            if( j == 0 )
            {
                fRunPara->getListOfExcludedSkyRegions( &iAnasumInputFile, fRunPara->fRunList[j].fRunOn );
            }
            
            // close input file
            iAnasumInputFile.Close();
        }
    }
    //////////////////////////////////////////////////////////////////
    // STANDARD ANALYSIS (analysis all individual runs; sequentiell)
    else
    {
        for( unsigned int j = 0; j < fRunPara->fRunList.size(); j++ )
        {
            sprintf( i_temp, "run_%d", fRunPara->fRunList[j].fRunOn );
            sprintf( i_title, "results for run %d", fRunPara->fRunList[j].fRunOn );
            fOPfile->cd();
            TDirectory* iDir = fOPfile->GetDirectory( i_temp );
            // update analysis
            if( fAnalysisRunMode == 1 && iDir )
            {
                cout << "Found directory for run " << fRunPara->fRunList[j].fRunOn << endl;
                fOldRunList.insert( fRunPara->fRunList[j].fRunOn );
                fRunDir.push_back( iDir );
                fStereoRunDir.push_back( iDir->GetDirectory( "stereo" ) );
            }
            else
            {
                fRunDir.push_back( fOPfile->mkdir( i_temp, i_title ) );
                if( !fRunDir.back() )
                {
                    cout << "VAnaSum::initialize error creating run directory " << i_temp << " in output file " << fOPfile->GetName() << endl;
                    cout << "(run " << fRunPara->fRunList[j].fRunOn << ")" << endl;
                    exit( EXIT_FAILURE );
                }
                
                sprintf( i_temp, "stereo" );
                sprintf( i_title, "stereo results for run %d", fRunPara->fRunList[j].fRunOn );
                fRunDir.back()->cd();
                fStereoRunDir.push_back( fRunDir.back()->mkdir( i_temp, i_title ) );
                
                if( !fStereoRunDir.back() )
                {
                    cout << "VAnaSum::initialize error creating stereo run directory ";
                    cout << i_temp << " in output file " << fOPfile->GetName() << endl;
                    cout << "(run " << fRunPara->fRunList[j].fRunOn << ")" << endl;
                    exit( EXIT_FAILURE );
                }
            }
        }
        
        //////////////////////////////////////////////////
        // set up data chains (all runs combined)
        // get wobble offsets
        // get run duration
        fTotalExposureOn = 0.;
        cout << endl;
        cout << "-----------------------------------------------------------------------(3)" << endl;
        for( unsigned int j = 0; j < fRunPara->fRunList.size(); j++ )
        {
            sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), fRunPara->fRunList[j].fRunOn, fSuffix.c_str() );
            cout << "Chaining run " << j << " of " << fRunPara->fRunList.size() << " runs with source data: " << i_temp << endl;
            
            // get azimuth range
            double azmin, azmax = 0.;
            fRunAzMeanOn[fRunPara->fRunList[j].fRunOn] = getAzRange( fRunPara->fRunList[j].fRunOn, "data", azmin, azmax );
            fRunAzMinOn[fRunPara->fRunList[j].fRunOn] = azmin;
            fRunAzMaxOn[fRunPara->fRunList[j].fRunOn] = azmax;
            
            // get mean noise level
            fRunPedVarsOn[fRunPara->fRunList[j].fRunOn] = getNoiseLevel( fRunPara->fRunList[j].fRunOn );
        }
        fTotalExposureOff = 0.;
        
        cout << "-----------------------------------------------------------------------(4)" << endl;
        for( unsigned int j = 0; j < fRunPara->fRunList.size(); j++ )
        {
            sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), fRunPara->fRunList[j].fRunOff, fSuffix.c_str() );
            cout << "Chaining run " << j << " of " << fRunPara->fRunList.size() << " runs with background data: " << i_temp << endl;
            
            // get azimuth range
            double azmin, azmax = 0.;
            fRunAzMeanOff[fRunPara->fRunList[j].fRunOff] = getAzRange( fRunPara->fRunList[j].fRunOff, "data", azmin, azmax );
            fRunAzMinOff[fRunPara->fRunList[j].fRunOff] = azmin;
            fRunAzMaxOff[fRunPara->fRunList[j].fRunOff] = azmax;
            
            // get noise level
            fRunPedVarsOff[fRunPara->fRunList[j].fRunOff] = getNoiseLevel( fRunPara->fRunList[j].fRunOff );
            
            // cp some information over to anasum file from mscw file
            if( j < fStereoRunDir.size() && fStereoRunDir[j] )
            {
                fStereoRunDir[j]->cd();
                char i_temp1[2000];
                sprintf( i_temp1, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), fRunPara->fRunList[j].fRunOn, fSuffix.c_str() );
                TFile* oldfile = new TFile( i_temp1 );
                if( !oldfile->IsZombie() )
                {
                    // copy TTree telconfig to anasum.root file
                    TTree* iTree = ( TTree* )oldfile->Get( "telconfig" );
                    if( iTree )
                    {
                        fStereoRunDir[j]->cd();
                        TTree* newtree = iTree->CloneTree();
                        if( newtree )
                        {
                            newtree->Write();
                            delete newtree;
                        }
                    }
                    // copy TTree pointingDataReduced to anasum.root file
                    TTree* jTree = ( TTree* )oldfile->Get( "pointingDataReduced" );
                    if( jTree )
                    {
                        fStereoRunDir[j]->cd();
                        TTree* newtreej = jTree->CloneTree();
                        if( newtreej )
                        {
                            newtreej->Write();
                            delete newtreej;
                        }
                    }
                    
                    // copy TTree deadPixelRegistry to anasum.root file
                    TTree* kTree = ( TTree* )oldfile->Get( "deadPixelRegistry" ) ;
                    if( kTree )
                    {
                        fStereoRunDir[j]->cd();
                        TTree* newtreek = kTree->CloneTree();
                        if( newtreek )
                        {
                            newtreek->Write();
                            delete newtreek;
                        }
                    }
                }
                // copy VEvndispRunParameter 'runparameterV2' to anasum.root file
                fStereoRunDir[j]->cd();
                VEvndispRunParameter* evnrunpar = ( VEvndispRunParameter* )oldfile->Get( "runparameterV2" );
                evnrunpar->Write();
                delete oldfile;
            }
            
        }
    }
    cout << "-----------------------------------------------------------------------(5)" << endl;
    
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // VStereoAnalysis will check validity of the data trees.
    fStereoOn = new VStereoAnalysis( true, "on", fRunPara, fStereoRunDir, fStereoTotalDir, fDatadir, iRandomSeed, ( fAnalysisRunMode == 1 ) );
    fRunExposureOn = fStereoOn->getRunDuration(); // Default exposure is run duration
    fStereoOff = new VStereoAnalysis( false, "off", fRunPara, fStereoRunDir, fStereoTotalDir, fDatadir, iRandomSeed, ( fAnalysisRunMode == 1 ) );
    // Default exposure is run duration
    fRunExposureOff = fStereoOff->getRunDuration();
    
    // rate plots and run summary
    // (rate plot times are relevant for ON runs only)
    fRunSummary = new VRunSummary();
    
    fMeanRawRateOn = 0.;
    fMeanRawRateOff = 0.;
    fMeanElevationOn = 0.;
    fMeanElevationOff = 0.;
    fNMeanElevation = 0.;
    fMeanAzimuthOn = 0.;
    fMeanAzimuthOff = 0.;
    fMeanDeadTimeOn = 0.;
    fMeanDeadTimeOff = 0.;
    fMeanPedVarsOn = 0.;
    fMeanPedVarsOff = 0.;
    
    //////////////////////////////////////////////////////////
    // fill runsummary (running from anasum, new total only)
    if( fAnalysisRunMode == 1 )
    {
        fRunSummary->fill( fDatadir, fTotalDirName, fRunPara->fRunList );
        
        fTotalExposureOn = fRunSummary->fTotalExposureOn;
        fTotalExposureOff = fRunSummary->fTotalExposureOff;
        fRunExposureOn  = fRunSummary->f_exposureOn;
        fRunExposureOff = fRunSummary->f_exposureOff;
        fStereoOn->setRunMJD( fRunSummary->fRunMJD );
        fStereoOff->setRunMJD( fRunSummary->fRunMJD );
        fMeanAzimuthOn  = fRunSummary->fMeanAzimuthOn;
        fMeanAzimuthOff = fRunSummary->fMeanAzimuthOff;
        fMeanElevationOn = fRunSummary->fMeanElevationOn;
        fMeanElevationOff = fRunSummary->fMeanElevationOff;
        fNMeanElevation = fRunSummary->fNMeanElevation;
        fMeanDeadTimeOn = fRunSummary->fMeanDeadTimeOn;
        fMeanDeadTimeOff = fRunSummary->fMeanDeadTimeOff;
        fMeanRawRateOn  = fRunSummary->fMeanRawRateOn;
        fMeanRawRateOff = fRunSummary->fMeanRawRateOff;
        fMeanPedVarsOn  = fRunSummary->fMeanPedVarsOn;
        fMeanPedVarsOff = fRunSummary->fMeanPedVarsOff;
        if( fStereoOn )
        {
            fStereoOn->setRunExposure( fRunExposureOn );
        }
        if( fStereoOff )
        {
            fStereoOff->setRunExposure( fRunExposureOff );
        }
    }
}


/*!
 *   do stereo analysis
 *
 *   loop first over all runs in run list, then do combined analysis
 *
 */
void VAnaSum::doStereoAnalysis()
{
    if( fAnalysisRunMode != 1 )
    {
        // do stereo analysis for each run
        for( unsigned int i = 0; i < fRunPara->fRunList.size(); i++ )
        {
            cout << endl;
            cout << "---------------------------------------" << endl;
            fRunPara->printStereoParameter( i );
            
            // analyze run
            doStereoAnalysis( i, fRunPara->fRunList[i].fRunOn, fRunPara->fRunList[i].fRunOff, fStereoRunDir[i] );
            
            cout << "---------------------------------------" << endl;
        }
    }
    fStereoTotalDir->cd();
    
    // now combine all runs to give 'total' results
    cout << endl;
    cout << "---------------------------------------" << endl;
    cout << "Stereo analysis for all runs:" << endl;
    
    doStereoAnalysis( -1, -1, -1, fStereoTotalDir );
    
    cout << "---------------------------------------" << endl;
    
    // print table with all results
    if( fRunSummary )
    {
        fRunSummary->print();
    }
}


/*
 *  run stereo analysis for the given run (sky map and histogram filling, on-off)
 *
 *  onrun == -1 : combined analysis
 *
 */

void VAnaSum::doStereoAnalysis( int icounter, int onrun, int offrun, TDirectory* idir )
{
    if( !idir->cd() )
    {
        cout << "VAnaSum::doStereoAnalysis error, directory not found " << endl;
    }
    
    ////////////////////////////////////////////////////////////
    // fill on and off histograms and sky maps
    fStereoOn->fillHistograms( icounter, onrun, fRunAzMinOn[onrun], fRunAzMaxOn[onrun], fRunPedVarsOn[onrun] );
    fStereoOff->fillHistograms( icounter, offrun, fRunAzMinOff[offrun], fRunAzMaxOff[offrun], fRunPedVarsOff[offrun] );
    
    ////////////////////////////////////////////////////////////
    // get and print exposures
    double iexp_on    = fTotalExposureOn;
    double iexp_off   = fTotalExposureOff;
    if( onrun != -1 )
    {
        iexp_on = fStereoOn->getEffectiveExposure( onrun );
        fRunExposureOn[onrun] = iexp_on;
        fTotalExposureOn += iexp_on;
        
        iexp_off = fStereoOff->getEffectiveExposure( offrun );
        fRunExposureOff[offrun] = iexp_off;
        fTotalExposureOff += iexp_off;
    }
    
    if( onrun != -1 && offrun != -1 )
    {
        cout << endl << "Mean properties for this pair: ON=" << onrun;
        cout << ", OFF=" << offrun << " -----------------------------" << endl;
    }
    cout << "\t Exposure ON=" << iexp_on << " secs (" << iexp_on / 60. << " min)";
    cout <<               "  OFF=" << iexp_off << " secs (" << iexp_off / 60. << " min),";
    cout << "\t Az range ON [" << fRunAzMinOn[onrun] << "," << fRunAzMaxOn[onrun] << "],";
    cout << " OFF [" << fRunAzMinOff[offrun] << "," << fRunAzMaxOff[offrun] << "]" << endl;
    
    ////////////////////////////////////////////////////////////
    // calculate and print mean elevation and raw rate
    if( onrun != -1 )
    {
        cout << "\t mean elevation: " << fStereoOn->getMeanElevation() << " (ON), " << fStereoOff->getMeanElevation() << " (OFF)" << endl;
        cout << "\t mean azimuth: " << fRunAzMeanOn[onrun] << " (ON), " << fRunAzMeanOff[offrun] << " (OFF)" << endl;
        if( fStereoOn->getMeanElevation() > 0. )
        {
            fMeanElevationOn += fStereoOn->getMeanElevation();
            fMeanAzimuthOn += fStereoOn->getMeanAzimuth();
        }
        else if( fStereoOff->getMeanElevation() > 0 )
        {
            fMeanElevationOn += fStereoOff->getMeanElevation();
            fMeanAzimuthOn += fStereoOff->getMeanAzimuth();
        }
        fMeanElevationOff += fStereoOff->getMeanElevation();
        fMeanAzimuthOff += fStereoOff->getMeanAzimuth();
        fNMeanElevation++;
        if( iexp_on > 0. && iexp_off > 0. )
        {
            cout << "\t trigger rate : " << fStereoOn->getRawRate() / iexp_on << " Hz (ON), ";
            cout << fStereoOff->getRawRate() / iexp_off << " Hz (Off)" << endl;
            fMeanRawRateOn += fStereoOn->getRawRate() / iexp_on;
            fMeanRawRateOff += fStereoOff->getRawRate() / iexp_off;
        }
        fMeanPedVarsOn += fRunPedVarsOn[onrun];
        fMeanPedVarsOff += fRunPedVarsOff[offrun];
    }
    else
    {
        if( fNMeanElevation > 0. )
        {
            fMeanElevationOn /= fNMeanElevation;
            fMeanElevationOff /= fNMeanElevation;
            fMeanAzimuthOn /= fNMeanElevation;
            if( fMeanAzimuthOn > 180. )
            {
                fMeanAzimuthOn -= 360.;
            }
            fMeanAzimuthOff /= fNMeanElevation;
            if( fMeanAzimuthOff > 180. )
            {
                fMeanAzimuthOff -= 360.;
            }
            fMeanRawRateOn /= fNMeanElevation;
            fMeanRawRateOff /= fNMeanElevation;
            fMeanPedVarsOn /= fNMeanElevation;
            fMeanPedVarsOff /= fNMeanElevation;
        }
        cout << "\t mean elevation: " << fMeanElevationOn << " (ON), " << fMeanElevationOff << " (OFF)" << endl;
        cout << "\t trigger rate : " << fMeanRawRateOn << " Hz (ON), " << fMeanRawRateOff << " Hz (Off)" << endl;
        cout << "\t mean pedvars: " << fMeanPedVarsOn << " (ON), " << fMeanPedVarsOff << " (OFF)" << endl;
    }
    
    ////////////////////////////////////////////////////////////
    // create alpha histogram for significance calculations
    // (called for correlated and uncorrelated histograms)
    fStereoOff->scaleAlpha( fStereoOn->getAlpha(), false );
    fStereoOff->scaleAlpha( fStereoOn->getAlphaUC(), true );
    
    ////////////////////////////////////////////////////////////
    // ON / OFF Analysis
    ////////////////////////////////////////////////////////////
    VOnOff* fstereo_onoff = new VOnOff();
    
    // normalization at target position
    double i_norm_alpha = fStereoOff->getAlphaNorm()->GetBinContent( fStereoOff->getAlphaNorm()->GetXaxis()->FindBin( -1.*fRunPara->fTargetShiftWest ),
                          fStereoOff->getAlphaNorm()->GetYaxis()->FindBin( -1.*fRunPara->fTargetShiftNorth ) );
                          
    // on-off for 1D histograms
    fstereo_onoff->doOnOffforParameterHistograms( fStereoOn->getParameterHistograms(), fStereoOff->getParameterHistograms(), i_norm_alpha, ( onrun == -1 ) );
    
    // on-off for correlated maps
    fstereo_onoff->doOnOffforSkyHistograms( fStereoOn->getSkyHistograms( false ), fStereoOff->getSkyHistograms( false ), fStereoOff->getAlphaNorm() );
    // on-off for uncorrelated maps
    fstereo_onoff->doOnOffforSkyHistograms( fStereoOn->getSkyHistograms( true ), fStereoOff->getSkyHistograms( true ), fStereoOff->getAlphaNormUC() );
    
    // Li & Ma significance maps and
    // print out maximum in maps
    cout << "\t Maximum in CORRELATED maps: " << endl;
    TH2D* hStSig = ( TH2D* )fstereo_onoff->do2DSignificance(
                       fStereoOn->getStereoSkyMap(),
                       fStereoOff->getStereoSkyMap(),
                       fStereoOff->getAlphaNorm() );
    cout << "\t Maximum in UNCORRELATED maps: " << endl;
    TH2D* hStSigUC = ( TH2D* )fstereo_onoff->do2DSignificance(
                         fStereoOn->getStereoSkyMapUC(),
                         fStereoOff->getStereoSkyMapUC(),
                         fStereoOff->getAlphaNormUC() );
                         
    ////////////////////////////////////////////////////////////
    // calulate significance in source bin
    
    // number of on events
    double i_nevts_on = fStereoOn->getStereoSkyMap()->GetBinContent( fStereoOn->getStereoSkyMap()->GetXaxis()->FindBin( -1.*fRunPara->fTargetShiftWest ),
                        fStereoOn->getStereoSkyMap()->GetYaxis()->FindBin( -1.*fRunPara->fTargetShiftNorth ) );
    // number of off events
    double i_nevts_off = fStereoOff->getStereoSkyMap()->GetBinContent( fStereoOff->getStereoSkyMap()->GetXaxis()->FindBin( -1.*fRunPara->fTargetShiftWest ),
                         fStereoOff->getStereoSkyMap()->GetYaxis()->FindBin( -1.*fRunPara->fTargetShiftNorth ) );
                         
    double i_sig = VStatistics::calcSignificance( i_nevts_on, i_nevts_off, i_norm_alpha );
    double i_rate = 0.;
    double i_rateE = 0.;
    double i_rateOFF = 0.;
    if( iexp_on > 0. && iexp_off > 0. )
    {
        i_rate = ( i_nevts_on - i_norm_alpha * i_nevts_off ) * 60. / iexp_on;       // rates in 1/min
        i_rateOFF = i_norm_alpha * i_nevts_off * 60. / iexp_off;                    // rates in 1/min
        i_rateE = sqrt( i_nevts_on + i_norm_alpha * i_norm_alpha * i_nevts_off ) * 60. / iexp_on;
    }
    
    cout << endl;
    cout << "\t ---------------------------- " << endl;
    cout << "\t RESULTS FOR SOURCE POSITION: " << endl;
    cout << "\t ---------------------------- " << endl;
    cout << "\t ON:" << i_nevts_on << "  OFF:" << setprecision( 4 ) << i_nevts_off* i_norm_alpha << " (";
    cout << "off " << i_nevts_off << ", alpha=" << i_norm_alpha << ")" << endl;
    cout << "\t " << setprecision( 4 ) <<  i_sig << " Sigma  " << i_rate << "+/-" << i_rateE << " gammas/min" << endl;
    cout << "\t background rate: " << i_rateOFF << " CR/min" << endl;
    cout << endl;
    
    // calculate q-factors
    fstereo_onoff->doQfactors( fStereoOn->getParameterHistograms(), fStereoOff->getParameterHistograms(), 1. );
    
    // fill run summary tree
    fillRunSummary( onrun, offrun, iexp_on, iexp_off, i_nevts_on, i_nevts_off, i_norm_alpha, i_sig, i_rate, i_rateOFF, fstereo_onoff );
    
    /////////////////////////////////////////////////////////
    // finalize and write everything to disk
    idir->cd();
    fstereo_onoff->fill1DSignificanceHistogram();
    fstereo_onoff->writeHistograms( hStSig, hStSigUC );
    fStereoOn->writeHistograms( true );
    fStereoOff->writeHistograms( false );
    if( onrun != -1 )
    {
        fStereoOff->writeDebugHistograms();
    }
    if( onrun == -1 )
    {
        // write run summary to disk
        fRunSummary->write();
    }
    // write list of excluded regions to disk (as a tree)
    fRunPara->writeListOfExcludedSkyRegions( onrun );
    
    // close data files, etc.
    fStereoOn->terminate();
    fStereoOff->terminate();
    
    // clean up
    delete fstereo_onoff;
}

/*!
 *   from http://root.cern.ch/phpBB2/viewtopic.php?t=2789
 *
 */
void VAnaSum::copyDirectory( TDirectory* source )
{
    //copy all objects and subdirs of directory source as a subdir of the current directory
    TDirectory* savdir = gDirectory;
    TDirectory* adir = savdir->mkdir( source->GetName() );
    if( !adir )
    {
        cout << "VAnaSum::copyDirectory error creating directory " << source->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    adir->cd();
    //loop on all entries of this directory
    TKey* key;
    TIter nextkey( source->GetListOfKeys() );
    while( ( key = ( TKey* )nextkey() ) )
    {
        const char* classname = key->GetClassName();
        TClass* cl = gROOT->GetClass( classname );
        if( !cl )
        {
            continue;
        }
        if( cl->InheritsFrom( "TDirectory" ) )
        {
            source->cd( key->GetName() );
            TDirectory* subdir = gDirectory;
            adir->cd();
            copyDirectory( subdir );
            adir->cd();
        }
        else if( cl->InheritsFrom( "TTree" ) )
        {
            if( !source->Get( "tRE" ) )
            {
                TTree* T = ( TTree* )source->Get( key->GetName() );
                adir->cd();
                if( T )
                {
                    TTree* newT = T->CloneTree();
                    if( newT )
                    {
                        newT->Write();
                        delete newT;
                    }
                    delete T;
                }
            }
        }
        else
        {
            source->cd();
            TObject* obj = key->ReadObj();
            if( obj )
            {
                adir->cd();
                obj->Write();
                delete obj;
            }
        }
    }
    adir->SaveSelf( kTRUE );
    savdir->cd();
}

/*!

   fill run summary tree

   - one entry per run
   - one entry for combined analysis

*/
void VAnaSum::fillRunSummary( int onrun, int offrun, double iexp_on, double iexp_off,
                              double i_nevts_on, double i_nevts_off, double i_norm_alpha,
                              double i_sig, double i_rate, double i_rateOFF, VOnOff* fstereo_onoff )
{
    if( !fRunSummary )
    {
        return;
    }
    
    // fill results tree
    fRunSummary->runOn = onrun;
    fRunSummary->runOff = offrun;
    if( onrun != -1 )
    {
        fRunSummary->MJDOn           = fStereoOn->getMJD( onrun );
        fRunSummary->MJDOn_runStart  = fStereoOn->getMJDStart( onrun );
        fRunSummary->MJDOn_runStopp  = fStereoOn->getMJDStopp( onrun );
        fRunSummary->RunDurationOn   = fStereoOn->getRunDuration( onrun );
        fRunSummary->MJDOff          = fStereoOff->getMJD( offrun );
        fRunSummary->MJDOff_runStart = fStereoOn->getMJDStart( offrun );
        fRunSummary->MJDOff_runStopp = fStereoOn->getMJDStopp( offrun );
        fRunSummary->RunDurationOff  = fStereoOn->getRunDuration( offrun );
        fRunSummary->fWobbleNorth = fStereoOn->getWobbleNorth();
        fRunSummary->fWobbleWest = fStereoOn->getWobbleWest();
        fRunSummary->fNTel = fRunPara->fMapRunList[onrun].fNTel;
        fRunSummary->elevationOn = fStereoOn->getMeanElevation();
        fRunSummary->elevationOff = fStereoOff->getMeanElevation();
        fRunSummary->azimuthOn = fRunAzMeanOn[onrun];
        fRunSummary->azimuthOff = fRunAzMeanOff[offrun];
        for( unsigned int p = 0; p < fRunPara->fRunList.size(); p++ )
        {
            stringstream iTelCombination;
            if( fRunPara->fRunList[p].fRunOn == onrun )
            {
                for( unsigned int t = 0; t < fRunPara->fRunList[p].fTelToAnalyze.size(); t++ )
                {
                    iTelCombination << "T";
                    iTelCombination << fRunPara->fRunList[p].fTelToAnalyze[t] + 1;
                }
                if( iTelCombination.str().size() < 300 )
                {
                    sprintf( fRunSummary->fTelList, "%s", iTelCombination.str().c_str() );
                }
                else
                {
                    sprintf( fRunSummary->fTelList, "%s", iTelCombination.str().substr( 0, 299 ).c_str() );
                }
            }
        }
        fRunSummary->fTheta2Max = fRunPara->fMapRunList[onrun].fSourceRadius;
    }
    else
    {
        fRunSummary->MJDOn = 0.;    //Could make this the mean MJD of all ON runs included in the summed analysis
        fRunSummary->MJDOn_runStart = 0.;
        fRunSummary->MJDOn_runStopp = 0.;
        fRunSummary->RunDurationOn  = 0.;
        fRunSummary->MJDOff = 0.;
        fRunSummary->MJDOff_runStart = 0.;
        fRunSummary->MJDOff_runStopp = 0.;
        fRunSummary->RunDurationOff  = 0.;
        fRunSummary->fWobbleNorth = 0.;
        fRunSummary->fWobbleWest = 0.;
        fRunSummary->fNTel = 0;
        sprintf( fRunSummary->fTelList, "NOTSET" );
        fRunSummary->elevationOn = fMeanElevationOn;
        fRunSummary->elevationOff = fMeanElevationOff;
        fRunSummary->azimuthOn = fMeanAzimuthOn;
        fRunSummary->azimuthOff = fMeanAzimuthOff;
    }
    if( onrun != -1 )
    {
        if( fRunPara->fMapRunList[onrun].fTarget.size() < 300 )
        {
            sprintf( fRunSummary->fTargetName, "%s", fRunPara->fMapRunList[onrun].fTarget.c_str() );
        }
        else
        {
            sprintf( fRunSummary->fTargetName, "%s", fRunPara->fMapRunList[onrun].fTarget.substr( 0, 299 ).c_str() );
        }
    }
    if( onrun != -1 && fRunPara->fMapRunList.find( onrun ) != fRunPara->fMapRunList.end() )
    {
        fRunSummary->fTargetRA = fRunPara->fMapRunList[onrun].fTargetRA;
        fRunSummary->fTargetDec = fRunPara->fMapRunList[onrun].fTargetDec;
        fRunSummary->fTargetRAJ2000 = fRunPara->fMapRunList[onrun].fTargetRAJ2000;
        fRunSummary->fTargetDecJ2000 = fRunPara->fMapRunList[onrun].fTargetDecJ2000;
        fRunSummary->fTargetShiftWest = fRunPara->fMapRunList[onrun].fTargetShiftWest;
        fRunSummary->fTargetShiftNorth = fRunPara->fMapRunList[onrun].fTargetShiftNorth;
    }
    else
    {
        fRunSummary->fTargetRA = 0.;
        fRunSummary->fTargetDec = 0.;
        fRunSummary->fTargetRAJ2000 = 0.;
        fRunSummary->fTargetDecJ2000 = 0.;
        fRunSummary->fTargetShiftWest = 0.;
        fRunSummary->fTargetShiftNorth = 0.;
    }
    fRunSummary->fSkyMapCentreRAJ2000 = fRunPara->fSkyMapCentreRAJ2000;
    fRunSummary->fSkyMapCentreDecJ2000 = fRunPara->fSkyMapCentreDecJ2000;
    fRunSummary->fTargetShiftRAJ2000 = fRunPara->fTargetShiftRAJ2000;
    fRunSummary->fTargetShiftDecJ2000 = fRunPara->fTargetShiftDecJ2000;
    fRunSummary->tOn = iexp_on;
    fRunSummary->tOff = iexp_off;
    if( onrun != -1 )
    {
        if( iexp_on > 0. )
        {
            fRunSummary->RawRateOn = fStereoOn->getRawRate() / iexp_on;
        }
        else
        {
            fRunSummary->RawRateOn = 0.;
        }
        if( iexp_off > 0. )
        {
            fRunSummary->RawRateOff = fStereoOff->getRawRate() / iexp_off;
        }
        else
        {
            fRunSummary->RawRateOff = 0.;
        }
    }
    else
    {
        fRunSummary->RawRateOn = fMeanRawRateOn;
        fRunSummary->RawRateOff = fMeanRawRateOff;
    }
    if( onrun != -1 )
    {
        fRunSummary->pedvarsOn = fRunPedVarsOn[onrun];
        fRunSummary->pedvarsOff = fRunPedVarsOff[offrun];
    }
    else
    {
        fRunSummary->pedvarsOn = fMeanPedVarsOn;
        fRunSummary->pedvarsOff = fMeanPedVarsOff;
    }
    fRunSummary->NOn = i_nevts_on;
    fRunSummary->NOff = i_nevts_off;
    fRunSummary->NOffNorm = i_nevts_off * i_norm_alpha;
    fRunSummary->OffNorm = i_norm_alpha;
    fRunSummary->Signi = i_sig;
    fRunSummary->Rate = i_rate;
    double i_tnorm = 1.;
    if( fRunSummary->tOff > 0. )
    {
        i_tnorm = fRunSummary->tOn / fRunSummary->tOff;
    }
    if( fRunSummary->tOn > 0. )
    {
        fRunSummary->RateE = sqrt( i_nevts_on + i_tnorm * i_tnorm * i_norm_alpha * i_norm_alpha * i_nevts_off ) / fRunSummary->tOn * 60.;
    }
    else
    {
        fRunSummary->RateE = 0.;
    }
    fRunSummary->RateOff = i_rateOFF;
    if( fRunSummary->tOff > 0. )
    {
        fRunSummary->RateOffE = sqrt( i_tnorm * i_tnorm * i_norm_alpha * i_norm_alpha * i_nevts_off ) / fRunSummary->tOff * 60.;
    }
    else
    {
        fRunSummary->RateOffE = 0.;
    }
    if( onrun != -1 )
    {
        fRunSummary->DeadTimeFracOn = fStereoOn->getDeadTimeFraction();
        fMeanDeadTimeOn += fStereoOn->getDeadTimeFraction() * iexp_on;
        fRunSummary->DeadTimeFracOff = fStereoOff->getDeadTimeFraction();
        fMeanDeadTimeOff += fStereoOff->getDeadTimeFraction() * iexp_off;
    }
    else
    {
        if( iexp_on > 0. )
        {
            fRunSummary->DeadTimeFracOn = fMeanDeadTimeOn / iexp_on;
        }
        if( iexp_off > 0. )
        {
            fRunSummary->DeadTimeFracOff = fMeanDeadTimeOff / iexp_off;
        }
    }
    if( fstereo_onoff )
    {
        fRunSummary->MaxSigni = fstereo_onoff->getMaxSigma();
        fRunSummary->MaxSigniX = fstereo_onoff->getMaxSigmaX();
        fRunSummary->MaxSigniY = fstereo_onoff->getMaxSigmaY();
    }
    fRunSummary->fill();
}


void VAnaSum::terminate()
{
    if( fOPfile )
    {
        fOPfile->Close();
    }
}


/*!
 *
 * get azimuth range of a run
 *
 *
 */
double VAnaSum::getAzRange( int i_run, string i_treename, double& azmin, double& azmax )
{
    double azmean = 0.;
    
    azmin = 1.e3;
    azmax = -1.e3;
    
    char i_temp[2000];
    sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), i_run, fSuffix.c_str() );
    TFile* i_f = new TFile( i_temp );
    if( i_f->IsZombie() )
    {
        cout << "VAnaSum::getAZRange fatal error: file not found, " << i_temp << endl;
        exit( EXIT_FAILURE );
    }
    
    TTree* i_tree = ( TTree* )i_f->Get( i_treename.c_str() );
    if( !i_tree )
    {
        cout << "VAnaSum::getAZRange tree not found " << i_treename << endl;
        exit( EXIT_FAILURE );
    }
    Float_t ArrayPointing_Azimuth;
    i_tree->SetBranchAddress( "ArrayPointing_Azimuth", &ArrayPointing_Azimuth );
    
    i_tree->GetEntry( 0 );
    azmin = ArrayPointing_Azimuth;
    i_tree->GetEntry( i_tree->GetEntries() - 1 );
    azmax = ArrayPointing_Azimuth;
    
    if( azmin > 180. )
    {
        azmin -= 360.;
    }
    if( azmax > 180. )
    {
        azmax -= 360.;
    }
    cout << "\t azimuth range: [" << azmin << "," << azmax << "]" << endl;
    
    // calculate mean az
    // mean azimuth angle
    if( azmin > 120. && azmax < -120. )
    {
        azmax += 360.;
    }
    else if( azmin < -150. && azmax > 120. )
    {
        azmin += 360.;
    }
    
    azmean = 0.5 * ( azmin + azmax );
    if( azmean > 180. )
    {
        azmean -= 360.;
    }
    
    i_f->Close();
    delete i_f;
    
    return azmean;
}

/*!

    get mean noise level for a run

*/
double VAnaSum::getNoiseLevel( int i_run )
{
    char i_temp[2000];
    double ipedv = -1.;
    
    sprintf( i_temp, "%s%s%d%s", fDatadir.c_str(), fPrefix.c_str(), i_run, fSuffix.c_str() );
    TFile* i_f = new TFile( i_temp );
    if( i_f->IsZombie() )
    {
        cout << "VAnaSum::getNoiseLevel fatal error: file not found, " << i_temp << endl;
        exit( EXIT_FAILURE );
    }
    VTableLookupRunParameter* fR = ( VTableLookupRunParameter* )i_f->Get( "TLRunParameter" );
    if( fR )
    {
        ipedv = fR->meanpedvars;
    }
    else
    {
        ipedv = -1.;
    }
    i_f->Close();
    cout << "\t mean pedestal variations in run " << i_run << ": " << ipedv << endl;
    return ipedv;
}

