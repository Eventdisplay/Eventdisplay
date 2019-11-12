/*! \class VDST
    \brief writes data summary files with sums and times for each pixel

    output is after pedestal substraction, gain and toffset correction


*/

#include <VDST.h>

VDST::VDST( bool iMode, bool iMC )
{
    if( getDebugFlag() )
    {
        cout << "VDST::VDST()" << endl;
    }
    // laser run as source file? (Hardcoded, should be in VEvndispRunParameter)
    fBLaser = false;
    fRaw = false;
    
    // initalize flag (true after first event)
    fDSTini = false;
    fDSTfile = 0;
    // no dst output, don't do anything
    if( !iMode )
    {
        return;
    }
    // open dst output file
    fDSTfile = new TFile( getRunParameter()->fdstfile.c_str(), "RECREATE" );
    if( fDSTfile->IsZombie() )
    {
        cout << "VDST::VDST error while create dst file" << endl;
        exit( EXIT_FAILURE );
    }
    
    initDSTTree( true, getRunParameter()->fdstcalibration );
    setMC( iMC );
    
    fVImageCleaning = new VImageCleaning( getData() );
}

VDST::~VDST()
{
    if( fDSTfile && !fDSTfile->IsZombie() )
    {
        fDSTfile->Close();
    }
    if( fVImageCleaning )
    {
        delete fVImageCleaning;
    }
}


void VDST::initialize()
{
    initializeDataReader();
}

/*
    analyse this event and fill the DST tree
*/
void VDST::fill()
{
    if( getDebugFlag() )
    {
        cout << "VDST::fill()" << endl;
    }
    if( fDST_tree == 0 )
    {
        return;
    }
    // init base analyzer
    if( !fDSTini )
    {
        for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
        {
            setTelID( getTeltoAna()[i] );
            // set summation vectors
            setAnaData();
            // find dead channels
            findDeadChans( false, true );
            findDeadChans( true, true );
            
            // set special channels
            setSpecialChannels();
        }
        fDSTini = true;
    }
    // get geometry
    fDSTntel = getNTel();
    // test geometry
    if( fDSTntel > VDST_MAXTELESCOPES || ( int )getNChannels() > VDST_MAXCHANNELS )
    {
        cout << "void VDST::fill(): error: too many cameras or pmts, maximum is " << VDST_MAXTELESCOPES << "/" << VDST_MAXCHANNELS << endl;
        exit( -1 );
    }
    // ignore following code
    // (GM) XXXX check muon/laser list
    /*     int iZ = 0;
         for( unsigned int k = 0; k < vEventList.size(); k++ )
         {
            iZ++;
            if( vEventList[k] != getEventNumber() ) continue;
        else break;
         }
         if( iZ > 12 ) return;  */
    // (GM) XXXX end check muon/laser list
    // get basic event data
    fDSTrunnumber = getRunNumber();
    fDSTeventnumber = getEventNumber();
    fDSTeventtype = ( unsigned int )( getReader()->getATEventType() );
    fDSTgps0 = getReader()->getGPS0();
    fDSTgps1 = getReader()->getGPS1();
    fDSTgps2 = getReader()->getGPS2();
    fDSTgps3 = getReader()->getGPS3();
    fDSTgps4 = getReader()->getGPS4();
    fDSTgpsyear = ( unsigned int )getReader()->getATGPSYear() + 2000;
    fDSTATgpsyear = ( unsigned int )getReader()->getATGPSYear();
    // update pointing
    // get pointing and local trigger (stored as unsigned long)
    bitset<8 * sizeof( unsigned long ) > i_localTrigger;
    fDSTNTrig = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( isTeltoAna( i ) )
        {
            // local trigger
            if( fReader->hasLocalTrigger( i ) )
            {
                i_localTrigger.set( i, 1 );
                fDSTLTrig_list[fDSTNTrig] = i;
                fDSTNTrig++;
            }
            else
            {
                fDSTLTrig_list[i] = 0;
            }
            // pointing
            getPointing()[i]->setTelPointing( getEventMJD(), getEventTime() );
            if( i < fReader->getTelAzimuth().size() )
            {
                fDSTpointAzimuth[i]   = fReader->getTelAzimuth()[i];
            }
            else if( i < getPointing().size() && getPointing()[i] )
            {
                fDSTpointAzimuth[i]   = getPointing()[i]->getTelAzimuth();
            }
            if( i < fReader->getTelElevation().size() )
            {
                fDSTpointElevation[i] = fReader->getTelElevation()[i];
            }
            else if( i < getPointing().size() && getPointing()[i] )
            {
                fDSTpointElevation[i]   = getPointing()[i]->getTelElevation();
            }
        }
    }
    fDSTLTrig = i_localTrigger.to_ulong();
    // fill DST tree only for triggered events
    if( fMC && getRunParameter()->fIsMC == 2 && fDSTLTrig == 0 )
    {
        return;
    }
    // get the Monte Carlo data
    if( fMC )
    {
        fDSTprimary = fReader->getMC_primary();
        fDSTenergy = ( float )fReader->getMC_energy();
        fDSTxcore = ( float )fReader->getMC_X();
        fDSTycore = ( float )fReader->getMC_Y();
        fDSTze = ( float )fReader->getMC_Ze();
        fDSTaz = ( float )fReader->getMC_Az();
        fDSTTel_xoff = ( float )fReader->getMC_Xoffset();
        fDSTTel_yoff = ( float )fReader->getMC_Yoffset();
    }
    // check if there was an array trigger, otherwise skip event
    if( fReader->isMC() && fDSTLTrig == 0 )
    {
        resetDataVectors();
        fDST_tree->Fill();
        return;
    }
    
    // temporary variable
    vector< float >  t_PulseTimingTemp;
    
    // ntel is always to total number of telescopes in the DST file
    fDSTntel_data = getNTel();
    for( unsigned int i = 0; i < fDSTntel_data; i++ )
    {
        fDSTtel_data[i] = i;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loop over all telescope and analyse channel by channel
    // get sums and toffsets
    int intubes = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        double i_total = 0;
        if( isTeltoAna( i ) )
        {
            intubes = 0;
            setTelID( i );
            
            // set number of samples
            setNSamples( fReader->getNumSamples() );
            // set hilo
            fillHiLo();
            // find dead channels (time dependent)
            //	findDeadChans( false, false );
            //	findDeadChans( true, false );
            // integrate FADC traces -> calculate integrated charges and pulse timing
            calcTZerosSums( getSumFirst(), getSumFirst() + getSumWindow(), getTraceIntegrationMethod() );
            // image cleaning if image threshold is > 0.
            if( getImageCleaningParameter() && getImageCleaningParameter()->fimagethresh > 0. )
            {
                if( fVImageCleaning )
                {
                    fVImageCleaning->cleanImagePedvars( getImageCleaningParameter() );
                }
                if( i < getRunParameter()->fDoublePass.size() && getRunParameter()->fDoublePass[i] )
                {
                    calcSecondTZerosSums();
                    if( fVImageCleaning )
                    {
                        fVImageCleaning->cleanImagePedvars( getImageCleaningParameter() );
                    }
                }
            }
            else
            {
                setImage( true );
                setBorder( false );
            }
            // calculate total sum over all channels
            for( unsigned int j = 0; j < getNChannels(); j++ )
            {
                i_total += getSums()[j];
            }
        }
        // fill pulse timing level
        if( fRunPar->fpulsetiminglevels.size() < getDSTpulsetiminglevelsN() )
        {
            for( unsigned int t = 0; t < fRunPar->fpulsetiminglevels.size(); t++ )
            {
                fDSTpulsetiminglevels[i][t] = fRunPar->fpulsetiminglevels[t];
            }
        }
        
        // fill dst arrays
        for( unsigned int j = 0; j < getNChannels(); j++ )
        {
            fDSTChan[i][j] = j;
            
            // fill values for this event, this telescope and this pixel if:
            //    i) this is a valid laser event
            //           or
            //   ii) fRunPar->fdstwriteallpixel is set to true
            //           or
            //  iii) this is a border or image pixel
            if( isTeltoAna( i ) && ( ( !fBLaser || ( i_total > fRunPar->fLaserSumMin ) )
                                     && ( fRunPar->fdstwriteallpixel || getBorder()[j] || getImage()[j] ) ) )
            {
                intubes++;
                ///////////////////////////////////////////////
                // for laser only
                if( fBLaser )
                {
                    int corrfirst = TMath::Nint( getTZeros()[j] ) - 3;
                    
                    fDSTpedestal[i][j] = ( float )getPeds( getHiLo()[j] )[j];
                    fDSTsums[i][j] = ( float )fTraceHandler->getTraceSum( corrfirst, corrfirst + getSumWindow(), false );
                    fDSTsums2[i][j] = fDSTsums[i][j];
                    // ignore dead low gain channels
                    fDSTdead[i][j] = ( unsigned int )getDead()[j];
                    fDSTsumwindow[i][j] = getCurrentSumWindow()[j];
                    fDSTsumfirst[i][j] = getTCorrectedSumFirst()[j];
                    // fill pulse timing
                    if( fRunPar->fpulsetiminglevels.size() < getDSTpulsetiminglevelsN() )
                    {
                        t_PulseTimingTemp = fTraceHandler->getPulseTiming( corrfirst, corrfirst + getSumWindow(), 0, getNSamples() );
                        for( unsigned int t = 0; t < fRunPar->fpulsetiminglevels.size(); t++ )
                        {
                            fDSTpulsetiming[i][t][j] = t_PulseTimingTemp[t];
                        }
                    }
                    double i_max = 0.;
                    unsigned int maxpos = 0;
                    unsigned int n255 = 0;
                    fTraceHandler->getTraceMax( corrfirst, corrfirst + getSumWindow(), i_max, maxpos, n255 );
                    if( maxpos != 99999 )
                    {
                        fDSTMax[i][j] = ( short )i_max;
                        fDSTRawMax[i][j] = ( short )( i_max + getPeds( getHiLo()[j] )[j] );
                    }
                    else
                    {
                        fDSTMax[i][j] = 0;
                        fDSTRawMax[i][j] = 0;
                    }
                }
                ///////////////////////////////////////////////
                // normal data
                else
                {
                    fDSTpedestal[i][j] = ( float )getPeds( getHiLo()[j] )[j];
                    fDSTsums[i][j] = ( float )getSums()[j];
                    fDSTsums2[i][j] = ( float )getSums2()[j];
                    //set channel dead if it is ( dead in low gain AND in low gain ) OR (dead in high gain AND in high gain )
                    fDSTdead[i][j] = ( unsigned int )getDead( getHiLo()[j] )[j];
                    fDSTsumwindow[i][j] = getCurrentSumWindow()[j];
                    fDSTsumfirst[i][j] = getTCorrectedSumFirst()[j];
                    if( fRunPar->fpulsetiminglevels.size() < getDSTpulsetiminglevelsN() )
                    {
                        for( unsigned int t = 0; t < fRunPar->fpulsetiminglevels.size(); t++ )
                        {
                            fDSTpulsetiming[i][t][j] = ( float )getPulseTiming()[t][j];
                        }
                    }
                    fDSTTraceWidth[i][j] = ( float )getTraceWidth()[j];
                    fDSTMax[i][j] = ( short )getTraceMax()[j];
                    fDSTN255[i][j] = getTraceN255()[j];
                    fDSTRawMax[i][j] = ( short )getTraceRawMax()[j];
                }
                fDSTHiLo[i][j] = ( unsigned int )getHiLo()[j];
            }
            else
            {
                fDSTpedestal[i][j] = 0.;
                fDSTsums[i][j] = 0.;
                fDSTsums2[i][j] = 0.;
                fDSTsumwindow[i][j] = 0;
                fDSTsumfirst[i][j] = 0;
                if( fRunPar->fpulsetiminglevels.size() < getDSTpulsetiminglevelsN() )
                {
                    for( unsigned int t = 0; t < fRunPar->fpulsetiminglevels.size(); t++ )
                    {
                        fDSTpulsetiming[i][t][j] = 0.;
                    }
                }
                fDSTMax[i][j] = 0;
                fDSTRawMax[i][j] = 0;
                fDSTHiLo[i][j] = 0;
                fDSTN255[i][j] = 0;
            }
        }
        if( fMC )
        {
            fDSTLTtime[i] = ( float )getReader()->getLocalTriggerTime( i );
            fDSTLDTtime[i] = ( float )getReader()->getLocalDelayedTriggerTime( i );
            fDSTL2TrigType[i] = ( int )getReader()->getLocalTriggerType( i );
        }
    }
    // only write events with more than NN ntubes to disk (for one telescope only)
    if( getNTel() == 1 && intubes < getRunParameter()->fdstminntubes )
    {
        return;
    }
    fDST_tree->Fill();
}


void VDST::terminate()
{
    if( getDebugFlag() )
    {
        cout << "VDST::terminate()" << endl;
    }
    // now write everything to disk
    if( fDSTfile )
    {
        if( fDSTfile->cd() )
        {
            cout << "writing data summary tree to " << fDSTfile->GetName() << endl;
            cout << "\t total number of events in dst tree: " << fDST_tree->GetEntries() << endl;
            fDST_tree->Write();
            getRunParameter()->Write();
            // write detector tree
            if( getDetectorTree() )
            {
                if( fDebug )
                {
                    cout << "\t writing detector tree: " << getDetectorTree()->GetName() << endl;
                }
                getDetectorTree()->Write();
            }
            if( fDebug )
            {
                cout << "\t writing mean pulses " << endl;
            }
            // write pulse shape histograms
            if( fDSTfile->mkdir( "meanPulses" )->cd() )
            {
                for( unsigned int i = 0; i < getNTel(); i++ )
                {
                    setTelID( i );
                    if( getMeanPulseHistograms() )
                    {
                        getMeanPulseHistograms()->Write();
                    }
                }
            }
            // write pulse sum histograms
            if( fDSTfile->mkdir( "pulseSums" )->cd() )
            {
                for( unsigned int i = 0; i < getNTel(); i++ )
                {
                    setTelID( i );
                    if( getIntegratedChargeHistograms() )
                    {
                        getIntegratedChargeHistograms()->Write();
                    }
                }
            }
            // write calibration data
            writeCalibrationData();
        }
        fDSTfile->Close();
    }
}

bool VDST::writeCalibrationData()
{
    if( fDebug )
    {
        cout << "VDST::writeCalibrationData()" << endl;
    }
    if( fDSTfile )
    {
        fDSTfile->cd();
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // DST analysis
    // (same code as in c_VDST)

    // max channels very different for CTA and VERITAS; adjust
    float *fPedvar_high = 0;
    float *fPedvar_low = 0;
    if( fRunPar->getObservatory() == "VERITAS" )
    {
         fPedvar_high = new float[500*120];
         fPedvar_low= new float[500*120];
    }
    else
    {
         fPedvar_high = new float[VDST_MAXCHANNELS * VDST_MAXSUMWINDOW];
         fPedvar_low= new float[VDST_MAXCHANNELS * VDST_MAXSUMWINDOW];
    }

    int fTelID = 0;
    unsigned int nPixel = 0;
    // integration window
    unsigned int fnum_sumwindow = 1;
    unsigned int fsumwindow[VDST_MAXSUMWINDOW];
    float fPed_high[VDST_MAXCHANNELS];
    float fPed_low[VDST_MAXCHANNELS];
    float fConv_high[VDST_MAXCHANNELS];
    float fConv_low[VDST_MAXCHANNELS];
    float fTZero[VDST_MAXCHANNELS];
    
    for( unsigned int i = 0; i < VDST_MAXCHANNELS; i++ )
    {
        fPed_high[i] = 0.;
        fPed_low[i] = 0.;
        fConv_high[i] = 0.;
        fConv_low[i] = 0.;
        fTZero[i] = -999.;
        for( unsigned int j = 0; j < VDST_MAXSUMWINDOW; j++ )
        {
            fPedvar_high[i * VDST_MAXSUMWINDOW + j] = 0;
            fPedvar_low[i * VDST_MAXSUMWINDOW + j] = 0;
        }
    }
    
    TTree* t = new TTree( "calibration", "calibration data" );
    
    char hname[200];
    t->Branch( "TelID", &fTelID, "TelID/I" );
    t->Branch( "NPixel", &nPixel, "NPixel/i" );
    if( t->GetBranchStatus( "num_sumwindow" ) )
    {
        t->Branch( "num_sumwindow", &fnum_sumwindow, "num_sumwindow/i" );
        t->Branch( "sumwindow", fsumwindow, "sumwindow[num_sumwindow]/i" );
    }
    else
    {
        fnum_sumwindow = 0;
    }
    t->Branch( "ped_high", fPed_high, "ped_high[NPixel]/F" );
    sprintf( hname, "pedvar_high[%d]/F", VDST_MAXCHANNELS * VDST_MAXSUMWINDOW );
    t->Branch( "pedvar_high", fPedvar_high, hname );
    t->Branch( "ped_low", fPed_low, "ped_low[NPixel]/F" );
    sprintf( hname, "pedvar_low[%d]/F", VDST_MAXCHANNELS * VDST_MAXSUMWINDOW );
     t->Branch( "pedvar_low", fPedvar_low, hname );
    t->Branch( "conv_high", fConv_high, "conv_high[NPixel]/F" );
    t->Branch( "conv_low", fConv_low, "conv_low[NPixel]/F" );
    t->Branch( "tzero", fTZero, "tzero[NPixel]/F" );

    fnum_sumwindow = getRunParameter()->fCalibrationSumWindow;
    for( unsigned int i = 0; i < ( unsigned int )getRunParameter()->fCalibrationSumWindow; i++ )
    {
        fsumwindow[i] = i + 1;
    }
    
    for( unsigned int itel = 0; itel <  getNTel(); itel++ )
    {
        setTelID( itel );
        fTelID = getTelID();
        
        // correct number of samples
        if( getNSamples() < fnum_sumwindow )
        {
            fnum_sumwindow = getNSamples();
        }
        
        nPixel = ( unsigned int )getNChannels();
        if( VDST_MAXCHANNELS < nPixel )
        {
            cout << "DST_fillCalibrationTree error: number of pixels (" << nPixel << ") exeeds allowed range (" << VDST_MAXCHANNELS << ")" << endl;
            cout << "\t adjust arrays..." << endl;
            return false;
        }
        for( unsigned int p = 0; p < nPixel; p++ )
        {
            fPed_high[p] = getPeds()[p];
            
            for( unsigned int i = 0; i < ( unsigned int )fnum_sumwindow; i++ )
            {
                fPedvar_high[p * VDST_MAXSUMWINDOW + i] = getPedvars( i + 1 )[p];
                fPedvar_low[p * VDST_MAXSUMWINDOW + i] = getPedvars( i + 1, true )[p];
            }
            fPed_low[p] = getPeds( true )[p];
            fConv_high[p] = 1.;
            fConv_low[p] = 1.;
            if( p < getAverageTZeros().size() )
            {
                fTZero[p] = getAverageTZeros()[p];
            }
            else
            {
               fTZero[p] = 0.;
            }
        }
        
        t->Fill();
    }
    t->Write();
    
    if( fDebug )
    {
        cout << "END VDST::writeCalibrationData()" << endl;
    }

    delete [] fPedvar_high;
    delete [] fPedvar_low;
    
    return true;
}
