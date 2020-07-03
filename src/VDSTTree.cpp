/*! \class VDSTTree
    \brief writes data summary files with sums and times for each pixel

    output is after pedestal substraction, gain and toffset correction


*/

#include <VDSTTree.h>

VDSTTree::VDSTTree()
{
    // source data monte carlo?
    fMC = true;
    fFullTree = false;
    fDST_ADC_set = false;
    
    setFADC();
    setFillPELeaf();
    setFillPeakADFLeaf();
    
    fMCtree = 0;
    fDST_tree = 0;
    fDST_conf = 0;
    
    // initialize
    fTelescopeCounter_temp = -1;
    fDSTLTrig = 0;
    fDSTNTrig = 0;
    fDSTgpsyear = 0;
    fDSTATgpsyear = 0;
    fDSTgps0 = 0;
    fDSTgps1 = 0;
    fDSTgps2 = 0;
    fDSTgps3 = 0;
    fDSTgps4 = 0;
    fDSTntel = 0;
    fDSTntel_data = 0;
    
    fDSTrunnumber = 0;
    fDSTeventnumber = 0;
    fDSTeventtype = 0;
    fDSTprimary = 0;
    fDSTenergy = 0.;
    fDSTxcore = 0.;
    fDSTycore = 0.;
    fDSTze = 0;
    fDSTaz = 0;
    fDSTTel_xoff = 0;
    fDSTTel_yoff = 0;
    
    fDSTMeanPulseTimingMinLightLevel = 5.;
    for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
    {
        fDSTMeanPulseTimingHistogram[i] = 0;
        for( unsigned int j = 0; j < VDST_MAXCHANNELS; j++ )
        {
            fDSTMeanPulseTiming[i][j] = 0.;
            fDSTMeanPulseTiming_N[i][j] = 0.;
        }
    }

    resetDataVectors();
    
}


bool VDSTTree::initMCTree()
{
    fMCtree = new TTree( "mc", "mc events" );
    fMCtree->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    fMCtree->Branch( "runNumber", &fDSTrunnumber, "runNumber/i" );
    fMCtree->Branch( "eventNumber", &fDSTeventnumber, "eventNumber/i" );
    fMCtree->Branch( "MCprim", &fDSTprimary, "MCprimary/s" );
    fMCtree->Branch( "MCe0", &fDSTenergy, "MCenergy/F" );
    fMCtree->Branch( "MCxcore", &fDSTxcore, "MCxcore/F" );
    fMCtree->Branch( "MCycore", &fDSTycore, "MCycore/F" );
    fMCtree->Branch( "MCze", &fDSTze, "MCze/F" );
    fMCtree->Branch( "MCaz", &fDSTaz, "MCaz/F" );
    fMCtree->Branch( "MCxoff", &fDSTTel_xoff, "MCxoff/F" );
    fMCtree->Branch( "MCyoff", &fDSTTel_yoff, "MCyoff/F" );
    
    return true;
}

/*
   init DST tree for writing
*/
bool VDSTTree::initDSTTree( bool iFullTree, bool iCalibrationTree )
{
    char tname[1000];
    // don't allow full and calibration mode simultaneously
    if( iCalibrationTree )
    {
        iFullTree = false;
    }
    
    // DST tree definition
    fDST_tree = new TTree( "dst", "data summary tree" );
    fDST_tree->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    fDST_tree->Branch( "runNumber", &fDSTrunnumber, "runNumber/i" );
    fDST_tree->Branch( "eventNumber", &fDSTeventnumber, "eventNumber/i" );
    fDST_tree->Branch( "eventtype", &fDSTeventtype, "eventtype/i" );
    if( iFullTree && !iCalibrationTree )
    {
        fDST_tree->Branch( "gpsyear", &fDSTgpsyear, "gpsyear/i" );
        fDST_tree->Branch( "ATgpsyear", &fDSTATgpsyear, "ATgpsyear/i" );
        fDST_tree->Branch( "gps0", &fDSTgps0, "gps0/i" );
        fDST_tree->Branch( "gps1", &fDSTgps1, "gps1/i" );
        fDST_tree->Branch( "gps2", &fDSTgps2, "gps2/i" );
        fDST_tree->Branch( "gps3", &fDSTgps3, "gps3/i" );
        fDST_tree->Branch( "gps4", &fDSTgps4, "gps4/i" );
    }
    fDST_tree->Branch( "ntel", &fDSTntel, "ntel/i" );
    if( !iCalibrationTree )
    {
        // all following array are filled for all telescopes
        fDST_tree->Branch( "Paz", &fDSTpointAzimuth, "Paz[ntel]/F" );
        fDST_tree->Branch( "Pel", fDSTpointElevation, "Pel[ntel]/F" );
        // following arrays are filled only for all triggered telescopes
        fDST_tree->Branch( "ltrig", &fDSTLTrig, "ltrig/i" );
        fDST_tree->Branch( "ntrig", &fDSTNTrig, "ntrig/i" );
        fDST_tree->Branch( "ltrig_list", fDSTLTrig_list, "ltrig_list[ntrig]/i" );
        if( fMC )
        {
            fDST_tree->Branch( "LTtime", fDSTLTtime, "LTtime[ntrig]/F" );
            if( iFullTree )
            {
                fDST_tree->Branch( "LDTtime", fDSTLDTtime, "LDTtime[ntrig]/F" );
            }
            fDST_tree->Branch( "L2TrigType", fDSTL2TrigType, "L2TrigType[ntrig]/s" );
        }
    }
    fDST_tree->Branch( "ntel_data", &fDSTntel_data, "ntel_data/i" );
    // following arrays are filed only for all telescopes with data
    if( !iCalibrationTree )
    {
        fDST_tree->Branch( "tel_data", fDSTtel_data, "tel_data[ntel_data]/i" );
        fDST_tree->Branch( "tel_zero_suppression", fDSTTelescopeZeroSupression, "tel_zero_suppression[ntel_data]/s" );
        sprintf( tname, "chan[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "chan", fDSTChan, tname );
        sprintf( tname, "recorded[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "recorded", fDSTRecord, tname );
        fDST_tree->Branch( "nL1trig", fDSTnL1trig, "nL1trig[ntel_data]/s" );
        sprintf( tname, "L1trig[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "L1trig", fDSTL1trig, tname );
    }
    sprintf( tname, "ped[ntel_data][%d]/F", VDST_MAXCHANNELS );
    fDST_tree->Branch( "ped", fDSTpedestal, tname );
    sprintf( tname, "sum[ntel_data][%d]/F", VDST_MAXCHANNELS );
    fDST_tree->Branch( "sum", fDSTsums, tname );
    sprintf( tname, "sum2[ntel_data][%d]/F", VDST_MAXCHANNELS );
    fDST_tree->Branch( "sum2", fDSTsums2, tname );
    sprintf( tname, "dead[ntel_data][%d]/s", VDST_MAXCHANNELS );
    fDST_tree->Branch( "dead", fDSTdead, tname );
    if( !iCalibrationTree )
    {
        sprintf( tname, "zerosuppressed[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "zerosuppressed", fDSTZeroSuppressed, tname );
    }
    sprintf( tname, "sumwindow[ntel_data][%d]/s", VDST_MAXCHANNELS );
    fDST_tree->Branch( "sumwindow", fDSTsumwindow, tname );
    sprintf( tname, "sumfirst[ntel_data][%d]/s", VDST_MAXCHANNELS );
    fDST_tree->Branch( "sumfirst", fDSTsumfirst, tname );
    if( !iCalibrationTree )
    {
        sprintf( tname, "tzero[ntel_data][%d]/F", VDST_MAXCHANNELS );
        fDST_tree->Branch( "tzero", fDSTt0, tname );
        sprintf( tname, "Width[ntel_data][%d]/F", VDST_MAXCHANNELS );
        fDST_tree->Branch( "Width", fDSTTraceWidth, tname );
        // timing levels
        sprintf( tname, "pulsetiminglevel[ntel_data][%d]/F", VDST_MAXTIMINGLEVELS );
        fDST_tree->Branch( "pulsetiminglevel", fDSTpulsetiminglevels, tname );
        sprintf( tname, "pulsetiming[ntel_data][%d][%d]/F", VDST_MAXTIMINGLEVELS, VDST_MAXCHANNELS );
        fDST_tree->Branch( "pulsetiming", fDSTpulsetiming, tname );
        sprintf( tname, "Max[ntel_data][%d]/S", VDST_MAXCHANNELS );
        fDST_tree->Branch( "Max", fDSTMax, tname );
    }
    if( iFullTree )
    {
        sprintf( tname, "RawMax[ntel_data][%d]/S", VDST_MAXCHANNELS );
        fDST_tree->Branch( "RawMax", fDSTRawMax, tname );
    }
    sprintf( tname, "HiLo[ntel_data][%d]/s", VDST_MAXCHANNELS );
    fDST_tree->Branch( "HiLo", fDSTHiLo, tname );
    sprintf( tname, "N255[ntel_data][%d]/s", VDST_MAXCHANNELS );
    fDST_tree->Branch( "N255", fDSTN255, tname );
    // FADC trace
    fDST_tree->Branch( "numSamples", fDSTnumSamples, "numSamples[ntel_data]/s" );
    if( fReadWriteFADC )
    {
        sprintf( tname, "FADC[ntel_data][%d][%d]/s", VDST_MAXSUMWINDOW, VDST_MAXCHANNELS );
        fDST_tree->Branch( "Trace", fDSTtrace, tname );
    }
    
    //PhotoElectrons
    if( fMC && !iCalibrationTree )//to be changed to something like fReadPE...
    {
        sprintf( tname, "Pe[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "Pe", fDSTPe, tname );
        sprintf( tname, "ADC_HG[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "ADC_HG", fDSTPadcHG, tname );
        sprintf( tname, "ADC_LG[ntel_data][%d]/s", VDST_MAXCHANNELS );
        fDST_tree->Branch( "ADC_LG", fDSTPadcLG, tname );
    }
    
    // MC block
    if( fMC )
    {
        fDST_tree->Branch( "MCprim", &fDSTprimary, "MCprimary/s" );
        fDST_tree->Branch( "MCe0", &fDSTenergy, "MCenergy/F" );
        fDST_tree->Branch( "MCxcore", &fDSTxcore, "MCxcore/F" );
        fDST_tree->Branch( "MCycore", &fDSTycore, "MCycore/F" );
        fDST_tree->Branch( "MCze", &fDSTze, "MCze/F" );
        fDST_tree->Branch( "MCaz", &fDSTaz, "MCaz/F" );
        fDST_tree->Branch( "MCxoff", &fDSTTel_xoff, "MCxoff/F" );
        fDST_tree->Branch( "MCyoff", &fDSTTel_yoff, "MCyoff/F" );
    }
    resetDataVectors();
    
    return true;
}


void VDSTTree::resetDataVectors( unsigned int iCH, unsigned int iMaxNTel, unsigned int iMaxPrevNTel, unsigned int iMaxNChannels,
                                 unsigned int iMaxNTimingLevels, unsigned int iMaxNSamples, bool iTriggerReset, bool iIsCTADST )
{
    // reset the data vectors
    if( iMaxNTel >= VDST_MAXTELESCOPES )
    {
        iMaxNTel = VDST_MAXTELESCOPES;
    }
    if( iMaxPrevNTel >= VDST_MAXTELESCOPES || iMaxPrevNTel == 0 )
    {
        iMaxPrevNTel = VDST_MAXTELESCOPES;
    }
    if( iMaxNChannels >= VDST_MAXCHANNELS )
    {
        iMaxNChannels = VDST_MAXCHANNELS;
    }
    if( iMaxNTimingLevels >= VDST_MAXTIMINGLEVELS )
    {
        iMaxNTimingLevels = VDST_MAXTIMINGLEVELS;
    }
    if( iMaxNSamples >= VDST_MAXSUMWINDOW )
    {
        iMaxNSamples = VDST_MAXSUMWINDOW;
    }
    
    // reset trigger data
    fDSTLTrig = 0;
    fDSTNTrig = 0;
    for( unsigned int i = 0; i < iMaxNTel; i++ )
    {
        fDSTLTrig_list[i] = 0;
        fDSTnL1trig[i] = 0;
        fDSTLTtime[i] = 0.;
        fDSTL2TrigType[i] = 0;
    }
    // return if previous event changed trigger variables only
    if( iTriggerReset )
    {
        return;
    }
    
    // loop over telescopes
    for( unsigned int i = 0; i < iMaxNTel; i++ )
    {
        fDSTpointAzimuth[i] = 0.;
        fDSTpointElevation[i] = 0.;
        fDSTpointTrackingKnown[i] = 0;
    }
    
    // for all non-CTA DSTs
    if( !iIsCTADST )
    {
        for( unsigned int i = 0; i < iMaxNTel; i++ )
        {
            fDSTLDTtime[i] = 0.;
            for( unsigned int j = 0; j < iMaxNChannels; j++ )
            {
                fDSTt0[i][j] = 0.;
                fDSTRawMax[i][j] = 0;
                fDSTTraceWidth[i][j] = 0.;
                fDSTN255[i][j] = 0;
            }
        }
    }
    
    // loop over all telescopes with data in the previous event
    for( unsigned int i = 0; i < iMaxPrevNTel; i++ )
    {
        fDSTtel_data[i] = 0;
        fDSTTelescopeZeroSupression[i] = 0;
        fDSTnumSamples[i] = 0;
    }
    
    for( unsigned int i = 0; i < iMaxPrevNTel; i++ )
    {
        for( unsigned int j = 0; j < iMaxNChannels; j++ )
        {
            if( iCH != 0 )
            {
                fDSTChan[i][j] = iCH;
            }
            else
            {
                fDSTChan[i][j] = j;
            }
        }
    }
    
    memset( fDSTdead, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTdead[0][0] ) );
    memset( fDSTZeroSuppressed, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTZeroSuppressed[0][0] ) );
    memset( fDSTL1trig, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTL1trig[0][0] ) );
    // FADC mode
    if( fReadWriteFADC )
    {
        memset( fDSTtrace, 0, VDST_MAXTELESCOPES * VDST_MAXSUMWINDOW * VDST_MAXCHANNELS * sizeof( fDSTtrace[0][0][0] ) );
    }
    // QADC mode
    else
    {
        std::fill( &fDSTRecord[0][0], &fDSTRecord[0][0] + VDST_MAXTELESCOPES * VDST_MAXCHANNELS, 1 );
        std::fill( &fDSTpedestal[0][0], &fDSTpedestal[0][0] + VDST_MAXTELESCOPES * VDST_MAXCHANNELS, 0. );
        std::fill( &fDSTsums[0][0], &fDSTsums[0][0] + VDST_MAXTELESCOPES * VDST_MAXCHANNELS, 0. );
        std::fill( &fDSTsums2[0][0], &fDSTsums2[0][0] + VDST_MAXTELESCOPES * VDST_MAXCHANNELS, 0. );
        std::fill( &fDSTpulsetiming[0][0][0], &fDSTpulsetiming[0][0][0] + VDST_MAXTELESCOPES * VDST_MAXCHANNELS * VDST_MAXTIMINGLEVELS, 0. );
        std::fill( &fDSTpulsetiminglevels[0][0], &fDSTpulsetiminglevels[0][0] + VDST_MAXTELESCOPES * VDST_MAXTIMINGLEVELS, 0. );
        memset( fDSTsumwindow, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTsumwindow[0][0] ) );
        memset( fDSTsumfirst, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTsumfirst[0][0] ) );
        memset( fDSTMax, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTMax[0][0] ) );
        if( fFillPeakADC )
        {
            memset( fDSTPadcHG, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTPadcHG[0][0] ) );
            memset( fDSTPadcLG, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTPadcLG[0][0] ) );
        }
    }
    // write PEs
    if( fFillPELeaf )
    {
        memset( fDSTPe, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTPe[0][0] ) );
    }
    memset( fDSTHiLo, 0, VDST_MAXTELESCOPES * VDST_MAXCHANNELS * sizeof( fDSTHiLo[0][0] ) );
    
}

/*
     init DST tree for reading
*/
bool VDSTTree::initDSTTree( TTree* t, TTree* c )
{
    fDST_tree = t;
    fDST_conf = c;
    
    // inform about empty trees
    if( fDST_tree && fDST_tree->GetEntries() == 0 )
    {
        cout << "DST tree found: no entries" << endl;
    }
    
    fDST_vlist_of_telescopes.clear();
    unsigned int iNChannels = 0;
    int fTelID = 0;
    fDST_conf->SetBranchAddress( "NPixel", &iNChannels );
    fDST_conf->SetBranchAddress( "TelID", &fTelID );
    fDSTntel = ( unsigned int )fDST_conf->GetEntries();
    for( int i = 0; i < fDST_conf->GetEntries(); i++ )
    {
        fDST_conf->GetEntry( i );
        
        fDSTnchannel[i] = iNChannels;
        fDST_vlist_of_telescopes.push_back( fTelID );
    }
    
    if( fDST_tree->GetBranch( "gpsyear" ) )
    {
        fFullTree = true;
    }
    if( fDST_tree->GetBranch( "MCe0" ) )
    {
        fMC = true;
    }
    
    fDST_tree->SetBranchAddress( "runNumber", &fDSTrunnumber );
    fDST_tree->SetBranchAddress( "eventNumber", &fDSTeventnumber );
    fDST_tree->SetBranchAddress( "ntel", &fDSTntel );
    if( fFullTree )
    {
        fDST_tree->SetBranchAddress( "gpsyear", &fDSTgpsyear );
        fDST_tree->SetBranchAddress( "ATgpsyear", &fDSTATgpsyear );
        fDST_tree->SetBranchAddress( "gps0", &fDSTgps0 );
        fDST_tree->SetBranchAddress( "gps1", &fDSTgps1 );
        fDST_tree->SetBranchAddress( "gps2", &fDSTgps2 );
    }
    fDST_tree->SetBranchAddress( "Paz", fDSTpointAzimuth );
    fDST_tree->SetBranchAddress( "Pel", fDSTpointElevation );
    fDST_tree->SetBranchAddress( "ltrig", &fDSTLTrig );
    fDST_tree->SetBranchAddress( "ntrig", &fDSTNTrig );
    fDST_tree->SetBranchAddress( "ltrig_list", fDSTLTrig_list );
    if( fMC )
    {
        fDST_tree->SetBranchAddress( "LTtime", fDSTLTtime );
        if( fFullTree )
        {
            fDST_tree->SetBranchAddress( "LDTtime", fDSTLDTtime );
        }
        fDST_tree->SetBranchAddress( "L2TrigType", fDSTL2TrigType );
    }
    fDST_tree->SetBranchAddress( "ntel_data", &fDSTntel_data );
    fDST_tree->SetBranchAddress( "tel_data", fDSTtel_data );
    if( fDST_tree->GetBranchStatus( "tel_zero_suppression" ) )
    {
        fDST_tree->SetBranchAddress( "tel_zero_suppression", fDSTTelescopeZeroSupression );
    }
    if( fDST_tree->GetBranchStatus( "recorded" ) )
    {
        fDST_tree->SetBranchAddress( "recorded", fDSTRecord );
    }
    fDST_tree->SetBranchAddress( "nL1trig", fDSTnL1trig );
    fDST_tree->SetBranchAddress( "L1trig", fDSTL1trig );
    if( fDST_tree->GetBranchStatus( "ped" ) )
    {
        fDST_tree->SetBranchAddress( "ped", fDSTpedestal );
    }
    
    fDST_tree->SetBranchAddress( "sum", fDSTsums );
    if( fDST_tree->GetBranchStatus( "sum2" ) )
    {
        fDST_tree->SetBranchAddress( "sum2", fDSTsums2 );
    }
    else
    {
        fDST_tree->SetBranchAddress( "sum", fDSTsums );
    }
    fDST_tree->SetBranchAddress( "dead", fDSTdead );
    if( fDST_tree->GetBranchStatus( "zerosuppressed" ) )
    {
        fDST_tree->SetBranchAddress( "zerosuppressed", fDSTZeroSuppressed );
    }
    fDST_tree->SetBranchAddress( "tzero", fDSTt0 );
    fDST_tree->SetBranchAddress( "Width", fDSTTraceWidth );
    if( fDST_tree->GetBranchStatus( "Trace" ) )
    {
        fDST_tree->SetBranchAddress( "Trace", fDSTtrace );
        setFADC( true );
    }
    if( fDST_tree->GetBranchStatus( "Pe" ) )
    {
        fDST_tree->SetBranchAddress( "Pe", fDSTPe );
        setMC( true );
    }
    if( fDST_tree->GetBranchStatus( "numSamples" ) )
    {
        fDST_tree->SetBranchAddress( "numSamples", fDSTnumSamples );
    }
    fDST_tree->SetBranchAddress( "pulsetiminglevel", fDSTpulsetiminglevels );
    fDST_tree->SetBranchAddress( "pulsetiming", fDSTpulsetiming );
    fDST_tree->SetBranchAddress( "Max", fDSTMax );
    if( fDST_tree->GetBranchStatus( "ADC_HG" )
            && fDST_tree->GetBranchStatus( "ADC_LG" ) )
    {
        fDST_tree->SetBranchAddress( "ADC_HG", fDSTPadcHG );
        fDST_tree->SetBranchAddress( "ADC_LG", fDSTPadcLG );
        //                fDST_ADC_set = true;
    }
    if( fFullTree )
    {
        fDST_tree->SetBranchAddress( "RawMax", fDSTRawMax );
    }
    fDST_tree->SetBranchAddress( "HiLo", fDSTHiLo );
    if( fMC )
    {
        fDST_tree->SetBranchAddress( "MCprim", &fDSTprimary );
        fDST_tree->SetBranchAddress( "MCe0", &fDSTenergy );
        fDST_tree->SetBranchAddress( "MCxcore", &fDSTxcore );
        fDST_tree->SetBranchAddress( "MCycore", &fDSTycore );
        fDST_tree->SetBranchAddress( "MCze", &fDSTze );
        fDST_tree->SetBranchAddress( "MCaz", &fDSTaz );
        fDST_tree->SetBranchAddress( "MCxoff", &fDSTTel_xoff );
        fDST_tree->SetBranchAddress( "MCyoff", &fDSTTel_yoff );
    }
    
    // reset data vectors
    resetDataVectors();
    
    return fMC;
}


int VDSTTree::hasData( int iTelID )
{
    if( iTelID < 0 )
    {
        cout << "\tB" << endl;
        return -1;
    }
    if( iTelID >= ( int )fDST_vlist_of_telescopes.size() )
    {
        return -1;
    }
    
    for( unsigned int j = 0; j < fDSTntel_data; j++ )
    {
        if( fDST_vlist_of_telescopes[iTelID] == fDSTtel_data[j] )
        {
            return j;
        }
    }
    return -1;
}

int VDSTTree::getDSTTelescopeNumber( unsigned int iTelHyperArray_ID )
{
    for( unsigned int i = 0; i < fDST_vlist_of_telescopes.size(); i++ )
    {
        if( fDST_vlist_of_telescopes[i] == iTelHyperArray_ID )
        {
            return ( int )i;
        }
    }
    return -1;
}


int VDSTTree::hasLocalTrigger( int iTelID )
{
    if( iTelID < 0 )
    {
        return -1;
    }
    if( iTelID >= ( int )fDST_vlist_of_telescopes.size() )
    {
        return -1;
    }
    
    for( unsigned int j = 0; j < fDSTNTrig; j++ )
    {
        if( fDST_vlist_of_telescopes[iTelID] == fDSTLTrig_list[j] )
        {
            return j;
        }
    }
    return -1;
}


double VDSTTree::getDSTPedestal( int iChannelID, bool iPrint )
{
    if( iPrint )
    {
        cout << "\t getDSTPedestal " << fTelescopeCounter_temp << "\t channel " << iChannelID;
        cout << "\t max channel " << VDST_MAXCHANNELS << endl;
        cout << "\t pedestal " << fDSTpedestal[fTelescopeCounter_temp][iChannelID];
        cout << endl;
    }
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else
    {
        return fDSTpedestal[fTelescopeCounter_temp][iChannelID];
    }
    return 0;
}

/*
 * return sum from DST tree
 *
 * in case the adc values are stored directly: return adc value
 * (might be better to add a flag for this)
 *
 */

double VDSTTree::getDSTSums( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else if( fDST_ADC_set )
    {
        if( getDSTHiLo( iChannelID ) )
        {
            return ( double )fDSTPadcHG[fTelescopeCounter_temp][iChannelID];
        }
        else
        {
            return ( double )fDSTPadcLG[fTelescopeCounter_temp][iChannelID];
        }
    }
    else
    {
        return fDSTsums[fTelescopeCounter_temp][iChannelID];
    }
    return 0;
}

unsigned short int  VDSTTree::getDSTPe( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else
    {
        return fDSTPe[fTelescopeCounter_temp][iChannelID];
    }
    return 0;
}

unsigned short int VDSTTree::getDSTPadcHG( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getDSTPadcHG( iChannelID );
}

unsigned short int VDSTTree::getDSTPadcHG( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0;
    }
    return ( unsigned short int )( fDSTPadcHG[fTelescopeCounter_temp][iChannelID] );
}
unsigned short int VDSTTree::getDSTPadcLG( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getDSTPadcLG( iChannelID );
}

unsigned short int VDSTTree::getDSTPadcLG( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0;
    }
    
    return ( unsigned short int )( fDSTPadcLG[fTelescopeCounter_temp][iChannelID] );
}

double VDSTTree::getDSTMax( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getDSTMax( iChannelID );
}

double VDSTTree::getDSTMax( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else
    {
        return ( double )( fDSTMax[fTelescopeCounter_temp][iChannelID] );
    }
    
    return 0;
}

double VDSTTree::getDSTRawMax( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getDSTRawMax( iChannelID );
}

double VDSTTree::getDSTRawMax( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else
    {
        return ( double )( fDSTRawMax[fTelescopeCounter_temp][iChannelID] ) / 100.;
    }
    
    return 0;
}


double VDSTTree::getDSTWidth( int iTelID, int iChannelID )
{
    int iTel = hasData( iTelID );
    
    if( iTel < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else
    {
        return ( double )( fDSTTraceWidth[iTel][iChannelID] );
    }
    
    return 0;
}


double VDSTTree::getDSTTZeros( int iTelID, int iChannelID )
{
    int iTel = hasData( iTelID );
    
    if( iTel < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else
    {
        return ( double )( fDSTt0[iTel][iChannelID] );
    }
    
    return 0;
}

double VDSTTree::getDSTpulsetiming( int iTelID, int iChannelID, int iTimingLevelN )
{
    setTelCounter( iTelID );
    return getDSTpulsetiming( iChannelID, iTimingLevelN );
}

double VDSTTree::getDSTpulsetiming( int iChannelID, int iTimingLevelN )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0.;
    }
    else if( iTimingLevelN < VDST_MAXTIMINGLEVELS )
    {
        return ( double )( fDSTpulsetiming[fTelescopeCounter_temp][iTimingLevelN][iChannelID] );
    }
    
    return 0;
}

unsigned int VDSTTree::getDSTDead( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getDSTDead( iChannelID );
}

unsigned int VDSTTree::getDSTDead( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0;
    }
    else
    {
        return fDSTdead[fTelescopeCounter_temp][iChannelID];
    }
    
    return 0;
}

unsigned short int VDSTTree::getZeroSupppressed( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getZeroSupppressed( iChannelID );
}

unsigned short int VDSTTree::getZeroSupppressed( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0;
    }
    else
    {
        return fDSTZeroSuppressed[fTelescopeCounter_temp][iChannelID];
    }
    
    // default return value: not suppressed
    return 0;
}

unsigned short int VDSTTree::getDSTNumSample( unsigned int iTelID )
{
    int iTel = hasData( iTelID );
    
    if( iTel >= 0 )
    {
        return fDSTnumSamples[iTelID];
    }
    
    return 0;
}

UShort_t VDSTTree::getDSTHiLo( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getDSTHiLo( iChannelID );
}

UShort_t VDSTTree::getDSTHiLo( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0;
    }
    else
    {
        return fDSTHiLo[fTelescopeCounter_temp][iChannelID];
    }
    
    return 0;
}

unsigned short int VDSTTree::getDSTTrace( unsigned int iTelID, unsigned int iChannelID, unsigned short int iSample )
{
    setTelCounter( iTelID );
    return getDSTTrace( iChannelID, iSample );
}

unsigned short int VDSTTree::getDSTTrace( unsigned int iChannelID, unsigned short int iSample )
{
    if( fTelescopeCounter_temp < 0 )
    {
        return 3;
    }
    
    return fDSTtrace[fTelescopeCounter_temp][iSample][iChannelID];
}

unsigned int VDSTTree::getTrigL1( int iTelID, int iChannelID )
{
    setTelCounter( iTelID );
    return getTrigL1( iChannelID );
}

unsigned int VDSTTree::getTrigL1( int iChannelID )
{
    if( fTelescopeCounter_temp < 0 && iChannelID < VDST_MAXCHANNELS )
    {
        return 0;
    }
    else
    {
        return fDSTL1trig[fTelescopeCounter_temp][iChannelID];
    }
    
    return 0;
}


bool VDSTTree::getDSTLocalTrigger( int iTelID )
{
    int iTel = hasLocalTrigger( iTelID );
    if( iTel < 0 )
    {
        return false;
    }
    
    return true;
}


float VDSTTree::getDSTLocalTriggerTime( int iTelID )
{
    int iTel = hasLocalTrigger( iTelID );
    if( iTel < 0 )
    {
        return 0;
    }
    
    return fDSTLTtime[iTel];
}

unsigned short int VDSTTree::getDSTL2TriggerType( int iTelID )
{
    int iTel = hasLocalTrigger( iTelID );
    if( iTel < 0 )
    {
        return 0;
    }
    
    return fDSTL2TrigType[iTel];
}

float VDSTTree::getDSTLocalDelayedTriggerTime( int iTelID )
{
    int iTel = hasLocalTrigger( iTelID );
    if( iTel < 0 )
    {
        return 0;
    }
    
    return fDSTLDTtime[iTel];
}

/*

     reading array configuration from an ASCII file
     (simple text file)

     expect maximum five colums (only four are read here)

     <TelID>  <FOV(in [deg])>  <Dynamic range (bit)>  <Raw/calibrated sum> <CTAO telescope ID>

    Note: Raw/calibrated sum is important for analyses using no trace integration

*/
map< unsigned int, VDSTTelescopeConfiguration> VDSTTree::readArrayConfig( string iFile )
{
    fDST_list_of_telescopes.clear();
    fDST_vlist_of_telescopes.clear();
    if( iFile.size() == 0 )
    {
        return fDST_list_of_telescopes;
    }
    
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VDSTTree::readArrayConfig error: file not found, " << iFile << endl;
        return fDST_list_of_telescopes;
    }
    string iLine;
    string iT1;
    string iT2;
    cout << "reading sub-array configuration from " << iFile << endl;
    
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            // lines with comments
            if( iLine.substr( 0, 1 ) == "#" )
            {
                continue;
            }
            istringstream is_stream( iLine );
            is_stream >> iT1;
            // FOV
            fDST_list_of_telescopes[atoi( iT1.c_str() )].FOV = 20.5;
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> fDST_list_of_telescopes[atoi( iT1.c_str() )].FOV;
            }
            // Dynamic range (in cfg file in BIT)
            fDST_list_of_telescopes[atoi( iT1.c_str() )].DynamicRange = -1.;
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iT2;
                fDST_list_of_telescopes[atoi( iT1.c_str() )].DynamicRange = TMath::Power( 2., atof( iT2.c_str() ) );
            }
            // use RAW sum or calibrated sum (default)
            fDST_list_of_telescopes[atoi( iT1.c_str() )].RAWsum = false;
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iT2;
                fDST_list_of_telescopes[atoi( iT1.c_str() )].RAWsum = ( bool )atoi( iT2.c_str() );
            }
            // read telescope name
            fDST_list_of_telescopes[atoi( iT1.c_str() )].TelescopeName = "";
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iT2;
                fDST_list_of_telescopes[atoi( iT1.c_str() )].TelescopeName = iT2;
            }
            fDST_vlist_of_telescopes.push_back( atoi( iT1.c_str() ) );
        }
    }
    is.close();
    
    map< unsigned int, VDSTTelescopeConfiguration >::iterator iter;
    unsigned int z = 0;
    for( iter = fDST_list_of_telescopes.begin(); iter != fDST_list_of_telescopes.end(); ++iter )
    {
        cout << "\t Telescope ID " << iter->first << "  FOV " << iter->second.FOV;
        cout << "  Dyn " << iter->second.DynamicRange << " RAW " << iter->second.RAWsum;
        if(  iter->second.TelescopeName.size() > 0 )
        {
            cout << " Name " << iter->second.TelescopeName;
        }
        cout << " (#" << z << ")" << endl;
        z++;
    }
    
    return fDST_list_of_telescopes;
}


float VDSTTree::getDSTMCaz()
{
    return fDSTaz;
}

unsigned int VDSTTree::getDSTpulsetiminglevelsN()
{
    return ( unsigned int )( VDST_MAXTIMINGLEVELS );
}

int VDSTTree::setTelCounter( int iTelID )
{
    fTelescopeCounter_temp = hasData( iTelID );
    
    return fTelescopeCounter_temp;
}

void VDSTTree::fillDSTMeanPulseTiming( unsigned int iTelID, unsigned int iChannelID, double iTime, int iNSamples )
{
    int iTel = getDSTTelescopeNumber( iTelID + 1 );
    if( iTel < 0 || iNSamples == 0 )
    {
        return;
    }
    
    if( iTel < VDST_MAXTELESCOPES && iChannelID < VDST_MAXCHANNELS )
    {
        fDSTMeanPulseTiming[iTel][iChannelID] += iTime;
        fDSTMeanPulseTiming_N[iTel][iChannelID]++;
        if( !fDSTMeanPulseTimingHistogram[iTel] )
        {
            char hname[200];
            sprintf( hname, "hPT_%d", iTel );
            fDSTMeanPulseTimingHistogram[iTel] = new TH1F( hname, "", 100 * iNSamples, 0., ( float )iNSamples );
        }
        fDSTMeanPulseTimingHistogram[iTel]->Fill( iTime );
    }
}

double VDSTTree::getDSTMeanPulseTiming( unsigned int iTelID, unsigned int iChannelID )
{
    int iTel = getDSTTelescopeNumber( iTelID + 1 );
    if( iTel < 0. )
    {
        return -9999.;
    }
    
    if( iTel < VDST_MAXTELESCOPES && iChannelID < VDST_MAXCHANNELS )
    {
        if( fDSTMeanPulseTiming_N[iTel][iChannelID] > 0. )
        {
            return fDSTMeanPulseTiming[iTel][iChannelID] / fDSTMeanPulseTiming_N[iTel][iChannelID];
        }
        else
        {
            return -9999.;
        }
    }
    
    return 0.;
}


double VDSTTree::getDSTMedianPulseTimingPerTelescope( unsigned int iTelID )
{
    int iTel = getDSTTelescopeNumber( iTelID + 1 );
    if( iTel < 0. )
    {
        return -9999.;
    }
    
    if( iTel < VDST_MAXTELESCOPES )
    {
        if( fDSTMeanPulseTimingHistogram[iTel] )
        {
            double i_a[] = { 0.5 };
            double i_b[] = { 0.0 };
            fDSTMeanPulseTimingHistogram[iTel]->GetQuantiles( 1, i_b, i_a );
            return i_b[0];
        }
        else
        {
            return -9999.;
        }
    }
    
    return 0.;
}


double VDSTTree::getDSTMeanPulseTimingPerTelescope( unsigned int iTelID )
{
    int iTel = getDSTTelescopeNumber( iTelID + 1 );
    if( iTel < 0. )
    {
        return -9999.;
    }
    
    if( iTel < VDST_MAXTELESCOPES )
    {
        if( fDSTMeanPulseTimingHistogram[iTel] )
        {
            return fDSTMeanPulseTimingHistogram[iTel]->GetMean();
        }
        else
        {
            return -9999.;
        }
    }
    return -9999.;
}

double VDSTTree::getDSTRMSPulseTimingPerTelescope( unsigned int iTelID )
{
    int iTel = getDSTTelescopeNumber( iTelID + 1 );
    if( iTel < 0. )
    {
        return -9999.;
    }
    
    if( iTel < VDST_MAXTELESCOPES )
    {
        if( fDSTMeanPulseTimingHistogram[iTel] )
        {
            return fDSTMeanPulseTimingHistogram[iTel]->GetRMS();
        }
        else
        {
            return -9999.;
        }
    }
    return -9999.;
}

double VDSTTree::getDSTNEventsPulseTimingPerTelescope( unsigned int iTelID )
{
    int iTel = getDSTTelescopeNumber( iTelID + 1 );
    if( iTel < 0. )
    {
        return -9999.;
    }
    
    if( iTel < VDST_MAXTELESCOPES )
    {
        if( fDSTMeanPulseTimingHistogram[iTel] )
        {
            return fDSTMeanPulseTimingHistogram[iTel]->GetEntries();
        }
        else
        {
            return -9999.;
        }
    }
    return -9999.;
}
