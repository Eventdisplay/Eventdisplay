/*! \class VCalibrationData
     \brief all calibration data is stored here

*/

#include <VCalibrationData.h>

VCalibrationData::VCalibrationData( unsigned int iTel, string iPedfile, string iGainfile, string iTofffile,
                                    string iPedLowGainfile, string iGainLowGainfile, string iToffLowGainfile,
                                    string iLowGainMultfile, string iTzerofile, string iTzeroLowGainfile,
                                    string iObservatory )
{
    fUsePedestalsInTimeSlices = true;
    fLowGainUsePedestalsInTimeSlices = true;
    
    fTelID = iTel;
    
    fFileName.push_back( iPedfile );
    fHistoName.push_back( "pedestal" );
    fFileName.push_back( iGainfile );
    fHistoName.push_back( "gain" );
    fFileName.push_back( iTofffile );
    fHistoName.push_back( "toff" );
    fFileName.push_back( iPedLowGainfile );
    fHistoName.push_back( "pedestal_lowGain" );
    fFileName.push_back( iGainLowGainfile );
    fHistoName.push_back( "gain_lowGain" );
    fFileName.push_back( iToffLowGainfile );
    fHistoName.push_back( "toff_lowGain" );
    fFileName.push_back( iLowGainMultfile );
    fHistoName.push_back( "lowgain_mult" );
    fFileName.push_back( iTzerofile );
    fHistoName.push_back( "tzero" );
    fFileName.push_back( iTzeroLowGainfile );
    fHistoName.push_back( "tzero_lowGain" );
    
    for( unsigned int i = 0; i < fFileName.size(); i++ )
    {
        fFile.push_back( 0 );
    }
    
    
    fPedFromPLine = false;
    
    fReader = 0;
    
    fBoolLowGainPedestals = false;
    fBoolLowGainGains = false;
    fBoolLowGainTOff = false;
    
    fSumWindow = 0;
    
    fTS_ped_temp_time = 0.;
    fLowGainMultiplier_Trace = 0;
    
    setAverageTZero( 0., true );
    setAverageTZero( 0., false );
    
    // summary histograms
    
    hisList = new TList();
    
    char hname[200];
    char hname_var[200];
    for( unsigned int i = 0; i < fFile.size(); i++ )
    {
        sprintf( hname, "h%s_%d", fHistoName[i].c_str(), iTel + 1 );
        sprintf( hname_var, "h%s_var_%d", fHistoName[i].c_str(), iTel + 1 );
        if( i == C_PED || i == C_PEDLOW )
        {
            if( iObservatory == "CTA" )
            {
                fHisto_mean.push_back( new TH1F( hname, "", 800, 0., 200. ) );
            }
            else
            {
                fHisto_mean.push_back( new TH1F( hname, "", 200, 0., 50. ) );
            }
            fHisto_mean.back()->SetXTitle( "mean pedestal [dc]" );
            if( iObservatory == "CTA" )
            {
                fHisto_variance.push_back( new TH1F( hname_var, "", 800, 0., 200. ) );
            }
            else
            {
                fHisto_variance.push_back( new TH1F( hname_var, "", 200, 0., 50. ) );
            }
            fHisto_variance.back()->SetXTitle( "pedestal variance [dc]" );
        }
        else if( i == C_TOFF || i == C_TOFFLOW )
        {
            fHisto_mean.push_back( new TH1F( hname, "", 200, -20., 20. ) );
            fHisto_mean.back()->SetXTitle( "time offset [sample]" );
            fHisto_variance.push_back( new TH1F( hname_var, "", 200, 0., 10. ) );
            fHisto_variance.back()->SetXTitle( "time offset variance [sample]" );
        }
        else if( i == C_GAIN || i == C_GAINLOW )
        {
            fHisto_mean.push_back( new TH1F( hname, "", 200, 0., 3. ) );
            fHisto_mean.back()->SetXTitle( "relative gain" );
            fHisto_variance.push_back( new TH1F( hname_var, "", 200, 0., 2. ) );
            fHisto_variance.back()->SetXTitle( "relative gain variance" );
        }
        else if( i == C_LOWGAIN )
        {
            fHisto_mean.push_back( new TH1F( hname, "", 200, 0., 20. ) );
            fHisto_mean.back()->SetXTitle( "low gain multiplier" );
            fHisto_variance.push_back( new TH1F( hname_var, "", 200, 0., 20. ) );
            fHisto_variance.back()->SetXTitle( "low gain multiplier variance" );
        }
        else if( i == C_TZERO || i == C_TZEROLOW )
        {
            fHisto_mean.push_back( new TH1F( hname, "", 1000, 0., 20. ) );
            fHisto_mean.back()->SetXTitle( "average pulse time t_{0} [sample]" );
            fHisto_variance.push_back( new TH1F( hname_var, "", 200, 0., 20. ) );
            fHisto_variance.back()->SetXTitle( "average pulse time t_{0} (variance) [sample]" );
        }
        
        hisList->Add( fHisto_mean.back() );
        hisList->Add( fHisto_variance.back() );
    }
    
    
}


void VCalibrationData::initialize( unsigned int i_channel, unsigned int nSamples, bool iUsePedestalsInTimeSlices,
                                   bool iLowGainUsePedestalsInTimeSlices, bool iPedsFromPLine, bool iReadCalibDB,
                                   bool i_isDSTMC, bool iDebug, int iRunMode, bool isTeltoAna )
{
    if( iDebug )
    {
        cout << "VCalibrationData::initialize " << i_channel << "\t" << fTelID << endl;
    }
    
    fPedFromPLine = iPedsFromPLine;
    fUsePedestalsInTimeSlices = iUsePedestalsInTimeSlices;
    fLowGainUsePedestalsInTimeSlices = iLowGainUsePedestalsInTimeSlices;
    
    for( unsigned int i = 0; i < fHisto_mean.size(); i++ )
    {
        if( fHisto_mean[i] )
        {
            fHisto_mean[i]->Reset();
        }
    }
    for( unsigned int i = 0; i < fHisto_variance.size(); i++ )
    {
        if( fHisto_variance[i] )
        {
            fHisto_variance[i]->Reset();
        }
    }
    
    char c_name[2000];
    for( unsigned int i = 0; i < fFileName.size(); i++ )
    {
        if( fFileName[i].size() > 0 && !fPedFromPLine && !i_isDSTMC )
        {
            // readcalibDB: gains and toffs are read from the VDB
            if( fFileName[i].find( "gain" ) != string::npos && iReadCalibDB )
            {
                fFile[i] = 0;
            }
            else if( fFileName[i].find( "toff" ) != string::npos && iReadCalibDB )
            {
                fFile[i] = 0;
            }
            // tzero calculation (file is opened later)
            else if( fFileName[i].find( "tzero" ) != string::npos && iRunMode == 7 )
            {
                fFile[i] = 0;
            }
            else
            {
                if( isTeltoAna && iRunMode != 2 )
                {
                    sprintf( c_name, "%s.root", fFileName[i].c_str() );
                    fFile[i] = new TFile( c_name, "READ" );
                    if( fFile[i]->IsZombie() )
                    {
                        fFile[i] = 0;
                    }
                }
                else
                {
                    fFile[i] = 0;
                }
            }
        }
        else
        {
            fFile[i] = 0;
        }
    }
    
    fFADCStopOffsets.resize( i_channel, 0. );
    fChannelStatus.resize( i_channel, 1 );
    
    valarray< double > itemp_ped;
    valarray< valarray< double > > itemp_pedvars;
    itemp_pedvars.resize( nSamples + 1, itemp_ped );
    
    // (time dependent pedestal vectors are initialzed in VCalibrator)
    
    // high gain channels
    fPeds.resize( i_channel, 0. );
    fPeds_perEvent.resize( i_channel, 0. );
    fVPedvars.resize( nSamples + 1, sqrt( fPeds ) );
    fPedrms.resize( i_channel, 0. );
    fVmeanPedvars.resize( nSamples + 1, 0. );
    fVmeanRMSPedvars.resize( nSamples + 1, 0. );
    
    fLowGainPedsrms.resize( i_channel, 0. );
    
    fGains.resize( i_channel, 0. );
    fGains_DefaultSetting.resize( i_channel, true );
    fGainvars.resize( i_channel, 0. );
    
    fTOffsets.resize( i_channel, 0. );
    fTOffsetvars.resize( i_channel, 0. );
    
    fAverageTzero.resize( i_channel, 0. );
    fAverageTzerovars.resize( i_channel, 0. );
    
    fFADCtoPhe.resize( i_channel, 1. );
    
    // low gain channels
    fLowGainPeds.resize( i_channel, 0. );
    fVLowGainPedvars.resize( nSamples + 1, fLowGainPeds );
    fVmeanLowGainPedvars.resize( nSamples + 1, 0. );
    fVmeanRMSLowGainPedvars.resize( nSamples + 1, 0. );
    
    fLowGainGains.resize( i_channel, 0. );
    fLowGainGains_DefaultSetting.resize( i_channel, 0. );
    fLowGainGainvars.resize( i_channel, 0. );
    
    fLowGainTOffsets.resize( i_channel, 0. );
    fLowGainTOffsetvars.resize( i_channel, 0. );
    fLowGainAverageTzero.resize( i_channel, 0. );
    fLowGainAverageTzerovars.resize( i_channel, 0. );
    
    fLowGainFADCtoPhe.resize( i_channel, 1. );
    
    itemp_ped.resize( i_channel, 1. );
    // low gain multiplier settings
    fLowGainMultiplier_Trace = 0;
    fLowGainMultiplier_Camera.resize( i_channel, 0. );
    
    if( iDebug )
    {
        cout << "END VCalibrationData::initialize " << endl;
    }
}

TH1F* VCalibrationData::getHistogram( unsigned int iTel, unsigned int iChannel, unsigned int iWindowSize, VCalibrationData::E_PEDTYPE iType, ULong64_t iTelType )
{
    char iHName[200];
    std::ostringstream iSname;
    
    if( iType == C_PED || iType == C_PEDLOW )
    {
        sprintf( iHName, "distributions/hped_%d_%d_%d", iTel, iWindowSize, iChannel );
        if( iType < ( int )fFile.size() && fFile[iType] )
        {
            if( !fFile[iType]->Get( iHName ) )
            {
                sprintf( iHName, "distributions_%d/hped_%d_%d_%d", iTel + 1, iTel, iWindowSize, iChannel );
            }
            if( !fFile[iType]->Get( iHName ) && iTelType != 99999 )
            {
                iSname << "distributions_" << iTelType << "/hped_" << iTelType << "_" << iWindowSize << "_" << iChannel;
                sprintf( iHName, "%s", iSname.str().c_str() );
            }
        }
    }
    else if( iType == C_TOFF || iType == C_TOFFLOW )
    {
        sprintf( iHName, "htoff_%d", iChannel );
    }
    else if( iType == C_TZERO || iType == C_TZEROLOW )
    {
        sprintf( iHName, "htzero_%d_%d", iTel + 1, iChannel );
    }
    else if( iType == C_GAIN || iType == C_GAINLOW )
    {
        sprintf( iHName, "hgain_%d", iChannel );
    }
    if( iType < ( int )fFile.size() && fFile[iType] )
    {
        return ( TH1F* )fFile[iType]->Get( iHName );
    }
    
    return 0;
}


TH1F* VCalibrationData::getHistoPed( unsigned int iTel, unsigned int iChannel, unsigned int iWindowSize, bool iLowGain, ULong64_t iTelType )
{
    if( fReader && fPedFromPLine )
    {
        if( fReader->getDataFormat() == "grisu" )
        {
            return fReader->getPedHisto( iTel, iChannel );
        }
        else
        {
            return 0;
        }
    }
    if( iLowGain )
    {
        return getHistogram( iTel, iChannel, iWindowSize, C_PEDLOW, iTelType );
    }
    
    return getHistogram( iTel, iChannel, iWindowSize, C_PED, iTelType );
}


TH1F* VCalibrationData::getHistoGain( unsigned int iTel, unsigned int iChannel, bool iLowGain )
{
    if( iLowGain )
    {
        return getHistogram( iTel, iChannel, 0, C_GAINLOW );
    }
    
    return getHistogram( iTel, iChannel, 0, C_GAIN );
}


TH1F* VCalibrationData::getHistoToff( unsigned int iTel, unsigned int iChannel, bool iLowGain )
{
    if( iLowGain )
    {
        return getHistogram( iTel, iChannel, 0, C_TOFFLOW );
    }
    
    return getHistogram( iTel, iChannel, 0, C_TOFF );
}

TH1F* VCalibrationData::getHistoAverageTzero( unsigned int iTel, unsigned int iChannel, bool iLowGain )
{
    if( iLowGain )
    {
        return getHistogram( iTel, iChannel, 0, C_TZEROLOW );
    }
    
    return getHistogram( iTel, iChannel, 0, C_TZERO );
}

TH1F*  VCalibrationData::getHistoDist( int iType, bool iDist )
{
    if( !iDist )
    {
        if( iType < ( int )fHisto_variance.size() )
        {
            return fHisto_variance[iType];
        }
    }
    else
    {
        if( iType < ( int )fHisto_mean.size() )
        {
            return fHisto_mean[iType];
        }
    }
    return 0;
}

double VCalibrationData::getAverageTZero( bool iLowGain )
{
    if( iLowGain )
    {
        return fAverageTZero_lowgain;
    }
    
    return fAverageTZero_highgain;
}

void  VCalibrationData::setAverageTZero( double iAverageTzero, bool iLowGain )
{
    if( iLowGain )
    {
        fAverageTZero_lowgain = iAverageTzero;
    }
    
    fAverageTZero_highgain = iAverageTzero;
}


bool VCalibrationData::terminate( vector< unsigned int > iDead, vector< unsigned int > iDeadLow, unsigned int iTraceIntegrationMethod, bool iDST )
{
    TDirectory* iDir = gDirectory;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // no DST analysis
    if( !iDST )
    {
    
        // fill tree with used calibration data
        const unsigned int iMAXSUMWINDOWS = 1000;
        const unsigned int iMAXDEFWINDOWS = 10;
        char hname[200];
        char htitle[200];
        sprintf( hname, "calib_%d", fTelID + 1 );
        sprintf( htitle, "calibration data (Telescope %d)", fTelID + 1 );
        int ipix = 0;
        double iped = 0.;
        double ipedvar = 0.;
        double igain = 0.;
        double igainvar = 0.;
        double itoff = 0.;
        double itoffvar = 0.;
        double itzero = 0.;
        double itzerovar = 0.;
        double ipedlowgain = 0.;
        double ipedvarlowgain = 0.;
        double iFADCtoPhe = 0.;
        double iFADCtoPhelowgain = 0.;
        double igainlowgain = 0.;
        double ilowgainmultiplier_trace = 0;
        double ilowgainmultiplier_sum[iMAXDEFWINDOWS][iMAXSUMWINDOWS];
        unsigned int inlowgaindefaultsumwindows = 0;
        double ilowgaindefaultsumwindows[iMAXDEFWINDOWS];
        double igainvarlowgain = 0.;
        double itofflowgain = 0.;
        double itoffvarlowgain = 0.;
        double itzerolowgain = 0.;
        double itzerovarlowgain = 0.;
        unsigned int nsumwindows = 0;
        unsigned int sumwindow = 0;
        unsigned int iintegrationMethod = iTraceIntegrationMethod;
        int istat = 0;
        int istatLow = 0;
        double pedvarV[iMAXSUMWINDOWS];
        double pedvarLowGainV[iMAXSUMWINDOWS];
        for( unsigned int j = 0; j < iMAXDEFWINDOWS; j++ )
        {
            ilowgaindefaultsumwindows[j] = 0.;
        }
        for( unsigned int i = 0; i < iMAXSUMWINDOWS; i++ )
        {
            pedvarV[i] = 0.;
            pedvarLowGainV[i] = 0.;
            for( unsigned int j = 0; j < iMAXDEFWINDOWS; j++ )
            {
                ilowgainmultiplier_sum[j][i] = 0;
            }
        }
        TTree* fCalibrationTree = new TTree( hname, htitle );
        fCalibrationTree->Branch( "pix", &ipix, "pix/I" );
        fCalibrationTree->Branch( "ped", &iped, "ped/D" );
        fCalibrationTree->Branch( "state", &istat, "state/I" );
        fCalibrationTree->Branch( "stateLow", &istatLow, "stateLow/I" );
        fCalibrationTree->Branch( "sumwindow", &sumwindow, "sumwindow/i" );
        fCalibrationTree->Branch( "pedvar", &ipedvar, "pedvar/D" );
        fCalibrationTree->Branch( "nsumwindows", &nsumwindows, "nsumwindows/i" );
        fCalibrationTree->Branch( "pedvarV", pedvarV, "pedvarV[nsumwindows]/D" );
        fCalibrationTree->Branch( "gain", &igain, "gain/D" );
        fCalibrationTree->Branch( "gainvar", &igainvar, "gainvar/D" );
        fCalibrationTree->Branch( "FADCtoPhe", &iFADCtoPhe, "FADCtoPhe/D" );
        fCalibrationTree->Branch( "toff", &itoff, "toff/D" );
        fCalibrationTree->Branch( "toffvar", &itoffvar, "toffvar/D" );
        fCalibrationTree->Branch( "tzero", &itzero, "tzero/D" );
        fCalibrationTree->Branch( "tzerovar", &itzerovar, "tzerovar/D" );
        fCalibrationTree->Branch( "pedLowGain", &ipedlowgain, "pedLowGain/D" );
        fCalibrationTree->Branch( "pedvarLowGain", &ipedvarlowgain, "pedvarLowGain/D" );
        fCalibrationTree->Branch( "pedvarVLowGain", pedvarLowGainV, "pedvarVLowGain[nsumwindows]/D" );
        fCalibrationTree->Branch( "gainLowGain", &igainlowgain, "gainLowGain/D" );
        fCalibrationTree->Branch( "gainvarLowGain", &igainvarlowgain, "gainvarLowGain/D" );
        fCalibrationTree->Branch( "FADCtoPheLowGain", &iFADCtoPhelowgain, "FADCtoPheLowGain/D" );
        fCalibrationTree->Branch( "toffLowGain", &itofflowgain, "toffLowGain/D" );
        fCalibrationTree->Branch( "toffvarLowGain", &itoffvarlowgain, "toffvarLowGain/D" );
        fCalibrationTree->Branch( "tzeroLowGain", &itzerolowgain, "tzeroLowGain/D" );
        fCalibrationTree->Branch( "tzerovarLowGain", &itzerovarlowgain, "tzerovarLowGain/D" );
        fCalibrationTree->Branch( "lowgainmultiplier_trace", &ilowgainmultiplier_trace, "lowgainmultiplier_trace/D" );
        fCalibrationTree->Branch( "integrationMethod", &iintegrationMethod, "integrationMethod/I" );
        fCalibrationTree->Branch( "lowgainmultiplier_sum", &ilowgainmultiplier_sum, "lowgainmultiplier_sum[10][1000]/D" );
        fCalibrationTree->Branch( "nlowgaindefaultsumwindows", &inlowgaindefaultsumwindows, "nlowgaindefaultsumwindows/i" );
        fCalibrationTree->Branch( "lowgaindefaultsumwindows", &ilowgaindefaultsumwindows, "lowgaindefaultsumwindows[nlowgaindefaultsumwindows]/D" );
        
        
        if( fPeds.size() == fPedrms.size() && fPeds.size() == fGains.size() && fPeds.size() == fGainvars.size()
                && fPeds.size() == fTOffsets.size() && fPeds.size() == fTOffsetvars.size() && fPeds.size() == fLowGainPeds.size()
                && fPeds.size() == fFADCtoPhe.size() && fPeds.size() == fLowGainFADCtoPhe.size() )
        {
            for( unsigned int i = 0; i < fPeds.size(); i++ )
            {
                ipix = ( int )i;
                if( i < fPeds.size() )
                {
                    iped = fPeds[i];
                }
                if( i < iDead.size() )
                {
                    istat = iDead[i];
                }
                if( i < iDeadLow.size() )
                {
                    istatLow = iDeadLow[i];
                }
                sumwindow = fSumWindow;
                if( sumwindow < fVPedvars.size() )
                {
                    ipedvar = fVPedvars[sumwindow][i];
                }
                else
                {
                    ipedvar = 0.;
                }
                
                nsumwindows = fVPedvars.size();
                if( nsumwindows < iMAXSUMWINDOWS )
                {
                    for( unsigned int s = 0; s < fVPedvars.size(); s++ )
                    {
                        pedvarV[s] = fVPedvars[s][i];
                    }
                }
                if( i < fGains.size() )
                {
                    igain = fGains[i];
                }
                if( i < fGains_DefaultSetting.size() )
                {
                    fGains_DefaultSetting[i] = false;
                }
                if( i < fGainvars.size() )
                {
                    igainvar = fGainvars[i];
                }
                if( i < fFADCtoPhe.size() )
                {
                    iFADCtoPhe = fFADCtoPhe[i];
                }
                if( i < fTOffsets.size() )
                {
                    itoff = fTOffsets[i];
                }
                if( i < fTOffsetvars.size() )
                {
                    itoffvar = fTOffsetvars[i];
                }
                if( i < fAverageTzero.size() )
                {
                    itzero = fAverageTzero[i];
                }
                if( i < fAverageTzerovars.size() )
                {
                    itzerovar = fAverageTzerovars[i];
                }
                if( i < fLowGainPeds.size() )
                {
                    ipedlowgain = fLowGainPeds[i];
                }
                if( sumwindow < fVLowGainPedvars.size() && i < fVLowGainPedvars[sumwindow].size() )
                {
                    ipedvarlowgain = fVLowGainPedvars[sumwindow][i];
                }
                else
                {
                    ipedvarlowgain = 0.;
                }
                if( nsumwindows < iMAXSUMWINDOWS )
                {
                    for( unsigned int s = 0; s < fVLowGainPedvars.size(); s++ )
                    {
                        if( i < fVLowGainPedvars[s].size() )
                        {
                            pedvarLowGainV[s] = fVLowGainPedvars[s][i];
                        }
                    }
                }
                
                ilowgainmultiplier_trace = fLowGainMultiplier_Trace;
                inlowgaindefaultsumwindows = 0;
                
                for( unsigned int k = 0; k < getLowGainDefaultSumWindows().size() && inlowgaindefaultsumwindows < iMAXDEFWINDOWS; k++ )
                {
                    if( getLowGainDefaultSumWindows()[k].first == iTraceIntegrationMethod )
                    {
                        ilowgaindefaultsumwindows[inlowgaindefaultsumwindows] = getLowGainDefaultSumWindows()[k].second ;
                        for( unsigned int j = 1; j <= nsumwindows && j < iMAXSUMWINDOWS; j++ )
                        {
                            ilowgainmultiplier_sum[inlowgaindefaultsumwindows][j] = getLowGainMultiplier_Sum( iTraceIntegrationMethod, getLowGainDefaultSumWindows()[k].second, j );
                        }
                        inlowgaindefaultsumwindows++;
                    }
                }
                
                
                if( i < fLowGainGains.size() )
                {
                    igainlowgain = fLowGainGains[i];
                }
                if( i < fLowGainGains_DefaultSetting.size() )
                {
                    fLowGainGains_DefaultSetting[i] = false;
                }
                if( i < fLowGainGainvars.size() )
                {
                    igainvarlowgain = fLowGainGainvars[i];
                }
                if( i < fLowGainFADCtoPhe.size() )
                {
                    iFADCtoPhelowgain = fLowGainFADCtoPhe[i];
                }
                if( i < fLowGainTOffsets.size() )
                {
                    itofflowgain = fLowGainTOffsets[i];
                }
                if( i < fLowGainTOffsetvars.size() )
                {
                    itoffvarlowgain = fLowGainTOffsetvars[i];
                }
                if( i < fLowGainAverageTzero.size() )
                {
                    itzerolowgain = fLowGainAverageTzero[i];
                }
                if( i < fLowGainAverageTzerovars.size() )
                {
                    itzerovarlowgain = fLowGainAverageTzerovars[i];
                }
                
                fCalibrationTree->Fill();
            }
        }
        fCalibrationTree->Write();
        // do not write histograms -> same information as in fCalibrationTree
        //       hisList->Write();
    }
    
    iDir->cd();
    
    return true;
}


unsigned int VCalibrationData::getTSTimeIndex( double iTime, unsigned int& i1, unsigned int& i2, double& ifrac1, double& ifrac2 )
{
    unsigned iTS_size = getTimeTS_vector().size();
    
    i1 = 0;
    i2 = 0;
    ifrac1 = 1.0;
    ifrac2 = 0.;
    
    if( iTS_size > 0 )
    {
        // time is before first time point
        if( iTime < getTimeTS_vector()[0] )
        {
            i1 = 0;
            i2 = 0;
            ifrac1 = 1.;
            ifrac2 = 0.;
        }
        // time is after first time point
        else if( iTime > getTimeTS_vector()[iTS_size - 1] )
        {
            i1 = iTS_size - 1;
            i2 = iTS_size - 1;
            ifrac1 = 0.;
            ifrac2 = 1.;
        }
        // time is in between
        else
        {
            // get index for closest time point
            unsigned int iIndex = 0;
            double iMaxDiff = 1.e15;
            // time is smaller than closest time point
            for( unsigned int i = 0; i < iTS_size; i++ )
            {
                if( fabs( getTimeTS_vector()[i] - iTime ) < iMaxDiff )
                {
                    iMaxDiff = fabs( getTimeTS_vector()[i] - iTime );
                    iIndex = i;
                }
            }
            // get indexes and fractions
            i1 = iIndex;
            if( iTime - getTimeTS_vector()[iIndex] > 0 )
            {
                if( iIndex < getTimeTS_vector().size() - 1 )
                {
                    i2 = iIndex + 1;
                }
                else
                {
                    i2 = getTimeTS_vector().size() - 1;
                }
            }
            else
            {
                if( iIndex > 0 )
                {
                    i2 = iIndex - 1;
                }
                else
                {
                    i2 = 0;
                }
            }
            
            double id = getTimeTS_vector()[i1] - getTimeTS_vector()[i2];
            if( id != 0. )
            {
                ifrac1 = 1. - ( getTimeTS_vector()[i1] - iTime ) / id;
                ifrac2 = 1. - ( iTime - getTimeTS_vector()[i2] ) / id;
            }
            else
            {
                ifrac1 = 1.;
                ifrac2 = 0.;
            }
        }
    }
    return 0;
}

valarray<double>& VCalibrationData::getPeds( bool iLowGain, double iTime )
{
    // return time dependent pedestals
    if( usePedestalsInTimeSlices( iLowGain ) && iTime > 0. && getPedsTS_vector( iLowGain ).size() > 0. )
    {
        if( fabs( iTime - fTS_ped_temp_time ) < 1.e-3 )
        {
            fTS_ped_temp_time = iTime;
            return fTS_ped_temp;
        }
        fTS_ped_temp_time = iTime;
        
        unsigned int i1 = 0;
        unsigned int i2 = 0;
        double ifrac1 = 0.;
        double ifrac2 = 0.;
        getTSTimeIndex( iTime, i1, i2, ifrac1, ifrac2 );
        
        if( fTS_ped_temp.size() != getPedsTS_vector( iLowGain )[i1].size() )
        {
            fTS_ped_temp.resize( getPedsTS_vector( iLowGain )[i1].size(), 0. );
        }
        
        for( unsigned int i = 0; i < fTS_ped_temp.size(); i++ )
        {
            fTS_ped_temp[i] = ( ifrac1 * getPedsTS_vector( iLowGain )[i1][i] + ifrac2 * getPedsTS_vector( iLowGain )[i2][i] );
        }
        
        return fTS_ped_temp;
    }
    
    if( iLowGain )
    {
        return fLowGainPeds;
    }
    // return pedestal value from DST file
    // (only if at least one event is already read from the DST file)
    if( fReader && fReader->isDST() && fReader->getDSTTreeEvent() > 0 )
    {
        if( fReader->getPedestal().size() == fPeds_perEvent.size() )
        {
            return fReader->getPedestal();
        }
        else
        {
            return fPeds;
        }
    }
    
    return fPeds;
}

valarray<double>& VCalibrationData::getPedvars( bool iLowGain, unsigned int iSW, double iTime )
{
    //////////////////////////////////////////////////////////
    // pedvars in time slices
    //////////////////////////////////////////////////////////
    if( usePedestalsInTimeSlices( iLowGain ) && iTime > 0. && getPedvarsVTS_vector( iLowGain ).size() > 0. )
    {
        // check validity of summation window
        if( iSW > 0 && iSW - 1 < getPedvarsVTS_vector( iLowGain )[0].size() )
        {
            // pedvars is available for this time (maybe filled by a previous call for another channel)
            if( fTS_pedvar_temp_time.find( iSW ) != fTS_pedvar_temp_time.end() && fabs( iTime - fTS_pedvar_temp_time[iSW] ) < 1.e-3 )
            {
                fTS_pedvar_temp_time[iSW] = iTime;
                return fTS_pedvar_temp[iSW];
            }
            // no search for a pedvars for given time
            fTS_pedvar_temp_time[iSW] = iTime;
            
            // get time index
            unsigned int i1 = 0;
            unsigned int i2 = 0;
            double ifrac1 = 0.;
            double ifrac2 = 0.;
            getTSTimeIndex( iTime, i1, i2, ifrac1, ifrac2 );
            
            if( fTS_pedvar_temp[iSW].size() != getPedvarsVTS_vector( iLowGain )[i1][iSW - 1].size() )
            {
                fTS_pedvar_temp[iSW].resize( getPedvarsVTS_vector( iLowGain )[i1][iSW - 1].size(), 0. );
            }
            
            // loop over all channels and calculate pedvars for this time (weighted mean between time bins)
            for( unsigned int i = 0; i < fTS_pedvar_temp[iSW].size(); i++ )
            {
                fTS_pedvar_temp[iSW][i] = ( ifrac1 * getPedvarsVTS_vector( iLowGain )[i1][iSW - 1][i] + ifrac2 * getPedvarsVTS_vector( iLowGain )[i2][iSW - 1][i] );
            }
            
            return fTS_pedvar_temp[iSW];
        }
        else
        {
            cout << "VCalibrationData::getPedvars: invalid summation window for pedestal variations in time slices: " << iSW << " ";
            if( getPedvarsVTS_vector( iLowGain ).size() > 0 )
            {
                cout << getPedvarsVTS_vector( iLowGain )[0].size();
            }
            cout << endl;
            exit( EXIT_FAILURE );
        }
    }
    //////////////////////////////////////////////////////////
    // time independent pedestal variations
    // most DST styles are without time-dependent pedestals -> ped/pedvars vector have length 1
    if( iLowGain )
    {
        if( iSW < fVLowGainPedvars.size() )
        {
            return fVLowGainPedvars[iSW];
        }
        else if( fVLowGainPedvars.size() > 0 )
        {
            return fVLowGainPedvars[0];
        }
    }
    else
    {
        if( iSW < fVPedvars.size() )
        {
            return fVPedvars[iSW];
        }
        else if( fVPedvars.size() > 0 )
        {
            return fVPedvars[0];
        }
    }
    
    // should never arrive here
    return fValArrayDouble;
}


void VCalibrationData::setPeds( unsigned int iChannel, double iPed, bool iLowGain )
{
    if( iLowGain && iChannel < fLowGainPeds.size() )
    {
        fLowGainPeds[iChannel] = iPed;
    }
    
    if( iChannel < fPeds.size() )
    {
        fPeds[iChannel] = iPed;
    }
    
}


double VCalibrationData::getmeanPedvars( bool iLowGain, unsigned int iSW )
{
    if( !iLowGain && iSW < fVmeanPedvars.size() )
    {
        return fVmeanPedvars[iSW];
    }
    
    if( iLowGain && iSW < fVmeanLowGainPedvars.size() )
    {
        return fVmeanLowGainPedvars[iSW];
    }
    
    return 0.;
}


double VCalibrationData::getmeanRMSPedvars( bool iLowGain, unsigned int iSW )
{
    if( !iLowGain && iSW < fVmeanRMSPedvars.size() )
    {
        return fVmeanRMSPedvars[iSW];
    }
    
    if( iLowGain && iSW < fVmeanRMSLowGainPedvars.size() )
    {
        return fVmeanRMSLowGainPedvars[iSW];
    }
    
    return 0.;
}


void VCalibrationData::getmeanPedvars( double& imean, double& irms, bool iLowGain, unsigned int iSW, double iTime )
{
    imean = 0.;
    irms  = 0.;
    
    double its_n = 0.;
    double its_sum = 0.;
    double its_sum2 = 0.;
    
    double iMinPedVars = 2.5;
    if( iSW < 8 )
    {
        iMinPedVars = 0.5;
    }
    
    valarray< double > ipedvar = getPedvars( iLowGain, iSW, iTime );
    
    for( unsigned int i = 0; i < ipedvar.size(); i++ )
    {
        if( ipedvar[i] > iMinPedVars )
        {
            its_n++;
            its_sum += ipedvar[i];
            its_sum2 += ipedvar[i] * ipedvar[i];
        }
    }
    if( its_n > 0. )
    {
        imean = its_sum / its_n;
    }
    if( its_n > 1. )
    {
        irms = sqrt( 1. / ( its_n - 1. ) * ( its_sum2 - 1. / its_n * its_sum * its_sum ) );
    }
}



double VCalibrationData::getLowGainMultiplier_Sum( unsigned int iMethod, int iSumWindow, int jSumWindow )
{
    std::pair<int, int> temp1( iSumWindow, jSumWindow );
    std::pair< unsigned int, std::pair<int, int> > temp2( iMethod, temp1 );
    if( fLowGainMultiplier_Sum.count( temp2 ) == 0 )
    {
        return 0;
    }
    else
    {
        return fLowGainMultiplier_Sum.at( temp2 );
    }
}
/*
The LowGainSumCorrection is used to correct the integrated charge, calculated from the trace * fLowGainMultiplier_Trace, with the low gain multiplier for a given sum window/nominal sum window. It returns 1 (no correction needed) if the low gain multiplier is 0 ( i.e. non-existent) .
*/
double VCalibrationData::getLowGainSumCorrection( unsigned int iMethod, int iSumWindow, int jSumWindow, bool HiLo )
{
    if( !HiLo )
    {
        return 1.0;
    }
    if( getLowGainMultiplier_Sum( iMethod, iSumWindow, jSumWindow ) == 0 || getLowGainMultiplier_Trace() == 0 )
    {
        return 1.0;
    }
    return getLowGainMultiplier_Sum( iMethod, iSumWindow, jSumWindow ) / fLowGainMultiplier_Trace ;
}

bool VCalibrationData::setLowGainMultiplier_Sum( unsigned int iMethod, int iSumWindow, int jSumWindow , double lmult )
{
    std::pair<int, int> temp1( iSumWindow, jSumWindow );
    std::pair< unsigned int, std::pair<int, int> > temp2( iMethod, temp1 );
    fLowGainMultiplier_Sum[temp2] = lmult;
    return true;
}


/*!
    get smallest non-zero pedestal value
*/
double VCalibrationData::getPed_min( bool iLowGain )
{
    double i_min = 1.e6;
    for( unsigned int i = 0; i < getPeds( iLowGain ).size(); i++ )
    {
        if( getPeds( iLowGain )[i] > 0. && getPeds( iLowGain )[i] < i_min )
        {
            i_min = getPeds( iLowGain )[i];
        }
    }
    return i_min;
}

/*!
    get largest non-zero pedestal value
*/
double VCalibrationData::getPed_max( bool iLowGain )
{
    double i_max = 0.;
    for( unsigned int i = 0; i < getPeds( iLowGain ).size(); i++ )
    {
        if( getPeds( iLowGain ) [i] > 0. && getPeds( iLowGain )[i] > i_max )
        {
            i_max = getPeds( iLowGain )[i];
        }
    }
    return i_max;
}

void VCalibrationData::recoverLowGainPedestals()
{
    for( unsigned int i = 0; i < getNSummationWindows(); i++ )
    {
        double iM = 0.;
        double iN = 0.;
        for( unsigned int j = 0; j < getPedvars( true, i ).size(); j++ )
        {
            if( getPedvars( true, i )[j] > 0.1 )
            {
                iM += getPedvars( true, i )[j];
                iN++;
            }
        }
        // recover pedvars
        if( iN > 0. )
        {
            iM /= iN;
        }
        if( iM > 0. )
        {
            for( unsigned int j = 0; j < getPedvars( true, i ).size(); j++ )
            {
                if( getPedvars( true, i )[j] < 0.1 )
                {
                    getPedvars( true, i )[j] = iM;
                }
            }
        }
    }
}

/*
 * get IPR graph for NN image cleaning
 *
 */
TGraphErrors* VCalibrationData::getIPRGraph( unsigned int iSumWindow, bool iMakeNewGraph )
{
    if( fGraphIPRGraph.find( iSumWindow ) != fGraphIPRGraph.end() && fGraphIPRGraph[iSumWindow] )
    {
        return fGraphIPRGraph[iSumWindow];
    }
    else if( iMakeNewGraph )
    {
        fGraphIPRGraph[iSumWindow] = new TGraphErrors( 1 );
        fGraphIPRGraph[iSumWindow]->SetTitle( "" );
        char hname[200];
        sprintf( hname, "IRPFGraph_TelID%d_SumWindow%d", fTelID, iSumWindow );
        fGraphIPRGraph[iSumWindow]->SetName( hname );
        return fGraphIPRGraph[iSumWindow];
    }
    return 0;
}

void VCalibrationData::setIPRGraph( unsigned int iSumWindow, TGraphErrors* g )
{
    fGraphIPRGraph[iSumWindow] = g;
}

/*
 *  get average dc to pe per telescope
 *
 */
double VCalibrationData::getTelescopeAverageFADCtoPhe( bool iLowGain )
{
    double n = 0.;
    double i_n = 0.;
    if( iLowGain )
    {
        for( unsigned int i = 0; i < fLowGainFADCtoPhe.size(); i++ )
        {
            if( fLowGainFADCtoPhe[i] > 0. )
            {
                i_n += fLowGainFADCtoPhe[i];
                n++;
            }
        }
    }
    else
    {
        for( unsigned int i = 0; i < fFADCtoPhe.size(); i++ )
        {
            if( fFADCtoPhe[i] > 0. )
            {
                i_n += fFADCtoPhe[i];
                n++;
            }
        }
    }
    if( n > 0. )
    {
        return i_n / n;
    }
    return 0.;
}
