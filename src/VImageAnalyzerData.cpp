/*! \class VImageAnalyzerData
     \brief all data concerning the channel analysis is stored in this class

*/

#include <VImageAnalyzerData.h>

VImageAnalyzerData::VImageAnalyzerData( unsigned int iTelID, unsigned int iShortTree,
                                        bool bCalibration, bool bWriteImagePixelList )
{
    fTelID = iTelID;
    if( !bCalibration )
    {
        fAnaHistos = new VImageAnalyzerHistograms( iTelID );
        fAnaHistos->init();
    }
    
    fFillMeanTraces = false;
    fFillPulseSum = false;
    
    if( !bCalibration )
    {
        fImageParameter = new VImageParameter( iShortTree, bWriteImagePixelList );
        fImageParameterLogL = new VImageParameter( iShortTree, bWriteImagePixelList );
    }
    
    // initialize time since run start
    fTimeSinceRunStart = -1.;
    fTimeRunStart = 0.;
    
    // random generator for setting randomly channels dead
    fRandomMakeDeadChannelsSeed = 0;
    fRandomMakeDeadChannels = new TRandom3( fRandomMakeDeadChannelsSeed );
    
    fNDead = 0;
    fLowGainNDead = 0;
    
    hMeanPulses = 0;
    hPulseSum = 0;
    
    fSpecialChannel = 0;
    
        fpulsetiming_triggertime_index = 9999;
    fpulsetiming_tzero_index = 9999;
    fpulsetiming_width_index = 9999;
    
    setTraceIntegrationMethod();
}


void VImageAnalyzerData::initialize( unsigned int iChannels, unsigned int iMaxChannel, bool iDebug,
                                     int iseed, unsigned int iSamples, unsigned int ipulsetiminglevel,
			             unsigned int ipulsetiming_tzero_index, unsigned int ipulsetiming_width_index,
                                     unsigned int ipulsetiming_triggertime_index )
{
    if( iDebug )
    {
        cout << "VImageAnalyzerData::initialize" << endl;
    }
    fNChannels = iChannels;
    fMaxChannels = iMaxChannel;
    fNSamples = iSamples;
    fNDead = 0;
    fpulsetiming_tzero_index = ipulsetiming_tzero_index;
    fpulsetiming_width_index = ipulsetiming_width_index;
        fpulsetiming_triggertime_index = ipulsetiming_triggertime_index;
    fLowGainNDead = 0;
    fDead.resize( iChannels, 0 );
    fMasked.resize( iMaxChannel, 0 );
    fDeadRecovered.resize( iChannels, false );
    fDeadUI.resize( iChannels, 0 );
    fLowGainDead.resize( iChannels, 0 );
    fLowGainDeadRecovered.resize( iChannels, false );
    fLowGainDeadUI.resize( iChannels, 0 );
    fTemplateMu.resize( iChannels, 0. );
    fModel3DMu.resize( iChannels, 0. );
    fModel3DClean.resize( iChannels, false );
    fSums.resize( iChannels, 0. );
    fSums2.resize( iChannels, 0. );
    fHiLo.resize( iChannels, false );
    fZeroSuppressed.resize( iChannels, 0 );
    fLLEst.resize( iChannels, false );
    fPulseTimingUncorrected.resize( ipulsetiminglevel, fSums );
    fPulseTimingCorrected.resize( ipulsetiminglevel, fSums );
    fPulseTimingAverageTime.resize( iChannels, 0. );
    fPulseTimingAverageTimeCorrected.resize( iChannels, 0. );
    fTCorrectedSumLast.resize( iChannels, 0 );
    fTCorrectedSumFirst.resize( iChannels, 0 );
    fTCorrectedSum2Last.resize( iChannels, 0 );
    fTCorrectedSum2First.resize( iChannels, 0 );
    fCurrentSummationWindow.resize( iChannels, 0 );
    fCurrentSummationWindow_2.resize( iChannels, 0 );
    fImage.resize( iChannels, false );
    fBorder.resize( iChannels, false );
    fTrigger.resize( iChannels, false );
    fBrightNonImage.resize( iChannels, false );
    fImageBorderNeighbour.resize( iChannels, false );
    fBorderBorderNeighbour.resize( iChannels, false );
    fTraceWidth.resize( iChannels, 0. );
    fTraceMax.resize( iChannels, 0. );
    fTraceN255.resize( iChannels, 0 );
    fRawTraceMax.resize( iChannels, 0. );
    fImageUser.resize( iChannels, 0 );
    
    fClusterID.resize( iChannels, 0 );
    fClusterNpix.resize( iChannels, 0 );
    fClusterSize.resize( iChannels, 0. );
    fClusterTime.resize( iChannels, 0. );
    fClusterCenx.resize( iChannels, 0. );
    fClusterCeny.resize( iChannels, 0. );
    fncluster_cleaned = 0;
    fncluster_uncleaned = 0;
    
    fCorrelationCoefficient.resize( iChannels, false );
    
    fFADCstopTZero.resize( 4, 0. );
    fFADCstopSum.resize( 4, 0. );
    
    fRandomMakeDeadChannelsSeed = iseed;
    fRandomMakeDeadChannels->SetSeed( fRandomMakeDeadChannelsSeed );
}


void VImageAnalyzerData::initializeMeanPulseHistograms()
{
    fFillMeanTraces = true;
    
    // set mean pulse histograms
    hMeanPulses = new TList();
    char hname[200];
    char htitle[200];
    for( unsigned int j = 0; j < fNChannels; j++ )
    {
        sprintf( hname, "hPulseHigh_%d_%u", fTelID + 1, j );
        sprintf( htitle, "mean pulse, high gain (tel %d, channel %u)", fTelID + 1, j );
        hMeanPulseHigh.push_back( new TProfile2D( hname, htitle, 50, 0., 5., 2 * fNSamples, -( double )fNSamples, ( double )fNSamples ) );
        hMeanPulseHigh.back()->SetXTitle( "log_{10} integrated charge" );
        hMeanPulseHigh.back()->SetYTitle( "sample (time corrected)" );
        hMeanPulses->Add( hMeanPulseHigh.back() );
        sprintf( hname, "hPulseLow_%d_%u", fTelID + 1, j );
        sprintf( htitle, "mean pulse, low gain (tel %d, channel %u)", fTelID + 1, j );
        hMeanPulseLow.push_back( new TProfile2D( hname, htitle, 50, 0., 5., 2 * fNSamples, -( double )fNSamples, ( double )fNSamples ) );
        hMeanPulseLow.back()->SetXTitle( "log_{10} integrated charge" );
        hMeanPulseLow.back()->SetYTitle( "sample (time corrected)" );
        hMeanPulses->Add( hMeanPulseLow.back() );
    }
}


void VImageAnalyzerData::initializeIntegratedChargeHistograms()
{
    fFillPulseSum = true;
    
    hPulseSum = new TList();
    char hname[200];
    char htitle[200];
    for( unsigned int j = 0; j < fNChannels; j++ )
    {
        sprintf( hname, "hSumHigh_%d_%u", fTelID + 1, j );
        sprintf( htitle, "integrated charge, high gain (tel %d, channel %u)", fTelID + 1, j );
        hPulseSumHigh.push_back( new TH1F( hname, htitle, 200, 0., 6. ) );
        hPulseSumHigh.back()->SetXTitle( "log_{10} integrated charge" );
        hPulseSum->Add( hPulseSumHigh.back() );
        
        sprintf( hname, "hSumLow_%d_%u", fTelID + 1, j );
        sprintf( htitle, "integrated charge, low gain (tel %d, channel %u)", fTelID + 1, j );
        hPulseSumLow.push_back( new TH1F( hname, htitle, 200, 0., 6. ) );
        hPulseSumLow.back()->SetXTitle( "log_{10} integrated charge" );
        hPulseSum->Add( hPulseSumLow.back() );
    }
}


void VImageAnalyzerData::setTrace( unsigned int iChannel, vector< double > fT, bool iHiLo, double iPeds )
{
    if( fFillMeanTraces && iChannel < hMeanPulseLow.size() && iChannel < hMeanPulseHigh.size() )
    {
        // expect that trace sum is set before!!!
        if( fSums[iChannel] > 0. && getTZeros( true )[iChannel] > 0 )
        {
            for( unsigned t = 0; t < fT.size(); t++ )
            {
                if( iHiLo )
                {
                    hMeanPulseLow[iChannel]->Fill( log10( fSums[iChannel] ), ( double )( t ) - getTZeros( true )[iChannel], fT[t] - iPeds );
                }
                else
                {
                    hMeanPulseHigh[iChannel]->Fill( log10( fSums[iChannel] ), ( double )( t ) - getTZeros( true )[iChannel], fT[t] - iPeds );
                }
            }
        }
    }
}


void VImageAnalyzerData::fillPulseSum( unsigned int iChannel, double iS, bool iHiLo )
{
    if( iS > 0. )
    {
        if( !iHiLo && iChannel < hPulseSumHigh.size() )
        {
            hPulseSumHigh[iChannel]->Fill( log10( iS ) );
        }
        if( iHiLo && iChannel < hPulseSumLow.size() )
        {
            hPulseSumLow[iChannel]->Fill( log10( iS ) );
        }
    }
}

vector<unsigned int>& VImageAnalyzerData::getFADCstopTrigChannelID()
{
    if( fSpecialChannel )
    {
        return fSpecialChannel->getFADCstopTrigChannelID();
    }
    
    return iDummyVectorUI;
}

bool VImageAnalyzerData::readSpecialChannels( int iRunNumber, string iFile, string iDirectory )
{
    if( fSpecialChannel )
    {
        fSpecialChannel->reset();
    }
    else
    {
        fSpecialChannel = new VSpecialChannel( fTelID );
    }
    fSpecialChannel->readSpecialChannels( iRunNumber, iFile, iDirectory );
    
    return !fSpecialChannel->isZombie();
}

double VImageAnalyzerData::getHIGHQE_gainfactor( unsigned int iChannel )
{
    if( fSpecialChannel )
    {
        return fSpecialChannel->getHIGHQE_gainfactor( iChannel );
    }
    
    return -1.;
}

valarray<double>& VImageAnalyzerData::getTraceAverageTime( bool iCorrected )
{
     if( iCorrected ) return fPulseTimingAverageTimeCorrected;

     return fPulseTimingAverageTime;
}

valarray<double>& VImageAnalyzerData::getTTrigger()
{
	// return tzero according to pulse timing vector
	if( fpulsetiming_triggertime_index < fPulseTimingCorrected.size() )
	{
		return fPulseTimingCorrected[fpulsetiming_triggertime_index];
	}
	
	// this is a serious problem and should never happen
	cout << "VImageAnalyzerData::getTTrigger error: trigger time index out of range " << endl;
	cout << "\t" << fpulsetiming_triggertime_index << "\t" << fPulseTimingCorrected.size() << endl;
	exit( EXIT_FAILURE );
	
	return fPulseTimingCorrected[0]; // should never happen
}


valarray<double>& VImageAnalyzerData::getTZeros( bool iCorrected )
{
    // return tzero according to pulse timing vector
    if( iCorrected && fpulsetiming_tzero_index < fPulseTimingCorrected.size() )
    {
        return fPulseTimingCorrected[fpulsetiming_tzero_index];
    }
    else if( fpulsetiming_tzero_index < fPulseTimingUncorrected.size() )
    {
        return fPulseTimingUncorrected[fpulsetiming_tzero_index];
    }
    
    // this is a serious problem and should never happen
    cout << "VImageAnalyzerData::getTZeros error: tzero index out of range" << endl;
    cout << "\t" << fpulsetiming_tzero_index << "\t" << fPulseTimingCorrected.size() << "\t" << fPulseTimingUncorrected.size() << endl;
    exit( EXIT_FAILURE );
    
    return fPulseTimingUncorrected[0]; // should never happen
}

valarray<double>& VImageAnalyzerData::getTraceWidth( bool iCorrected )
{
    if( iCorrected )
    {
        if( fpulsetiming_width_index < fPulseTimingCorrected.size() )
        {
            return fPulseTimingCorrected[fpulsetiming_width_index];
        }
    }
    else
    {
        if( fpulsetiming_width_index < fPulseTimingUncorrected.size() )
        {
            return fPulseTimingUncorrected[fpulsetiming_width_index];
        }
    }
    
    // this is a serious problem and should never happen
    cout << "VImageAnalyzerData::getTraceWidth error: tzero or width index out of range" << endl;
    cout << fpulsetiming_width_index << "\t";
    cout << fPulseTimingCorrected.size() << "\t" << fPulseTimingUncorrected.size() << endl;
    exit( -1 );
    
    return fPulseTimingUncorrected[0]; // should never happen
}
