/*! \class VBaseRawDataReader

    contains common elements for VRawDataReader and VBFDataReader


*/

#include <VRawDataReader.h>

VBaseRawDataReader::VBaseRawDataReader( string sourcefile, int isourcetype, unsigned int iNTel, bool iDebug )
{
    fDebug = iDebug;
    if( fDebug )
    {
        cout << "VBaseRawDataReader::VBaseRawDataReader" << endl;
    }
    if( fDebug )
    {
        setSimuDebugFlag();
    }
    fDataFormat = "rawdata";
    fDataFormatNum = 0;
    fEventNumber = 0;
    fNTel = iNTel;
    fTelID = 0;
    fMonteCarloHeader = 0;
    // noise file reader values for external noise files
    fNoiseFileReader = 0;
    fNoiseFilePedestal = 0;
    fNoiseFileFADCRange = 250;
    
    // additional Gaussian noise
    finjectGaussianNoise = -1.;
    fRandomInjectGaussianNoise = 0;
    
    // source types
    if( isourcetype == 2 )
    {
        fDataFormat = "MCvbf";
        fDataFormatNum = 2;
    }
    else if( isourcetype == 3 )
    {
        fDataFormat = "Rawvbf";
        fDataFormatNum = 3;
    }
    fSourceFileName = sourcefile;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fEvent.push_back( 0 );
    }
}


void VBaseRawDataReader::setDebug( bool iDebug )
{
    fDebug = iDebug;
    if( fDebug )
    {
        setSimuDebugFlag();
    }
}


bool VBaseRawDataReader::isMC()
{
    if( fDataFormat == "MCvbf" )
    {
        return true;
    }
    
    return false;
}


VBaseRawDataReader::~VBaseRawDataReader()
{
    if( fDebug )
    {
        cout << "VBaseRawDataReader::~VBaseRawDataReader()" << endl;
    }
    /*   delete fEvent; */
}


/*!
    1 = rawdata
    3 = MCvbf
*/
unsigned int VBaseRawDataReader::getDataFormatNum()
{
    return fDataFormatNum;
}


bool VBaseRawDataReader::setTelescopeID( unsigned int iTelID )
{
    if( iTelID < fNTel )
    {
        fTelID = iTelID;
    }
    else
    {
        fTelID = 0;
        return false;
    }
    return true;
}


uint32_t VBaseRawDataReader::getHitID( uint32_t i )
{
    if( !fEvent[fTelID] )
    {
        return 0;
    }
    uint32_t t_hitid = 0;
    try
    {
        t_hitid = fEvent[fTelID]->getHitID( i );
    }
    catch( VIndexOutOfBoundsException& ex )
    {
        cout << "VRawDataReader::getHitID error " << ex << endl;
        return 0;
    }
    return t_hitid;
}


uint16_t VBaseRawDataReader::getMaxChannels()
{
    if( !fEvent[fTelID] )
    {
        return 499;
    }
    if( fTelID < fEvent.size() )
    {
        if( fEvent[fTelID] )
        {
            return fEvent[fTelID]->getMaxChannels();
        }
    }
    else
    {
        cout << "VBaseRawDataReader::getMaxChannels(): problem with event size " << fTelID << "\t" << fEvent.size() << endl;
    }
    return 499;
}


uint16_t VBaseRawDataReader::getNumSamples()
{
    if( !fEvent[fTelID] )
    {
        return 0;
    }
    if( fTelID < fEvent.size() )
    {
        if( fEvent[fTelID] )
        {
            return fEvent[fTelID]->getNumSamples();
        }
    }
    return 0;
}


uint8_t VBaseRawDataReader::getNewEventType( unsigned int itelID )
{
    if( !fEvent[itelID] )
    {
        return 0;
    }
    
    if( itelID < fEvent.size() )
    {
        return fEvent[itelID]->getEventType().getBestNewStyleCode();
    }
    
    return 0;
}


uint8_t VBaseRawDataReader::getEventType()
{
    if( fEvent[fTelID] )
    {
        return fEvent[fTelID]->getEventType().getBestOldStyleCode();
    }
    
    return 0;
}


uint8_t VBaseRawDataReader::getNewEventType()
{
    if( fEvent[fTelID] )
    {
        return fEvent[fTelID]->getEventType().getBestNewStyleCode();
    }
    
    return 0;
}


int VBaseRawDataReader::getNumberofFullTrigger()
{
    int z = 0;
    if( fEvent[fTelID] )
    {
        unsigned int i_maxchannel = getMaxChannels();
        for( unsigned i = 0; i < i_maxchannel; i++ )
        {
            if( fEvent[fTelID]->getTriggerBit( i ) )
            {
                z++;
            }
        }
    }
    return z;
}


uint16_t VBaseRawDataReader::getGPSYear()
{
    if( fEvent[fTelID] )
    {
        if( getRunNumber() < 33242 || getRunNumber() > 33253 )
        {
            return fEvent[fTelID]->getGPSYear();
        }
        else
        {
            return fEvent[fTelID]->getGPSYear() + 1;
        }
    }
    return 50;
}


bool VBaseRawDataReader::wasLossyCompressed()
{
#ifdef VBF_027
    if( fEvent[fTelID] )
    {
        return fEvent[fTelID]->wasLossyCompressed();
    }
    else
    {
        return false;
    }
#endif
    
    return false;
}


bool VBaseRawDataReader::getHiLo( uint32_t i )
{
    if( fEvent[fTelID] )
    {
        try
        {
            return fEvent[fTelID]->getHiLo( i );
        }
        catch( VException& e )
        {
            cout << "EXCEPTION " << i << endl;
            return false;
        }
    }
    
    return false;
}

/*
 * inject Gaussian noise
 *
 * set std dev in units of dc
 *
 */

void VBaseRawDataReader::injectGaussianNoise( double injectGaussianNoise,  UInt_t seed )
{
    finjectGaussianNoise = injectGaussianNoise;
    
    if( !fRandomInjectGaussianNoise )
    {
        fRandomInjectGaussianNoise = new TRandom3( seed );
    }
}

/*
 *
 * initialise trace noise generator
 *
 * reading from external noise library
 *
 */
bool VBaseRawDataReader::initTraceNoiseGenerator( unsigned int iType, string iT, VDetectorGeometry* iD, vector<int> iSW,
        bool iDebug, int iseed, double iDefaultPed, vector<double> iFADCCorrect )
{
    if( fDebug )
    {
        cout << "VBaseRawDataReader::initTraceNoiseGenerator " << endl;
    }
    fNoiseFileReader = new VNoiseFileReader( iType, iT );
    fNoiseFilePedestal = ( uint8_t )iDefaultPed;
    if( iD && iD->getFADCRange() < 255 )
    {
        fNoiseFileFADCRange = ( uint8_t )iD->getFADCRange();
    }
    
    // preliminary: use value from Telescope 1 for all telescopes
    double iCorrection = 1.;
    if( iFADCCorrect.size() > 0 && iFADCCorrect[0] > 0. )
    {
        iCorrection = 1.;
        if( TMath::Abs( iCorrection - 1. ) > 0.01 && iCorrection > 0. )
        {
            cout << "init trace noise generator: use gain correction from telescope 1 for all telescope (" << 1. / iCorrection << ")" << endl;
        }
    }
    bool iB = fNoiseFileReader->init( iD, iD->getNumTelescopes(), iSW, iDebug, iseed, iCorrection );
    fNoiseFileReader->setDefaultGrisuPed( iDefaultPed );
    
    return iB;
}

/*
 * is this channel zero suppressed:
 *
 * true: channel is zero suppressed
 */
bool VBaseRawDataReader::isZeroSuppressed( unsigned int channel )
{
    pair< bool, uint32_t > i_hitIndexPair = getChannelHitIndex( channel );
    
    if( !i_hitIndexPair.first )
    {
        return true;
    }
    
    return false;
}

/*
 * return type of zero suppression
 *
 * 0 = not surpressed
 * bit 0: charge suppressed (but samples available)
 * bit 1: samples suppressed (but charge available)
 * bit 0+1: completely suppressed
 *
 * for VERITAS raw data: charge and samples are suppressed
 */
unsigned short int VBaseRawDataReader::getZeroSuppressionFlag( unsigned int channel )
{
    pair< bool, uint32_t > i_hitIndexPair = getChannelHitIndex( channel );
    
    if( !i_hitIndexPair.first )
    {
        return 3;
    }
    
    return 0;
}


uint8_t VBaseRawDataReader::getSample( unsigned channel, unsigned sample, bool iNewNoiseTrace )
{
    uint8_t iSampleValue = 0;
    try
    {
        if( fEvent[fTelID] )
        {
            iSampleValue = fEvent[fTelID]->getSample( channel, sample );
        }
    }
    catch( ... )
    {
        cout << "VBaseRawDataReader::getSample error: failed for channel " << channel << " and sample " << sample << endl;
        return 0;
    }
    
    // add noise from external noise library to traces
    // (e.g. VTS grisu MC are simulated without noise, noise is added here to the samples)
    if( fNoiseFileReader && !getHiLo( channel ) )
    {
        uint8_t iNoiseSampleValue = fNoiseFileReader->getNoiseSample( fTelID, channel, sample, iNewNoiseTrace );
        if( iSampleValue > iNoiseSampleValue && iSampleValue > fNoiseFileFADCRange - iNoiseSampleValue + fNoiseFilePedestal )
        {
            return fNoiseFileFADCRange;
        }
        return iSampleValue + iNoiseSampleValue - fNoiseFilePedestal;
    }
    // add gaussian noise
    if( finjectGaussianNoise > 0. && fRandomInjectGaussianNoise && !getHiLo( channel ) )
    {
        double iNoiseGaus = fRandomInjectGaussianNoise->Gaus(  0., finjectGaussianNoise );
        if( iSampleValue + iNoiseGaus < 256 )
        {
            return iSampleValue + iNoiseGaus;
        }
        else
        {
            return iSampleValue;
        }
    }
    return iSampleValue;
}


std::vector< uint8_t > VBaseRawDataReader::getSamplesVec()
{
    // standard way
    if( !fNoiseFileReader )
    {
        if( fEvent[fTelID] )
        {
            return fEvent[fTelID]->getSamplesVec();
        }
        else
        {
            return fDummyUint8V;
        }
    }
    else
    {
        if( fEvent[fTelID] )
        {
            std::vector< uint8_t > i_temp = fEvent[fTelID]->getSamplesVec();
            std::vector< uint8_t > i_pedV = fNoiseFileReader->getNoiseVec( fTelID, fHitID );
            if( i_temp.size() == i_pedV.size() )
            {
                for( unsigned int i = 0; i < i_temp.size(); i++ )
                {
                    if( i_temp[i] > fNoiseFileFADCRange - i_pedV[i] + fNoiseFilePedestal )
                    {
                        i_temp[i] = fNoiseFileFADCRange;
                    }
                    else
                    {
                        i_temp[i] += i_pedV[i] - fNoiseFilePedestal;
                    }
                }
            }
            return i_temp;
        }
    }
    return fDummyUint8V;
}


void VBaseRawDataReader::selectHitChan( uint32_t i )
{
    fHitID = i;
    
    if( fEvent[fTelID] )
    {
        return fEvent[fTelID]->selectHitChan( i );
    }
    
    return;
}


void VBaseRawDataReader::setSumWindow( unsigned int iTelID, int isw )
{
    if( fNoiseFileReader )
    {
        fNoiseFileReader->setSumWindow( iTelID, isw );
    }
}


valarray<double>& VBaseRawDataReader::getPeds()
{
    if( fNoiseFileReader )
    {
        return fNoiseFileReader->getPeds();
    }
    
    return v;
}


valarray<double>& VBaseRawDataReader::getPedvars()
{
    if( fNoiseFileReader )
    {
        return fNoiseFileReader->getPedvars();
    }
    
    return v;
}


vector< valarray<double> >& VBaseRawDataReader::getPedvarsAllSumWindows()
{
    if( fNoiseFileReader )
    {
        return fNoiseFileReader->getPedvarsAllSumWindows();
    }
    
    return vv;
}


valarray<double>& VBaseRawDataReader::getPedRMS()
{
    if( fNoiseFileReader )
    {
        return fNoiseFileReader->getPedRMS();
    }
    
    return v;
}

bool VBaseRawDataReader::hasFADCTrace()
{
    if( fTelID < fEvent.size() && fEvent[fTelID] )
    {
        if( fEvent[fTelID]->getNumSamples() == 0 )
        {
            return false;
        }
    }
    return true;
}

valarray< double >& VBaseRawDataReader::getSums( unsigned int iNChannel )
{
    if( iNChannel != 99999 && fSums.size() != iNChannel )
    {
        fSums.resize( iNChannel );
    }
    
    if( fTelID < fEvent.size() && fEvent[fTelID] )
    {
        fSums = 0.;
        for( unsigned int i = 0; i < fEvent[fTelID]->getNumChannelsHit(); i++ )
        {
            unsigned int i_channelHitID = getHitID( i );
            if( i_channelHitID < fSums.size() )
            {
                fSums[i_channelHitID] = fEvent[fTelID]->getCharge( i );
            }
        }
    }
    else
    {
        fSums = 0.;
    }
    return fSums;
}

valarray< double >& VBaseRawDataReader::getTraceMax( unsigned int iNChannel )
{
    if( iNChannel != 99999 && fTraceMax.size() != iNChannel )
    {
        fTraceMax.resize( iNChannel );
    }
    fTraceMax = 0.;
    
    return fTraceMax;
}

vector< valarray< double > >& VBaseRawDataReader::getTracePulseTiming( unsigned int iNChannel )
{
    if( fTelID < fEvent.size() && fEvent[fTelID] )
    {
        // check only first entry (anyway a dummy vector)
        if( fTracePulseTiming.size() == VDST_MAXTIMINGLEVELS && fTracePulseTiming[0].size() == iNChannel )
        {
            return fTracePulseTiming;
        }
        valarray< double > iTemp( 0., iNChannel );
        fTracePulseTiming.clear();
        for( unsigned int t = 0; t < VDST_MAXTIMINGLEVELS; t++ )
        {
            fTracePulseTiming.push_back( iTemp );
        }
    }
    return fTracePulseTiming;
    
}
