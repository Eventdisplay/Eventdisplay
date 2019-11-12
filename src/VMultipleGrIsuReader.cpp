/*! \class VMultipleGrIsuReader
    \brief read grisu simulation data from several files (one per telescope)


*/

#include "VMultipleGrIsuReader.h"

VMultipleGrIsuReader::VMultipleGrIsuReader( unsigned int nFiles, vector< unsigned int > iTelToAna, bool iDebug )
{
    fDebug = iDebug;
    
    fNFiles = nFiles;
    fTelescopeID = 0;
    fFileID = 0;
    fTeltoAna = iTelToAna;
    
    setEventStatus( 1 );
    
    for( unsigned int i = 0; i < fNFiles; i++ )
    {
        fSelectedTelescope.push_back( true );
        fLocalTrigger.push_back( false );
    }
    fNoiseFileReader = 0;
    // initialize pointing vectors
    fTelElevation.assign( fNFiles, 0. );
    fTelAzimuth.assign( fNFiles, 0. );
}


bool VMultipleGrIsuReader::init( VDetectorGeometry* iD, string i_sourcefile, vector< int > i_sumwindow, int i_telnumberoffset, int i_sampleoffset, double ifadcscale, int iseed, string iExPedFile, bool iSingleExternalPedFile, double iDefaultPed )
{
    char hname[2000];
    bool iB = true;
    
    map< unsigned int, unsigned int > iTelGrisu = iD->getTelIDGrisu();
    
    //////////////////////////////////////////////////////////////////
    // make file names
    cout << "VMultipleGrIsuReader::init source file: " << fNFiles << endl;
    for( unsigned int i = 0; i < fNFiles; i++ )
    {
        // check if this telescope should be analyzed
        bool bAna = false;
        for( unsigned int t = 0; t < fTeltoAna.size(); t++ ) if( fTeltoAna[t] == i )
            {
                bAna = true;
            }
            
        // make file names for source file
        if( bAna )
        {
            // get telescope ID mixing map
            // (this is for the case when N telescope are simulated in corsika and grisu, but
            //  eventdisplay should analyse a subset only)
            if( iTelGrisu.find( i ) != iTelGrisu.end() )
            {
                sprintf( hname, "%s.T%d", i_sourcefile.c_str(), iTelGrisu[i] + 1 );
            }
            else
            {
                sprintf( hname, "%s.T%d", i_sourcefile.c_str(), i + 1 );
            }
            cout << "VMultipleGrIsuReader::init sourcefile for telescope " << i + 1 << ": " << hname << endl;
            fSourceFileName.push_back( hname );
        }
        else
        {
            fSourceFileName.push_back( "" );
        }
    }
    
    /////////////////////////////////////////////////////////
    // create noise file reader
    // (one single noise file for all grisu files)
    if( iSingleExternalPedFile )
    {
        fNoiseFileReader = new VNoiseFileReader( 0, iExPedFile );
        fNoiseFileReader->init( iD, 1, i_sumwindow, fDebug, iseed );
        fNoiseFileReader->setDefaultGrisuPed( iDefaultPed );
    }
    /////////////////////////////////////////////////////////
    // create all the readers
    for( unsigned int i = 0; i < fNFiles; i++ )
    {
        if( fSourceFileName[i].size() > 0 )
        {
            if( !iSingleExternalPedFile )
            {
                fReader.push_back( new VGrIsuReader( iD, 1, fSourceFileName[i], i_sumwindow, i_telnumberoffset, i_sampleoffset, ifadcscale, fDebug, iseed, iExPedFile ) );
            }
            else
            {
                fReader.push_back( new VGrIsuReader( iD, 1, fSourceFileName[i], i_sumwindow, i_telnumberoffset, i_sampleoffset, ifadcscale, fDebug, iseed, "" ) );
                // fill random pedestals for all the readers
                fillRandomPeds( fReader.back(), iseed * 3 );
            }
            // reading from multiple files
            fReader.back()->setMultiGrIsuReader( true );
            // one telescope per file
            fReader.back()->setTelescopeID( 0 );
        }
        else
        {
            fReader.push_back( 0 );
        }
    }
    return iB;
}


void VMultipleGrIsuReader::fillRandomPeds( VGrIsuReader* g, int iseed )
{
    if( ! g )
    {
        return;
    }
    if( !fNoiseFileReader )
    {
        return;
    }
    
    cout << "\t VMultipleGrIsuReader: randomizing pedestals" << endl;
    
    TRandom3 iRandom( iseed );
    
    if( fNoiseFileReader->getFullNoiseVec().size() > 0 && fNoiseFileReader->getFullNoiseVec()[0].size() > 0 )
    {
        g->assignGrisuPeds( fNoiseFileReader->getFullNoiseVec()[0][0].size() );
    }
    
    // randomize all values
    if( fNoiseFileReader->getPeds().size() == g->getPeds().size() && g->getPedvars().size() == fNoiseFileReader->getPedvars().size() && g->getPedvarsAllSumWindows().size() == fNoiseFileReader->getPedvarsAllSumWindows().size() && g->getPedRMS().size() == fNoiseFileReader->getPedRMS().size() )
    {
        for( unsigned int i = 0; i < fNoiseFileReader->getPeds().size(); i++ )
        {
            int f = iRandom.Integer( fNoiseFileReader->getPeds().size() );
            
            g->getPeds()[i] = fNoiseFileReader->getPeds()[f];
            g->getPedvars()[i] = fNoiseFileReader->getPedvars()[f];
            for( unsigned int w = 0; w < g->getPedvarsAllSumWindows().size(); w++ )
            {
                g->getPedvarsAllSumWindows()[w][i] = fNoiseFileReader->getPedvarsAllSumWindows()[w][f];
            }
            g->getPedRMS()[i] = fNoiseFileReader->getPedRMS()[i];
            
            for( unsigned int b = 0; b < fNoiseFileReader->getFullNoiseVec( 0, f ).size(); b++ )
            {
                g->getFullNoiseVec( 0, i )[b] = fNoiseFileReader->getFullNoiseVec( 0, f )[b];
            }
        }
    }
    else
    {
        cout << "ERROR: VMultipleGrIsuReader::fillRandomPeds, pedestal vectors of different lengths" << endl;
        return;
    }
}


bool VMultipleGrIsuReader::checkTelescopeID( unsigned int iT )
{
    if( iT < fNFiles && iT < fReader.size() )
    {
        return true;
    }
    
    cout << "VMultipleGrIsuReader::checkTelescopeID; ERROR: telescope ID out of range: " << iT << " " << fNFiles << " " << fReader.size() << endl;
    
    return false;
}


VGrIsuReader* VMultipleGrIsuReader::getReader()
{
    if( fReader.size() > 0 && fReader[fTelescopeID] )
    {
        return fReader[fTelescopeID];
    }
    
    return 0;
}


string VMultipleGrIsuReader::getDataFormat()
{
    if( getReader() )
    {
        return getReader()->getDataFormat();
    }
    
    string a;
    return a;
}


string VMultipleGrIsuReader::getSourceFileName()
{
    if( fSourceFileName.size() > 0 )
    {
        return fSourceFileName[fTelescopeID];
    }
    
    string a;
    return a;
}


std::pair< bool, uint32_t > VMultipleGrIsuReader::getChannelHitIndex( uint32_t t )
{
    if( getReader() )
    {
        return getReader()->getChannelHitIndex( t );
    }
    
    pair< bool, uint32_t > a;
    return a;
}


uint32_t VMultipleGrIsuReader::getEventNumber()
{
    if( getReader() )
    {
        return getReader()->getEventNumber();
    }
    
    return 0;
}


// don't know the difference between hit and trig
std::vector< bool > VMultipleGrIsuReader::getFullHitVec()
{
    if( getReader() )
    {
        return getReader()->getFullHitVec();
    }
    
    vector< bool > a;
    return a;
}


std::vector< bool > VMultipleGrIsuReader::getFullTrigVec()
{
    if( getReader() )
    {
        return getReader()->getFullTrigVec();
    }
    
    vector< bool > a;
    return a;
}


int VMultipleGrIsuReader::getNumberofFullTrigger()
{
    if( getReader() )
    {
        return getReader()->getNumberofFullTrigger();
    }
    
    return 0;
}


std::vector< int > VMultipleGrIsuReader::getFullAnaVec()
{
    if( getReader() )
    {
        return getReader()->getFullAnaVec();
    }
    
    vector< int > a;
    return a;
}


uint8_t VMultipleGrIsuReader::getEventType()
{
    if( getReader() )
    {
        return getReader()->getEventType();
    }
    
    return 0;
}


uint8_t VMultipleGrIsuReader::getATEventType()
{
    if( getReader() )
    {
        return getReader()->getATEventType();
    }
    
    return 0;
}


uint32_t VMultipleGrIsuReader::getHitID( uint32_t t )
{
    if( getReader() )
    {
        return getReader()->getHitID( t );
    }
    
    return 0;
}


bool VMultipleGrIsuReader::getHiLo( uint32_t i )
{
    if( getReader() )
    {
        return getReader()->getHiLo( i );
    }
    
    return 0;
}


uint16_t VMultipleGrIsuReader::getMaxChannels()
{
    if( getReader() )
    {
        return getReader()->getMaxChannels();
    }
    
    return 0;
}


std::vector< uint8_t > VMultipleGrIsuReader::getNoiseVec( unsigned int iTel, uint32_t iHitID )
{
    if( checkTelescopeID( iTel ) && fReader[iTel] )
    {
        return fReader[iTel]->getNoiseVec( 0, iHitID );
    }
    
    std::vector< uint8_t > a;
    return a;
}


uint16_t VMultipleGrIsuReader::getNumChannelsHit()
{
    if( getReader() )
    {
        return getReader()->getNumChannelsHit();
    }
    
    return 0;
}


uint16_t VMultipleGrIsuReader::getNumSamples()
{
    if( getReader() )
    {
        return getReader()->getNumSamples();
    }
    
    return 0;
}


TH1F* VMultipleGrIsuReader::getPedHisto( unsigned int iTel, unsigned int ichannel )
{
    if( checkTelescopeID( iTel ) && fReader[iTel] )
    {
        return fReader[iTel]->getPedHisto( 0, ichannel );
    }
    
    return 0;
}


std::valarray< double >& VMultipleGrIsuReader::getPeds()
{
    if( getReader() )
    {
        return getReader()->getPeds();
    }
    
    return vvv_valarray;
}


std::valarray< double >& VMultipleGrIsuReader::getPedvars()
{
    if( getReader() )
    {
        return getReader()->getPedvars();
    }
    
    return vvv_valarray;
}


std::vector< valarray<double> >& VMultipleGrIsuReader::getPedvarsAllSumWindows()
{
    if( getReader() )
    {
        return getReader()->getPedvarsAllSumWindows();
    }
    
    return vvv_v_vvv_valarray;
}


std::valarray< double >& VMultipleGrIsuReader::getPedRMS()
{
    if( getReader() )
    {
        return getReader()->getPedRMS();
    }
    
    return vvv_valarray;
}


std::vector< uint8_t > VMultipleGrIsuReader::getSamplesVec()
{
    if( getReader() )
    {
        return getReader()->getSamplesVec();
    }
    
    std::vector< uint8_t > a;
    return a;
}


std::vector< double > VMultipleGrIsuReader::getTelElevation()
{
    return fTelElevation;
}


std::vector< double > VMultipleGrIsuReader:: getTelAzimuth()
{
    return fTelAzimuth;
}


void VMultipleGrIsuReader::selectHitChan( uint32_t t )
{
    if( getReader() )
    {
        getReader()->selectHitChan( t );
    }
}


void VMultipleGrIsuReader::setPedestalEventMode( int i_Nped )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setPedestalEventMode( i_Nped );
        }
    }
}


void VMultipleGrIsuReader::setRandomDead( int iNC, int iNB )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setRandomDead( iNC, iNB );
        }
    }
}


void VMultipleGrIsuReader::setSumWindow( vector< int > i_sum )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setSumWindow( i_sum );
        }
    }
}


void VMultipleGrIsuReader::setSumWindow( unsigned int iTelID, int isw )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setSumWindow( iTelID, isw );
        }
    }
}


void VMultipleGrIsuReader::setTelNumberOffset( int iOff )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setTelNumberOffset( iOff );
        }
    }
}


bool VMultipleGrIsuReader::setTelescopeID( unsigned int t )
{
    fTelescopeID = t;
    return checkTelescopeID( t );
}


void VMultipleGrIsuReader::setTraceFile( string i_Trf )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setTraceFile( i_Trf );
        }
    }
}


void VMultipleGrIsuReader::setTrigger( vector<bool> iImage, vector<bool> iBorder )
{
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fReader[i]->setTrigger( iImage, iBorder );
        }
    }
}


vector< bool >& VMultipleGrIsuReader::getLocalTrigger()
{
    return fLocalTrigger;
}


float VMultipleGrIsuReader::getLocalTriggerTime( unsigned int iTel )
{
    if( checkTelescopeID( iTel ) && fReader[iTel] )
    {
        return fReader[iTel]->getLocalTriggerTime( 0 );
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getLocalDelayedTriggerTime( unsigned int iTel )
{
    if( checkTelescopeID( iTel ) && fReader[iTel] )
    {
        return fReader[iTel]->getLocalDelayedTriggerTime( 0 );
    }
    
    return 0.;
}


unsigned int VMultipleGrIsuReader::getNTelLocalTrigger()
{
    if( getReader() )
    {
        return getReader()->getNTelLocalTrigger();
    }
    
    return 0;
}


bool VMultipleGrIsuReader::hasArrayTrigger()
{
    // require at least two telescopes with local trigger
    unsigned int ntrig = 0;
    
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        if( fReader[i] && fReader[i]->hasLocalTrigger( 0 ) )
        {
            ntrig++;
        }
        if( ntrig > 1 )
        {
            return true;
        }
    }
    
    return false;
}


bool VMultipleGrIsuReader::hasLocalTrigger( unsigned int iTel )
{
    if( checkTelescopeID( iTel ) && fReader[iTel] )
    {
        return fReader[iTel]->hasLocalTrigger( 0 );
    }
    
    return false;
}


double VMultipleGrIsuReader::getXimpactrot()
{
    if( getReader() )
    {
        return getReader()->getXimpactrot();
    }
    
    return 0.;
}


double VMultipleGrIsuReader::getYimpactrot()
{
    if( getReader() )
    {
        return getReader()->getYimpactrot();
    }
    
    return 0.;
}


void VMultipleGrIsuReader::setDefaultPed( double iD )
{
    if( getReader() )
    {
        return getReader()->setDefaultPed( iD );
    }
}


bool VMultipleGrIsuReader::getNextEvent()
{
    bool iB = true;
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        /////////////////////////////////////////////////////////////////
        // get next event from each file
        bool iC = false;
        if( fReader[i] && fSelectedTelescope[i] )
        {
            iC = fReader[i]->getNextEvent();
        }
        iB = ( iC && iB );
        /////////////////////////////////////////////////////////////////
        // set the local trigger for this telescope
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fLocalTrigger[i] = fReader[i]->getLocalTrigger()[0];
        }
        else
        {
            fLocalTrigger[i] = 0;
        }
        /////////////////////////////////////////////////////////////////
        // get telescope poiting
        if( fReader[i] && fSelectedTelescope[i] )
        {
            fTelElevation[i] = fReader[i]->getTelElevation()[0];
            fTelAzimuth[i]   = fReader[i]->getTelAzimuth()[0];
        }
        else
        {
            fTelElevation[i] = 0.;
            fTelAzimuth[i]   = 0.;
        }
    }
    /////////////////////////////////////////////////////////////
    // set eventstatuts to 999 for the case of no success
    if( !iB )
    {
        setEventStatus( 999 );
    }
    return iB;
}


bool VMultipleGrIsuReader::getNextPedestalEvent()
{
    bool iB = true;
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        bool iC = false;
        if( fReader[i] && fSelectedTelescope[i] )
        {
            iC = fReader[i]->getNextPedestalEvent();
        }
        iB = ( iC && iB );
    }
    return iB;
}


bool VMultipleGrIsuReader::getNextShowerEvent()
{
    bool iB = true;
    for( unsigned int i = 0; i < fReader.size(); i++ )
    {
        bool iC = false;
        if( fReader[i] && fSelectedTelescope[i] )
        {
            iC = fReader[i]->getNextShowerEvent();
        }
        iB = ( iC && iB );
    }
    return iB;
}


int VMultipleGrIsuReader::getMC_primary()
{
    if( getReader() )
    {
        return getReader()->getMC_primary();
    }
    
    return 0;
}


float VMultipleGrIsuReader::getMC_energy()
{
    if( getReader() )
    {
        return getReader()->getMC_energy();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_X()
{
    if( getReader() )
    {
        return getReader()->getMC_X();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Y()
{
    if( getReader() )
    {
        return getReader()->getMC_Y();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Xcos()
{
    if( getReader() )
    {
        return getReader()->getMC_Xcos();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Ycos()
{
    if( getReader() )
    {
        return getReader()->getMC_Ycos();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Ze()
{
    if( getReader() )
    {
        return getReader()->getMC_Ze();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Az()
{
    if( getReader() )
    {
        return getReader()->getMC_Az();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Xoffset()
{
    if( getReader() )
    {
        return getReader()->getMC_Xoffset();
    }
    
    return 0.;
}


float VMultipleGrIsuReader::getMC_Yoffset()
{
    if( getReader() )
    {
        return getReader()->getMC_Yoffset();
    }
    
    return 0.;
}
