/*! \class VBFDataReader

    steering class for reading vbf data format


*/

#include <VBFDataReader.h>

VBFDataReader::VBFDataReader( string sourcefile, int isourcetype, unsigned int iNTel, bool iDebug, unsigned int iPrintDetectorConfig ):
    VBaseRawDataReader( sourcefile, isourcetype, iNTel, iDebug ),
    pack( NULL ),
    reader( sourcefile, false, true ),
    index( 0 )
{
    at = 0;
    ae = 0;
    fArrayTrigger = false;
    fNIncompleteEvent.assign( iNTel, 0 );
    setDebug( iDebug );
    fPrintDetectorConfig = iPrintDetectorConfig;
}


VBFDataReader::~VBFDataReader()
{
    if( fDebug )
    {
        cout << "VBFDataReader::~VBFDataReader()" << endl;
    }
    if( pack != NULL )
    {
        delete pack;
    }
}


bool VBFDataReader::getNextEvent()
{
    if( fDebug )
    {
        cout << "VBFDataReader::getNextEvent() " << endl;
    }
    bool bSimulations = false;
    
    try
    {
        if( fDebug )
        {
            cout << "bool VBFRawDataReader::getNextEvent() " << index << endl;
        }
        for( ;; )
        {
            if( !reader.hasPacket( index ) )
            {
                setEventStatus( 999 );
                return false;
            }
            VPacket* old_pack = pack;
            try
            {
                pack = reader.readPacket( index );
            }
            catch( const std::exception& e )
            {
                std::cout << "VBFDataReader::getNextEvent: exception while reading file: "
                          << e.what() << std::endl;
                pack = old_pack;
                setEventStatus( 0 );
                return false;
            }
            delete old_pack;
            if( fDebug )
            {
                cout << "\t VBFRawDataReader::getNextEvent(): index " << index << endl;
            }
            fEventNumber = index - 1;
            index++;
            
            // check if this is a simulation header
            if( pack->hasSimulationHeader() )
            {
                printSimulationHeader( pack, fPrintDetectorConfig );
                fMonteCarloHeader = fillSimulationHeader( pack );
            }
            
            if( setSimulationData( pack ) )
            {
                fDataFormat = "MCvbf";
                fTelElevation.clear();
                fTelAzimuth.clear();
                for( unsigned int i = 0; i < fNTel; i++ )
                {
                    fTelElevation.push_back( getSMC_TelPointing_Elevation() );
                    fTelAzimuth.push_back( getSMC_TelPointing_Azimuth() );
                }
                // get eventnumber from simulation package
                fEventNumber = getSMC_eventNumber();
                bSimulations = true;
            }
            
            bool gotOneEv = false;
            fArrayTrigger = false;
            
            if( pack->hasArrayEvent() )
            {
                ae = pack->getArrayEvent();
                if( fDebug )
                {
                    cout << "\t VBFRawDataReader::getNextEvent(): hasArrayEvent ";
                    cout << ", has event number: " << ae->hasEventNumber() << " ";
                    if( ae->hasEventNumber() )
                    {
                        cout << ", array package eventnumber " << ae->getEventNumber();
                    }
                    cout << ", counting event number " << fEventNumber;
                    cout << endl;
                }
                
                setEventStatus( 0 );
                bitset< 8 * sizeof( ULong64_t ) > ib;
                
                // check if event for requested telescopes are available
                bool bComplete = true;
                for( unsigned int i = 0; i < getTeltoAna().size(); i++ )
                {
                    if( getTeltoAna()[i] < ae->getPresentTelescopes().size() )
                    {
                        if( !ae->getPresentTelescopes()[getTeltoAna()[i]] )
                        {
                            fNIncompleteEvent[getTeltoAna()[i]]++;
                            bComplete = false;
                            if( fNIncompleteEvent[getTeltoAna()[i]] < 25 && !bSimulations )
                            {
                                cout << "VBFDataReader::getNextEvent(): missing telescope event for telescope " << getTeltoAna()[i] + 1 << endl;
                                cout << "\t ----  skipping event " << ae->getEventNumber() << " ----" << endl;
                            }
                            else if( fNIncompleteEvent[getTeltoAna()[i]] == 25 && !bSimulations )
                            {
                                cout << "VBFDataReader::getNextEvent(): missing telescope event for telescope " << getTeltoAna()[i] + 1 << endl;
                                cout << "\t ---------------------------------------------------------------" << endl;
                                cout << "\t ---   more than 25 missing telescope events for telescope " << getTeltoAna()[i] + 1 << " ----" << endl;
                                cout << "\t ---------------------------------------------------------------" << endl;
                            }
                        }
                        else
                        {
                            if( getTeltoAna()[i] < ib.size() )
                            {
                                ib.set( getTeltoAna()[i], 1 );
                            }
                        }
                    }
                    setEventStatus( ( unsigned int )ib.to_ulong() );
                }
                if( fDebug )
                {
                    cout << "\t VBFRawDataReader::getNextEvent(): hasArrayEvent, complete? " << ib.to_ulong() << " " << bComplete << endl;
                }
                if( !bComplete )
                {
                    return false;
                }
                if( ae )
                {
                    for( unsigned int i = 0; i < ae->getNumEvents(); ++i )
                    {
                        VEvent* ev = ae->getEvent( i );
                        fEvent[ev->getNodeNumber()] = ev;
                        gotOneEv = true;
                    }
                }
                else
                {
                    for( unsigned int i = 0; i < fEvent.size(); i++ )
                    {
                        fEvent[i] = 0;
                    }
                    return false;
                }
                
                if( ae->hasTrigger() )
                {
                    fArrayTrigger = true;
                    if( fDebug )
                    {
                        cout << "\t VBFRawDataReader::getNextEvent(): hasTrigger" << endl;
                    }
                    at = ae->getTrigger();
                }
                else
                {
                    fArrayTrigger = false;
                    at = NULL;
                }
            }
            else
            {
                for( unsigned int i = 0; i < fEvent.size(); i++ )
                {
                    fEvent[i] = 0;
                }
                ae = 0;
                at = 0;
            }
            if( fDebug )
            {
                cout << "\t VBFRawDataReader::getNextEvent(): hasArrayEvent, event found? " << gotOneEv << " ";
                if( ae )
                {
                    cout << ae->getNumEvents();
                }
                cout << "\t" << gotOneEv;
                cout << endl;
            }
            
            if( !gotOneEv && !bSimulations )
            {
                continue;
            }
            
            if( fDebug )
            {
                cout << "bool VBFRawDataReader::getNextEvent() (end)" << endl;
            }
            
            return true;
        }
    }
    catch( const std::exception& e )
    {
        std::cout << "unexpected exception: " << e.what() << std::endl;
        setEventStatus( 999 );
        return false;
    }
}


unsigned int VBFDataReader::getNTel()
{
    unsigned int z = 0;
    unsigned int t = reader.getConfigMask().size();
    for( unsigned int i = 0; i < t; i++ )
    {
        if( reader.getConfigMask()[i] )
        {
            z++;
        }
    }
    return z;
}


vector< bool >& VBFDataReader::getLocalTrigger()
{
    if( getArrayEvent() )
    {
        ib_temp = getArrayEvent()->getExpectedTelescopes();
    }
    else
    {
        ib_temp.clear();
    }
    
    /* -------
        additonal test to check local trigger vector with eventtype
        but: has to be used differently for pedestal runs, etc...
             needs some more work
    // test if vector is nonzero
       bool bZero = true;
       for( unsigned int i = 0; i < getNTel(); i++ )
       {
          if( ib_temp[i] )
          {
             bZero = false;
    }
    }
    if( bZero )
    {
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
    if( getNewEventType( i ) == 1 ) ib_temp[i] = true;
    }
    }
    ---- */
    
    return ib_temp;
}


bool VBFDataReader::hasLocalTrigger( unsigned int iTel )
{
    if( fDebug )
    {
        cout << "VBFDataReader::hasLocalTrigger" << endl;
    }
    if( iTel < getLocalTrigger().size() )
    {
        return getLocalTrigger()[iTel];
    }
    return false;
}


unsigned int VBFDataReader::getNTelLocalTrigger()
{
    unsigned int iNTrig = 0;
    for( unsigned int i = 0; i < getLocalTrigger().size(); i++ )
    {
        if( getLocalTrigger()[i] )
        {
            iNTrig++;
        }
    }
    return iNTrig;
}


uint16_t VBFDataReader::getATGPSYear()
{
    if( at == NULL )
    {
        return getGPSYear();
    }
    else
    {
        // horrible fudge for a period where GPS is off by a year (VERITAS)
        if( at->getRunNumber() > 33242 && at->getRunNumber() < 33254 )
        {
            return at->getGPSYear() + 1;
        }
        return at->getGPSYear();
    }
}

uint16_t VBFDataReader::getNumSamples()
{
    if( fTelID < fEvent.size() && fEvent[fTelID] )
    {
        return fEvent[fTelID]->getNumSamples();
    }
    
    return 0;     // note: VEventLoop::initEventLoop depends on this
}

bool VBFDataReader::hasArrayTrigger()
{
    return fArrayTrigger;
}
