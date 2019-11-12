/*! \class VRawDataReader

    steering class for reading raw data format



*/

#include <VRawDataReader.h>

VRawDataReader::VRawDataReader( string sourcefile,
                                int isourcetype,
                                unsigned int iNTel,
                                bool iDebug ):
    VBaseRawDataReader( sourcefile, isourcetype, iNTel, iDebug )
{
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fFileRead.push_back( new VRawDataFileRead( sourcefile, i ) );
        fEvent.push_back( new VRawEventData() );
    }
}


VRawDataReader::~VRawDataReader()
{
    if( fDebug )
    {
        cout << "VRawDataReader::~VRawDataReader()" << endl;
    }
    /*   delete fFileRead;
       delete fEvent; */
}


bool VRawDataReader::getNextEvent()
{
    if( fDebug )
    {
        cout << "bool VRawDataReader::getNextEvent()" << endl;
    }
    bool iStatus = true;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        try
        {
            iStatus = fFileRead[i]->getNextEvent( *fEvent[i] );
            if( fDebug )
            {
                cout << "bool VRawDataReader::getNextEvent() status: " << iStatus << endl;
            }
        }
        catch( ... )
        {
            std::cout << "VRawDataReader::getNextEvent: exception while reading file" << std::endl;
            return false;
        }
        
        if( !iStatus )
        {
            return false;
        }
        else
        {
            //	 fEventParser[i]->setEvent(fEvent[i]);
            if( fDataFormat == "MCvbf" )
            {
                // ugly, but same data for all telescopes
                iStatus = setSimulationData( fEvent[i] );
                if( !iStatus )
                {
                    std::cout << "VRawDataReader::getNextEvent: no simulation bank found" << std::endl;
                }
            }
        }
    }
    if( fDebug )
    {
        cout << "bool VRawDataReader::getNextEvent() (end)" << endl;
    }
    
    return true;
}
