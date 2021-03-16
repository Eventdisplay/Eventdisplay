//! VSpecialChannel read and administrate list of special channels (e.g. L2 channels)

#ifndef VSpecialChannel_H
#define VSpecialChannel_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class VSpecialChannel
{
    private:
    
        bool  fDebug;
        bool  fIsZombie;               // state of this class: was reading of data file successfull?
        
        unsigned int fTelID;           // telescope identifier
        string fSpecialChannelFile;    // file name
        
        vector< unsigned int >      fFADCstopTrigChannelID;      // L2 channels fed into the FADCs crates
        map< unsigned int, double > fHIGHQE_gainfactor;          // relative gain of HIGHQE channels to ordinary channels
        map< unsigned int, unsigned int > fChannelStatus;              // channel status (0 for switched off)
        
        
    public:
    
        VSpecialChannel() {}
        VSpecialChannel( unsigned int iTelID );      // telescope ID: T1=0
        ~VSpecialChannel() {}
        
        unsigned int            getTelID()
        {
            return fTelID;
        }
        bool                    isZombie()
        {
            return fIsZombie;
        }
        unsigned int            getChannelStatus( unsigned int );
        vector< unsigned int >& getFADCstopTrigChannelID()
        {
            return fFADCstopTrigChannelID;
        }
        double                  getHIGHQE_gainfactor( unsigned int );
        bool                    readSpecialChannels( int iRun, string iFile, string iDirectory );
        bool                    readThroughput( string iEpoch, string iFile, string iDirectory, unsigned int iNChannel );
        void                    reset();
        void                    setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        void                    setTelID( unsigned int iTelID )
        {
            fTelID = iTelID;
        }
};

#endif
