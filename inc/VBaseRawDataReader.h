//! VBaseRawDataReader contains common elements of VRawDataReader and VBFDataReader

#ifndef VBASEDATAREADER_H
#define VBASEDATAREADER_H

#include "VDetectorGeometry.h"
#include "VMonteCarloRunHeader.h"
#include "VNoiseFileReader.h"
#include "VRawEventParser.h"
#include "VRawDataExceptions.h"
#include "VVirtualDataReader.h"
#include "VSimulationDataReader.h"

#include "TRandom3.h"

#include <iostream>
#include <stdint.h>
#include <string>
#include <valarray>
#include <vector>

using namespace std;

class VBaseRawDataReader : public VVirtualDataReader, public VSimulationDataReader
{
    protected:
        bool fDebug;
        vector< VEvent* > fEvent;
        uint32_t        fEventNumber;
        unsigned int    fDataFormatNum;
        string          fDataFormat;
        string          fSourceFileName;
        unsigned int    fNTel;
        unsigned int    fTelID;
        uint32_t        fHitID;

        vector< double > fTelElevation;
        vector< double > fTelAzimuth;

        vector< int > fNIncompleteEvent;

        vector< bool > fDummyBoolV;
        vector< uint8_t > fDummyUint8V;
        std::pair< bool, uint32_t > fDummyPair;

        VNoiseFileReader* fNoiseFileReader;
        uint8_t           fNoiseFilePedestal;
        uint8_t           fNoiseFileFADCRange;

        double            finjectGaussianNoise;
        TRandom3*         fRandomInjectGaussianNoise;

        vector< uint16_t > fDefaultMaxNChannels;

        // trace amplitude correction
        vector< float > fTraceAmplitudeCorrectionS;
        vector< float > fTraceAmplitudeCorrectionG;

        VMonteCarloRunHeader* fMonteCarloHeader;

        // QADC values
        std::valarray<double> fSums;
        std::valarray<double> fTraceMax;
        std::vector< valarray< double > > fTracePulseTiming;

        // placeholders
        std::valarray<double> v;
        std::vector< std::valarray<double> > vv;

    public:
        VBaseRawDataReader( string,
                            int isourcetype,
                            unsigned int iNTel,
                            bool iDebug );
        virtual ~VBaseRawDataReader();

        string                      getDataFormat()
        {
            return fDataFormat;
        }
        unsigned int                getDataFormatNum();
        string                      getSourceFileName()
        {
            return fSourceFileName;
        }

        // raweventparser
        virtual bool                hasAT()
        {
            return false;
        }
        std::pair< bool, uint32_t > getChannelHitIndex( uint32_t i )
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getChannelHitIndex( i );
            }
            else
            {
                return fDummyPair;
            }
        }
        uint32_t                    getRunNumber()
        {
            return 0;
        }
        uint32_t                    getEventNumber()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getEventNumber();
            }
            else
            {
                return fEventNumber;
            }
        }
        uint8_t                     getNewEventType();
        uint8_t                     getNewEventType( unsigned int itel );
        uint8_t                     getEventType();
        uint8_t                     getATEventType()
        {
            return getEventType();
        }
        bool                        getHitBit( uint32_t i )
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getHitBit( i );
            }
            else
            {
                return true;
            }
        }
        std::vector< bool >         getFullHitVec()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getFullHitVec();
            }
            else
            {
                return fDummyBoolV;
            }
        }
        std::vector< bool >         getFullTrigVec()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getFullTrigVec();
            }
            else
            {
                return fDummyBoolV;
            }
        }
        int                         getNumberofFullTrigger();
        uint32_t                    getGPS0()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getGPSTime()[0];
            }
            else
            {
                return 0;
            }
        }
        uint32_t                    getGPS1()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getGPSTime()[1];
            }
            else
            {
                return 0;
            }
        }
        uint32_t                    getGPS2()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getGPSTime()[2];
            }
            else
            {
                return 0;
            }
        }
        uint32_t                    getGPS3()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getGPSTime()[3];
            }
            else
            {
                return 0;
            }
        }
        uint32_t                    getGPS4()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getGPSTime()[4];
            }
            else
            {
                return 0;
            }
        }
        uint16_t                    getGPSYear();
        virtual uint16_t            getATGPSYear()
        {
            return getGPSYear();
        }
        uint16_t                    getMaxChannels();
        vector< int >               getNIncompleteEvents()
        {
            return fNIncompleteEvent;
        }
        uint16_t                    getNumChannelsHit()
        {
            if( fEvent[fTelID] )
            {
                return fEvent[fTelID]->getNumChannelsHit();
            }
            else
            {
                return 0;
            }
        }
        uint16_t                    getNumSamples();
        unsigned int                getNTel()
        {
            return 0;
        }
        uint8_t                     getSample( unsigned channel, unsigned sample, bool iNewNoiseTrace = true );
        std::vector< uint8_t >      getSamplesVec();
        uint32_t                    getHitID( uint32_t i );
        bool                        getHiLo( uint32_t i );
        unsigned int                getTelescopeID()
        {
            return fTelID;
        }
        bool                        hasFADCTrace();
        void                        selectHitChan( uint32_t i );
        bool isZeroSuppressed( unsigned int iChannel );
        unsigned short int getZeroSuppressionFlag( unsigned int iChannel );
        bool                        wasLossyCompressed();

        // rawfile
        virtual bool                getNextEvent() = 0;
        bool                        setTelescopeID( unsigned int iTelID );

        bool                       isMC();
        int                        getMC_primary()
        {
            return getSMC_primary();
        }
        float                      getMC_energy()
        {
            return getSMC_energy();
        }
        float                      getMC_X()
        {
            return getSMC_X();
        }
        float                      getMC_Y()
        {
            return getSMC_Y();
        }
        float                      getMC_Xcos()
        {
            return getSMC_Xcos();
        }
        float                      getMC_Ycos()
        {
            return getSMC_Ycos();
        }
        float                      getMC_Ze()
        {
            return getSMC_Ze();
        }
        float                      getMC_Az()
        {
            return getSMC_Az();
        }
        float getMCFirstInteractionHeight()
        {
            return getSMCFirstInteractionHeight() ;
        }
        float getMCFirstInteractionDepth()
        {
            return getSMCFirstInteractionDepth();
        }
        int getMCCorsikaRunID()
        {
            return getSMCCorsikaRunID() ;
        }
        int getMCCorsikaShowerID()
        {
            return getSMCCorsikaShowerID();
        }

        // sky to camera coordinates
        float                      getMC_Xoffset()
        {
            return -1. * getSMC_Xoffset();
        }
        float                      getMC_Yoffset()
        {
            return -1. * getSMC_Yoffset();
        }

        // return QADC values
        valarray< double >&        getSums( unsigned int iNChannel = 99999 );
        valarray< double >&        getTraceMax( unsigned int iNChannel = 99999 );
        valarray< double >&        getTraceRawMax( unsigned int iNChannel = 99999 )
        {
            return getTraceMax( iNChannel );
        }
        vector< valarray< double > >& getTracePulseTiming( unsigned int iNChannel = 99999 );

        // noise file related stuff
        valarray<double>&          getPeds();
        valarray<double>&          getPedvars();
        vector< valarray<double> >& getPedvarsAllSumWindows();
        valarray<double>&          getPedRMS();
        bool                       initTraceNoiseGenerator( unsigned int, string, VDetectorGeometry*, vector<int>, bool, int, double, vector<double> );
        void                       injectGaussianNoise( double injectGaussianNoise, UInt_t seed = 0 );
        bool                       initThroughputCorrection( double , vector< float >, vector< float > );
        void                       setSumWindow( unsigned int iTelID, int isw );

        vector< double >           getTelElevation()
        {
            return fTelElevation;
        }
        vector< double >           getTelAzimuth()
        {
            return fTelAzimuth;
        }

        void                       setDebug( bool iDebug = false );

        VMonteCarloRunHeader* getMonteCarloHeader()
        {
            return fMonteCarloHeader;
        }

        void setDefaultMaxNChannels( vector< uint16_t > i_nChannels )
        {
            fDefaultMaxNChannels = i_nChannels;
        }
};
#endif
