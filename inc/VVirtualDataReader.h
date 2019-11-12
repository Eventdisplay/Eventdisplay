//! VEvndispDataReader  wrapper classe for data reading

#ifndef VVIRTUALDATAREADER_H
#define VVIRTUALDATAREADER_H

#include "VMonteCarloRunHeader.h"
#ifndef NOVBF
#include <VRawDataExceptions.h>
#include <VRawEventData.h>
#else
#include <VDummyArrayTrigger.h>
#endif

#include <iostream>
#include <stdint.h>
#include <string>
#include <valarray>
#include <vector>

#include "TH1F.h"
#include "TTree.h"

using namespace std;

class VVirtualDataReader
{
    private:
        std::valarray<double> v;
        std::vector< unsigned int > vUI;
        std::vector< std::valarray<double> > vv;
        std::vector<int> b;
        std::vector<bool> f;
        std::vector<double> d;
        std::vector< uint16_t > iSampleVec16bit;
        
        string           fSourceFileName;
        vector< unsigned int > fTeltoAna;
        unsigned int fEventStatus;
        uint16_t     fNumSamplesTemp;
        
    public:
        VVirtualDataReader();
        //        VVirtualDataReader() throw (VFileException);
        virtual                            ~VVirtualDataReader() {}
        
        //!< get data format (MC/rawdata/MCvbf/Rawvbf)
        virtual string                      getDataFormat()
        {
            return "virtual";
        }
        //!< data format as integer value (rawdata=0, grisu MC=1,MCvbf=2,Rawvbf=3,DST=4)
        virtual unsigned int                getDataFormatNum()
        {
            return 0;
        }
        virtual vector< int >               getNIncompleteEvents()
        {
            return b;
        }
        virtual unsigned int                getNumTelescopes()
        {
            return 0;
        }
        virtual string                      getSourceFileName()
        {
            return fSourceFileName;
        }
        virtual unsigned int                getTelescopeID()
        {
            return 0;
        }
        //!< set summation window
        virtual void                        setSumWindow( vector< int > i_sum )
        {
        }
        virtual void                        setSumWindow( unsigned int iTelID, int isw ) { }
        //!< read peds
        virtual void                        readPeds( unsigned int i )
        {
        }
        virtual bool                        setTelescopeID( unsigned int a )
        {
            return true;
        }
        // raweventparser
        virtual std::pair< bool, uint32_t > getChannelHitIndex( uint32_t ) = 0;
        virtual uint32_t                    getEventNumber() = 0;
        virtual uint8_t                     getNewEventType()
        {
            return 0;
        }
        virtual uint8_t                     getEventType() = 0;
        unsigned int                        getEventStatus()
        {
            return fEventStatus;
        }
        virtual uint8_t                     getATEventType() = 0;
        virtual uint32_t                    getRunNumber() = 0;
        virtual unsigned int                getDSTTreeEvent()
        {
            return 0;
        }
        virtual bool                        getHitBit( uint32_t )
        {
            return true;
        }
        virtual std::vector< bool >          getFullHitVec() = 0;
        virtual std::vector< bool >          getFullTrigVec() = 0;
        virtual int                         getNumberofFullTrigger() = 0;
        virtual std::vector< int >          getFullAnaVec()
        {
            return b;
        }
        virtual uint32_t                    getGPS0() = 0;
        virtual uint32_t                    getGPS1() = 0;
        virtual uint32_t                    getGPS2() = 0;
        virtual uint32_t                    getGPS3() = 0;
        virtual uint32_t                    getGPS4() = 0;
        virtual uint16_t                    getGPSYear() = 0;
        virtual uint16_t                    getATGPSYear() = 0;
        virtual uint32_t                    getHitID( uint32_t i ) = 0;
        virtual uint16_t                    getMaxChannels() = 0;
        virtual VMonteCarloRunHeader*       getMonteCarloHeader() = 0;
        virtual uint16_t                    getNumChannelsHit() = 0;
        virtual uint16_t                    getNumSamples() = 0;
        virtual unsigned int                getNTel() = 0;
        virtual bool                        getHiLo( uint32_t i ) = 0;
        virtual TH1F*                       getPedHisto( unsigned int a, unsigned int b )
        {
            return 0;
        }
        virtual valarray<double>&           getPeds()
        {
            return v;
        }
        virtual valarray<double>&           getPedvars()
        {
            return v;
        }
        virtual vector< valarray<double> >& getPedvarsAllSumWindows()
        {
            return vv;
        }
        virtual valarray<double>&           getPedRMS()
        {
            return v;
        }
        vector< unsigned int >&             getTeltoAna()
        {
            return fTeltoAna;
        }
        virtual double                      getPedestal( unsigned int channel )
        {
            return 0.;
        }
        virtual uint8_t                     getSample( unsigned channel, unsigned sample, bool iNewNoiseTrace = true )
        {
            return 3;
        }
        virtual std::vector< uint8_t >      getSamplesVec() = 0;
        virtual uint16_t                    getSample16Bit( unsigned channel, unsigned sample, bool iNewNoiseTrace = true )
        {
            return 3;
        }
        double                              getSample_double( unsigned channel, unsigned sample, bool iNewNoiseTrace = true );
        virtual std::vector< uint16_t >     getSamplesVec16Bit()
        {
            return iSampleVec16bit;
        }
        virtual void                        selectHitChan( uint32_t ) = 0;
        void                                setNumSamples( unsigned int iT, uint16_t iS )
        {
            fNumSamplesTemp = iS;
        }
        void                                setTeltoAna( vector< unsigned int > iT )
        {
            fTeltoAna = iT;
        }
        void                                setEventStatus( unsigned int iS )
        {
            fEventStatus = iS;
        }
        virtual bool                        wasLossyCompressed() = 0;
        
        // rawfile
        virtual bool                        getNextEvent() = 0;
        virtual bool                        getPrevEvent()
        {
            return true;
        }
#ifndef NOVBF
        virtual VArrayTrigger*              getArrayTrigger()
        {
            return 0;
        }
#else
        virtual VDummyArrayTrigger*         getArrayTrigger()
        {
            return 0;
        }
#endif
        // telescopes pointing
        virtual std::vector< double >       getTelElevation()
        {
            return d;
        }
        virtual std::vector< double >       getTelAzimuth()
        {
            return d;
        }
        
        // MC
        virtual void                       setTrigger( vector<bool> iImage, vector<bool> iBorder ) {}
        virtual bool                       isMC() //!< is data MC?
        {
            return false;
        }
        //!< is data MC?
        virtual bool                       isGrisuMC()
        {
            return false;
        }
        bool                       isDST()
        {
            if( getDataFormatNum() == 4 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        bool                       isPEFormat()
        {
            if( getDataFormatNum() == 6 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        bool                       isVBFMC()
        {
            if( getDataFormatNum() == 2 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        virtual TTree*                     getMCTree()
        {
            return 0;
        }
        //!< MC primary type
        virtual int                        getMC_primary()
        {
            return 0;
        }
        //!< MC primary energy
        virtual float                      getMC_energy()
        {
            return 0.;
        }
        //!< MC x-coordinate of impact point on ground plane
        virtual float                      getMC_X()
        {
            return 0.;
        }
        //!< MC y-coordinate of impact point on ground plane
        virtual float                      getMC_Y()
        {
            return 0.;
        }
        //!< MC x direction cosine of primary in ground coordinate system
        virtual float                      getMC_Xcos()
        {
            return 0.;
        }
        //!< MC y direction cosine of primary in ground coordinate system
        virtual float                      getMC_Ycos()
        {
            return 0.;
        }
        //!< MC zenith angle of primary in ground coordinate system
        virtual float                      getMC_Ze()
        {
            return 0.;
        }
        //!< MC azimuth angle of primary in ground coordinate system
        virtual float                      getMC_Az()
        {
            return 0.;
        }
        //!< MC x coordinate of source location in degrees
        virtual float                      getMC_Xoffset()
        {
            return  0.;
        }
        //!< MC x coordinate of source location in degrees
        virtual float                      getMC_Yoffset()
        {
            return  0.;
        }
        virtual float getMCFirstInteractionHeight()
        {
            return 0 ;
        }
        virtual float getMCFirstInteractionDepth()
        {
            return 0 ;
        }
        virtual int getMCCorsikaRunID()
        {
            return 0 ;
        }
        virtual int getMCCorsikaShowerID()
        {
            return 0 ;
        }
        
        virtual bool hasFADCTrace()
        {
            return true;
        }
        virtual bool isZeroSuppressed( unsigned int iChannel )
        {
            return false;
        }
        // flag describing zero suppression:
        // 0 = not surpressed
        // bit 0: charge suppressed
        // bit 1: samples suppressed
        virtual unsigned short int getZeroSuppressionFlag( unsigned int iChannel )
        {
            return 0;
        }
        virtual bool has16Bit()
        {
            return false;
        }
        
        //!< has this event an array trigger or not
        virtual bool                       hasArrayTrigger()
        {
            return true;
        }
        //!< return number of telescopes with local trigger
        virtual unsigned int               getNTelLocalTrigger()
        {
            return 0;
        }
        virtual std::vector< bool >&       getLocalTrigger()
        {
            return f;
        }
        virtual bool                       hasLocalTrigger( unsigned int iTel )
        {
            return false;
        }
        virtual float                      getLocalTriggerTime( unsigned int iTel )
        {
            return -999.;
        }
        virtual float                      getLocalDelayedTriggerTime( unsigned int iTel )
        {
            return -999.;
        }
        virtual unsigned short int         getLocalTriggerType( unsigned int iTel )
        {
            return 0;
        }
        
        // DST returns
        virtual std::valarray< double >&     getPedestal()
        {
            return v;
        }
        virtual std::valarray< double >&     getPE( unsigned int iNChannel = 99999 )
        {
            return v;
        }
        virtual std::valarray< double >&     getSums( unsigned int iNChannel = 99999 )
        {
            return v;
        }
        virtual std::valarray< double >&     getTZeros( unsigned int iNChannel = 99999 )
        {
            return v;
        }
        virtual std::valarray< double >&     getTraceMax( unsigned int iNChannel = 99999 )
        {
            return v;
        }
        virtual std::valarray< double >&     getTraceRawMax( unsigned int iNChannel = 99999 )
        {
            return v;
        }
        virtual std::valarray< double >&     getTraceWidth( unsigned int iNChannel = 99999 )
        {
            return v;
        }
        virtual std::vector< valarray< double > >& getTracePulseTiming( unsigned int iNChannel = 99999 )
        {
            return vv;
        }
        virtual std::vector< unsigned int >& getDead()
        {
            return vUI;
        }
        
        virtual void setPerformFADCAnalysis( unsigned int iTel, bool iB )
        {
            iB = false;
        }
        
};
#endif
