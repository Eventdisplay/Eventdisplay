//! VMultipleGrIsuReader  read several grisu files at the same time

#ifndef VMultipleGrIsuReader_H
#define VMultipleGrIsuReader_H

#include "VDetectorGeometry.h"
#include "VGrIsuReader.h"
#include "VMonteCarloRunHeader.h"
#include "VNoiseFileReader.h"

#include "TH1F.h"
#include "TRandom3.h"

#include <iostream>
#include <map>
#include <string>
#include <valarray>
#include <vector>

using namespace std;

class VMultipleGrIsuReader : public VVirtualDataReader
{
    private:
    
        bool fDebug;
        
        unsigned int fNFiles;                     // number of files to be open
        string fSelectedTelescopes;
        
        vector< unsigned int > fTeltoAna;
        vector< bool > fSelectedTelescope;
        vector< string > fSourceFileName;
        vector< VGrIsuReader* > fReader;
        
        VNoiseFileReader* fNoiseFileReader;
        
        VDetectorGeometry* fDetectorGeo;
        
        unsigned int fTelescopeID;
        unsigned int fFileID;
        
        vector< bool > fLocalTrigger;
        
        // telescope pointing
        vector< double > fTelElevation;           //!< telescope pointing, elevation [deg]
        vector< double > fTelAzimuth;             //!< telescope pointing, azimuth [deg]
        
        ///////////////////////////////////////////////////////////////////////////
        // needed for good return values only
        std::vector< bool > vv_bool;
        std::valarray< double > vvv_valarray;
        std::vector< valarray<double> > vvv_v_vvv_valarray;
        ///////////////////////////////////////////////////////////////////////////
        
        bool          checkTelescopeID( unsigned int );
        void          fillRandomPeds( VGrIsuReader* g, int iseed );
        VGrIsuReader* getReader();
        
    public:
    
        VMultipleGrIsuReader( unsigned int nFiles, vector< unsigned int > iteltoana, bool iDebug );
        ~VMultipleGrIsuReader() {}
        
        bool                        init( VDetectorGeometry*, string i_sourcefile, vector< int > i_sumwindow, int i_telnumberoffset, int i_sampleoffset, double ifadcscale, int iseed, string iExPedFile = "", bool iSingleExternalPedFile = true, double iDefaultPed = 20. );
        
        string                      getDataFormat();
        unsigned int                getDataFormatNum()
        {
            return 1;
        }
        string                      getSourceFileName();
        // raweventparser
        //!<
        std::pair< bool, uint32_t > getChannelHitIndex( uint32_t );
        uint32_t                    getRunNumber()
        {
            return 0;
        }
        uint32_t                    getEventNumber();
        // don't know the difference between hit and trig
        std::vector< bool >         getFullHitVec();
        std::vector< bool >         getFullTrigVec();
        int                         getNumberofFullTrigger();
        std::vector< int >          getFullAnaVec();
        uint8_t                     getEventType();
        uint8_t                     getATEventType();
        uint32_t                    getGPS0()     //!< no MC time -> returns 0
        {
            return 0;
        }
        uint32_t                    getGPS1()     //!< no MC time -> returns 0
        {
            return 0;
        }
        uint32_t                    getGPS2()     //!< no MC time -> returns 0
        {
            return 0;
        }
        uint32_t                    getGPS3()     //!< no MC time -> returns 0
        {
            return 0;
        }
        uint32_t                    getGPS4()     //!< no MC time -> returns 0
        {
            return 0;
        }
        uint16_t                    getGPSYear()  //!< no MC time -> returns 0
        {
            return 0;
        }
        uint16_t                    getATGPSYear()//!< no MC time -> returns 0
        {
            return 0;
        }
        //!<
        uint32_t                    getHitID( uint32_t );
        bool                        getHiLo( uint32_t i );
        uint16_t                    getMaxChannels();
        std::vector< uint8_t >      getNoiseVec( unsigned int iTel, uint32_t iHitID );
        uint16_t                    getNumChannelsHit();
        uint16_t                    getNumSamples();
        unsigned int                getNTel()     //!< return number of telescopes
        {
            return fDetectorGeo->getNumTelescopes();
        }
        //!< return number of telescopes
        uint32_t                    getNumTelescopes()
        {
            return fDetectorGeo->getNumTelescopes();
        }
        TH1F*                       getPedHisto( unsigned int itel, unsigned int ichannel );
        std::valarray< double >&    getPeds();
        std::valarray< double >&    getPedvars();
        std::vector< valarray<double> >& getPedvarsAllSumWindows();
        std::valarray< double >&    getPedRMS();
        //!< return FADC samples vector for current hit channel
        std::vector< uint8_t >      getSamplesVec();
        std::vector< double >       getTelElevation();
        std::vector< double >       getTelAzimuth();
        //!< get selected telescope number
        unsigned int                getTelescopeID()
        {
            return fTelescopeID;
        }
        //!< select hit channel
        void                        selectHitChan( uint32_t );
        void                        setPedestalEventMode( int i_Nped );
        void                        setRandomDead( int iNC, int iNB );
        void                        setSumWindow( vector< int > i_sum );
        void                        setSumWindow( unsigned int iTelID, int isw );
        void                        setTelNumberOffset( int iOff );
        //!< select telescope
        bool                        setTelescopeID( unsigned int );
        void                        setTraceFile( string i_Trf );
        //!< set trigger values
        void                        setTrigger( vector<bool> iImage, vector<bool> iBorder );
        vector< bool >&             getLocalTrigger();
        float                       getLocalTriggerTime( unsigned int iTel );
        float                       getLocalDelayedTriggerTime( unsigned int iTel );
        //!< return number of telescopes with local trigger
        unsigned int                getNTelLocalTrigger();
        bool                        hasArrayTrigger();
        //!< check if telescope has local trigger
        bool                        hasLocalTrigger( unsigned int iTel );
        double                      getXimpactrot();
        double                      getYimpactrot();
        void                        setDefaultPed( double iD );
        
        bool                        wasLossyCompressed()
        {
            return false;
        }
        
        // rawfile
        //!< read in next event
        bool                        getNextEvent();
        //!< create next pedestal event
        bool                        getNextPedestalEvent();
        //!< read next shower event from file
        bool                        getNextShowerEvent();
        
        // MC
        bool                       isMC()         //!< GrIsu type data is always MC
        {
            return true;
        }
        bool                       isGrisuMC()
        {
            return true;
        }
        int                        getMC_primary();
        float                      getMC_energy();
        float                      getMC_X();
        float                      getMC_Y();
        float                      getMC_Xcos();
        float                      getMC_Ycos();
        float                      getMC_Ze();
        float                      getMC_Az();
        float                      getMC_Xoffset();
        float                      getMC_Yoffset();
        
        VMonteCarloRunHeader* getMonteCarloHeader()
        {
            return 0;
        }
        void setPerformFADCAnalysis( unsigned int iTel, bool iB )
        {
            ;
        }
};
#endif
