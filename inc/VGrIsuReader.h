//! VGrIsuReader steering class for GrIsu MC output format

#ifndef VGRISUREADER_H
#define VGRISUREADER_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"

#include "VDetectorGeometry.h"
#include "VMonteCarloRunHeader.h"
#include "VSkyCoordinatesUtilities.h"
#include "VVirtualDataReader.h"

using namespace std;

class VGrIsuReader : public VVirtualDataReader
{
    private:
        bool fDebug;
        bool fPedestalMode;                       //!< if true, create pedestal events, otherwise read shower events from file
        unsigned int  fNPedestalEvents;           //!< total number of pedestal events to be generated
        VDetectorGeometry* fDetectorGeo;          //!< telescope and camera geometry

        bool fMultiGrIsuReader;
        bool fIgnoreCFGVersions;

        double degrad;

        streampos sp;                             //!< current position in data file
        string fDataFormat;

        string fSourceFileName;
        bool fBZipped;                            //!< source file is zipped
        bool fBZipped2;                           //!< source file is zipped
        ifstream is;
        string is_Temp;
        bool fLastWithData;
        uint32_t fEventNumberPreviousEvent;
        uint32_t fEventNumber;
        vector< int > fSumWindow;

        string fExternalPedFile;
        ifstream is_ped;

        unsigned int fNTel;

        unsigned int fTelescopeID;                //!< telescope ID
        vector< uint32_t > fMaxChannels;          //!< number of channels (telescope vector)
        vector< uint16_t > fNumSamples;
        vector< uint32_t > fSelectedHitChan;
        vector< uint16_t > fNumChannelsHit;
        int fTelNumberOffset;                     //!< offset in telescope numbering (from different grisu versions)
        int fSampleOffset;                        //!< offset in FADC sample reading
        double fFADCScale;

        uint8_t  fEventType;                      //!< event type (ignore/observing/pedestal/laser=0/1/2/3)
        vector< std::valarray< double > > fPeds;  //!< pedestal vector  (case of trace library: tube numbering of real data file)
        //!< pedestal variance (case of trace library: tube numbering of real data file)
        vector< std::valarray< double > > fPedvars;
        //!< pedestal variance for all summation windows
        vector< vector< std::valarray< double > > > fVPedvars;
        vector< std::valarray< double > > fPedRMS;//!< pedestal RMS
        uint8_t fdefaultPedUI;
        double fdefaultPed;                       //!< default values, if no pedestals are in the GrIsu input file
        double fdefaultPedvars;                   //!< default values, if no pedestals are in the GrIsu input file
        vector< vector< TH1F* > > fhPeds;
        //!< vector with pedestals from GrIsu file ("P" line) [Tel][Channel][Val]
        vector< vector< vector< uint8_t > > > fGrIsuPeds;
        //!< size of fGrIsuPeds
        vector< vector< unsigned int > > fGrIsuPedsN;

        unsigned int fNoiseTraceStart;
        bool bNoiseTraceFilled;

        // trigger values (vectors [NTelescope])
        int fArrayTrigger;
        vector< bool > fLocalTrigger;             //! bit coded integer, 1 for local trigger
        vector< float > fLocalTriggerTime;
        vector< float > fLocalDelayedTriggerTime;
        double fXimpactrot;                       //!< impact point in shower coordinates
        double fYimpactrot;                       //!< impact point in shower coordinates

        TRandom3* fRandomGen;                     //!< random generator for pedestal generation

        vector< std::vector< bool > > fHiLo;
        vector< std::vector< bool > > fFullHitVec;
        vector< std::vector< bool > > fFullTrigVec;
        vector< int > fNumberofFullTrigger;
        vector< std::vector< int > > fFullAnaVec;
        vector< std::vector< unsigned int> > fPixelConvertVecM;
        vector< std::vector< unsigned int> > fPixelConvertVecR;

        //!< vector with sample vectors of current event
        vector< vector< std::vector< uint8_t > > > fSamplesVec;
        vector< uint8_t > ftempSample;
        vector< int > ftempSampleUI;
        vector< float > ftempSampleFL;

        vector< uint8_t >  vv8;
        vector< vector< uint8_t > > v8;           //!< ignore

        // traces for background
        TFile* fTraceFile;                        //!< file for background trace library
        string fTraceFileName;
        vector< TTree* > fTraceTree;              //!< vector for background trace trees (one per telescope, preliminary)
        short int fTrace[500][64];                //!< background trace array

        //  random dead channels
        int fMCNdead;                             //!< number of pixels set randomly dead
        int fMCNdeadboard;                        //!< number of boards set randomly dead (10 dead pixels in a row)

        // MC values
        int fMC_primary;                          //!< MC primary type
        float  fMC_energy;                        //!< MC primary energy in TeV
        float  fMC_X;                             //!< MC x-coordinate of impact point on ground plane
        float  fMC_Y;                             //!< MC y-coordinate of impact point on ground plane
        float  fMC_Xcos;                          //!< MC x direction cosine of primary in ground coordinate system
        float  fMC_Ycos;                          //!< MC y direction cosine of primary in ground coordinate system
        float  fMC_Ze;                            //!< MC zenith angle of primary in ground coordinate system
        float  fMC_Az;                            //!< MC azimuth angle of primary in ground coordinate system
        float  fMC_Xoffset;                       //!< MC x coordinate of source location in degrees
        float  fMC_Yoffset;                       //!< MC y coordinate of source location in degrees

        // telescope pointing
        vector< double > fTelElevation;           //!< telescope pointing, elevation [deg]
        vector< double > fTelAzimuth;             //!< telescope pointing, azimuth [deg]

        unsigned int fGrisuVersion;

        void fillBackgroundfromTraceLibrary();
        void fillBackgroundfromPlines();
        unsigned int  getGrisuVersion();
        void openDataFile( bool iPeds = false );
        void closeDataFile( bool iPeds = false );
        void initVectors();                       //!< initialize the data vectors
        bool openTraceLibrary();                  //!< open root file for background trace library
        bool readCamera( int iTel );              //!< read trigger and analysis pixel from camera file
        void readPeds( unsigned int i = 0 );
        void readPedsfromLibrary( unsigned int iTel );
        void readPedsfromPlines();
        bool readPixelMix( int iTel );
        void resetEvent();

    public:
        VGrIsuReader( VDetectorGeometry*, unsigned int iNTel, string iExPedFile, vector< int > i_sumwindow, bool iDebug, int iseed, double ifadcscale = 1. );
        VGrIsuReader( VDetectorGeometry*, unsigned int iNTel, string i_sourcefile, vector< int > i_sumwindow, int i_telnumberoffset, int i_sampleoffset, double ifadcscale, bool iDebug, int iseed, string iExPedFile = "", bool iIgnoreCFGFiles = false );
        ~VGrIsuReader();
        string                      getDataFormat()
        {
            return fDataFormat;
        }
        unsigned int                getDataFormatNum()
        {
            return 1;
        }
        string                      getSourceFileName()
        {
            return fSourceFileName;
        }
        // raweventparser
        //!<
        std::pair< bool, uint32_t > getChannelHitIndex( uint32_t );
        uint32_t                    getRunNumber()
        {
            return 0;
        }
        //!< get event number
        uint32_t                    getEventNumber()
        {
            return fEventNumber;
        }
        // don't know the difference between hit and trig
        //!< get hit vector
        std::vector< bool >         getFullHitVec()
        {
            return fFullHitVec[fTelescopeID];
        }
        //!< get triggered channels
        std::vector< bool >         getFullTrigVec()
        {
            return fFullTrigVec[fTelescopeID];
        }
        int                         getNumberofFullTrigger()
        {
            return fNumberofFullTrigger[fTelescopeID];
        }
        //!< get pixels for analysis
        std::vector< int >          getFullAnaVec()
        {
            return fFullAnaVec[fTelescopeID];
        }
        uint8_t                     getEventType()
        {
            return fEventType;
        }
        uint8_t                     getATEventType()
        {
            return fEventType;
        }
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
        //!< return maximum number of channels
        uint16_t                    getMaxChannels()
        {
            return fMaxChannels[fTelescopeID];
        }
        vector< vector< vector< uint8_t > > > getFullNoiseVec()
        {
            return fGrIsuPeds;
        }
        vector< vector< uint8_t > >& getFullNoiseVec( unsigned int iTel );
        vector< uint8_t >&          getFullNoiseVec( unsigned int iTel, unsigned int iChannel );
        uint8_t                     getNoiseSample( unsigned int iTel, uint32_t iHitID, unsigned int iSample, bool iNewTrace = true );
        std::vector< uint8_t >&     getNoiseVec( unsigned int iTel, uint32_t iHitID, bool iNewTrace = true );
        //!< return number of hit channels
        uint16_t                    getNumChannelsHit()
        {
            return fNumChannelsHit[fTelescopeID];
        }
        //!< return number of samples in FADC trace
        uint16_t                    getNumSamples()
        {
            return fNumSamples[fTelescopeID];
        }
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
        std::valarray< double >&    getPeds()
        {
            return fPeds[fTelescopeID];
        }
        std::valarray< double >&    getPedvars()
        {
            return fPedvars[fTelescopeID];
        }
        std::vector< valarray<double> >& getPedvarsAllSumWindows()
        {
            return fVPedvars[fTelescopeID];
        }
        std::valarray< double >&    getPedRMS()
        {
            return fPedRMS[fTelescopeID];
        }
        uint8_t                     getSample( unsigned channel, unsigned sample, bool iNewNoiseTrace = true );
        //!< return FADC samples vector for current hit channel
        std::vector< uint8_t >      getSamplesVec();
        std::vector< double >       getTelElevation()
        {
            return fTelElevation;
        }
        std::vector< double >       getTelAzimuth()
        {
            return fTelAzimuth;
        }
        //!< get selected telescope number
        unsigned int                getTelescopeID()
        {
            return fTelescopeID;
        }
        //!< select hit channel
        void                        selectHitChan( uint32_t );
        void                        setIgnoreCFGVersions( bool iB = false )
        {
            fIgnoreCFGVersions = iB;
        }
        void                        setPedestalEventMode( int i_Nped )
        {
            fPedestalMode = true;
            fNPedestalEvents = i_Nped;
        }
        void                        setRandomDead( int iNC, int iNB );
        //!< set summation window
        void                        setSumWindow( vector< int > i_sum )
        {
            fSumWindow = i_sum;
        }
        void                        setSumWindow( unsigned int iTelID, int isw )
        {
            if( iTelID < fSumWindow.size() )
            {
                fSumWindow[iTelID] = isw;
            }
        }
        //!< set telescope numbering offset
        void                        setTelNumberOffset( int iOff )
        {
            fTelNumberOffset = iOff;
        }
        //!< select telescope
        bool                        setTelescopeID( unsigned int );
        void                        setTraceFile( string i_Trf );
        //!< set trigger values
        void                        setTrigger( vector<bool> iImage, vector<bool> iBorder );
        //!< true if there was a local trigger
        vector< bool >&             getLocalTrigger()
        {
            return fLocalTrigger;
        }
        float                       getLocalTriggerTime( unsigned int iTel )
        {
            return fLocalTriggerTime[iTel];
        }
        float                       getLocalDelayedTriggerTime( unsigned int iTel )
        {
            return fLocalDelayedTriggerTime[iTel];
        }
        //!< return number of telescopes with local trigger
        unsigned int                getNTelLocalTrigger();
        bool                        hasArrayTrigger()
        {
            return !( fArrayTrigger == 0 );
        }
        //!< check if telescope has local trigger
        bool                        hasLocalTrigger( unsigned int iTel );
        double                      getXimpactrot()
        {
            return fXimpactrot;
        }
        double                      getYimpactrot()
        {
            return fYimpactrot;
        }
        void                        setDefaultPed( double iD );

        void                        assignGrisuPeds( unsigned int i = 0 );
        void                        setMultiGrIsuReader( bool iB = true )
        {
            fMultiGrIsuReader = iB;
        }

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
        int                        getMC_primary()//!< MC primary type
        {
            return fMC_primary;
        }
        float                      getMC_energy() //!< MC primary energy
        {
            return fMC_energy;
        }
        float                      getMC_X()      //!< MC x-coordinate of impact point on ground plane
        {
            return fMC_X;
        }
        float                      getMC_Y()      //!< MC y-coordinate of impact point on ground plane
        {
            return fMC_Y;
        }
        float                      getMC_Xcos()   //!< MC x direction cosine of primary in ground coordinate system
        {
            return fMC_Xcos;
        }
        float                      getMC_Ycos()   //!< MC y direction cosine of primary in ground coordinate system
        {
            return fMC_Ycos;
        }
        float                      getMC_Ze()     //!< MC zenith angle of primary in ground coordinate system [rad]
        {
            return fMC_Ze * degrad;
        }
        float                      getMC_Az()     //!< MC azimuth angle of primary in ground coordinate system [rad]
        {
            return fMC_Az * degrad;
        }
        float                      getMC_Xoffset()//!< MC x coordinate of source location in degrees
        {
            return  fMC_Xoffset;
        }
        float                      getMC_Yoffset()//!< MC x coordinate of source location in degrees
        {
            return  fMC_Yoffset;
        }
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
