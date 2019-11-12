//! VDSTReader reader for data summary source files

#ifndef VDSTREADER_H
#define VDSTREADER_H

#include <VGlobalRunParameter.h>
#include <VDSTTree.h>
#include <VVirtualDataReader.h>

#include "TFile.h"
#include "TTree.h"

#include <bitset>
#include <iostream>
#include <valarray>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////
// MAXIMUM NUMBERS OF TELESCOPES AND CHANNELS ARE DEFINED IN EVNDISP_definition.h
///////////////////////////////////////////////////////////////////////////////////

class VDSTReader : public VVirtualDataReader
{
    private:
        bool fDebug;
        string fSourceFileName;
        TFile* fDSTfile;                          //!< dst source file
        unsigned int fDSTtreeEvent;               //!< tree event number
        VDSTTree* fDSTTree;                       //!< data tree
        bool   fMC;                               //!< source data is Monte Carlo
        
        vector< bool > fPerformFADCAnalysis;              //!< look at FADC traces
        
        unsigned int fTelID;
        unsigned int fNTelescopes;
        vector< uint16_t >     fNumSamples;
        unsigned int fSelectedHitChannel;
        vector< unsigned int > fNChannel;
        vector< valarray< double > > fPedestal;
        vector< valarray< double > > fSums;
        vector< valarray< double > > fPe;
        vector< vector< valarray< double > > > fTracePulseTiming;
        vector< vector < unsigned int > > fDead;
        vector< valarray< double > > fTraceMax;
        vector< valarray< double > > fRawTraceMax;
        vector< float > fLTtime;
        vector< float > fLDTtime;
        vector< double > fTelAzimuth;
        vector< double > fTelElevation;
        
        vector< vector< bool > > fFullHitVec;
        vector< vector< bool > > fFullTrigVec;
        vector< vector < bool > > fHiLo;
        vector< int > fNumberofFullTrigger;
        
        vector< uint8_t > fDummySample;
        vector< uint16_t > fDummySample16Bit;
        vector< vector< vector< uint16_t > > > fFADCTrace;
        
        vector< bool > fDSTvltrig;
        vector< unsigned short int > fDSTl2trig_type;
        
        bool init();                              //!< open source file and init tree
        
    public:
        VDSTReader( string isourcefile, bool iMC, int iNTel, bool iDebug );
        ~VDSTReader() {}
        bool isZeroSuppressed( unsigned int iChannel );
        unsigned short int getZeroSuppressionFlag( unsigned int iChannel );
        std::pair< bool, uint32_t > getChannelHitIndex( uint32_t i_channel );
        string                      getDataFormat()
        {
            return "DST";
        }
        unsigned int                getDataFormatNum()
        {
            return 4;
        }
        vector<unsigned int>&       getDead()
        {
            return fDead[fTelID];
        }
        uint32_t                    getEventNumber()
        {
            return fDSTTree->getDSTEventNumber();
        }
        uint8_t                     getEventType()
        {
            return fDSTTree->getDSTEventType();
        }
        uint8_t                     getATEventType()
        {
            return fDSTTree->getDSTEventType();
        }
        unsigned int   getDSTTreeEvent()
        {
            return fDSTtreeEvent;
        }
        vector< bool >              getFullHitVec()
        {
            return fFullHitVec[fTelID];
        }
        vector< bool >              getFullTrigVec()
        {
            return fFullTrigVec[fTelID];
        }
        int                         getNumberofFullTrigger()
        {
            return fNumberofFullTrigger[fTelID];
        }
        uint32_t                    getGPS0()
        {
            return fDSTTree->getDSTGPS0();
        }
        uint32_t                    getGPS1()
        {
            return fDSTTree->getDSTGPS1();
        }
        uint32_t                    getGPS2()
        {
            return fDSTTree->getDSTGPS2();
        }
        uint32_t                    getGPS3()
        {
            return fDSTTree->getDSTGPS3();
        }
        uint32_t                    getGPS4()
        {
            return fDSTTree->getDSTGPS4();
        }
        uint16_t                    getGPSYear()
        {
            return fDSTTree->getDSTGPSYear();
        }
        uint16_t                    getATGPSYear()
        {
            return fDSTTree->getDSTATGPSYear();
        }
        uint32_t                    getHitID( uint32_t );
        bool                        getHiLo( uint32_t i )
        {
            if( i < fHiLo[fTelID].size() )
            {
                return fHiLo[fTelID][i];
            }
            else
            {
                return 0;
            }
        }
        vector< bool >&             getLocalTrigger()
        {
            return fDSTvltrig;
        }
        unsigned short int          getLocalTriggerType( unsigned int iTelID )
        {
            if( iTelID < fDSTl2trig_type.size() )
            {
                return fDSTl2trig_type[iTelID];
            }
            else
            {
                return 0;
            }
        }
        float                       getLocalTriggerTime( unsigned int iTel )
        {
            if( iTel < fLTtime.size() )
            {
                return fLTtime[iTel];
            }
            else
            {
                return -999.;
            }
        }
        float                       getLocalDelayedTriggerTime( unsigned int iTel )
        {
            if( iTel < fLDTtime.size() )
            {
                return fLDTtime[iTel];
            }
            else
            {
                return -999;
            }
        }
        uint16_t                    getMaxChannels()
        {
            return fNChannel[fTelID];
        }
        uint16_t                    getNumChannelsHit()
        {
            return fNChannel[fTelID];
        }
        uint16_t                    getNumSamples()
        {
            if( fTelID < fNumSamples.size() )
            {
                return fNumSamples[fTelID];
            }
            else
            {
                return 0;
            }
        }
        TTree*                      getMCTree();
        int                         getMC_primary()
        {
            return fDSTTree->getDSTMCPrimary();    //!< MC primary type
        }
        float                       getMC_energy()                  //!< MC primary energy
        {
            return fDSTTree->getDSTMCEnergy();
        }
        float     getMC_X()                       //!< MC x-coordinate of impact point on ground plane
        {
            return fDSTTree->getDSTMCxcore();
        }
        float     getMC_Y()                       //!< MC y-coordinate of impact point on ground plane
        {
            return fDSTTree->getDSTMCycore();
        }
        float     getMC_Xcos()                    //!< MC x direction cosine of primary in ground coordinate system
        {
            return fDSTTree->getDSTMCxcos();
        }
        float     getMC_Ycos()                    //!< MC y direction cosine of primary in ground coordinate system
        {
            return fDSTTree->getDSTMCycos();
        }
        float     getMC_Ze()                      //!< MC zenith angle of primary
        {
            return fDSTTree->getDSTMCze();
        }
        float     getMC_Az()                      //!< MC azimuth angle of primary
        {
            return fDSTTree->getDSTMCaz();
        }
        float     getMC_Xoffset()                 //!< MC x coordinate of source location in degrees
        {
            return fDSTTree->getDSTMCxoff();
        }
        float     getMC_Yoffset()                 //!< MC x coordinate of source location in degrees
        {
            return fDSTTree->getDSTMCyoff();
        }
        bool         getNextEvent();
        VMonteCarloRunHeader*         getMonteCarloHeader();
        unsigned int                  getNumTelescopes()
        {
            return fNTelescopes;
        }
        unsigned int                  getNTelLocalTrigger()
        {
            return fDSTTree->getDSTNLocalTrigger();
        }
        unsigned int                  getNTel()
        {
            return fNTelescopes;
        }
        uint32_t                      getRunNumber()
        {
            return fDSTTree->getDSTRunNumber();
        }
        valarray< double >&           getPedestal()
        {
            return fPedestal[fTelID];
        }
        vector< uint8_t >             getSamplesVec();
        uint8_t                       getSample( unsigned channel, unsigned sample, bool iNewNoiseTrace = true );
        vector< uint16_t >            getSamplesVec16Bit();
        uint16_t                      getSample16Bit( unsigned channel, unsigned sample, bool iNewNoiseTrace = true );
        valarray< double >&           getSums( unsigned int iNChannel = 99999 )
        {
            return fSums[fTelID];
        }
        valarray< double >&           getPE( unsigned int iNChannel = 99999 )
        {
            return fPe[fTelID];
        }
        string                        getSourceFileName()
        {
            return fSourceFileName;
        }
        vector< double >              getTelAzimuth()
        {
            return fTelAzimuth;
        }
        vector< double >              getTelElevation()
        {
            return fTelElevation;
        }
        unsigned int                  getTelescopeID()
        {
            return fTelID;
        }
        valarray< double >&           getTraceMax( unsigned int iNChannel = 99999 )
        {
            return fTraceMax[fTelID];
        }
        valarray< double >&           getTraceRawMax( unsigned int iNChannel = 99999 )
        {
            return fRawTraceMax[fTelID];
        }
        vector< valarray< double > >& getTracePulseTiming( unsigned int iNChannel = 99999 )
        {
            return fTracePulseTiming[fTelID];
        }
        bool      has16Bit()
        {
            return true;
        }
        bool      hasFADCTrace()
        {
            if( fDSTTree )
            {
                return fDSTTree->getFADC();
            }
            else
            {
                return false;
            }
        }
        bool      hasLocalTrigger( unsigned int iTel )
        {
            if( fDSTTree->hasLocalTrigger( iTel ) < 0 )
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        bool      isMC()
        {
            return fMC;
        }
        void      selectHitChan( uint32_t hit )
        {
            fSelectedHitChannel = hit;
        }
        void      setNumSamples( unsigned int iTelID, uint16_t iS )
        {
            if( iTelID < fNumSamples.size() )
            {
                fNumSamples[iTelID] = iS;
            }
        }
        void      setPerformFADCAnalysis( unsigned int iTelID, bool iB )
        {
            if( iTelID < fPerformFADCAnalysis.size() )
            {
                fPerformFADCAnalysis[iTelID] = iB;
            }
        }
        bool      setTelescopeID( unsigned int );
        void      setTrigger( vector<bool> iImage, vector<bool> iBorder );          //!< set trigger values
        bool      wasLossyCompressed()
        {
            return false;
        }
};
#endif
