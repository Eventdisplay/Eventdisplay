//! VPEReader reader for data summary source files

#ifndef VPEReader_H
#define VPEReader_H

#include "TFile.h"
#include "TTree.h"

#include <bitset>
#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>

#include "VDetectorGeometry.h"
#include "VMonteCarloRunHeader.h"
#include "VPETree.h"
#include "VSkyCoordinatesUtilities.h"
#include "VVirtualDataReader.h"

using namespace std;

class VPEReader : public VVirtualDataReader
{
    private:
    
        double degrad;
        
        bool fDebug;
        string fSourceFileName;
        
        VDetectorGeometry* fDetectorGeometry;
        
        TFile* fPE_file;                          //!< pe source file
        unsigned int fPE_treeEvent;               //!< tree event number
        vector< VPETree* > fPE_Tree;              //!< data tree
        bool   fMC;                               //!< source data is Monte Carlo
        
        unsigned int fTelID;
        unsigned int fNTelescopes;
        vector< unsigned int > fTeltoAna;
        unsigned int fSelectedHitChannel;
        vector< unsigned int > fNChannel;
        vector< valarray< double > > fSums;
        vector< vector< valarray< double > > > fTracePulseTiming;
        vector< float > fLTtime;
        vector< float > fLDTtime;
        vector< double > fTelAzimuth;
        vector< double > fTelElevation;
        
        vector< vector< bool > > fFullHitVec;
        vector< vector< bool > > fFullTrigVec;
        vector< vector < bool > > fHiLo;
        vector< int > fNumberofFullTrigger;
        vector< uint8_t > fDummySample;
        
        // tree variables
        unsigned int fPE_runnumber;
        unsigned int fPE_eventnumber;
        unsigned int fPE_eventtype;
        
        // MC parameters
        int fPE_primary;
        float fPE_energy;
        float fPE_xcore;
        float fPE_ycore;
        float fPE_az;
        float fPE_ze;
        float fPE_xcos;
        float fPE_ycos;
        float fPE_Tel_xoff;
        float fPE_Tel_yoff;
        
        bool           fArrayTrigger;
        unsigned int   fNLocalTrigger;
        vector< bool > fLocalTrigger;
        
        bool init();                              //!< open source file and init tree
        
    public:
        VPEReader( string isourcefile, vector< unsigned int > iTeltoAna, VDetectorGeometry* iD, bool iDebug );
        ~VPEReader();
        std::pair< bool, uint32_t > getChannelHitIndex( uint32_t i_channel );
        string    getDataFormat()
        {
            return "PE";
        }
        unsigned int getDataFormatNum()
        {
            return 6;
        }
        uint32_t  getEventNumber()
        {
            return fPE_eventnumber;
        }
        uint8_t   getEventType()
        {
            return fPE_eventtype;
        }
        uint8_t   getATEventType()
        {
            return fPE_eventtype;
        }
        vector< bool >    getFullHitVec()
        {
            return fFullHitVec[fTelID];
        }
        vector< bool >    getFullTrigVec()
        {
            return fFullTrigVec[fTelID];
        }
        int                    getNumberofFullTrigger()
        {
            return fNumberofFullTrigger[fTelID];
        }
        uint32_t  getGPS0()
        {
            return 0;
        }
        uint32_t  getGPS1()
        {
            return 0;
        }
        uint32_t  getGPS2()
        {
            return 0;
        }
        uint32_t  getGPS3()
        {
            return 0;
        }
        uint32_t  getGPS4()
        {
            return 0;
        }
        uint16_t  getGPSYear()
        {
            return 0;
        }
        uint16_t  getATGPSYear()
        {
            return 0;
        }
        uint32_t  getHitID( uint32_t );
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
            return fLocalTrigger;
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
        uint16_t  getMaxChannels()
        {
            return fNChannel[fTelID];
        }
        uint16_t  getNumChannelsHit()             //!< preli
        {
            return fNChannel[fTelID];
        }
        uint16_t  getNumSamples()                 //!< no samples in pe file!
        {
            return 0;
        }
        int       getMC_primary()                 //!< MC primary type
        {
            return fPE_primary;
        }
        float     getMC_energy()                  //!< MC primary energy
        {
            return fPE_energy;
        }
        float     getMC_X()                       //!< MC x-coordinate of impact point on ground plane
        {
            return fPE_xcore;
        }
        float     getMC_Y()                       //!< MC y-coordinate of impact point on ground plane
        {
            return fPE_ycore;
        }
        float     getMC_Xcos()                    //!< MC x direction cosine of primary in ground coordinate system
        {
            return fPE_xcos;
        }
        float     getMC_Ycos()                    //!< MC y direction cosine of primary in ground coordinate system
        {
            return fPE_ycos;
        }
        float     getMC_Ze()                      //!< MC zenith angle of primary
        {
            return fPE_ze;
        }
        float     getMC_Az()                      //!< MC azimuth angle of primary
        {
            return fPE_az;
        }
        float     getMC_Xoffset()                 //!< MC x coordinate of source location in degrees
        {
            return  fPE_Tel_xoff;
        }
        float     getMC_Yoffset()                 //!< MC x coordinate of source location in degrees
        {
            return  fPE_Tel_yoff;
        }
        bool      getNextEvent();
        unsigned int getNumTelescopes()
        {
            return fNTelescopes;
        }
        unsigned int getNTel()
        {
            return fNTelescopes;
        }
        unsigned int getNTelLocalTrigger()
        {
            return fNLocalTrigger;
        }
        uint32_t     getRunNumber()
        {
            return fPE_runnumber;
        }
        vector< uint8_t >   getSamplesVec()
        {
            return fDummySample;
        }
        valarray< double >& getSums( unsigned int iNChannel = 99999 )
        {
            return fSums[fTelID];
        }
        string    getSourceFileName()
        {
            return fSourceFileName;
        }
        vector< double >  getTelAzimuth();
        vector< double >  getTelElevation();
        unsigned int getTelescopeID()
        {
            return fTelID;
        }
        vector< valarray< double > >& getTracePulseTiming( unsigned int iNChannel = 99999 )
        {
            return fTracePulseTiming[fTelID];
        }
        bool      hasArrayTrigger()
        {
            return fArrayTrigger;
        }
        bool      hasFADCTrace()
        {
            return false;
        }
        bool      hasLocalTrigger( unsigned int iTel )
        {
            if( iTel < fLocalTrigger.size() )
            {
                return fLocalTrigger[iTel];
            }
            else
            {
                return false;
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
        bool      setTelescopeID( unsigned int );
        //!< set trigger values
        void      setTrigger( vector<bool> iImage, vector<bool> iBorder );
        bool      wasLossyCompressed()
        {
            return false;
        }
        
        VMonteCarloRunHeader* getMonteCarloHeader()
        {
            return 0;
        }
        void setPerformFADCAnalysis( unsigned int iTel, bool iB )
        {
            iB = false;
        }
};
#endif
