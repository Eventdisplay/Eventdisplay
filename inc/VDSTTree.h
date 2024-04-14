//! VDSTTree   data summarizer class (stores pixel sums and times in a tree)
#ifndef VDSTTree_H
#define VDSTTree_H

#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>

#include "VGlobalRunParameter.h"

////////////////////////////////////////////////////////////////////////////////////////
// MAXIMUM NUMBERS OF TELESCOPES AND CHANNELS ARE DEFINED IN inc/VGlobalRunParameter.h
////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

// basic telescope configuration read from the configuration files
// (usually located in $CTA_EVNDISP_AUX_DIR/DetectorGeometry)
struct VDSTTelescopeConfiguration
{
    float FOV;
    float DynamicRange;
    bool  RAWsum;
    string TelescopeName;
};

class VDSTTree
{
    public:

        TTree* fDST_tree;
        TTree* fDST_conf;
        TTree* fMCtree;
        // [telID] = FOV
        map< unsigned int, VDSTTelescopeConfiguration > fDST_list_of_telescopes;
        vector< unsigned int > fDST_vlist_of_telescopes;

        bool fMC;
        bool fFullTree;
        bool fDST_ADC_set;

        // temporary telescope counter
        int fTelescopeCounter_temp;

        unsigned int fDSTnchannel[VDST_MAXTELESCOPES];

        unsigned int fDSTrunnumber;
        unsigned int fDSTeventnumber;
        unsigned int fDSTeventtype;
        unsigned int fDSTgps0;
        unsigned int fDSTgps1;
        unsigned int fDSTgps2;
        unsigned int fDSTgps3;
        unsigned int fDSTgps4;
        unsigned int fDSTgpsyear;
        unsigned int fDSTATgpsyear;
        // triggered telescopes
        unsigned int fDSTLTrig;
        unsigned int fDSTNTrig;
        unsigned int fDSTLTrig_list[VDST_MAXTELESCOPES];
        // maximum number of telescopes is VDST_MAXTELESCOPES
        // maximum number of channels per camera VDST_MAXCHANNELS
        unsigned int fDSTntel;
        unsigned int fDSTntel_data;
        unsigned int fDSTtel_data[VDST_MAXTELESCOPES];
        float        fDSTpointAzimuth[VDST_MAXTELESCOPES];
        float        fDSTpointElevation[VDST_MAXTELESCOPES];
        int          fDSTpointTrackingKnown[VDST_MAXTELESCOPES];

        unsigned short int fDSTChan[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];

        // data recording parameters
        unsigned short int fDSTRecord[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int fDSTTelescopeZeroSupression[VDST_MAXTELESCOPES];

        // adc parameters
        float        fDSTpedestal[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float        fDSTsums[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];                // integrated charge
        float        fDSTsums2[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];                // integrated charge
        unsigned short int fDSTdead[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int fDSTZeroSuppressed[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int fDSTsumwindow[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int fDSTsumfirst[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float        fDSTt0[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        short        fDSTMax[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int    fDSTPadcHG[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int    fDSTPadcLG[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        // assume that all pulse timing levels are the same for all channels in a telescope
        float        fDSTpulsetiminglevels[VDST_MAXTELESCOPES][VDST_MAXTIMINGLEVELS];
        float        fDSTpulsetiming[VDST_MAXTELESCOPES][VDST_MAXTIMINGLEVELS][VDST_MAXCHANNELS];
        short int    fDSTRawMax[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float        fDSTLTtime[VDST_MAXTELESCOPES];
        float        fDSTLDTtime[VDST_MAXTELESCOPES];
        unsigned short int fDSTL2TrigType[VDST_MAXTELESCOPES];
        unsigned short int fDSTHiLo[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int fDSTN255[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        unsigned short int fDSTnL1trig[VDST_MAXTELESCOPES];
        unsigned short int fDSTL1trig[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        //////////////////////////////////////////////////////////////////////////////////////
        // FADC traces
        bool               fReadWriteFADC;
        unsigned short int fDSTnumSamples[VDST_MAXTELESCOPES];
        unsigned short int fDSTtrace[VDST_MAXTELESCOPES][VDST_MAXSUMWINDOW][VDST_MAXCHANNELS];
        //////////////////////////////////////////////////////////////////////////////////////
        // photoelectrons
        bool  fFillPELeaf;
        unsigned short int fDSTPe[VDST_MAXTELESCOPES][VDST_MAXCHANNELS]; // sum of Che pe in each pixel
        // peak ADC values
        bool  fFillPeakADC;

        // mean pulse timing
        float fDSTMeanPulseTiming[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float fDSTMeanPulseTiming_N[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float fDSTMeanPulseTimingMinLightLevel;
        TH1F* fDSTMeanPulseTimingHistogram[VDST_MAXTELESCOPES];
        float fDSTTraceWidth[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        //////////////////////////////////////////////////////////////////////////////////////
        // MC parameters
        unsigned short int fDSTprimary;
        float fDSTenergy;
        float fDSTxcore;
        float fDSTycore;
        float fDSTaz;
        float fDSTze;
        float fDSTTel_xoff;
        float fDSTTel_yoff;

        //////////////////////////////////////////////////////////////////////////////////////
        VDSTTree();
        ~VDSTTree() {}
        map< unsigned int, VDSTTelescopeConfiguration> getArrayConfig()
        {
            return fDST_list_of_telescopes;
        }
        bool  getFADC()
        {
            return fReadWriteFADC;
        }
        TTree* getDSTTree()
        {
            return fDST_tree;
        }
        TTree* getMCTree()
        {
            return fMCtree;
        }
        bool isMC()
        {
            return fMC;
        }
        bool initDSTTree( bool iFullTree = false, bool iCalibrationTree = false );
        bool initDSTTree( TTree* t, TTree* c );
        bool initMCTree();
        map< unsigned int, VDSTTelescopeConfiguration > readArrayConfig( string );
        map< unsigned int, unsigned int > readTelescopeTypeList( string );
        void resetDataVectors( unsigned int iCH = 0,
                               unsigned int iMaxNTel = VDST_MAXTELESCOPES,
                               unsigned int iMaxPrevNTel = VDST_MAXTELESCOPES,
                               unsigned int iMaxNChannels = VDST_MAXCHANNELS,
                               unsigned int iMaxNTimingLevels = VDST_MAXTIMINGLEVELS,
                               unsigned int iMaxNSamples = VDST_MAXSUMWINDOW,
                               bool iTriggerReset = false, bool iIsCTADST = false );
        void setFADC( bool iFADC = true )
        {
            fReadWriteFADC = iFADC;
        }
        void setFillPeakADFLeaf( bool iPADC = false )
        {
            fFillPeakADC = iPADC;
        }
        void setFillPELeaf( bool iPE = false )
        {
            fFillPELeaf = iPE;
        }
        void setMC( bool iMC = true )
        {
            fMC = iMC;
        }

        // getters for all variables
        uint32_t     getDSTRunNumber()
        {
            return fDSTrunnumber;
        }
        uint32_t     getDSTEventNumber()
        {
            return fDSTeventnumber;
        }
        uint8_t      getDSTEventType()
        {
            return fDSTeventtype;
        }
        uint32_t     getDSTGPS0()
        {
            return fDSTgps0;
        }
        uint32_t     getDSTGPS1()
        {
            return fDSTgps1;
        }
        uint32_t     getDSTGPS2()
        {
            return fDSTgps2;
        }
        uint32_t     getDSTGPS3()
        {
            return fDSTgps3;
        }
        uint32_t     getDSTGPS4()
        {
            return fDSTgps4;
        }
        uint16_t     getDSTGPSYear()
        {
            return fDSTgpsyear;
        }
        uint16_t     getDSTATGPSYear()
        {
            return fDSTATgpsyear;
        }
        unsigned int getDSTNumTelescopes()
        {
            return fDSTntel;
        }
        unsigned int getDSTNTel()
        {
            return fDSTntel;
        }
        unsigned int getDSTNChannels( unsigned int iTel )
        {
            if( iTel < getDSTNTel() )
            {
                return fDSTnchannel[iTel];
            }
            else
            {
                return 0;
            }
        }
        float        getDSTTelAzimuth( unsigned int iTel )
        {
            if( iTel < getDSTNTel() )
            {
                return fDSTpointAzimuth[iTel];
            }
            else
            {
                return 0.;
            }
        }
        float        getDSTTelElevation( unsigned int iTel )
        {
            if( iTel < getDSTNTel() )
            {
                return fDSTpointElevation[iTel];
            }
            else
            {
                return 0.;
            }
        }
        unsigned int getDSTNLocalTrigger()
        {
            return fDSTNTrig;
        }
        bool         getDSTLocalTrigger( int iTelID );
        float        getDSTLocalTriggerTime( int iTelID );
        float        getDSTLocalDelayedTriggerTime( int iTelID );
        unsigned short int getDSTL2TriggerType( int iTelID );

        double       getDSTPedestal( int iChannelID, bool iPrint = false );
        double       getDSTSums( int iChannelID );
        unsigned short int       getDSTPe( int iChannelID );
        unsigned short int       getDSTPe( int iTelID, int iChannelID );
        unsigned short int       getDSTPadcHG( int iChannelID );
        unsigned short int       getDSTPadcHG( int iTelID, int iChannelID );
        unsigned short int       getDSTPadcLG( int iChannelID );
        unsigned short int       getDSTPadcLG( int iTelID, int iChannelID );
        double       getDSTMax( int iChannelID );
        double       getDSTMax( int iTelID, int iChannelID );
        double       getDSTRawMax( int iChannelID );
        double       getDSTRawMax( int iTelID, int iChannelID );
        double       getDSTWidth( int iTelID, int iChannelID );
        double       getDSTTZeros( int iTelID, int iChannelID );
        unsigned int getDSTpulsetiminglevelsN();
        double       getDSTpulsetiming( int iChannelID, int iTimingLevelN );
        double       getDSTpulsetiming( int iTelID, int iChannelID, int iTimingLevelN );
        unsigned int getDSTDead( int iChannelID );
        unsigned int getDSTDead( int iTelID, int iChannelID );
        unsigned short int getZeroSupppressed( int iChannelID );
        unsigned short int getZeroSupppressed( int iTelID, int iChannelID );
        UShort_t     getDSTHiLo( int iChannelID );
        UShort_t     getDSTHiLo( int iTelID, int iChannelID );
        unsigned int getNTrigL1( unsigned int iTelID )
        {
            if( iTelID < getDSTNTel() )
            {
                return fDSTnL1trig[iTelID];
            }
            else
            {
                return 0;
            }
        }
        unsigned int getTrigL1( int iChannelID );
        unsigned int getTrigL1( int iTelID, int iChannelID );

        unsigned short int getDSTNumSample( unsigned int iTelID );
        unsigned short int getDSTTrace( unsigned int iChannelID, unsigned short int iSample );
        unsigned short int getDSTTrace( unsigned int iTelID, unsigned int iChannelID, unsigned short int iSample );

        unsigned short int getDSTMCPrimary()
        {
            return  fDSTprimary;
        }
        float        getDSTMCEnergy()             // [TeV]
        {
            return fDSTenergy;
        }
        // VERITAS coordinate system: x -> east, y-> north
        float        getDSTMCxcore()              // [m]
        {
            return fDSTxcore;
        }
        float        getDSTMCycore()              // [m]
        {
            return fDSTycore;
        }
        float        getDSTMCaz();                // [deg]
        float        getDSTMCze()                 // [deg]
        {
            return fDSTze;
        }
        float        getDSTMCxcos()
        {
            return TMath::Cos( fDSTze * atan( 1. ) / 45. ) * TMath::Cos( fDSTaz * atan( 1. ) / 45. );
        }
        float        getDSTMCycos()
        {
            return TMath::Cos( fDSTze * atan( 1. ) / 45. ) * TMath::Sin( fDSTaz * atan( 1. ) / 45. );
        }
        float        getDSTMCxoff()               // offset of source from camera centre (in camera coordinates) [deg]
        {
            return fDSTTel_xoff;
        }
        float        getDSTMCyoff()               // offset of source from camera centre (in camera coordinates) [deg]
        {
            return fDSTTel_yoff;
        }
        int          getDSTTelescopeNumber( unsigned int iTelHyperArray_ID );

        void         fillDSTMeanPulseTiming( unsigned int iTelID, unsigned int iChannelID, double iTime, int iNSamples = 0 );
        double       getDSTMeanPulseTimingPerTelescope( unsigned int iTelID );
        double       getDSTMedianPulseTimingPerTelescope( unsigned int iTelID );
        double       getDSTRMSPulseTimingPerTelescope( unsigned int iTelID );
        double       getDSTNEventsPulseTimingPerTelescope( unsigned int iTelID );
        double       getDSTMeanPulseTiming( unsigned int iTelID, unsigned int iChannelID );
        double       getDSTMeanPulseTimingMinLightLevel()
        {
            return fDSTMeanPulseTimingMinLightLevel;
        }

        int          hasLocalTrigger( int iTelID );
        int          hasData( int iTelID );

        int         setTelCounter( int iTelID );
};
#endif
