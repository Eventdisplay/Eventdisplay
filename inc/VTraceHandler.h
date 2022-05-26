//! VTraceHandler

#ifndef VTRACE_H
#define VTRACE_H

#include <fstream>
#include <iostream>
#include <stdint.h>
#include <vector>

#include "VVirtualDataReader.h"

using namespace std;

class VTraceHandler
{
    private:
        //!< linear interpolation
        double getLinInterpol( double y5, int x1, double y1, int x2, double y2 );
        
    protected:
        unsigned int      fTraceIntegrationMethod; //   set trace integration method (see setter in source file for definition)
        vector< double >  fpTrace;                //!< the FADC trace
        unsigned int      fpulsetiming_maxPV;
        unsigned int      fpulsetiminglevels_size;
        vector< float >   fpulsetiminglevels;     //!< levels in fraction of maximum for pulse timing calculation
        vector< float >   fpulsetiming;           //!< pulse timing vector (refilled for each call)
        bool              fFindPulseTiming;       //!< true if pulse timing finder was successfull
        unsigned int fpTrazeSize;                 //!< length of the FADC trace
        double fPed;                              //!< Ped value for this trace
        double fTraceAverageTime;                 //!< average time for this trace
        unsigned int fSumWindowFirst;
        unsigned int fSumWindowLast;
        unsigned int fChanID;                     //!< channel ID
        bool fHiLo;                               //!< hilo bit
        int fDynamicRange;                        //!< dynamic range
        int fMaxThreshold;
        unsigned int fMC_FADCTraceStart;          // start of FADC trace (in case the simulated trace is longer than needed)
        bool     kIPRmeasure;                     // if signal extractor is in IPR measurements mode
        
        // digital filter parameters
        unsigned int fDF_method;
        unsigned int fDF_upsample;
        float        fDF_polezero;
        double       fDF_tracemax;
        int          PzpsaSmoothUpsampleFloat( int n, int us, float* ip, float bl, float pz, float* op, float* max, int* at );
        
        bool     apply_digitalFilter();
        bool     apply_lowgain( double );
        double   calculateTraceSum_fixedWindow( unsigned int , unsigned int, bool );
        double   calculateTraceSum_slidingWindow( unsigned int iSearchStart, unsigned int iSearchEnd, int iIntegrationWindow, bool fRaw );
                                      
        void     reset();
        
    public:
        VTraceHandler();
        virtual ~VTraceHandler() {};
        
        virtual void setTrace( vector< uint8_t >, double, unsigned int, double iHiLo = -1. ); //!< pass the trace values (with hilo)
        virtual void setTrace( vector< uint16_t >, double, unsigned int, double iHilo = -1. ); //!< pass the trace values (with hilo)
        virtual void setTrace( VVirtualDataReader* iReader, unsigned int iNSamples, double ped,
                               unsigned int iChanID, unsigned int iHitID, double iHilo = -1. );
        vector< double >& getTrace()
        {
            return fpTrace;
        }
        unsigned int getTraceLength()
        {
            return fpTrazeSize;
        }
        double getTraceAverageTime()
        {
            return fTraceAverageTime;
        }
        unsigned int getTraceIntegrationFirst()
        {
            return fSumWindowFirst;
        }
        unsigned int getTraceIntegrationLast()
        {
            return fSumWindowLast;
        }
        vector<float>  getFADCTiming( unsigned int fFirst, unsigned int fLast, bool debug = false );
        virtual double getTraceSum( unsigned int iSumWindowFirst,
                                    unsigned int iSumWindowLast,
                                    bool iRaw,
                                    unsigned int iTraceIntegrationMethod = 9999,
                                    bool iForceWindowStart = false,
                                    unsigned int iSlidingWindowLast = 9999 );
        virtual vector< float >& getPulseTiming( unsigned int fFirst, unsigned int fLast, 
                                                 unsigned int fTFirst, unsigned int fTLast,
                                                 bool iReverseSearchinLowGain = false );
        virtual bool   getPulseTimingStatus()
        {
            return fFindPulseTiming;
        }
        unsigned int   getTraceIntegrationMethod()
        {
            return fTraceIntegrationMethod;
        }
        virtual double getTraceMax();
        virtual double getTraceMax( unsigned int& n255, unsigned int& maxpos );
        virtual void   getTraceMax( unsigned int fFirst, unsigned int fLast,
                                    double& max, unsigned int& maxpos,
                                    unsigned int& n255, bool iReverseSearchinLowGain = false );
        void    setDigitalFilterParameters( unsigned int iMethod = 0, unsigned int iUpSample = 4, float iPoleZero = 0.75 )
        {
            fDF_method   = iMethod;
            fDF_upsample = iUpSample;
            fDF_polezero = iPoleZero;
        }
        void    setDynamicRange( int iD )
        {
            fDynamicRange = iD;
        }
        void    setMaxThreshold( int iD )
        {
            fMaxThreshold = iD;
        }
        void    setMC_FADCTraceStart( unsigned int iO = 0 )
        {
            fMC_FADCTraceStart = iO;
        }
        void    setIPRmeasure( bool iIPRmeasure = true )
        {
            kIPRmeasure = iIPRmeasure;
        }
        void    setPulseTimingLevels( vector< float > iP );
        bool    setTraceIntegrationmethod( unsigned int iT = 1 );
};
#endif
