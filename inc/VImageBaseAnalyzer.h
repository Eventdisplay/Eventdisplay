//! VImageBaseAnalyzer   basic image analyzer routine (used by VImageAnalyzer, VDST, and VCalibrator)

#ifndef VIMAGEBASEANALYZER_H
#define VIMAGEBASEANALYZER_H

#include <VEvndispData.h>

#include "TTree.h"

#include <iostream>
#include <stdexcept>
#include <valarray>
#include <vector>

using namespace std;

class VImageBaseAnalyzer : public VEvndispData
{
    protected:
        vector<bool> fCalibrated;                 //!  true = calibration is done
        bool fRaw;
        
        void calcSecondTZerosSums();
        void calcTZeros( int , int );
        void calcTZerosSums( int, int, unsigned int );
        unsigned int getDynamicSummationWindow( unsigned int chanID );
        int  getFADCTraceIntegrationPosition( int iPos );
        void FADCStopCorrect();
        bool setSpecialChannels();
        void timingCorrect();
        TTree* makeDeadChannelTree();
        
        void initializeTrace( bool iMakingPeds, unsigned int i_channelHitID,
                              unsigned int i, unsigned int iTraceIntegrationMethod );
                              
    public:
        VImageBaseAnalyzer() {}
        ~VImageBaseAnalyzer() {}
        
        void           calcSums( int iFirst , int iLast, bool iMakingPeds, bool iLowGainOnly = false, unsigned int iTraceIntegrationMethod = 9999 );
        unsigned int   fillHiLo();                          //!< fill hi/low gain vector
        int            fillSaturatedChannels();
        unsigned int   fillZeroSuppressed();
        void           findDeadChans( bool iLowGain = false, bool iFirst = true );
        void           gainCorrect();
};
#endif
