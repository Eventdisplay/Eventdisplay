//! VDeadTime dead time calculator

#ifndef VDEADTIME_H
#define VDEADTIME_H

#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TNamed.h"

#include <iostream>
#include <vector>

using namespace std;

class VDeadTime : public TNamed
{
    private:
    
        bool bIsOn;
        
        double fDeadTimeMiss;
        double fDeadTimeMS;
        double fDeadTimeFrac;
        
        double fTFitMin;
        double fTFitMax;
        
        double ft0;
        double fRunStart;
        TH1D* hTimeDiff;
        TH1D* hTimeDiffLog;
        TH2D* hTimeDiff2D;
        TF1*  hFTimeDiff;
        TGraphErrors* hgDeadTime;
        TH1D* hNEventTime;
        TList* hisList;
        
        // scalar related histograms
        TH1D* hScalarClock;
        TH1D* hScalarBusy;
        TH1D* hScalarDeadTimeFraction;
        bool         fFirstScalarEvent;
        unsigned int fScalarClockLastEvent;
        unsigned int fScalarBusyLastEvent;
        double       fScalarTimeOfFirstEvent;
        double       fScalarDeadTimeFrac;
        double       fScalarDeadTimeChi2;
        
        unsigned int fDeadTimeFrac_status;
        double       fDeadTimeConsistencyCheck_allowedDifference;      // allowed difference between the two methods
        
        double  calculateDeadTimeFromScalars();
        double  calculateDeadTimeFromTimeDifferences();
        
    public:
    
        VDeadTime( bool iIsOn = true );
        ~VDeadTime() {}
        double calculateDeadTime();
        bool   checkStatus();
        void   defineHistograms( float iRunDuration = 0., bool iNoWarning = false );
        double fillTimeDifferenceHistograms( double time );
        double fillDeadTime( double time, unsigned int* tenMHzClock = 0 );
        void   fillTenMHzClockArray( double time, unsigned int* tenMHzClock );
        double getScalarDeadTimeFraction()
        {
            return fScalarDeadTimeFrac ;
        }
        double getDeadTimeMS()
        {
            return fDeadTimeMS;
        }
		double getDeadTimeFraction( double iT_run_s = -99., bool iTimeDiffMethod = false, bool iCheckForConsistentDeadTime = true );
		double getDeadTimeFraction( vector< bool > iMask, bool iTimeDiffMethod = false, bool iCheckForConsistentDeadTime = true );
        TList* getDeadTimeHistograms();
        unsigned int getDeadTimeFraction_status()
        {
            return fDeadTimeFrac_status;   // 0=failed, 1=success
        }
        void   printDeadTime();
        bool   readHistograms( TDirectoryFile* iDir );
        void   reset();
        void   setDeadTimeConsistencyCheck_allowedDifference( double iDiff = 0.2 )
        {
            fDeadTimeConsistencyCheck_allowedDifference = iDiff;
        }
        void   writeHistograms( bool iDebug_IO = false );
        
        ClassDef( VDeadTime, 2 ) ;
};
#endif
