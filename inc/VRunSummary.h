//! VRunSummary  class for rate/significance results from each run

#ifndef VRUNSUMMARY_H
#define VRUNSUMMARY_H

#include "CRunSummary.h"
#include "VAnaSumRunParameter.h"
#include "VSkyCoordinatesUtilities.h"

#include <iomanip>
#include <iostream>
#include <string>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TMinuit.h"
#include "TTree.h"

using namespace std;

class VRunSummary
{
    private:
        TTree* fRunSummaryTree;
        
        void fitEnergyHistogram( TH1D*, double minE, double maxE );
        void init();
        bool initTree();
        bool setBranches();
        
    public:
    
        int runOn;
        double MJDOn;              // Default value is mid-point; modified to mean time of accepted events by fillHistograms
        double MJDOn_runStart;
        double MJDOn_runStopp;
        double RunDurationOn;
        Char_t fTargetName[300];
        
        int runOff;
        double MJDOff;
        double MJDOff_runStart;
        double MJDOff_runStopp;
        double RunDurationOff;
        
        double fTargetRA;
        double fTargetDec;
        double fTargetRAJ2000;
        double fTargetDecJ2000;
        double fSkyMapCentreRAJ2000;
        double fSkyMapCentreDecJ2000;
        double fTargetShiftRAJ2000;
        double fTargetShiftDecJ2000;
        double fTargetShiftNorth;
        double fTargetShiftWest;
        double fWobbleNorth;
        double fWobbleWest;
        unsigned int fNTel;
        Char_t fTelList[300];
        double elevationOn;
        double azimuthOn;
        double elevationOff;
        double azimuthOff;
        double fTheta2Max;
        double RawRateOn;
        double RawRateOff;
        double pedvarsOn;
        double pedvarsOff;
        double NOn;
        double NOff;
        double NOffNorm;
        double OffNorm;
        double Signi;
        double Rate;
        double RateE;
        double RateOff;
        double RateOffE;
        double DeadTimeFracOn;
        double DeadTimeFracOff;
        double MaxSigni;
        double MaxSigniX;
        double MaxSigniY;
        
        ////////////////
        // variables used for total analysis only
        //
        map< int, double >  fRunMJD;
        double fTotalExposureOn;
        double fTotalExposureOff;
        map< int, double > f_exposureOn;
        map< int, double > f_exposureOff;
        double fMeanElevationOn;
        double fMeanElevationOff;
        double fMeanAzimuthOn;
        double fMeanAzimuthOff;
        double fNMeanElevation;
        double fMeanDeadTimeOn;
        double fMeanDeadTimeOff;
        double fMeanRawRateOn;
        double fMeanRawRateOff;
        double fMeanPedVarsOn;
        double fMeanPedVarsOff;
        
        double fTotTargetRA;
        double fTotTargetDec;
        double fTotTargetRAJ2000;
        double fTotTargetDecJ2000;
        
        VRunSummary();
        ~VRunSummary() {}
        void fill();
        bool fill( string, string, vector< VAnaSumRunParameterDataClass > );
        void print();
        void write();
};
#endif
