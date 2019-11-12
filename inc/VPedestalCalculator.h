//! VPedestalCalculator pedestal calculation in time slices

#ifndef VPEDESTALCALCULATOR_H
#define VPEDESTALCALCULATOR_H

#include <VImageBaseAnalyzer.h>
#include <VGlobalRunParameter.h>
#include <VSkyCoordinatesUtilities.h>

#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VPedestalCalculator : public VImageBaseAnalyzer
{
    private:
        bool fDebug;
        
        bool bCalibrationRun;
        
        double fLengthofTimeSlice;
        int fSumFirst;
        int fSumWindow;
        
        int runNumber;
        int MJD;
        double time;
        unsigned int fNPixel;
        double x[VDST_MAXCHANNELS];
        double y[VDST_MAXCHANNELS];
        double xRot[VDST_MAXCHANNELS];
        double yRot[VDST_MAXCHANNELS];
        
        vector< TTree* > fTree;
        
        vector< vector< vector< float > > > fpedcal_n;
        vector< vector< vector< float > > > fpedcal_mean;
        vector< vector< vector< float > > > fpedcal_mean2;
        
        vector< vector< float > > v_temp_pedEntries;
        vector< vector< float > > v_temp_ped;
        vector< vector< float > > v_temp_pedvar;
        
        // timing vector
        vector< double > fTimeVec;
        
        double adjustTimeSliceLength( double iLengthofTimeSlice, double iRunStartTime, double iRunStoppTime );
        void fillTimeSlice( unsigned int );
        void reset();
        
    public:
    
        vector< vector< int > > v_MJD;            //! [telid][time slice]
        vector< vector< double > > v_time;        //! [telid][time slice]
        //! [telid][time slice][npixel][summation window]
        vector< vector< vector< vector< float > > > > v_pedEntries;
        //! [telid][time slice][npixel][summation window]
        vector< vector< vector< vector< float > > > > v_ped;
        //! [telid][time slice][npixel][summation window]
        vector< vector< vector< vector< float > > > > v_pedvar;
        
        VPedestalCalculator();
        ~VPedestalCalculator() {}
        
        void doAnalysis( bool iLowGain = false );
        vector< TTree* > getPedestalTree()
        {
            return fTree;
        }
        bool initialize();
        bool initialize( bool ibCalibrationRun, unsigned int iNPixel, double iLengthofTimeSlice, int iSumFirst, int iSumWindow,
                         double iRunStartTime = -99., double iRunStoppTime = -99. );
        void terminate( bool iWrite = true, bool bDebug_IO = false );
};
#endif
