#ifndef VLowGainCalibrator_h
#define VLowGainCalibrator_h

#include "VGlobalRunParameter.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TMarker.h"

#include "TSystem.h"
#include "TH1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TText.h"
#include "TString.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"


using namespace std;

class VLowGainCalibrator
{

        //	private:
    public:
    
        enum status { GOOD, NO_POINTS, BAD_CHI2, NOT_PROPORTIONAL };
        bool fDEBUG;
        
        
        const static int fNTel = 4;
        const static int fNPix = 500;
        int fWindow;
        
        int fRun;
        
        TTree* fDsttree;
        
        int fChanMon_start;
        int fChanMon_stop;
        int fChan_start;
        int fChan_stop;
        
        int fNLiveMonitor_min;		//min number of good monitor channels to use the events
        bool fUseMedian;		//use median or mean of monitor charges in one event
        double fSumMonitor_min;		//min sum to be considered a monitor channel
        double fLightLevelWidth;	//Width of light levels (in terms of sigma)
        
        
        double fFitPure_min;		//min. purity of a light level to be considered for the fit. Eg 0.8 -> 80% of events have to be in high gain/low gain.
        double fFitNPoints_min;		//min. number of points for the fit (high/low is fitted separately)
        double fFitRSat_max;		//max. ratio of saturated high gain pixels in level.
        double fFitProb_min;		//min. fit probability assuming chi2 distribution.
        double fFitB_max;		//max of abs. value of b (if fitting a*x+b, not used currently).
        double fFitNEvents_min;		//min number of events in one level to make a data point for the fit.
        
        TH1D* fMonitorChargeHist[fNTel];
        unsigned int fNLightLevels[fNTel];
        
        vector<double> fLightLevelMean[fNTel];
        vector<double> fLightLevelSigma[fNTel];
        vector<double> fLightLevelMeanError[fNTel];
        vector<double> fLightLevelSigmaError[fNTel];
        
        vector<int> fN[fNTel][fNPix][2];	//number of entries per level, per tel, per pixel, per high/low gain.
        vector<int> fNSat[fNTel][fNPix][2];	//number of saturated traces per level, per tel, per pixel, per high/low gain.
        
        vector<double> fX[fNTel][fNPix][2];
        vector<double> fY[fNTel][fNPix][2];
        vector<double> fX2[fNTel][fNPix][2];
        vector<double> fY2[fNTel][fNPix][2];
        vector<double> fXY[fNTel][fNPix][2];
        
        double mHi[fNTel][fNPix] ;
        double mLo[fNTel][fNPix] ;
        
        
        TFile* fOutfile[fNTel];
        TFile* fDebugfile[fNTel];
        
        TFile* fInfile;
        
        TTree* fOuttree[fNTel];
        int fTree_Channel;
        double fTree_m[2];
        double fTree_mErr[2];
        double fTree_chi2[2];
        int fTree_ndf[2];
        status fTree_status[2];
        
        TTree* fDebugtree[fNTel];
        int fTree_eventNumber;
        double fTree_QMon;
        double fTree_QMonMean;
        int fTree_level;
        double fTree_Q;
        int fTree_hilo;
        int fTree_RawMax;
        int fTree_tel;
        int fTree_run;
        double fTree_TZero;
        
        vector<double> fLMult;
        vector<int> fDebugChannels;
        
        bool isNewPixel( int tel, int iChan );
        
        
        unsigned int eventNumber;
        float sumHiLo[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float sumMonitor[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        UShort_t sumfirst[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        UShort_t HiLo[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        UShort_t sumwindow[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        UShort_t dead[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        Short_t RawMax[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
        float pulsetiming[VDST_MAXTELESCOPES][VDST_MAXTIMINGLEVELS][VDST_MAXCHANNELS];
        
        double calcMonitorCharge( int tel, int ientry = -1 );
        double calcMeanMonitorCharge( int tel , int ientry = -1 );
        double calcMedianMonitorCharge( int tel , int ientry = -1 );
        
        bool fIsOk;
        
        unsigned int fMinDSTEvent;
        unsigned int fMaxDSTEvent;
        
        
        //	public:
        VLowGainCalibrator( int run, int sw, bool isInnerHigh, bool sw2monitor = false, TString dir = "./", TString outdir = "" );
        virtual ~VLowGainCalibrator();
        
        void setDebug( bool debug = true )
        {
            fDEBUG = debug ;
        }
        void resetLightLevels( );
        void resetLightLevels( int tel );
        void setMonitorChargeOptions( int nLive_min = 100, double sum_min = -100, bool useMedian = true, double width = 2 ) ;
        void setFitOptions( int n_min = 2, double pure_min = 0.8, double sat_max = 0.8, double nevents_min = 0, double prob_min = 0.01, double b_max = 2.0 );
        
        bool makeMonitorChargeHists( );
        int  checkLightLevels( int tel );
        bool calculateMeanCharges();
        
        bool doTheFit( ) ;
        
        bool terminate( );
        void findLightLevels( int tel, int iPeakSignificance = 2 , bool iDraw = true );
        bool findLightLevels( bool iDraw = true );
        void setLowGainMultiplierUsedInDST( double lmult = 6.0 )
        {
            fLMult.assign( fNTel, lmult );
        }
        
        void setDebugChannels( vector<int> channels )
        {
            fDebugChannels = channels ;
        }
        void addDebugChannel( int channel )
        {
            fDebugChannels.push_back( channel );
        }
        void setAllDebugChannels()
        {
            fDebugChannels.clear();
            for( int i = fChan_start; i < fChan_stop; i++ )
            {
                fDebugChannels.push_back( i );
            }
        }
        bool isDebugChannel( int channel )
        {
            return std::find( fDebugChannels.begin(), fDebugChannels.end(), channel ) != fDebugChannels.end();
        }
        
        void setDSTEventLimits( unsigned int a, unsigned int b )
        {
            fMinDSTEvent = a;
            fMaxDSTEvent = b;
        }
        
        ClassDef( VLowGainCalibrator, 5 );
        
};

#endif
