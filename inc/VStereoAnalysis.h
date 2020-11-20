//! VStereoAnalysis    class for producing one dimensional histograms from parameterized VERITAS data.

#ifndef VStereoAnalysis_H
#define VStereoAnalysis_H

#include "CData.h"
#include "VGammaHadronCuts.h"
#include "VAnaSumRunParameter.h"
#include "VAstronometry.h"
#include "VEvndispRunParameter.h"
#include "VTimeMask.h"
#include "VDeadTime.h"
#include "VEffectiveAreaCalculator.h"
#include "VExclusionRegions.h"
#include "VStereoHistograms.h"
#include "VStereoMaps.h"
#include "VSkyCoordinates.h"
#include "VSkyCoordinatesUtilities.h"

#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TParameter.h"
#include "TObject.h"
#include "TNamed.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VStereoAnalysis
{
    public:
    
        VStereoAnalysis( bool isOnrun, string i_hsuffix, VAnaSumRunParameter* irunpara,
                         vector< TDirectory* > iRDir, TDirectory* iDir, string iDataDir, int iRandomSeed, bool iTotalAnalysisOnly );
        ~VStereoAnalysis();
        double fillHistograms( int icounter, int irun, double AzMin = -1.e3, double AzMax = 1.e3, double iPedVar = -1. );
        TH2D*  getAlpha();
        TH2D*  getAlphaUC();
        TH2D*  getAlphaNorm();
        TH2D*  getAlphaNormUC();
        double getDeadTimeFraction();
        double getEffectiveExposure( int i_run )
        {
            return ( fRunExposure.find( i_run ) != fRunExposure.end() ? fRunExposure[i_run] : 0. );
        }
        TList* getEnergyHistograms();
        TList* getHisList();
        double getMeanAzimuth()
        {
            return fMeanAzimuth;
        }
        double getMeanElevation()
        {
            return fMeanElevation;
        }
        TH1D*  getMeanSignalBackgroundAreaRatio();
        TH1D* getMeanSignalBackgroundAreaRatioUC();
        double getMJD( int i_run )
        {
            return ( fRunMJD.find( i_run ) != fRunMJD.end() ? fRunMJD[i_run] : 0. );
        }
        double getMJDStart( int i_run )
        {
            return ( fRunMJDStart.find( i_run ) != fRunMJDStart.end() ? fRunMJDStart[i_run] : 0. );
        }
        double getMJDStopp( int i_run )
        {
            return ( fRunMJDStopp.find( i_run ) != fRunMJDStopp.end() ? fRunMJDStopp[i_run] : 0. );
        }
        TList* getParameterHistograms();
        double getRawRate();                      //! return number of entries in rate histograms
        vector< double > getRateCounts();
        vector< double > getRateTime();
        vector< double > getRateTimeIntervall();
        map< int, double > getRunMJD() const
        {
            return fRunMJD;
        }
        map< int, double > getRunDuration() const
        {
            return fRunDuration;
        }
        double getRunDuration( int i_run )
        {
            return ( fRunDuration.find( i_run ) != fRunDuration.end() ? fRunDuration[i_run] : 0. );
        }
        TList* getSkyHistograms( bool bUC );
        TH2D*  getStereoSkyMap();
        TH2D*  getStereoSkyMapUC();
        TH1D*  getTheta2();
        double getWobbleNorth();
        double getWobbleWest();
        TTree* getTreeWithSelectedEvents()
        {
            return fTreeSelectedEvents;
        }
        void   scaleAlpha( TH2D* halpha_on, bool bUC );
        void   setCuts( VAnaSumRunParameterDataClass iL, int irun );
        bool   setExclusionRegions( double i_raC, double i_decC );
        void   setNoSkyPlots( bool iS )
        {
            fNoSkyPlots = iS;
        }
        void   setRunExposure( map< int, double > iExpl )
        {
            fRunExposure = iExpl;
        }
        void   setRunMJD( map< int, double > iRunMJD )
        {
            fRunMJD = iRunMJD;
        }
        void   setRunTimes();
        string setRunTimes( CData* iData );
        bool   terminate();
        void   writeDebugHistograms();
        void   writeHistograms( bool bOn );
        
        
    private:
    
        bool fDebug;
        bool bIsGamma;
        bool bTotalAnalysisOnly;
        
        bool fIsOn;
        bool fNoSkyPlots;                         //! do full sky plots (if false, analysed source region only)
        
        VAnaSumRunParameter* fRunPara;
        
        TGraphAsymmErrors* gMeanEffectiveArea;
        TGraph2DErrors*    gTimeBinnedMeanEffectiveArea;
        TGraphErrors*      gMeanEsys_MC;
        
        TH2D* hMeanResponseMatrix;
        TGraphAsymmErrors* gMeanEffectiveAreaMC;

        VStereoMaps* fMap;
        VStereoMaps* fMapUC;
        int fHisCounter;
        map< int, double > fRunMJDStart;
        map< int, double > fRunMJDStopp;
        map< int, double > fRunMJD;               // Default value is mid-point; modified to mean time of accepted events by fillHistograms
        // If fRunMJD is defined from a VRunSummary.fRunMJD it will contain extra ON/OFF runs
        // because the VRunSummary.fRunMJD is a union of all the runs:
        // always access via the run index to be sure of getting the correct run
        map< int, double > fRunDuration;          // Raw run length from data tree
        map< int, double > fRunExposure;          // Open portion of time mask
        vector< VStereoHistograms* > fHisto;
        VStereoHistograms* fHistoTot;
        
        // tree containing events
        TTree* fTreeSelectedEvents;
        int    fTreeSelected_runNumber;
        int    fTreeSelected_eventNumber;
        int    fTreeSelected_MJD;
        double fTreeSelected_Time;
        int    fTreeSelected_NImages;
		ULong64_t fTreeSelected_ImgSel;
        double fTreeSelected_theta2;
        double fTreeSelected_Xoff;
        double fTreeSelected_Yoff;
        double fTreeSelected_Xoff_derot;
        double fTreeSelected_Yoff_derot;
        double fTreeSelected_Xcore;
        double fTreeSelected_Ycore;
        double fTreeSelected_MSCW;
        double fTreeSelected_MSCL;
        double fTreeSelected_MWR;
        double fTreeSelected_MLR;
        double fTreeSelected_Erec;
        double fTreeSelected_EChi2;
        double fTreeSelected_dE;
        float  fTreeSelected_EmissionHeight;
        float  fTreeSelected_EmissionHeightChi2;
        double fTreeSelected_SizeSecondMax;
        double fTreeSelected_Az;
        double fTreeSelected_El;
        UInt_t fTreeSelected_IsGamma;
        UInt_t fTreeSelected_DirectionCut;
        
		TTree* fDL3EventTree;
		int     fDL3EventTree_runNumber;
		int     fDL3EventTree_eventNumber;
		double  fDL3EventTree_Time;
		int     fDL3EventTree_MJD;
		double  fDL3EventTree_Xoff;
		double  fDL3EventTree_Yoff;
		double  fDL3EventTree_Xderot;
		double  fDL3EventTree_Yderot;
		double  fDL3EventTree_RA;
		double  fDL3EventTree_DEC;
		double  fDL3EventTree_Erec;
		double  fDL3EventTree_Erec_Err;
		double  fDL3EventTree_dE;
		double  fDL3EventTree_Xcore;
		double  fDL3EventTree_Ycore;
		int     fDL3EventTree_NImages;
		ULong64_t fDL3EventTree_ImgSel;
		double  fDL3EventTree_MSCW;
		double  fDL3EventTree_MSCL;
		double  fDL3EventTree_Az ;
		double  fDL3EventTree_El ;
		double  fDL3EventTree_EmissionHeight ;
		double  fDL3EventTree_Acceptance ;
		VRadialAcceptance* fDL3_Acceptance;
        
        double  fDeadTimeStorage ;
        //double fullMJD ;
        VSkyCoordinates* fVsky ;  // for RADec to AzimElev conversion
        
        double fTreeSelected_MVA;
        
        double fTotCount;
        
        map < int, double > f_t_in_s_min;
        map < int, double > f_t_in_s_max;
        double fMeanAzimuth;
        double fMeanElevation;
        double fNMeanElevation;
        
        CData* fDataRun;
        TTree* fDataRunTree;
        TFile* fDataFile;
        string fInstrumentEpoch;
        vector< unsigned int > fTelToAnalyze;
        
        vector< VSkyCoordinates* > fAstro;        //!< Astronomical source parameters for this analysis
        VGammaHadronCuts* fCuts;                  //!< Parameter Cuts
        VTimeMask* fTimeMask;                     //!< Time Cuts
        
        // dead time calculators
        vector< VDeadTime* > fDeadTime;
        
        // rate counters
        vector< vector< double > > fRateCounts;
        vector< vector< double > > fRateTime;
        vector< vector< double > > fRateTimeIntervall;
        vector< double > fRateCountsTot;
        vector< double > fRateTimeTot;
        vector< double > fRateTimeIntervallTot;
        
        // directories
        TDirectory* fDirTot;
        vector< TDirectory* > fDirTotRun;
        
        double combineHistograms();
        void   defineAstroSource();
        bool   closeDataFile();
        CData* getDataFromFile( int i_runNumber );
        
        void fill_TreeWithSelectedEvents( CData*, double, double, double, bool );
        bool init_TreeWithSelectedEvents( int, bool );
        void reset_TreeWithSelectedEvents();
        
		void fill_DL3Tree( CData* c , 
                                   double i_xderot, double i_yderot, 
                                   unsigned int icounter, double i_UTC );
		bool init_DL3Tree( int irun, int icounter );
		void write_DL3Tree();
        
        // derotation and J2000
        void getDerotatedCoordinates( unsigned int, CData* iData, double& x_derot, double& y_derot );
        
        int  getDataRunNumber() const;            // Check for existence of fDataRun and try to retrieve run number from first entry of the tree
};
#endif

