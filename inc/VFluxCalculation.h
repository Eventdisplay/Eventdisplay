//! VFluxCalculation calculate fluxes and upper flux limits
#ifndef VFLUXCALCULATION_H
#define VFLUXCALCULATION_H

#include "TDirectory.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TH2D.h"
#include "TMath.h"
#include "TObject.h"
#include "TTree.h"

#include "CRunSummary.h"

#include "VFluxAndLightCurveUtilities.h"
#include "VFluxDataPoint.h"
#include "VOrbitalPhaseData.h"
#include "VStatistics.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

class VFluxCalculation : public TObject
{
    private:
        bool fDebug;
        
        // input data
        vector< TFile* > fFile;                   //!< list of anasum data files
        bool bZombie;                             //!< no file or invalid file connected
        CRunSummary* fData;                       //!< run summary tree
        
        // data vectors used internally
        vector< VFluxDataPoint > fFluxData_perRun;      //! data vector for analysis results (one entry per run) [runnumber]
        vector< VFluxDataPoint > fFluxData_perTimeBin;  //! data vector for analysis results (one entry per time bin) [time bin]
        
        
        bool fUseIntraRunBins;
        bool fUseRunWiseBins;
        
        // data vectors with results of analysis
        double fRequestedTimeBins_MJD_start_value;	//MJD value at which to start the time bins. if 0, start at the start of the first run. If <0, start at the start of the run always.
        double fRequestedTimeBinWidth_s;                //! width of time bins in the analysis
        //! (note: restricted to multiple of time bins used in anasum or to
        //!        time bins much larger than a typical run length (e.g. a day))
        //! (use -1 for run-wise analysis)
        vector< double > fRequestedTimeBins_MJD_start;  //! vector with time periods to be used in flux analysis [MJD]
        vector< double > fRequestedTimeBins_MJD_stopp;  //! vector with time periods to be used in flux analysis [MJD]
        
        vector< double > fRequestedPhaseBins_start;  //! vector with time periods to be used in flux analysis [binary phase]
        vector< double > fRequestedPhaseBins_stopp;  //! vector with time periods to be used in flux analysis [binary phase]
        
        vector< VFluxDataPoint > fFluxDataVector;       //! data vector after rebinning in time [time bin] or phase
        VFluxDataPoint           fFluxDataCombined;     //! combined flux data point
        
        // maximum energy for integration of spectra
        double fMaxSave_MCEnergy_TeV;
        // spectral parameters for flux calculation (assuming power law)
        double fMinEnergy_TeV;                        //!< calculate flux limit above this energy [TeV]
        double fMaxEnergy_TeV;                        //!< maximum energy to be taken into account [TeV]
        double fE0_TeV;                               //!< flux normalization energy [TeV]
        double fSpectralIndex;                        //!< assumed spectral index
        
        // significance and upper flux limit parameters
        int    fLiMaEqu;
        double fThresholdSignificance;
        double fMinEvents;
        double fUpperLimit;
        int    fUpperLimitMethod;
        bool   fBoundedLimits;
        int    fCalculateExpectedLimitsN;		//number of loops for calculating expected limits. If n<=0, no expected limits will be calculated
        // orbital phase data
        VOrbitalPhaseData fOrbitalPhaseData;
        
        // private functions
        bool   calculateCombinedFluxes();
        bool   calculateFluxes();
        bool   calculateSignificancesAndUpperLimits();
        void   closeFiles();
        bool   getIntegralEffectiveArea();
        bool   getNumberOfEventsinEnergyInterval();
        bool   getNumberOfEventsinEnergyInterval( TFile* iFile );
        bool   openAnasumDataFile( string ifile );
        bool   openAnasumDataFile( vector< string > ifile );
        void   reset();
        void   resetRunList();
        
    public:
    
        VFluxCalculation();
        VFluxCalculation( string iDataFile,
                          int iRunMin = 0, int iRunMax = 100000,
                          double iMJDMin = -99., double iMJDMax = -99.,
                          double iFluxMultiplier = 1.,
                          bool iDebug = false );
        VFluxCalculation( vector< string > iFile_vector,
                          int iRunMin = 0, int iRunMax = 100000,
                          double iMJDMin = -99., double iMJDMax = -99.,
                          bool iDebug = false );
        ~VFluxCalculation();
        
        bool   calculateIntegralFlux( double iMinEnergy_TeV = 1. );
        bool   IsZombie()
        {
            return bZombie;
        }
        const VFluxDataPoint*    getFluxDataCombined()
        {
            return &fFluxDataCombined;
        }
        vector< VFluxDataPoint > getFluxDataVector()
        {
            return fFluxDataVector;
        }
        vector< double >* getRequestedTimeBins_MJD_start()
        {
            return &fRequestedTimeBins_MJD_start;
        }
        vector< double >* getRequestedTimeBins_MJD_stopp()
        {
            return &fRequestedTimeBins_MJD_stopp;
        }
        vector< double >* getRequestedPhaseBins_start()
        {
            return &fRequestedPhaseBins_start;
        }
        vector< double >* getRequestedPhaseBins_stopp()
        {
            return &fRequestedPhaseBins_stopp;
        }
        
        const VFluxDataPoint*    getFluxDataPerRun( int iRunNumber );
        unsigned int  loadFluxDataVectorFromAsciiFile( string iAsciiFile,
                double iFluxMultiplier = 1.,
                double iMJDMin = -99., double iMJDMax = -99. );
        unsigned int  loadRunListFromAnasumDataFile( int iRunMin = 0, int iRunMax = 100000,
                double iMJDMin = -99., double iMJDMax = -99. );
        void          printResults( ostream& output = std::cout );
        void          printResultSummary( ostream& output = std::cout );
        void          printRunList( ostream& output = std::cout );
        void          printTimeBinWidth();
        void          setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        void          setMaxSaveMC_Energy_TeV( double iEnergy_TeV = 300. )
        {
            fMaxSave_MCEnergy_TeV = iEnergy_TeV;
        }
        void          setPhaseFoldingValues( double iZeroPhase_MJD = -99., double iOrbit_Days = -99.,
                                             double iOrbitError_low_Days = 0., double iOrbitError_up_Days = 0. );
        void          setSpectralParameters( double iMinEnergy_TeV = 0., double E0_TeV = 1.,
                                             double alpha = -2.5, double iMaxEnergy_TeV = 300. );
        void          setSignificanceParameters( double iThresholdSignificance = 2., double iMinEvents = 5,
                double iUpperLimit = 0.99, int iUpperlimitMethod = 5, int iLiMaEqu = 17,
                bool iBoundedLimits = true );
        bool          setTimeBinVector( string iTimeBinFile, bool iPrint = false );
        void          setTimeBinVector( vector< double > iRequestedTimeBins_MJD_start, vector< double > iRequestedTimeBins_MJD_stopp );
        void          setPhaseBinVector( vector< double > iRequestedPhaseBins_start, vector< double > iRequestedPhaseBins_stopp );
        void 	      setPhaseBins( int nBins, double iStart = 0 );
        void	      setCalculateExpectedLimitsN( int n = 10000 )
        {
            fCalculateExpectedLimitsN = n ;
        }
        void 	      calculateExpectedLimitCombinedOnly( int n = 10000 );
        void	      setRunwiseLightCurve();
        void	      setTimeBinsInDays( double nNights = 1, double iMJD_start = 0 );
        void	      setTimeBinsInSeconds( double seconds = 0 );
        void	      clearTimeBinVectors();
        void	      clearPhaseBinVectors();
        
        ClassDef( VFluxCalculation, 18 );
};
#endif
