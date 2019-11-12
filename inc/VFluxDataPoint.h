//! data element and algorithm for flux calculation / light curve
#ifndef VFluxData_H
#define VFluxData_H

#include "TMath.h"
#include "VFluxAndLightCurveUtilities.h"
#include "VOrbitalPhaseData.h"
#include "VStatistics.h"
#include "TRandom3.h"
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

using namespace std;

class VFluxDataPoint
{
    private:
    
        // significance and upper flux limit parameters
        int    fLiMaEqu;
        double fThresholdSignificance;
        double fMinEvents;
        double fUpperLimit;
        int    fUpperLimitMethod;
        bool   fBoundedLimits;
        // variables used for operator+
        double fAlpha_U;             //!< temporary result used for alpha calculation
        double fAlpha_L;             //!< temporary result used for alpha calculation
        double fTotFluxSum;          //!< temporary result used for flux calculation
        double fTotFluxWeight;       //!< temporary result used for flux calculation
        double fTotZeSum;            //!< temporary result used for mean zenith
        double fTotWobbleOffset;     //!< temporary result used for mean wobble offset
        double fTotPedvars;          //!< temporary result used for mean pedvar
        double fTotEffArea;          //!< temporary result used for mean effective area
        
        unsigned int fDebugCounter;
        
    public:
    
        // more (a lot of) printouts
        bool   fDebug;
        
        // general info
        string fName;                //!< general name given
        string fDataFileName;        //!< name of anasum file used for this analysis
        bool   fIsCombinedFluxData;  //!< data point is a combination of several
        
        // run information
        int    fRunNumber;           //!< run number
        double fMJD_RunStart;        //!< MJD at start of run [days]
        double fMJD_RunStop;         //!< MJD at end of run [days]
        set< int > fRunNumber_list;  //!< list of of all runs used for this data point
        
        // time bin data
        double fMJD_Start;           //!< MJD at start of time bin (equivalent to fMJD_RunStart for run-wise calculation) [days]
        double fMJD_Stop;            //!< MJD at end of time bin (equivalent to fMJD_RunStart for run-wise calculation) [days]
        double fMJD;                 //!< MJD (centre of time bin) [days]
        double fOrbitalPhase_Start;  //!< orbital phase at start of time bin (in case orbital phase parameters are given)
        double fOrbitalPhase_Stop;   //!< orbital phase at end of time bin (in case orbital phase parameters are given)
        double fOrbitalPhase;        //!< orbital phse (centre of time bin, in case orbital phase parameters are given)
        double fSecondsIntoRun;      //!< seconds into the run (centre of time bin) [s]
        double fTimeBinDuration_sec; //!< length of time bin (equivalent to length of run for run-wise calculation) [s]
        bool   fTimeMask_open;       //!< this time bin might be excluded from analysis (false = good data bin)
        
        // life time data
        double fExposure;            //!< life time [s]
        double fExposure_deadTimeCorrected;  //! dead time corrected life time [s]
        double fDeadTime;            //!< dead time fraction
        
        double fZe;                  //!< mean zenith angle [deg]
        double fWobbleOffset;        //!< wobble offset [deg]
        double fPedvars;             //!< pedvars
        
        // events per time bin
        double fNon;                 //!< N_on
        double fNoff;                //!< N_off
        double fAlpha;               //!< alpha
        double fSignificance;        //!< significance
        bool   fSignificantDataPoint; //!< data point passes all significances/nevents requirements
        
        // events per time time: signal calculation
        double fNdiff;               //!< N_diff = N_on - alpha * N_off
        double fNdiffE;              //!< Poisson error in N_on - alpha * N_off
        double fNdiff_Rolke;         //!< estimated diff after Rolke et al
        double fCI_lo_1sigma;        //!< counts: lower value of 1 sigma confidence interval
        double fCI_up_1sigma;        //!< counts: upper value of 1 sigma confidence interval
        double fCI_lo_3sigma;        //!< counts: lower value of 3 sigma confidence interval
        double fCI_up_3sigma;        //!< counts: upper value of 3 sigma confidence interval
        double fUL;                  //!< counts: upper flux limit in events (-99. if not set)
        
        int    fCalculateExpectedLimitsN ;
        double fUL_Expected;
        double fUL_Expected_lo_1sigma;
        double fUL_Expected_up_1sigma;
        double fUL_Expected_lo_2sigma;
        double fUL_Expected_up_2sigma;
        
        double fFluxUL_Expected;
        double fFluxUL_Expected_lo_1sigma;
        double fFluxUL_Expected_up_1sigma;
        double fFluxUL_Expected_lo_2sigma;
        double fFluxUL_Expected_up_2sigma;
        
        double fFluxConstantUL_Expected;
        double fFluxConstantUL_Expected_lo_1sigma;
        double fFluxConstantUL_Expected_up_1sigma;
        double fFluxConstantUL_Expected_lo_2sigma;
        double fFluxConstantUL_Expected_up_2sigma;
        
        
        // rate calculation
        double fRate;                //!< gamma-ray rate in [1/min]
        double fRateE;               //!< error in gamma-ray rate in [1/min]
        double fRate_Rolke;          //!< gamma-ray rate in [1/min] (after Rolke et al)
        double fRate_lo_1sigma;      //!< rate: lower value of 1 sigma confidence interval
        double fRate_up_1sigma;      //!< rate: upper value of 1 sigma confidence interval
        
        // flux calculation
        double fEffArea_cm2;         //!< integrated effective area [cm^2]
        double fWeightedEffArea_cm2;         //!< normalized effective area [cm^2]
        double fFluxConstant;        //!< flux constant or upper flux limit constant
        double fFluxConstantE;       //!< error in flux constant or upper flux limit constant
        double fFlux;                //!< flux
        double fFluxE;               //!< flux: Poissonian flux in error
        double fFlux_Rolke;          //!< flux using diff after Rolke et al
        double fFluxCI_1sigma;       //!< flux: mean 1 sigma confidence interval
        double fFluxCI_lo_1sigma;    //!< flux: lower value of 1 sigma confidence interval
        double fFluxCI_up_1sigma;    //!< flux: upper value of 1 sigma confidence interval
        double fFluxCI_lo_3sigma;    //!< flux: lower value of 3 sigma confidence interval
        double fFluxCI_up_3sigma;    //!< flux: upper value of 3 sigma confidence interval
        double fFluxUL;              //!< upper flux limit
        string fFluxUnitString;      //!< unit string (used in axis for plotting)
        
        // spectral parameters used for flux calculation
        double fMinEnergy_TeV;
        double fMaxEnergy_TeV;
        double fE0_TeV;
        double fSpectralIndex;
        
        // orbital phase data
        VOrbitalPhaseData fOrbitalPhaseData;
        
        VFluxDataPoint( string iName = "flux_data_point" );
        virtual ~VFluxDataPoint() {}
        
        void calculateFlux();
        bool calculateIntegralEffectiveArea( vector< double > energy_axis, vector< double > effArea );
        void calculateOrbitalPhaseData();
        void calculateOrbitalPhaseData( VOrbitalPhaseData iOrbitalPhaseData );
        void calculateSignificancesAndUpperLimits();
        bool isCombinedDataElement()
        {
            return fIsCombinedFluxData;
        }
        bool isSignificantDataPoint()
        {
            return fSignificantDataPoint;
        }
        bool isTimeInsideRun( double iT_seconds );
        bool hasOrbitalPhases();
        void printRunData( ostream& output = std::cout );
        void printResults( ostream& output = std::cout );
        void printResultSummary( ostream& output = std::cout );
        void reset();
        void resetFluxValues();
        void setCombinedDataElement()
        {
            fIsCombinedFluxData = true;
        }
        void setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        void setSpectralParameters( double iMinEnergy_TeV = 0., double E0 = 1., double alpha = -2.5, double iMaxEnergy_TeV = 300. );
        void setSignificanceParameters( double iThresholdSignificance = 3., double iMinEvents = 5, double iUpperLimit = 0.99,
                                        int iUpperlimitMethod = 5, int iLiMaEqu = 17, bool iBounded = true );
                                        
        void calculateExpectedLimit( );
        
        VFluxDataPoint operator+( VFluxDataPoint& s );
        VFluxDataPoint& operator+=( VFluxDataPoint& s );
        void setCalculateExpectedLimitsN( int n = 100000 )
        {
            fCalculateExpectedLimitsN = n ;
        }
        bool operator<( const VFluxDataPoint& s1 ) const;
        
        ClassDef( VFluxDataPoint, 6 );
};

#endif
