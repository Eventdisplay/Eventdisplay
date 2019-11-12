//! VDifferentialFluxData data and converter class for differential flux values

#ifndef VDifferentialFluxData_H
#define VDifferentialFluxData_H

#include "TMath.h"
#include "TObject.h"

#include "VFluxAndLightCurveUtilities.h"

#include <iomanip>
#include <iostream>

using namespace std;

class VDifferentialFluxData : public TObject
{
    public:
    
        double MJD_min;                           // time range of observations
        double MJD_max;                           // time range of observations
        double Energy;                            // [TeV]
        double Energy_lowEdge;                    // lower bound of energy bin [TeV]
        double Energy_upEdge;                     // upper bound of energy bin [TeV]
        unsigned int Energy_lowEdge_bin;
        unsigned int Energy_upEdge_bin;
        double Energy_Hz;                         // [Hz]
        double EnergyWeightedMean;                // [TeV]
        double dE;                                // size of energy bin (in TeV)
        double EffectiveArea;                     // effective area in cm^2
        double DifferentialFlux;                  // [1/cm2/s/TeV]
        double DifferentialFluxError;             // [1/cm2/s/TeV]
        double DifferentialFluxError_up;          // [1/cm2/s/TeV]
        double DifferentialFluxError_low;         // [1/cm2/s/TeV]
        double DifferentialFlux_vFv;              // vF_v [ergs/cm2/s]
        double DifferentialFluxError_vFv;         // error in vF_v [ergs/cm2/s]
        double DifferentialFluxError_up_vFv;      // error in vF_v [ergs/cm2/s]
        double DifferentialFluxError_low_vFv;     // error in vF_v [ergs/cm2/s]
        double DifferentialFluxWeighting;         // weight used for the combination of fluxes
        double ObsTime;                           // observation time [s] (dead time corrected)
        double ExposureTime;                      // exposure [s]
        double NOn;
        double NOn_error;
        double NOff;
        double NOff_error;
        double NOff_alpha;
        double Significance;
        
        VDifferentialFluxData();
        ~VDifferentialFluxData() {}
        void fillEvent( double iMinMJD = 0., double iMaxMJD = 1.e14 );  // calculate vFv etc., fill MJDs, etc.
        void print( bool bSED = false );
        void printClean( bool bSED = false, bool b_PrintError_low_up = true );
        void resetCountsAndFluxes();
        
        double nuFnu( double F, double gamma, double e1, double e2, double e3 = -9999. );
        
        ClassDef( VDifferentialFluxData, 10 );
};
#endif
