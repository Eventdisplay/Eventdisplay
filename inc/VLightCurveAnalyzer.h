//! VLightCurveAnalyzer analyze a light curve
#ifndef VLightCurveAnalyzer_H
#define VLightCurveAnalyzer_H

#include "TMath.h"

#include "VFluxAndLightCurveUtilities.h"
#include "VFluxDataPoint.h"

#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

class VLightCurveAnalyzer
{
    private:
    
        bool fDebug;
        
        // central data object
        vector< VFluxDataPoint > fFluxDataVector;
        
    public:
    
        VLightCurveAnalyzer();
        VLightCurveAnalyzer( vector< VFluxDataPoint > iDataVector );
        ~VLightCurveAnalyzer() {}
        
        void setDataVector( vector< VFluxDataPoint > iDataVector );
        void setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        
        double  get_Flux_Mean();
        double  get_Flux_Variance();
        double  get_Flux_VariabilityIndex( double iSystematicFraction = 1. );
};

#endif
