//! VMathsandFunctions collects math utilities and external function to be used by TF1/TF2/etc

#ifndef VMathsandFunctions_h
#define VMathsandFunctions_h

#include "VAnalysisUtilities.h"
#include "VPlotUtilities.h"

#include "TMath.h"

using namespace std;

class VDouble_gauss : public VAnalysisUtilities, public VPlotUtilities
{
    public:
    
        double operator()( double* x, double* par );
        
        ClassDef( VDouble_gauss, 1 );
        
};

class VFun_gauss : public VAnalysisUtilities, public VPlotUtilities
{
    public:
    
        double operator()( double* x, double* par );
        
        ClassDef( VFun_gauss, 1 );
};

namespace VMathsandFunctions
{
    double getMeanEnergyInBin( unsigned int iMethod, double e_min_log10, double e_max_log10, double iSpectralIndex );
    double getBaryCentricMeanEnergy( double e_min_log10, double e_max_log10, double iSpectralIndex );
    double getMeanEnergy( double e_min_log10, double e_max_log10 );
    double getRatioError( double x1, double x2, double ex1, double ex2 );   // x1 / x2
    double getSpectralWeightedMeanEnergy( double e_min_log10, double e_max_log10, double iSpectralIndex );
}

#endif

