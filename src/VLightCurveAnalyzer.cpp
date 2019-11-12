/*  \class VLightCurveAnalyzer
    \brief analyze a light curve

    Note: take care which statistical methods you need to use:

    - fFlux/fFluxE: Poisson
    - fFlux_Rolke/fRate_lo_1sigma/fRate_up_1sigma: Rolke et al

*/

#include "VLightCurveAnalyzer.h"

VLightCurveAnalyzer::VLightCurveAnalyzer()
{
    fDebug = false;
    
}

VLightCurveAnalyzer::VLightCurveAnalyzer( vector< VFluxDataPoint > iDataVector )
{
    fDebug = false;
    
    setDataVector( iDataVector );
}

void VLightCurveAnalyzer::setDataVector( vector< VFluxDataPoint > iDataVector )
{
    fFluxDataVector = iDataVector;
}

/*

    calculate mean flux

*/
double VLightCurveAnalyzer::get_Flux_Mean()
{
    double iMean = 0.;
    double iNN = 0.;
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        iMean += fFluxDataVector[i].fFlux;
        iNN++;
    }
    if( iNN > 0. )
    {
        return iMean / iNN;
    }
    
    return -1.e99;
}

/*

    calculate the variance of the light curve

*/
double VLightCurveAnalyzer::get_Flux_Variance()
{
    double Sx = 0.;
    double Sxx = 0.;
    double iNN = 0.;
    
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        Sx  += fFluxDataVector[i].fFlux;
        Sxx += fFluxDataVector[i].fFlux * fFluxDataVector[i].fFlux;
        iNN++;
    }
    
    if( iNN > 1. )
    {
        return ( 1. / ( iNN - 1. ) * ( Sxx - 1. / iNN * Sx * Sx ) );
    }
    
    return 0.;
}


/*

   calculate variability index using a Chi2 criterion as described
   in section 4.5 in ApJS 188, 405, 2010 (1st FERMI LAT catalogue paper)

*/
double VLightCurveAnalyzer::get_Flux_VariabilityIndex( double iSystematicFraction )
{
    vector< double > F;
    vector< double > sigmaF;
    vector< double > w;
    double w_sum = 0.;
    double w_F = 0.;
    
    // fill vectors and calculate weights
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        F.push_back( fFluxDataVector[i].fFlux );
        sigmaF.push_back( fFluxDataVector[i].fFluxCI_1sigma );
        
        w.push_back( 1. / ( sigmaF.back()*sigmaF.back() + ( iSystematicFraction * F.back() ) * ( iSystematicFraction * F.back() ) ) );
        w_sum += w.back();
        w_F += F.back() * w.back();
    }
    
    // mean weighted flux
    double F_mean = 0.;
    if( w_sum > 0. )
    {
        F_mean = w_F / w_sum;
    }
    
    // calculate variability index
    double V = 0.;
    for( unsigned int j = 0; j < F.size(); j++ )
    {
        V += w[j] * ( F[j] - F_mean ) * ( F[j] - F_mean );
    }
    
    return V;
}

