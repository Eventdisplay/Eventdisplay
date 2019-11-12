/*! \class VMathsandFunctions
    \brief collects math utilities and external function to be used by TF1/TF2/etc

*/

#include "VMathsandFunctions.h"

ClassImp( VDouble_gauss )
ClassImp( VFun_gauss )

///////////////////////////////////////////////////////////////////////////////
/*
    1D two-sided normal distribution
*/
double VDouble_gauss::operator()( double* x, double* par )
{

    double norm = par[0];
    double mean = par[1];
    double sigma1 = par[2];
    double sigma2 = par[3];
    
    if( sigma1 <= 0. || sigma2 <= 0. )
    {
        return 0.;
    }
    
    if( x[0] < mean )
    {
        return norm * TMath::Exp( -0.5 * ( x[0] - mean ) * ( x[0] - mean ) / sigma1 / sigma1 );
    }
    return norm * TMath::Exp( -0.5 * ( x[0] - mean ) * ( x[0] - mean ) / sigma2 / sigma2 );
    
}


///////////////////////////////////////////////////////////////////////////////
/*
     1D Normal distributions
*/
double VFun_gauss::operator()( double* x, double* par )
{

    return par[0] * TMath::Exp( -1 * x[0] * x[0] / 2. );
    
}


/*!
    calculate spectral weighted mean energy in an energy bin

    e_min_log10:  lower bin limit (in log10)
    e_max_log10:  upper bin limit (in log10)
    g:  spectral index assumed (negative value)

    return value: spectral weighted mean energy in logarithmic values (log10)

    see as well Lafferty and Wyatt (1995), NIMA 355, 541
*/
double VMathsandFunctions::getSpectralWeightedMeanEnergy( double e_min_log10, double e_max_log10, double iSpectralIndex )
{
    double xL = TMath::Power( 10., e_min_log10 );
    double xU = TMath::Power( 10., e_max_log10 );
    
    // return mean (log value, mean on the linear scale)
    if( iSpectralIndex == 0. )
    {
        return TMath::Log10( 0.5 * ( xL + xU ) );
    }
    
    // upper limit of bin should be larger than lower limit
    if( xU <= xL )
    {
        return e_min_log10;
    }
    
    // calculate weighted mean
    double xM = 0.;
    xM = 1. / ( iSpectralIndex + 1. ) / ( xU - xL ) * ( TMath::Power( xU, iSpectralIndex + 1. ) - TMath::Power( xL, iSpectralIndex + 1. ) );
    xM = TMath::Log10( xM ) / iSpectralIndex;
    
    return xM;
}

/*

   calculate barycentre of a bin

    e_min_log10:  lower bin limit (in log10)
    e_max_log10:  upper bin limit (in log10)
    g:  spectral index assumed (negative value)

    return value: spectral weighted mean energy in logarithmic values (log10)
*/
double VMathsandFunctions::getBaryCentricMeanEnergy( double e_min_log10, double e_max_log10, double iSpectralIndex )
{
    double xL = TMath::Power( 10., e_min_log10 );
    double xU = TMath::Power( 10., e_max_log10 );
    
    double xM = -99.e99;
    // normalisation
    xM  = ( iSpectralIndex + 1. ) / ( TMath::Power( xU, iSpectralIndex + 1. ) - TMath::Power( xL, iSpectralIndex + 1. ) );
    // integration
    xM *=  1. / ( iSpectralIndex + 2. ) * ( TMath::Power( xU, iSpectralIndex + 2. ) - TMath::Power( xL, iSpectralIndex + 2. ) );
    
    if( xM > 0. )
    {
        xM = log10( xM );
    }
    
    return xM;
}

/*

    calculate mean energy of a bin (log bin)

    e_min_log10:  lower bin limit (in log10)
    e_max_log10:  upper bin limit (in log10)
*/
double VMathsandFunctions::getMeanEnergy( double e_min_log10, double e_max_log10 )
{
    double xL = TMath::Power( 10., e_min_log10 );
    double xU = TMath::Power( 10., e_max_log10 );
    
    return TMath::Log10( 0.5 * ( xL + xU ) );
}

double VMathsandFunctions::getMeanEnergyInBin( unsigned int iMethod, double e_min_log10, double e_max_log10, double iSpectralIndex )
{
    if( iMethod == 0 )
    {
        return getMeanEnergy( e_min_log10, e_max_log10 );
    }
    else if( iMethod == 1 )
    {
        return getBaryCentricMeanEnergy( e_min_log10, e_max_log10, iSpectralIndex );
    }
    else if( iMethod == 2 )
    {
        return getSpectralWeightedMeanEnergy( e_min_log10, e_max_log10, iSpectralIndex );
    }
    
    return -1.e99;
}

/*

  simple error propagation for f=x1/x2

*/
double VMathsandFunctions::getRatioError( double x1, double x2, double ex1, double ex2 )
{
    if( x2 != 0. )
    {
        return TMath::Sqrt( 1. / x2 / x2 * ex1 * ex1 + x1 * x1 / x2 / x2 / x2 / x2 * ex2 * ex2 );
    }
    
    return 0.;
}
