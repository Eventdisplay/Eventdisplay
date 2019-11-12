/*  VMedianCalculator
 *
 *  efficient median calculator
 *
 *  data set of up to nDim_exact elements:     precise median calculation
 *  data set of more than nDim_exact elements: approximation
 *
 */

#include "VMedianCalculator.h"

VMedianCalculator::VMedianCalculator()
{
    reset();
}

void VMedianCalculator::reset()
{
    n_counter = 0.;
    x.clear();
    
    mean_x  = 0.;
    mean_xx = 0.;
    
    for( unsigned int i = 0; i < 3; i++ )
    {
        quantiles[i] = 0.;
    }
    prob[0] = 0.16;
    prob[1] = 0.50;
    prob[2] = 0.84;
    
    setNExact( 10000 );
    setEta();
    
    x.resize( 1000 );
}

void VMedianCalculator::fill( double ivalue )
{
    // approximation: use just first nDim_exact events
    if( n_counter < nDim_exact )
    {
        x.push_back( ivalue );
    }
    
    // mean and rms
    mean_x  += ivalue;
    mean_xx += ivalue * ivalue;
    
    n_counter++;
}

double VMedianCalculator::getMean()
{
    if( n_counter > 0 )
    {
        return mean_x / ( ( float )n_counter );
    }
    
    return 0.;
}

double VMedianCalculator::getRMS()
{
    if( n_counter > 1 )
    {
        return 1. / ( ( float( n_counter ) - 1. ) * ( mean_xx - mean_x * mean_x ) );
    }
    
    return 0.;
}

double VMedianCalculator::getMedian( float& medianwidth, int& N )
{
    double* i_x = &x[0];
    int i_n = ( int )x.size();
    double i_a[] = { prob[0], prob[1], prob[2] };
    double i_b[] = { 0.0,  0.0, 0.0  };
    TMath::Quantiles( i_n, 3, i_x, i_b, i_a, kFALSE );
    
    medianwidth = float( i_b[2] - i_b[0] );
    N = n_counter;
    
    return i_b[1];
}

double VMedianCalculator::getMedianWidth()
{
    double* i_x = &x[0];
    int i_n = ( int )x.size();
    double i_a[] = { prob[0], prob[1], prob[2] };
    double i_b[] = { 0.0,  0.0, 0.0  };
    TMath::Quantiles( i_n, 3, i_x, i_b, i_a, kFALSE );
    return ( i_b[2] - i_b[0] );
}
