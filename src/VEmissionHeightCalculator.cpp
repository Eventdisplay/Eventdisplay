/*! \class VEmissionHeightCalculator
    \brief  calculate mean height of maximum emission (in shower coordinates)

*/

#include "VEmissionHeightCalculator.h"

VEmissionHeightCalculator::VEmissionHeightCalculator()
{
    fDebug = false;
    
    fNTel = 0;
    fNTelPairs = 0;
    fEmissionHeight = 0.;
    fEmissionHeightChi2 = 0.;
    fEmissionHeightT.assign( 1000, -99. );
}


VEmissionHeightCalculator::~VEmissionHeightCalculator()
{
}


double VEmissionHeightCalculator::getEmissionHeight( double* cen_x, double* cen_y, double* size, double az, double el )
{
    double iEmissionHeight = 0.;
    double iEmissionHeightWeight = 0.;
    double iEmissionHeightWeightTemp = 0.;
    double iEmissionHeight2 = 0.;
    double iEmissionHeightTemp = 0.;
    double iNEM_pairs = 0.;
    
    double fTelescopeDistanceSC = 0.;
    double fImageDistance = 0.;
    // counter for telescope pairs
    int nTPair = 0;
    
    // reset emission heights
    fEmissionHeight = 0.;
    fEmissionHeightChi2 = 0.;
    // loop over all telescope pairs
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        for( unsigned int j = i; j < fNTel; j++ )
        {
            // require elevation > 0. (test if telescope is present in analysis)
            if( i != j )
            {
                if( el > 0. && size[i] > 0. && size[j] > 0. )
                {
                    // get tangens of distance between the two image centroids
                    fImageDistance = TMath::Tan( imageDistance( cen_x[i], cen_x[j], cen_y[i], cen_y[j] ) / TMath::RadToDeg() );
                    if( fImageDistance > 0. )
                    {
                        // get distance between the two telescopes in shower coordinates
                        fTelescopeDistanceSC = getTelescopeDistanceSC( i, j, az, el );
                        // calculate emission height [km]
                        iEmissionHeightTemp = fTelescopeDistanceSC / fImageDistance / 1.e3;
                        // weight for pairwise emission height calculation
                        iEmissionHeightWeightTemp = 1. / ( ( 1. / log10( size[i] ) ) + ( 1. / log10( size[j] ) ) );
                        iEmissionHeightWeight    += iEmissionHeightWeightTemp;
                        iEmissionHeight          += iEmissionHeightTemp * iEmissionHeightWeightTemp;
                        iEmissionHeight2         += iEmissionHeightTemp * iEmissionHeightTemp * iEmissionHeightWeightTemp;
                        iNEM_pairs++;
                        if( nTPair < 1000 )
                        {
                            fEmissionHeightT[nTPair] = iEmissionHeightTemp;
                        }
                    }
                    nTPair++;
                }
            }
        }
    }
    // return mean correction factor (mean of values from all telescope pairs)
    fNTelPairs = nTPair;
    if( iNEM_pairs > 0. && iEmissionHeightWeight > 0. )
    {
        // calculate mean emission height
        fEmissionHeight = iEmissionHeight / iEmissionHeightWeight;
        iEmissionHeight2 /= iEmissionHeightWeight;
        if( iNEM_pairs > 1. )
        {
            fEmissionHeightChi2 = sqrt( 1. / ( iNEM_pairs - 1. ) * ( iEmissionHeight2 - fEmissionHeight * fEmissionHeight ) );
        }
        else
        {
            fEmissionHeightChi2 = 0.;
        }
    }
    else
    {
        fEmissionHeight = 0.;
    }
    
    return fEmissionHeight;
}

void VEmissionHeightCalculator::setTelescopePositions( vector< float > x, vector< float > y, vector< float > z )
{
    fNTel = x.size();
    if( x.size() != y.size()
            || x.size() != z.size() )
    {
        return;
    }
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fTelX.push_back( x[i] );
        fTelY.push_back( y[i] );
        fTelZ.push_back( z[i] );
        if( fDebug )
        {
            cout << "VEmissionHeightCalculator::setTelescopePositions " << i + 1 << "\t" << fTelX.back() << "\t" << fTelY.back() << "\t" << fTelZ.back() << endl;
        }
    }
}


void VEmissionHeightCalculator::setTelescopePositions( unsigned int ntel, double* x, double* y, double* z )
{
    fNTel = ntel;
    for( unsigned int i = 0; i < ntel; i++ )
    {
        fTelX.push_back( x[i] );
        fTelY.push_back( y[i] );
        fTelZ.push_back( z[i] );
        if( fDebug )
        {
            cout << "VEmissionHeightCalculator::setTelescopePositions " << i + 1 << "\t" << fTelX.back() << "\t" << fTelY.back() << "\t" << fTelZ.back() << endl;
        }
    }
}


/*!
    get distance between telescopes in shower coordinates
*/
double VEmissionHeightCalculator::getTelescopeDistanceSC( unsigned int iTel1, unsigned int iTel2, double az, double z )
{
    if( iTel1 >= fTelX.size() || iTel2 >= fTelX.size() )
    {
        cout << "VEmissionHeightCalculator::getTelescopeDistanceSC error: telescope identifier out of range: " << fTelX.size() << "\t" << iTel1 << "\t" << iTel2 << endl;
        return -999.;
    }
    
    double t1[3], t2[3];
    
    az /= TMath::RadToDeg();
    z  /= TMath::RadToDeg();
    
    t1[0] = fTelX[iTel1];
    t1[1] = fTelY[iTel1];
    t1[2] = fTelZ[iTel1];
    
    t2[0] = fTelX[iTel2];
    t2[1] = fTelY[iTel2];
    t2[2] = fTelZ[iTel2];
    
    return VUtilities::line_point_distance( t1[0], t1[1], t1[2], 90. - z * TMath::RadToDeg(), az * TMath::RadToDeg(), t2[0], t2[1], t2[2] );
}


/*!
    calculate distance between two images in the camera
*/
double VEmissionHeightCalculator::imageDistance( double c1x, double c2x, double c1y, double c2y )
{
    return sqrt( ( c1x - c2x ) * ( c1x - c2x ) + ( c1y - c2y ) * ( c1y - c2y ) );
}

