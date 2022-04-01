/* \file VSimpleStereoReconstructor.cpp
 * \brief a simple core and direction reconstruction
 *
 * should not be used for any VERITAS data analysis
 */


#include "VSimpleStereoReconstructor.h"

VSimpleStereoReconstructor::VSimpleStereoReconstructor()
{
    initialize();
    fTelElevation = 0.;
    fTelAzimuth   = 0.;
    reset();
}


void VSimpleStereoReconstructor::initialize( unsigned int iNImages_min,
        float iAxesAngles_min )
{
    fNImages_min    = iNImages_min;
    // require at least stereo events
    if( fNImages_min < 2 )
    {
        fNImages_min = 2;
    }
    fAxesAngles_min = iAxesAngles_min;
}

void VSimpleStereoReconstructor::reset()
{
    
    fiangdiff = 0.;
    fShower_Xoffset = -9999.;
    fShower_Yoffset = -9999.;
    fShower_stdS = -9999.;
    fShower_Chi2 = -9999.;
    fShower_Ze = -9999.;
    fShower_Az = -9999.;
    fShower_DispDiff = -9999.;
    
    fShower_Xcore = -9999.;
    fShower_Ycore = -9999.;
    fShower_stdP = -9999.;
    fmean_iangdiff = 0.;
}


/*!
      reconstruction of shower direction and core

      Hofmann et al 1999, Method 1 (HEGRA method)

      shower direction by intersection of image axes
      shower core by intersection of lines connecting reconstruced shower
      direction and image centroids

      corresponds to rcs_method4 in VArrayAnalyzer

      should be used for MC only

*/
bool VSimpleStereoReconstructor::reconstruct_direction_and_core( unsigned int i_ntel,
        double iArrayElevation,
        double iArrayAzimuth,
        double* iTelX,
        double* iTelY,
        double* iTelZ,
        double* img_size,
        double* img_cen_x,
        double* img_cen_y,
        double* img_cosphi,
        double* img_sinphi,
        double* img_width,
        double* img_length,
        double* img_weight )
{
    // telescope pointings
    fTelElevation = iArrayElevation;
    fTelAzimuth   = iArrayAzimuth;

    // make sure that all data arrays exist
    if( !img_size || !img_cen_x || !img_cen_y
            || !img_cosphi || !img_sinphi
            || !img_width || !img_length
            || !img_weight )
    {
        reset();
        return false;
    }
    
    float xs = 0.;
    float ys = 0.;
    
    // fill data vectors for direction reconstruction
    // (vectors are refilled for core reconstruction)
    vector< float > m;
    vector< float > x;
    vector< float > y;
    vector< float > w;
    vector< float > l;
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        if( img_size[i] > 0. )
        {
            w.push_back( img_size[i] * img_weight[i] );
            x.push_back( img_cen_x[i] );
            y.push_back( img_cen_y[i] );
            // in VArrayAnalyzer, we do a recalculatePhi. Is this needed (for LL)?
            // (not needed, but there will be a very small (<1.e-5) number of showers
            // with different phi values (missing accuracy in conversion from float
            // to double)
            if( img_cosphi[i] != 0. )
            {
                m.push_back( img_sinphi[i] / img_cosphi[i] );
            }
            else
            {
                m.push_back( 1.e9 );
            }
            if( img_length[i] > 0. )
            {
                l.push_back( img_width[i] / img_length[i] );
            }
            else
            {
                l.push_back( 1. );
            }
        }
    }
    // are there enough images the run an array analysis
    if( w.size() < fNImages_min )
    {
        reset();
        return false;
    }
    
    // don't do anything if angle between image axis is too small (for 2 images only)
    if( w.size() == 2 )
    {
        fiangdiff = -1.*fabs( atan( m[0] ) - atan( m[1] ) ) * TMath::RadToDeg();
    }
    else
    {
        fiangdiff = 0.;
    }
    
    ///////////////////////////////
    // direction reconstruction
    ////////////////////////////////////////////////
    // Hofmann et al 1999, Method 1 (HEGRA method)
    // (modified weights)
    
    float itotweight = 0.;
    float iweight = 1.;
    float ixs = 0.;
    float iys = 0.;
    float iangdiff = 0.;
    float b1 = 0.;
    float b2 = 0.;
    vector< float > v_xs;
    vector< float > v_ys;
    fmean_iangdiff = 0.;
    float fmean_iangdiffN = 0.;
    
    for( unsigned int ii = 0; ii < m.size(); ii++ )
    {
        for( unsigned int jj = 1; jj < m.size(); jj++ )
        {
            if( ii >= jj )
            {
                continue;
            }
            
            // check minimum angle between image lines; ignore if too small
            iangdiff = fabs( atan( m[jj] ) - atan( m[ii] ) );
            if( iangdiff < fAxesAngles_min * TMath::DegToRad() ||
                    fabs( 180. * TMath::DegToRad() - iangdiff ) < fAxesAngles_min * TMath::DegToRad() )
            {
                continue;
            }
            // mean angle between images
            if( iangdiff < 90. * TMath::DegToRad() )
            {
                fmean_iangdiff += iangdiff * TMath::RadToDeg();
            }
            else
            {
                fmean_iangdiff += ( 180. - iangdiff * TMath::RadToDeg() );
            }
            fmean_iangdiffN++;
            
            // weight is sin of angle between image lines
            iangdiff = fabs( sin( fabs( atan( m[jj] ) - atan( m[ii] ) ) ) );
            
            b1 = y[ii] - m[ii] * x[ii];
            b2 = y[jj] - m[jj] * x[jj];
            
            // line intersection
            if( m[ii] != m[jj] )
            {
                xs = ( b2 - b1 )  / ( m[ii] - m[jj] );
            }
            else
            {
                xs = 0.;
            }
            ys = m[ii] * xs + b1;
            
            iweight  = 1. / ( 1. / w[ii] + 1. / w[jj] ); // weight 1: size of images
            iweight *= ( 1. - l[ii] ) * ( 1. - l[jj] ); // weight 2: elongation of images (width/length)
            iweight *= iangdiff;                      // weight 3: angular differences between the two image axis
            iweight *= iweight;                       // use squared value
            
            ixs += xs * iweight;
            iys += ys * iweight;
            itotweight += iweight;
            
            v_xs.push_back( xs );
            v_ys.push_back( ys );
        }
    }
    // average difference between image pairs
    if( fmean_iangdiffN > 0. )
    {
        fmean_iangdiff /= fmean_iangdiffN;
    }
    else
    {
        fmean_iangdiff = 0.;
    }
    if( w.size() > 2 )
    {
        fiangdiff = fmean_iangdiff;
    }
    // check validity of weight
    if( itotweight > 0. )
    {
        ixs /= itotweight;
        iys /= itotweight;
        fShower_Xoffset = ixs;
        fShower_Yoffset = iys;
        // calculate dispdiff
        // (this is not exactly dispdiff, but
        //  an equivalent measure comparable to dispdiff)
        fShower_DispDiff = 0.;
        float z = 0;
        for( unsigned n = 0; n < v_xs.size(); n++ )
        {
            for( unsigned m = 0; m < v_xs.size(); m++ )
            {
                if( n > m )
                {
                    fShower_DispDiff += ( v_xs[n] - v_xs[m] ) * ( v_xs[n] - v_xs[m] );
                    fShower_DispDiff += ( v_ys[n] - v_ys[m] ) * ( v_ys[n] - v_ys[m] );
                    z++;
                }
            }
        }
        if( z > 0. )
        {
            fShower_DispDiff /= z;
        }
    }
    else
    {
        fShower_Xoffset = -99999.;
        fShower_Yoffset = -99999.;
        fShower_DispDiff = -999999.;
    }
    
    // fill correct shower direction
    // do not continue with core reconstruction in case there is no valid
    // direction reconstruction
    // (y sign flip!)
    if( !fillShowerDirection( fShower_Xoffset, -1.*fShower_Yoffset ) )
    {
        return false;
    }
    
    ////////////////////////////////////////////////
    // core reconstruction
    ////////////////////////////////////////////////
    
    // calculated telescope positions in shower coordinates
    float i_xcos = sin( ( 90. - fTelElevation ) / TMath::RadToDeg() ) * sin( ( fTelAzimuth - 180. ) / TMath::RadToDeg() );
    float i_ycos = sin( ( 90. - fTelElevation ) / TMath::RadToDeg() ) * cos( ( fTelAzimuth - 180. ) / TMath::RadToDeg() );
    float i_xrot, i_yrot, i_zrot = 0.;
    
    float ximp = 0.;
    float yimp = 0.;
    float stdp = 0.;
    
    float i_cenx = 0.;
    float i_ceny = 0.;
    
    x.clear();
    y.clear();
    m.clear();
    w.clear();
    
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        if( img_size[i] > 0. && img_length[i] > 0. )
        {
            // telescope coordinates
            // shower coordinates (telecope pointing)
            tel_impact( i_xcos, i_ycos, iTelX[i], iTelY[i], iTelZ[i], &i_xrot, &i_yrot, &i_zrot, false );
            // shower coordinates (shower direction)
            // x.push_back( iTelX[i] - ixs / TMath::RadToDeg() * iTelZ[i] );
            // y.push_back( iTelY[i] - iys / TMath::RadToDeg() * iTelZ[i] );
            x.push_back( i_xrot - ixs / TMath::RadToDeg() * i_zrot );
            y.push_back( i_yrot - iys / TMath::RadToDeg() * i_zrot );
            
            // gradient of image
            i_cenx = img_cen_x[i] - ixs;
            i_ceny = img_cen_y[i] - iys;
            if( i_cenx != 0. )
            {
                m.push_back( -1. * i_ceny / i_cenx );
            }
            else
            {
                m.push_back( 1.e9 );
            }
            // image weight
            iweight = img_size[i];
            iweight *= ( 1. - img_width[i] / img_length[i] );
            w.push_back( iweight * iweight );
        }
    }
    
    // Now call perpendicular_distance for the fit, returning ximp and yimp
    rcs_perpendicular_fit( x, y, w, m, ( int )w.size(), &ximp, &yimp, &stdp );
    
    // return to ground coordinates
    fillShowerCore( ximp, yimp );
    fShower_stdP = stdp;
    
    fShower_Chi2 = 0.;
    
    return true;
}

/*
 * check validity outcome of direction reconstruction
 *
 */
bool VSimpleStereoReconstructor::fillShowerDirection( float xoff, float yoff )
{
    if( TMath::IsNaN( yoff ) || TMath::IsNaN( yoff )
            || xoff < -9998. || yoff < -9998. || yoff > 9999.5 )
    {
        reset();
        return false;
    }
    fShower_Xoffset = xoff;
    fShower_Yoffset = yoff;
    
    // ze / az
    double ze = 0.;
    double az = 0.;
    VAstronometry::vlaDtp2s( -1.* fShower_Xoffset*TMath::DegToRad(), 
                                  fShower_Yoffset*TMath::DegToRad(),
                                  fTelAzimuth * TMath::DegToRad(),
                                  fTelElevation * TMath::DegToRad(),
                                   &az, &ze );
    az *= TMath::RadToDeg();
    ze = 90. - ze * TMath::RadToDeg();
            
    if( TMath::IsNaN( ze ) )
    {
        fShower_Ze = -99999.;
    }
    fShower_Ze = ze;
    fShower_Az = VAstronometry::vlaDranrm( az * TMath::DegToRad() ) * TMath::RadToDeg();

    return true;
}

/*
 * calculate shower core in ground coordinates and
 * check validity of core reconstruction results
 */
bool VSimpleStereoReconstructor::fillShowerCore( float ximp, float yimp )
{
    // check validity
    if( !isnormal( ximp ) || !isnormal( yimp ) )
    {
        ximp = -99999.;
        yimp = -99999.;
        return false;
    }
    // reconstructed shower core in ground coordinates
    float i_xcos = 0.;
    float i_ycos = 0.;
    float zimp = 0.;
    float igz = 0.;
    // calculate z in shower coordinates (for z=0 in ground coordinates)
    if( fShower_Ze != 0. )
    {
        zimp = yimp / tan( ( 90. - fShower_Ze ) * TMath::DegToRad() );
    }
    // calculate direction cosinii
    // taking telescope plane as reference plane.
    i_xcos = sin( ( 90. - fTelElevation ) * TMath::DegToRad() )
             * sin( ( fTelAzimuth - 180. ) * TMath::DegToRad() );
    if( fabs( i_xcos ) < 1.e-7 )
    {
        i_xcos = 0.;
    }
    i_ycos = sin( ( 90. - fTelElevation ) * TMath::DegToRad() )
             * cos( ( fTelAzimuth - 180. ) * TMath::DegToRad() );
    if( fabs( i_ycos ) < 1.e-7 )
    {
        i_ycos = 0.;
    }
    tel_impact( i_xcos, i_ycos, ximp, yimp, zimp, &fShower_Xcore, &fShower_Ycore, &igz, true );
    if( isinf( fShower_Xcore ) || isinf( fShower_Ycore )
            || TMath::IsNaN( fShower_Xcore ) || TMath::IsNaN( fShower_Ycore ) )
    {
        fShower_Xcore = -99999;
        fShower_Ycore = -99999;
    }
    return true;
}
