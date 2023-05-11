/*! \class VStereoMaps
 *  fill sky maps and calculate on and off events
 *
 *
 *  Random generator is used for:
 *  =============================
 *
 *  RING BACKGROUND MODEL
 *
 *   - position of test source in bin i,j in sky map is
 *     taken randomly from this bin
 *    - position of test source in bin i,j in sky map is
 *      taken randomly from this bin
 *
 *  REFLECTED REGION MODEL
 *
 *    - definition of off source region: if more than fRE_nMaxoffsource off source
 *      region would fit into the sky plot, choose fRE_nMaxoffsource regions
 *      randomly from the available regions (this is optional)
 *    - position of test source in bin i,j in sky map is
 *      taken randomly from this bin
 *
 *
 */

#include "VStereoMaps.h"

VStereoMaps::VStereoMaps( bool iuc, int iRandomSeed, bool iTMPL_RE_nMaxoffsource )
{
    fData = 0;
    
    bUncorrelatedSkyMaps = iuc;
    fNoSkyPlots = false;
    
    fRM_file = 0;
    fInitRun = 0;
    
    fRandom = new TRandom3( iRandomSeed );
    
    fAcceptance = 0;
    
    fTMPL_RE_nMaxoffsource = iTMPL_RE_nMaxoffsource;
    
    hmap_stereo = 0;
    hmap_alpha = 0;
    hmap_ratio = 0;
    
    // diagnostic histograms for reflected region analysis
    hAuxHisList = 0;
    
    fTheta2Cut_Max = 0.;
    
    fTheta2_length = 0;
    fTheta2.assign( 100, 0. );
    fTheta2_weight.assign( 100, 0. );
    fTheta2_weightREonly.assign( 100, 0. );
    fTheta2_All.assign( 25, 99.0 );
    fTargetShiftNorth = 0;
    fTargetShiftWest = 0.;
    
    hAux_theta2On = 0;
    hAux_theta2Off = 0;
    hAux_theta2Ratio = 0;
}

VStereoMaps::~VStereoMaps()
{
    if( fRandom )
    {
        delete fRandom;
    }
    if( hAux_theta2On )
    {
        delete hAux_theta2On;
    }
    if( hAux_theta2Off )
    {
        delete hAux_theta2Off;
    }
    if( hAux_theta2Ratio )
    {
        delete hAux_theta2Ratio;
    }
    if( hRE_regions )
    {
        delete hRE_regions;
    }
    if( fAcceptance )
    {
        delete fAcceptance;
    }
}

void VStereoMaps::setTargetShift( double iW, double iN )
{
    fTargetShiftWest = iW;
    fTargetShiftNorth = iN;
}

/*
    copy relevant parameters for exclusion regions from the
    anasum runparameters into vectors
*/
void VStereoMaps::setRegionToExclude( vector< VListOfExclusionRegions* > iF )
{
    fListOfExclusionRegions = iF;
}

void VStereoMaps::setHistograms( TH2D* i_map_stereo, TH2D* i_map_alpha, TH1D* i_map_ratio )
{

    hmap_stereo = i_map_stereo;
    hmap_alpha = i_map_alpha;
    hmap_ratio = i_map_ratio;
}


void VStereoMaps::setRunList( VAnaSumRunParameterDataClass iL )
{
    fRunList = iL;
    fRunList.fWobbleNorthMod *= -1.;
}

/*

    fill an ON or OFF event into a sky map

*/
bool VStereoMaps::fill( bool is_on, double x_sky, double y_sky, double theta2Cut_max, int irun, bool i_isGamma, double& i_theta2 )
{
    // maximum allow theta2 cut is the size of the source region
    fTheta2Cut_Max = theta2Cut_max;
    if( fTheta2Cut_Max > fRunList.fSourceRadius )
    {
        fTheta2Cut_Max = fRunList.fSourceRadius;
    }
    
    if( is_on )
    {
        return fillOn( x_sky, y_sky, irun, i_isGamma, i_theta2 );
    }
    else
    {
        return fillOff( x_sky, y_sky, irun, i_isGamma, i_theta2 );
    }
    
    return false;
}


/*

    fill ON event into a sky map

*/
bool VStereoMaps::fillOn( double x, double y, int irun, bool i_isGamma, double& i_theta2 )
{
    double i_weight = 1.;
    double i_MeanSignalBackgroundAreaRatio = 1.;
    
    // calculate theta2: - relative wobble position + shift in derotated sky coordinates
    double theta2  = ( x - fRunList.fWobbleWestMod  + fTargetShiftWest ) * ( x - fRunList.fWobbleWestMod  + fTargetShiftWest );
    theta2 += ( y - fRunList.fWobbleNorthMod + fTargetShiftNorth ) * ( y - fRunList.fWobbleNorthMod + fTargetShiftNorth );
    
    // calculate on to off area ratio due to energy dependent theta2 cut
    // (this is additional to the alpha values)
    if( fRunList.fBackgroundModel == eREFLECTEDREGION )
    {
        if( fRunList.fSourceRadius > 0. )
        {
            i_MeanSignalBackgroundAreaRatio = fTheta2Cut_Max / fRunList.fSourceRadius;
        }
    }
    else if( fRunList.fBackgroundModel == eRINGMODEL )
    {
        if( fRunList.fSourceRadius > 0. && fTheta2Cut_Max > 0. )
        {
            double iR1 = fRunList.fRM_RingRadius * fRunList.fRM_RingWidth / fRunList.fSourceRadius;
            double iR2 = fRunList.fRM_RingRadius * fRunList.fRM_RingWidth / fTheta2Cut_Max;
            if( iR2 > 0. )
            {
                i_MeanSignalBackgroundAreaRatio = iR1 / iR2;
            }
        }
    }
    
    
    if( i_isGamma )
    {
    
        //////////////////////////////////
        // initialisations
        if( irun != fInitRun )
        {
            cout << "\t filling ON events...";
            if( !bUncorrelatedSkyMaps )
            {
                cout << " (correlated maps)" << endl;
            }
            else
            {
                cout << " (uncorrelated maps)" << endl;
            }
            fInitRun = irun;
            
            //////////////////////////////////
            // REFLECTED REGION MODEL
            //////////////////////////////////
            if( fRunList.fBackgroundModel == eREFLECTEDREGION )
            {
                defineAcceptance();
            }
            ////////////////////////////////////////
            // RINGMODEL
            ////////////////////////////////////////
            else if( fRunList.fBackgroundModel == eRINGMODEL )
            {
                // at the beginning of each run, get acceptance curve
                if( !initialize_RingBackgroundModel( true ) )
                {
                    return false;
                }
            }                                     // irun != fInitRun
        }
        
        // this cuts events away further than the maximum distance
        if( fAcceptance && fAcceptance->isExcludedfromSource( x, y ) )
        {
            return false;
        }
        // fill 2D stereo Maps for all models (ON map)
        makeTwoDStereo_BoxSmooth( x - fRunList.fWobbleWestMod, y - fRunList.fWobbleNorthMod, i_weight, sqrt( fTheta2Cut_Max ), i_MeanSignalBackgroundAreaRatio );
    }
    
    // direction cut
    // (theta2 cut might be energy dependent)
    if( theta2 < fTheta2Cut_Max )
    {
        i_theta2 = theta2;
        return true;
    }
    
    return false;
}

/*

   fill OFF events into a sky map

*/
bool VStereoMaps::fillOff( double x, double y, int irun, bool i_isGamma, double& i_theta2 )
{
    // printout at beginning of each run
    if( irun != fInitRun )
    {
        cout << "\t filling OFF events...";
        if( !bUncorrelatedSkyMaps )
        {
            cout << " (correlated maps)" << endl;
        }
        else
        {
            cout << " (uncorrelated maps)" << endl;
        }
        if( fRunList.fBackgroundModel == eONOFF )
        {
            fInitRun = irun;
        }
    }
    ///////////////////////////////////////
    // ON/OFF analysis
    // (same method as filling ON runs)
    
    if( fRunList.fBackgroundModel == eONOFF )
    {
        return fillOn( x, y, irun, i_isGamma, i_theta2 );
    }
    
    ///////////////////////////////////////
    
    ///////////////////////////////////////
    // RINGMODEL
    
    else if( fRunList.fBackgroundModel == eRINGMODEL )
    {
        return fill_RingBackgroundModel( x, y, irun, i_isGamma );
    }
    
    ///////////////////////////////////////
    
    ///////////////////////////////////////
    //  REFLECTED REGION MODEL
    
    else if( fRunList.fBackgroundModel == eREFLECTEDREGION )
    {
        return fill_ReflectedRegionModel( x, y, irun, i_isGamma, i_theta2 );
    }
    
    ///////////////////////////////////////
    
    ///////////////////////////////////////
    // TEMPLATE MODEL
    // (work in progress, not completely implemeted yet)
    else if( fRunList.fBackgroundModel == eTEMPLATE )
    {
        // apply new shape cuts
        i_isGamma = true;
        if( fData->MSCW < fRunList.fTE_mscw_min || fData->MSCW > fRunList.fTE_mscw_max )
        {
            i_isGamma = false;
        }
        if( fData->MSCL < fRunList.fTE_mscl_min || fData->MSCL > fRunList.fTE_mscl_max )
        {
            i_isGamma = false;
        }
        return fillOn( x, y, irun, i_isGamma, i_theta2 );
    }
    
    return false;
}

/*

    fill (in general ON) event into a sky map

    simple box smoothing

    i_MeanSignalBackgroundAreaRatio = ratio of signal to background area (e.g. for energy dependent theta2 cut)

*/
void VStereoMaps::makeTwoDStereo_BoxSmooth( double i_xderot, double i_yderot, double i_weight, double thetaCutMax, double i_MeanSignalBackgroundAreaRatio )
{
    // fill uncorrelated maps
    // (only events which fall directly in this bin are counted)
    if( bUncorrelatedSkyMaps )
    {
        hmap_stereo->Fill( i_xderot, i_yderot );
        hmap_alpha->Fill( i_xderot, i_yderot, i_weight );
        if( hmap_ratio )
        {
            hmap_ratio->Fill( i_MeanSignalBackgroundAreaRatio );
        }
        return;
    }
    ///////////////////////////////////
    // fill correlated maps
    
    // Constructs a 2D skymap of reconstructed source location on the camera plane
    // smoothed by a box (i.e accept all events within a given radius)
    // thetaCutMax is the radius of the smoothing box
    // (and at the same time sqrt(theta2Max)
    
    // bin numbers for current event
    int i_x = hmap_stereo->GetXaxis()->FindBin( i_xderot );
    int i_y = hmap_stereo->GetYaxis()->FindBin( i_yderot );
    
    // get upper and lower limit for loop (add 2 bins to be on the save side)
    int fn_r0X = int( thetaCutMax / hmap_stereo->GetXaxis()->GetBinWidth( 2 ) ) + 2;
    int fn_r0Y = int( thetaCutMax / hmap_stereo->GetYaxis()->GetBinWidth( 2 ) ) + 2;
    
    int ix_start = i_x - fn_r0X;
    if( ix_start < 0 )
    {
        ix_start = 0;
    }
    int ix_stopp = i_x + fn_r0X;
    if( ix_stopp > hmap_stereo->GetNbinsX() )
    {
        ix_stopp = hmap_stereo->GetNbinsX();
    }
    
    int iy_start = i_y - fn_r0Y;
    if( iy_start < 0 )
    {
        iy_start = 0;
    }
    int iy_stopp = i_y + fn_r0Y;
    if( iy_stopp > hmap_stereo->GetNbinsY() )
    {
        iy_stopp = hmap_stereo->GetNbinsY();
    }
    
    double i_xbin = 0.;
    double i_ybin = 0.;
    double i_r = 0.;
    
    for( int i = ix_start; i < ix_stopp; i++ )
    {
        for( int j = iy_start; j < iy_stopp; j++ )
        {
            i_xbin = hmap_stereo->GetXaxis()->GetBinCenter( i + 1 );
            i_ybin = hmap_stereo->GetYaxis()->GetBinCenter( j + 1 );
            // test if this position is inside maximum accepted distance from camera center
            if( sqrt( ( i_xbin + fRunList.fWobbleWestMod ) * ( i_xbin + fRunList.fWobbleWestMod ) +
                      ( i_ybin + fRunList.fWobbleNorthMod ) * ( i_ybin + fRunList.fWobbleNorthMod ) ) > fRunList.fmaxradius )
            {
                continue;
            }
            
            // theta2 cut
            i_r = sqrt( ( i_xderot - i_xbin ) * ( i_xderot - i_xbin ) + ( i_yderot - i_ybin ) * ( i_yderot - i_ybin ) );
            if( i_r <= thetaCutMax )
            {
                hmap_stereo->Fill( i_xbin, i_ybin );
                hmap_alpha->Fill( i_xbin, i_ybin, i_weight );
                if( hmap_ratio )
                {
                    hmap_ratio->Fill( i_MeanSignalBackgroundAreaRatio );
                }
            }
        }
    }
}


/*!
 *
 * return weighting for target bin
 *
 */
void VStereoMaps::finalize( bool iIsOn, double OnOff_Alpha )
{
    //  if there is one run in on/off, assume that for all runs
    
    ///////////////////////////////////////////
    // ONOFF
    if( fRunList.fBackgroundModel == eONOFF )
    {
        if( OnOff_Alpha <= 0 )
        {
            cout << "WARNING: void VStereoMaps::finalize( bool iIsOn, double OnOff_Alpha)" << endl;
            cout << "Open probability space for data is less then or equal to zero" << endl;
            cout << "This is probably unintended and results from an error in the" << endl;
            cout << "RFCutLowerVals/RFCutUpperVals values in the cut file" << endl;
            cout << "(" << iIsOn << ", " << OnOff_Alpha << ")" << endl;
        }
        
        for( int i = 1; i <= hmap_alpha->GetNbinsX(); i++ )
        {
            for( int j = 1; j <= hmap_alpha->GetNbinsY(); j++ )
            {
                hmap_alpha->SetBinContent( i, j, OnOff_Alpha );
            }
        }
    }
    
    ///////////////////////////////////////////
    // RINGMODEL
    
    else if( fRunList.fBackgroundModel == eRINGMODEL )
    {
        RM_getAlpha( iIsOn );
    }
    
    ///////////////////////////////////////////
    // REFLECTED REGION MODEL
    
    else if( fRunList.fBackgroundModel == eREFLECTEDREGION )
    {
        RE_getAlpha( iIsOn );
    }
    
    // remove stuff we don't need anymore
    cleanup();
}


/*!
 *
 *   calculate accepted weighted area of source or background region
 *
 *   following Funk 2005, but adapted to case here
 *
 */
void VStereoMaps::RM_getAlpha( bool iIsOn )
{
    cout << "\t finalizing ring model alpha plots ";
    if( bUncorrelatedSkyMaps )
    {
        cout << "(uncorrelated maps)";
    }
    else
    {
        cout << "(correlated maps)";
    }
    cout << endl;
    
    // check acceptance
    if( !fAcceptance )
    {
        cout << "VStereoMaps::RM_getAlpha error: no acceptance map" << endl;
        return;
    }
    
    double iRingWidth = fRunList.fRM_RingWidth;
    
    // calculate ring radii
    // outer ring radius is fRM_RingRadius+fRM_RingWidth/2.
    // inner ring radius is fRM_RingRadius-fRM_RingWidth/2.
    
    double i_rU = fRunList.fRM_RingRadius + iRingWidth / 2.;
    double i_rL = fRunList.fRM_RingRadius - iRingWidth / 2.;
    
    double i_rS = sqrt( fRunList.fSourceRadius );
    
    // check if background region does not overlap with source region
    if( i_rL < i_rS )
    {
        cout << "VStereoMaps::RM_getAlpha error: background region overlaps with source region ";
        cout << "inner ring radius " << i_rL << " < source radius " << i_rS << endl;
        cout << "Check your RBM parameters! " << endl;
        return;
    }
    
    double x = 0.;
    double x_w = hmap_stereo->GetXaxis()->GetBinWidth( 2 );
    int ix_start = 0;
    int ix_stopp = 0;
    double y = 0.;
    double y_w = hmap_stereo->GetYaxis()->GetBinWidth( 2 );
    int iy_start = 0;
    int iy_stopp = 0;
    double cx = 0.;
    double cy = 0.;
    double cr = 0.;
    // expected number of bins in source region
    // (assume also that x and y binning is the same;
    //  this assumption is not critical)
    double iNB_expected = 1.;
    if( x_w > 0. )
    {
        iNB_expected = i_rS * i_rS * TMath::Pi() / x_w / x_w;
    }
    if( iNB_expected <= 0. )
    {
        iNB_expected = 1.;
    }
    
    double i_acc = 0.;
    
    // all calculations are with camera center at (0,0), but alpha histograms are
    // filled relative to source center at (0,0)
    int i_xoff = TMath::Nint( fRunList.fWobbleWestMod / x_w );
    int j_yoff = TMath::Nint( fRunList.fWobbleNorthMod / y_w );
    
    // loop over all possible source positions,
    // calculate alpha for source or ring positions
    
    // get boundaries for loop over ring area (square which encloses ring)
    // get loop boundaries (the square on the map which encloses the ring)
    
    // lower edge of bins
    double i_xLE = 0.;
    double i_yLE = 0.;
    
    int i_nbinsXstart = 1;
    int i_nbinsX = hmap_alpha->GetNbinsX();
    int i_nbinsXstopp = i_nbinsX;
    int i_nbinsYstart = 1;
    int i_nbinsY = hmap_alpha->GetNbinsY();
    int i_nbinsYstopp = i_nbinsY;
    if( fNoSkyPlots )
    {
        i_nbinsXstart = hmap_stereo->GetXaxis()->FindBin( fRunList.fWobbleWestMod );
        i_nbinsXstopp = hmap_stereo->GetXaxis()->FindBin( fRunList.fWobbleWestMod );
        i_nbinsYstart = hmap_stereo->GetYaxis()->FindBin( fRunList.fWobbleNorthMod );
        i_nbinsYstopp = hmap_stereo->GetYaxis()->FindBin( fRunList.fWobbleNorthMod );
    }
    
    // x,y (or i,j) is the position of the test source
    for( int i = i_nbinsXstart; i <= i_nbinsXstopp; i++ )
    {
        i_xLE = hmap_alpha->GetXaxis()->GetBinLowEdge( i );
        x = fRandom->Uniform( i_xLE, i_xLE + x_w );
        
        // define area which has to be test for this test source (x)
        ix_start = i - ( int )( i_rU / x_w + 1 );
        if( ix_start < 1 )
        {
            ix_start = 1;
        }
        ix_stopp = i + ( int )( i_rU / x_w + 1 );
        if( ix_stopp > i_nbinsX )
        {
            ix_stopp = i_nbinsX;
        }
        
        for( int j = i_nbinsYstart; j <= i_nbinsYstopp; j++ )
        {
            i_yLE = hmap_alpha->GetYaxis()->GetBinLowEdge( j );
            y = fRandom->Uniform( i_yLE, i_yLE + y_w );
            
            // check if events is in fiducial area
            if( sqrt( x * x + y * y ) > fRunList.fmaxradius )
            {
                hmap_alpha->SetBinContent( i - i_xoff, j - j_yoff, 0. );
                continue;
            }
            
            // define area which has to be test for this test source (y)
            iy_start = j - ( int )( i_rU / y_w + 1 );
            if( iy_start < 1 )
            {
                iy_start = 1;
            }
            iy_stopp = j + ( int )( i_rU / y_w + 1 );
            if( iy_stopp > i_nbinsY )
            {
                iy_stopp = i_nbinsY;
            }
            
            // reset acceptance
            i_acc = 0.;
            
            // loop over all possible ring positions
            for( int ci = ix_start; ci <= ix_stopp; ci++ )
            {
                i_xLE = hmap_alpha->GetXaxis()->GetBinLowEdge( ci );
                
                for( int cj = iy_start; cj <= iy_stopp; cj++ )
                {
                    i_yLE = hmap_alpha->GetYaxis()->GetBinLowEdge( cj );
                    
                    // cx, cy are camera coordinates relative to the camera center
                    cx = fRandom->Uniform( i_xLE, i_xLE + x_w );
                    cy = fRandom->Uniform( i_yLE, i_yLE + y_w );
                    
                    cr = ( cx - x ) * ( cx - x ) + ( cy - y ) * ( cy - y );
                    
                    //////////////////////////////////////////////
                    // get sum of acceptance for on region
                    if( iIsOn )
                    {
                        // check is event is in fiducial area
                        if( fAcceptance->isExcludedfromSource( cx, cy ) )
                        {
                            continue;
                        }
                        // check if event is in 'on' region
                        if( bUncorrelatedSkyMaps )
                        {
                            // check if event is in bin
                            if( ci == i && cj == j )
                            {
                                i_acc += fAcceptance->getAcceptance( cx, cy );
                            }
                        }
                        else
                        {
                            // check if event is inside source region
                            if( cr < i_rS * i_rS )
                            {
                                i_acc += fAcceptance->getAcceptance( cx, cy );
                            }
                        }
                    }
                    //////////////////////////////////////////////
                    // get acceptance for OFF regions
                    else
                    {
                        // don't use events in source region
                        if( fAcceptance->isExcludedfromBackground( cx, cy ) )
                        {
                            continue;
                        }
                        // check if event is inside ring
                        if( cr > i_rL * i_rL && cr < i_rU * i_rU )
                        {
                            i_acc += fAcceptance->getAcceptance( cx, cy );
                        }
                    }
                }
            }
            // fill the alpha map (scaled to acceptance 1)
            hmap_alpha->SetBinContent( i - i_xoff, j - j_yoff, i_acc / iNB_expected );
        }
    }
}


/*!

  normalization map for reflected region model


*/
void VStereoMaps::RE_getAlpha( bool iIsOn )
{
    for( int i = 1; i <= hmap_alpha->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= hmap_alpha->GetNbinsY(); j++ )
        {
            if( iIsOn )
            {
                hmap_alpha->SetBinContent( i, j, 1. );
            }
            else
            {
                // correlated maps
                if( !bUncorrelatedSkyMaps && i < ( int )fRE_off.size() && j < ( int )fRE_off[i].size() )
                {
                    hmap_alpha->SetBinContent( hmap_alpha->GetXaxis()->FindBin( hmap_alpha->GetXaxis()->GetBinCenter( i ) - fRunList.fWobbleWestMod ),
                                               hmap_alpha->GetYaxis()->FindBin( hmap_alpha->GetYaxis()->GetBinCenter( j ) - fRunList.fWobbleNorthMod ),
                                               ( double )fRE_off[i][j].noff );
                }
                // uncorrelated maps
                else if( hmap_stereo->GetBinContent( i, j ) > 0. )
                {
                    hmap_alpha->SetBinContent( i, j, hmap_alpha->GetBinContent( i, j ) / hmap_stereo->GetBinContent( i, j ) );
                }
                else
                {
                    hmap_alpha->SetBinContent( i, j, 0. );
                }
            }
        }
    }
}

/*

       main routine for reflected region model

*/
bool VStereoMaps::fill_ReflectedRegionModel( double x, double y, int irun, bool i_isGamma, double& i_theta2 )
{
    // at the beginning of each run, set up off regions
    if( irun != fInitRun )
    {
        fInitRun = irun;
        if( !initialize_ReflectedRegionModel() )
        {
            // intialization failed, return 0
            return false;
        }
        
        int i_nbinsX = hmap_stereo->GetNbinsX();
        int i_nbinsY = hmap_stereo->GetNbinsY();
        
        f_RE_binXW = hmap_stereo->GetXaxis()->GetBinWidth( 2 );
        f_RE_binYW = hmap_stereo->GetYaxis()->GetBinWidth( 2 );
        
        // only loop over areas in map which are inside the fiducal area
        // (assume box which encompasses fiducal area)
        f_RE_xstart = hmap_stereo->GetXaxis()->FindBin( -1. * fRunList.fmaxradius );
        if( f_RE_xstart < 1 )
        {
            f_RE_xstart = 1;
        }
        f_RE_xstopp = hmap_stereo->GetXaxis()->FindBin( fRunList.fmaxradius );
        if( f_RE_xstopp > i_nbinsX )
        {
            f_RE_xstopp = i_nbinsX;
        }
        f_RE_ystart = hmap_stereo->GetYaxis()->FindBin( -1. * fRunList.fmaxradius );
        if( f_RE_ystart < 1 )
        {
            f_RE_ystart = 1;
        }
        f_RE_ystopp = hmap_stereo->GetYaxis()->FindBin( fRunList.fmaxradius );
        if( f_RE_ystopp > i_nbinsY )
        {
            f_RE_ystopp = i_nbinsY;
        }
        
        // normalisation for uncorrelated plots
        // (off regions are still of source region, this is not the bin size)
        // (TTTH2)
        // add here energy dependent norm for theta2 cut
        f_RE_AreaNorm = 1.;
        if( bUncorrelatedSkyMaps && fRunList.fSourceRadius > 0. )
        {
            f_RE_AreaNorm = hmap_stereo->GetXaxis()->GetBinWidth( 2 ) * hmap_stereo->GetYaxis()->GetBinWidth( 2 ) / TMath::Pi() / fRunList.fSourceRadius;
            if( f_RE_AreaNorm > 0. )
            {
                f_RE_AreaNorm = 1. / f_RE_AreaNorm;
            }
        }
        
        // bin of source direction
        f_RE_WW = hmap_stereo->GetXaxis()->FindBin( fRunList.fWobbleWestMod - fTargetShiftWest );
        f_RE_WN = hmap_stereo->GetYaxis()->FindBin( fRunList.fWobbleNorthMod - fTargetShiftNorth );
    }
    
    ///////////////////////////
    // filling
    
    // first check if (x,y) is inside the fiducal area in the camera
    // (fiducal area is defined as distance to center + ringradius < cameraradius)
    //
    double i_evDist = sqrt( x * x + y * y );
    if( i_evDist > fRunList.fmaxradius )
    {
        return false;
    }
    
    // now loop over the whole map to check if events is in one of the off regions
    if( i_isGamma )
    {
        double i_cx = 0.;
        double i_cy = 0.;
        double i_binDist = 0.;
        unsigned int i_nr = 0;
        
        for( int i = f_RE_xstart; i <= f_RE_xstopp; i++ )
        {
            i_cx =  hmap_stereo->GetXaxis()->GetBinCenter( i );
            
            for( int j = f_RE_ystart; j <= f_RE_ystopp; j++ )
            {
                i_cy =  hmap_stereo->GetYaxis()->GetBinCenter( j );
                
                // check if event is in the same ring as this bin (all off regions are in a ring around the camera center)
                i_binDist = sqrt( i_cx * i_cx + i_cy * i_cy );
                
                if( i_evDist > i_binDist + fRE_roffTemp )
                {
                    continue;
                }
                if( i_evDist < i_binDist - fRE_roffTemp )
                {
                    continue;
                }
                
                // loop over all off regions for this bin
                i_nr = fRE_off[i][j].xoff.size();
                if( fRE_off[i][j].noff == 0 )
                {
                    i_nr = 0;
                }
                
                for( unsigned int p = 0; p < i_nr; p++ )
                {
                    // apply theta2 cut in background region
                    double theta2 = ( x - fRE_off[i][j].xoff[p] ) * ( x - fRE_off[i][j].xoff[p] )
                                    + ( y - fRE_off[i][j].yoff[p] ) * ( y - fRE_off[i][j].yoff[p] );
                                    
                    if( theta2 < fRE_off[i][j].roff[p]*fRE_off[i][j].roff[p] )
                    {
                        i_theta2 = theta2;
                        hmap_stereo->Fill( i_cx - fRunList.fWobbleWestMod, i_cy - fRunList.fWobbleNorthMod );
                        hmap_alpha->Fill( i_cx - fRunList.fWobbleWestMod, i_cy - fRunList.fWobbleNorthMod, ( double )fRE_off[i][j].noff * f_RE_AreaNorm );
                    }
                }
            }
        }
    }
    
    bool is_inside = false;
    // return value for source position (weight to calculate MSCW/MSCL etc. plots)
    if( f_RE_WW < ( int )fRE_off.size() && f_RE_WN < ( int )fRE_off[f_RE_WW].size() )
    {
        for( unsigned int p = 0; p < fRE_off[f_RE_WW][f_RE_WN].xoff.size(); p++ )
        {
            if( p < 25 )
            {
                fTheta2_All[p] = ( x - fRE_off[f_RE_WW][f_RE_WN].xoff[p] ) * ( x - fRE_off[f_RE_WW][f_RE_WN].xoff[p] ) +
                                 ( y - fRE_off[f_RE_WW][f_RE_WN].yoff[p] ) * ( y - fRE_off[f_RE_WW][f_RE_WN].yoff[p] );
            }
            
            if( ( x - fRE_off[f_RE_WW][f_RE_WN].xoff[p] ) * ( x - fRE_off[f_RE_WW][f_RE_WN].xoff[p] ) +
                    ( y - fRE_off[f_RE_WW][f_RE_WN].yoff[p] ) * ( y - fRE_off[f_RE_WW][f_RE_WN].yoff[p] ) <
                    fRE_off[f_RE_WW][f_RE_WN].roff[p]*fRE_off[f_RE_WW][f_RE_WN].roff[p] )
            {
                double t2temp = fTheta2_All[0];
                fTheta2_All[0] = ( x - fRE_off[f_RE_WW][f_RE_WN].xoff[p] ) * ( x - fRE_off[f_RE_WW][f_RE_WN].xoff[p] ) +
                                 ( y - fRE_off[f_RE_WW][f_RE_WN].yoff[p] ) * ( y - fRE_off[f_RE_WW][f_RE_WN].yoff[p] );
                i_theta2 = fTheta2_All[0];
                fTheta2_All[p] = t2temp;
                is_inside = true;
            }
            
        }
    }
    
    return is_inside;
}


/*!
 *   calculate number and positions of background regions
 *
 *   store this in 2D vector of sRE_REGIONS
 *
 *   x,y not rotated to source position
 */
bool VStereoMaps::initialize_ReflectedRegionModel()
{
    cout << "\t initialize REFLECTED REGION MODEL for run " << fInitRun;
    if( bUncorrelatedSkyMaps )
    {
        cout << " (uncorrelated plots)";
    }
    cout << endl;
    
    //  control histograms and delete the one from the previous runs
    initialize_ReflectedRegionHistograms();
    
    // reflected region variables
    double x = 0.;
    double y = 0.;
    int x_bin = 0;
    int y_bin = 0;
    double x_wobble = 0.;
    double y_wobble = 0.;
    int x_bin_wobble = 0;
    int y_bin_wobble = 0;
    double r = 0.;                                // size of source region
    // reflected regions (all in wobble shifted coordinates)
    int n_r = 0;
    double x_re[1000];
    double y_re[1000];
    double r_re[1000];
    // exclusion regions
    int n_ex = 0;
    double x_ex[1000];
    double y_ex[1000];
    double r1_ex[1000];
    double r2_ex[1000];
    double ang_ex[1000];
    for( unsigned int i = 0; i < 1000; i++ )
    {
        x_re[i] = 0.;
        y_re[i] = 0.;
        r_re[i] = 0.;
        x_ex[i] = 0.;
        y_ex[i] = 0.;
        r1_ex[i] = 0.;
        r2_ex[i] = 0.;
        ang_ex[i] = 0;
    }
    
    // tree with all reflected regions (only for correlated maps)
    // (this tree is written to the output root file in the directory run_XXX/stereo/debug)
    if( !bUncorrelatedSkyMaps )
    {
        hRE_regions = new TTree( "tRE", "reflected regions" );
        hRE_regions->Branch( "x", &x, "x/D" );
        hRE_regions->Branch( "y", &y, "y/D" );
        hRE_regions->Branch( "x_bin", &x_bin, "x_bin/I" );
        hRE_regions->Branch( "y_bin", &y_bin, "y_bin/I" );
        hRE_regions->Branch( "x_wobble", &x_wobble, "x_wobble/D" );
        hRE_regions->Branch( "y_wobble", &y_wobble, "y_wobble/D" );
        hRE_regions->Branch( "r", &r, "r/D" );
        hRE_regions->Branch( "x_bin_wobble", &x_bin_wobble, "x_bin_wobble/I" );
        hRE_regions->Branch( "y_bin_wobble", &y_bin_wobble, "y_bin_wobble/I" );
        hRE_regions->Branch( "n_r", &n_r, "n_r/I" );
        hRE_regions->Branch( "x_re", x_re, "x_re[n_r]/D" );
        hRE_regions->Branch( "y_re", y_re, "y_re[n_r]/D" );
        hRE_regions->Branch( "r_re", r_re, "r_re[n_r]/D" );
        hRE_regions->Branch( "n_ex", &n_ex, "n_ex/I" );
        hRE_regions->Branch( "x_ex", x_ex, "x_ex[n_ex]/D" );
        hRE_regions->Branch( "y_ex", y_ex, "y_ex[n_ex]/D" );
        hRE_regions->Branch( "r1_ex", r1_ex, "r1_ex[n_ex]/D" );
        hRE_regions->Branch( "r2_ex", r2_ex, "r2_ex[n_ex]/D" );
        hRE_regions->Branch( "ang_ex", ang_ex, "ang_ex[n_ex]/D" );
        if( hAuxHisList )
        {
            hAuxHisList->Add( hRE_regions );
        }
        else
        {
            cout << "VStereoMaps::initialize_ReflectedRegionModel() warning: aux list does not exist"  << endl;
        }
    }
    
    n_ex = ( int )fListOfExclusionRegions.size();
    // (don't expect more than a thousand exclusion regions)
    if( n_ex >= 1000 )
    {
        n_ex = 1000;
    }
    // copy exclusion region positions to tree variables
    for( unsigned int i = 0; i < fListOfExclusionRegions.size(); i++ )
    {
        if( fListOfExclusionRegions[i] )
        {
            x_ex[i]   = fListOfExclusionRegions[i]->fExcludeFromBackground_West;
            y_ex[i]   = fListOfExclusionRegions[i]->fExcludeFromBackground_North;
            r1_ex[i]  = fListOfExclusionRegions[i]->fExcludeFromBackground_Radius1;
            r2_ex[i]  = fListOfExclusionRegions[i]->fExcludeFromBackground_Radius2;
            ang_ex[i] = fListOfExclusionRegions[i]->fExcludeFromBackground_RotAngle;
        }
    }
    
    // off region parameters
    sRE_REGIONS i_off;
    i_off.noff = 0;
    
    // empty vector
    vector< double > i_Dempty;
    
    // source extension (equal to radius of off regions)
    fRE_roffTemp = sqrt( fRunList.fSourceRadius );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // calculate off regions (following Zufelde, 2005, p.28)
    //
    //  for all bins in stereo maps
    //
    //    number of off regions depends on distance of bin to camera center (how many areas fit in?)
    //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // set number of bins to be covered
    int i_nbinsX_min = 0;
    int i_nbinsX = hmap_stereo->GetNbinsX();
    int i_nbinsY_min = 0;
    int i_nbinsY = hmap_stereo->GetNbinsY();
    // check if full sky analysis is required
    if( fNoSkyPlots )
    {
        int x_bin_wobble = hmap_stereo->GetXaxis()->FindBin( fRunList.fWobbleWestMod );
        i_nbinsX_min = x_bin_wobble;
        i_nbinsX     = x_bin_wobble;
        
        int y_bin_wobble = hmap_stereo->GetYaxis()->FindBin( fRunList.fWobbleNorthMod );
        i_nbinsY_min = y_bin_wobble;
        i_nbinsY     = y_bin_wobble;
    }
    
    // set up 2D vector of off source parameters ([n_x][n_y])
    fRE_off.clear();
    vector< sRE_REGIONS > i_x_off;
    for( int j = 0; j <= i_nbinsY; j++ )
    {
        i_x_off.push_back( i_off );
        i_x_off.back().noff = 0;
    }
    for( int i = 0; i <= i_nbinsX; i++ )
    {
        fRE_off.push_back( i_x_off );
    }
    
    // distance of bin to camera center
    double ids = 0.;
    
    // vectors with off source positions
    vector< double > r_off;
    vector< double > x_off;
    vector< double > y_off;
    // vectors with off source positions (temporary)
    vector< double > r_offTemp;
    vector< double > x_offTemp;
    vector< double > y_offTemp;
    
    // calculate maximum number of off source regions possible
    // (maximum number of off source regions at the edge of the camera)
    int n_max_RE = 0;
    if( fRunList.fmaxradius > 0. )
    {
        n_max_RE = ( int )( TMath::Pi() / ( asin( fRE_roffTemp / fRunList.fmaxradius ) ) );
    }
    else
    {
        cout << "VStereoMaps::initialize_ReflectedRegionModel() warning: max camera radius is zero: ";
        cout << fRunList.fmaxradius << endl;
        return false;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // loop over all bins in the stereo maps
    //
    // coordinate system: relative to camera centre
    // (shifted at the very end to sky map centred coordinate system)
    //////////////////////////////////////////////////////////////////////////////////////
    for( int i = i_nbinsX_min; i <= i_nbinsX; i++ )
    {
        x = hmap_stereo->GetXaxis()->GetBinCenter( i );
        if( TMath::Abs( x ) < 1.e-5 )
        {
            x = 0.;
        }
        
        for( int j = i_nbinsY_min; j <= i_nbinsY; j++ )
        {
            y = hmap_stereo->GetYaxis()->GetBinCenter( j );
            if( TMath::Abs( y ) < 1.e-5 )
            {
                y = 0.;
            }
            
            // distance of this bin from camera center
            ids = sqrt( x * x + y * y );
            
            // initialise list of reflected regions for this bin
            n_r = 0;
            r_off.clear();
            x_off.clear();
            y_off.clear();
            fRE_off[i][j].roff = r_off;
            fRE_off[i][j].xoff = x_off;
            fRE_off[i][j].yoff = y_off;
            fRE_off[i][j].noff = 0;
            
            // bin is outside confidence region (distance of bin + off source radius)
            // and bin is too close to center of camera
            if( ids < ( fRunList.fmaxradius ) && ids > fRE_roffTemp && fRE_roffTemp != 0. )
            {
                // angular size of the on region seen from the observation position
                double w = asin( fRE_roffTemp / ids );
                // number of off source regions
                if( w > 0 )
                {
                    n_r = ( int )( TMath::Pi() / w );
                }
                else
                {
                    n_r = 0;
                }
                // set maximum number of off source region possible for this particular distance to camera center
                n_max_RE = n_r;
                
                // test if there are enough off source regions
                if( n_r >= fRunList.fRE_nMinoffsource )
                {
                    // phi angle of on position
                    double phi_0 = atan2( y, x );
                    
                    // off positions
                    double phi_i = 0.;
                    double x_t = 0.;
                    double y_t = 0.;
                    
                    // try to fit at least fRunList.fRE_nMinoffsource into the available space, work by trial and error
                    // (rotate source start of source regions in XX degrees steps)
                    // (step size is a trade off between accuracy and speed)
                    for( double t = 0.; t < 360.; t += 1.0 )
                    {
                        // reset all previously filled vector
                        r_offTemp.clear();
                        x_offTemp.clear();
                        y_offTemp.clear();
                        // loop over all possible off source positions
                        for( int p = 0; p < n_r; p++ )
                        {
                            // get off-source positions
                            phi_i = phi_0 + TMath::Pi() + ( 2 * p + 1 - n_r ) * w;
                            phi_i += ( ( double )( t ) ) / TMath::RadToDeg();
                            x_t = ids * cos( phi_i );
                            y_t = ids * sin( phi_i );
                            
                            // check if off source region is not included in this off position
                            // (require to be at least 2.0*times theta2 circle away)
                            if( ( x_t - x ) * ( x_t - x ) + ( y_t - y ) * ( y_t - y )
                                    < ( 2. + fRunList.fRE_distanceSourceOff ) * ( 2. + fRunList.fRE_distanceSourceOff )*fRunList.fSourceRadius )
                            {
                                continue;
                            }
                            
                            // check if real source region is not included in this off position
                            // check if off source region is not excluded from background
                            bool bExclude = false;
                            for( unsigned int ex = 0; ex < fListOfExclusionRegions.size(); ex++ )
                            {
                                if( !fListOfExclusionRegions[ex] )
                                {
                                    continue;
                                }
                                if( fListOfExclusionRegions[ex]->isInsideExclusionRegion( x_t, y_t, fRE_roffTemp ) )
                                {
                                    bExclude = true;
                                }
                            }
                            if( bExclude )
                            {
                                continue;
                            }
                            
                            // fill an off region
                            r_offTemp.push_back( fRE_roffTemp );
                            x_offTemp.push_back( x_t );
                            y_offTemp.push_back( y_t );
                        }
                        // test if this configuration has more off source regions than a previous one
                        if( x_offTemp.size() > x_off.size() )
                        {
                            r_off = r_offTemp;
                            x_off = x_offTemp;
                            y_off = y_offTemp;
                        }
                        // test if there are enough off source regions
                        if( ( int )x_off.size() >= fRunList.fRE_nMaxoffsource || ( x_off.size() > 10 && ( int )x_off.size() >= n_max_RE - 5 ) )
                        {
                            break;
                        }
                    }
                    n_r = ( int )x_off.size();
                    
                    // check maximum number of sources, if too many, remove some (randomly choosen)
                    // (this might introduce a gradient across the sky maps)
                    if( fTMPL_RE_nMaxoffsource )
                    {
                        if( n_r > fRunList.fRE_nMaxoffsource )
                        {
                            unsigned int remo_iter = 0;
                            while( ( int )r_off.size() > fRunList.fRE_nMaxoffsource )
                            {
                                remo_iter = ( unsigned int )fRandom->Integer( ( int )r_off.size() );
                                r_off.erase( r_off.begin() + remo_iter );
                                x_off.erase( x_off.begin() + remo_iter );
                                y_off.erase( y_off.begin() + remo_iter );
                            }
                        }
                    }
                    // remove those regions which are furthest away
                    // (default)
                    else
                    {
                        double i_dist = 0.;
                        while( ( int )r_off.size() > fRunList.fRE_nMaxoffsource )
                        {
                            double i_max = 0.;
                            unsigned int i_maxtt = 99999;
                            unsigned int tt = 0;
                            for( tt = 0; tt < r_off.size(); tt++ )
                            {
                                i_dist = sqrt( ( x_off[tt] - x ) * ( x_off[tt] - x ) + ( y_off[tt] - y ) * ( y_off[tt] - y ) );
                                if( i_dist > i_max )
                                {
                                    i_max = i_dist;
                                    i_maxtt = tt;
                                }
                            }
                            if( i_maxtt != 99999 && i_maxtt < r_off.size() )
                            {
                                r_off.erase( r_off.begin() + i_maxtt );
                                x_off.erase( x_off.begin() + i_maxtt );
                                y_off.erase( y_off.begin() + i_maxtt );
                            }
                        }
                    }
                    
                    // number of off source regions
                    n_r = ( int )x_off.size();
                    
                    if( hRE_NRegions )
                    {
                        hRE_NRegions->SetBinContent( i, j, n_r );
                    }
                }
            }
            ////////////////////////////////////////
            // now fill off source regions
            
            // not enough off source regions
            if( n_r < fRunList.fRE_nMinoffsource )
            {
                // set n_r to 0 (otherwise problems filling debug tree)
                n_r = 0;
                fRE_off[i][j].roff = i_Dempty;
                fRE_off[i][j].xoff = i_Dempty;
                fRE_off[i][j].yoff = i_Dempty;
                fRE_off[i][j].noff = 0;
            }
            // valid number of off source regions
            else
            {
                fRE_off[i][j].roff = r_off;
                fRE_off[i][j].xoff = x_off;
                fRE_off[i][j].yoff = y_off;
                fRE_off[i][j].noff = n_r;
            }
            
            // fill tree with reflected regions (only for correlated maps)
            if( !bUncorrelatedSkyMaps && hRE_regions )
            {
                x_bin = i;
                y_bin = j;
                x_wobble = x - fRunList.fWobbleWestMod;
                y_wobble = y - fRunList.fWobbleNorthMod;
                x_bin_wobble = hmap_stereo->GetXaxis()->FindBin( x_wobble );
                y_bin_wobble = hmap_stereo->GetYaxis()->FindBin( y_wobble );
                r = fRE_roffTemp;
                
                if( n_r > 0 )
                {
                    // all coordinates are wobble shifted (from unwobbled to wobbled)
                    for( int p = 0; p < n_r; p++ )
                    {
                        if( p > 999 )
                        {
                            continue;
                        }
                        x_re[p] = fRE_off[i][j].xoff[p] - fRunList.fWobbleWestMod;
                        y_re[p] = fRE_off[i][j].yoff[p] - fRunList.fWobbleNorthMod;
                        r_re[p] = fRE_off[i][j].roff[p];
                    }
                    hRE_regions->Fill();
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    cout << "\t\t ....reflected regions initialized" << endl;
    
    return true;
}


void VStereoMaps::initialize_theta2()
{
    // get available phase space for theta2 plots
    
    if( !bUncorrelatedSkyMaps )
    {
        char hname[200];
        sprintf( hname, "hAux_theta2On" );
        hAux_theta2On = new TH1D( hname, "", 200, 0., 2. );
        hAux_theta2On->SetXTitle( "#theta^{2} [deg^{2}]" );
        hAux_theta2On->SetYTitle( "entries" );
        if( hAuxHisList )
        {
            hAuxHisList->Add( hAux_theta2On );
        }
        
        sprintf( hname, "hAux_theta2Off" );
        hAux_theta2Off = new TH1D( hname, "", 200, 0., 2. );
        hAux_theta2Off->SetXTitle( "#theta^{2} [deg^{2}]" );
        hAux_theta2Off->SetYTitle( "entries" );
        hAux_theta2Off->SetLineColor( 2 );
        if( hAuxHisList )
        {
            hAuxHisList->Add( hAux_theta2Off );
        }
        
        sprintf( hname, "hAux_theta2Ratio" );
        hAux_theta2Ratio = new TH1D( hname, "", 200, 0., 2. );
        hAux_theta2Ratio->SetXTitle( "#theta^{2} [deg^{2}]" );
        hAux_theta2Ratio->SetYTitle( "ratio (off/on)" );
        if( hAuxHisList )
        {
            hAuxHisList->Add( hAux_theta2Ratio );
        }
        
        int nxybin = 5000;
        double i_xmin = hmap_stereo->GetXaxis()->GetXmin();
        
        double i_xybinW = fabs( 2 * i_xmin ) / ( double )nxybin;
        
        if( !fAcceptance )
        {
            fAcceptance = new VRadialAcceptance( fRunList.fAcceptanceFile );
            fAcceptance->Set2DAcceptanceMode( fRunList.f2DAcceptanceMode ) ;
        }
        
        double x = 0.;
        double y = 0.;
        
        // loop over all bins in stereo maps
        // (take only every second bin)
        double i_tAcc = 1.;
        for( int i = 0; i < nxybin; i += 2 )
        {
            x = i_xmin + ( double )i * i_xybinW + fRandom->Uniform( -1.*i_xybinW / 2., i_xybinW / 2. );
            
            for( int j = 0; j < nxybin; j += 2 )
            {
                y = i_xmin + ( double )j * i_xybinW + fRandom->Uniform( -1.*i_xybinW / 2., i_xybinW / 2. );
                
                i_tAcc = fAcceptance->getAcceptance( x, y );
                
                calculateTheta2( true, x, y );
                for( unsigned int t = 0; t < getTheta2_length(); t++ )
                {
                    hAux_theta2On->Fill( getTheta2()[t], getTheta2_weigthREonly()[t] * i_tAcc );
                }
                
                calculateTheta2( false, x, y );
                for( unsigned int t = 0; t < getTheta2_length(); t++ )
                {
                    hAux_theta2Off->Fill( getTheta2()[t], getTheta2_weigthREonly()[t] * i_tAcc );
                }
            }
        }
        // calculate the ratio off/on
        hAux_theta2Ratio->Divide( hAux_theta2Off, hAux_theta2On );
    }
    else
    {
        hAux_theta2On = 0;
        hAux_theta2Off = 0;
        hAux_theta2Ratio = 0;
    }
}


double VStereoMaps::phiInt( double iphi )
{
    if( iphi > 2 * TMath::Pi() )
    {
        int nph = ( int )( iphi / ( 2 * TMath::Pi() ) );
        iphi = iphi - ( double )nph * 2 * TMath::Pi();
    }
    return iphi;
}


/*

    fill events for ring background model background estimation ('OFF' events)

*/
bool VStereoMaps::fill_RingBackgroundModel( double x, double y, int irun, bool i_isGamma )
{
    // at the beginning of each run, get a acceptance curve
    if( irun != fInitRun )
    {
        fInitRun = irun;
        if( !initialize_RingBackgroundModel( false ) )
        {
            return false;
        }
    }
    
    // check if (x,y) is inside the fiducal area in the camera
    // check if event is inside an excluded region
    if( fAcceptance->isExcludedfromBackground( x, y ) )
    {
        return false;
    }
    
    // calculate ring radii
    // outer ring radius is fRM_RingRadius+fRM_RingWidth/2.
    // inner ring radius is fRM_RingRadius-fRM_RingWidth/2.
    
    double iRingWidth = fRunList.fRM_RingWidth;
    double i_rU = fRunList.fRM_RingRadius + iRingWidth / 2.;
    double i_rL = fRunList.fRM_RingRadius - iRingWidth / 2.;
    
    // get loop boundaries (the square on the map which encloses the ring)
    
    // coordinates of current shower on the map
    int i_x = hmap_stereo->GetXaxis()->FindBin( x );
    int i_y = hmap_stereo->GetYaxis()->FindBin( y );
    
    int ix_start = i_x - ( int )( i_rU / hmap_stereo->GetXaxis()->GetBinWidth( 2 ) + 1 );
    if( ix_start < 0 )
    {
        ix_start = 0;
    }
    int ix_stopp = i_x + ( int )( i_rU / hmap_stereo->GetXaxis()->GetBinWidth( 2 ) + 1 );
    if( ix_stopp >= hmap_stereo->GetNbinsX() )
    {
        ix_stopp = hmap_stereo->GetNbinsX();
    }
    
    int iy_start = i_y - ( int )( i_rU / hmap_stereo->GetYaxis()->GetBinWidth( 2 ) + 1 );
    if( iy_start < 0 )
    {
        iy_start = 0;
    }
    int iy_stopp = i_y + ( int )( i_rU / hmap_stereo->GetYaxis()->GetBinWidth( 2 ) + 1 );
    if( iy_stopp >= hmap_stereo->GetNbinsY() )
    {
        iy_stopp = hmap_stereo->GetNbinsY();
    }
    
    // bin center of current bin
    double i_cx = 0.;
    double i_cy = 0.;
    double i_cr = 0.;
    
    // now loop over the interesting region on the map
    // test if event is in any of these rings
    if( i_isGamma )
    {
        // loop over box with side length ringradius + ringwidth
        for( int i = ix_start; i <= ix_stopp; i++ )
        {
            for( int j = iy_start; j <= iy_stopp; j++ )
            {
                i_cx = hmap_stereo->GetXaxis()->GetBinCenter( i );
                i_cy = hmap_stereo->GetYaxis()->GetBinCenter( j );
                // check if bin is inside fiducial area
                if( sqrt( i_cx * i_cx + i_cy * i_cy ) > fRunList.fmaxradius )
                {
                    continue;
                }
                
                // check if bin is inside the ring
                i_cr = ( ( i_cx - x ) * ( i_cx - x ) + ( i_cy - y ) * ( i_cy - y ) );
                if( i_cr < i_rU * i_rU && i_cr > i_rL * i_rL )
                {
                    hmap_stereo->Fill( i_cx - fRunList.fWobbleWestMod, i_cy - fRunList.fWobbleNorthMod );
                }
            }
        }
    }
    
    // determine if event is in off region of source region, i.e. inside the ring
    i_cr = ( x - fRunList.fWobbleWestMod ) * ( x - fRunList.fWobbleWestMod ) + ( y - fRunList.fWobbleNorthMod ) * ( y - fRunList.fWobbleNorthMod );
    if( i_cr < i_rU * i_rU && i_cr > i_rL * i_rL )
    {
        return true;
    }
    
    return false;
}


bool VStereoMaps::initialize_RingBackgroundModel( bool iIsOn )
{
    defineAcceptance();
    
    if( !iIsOn )
    {
        initialize_Histograms();
        if( !bUncorrelatedSkyMaps )
        {
            initialize_theta2();
        }
    }
    
    return true;
}


/*!
 *   get theta2 for different model
 *
 *   for ON:  returns vector of length 1 with theta2 to source position
 *   for OFF:
 *            on-off: returns vector of length 1 with theta2 to source position
 *            Reflected regions: return vector of length N_regions with theta2 to center of regions
 *
 */
void VStereoMaps::calculateTheta2( bool isOn, double x, double y )
{
    double theta2 = 0.;
    fTheta2_length = 0;
    double i_phaseSpaceCorrection = 1.;
    
    double iWeight = 1.;
    
    // ON
    if( isOn )
    {
        theta2  = ( x - fRunList.fWobbleWestMod + fTargetShiftWest ) * ( x - fRunList.fWobbleWestMod + fTargetShiftWest );
        theta2 += ( y - fRunList.fWobbleNorthMod + fTargetShiftNorth ) * ( y - fRunList.fWobbleNorthMod + fTargetShiftNorth );
        
        if( fTheta2.size() > 0 && fTheta2_weight.size() > 0 && fTheta2_weightREonly.size() )
        {
            if( fAcceptance )
            {
                iWeight = fAcceptance->getAcceptance( x, y );
            }
            fTheta2[0] = theta2;
            fTheta2_weight[0] = iWeight;
            fTheta2_weightREonly[0] = iWeight;
            fTheta2_length = 1;
        }
        else
        {
            fTheta2_length = 0;
        }
    }
    // theta2 from background
    else
    {
    
        ///////////////////////////////
        // theta2 for ONOFF
        if( fRunList.fBackgroundModel == eONOFF )
        {
            double x_off = fRunList.fWobbleWestMod - fTargetShiftWest;
            double y_off = fRunList.fWobbleNorthMod - fTargetShiftNorth;
            
            theta2  = ( x - x_off ) * ( x - x_off );
            theta2 += ( y - y_off ) * ( y - y_off );
            
            if( fTheta2.size() > 0 && fTheta2_weight.size() > 0 )
            {
                fTheta2[0] = theta2;
                if( fAcceptance )
                {
                    iWeight = fAcceptance->getAcceptance( x, y );
                }
                fTheta2_weight[0] = iWeight;
                fTheta2_length = 1;
            }
            else
            {
                fTheta2_length = 0;
            }
        }
        ///////////////////////////////////////////////
        // RINGMODEL and REFLECTED REGION MODEL
        //
        // theta2 plots are ony filled for illustrations, no significances can be calculated
        // (background is only sampled once)
        //
        // background events come from theta2 region opposite to the potential source region
        //
        else if( fRunList.fBackgroundModel == eRINGMODEL || fRunList.fBackgroundModel == eREFLECTEDREGION )
        {
            theta2  = ( x - fRunList.fWobbleWestMod + fTargetShiftWest ) * ( x - fRunList.fWobbleWestMod + fTargetShiftWest );
            theta2 += ( y - fRunList.fWobbleNorthMod + fTargetShiftNorth ) * ( y - fRunList.fWobbleNorthMod + fTargetShiftNorth );
            // use larger off region distance (wobble offset)
            // (if this region is too small, theta2 plots will be contaminated by missreconstructed gammas for strong sources)
            if( theta2 < fRunList.fWobbleWestMod * fRunList.fWobbleWestMod + fRunList.fWobbleNorthMod * fRunList.fWobbleNorthMod )
            {
                fTheta2_length = 0;
                return;
            }
            theta2  = ( x + fRunList.fWobbleWestMod - fTargetShiftWest ) * ( x + fRunList.fWobbleWestMod - fTargetShiftWest );
            theta2 += ( y + fRunList.fWobbleNorthMod - fTargetShiftNorth ) * ( y + fRunList.fWobbleNorthMod - fTargetShiftNorth );
            fTheta2[0] = theta2;
            i_phaseSpaceCorrection = 1.;
            if( hAux_theta2Ratio )
            {
                i_phaseSpaceCorrection = hAux_theta2Ratio->GetBinContent( hAux_theta2Ratio->FindBin( theta2 ) );
            }
            if( i_phaseSpaceCorrection > 0. )
            {
                i_phaseSpaceCorrection = 1. / i_phaseSpaceCorrection;
            }
            
            if( fAcceptance )
            {
                iWeight = fAcceptance->getAcceptance( x, y );
            }
            fTheta2_weight[0] = i_phaseSpaceCorrection * iWeight;;
            fTheta2_weightREonly[0] = iWeight;
            
            fTheta2_length = 1;
        }
    }
}


/*!
 *
 * diagnostic histograms for reflected region analysis
 */
void VStereoMaps::initialize_ReflectedRegionHistograms()
{
    initialize_Histograms();
    if( !hmap_stereo )
    {
        hRE_NRegions = 0;
        return;
    }
    else
    {
        hRE_NRegions = ( TH2D* )hmap_stereo->Clone( "hmap_stereo" );
        hRE_NRegions->SetTitle( "number of reflected regions" );
        hRE_NRegions->SetStats( 0 );
        if( hAuxHisList )
        {
            hAuxHisList->Add( hRE_NRegions );
        }
    }
    defineAcceptance();
    initialize_theta2();
}


bool VStereoMaps::initialize_Histograms()
{
    if( hAuxHisList )
    {
        hAuxHisList->Delete();
    }
    else
    {
        hAuxHisList = new TList();
    }
    
    return true;
}


/*!
 *   called after each run
 */
void VStereoMaps::cleanup()
{
    // clean sky plots
    // - removing bins which are not completely in fiducial area circle
    // (sky plot and alpha plot)
    double x = 0.;
    double y = 0.;
    double r = 0.;
    double w = hmap_alpha->GetXaxis()->GetBinWidth( 2 );
    
    for( int i = 1; i <= hmap_alpha->GetNbinsX(); i++ )
    {
        x = hmap_alpha->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= hmap_alpha->GetNbinsY(); j++ )
        {
            y = hmap_alpha->GetYaxis()->GetBinCenter( j );
            
            r = ( x + fRunList.fWobbleWestMod ) * ( x + fRunList.fWobbleWestMod ) + ( y + fRunList.fWobbleNorthMod ) * ( y + fRunList.fWobbleNorthMod );
            
            if( r > ( fRunList.fmaxradius - w ) * ( fRunList.fmaxradius - w ) )
            {
                hmap_alpha->SetBinContent( i, j, 0. );
                hmap_stereo->SetBinContent( i, j, 0. );
            }
        }
    }
}


/////////////////////////////////////////////////

bool VStereoMaps::defineAcceptance()
{
    TDirectory* iDir = gDirectory;
    // define acceptance
    if( fAcceptance )
    {
        delete fAcceptance;
    }
    int iRun = -1;
    if( fRunList.fRunWiseRadialAcceptance )
    {
        iRun = fRunList.fRunOff;
    }
    fAcceptance = new VRadialAcceptance( fRunList.fAcceptanceFile, iRun );
    fAcceptance->Set2DAcceptanceMode( fRunList.f2DAcceptanceMode ) ;
    if( !fAcceptance )
    {
        cout << "VStereoMaps::defineAcceptance: ERROR, no acceptance curves found" << endl;
        return false;
    }
    iDir->cd();
    // set exclusion regions
    fAcceptance->setRegionToExcludeAcceptance( fListOfExclusionRegions );
    
    // set the position of the potential gamma-ray source
    // (in de-rotated camera coordinates)
    fAcceptance->setSource( fRunList.fWobbleWestMod + fTargetShiftWest,
                            fRunList.fWobbleNorthMod + fTargetShiftNorth,
                            sqrt( fRunList.fSourceRadius ), fRunList.fmaxradius );
                            
    return true;
}
