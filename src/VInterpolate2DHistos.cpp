/*! \class VInterpolate2DHistos
 *  \brief interpolate empty bins in 2D histograms
 *
 *
 *
 */

#include "VInterpolate2DHistos.h"

VInterpolate2DHistos::VInterpolate2DHistos( int iseed )
{
    fRandom = new TRandom3( iseed );
}

/*
 *   fill empty bins by averaging over surrounding bins
 *
 *   search width is given by iwidth
 *
 *   interpolation can be done several times (maxiter)
 *
 *   **existing bins are not modified**
 *
 */
TH2F* VInterpolate2DHistos::doSimpleInterpolation( TH2F* h, string iname, int iwidth, int maxiter, bool bError,
        TH2F* hNevents, int iMinEvents )
{
    if( !h )
    {
        return 0;
    }
    
    char hname[600];
    sprintf( hname, "%s_%s", h->GetName(), iname.c_str() );
    TH2F* hs = ( TH2F* )h->Clone( hname );
    
    double z = 0.;
    double imean = 0.;
    double ierror = 0.;
    // iterate several times over histogram
    for( int k = 0; k < maxiter; k++ )
    {
        TH2F* htemp = ( TH2F* )hs->Clone();
        for( int i = 1; i <= hs->GetNbinsX(); i++ )
        {
            for( int j = 1; j <= hs->GetNbinsY(); j++ )
            {
                // check number of events
                if( hNevents && hNevents->GetBinContent( i, j ) < iMinEvents )
                {
                    continue;
                }
                else if( hs->GetBinContent( i, j ) <= 0. )
                {
                    continue;
                }
                ///////////////////////////////////////
                // get mean of all bins around this bin
                imean = 0.;
                ierror = 0.;
                z = 0.;
                for( int ii = i - iwidth; ii <= i + iwidth; ii++ )
                {
                    for( int jj = j - iwidth; jj <= j + iwidth; jj++ )
                    {
                        if( ii > 0 && jj > 0 && htemp->GetBinContent( ii, jj ) > 0. )
                        {
                            if( !bError )
                            {
                                imean += htemp->GetBinContent( ii, jj );
                            }
                            else
                            {
                                imean += htemp->GetBinContent( ii, jj ) * htemp->GetBinContent( ii, jj );
                            }
                            ierror += htemp->GetBinError( ii, jj ) * htemp->GetBinError( ii, jj );
                            z++;
                        }
                    }
                }
                if( z > 0. )
                {
                    if( !bError )
                    {
                        hs->SetBinContent( i, j, imean / z );
                    }
                    else
                    {
                        hs->SetBinContent( i, j, sqrt( imean ) / z );
                    }
                    hs->SetBinError( i, j, sqrt( ierror ) / z );
                }
            }
        }
        delete htemp;
    }
    return hs;
}

/*
 * inter/extrapolation function made for lookup tables
 *
 * use iname=="fitexpo" for energy lookup tables
 *
 * use iname=="fitpol3" for mscw/mscl lookup tables
 *
 * use iname=="fitpol2" for mscw/mscl lookup tables
 *
 * do fits in slices in X
 */
TH2F* VInterpolate2DHistos::doLogLinearExtrapolation( TH2F* h, string iname, TH2F* hNevents, int iMinEvents )
{
    if( !h )
    {
        return 0;
    }
    
    char hname[600];
    sprintf( hname, "%s_%s", h->GetName(), iname.c_str() );
    TH2F* hs = ( TH2F* )h->Clone( hname );
    if( !hs )
    {
        return 0;
    }
    
    // loop over all bin in distance R
    for( int i = 1; i <= hs->GetNbinsY(); i++ )
    {
        TH1D* h = hs->ProjectionX( "d", i, i );
        if( h )
        {
            // count number of valid bins
            int z = 0;
            for( int b = 1; b <= h->GetNbinsX(); b++ )
            {
                if( h->GetBinContent( b ) > 0. )
                {
                    z++;
                }
            }
            // require at least four points for fit
            if( z > 3 )
            {
                /////////////////////////////////////
                // error calculation for fitted bins
                double i_Eav = 0.;
                double i_n = 0.;
                // calculate relative average error
                for( int b = 1; b <= hs->GetNbinsX(); b++ )
                {
                    if( hs->GetBinError( b, i ) > 0. )
                    {
                        i_Eav += hs->GetBinError( b, i ) / hs->GetBinContent( b, i );
                        i_n++;
                    }
                }
                if( i_n > 0. )
                {
                    i_Eav /= i_n;
                }
                /////////////////////////////////////
                // Fit and fill empty bins
                string iFitFunction = "";
                if( iname.find( "expo" ) != string::npos )
                {
                    iFitFunction = "expo";
                }
                else if( iname.find( "pol3" ) != string::npos )
                {
                    iFitFunction = "pol3";
                }
                else if( iname.find( "pol2" ) != string::npos )
                {
                    iFitFunction = "pol2";
                }
                h->Fit( iFitFunction.c_str(), "Q0" );
                TF1* f = h->GetFunction( iFitFunction.c_str() );
                
                if( f )
                {
                    // fill histogram with fitted values
                    // this is done only for those bins with less
                    // events than the given number
                    for( int b = 1; b <= hs->GetNbinsX(); b++ )
                    {
                        if( hNevents && hNevents->GetBinContent( b, i ) <= iMinEvents )
                        {
                            hs->SetBinContent( b, i, f->Eval( hs->GetXaxis()->GetBinCenter( b ) ) );
                            // error is set using the average relative error
                            if( i_Eav > 0. )
                            {
                                hs->SetBinError( b, i, i_Eav * f->Eval( hs->GetXaxis()->GetBinCenter( b ) ) );
                            }
                        }
                    }
                }
            }
            delete h;
        }
    }
    return hs;
}

/*
 *  Gaussian smoothing / interpolation
 *
 *
 *   **THIS IS UNTESTED AND SHOULD NOT BE USED**
*/
TH2F* VInterpolate2DHistos::doGaussianInterpolation( TH2F* h, string iname, TH2F* hNevents, int nGausN, double nWidth )
{
    if( !h || !hNevents )
    {
        return 0;
    }
    
    char hname[600];
    sprintf( hname, "%s_%s", h->GetName(), iname.c_str() );
    TH2F* hs = ( TH2F* )h->Clone( hname );
    
    sprintf( hname, "%s_%s_%s", h->GetName(), iname.c_str(), "temp" );
    TProfile2D h2D( hname, "", hs->GetNbinsX(), hs->GetXaxis()->GetXmin(), hs->GetXaxis()->GetXmax(),
                    hs->GetNbinsY(), hs->GetYaxis()->GetXmin(), hs->GetYaxis()->GetXmax(),
                    -100., 100. );
                    
    double x = 0.;
    double y = 0.;
    double xc = 0.;
    double yc = 0.;
    double z = 0.;
    
    double xWidth = hs->GetXaxis()->GetBinWidth( 1 );
    double yWidth = hs->GetYaxis()->GetBinWidth( 1 );
    
    for( int i = 1; i <= hs->GetNbinsX(); i++ )
    {
        xc = hs->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= hs->GetNbinsY(); j++ )
        {
            yc = hs->GetYaxis()->GetBinCenter( j );
            z  = hs->GetBinContent( i, j );
            int n  = ( int )( hNevents->GetBinContent( i, j ) * nGausN );
            
            // loop over all bin entries
            for( int k = 0; k < n; k++ )
            {
                x = fRandom->Gaus( xc, xWidth * nWidth / 2. );
                y = fRandom->Gaus( yc, yWidth * nWidth / 2. );
                h2D.Fill( x, y, z );
            }
        }
    }
    
    // copy histograms
    for( int i = 1; i <= hs->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= hs->GetNbinsY(); j++ )
        {
            hs->SetBinContent( i, j, h2D.GetBinContent( i, j ) );
        }
    }
    return hs;
}
