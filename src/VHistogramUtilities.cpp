/*! \class VHistogramUtilities
    \brief utility class to manipulate histograms

*/

#include "VHistogramUtilities.h"

VHistogramUtilities::VHistogramUtilities()
{
    setDebug( false );
}



/*
    calculate residuals between a 1D histogram and a function

    input:

    string iname :  histogram name for return histogram
    TH1 *h :        input 1D histogram
    TF1 *f :        input 1D function

    return:

    TH1D* :         residual
*/
TH1D* VHistogramUtilities::get_ResidualHistogram_from_TF1( string iname, TH1* h, TF1* f )
{
    if( iname.size() == 0 )
    {
        return 0;
    }
    if( !h || !f )
    {
        return 0;
    }
    
    TH1D* hDiff = new TH1D( iname.c_str(), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
    
    for( int i = 1; i <= hDiff->GetNbinsX(); i++ )
    {
        if( h->GetBinContent( i ) )
        {
            hDiff->SetBinContent( i, ( h->GetBinContent( i ) - f->Eval( h->GetBinCenter( i ) ) ) / h->GetBinContent( i ) );
            hDiff->SetBinError( i, h->GetBinError( i ) / h->GetBinContent( i ) );
        }
    }
    
    return hDiff;
}

/*

    convert a graph to a histogram

    in cases that there are more points per bin, average over these values

*/
void VHistogramUtilities::fill_Graph_in_Histogram( TGraphAsymmErrors* g, TH1* h, bool bRequireYGtZero )
{
    if( !g || !h )
    {
        return;
    }
    
    double x = 0.;
    double y = 0.;
    
    string hname = h->GetName();
    hname += "_P";
    TProfile* iP = new TProfile( hname.c_str(), "", h->GetNbinsX(),
                                 h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
    hname += "_E";
    TProfile* iE = new TProfile( hname.c_str(), "", h->GetNbinsX(),
                                 h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
                                 
    for( int j = 0; j < g->GetN(); j++ )
    {
        g->GetPoint( j, x, y );
        
        if( bRequireYGtZero && y <= 0. )
        {
            continue;
        }
        
        iP->Fill( x, y );
        iE->Fill( x, 0.5 * ( g->GetErrorYlow( j ) + g->GetErrorYhigh( j ) ) );
    }
    
    for( int i = 0; i <= h->GetNbinsX(); i++ )
    {
        h->SetBinContent( i, iP->GetBinContent( i ) );
        h->SetBinError( i, iE->GetBinContent( i ) );
    }
    
    delete iP;
    delete iE;
}



/*
    get systematic error in reconstruction from a 2D histogram
    (essentially a profile with mean or median)

    options for iMeanType:

       mean   :  error bars are the error of the mean
       meanS  :  error bars are the width of the distribution
       median :  median value

*/
TGraphErrors* VHistogramUtilities::get_Profile_from_TH2D( TH2D* iP, TGraphErrors* g, string iMeanType, int rbin, double iXaxisValue, double iMinusValue )
{
    if( !iP )
    {
        return 0;
    }
    
    if( g == 0 )
    {
        g = new TGraphErrors( 1 );
    }
    
    int zz = 0;
    //////////////////////////////////////////////////////////////////////
    // median
    if( iMeanType == "median" )
    {
        double i_a[] = { 0.32, 0.5, 0.68 };
        double i_b[] = { 0.0,  0.0, 0.0  };
        
        for( int b = 1; b <= iP->GetNbinsX(); b++ )
        {
            if( iP->GetXaxis()->GetBinCenter( b ) < iXaxisValue )
            {
                continue;
            }
            
            string hname = iP->GetName();
            hname += "_AA";
            TH1D* h = ( TH1D* )iP->ProjectionY( hname.c_str(), b, b );
            if( h && h->GetEntries() > 3. )
            {
                h->GetQuantiles( 3, i_b, i_a );
                
                g->SetPoint( zz, iP->GetXaxis()->GetBinCenter( b ), i_b[1] - iMinusValue );
                g->SetPointError( zz, 0., ( i_b[0] + i_b[2] ) / 2. / TMath::Sqrt( h->GetEntries() ) );
                zz++;
                delete h;
            }
        }
    }
    //////////////////////////////////////////////////////////////////////
    // mean
    else
    {
        for( int b = 1; b <= iP->GetNbinsX(); b++ )
        {
            if( iP->GetXaxis()->GetBinCenter( b ) < iXaxisValue )
            {
                continue;
            }
            
            TH1D* h = ( TH1D* )iP->ProjectionY( "a", b, b );
            if( h && h->GetEntries() > 3. )
            {
                g->SetPoint( zz, iP->GetXaxis()->GetBinCenter( b ), h->GetMean() - iMinusValue );
                if( iMeanType == "meanS" )
                {
                    g->SetPointError( zz, 0., h->GetRMS() );
                }
                else
                {
                    g->SetPointError( zz, 0., h->GetRMS() / TMath::Sqrt( h->GetEntries() ) );
                }
                zz++;
            }
            if( h )
            {
                delete h;
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////
    // rebinning
    if( rbin == 2 )
    {
        int z = 0;
        for( int b = 0; b < g->GetN() - 1; b = b + 2 )
        {
            double x1, y1, x2, y2 ;
            
            g->GetPoint( b, x1, y1 );
            g->GetPoint( b + 1, x2, y2 );
            g->SetPoint( z, ( x1 + x2 ) / 2., ( y1 + y2 ) / 2. );
            z++;
        }
        g->Set( ( int )( g->GetN() / rbin ) );
    }
    return g;
}

/*

    interpolate a response matrix

    This is used in the CTA sensitivity estimations, where the response matrixes calculated from
    proton simulation are often not well defined

*/

TH2* VHistogramUtilities::interpolateResponseMatrix( TH2* hResponseMatrix, string iNewHistoName )
{
    if( !hResponseMatrix )
    {
        return 0;
    }
    
    // clone input histogram
    char hname[600];
    if( iNewHistoName.size() == 0 )
    {
        sprintf( hname, "%s_IRM", hResponseMatrix->GetName() );
    }
    else
    {
        sprintf( hname, "%s", iNewHistoName.c_str() );
    }
    TH2* iHInter = ( TH2* )hResponseMatrix->Clone( hname );
    
    // calculate median energy for each y-bin (true energy axis)
    double xq[] = { 0.16, 0.50, 0.84 };
    double yq[] = { 0.0,  0.0, 0.0  };
    TH1D hOff( "hMeanOffset", "", 100, -2., 2. );
    vector< int > iYBinToInterpolate;
    int iyBinFirstToInterpolate = -1;
    for( int i = 1; i <= hResponseMatrix->GetNbinsY(); i++ )
    {
        TH1F* h = ( TH1F* )hResponseMatrix->ProjectionX( "p_y", i, i );
        if( h )
        {
            if( h->GetEntries() > 0 )
            {
                h->GetQuantiles( 3, yq, xq );
            }
            // count number of filled bins
            unsigned int i_filled = 0;
            for( int j = 1; j <= h->GetNbinsX(); j++ )
            {
                if( h->GetBinContent( j ) > 0 )
                {
                    i_filled++;
                }
            }
            if( i_filled <= 3 )
            {
                if( h->GetEntries() > 0 )
                {
                    hOff.Fill( yq[1] - hResponseMatrix->GetYaxis()->GetBinCenter( i ) );
                }
                iYBinToInterpolate.push_back( i );
            }
            else
            {
                if( iyBinFirstToInterpolate < 0 && i > 1 )
                {
                    iyBinFirstToInterpolate = i;
                }
            }
        }
    }
    hOff.GetQuantiles( 3, yq, xq );
    //////////////////////////////////
    // interpolate the missing bins
    double y = 0;
    double x = 0.;
    // rms shouldn't be larger than bin width
    double rms = yq[2] - yq[1];
    if( rms > iHInter->GetXaxis()->GetBinWidth( 1 ) )
    {
        rms = iHInter->GetXaxis()->GetBinWidth( 1 );
    }
    for( unsigned int i = 0; i < iYBinToInterpolate.size(); i++ )
    {
        // don't interpolate the smallest energies
        if( iYBinToInterpolate[i] <= iyBinFirstToInterpolate )
        {
            continue;
        }
        // fill a normal distribution at the median position
        y = iHInter->GetYaxis()->GetBinCenter( iYBinToInterpolate[i] );
        for( unsigned int j = 0; j < 100; j++ )
        {
            x = gRandom->Gaus( y + yq[1], rms );
            iHInter->Fill( x, y );
        }
    }
    normalizeTH2D_y( iHInter );
    
    return iHInter;
}

/*

   normalize a matrix in each y-row

*/
bool VHistogramUtilities::normalizeTH2D_y( TH2* h )
{
    if( !h )
    {
        return false;
    }
    
    double i_sum = 0.;
    for( int i = 1; i <= h->GetNbinsY(); i++ )
    {
        i_sum = 0.;
        for( int j = 1; j <= h->GetNbinsX(); j++ )
        {
            i_sum += h->GetBinContent( j, i );
        }
        if( i_sum > 0. )
        {
            for( int j = 1; j <= h->GetNbinsX(); j++ )
            {
                h->SetBinContent( j, i, h->GetBinContent( j, i ) / i_sum );
            }
        }
    }
    return true;
}

TH1D* VHistogramUtilities::get_Cumulative_Histogram( TH1D* iH_in, bool iNormalize, bool iLeft_to_right, double i_bin_value, double i_min_value )
{
    if( !iH_in )
    {
        return 0;
    }
    
    // clone input histogram
    char hname[600];
    sprintf( hname, "%s_CUMU", iH_in->GetName() );
    TH1D* iH_out = ( TH1D* )iH_in->Clone( hname );
    iH_out->Reset();
    
    float z = iH_out->GetNbinsX();
    
    if( iLeft_to_right )
    {
        iH_out->SetBinContent( 1, iH_in->GetBinContent( 1 ) );
        // loop over all bins
        for( int i = 2; i <= iH_in->GetNbinsX(); i++ )
        {
            if( iH_in->GetBinCenter( i ) > i_bin_value )
            {
                z = i - 1;
                break;
            }
            // avoid overshooting of from e.g. negative excesses
            if( iH_in->GetBinContent( i ) > i_min_value )
            {
                iH_out->SetBinContent( i, iH_in->GetBinContent( i ) + iH_out->GetBinContent( i - 1 ) );
            }
            else
            {
                iH_out->SetBinContent( i, iH_out->GetBinContent( i - 1 ) );
            }
        }
    }
    else
    {
        iH_out->SetBinContent( iH_in->GetNbinsX(), iH_in->GetBinContent( iH_in->GetNbinsX() ) );
        // loop over all bins
        for( int i = iH_in->GetNbinsX() - 1; i >= 1; i-- )
        {
            if( iH_in->GetBinCenter( i ) < i_bin_value )
            {
                z = i + 1;
                break;
            }
            // avoid overshooting of from e.g. negative excesses
            if( iH_in->GetBinContent( i ) > i_min_value )
            {
                iH_out->SetBinContent( i, iH_in->GetBinContent( i ) + iH_out->GetBinContent( i + 1 ) );
            }
            else
            {
                iH_out->SetBinContent( i, iH_out->GetBinContent( i + 1 ) );
            }
        }
    }
    if( iNormalize && iH_out->GetBinContent( z ) > 0. )
    {
        iH_out->Scale( 1. / iH_out->GetBinContent( z ) );
    }
    
    iH_out->SetMaximum( iH_out->GetMaximum() * 1.2 );
    
    return iH_out;
}

/*

regioncode:  extra specifier for additional region cuts
             'a' - exclude x<0 (only use top    half of skymap)
             'b' - exclude x>0 (only use bottom half of skymap)

*/
TH1D* VHistogramUtilities::get_Bin_Distribution( TH2D* h, int ion, double rmax, double rSource, bool iDiff, TH2D* hTest,
        int iExcN, float* iExcX, float* iExcY, float* iExcR1, float* iExcR2, float* iExcTheta, string regioncode )
{
    if( !h )
    {
        return 0;
    }
    
    // no distance cut for rmax < 0
    if( rmax < 0. )
    {
        rmax = 1.e6;
    }
    
    /////////////////////////////////
    // definition of 1D histogram
    int nbin = 100;
    double xmin = -7.1;
    double xmax = 10.9;
    
    char regioncode_histname[100] = "" ;
    if( regioncode.length() > 0 )
    {
        sprintf( regioncode_histname, "_%s", regioncode.c_str() ) ;
    }
    
    
    char hname[200];
    if( iDiff )
    {
        sprintf( hname, "hdiff1D_%d_%f%s", ion, rSource, regioncode_histname );
    }
    else
    {
        sprintf( hname, "hsig1D_%d_%f_%d%s", ion, rSource, iExcN, regioncode_histname );
    }
    
    TH1D* h1D = 0;
    if( gDirectory->Get( hname ) )
    {
        h1D = ( TH1D* )gDirectory->Get( hname );
        h1D->Reset();
    }
    else
    {
        h1D = new TH1D( hname, "", nbin, xmin, xmax );
    }
    
    // setup regioncode flags
    string regioncode_a = "a" ;
    bool   regioncodeflag_a = false ;
    if( regioncode_a.compare( regioncode ) == 0 )
    {
        regioncodeflag_a = true ;
    }
    string regioncode_b = "b" ;
    bool   regioncodeflag_b = false ;
    if( regioncode_b.compare( regioncode ) == 0 )
    {
        regioncodeflag_b = true ;
    }
    
    ////////////////////////////////////////////
    // get the distribution and fill histogram
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        double x_r = h->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            // exclude bins with no on counts
            if( hTest && hTest->GetBinContent( i, j ) == 0. )
            {
                continue;
            }
            
            double y_r = h->GetYaxis()->GetBinCenter( j );
            // radius cut
            if( TMath::Sqrt( x_r * x_r + y_r * y_r ) > rmax )
            {
                continue;
            }
            
            // regioncode cuts
            // regioncode = "a" , exclude x<0
            if( regioncodeflag_a )
            {
                if( x_r < 0 )
                {
                    continue ;
                }
            }
            // regioncode = "b" , exclude x>0
            else if( regioncodeflag_b )
            {
                if( x_r > 0 )
                {
                    continue ;
                }
            }
            
            // exclusion regions
            bool iBBreak = false;
            for( int e = 0; e < iExcN; e++ )
            {
                //if( ( x_r - iExcX[e] ) * ( x_r - iExcX[e] ) / ( iExcR1[e]*iExcR1[e] ) + ( y_r - iExcY[e] ) * ( y_r - iExcY[e] ) / ( iExcR2[e]*iExcR2[e] ) < 1. )
                if( TMath::Power( ( ( x_r - iExcX[e] ) * TMath::Cos( iExcTheta[e] * TMath::DegToRad() ) + ( y_r - iExcY[e] ) * TMath::Sin( iExcTheta[e] * TMath::DegToRad() ) ) / iExcR1[e], 2 ) + TMath::Power( ( ( x_r - iExcX[e] ) * TMath::Sin( iExcTheta[e] * TMath::DegToRad() ) - ( y_r - iExcY[e] ) * TMath::Cos( iExcTheta[e] * TMath::DegToRad() ) ) / iExcR2[e] , 2 ) < 1. )
                {
                    iBBreak = true;
                    break;
                }
            }
            if( iBBreak )
            {
                continue;
            }
            // exclude bins around center of sky map
            if( TMath::Sqrt( x_r * x_r + y_r * y_r ) > rSource )
            {
                h1D->Fill( h->GetBinContent( i, j ) );
            }
        }
    }
    return h1D;
}

bool VHistogramUtilities::get_Graph_from_Histogram( TH1* h, TGraphErrors* g, bool bIgnoreErrors, double iMinBinContent,
        double iXmin, double iXmax )
{
    if( !h || !g )
    {
        return false;
    }
    
    unsigned int z = 0;
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        if( h->GetXaxis()->GetBinCenter( i ) < iXmin )
        {
            continue;
        }
        if( h->GetXaxis()->GetBinCenter( i ) > iXmax )
        {
            continue;
        }
        if( h->GetBinContent( i ) > iMinBinContent )
        {
            g->SetPoint( z, h->GetXaxis()->GetBinCenter( i ), h->GetBinContent( i ) );
            if( bIgnoreErrors )
            {
                g->SetPointError( z, 0., 0. );
            }
            else
            {
                g->SetPointError( z, 0., h->GetBinError( i ) );
            }
            z++;
        }
    }
    return true;
}

bool VHistogramUtilities::get_Graph_from_Histogram( TH1* h, TGraphAsymmErrors* g, bool bIgnoreErrors, bool bLinXaxis,
        double iCutUnrealisticErrors, double iXmin, double iXmax )
{
    if( !h || !g )
    {
        return false;
    }
    
    unsigned int z = 0;
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        if( h->GetXaxis()->GetBinCenter( i ) < iXmin )
        {
            continue;
        }
        if( h->GetXaxis()->GetBinCenter( i ) > iXmax )
        {
            continue;
        }
        if( h->GetBinContent( i ) > 0. )
        {
            if( bIgnoreErrors )
            {
                g->SetPointEYlow( z, 0. );
                g->SetPointEYhigh( z, 0 );
            }
            else
            {
                // remove unrealistic errors (e.g. error is iCutUnrealisticErrorsx the bin content)
                if( iCutUnrealisticErrors > 0. && h->GetBinContent( i ) && h->GetBinError( i ) / h->GetBinContent( i ) > iCutUnrealisticErrors )
                {
                    g->SetPointEYlow( z, 0. );
                    g->SetPointEYhigh( z, 0. );
                }
                else
                {
                    g->SetPointEYlow( z, h->GetBinError( i ) );
                    g->SetPointEYhigh( z, h->GetBinError( i ) );
                }
            }
            if( !bLinXaxis )
            {
                g->SetPoint( z, h->GetXaxis()->GetBinCenter( i ), h->GetBinContent( i ) );
                g->SetPointEXlow( z, h->GetXaxis()->GetBinCenter( i ) - h->GetXaxis()->GetBinLowEdge( i ) );
                g->SetPointEXhigh( z, h->GetXaxis()->GetBinLowEdge( i ) + h->GetXaxis()->GetBinWidth( i ) - h->GetXaxis()->GetBinCenter( i ) );
            }
            else
            {
                if( h->GetXaxis()->GetBinCenter( i ) > 0. )
                {
                    g->SetPoint( z, TMath::Log10( h->GetXaxis()->GetBinCenter( i ) ), h->GetBinContent( i ) );
                    g->SetPointEXlow( z, TMath::Log10( h->GetXaxis()->GetBinCenter( i ) - h->GetXaxis()->GetBinLowEdge( i ) ) );
                    g->SetPointEXhigh( z, TMath::Log10( h->GetXaxis()->GetBinLowEdge( i ) + h->GetXaxis()->GetBinWidth( i ) - h->GetXaxis()->GetBinCenter( i ) ) );
                }
                else
                {
                    continue;
                }
            }
            z++;
        }
    }
    return true;
}

/*

    get histograms from CTA WP Phys sensitivity files

    not clear if this is the right place, but used by VInstrumentResponseFunctionReader and VSensitivityCalculator

*/
TH1F* VHistogramUtilities::get_CTA_IRF_Histograms( string iHistogramName, double iCameraOffset )
{
    TH1F* h = 0;
    TH2F* h2D = 0;
    char hname[200];
    // gamma-ray effective area
    // get on-axis results
    if( iCameraOffset <= 1.e-2 )
    {
        h = ( TH1F* )gDirectory->Get( iHistogramName.c_str() );
    }
    // get off-axis results
    else
    {
        sprintf( hname, "%s_offaxis", iHistogramName.c_str() );
        h2D = ( TH2F* )gDirectory->Get( hname );
        if( !h2D )
        {
            return 0;
        }
        sprintf( hname, "%s_px", h2D->GetName() );
        h = ( TH1F* )h2D->ProjectionX( hname, h2D->GetYaxis()->FindBin( iCameraOffset ), h2D->GetYaxis()->FindBin( iCameraOffset ) );
        if( !h )
        {
            return 0;
        }
    }
    
    return h;
}

/*

   get histograms from CTA WP Phys sensitivity files

   2D version: return profile

*/

TH1F* VHistogramUtilities::get_CTA_IRF_Histograms_from2D( string iHistogramName, double iSummand )
{

    TH2F* h2D = ( TH2F* )gDirectory->Get( iHistogramName.c_str() );
    if( h2D )
    {
        TH1F* h = ( TH1F* )h2D->ProfileX();
        if( h && iSummand != 0. )
        {
            for( int i = 1; i <= h->GetNbinsX(); i++ )
            {
                h->SetBinContent( i, h->GetBinContent( i ) + iSummand );
            }
            return h;
        }
    }
    
    return 0;
}

int VHistogramUtilities::findBinInGraph( TGraph* g, double x )
{
    if( !g )
    {
        return -1;
    }
    
    double i_x = 0.;
    double i_x_low = 0.;
    double i_x_high = 0.;
    double i_y = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, i_x, i_y );
        i_x_low  = g->GetErrorXlow( i );
        i_x_high = g->GetErrorXhigh( i );
        
        if( x > i_x - i_x_low && x <= i_x + i_x_high )
        {
            return i;
        }
    }
    
    return -1;
}

TH1* VHistogramUtilities::normalizeTH1( TH1* h, bool iIntegral )
{
    if( !h )
    {
        return 0;
    }
    
    double iSum = 0.;
    double iN = 0.;
    
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        iSum += h->GetBinContent( i );
        iN++;
    }
    
    double f = 1.;
    if( iSum > 0. )
    {
        if( iIntegral )
        {
            f =  1. / iSum;
        }
        else
        {
            f = iN / iSum;
        }
    }
    
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        h->SetBinError( i, h->GetBinError( i ) * f );
        h->SetBinContent( i, h->GetBinContent( i ) * f );
    }
    
    return h;
}

bool VHistogramUtilities::divide( TGraphAsymmErrors* g, TGraphAsymmErrors* g1, TGraph* g2 )
{
    if( !g || !g1 || !g2 )
    {
        return false;
    }
    
    double x1 = 0.;
    double y1 = 0.;
    
    int z = 0;
    for( int i = 0; i < g1->GetN(); i++ )
    {
        g1->GetPoint( i, x1, y1 );
        
        double y2 = g2->Eval( x1 );
        
        if( y1 != 0. )
        {
            g->SetPoint( z, x1, y2 / y1 );
            
            g->SetPointError( z, g1->GetErrorXlow( i ), g1->GetErrorXhigh( i ),
                              y2 * 0.5 * ( g1->GetErrorYhigh( i ) + g1->GetErrorYlow( i ) ) / y1 / y1,
                              y2 * 0.5 * ( g1->GetErrorYhigh( i ) + g1->GetErrorYlow( i ) ) / y1 / y1 );
            z++;
        }
    }
    
    return true;
}

/*

     divide two graphs by each other (take also care of errors)

*/
bool VHistogramUtilities::divide( TGraphAsymmErrors* g, TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, double epsilon )
{
    if( !g || !g1 || !g2 )
    {
        return false;
    }
    
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    
    int z = 0;
    for( int i = 0; i < g1->GetN(); i++ )
    {
        g1->GetPoint( i, x1, y1 );
        
        for( int j = 0; j < g2->GetN(); j++ )
        {
            g2->GetPoint( j, x2, y2 );
            
            if( TMath::Abs( x1 - x2 ) < epsilon && y1 != 0. && y2 != 0. )
            {
                double y = y2 / y1;
                g->SetPoint( z, x1, y );
                
                
                // double iErr1 = 1. / y1 * 0.5 * ( g1->GetErrorYhigh( i ) + g1->GetErrorYlow( i ) );
                double iErr1 = y2 / y1 / y1 * 0.5 * ( g1->GetErrorYhigh( i ) + g1->GetErrorYlow( i ) );
                // double iErr2 = 1. / y2 * 0.5 * ( g2->GetErrorYhigh( j ) + g2->GetErrorYlow( j ) );
                double iErr2 = 0.;
                if( g2->GetErrorYhigh( j ) > 0. && g2->GetErrorYlow( j ) > 0. )
                {
                    iErr2 = 1. / y1 * 0.5 * ( g2->GetErrorYhigh( j ) + g2->GetErrorYlow( j ) );
                }
                
                double yErr = y * sqrt( iErr1 * iErr1 + iErr2 * iErr2 );
                
                if( yErr > y )
                {
                    continue;
                }
                
                g->SetPointError( z, g1->GetErrorXlow( i ), g1->GetErrorXhigh( i ), yErr, yErr );
                z++;
                break;
            }
        }
    }
    
    return true;
}

/*!

    warning: not clear if '0' is a good return value for a bad value

*/
double VHistogramUtilities::interpolateTH2D( TH2* h, double x, double y )
{
    if( !h )
    {
        return 0.;
    }
    
    if( x >= h->GetXaxis()->GetXmin() && x <= h->GetXaxis()->GetXmax()
            && y >= h->GetYaxis()->GetXmin() && y <= h->GetYaxis()->GetXmax() )
    {
        return h->Interpolate( x, y );
    }
    
    return 0.;
}

TH2D* VHistogramUtilities::calculateContainmentDistance( TH2D* h, string inewHistogramName )
{
    if( !h )
    {
        return 0;
    }
    
    TH2D* hNew = new TH2D( inewHistogramName.c_str(), h->GetTitle(), h->GetNbinsX(),
                           h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                           100, 0., 1. );
    hNew->SetXTitle( h->GetXaxis()->GetTitle() );
    hNew->SetYTitle( "containment" );
    hNew->SetZTitle( "angular resolution [deg]" );
    
    // loop over all energy bins
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        double iTot = 0.;
        for( int j = 1; j < h->GetNbinsY(); j++ )
        {
            iTot += h->GetBinContent( i, j );
        }
        if( iTot < 1.e-12 )
        {
            continue;
        }
        // integrate from high to low y values
        double iSum = 0.;
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            iSum += h->GetBinContent( i, j );
            hNew->SetBinContent( i, hNew->GetYaxis()->FindBin( iSum / iTot ), h->GetYaxis()->GetBinCenter( j ) );
        }
    }
    return hNew;
}


/*

   reduce size of histogram by removing all irrelevant bins (empty bins)

*/
TH2F* VHistogramUtilities::reduce2DHistogramSize( TH2* h, string inewHistogramName )
{
    if( !h )
    {
        return 0;
    }
    
    int nBinX_min = 99999;
    int nBinX_max = -1;
    int nBinY_min = 99999;
    int nBinY_max = -1;
    
    // get minima and maxima
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            if( h->GetBinContent( i, j ) > 0. )
            {
                if( i < nBinX_min )
                {
                    nBinX_min = i;
                }
                if( j < nBinY_min )
                {
                    nBinY_min = j;
                }
                if( i > nBinX_max )
                {
                    nBinX_max = i;
                }
                if( j > nBinY_max )
                {
                    nBinY_max = j;
                }
            }
        }
    }
    if( nBinX_max < 0 )
    {
        nBinX_min = 1;
        nBinX_max = h->GetNbinsX();
        nBinY_min = 1;
        nBinY_max = h->GetNbinsY();
    }
    nBinY_min = 1;
    // create new histogram with reduced binning
    float xmin = h->GetXaxis()->GetBinLowEdge( 1 );
    if( nBinX_min > 0 )
    {
        xmin = h->GetXaxis()->GetBinLowEdge( nBinX_min );
    }
    float xmax = h->GetXaxis()->GetBinLowEdge( 1 ) +  h->GetXaxis()->GetBinWidth( 1 );
    if( nBinX_max > 0 )
    {
        xmax = h->GetXaxis()->GetBinLowEdge( nBinX_max ) + h->GetXaxis()->GetBinWidth( nBinX_max );
    }
    float ymin = h->GetYaxis()->GetBinLowEdge( 1 );
    float ymax = h->GetYaxis()->GetBinLowEdge( 1 ) +  h->GetYaxis()->GetBinWidth( 1 );
    if( ymax > 0 )
    {
        ymax = h->GetYaxis()->GetBinLowEdge( nBinY_max ) + h->GetYaxis()->GetBinWidth( nBinY_max );
    }
    
    TH2F* hNew = new TH2F( inewHistogramName.c_str(), "", nBinX_max - nBinX_min + 1, xmin, xmax, nBinY_max, ymin, ymax );
    hNew->SetXTitle( h->GetXaxis()->GetTitle() );
    hNew->SetYTitle( h->GetYaxis()->GetTitle() );
    hNew->SetZTitle( h->GetZaxis()->GetTitle() );
    hNew->Sumw2();
    
    for( int i = nBinX_min; i <= nBinX_max; i++ )
    {
        for( int j = nBinY_min; j <= nBinY_max; j++ )
        {
            if( h->GetBinContent( i, j ) > 0. )
            {
                hNew->SetBinContent( i - nBinX_min + 1, j - nBinY_min + 1, h->GetBinContent( i, j ) );
                hNew->SetBinError( i - nBinX_min + 1, j - nBinY_min + 1, h->GetBinError( i, j ) );
            }
        }
    }
    
    return hNew;
}
