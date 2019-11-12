/*! \class VPlotUtilities


*/

#include "VPlotUtilities.h"

VPlotUtilities::VPlotUtilities()
{
    setCanvasSize();
    setPlottingStyle();
    setColorAxisPalette();
    setColorAxisDataVector_minmax();
    fColorAxis_axis = 0;
}


VPlotUtilities::~VPlotUtilities()
{

}


/*
    xmin/xmax on linear axis scale
    xmin == xmax: get values from histogram
*/
void VPlotUtilities::plot_nullHistogram( TPad* c, TH1* h, bool bLogX, bool bLogY, double yTitleOffset, double xmin, double xmax )
{
    if( !c || !h )
    {
        return;
    }
    TPad* iPad = ( TPad* )c;
    
    if( TMath::Abs( xmax - xmin ) < 1.e-4 )
    {
        xmin = TMath::Power( 10., h->GetXaxis()->GetXmin() );
        xmax = TMath::Power( 10., h->GetXaxis()->GetXmax() );
    }
    
    if( !bLogX )
    {
        h->GetYaxis()->SetTitleOffset( yTitleOffset );
        h->Draw();
    }
    else
    {
        h->Draw( "AH" );
        iPad->Update();
        
        // X-axis
        TGaxis* x1 = new TGaxis( iPad->GetUxmin(), iPad->GetUymin(), iPad->GetUxmax(), iPad->GetUymin(), xmin, xmax, 510, "G" );
        string xtitle = h->GetXaxis()->GetTitle();
        if( xtitle.size() > 0 && xtitle.find( "log_{10}" ) != string::npos && xtitle.find( "log_{10}" ) != xtitle.size() )
        {
            x1->SetTitle( xtitle.substr( xtitle.find( "log_{10}" ) + 9, xtitle.size() ).c_str() );
        }
        else
        {
            x1->SetTitle( h->GetXaxis()->GetTitle() );
        }
        x1->SetLabelOffset( 0.003 );
        x1->SetTitleOffset( 1.25 );
        x1->Draw();
        
        TGaxis* x2 = new TGaxis( iPad->GetUxmin(), iPad->GetUymax(), iPad->GetUxmax(), iPad->GetUymax(), xmin, xmax, 510, "-UG" );
        x2->Draw();
        
        // y-axis
        if( bLogY )
        {
            TGaxis* y1 = new TGaxis( iPad->GetUxmin(), iPad->GetUymin(), iPad->GetUxmin(), iPad->GetUymax(), h->GetMinimum(), h->GetMaximum(), 510, "G" );
            y1->SetTitle( h->GetYaxis()->GetTitle() );
            y1->SetTitleOffset( yTitleOffset );
            y1->Draw();
            return;
            
            TGaxis* y2 = new TGaxis( iPad->GetUxmax(), iPad->GetUymin(), iPad->GetUxmax(), iPad->GetUymax(), h->GetMinimum(), h->GetMaximum(), 510, "+GU" );
            y2->Draw();
        }
        else
        {
            TGaxis* y1 = new TGaxis( iPad->GetUxmin(), iPad->GetUymin(), iPad->GetUxmin(), iPad->GetUymax(), h->GetMinimum(), h->GetMaximum(), 510, "" );
            y1->SetTitle( h->GetYaxis()->GetTitle() );
            y1->SetTitleOffset( yTitleOffset );
            y1->Draw();
            
            TGaxis* y2 = new TGaxis( iPad->GetUxmax(), iPad->GetUymin(), iPad->GetUxmax(), iPad->GetUymax(), h->GetMinimum(), h->GetMaximum(), 510, "+U" );
            y2->SetTitleOffset( yTitleOffset );
            y2->Draw();
        }
    }
}


void VPlotUtilities::setHistogramPlottingStyle( TH1* h )
{
    setHistogramPlottingStyle( h, fPlottingColor, fPlottingLineWidth, fPlottingMarkerSize, fPlottingMarkerStyle, 1, fPlottingFillStyle );
}


void VPlotUtilities::setHistogramPlottingStyle( TH1* his, int icolor, double iwidth, double isize, int imarker, int irebin, int iFillStyle )
{
    if( !his )
    {
        return;
    }
    
    his->SetLineColor( icolor );
    his->SetMarkerColor( icolor );
    his->SetLineWidth( ( Width_t )iwidth );
    his->SetMarkerSize( isize );
    if( imarker != 0 )
    {
        his->SetMarkerStyle( imarker );
    }
    if( iFillStyle >= 0 )
    {
        his->SetFillStyle( iFillStyle );
        his->SetFillColor( icolor );
    }
    his->SetStats( 0 );
    his->GetYaxis()->SetTitleOffset( 1.3 );
    if( irebin != 1 )
    {
        his->Rebin( irebin );
    }
}

void VPlotUtilities::setFunctionPlottingStyle( TF1* fun )
{
    setFunctionPlottingStyle( fun, fPlottingColor, fPlottingLineWidth, fPlottingMarkerSize, fPlottingMarkerStyle, fPlottingFillStyle );
}

void VPlotUtilities::setFunctionPlottingStyle( TF1* fun, int icolor, double iwidth, double isize, int imarker, int iFillStyle )
{
    if( !fun )
    {
        return;
    }
    
    fun->SetLineColor( icolor );
    fun->SetMarkerColor( icolor );
    fun->SetLineWidth( ( Width_t )iwidth );
    fun->SetMarkerSize( isize );
    if( imarker != 0 )
    {
        fun->SetMarkerStyle( imarker );
    }
    if( iFillStyle >= 0 )
    {
        fun->SetFillStyle( iFillStyle );
        fun->SetFillColor( icolor );
    }
}

void VPlotUtilities::setHistogramPlottingStyle( TH2D* his, double iminF )
{

    if( !his )
    {
        return;
    }
    
    his->SetStats( 0 );
    
    if( iminF > -10. )
    {
        double idiffmin = 99999;
        for( int i = 1; i <= his->GetNbinsX(); i++ )
        {
            for( int j = 1; j <= his->GetNbinsY(); j++ )
            {
                if( his->GetBinContent( i, j ) > -9999. &&  his->GetBinContent( i, j ) < idiffmin )
                {
                    idiffmin = his->GetBinContent( i, j );
                }
            }
        }
        his->SetMinimum( idiffmin * iminF );
    }
    if( his->GetMinimum() == his->GetMaximum() )
    {
        his->SetMinimum( his->GetMinimum() - 0.1 );
        his->SetMaximum( his->GetMaximum() + 0.1 );
    }
    
}

void VPlotUtilities::setArrowPlottingStyle( TArrow* a, int icolor, double iwidth, int iLineStyle )
{
    if( !a )
    {
        return;
    }
    
    a->SetLineColor( icolor );
    a->SetLineWidth( ( Width_t )iwidth );
    a->SetLineStyle( ( Style_t )iLineStyle );
    a->SetFillColor( icolor );
}


void VPlotUtilities::setArrowPlottingStyle( TArrow* a )
{
    if( !a )
    {
        return;
    }
    
    a->SetLineColor( fPlottingColor );
    a->SetLineWidth( ( Width_t )fPlottingLineWidth );
    a->SetLineStyle( ( Style_t )fPlottingLineStyle );
    a->SetFillColor( fPlottingColor );
    a->SetFillStyle( fPlottingFillStyle );
}


void VPlotUtilities::setGraphPlottingStyle( TGraph* g )
{
    setGraphPlottingStyle( g, fPlottingColor, fPlottingLineWidth, fPlottingMarkerStyle, fPlottingMarkerSize, fPlottingFillStyle, fPlottingLineStyle );
}


void VPlotUtilities::setGraphPlottingStyle( TGraph* g, int icolor, double iwidth, int imarker, double isize, int ifillstyle, int iLineStyle )
{
    if( !g )
    {
        return;
    }
    
    g->SetTitle( "" );
    g->SetLineColor( icolor );
    g->SetLineWidth( ( Width_t )iwidth );
    g->SetLineStyle( ( Style_t )iLineStyle );
    g->SetMarkerColor( icolor );
    g->SetMarkerSize( isize );
    g->SetMarkerStyle( imarker );
    g->SetFillColor( icolor );
    g->SetFillStyle( ifillstyle );
}


void VPlotUtilities::default_settings()
{
    gStyle->SetTitleOffset( 1.2, "Y" );
    gStyle->SetPadGridX( 0 );
    gStyle->SetPadGridY( 0 );
    gROOT->SetStyle( "Plain" );
    gStyle->SetPalette( 55 );
    setPlotHistogramTitle();
}


void VPlotUtilities::setColorAxisPalette( int palette, int ncolors )
{
    gStyle->SetPalette( palette );
    if( ncolors > 0 )
    {
        gStyle->SetNumberContours( ncolors );
    }
    fColorAxis_ncolor = gStyle->GetNumberOfColors();
    fColorAxis_ncont  = gStyle->GetNumberContours();
}


void VPlotUtilities::setColorAxisDataVector_minmax( double imin, double imax )
{
    fColorAxis_vmin = imin;
    fColorAxis_vmax = imax;
}


int VPlotUtilities::getColorAxisColor( double iV )
{
    double wls  = fColorAxis_vmax - fColorAxis_vmin;
    double scale = double( fColorAxis_ncont ) / wls;
    
    int color = int( 0.01 + ( iV - fColorAxis_vmin ) * scale );
    
    return gStyle->GetColorPalette( int( ( color + 0.99 ) * float( fColorAxis_ncolor - 1 ) / float( fColorAxis_ncont ) ) );
}


TGaxis* VPlotUtilities::getColorAxisAxis( double x1, double x2, double y1, double y2, string AxisTitle, Int_t ndiv, string iOption )
{
    double wls  = fColorAxis_vmax - fColorAxis_vmin;
    double scale = double( fColorAxis_ncont ) / wls;
    
    int color, theColor;
    double w1, w2;
    double ymin = 0.;
    double ymax = 0.;
    
    // plot axis on the right side
    TBox* iBox = new TBox();
    for( int i = 0; i < fColorAxis_ncont; i++ )
    {
        w1 = fColorAxis_vmin + wls / fColorAxis_ncont * i;
        if( w1 < fColorAxis_vmin )
        {
            w1 = fColorAxis_vmin;
        }
        w2 = fColorAxis_vmax;
        if( i < fColorAxis_ncont - 1 )
        {
            w2 = fColorAxis_vmin + wls / fColorAxis_ncont * ( i + 1 );
        }
        if( w2 <= fColorAxis_vmin )
        {
            continue;
        }
        ymin = y1 + ( w1 - fColorAxis_vmin ) * ( y2 - y1 ) / wls;
        ymax = y1 + ( w2 - fColorAxis_vmin ) * ( y2 - y1 ) / wls;
        color = int( 0.01 + ( w1 - fColorAxis_vmin ) * scale );
        theColor = int( ( color + 0.99 ) * float( fColorAxis_ncolor ) / float( fColorAxis_ncont ) );
        iBox->SetFillColor( gStyle->GetColorPalette( theColor ) );
        iBox->DrawBox( x1, ymin, x2, ymax );
    }
    fColorAxis_axis = new TGaxis( x2, y1, x2, y2, fColorAxis_vmin, fColorAxis_vmax, ndiv, iOption.c_str() );
    fColorAxis_axis->SetTitle( AxisTitle.c_str() );
    fColorAxis_axis->SetLabelSize( 0.02 );
    fColorAxis_axis->SetTitleSize( 0.02 );
    fColorAxis_axis->SetTitleOffset( 1.4 );
    fColorAxis_axis->Draw();
    
    return fColorAxis_axis;
}


void VPlotUtilities::setPadMargins( TCanvas* c, int nPads, double lM, double rM )
{
    if( !c )
    {
        return;
    }
    
    if( nPads > 1 )
    {
        for( int i = 0; i < nPads; i++ )
        {
            if( c->GetPad( i + 1 ) )
            {
                if( lM > -999. )
                {
                    c->GetPad( i + 1 )->SetLeftMargin( lM );
                }
                if( rM > -999. )
                {
                    c->GetPad( i + 1 )->SetRightMargin( rM );
                }
            }
            else
            {
                cout << "pad not found: " << c->GetName() << "\t" << i + 1 << endl;
            }
        }
    }
    else
    {
        if( lM > -999. )
        {
            c->SetLeftMargin( lM );
        }
        if( rM > -999. )
        {
            c->SetRightMargin( rM );
        }
    }
}


void VPlotUtilities::setTitles( TH1* his, string iname, string ititle, string ytitle )
{
    string itemp;
    // set names and titles
    itemp = his->GetName();
    itemp.replace( itemp.rfind( "on" ), 2, iname );
    his->SetName( itemp.c_str() );
    itemp = his->GetTitle();
    if( itemp.rfind( "(on)" ) < itemp.size() )
    {
        itemp.replace( itemp.rfind( "(on)" ), 4, ititle );
    }
    else
    {
        itemp += ititle;
    }
    his->SetTitle( itemp.c_str() );
    if( ytitle.size() > 0 )
    {
        his->GetYaxis()->SetTitle( ytitle.c_str() );
    }
}


void VPlotUtilities::setBlackAndWhitePalette()
{
    const int fBWNum = 30;
    int fBWPalette[fBWNum];
    int iOff = 0;
    int iCol = 9080;
    for( int i = 0; i < fBWNum; i++ )
    {
        TColor* color = new TColor( iCol + i, 1 - float( i + iOff ) / ( fBWNum + iOff ), 1 - float( i + iOff ) / ( fBWNum + iOff ), 1 - float( i + iOff ) / ( fBWNum + iOff ), "" );
        color->GetNumber();
        fBWPalette[i] = iCol + i;
    }
    gStyle->SetPalette( fBWNum, fBWPalette );
    gStyle->SetNumberContours( 20 );
}


TH2D* VPlotUtilities::removeOuterRing( TH2D* h, double r, double ivalue )
{
    if( !h )
    {
        return 0;
    }
    
    r *= r;
    
    double x = 0.;
    double y = 0.;
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        x = h->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            y = h->GetYaxis()->GetBinCenter( j );
            
            if( x * x + y * y > r )
            {
                h->SetBinContent( i, j, ivalue );
            }
            //       if( x*x + (y-0.5)*(y-0.5) > r ) h->SetBinContent( i, j, ivalue );
        }
    }
    
    return h;
}

unsigned int VPlotUtilities::setPlottingAxis( string iName, string iAxis, bool iLog, double xmin, double xmax, string iAxisTitle )
{
    if( fPlottingAxisData.find( iName ) == fPlottingAxisData.end() )
    {
        fPlottingAxisData[iName] = new VPlottingAxisData();
    }
    
    fPlottingAxisData[iName]->fName  = iName;
    fPlottingAxisData[iName]->fAxis  = iAxis;
    fPlottingAxisData[iName]->fAxisTitle = iAxisTitle;
    fPlottingAxisData[iName]->fLogAxis = iLog;
    fPlottingAxisData[iName]->fMinValue = xmin;
    fPlottingAxisData[iName]->fMaxValue = xmax;
    
    return fPlottingAxisData.size();
}

VPlottingAxisData* VPlotUtilities::getPlottingAxis( string iName )
{
    if( fPlottingAxisData.find( iName ) != fPlottingAxisData.end() )
    {
        return fPlottingAxisData[iName];
    }
    
    return 0;
}

unsigned int VPlotUtilities::listPlottingAxis()
{
    map< string, VPlottingAxisData* >::iterator i_iter;
    
    for( i_iter = fPlottingAxisData.begin(); i_iter != fPlottingAxisData.end(); i_iter++ )
    {
        cout << ( *i_iter ).first << ":\t\t";
        if( ( *i_iter ).second )
        {
            cout << " [" << ( *i_iter ).second->fMinValue;
            cout << ", " << ( *i_iter ).second->fMaxValue;
            cout << "] " << ( *i_iter ).second->fLogAxis;
            cout << "\t " << ( *i_iter ).second->fAxis;
            cout << "title: " << ( *i_iter ).second->fAxisTitle;
        }
        cout << endl;
    }
    
    return fPlottingAxisData.size();
}

void VPlotUtilities::plotHistogramTitle( TH1* h )
{
    if( !h )
    {
        return;
    }
    
    if( fPlotHistogramTitle.size() > 0 )
    {
        TText* iT = new TText( -4., 1.e3, fPlotHistogramTitle.c_str() );
        iT->SetNDC();
        if( fPlotHistogramTitle_x > 0. && fPlotHistogramTitle_y > 0. )
        {
            iT->SetX( fPlotHistogramTitle_x );
            iT->SetY( fPlotHistogramTitle_y );
        }
        else
        {
            iT->SetX( 0.2 );
            iT->SetY( 0.8 );
        }
        iT->Draw();
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

VPlottingAxisData::VPlottingAxisData()
{
    fName = "";
    fAxisTitle = "";
    fAxis = "";
    fLogAxis = false;
    fMinValue = 0.;
    fMaxValue = 0.;
}
