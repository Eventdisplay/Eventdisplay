/*! \class VPlotCompareDataWithMC
 *
 * compare data with MC
 *
 */

#include "VPlotCompareDataWithMC.h"

VPlotCompareDataWithMC::VPlotCompareDataWithMC( string iFile )
{
    setDebug();
    setPosterPlot();
    setNTel();
    setPrintName();
    setRelativePlotRange();
    
    fDataFileName = iFile;
    fDataFile = 0;
    if( fDataFileName.size() > 0 )
    {
        openDataFile( fDataFileName );
    }
    
    gStyle->SetPadGridX( 0 );
    gStyle->SetPadGridY( 0 );
    gStyle->SetPalette( 1 );
}

void VPlotCompareDataWithMC::help()
{
    cout << endl;
    cout << "compare image and shower parameter distribution of simulations and on/off data" << endl;
    cout << "------------------------------------------------------------------------------" << endl;
	cout << endl << endl;
    cout << "shower parameter distributions:  stereo_parameter()  " << endl << endl;
    cout << "mscw/mscl energy dependent:      msc_vs_energy_plots()  " << endl << endl;
    cout << "mwr/mlr energy dependent:        mwr_vs_energy_plots()  " << endl << endl;
    cout << "multiplicity plots:              multiplicity_plots() " << endl << endl;
    cout << "emission height:                 emission_height()" << endl << endl;
    cout << "core plots:                      core_plots()" << endl << endl;
    cout << "core distance plots:             distance_plots()" << endl << endl;
    cout << "centroid plots:                  centroids()" << endl << endl;
    cout << "image parameter distributions:   single_telescope()" << endl << endl;
    cout << "width/length energy dependent:   widthlength_vs_energy_plots()" << endl << endl;
    cout << "erecratio energy dependent:      erecRatio_vs_energy_plots()" << endl << endl;
    cout << "mva energy dependent:            mva_vs_energy_plots()" << endl << endl;
    cout << "mva parameter distribution:      mva_parameter()" << endl << endl;
    cout << endl;
}


void VPlotCompareDataWithMC::setAxisTitles( TH2D* h, string iS, int iTel )
{
    if( !h )
    {
        return;
    }
    char htit[200];
    sprintf( htit, "%s - Telescope %d (%s)", h->GetXaxis()->GetTitle(), iTel, iS.c_str() );
    h->SetXTitle( htit );
    sprintf( htit, "%s - Telescope %d (%s)", h->GetYaxis()->GetTitle(), iTel, iS.c_str() );
    h->SetYTitle( htit );
}


bool VPlotCompareDataWithMC::openDataFile( string ifile )
{
    fDataFileName = ifile;
    fDataFile = new TFile( fDataFileName.c_str() );
    if( fDataFile->IsZombie() )
    {
        cout << "VPlotCompareDataWithMC::openDataFile() error opening " << fDataFileName << endl;
        return false;
    }
    
    return true;
}

void VPlotCompareDataWithMC::drawMatchingTests( TH1D* h1, TH1D* h2, double x_min, double x_max )
{
    if( !h1 || !h2 )
    {
        return;
    }
    
    if( x_min < -998. )
    {
        x_min = h1->GetXaxis()->GetXmin();
    }
    if( x_max < -998. )
    {
        x_max = h1->GetXaxis()->GetXmax();
    }
    
    // calculate matching of distributions
    double KSProb = h1->KolmogorovTest( h2 );
    double KSSig  = TMath::ErfInverse( 1. - KSProb ) * TMath::Sqrt( 2. );
    double Chi2Prob = h2->Chi2Test( h1, "WW" );
    double Chi2Sig = TMath::ErfInverse( 1. - Chi2Prob ) * TMath::Sqrt( 2. );
    
    char hname[200];
    sprintf( hname, "mean (MC): %.2f#pm %.2f, mean (data): %.2f#pm %.2f", h1->GetMean(), h1->GetRMS(), h2->GetMean(), h2->GetRMS() );
    TLatex* iM = new TLatex( x_min + 0.1 * ( x_max - x_min ), 0.84 * h1->GetMaximum(), hname );
    iM->SetTextSize( iM->GetTextSize() * 0.6 );
    iM->Draw();
    sprintf( hname, "KS-test | P = %1.2e (%1.1f #sigma)", KSProb, KSSig );
    TLatex* iK = new TLatex( x_min + 0.1 * ( x_max - x_min ), 0.78 * h1->GetMaximum(), hname );
    iK->SetTextSize( iK->GetTextSize() * 0.6 );
    iK->Draw();
    sprintf( hname, "Chi2 | P = %1.2e (%1.1f #sigma)", Chi2Prob, Chi2Sig );
    TLatex* iC = new TLatex( x_min + 0.1 * ( x_max - x_min ), 0.72 * h1->GetMaximum(), hname );
    iC->SetTextSize( iC->GetTextSize() * 0.6 );
    iC->Draw();
}


void VPlotCompareDataWithMC::plotLegend( TH1D* hsims, TH1D* hdiff, double x0 )
{
    if( !hsims || !hdiff )
    {
        return;
    }
    
    TLegend* iLegend = new TLegend( x0, 0.68, 0.85 + ( x0 - 0.5 ), 0.85 );
    iLegend->AddEntry( hsims, "simulations", "pl" );
    iLegend->AddEntry( hdiff, "On-Off", "pl" );
    iLegend->Draw();
}

/*!
  get scale factor between simulations and data

  bContents = 1:   scale to same contents
  bContents = 2:   scale to same maximum value
  bContents = 3:   scale to same maximum value (three bins around maximum)
*/
void VPlotCompareDataWithMC::getScaling( TH1D* h_sims, TH1D* h_diff, double& s_sims, double& s_diff,
        int bContents, double xmin, double xmax )
{
    if( !h_sims || !h_diff )
    {
        return;
    }
    double z = 0.;
    ////////////////////////////////////
    // scale to same contents (integral)
    if( bContents == 1 )
    {
        int i_min = 1;
        int i_max = h_sims->GetNbinsX();
        if( xmin > -9998 )
        {
            i_min = h_sims->GetXaxis()->FindBin( xmin );
        }
        if( xmax <  9998 )
        {
            i_max = h_sims->GetXaxis()->FindBin( xmax );
        }
        for( int i = i_min; i <= i_max; i++ )
        {
            z += h_sims->GetBinContent( i );
        }
        s_sims = z;
        z = 0;
        i_min = 1;
        i_max = h_diff->GetNbinsX();
        if( xmin > -9998 )
        {
            i_min = h_diff->GetXaxis()->FindBin( xmin );
        }
        if( xmax <  9998 )
        {
            i_max = h_diff->GetXaxis()->FindBin( xmax );
        }
        for( int i = i_min; i <= i_max; i++ )
        {
            if( h_diff->GetBinContent( i ) > 0 )
            {
                z += h_diff->GetBinContent( i );
            }
        }
        
        s_diff = 1.;
        cout << " Bin Content:  data: " << z << "\t sims: " << s_sims << endl;
        if( s_sims > 0. )
        {
            s_sims = z / s_sims;
        }
    }
    //////////////////////////////////
    // scale to same maximum
    else if( bContents == 2 )
    {
        s_sims = h_sims->GetMaximum();
        z      = h_diff->GetMaximum();
        cout << h_sims->GetName() << " Maximum : data:" << z << "\t sims: " << s_sims << endl;
        if( s_sims > 0. )
        {
            s_sims = z / s_sims;
        }
        s_diff = 1.;
    }
    //////////////////////////////////
    // scale to peak (three bins around maximum)
    else if( bContents == 3 )
    {
        int imaxbin = h_sims->GetMaximumBin();
        s_sims = h_sims->GetBinContent( imaxbin );
        if( imaxbin > 1 )
        {
            s_sims += h_sims->GetBinContent( imaxbin - 1 );
        }
        if( imaxbin < h_sims->GetNbinsX() )
        {
            s_sims += h_sims->GetBinContent( imaxbin + 1 );
        }
        
        imaxbin = h_diff->GetMaximumBin();
        z = h_diff->GetBinContent( imaxbin );
        if( imaxbin > 1 )
        {
            z += h_diff->GetBinContent( imaxbin - 1 );
        }
        if( imaxbin < h_diff->GetNbinsX() )
        {
            z += h_diff->GetBinContent( imaxbin + 1 );
        }
        if( s_sims > 0. )
        {
            s_sims = z / s_sims;
        }
        s_diff = 1.;
        cout << h_sims->GetName() << " scale to 3 bins around maximum " << endl;
    }
    
    // make sure that results are positiv
    if( s_sims < 0. )
    {
        s_sims *= -1.;
        s_diff *= -1.;
    }
    cout << "Scaling: SIMS " << s_sims << "\t" << z << "\t DIFF " << s_diff << endl;
}

void VPlotCompareDataWithMC::getScaling( double& s_sims, double& s_diff, string his,
        int bContents, double xmin, double xmax )
{
    if( fDataFile == 0 )
    {
        cout << "NO SCALING POSSIBLE, no file" << endl;
        s_sims = 1.;
        s_diff = 1.;
        return;
    }
    cout << "scale on histograms " << his << ", scale to";
    if( bContents == 1 )
    {
        cout << " histogram contents" << endl;
    }
    else if( bContents == 2 )
    {
        cout << " histogram maximum" << endl;
    }
    else if( bContents == 3 )
    {
        cout << " histogram maximum (peak)" << endl;
    }
    
    char hname[200];
    sprintf( hname, "h%s_SIMS", his.c_str() );
    TH1D* h_sims = ( TH1D* )fDataFile->Get( hname );
    sprintf( hname, "h%s_DIFF", his.c_str() );
    TH1D* h_diff = ( TH1D* )fDataFile->Get( hname );
    if( !h_sims || !h_diff )
    {
        cout << "NO SCALING POSSIBLE, no histograms " << h_sims << " " << h_diff << "\t( " << hname << ")" << endl;
        s_sims = 1.;
        s_diff = 1.;
        return;
    }
    
    getScaling( h_sims, h_diff, s_sims, s_diff, bContents, xmin, xmax );
}


void VPlotCompareDataWithMC::setHistogramAtt( TH2D* his, double imin )
{
    if( !his )
    {
        return;
    }
    
    if( imin > -900. )
    {
        his->SetMinimum( imin );
    }
    his->SetStats( 0 );
    his->GetYaxis()->SetTitleOffset( 1.3 );
}

void VPlotCompareDataWithMC::setHistogramAtt( TH1D* his, int icolor, double iwidth, double isize, int imarker, int irebin, double iTitleOffset )
{
    his->SetLineColor( icolor );
    his->SetMarkerColor( icolor );
    his->SetLineWidth( ( Width_t )iwidth );
    his->SetMarkerSize( isize );
    if( imarker != 0 )
    {
        his->SetMarkerStyle( imarker );
    }
    his->SetStats( 0 );
    his->GetYaxis()->SetTitleOffset( iTitleOffset );
    his->GetXaxis()->SetTitleSize( 0.05 );
    his->GetYaxis()->SetTitleSize( 0.05 );
    if( irebin != 1 )
    {
        his->Rebin( irebin );
    }
}

TF1* VPlotCompareDataWithMC::do_theta2Fit( TH1D* h, int icolor, int istyle )
{
    if( !h )
    {
        return 0;
    }
    
    cout << "now fitting " << h->GetName() << endl;
    
    char hname[200];
    sprintf( hname, "fT2_%s", h->GetName() );
    TF1* fTheta2 = new TF1( hname, "[0]*([1]/2./[2]/[2]*TMath::Exp(-x/2./[2]/[2])+(1-[1])/2./[3]/[3]*TMath::Exp(-x/2./[3]/[3]))", 0., 0.3 );
    fTheta2->SetParameter( 0, 5. );
    fTheta2->SetParameter( 1, 0.5 );
    fTheta2->SetParLimits( 1, 0., 1. );
    fTheta2->SetParameter( 2, 0.03 );
    fTheta2->SetParLimits( 2, 0., 1. );
    fTheta2->SetParameter( 3, 0.09 );
    fTheta2->SetParLimits( 3, 0., 1. );
    fTheta2->SetLineColor( icolor );
    fTheta2->SetLineStyle( istyle );
    
    h->Fit( fTheta2, "REM0" );
    //  cout << "CHI2 " << fTheta2->GetChisquare()/htheta2_diff->GetXaxis()->GetXmax() << endl;
    
    cout << endl;
    
    return fTheta2;
}

void VPlotCompareDataWithMC::plotCummulativePlot( TH1D* h1, TH1D* h2, double xmin, double xmax, double iSystematicCutCheck, int iTelescope, bool iLeftToRight, double iBinValue )
{
    if( !h1 || !h2 )
    {
        return;
    }
    
    if( xmin < -998. )
    {
        xmin = h1->GetXaxis()->GetXmin();
    }
    if( xmax < -998. )
    {
        xmax = h1->GetXaxis()->GetXmax();
    }
    TH1D* hCumu_1 = VHistogramUtilities::get_Cumulative_Histogram( h1, true, iLeftToRight, iBinValue );
    TH1D* hCumu_2 = VHistogramUtilities::get_Cumulative_Histogram( h2, true, iLeftToRight, iBinValue );
    if( !hCumu_1 || !hCumu_2 )
    {
        return;
    }
    
    hCumu_1->SetMaximum( 1.15 );
    hCumu_1->SetMinimum( 0. );
    hCumu_1->SetAxisRange( xmin, xmax );
    hCumu_1->SetYTitle( "cumulative distribution" );
    hCumu_1->SetLineWidth( 2 );
    hCumu_2->SetLineWidth( 2 );
    
    if( iTelescope > 0 )
    {
        char hname[500];
        sprintf( hname, "%s (T%d)", hCumu_1->GetXaxis()->GetTitle(), iTelescope );
        hCumu_1->SetXTitle( hname );
    }
    
    hCumu_1->Draw();
    hCumu_2->Draw( "sames" );
    
    TLine* iL = new TLine( xmin, 1., xmax, 1. );
    iL->SetLineStyle( 2 );
    iL->Draw();
    
    // calculate systematics a certain value
    if( iSystematicCutCheck > xmin && iSystematicCutCheck < xmax )
    {
        TLine* iLS = new TLine( iSystematicCutCheck, 0., iSystematicCutCheck, 1.15 );
        iLS->Draw();
        double iSys = TMath::Abs( hCumu_1->GetBinContent( hCumu_1->FindBin( iSystematicCutCheck ) )
                                  - hCumu_2->GetBinContent( hCumu_2->FindBin( iSystematicCutCheck ) ) );
        char hname[200];
        sprintf( hname, "Sys.in flux estimates: %.1f%%", iSys * 100. );
        TLatex* iM = new TLatex( xmin + 0.1 * ( xmax - xmin ), 0.88 * hCumu_1->GetMaximum(), hname );
        iM->Draw();
    }
}

void VPlotCompareDataWithMC::plotRelativePlot( TH1D* h1, TH1D* h2, double xmin, double xmax, int iTelescope )
{
    if( !h1 || !h2 )
    {
        return;
    }
    
    if( xmin < -998. )
    {
        xmin = h1->GetXaxis()->GetXmin();
    }
    if( xmax < -998. )
    {
        xmax = h1->GetXaxis()->GetXmax();
    }
    
    char hname[200];
    sprintf( hname, "%s_rel", h1->GetName() );
    TH1D* hRel = ( TH1D* )h1->Clone( hname );
    hRel->Divide( h2 );
    
    setHistogramAtt( hRel, 1, 1, 1, 20, 1 );
    hRel->SetAxisRange( xmin, xmax );
    hRel->SetYTitle( "sims/data" );
    hRel->SetMinimum( fRelatePlotRange_min );
    hRel->SetMaximum( fRelatePlotRange_max );
    
    if( iTelescope > 0 )
    {
        char hname[500];
        sprintf( hname, "%s (T%d)", hRel->GetXaxis()->GetTitle(), iTelescope );
        hRel->SetXTitle( hname );
    }
    
    hRel->Draw();
    
    TLine* iL = new TLine( xmin, 1., xmax, 1. );
    iL->SetLineStyle( 2 );
    iL->Draw();
}

TCanvas* VPlotCompareDataWithMC::plotRelativePlots( char* i_CanvasName, char* i_CanvasTitle, TH1D* h1, TH1D* h2, double xmin, double xmax )
{
    if( !h1 || !h2 )
    {
        return 0;
    }
    
    char hname[600];
    char htitle[600];
    sprintf( hname, "rel_%s", i_CanvasName );
    sprintf( htitle, "%s rel. plot", i_CanvasTitle );
    TCanvas* cRel = new TCanvas( hname, htitle, 200, 200, 400, 400 );
    cRel->SetGridx( 0 );
    cRel->SetGridy( 0 );
    cRel->SetLeftMargin( 0.13 );
    
    sprintf( hname, "rel_%s", h1->GetName() );
    TH1D* hR = ( TH1D* )h1->Clone( hname );
    setHistogramAtt( hR, 1, 3, 1, 20, 1 );
    hR->Divide( h2 );
    hR->SetMinimum( fRelatePlotRange_min );
    hR->SetMaximum( fRelatePlotRange_max );
    hR->GetYaxis()->SetTitle( "sims / data" );
    hR->SetAxisRange( xmin, xmax );
    
    hR->Draw();
    
    TLine* iL = new TLine( xmin, 1., xmax, 1. );
    iL->SetLineStyle( 2 );
    iL->Draw();
    
    return cRel;
}

/*
 * compare telescope multiplicities
 *
 */
void VPlotCompareDataWithMC::multiplicity_plots()
{
    if( !fDataFile )
    {
        return;
    }
    
    // get the scaling between simulations and data
    double s_sims = 1.;
    double s_diff = 1.;
    getScaling( s_sims, s_diff, "NImages", 1 );
    
    char hname[600];
    char htitle[600];
    
    // canvases
    sprintf( hname, "cTrigger_%s", fDataFileName.c_str() );
    sprintf( htitle, "multiplicity plots (%s)", fDataFileName.c_str() );
    TCanvas* cTrigger = new TCanvas( hname, htitle, 10, 10, 800, 400 );
    cTrigger->SetGridx( 0 );
    cTrigger->SetGridy( 0 );
    cTrigger->Divide( 2, 1 );
    
    sprintf( hname, "cTriggerRel_%s", fDataFileName.c_str() );
    sprintf( htitle, "multiplicity plots (relative dist., %s)", fDataFileName.c_str() );
    TCanvas* cTriggerRel = new TCanvas( hname, htitle, 410, 10, 800, 400 );
    cTriggerRel->SetGridx( 0 );
    cTriggerRel->SetGridy( 0 );
    cTriggerRel->Divide( 2, 1 );
    
    TH1D* hNImages_SIMS = ( TH1D* )fDataFile->Get( "hNImages_SIMS" );
    TH1D* hNImages_DIFF = ( TH1D* )fDataFile->Get( "hNImages_DIFF" );
    TH1D* hImgSel_SIMS = ( TH1D* )fDataFile->Get( "hImgSel_SIMS" );
    TH1D* hImgSel_DIFF = ( TH1D* )fDataFile->Get( "hImgSel_DIFF" );
    
    if( !hNImages_SIMS || !hNImages_DIFF || !hImgSel_SIMS || !hImgSel_DIFF )
    {
        cout << hNImages_SIMS << "\t" << hNImages_DIFF << "\t" << hImgSel_SIMS << "\t" << hImgSel_DIFF << endl;
        return;
    }
    
    setHistogramAtt( hNImages_SIMS, 2, 3, 1, 20, 1 );
    setHistogramAtt( hNImages_DIFF, 1, 3, 1, 21, 1 );
    if( hNImages_SIMS->GetEntries() > 0 )
    {
        hNImages_SIMS->Scale( s_sims );
    }
    if( hNImages_DIFF->GetEntries() > 0 )
    {
        hNImages_DIFF->Scale( s_diff );
    }
    
    cTrigger->cd( 1 );
    gPad->SetGridx( 0 );
    gPad->SetGridy( 0 );
    gPad->SetLeftMargin( 0.13 );
    hNImages_SIMS->SetMaximum( hNImages_SIMS->GetMaximum() * 1.2 );
    hNImages_SIMS->Draw();
    hNImages_SIMS->SetYTitle( "number of showers [a.u.]" );
    hNImages_DIFF->Draw( "same" );
    
    setHistogramAtt( hImgSel_SIMS, 2, 3, 1, 20, 1 );
    setHistogramAtt( hImgSel_DIFF, 1, 3, 1, 21, 1 );
    if( hImgSel_SIMS->GetEntries() > 0 )
    {
        hImgSel_SIMS->Scale( s_sims );
    }
    if( hImgSel_DIFF->GetEntries() > 0 )
    {
        hImgSel_DIFF->Scale( s_diff );
    }
    
    cTrigger->cd( 2 );
    gPad->SetGridx( 0 );
    gPad->SetGridy( 0 );
    gPad->SetLeftMargin( 0.13 );
    hImgSel_SIMS->SetMaximum( hImgSel_SIMS->GetMaximum() * 1.2 );
    hImgSel_SIMS->Draw();
    hImgSel_SIMS->SetYTitle( "number of showers [a.u.]" );
    hImgSel_DIFF->Draw( "same" );
    //
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        char hname[200];
        sprintf( hname, "%s-Trigger.pdf", fPrintName.c_str() );
        cTrigger->Print( hname );
    }
    
    
    // relative plots
    cTriggerRel->cd( 1 );
    gPad->SetLeftMargin( 0.13 );
    plotRelativePlot( hNImages_SIMS, hNImages_DIFF );
    
    cTriggerRel->cd( 2 );
    gPad->SetLeftMargin( 0.13 );
    plotRelativePlot( hImgSel_SIMS, hImgSel_DIFF );
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        char hname[200];
        sprintf( hname, "%s-TriggerRel.pdf", fPrintName.c_str() );
        cTriggerRel->Print( hname );
    }
    
}

/*

  plot three single canvases for a certain histogram

*/
TCanvas* VPlotCompareDataWithMC::plot_singleCanvas( string iHistoName, string iCanvasTitle, double iHistoXAxisMax, string iScalingVariable )
{
    if( !fDataFile )
    {
        return 0;
    }
    
    char hname[600];
    char htitle[600];
    
    // get the scaling between simulations and data
    double s_sims = 1.;
    double s_diff = 1.;
    getScaling( s_sims, s_diff, iScalingVariable, 1, -9999., iHistoXAxisMax );
    
    sprintf( hname, "c%s_%s", iHistoName.c_str(), fDataFileName.c_str() );
    sprintf( htitle, "%s (%s)", iCanvasTitle.c_str(), fDataFileName.c_str() );
    TCanvas* cEMH = new TCanvas( hname, htitle, 610, 410, 900, 600 );
    cEMH->SetGridx( 0 );
    cEMH->SetGridy( 0 );
    
    sprintf( hname, "%s_SIMS", iHistoName.c_str() );
    TH1D* hHistogram_SIMS = ( TH1D* )fDataFile->Get( hname );
    sprintf( hname, "%s_DIFF", iHistoName.c_str() );
    TH1D* hHistogram_DIFF = ( TH1D* )fDataFile->Get( hname );
    
    if( !hHistogram_SIMS || !hHistogram_DIFF )
    {
        return 0;
    }
    
    setHistogramAtt( hHistogram_SIMS, 2, 3, 1, 20, 1 );
    setHistogramAtt( hHistogram_DIFF, 1, 3, 1, 21, 1 );
    cout << "SCALE " << s_sims << "\t" << s_diff << endl;
    if( hHistogram_SIMS->GetEntries() > 0 )
    {
        hHistogram_SIMS->Scale( s_sims );
    }
    if( hHistogram_DIFF->GetEntries() > 0 && TMath::Abs( s_sims - 1. ) > 1.e-4 )
    {
        hHistogram_DIFF->Scale( s_diff );
    }
    
    hHistogram_SIMS->SetMaximum( hHistogram_SIMS->GetMaximum() * 1.4 );
    hHistogram_SIMS->SetAxisRange( 0., iHistoXAxisMax );
    hHistogram_SIMS->Draw();
    hHistogram_DIFF->Draw( "same" );
    
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        sprintf( hname, "%s-%s.pdf", fPrintName.c_str(), iHistoName.c_str() );
        cEMH->Print( hname );
    }
    
    TCanvas* c = plotRelativePlots( hname, htitle, hHistogram_SIMS, hHistogram_DIFF, 0., iHistoXAxisMax );
    if( c && fPrintName.size() > 0 )
    {
        sprintf( hname, "%s-%sRel.pdf", fPrintName.c_str(), iHistoName.c_str() );
        c->Print( hname );
    }
    return cEMH;
}

TCanvas* VPlotCompareDataWithMC::emission_height( double iEmissionHeightMax )
{
    return plot_singleCanvas( "hEmissionHeight", "emission height", iEmissionHeightMax, "EmissionHeight" );
}

void VPlotCompareDataWithMC::widthlength_vs_energy_plots( int iTelescope, int iRebin, double xmin, double xmax )
{
    if( !fDataFile )
    {
        return;
    }
    
    plot_energyDependentDistributions( "length", iRebin, xmin, xmax , "REL", iTelescope );
    plot_energyDependentDistributions( "width", iRebin, xmin, xmax , "REL", iTelescope );
    plot_energyDependentDistributions( "length", iRebin, xmin, xmax, "SIMSDIFF", iTelescope );
    plot_energyDependentDistributions( "width", iRebin, xmin, xmax, "SIMSDIFF", iTelescope );
    
}


/*

    plot mscw and mscl energy dependent

*/
void VPlotCompareDataWithMC::msc_vs_energy_plots( int iRebin, double xmin, double xmax, double iSystematicCutCheck )
{
    if( !fDataFile )
    {
        return;
    }
    
    plot_energyDependentDistributions( "MSCW", iRebin, xmin, xmax , "CUMU", 0, iSystematicCutCheck );
    plot_energyDependentDistributions( "MSCL", iRebin, xmin, xmax , "CUMU", 0, iSystematicCutCheck );
    plot_energyDependentDistributions( "MSCL", iRebin, xmin, xmax , "REL" );
    plot_energyDependentDistributions( "MSCW", iRebin, xmin, xmax , "REL" );
    plot_energyDependentDistributions( "MSCL", iRebin, xmin, xmax );
    plot_energyDependentDistributions( "MSCW", iRebin, xmin, xmax );
}

/*

    plot mva energy dependent

*/
void VPlotCompareDataWithMC::mva_vs_energy_plots( int iRebin, double xmin, double xmax, double iSystematicCutCheck )
{
    if( !fDataFile )
    {
        return;
    }
    
    plot_energyDependentDistributions( "MVA", iRebin, xmin, xmax , "CUMU", 0, iSystematicCutCheck );
    plot_energyDependentDistributions( "MVA", iRebin, xmin, xmax , "REL" );
    plot_energyDependentDistributions( "MVA", iRebin, xmin, xmax );
}

/*
 *
 * plot erecratio vs energy
 *
 */
void VPlotCompareDataWithMC::erecRatio_vs_energy_plots( int iTelescope, int iRebin, double xmin, double xmax )
{
    if( !fDataFile )
    {
        return;
    }
    
    plot_energyDependentDistributions( "erecratio", iRebin, xmin, xmax, "REL", iTelescope );
    plot_energyDependentDistributions( "erecratio", iRebin, xmin, xmax, "SIMSDIFF", iTelescope );
}

/*
 * energy dependent plots
 *
 */
TCanvas* VPlotCompareDataWithMC::plot_energyDependentDistributions( string iVariable, int iRebin, double x_min, double x_max,
        string iPlot, int iTelescope, double iSystematicCutCheck,
        string iXVariable, double y_min, bool iPlotLogY )
{
    if( !fDataFile )
    {
        return 0;
    }
    
    double KSProb = 0;
    double KSSig = 0;
    double Chi2Prob = 0.;
    double Chi2Sig = 0.;
    
    // get the scaling between simulations and data
    double s_sims = 1.;
    double s_diff = 1.;
    
    double error_sims = 0.;
    double error_diff = 0.;
    
    char hname[600];
    char htitle[600];
    if( iTelescope > 0 )
    {
        sprintf( hname, "c_%s_%s_%s_%s_%d", iVariable.c_str(), iXVariable.c_str(), fDataFile->GetName(), iPlot.c_str(), iTelescope );
        sprintf( htitle, "%s vs %s T%d (%s) %s", iVariable.c_str(), iXVariable.c_str(), iTelescope, fDataFile->GetName(), iPlot.c_str() );
    }
    else
    {
        sprintf( hname, "c_%s_%s_%s_%s", iVariable.c_str(), iXVariable.c_str(), fDataFile->GetName(), iPlot.c_str() );
        sprintf( htitle, "%s vs %s (%s) %s", iVariable.c_str(), iXVariable.c_str(), fDataFile->GetName(), iPlot.c_str() );
    }
    TCanvas* c_MS = new TCanvas( hname, htitle, 100, 10, 900, 600 );
    c_MS->SetGridx( 0 );
    c_MS->SetGridy( 0 );
    c_MS->Divide( 3, 2 );
    
    if( iTelescope > 0 )
    {
        sprintf( hname, "h%s%s_%d_SIMS", iVariable.c_str(), iXVariable.c_str(), iTelescope );
    }
    else
    {
        sprintf( hname, "h%s%s_SIMS", iVariable.c_str(), iXVariable.c_str() );
    }
    TH2D* h_sims = ( TH2D* )fDataFile->Get( hname );
    if( iTelescope > 0 )
    {
        sprintf( hname, "h%s%s_%d_DIFF", iVariable.c_str(), iXVariable.c_str(), iTelescope );
    }
    else
    {
        sprintf( hname, "h%s%s_DIFF", iVariable.c_str(), iXVariable.c_str() );
    }
    TH2D* h_diff = ( TH2D* )fDataFile->Get( hname );
    
    if( !h_sims || !h_diff )
    {
        return 0;
    }
    
    // loop over all bins in energy
    for( int i = 1; i <= h_sims->GetXaxis()->GetNbins(); i++ )
    {
        if( iTelescope > 0 )
        {
            sprintf( hname, "h%s%s_%d_SIMS_%d", iVariable.c_str(), iXVariable.c_str(), iTelescope, i );
        }
        else
        {
            sprintf( hname, "h%s%s_SIMS_%d", iVariable.c_str(), iXVariable.c_str(), i );
        }
        TH1D* hSims = h_sims->ProjectionY( hname, i, i );
        setHistogramAtt( hSims, 2, 1, 1, 20, iRebin );
        if( iTelescope > 0 )
        {
            sprintf( hname, "h_%s%s_%d_diff_%d", iVariable.c_str(), iXVariable.c_str(), iTelescope, i );
        }
        else
        {
            sprintf( hname, "h_%s%s_diff_%d", iVariable.c_str(), iXVariable.c_str(), i );
        }
        TH1D* hDiff = h_diff->ProjectionY( hname, i, i );
        setHistogramAtt( hDiff, 1, 1, 1, 21, iRebin );
        
        getScaling( hSims, hDiff, s_sims, s_diff, 3 );
        if( hSims->GetEntries() > 0 )
        {
            hSims->Scale( s_sims );
        }
        if( hSims->GetEntries() > 0 )
        {
            hDiff->Scale( s_diff );
        }
        
        for( int j = 0; j < hSims->GetNbinsX(); j++ )
        {
            error_sims = hSims->GetBinError( j ) * s_sims;
            hSims->SetBinError( j, error_sims );
            error_diff = hDiff->GetBinError( j ) * s_diff;
            hDiff->SetBinError( j, error_diff );
        }
        
        hSims->SetAxisRange( x_min, x_max );
        hDiff->SetAxisRange( x_min, x_max );
        if( !iPlotLogY )
        {
            hSims->SetMaximum( hSims->GetMaximum() * 1.8 );
        }
        else
        {
            hSims->SetMaximum( TMath::Exp( TMath::Log( hSims->GetMaximum() ) * 1.8 ) );
        }
        if( y_min > 0. )
        {
            hSims->SetMinimum( y_min );
        }
        else
        {
            hSims->SetMinimum( 0. );
        }
        
        // calculate matching of distributions
        KSProb = hSims->KolmogorovTest( hDiff );
        KSSig  = TMath::ErfInverse( 1. - KSProb ) * TMath::Sqrt( 2. );
        Chi2Prob = hDiff->Chi2Test( hSims, "WW" );
        Chi2Sig = TMath::ErfInverse( 1. - Chi2Prob ) * TMath::Sqrt( 2. );
        
        // draw histograms
        TPad* g = ( TPad* )c_MS->cd( i );
        if( g )
        {
            g->SetRightMargin( 0.01 );
            g->SetLeftMargin( 0.13 );
            g->SetBottomMargin( 0.13 );
            if( iPlotLogY )
            {
                g->SetLogy( true );
            }
        }
        if( iPlot == "SIMSDIFF" )
        {
            if( iTelescope > 0 )
            {
                sprintf( hname, "%s (T%d)", hSims->GetXaxis()->GetTitle(), iTelescope );
                hSims->SetXTitle( hname );
            }
            hSims->Draw();
            hDiff->Draw( "same" );
            
            float x = 0.;
            if( iVariable.find( "MW" ) != string::npos || iVariable.find( "ML" ) != string::npos )
            {
                x = 1.;
            }
            TLine* lLine = new TLine( x, hSims->GetMinimum(), x, hSims->GetMaximum() );
            lLine->SetLineStyle( 2 );
            lLine->Draw();
        }
        else if( iPlot == "REL" )
        {
            plotRelativePlot( hSims, hDiff, x_min, x_max, iTelescope );
        }
        else if( iPlot == "CUMU" )
        {
            if( iVariable == "MVA" )
            {
                plotCummulativePlot( hSims, hDiff, x_min, x_max, iSystematicCutCheck, iTelescope, false, -0.95 );
            }
            else
            {
                plotCummulativePlot( hSims, hDiff, x_min, x_max, iSystematicCutCheck, iTelescope );
            }
        }
        
        
        if( iXVariable == "Erec" )
        {
            sprintf( hname, "%.1f < log_{10} E_{rec} < %.1f", h_sims->GetXaxis()->GetBinLowEdge( i ),  h_sims->GetXaxis()->GetBinUpEdge( i ) );
        }
        else if( iXVariable == "ntubes" )
        {
            sprintf( hname, "%.1f < ntubes < %.1f", h_sims->GetXaxis()->GetBinLowEdge( i ),  h_sims->GetXaxis()->GetBinUpEdge( i ) );
        }
        else if( iXVariable == "size" )
        {
            sprintf( hname, "%.1f < log_{10} size < %.1f", h_sims->GetXaxis()->GetBinLowEdge( i ),  h_sims->GetXaxis()->GetBinUpEdge( i ) );
        }
        else if( iXVariable == "sizeHG" )
        {
            sprintf( hname, "%.1f < log_{10} sizeHG < %.1f", h_sims->GetXaxis()->GetBinLowEdge( i ),  h_sims->GetXaxis()->GetBinUpEdge( i ) );
        }
        else if( iXVariable == "sizeLG" )
        {
            sprintf( hname, "%.1f < log_{10} sizeLG < %.1f", h_sims->GetXaxis()->GetBinLowEdge( i ),  h_sims->GetXaxis()->GetBinUpEdge( i ) );
        }
        
        
        TLatex* iT = new TLatex;
        iT->SetTextSize( iT->GetTextSize() * 0.6 );
        iT->DrawLatexNDC( 0.22, 0.84, hname );
        
        
        //don't print the rest of the stuff if we are plotting relative or cumulative distributions.
        if( iPlot == "SIMSDIFF" )
        {
            sprintf( hname, "mean (MC): %.2f#pm %.2f, mean (data): %.2f#pm %.2f", hSims->GetMean(), hSims->GetRMS(), hDiff->GetMean(), hDiff->GetRMS() );
            iT->DrawLatexNDC( 0.22, 0.80, hname );
            sprintf( hname, "KS-test | P = %1.2e (%1.1f #sigma)", KSProb, KSSig );
            iT->DrawLatexNDC( 0.22, 0.76, hname );
            sprintf( hname, "Chi2 | P = %1.2e (%1.1f #sigma)", Chi2Prob, Chi2Sig );
            iT->DrawLatexNDC( 0.22, 0.72, hname );
        }
        
    }
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        if( iTelescope > 0 )
        {
            sprintf( hname, "%s-%s-%s-T%d-%s.pdf", fPrintName.c_str(), iXVariable.c_str(), iVariable.c_str(), iTelescope, iPlot.c_str() );
        }
        else
        {
            sprintf( hname, "%s-%s-%s-%s.pdf", fPrintName.c_str(), iXVariable.c_str(), iVariable.c_str(), iPlot.c_str() );
        }
        c_MS->Print( hname );
    }
    
    return c_MS;
}

/*
 *
 * plot stereo parameter
 *
 */
TCanvas* VPlotCompareDataWithMC::stereo_parameter()
{
    if( !fDataFile )
    {
        return 0;
    }
    
    // two canvases, one with sims and diff, one with on/off histograms
    
    char hname[600];
    char htitle[600];
    TCanvas* cOO = new TCanvas( "cOO", "relative plots (stereo parameter)", 100, 10, 900, 600 );
    cOO->SetGridx( 0 );
    cOO->SetGridy( 0 );
    cOO->Divide( 2, 2 );
    
    sprintf( hname, "cSD_%s", fDataFileName.c_str() );
    sprintf( htitle, "sims/diff (stereo parameter, %s)", fDataFileName.c_str() );
    TCanvas* cSD = new TCanvas( hname, htitle, 10, 10, 900, 600 );
    cSD->SetGridx( 0 );
    cSD->SetGridy( 0 );
    cSD->Divide( 2, 2 );
    
    // theta2 < 0.02
    //
    
    TH1D* ht2_sims = ( TH1D* )fDataFile->Get( "htheta2_SIMS" );
    setHistogramAtt( ht2_sims, 2, 1, 1, 20, 1 );
    ht2_sims->SetYTitle( "number of shower [a.u.]" );
    
    TH1D* ht2_diff = ( TH1D* )fDataFile->Get( "htheta2_DIFF" );
    setHistogramAtt( ht2_diff, 1, 1, 1, 25, 1 );
    
    
    // get the scaling between simulations and data
    double s_sims = 1.;
    double s_diff = 1.;
    getScaling( s_sims, s_diff, "theta2", 1 );
    
    ht2_sims->Scale( s_sims );
    ht2_diff->Scale( s_diff );
    
    // plot everything
    
    cSD->cd( 1 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    ht2_sims->GetYaxis()->SetTitleOffset( 1.5 );
    ht2_sims->SetAxisRange( 0., 0.05 );
    ht2_sims->Draw();
    ht2_diff->Draw( "same" );
    
    drawMatchingTests( ht2_sims, ht2_diff, 0., 0.05 );
    
    cOO->cd( 1 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    plotRelativePlot( ht2_sims, ht2_diff, 0., 0.05 );
    
    // theta2 < 0.02
    //
    
    TH1D* hlt2_sims = ( TH1D* )fDataFile->Get( "hltheta2_SIMS" );
    setHistogramAtt( hlt2_sims, 2, 3, 1, 20, 1 );
    hlt2_sims->SetYTitle( "number of shower [a.u.]" );
    
    TH1D* hlt2_diff = ( TH1D* )fDataFile->Get( "hltheta2_DIFF" );
    setHistogramAtt( hlt2_diff, 1, 3, 1, 25, 1 );
    
    getScaling( s_sims, s_diff, "ltheta2", 2 );
    hlt2_sims->Scale( s_sims );
    hlt2_diff->Scale( s_diff );
    
    hlt2_sims->SetAxisRange( -5., -1. );
    
    cSD->cd( 2 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    hlt2_sims->GetYaxis()->SetTitleOffset( 1.5 );
    hlt2_sims->SetMaximum( hlt2_sims->GetMaximum() * 1.3 );
    hlt2_sims->Draw();
    hlt2_diff->Draw( "same" );
    
    drawMatchingTests( hlt2_sims, hlt2_diff, -5., 1. );
    
    cOO->cd( 2 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    plotRelativePlot( hlt2_sims, hlt2_diff, -5., 1. );
    
    // MSCW
    //
    
    TH1D* hmscw_sims = ( TH1D* )fDataFile->Get( "hMSCW_SIMS" );
    setHistogramAtt( hmscw_sims, 2, 1, 1, 20, 2 );
    hmscw_sims->SetYTitle( "number of shower [a.u.]" );
    
    TH1D* hmscw_diff = ( TH1D* )fDataFile->Get( "hMSCW_DIFF" );
    setHistogramAtt( hmscw_diff, 1, 1, 1, 25, 2 );
    
    hmscw_sims->SetAxisRange( -1., 1. );
    getScaling( s_sims, s_diff, "MSCW", 2, -0.5, 0.5 );
    if( hmscw_sims->GetEntries() > 0 )
    {
        hmscw_sims->Scale( s_sims );
    }
    if( hmscw_diff->GetEntries() > 0 )
    {
        hmscw_diff->Scale( s_diff );
    }
    
    cSD->cd( 3 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    hmscw_sims->GetYaxis()->SetTitleOffset( 1.5 );
    hmscw_sims->SetMaximum( hmscw_sims->GetMaximum() * 1.5 );
    hmscw_sims->Draw();
    hmscw_diff->Draw( "same" );
    TLine* lmscw = new TLine( 0., hmscw_sims->GetMinimum(), 0., hmscw_sims->GetMaximum() );
    lmscw->SetLineStyle( 2 );
    lmscw->Draw();
    
    drawMatchingTests( hmscw_sims, hmscw_diff, -1., 1. );
    
    cOO->cd( 3 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    plotRelativePlot( hmscw_sims, hmscw_diff, -1., 1. );
    
    // MSCL
    //
    
    TH1D* hmscl_sims = ( TH1D* )fDataFile->Get( "hMSCL_SIMS" );
    setHistogramAtt( hmscl_sims, 2, 3, 1, 21, 2 );
    hmscl_sims->SetYTitle( "number of shower [a.u.]" );
    
    TH1D* hmscl_diff = ( TH1D* )fDataFile->Get( "hMSCL_DIFF" );
    setHistogramAtt( hmscl_diff, 1, 3, 1, 25, 2 );
    hmscl_diff->SetLineWidth( 3 );
    hmscl_diff->SetStats( 0 );
    
    hmscl_sims->SetAxisRange( -1., 1. );
    getScaling( s_sims, s_diff, "MSCL", 1, -0.75, 0.75 );
    if( hmscl_sims->GetEntries() > 0 )
    {
        hmscl_sims->Scale( s_sims );
    }
    if( hmscl_diff->GetEntries() > 0 )
    {
        hmscl_diff->Scale( s_diff );
    }
    
    cSD->cd( 4 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    hmscl_sims->GetYaxis()->SetTitleOffset( 1.5 );
    hmscl_sims->SetMaximum( hmscl_sims->GetMaximum() * 1.5 );
    hmscl_sims->Draw();
    hmscl_diff->Draw( "same" );
    TLine* lmscl = new TLine( 0., hmscl_sims->GetMinimum(), 0., hmscl_sims->GetMaximum() );
    lmscl->SetLineStyle( 2 );
    lmscl->Draw();
    drawMatchingTests( hmscl_sims, hmscl_diff, -1., 1. );
    
    cOO->cd( 4 );
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    plotRelativePlot( hmscl_sims, hmscl_diff, -1., 1. );
    
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        char hn[200];
        sprintf( hn, "%s-Stereo.pdf", fPrintName.c_str() );
        cSD->Print( hn );
        sprintf( hn, "%s-StereoRel.pdf", fPrintName.c_str() );
        cOO->Print( hn );
    }
    return cSD;
}

/*
*
* plot mva parameter
*
*/

void VPlotCompareDataWithMC::mva_parameter()
{
    if( !fDataFile )
    {
        return;
    }
    
    char hname[600];
    char htitle[600];
    TCanvas* cOO = new TCanvas( "cOO", "relative plots (mva parameter)", 100, 10, 600, 600 );
    cOO->SetGridx( 0 );
    cOO->SetGridy( 0 );
    
    sprintf( hname, "cSD_%s", fDataFileName.c_str() );
    sprintf( htitle, "sims/diff (mva parameter, %s)", fDataFileName.c_str() );
    TCanvas* cSD = new TCanvas( hname, htitle, 10, 10, 600, 600 );
    cSD->SetGridx( 0 );
    cSD->SetGridy( 0 );
    
    TH1D* hmva_sims = ( TH1D* )fDataFile->Get( "hMVA_SIMS" );
    setHistogramAtt( hmva_sims, 2, 3, 1, 20, 2 );
    hmva_sims->SetYTitle( "number of events per bin" );
    
    TH1D* hmva_diff = ( TH1D* )fDataFile->Get( "hMVA_DIFF" );
    setHistogramAtt( hmva_diff, 1, 3, 1, 21, 2 );
    hmva_diff->SetLineWidth( 3 );
    hmva_diff->SetStats( 0 );
    
    // scaling of the histograms
    double s_sims = 1.;
    double s_diff = 1.;
    
    double error_sims = 0.;
    double error_diff = 0.;
    
    hmva_sims->SetAxisRange( -1., 1. );
    getScaling( s_sims, s_diff, "MVA", 1, -0.75, 0.75 );
    if( hmva_sims->GetEntries() > 0 )
    {
        hmva_sims->Scale( s_sims );
    }
    if( hmva_diff->GetEntries() > 0 )
    {
        hmva_diff->Scale( s_diff );
    }
    
    for( int j = 0; j < hmva_sims->GetNbinsX(); j++ )
    {
        error_sims = hmva_sims->GetBinError( j ) * s_sims;
        hmva_sims->SetBinError( j, error_sims );
        error_diff = hmva_diff->GetBinError( j ) * s_diff;
        hmva_diff->SetBinError( j, error_diff );
    }
    
    cSD->cd();
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    hmva_sims->GetYaxis()->SetTitleOffset( 1.3 );
    hmva_sims->SetMaximum( hmva_sims->GetMaximum() * 1.5 );
    hmva_sims->Draw();
    hmva_diff->Draw( "same" );
    TLine* lmva = new TLine( 0., hmva_sims->GetMinimum(), 0., hmva_sims->GetMaximum() );
    lmva->SetLineStyle( 2 );
    lmva->Draw();
    drawMatchingTests( hmva_sims, hmva_diff, -1., 1. );
    
    cOO->cd();
    gPad->SetLeftMargin( 0.14 );
    gPad->SetBottomMargin( 0.12 );
    plotRelativePlot( hmva_sims, hmva_diff, -1., 1. );
    
    if( fPrintName.size() > 0 )
    {
        char hn[200];
        sprintf( hn, "%s-Stereo.pdf", fPrintName.c_str() );
        cSD->Print( hn );
        sprintf( hn, "%s-StereoRel.pdf", fPrintName.c_str() );
        cOO->Print( hn );
    }
}

/*
 *
 * core parameter plots
 *
 * horrible code...programmer doesn't know what a loop is
*/

TCanvas* VPlotCompareDataWithMC::core_plots( int iRebin, int iScaling )
{
    if( !fDataFile )
    {
        return 0;
    }
    
    TCanvas* cOCore = new TCanvas( "cOCore", "on/off (core positions)", 100, 10, 1200, 600 );
    cOCore->SetGridx( 0 );
    cOCore->SetGridy( 0 );
    cOCore->Divide( 3, 2 );
    
    TCanvas* cSCore = new TCanvas( "cSCore", "sims & diff (core positions)", 10, 10, 1200, 600 );
    cSCore->SetGridx( 0 );
    cSCore->SetGridy( 0 );
    cSCore->Divide( 3, 2 );
    
    TCanvas* cSCoreRel = new TCanvas( "cSCoreRel", "sims/diff (core positions)", 10, 110, 600, 630 );
    cSCoreRel->SetGridx( 0 );
    cSCoreRel->SetGridy( 0 );
    cSCoreRel->Divide( 2, 2 );
    
    // xcore
    //
    TH1D* hXcore_sims = ( TH1D* )fDataFile->Get( "hXcore_SIMS" );
    setHistogramAtt( hXcore_sims, 2, 1, 1, 20, iRebin );
    
    TH1D* hXcore_on = ( TH1D* )fDataFile->Get( "hXcore_ON" );
    setHistogramAtt( hXcore_on, 3, 1, 1, 20, iRebin );
    
    TH1D* hXcore_off = ( TH1D* )fDataFile->Get( "hXcore_OFF" );
    setHistogramAtt( hXcore_off, 4, 1, 1, 20, iRebin );
    
    TH1D* hXcore_diff = ( TH1D* )fDataFile->Get( "hXcore_DIFF" );
    setHistogramAtt( hXcore_diff, 1, 1, 1, 21, iRebin );
    
    hXcore_sims->SetAxisRange( -250., 250. );
    hXcore_on->SetAxisRange( -250., 250. );
    
    double nSims = 0.;
    double nDiff = 0.;
    
    getScaling( nSims, nDiff, "Xcore", iScaling );
    hXcore_diff->Scale( nDiff );
    hXcore_sims->Scale( nSims );
    
    cSCore->cd( 1 );
    gPad->SetLeftMargin( 0.13 );
    hXcore_sims->SetMaximum( hXcore_sims->GetMaximum() * 1.8 );
    hXcore_sims->Draw();
    hXcore_diff->Draw( "same" );
    plotLegend( hXcore_sims, hXcore_diff, 0.15 );
    
    cSCoreRel->cd( 1 );
    gPad->SetLeftMargin( 0.13 );
    plotRelativePlot( hXcore_sims, hXcore_diff, -250., 250. );
    
    cOCore->cd( 1 );
    gPad->SetLeftMargin( 0.13 );
    hXcore_on->Draw();
    hXcore_off->Draw( "same" );
    
    // Ycore
    //
    
    TH1D* hYcore_sims = ( TH1D* )fDataFile->Get( "hYcore_SIMS" );
    setHistogramAtt( hYcore_sims, 2, 1, 1, 20, iRebin );
    hYcore_sims->SetMaximum( hYcore_sims->GetMaximum() * 1.5 );
    
    TH1D* hYcore_on = ( TH1D* )fDataFile->Get( "hYcore_ON" );
    setHistogramAtt( hYcore_on, 3, 1, 1, 20, iRebin );
    
    TH1D* hYcore_off = ( TH1D* )fDataFile->Get( "hYcore_OFF" );
    setHistogramAtt( hYcore_off, 4, 1, 1, 20, iRebin );
    
    TH1D* hYcore_diff = ( TH1D* )fDataFile->Get( "hYcore_DIFF" );
    setHistogramAtt( hYcore_diff, 1, 1, 1, 21, iRebin );
    
    hYcore_sims->SetAxisRange( -250., 250. );
    hYcore_on->SetAxisRange( -250., 250. );
    
    getScaling( nSims, nDiff, "Ycore", iScaling );
    hYcore_diff->Scale( nDiff );
    hYcore_sims->Scale( nSims );
    
    cSCore->cd( 2 );
    hYcore_sims->SetMaximum( hYcore_sims->GetMaximum() * 1.8 );
    hYcore_sims->Draw();
    hYcore_diff->Draw( "same" );
    plotLegend( hYcore_sims, hYcore_diff, 0.13 );
    
    cSCoreRel->cd( 2 );
    plotRelativePlot( hYcore_sims, hYcore_diff, -250., 250. );
    
    cOCore->cd( 2 );
    hYcore_on->Draw();
    hYcore_off->Draw( "same" );
    
    // az & elevation
    cSCore->cd( 3 );
    TH1D* hArrayEl_sims = ( TH1D* )fDataFile->Get( "hArrayEl_SIMS" );
    TH1D* hArrayEl_on = ( TH1D* )fDataFile->Get( "hArrayEl_ON" );
    TH1D* hArrayEl_off = ( TH1D* )fDataFile->Get( "hArrayEl_OFF" );
    TH1D* hArrayEl_diff = ( TH1D* )fDataFile->Get( "hArrayEl_DIFF" );
    if( hArrayEl_sims && hArrayEl_on && hArrayEl_off && hArrayEl_diff )
    {
        setHistogramAtt( hArrayEl_sims, 2, 1, 1, 20, iRebin );
        setHistogramAtt( hArrayEl_on, 3, 1, 1, 20, iRebin );
        setHistogramAtt( hArrayEl_off, 4, 1, 1, 20, iRebin );
        setHistogramAtt( hArrayEl_diff, 1, 1, 1, 21, iRebin );
        hArrayEl_sims->SetAxisRange( 40., 90. );
        hArrayEl_on->SetAxisRange( 40., 90. );
        getScaling( nSims, nDiff, "ArrayEl", iScaling );
        hArrayEl_diff->Scale( nDiff );
        hArrayEl_sims->Scale( nSims );
        
        cSCore->cd( 3 );
        hArrayEl_sims->Draw();
        hArrayEl_diff->Draw( "same" );
        plotLegend( hArrayEl_sims, hArrayEl_diff, 0.13 );
        
        cSCoreRel->cd( 3 );
        plotRelativePlot( hArrayEl_sims, hArrayEl_diff, 40., 90. );
        
        cOCore->cd( 3 );
        hArrayEl_off->Draw();
        hArrayEl_on->Draw( "same" );
    }
    cSCore->cd( 6 );
    TH1D* hArrayAz_sims = ( TH1D* )fDataFile->Get( "hArrayAz_SIMS" );
    TH1D* hArrayAz_on = ( TH1D* )fDataFile->Get( "hArrayAz_ON" );
    TH1D* hArrayAz_off = ( TH1D* )fDataFile->Get( "hArrayAz_OFF" );
    TH1D* hArrayAz_diff = ( TH1D* )fDataFile->Get( "hArrayAz_DIFF" );
    if( hArrayAz_sims && hArrayAz_on && hArrayAz_off && hArrayAz_diff )
    {
        setHistogramAtt( hArrayAz_sims, 2, 1, 1, 20, 1 );
        setHistogramAtt( hArrayAz_on, 3, 1, 1, 20, 1 );
        setHistogramAtt( hArrayAz_off, 4, 1, 1, 20, 1 );
        setHistogramAtt( hArrayAz_diff, 1, 1, 1, 21, 1 );
        getScaling( nSims, nDiff, "ArrayAz", iScaling );
        hArrayAz_diff->Scale( nDiff );
        hArrayAz_sims->Scale( nSims );
        
        cSCore->cd( 6 );
        hArrayAz_sims->SetMaximum( hArrayAz_sims->GetMaximum() * 1.8 );
        hArrayAz_sims->Draw();
        hArrayAz_diff->Draw( "same" );
        plotLegend( hArrayAz_sims, hArrayAz_diff, 0.13 );
        
        cSCoreRel->cd( 4 );
        plotRelativePlot( hArrayAz_sims, hArrayAz_diff, 0., 360. );
        
        cOCore->cd( 6 );
        hArrayAz_off->Draw();
        hArrayAz_on->Draw( "same" );
    }
    
    // //////////////////////
    // XY plot
    
    cOCore->cd( 4 );
    TH2D* hXYcore_on = ( TH2D* )fDataFile->Get( "hXYcore_ON" );
    setHistogramAtt( hXYcore_on, 1. );
    hXYcore_on->SetXTitle( "core position X (ON) [m]" );
    hXYcore_on->SetYTitle( "core position Y (ON) [m]" );
    
    hXYcore_on->Draw( "colz" );
    
    cOCore->cd( 5 );
    TH2D* hXYcore_off = ( TH2D* )fDataFile->Get( "hXYcore_OFF" );
    setHistogramAtt( hXYcore_off, 1. );
    hXYcore_off->SetXTitle( "core position X (OFF) [m]" );
    hXYcore_off->SetYTitle( "core position Y (OFF) [m]" );
    
    hXYcore_off->Draw( "colz" );
    
    cSCore->cd( 4 );
    TH2D* hXYcore_diff = ( TH2D* )fDataFile->Get( "hXYcore_DIFF" );
    setHistogramAtt( hXYcore_diff, 1. );
    hXYcore_diff->SetXTitle( "core position X (ON-OFF) [m]" );
    hXYcore_diff->SetYTitle( "core position Y (ON-OFF) [m]" );
    
    hXYcore_diff->Draw( "colz" );
    
    cSCore->cd( 5 );
    TH2D* hXYcore_sims = ( TH2D* )fDataFile->Get( "hXYcore_SIMS" );
    setHistogramAtt( hXYcore_sims, 1. );
    
    hXYcore_sims->SetXTitle( "core position X (SIMS) [m]" );
    hXYcore_sims->SetYTitle( "core position Y (SIMS) [m]" );
    hXYcore_sims->Draw( "colz" );
    
    
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        char hn[200];
        sprintf( hn, "%s-Core.pdf", fPrintName.c_str() );
        cSCore->Print( hn );
        sprintf( hn, "%s-OnOffCore.pdf", fPrintName.c_str() );
        cOCore->Print( hn );
        sprintf( hn, "%s-CoreRel.pdf", fPrintName.c_str() );
        cSCoreRel->Print( hn );
    }
    return cSCore;
}

void VPlotCompareDataWithMC::centroids()
{

    if( !fDataFile )
    {
        return;
    }
    
    TCanvas* cOCentro = new TCanvas( "cOCentro", "on/off (centroids plots)", 100, 10, 600, fNTel * 300 );
    cOCentro->SetGridx( 0 );
    cOCentro->SetGridy( 0 );
    cOCentro->Divide( 2, fNTel );
    
    TCanvas* cSCentro = new TCanvas( "cSCentro", "sims/diff (centroids plots)", 10, 10, 600, fNTel * 300 );
    cSCentro->SetGridx( 0 );
    cSCentro->SetGridy( 0 );
    cSCentro->Divide( 2, fNTel );
    
    int iC = 1;
    
    TH2D* hCenXY_sims[200];
    TH2D* hCenXY_diff[200];
    
    TH2D* hCenXY_off[200];
    TH2D* hCenXY_on[200];
    
    char hname[200];
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        sprintf( hname, "hcen_xy%u_SIMS", i + 1 );
        hCenXY_sims[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hCenXY_sims[i], -999. );
        setAxisTitles( hCenXY_sims[i], "sims", i + 1 );
        
        sprintf( hname, "hcen_xy%u_DIFF", i + 1 );
        hCenXY_diff[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hCenXY_diff[i], 0.01 );
        setAxisTitles( hCenXY_diff[i], "on-off", i + 1 );
        
        sprintf( hname, "hcen_xy%u_ON", i + 1 );
        hCenXY_on[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hCenXY_on[i], -999. );
        setAxisTitles( hCenXY_on[i], "on", i + 1 );
        
        sprintf( hname, "hcen_xy%u_OFF", i + 1 );
        hCenXY_off[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hCenXY_off[i], -999. );
        setAxisTitles( hCenXY_off[i], "off", i + 1 );
        
        cSCentro->cd( iC );
        hCenXY_sims[i]->Draw( "contz" );
        
        cOCentro->cd( iC );
        hCenXY_on[i]->Draw( "contz" );
        
        iC++;
        cSCentro->cd( iC );
        hCenXY_diff[i]->Draw( "contz" );
        
        cOCentro->cd( iC );
        hCenXY_off[i]->Draw( "contz" );
        
        iC++;
    }
}

/*

  plot distance from telescope to shower core

*/
TCanvas* VPlotCompareDataWithMC::distance_plots()
{
    if( !fDataFile )
    {
        return 0;
    }
    
    TCanvas* cODist = new TCanvas( "cODist", "on/off (distance plots)", 100, 10, 900, fNTel * 150 );
    cODist->SetGridx( 0 );
    cODist->SetGridy( 0 );
    cODist->Divide( 3, fNTel );
    
    TCanvas* cSDist = new TCanvas( "cSDist", "sims/diff (distance plots)", 10, 10, 900, fNTel * 150 );
    cSDist->SetGridx( 0 );
    cSDist->SetGridy( 0 );
    cSDist->Divide( 4, fNTel );
    
    char hname[200];
    char htitle[200];
    
    TCanvas* cdistsim[100];
    for( int i = 0; i < 100; i++ )
    {
        cdistsim[i] = 0;
    }
    if( fNTel > 100 )
    {
        cout << "too many telescopes ..." << endl;
        return 0;
    }
    if( fPlotPoster )
    {
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            sprintf( hname, "cdistsim_%d", i );
            sprintf( htitle, "distance plot (Telescope %d)", i + 1 );
            cdistsim[i] = new TCanvas( hname, htitle, 10 + i * 200, 10 + i * 200, 400, 400 );
            cdistsim[i]->SetGridx( 0 );
            cdistsim[i]->SetGridy( 0 );
        }
    }
    
    // get histogram scaling
    double s_sims = 1.;
    double s_diff = 1.;
    
    int iCsi = 1;
    int iCoo = 1;
    
    TH1D* hR_sims[200];
    TH1D* hR_diff[200];
    TH1D* hR_on[200];
    TH1D* hR_off[200];
    
    TH2D* hdistR_sims[200];
    TH2D* hdistR_diff[200];
    TH2D* hdistR_on[200];
    TH2D* hdistR_off[200];
    
    TH1D* hrel = 0;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        // R
        sprintf( hname, "hr_%u_SIMS", i + 1 );
        hR_sims[i] = ( TH1D* )fDataFile->Get( hname );
        setHistogramAtt( hR_sims[i], 2, 1, 0.5, 20, 1 );
        hR_sims[i]->SetMaximum( hR_sims[i]->GetMaximum() * 1.3 );
        hR_sims[i]->SetYTitle( "number of shower [a.u.]" );
        
        sprintf( hname, "hr_%u_DIFF", i + 1 );
        hR_diff[i] = ( TH1D* )fDataFile->Get( hname );
        setHistogramAtt( hR_diff[i], 1, 1, 0.5, 21, 1 );
        
        sprintf( hname, "hr_%u_ON", i + 1 );
        hR_on[i] = ( TH1D* )fDataFile->Get( hname );
        setHistogramAtt( hR_on[i], 3, 1, 0.5, 20, 1 );
        
        sprintf( hname, "hr_%u_OFF", i + 1 );
        hR_off[i] = ( TH1D* )fDataFile->Get( hname );
        setHistogramAtt( hR_off[i], 4, 1, 0.5, 21, 1 );
        
        
        sprintf( hname, "r_%u", i + 1 );
        getScaling( s_sims, s_diff, hname, 1 );
        if( hR_sims[i]->GetEntries() > 0 )
        {
            hR_sims[i]->Scale( s_sims );
        }
        if( hR_diff[i]->GetEntries() > 0 )
        {
            hR_diff[i]->Scale( s_diff );
        }
        
        cSDist->cd( iCsi );
        gPad->SetLeftMargin( 0.13 );
        hR_sims[i]->Draw( "e" );
        hR_diff[i]->Draw( "e same" );
        
        cODist->cd( iCoo );
        gPad->SetLeftMargin( 0.13 );
        hR_on[i]->Draw( "e" );
        hR_off[i]->Draw( "e same" );
        
        iCsi++;
        iCoo++;
        
        // relative plots
        if( hR_sims[i] && hR_diff[i] )
        {
            sprintf( hname, "hR_RE_%u", i );
            hrel = ( TH1D* )hR_sims[i]->Clone( hname );
            hrel->Divide( hR_diff[i] );
            hrel->SetYTitle( "sims/data" );
            hrel->SetMinimum( fRelatePlotRange_min );
            hrel->SetMaximum( fRelatePlotRange_max );
            setHistogramAtt( hrel, 1, 1, 0.5, 21, 1 );
            cSDist->cd( iCsi );
            gPad->SetLeftMargin( 0.13 );
            hrel->Draw( "e" );
            TLine* iLine = new TLine( hrel->GetXaxis()->GetXmin(), 1., hrel->GetXaxis()->GetXmax(), 1. );
            iLine->SetLineStyle( 2 );
            iLine->Draw();
        }
        iCsi++;
        
        // distR
        //
        sprintf( hname, "hdistR%u_SIMS", i + 1 );
        hdistR_sims[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hdistR_sims[i], -999. );
        setAxisTitles( hdistR_sims[i], "sims", i + 1 );
        hdistR_sims[i]->SetAxisRange( 0., 1.5, "Y" );
        
        sprintf( hname, "hdistR%u_DIFF", i + 1 );
        hdistR_diff[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hdistR_diff[i], 0.001 );
        setAxisTitles( hdistR_diff[i], "on-off", i + 1 );
        hdistR_diff[i]->SetAxisRange( 0., 1.5, "Y" );
        
        sprintf( hname, "hdistR%u_ON", i + 1 );
        hdistR_on[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hdistR_on[i], -999. );
        setAxisTitles( hdistR_on[i], "on", i + 1 );
        hdistR_on[i]->SetAxisRange( 0., 1.5, "Y" );
        
        sprintf( hname, "hdistR%u_OFF", i + 1 );
        hdistR_off[i] = ( TH2D* )fDataFile->Get( hname );
        setHistogramAtt( hdistR_off[i], -999. );
        setAxisTitles( hdistR_off[i], "off", i + 1 );
        hdistR_off[i]->SetAxisRange( 0., 1.5, "Y" );
        
        // draw everything
        cSDist->cd( iCsi );
        hdistR_sims[i]->Draw( "contz" );
        
        cODist->cd( iCoo );
        hdistR_on[i]->Draw( "contz" );
        
        iCsi++;
        iCoo++;
        
        cSDist->cd( iCsi );
        hdistR_diff[i]->Draw( "contz" );
        
        cODist->cd( iCoo );
        hdistR_off[i]->Draw( "contz" );
        
        iCsi++;
        iCoo++;
    }
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        char hn[200];
        sprintf( hn, "%s-Dist.pdf", fPrintName.c_str() );
        cSDist->Print( hn );
        sprintf( hn, "%s-DistRel.pdf", fPrintName.c_str() );
        cODist->Print( hn );
    }
    return cSDist;
}

/*
 * compare single telescope (image) parameters
 *
 * T1 = 1, T2 = 2
 *
 * plot is SIMSDIFF, ONOFF, REL
 *
 */

TCanvas* VPlotCompareDataWithMC::single_telescope( int telid )
{
    TCanvas* c = 0;
    if( telid > 0 )
    {
        c = single_telescope( telid, "REL" );
        single_telescope( telid, "SIMSDIFF" );
    }
    else
    {
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            single_telescope( i + 1, "REL" );
            single_telescope( i + 1, "SIMSDIFF" );
        }
    }
    return c;
}

TCanvas* VPlotCompareDataWithMC::single_telescope( int telid, string iPlot, bool iOneCanvas,
        int iScalingMethod, int i_rebin )
{

    if( iPlot != "SIMSDIFF" && iPlot != "ONOFF" && iPlot != "REL" )
    {
        cout << "error: unknown plotting mode (allowed are SIMSDIFF, ONOFF, REL)" << endl;
        return 0;
    }
    if( !fDataFile )
    {
        return 0;
    }
    
    double KSProb = 0;
    double KSSig = 0;
    char text[1000];
    // scaling factor
    double s_sims = 1.;
    double s_diff = 1.;
    char htitle[600];
    sprintf( htitle, "width_%d", telid );
    getScaling( s_sims, s_diff, htitle, iScalingMethod );
    
    //////////////////////////////////////
    // histogram names to be plotted
    vector< string > hname;
    vector< int >    f_rebin;
    vector< int >   f_logy;
    vector< double > f_x_min;
    vector< double > f_x_max;
    hname.push_back( "width" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 0.30 );
    hname.push_back( "length" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 0.50 );
    hname.push_back( "dist" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 2.10 );
    hname.push_back( "size" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 1 );
    f_x_min.push_back( 1.5 );
    f_x_max.push_back( 6.50 );
    hname.push_back( "sizeHG" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 1 );
    f_x_min.push_back( 1.5 );
    f_x_max.push_back( 6.50 );
    hname.push_back( "sizeLG" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 1 );
    f_x_min.push_back( 1.5 );
    f_x_max.push_back( 6.5 );
    hname.push_back( "fraclow" );
    f_rebin.push_back( 1 );
    f_logy.push_back( 1 );
    f_x_min.push_back( 0.0 );
    f_x_max.push_back( 1. );
    hname.push_back( "nlowgain" );
    f_rebin.push_back( 1 );
    f_logy.push_back( 1 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 40. );
    hname.push_back( "los" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 40. );
    hname.push_back( "asym" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( -2.0 );
    f_x_max.push_back( 2.0 );
    hname.push_back( "cen_x" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( -2.0 );
    f_x_max.push_back( 2.0 );
    hname.push_back( "cen_y" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( -2.0 );
    f_x_max.push_back( 2.0 );
    hname.push_back( "ntubes" );
    f_rebin.push_back( 1 );
    f_logy.push_back( 1 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 90. );
    hname.push_back( "mwrt" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 0.5 );
    f_x_max.push_back( 1.5 );
    hname.push_back( "mltt" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 0.5 );
    f_x_max.push_back( 1.5 );
    hname.push_back( "loss" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 1 );
    f_x_min.push_back( 0. );
    f_x_max.push_back( 0.25 );
    hname.push_back( "tgrad_x" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( -7.5 );
    f_x_max.push_back( 7.5 );
    hname.push_back( "pedvarT" );
    f_rebin.push_back( i_rebin );
    f_logy.push_back( 0 );
    f_x_min.push_back( 5. );
    f_x_max.push_back( 10. );
    
    // loop over all histograms and plot them
    char hn[600];
    char cn[600];
    char ct[600];
    
    TH1D* hsims = 0;
    TH1D* hdiff = 0;
    TH1D* hon = 0;
    TH1D* hoff = 0;
    TH1D* hrel = 0;
    
    TCanvas* hc = 0;
    // canvas for all in one
    if( iOneCanvas )
    {
        sprintf( cn, "image parameter comparision (telescope %d, file %s, %s)", telid, fDataFileName.c_str(), iPlot.c_str() );
        sprintf( ct, "cimage_%d_%s_%s", telid, iPlot.c_str(), fDataFileName.c_str() );
        hc = new TCanvas( ct, cn, 10, 10, 1300, 800 );
        hc->SetGridx( 0 );
        hc->SetGridy( 0 );
        hc->Divide( 6, 3 );
    }
    TLegend* iL = 0;
    
    /////////////////////////////////////////////////
    // loop over all histograms and plot them
    for( unsigned int j = 0; j < hname.size(); j++ )
    {
        sprintf( hn, "h%s_%d_SIMS", hname[j].c_str(), telid );
        hsims = ( TH1D* )fDataFile->Get( hn );
        if( !hsims )
        {
            cout << "sims histogram not found " << hn << endl;
            continue;
        }
        sprintf( hn, "h%s_%d_DIFF", hname[j].c_str(), telid );
        hdiff = ( TH1D* )fDataFile->Get( hn );
        if( !hdiff )
        {
            cout << "diff histogram not found " << hn << endl;
            continue;
        }
        sprintf( hn, "h%s_%d_ON", hname[j].c_str(), telid );
        hon = ( TH1D* )fDataFile->Get( hn );
        if( !hon )
        {
            cout << "on histogram not found " << hn << endl;
            continue;
        }
        sprintf( hn, "h%s_%d_OFF", hname[j].c_str(), telid );
        hoff = ( TH1D* )fDataFile->Get( hn );
        if( !hoff )
        {
            cout << "off histogram not found " << hn << endl;
            continue;
        }
        // rebin histograms
        hsims->Rebin( f_rebin[j] );
        hdiff->Rebin( f_rebin[j] );
        hon->Rebin( f_rebin[j] );
        hoff->Rebin( f_rebin[j] );
        sprintf( htitle, "%s_%d", hname[j].c_str(), telid );
        getScaling( s_sims, s_diff, htitle, iScalingMethod );
        // normalize sims histograms to data histograms
        hsims->Scale( s_sims );
        if( TMath::Abs( s_diff - 1. ) > 1.e-2 )
        {
            hdiff->Scale( s_diff );
        }
        // relative histograms
        sprintf( hn, "h%s_%d_RE", hname[j].c_str(), telid );
        hrel = ( TH1D* )hsims->Clone( hn );
        hrel->Divide( hdiff );
        
        if( iOneCanvas )
        {
            TPad* g = ( TPad* )hc->cd( j + 1 );
            g->SetGridx( 0 );
            g->SetGridy( 0 );
            g->SetRightMargin( 0.01 );
            g->SetLeftMargin( 0.13 );
            g->SetBottomMargin( 0.13 );
            if( iPlot != "REL" )
            {
                g->SetLogy( f_logy[j] );
            }
        }
        else
        {
            sprintf( cn, "c%d%s%s", telid, iPlot.c_str(), hname[j].c_str() );
            sprintf( ct, "%s (telescope %d, %s)", hname[j].c_str(), telid, iPlot.c_str() );
            int xmax = 600;
            if( fPlotPoster )
            {
                xmax = 400;
            }
            hc = new TCanvas( cn, ct, 10, 100, xmax, 400 );
            hc->SetGridx( 0 );
            hc->SetGridy( 0 );
            hc->SetLogy( f_logy[j] );
            hc->cd();
        }
        //		iL = new TLegend( 0.80 , 0.80, 1.05, 1.05 );
        
        double iTitleOffset = 1.3;
        if( !iOneCanvas )
        {
            iTitleOffset = 1.;
        }
        
        setHistogramAtt( hsims, 2, 1, 0.5, 20, 1, iTitleOffset );
        if( fPlotPoster )
        {
            setHistogramAtt( hsims, 2, 3, 2, 20, 1, iTitleOffset );
        }
        setHistogramAtt( hdiff, 1, 1, 0.5, 21, 1, iTitleOffset );
        if( fPlotPoster )
        {
            setHistogramAtt( hdiff, 1, 3, 2, 21, 1, iTitleOffset );
        }
        setHistogramAtt( hon, 3, 1, 1, 20, 1 );
        if( fPlotPoster )
        {
            setHistogramAtt( hon, 2, 3, 1, 20, 1 );
        }
        setHistogramAtt( hoff, 4, 1, 1, 21, 1 );
        if( fPlotPoster )
        {
            setHistogramAtt( hoff, 1, 3, 1, 21, 1 );
        }
        setHistogramAtt( hrel, 9, 1, 1, 20, 1 );
        if( fPlotPoster )
        {
            setHistogramAtt( hrel, 1, 3, 1, 21, 1 );
        }
        
        hdiff->SetYTitle( "number of shower [a.u.]" );
        if( !f_logy[j] && hdiff->GetMinimum() < -5. )
        {
            hdiff->SetMinimum( -5. );
        }
        if( f_logy[j] )
        {
            hdiff->SetMaximum( hdiff->GetBinContent( hdiff->GetMaximumBin() ) * 1.5 );
        }
        else
        {
            hdiff->SetMaximum( hdiff->GetBinContent( hdiff->GetMaximumBin() ) * 1.1 );
        }
        hrel->SetYTitle( "sims/data" );
        hrel->SetMinimum( fRelatePlotRange_min );
        hrel->SetMaximum( fRelatePlotRange_max );
        
        if( hdiff->GetXaxis()->GetXmin() > f_x_min[j] )
        {
            f_x_min[j] = hdiff->GetXaxis()->GetXmin();
        }
        if( hdiff->GetXaxis()->GetXmax() < f_x_max[j] )
        {
            f_x_max[j] = hdiff->GetXaxis()->GetXmax();
        }
        hdiff->SetAxisRange( f_x_min[j], f_x_max[j] );
        hsims->SetAxisRange( f_x_min[j], f_x_max[j] );
        hon->SetAxisRange( f_x_min[j], f_x_max[j] );
        hoff->SetAxisRange( f_x_min[j], f_x_max[j] );
        hrel->SetAxisRange( f_x_min[j], f_x_max[j] );
        
        ////////////////////////////////////////////////
        // difference plots
        if( iPlot == "SIMSDIFF" )
        {
            hdiff->Draw( "cle" );
            hsims->Draw( "cle same" );
            if( iL )
            {
                sprintf( cn, "telescope %d", telid );
                iL->AddEntry( hdiff, cn, "pl" );
                sprintf( cn, "simulations" );
                iL->AddEntry( hsims, cn, "pl" );
            }
            
            if( !gPad->GetLogy() )
            {
                TLine* iL0 = new TLine( f_x_min[j], 0., f_x_max[j], 0. );
                iL0->SetLineStyle( 2 );
                iL0->Draw();
            }
        }
        ////////////////////////////////////////////////
        // on/off plots
        else if( iPlot == "ONOFF" )
        {
            hon->Draw( "cle" );
            hoff->Draw( "cle same" );
            if( iL )
            {
                sprintf( cn, "on telescope %d", telid );
                iL->AddEntry( hon, cn, "pl" );
                sprintf( cn, "off telescope %d", telid );
                iL->AddEntry( hoff, cn, "pl" );
            }
        }
        ////////////////////////////////////////////////
        // relative plots
        else if( iPlot == "REL" )
        {
            hrel->Draw( "cle" );
        }
        if( !fPlotPoster && iPlot != "REL" )
        {
            if( iL )
            {
                iL->Draw();
            }
        }
        
        // line for mscwt and msclt histograms
        if( iPlot != "REL" )
        {
            if( hname[j] == "hmwrt" || hname[j] == "mltt" )
            {
                TLine* iLine = new TLine( 1., hdiff->GetMinimum(), 1., hdiff->GetMaximum() );
                iLine->SetLineStyle( 2 );
                iLine->Draw();
            }
        }
        // line at 1 for relative plots
        else
        {
            TLine* iLine = new TLine( f_x_min[j], 1., f_x_max[j], 1. );
            iLine->SetLineStyle( 2 );
            iLine->Draw();
        }
        if( iPlot == "SIMSDIFF" || iPlot == "REL" )
        {
            // calculate probabilities of agreement
            KSProb = hsims->KolmogorovTest( hdiff );
            KSSig = TMath::ErfInverse( 1 - KSProb ) * TMath::Sqrt( 2 );
            if( KSProb != 0 )
            {
                sprintf( text, "%s (TEL%d) | KS P = %1.2e (%1.1f #sigma)", hname[j].c_str(), telid, KSProb, KSSig );
            }
            else
            {
                sprintf( text, "%s (TEL%d) | KS P = %1.2e (#infty #sigma)", hname[j].c_str(), telid, KSProb );
            }
            cout << text << endl;
            
            TLatex* iT = new TLatex();
            iT->SetText( 0.11, 0.92, text );
            iT->SetNDC();
            iT->Draw();
        }
    }
    // print the canvas
    if( fPrintName.size() > 0 )
    {
        sprintf( hn, "%s-SINGLET-%s-T%d.pdf", fPrintName.c_str(), iPlot.c_str(), telid );
        hc->Print( hn );
    }
    return hc;
}

void VPlotCompareDataWithMC::msc_plots( char* offFile, char* helium, char* proton, double xmin, double xmax, string ivar )
{
    char hname[200];
    
    if( !fDataFile )
    {
        return;
    }
    
    TCanvas* cMSCsim = 0;
    sprintf( hname, "c%sSim", ivar.c_str() );
    cMSCsim = new TCanvas( hname, hname, 450, 510, 400, 400 );
    cMSCsim->SetGridx( 0 );
    cMSCsim->SetGridy( 0 );
    
    sprintf( hname, "h%s_SIMS", ivar.c_str() );
    TH1D* hMSC_sims = ( TH1D* )fDataFile->Get( hname );
    setHistogramAtt( hMSC_sims, 2, 1, 1, 24, 1 );
    hMSC_sims->SetYTitle( "number of showers [a.u.]" );
    
    sprintf( hname, "h%s_DIFF", ivar.c_str() );
    TH1D* hMSC_diff = ( TH1D* )fDataFile->Get( hname );
    setHistogramAtt( hMSC_diff, 1, 1, 1, 21, 1 );
    
    sprintf( hname, "h%s_OFF", ivar.c_str() );
    TH1D* hMSC_off = ( TH1D* )fDataFile->Get( hname );
    setHistogramAtt( hMSC_off, 4, 1, 1, 21, 1 );
    hMSC_off->Rebin( 2 );
    
    // get the scaling between simulations and data
    double s_sims = 1.;
    double s_diff = 1.;
    getScaling( s_sims, s_diff, ivar.c_str(), true );
    hMSC_sims->SetAxisRange( xmin, xmax );
    if( hMSC_sims->GetEntries() > 0 )
    {
        hMSC_sims->Scale( s_sims );
    }
    if( hMSC_diff->GetEntries() > 0 )
    {
        hMSC_diff->Scale( s_diff * 1.15 );
    }
    if( hMSC_off->GetMaximum() > 0 )
    {
        hMSC_off->Scale( hMSC_diff->GetMaximum() / hMSC_off->GetMaximum() / 2. );
    }
    
    cMSCsim->cd();
    hMSC_sims->SetMaximum( hMSC_sims->GetMaximum() * 1.3 );
    if( ivar == "MSCW" )
    {
        hMSC_sims->SetXTitle( "mean scaled width" );
    }
    else
    {
        hMSC_sims->SetXTitle( "mean scaled length" );
    }
    hMSC_sims->Draw( "e" );
    TBox* iB = new TBox( -1.5, 0., 0.5,  hMSC_sims->GetMaximum() );
    iB->SetFillColor( 5 );
    iB->SetFillStyle( 1001 );
    iB->SetFillColor( 18 );
    iB->Draw();
    hMSC_sims->Draw( "axis same" );
    hMSC_sims->Draw( "e same" );
    hMSC_diff->Draw( "e same" );
    //   hMSC_off->Draw( "e same" );
    
    // get off histogram from anasum output file
    if( offFile )
    {
        TFile* fOff = new TFile( offFile );
        if( !fOff->IsZombie() )
        {
            fOff->cd( "total/stereo/stereoParameterHistograms/" );
            if( ivar == "MSCW" )
            {
                sprintf( hname, "hmscw_off" );
            }
            else if( ivar == "MSCL" )
            {
                sprintf( hname, "hmscl_off" );
            }
            TH1D* hMSC_off_ana = ( TH1D* )gDirectory->Get( hname );
            if( !hMSC_off_ana )
            {
                cout << "histogram not found: " << hname << endl;
                return;
            }
            
            hMSC_off_ana->Rebin( 2 );
            hMSC_off_ana->SetLineColor( 4 );
            hMSC_off_ana->SetMarkerColor( 4 );
            hMSC_off_ana->SetLineWidth( 2 );
            if( hMSC_off_ana->GetMaximum() > 0 )
            {
                hMSC_off_ana->Scale( hMSC_diff->GetMaximum() / hMSC_off_ana->GetMaximum() / 2. );
            }
            hMSC_off_ana->Draw( "same" );
        }
    }
    if( helium && proton )
    {
        int nbins = 50;
        TH1D* hHelium = 0;
        TH1D* hProton = 0;
        TH1D* hCR = new TH1D( "hCR", "", nbins, -5., 15. );
        hCR->SetLineWidth( 2 );
        hCR->SetLineColor( 3 );
        hCR->SetMarkerColor( 3 );
        hCR->Sumw2();
        hCR->SetStats( 0 );
        TFile* fH = new TFile( helium );
        if( !fH->IsZombie() )
        {
            TTree* t = ( TTree* )gDirectory->Get( "data" );
            hHelium = new TH1D( "hHelium", "", nbins, -5., 15. );
            hHelium->Sumw2();
            sprintf( hname, "%s > -99. &&Chi2>=0&&(Xoff*Xoff+Yoff*Yoff)<0.1", ivar.c_str() );
            t->Project( hHelium->GetName(), ivar.c_str(), hname );
        }
        hHelium->SetStats( 0 );
        TFile* fP = new TFile( proton );
        if( !fP->IsZombie() )
        {
            TTree* t = ( TTree* )gDirectory->Get( "data" );
            hProton = new TH1D( "hProton", "", nbins, -5., 15. );
            hProton->Sumw2();
            sprintf( hname, "%s > -99. &&Chi2>=0&&(Xoff*Xoff+Yoff*Yoff)<0.1", ivar.c_str() );
            t->Project( hProton->GetName(), ivar.c_str(), hname );
        }
        // now add everything up assuming 70% proton and 30% helium
        // N_HE = 2.55e7, N_p = 6.6e7
        double iHEScale = 6.6e7 / 2.55e7 * 0.3;
        cout << "scaling helium histogram by " << iHEScale << endl;
        hCR->Add( hHelium, hProton, iHEScale, 1. );
        
        if( hCR->GetMaximum() > 0. )
        {
            hCR->Scale( hMSC_diff->GetMaximum() / hCR->GetMaximum() / 2. * 1.05 );
        }
        hCR->Draw( "same" );
    }
    TLegend* iLegend = new TLegend( 0.4, 0.60, 0.85, 0.85 );
    iLegend->AddEntry( hMSC_sims, "simulations (#gamma-rays)", "pl" );
    iLegend->AddEntry( hMSC_diff, "On-Off (Crab data)", "pl" );
    //  iLegend->AddEntry( hCR, "simulations (Cosmic-rays)", "l" );
    //  iLegend->AddEntry( hMSC_off_ana, "data (Cosmic-rays)", "l" );
    iLegend->Draw();
}


/*

    plot mean width plots (energy dependent)
*/
void VPlotCompareDataWithMC::mwr_vs_energy_plots( int iRebin, double xmin, double xmax, double iSystematicCutCheck )
{
    if( !fDataFile )
    {
        return;
    }
    
    plot_energyDependentDistributions( "MLR", iRebin, xmin, xmax, "CUMU", 0, iSystematicCutCheck );
    plot_energyDependentDistributions( "MWR", iRebin, xmin, xmax, "CUMU", 0, iSystematicCutCheck );
    plot_energyDependentDistributions( "MLR", iRebin, xmin, xmax, "REL" );
    plot_energyDependentDistributions( "MWR", iRebin, xmin, xmax, "REL" );
    plot_energyDependentDistributions( "MLR", iRebin, xmin, xmax );
    plot_energyDependentDistributions( "MWR", iRebin, xmin, xmax );
    
    return;
}

/*
 * plot everything
 *
 */

void VPlotCompareDataWithMC::plot( string iPrintName )
{

    TCanvas* c = new TCanvas( "cP", "empty page", 10, 10, 600, 400 );
    c->Draw();
    char hname[400];
    sprintf( hname, "%s.pdf(", iPrintName.c_str() );
    c->Print( hname );
    sprintf( hname, "%s.pdf", iPrintName.c_str() );
    TCanvas* cP = 0;
    cP = plot_energyDependentDistributions( "MSCW", 2, -1.5, 1.5, "SIMSDIFF" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MSCW", 1, -1.5, 1.5, "CUMU" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MSCL", 2, -1.5, 1.5, "SIMSDIFF" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MSCL", 1, -1.5, 1.5, "CUMU" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MWR", 2, 0.7, 1.3, "SIMSDIFF" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MWR", 1, 0.7, 1.3, "CUMU" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MWL", 1, 0.7, 1.3, "SIMSDIFF" );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "MWL", 1, 0.7, 1.3, "CUMU" );
    if( cP )
    {
        cP->Print( hname );
    }
    // telescope dependent plots
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = single_telescope( i, "SIMSDIFF" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = single_telescope( i, "REL" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "width", 1, 0.03, 0.16, "SIMSDIFF", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "width", 1, 0.03, 0.16, "SIMSDIFF", i, -99., "size" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "width", 1, 0.03, 0.16, "SIMSDIFF", i, -99., "sizeHG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "width", 1, 0.03, 0.16, "SIMSDIFF", i, -99., "sizeLG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "width", 1, 0.03, 0.16, "SIMSDIFF", i, -99., "ntubes" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "length", 2, 0.05, 0.75, "SIMSDIFF", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "length", 2, 0.05, 0.75, "SIMSDIFF", i, -99., "size" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "length", 2, 0.05, 0.75, "SIMSDIFF", i, -99., "sizeHG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "length", 2, 0.05, 0.75, "SIMSDIFF", i, -99., "sizeLG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "length", 2, 0.05, 0.75, "SIMSDIFF", i, -99., "ntubes" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "mwrt", 1, 0.7, 1.4, "SIMSDIFF", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "mwrt", 1, 0.7, 1.4, "SIMSDIFF", i, -99., "size" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "mwrt", 1, 0.7, 1.4, "SIMSDIFF", i, -99., "sizeHG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "mwrt", 1, 0.7, 1.4, "SIMSDIFF", i, -99., "sizeLG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "mwrt", 1, 0.7, 1.4, "SIMSDIFF", i, -99., "ntubes" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "erecratio", 1, 0.5, 1.5, "SIMSDIFF", i );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "erecratio", 1, 0.5, 1.5, "REL", i );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "erecratio", 1, 0.5, 1.5, "SIMSDIFF", i, -99., "size" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "erecratio", 1, 0.5, 1.5, "SIMSDIFF", i, -99., "sizeHG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "erecratio", 1, 0.5, 1.5, "SIMSDIFF", i, -99., "sizeLG" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "erecratio", 1, 0.5, 1.5, "SIMSDIFF", i, -99., "ntubes" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "size", 1, 2., 5., "SIMSDIFF", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "size", 1, 2., 5., "REL", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "sizeHG", 1, 2., 5., "SIMSDIFF", i, -99., "Erec", 1, true );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "sizeHG", 1, 2., 5., "REL", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "sizeLG", 2, 3., 5., "SIMSDIFF", i, -99., "Erec", 1, true );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "sizeLG", 2, 3., 5., "REL", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "ntubes", 2, 0.0, 90., "SIMSDIFF", i, -99., "Erec", 1, true );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "ntubes", 2, 0.0, 90., "REL", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "nlowgain", 2, 0.0, 40., "SIMSDIFF", i, -99., "Erec", 1, true );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "nlowgain", 2, 0.0, 40., "REL", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "fraclow", 4, 0.0, 1.0, "SIMSDIFF", i, -99., "Erec", 0.2, true );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "fraclow", 4, 0.0, 1.0, "REL", i, -99., "Erec" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    for( unsigned int i = 1; i <= 4; i++ )
    {
        cP = plot_energyDependentDistributions( "r", 4, 0., 400., "SIMSDIFF", i, -99., "core distance" );
        if( cP )
        {
            cP->Print( hname );
        }
    }
    cP = plot_energyDependentDistributions( "ltheta2", 1, -5., -1., "SIMSDIFF", 0 );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = plot_energyDependentDistributions( "theta2", 1, 0., 0.04, "SIMSDIFF", 0 );
    if( cP )
    {
        cP->Print( hname );
    }
    cP = stereo_parameter();
    if( cP )
    {
        cP->Print( hname );
    }
    cP = core_plots();
    if( cP )
    {
        cP->Print( hname );
    }
    cP = distance_plots();
    if( cP )
    {
        cP->Print( hname );
    }
    cP = emission_height();
    if( cP )
    {
        cP->Print( hname );
    }
    
    sprintf( hname, "%s.pdf)", iPrintName.c_str() );
    c->Print( hname );
    
}
