/*  \class VPlotRadialAcceptance
    \brief plot radial acceptance curves

*/

#include "VPlotRadialAcceptance.h"

VPlotRadialAcceptance::VPlotRadialAcceptance( string iFile, int iAz )
{
    fDebug = false;
    
    setName( "radial acceptance" );
    
    setAxisRange();
    
    fAcceptanceFile = 0;
    fAcceptanceHisto = 0;
    fAcceptanceHisto2D_rot = 0;
    fAcceptanceHisto2D_rot_normalized = false;
    fAcceptanceHisto2D_derot = 0;
    fAcceptanceHisto2D_derot_normalized = false;
    fAcceptanceFunction = 0;
    hPhiDist = 0;
    hPhiDistDeRot = 0;
    
    if( iFile.size() > 0 )
    {
        openAcceptanceFile( iFile, iAz );
    }
}

/*
 * open radial acceptance file and read histograms / functions
 *
 */
bool VPlotRadialAcceptance::openAcceptanceFile( string iFile, int iAzBin )
{
    ////////////////////////////////
    // open acceptance file
    fAcceptanceFile = new TFile( iFile.c_str() );
    if( fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::addAcceptanceFile: error adding acceptance file" << endl;
        cout << iFile << endl;
        return false;
    }
    char hname[200];
    if( iAzBin >= 0 )
    {
        sprintf( hname, "az_%d", iAzBin );
        if( !fAcceptanceFile->cd( hname ) )
        {
            cout << "VPlotRadialAcceptance::openAcceptanceFile error, directory for the following az bin not found: " << iAzBin << endl;
            return false;
        }
    }
    ////////////////////////////////
    // read acceptance histograms from file
    //
    // Note that some of the histograms might not be available in the
    // radial acceptance file
    cout << "reading acceptance histograms from " << fAcceptanceFile->GetName() << endl;
    fAcceptanceHisto = ( TH1F* )fAcceptanceFile->Get( "hRadialAcceptance" );
    if( !fAcceptanceHisto )
    {
        fAcceptanceHisto = ( TH1F* )fAcceptanceFile->Get( "hAccZe_0" );
        if( !fAcceptanceHisto )
        {
            cout << "VPlotRadialAcceptance::addAcceptanceFile: error finding acceptance histogram" << endl;
            cout << hname << endl;
            return false;
        }
    }
    // read acceptance fit function from file
    fAcceptanceFunction = ( TF1* )fAcceptanceFile->Get( "fRadialAcceptance" );
    if( !fAcceptanceFunction )
    {
        fAcceptanceFunction = ( TF1* )fAcceptanceFile->Get( "fAccZe_0" );
        if( !fAcceptanceFunction )
        {
            cout << "VPlotRadialAcceptance::addAcceptanceFile: error finding acceptance fit function" << endl;
            cout << hname << endl;
            return false;
        }
    }
    sprintf( hname, "hXYAccTotDeRot" );
    fAcceptanceHisto2D_derot_normalized = false;
    fAcceptanceHisto2D_derot = ( TH2F* )fAcceptanceFile->Get( hname );
    
    sprintf( hname, "hXYaccTot" );
    fAcceptanceHisto2D_rot = ( TH2F* )fAcceptanceFile->Get( hname );
    fAcceptanceHisto2D_rot_normalized = false;
    
    // read AZ dependent acceptance histograms from file and fit functions from file
    fAcceptancePhiHisto.clear();
    fAcceptancePhiFitFunction.clear();
    fAcceptancePhiHistoDeRot.clear();
    fAcceptancePhiFitFunctionDeRot.clear();
    for( int i = 0; i < 16; i++ )
    {
        sprintf( hname, "hAccPhi_%d", i );
        fAcceptancePhiHisto.push_back( ( TH1F* )fAcceptanceFile->Get( hname ) );
        sprintf( hname, "fAccPhi_%d", i );
        fAcceptancePhiFitFunction.push_back( ( TF1* )fAcceptanceFile->Get( hname ) );
        sprintf( hname, "hAccPhiDerot_%d", i );
        fAcceptancePhiHistoDeRot.push_back( ( TH1F* )fAcceptanceFile->Get( hname ) );
        sprintf( hname, "fAccPhiDerot_%d", i );
        fAcceptancePhiFitFunctionDeRot.push_back( ( TF1* )fAcceptanceFile->Get( hname ) );
    }
    // read AZ distributions
    hPhiDist = ( TH1F* )fAcceptanceFile->Get( "hPhiDist" );
    if( hPhiDist )
    {
        hPhiDist->Rebin( 5 );
    }
    hPhiDistDeRot = ( TH1F* )fAcceptanceFile->Get( "hPhiDistDeRot" );
    if( hPhiDistDeRot )
    {
        hPhiDistDeRot->Rebin( 5 );
    }
    
    return true;
}

/*

    plot acceptance curves

*/
TCanvas* VPlotRadialAcceptance::plotRadialAcceptance( TCanvas* cX , int iColor )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plot() error: data missing";
        return 0;
    }
    
    bool bPlotSame = false;
    // canvas
    if( cX )
    {
        cX->cd();
        bPlotSame = true;
    }
    else
    {
        char hname[2000];
        char htitle[2000];
        sprintf( hname, "cAcceptance_%s", VUtilities::removeSpaces( fName ).c_str() );
        sprintf( htitle, "%s", fName.c_str() );
        cX = new TCanvas( hname, htitle, 10, 10, 600, 600 );
        cX->SetGridx( 0 );
        cX->SetGridy( 0 );
    }
    
    // plot histograms
    if( fAcceptanceHisto )
    {
        fAcceptanceHisto->SetMinimum( fAxis_y_min );
        fAcceptanceHisto->SetMaximum( fAxis_y_max );
        fAcceptanceHisto->SetAxisRange( fAxis_x_min, fAxis_x_max );
        fAcceptanceHisto->SetTitle( "" );
        setHistogramPlottingStyle( fAcceptanceHisto, iColor, 1., 1.0, 20 );
        if( bPlotSame )
        {
            fAcceptanceHisto->Draw( "e same" );
        }
        else
        {
            fAcceptanceHisto->Draw( "e" );
        }
        fAcceptanceHisto->GetYaxis()->SetTitleOffset( 1.2 );
    }
    if( fAcceptanceFunction )
    {
        setFunctionPlottingStyle( fAcceptanceFunction, iColor );
        fAcceptanceFunction->Draw( "same" );
    }
    
    return cX;
}

/*
    plot radial acceptances per run

*/
bool VPlotRadialAcceptance::plotRadialAcceptance_perRun( bool iPlotTotalAcceptanceFunction )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plotRadialAcceptance_perRun() error:";
        cout << " cannot read radial acceptance file" << endl;
        return false;
    }
    vector< TH1F* > h1D;
    vector< TF1* >  h1F;
    
    TIter next( fAcceptanceFile->GetListOfKeys() );
    while( TObject* obj = next() )
    {
        string iName = obj->GetName();
        if( iName.find( "hRadialAcceptance_perRun_" ) != string::npos
                && iName.find( "Fit" ) == string::npos )
        {
            h1D.push_back( ( TH1F* )fAcceptanceFile->Get( iName.c_str() ) );
        }
        else if( iName.find( "fRadialAcceptance_perRun_" ) != string::npos )
        {
            h1F.push_back( ( TF1* )fAcceptanceFile->Get( iName.c_str() ) );
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // normalize and plot everything
    float nentries = ( float )h1D.size();
    
    int nx = ( int )( sqrt( nentries ) );
    int ny = nentries / ( int )( sqrt( nentries ) );
    if( nx * ny < nentries )
    {
        ny += 1;
    }
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "cRad_perRun" );
    sprintf( htitle, "acceptances per run" );
    TCanvas* cRad_perRun = new TCanvas( hname, htitle, 10, 10, 800, 800 );
    cRad_perRun->Divide( nx, ny );
    cRad_perRun->Draw();
    
    for( unsigned int i = 0; i < h1D.size(); i++ )
    {
        if( !h1D[i] )
        {
            continue;
        }
        TPad* i_C = ( TPad* )cRad_perRun->cd( i + 1 );
        i_C->SetGridx( 0 );
        i_C->SetGridy( 0 );
        
        h1D[i]->SetMinimum( fAxis_y_min );
        h1D[i]->SetMaximum( fAxis_y_max );
        h1D[i]->SetAxisRange( fAxis_x_min, fAxis_x_max );
        h1D[i]->Draw();
        
        if( h1F[i] )
        {
            h1F[i]->Draw( "same" );
        }
        if( iPlotTotalAcceptanceFunction && fAcceptanceFunction )
        {
            fAcceptanceFunction->SetLineColor( 1 );
            fAcceptanceFunction->SetLineStyle( 2 );
            fAcceptanceFunction->Draw( "same" );
        }
    }
    return true;
}

/*

    plot residuals between fit function and measured histogram

*/
TCanvas* VPlotRadialAcceptance::plotResiduals( TCanvas* cX, double i_res_min, double i_res_max, bool iPlotChi2 )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plotResiduals() error: data missing";
        return 0;
    }
    
    // canvas
    char hname[2000];
    if( cX )
    {
        cX->cd();
    }
    else
    {
        char htitle[2000];
        sprintf( hname, "cAcceptanceResiduals_%s", VUtilities::removeSpaces( fName ).c_str() );
        sprintf( htitle, "%s (residuals)", fName.c_str() );
        cX = new TCanvas( hname, htitle, 420, 10, 600, 600 );
        cX->SetGridx( 0 );
        cX->SetGridy( 0 );
    }
    
    if( fAcceptanceHisto && fAcceptanceFunction )
    {
        sprintf( hname, "%s_residual", fAcceptanceHisto->GetName() );
        TH1D* hRes = VHistogramUtilities::get_ResidualHistogram_from_TF1( hname, fAcceptanceHisto, fAcceptanceFunction );
        
        if( hRes )
        {
            setHistogramPlottingStyle( hRes );
            hRes->SetTitle( "" );
            hRes->SetMinimum( i_res_min );
            hRes->SetMaximum( i_res_max );
            hRes->SetAxisRange( fAxis_x_min, fAxis_x_max );
            hRes->SetXTitle( "distance to camera centre [deg]" );
            hRes->SetYTitle( "residial fit - measurement" );
            hRes->Draw();
            TLine* iL = new TLine( hRes->GetXaxis()->GetXmin(), 0., fAxis_x_max, 0. );
            iL->Draw();
            
            if( iPlotChi2 )
            {
                double sum2 = 0.;
                int n = 0;
                for( int i = 1; i <= fAcceptanceHisto->GetNbinsX(); i++ )
                {
                    if( fAcceptanceHisto->GetBinContent( i ) > 0. && fAcceptanceHisto->GetBinError( i ) > 0. )
                    {
                        sum2 += ( fAcceptanceHisto->GetBinContent( i ) - fAcceptanceFunction->Eval( fAcceptanceHisto->GetBinCenter( i ) ) )
                                * ( fAcceptanceHisto->GetBinContent( i ) - fAcceptanceFunction->Eval( fAcceptanceHisto->GetBinCenter( i ) ) )
                                / fAcceptanceHisto->GetBinError( i ) / fAcceptanceHisto->GetBinError( i );
                        n++;
                    }
                }
                sprintf( hname, "Fit Chi2/N: %.2f/%d", sum2, n );
                TText* iT = new TText( 0.5 * fAxis_x_max, 0.7 * i_res_max, hname );
                iT->Draw();
            }
        }
    }
    
    return cX;
}

/*
    service function for the 2D radial acceptance plotting

*/
void VPlotRadialAcceptance::plotRadialAcceptance2D_plotter( TH2F* h, bool iNorm, double iNormalizationRadius,
        bool iSmooth, float iMaxRadius,
        TPad* c_h2D, TPad* c_h2Dresidual )
{
    if( !h )
    {
        return;
    }
    char hname[2000];
    
    //////////////////////////////
    // normalize to centre bins
    if( !iNorm )
    {
        float i_nBins = 0;
        float i_nCounts = 0;
        for( int i = 1; i <= h->GetNbinsX(); i++ )
        {
            float x = h->GetXaxis()->GetBinCenter( i );
            for( int j = 1; j <= h->GetNbinsY(); j++ )
            {
                float y = h->GetYaxis()->GetBinCenter( j );
                
                if( x * x + y * y < iNormalizationRadius * iNormalizationRadius )
                {
                    i_nBins++;
                    i_nCounts += h->GetBinContent( i, j );
                }
                
                if( x * x + y * y > iMaxRadius * iMaxRadius )
                {
                    h->SetBinContent( i, j, 0. );
                }
            }
        }
        if( i_nCounts > 0. )
        {
            h->Scale( i_nBins / i_nCounts );
            
            cout << "Normalization: " << i_nBins / i_nCounts << endl;
        }
        h->SetMaximum( 1.2 );
        
        cout << "Normalization bins: " << i_nBins << endl;
        cout << "Normalization counts: " << i_nCounts << endl;
    }
    /////////////////////////////////////////////////////////////////////
    // smooth histogram
    if( iSmooth )
    {
        h->Smooth( 1 );
    }
    
    ///////////////////////////////////////////////////////////////
    // calculate residuals to radial acceptance
    sprintf( hname, "h2DResidual%s", h->GetName() );
    TH2F* h2DResidual = ( TH2F* )h->Clone( hname );
    h2DResidual->Reset();
    float i_plotMin = 50.;
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        float x = h->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            float y = h->GetYaxis()->GetBinCenter( j );
            
            if( fAcceptanceFunction && fAcceptanceFunction->Eval( sqrt( x * x + y * y ) ) )
            {
                float iRatio =  h->GetBinContent( i, j ) / fAcceptanceFunction->Eval( sqrt( x * x + y * y ) );
                h2DResidual->SetBinContent( i, j, iRatio );
                if( iRatio < i_plotMin )
                {
                    i_plotMin = iRatio;
                }
            }
            else
            {
                h2DResidual->SetBinContent( i, j, 0. );
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////
    // plot everything
    if( c_h2D && c_h2D->cd() )
    {
        c_h2D->SetGridx( 0 );
        c_h2D->SetGridy( 0 );
        h->SetTitle( "" );
        h->SetStats( 0 );
        h->Draw( "colz" );
    }
    if( c_h2Dresidual && c_h2Dresidual->cd() )
    {
        c_h2Dresidual->SetGridx( 0 );
        c_h2Dresidual->SetGridy( 0 );
        h2DResidual->SetTitle( "" );
        h2DResidual->SetStats( 0 );
        h2DResidual->SetMaximum( 1.3 );
        h2DResidual->SetMinimum( 0.7 );
        h2DResidual->Draw( "colz" );
    }
}

/*
 *
 * plot 2D acceptances (combined)
 *
 */
void VPlotRadialAcceptance::plotRadialAcceptance2D( bool iRotated, double iNormalizationRadius, bool iSmooth, float iMaxRadius )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plot() error: data missing";
        return;
    }
    TH2F* h = 0;
    bool iNorm = false;
    if( iRotated )
    {
        h = fAcceptanceHisto2D_rot;
        iNorm = fAcceptanceHisto2D_rot_normalized;
    }
    else
    {
        h = fAcceptanceHisto2D_derot;
        iNorm = fAcceptanceHisto2D_derot_normalized;
    }
    
    char hname[2000];
    char htitle[2000];
    sprintf( hname, "cAcceptance2D_%s_%d", VUtilities::removeSpaces( fName ).c_str(), iRotated );
    if( iRotated )
    {
        sprintf( htitle, "%s (2D, rotated)", fName.c_str() );
    }
    else
    {
        sprintf( htitle, "%s (2D, derotated)", fName.c_str() );
    }
    TCanvas* cXY = new TCanvas( hname, htitle, 10, 10, 600, 600 );
    cXY->SetGridx( 0 );
    cXY->SetGridy( 0 );
    
    
    sprintf( hname, "cAcceptanceResidual2D_%s_%d", VUtilities::removeSpaces( fName ).c_str(), iRotated );
    if( iRotated )
    {
        sprintf( htitle, "%s ratio to radial acceptance fit (2D, rotated)", fName.c_str() );
    }
    else
    {
        sprintf( htitle, "%s ratio to radial acceptance fit (2D, derotated)", fName.c_str() );
    }
    TCanvas* cXYR = new TCanvas( hname, htitle, 610, 10, 600, 600 );
    cXYR->SetGridx( 0 );
    cXYR->SetGridy( 0 );
    
    plotRadialAcceptance2D_plotter( h, iNorm, iNormalizationRadius, iSmooth, iMaxRadius,
                                    ( TPad* )cXY, ( TPad* )cXYR );
                                    
}

void VPlotRadialAcceptance::setAxisRange( double x_min, double x_max, double y_min, double y_max )
{
    fAxis_x_min = x_min;
    fAxis_x_max = x_max;
    fAxis_y_min = y_min;
    fAxis_y_max = y_max;
}


TCanvas* VPlotRadialAcceptance::plotPhiDependentRadialAcceptances( TCanvas* cX, int iIterator, bool iDeRot )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plotPhiDependentRadialAcceptances() error: data missing";
        return 0;
    }
    
    bool bPlotSame = false;
    // canvas
    if( cX )
    {
        cX->cd();
        bPlotSame = true;
    }
    else
    {
        char hname[2000];
        sprintf( hname, "cAcceptancePhi_%s_%d", VUtilities::removeSpaces( fName ).c_str(), ( int )iDeRot );
        string iTitle = fName;
        if( iDeRot )
        {
            fName += "(Phi dependent, derotated)";
        }
        else
        {
            fName += "(Phi dependent)";
        }
        cX = new TCanvas( hname, fName.c_str(), 60, 610, 600, 600 );
        cX->SetGridx( 0 );
        cX->SetGridy( 0 );
    }
    
    vector< TH1F* > iHisto;
    vector< TF1* > iF1;
    if( !iDeRot )
    {
        iHisto = fAcceptancePhiHisto;
        iF1    = fAcceptancePhiFitFunction;
    }
    else
    {
        iHisto = fAcceptancePhiHistoDeRot;
        iF1    = fAcceptancePhiFitFunctionDeRot;
    }
    
    // plot all histograms and plot them
    int i_color = 1;
    for( unsigned int i = 0; i < iHisto.size(); i += iIterator )
    {
        if( i < iHisto.size() && iHisto[i] )
        {
            cout << "Phi: " << iHisto[i]->GetTitle() << "\t color " << i_color << endl;
            iHisto[i]->SetMinimum( fAxis_y_min );
            iHisto[i]->SetMaximum( fAxis_y_max );
            iHisto[i]->SetAxisRange( fAxis_x_min, fAxis_x_max );
            iHisto[i]->SetTitle( "" );
            setHistogramPlottingStyle( iHisto[i], i_color, 1., 1.5, 20 );
            if( bPlotSame )
            {
                iHisto[i]->Draw( "e same" );
            }
            else
            {
                iHisto[i]->Draw( "e" );
                bPlotSame = true;
            }
            iHisto[i]->GetYaxis()->SetTitleOffset( 1.2 );
            if( i < iF1.size() && iF1[i] )
            {
                setFunctionPlottingStyle( iF1[i], i_color, 1. );
                iF1[i]->SetLineStyle( 2 );
                iF1[i]->Draw( "same" );
            }
            
            i_color++;
        }
    }
    if( fAcceptanceFunction && bPlotSame )
    {
        setFunctionPlottingStyle( fAcceptanceFunction, 1, 2 );
        fAcceptanceFunction->Draw( "same" );
    }
    
    return cX;
}

TCanvas*  VPlotRadialAcceptance::plotPhiDistributions( TCanvas* cX, int iColor )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plotPhiDistributions() error: data missing";
        return 0;
    }
    
    bool bPlotSame = false;
    // canvas
    if( cX )
    {
        cX->cd();
        bPlotSame = true;
    }
    else
    {
        char hname[2000];
        sprintf( hname, "cPhiDistribution%s", VUtilities::removeSpaces( fName ).c_str() );
        string iTitle = fName;
        fName += ", Phi distribution";
        cX = new TCanvas( hname, fName.c_str(), 60, 610, 600, 600 );
        cX->SetGridx( 0 );
        cX->SetGridy( 0 );
    }
    
    if( hPhiDist )
    {
        setHistogramPlottingStyle( hPhiDist, iColor );
        hPhiDist->SetLineWidth( 2 );
        hPhiDist = ( TH1F* )VHistogramUtilities::normalizeTH1( ( TH1* )hPhiDist, false );
        hPhiDist->GetXaxis()->SetTitle( "azimuth (camera) [deg]" );
        hPhiDist->SetMinimum( 0 );
        if( bPlotSame )
        {
            hPhiDist->Draw( "he same" );
        }
        else
        {
            hPhiDist->Draw( "he" );
        }
        
        TLine* iL = new TLine( hPhiDist->GetXaxis()->GetXmin(), 1., hPhiDist->GetXaxis()->GetXmax(), 1. );
        iL->SetLineStyle( 2 );
        iL->Draw();
    }
    return cX;
}

/*
    plot 2D acceptances and relative acceptances for each run

*/
bool VPlotRadialAcceptance::plotRadialAcceptance2D_perRun( double iNormalizationRadius, bool iSmooth, float iMaxRadius )
{
    if( !fAcceptanceFile || fAcceptanceFile->IsZombie() )
    {
        cout << "VPlotRadialAcceptance::plotRadialAcceptance2D_perRun() error:";
        cout << " cannot read radial acceptance file" << endl;
        return false;
    }
    vector< TH2F* > h2D;
    
    TIter next( fAcceptanceFile->GetListOfKeys() );
    while( TObject* obj = next() )
    {
        string iName = obj->GetName();
        if( iName.find( "hXYAccRun_" ) != string::npos )
        {
            h2D.push_back( ( TH2F* )fAcceptanceFile->Get( iName.c_str() ) );
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // normalize and plot everything
    float nentries = ( float )h2D.size();
    
    int nx = ( int )( sqrt( nentries ) );
    int ny = nentries / ( int )( sqrt( nentries ) );
    if( nx * ny < nentries )
    {
        ny += 1;
    }
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "cAllMaps" );
    sprintf( htitle, "all maps" );
    TCanvas* allMaps = new TCanvas( hname, htitle, 10, 10, 800, 800 );
    allMaps->Divide( nx, ny );
    allMaps->Draw();
    
    sprintf( hname, "cAllMapsR" );
    sprintf( htitle, "all maps (residuals)" );
    TCanvas* allMapsR = new TCanvas( hname, htitle, 110, 10, 800, 800 );
    allMapsR->Divide( nx, ny );
    
    for( unsigned int i = 0; i < h2D.size(); i++ )
    {
        if( h2D[i] )
        {
            plotRadialAcceptance2D_plotter( h2D[i],
                                            false, iNormalizationRadius,
                                            iSmooth, iMaxRadius,
                                            ( TPad* )allMaps->cd( i + 1 ),
                                            ( TPad* )allMapsR->cd( i + 1 ) );
        }
    }
    
    
    
    
    return true;
}
