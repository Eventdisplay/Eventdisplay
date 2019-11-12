/*  \file VPlotVERITASPerformance
 *  make plots for VERITAS public performance page
 *
 *  note: stil many hardwired options at this point
*/


#include "VPlotVERITASPerformance.h"

VPlotVERITASPerformanceEpochData::VPlotVERITASPerformanceEpochData()
{
    fName = "";
    fTitle = "";
    fPlotThisEpoch = false;
    fColor = 1;
    fLineStyle = 1;
    fCutsIRFFile = "";
    fZe = 20.;
    fNoiseLevel = 0;
    fGammaRayRate = 0.;
    fBckRate = 0.;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

VPlotVERITASPerformance::VPlotVERITASPerformance()
{
    setPlotEnergyRange();
    
    fPlotter = new VPlotInstrumentResponseFunction();
    resetSettings();
    
}

bool VPlotVERITASPerformance::setDataSet( bool iV6Only, string iCut )
{
    if( !fPlotter )
    {
        return false;
    }
    // clear data vectors and plotter
    fEpochData.clear();
    fPlotter->resetInstrumentResponseData();
    
    //////////////////////////////////////////
    // set data for different cuts and epochs
    
    // moderate cuts
    if( iCut.find( "moderate" ) != string::npos )
    {
        if( !iV6Only )
        {
            fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
            fEpochData.back()->fName = "V4-moderate";
            fEpochData.back()->fTitle = "Configuration 2007-2009 (V4)";
            fEpochData.back()->fPlotThisEpoch = true;
            fEpochData.back()->fColor = 8;
            fEpochData.back()->fCutsIRFFile = "IRF-GRISU-SW6-moderate2Tel-V4.root";
            fEpochData.back()->fNoiseLevel = 150;
            
            fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
            fEpochData.back()->fName = "V5-moderate";
            fEpochData.back()->fTitle = "Configuration 2009-2012 (V5)";
            fEpochData.back()->fPlotThisEpoch = true;
            fEpochData.back()->fColor = 46;
            fEpochData.back()->fCutsIRFFile = "IRF-GRISU-SW6-moderate2Tel-V5.root";
            fEpochData.back()->fNoiseLevel = 150;
            fEpochData.back()->fGammaRayRate = 6.4;
            fEpochData.back()->fBckRate = 0.21;
        }
        
        fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
        fEpochData.back()->fName = "V6-moderate";
        fEpochData.back()->fTitle = "Configuration 2012-today (V6)";
        fEpochData.back()->fPlotThisEpoch = true;
        fEpochData.back()->fColor = 1;
        fEpochData.back()->fCutsIRFFile = "IRF-CARE-SW18-moderate2Tel-V6.root";
        fEpochData.back()->fNoiseLevel = 170;
        fEpochData.back()->fGammaRayRate = 6.8;
        fEpochData.back()->fBckRate = 0.22;
    }
    else if( iCut.find( "V6V5-operationModes" ) != string::npos )
    {
    
        fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
        fEpochData.back()->fName = "V4-soft";
        fEpochData.back()->fTitle = "Configuration 2007-2009 (V4): soft cuts";
        fEpochData.back()->fPlotThisEpoch = true;
        fEpochData.back()->fColor = 8;
        fEpochData.back()->fCutsIRFFile = "IRF-GRISU-SW6-Soft2Tel-V4.root";
        fEpochData.back()->fNoiseLevel = 150;
        
        fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
        fEpochData.back()->fName = "V5-soft";
        fEpochData.back()->fTitle = "Configuration 2009-2012 (V5): soft cuts";
        fEpochData.back()->fPlotThisEpoch = true;
        fEpochData.back()->fColor = 46;
        fEpochData.back()->fCutsIRFFile = "IRF-GRISU-SW6-Soft2Tel-V5.root";
        fEpochData.back()->fNoiseLevel = 150;
        
        fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
        fEpochData.back()->fName = "V6-soft";
        fEpochData.back()->fTitle = "Configuration 2012-today (V6): soft cuts";
        fEpochData.back()->fPlotThisEpoch = true;
        fEpochData.back()->fColor = 1;
        fEpochData.back()->fCutsIRFFile = "IRF-CARE-SW18-Soft2Tel-V6.root";
        fEpochData.back()->fNoiseLevel = 170;
        
    }
    if( iCut.find( "RV" ) != string::npos )
    {
        fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
        fEpochData.back()->fName = "V6-redHV";
        fEpochData.back()->fTitle = "Reduced HV operations (V6): soft cuts";
        fEpochData.back()->fPlotThisEpoch = true;
        fEpochData.back()->fColor = 9;
        fEpochData.back()->fLineStyle = 2;
        fEpochData.back()->fCutsIRFFile = "effArea-v470-auxv01-CARE-Cut-NTel2-PointSource-Soft-GEO-V6-ATM21-redHV-T1234.root";
        fEpochData.back()->fNoiseLevel = 300;
    }
    if( iCut.find( "UV" ) != string::npos )
    {
        fEpochData.push_back( new VPlotVERITASPerformanceEpochData() );
        fEpochData.back()->fName = "V6-UVFilter";
        fEpochData.back()->fTitle = "UV filter operations (V6): soft cuts";
        fEpochData.back()->fPlotThisEpoch = true;
        fEpochData.back()->fColor = 41;
        fEpochData.back()->fLineStyle = 2;
        fEpochData.back()->fCutsIRFFile = "effArea-v470-auxv01-CARE-Cut-NTel2-PointSource-Soft-GEO-V6-ATM21-UV-T1234.root";
        fEpochData.back()->fNoiseLevel = 230;
        fEpochData.back()->fZe = 19.;
    }
    
    
    for( unsigned int i = 0; i < fEpochData.size(); i++ )
    {
        fEpochData[i]->fPlotThisEpoch = true;
        if( !fPlotter->addInstrumentResponseData( fEpochData[i]->fCutsIRFFile, fEpochData[i]->fZe,
                0.5, 0, 1.5,
                fEpochData[i]->fNoiseLevel, "A_MC",
                fEpochData[i]->fColor, fEpochData[i]->fLineStyle,
                20, 1, 0.2, 1.e6, 4., "l", 0.20 ) )
        {
            return false;
        }
    }
    
    return true;
}

void VPlotVERITASPerformance::resetSettings()
{
    fPlotter->setPlottingAxis( "energy_Lin", "X", true, fPlotEnergy_TeV_min, fPlotEnergy_TeV_max, "energy [TeV]" );
    fPlotter->setPlottingLogEnergyAxis( true );
}

void VPlotVERITASPerformance::setPlotAllEpochs()
{
    for( unsigned int i = 0; i < fEpochData.size(); i++ )
    {
        fEpochData[i]->fPlotThisEpoch = true;
    }
}

/*

    plot energy resolution for V4, V5 and V6

*/
void VPlotVERITASPerformance::plotEnergyResolutionForDifferentEpochs( double iMax )
{
    if( !fPlotter )
    {
        return;
    }
    
    TCanvas* c = fPlotter->plotEnergyResolution( 0., iMax, 0, true );
    
    TLegend* iL = getLegend();
    if( iL )
    {
        iL->SetHeader( "VERITAS Energy Resolution" );
        iL->Draw();
    }
    plotVersionNumber( c );
    printCanvas( c, "EnergyResolutionEpochs" );
}

/*

    plot angular resolution for V4, V5 and V6

*/
void VPlotVERITASPerformance::plotAngularResolutionForDifferentEpochs( double iMaxTheta )
{
    if( !fPlotter )
    {
        return;
    }
    
    TCanvas* c = fPlotter->plotAngularResolution( "energy", "68", 0.005, iMaxTheta );
    c->SetName( "AngularResolutionVsEnergy" );
    TGraphErrors* iG = fPlotter->getLastPlottedGraph();
    if( iG )
    {
        iG->RemovePoint( 0 );
    }
    TLegend* iL = getLegend();
    if( iL )
    {
        iL->SetHeader( "VERITAS Angular Resolution" );
        iL->Draw();
    }
    plotVersionNumber( c );
    printCanvas( c, "AngularResolutionEpochs" );
}

/*

    plot and return a legend

*/
TLegend* VPlotVERITASPerformance::getLegend( bool iLeft, bool iBottom )
{
    float iY = 0.7;
    float iY_d = 0.15;
    float iX = 0.25;
    float iX_d = 0.62;
    if( iLeft )
    {
        iX = 0.19;
        iX_d = 0.41;
    }
    if( iBottom )
    {
        iY = 0.15;
        iY_d = 0.2;
        iX = 0.30;
    }
    TLegend* iL = new TLegend( iX, iY, iX + iX_d, iY + iY_d );
    iL->SetFillColor( 10 );
    for( unsigned int i = 0; i < fEpochData.size(); i++ )
    {
        TGraph* iG = new TGraph( 1 );
        iG->SetLineWidth( 3 );
        iG->SetLineStyle( fEpochData[i]->fLineStyle );
        iG->SetLineColor( fEpochData[i]->fColor );
        iL->AddEntry( iG, fEpochData[i]->fTitle.c_str(), "L" );
    }
    return iL;
}

/*

    plot angular resolution for V6 vs HAWC and Fermi LAT performance

*/
void VPlotVERITASPerformance::plotAngularResolutionInComparison( bool iLatPass8 )
{

    if( !fPlotter )
    {
        return;
    }
    
    setPlotEnergyRange( 0.01, 50. );
    resetSettings();
    
    TCanvas* c = fPlotter->plotAngularResolution( "energy", "68", 0.005, 0.35 );
    
    ///////////////////////////////////////////////////////////////
    // read text file for data from other instruments
    
    double x = 0.;
    double y = 0.;
    int z = 0;
    // HAWC
    TGraph* hawc = new TGraph( "Instrument/HAWC_AngRes68.dat" );
    if( hawc )
    {
        for( int i = 0; i < hawc->GetN(); i++ )
        {
            hawc->GetPoint( i, x, y );
            if( x > 0. )
            {
                hawc->SetPoint( z, log10( x ), y );
                z++;
            }
        }
        hawc->SetLineWidth( 3 );
        hawc->SetLineColor( 30 );
        hawc->SetLineStyle( 4 );
        hawc->Draw( "c" );
    }
    // LAT Pass 8
    string iLatFile = "Instrument/LAT_pass8_AngRes68.dat";
    if( !iLatPass8 )
    {
        iLatFile = "Instrument/LAT_pass7v15_AngRes68.dat";
    }
    TGraph* latpx = new TGraph( iLatFile.c_str() );
    if( latpx )
    {
        z = 0;
        for( int i = 0; i < latpx->GetN(); i++ )
        {
            latpx->GetPoint( i, x, y );
            if( x > 0. )
            {
                latpx->SetPoint( z, log10( x ), y );
                z++;
            }
        }
        latpx->Print();
        latpx->SetLineWidth( 3 );
        latpx->SetLineStyle( 5 );
        latpx->SetLineColor( 46 );
        latpx->Draw( "c" );
    }
    // MAGIC
    TGraph* magic = new TGraph( "Instrument/MAGIC_AngRes68.dat" );
    if( magic )
    {
        z = 0;
        for( int i = 0; i < magic->GetN(); i++ )
        {
            magic->GetPoint( i, x, y );
            if( x > 0. )
            {
                magic->SetPoint( z, log10( x / 1.e3 ), y );
                z++;
            }
        }
        magic->SetLineStyle( 5 );
        magic->SetLineColor( fOtherInstrumentColor );
        //          magic->Draw( "c" );
    }
    
    TLegend* iL = getLegend( true );
    if( iL )
    {
        iL->SetHeader( "VERITAS Angular Resolution" );
        if( latpx && !iLatPass8 )
        {
            iL->AddEntry( latpx, "Fermi LAT P7_V15", "L" );
        }
        else if( latpx )
        {
            iL->AddEntry( latpx, "LAT Pass 8", "L" );
        }
        
        if( hawc )
        {
            iL->AddEntry( hawc, "HAWC (Abeysekara et al 2013)", "L" );
        }
        iL->Draw();
    }
    
    plotVersionNumber( c );
    printCanvas( c, "AngularResolution" );
}

/*

    plot effective area

*/
void VPlotVERITASPerformance::plotEffectiveArea( float iEnergy_TeV_Max, bool iPlotOtherInstruments )
{

    if( !fPlotter )
    {
        return;
    }
    
    setPlotEnergyRange( 0.03, iEnergy_TeV_Max );
    resetSettings();
    
    TCanvas* c = fPlotter->plotEffectiveArea( 10., 5.e5, 0, false );
    
    TLegend* iL = getLegend( false, true );
    
    // MAGIC effective area (ApJ 2014)
    if( iPlotOtherInstruments )
    {
        TGraph* magic = new TGraph( "Instrument/MAGIC_EffectiveArea.dat" );
        if( magic )
        {
            magic->SetLineStyle( 5 );
            magic->SetLineColor( 6 );
            magic->Draw( "c" );
            iL->AddEntry( magic, "MAGIC (ApJ 2014)", "l" );
        }
    }
    // draw header
    if( iL )
    {
        iL->SetHeader( "VERITAS Effective Area" );
        iL->Draw();
    }
    plotVersionNumber( c );
    printCanvas( c, "EffectiveArea" );
    
}

/*

    plot sensitivity vs time

    (calculate from gamma-ray and background rates)
*/
void VPlotVERITASPerformance::plotSensitivity_vs_time( int iPlotID, bool iGuidingLines )
{
    VSensitivityCalculator* i_Sens = new VSensitivityCalculator();
    i_Sens->setSignificanceParameter( 5., 0. );
    i_Sens->setSourceStrengthRange_CU( 1.e-3, 5. );
    
    vector< double > s;
    s.push_back( 0.01 );
    s.push_back( 0.02 );
    s.push_back( 0.05 );
    s.push_back( 0.10 );
    i_Sens->setSourceStrengthVector_CU( s );
    
    // QQQ  should be selected from data vector!
    i_Sens->addDataSet( 6.8, 0.22, 1. / 10., "V6-Moderate" );
    i_Sens->addDataSet( 6.4, 0.21, 1. / 10., "V5-Moderate" );
    TCanvas* c = i_Sens->plotObservationTimevsFlux( iPlotID, 0, 2, 2, iGuidingLines );
    
    i_Sens->list_sensitivity( iPlotID );
    
    TPaveText* iT = new TPaveText( 0.42, 0.75, 0.85, 0.85, "NDC" );
    iT->SetFillColor( 10 );
    iT->AddText( "VERITAS 2012-today (V6)" );
    iT->AddText( "E > 240 GeV" );
    iT->Draw();
    
    plotVersionNumber( c );
    printCanvas( c, "ObservingTimeVsFlux" );
    
}

void VPlotVERITASPerformance::printCanvas( TCanvas* c, string iTitle )
{
    if( c && fPlotPrintName.size() > 0 )
    {
        c->Update();
        
        char htitle[500];
        sprintf( htitle, "%s-%s.pdf", fPlotPrintName.c_str(), iTitle.c_str() );
        c->Print( htitle );
    }
}

void VPlotVERITASPerformance::plotVersionNumber( TCanvas* c )
{
    if( c )
    {
        c->cd();
        
        if( fVersionNumber.size() > 0 )
        {
            TText* iV = new TText();
            iV->SetNDC( true );
            iV->SetTextFont( 42 );
            iV->SetTextAngle( 90. );
            iV->SetTextSize( iV->GetTextSize() * 0.5 );
            iV->DrawTextNDC( 0.96, 0.2, fVersionNumber.c_str() );
        }
    }
}

/*

   plot all plots for the performance page

*/
void VPlotVERITASPerformance::plotPlotsForPerformancePage( string iPrintName, string iVersionNumber )
{
    setPrintName( iPrintName );
    setVersionNumber( iVersionNumber );
    
    ////////////////////////////////////
    // V6 sensitivity (moderate cuts)
    plotSensitivity_vs_time( 0, true );
    
    ////////////////////////////////////
    // V4/V5/V6
    setDataSet( false );
    
    // angular resolution
    plotAngularResolutionForDifferentEpochs( 0.25 );
    
    // energy resolution
    plotEnergyResolutionForDifferentEpochs();
    
    // angular resolutions in comparison with other instruments
    plotAngularResolutionInComparison( false );
    
    ////////////////////////////////////
    // V4/V5/V6 and red HV soft
    setDataSet( true, "V6V5-operationModes RV UV" );
    
    // effective area
    plotEffectiveArea();
    
}

TGraphAsymmErrors* VPlotVERITASPerformance::readDifferentialSensitivity( string iFile )
{
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        return 0;
    }
    TGraphAsymmErrors* iG = new TGraphAsymmErrors( 1 );
    double x = 0.;
    double y = 0.;
    double y_l = 0.;
    double y_u = 0.;
    int z = 0;
    while( is >> x >> y >> y_l >> y_u )
    {
        iG->SetPoint( z, x, y );
        iG->SetPointEYhigh( z, y_u );
        iG->SetPointEYlow( z, y_l );
        z++;
    }
    return iG;
}

void VPlotVERITASPerformance::plotDifferentialSensitivity( double iE_min, double iE_max )
{

    TCanvas* c = new TCanvas( "cDiffSens", "differential sensitivity (CU)", 10, 10, 800, 600 );
    c->SetGridx( 0 );
    c->SetGridy( 0 );
    c->SetLogx( 1 );
    c->SetLogy( 1 );
    c->Draw();
    
    TH1D* hnull = new TH1D( "hDiffSens", "", 100, iE_min, iE_max );
    hnull->SetStats( 0 );
    hnull->SetXTitle( "energy [TeV]" );
    hnull->SetYTitle( "Source strength [Crab units]" );
    hnull->SetMinimum( 5.e-3 );
    hnull->SetMaximum( 1.2e-0 );
    
    hnull->Draw();
    
    TGraphAsymmErrors* iGSoft = readDifferentialSensitivity( "Instrument/VERITAS-DiffSens-SoftCut.dat" );
    if( iGSoft )
    {
        iGSoft->SetFillColor( 4 );
        iGSoft->Draw( "3" );
    }
    TGraphAsymmErrors* iGMedium = readDifferentialSensitivity( "Instrument/VERITAS-DiffSens-MediumCut.dat" );
    if( iGMedium )
    {
        iGMedium->SetFillColor( 3 );
        iGMedium->Draw( "3" );
    }
    TGraphAsymmErrors* iGHard = readDifferentialSensitivity( "Instrument/VERITAS-DiffSens-HardCut.dat" );
    if( iGHard )
    {
        iGHard->SetFillColor( 2 );
        iGHard->Draw( "3" );
    }
}


