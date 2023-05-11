/*  \class VLightCurvePlotter
    \brief plot light curves

*/

#include "VLightCurvePlotter.h"

VLightCurvePlotter::VLightCurvePlotter()
{
    fDebug = false;
    
    reset();
}

VLightCurvePlotter::VLightCurvePlotter( vector< VFluxDataPoint > iDataVector )
{
    fDebug = false;
    
    reset();
    
    setDataVector( iDataVector );
}

void VLightCurvePlotter::reset()
{
    fName          = "";
    fMinEnergy_TeV = -99.;
    fMaxEnergy_TeV = 1.e10;
    fCanvasLightCurve = 0;
    setLightCurveTimeAxis();
    setLightCurveFluxAxis();
    fIgnoreUpperLimits = false;
    setPlotObservingIntervals();
    setLightCurveTimeAxis_to_OrbitalPhase();
}

void VLightCurvePlotter::setDataVector( vector< VFluxDataPoint > iDataVector )
{
    fFluxDataVector = iDataVector;
    
    if( fFluxDataVector.size() > 0 )
    {
        fName          = fFluxDataVector[0].fName;
        fMinEnergy_TeV = fFluxDataVector[0].fMinEnergy_TeV;
        fMaxEnergy_TeV = fFluxDataVector[0].fMaxEnergy_TeV;
        fOrbitalPhaseData = fFluxDataVector[0].fOrbitalPhaseData;
    }
    
    setLightCurveTimeAxis();
    setLightCurveFluxAxis();
}


TGraphAsymmErrors* VLightCurvePlotter::plotFluxes_vs_Variable( string iVariable, string iAxisTitle )
{
    char hname[200];
    char htitle[200];
    sprintf( htitle, "c%s", iVariable.c_str() );
    sprintf( hname, "fluxes vs %s", iAxisTitle.c_str() );
    
    TCanvas* cFCanvas = new TCanvas( htitle, hname, 400, 10, 600, 400 );
    cFCanvas->SetGridx( 0 );
    cFCanvas->SetGridy( 0 );
    cFCanvas->Draw();
    
    TGraphAsymmErrors* iFluxGraph = new TGraphAsymmErrors( ( int )fFluxDataVector.size() );
    iFluxGraph->SetTitle( "" );
    iFluxGraph->SetMarkerStyle( 24 );
    iFluxGraph->SetMarkerSize( 1. );
    iFluxGraph->SetLineWidth( 2 );
    
    int z = 0;
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        if( fFluxDataVector[i].fMJD > 10 )
        {
            if( iVariable == "fPedVars" )
            {
                iFluxGraph->SetPoint( z, fFluxDataVector[i].fPedvars, fFluxDataVector[i].fFlux_Rolke );
            }
            else if( iVariable == "fWobbleOffset" )
            {
                iFluxGraph->SetPoint( z, fFluxDataVector[i].fWobbleOffset, fFluxDataVector[i].fFlux_Rolke );
            }
            else if( iVariable == "fElevation" )
            {
                iFluxGraph->SetPoint( z, 90. - fFluxDataVector[i].fZe, fFluxDataVector[i].fFlux_Rolke );
            }
            iFluxGraph->SetPointEYhigh( z, fFluxDataVector[i].fFluxCI_up_1sigma );
            iFluxGraph->SetPointEYlow( z, fFluxDataVector[i].fFluxCI_lo_1sigma );
            z++;
        }
    }
    iFluxGraph->Draw( "ap" );
    iFluxGraph->GetHistogram()->SetXTitle( iAxisTitle.c_str() );
    // get flux axis string from light curve data
    if( fFluxDataVector.size() > 0 )
    {
        fPlotting_FluxAxis_title = fFluxDataVector[0].fFluxUnitString;
    }
    iFluxGraph->GetHistogram()->SetYTitle( getLightCurveFluxAxisTitle().c_str() );
    
    // draw a line at zero
    if( iFluxGraph->GetHistogram()->GetYaxis()->GetXmin() < 0. )
    {
        TLine* iL2 = new TLine( iFluxGraph->GetHistogram()->GetXaxis()->GetXmin(), 0., iFluxGraph->GetHistogram()->GetXaxis()->GetXmax(), 0. );
        iL2->SetLineStyle( 2 );
        iL2->Draw();
    }
    
    return iFluxGraph;
}

/*

   plot fluxes vs pedestal variations
   (check for systematics)

*/
TGraphAsymmErrors* VLightCurvePlotter::plotFluxes_vs_pedvars()
{
    return plotFluxes_vs_Variable( "fPedVars", "mean pedvars" );
}


/*

   plot fluxes vs wobble offsets
   (check for systematics)

*/
TGraphAsymmErrors* VLightCurvePlotter::plotFluxes_vs_wobbleOffset()
{
    return plotFluxes_vs_Variable( "fWobbleOffset", "wobble offset [deg]" );
}

/*

   plot fluxes vs elevation
   (check for systematics)

*/
TGraphAsymmErrors* VLightCurvePlotter::plotFluxes_vs_elevation()
{
    return plotFluxes_vs_Variable( "fElevation", "elevation [deg]" );
}

/*

    set min/max for time axis

    units of time axis migth be MJD or Phase

*/
void VLightCurvePlotter::setLightCurveTimeAxis( double iXmin, double iXmax, double iPlottingTimeAxis_tolerance_dayfraction, double iPlotting_MJD_offset_days )
{
    fPlotting_TimeAxis_min = iXmin;
    fPlotting_TimeAxis_max = iXmax;
    fPlottingMJD_tolerance_dayfraction = iPlottingTimeAxis_tolerance_dayfraction;
    fPlotting_MJD_offset_days = iPlotting_MJD_offset_days;
    
    // default values are start and end of light curves
    bool iAddTolerance = false;
    if( fPlotting_TimeAxis_min < 0 && fFluxDataVector.size() > 0 )
    {
        VFluxDataPoint i_temp = *min_element( fFluxDataVector.begin(), fFluxDataVector.end() );
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            fPlotting_TimeAxis_min = fOrbitalPhaseData.getOrbitalPhase( i_temp.fMJD );
        }
        else
        {
            fPlotting_TimeAxis_min = i_temp.fMJD;
        }
        iAddTolerance = true;
    }
    if( fPlotting_TimeAxis_max < 0 && fFluxDataVector.size() > 0 )
    {
        VFluxDataPoint i_temp = *max_element( fFluxDataVector.begin(), fFluxDataVector.end() );
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            fPlotting_TimeAxis_max = fOrbitalPhaseData.getOrbitalPhase( i_temp.fMJD );
        }
        else
        {
            fPlotting_TimeAxis_max = i_temp.fMJD;
        }
        iAddTolerance = true;
    }
    if( iAddTolerance )
    {
        double i_diff = fPlotting_TimeAxis_max - fPlotting_TimeAxis_min;
        fPlotting_TimeAxis_min -= fPlottingMJD_tolerance_dayfraction * i_diff;
        fPlotting_TimeAxis_max += fPlottingMJD_tolerance_dayfraction * i_diff;
    }
    
    // offset used for plotting
    fPlotting_TimeAxis_min -= fPlotting_MJD_offset_days;
    fPlotting_TimeAxis_max -= fPlotting_MJD_offset_days;
    
    // make sure that the values are not identical
    if( TMath::Abs( fPlotting_TimeAxis_max - fPlotting_TimeAxis_min ) < 1.e-4 )
    {
        fPlotting_TimeAxis_min -= 0.5;
        fPlotting_TimeAxis_max += 0.5;
    }
    
}

/*

    set min/max for flux axis

    QQQQ cannot deal yet with flux values < 0

*/
void VLightCurvePlotter::setLightCurveFluxAxis( double iYmin, double iYmax, string iAxisTitle )
{
    fPlotting_Flux_min = iYmin;
    fPlotting_Flux_max = iYmax;
    fPlotting_FluxAxis_title = iAxisTitle;
    
    // default values are read from data vector
    if( fFluxDataVector.size() > 0 )
    {
        double i_max = -1.e10;
        double i_min =  1.e10;
        double i_f_min = 1.e10;
        double i_f_max = -1.e10;
        for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
        {
            if( fFluxDataVector[i].isSignificantDataPoint() || fIgnoreUpperLimits )
            {
                i_f_min = fFluxDataVector[i].fFlux_Rolke - fFluxDataVector[i].fFluxCI_lo_1sigma;
                i_f_max = fFluxDataVector[i].fFlux_Rolke + fFluxDataVector[i].fFluxCI_up_1sigma;
            }
            else if( fFluxDataVector[i].fFluxUL > 0. )
            {
                i_f_min = fFluxDataVector[i].fFluxUL;
                i_f_max = fFluxDataVector[i].fFluxUL;
            }
            else
            {
                continue;
            }
            if( i_f_max > i_max )
            {
                i_max = i_f_max;
            }
            if( i_f_min < i_min )
            {
                i_min = i_f_min;
            }
        }
        // check if values are preset
        if( fPlotting_Flux_min < -8.e10 )
        {
            fPlotting_Flux_min = 0.8 * i_min;
        }
        if( fPlotting_Flux_max < -8.e10 )
        {
            fPlotting_Flux_max = 1.2 * i_max;
        }
        // make sure that the values are not identical
        if( TMath::Abs( fPlotting_Flux_max - fPlotting_Flux_min ) < 1.e-55 )
        {
            fPlotting_Flux_min *= 0.9;
            fPlotting_Flux_max *= 1.1;
        }
    }
}



/*

      plot a light curve for the given data

      TCanvas* iCanvasLightCurve      	0 for a new canvas, or pointer to an existing canvas
      string iCanvasName              	canvas name (for new canvas)

*/
TCanvas* VLightCurvePlotter::plotLightCurve( TCanvas* iCanvasLightCurve, string iCanvasName, string iPlottingOption, double iMaxMJDError )
{
    char hname[800];
    char htitle[5000];
    
    TH1D* hLightCurve = 0;
    
    ////////////////////////////
    // set canvas
    if( !iCanvasLightCurve )
    {
        sprintf( hname, "%s", iCanvasName.c_str() );
        if( fName.size() > 0 )
        {
            sprintf( htitle, "%s", fName.c_str() );
        }
        else
        {
            sprintf( htitle, "light curve" );
        }
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            sprintf( hname, "%s_%d_%d", iCanvasName.c_str(), ( int )fOrbitalPhaseData.fZeroPhase_MJD, ( int )fOrbitalPhaseData.fOrbit_days );
            sprintf( htitle, "%s, T_{0}=%.1f, period=%.1f days", fName.c_str(), fOrbitalPhaseData.fZeroPhase_MJD, fOrbitalPhaseData.fOrbit_days );
        }
        
        fCanvasLightCurve = new TCanvas( hname, htitle, 10, 10, getPlottingCanvasX(), getPlottingCanvasY() );
        fCanvasLightCurve->SetGridx( 0 );
        fCanvasLightCurve->SetGridy( 0 );
        
        ///////////////////////////////////
        // initialize null histograms
        double i_xmin = fPlotting_TimeAxis_min;
        double i_xmax = fPlotting_TimeAxis_max;
        sprintf( hname, "MJD" );
        if( fPlotting_MJD_offset_days > 0. )
        {
            sprintf( hname, "MJD - %.1f days", fPlotting_MJD_offset_days );
        }
        sprintf( htitle, "hLightCurve" );
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            i_xmin = 0.;
            i_xmax = 1.;
            sprintf( hname, "phase" );
            sprintf( htitle, "hLightCurve_%d_%d", ( int )fOrbitalPhaseData.fZeroPhase_MJD, ( int )fOrbitalPhaseData.fOrbit_days );
        }
        
        hLightCurve = new TH1D( htitle, "", 100, i_xmin, i_xmax );
        hLightCurve->SetStats( 0 );
        hLightCurve->SetXTitle( hname );
        hLightCurve->GetXaxis()->CenterTitle( true );
        // get flux axis string from light curve data
        if( fFluxDataVector.size() > 0 )
        {
            fPlotting_FluxAxis_title = fFluxDataVector[0].fFluxUnitString;
        }
        hLightCurve->SetYTitle( getLightCurveFluxAxisTitle().c_str() );
        hLightCurve->SetMinimum( getLightCurveFluxAxisRange_Min() );
        hLightCurve->SetMaximum( getLightCurveFluxAxisRange_Max() );
        hLightCurve->Draw( "" );
        hLightCurve->Draw( "AH" );
        
        cout << "Plotting range: ";
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            cout << i_xmin << " < orbital phase < " << i_xmax;
        }
        else
        {
            cout << i_xmin << " < time < " << i_xmax;
        }
        cout << ", " << hLightCurve->GetMinimum() << " < flux < " << hLightCurve->GetMaximum() << endl;
        
        plot_nullHistogram( fCanvasLightCurve, hLightCurve, false, true, 1.2, fPlotting_TimeAxis_min, fPlotting_TimeAxis_max );
        
    }
    // canvas exists - get histogram in canvas
    else
    {
        fCanvasLightCurve = iCanvasLightCurve;
        fCanvasLightCurve->cd();
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            sprintf( htitle, "hLightCurve_%d_%d", ( int )fOrbitalPhaseData.fZeroPhase_MJD, ( int )fOrbitalPhaseData.fOrbit_days );
        }
        else
        {
            sprintf( htitle, "hLightCurve" );
        }
        hLightCurve = ( TH1D* )fCanvasLightCurve->GetListOfPrimitives()->FindObject( htitle );
        if( !hLightCurve )
        {
            hLightCurve = ( TH1D* )fCanvasLightCurve->GetListOfPrimitives()->FindObject( "hLightCurve" );
        }
        if( !hLightCurve )
        {
            cout << "VLightCurvePlotter::plot: no light curve histogram found with name " << htitle << endl;
            return fCanvasLightCurve;
        }
    }
    
    ////////////////////////////////////////////////////////////////////////
    // flux and upper flux plotting
    
    fLightCurveGraph = new TGraphAsymmErrors( 0 );
    
    // loop over all values in flux data vector
    unsigned int z = 0;
    for( unsigned int i = 0; i < fFluxDataVector.size(); i++ )
    {
        // check if this is a good time fime
        if( fFluxDataVector[i].fTimeMask_open == true )
        {
            continue;
        }
        // x-axis: MJD or phase
        double iTimeAxisValue_mean     =  fFluxDataVector[i].fMJD - fPlotting_MJD_offset_days;
        double iTimeAxisValue_error_lo = fFluxDataVector[i].fMJD - fFluxDataVector[i].fMJD_Start;
        double iTimeAxisValue_error_hi = fFluxDataVector[i].fMJD_Stop - fFluxDataVector[i].fMJD;
        if( fLightCurveTimeAxis_is_OrbitalPhase )
        {
            iTimeAxisValue_mean = fFluxDataVector[i].fOrbitalPhase;
            iTimeAxisValue_error_lo = fFluxDataVector[i].fOrbitalPhase - fFluxDataVector[i].fOrbitalPhase_Start;
            if( iTimeAxisValue_error_lo < 0. )
            {
                iTimeAxisValue_error_lo = 0.;
            }
            iTimeAxisValue_error_hi = fFluxDataVector[i].fOrbitalPhase_Stop - fFluxDataVector[i].fOrbitalPhase;
            // not really true, but probably good enough
            if( iTimeAxisValue_error_hi + fFluxDataVector[i].fOrbitalPhase > 1. )
            {
                iTimeAxisValue_error_hi = 1. - fFluxDataVector[i].fOrbitalPhase;
            }
        }
        // don't plot 'error' in MJD
        if( !fLightCurveTimeAxis_is_OrbitalPhase && iMaxMJDError > 0. && iTimeAxisValue_error_lo > iMaxMJDError )
        {
            iTimeAxisValue_error_lo = 0.;
        }
        // don't plot 'error' in MJD
        if( !fLightCurveTimeAxis_is_OrbitalPhase && iMaxMJDError > 0. && iTimeAxisValue_error_hi > iMaxMJDError )
        {
            iTimeAxisValue_error_hi = 0.;
        }
        
        /////////////
        // plot fluxes or confidence intervals
        if( fFluxDataVector[i].isSignificantDataPoint() || fIgnoreUpperLimits )
        {
            fLightCurveGraph->SetPoint( z, iTimeAxisValue_mean, fFluxDataVector[i].fFlux_Rolke );
            if( fPlot_observing_intervals_as_errors )
            {
                fLightCurveGraph->SetPointEXlow( z, iTimeAxisValue_error_lo );
                fLightCurveGraph->SetPointEXhigh( z, iTimeAxisValue_error_hi );
            }
            else
            {
                fLightCurveGraph->SetPointEXlow( z, 0. );
                fLightCurveGraph->SetPointEXhigh( z, 0. );
            }
            fLightCurveGraph->SetPointEYlow( z, fFluxDataVector[i].fFluxCI_lo_1sigma );
            fLightCurveGraph->SetPointEYhigh( z, fFluxDataVector[i].fFluxCI_up_1sigma );
            z++;
        }
        /////////////
        // plot upper flux limits
        else if( fFluxDataVector[i].fFluxUL > 0. )
        {
            TArrow* fUL = new TArrow( iTimeAxisValue_mean, fFluxDataVector[i].fFluxUL, iTimeAxisValue_mean, fFluxDataVector[i].fFluxUL - 0.05 * hLightCurve->GetMaximum(), 0.01, "|-|>" );
            setArrowPlottingStyle( fUL );
            fUL->Draw();
        }
    }
    if( fLightCurveGraph->GetN() > 0 )
    {
        setGraphPlottingStyle( ( TGraph* )fLightCurveGraph );
        fLightCurveGraph->Draw( iPlottingOption.c_str() );
    }
    
    return fCanvasLightCurve;
}


/*

    get title of flux axis (depend on energy ranges)

*/
string VLightCurvePlotter::getLightCurveFluxAxisTitle()
{

    if( fPlotting_FluxAxis_title == "not_set" )
    {
        char hname[1000];
        if( fMaxEnergy_TeV < 0. || fMaxEnergy_TeV > 1.e2 - 1. )
        {
            // determine number of decimal places (do not allow more than three)
            if( ( int )( fMinEnergy_TeV * 10. ) % 10 == 0 && ( int )( fMinEnergy_TeV * 100. ) % 10 == 0 && ( int )( fMinEnergy_TeV * 1000. ) % 10 == 0 )
            {
                sprintf( hname, "Flux (E>%.0f TeV) [cm^{-2}s^{-1}]", fMinEnergy_TeV );
            }
            else if( ( int )( fMinEnergy_TeV * 100. ) % 10 == 0 && ( int )( fMinEnergy_TeV * 1000. ) % 10 == 0 )
            {
                sprintf( hname, "Flux (E>%.1f TeV) [cm^{-2}s^{-1}]", fMinEnergy_TeV );
            }
            else if( ( int )( fMinEnergy_TeV * 1000. ) % 10 == 0 )
            {
                sprintf( hname, "Flux (E>%.2f TeV) [cm^{-2}s^{-1}]", fMinEnergy_TeV );
            }
            else
            {
                sprintf( hname, "Flux (E>%.3f TeV) [cm^{-2}s^{-1}]", fMinEnergy_TeV );
            }
        }
        else if( fMaxEnergy_TeV > 0. && fMaxEnergy_TeV < 1.e2 )
        {
            // VHE data
            if( fMinEnergy_TeV > 1.e-2 )
            {
                sprintf( hname, "F(%d GeV < E < %0.1f TeV) [cm^{-2}s^{-1}]", ( int )( fMinEnergy_TeV * 1.e3 ), fMaxEnergy_TeV );
            }
            // X-ray data
            else if( fMinEnergy_TeV < 1.e-5 )
            {
                sprintf( hname, "F(%.1f keV < E < %0.1f keV) [cm^{-2}s^{-1}]", fMinEnergy_TeV * 1.e6, fMaxEnergy_TeV * 1.e6 );
            }
        }
        
        fPlotting_FluxAxis_title = hname;
    }
    
    return fPlotting_FluxAxis_title;
}

/*


*/
void VLightCurvePlotter::setLightCurveTimeAxis_to_OrbitalPhase( bool iOrbitalPhase )
{
    fLightCurveTimeAxis_is_OrbitalPhase = iOrbitalPhase;
    setLightCurveTimeAxis();
}
