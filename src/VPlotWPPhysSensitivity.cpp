/*! \class VPlotWPPhysSensitivity

     plotting class for CTA Sensitivities and Performances

*/

#include "VPlotWPPhysSensitivity.h"

VPlotWPPhysSensitivity::VPlotWPPhysSensitivity()
{
    fIRF = 0;
    bPlotNoLegend = false;
    setEnergyRange_Lin_TeV();
    setCrabSpectraFile();
    setPlotCTARequirements();
    setPlotCrabLines();
    setMaximumAllowedEnergyBias();
    setCurrentInstrumentFile();
    setCurrentInstrumentPlotVector();
    
    fSensitivityFOM = -1.;
    fSensitivityFOM_error = -1.;
    
    fUseIntegratedSensitivityForOffAxisPlots = false;
    
    setNorthSouthComparision( false );
}

VPlotWPPhysSensitivity::~VPlotWPPhysSensitivity()
{
    if( fIRF )
    {
        delete fIRF;
    }
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] )
        {
            delete fData[i];
        }
    }
}

void VPlotWPPhysSensitivity::reset()
{
    fData.clear();
}

bool VPlotWPPhysSensitivity::addDataSet( VSiteData* iData )
{
    fData.push_back( new VSiteData() );
    fData.back()->fSiteName = iData->fSiteName;
    fData.back()->fArray = iData->fArray;
    fData.back()->fObservationTime_s = iData->fObservationTime_s;
    fData.back()->fCameraOffset_deg = iData->fCameraOffset_deg;
    fData.back()->fSiteFile_Emin = iData->fSiteFile_Emin;
    fData.back()->fSiteFile_Emax = iData->fSiteFile_Emax;
    
    fData.back()->fPlottingColor = iData->fPlottingColor;
    fData.back()->fPlottingLineStyle = iData->fPlottingLineStyle;
    fData.back()->fPlottingFillStyle = iData->fPlottingFillStyle;
    fData.back()->fPlottingMarkerStyle = iData->fPlottingMarkerStyle;
    fData.back()->fLegend = iData->fLegend;
    
    return true;
}

bool VPlotWPPhysSensitivity::addDataSet( string iAnalysis, string iSubArray,
        double iObservationTime_s, double iOffset_deg,
        string iLegend, int iColor, int iLineStyle, int iFillStyle )
{
    VSiteData i_temp;
    i_temp.fSiteName = iAnalysis;
    i_temp.fArray.push_back( iSubArray );
    i_temp.fObservationTime_s.push_back( iObservationTime_s );
    i_temp.fCameraOffset_deg.push_back( iOffset_deg );
    
    i_temp.fPlottingColor.push_back( iColor );
    i_temp.fPlottingLineStyle.push_back( iLineStyle );
    i_temp.fPlottingFillStyle.push_back( iFillStyle );
    i_temp.fPlottingMarkerStyle.push_back( 21 );
    i_temp.fLegend.push_back( iLegend );
    
    return addDataSet( &i_temp );
}

/*
    loop over list of sites read from a plotlist file
    and fill the data set vector

*/
bool VPlotWPPhysSensitivity::addDataSets( string iDataSettxtFile, string iDirectionString )
{
    std::cout << "VPlotWPPhysSensitivity::addDataSets " << std::endl;
    unsigned int z_site = 0;
    for( ;; )
    {
        fData.push_back( new VSiteData() );
        if( !fData.back()->addDataSet( iDataSettxtFile, z_site, iDirectionString ) )
        {
            fData.pop_back();
            break;
        }
        z_site++;
    }
    
    if( fData.size() == 0 )
    {
        return false;
    }
    
    return true;
}


/*

   plot IRFs for all data sets

*/
bool VPlotWPPhysSensitivity::plotIRF( string iPrint,
                                      double iEffAreaMin, double iEffAreaMax,
                                      double iAngularResolutionMin, double iAngularResolutionMax,
                                      double iEnergyResolutionMin, double iEnergyResolutionMax,
                                      TPad* iEffAreaPad, TPad* iAngResPad, TPad* iEResPad,
                                      bool iPlotEnergyBias, bool iLogAngRes )
{
    fIRF = new VPlotInstrumentResponseFunction();
    
    fIRF->setCanvasSize( 400, 400 );
    fIRF->setPlottingAxis( "energy_Lin", "X", true, fMinEnergy_TeV, fMaxEnergy_TeV, "energy [TeV]" );
    fIRF->setPlottingAxis( "effarea_Lin", "X", true, iEffAreaMin, iEffAreaMax );
    fIRF->setPlottingAxis( "energyresolution_Lin", "X", false, 0., 0.7 );
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] )
        {
            for( unsigned int j = 0; j < fData[i]->fSiteFileName.size(); j++ )
            {
                if( fData[i]->fSiteFile_exists[j] )
                {
                    fIRF->addInstrumentResponseData( fData[i]->fSiteFileName[j], 20., fData[i]->fCameraOffset_deg[j],
                                                     0, 2.4, 200, "A_MC", fData[i]->fPlottingColor[j],
                                                     fData[i]->fPlottingLineStyle[j],
                                                     fData[i]->fPlottingMarkerStyle[j], 0.75,
                                                     fData[i]->fSiteFile_Emin[j], fData[i]->fSiteFile_Emax[j] );
                }
            }
        }
    }
    
    char hname[2000];
    
    ////////////////////////////
    // effective areas
    TCanvas* c = fIRF->plotEffectiveArea( iEffAreaMin, iEffAreaMax, iEffAreaPad );
    // plotLegend( c, true );
    if( fRequirementsString.size() > 0
            && fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
    {
        fPlotCTARequirements[0]->plotRequirement_EffectiveArea( c );
    }
    if( iPrint.size() > 0 )
    {
        sprintf( hname, "%s-EffArea.pdf", iPrint.c_str() );
        if( c )
        {
            c->Print( hname );
        }
    }
    ////////////////////////////
    // angular resolution (68%)
    if( iLogAngRes )
    {
        fIRF->setPlottingAxis( "angularesolution_Log", "X", true, 0.009 );
        c = fIRF->plotAngularResolution( "energy", "68", iAngularResolutionMin, iAngularResolutionMax, iAngResPad, true );
    }
    else
    {
        fIRF->setPlottingAxis( "angularesolution_Lin", "X", false, 0.002 );
        c = fIRF->plotAngularResolution( "energy", "68", iAngularResolutionMin, iAngularResolutionMax, iAngResPad, false );
    }
    // plotLegend( c, false );
    if( fRequirementsString.size() > 0
            && fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
    {
        fPlotCTARequirements[0]->plotRequirement_AngularResolution( c );
    }
    if( iPrint.size() > 0 )
    {
        sprintf( hname, "%s-AngRes.pdf", iPrint.c_str() );
        if( c )
        {
            c->Print( hname );
        }
    }
    // angular resolution (80%)
    // c = fIRF->plotAngularResolution( "energy", "80", -1.e99, 0 );
    // plotLegend( c, false );
    ////////////////////////////
    // energy resolution
    c = fIRF->plotEnergyResolution( iEnergyResolutionMin, iEnergyResolutionMax, iEResPad, false );
    // plotLegend( c, false );
    if( fRequirementsString.size() > 0
            && fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
    {
        fPlotCTARequirements[0]->plotRequirement_EnergyResolution( c );
    }
    if( iPrint.size() > 0 )
    {
        sprintf( hname, "%s-ERes.pdf", iPrint.c_str() );
        if( c )
        {
            c->Print( hname );
        }
    }
    // energy bias
    if( iPlotEnergyBias )
    {
        c = fIRF->plotEnergyReconstructionBias( "mean", 0.8, 1.5 );
        // plotLegend( c, false );
        if( iPrint.size() > 0 )
        {
            sprintf( hname, "%s-EBias.pdf", iPrint.c_str() );
            if( c )
            {
                c->Print( hname );
            }
        }
    }
    
    return true;
}

void VPlotWPPhysSensitivity::initialProjectedSensitivityPlots( bool iIncludeLowestEnergy )
{
    fProjectionEnergy_min_logTeV.clear();
    fProjectionEnergy_max_logTeV.clear();
    // (hard coded energies here...not good)
    if( !fUseIntegratedSensitivityForOffAxisPlots )
    {
        //		fProjectionEnergy_min_logTeV.push_back( log10( 60.0 ) );
        //		fProjectionEnergy_max_logTeV.push_back( log10( 70.0 ) );
        //		fProjectionEnergy_min_logTeV.push_back( log10( 30.0 ) );
        //		fProjectionEnergy_max_logTeV.push_back( log10( 40.0 ) );
        fProjectionEnergy_min_logTeV.push_back( log10( 50. ) );
        fProjectionEnergy_max_logTeV.push_back( log10( 100. ) );
        fProjectionEnergy_min_logTeV.push_back( log10( 5. ) );
        fProjectionEnergy_max_logTeV.push_back( log10( 10. ) );
        /*fProjectionEnergy_min_logTeV.push_back( log10( 0.5 ) );
        fProjectionEnergy_max_logTeV.push_back( log10( 0.8 ) );
        if( iIncludeLowestEnergy )
        {
            fProjectionEnergy_min_logTeV.push_back( log10( 0.05 ) );
            fProjectionEnergy_max_logTeV.push_back( log10( 0.08 ) );
        } */
    }
    // integrated sensitivity
    else
    {
        fProjectionEnergy_min_logTeV.push_back( log10( 0.03 ) ); // choose 128 GeV the be at the lower end of the corresponding bin on the log axis)
        fProjectionEnergy_max_logTeV.push_back( log10( 0.03 ) );
        //		fProjectionEnergy_min_logTeV.push_back( log10( 0.030 ) ); // choose 128 GeV the be at the lower end of the corresponding bin on the log axis)
        //		fProjectionEnergy_max_logTeV.push_back( log10( 0.030 ) );
        //		fProjectionEnergy_min_logTeV.push_back( log10( 0.032 ) ); // choose 128 GeV the be at the lower end of the corresponding bin on the log axis)
        //		fProjectionEnergy_max_logTeV.push_back( log10( 0.032 ) );
    }
    // graphs
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fProjectionSensitivityvsCameraOffset.find( fData[i]->fReferenceSiteName ) == fProjectionSensitivityvsCameraOffset.end() )
        {
            char hname[200];
            for( unsigned int j = 0; j < fProjectionEnergy_min_logTeV.size(); j++ )
            {
                fProjectionSensitivityvsCameraOffset[fData[i]->fReferenceSiteName].push_back( new TGraphAsymmErrors( ) );
                sprintf( hname, "fProjectionSensitivityvsCameraOffset_%s_%d", fData[i]->fReferenceSiteName.c_str(), j );
                fProjectionSensitivityvsCameraOffset[fData[i]->fReferenceSiteName].back()->SetName( hname );
                fProjectionSensitivityvsCameraOffset[fData[i]->fReferenceSiteName].back()->SetTitle( "" );
                
                setGraphPlottingStyle( fProjectionSensitivityvsCameraOffset[fData[i]->fReferenceSiteName].back(), i + 1, 1., 20, 1.5 );
            }
        }
    }
}

/*

     off-axis sensitivity plots

     input graph (g) with sensitivities
     (can be differential or integral sensitivity)

*/
void VPlotWPPhysSensitivity::fillProjectedSensitivityPlot( unsigned int iDataSet, TGraphAsymmErrors* g )
{
    if( !g )
    {
        return;
    }
    if( iDataSet >= fData.size() )
    {
        return;
    }
    
    if( fProjectionSensitivityvsCameraOffset.find( fData[iDataSet]->fReferenceSiteName ) == fProjectionSensitivityvsCameraOffset.end() )
    {
        return;
    }
    
    if( g )
    {
        VHistogramUtilities h;
        //////////////////////////////
        // loop over all energy bins
        for( unsigned int i = 0; i < fProjectionEnergy_min_logTeV.size(); i++ )
        {
            if( i < fProjectionSensitivityvsCameraOffset[fData[iDataSet]->fReferenceSiteName].size()
                    && fProjectionSensitivityvsCameraOffset[fData[iDataSet]->fReferenceSiteName][i] )
            {
                int z = fProjectionSensitivityvsCameraOffset[fData[iDataSet]->fReferenceSiteName][i]->GetN();
                int iBin_min = h.findBinInGraph( g, fProjectionEnergy_min_logTeV[i] );
                int iBin_max = h.findBinInGraph( g, fProjectionEnergy_max_logTeV[i] );
                double i_m = 0.;
                double i_m_lE = 0.;
                double i_m_hE = 0.;
                double i_mz = 0.;
                // average over energies
                for( int b = iBin_min; b <= iBin_max; b++ )
                {
                    if( b >= 0 )
                    {
                        // check that this is not below or above the highest bin
                        double x = 0.;
                        double y = 0.;
                        g->GetPoint( b, x, y );
                        
                        i_m += y;
                        i_m_lE += g->GetErrorYlow( b ) * g->GetErrorYlow( b );
                        i_m_hE += g->GetErrorYhigh( b ) * g->GetErrorYhigh( b );
                        i_mz++;
                    }
                }
                if( i_mz > 0 && i_m > 0. && fData[iDataSet]->fCameraOffset_deg.size() > 0 )
                {
                    fProjectionSensitivityvsCameraOffset[fData[iDataSet]->fReferenceSiteName][i]->SetPoint( z, fData[iDataSet]->fCameraOffset_deg[0], i_mz / i_m );
                    fProjectionSensitivityvsCameraOffset[fData[iDataSet]->fReferenceSiteName][i]->SetPointEYhigh( z, sqrt( i_m_hE / i_mz ) );
                    fProjectionSensitivityvsCameraOffset[fData[iDataSet]->fReferenceSiteName][i]->SetPointEYlow( z, sqrt( i_m_lE / i_mz ) );
                }
            }
        }
    }
}

/*

    plot relative off-axis sensitivity

*/
TCanvas* VPlotWPPhysSensitivity::plotProjectedSensitivities( TCanvas* c, double iMaxOffSet, int iColor )
{
    bool bUpsideDown = false;
    
    TCanvas* cC = 0;
    if( c )
    {
        cC = c;
    }
    else
    {
        cC = new TCanvas( "cCProjection", "projected off-axis sensitivities", 800, 100, 600, 600 );
        cC->SetGridx( 0 );
        cC->SetGridy( 0 );
        cC->SetLeftMargin( 0.13 );
        cC->Draw();
        TH1D* hnull = new TH1D( "hnull", "", 10, 0.0, iMaxOffSet + 0.5 );
        hnull->SetStats( 0 );
        hnull->SetXTitle( "angular distance to field of view center [deg]" );
        hnull->SetYTitle( "Sensitivity relative to field of view center" );
        if( bUpsideDown )
        {
            hnull->SetMinimum( 0.51 );
            hnull->SetMaximum( 20. );
        }
        else
        {
            hnull->SetMinimum( 0.01 );
            hnull->SetMaximum( 1.1 );
        }
        hnull->Draw();
        TLine* iLi = new TLine( hnull->GetXaxis()->GetXmin(), 1., hnull->GetXaxis()->GetXmax(), 1. );
        iLi->Draw();
        TLine* iL2 = new TLine( hnull->GetXaxis()->GetXmin(), 0.5, hnull->GetXaxis()->GetXmax(), 0.5 );
        iL2->SetLineStyle( 2 );
        iL2->Draw();
    }
    cC->cd();
    
    TLegend* iL = new TLegend( 0.55, 0.65, 0.88, 0.88 );
    char hname[200];
    
    map< string, vector< TGraphAsymmErrors* > >::iterator i_fProjectionSensitivityvsCameraOffset_iter;
    
    for( i_fProjectionSensitivityvsCameraOffset_iter = fProjectionSensitivityvsCameraOffset.begin();
            i_fProjectionSensitivityvsCameraOffset_iter != fProjectionSensitivityvsCameraOffset.end();
            i_fProjectionSensitivityvsCameraOffset_iter++ )
    {
        for( unsigned int i = 0; i < fProjectionEnergy_min_logTeV.size(); i++ )
        {
            if( i < i_fProjectionSensitivityvsCameraOffset_iter->second.size()
                    && i_fProjectionSensitivityvsCameraOffset_iter->second[i] )
            {
                TGraphAsymmErrors* iGraph = i_fProjectionSensitivityvsCameraOffset_iter->second[i];
                cout << "Projected sensitivity" << endl;
                iGraph->Print();
                cout << "---------------------" << endl;
                // normalize graphs to average of first two points
                double x = 0.;
                double y = 0.;
                double y_norm = 0.;
                iGraph->GetPoint( 0, x, y_norm );
                iGraph->GetPoint( 0, x, y );          // NOTE!!!! (normalized to first point)
                y_norm = 0.5 * ( y + y_norm );
                double y_normE_low  = sqrt( iGraph->GetErrorYlow( 0 ) * iGraph->GetErrorYlow( 0 ) );
                double y_normE_high = sqrt( iGraph->GetErrorYhigh( 0 ) * iGraph->GetErrorYhigh( 0 ) );
                
                for( int p = iGraph->GetN(); p > 0; p-- )
                {
                    iGraph->GetPoint( p - 1, x, y );
                    if( y_norm > 0. )
                    {
                        if( bUpsideDown )
                        {
                            iGraph->SetPoint( p, x, y_norm / y );
                        }
                        else
                        {
                            iGraph->SetPoint( p, x, y / y_norm );
                        }
                    }
                    else
                    {
                        iGraph->SetPoint( p, x, 0. );
                    }
                    iGraph->SetPointEYhigh( p, VMathsandFunctions::getRatioError( y, y_norm, iGraph->GetErrorYhigh( p ), y_normE_high ) );
                    iGraph->SetPointEYlow( p, VMathsandFunctions::getRatioError( y, y_norm, iGraph->GetErrorYlow( p ), y_normE_low ) );
                    if( y_norm > 0. )
                    {
                        cout << x << "\t" << y / y_norm << "\t" << endl;
                    }
                }
                iGraph->SetPoint( 0, 0., 1. );
                setGraphPlottingStyle( iGraph, i + 1, 2., 21, 1.5 );
                iGraph->Draw( "pl" );
                if( fProjectionEnergy_min_logTeV[i] >= 0. )
                {
                    sprintf( hname, "%d-%d TeV", ( int )( TMath::Power( 10., fProjectionEnergy_min_logTeV[i] ) + 0.5 ),
                             ( int )( TMath::Power( 10., fProjectionEnergy_max_logTeV[i] ) + 0.5 ) );
                }
                else
                {
                    sprintf( hname, "%.2f-%.2f TeV", TMath::Power( 10., fProjectionEnergy_min_logTeV[i] ),
                             TMath::Power( 10., fProjectionEnergy_max_logTeV[i] ) );
                }
                iL->AddEntry( iGraph, hname, "lp" );
            }
        }
    }
    // plot the legend only if there is more than one energy bin
    /*if( iColor < 0 && fProjectionEnergy_min_logTeV.size() > 1 )
    {
        iL->Draw();
    }
    // plot some text
    TText* iT = new TText( 1., 0.1, "Off-axis sensitivity" );
    iT->Draw(); */
    
    return cC;
}

/*

    plot ratio of sensitivity curves

    - to required sensitivity
    - to goal sensitivity
    - to another sensitivity curve

    iRatioSelector = 0: ratio to requirement
    iRatioSelector = 1: ratio to goal
    iRatioSelector !0, 1 ratio to first sensitivity in list
*/
bool VPlotWPPhysSensitivity::plotSensitivityRatio( string iPrint,
        double ymin, double ymax,
        unsigned int iRatioSelector,
        TPad* iSensRatio,
        unsigned int iRatioGraphCounter )
{
    char hname[200];
    char htitle[200];
    sprintf( hname, "cSensRatio_%d", iRatioSelector );
    if( iRatioSelector == 1 )
    {
        sprintf( htitle, "ratio of sensitivities (to goal)" );
    }
    else if( iRatioSelector == 0 )
    {
        sprintf( htitle, "ratio of sensitivities (to requirement)" );
    }
    else
    {
        sprintf( htitle, "ratio of sensitivities" );
    }
    TCanvas* cSensRatio = 0;
    if( !iSensRatio )
    {
        cSensRatio = new TCanvas( hname, htitle, 200, 200, 900, 600 );
        cSensRatio->SetGridy( 0 );
        cSensRatio->SetGridx( 0 );
    }
    else
    {
        cSensRatio = ( TCanvas* )iSensRatio;
    }
    cSensRatio->cd();
    
    // plot null histogram
    sprintf( hname, "hSensRatio_%d", iRatioSelector );
    TH1D* hNull = new TH1D( hname, "", 10, log10( fMinEnergy_TeV ), log10( fMaxEnergy_TeV ) );
    hNull->SetStats( 0 );
    hNull->SetXTitle( "log_{10} energy [TeV]" );
    if( iRatioSelector == 1 )
    {
        hNull->SetYTitle( "Sensitivity ratio to goal sensitivity" );
    }
    else if( iRatioSelector == 0 )
    {
        hNull->SetYTitle( "Sensitivity ratio to required sensitivity" );
    }
    else
    {
        hNull->SetYTitle( "Sensitivity ratio" );
    }
    hNull->SetMinimum( ymin );
    hNull->SetMaximum( ymax );
    hNull->Draw( "" );
    hNull->Draw( "AH" );
    
    plot_nullHistogram( cSensRatio, hNull, true , false, 1.2, fMinEnergy_TeV, fMaxEnergy_TeV );
    
    TLine* iL = new TLine( log10( fMinEnergy_TeV ), 1., log10( fMaxEnergy_TeV ), 1. );
    iL->SetLineStyle( 2 );
    iL->SetLineWidth( 3 );
    iL->Draw();
    
    // North/South comparision with Legend
    TLegend* iLL = 0;
    if( fNorthSouthComparision )
    {
        if( iRatioGraphCounter == 9998 )
        {
            iLL = new TLegend( 0.12, 0.75, 0.88, 0.88 );
        }
        else
        {
            iLL = new TLegend( 0.12, 0.15, 0.88, 0.3 );
        }
        //iLL->SetFillColorAlpha( 0, 0.1 );
    }
    // PPUT calculation
    VPPUTValues i_PPUT;
        
    
    ////////////////////////////////////////////////////////////////////
    // loop over all data sets and divide it by the selected one
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( !fData[i] )
        {
            continue;
        }
        if( iRatioSelector > 1 && i == iRatioGraphCounter )
        {
            continue;
        }
        
        ///////////////////////////////////////////////////////////////////////
        // select graph to which the sensitivity ratio is calculated to
        
        TGraphAsymmErrors* gRelG = 0;
        // relative to requirement
        if( !fNorthSouthComparision && iRatioGraphCounter == 998 && fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
            //                if( !fNorthSouthComparision && fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
        {
            gRelG = ( TGraphAsymmErrors* )fPlotCTARequirements[0]->getRequiredDifferentalSensitivity();
        }
        // relative to goal
        else if( !fNorthSouthComparision && iRatioGraphCounter == 998 && fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
        {
            gRelG = 0;
        }
        // hard wired values for NS comparision
        else if( fNorthSouthComparision )
        {
            // do not plot first and third ratio
            if( i == 0 || i == 2 )
            {
                continue;
            }
            if( fData.size() == 4 )
            {
                if( iRatioGraphCounter == 9998 && i == 1 && fPlotCTARequirements.size() > 1 && fPlotCTARequirements[1] )
                {
                    gRelG = ( TGraphAsymmErrors* )fPlotCTARequirements[1]->getRequiredDifferentalSensitivity();
                }
                else if( iRatioGraphCounter == 9998 && i == 3 && fPlotCTARequirements.size() > 1 && fPlotCTARequirements[1] )
                {
                    gRelG = ( TGraphAsymmErrors* )fPlotCTARequirements[0]->getRequiredDifferentalSensitivity();
                }
                else
                {
                    if( i == 1 && fData[0] && fData[0]->fGraphSensitivity[0] )
                    {
                        gRelG = fData[0]->fGraphSensitivity[0];
                    }
                    else if( i == 3 && fData[2] && fData[2]->fGraphSensitivity[0] )
                    {
                        gRelG = fData[2]->fGraphSensitivity[0];
                    }
                }
            }
        }
        // relative to another graph
        else
        {
            // prefer interpolated graph
            if( fData.size() > iRatioGraphCounter
                    && fData[iRatioGraphCounter]
                    && fData[iRatioGraphCounter]->fGraphSensitivityInterPolated )
            {
                gRelG = fData[iRatioGraphCounter]->fGraphSensitivityInterPolated;
            }
            else if( fData.size() > iRatioGraphCounter
                     && fData[iRatioGraphCounter]
                     && fData[iRatioGraphCounter]->fGraphSensitivity[0] )
            {
                gRelG = fData[iRatioGraphCounter]->fGraphSensitivity[0];
            }
        }
        if( !gRelG )
        {
            return false;
        }
        if( fNorthSouthComparision )
        {
            iL->SetLineColor( 1 );
        }
        else
        {
            iL->SetLineColor( gRelG->GetLineColor() );
        }

        for( unsigned int j = 0; j < fData[i]->fGraphSensitivity.size(); j++ )
        {
            //
            if( fData[i]->fGraphSensitivity[j] )
            {
                TGraphAsymmErrors* g = new TGraphAsymmErrors();
                if( gRelG )
                {
                    VHistogramUtilities::divide( g, fData[i]->fGraphSensitivity[j], gRelG );
                }
                if( g->GetN() > 0 )
                {
                    setGraphPlottingStyle( g, fData[i]->fGraphSensitivity[j]->GetMarkerColor(), 1.,
                                           fData[i]->fPlottingMarkerStyle[j], 1.1,
                                           fData[i]->fPlottingFillStyle[j],
                                           fData[i]->fPlottingLineStyle[j] );
                    g->Draw( "p" );
                    if( fData[i]->fLegend.size() > 0 )
                    {
                        i_PPUT.add( fData[i]->fLegend[0], g );
                    }
                    
                    if( fNorthSouthComparision && fData[i]->fLegend.size() > 0 )
                    {
                        ostringstream iTitle;
                        iTitle << fData[i]->fLegend[0] << " relative to ";
                        if( iRatioGraphCounter == 9998 )
                        {
                            iTitle << "required sensitivity";
                        }
                        else if( i == 1 && fData[0]->fLegend.size() > 0 )
                        {
                            iTitle << fData[0]->fLegend[0];
                        }
                        else if( i == 3 && fData[2]->fLegend.size() > 0 )
                        {
                            iTitle << fData[2]->fLegend[0];
                        }
                        iLL->AddEntry( g, iTitle.str().c_str(), "pl" );
                    }
                }
            }
        }
        // TMPTMP add current instruments
        vector< TGraph* > iCurrentInstruments = plotCurrentInstruments( 0 );
        for( unsigned int i = 0; i < iCurrentInstruments.size(); i++ )
        {
            //TMPTMP
            //                    continue;
            if( !iCurrentInstruments[i] )
            {
                continue;
            }
            if( gRelG )
            {
                TGraphAsymmErrors* g = new TGraphAsymmErrors();
                g->SetLineColor( iCurrentInstruments[i]->GetLineColor() );
                g->SetLineStyle( iCurrentInstruments[i]->GetLineStyle() );
                g->SetLineWidth( iCurrentInstruments[i]->GetLineWidth() );
                g->SetLineWidth( 2 );
                int z = 0;
                double x = 0.;
                double y = 0.;
                double yc = 0.;
                // get man and min energy of current instrument curve
                double x_min = 0.;
                double x_max = 0.;
                iCurrentInstruments[i]->GetPoint( 0, x_min, y );
                iCurrentInstruments[i]->GetPoint( iCurrentInstruments[i]->GetN() - 1, x_max, y );
                for( int p = 0; p < gRelG->GetN(); p++ )
                {
                    gRelG->GetPoint( p, x, y );
                    if( x < x_min || x > x_max )
                    {
                        continue;
                    }
                    yc = iCurrentInstruments[i]->Eval( x );
                    if( yc > 0. )
                    {
                        g->SetPoint( z, x, y / yc );
                        z++;
                    }
                }
                if( z > 0 )
                {
                    g->Draw( "l" );
                }
            }
        }
    }
    i_PPUT.printLatexTable();
    if( cSensRatio )
    {
        // TEMPTEMPTEMP
       	// plotLegend( cSensRatio, false, false );
        
        if( iLL )
        {
            iLL->Draw();
        }
        // print results
        if( iPrint.size() > 0 )
        {
            char hname[2000];
            char hnameC[2000];
            if( iRatioSelector == 1 )
            {
                sprintf( hname, "%s-SensitivityRatioGoal.pdf", iPrint.c_str() );
                sprintf( hnameC, "%s_SensitivityRatioGoal.root", iPrint.c_str() );
            }
            else if( iRatioSelector == 0 )
            {
                sprintf( hname, "%s-SensitivityRatio.pdf", iPrint.c_str() );
                sprintf( hnameC, "%s_SensitivityRatio.root", iPrint.c_str() );
            }
            else
            {
                sprintf( hnameC, "%s_SensitivityRatio_%d.root", iPrint.c_str(), iRatioSelector );
                sprintf( hname, "%s-SensitivityRatio_%d.pdf", iPrint.c_str(), iRatioSelector );
            }
            if( cSensRatio )
            {
                cSensRatio->Print( hname );
                cSensRatio->Print( hnameC );
            }
        }
    }
    
    return true;
}


/*

    plot sensitivities and data rates for different data sets

*/
bool VPlotWPPhysSensitivity::plotSensitivity( string iPrint,
        double iMinSensitivity,
        double iMaxSensitivity,
        string iUnit,
        TPad* iSensitivityPad, TPad* iBckPad )
{
    TCanvas* cIntSens = 0;
    TCanvas* cSens = 0;
    TCanvas* cSensInter = ( TCanvas* )iSensitivityPad;
    TCanvas* cBck = ( TCanvas* )iBckPad;
    
    bool iIncludeLowestEnergy = true;
    
    initialProjectedSensitivityPlots( iIncludeLowestEnergy );
    
    ////////////////////////////////////////////////////////
    // loop over all data sets
    unsigned int z = 0;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( !fData[i] || !fData[i]->checkIntegrity() )
        {
            continue;
        }
        // loop over all sites
        for( unsigned j = 0; j < fData[i]->fSiteFileName.size(); j++ )
        {
            cout << "NOW AT " << fData[i]->fSiteFileName[j] << ", " << fData[i]->fCameraOffset_deg[j] << endl;
            VSensitivityCalculator* a = new VSensitivityCalculator();
            if( fMaximumAllowedEnergyBias > 0. )
            {
                cout << "\t max allowed energy bias: " << fMaximumAllowedEnergyBias << endl;
            }
            a->setMaximumEnergyBias( fMaximumAllowedEnergyBias );
            a->setMonteCarloParametersCTA_MC( fData[i]->fSiteFileName[j], fData[i]->fCameraOffset_deg[j], fCrabSpectraFile, fCrabSpectraID );
            a->setPlotCrabLines( bPlotCrabLines );
            a->setEnergyRange_Lin( fMinEnergy_TeV, fMaxEnergy_TeV );
            a->setPlotCanvasSize( 900, 600 );
            a->setPlottingStyle( fData[i]->fPlottingColor[j], fData[i]->fPlottingLineStyle[j], 2., 21, 1., fData[i]->fPlottingFillStyle[j] );
            if( iUnit == "ENERGY" )
            {
                a->setFluxRange_ENERG( iMinSensitivity, iMaxSensitivity );
            }
            else if( iUnit == "CU" )
            {
                a->setFluxRange_CU( iMinSensitivity, iMaxSensitivity );
            }
            TCanvas* c_temp = 0;
            // cSens name = cSensitivity (default from VSensitivityCalculator::fPlot_CanvasName )
            c_temp = a->plotDifferentialSensitivityvsEnergyFromCrabSpectrum( cSens, "CTA-PHYS", iUnit, 0.2,
                     fData[i]->fSiteFile_Emin[j],
                     fData[i]->fSiteFile_Emax[j] );
            if( c_temp )
            {
                cSens = c_temp;
            }
            // integral sensitivity
            /*TCanvas *c_tempInter = 0;
            // cIntSens name (changing VSensitivityCalculator::fPlot_CanvasName default value)
            a->setPlotCanvasName("cIntegratedSensitivity", "integral sensitivity vs energy" );
            a->setFluxRange_ENERG(  2.5e-14, 2.5e-11 );
            c_tempInter = a->plotIntegralSensitivityvsEnergyFromCrabSpectrum( cIntSens, "CTA-PHYS", iUnit,
              		                                                  fData[i]->fSiteFile_Emin[j], fData[i]->fSiteFile_Emax[j] );
            if( c_tempInter ) cIntSens = c_tempInter; */
            c_temp = a->plotSignalBackgroundRates( cBck, ( z == 0 ), 2.e-7, 14. ); // plot also protons and electrons
            if( c_temp )
            {
                cBck = c_temp;
            }
            fData[i]->fGraphSensitivity[j] = ( TGraphAsymmErrors* )a->getSensitivityGraph();
            z++;
        }
        ////////////////////////////////////////////////////////
        // plot a second window with interpolated sensitivities
        TGraphAsymmErrors* iGraphSensitivity = 0;
        if( fData[i]->fSiteFileName.size() > 1 )
        {
            iGraphSensitivity = fData[i]->getCombinedSensitivityGraph( true, "" );
        }
        else if( fData[i]->fGraphSensitivity.size() > 0 )
        {
            iGraphSensitivity = fData[i]->fGraphSensitivity[0];
        }
        fData[i]->fGraphSensitivityInterPolated = iGraphSensitivity;
        if( iGraphSensitivity )
        {
            VSensitivityCalculator* b = new VSensitivityCalculator();
            b->setMonteCarloParametersCTA_MC( "", fData[i]->fCameraOffset_deg[0], fCrabSpectraFile, fCrabSpectraID );
            b->setSensitivityGraph( iGraphSensitivity );
            b->setPlotCanvasSize( 900, 600 );
            b->setPlotCanvasName( "cSensitivityInterpolated", "sensitivity vs energy (interpolated)" );
            b->setFluxRange_ENERG( iMinSensitivity, iMaxSensitivity );
            b->setEnergyRange_Lin( fMinEnergy_TeV, fMaxEnergy_TeV );
            b->setPlotCrabLines( bPlotCrabLines );
            b->setMaximumEnergyBias( fMaximumAllowedEnergyBias );
            if( fData[i]->fPlottingColor.size() > 0 )
            {
                b->setPlottingStyle( fData[i]->fPlottingColor[0], fData[i]->fPlottingLineStyle[0], 2.,
                                     fData[i]->fPlottingMarkerStyle[0], 1., fData[i]->fPlottingFillStyle[0] );
            }
            cSensInter = b->plotSensitivityvsEnergyFromCrabSpectrum( cSensInter, iUnit, 0.2, ( i == 0 ) );
            if( fRequirementsString.size() > 0 )
            {
                for( unsigned int p = 0; p < fPlotCTARequirements.size(); p++ )
                {
                    if( fPlotCTARequirements[p] )
                    {
                        fPlotCTARequirements[p]->plotRequirement_DifferentialSensitivity( cSensInter );
                    }
                }
            }
            ///// TMP TMP plot other instruments
            plotCurrentInstruments( cSensInter );
            ///// (END) TMP TMP plot other instruments
            ///// TMP TMP systematics
            // plot systematic line around sensitivity curve
            // plotSensitivitySystematicUncertainties( cSensInter, iGraphSensitivity );
            ///// (END) TMP TMP plot systematics
            iGraphSensitivity->Draw( "p" );

            // plot again the sensitivity graph on top of everything

            //
            // (end of general sensitivity plotting)
            // off-axis sensitivities
            if( fUseIntegratedSensitivityForOffAxisPlots )
            {
                cout << "\t using integrated sensitivity for off-axis sensitivity" << endl;
                TGraphAsymmErrors* iGraphIntegratedSensitivity = fData[i]->getCombinedSensitivityGraph( true, "", true );
                fillProjectedSensitivityPlot( i, iGraphIntegratedSensitivity );
            }
            else
            {
                cout << "\t using differential sensitivity for off-axis sensitivity" << endl;
                fillProjectedSensitivityPlot( i, iGraphSensitivity );
            }
        }
    }
    /////////////////////////////
    // print results
    if( cIntSens )
    {
        plotLegend( cIntSens, false );
        if( iPrint.size() > 0 )
        {
            char hname[2000];
            sprintf( hname, "%s-IntegratedSensitivity.pdf", iPrint.c_str() );
            cIntSens->Print( hname );
        }
    }
    if( cSens )
    {
        plotLegend( cSens, false );
        if( iPrint.size() > 0 )
        {
            char hname[2000];
            sprintf( hname, "%s-Sensitivity.pdf", iPrint.c_str() );
            cSens->Print( hname );
        }
    }
    if( cSensInter )
    {
        plotLegend( cSensInter, false );
        if( iPrint.size() > 0 )
        {
            char hname[2000];
            sprintf( hname, "%s-SensitivityInter.pdf", iPrint.c_str() );
            cSensInter->Print( hname );
        }
    }
    if( cBck )
    {
        //plotLegend( cBck, false );
        if( iPrint.size() > 0 )
        {
            char hname[2000];
            sprintf( hname, "%s-BRates.pdf", iPrint.c_str() );
            if( cBck )
            {
                cBck->Print( hname );
            }
        }
    }
    
    return true;
}

/*
 * plot a legend with the names of the different array layouts
 *
 */
bool VPlotWPPhysSensitivity::plotLegend( TCanvas* c, bool iDown, bool iLeft, bool iAddFirst )
{
    if( bPlotNoLegend )
    {
         return true;
    }
    if( !c )
    {
        return false;
    }
    c->cd();
    
    double x = 0.2 + 0.35;
    // TMPTMP
    x = 0.1 + 0.35;
    // (END) TMPTMP
    if( iLeft )
    {
        x = 0.15;
    }
    double y = 0.60;
    if( iDown )
    {
        y -= 0.50;
        if( y <= 0.1 )
        {
            y = 0.13;
        }
    }
    double y_yp = y + 0.18;
    if( fData.size() == 2 )
    {
        y_yp -= 0.1;
    }
    // TMPTMP
    double x_p = 0.3;
    x_p = 0.4;
    y_yp += 0.1;
    // (END) TMPTMP
    TLegend* iL = new TLegend( x, y, x + x_p, y_yp );
    iL->SetBorderSize( 0 );
    //iL->SetFillColorAlpha( 0, 0.1 );
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( i == 0 && !iAddFirst )
        {
            continue;
        }
        for( unsigned int j = 0; j < fData[i]->fSiteFile_exists.size(); j++ )
        {
            if( fData[i]->fSiteFile_exists[j] && fData[i]->fLegend.size() > 0 )
            {
                TGraph* g = new TGraph( 1 );
                g->SetLineColor( fData[i]->fPlottingColor[j] );
                g->SetLineStyle( fData[i]->fPlottingLineStyle[j] );
                g->SetFillStyle( fData[i]->fPlottingFillStyle[j] );
                g->SetFillColor( fData[i]->fPlottingColor[j] );
                g->SetMarkerColor( fData[i]->fPlottingColor[j] );
                g->SetMarkerStyle( fData[i]->fPlottingMarkerStyle[j] );

//                g->SetMarkerStyle( 1 );
                if( fData[i]->fLegend[j].size() > 0 && fData[i]->fLegend[j].find( "NO_LEGEND" ) == string::npos )
                {
                    iL->AddEntry( g, fData[i]->fLegend[j].c_str(), "pl" );
                }
            }
        }
    }
    string iPadName = c->GetName();
    if( fPlotCTARequirements.size() > 0 && fPlotCTARequirements[0] )
    {
        // add requirements line to the tlegend
        TGraph* i_g_req = new TGraph( 1 );
        i_g_req->SetLineColor( 2 );
        i_g_req->SetLineStyle( 2 );
        TGraph* i_g_goal = new TGraph( 1 );
        i_g_goal->SetLineColor( 3 );
        i_g_goal->SetLineStyle( 2 );
        
        // no requirements for effective areas
        if( iPadName.find( "cEA_EFF" ) == string::npos )
        {
            ostringstream i_Title;
            if( fPlotCTARequirements.size() == 1 )
            {
                i_Title << " Requirement (" << fPlotCTARequirements[0]->getTitle() << ")";
            }
            else
            {
                i_Title << " Requirement";
            }
            iL->AddEntry( i_g_req, i_Title.str().c_str(), "pl" );
        }
    }
    
    // draw legend
    if( iL->GetNRows() > 0 )
    {
        iL->Draw();
    }
    return true;
}

/*
 * plot requirements
 * (updated requirements from Dec 2017)
 *
*/
bool VPlotWPPhysSensitivity::setPlotCTARequirements( string iRequirements, 
        float fRequirementsScalingFactor, 
        double iRequirementsLineWidth,
        bool iRequirementsSystematics )
{
    fRequirementsString = iRequirements;

    if( fRequirementsString.size() == 0 )
    {
        return false;
    }
    
    fPlotCTARequirements.push_back( new VCTARequirements() );
    fPlotCTARequirements.back()->setRequirementsGraphLineWidth( iRequirementsLineWidth );
    fPlotCTARequirements.back()->setRequirementsPlotSystematics( iRequirementsSystematics );
    
    return fPlotCTARequirements.back()->setRequirement( iRequirements, fRequirementsScalingFactor );
}

/*
 * set options for North/South Comparision
*/
void VPlotWPPhysSensitivity::setNorthSouthComparision( bool iNS )
{
    fNorthSouthComparision = iNS;
}

/*
 * get systematic error as function of log_10 energy (TeV)
 *
 */

double VPlotWPPhysSensitivity::getSensitivitySystematicUncertaintiesFactor( double iLe )
{
    if( iLe < log10( 0.2 ) )
    {
        return 0.4;
    }

    return 0.3;
}

/*
 * add systematic uncertainties as 'bracket' errors
 *
 * values are hardwired
 *
 */
void VPlotWPPhysSensitivity::plotSensitivitySystematicUncertainties(
        TCanvas* c, TGraphAsymmErrors* iGraph )
{
    if( c )
    {
        c->cd();
    }

    double fsys = 0.3;

    double x = 0.;
    double y = 0.;
    TGraphAsymmErrors *iGraphSys = new TGraphAsymmErrors( 1 );
    iGraphSys->SetFillColor( 16 );
    iGraphSys->SetFillStyle( 3144 );
    TGraphAsymmErrors *iGraphSysBrackets = new TGraphAsymmErrors( 1 );
    iGraphSysBrackets->SetLineColor( iGraph->GetLineColor() );
    int z = 0;
    for( int i = 0; i < iGraph->GetN(); i++ )
    {
        fsys = getSensitivitySystematicUncertaintiesFactor( x );
        iGraph->GetPoint( i, x, y );
        iGraphSysBrackets->SetPoint( z, x, y );
        iGraphSysBrackets->SetPointEYhigh( z, y * fsys );
        iGraphSysBrackets->SetPointEYlow( z, y * fsys );

        // (below temporary)
        if( i == 0 )
        {
            iGraph->GetPoint( i, x, y );
            y = iGraph->Eval( x - iGraph->GetErrorXlow( i ) );
            iGraphSys->SetPoint( z, x - iGraph->GetErrorXlow( i ), y );
            iGraphSys->SetPointEYhigh( z, y * fsys );
            iGraphSys->SetPointEYlow( z, y * fsys  );
            z++;
        }
        iGraph->GetPoint( i, x, y );
        iGraphSys->SetPoint( z, x, y );
        iGraphSys->SetPointEYhigh( z, y * fsys );
        iGraphSys->SetPointEYlow( z, y * fsys  );
        z++;
        y = iGraph->Eval( x + iGraph->GetErrorXlow( i ) );
        iGraphSys->SetPoint( z, x + iGraph->GetErrorXhigh( i ), y );
        iGraphSys->SetPointEYhigh( z, y * fsys );
        iGraphSys->SetPointEYlow( z, y * fsys  );
        z++;

    }  
    //iGraphSys->Draw( "3" );

    iGraphSysBrackets->Draw( "[]" );
    

}

/*
 * default list of instruments to be plotted
 */
void VPlotWPPhysSensitivity::setCurrentInstrumentPlotVector()
{
    fCurrentInstrumentVector.clear();
    fCurrentInstrumentVector.push_back( "LATp8_0_0" );
    fCurrentInstrumentVector.push_back( "LATp8_120_145" );
    fCurrentInstrumentVector.push_back( "MAGIC" );
    fCurrentInstrumentVector.push_back( "VERITAS" );
    fCurrentInstrumentVector.push_back( "HESS" );
    fCurrentInstrumentVector.push_back( "HAWC5y" );
    fCurrentInstrumentVector.push_back( "HAWC1y" );
    fCurrentInstrumentVector.push_back( "LHAASO1y" );
    fCurrentInstrumentVector.push_back( "SWGO1y" );
    fCurrentInstrumentVector.push_back( "SWGO5y" );
}

/*
 * temporary routing
 *
 * plot sensitivities from current instruments
 */

vector< TGraph* > VPlotWPPhysSensitivity::plotCurrentInstruments( TCanvas* c )
{
    if( c )
    {
        c->cd();
    }
    
    vector< TGraph* > iG;
    
    TFile* iF = new TFile( fCurrentInstrumentRootFile.c_str() );
    if( iF->IsZombie() )
    {
        return iG;
    }
    cout << "reading current instrument performances from : " << fCurrentInstrumentRootFile << endl;
    
    for( unsigned int i = 0; i < fCurrentInstrumentVector.size(); i++ )
    {
        iG.push_back( ( TGraph* )iF->Get( fCurrentInstrumentVector[i].c_str() ) );
        if( iG.back() )
        {
            iG.back()->SetLineWidth( 2 );
            // scale Fermi sensitivity by 25% (4 vs 5 bins per decade)
            // (2017-09-29 GM: already done for the current root file)
            /*            if( fCurrentInstrumentVector[i] == "LATp8_0_0" || fCurrentInstrumentVector[i] == "LATp8_120_145" )
                        {
                            double x = 0.;
                            double y = 0.;
                            for( int i = 0; i < iG.back()->GetN(); i++ )
                            {
                                iG.back()->GetPoint( i, x, y );
                                iG.back()->SetPoint( i, x, y*1.25 );
                            }
                        } */
            if( c )
            {
                iG.back()->Draw( "l" );
            }
        }
/*        if( c )
        {
            string i_text = fCurrentInstrumentVector[i] + "_text";
            TText* i_t = ( TText* )iF->Get( i_text.c_str() );
            // TMPTMP: plot only KSP...
            if( i_t && fPlotCTARequirementsID == 3 )
            {
                i_t->Draw();
           } 
        } */
    }
    return iG;
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void VPPUTValues::add( string iName, TGraph *iG )
{

    // maximum of energy range is site dependent
    // (hard coded):
    //    - South: 100 TeV
    //    - North: 50 TeV
    double iMaxE = 110.;
    // iMaxE = 50.;

    fSetName.push_back( iName );
    fPPUT.push_back( getPPUT( iG, false, log10( 0.02 ), log10( iMaxE ) ) );
    fPPUTError.push_back( getPPUT( iG, true, log10( 0.02 ), log10( iMaxE ) ) );

    // low energy PPUT from 20 GeV to 150 GeV 
    // (energy range for which LSTs much fulfill full system requirements)
    fLowEPPUT.push_back( getPPUT( iG, false, log10( 0.02 ), log10( 0.11 ) ) );
    fLowEPPUTError.push_back( getPPUT( iG, true, log10( 0.02 ), log10( 0.11 ) ) );

    // medium energy PPUT 
    fMidEPPUT.push_back( getPPUT( iG, false, log10( 0.09 ), log10( 11. ) ) );
    fMidEPPUTError.push_back( getPPUT( iG, true, log10( 0.09 ), log10( 11. ) ) );

    // high energy PPUT 
    fHighEPPUT.push_back( getPPUT( iG, false, log10( 5. ), log10( iMaxE ) ) );
    fHighEPPUTError.push_back( getPPUT( iG, true, log10( 5. ), log10( iMaxE ) ) );

    // no LST PPUT
    fnoLoEPPUT.push_back( getPPUT( iG, false, log10( 0.09 ), log10( iMaxE ) ) );
    fnoLoEPPUTError.push_back( getPPUT( iG, true, log10( 0.09 ), log10( iMaxE ) ) );
}

double VPPUTValues::getPPUT( TGraph *iG, bool iError, double ilogEMin, double ilogEMax )
{
    if( iG )
    {
        double x = 0.;
        double y = 0.;
        double dy = 0.;
        double iFOM = 1.;
        double iFOMerror = 0.;
        double n = 0.;
        for( int i = 0; i < iG->GetN(); i++ )
        {
            iG->GetPoint( i, x, y );
            dy = iG->GetErrorYlow( i );
            // if( y > 1.e-3 && x > ilogEMin && x < ilogEMax )
            // allow PPUT to be zero if there is no sensitivity
            // in the required energy range
            if( (x - ilogEMin ) > 1.e-3 && ( ilogEMax - x ) > 1.e-3 )
            {
                iFOM *= y;
                iFOMerror += dy;
                n++;
            }
            // check if first point is close to the required minimal energy
            if( i == 0 )
            {
                if( x - ilogEMin > 0.1 )
                {
                    return 0.;
                }
                // HARD FIX to remove 20 GeV weirdness
                if( ilogEMin < -1.65 && y < 0.15 )
                {
                    return 0.;
                }
            }
        }
        if( n > 0. )
        {
            if( iError )
            {
                iFOMerror = TMath::Power( iFOM, 1./n ) * sqrt( iFOMerror );
                iFOMerror = 1./n * TMath::Power( TMath::Power( iFOM, 1./n ), 1./n-1) * iFOMerror;
                return iFOMerror;
            }
            else
            {
                return TMath::Power( iFOM, 1./n );
            }
        }

     }
    return 0.;
}

void VPPUTValues::print()
{
    for( unsigned int i = 0; i < fSetName.size(); i++ )
    {
        cout << "PPUT: ";
        cout << fSetName[i] << ": ";
        cout << fPPUT[i] << " +- " << fPPUTError[i] << endl;
    }
}

void VPPUTValues::printLatexTable()
{

    for( unsigned int i = 0; i < fSetName.size(); i++ )
    {
        cout << "\\hyperref[SensMult-Array-" + VUtilities::removeSpaces(fSetName[i]);
        cout << "]{" + VUtilities::removeSpaces(fSetName[i]) << "} & ";
        cout << fixed << setfill('0') << setw(4);
        if( fPPUT[i] > 1.e-5 )
        {
            cout << setprecision( 2 ) << fPPUT[i] << " $\\pm$ " << setprecision( 2 ) << fPPUTError[i] << " & ";
        }
        else
        {
            cout << " - & ";
        }
        cout << setprecision( 2 ) << fnoLoEPPUT[i] << " $\\pm$ " << setprecision( 2 ) << fnoLoEPPUTError[i] << " & ";
        if( fPPUT[i] > 1.e-5 )
        {
            cout << setprecision( 2 ) << fLowEPPUT[i] << " $\\pm$ " << setprecision( 2 ) << fLowEPPUTError[i] << " & ";
        }
        else
        {
            cout << " - & ";
        }
        cout << setprecision( 2 ) << fMidEPPUT[i] << " $\\pm$ " << setprecision( 2 ) << fMidEPPUTError[i] << " & ";
        cout << setprecision( 2 ) << fHighEPPUT[i] << " $\\pm$ " << setprecision( 2 ) << fHighEPPUTError[i]; 
        cout << "\\\\" << endl;
    }
}
