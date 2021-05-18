/*!  VWPPhysSensitivityPlotsMaker

     plot sets of sensitivities for CTA/VTS

     Example (executed on the root commane line):

     WPPhysSensitivityPlotsMaker a;
     a.compareDataSets("DataSets.list");

     with < DataSets.list > see example files

*/

#include "VWPPhysSensitivityPlotsMaker.h"

VWPPhysSensitivityPlotsMaker::VWPPhysSensitivityPlotsMaker()
{

    vector< double > iOffAxisValue;
    iOffAxisValue.push_back( 0.5 );
    iOffAxisValue.push_back( 1.5 );
    iOffAxisValue.push_back( 2.5 );
    iOffAxisValue.push_back( 3.25 );
    iOffAxisValue.push_back( 3.75 );
    //    iOffAxisValue.push_back( 4.25 );
    //    iOffAxisValue.push_back( 4.75 );
    
    cout << "VWPPhysSensitivityPlotsMaker: hardwired offsets from camera center: ";
    for( unsigned int i = 0; i < iOffAxisValue.size(); i++ )
    {
        cout << iOffAxisValue[i] << ", ";
    }
    cout << " [deg]" << endl;
    
    setOffAxisAngle( iOffAxisValue );
    setEnergyRange_Lin_TeV();
    setObservingTime();
    setAxisUnits();
    setSensitivityRatioLimits();
    setIRFRatioLimits();
    setEffectiveAreaLimits();
    setResolutionLimits();
    setPrintingOptions();
    setPlotRequirements();
    setPlotCrabLines();
    setPlotRequirementsSystematics();
    setMaximumAllowedEnergyBias();
    setPlotAllInOneCanvasStyle();
    
    fPlotAllInOneCanvas = 0;
    fPlotProjectedSensitivity = 0;
    fSensitivityPad = 0;
    fSensitivityTitlePad = 0;
    fSensitivityRatioPad = 0;
    fEffAreaPad = 0;
    fBckRatesPad = 0;
    fERes = 0;
    fAngRes = 0;
    fPlotCurrentInstrumentVectorLabel = false;

    setPlotNoLegend();
    
}

/*
 * set requirements
 *
 * e.g. South-50h, North-50h, South-30m-LST, ...
*/
void VWPPhysSensitivityPlotsMaker::setPlotRequirements( string iRequirement, float iRequirementScaling, double iRequirementsLineWidth )
{
    fPlotCTARequirementsString = iRequirement;
    fRequirementsScalingFactor = iRequirementScaling;
    fRequirementsLineWidth = iRequirementsLineWidth;
}

void VWPPhysSensitivityPlotsMaker::plotAllInOneCanvas( bool iCanvasBatch )
{
    fPlotAllInOneCanvas = new TCanvas( "cSensitivityFullCanvas", "sensitivity", 1, 1, 1400, 845 );
    fPlotAllInOneCanvas->SetTopMargin( 0.03 );
    fPlotAllInOneCanvas->SetGridx( 0 );
    fPlotAllInOneCanvas->SetGridy( 0 );
    if( iCanvasBatch )
    {
        fPlotAllInOneCanvas->SetBatch();
    }
    
    // Canvas for eff ratio, background ratio, angular resolution ratio
    if( fPlotAllInOneCanvasStyle == "AllRatios" )
    {
        fSensitivityTitlePad = new TPad( "cSensitivityTitlePad", "title", 0.01, 0.99, 0.68, 1.00 );
        fSensitivityPad = new TPad( "cSensitivityPad", "sensitivity", 0.01, 0.34, 0.60, 0.99 );
        fSensitivityRatioPad = new TPad( "cSensitivityPadRatio", "sensitivity ratio", 0.02, 0.01, 0.40, 0.34 );
        fEffAreaPad = new TPad( "cEffAreaPad", "effective area", 0.60, 0.67, 0.80, 0.99 );
        fEffAreaRatioPad = new TPad( "cEffAreaRatioPad", "effective area", 0.80, 0.67, 0.99, 0.99 );
        fBckRatesPad = new TPad( "cBckRatesPad", "background rates", 0.60, 0.34, 0.80, 0.67 );
        fBckRatesPadRatio = new TPad( "cBckRatesPadRatio", "background rates", 0.80, 0.34, 0.99, 0.67 );
        fAngRes = new TPad( "cAngRes", "angular resolution", 0.60, 0.01, 0.80, 0.34 );
        fAngResRatio = new TPad( "cAngResRatio", "angular resolution", 0.80, 0.01, 0.99, 0.34 );
        fERes = new TPad( "cERes", "energy resolution", 0.40, 0.01, 0.60, 0.34 );
    }
    // Six canvas plot (old style, as used before Feb 2021)
    else if( fPlotAllInOneCanvasStyle == "SensitivityRatioOnly" )
    {
        fSensitivityTitlePad = new TPad( "cSensitivityTitlePad", "title",  0.01, 0.95, 0.68, 0.99 );
        fSensitivityPad = new TPad( "cSensitivityPad", "sensitivity", 0.01, 0.34, 0.68, 0.95 );
        fSensitivityRatioPad = new TPad( "cSensitivityPadRatio", "sensitivity ratio", 0.02, 0.01, 0.38, 0.34 );
        fEffAreaPad = new TPad( "cEffAreaPad", "effective area", 0.69, 0.67, 0.99, 0.99 );
        fEffAreaRatioPad = 0;
        fBckRatesPad = new TPad( "cBckRatesPad", "background rates", 0.69, 0.34, 0.99, 0.67 );
        fBckRatesPadRatio = 0;
        fAngRes = new TPad( "cAngRes", "angular resolution", 0.69, 0.01, 0.99, 0.34 );
        fAngResRatio = 0;
        fERes = new TPad( "cERes", "energy resolution", 0.40, 0.01, 0.68, 0.34 );
    }


    if( fSensitivityTitlePad )
    {
        fSensitivityTitlePad->SetLeftMargin( 0.15 );
        fSensitivityTitlePad->SetRightMargin( 0.03 );
        fSensitivityTitlePad->SetTopMargin( 0.03 );
        fSensitivityTitlePad->SetBottomMargin( 0.03 );
        fSensitivityTitlePad->Draw();
    }
    if( fSensitivityPad )
    {
        fSensitivityPad->SetLeftMargin( 0.15 );
        fSensitivityPad->SetRightMargin( 0.03 );
        fSensitivityPad->Draw();
    }
    if( fEffAreaPad )
    {
        fEffAreaPad->SetLeftMargin( 0.15 );
        fEffAreaPad->SetRightMargin( 0.05 );
        fEffAreaPad->Draw();
    }
    if( fEffAreaRatioPad )
    {
        fEffAreaRatioPad->SetLeftMargin( 0.15 );
        fEffAreaRatioPad->SetRightMargin( 0.05 );
        fEffAreaRatioPad->Draw();
    }
    if( fBckRatesPad )
    {
        fBckRatesPad->SetLeftMargin( 0.15 );
        fBckRatesPad->SetRightMargin( 0.05 );
        fBckRatesPad->Draw();
    }
    if( fBckRatesPadRatio )
    {
        fBckRatesPadRatio->SetLeftMargin( 0.15 );
        fBckRatesPadRatio->SetRightMargin( 0.05 );
        fBckRatesPadRatio->Draw();
    }
    if( fSensitivityRatioPad )
    {
        fSensitivityRatioPad->SetBottomMargin( 0.13 );
        fSensitivityRatioPad->SetTopMargin( 0.05 );
        fSensitivityRatioPad->SetRightMargin( 0.05 );
        fSensitivityRatioPad->Draw();
    }
    if( fAngRes )
    {
        fAngRes->SetLeftMargin( 0.13 );
        fAngRes->SetRightMargin( 0.05 );
        fAngRes->SetBottomMargin( 0.13 );
        fAngRes->SetTopMargin( 0.05 );
        fAngRes->Draw();
    }
    if( fAngResRatio )
    {
        fAngResRatio->SetLeftMargin( 0.13 );
        fAngResRatio->SetRightMargin( 0.05 );
        fAngResRatio->SetBottomMargin( 0.13 );
        fAngResRatio->SetTopMargin( 0.05 );
        fAngResRatio->Draw();
    }
    if( fERes )
    {
        fERes->SetLeftMargin( 0.15 );
        fERes->SetRightMargin( 0.05 );
        fERes->SetBottomMargin( 0.13 );
        fERes->SetTopMargin( 0.05 );
        fERes->Draw();
    }
}

/*
   plot ratios of graphs in a canvas

*/
void VWPPhysSensitivityPlotsMaker::plotRatioPlot( TPad *corg,
                                                  TPad *cplot,
                                                  double ymin,
                                                  double ymax,
	                                          bool revert_ratio )
{
    if( !corg ) return;
    if( !cplot ) return;
    cplot->cd();

    unsigned int z = 0;
    TGraphAsymmErrors *iNullGraph = 0;

    TIter next(corg->GetListOfPrimitives());
    while (TObject *obj = next() )
    {
        string objname = obj->GetName();
        string objclass = obj->ClassName();
	string iYTitle;
        if( objclass == "TH1D" || objclass == "TH1F" )
        {
            TH1F *h = (TH1F*)obj;
            TH1F *hR = new TH1F( (objname+"ratio").c_str(), "",
                                  h->GetXaxis()->GetNbins(),
                                  h->GetXaxis()->GetXmin(), 
                                  h->GetXaxis()->GetXmax() );
            hR->SetXTitle( h->GetXaxis()->GetTitle() );
            iYTitle = h->GetYaxis()->GetTitle();
            iYTitle = iYTitle.substr( 0, iYTitle.find( "[" ) );
            hR->SetYTitle( (iYTitle + "ratio" ).c_str() );
            hR->SetMinimum( ymin );
            hR->SetMaximum( ymax );
            hR->SetStats( 0 );
            hR->Draw( "" );
            hR->Draw( "AH" );
            plot_nullHistogram( cplot, hR, 
                                true, false, 1.2, 
                                TMath::Power( 10., h->GetXaxis()->GetXmin() ),
                                TMath::Power( 10., h->GetXaxis()->GetXmax() ) );
            TLine *iL = new TLine( h->GetXaxis()->GetXmin(), 1., 
                                   h->GetXaxis()->GetXmax(), 1. );
            iL->SetLineStyle( 2 );
            iL->Draw();
        }
        else if( objclass == "TGraphAsymmErrors" ||
                 objclass == "TGraphErrors" )
        {
            TGraphAsymmErrors *r = (TGraphAsymmErrors*)obj;
            if( z == 0 )
            {
                iNullGraph = r;
            }
            else if( r->GetN() > 1 )
            {
                TGraphAsymmErrors *n = new TGraphAsymmErrors( 1 );
                bool divide_success = false;
                if( revert_ratio )
                {
                    divide_success = VHistogramUtilities::divide( n,
                                          r,
                                          iNullGraph,
                                          1.e-3,
                                          true );
                        }
                else
                {
                    divide_success = VHistogramUtilities::divide( n,
                                          iNullGraph,
                                          r,
                                          1.e-3,
                                          true );
                        }
                if( divide_success )
                {
                    n->SetLineColor( r->GetLineColor() );
                    n->SetMarkerColor( r->GetMarkerColor() );
                    n->SetMarkerStyle( r->GetMarkerStyle() );
                    n->SetMarkerSize( r->GetMarkerSize() );
                    n->Draw( "p" );
                }
            }
            z++;
        }
    }
}

/*
 * make IRF comparison plots
 *
 * read IRFs from root files
 *
 */
void VWPPhysSensitivityPlotsMaker::compareDataSets( string iDataSetFile,
        string iDirectionString,
        bool iUseIntegratedSensitivityForOffAxisPlots,
        unsigned int iRatioCounter,
        string iTitleText )
{
    VPlotWPPhysSensitivity a;
    a.setPlotNoLegend( bPlotNoLegend );
    a.setPlotCrabLines( bPlotCrabLines );
    if( fCurrentInstrumentVector.size() > 0 )
    {
        a.setCurrentInstrumentPlotVector( fCurrentInstrumentVector,
		                          fPlotCurrentInstrumentVectorLabel );
    }
    if( iRatioCounter == 9999 )
    {
        a.setNorthSouthComparision( true );
        setSensitivityRatioLimits( 0.01, 1.1 );
    }
    else if( iRatioCounter == 9998 )
    {
        a.setNorthSouthComparision( true );
        setSensitivityRatioLimits( 0. , 3.7 );
    }
    //	a.setUseIntegratedSensitivityForOffAxisPlots( iUseIntegratedSensitivityForOffAxisPlots );
    a.setPlotCTARequirements( fPlotCTARequirementsString, 
            fRequirementsScalingFactor, 
            fRequirementsLineWidth,
            fRequirementsSystematics );
    // (outdated, orginal requirements plotting)
    a.setEnergyRange_Lin_TeV( fMinEnergy_TeV, fMaxEnergy_TeV );
    if( !a.addDataSets( iDataSetFile, iDirectionString ) )
    {
        return;
    }
    a.plotIRF( fPrintingOptions,
               fEffArea_min, fEffArea_max,
               fAngularResolution_min, fAngularResolution_max,
               fEnergyResolution_min, fEnergyResolution_max,
               fEffAreaPad, fAngRes, fERes, false,
               fPlotAngResLogY );
    a.plotSensitivity( fPrintingOptions, fSensitivity_min, fSensitivity_max, fSensitivity_Unit, fSensitivityPad, fBckRatesPad );
    //fPlotProjectedSensitivity = a.plotProjectedSensitivities( 0, 5. );
    if( iRatioCounter == 998 )
    {
        a.plotSensitivityRatio( fPrintingOptions, fSensitivityRatio_min, fSensitivityRatio_max, 0, fSensitivityRatioPad, iRatioCounter );
    }
    else
    {
        a.plotSensitivityRatio( fPrintingOptions, fSensitivityRatio_min, fSensitivityRatio_max, 2, fSensitivityRatioPad, iRatioCounter );
    }
    // background ratio
    if( fBckRatesPad )
    {
        plotRatioPlot( fBckRatesPad, fBckRatesPadRatio, fBackgroundRatio_min, fBackgroundRatio_max );
    }
    if( fAngRes )
    {
        plotRatioPlot( fAngRes, fAngResRatio, fAngResRatio_min, fAngResRatio_max, true );
    }
    if( fEffAreaPad )
    {
        plotRatioPlot( fEffAreaPad, fEffAreaRatioPad, fEffAreaRatio_min, fEffAreaRatio_max );
    }
    
    if( iTitleText.size() > 0 && fSensitivityTitlePad )
    {
        fSensitivityTitlePad->cd();
        
        TText* iTTitleText = new TText( 0.01, 0.2, iTitleText.c_str() );
        //            iTTitleText->SetTextAngle( 90 );
        iTTitleText->SetTextSize( 0.8 );
        iTTitleText->SetTextAlign( 11 );
        iTTitleText->SetTextColor( 1 );
        iTTitleText->SetTextFont( 42 );
        iTTitleText->Draw();
    }
    
}

void VWPPhysSensitivityPlotsMaker::printPlotCTARequirementsIDs()
{
    cout << "requirements IDs: " << endl;
    cout << "0 \t South, 50h, FOM for energy range [0.03,100]" << endl;
    cout << "1 \t South, 5h, FOM for energy range [0.03,100]" << endl;
    cout << "2 \t South, 0.5h, FOM for energy range [0.03,100]" << endl;
    cout << "3 \t North, 50h, FOM for energy range [0.03,20]" << endl;
    cout << "4 \t North, 5h, FOM for energy range [0.03,20]" << endl;
    cout << "5 \t North, 0.5h, FOM for energy range [0.03,20]" << endl;
    cout << endl;
}

void VWPPhysSensitivityPlotsMaker::compareOffAxisSensitivities( string iSubArray, string iDataSet )
{
    vector< string > iD;
    iD.push_back( iDataSet );
    compareOffAxisSensitivities( iSubArray, iD );
}

/*

   compare off axis sensitivities

   note that some values are hardwired in plotProjectedSensitivities

*/
void VWPPhysSensitivityPlotsMaker::compareOffAxisSensitivities( string iSubArray, vector< string > iDataSet )
{
    if( iSubArray.size() > 0 )
    {
        fListOfArrays.clear();
        fListOfArrays.push_back( iSubArray );
    }
    if( iDataSet.size() > 0 )
    {
        fListofDataSets.clear();
        fListofDataSets = iDataSet;
    }
    cout << "Compare " << fListOfArrays.size() << " array(s) in " << fListofDataSets.size() << " data set(s)" << endl;
    
    //TCanvas* c = 0;
    for( unsigned int j = 0; j < fListofDataSets.size(); j++ )
    {
        for( unsigned int i = 0; i < fListOfArrays.size(); i++ )
        {
            VPlotWPPhysSensitivity a;
            a.setEnergyRange_Lin_TeV( fMinEnergy_TeV, fMaxEnergy_TeV );
            for( unsigned int k = 0; k < fOffAxisAngle.size(); k++ )
            {
                cout << "\t adding " << fListofDataSets[j] << ", " << fListOfArrays[i];
                cout << ", " << fObservingTime_s << ", " << fOffAxisAngle[k] << ", " << k << ", " << j << endl;
                a.addDataSet( fListofDataSets[j], fListOfArrays[i], fObservingTime_s, fOffAxisAngle[k], "", k + 1, j + 1 );
            }
            string iP = "";
            if( fPrintingOptions.size() > 0 )
            {
                iP += fPrintingOptions + "-" + fListOfArrays[i];
            }
            a.plotIRF( iP );
            a.plotSensitivity( iP, fSensitivity_min, fSensitivity_max, fSensitivity_Unit );
            
            // c = a.plotProjectedSensitivities( c, fOffAxisAngle.back(), j + 1 );
        }
    }
}

/*
 * set the correct flux axis units for the given observation time
 *
 * known values: 50h, 5h, 0.5h
 */
void VWPPhysSensitivityPlotsMaker::setAxisUnits( string iObservationTime )
{
    if( iObservationTime == "0.5h" || iObservationTime == "30m" )
    {
        cout << "Setting axis units for 0.5h observations" << endl;
        setAxisUnits( 1.e-12, 5.5e-9 );
    }
    else if( iObservationTime == "5h" || iObservationTime == "05h" )
    {
        cout << "Setting axis units for 5h observations" << endl;
        setAxisUnits( 2.e-13, 2.5e-10 );
    }
    else if( iObservationTime == "100s" )
    {
        cout << "Setting axis units for 100s observations" << endl;
        setAxisUnits( 1.e-11, 2.5e-9 );
    }
    else
    {
        cout << "Setting axis units for 50h observations" << endl;
        setAxisUnits();
    }
}

void VWPPhysSensitivityPlotsMaker::setAxisUnits( double iMinSensitivity, double iMaxSensitivity, string iUnit )
{
    fSensitivity_min = iMinSensitivity;
    fSensitivity_max = iMaxSensitivity;
    fSensitivity_Unit = iUnit;
}

void VWPPhysSensitivityPlotsMaker::resetVectors()
{
    fListOfArrays.clear();
    fListofDataSets.clear();
    fOffAxisAngle.clear();
}


bool VWPPhysSensitivityPlotsMaker::writeTexFileBody( string iTexFile, string iTexFileTitle )
{
    // tex file
    cout << "Writing tex file: " << iTexFile << endl;
    ofstream os;
    os.open( iTexFile.c_str() );
    if( !os )
    {
        cout << "VWPPhysSensitivityPlotsMaker::writeTexFileBody: failed writing to " << iTexFile << endl;
        return false;
    }
    
    // intro
    os << "\\documentclass[11pt]{scrartcl}" << endl;
    os << "\\usepackage[a4paper,landscape,scale=0.9]{geometry}" << endl;
    os << "\\usepackage{graphicx}" << endl;
    os << "\\usepackage{epstopdf}" << endl;
    os << "\\usepackage[pdftex,colorlinks=true,bookmarks=false,bookmarksopen=false]{hyperref}" << endl;
    
    os << "\\title{CTA sensitivities with EVNDISP \\\\ ";
    os << iTexFileTitle << "}" << endl;
    os << "\\author{Gernot Maier \\\\ DESY}" << endl;
    os << "\\date{\\today}" << endl;
    
    os << "\\begin{document}" << endl;
    os << "\\maketitle" << endl;
    os << "\\tableofcontents" << endl;
    
    os << "\\newpage" << endl;
    
    // images
    if( fListOfArrays.size() > 0 )
    {
        for( unsigned int i = 0; i < fListOfArrays.size(); i++ )
        {
            os << "\\begin{figure}" << endl;
            os << "\\centering\\includegraphics[width=0.69\\linewidth]{" << fPrintingOptions << "-" << fListOfArrays[i] << "-Sensitivity.pdf}" << endl;
            os << "\\centering\\includegraphics[width=0.29\\linewidth]{ArrayLayout-" << fListOfArrays[i] << ".pdf}" << endl;
            os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-" << fListOfArrays[i] << "-EffArea.pdf}" << endl;
            os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-" << fListOfArrays[i] << "-BRates.pdf}" << endl;
            os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-" << fListOfArrays[i] << "-AngRes.pdf}" << endl;
            os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-" << fListOfArrays[i] << "-ERes.pdf}" << endl;
            os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-" << fListOfArrays[i] << "-EBias.pdf}" << endl;
            os << "\\end{figure}" << endl;
            os << "\\clearpage" << endl;
            os << endl << endl;
        }
    }
    else
    {
        os << "\\begin{figure}" << endl;
        os << "\\centering\\includegraphics[width=0.69\\linewidth]{" << fPrintingOptions << "-Sensitivity.pdf}" << endl;
        os << "\\centering\\includegraphics[width=0.29\\linewidth]{ArrayLayout-" << fPrintingOptions << ".pdf}" << endl;
        os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-EffArea.pdf}" << endl;
        os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-BRates.pdf}" << endl;
        os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-AngRes.pdf}" << endl;
        os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-ERes.pdf}" << endl;
        os << "\\centering\\includegraphics[width=0.18\\linewidth]{" << fPrintingOptions << "-EBias.pdf}" << endl;
        os << "\\end{figure}" << endl;
        os << "\\clearpage" << endl;
        os << endl << endl;
    }
    // tex file closing
    os << endl;
    os << "\\end{document}" << endl;
    
    os.close();
    
    return true;
}

void VWPPhysSensitivityPlotsMaker::setLSTSettings()
{
    setEnergyRange_Lin_TeV( 0.01, 0.9 );
    setResolutionLimits( 0.5, 0.4 );
    setSensitivityRatioLimits( 0.4, 1.2 );
    setEffectiveAreaLimits( 5000., 4.e6 );
}
