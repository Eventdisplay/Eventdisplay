/*! \file  VPlotInstrumentResponseFunction.cpp
    \brief effective area and IRF plotter

    to get last resolution graph: getLastPlottedGraph()

*/

#include "VPlotInstrumentResponseFunction.h"

VPlotInstrumentResponseFunction::VPlotInstrumentResponseFunction()
{
    fDebug = false;
    
    fName = "EA";
    
    fTF1_fitResolution = 0;
    gLastPlottedGraph = 0;
    
    
    setLegendParameters();
    setResolutionFitting();
    setCanvasSize();
    
    setPlottingDefaults();
    
}

VPlotInstrumentResponseFunction::~VPlotInstrumentResponseFunction()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] )
        {
            delete fData[i];
        }
    }
}

void VPlotInstrumentResponseFunction::setPlottingDefaults()
{
    setPlottingAxis( "energy_Lin", "X", false, 0.005, 200., "energy [TeV]" );
    setPlottingAxis( "distance_Lin", "X", false, 0., 500., "distance [m]" );
    setPlottingAxis( "nimages_Lin", "X", false, 0., 50., "number of images" );
    
    setPlottingAxis( "effarea_Lin", "Y", true, 1.0, 5.e7, "effective area [m^{2}]" );
    setPlottingAxis( "angularesolution_Lin", "Y", false, 0.001, 0.35, "angular resolution [deg]" );
    setPlottingAxis( "angularesolution_Log", "Y", true, 0.001, 1.1, "angular resolution [deg]" );
    setPlottingAxis( "coreresolution_Lin", "Y", false, 0., 40.0, "core resolution [m]" );
    setPlottingAxis( "energyresolution_Lin", "Y", false, 0., 0.40, "energy resolution" );
    
}

/*

   read instrument response list of files and load them for e.g. plotting

   this replaces several VPlotInstrumentResponseFunction::addInstrumentResponseData( filename, ...)

   the format of the file list is:

    * 1 myEffectiveArea1.root 20. 16. 0.5 2.4 200 A_MC p 1 1 20 Number 1
    * 1 myEffectiveArea2.root  20. 16. 0.5 2.4 200 A_MC p 2 1 21 Number 2
    * ID FILENAME (zenith) (azimuth bin) (wobble offset) (spectral index) (noise level) (A_MC/A_REC) (plotting style) (color) (line style) (marker style) (name)
*/
bool VPlotInstrumentResponseFunction::addInstrumentResponseData( int iDataID, string iFileList )
{
    ifstream is;
    is.open( iFileList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VPlotInstrumentResponseFunction::addInstrumentResponseData: error opening input file list: " << iFileList << endl;
        return false;
    }
    string is_line;
    string temp;
    while( getline( is, is_line ) )
    {
        istringstream is_stream( is_line );
        is_stream >> temp;
        if( temp != "*" )
        {
            continue;
        }
        is_stream >> temp;
        // check set number
        if( atoi( temp.c_str() ) != iDataID )
        {
            continue;
        }
        // read this line
        fData.push_back( new VInstrumentResponseFunctionReader() );
        fData.back()->setDebug( fDebug );
        fData.back()->setPlottingStyle( fData.size() + 1, 1, 2., 20 + fData.size() );
        fData.back()->fillData( is_line, iDataID );
    }
    is.close();
    
    listDataSets();
    
    return true;
}

bool VPlotInstrumentResponseFunction::addInstrumentResponseData( string iFile, string iA_MC )
{
    return addInstrumentResponseData( iFile, 20., 0.0, 0, 0, 200, iA_MC );
}

bool VPlotInstrumentResponseFunction::addInstrumentResponseData( string iFile, double iZe, double iWoff,
        int iAzBin, double iIndex, int iNoise, string iA_MC,
        int iColor, int iLineStyle,
        int iMarkerStyle, float iMarkerSize,
        float iEmin_linTeV, float iEmax_linTeV,
        float iLineWidth, string iPlotOption,
        float iEnergy_TeV_Eres_oneSided, string iLegend )
{
    // read IRF
    VInstrumentResponseFunctionReader* iTempIRFReader = new VInstrumentResponseFunctionReader();
    iTempIRFReader->setDebug( fDebug );
    if( iColor < 0 )
    {
        iColor = fData.size() + 1;
    }
    if( iLineStyle < 0 )
    {
        iLineStyle = 1;
    }
    if( iMarkerStyle < 0 )
    {
        iMarkerStyle =  20 + fData.size();
    }
    if( iMarkerSize < 0 )
    {
        iMarkerSize  = 2.;
    }
    iTempIRFReader->setPlottingStyle( iColor, iLineStyle, iLineWidth, iMarkerStyle, iMarkerSize );
    if( iPlotOption.size() > 0 )
    {
        iTempIRFReader->setPlotOption( iPlotOption );
    }
    iTempIRFReader->setEnergyRange( iEmin_linTeV, iEmax_linTeV );
    iTempIRFReader->setEnergyResolutionMethod( 1, false, iEnergy_TeV_Eres_oneSided );
    if( !iTempIRFReader->fillData( iFile, iZe, iWoff, iAzBin, iIndex, iNoise, iA_MC ) )
    {
        cout << "VPlotInstrumentResponseFunction::addInstrumentResponseData() error filling effective area / IRF data" << endl;
        return false;
    }
    iTempIRFReader->fLegend = iLegend;
    fData.push_back( iTempIRFReader );
    
    listDataSets();
    
    return true;
}

bool VPlotInstrumentResponseFunction::removeInstrumentResponseData( int iDataSetID )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return false;
    }
    
    fData.erase( fData.begin() + iDataSetID );
    
    listDataSets();
    
    return true;
}

void VPlotInstrumentResponseFunction::resetInstrumentResponseData()
{
    fData.clear();   // this is not a clean way to get rid of the data -> fix
}

TCanvas* VPlotInstrumentResponseFunction::plotEffectiveArea( double iEffAreaMin_m2, double iEffAreaMax_m2, TPad* iEffAreaPad, bool iSmooth )
{
    if( fData.size() == 0 )
    {
        return 0;
    }
    
    // set min/maximum value in effective area axis
    if( iEffAreaMin_m2 > 0. )
    {
        getPlottingAxis( "effarea_Lin" )->fMinValue = iEffAreaMin_m2;
    }
    if( iEffAreaMax_m2 > 0. )
    {
        getPlottingAxis( "effarea_Lin" )->fMaxValue = iEffAreaMax_m2;
    }
    
    char hname[200];
    
    TCanvas* iEffectiveAreaPlottingCanvas = 0;
    if( iEffAreaPad )
    {
        iEffectiveAreaPlottingCanvas = ( TCanvas* )iEffAreaPad;
    }
    else
    {
        sprintf( hname, "cEA_EFF" );
        iEffectiveAreaPlottingCanvas = new TCanvas( hname, "effective area", 10, 10, fCanvasSize_X, fCanvasSize_Y );
        iEffectiveAreaPlottingCanvas->SetGridx( 0 );
        iEffectiveAreaPlottingCanvas->SetGridy( 0 );
        iEffectiveAreaPlottingCanvas->SetLeftMargin( 0.15 );
        iEffectiveAreaPlottingCanvas->SetRightMargin( 0.07 );
    }
    if( !iEffectiveAreaPlottingCanvas )
    {
        return 0;
    }
    iEffectiveAreaPlottingCanvas->cd();
    TLegend* iLegend = makeLegend();
    
    if( !iEffectiveAreaPlottingCanvas->GetListOfPrimitives()->FindObject( "heff" ) )
    {
        TH1D* heff = new TH1D( "heff", "", 100, log10( getPlottingAxis( "energy_Lin" )->fMinValue ), log10( getPlottingAxis( "energy_Lin" )->fMaxValue ) );
        heff->SetStats( 0 );
        heff->SetXTitle( "log_{10} energy [TeV]" );
        heff->SetYTitle( "effective area [m^{2}]" );
        heff->SetMinimum( getPlottingAxis( "effarea_Lin" )->fMinValue );
        heff->SetMaximum( getPlottingAxis( "effarea_Lin" )->fMaxValue );
        heff->Draw( "" );
        heff->Draw( "AH" );
        
        plot_nullHistogram( iEffectiveAreaPlottingCanvas, heff, getPlottingAxis( "energy_Lin" )->fLogAxis,
                            getPlottingAxis( "effarea_Lin" )->fLogAxis, 1.3,
                            getPlottingAxis( "energy_Lin" )->fMinValue, getPlottingAxis( "energy_Lin" )->fMaxValue );
    }
    
    int z = 0;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        TGraphAsymmErrors* g = 0;
        
        if( fData[i]->fA_MC == "A_MC" )
        {
            g = fData[i]->gEffArea_MC;
        }
        else if( fData[i]->fA_MC == "A_MCnoTH2" || fData[i]->fA_MC == "A_MC_noTH2" )
        {
            g = fData[i]->gEffArea_MCnoTh2;
        }
        else if( fData[i]->fA_MC == "A_RECnoTH2" || fData[i]->fA_MC == "A_REC_noTH2" )
        {
            g = fData[i]->gEffArea_RecnoTh2;
        }
        else
        {
            g = fData[i]->gEffArea_Rec;
        }
        
        if( !g )
        {
            cout << "VPlotInstrumentResponseFunction::plotEffectiveArea() warning: no graph (" << fData[i]->fA_MC << ")";
            cout << " found for data set " << i << endl;
            continue;
        }
        
        if( g->GetN() > 0. )
        {
            if( fDebug )
            {
                g->Print();
            }
            if( iSmooth )
            {
                // smooth this graph using a Kernel of width 0.4
                TGraphSmooth* Smooth = new TGraphSmooth( "s" );
                TGraph* Smooth_graph = Smooth->SmoothKern( g, "normal", 0.4, 200 );
                if( Smooth_graph )
                {
                    Smooth_graph->SetLineStyle( g->GetLineStyle() );
                    Smooth_graph->SetLineColor( g->GetLineColor() );
                    Smooth_graph->SetLineWidth( g->GetLineWidth() );
                    Smooth_graph->SetMarkerStyle( g->GetMarkerStyle() );
                    Smooth_graph->SetMarkerColor( g->GetMarkerColor() );
                    Smooth_graph->SetMarkerSize( g->GetMarkerSize() );
                    Smooth_graph->Draw( fData[i]->fPlotOption.c_str() );
                    if( iLegend )
                    {
                        iLegend->AddEntry( Smooth_graph, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                    }
                }
            }
            else
            {
                g->Draw( fData[i]->fPlotOption.c_str() );
                if( iLegend )
                {
                    iLegend->AddEntry( g, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
            }
            z++;
        }
    }
    if( z > 0 && getPlottingAxis( "effarea_Lin" )->fLogAxis )
    {
        iEffectiveAreaPlottingCanvas->SetLogy( 1 );
    }
    if( iLegend )
    {
        iLegend->Draw();
    }
    return iEffectiveAreaPlottingCanvas;
}

TCanvas* VPlotInstrumentResponseFunction::plotWeightedRate()
{
    if( fData.size() == 0 )
    {
        return 0;
    }
    
    char hname[200];
    
    sprintf( hname, "cEA_WR" );
    TCanvas* iWeightedRatePlottingCanvas = new TCanvas( hname, "rate", 10, 10, fCanvasSize_X, fCanvasSize_Y );
    iWeightedRatePlottingCanvas->SetGridx( 0 );
    iWeightedRatePlottingCanvas->SetGridy( 0 );
    iWeightedRatePlottingCanvas->SetLeftMargin( 0.15 );
    iWeightedRatePlottingCanvas->SetRightMargin( 0.07 );
    iWeightedRatePlottingCanvas->cd();
    TLegend* iLegend = makeLegend();
    
    TH1D* hWeightedRate = new TH1D( "hWeightedRate", "", 100, log10( getPlottingAxis( "energy_Lin" )->fMinValue ),
                                    log10( getPlottingAxis( "energy_Lin" )->fMaxValue ) );
    hWeightedRate->SetStats( 0 );
    hWeightedRate->SetXTitle( "log_{10} energy [TeV]" );
    hWeightedRate->SetYTitle( "rate [1/min]" );
    if( fData.size() > 0 && fData[0]->hWeightedRate )
    {
        hWeightedRate->SetMaximum( fData[0]->hWeightedRate->GetMaximum() * 1.5 );
    }
    hWeightedRate->Draw( "" );
    hWeightedRate->Draw( "AH" );
    
    plot_nullHistogram( iWeightedRatePlottingCanvas, hWeightedRate, getPlottingAxis( "energy_Lin" )->fLogAxis,
                        getPlottingAxis( "effarea_Lin" )->fLogAxis, 1.3,
                        getPlottingAxis( "energy_Lin" )->fMinValue, getPlottingAxis( "energy_Lin" )->fMaxValue );
                        
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] && fData[i]->hWeightedRate )
        {
            setHistogramPlottingStyle( fData[i]->hWeightedRate, i + 1, 2. );
            fData[i]->hWeightedRate->Draw( "same" );
            if( iLegend )
            {
                iLegend->AddEntry( fData[i]->hWeightedRate, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
            }
        }
    }
    if( iLegend )
    {
        iLegend->Draw();
    }
    
    // return plotting canvas
    return iWeightedRatePlottingCanvas;
}

bool VPlotInstrumentResponseFunction::checkDataSetID( unsigned int iDataSetID )
{
    if( iDataSetID >= fData.size() )
    {
        cout << "Error: data set ID out of range. Should be <" << fData.size() << endl;
        return false;
    }
    return true;
}

void VPlotInstrumentResponseFunction::plotCutEfficiency( unsigned int iDataSetID )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return;
    }
    
    char hname[200];
    
    sprintf( hname, "cEA_cuteff_%d", iDataSetID );
    TCanvas* iCutEfficencyPlottingCanvas = new TCanvas( hname, "cut efficiency", 10, 10, fCanvasSize_X, fCanvasSize_Y );
    iCutEfficencyPlottingCanvas->SetGridx( 0 );
    iCutEfficencyPlottingCanvas->SetGridy( 0 );
    
    sprintf( hname, "hceff_%d", iDataSetID );
    TH1D* hceff = new TH1D( hname, "", 100, log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ) );
    hceff->SetStats( 0 );
    hceff->SetXTitle( "log_{10} energy [TeV]" );
    hceff->SetYTitle( "cut efficiency" );
    hceff->SetMinimum( 1. );
    if( fData[iDataSetID]->hCutEfficiency.size() > 0 && fData[iDataSetID]->hCutEfficiency[0] )
    {
        hceff->SetMaximum( fData[iDataSetID]->hCutEfficiency[0]->GetMaximum() * 1.5 );
        cout << "Maximum in 0 histogram for cut efficiency: " <<  fData[iDataSetID]->hCutEfficiency[0]->GetMaximum() << endl;
    }
    hceff->Draw( "" );
    hceff->Draw( "AH" );

    plot_nullHistogram( iCutEfficencyPlottingCanvas, hceff, getPlottingAxis( "energy_Lin" )->fLogAxis, true,
                        hceff->GetYaxis()->GetTitleOffset(), getPlottingAxis( "energy_Lin" )->fMinValue, getPlottingAxis( "energy_Lin" ) ->fMaxValue );
                        
    int z = 0;
    for( unsigned int i = 0; i < fData[iDataSetID]->hCutEfficiency.size(); i++ )
    {
        if( fData[iDataSetID]->hCutEfficiency[i] )
        {
            fData[iDataSetID]->hCutEfficiency[i]->Draw( "same" );
            cout << i + 1 << "\t" << fData[iDataSetID]->hCutEfficiency[i]->GetName();
            cout << " (color: " << fData[iDataSetID]->hCutEfficiency[i]->GetMarkerColor();
            cout << ", marker: " << fData[iDataSetID]->hCutEfficiency[i]->GetMarkerStyle() << ")";
            cout << "; " << fData[iDataSetID]->hCutEfficiency[i]->GetEntries();
            cout << endl;
            z++;
        }
    }
    if( z > 0 )
    {
        iCutEfficencyPlottingCanvas->SetLogy( 1 );
    }
}

TCanvas* VPlotInstrumentResponseFunction::plotEnergyReconstructionBias2D( unsigned int iDataSetID,
        double iYmin, double iYmax,
        bool iNoDirectionCut )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return 0;
    }
    
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "cREA_Eerr_%d", iDataSetID );
    sprintf( htitle, "relative error in energy reconstruction (%d)", iDataSetID );
    TCanvas* iEnergyReconstructionErrorCanvas = new TCanvas( hname, htitle, 610, 10, fCanvasSize_X, fCanvasSize_Y );
    iEnergyReconstructionErrorCanvas->SetGridx( 0 );
    iEnergyReconstructionErrorCanvas->SetGridy( 0 );
    iEnergyReconstructionErrorCanvas->Draw();
    
    TH2D* h2 = 0;
    if( iNoDirectionCut )
    {
        h2 = fData[iDataSetID]->hEsysMCRelative2DNoDirectionCut;
    }
    else
    {
        h2 = fData[iDataSetID]->hEsysMCRelative2D;
    }
    
    if( h2 )
    {
        h2->SetTitle( "" );
        h2->GetYaxis()->SetTitleOffset( 1.2 );
        h2->SetStats( 0 );
        h2->SetAxisRange( iYmin, iYmax, "Y" );
        if( log10( getPlottingAxis( "energy_Lin" )->fMinValue ) < h2->GetXaxis()->GetXmin() )
        {
            getPlottingAxis( "energy_Lin" )->fMinValue = TMath::Power( 10., h2->GetXaxis()->GetXmin() );
        }
        if( log10( getPlottingAxis( "energy_Lin" )->fMaxValue ) > h2->GetXaxis()->GetXmax() )
        {
            getPlottingAxis( "energy_Lin" )->fMaxValue = TMath::Power( 10., h2->GetXaxis()->GetXmax() );
        }
        
        h2->SetAxisRange( log10( getPlottingAxis( "energy_Lin" )->fMinValue ),
                          log10( getPlottingAxis( "energy_Lin" )->fMaxValue ), "X" );
        if( h2->GetEntries() > 0. )
        {
            iEnergyReconstructionErrorCanvas->SetLogz( 1 );
        }
        h2->Draw( "colz" );
        // line at 1
        TLine* iL = new TLine( log10( getPlottingAxis( "energy_Lin" )->fMinValue ), 1., log10( getPlottingAxis( "energy_Lin" )->fMaxValue ), 1. );
        iL->SetLineStyle( 2 );
        iL->Draw();
    }
    else
    {
        cout << "No reconstruction bias histogram found" << endl;
    }
    return iEnergyReconstructionErrorCanvas;
}


void VPlotInstrumentResponseFunction::plotEnergyReconstructionLogBias2D( unsigned int iDataSetID, string iM, double iYmin, double iYmax )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return;
    }
    
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "cEA_Eerr_%d", iDataSetID );
    sprintf( htitle, "error in energy reconstruction (%d)", iDataSetID );
    TCanvas* iEnergyReconstructionErrorCanvas = new TCanvas( hname, htitle, 610, 10, fCanvasSize_X, fCanvasSize_Y );
    iEnergyReconstructionErrorCanvas->SetGridx( 0 );
    iEnergyReconstructionErrorCanvas->SetGridy( 0 );
    iEnergyReconstructionErrorCanvas->Draw();
    
    if( fData[iDataSetID]->hEsys )
    {
        fData[iDataSetID]->hEsys->SetTitle( "" );
        fData[iDataSetID]->hEsys->GetYaxis()->SetTitleOffset( 1.2 );
        fData[iDataSetID]->hEsys->SetStats( 0 );
        fData[iDataSetID]->hEsys->SetAxisRange( iYmin, iYmax, "Y" );
        fData[iDataSetID]->hEsys->SetAxisRange( log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ), "X" );
        if( fData[iDataSetID]->hEsys->GetEntries() > 0. )
        {
            iEnergyReconstructionErrorCanvas->SetLogz( 1 );
        }
        fData[iDataSetID]->hEsys->Draw( "colz" );
        // plot energy systematics
        if( iM == "mean"     && fData[iDataSetID]->gEnergyLogBias_Mean )
        {
            fData[iDataSetID]->gEnergyLogBias_Mean->Draw( "p" );
        }
        if( iM == "median"   && fData[iDataSetID]->gEnergyLogBias_Median )
        {
            fData[iDataSetID]->gEnergyLogBias_Median->Draw( "p" );
        }
        
        // line at 0
        TLine* iL = new TLine( log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), 0., log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ), 0. );
        iL->SetLineStyle( 2 );
        iL->Draw();
    }
}

void VPlotInstrumentResponseFunction::printResponseMatrixTypes()
{
    cout << "Response matrix types: " << endl;
    cout << "\t default" << endl;
    cout << "\t QC" << endl;
    cout << "\t noTheta2Cut" << endl;
}

/*
 *
 * plot energy migration matrix
 *
 * iMatrixType = 'default'
 * iMatrixType = 'QC'
 * iMatrixType = 'noTheta2Cut'
 *
 */
void VPlotInstrumentResponseFunction::plotEnergyReconstructionMatrix( unsigned int iDataSetID, bool bFineBinning,
        string iMatrixType, bool bInterPol,
        bool bPlotMedian )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return;
    }
    
    ostringstream hname;
    hname << "cEA_Ematrix_" << iDataSetID << "_" << bFineBinning << "_" << iMatrixType << "_" << bInterPol;
    ostringstream htitle;
    htitle << "energy reconstruction matrix (" << iDataSetID << "," << bInterPol << ")";
    if( bFineBinning )
    {
        htitle << " (fine binning)";
    }
    if( iMatrixType != "default" )
    {
	htitle << ", " << iMatrixType;
    }
    TCanvas* iEnergyReconstructionMatrixCanvas = new TCanvas( hname.str().c_str(), htitle.str().c_str(), 610, 10, fCanvasSize_X, fCanvasSize_Y );
    iEnergyReconstructionMatrixCanvas->SetGridx( 0 );
    iEnergyReconstructionMatrixCanvas->SetGridy( 0 );
    iEnergyReconstructionMatrixCanvas->SetLeftMargin( 0.11 );
    iEnergyReconstructionMatrixCanvas->SetRightMargin( 0.13 );
    
    TH2D* i_hRecMatrix = 0;
    if( iMatrixType == "default" )
    {
        if( bFineBinning )
        {
            i_hRecMatrix = fData[iDataSetID]->hERecMatrix;
        }
        else
        {
            i_hRecMatrix = fData[iDataSetID]->hERecMatrixCoarse;
        }
    }
    else if( iMatrixType == "QC" )
    {
        if( bFineBinning )
        {
            i_hRecMatrix = fData[iDataSetID]->hERecMatrixQC;
        }
        else
        {
            i_hRecMatrix = fData[iDataSetID]->hERecMatrixCoarseQC;
        }
    }
    else if( iMatrixType == "noTheta2Cut" )
    {
        if( bFineBinning )
        {
            i_hRecMatrix = fData[iDataSetID]->hERecMatrixNoDirectionCuts;
        }
        else
        {
            i_hRecMatrix = fData[iDataSetID]->hERecMatrixCoarseNoDirectionCuts;
        }
    }
    if( !i_hRecMatrix )
    {
        return;
    }
    
    i_hRecMatrix->SetTitle( "" );
    i_hRecMatrix->GetYaxis()->SetTitleOffset( 1.2 );
    i_hRecMatrix->SetStats( 0 );
    if( i_hRecMatrix->GetEntries() > 0. )
    {
        iEnergyReconstructionMatrixCanvas->SetLogz( 1 );
    }
    i_hRecMatrix->SetXTitle( "log_{10} energy_{rec} [TeV]" );
    i_hRecMatrix->SetYTitle( "log_{10} energy_{MC} [TeV]" );
    i_hRecMatrix->SetAxisRange( log10( getPlottingAxis( "energy_Lin" )->fMinValue ), log10( getPlottingAxis( "energy_Lin" )->fMaxValue ), "X" );
    i_hRecMatrix->SetAxisRange( log10( getPlottingAxis( "energy_Lin" )->fMinValue ), log10( getPlottingAxis( "energy_Lin" )->fMaxValue ), "Y" );
    if( i_hRecMatrix->GetMaximum() < 1.5 )
    {
        i_hRecMatrix->SetMaximum( 1.5 );
    }
    if( !bInterPol )
    {
        i_hRecMatrix->Draw( "colz" );
    }
    else
    {
        TH2D* i_InterPol = ( TH2D* )VHistogramUtilities::interpolateResponseMatrix( i_hRecMatrix );
        if( i_InterPol )
        {
            i_InterPol->Draw( "colz" );
        }
    }
    // diagonal
    TLine* iL = new TLine( log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ),
                           log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ),
                           log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ),
                           log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ) );
    iL->SetLineStyle( 2 );
    iL->Draw();
    
    if( !bPlotMedian )
    {
        return;
    }
    
    // plot median
    double xq[3];
    double yq[] = { 0.0,  0.0, 0.0  };
    TGraphAsymmErrors* gReco = new TGraphAsymmErrors( 1 );
    TH1D hOff( "hMeanOffset", "", 100, -2., 2. );
    int zReco = 0;
    for( int i = 1; i <= i_hRecMatrix->GetNbinsY(); i++ )
    {
        TH1F* h = ( TH1F* )i_hRecMatrix->ProjectionX( "p_y", i, i );
        if( h && h->GetEntries() > 0 )
        {
            xq[0] = 0.16;
            xq[1] = 0.50;
            xq[2] = 0.84;
            h->GetQuantiles( 3, yq, xq );
            // count number of filled bins
            unsigned int i_filled = 0;
            for( int j = 1; j <= h->GetNbinsX(); j++ )
            {
                if( h->GetBinContent( j ) > 0 )
                {
                    i_filled++;
                }
            }
            if( i_filled <= 5 )
            {
                hOff.Fill( yq[1] - i_hRecMatrix->GetYaxis()->GetBinCenter( i ) );
                gReco->SetPoint( zReco, 0., i_hRecMatrix->GetYaxis()->GetBinCenter( i ) );
            }
            else
            {
                gReco->SetPoint( zReco, yq[1], i_hRecMatrix->GetYaxis()->GetBinCenter( i ) );
                gReco->SetPointError( zReco, yq[1] - yq[0], yq[2] - yq[1], 0., 0. );
            }
            zReco++;
        }
    }
    // fill missing points
    double x = 0;
    double y = 0;
    hOff.GetQuantiles( 3, yq, xq );
    for( int i = 0; i < gReco->GetN(); i++ )
    {
        gReco->GetPoint( i, x, y );
        if( TMath::Abs( x ) < 1.e-7 )
        {
            gReco->SetPoint( i, y + yq[1], y );
            gReco->SetPointError( i, yq[2] - yq[1], yq[2] - yq[1], 0., 0. );
        }
    }
    gReco->SetLineColor( 4 );
    gReco->SetMarkerColor( 4 );
    gReco->SetMarkerStyle( 28 );
    gReco->Draw( "pl" );
    
}


TCanvas* VPlotInstrumentResponseFunction::plotCutEfficiencyRatio( unsigned int iDataSetID, unsigned int iCutID,
        double iPlotMinimum, double iPlotMaximum )
{
    if( !checkDataSetID( iDataSetID ) && iDataSetID < 999 )
    {
        return 0;
    }
    
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "cEA_cuteffratio_%d_%d", iDataSetID, iCutID );
    sprintf( htitle, "cut efficiency ratio (%d, %d)", iDataSetID, iCutID );
    TCanvas* iCutEfficencyRatioPlottingCanvas = new TCanvas( hname, htitle, 10, 10, 600, 600 );
    iCutEfficencyRatioPlottingCanvas->SetGridx( 1 );
    iCutEfficencyRatioPlottingCanvas->SetGridy( 1 );
    
    sprintf( hname, "hceffratio_%d", iDataSetID );
    TH1D* hceff = new TH1D( hname, "", 100, log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ) );
    hceff->SetStats( 0 );
    hceff->SetXTitle( "log_{10} energy [TeV]" );
    hceff->SetYTitle( "cut efficiency (ratio)" );
    hceff->SetMinimum( iPlotMinimum );
    hceff->SetMaximum( iPlotMaximum );
    hceff->Draw( "" );
    hceff->Draw( "AH" );
    
    plot_nullHistogram( iCutEfficencyRatioPlottingCanvas, hceff, getPlottingAxis( "energy_Lin" )->fLogAxis, false, hceff->GetYaxis()->GetTitleOffset(), getPlottingAxis( "energy_Lin" ) ->fMinValue, getPlottingAxis( "energy_Lin" ) ->fMaxValue );
    
    if( iDataSetID < 999 )
    {
        for( unsigned int i = 0; i < fData[iDataSetID]->hCutEfficiencyRelativePlots.size(); i++ )
        {
            if( fData[iDataSetID]->hCutEfficiencyRelativePlots[i] )
            {
                fData[iDataSetID]->hCutEfficiencyRelativePlots[i]->Draw( "same" );
                cout << i + 1 << "\t" << fData[iDataSetID]->hCutEfficiencyRelativePlots[i]->GetName();
                cout << " (color: " << fData[iDataSetID]->hCutEfficiencyRelativePlots[i]->GetMarkerColor();
                cout << ", marker: " << fData[iDataSetID]->hCutEfficiencyRelativePlots[i]->GetMarkerStyle() << ")";
                cout << endl;
            }
        }
    }
    else
    {
        for( unsigned int i = 0; i < fData.size(); i++ )
        {
            if( iCutID < fData[i]->hCutEfficiencyRelativePlots.size() && fData[i]->hCutEfficiencyRelativePlots[iCutID] )
            {
                fData[i]->hCutEfficiencyRelativePlots[iCutID]->SetMarkerStyle( 20 + i );
                fData[i]->hCutEfficiencyRelativePlots[iCutID]->Draw( "same" );
                cout << i + 1 << "\t" << fData[i]->hCutEfficiencyRelativePlots[iCutID]->GetName() << endl;
                cout << " (color: " << fData[i]->hCutEfficiencyRelativePlots[iCutID]->GetMarkerColor() << ")" << endl;
            }
        }
    }

    return iCutEfficencyRatioPlottingCanvas;
    
}

void VPlotInstrumentResponseFunction::listDataSets()
{
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        cout << i << "\t";
        cout << fData[i]->isZombie() << "\t";
        cout << fData[i]->fFile << "\t";
        cout << endl;
    }
}

unsigned int VPlotInstrumentResponseFunction::getNumberOfGoodDataSets()
{
    unsigned int z = 0;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i] && !fData[i]->isZombie() )
        {
            z++;
        }
    }
    return z;
}

TCanvas* VPlotInstrumentResponseFunction::plotEffectiveAreaRatio( unsigned int iDataSetID, double ymin, double ymax )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return 0;
    }
    
    char hname[200];
    
    sprintf( hname, "cEA_EFFRatio" );
    TCanvas* iEffectiveAreaRatioPlottingCanvas = new TCanvas( hname, "effective area ratio", 10, 10, fCanvasSize_X, fCanvasSize_Y );
    iEffectiveAreaRatioPlottingCanvas->SetGridx( 0 );
    iEffectiveAreaRatioPlottingCanvas->SetGridy( 0 );
    iEffectiveAreaRatioPlottingCanvas->cd();
    TLegend* iLegend = makeLegend();
    
    TH1D* heffR = new TH1D( "heffR", "", 100, log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ) );
    heffR->SetStats( 0 );
    heffR->SetXTitle( "log_{10} energy [TeV]" );
    heffR->SetYTitle( "ratio of effective areas" );
    heffR->SetMinimum( ymin );
    heffR->SetMaximum( ymax );
    heffR->Draw( "" );
    heffR->Draw( "AH" );
    
    plot_nullHistogram( iEffectiveAreaRatioPlottingCanvas, heffR, getPlottingAxis( "energy_Lin" )->fLogAxis, false, 1.3,
                        getPlottingAxis( "energy_Lin" ) ->fMinValue, getPlottingAxis( "energy_Lin" ) ->fMaxValue );
                        
    TLine* iL = new TLine( log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), 1., log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ), 1. );
    iL->SetLineWidth( 2 );
    iL->SetLineStyle( 1 );
    iL->Draw();
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( i == iDataSetID )
        {
            continue;
        }
        TGraphAsymmErrors* g = 0;
        
        if( fData[i]->fA_MC == "A_MC" )
        {
            fData[i]->calculateEffectiveAreaRatios( fData[iDataSetID]->gEffArea_MC );
            g = fData[i]->gEffArea_MC_Ratio;
        }
        else if( fData[i]->fA_MC == "A_MCnoTH2" || fData[i]->fA_MC == "A_MC_noTH2" )
        {
            fData[i]->calculateEffectiveAreaRatios( fData[iDataSetID]->gEffArea_MCnoTh2 );
            g = fData[i]->gEffArea_MCnoTh2_Ratio;
        }
        else if( fData[i]->fA_MC == "A_RECnoTH2" || fData[i]->fA_MC == "A_REC_noTH2" )
        {
            fData[i]->calculateEffectiveAreaRatios( fData[iDataSetID]->gEffArea_RecnoTh2 );
            g = fData[i]->gEffArea_RecnoTh2_Ratio;
        }
        else
        {
            fData[i]->calculateEffectiveAreaRatios( fData[iDataSetID]->gEffArea_Rec );
            g = fData[i]->gEffArea_Rec_Ratio;
        }
        
        if( !g )
        {
            continue;
        }
        if( iLegend )
        {
            iLegend->AddEntry( g, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
        }
        g->Draw( fData[i]->fPlotOption.c_str() );
    }
    if( iLegend )
    {
        iLegend->Draw();
    }
    return iEffectiveAreaRatioPlottingCanvas;
}

/*

    plot the energy resolution

    note the different ways to calculate the energy resolution

*/
TCanvas* VPlotInstrumentResponseFunction::plotEnergyResolution( double ymin, double ymax,
        TPad* iResolutionPad,
        bool iSmoothEnergyResolution,
        string iMatrixType )
{
    if( fDebug )
    {
        cout << "VPlotInstrumentResponseFunction::plotEnergyResolution " << ymax << endl;
    }
    
    // canvas
    char hname[200];
    sprintf( hname, "cEA_energyResolution_%s", iMatrixType.c_str() );
    TCanvas* iEnergyResolutionPlottingCanvas = 0;
    if( iResolutionPad )
    {
        iEnergyResolutionPlottingCanvas = ( TCanvas* )iResolutionPad;
    }
    else
    {
        char htitle[200];
        sprintf( htitle, "energy resolution (%s)", iMatrixType.c_str() );
        iEnergyResolutionPlottingCanvas = new TCanvas( hname, htitle, 10, 10, fCanvasSize_X, fCanvasSize_Y );
        iEnergyResolutionPlottingCanvas->SetGridx( 0 );
        iEnergyResolutionPlottingCanvas->SetGridy( 0 );
        iEnergyResolutionPlottingCanvas->SetLeftMargin( 0.15 );
        iEnergyResolutionPlottingCanvas->SetRightMargin( 0.07 );
    }
    if( !iEnergyResolutionPlottingCanvas )
    {
        return 0;
    }
    iEnergyResolutionPlottingCanvas->cd();
    TLegend* iLegend = makeLegend() ;
    
    // plotting frame
    TH1D* he0 = 0;
    if( !iEnergyResolutionPlottingCanvas->GetListOfPrimitives()->FindObject( "he0" ) )
    {
        if( getPlottingAxis( "energy_Lin" )->fMinValue > 0. && getPlottingAxis( "energy_Lin" )->fMaxValue > 0. )
        {
            sprintf( hname, "he0%s", iMatrixType.c_str() );
            he0 = new TH1D( hname, "", 100, log10( getPlottingAxis( "energy_Lin" )->fMinValue ),
                            log10( getPlottingAxis( "energy_Lin" )->fMaxValue ) );
            he0->SetStats( 0 );
            he0->SetXTitle( "log_{10} energy [TeV]" );
            he0->SetYTitle( "energy resolution" );
            he0->SetMinimum( ymin );
            he0->SetMaximum( ymax );
        }
        else
        {
            cout << "VPlotInstrumentResponseFunction::plotEnergyResolution error: negative energy axis: ";
            cout << getPlottingAxis( "energy_Lin" )->fMinValue << "\t";
            cout << getPlottingAxis( "energy_Lin" )->fMaxValue << endl;
            return 0;
        }
        plot_nullHistogram( iEnergyResolutionPlottingCanvas, he0, getPlottingAxis( "energy_Lin" )->fLogAxis,
                            false, he0->GetYaxis()->GetTitleOffset() * 1.3,
                            getPlottingAxis( "energy_Lin" )->fMinValue, getPlottingAxis( "energy_Lin" )->fMaxValue );
    }
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        TGraphErrors* iEnergyResolution = fData[i]->gEnergyResolution;
        if( iMatrixType == "QC" )
        {
            iEnergyResolution = fData[i]->gEnergyResolutionQC;
        }
        else if( iMatrixType == "noTheta2Cut" )
        {
            iEnergyResolution = fData[i]->gEnergyResolutionNoDirectionCuts;
        }
        if( iEnergyResolution )
        {
            // smooth the energy resolution graph to remove wiggles due to
            // local fluctuations
            if( iSmoothEnergyResolution )
            {
                // smooth this graph using a Kernel of width 0.4
                TGraphSmooth* Eres_Smooth = new TGraphSmooth( "s" );
                TGraph* Eres_Smooth_graph = Eres_Smooth->SmoothKern( iEnergyResolution, "normal", 0.4, 200 );
                if( Eres_Smooth_graph )
                {
                    Eres_Smooth_graph->SetLineStyle( iEnergyResolution->GetLineStyle() );
                    Eres_Smooth_graph->SetLineColor( iEnergyResolution->GetLineColor() );
                    Eres_Smooth_graph->SetLineWidth( iEnergyResolution->GetLineWidth() );
                    Eres_Smooth_graph->SetMarkerStyle( iEnergyResolution->GetMarkerStyle() );
                    Eres_Smooth_graph->SetMarkerColor( iEnergyResolution->GetMarkerColor() );
                    Eres_Smooth_graph->SetMarkerSize( iEnergyResolution->GetMarkerSize() );
                    Eres_Smooth_graph->Draw( fData[i]->fPlotOption.c_str() );
                    if( iLegend )
                    {
                        iLegend->AddEntry( Eres_Smooth_graph, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                    }
                }
            }
            else
            {
                iEnergyResolution->Draw( fData[i]->fPlotOption.c_str() );
                if( iLegend )
                {
                    iLegend->AddEntry( iEnergyResolution,
                                       fData[i]->fLegend.c_str(),
                                       fData[i]->fPlotOption.c_str() );
                }
            }
        }
    }
    if( iLegend )
    {
        iLegend->Draw();
    }
    
    return iEnergyResolutionPlottingCanvas;
}

void VPlotInstrumentResponseFunction::plotEnergySpectra( bool iWeighted, double iYMax, int iRebin )
{
    char hname[200];
    char htitle[200];
    
    if( iWeighted )
    {
        sprintf( hname, "cEA_energy" );
        sprintf( htitle, "energy spectra (spectral weighted)" );
    }
    else
    {
        sprintf( hname, "cEA_energyUW" );
        sprintf( htitle, "energy spectra (not spectral weighted)" );
    }
    TCanvas* iEnergySpectraPlottingCanvas = new TCanvas( hname, htitle, 10, 10, fCanvasSize_X, fCanvasSize_Y );
    iEnergySpectraPlottingCanvas->SetGridx( 0 );
    iEnergySpectraPlottingCanvas->SetGridy( 0 );
    iEnergySpectraPlottingCanvas->SetLeftMargin( 0.15 );
    iEnergySpectraPlottingCanvas->SetRightMargin( 0.07 );
    iEnergySpectraPlottingCanvas->cd();
    TLegend* iLegend = makeLegend();
    if( iWeighted )
    {
        sprintf( hname, "he0" );
    }
    else
    {
        sprintf( hname, "he0UW" );
    }
    TH1D* he0 = new TH1D( hname, "", 100, log10( getPlottingAxis( "energy_Lin" )->fMinValue ),
                          log10( getPlottingAxis( "energy_Lin" )->fMaxValue ) );
    he0->SetStats( 0 );
    he0->SetXTitle( "log_{10} energy [TeV]" );
    he0->SetYTitle( "number of events/bin" );
    he0->SetMinimum( 0.5 );
    if( fData.size() > 0 && fData[0]->hEmc && iYMax < 0. )
    {
        he0->SetMaximum( fData[0]->hEmc->GetMaximum() * 1.5 );
    }
    else
    {
        he0->SetMaximum( iYMax );
    }
    he0->Draw( "" );
    he0->Draw( "AH" );
    
    plot_nullHistogram( iEnergySpectraPlottingCanvas, he0, getPlottingAxis( "energy_Lin" )->fLogAxis, true, he0->GetYaxis()->GetTitleOffset() * 1.3,
                        getPlottingAxis( "energy_Lin" )->fMinValue, getPlottingAxis( "energy_Lin" )->fMaxValue );
                        
                        
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        bool addEntry = true;
        if( fData[i]->hEmc )
        {
            fData[i]->hEmc->Draw( "same" );
            if( iLegend && addEntry )
            {
                iLegend->AddEntry( fData[i]->hEmc, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
            }
            addEntry = false;
        }
        if( iWeighted )
        {
            if( fData[i]->hEcut )
            {
                fData[i]->hEcut->Draw( "same" );
                if( iLegend && addEntry )
                {
                    iLegend->AddEntry( fData[i]->hEcut, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                addEntry = false;
                
            }
            if( fData[i]->hEcut_rec )
            {
                fData[i]->hEcut_rec->Draw( "same" );
                if( iLegend && addEntry )
                {
                    iLegend->AddEntry( fData[i]->hEcut_rec, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                addEntry = false;
                
            }
        }
        else
        {
            if( fData[i]->hEcutUW )
            {
                if( iRebin > 1 )
                {
                    fData[i]->hEcutUW->Rebin( iRebin );
                    
                }
                fData[i]->hEcutUW->Draw( "same" );
                if( iLegend && addEntry )
                {
                    iLegend->AddEntry( fData[i]->hEcutUW, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                addEntry = false;
                
            }
            if( fData[i]->hEcut_recUW )
            {
                if( iRebin > 1 )
                {
                    fData[i]->hEcut_recUW->Rebin( iRebin );
                }
                fData[i]->hEcut_recUW->Draw( "same" );
                if( iLegend && addEntry )
                {
                    iLegend->AddEntry( fData[i]->hEcut_recUW, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                addEntry = false;
            }
        }
    }
    for( unsigned int i = 0; i < fData[0]->hCutEfficiency.size(); i++ )
    {
        if( fData[0]->hCutEfficiency[i] )
        {
            if( iRebin > 1 )
            {
                fData[0]->hCutEfficiency[i]->Rebin( iRebin );
            }
            fData[0]->hCutEfficiency[i]->Draw( "same" );
            cout << i + 1 << "\t" << fData[0]->hCutEfficiency[i]->GetName();
            cout << " (color: " << fData[0]->hCutEfficiency[i]->GetMarkerColor();
            cout << ", marker: " << fData[0]->hCutEfficiency[i]->GetMarkerStyle() << ")";
            cout << endl;
        }
    }
    if( iLegend )
    {
        iLegend->Draw();
    }
}

void VPlotInstrumentResponseFunction::plotEnergyReconstructionLogBias( string iM, double ymin, double ymax )
{
    plotEnergyReconstructionBias( iM, ymin, ymax, true );
}

TCanvas* VPlotInstrumentResponseFunction::plotEnergyReconstructionBias( string iM, double ymin, double ymax,
        bool iLogBias, TCanvas* iEnergySystematicsPlottingCanvas )
{
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "cEA_energy_bias_%d", ( int )iLogBias );
    if( iLogBias )
    {
        sprintf( htitle, "log energy bias" );
    }
    else
    {
        sprintf( htitle, "energy bias" );
    }
    TLegend* iLegend = 0;
    if( !iEnergySystematicsPlottingCanvas )
    {
        iEnergySystematicsPlottingCanvas = new TCanvas( hname, htitle, 10, 10, fCanvasSize_X, fCanvasSize_Y );
        iEnergySystematicsPlottingCanvas->SetGridx( 0 );
        iEnergySystematicsPlottingCanvas->SetGridy( 0 );
        iEnergySystematicsPlottingCanvas->SetLeftMargin( 0.15 );
        iEnergySystematicsPlottingCanvas->SetRightMargin( 0.07 );
        iEnergySystematicsPlottingCanvas->cd();
        iLegend = makeLegend();
        
        sprintf( hname, "he0_sys" );
        if( iLogBias )
        {
            sprintf( hname, "he0_sysL" );
        }
        TH1D* he0_sys = new TH1D( hname, "", 100, log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ),
                                  log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ) );
        he0_sys->SetStats( 0 );
        he0_sys->SetXTitle( "log_{10} energy [TeV]" );
        if( iLogBias )
        {
            he0_sys->SetYTitle( "energy bias (log_{10} E_{rec}/E_{MC})" );
        }
        else
        {
            he0_sys->SetYTitle( "energy bias E_{rec}/E_{MC}" );
        }
        he0_sys->SetMinimum( ymin );
        he0_sys->SetMaximum( ymax );
        he0_sys->Draw( "" );
        he0_sys->Draw( "AH" );
        
        plot_nullHistogram( iEnergySystematicsPlottingCanvas, he0_sys, getPlottingAxis( "energy_Lin" )->fLogAxis,
                            false, he0_sys->GetYaxis()->GetTitleOffset() * 1.3,
                            getPlottingAxis( "energy_Lin" )->fMinValue, getPlottingAxis( "energy_Lin" ) ->fMaxValue );
        if( !iLogBias )
        {
            TLine* iL = new TLine( log10( getPlottingAxis( "energy_Lin" ) ->fMinValue ), 1.,
                                   log10( getPlottingAxis( "energy_Lin" ) ->fMaxValue ), 1. );
            iL->SetLineStyle( 2 );
            iL->Draw();
        }
    }
    else
    {
        iEnergySystematicsPlottingCanvas->cd();
    }
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( iLogBias )
        {
            if( iM == "mean" && fData[i]->gEnergyLogBias_Mean )
            {
                fData[i]->gEnergyLogBias_Mean->Print();
                fData[i]->gEnergyLogBias_Mean->Draw( "p" );
                if( iLegend )
                {
                    iLegend->AddEntry( fData[i]->gEnergyLogBias_Mean, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                if( fDebug )
                {
                    fData[i]->gEnergyLogBias_Mean->Print();
                }
            }
            else if( iM == "median" && fData[i]->gEnergyLogBias_Median )
            {
                fData[i]->gEnergyLogBias_Median->Draw( "p" );
                if( iLegend )
                {
                    iLegend->AddEntry( fData[i]->gEnergyLogBias_Median, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                if( fDebug )
                {
                    fData[i]->gEnergyLogBias_Median->Print();
                }
            }
            else
            {
                cout << "no (log) graph found (" << iM << ")" << endl;
            }
        }
        else
        {
            if( iM == "mean" && fData[i]->gEnergyBias_Mean )
            {
                fData[i]->gEnergyBias_Mean->Print();
                fData[i]->gEnergyBias_Mean->Draw( "pl" );
                if( iLegend )
                {
                    iLegend->AddEntry( fData[i]->gEnergyBias_Mean, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                if( fDebug )
                {
                    fData[i]->gEnergyBias_Mean->Print();
                }
            }
            else if( iM == "median" && fData[i]->gEnergyBias_Median )
            {
                fData[i]->gEnergyBias_Median->Draw( "p" );
                if( iLegend )
                {
                    iLegend->AddEntry( fData[i]->gEnergyBias_Median, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
                }
                if( fDebug )
                {
                    fData[i]->gEnergyBias_Median->Print();
                }
            }
            else
            {
                cout << "no (lin) graph found (" << iM << ")" << endl;
            }
        }
    }
    if( iLegend )
    {
        iLegend->Draw() ;
    }
    return iEnergySystematicsPlottingCanvas;
}

TCanvas* VPlotInstrumentResponseFunction::plotAngularResolution2D( unsigned int iDataSetID, string iXaxis, string iProbabilityString, double iEnergySlice_TeV )
{
    string iResolutionTreeName = "t_angular_resolution";
    if( iProbabilityString != "68" )
    {
        iResolutionTreeName += "_0" + iProbabilityString + "p";
    }
    return plotResolution2D( iDataSetID, "angres" + iProbabilityString, "angular resolution vs " + iXaxis,
                             "angular resolution (" + iProbabilityString + "%) [deg]",
                             getPlottingAxis( "angularesolution_Lin" )->fMinValue, getPlottingAxis( "angularesolution_Lin" )->fMaxValue,
                             iResolutionTreeName, iXaxis, iEnergySlice_TeV );
}

TCanvas* VPlotInstrumentResponseFunction::plotAngularResolution( string iXaxis, string iProbabilityString,
        double iMin, double iMax,
        TPad* iResolutionPad, bool iLogY )
{
    string iResolutionTreeName = "t_angular_resolution";
    if( iMax > 0. )
    {
        getPlottingAxis( "angularesolution_Lin" )->fMaxValue = iMax;
        getPlottingAxis( "angularesolution_Log" )->fMaxValue = iMax;
    }
    else
    {
        if( iLogY )
        {
            iMax = getPlottingAxis( "angularesolution_Log" )->fMaxValue;
        }
        else
        {
            iMax = getPlottingAxis( "angularesolution_Lin" )->fMaxValue;
        }
    }
    if( iMin > 0. )
    {
        getPlottingAxis( "angularesolution_Lin" )->fMinValue = iMin;
        getPlottingAxis( "angularesolution_Log" )->fMinValue = iMin;
    }
    else
    {
        if( iLogY )
        {
            iMin = getPlottingAxis( "angularesolution_Log" )->fMinValue;
        }
        else
        {
            iMin = getPlottingAxis( "angularesolution_Lin" )->fMinValue;
        }
    }
    if( iProbabilityString != "68" )
    {
        iResolutionTreeName += "_0" + iProbabilityString + "p";
    }
    
    return plotResolution( "angres"  + iProbabilityString, "angular resolution vs " + iXaxis + " (" + iProbabilityString + "%)",
                           "angular resolution [deg]", iMin, iMax,
                           iResolutionTreeName, iXaxis, iResolutionPad, iLogY );
}

TCanvas* VPlotInstrumentResponseFunction::plotCoreResolution( string iXaxis, double iMax )
{
    if( iMax > 0. )
    {
        getPlottingAxis( "coreresolution_Lin" )->fMaxValue = iMax;
    }
    return plotResolution( "coreres", "core resolution vs " + iXaxis, "core resolution [m]",
                           getPlottingAxis( "coreresolution_Lin" )->fMinValue,
                           getPlottingAxis( "coreresolution_Lin" )->fMaxValue, "t_core_resolution", iXaxis );
}

TCanvas* VPlotInstrumentResponseFunction::plotCoreResolution2D( unsigned int iDataSetID, string iXaxis )
{
    return plotResolution2D( iDataSetID, "coreres", "core resolution vs " + iXaxis, "core resolution [m]",
                             getPlottingAxis( "coreresolution_Lin" )->fMinValue,
                             getPlottingAxis( "coreresolution_Lin" )->fMaxValue, "t_core_resolution", iXaxis );
}

/*
TCanvas* VPlotInstrumentResponseFunction::plotEnergyResolution( string iXaxis )
{
	iXaxis = "bubabub";
	return plotResolution( "energyres", "energy resolution vs energy", "energy resolution",
						   getPlottingAxis( "energyresolution_Lin" )->fMinValue,
						   getPlottingAxis( "energyresolution_Lin" )->fMaxValue, "t_energy_resolution", "energy" );
} */

TCanvas* VPlotInstrumentResponseFunction::plotEnergyResolution2D( unsigned int iDataSetID )
{
    return plotResolution2D( iDataSetID, "energyres", "energy resolution vs energy", "energy resolution",
                             getPlottingAxis( "energyresolution_Lin" )->fMinValue,
                             getPlottingAxis( "energyresolution_Lin" )->fMaxValue, "t_energy_resolution", "energy" );
}

TCanvas* VPlotInstrumentResponseFunction::plotResolution2D( unsigned int iDataSetID, string iName,
        string iCanvasTitle, string iYTitle,
        double iYmin, double iYmax,
        string iResolutionTreeName, string iXaxis,
        double iEnergySlice_TeV )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return 0;
    }
    
    if( fData.size() == 0 )
    {
        return 0;
    }
    
    char hname[800];
    char htitle[800];
    
    unsigned int i_Plotting_Selector = 0;
    double i_Plotting_Min = 0.;
    double i_Plotting_Max = 0.;
    double i_Plotting_L_Min = 0.;
    double i_Plotting_L_Max = 0.;
    bool   i_Plotting_log = false;
    string iXTitle = "";
    if( iXaxis == "energy" )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_DIFF;
    }
    else if( iXaxis == "nimages" )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_NIMAG;
    }
    else if( iXaxis == "distance" )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_DIST;
    }
    string iPlottingAxis = iXaxis + "_Lin";
    if( getPlottingAxis( iPlottingAxis ) )
    {
        i_Plotting_Min = getPlottingAxis( iPlottingAxis )->fMinValue;
        i_Plotting_Max = getPlottingAxis( iPlottingAxis )->fMaxValue;
        if( iXaxis == "energy" )
        {
            if( getPlottingAxis( iPlottingAxis )->fMinValue > 0. )
            {
                i_Plotting_L_Min = log10( getPlottingAxis( iPlottingAxis )->fMinValue );
            }
            else
            {
                cout << "VPlotInstrumentResponseFunction::plotResolution2D log min axis range <=0: " << getPlottingAxis( iPlottingAxis )->fMinValue;
                return 0;
            }
            if( getPlottingAxis( iPlottingAxis )->fMaxValue > 0. )
            {
                i_Plotting_L_Max = log10( getPlottingAxis( iPlottingAxis )->fMaxValue );
            }
            else
            {
                cout << "VPlotInstrumentResponseFunction::plotResolution2D log max axis range <=0: " << getPlottingAxis( iPlottingAxis )->fMaxValue;
                return 0;
            }
        }
        else
        {
            i_Plotting_L_Min = i_Plotting_Min;
            i_Plotting_L_Max = i_Plotting_Max;
        }
        i_Plotting_log = getPlottingAxis( iPlottingAxis )->fLogAxis;
        iXTitle = getPlottingAxis( iPlottingAxis )->fAxisTitle;
    }
    else
    {
        cout << "VPlotInstrumentResponseFunction::plotResolution2D:: X-axis not found" << endl;
        cout << "(available X-axes: energy, nimages, distance)" << endl;
        return 0;
    }
    if( fDebug )
    {
        cout << "Axis range: ";
        cout << "X " << i_Plotting_L_Min << "\t" << i_Plotting_L_Max << "\t" << i_Plotting_log;
        cout << ", Y " << iYmin << "\t" << iYmax;
        cout << endl;
    }
    
    // create canvas
    sprintf( hname, "c2D%s_%s_%d_%d", iName.c_str(), iXaxis.c_str(), iDataSetID, ( int )iEnergySlice_TeV );
    sprintf( htitle, "%s (data set %d)", iCanvasTitle.c_str(), iDataSetID );
    TCanvas* iResolutionPlottingCanvas = new TCanvas( hname, htitle, 210, 10, fCanvasSize_X, fCanvasSize_Y );
    iResolutionPlottingCanvas->SetGridx( 0 );
    iResolutionPlottingCanvas->SetGridy( 0 );
    iResolutionPlottingCanvas->SetLeftMargin( 0.13 );
    iResolutionPlottingCanvas->SetRightMargin( 0.13 );
    
    // get 2D histo
    TH2D* h = 0;
    for( unsigned int j = 0; j < fData[iDataSetID]->fIRF_TreeNames.size(); j++ )
    {
        if( fData[iDataSetID]->fIRF_TreeNames[j] == iResolutionTreeName )
        {
            if( j < fData[iDataSetID]->fIRF_Data.size() && fData[iDataSetID]->fIRF_Data[j] &&
                    i_Plotting_Selector < fData[iDataSetID]->fIRF_Data[j]->f2DHisto.size() )
            {
                h = fData[iDataSetID]->fIRF_Data[j]->f2DHisto[i_Plotting_Selector];
            }
            // plot everything
            if( h )
            {
                // get containment probability
                if( j < fData[iDataSetID]->fIRF_Data.size() && fData[iDataSetID]->fIRF_Data[j]
                        && i_Plotting_Selector < fData[iDataSetID]->fIRF_Data[j]->fContainmentProbability.size() )
                {
                    sprintf( hname, "%s (%d%%)", h->GetYaxis()->GetTitle(),
                             ( int )( fData[iDataSetID]->fIRF_Data[j]->fContainmentProbability[i_Plotting_Selector] * 100. ) );
                    h->SetYTitle( hname );
                }
                setHistogramPlottingStyle( h, -99. );
                h->SetAxisRange( i_Plotting_L_Min, i_Plotting_L_Max, "X" );
                h->SetAxisRange( iYmin, iYmax, "Y" );
                h->SetYTitle( iYTitle.c_str() );
                // plot 2D histogram
                if( iEnergySlice_TeV < 0. )
                {
                    h->GetYaxis()->SetTitleOffset( 1.5 );
                    h->Draw( "colz" );
                }
                // plot a slice of the 2D histogram
                else
                {
                    h->SetAxisRange( -1., -1., "Y" );
                    sprintf( hname, "%s_%d_%d", h->GetName(), ( int )iEnergySlice_TeV, h->GetXaxis()->FindBin( log10( iEnergySlice_TeV ) ) );
                    TH1D* h1D = h->ProjectionY( hname, h->GetXaxis()->FindBin( log10( iEnergySlice_TeV ) ), h->GetXaxis()->FindBin( log10( iEnergySlice_TeV ) ) );
                    if( h1D )
                    {
                        h1D->GetXaxis()->SetTitleOffset( 1.2 );
                        //		      h1D = get_Cumulative_Histogram( h1D, true, true );
                        h1D->Draw();
                    }
                }
                if( fDebug )
                {
                    cout << "HISTOGRAM " << h->GetName() << endl;
                    cout << "\t entries: " << h->GetEntries() << endl;
                }
                break;
            }
        }
    }
    if( !h )
    {
        cout << "VPlotInstrumentResponseFunction::plotResolution2D() warning: no histogram found for data set " << iDataSetID << endl;
    }
    
    return iResolutionPlottingCanvas;
}

/*!

   plot e.g. angular resolution

*/
TCanvas*  VPlotInstrumentResponseFunction::plotResolution( string iName, string iCanvasTitle, string iYTitle,
        double iYmin, double iYmax,
        string iResolutionTreeName,
        string iXaxis, TPad* iResolutionPad,
        bool iYLog )
{
    if( fData.size() == 0 )
    {
        return 0;
    }
    
    char hname[200];
    
    unsigned int i_Plotting_Selector = 0;
    // plotting axis
    double i_Plotting_Min = 0.;
    double i_Plotting_Max = 0.;
    double i_Plotting_L_Min = 0.;
    double i_Plotting_L_Max = 0.;
    bool   i_Plotting_log = false;
    string iXTitle = "";
    if( iXaxis == "energy" )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_DIFF;
    }
    else if( iXaxis == "nimages" )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_NIMAG;
    }
    else if( iXaxis == "distance" )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_DIST;
    }
    string iPlottingAxis = iXaxis + "_Lin";
    if( getPlottingAxis( iPlottingAxis ) )
    {
        i_Plotting_Min = getPlottingAxis( iPlottingAxis )->fMinValue;
        i_Plotting_Max = getPlottingAxis( iPlottingAxis )->fMaxValue;
        if( iXaxis == "energy" )
        {
            if( getPlottingAxis( iPlottingAxis )->fMinValue > 0. )
            {
                i_Plotting_L_Min = log10( getPlottingAxis( iPlottingAxis )->fMinValue );
            }
            else
            {
                cout << "VPlotInstrumentResponseFunction::plotResolution log min axis range <=0: " << getPlottingAxis( iPlottingAxis )->fMinValue;
                return 0;
            }
            if( getPlottingAxis( iPlottingAxis )->fMaxValue > 0. )
            {
                i_Plotting_L_Max = log10( getPlottingAxis( iPlottingAxis )->fMaxValue );
            }
            else
            {
                cout << "VPlotInstrumentResponseFunction::plotResolution log max axis range <=0: " << getPlottingAxis( iPlottingAxis )->fMaxValue;
                return 0;
            }
            iXTitle = "log_{10} energy [TeV]";
        }
        else
        {
            i_Plotting_L_Min = i_Plotting_Min;
            i_Plotting_L_Max = i_Plotting_Max;
            iXTitle = getPlottingAxis( iPlottingAxis )->fAxisTitle;
        }
        i_Plotting_log = getPlottingAxis( iPlottingAxis )->fLogAxis;
    }
    else
    {
        cout << "VPlotInstrumentResponseFunction::plotResolution:: X-axis not found" << endl;
        cout << "(available X-axes: energy, nimages, distance)" << endl;
        return 0;
    }
    if( fDebug )
    {
        cout << "Axis range: ";
        cout << "X " << i_Plotting_L_Min << "\t" << i_Plotting_L_Max << "\t" << i_Plotting_log;
        cout << ", Y " << iYmin << "\t" << iYmax;
        cout << endl;
    }
    
    // create canvas
    TCanvas* iResolutionPlottingCanvas = 0;
    bool i_renewPlot = true;
    if( iResolutionPad )
    {
        iResolutionPlottingCanvas = ( TCanvas* )iResolutionPad;
        i_renewPlot = false;
    }
    else
    {
        sprintf( hname, "c%s_%s", iName.c_str(), iXaxis.c_str() );
        iResolutionPlottingCanvas = new TCanvas( hname, iCanvasTitle.c_str(), 210, 10, fCanvasSize_X, fCanvasSize_Y );
        iResolutionPlottingCanvas->SetGridx( 0 );
        iResolutionPlottingCanvas->SetGridy( 0 );
        iResolutionPlottingCanvas->SetLeftMargin( 0.15 );
        iResolutionPlottingCanvas->SetRightMargin( 0.07 );
    }
    if( !iResolutionPlottingCanvas )
    {
        return 0;
    }
    iResolutionPlottingCanvas->cd();
    TLegend* iLegend = makeLegend();
    
    sprintf( hname, "har_%s_%s", iName.c_str(), iXaxis.c_str() );
    TH1D* har = 0;
    if( !iResolutionPlottingCanvas->GetListOfPrimitives()->FindObject( hname ) )
    {
        har = new TH1D( hname, "", 100, i_Plotting_L_Min, i_Plotting_L_Max );
        har->SetStats( 0 );
        har->SetXTitle( iXTitle.c_str() );
        har->SetYTitle( iYTitle.c_str() );
        har->SetMinimum( iYmin );
        har->SetMaximum( iYmax );
        har->Draw( "" );
        har->Draw( "AH" );
        
        plot_nullHistogram( iResolutionPlottingCanvas, har, i_Plotting_log, false, 1.6, i_Plotting_Min, i_Plotting_Max );
    }
    else
    {
        har = ( TH1D* )iResolutionPlottingCanvas->GetListOfPrimitives()->FindObject( hname );
    }
    
    // get resolution graphs for the whole data sample
    int z = 0;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        // get graph
        TGraphErrors* g = 0;
        for( unsigned int j = 0; j < fData[i]->fIRF_TreeNames.size(); j++ )
        {
            if( fData[i]->fIRF_TreeNames[j] == iResolutionTreeName )
            {
                if( j < fData[i]->fIRF_Data.size()
                        && fData[i]->fIRF_Data[j]
                        && i_Plotting_Selector < fData[i]->fIRF_Data[j]->fResolutionGraph.size() )
                {
                    // try to get EVNDISP resolution graph
                    g = fData[i]->fIRF_Data[j]->fResolutionGraph[i_Plotting_Selector];
                }
                // get containment probability (for first data set only)
                if( g && i == 0 )
                {
                    if( j < fData[i]->fIRF_Data.size() && fData[i]->fIRF_Data[j]
                            && i_Plotting_Selector < fData[i]->fIRF_Data[j]->fContainmentProbability.size() )
                    {
                        sprintf( hname, "%s (%d%%)", har->GetYaxis()->GetTitle(),
                                 ( int )( fData[i]->fIRF_Data[j]->fContainmentProbability[i_Plotting_Selector] * 100. ) );
                        har->SetYTitle( hname );
                        if( i_renewPlot )
                        {
                            plot_nullHistogram( iResolutionPlottingCanvas, har, i_Plotting_log, false, 1.6, i_Plotting_Min, i_Plotting_Max );
                        }
                    }
                }
                else if( fDebug )
                {
                    cout << "VPlotInstrumentResponseFunction::plotResolution() no graph found" << endl;
                }
            }
        }
        // plot null histogram
        if( i == 0 && i_renewPlot )
        {
            plot_nullHistogram( iResolutionPlottingCanvas, har, i_Plotting_log, iYLog, 1.6, i_Plotting_Min, i_Plotting_Max );
        }
        /////////////////////////////////////////////////////
        //// *** important for CTA is the following *** /////
        /////////////////////////////////////////////////////
        if( !g )
        {
            // try to get CTA resolution graph
            if( iName.find( "angres" ) != string::npos )
            {
                if( iName.find( "80" ) != string::npos )
                {
                    g = fData[i]->gAngularResolution80;
                }
                else if( iName.find( "95" ) != string::npos )
                {
                    g = fData[i]->gAngularResolution95;
                }
                else
                {
                    g = fData[i]->gAngularResolution;
                }
            }
            if( !g )
            {
                cout << "VPlotInstrumentResponseFunction::plotResolution() warning: no graph found for data set " << i << endl;
                continue;
            }
        }
        else if( g->GetN() == 0 )
        {
            cout << "VPlotInstrumentResponseFunction::plotResolution() warning: graph without points in data set " << i << endl;
            cout << "(" << fData[i]->fIRF_TreeNames.size() << ", " << iResolutionTreeName << ")" << endl;
            cout << "(";
            for( unsigned int j = 0; j < fData[i]->fIRF_TreeNames.size(); j++ )
            {
                cout << fData[i]->fIRF_TreeNames[j] << " ";
            }
            cout << ")" << endl;
            continue;
        }
        // draw the resolution graph
        if( g )
        {
            if( fDebug )
            {
                g->Print();
            }
            g->Draw( fData[i]->fPlotOption.c_str() );
            if( fFunction_fitResolution.size() > 0 )
            {
                fitResolution( g );
            }
            if( iLegend )
            {
                iLegend->AddEntry( g, fData[i]->fLegend.c_str(), fData[i]->fPlotOption.c_str() );
            }
            z++;
        }
        gLastPlottedGraph = g;
    }
    if( iLegend )
    {
        iLegend->Draw();
    }
    if( iYLog )
    {
        iResolutionPlottingCanvas->SetLogy( 1 );
    }
    return iResolutionPlottingCanvas;
}

bool VPlotInstrumentResponseFunction::setResolutionFitting( string iFitFunction, double iFitXmin, double iFitXmax )
{
    fFunction_fitResolution = iFitFunction;
    fXmin_fitResolution = iFitXmin;
    fXmax_fitResolution = iFitXmax;
    
    return true;
}

bool VPlotInstrumentResponseFunction::fitResolution( TGraphErrors* g )
{
    if( !g )
    {
        return false;
    }
    
    char hname[400];
    sprintf( hname, "fitResolution" );
    fTF1_fitResolution = new TF1( hname, fFunction_fitResolution.c_str(), fXmin_fitResolution, fXmax_fitResolution );
    
    g->Fit( fTF1_fitResolution, "R" );
    
    return true;
}

bool VPlotInstrumentResponseFunction::write_fitResolutionFunction( string iOutName, string iName )
{
    if( !fTF1_fitResolution )
    {
        cout << "PlotInstrumentResponseFunction::write_fitResolutionFunction error: no function defined" << endl;
        return false;
    }
    TFile f( iOutName.c_str(), "RECREATE" );
    if( f.IsZombie() )
    {
        cout << "PlotInstrumentResponseFunction::write_fitResolutionFunction error opening output file: " << iOutName << endl;
        return false;
    }
    if( iName.size() > 0 )
    {
        fTF1_fitResolution->SetName( iName.c_str() );
    }
    fTF1_fitResolution->Write();
    f.Close();
    
    return true;
}

TH1D* VPlotInstrumentResponseFunction::getTheta2orThetaHistogram( unsigned int iDataSetID, double i_Energy_TeV_lin,
        bool iTheta2 )
{
    if( !checkDataSetID( iDataSetID ) )
    {
        return 0;
    }
    
    if( fData.size() == 0 )
    {
        return 0;
    }
    
    string iResolutionTreeName = "t_angular_resolution";
    
    // check if theta or theta2 histograms should be used
    unsigned int i_Plotting_Selector = VInstrumentResponseFunctionData::E_DIFF2;
    if( !iTheta2 )
    {
        i_Plotting_Selector = VInstrumentResponseFunctionData::E_DIFF;
    }
    
    
    // get 2D histo
    TH2D* h2D = 0;
    TH1D* h1D = 0;
    for( unsigned int j = 0; j < fData[iDataSetID]->fIRF_TreeNames.size(); j++ )
    {
        if( fData[iDataSetID]->fIRF_TreeNames[j] == iResolutionTreeName )
        {
            if( j < fData[iDataSetID]->fIRF_Data.size() && fData[iDataSetID]->fIRF_Data[j] && i_Plotting_Selector < fData[iDataSetID]->fIRF_Data[j]->f2DHisto.size() )
            {
                h2D = fData[iDataSetID]->fIRF_Data[j]->f2DHisto[i_Plotting_Selector];
                if( h2D )
                {
                    if( i_Energy_TeV_lin > 0. )
                    {
                        char iname[600];
                        sprintf( iname, "%s_%d_%f", h2D->GetName(), iDataSetID, i_Energy_TeV_lin );
                        h1D = h2D->ProjectionY( iname, h2D->GetXaxis()->FindBin( log10( i_Energy_TeV_lin ) ),
                                                h2D->GetXaxis()->FindBin( log10( i_Energy_TeV_lin ) ) );
                        if( iTheta2 )
                        {
                            setHistogramPlottingStyle( h1D, iDataSetID + 1, 3. );
                        }
                        else
                        {
                            setHistogramPlottingStyle( h1D, iDataSetID + 1, 1. );
                        }
                    }
                    else
                    {
                        h1D = 0;
                    }
                }
            }
        }
    }
    
    return h1D;
}

TCanvas* VPlotInstrumentResponseFunction::plotTheta2( double iTheta2AxisMax, bool iCumulative )
{
    vector< double > i_temp_vector;
    return plotPSF( i_temp_vector, iTheta2AxisMax, iCumulative, true );
}

TCanvas* VPlotInstrumentResponseFunction::plotTheta2( vector< double > i_Energy_TeV_lin, double iTheta2AxisMax, bool iCumulative )
{
    return plotPSF( i_Energy_TeV_lin, iTheta2AxisMax, iCumulative, true );
}

TCanvas* VPlotInstrumentResponseFunction::plotTheta( double iTheta2AxisMax, bool iCumulative )
{
    vector< double > i_temp_vector;
    return plotPSF( i_temp_vector, iTheta2AxisMax, iCumulative, false );
}

TCanvas* VPlotInstrumentResponseFunction::plotTheta( vector< double > i_Energy_TeV_lin, double iTheta2AxisMax, bool iCumulative )
{
    return plotPSF( i_Energy_TeV_lin, iTheta2AxisMax, iCumulative, false );
}


TCanvas* VPlotInstrumentResponseFunction::plotPSF( vector< double > i_Energy_TeV_lin, double iTheta2AxisMax,
        bool iCumulative, bool iPlotTheta2 )
{

    if( i_Energy_TeV_lin.size() == 0 )
    {
        i_Energy_TeV_lin.push_back( 0.2 );
        i_Energy_TeV_lin.push_back( 0.5 );
        i_Energy_TeV_lin.push_back( 1.0 );
        i_Energy_TeV_lin.push_back( 5.0 );
    }
    
    char hname[600];
    if( iPlotTheta2 )
    {
        sprintf( hname, "Theta2_ID_%d", iCumulative );
    }
    else
    {
        sprintf( hname, "Theta_ID_%d", iCumulative );
    }
    TCanvas* c = new TCanvas( hname, hname, 10, 10, 600, 600 );
    c->Divide( TMath::Nint( sqrt( i_Energy_TeV_lin.size() ) ), TMath::Nint( sqrt( i_Energy_TeV_lin.size() ) ) );
    for( unsigned int j = 0; j < i_Energy_TeV_lin.size(); j++ )
    {
        c->cd( j + 1 );
        gPad->SetGridx( 0 );
        gPad->SetGridy( 0 );
        // histogram frame
        if( iPlotTheta2 )
        {
            sprintf( hname, "hTheta2_ID_%d_%d", iCumulative, j );
        }
        else
        {
            sprintf( hname, "hTheta_ID_%d_%d", iCumulative, j );
        }
        TH1D* hnull = new TH1D( hname, "", 100, 0., iTheta2AxisMax );
        if( iPlotTheta2 )
        {
            hnull->SetXTitle( "#Theta^{2}" );
        }
        else
        {
            hnull->SetXTitle( "#Theta" );
        }
        hnull->SetMaximum( 1.1 );
        hnull->SetStats( 0 );
        plot_nullHistogram( ( TPad* )gPad, hnull, false, false, 1.3, 0., iTheta2AxisMax );
        hnull->GetXaxis()->SetNdivisions( 505 );
        hnull->Draw();
        
        if( i_Energy_TeV_lin[j] < 0.1 )
        {
            sprintf( hname, "%.2f TeV", i_Energy_TeV_lin[j] );
        }
        else
        {
            sprintf( hname, "%.1f TeV", i_Energy_TeV_lin[j] );
        }
        TText* iT = new TText( iTheta2AxisMax * 0.6, hnull->GetMaximum() * 0.5, hname );
        iT->Draw();
        
    }
    
    // loop over all data sets
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        for( unsigned int j = 0; j < i_Energy_TeV_lin.size(); j++ )
        {
            c->cd( j + 1 );
            gPad->SetGridx( 0 );
            gPad->SetGridy( 0 );
            
            TH1D* h = getTheta2orThetaHistogram( i, i_Energy_TeV_lin[j], iPlotTheta2 );
            if( h )
            {
                TH1D* hCumu = get_Cumulative_Histogram( h, true, true );
                if( iCumulative )
                {
                    hCumu->Draw( "same" );
                }
                else
                {
                    // rebin
                    if( iPlotTheta2 )
                    {
                        h->Rebin( 2 );
                    }
                    else
                    {
                        if( h->GetNbinsX() % 4  == 0 )
                        {
                            h->Rebin( 4 );
                        }
                        else
                        {
                            h->Rebin( 5 );
                        }
                    }
                    // normalize
                    if( h->GetMaximum() > 0. )
                    {
                        h->Scale( 1. / h->GetMaximum() );
                    }
                    h->Draw( "same" );
                }
                // get 68% per value
                double x68 = hCumu->GetXaxis()->GetBinCenter( hCumu->FindFirstBinAbove( 0.68 ) );
                TLine* iL68 = new TLine( x68, 0., x68, 1.1 );
                iL68->SetLineStyle( 2 );
                iL68->SetLineColor( hCumu->GetLineColor() );
                iL68->Draw();
                cout << "68% value at " << i_Energy_TeV_lin[j] << " TeV: ";
                cout << sqrt( x68 ) << " deg" << endl;
                // set line at '1' for cumulative histograms
                if( iCumulative )
                {
                    TLine* iL = new TLine( h->GetXaxis()->GetXmin(), 1., iTheta2AxisMax, 1. );
                    iL->SetLineStyle( 2 );
                    iL->Draw();
                }
            }
        }
    }
    return c;
}


TLegend* VPlotInstrumentResponseFunction::makeLegend()
{

    TLegend* iLegend = 0;
    if( fLegendX1 > -99 )
    {
        iLegend = new TLegend( fLegendX1, fLegendY1, fLegendX2, fLegendY2, fLegendHeader.Data(), fLegendOption.Data() );
        if( iLegend )
        {
            iLegend->SetFillColor( kWhite );
            iLegend->SetNColumns( fLegendNCol );
            
            //center header
            TLegendEntry* header = ( TLegendEntry* )iLegend->GetListOfPrimitives()->First();
            if( header )
            {
                header->SetTextAlign( 22 );
            }
        }
    }
    return iLegend;
}

