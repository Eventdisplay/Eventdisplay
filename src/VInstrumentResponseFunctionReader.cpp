/*! \file VInstrumentResponseFunctionReader

    data class for effective area plotting

*/

#include "VInstrumentResponseFunctionReader.h"

VInstrumentResponseFunctionReader::VInstrumentResponseFunctionReader()
{
    fIsZombie = true;
    fDebug = false;
    
    fGammaHadronCuts_directionCut_selector = 0;
    
    fFile = "";
    fA_MC = "A_MC";
    fZe = 0.;
    fWoff = 0.;
    fAzbin = 0;
    fIndex = 0.;
    fNoise = 0;
    setPlotOption();
    fColor = 1;
    fLineStyle = 1;
    fMarkerStyle = 20;
    fLegend = "";
    
    gEffArea_MC = 0;
    gEffArea_MC_Ratio = 0;
    gEffArea_Rec = 0;
    gEffArea_Rec_Ratio = 0;
    gEffArea_MCnoTh2 = 0;
    gEffArea_MCnoTh2_Ratio = 0;
    gEffArea_RecnoTh2 = 0;
    gEffArea_RecnoTh2_Ratio = 0;
    
    hEmc = 0;
    hEcut = 0;
    hEcut_rec = 0;
    hEcutUW = 0;
    hEcut_recUW = 0;
    hEsys = 0;
    hERecMatrix = 0;
    hERecMatrixCoarse = 0;
    hERecMatrixNoDirectionCuts = 0;
    hERecMatrixCoarseNoDirectionCuts = 0;
    hEsysMCRelative = 0;
    hEsysMCRelative2D = 0;
    hEsysMCRelative2DNoDirectionCut = 0;
    gEnergyResolution = 0;
    gEnergyResolutionQC = 0;
    gEnergyResolutionNoDirectionCuts = 0;
    gEnergyBias_Mean = 0;
    gEnergyBias_Median = 0;
    gEnergyLogBias_Mean = 0;
    gEnergyLogBias_Median = 0;
    gAngularResolution = 0;
    gAngularResolution80 = 0;
    gAngularResolution95 = 0;
    h2DAngularPSF = 0;
    h2DAngularPSFEMC = 0;
    gEffArea_Recp80 = 0;
    hWeightedRate = 0;
    hWeightedRate005 = 0;
    
    // methods for energy resolution
    fEnergyResolutionMethod = 1;
    fEnergyXaxisIsEtrue = false;
    fEnergy_TeV_Eres_oneSided = -99.;
    
    setEnergyRange();
    
    initializeIRFData();
}

bool VInstrumentResponseFunctionReader::initializeIRFData()
{
    // read angular and core resolution
    fIRF_TreeNames.push_back( "t_angular_resolution" );
    fIRF_TreeNames.push_back( "t_angular_resolution_080p" );
    fIRF_TreeNames.push_back( "t_angular_resolution_095p" );
    fIRF_TreeNames.push_back( "t_core_resolution" );
    fIRF_TreeNames.push_back( "t_energy_resolution" );
    
    for( unsigned int i = 0; i < fIRF_TreeNames.size(); i++ )
    {
        vector< TH2D* > iTemp;
        fIRF_Data.push_back( 0 );
    }
    
    return true;
}

bool VInstrumentResponseFunctionReader::fillData( string iDataLine, int iDataID )
{
    if( iDataLine.size() == 0 )
    {
        fIsZombie = true;
        return false;
    }
    string temp;
    
    istringstream is_stream( iDataLine );
    is_stream >> temp;
    if( temp != "*" )
    {
        return false;
    }
    is_stream >> temp;
    // check set number
    if( atoi( temp.c_str() ) != iDataID )
    {
        return false;
    }
    // read this line
    is_stream >> fFile;
    is_stream >> temp;
    fZe = atof( temp.c_str() );
    is_stream >> temp;
    fAzbin = atoi( temp.c_str() );
    is_stream >> temp;
    fWoff = atof( temp.c_str() );
    is_stream >> temp;
    fIndex = atof( temp.c_str() );
    is_stream >> temp;
    fNoise = atoi( temp.c_str() );
    is_stream >> fA_MC;
    
    // plotting options
    is_stream >> fPlotOption;
    is_stream >> temp;
    fPlottingColor = atoi( temp.c_str() );
    is_stream >> temp;
    fPlottingLineStyle = atoi( temp.c_str() );
    is_stream >> temp;
    fPlottingMarkerStyle = atoi( temp.c_str() );
    fLegend = is_stream.str();
    fLegend = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
    
    return fillData();
}

bool VInstrumentResponseFunctionReader::fillData( string iFile, double iZe, double iWoff, int iAzBin, double iIndex, int iNoise, string iA_MC )
{
    fFile = iFile;
    fZe = iZe;
    fWoff = iWoff;
    fAzbin = iAzBin;
    fIndex = iIndex;
    fNoise = iNoise;
    fA_MC = iA_MC;
    
    return fillData();
}

bool VInstrumentResponseFunctionReader::fillData()
{
    // read in all the necessary data from the effective area tree
    
    if( !getDataFromFile() )
    {
        fIsZombie = true;
        return false;
    }
    
    // calculate ratio of cut - efficiencies
    if( !calculateCutEfficiencies() )
    {
        return false;
    }
    
    fIsZombie = false;
    
    return true;
}

/*

    read response functions from CTA file

    see http://www.cta-observatory.org/ctawpcwiki/index.php/WP_MC for documentation

*/
bool VInstrumentResponseFunctionReader::getDataFromCTAFile()
{
    bool bLinX = false;  // energy axis for effective areas
    
    TH1F* h = 0;
    // gamma-ray effective area vs reconstruction energy
    h = get_CTA_IRF_Histograms( "EffectiveArea", fWoff );
    if( !h )
    {
        h = get_CTA_IRF_Histograms( "harea_gamma", fWoff );
        if( h )
        {
            bLinX = true;
        }
    }
    if( h )
    {
        gEffArea_Rec = new TGraphAsymmErrors( 1 );
        gEffArea_Rec->SetName( "gEffArea_Rec" );
        setGraphPlottingStyle( gEffArea_Rec );
        get_Graph_from_Histogram( h, gEffArea_Rec, false, bLinX );
    }
    // gamma-ray effective area vs true energy
    h = get_CTA_IRF_Histograms( "EffectiveAreaEtrue", fWoff );
    if( !h )
    {
        h = get_CTA_IRF_Histograms( "harea_gamma", fWoff );
        if( h )
        {
            bLinX = true;
        }
        else
        {
            h = get_CTA_IRF_Histograms( "EffectiveArea", fWoff );
        }
    }
    if( h )
    {
        gEffArea_MC = new TGraphAsymmErrors( 1 );
        gEffArea_MC->SetName( "gEffArea_MC" );
        setGraphPlottingStyle( gEffArea_MC );
        get_Graph_from_Histogram( h, gEffArea_MC, false, bLinX, -1., log10( fEnergyLinTeV_min ), log10( fEnergyLinTeV_max ) );
    }
    else
    {
        gEffArea_MC = 0;
    }
    
    ///////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////
    // name and axis units are not consistent in the CTA files!!!
    ///////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////
    // energy resolution
    h = 0;
    gEnergyResolution = new TGraphErrors( 1 );
    h = ( TH1F* )get_CTA_IRF_Histograms( "ERes", fWoff );
    if( !h )
    {
        h = ( TH1F* )get_CTA_IRF_Histograms( "EnResol_RMS", fWoff );
    }
    if( h )
    {
        get_Graph_from_Histogram( h, gEnergyResolution, true, 0., log10( fEnergyLinTeV_min ), log10( fEnergyLinTeV_max ) );
        setGraphPlottingStyle( gEnergyResolution );
    }
    ///////////////////////////////////////////////////////////////
    // energy bias
    h = 0;
    gEnergyBias_Mean = new TGraphErrors( 1 );
    h = ( TH1F* )get_CTA_IRF_Histograms( "Ebias", fWoff );
    if( !h )
    {
        h = ( TH1F* )get_CTA_IRF_Histograms_from2D( "EestOverEtrue", -1. );
    }
    if( h )
    {
        get_Graph_from_Histogram( h, gEnergyBias_Mean, true, -100., log10( fEnergyLinTeV_min ), log10( fEnergyLinTeV_max ) );
        setGraphPlottingStyle( gEnergyBias_Mean );
    }
    ///////////////////////////////////////////////////////////////
    // angular resolution graphs
    
    gAngularResolution   = getAngularResolutionGraphs( "AngRes", "AngResolution68" );
    gAngularResolution80 = getAngularResolutionGraphs( "AngRes80", "AngResolution80" );
    gAngularResolution95 = getAngularResolutionGraphs( "AngRes95", "AngResolution95" );
    
    return true;
}

/*!

    read IRF from a root file

    might be a

    - evndisp IRF file (produced with makeEffectiveArea)
    - CTA WP MC response file

*/
bool VInstrumentResponseFunctionReader::getDataFromFile()
{
    if( fDebug )
    {
        cout << "VInstrumentResponseFunctionReader::getDataFromFile" << endl;
    }
    
    TFile* iFile = new TFile( fFile.c_str() );
    if( iFile->IsZombie() )
    {
        return false;
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // read effective areas
    //
    TTree* t = ( TTree* )iFile->Get( "fEffArea" );
    // not effective area tree - is this a CTA file?
    if( !t )
    {
        // try to read CTA file
        // (WP Phys style)
        if( !getDataFromCTAFile() )
        {
            cout << "Error: effective area histogram not found in CTA-style file" << endl;
            return false;
        }
        else
        {
            if( iFile )
            {
                iFile->Close();
            }
            return true;
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////
    // read IRFs from a EVNDISP effective area (response) file
    ////////////////////////////////////////////////////////////////////
    CEffArea* c = new CEffArea( t );
    
    bool bFound = false;
    for( int j = 0; j < c->fChain->GetEntries(); j++ )
    {
        c->GetEntry( j );
        
        if( fDebug )
        {
            cout << "VInstrumentResponseFunctionReader::getDataFromFile: reading event " << j << endl;
        }
        
        ////////////////////////////////////////////////////////////////
        // find right entry according to ze, az, woff, etc.
        // ignore all checks if there is only one entry in this tree
        if( c->fChain->GetEntries() > 1 )
        {
            // azimuth
            if( fDebug )
            {
                cout << "\t AZ: tree entry " << j << "\t az " << c->az << "\t az bin " << fAzbin << endl;
            }
            if( c->az != fAzbin )
            {
                continue;
            }
            // spectral index
            if( fDebug )
            {
                cout << "\t Index: " << j << "\t" << c->index << "\t " << fIndex << endl;
            }
            if( TMath::Abs( c->index - fIndex ) > 0.05 )
            {
                continue;
            }
            // wobble offset
            if( fDebug )
            {
                cout << "\t Woff: " << j << "\t" << c->Woff << "\t" << fWoff << endl;
            }
            if( TMath::Abs( c->Woff - fWoff ) > 0.05 )
            {
                continue;
            }
            // noise level
            if( fDebug )
            {
                cout << "\t Noise: " << j << "\t" << c->noise << "\t" << fNoise << endl;
            }
            if( c->noise != fNoise )
            {
                continue;
            }
            // zenith angle
            if( fDebug )
            {
                cout << "\t Ze: " << j << "\t" << c->ze << "\t" << fZe << endl;
            }
            if( TMath::Abs( c->ze - fZe ) > 3. )
            {
                continue;
            }
        }
        cout << "\t FOUND EFFECTIVE AREA (entry " << j << ")" << endl;
        bFound = true;
        
        // get effective areas (these are effective area after the usual angular cut)
        if( c->gEffAreaMC )
        {
            gEffArea_MC  = ( TGraphAsymmErrors* )c->gEffAreaMC->Clone();
            setGraphPlottingStyle( gEffArea_MC );
        }
        else if( c->nbins > 0 )
        {
            gEffArea_MC  = new TGraphAsymmErrors( 1 );
            unsigned int z = 0;
            for( int k = 0; k < c->nbins; k++ )
            {
                if( c->eff[k] > 0. )
                {
                    gEffArea_MC->SetPoint( z, c->e0[k], c->eff[k] );
                    gEffArea_MC->SetPointEYlow( z, c->eff_error[k] );
                    gEffArea_MC->SetPointEYhigh( z, c->eff_error[k] );
                    z++;
                }
            }
            if( gEffArea_MC )
            {
                setGraphPlottingStyle( gEffArea_MC );
            }
        }
        else
        {
            gEffArea_MC = 0;
        }
        if( c->gEffAreaRec )
        {
            gEffArea_Rec = ( TGraphAsymmErrors* )c->gEffAreaRec->Clone();
            if( gEffArea_Rec )
            {
                setGraphPlottingStyle( gEffArea_Rec );
            }
        }
        else if( c->nbins > 0 )
        {
            gEffArea_Rec  = new TGraphAsymmErrors( 1 );
            unsigned int z = 0;
            for( int k = 0; k < c->nbins; k++ )
            {
                if( c->Rec_eff[k] > 0. )
                {
                    gEffArea_Rec->SetPoint( z, c->e0[k], c->Rec_eff[k] );
                    gEffArea_Rec->SetPointEYlow( z, c->Rec_eff_error[k] );
                    gEffArea_Rec->SetPointEYhigh( z, c->Rec_eff_error[k] );
                    z++;
                }
            }
            setGraphPlottingStyle( gEffArea_Rec );
        }
        else
        {
            gEffArea_Rec = 0;
        }
        // effective areas without direction cut
        if( c->gEffAreaNoTh2MC )
        {
            gEffArea_MCnoTh2 = c->gEffAreaNoTh2MC;
            if( gEffArea_MCnoTh2 )
            {
                setGraphPlottingStyle( gEffArea_MCnoTh2 );
            }
        }
        else
        {
            gEffArea_MCnoTh2 = 0;
        }
        if( c->gEffAreaNoTh2Rec )
        {
            gEffArea_RecnoTh2 = c->gEffAreaNoTh2Rec;
            if( gEffArea_RecnoTh2 )
            {
                setGraphPlottingStyle( gEffArea_RecnoTh2 );
            }
        }
        else
        {
            gEffArea_RecnoTh2 = 0;
        }
        // get energy counting histogram
        if( c->hEmc )
        {
            hEmc = ( TH1D* )c->hEmc->Clone();
            setHistogramPlottingStyle( hEmc );
        }
        if( c->hEcut )
        {
            hEcut = ( TH1D* )c->hEcut->Clone();
            setHistogramPlottingStyle( hEcut );
        }
        if( c->hEcutRec )
        {
            hEcut_rec = ( TH1D* )c->hEcutRec->Clone();
            setHistogramPlottingStyle( hEcut_rec );
            hEcut_rec->SetMarkerStyle( hEcut->GetMarkerStyle() + 4 );
        }
        if( c->hEcutUW )
        {
            hEcutUW = ( TH1D* )c->hEcutUW->Clone();
            setHistogramPlottingStyle( hEcutUW );
        }
        if( c->hEcutRecUW )
        {
            hEcut_recUW = ( TH1D* )c->hEcutRecUW->Clone();
            setHistogramPlottingStyle( hEcut_recUW );
            hEcut_recUW->SetMarkerStyle( hEcutUW->GetMarkerStyle() + 4 );
        }
        // get energy reconstruction matrix
        hERecMatrix = ( TH2D* )c->hResponseMatrixFine;
        hERecMatrixCoarse = ( TH2D* )c->hResponseMatrix;
        hERecMatrixQC = ( TH2D* )c->hResponseMatrixFineQC;
        hERecMatrixCoarseQC = ( TH2D* )c->hResponseMatrixQC;
        hERecMatrixNoDirectionCuts = ( TH2D* )c->hResponseMatrixFineNoDirectionCuts;
        hERecMatrixCoarseNoDirectionCuts = ( TH2D* )c->hResponseMatrixNoDirectionCuts;
        // get error in energy reconstruction
        hEsys = ( TH2D* )c->hEsys2D;
        // erec/emc
        hEsysMCRelative = ( TProfile* )c->hEsysMCRelative;
        hEsysMCRelative2D = ( TH2D* )c->hEsysMCRelative2D;
        hEsysMCRelative2DNoDirectionCut = ( TH2D* )c->hEsysMCRelative2DNoDirectionCut;
        // get energy resolution (!!)
        //       getEnergyResolutionPlot( (TProfile*)c->hEsysMCRelative );
        // energy resolution calculation as 68% value
        //       getEnergyResolutionPlot68( (TH2D*)c->hEsysMCRelative2D );
        if( fEnergyResolutionMethod == 0 )
        {
            // energy resolution is RMS
            // (this is non standard)
            getEnergyResolutionPlot( ( TH2D* )c->hEsysMCRelativeRMS );
        }
        // Fermi LAT method (default)
        else if( fEnergyResolutionMethod == 1 )
        {
            gEnergyResolution = getEnergyResolutionMPropInterval( ( TH2D* )c->hResponseMatrixFine, fEnergyXaxisIsEtrue, gEnergyResolution );
            gEnergyResolutionQC = getEnergyResolutionMPropInterval( ( TH2D* )c->hResponseMatrixFineQC, fEnergyXaxisIsEtrue, gEnergyResolutionQC );
            gEnergyResolutionNoDirectionCuts = getEnergyResolutionMPropInterval( ( TH2D* )c->hResponseMatrixFineNoDirectionCuts, fEnergyXaxisIsEtrue, gEnergyResolutionNoDirectionCuts );
        }
        setGraphPlottingStyle( gEnergyResolution );
        setGraphPlottingStyle( gEnergyResolutionQC );
        setGraphPlottingStyle( gEnergyResolutionNoDirectionCuts );
        // get energy bias
        gEnergyBias_Mean = get_Profile_from_TH2D( ( TH2D* )c->hEsysMCRelativeRMS, 0, "mean", 1, -10., -1. );
        setGraphPlottingStyle( gEnergyBias_Mean );
        gEnergyBias_Median = get_Profile_from_TH2D( ( TH2D* )c->hEsysMCRelativeRMS, 0, "median", 1, -10., -1. );
        setGraphPlottingStyle( gEnergyBias_Median, 1, 1., 7 );
        gEnergyLogBias_Mean = get_Profile_from_TH2D( ( TH2D* )c->hEsys2D, 0, "mean", 1, -10. );
        setGraphPlottingStyle( gEnergyLogBias_Mean, 1, 1., 7 );
        gEnergyLogBias_Median = get_Profile_from_TH2D( ( TH2D* )c->hEsys2D, 0, "median", 1, -10. );
        setGraphPlottingStyle( gEnergyLogBias_Median );
        // get rate histograms
        hWeightedRate = ( TH1D* )c->hWeightedRate;
        hWeightedRate005 = ( TH1D* )c->hWeightedRate005;
        // get cut efficiencies
        if( c->hhEcutTrigger )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutTrigger->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        if( c->hhEcutFiducialArea )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutFiducialArea->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        if( c->hhEcutStereoQuality )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutStereoQuality->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        if( c->hhEcutTelType )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutTelType->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        if( c->hhEcutDirection )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutDirection->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        if( c->hhEcutEnergyReconstruction )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutEnergyReconstruction->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        if( c->hhEcutGammaHadron )
        {
            hCutEfficiency.push_back( ( TH1D* )c->hhEcutGammaHadron->Clone() );
        }
        else
        {
            hCutEfficiency.push_back( 0 );
        }
        for( unsigned int i = 0; i < hCutEfficiency.size(); i++ )
        {
            setHistogramPlottingStyle( hCutEfficiency[i], i + 1, 1., 1.5, fPlottingMarkerStyle );
        }
        
        break;
    }
    //////////////////////////////////////////////////////////////
    // read resolution files
    
    for( unsigned int i = 0; i < fIRF_TreeNames.size(); i++ )
    {
        cout << "reading IRF " << fIRF_TreeNames[i].c_str() << endl;
        TTree* t_a = ( TTree* )iFile->Get( fIRF_TreeNames[i].c_str() );
        fIRF_Data[i] = getIRFFromFile( t_a );
        if( !fIRF_Data[i] )
        {
            cout << " ...not found (index " << i << ")" << endl;
        }
        else
        {
            cout << " ...found" << endl;
            bFound = true;
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // now fill effective areas for 80% containment radius
    //
    // for this: read theta2 cut, get containment radius for this theta2 cut,
    //           then scale the effective areas accordingly
    
    // read gamma/hadron cuts from disk
    VGammaHadronCuts* i_cuts = ( VGammaHadronCuts* )iFile->Get( "GammaHadronCuts" );
    if( i_cuts && gEffArea_Rec )
    {
        // check if theta2 graph was used
        fGammaHadronCuts_directionCut_selector = i_cuts->getDirectionCutSelector();
        
        /////////////////////////////////
        // read theta2 from disk
        
        double i_x = 0;
        double i_y = 0.;
        int z = 0;
        for( unsigned int i = 0; i < fIRF_TreeNames.size(); i++ )
        {
            if( fIRF_TreeNames[i] == "t_angular_resolution" && fIRF_Data[i] )
            {
                ///////////////////////////////////////
                // 2D angular difference histogram
                // (difference between reconstruction and true direction vs energy)
                TH2D* hRes = ( TH2D* )fIRF_Data[i]->f2DHisto[VInstrumentResponseFunctionData::E_DIFF];
                if( hRes )
                {
                    z = 0;
                    gEffArea_Recp80 = new TGraphAsymmErrors( 1 );
                    setGraphPlottingStyle( gEffArea_Recp80, 2 );
                    // loop over all energies (x-axis of hRes) and get containment level for given theta2 cut
                    float i_ffactor_80p = 1.;
                    for( int p = 0; p < gEffArea_Rec->GetN(); p++ )
                    {
                        gEffArea_Rec->GetPoint( p, i_x, i_y );
                        i_ffactor_80p = i_cuts->getTheta2Cut_max( TMath::Power( 10., i_x ) );
                        if( i_ffactor_80p > 0. )
                        {
                            i_ffactor_80p = sqrt( i_ffactor_80p );
                        }
                        // integrate and calculate contaiment level for this bin in log10 energy
                        int i_bin = hRes->GetXaxis()->FindBin( i_x );
                        double iTot = 0.;
                        for( int j = 1; j <= hRes->GetNbinsY(); j++ )
                        {
                            iTot += hRes->GetBinContent( i_bin, j );
                        }
                        if( iTot < 1.e-12 )
                        {
                            continue;
                        }
                        double iSum = 0.;
                        for( int j = 1; j <= hRes->GetNbinsY(); j++ )
                        {
                            iSum += hRes->GetBinContent( i_bin, j );
                            if( hRes->GetYaxis()->GetBinCenter( j ) > i_ffactor_80p )
                            {
                                if( iSum / iTot > 0. )
                                {
                                    gEffArea_Recp80->SetPoint( z, i_x, i_y * 0.8 / iSum * iTot );
                                    gEffArea_Recp80->SetPointEYhigh( z, gEffArea_Rec->GetErrorYhigh( p )  * 0.8 / iSum * iTot );
                                    gEffArea_Recp80->SetPointEYlow( z, gEffArea_Rec->GetErrorYlow( p )  * 0.8 / iSum * iTot );
                                    z++;
                                }
                                break;
                            }
                        }
                    }
                }
                break; // found the angular resolution tree
            }
        }
    }
    else
    {
        fGammaHadronCuts_directionCut_selector = 0;
    }
    
    if( !bFound )
    {
        return false;
    }
    
    return true;
}



VInstrumentResponseFunctionData* VInstrumentResponseFunctionReader::getIRFFromFile( TTree* t )
{
    if( !t )
    {
        return 0;
    }
    
    VInstrumentResponseFunctionData* c = new VInstrumentResponseFunctionData();
    TBranch* br = t->GetBranch( "IRF" );
    br->SetAddress( &c );
    
    for( int j = 0; j < t->GetEntries(); j++ )
    {
        t->GetEntry( j );
        
        if( fDebug )
        {
            cout << "VInstrumentResponseFunctionReader::getDataFromFile (resolution data): reading event " << j << endl;
        }
        
        // check that there is data for this tree entry
        if( !c )
        {
            continue;
        }
        
        // ignore all values if there is only one entry in this tree
        if( t->GetEntries() > 1 )
        {
            // azimuth
            if( fDebug )
            {
                cout << "IRF AZ: " << j << ", found: " << c->fAz_bin << ", searched for: " << fAzbin << endl;
            }
            if( c->fAz_bin != fAzbin )
            {
                continue;
            }
            // spectral index
            if( fDebug )
            {
                cout << "IRF Index: " << j << ", found: " << c->fSpectralIndex << ", searched for: " << fIndex << endl;
            }
            if( TMath::Abs( c->fSpectralIndex - fIndex ) > 0.05 )
            {
                continue;
            }
            // wobble offset
            if( fDebug )
            {
                cout << "IRF Woff: " << j << ", found: " << c->fWobble << ", searched for: " << fWoff << endl;
            }
            if( TMath::Abs( c->fWobble - fWoff ) > 0.05 )
            {
                continue;
            }
            // noise level
            if( fDebug )
            {
                cout << "IRF Noise: " << j << ", found: " << c->fNoise << ", searched for: " << fNoise << endl;
            }
            if( c->fNoise != fNoise )
            {
                continue;
            }
            // zenith angle
            if( fDebug )
            {
                cout << "IRF Ze: " << j << ", found: " << c->fZe << ", searched for: " << fZe << endl;
            }
            if( TMath::Abs( c->fZe - fZe ) > 3. )
            {
                continue;
            }
        }
        if( c && c->fResolutionGraph.size() > 0 )
        {
            for( unsigned int r = 0; r < c->fResolutionGraph.size(); r++ )
            {
                if( c->fResolutionGraph[r] )
                {
                    setGraphPlottingStyle( c->fResolutionGraph[r] );
                }
            }
        }
        
        return ( VInstrumentResponseFunctionData* )c->Clone();
        
    }
    
    return 0;
}

TGraphErrors* VInstrumentResponseFunctionReader::getResolutionGraph( TGraphErrors* g )
{
    g = new TGraphErrors( 1 );
    g->SetMarkerStyle( 20 );
    g->SetMarkerSize( 2 );
    g->SetLineWidth( 2 );
    g->SetTitle( "" );
    g->SetName( "" );
    setGraphPlottingStyle( g );
    
    return g;
}

/*

      calculate energy resolution as the 68% interval around the
      most probably or median value for the energy reconstruction

      adapted from Fermi LAT, see Ackermann et al 2012, ApJS 203:4, Figure 67 (p. 54)

*/
TGraphErrors* VInstrumentResponseFunctionReader::getEnergyResolutionMPropInterval( TH2D* migmatrix,
        bool bXaxisIsEtrue,
        TGraphErrors* iEnergyResolution )
{
    ////////////////////////////////////
    // hard wired parameters
    
    // half width of 68%
    double fMProp_fContainment = 0.68 / 2.;
    // bootstrap to calculate errors on energy resolution
    // (don't bootstart if nBootStrap = 1 )
    unsigned int fMProp_nBootStrap = 100;
    // minimum number of events required per distribution
    // to be included in graph
    int fMProp_minEvents = 10;
    // maximum RMS to value ratio to be included in graph
    double fMProp_maxerror = 0.3;
    // use median or most probable?
    // (median is more stable; even for low statistics bin)
    bool bUseMedian = true;   // median is the default for e.g. prod3
    // energy resolution is averaged over these number of bins in the
    // migration matrix
    int i_zoffset_bins = 20;
    
    if( bUseMedian )
    {
        cout << "VInstrumentResponseFunctionReader::getEnergyResolutionMPropInterval:";
        cout << "use median of energy distributions" << endl;
    }
    else
    {
        cout << "VInstrumentResponseFunctionReader::getEnergyResolutionMPropInterval:";
        cout << "use mop of energy distributions" << endl;
    }
    cout << "\t minevents: " << fMProp_minEvents;
    cout << ", maxErrorRatio: " << fMProp_maxerror;
    cout << ", bootstrap loops: " << fMProp_nBootStrap;
    cout << endl;
    
    if( !migmatrix )
    {
        return iEnergyResolution;
    }
    // get a (empty) graph
    iEnergyResolution = getResolutionGraph( iEnergyResolution );
    if( !iEnergyResolution )
    {
        return iEnergyResolution;
    }
    
    char hname[200];
    
    // temporary variables
    int z = 0;
    int i_z = 0;
    float eres = 0.;
    float eres_temp = 0.;
    float eres_tempRMS = 0.;
    double energy = 0.;
    double energy_mop_med = 0.;
    
    TH1D i_h_zoffsetsHisto( "i_h_zoffsetsHisto", "", 1000., 0., 1. );
    
    int nbins = 0;
    bXaxisIsEtrue = true;
    bXaxisIsEtrue = false;
    if( bXaxisIsEtrue )
    {
        nbins = migmatrix->GetNbinsY();
        cout << "\t Migration matrix: number of bins in true energy: " << nbins << endl;
    }
    else
    {
        nbins = migmatrix->GetNbinsX();
        cout << "\t Migration matrix: number of bins in reconstructed energy: " << nbins << endl;
    }
    // MigMatrix is not fine binned: reduce averaging to 2
    if( nbins < 100 )
    {
        i_zoffset_bins = 2;
    }
    
    ///////////////////////////////////////////////////////////////////////
    // loop over axis in true energy of migration matrix
    for( int i = 2; i < nbins; i++ )
    {
        TH1F* hDist = 0;
        // fixed E_true, energy resolution from E_rec distribution
        if( bXaxisIsEtrue )
        {
            hDist = ( TH1F* )migmatrix->ProjectionX( "hx", i, i );
            energy = migmatrix->GetYaxis()->GetBinCenter( i );
        }
        // fixed E_rec, energy resolution from E_true distribution
        else
        {
            hDist = ( TH1F* )migmatrix->ProjectionY( "hy", i, i );
            energy = migmatrix->GetXaxis()->GetBinCenter( i );
        }

        // require at least a few entries
        if( hDist->GetEntries() < fMProp_minEvents )
        {
            continue;
        }
        
        ///////////////////////////////////
        // Fermi LAT method
        // calculates half-width of +-34% interval around
        // most probably or median energy
        
        // (normalized) cumulative distribution
        sprintf( hname, "%s_CUMU", hDist->GetName() );
        TH1F* hcum = ( TH1F* )hDist->Clone( hname );
        
        // bootstrap to calculate errors on energy resolution
        // (don't bootstart if nBootStrap = 1 )
        for( unsigned int b = 0; b < fMProp_nBootStrap; b++ )
        {
            // bootstrap only if more than one cycle is
            // required
            if( fMProp_nBootStrap > 1 )
            {
                hcum->Reset();
                hcum->FillRandom( hDist, hDist->GetEntries() );
            }
            
            // cumulative distribution
            for( int j = 2; j <= hcum->GetNbinsX(); j++ )
            {
                hcum->SetBinContent( j, hcum->GetBinContent( j ) + hcum->GetBinContent( j - 1 ) );
            }
            if( hcum->GetMaximum() > 0. )
            {
                hcum->Scale( 1. / hcum->GetMaximum() );
            }
            
            // most probable energy
            double energy_mprob = hDist->GetXaxis()->GetBinCenter( hDist->GetMaximumBin() );
            // median energy
            double energy_median = 0.5 * ( hcum->GetXaxis()->GetBinCenter( hcum->FindFirstBinAbove( 0.5 ) )
                                           + hcum->GetXaxis()->GetBinCenter( hcum->FindFirstBinAbove( 0.5 ) - 1 ) );
            // use median or mprob?
            // (median is more stable for small statistic
            //  bins)
            if( bUseMedian )
            {
                energy_mop_med = energy_median;
            }
            else
            {
                energy_mop_med = energy_mprob;
            }
            
            // determine energy resolution:
            // find 34% interval around energy value
            float prob_mprob = hcum->GetBinContent( hcum->GetXaxis()->FindBin( energy_mop_med ) );
            if( bUseMedian )
            {
                prob_mprob = 0.5;
            }
            
            // energy of med -34%
            float prob_mprob_minus = prob_mprob - fMProp_fContainment;
            if( prob_mprob_minus < 0. )
            {
                prob_mprob_minus = 0.;
            }
            int prob_minus = hcum->FindFirstBinAbove( prob_mprob_minus );
            double energy_prob_minus = hcum->GetXaxis()->GetBinCenter( prob_minus );
            if( prob_minus > 2 )
            {
                energy_prob_minus = 0.5 * ( hcum->GetXaxis()->GetBinCenter( prob_minus ) + hcum->GetXaxis()->GetBinCenter( prob_minus - 1 ) );
            }
            // energy of med +34%
            float prob_mprob_plus  = prob_mprob + fMProp_fContainment;
            if( prob_mprob_plus > 1. )
            {
                prob_mprob_plus = 1. - 1.e-3;
            }
            int prob_plus  = hcum->FindFirstBinAbove( prob_mprob_plus ); // TMPTMP -1 ?
            double energy_prob_plus = hcum->GetXaxis()->GetBinCenter( prob_plus );
            if( prob_plus > 2 )
            {
                energy_prob_plus = 0.5 * ( hcum->GetXaxis()->GetBinCenter( prob_plus ) + hcum->GetXaxis()->GetBinCenter( prob_plus - 1 ) );
            }
            // energy resolution
            eres = TMath::Power( 10., energy_prob_plus ) - TMath::Power( 10., energy_prob_minus );
            eres /= TMath::Power( 10., energy_mop_med );
            eres *= 0.5;  // half width
            // use one sided interval in the threshold area
            if( fEnergy_TeV_Eres_oneSided > 0. && energy < log10( fEnergy_TeV_Eres_oneSided )
                    && eres > 1.e-4 )
            {
                cout << "Double sided log10(Erec) at " << energy << " TeV: " << eres * 100. << "%" << endl;
                eres = TMath::Power( 10., hcum->GetXaxis()->GetBinCenter( prob_plus ) ) - TMath::Power( 10., energy_mop_med );
                eres /= TMath::Power( 10., energy_mop_med );
                cout << "One sided Erec at " << energy << " TeV: " << eres * 100. << "%" << endl;
            }
            
            i_h_zoffsetsHisto.Fill( eres );
        }
        i_z++;
        
        // average of i_zoffset_bins bins
        if( i_z == i_zoffset_bins )
        {
            double energy_graph_midpoint = migmatrix->GetXaxis()->GetBinCenter( i - 0.5 * i_zoffset_bins );
            eres_temp    = i_h_zoffsetsHisto.GetMean();
            eres_tempRMS = i_h_zoffsetsHisto.GetRMS();
            
            if( eres_temp > 0 && eres_tempRMS / eres_temp < fMProp_maxerror )
            {
                iEnergyResolution->SetPoint( z, energy_graph_midpoint, eres_temp );
                iEnergyResolution->SetPointError( z, 0., eres_tempRMS );
                z++;
            }
            
            i_h_zoffsetsHisto.Reset();
            i_z = 0;
        }
        if( hcum )
        {
            delete hcum;
        }
    }
    
    return iEnergyResolution;
}

/*
 * energy resolution calculation
 *
 * old and outdated method - do not use
 *
 */
void VInstrumentResponseFunctionReader::getEnergyResolutionPlot68( TH2D* iP, double iReferenceValue )
{
    if( !iP )
    {
        gEnergyResolution = 0;
        return;
    }
    gEnergyResolution = getResolutionGraph( gEnergyResolution );
    
    
    int zz = 0;
    double e_res = 0.;
    for( int b = 1; b <= iP->GetNbinsX(); b++ )
    {
        TH1D* h = iP->ProjectionY( "p_x", b, b + 1 );
        if( h && h->GetEntries() > 3. )
        {
            // calculate quantiles
            double xq[3];
            double yq[3];
            /*	    xq[0] = 0.5-0.6826895/2.;
            	    xq[1] = 0.5;
            	    xq[2] = 0.5+0.6826895/2.;
            	    h->GetQuantiles( 3, yq, xq );
                        if( iP->GetXaxis()->GetBinCenter( b ) < iMinEnergy ) continue;
            // +-1 sigma around median
            	    e_res = (yq[2]-yq[0])*0.5; */
            // 68% distribution around 1 (bb_ref, expected value)
            TH1D hh( "h", "", h->GetNbinsX(), 0., h->GetXaxis()->GetXmax() - 1. );
            double bb_ref = iReferenceValue;
            // < -998: relative to mean
            if( iReferenceValue < -998. )
            {
                bb_ref = h->GetMean();
            }
            // >  998: relative to median
            else if( iReferenceValue > 998. )
            {
                xq[0] = 0.50;
                h->GetQuantiles( 1, yq, xq );
                bb_ref = yq[0];
            }
            // fill 1D histogram before integration
            for( int bb = 1; bb <= h->GetNbinsX(); bb++ )
            {
                if( h->GetBinCenter( bb ) < bb_ref )
                {
                    hh.Fill( bb_ref - h->GetBinCenter( bb ), h->GetBinContent( bb ) );
                }
                else
                {
                    hh.Fill( h->GetBinCenter( bb ) - bb_ref, h->GetBinContent( bb ) );
                }
            }
            xq[0] = 0.68;
            hh.GetQuantiles( 1, yq, xq );
            e_res = yq[0];
            gEnergyResolution->SetPoint( zz, iP->GetXaxis()->GetBinCenter( b ), e_res );
            if( h->GetEntries() > 1. )
            {
                gEnergyResolution->SetPointError( zz, 0., h->GetRMS() / sqrt( h->GetEntries() - 1. ) );
            }
            else
            {
                gEnergyResolution->SetPointError( zz, 0., 0. );
            }
            zz++;
        }
    }
    return;
}

void VInstrumentResponseFunctionReader::getEnergyResolutionPlot( TH2D* iP, double iMinEnergy )
{
    if( !iP )
    {
        gEnergyResolution = 0;
        return;
    }
    
    gEnergyResolution = new TGraphErrors( 1 );
    gEnergyResolution->SetMarkerStyle( 20 );
    gEnergyResolution->SetMarkerSize( 2 );
    gEnergyResolution->SetLineWidth( 2 );
    gEnergyResolution->SetTitle( "" );
    gEnergyResolution->SetName( "" );
    setGraphPlottingStyle( gEnergyResolution );
    
    int zz = 0;
    for( int b = 1; b <= iP->GetNbinsX(); b++ )
    {
        TH1D* h = iP->ProjectionY( "p_x", b, b );
        if( h && h->GetEntries() > 10. )
        {
            if( iP->GetXaxis()->GetBinCenter( b ) < iMinEnergy )
            {
                continue;
            }
            gEnergyResolution->SetPoint( zz, iP->GetXaxis()->GetBinCenter( b ), h->GetRMS() );
            gEnergyResolution->SetPointError( zz, 0., h->GetRMS() / sqrt( h->GetEntries() - 1. ) );
            zz++;
        }
    }
    return;
}


void VInstrumentResponseFunctionReader::getEnergyResolutionPlot( TProfile* iP, int i_rebin, double iMinEnergy )
{
    if( !iP )
    {
        gEnergyResolution = 0;
        return;
    }
    
    iP->Rebin( i_rebin );
    
    gEnergyResolution = new TGraphErrors( 1 );
    gEnergyResolution->SetMarkerStyle( 20 );
    gEnergyResolution->SetMarkerSize( 2 );
    gEnergyResolution->SetLineWidth( 2 );
    gEnergyResolution->SetTitle( "" );
    gEnergyResolution->SetName( "" );
    setGraphPlottingStyle( gEnergyResolution );
    
    string iErrorOption = iP->GetErrorOption();
    
    int zz = 0;
    for( int b = 1; b <= iP->GetNbinsX(); b++ )
    {
        if( iP->GetBinEntries( b ) > 3. )
        {
            if( iP->GetXaxis()->GetBinCenter( b ) < iMinEnergy )
            {
                continue;
            }
            if( iErrorOption == "s" )
            {
                gEnergyResolution->SetPoint( zz, iP->GetXaxis()->GetBinCenter( b ), iP->GetBinError( b ) );
                if( iP->GetBinEntries( b ) > 0. )
                {
                    gEnergyResolution->SetPointError( zz, 0., iP->GetBinError( b ) / sqrt( iP->GetBinEntries( b ) - 1. ) );
                }
                else
                {
                    gEnergyResolution->SetPointError( zz, 0., 0. );
                }
            }
            else
            {
                gEnergyResolution->SetPoint( zz, iP->GetXaxis()->GetBinCenter( b ), iP->GetBinError( b )*sqrt( iP->GetBinEntries( b ) - 1. ) );
                gEnergyResolution->SetPointError( zz, 0., iP->GetBinError( b ) );
            }
            zz++;
        }
    }
    return;
}


bool VInstrumentResponseFunctionReader::calculateCutEfficiencies()
{
    char hname[200];
    // create histograms
    for( unsigned int i = 0; i < hCutEfficiency.size(); i++ )
    {
        if( hCutEfficiency[i] )
        {
            sprintf( hname, "%s_R", hCutEfficiency[i]->GetName() );
            hCutEfficiencyRelativePlots.push_back( ( TH1D* )hCutEfficiency[i]->Clone( hname ) );
            hCutEfficiencyRelativePlots.back()->Reset();
        }
        else
        {
            hCutEfficiencyRelativePlots.push_back( 0 );
        }
    }
    // calculate relative plots
    for( unsigned int i = 0; i < hCutEfficiencyRelativePlots.size(); i++ )
    {
        if( hCutEfficiencyRelativePlots[i] )
        {
            hCutEfficiencyRelativePlots[i]->Divide( hCutEfficiency[i], hCutEfficiency[0] );
        }
    }
    
    return true;
}

bool VInstrumentResponseFunctionReader::calculateEffectiveAreaRatios( TGraphAsymmErrors* g0 )
{
    if( !g0 )
    {
        return false;
    }
    
    gEffArea_MC_Ratio       = calculateEffectiveAreaRatios( g0, gEffArea_MC );
    gEffArea_Rec_Ratio      = calculateEffectiveAreaRatios( g0, gEffArea_Rec );
    gEffArea_MCnoTh2_Ratio  = calculateEffectiveAreaRatios( g0, gEffArea_MCnoTh2 );
    gEffArea_RecnoTh2_Ratio = calculateEffectiveAreaRatios( g0, gEffArea_RecnoTh2 );
    
    return true;
}

TGraphAsymmErrors* VInstrumentResponseFunctionReader::calculateEffectiveAreaRatios( TGraphAsymmErrors* g0, TGraphAsymmErrors* g1 )
{
    if( !g0 || !g1 )
    {
        return 0;
    }
    
    TGraphAsymmErrors* g = new TGraphAsymmErrors( 1 );
    
    double e = 0.;
    int z = 0;
    
    double x0 = 0.;
    double y0 = 0.;
    double y0_U = 0.;
    double y0_L = 0.;
    
    double x1 = 0.;
    double y1 = 0.;
    double y1_U = 0.;
    double y1_L = 0.;
    
    for( int f = 0; f < g1->GetN(); f++ )
    {
        g1->GetPoint( f, x1, y1 );
        y1_U = g1->GetErrorYhigh( f );
        y1_L = g1->GetErrorYlow( f );
        
        for( int k = 0; k < g0->GetN(); k++ )
        {
            // needed for correct error bars
            g0->GetPoint( k, x0, y0 );
            //x0 = x1;
            //y0 = g0->Eval( x0 );
            y0_U = g0->GetErrorYhigh( k );
            y0_L = g0->GetErrorYlow( k );
            
            if( y0 > 0. )
            {
                if( TMath::Abs( x0 - x1 ) < 0.001 )
                {
                    g->SetPoint( z, x0, y1 / y0 );
                    
                    // note: standard error propagation not correct in this case
                    e = y1_U * y1_U / y0 / y0 + y1 * y1 / y0 / y0 / y0 / y0 * y0_U * y0_U;
                    g->SetPointEYhigh( z, sqrt( e ) );
                    // (Preli)
                    //g->SetPointEYhigh( z, 0. );
                    
                    e = y1_L * y1_L / y0 / y0 + y1 * y1 / y0 / y0 / y0 / y0 * y0_L * y0_L;
                    g->SetPointEYlow( z, sqrt( e ) );
                    // (Preli)
                    //g->SetPointEYlow( z, 0. );
                    
                    z++;
                }
            }
        }
    }
    setGraphPlottingStyle( g );
    
    return g;
}

/*
 *
 * copy effective area values into TH1F histograms
 *
 */

bool VInstrumentResponseFunctionReader::fillEffectiveAreasHistograms( TH1F* hEffRec,
        string iContainmentRadius,
        TH1F* hEff_MC )
{
    if( !hEffRec )
    {
        return false;
    }
    TGraphAsymmErrors* g = gEffArea_Rec;
    
    // make sure that 80% containment radius effective areas are available
    if( iContainmentRadius.size() > 0 && iContainmentRadius == "80" )
    {
        if( gEffArea_Recp80 )
        {
            cout << "VInstrumentResponseFunctionReader::fillEffectiveAreasHistograms: ";
            cout << "found scaled effective areas for 80\% case " << endl;
            g = gEffArea_Recp80;
        }
        else
        {
            return false;
        }
    }
    else if( iContainmentRadius.size() > 0 && iContainmentRadius == "noTH2Cut" )
    {
        if( gEffArea_RecnoTh2 )
        {
            cout << "VInstrumentResponseFunctionReader::fillEffectiveAreasHistograms: ";
            cout << "found effective area without theta2 cut " << endl;
            g = gEffArea_RecnoTh2;
        }
    }
    // fill histogram
    if( g )
    {
        VHistogramUtilities::fill_Graph_in_Histogram( g, hEffRec, false );
    }
    else
    {
        cout << "VInstrumentResponseFunctionReader::fillEffectiveAreasHistograms() warning: " << endl;
        cout << "\t no effective are graph found" << endl;
        return false;
    }
    if( gEffArea_MC && hEff_MC )
    {
        g = 0;
        if( iContainmentRadius.size() > 0 && iContainmentRadius == "noTH2Cut" && gEffArea_MCnoTh2 )
        {
            g = gEffArea_MCnoTh2;
        }
        else if( gEffArea_MC )
        {
            g = gEffArea_MC;
        }
        if( g )
        {
            VHistogramUtilities::fill_Graph_in_Histogram( g, hEff_MC, false );
        }
    }
    
    return true;
}

bool VInstrumentResponseFunctionReader::fillBiasHistograms( TH1F* h, string iMeanOrMedian )
{
    if( !h )
    {
        return false;
    }
    
    TGraphErrors* g = 0;
    if( iMeanOrMedian == "median" )
    {
        g = gEnergyBias_Median;
    }
    else if( iMeanOrMedian == "mean" )
    {
        g = gEnergyBias_Mean;
    }
    else
    {
        cout << "VInstrumentResponseFunctionReader::fillBiasHistograms warning: unknown string: " << iMeanOrMedian << endl;
        cout << "\t allowed values: mean or median" << endl;
        return false;
    }
    if( !g )
    {
        cout << "VInstrumentResponseFunctionReader::fillBiasHistograms warning: no bias graph found" << endl;
        return false;
    }
    if( g->GetN() < 2 )
    {
        cout << "VInstrumentResponseFunctionReader::fillBiasHistograms warning: bias graph with no points" << endl;
        return false;
    }
    // reset histogram binning
    double x_axis[g->GetN() + 1];
    // Obs: assume fixed binning:
    double i_binWidth = 0.5 * ( g->GetX()[1] - g->GetX()[0] );
    for( int i = 0; i < g->GetN(); i++ )
    {
        x_axis[i] = g->GetX()[i] - i_binWidth;
    }
    x_axis[g->GetN()] = g->GetX()[g->GetN() - 1] + i_binWidth;
    h->SetBins( g->GetN(), x_axis );
    // fill histogram
    double x = 0.;
    double y = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        g->GetPoint( i, x, y );
        h->SetBinContent( i + 1, y );
        h->SetBinError( i + 1, g->GetErrorY( i ) );
    };
    
    return true;
}

/*
 * read angular resolution graph from file
 *
*/
TGraphErrors* VInstrumentResponseFunctionReader::getAngularResolutionGraphs( string iHName, string iHNameOldStyle )
{
    TGraphErrors* i_g = new TGraphErrors( 1 );
    TH1F* h = ( TH1F* )get_CTA_IRF_Histograms( iHName.c_str(), fWoff );
    // try Paris style file
    if( !h )
    {
        h = ( TH1F* )get_CTA_IRF_Histograms( iHNameOldStyle.c_str(), fWoff );
    }
    if( h )
    {
        get_Graph_from_Histogram( h, i_g, true, 0., log10( fEnergyLinTeV_min ), log10( fEnergyLinTeV_max ) );   // ignore errors in resolution graph
        setGraphPlottingStyle( i_g );
    }
    return i_g;
}


/*!

    calculate threshold (usually 68%) reconstruction accuracy from 2D histogram

    iHistogram is a 2D histogram (e.g. angular difference vs energy)

    errors are calculated by bootstrapping

*/
TGraphErrors*  VInstrumentResponseFunctionReader::getAngularResolution( TH2D* iHistogram,
        double iContainmentProbability,
        int iMinRequiredEvents )
{
    ////////////////////////////////////
    // hard wired parameters
    
    // bootstrap to calculate errors on energy resolution
    // (don't bootstart if nBootStrap = 1 )
    unsigned int fMProp_nBootStrap = 100;
    
    // require that the number of events does not change
    // significantly between neighbouring bins
    double iMinBinFraction = 0.05;
    
    if( !iHistogram )
    {
        cout << "VInstrumentResponseFunctionReader::getAngularResolution() error: no 2D histogram given" << endl;
        return 0;
    }
    
    // get an empty graph
    TGraphErrors* g = 0;
    g = getResolutionGraph( g );
    if( !g )
    {
        cout << "VInstrumentResponseFunctionReader::getAngularResolution() error: failure defining graph" << endl;
        return 0;
    }
    
    // containment probability is given in percent
    iContainmentProbability /= 100.;
    
    // temporary variables
    vector< double > vEnergy;
    vector< double > vRes;
    vector< double > vResE;
    vector< double > vResN;   // number of entries
    TH1D* iTemp = 0;
    char iname[800];
    
    TH1D i_h_nbootstrap( "i_h_nbootstrap", "", 1000., 0., 1. );
    
    double i_energy = 0.;
    
    //////////////////////////////////////////////////////////////////////////////
    // loop over all energy bins and project each bin into a TH1D
    for( int i = 1; i <= iHistogram->GetNbinsX(); i++ )
    {
        i_h_nbootstrap.Reset();
        
        // define temporary histogram and fill with projection
        sprintf( iname, "iH_%d", i );
        iTemp = iHistogram->ProjectionY( iname, i, i );
        if( !iTemp )
        {
            continue;
        }
        i_energy = iHistogram->GetXaxis()->GetBinCenter( i );
        
        double i_entries = iTemp->GetEntries();
        i_entries -= iTemp->GetBinContent( 0 );
        i_entries -= iTemp->GetBinContent( iTemp->GetNbinsX() + 1 );
        
        if( i_entries < iMinRequiredEvents )
        {
            continue;
        }
        
        TH1D* hBootS = ( TH1D* )iTemp->Clone( "BootstrapClone" );
        
        // bootstrap to calculate errors on energy resolution
        // (don't bootstart if nBootStrap = 1 )
        for( unsigned int b = 0; b < fMProp_nBootStrap; b++ )
        {
            // bootstrap only if more than one cycle is
            // required
            if( fMProp_nBootStrap > 1 )
            {
                hBootS->Reset();
                hBootS->FillRandom( iTemp, i_entries );
            }
            
            //////////////////////////////////////////////////////////
            // calculate containment
            double iTotSum = hBootS->GetEntries();
            if( iTotSum  > 0. )
            {
                double iTempSum = 0.;
                for( int j = 1; j <= hBootS->GetNbinsX(); j++ )
                {
                    iTempSum += hBootS->GetBinContent( j );
                    if( iTempSum / iTotSum  > iContainmentProbability )
                    {
                        i_h_nbootstrap.Fill( hBootS->GetBinCenter( j ) );
                        break;
                    }
                }
            }
        }
        // require a minimum number of events to calculate
        // a good containment radius and RMS
        if( iTemp->GetEntries() > iMinRequiredEvents )
        {
            vEnergy.push_back( i_energy );
            vRes.push_back( i_h_nbootstrap.GetMean() );
            vResE.push_back( i_h_nbootstrap.GetRMS() );
            vResN.push_back( iTemp->GetEntries() );
        }
        if( iTemp )
        {
            delete iTemp;
        }
    }
    ////////////////////////////////////
    // check results for consistency
    //
    // require that each bin contains a
    // number of events which exceeds 5%
    // of the neighbouring bins
    
    for( unsigned i = 0; i < vEnergy.size(); i++ )
    {
        // left bin
        if( i > 0 && vResN[i - 1] > 0. )
        {
            if( vResN[i] / vResN[i - 1] < iMinBinFraction )
            {
                vResN[i] = 0.;
            }
        }
        // right bin
        if( vEnergy.size() > 1 && i < vEnergy.size() - 2 && vResN[i + 1] > 0. )
        {
            if( vResN[i] / vResN[i + 1] < iMinBinFraction )
            {
                vResN[i] = 0.;
            }
        }
    }
    
    // fill graph with angular resolution
    int z = 0;
    for( unsigned i = 0; i < vEnergy.size(); i++ )
    {
        if( vResN[i] > 0. )
        {
            g->SetPoint( z, vEnergy[i], vRes[i] );
            g->SetPointError( z, 0., vResE[i] );
            z++;
        }
    }
    
    return g;
}


/*
 * fill resolution into a histogram (h) for a given containment radius
 *
 * different ways to obtain resolution:
 * - read from file
 * - recalculate from 2D histograms
 *
 *  used for angular energy resolution
 *  (depends on iResolutionTreeName (t_energy_resolution or t_angular_resolution))
 *
 */
bool VInstrumentResponseFunctionReader::fillResolutionHistogram( TH1F* h, string iContainmentRadius,
        string iResolutionTreeName, bool iEnergyAxis_MC )
{
    if( !h )
    {
        return false;
    }
    
    // inconsistent naming: trees 68% trees have no suffix,
    // all other trees have the containment probability in tree
    // name
    if( iContainmentRadius != "68" )
    {
        iResolutionTreeName += "_0" + iContainmentRadius + "p";
    }
    
    for( unsigned int j = 0; j < fIRF_TreeNames.size(); j++ )
    {
        if( fIRF_TreeNames[j] == iResolutionTreeName )
        {
            if( j < fIRF_Data.size() && fIRF_Data[j] )
            {
                TGraphErrors* g = 0;
                if( iResolutionTreeName == "t_energy_resolution" )
                {
                    g = gEnergyResolution;
                }
                else if( iResolutionTreeName.find( "t_angular_resolution" ) != string::npos )
                {
                    cout << "VInstrumentResponseFunctionReader::fillResolutionHistogram ";
                    cout << "filling angular resolution (" << iContainmentRadius << "\% containment";
                    if( iEnergyAxis_MC )
                    {
                        cout << ", energy axis is true energy)";
                    }
                    else
                    {
                        cout << ", energy axis is reconstructed energy)";
                    }
                    cout << endl;
                    if( iEnergyAxis_MC )
                    {
                        g = getAngularResolution( fIRF_Data[j]->f2DHisto[VInstrumentResponseFunctionData::E_DIFF_MC],
                                                  atof( iContainmentRadius.c_str() ) );
                    }
                    else
                    {
                        g = getAngularResolution( fIRF_Data[j]->f2DHisto[VInstrumentResponseFunctionData::E_DIFF],
                                                  atof( iContainmentRadius.c_str() ) );
                    }
                    if( iContainmentRadius == "68" )
                    {
                        gAngularResolution = g;
                        // clone the 2D angular difference histogram to be used later
                        // (only done once for 68% containment radius, at it is the same
                        // for all raddii)
                        h2DAngularPSF = fIRF_Data[j]->f2DHisto[VInstrumentResponseFunctionData::E_DIFF];
                        if( VInstrumentResponseFunctionData::E_DIFF_MC < fIRF_Data[j]->f2DHisto.size() )
                        {
                            h2DAngularPSFEMC = fIRF_Data[j]->f2DHisto[VInstrumentResponseFunctionData::E_DIFF_MC];
                        }
                        else
                        {
                            h2DAngularPSFEMC = 0;
                        }
                    }
                    else if( iContainmentRadius == "80" )
                    {
                        gAngularResolution80 = g;
                    }
                    else if( iContainmentRadius == "95" )
                    {
                        gAngularResolution95 = g;
                    }
                }
                else
                {
                    g = fIRF_Data[j]->fResolutionGraph[VInstrumentResponseFunctionData::E_DIFF];
                    cout << "VInstrumentResponseFunctionReader::fillResolutionHistogram reading ";
                    cout << " resolution graph (" << iResolutionTreeName << ") from file " << endl;
                }
                
                if( g )
                {
                    double x = 0.;
                    double y = 0.;
                    double yE = 0.;
                    for( int i = 0; i < g->GetN(); i++ )
                    {
                        g->GetPoint( i, x, y );
                        yE = g->GetErrorY( i );
                        if( y > 0. )
                        {
                            h->SetBinContent( h->FindBin( x ), y );
                            h->SetBinError( h->FindBin( x ), yE );
                        }
                    }
                }
            }
        }
    }
    return true;
}
