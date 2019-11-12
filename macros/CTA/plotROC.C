/*  plot background rejection vs signal efficiency
 *  from the TMVA results files
 *
 *  requires .L $EVNDISPSYS/lib/libVAnaSum.so
 *  e.g.
 *  .L $EVNDISPSYS/lib/libVAnaSum.so
 *  .L plotROC.C
 *   plotROC( "$CTA_USER_DATA_DIR/TMVA/", 0, 8 );
 */

#include <sstream>
#include <string>
#include <vector>

TGraph* getROCGraph( string iBDTFile )
{
    // open MVA file
    TFile* iFile = new TFile( iBDTFile.c_str() );
    if( iFile->IsZombie() )
    {
        cout << "Error: file not found: " << iFile << endl;
        return 0;
    }
    if( !iFile->cd( "Method_BDT/BDT_0/" ) )
    {
        cout << "Error: directories with BDT histograms not found" << endl;
        return 0;
    }
    TH1D*  effS = ( TH1D* )gDirectory->Get( "MVA_BDT_0_effS" );
    TH1D*  effB = ( TH1D* )gDirectory->Get( "MVA_BDT_0_effB" );
    if( !effS || !effB )
    {
        cout << "Error: histograms not found" << endl;
        cout << effS << "\t" << effB << endl;
        return 0;
    }
    TGraph* iROC = new TGraph( 1 );
    iROC->SetMarkerStyle( 7 );
    iROC->SetTitle( "" );
    for( int i = 1; i <= effS->GetNbinsX(); i++ )
    {
        iROC->SetPoint( i, effS->GetBinContent( i ), 1. - effB->GetBinContent( i ) );
    }
    
    return iROC;
}

/*
 * get minimum / maximum energy
 *
 */
double getROCEnergy( string iBDTFile, bool iMin )
{
    // open MVA file
    TFile* iFile = new TFile( iBDTFile.c_str() );
    if( iFile->IsZombie() )
    {
        cout << "Error: file not found: " << iFile << endl;
        return 0.;
    }
    // get energy range from title
    ostringstream iTitle;
    VTMVARunDataEnergyCut* iE = ( VTMVARunDataEnergyCut* )iFile->Get( "fDataEnergyCut" );
    if( iE )
    {
        if( iMin )
        {
            return iE->fEnergyCut_Log10TeV_min;
        }
        else
        {
            return iE->fEnergyCut_Log10TeV_max;
        }
    }
    
    return 0.;
}

/*
 * plot a series of BDT ROC curves
 *
 * (function called by user)
 */
void plotROC( string iDirectory, int iFile_min = 0, int iFile_max = 8 )
{
    vector< TGraph* > fROC_Graph;
    vector< double > fROC_E_min;
    vector< double > fROC_E_max;
    
    for( unsigned int i = iFile_min; i <= iFile_max; i++ )
    {
        ostringstream iFile;
        iFile << iDirectory << "/BDT_" << i << ".root";
        fROC_Graph.push_back( getROCGraph( iFile.str() ) );
        fROC_E_min.push_back( getROCEnergy( iFile.str(), true ) );
        fROC_E_max.push_back( getROCEnergy( iFile.str(), false ) );
    }
    
    // graph with AUC vs energy
    TGraph* fAUC = new TGraph( 1 );
    fAUC->SetMarkerStyle( 21 );
    fAUC->SetTitle( "" );
    
    TCanvas* cROC = new TCanvas( "cROC", "ROC", 10, 10, 600, 400 );
    cROC->SetGridx( 0 );
    cROC->SetGridy( 0 );
    cROC->Draw();
    
    TLegend* iL = new TLegend( 0.15, 0.15, 0.55, 0.55 );
    
    unsigned int z = 0;
    for( unsigned int i = 0; i < fROC_Graph.size(); i++ )
    {
        if( !fROC_Graph[i] )
        {
            continue;
        }
        fROC_Graph[i]->SetLineColor( i + 1 );
        fROC_Graph[i]->SetMarkerColor( i + 1 );
        ostringstream iLL;
        iLL << fROC_E_min[i] << " < log10( E/TeV ) <" << fROC_E_max[i];
        iL->AddEntry( fROC_Graph[i], iLL.str().c_str(), "l" );
        if( z == 0 )
        {
            fROC_Graph[i]->Draw( "apl" );
            fROC_Graph[i]->GetHistogram()->GetXaxis()->SetTitle( "signal efficiency" );
            fROC_Graph[i]->GetHistogram()->GetYaxis()->SetTitle( "background rejection" );
        }
        else
        {
            fROC_Graph[i]->Draw( "pl" );
        }
        fAUC->SetPoint( z, TMath::Power( 10., 0.5 * ( fROC_E_min[i] + fROC_E_max[i] ) ), fROC_Graph[i]->Integral() );
        z++;
    }
    iL->Draw();
    
    // AUC
    
    TCanvas* cAUC = new TCanvas( "cAUC", "AUC", 610, 10, 600, 400 );
    cAUC->SetGridx( 0 );
    cAUC->SetGridy( 0 );
    cAUC->SetLogx( 1 );
    cAUC->Draw();
    
    fAUC->Draw( "ap" );
    fAUC->GetHistogram()->GetXaxis()->SetTitle( "Energy [TeV]" );
    fAUC->GetHistogram()->GetYaxis()->SetTitle( "Area under ROC" );
    fAUC->Print();
}

