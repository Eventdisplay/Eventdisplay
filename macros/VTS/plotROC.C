/*  plot background rejection vs signal efficiency
 *  from the TMVA results files
 *
 *  requires .L $EVNDISPSYS/lib/libVAnaSum.so
 *  e.g.
 *  .L $EVNDISPSYS/lib/libVAnaSum.so
 *  .L plotROC.C
 *  plotROC( "$VERITAS_USER_DATA_DIR/TMVA/", 0, 4, 0, "mva", "BDT", 10, 0.6 );
 */

#include <sstream>
#include <string>
#include <vector>

TGraph* getROCGraph( string iBDTFile, unsigned int iMethodN, string iMethodName, double minSignal = 0. )
{
    // open MVA file
    TFile* iFile = new TFile( iBDTFile.c_str() );
    if( iFile->IsZombie() )
    {
        cout << "Error: file not found: " << iFile << endl;
        return 0;
    }
    char hname[200];
    sprintf( hname, "Method_BDT_%d/BDT_%d/", iMethodN, iMethodN );
    if( !iFile->cd( hname ) )
    {
        cout << "Error: directories with BDT histograms not found:" << hname << endl;
        return 0;
    }
    sprintf( hname, "MVA_BDT_%d_effS", iMethodN );
    TH1D*  effS = ( TH1D* )gDirectory->Get( hname );
    sprintf( hname, "MVA_BDT_%d_effB", iMethodN );
    TH1D*  effB = ( TH1D* )gDirectory->Get( hname );
    if( !effS || !effB )
    {
        cout << "Error: histograms not found" << endl;
        cout << effS << "\t" << effB << endl;
        return 0;
    }
    TGraph* iROC = new TGraph( 1 );
    iROC->SetMarkerStyle( 7 );
    iROC->SetTitle( "" );
    int z = 0;
    for( int i = 1; i <= effS->GetNbinsX(); i++ )
    {
        if( effS->GetBinContent( i ) > minSignal && effS->GetBinContent( i ) > 0. &&  1. - effB->GetBinContent( i ) > 0. )
        {
            iROC->SetPoint( z, effS->GetBinContent( i ), 1. - effB->GetBinContent( i ) );
            z++;
        }
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
void plotROC( string iDirectory, int iEnergy_min = 0, int iEnergy_max = 8,
              int iZe = -1,
              string iMethodFileName = "mva",
              string iMethodName = "BDT", unsigned int iNMethods = 1,
              double minSignal = 0. )
{
    vector< TGraph* > fROC_Graph;
    vector< double > fROC_E_min;
    vector< double > fROC_E_max;
    
    vector< string > fROC_FileName;
    
    for( unsigned int i = iEnergy_min; i <= iEnergy_max; i++ )
    {
        for( int j = iZe; j < iZe + 1; j++ )
        {
            ostringstream iFile;
            iFile << iDirectory << "/" << iMethodFileName;
            iFile << "_" << i << "_" << j;
            iFile << ".root";
            fROC_FileName.push_back( iFile.str() );
        }
    }
    
    // get energies
    for( unsigned int i = 0; i < fROC_FileName.size(); i++ )
    {
        fROC_E_min.push_back( getROCEnergy( fROC_FileName[i].c_str(), true ) );
        fROC_E_max.push_back( getROCEnergy( fROC_FileName[i].c_str(), false ) );
    }
    
    TCanvas* cROC = new TCanvas( "cROC", "ROC", 10, 10, 900, 600 );
    cROC->SetGridx( 0 );
    cROC->SetGridy( 0 );
    cROC->Divide( 2, TMath::Ceil( ( iEnergy_max - iEnergy_min + 1 ) / 2 + 0.5 ) );
    cROC->Draw();
    
    
    // loop over all mvas and retrieve ROC
    for( unsigned int i = 0; i < fROC_FileName.size(); i++ )
    {
        cROC->cd( i + 1 );
        
        for( unsigned m = 0; m < iNMethods; m++ )
        {
            cout << "Reading from " << fROC_FileName[i] << ", method " << m << endl;
            TGraph* iG = getROCGraph( fROC_FileName[i].c_str(), m, iMethodName, minSignal );
            if( m + 1 != 10 )
            {
                iG->SetLineColor( m + 1 );
                iG->SetMarkerColor( m + 1 );
            }
            else
            {
                iG->SetLineColor( 11 );
                iG->SetMarkerColor( 11 );
            }
            cout << i << "\t" << iMethodName << "\t" << m << endl;
            
            if( m == 0 )
            {
                iG->Draw( "apl" );
                iG->GetHistogram()->GetXaxis()->SetTitle( "signal efficiency" );
                iG->GetHistogram()->GetYaxis()->SetTitle( "background rejection" );
            }
            else
            {
                iG->Draw( "pl" );
            }
        }
    }
    
    // plot FAUC
    TCanvas* cFAUC = new TCanvas( "cFAUC", "FAUC", 300, 300, 600, 600 );
    cFAUC->SetLogx( 1 );
    cFAUC->Draw();
    
    for( unsigned m = 0; m < iNMethods; m++ )
    {
        TGraph* iF = new TGraph( 1 );
        iF->SetTitle( "" );
        
        for( unsigned int i = 0; i < fROC_FileName.size(); i++ )
        {
            TGraph* iG = getROCGraph( fROC_FileName[i].c_str(), m, iMethodName, minSignal );
            
            iF->SetPoint( i, TMath::Power( 10., 0.5 * ( fROC_E_min[i] + fROC_E_max[i] ) ), iG->Integral() );
        }
        cout << "_______________________________ " << m << endl;
        iF->Print();
        if( m + 1 != 10 )
        {
            iF->SetLineColor( m + 1 );
            iF->SetMarkerColor( m + 1 );
        }
        else
        {
            iF->SetLineColor( 11 );
            iF->SetMarkerColor( 11 );
        }
        iF->SetMarkerStyle( 20 );
        iF->SetMarkerSize( 2. );
        if( m == 0 )
        {
            iF->Draw( "apl" );
            iF->GetHistogram()->GetXaxis()->SetTitle( "Energy [TeV]" );
            iF->GetHistogram()->GetYaxis()->SetTitle( "Area under ROC" );
        }
        else
        {
            iF->Draw( "pl" );
        }
    }
    
}

