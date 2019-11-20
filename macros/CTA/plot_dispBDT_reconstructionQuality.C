/*
 * plot reconstruction quality from disp training
 * (telescopes types are set for CTA,
 *  change for other telescopes)
 *
 *  check for overtraining:
 *  plot disp/dispCore/energy resolution of testing
 *  and training tree
 *
 *  reads test/train trees produced
 *  by TMVA training
 *  trainTMVAforAngularReconstruction
 *
 *  execute with .L plot_dispBDT_reconstructionQuality.C++
 *
 *  expect a series of BDT training results in
 *  directories with consecutive numbering
 *
*/

#include "TCanvas.h"
#include "TCollection.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObject.h"
#include "TROOT.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

/*
 * plot distribution of BDT training variables
 *
 * show distributions for events with
 * large deviations
 *
 */
void plotVariables( TTree* t,
                    string iDispType,
                    string iTelType )
{
    if( !t )
    {
        return;
    }
    
    double fEMin = 10.;
    double fEMax = 300.;
    
    /////////////////////////
    // variables and their
    // hardwired min/max
    
    vector< string > fVarType;
    vector< float > fVarMin;
    vector< float > fVarMax;
    // size needs to be the first variable
    fVarType.push_back( "size" );
    fVarMin.push_back( 2. );
    fVarMax.push_back( 10. );
    // disp needs to be the second variable for BDTdisp
    if( iDispType == "BDTDisp" )
    {
        fVarType.push_back( "disp" );
        fVarMin.push_back( 0. );
        fVarMax.push_back( 20. );
    }
    fVarType.push_back( "width" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 2. );
    fVarType.push_back( "length" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 10. );
    fVarType.push_back( "wol" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 1. );
    fVarType.push_back( "tgrad_x_T_tgrad_x" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 5. );
    fVarType.push_back( "cross" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 25. );
    fVarType.push_back( "asym" );
    fVarMin.push_back( -2.5 );
    fVarMax.push_back( 2.5 );
    fVarType.push_back( "loss" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 0.3 );
    if( iDispType == "BDTDispEnergy" )
    {
        fVarType.push_back( "EHeight" );
        fVarMin.push_back( 0. );
        fVarMax.push_back( 100. );
        fVarType.push_back( "Rcore" );
        fVarMin.push_back( 0. );
        fVarMax.push_back( 1000. );
    }
    else
    {
        fVarType.push_back( "dist" );
        fVarMin.push_back( 0. );
        fVarMax.push_back( 4. );
        fVarType.push_back( "fui" );
        fVarMin.push_back( 0. );
        fVarMax.push_back( 1. );
    }
    fVarType.push_back( "Rcore-MCrcore" );
    fVarMin.push_back( -500. );
    fVarMax.push_back( 500. );
    
    vector< float > fVar( fVarType.size(), 0. );
    
    vector< TH1D* > fHAll;
    vector< TH1D* > fHDivU;
    vector< TH1D* > fHDivD;
    
    float MCe0 = 0.;
    float MCrcore = 0.;
    float disp_true = 0.;
    t->SetBranchAddress( "MCe0", &MCe0 );
    t->SetBranchAddress( "MCrcore", &MCrcore );
    
    ostringstream BDTvar;
    BDTvar << "BDT_" << iTelType;
    t->SetBranchAddress( BDTvar.str().c_str(), &disp_true );
    
    // trees and histograms
    for( unsigned int v = 0; v < fVarType.size(); v++ )
    {
        ostringstream iHName;
        ostringstream iHNameU;
        ostringstream iHNameD;
        if( fVarType[v].find( "MC" ) == string::npos )
        {
            t->SetBranchAddress( fVarType[v].c_str(), &fVar[v] );
            
            iHName << "h" << iDispType << iTelType << "_" << fVarType[v];
            iHNameU << "hDivU" << iDispType << iTelType << "_" << fVarType[v];
            iHNameD << "hDivD" << iDispType << iTelType << "_" << fVarType[v];
        }
        else
        {
            iHName << "h" << iDispType << iTelType << "_" << v;
            iHNameU << "hDivU" << iDispType << iTelType << "_" << v;
            iHNameD << "hDivD" << iDispType << iTelType << "_" << v;
        }
        
        fHAll.push_back( new TH1D( iHName.str().c_str(), "", 100, fVarMin[v], fVarMax[v] ) );
        fHAll.back()->SetLineColor( 1 );
        fHAll.back()->SetLineWidth( 2 );
        fHAll.back()->SetStats( 0 );
        fHAll.back()->SetXTitle( fVarType[v].c_str() );
        
        fHDivU.push_back( new TH1D( iHNameU.str().c_str(), "", 100, fVarMin[v], fVarMax[v] ) );
        fHDivU.back()->SetLineColor( 801 );
        fHDivU.back()->SetLineWidth( 2 );
        fHDivU.back()->SetStats( 0 );
        fHDivU.back()->SetXTitle( fVarType[v].c_str() );
        fHDivD.push_back( new TH1D( iHNameD.str().c_str(), "", 100, fVarMin[v], fVarMax[v] ) );
        fHDivD.back()->SetLineColor( 633 );
        fHDivD.back()->SetLineWidth( 2 );
        fHDivD.back()->SetStats( 0 );
        fHDivD.back()->SetXTitle( fVarType[v].c_str() );
    }
    
    ///////////////////////////////////////////////////////
    // loop over all entries and fill the histograms
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        if( MCe0 < fEMin || MCe0 > fEMax )
        {
            continue;
        }
        
        for( unsigned int v = 0; v < fHAll.size(); v++ )
        {
            if( v == fVar.size() - 1 )
            {
                fVar[v] = fVar[v - 1] - MCrcore;
            }
            
            fHAll[v]->Fill( fVar[v] );
            
            // plot all misreconstructed events
            if( iDispType == "BDTDispEnergy" )
            {
                if( TMath::Power( 10., disp_true * fVar[0] ) / MCe0 > 1.5 )
                {
                    fHDivU[v]->Fill( fVar[v] );
                }
                if( TMath::Power( 10., disp_true * fVar[0] ) / MCe0 < 0.5 )
                {
                    fHDivD[v]->Fill( fVar[v] );
                }
            }
            else if( iDispType == "BDTDisp" )
            {
                if( disp_true > 0. && TMath::Abs( fVar[1]/disp_true ) > 1.25 )
                {
                    fHDivU[v]->Fill( fVar[v] );
                }
                if( disp_true > 0. && TMath::Abs( fVar[1]/disp_true ) < 0.75 )
                {
                    fHDivD[v]->Fill( fVar[v] );
                }
            }
        }
    }
    ostringstream iCanvasName;
    iCanvasName << "c" << iTelType << "_" << iDispType;
    TCanvas* c = new TCanvas( iCanvasName.str().c_str(), "", 100, 100, 1600, 1600 );
    c->Divide( 4, 3 );
    
    ///////////////////////////////////////////////////////
    // plot everything
    for( unsigned int v = 0; v < fVarType.size(); v++ )
    {
        TPad* p = ( TPad* )c->cd( v + 1 );
        p->SetGridx( 0 );
        p->SetGridy( 0 );
        p->SetLogy( 1 );
        p->Draw();
        
        fHAll[v]->SetMinimum( 1 );
        fHAll[v]->Draw();
        fHDivU[v]->Draw( "same" );
        fHDivD[v]->Draw( "same" );
    }
    ostringstream iPrintName;
    iPrintName << "VAR" << iTelType << "-" << iDispType << ".pdf";
    c->Print( iPrintName.str().c_str() );
}

/*
 * calculate and plot disp resolution
 * for the given eventdisp
 */
double plotAngularResolution( string iDispType, TTree* t, int iColor, string iTelType, bool bFirst,
                            unsigned int iCounter, string iName, unsigned int iMaxNDir,
                            bool iPrint = false )
{
    if( !t )
    {
        return 0.; 
    }
    // 2D histogram with error
    ostringstream iMigMatrix;
    iMigMatrix << "MigMatrixDisp_" << iDispType << iTelType << "_" << iCounter;
    TH2F* migmatrix = 0;
    if( iDispType == "BDTDisp" )
    {
        migmatrix = new TH2F( iMigMatrix.str().c_str(), "", 22, -2., 2.4, 500, 0., 5. );
    }
    else
    {
        migmatrix = new TH2F( iMigMatrix.str().c_str(), "", 22, -2., 2.4, 500, 0., 500. );
    }
    
    float MCe0 = 0.;
    float MCrcore = 0.;
    float disp = 0.;
    float dispBDT = 0.;
    
    t->SetBranchAddress( "MCe0", &MCe0 );
    if( iDispType == "BDTDisp" )
    {
        t->SetBranchAddress( "disp", &disp );
    }
    else
    {
        t->SetBranchAddress( "dispCore", &disp );
    }
    ostringstream iBDTV;
    iBDTV << "BDT_" << iTelType;
    if( t->GetBranchStatus( iBDTV.str().c_str() ) )
    {
        t->SetBranchAddress( iBDTV.str().c_str(), &dispBDT );
    }
    else
    {
        return 0.;
    }
    cout << "Filling migration matrix for  " << iBDTV.str() << endl;
    
    // fill migmatrix
    for( int e = 0; e < t->GetEntries(); e++ )
    {
        t->GetEntry( e );
        
        migmatrix->Fill( log10( MCe0 ), TMath::Abs( disp - dispBDT ) );
    }
    TCanvas* cDispres = 0;
    ostringstream iCanvasName;
    iCanvasName << "dispResolution" << iTelType;
    if( bFirst )
    {
        cDispres = new TCanvas( iCanvasName.str().c_str(), "#Delta disp (68\% containment)", 20, 20, 700, 700 );
        cDispres->SetLeftMargin( 0.13 );
    }
    else
    {
        cDispres = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( iCanvasName.str().c_str() );
    }
    if( !cDispres )
    {
        return 0.;
    }
    cDispres->cd();
    char hname[600];
    
    TGraphErrors* Dres = new TGraphErrors( 1 );
    Dres->SetMarkerStyle( 20 );
    int z = 0;
    
    //////////////////////////////////////////////////////////////////////////////////
    // calculate 68% value
    
    double xq[] = { 0.5,  0.68, 0.95 };
    double yq[] = { 0.0,  0.0, 0.0  };
    
    cout << "----------- " << iCounter << "  " << iName << " -------------" << endl;
    for( int i = 1; i <= migmatrix->GetNbinsX(); i++ )
    {
        ostringstream iH;
        iH << "hDisp" << i;
        TH1F* h = ( TH1F* )migmatrix->ProjectionY( iH.str().c_str(), i, i );
        if( h->GetEntries() > 10 )
        {
            h->GetQuantiles( 3, yq, xq );
            
            Dres->SetPoint( z, migmatrix->GetXaxis()->GetBinCenter( i ), yq[1] );
            z++;
            cout << iDispType << " resolution at  " << TMath::Power( 10., migmatrix->GetXaxis()->GetBinCenter( i ) )  << " TeV: ";
            cout << yq[1] << endl;
        }
    }
    cDispres->cd();
    if( iMaxNDir < 10 )
    {
        Dres->SetLineColor( iColor );
        Dres->SetLineWidth( 3. );
    }
    
    if( bFirst )
    {
        ostringstream iHEres;
        iHEres << "HDrec_" << iTelType;
        TH1D* h = new TH1D( iHEres.str().c_str(), "", 100, -2., 2.5 );
        if( iDispType == "BDTDisp" )
        {
            h->SetMaximum( 0.5 );
        }
        else
        {
            h->SetMaximum( 100. );
        }
        h->SetStats( 0 );
        h->SetXTitle( "log_{10} Energy/TeV" );
        ostringstream iTEres;
        iTEres << "#Delta " << iDispType << " (68\% containment)";
        h->SetYTitle( iTEres.str().c_str() );
        h->GetYaxis()->SetTitleOffset( 1.4 );
        h->Draw();
    }
    Dres->Draw( "l" );


    if( iPrint )
    {
        ostringstream iPrintName;
        iPrintName << "BDTdisp-Disp-" << iTelType << "-" << iName << ".pdf";
        cDispres->Print( iPrintName.str().c_str() );
    }
    
    // calculate figure of merit for energy resolution
    double x = 0.;
    double y = 0.;
    double iMetric = 0.;
    double iMetricN = 0.;
    for( int i = 0; i < Dres->GetN(); i++ )
    {
         Dres->GetPoint( i, x, y );
         if( x > -1.5 && x < 1.5 )
         {
              iMetric += y;
              iMetricN ++;
         }
    }
    if( iMetricN > 0. )
    {
         return iMetric / iMetricN;
    } 
    return 0.; 
}

/*
 * calculate and plot energy resolution
 * for the given events
 *
*/
double plotEnergyResolution( TTree* t, int iColor, string iTelType, bool bFirst,
                           int bPlotMigrationMatrix, unsigned int iCounter,
                           string iName, unsigned int iMaxNDir, 
                           bool iPrint = false )
{
    if( !t )
    {
        return 0.;
    }
    // hardwired parameters
    double fEnergy_TeV_Eres_min = 0.02;
    bool bUseMedian = true;
    bool bXaxisIsEtrue = true;
    
    // migration matrix
    ostringstream iMigMatrix;
    iMigMatrix << "MigMatrix_" << iTelType << "_" << iCounter;
    TH2F* migmatrix = new TH2F( iMigMatrix.str().c_str(), "", 500, -2., 2.5, 500, -2., 2.5 );
    
    float size = 0.;
    float MCe0 = 0.;
    float dispBDT = 0.;
    
    t->SetBranchAddress( "size", &size );
    t->SetBranchAddress( "MCe0", &MCe0 );
    ostringstream iBDTV;
    iBDTV << "BDT_" << iTelType;
    if( t->GetBranchStatus( iBDTV.str().c_str() ) )
    {
        t->SetBranchAddress( iBDTV.str().c_str(), &dispBDT );
    }
    else
    {
        return 0.;
    }
    
    // fill migmatrix
    for( int e = 0; e < t->GetEntries(); e++ )
    {
        t->GetEntry( e );
        
        migmatrix->Fill( dispBDT * size, log10( MCe0 ) );
    }
    
    // plot migration matrix
    if( bPlotMigrationMatrix )
    {
        cout << "Plotting migration matrix for " << iCounter << endl;
        TCanvas* c_mig = new TCanvas( "migrationMatrix", "energy_{MC} (TeV)", 10, 10, 600, 600 );
        migmatrix->SetStats( 0 );
        migmatrix->Draw( "colz" );
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    // calculate energy resolution as the 68% interval around the most probably value
    // for the energy reconstruction
    
    TCanvas* cEres = 0;
    ostringstream iCanvasName;
    iCanvasName << "energyResolution" << iTelType;
    if( bFirst )
    {
        cEres = new TCanvas( iCanvasName.str().c_str(), "#Delta E/E (68\% containment)", 20, 20, 700, 700 );
        cEres->SetLeftMargin( 0.13 );
    }
    else
    {
        cEres = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( iCanvasName.str().c_str() );
    }
    if( !cEres )
    {
        return 0.;
    }
    cEres->cd();
    char hname[600];
    
    TGraphErrors* Eres = new TGraphErrors( 1 );
    Eres->SetMarkerStyle( 20 );
    int z = 0;
    int i_z = 0;
    int zC = 1;
    float eres = 0.;
    float eres_temp = 0.;
    float eres_tempRMS = 0.;
    double energy = 0.;
    double energy_mop_med = 0.;
    
    int nbins = 0;
    if( bXaxisIsEtrue )
    {
        nbins = migmatrix->GetNbinsY();
    }
    else
    {
        nbins = migmatrix->GetNbinsX();
    }
    
    // energy resolution is averaged over these number of bins in the
    // migration matrix
    int i_zoffset_bins = 20;
    TH1D i_h_zoffsetsHisto( "i_h_zoffsetsHisto", "", 1000., 0., 1. );
    cout << "----------- " << iCounter << "  " << iName << " -------------" << endl;
    
    ///////////////////////////////////////////////////////////////////////
    // loop over axis in energy of migration matrix
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
        if( fEnergy_TeV_Eres_min > 0. && energy < log10( fEnergy_TeV_Eres_min ) )
        {
            continue;
        }
        
        // require at least a few entries
        if( hDist->GetEntries() < 10 )
        {
            continue;
        }
        
        ///////////////////////////////////
        // Fermi LAT method
        // half-width of +-34% intervall around most probably energy
        
        // most probable energy
        double energy_mprob = hDist->GetXaxis()->GetBinCenter( hDist->GetMaximumBin() );
        
        // cummulative distribution
        sprintf( hname, "%s_CUMU", hDist->GetName() );
        TH1F* hcum = ( TH1F* )hDist->Clone( hname );
        hcum->Reset();
        
        hcum->SetBinContent( 1, hDist->GetBinContent( 1 ) );
        for( int j = 2; j <= hDist->GetNbinsX(); j++ )
        {
            hcum->SetBinContent( j, hDist->GetBinContent( j ) + hcum->GetBinContent( j - 1 ) );
        }
        if( hcum->GetMaximum() > 0. )
        {
            hcum->Scale( 1. / hcum->GetMaximum() );
        }
        // get energy of median
        double energy_median = 0.5 * ( hcum->GetXaxis()->GetBinCenter( hcum->FindFirstBinAbove( 0.5 ) )
                                       + hcum->GetXaxis()->GetBinCenter( hcum->FindFirstBinAbove( 0.5 ) - 1 ) );
        // use median or mprob?
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
        float prob_mprob_minus = prob_mprob - 0.34;
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
        float prob_mprob_plus  = prob_mprob + 0.34;
        if( prob_mprob_plus > 1. )
        {
            prob_mprob_plus = 1. - 1.e-3;
        }
        int prob_plus  = hcum->FindFirstBinAbove( prob_mprob_plus );
        double energy_prob_plus = hcum->GetXaxis()->GetBinCenter( prob_plus );
        if( prob_plus > 2 )
        {
            energy_prob_plus = 0.5 * ( hcum->GetXaxis()->GetBinCenter( prob_plus ) + hcum->GetXaxis()->GetBinCenter( prob_plus - 1 ) );
        }
        // energy resolution
        eres = TMath::Power( 10., energy_prob_plus ) - TMath::Power( 10., energy_prob_minus );
        eres /= TMath::Power( 10., energy_mop_med );
        eres *= 0.5;  // half width
        
        eres_temp += eres;
        i_h_zoffsetsHisto.Fill( eres );
        i_z++;
        
        
        // average of i_zoffset_bins bins
        if( i_z == i_zoffset_bins )
        {
            double energy_graph_midpoint = migmatrix->GetXaxis()->GetBinCenter( i - 0.5 * i_zoffset_bins );
            eres_temp    = i_h_zoffsetsHisto.GetMean();
            eres_tempRMS = i_h_zoffsetsHisto.GetRMS();
            
            Eres->SetPoint( z, energy_graph_midpoint, eres_temp );
            Eres->SetPointError( z, 0., eres_tempRMS );
            
            cout << "Energy resolution at " << TMath::Power( 10., energy_graph_midpoint ) << " TeV: ( ";
            cout << eres_temp * 100. << " +- " << eres_tempRMS * 100. << " )%";
            cout << " (bin [" << i - i_zoffset_bins << ", " << i << "]";
            cout << endl;
            z++;
            i_h_zoffsetsHisto.Reset();
            i_z = 0;
        }
        
    }
    cEres->cd();
    if( iMaxNDir < 10 )
    {
        Eres->SetLineColor( iColor );
        Eres->SetLineWidth( 3. );
    }
    
    if( bFirst )
    {
        ostringstream iHEres;
        iHEres << "HErec_" << iTelType;
        TH1D* h = new TH1D( iHEres.str().c_str(), "", 100, -2., 2.5 );
        h->SetStats( 0 );
        h->SetXTitle( "log_{10} Energy/TeV" );
        h->SetYTitle( "#Delta E/E (68\% containment)" );
        h->GetYaxis()->SetTitleOffset( 1.4 );
        h->SetMaximum( 0.5 );
        h->Draw();
    }
    Eres->Draw( "l" );


    if( iPrint )
    {
        ostringstream iPrintName;
        iPrintName << "BDTdisp-DispEnergy-" << iTelType << "-" << iName << ".pdf";
        cEres->Print( iPrintName.str().c_str() );
    }
    
    // calculate figure of merit for energy resolution

    // calculate figure of merit for energy resolution
    double x = 0.;
    double y = 0.;
    double iMetric = 0.;
    double iMetricN = 0.;
    for( int i = 0; i < Eres->GetN(); i++ )
    {
         Eres->GetPoint( i, x, y );
         if( x > -1.5 && x < 1.5 )
         {
              iMetric += y;
              iMetricN ++;
         }
    }
    if( iMetricN > 0. )
    {
         return iMetric / iMetricN;
    }
    return 0.;
}

double plotResolution( string iDispType, TTree* t, int iColor, string iTelType, bool bFirst,
                     int bPlotMigrationMatrix, unsigned int iCounter,
                     string iName, unsigned int iMaxNDir = 0, bool iPrint = false )
{
    if( iDispType == "BDTDispEnergy" )
    {
        return plotEnergyResolution( t, iColor, iTelType, bFirst, bPlotMigrationMatrix, iCounter, iName, iMaxNDir, iPrint );
    }
    else if( iDispType == "BDTDisp" || iDispType == "BDTDispCore" )
    {
        return plotAngularResolution( iDispType, t, iColor, iTelType, bFirst, iCounter, iName, iMaxNDir, iPrint );
    }

    return 0.;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// only function called by user!

/*
 *   plot reconstruction quality for DISP reconstruction
 *
 *   at this point: only for dispEnergy
 *
 *
 *  fTelType.push_back( "138704810" ); //  LST
     fTelType.push_back( "10408618" ); //   FlashCam MST
     fTelType.push_back( "10408418" ); //   NectarCam MST
     fTelType.push_back( "201511619" ); //  ASTRI SST
     fTelType.push_back( "201309316" ); //  GC-SST
     fTelType.push_back( "909924" ); //     DC-SST
     fTelType.push_back( "207308707" ); //  SC-MST


     different plotting modes are:
     "Resolution"
     "TrainingVariables" (add a number in iPlotSingleHistograms to say which one should be plotted)
     "OverTraining"

*/
void plot_dispBDT_reconstructionQuality( string iDispType = "BDTDispEnergy",
        string iTelType = "138704810",
        string iDirectoryWithTrees = "./",
        string iPlottingMode = "Resolution",
        int iPlotSingleHistograms = 0,
        unsigned int nMaxDispDir = 363 )
{
    ////////////////////////
    // directories with BDT training data to be plotted
    vector< string > fDataDir;
    vector< string > fDataName;
    
    // loop over all training directories
    ostringstream iFileName;
    
    /*iFileName << "BDTdisp.Nb.3AL4-BN15.M69/" << iDispType << "/0deg/";
    fDataDir.push_back( iFileName.str() );
    iFileName.str( "" );
    fDataName.push_back( "M69" ); */

    if( iPlottingMode == "Resolution" )
    {
        for( unsigned int i = 1; i < nMaxDispDir; i++ )
         {
             iFileName.str( "" );
             iFileName << "BDTdisp.Nb.3AL4-BN15.M" << i << "/" << iDispType << "/0deg/";
             fDataDir.push_back( iFileName.str() );

             iFileName.str( "" );
             iFileName << "BDTdisp.Nb.3AL4-BN15.M" << i;
             fDataName.push_back( iFileName.str() );
         }
     }
     else
     {
         iFileName.str( "" );
         iFileName << "BDTdisp.Nb.3AL4-BN15.T" << iPlotSingleHistograms << "/" << iDispType << "/0deg/";
         fDataDir.push_back( iFileName.str() );

         iFileName.str( "" );
         iFileName << "BDTdisp.Nb.3AL4-BN15.T" << iPlotSingleHistograms;
         fDataName.push_back( iFileName.str() );
     }

     // figure of merrit to find best resolution
     map< double, unsigned int > fResolutionMetric;
    
    
    ///////////////////////////////////
    // get all test and training trees
    vector< TTree* > fTestTree;
    vector< TTree* > fTrainingTree;
    vector< string > fName;
    unsigned int z = 0;
    for( unsigned int t = 0; t < fDataDir.size(); t++ )
    {
        ostringstream iFileName;
        iFileName << iDirectoryWithTrees << fDataDir[t] << "/" << iDispType << "_" << iTelType << ".tmva.root";
        cout << "opening " << iFileName.str() << endl;
        TFile* f = new TFile( iFileName.str().c_str() );
        if( f->IsZombie() )
        {
            fTestTree.push_back( 0 );
            continue;
        }
        if( f->Get( "TestTree" ) )
        {
            fTestTree.push_back( ( TTree* )f->Get( "TestTree" ) );
            fTestTree.back()->SetLineColor( fTestTree.size() );
            fTestTree.back()->SetLineWidth( 2 );
            fName.push_back( fDataName[t] );
            cout << "Directory " << fDataDir[t] << ": ID " << fTestTree.size() - 1 << endl;
            z++;
            if( f->Get( "TrainTree" ) )
            {
                fTrainingTree.push_back( ( TTree* )f->Get( "TrainTree" ) );
                fTrainingTree.back()->SetLineColor( fTestTree.size() );
                fTrainingTree.back()->SetLineWidth( 2 );
            }
            else
            {
                fTrainingTree.push_back( 0 );
            }
        }
        else
        {
             fTestTree.push_back( 0 );
        }
    }
    cout << "Found " << z << " test trees" << endl;
    
    if( iPlottingMode == "Resolution" )
    {
        ////////////////////////////////////////////////////////////
        // plot resolution
        for( unsigned int i = 0; i < fTestTree.size(); i++ )
        {
            if( fTestTree[i] != 0 )
            {
                double i_metric = plotResolution( iDispType, fTestTree[i], i + 1, iTelType, ( i == 0 ), 
                             ( ( int )i == iPlotSingleHistograms ), i, fDataName[i], nMaxDispDir );
                fResolutionMetric[i_metric] = i;
            }
        }
        // plot metric (sorted)
        cout << "Resolution metrics: " << endl;
        unsigned int z = 0;
        for( map< double, unsigned int >::iterator it = fResolutionMetric.begin(); 
               it != fResolutionMetric.end(); ++it )
        {
             if( it->first > 0 )
             {
                 cout << "metrics " << it->second+1 << "\t\t" << it->first << endl;
                 if( z < 20 )
                 {
                        plotResolution( iDispType, fTestTree[it->second], z+29, iTelType, false,
                        false, it->second, fDataName[it->second], 1 );
                        z++;
                 }
             }
        } 
        z = 0;
        for( map< double, unsigned int >::iterator it = fResolutionMetric.begin(); 
               it != fResolutionMetric.end(); ++it )
        {
             if( it->first > 0 && z < 20 )
             {
                 cout << "T20 metrics " << it->second+1 << "\t\t" << it->first << endl;
             }
             z++;
        }
    }
    else if( iPlottingMode == "TrainingVariables" )
    {
        ////////////////////////////////////////////////////////////
        // plot individual variables and the dependency of the
        // enerygy resolution for them
        if( iPlotSingleHistograms >= 0
             && fTestTree.size() > 0
             && fTestTree[0] )
        {
            plotVariables( fTestTree[0], iDispType, iTelType );
        }
    }
    else if( iPlottingMode == "OverTraining" )
    {
        ////////////////////////////////////////////////////////////
        // plot resolution for Test and Training tree
        // this is a check for over training
        //
        if( iPlotSingleHistograms >= 0
                && fTestTree.size() > 0
                && fTestTree[0] 
                && fTrainingTree.size() > 0 ) 
        {
            plotResolution( iDispType, fTestTree[0], 1, iTelType, true, false, 0, fName[0], 1 );
            if( fTrainingTree[0] )
            {
                fTrainingTree[0]->SetLineColor( 2 );
                plotResolution( iDispType, fTrainingTree[0], 2, iTelType, false, false, 0, fName[0], 1, true );
            }
            else
            {
                cout << "no training tree found" << endl;
            }
            cout << "Black: test tree" << endl;
            cout << "Red: training tree" << endl;
        }
    }
    else
    {
        cout << "Warning: unknown plotting mode" << endl;
    }
    
}

void plotMVA( int iMVAID, string iSite = "North" )
{
      if( iSite == "North" )
      {
           plot_dispBDT_reconstructionQuality( "BDTDisp", "10408418", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDispEnergy", "10408418", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDisp", "138704810", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDispEnergy", "138704810", "./", "OverTraining", iMVAID, 1);
      }
      else
      {
           plot_dispBDT_reconstructionQuality( "BDTDisp", "10408618", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDispEnergy", "10408618", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDisp", "138704810", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDispEnergy", "138704810", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDisp", "201309316", "./", "OverTraining", iMVAID, 1);
           plot_dispBDT_reconstructionQuality( "BDTDispEnergy", "201309316", "./", "OverTraining", iMVAID, 1);
      }
}
    

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MVA
{
    public:
        unsigned int fID;
        string fName;
        float fEres;
        float fEresError;
        TGraphErrors* fEresGraph;
        
        MVA( unsigned int iColor = 1 );
        ~MVA() {}
};

MVA::MVA( unsigned int iColor )
{
    fEresGraph = new TGraphErrors( 1 );
    fEresGraph->SetTitle( "" );
    fEresGraph->SetMarkerStyle( 24 );
    fEresGraph->SetMarkerColor( iColor );
    fEresGraph->SetLineColor( iColor );
}

/*
 * use the printout of plot_dispBDT_reconstructionQuality.C
 * and print a selected number of resolution plots
 */
void plotResolutionFromASCIIFile( string iFile, float iEnergy, float EresMax = 0.3,
                                  unsigned int iReferenceID = 99999, double iMax = 0.5 )
{
    vector< MVA* > fData;
    
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        return;
    }
    
    string is_line;
    string itemp;
    unsigned int iID = 0;
    unsigned int z = 1;
    unsigned int t = 0;
    
    while( getline( is, is_line ) )
    {
        // new data set
        if( is_line.find( "---------" ) != string::npos )
        {
            stringstream is_stream( is_line );
            is_stream >> itemp;
            fData.push_back( new MVA( z ) );
            is_stream >> fData.back()->fID;
            is_stream >> fData.back()->fName;
            z++;
            t = 0;
        }
        else if( is_line.find( "resolution" ) != string::npos )
        {
            float e = atof( is_line.substr( is_line.find( "at" ) + 3, 6 ).c_str() );
            stringstream is_stream( is_line.substr( is_line.find( "TeV:" ) + 4 ) );
            float eres = 0.;
            float eresError = 0.;
            if( is_line.find( "Energy resolution" ) != string::npos )
            {
                is_stream >> itemp;
                is_stream >> eres;
                is_stream >> itemp;
                is_stream >> eresError;
                eres /= 100.;
                eresError /= 1000.;
            }
            else
            {
                is_stream >> eres;
                eresError = 0.;
            }
            if( eres > 1.e-3 )
            {
                fData.back()->fEresGraph->SetPoint( t, log10( e ), eres );
                fData.back()->fEresGraph->SetPointError( t, 0., eresError );
                t++;
            }
            if( TMath::Abs( e / iEnergy - 1. ) < 0.05 )
            {
                fData.back()->fEres = eres;
                fData.back()->fEresError = eresError;
            }
        }
    }
    is.close();
    
    cout << "Total number of data points: " << fData.size() << endl;
    
    TGraphErrors* fEres = new TGraphErrors( 1 );
    fEres->SetTitle( "" );
    fEres->SetMarkerStyle( 21 );
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fEres->SetPoint( i, fData[i]->fID, fData[i]->fEres );
        fEres->SetPointError( i, 0., fData[i]->fEresError );
    }
    fEres->Print();
    
    TCanvas* c = new TCanvas( "cE", "resolution vs MVA training ID", 10, 10, 1000, 600 );
    c->SetGridx( 0 );
    c->SetGridy( 0 );
    c->SetLeftMargin( 0.13 );
    c->Draw();
    
    fEres->Draw( "ap" );
    fEres->GetHistogram()->SetXTitle( "MVA ID" );
    stringstream ytitle;
    ytitle << "#Delta E/E (68\% containment) at " << iEnergy << " TeV";
    fEres->GetHistogram()->SetYTitle( ytitle.str().c_str() );
    fEres->GetHistogram()->GetYaxis()->SetTitleOffset( 1.4 );
    
    // plot those resolution values which are better than a certain value
    TCanvas* cM = new TCanvas( "cEM", "resolution", 100, 10, 1000, 600 );
    cM->SetGridx( 0 );
    cM->SetGridy( 0 );
    cM->SetLeftMargin( 0.13 );
    cM->Draw();
    
    TH1D* h = new TH1D( "eres", "", 100, -2., 2.5 );
    h->SetStats( 0 );
    h->SetXTitle( "log_{10} Energy/TeV" );
    h->SetYTitle( "#Delta E/E (68\% containment)" );
    h->GetYaxis()->SetTitleOffset( 1.4 );
    h->SetMaximum( iMax );
    h->Draw();
    
    z = 0;
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( fData[i]->fEres < EresMax )
        {
            cout << "###########################" << endl;
            cout << "MVA ID " << fData[i]->fID << "  " << fData[i]->fName << endl;
            fData[i]->fEresGraph->Print();
            fData[i]->fEresGraph->SetLineColor( z + 1 );
            fData[i]->fEresGraph->SetMarkerColor( z + 1 );
            fData[i]->fEresGraph->Draw( "lp" );
            z++;
        }
    }
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        if( iReferenceID == i )
        {
            cout << "Plotting reference ID " << i << endl;
            cout << "\t MVA ID " << fData[i]->fID << "  " << fData[i]->fName << endl;
            fData[i]->fEresGraph->SetLineWidth( 3. );
            fData[i]->fEresGraph->SetLineColor( 1 );
            fData[i]->fEresGraph->SetMarkerColor( 1. );
            fData[i]->fEresGraph->SetMarkerStyle( 20 );
            fData[i]->fEresGraph->Draw( "lp" );
        }
        else if( i == 0 )
        {
            fData[i]->fEresGraph->SetLineWidth( 3. );
            fData[i]->fEresGraph->SetLineColor( 2 );
            fData[i]->fEresGraph->SetMarkerColor( 2 );
            fData[i]->fEresGraph->SetMarkerStyle( 20 );
            fData[i]->fEresGraph->SetLineStyle( 2 );
            fData[i]->fEresGraph->Draw( "lp" );
        }
    }
    
    cM->Print( "dispMVA.pdf" );
}

