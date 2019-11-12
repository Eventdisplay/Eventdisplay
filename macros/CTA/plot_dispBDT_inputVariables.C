/*
 * plot input variables for disp training
 * (telescopes types are set for CTA,
 *  change for other telescopes)
 *
 *  reads data trees produced
 *  for BDT training from
 *  trainTMVAforAngularReconstruction
 *
 *  execute with .L plot_dispBDT_inputVariables.C++
 *
*/

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void plot_dispBDT_inputVariables( string iDirectoryWithTrees = "./" )
{
    vector< string > fTelType;
    fTelType.push_back( "138704810" ); //  LST
    fTelType.push_back( "10408618" ); //   FlashCam MST
    fTelType.push_back( "10408418" ); //   NectarCam MST
    fTelType.push_back( "201511619" ); //  ASTRI SST
    fTelType.push_back( "201309316" ); //  GC-SST
    fTelType.push_back( "909924" ); //     DC-SST
    fTelType.push_back( "207308707" ); //  SC-MST
    
    // data trees
    vector< TTree* > fData;
    unsigned int z = 0;
    for( unsigned int t = 0; t < fTelType.size(); t++ )
    {
        ostringstream iFileName;
        iFileName << iDirectoryWithTrees << "/BDTDispEnergy_" << fTelType[t] << ".root";
        TFile* f = new TFile( iFileName.str().c_str() );
        if( f->IsZombie() )
        {
            continue;
        }
        ostringstream iTreeName;
        iTreeName << "dispTree_" << fTelType[t];
        if( f->Get( iTreeName.str().c_str() ) )
        {
            fData.push_back( ( TTree* )f->Get( iTreeName.str().c_str() ) );
            fData.back()->SetLineColor( z + 1 );
            fData.back()->SetLineWidth( 2 );
            cout << "Telescope type " << fTelType[t] << ": color " << z + 1 << endl;
            z++;
        }
    }
    cout << "Found " << z << " data trees" << endl;
    
    
    vector< string > fVarType;
    vector< float > fVarMin;
    vector< float > fVarMax;
    fVarType.push_back( "width" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 2. );
    fVarType.push_back( "length" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 10. );
    fVarType.push_back( "wol" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 1. );
    fVarType.push_back( "size" );
    fVarMin.push_back( 2. );
    fVarMax.push_back( 10. );
    fVarType.push_back( "tgrad_x*tgrad_x" );
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
    fVarType.push_back( "EHeight" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 100. );
    fVarType.push_back( "Rcore" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 3500. );
    fVarType.push_back( "disp" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 20. );
    fVarType.push_back( "dispEnergy" );
    fVarMin.push_back( -1. );
    fVarMax.push_back( 1. );
    fVarType.push_back( "dispCore" );
    fVarMin.push_back( 0. );
    fVarMax.push_back( 3500. );
    
    // loop over all variables and plot them
    for( unsigned int v = 0; v < fVarType.size(); v++ )
    {
        cout << "Filling and plotting " << fVarType[v] << endl;
        ostringstream iCanvasName;
        iCanvasName << "c" << fVarType[v];
        TCanvas* c = new TCanvas( iCanvasName.str().c_str(), fVarType[v].c_str(), 10, 10 + 10 * v, 600, 600 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetLogy( 1 );
        c->Draw();
        
        for( unsigned int t = 0; t < fData.size(); t++ )
        {
            string iDrawOption = "";
            if( t != 0 )
            {
                iDrawOption = "same";
            }
            ostringstream iHName;
            iHName << "h" << fData[t]->GetName() << "_" << fVarType[v];
            TH1D* h = new TH1D( iHName.str().c_str(), "", 100, fVarMin[v], fVarMax[v] );
            h->SetLineColor( fData[t]->GetLineColor() );
            h->SetLineWidth( 2 );
            h->SetStats( 0 );
            h->SetXTitle( fVarType[v].c_str() );
            
            fData[t]->Project( iHName.str().c_str(), fVarType[v].c_str(),
                               "size>1.&&ntubes>4.&&length>0.&&loss<0.2" );
            h->Draw( iDrawOption.c_str() );
            cout << "\t" << iHName.str() << endl;
            cout << "\t\t Overflow: " << h->GetBinContent( h->GetNbinsX() + 1 ) << endl;
            cout << "\t\t Underflow: " << h->GetBinContent( 0 ) << endl;
        }
        ostringstream iPrintName;
        iPrintName << fVarType[v] << ".pdf";
        c->Print( iPrintName.str().c_str() );
    }
}


