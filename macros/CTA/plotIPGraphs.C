/*  plot IPR graphs for all CTA telescope types
 *
 *  compare results from two IPR files
 *
 *  use with .L plotIPR.C++
 *
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"

/*
 *  plot IPR graphs (from two or one file)
 *
*/
void plotIPRGraphs( string iFile1, string iFile2 = "", float plot_min = 20., float plot_max = 9.e7 )
{
    // root files with IPR graphs (...ped.root)
    vector< string > fFile;
    fFile.push_back( iFile1 );
    if( iFile2.size() > 0 )
    {
        fFile.push_back( iFile2 );
    }
    
    // list of telescope types and the corresponding sum windows
    vector< ULong64_t > fTelType;
    vector< unsigned int > fSumWindow;
    vector< string > fTelTypeName;
    
    fTelTypeName.push_back( "LST" );
    fTelType.push_back( 138704810 );
    fSumWindow.push_back( 4 );
    
    fTelTypeName.push_back( "FlashCam MST" );
    fTelType.push_back( 10408618 );
    fSumWindow.push_back( 2 );
    fTelTypeName.push_back( "NectarCam MST" );
    fTelType.push_back( 10408418 );
    fSumWindow.push_back( 6 );
    fTelTypeName.push_back( "SC-MST" );
    fTelType.push_back( 207308707 );
    fSumWindow.push_back( 6 );
    
    fTelTypeName.push_back( "ASTRI SST" );
    fTelType.push_back( 201511619 );
    fSumWindow.push_back( 10 );
    fTelTypeName.push_back( "GC-SST" );
    fTelType.push_back( 201309316 );
    fSumWindow.push_back( 10 );
    fTelTypeName.push_back( "DC-SST" );
    fTelType.push_back( 909924 );
    fSumWindow.push_back( 4 );
    
    vector< bool > fPlotted( fTelType.size(), false );
    
    
    for( unsigned int i = 0; i < fTelType.size(); i++ )
    {
        std::ostringstream iGraphName;
        iGraphName << "IPRcharge_TelType" << fTelType[i] << "_SW" << fSumWindow[i];
        
        std::ostringstream iGraphTitle;
        iGraphTitle << fTelTypeName[i] << " (" << fTelType[i] << ", summation window " << fSumWindow[i] << ")";
        
        std::ostringstream iCanvasName;
        iCanvasName << "cIPRcharge_TelType" << fTelType[i] << "_SW" << fSumWindow[i];
        
        std::ostringstream iCanvasTitle;
        iCanvasTitle << "IPRcharge graph for " << fTelType[i] << " (summation window " << fSumWindow[i] << ")";
        
        TCanvas* c = new TCanvas( iCanvasName.str().c_str(), iCanvasTitle.str().c_str(), 50 + i * 100, 50 + i * 10, 600, 400 );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        c->SetLogy( 1 );
        c->Draw();
        
        for( unsigned int f = 0; f < fFile.size(); f++ )
        {
            TFile* fF = new TFile( fFile[f].c_str() );
            if( fF->IsZombie() )
            {
                continue;
            }
            TGraphErrors* ipr = ( TGraphErrors* )fF->Get( iGraphName.str().c_str() );
            if( ipr )
            {
                ipr->SetLineColor( f + 1 );
                ipr->SetLineWidth( 2 );
                ipr->SetMinimum( plot_min );
                ipr->SetMaximum( plot_max );
                ipr->SetTitle( iGraphTitle.str().c_str() );
                if( !fPlotted[i] )
                {
                    ipr->Draw();
                    fPlotted[i] = true;
                }
                else
                {
                    ipr->Draw( "same" );
                }
            }
        }
        
        // print everything
        std::ostringstream iPrintName;
        iPrintName << "IPRComparison_TelType_" << fTelType[i] << ".pdf";
        c->Print( iPrintName.str().c_str() );
    }
}


/*
 * merge IPR graphs from several files
 *
 * paths need to be added by hand
 *
 */
void mergeIPRGraphs( string iOutFile )
{

    // list of telescope types and the corresponding sum windows
    vector< ULong64_t > fTelType;
    vector< unsigned int > fSumWindow_min;
    vector< unsigned int > fSumWindow_max;
    vector< string > fTelTypeName;
    vector< string > fIPRFileName;
    
    fTelTypeName.push_back( "LST" );
    fTelType.push_back( 138704810 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 10 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20160603/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranal.ped.root" );
    
    fTelTypeName.push_back( "FlashCam MST" );
    fTelType.push_back( 10408618 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 10 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20160603/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranal.ped.root" );
    
    fTelTypeName.push_back( "NectarCam MST" );
    fTelType.push_back( 10408418 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 10 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20160603/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranal.ped.root" );
    
    fTelTypeName.push_back( "GC-SST" );
    fTelType.push_back( 201309316 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 10 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20160603/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranal.ped.root" );
    
    fTelTypeName.push_back( "DC-SST" );
    fTelType.push_back( 909924 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 10 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20160603/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranal.ped.root" );
    
    fTelTypeName.push_back( "ASTRI SST" );
    fTelType.push_back( 201511619 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 20 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20171029/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranalv02120ns.ped.root" );
    
    
    fTelTypeName.push_back( "SC-MST" );
    fTelType.push_back( 207308707 );
    fSumWindow_min.push_back( 1 );
    fSumWindow_max.push_back( 10 );
    fIPRFileName.push_back( "$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.Paranal-20171129/gamma_20deg_180deg_run5___cta-prod3-demo_desert-2150m-Paranal.SCT.ped.root" );
    
    
    // output (merged file)
    TFile* fO = new TFile( iOutFile.c_str(), "recreate" );
    if( fO->IsZombie() )
    {
        cout << "Error opening output file " << iOutFile << endl;
        return;
    }
    cout << "Output (merged) file: " << iOutFile << endl;
    
    /////////////////////////////////////////////////////////////
    // loop over all telescope types and copy files over
    for( unsigned int i = 0; i < fTelType.size(); i++ )
    {
    
        // source file
        TFile* iS = new TFile( fIPRFileName[i].c_str() );
        if( iS->IsZombie() )
        {
            cout << "Error opening source file " << fIPRFileName[i].c_str() << endl;
            continue;
        }
        cout << "Reading " << fTelType[i] << " from " <<  fIPRFileName[i].c_str() << endl;
        for( unsigned int sw = fSumWindow_min[i]; sw <= fSumWindow_max[i]; sw++ )
        {
            // ipr graph
            std::ostringstream iGraphName;
            iGraphName << "IPRcharge_TelType" << fTelType[i] << "_SW" << sw;
            TGraphErrors* ipr = ( TGraphErrors* )iS->Get( iGraphName.str().c_str() );
            if( ipr )
            {
                fO->cd();
                ipr->Write();
            }
        }
        // pedestal tree
        iS->cd();
        std::ostringstream iTreenName;
        iTreenName << "tPeds_" << fTelType[i];
        TTree* itp = ( TTree* )iS->Get( iTreenName.str().c_str() );
        if( itp )
        {
            fO->cd();
            TTree* itpN = ( TTree* )itp->CloneTree();
            itpN->Write();
        }
        iS->Close();
    }
}


