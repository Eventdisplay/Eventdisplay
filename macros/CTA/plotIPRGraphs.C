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
 * returns TelTypeID for a given telescope type
 *
 * telescope type ID changes with productions
 * loop here through all possibilities starting
 * from the most recent ones
 *
 */
ULong64_t find_tel_type( TFile *f, string iTelType, unsigned int i_sum_window )
{
    if( !f || f->IsZombie() )
    {
        return 0;
    }
    vector< ULong64_t > fTelType;
    if( iTelType == "LST" )
    {
        fTelType.push_back( 138704810 );
    }
    else if( iTelType == "FlashCam-MST" )
    {
        fTelType.push_back( 10408618 );
    }
    else if( iTelType == "NectarCam-MST" )
    {
        fTelType.push_back( 10408418 );
        fTelType.push_back( 10608418 );
    }
    else if( iTelType == "SC-MST" )
    {
        fTelType.push_back( 207308707 );
        fTelType.push_back( 205008707 );
    }
    else if( iTelType == "SST" )
    {
        fTelType.push_back( 201409917 );
        fTelType.push_back( 201109917 );
        fTelType.push_back( 201109916 );
    }
    else if( iTelType == "ASTRI-SST" )
    {
        fTelType.push_back( 201511619 );
    }
    else if( iTelType == "GC-SST" )
    {
        fTelType.push_back( 201309316 );
    }
    else if( iTelType == "DC-SST" )
    {
        fTelType.push_back( 909924 );
    }

    for( unsigned int i = 0; i < fTelType.size(); i++ )
    {
        std::ostringstream iGraphName;
        iGraphName << "IPRcharge_TelType" << fTelType[i] << "_SW" << i_sum_window;
        if( f->Get( iGraphName.str().c_str() ) )
        {
            return fTelType[i];
        }
    }
    return 0;
}

/*
 *  plot IPR graphs (from two or one file)
 *
*/
void plotIPRGraphs( string iFile1, string iFile2 = "", 
        float plot_min = 20., float plot_max = 9.e7,
        string Production = "PROD6" )
{
    // root files with IPR graphs (...ped.root)
    vector< string > fFile;
    fFile.push_back( iFile1 );
    if( iFile2.size() > 0 )
    {
        fFile.push_back( iFile2 );
    }
    
    // list of telescope types and the corresponding sum windows
    vector< unsigned int > fSumWindow;
    vector< string > fTelTypeName;
    
    fTelTypeName.push_back( "LST" );
    fSumWindow.push_back( 2 );
    
    fTelTypeName.push_back( "FlashCam-MST" );
    fSumWindow.push_back( 2 );
    fTelTypeName.push_back( "NectarCam-MST" );
    fSumWindow.push_back( 4 );
    fTelTypeName.push_back( "SC-MST" );
    fSumWindow.push_back( 4 );

    if( Production == "PROD3" )
    {
        fTelTypeName.push_back( "ASTRI-SST" );
        fSumWindow.push_back( 10 );
        fTelTypeName.push_back( "GC-SST" );
        fSumWindow.push_back( 10 );
        fTelTypeName.push_back( "DC-SST" );
        fSumWindow.push_back( 4 );
    }
    else
    {
        fTelTypeName.push_back( "SST" );
        fSumWindow.push_back( 6 );
    }
    
    vector< bool > fPlotted( fTelTypeName.size(), false );
    
    for( unsigned int i = 0; i < fTelTypeName.size(); i++ )
    {
        
        std::ostringstream iCanvasName;
        iCanvasName << "cIPRcharge_TelType" << fTelTypeName[i] << "_SW" << fSumWindow[i];
        
        std::ostringstream iCanvasTitle;
        iCanvasTitle << "IPRcharge graph for " << fTelTypeName[i] << " (summation window " << fSumWindow[i] << ")";
        
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
            ULong64_t i_tel_type = find_tel_type( fF, fTelTypeName[i], fSumWindow[i] );
            std::ostringstream iGraphName;
            iGraphName << "IPRcharge_TelType" << i_tel_type << "_SW" << fSumWindow[i];
            std::ostringstream iGraphTitle;
            iGraphTitle << fTelTypeName[i] << " (" << i_tel_type << ", summation window " << fSumWindow[i] << ")";
            
            TGraphErrors* ipr = ( TGraphErrors* )fF->Get( iGraphName.str().c_str() );
            if( ipr )
            {
/*                if( i_tel_type == 201109916 )
                {
                    double x, y;
                    for( int p = 0; p < ipr->GetN(); p++ )
                    {
                        ipr->GetPoint( p, x, y );
                        ipr->SetPoint( p, x*6.74/2.77, y );
                        ipr->GetPoint( p, x, y );
                    }
                } */
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
        iPrintName << "IPRComparison_TelType_" << fTelTypeName[i] << ".pdf";
        c->Print( iPrintName.str().c_str() );
    }
}
