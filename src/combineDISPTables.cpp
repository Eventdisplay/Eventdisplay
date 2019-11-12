/*! \file combineDISPTables.cpp
    \brief merge different DISP tables into one file


*/

#include "VGlobalRunParameter.h"
#include "VDispTableReader.h"

#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


int main( int argc, char* argv[] )
{

    cout << endl;
    cout << "combineDISPTables " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "=======================================================" << endl;
    
    if( argc != 3 && argc != 4 )
    {
        cout << "combineDISPTables " << "<file list> <output root file> [options]" << endl;
        cout << endl;
        cout << endl;
        exit( 0 );
    }
    
    string fFileList = argv[1];
    string fOutFile  = argv[2];
    
    vector< string > fInputFile;
    
    ifstream is;
    is.open( fFileList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error: file list not found" << endl;
        exit( -1 );
    }
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            fInputFile.push_back( is_line );
        }
    }
    is.close();
    
    // new disp table (merger product)
    TFile* fTot = new TFile( fOutFile.c_str(), "RECREATE" );
    if( fTot->IsZombie() )
    {
        cout << "error creating output file: " << fOutFile << endl;
        exit( -1 );
    }
    VDispTableReader* fData = new VDispTableReader();
    fData->SetName( "dispTable" );
    fData->initialize( false );
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loop over all files
    for( unsigned int i = 0; i < fInputFile.size(); i++ )
    {
        cout << "now reading " << fInputFile[i] << endl;
        TFile f( fInputFile[i].c_str() );
        if( f.IsZombie() )
        {
            cout << "error reading input file " << fInputFile[i] << endl;
            exit( -1 );
        }
        // read disp table from this file
        VDispTableReader* a = ( VDispTableReader* )f.Get( "dispTable" );
        if( !a )
        {
            cout << "error: no disp table found in " << fInputFile[i] << endl;
            exit( -1 );
        }
        a->initialize( true );
        TTree* t = a->getTree();
        if( !t )
        {
            cout << "error: no tree in disp table in " << fInputFile[i] << endl;
            exit( -1 );
        }
        // loop over all file entries
        for( int j = 0; j < t->GetEntries(); j++ )
        {
            t->GetEntry( j );
            
            cout << "\t filling entry " << j << "\t" << a->ze << "\t" << a->az_min << "\t" << a->az_max << "\t" << a->woff << "\t" << a->pedvar << endl;
            
            fData->fill( a->ze, a->az_bin, a->az_min, a->az_max, a->woff, a->pedvar, ( TH2* )a->h2D_DispTable, ( TH2* )a->h2D_DispTableN, ( TH2* )a->h2D_DispPhiTable, ( TH2* )a->h2D_DispMissTable, ( TH3* )a->h3D_DispTable, ( TH3* )a->h3D_DispTableN, ( TH3* )a->h3D_DispPhiTable, ( TH3* )a->h3D_DispMissTable );
            
            if( j == 0 )
            {
                fData->fWidthScaleParameter = a->fWidthScaleParameter;
                fData->fLengthScaleParameter = a->fLengthScaleParameter;
            }
        }
        f.Close();
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fTot->cd();
    fData->Write();
    cout << "==================================================================================" << endl;
    cout << "==================================================================================" << endl;
    fData->print( true );
    fTot->Close();
}


