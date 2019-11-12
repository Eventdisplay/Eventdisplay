/*!  \file printDISPTables.cpp
     \brief read tables for disp angular reconstruction

     ** PRELIMINARY **


*/

#include "VGlobalRunParameter.h"
#include "TFile.h"

#include "VDispTableReader.h"

#include <iostream>
#include <string>

using namespace std;

int main( int argc, char* argv[] )
{
    cout << endl;
    cout << "printDISPTables " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "=======================================================" << endl;
    if( argc != 2 )
    {
        cout << "printDISPTables <table file> " << endl;
        exit( 0 );
    }
    TFile* f = new TFile( argv[1] );
    if( f->IsZombie() )
    {
        exit( 0 );
    }
    
    VDispTableReader* fData = ( VDispTableReader* )f->Get( "dispTable" );
    if( fData )
    {
        fData->initialize( true );
        fData->print( true );
    }
    f->Close();
}
