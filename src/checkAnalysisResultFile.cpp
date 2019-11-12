/*
 *  check if an analysis file was correctly closed and
 *  that it contains all trees
*/

#include "TFile.h"
#include "TTree.h"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

int main( int argc, char* argv[] )
{
    cout << argc << endl;
    if( argc != 2 && argc != 3 )
    {
        cout << "./checkAnalysisResultFile <root file> [ntel]" << endl;
        cout << endl;
        exit( -1 );
    }
    string iFile = argv[1];
    int ntel = 4;
    if( argc == 3 )
    {
        ntel = atoi( argv[2] );
    }
    
    TFile iF( iFile.c_str() );
    if( iF.IsZombie() )
    {
        exit( EXIT_FAILURE );
    }
    TTree* s = ( TTree* )iF.Get( "showerpars" );
    if( !s )
    {
        cout << "missing showerpars tree in " << iFile << endl;
        exit( EXIT_FAILURE );
    }
    char hname[200];
    for( int i = 0; i < ntel; i++ )
    {
        sprintf( hname, "Tel_%d/tpars", i + 1 );
        TTree* t = ( TTree* )iF.Get( hname );
        if( !t )
        {
            cout << "missing " << hname << " tree in " << iFile << endl;
            exit( EXIT_FAILURE );
        }
    }
    vector< string > i_objects;
    i_objects.push_back( "runparameterV2" );
    i_objects.push_back( "telconfig" );
    i_objects.push_back( "MC_runheader" );
    i_objects.push_back( "MChistos" );
    i_objects.push_back( "EvndispReconstructionParameter" );
    i_objects.push_back( "pointingData" );
    for( unsigned int i = 0; i < i_objects.size(); i++ )
    {
        TNamed* r = ( TNamed* )iF.Get( i_objects[i].c_str() );
        if( !r )
        {
            cout << "missing object " << i_objects[i] << " in " << iFile << endl;
            exit( EXIT_FAILURE );
        }
    }
    if( iF.TestBit( TFile::kRecovered ) )
    {
        cout << "file has been recoverd " << iFile << endl;
    }
    iF.Close();
    
    cout << "file sucess " << iFile << endl;
    exit( EXIT_SUCCESS );
}
