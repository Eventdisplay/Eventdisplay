#include <iostream>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"

using namespace std;

int main( int argc, char* argv[] )
{
    if( argc != 3 )
    {
        printf( "Error, must give input filename as first argument and the treename as second argument, exiting...\n" );
        exit( 1 );
    }
    
    char inputfile[200] = "" ;
    sprintf( inputfile, "%s", argv[1] ) ;
    
    char treename[200] = "" ;
    sprintf( treename, "%s/%s", inputfile, argv[2] ) ;
    
    TFile* inpfile = new TFile( inputfile, "READ" ) ;
    TChain* dataChain = new TChain() ;
    
    dataChain->Add( treename ) ;
    printf( "file %s\n", inputfile ) ;
    printf( "treename %s\n", argv[2] ) ;
    printf( "nevents %d\n", ( int )( dataChain->GetEntries() ) ) ;
    
    int eventNumber = 0 ;
    dataChain->SetBranchAddress( "eventNumber", &eventNumber ) ;
    dataChain->GetEntry( dataChain->GetEntries() - 1 ) ;
    int maxeventnumber = eventNumber ;
    
    printf( "maxeventnumber %d\n", maxeventnumber ) ;
    
    inpfile->Close();
}
