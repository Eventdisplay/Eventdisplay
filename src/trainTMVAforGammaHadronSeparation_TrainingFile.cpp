/*  \file  trainTMVAforGammaHadronSeparation_TrainingFile.cpp
    \brief make background training files for TMVA using MC and randomize theta2

    Note: no wobble offsets are implemented yet

*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

/*

   get list of file names (full paths) from parameter file

*/
vector< string > getListofFiles( string iFileName )
{
    ifstream is;
    is.open( iFileName.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, no input file list found: " << iFileName << endl;
        exit( -1 );
        
    }
    cout << "list of files from " << iFileName << endl;
    cout << endl;
    string is_line;
    string temp;
    vector< string > iFileList;
    
    while( getline( is, is_line ) )
    {
        stringstream is_stream( is_line );
        is_stream >> temp;
        if( temp != "*" )
        {
            continue;
        }
        is_stream >> temp;
        if( temp != "BACKGROUNDFILE" )
        {
            continue;
        }
        is_stream >> temp;
        iFileList.push_back( temp );
    }
    is.close();
    
    return iFileList;
}


/*

   main...

*/
int main( int argc, char* argv[] )
{
    if( argc != 5 )
    {
        cout << endl;
        cout << "trainTMVAforGammaHadronSeparation_TrainingFile <run parameter file> <output directory for mscw_energy file> <maximum theta2 value [deg^2]> <seed>" << endl;
        cout << endl;
        cout << "<run parameter file>";
        cout << "parameter file with list of background mscw file (same format as in trainTMVAforGammaHadronSeparation)" << endl;
        cout << "<seed>                 seed for random generator (default: 0)" << endl;
        cout << endl;
        cout << "Note: no wobble offsets are implemented yet" << endl;
        cout << "Note: you might want to have some additional cuts (e.g. distance to the camera center) hardwired" << endl;
        cout << endl;
        exit( 0 );
    }
    
    cout << endl;
    cout << "trainTMVAforGammaHadronSeparation_TrainingFile" << endl;
    cout << "====================================" << endl;
    cout << endl;
    
    double fMaxTheta = TMath::Sqrt( atof( argv[3] ) );
    UInt_t fRandomSeed = ( UInt_t )( atoi( argv[4] ) );
    
    cout << "replacing Xoff/Yoff by random value inside a circle of radius " << fMaxTheta << " deg" << endl;
    cout << "replace MCxoff/MCyoff by 0,0 " << endl;
    
    // random generator
    TRandom3* fRandom = new TRandom3( fRandomSeed );
    
    // get list of files to change
    vector< string > iFileList = getListofFiles( argv[1] );
    string i_outputDir = argv[2];
    
    // loop over all input files
    for( unsigned int f = 0; f < iFileList.size(); f++ )
    {
    
        // input data file
        string iIfile = iFileList[f];
        
        TFile* fInput = new TFile( iIfile.c_str() );
        if( fInput->IsZombie() )
        {
            cout << "error reading input file " << iIfile << endl;
            cout << "exit..." << endl;
            exit( -1 );
        }
        // input data tree
        TTree* fInputData = ( TTree* )fInput->Get( "data" );
        if( !fInputData )
        {
            cout << "error reading input tree from " << iIfile << endl;
            exit( -1 );
        }
        double fXoff = 0.;
        double fYoff = 0.;
        double fMCXoff = 0.;
        double fMCYoff = 0.;
        fInputData->SetBranchAddress( "Xoff", &fXoff );
        fInputData->SetBranchAddress( "Yoff", &fYoff );
        fInputData->SetBranchAddress( "MCxoff", &fMCXoff );
        fInputData->SetBranchAddress( "MCyoff", &fMCYoff );
        Int_t nentries = ( Int_t )fInputData->GetEntries();
        
        // output data file
        string iOfile = i_outputDir + "/" + gSystem->BaseName( iFileList[f].c_str() );
        if( iIfile == iOfile )
        {
            cout << "error: input and output file is identical" << endl;
            exit( -1 );
        }
        TFile* fOutput = new TFile( iOfile.c_str(), "RECREATE" );
        if( fOutput->IsZombie() )
        {
            cout << "error opening output file " << iOfile << endl;
            exit( -1 );
        }
        
        // output data tree
        TTree* fOutputData = fInputData->CloneTree( 0 );
        
        double i_theta = 0.;
        double i_r = 0.;
        
        // loop over all entries
        cout << "Looping over " << nentries << " events in " << fInput->GetName() << endl;
        for( Int_t i = 0; i < nentries; i++ )
        {
            fInputData->GetEntry( i );
            
            i_r = TMath::Sqrt( fRandom->Uniform( 1. ) );
            i_theta = fRandom->Uniform( TMath::TwoPi() );
            fXoff = fMaxTheta * i_r * TMath::Cos( i_theta );
            fYoff = fMaxTheta * i_r * TMath::Sin( i_theta );
            
            fMCXoff = 0.;
            fMCYoff = 0.;
            
            fOutputData->Fill();
        }
        cout << "writing new tree to " << fOutput->GetName() << endl;
        fOutputData->AutoSave();
        
        delete fInput;
        delete fOutput;
    }  // end of loop over all files
}
