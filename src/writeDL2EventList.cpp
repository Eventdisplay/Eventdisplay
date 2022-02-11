/*! \file writeDL2 eventlist
 * 
 * write tree with DL2-level data:
 * reconstructed events including g/h separation parameters
 *  
 * input is a file list of mscw_energy output file from simulations
 * or from data
 *
 */

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TTree.h"

#include "CData.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "VDL2Writer.h"
#include "VGlobalRunParameter.h"
#include "VMonteCarloRunHeader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

VEffectiveAreaCalculatorMCHistograms* copyMCHistograms( TChain* c );

void printHelp( int argc, char* argv[] )
{
    // print version only
    if( argc == 2 )
    {
        string fCommandLine = argv[1];
        if( fCommandLine == "-v" || fCommandLine == "--version" )
        {
            VGlobalRunParameter iRunPara;
            cout << iRunPara.getEVNDISP_VERSION() << endl;
            exit( EXIT_SUCCESS );
        }
    }
    cout << endl;
    cout << "writeDL2EventList " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "-----------------------------" << endl;
    
    /////////////////////////////////////////////////////////////////
    // read command line parameters
    if( argc != 3 )
    {
        cout << "DL2 event list writer" << endl;
        cout << endl;
        cout << "./writeDL2EventList <config file> <output file>" << endl;
        exit( EXIT_SUCCESS );
    }
}

bool writeAdditionalObjects( TChain *c, 
                       TFile *iOutFile )
{
    if( !c || !iOutFile ) return false;

    // writing Monte Carlo header to disk
    TFile* iF = (TFile*)c->GetFile();
    if( !iF ) return false;
    VMonteCarloRunHeader* iMC = (VMonteCarloRunHeader*)iF->Get( "MC_runheader" );
    TTree *telconfig = (TTree*)iF->Get( "telconfig" );
    iOutFile->cd();
    if( iMC ) iMC->Write();
    if( telconfig ) telconfig->CloneTree()->Write();
    return true;
}


void fillDL2EventList( string fConfigFile,
                       TFile *fOutputfile )
{
    // DL2 writer
    VDL2Writer fDL2Writer( fConfigFile );
    
    // load input data chain
    TChain* c = new TChain( "data" );
    for( unsigned int i = 0; i < fDL2Writer.getDataFile().size(); i++ )
    {
        cout << "Addding " << fDL2Writer.getDataFile()[i] << endl;
        if( !c->Add( fDL2Writer.getDataFile()[i].c_str(), -1 ) )
        {
            cout << "Error while trying to add mscw data tree from file ";
            cout << fDL2Writer.getDataFile()[i] << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
    }
    CData d( c, true, true );
    
    // MC histograms
    copyMCHistograms( c );

    // fill DL2 trees
    fOutputfile->cd();
    fDL2Writer.fill( &d );
    
    // write results to disk
    fDL2Writer.writeDataTrees( fOutputfile, c );

    // copy additional useful objects
    writeAdditionalObjects( c, fOutputfile );
}
    

//////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{
    printHelp( argc, argv );
    
    string fConfigFile = argv[1];
    string fOutputfileName = argv[2];

    // output root file with DL2 event list
    TFile* fOutputfile = new TFile( fOutputfileName.c_str(), "RECREATE" );
    if( fOutputfile->IsZombie() )
    {
        cout << "Error in opening output file: " << fOutputfile->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }

    fillDL2EventList( fConfigFile, 
                      fOutputfile );
    
    fOutputfile->Close();
    cout << "end..." << endl;
}

/*
 * return MC histograms from mscw_energy files
 *
 * - loop over all mscw_energy files
 * - add MC histograms
 */
VEffectiveAreaCalculatorMCHistograms* copyMCHistograms( TChain* c )
{
    VEffectiveAreaCalculatorMCHistograms* iMC_his = 0;
    if( c )
    {
        // loop over all files and add MC histograms
        TObjArray* fileElements = c->GetListOfFiles();
        TChainElement* chEl = 0;
        TIter next( fileElements );
        unsigned int z = 0;
        while( ( chEl = ( TChainElement* )next() ) )
        {
            TFile* ifInput = new TFile( chEl->GetTitle() );
            if( !ifInput->IsZombie() )
            {
                if( z == 0 )
                {
                    iMC_his = ( VEffectiveAreaCalculatorMCHistograms* )ifInput->Get( "MChistos" );
                }
                else
                {
                    if( iMC_his )
                    {
                        iMC_his->add( ( VEffectiveAreaCalculatorMCHistograms* )ifInput->Get( "MChistos" ) );
                    }
                    ifInput->Close();
                }
                z++;
            }
        }
    }
    return iMC_his;
}

