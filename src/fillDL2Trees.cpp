/*! \file  fill DL2 trees
 *  
 *  applying gamma / hadron separation (e.g., TMVA)
 *
 *   input is a file list of mscw_energy output file from gamma-ray simulations
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
bool writeMCRunHeader( TChain *c, TFile *iOutFile );

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
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
    cout << "fillDL2Trees " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    
    /////////////////////////////////////////////////////////////////
    // read command line parameters
    if( argc != 3 )
    {
        cout << endl;
        cout << "./fillDL2Trees <config file> <output file>" << endl;
        exit( EXIT_SUCCESS );
    }
    string fConfigFile = argv[1];
    string fOutputfileName = argv[2];

    // open output file and write results to dist
    TFile* fOutputfile = new TFile( fOutputfileName.c_str(), "RECREATE" );
    if( fOutputfile->IsZombie() )
    {
        cout << "Error in opening output file: " << fOutputfile->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    // DL2 writer
    VDL2Writer fDL2Writer( fConfigFile );
    
    // load data chain
    TChain* c = new TChain( "data" );
    if( !c->Add( fDL2Writer.getDataFile().c_str(), -1 ) )
    {
        cout << "Error while trying to add mscw data tree from file ";
        cout << fDL2Writer.getDataFile() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    CData d( c, true, true );
    
    // MC histograms
    VEffectiveAreaCalculatorMCHistograms *fMC_histo = copyMCHistograms( c );
    if( fMC_histo )
    {
        fMC_histo->print();
    }
    else
    {
        cout << "Warning: failed reading MC histograms" << endl;
    }
    
    // fill DL2 trees
    fOutputfile->cd();
    fDL2Writer.fill( &d );
    
    // write results to disk
    if( fDL2Writer.getEventCutDataTree() )
    {
        cout << "writing event data trees: (";
        cout << fDL2Writer.getEventCutDataTree()->GetName();
        cout << ") to " << fOutputfile->GetName() << endl;
        fOutputfile->cd();
        if( fDL2Writer.getEventCutDataTree() )
        {
            fDL2Writer.getEventCutDataTree()->Write();
        }
        if( c )
        {
            c->Merge(fOutputfile, 0, "keep" );
        }
    }
    writeMCRunHeader( c, fOutputfile );
    
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

bool writeMCRunHeader( TChain *c, 
                       TFile *iOutFile )
{
    if( !c || !iOutFile ) return false;

    // writing monte carlo header to disk
    TFile* iF = (TFile*)c->GetFile();
    if( !iF ) return false;
    VMonteCarloRunHeader* iMC = (VMonteCarloRunHeader*)iF->Get( "MC_runheader" );
    iOutFile->cd();
    if( iMC )
    {
        iMC->Write();
        return true;
    }
    return false;
}

