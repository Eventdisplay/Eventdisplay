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

#include "VGlobalRunParameter.h"
#include "CData.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "VGammaHadronCuts.h"
#include "VDL2Writer.h"
#include "VInstrumentResponseFunctionRunParameter.h"
#include "VMonteCarloRunHeader.h"
#include "VTableLookupRunParameter.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

VEffectiveAreaCalculatorMCHistograms* copyMCHistograms( TChain* c );

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
        cout << "./fillDL2Trees" << endl;
        exit( EXIT_SUCCESS );
    }
    string fOutputfileName = argv[2];
    
    /////////////////////////////////////////////////////////////////
    // read run parameters from file
    VInstrumentResponseFunctionRunParameter* fRunPara = new VInstrumentResponseFunctionRunParameter();
    fRunPara->SetName( "fillDL2Tree_runparameter" );
    if( !fRunPara->readRunParameterFromTextFile( argv[1] ) )
    {
        cout << "error reading runparameters from text file" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fRunPara->print();
    
    /////////////////////////////////////////////////////////////////
    // open output file and write results to dist
    TFile* fOutputfile = new TFile( fOutputfileName.c_str(), "RECREATE" );
    if( fOutputfile->IsZombie() )
    {
        cout << "Error in opening output file: " << fOutputfile->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    /////////////////////////////////////////////////////////////////
    // gamma/hadron cuts
    VGammaHadronCuts* fCuts = new VGammaHadronCuts();
    fCuts->initialize();
    fCuts->setNTel( fRunPara->telconfig_ntel, 
                    fRunPara->telconfig_arraycentre_X, 
                    fRunPara->telconfig_arraycentre_Y );
    fCuts->setInstrumentEpoch( fRunPara->getInstrumentEpoch( true ) );
    fCuts->setTelToAnalyze( fRunPara->fTelToAnalyse );
    fCuts->setReconstructionType( fRunPara->fReconstructionType );
    if( !fCuts->readCuts( fRunPara->fCutFileName, 2 ) )
    {
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE ) ;
    }
    fRunPara->fGammaHadronCutSelector = fCuts->getGammaHadronCutSelector();
    fRunPara->fDirectionCutSelector   = fCuts->getDirectionCutSelector();
    fCuts->initializeCuts( -1, fRunPara->fGammaHadronProbabilityFile );
    fCuts->printCutSummary();
    
    /////////////////////////////////////////////////////////////////
    // read MC header (might not be there, no problem; but depend on right input in runparameter file)
    VMonteCarloRunHeader* iMonteCarloHeader = fRunPara->readMCRunHeader();
    
    /////////////////////////////////////////////////////////////////////////////
    // DL2 writer
    VDL2Writer fDL2Writer( fRunPara, fCuts );
    
    /////////////////////////////////////////////////////////////////////////////
    // set effective area Monte Carlo histogram class
    VEffectiveAreaCalculatorMCHistograms* fMC_histo = 0;
    
    /////////////////////////////////////////////////////////////////////////////
    // load data chain
    TChain* c = new TChain( "data" );
    if( !c->Add( fRunPara->fdatafile.c_str(), -1 ) )
    {
        cout << "Error while trying to add mscw data tree from file " << fRunPara->fdatafile  << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    CData d( c, true, true );
    fCuts->setDataTree( &d );
    d.setReconstructionType( fCuts->fReconstructionType );
    
    //////////////////////////////////////////////////////////////////////////////
    // MC histograms
    if( fRunPara->fFillingMode != 1 && fRunPara->fFillingMode != 2 )
    {
        fMC_histo = copyMCHistograms( c );
        if( fMC_histo )
        {
            fMC_histo->matchDataVectors( fRunPara->fAzMin, 
                                         fRunPara->fAzMax, 
                                         fRunPara->fSpectralIndex );
            fMC_histo->print();
        }
        else
        {
            cout << "Warning: failed reading MC histograms" << endl;
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////
    // fill DL2 trees
    fOutputfile->cd();
    fDL2Writer.fill( &d, fRunPara->fEnergyReconstructionMethod );
    
    
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
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
        if( fCuts->getMVACutGraphs().size() > 0 )
        {
              cout << "writing MVA cut graphs to " << fOutputfile->GetName() << endl;
              fOutputfile->cd();
              if( fOutputfile->mkdir( "mvaGraphs" ) )
              {
                  fOutputfile->cd( "mvaGraphs" );   
                  for( unsigned int i = 0; i < fCuts->getMVACutGraphs().size(); i++ )
                  {
                      if( fCuts->getMVACutGraphs()[i] )
                      {
                          fCuts->getMVACutGraphs()[i]->Write();
                      }
                  }
              }
              fOutputfile->cd();
        }
    // writing cuts to disk
    if( fCuts )
    {
        fCuts->terminate( ( fRunPara->fEffArea_short_writing || fRunPara->fFillingMode == 3 ) );
    }
    // writing monte carlo header to disk
    if( iMonteCarloHeader )
    {
        iMonteCarloHeader->Write();
    }
    // write run parameters to disk
    if( fRunPara )
    {
        fRunPara->Write();
    }
    
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


