/*  calculate MVA value and write a new tree with one entry
 *  per input event
 *
 *  input data is a mscw root file
 *
 *  number of energy and ze bins fixed for typical
 *  CTA values
 *
 *
 */


#include "TChain.h"
#include "TTree.h"

#include "CData.h"
#include "VEvndispRunParameter.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "VTMVAEvaluator.h"

#include <iostream>
#include <string>

using namespace std;

///////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
    if( argc != 4 )
    {
        cout << "./writeEventListTMVA <input mscw file> <output root file> <BDT location>" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    cout << "./writeEventListTMVA" << endl;
    
    
    // input MSCW file
    string fInputMSWFile = argv[1];
    
    // output file
    string fOutputFile = argv[2];
    
    // TMVA files
    string fTMVAFile = argv[3];
    
    // (end of user input)
    
    ////////////////////////////////////////////////////////
    // (hardwired values)
    unsigned int fTMVAWeightFileIndex_Emin = 0;
    unsigned int fTMVAWeightFileIndex_Emax = 8;
    unsigned int fTMVAWeightFileIndex_Zmin = 0;
    unsigned int fTMVAWeightFileIndex_Zmax = 0;
    double fTMVAEnergy_StepSize = 0.2;
    string fInstrumentEpoch = "CTA";
    
    /////////////////////////////////////////
    // open input data files
    TChain* cD = new TChain( "data" );
    if( !cD->Add( fInputMSWFile.c_str() ) )
    {
        cout << "Error opening input file " << fInputMSWFile << endl;
        exit( EXIT_FAILURE );
    }
    CData* fData = new CData( cD, true, true );
    cout << "Reading data from " << fInputMSWFile << endl;
    cout << "Total number of events: " << fData->fChain->GetEntries() << endl;
    
    /////////////////////////////////////////
    // open output data file
    TFile* fOutputData = new TFile( fOutputFile.c_str(), "RECREATE" );
    if( fOutputData->IsZombie() )
    {
        cout << "Error creating output file " << fOutputFile << endl;
        exit( EXIT_FAILURE );
    }
    TTree* fTMVA = new TTree( "BDTEventList", "BDT results" );
    double BDT = 0.;
    fTMVA->Branch( "BDT", &BDT, "BDT/D" );
    
    VTMVAEvaluator* fTMVAEvaluator = new VTMVAEvaluator();
    if( !fTMVAEvaluator->initializeWeightFiles( fTMVAFile,
            fTMVAWeightFileIndex_Emin, fTMVAWeightFileIndex_Emax,
            fTMVAWeightFileIndex_Zmin, fTMVAWeightFileIndex_Zmax,
            fTMVAEnergy_StepSize, fInstrumentEpoch ) )
    {
        cout << "Error initializing weight files" << endl;
        exit( EXIT_FAILURE );
    }
    
    fTMVAEvaluator->initializeDataStrutures( fData );
    
    /////////////////////////////////////////
    // loop over all events
    int nentries = fData->fChain->GetEntries();
    cout << endl << endl;
    cout << "loop over " << nentries << " entries " << endl;
    for( int i = 0; i < nentries; i++ )
    {
        fData->GetEntry( i );
        
        fTMVAEvaluator->evaluate();
        BDT = fTMVAEvaluator->getTMVA_EvaluationResult();
        fTMVA->Fill();
    }
    
    cout << "Writing event list (" << fTMVA->GetEntries() << " entries) to " << fOutputData->GetName() << endl;
    fOutputData->cd();
    fTMVA->Write();
    
    //////////////////////////////////////////////
    // get histogram with thrown events
    TFile* fTE = new TFile( fInputMSWFile.c_str() );
    if( fTE->IsZombie() )
    {
        cout << "Error opning mscw file " << fInputMSWFile << endl;
        cout << "Exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    VEffectiveAreaCalculatorMCHistograms* iMChistos = ( VEffectiveAreaCalculatorMCHistograms* )fTE->Get( "MChistos" );
    if( iMChistos && iMChistos->getHistogram_Emc( 8, 0 ) )
    {
        TH1D* h = ( TH1D* )iMChistos->getHistogram_Emc( 8, 0 )->Clone();
        if( h )
        {
            fOutputData->cd();
            cout << "Writing histogram with thrown events" << endl;
            h->Write();
        }
    }
    else
    {
        cout << "Error: cannot find MC histograms" << endl;
    }
    
    fOutputData->Close();
    
    cout << "(end of program)" << endl;
}


