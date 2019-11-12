/*! \file VTS.calculateCrabRateFromMC.cpp
    \brief calculate gamma ray rate from Crab Nebula with effective areas from MC


*/

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"

#include "CEffArea.h"
#include "VGlobalRunParameter.h"
#include "VMonteCarloRateCalculator.h"
#include "VEnergySpectrumfromLiterature.h"

#include <iostream>
#include <vector>


using namespace std;

int main( int argc, char* argv[] )
{
    const double index = 2.4;
    bool bDebug = false;
    
    cout << endl;
    cout << "VTS.calculateCrabRateFromMC (" << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "--------------------------------" << endl;
    if( argc != 4 )
    {
        cout << "VTS.calculateCrabRateFromMC <effective area file> <outputfile> <energy threshold [TeV]>" << endl;
        cout << endl;
        exit( 0 );
    }
    cout << endl;
    cout << "selecting effective areas with spectral index " << index << endl;
    
    string ieff = argv[1];
    string ioffile = argv[2];
    // energy threshold
    double fEnergyThreshold = atof( argv[3] );
    cout << "energy threshold: " << fEnergyThreshold << " TeV" << endl;
    if( fEnergyThreshold < 1.e-2 )
    {
        fEnergyThreshold = 1.e-2;    // take care of log(fEnergyThreshold)
    }
    
    double ze = 0.;
    int az = 0;
    double Woff = 0.;
    int noise = 0;
    double pedvar = 0.;
    unsigned int nrates = 0;
    double MCrate[1000];
    
    TFile* f1 = new TFile( ioffile.c_str(), "RECREATE" );
    if( f1->IsZombie() )
    {
        cout << "error opening output file: " << ioffile << endl;
        exit( 0 );
    }
    TTree* fMC = new TTree( "fMCRate", "MC rate predictions" );
    fMC->Branch( "ze", &ze, "ze/D" );
    fMC->Branch( "az", &az, "az/I" );
    fMC->Branch( "Woff", &Woff, "Woff/D" );
    fMC->Branch( "noise", &noise, "noise/I" );
    fMC->Branch( "pedvar", &pedvar, "pedvar/D" );
    fMC->Branch( "nrates", &nrates, "nrates/i" );
    fMC->Branch( "MCrate", MCrate, "MCrate[nrates]/D" );
    
    TFile* f = new TFile( ieff.c_str() );
    if( f->IsZombie() )
    {
        cout << "error opening file with effective areas: " << ieff << endl;
        exit( 0 );
    }
    TTree* t = ( TTree* )gDirectory->Get( "fEffArea" );
    if( !t )
    {
        exit( 0 );
    }
    CEffArea* c = new CEffArea( t );
    
    // Crab Nebula Spectra
    vector< unsigned int > fID;
    fID.push_back( 1 );                           // Whipple
    fID.push_back( 5 );                           // HEGRA
    fID.push_back( 2 );                           // H.E.S.S.
    fID.push_back( 6 );                           // MAGIC
    vector< VEnergySpectrumfromLiterature* > fESpecFun;
    char hname[2000];
    VGlobalRunParameter* fRunParameter = new VGlobalRunParameter();
    sprintf( hname, "%s/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat", fRunParameter->getDirectory_EVNDISPAnaData().c_str() );
    VEnergySpectrumfromLiterature* fESpec = new VEnergySpectrumfromLiterature( hname );
    if( fESpec->isZombie() )
    {
        exit( -1 );
    }
    for( unsigned int i = 0; i < fID.size(); i++ )
    {
        fESpecFun.push_back( fESpec );
    }
    
    // error checks
    if( fESpecFun.size() != fID.size() )
    {
        cout << "Error: vector mismatch - check code" << endl;
        exit( -1 );
    }
    
    // print everything
    for( unsigned int i = 0; i < fESpecFun.size(); i++ )
    {
        if( fESpecFun[i] )
        {
            fESpecFun[i]->listValues( fID[i] );
        }
    }
    
    // rate calculator
    VMonteCarloRateCalculator* fMCR = new VMonteCarloRateCalculator();
    
    cout << endl;
    cout << "reading " << c->fChain->GetEntries() << " effective areas " << endl;
    
    // loop over all effective areas and calculate expected rates
    unsigned int iN = c->fChain->GetEntries();
    for( unsigned int i = 0; i < iN; i++ )
    {
        c->GetEntry( i );
        
        if( TMath::Abs( c->index - index ) > 1.e-2 && c->fChain->GetEntries() != 1 )
        {
            continue;
        }
        if( i % 5000  == 0 )
        {
            cout << "now at entry " << i << endl;
        }
        
        ze     = c->ze;
        az     = c->az;
        Woff   = c->Woff;
        noise  = c->noise;
        pedvar = c->pedvar;
        
        vector< double > fenergy;
        vector< double > feffectivearea;
        for( int e = 0; e < c->nbins; e++ )
        {
            fenergy.push_back( c->e0[e] );
            feffectivearea.push_back( c->eff[e] );
        }
        
        nrates = fESpecFun.size();
        // hardwire Whipple spectrum (much faster than outer functino call)
        // nrates = 1;
        // MCrate[0] = fMCR->getMonteCarloRate( fenergy, feffectivearea, -2.440, 3.250e-11, 1., fEnergyThreshold, 1.e7, bDebug );
        for( unsigned int t = 0; t < fESpecFun.size(); t++ )
        {
            cout << "XXXX " << fESpecFun[t]->getPowerLaw_Index( fID[t] ) << "\t" << fESpecFun[t]->getPowerLaw_FluxConstant_at1TeV( fID[t] ) * 1.e11 << endl;
            MCrate[t] = fMCR->getMonteCarloRate( fenergy, feffectivearea,
                                                 fESpecFun[t]->getPowerLaw_Index( fID[t] ),
                                                 fESpecFun[t]->getPowerLaw_FluxConstant_at1TeV( fID[t] ),
                                                 1.,
                                                 fEnergyThreshold, 1.e7, bDebug );
            if( bDebug )
            {
                cout << "MC Rate " << ze << "\t" << MCrate[t] << endl;
            }
        }
        fMC->Fill();
    }
    cout << endl;
    cout << "writing results to " << f1->GetName() << endl;
    f1->cd();
    fMC->Write();
    f1->Close();
}
