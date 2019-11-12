/* writeParticleRateFilesForTMVA

   write files with particle number spectra for on (gamma) and off (protons+electrons) counts

   files are needed e.g. for setting the optimal cut value for TMVA cuts

   Input:
       * Combined anasum files for each zenith angle bin

   Output:
       * root file with signal and background rates per zenith angle and energy bin

*/

#include "VTMVARunData.h"
#include "VInstrumentResponseFunctionReader.h"
#include "VWPPhysSensitivityFile.h"

#include <iostream>
#include <string>
#include <vector>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"

using namespace std;

/*
 *  calculate signal rates from effective areas
 *
 */
void readRateGraphsFromEffectiveAreaFile( string iEffAreaDirectory,
        string iEffAreaFile,
        int i_noise,
        TGraph2DErrors* tRatePerEnergySignal )
{
    if( !tRatePerEnergySignal )
    {
        cout << "Error reading rate histograms from effective area files" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    cout << endl << "reading rate graphs from effective area files " << endl << endl;
    ///////////
    // zenith bins from simulations
    // (should match simulations)
    vector< int > i_ze;
    i_ze.push_back( 0 );
    i_ze.push_back( 20 );
    i_ze.push_back( 30 );
    i_ze.push_back( 35 );
    i_ze.push_back( 40 );
    i_ze.push_back( 45 );
    i_ze.push_back( 50 );
    i_ze.push_back( 55 );
    i_ze.push_back( 60 );
    i_ze.push_back( 65 ); 
    // az bin 8: averaged over all azimuth angles
    int i_Az = 8;
    // spectral index (does not matter, weighted histograms are
    // the same for all spectral indexes)
    double i_index = 2.4;
    // wobble offset 0.5
    double i_woff = 0.5;
    
    int m = 0;
    for( unsigned int i = 0; i < i_ze.size(); i++ )
    {
        string iFile = iEffAreaDirectory + "/" + iEffAreaFile;;
        stringstream i_zeString;
        if( i_ze[i] < 10 )
        {
            i_zeString << "Ze0" << i_ze[i] << "deg";
        }
        else
        {
            i_zeString << "Ze" << i_ze[i] << "deg";
        }
        if( iFile.find( "ZENITANGLE" ) != string::npos )
        {
            iFile.replace( iFile.find( "ZENITANGLE" ), std::string( "ZENITANGLE" ).size(), i_zeString.str() );
        }
        stringstream i_NoiseString;
        i_NoiseString << i_noise;
        if( iFile.find( "NOISE" ) != string::npos )
        {
            iFile.replace( iFile.find( "NOISE" ), std::string( "NOISE" ).size(), i_NoiseString.str() );
        }
        cout << "Openening effective area file: " << iFile << endl;
        VInstrumentResponseFunctionReader iTempIRFReader;
        if( iTempIRFReader.fillData( iFile, i_ze[i], i_woff, i_Az, i_index, i_noise, "A_REC" ) )
        {
            cout << "found data histogram with weights " << endl;
            // get weighted rate histogram and fill the corresponding graph
            if( iTempIRFReader.hWeightedRate005 )
            {
                for( int b = 1; b <= iTempIRFReader.hWeightedRate005->GetNbinsX(); b++ )
                {
                    // dE for differential calculation
                    double iE_min = TMath::Power( 10., iTempIRFReader.hWeightedRate005->GetBinLowEdge( b ) );
                    double iE_max = TMath::Power( 10., iTempIRFReader.hWeightedRate005->GetBinLowEdge( b )
                                                    + iTempIRFReader.hWeightedRate005->GetBinWidth( b ) );
                    double dE = iE_max - iE_min;
                    if( iTempIRFReader.hWeightedRate005->GetBinContent( b ) > 0 && dE > 0. )
                    {
                        // (note: unit of hWeightedRate005 is per min)
                        tRatePerEnergySignal->SetPoint( m, iTempIRFReader.hWeightedRate005->GetXaxis()->GetBinCenter( b ),
                                                        1. / cos( i_ze[i]*TMath::DegToRad() ),
                                                        iTempIRFReader.hWeightedRate005->GetBinContent( b ) / dE / 60. );
                        tRatePerEnergySignal->SetPointError( m, 0., 0., 0. );
                        m++;
                    }
                }
            }
        }
    }
}

/*
 *  read on/off/signal rate from an anasum file
 *
 *  note: on and signal rates are filled only if iReadMC is false
 *
 *  for on/signal rates, a gamma-ray source is required
 *
 */

void readRateGraphsFromAnasumFile( string iDataDir,
                                   TGraph2DErrors* tRatePerEnergyON,
                                   TGraph2DErrors* tRatePerEnergyOFF,
                                   TGraph2DErrors* tRatePerEnergySignal, bool iReadMC )
{
    if( !tRatePerEnergyON || !tRatePerEnergyOFF || !tRatePerEnergySignal )
    {
        cout << "Error getting rate graphs" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    Int_t runOn = 0;
    Double_t elevationOff = 0.;
    Double_t tOn = 0.;
    Double_t alphaNorm = 0.;
    
    /* Name and directory to the COMBINED anasum root file */
    TFile* f1 = new TFile( iDataDir.c_str(), "OPEN" );
    if( f1->IsZombie() )
    {
        cout << "Error opening anasum file: " << endl;
        cout << "\t" << iDataDir << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    cout << "Reading rate graphs from " << f1->GetName() << endl;
    
    char histnameOn[256], histnameOff[256];
    TTree* RunSum = ( TTree* )f1->Get( "total_1/stereo/tRunSummary" );
    if( !RunSum )
    {
        cout << "Error reading run summary tree from anasum file: " << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    RunSum->SetBranchAddress( "runOn", &runOn );
    RunSum->SetBranchAddress( "runOn", &runOn );
    RunSum->SetBranchAddress( "OffNorm", &alphaNorm );
    RunSum->SetBranchAddress( "elevationOff", &elevationOff );
    RunSum->SetBranchAddress( "tOn", &tOn );
    double iZenithperrun = 0.;
    double iSecantperrrun = 0.;
    
    double new_on = 0.;
    double new_off = 0.;
    double new_on_error = 0.;
    double new_off_error = 0.;
    double new_excess = 0.;
    
    unsigned int iBin = 0;
    double iBinWidth = 1.;
    
    // -----------------------------------------------------------------
    //
    // loop over the number of runs in the combined anasum root file
    //
    // check if the zenith angle is within a range of 5 degrees
    // (zenith angle intervals are used for an easier optimization)
    //
    // -----------------------------------------------------------------
    
    
    // counter
    int m = 0;
    ///////////////////////////
    // loop over all runs
    for( unsigned int i = 0; i < ( RunSum->GetEntries() - 1 ); i++ )
    {
    
        RunSum->GetEntry( i );
        // zenith angle and secant of zenith angle
        iZenithperrun = 90.0 - elevationOff;
        if( iZenithperrun > 89. )
        {
            iZenithperrun = 89.;
        }
        iSecantperrrun = 1. / cos( iZenithperrun * TMath::DegToRad() );
        
        // read counting histograms from disk
        sprintf( histnameOn, "run_%d/stereo/energyHistograms/herecCounts_on", runOn );
        TH1D* herecCounts_on = ( TH1D* )f1->Get( histnameOn );
        if( !herecCounts_on )
        {
            cout << "Error: on counting histogram not found" << endl;
        }
        sprintf( histnameOff, "run_%d/stereo/energyHistograms/herecCounts_off", runOn );
        TH1D* herecCounts_off = ( TH1D* )f1->Get( histnameOff );
        if( !herecCounts_off )
        {
            cout << "Error: off counting histogram not found" << endl;
        }
        
        cout << "Run Number: " << runOn;
        cout << ", Elevation angle = " << 90 - iZenithperrun;
        cout << ", Zenith angle = " << iZenithperrun;
        cout << ", Secant = " << iSecantperrrun;
        cout << ", Run Length = " << tOn << "s" << endl;
        
        if( tOn <= 0. )
        {
            continue;
        }
        ////////////////////////////////
        // get the bin content of the ON and OFF events per second and scale the off-counts with the alpha factor
        // assume same binning in on and off histograms
        iBin = herecCounts_on->GetNbinsX();
        for( unsigned int k = 1; k <= iBin; k++ )
        {
            // note! Rate graphs are dN/dE
            // unit is here dE (not d (log10 E))
            iBinWidth  = TMath::Power( 10., herecCounts_on->GetBinLowEdge( k )
                                       + herecCounts_on->GetBinWidth( k ) );
            iBinWidth -= TMath::Power( 10., herecCounts_on->GetBinLowEdge( k ) );
            if( iBinWidth <= 0. )
            {
                continue;
            }
            
            new_on        = ( herecCounts_on->GetBinContent( k ) ) / tOn / iBinWidth;
            new_on_error  = ( herecCounts_on->GetBinError( k ) ) / tOn / iBinWidth;
            
            new_off       = ( herecCounts_off->GetBinContent( k ) ) * alphaNorm / tOn / iBinWidth;
            new_off_error = ( herecCounts_off->GetBinError( k ) ) * alphaNorm / tOn / iBinWidth;
            
            if( !iReadMC )
            {
                /* ON rate (signal and background) */
                tRatePerEnergyON->SetPoint( m, herecCounts_on->GetBinCenter( k ), iSecantperrrun, new_on );
                tRatePerEnergyON->SetPointError( m, 0., 0., new_on_error );
                
                /* Signal rate (ON minus OFF) */
                new_excess = new_on - new_off;
                if( new_excess < 0. )
                {
                    new_excess = 0.;
                }
                tRatePerEnergySignal->SetPoint( m, herecCounts_on->GetBinCenter( k ), iSecantperrrun, new_excess );
                tRatePerEnergySignal->SetPointError( m, 0., 0., sqrt( new_on_error * new_on_error
                                                     + new_off_error * new_off_error ) );
            }
            
            /* OFF rate */
            tRatePerEnergyOFF->SetPoint( m, herecCounts_off->GetBinCenter( k ), iSecantperrrun, new_off );
            tRatePerEnergyOFF->SetPointError( m, 0., 0., new_off_error );
            
            m++;
        }
    }
}


/*

    fill graphs with particle numbers for on (gamma) and off (protons+electrons) counts

    files are needed e.g. for setting the optimal cut value for TMVA cuts

*/
int main( int argc, char* argv[] )
{
    cout << endl;
    cout << "writeParticleRateFilesForTMVA (version ";
    cout << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "=================================================" << endl;
    cout << endl ;
    if( argc < 3 )
    {
        cout << "writeParticleRateFilesForTMVA <combined anasum root file> <outputfile> ";
        cout << "[directory with effective areas] [effective area file name] [Noise level]" << endl;
        cout << endl;
        cout << "       write signal and backround rates vs energy and zenith angle" << endl;
        cout << endl;
        cout << "       use data for signal and background: only give <combined anasum root file> " << endl;
        cout << "       use data for background and MC for signal: give additionally the effective area file name" << endl;
        cout << endl;
        cout << "       [effective area file name]  zenith angle and noise value are added to the name by" << endl;
        cout << "                                   search and replacement: e.g. effArea-ZENITANGLE-NOISE-SoftCut.root" << endl;
        cout << endl << endl;
        exit( EXIT_SUCCESS );
    }
    
    
    string iDataDir = argv[1];
    string iOutFil = argv[2];
    string iEffAreaDir = "";
    string iEffAreaFile = "";
    int    iNoiseLevel = 0;
    if( argc == 6 )
    {
        iEffAreaDir = argv[3];
        iEffAreaFile = argv[4];
        iNoiseLevel = atoi( argv[5] );
    }
    
    /////////////////////////////////////////////////////////////////
    // 2D graphs which contain all results
    TGraph2DErrors* tRatePerEnergyON     = new TGraph2DErrors();
    tRatePerEnergyON->SetTitle( "differential signal+background rate (1/s/TeV); Log_{10}E [TeV]; Secant(zenith angle)" );
    tRatePerEnergyON->SetMarkerStyle( 20 );
    TGraph2DErrors* tRatePerEnergyOFF    = new TGraph2DErrors();
    tRatePerEnergyOFF->SetTitle( "differential background rate (1/s/TeV); Log_{10}E [TeV]; Secant(zenith angle)" );
    tRatePerEnergyOFF->SetMarkerStyle( 20 );
    TGraph2DErrors* tRatePerEnergySignal = new TGraph2DErrors();
    tRatePerEnergySignal->SetTitle( "differential signal rate (1/s/TeV); Log_{10}E [TeV]; Secant(zenith angle)" );
    tRatePerEnergySignal->SetMarkerStyle( 20 );
    
    /////////////////////////////////////////////////////////////////
    // rate graphs from anasum file
    readRateGraphsFromAnasumFile( iDataDir, tRatePerEnergyON, 
                                  tRatePerEnergyOFF, 
                                  tRatePerEnergySignal, 
                                  ( argc > 3 ) );
    
    /////////////////////////////////////////////////////////////////
    // signal rate graphs from effective area files
    if( argc == 6 )
    {
        readRateGraphsFromEffectiveAreaFile( iEffAreaDir, iEffAreaFile, iNoiseLevel, tRatePerEnergySignal );
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // write particle number file
    
    cout << endl;
    cout << "writing rate graphs to " << iOutFil << endl;
    
    TFile* fRateFile = new TFile( iOutFil.c_str(), "RECREATE" );
    if( fRateFile->IsZombie() )
    {
        cout << "Error writing rate graph file: " << iOutFil << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    // write graphs of ON rate to file
    tRatePerEnergyON->Write( "gONRate" );
    
    // write graphs of OFF rate to file
    tRatePerEnergyOFF->Write( "gBGRate" );
    
    // write graphs of signal rate to file
    tRatePerEnergySignal->Write( "gSignalRate" );
    
    fRateFile->Close();
    
    cout << "\t finished (success)" << endl;
    
    return 0;
}
