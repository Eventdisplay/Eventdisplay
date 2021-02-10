/*   macro to plot average and maximum
 *   number of pe as function of energy
 *   and core distance
 *
 *   input are dst.root files produced with the converter
 *   using the -pe flag
 *
 */

#include "TCanvas.h"
#include "TChain.h"
#include "TH2D.h"
#include "TH1F.h"


#include <iostream>
#include <map>
#include <string>
#include <vector>

void plot( string iFileDSTFileName )
{
    TChain* fData = new TChain( "dst" );
    if( !fData )
    {
        return;
    }
    fData->Add( iFileDSTFileName.c_str() );
    const unsigned int VDST_MAXTELESCOPES = 100;
    const unsigned int VDST_MAXCHANNELS = 12000;
    unsigned int ntel_data = 0;
    unsigned int tel_data[VDST_MAXTELESCOPES];
    unsigned short int fPe[VDST_MAXTELESCOPES][VDST_MAXCHANNELS];
    float MCe0 = 0.;
    float MCxcore = 0.;
    float MCycore = 0.;
    
    fData->SetBranchAddress( "ntel_data", &ntel_data );
    fData->SetBranchAddress( "tel_data", &tel_data );
    fData->SetBranchAddress( "Pe", fPe );
    fData->SetBranchAddress( "MCe0", &MCe0 );
    fData->SetBranchAddress( "MCxcore", &MCxcore );
    fData->SetBranchAddress( "MCycore", &MCycore );
    
    // read number of channels per telescope and pixel positions
    // (note this assume that the dst files are for one telescope type only)
    TChain* ftelconfig = new TChain( "telconfig" );
    ftelconfig->Add( iFileDSTFileName.c_str() );
    if( !fData )
    {
        return;
    }
    map< int,  unsigned int > fNPixel;
    map< int,  float > fTelX;
    map< int, float > fTelY;
    int TelID = 0;
    unsigned int NPixel = 0;
    float TelX = 0.;
    float TelY = 0.;
    ftelconfig->SetBranchAddress( "TelID", &TelID );
    ftelconfig->SetBranchAddress( "TelX", &TelX );
    ftelconfig->SetBranchAddress( "TelY", &TelY );
    ftelconfig->SetBranchAddress( "NPixel", &NPixel );
    for( unsigned int t = 0; t < ftelconfig->GetEntries(); t++ )
    {
        ftelconfig->GetEntry( t );
        fNPixel[TelID] = NPixel;
        fTelX[TelID] =  TelX;
        fTelY[TelID] = TelY;
    }
    cout << "Number of telescopes: " << fNPixel.size() << endl;
    cout << "Number of pixels per telesope: " << NPixel << endl;
    
    ///////////////////////////////////////////////
    // histogram definitions
    const unsigned int nEbins = 25;
    const unsigned int nRbins = 20;
    TH2D* hMeanPe = new TH2D( "hMeanPE", "", nEbins, -2., 3, nRbins, 0., 2000. );
    hMeanPe->SetXTitle( "log_{10} E_{MC}/TeV" );
    hMeanPe->SetYTitle( "impact distance (m)" );
    hMeanPe->SetStats( 0 );
    TH2D* hMaxPe = new TH2D( "hMaxPE", "", nEbins, -2., 3, nRbins, 0., 2000. );
    hMaxPe->SetXTitle( "log_{10} E_{MC}/TeV" );
    hMaxPe->SetYTitle( "impact distance (m)" );
    hMaxPe->SetStats( 0 );
    
    // 1D histograms
    vector< vector< TH1F* > > fVPe;
    char hname[200];
    cout << "Initializing " << nEbins* nRbins << " histograms" << endl;
    for( unsigned int e = 0; e < nEbins; e++ )
    {
        vector< TH1F*> tempfPE;
        for( unsigned r = 0; r < nRbins; r++ )
        {
            sprintf( hname, "hPE_%d_%d", e, r );
            tempfPE.push_back( new TH1F( hname, "", 100, 0., 2000. ) );
        }
        fVPe.push_back( tempfPE );
    }
    
    // loop over all entries in histogram
    cout << "filling histograms from " << fData->GetEntries() << " entries" << endl;
    int nEbin = 0;
    int nRbin = 0;
    float R = 0.;
    for( int i = 0; i < fData->GetEntries(); i++ )
    {
        fData->GetEntry( i );
        
        nEbin = hMeanPe->GetXaxis()->FindBin( log10( MCe0 ) ) - 1;
        
        for( unsigned int n = 0; n < ntel_data; n++ )
        {
            if( fTelX.find( tel_data[n] ) != fTelX.end() )
            {
                R = sqrt( ( MCxcore - fTelX[tel_data[n]] ) * ( MCxcore - fTelX[tel_data[n]] ) + ( MCycore - fTelY[tel_data[n]] ) * ( MCycore - fTelY[tel_data[n]] ) );
                
                nRbin = hMeanPe->GetYaxis()->FindBin( R ) - 1;
                
                if( nEbin < ( int )fVPe.size() && nRbin < ( int )fVPe[nEbin].size() )
                {
                    for( unsigned int p = 0; p < fNPixel[tel_data[n]]; p++ )
                    {
                        if( fPe[n][p] > 0 )
                        {
                            fVPe[nEbin][nRbin]->Fill( fPe[n][p] );
                        }
                    }
                }
            }
        }
        
    }
    
    // fill mean and maximum
    for( unsigned int e = 0; e < fVPe.size(); e++ )
    {
        for( unsigned int r = 0; r < fVPe[e].size(); r++ )
        {
            hMeanPe->SetBinContent( e + 1, r + 1, fVPe[e][r]->GetMean() );
            for( int b = fVPe[e][r]->GetNbinsX(); b > 0; b-- )
            {
                if( fVPe[e][r]->GetBinContent( b ) > 0 )
                {
                    hMaxPe->SetBinContent( e + 1, r + 1, fVPe[e][r]->GetXaxis()->GetBinCenter( b ) );
                    break;
                }
            }
        }
    }
    
    
    // plotting
    //
    TCanvas* cCMean = new TCanvas( "cCMean", "mean pe", 10, 10, 800, 600 );
    cCMean->SetGridx( 0 );
    cCMean->SetGridy( 0 );
    hMeanPe->Draw( "colz" );
    
    TCanvas* cCMax = new TCanvas( "cCMax", "Max pe", 910, 10, 800, 600 );
    cCMax->SetGridx( 0 );
    cCMax->SetGridy( 0 );
    hMaxPe->Draw( "colz" );
}


