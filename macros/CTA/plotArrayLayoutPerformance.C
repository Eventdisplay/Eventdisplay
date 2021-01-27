/*   \file plotArrayLayoutPerformance.C

     plot performance (sensitivity, angular resolution)
     in energy bins vs array layout and scaling

     input is a tree produced with writeCTAWPPhysSensitivityTree

     - optimized for prod3
     - some hard wired parameters


     todotodo
     // add all north layouts
     // different stagins, all in one?
*/


#include "cIRFData.C"
#include "../../inc/VCTASensitivityRequirements.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TTree.h"

using namespace std;

class CPerformance
{
    private:
    
        // base performance for PPUT calculation
        unsigned int fBasePerformance;
        TH1F*        fBasePerformanceHistogram;
        TGraph*      fReqAngularResolution;
        TGraph*      fReqEnergyResolution;
        
        void    fillPPUTRanking();
        double  getBasePerformance( double iE );
        double  getBasePerformanceAngularResolution( double iE );
        double  getBasePerformanceEnergyResolution( double iE );
        
    public:
    
        string fName;
        vector< string > fArrayLayout;
        bool   fSouthSite;
        bool   fKSP;
        
        vector< vector< string > > fArrayLayoutLabels;
        vector< TGraphErrors* > fSensitivity;
        vector< vector< vector< double > > > fSensitivityPPUT;
        vector< vector< vector< double > > > fSensitivityPPUT_E;
        vector< vector< vector< double > > > fSensitivityPPUT_N;
        vector< vector< vector< int > > > fSensitivityPPUTRanking;
        
        vector< TGraphErrors* > fAngRes;
        vector< vector< vector< double > > > fAngularResolution;
        vector< TGraphErrors* > fERes;
        vector< vector< vector< double > > > fEnergyResolution;
        bool         fPlotResolutionRequirement;
        
        vector< unsigned int >  fCounter;
        unsigned int fMSTCounter;
        unsigned int fSSTCounter;
        int fObservingTime_s;
        float fDistanceCameraCentre_deg;
        string fPointingDirection;
        
        // plotting
        vector< unsigned int > fEnergyBins_min;
        vector< unsigned int > fEnergyBins_max;
        vector< TText* > fEnergyLabels;
        int fColor;
        
        CPerformance( string iName, bool iSouth, vector< string > iA, unsigned int iLowEnergyBin,
                      unsigned int iNScalings, int iColor = 1, int iMarkerStyle = 24, bool iKSP = false );
        ~CPerformance();
        
        void addSensitivityPoint( unsigned int iArrayID, unsigned int iScalings,
                                  unsigned int iEbin, double iE_log10TeV,
                                  double iDiffSens, double iDiffSensError = 0.,
                                  double iAngres = 0., double iEres = 0.,
                                  int iTreeEntry = 0 );
                                  
        unsigned int getArrayID( string iAL, int iObsTime, double iOffset_deg );
        
        void setCameraOffset( float f )
        {
            fDistanceCameraCentre_deg = f;
        }
        void setEnergyBins( unsigned int iLowEnergyOffset = 0 );
        void setEnergyBins( vector< unsigned int > iEnergyBins_min, vector< unsigned int > iEnergyBins_max, vector< TText* > iEnergyLabels )
        {
            fEnergyBins_min = iEnergyBins_min;
            fEnergyBins_max = iEnergyBins_max;
            fEnergyLabels = iEnergyLabels;
        }
        void setObservingTime( int iObservingTime )
        {
            fObservingTime_s = iObservingTime;
        }
        bool setPerformanceBase( unsigned int iPerformanceID, string iWPPhysFileName = "" );
        void setPointingDirection( string iPointingDirection )
        {
            fPointingDirection = iPointingDirection;
        }
        void setPlotResolutionRequirement()
        {
            fPlotResolutionRequirement = true;
        }
        void setKSP( bool iKSP = false )
        {
            fKSP = iKSP;
        }
        void setSouth( bool iSouth = true )
        {
            fSouthSite = iSouth;
        }
        void terminate();
};

CPerformance::CPerformance( string iName, bool iSouth, vector< string > iA, unsigned int iLowEnergyBin,
                            unsigned int iNScalings, int iColor, int iMarkerStyle, bool iKSP )
{
    fName = iName;
    fMSTCounter = 1;
    fSSTCounter = 3;
    fObservingTime_s = 50 * 3600; // in [s]; default = 50 h
    fDistanceCameraCentre_deg = 2.5;
    fArrayLayout = iA;
    fColor = iColor;
    fPlotResolutionRequirement = false;
    fPointingDirection = "";
    setSouth( iSouth );
    setKSP( iKSP );
    
    // note: setSouth and setKSP must be called before setting the energy bins;
    setEnergyBins( iLowEnergyBin );
    
    fBasePerformance = 0;
    fBasePerformanceHistogram = 0;
    
    // energy bins
    for( unsigned int i = 0; i < fEnergyBins_min.size(); i++ )
    {
        vector< double > i1( fArrayLayout.size(), 1. );
        vector< double > i0( fArrayLayout.size(), 0. );
        vector< int > ii( fArrayLayout.size(), 0 );
        vector< vector< double > > ii1;
        vector< vector< double > > ii0;
        vector< vector< int > > iii;
        
        // second dimension: array scalings
        for( unsigned int s = 0; s < iNScalings; s++ )
        {
            // third dimension: array layouts
            ii1.push_back( i1 );
            ii0.push_back( i0 );
            iii.push_back( ii );
            
        }
        fSensitivity.push_back( new TGraphErrors( 1 ) );
        fSensitivity.back()->SetMarkerStyle( iMarkerStyle );
        if( iMarkerStyle != 21 )
        {
            fSensitivity.back()->SetMarkerSize( 0.8 );
        }
        fSensitivity.back()->SetTitle( "" );
        fSensitivity.back()->SetMarkerColor( iColor );
        fSensitivity.back()->SetLineColor( iColor );
        fSensitivityPPUT.push_back( ii1 );
        fSensitivityPPUT_E.push_back( ii0 );
        fSensitivityPPUT_N.push_back( ii0 );
        fSensitivityPPUTRanking.push_back( iii );
        fAngularResolution.push_back( ii1 );
        fEnergyResolution.push_back( ii1 );
        
        fAngRes.push_back( new TGraphErrors( 1 ) );
        fAngRes.back()->SetMarkerStyle( 5 );
        fAngRes.back()->SetMarkerSize( 0.8 );
        fAngRes.back()->SetTitle( "" );
        fAngRes.back()->SetMarkerColor( iColor );
        fAngRes.back()->SetLineColor( iColor );
        
        fERes.push_back( new TGraphErrors( 1 ) );
        fERes.back()->SetMarkerStyle( 2 );
        fERes.back()->SetMarkerSize( 0.8 );
        fERes.back()->SetTitle( "" );
        fERes.back()->SetMarkerColor( iColor );
        fERes.back()->SetLineColor( iColor );
        
        fCounter.push_back( 0 );
    }
    
    cout << "Initialize performance class for " << fName << ": " << endl;
    cout << "\t # of energy bins: " << fEnergyBins_min.size() << "(" << fSensitivityPPUT.size() << ")" << endl;
    cout << "\t # of array layouts: " << fArrayLayout.size() << "(" << fSensitivityPPUT[0][0].size() << ")" << endl;
    cout << "\t # of scalings: " << iNScalings << "(" << fSensitivityPPUT[0].size() << ")" << endl;
    
    // angular resolution requirement
    //
    // from MAN-PO/121004, version 2.5, June 14, 2013
    // from Figure 4
    fReqAngularResolution = new TGraph( 18 );
    fReqAngularResolution->SetPoint( 0, log10( 0.024 ), 0.394 );
    fReqAngularResolution->SetPoint( 1, log10( 0.028 ), 0.317 );
    fReqAngularResolution->SetPoint( 2, log10( 0.041 ), 0.220 );
    fReqAngularResolution->SetPoint( 3, log10( 0.059 ), 0.171 );
    fReqAngularResolution->SetPoint( 4, log10( 0.102 ), 0.127 );
    fReqAngularResolution->SetPoint( 5, log10( 0.167 ), 0.104 );
    fReqAngularResolution->SetPoint( 6, log10( 0.370 ), 0.081 );
    fReqAngularResolution->SetPoint( 7, log10( 1.147 ), 0.060 );
    fReqAngularResolution->SetPoint( 8, log10( 1.674 ), 0.054 );
    fReqAngularResolution->SetPoint( 9, log10( 2.905 ), 0.048 );
    fReqAngularResolution->SetPoint( 10, log10( 6.192 ), 0.041 );
    fReqAngularResolution->SetPoint( 11, log10( 9.638 ), 0.038 );
    fReqAngularResolution->SetPoint( 12, log10( 16.568 ), 0.034 );
    fReqAngularResolution->SetPoint( 13, log10( 33.101 ), 0.030 );
    fReqAngularResolution->SetPoint( 14, log10( 51.349 ), 0.028 );
    fReqAngularResolution->SetPoint( 15, log10( 88.158 ), 0.026 );
    fReqAngularResolution->SetPoint( 16, log10( 122.198 ), 0.025 );
    fReqAngularResolution->SetPoint( 17, log10( 256.706 ), 0.023 );
    
    // energy resolution requirement
    // from Figure 5
    fReqEnergyResolution = new TGraph( 4 );
    fReqEnergyResolution->SetPoint( 0, log10( 0.03 ), 0.40 );
    fReqEnergyResolution->SetPoint( 1, log10( 0.10 ), 0.20 );
    fReqEnergyResolution->SetPoint( 2, log10( 1.00 ), 0.10 );
    fReqEnergyResolution->SetPoint( 3, log10( 10.00 ), 0.10 );
}

CPerformance::~CPerformance()
{
    for( unsigned int i = 0; i < fSensitivity.size(); i++ )
    {
        if( fSensitivity[i] )
        {
            delete fSensitivity[i];
        }
    }
}

/*
 * set energy bins for PPUT calculation
 *
 * different for South and North
 *
*/
void CPerformance::setEnergyBins( unsigned int iLowEnergyBin )
{
    fEnergyBins_min.clear();
    fEnergyBins_max.clear();
    fEnergyLabels.clear();
    
    // TText label positions
    double i_Ty = 0.85;
    double i_Tx = 0.15;
    
    // lowest energy bin (usually 30 GeV?)
    unsigned int iLowestEnergyBin = 2;
    
    /////////////////////////////////////////////////////
    // South
    if( fSouthSite )
    {
        if( !fKSP )
        {
            fEnergyBins_min.push_back( iLowestEnergyBin + iLowEnergyBin );    // [31.6,50.1] GeV
            fEnergyBins_max.push_back( 20 );   // [80,125] TeV
            if( iLowEnergyBin == 1 )
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT E = [50 GeV, 125 TeV]" ) );
            }
            else
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT E = [20 GeV, 125 TeV]" ) );
            }
            
            fEnergyBins_min.push_back( 6 );    // [125] GeV
            fEnergyBins_max.push_back( 18 );   // [80,] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT MST E = [125 GeV, 50 TeV]" ) );
            
            fEnergyBins_min.push_back( iLowestEnergyBin + iLowEnergyBin );  // [31.6,39.8] GeV
            fEnergyBins_max.push_back( 5 );   //  [80,125] GeV
            if( iLowEnergyBin == 1 )
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Low E = [50 GeV, 125 GeV]" ) );
            }
            else
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Low E = [20 GeV, 125 GeV]" ) );
            }
            
            fEnergyBins_min.push_back( 5 );   //  [80,125] GeV
            fEnergyBins_max.push_back( 10 );  //  [0.8, 1.25] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Mid E = [80 GeV, 1.25 TeV]" ) );
            
            fEnergyBins_min.push_back( 11 );  //  [1.25, 2.] TeV
            fEnergyBins_max.push_back( 15 );  //  [8, 12.5] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT High E = [1.25 TeV, 12.5 TeV]" ) );
            
            fEnergyBins_min.push_back( 15 );  //  [8, 12.5] TeV
            fEnergyBins_max.push_back( 20 );   // [80,125] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Very-high E = [8 TeV, 125 TeV]" ) );
        }
        // KSP studies
        else
        {
            // full energy range
            fEnergyBins_min.push_back( iLowestEnergyBin + iLowEnergyBin );    // [31.6,50.1] GeV
            fEnergyBins_max.push_back( 20 );   // [80,125] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT E = [20 GeV, 125 TeV]" ) );
            
            // KSP 1
            
            // KSP 2
            fEnergyBins_min.push_back( 5 );
            fEnergyBins_max.push_back( 20 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 2 E = [80 GeV, 125 TeV]" ) );
            
            // KSP 3
            fEnergyBins_min.push_back( 7 );
            fEnergyBins_max.push_back( 18 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 3 E = [200 GeV, 50 TeV]" ) );
            
            // KSP 4
            fEnergyBins_min.push_back( 4 );
            fEnergyBins_max.push_back( 15 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 4 E = [50 GeV, 12.5 TeV]" ) );
            
            // KSP 5
            fEnergyBins_min.push_back( 3 );
            fEnergyBins_max.push_back( 20 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 5 E = [20 GeV, 125 TeV]" ) );
            
            // KSP 6
            fEnergyBins_min.push_back( 5 );
            fEnergyBins_max.push_back( 20 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 6 E = [80 GeV, 125 TeV]" ) );
            
            // KSP 7
            fEnergyBins_min.push_back( 5 );
            fEnergyBins_max.push_back( 16 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP1 7 E = [80 GeV, 20 TeV]" ) );
            
            // KSP 8
            fEnergyBins_min.push_back( 2 );
            fEnergyBins_max.push_back( 18 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 8 E = [20 GeV, 50 TeV]" ) );
        }
    }
    /////////////////////////////////////////////////////
    // North
    else
    {
        if( !fKSP )
        {
            fEnergyBins_min.push_back( iLowestEnergyBin + iLowEnergyBin );    // [31.6,50.1] GeV
            fEnergyBins_max.push_back( 18 );   // [80,125] TeV
            if( iLowEnergyBin == 1 )
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT E = [50 GeV, 50 TeV]" ) );
            }
            else
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT E = [20 GeV, 50 TeV]" ) );
            }
            
            fEnergyBins_min.push_back( 6 );    // [125] GeV
            fEnergyBins_max.push_back( 18 );   // [80,] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT MST E = [125 GeV, 50 TeV]" ) );
            
            fEnergyBins_min.push_back( iLowestEnergyBin + iLowEnergyBin );  // [31.6,39.8] GeV
            fEnergyBins_max.push_back( 5 );   //  [80,125] GeV
            if( iLowEnergyBin == 1 )
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Low E = [50 GeV, 125 GeV]" ) );
            }
            else
            {
                fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Low E = [20 GeV, 125 GeV]" ) );
            }
            
            fEnergyBins_min.push_back( 5 );   //  [80,125] GeV
            fEnergyBins_max.push_back( 10 );  //  [0.8, 1.25] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Mid E = [80 GeV, 1.25 TeV]" ) );
            
            fEnergyBins_min.push_back( 11 );  //  [1.25, 2.] TeV
            fEnergyBins_max.push_back( 15 );  //  [8, 12.5] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT High E = [1.25 TeV, 12.5 TeV]" ) );
            
            fEnergyBins_min.push_back( 15 );  //  [8, 12.5] TeV
            fEnergyBins_max.push_back( 18 );   // [80,125] TeV
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT Very-high E = [8 TeV, 50 TeV]" ) );
        }
        else
        {
            // KSP 1
            
            // KSP 2
            fEnergyBins_min.push_back( 5 );
            fEnergyBins_max.push_back( 16 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 2 E = [80 GeV, 20 TeV]" ) );
            
            // KSP 3
            /*            fEnergyBins_min.push_back( 7 );
                        fEnergyBins_max.push_back( 18 );   //
                        fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 3 E = [200 GeV, 50 TeV]" ) ); */
            
            // KSP 4
            fEnergyBins_min.push_back( 4 );
            fEnergyBins_max.push_back( 15 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 4 E = [50 GeV, 12.5 TeV]" ) );
            
            // KSP 5
            fEnergyBins_min.push_back( 3 );
            fEnergyBins_max.push_back( 20 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 5 E = [20 GeV, 125 TeV]" ) );
            
            /*            // KSP 6
                        fEnergyBins_min.push_back( 5 );
                        fEnergyBins_max.push_back( 20 );   //
                        fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 6 E = [80 GeV, 125 TeV]" ) ); */
            
            // KSP 7
            fEnergyBins_min.push_back( 5 );
            fEnergyBins_max.push_back( 16 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP1 7 E = [80 GeV, 20 TeV]" ) );
            
            // KSP 8
            fEnergyBins_min.push_back( 2 );
            fEnergyBins_max.push_back( 16 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 8 E = [20 GeV, 20 TeV]" ) );
            
            // KSP 9
            fEnergyBins_min.push_back( 2 );
            fEnergyBins_max.push_back( 15 );   //
            fEnergyLabels.push_back( new TText( i_Tx, i_Ty, "PPUT KSP 9 E = [20 GeV, 12.5 TeV]" ) );
            
        }
    }
}

/*
    search for first matching array layout and return ID
*/
unsigned int CPerformance::getArrayID( string iA, int iObsTime, double iOffset_deg )
{
    for( unsigned int i = 0; i < fArrayLayout.size(); i++ )
    {
        if( iA.find( fArrayLayout[i] ) != string::npos )
        {
            if( iObsTime == fObservingTime_s )
            {
                if( TMath::Abs( iOffset_deg - fDistanceCameraCentre_deg ) < 1.e-3 )
                {
                    return i;
                }
            }
        }
    }
    
    return 9999;
}


bool CPerformance::setPerformanceBase( unsigned int iPerformanceID, string iWPPhysFileName )
{
    fBasePerformance = iPerformanceID;
    if( fBasePerformance < 3 )
    {
        return true;
    }
    // read a DiffSens histogram as base performance histogram from file
    TFile* iF = new TFile( iWPPhysFileName.c_str() );
    if( iF->IsZombie() )
    {
        cout << "Error opening file with base performance histogram " << iF->GetName() << endl;
        return false;
    }
    if( fDistanceCameraCentre_deg < 1.e-2 )
    {
        fBasePerformanceHistogram = ( TH1F* )iF->Get( "DiffSens" );
    }
    else
    {
        TH2F* h = ( TH2F* )iF->Get( "DiffSens_offaxis" );
        if( h )
        {
            fBasePerformanceHistogram = new TH1F( "DiffSens", "", h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
            int iy = h->GetYaxis()->FindBin( fDistanceCameraCentre_deg );
            for( int i = 1; i <= h->GetXaxis()->GetNbins(); i++ )
            {
                fBasePerformanceHistogram->SetBinContent( i, h->GetBinContent( i, iy ) );
                fBasePerformanceHistogram->SetBinError( i, h->GetBinError( i, iy ) );
            }
        }
    }
    if( !fBasePerformanceHistogram )
    {
        cout << "Error: base performance histogram not found in " << iF->GetName() << endl;
        return false;
    }
    cout << "Performance base file " << iWPPhysFileName << endl;
    
    return true;
}

/*
   get energy resolution requirement


*/

double CPerformance::getBasePerformanceEnergyResolution( double iE )
{
    if( fReqEnergyResolution )
    {
        return fReqEnergyResolution->Eval( iE );
    }
    
    return 0.;
}


/*
   get angular resolution requirement


*/

double CPerformance::getBasePerformanceAngularResolution( double iE )
{
    if( fReqAngularResolution )
    {
        return fReqAngularResolution->Eval( iE );
    }
    
    return 0.;
}

/*

   get the baseline performance needed for the PPUT calculation

   fBasePerformance = 0: use required sensitivity (South)
   fBasePerformance = 1: use required sensitivity (North)

   baseline performance might be relative to a given histogram

*/

double CPerformance::getBasePerformance( double iE )
{
    if( fBasePerformance == 0 )
    {
        return VCTASensitivityRequirements::Flux_req50_E2erg_south( TMath::Power( 10., iE ) );
    }
    else if( fBasePerformance == 1 )
    {
        return VCTASensitivityRequirements::Flux_req50_E2erg_north( TMath::Power( 10., iE ) );
    }
    else
    {
        if( fBasePerformanceHistogram )
        {
            return fBasePerformanceHistogram->GetBinContent( fBasePerformanceHistogram->FindBin( iE ) );
        }
    }
    
    
    return 1.;
}

/*

   add a new sensitivity point to PPUT term

*/

void CPerformance::addSensitivityPoint( unsigned int iArrayID, unsigned int iScaling,
                                        unsigned int iEbin, double iE,
                                        double iDiffSens, double iDiffSensError,
                                        double iAngRes, double iERes,
                                        int iTreeEntry )
{
    // required sensitivity
    double req = getBasePerformance( iE );
    
    if( req < 1.e-30 )
    {
        cout << "WARNING: required sensitivity too small at " << TMath::Power( 10., iE ) << " TeV: " << req;
        cout << " (tree entry " << iTreeEntry << ")" << endl;
    }
    if( iDiffSens < 1.e-30 )
    {
        cout << "WARNING: sensitivity too small at " << TMath::Power( 10., iE ) << " TeV (PPUT Ebin " << iEbin << "): ";
        cout << iDiffSens;
        cout << " (tree entry " << iTreeEntry << ", " << iArrayID << ", " << iScaling << ")" << endl;
    }
    
    // required angular resolution
    double req_angres = getBasePerformanceAngularResolution( iE );
    // required energy resolution
    double req_eres = getBasePerformanceEnergyResolution( iE );
    
    /////////////////////////////////////////////////////////
    // require
    // - sensitivity > 0
    // - required sensitivity > 0
    // **otherwise PPUT will be zero!**
    if( iEbin < fSensitivityPPUT.size()
            && iScaling < fSensitivityPPUT[iEbin].size()
            && iArrayID < fSensitivityPPUT[iEbin][iScaling].size()
            && iDiffSens > 0.
            && req > 0. )
    {
        // add to PPUT term
        fSensitivityPPUT[iEbin][iScaling][iArrayID] *= req / iDiffSens;
        fSensitivityPPUT_N[iEbin][iScaling][iArrayID]++;
        
        fSensitivityPPUT_E[iEbin][iScaling][iArrayID] += iDiffSensError * iDiffSensError * req * req / iDiffSens / iDiffSens / iDiffSens / iDiffSens;
        
        // check angular resolution
        // (negative if requirement is not met)
        if( iAngRes > req_angres )
        {
            fAngularResolution[iEbin][iScaling][iArrayID] = -1.;
        }
        else if( fAngularResolution[iEbin][iScaling][iArrayID] > 0. )
        {
            fAngularResolution[iEbin][iScaling][iArrayID] = 1.6;
        }
        // check energy resolution
        if( iERes > req_eres )
        {
            fEnergyResolution[iEbin][iScaling][iArrayID] = -1.;
        }
        else if( fEnergyResolution[iEbin][iScaling][iArrayID] > 0. )
        {
            fEnergyResolution[iEbin][iScaling][iArrayID] = 0.35;
        }
    }
    else
    {
        // fudge: allow last bin at 100 TeV to be zero
        // (migh happen for 0.5h sensitivity files)
        if( TMath::Abs( iE - 2. ) > 1.e-5 )
        {
            fSensitivityPPUT[iEbin][iScaling][iArrayID] = 0.;
            fSensitivityPPUT_N[iEbin][iScaling][iArrayID] = 0.;
            fSensitivityPPUT_E[iEbin][iScaling][iArrayID] = 0.;
        }
    }
}

/*
 * fill a ranking for each energy bin of the PPUTs
 *
 */

void CPerformance::fillPPUTRanking()
{
    // energy bins
    for( unsigned int i = 0; i < fSensitivityPPUT.size(); i++ )
    {
        // pair consist of index for array scaling / array layouts
        multimap< double, pair< unsigned int, unsigned int > > iT;
        // array scaling
        for( unsigned int j = 0; j < fSensitivityPPUT[i].size(); j++ )
        {
            // array layouts
            for( unsigned int k = 0; k < fSensitivityPPUT[i][j].size(); k++ )
            {
                pair< unsigned int, unsigned int > iP;
                iP.first = j;
                iP.second = k;
                if( fArrayLayout[k].find( "S.3HB8" ) != string::npos
                        && fArrayLayout[k].find( "S.3HB8-NG" ) == string::npos )
                {
                    iT.insert( std::pair< double, pair< unsigned int, unsigned int > >( 0., iP ) );
                }
                else if( TMath::Abs( fSensitivityPPUT_E[i][iP.first][iP.second] ) > 1.e-4 )
                {
                    iT.insert( std::pair< double, pair< unsigned int, unsigned int > >( fSensitivityPPUT[i][j][k], iP ) );
                }
                else
                {
                    iT.insert( std::pair< double, pair< unsigned int, unsigned int > >( 0., iP ) );
                }
            }
        }
        // max PPUT
        std::multimap< double, pair< unsigned int, unsigned int > >::reverse_iterator it_RT = iT.rbegin();
        double iPPUT_max = it_RT->first;
        // now fill ranking
        // within 5%  = 2
        // within 10% = 3
        // within 25% = 4
        for( std::multimap< double, pair< unsigned int, unsigned int > >::iterator it_iT = iT.begin();
                it_iT != iT.end(); ++it_iT )
        {
            pair< unsigned int, unsigned int > iP = ( *it_iT ).second;
            if( iP.first < fSensitivityPPUTRanking[i].size()
                    && iP.second < fSensitivityPPUTRanking[i][iP.first].size() )
            {
                fSensitivityPPUTRanking[i][iP.first][iP.second] = 10;
                if( ( *it_iT ).first > 0.75 * iPPUT_max )
                {
                    fSensitivityPPUTRanking[i][iP.first][iP.second] = 4;
                }
                if( ( *it_iT ).first > 0.90 * iPPUT_max )
                {
                    fSensitivityPPUTRanking[i][iP.first][iP.second] = 3;
                }
                if( ( *it_iT ).first > 0.95 * iPPUT_max )
                {
                    fSensitivityPPUTRanking[i][iP.first][iP.second] = 2;
                }
                if( ( *it_iT ).first > 0.999 * iPPUT_max )
                {
                    fSensitivityPPUTRanking[i][iP.first][iP.second] = 1;
                }
            }
        }
    }
}

/*
 * perform final steps of analysis:
 * - PPUT calculation
 * - PPUT ranking
 *
 */
void CPerformance::terminate()
{
    char hname[200];
    fArrayLayoutLabels.clear();
    
    // energy bins
    for( unsigned int i = 0; i < fSensitivityPPUT.size(); i++ )
    {
        vector< string > iArrayLayoutLabel;
        
        int z = 0;
        int zbin = 0;
        // array scaling
        for( unsigned int j = 0; j < fSensitivityPPUT[i].size(); j++ )
        {
            // array layouts
            for( unsigned int k = 0; k < fSensitivityPPUT[i][j].size(); k++ )
            {
                /////////////////////////////////
                cout << "PPUT for array ";
                if( k < fArrayLayout.size() )
                {
                    cout << fArrayLayout[k] << ", scaling " << j + 1;
                    cout << ", energy bin " << i;
                    cout << " (averaged over " << fSensitivityPPUT_N[i][j][k] << " points)\t";
                }
                if( fSensitivityPPUT_N[i][j][k] > 0. )
                {
                    fSensitivityPPUT[i][j][k] = TMath::Power( fSensitivityPPUT[i][j][k], 1. / fSensitivityPPUT_N[i][j][k] );
                    
                    if( fSensitivityPPUT_N[i][j][k] > 2. )
                    {
                        fSensitivityPPUT_E[i][j][k] = fSensitivityPPUT[i][j][k] * sqrt( fSensitivityPPUT_E[i][j][k] );
                        fSensitivityPPUT_E[i][j][k] = 1. / fSensitivityPPUT_N[i][j][k]
                                                      * TMath::Power( fSensitivityPPUT[i][j][k], 1. / fSensitivityPPUT_N[i][j][k] - 1. )
                                                      * fSensitivityPPUT_E[i][j][k];
                    }
                    else
                    {
                        fSensitivityPPUT_E[i][j][k] = 0.;
                        
                    }
                    if( !isnan( fSensitivityPPUT_E[i][j][k] ) )
                    {
                        cout << "  PPUT: " << fSensitivityPPUT[i][j][k] << " +- " << fSensitivityPPUT_E[i][j][k];
                        
                        cout << endl;
                        
                        fSensitivity[i]->SetPoint( z, zbin + 1, fSensitivityPPUT[i][j][k] );
                        fSensitivity[i]->SetPointError( z, 0., fSensitivityPPUT_E[i][j][k] );
                        
                        fAngRes[i]->SetPoint( z, zbin + 1, fAngularResolution[i][j][k] );
                        fERes[i]->SetPoint( z, zbin + 1, fEnergyResolution[i][j][k] );
                    }
                    else
                    {
                        cout << "   PPUT_E NAN " << endl;
                        fSensitivity[i]->SetPoint( z, zbin + 1, 0. );
                        fSensitivity[i]->SetPointError( z, 0., 0. );
                        
                        fAngRes[i]->SetPoint( z, zbin + 1, -1. );
                        fERes[i]->SetPoint( z, zbin + 1, -1. );
                    }
                    fSensitivityPPUT_N[i][j][k] = -1.;
                }
                else
                {
                    fSensitivity[i]->SetPoint( z, zbin + 1, 0. );
                    fSensitivity[i]->SetPointError( z, 0., 0. );
                    
                    fAngRes[i]->SetPoint( z, zbin + 1, -1. );
                    fERes[i]->SetPoint( z, zbin + 1, -1. );
                    cout << endl;
                }
                // set array labels
                if( fSouthSite )
                {
                    sprintf( hname, "%s-%d", fArrayLayout[k].c_str(), j + 1 );
                }
                else
                {
                    sprintf( hname, "%s", fArrayLayout[k].c_str() );
                }
                string iT = hname;
                iArrayLayoutLabel.push_back( iT.substr( 2, string::npos ) );
                z++;
                zbin++;
            }
        }
        fArrayLayoutLabels.push_back( iArrayLayoutLabel );
    }
    fillPPUTRanking();
}


////////////////////////////////////////////////////////////////////////
/*
      plotting types:

      - "onAxis-N" onaxis sensitivities for NectarCam combinations (for the hardcoded observing time)
      - "onAxis-F" onaxis sensitivities for FlashCam combinations (for the hardcoded observing time)
*/
CPerformance* fillData( string iName, string iDataFile, string iPerformanceBaseFile,
                        vector< string > iArrayLayout, unsigned int iScaling,
                        int iColor, int iMarker, int iLowEnergyBin = 0,
                        float iCameraOffset = 0., string iPointingDirection = "",
                        int iObservingTime = 50 * 3600, bool iSouth = true, bool iKSP = false )
{
    // open file with sensitivity tree
    TFile* f = new TFile( iDataFile.c_str() );
    if( f->IsZombie() )
    {
        return 0;
    }
    
    cout << "reading sensitivity values from " << f->GetName() << endl;
    
    // sensitivity tree
    TTree* t = ( TTree* )f->Get( "IRFData" );
    if( !t )
    {
        return 0;
    }
    
    // performance class definition
    // (no checks if telescope counters exist: take care!)
    CPerformance* fPerformance = new CPerformance( iName, iSouth, iArrayLayout, iLowEnergyBin, iScaling, iColor, iMarker, iKSP );
    fPerformance->setCameraOffset( iCameraOffset );
    fPerformance->setObservingTime( iObservingTime );
    fPerformance->setPerformanceBase( 3, iPerformanceBaseFile );
    fPerformance->setPointingDirection( iPointingDirection );
    
    // tree access
    cIRFData* fData = new cIRFData( t );
    
    /////////////////////////
    // loop over all array layouts
    // search for all array layouts with similar names
    for( int i = 0; i < fData->fChain->GetEntries(); i++ )
    {
        fData->GetEntry( i );
        
        unsigned int iArrayID = fPerformance->getArrayID( fData->Array, fData->ObsTime_s, fData->Offset_deg );
        if( iArrayID == 9999 )
        {
            continue;
        }
        cout << "ARRAY " << fData->Array << "\t" << fData->Scaling << " (array id " << iArrayID << ")";
        cout << " [tree entry " << i << "] " << endl;
        cout << "--------------------------------" << endl;
        
        for( unsigned int e = 0; e < fPerformance->fEnergyBins_min.size(); e++ )
        {
            for( unsigned int ee = fPerformance->fEnergyBins_min[e]; ee <= fPerformance->fEnergyBins_max[e]; ee++ )
            {
                // differential sensitivity
                fPerformance->addSensitivityPoint( iArrayID, fData->Scaling - 1, e, fData->Energy_logTeV[ee],
                                                   fData->DiffSens[ee], fData->DiffSensError[ee],
                                                   fData->AngRes[ee], fData->ERes[ee], i );
            }
        }
    }
    // calculate PPUTs
    fPerformance->terminate();
    
    f->Close();
    
    return fPerformance;
}

/*
 * check if site is South or North
 */
bool isSouth( string iSite )
{
    if( iSite.find( "lapalma" ) != string::npos )
    {
        return false;
    }
    return true;
}

/*
 * write a latex table with all the results: final lines
 */
void writeLatexFileFooter( ofstream& os )
{
    os << "\\end{document}" << endl;
}

/*
 * write a latex table with all the results: final lines
 */
void writeLatexFileHeader( ofstream& os )
{
    // intro
    os << "\\documentclass[8pt]{scrartcl}" << endl;
    os << "\\usepackage[a4paper,landscape,scale=0.9]{geometry}" << endl;
    os << "\\usepackage{graphicx}" << endl;
    os << "\\usepackage{epstopdf}" << endl;
    os << "\\usepackage{longtable}" << endl;
    os << "\\usepackage{color}" << endl;
    os << "\\usepackage[pdftex,colorlinks=true,bookmarks=false,bookmarksopen=false]{hyperref}" << endl;
    
    os << "\\begin{document}" << endl;
}

/*
 * write a latex table and a txt table with all the results
 * (this function writes a single tex line)
*/
bool writeLatexTable( vector< CPerformance* > fPerformance, ostream& os, ostream& os_txt )
{

    for( unsigned int c = 0; c < fPerformance.size(); c++ )
    {
        os << "\\section*{" << fPerformance[c]->fName << "}" << endl;
        os << "\\begin{longtable}{l|c";
        if( fPerformance.size() > 0 )
        {
            for( unsigned int i = 0; i < fPerformance[0]->fEnergyLabels.size(); i++ )
            {
                os << "|c";
            }
        }
        os << "}" << endl;
        os_txt << fPerformance[c]->fName << endl;
        cout << "========================================================" << endl;
        cout << endl;
        // table title
        // (first latex then txt)
        os << "Array ";
        os << " & Scaling ";
        if( fPerformance.size() > 0 )
        {
            for( unsigned int i = 0; i < fPerformance[0]->fEnergyLabels.size(); i++ )
            {
                os << " & PPUT";
            }
        }
        os << "\\\\" << endl;
        
        // table title (line two)
        os << "Layout ";
        os << " & ";
        if( fPerformance.size() > 0 )
        {
            for( unsigned int i = 0; i < fPerformance[0]->fEnergyLabels.size(); i++ )
            {
                if( fPerformance[0]->fEnergyLabels[i] )
                {
                    string iT = fPerformance[0]->fEnergyLabels[i]->GetTitle();
                    os << " & " << iT.substr( iT.find( "=" ) + 1, string::npos );
                }
            }
        }
        os << "\\\\" << endl;
        os << "\\hline" << endl;
        os << "\\hline" << endl;
        // assume same amount of array layouts and scalings in all energy bins
        if( fPerformance[c]->fSensitivityPPUT.size() < 1 )
        {
            continue;
        }
        if( fPerformance[c]->fSensitivityPPUT[0].size() < 1 )
        {
            continue;
        }
        //////////////
        // table body
        
        // array scaling
        for( unsigned int j = 0; j < fPerformance[c]->fSensitivityPPUT[0].size(); j++ )
        {
            // array layouts
            for( unsigned int k = 0; k < fPerformance[c]->fSensitivityPPUT[0][j].size(); k++ )
            {
                // skip all HB8 and HF8 results which are not for scaling 1
                if( ( fPerformance[c]->fArrayLayout[k].find( "S.3HB8" ) != string::npos
                        || fPerformance[c]->fArrayLayout[k].find( "S.3HF8" ) != string::npos )
                        && j != 0 )
                {
                    continue;
                }
                // skip all HB8 and HF8 results which are not for NG arrays
                if( fPerformance[c]->fArrayLayout[k].find( "S.3HB8" ) != string::npos
                        && fPerformance[c]->fArrayLayout[k].find( "S.3HB8-NG" ) == string::npos )
                {
                    continue;
                }
                if( fPerformance[c]->fArrayLayout[k].find( "S.3HF8" ) != string::npos
                        && fPerformance[c]->fArrayLayout[k].find( "S.3HF8-NG" ) == string::npos )
                {
                    continue;
                }
                os << fPerformance[c]->fArrayLayout[k] << " & " << j + 1 << " ";
                os_txt << fPerformance[c]->fArrayLayout[k] << " \t " << j + 1;
                os_txt << "\t" << fPerformance[c]->fPointingDirection;
                
                // energy bins
                // (assume same number of energy bins everywhere!)
                for( unsigned int i = 0; i < fPerformance[c]->fSensitivityPPUT.size(); i++ )
                {
                    //                          if( TMath::Abs( fPerformance[c]->fSensitivityPPUT_E[i][j][k] - 0. ) < 1.e-5 )
                    if( isnan( fPerformance[c]->fSensitivityPPUT_E[i][j][k] )
                            ||  TMath::Abs( fPerformance[c]->fSensitivityPPUT_N[i][j][k] ) < 0.5 )
                    {
                        os_txt << "\t" << fixed << setprecision( 2 );
                        os_txt << 0.0;
                        os << " & " << fixed << setprecision( 2 ) << 0.0;
                        continue;
                    }
                    // 1%
                    if( fPerformance[c]->fSensitivityPPUTRanking[i][j][k] == 1 )
                    {
                        os << " & {\\color{green}\\textbf *" << fixed << setprecision( 2 );
                        os << fPerformance[c]->fSensitivityPPUT[i][j][k] << "*}";
                    }
                    // 5%
                    else if( fPerformance[c]->fSensitivityPPUTRanking[i][j][k] == 2 )
                    {
                        os << " & {\\color{green}\\textbf " << fixed << setprecision( 2 );
                        os << fPerformance[c]->fSensitivityPPUT[i][j][k] << "}";
                    }
                    // 10%
                    else if( fPerformance[c]->fSensitivityPPUTRanking[i][j][k] == 3 )
                    {
                        os << " & {\\color{magenta} " << fixed << setprecision( 2 );
                        os << fPerformance[c]->fSensitivityPPUT[i][j][k] << "}";
                    }
                    else
                    {
                        os << " & " << fixed << setprecision( 2 );
                        os << fPerformance[c]->fSensitivityPPUT[i][j][k];
                    }
                    os_txt << "\t" << fixed << setprecision( 2 );
                    os_txt << fPerformance[c]->fSensitivityPPUT[i][j][k];
                }
                os << " \\\\" << endl;
                os_txt << endl;
            }
        }
        os << "\\end{longtable}" << endl;
        os << "Colors: {\\color{green} PPUT within 5\\% of max PPUT (max PPUT is indicated by **)}, {\\color{magenta} PPUT within 10\\% of max PPUT}." << endl;
        os << "Note: S.HB8 results are only available for S.3HB-NG-1 (no averaging over different telescope/camera types)" << endl;
        os << "\\newpage" << endl;
    }
    
    
    return true;
}

/*
    draw performances and resolution curves

*/
TCanvas* drawPerformanceCurves( vector< CPerformance* > fPerformance, string iCanvasTitle, string iTelescopeCombination )
{
    //////////////////////////////
    // draw everything
    
    char hname[200];
    
    // A4     TCanvas *cC = new TCanvas( "cC", "sensitivity vs array layout", 1, 1, 297*3, 210*3 );
    sprintf( hname, "cC%s", iTelescopeCombination.c_str() );
    TCanvas* cC = new TCanvas( hname, "sensitivity vs array layout", 1, 1, ( int )( 1440 * 0.8 ), ( int )( 900 * 0.8 ) );
    cC->Divide( 2, 3, 0.02 );
    cC->Draw();
    
    // get minimum and maximum relative sensitivity
    double i_diffsens_min = 1.e9;
    double i_diffsens_max = -1.e9;
    for( unsigned int i = 0; i < fPerformance.size(); i++ )
    {
        for( unsigned int p = 0; p < fPerformance[i]->fSensitivity.size(); p++ )
        {
            if( fPerformance[i]->fSensitivity[p] )
            {
                if( fPerformance[i]->fSensitivity[p]->GetMinimum() < i_diffsens_min )
                {
                    i_diffsens_min = fPerformance[i]->fSensitivity[p]->GetMinimum() < i_diffsens_min;
                }
                if( fPerformance[i]->fSensitivity[p]->GetMaximum() > i_diffsens_max )
                {
                    i_diffsens_max = fPerformance[i]->fSensitivity[p]->GetMaximum() < i_diffsens_max;
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////
    // manually set min / max settings
    vector< double > i_DiffSens_min;
    vector< double > i_DiffSens_max;
    
    // plot
    for( unsigned int i = 0; i < fPerformance.size(); i++ )
    {
        for( unsigned int e = 0; e < fPerformance[i]->fSensitivityPPUT.size(); e++ )
        {
            TPad* p = ( TPad* )cC->cd( e + 1 );
            
            // sensitivity
            if( i != 0 )
            {
                fPerformance[i]->fSensitivity[e]->Draw( "p" );
            }
            else
            {
                p->SetTopMargin( 0.01 );
                p->SetRightMargin( 0.03 );
                p->SetBottomMargin( 0.30 );
                p->SetGridx( 0 );
                p->SetGridy( 0 );
                
                fPerformance[i]->fSensitivity[e]->SetMinimum( i_diffsens_min * 0.95 );
                fPerformance[i]->fSensitivity[e]->SetMaximum( i_diffsens_max * 1.05 );
                fPerformance[i]->fSensitivity[e]->SetMinimum( 0.25 );
                fPerformance[i]->fSensitivity[e]->SetMaximum( 1.65 );
                fPerformance[i]->fSensitivity[e]->Draw( "ap" );
                
                // set bin labels (from first histogram)
                TH1F* h = fPerformance[i]->fSensitivity[e]->GetHistogram();
                if( h )
                {
                    h->GetXaxis()->SetLabelSize( 0.08 );
                    for( int b = 0; b < fPerformance[i]->fSensitivity[e]->GetN(); b++ )
                    {
                        if( e < fPerformance[i]->fArrayLayoutLabels.size()
                                &&  b < ( int )fPerformance[i]->fArrayLayoutLabels[e].size() )
                        {
                            h->GetXaxis()->SetBinLabel( h->GetXaxis()->FindBin( b + 1 ), fPerformance[i]->fArrayLayoutLabels[e][b].c_str() );
                        }
                    }
                    sprintf( hname, "PPUT (relative to S.3HB1-%s-3, Av PT)", iTelescopeCombination.c_str() );
                    h->GetYaxis()->SetTitle( hname );
                    h->GetYaxis()->SetTitleOffset( 1.1 );
                }
                if( e < fPerformance[i]->fEnergyLabels.size() )
                {
                    fPerformance[i]->fEnergyLabels[e]->SetTextSize( fPerformance[i]->fEnergyLabels[e]->GetTextSize() * 1.5 );
                    fPerformance[i]->fEnergyLabels[e]->SetNDC();
                    fPerformance[i]->fEnergyLabels[e]->Draw();
                }
                
                TLine* iLP = new TLine( h->GetXaxis()->GetXmin(), 1., h->GetXaxis()->GetXmax(), 1. );
                iLP->SetLineStyle( 3 );
                iLP->Draw();
            }
            if( fPerformance[i]->fPlotResolutionRequirement )
            {
                fPerformance[i]->fAngRes[e]->Draw( "p" );
                fPerformance[i]->fERes[e]->Draw( "p" );
            }
        }
    }
    
    cC->Update();
    cC->cd();
    
    if( iCanvasTitle.size() > 0 )
    {
        TText* iTit = new TText( 0.505, 0.1, iCanvasTitle.c_str() );
        iTit->SetTextAngle( 90. );
        iTit->Draw();
    }
    
    cC->Update();
    return cC;
}

/*
 * read list of arrays from a text file
 *
 */
vector< string > readArrayList( string iArrayList, string iTelescopeCombination, bool iSouth = true )
{
    vector< string > iA;
    
    // read list of arrays
    ifstream is;
    is.open( iArrayList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error, array list file not found: " << endl;
        cout << iArrayList << endl;
        return iA;
    }
    
    string is_line;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 && is_line.substr( 0, 1 ) != "#" )
        {
            if( iSouth )
            {
                iA.push_back( is_line );
            }
            else
            {
                ostringstream iAs;
                // five scalings
                for( unsigned int s = 1; s <= 5; s++ )
                {
                    iAs.str( "" );
                    iAs << is_line << "-" << s;
                    iA.push_back( iAs.str() );
                }
            }
        }
    }
    
    is.close();
    
    // add telescope combination
    for( unsigned int i = 0; i < iA.size(); i++ )
    {
        if( iTelescopeCombination.size() > 0 && iA[i].find( "-" ) == string::npos )
        {
            iA[i] = iA[i] + "-" + iTelescopeCombination;
        }
    }
    
    return iA;
}

string getSetNameString( float iCameraOffset = 0., int iObservingTime = 50 * 3600 )
{
    ostringstream iA;
    
    iA << fixed << setprecision( 1 ) << iObservingTime / 3600. << "h";
    if( iCameraOffset < 1.e-2 )
    {
        iA << "-onAxis";
    }
    else
    {
        iA << fixed << setprecision( 1 ) << iCameraOffset << "deg";
    }
    
    return iA.str();
}

/*
 *   main plotting function
 *
 *   iArrayList = list of array layouts
 *
 */
TCanvas* plot( ofstream& os, ofstream& os_txt, string iArrayList, string iSite = "paranal-20deg", string iCanvasTitle = "",
               string iTelescopeCombination = "ND", string iTelescopeBaseLineCombination = "ND",
               int iSubArrayID = 0, float iCameraOffset = 0., int iObservingTime = 50 * 3600, bool iKSP = false )
{
    // get list of arrays
    vector< string > iArrayLayout = readArrayList( iArrayList, iTelescopeCombination, isSouth( iSite ) );
    
    // data directory
    string iDDIR = "/Users/maierg/Experiments/CTA/Analysis/sensitivity/data/";
    
    vector< CPerformance* > fPerformance;
    ostringstream iIRFFile;
    ostringstream iLatexSectionTitle;
    ///////////////////////////////////////////////////
    // Parnal 20 deg full arrays
    if( iSite == "paranal-20deg" )
    {
        // Performance base file
        ostringstream iPerformanceBaseFile;
        iPerformanceBaseFile << iDDIR;
        iPerformanceBaseFile << "/S.20deg.20160313/DESY.d20160304.V5.ID";
        iPerformanceBaseFile << iSubArrayID;
        iPerformanceBaseFile << "NIM2.prod3-paranalp05-NN.S.3HB8-";
        iPerformanceBaseFile << iTelescopeBaseLineCombination << "-1.";
        iPerformanceBaseFile << iObservingTime << "s.root";
        
        // Average Pointing
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-paranalp05-NN-ID" << iSubArrayID << "-Average-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "Paranal (20deg zenith, Average pointing ";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << "): " << iCanvasTitle;
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 5, 1, 21, 0, iCameraOffset, "Average",
                                          iObservingTime, isSouth( iSite ), iKSP ) );
        if( !fPerformance.back() )
        {
            return 0;
        }
        fPerformance.back()->setPlotResolutionRequirement();
        
        // Pointing North
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-paranalp05-NN-ID" << iSubArrayID << "-0deg-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "Paranal (20deg zenith, Pointing North";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << "): " << iCanvasTitle;
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 5, 2, 24, 0, iCameraOffset, "North",
                                          iObservingTime, isSouth( iSite ) ) );
                                          
        // Pointing South
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-paranalp05-NN-ID" << iSubArrayID << "-180deg-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "Paranal (20deg zenith, Pointing South";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << "): " << iCanvasTitle;
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 5, 4, 24, 0, iCameraOffset, "South",
                                          iObservingTime, isSouth( iSite ) ) );
    }
    ///////////////////////////////////////////////////
    // Parnal 40 deg full arrays
    else if( iSite == "paranal-40deg" )
    {
        // Performance base file
        ostringstream iPerformanceBaseFile;
        iPerformanceBaseFile << iDDIR;
        iPerformanceBaseFile << "/S.40deg.20160313/DESY.d20160304.V5.ID";
        iPerformanceBaseFile << iSubArrayID;
        iPerformanceBaseFile << "NIM2.prod3-paranalp05-40deg-NN.S.3HB1-";
        iPerformanceBaseFile << iTelescopeBaseLineCombination << "-3.";
        iPerformanceBaseFile << iObservingTime << "s.root";
        
        // Average Pointing
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-paranalp05-40deg-NN-ID" << iSubArrayID << "-Average-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "Paranal (40deg zenith, Average pointing";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << "): " << iCanvasTitle;
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 5, 1, 21, 1, iCameraOffset, "Average",
                                          iObservingTime, isSouth( iSite ) ) );
        fPerformance.back()->setPlotResolutionRequirement();
        
        // Pointing North
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-paranalp05-40deg-NN-ID" << iSubArrayID << "-0deg-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "Paranal (40deg zenith, Pointing North";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << "): " << iCanvasTitle;
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 5, 2, 24, 1, iCameraOffset, "North",
                                          iObservingTime, isSouth( iSite ) ) );
                                          
        // Pointing South
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-paranalp05-40deg-NN-ID" << iSubArrayID << "-180deg-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "Paranal (40deg zenith, Pointing South";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << "): " << iCanvasTitle;
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 5, 4, 24, 1, iCameraOffset, "South",
                                          iObservingTime, isSouth( iSite ) ) );
    }
    ///////////////////////////////////////////////////
    // La Palma 20 deg
    else if( iSite == "lapalma-20deg" )
    {
        // Performance base file
        ostringstream iPerformanceBaseFile;
        iPerformanceBaseFile << iDDIR;
        iPerformanceBaseFile << "/N.20deg.20160313/DESY.d20160304.V5.ID";
        iPerformanceBaseFile << iSubArrayID;
        iPerformanceBaseFile << "NIM2.prod3-LaPalmap05-NN.N.3AL4M15-3-F.";
        iPerformanceBaseFile << iObservingTime << "s.root";
        
        // Average Pointing
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-LaPalmap05-NN-ID" << iSubArrayID << "-Average-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "La Palma (20deg zenith, Average pointing ";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << ")";
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 1, 1, 21, 0, iCameraOffset, "Average",
                                          iObservingTime, isSouth( iSite ) ) );
        if( !fPerformance.back() )
        {
            return 0;
        }
        fPerformance.back()->setPlotResolutionRequirement();
        
        // Pointing North
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-LaPalmap05-NN-ID" << iSubArrayID << "-0deg-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "La Palma (20deg zenith, Pointing North";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << ")";
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 1, 2, 24, 0, iCameraOffset, "North",
                                          iObservingTime, isSouth( iSite ) ) );
                                          
        // Pointing South
        iIRFFile.str( "" );
        iIRFFile << iDDIR << "NS.SensitivityTrees.20160313/";
        iIRFFile << "prod3-LaPalmap05-NN-ID" << iSubArrayID << "-180deg-sensitivityTree.root";
        
        iLatexSectionTitle.str( "" );
        iLatexSectionTitle <<  "La Palma (20deg zenith, Pointing South";
        iLatexSectionTitle << getSetNameString( iCameraOffset, iObservingTime ) << ")";
        fPerformance.push_back( fillData( iLatexSectionTitle.str(), iIRFFile.str(), iPerformanceBaseFile.str(),
                                          iArrayLayout, 1, 4, 24, 0, iCameraOffset, "South",
                                          iObservingTime, isSouth( iSite ) ) );
    }
    
    writeLatexTable( fPerformance, os, os_txt );
    
    return drawPerformanceCurves( fPerformance, iCanvasTitle, iTelescopeBaseLineCombination );
    
}

/*
 * return a string with the Canvas title printed in the centre of the Canvas
 *
 * allowed  layout configurations:
 * - baseline
 * - descoped
 * - staging1, staging2, .., staging7
*/
string getCanvasTitle( string layout_configurations, string telescope_combinations, bool isSouth = true )
{
    std::ostringstream i_title;
    i_title << layout_configurations << " ";
    if( telescope_combinations.size() > 0 )
    {
        if( telescope_combinations.substr( 0, 1 ) == "N" )
        {
            i_title << "[MST-N],";
        }
        else if( telescope_combinations.substr( 0, 1 ) == "F" )
        {
            i_title << "[MST-F],";
        }
    }
    else if( telescope_combinations.size() == 0 )
    {
        i_title << "[MST-N,MST-F], [1mDC,GCT]";
    }
    if( telescope_combinations.size() > 1 )
    {
        if( telescope_combinations.substr( 1, 1 ) == "G" )
        {
            i_title << "[GCT]";
        }
        else if( telescope_combinations.substr( 1, 1 ) == "D" )
        {
            i_title << "[1mDC]";
        }
    }
    else if( telescope_combinations.size() == 1 )
    {
        if( isSouth )
        {
            i_title << "[1mDC,GCT]";
        }
    }
    
    return i_title.str();
}

/*
 * get name of subarray for title pages
 */
string getSubArrayName( int iSubArrayID )
{
    if( iSubArrayID == 1 )
    {
        return "LST arrays";
    }
    else if( iSubArrayID == 2 )
    {
        return "MST arrays";
    }
    else if( iSubArrayID == 3 )
    {
        return "SST arrays";
    }
    else if( iSubArrayID == 4 )
    {
        return "MST+SST arrays";
    }
    return "";
}


/*
 * print the last page for the pdf document
 *
*/
void print_finalPage( string iPrint )
{
    char hname[200];
    if( iPrint.size() > 0 )
    {
        TCanvas* cEmpty = new TCanvas( "b", "b", 1, 1, ( int )( 1440 * 0.8 ), ( int )( 900 * 0.8 ) );
        cEmpty->Draw();
        sprintf( hname, "%s.pdf)", iPrint.c_str() );
        cEmpty->Print( hname );
    }
}

/*
 * print the first page for the pdf document
 *
 * contains title text and date
 *
 */
void print_titlePage( string iPrint, string iText, int iSubArrayID, float iCameraOffset, int iObservingTime )
{
    char hname[200];
    if( iPrint.size() > 0 )
    {
        TCanvas* cEmpty = 0;
        cEmpty = new TCanvas( "a", "b", 1, 1, ( int )( 1440 * 0.8 ), ( int )( 900 * 0.8 ) );
        cEmpty->Draw();
        
        TPaveText* iPT = new TPaveText( 0.3, 0.5, 0.7, 0.9, "NB" );
        iPT->SetFillStyle( 0 );
        TText* t = iPT->AddText( iText.c_str() );
        t->SetTextColor( 50 );
        t->SetTextSize( t->GetTextSize() * 2 );
        t = iPT->AddText( "PPUT Summary Plots" );
        t->SetTextColor( 50 );
        t->SetTextSize( t->GetTextSize() * 2 );
        t->SetTextColor( 50 );
        if( iCameraOffset < 1.e-2 )
        {
            sprintf( hname, "on axis %s %.1fh", getSubArrayName( iSubArrayID ).c_str(), ( float )iObservingTime / 3600. );
        }
        else
        {
            sprintf( hname, "offaxis %.1f deg %s %.1fh", iCameraOffset, getSubArrayName( iSubArrayID ).c_str(), ( float )iObservingTime / 3600. );
        }
        t = iPT->AddText( hname );
        t->SetTextColor( 50 );
        t->SetTextFont( 52 );
        t = iPT->AddText( "DESY analysis - 2016-04-18" );
        t->SetTextColor( 50 );
        t->SetTextFont( 52 );
        
        iPT->Draw();
        
        // explanation
        TPaveText* iPE = new TPaveText( 0.3, 0.1, 0.7, 0.5, "NP" );
        iPE->SetTextFont( 52 );
        iPE->SetFillStyle( 0 );
        t = iPE->AddText( "PPUT calculation relative to the sensitivity" );
        t = iPE->AddText( "of the corresponding S.3HB1-3 array" );
        t = iPE->AddText( "and for the given energy range." );
        t = iPE->AddText( "Black markers: average pointing" );
        t = iPE->AddText( "Red markers: pointing north" );
        t->SetTextColor( 2 );
        t = iPE->AddText( "Blue markers: pointing south" );
        t->SetTextColor( 4 );
        t = iPE->AddText( "'x' markers: angular resolution requirement met" );
        t = iPE->AddText( "'+' markers: energy resolution requirement met" );
        
        iPE->Draw();
        
        sprintf( hname, "%s.pdf(", iPrint.c_str() );
        cEmpty->Print( hname );
    }
}

/*
 * main plotting function
 *
 * works for paranal 20 and 40 deg
 *
*/
void plotParanal( string iSite = "paranal-20deg", int iSubArrayID = 0,
                  float iCameraOffset = 0., int iObservingTime = 50 * 3600,
                  bool iShort = false )
{
    TCanvas* c = 0;
    
    vector< string > layout_configurations;
    //     layout_configurations.push_back( "baseline" );
    layout_configurations.push_back( "KSPSouth" );
    bool iKSP = true;
    if( !iShort )
    {
        if( isSouth( iSite ) )
        {
            layout_configurations.push_back( "descoped" );
        }
        layout_configurations.push_back( "staging1" );
        layout_configurations.push_back( "staging2" );
        layout_configurations.push_back( "staging3" );
        layout_configurations.push_back( "staging4" );
        if( isSouth( iSite ) )
        {
            layout_configurations.push_back( "staging5" );
            layout_configurations.push_back( "staging6" );
            layout_configurations.push_back( "staging7" );
        }
    }
    
    // baseline combinations: PPUTs are calculated relative to this value
    vector< string > telescope_combinations;
    vector< string > telescope_baselinecombinations;
    // Paranal Site
    if( isSouth( iSite ) )
    {
        telescope_combinations.push_back( "NG" );
        telescope_baselinecombinations.push_back( "NG" );
        if( !iShort )
        {
            telescope_combinations.push_back( "ND" );
            telescope_baselinecombinations.push_back( "ND" );
            telescope_combinations.push_back( "FG" );
            telescope_baselinecombinations.push_back( "FG" );
            telescope_combinations.push_back( "FD" );
            telescope_baselinecombinations.push_back( "FD" );
            telescope_combinations.push_back( "N" );
            telescope_baselinecombinations.push_back( "NG" );
            telescope_combinations.push_back( "F" );
            telescope_baselinecombinations.push_back( "FG" );
            telescope_combinations.push_back( "" );
            telescope_baselinecombinations.push_back( "FG" );
        }
    }
    // LaPalma Site
    else
    {
        telescope_combinations.push_back( "F" );
        telescope_baselinecombinations.push_back( "F" );
    }
    
    // basename for all output files
    ostringstream iOutputFileName;
    iOutputFileName << iSite << "-PPUT-ID" << iSubArrayID << "ObsTime";
    iOutputFileName << getSetNameString( iCameraOffset, iObservingTime );
    
    // open a latex file for the summary tables
    ostringstream iLatexFileName;
    iLatexFileName << iOutputFileName.str() << "Table.tex";
    ofstream os;
    os.open( iLatexFileName.str().c_str() );
    if( !os )
    {
        cout << "failed opening latex file: " << iLatexFileName.str() << endl;
        return;
    }
    writeLatexFileHeader( os );
    // open a text file for summary tables
    ostringstream iTXTFileName;
    iTXTFileName << iOutputFileName.str() << "Table.txt";
    ofstream os_txt;
    os_txt.open( iTXTFileName.str().c_str() );
    if( !os_txt )
    {
        cout << "failed opening latex file: " << iTXTFileName.str() << endl;
        return;
    }
    // open a pdf file for printing
    // (all plots are printed in a single pdf file)
    ostringstream iPDFName;
    iPDFName << iOutputFileName.str();
    print_titlePage( iPDFName.str(), iSite, iSubArrayID, iCameraOffset, iObservingTime );
    
    ///////////////////////////////////////////////////
    // loop over all configuration and print results
    for( unsigned int i = 0; i < layout_configurations.size(); i++ )
    {
        for( unsigned int j = 0; j < telescope_combinations.size(); j++ )
        {
            std::ostringstream i_list;
            if( isSouth( iSite ) )
            {
                i_list << "/Users/maierg/Experiments/CTA/Analysis/sensitivity/ArrayLayoutAnalysis/subArray.prod3.";
            }
            else
            {
                i_list << "/Users/maierg/Experiments/CTA/Analysis/sensitivity/ArrayLayoutAnalysis/subArray.prod3N.";
            }
            i_list << layout_configurations[i] << ".list";
            
            c = plot( os, os_txt, i_list.str(), iSite, getCanvasTitle( layout_configurations[i], telescope_combinations[j], isSouth( iSite ) ),
                      telescope_combinations[j], telescope_baselinecombinations[j], iSubArrayID, iCameraOffset, iObservingTime, iKSP );
            if( c )
            {
                string i_print = iPDFName.str() + ".pdf";
                c->Print( i_print.c_str() );
            }
            else
            {
                return;
            }
        }
    }
    
    // summmary latex file: closing lines
    writeLatexFileFooter( os );
    // print final page of pdf document
    print_finalPage( iPDFName.str() );
    os.close();
    os_txt.close();
}

/*
 * plot PPUT for the different camera configurations in one plot
 *
 */
void plot_PPUT_CameraConfigurations( string iPPUT_Txt_file )
{
    vector< string > array_layout;
    array_layout.push_back( "HB1" );
    array_layout.push_back( "HI1" );
    vector< string > telescope_combinations;
    vector< string > label_name;
    vector< int > color_ID;
    telescope_combinations.push_back( "NG" );
    label_name.push_back( "NectarCam/GCT" );
    color_ID.push_back( 8 );
    telescope_combinations.push_back( "ND" );
    label_name.push_back( "NectarCam/1m-DC" );
    color_ID.push_back( 9 );
    telescope_combinations.push_back( "FG" );
    label_name.push_back( "FlashCam/GCT" );
    color_ID.push_back( 20 );
    telescope_combinations.push_back( "FD" );
    label_name.push_back( "FlashCam/1m-DC" );
    color_ID.push_back( 46 );
    
    vector< TGraphErrors* > fPPUT;
    vector< int > z;   // simple counter
    vector< float > x_ffset;
    for( unsigned int i = 0; i < telescope_combinations.size(); i++ )
    {
        fPPUT.push_back( new TGraphErrors( 1 ) );
        fPPUT.back()->SetMarkerStyle( 20 + i );
        fPPUT.back()->SetMarkerSize( 1.5 );
        fPPUT.back()->SetMarkerColor( color_ID[i] );
        fPPUT.back()->SetName( "" );
        z.push_back( 1 );
        fPPUT.back()->SetPoint( 0, 0., 1.e5 );
        x_ffset.push_back( -0.2 + 0.1 * i );
    }
    
    // open txt file
    ifstream is;
    is.open( iPPUT_Txt_file.c_str() );
    if( !is )
    {
        cout << "Error opening txt file " << iPPUT_Txt_file << endl;
        return;
    }
    
    ///////////////////////
    // read txt file
    string i_line;
    string t;
    while( getline( is, i_line ) )
    {
        istringstream is_stream( i_line );
        
        if( i_line.size() > 0 && i_line.substr( 0, 1 ) != "S" )
        {
            continue;
        }
        is_stream >> t;
        for( unsigned int i = 0; i < telescope_combinations.size(); i++ )
        {
            for( unsigned int a = 0; a < array_layout.size(); a++ )
            {
                string iA = "S.3" + array_layout[a] + "-" + telescope_combinations[i];
                if( iA != t )
                {
                    continue;
                }
                is_stream >> t;
                // scaling
                for( int s = 0; s < 5; s++ )
                {
                    if( atoi( t.c_str() ) != s + 1 )
                    {
                        continue;
                    }
                    // pointing
                    is_stream >> t;
                    if( t != "Average" )
                    {
                        continue;
                    }
                    
                    is_stream >> t;
                    
                    // a new point to the TGraph
                    fPPUT[i]->SetPoint( z[i], z[i] + x_ffset[i], atof( t.c_str() ) );
                    
                    z[i]++;
                    
                    cout << iA << "-" << s + 1 << ": " << t << " (" << z[i] << ")" << endl;
                }
            }
        }
    }
    is.close();
    
    // now draw everything
    TCanvas* cPPUT = new TCanvas( "cPPUT", "PPUT vs telescope/camera type", 10, 10, 800, 400 );
    cPPUT->SetGridx( 0 );
    cPPUT->SetGridy( 0 );
    cPPUT->SetTopMargin( 0.05 );
    cPPUT->SetRightMargin( 0.03 );
    cPPUT->SetBottomMargin( 0.25 );
    
    TLegend* iL = new TLegend( 0.2, 0.7, 0.4, 0.9 );
    
    for( unsigned int i = 0; i < fPPUT.size(); i++ )
    {
        if( !fPPUT[i] )
        {
            continue;
        }
        
        if( i == 0 )
        {
            fPPUT[i]->SetMinimum( 0.80 );
            fPPUT[i]->SetMaximum( 1.20 );
            fPPUT[i]->SetTitle( "" );
            fPPUT[i]->SetName( "" );
            fPPUT[i]->Draw( "ap" );
            
            // set bin labels (from first histogram)
            TH1F* h = fPPUT[i]->GetHistogram();
            if( h )
            {
                h->GetXaxis()->SetLabelSize( 0.08 );
                unsigned int b = 1;
                for( unsigned int s = 0; s < 5; s++ )
                {
                    for( unsigned int i = 0; i < array_layout.size(); i++ )
                    {
                        ostringstream i_BinLabel;
                        i_BinLabel << "S.3" << array_layout[i] << "-" << s + 1;
                        h->GetXaxis()->SetBinLabel( h->GetXaxis()->FindBin( b ), i_BinLabel.str().c_str() );
                        b++;
                    }
                }
                h->GetXaxis()->SetLabelSize( 0.07 );
                h->GetYaxis()->SetTitle( "PPUT [20 GeV - 125 TeV] (relative to S.3HB1-3)" );
                h->GetYaxis()->SetTitleOffset( 1.1 );
                TLine* iLP = new TLine( h->GetXaxis()->GetXmin(), 1., h->GetXaxis()->GetXmax(), 1. );
                iLP->SetLineStyle( 3 );
                iLP->Draw();
            }
        }
        else
        {
            fPPUT[i]->Draw( "p" );
        }
        
        iL->AddEntry( fPPUT[i], label_name[i].c_str(), "p" );
    }
    iL->Draw();
}


/*
 * consistency check of results in performance tree
 *
 * print arrays-scalings with unusual values or number of entries
 * in sensitivity / angular / energy resolution histograms
 *
 * Input:
 * - IRF data tree file name
 * - threshold on minimum number of entries for printout
 *
 * Hardwired number of energy bins
 */
void checkIRFDataFile( string iIRFDataTree, int iPrintTh = 1000, int ObsTime_s = 180000, float Offset_deg = 0. )
{
    TFile* f = new TFile( iIRFDataTree.c_str() );
    if( f->IsZombie() )
    {
        cout << "Error: opening file " << iIRFDataTree << endl;
        return;
    }
    TTree* t = ( TTree* )f->Get( "IRFData" );
    if( !t )
    {
        cout << "IRFData tree not found" << endl;
        return;
    }
    cIRFData* fData = new cIRFData( t );
    
    for( int i = 0; i < fData->fChain->GetEntries(); i++ )
    {
        fData->GetEntry( i );
        
        if( fData->ObsTime_s != ObsTime_s )
        {
            continue;
        }
        if( TMath::Abs( fData->Offset_deg - Offset_deg ) > 1.e-2 )
        {
            continue;
        }
        
        // count valid sensitivity points
        int z = 0;
        for( int j = 0; j < 21; j++ )
        {
            if( fData->DiffSens[j] > 0. )
            {
                z++;
            }
        }
        if( z < iPrintTh )
        {
            cout << "DiffSens: " << fData->Array << ": " << z << endl;
        }
    }
    cout << "(checked " << fData->fChain->GetEntries() << " entries)" << endl;
}
