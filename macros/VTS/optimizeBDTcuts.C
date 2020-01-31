/* \file optimizeBDTcuts.C
 *
 * \brief optimize BDT cuts
 *
 *  
 */

R__ADD_LIBRARY_PATH($EVNDISPSYS/lib)
R__LOAD_LIBRARY(libVAnaSum.so)

#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"

using namespace std;

void help()
{
    cout << endl;
    cout << "optimize the BDT cuts" << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "(requires full path (no ENVIRONMENTAL variables in path name)" << endl;
    cout << endl;
}

void optimizeBDTcuts( string particleraterootfile, 
                      const string weightFileDir, string weightFileName = "mva",
                      string MVAName = "BDT", unsigned int MVACounter = 0,
                      double observing_time_h = 20., 
                      int weightFileIndex_Emin = 0, int weightFileIndex_Emax = 2,
                      int weightFileIndex_Zmin = 0, int weightFileIndex_Zmax = 2,
                      double significance = 3., int min_events = 3,
                      bool iPlotEfficiencyPlots = false, bool iPlotOptimisationResults = true,
                      string iWriteTMACuts = "" )
{
    gSystem->Load( "$EVNDISPSYS/lib/libVAnaSum.so" ) ;

    // fixed parameters
    double min_sourcestrength_CU = 0.00001;
    // conversion parameter for particle rates and
    // observation time (h-->s)
    double timeparticlerate = 3600.;
    // use energy bins defined for training
    // (note that all other options will 
    // lead to a misassignment of xml files 
    // to energy bins!)
    double energyStepSize = 0.25;
    // typical alpha
    double min_backgroundrateratio = 1. / 5.;
    double min_backgroundevents = 0.;
    double signalefficiency = 0.90;
    // systematic cut
    double iSignaltoMinBackgroundRateRatio = 0.;
    
    VTMVAEvaluator a;
    a.setTMVAMethod( MVAName, MVACounter );
    a.setPrintPlotting( false );
    a.setPlotEfficiencyPlotsPerBin( iPlotEfficiencyPlots );
    a.setParticleNumberFile( particleraterootfile, timeparticlerate );

    a.setSensitivityOptimizationParameters( significance, min_events, observing_time_h, 
                            min_backgroundrateratio, min_backgroundevents, iSignaltoMinBackgroundRateRatio );
    a.setSensitivityOptimizationFixedSignalEfficiency( signalefficiency );
    a.setSensitivityOptimizationSourceStrength( min_sourcestrength_CU );
    // default: smoothing happens in makeEffectiveArea
    a.setSmoothAndInterpolateMVAValues( false );


    ostringstream iFullWeightFileName;
    iFullWeightFileName << weightFileDir << "/" << weightFileName;
    a.initializeWeightFiles( iFullWeightFileName.str().c_str(), 
                             weightFileIndex_Emin, weightFileIndex_Emax, 
                             weightFileIndex_Zmin, weightFileIndex_Zmax, 
                             energyStepSize, "VTS", "UseAveragedCounts" );
    // (do not change "UseAveragedCounts", interpolation is less robust)

    if( iPlotOptimisationResults )
    {
        a.plotSignalAndBackgroundEfficiencies( false, 1.e-2, -0.2, 0.4 );
    }
    
    a.printOptimizedMVACutValues();
    if( iWriteTMACuts.size() > 0 )
    {
        a.writeOptimizedMVACutValues( iWriteTMACuts );
    }

}

/* 
 * compare optimisation results from different MVAs
 *
 */
void plotCompare( unsigned int iZeBin = 0 )
{
     vector< string > iFileName;
     iFileName.push_back( "TL5035MA20_ID0_NTel2-PointSource-Soft-MVA-Preselection_RE_N8_MC_MVA.default.root" );
     /*iFileName.push_back( "TL5035MA20_ID0_NTel2-PointSource-Soft-MVA-Preselection_RE_N8_MC_MVA.noEnergyBins.root" );
     iFileName.push_back( "TL5035MA20_ID0_NTel2-PointSource-Soft-MVA-Preselection_RE_N8_MC_MVA.lessEnergyBins.root" );
     iFileName.push_back( "TL5035MA20_ID0_NTel2-PointSource-Soft-MVA-Preselection_RE_N8_MC_MVA.moreEnergyBins.root" );  */

     TCanvas *cTMVA = new TCanvas( "cTMVA", "TMVA optimised cut", 10, 10, 800, 800 );
     cTMVA->SetGridx( 0 );
     cTMVA->SetGridy( 0 );
     TCanvas *cSEff = new TCanvas( "cSEff", "signal efficiency", 810, 10, 800, 800 );
     cSEff->SetGridx( 0 );
     cSEff->SetGridy( 0 );
     TCanvas *cBEff = new TCanvas( "cBEff", "background efficiency", 810, 810, 800, 800 );
     cBEff->SetGridx( 0 );
     cBEff->SetGridy( 0 );
     cBEff->SetLogy( 1 );

     char hname[200];
     for( unsigned int i = 0; i < iFileName.size(); i++ )
     {
         TFile *iF = new TFile( iFileName[i].c_str() );
         if( iF->IsZombie() )
         {
             continue;
         }
         cout << "reading " << iF->GetName() << endl;
         cTMVA->cd();
         sprintf( hname, "TMVACutValue_ze%u", iZeBin );
         TGraphAsymmErrors *iG = (TGraphAsymmErrors*)iF->Get( hname );
         if( iG )
         {
             iG->SetTitle( "" );
             iG->SetMarkerColor( i + 1 );
             iG->SetLineColor( i + 1 );
             iG->SetMarkerStyle( 20 + i );
             iG->SetMinimum( -0.1 );
             iG->SetMaximum( 0.5 );
             if( i == 0 )
             {
                 iG->Draw( "ap" );
                 iG->GetHistogram()->SetXTitle( "log_{10} energy (TeV)" );
                 iG->GetHistogram()->SetYTitle( "MVA cut value" );
             }
             else
             {
                iG->Draw( "p" );
             }
         }
         else
         {
             cout << "MVA cut graph not found " << iF->GetName() << endl;
             cout << "\t " << endl;
         }
         cSEff->cd();
         sprintf( hname, "SignalEfficiency_ze%u", iZeBin );
         iG = (TGraphAsymmErrors*)iF->Get( hname );
         if( iG )
         {
             iG->SetTitle( "" );
             iG->SetMarkerColor( i + 1 );
             iG->SetLineColor( i + 1 );
             iG->SetMarkerStyle( 20 + i );
             iG->SetMinimum( 0. );
             iG->SetMaximum( 1. );
             if( i == 0 )
             {
                 iG->Draw( "apl" );
                 iG->GetHistogram()->SetXTitle( "log_{10} energy (TeV)" );
                 iG->GetHistogram()->SetYTitle( "signal efficiency" );
             }
             else
             {
                iG->Draw( "pl" );
             }
         }
         cBEff->cd();
         sprintf( hname, "BackgroundEfficiency_ze%u", iZeBin );
         iG = (TGraphAsymmErrors*)iF->Get( hname );
         if( iG )
         {
             iG->SetTitle( "" );
             iG->SetMarkerColor( i + 1 );
             iG->SetLineColor( i + 1 );
             iG->SetMarkerStyle( 20 + i );
             iG->SetMinimum( 5.e-5 );
             iG->SetMaximum( 1. );
             if( i == 0 )
             {
                 iG->Draw( "apl" );
                 iG->GetHistogram()->SetXTitle( "log_{10} energy (TeV)" );
                 iG->GetHistogram()->SetYTitle( "background efficiency" );
             }
             else
             {
                iG->Draw( "pl" );
             }
         }
         iF->Close();
     }
}

/*  smooth MVA eval
 *
 */
void smoothMVA( string iFileName, string iGraphName = "TMVACutValue_ze",
                int iZeBin = 0, 
                double dEres = 0.25, int iMVAiter = 1,
                double energy_constant_mva = 1.,
                double acceptMax = 0.5,
                double plotMax = 0.3, double plotMin = 0. )
{

     TCanvas *cTMVA = new TCanvas( "cTMVA", "TMVA optimised cut", 10, 10, 800, 800 );
     cTMVA->SetGridx( 0 );
     cTMVA->SetGridy( 0 );

     TFile *iF = new TFile( iFileName.c_str() );
     if( iF->IsZombie() )
     {
         return;
     }
     char hname[200];
     sprintf( hname, "%s%d", iGraphName.c_str(), iZeBin );
     TGraphAsymmErrors *iG = (TGraphAsymmErrors*)iF->Get( hname );
     iG->SetTitle( "" );
     iG->SetMarkerColor( 1 );
     iG->SetLineColor( 1 );
     iG->SetMarkerStyle( 20 );
     iG->SetMinimum( plotMin );
     iG->SetMaximum( plotMax );
     iG->Draw( "ap" );
     iG->GetHistogram()->SetXTitle( "log_{10} energy (TeV)" );
     iG->GetHistogram()->SetYTitle( "mva optimised cut" );

     double x = 0.;
     double y = 0.;
     int z = 0;
     // first point is lowest energy
     iG->GetPoint( 0, x, y );
     double e_min = x - iG->GetErrorXlow( 0 );
     iG->GetPoint( iG->GetN()-1, x, y );
     double e_max = x + iG->GetErrorXhigh( iG->GetN()-1 );

     // TGraph::Eval
     TProfile *iSmoothedEval = new TProfile( "smoothedEval", "", 100., e_min, e_max, -1., 1. );
     iSmoothedEval->SetLineColor(2);
     // TGraph::Eval + Gaussian bin width
     TProfile *iSmoothedEvalGauss = new TProfile( "smoothedEvalGauss", "", 100., e_min, e_max, -1., 1. );
     iSmoothedEvalGauss->SetLineColor(3);

     double iV = 0.;
     double e_log10_G = 0.;
/*     for( unsigned int i = 0; i < 1000; i++ )
     {
          double e_log10 = gRandom->Uniform( e_min, e_max );

          iSmoothedEval->Fill( e_log10, iG->Eval( e_log10 ) );

          double iM = 0.;
          double iN = 0.;
          for( unsigned int j = 0; j < 1000; j++ )
          {
              iV = iG->Eval( e_log10 + gRandom->Gaus( 0., dEres ) );
              if( iV < acceptMax )
              {
                  iM += iV;
                  iN++;
              }
          }
          if( iN > 0. )
          {
              iSmoothedEvalGauss->Fill( e_log10, iM / iN );
          }

     }
     iSmoothedEval->Draw( "same" );
     iSmoothedEvalGauss->Draw( "same" ); */

     //////////////////////////////////
     // Gaussian smoothing
     TGraph *iSmoothedGaussian = new TGraph( 1 );
     iSmoothedGaussian->SetLineColor( 30 );
     iSmoothedGaussian->SetMarkerColor( 30 );
     iSmoothedGaussian->SetMarkerStyle( 22 );
     z = 0;
     for( unsigned int i = 0; i < 1000; i++ )
     {
          double e_log10 = e_min + i * (e_max - e_min)/1000.;
          double iM = 0.;
          double iN = 0.;
          for( unsigned int j = 0; j < 1000; j++ )
          {
              e_log10_G = e_log10 + gRandom->Gaus( 0., dEres );
              // make sure that values are in the given energy range
              // (extrapolation might otherwise result in very wrong
              // values)
              if( e_log10_G > e_min )
              {
                  iV = iG->Eval( e_log10_G );
                  // accept only cut values below a maximum
                  // (optimisation might give in the threshold region
                  // very high values)
                  if( iV < acceptMax )
                  {
                      iM += iV;
                      iN++;
                  }
              }
          }
          if( iN > 0. )
          {
             if( e_log10 > energy_constant_mva )
             {
                 iSmoothedGaussian->SetPoint( z, e_log10,
                          iSmoothedGaussian->Eval( e_min + (i-1) * (e_max - e_min)/1000. ) );
             }
             else
             {
                 iSmoothedGaussian->SetPoint( z, e_log10, iM / iN );
             }
             z++;
          }
     }
     iSmoothedGaussian->Draw( "pl" );
     cout << "Green color: smoothed Gaussian with " << dEres << " width" << endl;

     // moving average
     TGraph *iMovingAverage = new TGraph( 1 );
     iMovingAverage->SetMarkerStyle( 21 );
     iMovingAverage->SetMarkerColor( 4 );
     z = 0;
     for( int i = 0; i < iG->GetN(); i++ )
     {
          double iM = 0.;
          double iN = 0.;
          for( int j = 0; j < iG->GetN(); j++ )
          {
              if( j >= i - iMVAiter && j <= i + iMVAiter )
              {
                    iG->GetPoint( j, x, y );
                    double w = 1.;
                    if( TMath::Abs( j - i ) + 1 > 0. )
                    {
                         w = 1./(TMath::Abs( j - i )+1);
                    }
                    iM += y * w;
                    iN += w;
              }
          }
          if( iN > 0. )
          {
              iG->GetPoint( i, x, y );
              iMovingAverage->SetPoint( z, x, iM / iN );
              z++;
          }
      }
      iMovingAverage->Draw( "pl" );
      cout << "Blue color: moving average" << endl;
}
