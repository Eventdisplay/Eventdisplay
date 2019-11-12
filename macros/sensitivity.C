/* \file plot_sensitivity.C
   \brief example script of how to plot integral or differential sensitivities

   do help() for help...

   plot sensitivities from Crab Nebula data or simulations(CTA or VERITAS)

   Author: Gernot Maier, Heike Prokoph

*/

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void getPlottingData( bool bIntegral, string bUnit, string iObservatory );
void plotSensitivity( char* iData_anasumFile1, char* iData_anasumFile2, bool bIntegral,
                      char* iMC_Gamma, char* iMC_Proton, char* iMC_Helium, char* iMC_Electron,
                      string iFluxUnit, unsigned int iCrabSpec_ID, string iObservatory, double iObservingTime_h );

void plotDebugComparisionPlots( string iFile, int iColor );

/*

*/
void help()
{
    cout << endl;
    cout << "plot sensitivity from Crab Nebula data (VTS) or MC effective areas files (VTS and CTA)" << endl;
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "plot differential sensitivity: " << endl;
    cout << "------------------------------"  << endl;
    cout << endl;
    cout << "plotDifferentialSensitivity( string fFluxUnit, char *ifile1 = 0, char *ifile2 = 0, char *iMC_Gamma = 0, char *iMC_Proton = 0, ... ) " << endl;
    cout << endl;
    cout << "flux units: " << endl;
    cout << "\t PFLUX \t Flux Sensitivity [cm^{-2} s^{-1} TeV^{-1}]" << endl;
    cout << "\t ENERG \t E^{2} x Flux Sensitivity [erg cm^{-2} s^{-1}]" << endl;
    cout << "\t CU \t Flux Sensitivity in Crab Units [C.U.]" << endl;
}

/*

   plotting data (axes values, sensitivity text files, etc)

*/
struct cSensitivityPlottingData
{
    bool bSet;
    bool bIntegral;
    
    string fESpecDataFile_CrabNebula;
    string fESpecDataFile_CosmicRays;
    double fPlotting_flux_min;             // note: can have different units
    double fPlotting_flux_max;
    double fPlotting_energy_min_TeV;
    double fPlotting_energy_max_TeV;
    vector< string > fSensitivityvsEnergyFromTextTFile;
    vector< string > fSensitivityvsEnergyFromTextTFile_LegendTitles;
    vector< int >    fSensitivityvsEnergyFromTextTFile_LineColor;
    vector< int >    fSensitivityvsEnergyFromTextTFile_LineStyle;
    
};

cSensitivityPlottingData fPD;

/*!

     get plotting parameters
     (axis dimensions, etc)

*/
void getPlottingData( bool bIntegral, string bUnit, string iObservatory )
{
    fPD.bSet = false;
    fPD.fSensitivityvsEnergyFromTextTFile.clear();
    fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.clear();
    fPD.fSensitivityvsEnergyFromTextTFile_LineColor.clear();
    fPD.fSensitivityvsEnergyFromTextTFile_LineStyle.clear();
    
    // Crab Nebula and Cosmic ray fluxes
    fPD.fESpecDataFile_CrabNebula = "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat";
    fPD.fESpecDataFile_CosmicRays = "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CR.dat";
    
    // energy axis (min and max in TeV)
    fPD.fPlotting_energy_min_TeV = 0.01;
    fPD.fPlotting_energy_max_TeV = 250.00;
    if( iObservatory == "VTS" || iObservatory == "VERITAS" )
    {
        fPD.fPlotting_energy_min_TeV = 0.06;
        fPD.fPlotting_energy_max_TeV = 15.00;
    }
    
    //////////////////
    // integral flux
    cout << "FILLING " << bIntegral << "\t" << bUnit << endl;
    if( bIntegral )
    {
        if( bUnit == "PFLUX" )
        {
            fPD.bIntegral = true;
            fPD.fPlotting_flux_min = 8.e-16;
            fPD.fPlotting_flux_max = 2.e-10;
            fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/HESS_IntegralSensitivity_Achieved_Moriond2009.txt" );
            fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "HESS" );
            fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/MAGIC_IntegralSensitivity_Data_Colin2010_PFLUX.txt" );
            fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "MAGIC" );
            fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/MAGICII_IntegralSensitivity_MC_Colin2010_PFLUX.txt" );
            fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "MAGIC II" );
            
            fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/GLAST5Y_IntegralSensitivity_LOI_2009_PFLUX.txt" );
            fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "Fermi LAT 5y" );
            
            fPD.bSet = true;
        }
        else if( bUnit == "CU" )
        {
            fPD.bIntegral = true;
            fPD.fPlotting_flux_min = 0.0002;
            fPD.fPlotting_flux_max = 0.85;
            
            fPD.bSet = true;
        }
        else if( bUnit == "ENERGY" )
        {
            fPD.bIntegral = true;
            fPD.fPlotting_flux_min = 1.e-14;
            fPD.fPlotting_flux_max = 1.e-11;
            
            fPD.bSet = true;
        }
    }
    /////////////////////////////////////////////////
    // differential flux
    else
    {
        if( bUnit == "PFLUX" )
        {
            fPD.bIntegral = false;
            fPD.fPlotting_flux_min = 2.e-15;
            fPD.fPlotting_flux_max = 8.e-09;
            if( iObservatory == "VTS" || iObservatory == "VERITAS" )
            {
                fPD.fPlotting_flux_min = 1.e-15;
                fPD.fPlotting_flux_max = 2.e-08;
                fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$VERITAS_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/MAGIC_DifferentialSensitivity_020dE_1108.1477_PFLUX.txt" );
                fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "MAGIC stereo" );
            }
            else
            {
                fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/CTA_Typical_DifferentialSensitivity_020dE_LOI_2009_PFLUX.txt" );
                fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "CTA_LOI" );
            }
            
            fPD.bSet = true;
        }
        else if( bUnit == "CU" )
        {
            if( iObservatory == "VTS" || iObservatory == "VERITAS" )
            {
                fPD.fPlotting_flux_min = 0.002;
                fPD.fPlotting_flux_max = 0.85;
                fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$VERITAS_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/MAGIC_DifferentialSensitivity_020dE_1108.1477_CU.txt" );
                fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "MAGIC stereo" );
            }
            else
            {
                fPD.fPlotting_flux_min = 0.001;
                fPD.fPlotting_flux_max = 1.;
                fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/CTA_Typical_DifferentialSensitivity_020dE_Zurich2009_CU.txt" );
                fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "KB_Zurich2009" );
            }
            
            fPD.bSet = true;
        }
        else if( bUnit == "ENERGY" )
        {
            fPD.fSensitivityvsEnergyFromTextTFile.push_back( "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/sensitivity/CTA_Typical_DifferentialSensitivity_020dE_Zurich2009_CU.txt" );
            fPD.fSensitivityvsEnergyFromTextTFile_LegendTitles.push_back( "KB_Zurich2009" );
            fPD.fPlotting_flux_min = 1.e-14;
            fPD.fPlotting_flux_max = 2.e-10;
            
            if( iObservatory == "VTS" || iObservatory == "VERITAS" )
            {
                fPD.fPlotting_flux_min = 2.e-13;
                fPD.fPlotting_flux_max = 2.e-10;
            }
            
            fPD.bSet = true;
        }
    }
    
    // marker colors and line styles
    for( unsigned int i = 0; i < fPD.fSensitivityvsEnergyFromTextTFile.size(); i++ )
    {
        fPD.fSensitivityvsEnergyFromTextTFile_LineColor.push_back( 6 + i );
        fPD.fSensitivityvsEnergyFromTextTFile_LineStyle.push_back( 6 );
    }
}

////////////////////////////////////////////////////////////////////////////////////
/*

   plot integral sensitivity from measured data (using Crab spectra) and MC

*/
////////////////////////////////////////////////////////////////////////////////////
void plotIntegralSensitivity( string iFluxUnit = "PFLUX",
                              char* iDatafile1 = 0, char* iDatafile2 = 0,
                              char* iMC_Gamma = 0, char* iMC_Proton = 0, char* iMC_Helium = 0, char* iMC_Electron = 0,
                              unsigned int iCrabSpec_ID = 5, string iObservatory = "CTA", double iObservingTime_h = 50. )
{
    plotSensitivity( iDatafile1, iDatafile2, true,
                     iMC_Gamma, iMC_Proton, iMC_Helium, iMC_Electron, iFluxUnit, iCrabSpec_ID, iObservatory, iObservingTime_h );
}

////////////////////////////////////////////////////////////////////////////////////
/*

   plot differential sensitivity from measured data (using Crab spectra) and MC


*/
////////////////////////////////////////////////////////////////////////////////////
void plotDifferentialSensitivity( string iFluxUnit = "ENERGY",
                                  char* iDatafile1 = 0, char* iDatafile2 = 0,
                                  char* iMC_Gamma = 0, char* iMC_Proton = 0, char* iMC_Helium = 0, char* iMC_Electron = 0,
                                  unsigned int iCrabSpec_ID = 5, string iObservatory = "CTA",
                                  double iObservingTime_h = 50. )
{
    plotSensitivity( iDatafile1, iDatafile2, false,
                     iMC_Gamma, iMC_Proton, iMC_Helium, iMC_Electron, iFluxUnit, iCrabSpec_ID, iObservatory, iObservingTime_h );
}


////////////////////////////////////////////////////////////////////////////////////
/*

   plot sensitivity

   this function should not be called by the user. Use plotIntegralSensitivity() or plotDifferentialSensitivity()

   iData_anasumFile1        anasum file from Crab Nebula analysis
   iData_anasumFile2        anasum file from Crab Nebula analysis
   iMC_Gamma                effective area file for primary gamma rays
   iMC_Proton               effective area file for primary protons
   iMC_Helium               effective area file for primary helium
   iMC_Electron             effective area file for primary electron
   iFluxUnit                flux units (see help() )
   iCrabSpec_ID             Crab Nebula spectrum ID (read from text file with spectral parameters)

*/
////////////////////////////////////////////////////////////////////////////////////
void plotSensitivity( char* iData_anasumFile1, char* iData_anasumFile2, bool bIntegral,
                      char* iMC_Gamma, char* iMC_Proton, char* iMC_Helium, char* iMC_Electron,
                      string iFluxUnit, unsigned int iCrabSpec_ID, string iObservatory,
                      double iObservingTime_h )
{

    // get values for plotting
    getPlottingData( bIntegral, iFluxUnit, iObservatory );
    if( !fPD.bSet )
    {
        cout << "plotSensitivity error: no plotting datat set for " << iFluxUnit << "\t" << bIntegral << endl;
        return;
    }
    
    // sensitivity plotter
    VSensitivityCalculator a;
    a.setEnergyRange_Lin( fPD.fPlotting_energy_min_TeV, fPD.fPlotting_energy_max_TeV );           // x-axis range: in TeV
    // Monte Carlo only
    if( iFluxUnit == "PFLUX" )
    {
        a.setFluxRange_PFLUX( fPD.fPlotting_flux_min, fPD.fPlotting_flux_max );           // y-axis range: in 1/cm2/s
    }
    else if( iFluxUnit == "ENERGY" )
    {
        a.setFluxRange_ENERG( fPD.fPlotting_flux_min, fPD.fPlotting_flux_max );
    }
    else
    {
        a.setFluxRange_CU( fPD.fPlotting_flux_min, fPD.fPlotting_flux_max );
    }
    a.setPlotCanvasSize( 900, 600 );                                                              // size of canvases
    
    // set Crab Nebula spectrum used for relative flux calculation
    if( !a.setEnergySpectrumfromLiterature( fPD.fESpecDataFile_CrabNebula, iCrabSpec_ID ) )
    {
        return;
    }
    
    // lots of debug output
    a.setDebug( false );
    
    a.setSignificanceParameter( 5., 10., iObservingTime_h, 0.05, 0.2 );
    
    //////////////////////////////////////////////////////////////////
    // plot sensitivity canvas
    TCanvas* c = a.plotCanvas_SensitivityvsEnergy( iFluxUnit, bIntegral );
    if( !c )
    {
        return;
    }
    
    //////////////////////////////////////////////////////////////////
    // plot sensitivities from text files
    for( unsigned int i = 0; i < fPD.fSensitivityvsEnergyFromTextTFile.size(); i++ )
    {
        a.plotSensitivityvsEnergyFromTextTFile( c, fPD.fSensitivityvsEnergyFromTextTFile[i],
                                                fPD.fSensitivityvsEnergyFromTextTFile_LineColor[i], 3, fPD.fSensitivityvsEnergyFromTextTFile_LineStyle[i],
                                                iFluxUnit );
    }
    a.setPlottingStyle( 1, 1, 2, 20, 2., 3002 );
    
    //////////////////////////////////////////////////////////////////////
    // calculate sensitivities from data taken towards the Crab Nebula
    //////////////////////////////////////////////////////////////////////
    if( iData_anasumFile1 != 0 )
    {
        if( bIntegral )
        {
            a.plotIntegralSensitivityvsEnergyFromCrabSpectrum( c, iData_anasumFile1, 1, iFluxUnit, 0.05, 100. );
        }
        else
        {
            a.plotDifferentialSensitivityvsEnergyFromCrabSpectrum( c, iData_anasumFile1, 1, iFluxUnit, 0.2, 0.05, 8. );
        }
    }
    // plot second file
    if( iData_anasumFile2 != 0 )
    {
        a.setPlottingStyle( 3, 1, 2, 20, 2., 3002 );
        if( bIntegral )
        {
            a.plotIntegralSensitivityvsEnergyFromCrabSpectrum( c, iData_anasumFile2, 3, iFluxUnit, 0.05, 100. );
        }
        else
        {
            a.plotDifferentialSensitivityvsEnergyFromCrabSpectrum( c, iData_anasumFile2, 3, iFluxUnit, 0.2, 0.05, 100. );
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////
    // calculate sensitivities from Monte Carlo effective areas
    ////////////////////////////////////////////////////////////////////////
    if( iMC_Gamma && iMC_Proton )
    {
        VSensitivityCalculator bb;
        bb.setDebug( false );             // creates lots of output
        // set Crab Nebula spectrum
        bb.setEnergySpectrumfromLiterature( fPD.fESpecDataFile_CrabNebula, iCrabSpec_ID );
        // draw some debugging plots
        char hname[600];
        sprintf( hname, "plotDebug_%d", bIntegral );
        bb.setPlotDebug( hname );
        // energy range to be plotted
        bb.setEnergyRange_Lin( 0.01, 250. );
        // significance parameters
        bb.setSignificanceParameter( 5., 10., iObservingTime_h, 0.05, 0.2 );
        // set colors different in case data is plotted
        if( iData_anasumFile1 )
        {
            bb.setPlottingStyle( 4, 1, 2, 20, 2., 3002 );
        }
        
        //////////////////////////////////////////////////////////////////////////
        // select bins and index from gamma and proton effective area files
        // CTA
        int i_Azbin_gamma = 0;
        double i_index_gamma = 2.5;
        int i_noise_gamma = 250;
        double i_woff_gamma = 0.;
        
        int i_Azbin_proton = 0;
        double i_index_proton = 2.6;
        int i_noise_proton = 250;
        double i_woff_proton = 0.;
        
        int i_Azbin_electron = 0;
        double i_index_electron = 3.0;
        int i_noise_electron = 250;
        double i_woff_electron = 0.;
        
        if( iObservatory == "VTS" || iObservatory == "VERITAS" )
        {
            i_Azbin_gamma = 0;
            i_index_gamma = 2.4;
            i_noise_gamma = 130;
            i_noise_gamma = 200;
            i_woff_gamma = 0.5;
            
            i_Azbin_proton = 0;
            i_index_proton = 2.6;
            
            i_noise_proton = 130;
            i_noise_proton = 200;
            i_woff_proton = 0.;
            
            i_Azbin_electron = 0;
            i_index_electron = 3.0;
            
            i_noise_electron = 130;
            i_noise_electron = 200;
            
            i_woff_electron = 0.;
        }
        cout << "SETTING EFFECTIVE AREA SEARCH VALUES TO " << iObservatory << endl;
        //////////////////////////////////////////////////////////////////////////
        
        // gammas
        bb.setMonteCarloParameters( 1, fPD.fESpecDataFile_CrabNebula, iCrabSpec_ID, iMC_Gamma, 20.,
                                    i_Azbin_gamma, i_woff_gamma, i_noise_gamma, i_index_gamma );
        // protons
        bb.setMonteCarloParameters( 14, fPD.fESpecDataFile_CosmicRays, 0, iMC_Proton, 20.,
                                    i_Azbin_proton, i_woff_proton, i_noise_proton, i_index_proton );
        // helium (spectral index?)
        if( iMC_Helium )
        {
            bb.setMonteCarloParameters( 402, fPD.fESpecDataFile_CosmicRays, 1, iMC_Helium, 20., 0, 0.0, 200, 2.0 );
        }
        // electrons (spectral index?)
        if( iMC_Electron )
        {
            bb.setMonteCarloParameters( 2, fPD.fESpecDataFile_CosmicRays, 2, iMC_Electron, 20.,
                                        i_Azbin_electron, i_woff_electron, i_noise_electron, i_index_electron );
        }
        
        // energy range determined by looking at number of noff events (need off events to determine sensitivity)
        if( bIntegral )
        {
            bb.plotIntegralSensitivityvsEnergyFromCrabSpectrum( c, "MC", 1, iFluxUnit, 0.01, 500. );
        }
        else
        {
            bb.plotDifferentialSensitivityvsEnergyFromCrabSpectrum( c, "MC", 1, iFluxUnit, 0.2, 0.01 );
        }
        
        bb.plotSensitivityLimitations( c );
    }
    return;
    // plot different limitations in sensitivity calculation
    a.plotSensitivityLimitations( c );
    
}

////////////////////////////////////////////////////////////////////////////
/*

  plot sensitivity vs time for a given gamma-ray and background rate from the Crab)

*/
void plotSensitivityVStime()
{
    VSensitivityCalculator a;
    // V330, theta2 < 0.15, size > 500, mscw/mscl<0.5 (May 2008)
    a.addDataSet( 4.97, 0.12, 0.09, "VERITAS (Autumn 2009)" );
    // from Crab data set Autumn 2009. LL R4 theta2 < 0.08 (mscw<0.35 && mscl<0.7; April 2010)
    a.addDataSet( 5.19, 0.27, 0.09, "VERITAS (Spring 2009)" );
    
    vector< double > s;
    s.push_back( 0.01 );
    s.push_back( 0.05 );
    s.push_back( 0.30 );
    a.list_sensitivity( 0 );
    a.setSourceStrengthVector_CU( s );
    
    TCanvas* c = a.plotObservationTimevsFlux( 0, 0, 2 );
    //    a.plotObservationTimevsFlux( 1, c, 1, 2 );
    a.plotObservationTimevsFlux( 0, c, 2 );
}

/*

   plot some debug plots for comparision

*/
void plotDebugComparisionPlots( string iFileName, int iColor )
{
    TH1F* hGammaEffArea = 0;
    TH1F* hBGRate = 0;
    TH1F* hDiffSens = 0;
    TCanvas* c = 0;
    
    TFile* iFile = new TFile( iFileName.c_str() );
    if( iFile->IsZombie() )
    {
        return;
    }
    
    hGammaEffArea = ( TH1F* )iFile->Get( "EffectiveArea" );
    hBGRate       = ( TH1F* )iFile->Get( "BGRate" );
    hDiffSens     = ( TH1F* )iFile->Get( "DiffSens" );
    
    // effective area canvas
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "plotDebug_0_1" );
    if( c && hGammaEffArea )
    {
        c->cd();
        hGammaEffArea->SetLineWidth( 2 );
        hGammaEffArea->SetLineColor( iColor );
        hGammaEffArea->SetMarkerColor( iColor );
        hGammaEffArea->Draw( "same" );
    }
    
    // background rates
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "plotDebug_0_0" );
    if( c && hBGRate )
    {
        c->cd();
        hBGRate->Scale( 60. );
        hBGRate->SetLineWidth( 2 );
        hBGRate->SetLineColor( iColor );
        hBGRate->SetMarkerColor( iColor );
        hBGRate->Draw( "same" );
    }
    
    // sensitivity
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "iCanvas" );
    if( c && hDiffSens )
    {
        c->cd();
        hDiffSens->SetLineWidth( 2 );
        hDiffSens->SetLineColor( iColor );
        hDiffSens->SetMarkerColor( iColor );
        hDiffSens->Draw( "same" );
    }
}

/*

    write files with particle number spectra for on (gamma) and off (protons+electrons) counts

    files are needed e.g. for setting the optimal cut value for TMVA cuts

    use writeAllParticleNumberFiles() for writing the files for all sub arrays

*/
void writeParticleNumberFile( char* iMC_Gamma = 0, char* iMC_Proton = 0, char* iMC_Electron = 0,
                              unsigned int iCrabSpec_ID = 5, string iParticleNumberFile = "particleNumbers.tmp.root",
                              string iObservatory = "CTA" )
{
    string iESpecDataFile_CrabNebula = "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat";
    string iESpecDataFile_CosmicRays = "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CR.dat";
    
    if( iMC_Gamma && iMC_Proton )
    {
        VSensitivityCalculator b;
        b.setDebug( false );             // creates lots of output
        // set Crab Nebula spectrum
        b.setEnergySpectrumfromLiterature( iESpecDataFile_CrabNebula, iCrabSpec_ID );
        // draw some debugging plots
        b.setWriteParticleNumberFile( iParticleNumberFile );
        // CTA
        int i_Azbin_gamma = 0;
        double i_index_gamma = 2.5;
        int i_noise_gamma = 250;
        double i_woff_gamma = 0.;
        
        int i_Azbin_proton = 0;
        double i_index_proton = 2.6;
        int i_noise_proton = 250;
        double i_woff_proton = 0.;
        cout << "SETTING EFFECTIVE AREA SEARCH VALUES TO " << iObservatory << endl;
        //////////////////////////////////////////////////////////////////////////
        
        // gammas
        b.setMonteCarloParameters( 1, iESpecDataFile_CrabNebula, iCrabSpec_ID, iMC_Gamma, 20.,
                                   i_Azbin_gamma, i_woff_gamma, i_noise_gamma, i_index_gamma );
        // protons
        b.setMonteCarloParameters( 14, iESpecDataFile_CosmicRays, 0, iMC_Proton, 20.,
                                   i_Azbin_proton, i_woff_proton, i_noise_proton, i_index_proton );
        // electrons (spectral index?)
        if( iMC_Electron )
        {
            b.setMonteCarloParameters( 2, iESpecDataFile_CosmicRays, 2, iMC_Electron, 20., 0, 0.0, 250, 3.0 );
        }
        b.calculateParticleNumberGraphs_MC( 0.2 );
    }
}

/*

    write files with particle number spectra for on (gamma) and off (protons+electrons) counts

    files are needed e.g. for setting the optimal cut value for TMVA cuts

    (for all sub arrays)

*/
void writeAllParticleNumberFiles( char* iSubArrayFile = 0,
                                  int iOffSetCounter = 8, int iRecID = 0,
                                  char* iMC_Gamma_onSource = "gamma_onSource", char* iMC_Gamma_cone10 = "gamma_cone10",
                                  char* iMC_Proton = "proton", char* iMC_Electron = "electron" )
{
    vector< string > SubArray;
    ifstream is;
    is.open( gSystem->ExpandPathName( iSubArrayFile ), ifstream::in );
    if( !is )
    {
        cout << "error opening subarray file " << iSubArrayFile << endl;
        return;
    }
    
    string is_line;
    while( getline( is, is_line ) )
    {
        SubArray.push_back( is_line );
    }
    
    char iGamma[200];
    char iProton[200];
    char iElectron[200];
    char iParticleNumberFile[200];
    
    for( unsigned int i = 0; i < SubArray.size(); i++ )
    {
        cout << "STARTING SUBARRAY " << SubArray[i] << endl;
        
        if( iMC_Gamma_onSource )
        {
            sprintf( iParticleNumberFile, "ParticleNumbers.%s.00.root", SubArray[i].c_str() );
            sprintf( iGamma, "%s.%s_ID%d.eff-%d.root", iMC_Gamma_onSource, SubArray[i].c_str(), iRecID, 0 );
            sprintf( iProton, "%s.%s_ID%d.eff-%d.root", iMC_Proton, SubArray[i].c_str(), iRecID, 0 );
            if( iMC_Electron )
            {
                sprintf( iElectron, "%s.%s_ID%d.eff-%d.root", iMC_Electron, SubArray[i].c_str(), iRecID, 0 );
                writeParticleNumberFile( iGamma, iProton, iElectron, 6, iParticleNumberFile );
            }
            else
            {
                writeParticleNumberFile( iGamma, iProton, 0, 6, iParticleNumberFile );
            }
        }
        
        // offset files
        for( int j = 0; j < iOffSetCounter; j++ ) // use first bin on source particle file
        {
            sprintf( iGamma, "%s.%s_ID%d.eff-%d.root", iMC_Gamma_cone10, SubArray[i].c_str(), iRecID, j );
            sprintf( iProton, "%s.%s_ID%d.eff-%d.root", iMC_Proton, SubArray[i].c_str(), iRecID, j );
            sprintf( iElectron, "%s.%s_ID%d.eff-%d.root", iMC_Electron, SubArray[i].c_str(), iRecID, j );
            
            sprintf( iParticleNumberFile, "ParticleNumbers.%s.%d.root", SubArray[i].c_str(), j );
            
            writeParticleNumberFile( iGamma, iProton, iElectron, 6, iParticleNumberFile );
        }
    }
}

/*

   plot all IRFs and sensitivity in comparision

*/

void plotIRF( string iSubArray, bool iSensitivity, string iObservingTime_H, string iWPPhysSensitivityDirectory )
{
    string iGammaEffArea    = "gamma_onSource." + iSubArray + "_ID0.eff-0.root";
    string iProtonEffArea   = "proton." + iSubArray + "_ID0.eff-0.root";
    string iElectronEffArea = "electron." + iSubArray + "_ID0.eff-0.root";
    
    string iIFAE = iWPPhysSensitivityDirectory + "IFAEPerformanceBCDEINANB_Nov2011/Subarray";
    iIFAE += iSubArray + "_IFAE_" + iObservingTime_H + "hours_20111109.root";
    string iISDC = iWPPhysSensitivityDirectory + "ISDC/ISDC_2000m_KonradB_optimal_";
    iISDC += iSubArray;
    if( iObservingTime_H == "50" )
    {
        iISDC += "_50.0";
    }
    iISDC += "h_20deg_20110615.root";
    
    VPlotInstrumentResponseFunction a;
    
    a.setPlottingAxis( "energy_Lin", "X", true, 0.01, 200 );
    a.setPlottingAxis( "effarea_Lin", "X", true, 50., 5.e7 );
    a.setPlottingAxis( "energyresolution_Lin", "X", false, 0., 0.7 );
    
    a.addInstrumentResponseData( iGammaEffArea );
    a.addInstrumentResponseData( iIFAE );
    a.addInstrumentResponseData( iISDC );
    
    a.plotEffectiveArea();
    a.plotAngularResolution();
    a.plotEnergyResolution( 0.7 );
    
    if( iSensitivity )
    {
        //        plotDifferentialSensitivity( "ENERGY", 0, 0, iGammaEffArea, iProtonEffArea, 0, iElectronEffArea );
        plotDebugComparisionPlots( iIFAE, 2 );
        plotDebugComparisionPlots( iISDC, 3 );
    }
}

