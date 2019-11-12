//! VSensitivityCalculator

#ifndef VSensitivityCalculator_H
#define VSensitivityCalculator_H

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "CEffArea.h"
#include "VGammaHadronCuts.h"
#include "VEnergySpectrum.h"
#include "VEnergySpectrumfromLiterature.h"
#include "VHistogramUtilities.h"
#include "VMonteCarloRateCalculator.h"
#include "VMonteCarloRunHeader.h"
#include "VStatistics.h"
#include "VTMVAEvaluator.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct sSensitivityData
{
    string fName;
    double fSignal;
    double fBackground;
    double fAlpha;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// data class for response functions for a give primary
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class VSensitivityCalculatorDataResponseFunctions
{
    public:
    
        string fName;
        unsigned int fParticleID;                              //              particle ID (1=gamma,2=electron,14=proton,402=helium)
        string fSpectralParameterFile;                         //              data file (+path) with spectral information
        //              (e.g. Crab Nebula or proton spectrum)
        unsigned int fSpectralParameterID;                     //              ID of spectrum in given file fSpectralParameterFile
        string fEffectiveAreaFile;                             //              file name (+path) of effective area file
        double ze;                                             // [deg]        zenith angle of selected effective area
        int az;                                                //              azimuth bin of selected effective area
        double woff;                                           // [deg]        wobble offset of selected effective area
        int noise;                                             // [GrISU]      noise level of selected effective area
        double index;                                          //              spectral index of selected effective area
        double theta2_min;                                     // [deg^2]      direction cut (if energy independent)
        double theta2_MCScatterAngle;                          // [deg^2]      scattering angle^2 of primary direction in CORSIKA
        //              (e.g. 10deg in most CTA simulations)
        TGraph* gSolidAngle_DirectionCut_vs_EnergylgTeV;       // [sr, lg TeV] solid angle of direction cut (as function of energy)
        TGraph* gTheta2Cuts_vsEnergylgTeV;                     // [deg, TeV]   theta2 cut vs energy
        double SolidAngle_MCScatterAngle;                      // [sr]         scattering solid angle of primary direction in CORSIKA
        //              (e.g. 10deg in most CTA simulations)
        double alpha;
        int effArea_Ebins;
        double effArea_Emin;
        double effArea_Emax;
        double energy_min_log;
        double energy_max_log;
        vector< double > energy;
        vector< double > energy_lowEdge;
        vector< double > energy_upEdge;
        vector< double > effArea;
        vector< double > effArea_error;
        TH2D* hResponseMatrix;
        TH1D* hWeightedRate;
        
        VSensitivityCalculatorDataResponseFunctions();
        ~VSensitivityCalculatorDataResponseFunctions();
        
        double getSolidAngle_DirectionCut( double e );         // get solid angle for direction cut ( energy in lg10 [TeV] )
        void   initializeHistograms( string );
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class VSensitivityCalculator : public TObject, public VPlotUtilities, public VHistogramUtilities
{
    private:
    
        double fConstant_Flux_To_Ergs;
        
        bool fDebug;
        
        // source strength vector (to plot and to list)
        vector< double > fSourceStrength;
        double fSourceStrength_min;
        double fSourceStrength_max;
        double fSourceStrength_step;
        
        // energy spectrum from literature
        VEnergySpectrumfromLiterature* fEnergySpectrumfromLiterature;
        unsigned int fEnergySpectrumfromLiterature_ID;
        vector< TGraph* > fCrabFlux_SourceStrength;
        
        // TMVA cut optimization results
        VTMVAEvaluatorResults* fTMVAEvaluatorResults;
        bool   fRequireCutsToBeOptimized;
        
        // observing time
        double fObservationTime_h;                  // [h]
        double fObservationTime_min;
        double fObservationTime_max;
        int    fObservationTime_steps;
        
        // significance calculation
        unsigned int fLiAndMaEqu;
        double       fSignificance_min;
        double       fEvents_min;
        double       fBackgroundEvents_min;
        double       fMinBackgroundRateRatio_min;
        double       fAlpha;
        bool         fSetEvents_minCutOnly;
        
        map< int, double > fSignificanceLimited;             // sensitivity is significance limited at this energy (linear) [GeV] !!!
        map< int, double > fMinEventsLimited;                // sensitivity is limited by minimum number of events at this energy [GeV] !!!
        map< int, double > fMinBackgroundEventsLimited;      // sensitivity is limited by minimum number of background events at this energy [GeV] !!!
        map< int, double > fMinNoBackground;                 // no background events
        
        double fEnergy_min_Log;
        double fEnergy_max_Log;
        double fEnergy_dE_log10;
        double fMaxEnergyBias;
        
        // data vectors for MC and int/diff sensitivity calculation
        map< unsigned int, VSensitivityCalculatorDataResponseFunctions* > fMC_Data;                        //! [particle ID]
        map< unsigned int, VSensitivityCalculatorDataResponseFunctions* >::iterator fMC_Data_iterator;     //! [particle ID]
        VSensitivityCalculatorDataResponseFunctions fMC_GammaData;
        VSensitivityCalculatorDataResponseFunctions fMC_ProtonData;
        double fMC_BackgroundMissingParticleFraction;     //! this is the fraction of background missing due to missing He/.. simulations
        
        string fMCCTA_File;
        double fMCCTA_cameraoffset_deg;
        
        unsigned int    fCurrentDataSet;
        vector< sSensitivityData > fData;
        map< unsigned int, TGraph* > fGraphObsvsTime;
        
        // plotting values
        TH1D*  hnull;
        string fPlot_CanvasName;
        string fPlot_CanvasTitle;
        int    fPlot_CanvasSize_x;
        int    fPlot_CanvasSize_y;
        double fPlot_flux_PFLUX_min;
        double fPlot_flux_PFLUX_max;
        double fPlot_flux_ENERG_min;
        double fPlot_flux_ENERG_max;
        double fPlot_flux_CU_min;
        double fPlot_flux_CU_max;
        
        // sensitivity and rate graph
        TGraphAsymmErrors* gSensitivityvsEnergy;
        TGraphAsymmErrors* gSignalRate;
        TGraphAsymmErrors* gBGRate;
        TGraphAsymmErrors* gProtonRate;
        TGraphAsymmErrors* gElectronRate;
        
        // plotting debug stuff
        vector< TCanvas* > cPlotDebug;
        string fPlotDebugName;
        string fDebugParticleNumberFile;                // write non/noff to disk
        // values of Crab fluxes to be plotted as lines into the sensitivity vs energy graph (in Crab Units)
        vector< double > fPlottingCrabFlux_CU;
        
        // private functions
        void       applyEnergyBiasCut( TH1F*, TGraphAsymmErrors*, double iMaxEBias );
        bool       checkCutOptimization( double iEnergy, bool iPrint = false );
        bool       checkDataSet( unsigned int iD, string iName );
        bool       checkUnits( string iUnit );
        bool       fillSensitivityHistogramfromGraph( TGraph* g, TH1F* h, double iScale = 1.,
                bool bFillErrors = true, bool bSQRT = false );
        bool       fillSensitivityHistogramfromMap( map< int, double > m, TH1F* h );
        TGraph*           getCrabSpectrum( bool bIntegralSpectrum,  string bUnit = "CU", bool bReset = true, bool iUseMAGICCurvedCrabSpectrum = false );
        vector< TGraph* > getCrabSpectrum( vector< double > i_fCrabFlux, bool bIntegralSpectrum,  string bUnit = "CU", bool bReset = true, bool iUseMAGICCurvedCrabSpectrum = false );
        
        vector< VDifferentialFluxData > getDifferentFluxVectorfromData( string iAnasumCrabFile, double dE_Log10, double& iNorm );
        vector< VDifferentialFluxData > getDifferentialFluxVectorfromMC( double dE_Log10, double& iNorm );
        vector< VDifferentialFluxData > getDifferentialFluxVectorfromCTA_MC( double& iNorm );
        vector< VDifferentialFluxData > getDifferentialFluxVectorfromMC_ErrorMessage( string );
        
        void       fillBackgroundParticleNumbers( vector< VDifferentialFluxData > iDifferentialFlux,
                map< unsigned int, vector< double > > i_flux_NOff,
                map< unsigned int, vector< double > > i_flux_NOff_error );
        bool       fillBackroundvsSquareDegree( TGraphAsymmErrors* i_R, TH1F* i_H );
        void       fillParticleNumbersGraphs( vector< VDifferentialFluxData > iDifferentialFlux, double alpha, double dE_Log10 = 1. );
        bool       getMonteCarlo_EffectiveArea( VSensitivityCalculatorDataResponseFunctions* iMCPara );
        double     getMonteCarloRateFromWeightedRateHistogram( double iE_low, double iE_up, bool iRateError, TH1D* iWeightedRateHistogram );
        
        
        TGraphAsymmErrors* getSensitivityGraphFromWPPhysFile( string bUnit,
                                                              double iEnergyMin_TeV_lin, 
                                                              double iEnergyMax_TeV_lin, 
                                                              double dE_Log10 );
        void       plot_guidingLines( double x, TGraph* g, bool iHours );
        void       plotEffectiveArea();
        void       plotDebugPlotsBackgroundParticleNumbers( vector< VDifferentialFluxData > iDifferentialFlux,
                map< unsigned int, vector< double > > i_flux_NOff,
                map< unsigned int, vector< double > > i_flux_NOff_error );
        void       plotDebugPlotsParticleNumbers();
        void       prepareDebugPlots();
        void       smoothenSensitivityGraph( TGraphAsymmErrors*, TH1F* );
        
    public:
    
        VSensitivityCalculator();
        ~VSensitivityCalculator() {}
        
        unsigned int  addDataSet( double iGammaRayRate, double iBackGroundRate, double iAlpha, string iName = "" );
        unsigned int  addDataSetOnOff( double iOnRate, double iOffRate, double iAlpha, string iName = "" );
        double   calculateObservationTimevsFlux( unsigned int iD );
        bool     calculateParticleNumberGraphs_MC( double dE_Log10 );
        bool     fillSensitivityHistograms( TH1F* iSensitivity = 0, TH1F* iBGRate = 0, TH1F* iBGRateSqDeg = 0,
                                            TH1F* iProtonRate = 0, TH1F* iProtonRateSqDeg = 0,
                                            TH1F* iElectronRate = 0, TH1F* iElectronRateSqDeg = 0,
                                            bool iHE = false,
                                            TH1F* iTheta2 = 0, TH1F* iTheta = 0 );
        bool     fillSensitivityLimitsHistograms( vector<TH1F*>& h );
        bool     getDebug()
        {
            return fDebug;
        }
        double   getSensitivity( unsigned int iD, double energy = -1., unsigned int iFillStatistics = 0 );
        double   getSensitivity( double iSignal, double iBackground, double iAlpha, double energy = -1., unsigned int iFillStatistics = 0 );
        TGraphAsymmErrors*  getSensitivityGraph()
        {
            return gSensitivityvsEnergy;
        }
        TGraph*  getObservationTimevsFluxGraph( unsigned int ID = 0 );
        unsigned int  listDataSets();
        void     listUnits();
        void     list_sensitivity( unsigned int iD = 0 );
        TCanvas* plotCanvas_SensitivityvsEnergy( string bUnit, bool bIntegralSensitivity, bool bNewCanvas = false );
        TCanvas* plotObservationTimevsFlux( unsigned int iD = 0, TCanvas* c = 0, int iLineColor = 1, double iLineWidth = 4., bool bGuidingLines = true );
        void     plotObservationTimevsFluxFromTextFile( TCanvas* c, string iTextFile, int iLineColor = 4, double iLineWidth = 1., int iLineStyle = 2 );
        TCanvas* plotSignalBackgroundRates( TCanvas* c = 0, bool bPlotParticleBackground = false, double iRateMinimum = 2.e-7, double iRateMaximum = 1.e-1 );
        void     plotSensitivityLimitations( TCanvas* cSensitivity = 0, double iY = -9999. );
        TCanvas* plotSensitivityvsEnergyFromTextTFile( TCanvas* c, string iTextFile,
                int iColor = 1, double iLineWidth = 2, int iLineStyle = 2,
                string bUnit = "CU", bool bIntegralSensitivity = true,
                double iMultFactor = 1 );
        bool     calculateSensitivityvsEnergyFromCrabSpectrum( string iAnasumCrabFile, string bUnit = "CU",
                double dE_Log10 = 0.20,
                double iEnergyMin = 0.01, double iEnergyMax = 1.e6 );
        TCanvas* plotDifferentialSensitivityvsEnergyFromCrabSpectrum( TCanvas* c, string iAnasumCrabFile,
                string bUnit = "CU",
                double dE_Log10 = 0.20,
                double iEnergyMin_TeV_lin = 0.01, double iEnergyMax_TeV_lin = 1.e6 );
        TCanvas* plotIntegralSensitivityvsEnergyFromCrabSpectrum( TCanvas* c, string iAnasumCrabFile,
                string bUnit = "CU",
                double iEnergyMin_TeV_lin = 0.01, double iEnergyMax_TeV_lin = 1.e6 );
        TCanvas*   plotSensitivityvsEnergyFromCrabSpectrum( TCanvas* c, string bUnit = "CU", double dE_Log10 = 0.25, bool bNewCanvas = false );
        void     plotSignificanceParameters( TCanvas* cSensitivity = 0 );
        bool     printSensitivity();
        bool     removeDataSet( unsigned int iD );
        void     reset();
        void     setBackgroundMissingParticleFraction( double iB = 0. )
        {
            fMC_BackgroundMissingParticleFraction = iB;
        }
        void     setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        void     setPlotCanvasSize( int x = 600, int y = 600 )
        {
            fPlot_CanvasSize_x = x;
            fPlot_CanvasSize_y = y;
        }
        void     setPlotCanvasName( string iN = "cSensitivity", string iT = "sensitivity vs energy" )
        {
            fPlot_CanvasName = iN;
            fPlot_CanvasTitle = iT;
        }
        void     setPlotCrabFluxLineValues( vector< double > iF )
        {
            fPlottingCrabFlux_CU = iF;
        }
        void     setPlotDebug( string iName )
        {
            fPlotDebugName = iName;
        }
        bool     setCurrentDataSet( unsigned int iD );
        void     setEnergyRange_Log( double iLMin = -1.9, double iLMax = 3.1 );
        void     setEnergyRange_Lin( double iLMin = TMath::Power( 10., -1.9 ), double iLMax = TMath::Power( 10., 3.1 ) );
        bool     setEnergySpectrumfromLiterature( string iFile, unsigned int iID = 1 );
        void     setEnergySpectrumfromLiterature_ID( unsigned int iID )
        {
            fEnergySpectrumfromLiterature_ID = iID;
        }
        void     setFluxRange_PFLUX( double iMin = 1.e-15, double iMax = 5.e-11 )
        {
            fPlot_flux_PFLUX_min = iMin;
            fPlot_flux_PFLUX_max = iMax;
        }
        void     setFluxRange_ENERG( double iMin = 4.e-14, double iMax = 2.5e-10 )
        {
            fPlot_flux_ENERG_min = iMin;
            fPlot_flux_ENERG_max = iMax;
        }
        void     setFluxRange_CU( double iMin = 1.e-4, double iMax = 10. )
        {
            fPlot_flux_CU_min = iMin;
            fPlot_flux_CU_max = iMax;
        }
        void     setMaximumEnergyBias( double iEmax_bias = 1.e5 )
        {
            fMaxEnergyBias = iEmax_bias;
        }
        bool     setMonteCarloParametersCTA_MC( string iCTA_MCFile, double iMCCTA_cameraoffset_deg, string iSpectralParameterFile, unsigned int iSpectralParameterID );
        void     setMonteCarloParameters( unsigned int iParticleID, string iSpectralParameterFile, unsigned int iSpectralParameterID,
                                          string iGammaEffectiveAreaFile, double ze = 20.,
                                          int az = 0, double woff = 0.5, int noise = 150, double index = 2.5,
                                          double iEnergy_min_log = -10., double iEnergy_max_log = 10., string bUnit = "CU" );
        void     setObservationTimeRange( double iObs_min = 0.5e-3, double iObs_max = 5.e4, int iObs_steps = 1000 );    // hours
        void     setRequireCutsToBeOptimized( bool iB = true )
        {
            fRequireCutsToBeOptimized = iB;
        }
        void     setSensitivityGraph( TGraphAsymmErrors* g )
        {
            gSensitivityvsEnergy = g;
        }
        void     setSignificanceParameter( double iSignificance = 5., double iMinEvents = 10., double iObservationTime = 50.,
                                           double iMinBackgroundRateRatio = 0.05, double alpha = 0.2 );
        void     setSourceStrengthRange_CU( double iMin = 0.01, double iMax = 1.5, double iStep = 0.005, bool iLog = false );
        void     setSourceStrengthVector_CU();
        void     setSourceStrengthVector_CU( vector< double > );
        void     setWriteParticleNumberFile( string iFile )
        {
            fDebugParticleNumberFile = iFile;
        }
        
        ClassDef( VSensitivityCalculator, 27 );
};
#endif
