//! VEnergySpectrum analyse and plot energy spectra

#ifndef VEnergySpectrum_h
#define VEnergySpectrum_h

#include "VAnalysisUtilities.h"
#include "VDifferentialFluxData.h"
#include "VEnergySpectrumfromLiterature.h"
#include "VEnergyThreshold.h"
#include "VHistogramUtilities.h"
#include "VMathsandFunctions.h"
#include "VPlotUtilities.h"
#include "VSpectralFitter.h"
#include "VStatistics.h"

#include "TArrow.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRolke.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VEnergySpectrum : public VAnalysisUtilities, public VPlotUtilities
{
    private:
    
        bool  bZombie;                            //! status of files
        bool  bCombineRuns;                       //! rerun combine runs
        bool  bAsciiDataFile;
        
        int  fTotalRun;                           // total number of runs found and used
        
        string fDataSetName;                      //! name for histogram, canvases, etc. (no spaces, special characters)
        
        // vector with differential flux values and upper limits
        vector< vector< VDifferentialFluxData > > fDifferentialFlux_perRun;
        vector< VDifferentialFluxData > fDifferentialFlux;
        
        double fAnalysisMinEnergy;                // linear energy axis [TeV]
        double fAnalysisMaxEnergy;                // linear energy axis [TeV]
        double fAnalysisEnergyBinning;            // log energy axis [log10(TeV)]
        bool   fAnalysisHistogramAddingUseLowEdge;
        
        unsigned int fAnalysisEnergyThresholdDefinition; // 0 = no energy threshold, 1 = XX % systematic in energy reconstruction,
        // 2 = fraction of maximum effective area, 3 = user defined
        double fAnalysisMaxEnergySystematic;      // maximum allowed systematic error in energy reconstruction
        double fAnalysisMaxEffectiveAreaFraction; // energy threshold is defined as fraction of maximum effective area
        double fAnalysisEnergyThreshold_min_absolute_TeV;
        
        double fAnalysisSignificance;
        double fAnalysisMinEvents;
        double fAnalysisUpperLimits;
        int    fAnalysisLiAndMaEquation;
        int    fAnalysisUpperLimitAlgorithm;
        string fErrorCalculationMethod;
        bool   fAnalysisForFit;
        
        // spectral fit variables
        VSpectralFitter* fSpectralFitter;
        int    fSpectralFitFunction;              // 0: power law
        // Default: 1 TeV
        double fSpectralFitFluxNormalisationEnergy;
        double fSpectralFitEnergy_min;
        double fSpectralFitEnergy_max;
        
        double fOffsetDistance;                  // offset distance to get correct counting histograms
        
        // plotting variables
        TCanvas* fPlottingCanvas;
        // spectral weight to calculate bin centers
        double fPlottingSpectralWeightForBinCenter;
        double fPlottingYaxisMin;
        double fPlottingYaxisMax;
        double fPlottingMultiplierIndex;
        double fPlottingMinEnergy;                // linear energy axis [TeV]
        double fPlottingMaxEnergy;                // linear energy axis [TeV]
        bool   fPlottingLogEnergyAxis;            // plot log or lin values on energy axis (default=true)
        bool   fPlottingUpperLimits;              // plot upper limits
        TGraphAsymmErrors* gEnergySpectrum;
        TGraphAsymmErrors* gEnergySpectrumFitResiduals;
        
        // rebinning of energy spectra
        TH1D* nRebinner;
        TH1D* nOriginalBinner;
        
        bool   bUseRebinner;
        vector< int > newBinningGroupings;
        
        int    fPlottingEnergySpectrumFitLineColor;
        int    fPlottingEnergySpectrumFitLineStyle;
        float  fPlottingEnergySpectrumFitLineWidth;
        
        TF1*   fEnergySpectrumFit;
        
        void   addValueToHistogram( TH1* h, double iTObs, double iEThreshold_log10TeV );
        void   addValueToHistogram( TH1* h, TGraph* g, double iTObs, double iEThreshold_log10TeV );
        void   addHistogram( TH1* h1, TH1* h2, double iEThreshold );
        int    getRebinningGrouping( TH1* h, double iNewBinWidth );
        bool   openAsciiFile( string iFile );
        void   rebinEnergySpectrum( TH1D* h, double iER );
        TH1*   setVariableBinning( TH1* a );
        void   setOriginalBinner( TH1* a );
        
        // histograms
        TH1D* hErecCountsOn;
        TH1D* hErecCountsOff;
        TH1D* hErecTotalTime;
        TH1D* hErecTotalTimeDeadTimeCorrected;
        TH1D* hEffArea;
        
        // definition of energy value in a log bin
        unsigned int fEnergyInBinDefinition;    // 0 = mean of (lin) energy, 1 = barycentric mean, 2 = spectral weighted mean
        
        // energy thresholds
        string fEnergyThresholdFileName;
        double fEnergyThresholdFixedValue;
        
        // total numbers
        double fTotalObservationTime;
        double fTotalObservationTimeDeadTimeCorrected;
        double fTotalNormalisationFactor;
        
        // initialize all variables
        void               initializeRunVariables();
        // calculate fluxes and errors
        void               calculateDifferentialFluxes();
        void               calculateDifferentialFluxes_ErrorsAndSignificances( VDifferentialFluxData& iFlux, bool iCombinedAnalysis, double iNormalization );
        // calculat fluxes per run
        vector< VDifferentialFluxData > calculateDifferentialFluxesPerRun( TH1D* i_hhErecCountsOn, TH1D* i_hErecCountsOff, double iAlpha,
                double iObsTimeDeadTimeCorrected, double iObsTime,
                TGraphErrors* iGEff, double iEnergyThreshold_log10TeV );
        // fill and plot energy spectrum graph
        TGraphAsymmErrors* plot_energySpectrum();
        
        /////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////
    public:
    
        VEnergySpectrum();
        VEnergySpectrum( string ifile, string iname = "E", int irun = -1, bool iSourceTypeIsAscii = false );
        ~VEnergySpectrum() {}
        
        double    calculateIntegralFluxFromFitFunction( double iMinEnergy_TeV, double iMaxEnergy_TeV = 1.e6 );
        bool      combineRuns();
        bool      combineRuns( vector< int > iRunList );
        vector< VDifferentialFluxData > getDifferentialFlux()
        {
            return fDifferentialFlux;
        }
        TH1D*     getEnergyCountingOnHistogram()
        {
            return hErecCountsOn;
        }
        TH1D*     getEnergyCountingOffHistogram()
        {
            return hErecCountsOff;
        }
        TGraphAsymmErrors* getEnergySpectrumGraph();
        TH1D* getTotalTimeHistogram( bool iDeadtimeCorrected = false )
        {
            if( iDeadtimeCorrected )
            {
                return hErecTotalTimeDeadTimeCorrected;
            }
            else
            {
                return hErecTotalTime;
            }
        }
        TF1*      getSpectralFitFunction()
        {
            if( fSpectralFitter )
            {
                return fSpectralFitter->getSpectralFitFunction();
            }
            else
            {
                return 0;
            }
        }
        double    getTotalNormalisationFactor()
        {
            return fTotalNormalisationFactor;
        }
        bool      isZombie()
        {
            return bZombie;
        }
        bool      openDataFile( string iFile, int irun = -1, bool iSourceTypeIsAscii = false );
        void      printDifferentialFluxes( bool bSED = false );
        void      printDifferentialFluxes_perRun( unsigned int iEnergyBin, bool bSED = false );
        void      printEnergyThreshold();
        TCanvas*  plot( TCanvas* c = 0 );
        TCanvas*  plotCrabNebulaSpectrum( double iPlottingMultiplierIndex = 0., double i_FitStart_TevLin = 0.1,
                                          double i_FitStop_TeVLin = 20., double i_EnergyBinningLog10 = 0.1,
                                          int i_SpectralFitFunction = 3 );
        TCanvas*  plotCountingHistograms( TCanvas* c = 0 );
        void      plotEnergyThresholds();
        void      plotEventNumbers( Double_t ts = 0.02 );
        void      plotFitValues( bool iForce = false );
        TCanvas*  plotMeanEffectiveArea( TCanvas* c = 0, double i_effMin = 0., double i_effMax = -99. );
        TCanvas*  plotResiduals( TCanvas* c = 0, TF1* f = 0 );
        TCanvas*  plotLifeTimevsEnergy( TCanvas* c = 0 );
        
        void printEnergyBins();
        void setAddHistogramParameters( bool iB = false )
        {
            fAnalysisHistogramAddingUseLowEdge = iB;
        }
        void setEnergyBinning( double iBin = 0.2 );
        bool setEnergyInBinDefinition( unsigned int iE = 2 );
        bool setErrorCalculationMethod( string iErrorMethod = "Rolke" );
        void setReBinningGroupings( vector< int > g )
        {
            newBinningGroupings = g;
        }
        void setRebinBool( bool i_rebin );
        void setEnergyRangeLinear( double xmin, double max );
        void setEnergyRangeLog( double xmin, double xmax );
        void setEnergyThreshold( double a )
        {
            fEnergyThresholdFixedValue = a;
        }
        void setEnergyThreshold( string a )
        {
            fEnergyThresholdFileName = a;
        }
        void setEnergyThresholdDefinition( unsigned int iDef = 1, double iSys = 0.1,
                                           double iMaxEff = 0.2, double iE_min_absolute_TeV = 0.07  )
        {
            fAnalysisEnergyThresholdDefinition = iDef;
            fAnalysisMaxEnergySystematic = iSys;
            fAnalysisMaxEffectiveAreaFraction = iMaxEff;
            fAnalysisEnergyThreshold_min_absolute_TeV = iE_min_absolute_TeV;
        }
        void setSignificanceParameters( double iSig = -5, double iMinOnEvents = -5., double iUpperLimit = 0.95, bool iFit = true,
                                        int iLiAndMa = 17, int iULAlgo = 5 );
                                        
        TF1* fitEnergySpectrum( string iname = "fit", bool bDraw = true );
        void setOffsetdistance( double iOff = -9999. )
        {
            fOffsetDistance = iOff;
        }
        void setSpectralFitFunction( int iD  = 0 )
        {
            fSpectralFitFunction = iD;
        }
        void setSpectralFitFluxNormalisationEnergy( double iS = 1. )
        {
            fSpectralFitFluxNormalisationEnergy = iS;
        }
        void setSpectralFitRangeLin( double xmin = 0.1, double xmax = 10. )
        {
            fSpectralFitEnergy_min = xmin;
            fSpectralFitEnergy_max = xmax;
        }
        void setSpectralFitRangeLog( double xmin = 0.1, double xmax = 10. )
        {
            fSpectralFitEnergy_min = TMath::Power( 10., xmin );
            fSpectralFitEnergy_max = TMath::Power( 10., xmax );
        }
        void setSpectralFitPlottingStyle( int iColor = 1, int iStyle = 1, float iWidth = 2. );
        
        TCanvas* getPlottingCanvas()
        {
            return fPlottingCanvas;
        }
        void setPlottingEnergyRangeLinear( double Emin_TeV = 0.05, double Emax_TeV = 20. )
        {
            fPlottingMinEnergy = Emin_TeV;
            fPlottingMaxEnergy = Emax_TeV;
        }
        void setPlottingEnergyRangeLog( double xmin = -1.2, double xmax = 1.3 )
        {
            fPlottingMinEnergy = TMath::Power( 10., xmin );
            fPlottingMaxEnergy = TMath::Power( 10., xmax );
        }
        void setPlottingLogEnergyAxis( bool iB = true )
        {
            fPlottingLogEnergyAxis = iB;
        }
        void setPlottingMultiplierIndex( double iS = 0. )
        {
            fPlottingMultiplierIndex = iS;
        }
        void setPlottingSpectralWeightForBinCenter( double iS = -2.5 )
        {
            fPlottingSpectralWeightForBinCenter = iS;
        }
        void setPlottingUpperLimits( bool iPlot = true )
        {
            fPlottingUpperLimits = iPlot;
        }
        void setPlottingYaxis( double iMin = 1.e-14, double iMax = 1.e-8 )
        {
            fPlottingYaxisMin = iMin;
            fPlottingYaxisMax = iMax;
        }
        double getIntegralFluxFromFitFunction( double, double );
        double getIntegralFluxErrorFromFitFunction( double, double );
        double getIntegralFluxFromHistogram( );
        double getIntegralFluxErrorFromHistogram( );
        int writeFluxInfoToTextFile( char* output_text_file ) ;
        int writeFitInfoToTextFile( char* output_text_file ) ;
        
        ClassDef( VEnergySpectrum, 17 );
};
#endif
