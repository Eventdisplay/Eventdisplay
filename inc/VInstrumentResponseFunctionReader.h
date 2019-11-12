//! VInstrumentResponseFunctionReader data class for reading response functions (effective area, angular resolution, etc )

#ifndef VInstrumentResponseFunctionReader_H
#define VInstrumentResponseFunctionReader_H

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"

#include "CEffArea.h"
#include "VAnalysisUtilities.h"
#include "VGammaHadronCuts.h"
#include "VHistogramUtilities.h"
#include "VInstrumentResponseFunctionData.h"
#include "VPlotUtilities.h"

using namespace std;

class VInstrumentResponseFunctionReader : public VAnalysisUtilities, public VPlotUtilities, public VHistogramUtilities
{
    private:
    
        bool   fIsZombie;
        bool   fDebug;
        
        unsigned int fEnergyResolutionMethod;   // 0 = RMS, 1 = Fermi Method (68% around mprop)
        bool         fEnergyXaxisIsEtrue;
        float        fEnergy_TeV_Eres_oneSided;
        
        int    fGammaHadronCuts_directionCut_selector;
        
        bool                             calculateCutEfficiencies();
        TGraphAsymmErrors*               calculateEffectiveAreaRatios( TGraphAsymmErrors* g0, TGraphAsymmErrors* g1 );
        TGraphErrors*                    getAngularResolutionGraphs( string iHName, string iHNameOldStyle );
        VInstrumentResponseFunctionData* getIRFFromFile( TTree* );
        bool                             getDataFromFile();
        bool                             getDataFromCTAFile();
        TGraphErrors*                    getResolutionGraph( TGraphErrors* );
        TGraphErrors*                    getAngularResolution( TH2D* iHistogram, double iContainmentProbability, int iMinRequiredEvents = 20 );
        TGraphErrors*                    getEnergyResolutionMPropInterval( TH2D* migmatrix,
                bool bXaxisIsEtrue = false,
                TGraphErrors* gEnergyResolutionNoDirectionCuts = 0 );
        void                             getEnergyResolutionPlot( TProfile* iP, int i_rebin = 2, double iMinEnergy = -10. );
        void                             getEnergyResolutionPlot( TH2D* iP, double iMinEnergy = -10. );
        void                             getEnergyResolutionPlot68( TH2D* iP, double iReferenceValue = -999. );
        
        bool                             initializeIRFData();
        
    public:
    
        string fFile;
        string fA_MC;
        
        //////////////////////////////////
        // conditions
        //////////////////////////////////
        double fZe;
        double fWoff;
        int    fAzbin;
        double fIndex;
        int    fNoise;
        
        float  fEnergyLinTeV_min;
        float  fEnergyLinTeV_max;
        
        //////////////////////////////////
        // data
        //////////////////////////////////
        // effective areas
        TGraphAsymmErrors* gEffArea_MC;
        TGraphAsymmErrors* gEffArea_Rec;
        TGraphAsymmErrors* gEffArea_Recp80;
        TGraphAsymmErrors* gEffArea_MCnoTh2;
        TGraphAsymmErrors* gEffArea_RecnoTh2;
        // effective area ratios
        TGraphAsymmErrors* gEffArea_MC_Ratio;
        TGraphAsymmErrors* gEffArea_Rec_Ratio;
        TGraphAsymmErrors* gEffArea_MCnoTh2_Ratio;
        TGraphAsymmErrors* gEffArea_RecnoTh2_Ratio;
        // energy spectra
        TH1D* hEmc;
        TH1D* hEcut;
        TH1D* hEcut_rec;
        TH1D* hEcutUW;
        TH1D* hEcut_recUW;
        // energy reconstruction matrix
        TH2D* hERecMatrix;
        TH2D* hERecMatrixCoarse;
        TH2D* hERecMatrixQC;
        TH2D* hERecMatrixCoarseQC;
        TH2D* hERecMatrixNoDirectionCuts;
        TH2D* hERecMatrixCoarseNoDirectionCuts;
        // e_rec/e_mc
        TH2D* hEsysMCRelative2D;
        TH2D* hEsysMCRelative2DNoDirectionCut;
        TProfile* hEsysMCRelative;
        // 2D energy error distribution
        TH2D* hEsys;
        // energy resolution
        TGraphErrors* gEnergyResolution;
        TGraphErrors* gEnergyResolutionQC;
        TGraphErrors* gEnergyResolutionNoDirectionCuts;
        // energy bias
        TGraphErrors* gEnergyBias_Mean;
        TGraphErrors* gEnergyBias_Median;
        TGraphErrors* gEnergyLogBias_Mean;
        TGraphErrors* gEnergyLogBias_Median;
        // angular resolution (filled for CTA only)
        TGraphErrors* gAngularResolution;               // 68% containment radius for angular resolution
        TGraphErrors* gAngularResolution80;             // 80% containment radius for angular resolution
        TGraphErrors* gAngularResolution95;             // 95% containment radius for angular resolution
        TH2D* h2DAngularPSF;                            // angular difference histogram
        TH2D* h2DAngularPSFEMC;                         // angular difference histogram (true energy axis)
        // cut efficiencies
        vector< TH1D* > hCutEfficiency;
        vector< TH1D* > hCutEfficiencyRelativePlots;
        // weight histograms
        TH1D* hWeightedRate;
        TH1D* hWeightedRate005;
        
        // resolution graphs
        vector< string >                           fIRF_TreeNames;
        vector< VInstrumentResponseFunctionData* > fIRF_Data;
        
        //////////////////////////////////
        // plotting
        //////////////////////////////////
        string fPlotOption;
        int    fColor;
        int    fLineStyle;
        int    fMarkerStyle;
        string fLegend;
        
        VInstrumentResponseFunctionReader();
        ~VInstrumentResponseFunctionReader() {}
        
        bool calculateEffectiveAreaRatios( TGraphAsymmErrors* g0 );
        bool fillBiasHistograms( TH1F* h = 0, string iMeanOrMedian = "mean" );
        bool fillEffectiveAreasHistograms( TH1F* h = 0, string iContainmentRadius = "", TH1F* hMC = 0 );
        bool fillResolutionHistogram( TH1F* h = 0, string iContainmentRadius = "68",
                                      string iResolutionTreeName = "t_angular_resolution",
                                      bool iEnergyAxis_MC = false );
        bool fillData();
        bool fillData( string iDataLine, int iDataID );
        bool fillData( string iFile, double iZe = 20., double iWoff = 0.5, int iAzBin = 0, double iIndex = 2.0, int iNoise = 200, string iA_MC = "A_MC" );
        TH2D* get2DAngularPSF( bool iEAxisMC = false )
        {
            if( iEAxisMC )
            {
                return h2DAngularPSFEMC;
            }
            return h2DAngularPSF;
        }
        TH2D* getRecvsMCEnergy( bool iNoDirectionCut = false )
        {
            if( iNoDirectionCut )
            {
                return hEsysMCRelative2DNoDirectionCut;
            }
            return hEsysMCRelative2D;
        }
        TH2D* getMigrationMatrix()
        {
            return hERecMatrix;
        }
        TH2D* getMigrationMatrixNoDirectionCuts()
        {
            return hERecMatrixNoDirectionCuts;
        }
        bool isZombie()
        {
            return fIsZombie;
        }
        void setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        void setEnergyRange( float iEmin_linTeV = 1.e-6, float iEmax_linTeV = 1.e6 )
        {
            fEnergyLinTeV_min = iEmin_linTeV;
            fEnergyLinTeV_max = iEmax_linTeV;
        }
        void setEnergyResolutionMethod( unsigned int iEmethod, bool iXaxis = true, float iEnergy_TeV_Eres_oneSided = -99. )
        {
            fEnergyResolutionMethod   = iEmethod;
            fEnergyXaxisIsEtrue       = iXaxis;
            fEnergy_TeV_Eres_oneSided = iEnergy_TeV_Eres_oneSided;
        }
        
        void setPlotOption( string iPlotOption = "pl" )
        {
            fPlotOption = iPlotOption;
        }
        
        ClassDef( VInstrumentResponseFunctionReader, 16 );
};


#endif
