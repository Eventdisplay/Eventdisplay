//! VLikelihoodFitter fit and plot spectral data using likelihood methods

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include "TRandom3.h"
#include <ctime>
#include <vector>
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TFeldmanCousins.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TMatrixD.h"
#include "TLine.h"
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/Factory.h>
#include "Math/ProbFunc.h"
#include "Fit/FitResult.h"
#include <TFitResult.h>
#include "VEnergySpectrum.h"
#include "VHistogramUtilities.h"

// #ifdef __CINT__

// #pragma link off all globals;
// #pragma link off all classes;
// #pragma link off all functions;

#ifndef __VLikelihoodFitter_H_
#define __VLikelihoodFitter_H_

using namespace std;


class VLikelihoodFitter : public VEnergySpectrum
{
    public:
    
    
        VLikelihoodFitter( string filename ) : VEnergySpectrum( filename )
        {
        
            cout << "Constructor\n";
            setEnergyBinning( 0.05 );
            combineRuns();
            setThresholdBias( 0.2 );
            fFitMin = -0.7;
            fFitMax = 1.3;
            fModel_intID = -999;
            setModel( 0 );
            setNormalisationEnergy( 1.0 );
            setBinWidth( 0.1 );
            fExcludeRun.clear();
            fMJD_Min = -999999;
            fMJD_Max = 999999;
            Index_Scan = 0;
            Index_Scan_Likelihood = 0;
            Norm_Scan = 0;
            Norm_Scan_Likelihood = 0;
            fConfidinceInterval = 0;
        }
        
        
        ~VLikelihoodFitter() {}
        
        bool initialize();
        
        
        
        
        
        void setFitMinMax( double i_min, double i_max )
        {
            fFitMin = i_min;
            fFitMax = i_max;
        }
        
        void setBinWidth( double i_BinWidth );
        
        void setBinningRecMC( int i_nRec, double* i_fRecBins, int i_nMC, double* i_fMCBins );
        
        
        void setMJDMinMax( double i_MJDMin, double i_MJDMax )
        {
            fMJD_Min = i_MJDMin;
            fMJD_Max = i_MJDMax;
        }
        
        
        void setThresholdBias( double i_thresh )
        {
            fThresholdBias = i_thresh;
        }
        
        vector <TGraphAsymmErrors*> getEffectiveAreasMC();
        
        
        vector <TGraphAsymmErrors*> getEffectiveAreasRec();
        
        
        vector <TH1D*> getCountingHistogramOn()
        {
            return fOnRebinnedHistograms;
        }
        
        vector <TH1D*> getCountingHistogramOff()
        {
            return fOffRebinnedHistograms;
        }
        
        TH2D* getResponseMatrix( int i )
        {
            return fResponseMatrixRebinned[i];
        }
        
        vector < vector <double> > getOnCounts()
        {
            return fOnCounts;
        }
        
        vector < vector <double> > getOffCounts()
        {
            return fOffCounts;
        }
        
        void setModel( int ID );
        
        void printRunInfo();
        TCanvas* getRunPlots( int i_Entry );
        TCanvas* getTotalCountsPlots();
        TCanvas* plotEnergyBias( int i_Entry );
        
        
        vector < vector <double> > getModelPredictedExcess( const double* parms, double i_eThresh );
        vector < vector <double> > getModelPredictedOff( const double* parms );
        double getLogLi( const double* parms );
        double getLogLiTotal( const double* parms );
        double getLogL0Total( const double* parms );
        TF1* getLikelihoodFit();
        TF1* getLikelihoodFitTotal();
        double getLogL0( const double* parms );
        double getChi2( const double* parms );
        double getChi2Total( const double* parms );
        double getNDF();
        
        
        void setNormalisationEnergy( double iNormEnergy )
        {
            E_Norm = iNormEnergy;
        }
        
        // Set Redshift
        void setRedShift( double i_Z );
        
        // Set the EBL Opacity
        int setEBLOpacity( TGraph2D* i_OpacityTable );
        
        // Setting the intrinsic Source Model
        void setIntrinsicModel( int i_ID );
        
        TGraphErrors* calculateConfidinceInterval( double* i_covmat, TF1* i_fitfunction, int i_model, int i_nparms );
        
        TGraphErrors* getConfidinceInterval()
        {
            return fConfidinceInterval;
        }
        TGraphAsymmErrors* getEnergySpectrum( TF1* iBestFit );
        // TGraphAsymmErrors* getEnergySpectrumTotal(TF1 *iBestFit);
        
        
        void printCountsInBins();
        void addExclusionDate( double i_MJDStart, double i_MJDStop );
        void excludeRun( int i_Run )
        {
            fExcludeRun.push_back( i_Run );
        }
        float* getIntergralFlux( double i_EMin, double i_EMax, TF1* i_Model, bool i_log = true );
        void setEBLOpacity( TGraph* i_EBLOpacity );
        
        
        double* getIndexScan()
        {
            return Index_Scan;
        }
        double* getIndexScanLikelihood()
        {
            return Index_Scan_Likelihood;
        }
        
        
        double* getNormScan()
        {
            return Norm_Scan;
        }
        double* getNormScanLikelihood()
        {
            return Norm_Scan_Likelihood;
        }
        
        
        double getRunwiseLogL( int i_run, const double* parms );
        double getRunwiseLogL0( int i_run, const double* parms );
        int getRunWiseFitInfo( int i_runnum, TF1* i_fit );
        void getRunwiseFitStatus( TF1* i_fit );
        
        double getBinTS( double i_on, double i_off, double i_alpha );
    private:
        TRandom3* fRandom;
        char hname[800];
        
        // Binning Information
        int nRec;
        double* fRecBins;
        int nMC;
        double* fMCBins;
        double fBinWidth;
        vector < double > fRecBinCentres;
        vector < double > fMCBinCentres;
        int nBinsFit;
        int nBinsFit_Total;
        int nBinsFit_runwise;
        double fThresholdBias;
        double E_Norm;
        double fMJD_Min;
        double fMJD_Max;
        vector < vector <double> > fExcludeMJD;
        vector <int> fExcludeRun;
        
        vector <double> fLastOn;
        vector <double> fLastOff;
        // Model Information
        int fModelID;
        int fModel_intID;
        TF1* fModel;
        TF1* fModel_int;
        int nParms;
        double fFitMin;
        double fFitMax;
        
        TGraphErrors* fConfidinceInterval;
        
        // Totals
        vector <int> fTotalOn;
        vector <int> fTotalOff;
        
        // Data Vectors
        vector <TH1D*> fOnRebinnedHistograms;
        vector <TH1D*> fOffRebinnedHistograms;
        vector <TH2D*> fResponseMatrixRebinned;
        vector <TGraphAsymmErrors*> fMeanEffectiveAreaRec;
        vector <TGraphAsymmErrors*> fMeanEffectiveAreaMC;
        
        
        vector < vector <double> > fOnCounts;
        vector < vector <double> > fOffCounts;
        vector < vector <double> > fModelPredictedExcess;
        vector < vector <double> > fModelPredictedOff;
        vector < vector <double> > fEnergyBias;
        
        vector < int > vLastOnCount;
        
        
        // EBL analysis
        double Z ;
        TGraph* fEBLOpacityGraph;
        // Accessing Raw Data
        vector <TH1D*> getCountingHistogramOnRaw();
        vector <TH1D*> getCountingHistogramOffRaw();
        vector <TH2D*> getResponseMatrixRaw();
        vector < vector <double> > getCounts( vector <TH1D*> i_hTemp );
        
        // Fit fuction that takes into account EBL attenuation
        double calculateIntrinsicSpectrum( Double_t* x, Double_t* parm );
        
        float* getSpectralPoint( double* parms , double BinMin, double BinMax, double iE_Norm, TF1* iBestFit );
        double* getSpectralPointTotal( double* parms , double BinMin, double BinMax, double iE_Norm );
        double getCrabFlux( double iF, double i_EMin, double i_EMax, double i_Gamma = -2.49 );
        double* Index_Scan;
        double* Index_Scan_Likelihood;
        double* Norm_Scan;
        double* Norm_Scan_Likelihood;
        
        
        
        
        
};

#endif
