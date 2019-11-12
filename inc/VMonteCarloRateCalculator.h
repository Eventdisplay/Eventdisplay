//! VMonteCarloRateCalculator calculate rates from Monte Carlo

#ifndef VMonteCarloRateCalculator_H
#define VMonteCarloRateCalculator_H

#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"

#include "CEffArea.h"
#include "VEnergySpectrumfromLiterature.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VMonteCarloRateCalculator : public VPlotUtilities
{
    private:
    
        TFile* fFile;
        TTree* fMCTree;
        double fze;
        int    faz;
        double fwoff;
        int    fnoise;
        unsigned int fnrates;
        double fMCrate[1000];
        
        vector< double > fenergy;
        vector< double > feffectiveArea;
        
        void getMinMaxRates( unsigned int n, double* r, double& i_mean, double& i_min, double& i_max );
        
    public:
        VMonteCarloRateCalculator();
        VMonteCarloRateCalculator( string iFile );
        virtual ~VMonteCarloRateCalculator() {}
        
        // fill rates
        double getMonteCarloRate( int nbins, double* e0, double* eff, double i_gamma, double i_phi, double iE0 = 1.,
                                  double iEMin = 1.e-20, double iEMax = 1.e20, double bDebug = false );
        double getMonteCarloRate( vector< double >& e, vector< double >& eff, double i_gamma, double i_phi, double iE0 = 1.,
                                  double iEMin = 1.e-20, double iEMax = 1.e20, bool bDebug = false );
        double getMonteCarloRate( int nbins, double* e0, double* eff, VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
                                  double iEMinBin, double iEMaxBin, bool bDebug = false );
        double getMonteCarloRate( vector< double > e, vector< double > eff, VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
                                  unsigned int iEMinBin, unsigned int iEMaxBin, TH2* iResponseMatrix = 0, bool bDebug = false );
        double getMonteCarloRate( vector< double > e, vector< double > eff, VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
                                  unsigned int iEMinBin, unsigned int iEMaxBin, double iEMin, double iEMax, TH2* iResponseMatrix = 0, bool bDebug = false );
        double getMonteCarloRate( vector< double > e, vector< double > eff, VEnergySpectrumfromLiterature* e_lit, unsigned int e_lit_ID,
                                  unsigned int iEMinBin, unsigned int iEMaxBin, double iEMin, double iEMax, TH2* iResponseMatrix,
                                  vector< double > e_gamma, bool bDebug );
                                  
        // read rates
        TCanvas* plot_MonteCarloRate_vs_wobbleOffsets( TCanvas* c = 0, double ze = 20., int az = 16, int noise = 200, string iPlottingOption = "3" );
        TGraphAsymmErrors* getMonteCarloRate_vs_wobbleOffsets( double ze = 20., int az = 16, int noise = 200 );
        
        ClassDef( VMonteCarloRateCalculator, 3 );
};
#endif
