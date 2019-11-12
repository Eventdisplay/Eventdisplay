//! VPlotMonteCarloQualityFactor fill and plot quality factors for MC

#ifndef VPlotMonteCarloQualityFactor_H
#define VPlotMonteCarloQualityFactor_H

#include "TCanvas.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLine.h"

#include "VPlotUtilities.h"
#include "CData.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VPlotMonteCarloQualityFactorData
{
    public:
    
        double fVar_max;
        double fVar_min;
        TH1D*  hSignal;
        TH1D*  hBackground;
        TH1D*  hQFactors_LowerCut;
        TH1D*  hQFactors_UpperCut;
        TGraphErrors* gQFactor_LowerCutE;
        TGraphErrors* gQFactor_UpperCutE;
        TGraphErrors* gQFactorMax_LowerCutE;
        TGraphErrors* gQFactorMax_UpperCutE;
        
        VPlotMonteCarloQualityFactorData();
        ~VPlotMonteCarloQualityFactorData() {}
};


class VPlotMonteCarloQualityFactor : public VPlotUtilities
{
    private:
    
        bool fDebug;
        
        double fEnergy_min;
        double fEnergy_max;
        double fMSCW_min;
        double fMSCW_max;
        double fMSCL_min;
        double fMSCL_max;
        
        CData* fSignalChain;
        CData* fBackgroundChain;
        
        map< string, VPlotMonteCarloQualityFactorData* > fData;
        
        void calculateQfactors();
        void fill( int iMaxNevents, CData* c, bool bSignal );
        void initializeHistograms();
        void resetHistograms();
        bool setDataChain( string iChain, bool bSignal );
        
    public:
    
        VPlotMonteCarloQualityFactor();
        ~VPlotMonteCarloQualityFactor() {}
        
        bool setBackgroundDataChain( string iChain );
        bool setSignalDataChain( string iChain );
        void setMSCCuts( double iMSCW_min = -2., double iMSCW_max = 0.5, double iMSCL_min = -2., double iMSCL_max = 0.5 );
        void setEnergyRange( double iEmin = -10., double iEmax = 10. )
        {
            fEnergy_min = iEmin;    // log10 [TeV]
            fEnergy_max = iEmax;
        }
        bool fill( int iMaxEvents = -1 );
        void fillEnergyDependence( int iMaxNevents = -1, double iEmin = -2., double iEmax = 2., double iEbin = 0.25 );
        void plot( bool iPrint = false );
};

#endif
