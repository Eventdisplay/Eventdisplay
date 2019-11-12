//! VEnergyThreshold  class to determine energy threshold

#ifndef VEnergyThreshold_H
#define VEnergyThreshold_H

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "TObject.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"

#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "VRunList.h"
#include "CEffArea.h"

using namespace std;

class VEnergyThreshold : public TObject
{
    private:
    
        bool fDebug;
        
        TFile* fEffAreaFile;
        CEffArea* fEffArea;
        
        TFile* fEnergyThresholdFile;
        
        double fEnergyThresholdFixedValue;
        string fEnergyThresholdFileName;
        
        TFile* fOutFile;
        TTree* fTreeEth;
        
        double fze;
        int fAzBin;
        double fXoff;
        double fYoff;
        double fWoff;
        int fTNoise;
        double fTNoisePE;
        double fTPedvar;
        double fSpectralIndex;
        double fdiffmax;
        double fesys_10p;
        double fesys_15p;
        double fesys_20p;
        double feffFract_05p;
        double feffFract_10p;
        double feffFract_20p;
        double feffFract_50p;
        double feffFract_90p;
        
        int fPlottingMarkerStyle;
        int fPlottingMarkerColor;
        float fPlottingMarkerSize;
        float fPlottingLineWidth;
        double fPlottingYmin;
        double fPlottingYmax;
        
        bool   openEnergyThresholdFile();
        bool   setUpThresholdTree();
        void   copyEntry();
        
        double interpolateEnergyThreshold( VRunList* );
        
    public:
    
        VEnergyThreshold();
        VEnergyThreshold( string ioutfilename );
        VEnergyThreshold( double iEthFixed, string iInFile = "" );
        ~VEnergyThreshold() {}
        bool   openEffectiveAreaFile( string ifile );
        bool   calculateEnergyThreshold( int nentries = -1 );
        double getEnergy_diffMax( TObject* h = 0, double index = -2.5 );
        double getEnergy_maxSystematic( TObject* h = 0, double iSys = 0.1, double iE_min_absolute_TeV = 0.07 );
        double getEnergy_MaxEffectiveAreaFraction( TObject* h = 0, double iFrac = 0.1 );
        double getEnergy_fixedValue()
        {
            return fEnergyThresholdFixedValue;
        }
        void plot_energyThresholds( string var = "E_diffmax", double ze = 20., double woff = 0.5,
                                    int noise = 150, double index = 2.4, int az = 16,
                                    bool bPlot = true, string plot_option = "p" );
        void setPlottingStyle( int iC = 1, int iS = 21, float iW = 2., float iL = 2. )
        {
            fPlottingMarkerStyle = iS;
            fPlottingMarkerColor = iC;
            fPlottingMarkerSize = iW;
            fPlottingLineWidth = ( Width_t )iL;
        }
        void setPlottingYaxis( double ymin = -1., double ymax = -1. )
        {
            fPlottingYmin = ymin;
            fPlottingYmax = ymax;
        }
        bool writeResults();
        
        double getEnergyThreshold( VRunList* );
        
        ClassDef( VEnergyThreshold, 3 );
};
#endif
