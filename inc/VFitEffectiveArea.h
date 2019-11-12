//! VFitEffectiveArea  effective area fitting

#ifndef VFITEFFECTIVEAREA_H
#define VFITEFFECTIVEAREA_H

#include "CEffArea.C"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

#include <iostream>
#include <string>

using namespace std;

// fit function
double fun_eff( double* e, double* p );
double fun_effLin( double* e, double* p );

class VFitEffectiveArea
{
    public:
    
        bool bZombie;
        
        string fInputFile;
        
        TF1* fEff;
        TF1* fEffLin;
        TList* hList;
        TGraphAsymmErrors* gEff;
        TGraphAsymmErrors* gEffLog;
        
        bool bAMC;
        int fAMC;
        
        TTree* fTree;
        double fZe;
        int fAz;
        double fAzMin;
        double fAzMax;
        double fIndex;
        double fFitMin;
        double fFitMax;
        
        string fName;
        string fTitle;
        
        TCanvas* cEffLog;
        TH1D* hnullLog;
        TCanvas* cEffLin;
        TH1D* hnullLin;
        double fCanvasMaximum;
        
        VFitEffectiveArea( string ifile );
        void fill();
        void fit( double xmin = -1., double xmax = 1. );
        void getEffectiveArea( double ze = 20., int az = 0, double index = 2.5, bool iAMC = false );
        TGraphAsymmErrors* getLogGraph();
        void makeCanvases();
        void resetFitParameters();
        void setCanvasMaximum( double iB )
        {
            fCanvasMaximum = iB;
        }
        void setLinearFitParameters();
        void setMarkers( TGraph* g, int iMarker, double iSize, int iColor );
        void smoothGraph( double iOut = 0.10, int iIter = 1 );
        void write( string ofile );
};
#endif
