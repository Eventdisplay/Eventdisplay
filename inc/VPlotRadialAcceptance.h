//! VPlotRadialAcceptance plot radial acceptance curves
//

#ifndef VPlotRadialAcceptance_H
#define VPlotRadialAcceptance_H

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TText.h"

#include <iostream>
#include <string>
#include <vector>

#include "VHistogramUtilities.h"
#include "VPlotUtilities.h"
#include "VUtilities.h"

using namespace std;

class VPlotRadialAcceptance : public VPlotUtilities
{
    private:
    
        bool fDebug;
        
        string fName;
        
        TFile* fAcceptanceFile;
        TH1F*  fAcceptanceHisto;
        TH2F*  fAcceptanceHisto2D_rot;
        bool   fAcceptanceHisto2D_rot_normalized;
        TH2F*  fAcceptanceHisto2D_derot;
        bool   fAcceptanceHisto2D_derot_normalized;
        TF1*   fAcceptanceFunction;
        vector< TH1F* > fAcceptancePhiHisto;
        vector< TF1* >  fAcceptancePhiFitFunction;
        vector< TH1F* > fAcceptancePhiHistoDeRot;
        vector< TF1* >  fAcceptancePhiFitFunctionDeRot;
        TH1F*  hPhiDist;
        TH1F*  hPhiDistDeRot;
        
        double fAxis_x_min;
        double fAxis_x_max;
        double fAxis_y_min;
        double fAxis_y_max;
        
        void     plotRadialAcceptance2D_plotter( TH2F* h, bool iRotated,
                double iNormalizationRadius = 0.35,
                bool iSmooth = false, float iMaxRadius = 2.,
                TPad* c_h2D = 0, TPad* c_h2Dresidual = 0 );
                
    public:
    
        VPlotRadialAcceptance( string iFile = "", int iAzBin = -1 );
        ~VPlotRadialAcceptance() {}
        
        TF1*  getAcceptanceFunction()
        {
            return fAcceptanceFunction;
        }
        TH1F* getAcceptanceHisto()
        {
            return fAcceptanceHisto;
        }
        
        //TCanvas* plotRadialAcceptance( TCanvas* cX = 0 );
        TCanvas* plotRadialAcceptance( TCanvas* cX = 0 , int iColor = 1 );
        bool     plotRadialAcceptance_perRun( bool iPlotTotalAcceptanceFunction = false );
        void     plotRadialAcceptance2D( bool iRotated = false, double iNormalizationRadius = 0.35, bool iSmooth = false, float iMaxRadius = 2. );
        bool     plotRadialAcceptance2D_perRun( double iNormalizationRadius = 0.35, bool iSmooth = false, float iMaxRadius = 2. );
        TCanvas* plotPhiDependentRadialAcceptances( TCanvas* cX = 0, int iIterator = 4, bool iDeRot = false );
        TCanvas* plotPhiDistributions( TCanvas* cX = 0, int iColor = 1 );
        TCanvas* plotResiduals( TCanvas* cX = 0, double i_res_min = -0.5, double i_res_max = 0.5, bool iDrawChi2 = true );
        bool     openAcceptanceFile( string iFile, int iAzBin = -1 );
        void     setAxisRange( double x_min = 0., double x_max = 2.5, double y_min = 0., double y_max = 1.5 );
        void     setName( string iName )
        {
            fName = iName;
        }
        
        ClassDef( VPlotRadialAcceptance, 1 );
};

#endif
