//! VHistogramUtilities.h  utility class to manipulate histograms

#ifndef VHistogramUtilities_H
#define VHistogramUtilities_H

#include "TDirectory.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TMath.h"
#include "TProfile.h"
#include "TRandom.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VHistogramUtilities
{
    private:
    
        bool fDebug;
        
        
    public:
    
        VHistogramUtilities();
        ~VHistogramUtilities() {}
        
        static bool      divide( TGraphAsymmErrors* g, TGraphAsymmErrors* g1, TGraph* g2 );
        static bool      divide( TGraphAsymmErrors* g, TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, double epsilon = 1.e-3 );
        int              findBinInGraph( TGraph* g, double x );
        void             fill_Graph_in_Histogram( TGraphAsymmErrors* g, TH1* h, bool bRequireYGtZero = false );
        TH1D*            get_Bin_Distribution( TH2D* h, int ion, double rmax, double rSource, bool iDiff, TH2D* hTest,
                                               int iN = 0, float* x = 0, float* y = 0, float* r1 = 0, float* r2 = 0, float* theta = 0, string regioncode = "" );
        TH1D*            get_Cumulative_Histogram( TH1D* iH_in, bool iNormalize, bool iLeft_to_right, double i_binvalue = 1.e30, double i_min_value = 0. );
        bool             get_Graph_from_Histogram( TH1* h, TGraphErrors* g, bool bIgnoreErrors = false, double iMinBinContent = 0., double iXmin = -1.e50, double iXmax = 1.e50 );
        bool             get_Graph_from_Histogram( TH1* h, TGraphAsymmErrors* g, bool bIgnoreErrors = false, bool bLinXaxis = false, double iCutUnrealisticErrors = -1., double iXmin = -1.e50, double iXmax = 1.e50 );
        TGraphErrors*    get_Profile_from_TH2D( TH2D* iH, TGraphErrors* g = 0, string iM = "median", int i_rebin = 2,
                                                double iXaxisValue = -10., double iMinusValue = 0. );
        static TH1D*     get_ResidualHistogram_from_TF1( string iname = "", TH1* h = 0, TF1* f = 0 );
        
        TH1F*            get_CTA_IRF_Histograms( string iHistogramName, double iCameraOffset );
        TH1F*            get_CTA_IRF_Histograms_from2D( string iHistogramName,  double iSummand = 0. );
        
        static TH1*      normalizeTH1( TH1* h, bool iIntegral );
        
        static TH2F*     reduce2DHistogramSize( TH2* h, string inewHistogramName );
        static TH2D*     calculateContainmentDistance( TH2D* h, string inewHistogramName );
        static double    interpolateTH2D( TH2* h, double x, double y );
        
        static TH2*      interpolateResponseMatrix( TH2*, string iNewHistoName = "" );
        static bool      normalizeTH2D_y( TH2* h );
        
        void             setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        
        ClassDefNV( VHistogramUtilities, 12 );
};

#endif

