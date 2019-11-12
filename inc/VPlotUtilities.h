//! VPlotUtilities some utilities to plot stuff

#ifndef VPlotUtilities_H
#define VPlotUtilities_H

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TMath.h"
#include "TObject.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TColor.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VPlottingAxisData
{
    public:
    
        string fName;
        string fAxis;           // X, Y, Z, ...
        string fAxisTitle;
        bool   fLogAxis;
        double fMinValue;
        double fMaxValue;
        
        VPlottingAxisData();
        virtual ~VPlottingAxisData() {}
        
        ClassDef( VPlottingAxisData, 1 );
};

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

class VPlotUtilities
{
    protected:
    
        int    fPlottingCanvasX;
        int    fPlottingCanvasY;
        int    fPlottingColor;
        int    fPlottingMarkerStyle;
        double fPlottingMarkerSize;
        int    fPlottingLineStyle;
        double fPlottingLineWidth;
        int    fPlottingFillStyle;
        
        int    fColorAxis_ncolor;
        int    fColorAxis_ncont;
        double fColorAxis_vmin;
        double fColorAxis_vmax;
        TGaxis* fColorAxis_axis;
        
        string fPlotHistogramTitle;
        double fPlotHistogramTitle_x;
        double fPlotHistogramTitle_y;
        
        map< string, VPlottingAxisData* > fPlottingAxisData;
        
        
    public:
    
        VPlotUtilities();
        virtual ~VPlotUtilities();
        void default_settings();
        void plot_nullHistogram( TPad* c, TH1* h, bool bLogX = false, bool bLogY = false, double yTitleOffset = 1., double xmin = 0., double xmax = 0. );
        void plotHistogramTitle( TH1* h );
        
        TGaxis* getColorAxisAxis( double x1 = 0.9, double x2 = 0.93, double y1 = 0.7, double y2 = 0.95, string AxisTitle = "", Int_t ndiv = 10, string iOption = "+L" );
        VPlottingAxisData* getPlottingAxis( string iName );
        int    getPlottingCanvasX()
        {
            return fPlottingCanvasX;
        }
        int    getPlottingCanvasY()
        {
            return fPlottingCanvasY;
        }
        int    getPlottingColor()
        {
            return fPlottingColor;
        }
        int    getPlottingMarkerStyle()
        {
            return fPlottingMarkerStyle;
        }
        double getPlottingMarkerSize()
        {
            return fPlottingMarkerSize;
        }
        int    getPlottingLineStyle()
        {
            return fPlottingLineStyle;
        }
        double getPlottingLineWidth()
        {
            return fPlottingLineWidth;
        }
        int    getPlottingFillStyle()
        {
            return fPlottingFillStyle;
        }
        int    getColorAxisColor( double iV );
        unsigned int listPlottingAxis();
        void   setCanvasSize( int x = 600, int y = 400 )
        {
            fPlottingCanvasX = x;
            fPlottingCanvasY = y;
        }
        void   setColorAxisDataVector_minmax( double imin = 0., double imax = 0. );
        void   setColorAxisPalette( int palette = 55, int ncolors = 100 );
        void   setArrowPlottingStyle( TArrow* iArr );
        void   setArrowPlottingStyle( TArrow* iArr, int icolor, double iwidth = 1., int iLineStyle = 1 );
        void   setFunctionPlottingStyle( TF1* iF );
        void   setFunctionPlottingStyle( TF1* iF, int icolor, double iwidth = 1., double isize = 1., int imarker = 1, int iFillStyle = 1 );
        void   setHistogramPlottingStyle( TH1* his );
        void   setHistogramPlottingStyle( TH1* his, int icolor, double iwidth = 1., double isize = 1., int imarker = 1, int irebin = 1, int iFillStyle = 1 );
        void   setHistogramPlottingStyle( TH2D* his, double iminF );
        void   setGraphPlottingStyle( TGraph* g );
        void   setGraphPlottingStyle( TGraph* g, int icolor, double iwidth = 1., int imarker = 20, double isize = 1., int iFillStyle = 0, int iLineStyle = 1 );
        unsigned int setPlottingAxis( string iName = "energy_Lin", string iAxis = "X", bool iLog = true, double xmin = 0.005, double xmax = 200., string iAxisTitle = "" );
        void   setPlottingStyle( int iColor = 1, int iStyle = 1, float iWidth = 2., int iMarkerStyle = 21,
                                 float iMarkerSize = 1., int iFillStyle = 1001 )
        {
            fPlottingColor = iColor;
            fPlottingLineStyle = iStyle;
            fPlottingLineWidth = iWidth;
            fPlottingMarkerStyle = iMarkerStyle;
            fPlottingMarkerSize = iMarkerSize;
            fPlottingFillStyle = iFillStyle;
        }
        void   setPadMargins( TCanvas* c, int nPads, double lM, double rM );
        void   setPlotHistogramTitle( string iS = "", double x = -1., double y = -1. )
        {
            fPlotHistogramTitle = iS;
            fPlotHistogramTitle_x = x;
            fPlotHistogramTitle_y = y;
        }
        void   setTitles( TH1* his, string iname, string ititle, string ytitle );
        void   setBlackAndWhitePalette();
        TH2D*  removeOuterRing( TH2D* h, double r, double ivalue );
        
        ClassDef( VPlotUtilities, 13 );
};
#endif
