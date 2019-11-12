//! VStereoReconstruction.h calculate and plot array and direction resolution

#ifndef VStereoReconstruction_h
#define VStereoReconstruction_h

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"

#include "VPlotUtilities.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class VStereoReconstructionData
{
    public:
    
        string  fTitle;
        int     fColor;
        int     fMarkerStyle;
        float   fMarkerSize;
        int     fLineStyle;
        Width_t fLineWidth;
        int     fFillStyle;
        string  fPlotStyle;
        
        TGraphAsymmErrors* gData;
        
        VStereoReconstructionData();
        ~VStereoReconstructionData() {}
        void   draw();
};


//###########################################################################

class VStereoReconstruction : public VPlotUtilities
{
    private:
    
        string fName;
        bool bDebug;
        
        // vector containing all the plotting data
        vector< VStereoReconstructionData* > fData;
        
        // plotting variables
        string fPlottingVariable;
        TCanvas* fPlottingCanvas;
        // spectral weight to calculate bin centers
        double fPlottingYaxisMin;
        double fPlottingYaxisMax;
        double fPlottingMinEnergy;                // linear energy axis [TeV]
        double fPlottingMaxEnergy;                // linear energy axis [TeV]
        bool   fPlottingLogEnergyAxis;            // plot log or lin values on energy axis (default=true)
        int    fPlottingCanvasSizeX;
        int    fPlottingCanvasSizeY;
        
    public:
    
        VStereoReconstruction();
        ~VStereoReconstruction() {}
        
        bool      addDataSet( string iTitle, string iFile, double emin_lin = -1., double emax_lin = -1., bool bEnergyAxis_linear_GeV = false, bool bResolutionAxis_arcmin = false );   // linear energy axis [TeV]
        bool      readDataSetsfromTextFile( string iFile, unsigned int iSet, bool bClearExistingDataSet = true );
        bool      removeDataSet( unsigned int iDataSet );
        bool      setPlottingAtt( unsigned int iDataSet, string iPlotStyle = "pc", int iColor = 1, int iMarkerStyle = 20, float iMarkerSize = 1., int iLineStyle = 1, float iLineWidth = 1., int iFillStyle = 1001 );
        
        TCanvas*  plot( TCanvas* c = 0 );
        void      plotLegend();
        void      setDebug( bool iB = false )
        {
            bDebug = iB;
        }
        void      setName( string a = "VStereoReconstruction" )
        {
            fName = a;
        }
        bool      setPlottingVariable( string iVar = "direction" );
        TCanvas*  getPlottingCanvas()
        {
            return fPlottingCanvas;
        }
        void      setPlottingCanvasSize( int iX = 600, int iY = 600 )
        {
            fPlottingCanvasSizeX = iX;
            fPlottingCanvasSizeY = iY;
        }
        void      setPlottingEnergyRangeLinear( double xmin = 0.01, double xmax = 200. )
        {
            fPlottingMinEnergy = xmin;
            fPlottingMaxEnergy = xmax;
        }
        void      setPlottingEnergyRangeLog( double xmin = -2.0, double xmax = 2.0 )
        {
            fPlottingMinEnergy = TMath::Power( 10., xmin );
            fPlottingMaxEnergy = TMath::Power( 10., xmax );
        }
        
        void      setPlottingLogEnergyAxis( bool iB = true )
        {
            fPlottingLogEnergyAxis = iB;
        }
        void      setPlottingYaxis( double iMin = 0., double iMax = 0.2 )
        {
            fPlottingYaxisMin = iMin;
            fPlottingYaxisMax = iMax;
        }
        
};


#endif
