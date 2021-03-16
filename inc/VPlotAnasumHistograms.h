//! VPlotAnasumHistograms  plot anasum results

#ifndef VPlotAnasumHistograms_h
#define VPlotAnasumHistograms_h

#include "VAstronometry.h"
#include "VAnalysisUtilities.h"
#include "VEnergyThreshold.h"
#include "VHistogramUtilities.h"
#include "VMathsandFunctions.h"
#include "VPlotUtilities.h"
#include "VSkyCoordinatesUtilities.h"
#include "VStatistics.h"
#include "CRunSummary.h"

#include "VStarCatalogue.h"

#include "TArrow.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TTimeStamp.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

struct sSource
{
    double fX;
    double fY;
    string fStarName;
};

class VPlotAnasumHistograms : public VAnalysisUtilities, public VPlotUtilities, public VHistogramUtilities
{
    private:
    
        string fAnasumDataFile;    //! name of input anasum file
        int    fRunNumber;         //! -1 all runs (run number for plotting can be set later)
        bool   fDebug;
        string fPlotMode;
        
        bool   fPlotCorrelated;    //! plot correlated sky plots
        float  fPlotDrawPSF;       //! plot PSF on top of sky maps (radius in deg)
        bool   fPlotUseHours;      //! plot hour/min/sec on sky maps axis
        int    fPlotZeroHours;     //! quick fix for hour axis (definition might depend on time zone)
        
        // some run Summary info needed for skymaps
        double fSkyMapCentreDecJ2000;
        double fSkyMapCentreRAJ2000;
        double fTargetShiftWest;
        double fTargetShiftNorth;
        
        // histograms
        TH1D* hmscw_on;
        TH1D* hmscw_off;
        TH1D* hmscw_diff;
        TH1D* hmscl_on;
        TH1D* hmscl_off;
        TH1D* hmscl_diff;
        
        TH1D* htheta2_on;
        TH1D* htheta2_off;
        TH1D* htheta2_diff;
        
        // helper functions
        TH1D* doQfactors( TH1D* hon, TH1D* hoff, TH1D* hdiff, bool bUpper, int iMethod, double iSourceStrength );
        TH2D* reflectXaxis( TH2D* h = 0, char* iNewName = 0 );
        
        /////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////
    public:
    
        VPlotAnasumHistograms();
        VPlotAnasumHistograms( string ifile, int ion = -1 );
        ~VPlotAnasumHistograms() {};
        
        void convert_derotated_RADECJ2000( double x = 0, double y = 0, double xerr = 0, double yerr = 0 );
        
        void drawPSF( TCanvas* c = 0, string iFile = 0, TH2D* h2 = 0, float iPSF = 0.1 );
        //   void fit_energy(double minE = -0.5, double maxE = 0.5 );
        bool openDataFile( string ifile );
        
        void help();                                                       // this will print all available functions
        
        
        // plotting functions
        void            plot_deadTimes();
        void            plot_mscPlots( int irebin = 2, double xmin = -2., double xmax = 4., string mscwfile = "" );
        void            plot_qualityHistograms( double iSourceStrength = 1., bool bUpper = true, int iMethod = 0 );
        TCanvas*        plot_theta2( double t2min = 0., double t2max = 0.3, int irbin = 5 );
        void            plot_theta2Correction();
        void            plot_UncorrelatedSkyPlots();
        void            plot_CorrelatedSkyPlots();
        void            plot_skyPlots_perRun( string iHistoName = "hmap_stereo_sig", double rmax = -1.,
                                              double zmin = -100., double zmax = -100., double rSource = 0.4,
                                              int iN = -1, unsigned int nstart = 0 );
        TCanvas*        plot_skyPlots( string iPlotMode = "colz", bool iSingleCanvases = false,
                                       double excess_min = -9999., double excess_max = -9999., double sig_min  = -999., double sig_max  = -999. );
        TCanvas*        plot_skyPlots_significance( bool iCorrelated = false, double rmax = 0., double zmin = -100., double zmax = -100., bool bPoster = false );
        TCanvas*        plot_on( bool iCorrelated = false, double rmax = -1. );
        TH1D*           plot_significanceDistributions( double rmax = 2.0, double rSource = 0.4,
                double xmin = -6.5, double xmax = 10.,
                TCanvas* c = 0, bool regioncodeflag = false );
        TCanvas*        plot_radec( int sPlot = 0, double rmax = -3., double zmin = -5., double zmax = -1000.,
                                    double xcenter = 0., double ycenter = 0.,
                                    bool bSlices = false, double fSliceXmin = -0.1, double fSliceXmax = 0.1, bool bProjX = true );
        vector<sSource> plot_catalogue( TCanvas* c, string iCatalogue = "Hipparcos_MAG8_1997.dat",
                                        double iMaxBrightness = 8.0, string iBand = "B", double iStarRadius = 0.30,
                                        int iColor = 1, int iLineStyle = 1, string hSkyMapName = "hmap_stereo_sig_REFLECTED",
                                        double iTextAngle = 45., int iMarkerStyle = 5 );
        void            plot_RBM_ring( double r, double iA, double t2, double iN );
        void            plot_reflectedRegions( TCanvas* iC, int i, int j, int iColor = 5 );
        void            plot_excludedRegions( TCanvas* c, int iLineColor = 6 );
        TH1D*           plot_triggerpattern( int ntel = 3, bool bPlot = true );
        TCanvas*	plot_cumulativeSignificance( bool doSqrtFit = true );
        
        void            setPlottingCorrelatedHistograms( bool iB = false )
        {
            fPlotCorrelated = iB;
        }
        void            setPlottingPSF_radius_deg( float iPSF_radius_deg = -1. )
        {
            fPlotDrawPSF = iPSF_radius_deg;
        }
        void            setPlottingUseHours( bool iB = true, int iZeroHours = 0 )
        {
            fPlotUseHours = iB;
            fPlotZeroHours = iZeroHours;
        }
        void            setDebugInfo( bool iB = false )
        {
            fDebug = iB;    // more debug output to screen
        }
        bool            setRunNumber( int iRun );                                      // select run for plotting
        
        ClassDef( VPlotAnasumHistograms, 16 );
};

#endif
