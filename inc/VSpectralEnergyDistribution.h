//! VSpectralEnergyDistribution

#ifndef VSpectralEnergyDistribution_H
#define VSpectralEnergyDistribution_H

#include "TArrow.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "VFluxAndLightCurveUtilities.h"

using namespace std;

/*
   effective wavelengths and zero-points for photometric bands
*/
struct sPhotoMetricBand
{
    string fBand;
    double fEffectiveWavelength_micron;
    double fF0_CIT_Jy;
    double fF0_UKIRT_Jy;
};


/*
   source specific galactic extinction correction
*/
struct sGalacticExtinction
{
    string fBand;
    double fCorrection;
};


/*

 */
struct sPhotonFlux
{
    string name;
    double MJD_min;
    double MJD_max;
    vector< double > energy_Hz;
    vector< double > energy_Hz_min;
    vector< double > energy_Hz_max;
    vector< double > energy_eV;
    vector< double > energy_eV_min;
    vector< double > energy_eV_max;
    vector< double > flux_ergscms;
    vector< double > flux_error_up_ergscms;
    vector< double > flux_error_down_ergscms;
    int Marker;
    int Color;
};

class VSpectralEnergyDistribution
{
    private:
    
        bool fDebug;
        
        string fName;                             // project name
        
        // constants
        vector< sPhotoMetricBand >      fPhotoMetricBand;
        vector< sGalacticExtinction >   fGalacticExtinction;
        
        // data sets
        // [dataset][data points]
        vector< vector< sPhotonFlux > > fSpectralFlux;
        
        // plot/print variables
        double fMJDMin;
        double fMJDMax;
        
        double fPlotting_EnergyRange_min_Hz;
        double fPlotting_EnergyRange_max_Hz;
        double fPlotting_FluxRange_min;
        double fPlotting_FluxRange_max;
        
        TCanvas* plotCanvas( int canvas_x = 600, int canvas_y = 600 );
        
    public:
    
        VSpectralEnergyDistribution( string name = "SED" );
        virtual ~VSpectralEnergyDistribution() {}
        
        TCanvas* plot( TCanvas* c = 0, int bLegend = false, int canvas_x = 900, int canvas_y = 600, bool bErrorX = true, bool bPlotName = false );
        double   getEffectiveWavelength( string iband, string iunit );
        double   getGalacticExtinctionCorrection( string iband );
        double   getFluxfromMagnitude( double magnitude, string band, string system = "CIT" );
        TGraph* plotModel( TCanvas* c, string ifile, int icolor = 1, int ilinestyle = 1, int ilinewidth = 2, bool isJyHz = false );
        TCanvas* plotPowerLaw( TCanvas* c, string iName, double iEMin_TeV, double iEMax_TeV,
                               double iNorm, double iGamma, double iNormEnergy = 1.,
                               bool bPlotButterfly = false, double iNormError = 0., double iGammaError = 0.,
                               int iLineColor = 1, int iLineStyle = 1 );
        void printASCII();
        bool readPhotoMetricBands( string ifile = "$OBS_EVNDISP_AUX_DIR/AstroData/Multiwavelengthdata/photometricBands.dat", bool iPrint = true );
        bool readGalacticExtinction( string ifile, bool iPrint = true );
        bool readDataFile( string name, string txtfile, double MJD_min, double MJD_max, bool bPrint = false, int imarker = 20, int icolor = 1 );
        TGraphErrors* readOpticalData( string name, string txtfile, string band, bool bPrint = false, int imarker = 20, int icolor = 1, bool bAverage = false, double iPlotMagnitudeMultiplier = 1., bool bCorrection = false, string icorfile = "" );
        bool readSED( string iname );
        bool readSwiftData( string name, string txtfile, double MJD_min, double MJD_max, bool bPrint = false, int imarker = 20, int icolor = 1 );
        bool readTeVEvndispData( string name, string txtfile, bool bPrint = false, int imarker = 20, int icolor = 1 );
        bool readFermiData( string name, string txtfile, double MJD_min, double MJD_max, bool bPrint = false, int imarker = 20, int icolor = 1, int iFormat = 0 );
        bool readXMMData( string name, string txtfile, double MJD_min, double MJD_max, bool bPrint = false, int imarker = 20, int icolor = 1, bool iModel = false, int iFormat = 0 );
        void setPlottingEnergyRange_Hz( double iMin = 1.e9, double iMax = 1.e28 )
        {
            fPlotting_EnergyRange_min_Hz = iMin;
            fPlotting_EnergyRange_max_Hz = iMax;
        }
        void setPlottingFluxRange( double iMin = 1.e-14, double iMax = 1.e-9 )
        {
            fPlotting_FluxRange_min = iMin;
            fPlotting_FluxRange_max = iMax;
        }
        void setTimeRange( double iMJDmin = 0., double iMJDmax = 1.e14 );
        
        ClassDef( VSpectralEnergyDistribution, 3 );
};
#endif
