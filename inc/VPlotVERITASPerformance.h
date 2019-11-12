//! VPlotVERITASPerformance  plot for VERITAS performance page

#ifndef VPlotVERITASPerformance_H
#define VPlotVERITASPerformance_H

#include <vector>

#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TText.h"

#include "VPlotInstrumentResponseFunction.h"
#include "VSensitivityCalculator.h"

using namespace std;

class VPlotVERITASPerformanceEpochData
{
    public:
    
        string fName;
        string fTitle;
        
        bool fPlotThisEpoch;
        int  fColor;
        int  fLineStyle;
        
        float fZe;
        int fNoiseLevel;
        
        float fGammaRayRate;
        float fBckRate;
        
        string fCutsIRFFile;
        
        VPlotVERITASPerformanceEpochData();
        ~VPlotVERITASPerformanceEpochData() {}
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class VPlotVERITASPerformance
{
    private:
    
        VPlotInstrumentResponseFunction* fPlotter;
        
        vector< VPlotVERITASPerformanceEpochData* > fEpochData;
        
        double fPlotEnergy_TeV_max;
        double fPlotEnergy_TeV_min;
        
        int    fOtherInstrumentColor;
        
        string fPlotPrintName;
        string fVersionNumber;
        
        TLegend* getLegend( bool iLeft = false, bool iBottom = false );
        void     printCanvas( TCanvas* c, string iTitle = "" );
        void     plotVersionNumber( TCanvas* c );
        TGraphAsymmErrors* readDifferentialSensitivity( string iFile );
        void     resetSettings();
        
    public:
    
        VPlotVERITASPerformance();
        ~VPlotVERITASPerformance() {}
        
        void plotAngularResolutionForDifferentEpochs( double iMaxTheta = 0.22 );
        void plotAngularResolutionInComparison( bool iLatPass8 = false );
        void plotDifferentialSensitivity( double iE_min = 0.05, double iE_max = 11. );
        void plotEffectiveArea( float iEnergy_TeV_Max = 50., bool iPlotOtherInstruments = false );
        void plotEnergyResolutionForDifferentEpochs( double iMax = 0.45 );
        void plotPlotsForPerformancePage( string iPrintName = "", string iVersionNumber = "" );
        void plotSensitivity_vs_time( int iPlotID = 0, bool iGuidingLines = true );
        
        
        void setPlotEnergyRange( double iPlotEnergy_TeV_min = 0.08, double iPlotEnergy_TeV_max = 25. )
        {
            fPlotEnergy_TeV_min = iPlotEnergy_TeV_min;
            fPlotEnergy_TeV_max = iPlotEnergy_TeV_max;
        }
        void setPlotOtherInstrumentColor( int iColor = 12 )
        {
            fOtherInstrumentColor = iColor;
        }
        void setVersionNumber( string iVersionText )
        {
            fVersionNumber = iVersionText;
        }
        void setPlotAllEpochs();
        void setPrintName( string iPrintName )
        {
            fPlotPrintName = iPrintName;
        }
        bool setDataSet( bool iV6Only = true, string iCut = "moderate" );
        
        
};

#endif
