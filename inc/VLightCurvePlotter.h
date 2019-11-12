//! VLightCurvePlotter plot light curves
#ifndef VLightCurvePlotter_H
#define VLightCurvePlotter_H

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"

#include "VFluxAndLightCurveUtilities.h"
#include "VFluxDataPoint.h"
#include "VOrbitalPhaseData.h"
#include "VPlotUtilities.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

class VLightCurvePlotter : public VPlotUtilities
{
    private:
    
        bool fDebug;
        
        string fName;            // name of data set
        double fMinEnergy_TeV;   // max and min energy of flux calculation (used for axis defintion)
        double fMaxEnergy_TeV;   // max and min energy of flux calculation (used for axis defintion)
        
        // central data object
        vector< VFluxDataPoint > fFluxDataVector;
        TGraphAsymmErrors* fLightCurveGraph;
        
        // orbital phase data element
        VOrbitalPhaseData fOrbitalPhaseData;
        
        // light-curve plotting
        TCanvas* fCanvasLightCurve;
        double   fPlotting_TimeAxis_min;
        double   fPlotting_TimeAxis_max;
        double   fPlottingMJD_tolerance_dayfraction;   //!< add this fractional amount to the MJD axis (for automatic determination of MJD axis range)
        double   fPlotting_MJD_offset_days;            //!< offset in days used for plotting
        bool     fLightCurveTimeAxis_is_OrbitalPhase;
        bool     fIgnoreUpperLimits;                   //!< don't plot any upper limits
        bool     fPlot_observing_intervals_as_errors; //!< plot observing intervals as x-errors
        
        double   fPlotting_Flux_min;
        double   fPlotting_Flux_max;
        string   fPlotting_FluxAxis_title;
        
        // private plotting functions
        TGraphAsymmErrors* plotFluxes_vs_Variable( string iVariable, string iAxisTitle );
        
        // private utility functions
        void               reset();
        
    public:
    
        VLightCurvePlotter();
        VLightCurvePlotter( vector< VFluxDataPoint > iDataVector );
        ~VLightCurvePlotter() {}
        
        // getters
        string             getLightCurveFluxAxisTitle();
        double             getLightCurveFluxAxisRange_Min()
        {
            return fPlotting_Flux_min;
        }
        double             getLightCurveFluxAxisRange_Max()
        {
            return fPlotting_Flux_max;
        }
        TGraphAsymmErrors* getLightCurveGraph()
        {
            return fLightCurveGraph;
        }
        
        // plotters
        TCanvas*           plotLightCurve( TCanvas* iCanvasLightCurve = 0, string iCanvasName = "cL", string iPlottingOption = "p", double iMaxMJDError = -1. );
        TGraphAsymmErrors* plotFluxes_vs_elevation();
        TGraphAsymmErrors* plotFluxes_vs_pedvars();
        TGraphAsymmErrors* plotFluxes_vs_wobbleOffset();
        TCanvas* plotPhaseDistribution( TCanvas* iCanvasPhaseDist = 0, string iCanvasName = "cPD", string iFluxState = "", int iColor = 1 );
        
        // setters
        void setDataVector( vector< VFluxDataPoint > iDataVector );
        void setDebug( bool iDebug = true )
        {
            fDebug = iDebug;
        }
        void  setIgnoreUpperLimits( bool iIgnore = true )
        {
            fIgnoreUpperLimits = iIgnore;
        }
        void setLightCurveFluxAxis( double iYmin = -9.e10, double iYmax = -9.e10, string iAxisTitle = "not_set" );
        void setLightCurveTimeAxis( double iXmin = -99., double iXmax = -99., double iPlottingTimeAxis_tolerance_dayfraction = 0.1, double iPlotting_MJD_offset_days = 0. );
        void setLightCurveTimeAxis_to_OrbitalPhase( bool iOrbitalPhase = false );
        void setPlotObservingIntervals( bool iPlot_observing_intervals_as_errors = true )
        {
            fPlot_observing_intervals_as_errors = iPlot_observing_intervals_as_errors;
        }
};

#endif
