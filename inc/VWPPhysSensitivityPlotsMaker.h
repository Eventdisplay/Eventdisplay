//! VWPPhysSensitivityPlotsMaker class to produce a tex files for CTA/VTS sensitivties

#ifndef VWPPhysSensitivityPlotsMaker_H
#define VWPPhysSensitivityPlotsMaker_H

#include "VPlotWPPhysSensitivity.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class VWPPhysSensitivityPlotsMaker : public VPlotUtilities
{
    private:
    
        double           fObservingTime_s;
        vector< string > fListofDataSets;
        vector< string > fListOfArrays;
        vector< double > fOffAxisAngle;
        
        double fSensitivity_min;
        double fSensitivity_max;
        double fSensitivityRatio_min;
        double fSensitivityRatio_max;
        string fSensitivity_Unit;
        double fMinEnergy_TeV;
        double fMaxEnergy_TeV;
        double fAngularResolution_min;
        double fAngularResolution_max;
        double fEnergyResolution_min;
        double fEnergyResolution_max;
        double fEffArea_min;
        double fEffArea_max;
        bool   fPlotAngResLogY;
        double fBackgroundRatio_min;
        double fBackgroundRatio_max;
        double fEffAreaRatio_min;
        double fEffAreaRatio_max;
        double fAngResRatio_min;
        double fAngResRatio_max;
        
        string fPrintingOptions;
        bool bPlotNoLegend;
        bool bPlotCrabLines;
        
        string fPlotCTARequirementsString;
        float  fRequirementsScalingFactor;
        int    fRequirementsLineWidth;
        bool   fRequirementsSystematics;

        vector< string > fCurrentInstrumentVector;
	    bool fPlotCurrentInstrumentVectorLabel;

        float  fMaximumAllowedEnergyBias;
 
        string fPlotAllInOneCanvasStyle;
        // projected off-axis sensitivity
        TCanvas* fPlotProjectedSensitivity;
        // pads for plotAllInOneCanvas()
        TCanvas* fPlotAllInOneCanvas;
        TPad* fSensitivityPad;
        TPad* fSensitivityTitlePad;
        TPad* fSensitivityRatioPad;
        TPad* fEffAreaPad;
        TPad* fEffAreaRatioPad;
        TPad* fBckRatesPad;
        TPad* fBckRatesPadRatio;
        TPad* fERes;
        TPad* fAngRes;
        TPad* fAngResRatio;
        
    public:
    
        VWPPhysSensitivityPlotsMaker();
        ~VWPPhysSensitivityPlotsMaker() {}
        
        void compareDataSets( string iDataSetFile, string iDirectionString = "", bool IntegratedSensitivityForOffAxisPlots = false,
                              unsigned int iRatioCounter = 0, string iTitlteText = "" );
        void compareOffAxisSensitivities( string iSubArray, vector< string > iDataSet );
        void compareOffAxisSensitivities( string iSubArray = "", string iDataSet = "" );
        TCanvas* getAllinOneCanvas()
        {
            return fPlotAllInOneCanvas;
        }
        TCanvas* getProjectedSensitivityCanvas()
        {
            return fPlotProjectedSensitivity;
        }
        void plotAllInOneCanvas( bool iCanvasBatch = false );
        void plotRatioPlot( TPad *corg, TPad *ratio, double ymin = 0., double ymax = 2. );
        void printPlotCTARequirementsIDs();
        void resetVectors();
        void setAxisUnits( string iObservationTime );   // set the correct y-axis scale for 50h, 5h, and 0.5h
        void setAxisUnits( double iMinSensitivity = 4.e-14, double iMaxSensitivity = 1.5e-10, string iUnit = "ENERGY" );
        void setCurrentInstrumentPlotVector( vector< string > iV, bool iPlotLabels = false )
        {
            fCurrentInstrumentVector = iV;
	        fPlotCurrentInstrumentVectorLabel = iPlotLabels;
        }
        void setEffectiveAreaLimits( double iEffArea_min = 5.e2, double iEffArea_max = 7.e6 )
        {
            fEffArea_min = iEffArea_min;
            fEffArea_max = iEffArea_max;
        }
        void setPlotCrabLines( bool iPlot = true )
        {
            bPlotCrabLines = iPlot;
        }
        void setPlotAllInOneCanvasStyle( string iStyle = "AllRatios" )
        {
             fPlotAllInOneCanvasStyle = iStyle;
        }
        void setResolutionLimits( double iAngularResolutionMax = 1.10, double iEnergyResolutionMax = 0.3,
                                  bool iAngresLogY = true,
                                  double iAngularResolutionMin = 0.8e-2, double iEnergyResolutionMin = 0. )
        {
            fAngularResolution_max = iAngularResolutionMax;
            fAngularResolution_min = iAngularResolutionMin;
            fEnergyResolution_max = iEnergyResolutionMax;
            fEnergyResolution_min = iEnergyResolutionMin;
            fPlotAngResLogY       = iAngresLogY;
        }
        void setSensitivityRatioLimits( double iRatio_min = 0.1, double iRatio_max = 1.5 )
        {
            fSensitivityRatio_min = iRatio_min;
            fSensitivityRatio_max = iRatio_max;
        }
        void setIRFRatioLimits( double iRatio_bck_min = 0., double iRatio_bck_max = 2.,
                                double iRatio_effarea_min = 0., double iRatio_effarea_max = 2.,
                                double iRatio_angres_min = 0., double iRatio_angres_max = 2. )
        {
            fBackgroundRatio_min = iRatio_bck_min;
            fBackgroundRatio_max = iRatio_bck_max;
            fEffAreaRatio_min = iRatio_effarea_min;
            fEffAreaRatio_max = iRatio_effarea_max;
            fAngResRatio_min = iRatio_angres_min;
            fAngResRatio_max = iRatio_angres_max;
        }

        void setEnergyRange_Lin_TeV( double iMinEnergy_TeV = 0.01, double iMaxEnergy_TeV = 200. )
        {
            fMinEnergy_TeV = iMinEnergy_TeV;
            fMaxEnergy_TeV = iMaxEnergy_TeV;
        }
        void setLSTSettings();
        void setMaximumAllowedEnergyBias( float emax_bias = -1. )
        {
            fMaximumAllowedEnergyBias = emax_bias;
        }
        void setObservingTime( double i_s = 180000. )
        {
            fObservingTime_s = i_s;
        }
        void setPrintingOptions( string iPrint = "" )
        {
            fPrintingOptions = iPrint;
        }
        void setPlotRequirements( string iRequirement = "", 
                float iRequirementsScalingFactor = 1., 
                double iLineWidth = 1.  );
        void setPlotRequirementsSystematics( bool iPlotSystematics = false )
        {
            fRequirementsSystematics = iPlotSystematics;
        }
        void setOffAxisAngle( vector< double > iA )
        {
            fOffAxisAngle = iA;
        }
        void setPlotNoLegend( bool iPlotNoLegend = false )
        {
            bPlotNoLegend = iPlotNoLegend;
        }
        bool writeTexFileBody( string iTexFile, string iTexFileTitle = "" );
};

#endif
