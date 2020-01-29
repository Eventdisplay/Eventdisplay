//! VPlotWPPhysSensitivity plot CTA style sensitivities

#ifndef VPlotWPPhysSensitivity_H
#define VPlotWPPhysSensitivity_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

#include "VHistogramUtilities.h"
#include "VMathsandFunctions.h"
#include "VCTARequirements.h"
#include "VCTASensitivityRequirements.h"
#include "VPlotInstrumentResponseFunction.h"
#include "VPlotUtilities.h"
#include "VSensitivityCalculator.h"
#include "VSiteData.h"
#include "VUtilities.h"

class VPPUTValues
{
    private:
        vector< string > fSetName;
        vector< double > fPPUT;
        vector< double > fPPUTError;

        vector< double > fLowEPPUT;
        vector< double > fLowEPPUTError;

        vector< double > fMidEPPUT;
        vector< double > fMidEPPUTError;

        vector< double > fHighEPPUT;
        vector< double > fHighEPPUTError;

        vector< double > fnoLoEPPUT;
        vector< double > fnoLoEPPUTError;

        double getPPUT( TGraph *iG, bool iError = false, double ilogEMin = -20., double ilogEMax = 20. );

    public:
        VPPUTValues() {}
        ~VPPUTValues() {}
        void add( string iName, TGraph *iG );
        void print();
        void printLatexTable();
};


class VPlotWPPhysSensitivity : public VPlotUtilities
{
    private:
    
        VPlotInstrumentResponseFunction* fIRF;
        
        vector< VSiteData* > fData;
        
        double fMinEnergy_TeV;
        double fMaxEnergy_TeV;
        string fCrabSpectraFile;
        unsigned int fCrabSpectraID;
        
        string fRequirementsString;
        vector< VCTARequirements* > fPlotCTARequirements;
        
        bool fUseIntegratedSensitivityForOffAxisPlots;
        bool fNorthSouthComparision;

        string fCurrentInstrumentRootFile;
        vector< string > fCurrentInstrumentVector;
        
        // FOM variables
        double fSensitivityFOM;
        double fSensitivityFOM_error;

        // energy bias treatment
        float fMaximumAllowedEnergyBias;
        
        // projected sensitvity plots
        vector< double >  fProjectionEnergy_min_logTeV;
        vector< double >  fProjectionEnergy_max_logTeV;
        map< string, vector< TGraphAsymmErrors* > > fProjectionSensitivityvsCameraOffset;

        bool bPlotNoLegend;
        bool bPlotCrabLines;
        
        void    fillProjectedSensitivityPlot( unsigned int i, TGraphAsymmErrors* g );
        void    initialProjectedSensitivityPlots( bool iIncludeLowestEnergy = true );
        vector< TGraph* >  plotCurrentInstruments( TCanvas* c );
        double  getSensitivitySystematicUncertaintiesFactor( double );
        void    plotSensitivitySystematicUncertainties( TCanvas* c, TGraphAsymmErrors* g );
        bool    plotLegend( TCanvas* c = 0, bool iDown = false, bool iLeft = false, bool iAddFirst = true );
        
    public:
    
        VPlotWPPhysSensitivity();
        ~VPlotWPPhysSensitivity();
        
        bool addDataSet( VSiteData* iData );
        bool addDataSet( string iAnalysis, string iSubArray = "E", double iObservationTime_s = 180000., double iOffset_deg = 0.0,
                         string iLegend = "", int iColor = 1, int iLineStyle = 1, int iFillStyle = 3001 );
        bool addDataSets( string iDataSettxtFile, string iDirectionString );
        double getSensitivityFOM()
        {
            return fSensitivityFOM;
        }
        double getSensitivityFOM_error()
        {
            return fSensitivityFOM_error;
        }
        vector< VSiteData* > getData()
        {
            return fData;
        }
        bool plotIRF( string iPrint = "",
                      double iEffAreaMin = 50., double iEffAreaMax = 5.e7,
                      double iAngularResolutionMin = 0.005, double iAngularResolutionMax = 0.30,
                      double iEnergyResolutionMin = 0., double iEnergyResolutionMax = 0.5,
                      TPad* iEffAreaPad = 0, TPad* iAngResPad = 0, TPad* iEResPad = 0, bool iPlotEnergyBias = true,
                      bool iLogAngRes = true );
        TCanvas* plotProjectedSensitivities( TCanvas*, double iMaxOffset, int iColor = -1 );
        bool plotSensitivity( string iPrint = "",
                              double iMinSensitivity = 4.e-14, double iMaxSensitivity = 2.5e-10,
                              string iUnit = "ENERGY",
                              TPad* iSensitivityPad = 0, TPad* iBckPad = 0 );
        bool plotSensitivityRatio( string iPrint,
                                   double ymin = 0.01, double ymax = 2.,
                                   unsigned int iRatioSelector = 0, TPad* iSensRatio = 0,
                                   unsigned int iRatioGraphCounter = 0 );
        void reset();
        void setCurrentInstrumentFile( string iCurrentInstrumentRootFile = "CurrentInstruments.root" )
        {
            fCurrentInstrumentRootFile = iCurrentInstrumentRootFile;
        }
        void setCurrentInstrumentPlotVector();
        void setCurrentInstrumentPlotVector( vector< string > iV )
        {
            fCurrentInstrumentVector = iV;
        }
        void setCrabSpectraFile( string iFile = "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat",
                                 unsigned int iSpectraID = 5 )
        {
            fCrabSpectraFile = iFile;
            fCrabSpectraID = iSpectraID;
        }
        void setEnergyRange_Lin_TeV( double iMinEnergy_TeV = 0.01, double iMaxEnergy_TeV = 200. )
        {
            fMinEnergy_TeV = iMinEnergy_TeV;
            fMaxEnergy_TeV = iMaxEnergy_TeV;
        }
        void setMaximumAllowedEnergyBias( float emax_bias = -1. )
        {
            fMaximumAllowedEnergyBias = emax_bias;
        }
        void setNorthSouthComparision( bool iNS = false );
        bool setPlotCTARequirements( string iRequirements = "", 
                float iRequirementsScalingFactor = 1., double iRequirementsLineWidth = 1.,
                bool iRequirementsSystematic = false );
        void setPlotNoLegend( bool iPlotNoLegend = false )
        {
             bPlotNoLegend = iPlotNoLegend;
        }
        void setPlotCrabLines( bool iPlot = true )
        {
            bPlotCrabLines = iPlot;
        }
        void setUseIntegratedSensitivityForOffAxisPlots( bool iUseIntegratedSensitivityForOffAxisPlots = false )
        {
            fUseIntegratedSensitivityForOffAxisPlots = iUseIntegratedSensitivityForOffAxisPlots;
        }
};

#endif
