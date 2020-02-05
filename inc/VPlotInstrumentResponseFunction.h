//! VPlotInstrumentResponseFunction effective area plotter

#ifndef VPlotInstrumentResponseFunction_H
#define VPlotInstrumentResponseFunction_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TText.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include "VHistogramUtilities.h"
#include "VInstrumentResponseFunctionReader.h"
#include "VPlotUtilities.h"

using namespace std;

class VPlotInstrumentResponseFunction : public VPlotUtilities, public VHistogramUtilities
{
    private:
    
        bool   fDebug;
        string fName;
        
        // effective area data
        vector< VInstrumentResponseFunctionReader* > fData;
        
        // resolution fitter
        string fFunction_fitResolution;
        double fXmin_fitResolution;
        double fXmax_fitResolution;
        TF1*   fTF1_fitResolution;
        
        // last plotted graph
        TGraphErrors* gLastPlottedGraph;
        
        // canvas size
        int    fCanvasSize_X;
        int    fCanvasSize_Y;
        
        // general resolution plotter
        TCanvas* plotPSF( vector< double > i_Energy_TeV_lin, double iTheta2AxisMax = 0.05,
                          bool iCumulative = false, bool iPlotTheta2 = false );
        TCanvas* plotResolution( string iName, string iCanvasTitle, string iYTitle,
                                 double iYmin, double iYmax, string iResolutionTreeName,
                                 string iXaxis, TPad* iResolutionPad = 0, bool iLogY = false );
        TCanvas* plotResolution2D( unsigned int iDataSetID, string iName, string iCanvasTitle, string iYTitle,
                                   double iYmin, double iYmax, string iResolutionTreeName, string iXaxis, double iEnergySlice_TeV = -1. );
        double fLegendX1, fLegendX2, fLegendY1, fLegendY2;
        int fLegendNCol;
        TString fLegendHeader, fLegendOption;
        TLegend* makeLegend();
        
    public:
    
        VPlotInstrumentResponseFunction();
        ~VPlotInstrumentResponseFunction();
        
        bool         addInstrumentResponseData( int iDataID, string iFileList );
        bool         addInstrumentResponseData( string iInstrumentResponseFile, string iA_MC );
        bool         addInstrumentResponseData( string iInstrumentResponseFile,
                                                double iZe = 20., double iWoff = 0.5, int iAzBin = 0,
                                                double iIndex = 2.4, int iNoise = 200, string iA_MC = "A_MC",
                                                int iColor = -99, int iLineStyle = -99,
                                                int iMarkerStyle = -99, float iMarkerSize = -99.,
                                                float iEmin_linTeV = 1.e-6, float iEmax_linTeV = 1.e6,
                                                float iLineWidth = 2., string iPlotOption = "",
                                                float iEnergy_TeV_Eres_oneSided = -99., string iLegend = "" );
        bool         checkDataSetID( unsigned int iDataSetID );
        bool         fitResolution( TGraphErrors* g );
        TGraphErrors* getLastPlottedGraph()
        {
            return gLastPlottedGraph;
        }
        unsigned int getNumberOfDataSets()
        {
            return fData.size();
        }
        unsigned int getNumberOfGoodDataSets();
        TH1D*        getTheta2orThetaHistogram( unsigned int iDataSetID = 0, double i_Energy_TeV_lin = 1., bool iTheta2 = true );
        void         listDataSets();
        TCanvas*     plotAngularResolution( string iXaxis = "energy", string iProbabilityString = "68",
                                            double iMin = -1.e99, double iMax = -1.e99,
                                            TPad* iResolutionPad = 0, bool iLogY = false );
        TCanvas*     plotAngularResolution2D( unsigned int iDataSetID = 0, string iXaxis = "energy", string iProbabilityString = "68", double iEnergySlice_TeV = -1. );
        TCanvas*     plotAngularResolutionContaintment( string iXaxis = "energy" );
        TCanvas*     plotCoreResolution( string iXaxis = "energy", double iMax = -1.e99 );
        TCanvas*     plotCoreResolution2D( unsigned int iDataSetID = 0, string iXaxis = "energy" );
        void         plotCutEfficiency( unsigned int iDataSetID = 0 );
        TCanvas*     plotCutEfficiencyRatio( unsigned int iDataSetID = 0, unsigned int iCutID = 999,
                                             double iPlotMinimum = 0., double iPlotMaximum = 1.2 );
        TCanvas*     plotEffectiveArea( double iEffAreaMin_m2 = -1., double iEffAreaMax_m2 = -1., TPad* iEffAreaPad = 0, bool iSmooth = false );
        TCanvas*     plotEffectiveAreaRatio( unsigned int iDataSetID = 0, double ymin = 0., double ymax = 1.3 );
        TCanvas*     plotEnergyReconstructionBias( string iM = "mean", double ymin = 0., double ymax = 3., bool iLogBias = false, TCanvas* c = 0 );
        TCanvas*     plotEnergyReconstructionBias2D( unsigned int iDataSetID = 0, double iYmin = 0., double iYmax = 3., bool iNoDirectionCut = false );
        void         plotEnergyReconstructionLogBias( string iM = "mean", double ymin = -1., double ymax = 1. );
        void         plotEnergyReconstructionLogBias2D( unsigned int iDataSetID = 0, string iM = "mean", double ymin = -0.8, double ymax = 0.8 );
        void         plotEnergyReconstructionMatrix( unsigned int iDataSetID = 0, bool bFineBinning = true,
                string iMatrixType = "default",
                bool bInterPol = false, bool bPlotMedian = false );
        TCanvas*     plotEnergyResolution( double ymin = 0., double ymax = 0.4,
                                           TPad* iResolutionPad = 0,
                                           bool iSmoothEnergyResolution = false,
                                           string iMatrixType = "default" );  // from makeEffectiveArea
        //		TCanvas*     plotEnergyResolution( string iXaxis );                                                // from IRF
        TCanvas*     plotEnergyResolution2D( unsigned int iDataSetID = 0 );
        void         plotEnergySpectra( bool iWeighted = true, double iYMax = -1., int iRebin = 1 );
        void         plotMigrationMatrix( unsigned int iDataSetID = 0, bool bFineBinning = true,
                                          string iMatrixType = "default",
                                          bool bInterPol = false, bool bPlotMedian = false )
        {
            plotEnergyReconstructionMatrix( iDataSetID, bFineBinning, iMatrixType, bInterPol, bPlotMedian );
        }
        TCanvas*     plotTheta( double iTheta2AxisMax = 0.05, bool iCumulative = false );
        TCanvas*     plotTheta( vector< double > i_Energy_TeV_lin, double iTheta2AxisMax = 0.05, bool iCumulative = false );
        TCanvas*     plotTheta2( double iTheta2AxisMax = 0.05, bool iCumulative = false );
        TCanvas*     plotTheta2( vector< double > i_Energy_TeV_lin, double iTheta2AxisMax = 0.05, bool iCumulative = false );
        TCanvas*     plotWeightedRate();
        bool         removeInstrumentResponseData( int iDataSetID );
        void         resetInstrumentResponseData();
        void         setCanvasSize( int iX = 600, int iY = 600 )
        {
            fCanvasSize_X = iX;
            fCanvasSize_Y = iY;
        }
        void         setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        bool         setResolutionFitting( string iFitFunction = "", double iFitXmin = -1., double iFitXmax = 2. );
        void         setPlottingDefaults();
        void         setPlottingLogEnergyAxis( bool iB = true )
        {
            if( getPlottingAxis( "energy_Lin" ) )
            {
                getPlottingAxis( "energy_Lin" )->fLogAxis = iB;
            }
        }
        bool         write_fitResolutionFunction( string iOutName, string iName = "" );
        vector< VInstrumentResponseFunctionReader* >* getData()
        {
            return &fData;
        }
        void setLegendParameters( double x1 = -99, double y1 = -99, double x2 = -99, double y2 = -99, int N = 1, TString header = "", TString opt = "NDC" )
        {
            fLegendX1 = x1;
            fLegendY1 = y1;
            fLegendX2 = x2;
            fLegendY2 = y2;
            fLegendNCol = N;
            fLegendHeader = header;
            fLegendOption = opt;
        }
        
};

#endif
