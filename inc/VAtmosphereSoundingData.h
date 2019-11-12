//! VAtmosphereSoundings read and analyse sounding data

#ifndef VAtmosphereSoundingData_H
#define VAtmosphereSoundingData_H

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#define VMAXNUMBEROFSOUNDINGPOINTS 1000


using namespace std;

class VAtmosphereSoundingData
{
    public:
    
        double MJD;
        int Year;
        int Month;
        int Day;
        double Hour;
        string Name;
        vector< double > fHeight_m;
        vector< double > fPressure_Pa;
        vector< double > fDensity_gcm3;
        vector< double > fThickness_gcm2;
        vector< double > fTemperature_K;
        vector< double > fDewPoint_K;
        vector< double > fRelativeHumidity;
        vector< double > fVaporMassDensity_gm3;
        vector< double > fMixingRatio_gkg;
        vector< double > fWindDirection_deg;
        vector< double > fWindSpeed_ms;
        vector< double > fIndexofRefraction;
        vector< double > fO2_cmkm;
        vector< double > fO3_cmkm;
        // vectors for errors
        vector< double > fRMS_Pressure_Pa;
        vector< double > fRMS_Density_gcm3;
        vector< double > fRMS_Temperature_K;
        vector< double > fRMS_RelHum;
        vector< double > fRMS_IndexofRefraction;
        
        int PlotColor;
        int PlotMarker;
        float PlotMarkerSize;
        int PlotLineStyle;
        Width_t PlotLineWidth;
        
        VAtmosphereSoundingData();
        ~VAtmosphereSoundingData() {}
        void setdefaultvalues( unsigned int iN );
        
        void makeGraphScaledDensity( );
        void makeGraphPressure( );
        void makeGraphHumidity( );
        void makeGraphTemperature( );
        void makeGraphIndex( );
        void makeGraphThickness( );
        bool write_2C1( string filename, vector<double>* ModtranHeights, double max_height = 150e3 );
        bool write_CORSIKA_UserProfile( unsigned int atmprofmodel, string iName );
        
        TGraphErrors* getDensityGraph()
        {
            if( !fGraphScaledDensityHeight )
            {
                makeGraphScaledDensity();
            }
            return fGraphScaledDensityHeight;
        }
        TGraphErrors* getPressureGraph()
        {
            if( !fGraphPressureHeight )
            {
                makeGraphPressure();
                
            }
            return fGraphPressureHeight;
        }
        TGraphErrors* getHumidityGraph()
        {
            if( !fGraphHumidityHeight )
            {
                makeGraphHumidity();
            }
            return fGraphHumidityHeight;
        }
        TGraphErrors* getTemperatureGraph()
        {
            if( !fGraphTemperatureHeight )
            {
                makeGraphTemperature();
            }
            return fGraphTemperatureHeight;
        }
        TGraphErrors* getIndexGraph()
        {
            if( !fGraphIndexHeight )
            {
                makeGraphIndex();
            }
            return fGraphIndexHeight;
        }
        TGraphErrors* getThicknessGraph()
        {
            if( !fGraphThicknessHeight )
            {
                makeGraphThickness();
            }
            return fGraphThicknessHeight;
        }
        TGraphErrors* getGraph( string which )
        {
            if( which == "density" )
            {
                return getDensityGraph();
            }
            else if( which == "temperature" )
            {
                return getTemperatureGraph() ;
            }
            else if( which == "pressure" )
            {
                return getPressureGraph() ;
            }
            else if( which == "humidity" )
            {
                return getHumidityGraph() ;
            }
            else if( which == "index" )
            {
                return getIndexGraph() ;
            }
            else if( which == "thickness" )
            {
                return getThicknessGraph() ;
            }
            else
            {
                return 0;
            }
        }
        
        TGraphErrors* fGraphScaledDensityHeight;
        TGraphErrors* fGraphPressureHeight;
        TGraphErrors* fGraphHumidityHeight;
        TGraphErrors* fGraphTemperatureHeight;
        TGraphErrors* fGraphIndexHeight;
        TGraphErrors* fGraphThicknessHeight;
        
        void setColor( int color );
        int getColor( )
        {
            return PlotColor;
        }
        
};

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////



#endif
