//! VAtmosphereSoundings read and analyse sounding data

#ifndef VAtmosphereSoundings_H
#define VAtmosphereSoundings_H

#include "VAstronometry.h"
#include "VAtmosphereSoundingData.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#define VMAXNUMBEROFSOUNDINGPOINTS 1000

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


class VAtmosphereSoundings
{
    private:
    
        bool             bDebug;
        
        vector< string > fListTextFile;
        
        // tree writing and reading
        TFile* fDataFile;
        TTree* fDataTree;
        
        double MJD;
        int Year;
        int Month;
        int Day;
        double Hour;
        unsigned int nPoints;
        double Pressure_Pa[VMAXNUMBEROFSOUNDINGPOINTS];
        double Height_m[VMAXNUMBEROFSOUNDINGPOINTS];
        double Density_gcm3[VMAXNUMBEROFSOUNDINGPOINTS];
        double Thickness_gcm2[VMAXNUMBEROFSOUNDINGPOINTS];
        double Temperature_K[VMAXNUMBEROFSOUNDINGPOINTS];
        double DewPoint_K[VMAXNUMBEROFSOUNDINGPOINTS];
        double RelativeHumidity[VMAXNUMBEROFSOUNDINGPOINTS];
        double VaporMassDensity_gm3[VMAXNUMBEROFSOUNDINGPOINTS];
        double MixingRatio_gkg[VMAXNUMBEROFSOUNDINGPOINTS];
        double WindDirection_deg[VMAXNUMBEROFSOUNDINGPOINTS];
        double WindSpeed_ms[VMAXNUMBEROFSOUNDINGPOINTS];
        double IndexofRefraction[VMAXNUMBEROFSOUNDINGPOINTS]; //new
        
        string           fTXTSearch_DataString;
        
        // data containers
        vector< VAtmosphereSoundingData* > fData;                         // profile from soundings data
        vector< VAtmosphereSoundingData* > fDataInterpol;
        vector< VAtmosphereSoundingData* > fDataCORSIKAMODTRAN;           // CORSIKA/MODTRAN profile
        vector< VAtmosphereSoundingData* > fDataUserProfile;              // user created profile
        vector< VAtmosphereSoundingData* > fAverageProfile;              // user created profile, henrikes method
        
        vector <double> fHeights;
        vector <double> fModtranHeights;
        
        // periods to calculate average over
        map< unsigned int, string > fPlottingPeriodFiles;
        map< unsigned int, vector< unsigned int > > fPlottingPeriodDates;
        
        // plotting
        bool   fPlotRelativePlots;
        bool   fBoolColorChange;
        vector< TLegend* > fPlottingLegend;
        bool   fPlottingLegendDraw;
        vector< TCanvas* > fCanvasProfile;
        vector< TCanvas* > fCanvas2D;
        string fPlottingPeriod;
        double fPlottingHeight_min;
        double fPlottingHeight_max;
        
        // observatory
        double fObservatoryLatitude;
        double fObservatoryHeight_km;
        
        unsigned int checkPlottingPeriodIdentifier( unsigned int );
        void   fillAtmosphericDensity();
        void   fillAtmosphericDensity( VAtmosphereSoundingData* );
        void   fillAtmosphericThickness();
        void   fillAtmosphericThickness( VAtmosphereSoundingData* );
        void   fillAtmosphericPressure();
        void   fillAtmosphericPressure( VAtmosphereSoundingData* );
        void   fillIndexofRefraction();
        void   fillO2();
        void   fillO3();
        void   fillWaterVaporDensity();
        void   fillWaterVaporDensity( VAtmosphereSoundingData* );
        double getInterpolation( double h, VAtmosphereSoundingData* iData, string iType );
        int    getMonth( string );
        unsigned int getHistogramIdentifier( unsigned int );
        vector< double > getDataVectorForUserAtmosphere( double iHeightMaxData, VAtmosphereSoundingData* iDataMonteCarlo, string iType );
        double getWaterVaporDensity( double T, double RH );
        double getWaterVaporMassDensity( double ATEMP );
        Color_t getSeasonColor( int iMonth );
        Color_t getSeasonColor( int counter, int atmo );
        TCanvas* plotCORSIKA( TCanvas* c, int iPlotID, vector< VAtmosphereSoundingData* > iData, double iHeightMin = 0., double iHeightMax = 120. );
        void   plotProfiles( unsigned int iYearStart, unsigned int iMonthStart, unsigned int iYearStop, unsigned int iMonthStop, bool b2D = false,
                             string iPlotOption = "", bool bSames = false );
        bool   readPlottingPeriodsFromTextFile( string );
        bool   readRootFile( unsigned npoints_min = 0 );
        
    public:
    
        VAtmosphereSoundings();
        VAtmosphereSoundings( string iRootFile, unsigned int npoints_min = 0 );
        ~VAtmosphereSoundings() {}
        bool     add_user_Atmosphere( unsigned int iIndexCORSIKAMODTRAN, double iHeightMaxData, string iName = "" );
        bool     readSoundingsFromTextFile( string iFileList );
        bool     readGDASFromTextFile( string iFileList );
        double   getAmosphericVaporPressure( double T );
        double   getDewPoint( double temperature, double relativeHumidity, int iMethod = 1 );
        void     list_datasets();
        void     list_datasets_CORSIKAMODTRAN();
        
        TCanvas* plotCORSIKA_Density_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 0, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_DewPoint_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 5, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_IndexofRefraction_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 1, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_Ozone_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 3, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_Pressure_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 6, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_RelativeHumidity_vs_Height( TCanvas* c, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 4, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_Temperature_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 2, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        TCanvas* plotCORSIKA_Thickness_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 7, fDataCORSIKAMODTRAN, iHeightMin, iHeightMax );
        }
        
        TCanvas* plotUserAtmosphere_Density_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 0, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_DewPoint_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 5, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_IndexofRefraction_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 1, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_Ozone_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 3, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_Pressure_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 6, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_RelativeHumidity_vs_Height( TCanvas* c, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 4, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_Temperature_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 2, fDataUserProfile, iHeightMin, iHeightMax );
        }
        TCanvas* plotUserAtmosphere_Thickness_vs_Height( TCanvas* c = 0, double iHeightMin = 0., double iHeightMax = 120. )
        {
            return plotCORSIKA( c, 7, fDataUserProfile, iHeightMin, iHeightMax );
        }
        
        void     plot2DProfiles( unsigned int iYearStart = 1980, unsigned int iMonthStart = 1,
                                 unsigned int iYearStop = 2020, unsigned int iMonthStop = 12 );
        void     plotAverages( unsigned int iYearStart = 1980, unsigned int iMonthStart = 1,
                               unsigned int iYearStop = 2020, unsigned int iMonthStop = 12,
                               string iPlotOption = "", bool iSames = false );
        void     plotAttributes_ColorChange( bool iB = true )
        {
            fBoolColorChange = iB;
        }
        void     plotAttributes_PlotLegend( bool iB = true )
        {
            fPlottingLegendDraw = iB;
        }
        int     read_CORSIKA_Atmosphere( string iFile, string iName = "", int iColor = 2, int iLineStyle = 1 );
        int     read_MODTRAN_Atmosphere( string iFile, string iName = "", int iColor = 2, int iLineStyle = 1 );
        bool     readSoundingsFromRootFile( string iRootFile, unsigned int npoints_min = 0 );
        void     setGeographicPosition( double iLatitude = 31.675, double iObsHeight_km = 1.27 )
        {
            fObservatoryLatitude = iLatitude;
            fObservatoryHeight_km = iObsHeight_km;
        }
        void     setPlottingPeriod( string iPeriod = "monthly" );
        void     setPlottingRangeHeight( double iHeightMin = 0., double iHeightMax = 30. )
        {
            fPlottingHeight_min = iHeightMin;    // in [km]
            fPlottingHeight_max = iHeightMax;
        }
        void     setPlottingRelativePlots( bool iB = false )
        {
            fPlotRelativePlots = iB;
        }
        bool     write_CORSIKA_UserProfile( unsigned int iMODTRANIndex, unsigned int atmprofmodel, string iName = "user profile" );
        bool     write_MODTRAN_UserProfile( unsigned int iIndexUserData, unsigned int defaultModel = 6, bool iWriteDewPoint = false );
        bool     writeRootFile( string iFile );
        
        int push_average_atmosphere( string name, vector<int>* years, vector<int>* months, vector<int>* days , vector<int>* hours  , vector<double>* mjds, unsigned int nMinPoints, int nMinFlights );
        bool isDateInRange( VAtmosphereSoundingData* Data, vector<int>* years, vector<int>* months, vector<int>* days, vector<int>* hours, vector<double>* mjds , unsigned int nMinPoints );
        bool isDateInRange( VAtmosphereSoundingData* Data, vector<double>* mjds, unsigned int nMinPoints );
        bool isDateInRange( VAtmosphereSoundingData* Data, double minMJD, double max_MJD , unsigned int nMinPoints );
        VAtmosphereSoundingData* makeDefaultAtmosphere( string season, string name = "noname", string opt = "", int year = 0 );
        VAtmosphereSoundingData* makeDefaultWinterAtmosphere( string name = "winter", string opt = "", int year = 0 ) ;
        VAtmosphereSoundingData* makeDefaultSummerAtmosphere( string name = "summer", string opt = "", int year = 0 ) ;
        VAtmosphereSoundingData* makeMeanMonthlyAtmosphere( int month, string name, string opt, int yearMin, int yearMax );
        VAtmosphereSoundingData* makeOneFlightAtmosphere( int year, int month, int day, int hour, string name, string opt );
        VAtmosphereSoundingData* makeMeanAtmosphereMJD( double minMJD, double maxMJD, string name, string opt );
        
        double interpolate( vector<double> raw, vector<double> raw_heights, vector<double>& result, string opt, double h );
        double safe_eval( TGraph* g, double h, string opt );
        
        VAtmosphereSoundingData* make_interpolated_atmosphere( VAtmosphereSoundingData* RawData );
        void make_interpolated_atmospheres();
        void setModtranHeights();
        vector<double> getModtranHeights()
        {
            return fModtranHeights;
        }
        vector<double> getHeights()
        {
            return fHeights;
        }
        void setHeights( vector<double> heights )
        {
            fHeights = heights;
        }
        //void setHeights( double min = 0.0, double max = 40000.0, double step = 100.0 )
        void setHeights( double min = 0.0, double max = 40000.0, double step = 1000.0 )
        {
            fHeights.clear();
            for( double h = min; h <= max; h += step )
            {
                fHeights.push_back( h );
            }
        }
        
        VAtmosphereSoundingData* getData( unsigned int i )
        {
            return i >= fData.size()		? 	0	:	fData.at( i );
        }
        VAtmosphereSoundingData* getDataInterpol( unsigned int i )
        {
            return i >= fDataInterpol.size()	? 	0	:	fDataInterpol.at( i );
        }
        VAtmosphereSoundingData* getDataMODTRAN( unsigned int i )
        {
            return i >= fDataCORSIKAMODTRAN.size()	? 	0	:	fDataCORSIKAMODTRAN.at( i );
        }
        VAtmosphereSoundingData* getDataUserProfile( unsigned int i )
        {
            return i >= fDataUserProfile.size()	? 	0	:	fDataUserProfile.at( i );
        }
        VAtmosphereSoundingData* getAverageProfile( unsigned int i )
        {
            return i >= fAverageProfile.size()	? 	0	:	fAverageProfile.at( i );
        }
        
        bool write_2C1( unsigned int iIndexAverageData, string filename, double max_height );
        
        TGraph* getResidualGraph( TGraph* data, TGraph* model , int color = 2 ) ;
        TCanvas* plot_season( vector<VAtmosphereSoundingData*> v, TString season_name, string value, TString outfileprefix ) ;
        TCanvas* plot_season( int year_start, int month_start, int day_start, int year_end, int month_end , int day_end, string value, int bWriteCorsika = -1 );
        TCanvas* plot_monthly( vector< int > year, vector< int > month, vector< int > day, double intervall_days, int offset_months, string value );
        
        int readEpochsAndAtmospheres( TString iDstart, double iMonthLength_days,
                                      string iEpochFile = "$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/VERITAS.Epochs.runparameter" );
                                      
                                      
};

#endif
