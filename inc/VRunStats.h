//! VRunStats get basic run data from DB

#ifndef VRUNSTATS_H
#define VRUNSTATS_H

#include "VAstronometry.h"
#include "VGlobalRunParameter.h"
#include "VDB_Connection.h"

#include <TFile.h>
#include <TMath.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <TTree.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class VDQMRunInfo
{
    public:
    
        unsigned int runNumber;
        int fieldVersion;
        
        double StartMJD;
        string SourceID;
        string ObsMode;
        int NTel;
        unsigned int TelBitMask;
        
        double WobbleOff_estim_deg;
        double WobbleDir_estim_deg;
        double El_mean_deg;
        double Az_mean_deg;
        double L3_runTime_s;
        double L3_lifeTime_s;
        double L3_rate;
        double Rate_mean;
        double Rate_RMS;
        double Rate_Chi2;
        
        double L3_threshold;
        
        double L2_threshold[4];
        double L2_threshold_mV[4];
        double L2_threshold_Spread_mV[4];
        double L2_throughput[4];
        
        double GPS_runtime_min;
        double GPS_lifetime_min;
        
        double MoonEl_mean_deg;
        double MoonPhase_mean;
        double MoonSeparation_deg[4];
        
        double RA_mean[4];
        double RA_RMS[4];
        double Dec_mean[4];
        double Dec_RMS[4];
        
        double L1_rate[4];
        double Median_current_uA[4];
        double L2_rate[4];
        double MissingEvents[4];
        double Muons_N[4];
        double Muons_MeanCharge[4];
        double Muons_RMSCharge[4];
        double PMT_GainLaser[4];
        
        double SuppressedChannels_Num[4];
        
        double PedVar_mean[4];
        
        int LaserRun[4];
        
        double FIR0_mean;
        double FIR0_RMS;
        double FIR_mean[4];
        double FIR_RMS[4];
        
        VDQMRunInfo();
        ~VDQMRunInfo();
        
};

class VDBAnalysisComments
{
    public:
    
        int runNumber;
        unsigned int status;
        float usable_duration;
        string comment;
        
        VDBAnalysisComments();
        ~VDBAnalysisComments() {}
};

class VDBDefaultHV
{
    public:
        double startMJD;
        double stopMJD;
        unsigned int telescope;
        unsigned int channel;
        double default_HV;
        
        VDBDefaultHV();
        ~VDBDefaultHV();
};

class VDBHV
{
    public:
        double startMJD;
        double stopMJD;
        unsigned int telescope;
        unsigned int channel;
        double HV;
        
        VDBHV();
        ~VDBHV();
};

class VDBCameraStatus
{
    public:
    
        double MJD;
        unsigned int telescope;
        
        double temp1;
        double temp2;
        double hum1;
        double hum2;
        double light1;
        double light2;
        
        VDBCameraStatus();
        ~VDBCameraStatus();
};

class VDBFIRData
{
    public:
    
        double MJD;
        unsigned int tel_ID;
        float ambient_temp;
        float radiant_sky_temp;
        float radiant_sky_temp_cor;
        
        VDBFIRData();
        ~VDBFIRData() {}
};


class VDBWeatherData
{
    public:
    
        double MJD;
        double AirTemperature;
        double AirPressure;
        double AirHumidity;
        double AirWindSpeed;
        double AirWindGust;
        double AirWindDirection;
        
        VDBWeatherData();
        ~VDBWeatherData() {}
};

class VDBRunData
{
    public:
    
        unsigned int runNumber;
        unsigned int ConfigMask;
        string SourceID;
        double StartMJD;
        double StopMJD;
        string start_time;
        string stop_time;
        double Duration;
        double RA;
        double Dec;
        double offsetRA;
        double offsetDec;
        double offset_distance;
        double offset_angle;
        double GalLong1958;
        double GalLat1958;
        
        vector< VDBFIRData* > fDBFIRData;
        vector< VDBWeatherData* > fDBWeatherData;
        vector< VDBCameraStatus* > fCameraStatus;
        
        VDBRunData();
        ~VDBRunData() {};
};

class VDBSourceData
{
    public:
    
        string SourceID;
        double fDec;
        double fRA;
        
        VDBSourceData();
        ~VDBSourceData() {}
};

/////////////////////////////////////////

class VRunStats : public TObject, public VGlobalRunParameter
{
    private:
    
        bool fDebug;
        
        bool fReadfullHVEvndispData;
        
        int fMax_run_id;
        int fMin_run_id;
        
        string fStartDate;
        string fStopDate;
        vector< VDBRunData* > fRunData;
        map< string, VDBSourceData* > fSourceData;
        vector< VDBWeatherData* > fWeatherData;
        vector< VDBFIRData* > fFIRData;
        vector< VDBCameraStatus* > fCameraStatus;
        vector< VDBAnalysisComments* > fAnalysisComments;
        vector< VDBDefaultHV* > fdefaultHVEvndispData;
        vector< VDBHV* > fHVEvndispData;
        vector< VDQMRunInfo* > fDQMData;
        
        void getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip );
        
        bool readDBAnalysisComments();
        bool readDBHVEvndispData( TSQLServer* f_db );
        bool readDB_default_HVEvndispData( TSQLServer* f_db );
        bool readDBCameraStatus( TSQLServer* f_db );
        bool readDBFIRInfo( TSQLServer* f_db );
        bool readDBRun_IDs( TSQLServer* f_db );
        bool readDBRunInfo( TSQLServer* f_db );
        bool readDBSourceInfo( TSQLServer* f_db );
        bool readDBWeatherInfo( TSQLServer* f_db );
        
        bool getDBWeatherData( TSQLServer* f_db, string iTimeStamp, double& iWindSpeed, double& WindGust, double& iWindDir, double& iTemp, double& iHumidity );
        bool getDBSourceCoordinates( TSQLServer* f_db, string iSource, double& iEVNTargetDec, double& iEVNTargetRA );
        
    public:
    
        VRunStats();
        ~VRunStats() {}
        
        void print();
        void printDQM( unsigned int iRun );
        bool readFromDB();
        bool readDQMData( string iDQMFile );
        void setReadFullHVEvndispData( bool iHVEvndispData = true )
        {
            fReadfullHVEvndispData = iHVEvndispData;
        }
        void setTimeRange( string iS = "2009-01-01", string iE = "2009-01-30" );
        bool writeRootFile( string ifile );
        
        ClassDef( VRunStats, 3 );
};
#endif
