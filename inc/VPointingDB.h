//! VPointingDB get pointing from DB

#ifndef VPOINTINGDB_H
#define VPOINTINGDB_H

#include "PointingMonitor.h"
#include "VAstronometry.h"
#include "VGlobalRunParameter.h"
#include "VTrackingCorrections.h"
#include "VDB_Connection.h"

#include <TMath.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class VPointingDB : public VGlobalRunParameter
{
    private:
        double degrad;
        
        double fObsLongitude;
        double fObsLatitude;
        
        bool fStatus;
        
        bool fGoodVPM;
        
        unsigned int fTelID;
        unsigned int fRunNumber;
        
        int fMJDRunStart;                         // run start time
        double fTimeRunStart;
        int fMJDRunStopp;                         // run stop time
        double fTimeRunStopp;
        
        unsigned int fMJD;
        double fTime;
        
        unsigned int fEventStatus;
        float fTelAzimuth;
        float fTelElevation;
        float fTelExpectedAzimuth;
        float fTelExpectedElevation;
        
        string fDBSourceName;
        float fDBTargetDec;
        float fDBTargetRA;
        float fDBWobbleNorth;
        float fDBWobbleEast;
        
        unsigned int fCounter;
        unsigned int fDBNrows;
        vector< unsigned int > fDBMJD;
        vector< double > fDBTime;
        vector< float > fDBTelElevationRaw;
        vector< float > fDBTelAzimuthRaw;
        vector< float > fDBTelElevation;
        vector< float > fDBTelAzimuth;
        vector< float > fDBTelExpectedElevation;
        vector< float > fDBTelExpectedAzimuth;
        
        int fNWarnings;
        
        
        VDB_Connection* fmy_connection;
        
        VTrackingCorrections* fTrackingCorrections;
        string fTPointCorrectionDate;
        
        bool getDBRunInfo();
        void getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip );
        void getDBSourceCoordinates( string iSource, float& iEVNTargetDec, float& iEVNTargetRA );
        bool readPointingFromDB();
        bool readPointingCalibratedVPMFromDB();
        bool readPointingUncalibratedVPMFromDB();
        bool readPointingFromVPMTextFile( string );
        
        void delete_myconnection()
        {
            if( fmy_connection )
            {
                delete fmy_connection;
                fmy_connection = 0;
            }
        }
        
    public:
    
        VPointingDB( unsigned int iTelID, unsigned int iRun );
        ~VPointingDB()
        {
            delete_myconnection();
        }
        
        
        bool   isGood()
        {
            return fStatus;
        }
        unsigned int getEventStatus()
        {
            return fEventStatus;
        }
        string getSourceNameDB()
        {
            return fDBSourceName;
        }
        float  getTargetDecDB()
        {
            return fDBTargetDec;
        }
        float  getTargetRADB()
        {
            return fDBTargetRA;
        }
        float  getTelExpectedAzimuthDB()
        {
            return fTelExpectedAzimuth;
        }
        float  getTelExpectedElevationDB()
        {
            return fTelExpectedElevation;
        }
        float  getTelAzimuthDB()
        {
            return fTelAzimuth;
        }
        float  getTelElevationDB()
        {
            return fTelElevation;
        }
        TTree* getTreePointingDB();
        float  getWobbleNorthDB()
        {
            return fDBWobbleNorth;
        }
        float  getWobbleEastDB()
        {
            return fDBWobbleEast;
        }
        unsigned int getTelID()
        {
            return fTelID;
        }
        bool   initialize( string iTPointCorrection, string iVPMDirectory, bool iVPMDB, bool iUncalibratedVPM );
        void   setObservatory( double iLongitude_deg = 0., double iLatitude_deg = 0. );
        bool   terminate();
        bool   updatePointing( int MJD, double iTime );
        void   getHorizonCoordinates( int MJD, double time, double dec, double ra, double& az, double& ze );
};
#endif
