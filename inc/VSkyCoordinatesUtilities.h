//! VSkyCoordinatesUtilities utilities for timing and coordinate transformations

#ifndef VSkyCoordinatesUtilities_H
#define VSkyCoordinatesUtilities_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"

#include "VAstronometry.h"
#include "VGlobalRunParameter.h"

namespace VSkyCoordinatesUtilities
{
    double addToMeanAzimuth( double iMean, double iAz );                   // mean azimuth calculation
    void   convert_derotatedCoordinates_to_J2000( double iMJD, double i_RA_J2000_deg, double i_DEC_J2000_deg, double& x, double& y );
    void   convert_derotatedCoordinates_to_J2000( int iMJD, double iTime, double i_TelAz, double iTelElevation, double& x, double& y );
    void   getCameraCentreCoordinates_J2000( double iMJD, double i_Target_RA_J200_deg, double i_Target_Dec_J2000_deg,
            double iNorth_deg, double i_East_deg,
            double& i_C_RA_J2000_deg, double& i_C_Dec_J2000_deg );
    void   getEquatorialCoordinates( int MJD, double time, double az, double ze, double& dec, double& ra );
    double getDerotationAngle( double MJD, double time, double iTelRA, double iTelDec, double iObservatoryLongitude, double iObservatoryLatitude );
    double getDerotationAngle( double i_UTC, double iTelRA, double iTelDec, double iObservatoryLongitude, double iObservatoryLatitude );
    double getRightAscension_inDegrees_fromHour( double h, double m, double s );
    double getHourAngle( double i_UTC, double iTelRA, double iObservatoryLongitude );
    double getSidereal( double i_UTC, double iObservatoryLongitude );
    double getDeclination_inDegrees_fromHour( double h, double m, double s );
    void   getHorizontalCoordinates( int MJD, double time, double dec_deg, double ra_deg, double& az_deg, double& ze_deg );
    double getMJD( int i_year, int i_month, int i_day );
    int    getMJD_from_SQLstring( string iSQLData, double& mjd, double& sec_of_day );
    double getTargetShiftWest( double iTargetRA_deg, double iTargetDec_deg, double ira_deg, double idec_deg );
    double getTargetShiftNorth( double iTargetRA_deg, double iTargetDec_deg, double ira_deg, double idec_deg );
    double getUTC( int i_mjd, double i_seconds );
    void   getWobbledDirection( double iNorth_deg, double iEast_deg, double idec_deg, double ira_deg, double& dec_W_deg, double& ra_W_deg );
    void   getWobbleOffset_in_RADec( double iNorth, double iEast, double idec, double ira, double& idiffdec, double& idiffra );
    void   precessTarget( double iMJD_end, double& ra_rad, double& dec_rad, double iMJD_start = 51544.5, bool bUnitIsDegrees = false );     // default start is J2000
    void   rotate( const double theta_rad, double& x, double& y );
}

#endif

