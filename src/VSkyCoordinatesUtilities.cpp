/*! \namespace VSkyCoordinatesUtilities
    \brief  utilities for timing and coordinate transformations

*/

#include "VSkyCoordinatesUtilities.h"

int VSkyCoordinatesUtilities::getMJD_from_SQLstring( string iSQLData, double& i_mjd, double& sec_of_day )
{
    if( iSQLData.size() < 19 )
    {
        i_mjd = 0;
        sec_of_day = 0;
        return -99;
    }
    int i_year  = atoi( iSQLData.substr( 0, 4 ).c_str() );
    int i_month = atoi( iSQLData.substr( 5, 2 ).c_str() );
    int i_day   = atoi( iSQLData.substr( 8, 2 ).c_str() );
    
    // MJD
    int i_stat = 0;
    if( i_year > 0 )
    {
        VAstronometry::vlaCldj( i_year, i_month, i_day, &i_mjd, &i_stat );
    }
    else
    {
        i_mjd = 0;
    }
    
    // fraction of day
    sec_of_day  = atof( iSQLData.substr( 11, 2 ).c_str() ) * 60.*60.;
    sec_of_day += atof( iSQLData.substr( 14, 2 ).c_str() ) * 60.;
    sec_of_day += atof( iSQLData.substr( 17, 2 ).c_str() );
    
    return i_stat;
}

double VSkyCoordinatesUtilities::getMJD( int i_year, int i_month, int i_day )
{
    double i_mjd = 0.;
    int i_stat = 0;
    VAstronometry::vlaCldj( i_year, i_month, i_day, &i_mjd, &i_stat );
    
    return i_mjd;
}


double VSkyCoordinatesUtilities::getUTC( int i_mjd, double i_seconds )
{
    //! and add the fractional day to the mjd to get UTC
    double i_utc = i_mjd + i_seconds / 24. / 60. / 60.;
    return i_utc;
}


void VSkyCoordinatesUtilities::rotate( const double theta_rad, double& x, double& y )
{
    const double s = sin( theta_rad );
    const double c = cos( theta_rad );
    const double _x = x * c - y * s;
    const double _y = y * c + x * s;
    x = _x;
    y = _y;
}


/*!
    calculate wobble offsets in ra,dec

    (all angles in degrees)

    based on wobble_offset() in wobble.cpp from SFegan
*/
void VSkyCoordinatesUtilities::getWobbleOffset_in_RADec( double iNorth, double iEast, double idec, double ira, double& idiffdec, double& idiffra )
{
    idec *= TMath::DegToRad();
    ira  *= TMath::DegToRad();
    
    double x = 0.;
    double y = 0.;
    double z = 1.;
    double theta_rad = sqrt( iNorth * iNorth + iEast * iEast ) * TMath::DegToRad();
    double phi_rad = -1.*atan2( iEast, iNorth );
    if( phi_rad < 0. )
    {
        phi_rad += TMath::TwoPi();
    }
    
    VSkyCoordinatesUtilities::rotate( -theta_rad, z, x );
    VSkyCoordinatesUtilities::rotate( phi_rad, y, x );
    // declination
    VSkyCoordinatesUtilities::rotate( TMath::PiOver2() - idec, z, x );
    idiffdec = ( atan2( z, sqrt( x * x + y * y ) ) - idec ) * TMath::RadToDeg();
    // right ascension
    idiffra = atan2( y, x );
    if( idiffra < 0. )
    {
        idiffra += TMath::TwoPi();
    }
    idiffra *= -1.*TMath::RadToDeg();
    
    if( TMath::Abs( idiffra ) < 1.e-9 )
    {
        idiffra  = 0.;
    }
    if( TMath::Abs( idiffdec ) < 1.e-9 )
    {
        idiffdec = 0.;
    }
}

/*
    get wobbled direction

    (all angles in degrees)

    based on wobble() in wobble.cpp from SFegan
*/
void VSkyCoordinatesUtilities::getWobbledDirection( double iNorth, double iEast, double idec, double ira, double& dec_W, double& ra_W )
{
    dec_W = idec * TMath::DegToRad();
    ra_W  = ira * TMath::DegToRad();
    
    double x = 0.;
    double y = 0.;
    double z = 1.;
    double theta_rad = sqrt( iNorth * iNorth + iEast * iEast ) * TMath::DegToRad();
    double phi_rad = -1.*atan2( iEast, iNorth );
    if( phi_rad < 0. )
    {
        phi_rad += TMath::TwoPi();
    }
    
    VSkyCoordinatesUtilities::rotate( -theta_rad, z, x );
    VSkyCoordinatesUtilities::rotate( phi_rad, y, x );
    VSkyCoordinatesUtilities::rotate( TMath::PiOver2() - dec_W, z, x );
    VSkyCoordinatesUtilities::rotate( ra_W, x, y );
    ra_W = atan2( y, x );
    if( ra_W < 0. )
    {
        ra_W += TMath::TwoPi();
    }
    dec_W = atan2( z, sqrt( x * x + y * y ) );
    
    dec_W *= TMath::RadToDeg();
    ra_W  *= TMath::RadToDeg();
}


/*

   starting (iMJD_start) and ending epoch (iMJD_end)

   [rad] or [deg] (depending on iUnitIsDeg)

*/
void VSkyCoordinatesUtilities::precessTarget( double iMJD_end, double& ra_rad, double& dec_rad, double iMJD_start, bool iUnitIsDeg )
{
    if( iUnitIsDeg )
    {
        ra_rad *= TMath::DegToRad();
        dec_rad *= TMath::DegToRad();
    }
    
    VAstronometry::vlaPreces( iMJD_start, iMJD_end, &ra_rad, &dec_rad );
    
    if( iUnitIsDeg )
    {
        ra_rad *= TMath::RadToDeg();
        dec_rad *= TMath::RadToDeg();
    }
}



/*!
    assume continuous increasing azimuth
*/
double VSkyCoordinatesUtilities::addToMeanAzimuth( double iMean, double iAz )
{
    if( iMean > 270. && iAz < 90. )
    {
        iMean += iAz + 360.;
    }
    else
    {
        iMean += iAz;
    }
    
    return iMean;
}


/*
     difference between to pointing directions in camera coordinates

     input coordinates in [deg]
*/
void VSkyCoordinatesUtilities::getDifferenceInCameraCoordinates( double tel_ze, double tel_az,
        double shower_ze,  double shower_az,
        float& x, float& y, float& z )
{
    // convert coordinates from [deg] to [rad]
    tel_az    *= TMath::DegToRad();
    shower_az *= TMath::DegToRad();
    double tel_el    = ( 90. - tel_ze ) * TMath::DegToRad();
    double shower_el = ( 90. - shower_ze ) * TMath::DegToRad();
    
    double cx = cos( shower_el ) * sin( shower_az );
    double cy = cos( shower_el ) * cos( shower_az );
    double cz = sin( shower_el );
    
    double i_temp = sin( tel_az ) * cx + cos( tel_az ) * cy;
    
    x = ( cos( tel_az ) * cx - sin( tel_az ) * cy ) * TMath::RadToDeg();
    z = ( cos( tel_el ) * i_temp + sin( tel_el ) * cz );
    y = ( -1.*sin( tel_el ) * i_temp + cos( tel_el ) * cz ) * TMath::RadToDeg();
    y *= -1.;
    
    if( fabs( x ) < 1.e-4 )
    {
        x = 0.;
    }
    if( fabs( y ) < 1.e-4 )
    {
        y = 0.;
    }
    if( fabs( z ) < 1.e-4 )
    {
        z = 0.;
    }
}

/*
    calculate shower direction from telescope pointing and reconstruction shower direction

     small angle approximation, assume small x,y (neglect z)
*/

void VSkyCoordinatesUtilities::getRotatedShowerDirection( double ze, double az, double y, double x, double& rze, double& raz )
{
    // get all directions in [rad]
    x /= TMath::RadToDeg();
    y /= ( -1.*TMath::RadToDeg() );
    // assume all telescopes point in same directions
    double el = ( 90. - ze ) / TMath::RadToDeg();
    az = az / TMath::RadToDeg();
    // these are the resulting directions
    
    double r = sqrt( 1. + x * x + y * y );
    double cx = x / r;
    double cy = 1. / r;
    double cz = y / r;
    
    // rotate telescope around elevation axis
    double ex = cx;
    double ey = cy * cos( el ) - cz * sin( el );
    double ez = cy * sin( el ) + cz * cos( el );
    // rotate around azimuth
    double rx, ry, rz;
    rx =     ex * cos( az ) + ey * sin( az );
    ry = -1.*ex * sin( az ) + ey * cos( az );
    rz = ez;
    // calculate new azimuth, zenith
    r = sqrt( rx * rx + ry * ry );
    // small value check
    if( fabs( r ) < 1.e-10 )
    {
        r = 0.;
    }
    if( fabs( rx ) < 1.e-10 )
    {
        rx = 0.;
    }
    if( fabs( ry ) < 1.e-10 )
    {
        ry = 0.;
    }
    if( fabs( rz ) < 1.e-10 )
    {
        rz = 0.;
    }
    
    if( r == 0. )
    {
        raz = az * TMath::RadToDeg();
    }
    else
    {
        raz = ( TMath::Pi() / 2. - atan2( ry, rx ) ) * TMath::RadToDeg();
        if( raz > 180. )
        {
            raz = -1.*( 360. - raz );
        }
        if( raz < -180. )
        {
            raz *= -1.;
        }
    }
    if( rz == 0. )
    {
        rze = 90. - el * TMath::RadToDeg();
    }
    else
    {
        rze = 90. - atan2( rz, r ) * TMath::RadToDeg();
    }
}


double VSkyCoordinatesUtilities::adjustAzimuthToRange( double az_deg )
{
    return VAstronometry::vlaDranrm( az_deg * TMath::DegToRad() ) * TMath::RadToDeg();
}

/*

   return values in [deg]

*/
double VSkyCoordinatesUtilities::getTargetShiftWest( double iTargetRA_deg, double iTargetDec_deg, double ira_deg, double idec_deg )
{
    double sep  = VAstronometry::vlaDsep( iTargetRA_deg * TMath::DegToRad(), iTargetDec_deg * TMath::DegToRad(),
                           ira_deg * TMath::DegToRad(), idec_deg * TMath::DegToRad() );
    double bear = VAstronometry::vlaDbear( iTargetRA_deg * TMath::DegToRad(), iTargetDec_deg * TMath::DegToRad(),
                            ira_deg * TMath::DegToRad(), idec_deg * TMath::DegToRad() );
                            
    double iShift = sep * sin( bear ) * TMath::RadToDeg();
    
    if( TMath::Abs( iShift ) < 1.e-8 )
    {
        iShift = 0.;
    }
    
    return iShift;
}


double VSkyCoordinatesUtilities::getTargetShiftNorth( double iTargetRA_deg, double iTargetDec_deg, double ira_deg, double idec_deg )
{
    double sep  = VAstronometry::vlaDsep( iTargetRA_deg * TMath::DegToRad(), iTargetDec_deg * TMath::DegToRad(),
                           ira_deg * TMath::DegToRad(), idec_deg * TMath::DegToRad() );
    double bear = VAstronometry::vlaDbear( iTargetRA_deg * TMath::DegToRad(), iTargetDec_deg * TMath::DegToRad(),
                            ira_deg * TMath::DegToRad(), idec_deg * TMath::DegToRad() );
                            
    double iShift = sep * cos( bear ) * TMath::RadToDeg();
    
    if( TMath::Abs( iShift ) < 1.e-8 )
    {
        iShift = 0.;
    }
    
    return iShift;
}
/*

   convert x,y in derotated coordinates of current epoch (iMJD) into x,y in J2000

   use Az/Elevation for pointing

*/
void VSkyCoordinatesUtilities::convert_derotatedCoordinates_to_J2000( int iMJD, double iTime,
        double iTelAz, double iTelEl,
        double& x, double& y )
{
    // get equatorial coordinates for current epoch
    double i_ra = 0.;
    double i_dec = 0.;
    getEquatorialCoordinates( iMJD, iTime, iTelAz, 90. - iTelEl, i_dec, i_ra );
    // precess cooridanates to J2000
    double i_raJ2000 = i_ra;
    double i_decJ2000 = i_dec;
    precessTarget( 51544., i_raJ2000, i_decJ2000, iMJD, true );
    
    // calculate wobble offset in ra/dec for current epoch
    double i_decDiff = 0.;
    double i_raDiff = 0.;
    getWobbleOffset_in_RADec( y, -x, i_dec, i_ra, i_decDiff, i_raDiff );
    if( i_raDiff < -180. )
    {
        i_raDiff += 360.;
    }
    double i_decWobble = i_dec + i_decDiff;
    double i_raWobble  = i_ra + i_raDiff;
    
    // correct for precession (from current epoch to J2000=MJD51544)
    precessTarget( 51544., i_raWobble, i_decWobble, iMJD, true );
    x = getTargetShiftWest( i_raJ2000, i_decJ2000, i_raWobble, i_decWobble ) * -1.;
    y = getTargetShiftNorth( i_raJ2000, i_decJ2000, i_raWobble, i_decWobble );
}

/*

   convert x,y in derotated coordinates of current epoch (iMJD) into x,y in J2000

   use RA/DEC in J2000 for pointing

*/
void VSkyCoordinatesUtilities::convert_derotatedCoordinates_to_J2000( double iMJD, double i_RA_J2000_deg, double i_DEC_J2000_deg, double& x, double& y )
{
    double i_ra = i_RA_J2000_deg * TMath::DegToRad();
    double i_dec = i_DEC_J2000_deg * TMath::DegToRad();
    precessTarget( iMJD, i_ra, i_dec );
    
    // calculate wobble offset in ra/dec for current epoch
    double i_decDiff = 0.;
    double i_raDiff = 0.;
    getWobbleOffset_in_RADec( y, -x, i_dec * TMath::RadToDeg(), i_ra * TMath::RadToDeg(), i_decDiff, i_raDiff );
    if( i_raDiff < -180. )
    {
        i_raDiff += 360.;
    }
    double i_decWobble = i_dec * TMath::RadToDeg() + i_decDiff;
    double i_raWobble  = i_ra * TMath::RadToDeg()  + i_raDiff;
    
    // correct for precession (from current epoch to J2000=MJD51544)
    precessTarget( 51544., i_raWobble, i_decWobble, iMJD, true );
    x = getTargetShiftWest( i_RA_J2000_deg, i_DEC_J2000_deg, i_raWobble, i_decWobble ) * -1.;
    y = getTargetShiftNorth( i_RA_J2000_deg, i_DEC_J2000_deg, i_raWobble, i_decWobble );
}

/*
   calculate coordinates of camera centre in J2000
*/
void VSkyCoordinatesUtilities::getCameraCentreCoordinates_J2000( double iMJD, double i_Target_RA_J200_deg, double i_Target_Dec_J2000_deg,
        double iNorth_deg, double i_East_deg,
        double& i_C_RA_J2000_deg, double& i_C_Dec_J2000_deg )
{
    double i_dec = i_Target_Dec_J2000_deg * TMath::DegToRad();
    double i_ra  = i_Target_RA_J200_deg * TMath::DegToRad();
    
    // precess to current epoch
    precessTarget( iMJD, i_ra, i_dec );
    
    // calculate wobble offset in current epoch
    double i_decDiff = 0.;
    double i_raDiff = 0.;
    getWobbleOffset_in_RADec( iNorth_deg, i_East_deg,
                              i_dec * TMath::RadToDeg(), i_ra * TMath::RadToDeg(),
                              i_decDiff, i_raDiff );
    if( i_raDiff < -180. )
    {
        i_raDiff += 360.;
    }
    
    i_C_RA_J2000_deg  = i_ra * TMath::RadToDeg() + i_raDiff;
    i_C_Dec_J2000_deg = i_dec * TMath::RadToDeg() + i_decDiff;
    
    // precess to J2000
    precessTarget( 51544., i_C_RA_J2000_deg, i_C_Dec_J2000_deg, iMJD, true );
}



/*!
  calculate angular distance between two directions

  input: all angles in radiant
  output: angular distance in degree

*/
double VSkyCoordinatesUtilities::angularDistance( double Az, double Ze, double Traz, double Trze )
{
    double value;
    
    value  = sin( Ze ) * sin( Trze ) * cos( ( Az - Traz ) );
    value += cos( Ze ) * cos( Trze );
    // limited accuracy results sometimes in values slightly larger than 1
    if( value > 1. )
    {
        value = 1.;
    }
    value = acos( value );
    value *= TMath::RadToDeg();
    
    return value;
}

/*
 *
 *  get equatorial from horizontal coordinates
 *
 *  all input/output angles in [deg]
 */
void VSkyCoordinatesUtilities::getEquatorialCoordinates( int MJD, double time, double az_deg, double ze_deg, double& dec_deg, double& ra_deg )
{
    // convert time to fraction of a day
    double iTime = time / 86400.;
    // transform coordinates
    double ha = 0.;
    VAstronometry::vlaDh2e( az_deg * TMath::DegToRad(), ( 90. - ze_deg ) * TMath::DegToRad(), VGlobalRunParameter::getObservatory_Latitude_deg() * TMath::DegToRad(), &ha, &dec_deg );
    // convert hour angle into ra
    // get Greenwich sideral time
    double iSid = VAstronometry::vlaGmsta( ( double )MJD, iTime );
    // calculate local sideral time
    iSid = iSid - VGlobalRunParameter::getObservatory_Longitude_deg() * TMath::DegToRad();
    // calculate right ascension
    ra_deg = VAstronometry::vlaDranrm( iSid - ha );
    // from [rad] to [deg]
    dec_deg *= TMath::RadToDeg();
    ra_deg  *= TMath::RadToDeg();
}

void VSkyCoordinatesUtilities::getHorizontalCoordinates( int MJD, double time, double dec_deg, double ra_deg, double& az_deg, double& ze_deg )
{
    // convert time to fraction of a day
    double iTime = time / 86400.;
    // get Greenwich sideral time
    double iSid = VAstronometry::vlaGmsta( ( double )MJD, iTime );
    // calculate local sideral time
    iSid = iSid - VGlobalRunParameter::getObservatory_Longitude_deg() * TMath::DegToRad();
    // calculate hour angle
    double ha = VAstronometry::vlaDranrm( iSid - ra_deg * TMath::DegToRad() );
    // get horizontal coordinates
    VAstronometry::vlaDe2h( ha, dec_deg * TMath::DegToRad(), VGlobalRunParameter::getObservatory_Latitude_deg() * TMath::DegToRad(), &az_deg, &ze_deg );
    // from [rad] to [deg]
    ze_deg = 90 - ze_deg * TMath::RadToDeg();
    az_deg *= TMath::RadToDeg();
}

double VSkyCoordinatesUtilities::getRightAscension_inDegrees_fromHour( double h, double m, double s )
{
    double i_RA_deg = h;
    
    i_RA_deg += m / 60.;
    i_RA_deg += s / 3600.;
    
    return i_RA_deg / 24. * 360.;
}

double VSkyCoordinatesUtilities::getDeclination_inDegrees_fromHour( double h, double m, double s )
{
    double i_DEC_deg = h;
    if( i_DEC_deg < 0. )
    {
        i_DEC_deg -= m / 60.;
        i_DEC_deg -= s / 3600.;
    }
    else
    {
        i_DEC_deg += m / 60.;
        i_DEC_deg += s / 3600.;
    }
    
    return i_DEC_deg;
}

/*
 * returns the camera derotation angle in radians
*/
double VSkyCoordinatesUtilities::getDerotationAngle( double i_UTC, double iTelRA, double iTelDec,
        double iObservatoryLongitude, double iObservatoryLatitude )
{
    return  -1.*VAstronometry::vlaPa( getHourAngle( i_UTC, iTelRA, iObservatoryLongitude ),
                    iTelDec, iObservatoryLatitude );
}

double VSkyCoordinatesUtilities::getDerotationAngle( double iMJD, double iTime, double iTelRA, double iTelDec,
        double iObservatoryLongitude, double iObservatoryLatitude )
{
    return getDerotationAngle( getUTC( iMJD, iTime ), iTelRA, iTelDec, iObservatoryLongitude, iObservatoryLatitude );
}


double VSkyCoordinatesUtilities::getSidereal( double i_UTC, double iObservatoryLongitude )
{
    return VAstronometry::vlaDranrm( VAstronometry::vlaGmst( i_UTC ) - iObservatoryLongitude );
}


double VSkyCoordinatesUtilities::getHourAngle( double i_UTC, double iTelRA, double iObservatoryLongitude )
{
    return VAstronometry::vlaDranrm( getSidereal( i_UTC, iObservatoryLongitude ) - iTelRA );
}

