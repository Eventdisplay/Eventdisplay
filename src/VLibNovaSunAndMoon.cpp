/*! \class VLibNovaSunAndMoon
    \brief get moon data (elevation, rise time, etc.)


    VERITAS basecamp location hardcoded
*/

#include "VLibNovaSunAndMoon.h"

VLibNovaSunAndMoon::VLibNovaSunAndMoon( double ilong, double ilat )
{
    reset();
    
    fObserver.lng = ilong;
    fObserver.lat = ilat;
}

VLibNovaSunAndMoon::VLibNovaSunAndMoon( int year, int month, int day, int hour, int minute, int second )
{
    reset();
    setDate( year, month, day, hour, minute, second );
}

VLibNovaSunAndMoon::VLibNovaSunAndMoon( double MJD )
{
    reset();
    setMJD( MJD );
}


void VLibNovaSunAndMoon::reset()
{
    // VERITAS basecamp (default values)
    fObserver.lat = 31.675;
    fObserver.lng = 249.0478;
    
    fJD = 0.;
}

double VLibNovaSunAndMoon::getMoonDisk()
{
    return ln_get_lunar_disk( fJD );
}

void VLibNovaSunAndMoon::getMoonElevationAzimuth( double& el, double& az )
{
    struct ln_equ_posn equ;
    struct ln_hrz_posn hrz;
    
    // get lunar position in equatorial coordinates
    ln_get_lunar_equ_coords( fJD, &equ );
    // transform to horizontal coordinates
    ln_get_hrz_from_equ( &equ, &fObserver, fJD, &hrz );
    
    el = hrz.alt;
    az = hrz.az;
    
}

double VLibNovaSunAndMoon::getMoonElevation()
{
    struct ln_equ_posn equ;
    struct ln_hrz_posn hrz;
    
    // get lunar position in equatorial coordinates
    ln_get_lunar_equ_coords( fJD, &equ );
    // transform to horizontal coordinates
    ln_get_hrz_from_equ( &equ, &fObserver, fJD, &hrz );
    
    return hrz.alt;
}

double VLibNovaSunAndMoon::getMoonAzimuth()
{
    struct ln_equ_posn equ;
    struct ln_hrz_posn hrz;
    
    // get lunar position in equatorial coordinates
    ln_get_lunar_equ_coords( fJD, &equ );
    // transform to horizontal coordinates
    ln_get_hrz_from_equ( &equ, &fObserver, fJD, &hrz );
    
    return hrz.az;
}

void VLibNovaSunAndMoon::setDate( int year, int month, int day, int hour, int minute, int second )
{
    struct ln_date date;
    date.years = year;
    date.months = month;
    date.days = day;
    date.hours = hour;
    date.minutes = minute;
    date.seconds = second;
    
    fJD = ln_get_julian_day( &date );
    cout << "Julian date: " << setprecision( 10 ) << fJD << endl;
}

void VLibNovaSunAndMoon::setMJD( double MJD )
{
    fJD = MJD + 2400000.5;
}

double VLibNovaSunAndMoon::getSunRiseTime()
{
    struct ln_rst_time rst;
    
    ln_get_solar_rst( fJD, &fObserver, &rst );
    
    return ( rst.rise - 2400000.5 );
}

double VLibNovaSunAndMoon::getSunSetTime()
{
    struct ln_rst_time rst;
    
    ln_get_solar_rst( fJD, &fObserver, &rst );
    
    return ( rst.set - 2400000.5 );
}

double VLibNovaSunAndMoon::getSunElevation()
{
    struct ln_equ_posn equ;
    struct ln_hrz_posn hrz;
    
    // get lunar position in equatorial coordinates
    ln_get_solar_equ_coords( fJD, &equ );
    // transform to horizontal coordinates
    ln_get_hrz_from_equ( &equ, &fObserver, fJD, &hrz );
    
    return hrz.alt;
}

