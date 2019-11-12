//! class VLibNovaSunAndMoon get moon data (elevation, rise time, etc.)
// libnova needed

#include <iomanip>
#include <iostream>

#include <libnova/lunar.h>
#include <libnova/julian_day.h>
#include <libnova/rise_set.h>
#include <libnova/solar.h>
#include <libnova/transform.h>

using namespace std;

class VLibNovaSunAndMoon
{
    private:
    
        struct ln_lnlat_posn fObserver;
        
        double fJD;
        
        void reset();
        
    public:
    
        VLibNovaSunAndMoon( double MJD );
        VLibNovaSunAndMoon( double ilong = 249.0478, double ilat = 31.675 );
        VLibNovaSunAndMoon( int year, int month, int day, int hour, int minute, int second );
        ~VLibNovaSunAndMoon() {}
        double getMJD()
        {
            return fJD + 2400000.5;
        }
        void   getMoonElevationAzimuth( double& el, double& az );
        double getMoonElevation();
        double getMoonAzimuth();
        double getMoonDisk();
        double getSunElevation();
        double getSunRiseTime();   //! (in MJD)
        double getSunSetTime();    //! (in MJD)
        void setDate( int year, int month, int day, int hour, int minute, int second );
        void setMJD( double MJD );
};

