//! VStar get altitude and azimuth for a star
// libnova package needed

#include <iostream>

#include <libnova/julian_day.h>
#include <libnova/transform.h>

using namespace std;

class VLibNovaStar
{
    private:
    
        struct ln_lnlat_posn fObserver;
        double fJD;
        
        struct ln_equ_posn fStar;
        
    public:
    
        VLibNovaStar( double ilong = 249.0478, double ilat = 31.675 );
        VLibNovaStar( double ra, double dec, double ilong = 249.0478, double ilat = 31.675 );
        void setStar( double ra, double dec );
        void getElevationAzimuth( double MJD, double& el, double& az );
        double getElevation( double MJD );
        double getAzimuth( double MJD );
};

