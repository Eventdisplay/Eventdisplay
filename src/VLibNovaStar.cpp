/* \class VLibNovaStar
   \brief get altitude and azimuth for a star


   VERITAS basecamp location hardcoded
*/

#include "VLibNovaStar.h"

VLibNovaStar::VLibNovaStar( double ilong, double ilat )
{
    fStar.ra = 0.;
    fStar.dec = 0.;
    
    fObserver.lng = ilong;
    fObserver.lat = ilat;
}

VLibNovaStar::VLibNovaStar( double ra, double dec, double ilong, double ilat )
{
    fStar.ra = ra;
    fStar.dec = dec;
    
    fObserver.lng = ilong;
    fObserver.lat = ilat;
}

void VLibNovaStar::getElevationAzimuth( double MJD, double& el, double& az )
{
    struct ln_hrz_posn hrz;
    
    ln_get_hrz_from_equ( &fStar, &fObserver, MJD + 2400000.5, &hrz );
    
    az = hrz.az;
    el = hrz.alt;
}


double VLibNovaStar::getElevation( double MJD )
{
    struct ln_hrz_posn hrz;
    
    ln_get_hrz_from_equ( &fStar, &fObserver, MJD + 2400000.5, &hrz );
    
    return hrz.alt;
}

double VLibNovaStar::getAzimuth( double MJD )
{
    struct ln_hrz_posn hrz;
    
    ln_get_hrz_from_equ( &fStar, &fObserver, MJD + 2400000.5, &hrz );
    
    return hrz.az;
}
