/* \class VAstronometry
 * \brief wrapper class for positional astronomy access
 * 
 * Allows to use
 *
 * - slalib
 * - sofa (preferred)
 *
*/
 
#include "VAstronometry.h"

VAstronometry::VAstronometry()
{
}

/*
 * Modified Julian Date to Gregorian year, month, day, and fraction of a day.
 *
*/
void VAstronometry::vlaDjcl( double djm, int* iy, int* im, int* id, double* fd, int* j )
{
#ifdef ASTROSLALIB
    slaDjcl( djm, iy, im, id, fd, j );
#elif ASTROSOFA
    *j = iauJd2cal( ASTRODJM0, djm, iy, im, id, fd );
#endif
}

/*
 * Gregorian calendar to Modified Julian Date.
 */
void VAstronometry::vlaCldj( int iy, int im, int id, double* djm, int* j )
{
#ifdef ASTROSLALIB
    slaCldj( iy, im, id, djm, j );
#elif ASTROSOFA
    double djm0 = 0.;
    *j = iauCal2jd( iy, im, id, &djm0, djm );
#endif
}

/*
 * Conversion from Universal Time to Greenwich mean sidereal time, with rounding errors minimized
*/
double VAstronometry::vlaGmsta( double date, double ut1 )
{
#ifdef ASTROSLALIB
    return slaGmsta( date, ut1 );
#elif ASTROSOFA
    date += DJM0;
    return iauGmst06( date, ut1, date, ut1 );
#endif
}

/*
 * Conversion from Universal Time to Sidereal Time.
*/
double VAstronometry::vlaGmst( double ut1 )
{
#ifdef ASTROSLALIB
    return slaGmst( ut1 );
#elif ASTROSOFA
    return iauGmst06( DJM0, ut1, DJM0, ut1 );
#endif
}

/*
 * Transformation from J2000.0 equatorial coordinates to IAU 1958 Galactic coordinates.
*/
void VAstronometry::vlaEqgal( double dr, double dd, double* dl, double* db )
{
#ifdef ASTROSLALIB
    slaEqgal( dr, dd, dl, db );
#elif ASTROSOFA
    iauIcrs2g( dr, dd, dl, db );
#endif
}

/*
 * Equatorial to horizon coordinates:  HA,Dec to Az,El
*/

void VAstronometry::vlaDe2h( double ha, double dec, double phi, double* az, double* el )
{
#ifdef ASTROSLALIB
    slaDe2h( ha, dec, phi, az, el );
#elif ASTROSOFA
    iauHd2ae( ha, dec, phi, az, el );
#endif
}


/*
 *  Precession
 *
 *  assume always FK5 precession
 *
 */
void VAstronometry::vlaPreces( double MJD_ep0, double MJD_ep1, double *ra, double *dc )
{
#ifdef ASTROSLALIB
    int  oy, om, od, j, ny, nd;
    double ofd, ep0, ep1;
    slaDjcl( MJD_ep0, &oy, &om, &od, &ofd, &j );
    slaClyd( oy, om, od, &ny, &nd, &j );
    ep0  = ny + nd / 365.25;
    slaDjcl( MJD_ep1, &oy, &om, &od, &ofd, &j );
    slaClyd( oy, om, od, &ny, &nd, &j );
    ep1  = ny + nd / 365.25;

    slaPreces( "FK5", ep0, ep1, ra, dc );

#elif ASTROSOFA
    
   // precession matrix
   double rot_prec[3][3];

   // days since year 2000
   double ep0_days_2000 = MJD_ep0 - DJM00;
   double ep1_days_2000 = MJD_ep1 - DJM00;

   /////////////////////////////////
   // Three cases in the following

   // MJD_ep0 is J2000
   if( ep0_days_2000 == 2000.0 )
   {
      iauPmat06( DJ00, ep1_days_2000, rot_prec );
   }
   // MJD_ep1 is J2000
   else if( ep1_days_2000 == 2000.0 )
   {
      iauPmat06( DJ00, ep0_days_2000, rot_prec );
      iauTr( rot_prec, rot_prec );
   }
   // all other cases
   else
   {
      double temp_tot_mat[3][3];
      iauPmat06( DJ00, ep0_days_2000, rot_prec );
      iauTr( rot_prec, rot_prec);
      iauPmat06( DJ00, ep1_days_2000, temp_tot_mat );
      iauRxr( rot_prec, temp_tot_mat, rot_prec );
   }
    double e1[3];
    double e2[3];

    // Convert spherical coordinates to Cartesian
    iauS2c( *ra, *dc, e1 );
    // apply precession matrix
    iauRxp( rot_prec, e1, e2 );
    // P-vector to spherical coordinates
    iauC2s( e2, ra, dc );
    *ra = iauAnp( *ra );
#endif
}

/*
 * Normalize angle into range 0-2 pi.
 */
double VAstronometry::vlaDranrm( double angle )
{
#ifdef ASTROSLALIB
    return slaDranrm( angle );
#elif ASTROSOFA
    return iauAnp( angle );
#endif
}

/*
 * Angle between two points on a sphere
*/
double VAstronometry::vlaDsep( double a1, double b1, double a2, double b2 )
{
#ifdef ASTROSLALIB
    return slaDsep( a1, b1, a2, b2 );
#elif ASTROSOFA
    return iauSeps( a1, b1, a2, b2 );
#endif
}

/*
 * Transform tangent plane coordinates into spherical
*/
void VAstronometry::vlaDtp2s( double xi, double eta, double raz, double decz,
        double* ra, double* dec )
{
#ifdef ASTROSLALIB
    slaDtp2s( xi, eta, raz, decz, ra, dec );
#elif ASTROSOFA
    iauTpsts( xi, eta, raz, decz, ra, dec );
#endif
}

/*
 *  Projection of spherical coordinates onto tangent plane ('gnomonic' projection - 'standard coordinates').
*/
void VAstronometry::vlaDs2tp( double ra, double dec, double raz, double decz,
        double* xi, double* eta, int* j )
{
#ifdef ASTROSLALIB
    slaDs2tp( ra, dec, raz, decz, xi, eta, j );
#elif ASTROSOFA
    *j = iauTpxes( ra, dec, raz, decz, xi, eta );
#endif
}

/*
 * Convert an angle in radians to hours, minutes, seconds
*/
void VAstronometry::vlaDr2tf( int ndp, double angle, char* sign, int ihmsf[4] )
{
#ifdef ASTROSLALIB
    slaDr2tf( ndp, angle, sign, ihmsf );
#elif ASTROSOFA
    iauA2tf( ndp, angle, sign, ihmsf );
#endif
}

/*
 * Bearing (position angle) of one point on a sphere relative to another.
*/

double VAstronometry::vlaDbear(double a1, double b1, double a2, double b2 )
{
#ifdef ASTROSLALIB
    return slaDbear( a1, b1, a2, b2 );
#elif ASTROSOFA
    return iauPas( a1, b1, a2, b2 );
#endif
}

/*
 *  Horizon to equatorial coordinates:  Az,El to HA,Dec
*/
void VAstronometry::vlaDh2e( double az, double el, double phi, double* ha, double* dec )
{
#ifdef ASTROSLALIB
    slaDh2e( az, el, phi, ha, dec );
#elif ASTROSOFA
    iauHd2ae( az, el, phi, ha, dec );
#endif
}

/*
 * HA, Dec to Parallactic Angle.
 */
double VAstronometry::vlaPa( double ha, double dec, double phi )
{
#ifdef ASTROSLALIB
    return slaPa( ha, dec, phi );
#elif ASTROSOFA
    return iauHd2pa( ha, dec, phi );
#endif
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
string VAstronometry::getAstronometryLibrary()
{
#ifdef ASTROSLALIB
    return "SLALIB";
#elif ASTROSOFA
    return "SOFA";
#endif
    return "NOTSET";
}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/* 
 * testing routines
 *
*/

void VAstronometry::test_vlaDjcl()
{
     double fd;
     int iy, im, id, j;

     VAstronometry::vlaDjcl( 58401.87500000, &iy, &im, &id, &fd, &j );
     cout << "VAstronometry::vlaDjcl " << getAstronometryLibrary();
     cout << " (2018, 10, 10, 0.853) \t\t";
     cout << iy << ", " << im << ", " << id << ", " << fd << ", " << endl;
}

void VAstronometry::test_vlaCldj()
{
    int j;
    double djm;
    VAstronometry::vlaCldj( 2018, 10, 10, &djm, &j );
    cout << "VAstronometry::vlaCldj " << getAstronometryLibrary();
    cout << " (58401) \t\t";
    cout << djm << "\t" << j << endl;
}

void VAstronometry::test_vlaPreces()
{
    double MJD_ep0 = 2451545.0 - 2400000.5;
    double MJD_ep1 = 58082.4;
    double ra = 83.6333 * DD2R;
    double dec = 22.0145 * DD2R;

    vlaPreces( MJD_ep0, MJD_ep1, &ra, &dec );

    cout << "VAstronometry::vlaPreces " << getAstronometryLibrary();
    cout << " (83.9027, 22.0253) \t\t";
    cout << setprecision(16) << ra * DR2D << "\t" << dec * DR2D << endl;
}

void VAstronometry::test_vlaGmsta()
{
    cout << "VAstronometry::vlaGmsta " << getAstronometryLibrary();
    cout << " (0.559428)  \t\t" << setprecision(16) << vlaGmsta(  50123.5, 0.2 );
    cout << endl;
}

void VAstronometry::test_vlaGmst()
{
    cout << "VAstronometry::vlaGmst " << getAstronometryLibrary();
    cout << " (0.559428)  \t\t" << setprecision(16) << vlaGmst( 0.2 );
    cout << endl;
}

void VAstronometry::test_vlaDranrm()
{
    cout << "VAstronometry::vlaDranrm " << getAstronometryLibrary();
    cout << " (6.08319) \t\t";
    cout << setprecision(16) << vlaDranrm( 2.*6.183185307179586477 ) << endl;
}

void VAstronometry::test_vlaDe2h()
{
    double az, el;
    VAstronometry::vlaDe2h( 1.1, 1.2, 0.3, &az, &el );
    cout << "VAstronometry::vlaDe2h " << getAstronometryLibrary();
    cout << setprecision(16) << " (" << 5.916889243730066194 * DR2D << ", " << 0.4472186304990486228 * DR2D << " ) \t\t";
    cout << az * DR2D << "\t" << el  * DR2D << endl;
}

void VAstronometry::test_vlaDh2e()
{
    double ha, dec;
    VAstronometry::vlaDh2e( 5.916889243730066194, 0.4472186304990486228, 0.3, &ha, &dec );
    cout << "VAstronometry::vlaDh2e " << getAstronometryLibrary();
    cout << " (1.1, 1.2 ) \t\t";
    cout << ha << ", " << dec << endl;
}


void VAstronometry::test_vlaEqgal()
{
    double b, l;
    VAstronometry::vlaEqgal( 83.633076 * DD2R, 22.014493 * DD2R, &b, &l );

    cout << "VAstronometry::vlaEqgal " << getAstronometryLibrary();
    cout << " (184.557451, -5.784369) \t\t";
    cout << setprecision(16) << b * DR2D << ", " << l * DR2D << endl;
}

void VAstronometry::test_vlaDsep()
{
    cout << "VAstronometry::vlaDsep " << getAstronometryLibrary();
    cout << " (2.346722016996998842) \t\t";
    cout << setprecision(16) << VAstronometry::vlaDsep( 1.0, 0.1, 0.2, -3. ) << endl;
}

void VAstronometry::test_vlaDtp2s()
{
    double ra, dec;
    VAstronometry::vlaDtp2s( -0.03, 0.07, 2.3, 1.5, &ra, &dec );
    cout << "VAstronometry::vlaDtp2s " << getAstronometryLibrary();
    cout << setprecision(16) << " (" << 0.7596127167359629775 * DR2D << ", " << 1.540864645109263028 * DR2D << ") \t\t";
    cout << ra * DR2D << ", " << dec * DR2D << endl;
}

void VAstronometry::test_vlaDs2tp()
{
    double xi, eta;
    int j;
    VAstronometry::vlaDs2tp( 0.7596127167359629775, 1.540864645109263028, 2.3, 1.5, &xi, &eta, &j );
    cout << "VAstronometry::vlaDs2tp " << getAstronometryLibrary();
    cout << setprecision(16) << " (-0.03, 0.07 ) \t\t" << xi << ", " << eta << ", " << j << endl;
}

void VAstronometry::test_vlaDr2tf()
{
    int ihmsf[4];
    char s;
    VAstronometry::vlaDr2tf( 4, -3.01234, &s, ihmsf );
    cout << "VAstronometry::vlaDr2tf " << getAstronometryLibrary();
    cout << " (-, 11, 30, 22, 6484) \t\t";
    cout << s << ", " << ihmsf[0] << ", " << ihmsf[1] << ", " << ihmsf[2] << ", " << ihmsf[3] << endl;
}

void VAstronometry::test_vlaDbear()
{
    cout << "VAstronometry::vlaDbear " << getAstronometryLibrary();
    cout << " (-2.724544922932270424) \t\t";
    cout << setprecision(16) << VAstronometry::vlaDbear( 1.0, 0.1, 0.2, -1. ) << endl;
}

void VAstronometry::test_vlaPa()
{
    cout << "VAstronometry::vlaPa " << getAstronometryLibrary();
    cout << " (1.906227428001995580) \t\t";
    cout << setprecision(16) << VAstronometry::vlaPa( 1.1, 1.2, 0.3 ) << endl;
}


void VAstronometry::test()
{
     test_vlaDjcl();
     test_vlaCldj();
     test_vlaDranrm();
     test_vlaDsep();
     test_vlaGmst();
     test_vlaGmsta();
     test_vlaGmst();
     test_vlaEqgal();
     test_vlaDe2h();
     test_vlaDh2e();
     test_vlaDtp2s();
     test_vlaDs2tp();
     test_vlaPreces();
     test_vlaDr2tf();
     test_vlaDbear();
     test_vlaPa();
}
