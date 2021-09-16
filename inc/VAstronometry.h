//! VAstronometry positional astronomy routines

#ifndef VAstronometry_H
#define VAstronometry_H

#include <iomanip>
#include <iostream>
#include <string>

#ifdef ASTROSLALIB
#include "VASlalib.h"
#elif ASTROSOFA
#include <sofa.h>
#include <sofam.h>
#endif

#define ASTRODJM0 (2400000.5)

using namespace std;

class VAstronometry
{
     private:

         static void test_vlaDjcl();
         static void test_vlaCldj();
         static void test_vlaPreces();
         static void test_vlaDranrm();
         static void test_vlaGmst();
         static void test_vlaGmsta();
         static void test_vlaEqgal();
         static void test_vlaDe2h();
         static void test_vlaDh2e();
         static void test_vlaDsep();
         static void test_vlaDs2tp();
         static void test_vlaDtp2s();
         static void test_vlaDr2tf();
         static void test_vlaDbear();
         static void test_vlaPa();

     public:

         VAstronometry();
         ~VAstronometry() {}

         static void vlaDjcl( double djm, int* iy, int* im, int* id, double* fd, int* j );
         static void vlaCldj( int iy, int im, int id, double* djm, int* j );
         static double vlaGmst( double ut1 );
         static double vlaGmsta( double date, double ut1 );

         static double vlaDranrm( double angle );
         static double vlaDsep( double a1, double b1, double a2, double b2 );
         static void vlaDs2tp( double ra, double dec, double raz, double decz,
                 double* xi, double* eta, int* j );
         static void vlaDtp2s( double xi, double eta, double raz, double decz,
                 double* ra, double* dec );
         static void vlaDr2tf( int ndp, double angle, char* sign, int ihmsf[4] );
         static double vlaDbear(double a1, double b1, double a2, double b2 );
         static double vlaPa( double ha, double dec, double phi );

         static void vlaEqgal( double dr, double dd, double* dl, double* db );
         static void vlaDe2h( double ha, double dec, double phi, double* az, double* el );
         static void vlaDh2e( double az, double el, double phi, double* ha, double* dec );
         static void vlaPreces( double MJD_ep0, double MJD_ep1, double *ra, double *dc );

         static string getAstronometryLibrary();
         static void test();
};

#endif
