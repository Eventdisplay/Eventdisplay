/*! \file printBinaryOrbitalPhase
    \brief print binary orbital phase


*/

#include <iomanip>
#include <iostream>
#include <stdlib.h>

#include "VAstronometry.h"

using namespace std;

int main( int argc, char* argv[] )
{
    if( argc != 6 && argc != 7 && argc != 8 )
    {
        cout << "printBinaryOrbitalPhase <year> <month> <day> <MJD> <Object> [orbit] [T0]" << endl;
        cout << endl;
        cout << "   set either year/month/day or MJD (set other values to -1)" << endl;
        cout << endl;
        cout << "   print orbit for a whole year: set month and day to -1 each" << endl;
        cout << endl;
        cout << "   <Object=\"LS I +61 303\", \"HESS J0632+057\", \"LS 5039\", \"PSR B1259–63\">" << endl;
        cout << endl;
        cout << "   [orbit]   orbital period in days (optional)" << endl;
        cout << "   [T0]      MJD 0" << endl;
        cout << endl;
        exit( 0 ) ;
    }
    int year = atoi( argv[1] );
    int month = atoi( argv[2] );
    int day = atoi( argv[3] );
    double djm = atoi( argv[4] );
    int j = 0;
    double fd = 0.;
    if( year < 0 && month < 0 && day < 0 && djm > 0 )
    {
        VAstronometry::vlaDjcl( djm, &year, &month, &day, &fd, &j );
    }
    bool bPrintWholeYear = false;
    if( month < 0 && day < 0 && djm < 0 )
    {
        month = 1;
        day = 1;
        bPrintWholeYear = true;
    }
    
    string object = argv[5];
    double orbit = 0.;
    double t0 = 0.;
    // LS I+61 303
    if( object == "LS I +61 303" )
    {
        orbit = 26.4960;
        t0    = 2443366.775 - 2400000.5;
    }
    // HESS J0632+057
    else if( object == "HESS J0632+057" )
    {
        orbit = 315.;
        t0    = 54857.;
    }
    // LS 5039
    else if( object == "LS 5039" )
    {
        orbit = 3.90603;
        t0    = 51942.59;
    }
    // PSR B1259–63 (needs to be verified)
    else if( object == "PSR B1259–63" )
    {
        orbit = 1237.;
        t0    = 53040.22;
    }
    else if( object == "IGR J00370+6122" )
    {
        orbit = 15.6627;
        t0    = 2450000 + 3000.04 - 2400000.5;
    }
    else if( object == "V662 Cas" )
    {
        orbit = 11.5983;
        t0    = 2400000 + 51825.3 - 2400000.5;
    }
    // user given timing parameters
    if( argc > 6 )
    {
        orbit = atof( argv[6] );
    }
    if( argc == 8 )
    {
        t0 = atof( argv[7] );
    }
    
    cout << "Orbital parameters for " << object << ": ";
    cout << "MJD_0 = " << t0;
    cout << ", Orbit [d]: " << orbit;
    cout << endl << endl;
    
    // get mjd for given date
    VAstronometry::vlaCldj( year, month, day, &djm, &j );
    if( j == 0 )
    {
        if( !bPrintWholeYear )
        {
            // calculate orbital phase
            cout << "Orbital phase of " << object << " for " << year << "-" << month << "-" << day << ", MJD " << djm << ": ";
            cout << endl;
            cout << ( djm - t0 ) / orbit - ( int )( ( djm - t0 ) / orbit ) << endl;
            double iP = ( djm - t0 ) / orbit - ( int )( ( djm - t0 ) / orbit );
            if( iP < 0. )
            {
                cout << 1. + iP << endl;
            }
            else
            {
                cout << iP << endl;
            }
        }
        else
        {
            float p = 0.;
            float p_old = p;
            double djm_old = djm;
            // 10 day interval
            for( int i = 0; i < 366; i++ )
            {
                p = ( djm  - t0 ) / orbit - ( int )( ( djm - t0 ) / orbit );
                if( p < p_old )
                {
                    VAstronometry::vlaDjcl( djm_old, &year, &month, &day, &fd, &j );
                    cout << "new orbit (MJD: " << djm_old << "): \t";
                    cout << year << "/" << month << "/" << day << endl;
                }
                VAstronometry::vlaDjcl( djm, &year, &month, &day, &fd, &j );
                cout << "Orbital phase of " << object << " for " << year << "-" << month << "-" << day << ", MJD " << djm << ": ";
                cout << ( djm - t0 ) / orbit - ( int )( ( djm - t0 ) / orbit ) << endl;
                p_old = p;
                djm_old = djm;
                djm++;
            }
        }
    }
    else
    {
        cout << "failed calculation of MJD" << endl;
        exit( -1 );
    }
    cout << endl;
}
