/*! \class VRunList
    \brief data class


*/

#include "VRunList.h"

#ifndef VRUNLIST_CPP
#define VRUNLIST_CPP


VRunList::VRunList()
{
    reset();
}


void VRunList::reset()
{
    runnumber = 0;
    MJD = 0.;
    tOn = 0.;
    deadTimeFraction = 0.;
    NOn = 0.;
    NOff = 0.;
    OffNorm = 0.;
    elevationOn = 0.;
    elevationOff = 0.;
    pedvarsOn = 0.;
    alpha = 1.;
    energyThreshold = 0.;
    TargetRAJ2000 = 0.;
    TargetDecJ2000 = 0.;
    phase = -1.e9;
}


void VRunList::print()
{
    cout << "Run: " << runnumber;
    cout << setprecision( 3 ) << " (MJD " << MJD << ")";
    cout << " tOn [min]: " << setprecision( 2 ) << fixed << tOn / 60.;
    cout << " (dead time " << setprecision( 2 ) << ( 1. - deadTimeFraction ) * 100. << "%) ";
    cout << " Non/NOff: " << setw( 8 ) << fixed << NOn << "/" << NOff;
    cout << setprecision( 3 ) << setw( 8 ) << " alpha: " << alpha;
    cout << " elevation [deg]: " << setw( 8 );
    if( elevationOn > 0. )
    {
        cout << elevationOn;
    }
    else
    {
        cout << elevationOff << " (OFF) ";
    }
    cout << " mean pedvar: " << pedvarsOn;
    if( energyThreshold > 0. )
    {
        cout << " energy threshold [GeV]: " << energyThreshold * 1.e3;
    }
    if( phase > 0. )
    {
        cout << " phase: " << phase;
    }
    cout << endl;
}

void VRunList::print( bool csv )
{
    if( csv )
    {
        cout << runnumber << ", " << setprecision( 6 ) << MJD;
        cout << ", " << setprecision( 2 ) << fixed << ( 1. - deadTimeFraction ) * 100.;
        cout << ", " << fixed << NOn << ", " << NOff;
        cout << setprecision( 3 ) << ", " << alpha << ", ";
        if( elevationOn > 0. )
        {
            cout << elevationOn;
        }
        else
        {
            cout << elevationOff;
        }
        cout << ", " << pedvarsOn;
        if( energyThreshold > 0. )
        {
            cout << "," << energyThreshold * 1.e3;
        }
        if( phase > 0. )
        {
            cout << ", " << phase;
        }
        cout << endl;
    }
    else
    {
        print();
    }
}

#endif
