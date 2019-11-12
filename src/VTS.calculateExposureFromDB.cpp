/*! file VTS.calculateExposureFromDB read exposure from DB

    (VERITAS only)


*/

#include "VGlobalRunParameter.h"
#include "VExposure.h"

#include <iostream>
#include <string>

int main( int argc, char* argv[] )
{
    cout << endl;
    cout << "VTS.calculateExposureFromDB (" << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << endl;
    if( argc != 4 )
    {
        cout << "./VTS.calculateExposureFromDB <db start time> <db stop time> <output root file>" << endl;
        cout << endl;
        cout << "example: VTS.calculateExposureFromDB 2009-01-09 2009-01-31 tt.root" << endl;
        cout << endl;
        exit( 0 );
    }
    
    string iA = argv[1];
    string iB = argv[2];
    string iO = argv[3];
    
    VExposure a;
    a.setTimeRange( iA, iB );
    a.readFromDB();
    a.writeRootFile( iO );
}
