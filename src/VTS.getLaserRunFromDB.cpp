/* VTS.getLaserRunFromDB.cpp  get laser run for a data run from DB

*/


#include <iostream>
#include <string>
#include <vector>

#include "VDBRunInfo.h"
#include "VGlobalRunParameter.h"

using namespace std;


int main( int argc, char* argv[] )
{
    if( argc != 3 )
    {
        cout << endl;
        cout << "./VTS.getLaserRunFromDB <telID> <data run number>" << endl;
        cout << endl;
        cout << "   read corresponding laser/flaser run number for a data run number from VTS data base" << endl;
        cout << "   (telescope ID start at 1 for T1)" << endl;
        cout << endl;
        exit( 0 );
    }
    int fTelID = atoi( argv[1] );
    int fRunNumber = atoi( argv[2] );
    const unsigned int fNTel = 4;
    if( fTelID > ( int )fNTel )
    {
        cout << "error: maximum number of telescopes are " << fNTel;
        cout << ", requested is laser run for telescope telescope " << fTelID << endl;
        exit( -1 );
    }
    else if( fTelID < 1 )
    {
        cout << "error: first telescope is T1 (1)" << endl;
        exit( -1 );
    }
    
    // get DBservers
    VGlobalRunParameter* fRunPara = new VGlobalRunParameter();
    
    // DB run info
    VDBRunInfo i_DBinfo( fRunNumber, fRunPara->getDBServer(), fNTel );
    if( i_DBinfo.isGood() )
    {
        vector< unsigned int > iL = i_DBinfo.getLaserRun();
        if( fTelID - 1 < ( int )iL.size() )
        {
            cout << iL[fTelID - 1] << endl;
        }
        else
        {
            cout << "0" << endl;
        }
    }
}

