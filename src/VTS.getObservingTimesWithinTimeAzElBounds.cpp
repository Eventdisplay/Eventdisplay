///////////////////////////////
//  Take a single runnumber, and print the elevation and azimuth
//  the telescope is pointing at once per second
#include "CData.h"
#include "VSkyCoordinates.h"
#include "VAstronometry.h"
#include "VDB_Connection.h"
#include "VGlobalRunParameter.h"
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <bitset>
// #include <algorithm>
// #include <functional>
// #include <cctype>
// #include <locale>
void getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip );

bool is_number( const std::string& s )
{
    std::string::const_iterator it = s.begin();
    while( it != s.end() && std::isdigit( *it ) )
    {
        ++it;
    }
    return !s.empty() && it == s.end() ;
}

using namespace std;

int main( int argc, char* argv[] )
{
    // parse argument
    if( argc != 7 )
    {
        cout << "VTS.getTelescopePointingWithinBounds:" << endl;
        cout << "  Find the blocks of time where the telescopes were within azimuth and elevation bounds." << endl;
        cout << endl;
        cout << "For Example:" << endl;
        cout << "  $ VTS.getTelescopePointingWithinBounds 2011-01-01 2015-01-01 70         90         210        280" << endl;
        cout << "  $ VTS.getTelescopePointingWithinBounds <min date> <max date> <min elev> <max elev> <min azim> <max azim>" << endl;
        cout << endl;
        cout << "Will output pairs of MJD times when the telescopes were within the given orientation bounds." << endl;
        cout << "  $ VTS.getTelescopePointingWithinBounds 2011-01-01 2015-01-01 70         90         210        280" << endl;
        cout << "  57065.5314298 57066.5314298" << endl;
        cout << "  57067.5314298 57069.5314298" << endl;
        cout << "  57071.5314298 57072.5314298" << endl;
        cout << endl;
        cout << "Arguments can be replaced with '-' to use its default value:" << endl;
        cout << "   default azimuth range  : 0-360" << endl;
        cout << "   default elevation range: 0-90" << endl;
        cout << "   default time range: 2011-01-01 to today" << endl;
        cout << endl;
        cout << "* there will be one line per ~2 seconds of observatory operation" << endl;
        cout << "* elevation is degrees above the local horizon" << endl;
        cout << "* azimuth is degrees from North(0deg) to East(90deg)" << endl;
        return 1 ;
    }
    
    int  yearstart = 2011 ;
    int  monthstart = 1    ;
    int  daystart  = 1    ;
    int  yearend   = 2015 ;
    int  monthend  = 1    ;
    int  dayend    = 1    ;
    int  azimlow   = 0    ;
    int  azimhigh  = 360  ;
    int  elevlow   = 0    ;
    int  elevhigh  = 90   ;
    if( strcmp( argv[1], "-" ) == 0 )
    {
        yearstart = 2011 ;
        monthstart = 1    ;
        daystart  = 1    ;
    }
    else
    {
        sscanf( argv[1], "%d-%d-%d", &yearstart, &monthstart, &daystart ) ;
    }
    if( strcmp( argv[2], "-" ) == 0 )
    {
        yearend = 2011 ;
        monthend = 1    ;
        dayend  = 1    ;
    }
    else
    {
        sscanf( argv[2], "%d-%d-%d", &yearend  , &monthend  , &dayend ) ;
    }
    
    if( strcmp( argv[3], "-" ) == 0 )
    {
        elevlow  = 0               ;
    }
    else
    {
        elevlow  = atoi( argv[3] ) ;
    }
    if( strcmp( argv[4], "-" ) == 0 )
    {
        elevhigh = 90              ;
    }
    else
    {
        elevhigh = atoi( argv[4] ) ;
    }
    if( strcmp( argv[5], "-" ) == 0 )
    {
        azimlow  = 0               ;
    }
    else
    {
        azimlow  = atoi( argv[5] ) ;
    }
    if( strcmp( argv[6], "-" ) == 0 )
    {
        azimhigh = 360             ;
    }
    else
    {
        azimhigh = atoi( argv[6] ) ;
    }
    
    printf( "read in arguments:\n" ) ;
    printf( "  starting date: %d-%d-%d\n", yearstart, monthstart, daystart ) ;
    printf( "  ending date  : %d-%d-%d\n", yearend  , monthend  , dayend ) ;
    printf( "  elev range: %d - %d\n", elevlow, elevhigh ) ;
    printf( "  azim range: %d - %d\n", azimlow, azimhigh ) ;
    
    VGlobalRunParameter* blah = new VGlobalRunParameter() ;
    //cout << " VGlobal->getDBServer(): " << blah->getDBServer() << endl;
    // start connection
    // fixed below to follow the party line
    string tmpdb_ver;
    tmpdb_ver = blah->getDBServer() + "/VERITAS" ;
    cout << "Using server " << tmpdb_ver << endl;
    VDB_Connection my_connection_ver( tmpdb_ver.c_str(), "readonly", "" ) ;
    if( !my_connection_ver.Get_Connection_Status() )
    {
        cout << "error connecting to db" << tmpdb_ver << endl;
        return -1;
    }
    string tmpdb_vof;
    tmpdb_vof = blah->getDBServer() + "/VOFFLINE" ;
    cout << "Using server " << tmpdb_vof << endl;
    VDB_Connection my_connection_vof( tmpdb_vof.c_str(), "readonly", "" ) ;
    if( !my_connection_vof.Get_Connection_Status() )
    {
        cout << "error connecting to db:" << tmpdb_vof << endl;
        return -1;
    }
    
    // calculate mjd's
    char   mjdstart_char[100] = ""  ;
    char   mjdend_char[100]   = ""  ;
    sprintf( mjdstart_char, "%4d-%02d-%02d 00:00:00", yearstart, monthstart, daystart ) ;
    sprintf( mjdend_char  , "%4d-%02d-%02d 00:00:00", yearend  , monthend  , dayend ) ;
    string mjdstart_str = mjdstart_char ;
    string mjdend_str   = mjdend_char   ;
    //printf("init:mjdstart:'%s'\n", mjdstart_str.c_str() ) ;
    //printf("init:mjdend  :'%s'\n", mjdend_str.c_str()   ) ;
    int mjdstart_day    = 0 ;
    int mjdend_day      = 0 ;
    double mjdstart_time = 0.0 ;
    double mjdend_time   = 0.0 ;
    getDBMJDTime( mjdstart_str, mjdstart_day, mjdstart_time, true );
    getDBMJDTime( mjdend_str  , mjdend_day  , mjdend_time  , true );
    double mjdstart = mjdstart_day + ( mjdstart_time / 86400.0 ) ;
    double mjdend   = mjdend_day   + ( mjdend_time   / 86400.0 ) ;
    //printf("mjdstart:%f\n", mjdstart ) ;
    //printf("mjdend  :%f\n", mjdend   ) ;
    
    // construct and submit query
    printf( "submitting database query, please wait....\n" );
    char c_query[1000] ;
    sprintf( c_query, "select mjd, ra, decl from tblPointing_Monitor_Telescope2_Calibrated_Pointing where mjd > %f and mjd < %f", mjdstart, mjdend ) ;
    if( !my_connection_vof.make_query( c_query ) )
    {
        cout << "Error, unable to get response from db for query \"" << c_query << "\"" << endl;
        cout << "Exiting..." << endl;
        return -1;
    }
    
    // loop over rows in response
    string itemp ;
    TSQLResult* db_res2 = my_connection_vof.Get_QueryResult();
    int fNRows = db_res2->GetRowCount();
    double imjd, ira, idecl ;
    vector<double> mjd  ;
    vector<int>    mjdDay  ;
    vector<double> mjdDayFraction  ;
    vector<double> mjdSecondsOfDay  ;
    vector<double> ra   ;
    vector<double> decl ;
    printf( "received %d database rows from tblPointing_Monitor_Telescope2_Calibrated_Pointing...\n", fNRows ) ;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res2->Next();
        if( !db_row )
        {
            break;
        }
        
        // time
        itemp = db_row->GetField( 0 ) ;
        sscanf( itemp.c_str(), "%lf", &imjd ) ;
        mjd.push_back( imjd );
        mjdDay.push_back( ( int )imjd ) ;
        mjdDayFraction.push_back( ( double )( imjd - ( int )imjd ) ) ;
        mjdSecondsOfDay.push_back( mjdDayFraction.back() * 86400.0 ) ;
        
        // ra
        itemp = db_row->GetField( 1 ) ;
        sscanf( itemp.c_str(), "%lf", &ira ) ;
        ra.push_back( ira );
        
        // decl
        itemp = db_row->GetField( 2 ) ;
        sscanf( itemp.c_str(), "%lf", &idecl ) ;
        decl.push_back( idecl );
        
    }
    
    // convert pointing data from ra/dec to elev/azi
    vector<double> selectedMJD ;
    double tmpelev = 0.0 ;
    double tmpazim = 0.0 ;
    int telescope  = 0   ; // 0-3 for T1-4
    float lat = blah->getObservatory_Latitude_deg() ;
    float lon = blah->getObservatory_Longitude_deg() ;
    VSkyCoordinates* vsky = new VSkyCoordinates() ;
    vsky->supressStdoutText( true ) ;
    vsky->setObservatory( lon, lat );
    
    vector<double> timezone_start;
    vector<double> timezone_end  ;
    bool withinBoundsLastRow = false ;
    bool withinBounds        = false ;
    double tdiff   = 0.0 ;
    int zone = 0 ;
    for( unsigned int i_row = 0 ; i_row < decl.size() ; i_row++ )
    {
        withinBounds = false ;
        // target's ra and decl in degrees
        vsky->setTargetJ2000( decl[i_row] * TMath::RadToDeg() , ra[i_row] * TMath::RadToDeg() ) ;
        // day that you're looking for elev and azimuth on
        vsky->precessTarget( mjdDay[i_row], telescope ) ;
        // calculate new param
        vsky->updatePointing( mjdDay[i_row], mjdSecondsOfDay[i_row] ) ;
        tmpelev = vsky->getTargetElevation() ;
        tmpazim = vsky->getTargetAzimuth()   ;
        
        // see if row is within elevation bounds
        if( tmpelev > elevlow && tmpelev < elevhigh )
        {
            // check to see how our azimuth bounds are oriented
            //  40-70  = span of 30  deg
            //  40-300 = span of 260 deg
            //  330-20 = span of 50  deg, and must be careful
            if( azimlow < azimhigh )
            {
                // then span is a normal one
                if( tmpazim > azimlow && tmpazim < azimhigh )
                {
                    withinBounds = true ;
                }
            }
            else
            {
                // then span crosses deg=0, and we need to be careful
                if( tmpazim > azimlow || tmpazim < azimhigh )
                {
                    withinBounds = true ;
                }
            }
        }
        
        // go through the logic of dealing with each row
        // deal with only getting 0 rows
        if( decl.size() == 0 )
        {
            printf( "no runs found, try wider bounds.\n" );
            return 1;
        }
        // deal with only getting 1 row, assume a time span of 2 seconds, centered on our only row
        else if( decl.size() == 1 )
        {
            zone += 1 ;
            printf( "Z%2d %12.6f only one row!\n", zone, mjd[i_row] );
            timezone_start.push_back( selectedMJD[0] - 0.00001157 ) ;
            timezone_end.push_back( selectedMJD[0] + 0.00001157 ) ;
        }
        else
        {
            // if its the first row and we're within the bounds, start a new time zone
            if( i_row == 0 )
            {
                if( withinBounds )
                {
                    zone += 1 ;
                    timezone_start.push_back( mjd[i_row] ) ;
                    printf( "Z%2d %12.6f new timezone - first row\n", zone, mjd[i_row] ) ;
                }
            }
            
            // if its the last row, close off the current timezone
            else if( i_row == mjd.size() - 1 )
            {
                if( withinBounds )
                {
                    timezone_end.push_back( mjd[i_row] ) ;
                    zone += 1 ;
                    printf( "Z%2d %12.6f closing timezone - last row\n"       , zone, mjd[i_row] ) ;
                }
            }
            
            // if this isn't the first or last rows...
            else
            {
                // check if we've entered a new time zone
                if( withinBounds == true  && withinBoundsLastRow == false )
                {
                    timezone_start.push_back( mjd[i_row] ) ;
                    zone += 1 ;
                    printf( "Z%2d %12.6f new timezone     - back within bounds\n", zone, mjd[i_row] );
                }
                // check if we've left a time zone
                else if( withinBounds == false && withinBoundsLastRow == true )
                {
                    timezone_end.push_back( mjd[i_row] ) ;
                    printf( "Z%2d %12.6f leaving timezone - leaving bounds\n", zone, mjd[i_row] ) ;
                }
                // check if the inter-row time difference is large enough to warrent a new timezone
                else if( withinBounds == true && withinBoundsLastRow == true )
                {
                    tdiff = mjd[i_row] - mjd[i_row - 1] ;
                    // check if the inter-row distance is high enough to warrent a new timezone
                    if( tdiff > 0.000035 )
                    {
                        printf( "Z%2d %12.6f leaving timezone - large time jump\n", zone, mjd[i_row] ) ;
                        zone += 1 ;
                        printf( "Z%2d %12.6f new timezone     - large time jump\n", zone, mjd[i_row] ) ;
                        timezone_end.push_back( mjd[i_row  ] );
                        timezone_start.push_back( mjd[i_row - 1] );
                    }
                }
            }
        }
        
        withinBoundsLastRow = withinBounds ;
    }
    
    printf( "found %d/%d time segments within bounds...\n", ( int )timezone_start.size(), ( int )timezone_end.size() ) ;
    if( timezone_start.size() != timezone_end.size() )
    {
        printf( "\nWarning, timezone vectors are not the same size!\n\n" ) ;
    }
    for( unsigned int i_row = 0 ; i_row < timezone_start.size() ; i_row++ )
    {
        printf( "MJDSEGMENT  %9f - %9f\n", timezone_start[i_row], timezone_end[i_row] ) ;
    }
    
    return 0 ;
}


void getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip )
{
    if( itemp.size() < 16 )
    {
        MJD = 0;
        Time = 0.;
        return;
    }
    if( bStrip )
    {
        itemp.replace( itemp.find( "-" ), 1, "" );
        itemp.replace( itemp.find( "-" ), 1, "" );
        itemp.replace( itemp.find( " " ), 1, "" );
        itemp.replace( itemp.find( ":" ), 1, "" );
        itemp.replace( itemp.find( ":" ), 1, "" );
    }
    int y, m, d, h, min, s, ms, l;
    double gMJD;
    
    // get y, m, d
    y = atoi( itemp.substr( 0, 4 ).c_str() );
    m = atoi( itemp.substr( 4, 2 ).c_str() );
    d = atoi( itemp.substr( 6, 2 ).c_str() );
    h = atoi( itemp.substr( 8, 2 ).c_str() );
    min = atoi( itemp.substr( 10, 2 ).c_str() );
    s = atoi( itemp.substr( 12, 2 ).c_str() );
    if( !bStrip )
    {
        ms = atoi( itemp.substr( 14, 3 ).c_str() );
    }
    else
    {
        ms = 0;
    }
    
    // calculate MJD
    VAstronometry::vlaCldj( y, m, d, &gMJD, &l );
    MJD = ( int )gMJD;
    Time = h * 60.*60. + min * 60. + s + ms / 1.e3;
}

