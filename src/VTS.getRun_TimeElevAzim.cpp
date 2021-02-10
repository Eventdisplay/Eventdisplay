///////////////////////////////
//  Take a single runnumber, and print the elevation and azimuth
//  the telescope is pointing at once per second
#include "CData.h"
#include "VAstronometry.h"
#include "VSkyCoordinates.h"
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

// remove any spaces at the beginning of the line (LeftTrim)
static inline std::string& ltrim( std::string& s )
{
    s.erase( s.begin(), std::find_if( s.begin(), s.end(), std::not1( std::ptr_fun<int, int>( std::isspace ) ) ) ) ;
    return s ;
}
static inline std::string& rchop( std::string& s )
{
    if( s.find_first_of( " " ) != std::string::npos )
    {
        s.erase( s.begin() + s.find_first_of( " " ), s.end() ) ;
    }
    return s ;
}
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
    if( argc != 3 )
    {
        cout << "Error, needs runnumber or runlist argument." << endl;
        cout << "For Example:" << endl;
        cout << "  $ VTS.getRun_TimeElevAzim -r 65765" << endl;
        cout << " or " << endl;
        cout << "  $ VTS.getRun_TimeElevAzim -f mysimplerunlist.dat" << endl;
        cout << endl;
        cout << "Outputdata will be in the format:" << endl;
        cout << "RUN##### MJD(days) Elevation(degs) Azimuth(degs)" << endl;
        cout << "RUN##### MJD(days) Elevation(degs) Azimuth(degs)" << endl;
        cout << "RUN##### MJD(days) Elevation(degs) Azimuth(degs)" << endl;
        cout << "RUN##### MJD(days) Elevation(degs) Azimuth(degs)" << endl;
        cout << endl;
        cout << "* there will be one line per ~2 seconds of observatory operation" << endl;
        cout << "* elevation is degrees above the local horizon" << endl;
        cout << "* azimuth is degrees from North(0deg) to East(90deg)" << endl;
        cout << "* to only get the data (and strip out the other text), do:" << endl;
        cout << "  $ VTS.getRun_TimeElevAzim -f mysimplerunlist.dat | grep -P \"^RUN\\d{5}\"" << endl;
        cout << "* or a specific run's data with:" << endl;
        cout << "  $ VTS.getRun_TimeElevAzim -f mysimplerunlist.dat | grep -P \"^RUN66565\"" << endl;
        cout << endl;
        return 1 ;
    }
    
    int IRUN  = 1;  // 1 is for runnumber
    int IFILE = 2; // 2 is for simple list of runnumbers
    //cout << "IRUN:  " << IRUN << endl;
    //cout << "IFILE: " << IFILE << endl;
    int  inputmode     = 0 ;
    int  c             = 0 ;
    int  inputrun      = 0 ;
    string inputfile ;
    while( ( c = getopt( argc, argv, "r:f:" ) ) != -1 )
    {
        switch( c )
        {
            case 'r':
                //cout << "-r detected" << endl;
                inputmode = IRUN ;
                cout << "optarg: " << optarg << endl;
                //if ( ! is_number( optarg ) ) {
                //	cout << "Error, argument for -r flag must be an integer, exiting..." << endl;
                //	return 1 ;
                //}
                inputrun  = atoi( optarg ) ;
                break ;
            case 'f':
                //cout << "-f detected" << endl;
                inputmode = IFILE ;
                inputfile = optarg ;
                cout << "Using runlist " << inputfile << " ..." << endl;
                break ;
            default:
                printf( "unknown flag '%c', exiting!\n", c ) ;
                return 1 ;
                break ;
        }
    }
    //cout << "inputmode: " << inputmode << endl;
    //cout << "inputrun : " << inputrun  << endl;
    //cout << "inputfile: " << inputfile << endl;
    //cout << endl;
    
    int runnumber = 0 ;
    vector<int> runlist ;
    if( inputmode == IRUN )
    {
        runlist.push_back( inputrun ) ;
    }
    if( inputmode == IFILE )
    {
        ifstream fin( inputfile.c_str() ) ;
        string   file_line ;
        if( fin.good() )
        {
            while( std::getline( fin, file_line ) )
            {
                //cout << "line:    '" << file_line << "'" << endl;
                ltrim( file_line ) ; // clean up the line
                //cout << "    ltrim:'" << file_line << "'" << endl;
                rchop( file_line ) ; // clean up the line
                //cout << "    rchop:'" << file_line << "'" << endl;
                if( is_number( file_line ) )
                {
                    runnumber = atoi( file_line.c_str() ) ;
                    if( runnumber < 9999 or runnumber > 99999 )
                    {
                        cout << "Error, '" << file_line << "' needs to be a valid runnumber (9999-99999)" << endl;
                        return 1 ;
                    }
                    else
                    {
                        runlist.push_back( runnumber ) ;
                    }
                }
                else
                {
                    cout << "Error, '" << file_line << "' needs to be a valid runnumber (9999-99999)" << endl;
                    return 1 ;
                }
            }
        }
        else
        {
            cout << "Error, file '" << inputfile << "' could not be read, exiting..." << endl;
            return 1;
        }
    }
    //cout << "runlist.size(): " << runlist.size() << endl;
    
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
        cout << "error connecting to db" << endl;
        return -1;
    }
    
    for( unsigned int i_run = 0 ; i_run < runlist.size() ; i_run++ )
    {
        cout << "Downloading observatory data for run " << runlist[i_run] << "..." << endl;
        runnumber = runlist[i_run] ;
        
        // get run start and finish, # of telescopes
        char c_query[1000] ;
        sprintf( c_query, "select run_id, data_start_time, data_end_time, config_mask from tblRun_Info where run_id = %d", runnumber ) ;
        //printf( "%s\n", query ) ;
        if( !my_connection_ver.make_query( c_query ) )
        {
            cout << "Error, unable to get response from db for query \"" << c_query << "\"" << endl;
            cout << "Exiting..." << endl;
            return -1;
        }
        
        //TSQLResult *db_res = f_db->Query( query );
        TSQLResult* db_res = my_connection_ver.Get_QueryResult();
        
        if( !db_res )
        {
            cout << "TSQLResult empty, exiting..." << endl ;
            exit( 1 ) ;
        }
        int fNRows = db_res->GetRowCount();
        string itemp;
        double MJDStart = 0.0, MJDEnd = 0.0 ;
        unsigned long int ImgSel = 0 ;
        for( int j = 0; j < fNRows; j++ )
        {
            TSQLRow* db_row = db_res->Next();
            if( !db_row )
            {
                break;
            }
            int iMJD = 0;
            double iTime = 0.;
            
            itemp = db_row->GetField( 1 );
            getDBMJDTime( itemp, iMJD, iTime, true );
            MJDStart = ( double )iMJD + iTime / 86400.0;
            //printf("date to MJD: %s : %f : sec %f\n", itemp.c_str(), MJDStart, (MJDStart-(int)MJDStart)*86400.0 ) ;
            
            itemp = db_row->GetField( 2 );
            getDBMJDTime( itemp, iMJD, iTime, true );
            MJDEnd = ( double )iMJD + iTime / 86400.0;
            //printf("date to MJD: %s : %f : sec %f\n", itemp.c_str(), MJDEnd, (MJDEnd-(int)MJDEnd)*86400.0 ) ;
            
            ImgSel = ( unsigned long int )db_row->GetField( 3 ) ;
        }
        
        //printf( "Run %d started at %f and ended at %f, using telescope combo %lu\n", runnumber, MJDStart, MJDEnd, ImgSel ) ;
        
        string tmpdb_vof;
        tmpdb_vof = blah->getDBServer() + "/VOFFLINE" ;
        VDB_Connection my_connection_vof( tmpdb_ver.c_str(), "readonly", "" ) ;
        if( !my_connection_vof.Get_Connection_Status() )
        {
            cout << "error connecting to db VOFFLINE" << endl;
            return -1;
        }
        
        // request pointing data
        vector <double> mjd  ;
        vector <int>    mjdDay ;
        vector <double> mjdDayFraction ;
        vector <double> mjdSecondsOfDay ;
        vector <double> ra   ;
        vector <double> decl ;
        bitset<8 * sizeof( ULong64_t )> a = ImgSel;
        int telescope = -1 ;
        for( unsigned int i_tel = 0 ; i_tel < 4 ; i_tel++ )
        {
            if( ! a.test( i_tel ) )
            {
                continue ;
            }
            sprintf( c_query, "select mjd, ra, decl from VOFFLINE.tblPointing_Monitor_Telescope%d_Calibrated_Pointing where mjd > %f and mjd < %f ;", i_tel, MJDStart, MJDEnd ) ;
            //printf( "%s\n", query) ;
            if( !my_connection_vof.make_query( c_query ) )
            {
                cout << "Error, unable to get response from db for query \"" << c_query << "\"" << endl;
                cout << "Exiting..." << endl;
                return -1;
            }
            
            //TSQLResult *db_res2 = f_db->Query( query );
            TSQLResult* db_res2 = my_connection_vof.Get_QueryResult();
            fNRows = db_res2->GetRowCount();
            //int count=0 ;
            //printf( "%d rows\n", fNRows ) ;
            double imjd, ira, idecl ;
            for( int j = 0; j < fNRows; j++ )
            {
                TSQLRow* db_row = db_res2->Next();
                if( !db_row )
                {
                    break;
                }
                //if ( count<8 ) { count++ ; printf( "Row %d  :  %s  :  %s  :  %s\n", j, db_row->GetField(0), db_row->GetField(1), db_row->GetField(2) ) ; }
                
                // time
                itemp = db_row->GetField( 0 ) ;
                sscanf( itemp.c_str(), "%lf", &imjd ) ;
                mjd.push_back( imjd );
                mjdDay.push_back( ( int )imjd ) ;
                mjdDayFraction.push_back( ( double )( imjd - ( int )imjd ) ) ;
                mjdSecondsOfDay.push_back( mjdDayFraction.back() * 86400.0 ) ;
                
                // ra decl
                itemp = db_row->GetField( 1 ) ;
                sscanf( itemp.c_str(), "%lf", &ira ) ;
                ra.push_back( ira );
                itemp = db_row->GetField( 2 ) ;
                sscanf( itemp.c_str(), "%lf", &idecl ) ;
                decl.push_back( idecl );
                //if ( count<8 ) { printf("     imjd:%.8f   ra:%f, decl:%f\n", imjd, ira, idecl ) ; }
                
            }
            
            // stop after first telescope, we only need one
            if( a.test( i_tel ) && fNRows > 0 )
            {
                telescope = i_tel ;
                break ;
            }
        }
        //f_db->Close() ;
        
        // convert pointing data from ra/dec to elev/azi
        //float lat = 31.675  ; // observatory coordinates
        float lat = blah->getObservatory_Latitude_deg() ;
        //float lon = 110.952 ;
        float lon = blah->getObservatory_Longitude_deg() ;
        VSkyCoordinates* vsky = new VSkyCoordinates() ;
        vsky->supressStdoutText( true ) ;
        vsky->setObservatory( lon, lat );
        for( unsigned int i_row = 0 ; i_row < decl.size() ; i_row++ )
        {
            // target's ra and decl in degrees
            vsky->setTargetJ2000( decl[i_row] * TMath::RadToDeg() , ra[i_row] * TMath::RadToDeg() ) ;
            // day that you're looking for elev and azimuth on
            vsky->precessTarget( mjdDay[i_row], telescope ) ;
            // calculate new param
            vsky->updatePointing( mjdDay[i_row], mjdSecondsOfDay[i_row] ) ;
            printf( "RUN%d %f %f %f\n", runnumber, mjd[i_row], vsky->getTargetElevation(), vsky->getTargetAzimuth() ) ;
            // not sure why mode=1 gives weird answers, mode=2 seems to match up with the values in the database-generated log
            // e.g. veritasm.sao.arizona.edu/DQM/cgi-bin/loggen/query_night?search=runid&format=html&dqm=1&dqm_inst=&4unid=54679
        }
        
        
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

