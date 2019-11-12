/**
 * \file
 * \brief Function definitions for Pointing Monitor interactions with the database [ VDB v4.1.0 (Jan 2010), adapted for TSQLServer (Oct 2010) ]
 *
 */

#include "PointingMonitor.h"

using namespace pointingmonitor;

PointingMonitor::PointingMonitor()
{

}

// local function to convert a date string to an integer by removing all non-number characters
uint32_t PointingMonitor::dateToUInt32( const string& dateStr )
{
    stringstream strBuf;
    for( uint32_t i = 0; i < dateStr.size(); ++i )
    {
        char c = dateStr[i];
        if( c >= '0' && c <= '9' )
        {
            strBuf << c;
        }
    }
    uint32_t dateInt = 0;
    strBuf >> dateInt;
    return dateInt;
}

/********************************************************************************/
/* Functions for converting between Modified Julian Date and calender date/time */
/********************************************************************************/

uint32_t PointingMonitor::mjdToDate( double mjd )
{
    time_t timeStamp = ( time_t )( mjd - 40587 ) * 86400;
    tm* timeStructure = gmtime( &timeStamp );
    uint64_t timeInt = timeStructure->tm_year + 1900;
    timeInt *= 100;
    timeInt += timeStructure->tm_mon + 1;
    timeInt *= 100;
    timeInt += timeStructure->tm_mday;
    return timeInt;
}

double PointingMonitor::dateTimeToMjd( uint64_t dateTime )
{
    tm timeStructure;
    timeStructure.tm_isdst = 0;
    timeStructure.tm_sec = dateTime % 100;
    dateTime /= 100;
    timeStructure.tm_min = dateTime % 100;
    dateTime /= 100;
    timeStructure.tm_hour = dateTime % 100;
    dateTime /= 100;
    timeStructure.tm_mday = dateTime % 100;
    dateTime /= 100;
    timeStructure.tm_mon = dateTime % 100 - 1;
    dateTime /= 100;
    timeStructure.tm_year = dateTime - 1900;
    time_t timeStamp = timegm( &timeStructure );
    double mjd = double( timeStamp ) / 86400 + 40587;
    return mjd;
}

/************************************************************************/
/* Functions for writing/reading camera parameters to/from the database */
/************************************************************************/

vector<pointingmonitor::CameraParameters> PointingMonitor::getCameraParametersList( uint32_t telescope_id, uint32_t date, uint32_t limit )
{
    vector<pointingmonitor::CameraParameters> parametersVec;
    
    // create database command string
    ostringstream strbuf;
    strbuf << "select * from tblPointing_Monitor_Camera_Parameters";
    strbuf << " where telescope_id=" << telescope_id;
    if( date != 0 )
    {
        strbuf << " and start_date<=" << date;
        strbuf << " and end_date>=" << date;
    }
    strbuf << " order by version desc";
    if( limit != 0 )
    {
        strbuf << " limit " << limit;
    }
    
    string fDBserver = getDBServer() + "VOFFLINE";
    //std::cout<<"PointingMonitor::getCameraParametersList "<<std::endl;
    VDB_Connection my_connection( fDBserver.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "PointingMonitor: failed to connect to database server" << endl;
        return parametersVec; // JG
    }
    if( !my_connection.make_query( strbuf.str().c_str() ) )
    {
        return parametersVec;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    int fNRows = db_res->GetRowCount();
    
    for( int i = 0; i < fNRows; i++ )
    {
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            continue;
        }
        
        // add new element to result vector and fill with numbers
        parametersVec.push_back( pointingmonitor::CameraParameters() );
        pointingmonitor::CameraParameters& parameters = parametersVec.back();
        
        parameters.isValid =           true;
        parameters.startDate =         dateToUInt32( db_row->GetField( 1 ) );
        parameters.endDate =           dateToUInt32( db_row->GetField( 2 ) );
        parameters.version =           atol( db_row->GetField( 3 ) );
        parameters.pmtRotation =       atof( db_row->GetField( 12 ) );
        parameters.pmtScale =          atof( db_row->GetField( 13 ) );
        parameters.skyCameraRotation = atof( db_row->GetField( 14 ) );
        parameters.skyCameraScale =    atof( db_row->GetField( 15 ) );
        parameters.referencePixelX =   atof( db_row->GetField( 16 ) );
        parameters.referencePixelY =   atof( db_row->GetField( 17 ) );
    }
    
    return parametersVec;
}


pointingmonitor::CameraParameters PointingMonitor::getCameraParameters( uint32_t telescope_id, uint32_t date )
{
    vector<pointingmonitor::CameraParameters> parametersVec
        = getCameraParametersList( telescope_id, date, 1 );
    if( parametersVec.size() == 0 )
    {
        return pointingmonitor::CameraParameters();
    }
    return parametersVec[0];
}


/*****************************************************************************/
/* Functions for writing/reading calibration parameters to/from the database */
/*****************************************************************************/

vector<pointingmonitor::CalibrationParameters> PointingMonitor::getCalibrationParametersList( uint32_t telescope_id, uint32_t date, uint32_t limit )
{
    vector<pointingmonitor::CalibrationParameters> parametersVec;
    
    // create database command string
    ostringstream strbuf;
    strbuf << "select * from tblPointing_Monitor_Calibration_Parameters";
    strbuf << " where telescope_id=" << telescope_id;
    if( date != 0 )
    {
        strbuf << " and start_date<=" << date;
        strbuf << " and end_date>=" << date;
    }
    strbuf << " order by version desc";
    if( limit != 0 )
    {
        strbuf << " limit " << limit;
    }
    
    string fDBserver = getDBServer() + "VOFFLINE";
    //std::cout<<"PointingMonitor::getCalibrationParametersList "<<std::endl;
    VDB_Connection my_connection( fDBserver.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "PointingMonitor: failed to connect to database server" << endl;
        return parametersVec; // JG
    }
    if( !my_connection.make_query( strbuf.str().c_str() ) )
    {
        return parametersVec;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    int fNRows = db_res->GetRowCount();
    
    for( int i = 0; i < fNRows; i++ )
    {
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            continue;
        }
        
        // add new element to result vector and fill with numbers
        parametersVec.push_back( pointingmonitor::CalibrationParameters() );
        pointingmonitor::CalibrationParameters& parameters = parametersVec.back();
        
        parameters.isValid =   true;
        parameters.startDate = dateToUInt32( db_row->GetField( 1 ) );
        parameters.endDate =   dateToUInt32( db_row->GetField( 2 ) );
        parameters.version =   atol( db_row->GetField( 3 ) );
        parameters.x0 =        atof( db_row->GetField( 4 ) );
        parameters.x1 =        atof( db_row->GetField( 5 ) );
        parameters.x2 =        atof( db_row->GetField( 6 ) );
        parameters.x3 =        atof( db_row->GetField( 7 ) );
        parameters.y0 =        atof( db_row->GetField( 8 ) );
        parameters.y1 =        atof( db_row->GetField( 9 ) );
        parameters.y2 =        atof( db_row->GetField( 10 ) );
        parameters.y3 =        atof( db_row->GetField( 11 ) );
        parameters.z0 =        atof( db_row->GetField( 12 ) );
        parameters.z1 =        atof( db_row->GetField( 13 ) );
        parameters.z2 =        atof( db_row->GetField( 14 ) );
        parameters.z3 =        atof( db_row->GetField( 15 ) );
        parameters.l0 =        atof( db_row->GetField( 16 ) );
        parameters.l1 =        atof( db_row->GetField( 17 ) );
        parameters.l2 =        atof( db_row->GetField( 18 ) );
        parameters.l3 =        atof( db_row->GetField( 19 ) );
    }
    
    return parametersVec;
}


pointingmonitor::CalibrationParameters PointingMonitor::getCalibrationParameters( uint32_t telescope_id, uint32_t date )
{
    vector<pointingmonitor::CalibrationParameters> parametersVec
        = getCalibrationParametersList( telescope_id, date, 1 );
    if( parametersVec.size() == 0 )
    {
        return pointingmonitor::CalibrationParameters();
    }
    return parametersVec[0];
}


/****************************************************************************/
/* Functions for writing/reading uncalibrated pointing to/from the database */
/****************************************************************************/

vector<pointingmonitor::UncalibratedPointing> PointingMonitor::getUncalibratedPointing( uint32_t telescope_id, double startmjd, double stopmjd )
{
    vector<pointingmonitor::UncalibratedPointing> pointingVec;
    
    // create database command string
    ostringstream strbuf;
    strbuf.precision( 13 );
    strbuf << "select * from tblPointing_Monitor_Telescope" << telescope_id << "_Pointing";
    strbuf << " where mjd<=" << stopmjd << " and mjd>=" << startmjd;
    
    string fDBserver = getDBServer() + "VOFFLINE";
    //std::cout<<"PointingMonitor::getUncalibratedPointing "<<std::endl;
    VDB_Connection my_connection( fDBserver.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "PointingMonitor: failed to connect to database server" << endl;
        return pointingVec; // JG
    }
    if( !my_connection.make_query( strbuf.str().c_str() ) )
    {
        return  pointingVec; // JG
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    int fNRows = db_res->GetRowCount();
    
    for( int i = 0; i < fNRows; i++ )
    {
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            continue;
        }
        
        // add new element to result vector and fill with numbers
        pointingVec.push_back( pointingmonitor::UncalibratedPointing() );
        pointingmonitor::UncalibratedPointing& pointing = pointingVec.back();
        
        pointing.mjd =       atof( db_row->GetField( 0 ) );
        pointing.ra =        atof( db_row->GetField( 1 ) );
        pointing.dec =       atof( db_row->GetField( 2 ) );
        pointing.rotation =  atof( db_row->GetField( 3 ) );
        pointing.elevation = atof( db_row->GetField( 4 ) );
        pointing.ledPosY =   atof( db_row->GetField( 5 ) );
    }
    
    return pointingVec;
}


/*****************************************************/
/* Functions for calculating the calibrated pointing */
/*****************************************************/

vector<pointingmonitor::CalibratedPointing> PointingMonitor::calibratedPointing( const vector<pointingmonitor::UncalibratedPointing>& pointingVec,
        const pointingmonitor::CameraParameters& camParameters,
        const pointingmonitor::CalibrationParameters& calParameters )
{
    vector<pointingmonitor::CalibratedPointing> calPointingVec;
    calPointingVec.reserve( pointingVec.size() );
    
    for( uint32_t i = 0; i < pointingVec.size(); ++i )
    {
    
        const pointingmonitor::UncalibratedPointing& pointing = pointingVec[i];
        
        // calculate pixel in sky camera that corresponds to central PMT
        float el = pointing.elevation - float( 60 * 3.14159265 / 180 );
        float centerX = calParameters.x0 + calParameters.x1 * el + calParameters.x2 * el * el
                        + calParameters.x3 * el * el * el;
        float ledy = calParameters.l0 + calParameters.l1 * el + calParameters.l2 * el * el
                     + calParameters.l3 * el * el * el;
        float yy = calParameters.y0 + calParameters.y1 * el + calParameters.y2 * el * el
                   + calParameters.y3 * el * el * el;
        float zz = calParameters.z0 + calParameters.z1 * el + calParameters.z2 * el * el
                   + calParameters.z3 * el * el * el;
        float centerY = ( pointing.ledPosY > 0 ) ? yy + ( pointing.ledPosY - ledy ) * zz : yy;
        
        // calculate offset from reference pixel in radians
        float dx = centerX - camParameters.referencePixelX;
        float dy = centerY - camParameters.referencePixelY;
        dx *= camParameters.skyCameraScale;
        dy *= camParameters.skyCameraScale;
        
        // calculate RA and Dec of central PMT using tangential projection
        float cosr, sinr, tx;
        cosr = cos( pointing.rotation );
        sinr = sin( pointing.rotation );
        tx = dx * cosr - dy * sinr;
        dy = dx * sinr + dy * cosr;
        dx = tx;
        
        float cos0, sin0, dect, rat;
        cos0 = cos( pointing.dec );
        sin0 = sin( pointing.dec );
        dect = cos0 - dy * sin0;
        rat = pointing.ra + atan2( dx, dect );
        dect = sqrt( dx * dx + dect * dect );
        dect = atan2( dy * cos0 + sin0, dect );
        
        // rotation angle calculation is not yet implemented
        // float rot = pointing.rotation - camParameters.skyCameraRotation;
        float rot = 0.0;
        
        calPointingVec.push_back( pointingmonitor::CalibratedPointing( pointing.mjd, rat, dect, rot ) );
        
    }
    
    return calPointingVec;
}

