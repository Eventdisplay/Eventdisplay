/*! \class VPointingDB
    \brief read pointing data from the database

*/

#include "VPointingDB.h"

VPointingDB::VPointingDB( unsigned int iTelID, unsigned int irun )
{
    fStatus = false;
    fRunNumber = irun;
    fTelID = iTelID;
    
    degrad = 45. / atan( 1. );
    
    setObservatory();
    
    fMJD = 0;
    fTime = 0.;
    fTelAzimuth = 0.;
    fTelElevation = 0.;
    fTelExpectedAzimuth = 0.;
    fTelExpectedElevation = 0.;
    fEventStatus = 0;
    
    fMJDRunStart = 0;
    fTimeRunStart = 0.;
    fMJDRunStopp = 0;
    fTimeRunStopp = 0.;
    
    fDBSourceName = "don't know";
    fDBTargetDec = 0.;
    fDBTargetRA = 0.;
    fDBWobbleNorth = 0.;
    fDBWobbleEast = 0.;
    
    fCounter = 0;
    fDBNrows = 0;
    
    fmy_connection = 0;
    
    fNWarnings = 0;
    
    fTrackingCorrections = 0;
}

bool VPointingDB::initialize( string iTPointCorrection, string iVPMDirectory, bool iVPMDB, bool iUncalibratedVPM )
{

    // setup tracking corrections (expert use!)
    if( iTPointCorrection.size() > 0 )
    {
        fTrackingCorrections = new VTrackingCorrections( fTelID );
        if( !fTrackingCorrections->readTrackingCorrectionsFromDB( iTPointCorrection ) )
        {
            cout << "VPointingDB: error while reading tracking correction from VERITAS database" << endl;
            cout << "             date wrong? (example for SQL date format: \"2007-10-05\")" << endl;
            exit( -1 );
        }
    }
    else
    {
        fTrackingCorrections = 0;
    }
    
    string iTempS  = getDBServer();
    iTempS        += "/VERITAS";
    
    //std::cout<<"VPointingDB::initialize "<<std::endl;
    fmy_connection = new VDB_Connection( iTempS.c_str(), "readonly", "" ) ;
    //std::cout<<"fmy_connection->Get_Connection_Status() "<<fmy_connection->Get_Connection_Status() <<std::endl;
    
    if( !fmy_connection->Get_Connection_Status() )
    {
        cout << "VPointingDB: failed to connect to database server: " << iTempS << endl;
        fStatus = false;
        exit( EXIT_FAILURE );
    }
    
    fStatus = getDBRunInfo();
    // read pointing from VPM text file
    if( iVPMDirectory.size() > 0 )
    {
        fStatus = readPointingFromVPMTextFile( iVPMDirectory );
    }
    // read pointing from calibrated pointing monitor (VPM) entries in DB
    else if( iVPMDB )
    {
        fGoodVPM = readPointingCalibratedVPMFromDB();
        // fall back to DB pointing if reading pointing monitor data failed
        if( !fGoodVPM )
        {
            cout << "VPointingDB warning: quality-selected VPM data not available, reverting to encoder data for telescope " << getTelID() + 1 << " for full duration of run " << fRunNumber << endl;
            fStatus = readPointingFromDB();
        }
        else
        {
            fStatus = fGoodVPM;
        }
    }
    // read pointing from uncalibrated pointing monitor (VPM) entries in DB
    else if( iUncalibratedVPM )
    {
        fGoodVPM = readPointingUncalibratedVPMFromDB();
        // fall back to DB pointing if reading pointing monitor data failed
        if( !fGoodVPM )
        {
            fStatus = readPointingFromDB();
        }
        else
        {
            fStatus = fGoodVPM;
        }
    }
    // read point from DB (no pointing monitor)
    else
    {
        fStatus = readPointingFromDB();
    }
    
    // Close the connection immediately after all the reading is complete
    // this is important - open connections to the DB cause DB replication problems.
    delete_myconnection();
    
    return fStatus;
}


bool VPointingDB::updatePointing( int iMJD, double iTime )
{
    fMJD = ( unsigned int )iMJD;
    fTime = iTime;
    
    // safety net
    if( fCounter > 3 )
    {
        fCounter -= 3;
    }
    else
    {
        fCounter  = 0;
    }
    
    bool iBreak = false;
    for( unsigned int i = fCounter; i < fDBNrows; i++ )
    {
        // two cases:
        //   i) DB time > event time
        //  ii) last event and time difference is less than 10 seconds
        if( fMJD == fDBMJD[i] && ( ( fDBTime[i] >= fTime ) || ( ( i + 1 ) == fDBNrows && fTime - fDBTime[i] < 10. ) ) )
        {
            // first event, no extrapolation possible
            if( i == 0 )
            {
                fTelAzimuth = fDBTelAzimuth[i];
                fTelElevation = fDBTelElevation[i];
                fTelExpectedAzimuth = fDBTelExpectedAzimuth[i];
                fTelExpectedElevation = fDBTelExpectedElevation[i];
            }
            else if( fabs( fDBTelElevation[i] ) < 5. )
            {
                fTelAzimuth = 0.;
                fTelElevation = 0.;
                fTelExpectedAzimuth = 0.;
                fTelExpectedElevation = 0.;
                fEventStatus = 3;
            }
            // extrapolate between this and previous measurement
            else
            {
                if( fDBTime[i] > fDBTime[i - 1] )
                {
                    /////////////////////////
                    // measured elevation and azimuth
                    fTelElevation =  fDBTelElevation[i]   * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                     + fDBTelElevation[i - 1] * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                    // check that we are not around 0 or around 360 degrees
                    if( fabs( fDBTelAzimuth[i] - fDBTelAzimuth[i - 1] ) < 1. )
                    {
                        fTelAzimuth = fDBTelAzimuth[i]   * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                      + fDBTelAzimuth[i - 1] * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                    }
                    else
                    {
                        if( fDBTelAzimuth[i] > 359. && fDBTelAzimuth[i - 1] < 1. )
                        {
                            fTelAzimuth = fDBTelAzimuth[i]         * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                          + ( fDBTelAzimuth[i - 1] + 360. ) * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                        }
                        else if( fDBTelAzimuth[i] < 1. && fDBTelAzimuth[i - 1] > 359. )
                        {
                            fTelAzimuth = ( fDBTelAzimuth[i] + 360. ) * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                          + fDBTelAzimuth[i - 1]     * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                        }
                        // this should not happen -> large gap in azimuth between two events
                        else
                        {
                            fTelAzimuth = fDBTelAzimuth[i];
                            fEventStatus = 1;
                        }
                        // be sure to be in [0,360]
                        if( fTelAzimuth > 360. )
                        {
                            fTelAzimuth -= 360.;
                        }
                    }
                    /////////////////////////
                    // target elevation and azimuth
                    fTelExpectedElevation = fDBTelExpectedElevation[i]   * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                            + fDBTelExpectedElevation[i - 1] * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                    // check that we are not around 0 or around 360 degrees
                    if( fabs( fDBTelExpectedAzimuth[i] - fDBTelExpectedAzimuth[i - 1] ) < 1. )
                    {
                        fTelExpectedAzimuth = fDBTelExpectedAzimuth[i]   * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                              + fDBTelExpectedAzimuth[i - 1] * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                    }
                    else
                    {
                        if( fDBTelExpectedAzimuth[i] > 359. && fDBTelExpectedAzimuth[i - 1] < 1. )
                        {
                            fTelExpectedAzimuth = fDBTelExpectedAzimuth[i]         * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                                  + ( fDBTelExpectedAzimuth[i - 1] + 360. ) * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                        }
                        else if( fDBTelExpectedAzimuth[i] < 1. && fDBTelExpectedAzimuth[i - 1] > 359. )
                        {
                            fTelExpectedAzimuth = ( fDBTelExpectedAzimuth[i] + 360. ) * ( 1. - ( fDBTime[i] - iTime ) / ( fDBTime[i] - fDBTime[i - 1] ) )
                                                  + fDBTelExpectedAzimuth[i - 1]     * ( 1. - ( iTime - fDBTime[i - 1] ) / ( fDBTime[i] - fDBTime[i - 1] ) );
                        }
                        // this should not happen -> large gap in azimuth between two events
                        else
                        {
                            fTelExpectedAzimuth = fDBTelExpectedAzimuth[i];
                            fEventStatus = 1;
                        }
                        // be sure to be in [0,360]
                        if( fTelExpectedAzimuth > 360. )
                        {
                            fTelExpectedAzimuth -= 360.;
                        }
                    }
                    ////////////////////////////
                }
                // MJD change
                else
                {
                    fTelAzimuth = fDBTelAzimuth[i];
                    fTelElevation = fDBTelElevation[i];
                    fTelExpectedAzimuth = fDBTelExpectedAzimuth[i];
                    fTelExpectedElevation = fDBTelExpectedElevation[i];
                    fEventStatus = 2;
                }
            }
            fCounter = i;
            iBreak = true;
            break;
        }
    }
    if( !iBreak )
    {
        fNWarnings++;
        if( fNWarnings <= 30 )
        {
            cout << "VPointingDB::updatePointing WARNING: no pointing found in DB vector for telescope " << fTelID + 1;
            cout << " and time (MJD,time): " << fMJD << "\t" << setprecision( 10 ) << fTime << endl;
            if( fDBMJD.size() > 0 )
            {
                cout << "\t first DB time is at " << fDBMJD[0] << ", " << setprecision( 10 ) << fDBTime[0];
                cout << ", last DB time is at " << fDBMJD[fDBMJD.size() - 1] << ", " << setprecision( 10 ) << fDBTime[fDBTime.size() - 1] << endl;
            }
            if( fNWarnings == 30 )
            {
                cout << "-------------- more than 30 warning on DB pointing vector  -------------------------" << endl;
                cout << "                  all following warning are suppressed " << endl;
            }
        }
        
        fTelAzimuth = 0.;
        fTelElevation = 0.;
        fTelExpectedAzimuth = 0.;
        fTelExpectedElevation = 0.;
        return false;
    }
    
    return true;
}


bool VPointingDB::terminate()
{

    //std::cout<<"VPointingDB::terminate "<<std::endl;
    delete_myconnection();
    return true;
}

void VPointingDB::getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip )
{
    if( itemp.size() < 16 )
    {
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
    
    try
    {
        int y, m, d, h, min, s, ms, l;
        double gMJD;
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
    catch( const std::out_of_range& oor )
    {
        MJD = 0;
        Time = 0.;
        cout << "VPointingDB::getDBMJDTime() error: DB time string too short " << itemp.size() << ", " << oor.what() << endl;
    }
}


bool VPointingDB::getDBRunInfo()
{

    if( !fmy_connection->Get_Connection_Status() )
    {
        return false;
    }
    
    char c_query[1000];
    
    sprintf( c_query, "select * from tblRun_Info where run_id=%d", fRunNumber );
    
    if( !fmy_connection->make_query( c_query ) )
    {
        fStatus = false;
        return false;
    }
    
    TSQLResult* db_res = fmy_connection->Get_QueryResult();
    if( !db_res )
    {
        return false;
    }
    
    TSQLRow* db_row = db_res->Next();
    if( !db_row )
    {
        return false;
    }
    
    if( db_row->GetField( 6 ) && db_row->GetField( 7 ) && db_row->GetField( 17 ) &&
            db_row->GetField( 18 ) && db_row->GetField( 19 ) )
    {
        string itemp = db_row->GetField( 6 );
        getDBMJDTime( itemp, fMJDRunStart, fTimeRunStart, true );
        itemp = db_row->GetField( 7 );
        getDBMJDTime( itemp, fMJDRunStopp, fTimeRunStopp, true );
        // add 1 min to be save
        fTimeRunStopp += 60.;
        
        fDBSourceName = db_row->GetField( 19 );
        getDBSourceCoordinates( fDBSourceName, fDBTargetDec, fDBTargetRA );
        
		float dist = atof( db_row->GetField( 17 ) );
		float angl = atof( db_row->GetField( 18 ) );
        fDBWobbleNorth = dist * cos( angl * TMath::DegToRad() );
        fDBWobbleEast = dist * sin( angl * TMath::DegToRad() );
    }
    else
    {
        return false;
    }
    
    return true;
}


bool VPointingDB::readPointingFromVPMTextFile( string iDirectory )
{
    if( iDirectory.size() < 1 )
    {
        return false;
    }
    
    char fname[1200];
    sprintf( fname, "%s/pointing_VPM.%d.t%d.dat", iDirectory.c_str(), fRunNumber, getTelID() + 1 );
    
    ifstream is;
    is.open( fname );
    if( !is )
    {
        cout << "VPointingDB: error opening pointing monitor data file: " << fname << endl;
        cout << "exiting...." << endl;
        exit( 0 );
    }
    string is_line;
    int nLines = 0;
    double iMJD = 0.;
    double iTemp = 0.;
    double iRA = 0.;
    double iDec = 0.;
    double iITime = 0.;
    double az = 0.;
    double ze = 0.;
    
    while( getline( is, is_line ) )
    {
        // This is the expected format (vpm v.0.9):
        // MJD REQ_AZ REQ_EL REQ_RA2000 REQ_DEC2000 POSITIONER_RA2000 POSITIONER_DEC2000 VPM_AZ VPM_EL VPM_RA2000 VPM_DEC2000
        // test line for completeness (expect 11 columns)
        int nC = 0;
        istringstream is_streamT( is_line );
        while( !(is_streamT>>std::ws).eof() )
        {
            is_streamT >> iTemp;
            nC++;
        };
        if( nC != 11 )
        {
            cout << "\t warning: incomplete line in pointing monitor file" << endl;
            continue;
        }
        // read a new line
        istringstream is_stream( is_line );
        is_stream >> iMJD;
        iITime = modf( iMJD, &iMJD );
        // ignore the following 8 columns
        for( int i = 0; i < 8; i++ )
        {
            is_stream >> iTemp;
        }
        is_stream >> iRA;
        is_stream >> iDec;
        if( iRA > 0. && iDec > 0. )
        {
            fDBMJD.push_back( ( unsigned int )( iMJD ) );
            fDBTime.push_back( iITime * 86400. );
            fDBTelElevationRaw.push_back( 0. );
            fDBTelAzimuthRaw.push_back( 0. );
            fDBTelRA.push_back( iRA * degrad );
            fDBTelDec.push_back( iDec * degrad );
            getHorizonCoordinates( fDBMJD.back(), fDBTime.back(), iDec * degrad, iRA * degrad, az, ze );
            fDBTelElevation.push_back( 90. - ze );
            fDBTelAzimuth.push_back( az );
            fDBTelExpectedElevation.push_back( 0. );
            fDBTelExpectedAzimuth.push_back( 0. );
        }
        nLines++;
    };
    cout << "Reading pointing data from VERITAS pointing monitor file for telescope " << getTelID() + 1;
    cout << ": found " << nLines << " rows in file " << fname << endl;
    if( nLines == 0 )
    {
        return false;
    }
    
    fDBNrows = fDBMJD.size();
    
    return true;
}

/*
    read calibrated (default) pointing monitor data from DB
*/
bool VPointingDB::readPointingCalibratedVPMFromDB()
{
    double startMJD = fMJDRunStart + ( fTimeRunStart / 86400. );
    double stopMJD = fMJDRunStopp + ( fTimeRunStopp / 86400. );
    
    // mysql query //
    char c_query[1000];
    sprintf( c_query, "SELECT mjd,ra,decl FROM tblPointing_Monitor_Telescope%d_Calibrated_Pointing WHERE mjd<=%.13f AND mjd>=%.13f", getTelID(), stopMJD, startMJD );
    
    // use VOFFLINE database for calibrated VPM //
    string iTempS  = getDBServer();
    iTempS += "/VOFFLINE";
    
    //std::cout<<"VPointingDB::readPointingCalibratedVPMFromDB "<<std::endl;
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "VPointingDB: failed to connect to database server: " << iTempS << endl;
        return false;
    }
    
    if( !my_connection.make_query( c_query ) )
    {
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    int fNRows = db_res->GetRowCount();
    cout << "Reading calibrated pointing monitor (VPM) data from database for telescope " << getTelID() + 1;
    cout << ": found " << fNRows << " rows in database" << endl;
    
    // get VPM quality flag from database
    char cflag_query[1000];
    sprintf( cflag_query, "SELECT vpm_config_mask FROM tblRun_Analysis_Comments WHERE run_id = %d", fRunNumber );
    if( !my_connection.make_query( cflag_query ) )
    {
        cout << "VPointingDB: error, missing VPM quality flag" << endl;
        return false;
    }
    TSQLResult* db_flag = my_connection.Get_QueryResult();
    
    
    TSQLRow* flag_row = db_flag->Next();
    if( !flag_row || !flag_row->GetField( 0 ) )
    {
        cout << "VPointingDB: error while reading VPM quality flag" << endl;
        return false;
    }
    int maskVPM = atoi( flag_row->GetField( 0 ) );
    
    if( getTelID() == 0 )      // T1 bad
    {
        if( maskVPM == 0 || maskVPM == 2 || maskVPM == 4 || maskVPM == 8 )
        {
            return false;
        }
        if( maskVPM == 6 || maskVPM == 10 || maskVPM == 12 || maskVPM == 14 )
        {
            return false;
        }
    }
    if( getTelID() == 1 )      // T2 bad
    {
        if( maskVPM == 0 || maskVPM == 1 || maskVPM == 4 || maskVPM == 8 )
        {
            return false;
        }
        if( maskVPM == 5 || maskVPM == 9 || maskVPM == 12 || maskVPM == 13 )
        {
            return false;
        }
    }
    if( getTelID() == 2 )      // T3 bad
    {
        if( maskVPM == 0 || maskVPM == 1 || maskVPM == 2 || maskVPM == 8 )
        {
            return false;
        }
        if( maskVPM == 3 || maskVPM == 9 || maskVPM == 10 || maskVPM == 11 )
        {
            return false;
        }
    }
    if( getTelID() == 3 )      // T4 bad
    {
        if( maskVPM == 0 || maskVPM == 1 || maskVPM == 2 || maskVPM == 4 )
        {
            return false;
        }
        if( maskVPM == 3 || maskVPM == 5 || maskVPM == 6 || maskVPM == 7 )
        {
            return false;
        }
    }
    
    // loop over all db entries
    string itemp;
    double iMJD;
    double iITime = 0.;
    double az = 0.;
    double ze = 0.;
    double iRA = 0.;
    double iDec = 0.;
    
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            continue;
        }
        
        iMJD = atof( db_row->GetField( 0 ) );
        iITime = modf( iMJD, &iMJD );
        fDBMJD.push_back( ( unsigned int )iMJD );
        fDBTime.push_back( iITime * 86400. );
        fDBTelElevationRaw.push_back( 0. );
        fDBTelAzimuthRaw.push_back( 0. );
        iRA = atof( db_row->GetField( 1 ) );
        iDec = atof( db_row->GetField( 2 ) );
        fDBTelRA.push_back( iRA * degrad );
        fDBTelDec.push_back( iDec * degrad );
        getHorizonCoordinates( fDBMJD.back(), fDBTime.back(), iDec * degrad, iRA * degrad, az, ze );
        fDBTelElevation.push_back( 90. - ze );
        fDBTelAzimuth.push_back( az );
        fDBTelExpectedElevation.push_back( 0. );
        fDBTelExpectedAzimuth.push_back( 0. );
    }
    
    fDBNrows = fDBMJD.size();
    
    return true;
}

/*
    read uncalibrated pointing monitor data from DB
*/
bool VPointingDB::readPointingUncalibratedVPMFromDB()
{
    double iMJD = 0.;
    double iITime = 0.;
    double iRA = 0.;
    double iDec = 0.;
    double az = 0.;
    double ze = 0.;
    
    // get date in year/month/day
    int year, month, day, j_status;
    double fracday;
    VAstronometry::vlaDjcl( ( double )fMJDRunStart, &year, &month, &day, &fracday, &j_status );
    year = year * 10000;
    month = month * 100;
    
    uint32_t telescopeID = getTelID();
    uint32_t obsDate = year + month + day;
    double startMJD = fMJDRunStart + ( fTimeRunStart / 86400. );
    double stopMJD = fMJDRunStopp + ( fTimeRunStopp / 86400. );
    
    vector< pointingmonitor::UncalibratedPointing > uncalibratedPointing;
    vector< pointingmonitor::CalibratedPointing > VPMcalibratedPointing;
    pointingmonitor::CameraParameters cameraParameters;
    pointingmonitor::CalibrationParameters calibrationParameters;
    
    pointingmonitor::PointingMonitor vpm;
    // read vector of uncalibrated pointings
    uncalibratedPointing = vpm.getUncalibratedPointing( telescopeID, startMJD, stopMJD );
    // read camera parameters
    cameraParameters = vpm.getCameraParameters( telescopeID, obsDate );
    // read calibration parameters
    calibrationParameters = vpm.getCalibrationParameters( telescopeID, obsDate );
    // do some checks and show some diagnostics output
    if( uncalibratedPointing.size() > 0 )
    {
        cout << "Found " << uncalibratedPointing.size() << " pointing monitor entries for telescope " << telescopeID + 1 << endl;
    }
    else
    {
        cout << "VPointingDB::readPointingMonitorFromDB warning: Unable to find any pointing monitor data for this telescope during the run," << endl;
        cout << "   reverting to encoder data for this telescope for full duration of the run";
        cout << " (telescope " << telescopeID + 1 << ")" << endl;
        
        return false;
    }
    
    if( cameraParameters.isValid && calibrationParameters.isValid )
    {
        VPMcalibratedPointing = vpm.calibratedPointing( uncalibratedPointing, cameraParameters, calibrationParameters );
    }
    
    // hard-wired limit (in arcsec) for offset between VPM and target position //
    // note: previously was 180, now loosened to 300
    double vpmlimit = 300.;
    
    double decoff = 0;
    double raoff = 0;
    
    for( uint32_t i = 0; i < VPMcalibratedPointing.size(); i++ )
    {
        decoff = fabs( ( 3600. * VPMcalibratedPointing[0].dec * degrad ) - ( 3600. * ( fDBTargetDec + fDBWobbleNorth ) ) );
        raoff = fabs( ( 3600. * VPMcalibratedPointing[0].ra * degrad ) - ( 3600. * ( fDBTargetRA + fDBWobbleEast ) ) );
        if( decoff > vpmlimit || raoff > vpmlimit )
        {
            cout << "VPointingDB::readPointingMonitorFromDB warning: For part of this run the pointing monitor data is off by more than ";
            cout << vpmlimit;
            cout << " arcsec from the target position for this telescope, reverting to encoder data for this telescope for full duration of the run" << endl;
            
            return false;
        }
        
    }
    
    // hard-wired time bin (in sec) to check for VPM data //
    int tbinwidth = 10;
    // hard-wired minimum fraction of run missing VPM data //
    double vpmrunfrac = 0.10;
    double timebin = ( fTimeRunStopp - fTimeRunStart ) / tbinwidth;
    
    double timelimit = timebin * vpmrunfrac;
    timebin = ( int )timebin;
    timelimit = ( int )timelimit;
    
    double starttime = startMJD * 86400.;
    int nbad = 0;
    
    for( int k = 0; k < timebin; k++ )
    {
        bool timegood = false;
        
        for( uint32_t i = 0; i < VPMcalibratedPointing.size(); i++ )
        {
        
            if( starttime < ( VPMcalibratedPointing[i].mjd * 86400 ) && ( starttime + tbinwidth ) > ( VPMcalibratedPointing[i].mjd * 86400 ) )
            {
                timegood = true;
                break;
            }
        }
        if( !timegood )
        {
            nbad += 1;
        }
        
        starttime += tbinwidth;
    }
    
    if( nbad > timelimit )
    {
        cout << "VPointingDB::readPointingMonitorFromDB warning: pointing monitor data is missing for more than ";
        cout << vpmrunfrac * 100.;
        cout << "% of the run for this telescope, reverting to encoder data for this telescope for full duration of the run" << endl;
        
        return false;
    }
    
    for( uint32_t i = 0; i < VPMcalibratedPointing.size(); i++ )
    {
        iMJD = VPMcalibratedPointing[i].mjd;
        iITime = modf( iMJD, &iMJD );
        fDBMJD.push_back( ( unsigned int )( iMJD ) );
        fDBTime.push_back( iITime * 86400. );
        fDBTelElevationRaw.push_back( 0. );
        fDBTelAzimuthRaw.push_back( 0. );
        
        iDec = VPMcalibratedPointing[i].dec;
        iRA = VPMcalibratedPointing[i].ra;
        getHorizonCoordinates( fDBMJD.back(), fDBTime.back(), iDec * degrad, iRA * degrad, az, ze );
        fDBTelElevation.push_back( 90. - ze );
        fDBTelAzimuth.push_back( az );
        fDBTelExpectedElevation.push_back( 0. );
        fDBTelExpectedAzimuth.push_back( 0. );
    }
    
    fDBNrows = fDBMJD.size();
    
    return true;
}

bool VPointingDB::readPointingFromDB()
{

    if( !fmy_connection->Get_Connection_Status() )
    {
        return false;
    }
    
    char iDate1[200];
    char iDate2[200];
    
    // get date in year/month/day
    int year, month, day, j_status;
    double fracday;
    VAstronometry::vlaDjcl( ( double )fMJDRunStart, &year, &month, &day, &fracday, &j_status );
    
    int hour = ( int )( fTimeRunStart / 3600. );
    int minute = ( int )( ( fTimeRunStart - hour * 3600. ) / 60. );
    int sec = ( int )( fTimeRunStart - hour * 3600. - minute * 60. );
    sprintf( iDate1, "%d%02d%02d%02d%02d%02d000", ( int )year, ( int )month, ( int )day, hour, minute, sec );
    
    VAstronometry::vlaDjcl( ( double )fMJDRunStopp, &year, &month, &day, &fracday, &j_status );
    hour = ( int )( fTimeRunStopp / 3600. );
    minute = ( int )( ( fTimeRunStopp - hour * 3600. ) / 60. );
    sec = ( int )( fTimeRunStopp - hour * 3600. - minute * 60. );
    sprintf( iDate2, "%d%02d%02d%02d%02d%02d000", ( int )year, ( int )month, ( int )day, hour, minute, sec );
    
    // mysql queries
    char c_query[1000];
    sprintf( c_query, "SELECT timestamp, elevation_raw, azimuth_raw, elevation_meas, azimuth_meas, elevation_target, azimuth_target FROM tblPositioner_Telescope%d_Status WHERE timestamp >= %s AND timestamp <= %s", getTelID(), iDate1, iDate2 );
    
    if( !fmy_connection->make_query( c_query ) )
    {
        return false;
    }
    TSQLResult* db_res = fmy_connection->Get_QueryResult();
    
    
    int fNRows = db_res->GetRowCount();
    
    // loop over all db entries
    int iMJD = 0;
    double iTime = 0.;
    
    cout << "Reading pointing data from VERITAS database for telescope " << getTelID() + 1;
    cout << ": found " << fNRows << " rows in database" << endl;
    
    string itemp;
    double el = 0.;
    double az = 0.;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            continue;
        }
        itemp = db_row->GetField( 0 );
        getDBMJDTime( itemp, iMJD, iTime, false );
        fDBMJD.push_back( ( unsigned int )iMJD );
        fDBTime.push_back( iTime );
        // reapply tracking corrections
        if( fTrackingCorrections )
        {
            // get elevation_raw, azimuth_raw from DB and reapply tracking corrections
            fTrackingCorrections->applyTrackingCorrections( atof( db_row->GetField( 1 ) ), atof( db_row->GetField( 2 ) ), el, az );
        }
        else
        {
            // get elevation_meas, azimuth_meas from DB
            el = atof( db_row->GetField( 3 ) );
            az = atof( db_row->GetField( 4 ) );
        }
        fDBTelElevationRaw.push_back( atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi() );
        fDBTelAzimuthRaw.push_back( atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi() );
        fDBTelElevation.push_back( el * 180. / TMath::Pi() );
        fDBTelAzimuth.push_back( az * 180. / TMath::Pi() );
        fDBTelRA.push_back( 0. );
        fDBTelDec.push_back( 0. );
        fDBTelExpectedElevation.push_back( atof( db_row->GetField( 5 ) ) * 180. / TMath::Pi() );
        fDBTelExpectedAzimuth.push_back( atof( db_row->GetField( 6 ) ) * 180. / TMath::Pi() );
    }
    fDBNrows = fDBMJD.size();
    
    return true;
}


void VPointingDB::getDBSourceCoordinates( string iSource, float& iEVNTargetDec, float& iEVNTargetRA )
{
    if( !fmy_connection->Get_Connection_Status() )
    {
        return;
    }
    
    char c_query[1000];
    
    sprintf( c_query, "select * from tblObserving_Sources where source_id like convert( _utf8 \'%s\' using latin1)", iSource.c_str() );
    
    if( !fmy_connection->make_query( c_query ) )
    {
        return;
    }
    TSQLResult* db_res = fmy_connection->Get_QueryResult();
    
    
    
    TSQLRow* db_row = db_res->Next();
    
    iEVNTargetDec = atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi();
    iEVNTargetRA = atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi();
    return;
}

/*
 * tree with pointing values as given by DB
 *
 */
TTree* VPointingDB::getTreePointingDB()
{
    char hname[200];
    char htitle[200];
    
    unsigned int MJD = 0;
    double Time = 0.;
    float iTelElR = 0.;
    float iTelAzR = 0.;
    float iTelEl = 0.;
    float iTelAz = 0.;
    float iTelElT = 0.;
    float iTelAzT = 0.;
    float iTelRA = 0.;
    float iTelDec = 0.;
    
    sprintf( hname, "db_pointing_%d", getTelID() + 1 );
    sprintf( htitle, "pointing from DB (Telescope %d)", getTelID() + 1 );
    TTree* tD = new TTree( hname, htitle );
    tD->Branch( "MJD", &MJD, "MJD/i" );
    tD->Branch( "Time", &Time, "Time/D" );
    // elevation / azimuth before application of tracking corrections
    tD->Branch( "ElRaw", &iTelElR, "ElRaw/F" );
    tD->Branch( "AzRaw", &iTelAzR, "AzRaw/F" );
    // elevation / azimuth after application of tracking corrections
    tD->Branch( "El", &iTelEl, "El/F" );
    tD->Branch( "Az", &iTelAz, "Az/F" );
    // expected elevation / azimuth (value requested by pointing software)
    tD->Branch( "ElTelExpected", &iTelElT, "ElTelExpected/F" );
    tD->Branch( "AzTelExpected", &iTelAzT, "AzTelExpected/F" );
    // ra/dec as read from DB
    tD->Branch("RA", &iTelRA, "RA/F" );
    tD->Branch("Dec", &iTelDec, "Dec/F" );
    
    for( unsigned int i = 0; i < fDBMJD.size(); i++ )
    {
        MJD = fDBMJD[i];
        Time = fDBTime[i];
        iTelElR = fDBTelElevationRaw[i];
        iTelAzR = fDBTelAzimuthRaw[i];
        iTelEl = fDBTelElevation[i];
        iTelAz = fDBTelAzimuth[i];
        iTelElT = fDBTelExpectedElevation[i];
        iTelAzT = fDBTelExpectedAzimuth[i];
        iTelRA = fDBTelRA[i];
        iTelDec = fDBTelDec[i];
        
        tD->Fill();
    }
    
    return tD;
}

/*
 *
 *
 */
void VPointingDB::getHorizonCoordinates( int MJD, double time, double decJ2000, double raJ2000, double& az, double& ze )
{
    raJ2000 /= degrad;
    decJ2000 /= degrad;
    
    // first precess target to current epoch
    VAstronometry::vlaPreces( 2451545.0 - 2400000.5, (double)(MJD+0.5), &raJ2000, &decJ2000 );
    
    // convert ra into hour angle
    double ha = 0.;
    double iTime = 0.;
    double iSid = 0.;
    // convert time to fraction of a day
    iTime = time / 86400.;
    // get Greenwich sideral time
    iSid = VAstronometry::vlaGmsta( ( double )MJD, iTime );
    // calculate local sideral time
    iSid = iSid - fObsLongitude;
    // calculate right ascension
    ha = VAstronometry::vlaDranrm( iSid - raJ2000 );
    // calculate equatorial coordinates
    double el = 0.;
    VAstronometry::vlaDe2h( ha, decJ2000, fObsLatitude, &az, &el );
    el *= degrad;
    
    ze = 90. - el;
    az *= degrad;
}

void VPointingDB::setObservatory( double iLongitude, double iLatitude )
{
    fObsLongitude = iLongitude / degrad;
    fObsLatitude = iLatitude / degrad;
}
