/*! \file VRunStats
    \brief get run and DQM data from DB

    Missing:

     * add remaining fields from tblRun_Info (run_type, observing_mode, etc.)
     * add mean HV values


*/

#include "VRunStats.h"


VDBRunData::VDBRunData()
{

    runNumber = 0;
    ConfigMask = 0;
    SourceID = "";
    StartMJD = 0.;
    StopMJD = 0.;
    start_time = "";
    stop_time = "";
    Duration = 0.;
    RA = 0.;
    Dec = 0.;
    offsetRA = 0.;
    offsetDec = 0.;
    offset_distance = 0.;
    offset_angle = 0.;
    GalLong1958 = 0.;
    GalLat1958 = 0.;
}

VDBAnalysisComments::VDBAnalysisComments()
{
    runNumber = 0;
    status = 0;
    usable_duration = 0.;
    comment = "";
}


VDBSourceData::VDBSourceData()
{
    SourceID = "";
    fDec = 0.;
    fRA = 0.;
}

VDBDefaultHV::VDBDefaultHV()
{
    startMJD = 0.;
    stopMJD = 0.;
    telescope = 0;
    channel = 0;
    default_HV = 0.;
}

VDBHV::VDBHV()
{
    startMJD = 0.;
    stopMJD = 0.;
    telescope = 0;
    channel = 0;
    HV = 0.;
}

VDBCameraStatus::VDBCameraStatus()
{
    MJD = 0.;
    telescope = 0;
    temp1 = 0.;
    temp2 = 0.;
    hum1 = 0.;
    hum2 = 0.;
    light1 = 0.;
    light2 = 0.;
}

VDBFIRData::VDBFIRData()
{
    MJD = 0.;

    tel_ID = 0;
    ambient_temp = 0.;
    radiant_sky_temp = 0.;
    radiant_sky_temp_cor = 0.;
}

VDBWeatherData::VDBWeatherData()
{
    MJD = 0.;

    AirTemperature = 0.;
    AirPressure = 0.;
    AirHumidity = 0.;
    AirWindSpeed = 0.;
    AirWindGust = 0.;
    AirWindDirection = 0.;
}


VDQMRunInfo::VDQMRunInfo()
{

    runNumber = 0;
    fieldVersion = 0;

    StartMJD = 0.;
    SourceID = "";
    ObsMode = "";
    WobbleOff_estim_deg = 0.;
    WobbleDir_estim_deg = 0.;
    El_mean_deg = 0.;
    Az_mean_deg = 0.;
    L3_runTime_s = 0.;
    L3_lifeTime_s = 0.;
    Rate_mean = 0.;
    Rate_RMS = 0.;
    Rate_Chi2 = 0.;

    L3_threshold = 0.;

    for( unsigned int i = 0; i < 4; i++ )
    {
        L1_rate[i] = 0.;
        Median_current_uA[i] = 0.;
        L2_threshold[i] = 0.;
        L2_threshold_mV[i] = 0.;
        L2_threshold_Spread_mV[i] = 0.;
        L2_throughput[i] = 0.;
    }

    GPS_runtime_min = 0.;
    GPS_lifetime_min = 0.;

    MoonEl_mean_deg = 0.;
    MoonPhase_mean = 0.;

    NTel = 0;
    TelBitMask = 0;

    FIR0_mean = 0.;
    FIR0_RMS = 0.;

    for( unsigned int i = 0; i < 4; i++ )
    {
        MoonSeparation_deg[i] = 0.;
        RA_mean[i] = 0.;
        RA_RMS[i] = 0.;
        Dec_mean[i] = 0.;
        Dec_RMS[i] = 0.;

        L2_rate[i] = 0.;
        MissingEvents[i] = 0.;
        Muons_N[i] = 0.;
        Muons_MeanCharge[i] = 0.;
        Muons_RMSCharge[i] = 0.;
        PMT_GainLaser[i] = 0.;

        SuppressedChannels_Num[i] = 0.;

        PedVar_mean[i] = 0.;

        LaserRun[i] = 0;
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VRunStats::VRunStats()
{
    fDebug = false;

    fReadfullHVEvndispData = false;

    fMax_run_id = 0;
    fMin_run_id = 0;

    setTimeRange( "2009-01-01", "2009-01-02" );
}


void VRunStats::setTimeRange( string iStart, string iStop )
{
    fStartDate = iStart;
    fStopDate  = iStop;
}

bool VRunStats::readDBAnalysisComments()
{
    string itemp = getDBServer() + "/VOFFLINE";

    //std::cout<<"VRunStats::readDBAnalysisComments "<<std::endl;
    VDB_Connection my_connection( itemp.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        return false;
    }

    char c_query[1000];
    sprintf( c_query, "select run_id , data_category   , status   , status_reason , tel_cut_mask , usable_duration , time_cut_mask , light_level , vpm_config_mask , authors  , comment from tblRun_Analysis_Comments where run_id >= %d and run_id <= %d", fMin_run_id, fMax_run_id );
    cout << c_query << endl;

    if( !my_connection.make_query( c_query ) )
    {
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();


    for( unsigned int i = 0; i < 5; i++ )
    {
        cout << i << "\t" << db_res->GetFieldName( i ) << endl;
    }

    int fNRows = db_res->GetRowCount();

    cout << "reading analysis comments from DB" << endl;

    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();

        if( !db_row )
        {
            break;
        }
        if( !db_row->GetField( 0 ) )
        {
            continue;
        }
        if( !db_row->GetField( 1 ) )
        {
            continue;
        }
        if( !db_row->GetField( 2 ) )
        {
            continue;
        }
        if( !db_row->GetField( 4 ) )
        {
            continue;
        }

        fAnalysisComments.push_back( new VDBAnalysisComments() );

        fAnalysisComments.back()->runNumber = atoi( db_row->GetField( 0 ) );
        // seg fault ??
        itemp = db_row->GetField( 1 );
        if( itemp == "good_run" )
        {
            fAnalysisComments.back()->status = 0;
        }
        else if( itemp == "needs_adjustments" )
        {
            fAnalysisComments.back()->status = 1;
        }
        else if( itemp == "minor_problems" )
        {
            fAnalysisComments.back()->status = 2;
        }
        else if( itemp == "major_problems" )
        {
            fAnalysisComments.back()->status = 3;
        }
        else if( itemp == "do_not_use" )
        {
            fAnalysisComments.back()->status = 4;
        }
        else if( itemp == "unknown" )
        {
            fAnalysisComments.back()->status = 5;
        }

        itemp = db_row->GetField( 2 );
        fAnalysisComments.back()->usable_duration  = atof( itemp.substr( 0, itemp.find( ":" ) ).c_str() ) * 3600.;
        fAnalysisComments.back()->usable_duration += atof( itemp.substr( itemp.find( ":" ) + 1, itemp.find( ":" ) ).c_str() ) * 60.;
        fAnalysisComments.back()->usable_duration += atof( itemp.substr( itemp.rfind( ":" ) + 1, itemp.size() ).c_str() );

        fAnalysisComments.back()->comment = db_row->GetField( 4 );
    }

    return true;
}


bool VRunStats::readFromDB()
{
    string itemp = getDBServer() + "/VERITAS";

    //std::cout<<"VRunStats::readFromDB "<<std::endl;
    VDB_Connection my_connection( itemp.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "error connecting to db" << endl;
        return false;
    }

    TSQLServer* f_db = my_connection.Get_ConnectionResult();

    if( !readDBRun_IDs( f_db ) )
    {
        return false;
    }

    if( !readDBSourceInfo( f_db ) )
    {
        return false;
    }

    if( !readDBWeatherInfo( f_db ) )
    {
        return false;
    }

    if( !readDB_default_HVEvndispData( f_db ) )
    {
        return false;
    }

    if( fReadfullHVEvndispData )
    {
        if( !readDBHVEvndispData( f_db ) )
        {
            return false;
        }
    }

    if( !readDBCameraStatus( f_db ) )
    {
        return false;
    }

    if( !readDBFIRInfo( f_db ) )
    {
        return false;
    }

    if( !readDBAnalysisComments() )
    {
        return false;
    }

    if( !readDBRunInfo( f_db ) )
    {
        return false;
    }



    return true;
}


bool VRunStats::readDBHVEvndispData( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    for( unsigned int i = 0; i < 4; i++ )
    {
        sprintf( c_query, "select * from tblHV_Telescope%d_Voltages where db_start_time >= \"%s 14:00:00\"  and db_end_time < \"%s 16:00:00\"", i,  fStartDate.c_str(), fStopDate.c_str() );
        cout << c_query << endl;

        TSQLResult* db_res = f_db->Query( c_query );
        if( !db_res )
        {
            return false;
        }

        int fNRows = db_res->GetRowCount();

        cout << "reading HVs for telescope " << i + 1 << " from DB" << endl;
        cout << "\t number of rows for HV read form DB: " << fNRows << endl;

        string itemp;
        for( int j = 0; j < fNRows; j++ )
        {
            TSQLRow* db_row = db_res->Next();

            if( !db_row )
            {
                break;
            }


            int iMJD = 0;
            double iTime = 0.;
            unsigned int iChannel = 0;

            itemp = db_row->GetField( 2 );
            iChannel = atoi( itemp.c_str() );
            if( iChannel < 500 )
            {

                fHVEvndispData.push_back( new VDBHV() );
                fHVEvndispData.back()->telescope = i;
                fHVEvndispData.back()->channel = iChannel;


                if( db_row->GetField( 0 ) )
                {
                    itemp = db_row->GetField( 0 );
                    getDBMJDTime( itemp, iMJD, iTime, true );
                    fHVEvndispData.back()->startMJD = ( double )iMJD + iTime / 86400;
                }
                else
                {
                    fHVEvndispData.back()->startMJD = 0.;
                }
                if( db_row->GetField( 1 ) )
                {
                    itemp = db_row->GetField( 1 );
                    getDBMJDTime( itemp, iMJD, iTime, true );
                    fHVEvndispData.back()->stopMJD = ( double )iMJD + iTime / 86400;
                }
                else
                {
                    fHVEvndispData.back()->stopMJD = 0.;
                }
                if( db_row->GetField( 3 ) )
                {
                    itemp = db_row->GetField( 3 );
                    fHVEvndispData.back()->HV = atof( itemp.c_str() );
                }
                else
                {
                    fHVEvndispData.back()->HV = 0.;
                }
            }
        }
    }
    return true;
}


bool VRunStats::readDB_default_HVEvndispData( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    for( unsigned int i = 0; i < 4; i++ )
    {
        sprintf( c_query, "select * from tblHV_Telescope%d_DefaultVoltage where db_start_time >= \"2005-01-01 14:00:00\" ", i );
        cout << c_query << endl;

        TSQLResult* db_res = f_db->Query( c_query );
        if( !db_res )
        {
            return false;
        }

        int fNRows = db_res->GetRowCount();
        //      for( unsigned int f = 0; f < db_res->GetFieldCount(); f++ ) cout << f << "\t" << db_res->GetFieldName( f ) << endl;

        cout << "reading HV default settings for telescope " << i + 1 << " from DB" << endl;
        cout << "\t number of rows for defaultHV read form DB: " << fNRows << endl;


        string itemp;
        for( int j = 0; j < fNRows; j++ )
        {
            TSQLRow* db_row = db_res->Next();

            if( !db_row )
            {
                break;
            }


            int iMJD = 0;
            double iTime = 0.;
            unsigned int iChannel = 0;

            itemp = db_row->GetField( 2 );
            iChannel = atoi( itemp.c_str() );
            if( iChannel < 500 )
            {

                fdefaultHVEvndispData.push_back( new VDBDefaultHV() );
                fdefaultHVEvndispData.back()->telescope = i;
                fdefaultHVEvndispData.back()->channel = iChannel;


                if( db_row->GetField( 0 ) )
                {
                    itemp = db_row->GetField( 0 );
                    getDBMJDTime( itemp, iMJD, iTime, true );
                    fdefaultHVEvndispData.back()->startMJD = ( double )iMJD + iTime / 86400;
                }
                else
                {
                    fdefaultHVEvndispData.back()->startMJD = 0.;
                }
                if( db_row->GetField( 1 ) )
                {
                    itemp = db_row->GetField( 1 );
                    getDBMJDTime( itemp, iMJD, iTime, true );
                    fdefaultHVEvndispData.back()->stopMJD = ( double )iMJD + iTime / 86400;
                }
                else
                {
                    fdefaultHVEvndispData.back()->stopMJD = 0.;
                }
                if( db_row->GetField( 3 ) )
                {
                    itemp = db_row->GetField( 3 );
                    fdefaultHVEvndispData.back()->default_HV = atof( itemp.c_str() );
                }
                else
                {
                    fdefaultHVEvndispData.back()->default_HV = 0.;
                }
            }
        }
    }
    return true;
}

bool VRunStats::readDBCameraStatus( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    // read run data
    sprintf( c_query, "select * from tblCamera_Status where timestamp >= \"%s 14:00:00\"  and timestamp < \"%s 16:00:00\"", fStartDate.c_str(), fStopDate.c_str() );
    cout << c_query << endl;

    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    for( unsigned int i = 0; i < 8; i++ )
    {
        cout << i << "\t" << db_res->GetFieldName( i ) << endl;
    }

    int fNRows = db_res->GetRowCount();

    cout << "reading camera status from DB" << endl;
    cout << "\t number of rows for camera status read from DB: " << fNRows << endl;

    string itemp;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();

        if( !db_row )
        {
            break;
        }

        fCameraStatus.push_back( new VDBCameraStatus() );

        int iMJD = 0;
        double iTime = 0.;

        itemp = db_row->GetField( 0 );
        getDBMJDTime( itemp, iMJD, iTime, true );
        fCameraStatus.back()->MJD = ( double )iMJD + iTime / 86400;

        fCameraStatus.back()->telescope = atoi( db_row->GetField( 1 ) );
        fCameraStatus.back()->temp1 = atof( db_row->GetField( 2 ) );
        fCameraStatus.back()->temp2 = atof( db_row->GetField( 3 ) );
        fCameraStatus.back()->hum1 = atof( db_row->GetField( 4 ) );
        fCameraStatus.back()->hum2 = atof( db_row->GetField( 5 ) );
        fCameraStatus.back()->light1 = atof( db_row->GetField( 6 ) );
        fCameraStatus.back()->light2 = atof( db_row->GetField( 7 ) );
    }

    return true;
}



bool VRunStats::readDBFIRInfo( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    // read run data
    sprintf( c_query, "select * from tblFIR_Pyrometer_Info where timestamp >= \"%s 14:00:00\"  and timestamp < \"%s 16:00:00\"", fStartDate.c_str(), fStopDate.c_str() );
    cout << c_query << endl;

    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    for( unsigned int i = 0; i < 4; i++ )
    {
        cout << i << "\t" << db_res->GetFieldName( i ) << endl;
    }

    int fNRows = db_res->GetRowCount();

    cout << "reading FIR data from DB" << endl;
    cout << "\t number of rows for FIR data read from DB: " << fNRows << endl;

    string itemp;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();

        if( !db_row )
        {
            break;
        }

        fFIRData.push_back( new VDBFIRData() );

        int iMJD = 0;
        double iTime = 0.;

        itemp = db_row->GetField( 0 );
        getDBMJDTime( itemp, iMJD, iTime, true );
        fFIRData.back()->MJD = ( double )iMJD + iTime / 86400;

        fFIRData.back()->tel_ID = atoi( db_row->GetField( 1 ) );
        fFIRData.back()->ambient_temp = atof( db_row->GetField( 2 ) );
        fFIRData.back()->radiant_sky_temp = atof( db_row->GetField( 3 ) );
        fFIRData.back()->radiant_sky_temp_cor = fFIRData.back()->radiant_sky_temp;
    }
    return true;
}


bool VRunStats::readDBWeatherInfo( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    // read run data
    sprintf( c_query, "select * from tblWeather_Status where timestamp >= \"%s 14:00:00\"  and timestamp < \"%s 16:00:00\"", fStartDate.c_str(), fStopDate.c_str() );
    cout << c_query << endl;

    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    for( unsigned int i = 0; i < 9; i++ )
    {
        cout << i << "\t" << db_res->GetFieldName( i ) << endl;
    }

    int fNRows = db_res->GetRowCount();

    cout << "reading weather info from DB" << endl;
    cout << "\t number of rows for weather info read form DB: " << fNRows << endl;

    string itemp;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();

        if( !db_row )
        {
            break;
        }

        fWeatherData.push_back( new VDBWeatherData() );

        int iMJD = 0;
        double iTime = 0.;

        itemp = db_row->GetField( 0 );
        getDBMJDTime( itemp, iMJD, iTime, true );
        fWeatherData.back()->MJD = ( double )iMJD + iTime / 86400;

        fWeatherData.back()->AirTemperature = atof( db_row->GetField( 6 ) );
        fWeatherData.back()->AirHumidity = atof( db_row->GetField( 7 ) );
        fWeatherData.back()->AirPressure = atof( db_row->GetField( 11 ) );
        fWeatherData.back()->AirWindDirection = atof( db_row->GetField( 5 ) );
        fWeatherData.back()->AirWindSpeed = atof( db_row->GetField( 2 ) );
        fWeatherData.back()->AirWindGust = atof( db_row->GetField( 3 ) );
    }

    return true;
}


bool VRunStats::readDBSourceInfo( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    // read run data
    sprintf( c_query, "select * from tblObserving_Sources" );
    cout << c_query << endl;

    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    int fNRows = db_res->GetRowCount();
    for( unsigned int i = 0; i < 3; i++ )
    {
        cout << i << "\t" << db_res->GetFieldName( i ) << endl;
    }

    cout << "reading source info from DB" << endl;
    cout << "\t number of sources read form DB: " << fNRows << endl;

    string itemp;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();

        itemp = db_row->GetField( 0 );
        fSourceData[itemp] = new VDBSourceData();
        fSourceData[itemp]->SourceID = itemp;
        fSourceData[itemp]->fDec = atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi();
        fSourceData[itemp]->fRA  = atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi();
    }
    return true;
}

bool VRunStats::readDBRun_IDs( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }
    char c_query[1000];

    cout << "reading run IDs from DB" << endl;

    // read run data
    sprintf( c_query, "select run_id from tblRun_Info where db_start_time >= \"%s.000000\"  and db_start_time < \"%s.160000\"", fStartDate.c_str(), fStopDate.c_str() );
    cout << c_query << endl;

    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    int fNRows = db_res->GetRowCount();

    cout << "\t number of runs read form DB: " << fNRows << endl;

    string itemp;

    for( int i = 0; i < fNRows; i++ )
    {
        TSQLRow* db_row = db_res->Next();

        if( !db_row )
        {
            break;
        }

        // get run number
        itemp = db_row->GetField( 0 );
        if( i == 0 )
        {
            fMin_run_id = atoi( itemp.c_str() );
        }
        if( i == fNRows - 1 )
        {
            fMax_run_id = atoi( itemp.c_str() );
        }
    }
    cout << "run_id range: " << fMin_run_id << "\t" << fMax_run_id << endl;
    return true;
}


bool VRunStats::readDBRunInfo( TSQLServer* f_db )
{
    if( !f_db )
    {
        return false;
    }
    char c_query[1000];

    cout << "reading run info from DB" << endl;

    // read run data
    sprintf( c_query, "select * from tblRun_Info where db_start_time >= \"%s.000000\"  and db_start_time < \"%s.160000\"", fStartDate.c_str(), fStopDate.c_str() );
    cout << c_query << endl;

    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    int fNRows = db_res->GetRowCount();
    for( unsigned int i = 0; i < 19; i++ )
    {
        cout << i << "\t" << db_res->GetFieldName( i ) << endl;
    }

    cout << "\t number of runs read form DB: " << fNRows << endl;

    string itemp;
    unsigned int j_start = 0;
    unsigned int k_start = 0;
    unsigned int cs_start = 0;

    for( int i = 0; i < fNRows; i++ )
    {
        TSQLRow* db_row = db_res->Next();

        if( !db_row )
        {
            break;
        }

        cout << "ROW " << i << endl;
        if( db_row->GetField( 0 ) )
        {
            cout << "\t RRRR " << db_row->GetField( 0 ) << endl;
        }

        // all fields should be defined
        if( !db_row->GetField( 19 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 19 );
        // don't use laser or charge injection runs
        if( itemp == "NOSOURCE" )
        {
            continue;
        }
        // check if this run is an observing run
        if( !db_row->GetField( 1 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 1 );
        if( itemp != "observing" )
        {
            continue;
        }
        // don't use aborted runs
        if( !db_row->GetField( 3 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 3 );
        if( itemp == "aborted" )
        {
            continue;
        }
        // don't use runs which are started only
        itemp = db_row->GetField( 3 );
        if( itemp == "started" )
        {
            continue;
        }
        // don't use runs which are defined only
        itemp = db_row->GetField( 3 );
        if( itemp == "defined" )
        {
            continue;
        }
        // don't use runs which are prepared only
        itemp = db_row->GetField( 3 );
        if( itemp == "prepared" )
        {
            continue;
        }

        fRunData.push_back( new VDBRunData() );
        VDBRunData* i_RunData = fRunData.back();

        // get source coordinates
        double iRa = 0.;
        double iDec = 0.;
        itemp = db_row->GetField( 19 );
        i_RunData->SourceID = itemp;
        if( fSourceData.find( itemp ) != fSourceData.end() && fSourceData[itemp] )
        {
            i_RunData->RA = fSourceData[itemp]->fRA;
            i_RunData->Dec = fSourceData[itemp]->fDec;
        }
        else
        {
            cout << "No coordinates available: " << db_row->GetField( 0 ) << " " << itemp << endl;
            continue;
        }
        i_RunData->offsetRA = atof( db_row->GetField( 15 ) ) * 180. / TMath::Pi();
        i_RunData->offsetDec = atof( db_row->GetField( 16 ) ) * 180. / TMath::Pi();
        i_RunData->ConfigMask = atoi( db_row->GetField( 10 ) );
        i_RunData->offset_distance = atof( db_row->GetField( 17 ) );
        i_RunData->offset_angle = atof( db_row->GetField( 18 ) );

        // get galactic coordinates
        double i_b = 0.;
        double i_l = 0.;
        iRa  += i_RunData->offsetRA;
        iDec += i_RunData->offsetDec;

        VAstronometry::vlaEqgal( iRa / 180. * TMath::Pi(), iDec / 180. * TMath::Pi(), &i_l, &i_b );
        i_RunData->GalLong1958 = i_l * 180. / TMath::Pi();
        i_RunData->GalLat1958  = i_b * 180. / TMath::Pi();

        // calculate MJD, etc.

        int iMJD = 0;
        double iTime1 = 0.;
        double iTime2 = 0.;

        itemp = db_row->GetField( 6 );
        getDBMJDTime( itemp, iMJD, iTime1, true );
        i_RunData->start_time = itemp;
        i_RunData->StartMJD = ( double )iMJD + iTime1 / 86400;
        if( fDebug )
        {
            cout << "T1 " << iMJD << " " << iTime1 << " ";
        }

        itemp = db_row->GetField( 7 );
        if( fDebug )
        {
            cout << itemp << " " << flush;
        }
        getDBMJDTime( itemp, iMJD, iTime2, true );
        if( fDebug )
        {
            cout << "T2 " << iMJD << " " << iTime2 << " ";
        }
        i_RunData->stop_time = itemp;
        i_RunData->StopMJD = ( double )iMJD + iTime2 / 86400;

        i_RunData->Duration = iTime2 - iTime1;

        // get run number
        itemp = db_row->GetField( 0 );
        if( fDebug )
        {
            cout << "RUN " << itemp << " " << flush;
        }
        i_RunData->runNumber = atoi( itemp.c_str() );

        cout << "Filling run " << i_RunData->runNumber << endl;

        // get camera status
        cout << "\t camera " << endl;
        for( unsigned int k = cs_start; k < fCameraStatus.size(); k++ )
        {
            if( !fCameraStatus[k] )
            {
                continue;
            }

            if( fCameraStatus[k]->MJD > i_RunData->StartMJD && fCameraStatus[k]->MJD < i_RunData->StopMJD )
            {
                i_RunData->fCameraStatus.push_back( fCameraStatus[k] );
                cs_start = k;
            }
        }


        // get FIR data
        cout << "\t FIR " << endl;
        for( unsigned int k = k_start; k < fFIRData.size(); k++ )
        {
            if( !fFIRData[k] )
            {
                continue;
            }

            if( fFIRData[k]->MJD > i_RunData->StartMJD && fFIRData[k]->MJD < i_RunData->StopMJD )
            {
                i_RunData->fDBFIRData.push_back( fFIRData[k] );
                k_start = k;
            }
        }

        // get weather info
        cout << "\t weather " << endl;
        for( unsigned int j = j_start; j < fWeatherData.size(); j++ )
        {
            if( !fWeatherData[j] )
            {
                continue;
            }

            if( fWeatherData[j]->MJD > i_RunData->StartMJD && fWeatherData[j]->MJD < i_RunData->StopMJD )
            {
                i_RunData->fDBWeatherData.push_back( fWeatherData[j] );
                j_start = j;
            }
        }
        cout << "\t\t end weather" << endl;

        if( fDebug )
        {
            cout << endl;
        }
    }

    cout << "total number of runs in runlist: " << fRunData.size() << endl;

    return true;
}


bool VRunStats::getDBSourceCoordinates( TSQLServer* f_db, string iSource, double& iEVNTargetDec, double& iEVNTargetRA )
{
    if( !f_db )
    {
        return false;
    }

    char c_query[1000];

    sprintf( c_query, "select * from tblObserving_Sources where source_id like convert( _utf8 \'%s\' using latin1)", iSource.c_str() );
    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    TSQLRow* db_row = db_res->Next();

    if( db_row &&  db_row->GetField( 2 ) && db_row->GetField( 1 ) )
    {
        iEVNTargetDec = atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi();
        iEVNTargetRA = atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi();
    }
    else
    {
        iEVNTargetDec = 0.;
        iEVNTargetRA = 0.;
        return false;
    }
    return true;
}


bool VRunStats::getDBWeatherData( TSQLServer* f_db, string iTimeStamp, double& iWindSpeed, double& iWindGust, double& iWindDir, double& iTemp, double& iHumidity )
{
    if( !f_db )
    {
        return false;
    }

    // weather data is saved once per minute
    iTimeStamp.replace( 17, 2, "00" );

    char c_query[1000];

    sprintf( c_query, "select * from tblWeather_Status where timestamp = \"%s\" ", iTimeStamp.c_str() );
    TSQLResult* db_res = f_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }

    TSQLRow* db_row = db_res->Next();

    if( db_row )
    {
        iWindSpeed = atof( db_row->GetField( 2 ) );
        iWindGust = atof( db_row->GetField( 3 ) );
        iTemp = atof( db_row->GetField( 6 ) );
        iHumidity = atof( db_row->GetField( 7 ) );
        iWindDir = atof( db_row->GetField( 5 ) );
    }

    return true;
}


void VRunStats::getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip )
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


void VRunStats::print()
{
    for( unsigned int i = 0; i < fRunData.size(); i++ )
    {
        if( fRunData[i] )
        {
            cout << fRunData[i]->runNumber << " " << fRunData[i]->start_time << " " << fRunData[i]->stop_time << " ";
            //	    << fRunData[i]->AirHumidity;
            //	    cout << " " << fRunData[i]->AirTemperature << " " << fRunData[i]->AirWindSpeed << " " << fRunData[i]->AirWindGust << endl;
            cout << endl;
        }
    }
}


bool VRunStats::writeRootFile( string ifile )
{
    if( ifile.size() <= 0 )
    {
        return false;
    }

    TFile f( ifile.c_str(), "RECREATE" );
    if( f.IsZombie() )
    {
        return false;
    }

    int iRun;
    Char_t iRunSourceID[300];
    unsigned int iConfigMask;
    double iRunStartMJD;
    double iRunStoppMJD;
    double iRunDuration;
    double iRunRA;
    double iRunDec;
    double iRunRAoffset;
    double iRunDecoffset;
    double iRunoffset_distance;
    double iRunoffset_angle;
    double l;
    double b;
    // weather data
    int nWeatherPoints;
    const int nWeatherPoints_MAX = 10000;
    double weatherMJD[nWeatherPoints_MAX];
    double wind[nWeatherPoints_MAX];
    double windGust[nWeatherPoints_MAX];
    double windDir[nWeatherPoints_MAX];
    double Humidl[nWeatherPoints_MAX];
    double temp[nWeatherPoints_MAX];
    double pressure[nWeatherPoints_MAX];
    // HV data
    double iHV_MJD_startMJD;
    double iHV_MJD_stopMJD;
    unsigned int iHV_telescope = 0;
    unsigned int iHV_channel = 0;
    double iHV_HV;

    // camera status
    const int nCSPoints_MAX = 10000;
    int nCameraStatusPoints = 0;
    double iCS_MJD[nCSPoints_MAX];
    unsigned int iCS_telescope[nCSPoints_MAX];
    double iCS_temp1[nCSPoints_MAX];
    double iCS_temp2[nCSPoints_MAX];
    double iCS_hum1[nCSPoints_MAX];
    double iCS_hum2[nCSPoints_MAX];
    double iCS_light1[nCSPoints_MAX];
    double iCS_light2[nCSPoints_MAX];

    // FIR data
    const int nFIRPoints_MAX = 10000;
    int FIR_NPoints = 0;
    unsigned int FIR_tel_ID[nFIRPoints_MAX];
    double FIR_MJD[nFIRPoints_MAX];
    double FIR_ambient_temp[nFIRPoints_MAX];
    double FIR_radiant_sky_temp[nFIRPoints_MAX];
    double FIR_radiant_sky_temp_cor[nFIRPoints_MAX];

    // DQM data
    int fieldVersion = 0;
    double StartMJD = 0.;
    int NTel = 0;
    unsigned int TelBitMask = 0;

    double WobbleOff_estim_deg = 0.;
    double WobbleDir_estim_deg = 0.;
    double El_mean_deg = 0.;
    double Az_mean_deg = 0.;
    double L3_runTime_s = 0.;
    double L3_lifeTime_s = 0.;
    double L3_rate = 0.;
    double Rate_mean = 0.;
    double Rate_RMS = 0.;
    double Rate_Chi2 = 0.;

    double L3_threshold = 0.;

    double L2_threshold[4];
    double L2_threshold_mV[4];
    double L2_threshold_Spread_mV[4];
    double L2_throughput[4];

    double GPS_runtime_min = 0.;
    double GPS_lifetime_min = 0.;

    double MoonEl_mean_deg = 0.;
    double MoonPhase_mean = 0.;
    double MoonSeparation_deg[4];

    double RA_mean[4];
    double RA_RMS[4];
    double Dec_mean[4];
    double Dec_RMS[4];

    double L1_rate[4];
    double Median_current_uA[4];
    double L2_rate[4];
    double MissingEvents[4];
    double Muons_N[4];
    double Muons_MeanCharge[4];
    double Muons_RMSCharge[4];
    double PMT_GainLaser[4];

    double SuppressedChannels_Num[4];

    double PedVar_mean[4];

    int LaserRun[4];

    double FIR0_mean = 0.;
    double FIR0_RMS = 0.;
    double FIR_mean[4];
    double FIR_RMS[4];
    for( unsigned int i = 0; i < 4; i++ )
    {
        L2_threshold[i] = 0.;
        L2_threshold_mV[i] = 0.;
        L2_threshold_Spread_mV[i] = 0.;
        L2_throughput [i] = 0.;
        MoonSeparation_deg[i] = 0.;
        RA_mean[i] = 0.;
        RA_RMS[i] = 0.;
        Dec_mean[i] = 0.;
        Dec_RMS[i] = 0.;
        L1_rate[i] = 0.;
        Median_current_uA[i] = 0.;
        L2_rate[i] = 0.;
        MissingEvents[i] = 0.;
        Muons_N[i] = 0.;
        Muons_MeanCharge[i] = 0.;
        Muons_RMSCharge[i] = 0.;
        PMT_GainLaser[i] = 0.;
        SuppressedChannels_Num[i] = 0.;
        PedVar_mean[i] = 0.;
        LaserRun[i] = 0;
        FIR_mean[i] = 0.;
        FIR_RMS[i] = 0.;
    }
    unsigned int iStatus = 0;
    double iUsableDuration = 0.;
    Char_t iComment[300];

    cout << "writing data to " << f.GetName() << endl;

    char hname[2000];

    sprintf( hname, "DB source table" );
    TTree s( "fSourceTable", hname );
    s.Branch( "sourceID", &iRunSourceID, "sourceID/C" );
    s.Branch( "ra", &iRunRA, "ra/D" );
    s.Branch( "dec", &iRunDec, "dec/D" );

    map< string, VDBSourceData* >::iterator i_iter;

    for( i_iter = fSourceData.begin(); i_iter != fSourceData.end(); ++i_iter )
    {
        VDBSourceData* i_temp = ( VDBSourceData* )i_iter->second;
        if( i_temp )
        {
            sprintf( iRunSourceID, "%s", i_temp->SourceID.c_str() );
            iRunRA = i_temp->fRA;
            iRunDec =  i_temp->fDec;

            s.Fill();
        }
    }
    s.Write();

    sprintf( hname, "DB Default HV" );
    TTree i_DHV( "fDefaultHV", hname );
    i_DHV.Branch( "startMJD", &iHV_MJD_startMJD, "startMJD/D" );
    i_DHV.Branch( "stopMJD", &iHV_MJD_stopMJD, "stopMJD/D" );
    i_DHV.Branch( "telescope", &iHV_telescope, "telescope/i" );
    i_DHV.Branch( "channel", &iHV_channel, "channel/i" );
    i_DHV.Branch( "defaultHV", &iHV_HV, "defaultHV/D" );

    for( unsigned int i = 0; i < fdefaultHVEvndispData.size(); i++ )
    {
        if( fdefaultHVEvndispData[i] )
        {
            iHV_MJD_startMJD = fdefaultHVEvndispData[i]->startMJD;
            iHV_MJD_stopMJD  = fdefaultHVEvndispData[i]->stopMJD;
            iHV_telescope    = fdefaultHVEvndispData[i]->telescope + 1;
            iHV_channel      = fdefaultHVEvndispData[i]->channel;
            iHV_HV   = fdefaultHVEvndispData[i]->default_HV;
            i_DHV.Fill();
        }
    }
    i_DHV.Write();

    sprintf( hname, "DB HV,  %s to %s", fStartDate.c_str(), fStopDate.c_str() );
    TTree i_MHV( "fHV", hname );
    i_MHV.Branch( "startMJD", &iHV_MJD_startMJD, "startMJD/D" );
    i_MHV.Branch( "stopMJD", &iHV_MJD_stopMJD, "stopMJD/D" );
    i_MHV.Branch( "telescope", &iHV_telescope, "telescope/i" );
    i_MHV.Branch( "channel", &iHV_channel, "channel/i" );
    i_MHV.Branch( "HV", &iHV_HV, "HV/D" );

    for( unsigned int i = 0; i < fHVEvndispData.size(); i++ )
    {
        if( fHVEvndispData[i] )
        {
            iHV_MJD_startMJD = fHVEvndispData[i]->startMJD;
            iHV_MJD_stopMJD  = fHVEvndispData[i]->stopMJD;
            iHV_telescope    = fHVEvndispData[i]->telescope + 1;
            iHV_channel      = fHVEvndispData[i]->channel;
            iHV_HV   = fHVEvndispData[i]->HV;
            i_MHV.Fill();
        }
    }
    i_MHV.Write();

    sprintf( hname, "DB FIR data,  %s to %s", fStartDate.c_str(), fStopDate.c_str() );
    TTree r( "fFIRTable", hname );
    r.Branch( "MJD", FIR_MJD, "MJD/D" );
    r.Branch( "TelID", FIR_tel_ID, "TelID/i" );
    r.Branch( "ambient_temp", FIR_ambient_temp, "ambient_temp/D" );
    r.Branch( "radiant_sky_temp", FIR_radiant_sky_temp, "radiant_sky_temp/D" );
    r.Branch( "radiant_sky_temp_cor", FIR_radiant_sky_temp_cor, "radiant_sky_temp_cor/D" );

    for( unsigned int i = 0; i < fFIRData.size(); i++ )
    {
        if( fFIRData[i] )
        {
            FIR_MJD[0] = fFIRData[i]->MJD;
            FIR_tel_ID[0] = fFIRData[i]->tel_ID;
            FIR_ambient_temp[0] = fFIRData[i]->ambient_temp;
            FIR_radiant_sky_temp[0] = fFIRData[i]->radiant_sky_temp;
            FIR_radiant_sky_temp_cor[0] = fFIRData[i]->radiant_sky_temp_cor;
            r.Fill();
        }
    }
    r.Write();

    sprintf( hname, "DB camera status,  %s to %s", fStartDate.c_str(), fStopDate.c_str() );
    TTree cs( "fCamerStatus", hname );
    cs.Branch( "MJD", iCS_MJD, "MJD/D" );
    cs.Branch( "telescope", iCS_telescope, "telescope/i" );
    cs.Branch( "temp1", iCS_temp1, "temp1/D" );
    cs.Branch( "temp2", iCS_temp2, "temp2/D" );
    cs.Branch( "hum1", iCS_hum1, "hum1/D" );
    cs.Branch( "hum2", iCS_hum2, "hum2/D" );
    cs.Branch( "light1", iCS_light1, "light1/D" );
    cs.Branch( "light2", iCS_light2, "light2/D" );

    for( unsigned int i = 0; i < fCameraStatus.size(); i++ )
    {
        if( fCameraStatus[i] )
        {
            iCS_MJD[0] = fCameraStatus[i]->MJD;
            iCS_telescope[0] = fCameraStatus[i]->telescope + 1;
            iCS_temp1[0] = fCameraStatus[i]->temp1;
            iCS_temp2[0] = fCameraStatus[i]->temp2;
            iCS_hum1[0] = fCameraStatus[i]->hum1;
            iCS_hum2[0] = fCameraStatus[i]->hum2;
            iCS_light1[0] = fCameraStatus[i]->light1;
            iCS_light2[0] = fCameraStatus[i]->light2;

            cs.Fill();
        }
    }
    cs.Write();


    sprintf( hname, "DB weather data,  %s to %s", fStartDate.c_str(), fStopDate.c_str() );
    TTree w( "fWeatherTable", hname );
    w.Branch( "MJD", weatherMJD, "MJD/D" );
    w.Branch( "temperature", temp, "temperature/D" );
    w.Branch( "pressure", pressure, "pressure/D" );
    w.Branch( "windSpeedAv", wind, "windSpeedAv/D" );
    w.Branch( "windSpeedMax", windGust, "windSpeedMax/D" );
    w.Branch( "windDirection", windDir, "windDirection/D" );
    w.Branch( "humidity", Humidl, "humidity/D" );

    for( unsigned int i = 0; i < fWeatherData.size(); i++ )
    {
        if( fWeatherData[i] )
        {
            weatherMJD[0] = fWeatherData[i]->MJD;
            temp[0] = fWeatherData[i]->AirTemperature;
            pressure[0] = fWeatherData[i]->AirPressure;
            wind[0] = fWeatherData[i]->AirWindSpeed;
            windGust[0] = fWeatherData[i]->AirWindGust;
            windDir[0] = fWeatherData[i]->AirWindDirection;
            Humidl[0] = fWeatherData[i]->AirHumidity;

            w.Fill();
        }
    }
    w.Write();

    // tree with all DQM information per run
    sprintf( hname, "DB entries, %s to %s", fStartDate.c_str(), fStopDate.c_str() );
    TTree t( "fRunTable", hname );
    t.Branch( "runNumber", &iRun, "runNumber/I" );
    t.Branch( "sourceID", &iRunSourceID, "sourceID/C" );
    t.Branch( "configmask", &iConfigMask, "configmask/i" );
    t.Branch( "MJDstart", &iRunStartMJD, "MJDstart/D" );
    t.Branch( "MJDstop", &iRunStoppMJD, "MJDstop/D" );
    t.Branch( "runLength", &iRunDuration, "runLength/D" );
    t.Branch( "ra", &iRunRA, "ra/D" );
    t.Branch( "dec", &iRunDec, "dec/D" );
    t.Branch( "raOffset", &iRunRAoffset, "raOffset/D" );
    t.Branch( "decOffset", &iRunDecoffset, "decOffset/D" );
    t.Branch( "offset_distance", &iRunoffset_distance, "offset_distance/D" );
    t.Branch( "offset_angle", &iRunoffset_angle, "offset_angle/D" );
    t.Branch( "l", &l, "l/D" );
    t.Branch( "b", &b, "b/D" );
    t.Branch( "nWeatherPoints", &nWeatherPoints, "nWeatherPoints/I" );
    t.Branch( "weatherMJD", weatherMJD, "weatherMJD[nWeatherPoints]/D" );
    t.Branch( "temperature", temp, "temperature[nWeatherPoints]/D" );
    t.Branch( "windSpeedAv", wind, "windSpeedAv[nWeatherPoints]/D" );
    t.Branch( "windSpeedMax", windGust, "windSpeedMax[nWeatherPoints]/D" );
    t.Branch( "windDirection", windDir, "windDirection[nWeatherPoints]/D" );
    t.Branch( "humidity", Humidl, "humidity[nWeatherPoints]/D" );

    t.Branch( "nFIRPoints", &FIR_NPoints, "nFIRPoints/I" );
    t.Branch( "FIR_tel_ID", FIR_tel_ID, "FIR_tel_ID[nFIRPoints]/i" );
    t.Branch( "FIR_MJD", FIR_MJD, "FIR_MJD[nFIRPoints]/D" );
    t.Branch( "FIR_ambient_temp", FIR_ambient_temp, "FIR_ambient_temp[nFIRPoints]/D" );
    t.Branch( "FIR_radiant_sky_temp", FIR_radiant_sky_temp, "FIR_radiant_sky_temp[nFIRPoints]/D" );
    t.Branch( "FIR_radiant_sky_temp_cor", FIR_radiant_sky_temp_cor, "FIR_radiant_sky_temp_cor[nFIRPoints]/D" );

    t.Branch( "nCameraStatusPoints", &nCameraStatusPoints, "nCameraStatusPoints/I" );
    t.Branch( "CameraStatus_tel_ID", iCS_telescope, "CameraStatus_tel_ID[nCameraStatusPoints]/I" );
    t.Branch( "CameraStatus_MJD", iCS_MJD, "CameraStatus_MJD[nCameraStatusPoints]/D" );
    t.Branch( "CameraStatus_temp1", iCS_temp1, "CameraStatus_temp1[nCameraStatusPoints]/D" );
    t.Branch( "CameraStatus_temp2", iCS_temp2, "CameraStatus_temp2[nCameraStatusPoints]/D" );
    t.Branch( "CameraStatus_hum1", iCS_hum1, "CameraStatus_hum1[nCameraStatusPoints]/D" );
    t.Branch( "CameraStatus_hum2", iCS_hum2, "CameraStatus_hum2[nCameraStatusPoints]/D" );
    t.Branch( "CameraStatus_light1", iCS_light1, "CameraStatus_light1[nCameraStatusPoints]/D" );
    t.Branch( "CameraStatus_light2", iCS_light2, "CameraStatus_light2[nCameraStatusPoints]/D" );

    t.Branch( "fieldVersion", &fieldVersion, "fieldVersion/I" );
    t.Branch( "StartMJD", &StartMJD, "StartMJD/D" );
    t.Branch( "NTel", &NTel, "NTel/I" );
    t.Branch( "TelBitMask", &TelBitMask, "TelBitMask/i" );
    t.Branch( "El_mean_deg", &El_mean_deg, "El_mean_deg/D" );
    t.Branch( "Az_mean_deg", &Az_mean_deg, "Az_mean_deg/D" );
    t.Branch( "offset_distance_estim_deg", &WobbleOff_estim_deg, "offset_distance_estim_deg/D" );
    t.Branch( "offset_angle_estim_deg", &WobbleDir_estim_deg, "offset_angle_estim_deg/D" );
    t.Branch( "L3_runTime_s", &L3_runTime_s, "L3_runTime_s/D" );
    t.Branch( "L3_lifeTime_s", &L3_lifeTime_s, "L3_lifeTime_s/D" );
    t.Branch( "L3_rate", &L3_rate, "L3_rate/D" );
    t.Branch( "Rate_mean", &Rate_mean, "Rate_mean/D" );
    t.Branch( "Rate_RMS", &Rate_RMS, "Rate_RMS/D" );
    t.Branch( "Rate_Chi2", &Rate_Chi2, "Rate_Chi2/D" );
    t.Branch( "L3_threshold", &L3_threshold, "L3_threshold/D" );
    t.Branch( "L2_threshold", L2_threshold, "L2_threshold[4]/D" );
    t.Branch( "L2_threshold_mV", L2_threshold_mV, "L2_threshold_mV[4]/D" );
    t.Branch( "L2_threshold_Spread_mV", L2_threshold_Spread_mV, "L2_threshold_Spread_mV[4]/D" );
    t.Branch( "L2_throughput", L2_throughput, "L2_throughput[4]/D" );
    t.Branch( "GPS_runtime_min", &GPS_runtime_min, "GPS_runtime_min/D" );
    t.Branch( "GPS_lifetime_min", &GPS_lifetime_min, "GPS_lifetime_min/D" );
    t.Branch( "MoonEl_mean_deg", &MoonEl_mean_deg, "MoonEl_mean_deg/D" );
    t.Branch( "MoonPhase_mean", &MoonPhase_mean, "MoonPhase_mean/D" );
    t.Branch( "MoonSeparation_deg", MoonSeparation_deg, "MoonSeparation_deg[4]/D" );
    t.Branch( "RA_mean", RA_mean, "RA_mean[4]/D" );
    t.Branch( "RA_RMS", RA_RMS, "RA_RMS[4]/D" );
    t.Branch( "Dec_mean", Dec_mean, "Dec_mean[4]/D" );
    t.Branch( "Dec_RMS", Dec_RMS, "Dec_RMS[4]/D" );
    t.Branch( "L1_rate", L1_rate, "L1_rate[4]/D" );
    t.Branch( "Median_current_uA", Median_current_uA, "Median_current_uA[4]/D" );
    t.Branch( "L2_rate", L2_rate, "L2_rate[4]/D" );
    t.Branch( "MissingEvents", MissingEvents, "MissingEvents[4]/D" );
    t.Branch( "Muons_N", Muons_N, "Muons_N[4]/D" );
    t.Branch( "Muons_MeanCharge", Muons_MeanCharge, "Muons_MeanCharge[4]/D" );
    t.Branch( "Muons_RMSCharge", Muons_RMSCharge, "Muons_RMSCharge[4]/D" );
    t.Branch( "PMT_GainLaser", PMT_GainLaser, "PMT_GainLaser[4]/D" );
    t.Branch( "SuppressedChannels_Num", SuppressedChannels_Num, "SuppressedChannels_Num[4]/D" );
    t.Branch( "PedVar_mean", PedVar_mean, "PedVar_mean[4]/D" );
    t.Branch( "LaserRun", LaserRun, "LaserRun[4]/I" );
    t.Branch( "FIR0_mean", &FIR0_mean, "FIR0_mean/D" );
    t.Branch( "FIR0_RMS", &FIR0_RMS, "FIR0_RMS/D" );
    t.Branch( "FIR_mean", FIR_mean, "FIR_mean[4]/D" );
    t.Branch( "FIR_RMS", FIR_RMS, "FIR_RMS[4]/D" );
    t.Branch( "status", &iStatus, "status/i" );
    t.Branch( "usable_duration", &iUsableDuration, "usable_duration/D" );
    t.Branch( "comment", &iComment, "comment/C" );

    unsigned int k_start = 0;
    /////////////////////
    // loop over all runs
    for( unsigned int i = 0; i < fRunData.size(); i++ )
    {
        iRun = fRunData[i]->runNumber;
        if( fRunData[i]->SourceID.size() < 300 )
        {
            sprintf( iRunSourceID, "%s", fRunData[i]->SourceID.c_str() );
        }
        else
        {
            sprintf( iRunSourceID, "%s", fRunData[i]->SourceID.substr( 0, 299 ).c_str() );
        }
        iConfigMask = fRunData[i]->ConfigMask;
        iRunStartMJD = fRunData[i]->StartMJD;
        iRunStoppMJD = fRunData[i]->StopMJD;
        iRunDuration = fRunData[i]->Duration;
        iRunRA = fRunData[i]->RA;
        iRunDec = fRunData[i]->Dec;
        iRunRAoffset = fRunData[i]->offsetRA;
        iRunDecoffset = fRunData[i]->offsetDec;
        iRunoffset_distance = fRunData[i]->offset_distance;
        iRunoffset_angle = fRunData[i]->offset_angle;
        l = fRunData[i]->GalLong1958;
        b = fRunData[i]->GalLat1958;

        // analysis comments
        iUsableDuration = iRunDuration;
        for( unsigned int j = 0; j < fAnalysisComments.size(); j++ )
        {
            if( iRun == fAnalysisComments[j]->runNumber )
            {
                iStatus = fAnalysisComments[j]->status;
                iUsableDuration = fAnalysisComments[j]->usable_duration;
                if( fAnalysisComments[j]->comment.size() < 300 )
                {
                    sprintf( iComment, "%s", fAnalysisComments[j]->comment.c_str() );
                }
                else
                {
                    sprintf( iComment, "%s", fAnalysisComments[j]->comment.substr( 0, 299 ).c_str() );
                }
                break;
            }
            else
            {
                iStatus = 0;
                sprintf( iComment, " " );
            }
        }

        // camera status
        nCameraStatusPoints = ( int )fRunData[i]->fCameraStatus.size();
        if( nCameraStatusPoints > nCSPoints_MAX )
        {
            continue;
        }
        for( int j = 0; j < nCameraStatusPoints; j++ )
        {
            if( fRunData[i]->fCameraStatus[j] )
            {
                iCS_MJD[j] = fRunData[i]->fCameraStatus[j]->MJD;
                iCS_telescope[j] = fRunData[i]->fCameraStatus[j]->telescope + 1;
                iCS_temp1[j] = fRunData[i]->fCameraStatus[j]->temp1;
                iCS_temp2[j] = fRunData[i]->fCameraStatus[j]->temp2;
                iCS_hum1[j] = fRunData[i]->fCameraStatus[j]->hum1;
                iCS_hum2[j] = fRunData[i]->fCameraStatus[j]->hum2;
                iCS_light1[j] = fRunData[i]->fCameraStatus[j]->light1;
                iCS_light2[j] = fRunData[i]->fCameraStatus[j]->light2;
            }
        }

        // FIR data
        FIR_NPoints = ( int )fRunData[i]->fDBFIRData.size();
        if( FIR_NPoints > nFIRPoints_MAX )
        {
            continue;
        }
        for( int j = 0; j < FIR_NPoints; j++ )
        {
            if( fRunData[i]->fDBFIRData[j] )
            {
                FIR_tel_ID[j] = fRunData[i]->fDBFIRData[j]->tel_ID;
                FIR_MJD[j] = fRunData[i]->fDBFIRData[j]->MJD;
                FIR_ambient_temp[j] = fRunData[i]->fDBFIRData[j]->ambient_temp;
                FIR_radiant_sky_temp[j] = fRunData[i]->fDBFIRData[j]->radiant_sky_temp;
                FIR_radiant_sky_temp_cor[j] = fRunData[i]->fDBFIRData[j]->radiant_sky_temp_cor;
            }
        }

        // weather data
        nWeatherPoints = ( int )fRunData[i]->fDBWeatherData.size();
        if( nWeatherPoints > nWeatherPoints_MAX )
        {
            continue;
        }

        for( int j = 0; j < nWeatherPoints; j++ )
        {
            if( fRunData[i]->fDBWeatherData[j] )
            {
                weatherMJD[j] = fRunData[i]->fDBWeatherData[j]->MJD;
                temp[j] = fRunData[i]->fDBWeatherData[j]->AirTemperature;
                pressure[i] = fRunData[i]->fDBWeatherData[j]->AirPressure;
                wind[j] = fRunData[i]->fDBWeatherData[j]->AirWindSpeed;
                windGust[j] = fRunData[i]->fDBWeatherData[j]->AirWindGust;
                windDir[j] = fRunData[i]->fDBWeatherData[j]->AirWindDirection;
                Humidl[j] = fRunData[i]->fDBWeatherData[j]->AirHumidity;
            }
        }

        for( unsigned int k = k_start; k < fDQMData.size(); k++ )
        {
            if( fDQMData[k] && fDQMData[k]->runNumber == fRunData[i]->runNumber )
            {
                fieldVersion = fDQMData[k]->fieldVersion;
                StartMJD = fDQMData[k]->StartMJD;
                NTel = fDQMData[k]->NTel;
                TelBitMask = fDQMData[k]->TelBitMask;
                WobbleOff_estim_deg = fDQMData[k]->WobbleOff_estim_deg;
                WobbleDir_estim_deg = fDQMData[k]->WobbleDir_estim_deg;
                El_mean_deg = fDQMData[k]->El_mean_deg;
                Az_mean_deg = fDQMData[k]->Az_mean_deg;
                L3_runTime_s = fDQMData[k]->L3_runTime_s;
                L3_lifeTime_s = fDQMData[k]->L3_lifeTime_s;
                L3_rate = fDQMData[k]->L3_rate;
                Rate_mean = fDQMData[k]->Rate_mean;
                Rate_RMS = fDQMData[k]->Rate_RMS;
                Rate_Chi2 = fDQMData[k]->Rate_Chi2;
                L3_threshold = fDQMData[k]->L3_threshold;
                GPS_runtime_min = fDQMData[k]->GPS_runtime_min;
                GPS_lifetime_min = fDQMData[k]->GPS_lifetime_min;
                MoonEl_mean_deg = fDQMData[k]->MoonEl_mean_deg;
                MoonPhase_mean = fDQMData[k]->MoonPhase_mean;
                for( unsigned int l = 0; l < 4; l++ )
                {
                    L2_threshold[l] = fDQMData[k]->L2_threshold[l];
                    L2_threshold_mV[l] = fDQMData[k]->L2_threshold_mV[l];
                    L2_threshold_Spread_mV[l] = fDQMData[k]->L2_threshold_Spread_mV[l];
                    L2_throughput[l] = fDQMData[k]->L2_throughput[l];
                    MoonSeparation_deg[l] = fDQMData[k]->MoonSeparation_deg[l];
                    RA_mean[l] = fDQMData[k]->RA_mean[l];
                    RA_RMS[l] = fDQMData[k]->RA_RMS[l];
                    Dec_mean[l] = fDQMData[k]->Dec_mean[l];
                    Dec_RMS[l] = fDQMData[k]->Dec_RMS[l];
                    L1_rate[l] = fDQMData[k]->L1_rate[l];
                    Median_current_uA[l]  = fDQMData[k]->Median_current_uA[l];
                    L2_rate[l] = fDQMData[k]->L2_rate[l];
                    MissingEvents[l] = fDQMData[k]->MissingEvents[l];
                    Muons_N[l] = fDQMData[k]->Muons_N[l];
                    Muons_MeanCharge[l] = fDQMData[k]->Muons_MeanCharge[l];
                    Muons_RMSCharge[l] = fDQMData[k]->Muons_RMSCharge[l];
                    PMT_GainLaser[l] = fDQMData[k]->PMT_GainLaser[l];
                    SuppressedChannels_Num[l] = fDQMData[k]->SuppressedChannels_Num[l];
                    PedVar_mean[l] = fDQMData[k]->PedVar_mean[l];
                    LaserRun[l] = fDQMData[k]->LaserRun[l];
                    FIR_mean[l] = fDQMData[k]->FIR_mean[l];
                    FIR_RMS[l] = fDQMData[k]->FIR_RMS[l];
                }
                FIR0_mean = fDQMData[k]->FIR0_mean;
                FIR0_RMS = fDQMData[k]->FIR0_RMS;

                k_start = k;
                break;
            }
            fieldVersion = -1;
        }

        t.Fill();
    }
    t.Write();

    f.Close();

    return true;
}


void VRunStats::printDQM( unsigned int iRun )
{
    for( unsigned int i = 0; i < fDQMData.size(); i++ )
    {
        if( fDQMData[i] && fDQMData[i]->runNumber == iRun )
        {
            cout << fDQMData[i]->runNumber << "\t";
            cout << fDQMData[i]->fieldVersion << "\t";
            cout << fDQMData[i]->StartMJD << "\t";
            cout << endl;
            break;
        }
    }
}


bool VRunStats::readDQMData( string iDQMFile )
{
    ifstream is;
    is.open( iDQMFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "VRunStats::readDQMData: DQM file not found: " << iDQMFile << endl;
        return false;
    }
    string iLine;
    string itemp;

    cout << "reading DQM file " << iDQMFile << endl;

    int z = 0;
    while( getline( is, iLine ) )
    {
        if( iLine.size() < 1 )
        {
            continue;
        }
        if( iLine.substr( 0, 1 ) == "#" )
        {
            continue;
        }

        fDQMData.push_back( new VDQMRunInfo() );

        istringstream is_stream( iLine );

        is_stream >> fDQMData.back()->runNumber;
        is_stream >> fDQMData.back()->fieldVersion;
        is_stream >> itemp;
        is_stream >> itemp;
        is_stream >> fDQMData.back()->StartMJD;
        is_stream >> fDQMData.back()->SourceID;
        is_stream >> fDQMData.back()->ObsMode;
        // get offset angle and distance
        if( fDQMData.back()->ObsMode.find( "wobble" ) != string::npos )
        {
            fDQMData.back()->WobbleOff_estim_deg = atof( fDQMData.back()->ObsMode.substr( 6,  fDQMData.back()->ObsMode.find( "@" ) ).c_str() );
            fDQMData.back()->WobbleDir_estim_deg = atof( fDQMData.back()->ObsMode.substr( fDQMData.back()->ObsMode.find( "@" ) + 1, fDQMData.back()->ObsMode.size() ).c_str() );
        }
        else
        {
            fDQMData.back()->WobbleOff_estim_deg = 0.;
            fDQMData.back()->WobbleDir_estim_deg = 0.;
        }
        is_stream >> fDQMData.back()->El_mean_deg;
        is_stream >> fDQMData.back()->Az_mean_deg;
        // 10
        is_stream >> fDQMData.back()->L3_runTime_s;
        is_stream >> fDQMData.back()->L3_lifeTime_s;
        is_stream >> fDQMData.back()->L3_rate;
        is_stream >> fDQMData.back()->Rate_mean;
        is_stream >> fDQMData.back()->Rate_RMS;
        is_stream >> fDQMData.back()->Rate_Chi2;
        is_stream >> fDQMData.back()->L3_threshold;
        is_stream >> fDQMData.back()->L2_threshold[0];
        is_stream >> fDQMData.back()->L2_threshold_mV[0];
        is_stream >> fDQMData.back()->L2_threshold_Spread_mV[0];
        // 20
        is_stream >> itemp;
        is_stream >> itemp;
        is_stream >> fDQMData.back()->L2_throughput[0];
        is_stream >> fDQMData.back()->L2_threshold[1];
        is_stream >> fDQMData.back()->L2_threshold_mV[1];
        is_stream >> fDQMData.back()->L2_threshold_Spread_mV[1];
        is_stream >> itemp;
        is_stream >> itemp;
        is_stream >> fDQMData.back()->L2_throughput[1];
        is_stream >> fDQMData.back()->L2_threshold[2];
        // 30
        is_stream >> fDQMData.back()->L2_threshold_mV[2];
        is_stream >> fDQMData.back()->L2_threshold_Spread_mV[2];
        is_stream >> itemp;
        is_stream >> itemp;
        is_stream >> fDQMData.back()->L2_throughput[2];
        is_stream >> fDQMData.back()->L2_threshold[3];
        is_stream >> fDQMData.back()->L2_threshold_mV[3];
        is_stream >> fDQMData.back()->L2_threshold_Spread_mV[3];
        is_stream >> itemp;
        is_stream >> itemp;
        // 40
        is_stream >> fDQMData.back()->L2_throughput[3];

        if( fDQMData.back()->fieldVersion < 1 )
        {
            continue;
        }

        is_stream >> fDQMData.back()->GPS_runtime_min;
        is_stream >> fDQMData.back()->GPS_lifetime_min;
        is_stream >> itemp;
        is_stream >> fDQMData.back()->MoonEl_mean_deg;
        is_stream >> fDQMData.back()->MoonPhase_mean;

        is_stream >> fDQMData.back()->RA_mean[0];
        is_stream >> fDQMData.back()->Dec_mean[0];
        is_stream >> fDQMData.back()->RA_RMS[0];
        is_stream >> fDQMData.back()->Dec_RMS[0];
        // 50
        is_stream >> fDQMData.back()->L1_rate[0];
        is_stream >> fDQMData.back()->Median_current_uA[0];
        is_stream >> fDQMData.back()->L2_rate[0];
        is_stream >> fDQMData.back()->MissingEvents[0];
        is_stream >> fDQMData.back()->Muons_N[0];
        is_stream >> fDQMData.back()->Muons_MeanCharge[0];
        is_stream >> fDQMData.back()->Muons_RMSCharge[0];
        is_stream >> fDQMData.back()->PMT_GainLaser[0];
        is_stream >> fDQMData.back()->MoonSeparation_deg[0];
        is_stream >> fDQMData.back()->SuppressedChannels_Num[0];
        // 60
        is_stream >> fDQMData.back()->RA_mean[1];
        is_stream >> fDQMData.back()->Dec_mean[1];
        is_stream >> fDQMData.back()->RA_RMS[1];
        is_stream >> fDQMData.back()->Dec_RMS[1];
        is_stream >> fDQMData.back()->L1_rate[1];
        is_stream >> fDQMData.back()->Median_current_uA[1];
        is_stream >> fDQMData.back()->L2_rate[1];
        is_stream >> fDQMData.back()->MissingEvents[1];
        is_stream >> fDQMData.back()->Muons_N[1];
        is_stream >> fDQMData.back()->Muons_MeanCharge[1];
        // 70
        is_stream >> fDQMData.back()->Muons_RMSCharge[1];
        is_stream >> fDQMData.back()->PMT_GainLaser[1];
        is_stream >> fDQMData.back()->MoonSeparation_deg[1];
        is_stream >> fDQMData.back()->SuppressedChannels_Num[1];
        is_stream >> fDQMData.back()->RA_mean[2];
        is_stream >> fDQMData.back()->Dec_mean[2];
        is_stream >> fDQMData.back()->RA_RMS[2];
        is_stream >> fDQMData.back()->Dec_RMS[2];
        is_stream >> fDQMData.back()->L1_rate[2];
        is_stream >> fDQMData.back()->Median_current_uA[2];
        // 80
        is_stream >> fDQMData.back()->L2_rate[2];
        is_stream >> fDQMData.back()->MissingEvents[2];
        is_stream >> fDQMData.back()->Muons_N[2];
        is_stream >> fDQMData.back()->Muons_MeanCharge[2];
        is_stream >> fDQMData.back()->Muons_RMSCharge[2];
        is_stream >> fDQMData.back()->PMT_GainLaser[2];
        is_stream >> fDQMData.back()->MoonSeparation_deg[2];
        is_stream >> fDQMData.back()->SuppressedChannels_Num[2];

        is_stream >> fDQMData.back()->RA_mean[3];
        is_stream >> fDQMData.back()->Dec_mean[3];
        // 90
        is_stream >> fDQMData.back()->RA_RMS[3];
        is_stream >> fDQMData.back()->Dec_RMS[3];
        is_stream >> fDQMData.back()->L1_rate[3];
        is_stream >> fDQMData.back()->Median_current_uA[3];
        is_stream >> fDQMData.back()->L2_rate[3];
        is_stream >> fDQMData.back()->MissingEvents[3];
        is_stream >> fDQMData.back()->Muons_N[3];
        is_stream >> fDQMData.back()->Muons_MeanCharge[3];
        is_stream >> fDQMData.back()->Muons_RMSCharge[3];
        is_stream >> fDQMData.back()->PMT_GainLaser[3];
        // 100
        is_stream >> fDQMData.back()->MoonSeparation_deg[3];
        is_stream >> fDQMData.back()->SuppressedChannels_Num[3];

        if( fDQMData.back()->fieldVersion < 2 )
        {
            continue;
        }

        is_stream >> fDQMData.back()->NTel;
        is_stream >> fDQMData.back()->TelBitMask;
        is_stream >> itemp;
        is_stream >> itemp;
        is_stream >> fDQMData.back()->FIR0_mean;
        is_stream >> fDQMData.back()->FIR0_RMS;

        is_stream >> fDQMData.back()->LaserRun[0];
        is_stream >> fDQMData.back()->PedVar_mean[0];
        // 110
        is_stream >> itemp;
        is_stream >> fDQMData.back()->FIR_mean[0];
        is_stream >> fDQMData.back()->FIR_RMS[0];
        is_stream >> fDQMData.back()->LaserRun[1];
        is_stream >> fDQMData.back()->PedVar_mean[1];
        is_stream >> itemp;
        is_stream >> fDQMData.back()->FIR_mean[1];
        is_stream >> fDQMData.back()->FIR_RMS[1];
        is_stream >> fDQMData.back()->LaserRun[2];
        is_stream >> fDQMData.back()->PedVar_mean[2];
        // 120
        is_stream >> itemp;
        is_stream >> fDQMData.back()->FIR_mean[3];
        is_stream >> fDQMData.back()->FIR_RMS[3];
        is_stream >> fDQMData.back()->LaserRun[3];
        is_stream >> fDQMData.back()->PedVar_mean[3];
        is_stream >> itemp;
        is_stream >> fDQMData.back()->FIR_mean[3];
        is_stream >> fDQMData.back()->FIR_RMS[3];

        for( unsigned int i = 0; i < 4; i++ )
        {
            if( TMath::Abs( fDQMData.back()->FIR_mean[i] ) < 1.e-5 )
            {
                fDQMData.back()->FIR_mean[i] = 0.;
            }
            if( TMath::Abs( fDQMData.back()->FIR_RMS[i] ) < 1.e-5 )
            {
                fDQMData.back()->FIR_RMS[i] = 0.;
            }
        }

        z++;
    }
    is.close();

    cout << "read " << z << " lines from DQM file" << endl;

    return true;
}
