/*! \class VDBRunInfo
    \brief get run info from DB

*/

#include "VDBRunInfo.h"

//
// irun      : runnumber to pull info for
// iDBserver : database url to get runinfo from, usually from VGlobalRunParameter->getDBServer()
// iNTel     : number of telescopes to look for?
VDBRunInfo::VDBRunInfo( int irun, string iDBserver, unsigned int iNTel )
{
    fRunNumber = irun;
    fDBStatus = false;
    
    fTargetName = "";;
    fTargetDec = -99.;
    fTargetRA = -99.;
    fWobbleNorth = 0.;
    fWobbleEast = 0.;
    fConfigMask = 0;
    //fConfigMaskDQM = 0;
    //fConfigMaskNew = 0;
    fTelToAna = 1234;
    fRunType = "";
    fObservingMode = "";
    fRunStatus = "";
    fWeather = "";
    fDataStoppTimeSQL = "";
    fDataStoppTime = 0;
    fDataStartTimeSQL = "";
    fDataStartTime = 0;
    fDataStartTimeHMS = 0;
    fDataStoppTimeHMS = 0;
    fDataStartTimeMJD = 0.;
    fDataStoppTimeMJD = 0.;
    fDuration = 0;
    
    readRunInfoFromDB( iDBserver );
    readRunDQM( iDBserver );
    getLaserRun( iDBserver, fRunNumber, iNTel );
}

void VDBRunInfo::readRunDQM( string iDBserver )
{

    int config_mask_new = readRunDQM( iDBserver, fRunNumber, getConfigMask() );
    fConfigMask = config_mask_new;
    if( fConfigMask == 1 )
    {
        fTelToAna = 1;
    }
    else if( fConfigMask == 2 )
    {
        fTelToAna = 2;
    }
    else if( fConfigMask == 3 )
    {
        fTelToAna = 12;
    }
    else if( fConfigMask == 4 )
    {
        fTelToAna = 3;
    }
    else if( fConfigMask == 5 )
    {
        fTelToAna = 13;
    }
    else if( fConfigMask == 6 )
    {
        fTelToAna = 23;
    }
    else if( fConfigMask == 7 )
    {
        fTelToAna = 123;
    }
    else if( fConfigMask == 8 )
    {
        fTelToAna = 4;
    }
    else if( fConfigMask == 9 )
    {
        fTelToAna = 14;
    }
    else if( fConfigMask == 10 )
    {
        fTelToAna = 24;
    }
    else if( fConfigMask == 11 )
    {
        fTelToAna = 124;
    }
    else if( fConfigMask == 12 )
    {
        fTelToAna = 34;
    }
    else if( fConfigMask == 13 )
    {
        fTelToAna = 134;
    }
    else if( fConfigMask == 14 )
    {
        fTelToAna = 234;
    }
    else if( fConfigMask == 15 )
    {
        fTelToAna = 1234;
    }
    
    return;
    
}

unsigned int VDBRunInfo::readRunDQM( string iDBserver, int run_number , unsigned int config_mask )
{

    unsigned int ConfigMaskDQM = 0;
    unsigned int ConfigMaskNew = 0;
    
    stringstream iTempS;
    iTempS << iDBserver << "/VOFFLINE";
    char c_query[1000];
    sprintf( c_query, "SELECT run_id , data_category   , status   , status_reason , tel_cut_mask , usable_duration , time_cut_mask , light_level , vpm_config_mask , authors  , comment from tblRun_Analysis_Comments where run_id=%d", run_number );
    
    //std::cout<<"VDBRunInfo::readRunDQM "<<std::endl;
    VDB_Connection my_connection( iTempS.str().c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        return config_mask;
    }
    if( !my_connection.make_query( c_query ) )
    {
        return config_mask;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    
    TSQLRow* db_row = db_res->Next();
    if( !db_row )
    {
        cout << "VDBRunInfo:readRunDQM:Info no row in VOFFLINE DB for run " << run_number << endl;
    }
    else
    {
    
        if( db_row->GetField( 4 ) )
        {
            ConfigMaskDQM = ( unsigned int )( atoi( db_row->GetField( 4 ) ) );
        }
        else
        {
            return config_mask;
        }
        
        // Check if the mask is 0
        if( ConfigMaskDQM == 0 )
        {
            return config_mask;
        }
        
        
        
        bitset<4> bitConfig( config_mask );
        bitset<4> bitDQM( ConfigMaskDQM );
        bitset<4> bitNDQM;
        
        for( int i = 0; i < ( int )bitDQM.size(); i++ )
        {
            bitNDQM.set( ( int )bitDQM.size() - i - 1, bitDQM.test( i ) );
        }
        
        bitNDQM = ~bitNDQM;
        bitset<4> bitNewConfig = bitConfig & bitNDQM;
        
        for( int i = 0; i < ( int )bitNewConfig.size(); i++ )
        {
            if( bitNewConfig.test( i ) )
            {
                ConfigMaskNew += ( unsigned int )pow( 2., i );
            }
        }
        
        
        config_mask = ConfigMaskNew;
        
    }
    
    
    
    return config_mask;
    
}

void VDBRunInfo::readRunInfoFromDB( string iDBserver )
{
    stringstream iTempS;
    iTempS << iDBserver << "/VERITAS";
    char c_query[1000];
    sprintf( c_query, "select * from tblRun_Info where run_id=%d", fRunNumber );
    
    //std::cout<<"VDBRunInfo::readRunInfoFromDB "<<std::endl;
    VDB_Connection my_connection( iTempS.str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        fDBStatus = false;
        return;
    }
    if( !my_connection.make_query( c_query ) )
    {
        fDBStatus = false;
        return;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    
    TSQLRow* db_row = db_res->Next();
    if( !db_row )
    {
        cout << "VDBRunInfo: failed reading a row from DB for run " << fRunNumber << endl;
        fDBStatus = false;
        return;
    }
    double iStarttime = 0.;
    double iStopptime = 0.;
    
    // get date
    if( db_row->GetField( 4 ) )
    {
        string iTemp = db_row->GetField( 4 );
        if( iTemp.size() > 8 )
        {
            fDBDate = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
        }
        else
        {
            fDBDate = 0;
        }
    }
    else
    {
        fDBDate = 0;
    }
    if( db_row->GetField( 6 ) )
    {
        string iTemp = db_row->GetField( 6 );
        if( iTemp.size() > 8 )
        {
            fDataStartTimeSQL = iTemp;
            iStarttime = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
            fDataStartTimeHMS = atoi( iTemp.substr( 11, 2 ).c_str() ) * 10000 + atoi( iTemp.substr( 14, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 17, 2 ).c_str() );
            fDataStartTime = atoi( iTemp.substr( 11, 2 ).c_str() ) * 60 * 60 + atoi( iTemp.substr( 14, 2 ).c_str() ) * 60 + atoi( iTemp.substr( 17, 2 ).c_str() );
        }
        else
        {
            iStarttime = 0;
            fDataStartTime = 0;
            fDataStartTimeHMS = 0;
            fDataStartTimeSQL = "";
        }
    }
    else
    {
        iStarttime = 0;
        fDataStartTime = 0;
        fDataStartTimeHMS = 0;
        fDataStartTimeSQL = "";
    }
    if( db_row->GetField( 7 ) )
    {
        string iTemp = db_row->GetField( 7 );
        if( iTemp.size() > 8 )
        {
            fDataStoppTimeSQL = iTemp;
            iStopptime = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
            fDataStoppTimeHMS = atoi( iTemp.substr( 11, 2 ).c_str() ) * 10000 + atoi( iTemp.substr( 14, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 17, 2 ).c_str() );
            fDataStoppTime = atoi( iTemp.substr( 11, 2 ).c_str() ) * 60 * 60 + atoi( iTemp.substr( 14, 2 ).c_str() ) * 60 + atoi( iTemp.substr( 17, 2 ).c_str() );
        }
        else
        {
            iStopptime = 0;
            fDataStoppTime = 0.;
            fDataStoppTimeHMS = 0;
            fDataStoppTimeSQL = "";
        }
    }
    else
    {
        iStopptime = 0;
        fDataStoppTime = 0;
        fDataStoppTimeHMS = 0;
        fDataStoppTimeSQL = "";
    }
    
    if( db_row->GetField( 8 ) )
    {
        string iTemp = db_row->GetField( 8 );
        if( iTemp.size() > 7 )
        {
            fDuration = atoi( iTemp.substr( 0, 1 ).c_str() ) * 3600 + atoi( iTemp.substr( 3, 4 ).c_str() ) * 60 + atoi( iTemp.substr( 6, 7 ).c_str() );
        }
        else
        {
            fDuration = 0;
        }
    }
    else
    {
        fDuration = 0;
    }
    if( TMath::Abs( fDuration < 1.e-4 ) )
    {
        double mjd = 0.;
        double isec_start = 0.;
        double isec_stopp = 0.;
        VSkyCoordinatesUtilities::getMJD_from_SQLstring( fDataStartTimeSQL, mjd, isec_start );
        VSkyCoordinatesUtilities::getMJD_from_SQLstring( fDataStoppTimeSQL, mjd, isec_stopp );
        fDuration = isec_stopp - isec_start;
    }
    
    if( db_row->GetField( 19 ) )
    {
        fTargetName = db_row->GetField( 19 );
    }
    else
    {
        fTargetName = "";
    }
    
    double imjd = 0.;
    int j = 0;
    int iy = iStarttime / 10000;
    int im = ( iStarttime - iy * 10000 ) / 100;
    int id = iStarttime - iy * 10000 - im * 100;
    if( iStarttime > 0 )
    {
        VAstronometry::vlaCldj( iy, im, id, &imjd, &j );
    }
    else
    {
        imjd = 0.;
    }
    fDataStartTimeMJD = imjd;
    iy = iStopptime / 10000;
    im = ( iStopptime - iy * 10000 ) / 100;
    id = iStopptime - iy * 10000 - im * 100;
    if( iStopptime > 0 )
    {
        VAstronometry::vlaCldj( iy, im, id, &imjd, &j );
    }
    else
    {
        imjd = 0.;
    }
    fDataStoppTimeMJD  = imjd;
    // calculate start and stop time
    
    if( db_row->GetField( 1 ) )
    {
        fRunType = db_row->GetField( 1 );
    }
    if( db_row->GetField( 2 ) )
    {
        fObservingMode = db_row->GetField( 2 );
    }
    if( db_row->GetField( 3 ) )
    {
        fRunStatus = db_row->GetField( 3 );
    }
    if( db_row->GetField( 9 ) )
    {
        fWeather = db_row->GetField( 9 );
    }
    
    float dist = 0.;
    if( db_row->GetField( 17 ) )
    {
        dist = atof( db_row->GetField( 17 ) );
    }
    float angl = 0.;
    if( db_row->GetField( 18 ) )
    {
        angl = atof( db_row->GetField( 18 ) );
    }
    fWobbleNorth = dist * cos( angl * TMath::DegToRad() );
    fWobbleEast = dist * sin( angl * TMath::DegToRad() );
    if( TMath::Abs( fWobbleNorth ) < 1.e-15 )
    {
        fWobbleNorth = 0.;
    }
    if( TMath::Abs( fWobbleEast ) < 1.e-15 )
    {
        fWobbleEast = 0.;
    }
    
    // get config mask
    if( db_row->GetField( 10 ) )
    {
        fConfigMask = ( unsigned int )( atoi( db_row->GetField( 10 ) ) );
    }
    else
    {
        fConfigMask = 0;
    }
    if( fConfigMask == 1 )
    {
        fTelToAna = 1;
    }
    else if( fConfigMask == 2 )
    {
        fTelToAna = 2;
    }
    else if( fConfigMask == 3 )
    {
        fTelToAna = 12;
    }
    else if( fConfigMask == 4 )
    {
        fTelToAna = 3;
    }
    else if( fConfigMask == 5 )
    {
        fTelToAna = 13;
    }
    else if( fConfigMask == 6 )
    {
        fTelToAna = 23;
    }
    else if( fConfigMask == 7 )
    {
        fTelToAna = 123;
    }
    else if( fConfigMask == 8 )
    {
        fTelToAna = 4;
    }
    else if( fConfigMask == 9 )
    {
        fTelToAna = 14;
    }
    else if( fConfigMask == 10 )
    {
        fTelToAna = 24;
    }
    else if( fConfigMask == 11 )
    {
        fTelToAna = 124;
    }
    else if( fConfigMask == 12 )
    {
        fTelToAna = 34;
    }
    else if( fConfigMask == 13 )
    {
        fTelToAna = 134;
    }
    else if( fConfigMask == 14 )
    {
        fTelToAna = 234;
    }
    else if( fConfigMask == 15 )
    {
        fTelToAna = 1234;
    }
    
    // get source coordinates
    sprintf( c_query, "select * from tblObserving_Sources where source_id like convert( _utf8 \'%s\' using latin1)", fTargetName.c_str() );
    if( !my_connection.make_query( c_query ) )
    {
        fDBStatus = false;
        return;
    }
    db_res = my_connection.Get_QueryResult();
    
    db_row = db_res->Next();
    if( !db_row )
    {
        cout << " VDBRunInfo::readRunInfoFromDB warning :  no source name in tblObserving_Sources like the name of the target for this run ("
             << fTargetName << ")" << endl;
    }
    else
    {
        fTargetDec = atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi();
        fTargetRA = atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi();
    }
    
    
    fDBStatus = true;
    return;
}


void VDBRunInfo::print()
{
    cout << "Run info from database for run " << fRunNumber << ":" << endl;
    cout << "Date: " << fDBDate << "(" << fDataStartTimeSQL << "," << fDataStoppTimeSQL << ")";
    cout << ", Duration: " << fDuration << " [s]";
    cout << ", " << fRunType << ", " << fObservingMode << ", " << fRunStatus;
    cout << ", Weather: " << fWeather << endl;
    cout << "Target: " << fTargetName;
    cout << " (ra,dec)(J2000)=(" << fTargetRA << ", " << fTargetDec << ")";
    cout << ", Wobble (N,E): " << fWobbleNorth << ", " << fWobbleEast << "), TelToAna: " << fTelToAna << endl;
    cout << "Laser/Flasher runs: ";
    for( unsigned int i = 0; i < getLaserRun().size(); i++ )
    {
        cout << "T" << i + 1 << ": " << getLaserRun()[i] << "   ";
    }
    cout << endl;
}
/*
mysql> describe tblRun_Analysis_Comments;
+-----------------+-----------------------------------------------------------------------------------------------+------+-----+---------+-------+
| Field           | Type                                                                                          | Null | Key | Default | Extra |
+-----------------+-----------------------------------------------------------------------------------------------+------+-----+---------+-------+
| run_id          | int(10) unsigned                                                                              |      | PRI | 0       |       |
| data_category   | enum('science','calibration','engineering','moonfilter','reducedhv','special')                | YES  |     | NULL    |       |
| status          | enum('good_run','needs_adjustments','minor_problems','major_problems','do_not_use','unknown') | YES  |     | NULL    |       |
| status_reason   | enum('weather','external','hardware','multiple','abort')                                      | YES  |     | NULL    |       |
| tel_cut_mask    | tinyint(3) unsigned                                                                           | YES  |     | NULL    |       |
| usable_duration | time                                                                                          | YES  |     | NULL    |       |
| time_cut_mask   | varchar(255)                                                                                  | YES  |     | NULL    |       |
| light_level     | enum('dark','low','moderate','high')                                                          | YES  |     | NULL    |       |
| vpm_config_mask | tinyint(3) unsigned                                                                           | YES  |     | NULL    |       |
| authors         | varchar(255)                                                                                  | YES  |     | NULL    |       |
| comment         | text                                                                                          | YES  |     | NULL    |       |
+-----------------+-----------------------------------------------------------------------------------------------+------+-----+---------+-------+
*/

string VDBRunInfo::readRunDQM_Assessment( string iDBserver, unsigned int iRunNumber )
{
    stringstream iTempS;
    iTempS << iDBserver << "/VOFFLINE";
    char c_query[1000];
    
    sprintf( c_query, "SELECT status from tblRun_Analysis_Comments where run_id=%d", iRunNumber );
    
    VDB_Connection my_connection( iTempS.str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "VDBRunInfo::readRunDQM_Assesment warning: Unable to get DB connection to " << iTempS.str() << endl;
        my_connection.Close_Connection();
        return "unknown";
    }
    if( !my_connection.make_query( c_query ) )
    {
        cout << "VDBRunInfo::readRunDQM_Assesment warning: Unable to make query " << c_query << endl;
        my_connection.Close_Connection();
        return "unknown";
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    string iStatus = "unknown";
    if( db_res->GetRowCount() > 0 )
    {
        TSQLRow* db_row = db_res->Next() ;
        if( !db_row )
        {
            cout << "VDBRunInfo::readRunDQM_Assessment warning: run " << iRunNumber << " not in DB (yet). " << endl;
            my_connection.Close_Connection();
            
        }
        else
        {
            if( !db_row->GetField( 0 ) )
            {
                cout << "VDBRunInfo::readRunDQM_Assessment warning: run " << iRunNumber;
                cout << " does not have a DQM status assigned in the DB." << endl;
                my_connection.Close_Connection();
            }
            else
            {
                iStatus =  db_row->GetField( 0 ) ;
                my_connection.Close_Connection();
            }
        }
        
    }
    return iStatus;
}


vector< unsigned int > VDBRunInfo::getLaserRun( string iDBserver, unsigned int iRunNumber, unsigned int iNTel )
{
    stringstream iTempS;
    iTempS << iDBserver << "/VERITAS";
    char c_query[1000];
    
    sprintf( c_query, "SELECT info.run_id, grp_cmt.excluded_telescopes, info.config_mask FROM tblRun_Info AS info, tblRun_Group AS grp, tblRun_GroupComment AS grp_cmt, (SELECT group_id FROM tblRun_Group WHERE run_id=%d) AS run_grp WHERE grp_cmt.group_id = run_grp.group_id AND grp_cmt.group_type='laser' AND grp_cmt.group_id=grp.group_id AND grp.run_id=info.run_id AND (info.run_type='flasher' OR info.run_type='laser')", iRunNumber );
    
    //std::cout<<"VDBRunInfo::getLaserRun "<<std::endl;
    VDB_Connection my_connection( iTempS.str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        fDBStatus = false;
        return fLaserRunID;
    }
    if( !my_connection.make_query( c_query ) )
    {
        fDBStatus = false;
        return fLaserRunID;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    vector< unsigned int > iLaserList;
    vector< unsigned int > iLaserConfigMask;
    vector< unsigned int > iLaserExclude;
    
    if( db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            iLaserList.push_back( atoi( db_row->GetField( 0 ) ) );
            iLaserExclude.push_back( atoi( db_row->GetField( 1 ) ) );
            iLaserConfigMask.push_back( atoi( db_row->GetField( 2 ) ) );
        }
        
    }
    else
    {
        cout << "WARNING: VDBRunInfo::getLaserRun() no laser run found for telescope " << iNTel << " and run " << iRunNumber << endl;
    }
    
    fLaserRunID.assign( iNTel, 0 );
    for( unsigned int t = 0; t < iNTel; t++ )
    {
        // check if this run is excluded from group
        // also check if telescope is within the config mask (taking DQM cuts into account)
        for( unsigned int i = 0; i < iLaserList.size(); i++ )
        {
            bitset< 8 > ibit( iLaserExclude[i] );
            unsigned int config_mask = readRunDQM( iDBserver, iLaserList[i], iLaserConfigMask[i] );
            bitset< 8 > ibit_mask( config_mask );
            if( !ibit.test( t ) && ibit_mask.test( t ) )
            {
            
                // make sure this run is not marked as 'do not use' by DQM
                string iStatus = readRunDQM_Assessment( iDBserver, iLaserList[i] );
                cout << "VDBRunInfo::getLaserRun(), T" << t + 1 << ": laser run " << iLaserList[i] << " has status " << iStatus << endl;
                
                if( iStatus == "do_not_use" )
                {
                    cout << "VDBRunInfo::getLaserRun() Warning: laser run " << iLaserList[i] << " was assessed as \"do not use\" by DQM.";
                    cout << " Issue needs to be investigated, use analysis results might be affected. " << endl;
                }
                fLaserRunID[t] = iLaserList[i];
            }
        }
    }
    
    return fLaserRunID;
}
