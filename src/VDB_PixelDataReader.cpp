/* \file  VDB_PixelDataReader
 * \brief reading pixel data from DB
 *
 */

#include "VDB_PixelDataReader.h"

VDB_PixelDataReader::VDB_PixelDataReader( vector< unsigned int > nPixel_per_telescope )
{
    setDebug( true );
    fDBStatus = true;
    
    fNPixel = nPixel_per_telescope;
    
    //////////////////////////////////////////////////////////////////
    // coding of different data vectors:
    //  ID = 0: L1 rates
    //  ID = 1: HV
    //  ID = 3: currents
    //  note: this is propagated into the data vector fPixelData
    fPixelDataType.push_back( "L1_Rates" );
    fPixelDataType.push_back( "HV" );
    fPixelDataType.push_back( "Currents" );
    fPixelDataType.push_back( "FADC_board" );
    fPixelDataType.push_back( "FADC_channel" );
    
    // define data vectors and histograms
    char htitle[800];
    char hname[800];
    for( unsigned int t = 0; t < fPixelDataType.size(); t++ )
    {
        vector< vector< VDB_PixelData* > > iD;
        vector< TH1F* > i_his;
        for( unsigned int i = 0; i < nPixel_per_telescope.size(); i++ )
        {
            // histograms
            sprintf( hname, "h_%s_tel%d", fPixelDataType[t].c_str(), i + 1 );
            sprintf( htitle, "%s (tel %d)", fPixelDataType[t].c_str(), i + 1 );
            i_his.push_back( new TH1F( hname, htitle, 200, 0., 1. ) );
            i_his.back()->SetXTitle( fPixelDataType[t].c_str() );
            // data
            vector< VDB_PixelData* > iP;
            for( unsigned int j = 0; j < nPixel_per_telescope[i]; j++ )
            {
                iP.push_back( new VDB_PixelData( fPixelDataType[t] ) );
                // preliminary: set sampling interval to 60s
                iP.back()->fTimeBinWidth_s = 60.;
            }
            iD.push_back( iP );
        }
        fPixelData.push_back( iD );
        fPixelData_histogram.push_back( i_his );
    }
    //////////////////////////////////////////////////////////////////
}

void VDB_PixelDataReader::print( string iDataType, unsigned int iTelID, unsigned int iPixel )
{
    for( unsigned int i = 0; i < fPixelDataType.size(); i++ )
    {
        if( fPixelDataType[i] == iDataType )
        {
            if( iTelID < fPixelData[i].size() && iPixel < fPixelData[i][iTelID].size() )
            {
                fPixelData[i][iTelID][iPixel]->print();
            }
            else
            {
                cout << "VDB_PixelDataReader::print() error: invalid telescope ID/pixel ID (";
                cout << fPixelDataType[i] << ", " << iTelID << "," << iPixel << ")" << endl;
                return;
            }
            return;
        }
    }
}

void VDB_PixelDataReader::print()
{
    for( unsigned int i = 0; i < fPixelDataType.size(); i++ )
    {
        cout << "VDB_PixelDataReader database: " << fPixelDataType[i] << ":";
        if( fPixelDataType[i].size() < 10 )
        {
            cout << "\t\t";
        }
        else
        {
            cout << "\t";
        }
        cout << "# tel=" << fPixelData[i].size();
        if( fPixelData[i].size() > 0 )
        {
            cout << ", # pixels: " << fPixelData[i][0].size();
            cout << ", # time stamps: " ;
            for( unsigned int iTel = 0; iTel < fPixelData[i].size(); iTel++ )
            {
                if( fPixelData[i][iTel].size() > 0 )
                {
                    cout << fPixelData[i][iTel][0]->fMJD.size() << " " ;
                }
            }
            cout << endl;
        }
    }
}

void VDB_PixelDataReader::fillDataRow( unsigned int iDataType, string iTimeStamp, int iTel, int iPixel, float iData )
{
    if( iDataType < fPixelData.size() )
    {
        if( ( unsigned int )iTel < fPixelData[iDataType].size() && ( unsigned int )iPixel < fPixelData[iDataType][iTel].size() )
        {
            double i_mjd = 0;
            double i_sec_of_day = 0;
            VSkyCoordinatesUtilities::getMJD_from_SQLstring( iTimeStamp, i_mjd, i_sec_of_day );
            fPixelData[iDataType][iTel][iPixel]->fMJD.push_back( i_mjd );
            fPixelData[iDataType][iTel][iPixel]->fsec_of_day.push_back( i_sec_of_day );
            fPixelData[iDataType][iTel][iPixel]->fData.push_back( iData );
            if( fDebug )
            {
                cout << "VDB_PixelDataReader::fillDataRow fill data type " << iDataType << ", time ";
                cout << iTimeStamp << ", T" << iTel + 1 << ", pixel " << iPixel << ": " << iData << endl;
            }
        }
    }
    else
    {
        cout << "VDB_PixelDataReader::fillDataRow(): invalid data type: " << iDataType << endl;
    }
}

/*
 * access the VERITAS DB and retrieve all values
 * (this is the only function where the DB is queried)
 *
 */
bool VDB_PixelDataReader::readFromDB( string iDBserver, unsigned int runNumber, string iDBStartTimeSQL, string fDBRunStoppTimeSQL )
{
    if( fDebug )
    {
        cout << "VDB_PixelDataReader::readFromDB: " << runNumber << endl;
    }
    
    stringstream iTempS;
    iTempS << iDBserver << "/VERITAS";
    char c_query[1000];
    VDB_Connection my_connection( iTempS.str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        fDBStatus = false;
        return false;
    }
    
    /////////////////////////////////////////////////////
    // read L1 rates
    
    ostringstream c_queryS;
    c_queryS << "select timestamp, telescope_id, pixel_id, rate from tblL1_TriggerInfo, tblRun_Info where timestamp";
    c_queryS << " >= tblRun_Info.data_start_time - INTERVAL 1 MINUTE AND timestamp ";
    c_queryS << " <=  tblRun_Info.data_end_time + INTERVAL 1 MINUTE AND tblRun_Info.run_id=";
    c_queryS << runNumber;
    
    if( !my_connection.make_query( c_queryS.str().c_str() ) )
    {
        fDBStatus = false;
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    if( db_res && db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( !db_row )
            {
                cout << "VDB_PixelDataReader::readFromDB: failed reading a DB row" << endl;
                fDBStatus = false;
                return false;
            }
            // mysql> describe tblL1_TriggerInfo;
            // +--------------+----------------------+------+-----+---------------------+-------+
            // | Field        | Type                 | Null | Key | Default             | Extra |
            // +--------------+----------------------+------+-----+---------------------+-------+
            // | timestamp    | datetime             |      | PRI | 0000-00-00 00:00:00 |       |
            // | run_id       | int(11)              |      | PRI | 0                   |       |
            // | telescope_id | tinyint(3) unsigned  |      | PRI | 0                   |       |
            // | pixel_id     | smallint(5) unsigned |      | PRI | 0                   |       |
            // | rate         | float                | YES  |     | NULL                |       |
            // +--------------+----------------------+------+-----+---------------------+-------+
            if( db_row->GetField( 0 ) &&  db_row->GetField( 1 ) && db_row->GetField( 2 ) && db_row->GetField( 3 ) )
            {
                fillDataRow( 0, db_row->GetField( 0 ), atoi( db_row->GetField( 1 ) ), atoi( db_row->GetField( 2 ) ), atof( db_row->GetField( 3 ) ) );
            }
        }
    }
    
    /////////////////////////////////////////////////////
    // read HV and currents measured
    
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        ostringstream c_queryS;
        c_queryS << "select * FROM tblHV_Telescope" << i << "_Status WHERE channel > 0 AND";
        c_queryS << " (db_start_time >= \"" << iDBStartTimeSQL << "\"";
        c_queryS << "- INTERVAL 1 MINUTE) AND (db_start_time <= \" " << fDBRunStoppTimeSQL << "\" )";
        
        if( !my_connection.make_query( c_queryS.str().c_str() ) )
        {
            fDBStatus = false;
            cout << "FAILED" << endl;
            return false;
        }
        TSQLResult* db_res = my_connection.Get_QueryResult();
        
        if( db_res && db_res->GetRowCount() > 0 )
        {
            while( TSQLRow* db_row = db_res->Next() )
            {
                if( !db_row )
                {
                    cout << "VDB_PixelDataReader::readFromDB: failed reading a DB row" << endl;
                    fDBStatus = false;
                    return false;
                }
                // +---------------+----------------------+------+-----+---------------------+-------+
                // | Field         | Type                 | Null | Key | Default             | Extra |
                // +---------------+----------------------+------+-----+---------------------+-------+
                // | db_start_time | datetime             |      | PRI | 0000-00-00 00:00:00 |       |
                // | db_end_time   | datetime             | YES  | MUL | NULL                |       |
                // | channel       | smallint(5) unsigned |      | PRI | 0                   |       |
                // | voltage_meas  | float                | YES  |     | NULL                |       |
                // | current_meas  | float                | YES  |     | NULL                |       |
                // +---------------+----------------------+------+-----+---------------------+-------+
                if( db_row->GetField( 1 ) && db_row->GetField( 2 ) && db_row->GetField( 3 ) )
                {
                    fillDataRow( 1, db_row->GetField( 1 ), i, atoi( db_row->GetField( 2 ) ) - 1, atof( db_row->GetField( 3 ) ) );
                }
                if( db_row->GetField( 1 ) && db_row->GetField( 2 ) && db_row->GetField( 4 ) )
                {
                    fillDataRow( 2, db_row->GetField( 1 ), i, atoi( db_row->GetField( 2 ) ) - 1, atof( db_row->GetField( 4 ) ) );
                }
            }
        }
    }
    
    
    //read FADC settings
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        //GetField( 0 ) = pixel_id (really channel, starts counting at 0)
        //GetField( 1 ) = fadc_id i.e. FADC module number
        //GetField( 2 ) = fadc_channel (per module)
        
        
        sprintf( c_query, "select c.pixel_id , s.fadc_id, c.fadc_channel from tblFADC_Slot_Relation as s, tblFADC_Channel_Relation as c where s.db_start_time < \"%s\" and c.db_start_time < \"%s\" and ( s.db_end_time IS NULL or s.db_end_time > \"%s\" ) and ( c.db_end_time IS NULL or c.db_end_time > \"%s\" ) and s.fadc_crate=c.fadc_crate and s.fadc_slot=c.fadc_slot and s.telescope_id=c.telescope_id and c.pixel_id IS NOT NULL and s.telescope_id=%d order by c.pixel_id ;", iDBStartTimeSQL.c_str(), iDBStartTimeSQL.c_str(), fDBRunStoppTimeSQL.c_str(), fDBRunStoppTimeSQL.c_str(), i );
        
        if( !my_connection.make_query( c_query ) )
        {
            fDBStatus = false;
            cout << "VDB_PixelDataReader::readFromDB: FAILED making query for reading FADC status" << endl;
            return false;
        }
        
        TSQLResult* db_res = my_connection.Get_QueryResult();
        
        if( db_res && db_res->GetRowCount() > 0 )
        {
            while( TSQLRow* db_row = db_res->Next() )
            {
                if( !db_row )
                {
                    cout << "VDB_PixelDataReader::readFromDB: failed reading a DB row" << endl;
                    fDBStatus = false;
                    return false;
                }
                
                if( db_row->GetField( 0 ) && db_row->GetField( 1 ) && db_row->GetField( 2 ) )
                {
                    fillDataRow( 3, iDBStartTimeSQL, i, atoi( db_row->GetField( 0 ) ) , atof( db_row->GetField( 1 ) ) );
                    fillDataRow( 4, iDBStartTimeSQL, i, atoi( db_row->GetField( 0 ) ) , atof( db_row->GetField( 2 ) ) );
                }
            }
        }
    }
    print();
    
    return true;
}

bool VDB_PixelDataReader::writeDataTree( unsigned int iTel )
{

    char htitle[200];
    char hname[200];
    float mjd = 0.;
    float time = 0.;
    float value[VDST_MAXCHANNELS];
    unsigned int nPixel = 0;
    for( unsigned int i = 0; i < fPixelData.size(); i++ )
    {
        if( iTel >= fPixelData[i].size() )
        {
            cout << "VDB_PixelDataReader::writeDataTree error: invalid data telescope number: ";
            cout << iTel << ", allowed is up to " << fPixelData[i].size() << endl;
            return false;
        }
        sprintf( htitle, "DB pixel data for %s (Tel %d)", fPixelDataType[i].c_str(), iTel + 1 );
        sprintf( hname, "dbpixeldata_%s", fPixelDataType[i].c_str() );
        TTree iTree( hname, htitle );
        iTree.Branch( "mjd", &mjd, "mjd/F" );
        iTree.Branch( "sec_of_day", &time, "sec_of_day/F" );
        iTree.Branch( "npixel", &nPixel, "npixel/i" );
        sprintf( hname, "%s", fPixelDataType[i].c_str() );
        sprintf( htitle, "%s[npixel]", fPixelDataType[i].c_str() );
        iTree.Branch( hname, value, htitle );
        
        if( fPixelData[i][iTel].size() > 0 )
        {
            nPixel = fPixelData[i][iTel].size();
            // loop over all MJDs (assuming the same for all pixel)
            for( unsigned int t = 0; t < fPixelData[i][iTel][0]->fMJD.size(); t++ )
            {
                mjd  = fPixelData[i][iTel][0]->fMJD[t];
                time = fPixelData[i][iTel][0]->fsec_of_day[t];
                for( unsigned int p = 0; p < fPixelData[i][iTel].size(); p++ )
                {
                    if( t < fPixelData[i][iTel][p]->fMJD.size() )
                    {
                        value[p] = fPixelData[i][iTel][p]->fData[t];
                    }
                    else
                    {
                        value[p] = -99.;
                    }
                }
                iTree.Fill();
            }
        }
        
        iTree.Write();
    }
    return true;
}

/*
 * return distribution of data values
 *
 */
TH1F* VDB_PixelDataReader::getDataHistogram( unsigned int iDataType, unsigned int iTel, int iMJD, float iTime )
{
    // fills fDummyReturnVector with data
    getDataVector( iDataType, iTel, iMJD, iTime );
    if( fDummyReturnVector.size() == 0 )
    {
        return 0;
    }
    
    // get min/maximum
    float v_min = *std::min_element( fDummyReturnVector.begin(), fDummyReturnVector.end() );
    float v_max = *std::max_element( fDummyReturnVector.begin(), fDummyReturnVector.end() );
    // make sure that min/max give a valid histogram
    if( v_min < 0. )
    {
        v_min = 0.;
    }
    if( TMath::Abs( v_max - v_min ) < 1.e-2 )
    {
        v_max = v_min + 1.;
    }
    if( fPixelData_histogram[iDataType][iTel] )
    {
        fPixelData_histogram[iDataType][iTel]->Reset();
        fPixelData_histogram[iDataType][iTel]->SetBins( fPixelData_histogram[iDataType][iTel]->GetNbinsX(), 0., v_max );
        for( unsigned int i = 0; i < fDummyReturnVector.size(); i++ )
        {
            fPixelData_histogram[iDataType][iTel]->Fill( fDummyReturnVector[i] );
        }
        return fPixelData_histogram[iDataType][iTel];
    }
    
    return 0;
}

/*
 * get a single value for a given telescope, channel, and time
 *
 * this function is very inefficient, don't call it in a loop!
 */
float VDB_PixelDataReader::getValue( unsigned int iDataType, unsigned int iTel, unsigned int iChannel, int iMJD, float iTime )
{
    // fills fDummyReturnVector with data
    getDataVector( iDataType, iTel, iMJD, iTime );
    if( fDummyReturnVector.size() == 0 )
    {
        return 0;
    }
    
    if( iChannel < fDummyReturnVector.size() )
    {
        return fDummyReturnVector[iChannel];
    }
    
    return 0.;
}


/*
 * return a data vector of size npixel
 *
 * different 'incomplete' returns:
 *
 * case 1: iDataType is an invalid data type:  zero size vector is returned
 * case 2: iTel is an invalid telescope number:  zero size vector is returned
 * case 3: iMJD is insided any known time interval: vector of size npixel filled with zeros is returned
 */
vector< float > VDB_PixelDataReader::getDataVector( unsigned int iDataType, unsigned int iTel, int iMJD, float iTime )
{
    fDummyReturnVector.clear();
    // is this a good data type ID?
    if( iDataType < fPixelData.size() )
    {
        // check telescope number
        if( iTel < fPixelData[iDataType].size() )
        {
            // reset return vector
            fDummyReturnVector.assign( fPixelData[iDataType][iTel].size(), 0. );
            // loop over all pixel
            for( unsigned int i = 0; i < fPixelData[iDataType][iTel].size(); i++ )
            {
                // loop over all pixel times to get correct time
                if( fPixelData[iDataType][iTel][i]->fMJD.size() > 1 )
                {
                    for( unsigned int t = 0; t < fPixelData[iDataType][iTel][i]->fMJD.size(); t++ )
                    {
                        // this should never happen
                        if( TMath::Abs( fPixelData[iDataType][iTel][i]->fMJD[t] - iMJD ) > 1.e-1 )
                        {
                            continue;
                        }
                        // first measurement is sometimes a bit late: look in a wide search window
                        double iTimeBinWidth = fPixelData[iDataType][iTel][i]->fTimeBinWidth_s;
                        if( t == 0 )
                        {
                            iTimeBinWidth *= 4.;
                        }
                        // check time and fill return vector
                        if( iTime >= fPixelData[iDataType][iTel][i]->fsec_of_day[t] - iTimeBinWidth
                                && iTime < fPixelData[iDataType][iTel][i]->fsec_of_day[t] )
                        {
                            fDummyReturnVector[i] = fPixelData[iDataType][iTel][i]->fData[t];
                            // first two and/or last time bin is sometimes not filled (run start/end and L1 rate interval star mismatch)
                            // use in this case the L1 rates from second bin
                            // (note: only for L1 rates)
                            if( iDataType == 0 && fDummyReturnVector[i] < 1.e-2 && fPixelData[iDataType][iTel][i]->fData.size() > 1 )
                            {
                                if( t == 0 || t == 1 )
                                {
                                    if( fPixelData[iDataType][iTel][i]->fData.size() > t + 1 )
                                    {
                                        fDummyReturnVector[i] = fPixelData[iDataType][iTel][i]->fData[t + 1];
                                    }
                                    if( fDummyReturnVector[i] < 1.e-2 && fPixelData[iDataType][iTel][i]->fData.size() > 3 )
                                    {
                                        fDummyReturnVector[i] = fPixelData[iDataType][iTel][i]->fData[t + 2];
                                    }
                                }
                                else if( t == fPixelData[iDataType][iTel][i]->fMJD.size() - 1 )
                                {
                                    fDummyReturnVector[i] = fPixelData[iDataType][iTel][i]->fData[t - 1];
                                }
                            }
                            break;
                        }
                    }
                }
                else if( fPixelData[iDataType][iTel][i]->fMJD.size() == 1 )
                {
                    fDummyReturnVector[i] = fPixelData[iDataType][iTel][i]->fData[0];
                    
                }
            }
        }
    }
    return fDummyReturnVector;
}

/*
 * return list of channels with data outside of the given range
 *
 */
vector< unsigned int > VDB_PixelDataReader::getDeadChannelList( unsigned int iDataType, unsigned int iTel, int iMJD,
        float iTime, float i_min, float i_max, bool bRMS )
{
    vector< unsigned int > i_channelList;
    // fills fDummyReturnVector with data
    getDataVector( iDataType, iTel, iMJD, iTime );
    if( fDummyReturnVector.size() == 0 )
    {
        return i_channelList;
    }
    
    // for relative differences, calculate mean and rms
    double i_mean = 0.;
    double i_rms = 0.;
    if( bRMS )
    {
        double i_mean2 = 0.;
        double i_n = 0.;
        for( unsigned int i = 0; i < fDummyReturnVector.size(); i++ )
        {
            if( fDummyReturnVector[i] > 0. )
            {
                i_mean  += fDummyReturnVector[i];
                i_mean2 += fDummyReturnVector[i] * fDummyReturnVector[i];
                i_n++;
            }
        }
        if( i_n > 1 )
        {
            i_rms = sqrt( ( 1. / ( i_n - 1. ) ) * ( i_mean2 - i_mean * i_mean / i_n ) );
            i_mean /= i_n;
        }
        
        i_min = i_mean - TMath::Abs( i_min ) * i_rms;
        i_max = i_mean + TMath::Abs( i_max ) * i_rms;
    }
    
    for( unsigned int i = 0; i < fDummyReturnVector.size(); i++ )
    {
        if( fDummyReturnVector[i] < i_min )
        {
            i_channelList.push_back( i );
        }
        else if( fDummyReturnVector[i] > i_max )
        {
            i_channelList.push_back( i );
        }
    }
    
    return i_channelList;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

VDB_PixelData::VDB_PixelData( string iDataType )
{
    fDataType = iDataType;
    fTimeBinWidth_s = 60.;
}

void VDB_PixelData::print()
{
    for( unsigned int i = 0; i < fMJD.size(); i++ )
    {
        cout << "\t Pixel: " << i + 1;
        cout << ", MJD " << fMJD[i];
        cout << ", sec of day " << fsec_of_day[i];
        cout << ", " << fDataType << ": " << fData[i];
        cout << endl;
    }
}
