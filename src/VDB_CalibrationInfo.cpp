/*! \class VDB_CalibrationInfo
    \brief read or get laser run calibration info from VOFFLINE DB

*/

#include "VDB_CalibrationInfo.h"
//======================== PUBLIC ================

/******************** CONSTRUCTORS ***************/
VDB_CalibrationInfo::VDB_CalibrationInfo()
{
    //..........................
    // TO BE COMPLETED maybe
    //..........................
    
    fread_flag = false;
    fcurrent_run = -111;
    fcurrent_tel = -111 ;
    fFile_to_write = "le_nom";
    freading_gain_or_toff = -111;
    fVOFFLINE_version_query = -111;
    fLOW_GAIN = -111;
    fServer  = "";
    
    fwrite_flag = false;
    fdofile_flag = false;
    fVOFFLINE_version_query = -1000;
    fquery_read = "";
    fNOW_DB = "demain";
    
    flong_char_query = 100000;
    fseparation = "|";
    
}
VDB_CalibrationInfo::VDB_CalibrationInfo( TString file_to_write_in_DB, TString DBserver, int low_gain )
{

    fwrite_flag = true;
    fServer  = DBserver;
    fFile_to_write = file_to_write_in_DB;
    fLOW_GAIN = low_gain;
    
    //--- global parameter not specified in the call
    fcurrent_run = -1000;
    fcurrent_tel = -1000 ;
    fNOW_DB = "demain";
    fVOFFLINE_version_query = -1000;
    fdofile_flag = false;
    fread_flag = false;
    fquery_read = "";
    
    flong_char_query = 100000;
    fseparation = "|";
}
VDB_CalibrationInfo::VDB_CalibrationInfo( int current_run, int current_tel, TString NOW_DB, int VOFFLINE_version_query, int LOW_GAIN )
{
    fdofile_flag = true;
    fcurrent_run = current_run;
    fcurrent_tel = current_tel ;
    fNOW_DB = NOW_DB;
    fVOFFLINE_version_query = VOFFLINE_version_query;
    fLOW_GAIN = LOW_GAIN;
    
    //--- global parameter not specified in the call
    fServer  = "";
    fFile_to_write = "";
    fwrite_flag = false;
    fread_flag = false;
    fquery_read = "";
    
    flong_char_query = 100000;
    fseparation = "|";
}
VDB_CalibrationInfo::VDB_CalibrationInfo( int laserrun , int tel , string name_out_file, int gain_or_toff, int VOFFLINE_version_query, int LOW_GAIN, TString DBserver )
{
    fread_flag = true;
    fcurrent_run = laserrun;
    fcurrent_tel = tel ;
    fFile_to_write = name_out_file;
    freading_gain_or_toff = gain_or_toff;
    fVOFFLINE_version_query = VOFFLINE_version_query;
    fLOW_GAIN = LOW_GAIN;
    fServer  = DBserver;
    
    //--- global parameter not specified in the call
    fwrite_flag = false;
    fdofile_flag = false;
    fVOFFLINE_version_query = -1000;
    fquery_read = "";
    fNOW_DB = "demain";
    flong_char_query = 100000;
    fseparation = "|";
    
}

/******************************* Writing function *******************************/
//---------------------------------------------------------------------------------------------------------------------------------------
//-- VDB_CalibrationInfo::DoFile_for_DBwriting                                                                                                         --
//-- write a file with the correct format in order to use write in VOFFLINE DB table: tblEventDisplay_Analysis_Calibration_Flasher     --
//-- public
//---------------------------------------------------------------------------------------------------------------------------------------
bool VDB_CalibrationInfo::DoFile_for_DBwriting( vector < double > Vchannel_gain, vector < double > Vmean_gain, vector < double > Vvar_gain, vector < double > Vchannel_toff, vector < double > Vmean_toff, vector < double > Vvar_toff, FILE*& wDB_file )
{

    bool writing_ok = false;
    if( !fdofile_flag )
    {
        std::cout << "ERROR: wrong constructor" << std::endl;
        std::cout << "use: VDB_CalibrationInfo::VDB_CalibrationInfo(int current_run,int current_tel,TString NOW_DB,int VOFFLINE_version_query,int LOW_GAIN)" << std::endl;
        std::cout << "if you want to call VDB_CalibrationInfo::DoFile_for_DBwriting " << std::endl;
        return writing_ok;
    }
    if( Vchannel_gain.size() != Vchannel_toff.size() )
    {
        printf( " ERROR: RUN %d tel %d \n number of channel in gain (%d) different from the number of channel in toff (%d) \n", fcurrent_run, fcurrent_tel, ( int ) Vchannel_gain.size(), ( int ) Vchannel_toff.size() );
        return writing_ok;
    }
    for( unsigned int i = 0; i < Vchannel_gain.size(); i++ )
    {
    
        if( Vchannel_gain[i] != Vchannel_toff[i] )
        {
            printf( " WARNING: RUN %d tel %d \n channel in gain (%f) different from channel in toff (%f) \n", fcurrent_run, fcurrent_tel, Vchannel_gain[i], Vchannel_toff[i] );
            return writing_ok;
        }
        else
        {
            fprintf( wDB_file, "%5d %s %d %s %3d %s %s %s %.4f %s %.4f %s %+.4f %s %+.4f %s %5d %s %5d \n", fcurrent_run, fseparation.Data(), fcurrent_tel, fseparation.Data(), ( int ) Vchannel_gain[i], fseparation.Data(), fNOW_DB.Data(), fseparation.Data(), Vmean_gain[i], fseparation.Data(), Vvar_gain[i], fseparation.Data(), Vmean_toff[i], fseparation.Data(), Vvar_toff[i], fseparation.Data(), fVOFFLINE_version_query, fseparation.Data(), fLOW_GAIN );
            
            
            writing_ok = true;
        }
        
    }
    
    return writing_ok;
}

/*

mysql> show columns from tblEventDisplay_Analysis_Calibration_Flasher;
+--------------------+-------------+------+-----+---------------------+-------+
| Field              | Type        | Null | Key | Default             | Extra |
+--------------------+-------------+------+-----+---------------------+-------+
| run_id             | int(11)     |      | PRI | 0                   |       |
| telescope          | tinyint(1)  |      | PRI | 0                   |       |
| channel_id         | smallint(5) |      | PRI | 0                   |       |
| update_time        | datetime    |      | PRI | 0000-00-00 00:00:00 |       |
| gain_mean          | float       |      |     | 0                   |       |
| gain_var           | float       |      |     | 0                   |       |
| toffset_mean       | float       |      |     | 0                   |       |
| toffset_var        | float       |      |     | 0                   |       |
| code_version       | smallint(5) |      | PRI | 0                   |       |
| high_low_gain_flag | smallint(5) |      |     | 0                   |       |
+--------------------+-------------+------+-----+---------------------+-------+

mysql> show columns from tblEventDisplay_Analysis_Status;
+-----------+------------------+------+-----+---------+-------+
| Field     | Type             | Null | Key | Default | Extra |
+-----------+------------------+------+-----+---------+-------+
| run_id    | int(10) unsigned | YES  | MUL | NULL    |       |
| version   | varchar(10)      | YES  |     | NULL    |       |
| laser_run | int(10) unsigned | YES  |     | NULL    |       |
| status    | int(10) unsigned | YES  |     | NULL    |       |
+-----------+------------------+------+-----+---------+-------+

*/

//--------------------------------------------------------------
//-- VDB_CalibrationInfo::write_inVOFFLINE_DB_from_file()                     --
//-- access VOFFLINE_DB for writing with an input file        --
//-- public
//--------------------------------------------------------------
void VDB_CalibrationInfo::write_inVOFFLINE_DB_from_file( string pass_word )
{

    if( !fwrite_flag )
    {
        std::cout << "ERROR: wrong constructor" << std::endl;
        std::cout << "use: VDB_CalibrationInfo::VDB_CalibrationInfo(TString file_to_write_in_DB,TString DBserver) " << std::endl;
        std::cout << "if you want to call  VDB_CalibrationInfo::write_inVOFFLINE_DB_from_file()" << std::endl;
        return;
    }
    
    //---- open the DB and check
    string iTempS;
    iTempS =  fServer + "/VOFFLINE?local";
    
    std::cout << "You are in the process of WRITING in VOFFLINE Data Base " << std::endl;
    std::cout << "in order to complete the process: " << std::endl;
    string the_password = pass_word;
    if( the_password == "" )
    {
        please_give_the_password();
    }
    
    //std::cout<<"VDB_CalibrationInfo::write_inVOFFLINE_DB_from_file "<<std::endl;
    VDB_Connection my_connection( iTempS.c_str(), "readwrite", the_password.c_str() ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        return;
    }
    string query = WriteQuery_to_write_in_DB();
    std::cout << "query " << query << std::endl;
    
    //---- do the query and check
    if( !my_connection.make_query( query.c_str() ) )
    {
        cout << "WARNING VDB_CalibrationInfo::write_inVOFFLINE_DB_from_file: failed to get something from the query " << endl;
        return;
    }
    
    return;
}

//---------------------------------------------
//-- VDB_CalibrationInfo::readVOFFLINE
//-- reading the VOFFLINE DB table
//-- put the result in a FILE
//-- public
//---------------------------------------------
bool VDB_CalibrationInfo::readVOFFLINE()
{
    std::cout << "Reading calibration information from VOFFLINE DB" << std::endl;
    
    if( !fread_flag )
    {
        std::cout << "ERROR: wrong constructor" << std::endl;
        std::cout << "use: VDB_CalibrationInfo::VDB_CalibrationInfo(int laserrun ,int tel ,string name_out_file,int gain_or_toff,int VOFFLINE_version_query,TString DBserver)" << std::endl;
        std::cout << "if you want to call VDB_CalibrationInfo::readVOFFLINE()" << std::endl;
        Vchannel.clear();
        Vmean.clear();
        Vvar.clear();
        return false;
    }
    
    //-- Create the query to read the DB
    Create_query_read();
    
    //-- read the DB and put the result in the vector
    if( Read_the_DB() )
    {
    
        //-- write the result of the reading in a FILE
        if( Vchannel.size() > 0 && fFile_to_write.Length() > 0 )
        {
        
            // make sure that directory exists
            gSystem->mkdir( gSystem->DirName( fFile_to_write.Data() ), true );
            FILE* file_output;
            file_output = fopen( fFile_to_write.Data(), "w" );
            for( unsigned int i = 0 ; i < Vchannel.size(); i++ )
            {
                fprintf( file_output, "%3d %.4f %.4f  \n", ( int ) Vchannel[i], Vmean[i], Vvar[i] );
            }
            fclose( file_output );
            return true;
        }
        else if( Vchannel.size() == 0 )
        {
            std::cout << "ERROR  VDB_CalibrationInfo::readVOFFLINE Vchannel.size()<1" << std::endl;
            
            return false;
        }
        return true;
        
    }
    else
    {
    
        std::cout << "ERROR reading went wrong " << std::endl;
        return false;
        
    }
    
}




//======================== PRIVATE ================


/**************** Writing in VOFFLINE DB functions ***********************/
//------------------------------------------------------------------------
//-- VDB_CalibrationInfo::please_give_the_password
//-- Interactive password asking (access the DB with readwrite rights)  --
//-- private
//------------------------------------------------------------------------
string VDB_CalibrationInfo::please_give_the_password()
{

    string password = "";
    char str [80];
    
    
    printf( "Please enter the password: " );
    int result_scan = scanf( "%s", str );
    
    if( result_scan > 0 )
    {
        password = str;
    }
    
    return password;
    
}
//------------------------------------------------------------------------
//-- VDB_CalibrationInfo::WriteQuery_to_write_in_DB()
//-- Write the query to copy in DB
//-- private
//------------------------------------------------------------------------
string VDB_CalibrationInfo::WriteQuery_to_write_in_DB()
{

    string query = "";
    
    char c_query[flong_char_query];// has to be long if we ask for a long run list in VDBSourceInfo
    sprintf( c_query, "LOAD DATA LOCAL INFILE '%s' INTO TABLE tblEventDisplay_Analysis_Calibration_Flasher FIELDS TERMINATED BY '%s' ;", fFile_to_write.Data(), fseparation.Data() );
    
    query = c_query;
    return query;
    
}

/******************** Reading VOFFLINE DB functions *******************************/

//---------------------------------------------------------------------------
//-- VDB_CalibrationInfo::Read_the_DB
//-- go read the DB with fquery_read
//-- fill the given vector with the result of the reading
//-- private
//
//---------------------------------------------------------------------------
bool VDB_CalibrationInfo::Read_the_DB()
{
    Vchannel.clear();
    Vmean.clear();
    Vvar.clear();
    
    //---- open the DB and check
    stringstream iTempS;
    iTempS << fServer << "/VOFFLINE";
    
    //std::cout<<"VDB_CalibrationInfo::Read_the_DB "<<std::endl;
    VDB_Connection my_connection( iTempS.str().c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "ERROR VDB_CalibrationInfo::Read_the_DB(): failed to connect to database server" << endl;
        cout << "\t server: " <<  fServer << endl;
        return false;
    }
    //---- do the query and check
    if( !my_connection.make_query( fquery_read.c_str() ) )
    {
        cout << "ERROR VDB_CalibrationInfo::Read_the_DB(): failed to get something from the query " << endl;
        cout << "ERROR laser run  " << fcurrent_run << " tel " << fcurrent_tel << " is not in the VOFFLINE DB (yet?)" << std::endl;
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    //---- read the query
    if( db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( !db_row )
            {
                cout << "WARNING VDB_CalibrationInfo::Read_the_DB(): failed reading a row from DB " << endl;
                cout << "ERROR laser run  " << fcurrent_run << " tel " << fcurrent_tel << " is not in the VOFFLINE DB (yet?)" << std::endl;
                return false;
            }
            
            Vchannel.push_back( ( unsigned int )atoi( db_row->GetField( 0 ) ) ) ;
            Vmean.push_back( atof( db_row->GetField( 1 ) ) ) ;
            Vvar.push_back( atof( db_row->GetField( 2 ) ) );
        }
    }
    else
    {
    
        cout << "ERROR laser run  " << fcurrent_run << " tel " << fcurrent_tel << " is not in the VOFFLINE DB (yet?)" << std::endl;
        return false;
        
    }
    // HARDWIRED TOTAL NUMBER OF CHANNELS
    if( Vchannel.size() < 499 )
    {
        cout << "ERROR laser run  " << fcurrent_run << " tel " << fcurrent_tel << " has " << Vchannel.size() << " channel filled in the DB. Should be 499" << std::endl;
        return false;
    }
    
    
    
    
    return true;
}
//---------------------------------------------------------------------------
//-- VDB_CalibrationInfo::Create_query_read()
//-- write the query to read gain OR toff from VOFFLINE DB table tblEventDisplay_Analysis_Calibration_Flasher
//-- the most recent entry is taken
//-- a given version
//---------------------------------------------------------------------------
void VDB_CalibrationInfo::Create_query_read()
{

    char c_query[flong_char_query];// has to be long if we ask for a long run list in VDBSourceInfo
    
    //--- gain only
    if( freading_gain_or_toff == 1 )
    {
    
        if( fVOFFLINE_version_query > 0 )
        {
            sprintf( c_query, " SELECT * FROM(SELECT tbl.channel_id, tbl.gain_mean, tbl.gain_var FROM tblEventDisplay_Analysis_Calibration_Flasher AS tbl WHERE tbl.telescope = %d AND  tbl.run_id = %d AND tbl.code_version = %d AND tbl.high_low_gain_flag = %d ORDER BY tbl.update_time DESC ) AS BIG_table GROUP BY channel_id ;", fcurrent_tel, fcurrent_run, fVOFFLINE_version_query, fLOW_GAIN );
        }
        else
        {
            sprintf( c_query, " SELECT * FROM(SELECT tbl.channel_id, tbl.gain_mean, tbl.gain_var FROM tblEventDisplay_Analysis_Calibration_Flasher AS tbl WHERE tbl.telescope = %d AND  tbl.run_id = %d AND tbl.high_low_gain_flag = %d ORDER BY tbl.update_time DESC ) AS BIG_table GROUP BY channel_id;", fcurrent_tel, fcurrent_run, fLOW_GAIN );
        }
        
    }
    //--- toff only
    else if( freading_gain_or_toff == 2 )
    {
    
        if( fVOFFLINE_version_query > 0 )
        {
            sprintf( c_query, " SELECT * FROM(SELECT tbl.channel_id, tbl.toffset_mean, tbl.toffset_var FROM tblEventDisplay_Analysis_Calibration_Flasher AS tbl WHERE tbl.telescope = %d AND  tbl.run_id = %d AND tbl.code_version = %d AND tbl.high_low_gain_flag = %d ORDER BY tbl.update_time DESC )AS BIG_table GROUP BY channel_id ;", fcurrent_tel, fcurrent_run, fVOFFLINE_version_query, fLOW_GAIN );
        }
        else
        {
            sprintf( c_query, " SELECT * FROM(SELECT tbl.channel_id, tbl.toffset_mean, tbl.toffset_var FROM tblEventDisplay_Analysis_Calibration_Flasher AS tbl WHERE tbl.telescope = %d AND  tbl.run_id = %d AND tbl.high_low_gain_flag = %d ORDER BY tbl.update_time DESC )AS BIG_table GROUP BY channel_id ;", fcurrent_tel, fcurrent_run, fLOW_GAIN );
        }
        
        
    }
    else
    {
        std::cout << "VDB_CalibrationInfo::Create_query_read() ERROR:  " << freading_gain_or_toff << " wrong value for freading_gain_or_toff (should be 0 or 1 or2 )" << std::endl;
    }
    
    fquery_read = c_query;
    
    return;
}

//http://dev.mysql.com/doc/refman/5.0/en/example-maximum-column-group-row.html
//SELECT *
//FROM (SELECT *
//FROM shop
//[WHERE conditions]
//ORDER BY price DESC) AS s
//GROUP BY article


// A priori ci-dessous c est la bonne maniere pour recuperer l entree la plus recente pour chaque channel
// ca devrait etre dans le mail a David williams, comme la ligne de commande que je vais probablement utiliser
/*
   SELECT *
   FROM(
         SELECT
	      tbl.channel_id,
	      tbl.gain_mean,
	      tbl.gain_var,
	      tbl.toffset_mean,
	      tbl.toffset_var
	 FROM
	      tblEventDisplay_Analysis_Calibration_Flasher AS tbl
	 WHERE
	      tbl.telescope = my_telescope
	      AND
	      tbl.run_id = my_run_id
	      AND
	      tbl.code_version = my_code_version
	      AND
	      tbl.high_low_gain_flag = my_high_low_gain_flag

	 ORDER BY tbl.update_time DESC

	)
AS BIG_table
    GROUP BY channel_id
*/




