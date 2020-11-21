/*! file updateDBlaserRUN

    (VERITAS only)

*/

#include <TApplication.h>
#include <TGClient.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include <string>
#include <vector>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <bitset>
#include <getopt.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>

#include "VGlobalRunParameter.h"
#include "VExposure.h"

#include "VDB_CalibrationInfo.h"
#include "VDB_Connection.h" // lucie

using namespace std;

//==============================================================================================================
// PARAMETERS DECLARATION
//==============================================================================================================
//--- command line
bool fold_calib = false;

bool fVOFFLINE_time_start_query_bool = false;
string fVOFFLINE_time_start_query = "";

bool fVOFFLINE_version_query_bool = false;
int fVOFFLINE_version_query = 1079;

TString fcalib_gain_tree_name = "tgain_";
TString fcalib_toff_tree_name = "ttoff_";

TString fcalib_gain_tree_name_old = "tGains_";
TString fcalib_toff_tree_name_old = "tToffs_";

int fLOW_GAIN = 0;

int forget_this_run = 63416;


//--- indicate if the function has succeeded
// TO BE INITIALIZED at the beginning of main
bool bool_all_is_downloaded = false;
bool bool_EVNDanalysis_prepared = false;
bool bool_EVNDanalysis_done = false;
bool bool_EVNDcalib_good = false;

//--- hard coded values
TString fServer = "";



int flong_char_query = 10000;
int fmax_number_tel = 4;
//lucie: this calib dir should be an option, if someone else use the script
TString fCalib_dir = "./";
// Password to write in the VOFF DB
string fpass_word = "";
//--- Tool-global-parameters
int fcurrent_run = -10;
int fcurrent_tel = -10;
int fcurrent_date = -10;

TString fNOW = "demain"; // initialized once, when doing get_new_laser_run_list() or  get_laser_run_info_from_DBs(unsigned int arg_run)
TString fNOW_DB = "demain aussi"; // initialized once, when doing get_new_laser_run_list() or  get_laser_run_info_from_DBs(unsigned int arg_run)
//========== Global Variable for the run information
unsigned int  f_VERITAS_DB_LaserRunNumber;//2016.
unsigned int  f_VERITAS_DB_LaserConfigMask;//2016.
unsigned int  f_VERITAS_DB_LaserExclTel;//2016.
unsigned int  f_VERITAS_DB_LaserMiss;//2016.
unsigned int  f_VERITAS_DB_LaserDate;//2016.

//==============================================================================================================
// FUNCTION DECLARATION
//==============================================================================================================
//-- Functions in the main
int parseOptions( int argc, char* argv[] );//2016. to be modified, need to add the run_number
void get_laser_run_info_from_DBs(); //2016.
TString prepare_CalibVOFF_writing(); // 2016.
void write_calib_DB( TString new_laser_run_list_name ); // 2016.
//-------------------------

//-- Active function
//------------------------------- used in get_new_laser_run_list():
unsigned long Check_telmissing_from_VOFFDB_for_one_run( unsigned int VERITAS_DB_LaserRunNumber_i_run, unsigned int VERITAS_DB_LaserConfigMask_i_run, unsigned int VERITAS_DB_LaserExclTel_i_run, vector < unsigned int > VOFFLINE_DB_LaserRunNumber_Tel );                                                                               //2016.
void check_run_calib( vector< int >  ListTel , vector< int >& list_of_valid_tel ); //2016 need.
bool test_gain_toff( TString file_root_name_gain, TString file_root_name_toff ); //2016 need.
bool construct_VOFFLINE_calibration_writing_file( FILE*& wDB_file );                             //2016 need.
//-------------------------------------------------------------------

//-- Accessing the DB
bool read_one_laserRUN_fromVERITAS_DB( unsigned int arg_run, unsigned int& VERITAS_DB_LaserRunNumber, unsigned int& VERITAS_DB_LaserConfigMask, unsigned int& VERITAS_DB_LaserExclTel, unsigned int& VERITAS_DB_LaserDate ); //2016.
bool read_one_laserRUN_fromVOFFLINE_DB( unsigned int arg_run, vector < unsigned int >& VOFFLINE_DB_LaserRunNumber_Telnum, vector < string >& VOFFLINE_DB_LaserDate_Telnum, vector < string >& VOFFLINE_DB_LaserVersion_Telnum, unsigned int Tel_num ); //2016.
string WriteQuery_to_get_one_LaserRun_fromVERITAS_DB( unsigned int arg_run ); //2016.
string WriteQuery_to_get_one_LaserRun_fromVOFFLINE_DB( unsigned int arg_run, unsigned int Tel_num ); //2016.
//=== tools functions ==========================
void Set_the_time();
//-- reading
bool read_calib_file( TString file_root_name, vector < double >& Vchannel, vector < double >& Vmean, vector < double >& Vvar, bool bgain );//2016 need.
//-- telescope mask related (2016 need)
vector<int> get_list_of_telescope( unsigned int tel_mask, unsigned int excluded_tel );
vector<int> get_list_of_telescope( unsigned int missing_tel_mask );
bool Tel_is_in_mask( Int_t Tel, Int_t config_mask );
bool Tel_is_excluded( Int_t Tel, Int_t config_mask, Int_t excluded_tel );
bool Incoherence_between_ConfigMask_and_excludedTel( Int_t Tel, Int_t config_mask, Int_t excluded_tel );
bool Tel_code_not_valid( Int_t Tel_code );
//-- reading DB time (could be integrated to VExposure) (2016 need)
int get_date_from_tblRun_Info_data_start_time( string iTemp );

//-- File existence and filing (2016 need)
bool does_file_exist( TString file_path_name );
bool is_file_empty( TString file_path_name );

//-- Reading TTree (2016 need)
TH1F* get_h_from_TTree( TFile* file_root, TString string_arbre, TString want, TString condition, TString histo_name, TString histo_bin, bool normalised );
bool fill_3V_from_TTree( TFile* file_root, TString string_arbre, TString want, TString condition, vector < double >& Varg1, vector < double >& Varg2, vector < double >& Varg3 );

//-- LOCAL
//-- constructing name
TString get_name_for_download( TString name );
TString get_name_run_to_be_analysed( TString name );
TString get_name_run_already_analysed( TString name );
TString get_name_calib_file( TString dir, TString suffixe, int tel_num, unsigned int RunNumber );
TString get_name_run_problem( TString name );
TString get_name_run_to_be_writen( TString name, bool long_table );
TString get_downloaded_laser_path( unsigned int RunNumber, unsigned int Date );


//==============================================================================================================
// FUNCTION DEFINITION
//==============================================================================================================


//------------------------------------
//-- main
//------------------------------------
int main( int argc, char* argv[] )
{

    //--- indicate if the function has succeeded
    bool_all_is_downloaded = false;
    bool_EVNDanalysis_prepared = false;
    bool_EVNDanalysis_done = false;
    bool_EVNDcalib_good = false;
    
    
    // 2016. to do: add the run number in the argument of the function
    //              make sure the option for the path to the directory where the calib results are is working
    parseOptions( argc, argv );
    
    
    
    // 2016 to do: would be great to read those name from the current EVNDISP version...
    std::cout << "REMINDER: Calibration checking, HARD coded name TTree read are: " << fcalib_gain_tree_name << " and " << fcalib_toff_tree_name << std::endl;
    // they are the same for the all run list, which have to be homogeneous in this respect
    // name could be change at the parse option level (but for this would be good to have option instroduced with more than one letter).
    
    
    //VGlobalRunParameter* blah = new VGlobalRunParameter() ;
    //fServer = blah->getDBServer();
    cout << "Using server " << fServer << endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "ATTENTION: make sure the Calib_dir is correctly defined (option -d). This is were the output of the .gain .ped ,toff will be " << std::endl;
    std::cout << "in this directory there should be one directory per telescope (Tel_1, Tel_2, ...) " << std::endl;
    std::cout << "Calib_dir  " << fCalib_dir << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    Set_the_time();
    // 2016. getting the info about the laser run from the DB
    std::cout << "Getting list of ALL laser flasher run from VERITAS DB " << std::endl;
    get_laser_run_info_from_DBs(); //2016. to be checked (to do)
    // 2016. Filing the following global variables
    //f_VERITAS_DB_LaserRunNumber;
    //f_VERITAS_DB_LaserConfigMask;
    //f_VERITAS_DB_LaserExclTel;
    //f_VERITAS_DB_LaserMiss;
    //f_VERITAS_DB_LaserDate;
    
    //=======================================
    // still need to do the get, some how, read the DB for the laser run I'm dealing with, because need the mask info and all that. All the info collected, I believe in the -g
    //=======================================
    /*
      Verification and writing, to be adapted
      - run by run bases
      - we know the run number
      - we know where the calibration info from EVNDISP analysis is
     */
    //=======================================
    
    //------------ test the output of the analysis and write the file to be used to fill the VOFFLine DB
    TString file_name_write_DB = prepare_CalibVOFF_writing(); // done: to be tested (seems ok, but problem with TTree named changed during an EVNDISP update, have to rerun all runs)
    //------------
    
    write_calib_DB( file_name_write_DB );
    
    
    // 2016 to do: delete the writing file file_name_write_DB;
    
    std::remove( file_name_write_DB.Data() );
    
    
    
    
    
}




//===============================================================================================
//==== MAIN FUNCTIONS
//===============================================================================================


//---------------------------------------------------------------------------
//-- 2016. get_laser_run_info_from_DBs()
//-- fill the global variable with the run info for the DBs
//-- it's only for the run we asked for, and this run might not be an official laser run.
//-- laser run present in VERITAS DB but not in VOFFLINE
//---------------------------------------------------------------------------
void get_laser_run_info_from_DBs()
{
    //-------------------------------------------------------------------
    //--- get the info for the laser run from VERITAS DB
    // 2016. need to modify this function. Get the same info, but for a given laser_run, regardless if it is an official laser run or not
    unsigned int VERITAS_DB_LaserRunNumber;
    unsigned int VERITAS_DB_LaserConfigMask;
    unsigned int VERITAS_DB_LaserExclTel;
    unsigned int VERITAS_DB_LaserDate;
    std::cout << "read laserRUN " << fcurrent_run << " fromVERITAS_DB " << std::endl;
    if( !read_one_laserRUN_fromVERITAS_DB( fcurrent_run, VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VERITAS_DB_LaserDate ) )
    {
        exit( 0 );
    }
    // 2016.
    
    
    //--- get the list of laser run from VOFFSET (for a given analysis version)
    vector < unsigned int >  VOFFLINE_DB_LaserRunNumber_Tel;
    vector < string >  VOFFLINE_DB_Laserdate_Tel;
    vector < string >  VOFFLINE_DB_Laserversion_Tel;
    //max number of telescope is hard coded here
    // and the oder in which the telescope are called is important
    if( !read_one_laserRUN_fromVOFFLINE_DB( fcurrent_run, VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 1 ) )
    {
        exit( 0 );
    }
    if( !read_one_laserRUN_fromVOFFLINE_DB( fcurrent_run, VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 2 ) )
    {
        exit( 0 );
    }
    if( !read_one_laserRUN_fromVOFFLINE_DB( fcurrent_run, VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 3 ) )
    {
        exit( 0 );
    }
    if( !read_one_laserRUN_fromVOFFLINE_DB( fcurrent_run, VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 4 ) )
    {
        exit( 0 );
    }
    
    
    //--- compare the list of laser run from VERITAS DB with the one from VOFFLINE
    std::cout << "Checking if this laser/flasher run is missing from VOFFLINE DB " << std::endl;
    //---- 2016. have a look, see if I need to modify this function. probably not, was in the loop
    unsigned long mask_missing_tel = Check_telmissing_from_VOFFDB_for_one_run( VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VOFFLINE_DB_LaserRunNumber_Tel );
    
    if( mask_missing_tel > 0 )
    {
        std::cout << " " << std::endl;		 //--- fill the global variables with: run mask exclu_tel mask_with_missing_tel_from_VOFFLINE_DB
        printf( "*** %d %2d %2d %2lu %d \n", VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, mask_missing_tel, VERITAS_DB_LaserDate );
        f_VERITAS_DB_LaserRunNumber = VERITAS_DB_LaserRunNumber;
        f_VERITAS_DB_LaserConfigMask = VERITAS_DB_LaserConfigMask;
        f_VERITAS_DB_LaserExclTel = VERITAS_DB_LaserExclTel;
        f_VERITAS_DB_LaserMiss = mask_missing_tel;
        f_VERITAS_DB_LaserDate = VERITAS_DB_LaserDate;
        std::cout << " " << std::endl;		 //--- fill the NEW_laser_list_File with: run mask exclu_tel mask_with_missing_tel_from_VOFFLINE_DB
    }
    else
    {
        f_VERITAS_DB_LaserRunNumber = 0;
        f_VERITAS_DB_LaserConfigMask = 0;
        f_VERITAS_DB_LaserExclTel = 0;
        f_VERITAS_DB_LaserMiss = 0;
        f_VERITAS_DB_LaserDate = 0;
    }
    
    return;
}




//---------------------------------------------------------------------------
//-- Look the output file of EVND calibration
//-- create a list of run ready to be copied in the DB
//-- create a file with info on all run
//---------------------------------------------------------------------------
TString prepare_CalibVOFF_writing()
{


    //-- Create File to write in DB, with valid run
    TString  wDB_file_name = "write_file";
    FILE* wDB_file;
    wDB_file = fopen( wDB_file_name.Data(), "w" );
    //-- Create File which collect all checking information on the run
    TString  pbRunFile_name = "pb_file";
    //---- 2016. Can keep it, even on the run by run basis. Or could simply print it out.
    
    //--- read File run list missing from VOFF DB and store information
    unsigned int VERITAS_DB_LaserRunNumber  = f_VERITAS_DB_LaserRunNumber;
    //int VERITAS_DB_LaserConfigMask = f_VERITAS_DB_LaserConfigMask;
    //unsigned int 	int VERITAS_DB_LaserExclTel    = f_VERITAS_DB_LaserExclTel;
    unsigned int VERITAS_DB_LaserMiss       = f_VERITAS_DB_LaserMiss;
    unsigned int VERITAS_DB_LaserDate       = f_VERITAS_DB_LaserDate;
    
    string  Summary_tel_list;
    string  Summary_calib_status;
    string  Summary_writing_status;
    //----------------
    
    
    //==== 2016 will only do once, what is in this loop
    
    if( fcurrent_run != ( int ) f_VERITAS_DB_LaserRunNumber )
    {
        std::cout << "PROBLEM run asked for " << fcurrent_run << " does not correspond to the run read in DB " << f_VERITAS_DB_LaserRunNumber << std::endl;
        exit( 0 );
    }
    fcurrent_date = VERITAS_DB_LaserDate; // 2016 don't really need that
    //--------------------
    // for summary
    string tel_list = "";
    string telgoodcalib = "-1";
    string telgoodwriting = "-1";
    //--------------------
    
    //----------------------------------------------------------------------------
    //-- Check calibration
    //----------------------------------------------------------------------------
    vector< int > list_tel_missing = get_list_of_telescope( VERITAS_DB_LaserMiss );
    vector< int > list_of_valid_tel;
    // 2016. Run wise function, probably don't have to change much
    check_run_calib( list_tel_missing, list_of_valid_tel );
    
    //--------------------------------------
    // for summary
    //-- result of the check ----- (just to have a look)
    if( list_of_valid_tel.size() > 0 )
    {
        if( list_tel_missing.size() < list_of_valid_tel.size() )
        {
            std::cout << "ERROR: BIG Problem !!!! check what the code is doing (RUN:" << fcurrent_run << ")" << std::endl;
            return wDB_file_name;
        }
        else if( list_tel_missing.size() == list_of_valid_tel.size() )
        {
            printf( "                       %5d run NO calib problem \n", ( int )fcurrent_run );
        }
        else if( list_tel_missing.size() > list_of_valid_tel.size() )
        {
            printf( "                       %5d run with calib problem for some telescopes\n", ( int ) fcurrent_run );
        }
    }
    else
    {
        printf( "                       %5d run with calib problem for ALL telescopes\n", ( int ) fcurrent_run );
    }
    //--------------------------------------
    
    //----------------------------------------------------------------------------
    //-- Construct Writing DB File if there is some telescopes with valid calib
    //----------------------------------------------------------------------------
    unsigned int counter = 0; // for summary
    if( list_of_valid_tel.size() > 0 )
    {
    
        telgoodcalib = "";// for summary
        
        //-- Loop on the telescope with valid calibration
        for( unsigned int i = 0; i < list_of_valid_tel.size() ; i++ )
        {
            fcurrent_tel = list_of_valid_tel[i];
            //---------------
            // for summary
            char ctelcal[10];
            sprintf( ctelcal, "%d", fcurrent_tel );
            string stelcal = ctelcal;
            telgoodcalib += stelcal;
            //---------------
            bool writing_ok = construct_VOFFLINE_calibration_writing_file( wDB_file );
            //-------------------------
            // for summary
            if( writing_ok )
            {
                if( counter > 0 )
                {
                    char ctel[10];
                    sprintf( ctel, "%d", fcurrent_tel );
                    string stel = ctel;
                    telgoodwriting += stel;
                    counter++;
                }
                else
                {
                    counter++;
                    char ctel[10];
                    sprintf( ctel, "%d", fcurrent_tel );
                    string stel = ctel;
                    telgoodwriting = stel;
                }
            }
            //-------------------------
        }
        //-- end loop tel
        //-------------------------
        // for summary
        if( list_of_valid_tel.size() > 0 )
        {
            if( list_tel_missing.size() < counter )
            {
                std::cout << "ERROR: BIG Problem writing !!!! check what the code is doing (RUN:" << fcurrent_run << ")" << std::endl;
                return wDB_file_name ;
            }
            else if( list_tel_missing.size() == counter )
            {
                printf( "                       %5d run with no problem writing\n", ( int ) fcurrent_run );
            }
            else if( list_tel_missing.size() > counter )
            {
                printf( "                       %5d run with some problems writing\n", ( int ) fcurrent_run );
            }
        }
        else
        {
            printf( "                       %5d run with problems writing\n", ( int ) fcurrent_run );
        }
        //-------------------------
    }
    //-- end constructing file
    //-------------------------
    //  for summary
    //-- get the list of telescope tested, for comparison
    for( unsigned int j = 0 ; j < list_tel_missing.size(); j++ )
    {
        char ctel[10];
        sprintf( ctel, "%d", list_tel_missing[j] );
        string stel = ctel;
        tel_list += stel;
    }
    Summary_tel_list =  tel_list ;
    Summary_calib_status =  telgoodcalib ;
    Summary_writing_status = telgoodwriting ;
    //-------------------------
    
    if( ( int ) fcurrent_run == ( int ) forget_this_run )
    {
        return wDB_file_name ;
    }
    
    TString calib_status = "BAD";
    TString writing_status = "BAD";
    if( Summary_tel_list.compare( Summary_calib_status ) == 0 )
    {
        calib_status = "GOOD";
    }
    if( Summary_tel_list.compare( Summary_writing_status ) == 0 )
    {
        writing_status = "GOOD";
    }
    
    TString Stellist = Summary_tel_list;
    TString Scalibstatus = Summary_calib_status;
    TString Swritingstatus = Summary_writing_status;
    
    printf( "d%8d %d %4s %4s %4s %4s %4s \n", VERITAS_DB_LaserDate, VERITAS_DB_LaserRunNumber, calib_status.Data(), writing_status.Data(), Stellist.Data(), Scalibstatus.Data(), Swritingstatus.Data() );
    //-----------------------------------------------------------------------------------------
    
    //-----------------------------------------------------------------------------------------
    //-- tells that the process of checking and writing the file is over (not that it was successful) AND that there is something to write
    if( !is_file_empty( wDB_file_name ) )
    {
        bool_EVNDcalib_good = true;
    }
    else
    {
        std::cout << "WRITING DB FILE is empty " << std::endl;
    }
    //-----------------------------------------------------------------------------------------
    
    
    fclose( wDB_file );
    
    return wDB_file_name;
}

//---------------------------------------------------------------------------
//-- write_calib_DB
//---------------------------------------------------------------------------
void write_calib_DB( TString wDB_file_name )
{
    VDB_CalibrationInfo* db_calib_info = new VDB_CalibrationInfo( wDB_file_name, fServer, fLOW_GAIN );
    db_calib_info->write_inVOFFLINE_DB_from_file( fpass_word );
    
    return;
}



//===========================================================================
//====== ACTIVE functions for look_at_evndisp_calibration_output



//---------------------------------------------------------------------------
//-- void construct_VOFFLINE_calibration_writing_file
//-- the list_of_valid_tel correspond to the list of telescope for which the calibration is ok according to check_run_calib
//-- independent function from checking calib output, clearer to see what is copied in the DB
//---------------------------------------------------------------------------
bool construct_VOFFLINE_calibration_writing_file( FILE*& wDB_file )
{


    bool construction_ok = false;
    
    
    printf( "writing ---Run %d Tel %4d \n", fcurrent_run, fcurrent_tel );
    //    fprintf(wDB_file, "--- Tel %d \n",fcurrent_tel);
    //-- get the name of the calibration root file
    TString file_gain_root = get_name_calib_file( fCalib_dir, "gain.root", fcurrent_tel, fcurrent_run );
    TString file_toff_root = get_name_calib_file( fCalib_dir, "toff.root", fcurrent_tel, fcurrent_run );
    
    vector < double > Vchannel_gain;
    vector < double > Vmean_gain;
    vector < double > Vvar_gain;
    vector < double > Vchannel_toff;
    vector < double > Vmean_toff;
    vector < double > Vvar_toff;
    
    
    if( read_calib_file( file_gain_root, Vchannel_gain, Vmean_gain, Vvar_gain, true )
            && read_calib_file( file_toff_root, Vchannel_toff, Vmean_toff, Vvar_toff, false ) )
    {
    
        VDB_CalibrationInfo* db_calib_info = new VDB_CalibrationInfo( fcurrent_run, fcurrent_tel, fNOW_DB, fVOFFLINE_version_query, fLOW_GAIN );
        construction_ok = db_calib_info->DoFile_for_DBwriting( Vchannel_gain, Vmean_gain, Vvar_gain, Vchannel_toff, Vmean_toff, Vvar_toff, wDB_file );
        //construction_ok = DoFile_for_DBwriting(Vchannel_gain,Vmean_gain,Vvar_gain,Vchannel_toff,Vmean_toff,Vvar_toff,wDB_file);
        delete db_calib_info;
    }
    
    
    return construction_ok;
}



//-------------------------------------------------------------
int parseOptions( int argc, char* argv[] )
{


    while( 1 )
    {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"cal_dir", required_argument, 0, 'd'},
            {"pass_word", required_argument, 0, 'p'},
            {"run_number", required_argument, 0, 'r'},
            {"server", required_argument, 0, 's'},
            {"time_start", required_argument, 0, 't'},
            {"req_version", required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };
        
        int option_index = 0;
        int c = getopt_long( argc, argv, "hd:p:r:s:t:v:", long_options, &option_index );
        if( argc == 1 )
        {
            c = 'h';
        }
        if( c == -1 )
        {
            break;
        }
        switch( c )
        {
            case 0:
                //                if(long_options[option_index].flag != 0)
                //                    break;
                //                printf ("option %s", long_options[option_index].name);
                //                if (optarg)
                //                    printf (" with arg %s", optarg);
                printf( "\n" );
                break;
            case 'h':
                // EVNDISP standard way to show the option without having to recompile the whole thing
                // should work together with a script that call the binary of this compiled file
                // the scripts has to call set observatory
                char* ENV;
                ENV = getenv( "OBS_EVNDISP_AUX_DIR" );
                char readme[500];
                sprintf( readme, "cat %s/ParameterFiles/EVNDISP.updateDBlaserRUN.runparameter", ENV );
                if( system( readme ) < 0 ) cout << "error in finding readme" << endl;
                exit( 0 );
                break;
            case 'd':
                fCalib_dir = optarg;
                break;
            case 'p':
                fpass_word = optarg;
                break;
            case 'r':
                sscanf( optarg, "%d", &fcurrent_run );
                break;
            case 's':
                fServer = optarg;
                break;
            case 't':
                fVOFFLINE_time_start_query_bool = true;
                fVOFFLINE_time_start_query = optarg; // year-month-day: 0000-00-00 OR year-month-day hour:minute:second 0000-00-00 00:00:00
                break;
            case 'v':
                fVOFFLINE_version_query_bool = true;
                sscanf( optarg, "%d", &fVOFFLINE_version_query );
                break;
            default:
                abort();
        }
    }
    
    std::cout << "fCalib_dir " << fCalib_dir << std::endl;
    std::cout << "fpass_word " << fpass_word << std::endl;
    std::cout << "fcurrent_run " << fcurrent_run << std::endl;
    std::cout << "fServer " << fServer << std::endl;
    std::cout << "fVOFFLINE_time_start_query_bool " << fVOFFLINE_time_start_query_bool << std::endl;
    std::cout << "fVOFFLINE_time_start_query " << fVOFFLINE_time_start_query << std::endl;
    std::cout << "fVOFFLINE_version_query_bool " << fVOFFLINE_version_query_bool << std::endl;
    std::cout << "fVOFFLINE_version_query_bool " << fVOFFLINE_version_query_bool << std::endl;
    
    
    
    
    return optind;
}

//===============================================================================================
//==== ACTIVE FUNCTIONS
//===============================================================================================


//---------------------------------------------------------------------------
//-- bool test_gain_toff
//-- gain and toff should be tested together, because need to compare the channel
//-- for example for dead pixel
//---------------------------------------------------------------------------
bool test_gain_toff( TString file_root_name_gain, TString file_root_name_toff )
{

    std::cout << "test_gain_toff , root file: " << file_root_name_gain << std::endl;
    
    bool good_calib = true;
    
    TString want = "";
    TString condition = "";
    TString histo_name = "name";
    TString histo_bin = "";
    
    std::vector<Int_t> V_pixel_zero_gain;
    std::vector<Int_t> V_pixel_zero_toff;
    
    //----------------------------------
    //-- test gain
    TFile* root_file_gain = new TFile( file_root_name_gain );
    
    if( root_file_gain->IsOpen() )
    {
    
        //   std::cout<<"root file is open "<<std::endl;
        
        TString tree_name = fcalib_gain_tree_name;//"tGains_";
        if( fcurrent_run < 63416 && fold_calib )
        {
            tree_name = fcalib_gain_tree_name_old;
        }
        
        tree_name += fcurrent_tel;
        //-- look at gain RMS
        want = "gain";
        histo_name = "h_gain";
        histo_bin = "(500, 0, 10)";
        TH1F* h_gain = get_h_from_TTree( root_file_gain, tree_name, want, condition, histo_name, histo_bin, false );
        
        if( h_gain->GetRMS() > 0.6 )
        {
            printf( "ERROR: gain RMS %f too large\n", h_gain->GetRMS() );
            good_calib = false;
        }
        //-- look at the number of channel with gain = 0
        want = "channel";
        condition = "gain == 0";
        histo_name = "h_channel";
        histo_bin = "(499, 0, 499)";
        TH1F* h_gz_channel = get_h_from_TTree( root_file_gain, tree_name, want, condition, histo_name, histo_bin, false );
        if( h_gz_channel->GetEntries() >= 499 )
        {
            printf( "ERROR: %4f channels with gain at zero \n", h_gz_channel->GetEntries() );
            good_calib = false;
        }
        else if( h_gz_channel->GetEntries() > 0 && h_gz_channel->GetEntries() < 499 )
        {
            //- list channel with gain = 0
            //-- loop on the pixel
            for( int ii = 1 ; ii <= 500; ii++ )
            {
                if( h_gz_channel->GetBinContent( ii ) > 0 )
                {
                    V_pixel_zero_gain.push_back( TMath::Nint( h_gz_channel->GetBinCenter( ii ) - 0.1 ) );
                }
            }
        }
        
        // close root file
        root_file_gain->Close(); // attention a cela, mes histo peuvent etre vide parce que je ferme le fichier root
    }
    else
    {
        printf( "ERROR: File %s does not exist \n", file_root_name_gain.Data() );
        good_calib = false;
    }
    //-------------------------------
    //-- test toff
    TFile* root_file_toff = new TFile( file_root_name_toff );
    if( root_file_toff->IsOpen() )
    {
    
        TString tree_name = fcalib_toff_tree_name;//"tToffs_";
        if( fcurrent_run < 63416 && fold_calib )
        {
            tree_name = fcalib_toff_tree_name_old;
        }
        tree_name += fcurrent_tel;
        //-- look at the number of channel with toff = 0
        want = "channel";
        condition = "toff == 0";
        histo_bin = "(499, 0, 499)";
        TH1F* h_tz_channel = get_h_from_TTree( root_file_toff, tree_name, want, condition, histo_name, histo_bin, false );
        if( h_tz_channel->GetEntries() >= 499 )
        {
            printf( "ERROR: %4f channels with toff at zero \n", h_tz_channel->GetEntries() );
            good_calib = false;
        }
        else if( h_tz_channel->GetEntries() > 0 && h_tz_channel->GetEntries() < 499 )
        {
            //- list channel with toff = 0
            //-- loop on the pixel
            for( int ii = 1 ; ii <= 500; ii++ )
            {
                if( h_tz_channel->GetBinContent( ii ) > 0 )
                {
                    V_pixel_zero_toff.push_back( TMath::Nint( h_tz_channel->GetBinCenter( ii ) - 0.1 ) );
                }
            }
            // we will compare with the number of channel with toff = 0 later
            // it's only a problem if the number of entries are not the same
            // not sure I understand this: if(h_gz_channel->GetBinContent(ii)>0)
        }
        
        // close root file
        root_file_toff->Close(); // attention a cela, mes histo peuvent etre vide parce que je ferme le fichier root
    }
    else
    {
        printf( "ERROR: File %s does not exist \n", file_root_name_toff.Data() );
        good_calib = false;
    }
    //----------------------- compare channel's gain and toff at zero
    if( V_pixel_zero_gain.size() > 0 || V_pixel_zero_toff.size() > 0 )
    {
        if( V_pixel_zero_gain.size() == V_pixel_zero_toff.size() )
        {
            printf( "Dead pixels: \n" );
            for( unsigned int j = 0 ; j < V_pixel_zero_gain.size(); j++ )
            {
                printf( " %d ", V_pixel_zero_gain[j] );
            }
            printf( " \n" );
            
        }
        else
        {
            if( TMath::Abs( ( Int_t ) V_pixel_zero_gain.size() - ( Int_t ) V_pixel_zero_toff.size() ) > 2 )
            {
                printf( "ERROR: NOT same number of pixel with zero gain (%4d) and zero toff (%4d) \n", ( int ) V_pixel_zero_gain.size(), ( int ) V_pixel_zero_toff.size() );
                good_calib = false;
            }
            else
            {
                printf( "WARNING: NOT same number of pixel with zero gain (%4d) and zero toff (%4d) \n", ( int ) V_pixel_zero_gain.size(), ( int ) V_pixel_zero_toff.size() );
                good_calib = true;
            }
            for( unsigned int ii = 0; ii < V_pixel_zero_gain.size(); ii++ )
            {
                printf( " %4d zero gain \n", V_pixel_zero_gain[ii] );
            }
            for( unsigned int ii = 0; ii < V_pixel_zero_toff.size(); ii++ )
            {
                printf( " %4d zero toff \n", V_pixel_zero_toff[ii] );
            }
        }
    }
    //-------------------------------------------------------------------
    return good_calib;
}



//---------------------------------------------------------------------------
//-- unsigned long Check_telmissing_from_VOFFDB
//---------------------------------------------------------------------------
//-- used in TString get_laser_run_info_from_DBs(unsigned int arg_run)
//-- not fit to be called on its own
//-------------------
//-- Check_telmissing_from_VOFFDB
//-- read VOFFLINE_DB_LaserRunNumber_Tel, which should be properly filled, as it is done in get_laser_run_info_from_DBs(unsigned int arg_run)
//-- Check if laser run in VERITAS DB are present in VOFFLINE DB
//-- return mask of telescopes missing from VOFFLINE DB (which should be there)
//---------------------------------------------------------------------------
unsigned long Check_telmissing_from_VOFFDB_for_one_run( unsigned int VERITAS_DB_LaserRunNumber_i_run, unsigned int VERITAS_DB_LaserConfigMask_i_run, unsigned int VERITAS_DB_LaserExclTel_i_run, vector < unsigned int > VOFFLINE_DB_LaserRunNumber_Tel )
{
    //	 std::cout<<"Check_telmissing_from_VOFFDB_for_one_run "<<std::endl;
    //	 printf ("-------- %d %d %d \n",VERITAS_DB_LaserRunNumber_i_run,VERITAS_DB_LaserConfigMask_i_run, VERITAS_DB_LaserExclTel_i_run);
    bitset<4> missing_tel; // with this we list the telescope missing from VOFFLINE DB
    missing_tel.reset();
    
    if( fcurrent_run != ( int ) VERITAS_DB_LaserRunNumber_i_run )
    {
        std::cout << "PROBLEM: the run given as argument (" << fcurrent_run << ") is not the one read in the DB (" << VERITAS_DB_LaserRunNumber_i_run << ")!!! " << std::endl;
        exit( 0 );
    }
    
    vector<int> list_of_tel = get_list_of_telescope( VERITAS_DB_LaserConfigMask_i_run, VERITAS_DB_LaserExclTel_i_run );
    
    
    //--- loop on the telescopes participating to the run
    for( unsigned i_tel = 0; i_tel < list_of_tel.size(); i_tel++ )
    {
        int tel_indice = list_of_tel[i_tel] - 1;
        std::cout << "Tel  " << list_of_tel[i_tel] << " tel_indice " << tel_indice << std::endl;
        //--- loop on the run in VOFFLINE DB
        bool run_in_voff_db = false;
        if( VERITAS_DB_LaserRunNumber_i_run == VOFFLINE_DB_LaserRunNumber_Tel[tel_indice] )
        {
            //		    std::cout<<"VOFFLINE == VERITAS "<<std::endl;
            run_in_voff_db = true;
            break;
        }
        
        //--- end loop on voff run for a given tel
        if( !run_in_voff_db )
        {
            std::cout << "RUN " << fcurrent_run << ": miss tel:" << missing_tel.set( tel_indice, 1 ) << " in VOFFLINE" << std::endl ;     // construction of the mask for the telescopes missing from VOFFLINE DB
        }
    }
    //--- end loop on telescope
    return missing_tel.to_ulong();
}



//---------------------------------------------------------------------------
//-- 2016. read_one_laserRUN_fromVOFFLINE_DB                                         --
//-- for each telescope, get the list of run in the data base              --
// dont need the first layer of vector, should get rid of them
//---------------------------------------------------------------------------
bool read_one_laserRUN_fromVOFFLINE_DB( unsigned int arg_run, vector < unsigned int >& VOFFLINE_DB_LaserRunNumber_Telnum, vector < string >& VOFFLINE_DB_LaserDate_Telnum, vector < string >& VOFFLINE_DB_LaserVersion_Telnum, unsigned int Tel_num )
{


    //	 std::cout<<"read_one_laserRUN_fromVOFFLINE_DB "<<std::endl;
    string query = WriteQuery_to_get_one_LaserRun_fromVOFFLINE_DB( arg_run, Tel_num );
    
    bool things_went_well = false;
    
    // DONE modified this function considering it should be for only one run!!!! first layer of vector to be removed, and make sure we have only one line
    // OK   actually. how do I get only one line? even for one run? (multiple version, each pixel???) => Grouped by run_id, and ordered by date
    // if more than one entry for a given telescope and given run, at smoe point, probably chose to keep only the most recent one. Have to find where this is done!!!!
    // probably after, when reading those vectors, only taking the first entrance for a given run_number/tel_number
    
    //---- open the DB and check
    string iTempS;
    iTempS =  fServer + "/VOFFLINE";
    //TSQLServer *f_db = TSQLServer::Connect( iTempS.c_str(), "readonly", "" );
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ; // lucie: DO NOT create a pointer. this way the connection is automatically closed when getting out of the function
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "ERROR read_one_laserRUN_fromVOFFLINE_DB: failed to connect to database server" << endl;
        cout << "\t server: " <<  fServer << endl;
        return things_went_well;
    }
    //---- do the query and check
    if( !my_connection.make_query( query.c_str() ) )
    {
        cout << "WARNING read_one_laserRUN_fromVOFFLINE_DB: failed to get something from the query " << endl;
        return things_went_well ;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    //---- read the query
    if( db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( atoi( db_row->GetField( 0 ) ) == forget_this_run )
            {
                return things_went_well ;
            }
            
            VOFFLINE_DB_LaserRunNumber_Telnum.push_back( atoi( db_row->GetField( 0 ) ) );
            string version = string( db_row->GetField( 1 ) );
            VOFFLINE_DB_LaserVersion_Telnum.push_back( version ) ;
            string date = string( db_row->GetField( 2 ) );
            VOFFLINE_DB_LaserDate_Telnum.push_back( date ) ;
        }
    }
    else
    {
        cout << "WARNING read_one_laserRUN_fromVOFFLINE_DB:  no laser run found " << endl;
        VOFFLINE_DB_LaserRunNumber_Telnum.push_back( -1 );
        VOFFLINE_DB_LaserVersion_Telnum.push_back( "toto" ) ;
        VOFFLINE_DB_LaserDate_Telnum.push_back( "now" ) ;
        return true ;
    }
    
    things_went_well = true;
    
    //-- close the DB
    //f_db->Close(); // lucie
    return things_went_well ;
}

//---------------------------------------------------------------------------
//-- 2016. read_one_laserRUN_fromVERITAS_DB                                          --
//-- get the list of all valid laser run in the data base                  --
//-- should be only one line to read, so, no need for vector
//---------------------------------------------------------------------------
//void read_one_laserRUN_fromVERITAS_DB(int arg_run, vector< unsigned int >& VERITAS_DB_LaserRunNumber, vector< unsigned int >& VERITAS_DB_LaserConfigMask, vector< unsigned int >& VERITAS_DB_LaserExclTel, vector< unsigned int >& VERITAS_DB_LaserDate )
bool read_one_laserRUN_fromVERITAS_DB( unsigned int arg_run, unsigned int& VERITAS_DB_LaserRunNumber, unsigned int& VERITAS_DB_LaserConfigMask, unsigned int& VERITAS_DB_LaserExclTel, unsigned int& VERITAS_DB_LaserDate )
{

    bool things_went_well = false;
    
    
    //---- open the DB and check
    string iTempS;
    iTempS =  fServer;
    iTempS += "/VERITAS";
    std::cout << "server:  " << iTempS << std::endl;
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ; // lucie: DO NOT create a pointer. this way the connection is automatically closed when getting out of the function
    if( !my_connection.Get_Connection_Status() )
    {
        //if( !f_db ){
        cout << "ERROR read_one_laserRUN_fromVERITAS_DB: failed to connect to database server" << endl;
        cout << "\t server: " <<  fServer << endl;
        return things_went_well;
    }
    string query = WriteQuery_to_get_one_LaserRun_fromVERITAS_DB( arg_run );
    
    //---- do the query and check
    if( !my_connection.make_query( query.c_str() ) )
    {
        return things_went_well;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    //---- read the query
    if( db_res->GetRowCount() > 0 )
    {
    
        if( db_res->GetRowCount() > 1 )
        {
            cout << " OOOOhhhhhhh more than one line out???? should not happen!!!!!" << std::endl;
            return things_went_well;
        }
        
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( db_row->GetField( 0 ) )
            {
                VERITAS_DB_LaserRunNumber =  atoi( db_row->GetField( 0 ) ) ;
            }
            else
            {
                std::cout << "WARNING: no run number  in VERITAS DB" << std::endl;
                VERITAS_DB_LaserRunNumber =  -1 ;
                return things_went_well;
            }
            if( db_row->GetField( 1 ) )
            {
                VERITAS_DB_LaserConfigMask = atoi( db_row->GetField( 1 ) ) ;
            }
            else
            {
                std::cout << "WARNING: no config Mask in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                VERITAS_DB_LaserConfigMask = -1 ;
                return things_went_well;
            }
            if( db_row->GetField( 2 ) )
            {
                VERITAS_DB_LaserExclTel = atoi( db_row->GetField( 2 ) ) ;
            }
            else
            {
                std::cout << "WARNING: no exclud Tel  in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                VERITAS_DB_LaserExclTel = -1 ;
                return things_went_well;
            }
            if( db_row->GetField( 3 ) )
            {
                VERITAS_DB_LaserDate = get_date_from_tblRun_Info_data_start_time( db_row->GetField( 3 ) ) ;
            }
            else
            {
                std::cout << "WARNING: no data start date in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                std::cout << "using db_start_time instead" << std::endl;
                if( db_row->GetField( 4 ) )
                {
                    VERITAS_DB_LaserDate = get_date_from_tblRun_Info_data_start_time( db_row->GetField( 4 ) ) ;
                }
                else
                {
                    std::cout << "WARNING: no db start date in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                    VERITAS_DB_LaserDate = -1 ;
                    return things_went_well;
                }
            }
        }
    }
    else
    {
        cout << "WARNING read_one_laserRUN_fromVERITAS_DB:  no laser run found " << endl;
        return things_went_well;
    }
    //-- close the DB
    //f_db->Close(); // lucie (closing automatically with the automatic deletion of my_connection)
    
    things_went_well = true;
    return things_went_well;
    
}

//---------------------------------------------------------------------------
//--2016. WriteQuery_to_get_one_LaserRun_info_fromVERITAS_DB() get out same info than WriteQuery_to_getLaserRun_fromVERITAS_DB()
//-- but for a given run number, regardless if it's an official flasher or not
//---------------------------------------------------------------------------
string WriteQuery_to_get_one_LaserRun_fromVERITAS_DB( unsigned int arg_run )
{



    //begining of the querry
    string query1 = "SELECT big_table.run_id, big_table.config_mask, big_table.excluded_telescopes, big_table.data_start_time , big_table.db_start_time FROM (SELECT Info.run_id, Info.run_type, Info.config_mask,Info.data_start_time,Info.db_start_time , grp_cmt.excluded_telescopes, grp_cmt.group_type, grp_cmt.group_id FROM tblRun_Info AS Info, tblRun_Group AS grp, tblRun_GroupComment AS grp_cmt WHERE ";
    // query a specific run number
    char c_query[flong_char_query];// has to be long if we ask for a long run list in VDBSourceInfo
    sprintf( c_query, " Info.run_id =  %u", arg_run );
    string run_query = c_query;
    // finish the query
    string query2 = " AND grp_cmt.group_id = grp.group_id AND grp.run_id = Info.run_id AND (Info.run_type = 'laser' OR Info.run_type = 'flasher') ORDER BY grp_cmt.excluded_telescopes ) AS big_table GROUP BY run_id,excluded_telescopes ; ";
    // built the full query
    string query = query1;
    query += run_query;
    query += query2;
    
    
    // if laser run belong to two different laser group with two different excluded_telescope, the smaller excluded_telescope is taken into account (so that the VOFFLINE DB is filled with for the maximum number of telescopes)
    //SELECT big_table.run_id, big_table.config_mask, big_table.excluded_telescopes big_table.data_start_time FROM
    //              (
    //               SELECT
    //                     Info.run_id,
    //                     Info.run_type,
    //                     Info.config_mask,
    //                     Info.data_start_time,
    //                     grp_cmt.excluded_telescopes,
    //                     grp_cmt.group_type,
    //                     grp_cmt.group_id
    //               FROM
    //                     tblRun_Info         AS Info,
    //                     tblRun_Group        AS grp,
    //                     tblRun_GroupComment AS grp_cmt
    //               WHERE
    //2016. added      Info.run_id = arg_run
    //2016. removed                     grp_cmt.group_type = 'laser'                               // because not sure the laser run belong to 'laser' group type if not official laser run
    //                 AND grp_cmt.group_id = grp.group_id                                         // need that or else too much lines
    //                 AND grp.run_id = Info.run_id                                                // need that or else too much lines
    //                 AND (Info.run_type = 'laser' OR Info.run_type = 'flasher')
    //
    //               ORDER BY
    //                     grp_cmt.excluded_telescopes
    //               )
    //         AS big_table
    //
    //         GROUP BY big_table.run_id;
    
    return query;
}

//---------------------------------------------------------------------------
//--2016. WriteQuery_to_get_one_LaserRun_fromVOFFLINE_DB()
//---------------------------------------------------------------------------
string WriteQuery_to_get_one_LaserRun_fromVOFFLINE_DB( unsigned int arg_run, unsigned int Tel_num )
{

    //begining of the querry
    string query = "SELECT big_table.run_id, big_table.code_version, big_table.update_time  FROM ( SELECT run_id, code_version, update_time  FROM tblEventDisplay_Analysis_Calibration_Flasher WHERE ";
    // select the telescope
    char c_query[flong_char_query];// has to be long if we ask for a long run list in VDBSourceInfo
    sprintf( c_query, "telescope = %u AND ", Tel_num );
    string tel_query = c_query;
    query += tel_query;
    sprintf( c_query, " run_id = %u", arg_run );
    tel_query = c_query;
    query += tel_query;
    
    // additional selection if requested in the command line
    if( fVOFFLINE_time_start_query_bool )
    {
        char c_time_query[flong_char_query];
        sprintf( c_time_query, "update_time > %s", fVOFFLINE_time_start_query.c_str() );
        string time_query = c_time_query;
        query += time_query;
        std::cout << "VOFFLINE Data Base: reading tblEventDisplay_Analysis_Calibration_Flasher for time > " << fVOFFLINE_time_start_query << std::endl;
    }
    if( fVOFFLINE_version_query_bool )
    {
        char c_version_query[flong_char_query];
        sprintf( c_version_query, "code_version = %d", fVOFFLINE_version_query );
        string version_query = c_version_query;
        query += version_query;
        std::cout << "VOFFLINE Data Base: reading tblEventDisplay_Analysis_Calibration_Flasher for version = " << fVOFFLINE_version_query << std::endl;
    }
    //end of the querry
    query += " ORDER BY update_time) as big_table GROUP BY big_table.run_id;";
    
    
    //SELECT big_table.run_id, big_table.code_version, big_table.update_time  FROM
    //             (
    //             SELECT run_id, code_version, update_time  FROM tblEventDisplay_Analysis_Calibration_Flasher
    //
    //            WHERE
    //                  telescope = %u
    //
    //                  run_id = %u
    //
    //                  update_time > %s
    //
    //                  code_version = %s
    //
    //            ORDER BY update_time
    //            )
    //            as big_table
    //
    //            GROUP BY big_table.run_id;
    
    
    
    return query;
}

//---------------------------------------------------------------------------
//-- get_date_from_tblRun_Info_data_start_time
//-- convert the date from data_start_time in table tblRun_Info in VERITAS DB
//-- to format year month day: yyyymmdd
//-- same way it is done in VExposure
//-- the formal has to be the one expected for download, because it give the name for the output directory for the downloaded file
//---------------------------------------------------------------------------
int get_date_from_tblRun_Info_data_start_time( string iTemp )
{

    int fDataStartTime = -1111;
    
    if( iTemp.size() > 8 )
    {
        fDataStartTime = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
    }
    
    return fDataStartTime;
    
}


//=== tools functions ==========================

void Set_the_time()
{
    //--- retrieve today's date
    time_t rawtime;
    struct tm* timeinfo;
    time( &rawtime );
    timeinfo = localtime( &rawtime );
    char today_now[1000]  ;
    sprintf( today_now, "%04d%02d%02d_%02dh%02dm%02ds", 1900 + timeinfo->tm_year, 1 + timeinfo->tm_mon, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec );
    fNOW = today_now;
    char today_DB[1000]  ;
    sprintf( today_DB, "%04d-%02d-%02d %02d:%02d:%02d", 1900 + timeinfo->tm_year, 1 + timeinfo->tm_mon, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec );
    std::cout << " today DB = " << today_DB << std::endl;
    fNOW_DB = today_DB;
    
    return;
    
}

//------------------------------------------------------------------------------------------
//-- bool read_calib_file
//-- put values in the caib file TTree in the tree given vector
//------------------------------------------------------------------------------------------
bool read_calib_file( TString file_root_name, vector < double >& Vchannel, vector < double >& Vmean, vector < double >& Vvar, bool bgain )
{

    bool reading_successfull = false;
    //-- name of the TTree and the leaf read
    TString tree_name =  "t";
    TString want_mean = "";
    if( bgain )
    {
        if( fcurrent_run < 63416 && fold_calib )
        {
            tree_name +=  "Gains_";
        }
        else
        {
            tree_name +=  "gain_";
        }
        tree_name += fcurrent_tel;
        want_mean = "gain";
    }
    else
    {
        if( fcurrent_run < 63416 && fold_calib )
        {
            tree_name +=  "Toffs_";
        }
        else
        {
            tree_name +=  "toff_";
        }
        tree_name += fcurrent_tel;
        want_mean = "toff";
    }
    TString want_var = want_mean;
    want_var += "var";
    
    
    TString want = "channel:";
    want += want_mean;
    want += ":";
    want += want_var;
    TString condition = "";
    
    //open root file
    TFile* root_file = new TFile( file_root_name );
    if( root_file->IsOpen() )
    {
        reading_successfull = fill_3V_from_TTree( root_file, tree_name, want, condition, Vchannel, Vmean, Vvar );
    }
    else
    {
        printf( "ERROR: File %s cannot be openned \n", file_root_name.Data() );
    }
    root_file->Close();
    
    return reading_successfull;
}


//===== test calibration
//---------------------------------------------------------------------------
//-- void check_run_calib
// for one run, on a list of good telescope ListTel (mask and not excluded), check the gain and toff root file
// create the list of telescope with correct calib output
//---------------------------------------------------------------------------
void check_run_calib( vector< int >  ListTel , vector< int >& list_of_valid_tel )
{

    // write introduction for the run
    printf( "-------------------------------------------------------------------------------- \n" );
    printf( "--- RUN %d      has %d telescope(s) to be tested \n", ( int ) fcurrent_run, ( int ) ListTel.size() );
    
    TString ltel_not_good = "";
    
    //-- Loop on the telescope
    for( unsigned int i = 0; i < ListTel.size() ; i++ )
    {
    
    
        fcurrent_tel = ListTel[i];
        printf( "--- Tel %d \n", fcurrent_tel );
        //-- get the name of the calibration root file to be checked
        TString file_gain_root = get_name_calib_file( fCalib_dir, "gain.root", ListTel[i], fcurrent_run );
        TString file_toff_root = get_name_calib_file( fCalib_dir, "toff.root", ListTel[i], fcurrent_run );
        //-- test the calib file for the run and a given telescope
        if( test_gain_toff( file_gain_root, file_toff_root ) )
        {
            list_of_valid_tel.push_back( ListTel[i] );
        }
        else
        {
            ltel_not_good += ListTel[i];
        }
        
        
        
    }
    //-- write conclu for the run
    printf( "--- RUN %d      has %d/%d telescope(s) with valid calibration: \n", fcurrent_run, ( int ) list_of_valid_tel.size(), ( int ) ListTel.size() );
    //  std::cout<<"--- 1 "<<std::endl;
    //   std::cout<<"--- 2 "<<std::endl;
    printf( "--- (Display: $EVNDISPSYS/scripts/VTS/VTS.EVNDISP.display %s %d ) \n", ltel_not_good.Data(), fcurrent_run );
    // std::cout<<"--- 3 "<<std::endl;
    
    if( list_of_valid_tel.size() > 0 )
    {
        for( unsigned int i = 0; i < list_of_valid_tel.size(); i++ )
        {
            printf( "---------  - %d %d \n", fcurrent_run, list_of_valid_tel[i] );
        }
    }
    if( list_of_valid_tel.size() < ListTel.size() )
    {
        printf( "--- Should run the Display: \n $EVNDISPSYS/scripts/VTS/VTS.EVNDISP.display %s %s \n (May have to find replace this laser run) \n", ltel_not_good.Data(), get_downloaded_laser_path( fcurrent_run, fcurrent_date ).Data() );
        
    }
    
    
    return;
}


//---------------------------------------------------------------------------
//-- TH1F* get_h
//-- return the histogram from the Drawing of the TTree
//---------------------------------------------------------------------------
TH1F* get_h_from_TTree( TFile* file_root, TString string_arbre, TString want, TString condition, TString histo_name, TString histo_bin, bool normalised )
{

    //    std::cout<<"get_h_from_TTree "<<std::endl;
    
    TTree* arbre_gain = ( TTree* ) file_root->Get( string_arbre.Data() )  ;
    
    // initial want, for example gain
    want += ">>hsqrt";
    want += histo_bin; // for exemple (500,0,10)
    
    //    std::cout<<"want "<<want<<std::endl;
    //   std::cout<<"condition "<<condition<<std::endl;
    arbre_gain->Draw( want, condition, "goff" );
    
    TH1F* h = ( TH1F* )gDirectory->Get( "hsqrt" );
    
    TH1F* hh = ( TH1F* ) h->Clone( histo_name );
    
    if( normalised )
    {
        hh->Scale( 1. / hh->GetEntries() );
    }
    return hh;
    
}

//---------------------------------------------------------------------------
//-- void fill_3V_from_TTree
//-- read the tree in the root file
//-- draw what you want (WARNING, you must want 3 things) under the given condition
//-- copy the drawn values in the 3 given vector
//---------------------------------------------------------------------------
bool fill_3V_from_TTree( TFile* file_root, TString string_arbre, TString want, TString condition, vector < double >& Varg1, vector < double >& Varg2, vector < double >& Varg3 )
{

    bool good_filling = false;
    TTree* arbre = ( TTree* ) file_root->Get( string_arbre.Data() );
    
    Long64_t draw_res = arbre->Draw( want, condition, "goff" );
    
    if( draw_res > 0 )
    {
        vector < double > Vec1( arbre->GetV1(), arbre->GetV1() + arbre->GetSelectedRows() ) ;
        vector < double > Vec2( arbre->GetV2(), arbre->GetV2() + arbre->GetSelectedRows() ) ;
        vector < double > Vec3( arbre->GetV3(), arbre->GetV3() + arbre->GetSelectedRows() ) ;
        
        Varg1 = Vec1;
        Varg2 = Vec2;
        Varg3 = Vec3;
        good_filling = true;
        
    }
    else if( draw_res == 0 )
    {
        std::cout << "fill_3V_from_TTree ERROR: Vector could not be filled because no row were selected" << std::endl;
        std::cout << "                          Check the TTree " << string_arbre << ", what you want:  " << want << " and your condition: " << condition << std::endl;
        
    }
    else
    {
        std::cout << "fill_3V_from_TTree ERROR: Erreur drawing the TTree " << string_arbre << std::endl;
    }
    
    
    return good_filling;
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
bool does_file_exist( TString file_path_name )
{
    bool file_exist = false;
    ifstream file;
    file.open( file_path_name );
    if( file.is_open() )
    {
        file_exist = true;
    }
    else
    {
        std::cout << "WARNING : File  " << file_path_name << " does not exist" << std::endl;
    }
    file.close();
    return file_exist;
}

bool is_file_empty( TString file_path_name )
{
    bool file_empty = true;
    ifstream file;
    file.open( file_path_name );
    if( file )
    {
        if( file.get() > 0 )
        {
            file_empty = false;
        }
        else
        {
            std::cout << "File " << file_path_name << " is empty " << std::endl; // not an error, if every thing is analysed, the file should be empty
        }
        
    }
    else
    {
        std::cout << "ERROR: File  " << file_path_name << " does not exist" << std::endl;
    }
    file.close();
    return file_empty;
}
//---------------------------------------------------------------------------
//-- Completion file name for download
//---------------------------------------------------------------------------
TString get_name_for_download( TString name )
{
    TString name_for_download = name;
    name_for_download += "_for_download";
    return name_for_download;
}
//---------------------------------------------------------------------------
//-- Completion file name for only run number
//---------------------------------------------------------------------------
TString get_name_run_to_be_analysed( TString name )
{
    TString name_only_run_number = name;
    name_only_run_number += "_ana_todo";
    return name_only_run_number;
}
//---------------------------------------------------------------------------
//--
//---------------------------------------------------------------------------
TString get_name_run_already_analysed( TString name )
{
    TString name_only_run_number = name;
    name_only_run_number += "_ana_done";
    return name_only_run_number;
}
//---------------------------------------------------------------------------
//--
//---------------------------------------------------------------------------
TString get_name_run_to_be_writen( TString name, bool long_table )
{
    TString name_wDB = name;
    if( long_table )
    {
        name_wDB += "_tblEventDisplay_Analysis_Calibration_Flasher";
    }
    else
    {
        name_wDB += "_tblEventDisplay_Analysis_Status";
    }
    return name_wDB;
}
//---------------------------------------------------------------------------
//--
//---------------------------------------------------------------------------
TString get_name_run_problem( TString name )
{
    TString name_pbCalib = name;
    name_pbCalib += "_pbCalib";
    return name_pbCalib;
}

//---------------------------------------------------------------------------
//-- built the name
//-- same procedure as in VExposure::downloadRunList()
//---------------------------------------------------------------------------
TString get_downloaded_laser_path( unsigned int RunNumber, unsigned int Date )
{

    char filename[800];
    char* ENVIR_VAR;
    ENVIR_VAR = getenv( "VERITAS_DATA_DIR" );
    sprintf( filename, "%s/data/d%d/%d.cvbf", ENVIR_VAR, Date, RunNumber );
    
    TString path =  filename;
    return path;
}

//---------------------------------------------------------------------------
//-- built the name
//---------------------------------------------------------------------------
TString get_name_calib_file( TString dir, TString suffixe, int tel_num, unsigned int RunNumber )
{

    TString file_name = dir;
    file_name += "/Tel_";
    file_name += tel_num;
    file_name += "/";
    file_name += RunNumber;
    file_name += ".";
    file_name += suffixe;
    
    return file_name;
}



//---------------------------------------------------------------------------
//-- get_list_of_telescope                                                  --
//-- return a vector with the number of the telescope participating to the run and not excluded from it
//---------------------------------------------------------------------------
vector<int> get_list_of_telescope( unsigned int config_mask, unsigned int excluded_tel )
{

    vector<int> list_of_telescope;
    
    if( Tel_code_not_valid( config_mask ) )
    {
        std::cout << "WARNING: config mask (" << config_mask << ") for run " << fcurrent_run << " is NOT VALID" << std::endl;
        return list_of_telescope;
    }
    if( Tel_code_not_valid( excluded_tel ) )
    {
        std::cout << "WARNING: excluded tel (" << excluded_tel << ") for run " << fcurrent_run << " is NOT VALID" << std::endl;
        return list_of_telescope;
    }
    
    for( int i = 0 ; i < fmax_number_tel ; i++ )
    {
        int Tel = i + 1;
        if( Incoherence_between_ConfigMask_and_excludedTel( Tel, config_mask, excluded_tel ) )
        {
            std::cout << "WARNING: Incoherence between  config mask (" << config_mask << ") and excluded tel (" << excluded_tel << ") for RUN " << fcurrent_run << std::endl;
            //return list_of_telescope;
        }
    }
    
    for( int i = 0 ; i < fmax_number_tel ; i++ )
    {
    
        int Tel = i + 1;
        if( Tel_is_in_mask( Tel, config_mask ) && !Tel_is_excluded( Tel, config_mask, excluded_tel ) )
        {
            list_of_telescope.push_back( Tel );
        }
    }
    
    return list_of_telescope;
    
}

//---------------------------------------------------------------------------
//-- get_list_of_telescope                                                 --
//-- return a vector with the telescope in the given mask                  --
//---------------------------------------------------------------------------
vector<int> get_list_of_telescope( unsigned int mask )
{

    vector<int> list_of_telescope;
    
    if( Tel_code_not_valid( mask ) )
    {
        std::cout << "WARNING: mask (" << mask << ") for run " << fcurrent_run << " is NOT VALID" << std::endl;
        return list_of_telescope;
    }
    
    for( int i = 0 ; i < fmax_number_tel ; i++ )
    {
        int Tel = i + 1;
        if( Tel_is_in_mask( Tel, mask ) )
        {
            list_of_telescope.push_back( Tel );
        }
    }
    
    return list_of_telescope;
    
}


//----------------------------------------------
//-- Tel_is_in_mask                           --
//-- return true if Tel is in mask            --
//----------------------------------------------
bool Tel_is_in_mask( Int_t Tel, Int_t config_mask )
{
    if( Tel_code_not_valid( config_mask ) )
    {
        std::cout << "WARNING: config_mask  " << config_mask << std::endl;
        return false;
    }
    Int_t Tel_test = Tel - 1;
    bitset< 8 > ibit_conf( config_mask );
    
    
    if( ibit_conf.test( Tel_test ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}
//----------------------------------------------
//-- Tel_is_excluded                          --
//-- return true if the telescope is excluded --
//----------------------------------------------
bool Tel_is_excluded( Int_t Tel, Int_t config_mask, Int_t excluded_tel )
{
    if( Tel_code_not_valid( config_mask ) || Tel_code_not_valid( excluded_tel ) )
    {
        std::cout << "WARNING.... " << std::endl;
        return true;
    }
    //if(Incoherence_between_ConfigMask_and_excludedTel(Tel,config_mask,excluded_tel)){ std::cout<<"WARNING.... "<<std::endl; return true;}
    
    Int_t Tel_test = Tel - 1;
    bitset< 8 > ibit_excl( excluded_tel );
    bitset< 8 > ibit_conf( config_mask );
    
    if( ibit_excl.test( Tel_test ) && ibit_conf.test( Tel_test ) )
    {
        return true;
    }
    else if( ibit_excl.test( Tel_test ) && !ibit_conf.test( Tel_test ) )
    {
        std::cout << "warning: Tel is excluded BUT REDONDANT INFORMATION because was not in configmask" << std::endl;
        return true;
    }
    else
    {
        return false;
    }
}
//---------------------------------------------------
//-- Incoherence_between_ConfigMask_and_excludedTel--
//-- return true if there is an incoherence        --
//---------------------------------------------------
bool Incoherence_between_ConfigMask_and_excludedTel( Int_t Tel, Int_t config_mask, Int_t excluded_tel )
{
    if( Tel_code_not_valid( config_mask ) || Tel_code_not_valid( excluded_tel ) )
    {
        std::cout << "WARNING... " << std::endl;
        return true;
    }
    Int_t Tel_test = Tel - 1;
    bitset< 4 > ibit_excl( excluded_tel );
    bitset< 4 > ibit_conf( config_mask );
    
    if( ibit_excl.test( Tel_test ) && !ibit_conf.test( Tel_test ) )
    {
        std::cout << "Tel " << Tel << " is excluded (" << excluded_tel << ") but is not in config_mask (" << config_mask << ")" << std::endl;
        return true;
    }
    else
    {
        //std::cout<<"No incoherence "<<std::endl;
        return false;
    }
}


//----------------------------------------------
//-- Tel_code_not_valid                       --
//-- return true if the tel code is not valid --
//----------------------------------------------
bool Tel_code_not_valid( Int_t Tel_code )
{
    Int_t Tel_max = fmax_number_tel;
    bitset< 8 > ibit( Tel_code );
    
    if( ibit.test( Tel_max ) )
    {
        std::cout << " WARNING: WRONG Tel_Code = " << Tel_code << " not possible if max number of Telescope = " << Tel_max << std::endl;
        return true;
    }
    else
    {
        return false;
    }
}

//===============================================


