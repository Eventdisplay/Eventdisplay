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
bool fwriteDB = false;
bool fget_list_new_laser_run = false;
bool fsave_list_laser_run = false;
bool fdownload_new_laser_run = false;
bool fprepareEVNDanalysis = false;
bool fdoEVNDanalysis = false;
bool ftestEVNDcalib = false;
bool freANA_ALL_RUN = false;
bool fold_calib = false;

bool fwrite_only = false;
TString fwrite_only_run_list = "";

bool finput_list_run = false;
TString flist_run = "";

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
//--- Tool-global-parameters
int fcurrent_run = -10;
int fcurrent_tel = -10;
int fcurrent_date = -10;

TString fNOW = "demain"; // initialized once, when doing get_new_laser_run_list()
TString fNOW_DB = "demain aussi"; // initialized once, when doing get_new_laser_run_list()

//==============================================================================================================
// FUNCTION DECLARATION
//==============================================================================================================
//-- Functions in the main
int parseOptions( int argc, char* argv[] );
TString get_new_laser_run_list();
void save_list_currently_in_DB( TString NOW, vector< unsigned int > VERITAS_DB_LaserRunNumber, vector< unsigned int > VERITAS_DB_LaserConfigMask, vector< unsigned int > VERITAS_DB_LaserExclTel, vector< unsigned int > VERITAS_DB_LaserDate, vector< vector < unsigned int > > VOFFLINE_DB_LaserRunNumber_Tel, vector< vector < string > > VOFFLINE_DB_Laserdate_Tel, vector< vector < string > > VOFFLINE_DB_Laserversion_Tel );
void download_new_laser_run( TString new_laser_run_list_name );
void prepareEVNDanalysis( TString new_laser_run_list_name );
void doEVNDanalysis( TString new_laser_run_list_name );
void prepare_CalibVOFF_writing( TString new_laser_run_list_name );
void write_calib_DB( TString new_laser_run_list_name );
//-------------------------

//-- Active function
//------------------------------- used in get_new_laser_run_list():
unsigned long Check_telmissing_from_VOFFDB( unsigned int VERITAS_DB_LaserRunNumber_i_run, unsigned int VERITAS_DB_LaserConfigMask_i_run, unsigned int VERITAS_DB_LaserExclTel_i_run, vector< vector < unsigned int > > VOFFLINE_DB_LaserRunNumber_Tel );
void check_run_calib( unsigned int RunNumber, int  mask, int excluded_tel_mask );
void check_run_calib( vector< int >  ListTel , vector< int >& list_of_valid_tel, FILE*& pbFile );
bool test_gain_toff( TString file_root_name_gain, TString file_root_name_toff , FILE*& pbFile );
bool construct_VOFFLINE_calibration_writing_file( FILE*& wDB_file );
//-------------------------------------------------------------------


//-- Accessing the DB
void read_laserRUN_fromVERITAS_DB( vector< unsigned int >& VERITAS_DB_LaserRunNumber, vector< unsigned int >& VERITAS_DB_LaserConfigMask, vector< unsigned int >& VERITAS_DB_LaserExclTel, vector< unsigned int >& VERITAS_DB_LaserDate );
void read_laserRUN_fromVOFFLINE_DB( vector < vector < unsigned int > >& VOFFLINE_DB_LaserRunNumber_Telnum, vector < vector < string > >& VOFFLINE_DB_LaserDate_Telnum, vector < vector < string > >& VOFFLINE_DB_LaserVersion_Telnum, unsigned int Tel_num );
string WriteQuery_to_getLaserRun_fromVERITAS_DB();
string WriteQuery_to_getLaserRun_fromVOFFLINE_DB( unsigned int Tel_num );

//=== tools functions ==========================
void Set_the_time();
//-- reading
void readfile_5_unsigned_int( TString file_name, vector< unsigned int >& VERITAS_DB_LaserRunNumber, vector< unsigned int >& VERITAS_DB_LaserConfigMask, vector< unsigned int >& VERITAS_DB_LaserExclTel, vector< unsigned int >& VERITAS_DB_LaserMiss, vector< unsigned int >& VERITAS_DB_LaserDate );
bool read_calib_file( TString file_root_name, vector < double >& Vchannel, vector < double >& Vmean, vector < double >& Vvar, bool bgain );
void create_list_for_download( TString file_name );
//-- telescope mask related
vector<int> get_list_of_telescope( unsigned int tel_mask, unsigned int excluded_tel );
vector<int> get_list_of_telescope( unsigned int missing_tel_mask );
bool Tel_is_in_mask( Int_t Tel, Int_t config_mask );
bool Tel_is_excluded( Int_t Tel, Int_t config_mask, Int_t excluded_tel );
bool Incoherence_between_ConfigMask_and_excludedTel( Int_t Tel, Int_t config_mask, Int_t excluded_tel );
bool Tel_code_not_valid( Int_t Tel_code );
//-- reading DB time (could be integrated to VExposure)
int get_date_from_tblRun_Info_data_start_time( string iTemp );

//-- File existence and filing
bool does_file_exist( TString file_path_name );
bool is_file_empty( TString file_path_name );

//-- Reading TTree
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
    //int main() {
    
    //--- indicate if the function has succeeded
    bool_all_is_downloaded = false;
    bool_EVNDanalysis_prepared = false;
    bool_EVNDanalysis_done = false;
    bool_EVNDcalib_good = false;
    
    parseOptions( argc, argv );
    if( ( fdownload_new_laser_run || fprepareEVNDanalysis ) && ( !fget_list_new_laser_run && !finput_list_run ) )
    {
        std::cout << "Need a run list" << std::endl;
        std::cout << "options:" << std::endl;
        std::cout << "-g (will get it automatically)" << std::endl;
        std::cout << "or" << std::endl;
        std::cout << "-l (provide your own run list with format: RunNumber MaskVERITAS_DB ExclVERITAS_DB Miss_fromVOFF_DB Date)" << std::endl;
        
        exit( 0 );
    }
    
    if( freANA_ALL_RUN )
    {
        std::cout << "WARNING: all RUN will be analysed event if the calibration already exist !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
    if( ftestEVNDcalib )
    {
        std::cout << "REMINDER: Calibration checking, HARD coded name TTree read are: " << fcalib_gain_tree_name << " and " << fcalib_toff_tree_name << std::endl;
    }
    // they are the same for the all run list, which have to be homogeneous in this respect
    // name could be change at the parse ooption level (but for this would be good to have option instroduced with more than one letter).
    
    std::cout << "======================================================================" << std::endl;
    VGlobalRunParameter* blah = new VGlobalRunParameter() ;
    fServer = blah->getDBServer();
    cout << "Using server " << fServer << endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "ATTENTION: make sure the Calib_dir is correctly defined (option -b). This is were the output of the .gain .ped ,toff will be " << std::endl;
    std::cout << "in this directory there should be one directory per telescope (Tel_1, Tel_2, ...) " << std::endl;
    std::cout << "Calib_dir  " << fCalib_dir << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    Set_the_time();
    //------------ getting the list of new laser run (present in the veritas DB but not in the VOFFLINE)
    TString new_laser_run_list_name = "";
    if( fget_list_new_laser_run )
    {
        // -g
        std::cout << "Getting list of ALL laser flasher run from VERITAS DB " << std::endl;
        new_laser_run_list_name = get_new_laser_run_list(); //working
    }
    else if( finput_list_run )
    {
        // -l
        new_laser_run_list_name = flist_run;
        std::cout << "manual input run list " << std::endl;
        std::cout << "BEWARE: " << std::endl;
        std::cout << "Structure has to be: " << std::endl;
        std::cout << "RunNumber MaskVERITAS_DB ExclVERITAS_DB Miss_fromVOFF_DB Date " << std::endl;
        // to do: add a function testing if the format of the run list is valid
    }
    else
    {
    
        std::cout << "ERROR: NEED a run list " << std::endl;
        std::cout << "use -g option or -l name_of_the_run_list " << std::endl;
        exit( 0 );
    }
    //------------ downloading the laser run
    // -d (need -g or -l)
    if( fdownload_new_laser_run )
    {
    
        if( finput_list_run )
        {
            create_list_for_download( new_laser_run_list_name );
        }
        download_new_laser_run( new_laser_run_list_name ); // working, needs to be on transfer
    }
    //------------ create the list of run to be analysed
    // -p (need -g or -l, not -d  we can prepare even if download is not complete, but -d should be done one time before)
    if( fprepareEVNDanalysis )
    {
        prepareEVNDanalysis( new_laser_run_list_name ); // working
    }
    //------------ launch the analysis
    // -e (if -e then -p should be automatic, not d because not on the same machine and if one run has problem to be transfered, then we still need to move on)
    if( fdoEVNDanalysis && bool_EVNDanalysis_prepared )
    {
        doEVNDanalysis( new_laser_run_list_name ); // working (seems to work even if we are on transfer)
    }
    //------------ test the output of the analysis and write the file to be used to fill the VOFFLine DB
    //-c (if -t, then -e and should be automatic and all the option implied but -e)
    if( ftestEVNDcalib && bool_EVNDanalysis_done )
    {
        prepare_CalibVOFF_writing( new_laser_run_list_name ); // done: to be tested (seems ok, but problem with TTree named changed during an EVNDISP update, have to rerun all runs)
    }
    //------------
    
    // -w
    if( fwriteDB && bool_EVNDcalib_good )
    {
        write_calib_DB( new_laser_run_list_name );
    }
    
    
    if( fwrite_only )
    {
        write_calib_DB( fwrite_only_run_list );
    }
    
    
}




//===============================================================================================
//==== MAIN FUNCTIONS
//===============================================================================================

//---------------------------------------------------------------------------
//-- get_new_laser_run_list                                                --
//-- create a file with list of new laser run                              --
//-- laser run present in VERITAS DB but not in VOFFLINE                   --
//---------------------------------------------------------------------------
TString get_new_laser_run_list()
{

    //--- Create File with the list of new laser run (to be updated in the VOFFLINE DB)
    
    //--- create the name of the file
    TString NEW_laser_list = "Laser_RunMaskExclMissDate_NEW_LIST_";
    NEW_laser_list += fNOW;
    //--- create the file
    FILE* NEW_laser_list_File;
    NEW_laser_list_File = fopen( NEW_laser_list.Data(), "w" );
    
    //--- Create file for downloading the laser run
    TString NEW_laser_list_for_download = get_name_for_download( NEW_laser_list );
    FILE* NEW_laser_list_File_for_download;
    if( fdownload_new_laser_run )
    {
        NEW_laser_list_File_for_download = fopen( NEW_laser_list_for_download.Data(), "w" );
    }
    
    //--- get the list of all laser run from VERITAS DB
    vector< unsigned int > VERITAS_DB_LaserRunNumber;
    vector< unsigned int > VERITAS_DB_LaserConfigMask;
    vector< unsigned int > VERITAS_DB_LaserExclTel;
    vector< unsigned int > VERITAS_DB_LaserDate;
    read_laserRUN_fromVERITAS_DB( VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VERITAS_DB_LaserDate );
    
    //--- get the list of laser run from VOFFSET (for a given analysis version)
    vector< vector < unsigned int > > VOFFLINE_DB_LaserRunNumber_Tel;
    vector< vector < string > > VOFFLINE_DB_Laserdate_Tel;
    vector< vector < string > > VOFFLINE_DB_Laserversion_Tel;
    //max number of telescope is hard coded here
    // and the oder in which the telescope are called is important
    read_laserRUN_fromVOFFLINE_DB( VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 1 );
    read_laserRUN_fromVOFFLINE_DB( VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 2 );
    read_laserRUN_fromVOFFLINE_DB( VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 3 );
    read_laserRUN_fromVOFFLINE_DB( VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel, 4 );
    
    //--- SAVING the list of run present in the DB
    if( fsave_list_laser_run )
    {
        save_list_currently_in_DB( fNOW, VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VERITAS_DB_LaserDate, VOFFLINE_DB_LaserRunNumber_Tel, VOFFLINE_DB_Laserdate_Tel, VOFFLINE_DB_Laserversion_Tel );
    }
    
    //--- compare the list of laser run from VERITAS DB with the one from VOFFLINE
    
    //--- loop on the laser run in VERITAS data base
    std::cout << "Checking wich laser/flasher run are missing from VOFFLINE DB " << std::endl;
    for( unsigned int i_run = 0 ; i_run < VERITAS_DB_LaserRunNumber.size(); i_run++ )
    {
        unsigned long mask_missing_tel = Check_telmissing_from_VOFFDB( VERITAS_DB_LaserRunNumber[i_run], VERITAS_DB_LaserConfigMask[i_run], VERITAS_DB_LaserExclTel[i_run], VOFFLINE_DB_LaserRunNumber_Tel );
        
        if( VERITAS_DB_LaserRunNumber[i_run] == 53139 )
        {
            printf( "** %d %2d %2d %2lu %d \n", VERITAS_DB_LaserRunNumber[i_run], VERITAS_DB_LaserConfigMask[i_run], VERITAS_DB_LaserExclTel[i_run], mask_missing_tel, VERITAS_DB_LaserDate[i_run] );
            
        }
        
        if( mask_missing_tel > 0 )
        {
            std::cout << " " << std::endl;		 //--- fill the NEW_laser_list_File with: run mask exclu_tel mask_with_missing_tel_from_VOFFLINE_DB
            printf( "*** %d %2d %2d %2lu %d \n", VERITAS_DB_LaserRunNumber[i_run], VERITAS_DB_LaserConfigMask[i_run], VERITAS_DB_LaserExclTel[i_run], mask_missing_tel, VERITAS_DB_LaserDate[i_run] );
            fprintf( NEW_laser_list_File, "%d %2d %2d %2lu %d \n", VERITAS_DB_LaserRunNumber[i_run], VERITAS_DB_LaserConfigMask[i_run], VERITAS_DB_LaserExclTel[i_run], mask_missing_tel, VERITAS_DB_LaserDate[i_run] );
            std::cout << " " << std::endl;		 //--- fill the NEW_laser_list_File with: run mask exclu_tel mask_with_missing_tel_from_VOFFLINE_DB
            if( fdownload_new_laser_run )
            {
                fprintf( NEW_laser_list_File_for_download, "%d %d \n", VERITAS_DB_LaserRunNumber[i_run], VERITAS_DB_LaserDate[i_run] );
            }
        }
    }
    
    //--- end loop on RUN in VERITAS DB
    //--- close the file
    fclose( NEW_laser_list_File );
    if( fdownload_new_laser_run )
    {
        fclose( NEW_laser_list_File_for_download );
    }
    return NEW_laser_list;
}



//---------------------------------------------------------------------------
//-- save_list_currently_in_DB                                             --
//---------------------------------------------------------------------------
void save_list_currently_in_DB( TString NOW, vector< unsigned int > VERITAS_DB_LaserRunNumber, vector< unsigned int > VERITAS_DB_LaserConfigMask, vector< unsigned int > VERITAS_DB_LaserExclTel, vector< unsigned int > VERITAS_DB_LaserDate, vector< vector < unsigned int > > VOFFLINE_DB_LaserRunNumber_Tel, vector< vector < string > > VOFFLINE_DB_Laserdate_Tel, vector< vector < string > > VOFFLINE_DB_Laserversion_Tel )
{

    //--- List of run from VERITAS DB
    TString All_laser_VERITAS_DB = "Laser_RunMaskExclDate_VERITAS_DB_";
    All_laser_VERITAS_DB += NOW ;
    FILE* All_laser_VERITAS_DB_File;
    All_laser_VERITAS_DB_File = fopen( All_laser_VERITAS_DB.Data(), "w" );
    for( unsigned int i_run = 0 ; i_run < VERITAS_DB_LaserRunNumber.size(); i_run++ )
    {
        fprintf( All_laser_VERITAS_DB_File, "%d %2d %2d %d \n", VERITAS_DB_LaserRunNumber[i_run], VERITAS_DB_LaserConfigMask[i_run], VERITAS_DB_LaserExclTel[i_run], VERITAS_DB_LaserDate[i_run] );
        
    }
    fclose( All_laser_VERITAS_DB_File );
    //--- list of run in the VOFFLINE DB
    for( unsigned i_tel = 0; i_tel < 4; i_tel++ ) // attention hard coded value for the number of telescope
    {
    
        int telescope =  i_tel + 1;
        TString All_laser_VOFFLINE_DB = "Laser_RunDateVers_VOFFLINE_DB_Tel";
        All_laser_VOFFLINE_DB += telescope;
        All_laser_VOFFLINE_DB += "_";
        if( fVOFFLINE_time_start_query_bool )
        {
            All_laser_VOFFLINE_DB += "uptime_" ;
            All_laser_VOFFLINE_DB += fVOFFLINE_time_start_query;
            All_laser_VOFFLINE_DB += "_" ;
        }
        if( fVOFFLINE_version_query_bool )
        {
            All_laser_VOFFLINE_DB += "Version_" ;
            All_laser_VOFFLINE_DB += fVOFFLINE_version_query;
            All_laser_VOFFLINE_DB += "_" ;
        }
        All_laser_VOFFLINE_DB += NOW ;
        FILE* All_laser_VOFFLINE_DB_File;
        All_laser_VOFFLINE_DB_File = fopen( All_laser_VOFFLINE_DB.Data(), "w" );
        //--- loop on the run in VOFFLINE DB
        for( unsigned int i_run_voff = 0; i_run_voff < VOFFLINE_DB_LaserRunNumber_Tel[i_tel].size(); i_run_voff++ )
        {
            fprintf( All_laser_VOFFLINE_DB_File, "%d %s %s\n", VOFFLINE_DB_LaserRunNumber_Tel[i_tel][i_run_voff], VOFFLINE_DB_Laserdate_Tel[i_tel][i_run_voff].c_str(), VOFFLINE_DB_Laserversion_Tel[i_tel][i_run_voff].c_str() );
        }
        fclose( All_laser_VOFFLINE_DB_File );
    }
    std::cout << " saving the list of all laser run from VERITAS DB (no repetition, smaller excluded telescope kept) " << std::endl;
    std::cout << " saving the list of all laser run from VOFFLINE DB (starting at a giving time and for a given version if requested) " << std::endl;
    return;
}

//---------------------------------------------------------------------------
//-- download_new_laser_run                                                --
//---------------------------------------------------------------------------
void download_new_laser_run( TString new_laser_run_list_name )
{


    VExposure a;
    a.setSelectLaser( 1 );
    a.readLaserRunDateListFromFile( ( string ) get_name_for_download( new_laser_run_list_name ) ); // run date yearmonthday
    a.downloadRunList();
    
    bool_all_is_downloaded = true;
    
    return;
}

//---------------------------------------------------------------------------
//-- prepareEVNDanalysis
//-- input: list run in VERITAS DB but not in VOFF
//-- output: list of run to be analysed
//-- output: list of run already analysed (to be used during the testing of the output of the analysis)
//-- checking the existence of the output of EVNDISP in output dir (which is hard coded)
//---------------------------------------------------------------------------
void prepareEVNDanalysis( TString new_laser_run_list_name )
{

    //-- Create File to list the run that have to be analysed
    TString  new_laser_run_list_name_for_EVNDISP_analyse_laser = get_name_run_to_be_analysed( new_laser_run_list_name );
    FILE* NEW_laser_list_File_for_EVNDISP_analyse_laser;
    NEW_laser_list_File_for_EVNDISP_analyse_laser = fopen( new_laser_run_list_name_for_EVNDISP_analyse_laser.Data(), "w" );
    //-- Create File to list the run already analysed
    TString  new_laser_run_list_name_EVNDISP_analyse_laser_DONE = get_name_run_already_analysed( new_laser_run_list_name );
    FILE* NEW_laser_list_File_EVNDISP_analyse_laser_DONE;
    NEW_laser_list_File_EVNDISP_analyse_laser_DONE = fopen( new_laser_run_list_name_EVNDISP_analyse_laser_DONE.Data(), "w" );
    
    // open the file listing the run missing from VOFFLINE DB
    ifstream file_new_laser_run_list_name;
    file_new_laser_run_list_name.open( new_laser_run_list_name );
    if( !file_new_laser_run_list_name )
    {
        std::cout << "ERROR: imposible to open file " << new_laser_run_list_name << std::endl;
        fclose( NEW_laser_list_File_for_EVNDISP_analyse_laser );
        fclose( NEW_laser_list_File_EVNDISP_analyse_laser_DONE );
        return;
    }
    
    
    unsigned int counter = 0;
    //--- read File run list missing from VOFF DB and store information
    vector< unsigned int > VERITAS_DB_LaserRunNumber;
    vector< unsigned int > VERITAS_DB_LaserConfigMask;
    vector< unsigned int > VERITAS_DB_LaserExclTel;
    vector< unsigned int > VERITAS_DB_LaserMiss;
    vector< unsigned int > VERITAS_DB_LaserDate;
    readfile_5_unsigned_int( new_laser_run_list_name, VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VERITAS_DB_LaserMiss, VERITAS_DB_LaserDate );
    
    //-- loop on the run missing from VOFFLINE DB
    for( unsigned int j = 0 ; j < VERITAS_DB_LaserRunNumber.size(); j++ )
    {
    
        // test if the file is on disk, other wise skip
        if( !does_file_exist( get_downloaded_laser_path( VERITAS_DB_LaserRunNumber[j], VERITAS_DB_LaserDate[j] ) ) )
        {
            std::cout << "***** Run " << VERITAS_DB_LaserRunNumber[j] << " cannot be launched" << std::endl;
            continue;
        }
        // if it's not, it is not an error, problem with the download, but we still want to prepare the analysis for the remaining run
        
        fcurrent_run = VERITAS_DB_LaserRunNumber[j];
        bool need_ana = false;
        
        vector<int> VList_of_telescope = get_list_of_telescope( VERITAS_DB_LaserMiss[j] );
        //-- loop on the telescope missing from VOFFLINE DB
        for( unsigned int i = 0 ; i < VList_of_telescope.size(); i++ )
        {
        
            //-- expected files:
            TString file_ped  = get_name_calib_file( fCalib_dir, "ped" , VList_of_telescope[i], fcurrent_run );
            TString file_gain = get_name_calib_file( fCalib_dir, "gain", VList_of_telescope[i], fcurrent_run );
            TString file_toff = get_name_calib_file( fCalib_dir, "toff", VList_of_telescope[i], fcurrent_run );
            
            //-- file must exist
            if( !does_file_exist( file_ped ) || !does_file_exist( file_gain ) || !does_file_exist( file_toff ) || freANA_ALL_RUN )
            {
                need_ana = true;
                counter++;
                break;
            }
            //-- note, emptiness of root file, a matter of testing if the analysis was alright, not to be done here
            
        }
        
        // just create a list with run missing something in the analysis.  Not for telescope wise analyse... (could do that...)
        if( need_ana )
        {
            fprintf( NEW_laser_list_File_for_EVNDISP_analyse_laser, "%d\n", fcurrent_run );
        }
        else
        {
            std::cout << fcurrent_run << " ana done" << std::endl;
            fprintf( NEW_laser_list_File_EVNDISP_analyse_laser_DONE, "%d %2d %2d %2u %d \n", VERITAS_DB_LaserRunNumber[j], VERITAS_DB_LaserConfigMask[j], VERITAS_DB_LaserExclTel[j], VERITAS_DB_LaserMiss[j], VERITAS_DB_LaserDate[j] );
        }
    }
    file_new_laser_run_list_name.close();
    fclose( NEW_laser_list_File_for_EVNDISP_analyse_laser );
    fclose( NEW_laser_list_File_EVNDISP_analyse_laser_DONE );
    
    
    bool_EVNDanalysis_prepared = true;
    
    std::cout << "****** " << counter << " run to be analysed " << std::endl;
    
    return;
}
//---------------------------------------------------------------------------
//-- Launch EVNDISP analysis on a list of laser run
//-- using script
//---------------------------------------------------------------------------
void doEVNDanalysis( TString new_laser_run_list_name )
{

    //-- get the name of the file to be passed to the script launching EVNDISP laser analysis
    TString run_list = get_name_run_to_be_analysed( new_laser_run_list_name );
    // test if the file exist
    if( !does_file_exist( run_list ) )
    {
        std::cout << "Error: File " << run_list << " does not exist" << std::endl;
        return;
    }
    
    //-- launch the analysis if run_list is not empty
    if( !is_file_empty( run_list ) )
    {
        char launching_EVNDISP_string_string[800];
        sprintf( launching_EVNDISP_string_string, "$EVNDISPSYS/scripts/VTS/VTS.EVNDISP.sub_analyse_laser.sh %s", run_list.Data() );
        std::cout << "COMMAND " << launching_EVNDISP_string_string << std::endl;
        if( system( launching_EVNDISP_string_string ) < 0 )
        {
            std::cout << "error launching " << launching_EVNDISP_string_string << endl;
        }
    }
    else
    {
        std::cout << "Nothing needs to be analysed " << std::endl;
        bool_EVNDanalysis_done = true;
    }
    
    return;
}

//---------------------------------------------------------------------------
//-- Look the output file of EVND calibration
//-- create a list of run ready to be copied in the DB
//-- create a file with info on all run
//---------------------------------------------------------------------------
void prepare_CalibVOFF_writing( TString new_laser_run_list_name )
{

    TString list_run_analysed = get_name_run_already_analysed( new_laser_run_list_name );
    
    //-- Create File to write in DB, with valid run
    TString  wDB_file_name = get_name_run_to_be_writen( new_laser_run_list_name, true );
    FILE* wDB_file;
    wDB_file = fopen( wDB_file_name.Data(), "w" );
    //-- Create File which collect all checking information on the run
    TString  pbRunFile_name = get_name_run_problem( new_laser_run_list_name );
    FILE* pbFile;
    pbFile = fopen( pbRunFile_name.Data(), "w" );
    
    
    
    //--- read File run list missing from VOFF DB and store information
    vector< unsigned int > VERITAS_DB_LaserRunNumber;
    vector< unsigned int > VERITAS_DB_LaserConfigMask;
    vector< unsigned int > VERITAS_DB_LaserExclTel;
    vector< unsigned int > VERITAS_DB_LaserMiss;
    vector< unsigned int > VERITAS_DB_LaserDate;
    readfile_5_unsigned_int( list_run_analysed, VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VERITAS_DB_LaserMiss, VERITAS_DB_LaserDate );
    
    printf( "--- List  %s  with %d runs to be tested \n", new_laser_run_list_name.Data(), ( int ) VERITAS_DB_LaserRunNumber.size() );
    fprintf( pbFile, "--- List  %s  with %d runs to be tested \n", new_laser_run_list_name.Data(), ( int ) VERITAS_DB_LaserRunNumber.size() );
    
    // lucie: to do (mais pas dans cette fonction)
    //global vision of the run list
    
    //-- run by run check
    //-- for summary
    vector < unsigned int > run_no_pb_calib;
    vector < unsigned int > run_some_tel_pb_calib;
    vector < unsigned int > run_pb_calib;
    vector < unsigned int > run_no_pb_writing;
    vector < unsigned int > run_some_pb_writing;
    vector < unsigned int > run_pb_writing;
    
    vector < string > Summary_tel_list;
    vector < string > Summary_calib_status;
    vector < string > Summary_writing_status;
    //----------------
    
    
    //---- Loop on Run
    std::cout << "star loop on run " << std::endl;
    for( unsigned int i = 0 ; i < VERITAS_DB_LaserRunNumber.size(); i++ )
    {
        fcurrent_run = VERITAS_DB_LaserRunNumber[i];
        fcurrent_date = VERITAS_DB_LaserDate[i];
        //--------------------
        // for summary
        string tel_list = "";
        string telgoodcalib = "-1";
        string telgoodwriting = "-1";
        //--------------------
        
        //----------------------------------------------------------------------------
        //-- Check calibration
        //----------------------------------------------------------------------------
        vector< int > list_tel_missing = get_list_of_telescope( VERITAS_DB_LaserMiss[i] );
        vector< int > list_of_valid_tel;
        check_run_calib( list_tel_missing, list_of_valid_tel, pbFile );
        
        
        //--------------------------------------
        // for summary
        //-- result of the check ----- (just to have a look)
        if( list_of_valid_tel.size() > 0 )
        {
            if( list_tel_missing.size() < list_of_valid_tel.size() )
            {
                std::cout << "ERROR: BIG Problem !!!! check what the code is doing (RUN:" << fcurrent_run << ")" << std::endl;
                return;
            }
            else if( list_tel_missing.size() == list_of_valid_tel.size() )
            {
                run_no_pb_calib.push_back( fcurrent_run );
            }
            else if( list_tel_missing.size() > list_of_valid_tel.size() )
            {
                run_some_tel_pb_calib.push_back( fcurrent_run );
            }
        }
        else
        {
            run_pb_calib.push_back( fcurrent_run );
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
                    return;
                }
                else if( list_tel_missing.size() == counter )
                {
                    run_no_pb_writing.push_back( fcurrent_run );
                }
                else if( list_tel_missing.size() > counter )
                {
                    run_some_pb_writing.push_back( fcurrent_run );
                }
            }
            else
            {
                run_pb_writing.push_back( fcurrent_run );
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
        Summary_tel_list.push_back( tel_list );
        Summary_calib_status.push_back( telgoodcalib );
        Summary_writing_status.push_back( telgoodwriting );
        //-------------------------
        
        
    }
    //---- end loop on run
    
    //-----------------------------------------------------------------------------------------
    printf( "------------------- SUMMARY --------------------------\n" );
    fprintf( pbFile, "------------------- SUMMARY --------------------------\n" );
    printf( "--- List  %s  with %d runs to be tested \n", new_laser_run_list_name.Data(), ( int ) VERITAS_DB_LaserRunNumber.size() );
    fprintf( pbFile, "--- List  %s  with %d runs to be tested \n", new_laser_run_list_name.Data(), ( int ) VERITAS_DB_LaserRunNumber.size() );
    
    fprintf( pbFile, "Total number of runs = %5d \n", ( int ) VERITAS_DB_LaserRunNumber.size() );
    fprintf( pbFile, "                       %5d runs NO calib problem \n", ( int ) run_no_pb_calib.size() );
    fprintf( pbFile, "                       %5d runs with calib problem for some telescopes\n", ( int ) run_some_tel_pb_calib.size() );
    fprintf( pbFile, "                       %5d runs with calib problem for ALL telescopes\n", ( int ) run_pb_calib.size() );
    fprintf( pbFile, "                       %5d runs with no problem writing\n", ( int ) run_no_pb_writing.size() );
    //loop on run
    for( unsigned int i = 0 ; i < VERITAS_DB_LaserRunNumber.size(); i++ )
    {
    
        if( ( int ) VERITAS_DB_LaserRunNumber[i] == ( int ) forget_this_run )
        {
            continue;
        }
        
        TString calib_status = "BAD";
        TString writing_status = "BAD";
        if( Summary_tel_list[i].compare( Summary_calib_status[i] ) == 0 )
        {
            calib_status = "GOOD";
        }
        if( Summary_tel_list[i].compare( Summary_writing_status[i] ) == 0 )
        {
            writing_status = "GOOD";
        }
        
        TString Stellist = Summary_tel_list[i];
        TString Scalibstatus = Summary_calib_status[i];
        TString Swritingstatus = Summary_writing_status[i];
        
        printf( "d%8d %d %4s %4s %4s %4s %4s \n", VERITAS_DB_LaserDate[i], VERITAS_DB_LaserRunNumber[i], calib_status.Data(), writing_status.Data(), Stellist.Data(), Scalibstatus.Data(), Swritingstatus.Data() );
        fprintf( pbFile, "d%8d %d %4s %4s %4s %4s %4s \n", VERITAS_DB_LaserDate[i], VERITAS_DB_LaserRunNumber[i], calib_status.Data(), writing_status.Data(), Stellist.Data(), Scalibstatus.Data(), Swritingstatus.Data() );
        
    }
    
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
    fclose( pbFile );
    
    return;
}

//---------------------------------------------------------------------------
//-- write_calib_DB
//---------------------------------------------------------------------------
void write_calib_DB( TString new_laser_run_list_name )
{

    TString  wDB_file_name = get_name_run_to_be_writen( new_laser_run_list_name, true );
    
    
    
    VDB_CalibrationInfo* db_calib_info = new VDB_CalibrationInfo( wDB_file_name, fServer, fLOW_GAIN );
    db_calib_info->write_inVOFFLINE_DB_from_file();
    
    
    //write_inVOFFLINE_DB_from_file(query_write_DB);
    
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
            {"get_list", no_argument, 0, 'g'},
            {"save_list", no_argument, 0, 's'},
            {"download", no_argument, 0, 'd'},
            {"prep_ana", no_argument, 0, 'p'},
            {"do_ana", no_argument, 0, 'e'},
            {"do_ana_all", no_argument, 0, 'a'},
            {"calib_test", no_argument, 0, 'c'},
            {"write_DB", no_argument, 0, 'w'},
            {"write_only", required_argument, 0, 'f'},
            {"run_list", required_argument, 0, 'l'},
            {"time_start", required_argument, 0, 't'},
            {"req_version", required_argument, 0, 'v'},
            {"cal_dir", required_argument, 0, 'b'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        int c = getopt_long( argc, argv, "hgsdpeacwf:l:t:v:b:", long_options, &option_index );
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
            case 'b':
                fCalib_dir = optarg;
            case 'h':
                // EVNDISP standard way to show the option without having to recompile the whole thing
                // should work together with a script that call the binary of this compiled file
                // the scripts has to call set observatory
                char* ENV;
                ENV = getenv( "OBS_EVNDISP_AUX_DIR" );
                char readme[500];
                sprintf( readme, "cat %s/ParameterFiles/EVNDISP.updateDBlaserRUN.runparameter", ENV );
                if( system( readme ) < 0 )
                {
                    cout << "error printing readme" << endl;
                }
                exit( 0 );
                break;
            case 'g':
                fget_list_new_laser_run = true;
                break;
            case 's':
                fsave_list_laser_run = true;
                break;
            case 'd':
                fdownload_new_laser_run = true;
                break;
            case 'p':
                fprepareEVNDanalysis = true;
                break;
            case 'e':
                fdoEVNDanalysis = true;
                fprepareEVNDanalysis = true;
                break;
            case 'a':
                fdoEVNDanalysis = true;
                freANA_ALL_RUN = true;
                break;
            case 'c':
                ftestEVNDcalib = true;
                fdoEVNDanalysis = true;
                fprepareEVNDanalysis = true;
                break;
            case 'w':
                fwriteDB = true;
                ftestEVNDcalib = true;
                fdoEVNDanalysis = true;
                fprepareEVNDanalysis = true;
                break;
            case 'f':
                // lucie warning = this is test mode. This option should not exist
                fwrite_only = true;
                fwrite_only_run_list = optarg;
            case 'l':
                finput_list_run = true;
                flist_run = optarg;
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
bool test_gain_toff( TString file_root_name_gain, TString file_root_name_toff , FILE*& pbFile )
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
            fprintf( pbFile, "ERROR: gain RMS %f too large\n", h_gain->GetRMS() );
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
            fprintf( pbFile, "ERROR: %4f channels with gain at zero \n \n", h_gz_channel->GetEntries() );
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
        fprintf( pbFile, "ERROR: File %s does not exist \n", file_root_name_gain.Data() );
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
            fprintf( pbFile, "ERROR: %4f channels with toff at zero \n \n", h_tz_channel->GetEntries() );
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
        fprintf( pbFile, "ERROR: File %s does not exist \n", file_root_name_toff.Data() );
        good_calib = false;
    }
    //----------------------- compare channel's gain and toff at zero
    if( V_pixel_zero_gain.size() > 0 || V_pixel_zero_toff.size() > 0 )
    {
        if( V_pixel_zero_gain.size() == V_pixel_zero_toff.size() )
        {
            printf( "Dead pixels: \n" );
            fprintf( pbFile, "Dead pixels: \n" );
            for( unsigned int j = 0 ; j < V_pixel_zero_gain.size(); j++ )
            {
                printf( " %d ", V_pixel_zero_gain[j] );
                fprintf( pbFile, " %d ", V_pixel_zero_gain[j] );
            }
            printf( " \n" );
            fprintf( pbFile, " \n " );
            
        }
        else
        {
            if( TMath::Abs( ( Int_t ) V_pixel_zero_gain.size() - ( Int_t ) V_pixel_zero_toff.size() ) > 2 )
            {
                printf( "ERROR: NOT same number of pixel with zero gain (%4d) and zero toff (%4d) \n", ( int ) V_pixel_zero_gain.size(), ( int ) V_pixel_zero_toff.size() );
                fprintf( pbFile, "ERROR: NOT same number of pixel with zero gain (%4d) and zero toff (%4d)  \n", ( int ) V_pixel_zero_gain.size(), ( int ) V_pixel_zero_toff.size() );
                good_calib = false;
            }
            else
            {
                printf( "WARNING: NOT same number of pixel with zero gain (%4d) and zero toff (%4d) \n", ( int ) V_pixel_zero_gain.size(), ( int ) V_pixel_zero_toff.size() );
                fprintf( pbFile, "WARNING: NOT same number of pixel with zero gain (%4d) and zero toff (%4d)  \n", ( int ) V_pixel_zero_gain.size(), ( int ) V_pixel_zero_toff.size() );
                good_calib = true;
            }
            for( unsigned int ii = 0; ii < V_pixel_zero_gain.size(); ii++ )
            {
                printf( " %4d zero gain \n", V_pixel_zero_gain[ii] );
                fprintf( pbFile, " %4d zero gain \n", V_pixel_zero_gain[ii] );
            }
            for( unsigned int ii = 0; ii < V_pixel_zero_toff.size(); ii++ )
            {
                printf( " %4d zero toff \n", V_pixel_zero_toff[ii] );
                fprintf( pbFile, " %4d zero toff \n", V_pixel_zero_toff[ii] );
            }
        }
    }
    //-------------------------------------------------------------------
    return good_calib;
}


//---------------------------------------------------------------------------
//-- unsigned long Check_telmissing_from_VOFFDB
//---------------------------------------------------------------------------
//-- used in TString get_new_laser_run_list()
//-- not fit to be called on its own
//-------------------
//-- Check_telmissing_from_VOFFDB
//-- read VOFFLINE_DB_LaserRunNumber_Tel, which should be properly filled, as it is done in get_new_laser_run_list()
//-- Check if laser run in VERITAS DB are present in VOFFLINE DB
//-- return mask of telescopes missing from VOFFLINE DB (which should be there)
//---------------------------------------------------------------------------
unsigned long Check_telmissing_from_VOFFDB( unsigned int VERITAS_DB_LaserRunNumber_i_run, unsigned int VERITAS_DB_LaserConfigMask_i_run, unsigned int VERITAS_DB_LaserExclTel_i_run, vector< vector < unsigned int > > VOFFLINE_DB_LaserRunNumber_Tel )
{
    //    printf ("-------- %d %d %d \n",VERITAS_DB_LaserRunNumber_i_run,VERITAS_DB_LaserConfigMask_i_run, VERITAS_DB_LaserExclTel_i_run);
    bitset<4> missing_tel; // with this we list the telescope missing from VOFFLINE DB
    missing_tel.reset();
    fcurrent_run = VERITAS_DB_LaserRunNumber_i_run;
    vector<int> list_of_tel = get_list_of_telescope( VERITAS_DB_LaserConfigMask_i_run, VERITAS_DB_LaserExclTel_i_run );
    //--- loop on the telescopes participating to the run
    for( unsigned i_tel = 0; i_tel < list_of_tel.size(); i_tel++ )
    {
        int tel_indice = list_of_tel[i_tel] - 1;
        //	    std::cout<<"Tel  "<<list_of_tel[i_tel]<<std::endl;
        //--- loop on the run in VOFFLINE DB
        bool run_in_voff_db = false;
        for( unsigned int i_run_voff = 0; i_run_voff < VOFFLINE_DB_LaserRunNumber_Tel[tel_indice].size(); i_run_voff++ )
        {
            if( VERITAS_DB_LaserRunNumber_i_run == VOFFLINE_DB_LaserRunNumber_Tel[tel_indice][i_run_voff] )
            {
                //		    std::cout<<"VOFFLINE == VERITAS "<<std::endl;
                run_in_voff_db = true;
                break;
            }
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
//-- read_laserRUN_fromVOFFLINE_DB                                         --
//-- for each telescope, get the list of run in the data base              --
//---------------------------------------------------------------------------
void read_laserRUN_fromVOFFLINE_DB( vector < vector < unsigned int > >& VOFFLINE_DB_LaserRunNumber_Telnum, vector < vector < string > >& VOFFLINE_DB_LaserDate_Telnum, vector < vector < string > >& VOFFLINE_DB_LaserVersion_Telnum, unsigned int Tel_num )
{
    string query = WriteQuery_to_getLaserRun_fromVOFFLINE_DB( Tel_num );
    
    vector< unsigned int >  i_VOFFLINE_DB_LaserRunNumber_Tel;
    vector< string >  i_VOFFLINE_DB_Laserdate_Tel;
    vector< string >  i_VOFFLINE_DB_Laserversion_Tel;
    
    //---- open the DB and check
    string iTempS;
    iTempS =  fServer + "/VOFFLINE";
    //TSQLServer *f_db = TSQLServer::Connect( iTempS.c_str(), "readonly", "" );
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ; // lucie: DO NOT create a pointer. this way the connection is automatically closed when getting out of the function
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "ERROR read_laserRUN_fromVOFFLINE_DB: failed to connect to database server" << endl;
        cout << "\t server: " <<  fServer << endl;
        VOFFLINE_DB_LaserRunNumber_Telnum.push_back( i_VOFFLINE_DB_LaserRunNumber_Tel );
        VOFFLINE_DB_LaserDate_Telnum.push_back( i_VOFFLINE_DB_Laserdate_Tel );
        VOFFLINE_DB_LaserVersion_Telnum.push_back( i_VOFFLINE_DB_Laserversion_Tel );
        return;
    }
    //---- do the query and check
    if( !my_connection.make_query( query.c_str() ) )
    {
        cout << "WARNING read_laserRUN_fromVOFFLINE_DB: failed to get something from the query " << endl;
        VOFFLINE_DB_LaserRunNumber_Telnum.push_back( i_VOFFLINE_DB_LaserRunNumber_Tel );
        VOFFLINE_DB_LaserDate_Telnum.push_back( i_VOFFLINE_DB_Laserdate_Tel );
        VOFFLINE_DB_LaserVersion_Telnum.push_back( i_VOFFLINE_DB_Laserversion_Tel );
        return;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    //---- read the query
    if( db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( atoi( db_row->GetField( 0 ) ) == forget_this_run )
            {
                continue;
            }
            
            i_VOFFLINE_DB_LaserRunNumber_Tel.push_back( atoi( db_row->GetField( 0 ) ) );
            string version = string( db_row->GetField( 1 ) );
            i_VOFFLINE_DB_Laserversion_Tel.push_back( version ) ;
            string date = string( db_row->GetField( 2 ) );
            i_VOFFLINE_DB_Laserdate_Tel.push_back( date ) ;
        }
    }
    else
    {
        cout << "WARNING read_laserRUN_fromVOFFLINE_DB:  no laser run found " << endl;
        
    }
    //-- close the DB
    //f_db->Close(); // lucie
    
    VOFFLINE_DB_LaserRunNumber_Telnum.push_back( i_VOFFLINE_DB_LaserRunNumber_Tel );
    VOFFLINE_DB_LaserDate_Telnum.push_back( i_VOFFLINE_DB_Laserdate_Tel );
    VOFFLINE_DB_LaserVersion_Telnum.push_back( i_VOFFLINE_DB_Laserversion_Tel );
    return;
}

//---------------------------------------------------------------------------
//-- read_laserRUN_fromVERITAS_DB                                          --
//-- get the list of all valid laser run in the data base                  --
//---------------------------------------------------------------------------
void read_laserRUN_fromVERITAS_DB( vector< unsigned int >& VERITAS_DB_LaserRunNumber, vector< unsigned int >& VERITAS_DB_LaserConfigMask, vector< unsigned int >& VERITAS_DB_LaserExclTel, vector< unsigned int >& VERITAS_DB_LaserDate )
{

    string query = WriteQuery_to_getLaserRun_fromVERITAS_DB();
    
    //---- open the DB and check
    string iTempS;
    iTempS =  fServer;
    iTempS += "/VERITAS";
    //TSQLServer *f_db = TSQLServer::Connect( iTempS.c_str(), "readonly", "" );
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ; // lucie: DO NOT create a pointer. this way the connection is automatically closed when getting out of the function
    if( !my_connection.Get_Connection_Status() )
    {
        //if( !f_db ){
        cout << "ERROR read_laserRUN_fromVERITAS_DB: failed to connect to database server" << endl;
        cout << "\t server: " <<  fServer << endl;
        return;
    }
    //---- do the query and check
    if( !my_connection.make_query( query.c_str() ) )
    {
        return;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    //---- read the query
    if( db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( db_row->GetField( 0 ) )
            {
                VERITAS_DB_LaserRunNumber.push_back( atoi( db_row->GetField( 0 ) ) );
            }
            else
            {
                std::cout << "WARNING: no run number  in VERITAS DB" << std::endl;
                VERITAS_DB_LaserRunNumber.push_back( -1 );
            }
            if( db_row->GetField( 1 ) )
            {
                VERITAS_DB_LaserConfigMask.push_back( atoi( db_row->GetField( 1 ) ) );
            }
            else
            {
                std::cout << "WARNING: no config Mask in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                VERITAS_DB_LaserConfigMask.push_back( -1 );
            }
            if( db_row->GetField( 2 ) )
            {
                VERITAS_DB_LaserExclTel.push_back( atoi( db_row->GetField( 2 ) ) );
            }
            else
            {
                std::cout << "WARNING: no exclud Tel  in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                VERITAS_DB_LaserExclTel.push_back( -1 );
            }
            if( db_row->GetField( 3 ) )
            {
                VERITAS_DB_LaserDate.push_back( get_date_from_tblRun_Info_data_start_time( db_row->GetField( 3 ) ) );
            }
            else
            {
                std::cout << "WARNING: no data start date in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                std::cout << "using db_start_time instead" << std::endl;
                if( db_row->GetField( 4 ) )
                {
                    VERITAS_DB_LaserDate.push_back( get_date_from_tblRun_Info_data_start_time( db_row->GetField( 4 ) ) );
                }
                else
                {
                    std::cout << "WARNING: no db start date in VERITAS DB for run " << atoi( db_row->GetField( 0 ) ) << std::endl;
                    VERITAS_DB_LaserDate.push_back( -1 );
                }
            }
        }
    }
    else
    {
        cout << "WARNING read_laserRUN_fromVERITAS_DB:  no laser run found " << endl;
    }
    //-- close the DB
    //f_db->Close(); // lucie (closing automatically with the automatic deletion of my_connection)
    
    std::cout << "************************************************************ " << VERITAS_DB_LaserExclTel.size() << " size vect excluded tel" << std::endl;
    
    return;
    
}

//---------------------------------------------------------------------------
//-- WriteQuery_to_getAllLaserRun_fromVERITAS_DB()
//---------------------------------------------------------------------------
string WriteQuery_to_getLaserRun_fromVERITAS_DB()
{

    char c_query[flong_char_query];// has to be long if we ask for a long run list in VDBSourceInfo
    
    //sprintf( c_query, "SELECT big_table.run_id, big_table.config_mask, big_table.excluded_telescopes, big_table.data_start_time , big_table.db_start_time FROM (SELECT Info.run_id, Info.run_type, Info.config_mask,Info.data_start_time,Info.db_start_time , grp_cmt.excluded_telescopes, grp_cmt.group_type, grp_cmt.group_id FROM tblRun_Info AS Info, tblRun_Group AS grp, tblRun_GroupComment AS grp_cmt WHERE  grp_cmt.group_type = 'laser' AND grp_cmt.group_id = grp.group_id AND grp.run_id = Info.run_id AND (Info.run_type = 'laser' OR Info.run_type = 'flasher') ORDER BY grp_cmt.excluded_telescopes ) AS big_table GROUP BY run_id;");
    sprintf( c_query, "SELECT big_table.run_id, big_table.config_mask, big_table.excluded_telescopes, big_table.data_start_time , big_table.db_start_time FROM (SELECT Info.run_id, Info.run_type, Info.config_mask,Info.data_start_time,Info.db_start_time , grp_cmt.excluded_telescopes, grp_cmt.group_type, grp_cmt.group_id FROM tblRun_Info AS Info, tblRun_Group AS grp, tblRun_GroupComment AS grp_cmt WHERE  grp_cmt.group_type = 'laser' AND grp_cmt.group_id = grp.group_id AND grp.run_id = Info.run_id AND (Info.run_type = 'laser' OR Info.run_type = 'flasher') ORDER BY grp_cmt.excluded_telescopes ) AS big_table GROUP BY run_id,excluded_telescopes ;" );
    
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
    //                     grp_cmt.group_type = 'laser'
    //                 AND grp_cmt.group_id = grp.group_id
    //                 AND grp.run_id = Info.run_id
    //                 AND (Info.run_type = 'laser' OR Info.run_type = 'flasher')
    //
    //               ORDER BY
    //                     grp_cmt.excluded_telescopes
    //               )
    //         AS big_table
    //
    //         GROUP BY big_table.run_id;
    
    string query = c_query;
    return query;
}
//---------------------------------------------------------------------------
//-- WriteQuery_to_getAllLaserRun_fromVOFFLINE_DB()
//---------------------------------------------------------------------------
string WriteQuery_to_getLaserRun_fromVOFFLINE_DB( unsigned int Tel_num )
{

    //begining of the querry
    string query = "SELECT big_table.run_id, big_table.code_version, big_table.update_time  FROM ( SELECT run_id, code_version, update_time  FROM tblEventDisplay_Analysis_Calibration_Flasher WHERE ";
    // select the telescope
    char c_query[flong_char_query];// has to be long if we ask for a long run list in VDBSourceInfo
    sprintf( c_query, "telescope = %u", Tel_num );
    string tel_query = c_query;
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

//---------------------------------------------------------------------------
//-- create_list_for_download
//-- read the file with all info from run missing from VOFFLINE DB
//-- create a file with the right structure to download those run
//---------------------------------------------------------------------------
void create_list_for_download( TString file_name )
{

    FILE* NEW_laser_list_File_for_download;
    TString NEW_laser_list_for_download = get_name_for_download( file_name );
    NEW_laser_list_File_for_download = fopen( NEW_laser_list_for_download.Data(), "w" );
    
    //--- read File run list missing from VOFF DB and store information
    vector< unsigned int > VERITAS_DB_LaserRunNumber;
    vector< unsigned int > VERITAS_DB_LaserConfigMask;
    vector< unsigned int > VERITAS_DB_LaserExclTel;
    vector< unsigned int > VERITAS_DB_LaserMiss;
    vector< unsigned int > VERITAS_DB_LaserDate;
    readfile_5_unsigned_int( file_name, VERITAS_DB_LaserRunNumber, VERITAS_DB_LaserConfigMask, VERITAS_DB_LaserExclTel, VERITAS_DB_LaserMiss, VERITAS_DB_LaserDate );
    
    for( unsigned int i = 0 ; i < VERITAS_DB_LaserRunNumber.size(); i++ )
    {
        fprintf( NEW_laser_list_File_for_download, "%d %d \n", VERITAS_DB_LaserRunNumber[i], VERITAS_DB_LaserDate[i] );
    }
    
    fclose( NEW_laser_list_File_for_download );
    return;
}
//---------------------------------------------------------------------------
//-- readfile_5_unsigned_int
//---------------------------------------------------------------------------
void readfile_5_unsigned_int( TString file_name, vector< unsigned int >& VERITAS_DB_LaserRunNumber, vector< unsigned int >& VERITAS_DB_LaserConfigMask, vector< unsigned int >& VERITAS_DB_LaserExclTel, vector< unsigned int >& VERITAS_DB_LaserMiss, vector< unsigned int >& VERITAS_DB_LaserDate )
{

    if( is_file_empty( file_name ) )
    {
        return;
    }
    else
    {
    
        ifstream file;
        file.open( file_name );
        
        unsigned RunNumber = -1;
        unsigned ConfigMask = -1;
        unsigned ExclTel = -1;
        unsigned Miss = -1;
        unsigned Date = -1;
        
        while( 1 )
        {
            file >> RunNumber >> ConfigMask >> ExclTel >> Miss >> Date;
            if( !file.good() )
            {
                break;
            }
            VERITAS_DB_LaserRunNumber.push_back( RunNumber );
            VERITAS_DB_LaserConfigMask.push_back( ConfigMask );
            VERITAS_DB_LaserExclTel.push_back( ExclTel );
            VERITAS_DB_LaserMiss.push_back( Miss );
            VERITAS_DB_LaserDate.push_back( Date );
        }
        file.close();
        return;
    }
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
// and write the pb in the pbFile
//---------------------------------------------------------------------------
void check_run_calib( vector< int >  ListTel , vector< int >& list_of_valid_tel, FILE*& pbFile )
{

    // write introduction for the run
    printf( "-------------------------------------------------------------------------------- \n" );
    fprintf( pbFile, "---------------------------------------------------------------------------------- \n" );
    printf( "--- RUN %d      has %d telescope(s) to be tested \n", ( int ) fcurrent_run, ( int ) ListTel.size() );
    fprintf( pbFile, "--- RUN %d      has %d telescope(s) to be tested \n", ( int ) fcurrent_run, ( int ) ListTel.size() );
    
    TString ltel_not_good = "";
    
    //-- Loop on the telescope
    for( unsigned int i = 0; i < ListTel.size() ; i++ )
    {
    
    
        fcurrent_tel = ListTel[i];
        printf( "--- Tel %d \n", fcurrent_tel );
        fprintf( pbFile, "--- Tel %d \n", fcurrent_tel );
        //-- get the name of the calibration root file to be checked
        TString file_gain_root = get_name_calib_file( fCalib_dir, "gain.root", ListTel[i], fcurrent_run );
        TString file_toff_root = get_name_calib_file( fCalib_dir, "toff.root", ListTel[i], fcurrent_run );
        //-- test the calib file for the run and a given telescope
        if( test_gain_toff( file_gain_root, file_toff_root, pbFile ) )
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
    fprintf( pbFile, "--- RUN %d      has %d/%d telescope(s) with valid calibration: \n", fcurrent_run, ( int ) list_of_valid_tel.size(), ( int ) ListTel.size() );
    //  std::cout<<"--- 1 "<<std::endl;
    fprintf( pbFile, "--- (Display: $EVNDISPSYS/scripts/VTS/VTS.EVNDISP.display %s %s ) \n", ltel_not_good.Data(), get_downloaded_laser_path( fcurrent_run, fcurrent_date ).Data() ); // bob done to be tested
    //   std::cout<<"--- 2 "<<std::endl;
    printf( "--- (Display: $EVNDISPSYS/scripts/VTS/VTS.EVNDISP.display %s %d ) \n", ltel_not_good.Data(), fcurrent_run );
    // std::cout<<"--- 3 "<<std::endl;
    
    if( list_of_valid_tel.size() > 0 )
    {
        for( unsigned int i = 0; i < list_of_valid_tel.size(); i++ )
        {
            printf( "---------  - %d %d \n", fcurrent_run, list_of_valid_tel[i] );
            fprintf( pbFile, "---------  - %d %d  \n", fcurrent_run, list_of_valid_tel[i] );
        }
    }
    if( list_of_valid_tel.size() < ListTel.size() )
    {
        printf( "--- Should run the Display: \n $EVNDISPSYS/scripts/VTS/VTS.EVNDISP.display %s %s \n (May have to find replace this laser run) \n", ltel_not_good.Data(), get_downloaded_laser_path( fcurrent_run, fcurrent_date ).Data() );
        fprintf( pbFile, "--- Should run the Display: \n $EVNDISPSYS/scripts/VTS/VTS.EVNDISP.display %s %s \n (May have to replace this laser run)\n", ltel_not_good.Data(), get_downloaded_laser_path( fcurrent_run, fcurrent_date ).Data() );
        
    }
    
    
    return;
}
//---------------------------------------------------------------------------
//-- independent use version
//---------------------------------------------------------------------------
void check_run_calib( unsigned int RunNumber, int  mask, int excluded_tel_mask )
{


    TString name = "Run_";
    name += RunNumber;
    name += ".check_calib";
    FILE* pbRunFile;
    pbRunFile = fopen( name.Data(), "w" );
    vector< int > ListTel = get_list_of_telescope( mask, excluded_tel_mask );
    vector< int > list_of_valid_tel;
    
    fcurrent_run = RunNumber;
    check_run_calib( ListTel , list_of_valid_tel, pbRunFile );
    
    fclose( pbRunFile );
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


