//! VDB_CalibrationInfo read or get laser run calibration info from VOFFLINE DB


#ifndef VDB_CALIBRATIONINFO_H
#define VDB_CALIBRATIONINFO_H

#include <TMath.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>
#include <TSystem.h>

#include <bitset>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "VDB_Connection.h"

using namespace std;

class VDB_CalibrationInfo
{
    protected:
    
        int fLowGain;
        TString fServer;
        TString fFile_to_write;
        bool fwrite_flag;
        bool fdofile_flag;
        bool fread_flag;
        int fcurrent_run;
        int fcurrent_tel;
        TString fNOW_DB;
        int fVOFFLINE_version_query;
        int fLOW_GAIN;
        
        int freading_gain_or_toff;//gain_or_toff : 1 = gain, 2=toff
        
        TString fseparation;
        string fquery_read;
        int flong_char_query;
        
        //------------------- Writing function related
        string please_give_the_password();//needed
        string WriteQuery_to_write_in_DB();//needed
        
        bool test_file_format( TString file_to_be_copied ); // security TO DO
        
        //------------------- reading function
        void Create_query_read();
        bool Read_the_DB();
        
        // data vectors
        vector < unsigned int > Vchannel;
        vector < double > Vmean;
        vector < double > Vvar;
        
        
        
    public:
        VDB_CalibrationInfo(); // default constructor
        
        //-- constructor to then call DoFile_for_DBwriting
        VDB_CalibrationInfo( int current_run, int current_tel, TString NOW_DB, int VOFFLINE_version_query, int LOW_GAIN );
        bool DoFile_for_DBwriting( vector < double > Vchannel_gain, vector < double > Vmean_gain, vector < double > Vvar_gain, vector < double > Vchannel_toff, vector < double > Vmean_toff, vector < double > Vvar_toff, FILE*& wDB_file );
        
        //-- constructor to then call write_inVOFFLINE_DB_from_file()
        VDB_CalibrationInfo( TString file_to_write_in_DB, TString DBserver, int low_gain ); // needed
        void write_inVOFFLINE_DB_from_file( string pass_word = "" );
        
        //-- constructor to then call read
        VDB_CalibrationInfo( int laserrun , int tel , string name_out_file, int gain_or_toff, int VOFFLINE_version_query, int LOW_GAIN, TString DBserver );
        bool readVOFFLINE();
        vector< unsigned int > getVectorChannelList()
        {
            return Vchannel;
        }
        vector< double > getVectorMean()
        {
            return Vmean;
        }
        vector< double > getVectorVariance()
        {
            return Vvar;
        }
        
        ~VDB_CalibrationInfo() {}
        
        
};

#endif





