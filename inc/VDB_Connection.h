//! VDB_Connection
// small class to connect to DB


#ifndef VDB_CONNECTION_H
#define VDB_CONNECTION_H

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

using namespace std;

class VDB_Connection
{
    protected:
    
        TSQLServer* f_db;
        TSQLResult* fdb_res;
        
        bool Connect();
        bool fDB_Connection_successfull;
        bool fDB_Query_successfull;
        
        string fDBserver;
        string fconnection_mode;
        string fconnection_option;
        
        int fMAX_PROCESS;
        int fNumb_Connection;
        
    public:
        VDB_Connection(); // default constructor
        
        VDB_Connection( string DBserver, string connection_mode, string connection_option );
        
        ~VDB_Connection()
        {
            if( f_db )
            {
                Close_Connection();
            }
            
            
        }
        
        bool make_query( const char* the_query );
        
        bool Get_Connection_Status()
        {
            return fDB_Connection_successfull;
        }
        bool Get_Query_Status()
        {
            return fDB_Query_successfull;
        }
        TSQLResult* Get_QueryResult()
        {
            return fdb_res;
        }
        TSQLServer* Get_ConnectionResult()
        {
            return f_db;
        }
        
        
        void Close_Connection()
        {
            if( f_db )
            {
                f_db->Close();
                f_db = 0;
            }
            return;
        }
        
        int  Get_Nb_Connection();
        
};

#endif





