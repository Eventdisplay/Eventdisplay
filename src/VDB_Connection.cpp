/*! \class VDB_Connection
    \brief connect to DB

*/

#include "VDB_Connection.h"
//======================== PUBLIC ================

/******************** CONSTRUCTORS ***************/
VDB_Connection::VDB_Connection()
{
    fMAX_PROCESS = 400;
    fNumb_Connection = -111;
    
    fDBserver = "";
    fconnection_mode = "";
    fconnection_option = "";
    fDB_Connection_successfull = false;
    fDB_Query_successfull = false;
    f_db = 0;
    
    
}

VDB_Connection::VDB_Connection( string DBserver, string connection_mode, string connection_option )
{

    fMAX_PROCESS = 400;
    fNumb_Connection = -111;
    
    fDBserver = DBserver;
    fconnection_mode = connection_mode;
    fconnection_option = connection_option;
    fDB_Connection_successfull = false;
    fDB_Query_successfull = false;
    f_db = 0;
    Connect();
    
}

bool VDB_Connection::Connect()
{

    // Connect
    f_db = TSQLServer::Connect( fDBserver.c_str(), fconnection_mode.c_str() , fconnection_option.c_str() );
    fDB_Connection_successfull = false;
    // Test the connection
    
    if( !f_db )
    {
        // connection failed
        cout << "VDB_Connection::Connect() info: failed to connect to database server, sleep for 60 s and try again..." << endl;
        gSystem->Sleep( 60000 );
        // try again
        f_db = TSQLServer::Connect( fDBserver.c_str(), fconnection_mode.c_str() , fconnection_option.c_str() );
        // give up
        if( !f_db )
        {
            cout << "VDB_Connection::Connect() Error: failed to connect to database server" << endl;
            cout << "\t server: " << fDBserver  << endl;
            Close_Connection();
            exit( EXIT_FAILURE );
            return fDB_Connection_successfull;
        }
        
    }
    else if( Get_Nb_Connection() > fMAX_PROCESS )
    {
        // too may connection are already opened
        Close_Connection();
        cout << "VDB_Connection::Connect() info: there were too many DB connections (" << fNumb_Connection << " > " << fMAX_PROCESS << "). Sleep for 60 s and try again..." << endl;
        gSystem->Sleep( 60000 );
        
        // try again
        unsigned int nTry = 3;
        for( unsigned int i = 1; i <= nTry; i++ )
        {
            f_db = TSQLServer::Connect( fDBserver.c_str(), fconnection_mode.c_str() , fconnection_option.c_str() );
            if( Get_Nb_Connection() > fMAX_PROCESS )
            {
                Close_Connection();
                cout << "VDB_Connection::Connect() Error: there were still too many DB connections (";
                cout << fNumb_Connection << " > " << fMAX_PROCESS << ")" << endl;
                if( i < nTry )
                {
                    cout << "\t Sleep for 60 s and try again..." << endl;
                }
                else
                {
                    cout << "\t server: " <<  fDBserver << endl;
                    cout << "\t still to many connections open - close and exiting now..." << endl;
                    exit( EXIT_FAILURE );
                }
            }
            return fDB_Connection_successfull;
        }
    }
    
    fDB_Connection_successfull = true;
    return fDB_Connection_successfull;
}

bool VDB_Connection::make_query( const char* the_query )
{

    //cout << "Query: " << the_query << endl;
    fDB_Query_successfull = false;
    
    if( !f_db )
    {
    
        return fDB_Query_successfull;
        
    }
    else
    {
    
        fdb_res = f_db->Query( the_query );
        
        if( !fdb_res )
        {
            std::cout << "VDB_Connection::make_query no result for query:  " << the_query << std::endl;
            return fDB_Query_successfull;
        }
        else
        {
            fDB_Query_successfull = true;
        }
        
    }
    
    return fDB_Query_successfull;
}

int  VDB_Connection::Get_Nb_Connection()
{

    if( f_db && f_db->Query( "show processlist" ) )
    {
        fNumb_Connection = f_db->Query( "show processlist" )->GetRowCount();
    }
    else
    {
        return -111;
    }
    
    return fNumb_Connection;
    
}


