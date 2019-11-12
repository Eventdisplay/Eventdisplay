/*! \class VDBTools
    \brief routines to access different tables of the DB


*/

#include "VDBTools.h"

///////////////////////////////////////////////////////////////////////////////////

VDB_ObservingSources_Data::VDB_ObservingSources_Data()
{
    SourceID = "";
    fDec = 0.;
    fRA = 0.;
    fEpoch = 0.;
}

///////////////////////////////////////////////////////////////////////////////////

VDB_ObservingSources::VDB_ObservingSources()
{
    setDebug( false );
}

bool VDB_ObservingSources::fill( TSQLServer* i_db )
{
    if( !i_db )
    {
        return false;
    }
    
    fVDB_ObservingSources_Data.clear();
    
    char c_query[1000];
    
    // read source data
    sprintf( c_query, "select * from tblObserving_Sources" );
    
    TSQLResult* db_res = i_db->Query( c_query );
    if( !db_res )
    {
        return false;
    }
    
    int fNRows = db_res->GetRowCount();
    
    cout << "reading source info from DB" << endl;
    cout << "\t number of sources read form DB: " << fNRows << endl;
    
    string itemp;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();
        
        if( !db_row->GetField( 0 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 0 );
        fVDB_ObservingSources_Data[itemp] = new VDB_ObservingSources_Data();
        fVDB_ObservingSources_Data[itemp]->SourceID = itemp;
        fVDB_ObservingSources_Data[itemp]->fDec = atof( db_row->GetField( 2 ) ) * 180. / TMath::Pi();
        fVDB_ObservingSources_Data[itemp]->fRA  = atof( db_row->GetField( 1 ) ) * 180. / TMath::Pi();
        fVDB_ObservingSources_Data[itemp]->fEpoch =  atof( db_row->GetField( 3 ) );
    }
    return true;
}

VDB_ObservingSources_Data* VDB_ObservingSources::get_ObservingSources_Data( string iSource )
{
    if( fVDB_ObservingSources_Data.find( iSource ) != fVDB_ObservingSources_Data.end() )
    {
        return fVDB_ObservingSources_Data[iSource];
    }
    
    cout << "VDB_ObservingSources::get_ObservingSources_Data: warning: source " << iSource << " not found" << endl;
    return 0;
}

void VDB_ObservingSources::list()
{
    map< string, VDB_ObservingSources_Data* >::iterator i_iter;
    
    for( i_iter = fVDB_ObservingSources_Data.begin(); i_iter != fVDB_ObservingSources_Data.end(); ++i_iter )
    {
        VDB_ObservingSources_Data* i_temp = ( VDB_ObservingSources_Data* )i_iter->second;
        if( i_temp )
        {
            cout << i_temp->SourceID << "\t";
            cout << i_temp->fDec << "\t";
            cout << i_temp->fRA << "\t";
            cout << i_temp->fEpoch << "\t";
            cout << endl;
        }
    }
}
