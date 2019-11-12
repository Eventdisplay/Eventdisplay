/*! \class VTrackingCorrections
    \brief read tracking corrections from DB and apply corrections to pointing

    use code from Steve Fegan to apply tracking corrections (thanks!)


*/

#include "VTrackingCorrections.h"

VTrackingCorrections::VTrackingCorrections( unsigned int iTelID )
{
    fStatus = false;
    fTelID = iTelID;
    
    correctionParameters = 0;
}


bool VTrackingCorrections::readTrackingCorrectionsFromDB( string iSQLDate )
{
    cout << "Reading T-point tracking corrections from VERITAS database for telescope " << fTelID + 1 << " at " << iSQLDate << endl;
    string iTemp = getDBServer() + "/VERITAS";
    
    //std::cout<<"VTrackingCorrections::readTrackingCorrectionsFromDB "<<std::endl;
    VDB_Connection my_connection( iTemp.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "VTrackingCorrections: failed to connect to database server" << endl;
        fStatus = false;
        return false;
    }
    char c_query[1000];
    
    sprintf( c_query, "select * from tblPositioner_Telescope%d_Corrections where db_start_time <= \"%s\" ", fTelID, iSQLDate.c_str() );
    
    if( !my_connection.make_query( c_query ) )
    {
        cout << "VTrackingCorrections: failed to find correct tables" << endl;
        cout << "\t" << c_query << endl;
        fStatus = false;
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    
    int fNRows = db_res->GetRowCount();
    
    TSQLRow* db_row = 0;
    for( int i = 0; i < fNRows; i++ )
    {
        db_row = db_res->Next();
    }
    if( db_row )
    {
        string itemp;
        cout << "0: " << db_row->GetField( 0 ) << ", ";
        if( db_row->GetField( 1 ) )
        {
            cout << "1: " << db_row->GetField( 1 ) << ", ";
        }
        for( int i = 2; i < 12; i++ )
        {
            cout << i << ": " << db_row->GetField( i ) << ", ";
        }
        cout << endl;
        for( int i = 12; i < 24; i++ )
        {
            cout << i << ": " << db_row->GetField( i ) << ", ";
        }
        cout << endl;
        correctionParameters = new CorrectionParameters();
        itemp = db_row->GetField( 2 );
        if( itemp == "enabled" )
        {
            correctionParameters->enable_offsets = true;
        }
        else
        {
            correctionParameters->enable_offsets = false;
        }
        itemp = db_row->GetField( 3 );
        if( itemp == "enabled" )
        {
            correctionParameters->enable_corrections = true;
        }
        else
        {
            correctionParameters->enable_corrections = false;
        }
        itemp = db_row->GetField( 14 );
        if( itemp == "enabled" )
        {
            correctionParameters->enable_vff = true;
        }
        else
        {
            correctionParameters->enable_vff = false;
        }
        correctionParameters->az_ratio = atof( db_row->GetField( 4 ) );
        correctionParameters->el_ratio = atof( db_row->GetField( 5 ) );
        correctionParameters->az_offset = atof( db_row->GetField( 6 ) );
        correctionParameters->el_offset = atof( db_row->GetField( 7 ) );
        correctionParameters->az_ns = atof( db_row->GetField( 8 ) );
        correctionParameters->az_ew = atof( db_row->GetField( 9 ) );
        correctionParameters->el_udew = atof( db_row->GetField( 10 ) );
        correctionParameters->fp_az = atof( db_row->GetField( 11 ) );
        correctionParameters->flex_el_A = atof( db_row->GetField( 12 ) );
        correctionParameters->flex_el_B = atof( db_row->GetField( 13 ) );
        correctionParameters->el_pos_vff_s = atof( db_row->GetField( 15 ) );
        correctionParameters->el_pos_vff_t = atof( db_row->GetField( 16 ) );
        correctionParameters->el_neg_vff_s = atof( db_row->GetField( 17 ) );
        correctionParameters->el_neg_vff_t = atof( db_row->GetField( 18 ) );
        correctionParameters->az_pos_vff_s = atof( db_row->GetField( 19 ) );
        correctionParameters->az_pos_vff_t = atof( db_row->GetField( 20 ) );
        correctionParameters->az_neg_vff_s = atof( db_row->GetField( 21 ) );
        correctionParameters->az_neg_vff_t = atof( db_row->GetField( 22 ) );
    }
    
    return true;
}


bool VTrackingCorrections::applyTrackingCorrections( double iElRaw, double iAzRaw, double& iElCorr, double& iAzCorr )
{
    double iEl = iElRaw;
    double iAz = iAzRaw;
    
    correctionParameters->undoAzElCorrections( iAz, iEl, true );
    
    iElCorr = iEl;
    iAzCorr = iAz;
    
    return fStatus;
}
