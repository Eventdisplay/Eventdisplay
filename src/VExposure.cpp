/*! \class VExposure
    retrieve a list of runs, fill and plot VERITAS exposure plots, fill elevation plots

    Note: some of the usage of this class is counter-intuitive


*/

#include "VExposure.h"


ClassImp( VExposure )

VExposure::VExposure( int nBinsL, int nBinsB )
{
    fDebug = false;
    fMakeRunList = false;
    
    bPlotElevationPlots = false;
    
    bPrintVerbose = false;
    bPrintTimeMask = false;
    
    fStartDate_SQL = "2007-09-01";
    fStopDate_SQL =  "2009-12-31";
    
    fAcceptance = 0;
    
    fPlotExtendedSources = true;
    fPlotSourceNames = true;
    fVDB_ObservingSources = 0;
    
    setMaximumIntegrationRadius();
    setCanvasSize();
    
    if( nBinsL > 0 && nBinsB > 0 )
    {
        double iBMin = -90.;
        double iBMax = 90.;
        
        fMapGal2D = new TH2D( "fMapGal2D", "", nBinsL * 2, -180., 180., nBinsB * 2, iBMin, iBMax );
        fMapGal2D->SetStats( 0 );
        fMapGal2D->SetXTitle( "Galactic longitude [deg]" );
        fMapGal2D->SetYTitle( "Galactic latitude [deg]" );
        
        fRadAccMapGal2D = new TH2D( "fRadAccMapGal2D", "", nBinsL * 2, -180., 180., nBinsB * 2, iBMin, iBMax );
        fRadAccMapGal2D->SetStats( 0 );
        fRadAccMapGal2D->SetXTitle( "Galactic longitude [deg]" );
        fRadAccMapGal2D->SetYTitle( "Galactic latitude [deg]" );
        
        fMapGal2D_aitoff = new TH2D( "fMapGal2D_aitoff", "", nBinsL, -180., 180., ( int )( nBinsB * 2. / 3. ), iBMin, iBMax );
        fMapGal2D_aitoff->SetStats( 0 );
        fMapGal2D_aitoff->SetXTitle( "Galactic longitude [deg]" );
        fMapGal2D_aitoff->SetYTitle( "Galactic latitude [deg]" );
        
        fRadAccMapGal2D_aitoff = new TH2D( "fRadAccMapGal2D_aitoff", "", nBinsL, -180., 180., ( int )( nBinsB * 2. / 3. ), iBMin, iBMax );
        fRadAccMapGal2D_aitoff->SetStats( 0 );
        fRadAccMapGal2D_aitoff->SetXTitle( "Galactic longitude [deg]" );
        fRadAccMapGal2D_aitoff->SetYTitle( "Galactic latitude [deg]" );
    }
    else
    {
        fMapGal2D = 0;
        fRadAccMapGal2D = 0;
        fMapGal2D_aitoff = 0;
        fRadAccMapGal2D_aitoff = 0;
    }
    
    fPlotVTSObjects = false;
    
    fDoCheckSums = false;
    
}


/*
   set more contours for the colz map
*/
void VExposure::set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
    gStyle->SetNumberContours( NCont );
}

void VExposure::setSelectLaser( int iSelectLaser )
{
    fSelectLaser = iSelectLaser;
}

void VExposure::setTimeRange( string iStart, string iStop )
{
    fStartDate_SQL = iStart;
    fStopDate_SQL  = iStop;
}

void VExposure::setSourceName( string sourceName )
{
    fTargetSourceName = sourceName;
}

void VExposure::setTelMinElevation( double iElevation )
{
    fTelMinElevation = iElevation;
}

void VExposure::setMinDuration( double iDuration )
{
    fMinDuration = iDuration;
}

bool VExposure::setPlannedObservation( vector<double> ra, vector<double> dec, vector<double> t )
{

    for( int unsigned i = 0; i < ra.size(); i++ )
    {
        double i_b = 0;
        double i_l = 0;
        
        VAstronometry::vlaEqgal( ra[i] / 180. * TMath::Pi(), dec[i] / 180. * TMath::Pi(), &i_l, &i_b );
        fRunGalLong1958.push_back( i_l * 180. / TMath::Pi() );
        fRunGalLat1958.push_back( i_b * 180. / TMath::Pi() );
        fRunDuration.push_back( t[i] );
        if( fDebug )
        {
            cout << " " << ra[i] << " " << dec[i] << " ";
        }
        if( fDebug )
        {
            cout << " " << fRunGalLong1958.back() << " " << fRunGalLat1958.back() << " " << flush;
        }
    }
    
    return true;
}

unsigned int VExposure::getLaserDate( unsigned int iRunNumber )
{

    string iTempS;
    iTempS = getDBServer() + "/VERITAS";
    
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "error connecting to db" << endl;
        return -1;
    }
    char c_query[1000];
    sprintf( c_query, "select * from tblRun_Info where run_id = %d", iRunNumber );
    if( !my_connection.make_query( c_query ) )
    {
        return -1;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    string itemp;
    TSQLRow* db_row = db_res->Next();
    if( !db_row )
    {
        return -1;
    }
    // get date
    if( db_row->GetField( 6 ) )
    {
        string iTemp = db_row->GetField( 6 );
        if( iTemp.size() > 8 )
        {
            fDataStartTime = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
        }
    }
    return fDataStartTime;
}

bool VExposure::readFromDB()
{
    string iTempS;
    iTempS = getDBServer() + "/VERITAS";
    
    //std::cout<<"VExposure::readFromDB() "<<std::endl;
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "error connecting to db" << endl;
        return false;
    }
    cout << "reading from " << iTempS << endl;
    char c_query[1000];
    
    if( !fMakeRunList )
    {
        sprintf( c_query, "select * from tblRun_Info where db_start_time >= \"%s.000000\"  and db_start_time < \"%s.160000\"", fStartDate_SQL.c_str(), fStopDate_SQL.c_str() );
    }
    else
    {
        sprintf( c_query, "select * from tblRun_Info where db_start_time >= \"%s.000000\"  and db_start_time < \"%s.160000\" and source_id like convert( _utf8 \'%s\' using latin1)", fStartDate_SQL.c_str(), fStopDate_SQL.c_str(), fTargetSourceName.c_str() );
    }
    cout << c_query << endl;
    
    if( !my_connection.make_query( c_query ) )
    {
        cout << "Returning false!" << endl;
        return false;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    int fNRows = db_res->GetRowCount();
    
    cout << "Number of runs read form db: " << fNRows << endl;
    
    string itemp;
    for( int j = 0; j < fNRows; j++ )
    {
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            break;
        }
        
        // all fields should be defined (check if field #19 is there)
        if( !db_row->GetField( 19 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 19 );
        // don't use laser or charge injection runs
        if( itemp == "NOSOURCE" )
        {
            continue;
        }
        
        if( !db_row->GetField( 1 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 1 );
        // check if this run is an observing run
        if( fObservingMode == "Normal" )
        {
            if( itemp != "observing" )
            {
                continue;
            }
        }
        // check if this run is an observing/redHV/UVFilter run
        else if( fObservingMode == "Special" )
        {
            if( itemp != "observing" && itemp != "obsFilter" && itemp != "obsLowHV" )
            {
                continue;
            }
        }
        
        fRunObsMode.push_back( itemp );
        
        // get run status
        if( !db_row->GetField( 3 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 3 );
        // don't use aborted runs
        if( itemp == "aborted" )
        {
            continue;
        }
        // don't use runs which are started only
        if( itemp == "started" )
        {
            continue;
        }
        // don't use runs which are defined only
        if( itemp == "defined" )
        {
            continue;
        }
        // don't use runs which are prepared only
        if( itemp == "prepared" )
        {
            continue;
        }
        
        // get date
        
        if( db_row->GetField( 6 ) )
        {
            string iTemp = db_row->GetField( 6 );
            if( iTemp.size() > 8 )
            {
                fDataStartTime = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
            }
        }
        
        if( fDebug )
        {
            cout << j << " " << db_row->GetField( 19 ) <<  " ";
        }
        if( fDebug )
        {
            cout << db_row->GetField( 1 ) << " " << db_row->GetField( 3 ) << " " << flush;
        }
        
        //////
        // get source coordinates
        double iRa = 0.;
        double iDec = 0.;
        itemp = db_row->GetField( 19 );
        if( !getDBSourceCoordinates( my_connection.Get_ConnectionResult(), itemp, iDec, iRa ) )
        {
            cout << "No coordinates available: " << db_row->GetField( 0 ) << " " << itemp << endl;
            continue;
        }
        
        fRunSourceID.push_back( itemp );
        fRunRA.push_back( iRa );
        fRunDec.push_back( iDec );
        
        fRunoffsetRA.push_back( atof( db_row->GetField( 15 ) ) * 180. / TMath::Pi() );
        fRunoffsetDec.push_back( atof( db_row->GetField( 16 ) ) * 180. / TMath::Pi() );
        fRunConfigMask.push_back( atoi( db_row->GetField( 10 ) ) );
        
        //////
        // get galactic coordinates
        double i_b = 0.;
        double i_l = 0.;
        iRa  += fRunoffsetRA.back();
        iDec += fRunoffsetDec.back();
        
        VAstronometry::vlaEqgal( iRa / 180. * TMath::Pi(), iDec / 180. * TMath::Pi(), &i_l, &i_b );
        fRunGalLong1958.push_back( i_l * 180. / TMath::Pi() );
        fRunGalLat1958.push_back( i_b * 180. / TMath::Pi() );
        if( fDebug )
        {
            cout << " " << iRa << " " << iDec << " ";
        }
        if( fDebug )
        {
            cout << " " << fRunGalLong1958.back() << " " << fRunGalLat1958.back() << " " << flush;
        }
        
        //////
        // calculate MJD, etc.
        int iMJD = 0;
        double iTime1 = 0.;
        double iTime2 = 0.;
        
        if( !db_row->GetField( 6 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 6 );
        getDBMJDTime( itemp, iMJD, iTime1, true );
        fRunStartMJD.push_back( ( double )iMJD + iTime1 / 86400 );
        if( fDebug )
        {
            cout << "T1 " << iMJD << " " << iTime1 << " ";
        }
        if( !db_row->GetField( 7 ) )
        {
            continue;
        }
        itemp = db_row->GetField( 7 );
        if( fDebug )
        {
            cout << itemp << " " << flush;
        }
        getDBMJDTime( itemp, iMJD, iTime2, true );
        if( fDebug )
        {
            cout << "T2 " << iMJD << " " << iTime2 << " ";
        }
        fRunStopMJD.push_back( ( double )iMJD + iTime2 / 86400 );
        fRunDuration.push_back( iTime2 - iTime1 );
        fRunDate.push_back( fDataStartTime );
        
        if( fDebug )
        {
            cout << fRunDuration.back() << " " << flush;
        }
        
        itemp = db_row->GetField( 0 );
        if( fDebug )
        {
            cout << "RUN " << itemp << " " << flush;
        }
        fRun.push_back( atoi( itemp.c_str() ) );
        
        //       itemp = db_row->GetField( 2 );
        fRunStatus.push_back( "" );
        
        // Wobble Stuff
        double fWobbleNorth = -9999.;
        double fWobbleEast = -9999.;
        
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
        if( fabs( angl ) < 0.1 )
        {
            fWobbleNorth = dist;
            fWobbleEast = 0.;
        }
        else if( fabs( angl - 90. ) < 0.1 )
        {
            fWobbleNorth = 0.;
            fWobbleEast = dist;
        }
        else if( fabs( angl - 180. ) < 0.1 )
        {
            fWobbleNorth = -1.*dist;
            fWobbleEast = 0.;
        }
        else if( fabs( angl - 270. ) < 0.1 )
        {
            fWobbleNorth = 0.;
            fWobbleEast = -1.*dist;
        }
        if( itemp  == "50308" )
        {
            fWobbleNorth = 0.;
            fWobbleEast = -0.5;
        }
        
        fWobbleN.push_back( fWobbleNorth );
        fWobbleE.push_back( fWobbleEast );
        
        //////
        // get local coordinates
        
        double ha = 0.;
        double iSid = 0.;
        double az, el;
        // get Greenwich sideral time
        iSid = VAstronometry::vlaGmsta( ( double )iMJD, ( iTime1 + iTime2 ) / 2. / 86400. );
        // calculate local sideral time
        iSid = iSid - VGlobalRunParameter::getObservatory_Longitude_deg() * TMath::DegToRad();
        // calculate hour angle
        ha = VAstronometry::vlaDranrm( iSid - iRa / 180. * TMath::Pi() );
        // get horizontal coordinates
        VAstronometry::vlaDe2h( ha, iDec / 180. * TMath::Pi(), VGlobalRunParameter::getObservatory_Latitude_deg() * TMath::DegToRad(), &az, &el );
        // fill vectors
        fRunTelElevation.push_back( el * 180. / TMath::Pi() );
        fRunTelAzimuth.push_back( az * 180. / TMath::Pi() );
        
        ////////////////////////////////////////////////////////////////////////
        if( fDebug )
        {
            cout << fRunStatus.back() << " " << fRun.back() << flush;
        }
        
        if( fDebug )
        {
            cout << endl;
        }
    }
    
    cout << "total number of runs in runlist: " << fRunStatus.size() << endl;
    
    return true;
}

bool VExposure::readFromDBList()
{
    string iTempS;
    iTempS = getDBServer() + "/VERITAS";
    
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        cout << "error connecting to db 2 " << iTempS << endl;
        return false;
    }
    cout << "Using external runlist. Reading from " << iTempS << endl;
    
    char c_query[1000];
    
    for( unsigned int i = 0; i < fRunDownloadList.size(); i++ )
    {
    
        sprintf( c_query, "select * from tblRun_Info where run_id=%d", fRunDownloadList[i] );
        if( fDebug )
        {
            cout << c_query << endl;
        }
        
        if( !my_connection.make_query( c_query ) )
        {
            cout << "VExposure::readFromDBList() error: could not connect to DB." << endl;
            return false;
        }
        TSQLResult* db_res =  my_connection.Get_QueryResult();
        
        string itemp;
        
        TSQLRow* db_row = db_res->Next();
        
        if( !db_row )
        {
        
            cout << "VExposure::readFromDBList() error: could not read run info for run " << fRunDownloadList[i] << endl;
            continue;
        }
        
        // get date
        
        if( db_row->GetField( 6 ) )
        {
            string iTemp = db_row->GetField( 6 );
            if( iTemp.size() > 8 )
            {
                fDataStartTime = atoi( iTemp.substr( 0, 4 ).c_str() ) * 10000 + atoi( iTemp.substr( 5, 2 ).c_str() ) * 100 + atoi( iTemp.substr( 8, 2 ).c_str() );
            }
            else
            {
                cout << "VExposure::readFromDBList() error: could not parse run date for run " <<  fRunDownloadList[i] << ". DB enstry for start time is " << iTemp << endl;
                continue;
            }
        }
        else
        {
            cout << "VExposure::readFromDBList() error: could not read start date for run " <<  fRunDownloadList[i] << endl;
            continue;
        }
        fRunDate.push_back( fDataStartTime );
        fRun.push_back( fRunDownloadList[i] );
        
        //////
        // get source coordinates
        double iRa = 0.;
        double iDec = 0.;
        itemp = "";
        
        if( db_row->GetField( 19 ) )
        {
            itemp = db_row->GetField( 19 );
            if( !getDBSourceCoordinates( my_connection.Get_ConnectionResult(), itemp, iDec, iRa ) )
            {
                cout << "No coordinates available: " << fRunDownloadList[i] << " " << itemp << endl;
            }
        }
        fRunSourceID.push_back( itemp );
        fRunRA.push_back( iRa );
        fRunDec.push_back( iDec );
        
        double ftemp = ( db_row->GetField( 15 ) != NULL  ? atof( db_row->GetField( 15 ) ) * 180. / TMath::Pi() : 0 );
        fRunoffsetRA.push_back( ftemp );
        
        ftemp = ( db_row->GetField( 16 ) != NULL ? 	 atof( db_row->GetField( 16 ) ) * 180. / TMath::Pi() : 0 );
        fRunoffsetDec.push_back( ftemp );
        
        ftemp = ( db_row->GetField( 10 ) != NULL ? 	atof( db_row->GetField( 10 ) ) * 180. / TMath::Pi() : 0 );
        fRunConfigMask.push_back( ftemp );
        
        //observing mode
        itemp = db_row->GetField( 1 ) ? db_row->GetField( 1 ) : "";
        fRunObsMode.push_back( itemp );
        
        //////
        // Wobble Stuff
        
        double fWobbleNorth = -9999.;
        double fWobbleEast = -9999.;
        
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
        if( fabs( angl ) < 0.1 )
        {
            fWobbleNorth = dist;
            fWobbleEast = 0.;
        }
        else if( fabs( angl - 90. ) < 0.1 )
        {
            fWobbleNorth = 0.;
            fWobbleEast = dist;
        }
        else if( fabs( angl - 180. ) < 0.1 )
        {
            fWobbleNorth = -1.*dist;
            fWobbleEast = 0.;
        }
        else if( fabs( angl - 270. ) < 0.1 )
        {
            fWobbleNorth = 0.;
            fWobbleEast = -1.*dist;
        }
        else if( TMath::Abs( dist ) < 1.e-2 )
        {
            fWobbleNorth = 0.;
            fWobbleEast  = 0.;
        }
        
        if( fRunDownloadList[i] == 50308 )
        {
            fWobbleNorth = 0.;
            fWobbleEast = -0.5;
        }
        
        fWobbleN.push_back( fWobbleNorth );
        fWobbleE.push_back( fWobbleEast );
        
        //////
        // get galactic coordinates
        double i_b = 0.;
        double i_l = 0.;
        iRa  += fRunoffsetRA.back();
        iDec += fRunoffsetDec.back();
        
        VAstronometry::vlaEqgal( iRa / 180. * TMath::Pi(), iDec / 180. * TMath::Pi(), &i_l, &i_b );
        fRunGalLong1958.push_back( i_l * 180. / TMath::Pi() );
        fRunGalLat1958.push_back( i_b * 180. / TMath::Pi() );
        if( fDebug )
        {
            cout << " " << iRa << " " << iDec << " ";
        }
        if( fDebug )
        {
            cout << " " << fRunGalLong1958.back() << " " << fRunGalLat1958.back() << " " << flush;
        }
        
        //////
        // calculate MJD, etc.
        int iMJD = 0;
        double iTime1 = 0.;
        double iTime2 = 0.;
        
        if( db_row->GetField( 6 ) )
        {
            itemp = db_row->GetField( 6 );
            getDBMJDTime( itemp, iMJD, iTime1, true );
            fRunStartMJD.push_back( ( double )iMJD + iTime1 / 86400 );
            if( fDebug )
            {
                cout << "T1 " << iMJD << " " << iTime1 << " ";
            }
        }
        else
        {
            fRunStartMJD.push_back( 0 );
        }
        
        
        if( db_row->GetField( 7 ) )
        {
            itemp = db_row->GetField( 7 );
            if( fDebug )
            {
                cout << itemp << " " << flush;
            }
            getDBMJDTime( itemp, iMJD, iTime2, true );
            if( fDebug )
            {
                cout << "T2 " << iMJD << " " << iTime2 << " ";
            }
            fRunStopMJD.push_back( ( double )iMJD + iTime2 / 86400 );
        }
        else
        {
            fRunStopMJD.push_back( 0 );
        }
        
        fRunDuration.push_back( iTime2 - iTime1 );
        
        
        if( fDebug )
        {
            cout << fRunDuration.back() << " " << flush;
        }
        
        
        //       itemp = db_row->GetField( 2 );
        fRunStatus.push_back( "" );
        
        //////
        // get local coordinates
        
        double ha = 0.;
        double iSid = 0.;
        double az, el;
        // get Greenwich sideral time
        iSid = VAstronometry::vlaGmsta( ( double )iMJD, ( iTime1 + iTime2 ) / 2. / 86400. );
        // calculate local sideral time
        iSid = iSid - VGlobalRunParameter::getObservatory_Longitude_deg() * TMath::DegToRad();
        // calculate hour angle
        ha = VAstronometry::vlaDranrm( iSid - iRa / 180. * TMath::Pi() );
        // get horizontal coordinates
        VAstronometry::vlaDe2h( ha, iDec / 180. * TMath::Pi(), VGlobalRunParameter::getObservatory_Latitude_deg() * TMath::DegToRad(), &az, &el );
        // fill vectors
        fRunTelElevation.push_back( el * 180. / TMath::Pi() );
        fRunTelAzimuth.push_back( az * 180. / TMath::Pi() );
        
        ////////////////////////////////////////////////////////////////////////
        if( fDebug )
        {
            cout << fRunStatus.back() << " " << fRun.back() << flush;
        }
        
        if( fDebug )
        {
            cout << endl;
        }
        
    }
    
    
    cout << "total number of runs in runlist: " << fRun.size() << endl;
    
    return true;
}


void VExposure::resetDataVectors()
{
    fRun.clear();
    fRunConfigMask.clear();
    fRunStartMJD.clear();
    fRunStopMJD.clear();
    fRunDuration.clear();
    fRunRA.clear();
    fRunDec.clear();
    fRunoffsetRA.clear();
    fRunoffsetDec.clear();
    fRunGalLong1958.clear();
    fRunGalLat1958.clear();
    fRunTelElevation.clear();
    fRunTelAzimuth.clear();
}


bool VExposure::readRootFile( string iRFile, double iMinMJD, double iMaxMJD )
{
    resetDataVectors();
    
    TFile f( iRFile.c_str() );
    if( f.IsZombie() )
    {
        return false;
    }
    
    int iRun = 0;
    Char_t iSourceID[300];
    unsigned int iConfigMask = 0;
    double iRunStartMJD = 0.;
    double iRunStoppMJD = 0.;
    double iRunDuration = 0.;
    double iRunRA = 0.;
    double iRunDec = 0.;
    double iRunRAoffset = 0.;
    double iRunDecoffset = 0.;
    double l = 0.;
    double b = 0.;
    double el = 0.;
    double az = 0.;
    
    TTree* t = ( TTree* )f.Get( "fRunTable" );
    if( !t )
    {
        return false;
    }
    
    t->SetBranchAddress( "runNumber", &iRun );
    t->SetBranchAddress( "sourceID", &iSourceID );
    t->SetBranchAddress( "configmask", &iConfigMask );
    t->SetBranchAddress( "MJDstart", &iRunStartMJD );
    t->SetBranchAddress( "MJDstop", &iRunStoppMJD );
    t->SetBranchAddress( "runLength", &iRunDuration );
    t->SetBranchAddress( "ra", &iRunRA );
    t->SetBranchAddress( "dec", &iRunDec );
    t->SetBranchAddress( "raOffset", &iRunRAoffset );
    t->SetBranchAddress( "decOffset", &iRunDecoffset );
    t->SetBranchAddress( "l", &l );
    t->SetBranchAddress( "b", &b );
    t->SetBranchAddress( "telEl", &el );
    t->SetBranchAddress( "telAz", &az );
    
    double i_tot = 0.;
    
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        
        if( iMinMJD > 0. && iRunStartMJD < iMinMJD )
        {
            continue;
        }
        if( iMaxMJD > 0. && iRunStartMJD > iMaxMJD )
        {
            continue;
        }
        
        fRun.push_back( iRun );
        fRunConfigMask.push_back( iConfigMask );
        fRunStartMJD.push_back( iRunStartMJD );
        fRunStopMJD.push_back( iRunStoppMJD );
        fRunDuration.push_back( iRunDuration );
        i_tot += iRunDuration;
        fRunSourceID.push_back( iSourceID );
        fRunRA.push_back( iRunRA );
        fRunDec.push_back( iRunDec );
        fRunoffsetRA.push_back( iRunRAoffset );
        fRunoffsetDec.push_back( iRunDecoffset );
        if( l > 180. )
        {
            l -= 360.;
        }
        fRunGalLong1958.push_back( l );
        fRunGalLat1958.push_back( b );
        fRunTelElevation.push_back( el );
        fRunTelAzimuth.push_back( az );
        
    }
    
    f.Close();
    
    cout << "total number of entries read in: " << fRunGalLat1958.size() << endl;
    cout << "total observing time: " << i_tot / 3600. << " [h]" << endl;
    
    return true;
}


bool VExposure::writeRootFile( string iOfile )
{
    if( iOfile.size() <= 0 )
    {
        return false;
    }
    
    int iRun;
    Char_t iRunSourceID[300];
    unsigned int iConfigMask;
    double iRunStartMJD;
    double iRunStoppMJD;
    double iRunDuration;
    double iRunRA;
    double iRunDec;
    double iRunRAoffset;
    double iRunDecoffset;
    double l;
    double b;
    double el;
    double az;
    
    TFile f( iOfile.c_str(), "RECREATE" );
    if( f.IsZombie() )
    {
        return false;
    }
    cout << "writing data to " << f.GetName() << endl;
    
    char hname[2000];
    sprintf( hname, "DB entries, %s to %s", fStartDate_SQL.c_str(), fStopDate_SQL.c_str() );
    TTree t( "fRunTable", hname );
    t.Branch( "runNumber", &iRun, "runNumber/I" );
    t.Branch( "sourceID", &iRunSourceID, "sourceID/C" );
    t.Branch( "configmask", &iConfigMask, "configmask/i" );
    t.Branch( "MJDstart", &iRunStartMJD, "MJDstart/D" );
    t.Branch( "MJDstop", &iRunStoppMJD, "MJDstop/D" );
    t.Branch( "runLength", &iRunDuration, "runLength/D" );
    t.Branch( "ra", &iRunRA, "ra/D" );
    t.Branch( "dec", &iRunDec, "dec/D" );
    t.Branch( "raOffset", &iRunRAoffset, "raOffset/D" );
    t.Branch( "decOffset", &iRunDecoffset, "decOffset/D" );
    t.Branch( "l", &l, "l/D" );
    t.Branch( "b", &b, "b/D" );
    t.Branch( "telEl", &el, "telEl/D" );
    t.Branch( "telAz", &az, "telAz/D" );
    
    for( unsigned int i = 0; i < fRun.size(); i++ )
    {
        iRun = fRun[i];
        if( fRunSourceID[i].size() < 300 )
        {
            sprintf( iRunSourceID, "%s", fRunSourceID[i].c_str() );
        }
        else
        {
            sprintf( iRunSourceID, "%s", fRunSourceID[i].substr( 0, 299 ).c_str() );
        }
        iConfigMask = fRunConfigMask[i];
        iRunStartMJD = fRunStartMJD[i];
        iRunStoppMJD = fRunStopMJD[i];
        iRunDuration = fRunDuration[i];
        iRunRA = fRunRA[i];
        iRunDec = fRunDec[i];
        iRunRAoffset = fRunoffsetRA[i];
        iRunDecoffset = fRunoffsetDec[i];
        l = fRunGalLong1958[i];
        b = fRunGalLat1958[i];
        el = fRunTelElevation[i];
        az = fRunTelAzimuth[i];
        
        t.Fill();
    }
    t.Write();
    
    f.Close();
    
    return true;
}

/*
    fill elevation vs galactic coordinates for the given observatory

    (note: preliminary, some numbers are hardwired)

*/
void VExposure::fillElevationPlot( int iYear, int iMonth, int ze_max_deg )
{
    fMapGal2D->Reset();
    fMapGal2D_aitoff->Reset();
    
    vector< int > i_MJD;
    if( iMonth < 0 )
    {
        for( unsigned int i = 0; i < 12; i++ )
        {
            if( i == 6 || i == 7 )
            {
                continue;
            }
            i_MJD.push_back( VSkyCoordinatesUtilities::getMJD( iYear, i + 1, 1 ) );
        }
    }
    else
    {
        i_MJD.push_back( VSkyCoordinatesUtilities::getMJD( iYear, iMonth, 1 ) );
    }
    
    double dec = 0.;
    double ra = 0.;
    double l = 0.;
    double b = 0.;
    
    cout << "Observatory coordinates: " << getObservatory_Longitude_deg() << "\t" << getObservatory_Latitude_deg() << endl;
    
    for( unsigned int i = 0; i < i_MJD.size(); i++ )
    {
        // observations from 7 pm to 6 am local time
        // (hardwired AZ times)
        for( int t = 120; t < 520; t++ )
        {
            double imjd = i_MJD[i];
            double time = 86400. / 1000.*( double )t;
            // zenith angle range
            for( int z = ze_max_deg; z > 0; z-- )
            {
                // azimuth range
                for( int a = 0; a < 360; a++ )
                {
                    VSkyCoordinatesUtilities::getEquatorialCoordinates( imjd, time, ( double )a, ( double )z, dec, ra );
                    VAstronometry::vlaEqgal( ra * TMath::DegToRad(), dec * TMath::DegToRad(), &l, &b );
                    
                    l *= TMath::RadToDeg();
                    if( l > 180 )
                    {
                        l -= 360.;
                    }
                    b *= TMath::RadToDeg();
                    
                    int i_l = fMapGal2D->GetXaxis()->FindBin( -1.*l );
                    int i_b = fMapGal2D->GetYaxis()->FindBin( b );
                    if( 90. - ( double )z > fMapGal2D->GetBinContent( i_l, i_b ) )
                    {
                        fMapGal2D->SetBinContent( i_l, i_b, 90. - ( double )z );
                    }
                }
            }
        }
    }
    
    /////////////////////////////////
    // calculate aitoff projection
    double al = 0.;
    double ab = 0.;
    double xl = 0.;
    double xb = 0.;
    int bl = 0;
    int bb = 0;
    
    for( int i = 1; i <= fMapGal2D->GetNbinsX(); i++ )
    {
        xl = fMapGal2D->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= fMapGal2D->GetNbinsY(); j++ )
        {
            xb = fMapGal2D->GetYaxis()->GetBinCenter( j );
            
            aitoff2xy( -1.*xl, xb, al, ab );
            
            bl = fMapGal2D_aitoff->GetXaxis()->FindBin( al );
            bb = fMapGal2D_aitoff->GetYaxis()->FindBin( ab );
            
            fMapGal2D_aitoff->SetBinContent( bl, bb, fMapGal2D->GetBinContent( i, j ) );
        }
    }
    gStyle->SetPalette( 1 );
    set_plot_style();
    
    bPlotElevationPlots = true;
}


void VExposure::fillExposureMap()
{
    cout << "fill exposure map (may take a while)" << endl;
    
    fMapGal2D->Reset();
    fRadAccMapGal2D->Reset();
    fMapGal2D_aitoff->Reset();
    fRadAccMapGal2D_aitoff->Reset();
    
    int    i_r_b = ( int )( fMaximumIntegrationRadius / fMapGal2D->GetYaxis()->GetBinWidth( 2 ) + 0.5 );
    
    double r_dist = 0.;
    
    //////////////////////////////////////////////////////////////////////
    // loop over all runs and fill the map
    for( unsigned int i = 0; i < fRunGalLong1958.size(); i++ )
    {
        double l = fRunGalLong1958[i];
        if( l > 180. )
        {
            l -= 360.;
        }
        int i_l = fMapGal2D->GetXaxis()->FindBin( l );
        int i_b = fMapGal2D->GetYaxis()->FindBin( fRunGalLat1958[i] );
        
        int i_b_start = i_b - i_r_b;
        if( i_b_start < 1 )
        {
            i_b_start = 1;
        }
        int i_b_stop = i_b + i_r_b;
        if( i_b_stop > fMapGal2D->GetNbinsY() )
        {
            i_b_stop = fMapGal2D->GetNbinsY();
        }
        
        for( int b = i_b_start; b <= i_b_stop; b++ )
        {
            double b_pos = fMapGal2D->GetYaxis()->GetBinCenter( b );
            
            // calculate extension in longitude
            int i_r_l = 1;
            if( cos( b_pos * TMath::Pi() / 180. ) > 0. )
            {
                i_r_l = ( int )( fMaximumIntegrationRadius / cos( b_pos * TMath::Pi() / 180. ) / fMapGal2D->GetXaxis()->GetBinWidth( 2 ) + 0.5 );
            }
            else
            {
                i_r_l = 0;
            }
            int i_l_start = i_l - i_r_l;
            if( i_l_start < 1 )
            {
                i_l_start = 1;
            }
            int i_l_stop = i_l + i_r_l;
            if( i_l_stop > fMapGal2D->GetNbinsX() )
            {
                i_l_stop = fMapGal2D->GetNbinsX();
            }
            
            for( int l = i_l_start; l <= i_l_stop; l++ )
            {
                double l_pos = fMapGal2D->GetXaxis()->GetBinCenter( l );
                
                r_dist = VAstronometry::vlaDsep( l_pos * TMath::Pi() / 180., b_pos * TMath::Pi() / 180., fRunGalLong1958[i] * TMath::Pi() / 180.,
                                  fRunGalLat1958[i] * TMath::Pi() / 180. ) * 180. / TMath::Pi();
                if( r_dist < fMaximumIntegrationRadius && fRunDuration[i] > 0. )
                {
                    // galactic longitudes are from 180. to -180.
                    l_pos *= -1.;
                    fMapGal2D->Fill( l_pos, b_pos, fRunDuration[i] / 3600. );
                    fRadAccMapGal2D->Fill( l_pos, b_pos, fRunDuration[i] / 3600. * getAcceptance( r_dist ) );
                }
            }
        }
    }
    cout << "entries " << fMapGal2D->GetEntries() << endl;
    
    /////////////////////////////////
    // calculate aitoff projection
    double al = 0.;
    double ab = 0.;
    double xl = 0.;
    double xb = 0.;
    int bl = 0;
    int bb = 0;
    
    for( int i = 1; i <= fMapGal2D->GetNbinsX(); i++ )
    {
        xl = fMapGal2D->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= fMapGal2D->GetNbinsY(); j++ )
        {
            xb = fMapGal2D->GetYaxis()->GetBinCenter( j );
            
            aitoff2xy( -1.*xl, xb, al, ab );
            
            bl = fMapGal2D_aitoff->GetXaxis()->FindBin( al );
            bb = fMapGal2D_aitoff->GetYaxis()->FindBin( ab );
            
            fMapGal2D_aitoff->SetBinContent( bl, bb, fMapGal2D->GetBinContent( i, j ) );
            fRadAccMapGal2D_aitoff->SetBinContent( bl, bb, fRadAccMapGal2D->GetBinContent( i, j ) );
        }
    }
    // now plot everything
    gStyle->SetPalette( 1 );
    set_plot_style();
    
}


TCanvas* VExposure::plot( double ilmin, double ilmax, double ibmin, double ibmax, unsigned int iReturnCanvas, int iPalette, TString opt )
{
    fPlottingCanvas.clear();
    fPlottingCanvas.push_back( plot2DGalactic( "cGal", "Exposure", 10, 10, fCanvasSize_x, fCanvasSize_y, fMapGal2D, ibmin, ibmax, ilmin, ilmax, false, iPalette, opt ) );
    fPlottingCanvas.push_back( plot2DGalactic( "cGal_aitoff", "Exposure (aitoff)", 10, 450, fCanvasSize_x, fCanvasSize_y, fMapGal2D_aitoff, ibmin, ibmax, ilmin, ilmax, true, iPalette, opt ) );
    if( fAcceptance )
    {
        fPlottingCanvas.push_back( plot2DGalactic( "cAccGal", "Acceptance corrected exposure", 650, 10, 600, 400,
                                   fRadAccMapGal2D, ibmin, ibmax, ilmin, ilmax, false, iPalette, opt ) );
        fPlottingCanvas.push_back( plot2DGalactic( "cAccGal_aitoff", "Acceptance corrected exposure (aitoff)", 650, 450, 600, 400,
                                   fRadAccMapGal2D_aitoff, ibmin, ibmax, ilmin, ilmax, true, iPalette, opt ) );
    }
    
    if( iReturnCanvas < fPlottingCanvas.size() )
    {
        return fPlottingCanvas[iReturnCanvas];
    }
    
    return 0;
}

/*
 * plots extension of HESS survey
 *
 * (note: this is the original HESS survey, not the extended one )
 *
 */
void VExposure::plot_HESSSkySurvey( TCanvas* c )
{
    if( !c )
    {
        return;
    }
    
    c->cd();
    
    // first phase
    TBox* ib1 = new TBox( -30., -3., 30., 3 );
    ib1->SetFillStyle( 0 );
    ib1->SetLineColor( 14 );
    ib1->SetLineStyle( 2 );
    ib1->Draw();
    TBox* ib2 = new TBox( 30., -3., 80., 3 );
    ib2->SetFillStyle( 0 );
    ib2->SetLineColor( 14 );
    ib2->SetLineStyle( 2 );
    ib2->Draw();
    TBox* ib3 = new TBox( -60., -3., -30., 3 );
    ib3->SetFillStyle( 0 );
    ib3->SetLineColor( 14 );
    ib3->SetLineStyle( 2 );
    ib3->Draw();
}


TCanvas* VExposure::plot2DGalactic( string iName, string iTitle, int ix, int iy, int iwx, int iwy, TH2* h,
                                    double ibmin, double ibmax, double ilmin, double ilmax, bool bAitoff, int iPalette, TString opt )
{
    if( !h )
    {
        return 0;
    }
    
    // canvas
    TCanvas* cGal = new TCanvas( iName.c_str(), iTitle.c_str(), ix, iy, iwx, iwy );
    cGal->SetGridx( 0 );
    cGal->SetGridy( 0 );
    cGal->SetRightMargin( 0.15 );
    cGal->Draw();
    
    // set axis ranges of galactic map
    h->SetAxisRange( ibmin, ibmax, "Y" );
    h->SetAxisRange( -1.*ilmax, -1.*ilmin, "X" );
    
    if( bPlotElevationPlots )
    {
        cGal->SetGridx( 1 );
        cGal->SetGridy( 1 );
        h->SetZTitle( "elevation [deg]" );
    }
    else
    {
        cGal->SetGridx( 0 );
        cGal->SetGridy( 0 );
        h->SetZTitle( "exposure [h]" );
    }
    gStyle->SetPalette( iPalette );
    h->Draw( opt );
    
    // plot axis
    
    TF1* IncValues = new TF1( "IncValues", "-x", ilmin, ilmax );
    IncValues->SetNpx( 500 );
    
    TGaxis* decAxis = new TGaxis( cGal->GetUxmin() - ilmax - h->GetXaxis()->GetBinWidth( 1 ) / 2., cGal->GetUymin() + ibmin - h->GetXaxis()->GetBinWidth( 1 ) / 2., cGal->GetUxmin() - ilmax - h->GetXaxis()->GetBinWidth( 1 ) / 2., cGal->GetUymin() + ibmax + h->GetXaxis()->GetBinWidth( 1 ) / 2., ibmin, ibmax, 505, "-" );
    decAxis->SetTitleFont( h->GetXaxis()->GetTitleFont() );
    decAxis->SetTitleSize( h->GetXaxis()->GetTitleSize() );
    decAxis->SetTitleOffset( h->GetXaxis()->GetTitleOffset() );
    decAxis->SetLabelFont( h->GetXaxis()->GetLabelFont() );
    decAxis->SetLabelSize( h->GetXaxis()->GetLabelSize() );
    decAxis->SetLabelColor( h->GetXaxis()->GetLabelColor() );
    decAxis->SetLabelOffset( 0.020 );
    decAxis->SetTextColor( h->GetXaxis()->GetTitleColor() );
    decAxis->SetLineColor( h->GetXaxis()->GetAxisColor() );
    decAxis->SetTitle( "galactic longitude [deg]" );
    //  decAxis->SetNdivisions( 510 );
    decAxis->SetTitle( "galactic latitude [deg]" );
    decAxis->Draw();
    
    //   TGaxis *raLowerAxis = new TGaxis( cGal->GetUxmin()-ilmax - h->GetYaxis()->GetBinWidth( 1 )/2., cGal->GetUymin()+ibmin-h->GetYaxis()->GetBinWidth(1)/2. , cGal->GetUxmin()-ilmin+h->GetYaxis()->GetBinWidth( 1 )/2., cGal->GetUymin()+ibmin -h->GetYaxis()->GetBinWidth(1)/2. , "IncValues", 4 );
    TGaxis* raLowerAxis = new TGaxis( cGal->GetUxmin() - ilmax - h->GetYaxis()->GetBinWidth( 1 ) / 2., cGal->GetUymin() + ibmin - h->GetYaxis()->GetBinWidth( 1 ) / 2. , cGal->GetUxmin() - ilmin + h->GetYaxis()->GetBinWidth( 1 ) / 2., cGal->GetUymin() + ibmin - h->GetYaxis()->GetBinWidth( 1 ) / 2. , "IncValues", 4 );
    raLowerAxis->SetTitleFont( h->GetXaxis()->GetTitleFont() );
    raLowerAxis->SetTitleSize( h->GetXaxis()->GetTitleSize() );
    raLowerAxis->SetTitleOffset( h->GetXaxis()->GetTitleOffset() );
    raLowerAxis->SetLabelFont( h->GetXaxis()->GetLabelFont() );
    raLowerAxis->SetLabelSize( h->GetXaxis()->GetLabelSize() );
    raLowerAxis->SetLabelColor( h->GetXaxis()->GetLabelColor() );
    raLowerAxis->SetLabelOffset( 0.011 );
    raLowerAxis->SetTextColor( h->GetXaxis()->GetTitleColor() );
    raLowerAxis->SetLineColor( h->GetXaxis()->GetAxisColor() );
    raLowerAxis->SetTitle( "galactic longitude [deg]" );
    //   raLowerAxis->SetNdivisions( 510 );
    raLowerAxis->Draw();
    
    for( unsigned int t = 0; t < fCatalogue.size(); t++ )
    {
        if( fCatalogue[t].size() > 0 )
        {
            analyseCatalogue( fCatalogue[t], ibmin, ibmax, ilmin, ilmax, h, bAitoff, fCatalogueMarkerStyle[t], fCatalogueMarkerColor[t],
                              fCatalogueTextAngle[t] );
        }
    }
    if( fPlotVTSObjects )
    {
        plotVTSObjects( bAitoff, ibmin, ibmax, ilmin, ilmax, 5, 1, 45., h );
    }
    
    // draw aitoff coordinate system
    if( bAitoff )
    {
        drawAitoffCoordinateSystem();
    }
    
    cGal->Update();
    
    return cGal;
}


/*
   copied from ROOT forum #6650
*/
void VExposure::drawAitoffCoordinateSystem()
{
    double radeg = TMath::Pi() / 180.;
    double la, lo, x, y, z;
    
    const int Nl = 11;                            // Number of drawn latitudes
    const int NL = 10;                            // Number of drawn longitudes
    int       M  = 60;
    
    TGraph*  latitudes[Nl];
    TGraph*  longitudes[NL];
    
    for( int j = 0; j < Nl; ++j )
    {
        latitudes[j] = new TGraph();
        
        la = -90. + 180. / ( Nl - 1 ) * j;
        for( int i = 0; i < M + 1; ++i )
        {
            lo = -180. + 360. / M * i;
            z  = sqrt( 1 + cos( la * radeg ) * cos( lo * radeg / 2. ) );
            x  = 180.*cos( la * radeg ) * sin( lo * radeg / 2. ) / z;
            y  = 90.*sin( la * radeg ) / z;
            latitudes[j]->SetPoint( i, x, y );
        }
    }
    for( int j = 0; j < Nl; ++j )
    {
        latitudes[j]->Draw( "l" );
    }
    
    for( int j = 0; j < NL; ++j )
    {
        longitudes[j] = new TGraph();
        lo = -180. + 360. / ( NL - 1 ) * j;
        for( int i = 0; i < M + 1; ++i )
        {
            la = -90. + 180. / M * i;
            z  = sqrt( 1 + cos( la * radeg ) * cos( lo * radeg / 2. ) );
            x  = 180.*cos( la * radeg ) * sin( lo * radeg / 2. ) / z;
            y  = 90.*sin( la * radeg ) / z;
            longitudes[j]->SetPoint( i, x, y );
        }
    }
    
    for( int j = 0; j < NL; ++j )
    {
        longitudes[j]->Draw( "l" );
    }
}

void VExposure::plotVTSObjects( bool bAitoff, double ibmin, double ibmax, double ilmin, double ilmax,
                                int iMarkerStyle, int iMarkerColor, double iTextAngle, TH2* h )
{
    double l = 0.;
    double b = 0.;
    double ra = 0.;
    double dec = 0.;
    string l_name = "";
    
    set< string > iObjects;
    for( unsigned int i = 0; i < fRunSourceID.size(); i++ )
    {
        if( iObjects.find( fRunSourceID[i] ) == iObjects.end() )
        {
            iObjects.insert( fRunSourceID[i] );
            
            ra  = fRunRA[i];
            dec = fRunDec[i];
            VAstronometry::vlaEqgal( ra / 180. * TMath::Pi(), dec / 180. * TMath::Pi(), &l, &b );
            
            plotObject( l * TMath::RadToDeg(), b * TMath::RadToDeg(), fRunSourceID[i], 0.,
                        ibmin, ibmax, ilmin, ilmax, h,
                        bAitoff, iMarkerStyle, iMarkerColor, iTextAngle );
        }
    }
}


void VExposure::analyseCatalogue( string iCatalogue, double ibmin, double ibmax, double ilmin, double ilmax, TH2* h,
                                  bool bAitoff, int iMarkerStyle, int iMarkerColor, double iTextAngle )
{
    VStarCatalogue* s = new VStarCatalogue();
    s->init( 54626., iCatalogue );
    //    s->printCatalogue();
    
    double l = 0.;
    double b = 0.;
    string l_name = "";
    
    
    cout << "total number of objects in catalogue: " << s->getNStar() << endl;
    for( unsigned int i = 0; i < s->getNStar(); i++ )
    {
        l = s->getStar_l( i );
        b = s->getStar_b( i );
        l_name = s->getStarName( i );
        
        if( iCatalogue == "VERITASDB" && l_name.substr( 0, 2 ) == "SS" )
        {
            continue;
        }
        if( iCatalogue == "VERITASDB" && l_name.substr( 0, 3 ) == "GRB" )
        {
            continue;
        }
        
        plotObject( l, b, l_name, s->getStarMajorDiameter( i ), ibmin, ibmax, ilmin, ilmax, h,
                    bAitoff, iMarkerStyle, iMarkerColor, iTextAngle );
    }
}

void VExposure::plotObject( double l, double b, string l_name, double iExtension,
                            double ibmin, double ibmax, double ilmin, double ilmax, TH2* h,
                            bool bAitoff, int iMarkerStyle, int iMarkerColor, double iTextAngle )
{
    double ib_range = ibmax - ibmin;
    double il_range = ilmax - ilmin;
    
    if( l > 180. )
    {
        l -= 360.;
    }
    if( l > ilmin && l < ilmax && b > ibmin && b < ibmax )
    {
        // aitoff projections?
        double ab = b;
        double al = -1.*l;
        if( bAitoff )
        {
            aitoff2xy( l, b, al, ab );
        }
        // extended?
        if( iExtension > 0. && fPlotExtendedSources )
        {
            TEllipse* m = new TEllipse( al, ab, iExtension, iExtension );
            m->SetFillStyle( 0 );
            m->SetLineColor( iMarkerColor );
            m->Draw();
        }
        else
        {
            TMarker* m = new TMarker( al, ab, 5 );
            m->SetMarkerColor( iMarkerColor );
            m->SetMarkerStyle( iMarkerStyle );
            m->SetMarkerSize( 2 );
            m->Draw();
        }
        
        if( h )
        {
            if( bAitoff )
            {
                aitoff2xy( l, b, al, ab );
            }
            else
            {
                al = -1.*l;
                ab = b;
            }
            int il = h->GetXaxis()->FindBin( al );
            int ib = h->GetYaxis()->FindBin( ab );
            if( h->GetBinContent( il, ib ) > 0. )
            {
                cout << h->GetBinContent( il, ib ) * 60.;
                cout << "  [min] exposure on object " << l_name;
                cout << "(l,b)=(" << l << ", " << b << ")\t";
                cout << endl;
            }
            if( fPlotSourceNames )
            {
                if( bAitoff )
                {
                    al = al + il_range * 0.01 * ( cos( iTextAngle * TMath::DegToRad() ) + sin( iTextAngle * TMath::DegToRad() ) );
                    ab = ab + ib_range * 0.01 * ( -1.*sin( iTextAngle * TMath::DegToRad() ) + cos( iTextAngle * TMath::DegToRad() ) );
                }
                else
                {
                    al = -1.*l + il_range * 0.01 * ( cos( iTextAngle * TMath::DegToRad() ) + sin( iTextAngle * TMath::DegToRad() ) );
                    ab = b + ib_range * 0.01 * ( -1.*sin( iTextAngle * TMath::DegToRad() ) + cos( iTextAngle * TMath::DegToRad() ) );
                }
                TText* t = new TText( al, ab, l_name.c_str() );
                t->SetTextColor( iMarkerColor );
                t->SetTextSize( t->GetTextSize() * 0.9 );
                t->SetTextAngle( iTextAngle );
                t->Draw();
            }
        }
    }
}


void VExposure::getDBMJDTime( string itemp, int& MJD, double& Time, bool bStrip )
{
    if( itemp.size() < 16 )
    {
        MJD = 0;
        Time = 0.;
        return;
    }
    if( bStrip )
    {
        itemp.replace( itemp.find( "-" ), 1, "" );
        itemp.replace( itemp.find( "-" ), 1, "" );
        itemp.replace( itemp.find( " " ), 1, "" );
        itemp.replace( itemp.find( ":" ), 1, "" );
        itemp.replace( itemp.find( ":" ), 1, "" );
    }
    int y, m, d, h, min, s, ms, l;
    double gMJD;
    // get y, m, d
    y = atoi( itemp.substr( 0, 4 ).c_str() );
    m = atoi( itemp.substr( 4, 2 ).c_str() );
    d = atoi( itemp.substr( 6, 2 ).c_str() );
    h = atoi( itemp.substr( 8, 2 ).c_str() );
    min = atoi( itemp.substr( 10, 2 ).c_str() );
    s = atoi( itemp.substr( 12, 2 ).c_str() );
    if( !bStrip )
    {
        ms = atoi( itemp.substr( 14, 3 ).c_str() );
    }
    else
    {
        ms = 0;
    }
    // calculate MJD
    VAstronometry::vlaCldj( y, m, d, &gMJD, &l );
    MJD = ( int )gMJD;
    Time = h * 60.*60. + min * 60. + s + ms / 1.e3;
}


bool VExposure::getDBSourceCoordinates( TSQLServer* f_db, string iSource, double& iEVNTargetDec, double& iEVNTargetRA )
{
    if( !f_db )
    {
        return false;
    }
    
    if( !fVDB_ObservingSources )
    {
        fVDB_ObservingSources = new VDB_ObservingSources();
        if( !fVDB_ObservingSources->fill( f_db ) )
        {
            return false;
        }
    }
    
    VDB_ObservingSources_Data* iTemp = fVDB_ObservingSources->get_ObservingSources_Data( iSource );
    if( iTemp )
    {
        iEVNTargetDec = iTemp->fDec;
        iEVNTargetRA  = iTemp->fRA;
        return true;
    }
    
    iEVNTargetDec = 0.;
    iEVNTargetRA = 0.;
    return false;
}


double VExposure::getAcceptance( double r )
{
    if( r < 0. || r > fAcceptance_MaxDistance )
    {
        return 0.;
    }
    
    double iacc = 1.;
    if( fAcceptance )
    {
        iacc = fAcceptance->Eval( r );
    }
    if( iacc < 0. )
    {
        iacc = 0.;
    }
    
    return iacc;
}


bool VExposure::readAcceptanceCurveFromFile( string iFile, double iAcceptance_MaxDistance )
{
    if( iFile.size() <= 0 )
    {
        return false;
    }
    
    cout << "reading acceptance curves from " << iFile << endl;
    TFile* fAcceptanceFile = new TFile( iFile.c_str() );
    if( fAcceptanceFile->IsZombie() )
    {
        return false;
    }
    
    fAcceptance = ( TF1* )fAcceptanceFile->Get( "fAccZe_0" );
    
    if( !fAcceptance )
    {
        cout << "Warning: acceptance curve not found in file: " << iFile << endl;
        return false;
    }
    
    fAcceptance_MaxDistance = iAcceptance_MaxDistance;
    
    return true;
}

/*
    apply DQM to a single run


    (still some hardwired numbers)
*/
bool VExposure::doDQM( unsigned int iID, double iMinDuration )
{
    //   if( iID >= (unsigned int)fRun.size() ) return false;
    
    // check run duration
    if( fRunDuration[iID] < iMinDuration )
    {
        return false;
    }
    
    // HARDWIRED: require at least three telescopes
    if( fRunConfigMask[iID] == 0 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 1 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 2 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 3 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 4 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 5 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 6 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 8 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 9 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 10 )
    {
        return false;
    }
    if( fRunConfigMask[iID] == 12 )
    {
        return false;
    }
    
    return true;
}


void VExposure::printListOfRuns( string iCatalogue, double iR, double iMinDuration, string iTeVCatalogue, double r_min, string iEventListFile )
{
    VStarCatalogue* s = new VStarCatalogue();
    s->init( 54626., iCatalogue );
    
    VStarCatalogue* tev = new VStarCatalogue();
    if( iTeVCatalogue.size() > 0 )
    {
        tev->init( 54626., iTeVCatalogue );
    }
    
    // file for event list
    TFile* fEventListfile = 0;
    TEventList* fEventList_tevcat = 0;
    TEventList* fEventList_inFOV = 0;
    TEventList* fEventList_inFOV_noTevcat = 0;
    if( iEventListFile.size() > 0 )
    {
        fEventListfile = new TFile( iEventListFile.c_str(), "RECREATE" );
        if( fEventListfile->IsZombie() )
        {
            return;
        }
        fEventList_tevcat = new TEventList( "tevcatList" );
        fEventList_inFOV = new TEventList( "veritas" );
        fEventList_inFOV_noTevcat = new TEventList( "V_noPointing" );
    }
    
    set< string > iVERITAS_targets;
    
    double r_centre = 0.;
    double r_VA_object = 0.;
    
    ////////////////////////////////////////////////////////
    // loop over all objects in catalogue
    cout << "total number of objects in catalogue: " << s->getNStar() << endl;
    for( unsigned int i = 0; i < s->getNStar(); i++ )
    {
        double r_tot = 0.;
        double r_tot_V5 = 0.;
        
        double r_centre_mean = 0.;
        double r_VA_object_mean = 0.;
        double r_N = 0.;
        double r_meanElevation = 0.;
        
        cout << "Testing now ";
        cout << s->getStarName( i );
        cout << " (l,b) = " << s->getStar_l( i ) << ", " << s->getStar_b( i );
        cout << " (ra,dec) = " << s->getStarRA2000( i ) << ", " << s->getStarDec2000( i );
        cout << endl;
        
        
        ////////////////////////////////////////////////////////
        // check if this catalogue object is in second catalogue
        int tev_select = -1;
        
        // loop over second catalogue
        for( unsigned int j = 0; j < tev->getNStar(); j++ )
        {
            double r = VAstronometry::vlaDsep( s->getStarRA2000( i ) * TMath::Pi() / 180., s->getStarDec2000( i ) * TMath::Pi() / 180., tev->getStarRA2000( j ) * TMath::Pi() / 180., tev->getStarDec2000( j ) * TMath::Pi() / 180. ) * 180. / TMath::Pi();
            if( r < r_min )
            {
                tev_select = ( int )j;
                cout << "tevCAT: " << tev->getStarName( j ) << endl;
                if( fEventList_tevcat )
                {
                    fEventList_tevcat->Enter( i );
                }
                break;
            }
        }
        if( tev_select > 0 )
        {
            cout << "\t OBJECT " << s->getStarName( i ) << " is nearby " << tev->getStarName( tev_select ) << endl;
            cout << endl;
        }
        
        ////////////////////////////////////////////////////////////////
        // now check if object is in FOV
        
        // loop over all runs
        r_tot = 0.;
        r_tot_V5 = 0.;
        iVERITAS_targets.clear();
        for( unsigned int j = 0; j < fRunRA.size(); j++ )
        {
            // calculate distance of catalogue object to camera center
            r_centre = VAstronometry::vlaDsep( s->getStarRA2000( i ) * TMath::Pi() / 180., s->getStarDec2000( i ) * TMath::Pi() / 180.,
                                ( fRunRA[j] + fRunoffsetRA[j] ) * TMath::Pi() / 180., ( fRunDec[j] + fRunoffsetDec[j] ) * TMath::Pi() / 180. ) * 180. / TMath::Pi();
            // do dqm
            if( !doDQM( j, iMinDuration ) )
            {
                continue;
            }
            
            // check if the catalogue object is close enough
            if( r_centre < iR )
            {
                // calculate distance of catalogue object to VERITAS object
                r_VA_object = VAstronometry::vlaDsep( s->getStarRA2000( i ) * TMath::Pi() / 180., s->getStarDec2000( i ) * TMath::Pi() / 180.,
                                       fRunRA[j] * TMath::Pi() / 180., fRunDec[j] * TMath::Pi() / 180. ) * 180. / TMath::Pi();
                // total time on object (all array configurations)
                r_tot += fRunDuration[j];
                // total time on object (new array configuration only)
                if( fRun[j] > 46642 )
                {
                    r_tot_V5 += fRunDuration[j];
                }
                // print some information about this run
                cout << "\tRUN " << fRun[j];
                if( fRunSourceID[j].size() > 0 )
                {
                    cout << "\t" << fRunSourceID[j];
                }
                cout << "\t MJD " << fRunStartMJD[j];
                cout << "\t CONFIGMASK " << fRunConfigMask[j];
                cout << "\t (ra,dec)=(" << fRunRA[j] << "," << fRunDec[j] << ")";
                cout << "\t DIST " << r_centre << " deg";
                cout << "\t DIST (VTS pointing) " << r_VA_object << " deg";
                cout << endl;
                // mean calculation
                r_centre_mean += r_centre;
                r_VA_object_mean += r_VA_object;
                r_meanElevation += fRunTelElevation[j];
                r_N++;
                iVERITAS_targets.insert( fRunSourceID[j] );
            }
        }
        // print summary only if some data was taken
        if( r_tot > 0. )
        {
            if( fEventList_inFOV )
            {
                fEventList_inFOV->Enter( i );
            }
            cout << "\t\t total time: " << r_tot / 3600. << " [h] (V5: " << r_tot_V5 / 3600. << ")" << endl;
            if( r_N > 0. )
            {
                cout << "\t\t mean distance to camera centre: " << r_centre_mean / r_N << " [deg]" << endl;
                cout << "\t\t mean distance to VERITAS object: " << r_VA_object_mean / r_N << " [deg]" << endl;
                cout << "\t\t mean elevation: " << r_meanElevation / r_N << " [deg]" << endl;
                if( fEventList_inFOV && r_VA_object_mean / r_N > 0.1 )
                {
                    fEventList_inFOV_noTevcat->Enter( i );
                }
            }
        }
    }
    if( fEventList_tevcat )
    {
        fEventList_tevcat->Write();
    }
    if( fEventList_inFOV )
    {
        fEventList_inFOV->Write();
    }
    if( fEventList_inFOV_noTevcat )
    {
        fEventList_inFOV_noTevcat->Write();
    }
    if( fEventListfile )
    {
        fEventListfile->Close();
    }
}

void VExposure::printShortRunList()
{
    for( unsigned int i = 0; i < fRunSourceID.size(); i++ )
    {
        cout << fRunSourceID[i] << "\t" << fRunDuration[i] << "\t" << fRun[i] << endl;
    }
}


void VExposure::printListOfRuns( double il, double ib, double iR, double iMinDuration, string iDQMfileList, string ofile, unsigned int iVerbose )
{
    double r_dist = 0.;
    double r_tot = 0.;
    
    // read list with runs which passed dqm
    vector< int > i_vDQMList;
    if( iDQMfileList.size() > 0 )
    {
        ifstream is;
        is.open( iDQMfileList.c_str(), ifstream::in );
        if( !is )
        {
            cout << "DQM file list not found: " << iDQMfileList << endl;
            return;
        }
        string is_line;
        while( getline( is, is_line ) )
        {
            if( atoi( is_line.c_str() ) > 0 )
            {
                i_vDQMList.push_back( atoi( is_line.c_str() ) );
            }
        }
        is.close();
    }
    ofstream os;
    if( ofile.size() > 0 )
    {
        os.open( ofile.c_str() );
        if( !os )
        {
            cout << "Error writing list of runs to " << ofile << endl;
            return;
        }
    }
    
    if( iVerbose == 1 )
    {
        cout << "Number of runs:  " << fRunGalLong1958.size() << endl;
    }
    else if( iVerbose == 2 && fRunGalLong1958.size() > 0 )
    {
        cout << "Number of runs: " << fRunGalLong1958.size() << endl;
    }
    for( unsigned int i = 0; i < fRunGalLong1958.size(); i++ )
    {
        r_dist = VAstronometry::vlaDsep( il * TMath::Pi() / 180., ib * TMath::Pi() / 180., fRunGalLong1958[i] * TMath::Pi() / 180., fRunGalLat1958[i] * TMath::Pi() / 180. ) * 180. / TMath::Pi();
        
        if( r_dist < iR && fRunDuration[i] > iMinDuration )
        {
            // check if this file passed the dqm
            bool bPassed = false;
            for( unsigned int q = 0; q < i_vDQMList.size(); q++ )
            {
                if( fRun[i] == i_vDQMList[q] )
                {
                    bPassed = true;
                    break;
                }
            }
            if( i_vDQMList.size() == 0 )
            {
                bPassed = true;
            }
            if( !bPassed )
            {
                if( iVerbose > 0 && i_vDQMList.size() > 0 )
                {
                    cout << fRun[i] << " FAILED DQM" << endl;
                }
                continue;
            }
            
            r_tot += fRunDuration[i];
            cout << "\t\tRUN " << fRun[i];
            if( fRunSourceID[i].size() > 0 )
            {
                cout << "\t" << fRunSourceID[i];
            }
            cout << "\t MJD " << fRunStartMJD[i];
            cout << "\t CONFIGMASK " << fRunConfigMask[i];
            cout << "\t DIST " << r_dist << " deg";
            cout << endl;
            if( os )
            {
                // do not write T1T2 runs to disk
                if( fRunConfigMask[i] != 3 )
                {
                    os << fRun[i] << endl;
                }
            }
        }
    }
    if( iVerbose > 0 )
    {
        cout << "\t\t total time: " << r_tot / 3600. << " [h]" << endl;
    }
    if( os )
    {
        os.close();
    }
    
}


/*
   copied from ROOT forum #6650

   (input in degrees)
*/
void VExposure::aitoff2xy( Double_t l, Double_t b, Double_t& Al, Double_t& Ab )
{
    Double_t x, y;
    
    Double_t alpha2 = ( l / 2 ) * TMath::DegToRad();
    Double_t delta  = b * TMath::DegToRad();
    Double_t r2     = TMath::Sqrt( 2. );
    Double_t f      = 2 * r2 / TMath::Pi();
    Double_t cdec   = TMath::Cos( delta );
    Double_t denom  = TMath::Sqrt( 1. + cdec * TMath::Cos( alpha2 ) );
    x      = cdec * TMath::Sin( alpha2 ) * 2.*r2 / denom;
    y      = TMath::Sin( delta ) * r2 / denom;
    x     *= TMath::RadToDeg() / f;
    y     *= TMath::RadToDeg() / f;
    x *= -1.;                                     // for a skymap swap left<->right
    Al = x;
    Ab = y;
}


void VExposure::plotMarker( double l, double b, double r, string iText, int iMarkerStyle, int iMarkerColor, int iMarkerSize, double iTextAngle )
{
    double ib_range = 180;
    double il_range = 360;
    
    if( fPlottingCanvas.size() == 0 )
    {
        return;
    }
    
    for( unsigned int i = 0; i < fPlottingCanvas.size(); i++ )
    {
        if( !fPlottingCanvas[i] )
        {
            continue;
        }
        
        fPlottingCanvas[i]->cd();
        bool bAitoff = false;
        double ab = b;
        double al = -1.*l;
        string i_cname = fPlottingCanvas[i]->GetName();
        if( i_cname.find( "aitoff" ) < i_cname.size() )
        {
            bAitoff = true;
            aitoff2xy( l, b, al, ab );
        }
        // extended?
        if( r > 0. )
        {
            TEllipse* m = new TEllipse( al, ab, r, r );
            m->SetFillStyle( 0 );
            m->Draw();
        }
        else
        {
            TMarker* m = new TMarker( al, ab, iMarkerStyle );
            m->SetMarkerColor( iMarkerColor );
            m->SetMarkerSize( iMarkerSize );
            m->Draw();
        }
        
        if( iText.size() > 0 )
        {
            if( bAitoff )
            {
                al = al + il_range * 0.01 * ( cos( iTextAngle * TMath::DegToRad() ) + sin( iTextAngle * TMath::DegToRad() ) );
                ab = ab + ib_range * 0.01 * ( -1.*sin( iTextAngle * TMath::DegToRad() ) + cos( iTextAngle * TMath::DegToRad() ) );
            }
            else
            {
                al = -1.*l + il_range * 0.01 * ( cos( iTextAngle * TMath::DegToRad() ) + sin( iTextAngle * TMath::DegToRad() ) );
                ab = b + ib_range * 0.01 * ( -1.*sin( iTextAngle * TMath::DegToRad() ) + cos( iTextAngle * TMath::DegToRad() ) );
            }
            TText* t = new TText( al, ab, iText.c_str() );
            t->SetTextColor( iMarkerColor );
            t->SetTextSize( t->GetTextSize() * 0.9 );
            t->SetTextAngle( iTextAngle );
            t->Draw();
            
        }
        fPlottingCanvas[i]->Update();
    }
}


void VExposure::addCatalogue( string iC, int iMarker, int iColor, double iAngle )
{
    fCatalogue.push_back( iC );
    fCatalogueMarkerColor.push_back( iColor );
    fCatalogueMarkerStyle.push_back( iMarker );
    fCatalogueTextAngle.push_back( iAngle );
}


bool VExposure::removeCataloge( unsigned int iB )
{
    if( iB < fCatalogue.size() )
    {
        fCatalogue.erase( fCatalogue.begin() + iB );
        fCatalogueMarkerColor.erase( fCatalogueMarkerColor.begin() + iB );
        fCatalogueMarkerStyle.erase( fCatalogueMarkerStyle.begin() + iB );
        fCatalogueTextAngle.erase( fCatalogueTextAngle.begin() + iB );
        return true;
    }
    return false;
}


void VExposure::listCatalogues()
{
    cout << "list of catalogues: " << endl;
    for( unsigned int i = 0; i < fCatalogue.size(); i++ )
    {
        cout << i << "\t" << fCatalogue[i] << "\t" << fCatalogueMarkerColor[i] << "\t";
        cout << fCatalogueMarkerStyle[i] << "\t" << fCatalogueTextAngle[i] << endl;
    }
}


void VExposure::plotTimeDifferencesbetweenRuns()
{
    fTimeDifferencesBetweenRuns = new TH1D( "hTimeDifferenceBetweenRuns", "", 1000, 0., 100. );
    fTimeDifferencesBetweenRuns->SetXTitle( "time difference [min]" );
    fTimeDifferencesBetweenRuns->SetYTitle( "# of runs" );
    fTimeDifferencesBetweenRuns->SetStats( 0 );
    
    if( fRunStartMJD.size() < 2 )
    {
        return;
    }
    
    double iDiff = 0.;
    for( unsigned int i = 0; i < fRunStartMJD.size() - 1; i++ )
    {
        iDiff  = fRunStartMJD[i + 1] - fRunStopMJD[i];
        iDiff *= 24. * 60.;
        
        fTimeDifferencesBetweenRuns->Fill( iDiff );
    }
    
    TCanvas* cTimeDiff = new TCanvas( "cTimeDiff", "time difference between runs", 10, 10, 400, 400 );
    cTimeDiff->SetGridx( 0 );
    cTimeDiff->SetGridy( 0 );
    if( fTimeDifferencesBetweenRuns->GetEntries() > 0. )
    {
        cTimeDiff->SetLogy( 1 );
    }
    
    fTimeDifferencesBetweenRuns->SetAxisRange( 0., 10. );
    fTimeDifferencesBetweenRuns->Draw();
}

void VExposure::printListOfRuns()
{

    int k = 0;
    int itemp = 0;
    double Total_Time = 0;
    vector< int > tcut;
    
    //////////////////////////////////
    // loop over all runs
    for( unsigned int j = 0; j < fRunRA.size(); j++ )
    {
        // total time on object (new array configuration only)
        if( fRunTelElevation[j] >= fTelMinElevation && fRunDuration[j] >= fMinDuration )
        {
            cout << fixed ;
            cout << "RUN " << fRun[j] << " (" << fRunObsMode[j] << ") ";
            cout << "\t" <<  fRunSourceID[j];
            cout << setprecision( 4 ) << "\tMJD " << fRunStartMJD[j];
            cout << "\tDate: " << fRunDate[j];
            cout << setprecision( 1 ) << "\tWobble:(N,E) (" << fWobbleN[j] << "," << fWobbleE[j] << ")";
            cout << "\tCONFIGMASK " << fRunConfigMask[j];
            cout << "\t(ra,dec)=(" << fRunRA[j] << "," << fRunDec[j] << ")";
            cout << "\tDuration[min] " << fRunDuration[j] / 60.;
            cout << "\t(El,Az) " << fRunTelElevation[j] << " " << fRunTelAzimuth[j];
            if( bPrintVerbose )
            {
                cout << "\t Category: " << fDataCat[j];
                cout << "\t Status: " << fStatus[j];
                cout << "\t Reason: " << fStatReason[j];
                cout << "\t Telescope Mask: " << fTelCutMask[j];
                cout << "\t Usable: " << fUsable[j];
                cout << "\t TimeCut: " << fTimeCutMask[j];
                cout << "\t Light Level: " << fLightLevel[j];
                if( j < fVPMcon.size() )
                {
                    cout << "\t VPM: " << fVPMcon[j];
                }
                cout << "\t Comment: " << fComment[j];
                cout << "\t Author: " << fAuthor[j];
            }
            else
            {
                cout << "\t " << fStatus[j];
            }
            
            if( fSelectLaser == 1 )
            {
                cout << "\tLaser: " << fRunLaserList[j][0];
                cout << " " << fRunLaserList[j][1];
                cout << " " << fRunLaserList[j][2];
                cout << " " << fRunLaserList[j][3];
            }
            cout << endl;
            
            k++;
            Total_Time += fRunDuration[j];
            
            fRunDownload.push_back( fRun[j] );
            fRunDownloadDate.push_back( fRunDate[j] );
            
            if( fSelectLaser == 1 )
            {
                for( unsigned int i = 0; i < 4 ; i++ )
                {
                    itemp = 0;
                    for( unsigned int m = 0 ; m < fLaserDownload.size() ; m++ )
                        if( fLaserDownload[m] == fRunLaserList[j][i] )
                        {
                            itemp = 1;
                        }
                    if( itemp == 0 )
                    {
                        fLaserDownload.push_back( fRunLaserList[j][i] );
                        fLaserDownloadDate.push_back( getLaserDate( fRunLaserList[j][i] ) );
                    }
                }
            }
            
        }
        
    }
    
    cout << endl;
    cout << "Selected " << k << " runs after cuts. A total of " << Total_Time / 60. / 60. << " [hrs]." << endl;
    cout << endl;
    
    if( fSelectLaser == 1 )
    {
    
        cout << endl;
        cout << "Printing Separate Set of Laser Run List (for easy copy n paste):" << endl;
        
        for( unsigned int i = 0 ; i < fLaserDownload.size() ; i++ )
        {
            cout << "\t" << fLaserDownloadDate[i] << "\t" << fLaserDownload[i] << endl;
        }
        cout << endl;
        
    }
    
    if( bPrintTimeMask )
    {
    
        int k = 0;
        double wordtemp = 0.;
        
        cout << endl;
        cout << "Time Masks For Eventdisplay:" << endl;
        cout << endl;
        
        for( unsigned int j = 0; j < fRunRA.size(); j++ )
        {
        
            if( fRunTelElevation[j] >= fTelMinElevation && fRunDuration[j] >= fMinDuration )
            {
            
                if( fTimeCutMask[j] != " " )
                {
                    stringstream stream1( fTimeCutMask[j] );
                    string word1;
                    
                    while( getline( stream1, word1, ',' ) )
                    {
                        stringstream stream2( word1 );
                        string word2;
                        
                        cout << "* " ;
                        cout << " " << fRun[j];
                        k = 0;
                        while( getline( stream2, word2, '/' ) )
                        {
                            if( k == 0 )
                            {
                                cout << " " << word2 ;
                                wordtemp = atof( word2.c_str() );
                            }
                            else if( k == 1 )
                            {
                                cout << " " << atof( word2.c_str() ) - wordtemp;
                            }
                            else
                            {
                                cout << "WARNING: In Time Mask Calculation." << endl;
                            }
                            k++;
                        }
                        cout << " 0" << endl;
                    }
                }
                
            }
            
        }
        
    }
    
}

void VExposure::outputAnasumRunlist( string fAnasumFile )
{

    FILE* anasumFile;
    
    string fatmo;
    string fepoch;
    string fconfig;
    string freplace;
    
    if( fAnasumFile != "" )
    {
        anasumFile = fopen( fAnasumFile.c_str(), "w" );
    }
    else
    {
        return;
    }
    
    cout << "ANSUM output (warning: OLD, USE WITH CAUTION): " << endl;
    
    // loop over all runs
    for( unsigned int j = 0; j < fRunRA.size(); j++ )
    {
        fepoch = "";
        fconfig = "";
        freplace = "REPLACE_";
        
        // check for atmospheric conditions (dates are set to full moon periods)
        fatmo = "ATM";
        // did not looked into such old data files yet (set to undefined atmosphere)
        if( fRunDate[j] < 20071026 )
        {
            fatmo += "XX_";
        }
        else if( fRunDate[j] > 20080420 && fRunDate[j] < 20081113 )
        {
            fatmo += "22_";
        }
        else if( fRunDate[j] > 20090509 && fRunDate[j] < 20091102 )
        {
            fatmo += "22_";
        }
        else if( fRunDate[j] > 20100428 && fRunDate[j] < 20101023 )
        {
            fatmo += "22_";
        }
        else if( fRunDate[j] > 20110418 && fRunDate[j] < 20111110 )
        {
            fatmo += "22_";
        }
        else if( fRunDate[j] > 20120506 && fRunDate[j] < 20121029 )
        {
            fatmo += "22_";
        }
        // TODO: not sure which of the two dates (still missing some data from Tucson) ???
        // else if( fRunDate[j] > 20130425 && fRunDate[j] < 20131019 ) fatmo += "22_";
        else if( fRunDate[j] > 20130425 && fRunDate[j] < 20131117 )
        {
            fatmo += "22_";
        }
        // did not looked into such new data files yet (set to undefined atmosphere)
        else if( fRunDate[j] > 20140401 )
        {
            fatmo += "XX_";
        }
        else
        {
            fatmo += "21_";
        }
        
        
        if( fRun[j] < 46642 )
        {
            fepoch += "V4_";
        }
        else if( fRun[j] > 46641 && fRun[j] < 63373 )
        {
            fepoch += "V5_";
        }
        else if( fRun[j] > 63372 )
        {
            fepoch += "V6_";
        }
        
        if( fRunConfigMask[j] == 15 )
        {
            fconfig = "1234";
        }
        else if( fRunConfigMask[j] == 7 )
        {
            fconfig = "123";
        }
        else if( fRunConfigMask[j] == 11 )
        {
            fconfig = "124";
        }
        else if( fRunConfigMask[j] == 13 )
        {
            fconfig = "134";
        }
        else if( fRunConfigMask[j] == 14 )
        {
            fconfig = "234";
        }
        else
        {
            fconfig = "2Tel";
        }
        
        freplace += fatmo;
        freplace += fepoch;
        freplace += fconfig;
        
        // total time on object (new array configuration only)
        if( fRunTelElevation[j] >= fTelMinElevation && fRunDuration[j] >= fMinDuration )
        {
        
            if( fStatus[j] == "good_run" )
            {
                cout << "* ";
            }
            cout << fRun[j] << " ";
            cout << fRun[j] << " ";
            cout << " 0 ";
            cout << freplace.c_str() << endl;
            
            if( fStatus[j] == "good_run" )
            {
                fprintf( anasumFile, "* %d %d %s\n", fRun[j], fRun[j], freplace.c_str() );
            }
            else
            {
                fprintf( anasumFile, "%d %d %s\n", fRun[j], fRun[j], freplace.c_str() );
            }
            
        }
        
    }
    
    
    if( fatmo == "ATMXX_" )
    {
        cout << endl;
        cout << "Warning: You are using either old or new data runs for which no transition dates have been set in VExposure::outputAnasumRunlist()" << endl;
        cout << "         Please check radiosonde measurements to find the transition date (they should optiminally be choosen during a full moon period) ..." << endl;
    }
    
    cout << endl;
    cout << "ANASUM output written to: " << fAnasumFile.c_str() << endl;
    cout << endl;
    
    fclose( anasumFile );
    
}

void VExposure::checkRunList()
{

    char filename[800];
    char filepath[800];
    
    char* ENVIR_VAR;
    
    ENVIR_VAR = getenv( "VERITAS_DATA_DIR" );
    
    cout << "Checking data runs: " << endl;
    
    for( unsigned int i = 0 ; i < fRunDownload.size(); i++ )
    {
    
        sprintf( filename, "%s/data/d%d/%d.cvbf", ENVIR_VAR, fRunDownloadDate[i], fRunDownload[i] );
        sprintf( filepath, "%s/data/d%d/", ENVIR_VAR, fRunDownloadDate[i] );
        
        ifstream datafile( filename );
        if( datafile.is_open() )
        {
            cout << filename << "\tOn Disk" << endl;
            checkMD5sum( fRunDownloadDate[i], fRunDownload[i] );
        }
        else
        {
            cout << filename << "\tNeeds Downloading" << endl;
        }
        
    }
    
    cout << endl;
    
    
    cout << "Checking flasher runs: " << endl;
    
    for( unsigned int i = 0 ; i < fLaserDownload.size(); i++ )
    {
    
        sprintf( filename, "%s/data/d%d/%d.cvbf", ENVIR_VAR, fLaserDownloadDate[i], fLaserDownload[i] );
        sprintf( filepath, "%s/data/d%d/", ENVIR_VAR, fLaserDownloadDate[i] );
        
        ifstream datafile( filename );
        if( datafile.is_open() )
        {
            cout << filename << "\tOn Disk" << endl;
            checkMD5sum( fLaserDownloadDate[i], fLaserDownload[i] );
        }
        else
        {
            cout << filename << "\tNeeds Downloading" << endl;
        }
        
    }
    
    cout << endl;
    
}

void VExposure::downloadRunList()
{

    char filename[800];
    char filepath[800];
    char dl_string[800];
    
    char* ENVIR_VAR;
    
    char mkdir_string[2000];
    char permission_string[2000];
    
    ENVIR_VAR = getenv( "VERITAS_DATA_DIR" );
    
    cout << "Download " << fRunDownload.size() << " data runs: " << endl;
    
    for( unsigned int i = 0 ; i < fRunDownload.size(); i++ )
    {
    
        sprintf( filename, "%s/data/d%d/%d.cvbf", ENVIR_VAR, fRunDownloadDate[i], fRunDownload[i] );
        sprintf( filepath, "%s/data/d%d/", ENVIR_VAR, fRunDownloadDate[i] );
        
        ifstream datafile( filename );
        if( datafile.is_open() )
        {
            cout << filename << " File Exists!" << endl;
            checkMD5sum( fRunDownloadDate[i], fRunDownload[i] );
        }
        else
        {
            ifstream paths( filepath );
            if( !paths.is_open() )
            {
                sprintf( mkdir_string, "mkdir %s", filepath );
                cout << "COMMAND: " << mkdir_string << endl;
                if( system( mkdir_string ) != 0 )
		{
		    cout << "error mkdir " << filepath << endl;
		    exit( EXIT_FAILURE );
                }
                sprintf( permission_string, "chmod g+sw %s", filepath );
                cout << "COMMAND: " << permission_string << endl;
		if( system( permission_string ) != 0 )
		{
		    cout << "error setting permissions" << endl;
		    exit( EXIT_FAILURE );
                }
            }
            
            if( system( "which bbftp" ) != 0 )
            {
                cout << "ERROR: \"which bbftp\" shows no match. Install bbftp and add to you $PATH." << endl;
                cout << "exiting ...." << endl;
                exit( EXIT_FAILURE );
            }
            // Then Download
            sprintf( dl_string, "bbftp -V -S -p 12 -u bbftp -e \"get /veritas/data/d%d/%d.cvbf %s/data/d%d/%d.cvbf\" %s", fRunDownloadDate[i], fRunDownload[i],
                     ENVIR_VAR, fRunDownloadDate[i], fRunDownload[i],
                     getRawDataServer().c_str() );
            sprintf( permission_string , "chmod g+w %s", filename );
            cout << dl_string << endl;
	    if( system( dl_string ) != 0 )
	    {
	        cout << "error downloading" << endl;
		cout << dl_string << endl;
		exit( EXIT_FAILURE );
            }
            cout << permission_string << endl;
	    if( system( permission_string ) != 0 )
	    {
	        cout << "error setting permissions" << endl;
		cout << permission_string << endl;
		exit( EXIT_FAILURE );
            }
            
            checkMD5sum( fRunDownloadDate[i], fRunDownload[i] );
        }
        
        
    }
    
    cout << endl;
    
    if( fSelectLaser == 1 )
    {
    
        cout << "Download " << fLaserDownload.size() << " flasher runs: " << endl;
        
        for( unsigned int i = 0 ; i < fLaserDownload.size(); i++ )
        {
        
            sprintf( filename, "%s/data/d%d/%d.cvbf", ENVIR_VAR, fLaserDownloadDate[i], fLaserDownload[i] );
            sprintf( filepath, "%s/data/d%d/", ENVIR_VAR, fLaserDownloadDate[i] );
            
            ifstream datafile( filename );
            if( datafile.is_open() )
            {
                cout << filename << " File Exists!" << endl;
                checkMD5sum( fLaserDownloadDate[i], fLaserDownload[i] );
            }
            else
            {
                ifstream paths( filepath );
                if( !paths.is_open() )
                {
                    sprintf( mkdir_string, "mkdir %s", filepath );
                    cout << "COMMAND: " << mkdir_string << endl;
                    if( system( mkdir_string ) != 0 )
		    { 
		        cout << "error making directory " << endl;
                    }
                    sprintf( permission_string, "chmod g+sw %s", filepath );
                    cout << "COMMAND: " << permission_string << endl;
                    if( system( permission_string ) != 0 )
		    { 
		        cout << "error setting permissions" << endl;
                    }
                }
                if( system( "which bbftp" ) != 0 )
                {
                    cout << "ERROR: \"which bbftp\" shows no match. Install bbftp and add to you $PATH." << endl;
                    cout << "exiting ...." << endl;
                    exit( EXIT_FAILURE );
                }
                // Then Download
                sprintf( dl_string, "bbftp -V -S -p 12 -u bbftp -e \"get /veritas/data/d%d/%d.cvbf %s/data/d%d/%d.cvbf\" %s", fLaserDownloadDate[i], fLaserDownload[i],
                         ENVIR_VAR, fLaserDownloadDate[i], fLaserDownload[i],
                         getRawDataServer().c_str() );
                sprintf( permission_string , "chmod g+w %s", filename );
                if( fLaserDownload[i] > 0 )
                {
                    cout << dl_string << endl;
                    if( system( dl_string ) != 0 )
		    {
                        cout << "error downloading" << endl;
                    }
                    cout << permission_string << endl;
                    if( system( permission_string ) != 0 )
		    {
                        cout << "error setting permission string" << endl;
                    }
                    checkMD5sum( fLaserDownloadDate[i], fLaserDownload[i] );
                }
                else
                {
                    cout << "Skipping Run " << fLaserDownload[i] << endl;
                }
                
            }
            
        }
        
    }
    
}

int VExposure::checkMD5sum( int date, int run, bool force_download )
{

    if( !fDoCheckSums )
    {
        return 0;
    }
    
    TString desy_sum = calcMD5sum( date, run );
    if( desy_sum == "" )
    {
        cout << "VExposure::checkMD5sum Warning: Unable to calculate checksum for run " << run << ", date " << date << endl;
        fRunsNoChecksum.push_back( run );
        return -1;
    }
    TString archive_sum = getArchiveMD5sum( date, run, force_download );
    if( archive_sum == "" )
    {
        cout << "VExposure::checkMD5sum Warning: Unable to download checksum for run " << run << ", date " << date << endl;
        fRunsNoChecksum.push_back( run );
        return -1;
    }
    
    if( archive_sum != desy_sum )
    {
        cout << "VExposure::checkMD5sum Warning: Checksum error for run " << run << ", date " << date << endl;
        cout << "Checksum from archive: " << archive_sum << ", calculated checksum: " << desy_sum << endl;
        fRunsBadChecksum.push_back( run );
        return 1;
    }
    cout << "VExposure::checkMD5sum: Run " << run << " survived the testsum check." << endl;
    cout << "Checksum from archive: " << archive_sum << ", calculated checksum: " << desy_sum << endl;
    fRunsGoodChecksum.push_back( run );
    return 0;
}

void VExposure::printChecksumSummary()
{
    if( !fDoCheckSums )
    {
        return;
    }
    
    unsigned int ntot = fRunsNoChecksum.size() + fRunsGoodChecksum.size() + fRunsBadChecksum.size() ;
    cout << "Checksum summary: " << ntot << " run(s); " << fRunsBadChecksum.size() << " bad, " << fRunsGoodChecksum.size() << " good, " << fRunsNoChecksum.size() << " could not be checked." << endl;
    
    if( fRunsGoodChecksum.size() == ntot && ntot > 0 )
    {
        cout << "Congratulations, all runs survived the checksum test. " << endl;
        
    }
    else
    {
        if( fRunsGoodChecksum.size() > 0 )
        {
            cout << "Got good checksums for " << fRunsGoodChecksum.size() << " run(s)."  << endl << "Good runs:\t" ;
            for( unsigned int i = 0; i <  fRunsGoodChecksum.size() ; i++ )
            {
                cout << fRunsGoodChecksum.at( i ) << "\t" ;
            }
            cout << endl << endl;
        }
        
        if( fRunsNoChecksum.size() > 0 )
        {
            cout << "Warning, checksums could not be compared for " << fRunsNoChecksum.size();
            cout << " run(s), see above for details. (Checksums are not available from the archive for old or very new runs.)"  << endl << "Not checked:\t" ;
            for( unsigned int i = 0; i <  fRunsNoChecksum.size() ; i++ )
            {
                cout << fRunsNoChecksum.at( i ) << "\t" ;
            }
            cout << endl << endl;
        }
        if( fRunsBadChecksum.size() > 0 )
        {
            cout << "Warning, wrong checksums for " << fRunsBadChecksum.size() << " run(s), see above for details. "  << endl << "Bad run(s): " ;
            for( unsigned int i = 0; i <  fRunsBadChecksum.size() ; i++ )
            {
                cout << fRunsBadChecksum.at( i ) << "\t" ;
            }
            cout << endl << endl;
        }
    }
    
}

TString VExposure::getArchiveMD5sum( int date, int run, bool force_download )
{

    bool attempt_download = true;
    
    //check if md5sum utility is available
    if( system( "which bbftp" ) != 0 )
    {
        cout << "VExposure::getArchiveMD5sum Error: \"which bbftp\" shows no match; will attempt to DL checksum file" << endl;
        attempt_download = false;
        if( force_download )
        {
            return "";
        }
    }
    
    char* ENVIR_VAR;
    ENVIR_VAR = getenv( "VERITAS_DATA_DIR" );
    
    TString basecamp_sumfilename = TString::Format( "%s/data/d%d/sumd%d", ENVIR_VAR, date, date );
    TString UTAH_sumfilename = TString::Format( "%s/data/d%d/CHPC_sumd%d", ENVIR_VAR, date, date );
    TString UCLA_sumfilename = TString::Format( "%s/data/d%d/UCLA_sumd%d", ENVIR_VAR, date, date );
    TString UCLA2_sumfilename = TString::Format( "%s/data/d%d/UCLA_sum", ENVIR_VAR, date );
    TString UCLA_new_sumfilename = TString::Format( "%s/data/d%d/UCLA_new_sum", ENVIR_VAR, date );
    //see elog http://veritash.sao.arizona.edu:8081/VERITAS-Operations/11544 . Some runs were reprocessed on disk at ucla to correct 2 rows of swapped signal cables.
    
    TString dlcommand = TString::Format( "bbftp -V -S -p 12 -u bbftp -e \"mget /veritas/data/d%d/*sum* %s/data/d%d/\" %s", date, ENVIR_VAR, date, getRawDataServer().c_str() );
    TString chmodcommand = TString::Format( "chmod g+w %s/data/d%d/*sum*", ENVIR_VAR, date );
    
    bool use_new_ucla_sumfile = ( run >= 69474 && run <= 69641 );
    
    TString checksum = "";
    //attempt to find checksum on disk
    if( !force_download )
    {
        if( use_new_ucla_sumfile )
        {
            checksum = readMD5sumFromFile( UCLA_new_sumfilename, run, !attempt_download );
        }
        else
        {
            checksum = readMD5sumFromFile( basecamp_sumfilename, run, !attempt_download );
            if( checksum == "" )
            {
                checksum = readMD5sumFromFile( UCLA_sumfilename, run, !attempt_download );
            }
            if( checksum == "" )
            {
                checksum = readMD5sumFromFile( UCLA2_sumfilename, run, !attempt_download );
            }
            if( checksum == "" )
            {
                checksum = readMD5sumFromFile( UTAH_sumfilename, run, !attempt_download );
            }
        }
        if( checksum != "" )
        {
            return checksum;
        }
        
    }
    
    if( attempt_download )
    {
        cout << dlcommand << endl;
        if( system( dlcommand.Data() ) != 0 )
        {
            cout << "VExposure::getArchiveMD5sum Error: Problem executing command : " << dlcommand << endl;
        }
	if( system( chmodcommand.Data() ) != 0 )
	{
	    cout << "error running " << dlcommand << endl;
        }
        if( use_new_ucla_sumfile )
        {
            checksum = readMD5sumFromFile( UCLA_new_sumfilename, run );
        }
        else
        {
            checksum = readMD5sumFromFile( basecamp_sumfilename, run );
            if( checksum == "" )
            {
                checksum = readMD5sumFromFile( UCLA_sumfilename, run );
            }
            if( checksum == "" )
            {
                checksum = readMD5sumFromFile( UCLA2_sumfilename, run );
            }
            if( checksum == "" )
            {
                checksum = readMD5sumFromFile( UTAH_sumfilename, run );
            }
        }
        if( checksum != "" )
        {
            return checksum;
        }
    }
    
    return "";
    
}

TString VExposure::readMD5sumFromFile( TString filename, int run, bool warn )
{

    ifstream intemp( filename.Data() );
    string linestring;
    TString iChecksum;
    TString iFile;
    
    TString search = TString::Format( "%d.cvbf", run );
    TString checksum = "";
    
    if( ! intemp.is_open() )
    {
        if( warn )
        {
            cout << "VExposure::readMD5sumFromFile Error: Unable to open file: " << filename << endl;
        }
        return "";
    }
    else
    {
        while( getline( intemp, linestring ) )
        {
            stringstream line( linestring );
            line >> iChecksum >> iFile;
            if( iFile.EndsWith( search.Data() ) )
            {
                checksum = iChecksum;
            }
        }
        intemp.close();
    }
    if( checksum == "" && warn )
    {
        cout <<   "VExposure::readMD5sumFromFile Error: Unable to find checksum for " << search << " in file: " << filename << endl;
    }
    return checksum;
    
}

TString VExposure::calcMD5sum( int date, int run )
{

    //check if md5sum utility is available
    if( system( "which md5sum" ) != 0 )
    {
        cout << "VExposure::calcMD5sum Error: \"which md5sum\" shows no match; will not calculate checksum." << endl;
        return "";
    }
    
    TString datafilename;
    TString resultfilename;
    TString tempfilename;
    TString filepath;
    char* ENVIR_VAR;
    
    ENVIR_VAR = getenv( "VERITAS_DATA_DIR" );
    
    datafilename.Form( "%s/data/d%d/%d.cvbf", ENVIR_VAR, date, run );
    resultfilename.Form( "%s/data/d%d/LOCAL_sumd%d", ENVIR_VAR, date, date );
    tempfilename.Form( "%s/data/d%d/TEMP_sum_%d_%d", ENVIR_VAR, date, run, ( int )time( 0 ) );
    filepath.Form( "%s/data/d%d/", ENVIR_VAR, date );
    
    //check if input file exists
    ifstream itemp( datafilename.Data() );
    if( ! itemp.good() )
    {
        itemp.close();
        cout << "VExposure::calcMD5sum Error: File " << datafilename.Data() << " not available. "  << endl;
        return "";
    }
    itemp.close();
    
    //now do the checksum.
    TString command = TString::Format( "md5sum %s | tr \"\\n\" \" \" > %s ; date >> %s", datafilename.Data(), tempfilename.Data(), tempfilename.Data() );
    cout << command << endl;
    if( system( command.Data() ) != 0 )
    {
        cout << "VExposure::calcMD5sum Error: Unable to execute command: " << command << endl;
        return "";
    }
    
    //read in checksum
    TString checksum = readMD5sumFromFile( tempfilename, run );
    
    if( checksum != "" )
    {
        command.Form( "cat %s >> %s", tempfilename.Data(), resultfilename.Data() );
        if( system( command.Data() ) == 0 )
        {
            command.Form( "rm %s", tempfilename.Data() );
	    if( system( command.Data() ) != 0 )
	    {
	        cout << "error calculating md5sum" << endl;
            }
            command.Form( "chmod g+w %s", resultfilename.Data() );
            if( system( command.Data() ) != 0 )
	    {
	        cout << "error setting permissions" << endl;
            }
        }
    }
    return checksum;
}

void VExposure::getLaserList()
{

    cout << "Reading Laser Runs for selected runs ....... " ;
    
    for( unsigned int i = 0; i < fRunRA.size(); i++ )
    {
        fRunLaserList.push_back( getLaserRun( fRun[i], 4 ) );
    }
    
    cout << "complete." << endl;
    cout << endl;
    
    return;
    
}

// lucie: This function is an old copy (no read of DQM mask) of the function with the same name in VDBRunInfo
// when used should be replaced by VDBRunInfo one.
// This function should disappear
vector< unsigned int > VExposure::getLaserRun( unsigned int iRunNumber, unsigned int iNTel )
{

    string iTempS;
    iTempS = getDBServer() + "/VERITAS";
    char c_query[1000];
    
    sprintf( c_query, "SELECT info.run_id, grp_cmt.excluded_telescopes FROM tblRun_Info AS info, tblRun_Group AS grp, tblRun_GroupComment AS grp_cmt, (SELECT group_id FROM tblRun_Group WHERE run_id=%d) AS run_grp WHERE grp_cmt.group_id = run_grp.group_id AND grp_cmt.group_type='laser' AND grp_cmt.group_id=grp.group_id AND grp.run_id=info.run_id AND (info.run_type='flasher' OR info.run_type='laser')", iRunNumber );
    
    //std::cout<<"VExposure::getLaserRun "<<std::endl;
    VDB_Connection my_connection( iTempS.c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        return fLaserRunID;
    }
    if( !my_connection.make_query( c_query ) )
    {
        return fLaserRunID;
    }
    TSQLResult* db_res = my_connection.Get_QueryResult();
    
    
    vector< unsigned int > iLaserList;
    vector< unsigned int > iLaserExclude;
    if( db_res->GetRowCount() > 0 )
    {
        while( TSQLRow* db_row = db_res->Next() )
        {
            if( !db_row )
            {
                cout << "VDBRunInfo: failed reading a row from DB for run " << iRunNumber << endl;
                //              fDBStatus = false;
                return fLaserRunID;
            }
            iLaserList.push_back( atoi( db_row->GetField( 0 ) ) );
            iLaserExclude.push_back( atoi( db_row->GetField( 1 ) ) );
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
        for( unsigned int i = 0; i < iLaserList.size(); i++ )
        {
            bitset< 8 > ibit( iLaserExclude[i] );
            if( !ibit.test( t ) )
            {
                fLaserRunID[t] = iLaserList[i];
            }
        }
    }
    
    return fLaserRunID;
}


void VExposure::setMakeRunList( bool iSet )
{

    if( iSet )
    {
        fMakeRunList = true;
    }
    
    return;
    
}

void VExposure::readRunListFromFile( string runlist )
{

    string is_line;
    ifstream inputfile;
    inputfile.open( runlist.c_str() );
    
    if( !inputfile )
    {
        cout << "ERROR: Input runlist not found: " << runlist << endl;
        exit( -1 );
    }
    
    while( getline( inputfile, is_line ) )
    {
        fRunDownloadList.push_back( atoi( is_line.c_str() ) );
    }
    
    inputfile.close();
    
}

void VExposure::readRunCommentsFromDB()
{

    stringstream iTempS;
    iTempS << getDBServer() << "/VOFFLINE";
    
    //std::cout<<"VExposure::readRunCommentsFromDB "<<std::endl;
    VDB_Connection my_connection( iTempS.str().c_str(), "readonly", "" ) ;
    if( !my_connection.Get_Connection_Status() )
    {
        return;
    }
    char c_query[1000];
    
    
    for( unsigned int i = 0; i < fRun.size(); i++ )
    {
    
        sprintf( c_query, "SELECT run_id , data_category   , status   , status_reason , tel_cut_mask , usable_duration , time_cut_mask , light_level , vpm_config_mask , authors  , comment from tblRun_Analysis_Comments where run_id=%d", fRun[i] );
        
        if( !my_connection.make_query( c_query ) )
        {
            return;
        }
        TSQLResult* db_res = my_connection.Get_QueryResult();
        
        TSQLRow* db_row = db_res->Next();
        if( !db_row )
        {
            if( fRun[i] > 46905 )
            {
                cout << "VDBExposure:Comments: failed reading a row from DB for run " << fRun[i] << endl;
            }
            fDataCat.push_back( " " );
            fStatus.push_back( " " );
            fStatReason.push_back( " " );
            fTelCutMask.push_back( " " );
            fUsable.push_back( " " );
            fTimeCutMask.push_back( " " );
            fLightLevel.push_back( " " );
            fVPMcon.push_back( " " );
            fAuthor.push_back( " " );
            fComment.push_back( " " );
        }
        else
        {
            if( db_row->GetField( 1 ) )
            {
                fDataCat.push_back( db_row->GetField( 1 ) );
            }
            else
            {
                fDataCat.push_back( " " );
            }
            if( db_row->GetField( 2 ) )
            {
                fStatus.push_back( db_row->GetField( 2 ) );
            }
            else
            {
                fStatus.push_back( " " );
            }
            if( db_row->GetField( 3 ) )
            {
                fStatReason.push_back( db_row->GetField( 3 ) );
            }
            else
            {
                fStatReason.push_back( " " );
            }
            if( db_row->GetField( 4 ) )
            {
                fTelCutMask.push_back( db_row->GetField( 4 ) );
            }
            else
            {
                fTelCutMask.push_back( " " );
            }
            if( db_row->GetField( 5 ) )
            {
                fUsable.push_back( db_row->GetField( 5 ) );
            }
            else
            {
                fUsable.push_back( " " );
            }
            if( db_row->GetField( 6 ) )
            {
                fTimeCutMask.push_back( db_row->GetField( 6 ) );
            }
            else
            {
                fTimeCutMask.push_back( " " );
            }
            if( db_row->GetField( 7 ) )
            {
                fLightLevel.push_back( db_row->GetField( 7 ) );
            }
            else
            {
                fLightLevel.push_back( " " );
            }
            if( db_row->GetField( 8 ) )
            {
                fVPMcon.push_back( db_row->GetField( 8 ) );
            }
            else
            {
                fVPMcon.push_back( " " );
            }
            if( db_row->GetField( 9 ) )
            {
                fAuthor.push_back( db_row->GetField( 9 ) );
            }
            else
            {
                fAuthor.push_back( " " );
            }
            if( db_row->GetField( 10 ) )
            {
                fComment.push_back( db_row->GetField( 10 ) );
            }
            else
            {
                fComment.push_back( " " );
            }
        }
    }
    
    return;
}
void VExposure::setPrintTimeMask( int iPrintTimeMask )
{
    if( iPrintTimeMask )
    {
        bPrintTimeMask = true;
    }
}
void VExposure::setPrintVerbose( int iPrintVerbose )
{
    if( iPrintVerbose )
    {
        bPrintVerbose = true;
    }
}

void VExposure::readLaserRunListFromFile( string runlist )
{

    string is_line;
    ifstream inputfile;
    inputfile.open( runlist.c_str() );
    
    if( !inputfile )
    {
        cout << "ERROR: Input LASER runlist not found: " << runlist << endl;
        exit( -1 );
    }
    
    while( getline( inputfile, is_line ) )
    {
        cout << atoi( is_line.c_str() ) << " " << getLaserDate( atoi( is_line.c_str() ) ) << endl;
        fLaserDownload.push_back( atoi( is_line.c_str() ) );
        fLaserDownloadDate.push_back( getLaserDate( atoi( is_line.c_str() ) ) );
    }
    
    inputfile.close();
    
}

void VExposure::readLaserRunDateListFromFile( string runlist )
{

    ifstream inputfile;
    inputfile.open( runlist.c_str() );
    
    if( !inputfile )
    {
        cout << "ERROR: Input LASER runlist not found: " << runlist << endl;
        exit( -1 );
    }
    
    int run = -1;
    int date_run = -1;
    
    while( 1 )
    {
        inputfile >> run >> date_run ;
        if( !inputfile.good() )
        {
            break;
        }
        fLaserDownload.push_back( run );
        fLaserDownloadDate.push_back( date_run );
    }
    
    inputfile.close();
    
}


void VExposure::setRunNumber( unsigned int number )
{

    fRunDownloadList.push_back( number );
    
}

void VExposure::setLaserNumber( unsigned int number )
{

    fLaserDownload.push_back( number );
    fLaserDownloadDate.push_back( getLaserDate( number ) );
    
}

void VExposure::setObservingMode( bool bObsMode )
{

    if( bObsMode )
    {
        fObservingMode = "Special";
    }
    else
    {
        fObservingMode = "Normal";
    }
    
}


void VExposure::setCanvasSize( double x, double y )
{
    fCanvasSize_x = x;
    fCanvasSize_y = y;
}
