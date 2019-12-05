/*! \class VPointing
    \brief get telescope pointing direction

*/

#include "VPointing.h"

VPointing::VPointing( unsigned int iTelID )
{
    fTelID = iTelID;
    
    fPointingTree = 0;
    fPointingDB = 0;
    fUseDB = false;
    fPointingType = 0;
    
    fTelAzimuthDB = 0.;
    fTelElevationDB = 0.;
    fNEventsWithNoDBPointing = 0;
    fPointingDB = 0;
    fEventStatus = 0;
    
    fPointingErrorX = 0.;
    fPointingErrorY = 0.;
    fMeanPointingErrorN = 0;
    fMeanPointingErrorX = 0.;
    fMeanPointingErrorY = 0.;
    fMeanPointingDistance = 0.;
    
    reset();
    setObservatory();
    
    initializePointingTree();
}


/*

    calculate azimuth and elevation of telescope pointing direction
    calculate azimuth and elevation of target

 */
void VPointing::setTelPointing( int MJD, double time, bool iUseDB, bool iFillPointingTree )
{
    fUseDB = iUseDB;
    
    // update pointing
    updatePointing( MJD, time );
    
    // telescope elevation/azimuth calculated from VERITAS DB entries
    if( fUseDB )
    {
        updatePointingfromDB( fMJD, fTime );
    }
    // calulation from source should always be successful
    else
    {
        fEventStatus = 1;
    }
    
    // now set the global elevation/azimuth to be used for the analysis
    if( fUseDB && fPointingType > 1 )
    {
        fTelElevation = fTelElevationDB;
        fTelAzimuth   = fTelAzimuthDB;
    }
    else
    {
        fTelElevation = fTelElevationCalculated;
        fTelAzimuth   = fTelAzimuthCalculated;
    }
    
    // fill pointing tree
    if( iFillPointingTree )
    {
        fillPointingTree();
    }
    
}

void VPointing::getPointingFromDB( int irun, string iTCorrection, string iVPMDirectory, bool iVPMDB, bool iUncalibratedVPM )
{
    fPointingType = 2;
    if( iVPMDB == true )
    {
        fPointingType = 6;    // read VPM data from VERITAS DB
    }
    else if( iUncalibratedVPM == true )
    {
        fPointingType = 5;    // read uncalibrated VPM data from VERITAS DB
    }
    else if( iVPMDirectory.size() > 0 )
    {
        fPointingType = 4;    // read VPM data from a text file
    }
    else if( iTCorrection.size() > 0 )
    {
        fPointingType = 3;    // read raw positioner data from VERITAS DB and apply tracking corrections
    }
    else
    {
        fPointingType = 2;    // read T-Point corrected positioner data from VERITAS DB
    }
    
#ifdef RUNWITHDB
    fPointingDB = new VPointingDB( fTelID, irun );
    fPointingDB->setObservatory( fObsLongitude * TMath::RadToDeg(), fObsLatitude * TMath::RadToDeg() );      // work in [deg]
    fPointingDB->initialize( iTCorrection, iVPMDirectory, iVPMDB, iUncalibratedVPM );
    if( !fPointingDB->isGood() )
    {
        cout << endl;
        cout << "FATAL ERROR: cannot connect to VERITAS database" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
#else
    fPointingDB = 0;
#endif
}

/*

   update pointing from VERITAS database

*/
bool VPointing::updatePointingfromDB( int MJD, double iTime )
{
    // do something if we read stuff from the db
    if( fPointingDB )
    {
        fPointingDB->updatePointing( MJD, iTime );
        
        // telescope pointings
        fTelAzimuthDB   = fPointingDB->getTelAzimuthDB();
        fTelElevationDB = fPointingDB->getTelElevationDB();
        fEventStatus    = fPointingDB->getEventStatus();
        
        if( fEventStatus != 3 )
        {
            // calculate pointing error in camera coordinates
            double iPx = 0.;
            double iPy = 0.;
            int j = 0;
            VAstronometry::vlaDs2tp( fTelAzimuthCalculated / TMath::RadToDeg(), fTelElevationCalculated / TMath::RadToDeg(),
                      fTelAzimuthDB / TMath::RadToDeg(), fTelElevationDB / TMath::RadToDeg(),
                      &iPx, &iPy, &j );
            if( j == 0 )
            {
                // azimuth from North to East
                fPointingErrorX = iPx * TMath::RadToDeg();
                // evndisp camera is upside down
                fPointingErrorY = iPy * TMath::RadToDeg();
            }
            // star not on tangent plane
            else
            {
                fPointingErrorX = 0.;
                fPointingErrorY = 0.;
                fEventStatus = 4;
            }
        }
        else
        {
            fPointingErrorX = 0.;
            fPointingErrorY = 0.;
            fNEventsWithNoDBPointing++;
        }
        fMeanPointingErrorN++;
        fMeanPointingErrorX += fPointingErrorX;
        fMeanPointingErrorY += fPointingErrorY;
        fMeanPointingDistance += sqrt( fPointingErrorX * fPointingErrorX + fPointingErrorY * fPointingErrorY );
    }
    if( fEventStatus != 3 && fEventStatus != 4 )
    {
        return true;
    }
    
    return false;
}

/*

   mainly check that the different information about pointing matches

   write results to disk

*/
void VPointing::terminate( bool i_isMC )
{
    // don't do anything for MC
    if( i_isMC )
    {
        return;
    }
    
    cout << "\t mean pointing mismatch between eventdisplay and DB for telescope " << getTelID() + 1 << ":  (x,y,r) [deg] ";
    if( fMeanPointingErrorN > 0 )
    {
        cout << fMeanPointingErrorX / ( double )fMeanPointingErrorN << ", " << fMeanPointingErrorY / ( double )fMeanPointingErrorN;
        cout << ", " << fMeanPointingDistance / ( double )fMeanPointingErrorN;
    }
    else
    {
        cout << "0., 0., 0.";
    }
    if( fNEventsWithNoDBPointing > 0 )
    {
        cout << ", number of events with no pointing information from the database: " << fNEventsWithNoDBPointing;
    }
    cout << endl;
    
    if( fMeanPointingErrorN > 0 && ( fMeanPointingDistance / ( double )fMeanPointingErrorN > 0.1 ) )
    {
        cout << "WARNING: LARGE MISMATCH BETWEEN EVENTDISPLAY AND DB POINTING DATA FOR TELESCOPE " << getTelID() + 1 << endl;
    }
    
    //  write results to disk
    if( fPointingDB )
    {
        TTree* t = fPointingDB->getTreePointingDB();
        if( t )
        {
            t->Write();
        }
        fPointingDB->terminate();
    }
    if( fPointingTree )
    {
        fPointingTree->Write();
    }
}


void VPointing::initializePointingTree()
{
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "pointing_%d", getTelID() + 1 );
    sprintf( htitle, "pointing (Telescope %d)", getTelID() + 1 );
    fPointingTree = new TTree( hname, htitle );
    fPointingTree->Branch( "MJD", &fMJD, "MJD/i" );
    fPointingTree->Branch( "Time", &fTime, "Time/D" );
    // elevation / azimuth for target
    // (calculated from target coordinates in run parameters)
    fPointingTree->Branch( "TargetAzimuth", &fTargetAzimuth, "TargetAzimuth/D" );
    fPointingTree->Branch( "TargetElevation", &fTargetElevation, "TargetElevation/D" );
    // ra / dec for target for J2000
    // (fixed target coordinates from  run parameters)
    fPointingTree->Branch( "TargetRAJ2000", &fTargetRAJ2000, "TargetRAJ2000/D" );
    fPointingTree->Branch( "TargetDecJ2000", &fTargetDecJ2000, "TargetDecJ2000/D" );
    // ra / dec for target for current epoch
    // (calcualted from target coordinates in run parameters)
    fPointingTree->Branch( "TargetRA", &fTargetRA, "TargetRA/D" );
    fPointingTree->Branch( "TargetDec", &fTargetDec, "TargetDec/D" );
    // elevation / azimuth for telescopes
    // (default: from DB (pointing monitior); otherwise calculate from target coordinates or vbf values)
    // (so default means that TelAzimuth = TelAzimuthDB and TelElevation = TelElevationDB)
    fPointingTree->Branch( "TelAzimuth", &fTelAzimuth, "TelAzimuth/D" );
    fPointingTree->Branch( "TelElevation", &fTelElevation, "TelElevation/D" );
    // flag indicated pointing type: see beginning of VPointing::getPointingFromDB() for the meaning of these flags
    fPointingTree->Branch( "PointingType", &fPointingType, "fPointingType/i" );
    // eleation / azimuth calculated from target coordinates and wobble offset
    fPointingTree->Branch( "TelAzimuthCalculated", &fTelAzimuthCalculated, "TelAzimuthCalculated/F" );
    fPointingTree->Branch( "TelElevationCalculated", &fTelElevationCalculated, "TelElevationCalculated/F" );
    // elevation / azimuth for telescopes from DB (pointing monitor)
    fPointingTree->Branch( "TelAzimuthDB", &fTelAzimuthDB, "TelAzimuthDB/F" );
    fPointingTree->Branch( "TelElevationDB", &fTelElevationDB, "TelElevationDB/F" );
    // difference between expected az/el from target coordinates and DB el/az
    fPointingTree->Branch( "PointingErrorX", &fPointingErrorX, "PointingErrorX/F" );
    fPointingTree->Branch( "PointingErrorY", &fPointingErrorY, "PointingErrorY/F" );
    // status of pointing results -> never used and need to be checked
    fPointingTree->Branch( "EventStatus", &fEventStatus, "EventStatus/i" );
    // telescope pointing directions in ra/dec
    fPointingTree->Branch( "TelRAJ2000", &fTelRAJ2000, "TelRAJ2000/D" );
    fPointingTree->Branch( "TelDecJ2000", &fTelDecJ2000, "TelDecJ2000/D" );
    fPointingTree->Branch( "TelRA", &fTelRA, "TelRA/D" );
    fPointingTree->Branch( "TelDec", &fTelDec, "TelDec/D" );
}

void VPointing::fillPointingTree()
{
    if( fPointingTree )
    {
        fPointingTree->Fill();
    }
}

/*

   set an artifical pointing error (e.g. from command line)

*/
void VPointing::setPointingError( double iX, double iY )
{
    fPointingType = 1;
    fPointingErrorX = iX;
    fPointingErrorY = iY;
    
    fMeanPointingErrorN = 1;
    fMeanPointingErrorX = fPointingErrorX;
    fMeanPointingErrorY = fPointingErrorY;
    fMeanPointingDistance = sqrt( fPointingErrorX * fPointingErrorX + fPointingErrorY * fPointingErrorY );
}
