/* \class VArrayPointing

   \brief data class for general pointing values valid for all telescopes

   e.g. reference pointing, target, e.g.



*/

#include "VArrayPointing.h"

VArrayPointing::VArrayPointing( bool bInitTree )
{
    setObservatory();
    reset();
    
    if( bInitTree )
    {
        initializePointingTree();
    }
}


void VArrayPointing::initializePointingTree()
{
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "pointingData" );
    sprintf( htitle, "pointing (Array)" );
    fPointingTree = new TTree( hname, htitle );
    fPointingTree->Branch( "MJD", &fMJD, "MJD/i" );
    fPointingTree->Branch( "Time", &fTime, "Time/D" );
    fPointingTree->Branch( "TargetAzimuth", &fTargetAzimuth, "TargetAzimuth/D" );
    fPointingTree->Branch( "TargetElevation", &fTargetElevation, "TargetElevation/D" );
    fPointingTree->Branch( "TargetRAJ2000", &fTargetRAJ2000, "TargetRAJ2000/D" );
    fPointingTree->Branch( "TargetDecJ2000", &fTargetDecJ2000, "TargetDecJ2000/D" );
    fPointingTree->Branch( "TargetRA", &fTargetRA, "TargetRA/D" );
    fPointingTree->Branch( "TargetDec", &fTargetDec, "TargetDec/D" );
    // TelAzimuth and TelAzimthCalculate should be the same (same for Elevation)
    fPointingTree->Branch( "TelAzimuth", &fTelAzimuth, "TelAzimuth/D" );
    fPointingTree->Branch( "TelElevation", &fTelElevation, "TelElevation/D" );
    fPointingTree->Branch( "TelAzimuthCalculated", &fTelAzimuthCalculated, "TelAzimuthCalculated/F" );
    fPointingTree->Branch( "TelElevationCalculated", &fTelElevationCalculated, "TelElevationCalculated/F" );
    fPointingTree->Branch( "TelRAJ2000", &fTelRAJ2000, "TelRAJ2000/D" );
    fPointingTree->Branch( "TelDecJ2000", &fTelDecJ2000, "TelDecJ2000/D" );
    fPointingTree->Branch( "TelRA", &fTelRA, "TelRA/D" );
    fPointingTree->Branch( "TelDec", &fTelDec, "TelDec/D" );
    
    
    sprintf( hname, "pointingDataReduced" );
    sprintf( htitle, "pointing information for the entire array, at 1 interpolated measurement per second" );
    fPntReduced = new TTree( hname, htitle );
    fPntReduced->Branch( "MJD",          &fMJD,          "MJD/i" );
    fPntReduced->Branch( "Time",         &fTime,         "Time/D" );
    fPntReduced->Branch( "TelAzimuth",   &fTelAzimuth,   "TelAzimuth/D" );
    fPntReduced->Branch( "TelElevation", &fTelElevation, "TelElevation/D" );
    fPntReduced->Branch( "TelRAJ2000",   &fTelRAJ2000,   "TelRAJ2000/D" );
    fPntReduced->Branch( "TelDecJ2000",  &fTelDecJ2000,  "TelDecJ2000/D" );
}


void VArrayPointing::fillPointingTree( bool iIsMC )
{
    if( fPointingTree && !iIsMC )
    {
        fPointingTree->Fill();
    }
}

void VArrayPointing::fillPntReduced()
{
    if( fPntReduced )
    {
        fPntReduced->Fill();
    }
}

void VArrayPointing::terminate( bool iDebug_IO, bool bIsMC )
{
    // full pointing tree only written to disk for data, not for MC
    if( fPointingTree && !bIsMC )
    {
        int i_nbytes = fPointingTree->Write();
        if( iDebug_IO )
        {
            cout <<  "WRITEDEBUG: pointing tree (nbytes " << i_nbytes << ")" << endl;
        }
    }
    if( fPntReduced )
    {
        int i_nbytes = fPntReduced->Write() ;
        if( iDebug_IO )
        {
            cout <<  "WRITEDEBUG: pointing tree (reduced) (nbytes " << i_nbytes << ")" << endl;
        }
    }
}

