//
// VDeadPixelOrganizer
// top-level class for taking in pixel information,
// converting it to a smaller format,
// and writing it to a TTree, which is passed up to
// the anasum-ctools event converter
//
// Nathan Kelley-Hoskins
// Aug 2014

#include "VDeadPixelOrganizer.h"
using namespace std ;

void VDeadPixelOrganizer::setupMap()
{
    //cout << coutprefix << "VDeadPixelOrganizer::setupMaps()" << endl;
    
    if( ntel <= 0 || npix <= 0 )
    {
        cout << coutprefix << "VDeadPixelOrganizer::setupMaps(): Error, must have ntel and npix > 0, exiting..." << endl;
        exit( 1 ) ;
    }
    
    int telid = 0 ;
    for( int i_tel = 0 ; i_tel < ntel ; i_tel++ )
    {
        telid = i_tel + 1 ;
        //cout << coutprefix << "adding telescope " << telid << endl;
        telescopeMap[telid] = VNTelescope() ;
        telescopeMap[telid].initialize( telid, npix, this ) ;
        //cout << telescopeMap[telid] << endl;
        
    }
    //cout << "NKH VDeadPixelOrganizer::VDeadPixelOrganizer(" << ntel<< "," << npix << "): telescopeMap.size()=" << telescopeMap.size() << endl;
    
}

void VNTelescope::initialize( int telid, int npix, VDeadPixelOrganizer* parentVDPO )
{
    setParentVDeadPixelOrganizerPointer( parentVDPO ) ;
    setTelid( telid ) ;
    setNpix( npix ) ;
    setupMap() ;
    
}

void VNTelescope::setupMap()
{
    //cout << coutprefix << "VNTelescope::setupMaps()" << endl;
    int pixid ;
    for( int i_pix = 0 ; i_pix < npix ; i_pix++ )
    {
        pixid = i_pix + 1 ;
        pixelMap[pixid] = VNPixel() ;
        pixelMap[pixid].initialize( pixid, this ) ;
    }
}

void VNPixel::initialize( int pixid, VNTelescope* parentVNT )
{
    setParentVNTelescopePointer( parentVNT ) ;
    setPixid( pixid ) ;
    setupPixelCoord( fTelescope->fOrganizer->fDetectorGeo ) ;
    setupGainMap() ;
}

void VNTelescope::setParentVDeadPixelOrganizerPointer( VDeadPixelOrganizer* vdpo )
{
    fOrganizer = vdpo ;
}

void VNPixel::setParentVNTelescopePointer( VNTelescope* ptr )
{
    fTelescope = ptr ;
}

void VNPixel::setPixid( int pid )
{
    pixid_evndisp  = pid - 1 ;
    pixid_dbtables = pid   ;
}


void VNGain::setParentVNPixelPointer( VNPixel* ptr )
{
    fPixel = ptr ;
}

void VNPixel::setupPixelCoord( VDetectorGeometry* fDetectorGeo )
{
    int telid = fTelescope->getTelid() - 1 ;
    pixelCoord.xcam = fTelescope->fOrganizer->fDetectorGeo->getXUnrotated( telid )[pixid_evndisp] ;
    pixelCoord.ycam = fTelescope->fOrganizer->fDetectorGeo->getYUnrotated( telid )[pixid_evndisp] ;
    //cout << coutprefix << "======== VNPIxel::setupPixelCoord() id=" << pixid_dbtables << " (" << pixelCoord.xcam << "," << pixelCoord.ycam << ")" << endl ;
}

// after
void VNPixel::setupGainMap()
{
    gainMap["low" ] = VNGain() ;
    gainMap["low" ].initialize( "low", this ) ;
    gainMap["high"] = VNGain() ;
    gainMap["high"].initialize( "high", this ) ;
}

void VNGain::initialize( string gainType, VNPixel* parentVNP )
{

    // handle the parent pointer and the pixel's camera coordinate
    setParentVNPixelPointer( parentVNP ) ;
    
    // figure out which gain we're using
    fGainName = gainType ;
    if( strcmp( "low", gainType.c_str() ) == 0 )
    {
        fIsLow    = true ;
        fGainChar = 'L'  ;
        //cout << coutprefix << "  detected as 'low'" << endl;
    }
    else if( strcmp( "high", gainType.c_str() ) == 0 )
    {
        fIsLow    = false ;
        fGainChar = 'H'   ;
        //cout << coutprefix << "  detected as 'high'" << endl;
    }
    else
    {
        cerr << coutprefix << "VNGain::VNGain() Error, gainType '" << gainType << "' unrecognized, exiting..." << endl;
        exit( 1 ) ;
    }
    
    // set up the first entry of stateHistory
    // starts at beginning of run, with state 0 (ok)
    // leave endtime empty until the next time updateState is called
    stateHistory.push_back( VNStateDuration() ) ;
    stateHistory.back().begtime = fPixel->fTelescope->fOrganizer->fRunStart ;
    stateHistory.back().state   = 0 ;
    
    // add the beginning
    
    
}

// I think c++ creates an empty copy of this class first
// when creating the maps, so in those, copies the pointers
// fPixel and fTelescope dont exist (fPixel=0x0, etc).
// So, if these exist, we know these are the real objects
// (and not the copies) so we can actually print info
// Otherwise, if we don't check fPixel and fTelescope,
// we get segfaults when trying to use these pointers.
// Really only matters in the constructer, I think.
bool VNGain::isRealObject()
{
    if( fPixel )
    {
        if( fPixel->fTelescope )
        {
            if( fPixel->fTelescope->fOrganizer )
            {
                return true ;
            }
        }
    }
    return false ;
    
}

bool VNPixel::isRealObject()
{
    if( fTelescope )
    {
        if( fTelescope->fOrganizer )
        {
            return true ;
        }
    }
    return false ;
}

bool VNTelescope::isRealObject()
{
    if( fOrganizer )
    {
        return true ;
    }
    return false ;
}

// tel: 1-4
// pix: 1-499
// mjd: day of the event
// time: seconds since beginning of mjd day
// lowGain = true means its the low gain's state
void VDeadPixelOrganizer::UpdatePixelState( int tel, int pix, bool lowGain, int mjd, double time, PixelStateInt state )
{

    char buff[200] ;
    sprintf( buff, "VDeadPixelOrganizer::UpdatePixelState( tel=%d, pix=%3d, gain=%s, mjd=%d, time=%f, state=%d)",
             tel, pix, lowGain ? "low" : "high" , mjd, time, state ) ;
    /*
    if ( tel <= tellimit && pix <= pixUlimit && pix >= pixLlimit ) {
    	cout << coutprefix << buff << endl;
    }
    */
    
    if( tel <= ntel )
    {
        if( pix <= telescopeMap[tel].getNpix() )
        {
            telescopeMap[tel].pixelMap[pix].gainMap[ lowGain ? "low" : "high" ].updateState( mjd, time, state ) ;
        }
        else
        {
            cerr << coutprefix << "VDeadPixelOrganizer::UpdatePixelState(): error, pixel '" << pix << "' doesn't exist in telescope '" << tel << "', exiting..." << endl;
        }
    }
    else
    {
        cerr << coutprefix << "VDeadPixelOrganizer::UpdatePixelState(): error, telescope '" << tel << "' doesn't exist, exiting..." << endl;
    }
}

void VDeadPixelOrganizer::setDetectorGeo( VDetectorGeometry* detgeo )
{
    fDetectorGeo = detgeo ;
}
VDetectorGeometry* VDeadPixelOrganizer::getDetectorGeo()
{
    return fDetectorGeo ;
}

VDeadPixelOrganizer::VDeadPixelOrganizer( int ntel, int npix, VDetectorGeometry* detGeo, int startMJD, double startSec, int endMJD, double endSec, string treename, int runNumber )
{
    setTreeName( treename ) ;
    setNtel( ntel ) ;
    setNpix( npix ) ;
    setDetectorGeo( detGeo ) ;
    fRunStart = VNTime( startMJD, startSec ) ;
    fRunEnd   = VNTime( endMJD,   endSec ) ;
    fPreviousEvent.fEventNumber = 0 ;
    fPreviousEvent.fEventTime   = fRunStart ;
    setRunNumber( runNumber ) ;
    setupMap() ;
    
    // writes copy of output tree to csv, in addition to the TTree
    fWriteCSV = true ;
}

void VDeadPixelOrganizer::printSummary()
{

    cout << endl;
    cout << coutprefix << " VDeadPixelOrganizer::printSummary():" << endl ;
    
    VNTelescope vnt ;
    VNPixel     vnp ;
    VNGain      vng ;
    MapOfVNTelescopes::const_iterator telIter ;
    for( telIter  = telescopeMap.begin() ; telIter != telescopeMap.end() ; ++telIter )
    {
    
        vnt = telIter->second ;
        MapOfVNPixels::const_iterator pixIter ;
        for( pixIter = vnt.pixelMap.begin() ; pixIter != vnt.pixelMap.end() ; ++pixIter )
        {
        
            vnp = pixIter->second ;
            /*
            MapOfVNGains::const_iterator gainIter ;
            for ( gainIter = vnp.gainMap.begin() ; gainIter != vnp.gainMap.end() ; ++gainIter ) {
            
                vng = gainIter->second ;
            
            	char buff[200] ;
            	sprintf( buff, "T%d P%3d G%c #%d", vnt.getTelid(), vnp.pixid_dbtables, vng.fGainChar, (int)vng.stateHistory.size()  ) ;
            	cout << coutprefix << "sum: " << buff << endl;
            
            
            }
            */
            vng = vnp.gainMap["high"] ;
            char buff[200] ;
            sprintf( buff, "T%d P%3d GL #%d", vnt.getTelid(), vnp.pixid_dbtables, ( int )vng.stateHistory.size() ) ;
            
            if( vnt.getTelid() <= tellimit && vnp.pixid_dbtables <= pixUlimit && vnp.pixid_dbtables >= pixLlimit )
            {
                cout << coutprefix << "sum: " << buff << endl;
            }
        }
    }
    cout << coutprefix << " VDeadPixelOrganizer::printSummary()" << endl ;
    cout << endl;
    
}

void VNGain::updateState( int mjd, double time, PixelStateInt state )
{
    //int telid = fPixel->fTelescope->getTelid() ;
    //int pixid = fPixel->pixid_dbtables ;
    
    char buff[200] ;
    sprintf( buff, "VNGain::updateState( tel=%d, pix=%d, gain=%s, mjd=%d, time=%f, state=%d )",
             fPixel->fTelescope->getTelid(),
             fPixel->pixid_dbtables,
             fGainName.c_str(), mjd, time, state ) ;
    /*
    if ( fPixel->fTelescope->getTelid() <= tellimit && fPixel->pixid_dbtables <= pixUlimit && fPixel->pixid_dbtables >= pixLlimit ) {
    	cout << coutprefix << buff << endl;
    	cout << coutprefix << "        updateState( seeing if new state(" << state << ") matches old state(" << stateHistory.back().state << ") needs to be added..." << endl;
    }
    */
    
    // if the last state doesn't match the current state,
    // we should end the last state's duration, and add a new one
    if( stateHistory.back().state != state )
    {
        /*
        sprintf( buff, "T%d P%3d G%c (%f,%f)", telid, pixid, fGainChar, fPixel->pixelCoord.xcam, fPixel->pixelCoord.ycam ) ;
        cout << coutprefix << "        updateState(" << buff << "): NEW STATE: state(" << state << ", mjd=" << mjd<< ", time=" << time << ") doesn't match previous state(" << stateHistory.back().state << "), will add new history row" << endl;
        */
        
        // if its different, then add it to the history
        int oldsize = stateHistory.size() ;
        
        // figure out the previous event's time was
        // (*not* the previous row's time)
        VNTime lastCheckedTime = fPixel->fTelescope->fOrganizer->fPreviousEvent.fEventTime ;
        
        // find halfway between the current events time and the previous event's time
        VNTime dividingTime    = lastCheckedTime.halfwayTime( VNTime( mjd, time ) ) ;
        
        // close off the last history row
        stateHistory.back().endtime = dividingTime ;
        
        // add a new history row for the new state
        stateHistory.push_back( VNStateDuration() ) ;
        stateHistory.back().begtime = dividingTime ;
        stateHistory.back().state   = state ;
        
        // output size change
        int newsize = stateHistory.size() ;
        if( fPixel->fTelescope->getTelid() <= tellimit && fPixel->pixid_dbtables <= pixUlimit && fPixel->pixid_dbtables >= pixLlimit )
        {
            cout << coutprefix << "        updateState(): added new row, histry had " << oldsize << " elements, now has " << newsize << " elements..." << endl;
        }
    }
    
    /*
    else
    {
    
    	if ( fPixel->fTelescope->getTelid() <= tellimit && fPixel->pixid_dbtables <= pixUlimit && fPixel->pixid_dbtables >= pixLlimit ) {
    		cout << coutprefix << "        updateState(): state(" << state << ") matches previous state, not doing anything..." << endl;
    	}
    }
    */
    
}

VNTime::VNTime( int mjd_days, double time_seconds )
{
    mjd  = mjd_days     ;
    time = time_seconds ;
}

// find point in time halfway between object's time and argument's time
VNTime VNTime::halfwayTime( const VNTime vnt )
{
    VNTime out = VNTime() ;
    if( mjd == vnt.mjd )
    {
        out.mjd = mjd ;
        out.time = ( int )( time + vnt.time ) / 2.0 ;
        return out ;
    }
    else
    {
        double objtime =     mjd + ( time / 86400.0 ) ;
        double argtime = vnt.mjd + ( vnt.time / 86400.0 ) ;
        
        double halftime = ( objtime + argtime ) / 2.0 ;
        
        // strip off the mantissa to get the mjd
        out.mjd  = ( int ) halftime ;
        
        // get the seconds from the mantissa
        out.time = ( halftime - out.mjd ) * 86400.0 ;
        return out ;
        
    }
    return out ;
}


void VDeadPixelOrganizer::updatePreviousEventInfo( int eventNumber, int eventMJD, double eventTime )
{
    fPreviousEvent.fEventNumber    = eventNumber ;
    fPreviousEvent.fEventTime.mjd  = eventMJD    ;
    fPreviousEvent.fEventTime.time = eventTime   ;
}

void VDeadPixelOrganizer::topOffRows()
{

    //cout << coutprefix << "closeAllOpen() fRunEnd:" << fRunEnd << endl;
    //char buff[200] ;
    VNTelescope* vnt ;
    VNPixel*      vnp ;
    VNGain*       vng ;
    MapOfVNTelescopes::iterator  telIter  ;
    MapOfVNPixels    ::iterator  pixIter  ;
    MapOfVNGains     ::iterator gainIter ;
    for( telIter  = telescopeMap.begin() ; telIter != telescopeMap.end() ; ++telIter )
    {
    
        vnt = &telIter->second ;
        for( pixIter = vnt->pixelMap.begin() ; pixIter != vnt->pixelMap.end() ; ++pixIter )
        {
        
            vnp = &pixIter->second ;
            for( gainIter = vnp->gainMap.begin() ; gainIter != vnp->gainMap.end() ; ++gainIter )
            {
            
                vng = &gainIter->second ;
                
                //sprintf( buff, "T%d P%3d G%c", vnt->getTelid(), vnp->getPixid(), vng->fGainChar ) ;
                // close off this gain's last duration in its state history
                //char buff2[200] = "" ;
                //sprintf( buff2, "closeAllOpen(%s) before: ", buff ) ;
                //cout << coutprefix << buff2 << vng->stateHistory.back().endtime << endl;
                vng->stateHistory.back().endtime = fRunEnd ;
                //cout << coutprefix << "                         after : " << vng->stateHistory.back().endtime << endl;
            }
        }
    }
}

//
void VDeadPixelOrganizer::organize()
{
    // start up ttree and storage vars
    fRegTree = new TTree( fTreeName.c_str(), "Tree of pixel states for this run" ) ;
    int telid, pixid, begtimeMJD, endtimeMJD, runnum;
    double begtimeSec, endtimeSec, pixXcam, pixYcam ;
    char gainChar ;
    PixelStateInt state ;
    
    // start
    FILE* pFile = 0 ;
    if( fWriteCSV )
    {
        pFile = fopen( "pixReg.csv", "w" ) ;
        fprintf( pFile, "runnumber, telid, pixid, gain,       xcam,       ycam, begMJD,         begSec, endMJD,         endSec, state\n" ) ;
    }
    
    // branches
    fRegTree->Branch( "runnumber" , &runnum     , "runnumber/I" ) ;
    fRegTree->Branch( "telid"     , &telid      , "telid/I" ) ;
    fRegTree->Branch( "pixid"     , &pixid      , "pixid/I" ) ;
    fRegTree->Branch( "gain"      , &gainChar   , "gainChar/B" ) ;
    fRegTree->Branch( "pixXcam"   , &pixXcam    , "pixXcam/D" ) ;
    fRegTree->Branch( "pixYcam"   , &pixYcam    , "pixYcam/D" ) ;
    fRegTree->Branch( "begTimeMJD", &begtimeMJD , "begTimeMJD/I" ) ;
    fRegTree->Branch( "begTimeSec", &begtimeSec , "begTimeSec/D" ) ;
    fRegTree->Branch( "endTimeMJD", &begtimeMJD , "endTimeMJD/I" ) ;
    fRegTree->Branch( "endTimeSec", &begtimeSec , "endTimeSec/D" ) ;
    fRegTree->Branch( "state"     , &state      , "state/s" ) ;
    
    // loop over all telescopes, pixels, and gains
    VNTelescope* vnt ;
    VNPixel*      vnp ;
    VNGain*       vng ;
    runnum = getRunNumber() ;
    
    // our loop iterators
    MapOfVNTelescopes::iterator  telIter ;
    MapOfVNPixels    ::iterator  pixIter ;
    MapOfVNGains     ::iterator gainIter ;
    vector<VNStateDuration>::iterator histIter ;
    
    // loop over each telescope
    for( telIter  = telescopeMap.begin() ; telIter != telescopeMap.end() ; ++telIter )
    {
    
        vnt   = &telIter->second ;
        telid = vnt->getTelid()  ;
        
        // loop over each pixel
        for( pixIter = vnt->pixelMap.begin() ; pixIter != vnt->pixelMap.end() ; ++pixIter )
        {
        
            vnp     = &pixIter->second     ;
            pixid   = vnp->pixid_dbtables  ;
            pixXcam = vnp->pixelCoord.xcam ;
            pixYcam = vnp->pixelCoord.ycam ;
            
            // loop over each gain type
            for( gainIter = vnp->gainMap.begin() ; gainIter != vnp->gainMap.end() ; ++gainIter )
            {
            
                vng      = &gainIter->second ;
                gainChar = vng->fGainChar    ;
                
                // loop over each row of history
                for( histIter = vng->stateHistory.begin() ; histIter != vng->stateHistory.end() ; ++histIter )
                {
                
                    // load the variables describing the state and its duration
                    state      = ( *histIter ).state        ;
                    begtimeMJD = ( *histIter ).begtime.mjd  ;
                    begtimeSec = ( *histIter ).begtime.time ;
                    endtimeMJD = ( *histIter ).endtime.mjd  ;
                    endtimeSec = ( *histIter ).endtime.time ;
                    
                    // add it to our tree if the pixel duration is not functional
                    if( state != 0 )
                    {
                        fRegTree->Fill() ;
                        
                        // write our extra csv file for double checking
                        if( fWriteCSV )
                        {
                            fprintf( pFile, "%9d, %5d, %5d, %4c, %10f, %10f, %6d, %14f, %6d, %14f, %d\n",
                                     runnum, telid, pixid, gainChar, pixXcam, pixYcam,
                                     begtimeMJD, begtimeSec, endtimeMJD, endtimeSec, state ) ;
                        }
                        
                    } // endif: state was fine(0)
                } // endfor: no more entries in the state history
            } // endfor: no more gains to loop over
        } // endfor: no more pixels to loop over
    } // endfor: no more telescopes to loop over
    
    cout << coutprefix << "fRegTree now has '" << fRegTree->GetEntries() << "' entries..." << endl;
    
    if( fWriteCSV )
    {
        fclose( pFile ) ;
    }
    
}


// do all the stuff to finish up the VDeadPixelOrganizer class's objectives
void VDeadPixelOrganizer::finalize()
{
    // close up all open durations in all histories
    topOffRows() ;
    
    // organize data structs into TTree
    organize() ;
    
    // write TTree to output evndisp.root file
    getTree()->Write() ;
}

