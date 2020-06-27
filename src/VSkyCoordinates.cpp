/*! \class VSkyCoordinates
    \brief handle telescope pointing direction and target positions

*/

#include "VSkyCoordinates.h"

VSkyCoordinates::VSkyCoordinates()
{
    fStarCatalogue = 0;
    
    reset();
    setObservatory();
}


VSkyCoordinates::~VSkyCoordinates()
{
    if( fStarCatalogue )
    {
        delete fStarCatalogue;
    }
}

void VSkyCoordinates::setObservatory( double iLongitude, double iLatitude )
{
    fObsLatitude  = iLatitude * TMath::DegToRad();
    fObsLongitude = iLongitude * TMath::DegToRad();
}

void VSkyCoordinates::reset()
{
    fMC = false;
    fSet = false;
    fPrecessed = false;
    fWobbleSet = false;
    
    fTelDec = 0.;
    fTelRA  = 0.;
    fTelDecJ2000 = 0.;
    fTelRAJ2000  = 0.;
    fTargetDec = 0.;
    fTargetRA  = 0.;
    fTargetDecJ2000 = 0.;
    fTargetRAJ2000  = 0.;
    fTelAzimuth = 0.;
    fTelElevation = 0.;
    fTelAzimuthCalculated = 0.;
    fTelElevationCalculated = 0.;
    
    fTargetAzimuth = 0.;
    fTargetElevation = 0.;
    fWobbleNorth = 0.;
    fWobbleEast = 0.;
    
    fMJD = 0;
    fTime = 0.;
    
    fSupressStdoutText = false ;
}

void VSkyCoordinates::precessTarget( int iMJD, int iTelID )
{
    if( !fPrecessed )
    {
        if( !fSupressStdoutText )
        {
            cout << "---------------------------------------------------------------------------------------------------------" << endl;
        }
        if( iTelID >= 0 )
        {
            if( !fSupressStdoutText )
            {
                cout << "Pointing telescope " << iTelID + 1 << endl;
            }
        }
        else
        {
            cout << "Array pointing " << endl;
        }
        if( !fSupressStdoutText )
        {
            cout << "\tPrecessing target ( " << getTargetName() << " ) from J2000 to MJD " << iMJD << endl;
            // TEMP
            cout << "\tJ2000   \t\t RA=" << setprecision( 6 ) << fTargetRA* TMath::RadToDeg() << " dec=" << fTargetDec* TMath::RadToDeg() << endl;
        }
        // ENDTEMP
        // precess target coordinates
        VSkyCoordinatesUtilities::precessTarget( iMJD, fTargetRA, fTargetDec );
        
        if( !fSupressStdoutText )
        {
            cout << "\tMJD " << iMJD;
            cout << "\t RA=" << fTargetRA* TMath::RadToDeg() << " dec=" << fTargetDec* TMath::RadToDeg() << endl;
        }
        // precess telescope coordinates
        VSkyCoordinatesUtilities::precessTarget( iMJD, fTelRA, fTelDec );
        
        fPrecessed = true;
    }
}


/*!

    set target coordinates in J2000

    \param iDec  [deg]
    \param iRA   [deg]
*/
bool VSkyCoordinates::setTargetJ2000( double iDec, double iRA )
{
    fTargetDecJ2000 = iDec / TMath::RadToDeg();
    fTargetRAJ2000  = iRA / TMath::RadToDeg();
    
    // unprecessed -> precess later in the analysis
    fTargetDec = iDec / TMath::RadToDeg();
    fTargetRA  = iRA / TMath::RadToDeg();
    
    fTelDecJ2000 = iDec / TMath::RadToDeg();
    fTelRAJ2000  = iRA / TMath::RadToDeg();
    // unprecessed -> precess later in the analysis
    fTelDec = iDec / TMath::RadToDeg();
    fTelRA  = iRA / TMath::RadToDeg();
    
    fSet = true;
    
    return true;
}

/*

    calculate azimuth and elevation of telescope pointing direction
    calculate azimuth and elevation of target
 */
void VSkyCoordinates::updatePointing( int MJD, double time )
{
    fMJD = ( unsigned int )MJD;
    fTime = time;
    
    double az = 0.;
    double el = 0.;
    
    // telescope elevation/azimuth calculated from source coordinates and time
    VSkyCoordinatesUtilities::getHorizontalCoordinates( MJD, time, fTelDec * TMath::RadToDeg(), fTelRA * TMath::RadToDeg(), az, el );
    el = 90. - el;
    fTelAzimuthCalculated   = ( float )az;
    fTelElevationCalculated = ( float )el;
    fTelElevation = fTelElevationCalculated;
    fTelAzimuth   = fTelAzimuthCalculated;
    
    // set target azimuth/elevation
    VSkyCoordinatesUtilities::getHorizontalCoordinates( MJD, time, fTargetDec * TMath::RadToDeg(), fTargetRA * TMath::RadToDeg(), fTargetAzimuth, fTargetElevation );
    fTargetElevation = 90. - fTargetElevation;
}

/*

    calculate right ascension / declination

    all angles (in/out) in [deg]

*/
void VSkyCoordinates::getEquatorialCoordinates( int MJD, double time, double az, double ze, double& dec, double& ra )
{
    if( fMC )
    {
        dec = fTelDec;
        ra  = fTelRA;
        return;
    }
    VSkyCoordinatesUtilities::getEquatorialCoordinates( MJD, time, az, ze, dec, ra );
}

/*

   add an offset in ra/dec

   this should happen before precession is applied

*/
bool VSkyCoordinates::setPointingOffset( double i_raOff_deg, double i_decOff_deg )
{
    if( isPrecessed() )
    {
        fTelRA = -9999.;
        fTelDec = -9999.;
        fTelRAJ2000 = -9999.;
        fTelDecJ2000 = -9999.;
        return false;
    }
    
    fTelRA      = fTargetRA + i_raOff_deg * TMath::DegToRad();
    fTelDec     = fTargetDec + i_decOff_deg * TMath::DegToRad();
    
    fTelRAJ2000      = fTargetRAJ2000 + i_raOff_deg * TMath::DegToRad();
    fTelDecJ2000     = fTargetDecJ2000 + i_decOff_deg * TMath::DegToRad();
    
    return true;
}


double VSkyCoordinates::derotateCoords( double i_UTC, double i_xin, double i_yin, double& i_xout, double& i_yout )
{
    double i_theta = VSkyCoordinatesUtilities::getDerotationAngle( i_UTC, fTelRA, fTelDec, fObsLongitude, fObsLatitude );
    i_xout = i_xin * cos( i_theta ) + i_yin * sin( i_theta );
    i_yout = i_yin * cos( i_theta ) - i_xin * sin( i_theta );
    return i_theta;
}

double VSkyCoordinates::getDerotationAngle( int i_mjd, double i_seconds )
{
    return VSkyCoordinatesUtilities::getDerotationAngle( VSkyCoordinatesUtilities::getUTC( i_mjd, i_seconds ),
            fTelRA, fTelDec, fObsLongitude, fObsLatitude );
}

double VSkyCoordinates::derotateCoords( int i_mjd, double i_seconds, double i_xin, double i_yin, double& i_xout, double& i_yout )
{
    double i_UTC = VSkyCoordinatesUtilities::getUTC( i_mjd, i_seconds );
    double i_theta = VSkyCoordinatesUtilities::getDerotationAngle( i_UTC, fTelRA, fTelDec, fObsLongitude, fObsLatitude );
    i_xout = i_xin * cos( i_theta ) + i_yin * sin( i_theta );
    i_yout = i_yin * cos( i_theta ) - i_xin * sin( i_theta );
    return i_theta;
}

double VSkyCoordinates::rotateCoords( int i_mjd, double i_seconds, double i_xin, double i_yin, double& i_xout, double& i_yout )
{
    double i_UTC = VSkyCoordinatesUtilities::getUTC( i_mjd, i_seconds );
    double i_theta = -1. * VSkyCoordinatesUtilities::getDerotationAngle( i_UTC, fTelRA, fTelDec, fObsLongitude, fObsLatitude );
    i_xout = i_xin * cos( i_theta ) + i_yin * sin( i_theta );
    i_yout = i_yin * cos( i_theta ) - i_xin * sin( i_theta );
    return i_theta;
}


/*!
 *
 *  should be called after precession, etc.
 *
 */
void VSkyCoordinates::setWobbleOffset( double iNorth, double iEast, int iTelID, int iMJD )
{
    fWobbleNorth = iNorth;
    fWobbleEast = iEast;
    
    if( !fWobbleSet )
    {
        double i_decDiff = 0.;
        double i_RADiff = 0.;
        VSkyCoordinatesUtilities::getWobbleOffset_in_RADec( iNorth, iEast, fTargetDec * TMath::RadToDeg(), fTargetRA * TMath::RadToDeg(), i_decDiff, i_RADiff );
        if( i_RADiff < -180. )
        {
            i_RADiff += 360.;
        }
        
        fTelRA  = fTargetRA + i_RADiff * TMath::DegToRad();
        fTelDec = fTargetDec + i_decDiff * TMath::DegToRad();
        
        if( iTelID >= 0 )
        {
            cout << "\tWobble mode, telescope " << iTelID + 1;
        }
        else
        {
            cout << "\tWobble mode, array ";
        }
        cout << " pointing to (ra,dec) = (" << fTelRA* TMath::RadToDeg() << ", " << fTelDec* TMath::RadToDeg() << ")";
        cout << ", (delta ra, delta dec) = (" << i_RADiff << ", " << i_decDiff << ")";
        cout << endl;
        
        // set J2000 telescope coordinates
        fTelRAJ2000 = fTelRA;
        fTelDecJ2000 = fTelDec;
        VSkyCoordinatesUtilities::precessTarget( 51544., fTelRAJ2000, fTelDecJ2000, iMJD );
        
        fWobbleSet = true;
    }
    else
    {
        cout << "VSkyCoordinates::setWobbleOffset warning, wobble offsets already set" << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
}

/*

   initial star catalogue and set FOV in catalogue

   results in a list of stars in the current FOV

   expect J2000 coordinates

*/
bool VSkyCoordinates::initStarCatalogue( string iCatalogueName, double iMJD,
        double xmin, double xmax, double ymin, double ymax,
        double iRASkyMapCentre_J2000, double iDecSkyMapCentre_J2000 )
{
    if( !fStarCatalogue )
    {
        fStarCatalogue = new VStarCatalogue();
    }
    if( fStarCatalogue )
    {
        double i_x = 0.;
        double i_y = 0.;
        if( fabs( xmin ) > fabs( xmax ) )
        {
            i_x = fabs( xmin );
        }
        else
        {
            i_x = fabs( xmax );
        }
        if( fabs( ymin ) > fabs( ymax ) )
        {
            i_y = fabs( ymin );
        }
        else
        {
            i_y = fabs( ymax );
        }
        if( !fStarCatalogue->init( iMJD, iCatalogueName ) )
        {
            cout << "Error reading star catalogue: " << iCatalogueName << endl;
            cout << "exiting..." << endl;
            return false;
        }
        fStarCatalogue->setFOV( iRASkyMapCentre_J2000, iDecSkyMapCentre_J2000, i_x, i_y, true );
    }
    
    return true;
}
