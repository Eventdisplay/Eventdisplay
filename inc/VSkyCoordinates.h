//! VSkyCoordinates get pointing direction of telescope

#ifndef VSKYCOORDINATES_H
#define VSKYCOORDINATES_H

#include "TMath.h"

#include <iomanip>
#include <iostream>

#include "VSkyCoordinatesUtilities.h"
#include "VStarCatalogue.h"

using namespace std;

class VSkyCoordinates
{
    protected:
    
        VStarCatalogue* fStarCatalogue;
        
        bool   fSet;                              //!< true if target or dec/ra is set
        bool   fPrecessed;                        //!< true if target position has been precessed
        bool   fMC;                               //!< true for Monte Carlo run
        
        unsigned int fMJD;
        double fTime;
        double fTelDec;                           //!< [rad] //! declination of pointing direction
        double fTelRA;                            //!< [rad] //! right ascension of pointing direction
        double fTelDecJ2000;                     //!< [rad] //! declination of pointing direction (J2000)
        double fTelRAJ2000;                      //!< [rad] //! right ascension of pointing direction (J2000)
        double fTargetDec;                        //!< [rad] //! declination of target
        double fTargetRA;                         //!< [rad] //! right ascension of target
        double fTargetDecJ2000;                        //!< [rad] //! declination of target (J2000)
        double fTargetRAJ2000;                         //!< [rad] //! right ascension of target (J2000)
        string fTargetName;
        double fTargetElevation;                  //!< [deg]
        double fTargetAzimuth;                    //!< [deg]
        bool   fWobbleSet;                        //!< be sure that wobble offset is only applied once
        double fWobbleNorth;                      //!< [deg] wobble offset north
        double fWobbleEast;                       //!< [deg] wobble offset east
        
        // telescope orientation
        double fTelAzimuth;                       //!< [deg]  return value to be used in the analysis
        double fTelElevation;                     //!< [deg]  return value to be used in the analysis
        float  fTelAzimuthCalculated;             //!< [deg]  elevation from source coordinates
        float  fTelElevationCalculated;           //!< [deg]  elevation from source coordinates
        
        double fObsLatitude;                      //!< [rad]
        double fObsLongitude;                     //!< [rad]
        double fSupressStdoutText ;
        
        void reset();
        
    public:
    
        VSkyCoordinates();
        ~VSkyCoordinates();
        
        double derotateCoords( int MJD, double time, double i_xin, double i_yin, double& i_xout, double& i_yout );
        double derotateCoords( double i_UTC, double i_xin, double i_yin, double& i_xout, double& i_yout );
        void   getEquatorialCoordinates( int MJD, double time, double az, double ze, double& dec, double& ra );
        string getTargetName()
        {
            return fTargetName;
        }
        double getDerotationAngle( int MJD, double time );
        void   getDerotatedShowerDirection( double ze, double az, float& y, float& x, double rze, double raz );
        double getTargetDec()
        {
            return fTargetDec * TMath::RadToDeg();
        }
        double getTargetRA()
        {
            return fTargetRA * TMath::RadToDeg();
        }
        double getTargetDecJ2000()
        {
            return fTargetDecJ2000 * TMath::RadToDeg();
        }
        double getTargetRAJ2000()
        {
            return fTargetRAJ2000 * TMath::RadToDeg();
        }
        double getTargetElevation()
        {
            return fTargetElevation;
        }
        double getTargetAzimuth()
        {
            return fTargetAzimuth;
        }
        double getTelAzimuth()
        {
            return fTelAzimuth;
        }
        double getTelElevation()
        {
            return fTelElevation;
        }
        double getTelDec()
        {
            return fTelDec;
        }
        double getTelRA()
        {
            return fTelRA;
        }
        double getTelLatitude()
        {
            return fObsLatitude * TMath::RadToDeg();
        }
        double getTelLongitude()
        {
            return fObsLongitude * TMath::RadToDeg();
        }
        VStarCatalogue* getStarCatalogue()
        {
            return fStarCatalogue;
        }
        double getWobbleNorth()
        {
            return fWobbleNorth;
        }
        double getWobbleEast()
        {
            return fWobbleEast;
        }
        unsigned int getMJD()
        {
            return fMJD ;
        }
        double getTime()
        {
            return fTime ;
        }
        bool   initStarCatalogue( string iCatalogueName, double iMJD, double xmin, double xmax, double ymin, double ymax,
                                  double iRASkyMapCentreJ2000, double iDecSkyMapCentreJ2000 );
        bool   isPrecessed()
        {
            return fPrecessed;
        }
        bool   isSet()
        {
            return fSet;
        }
        void   precessTarget( int iMJD, int iTelID = -1 );
        double rotateCoords( int i_mjd, double i_seconds, double i_xin, double i_yin, double& i_xout, double& i_yout );
        void   setMC()
        {
            fMC = true;
        }
        void   setObservatory( double iLongitude_deg = 0., double iLatitude_deg = 0. );
        bool   setPointingOffset( double i_raOff, double i_decOff );
        bool   setTargetJ2000( double iDec_deg, double iRA_deg );
        void   setTargetName( string iTargetName )
        {
            fTargetName = iTargetName;
        }
        void   setTelDec_deg( double iTelDec_deg )
        {
            fTelDec = iTelDec_deg * TMath::DegToRad();
        }
        void   setTelRA_deg( double iTelRA_deg )
        {
            fTelRA  = iTelRA_deg * TMath::DegToRad();
        }
        void   setTelAzimuth( double iTelAz )
        {
            fTelAzimuth = iTelAz;    //!< set telescope azimuth (e.g.for MC)
        }
        void   setTelElevation( double iTelEl )
        {
            fTelElevation = iTelEl;    //!< set telescope elevation (e.g. for MC)
        }
        void   setWobbleOffset( double iWobbleNorth, double iWobbleEast, int iTelID, int iMJD );
        void   updatePointing( int MJD, double time );
        
        // for hiding some text (so other output text shows up cleanly)
        // if setting = true, will hid some stdout text
        // if setting = false, class will behave normally
        void   supressStdoutText( bool setting )
        {
            fSupressStdoutText = setting ;
        }
};
#endif
