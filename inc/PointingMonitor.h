/**
 * \file
 * \brief Function declarations for Pointing Monitor interactions with the database
 *
 *
 * \author Dirk Pandel
 */

#ifndef POINTINGMONITOR_H
#define POINTINGMONITOR_H

#include "VGlobalRunParameter.h"
#include "VDB_Connection.h"

#include <stdint.h>
#include <vector>

#include <TMath.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSQLServer.h>

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>

#include <cmath>
#include <ctime>

using namespace std;

namespace pointingmonitor
{

    /// Structure holding the calibrated LED positions and the rotation and scale of the PMTs and the sky camera
    struct CameraParameters
    {
        bool isValid;                        ///< True if the object contains valid parameters
        uint32_t startDate;                  ///< Start date from which the parameters are valid (yyyymmdd)
        uint32_t endDate;                    ///< End date until which the parameters are valid (yyyymmdd)
        uint32_t version;                    ///< Version number of the parameters
        //JG/////std::vector<LEDPos> ledPositions;    ///< Vector of LED positions in the PMT coordinate system
        float pmtRotation;                   ///< Clockwise PMT rotation angle vs. the horizon (radians)
        float pmtScale;                      ///< Plate scale of the PMTs (radians/PMT)
        float skyCameraRotation;             ///< Clockwise rotation angle of the sky camera (radians)
        float skyCameraScale;                ///< Plate scale of the sky camera (radians/pixel)
        float referencePixelX;               ///< X-position of reference CCD pixel
        float referencePixelY;               ///< Y-position of reference CCD pixel
        
        /// Constructor: sets all members to 0
        CameraParameters() : isValid( false ), startDate( 0 ), endDate( 0 ), version( 0 ),
            pmtRotation( 0.0 ), pmtScale( 0.0 ),
            skyCameraRotation( 0.0 ), skyCameraScale( 0.0 ),
            referencePixelX( 0.0 ), referencePixelY( 0.0 ) {}
    };
    
    
    /// Structure holding the calibration parameters
    struct CalibrationParameters
    {
        bool isValid;                        ///< True if the object contains valid parameters
        uint32_t startDate;                  ///< Start date from which the parameters are valid (yyyymmdd)
        uint32_t endDate;                    ///< End date until which the parameters are valid (yyyymmdd)
        uint32_t version;                    ///< Version number of the parameters
        float x0;                            ///< Polynomial 1, 0th order
        float x1;                            ///< Polynomial 1, 1st order
        float x2;                            ///< Polynomial 1, 2nd order
        float x3;                            ///< Polynomial 1, 3rd order
        float y0;                            ///< Polynomial 2, 0th order
        float y1;                            ///< Polynomial 2, 1st order
        float y2;                            ///< Polynomial 2, 2nd order
        float y3;                            ///< Polynomial 2, 3rd order
        float z0;                            ///< Polynomial 3, 0th order
        float z1;                            ///< Polynomial 3, 1st order
        float z2;                            ///< Polynomial 3, 2nd order
        float z3;                            ///< Polynomial 3, 3rd order
        float l0;                            ///< Polynomial 4, 0th order
        float l1;                            ///< Polynomial 4, 1st order
        float l2;                            ///< Polynomial 4, 2nd order
        float l3;                            ///< Polynomial 4, 3rd order
        
        /// Constructor: sets all members to 0
        CalibrationParameters() : isValid( false ), startDate( 0 ), endDate( 0 ), version( 0 ),
            x0( 0.0 ), x1( 0.0 ), x2( 0.0 ), x3( 0.0 ),
            y0( 0.0 ), y1( 0.0 ), y2( 0.0 ), y3( 0.0 ),
            z0( 0.0 ), z1( 0.0 ), z2( 0.0 ), z3( 0.0 ),
            l0( 0.0 ), l1( 0.0 ), l2( 0.0 ), l3( 0.0 ) {}
    };
    
    
    /// Structure holding the MJD(UTC) time and the uncalibrated pointing
    struct UncalibratedPointing
    {
        double mjd;                          ///< Modified Julian Date
        float ra;                            ///< Right ascension (radians)
        float dec;                           ///< Declination (radians)
        float rotation;                      ///< Rotation angle (radians)
        float elevation;                     ///< Elevation from positioner (radians)
        float ledPosY;                       ///< Average y-position of LEDs (pixel)
        
        /// Constructor: sets all members to 0
        UncalibratedPointing() : mjd( 0.0 ), ra( 0.0 ), dec( 0.0 ), rotation( 0.0 ),
            elevation( 0.0 ), ledPosY( 0.0 ) {}
            
        /// Constructor: initializes all members
        UncalibratedPointing( double mjd, float ra, float dec, float rotation,
                              float elevation, float ledPosY )
            : mjd( mjd ), ra( ra ), dec( dec ), rotation( rotation ),
              elevation( elevation ), ledPosY( ledPosY ) {}
    };
    
    
    /// Structure holding the MJD(UTC) time and the calibrated pointing
    struct CalibratedPointing
    {
        double mjd;                          ///< Modified Julian Date
        float ra;                            ///< Right ascension (radians)
        float dec;                           ///< Rotation angle (radians)
        float rotation;                      ///< Rotation angle (radians)
        
        /// Constructor: sets all members to 0
        CalibratedPointing() : mjd( 0.0 ), ra( 0.0 ), dec( 0.0 ), rotation( 0.0 ) {}
        
        /// Constructor: initializes all members
        CalibratedPointing( double mjd, float ra, float dec, float rotation )
            : mjd( mjd ), ra( ra ), dec( dec ), rotation( rotation ) {}
    };
    
    
    
    /********************************************************************************/
    /* Functions for converting between Modified Julian Date and calender date/time */
    /********************************************************************************/
    
    class PointingMonitor : public VGlobalRunParameter
    {
        public:
        
            /**
             * \brief Convert Modified Julian Date to a date integer
             *
             * Converts MJD to a date integer of the form YYYYMMDD.
             * Fractions of a day are rounded down.
             * \param mjd          Modified Julian Date
             */
            uint32_t mjdToDate( double mjd );
            
            /**
             * \brief Convert a date/time integer to Modified Julian Date
             *
             * Converts a date/time integer to MJD.
             * \param dateTime     Date and time in the form YYYYMMDDhhmmss
             */
            double dateTimeToMjd( uint64_t dateTime );
            
            uint32_t dateToUInt32( const string& dateStr );
            
            /************************************************************************/
            /* Functions for writing/reading camera parameters to/from the database */
            /************************************************************************/
            
            /**
             * \brief Read sets of camera parameters for one telescope from the database
             *
             * Reads the calibrated LED positions and the rotation and scale
             * of the PMTs and the sky camera for a specified date from the
             * tblPointing_Monitor_Camera_Parameters table.
             * The function returns a vector of all database entries matching the date.
             * If date=0 or no date is given, all entries for a telescope are returned.
             * Entries are sorted in descending order of the version number.
             * \param telescope_id Telescope ID (0-3)
             * \param date         Date (yyyymmdd)
             * \param limit        Maximum number of entries to return (no limit if 0)
             */
            vector<CameraParameters>
            getCameraParametersList( uint32_t telescope_id, uint32_t date = 0, uint32_t limit = 0 );
            
            /**
             * \brief Read latest set of camera parameters for a specified date from the database
             *
             * Reads the calibrated LED positions and the rotation and scale
             * of the PMTs and the sky camera for a specified date from the
             * tblPointing_Monitor_Camera_Parameters table.
             * The function returns the database entry with the highest version number.
             * \param telescope_id Telescope ID (0-3)
             * \param date         Date (yyyymmdd)
             */
            CameraParameters getCameraParameters( uint32_t telescope_id, uint32_t date );
            
            
            /*****************************************************************************/
            /* Functions for writing/reading calibration parameters to/from the database */
            /*****************************************************************************/
            
            /**
             * \brief Read sets of calibration parameters for one telescope from the database
             *
             * Reads the pointing monitor calibration parameters for a specified date from the
             * tblPointing_Monitor_Calibration_Parameters table.
             * The function returns a vector of all database entries matching the date.
             * If date=0 or no date is given, all entries for a telescope are returned.
             * Entries are sorted in descending order of the version number.
             * \param telescope_id Telescope ID (0-3)
             * \param date         Date (yyyymmdd)
             * \param limit        Maximum number of entries to return (no limit if 0)
             */
            vector<CalibrationParameters>
            getCalibrationParametersList( uint32_t telescope_id, uint32_t date = 0, uint32_t limit = 0 );
            
            
            /**
             * \brief Read latest set of calibration parameters for a specified date from the database
             *
             * Reads the pointing monitor calibration parameters for a specified date from the
             * tblPointing_Monitor_Calibration_Parameters table.
             * The function returns the database entry with the highest version number.
             * \param telescope_id Telescope ID (0-3)
             * \param date         Date (yyyymmdd)
             */
            CalibrationParameters getCalibrationParameters( uint32_t telescope_id, uint32_t date );
            
            
            /****************************************************************************/
            /* Functions for writing/reading uncalibrated pointing to/from the database */
            /****************************************************************************/
            
            /**
             * \brief Read uncalibrated pointings for a given time interval from the database
             *
             * Reads a vector of UncalibratedPointing objects for a given time interval
             * from the tblPointing_Monitor_TelescopeX_Pointing tables.
             * \param telescope_id Telescope ID (0-3)
             * \param startmjd     Start time in MJD(UTC)
             * \param stopmjd      Stop time in MJD(UTC)
             */
            vector<UncalibratedPointing>
            getUncalibratedPointing( uint32_t telescope_id, double startmjd, double stopmjd );
            
            
            /*****************************************************/
            /* Functions for calculating the calibrated pointing */
            /*****************************************************/
            
            
            /**
             * \brief Calculates the calibrated pointing
             *
             * Calculates a vector of calibrated pointings
             * from a vector of uncalibrated pointings
             * and a set of camera and calibration parameters.
             * \param pointingVec    Vector of UncalibratedPointing objects
             * \param camParameters  CameraParameters object
             * \param calParameters  CalibrationParameters object
             */
            vector<CalibratedPointing>
            calibratedPointing( const vector<UncalibratedPointing>& pointingVec,
                                const CameraParameters& camParameters,
                                const CalibrationParameters& calParameters );
                                
            PointingMonitor();
            ~PointingMonitor() {}
    };
    
}   // end of namespace

#endif
