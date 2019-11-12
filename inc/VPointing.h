//! VPointing get pointing direction of telescope

#ifndef VPointing_H
#define VPointing_H

#include "TMath.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>

#include "VAstronometry.h"
#include "VPointingDB.h"
#include "VSkyCoordinates.h"
#include "VSkyCoordinatesUtilities.h"

using namespace std;

class VPointing : public VSkyCoordinates
{
    private:
    
        unsigned int fTelID;
        bool         fUseDB;                      //!< uses DB to calculate pointing directions
        
        unsigned int fPointingType;               //!< 0: pointing calculated from source coordinates (+wobble offsets)
        //!< 1: pointing calculated from source coordinates (+wobble offsets), added error from command line,
        //!< 2: read T-Point corrected positioner data from VERITAS DB
        //!< 3: read raw positioner data from VERITAS DB and apply tracking corrections
        //!< 4: from pointing monitor (text file)
        //!< 5: from pointing monitor (DB)
        unsigned int fEventStatus;
        unsigned int fNEventsWithNoDBPointing;
        
        float  fTelAzimuthDB;                     //!< [deg]  azimuth from VTS DB (from positioner or pointing monitor)
        float  fTelElevationDB;                   //!< [deg]  elevation from VTS DB (from positioner or pointing monitor)
        // difference between DB/VPM pointing and position
        float fPointingErrorX;                    //!< [deg]
        float fPointingErrorY;                    //!< [deg]
        unsigned int fMeanPointingErrorN;
        double fMeanPointingErrorX;               //!< [deg]
        double fMeanPointingErrorY;               //!< [deg]
        double fMeanPointingDistance;             //!< [deg]
        
        VPointingDB* fPointingDB;
        TTree* fPointingTree;
        
        // private functions
        void initializePointingTree();
        bool updatePointingfromDB( int, double );
        
    public:
    
        VPointing( unsigned int itelID );
        ~VPointing() {}
        
        void         fillPointingTree();
        float        getPointingErrorX()
        {
            return fPointingErrorX;
        }
        float        getPointingErrorY()
        {
            return fPointingErrorY;
        }
        unsigned int getPointingType()
        {
            return fPointingType;
        }
        void         getPointingFromDB( int irun, string iTCorrections, string iVPMDirectory, bool iVPMDB, bool iUncalibratedVPM );
        unsigned int getTelID()
        {
            return fTelID;
        }
        void         setPointingError( double, double );//!< Pointing error [deg]
        void         setTelPointing( int MJD, double time, bool iUseDB = false, bool iFillPointingTree = false );
        void         terminate( bool i_IsMC = false );
};
#endif
