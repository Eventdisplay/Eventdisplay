//
// VDeadPixelOrganizer
// top-level class for taking in pixel information,
// converting it to a smaller format,
// and writing it to a TTree, which is passed up to
// the anasum-ctools event converter
//
// Nathan Kelley-Hoskins
// Aug 2014

#include <ostream>
#include <iostream>
#include <map>
#include <vector>

#include <stdlib.h>
#include <stdio.h>

#include <TTree.h>

#include <VDetectorGeometry.h>
#include <VEvndispRunParameter.h>

using namespace std ;

#define coutprefix "NKH "
#define tellimit 2
#define pixLlimit 70
#define pixUlimit 90
class VDeadPixelOrganizer ;
class VNTelescope         ;
class VNPixel             ;

class VNTime
{
    public :
    
        // variables
        double mjd  ;
        double time ;
        
        // functions
        VNTime( int mjd_days = -1, double time_seconds = -1.0 ) ;
        VNTime halfwayTime( const VNTime vnt ) ;
        
        friend ostream& operator<<( ostream& os, const VNTime& vnt )
        {
            os << "MJD" << vnt.mjd << ":" << ( int )vnt.time << "s" ;
            return os ;
        }
        friend ostream& operator<<( ostream& os, const VNTime* vnt )
        {
            os << ( *vnt ) ;
            return os ;
        }
        
} ;

// Everywhere you see PixelStateInt, this refers to
// an unsigned int, as the pixel state may be changed
// to a bitmask or something more complicated later
typedef unsigned int PixelStateInt ;

// pixel state at an instant in time
class VNStateInstant
{
    public :
        VNTime        time  ;
        PixelStateInt state ;
} ;

// pixel state for a duration of time
class VNStateDuration
{
    public :
        VNTime        begtime ;
        VNTime        endtime ;
        PixelStateInt state   ; // 0=pixel is ok, any other means pixel was having problems
} ;

class VNCameraCoordinate
{
    public :
        double xcam ;  // x coordinate of pixel in camera coordinates (deg)
        double ycam ;  // y coordinate of pixel in camera coordinates (deg)
} ;

class VNEvent
{
    public :
        int    fEventNumber ;
        VNTime fEventTime   ;
} ;

class VNGain
{
    public :
        // variables
        bool   fIsLow ; // true if lowgain, false if highGain
        string fGainName ;
        char   fGainChar ;
        
        VNPixel* fPixel ;  // pointer to parent VNPixel
        vector<VNStateDuration> stateHistory ;
        
        // functions
        VNGain() {} ;
        void initialize( string gainType, VNPixel* parentVNP ) ;
        void updateState( int mjd, double time, PixelStateInt state ) ;
        void setParentVNPixelPointer( VNPixel* ptr ) ;
        bool isLowGain()
        {
            return  fIsLow ;
        } ;
        bool isHighGain()
        {
            return !fIsLow ;
        } ;
        bool isRealObject() ; // for testing if this is a real object or an empty copy
        
        // cout << VNGain
        friend ostream& operator<<( ostream& os, const VNGain& vng )
        {
            cout << coutprefix << "======== VNGain " << vng.fGainName << " : " << vng.stateHistory.size() << " rows of history" << endl;
            return os ;
        }
        
        // cout << * VNGain
        friend ostream& operator<<( ostream& os, const VNGain* vng )
        {
            os << ( *vng ) ;
            return os ;
        }
        
} ;

typedef map<string, VNGain> MapOfVNGains ;
class VNPixel
{
    public :
    
        // variables
        int pixid_evndisp  ; // 0-498
        int pixid_dbtables ; // 1-499
        VNCameraCoordinate pixelCoord ;
        MapOfVNGains gainMap ; // "high" or "low"
        VNTelescope* fTelescope ;  // this pixel's parent telescope
        
        // functions
        VNPixel() {} ;
        void initialize( int pixid, VNTelescope* parentVNT ) ;
        void setPixid( int pid ) ;
        int  getPixid( )
        {
            return pixid_dbtables ;
        } ;
        void setupPixelCoord( VDetectorGeometry* fDetectorGeo ) ;
        void setupGainMap() ;
        void setParentVNTelescopePointer( VNTelescope* ptr ) ;
        bool isRealObject() ;
        //void updateState( bool lowGain, int mjd, double ti, PixelStateInt st ) ;
        
        // cout << VNPixel
        friend ostream& operator<<( ostream& os, const VNPixel& vnp )
        {
            os << coutprefix << "====== VNPixel id=" << vnp.pixid_dbtables ;
            os << " ==========" << endl ;
            MapOfVNGains::const_iterator gainIter ;
            for( gainIter = vnp.gainMap.begin() ; gainIter != vnp.gainMap.end() ; ++gainIter )
            {
                os << gainIter->second ;
            }
            return os ;
        }
        
        // cout << * VNPixel
        friend ostream& operator<<( ostream& os, const VNPixel* vnp )
        {
            os << ( *vnp ) ;
            return os ;
        }
        
} ;

typedef map<int, VNPixel> MapOfVNPixels ;
class VNTelescope
{
    public :
    
        // variables
        MapOfVNPixels pixelMap ; // 1-499
        VDeadPixelOrganizer* fOrganizer ;  // pointer to parent class
        
        // functions
        VNTelescope() {} ;
        void initialize( int telid, int npix, VDeadPixelOrganizer* parentVDPO ) ;
        void setNpix( int npixels )
        {
            npix = npixels ;
        } ;
        int  getNpix( )
        {
            return npix    ;
        } ;
        void setTelid( int tid )
        {
            fTelID = tid   ;
        } ;
        int  getTelid( )
        {
            return fTelID  ;
        } ;
        void setParentVDeadPixelOrganizerPointer( VDeadPixelOrganizer* vdpo ) ;
        void setupMap() ;
        bool isRealObject() ;
        
        // cout << VNTelescope
        friend ostream& operator<<( ostream& os, const VNTelescope& vnt )
        {
            os << coutprefix << "==== VNTelescope id=" << vnt.fTelID << " ===========" << endl;
            MapOfVNPixels::const_iterator pixIter ;
            for( pixIter = vnt.pixelMap.begin() ; pixIter != vnt.pixelMap.end() ; ++pixIter )
            {
                os << pixIter->second ;
            }
            return os ;
        }
        
        // cout << *VNTelescope
        friend ostream& operator<<( ostream& os, const VNTelescope* vnt )
        {
            os << ( *vnt ) ;
            return os ;
        }
        
    private :
    
        // variables
        int fTelID ;
        int npix   ;
} ;

typedef map<int, VNTelescope> MapOfVNTelescopes ;
class VDeadPixelOrganizer
{
    public :
    
        // variables
        VNTime  fRunStart      ;
        VNTime  fRunEnd        ;
        VNEvent fPreviousEvent ; // previous event/time
        MapOfVNTelescopes   telescopeMap ;
        VDetectorGeometry* fDetectorGeo ;  // for getting pixel numbers and coordinates
        bool fWriteCSV ;
        
        // functions
        VDeadPixelOrganizer( int ntel = 0, int npix = 0, VDetectorGeometry* detGeo = 0, int startMJD = 0, double startSec = 0.0, int endMJD = 0, double endSec = 0.0, string treename = "deadPixReg", int runNumber = 0 ) ;
        void     setNtel( int    ntelescopes )
        {
            ntel = ntelescopes     ;
        } ;
        int      getNtel( )
        {
            return ntel            ;
        } ;
        void     setNpix( int    npixels )
        {
            npix = npixels         ;
        } ;
        int      getNpix( )
        {
            return npix            ;
        } ;
        void     setRunNumber( int    runnumber )
        {
            fRunNumber = runnumber ;
        } ;
        int      getRunNumber( )
        {
            return fRunNumber      ;
        } ;
        void     setTreeName( string treename )
        {
            fTreeName = treename   ;
        } ;
        string   getTreeName( )
        {
            return fTreeName       ;
        } ;
        TTree*   getTree( )
        {
            return fRegTree        ;
        } ;
        void                setDetectorGeo( VDetectorGeometry* detGeo ) ;
        VDetectorGeometry* getDetectorGeo( ) ;
        void updatePreviousEventInfo( int eventNumber, int eventMJD, double eventTime ) ;
        void setupMap() ; // call this to initialize object properly
        // call set_ntel() and set_npix() before this, though
        
        void UpdatePixelState( int tel, int pix, bool lowGain, int mjd, double time, PixelStateInt state ) ;
        void topOffRows() ;
        void printSummary() ;
        void organize() ;
        void finalize() ;     // call to write tree to file
        
        // cout << VDeadPixelOrganizer
        friend ostream& operator<<( ostream& os, const VDeadPixelOrganizer& dpo )
        {
            //os << "blaaah" << endl;
            os << coutprefix << "== VDeadPixelOrganizer ==========" << endl ;
            MapOfVNTelescopes::const_iterator telIter ;
            for( telIter  = dpo.telescopeMap.begin() ; telIter != dpo.telescopeMap.end() ; ++telIter )
            {
                os << coutprefix << telIter->second ;
            }
            
            return os ;
        }
        
        // cout << * VDeadPixelOrganizer
        friend ostream& operator<<( ostream& os, const VDeadPixelOrganizer* dpo )
        {
            os << ( *dpo ) ;
            return os ;
        }
        
    private :
    
        // variables
        int      ntel       ;
        int      npix       ;
        int      fRunNumber ;
        string   fTreeName  ;
        TTree*   fRegTree   ;
        
} ;

