/* VTimeMask.h
 * ---------------------------------------------
 * Class to facilitate rejecting events based on their time during the run.
 * ---------------------------------------------
 * Created: Mon  2 Mar 00:00:10 2009
 */
#ifndef VTIMEMASK_H
#  define VTIMEMASK_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <cmath>
#include <set>

#include <TROOT.h>
#include <TBits.h>
#include <TVector.h>
#include <TNamed.h>

#include "VSkyCoordinatesUtilities.h"

using namespace std;

class VTimeMask : public TNamed
{
    private:
        Int_t           run_id;                   // Run number
        Double_t        start_time;               // Start of run: decimal MJD in units of seconds
        string          mask_file;                // Location of file containing mask definitions
        
        vector< Bool_t >    mask;                 // One entry for each second of the run
        vector< UInt_t >    checked;              // Number of events checked against mask for each second
        vector< UInt_t >    accepted;             // Number of events allowed by the mask for each second
        vector< UInt_t >    counted;              // Number of events passing gamma and direction cuts
        
        Bool_t          set_success;              // If the mask was set from the mask_file
        Bool_t          override;                 // If we have encountered an error
        Int_t           outside_count;            // Number of events checked which weren't inside the mask time domain
        // i.e. those before the start or after the end
        
        void            initialise( Int_t run_number, Double_t run_startMJD_secs, Double_t run_endMJD_secs, string file_name );
        // Flake out in a dignified manner.
        void            problem( string reason, ostream& terminal = cout );
        
        // Internal time format is exact decimal MJD converted into seconds (<MJD>.<fraction of day> * secs_day)
        static const    Int_t           secs_day    = 24 * 60 * 60;
        // Create internal time format: double VAnaUtils::getUTC(int i_mjd, double i_seconds) returns <MJD>.<fraction of day>
        Double_t        exactMJD_secs( Int_t MJD, Double_t secsMid ) const
        {
            return secs_day * VSkyCoordinatesUtilities::getUTC( MJD, secsMid );
        }
        // Return integer MJD portion of internal time
        Int_t           timeMJD( Double_t internal_time ) const
        {
            return Int_t( floor( internal_time / secs_day ) );
        }
        // Return seconds since UTC midnight
        Double_t        time_secsMid( Double_t internal_time ) const
        {
            return internal_time - Double_t( secs_day ) * timeMJD( internal_time );
        }
        
        // How is the mask set at this internal time?
        Bool_t          checkMaskNow( Double_t now );
        // Make it up as we go along; collate statistics even if there was problem.
        Bool_t          loadMaskNow( Double_t now );
        
    public:
        VTimeMask();                              // Default constructor gives mask of zero size but reserves storage for 1200 seconds
        //~VTimeMask();	// Remove any ROOT objects with pointers
        // To properly initialise use one of the following constructors: default mask is open for all seconds
        // run_startMJD and run_endMJD are integer MJD; run_startTime and run_endTime are seconds since UTC midnight
        VTimeMask( Int_t run_number, Int_t run_startMJD, Double_t run_startTime, Int_t run_endMJD, Double_t run_endTime, string file_name = "" );
        // run_startUTC and run_endUTC are decimal MJD (<MJD>.<fraction of day>) following [double VAnaUtils::getUTC(int i_mjd, double i_seconds)]
        VTimeMask( Int_t run_number, Double_t run_startUTC, Double_t run_endUTC, string file_name = "" );
        
        // Intialise now
        // To read in the 'timeMask' directory from an anasum file,
        // you need to first do
        // VTimeMask * vtm = VTimeMask( runid, runstart, runend ) ;
        // vtm->readObjects( timeMaskDirectory ) ;
        // vtm->setMask( runid, runstart, runend ) ;
        Bool_t      setMask( Int_t run_number, Double_t run_startUTC, Double_t run_endUTC, string file_name = "" );
        Bool_t      setMask();                    // Retrieve user-defined mask from mask_file
        
        Bool_t      setMaskDuringPhaseCuts( Double_t eventUTC );
        
        
        // Is the event allowed to pass the mask?
        // if true, then that instant of time passes the cut
        // if false, does not pass cut.
        Bool_t      checkAgainstMask( Int_t eventMJD, Double_t eventTime )
        {
            return checkMaskNow( exactMJD_secs( eventMJD, eventTime ) );
        }
        Bool_t      checkAgainstMask( Double_t eventUTC )
        {
            return checkMaskNow( secs_day * eventUTC );
        }
        
        // Count events passing gamma and direction cuts
        void        countOn( Double_t eventUTC )
        {
            counted.at( Int_t( floor( secs_day * eventUTC - start_time ) ) )++;
        }
        
        
        Bool_t      getMaskStatus() const         // Has the mask encountered any errors; is it still active?
        {
            return !override;
        }
        Int_t       getRunID() const
        {
            return run_id;
        }
        Int_t       getMaskStartMJD() const
        {
            return timeMJD( start_time );
        }
        Double_t    getMaskStartTime() const
        {
            return time_secsMid( start_time );
        }
        Double_t    getMaskStartUTC() const
        {
            return start_time / secs_day;
        }
        UInt_t      getMaskSize() const           // The length of the run according to the mask
        {
            return mask.size();
        }
        vector< Bool_t > getMask() const
        {
            return mask;
        }
        UInt_t      getEffectiveDuration() const; // Number of open seconds in the mask
        UInt_t      getEventsTotal() const;       // Total of all events checked against mask
        UInt_t      getAcceptedTotal() const;     // Total of events allowed to pass the mask
        UInt_t      getCountedTotal() const;      // Total of events passing gamma and direction cuts
        void        resetEventTotals()
        {
            checked.assign( mask.size(), 0 );
            accepted.assign( mask.size(), 0 );
        }
        
        // Various takes on the mean time of the run; all in seconds since the beginning of the run.
        Double_t    getMeanTime_Run() const;      // Mid-point of the run
        Double_t    getMeanTime_Mask() const;     // Weighted by open seconds
        Double_t    getMeanTime_Events() const;   // Weighted by all events
        Double_t    getMeanTime_Accepted() const; // Weighted by events in open seconds
        // Various takes on the mean time of the run; all in <MJD>.<fraction of day>
        Double_t    getMeanUTC_Run() const
        {
            return ( start_time + getMeanTime_Run() ) / secs_day;
        }
        Double_t    getMeanUTC_Mask() const
        {
            return ( start_time + getMeanTime_Mask() ) / secs_day;
        }
        Double_t    getMeanUTC_Events() const
        {
            return ( start_time + getMeanTime_Events() ) / secs_day;
        }
        Double_t    getMeanUTC_Accepted() const
        {
            return ( start_time + getMeanTime_Accepted() ) / secs_day;
        }
        
        //Fill vectors with statistics used in rate plots, over intervals defined as near as possible to a user defined width
        void        getIntervalRates( vector< double >& event_count, vector< double >& interval_time, vector< double >& interval_size, double width ) const;
        
        // ROOT object versions of the internal arrays, to be saved to a ROOT file.
        void        writeObjects() const;
        const   TBits*      getMaskBits() const;  // Copy of the mask as a TBits object.
        // To check this version of the mask use:
        //	UInt_t GetNbits() const 			- for the size
        //	UInt_t CountBits(UInt_t startBit = 0) const	- for the number of open seconds
        //	Bool_t TestBitNumber( UInt_t second ) const 	- for one entry
        // TVector of all events checked against mask.
        const   TVector*    getCheckedVector() const;
        // TVector of events allowed to pass the mask.
        const   TVector*    getAcceptedVector() const;
        // Read stored ROOT versions of a time mask and fill initialised internal vectors.
        Bool_t      readObjects( TDirectory* iDir );
        
        // Display mean times determined according to open portions of the mask and optionally to event counts
        void        printMeanTime( Bool_t event_statistics = kFALSE, ostream& terminal = cout ) const;
        // Display summary of mask with optional event counts
        void        printMask( UInt_t interval_seconds = 60, Bool_t event_statistics = kFALSE, ostream& terminal = cout ) const;
        // display entire mask as a block of 0's and 1's, divided up into lines of 1 minute
        // much better for visually seeing the mask
        void        displayMask( ostream& terminal = cout ) ;
        
        vector<string> splitpath( string& str, set<char> delimiters ) ;
        
        string      getMaskFileName() ;
        
        ClassDef( VTimeMask, 5 ) ;
};
#endif                                            /* ifndef VTIMEMASK_H */
