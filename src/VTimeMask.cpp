/* VTimeMask.cpp
 * Created: Mon  2 Mar 00:11:01 2009
 *
 */

#include "VTimeMask.h"

using namespace std;

VTimeMask::VTimeMask()
{
    run_id                  = -1;
    start_time              = 0;
    mask_file               = "";
    
    mask.reserve( 1200 );
    checked.reserve( 1200 );
    accepted.reserve( 1200 );
    counted.reserve( 1200 );
    
    set_success             = kFALSE;
    override                = kTRUE;
    outside_count                   = 0;
}


VTimeMask::VTimeMask( Int_t run_number, Int_t run_startMJD, Double_t run_startTime, Int_t run_endMJD, Double_t run_endTime, string file_name )
{
    Double_t startMJD_secs  = exactMJD_secs( run_startMJD, run_startTime );
    Double_t endMJD_secs    = exactMJD_secs( run_endMJD, run_endTime );
    
    initialise( run_number, startMJD_secs, endMJD_secs, file_name );
}


VTimeMask::VTimeMask( Int_t run_number, Double_t run_startUTC, Double_t run_endUTC, string file_name )
{
    Double_t startMJD_secs  = secs_day * run_startUTC;
    Double_t endMJD_secs    = secs_day * run_endUTC;
    
    initialise( run_number, startMJD_secs, endMJD_secs, file_name );
}


void            VTimeMask::initialise( Int_t run_number, Double_t run_start, Double_t run_end, string file_name )
{
    run_id                                  = run_number;
    start_time                              = run_start;
    mask_file                               = ""; // Default filename
    if( !file_name.empty() )
    {
        mask_file = file_name;
    }
    set_success                             = kFALSE;
    override                                = kFALSE;
    outside_count                               = 0;
    
    UInt_t duration_seconds = 0;
    if( run_end - run_start > secs_day )
    {
        problem( "run duration greater than one day" );
    }
    else if( run_end < run_start )
    {
        problem( "end time before start time" );
    }
    else
    {
        duration_seconds = UInt_t( ceil( run_end - run_start ) );
    }
    
    mask.assign( duration_seconds, kTRUE );
    checked.assign( duration_seconds, 0 );
    accepted.assign( duration_seconds, 0 );
    counted.assign( duration_seconds, 0 );
    
    cout    << "\t time-mask initialised for run " << run_id
            << " for " << mask.size() << " seconds"
            << " starting at " << getMaskStartTime() << " secs on MJD " << getMaskStartMJD()
            << endl;
}


void            VTimeMask::problem( string reason, ostream& terminal )
{
    override = kTRUE;
    terminal    << "VTimeMask error: improperly initialised (" << reason << ")" << endl
                << "                 Will pass all events and take first event as start time." << endl;
    //terminal	<<"The story so far ..."<< endl;
    //printMask( 0, kTRUE );
}


Bool_t          VTimeMask::setMask( Int_t run_number, Double_t run_startUTC, Double_t run_endUTC, string file_name )
{
    Double_t startMJD_secs  = secs_day * run_startUTC;
    Double_t endMJD_secs    = secs_day * run_endUTC;
    
    initialise( run_number, startMJD_secs, endMJD_secs, file_name );
    
    return setMask();
}


Bool_t          VTimeMask::setMask()
{
    // Variables for user input
    Int_t run_number    = -1;
    Int_t mask_start    = -1;
    Int_t mask_width    = -1;
    Int_t mask_mode     =  0;
    
    // Variables for derived quantities
    Int_t mask_end      = -1;
    
    if( !mask_file.empty() )
    {
        ifstream settings_file( mask_file.c_str() );
        if( !settings_file.is_open() )
        {
            cout << "VTimeMask::setMask warning: failed to find time-mask file named \'" << mask_file << "\' (check your settings in runparameter.dat)." << endl;
            exit( EXIT_FAILURE );
        }
        else
        {
            cout << "\t setting time-mask from file " << mask_file << endl;
            // Retrieve each line and parse
            string settings_str = "";
            while( getline( settings_file, settings_str ) )
            {
                if( settings_str[0] != '*' )
                {
                    continue;
                }
                settings_str.erase( 0, 1 );
                istringstream settings_stream( settings_str );
                settings_stream >> run_number
                                >> mask_start
                                >> mask_width
                                >> mask_mode;
                                
                // Check sanity of user input (oh, those crazy users!)
                // Correct run number
                if( run_number != run_id )
                {
                    continue;
                }
                // If mask definition starts after run ends don't bother loading
                if( mask_start >= Int_t( mask.size() ) )
                {
                    continue;
                }                                 // If mask starts before run, crop the leading part
                else if( mask_start < 0 )
                {
                    mask_width = mask_width + mask_start;
                    mask_start = 0;
                }
                // Mask width should be finite
                if( mask_width < 1 )
                {
                    continue;
                }
                
                // Mask mode should represent a boolean
                if( mask_mode != 0 && mask_mode != 1 )
                {
                    continue;
                }
                
                cout    << "\t -> user requests " << mask_width
                        << " secs " << ( mask_mode == 0 ? "closed" : "open" )
                        << " starting at " << mask_start << " seconds ..." << endl;
                        
                // Derive limits of mask
                mask_end = mask_start + mask_width;
                mask_end = ( mask_end < Int_t( mask.size() ) ? mask_end : mask.size() );
                
                for( Int_t t = mask_start; t < mask_end ; t++ )
                {
                    mask.at( t ) = ( mask_mode == 0 ? kFALSE : kTRUE );
                }
                
                set_success = kTRUE;
                //printMask(100);
            }
        }
    }
    if( set_success )
    {
        cout << "\t -> mask is open for " << getEffectiveDuration() << " secs." << endl;
    }
    else
    {
        cout << "\t -> default mask is open for entire duration of run." << endl;
    }
    
    return set_success;
}


Bool_t VTimeMask::setMaskDuringPhaseCuts( Double_t now )
{

    now = secs_day * now;
    Int_t instant;
    if( abs( now - start_time ) < 0.1 )
    {
        instant = 0;
    }
    else
    {
        instant = Int_t( floor( now - start_time ) );
    }
    
    //  cout << "\t now start_time instant "<< now << " " << start_time << " " << now-start_time << " " << instant << endl;
    mask.at( instant ) = kFALSE;
    //  set_success = kTRUE
    //cout << "\t phase cuts are being applied, mask is open "<< endl;
    
    return kTRUE;
}



Bool_t VTimeMask::checkMaskNow( Double_t now )
{
    if( override )
    {
        return loadMaskNow( now );
    }
    
    Bool_t result = kFALSE;
    Int_t instant =  Int_t( floor( now - start_time ) );
    
    if( instant < 0 || instant >= Int_t( mask.size() ) )
    {
        outside_count++;
    }
    else
    {
        checked.at( instant )++;
        result = mask.at( instant );
        if( result )
        {
            accepted.at( instant )++;
        }
    }
    
    return result;
}


Bool_t          VTimeMask::loadMaskNow( Double_t now )
{
    // Let's start at the beginning
    if( mask.empty() )
    {
        start_time = now;
    }
    
    Int_t instant =  Int_t( floor( now - start_time ) );
    if( instant < 0 )
    {
        return kFALSE;
    }
    
    // Ensure mask can contain this instant
    if( mask.size() == UInt_t( instant ) )
    {
        mask.push_back( kTRUE );
        checked.push_back( 0 );
        accepted.push_back( 0 );
        counted.push_back( 0 );
    }
    else if( mask.size() < UInt_t( instant ) )
    {
        mask.resize( instant + 1, kTRUE );
        checked.resize( instant + 1, 0 );
        accepted.resize( instant + 1, 0 );
        counted.resize( instant + 1, 0 );
    }
    checked.at( instant )++;
    accepted.at( instant )++;
    
    return kTRUE;
}


UInt_t          VTimeMask::getEffectiveDuration() const
{
    UInt_t dur = 0;
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        if( mask.at( t ) )
        {
            dur++;
        }
    }
    return dur;
}


UInt_t          VTimeMask::getEventsTotal() const
{
    UInt_t count = 0;
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        count += checked.at( t );
    }
    return count;
}


UInt_t          VTimeMask::getAcceptedTotal() const
{
    UInt_t count = 0;
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        count += accepted.at( t );
    }
    return count;
}


UInt_t          VTimeMask::getCountedTotal() const
{
    UInt_t count = 0;
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        count += counted.at( t );
    }
    return count;
}


Double_t        VTimeMask::getMeanTime_Run() const
{
    Double_t    count       = Double_t( mask.size() );
    Double_t    sum         = 0;
    
    // Accumulate weighted sum and weighted count
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        sum     += ( t + 0.5 );                   // Weight by center of mask bin
    }
    
    // Calculate weighted mean
    return ( count > 0 ? sum / count : -1. );
}


Double_t        VTimeMask::getMeanTime_Mask() const
{
    Double_t    count       = Double_t( getEffectiveDuration() );
    Double_t    sum         = 0;
    
    // Accumulate weighted sum and weighted count
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        if( mask.at( t ) )
        {
            sum   += ( t + 0.5 );    // Weight by center of mask bin
        }
    }
    
    // Calculate weighted mean
    return ( count > 0 ? sum / count : -1. );
}


Double_t        VTimeMask::getMeanTime_Events() const
{
    Double_t    count       = Double_t( getEventsTotal() );
    Double_t    sum         = 0;
    
    // Accumulate weighted sum and weighted count
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        // Weight by center of mask bin
        sum         += ( t + 0.5 ) * checked.at( t );
    }
    
    // Calculate weighted mean
    return ( count > 0 ? sum / count : -1. );
}


Double_t        VTimeMask::getMeanTime_Accepted() const
{
    Double_t    count       = Double_t( getAcceptedTotal() );
    Double_t    sum         = 0;
    
    // Accumulate weighted sum and weighted count
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        // Weight by center of mask bin
        sum         += ( t + 0.5 ) * accepted.at( t );
    }
    
    // Calculate weighted mean
    return ( count > 0 ? sum / count : -1. );
}


void            VTimeMask::printMeanTime( Bool_t event_statistics, ostream& terminal ) const
{
    Double_t    all_res         = getMeanTime_Run();
    Double_t    open_res        = getMeanTime_Mask();
    Double_t    event_res       = 0.;
    if( event_statistics )
    {
        event_res       = getMeanTime_Events();
    }
    Double_t    accepted_res    = 0.;
    if( event_statistics )
    {
        accepted_res    = getMeanTime_Accepted();
    }
    
    // Accumulate output
    ostringstream   all_str( "" );
    ostringstream   open_str( "" );
    ostringstream   event_str( "" );
    ostringstream   accepted_str( "" );
    
    // Formatting for VStereoAnalysis output
    all_str << "\t ";
    open_str << "\t ";
    event_str << "\t ";
    accepted_str << "\t ";
    
    // Labels
    all_str << "Mean of run:\t\t\t";
    open_str << "Mean of mask:\t\t\t";
    event_str << "Mean of events:\t\t";
    accepted_str << "Mean of accepted events:\t";
    
    if( all_res > 0 )
    {
        all_str.width( 7 );
        all_str << all_res << "\t (" << getMaskStartTime() + all_res << ")";
    }
    else
    {
        all_str << "N/A (zero length)";
    }
    
    if( open_res > 0 )
    {
        open_str.width( 7 );
        open_str << open_res << "\t (" << getMaskStartTime() + open_res << ")";
    }
    else
    {
        open_str << "N/A (zero length)";
    }
    
    if( event_res > 0 )
    {
        event_str.width( 7 );
        event_str << event_res << "\t (" << getMaskStartTime() + event_res << ")";
    }
    else
    {
        event_str << "N/A (zero length)";
    }
    
    if( accepted_res > 0 )
    {
        accepted_str.width( 7 );
        accepted_str << accepted_res << "\t (" << getMaskStartTime() + accepted_res << ")";
    }
    else
    {
        accepted_str << "N/A (zero length)";
    }
    
    //Display on STDOUT or STDERR
    terminal << endl << "\t Mean times from time mask for run (start on MJD):\t " << run_id << "\t (" << getMaskStartTime() << " secs on " << getMaskStartMJD() << ")" << endl;
    terminal                    <<  all_str.str() << endl;
    if( event_statistics )
    {
        terminal        <<  event_str.str() << endl;
    }
    if( set_success )
    {
        terminal             <<  open_str.str() << endl;
    }
    if( event_statistics && set_success )
    {
        terminal     <<  accepted_str.str() << endl;
    }
    
    //terminal<< endl;
}


void            VTimeMask::printMask( UInt_t interval_seconds, Bool_t event_statistics, ostream& terminal ) const
{
    if( interval_seconds == 0 )
    {
        interval_seconds = 100;
    }
    Int_t column_width = event_statistics ? 5 : 4 ;
    
    // Set-up the rows of the table
    ostringstream   interval_str( "" );
    ostringstream   open_str( "" );
    ostringstream   openFrac_str( "" );
    ostringstream   checked_str( "" );
    ostringstream   accepted_str( "" );
    ostringstream   acceptedFrac_str( "" );
    
    openFrac_str.setf( ios_base::fixed, ios_base::floatfield );
    openFrac_str.precision( 2 );
    
    acceptedFrac_str.setf( ios_base::fixed, ios_base::floatfield );
    acceptedFrac_str.precision( 2 );
    
    // Formatting for VStereoAnalysis output
    interval_str << "\t ";
    open_str << "\t ";
    openFrac_str << "\t ";
    checked_str << "\t ";
    accepted_str << "\t ";
    acceptedFrac_str << "\t ";
    
    // Preface with totals for the entire run
    interval_str            << "Run (";
    interval_str.width( 4 );
    interval_str            << mask.size() << " secs)"               << " | ";
    open_str.width( 15 );
    open_str            << getEffectiveDuration()               << " | ";
    openFrac_str.width( 15 );
    if( mask.size() > 0 )
    {
        openFrac_str        << Double_t( getEffectiveDuration() ) / mask.size()   << " | ";
    }
    else
    {
        openFrac_str        << "N/A"                         << " | ";
    }
    checked_str.width( 15 );
    checked_str         << getEventsTotal()                     << " | ";
    accepted_str.width( 15 );
    accepted_str            << getAcceptedTotal()                   << " | ";
    acceptedFrac_str.width( 15 );
    if( getEventsTotal() > 0 )
    {
        acceptedFrac_str    << Double_t( getAcceptedTotal() ) / getEventsTotal()  << " | ";
    }
    else
    {
        acceptedFrac_str    << "N/A"                         << " | ";
    }
    
    //Row titles
    interval_str        << "Interval (" << interval_seconds << " secs)"    << "\t| ";
    open_str        << "Open (secs)\t"               << "\t| ";
    openFrac_str        << "Open Fraction"               << "\t| ";
    checked_str     << "Total Events\t"              << "\t| ";
    accepted_str        << "Accepted Events"             << "\t| ";
    acceptedFrac_str    << "Accepted Fraction"               << "\t| ";
    
    // Load the rows of the time-mask table
    UInt_t interval             = 0;
    UInt_t interval_count       = 0;
    Double_t open_count         = 0;
    Double_t checked_count      = 0;
    Double_t accepted_count     = 0;
    
    for( UInt_t t = 0; t < mask.size(); t++ )
    {
        // Accumulate statistics
        interval_count++;
        if( mask.at( t ) )
        {
            open_count++;
        }
        checked_count += checked.at( t );
        accepted_count += accepted.at( t );
        
        // Dump output for a completed interval and prepare for the next one
        if( interval_count == interval_seconds )
        {
            interval_str.width( column_width );
            interval_str            << interval                 << " | ";
            open_str.width( column_width );
            open_str            << open_count               << " | ";
            openFrac_str.width( column_width );
            openFrac_str            << open_count / interval_seconds      << " | ";
            checked_str.width( column_width );
            checked_str         << checked_count            << " | ";
            accepted_str.width( column_width );
            accepted_str            << accepted_count           << " | ";
            if( checked_count > 0 )
            {
                acceptedFrac_str.width( column_width );
                acceptedFrac_str    << accepted_count / checked_count     << " | ";
            }
            else
            {
                acceptedFrac_str.width( column_width );
                acceptedFrac_str    << " N/A"                << " | ";
            }
            
            interval_count  = 0;
            open_count      = 0;
            checked_count       = 0;
            accepted_count  = 0;
            interval++;
        }
        
    }
    // Dump output for any remaining partial interval
    if( interval_count != 0 )
    {
        interval_str.width( column_width );
        interval_str                << interval << " (" << interval_count << " sec)";
        open_str.width( column_width );
        open_str                << open_count;
        openFrac_str.width( column_width );
        openFrac_str                << open_count / interval_count;
        checked_str.width( column_width );
        checked_str             << checked_count;
        accepted_str.width( column_width );
        accepted_str                << accepted_count;
        if( checked_count > 0 )
        {
            acceptedFrac_str.width( column_width );
            acceptedFrac_str        << accepted_count / checked_count;
        }
        else
        {
            acceptedFrac_str.width( column_width );
            acceptedFrac_str        << " N/A";
        }
    }
    
    //Display on STDOUT or STDERR
    terminal    << endl << "\t Time mask for run: " << run_id << endl;
    
    // VStereoAnalysis takes run start as second event and run end as penultimate event (might be fixed sometime ... ?),
    // so under normal circumstances you can expect up to two events outside the mask domain.
    if( outside_count > 2 ) terminal   << "VTimeMask warning: " << outside_count << ( outside_count == 1 ? " event" : " events" )
                                           << " rejected because " << ( outside_count == 1 ? "it" : "they" )
                                           << " fell outside the time domain of the mask!" << endl;
                                           
    terminal    << interval_str.str()           << endl;
    if( set_success )
    {
        terminal << open_str.str()       << endl;
    }
    if( set_success )
    {
        terminal << openFrac_str.str()       << endl;
    }
    
    if( event_statistics )
    {
        terminal    << checked_str.str()            << endl;
        if( set_success )
        {
            terminal << accepted_str.str()       << endl;
        }
        if( set_success )
        {
            terminal << acceptedFrac_str.str()   << endl;
        }
    }
    
    //terminal	<< endl;
    
    return;
}


void        VTimeMask::getIntervalRates( vector< double >& event_count, vector< double >& interval_time, vector< double >& interval_size, double width ) const
{
    event_count.resize( 0 );
    interval_time.resize( 0 );
    interval_size.resize( 0 );
    
    double total = 0.;
    double mean = 0.;
    double live = 0.;
    
    double count = 0.;
    for( unsigned int t = 0; t < mask.size(); t++ )
    {
        total += counted.at( t );
        mean  += t * mask.at( t );
        live  += mask.at( t );
        
        count += 1.;
        
        if( count > width || t == mask.size() - 1 )
        {
            double default_meanUTC = ( start_time + ( interval_time.size() + 0.5 ) * width ) / secs_day;
            double adjusted_meanUTC = ( live > 0 ? ( start_time + mean / live ) / secs_day : 0 );
            
            event_count.push_back( total );
            interval_time.push_back( live > 0 ? adjusted_meanUTC : default_meanUTC );
            interval_size.push_back( live );
            
            total = 0.;
            mean = 0.;
            live = 0.;
            
            count = 0.;
        }
    }
}


void            VTimeMask::writeObjects() const
{
    TDirectory* iDir = gDirectory;
    TDirectory* wDir = 0;
    
    iDir->cd();
    wDir = ( TDirectory* )iDir->Get( "timeMask" );
    if( !wDir )
    {
        iDir->mkdir( "timeMask" )->cd();
    }
    else
    {
        wDir->cd();
    }
    
    // Mask objects 'newed'
    const   TBits*      iMaskBits   = getMaskBits();
    const   TVector*    iCheckedVector = getCheckedVector();
    const   TVector*    iAcceptedVector = getAcceptedVector();
    
    // Write mask objects
    iMaskBits->Write( "maskBits" );
    iCheckedVector->Write( "checkedEvtsVector" );
    iAcceptedVector->Write( "acceptedEvtsVector" );
    
    // remove all objects created with new in this class
    delete iMaskBits;
    delete iCheckedVector;
    delete iAcceptedVector;
    
    iDir->cd();
}


const TBits*    VTimeMask::getMaskBits() const
{
    TBits* bitsMask = new TBits( mask.size() );
    
    for( UInt_t i = 0; i < mask.size(); i++ )
    {
        bitsMask->SetBitNumber( i, mask.at( i ) );
    }
    
    return bitsMask;
}


const TVector*  VTimeMask::getCheckedVector() const
{
    float array[mask.size()];
    for( UInt_t i = 0; i < mask.size(); i++ )
    {
        array[i] = checked.at( i );
    }
    
    TVector* vectorChecked = new TVector( mask.size(), array );
    return vectorChecked;
}


const TVector*  VTimeMask::getAcceptedVector() const
{
    float array[mask.size()];
    for( UInt_t i = 0; i < mask.size(); i++ )
    {
        array[i] = accepted.at( i );
    }
    
    TVector* vectorAccepted = new TVector( mask.size(), array );
    return vectorAccepted;
}


Bool_t          VTimeMask::readObjects( TDirectory* iDir )
{
    if( !iDir )
    {
        return kFALSE;
    }
    
    cout << "Reading time mask from file " << iDir->GetPath() << endl;
    
    //Fetch objects
    const   TBits*      iMaskBits       = ( TBits* ) iDir->Get( "maskBits" );
    const   TVector*    iCheckedVector  = ( TVector* ) iDir->Get( "checkedEvtsVector" );
    const   TVector*    iAcceptedVector = ( TVector* ) iDir->Get( "acceptedEvtsVector" );
    
    //Parse into internal structures
    UInt_t maskSize = iMaskBits->GetNbits();
    
    mask.resize( maskSize, kTRUE );
    checked.resize( maskSize, 0 );
    accepted.resize( maskSize, 0 );
    
    for( UInt_t i = 0; i < maskSize; i++ )
    {
        mask.at( i )      = iMaskBits->TestBitNumber( i );
        // The following contains ROOT of a graphic nature: viewer discretion is advised.
        checked.at( i )   = UInt_t( ( iCheckedVector->GetMatrixArray() )[i] );
        accepted.at( i )  = UInt_t( ( iAcceptedVector->GetMatrixArray() )[i] );
    }
    
    // Does ROOT 'new' objects when it gets them from a file?
    // remove all objects created with new in this class
    //delete iMaskBits;
    //delete iCheckedVector;
    //delete iAcceptedVector;
    
    return kTRUE;
}


void VTimeMask::displayMask( ostream& terminal )
{
    double currMJD = 0.0 ;
    double startMJD = start_time / secs_day ;
    double stopMJD  = ( start_time + mask.size() ) / secs_day ;
    //terminal << "~~~~~~~~~~~~~~~~~~~" << endl;
    //terminal.precision(8) ;
    //terminal << "startMJD:" << startMJD << endl ;
    //terminal << "stopMJD :" << stopMJD  << endl ;
    // loop through every second of the run
    for( int i_dur = 0 ; i_dur < 10000 ; i_dur++ )
    {
        currMJD = startMJD + ( i_dur / ( 24.*60.*60. ) ) ;
        if( currMJD >= stopMJD )
        {
            break ;
        }
        // every minute since beginning of the run, start a new line
        if( i_dur % 60 == 0 )
        {
            terminal << endl ;
            terminal << " TIMEMASK RUN"    ;
            terminal.precision( 0 ) ;
            terminal << fixed << setw( 5 )  << run_id  ;
            terminal << " - MJD:" ;
            terminal.precision( 5 ) ;
            terminal << fixed << setw( 11 ) << currMJD ;
            terminal << " - sec:" ;
            terminal.precision( 0 ) ;
            terminal << fixed << setw( 4 )  << i_dur   ;
            terminal << " - "     ;
        }
        // display mask status at currMJD
        terminal << checkMaskNow( secs_day * currMJD ) ;
    }
    terminal << endl;
}

/*
 *  split up a path based on a few delimiting characters, into a vector of strings
 */
vector<string> VTimeMask::splitpath( string& str, std::set<char> delimiters )
{
    vector<string> result;
    
    char const* pch = str.c_str();
    char const* start = pch;
    for( ; *pch; ++pch )
    {
        if( delimiters.find( *pch ) != delimiters.end() )
        {
            if( start != pch )
            {
                string str( start, pch );
                result.push_back( str );
            }
            else
            {
                result.push_back( "" );
            }
            start = pch + 1;
        }
    }
    result.push_back( start );
    
    return result;
}

/*
 * get the mask file name
 * 'mask_file' by itself is hardcoded, so here we just get the basename, and splice it into the current $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/
 */
string VTimeMask::getMaskFileName()
{
    string EVNDISPSYS = gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" ) ;
    char delim_chars[] = { '/' } ;
    set<char> delims( delim_chars, delim_chars + 1 ) ;
    //set<char> delims{'/'};
    vector<string> path = splitpath( mask_file, delims );
    char relative_path[300] ;
    sprintf( relative_path, "%s/ParameterFiles/%s", EVNDISPSYS.c_str(), path.back().c_str() ) ;
    string relative_path_str( relative_path ) ;
    cout << "relative_path : " << relative_path_str << endl;
    return relative_path_str ;
    //return mask_file ;
}


/* Modeline for ViM {{{
* vim:set ts=4:
* vim600:fdm=marker fdl=0 fdc=3:
* }}} */
