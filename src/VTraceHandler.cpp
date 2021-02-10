/*! \class VTraceHandler
    \brief calculation of trace specific stuff like charge sums, getting the maximum, etc.


*/

#include "VTraceHandler.h"
#include "TMath.h"

VTraceHandler::VTraceHandler()
{
    fpTrace.assign( 64, 0. );
    fPed = 0.;
    fTraceAverageTime = 0.;
    fChanID = 0;
    fpTrazeSize = 0;
    fHiLo = false;
    fDynamicRange = 216;                 // 8bit FADC, switch to low-gain; used only for debugging info
    fMaxThreshold = 150;                 // used for trace max calculation in low gain
    // TMPTMPTMP
    //	fMaxThreshold = 0;                 // used for trace max calculation in low gain
    // TMPTMPTMP
    fMC_FADCTraceStart = 0;
    fpulsetiming_maxPV = 0;
    fpulsetiminglevels_size = 0;
    fFindPulseTiming = false;
    
    fTraceIntegrationMethod = 1;
    kIPRmeasure  = false;
    
    setDigitalFilterParameters();
    fDF_tracemax = 0.;
}

void VTraceHandler::reset()
{
    fTraceAverageTime = 0.;
    fSumWindowFirst = 0;
    fSumWindowLast  = 0;
    fHiLo = false;
}

/*
    analysis routines call this function

    (generally, this function is called)
*/
void VTraceHandler::setTrace( VVirtualDataReader* iReader, unsigned int iNSamples, double ped, unsigned int iChanID, unsigned int iHitID, double iHiLo )
{
    fPed = ped;
    // TMPTMPTMP
    //        fPed = 0.;
    // TMPTMPTMP
    fChanID = iChanID;
    
    reset();
    
    if( !iReader )
    {
        cout << "VTraceHandler::setTrace( VVirtualDataReader* iReader, double ped, unsigned int iChanID ): no reader set" << endl;
        return;
    }
    
    ///////////////////////////////////////
    // copy trace from raw data reader
    if( iNSamples != fpTrace.size() )
    {
        fpTrace.clear();
        for( unsigned int i = 0; i < iNSamples; i++ )
        {
            fpTrace.push_back( iReader->getSample_double( iHitID, i + fMC_FADCTraceStart, ( i == 0 ) ) );
        }
    }
    else for( unsigned int i = 0; i < iNSamples; i++ )
        {
            fpTrace[i] = iReader->getSample_double( iHitID, i + fMC_FADCTraceStart, ( i == 0 ) );
        }
        
    fpTrazeSize = fpTrace.size();
    
    ////////////////////////////
    // apply hi-lo gain ratio
    fHiLo = apply_lowgain( iHiLo );
    
    ////////////////////////////
    // apply digital filtering
    // (overwrites trace)
    if( fDF_method > 0 )
    {
        apply_digitalFilter();
    }
}

/*
 *  used only for time jitter calibration
 *
 */
void VTraceHandler::setTrace( vector<uint16_t> pTrace, double ped, unsigned int iChanID, double iHiLo )
{
    fPed = ped;
    fChanID = iChanID;
    reset();
    // copy trace
    unsigned int i_tsize = pTrace.size();
    if( i_tsize != fpTrace.size() )
    {
        fpTrace.clear();
        for( unsigned int i = 0; i < i_tsize; i++ )
        {
            fpTrace.push_back( ( double )pTrace[i] );
        }
    }
    else for( unsigned int i = 0; i < i_tsize; i++ )
        {
            fpTrace[i] = ( double )pTrace[i];
        }
        
    fpTrazeSize = fpTrace.size();
    fHiLo = apply_lowgain( iHiLo );
}

/*
 *  used only for time jitter calibration
 *
 */
void VTraceHandler::setTrace( vector<uint8_t> pTrace, double ped, unsigned int iChanID, double iHiLo )
{
    fPed = ped;
    fChanID = iChanID;
    reset();
    // copy trace
    unsigned int i_tsize = pTrace.size();
    if( i_tsize != fpTrace.size() )
    {
        fpTrace.clear();
        for( unsigned int i = 0; i < i_tsize; i++ )
        {
            fpTrace.push_back( ( double )pTrace[i] );
        }
    }
    else
    {
        for( unsigned int i = 0; i < i_tsize; i++ )
        {
            fpTrace[i] = ( double )pTrace[i];
        }
    }
    
    fpTrazeSize = fpTrace.size();
    fHiLo = apply_lowgain( iHiLo );
}


bool VTraceHandler::apply_lowgain( double iHiLo )
{
    // hilo switch is set
    if( iHiLo > 0. )
    {
        ///// TMPTMPTMP
        //                iHiLo = 1.;
        ///// TMPTMPTMP
        for( unsigned int i = 0; i < fpTrazeSize; i++ )
        {
            fpTrace[i]  = ( fpTrace[i] - fPed ) * iHiLo;
            fpTrace[i] += fPed;
        }
        return true;
    }
    return false;
}

/*
 * sum up FADC trace from fFirst to fLast
 *
 */
double VTraceHandler::calculateTraceSum_fixedWindow( unsigned int fFirst, unsigned int fLast, bool fRaw )
{
    double sum = 0.;
    double tcharge = 0.;
    fTraceAverageTime = 0.;
    for( unsigned int i = fFirst; i < fLast; i++ )
    {
        // require that trace is >0.
        // (CTA MC write trace values above a certain signal only)
        if( i < fpTrace.size() && fpTrace[i] > 0. )
        {
            if( !fRaw )
            {
                sum += fpTrace[i] - fPed;
                tcharge += ( i + 0.5 ) * ( fpTrace[i] - fPed );
            }
            else
            {
                sum += fpTrace[i];
                tcharge += ( i + 0.5 ) * fpTrace[i];
            }
        }
    }
    if( TMath::IsNaN( sum ) )
    {
        sum = 0.;
    }
    if( TMath::Abs( sum ) < 1.e-10 )
    {
        sum = 0.;
        fTraceAverageTime = 0.;
    }
    else
    {
        fTraceAverageTime = tcharge / sum;
    }
    return sum;
}

/*
 *  determines timing if L2 pulse to calculate crate timing jitter
 *
 *
*/
vector<float> VTraceHandler::getFADCTiming( unsigned int fFirst, unsigned int fLast, bool debug )
{

    if( fLast - fFirst <= 20 ) // small readout window -> don't bother with extra step
    {
        return getPulseTiming( fFirst, fLast, fFirst, fLast );
    }
    unsigned int i_start = fFirst;
    unsigned int i_stop = fLast;
    double trace_max = 0.;
    unsigned int n255 = 0;
    unsigned int maxpos = 0;
    getTraceMax( fFirst, fLast, trace_max, maxpos, n255 );
    
    bool have_first = false;
    bool have_second = false;
    
    float temp = 0;
    
    //find first bin above 50 dc & first bin after that where the trace goes down again
    for( unsigned int i = fFirst; i < fLast && !have_second; i++ )
    {
        if( !have_first && ( fpTrace[i] - fPed ) > 40 )
        {
            i_start = i;
            have_first = true;
        }
        if( have_first && ( fpTrace[i] - fPed ) < temp )
        {
            i_stop = i;
            have_second = true;
        }
        temp = fpTrace[i] - fPed;
    }
    
    if( !have_first && debug )
    {
        cout << "VTraceHandler::getFADCTiming()  Warning: coulnd't find bin with signal > 40 dc in range " << fFirst << " - " << fLast << endl;
    }
    if( i_start >= 4 )
    {
        i_start -= 4;
    }
    while( i_start < fFirst )
    {
        i_start++;
    }
    i_stop += 4;
    while( i_stop > fLast )
    {
        if( i_stop > 0 )
        {
            i_stop--;
        }
    }
    return getPulseTiming( i_start, i_stop, i_start, i_stop );
}

/*!
   calculate pulse timing
   (pulse times at different fraction of maximum)

   fFirst, fLast:   range where maximum is determined
   fTFirst, fTLast: range where timing parameters are determined
   iReverseSearchinLowGain: search time parameters from right to left
   (important some VERITAS runs with incomplete timing calibration)

*/

vector< float >& VTraceHandler::getPulseTiming(
           unsigned int fFirst, unsigned int fLast,
           unsigned int fTFirst, unsigned int fTLast,
           bool iReverseSearchinLowGain )
{
    // reset pulse timing vector
    for( unsigned int i = 0; i < fpulsetiming.size(); i++ )
    {
        fpulsetiming[i] = 0.;
    }
    unsigned int m_pos = 0;
    
    // by definition are there always an odd number of values -> centre value is 1
    double i_trace = 0.;
    
    // get pulse maximum
    double trace_max = 0.;
    unsigned int n255 = 0;
    unsigned int maxpos = 0;
    getTraceMax( fFirst, fLast, trace_max, maxpos, n255, iReverseSearchinLowGain );
    if( maxpos == 99999 )
    {
        return fpulsetiming;
    }
    fpulsetiming[fpulsetiming_maxPV] = ( float )maxpos + 0.5;
    
    // first half of the pulse
    // (loop backwards over pulse)
    bool bBreak = false;
    fFindPulseTiming = false;
    if( maxpos >= fpTrazeSize - 1 )
    {
        maxpos--;
    }
    for( unsigned int i = maxpos; i > fTFirst ; i-- )
    {
        i_trace = fpTrace[i] - fPed;
        // loop over all pulse level
        for( unsigned int m = 0; m < fpulsetiming_maxPV; m++ )
        {
            if( fpulsetiming_maxPV >= 1 + m )
            {
                m_pos = fpulsetiming_maxPV - 1 - m;
                if( m_pos < fpulsetiminglevels_size && fpulsetiming[m_pos] < 1.e-5 )
                {
                    if( i_trace < fpulsetiminglevels[m_pos] * trace_max )
                    {
                        fpulsetiming[m_pos] = getLinInterpol( fpulsetiminglevels[m_pos] * trace_max,
                                                              i, i_trace, i + 1, fpTrace[i + 1] - fPed );
                        if( m_pos == 0 )
                        {
                            fFindPulseTiming = true;
                            bBreak = true;
                        }
                    }
                }
            }
        }
        if( bBreak )
        {
            break;
        }
    }
    // second half of the pulse: contains pulse widths
    // (loop forwards over pulse)
    bBreak = false;
    if( maxpos > 0 )
    {
        for( unsigned int i = maxpos; i < fTLast; i++ )
        {
            i_trace = fpTrace[i] - fPed;
            // loop over all pulse level
            for( unsigned int m_pos = fpulsetiming_maxPV + 1; m_pos < fpulsetiminglevels_size; m_pos++ )
            {
                if( m_pos < fpulsetiminglevels_size && fpulsetiming[m_pos] < 1.e-5 )
                {
                    if( i_trace < fpulsetiminglevels[m_pos] * trace_max )
                    {
                        fpulsetiming[m_pos] = getLinInterpol( fpulsetiminglevels[m_pos] * trace_max, i, i_trace, i - 1, fpTrace[i - 1] - fPed );
                        fpulsetiming[m_pos] -= fpulsetiming[fpulsetiminglevels_size - m_pos - 1];
                        if( m_pos == fpulsetiminglevels_size - 1 )
                        {
                            bBreak = true;
                        }
                    }
                }
            }
            if( bBreak )
            {
                break;
            }
        }
    }
    
    return fpulsetiming;
}


/*
 *
 *   return trace value and position of trace maximum
 *
 *   trace value is pedestal subtracted
 *
 */
void VTraceHandler::getTraceMax(
        unsigned int fFirst, unsigned int fLast,
        double& tmax, unsigned int& maxpos,
        unsigned int& n255,
        bool iReverseSearchinLowGain )
{
    unsigned int nMax = ( unsigned int )( fDynamicRange * tmax );
    n255 = 0;
    // value at maximum
    tmax = -10000.;
    // position of maximum (in samples)
    maxpos = 99999;
    /////////////////////////////////////////////////////
    // determine maximum position in a high gain channel
    if( !fHiLo || !iReverseSearchinLowGain )
    {
        if( fFirst < fLast && fLast <= fpTrazeSize )
        {
            for( unsigned int i = fFirst; i < fLast; i++ )
            {
                if( fpTrace[i] > tmax )
                {
                    tmax   = fpTrace[i];
                    maxpos = i;
                }
            }
            tmax -= fPed;
        }
    }
    /////////////////////////////////////////////////////
    // low gain channel
    // (needs special treatment as end of the saturated high gain pulse is
    //  occassionally at the beginning of the readout window)
    else
    {
        if( fFirst < fLast && fLast <= fpTrazeSize )
        {
            // start search at end of the integration window
            for( unsigned int i = fLast; i > fFirst; i-- )
            {
                if( fpTrace[i - 1] < fMaxThreshold && tmax < 0. )
                {
                    continue;
                }
                if( fpTrace[i - 1] < fMaxThreshold / 2 )
                {
                    break;
                }
                if( fpTrace[i - 1] > tmax )
                {
                    tmax   = fpTrace[i - 1];
                    maxpos = i - 1;
                }
                // do some rough counting of saturation
                if( nMax > 0 && fpTrace[i - 1] > nMax )
                {
                    n255++;
                }
            }
            tmax -= fPed;
        }
    }
}


double VTraceHandler::getTraceMax()
{
    double tmax = 0.;
    unsigned int maxposInt = 0;
    unsigned int n255 = 0;
    getTraceMax( 0, fpTrazeSize, tmax, maxposInt, n255 );
    return tmax;
}

double VTraceHandler::getTraceMax( unsigned int& n255, unsigned int& maxposInt )
{
    double tmax = 0.;
    getTraceMax( 0, fpTrazeSize, tmax, maxposInt, n255 );
    return tmax;
}

double VTraceHandler::getLinInterpol( double y5, int x1, double y1, int x2, double y2 )
{
    double a = 0.;
    double b = 0.;
    if( x2 - x1 != 0 )
    {
        a = ( y2 - y1 ) / ( double )( x2 - x1 );
    }
    // shift by 0.5 to locate bin center
    if( a != 0. && x2 + 0.5 > 0 )
    {
        b = y2 - a * ( ( double )( x2 ) + 0.5 );
    }
    
    if( a != 0. )
    {
        return ( y5 - b ) / a;
    }
    
    return 0.;
}

void VTraceHandler::setPulseTimingLevels( vector< float > iP )
{
    fpulsetiminglevels = iP;
    fpulsetiminglevels_size = fpulsetiminglevels.size();;
    fpulsetiming.assign( fpulsetiminglevels_size, 0. );
    fpulsetiming_maxPV = ( fpulsetiminglevels_size - 1 ) / 2;
}

/*

    select trace integration method

    1: get trace sum between given integration range

    2: get maximum sum from sliding window method (use integration range to calculate integration window only)

*/
bool VTraceHandler::setTraceIntegrationmethod( unsigned int iT )
{
    // check method numbers
    if( iT > 5 )
    {
        return false;
    }
    
    fTraceIntegrationMethod = iT;
    
    return true;
}

/*
 *   trace integration
 *
 *   note that there are a number of different trace integration methods
 *
 */
double VTraceHandler::getTraceSum( unsigned int iSumWindowFirst, unsigned int iSumWindowLast, bool iRaw,
                                   unsigned int iTraceIntegrationMethod, bool iForceWindowStart,
                                   unsigned int iSlidingWindowLast )
{
    // set trace integration method
    if( iTraceIntegrationMethod < 9999 )
    {
        fTraceIntegrationMethod = iTraceIntegrationMethod;
    }
    
    // integrate from fFirst to fLast
    if( fTraceIntegrationMethod == 1 )
    {
        fSumWindowFirst = iSumWindowFirst;
        fSumWindowLast  = iSumWindowLast;
        return calculateTraceSum_fixedWindow( fSumWindowFirst, fSumWindowLast, iRaw );
    }
    // find maximum integral
    else if( fTraceIntegrationMethod == 2 )
    {
        if( !kIPRmeasure )
        {
            // special case: search over restricted window
            if( iForceWindowStart )
            {
                if( iSlidingWindowLast >= fpTrace.size() )
                {
                    iSlidingWindowLast = fpTrace.size();
                }
                return calculateTraceSum_slidingWindow( iSumWindowFirst,
                                                        iSlidingWindowLast,
                                                        iSumWindowLast - iSumWindowFirst,
                                                        iRaw );
            }
            // default: search over whole summation window
            else
            {
                return calculateTraceSum_slidingWindow( 0, fpTrace.size(), iSumWindowLast - iSumWindowFirst, iRaw );
            }
        }
        // IPR measurements from long trace file
        else
        {
            return calculateTraceSum_slidingWindow( 0.5 * fpTrace.size(),
                                                    fpTrace.size(),
                                                    iSumWindowLast - iSumWindowFirst,
                                                    iRaw );
        }
    }
    // return simple the trace maximum as trace sum
    else if( fTraceIntegrationMethod == 5 )
    {
        double peakamplitude = getTraceMax();
        double result = 0;
        if( iRaw )
        {
            result = peakamplitude + fPed;
        }
        else
        {
            result = peakamplitude;
        }
        return result;
    }
    // return digital filter results
    else if( fTraceIntegrationMethod == 6 )
    {
        return fDF_tracemax;
    }
    
    return 0.;
}


/*
 *
 * get maximum trace sum
 * (sliding window, search along trace for maximum sum)
 *
*/
double VTraceHandler::calculateTraceSum_slidingWindow( unsigned int iSearchStart,
        unsigned int iSearchEnd,
        int iIntegrationWindow,
        bool fRaw )
{
    unsigned int n = fpTrace.size();
    unsigned int window = iIntegrationWindow;
    unsigned int SearchEnd = iSearchEnd;
    if( ( n - window ) <= SearchEnd )
    {
        SearchEnd = n - window + 1;
    }
    int lolimit = 0;
    int uplimit = 0;
    float tcharge = 0.;
    double ampl = 0.;
    double charge = 0.;
    float ped = fPed;
    if( kIPRmeasure )
    {
        ped = 0.;
    }
    ////////////////////////////////////////
    // sample time and value (ped subracted)
    float muxBINS[n], FADC[n];
    for( unsigned int i = 1; i <= n; i++ )
    {
        muxBINS[i - 1] = i - 0.5;
        FADC[i - 1] = ( float )fpTrace.at( i - 1 ) - ped;
    }
    
    if( n == 0 )
    {
        fTraceAverageTime = muxBINS[1];
        return 0.;
    }
    
    ////////////////////////////////////////
    // special case for ped calculation
    if( fRaw )
    {
        for( unsigned int i = 0; i < ( unsigned int )iIntegrationWindow; i++ )
        {
            charge += ( float )fpTrace.at( n - 1 - i );
        }
        fTraceAverageTime = muxBINS[n - 1];
        fSumWindowFirst = n - iIntegrationWindow;
        fSumWindowLast  = n;
        
        return FADC[1];
    }
    
    ////////////////////////////////////////
    float xmax = 0.;
    for( unsigned int i = iSearchStart; i < int( window ) + iSearchStart; i++ )
    {
        xmax += FADC[i];
    }
    // extract charge
    for( unsigned int i = iSearchStart; i < SearchEnd; i++ )
    {
        if( charge < xmax )
        {
            charge = xmax;
            uplimit = i + int( window );
            if( uplimit > int( n ) )
            {
                uplimit = n;
            }
            lolimit = i;
            if( lolimit < 0 )
            {
                lolimit = 0;
            }
        }
        xmax = xmax - FADC[i] + FADC[i + window];
    }
    // arrival times *****************************************************
    for( int k = lolimit; k < uplimit; k++ )
    {
        tcharge += muxBINS[k] * FADC[k];
    }
    
    if( charge != 0. )
    {
        fTraceAverageTime = tcharge / charge;
    }
    if( fTraceAverageTime < iSearchStart )
    {
        fTraceAverageTime = 0.;
    }
    if( fTraceAverageTime > ( ( int )SearchEnd + ( int )window - 1 ) )
    {
        fTraceAverageTime = ( ( int )SearchEnd + ( int )window - 1 );
    }
    ampl -= fPed;
    
    fSumWindowFirst = lolimit;
    fSumWindowLast  = uplimit;
    
    return charge;
}

double VTraceHandler::getMaxSumAutoWindow( float AmplThresh, unsigned int iSearchStart, unsigned int iSearchEnd, unsigned int iIntegrationWindow, bool fRaw )
{
    //fRaw = false;
    unsigned int n = fpTrace.size();
    int saturflag = 0;
    int max = 0;
    unsigned int window = iIntegrationWindow;
    unsigned int windowcalib = 2 * iIntegrationWindow;
    unsigned int SearchEnd = iSearchEnd;
    if( ( n - window ) <= SearchEnd )
    {
        SearchEnd = n - window;
    }
    int lolimit = 0;
    int lolimitcal = 0;
    int uplimit = 0;
    int uplimitcal = 0;
    float tcharge = 0, tcharge2 = 0;
    float arrtime = 0, arrtime2 = 0;
    double charge = 0, charge2 = 0, ampl = 0.;
    
    float muxBINS[n], FADC[n];
    for( unsigned int i = 1; i <= n; i++ )
    {
        muxBINS[i - 1] = i - 0.5;
        FADC[i - 1] = ( float )fpTrace.at( i - 1 ) - fPed;
    }
    
    charge = 0, charge2 = 0;
    tcharge = 0, tcharge2 = 0;
    //maxbin   = LocMax(ampl);
    
    if( n == 0 )
    {
        return -1;
    }
    
    if( fRaw ) // special case for ped calculation
    {
        for( unsigned int i = 0; i < iIntegrationWindow; i++ )
        {
            //charge+=( float )fpTrace.at( i );
            charge += ( float )fpTrace.at( n - 1 - i );
        }
        //charge = FADC[1];
        arrtime = muxBINS[n - 1];
        fTraceAverageTime = arrtime;
        fSumWindowFirst = n - iIntegrationWindow;
        fSumWindowLast  = n;
        
        return charge;
    }
    
    float xmax = 0;
    for( unsigned int i = iSearchStart; i < int( window ) + iSearchStart; i++ )
    {
        xmax += FADC[i];
    }
    // extract charge for small window **********************************
    for( unsigned int i = iSearchStart; i < SearchEnd; i++ )
    {
        if( charge < xmax )
        {
            charge = xmax;
            max = i;
            uplimit = max + int( window );
            if( uplimit > int( n ) )
            {
                uplimit = n - 1;
            }
            lolimit = max;
            if( lolimit < 0 )
            {
                lolimit = 1;
            }
        }
        xmax = xmax - FADC[i] + FADC[i + window];
    }
    // extract charge for big window ************************************
    uplimitcal = int( uplimit ) + ( int( windowcalib ) - int( window ) ) / 2;
    if( uplimitcal > int( SearchEnd ) + ( int )window - 2 )
    {
        uplimitcal = ( int )SearchEnd + ( int )window - 1;
    }
    lolimitcal = int( lolimit ) - ( int( windowcalib ) - int( window ) ) / 2;
    if( lolimitcal < int( iSearchStart ) )
    {
        lolimitcal = ( int )iSearchStart;
    }
    // arrival times *****************************************************
    for( int k = lolimit; k < uplimit; k++ )
    {
        tcharge += muxBINS[k] * FADC[k];
    }
    for( int k = lolimitcal; k < uplimitcal; k++ )
    {
        tcharge2 += muxBINS[k] * FADC[k];
        charge2 += FADC[k];
    }
    // extract saturated charge (integrate everything to the end) ************************************
    if( saturflag )
    {
        tcharge = 0;
        charge = 0;
        charge2 = 0;
        uplimitcal = ( int )SearchEnd + ( int )window - 1;
        lolimitcal = lolimit - ( int( windowcalib ) - int( window ) ) / 2;
        if( lolimitcal < int( iSearchStart ) )
        {
            lolimitcal = iSearchStart;
        }
        for( int k = lolimitcal; k < uplimitcal; k++ )
        {
            tcharge += muxBINS[k] * FADC[k];
            charge2 += FADC[k];
        }
    }
    
    if( charge != 0 )
    {
        arrtime = tcharge / charge;
    }
    if( charge2 != 0 )
    {
        arrtime2 = tcharge2 / charge2;
    }
    if( arrtime < iSearchStart )
    {
        arrtime = 0;
    }
    if( arrtime > ( ( int )SearchEnd + ( int )window - 1 ) )
    {
        arrtime = ( ( int )SearchEnd + ( int )window - 1 );
    }
    if( arrtime2 < iSearchStart )
    {
        arrtime2 = 0;
    }
    if( arrtime2 > ( ( int )SearchEnd + ( int )window - 1 ) )
    {
        arrtime2 = ( ( int )SearchEnd + ( int )window - 1 );
    }
    ampl -= fPed;
    
    // extract charge for automatic integration window (extends window )
    int intwin = 0;
    int Start = lolimit;
    int Stop = uplimit;
    
    float chargeAutoWin = 0.;
    
    for( int i = lolimit; i > int( iSearchStart ); i-- )
    {
        if( FADC[i] >= AmplThresh )
        {
            intwin++;
            Start = i;
        }
        else
        {
            break;
        }
    }
    for( int i = uplimit; i < int( SearchEnd ); i++ )
    {
        if( FADC[i] >= AmplThresh )
        {
            intwin++;
            Stop = i;
        }
        else
        {
            break;
        }
    }
    
    for( int k = Start; k < Stop; k++ )
    {
        chargeAutoWin += FADC[k];
    }
    
    fTraceAverageTime = arrtime;
    fSumWindowFirst = Start;
    fSumWindowLast  = Stop;
    
    return charge;
}

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

bool VTraceHandler::apply_digitalFilter()
{
    // not opimized for speed, but try to use the FlashCam routines
    // in its orginal versions
    
    // input
    int n  = ( int )fpTrace.size();
    float* ip = new float[fpTrace.size()];
    int us = fDF_upsample;
    float bl = fPed;
    float pz = fDF_polezero;
    for( unsigned int t = 0; t < fpTrace.size(); t++ )
    {
        ip[t] = fpTrace[t];
    }
    
    // results
    float* op = new float[fpTrace.size()*fDF_upsample];
    float max = 0.;
    int   at  = 0;
    
    // filter
    PzpsaSmoothUpsampleFloat( n, us, ip, bl, pz, op, &max, &at );
    
    fpTrace.clear();
    for( int i = 0; i < n * us; i++ )
    {
        fpTrace.push_back( op[i] + fPed );
    }
    fpTrazeSize = fpTrace.size();
    
    fDF_tracemax    = max;
    fSumWindowFirst = ( unsigned int )at;
    fSumWindowLast  = fSumWindowFirst + 1;
    
    delete[] ip;
    delete[] op;
    
    return true;
}

/*
 *  digital filter code obtained by FlashCam
 *  (per email JiHi and GeHa)
 *
 */


/*=== Function ===============================================*/

int VTraceHandler::PzpsaSmoothUpsampleFloat
(
    int n, int us, float* ip, float bl, float pz,
    float* op, float* max, int* at
)

/*--- Description --------------------------------------------//

Upsample (expand the n input values to us samples each)
Subtract baseline bl and correct for a single pole
decay with the decay time pz and smooth the resulting trace
with two moving averages with a width of us.
The output is placed in array op and returns the new number
of samples (n*us).

//------------------------------------------------------------*/
{
    int   i, i1; 		// running indexes
    float v2, v1;		// the next and prev. input samples
    float sum1, sum2;	// the running sum of 1.st and 2.nd average
    float tmp;		// a temp var for intermediate copy
    float pzc2, pzc1;	// the next and prev. pz corrected value
    float* out1 = op;	// the out pointer of the first runsum
    float* out2 = op;	// the out pointer of the second runsum
    float mult = 1. / us / us;	// the multiplier to correct the two runsums
    float peakmax = -1e30;	// peak maximum
    int   peakat = 0; 	// peak position
    
    v1 = v2 = ( ip[0] - bl ) * mult;
    pzc2 = pzc1 = v2;
    sum1 = pzc2 * us;
    sum2 = sum1 * us;
    for( i = 0; i < us; i++ )
    {
        *out1++ = sum1;
    }
    for( i = 1; i < n; i++ )
    {
        v2 = ( ip[i] - bl ) * mult;
        pzc2 = ( v2 - v1 );
        v1 = v2;
        for( i1 = 0; i1 < us; i1++ )
        {
            sum1 += pzc2 - pzc1 * pz;
            *out1++ = sum1;
            tmp = *out2;
            *out2++ = sum2;
            if( sum2 > peakmax )
            {
                peakmax = sum2, peakat = out2 - op - 1;
            }
            sum2 += sum1 - tmp;
        }
        pzc1 = pzc2;
    }
    n *= us;
    for( v2 = op[n - 1]; out2 < ( op + n ); )
    {
        tmp = *out2;
        *out2++ = sum2;
        if( sum2 > peakmax )
        {
            peakmax = sum2, peakat = out2 - op - 1;
        }
        sum2 += v2 - tmp;
    }
    if( max )
    {
        *max = peakmax;
    }
    if( at )
    {
        *at = peakat;
    }
    return n;
}

