/** CTA.convert_hessio_to_VDST
 *  short a program to convert sim_telarray files (hessio) to EVNDISP DST format
 *
 *
 *  Author of skeleton (as part of the hessio distribution):  Konrad Bernloehr
 *
 *  Author of modifications for eventdisplay: Gernot Maier (DESY)
 */


#include "initial.h"
#include "io_basic.h"
#include "history.h"
#include "io_hess.h"
#include "fileopen.h"
#ifdef CTA_PROD2_TRGMASK
#include "io_trgmask.h"
#endif

#include <bitset>
#include <iostream>
#include <map>
#include <set>
#include <stdlib.h>
#include <string>
#include <vector>

#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TStopwatch.h>

#include "VGlobalRunParameter.h"
#include "VEvndispRunParameter.h"
#include "VDSTTree.h"
#include "VMonteCarloRunHeader.h"
#include "VAstronometry.h"
///////////////////////////////////////////////////////
// global variables
///////////////////////////////////////////////////////
// number of telescopes (checked)
unsigned int fNTel_checked = 0;
// maximum number of pixels for the current array configuration
unsigned int fGlobalMaxNumberofPixels = 0;
// maximum number of FADC samples for the current array configuration
unsigned int fGlobalMaxNumberofSamples = 0;
// maximum number of telescopes in previous event
unsigned int fGlobalNTelPreviousEvent = 0;
// set if this event does not pass the minimum number of telescope condition
bool fGlobalTriggerReset = false;
// map with telescope types (length of map is total number of telescopes)
map< int, ULong64_t > fTelescopeType;
// set with telescope types (length of map is total number of telescope types)
set< ULong64_t > fTelescopeTypeList;
// map wit simtel telescope IDs (ID, counter)
map< unsigned int, unsigned int > fTelescopeSimTelList;
#ifdef CTA_PROD2_TRGMASK
// trigger mask hash set
struct trgmask_hash_set* fTriggerMask_hash_set = 0;
#endif
// minimum number of photo electrons needed for trigger type 4
int fMinNumber_pe_per_telescope = 30;
// fill leaf with photoelectrons
bool fFillPELeaf = false;
// fill leaf with peak ADC values
bool fFillPeakADC = false;
// additional pedestal shift (e.g. for ASTRI analysis)
float fPedestalShift = 0.;
// Minimum energy cut
float fMinEnergy = 0.;
// Maximum energy cut
float fMaxEnergy = 1.e20;
// use this for writing of telconfig tree only
bool fWriteTelConfigTreeOnly = false;
///////////////////////////////////////////////////////

/*!
      hard wired timing configuration for pulse characterization
      (time of 100%, 20%, 50%.. of max pulse height)

      Changes must also be made to getTimingLevelIndex()

      1     	1
      0.2	2
      0.5	2
      0.8	2
      0.5	4
      0.2	4
      -40	5
*/
vector< float > getPulseTimingLevels()
{
    vector< float > t;
    t.push_back( 0.2 );
    t.push_back( 0.5 );
    t.push_back( 1.0 );
    t.push_back( 0.5 );
    t.push_back( 0.2 );
    
    return t;
}


unsigned int getTimingLevelIndex( unsigned int i )
{
    if( i == 0 )
    {
        return 2;    // 100%
    }
    else if( i == 1 )
    {
        return 0;    //  20%
    }
    else if( i == 2 )
    {
        return 1;    //  50%
    }
    else if( i == 4 )
    {
        return 3;    //  50%
    }
    else if( i == 5 )
    {
        return 4;    //  20%
    }
    
    return 9999;
}

/*
 * read plate scale stretching factor from external file
 *
 * given per telescope type
 *
 * usually a file called $CTA_EVNDISP_AUX_DIR//DetectorGeometry/...EffectiveFocalLength*
 *
 */
map< ULong64_t, float > readCameraScalingMap( string iFile )
{
    map< ULong64_t, float > iCameraScalingMap;
    
    ifstream is;
    is.open( iFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "readCameraScalingMap error: file not found, " << iFile << endl;
        return iCameraScalingMap;
    }
    string iLine;
    string iT1;
    cout << "reading camera scaling configuration from " << iFile << endl;
    
    ULong64_t teltype = 0;
    float scaling = 0.;
    
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 && iLine.substr( 0, 1 ) == "*" )
        {
            istringstream is_stream( iLine );
            is_stream >> iT1;
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> iT1;
                teltype = atoi( iT1.c_str() );
                if( !(is_stream>>std::ws).eof() )
                {
                    is_stream >> iT1;
                    scaling = atof( iT1.c_str() );
                    iCameraScalingMap[teltype] = scaling;
                    cout << "\t camera scaling for telescope type " << teltype << " :\t";
                    cout << 1. + iCameraScalingMap[teltype] << endl;
                }
            }
        }
    }
    is.close();
    
    return iCameraScalingMap;
}


/** The factor needed to transform from mean p.e. units to units of the single-p.e. peak:
    Depends on the collection efficiency, the asymmetry of the single p.e. amplitude
    distribution and the electronic noise added to the signals. */
#define CALIB_SCALE 0.92

/*
 *  FLAG_AMP_TMP       0: Use normal integrated amplitude.
 *                     1: Use integration around global peak position from
 *                        pulse shape analysis. May include all pixels or only selected.
 *                     2: Use integration around local peak position from
 *                        pulse shape analysis. Return 0 for pixels without
 *                        a fairly significant peak.
 */
#define FLAG_AMP_TMP 0

/* ---------------------- calibrate_pixel_amplitude ----------------------- */

/** Calibrate a single pixel amplitude, for cameras with two gains per pixel.
 *
 * @return Pixel amplitude in peak p.e. units.
 */

float calibrate_pixel_amplitude( AllHessData* hsdata, int itel, int ipix, int dummy, unsigned int& iLowGain );

float calibrate_pixel_amplitude( AllHessData* hsdata, int itel, int ipix, int dummy, unsigned int& iLowGain )
{
    int i = ipix, npix, significant, hg_known, lg_known;
    double npe, sig_hg, npe_hg, sig_lg, npe_lg;
    AdcData* raw;
    
    if( hsdata == NULL || itel < 0 || itel >= H_MAX_TEL )
    {
        return 0.;
    }
    npix = hsdata->camera_set[itel].num_pixels;
    if( ipix < 0 || ipix >= npix )
    {
        return 0.;
    }
    raw = hsdata->event.teldata[itel].raw;
    if( raw == NULL )
    {
        return 0.;
    }
    if( ! raw->known )
    {
        return 0.;
    }
    
    significant = hsdata->event.teldata[itel].raw->significant[i];
    
    hg_known = hsdata->event.teldata[itel].raw->adc_known[HI_GAIN][i];
    sig_hg = hg_known ? ( hsdata->event.teldata[itel].raw->adc_sum[HI_GAIN][i] -
                          hsdata->tel_moni[itel].pedestal[HI_GAIN][i] ) : 0.;
    npe_hg = sig_hg * hsdata->tel_lascal[itel].calib[HI_GAIN][i];
    
    lg_known = hsdata->event.teldata[itel].raw->adc_known[LO_GAIN][i];
    sig_lg = lg_known ? ( hsdata->event.teldata[itel].raw->adc_sum[LO_GAIN][i] -
                          hsdata->tel_moni[itel].pedestal[LO_GAIN][i] ) : 0.;
    npe_lg = sig_lg * hsdata->tel_lascal[itel].calib[LO_GAIN][i];
    
    iLowGain = false;
    if( !significant )
    {
        npe = 0.;
    }
    else if( hg_known && sig_hg < 1000000 && sig_hg > -1000 )
    {
        npe = npe_hg;
    }
    else
    {
        npe = npe_lg;
        iLowGain = true;
    }
    
    /* npe is in units of 'mean photo-electrons'. */
    /* We convert to experimentalist's */
    /* 'peak photo-electrons' now. */
    return ( float )( TMath::Nint( CALIB_SCALE * npe ) * 100. ) / 100.;
}

#include <signal.h>

void stop_signal_function( int isig );

static int interrupted;

/* ---------------------- stop_signal_function -------------------- */
/**
 *  Stop the program gracefully when it catches an INT or TERM signal.
 *
 *  @param isig  Signal number.
 *
 *  @return (none)
 */

void stop_signal_function( int isig )
{
    if( isig >= 0 )
    {
        fprintf( stderr, "Received signal %d\n", isig );
    }
    if( !interrupted )
        fprintf( stderr,
                 "Program stop signal. Stand by until current data block is finished.\n" );
                 
    interrupted = 1;
    
    signal( SIGINT, SIG_DFL );
    signal( SIGTERM, SIG_DFL );
}

/** Show program syntax */

static void syntax( char* program );

static void syntax( char* program )
{
    printf( "Syntax: %s [ options ] [ - | input_fname ... ]\n", program );
    printf( "Options:\n" );
    printf( "   -v               (More verbose output)\n" );
    printf( "   -q               (Much more quiet output)\n" );
    printf( "   -s               (Show data explained)\n" );
    printf( "   -S               (Show data explained, including raw data)\n" );
    printf( "   --history (-h)   (Show contents of history data block)\n" );
    printf( "   -i               (Ignore unknown data block types)\n" );
    printf( "   --max-events n   (Skip remaining data after so many triggered events.)\n" );
    printf( "   -a subarray file (list of telescopes to read with FOV.)\n" );
    printf( "   -o dst filename  (name of dst output file)\n" );
    printf( "   -f 0(QADC), 1(FADC), 2(MIXED) (write FADC samples to DST file;default=2)\n" );
    printf( "   -c pedfile.root  (file with pedestals and pedestal variances)\n" );
    printf( "   -t <trgmask directory> (directory with trigger mask files (corrections for Spring 2013 prod2 production)\n" );
    printf( "   -r on=1/off=0    (apply camera plate scaling for DC telescopes; default=1)\n" );
    printf( "   -rfile camera scaling file (read camera place scales from a file; default=none)\n" );
    printf( "   -NSB <float>     (scale pedvars by sqrt(...) of this value\n" );
    printf( "   -pe              (fill leaf with photoelectrons into DST tree (default: off)\n" );
    printf( "   -peakadc         (fill leaf with peakadc values into DST tree (default: off)\n" );
    printf( "   -pedshift float(dc) (apply additional pedestal shift to sum values (default: 0.; good value for ASTRI: 150.)\n" );
    printf( "   -minenergy float(TeV) (apply a minimum energy cut in TeV on events in sim_telarray file.)\n" );
    printf( "   -maxenergy float(TeV) (apply a maximum energy cut in TeV on events in sim_telarray file.)\n" );
    printf( "   -pedshift float(dc) (apply additional pedestal shift to sum values (default: 0.; good value for ASTRI: 150.)\n" );
    
    exit( EXIT_SUCCESS );
}

using namespace std;

/*

    read trigger mask from on external file

    this is the correction for the Spring 2013 prod2 production with wrong trigger settings

*/
bool read_trigger_mask( string trg_mask_file )
{
#ifdef CTA_PROD2_TRGMASK
    struct trgmask_set* tms = ( trgmask_set* )calloc( 1, sizeof( struct trgmask_set ) );
    if( fTriggerMask_hash_set )
    {
        free( fTriggerMask_hash_set );
    }
    fTriggerMask_hash_set = ( trgmask_hash_set* )calloc( 1, sizeof( struct trgmask_hash_set ) );
    
    IO_BUFFER* iobuf = allocate_io_buffer( 1000000L );
    if( iobuf == NULL )
    {
        cout << "read_trigger_mask(): error, cannot allocate I/O buffer" << endl;
        exit( EXIT_FAILURE );
    }
    iobuf->max_length = 200000000L;
    
    IO_ITEM_HEADER item_header;
    iobuf->input_file = fileopen( trg_mask_file.c_str(), "r" );
    if( iobuf->input_file != NULL )
    {
        cout << endl << "reading trigger masks from " << trg_mask_file << endl;
        unsigned int z = 0;
        for( ;; )
        {
            if( find_io_block( iobuf, &item_header ) != 0 )
            {
                break;
            }
            
            printf( "Found I/O block of type %ld\n", item_header.type );
            
            if( read_io_block( iobuf, &item_header ) != 0 )
            {
                break;
            }
            
            read_trgmask( iobuf, tms );
            
            trgmask_fill_hashed( tms, fTriggerMask_hash_set );
            
            z++;
        }
        if( z > 1 )
        {
            cout << "read_trigger_mask(): error, more than one iobuf - code cannot handle this yet" << endl;
            exit( EXIT_FAILURE );
        }
        fileclose( iobuf->input_file );
    }
    else
    {
        cout << "read_trigger_mask(): error, cannot open trigger mask file: " << endl;
        cout << "\t" << trg_mask_file << endl;
        return false;
    }
    
#endif
    return true;
}

/*
 * fill MC run headr (mainly with CORSIKA parameters, e.g. the scatter radius for the core positions)
 *
*/
bool DST_fillMCRunheader( VMonteCarloRunHeader* f, AllHessData* hsdata, bool iAddToFirstFile )
{
    // values from CORSIKA (from first file)
    if( !iAddToFirstFile )
    {
        f->runnumber = hsdata->run_header.run;
        f->shower_prog_id = hsdata->mc_run_header.shower_prog_id;
        f->shower_prog_vers = hsdata->mc_run_header.shower_prog_vers;
        f->detector_prog_id = hsdata->mc_run_header.detector_prog_id;
        f->detector_prog_vers = ( unsigned int )hsdata->mc_run_header.detector_prog_vers;
        f->obsheight = hsdata->mc_run_header.obsheight;
    }
    // TODO: add simulation dates from shower and detector simulation
    // (need a compiler flag for backwards compatibility)
    if( iAddToFirstFile )
    {
        f->num_showers += hsdata->mc_run_header.num_showers;
    }
    else
    {
        f->num_showers = hsdata->mc_run_header.num_showers;
    }
    f->num_use = hsdata->mc_run_header.num_use;
    f->core_pos_mode = hsdata->mc_run_header.core_pos_mode;
    f->core_range[0] = hsdata->mc_run_header.core_range[0];
    f->core_range[1] = hsdata->mc_run_header.core_range[1];
    f->az_range[0] = hsdata->mc_run_header.az_range[0];
    f->az_range[1] = hsdata->mc_run_header.az_range[1];
    f->alt_range[0] = hsdata->mc_run_header.alt_range[0];
    f->alt_range[1] = hsdata->mc_run_header.alt_range[1];
    f->diffuse = hsdata->mc_run_header.diffuse;
    f->viewcone[0] = hsdata->mc_run_header.viewcone[0];
    f->viewcone[1] = hsdata->mc_run_header.viewcone[1];
    if( iAddToFirstFile && f->E_range[0] != hsdata->mc_run_header.E_range[0] )
    {
        cout << "DST_fillMCRunheader(): warning, energy range differ, " << f->E_range[0] << ", " << hsdata->mc_run_header.E_range[0] << endl;
    }
    f->E_range[0] = hsdata->mc_run_header.E_range[0];
    if( iAddToFirstFile && f->E_range[1] != hsdata->mc_run_header.E_range[1] )
    {
        cout << "DST_fillMCRunheader(): warning, energy range differ, " << f->E_range[1] << ", " << hsdata->mc_run_header.E_range[1] << endl;
    }
    f->E_range[1] = hsdata->mc_run_header.E_range[1];
    if( iAddToFirstFile && f->spectral_index != hsdata->mc_run_header.spectral_index )
    {
        cout << "DST_fillMCRunheader(): warning, energy range differ, " << f->spectral_index << ", " << hsdata->mc_run_header.spectral_index << endl;
    }
    f->spectral_index = hsdata->mc_run_header.spectral_index;
    f->B_total = hsdata->mc_run_header.B_total;
    f->B_inclination = hsdata->mc_run_header.B_inclination;
    f->B_declination = hsdata->mc_run_header.B_declination;
    f->injection_height = hsdata->mc_run_header.injection_height;
    f->fixed_int_depth = hsdata->mc_run_header.fixed_int_depth;
    f->atmosphere = hsdata->mc_run_header.atmosphere;
    f->corsika_iact_options = hsdata->mc_run_header.corsika_iact_options;
    f->corsika_low_E_model = hsdata->mc_run_header.corsika_low_E_model;
    f->corsika_high_E_model = hsdata->mc_run_header.corsika_high_E_model;
    f->corsika_bunchsize = hsdata->mc_run_header.corsika_bunchsize;
    f->corsika_wlen_min = hsdata->mc_run_header.corsika_wlen_min;
    f->corsika_wlen_max = hsdata->mc_run_header.corsika_wlen_max;
    f->corsika_low_E_detail = hsdata->mc_run_header.corsika_low_E_detail;
    f->corsika_high_E_detail = hsdata->mc_run_header.corsika_high_E_detail;
    
    if( iAddToFirstFile )
    {
        f->combined_runHeader = true;
    }
    
    return true;
}

/*
 * fill MC event data (e.g. primary type, MC energy, etc)
 *
 * transform into VERITAS coordinate system
 */
bool DST_fillMCEvent( VDSTTree* fData, AllHessData* hsdata )
{
    // MC
    fData->fDSTeventnumber = hsdata->mc_event.event;
    fData->fDSTrunnumber = hsdata->run_header.run;
    fData->fDSTprimary = hsdata->mc_shower.primary_id;
    fData->fDSTenergy = hsdata->mc_shower.energy;
    fData->fDSTaz = hsdata->mc_shower.azimuth * TMath::RadToDeg();
    fData->fDSTze = 90. - hsdata->mc_shower.altitude * TMath::RadToDeg();
    // Observe: transform to VERITAS coordinate system
    // x: East, y: North
    fData->fDSTxcore = -1.* hsdata->mc_event.ycore;
    fData->fDSTycore = hsdata->mc_event.xcore;
    /////////////////////////////////////////////////////////////////////////////
    // calculate offset in camera coordinates from telescope and MC direction
    // (Note: assume that all telescope point into the same direction)
    double i_tel_el = hsdata->run_header.direction[1] * TMath::RadToDeg();
    double i_tel_az = hsdata->run_header.direction[0] * TMath::RadToDeg();
    double j_x, j_y = 0.;
    int j_j = 0;
    VAstronometry::vlaDs2tp( fData->fDSTaz *TMath::DegToRad(), (90.-fData->fDSTze)*TMath::DegToRad(),
                             i_tel_az * TMath::DegToRad(), i_tel_el * TMath::DegToRad(), 
                             &j_x, &j_y, &j_j );
    fData->fDSTTel_xoff = j_x * TMath::RadToDeg();
    fData->fDSTTel_yoff = -1.*j_y * TMath::RadToDeg();
    
    if( fData->fMCtree )
    {
        fData->fMCtree->Fill();
    }
    
    return true;
}

/*

    main function to read and convert data for each event

*/
bool DST_fillEvent( VDSTTree* fData, AllHessData* hsdata, map< unsigned int, VDSTTelescopeConfiguration > telescope_list,
                    unsigned int iWriteFADC )
{
    if( !fData || !hsdata )
    {
        return false;
    }
    
    // Use only high energy events
    if( hsdata->mc_shower.energy < fMinEnergy )
    {
        return true;
    }
    // Use only low energy events
    if( hsdata->mc_shower.energy > fMaxEnergy )
    {
        return true;
    }
    
    // ntel doesn't change from event to event, but needed by VDST
    fData->fDSTntel = ( unsigned short int )hsdata->run_header.ntel;
    if( telescope_list.size() == 0 )
    {
        for( unsigned int i = 0; i < fData->fDSTntel; i++ )
        {
            telescope_list[i + 1].FOV = 0.;
        }
    }
    fData->fDSTntel = telescope_list.size();
    
    // reset all arrays
    fData->resetDataVectors( 99999, fData->fDSTntel, fGlobalNTelPreviousEvent,
                             fGlobalMaxNumberofPixels,
                             getPulseTimingLevels().size(),
                             fGlobalMaxNumberofSamples,
                             fGlobalTriggerReset, true );
                             
    /////////////////////////////////////////////////////////////////////////////////////
    // event data
    fData->fDSTeventnumber = hsdata->mc_event.event;
    fData->fDSTrunnumber = hsdata->run_header.run;
    fData->fDSTeventtype = 0;
    fData->fDSTgpsyear = 2010;
    
    ////////////////////////////////////////////////
    // trigger data
    fData->fDSTNTrig = hsdata->event.central.num_teltrg;
    unsigned int i_ntel_trig = 0;
    bitset<8 * sizeof( unsigned long ) > i_localTrigger;
    // loop over all triggered events
    for( unsigned int t = 0; t < ( unsigned int )hsdata->event.central.num_teltrg; t++ )
    {
        if( telescope_list.find( hsdata->event.central.teltrg_list[t] ) != telescope_list.end() )
        {
            if( hsdata->event.central.teltrg_list[t] < ( int )i_localTrigger.size() )
            {
                i_localTrigger.set( hsdata->event.central.teltrg_list[t] - 1, true );
            }
            if( t < ( unsigned int )hsdata->run_header.ntel )
            {
                fData->fDSTLTrig_list[i_ntel_trig] = hsdata->event.central.teltrg_list[t];
                fData->fDSTLTtime[i_ntel_trig] = hsdata->event.central.teltrg_time[t];
                // read L2 trigger type
                // (filled for >= prod2 only;
                //  although not filled in prod3)
                // PROD2: bit1: majority
                //        bit2: analog sum
                //        bit3: digital sum
                //        bit4: pe trigger (set with fMinNumber_pe_per_telescope)
#if defined(CTA_PROD2) || defined(CTA_PROD3) || defined(CTA_MAX_SC) || defined(CTA_PROD3_MERGE)
                fData->fDSTL2TrigType[i_ntel_trig] = ( unsigned short int )hsdata->event.central.teltrg_type_mask[t];
#else
                fData->fDSTL2TrigType[i_ntel_trig] = 99;
#endif
                // trigger corrects (from log files - needed for Spring 2013 prod2 files)
#ifdef CTA_PROD2_TRGMASK
                if( fTriggerMask_hash_set )
                {
                    struct trgmask_entry* h_te = find_trgmask( fTriggerMask_hash_set, hsdata->mc_event.event,
                                                 hsdata->event.central.teltrg_list[t] );
                    if( h_te )
                    {
                        fData->fDSTL2TrigType[i_ntel_trig] = ( unsigned short int )h_te->trg_mask;
                    }
                    else
                    {
                        cout << "No trigger mask event found: event " << hsdata->mc_event.event;
                        cout << ", telescope " << hsdata->event.central.teltrg_list[t] << endl;
                        fData->fDSTL2TrigType[i_ntel_trig] = 99;
                    }
                }
#endif
                // add an additional artificial trigger flag on minimum number of pe
                if( hsdata->mc_event.mc_pe_list[t].npe > fMinNumber_pe_per_telescope && fData->fDSTL2TrigType[i_ntel_trig] != 99 )
                {
                    fData->fDSTL2TrigType[i_ntel_trig] += 8;
                }
            }
            i_ntel_trig++;
        }
    }
    fData->fDSTLTrig = i_localTrigger.to_ulong();
    fData->fDSTNTrig = i_ntel_trig;
    // check array trigger condition again
    if( fData->fDSTNTrig < ( unsigned int )hsdata->run_header.min_tel_trig )
    {
        // set trigger reset variable (needed for efficient resetting of DST arrays)
        fGlobalTriggerReset = true;
        return true;
    }
    fGlobalTriggerReset = false;
    
    /////////////////////////////////////////////////////////////////////////////////////
    // tracking data
    unsigned int z = 0;
    for( unsigned short int i = 0; i < ( unsigned short int )hsdata->run_header.ntel; i++ )
    {
        if( telescope_list.find( hsdata->event.trackdata[i].tel_id ) != telescope_list.end() )
        {
            fData->fDSTpointAzimuth[z]       = hsdata->event.trackdata[i].azimuth_raw * TMath::RadToDeg();
            fData->fDSTpointElevation[z]     = hsdata->event.trackdata[i].altitude_raw * TMath::RadToDeg();
            fData->fDSTpointTrackingKnown[z] = hsdata->event.trackdata[i].raw_known;
            z++;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // event data
    // (main loop)
    fData->fDSTntel_data = ( unsigned short int )hsdata->event.central.num_teldata;
    unsigned short int i_ntel_data = 0;
    
    // loop over all telescopes
    for( unsigned short int i = 0; i < fData->fDSTntel_data; i++ )
    {
        // check that this telescope is defined in the current array
        if( telescope_list.find( hsdata->event.central.teldata_list[i] ) == telescope_list.end() )
        {
            continue;
        }
        // get telescope ID
        fData->fDSTtel_data[i_ntel_data] = ( unsigned int )hsdata->event.central.teldata_list[i];
        unsigned int telID = fData->fDSTtel_data[i_ntel_data] - 1;
        
        if( fTelescopeSimTelList.find( telID + 1 ) != fTelescopeSimTelList.end() )
        {
            telID = fTelescopeSimTelList[telID + 1];
        }
        else
        {
            cout << "Error in DST_fillEvent: unknown telescope " << endl;
            cout << "\t telescope ID: " << telID + 1 << endl;
            cout << "\t telescope # in list: " << i << endl;
            cout << "\t total number of telescopes: " << hsdata->event.num_tel << endl;
            cout << "\t number of telescopes in this list: " << hsdata->event.central.num_teldata << endl;
            exit( EXIT_FAILURE );
        }
        // newer hessio versions need this statement for whatever reason
#if CTA_SC>1
        if( hsdata->event.teldata[telID].known == 0 )
        {
            continue;
        }
#endif
        ////////////////////////////////////////////////
        // get pixel data
        fData->fDSTTelescopeZeroSupression[i_ntel_data] = ( unsigned short int )hsdata->event.teldata[telID].raw->zero_sup_mode;
        fData->fDSTnumSamples[i_ntel_data]              = ( unsigned short int )hsdata->event.teldata[telID].raw->num_samples;
        
        // test of FADC samples are requested and then, if there are samples
        if( iWriteFADC == 1 && fData->fDSTnumSamples[i_ntel_data] == 0 )
        {
            cout << "Error in DST_fillEvent: sample length in hessio file is null for telescope ";
            cout << telID << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        
        // set maximum number of FADC samples (needed for efficient resetting of DST arrays)
        if( fData->fDSTnumSamples[i_ntel_data] > fGlobalMaxNumberofSamples )
        {
            fGlobalMaxNumberofSamples = fData->fDSTnumSamples[i_ntel_data];
        }
        
        // check number of pixel (change in inc/VGlobalRunParameters if not sufficient
        if( hsdata->camera_set[telID].num_pixels >= VDST_MAXCHANNELS )
        {
            cout << "ERROR: number of pixels too high: telescope " << telID + 1 << ": is: ";
            cout << hsdata->camera_set[telID].num_pixels << "\t max allowed: " << VDST_MAXCHANNELS << endl;
            exit( EXIT_FAILURE );
        }
        
        // get dynamic range from list of telescopes (for low gain)
        int i_lowGainSwitch = -1;
        if( telescope_list[fData->fDSTtel_data[i_ntel_data]].DynamicRange > 0 )
        {
            i_lowGainSwitch = ( int )( telescope_list[fData->fDSTtel_data[i_ntel_data]].DynamicRange * 0.99 );
        }
        
        /////////////////////////////////////
        // set maximum number of pixels (needed for efficient resetting of DST arrays)
        if( ( unsigned int )hsdata->camera_set[telID].num_pixels > fGlobalMaxNumberofPixels )
        {
            fGlobalMaxNumberofPixels = ( unsigned int )hsdata->camera_set[telID].num_pixels;
        }
        
        /////////////////////////////////////
        // loop over all pixel
        for( int p = 0; p < hsdata->camera_set[telID].num_pixels; p++ )
        {
            // reset low gain switch (filled later)
            unsigned int iLowGain = false;
            fData->fDSTHiLo[i_ntel_data][p] = iLowGain;
            fData->fDSTChan[i_ntel_data][p] = ( unsigned int )p;
            
            // channels not known are set dead
            fData->fDSTdead[i_ntel_data][p] = !( hsdata->event.teldata[telID].raw->adc_known[HI_GAIN][p] );
            
            ///////////////////////////////////
            // fill pe counting
            if( fFillPELeaf )
            {
                fData->fDSTPe[i_ntel_data][p] = hsdata->mc_event.mc_pe_list[telID].pe_count[p];
            }
            if( fFillPeakADC )
            {
                fData->fDSTPadcHG[i_ntel_data][p] = ( unsigned short int )( hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p] );
                fData->fDSTPadcLG[i_ntel_data][p] = ( unsigned short int )( hsdata->event.teldata[telID].raw->adc_sum[LO_GAIN][p] );
            }
            ////////////////////////////////////////////
            // zero suppression: channel recorded? Bit 0: sum, 1: samples
            // temporary: not clear how that works
            fData->fDSTZeroSuppressed[i_ntel_data][p] = 0;
            /*                        fData->fDSTZeroSuppressed[i_ntel_data][p] = hsdata->event.teldata[telID].raw->significant[p];
                                    // no suppression:     0
                                    // charge suppressed:  1
                                    // samples suppressed: 2
                                    // (tmp) in the following, either charge or samples are filled
                                    if( hsdata->event.teldata[telID].raw->significant[p] == 1 )
                                    {
                                        fData->fDSTZeroSuppressed[i_ntel_data][p] = 2;  // samples suppressed (charge available)
                                    }
                                    else if( hsdata->event.teldata[telID].raw->significant[p] > 1 )
                                    {
                                        fData->fDSTZeroSuppressed[i_ntel_data][p] = 1;  // charge suppressed (samples available)
                                    }
                                    else
                                    {
                                        fData->fDSTZeroSuppressed[i_ntel_data][p] = 3;  // charge and samples suppressed
                                    } */
            
            ///////////////////////
            // fill FADC trace
            unsigned int iTraceIsZero = 1;
            if( iWriteFADC )   // allowed values are 0/1/2
            {
                ///////////////////////////////////////
                // fill FADC trace
                iTraceIsZero = 0;
                iLowGain = HI_GAIN;
                // first bit = samples are available
                //                                if ( ((hsdata->event.teldata[telID].raw->adc_known[HI_GAIN][p]) & (1 << (1)))
                if( hsdata->event.teldata[telID].raw->adc_sample && hsdata->event.teldata[telID].raw->adc_sample[HI_GAIN]
                        && hsdata->event.teldata[telID].raw->num_samples > 0 )
                    // TMPTMP ZEROSUP                                 && hsdata->event.teldata[telID].raw->significant[p] > 1 )  // channel recorded? Bit 0: sum, 1: samples
                {
                    // LOW GAIN
                    // check if this trace is in low gain
                    // (simply check if HI_GAIN is above certain value)
                    if( i_lowGainSwitch > 0
                            && hsdata->event.teldata[telID].raw->num_gains == 2
                            && hsdata->event.teldata[telID].raw->adc_known[LO_GAIN][p] != 0 
                            && hsdata->event.teldata[telID].raw->adc_sample && hsdata->event.teldata[telID].raw->adc_sample[LO_GAIN] )
                    {
                        for( int t = 0; t < hsdata->event.teldata[telID].raw->num_samples; t++ )
                        {
                            // check if HI_GAIN trace value is above threshold, if yes, flip low gain
                            if( ( int )hsdata->event.teldata[telID].raw->adc_sample[HI_GAIN][p][t] > ( int )( i_lowGainSwitch ) )
                            {
                                iLowGain = LO_GAIN;
                            }
                        }
                    }
                    fData->fDSTHiLo[i_ntel_data][p] = ( bool )( iLowGain == LO_GAIN );
                    
                    // fill FADC trace into dst tree
                    for( int t = 0; t < hsdata->event.teldata[telID].raw->num_samples; t++ )
                    {
                        iTraceIsZero += hsdata->event.teldata[telID].raw->adc_sample[iLowGain][p][t];
                        fData->fDSTtrace[i_ntel_data][t][p] = hsdata->event.teldata[telID].raw->adc_sample[iLowGain][p][t];
                    }
                }
                /////////////////////////////
                // no FADC trace available
                else if( telescope_list.find( fData->fDSTtel_data[i_ntel_data] ) != telescope_list.end() )
                {
                    if( telescope_list[fData->fDSTtel_data[i_ntel_data]].RAWsum )
                        //                                       && hsdata->event.teldata[telID].raw->significant[p] > 0 )
                    {
                        fData->fDSTsums[i_ntel_data][p] = hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p]
                                                          - hsdata->tel_moni[telID].pedestal[HI_GAIN][p]
                                                          - fPedestalShift;
                    }
                    else
                    {
                        fData->fDSTsums[i_ntel_data][p] = calibrate_pixel_amplitude( hsdata, telID, p, FLAG_AMP_TMP, iLowGain );
                    }
                    fData->fDSTHiLo[i_ntel_data][p] = false;
                    
                    //  use raw sum for low-gain determination
                    if( i_lowGainSwitch > 0
                            && hsdata->event.teldata[telID].raw->num_gains == 2 )
                    {
                        if( ( int )hsdata->event.teldata[telID].raw->adc_sum[HI_GAIN][p] > i_lowGainSwitch )
                        {
                            if( telescope_list[fData->fDSTtel_data[i_ntel_data]].RAWsum )
                            {
                                fData->fDSTsums[i_ntel_data][p] = hsdata->event.teldata[telID].raw->adc_sum[LO_GAIN][p]
                                                                  - hsdata->tel_moni[telID].pedestal[LO_GAIN][p];
                            }
                            else
                            {
                                fData->fDSTsums[i_ntel_data][p] = calibrate_pixel_amplitude( hsdata, telID, p, FLAG_AMP_TMP, iLowGain );
                            }
                            // apply hi/lo multiplier
                            if( hsdata->tel_lascal[telID].calib[HI_GAIN][p] > 0. )
                            {
                                fData->fDSTsums[i_ntel_data][p] *= hsdata->tel_lascal[telID].calib[LO_GAIN][p] /
                                                                   hsdata->tel_lascal[telID].calib[HI_GAIN][p];
                            }
                            fData->fDSTHiLo[i_ntel_data][p] = true;
                        }
                    }
                }
                // pedestal
                if( hsdata->tel_moni[telID].num_ped_slices > 0 )
                {
                    fData->fDSTpedestal[i_ntel_data][p] = hsdata->tel_moni[telID].pedestal[iLowGain][p]
                                                          / ( double )( hsdata->tel_moni[telID].num_ped_slices );
                }
                // fill pulse timing level as set in sim_telarray: note, there are some hardwired levels here!
                for( int t = 0; t < hsdata->event.teldata[telID].pixtm->num_types; t++ )
                {
                    if( t < ( int )fData->getDSTpulsetiminglevelsN()
                            && getTimingLevelIndex( t ) < fData->getDSTpulsetiminglevelsN() )
                    {
                        fData->fDSTpulsetiminglevels[i_ntel_data][getTimingLevelIndex( t )] = hsdata->event.teldata[telID].pixtm->time_level[t];
                    }
                }
                ////////////////////////////////////////////////
                // mean pulse timing and pulse timing levels
                for( int t = 0; t < hsdata->event.teldata[telID].pixtm->num_types; t++ )
                {
                    // sim_telarray timing levels
                    if( t < ( int )fData->getDSTpulsetiminglevelsN() && getTimingLevelIndex( t ) < fData->getDSTpulsetiminglevelsN() )
                    {
                        fData->fDSTpulsetiming[i_ntel_data][getTimingLevelIndex( t )][p] = hsdata->event.teldata[telID].pixtm->timval[p][t];
                        if( getTimingLevelIndex( t ) == 1 && fData->fDSTsums[i_ntel_data][p] > fData->getDSTMeanPulseTimingMinLightLevel()
                                && hsdata->event.teldata[telID].pixtm->timval[p][t] > 1.e-1 )
                        {
                            fData->fillDSTMeanPulseTiming( telID, p, hsdata->event.teldata[telID].pixtm->peak_global,
                                                           hsdata->event.teldata[telID].raw->num_samples );
                        }
                    }
                }
                                // fill characteristic trigger time in pulse timing level vector after all entries before:
                                if( getPulseTimingLevels().size() < fData->getDSTpulsetiminglevelsN() )
                                {
                                        unsigned int iTriggerTimeIndex = getPulseTimingLevels().size();
                                        fData->fDSTpulsetiminglevels[i_ntel_data][iTriggerTimeIndex] = -1.;
                                        // read timing from trigger timing (for ASTRI - prod4)
                                        if( hsdata->event.teldata[telID].pixeltrg_time.known 
                                            && hsdata->event.teldata[telID].pixeltrg_time.num_times > 0 )
                                        {
                                              bool b_PixTimeFound = false;
                                              // current pixel is p - search for it in list
                                              for( int ip = 0; ip < hsdata->event.teldata[telID].pixeltrg_time.num_times; ip++ )
                                              {
                                                   if( hsdata->event.teldata[telID].pixeltrg_time.pixel_list[ip] == p )
                                                   {
                                                        fData->fDSTpulsetiming[i_ntel_data][iTriggerTimeIndex][p] = hsdata->event.teldata[telID].pixeltrg_time.pixel_time[ip];
                                                        fData->fDSTpulsetiming[i_ntel_data][iTriggerTimeIndex][p] *= hsdata->event.teldata[telID].pixeltrg_time.time_step;
                                                        b_PixTimeFound = true;
                                                        break;
                                                   }
                                              }
                                              if( !b_PixTimeFound )
                                              {
                                                   fData->fDSTpulsetiming[i_ntel_data][iTriggerTimeIndex][p] = -999.;
                                              }
                                        } 
                                  }

            }
            /////////////////////////////
            // fill QADC results
            else
            {
                fData->fDSTRecord[i_ntel_data][p] = hsdata->event.teldata[telID].raw->significant[p];
                // (low gain is not treated correctly here)
                fData->fDSTsums[i_ntel_data][p] = calibrate_pixel_amplitude( hsdata, telID, p, FLAG_AMP_TMP, iLowGain );
                fData->fDSTMax[i_ntel_data][p] = ( short )( hsdata->event.teldata[telID].pixtm->pulse_sum_loc[HI_GAIN][p] );
                if( FLAG_AMP_TMP > 0 )
                {
                    fData->fDSTsumwindow[i_ntel_data][p] = hsdata->event.teldata[telID].pixtm->after_peak;
                    fData->fDSTsumwindow[i_ntel_data][p] -= hsdata->event.teldata[telID].pixtm->before_peak;
                }
                else
                {
                    fData->fDSTsumwindow[i_ntel_data][p] = hsdata->pixel_set[telID].sum_bins;
                }
                // fill timing information
                for( int t = 0; t < hsdata->event.teldata[telID].pixtm->num_types; t++ )
                {
                    if( t < ( int )fData->getDSTpulsetiminglevelsN()
                            && getTimingLevelIndex( t ) < fData->getDSTpulsetiminglevelsN() )
                    {
                        fData->fDSTpulsetiming[i_ntel_data][getTimingLevelIndex( t )][p] = hsdata->event.teldata[telID].pixtm->timval[p][t];
                    }
                }
            }
        }
        //////////////////////////////
        // get local trigger data (this is not corrected for clipping!)
        unsigned int i_nL1trig = 0;
        for( int p = 0; p < hsdata->event.teldata[telID].trigger_pixels.pixels; p++ )
        {
            if( hsdata->event.teldata[telID].trigger_pixels.pixel_list[p] < VDST_MAXCHANNELS )
            {
                fData->fDSTL1trig[i_ntel_data][hsdata->event.teldata[telID].trigger_pixels.pixel_list[p]] = 1;
                i_nL1trig++;
            }
        }
        fData->fDSTnL1trig[i_ntel_data] = i_nL1trig;
        i_ntel_data++;  // telescope counter
    }
    fGlobalNTelPreviousEvent = i_ntel_data;
    fData->fDSTntel_data = i_ntel_data;
    
    ///////////////////
    // MC event data
    fData->fDSTprimary = hsdata->mc_shower.primary_id;
    fData->fDSTenergy = hsdata->mc_shower.energy;
    fData->fDSTaz = hsdata->mc_shower.azimuth * TMath::RadToDeg();
    fData->fDSTze = 90. - hsdata->mc_shower.altitude * TMath::RadToDeg();
    // Observe: transform to VERITAS coordinate system
    // x: East, y: North
    fData->fDSTxcore = -1.*hsdata->mc_event.ycore;
    fData->fDSTycore =     hsdata->mc_event.xcore;
    /////////////////////////////////////////////////////////////////////////////
    // calculate offset in camera coordinates from telescope and MC direction
    // (OBS: this assumes that all telescopes are pointing into the same direction)
    double i_tel_el = 0.;
    double i_tel_az = 0.;
    for( int itel = 0; itel < hsdata->run_header.ntel; itel++ )
    {
        if( fData->fDSTpointTrackingKnown[itel] > 0 )
        {
            i_tel_el = fData->fDSTpointElevation[itel];
            i_tel_az = fData->fDSTpointAzimuth[itel];
            break;
        }
    }
    double j_x, j_y = 0.;
    int j_j = 0;
    VAstronometry::vlaDs2tp( fData->fDSTaz *TMath::DegToRad(), (90.-fData->fDSTze)*TMath::DegToRad(),
                             i_tel_az * TMath::DegToRad(), i_tel_el * TMath::DegToRad(), 
                             &j_x, &j_y, &j_j );
    fData->fDSTTel_xoff = j_x * TMath::RadToDeg();;
    fData->fDSTTel_yoff = -1.*j_y * TMath::RadToDeg();;
    /////////////////////////////////////////////////////////////////////////////
    
    if( fData->fDST_tree )
    {
        fData->fDST_tree->Fill();
    }
    
    return true;
}

/*

   fill calibration tree

       e.g. pedestals, tzeros, gains, etc.

*/
TList* DST_fillCalibrationTree( VDSTTree* fData, AllHessData* hsdata,
                                map< unsigned int, VDSTTelescopeConfiguration> telescope_list,
                                string ipedfile, float iNSBScaling )
{
    if( !hsdata || fWriteTelConfigTreeOnly )
    {
        return 0;
    }
    cout << "filling calibration tree ";
    if( ipedfile.size() > 0 )
    {
        cout << " (using data from external ped file: " << ipedfile << ")";
    }
    else
    {
        cout << " (using data from simtelarray file)";
    }
    cout << endl;
    if( iNSBScaling <= 0. )
    {
        cout << "DST_fillCalibrationTree: invalid NSB scaling factor (should be >0): iNSBScaling" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    int fTelID = 0;
    
    // save current directory
    TDirectory* i_curDir = gDirectory;
    
    // number of pixels
    unsigned int nPixel = 0;
    // integration window
    unsigned int fnum_sumwindow = 1;
    unsigned int fsumwindow[VDST_MAXSUMWINDOW];
    float fPed_high[VDST_MAXCHANNELS];
    float* fPedvar_high = new float[VDST_MAXCHANNELS * VDST_MAXSUMWINDOW];
    float fPed_low[VDST_MAXCHANNELS];
    float* fPedvar_low = new float[VDST_MAXCHANNELS * VDST_MAXSUMWINDOW];
    float fConv_high[VDST_MAXCHANNELS];
    float fConv_low[VDST_MAXCHANNELS];
    float fTZero[VDST_MAXCHANNELS];
    float fTZeroMean = -999.;
    float fTZeroMedian = -999.;
    float fTZeroRMS = -999.;
    float fTZeroN = 0.;
    
    for( unsigned int i = 0; i < VDST_MAXCHANNELS; i++ )
    {
        fPed_high[i] = 0.;
        fPed_low[i] = 0.;
        fConv_high[i] = 0.;
        fConv_low[i] = 0.;
        fTZero[i] = -999.;
        for( unsigned int j = 0; j < VDST_MAXSUMWINDOW; j++ )
        {
            fPedvar_high[i * VDST_MAXSUMWINDOW + j] = 0.;
            fPedvar_low[i * VDST_MAXSUMWINDOW + j]  = 0.;
        }
    }
    // the following variables are needed for reading of external pedestal files
    float iT_sumwindow[VDST_MAXSUMWINDOW];
    float iT_pedvars[VDST_MAXSUMWINDOW];
    
    // list with root objects which are saved to file
    TList* hisList = new TList();
    
    // definition of calibration tree
    TTree* t = new TTree( "calibration", "calibration data" );
    hisList->Add( t );
    
    char hname[200];
    t->Branch( "TelID", &fTelID, "TelID/I" );
    t->Branch( "NPixel", &nPixel, "NPixel/i" );
    t->Branch( "num_sumwindow", &fnum_sumwindow, "num_sumwindow/i" );
    t->Branch( "sumwindow", fsumwindow, "sumwindow[num_sumwindow]/i" );
    t->Branch( "ped_high", fPed_high, "ped_high[NPixel]/F" );
    sprintf( hname, "pedvar_high[%d]/F", VDST_MAXCHANNELS * VDST_MAXSUMWINDOW );
    t->Branch( "pedvar_high", fPedvar_high, hname, VDST_MAXCHANNELS * VDST_MAXSUMWINDOW * 4 );
    t->Branch( "ped_low", fPed_low, "ped_low[NPixel]/F" );
    sprintf( hname, "pedvar_low[%d]/F", VDST_MAXCHANNELS * VDST_MAXSUMWINDOW );
    t->Branch( "pedvar_low", fPedvar_low, hname, VDST_MAXCHANNELS * VDST_MAXSUMWINDOW * 4 );
    t->Branch( "conv_high", fConv_high, "conv_high[NPixel]/F" );
    t->Branch( "conv_low", fConv_low, "conv_low[NPixel]/F" );
    t->Branch( "tzero", fTZero, "tzero[NPixel]/F" );
    t->Branch( "tzero_mean_tel", &fTZeroMean, "tzero_mean_tel/F" );
    t->Branch( "tzero_med_tel", &fTZeroMedian, "tzero_med_tel/F" );
    t->Branch( "tzero_rms_tel", &fTZeroRMS, "tzero_rms_tel/F" );
    t->Branch( "tzero_n_tel", &fTZeroN, "tzero_n_tel/F" );
    
    ///////////////////////////////////////////////////////
    // open external file with pedestals
    TFile* iPedFile = 0;
    if( ipedfile.size() > 0 )
    {
        iPedFile = new TFile( ipedfile.c_str(), "READ" );
        if( iPedFile->IsZombie() )
        {
            cout << "DST_fillCalibrationTree: error while opening external pedestal file: " << ipedfile << endl;
            exit( EXIT_FAILURE );
        }
    }
    
    ////////////////////////////////////////////////////////
    // loop over all telescopes and fill calibration trees
    for( int itel = 0; itel < hsdata->run_header.ntel; itel++ )
    {
        fTelID = hsdata->tel_moni[itel].tel_id;
        // select telescopes for this analysis
        if( telescope_list.size() == 0 || telescope_list.find( fTelID ) != telescope_list.end() )
        {
            cout << "\t filling calibration values for Telescope: " << itel << "\t" << fTelID;
            if( telescope_list[fTelID].TelescopeName.size() > 0 )
            {
                 cout << "\t" << telescope_list[fTelID].TelescopeName;
            }
            cout << " (FOV " << telescope_list[fTelID].FOV << " deg,";
            cout << " dynamic range: " << telescope_list[fTelID].DynamicRange;
            cout << ", RAWsum: " << telescope_list[fTelID].RAWsum;
            cout << ")" << endl;
            nPixel = ( unsigned int )hsdata->tel_moni[itel].num_pixels;
            if( VDST_MAXCHANNELS < nPixel )
            {
                cout << "DST_fillCalibrationTree error: number of pixels (" << nPixel;
                cout << ") exeeds allowed range (" << VDST_MAXCHANNELS << ")" << endl;
                cout << "\t adjust arrays..." << endl;
                exit( EXIT_FAILURE );
            }
            // pedestal (always taken from hessio file, as it might change from run to run)
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                if( hsdata->tel_moni[itel].num_ped_slices > 0. )
                {
                    fPed_high[p] = hsdata->tel_moni[itel].pedestal[HI_GAIN][p] / ( double )( hsdata->tel_moni[itel].num_ped_slices );
                    fPed_low[p] = hsdata->tel_moni[itel].pedestal[LO_GAIN][p] / ( double )( hsdata->tel_moni[itel].num_ped_slices );
                }
                else
                {
                    fPed_high[p] = 0.;
                    fPed_low[p] = 0.;
                }
            }
            
            ///////////////////////////////////////////////////////////////////////////////////
            // fill pedestal variances from external root file
            if( iPedFile && fTelescopeType.find( itel ) != fTelescopeType.end() )
            {
                std::ostringstream iSname;
                iSname << "tPeds_" << fTelescopeType[itel];
                TTree* iT = ( TTree* )iPedFile->Get( iSname.str().c_str() );
                if( !iT )
                {
                    cout << "DST_fillCalibrationTree error: pedestal tree not found for telescope ";
                    cout << itel << " (type " <<  fTelescopeType[itel] << ")" << endl;
                    continue;
                }
                // now copy values over to new tree
                iT->SetBranchAddress( "nsumwindows", &fnum_sumwindow );
                iT->SetBranchAddress( "sumwindow", iT_sumwindow );
                iT->SetBranchAddress( "pedvars", iT_pedvars );
                
                if( iT->GetEntries() < nPixel )
                {
                    cout << "DST_fillCalibrationTree error: number of pixels different in pedestal tree: ";
                    cout << nPixel << "\t" << iT->GetEntries() << endl;
                    return 0;
                }
                
                for( unsigned int p = 0; p < nPixel; p++ )
                {
                    iT->GetEntry( p );
                    
                    for( unsigned int w = 0; w < fnum_sumwindow; w++ )
                    {
                        fPedvar_high[p * VDST_MAXSUMWINDOW + w] = iT_pedvars[w] * sqrt( iNSBScaling );
                        fPedvar_low[p * VDST_MAXSUMWINDOW + w] = hsdata->tel_moni[itel].noise[LO_GAIN][p] * sqrt( iNSBScaling );
                    }
                }
                if( i_curDir )
                {
                    i_curDir->cd();
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            // read calibration data from simtel files
            // (not clear yet how to use this information later in the analysis)
            else
            {
                for( unsigned int p = 0; p < nPixel; p++ )
                {
                    fnum_sumwindow = 1;
                    fsumwindow[0] = ( unsigned int )hsdata->tel_moni[itel].num_ped_slices;
                    for( unsigned int w = 0; w < fnum_sumwindow; w++ )
                    {
                        fPedvar_high[p * VDST_MAXSUMWINDOW + w] = hsdata->tel_moni[itel].noise[HI_GAIN][p];
                        fPedvar_low[p * VDST_MAXSUMWINDOW + w] = hsdata->tel_moni[itel].noise[LO_GAIN][p];
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            // dc/pe conversion and pulse arrival time
            for( unsigned int p = 0; p < nPixel; p++ )
            {
                fConv_high[p] = hsdata->tel_lascal[itel].calib[HI_GAIN][p] * CALIB_SCALE;
                if( hsdata->tel_lascal[itel].known && hsdata->tel_lascal[itel].num_gains > 1 )
                {
                    fConv_low[p]  = hsdata->tel_lascal[itel].calib[LO_GAIN][p] * CALIB_SCALE;
                }
                else
                {
                    fConv_low[p] = -999.;
                }
                fTZero[p] = -999.;
                if( fData )
                {
                    fTZero[p] = fData->getDSTMeanPulseTiming( itel, p );
                }
            }
            if( fData )
            {
                fTZeroMean = fData->getDSTMeanPulseTimingPerTelescope( itel );
                fTZeroMedian = fData->getDSTMedianPulseTimingPerTelescope( itel );
                fTZeroRMS = fData->getDSTRMSPulseTimingPerTelescope( itel );
                fTZeroN = fData->getDSTNEventsPulseTimingPerTelescope( itel );
            }
            for( unsigned int w = 0; w < fnum_sumwindow; w++ )
            {
                fsumwindow[w] = ( unsigned int )( TMath::Nint( iT_sumwindow[w] ) );
            }
            t->Fill();
        }
    }
    // copy directories with pedestal distributions
    
    // loop over all telescope types
    //
    
    bool iFillPedDistributions = false;
    
    if( iFillPedDistributions )
    {
    
        cout << "copy directories with pedestal distributions..." << endl;
        set< ULong64_t >::iterator it;
        for( it = fTelescopeTypeList.begin(); it != fTelescopeTypeList.end(); ++it )
        {
            ULong64_t iTelescopeType = *it;
            std::ostringstream iDname;
            iDname << "distributions_" << iTelescopeType;
            if( iPedFile && iPedFile->cd( iDname.str().c_str() ) )
            {
                for( unsigned int w = 0; w < fnum_sumwindow; w++ )
                {
                
                    for( unsigned int p = 0; p < VDST_MAXCHANNELS; p++ )
                    {
                        std::ostringstream iHname;
                        iHname << "hpedPerTelescopeType_" << iTelescopeType << "_" << TMath::Nint( iT_sumwindow[w] ) << "_" << p ;
                        if( gDirectory->Get( iHname.str().c_str() ) )
                        {
                            hisList->Add( gDirectory->Get( iHname.str().c_str() ) );
                        }
                    }
                    /*
                     std::ostringstream iHname;
                     iHname << "hpedPerTelescopeType_" << iTelescopeType << "_" << TMath::Nint( iT_sumwindow[w] );
                     if( gDirectory->Get( iHname.str().c_str() ) )
                     {
                     hisList->Add( gDirectory->Get( iHname.str().c_str() ) );
                     }
                     */
                }
            }
            if( i_curDir )
            {
                i_curDir->cd();
            }
        }
    }
    
    // copy IPR graphs
    if( iPedFile )
    {
        i_curDir->cd();
        
        unsigned int iIPRgraphs = 0;
        
        // loop over all telescope types
        set< ULong64_t >::iterator it;
        for( it = fTelescopeTypeList.begin(); it != fTelescopeTypeList.end(); ++it )
        {
            ULong64_t iTelescopeType = *it;
            
            // loop over all summation windows (guess that there is a max)
            for( unsigned int w = 1; w < 100; w++ )
            {
                std::ostringstream iHname;
                iHname << "IPRcharge_TelType" << iTelescopeType << "_SW" << w;
                if( iPedFile->Get( iHname.str().c_str() ) )
                {
                    hisList->Add( iPedFile->Get( iHname.str().c_str() ) );
                    iIPRgraphs++;
                }
            }
        }
        if( iIPRgraphs > 0 )
        {
            cout << "found " << iIPRgraphs << " IPR graphs " << endl;
        }
    }
    
    
    delete [] fPedvar_high;
    delete [] fPedvar_low;
    
    return hisList;
}

/*
 * fill a tree with all information on the detector
 *
 */
TTree* DST_fill_detectorTree( AllHessData* hsdata, map< unsigned int, VDSTTelescopeConfiguration > telescope_list,
                              bool iApplyCameraScaling, map< ULong64_t, float > iCameraScalingMap )
{
    if( !hsdata )
    {
        return 0;
    }
    
    // HARDCODED: all large telescopes are parabolic
    //            all telescopes with a mirror area larger than this value are
    //            considered to be parabolic (in m^2)
    double fParabolic_mirrorArea = 380.;
    // HARDCODED: all telescopes with 2 mirrors only are SC telescopes
    int    fSC_number_of_mirrors = 2;
    cout << "Info:" << endl;
    cout << " assume that all telescopes with mirror area larger than " << fParabolic_mirrorArea;
    cout << " m^2 are of parabolic type" << endl;
    cout << " assume that all telescopes with " << fSC_number_of_mirrors;
    cout << " mirrors are Schwarzschild-Couder telescopes" << endl;
    
    // define tree
    int fTelID = 0;
    unsigned int fNTel = 0;
    Char_t fTelescopeName[300];
    float fTelxpos = 0.;
    float fTelypos = 0.;
    float fTelzpos = 0.;
    float fFocalLength = 0.;
	float fEffectiveFocalLength = 0.;
    float fCameraScaleFactor = 1.;
    float fCameraCentreOffset = 0.;
    float fCameraRotation = 0.;
    int   fCamCurvature = 0;
    unsigned int nPixel = 0;
    unsigned int nPixel_active = 0;
    unsigned int nSamples = 0;
    unsigned int fDynRange = 0;
    float Sample_time_slice = 0.;
    unsigned int nGains;
    float fHiLoScale = 0.;
    int   fHiLoThreshold = 0;
    float fHiLoOffset = 0.;
    const unsigned int fMaxPixel = 50000;
    float fXTubeMM[fMaxPixel];
    float fYTubeMM[fMaxPixel];
    float fRTubeMM[fMaxPixel];
    float fXTubeDeg[fMaxPixel];
    float fYTubeDeg[fMaxPixel];
    float fRTubeDeg[fMaxPixel];
    float fATubem2[fMaxPixel];
    unsigned int fPixelShape[fMaxPixel];
    int nDisabled = 0;
    int fTubeDisabled[fMaxPixel];
    float fMirrorArea = 0.;
    int fNMirrors = 0;
    float fFOV = 0.;
    ULong64_t fTelescope_type = 0;
    for( unsigned int i = 0; i < fMaxPixel; i++ )
    {
        fXTubeMM[i] = 0.;
        fYTubeMM[i] = 0.;
        fRTubeMM[i] = 0.;
        fXTubeDeg[i] = 0.;
        fYTubeDeg[i] = 0.;
        fRTubeDeg[i] = 0.;
        fTubeDisabled[i] = 0;
        fATubem2[i] = 0.;
    }
    
    TTree* fTreeDet = new TTree( "telconfig", "detector configuration" );
    
    fTreeDet->Branch( "NTel", &fNTel, "NTel/i" );
    fTreeDet->Branch( "TelID", &fTelID, "TelID/I" );
    fTreeDet->Branch( "TelescopeName", &fTelescopeName, "TelescopeName/C" );
    fTreeDet->Branch( "TelType", &fTelescope_type, "TelType/l" );
    fTreeDet->Branch( "TelX", &fTelxpos, "TelX/F" );
    fTreeDet->Branch( "TelY", &fTelypos, "TelY/F" );
    fTreeDet->Branch( "TelZ", &fTelzpos, "TelZ/F" );
    fTreeDet->Branch( "FocalLength", &fFocalLength, "FocalLength/F" );
	fTreeDet->Branch( "EffectiveFocalLength", &fEffectiveFocalLength, "EffectiveFocalLength/F" );
    fTreeDet->Branch( "FOV", &fFOV, "FOV/F" );
    fTreeDet->Branch( "CameraScaleFactor", &fCameraScaleFactor, "CameraScaleFactor/F" );
    fTreeDet->Branch( "CameraCentreOffset", &fCameraCentreOffset, "CameraCentreOffset/F" );
    fTreeDet->Branch( "CameraRotation", &fCameraRotation, "CameraRotation/F" );
    fTreeDet->Branch( "CamCurvature", &fCamCurvature, "CamCurvate/I" );
    fTreeDet->Branch( "NPixel", &nPixel, "NPixel/i" );
    fTreeDet->Branch( "NPixel_active", &nPixel_active, "NPixel_active/i" );
    fTreeDet->Branch( "NSamples", &nSamples, "NSamples/i" );
    fTreeDet->Branch( "Sample_time_slice", &Sample_time_slice, "Sample_time_slice/F" );
    fTreeDet->Branch( "NGains", &nGains, "NGains/i" );
    fTreeDet->Branch( "HiLoScale", &fHiLoScale, "HiLoScale/F" );
    fTreeDet->Branch( "HiLoThreshold", &fHiLoThreshold, "HiLoThreshold/I" );
    fTreeDet->Branch( "HiLoOffset", &fHiLoOffset, "HiLoOffset/F" );
    fTreeDet->Branch( "DynRange", &fDynRange, "DynRange/i" );
    fTreeDet->Branch( "XTubeMM", fXTubeMM, "XTubeMM[NPixel]/F" );
    fTreeDet->Branch( "YTubeMM", fYTubeMM, "YTubeMM[NPixel]/F" );
    fTreeDet->Branch( "RTubeMM", fRTubeMM, "RTubeMM[NPixel]/F" );
    fTreeDet->Branch( "XTubeDeg", fXTubeDeg, "XTubeDeg[NPixel]/F" );
    fTreeDet->Branch( "YTubeDeg", fYTubeDeg, "YTubeDeg[NPixel]/F" );
    fTreeDet->Branch( "RTubeDeg", fRTubeDeg, "RTubeDeg[NPixel]/F" );
    fTreeDet->Branch( "AreaTube_m2", fATubem2, "AreaTube_m2[NPixel]/F" );
    fTreeDet->Branch( "PixelShape", fPixelShape, "PixelShape[NPixel]/i" );
    fTreeDet->Branch( "NTubesOFF", &nDisabled, "NTubesOFF/I" );
    fTreeDet->Branch( "TubeOFF", fTubeDisabled, "TubeOFF[NPixel]/I" );
    fTreeDet->Branch( "NMirrors", &fNMirrors, "NMirrors/I" );
    fTreeDet->Branch( "MirrorArea", &fMirrorArea, "MirrorArea/F" );
    
    // check number of telescopes
    for( int itel = 0; itel <  hsdata->run_header.ntel; itel++ )
    {
        fTelID = hsdata->run_header.tel_id[itel];
        if( telescope_list.size() == 0 || telescope_list.find( fTelID ) != telescope_list.end() )
        {
            fNTel++;
        }
    }
    fNTel_checked = fNTel;
    
    ///////////////////////////////////////////////////////////
    // fill the tree with all information about the detector
    for( int itel = 0; itel <  hsdata->run_header.ntel; itel++ )
    {
        // actuall filling of telescope trees
        fTelID = hsdata->run_header.tel_id[itel];
        if( telescope_list.size() == 0 || telescope_list.find( fTelID ) != telescope_list.end() )
        {
            // Observe: transform to VERITAS coordinate system
            // x: East, y: North
            fTelxpos = -1.*hsdata->run_header.tel_pos[itel][1];
            fTelypos = hsdata->run_header.tel_pos[itel][0];
            fTelzpos = hsdata->run_header.tel_pos[itel][2];
            fFocalLength = hsdata->camera_set[itel].flen;
            // effective focal length (set only from prod4 on)
            if( hsdata->camera_set[itel].eff_flen > 0. )
            {
                 fEffectiveFocalLength = hsdata->camera_set[itel].eff_flen;
            }
            if( telescope_list[fTelID].DynamicRange > 0 )
            {
                fDynRange = ( unsigned int )telescope_list[fTelID].DynamicRange;
            }
            else
            {
                fDynRange = 0.;
            }
            fCameraScaleFactor = 1.;
            fCameraCentreOffset = 0.;
            fCameraRotation = -1.*hsdata->camera_set[itel].cam_rot * TMath::RadToDeg();
#if defined(CTA_PROD2)
            fCamCurvature = 0;
#else
            fCamCurvature = hsdata->camera_set[itel].curved_surface;
#endif
            fNMirrors = hsdata->camera_set[itel].num_mirrors;
            fMirrorArea = hsdata->camera_set[itel].mirror_area;
            
            nPixel = hsdata->camera_set[itel].num_pixels;
            nDisabled = hsdata->pixel_disabled[itel].num_HV_disabled;
            nSamples = hsdata->event.teldata[itel].raw->num_samples;
            Sample_time_slice = hsdata->pixel_set[itel].time_slice;
            nGains = hsdata->event.teldata[itel].raw->num_gains;
            fHiLoScale = hsdata->event.teldata[itel].raw->scale_hg8;
            fHiLoThreshold = hsdata->event.teldata[itel].raw->threshold;
            fHiLoOffset = hsdata->event.teldata[itel].raw->offset_hg8;
            
            // reset dead pixel list
            for( unsigned int p = 0; p < fMaxPixel; p++ )
            {
                fTubeDisabled[p] = 0;
            }
            
            nPixel_active = 0;
            float maxPix_dist = 0.;
            float pix_size = 0.;
            if( nPixel < fMaxPixel )
            {
                for( unsigned int p = 0; p < nPixel; p++ )
                {
                    fXTubeMM[p] = hsdata->camera_set[itel].xpix[p] * 1.e3;
                    fYTubeMM[p] = hsdata->camera_set[itel].ypix[p] * 1.e3;
#if defined(CTA_PROD3) || defined(CTA_MAX_SC) || defined(CTA_PROD3_MERGE)
                    // 0: circ., 1,3: hex, 2: square, -1: unknown
                    if( fPixelShape[p] >= 0 )
                    {
                        fPixelShape[p] = ( unsigned int )hsdata->camera_set[itel].pixel_shape[p];
                    }
                    else
#endif
                    {
                        fPixelShape[p] = 99;
                    }
#if defined(CTA_PROD3) || defined(CTA_MAX_SC) || defined(CTA_PROD3_MERGE)
                    // use as size the radius of the active area of the tube
                    if( hsdata->camera_set[itel].pixel_shape[p] == 0
                            || hsdata->camera_set[itel].pixel_shape[p] == 1 )
                    {
                        fRTubeMM[p] = sqrt( hsdata->camera_set[itel].area[p] / TMath::Pi() ) * 1.e3;
                    }
                    else
                    {
                        fRTubeMM[p] = sqrt( hsdata->camera_set[itel].area[p] ) * 1.e3 * 0.5;
                    }
#else
                    fRTubeMM[p] = sqrt( hsdata->camera_set[itel].area[p] / TMath::Pi() ) * 1.e3;
#endif
                    fATubem2[p] = hsdata->camera_set[itel].area[p];
                    if( p == 0 )
                    {
                        // prod2: pix_size in [deg] (diameter)
                        // (note that this is not entirely correct for the dual-mirror telescopes)
                        pix_size = atan2( ( double )hsdata->camera_set[itel].size[p], ( double )fFocalLength ) * TMath::RadToDeg();
                        
                    }
                    
                    // mm -> deg
                    fXTubeDeg[p] = atan2( ( double )fXTubeMM[p] / 1000., ( double )fFocalLength ) * TMath::RadToDeg();
                    fYTubeDeg[p] = atan2( ( double )fYTubeMM[p] / 1000., ( double )fFocalLength ) * TMath::RadToDeg();
                    fRTubeDeg[p] = atan2( ( double )fRTubeMM[p] / 1000., ( double )fFocalLength ) * TMath::RadToDeg();
                    
                    
                    float x2 = fXTubeDeg[p] * fXTubeDeg[p];
                    float y2 = fYTubeDeg[p] * fYTubeDeg[p];
                    if( sqrt( x2 + y2 ) * 2. > maxPix_dist )
                    {
                        maxPix_dist = sqrt( x2 + y2 ) * 2.;
                    }
                    // disable pixels which are too far out
                    if( telescope_list.size() != 0 && TMath::Abs( telescope_list[fTelID].FOV ) > 1.e-2 )
                    {
                        if( sqrt( x2 + y2 ) - fRTubeDeg[p] > telescope_list[fTelID].FOV * 0.5 )
                        {
                            fTubeDisabled[p] = 2;
                            nDisabled++;
                        }
                    }
                    if( fTubeDisabled[p] == 0 )
                    {
                        nPixel_active++;
                    }
                }
                // check for disabled pixels in sim_telarray
                for( int p_off = 0; p_off < hsdata->pixel_disabled[itel].num_HV_disabled; p_off++ )
                {
                    if( hsdata->pixel_disabled[itel].HV_disabled[p_off] < ( int )nPixel
                            && fTubeDisabled[hsdata->pixel_disabled[itel].HV_disabled[p_off]] == 0 )
                    {
                        fTubeDisabled[hsdata->pixel_disabled[itel].HV_disabled[p_off]] = 822;
                        nPixel_active--;
                    }
                }
            }
            if( telescope_list.size() == 0 || TMath::Abs( telescope_list[fTelID].FOV ) < 1.e5 )
            {
                fFOV = maxPix_dist;
            }
            else if( maxPix_dist < telescope_list[fTelID].FOV )
            {
                fFOV = maxPix_dist;
            }
            else
            {
                fFOV = telescope_list[fTelID].FOV;
            }
            // telescope name
            if( telescope_list.size() == 0 || telescope_list[fTelID].TelescopeName.size() == 0 )
            {
                sprintf( fTelescopeName, "Tel-%d", fTelID+1 );
            }
            else
            {
                sprintf( fTelescopeName, "%s", telescope_list[fTelID].TelescopeName.c_str() );
            }
            
            // telescope types
            fTelescope_type  = TMath::Nint( pix_size * 100. );
            fTelescope_type += TMath::Nint( fFOV * 10. ) * 100;
            fTelescope_type += TMath::Nint( fMirrorArea ) * 100 * 10 * 100;
            // all large telescopes are parabolic, all others are Davies-Cotton (hardwired)
            if( fMirrorArea > fParabolic_mirrorArea )
            {
                fTelescope_type += 100000000;
            }
            // Schwarzschild-Couder: check number of mirrors
            // (assumption is that SC telescope has 2 mirrors only)
            else if( fNMirrors == fSC_number_of_mirrors )
            {
                fTelescope_type += 200000000;
            }
            // Keep telescope IDs constant between prod3 and prod4 - 
            // despite changing parameters.

            // ASTRI
            if (fTelescope_type == 201411619)
            {
                fTelescope_type += 100000;
            }

            // TMP Prod4 fix for MSTs
            if (fTelescope_type == 10608618)
            {
                fTelescope_type = 10408618;
            }
            // (end of tmp)
            

            fTelescopeType[itel] = fTelescope_type;
            fTelescopeTypeList.insert( fTelescope_type );
            
            ////////////////////////////////////////////////
            // camera scaling (effective focal length)
            fCameraScaleFactor = 1.;
            if( fEffectiveFocalLength > 0. )
            {
                 fCameraScaleFactor = fFocalLength / fEffectiveFocalLength;
            }
            // scaling factor given from external file?
            else if( iCameraScalingMap.find( fTelescope_type ) != iCameraScalingMap.end() )
            {
                fCameraScaleFactor = 1. - iCameraScalingMap[fTelescope_type];
            }
            // hard wired value for DC telescopes
            // (for backwards compatibility,
            //  should not be used anymore)
            if( iApplyCameraScaling && fTelescope_type < 100000000 )
            {
                // assume that CTA MST is of intermediate-1.2 design
                if( TMath::Abs( fMirrorArea - 100. ) < 5. )
                {
                    // G.Hughes (2011/10/21): 2.79+-0.06 %
                    fCameraScaleFactor = 1. - 0.0279;
                }
            }
            
            fTreeDet->Fill();
        }
    }
    
    return fTreeDet;
}

/*

   main program

*/
int main( int argc, char** argv )
{
    // stop watch
    TStopwatch fStopWatch;
    fStopWatch.Start();
    
    IO_BUFFER* iobuf = NULL;
    IO_ITEM_HEADER item_header;
    const char* input_fname = NULL;
    int itel, rc = 0;
    int tel_id;
    int verbose = 0, ignore = 0, quiet = 0;
    int ntel_trg = 0, min_tel_trg = 0;
    int nev = 0, ntrg = 0;
    char* program = argv[0];
    int showdata = 0, showhistory = 0;
    size_t events = 0, max_events = 0;
    int iarg;
    
    fGlobalNTelPreviousEvent = 0;
    fGlobalMaxNumberofPixels = 0;
    fGlobalMaxNumberofSamples = 0;
    fGlobalTriggerReset = false;
    
    string config_file = "";             // file with list of telescopes
    string dst_file = "dst.root";        // output dst file
    string ped_file = "";                // file with pedestal and pedestal variances
    string trg_mask_dir = "";            // directory with trigger information (only for 2013 prod2 MC)
    unsigned int fWriteFADC = 2;         // fill FADC traces into converter (0=QADC, 1=FADC, 2=mix of both)
    bool   fApplyCameraScaling = true;   // apply camera plate scaling according for DC telescopes
    map< ULong64_t, float > fCameraScalingMap; // camera plate scale (per telescoe type; read from external file)
    float  fNSB_scaling = 1.;            // pedvar scaling due to higher NSB
    
    static AllHessData* hsdata;
    
    cout << endl;
    cout << "CTA.convert_hessio_to_VDST: A program to convert hessio data to EVNDISP DST files";
    cout << " (" << VGlobalRunParameter::fEVNDISP_VERSION << ")" << endl;
    cout << " (based on a skeleton program distributed with the hessio package)" << endl;
    cout << "=====================================================================" << endl;
    
    /* Show command line on output */
    if( getenv( "SHOWCOMMAND" ) != NULL )
    {
        for( iarg = 0; iarg < argc; iarg++ )
        {
            printf( "%s ", argv[iarg] );
        }
        printf( "\n" );
    }
    
    /* Catch INTerrupt and TERMinate signals to stop program */
    signal( SIGINT, stop_signal_function );
    signal( SIGTERM, stop_signal_function );
    
    if( argc < 2 )
    {
        syntax( program );
    }
    interrupted = 0;
    
    /* Check assumed limits with the ones compiled into the library. */
    H_CHECK_MAX();
    
    if( ( iobuf = allocate_io_buffer( 1000000L ) ) == NULL )
    {
        cout << "Cannot allocate I/O buffer" << endl;
        exit( EXIT_FAILURE );
    }
    iobuf->max_length = 100000000L;
    
    
    /* Command line options */
    while( argc > 1 )
    {
        if( strcmp( argv[1], "-i" ) == 0 )
        {
            ignore = 1;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-v" ) == 0 )
        {
            verbose = 1;
            quiet = 0;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-q" ) == 0 )
        {
            quiet = 1;
            verbose = 0;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-h" ) == 0 || strcmp( argv[1], "--history" ) == 0 )
        {
            showhistory = 1;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-s" ) == 0 )
        {
            showdata = 1;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-S" ) == 0 )
        {
            showdata = 1;
            putenv( ( char* )"PRINT_VERBOSE=1" );
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "--max-events" ) == 0 && argc > 2 )
        {
            max_events = atol( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-NSB" ) == 0 )
        {
            fNSB_scaling = atof( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-a" ) == 0 )
        {
            config_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-o" ) == 0 )
        {
            dst_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-t" ) == 0 )
        {
            trg_mask_dir = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-c" ) == 0 )
        {
            ped_file = argv[2];
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-f" ) == 0 )
        {
            fWriteFADC = atoi( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-pe" ) == 0 )
        {
            fFillPELeaf = true;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-peakadc" ) == 0 )
        {
            fFillPeakADC = true;
            argc--;
            argv++;
            continue;
        }
        else if( strcmp( argv[1], "-pedshift" ) == 0 )
        {
            fPedestalShift = atof( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-minenergy" ) == 0 )
        {
            fMinEnergy = atof( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "-maxenergy" ) == 0 )
        {
            fMaxEnergy = atof( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        // apply camera scaling for MSTs only
        else if( strcmp( argv[1], "-r" ) == 0 )
        {
            fApplyCameraScaling = atoi( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        // file with camera scaling values per telescope type
        else if( strcmp( argv[1], "-rfile" ) == 0 )
        {
            fCameraScalingMap = readCameraScalingMap( argv[2] );
            argc -= 2;
            argv += 2;
            continue;
        }
        else if( strcmp( argv[1], "--help" ) == 0 )
        {
            printf( "\nc_DST: A program to convert hessio data to EVNDISP DST files.\n\n" );
            syntax( program );
        }
        else if( argv[1][0] == '-' && argv[1][1] != '\0' )
        {
            printf( "Syntax error at '%s'\n", argv[1] );
            syntax( program );
        }
        else
        {
            break;
        }
    }
    
    if( verbose && !quiet )
    {
        showhistory = 1;
    }
    
    //////////////////////////////////////////////////////////////////
    // initialize eventdisplay dst and run header
    ///////////////////////////////////////////////////////////////////
    
    cout << endl << "NOTE: FIXED TIMING LEVELS READ FROM HESSIO FILE" << endl << endl;
    
    ///////////////////////////////////////////////////////////////////
    // open DST file
    TFile* fDSTfile = new TFile( dst_file.c_str(), "RECREATE" );
    if( fDSTfile->IsZombie() )
    {
        cout << "Error while opening DST output file: " << fDSTfile->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    cout << "DST tree will be written to " << dst_file << endl;
    if( fWriteFADC == 1 )
    {
        cout << "Writing FADC samples to DST file" << endl;
    }
    else if( fWriteFADC == 2 )
    {
        cout << "Writing FADC samples / QADC to DST file" << endl;
    }
    else
    {
        cout << "No FADC output" << endl;
    }
    if( fApplyCameraScaling )
    {
        cout << "Apply camera plate scaling for DC (intermediate) telescopes" << endl;
    }
    if( fMinNumber_pe_per_telescope > 0 )
    {
        cout << "Add trigger type #4 for events with >" << fMinNumber_pe_per_telescope << " pe per telescope" << endl;
    }
    if( fPedestalShift > 0. )
    {
        cout << "Apply pedestal shift by " << fPedestalShift << " dc ";
        cout << "(applied only to non-FADC data)" << endl;
    }
    if( fMinEnergy > 0. )
    {
        cout << "Apply minimum energy cut of " << fMinEnergy << " TeV " << endl;
    }
    if( fMaxEnergy < 1.e19 )
    {
        cout << "Apply maximum energy cut of " << fMaxEnergy << " TeV " << endl;
    }
    
    // new DST tree
    VDSTTree* fDST = new VDSTTree();
    fDST->setMC();
    map< unsigned int, VDSTTelescopeConfiguration > fTelescope_list = fDST->readArrayConfig( config_file );
    if( fTelescope_list.size() < 1 )
    {
        cout << "error reading array configuration file" << endl;
        cout << "(" << config_file << ")" << endl;
        exit( EXIT_FAILURE );
    }
    fDST->setFADC( fWriteFADC );
    fDST->setFillPELeaf( fFillPELeaf );
    if( fFillPELeaf )
    {
        cout << "Filling pe leaf " << endl;
    }
    fDST->setFillPeakADFLeaf( fFillPeakADC );
    if( !fWriteTelConfigTreeOnly )
    {
        fDST->initDSTTree( false );
        fDST->initMCTree();
    }
    
    // MC run header
    VMonteCarloRunHeader* fMC_header = new VMonteCarloRunHeader();
    fMC_header->SetName( "MC_runheader" );
    
    // long list of input file names
    string fRunPara_InputFileName = "";
    unsigned int fRunPara_NumInputFiles = 0;
    
    cout << "start proceeding with sim_telarray file" << endl;
    
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    // Now go over rest of the command line and read input file names
    while( argc > 1 || input_fname != NULL )
    {
        if( interrupted )
        {
            break;
        }
        if( argc > 1 )
        {
            if( argv[1][0] == '-' && argv[1][1] != '\0' )
            {
                syntax( program );
            }
            else
            {
                input_fname = argv[1];
                argc--;
                argv++;
            }
        }
        
        if( strcmp( input_fname , "-" ) == 0 )
        {
            iobuf->input_file = stdin;
        }
        else if( ( iobuf->input_file = fileopen( input_fname, READ_BINARY ) ) == NULL )
        {
            perror( input_fname );
            cout << "Cannot open input file." << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        cout << "opening simtel file " << input_fname << endl;
        
        fflush( stdout );
        fprintf( stderr, "%s\n", input_fname );
        string f_inputfilename = "";
        if( input_fname )
        {
            f_inputfilename = input_fname;
        }
        printf( "\nInput file '%s' has been opened.\n", input_fname );
        input_fname = NULL;
        
        fRunPara_InputFileName += f_inputfilename;
        fRunPara_NumInputFiles++;
        if( fRunPara_NumInputFiles > 0 )
        {
            fRunPara_InputFileName += ", ";
        }
        
        /////////////////////////////////////////////
        // read in trigger mask
        if( trg_mask_dir.size() > 0 )
        {
            string iFile = f_inputfilename;
            if( iFile.find_last_of( "/\\" ) != string::npos )
            {
                iFile = iFile.substr( iFile.find_last_of( "/\\" ) + 1,
                                      iFile.find( "simtel.gz" ) - iFile.find_last_of( "/\\" ) - 1 );
            }
            string trg_mask_file = trg_mask_dir + iFile + "trgmask.gz";
            read_trigger_mask( trg_mask_file );
        }
        
        
        /////////////////////////////////////////////
        // Loop over all data in the input data file
        for( ;; )
        {
            if( interrupted )
            {
                break;
            }
            
            /* Find and read the next block of data. */
            /* In case of problems with the data, just give up. */
            if( find_io_block( iobuf, &item_header ) != 0 )
            {
                break;
            }
            if( max_events > 0 && events >= max_events )
            {
                if( iobuf->input_file != stdin )
                {
                    break;
                }
                if( skip_io_block( iobuf, &item_header ) != 0 )
                {
                    break;
                }
                continue;
            }
            if( read_io_block( iobuf, &item_header ) != 0 )
            {
                break;
            }
            
            if( hsdata == NULL &&
                    item_header.type > IO_TYPE_HESS_RUNHEADER &&
                    item_header.type < IO_TYPE_HESS_RUNHEADER + 200 )
            {
                fprintf( stderr, "Trying to read event data before run header.\n" );
                fprintf( stderr, "Skipping this data block.\n" );
                continue;
            }
            
            //////////////////////////////////////////
            // check header types
            switch( ( int ) item_header.type )
            {
                /* =================================================== */
                case IO_TYPE_HESS_RUNHEADER:
                    /* Summary of a preceding run in the same file ? */
                    if( nev > 0 )
                    {
                        printf( "%d of %d events triggered.\n", ntrg, nev );
                    }
                    nev = ntrg = 0;
                    /* Structures might be allocated from previous run */
                    if( hsdata != NULL )
                    {
                        /* Free memory allocated inside ... */
                        for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
                        {
                            if( hsdata->event.teldata[itel].raw != NULL )
                            {
                                free( hsdata->event.teldata[itel].raw );
                                hsdata->event.teldata[itel].raw = NULL;
                            }
                            if( hsdata->event.teldata[itel].pixtm != NULL )
                            {
                                free( hsdata->event.teldata[itel].pixtm );
                                hsdata->event.teldata[itel].pixtm = NULL;
                            }
                            if( hsdata->event.teldata[itel].img != NULL )
                            {
                                free( hsdata->event.teldata[itel].img );
                                hsdata->event.teldata[itel].img = NULL;
                            }
                        }
                        /* Free main structure */
                        free( hsdata );
                        hsdata = NULL;
                        
                        /* Perhaps some cleaning needed in ROOT as well ... */
                        
                    }
                    hsdata = ( AllHessData* ) calloc( 1, sizeof( AllHessData ) );
                    if( ( rc = read_hess_runheader( iobuf, &hsdata->run_header ) ) < 0 )
                    {
                        cout << "Reading run header failed." << endl;
                        exit( EXIT_FAILURE );
                    }
                    if( !quiet )
                    {
                        printf( "Reading simulated data for %d telescope(s)\n", hsdata->run_header.ntel );
                    }
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_runheader(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_runheader( iobuf );
                    }
                    
                    cout << "Telescope configuration: " << endl;
                    for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
                    {
                        tel_id = hsdata->run_header.tel_id[itel];
                        hsdata->camera_set[itel].tel_id = tel_id;
                        fTelescopeSimTelList[tel_id] = itel;
                        cout << "\t initialize Telescope ID " << tel_id << " (is telescope # " << itel << ")" << endl;
                        hsdata->camera_org[itel].tel_id = tel_id;
                        hsdata->pixel_set[itel].tel_id = tel_id;
                        hsdata->pixel_disabled[itel].tel_id = tel_id;
                        hsdata->cam_soft_set[itel].tel_id = tel_id;
                        hsdata->tracking_set[itel].tel_id = tel_id;
                        hsdata->point_cor[itel].tel_id = tel_id;
                        hsdata->event.num_tel = hsdata->run_header.ntel;
                        hsdata->event.teldata[itel].tel_id = tel_id;
                        hsdata->event.trackdata[itel].tel_id = tel_id;
                        if( ( hsdata->event.teldata[itel].raw =
                                    ( AdcData* ) calloc( 1, sizeof( AdcData ) ) ) == NULL )
                        {
                            cout << "Not enough memory" << endl;
                            exit( EXIT_FAILURE );
                        }
                        hsdata->event.teldata[itel].raw->tel_id = tel_id;
                        if( ( hsdata->event.teldata[itel].pixtm =
                                    ( PixelTiming* ) calloc( 1, sizeof( PixelTiming ) ) ) == NULL )
                        {
                            cout << "Not enough memory" << endl;
                            exit( EXIT_FAILURE );
                        }
                        hsdata->event.teldata[itel].pixtm->tel_id = tel_id;
                        if( ( hsdata->event.teldata[itel].img =
                                    ( ImgData* ) calloc( 2, sizeof( ImgData ) ) ) == NULL )
                        {
                            cout << "Not enough memory" << endl;
                            exit( EXIT_FAILURE );
                        }
                        hsdata->event.teldata[itel].max_image_sets = 2;
                        hsdata->event.teldata[itel].img[0].tel_id = tel_id;
                        hsdata->event.teldata[itel].img[1].tel_id = tel_id;
                        hsdata->tel_moni[itel].tel_id = tel_id;
                        hsdata->tel_lascal[itel].tel_id = tel_id;
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MCRUNHEADER:
                    rc = read_hess_mcrunheader( iobuf, &hsdata->mc_run_header );
                    
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_mcrunheader(), rc = %d\n", rc );
                    }
                    //            if ( showdata )
                    print_hess_mcrunheader( iobuf );
                    
                    // fill EVNDISP DST run header
                    DST_fillMCRunheader( fMC_header, hsdata, ( fRunPara_NumInputFiles > 1 ) );
                    break;
                    
                /* =================================================== */
                case IO_TYPE_MC_INPUTCFG:
                {
                    struct linked_string corsika_inputs;
                    corsika_inputs.text = NULL;
                    corsika_inputs.next = NULL;
                    read_input_lines( iobuf, &corsika_inputs );
                    if( corsika_inputs.text != NULL )
                    {
                        struct linked_string* xl = NULL, *xln = NULL;
                        if( ! quiet )
                        {
                            printf( "\nCORSIKA was run with the following input lines:\n" );
                        }
                        for( xl = &corsika_inputs; xl != NULL; xl = xln )
                        {
                            if( ! quiet )
                            {
                                printf( "   %s\n", xl->text );
                            }
                            free( xl->text );
                            xl->text = NULL;
                            xln = xl->next;
                            xl->next = NULL;
                            if( xl != &corsika_inputs )
                            {
                                free( xl );
                            }
                        }
                    }
                }
                break;
                
                /* =================================================== */
                case 70: /* How sim_hessarray was run and how it was configured. */
                    if( showhistory )
                    {
                        list_history( iobuf, NULL );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CAMSETTINGS:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Camera settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_camsettings( iobuf, &hsdata->camera_set[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_camsettings(), rc = %d\n", rc );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CAMORGAN:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Camera organisation for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_camorgan( iobuf, &hsdata->camera_org[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_camorgan(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_camorgan( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_PIXELSET:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Pixel settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_pixelset( iobuf, &hsdata->pixel_set[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_pixelset(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_pixelset( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_PIXELDISABLE:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Pixel disable block for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_pixeldis( iobuf, &hsdata->pixel_disabled[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_pixeldis(), rc = %d\n", rc );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CAMSOFTSET:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Camera software settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_camsoftset( iobuf, &hsdata->cam_soft_set[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_camsoftset(), rc = %d\n", rc );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_POINTINGCOR:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Pointing correction for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_pointingcor( iobuf, &hsdata->point_cor[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_pointingco(), rc = %d\n", rc );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_TRACKSET:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Tracking settings for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_trackset( iobuf, &hsdata->tracking_set[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_trackset(), rc = %d\n", rc );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_EVENT:
                    rc = read_hess_event( iobuf, &hsdata->event, -1 );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_event(), rc = %d\n", rc );
                    }
                    events++;
                    if( showdata )
                    {
                        print_hess_event( iobuf );
                    }
                    /* Count number of telescopes (still) present in data and triggered */
                    ntel_trg = 0;
                    for( itel = 0; itel < hsdata->run_header.ntel; itel++ )
                    {
                        if( hsdata->event.teldata[itel].known )
                        {
                            /* If non-triggered telescopes record data (like HEGRA),
                               we may have to check the central trigger bit as well,
                               but ignore this for now. */
                            ntel_trg++;
                        }
                    }
                    if( hsdata->event.shower.known )
                    {
                        hsdata->event.shower.num_trg = ntel_trg;
                    }
                    if( ntel_trg < min_tel_trg )
                    {
                        continue;
                    }
                    ntrg++;
                    
                    // fill EVNDISP DST event
                    if( !fWriteTelConfigTreeOnly )
                    {
                        DST_fillEvent( fDST, hsdata, fTelescope_list, fWriteFADC );
                    }
                    
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_CALIBEVENT:
                {
                    int type = -1;
                    rc = read_hess_calib_event( iobuf, &hsdata->event, -1, &type );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_calib_event(), rc = %d, type=%d\n", rc, type );
                    }
                    if( !fWriteTelConfigTreeOnly )
                    {
                        DST_fillEvent( fDST, hsdata, fTelescope_list, fWriteFADC );
                    }
                }
                break;
                
                /* =================================================== */
                case IO_TYPE_HESS_MC_SHOWER:
                    rc = read_hess_mc_shower( iobuf, &hsdata->mc_shower );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_mc_shower(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_mc_shower( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MC_EVENT:
                    rc = read_hess_mc_event( iobuf, &hsdata->mc_event );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_mc_event(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_mc_event( iobuf );
                    }
                    
                    DST_fillMCEvent( fDST, hsdata );
                    
                    break;
                    
                /* =================================================== */
                case IO_TYPE_MC_TELARRAY:
                    if( hsdata && hsdata->run_header.ntel > 0 )
                    {
                        rc = read_hess_mc_phot( iobuf, &hsdata->mc_event );
                        if( verbose || rc != 0 )
                        {
                            printf( "read_hess_mc_phot(), rc = %d\n", rc );
                        }
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MC_PE_SUM:
                    rc = read_hess_mc_pe_sum( iobuf, &hsdata->mc_event.mc_pesum );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_mc_pe_sum(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_mc_pe_sum( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_TEL_MONI:
                    // Telescope ID among others in the header
                    tel_id = ( item_header.ident & 0xff ) |
                             ( ( item_header.ident & 0x3f000000 ) >> 16 );
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Telescope monitor block for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_tel_monitor( iobuf, &hsdata->tel_moni[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_tel_monitor(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_tel_monitor( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_LASCAL:
                    tel_id = item_header.ident; // Telescope ID is in the header
                    if( ( itel = find_tel_idx( tel_id ) ) < 0 )
                    {
                        char msg[256];
                        snprintf( msg, sizeof( msg ) - 1,
                                  "Laser/LED calibration for unknown telescope %d.", tel_id );
                        cout << msg << endl;
                        exit( EXIT_FAILURE );
                    }
                    rc = read_hess_laser_calib( iobuf, &hsdata->tel_lascal[itel] );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_laser_calib(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_laser_calib( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_RUNSTAT:
                    rc = read_hess_run_stat( iobuf, &hsdata->run_stat );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_run_stat(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_run_stat( iobuf );
                    }
                    break;
                    
                /* =================================================== */
                case IO_TYPE_HESS_MC_RUNSTAT:
                    rc = read_hess_mc_run_stat( iobuf, &hsdata->mc_run_stat );
                    if( verbose || rc != 0 )
                    {
                        printf( "read_hess_mc_run_stat(), rc = %d\n", rc );
                    }
                    if( showdata )
                    {
                        print_hess_mc_run_stat( iobuf );
                    }
                    break;
                    
                default:
                    if( !ignore )
                    {
                        fprintf( stderr, "Ignoring data block type %ld\n", item_header.type );
                    }
            }
        }
        
        if( iobuf->input_file != NULL && iobuf->input_file != stdin )
        {
            fileclose( iobuf->input_file );
        }
        iobuf->input_file = NULL;
        reset_io_block( iobuf );
        
        if( nev > 0 )
        {
            printf( "%d of %d events triggered\n", ntrg, nev );
        }
        
        if( hsdata != NULL )
        {
            hsdata->run_header.run = 0;
        }
    }
    // end while loop over all input files
    ////////////////////////////////////////////////////
    
    cout << "writing DST, MC and detector trees" << endl;
    if( fDST && fDST->getDSTTree() )
    {
        fDST->getDSTTree()->Write();
        cout << "\t (writing " << fDST->getDSTTree()->GetEntries() << " events)" << endl;
    }
    if( fDST && fDST->getMCTree() )
    {
        cout << "\t (writing " << fDST->getMCTree()->GetEntries() << " MC events)" << endl;
        fDST->getMCTree()->Write();
    }
    // writing detector tree
    TTree* i_detTree = DST_fill_detectorTree( hsdata, fTelescope_list, fApplyCameraScaling, fCameraScalingMap );
    if( i_detTree )
    {
        i_detTree->Write();
    }
    // writing Monte Carlo header
    if( fMC_header )
    {
        fMC_header->Write();
        fMC_header->print();
    }
    // writing calibration data
    TList* i_calibTree = DST_fillCalibrationTree( fDST, hsdata, fTelescope_list, ped_file, fNSB_scaling );
    if( fDSTfile )
    {
        fDSTfile->cd();
    }
    if( i_calibTree )
    {
        i_calibTree->Write();
    }
    ///////////////////////////////////////////////////
    // writing run parameters to dst file
    VEvndispRunParameter* fRunPara = new VEvndispRunParameter( false );
    fRunPara->SetName( "runparameterDST" );
    fRunPara->setSystemParameters();
    fRunPara->fEventDisplayUser = "CTA-DST";
    fRunPara->frunnumber = fDST->getDSTRunNumber();
    fRunPara->getObservatory() = "CTA";
    fRunPara->fsourcetype = 7;                                  // 7 is DST - MC
    fRunPara->fsourcefile = fRunPara_InputFileName;
    fRunPara->fIsMC = true;
    fRunPara->fCalibrationDataType = 0;                         // no pedvars available
    fRunPara->fFADCChargeUnit = "PE";
    fRunPara->fNTelescopes = fNTel_checked;
    if( fRunPara->fNTelescopes == 0 && i_detTree )
    {
        fRunPara->fNTelescopes = ( unsigned int )i_detTree->GetEntries();
    }
    fRunPara->fpulsetiminglevels = getPulseTimingLevels();
    fRunPara->setPulseZeroIndex();
    fRunPara->Write();
    ///////////////////////////////////////////////////
    
    if( fDSTfile )
    {
        fDSTfile->Close();
    }
    
    
    fStopWatch.Stop();
    fStopWatch.Print();
    
    cout << "MC events written to DST file: " << dst_file << endl;
    
    cout << "exit..." << endl;
    
    return 0;
}

