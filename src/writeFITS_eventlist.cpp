/*  writeFITS_eventlist.cpp
 *
 *  convert eventdisplay ROOT output to FITS format
 *
 *  based on evlio package written by K.Kosack and M.Raue
 *  (http://evlio.readthedocs.org/en/latest/)
 */


#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "stdint.h"

#include <iostream>

#include "FITSRecord.hh"

#include "CData.h"
#include "VEvndispRunParameter.h"
#include "VSkyCoordinates.h"
#include "VSkyCoordinatesUtilities.h"

using namespace std;

int writeTELARRAY( string iInputFile, string iOutputFile )
{

    // read templates from template directory
    char* TEMPLATE_DIR;
    TEMPLATE_DIR = getenv( "EVLIOSYS" );
    
    char templ[1200];
    sprintf( templ, "%s/../templates/evl/1.0.0/EventList.tpl", TEMPLATE_DIR );
    
    ///////////////////////////////
    // open and read data file
    TFile* mscwfile = new TFile( iInputFile.c_str() );
    if( mscwfile->IsZombie() )
    {
        cout << "error in reading mscw file: " << iInputFile << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    //////////////////////////////////////
    // header
    
    
    //////////////////////////////////////
    // array and telescope configuration
    
    TTree* fTreeDet = ( TTree* )mscwfile->Get( "telconfig" );
    if( !fTreeDet )
    {
        cout << "error while reading mscw file: no tree for telescope configuration (telconfig)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    
    // define tree
    UInt_t fNTel = 0;
    int fTelID = 0;
    float fTelxpos = 0.;
    float fTelypos = 0.;
    float fTelzpos = 0.;
    unsigned int fTelID_hyperArray = 0;
    ULong64_t fTelType = 0;
    float fFocalLength = 0.;
    float fFOV = 0.;
    float fCameraScaleFactor = 0.;
    float fCameraCentreOffset = 0.;
    float fCameraRotation = 0.;
    unsigned int nPixel = 0;
    unsigned int nSamples = 0;
    unsigned int nGains = 2;
    float fHiLoScale = 0.;
    int   fHiLoThreshold = 0;
    float fHiLoOffset = 0.;
    float fXTubeMM[VDST_MAXCHANNELS];
    float fYTubeMM[VDST_MAXCHANNELS];
    float fRTubeMM[VDST_MAXCHANNELS];
    float fXTubeDeg[VDST_MAXCHANNELS];
    float fYTubeDeg[VDST_MAXCHANNELS];
    float fRTubeDeg[VDST_MAXCHANNELS];
    float fATube_m2[VDST_MAXCHANNELS];
    int   fTubeOFF[VDST_MAXCHANNELS];
    float fMirrorArea = 0.;
    int   fNMirrors = 0;
    for( unsigned int i = 0; i < VDST_MAXCHANNELS; i++ )
    {
        fXTubeMM[i] = 0.;
        fYTubeMM[i] = 0.;
        fRTubeMM[i] = 0.;
        fXTubeDeg[i] = 0.;
        fYTubeDeg[i] = 0.;
        fRTubeDeg[i] = 0.;
        fATube_m2[i] = 0.;
        fTubeOFF[i] = 0;
    }
    
    fTreeDet->SetBranchAddress( "NTel", &fNTel );
    fTreeDet->SetBranchAddress( "TelID", &fTelID );
    fTreeDet->SetBranchAddress( "TelType", &fTelType );
    fTreeDet->SetBranchAddress( "TelX", &fTelxpos );
    fTreeDet->SetBranchAddress( "TelY", &fTelypos );
    fTreeDet->SetBranchAddress( "TelZ", &fTelzpos );
    fTreeDet->SetBranchAddress( "FocalLength", &fFocalLength );
    if( fTreeDet->GetBranchStatus( "FOV" ) )
    {
        fTreeDet->SetBranchAddress( "FOV", &fFOV );
    }
    fTreeDet->SetBranchAddress( "CameraScaleFactor", &fCameraScaleFactor );
    fTreeDet->SetBranchAddress( "CameraCentreOffset", &fCameraCentreOffset );
    fTreeDet->SetBranchAddress( "CameraRotation", &fCameraRotation );
    fTreeDet->SetBranchAddress( "NPixel", &nPixel );
    fTreeDet->SetBranchAddress( "NSamples", &nSamples );
    if( fTreeDet->GetBranchStatus( "NGains" ) )
    {
        fTreeDet->SetBranchAddress( "NGains", &nGains );
    }
    if( fTreeDet->GetBranchStatus( "HiLoScale" ) )
    {
        fTreeDet->SetBranchAddress( "HiLoScale", &fHiLoScale );
    }
    if( fTreeDet->GetBranchStatus( "HiLoThreshold" ) )
    {
        fTreeDet->SetBranchAddress( "HiLoThreshold", &fHiLoThreshold );
    }
    if( fTreeDet->GetBranchStatus( "HiLoOffset" ) )
    {
        fTreeDet->SetBranchAddress( "HiLoOffset", &fHiLoOffset );
    }
    fTreeDet->SetBranchAddress( "XTubeMM", fXTubeMM );
    fTreeDet->SetBranchAddress( "YTubeMM", fYTubeMM );
    fTreeDet->SetBranchAddress( "RTubeMM", fRTubeMM );
    fTreeDet->SetBranchAddress( "XTubeDeg", fXTubeDeg );
    fTreeDet->SetBranchAddress( "YTubeDeg", fYTubeDeg );
    fTreeDet->SetBranchAddress( "RTubeDeg", fRTubeDeg );
    if( fTreeDet->GetBranchStatus( "AreaTube_m2" ) )
    {
        fTreeDet->SetBranchAddress( "AreaTube_m2", fATube_m2 );
    }
    if( fTreeDet->GetBranchStatus( "TubeOFF" ) )
    {
        fTreeDet->SetBranchAddress( "TubeOFF", fTubeOFF );
    }
    if( fTreeDet->GetBranchStatus( "NMirrors" ) )
    {
        fTreeDet->SetBranchAddress( "NMirrors", &fNMirrors );
    }
    if( fTreeDet->GetBranchStatus( "MirrorArea" ) )
    {
        fTreeDet->SetBranchAddress( "MirrorArea", &fMirrorArea );
    }
    
    ///////////////////////////////////////////////////////////////////
    // FITS record: telescope array
    FITSRecord recTelArray( iOutputFile, templ, "TELARRAY" );
    recTelArray.setVerbose( 0 );
    
    int16_t telid = 0;
    const size_t NAMELEN = 10;
    //char telclass[NAMELEN];
    long long int telclass[NAMELEN] = {0};
    long long int subclass[NAMELEN] = {0};
    unsigned long telsubclass = 0;
    double pos_x = 0;
    double pos_y = 0;
    double pos_z = 0;
    double mirror_area = 110.0;
    double cam_area = 10.0;
    double focal_length = 1.0;
    double fieldofview = 3.5;
    double pix_size = 0.14;
    double pix_sep = fieldofview / pix_size;
    double npix = 0.;
    
    std::string telescop = "VERITAS";
    std::string array; //  name of array layout for reference
    int obs_id; //  name of observation this array corresponds to, or 0 if general
    int geolat; //  latitude of observatory
    int geolon; //  longitude of observatory
    int altitude; //  altitude of observatory (km)
    
    recTelArray.writeHeader( "TELESCOP", telescop );
    
    recTelArray.mapColumnToVar( "TELID", telid );
    recTelArray.mapColumnToVar( "TELCLASS", telclass );
    recTelArray.mapColumnToVar( "SUBCLASS", subclass );
    recTelArray.mapColumnToVar( "POSX", pos_x );
    recTelArray.mapColumnToVar( "POSY", pos_y );
    recTelArray.mapColumnToVar( "POSZ", pos_z );
    recTelArray.mapColumnToVar( "MIRAREA", mirror_area );
    recTelArray.mapColumnToVar( "CAMAREA", cam_area );
    recTelArray.mapColumnToVar( "FOCLEN", focal_length );
    recTelArray.mapColumnToVar( "FOV", fieldofview );
    recTelArray.mapColumnToVar( "N_PIX", npix );
    recTelArray.mapColumnToVar( "PIX_SIZE", pix_size );
    recTelArray.mapColumnToVar( "PIX_SEP", pix_sep );
    
    for( unsigned int i = 0; i < fTreeDet->GetEntries(); i++ )
    {
    
        fTreeDet->GetEntry( i );
        
        telid		= fTelID;
        telclass[i]      = ( long long int )fTelType;
        subclass[i]      = ( long long int )fTelType;
        pos_x   		= fTelxpos;
        pos_y   		= fTelypos;
        pos_z   		= fTelzpos;
        mirror_area      = fMirrorArea;
        cam_area         = fMirrorArea;
        focal_length     = fFocalLength;
        fieldofview      = fFOV;
        npix    		= nPixel;
        if( fRTubeDeg[0] > 0. )
        {
            pix_size = 2.*fRTubeDeg[0];
        }
        pix_sep 		= 2.*fRTubeDeg[0];
        
        recTelArray.write();
    }
    recTelArray.finishWriting();
    
    /////////////////////////////////////////////////////////////////////
    
    /*    int    event_id = -9999;
        bool flags = 0;
        int    multip = -9999;
        bool    telmask = 1;
        double ttime = -9999.0;
        double ra = -9999.0;
        double dec = -9999.0;
        double dir_err = -9999.0;
        double detx = -9999.0;
        double dety = -9999.0;
        double alt_pnt = -9999.0;
        double az_pnt = -9999.0;
        double alt = -9999.0;
        double az = -9999.0;
        double corex = -9999.0;
        double corey = -9999.0;
        double core_err = -9999.0;
        double xmax = -9999.0;
        double xmax_err = -9999.0;
        double energy = -9999.0;
        double energy_err = -9999.0;
    
        double mscw = -9999.;
        double mscw_err = -9999.;
        double mscl = -9999.;
        double mscl_err = -9999.;
    
        double i_UTC;
    
        rec2.mapColumnToVar( "EVENT_ID", event_id );
        rec2.mapColumnToVar( "FLAGS", flags );
        rec2.mapColumnToVar( "TIME", ttime );
        rec2.mapColumnToVar( "MULTIP", multip );
        rec2.mapColumnToVar( "TELMASK", telmask );
        rec2.mapColumnToVar( "RA", ra );
        rec2.mapColumnToVar( "DEC", dec );
        rec2.mapColumnToVar( "DIR_ERR", dir_err );
        rec2.mapColumnToVar( "DETX", detx );
        rec2.mapColumnToVar( "DETY", dety );
        rec2.mapColumnToVar( "ALT_PNT", alt_pnt );
        rec2.mapColumnToVar( "AZ_PNT", az_pnt );
        rec2.mapColumnToVar( "ALT", alt );
        rec2.mapColumnToVar( "AZ", az );
        rec2.mapColumnToVar( "COREX", corex );
        rec2.mapColumnToVar( "COREY", corey );
        rec2.mapColumnToVar( "CORE_ERR", core_err );
        rec2.mapColumnToVar( "XMAX", xmax );
        rec2.mapColumnToVar( "XMAX_ERR", xmax_err );
        rec2.mapColumnToVar( "ENERGY", energy );
        rec2.mapColumnToVar( "ENERGY_ERR", energy_err );
        rec2.mapColumnToVar( "HIL_MSW", mscw );
        rec2.mapColumnToVar( "HIL_MSW_ERR", mscw_err );
        rec2.mapColumnToVar( "HIL_MSL", mscl );
        rec2.mapColumnToVar( "HIL_MSL_ERR", mscl_err );
    
        double start = -9999.;
        double stop  =  9999.;
    
        rec3.mapColumnToVar( "START", start );
        rec3.mapColumnToVar( "STOP", stop );
    
        //strncpy( telclass, "VERITAS", NAMELEN );
    
    //    for( unsigned int i = 0; i < fTreeDet->GetEntries(); i++ )
    
    //////////////////////////////////////////////////////////
        //FITSRecord hdu( iOutputFile );
        FITSRecord rec2( iOutputFile, templ, "EVENTS" );
        FITSRecord rec3( iOutputFile, templ, "GTI" );
        TChain *chain = new TChain("data");
        chain->Add( iInputFile.c_str() );
        CData dchain( chain, true, 6, false );
    
        double decDiff = 0.;
        double raDiff = 0.;
        double tra = 83.6331;
        double tdec = 22.0145;
        double wobbleN = -99.;
        double wobbleE = -99.;
    
        string creator = "Hughes";
        //string telescop = "VERITAS";
        string telescop = "CTA";
        string date_obs = "1980-10-15";
        string time_obs = "00:00:00";
        string date_end = "1980-10-15";
        string time_end = "00:00:00";
        string timeunit = "s";
        string timesys = "TT";
        string timeref = "local";
        string tassign = "AZ";
        string object = "Crab";
        string radesys = "FK5";
        string obs_mode = "on";
    
        rec2.writeHeader( "CREATOR", creator );
    	cout << "111" << endl;
        rec2.writeHeader( "TELESCOP", telescop );
        rec2.writeHeader( "OBS_ID", 0 );
        rec2.writeHeader( "DATE_OBS", date_obs );
        rec2.writeHeader( "TIME_OBS", time_obs );
        rec2.writeHeader( "DATE_END", date_end );
        rec2.writeHeader( "TIME_END", time_end  );
        rec2.writeHeader( "TSTART", 0 );
        rec2.writeHeader( "TSTOP", 0 );
        rec2.writeHeader( "MJDREFI", 51910 );
        rec2.writeHeader( "MJDREFF", 7.4287037e-4 );
        rec2.writeHeader( "TIMEUNIT", timeunit );
        rec2.writeHeader( "TIMESYS", timesys );
        rec2.writeHeader( "TIMEREF", timeref );
        rec2.writeHeader( "TASSIGN", tassign );
        rec2.writeHeader( "TELAPSE", 0 );
        rec2.writeHeader( "ONTIME", 1200 );
        rec2.writeHeader( "LIVETIME", 1200 );
        rec2.writeHeader( "DEADC", 0 );
        rec2.writeHeader( "OBJECT", object );
        rec2.writeHeader( "RA_OBJ", tra );
        rec2.writeHeader( "DEC_OBJ",tdec );
        rec2.writeHeader( "RA_PNT", tra );
        rec2.writeHeader( "DEC_PNT", tdec );
        rec2.writeHeader( "ALT_PNT", 0 );
        rec2.writeHeader( "AZ_PNT", 0 );
        rec2.writeHeader( "RADECSYS", radesys );
        rec2.writeHeader( "EQUINOX", 2000 );
        rec2.writeHeader( "CONV_DEP", 0 );
        rec2.writeHeader( "CONV_RA", tra );
        rec2.writeHeader( "CONV_DEC", tdec );
        rec2.writeHeader( "OBS_MODE", obs_mode );
        rec2.writeHeader( "N_TELS", 0 );
        rec2.writeHeader( "TELLIST", 0 );
        rec2.writeHeader( "GEOLAT", 0 );
        rec2.writeHeader( "GEOLON", 0 );
        rec2.writeHeader( "ALTITUDE", 0 );
    
        for( unsigned int i = 0; i <  chain->GetEntries(); i++ )
        {
            chain->GetEntry(i);
    
    	event_id = dchain.runNumber;
    	flags    = false;
    	multip   = dchain.NImages;
    	telmask  = true;
    
            ttime    = dchain.Time;
    	detx     = dchain.Xoff_derot;
    	dety     = dchain.Yoff_derot;
    //	wobbleN  = dchain.WobbleN;
    //	wobbleE  = dchain.WobbleE;
    //XXXX
    	wobbleN  = 0.;
    	wobbleE  = 0.;
    
    	ra  = tra  + detx;
    	dec = tdec + dety;
    	ra  = detx;
    	dec = dety;
    	ra  += wobbleE;
    	dec += wobbleN;
    
    	ra  += tra;
    	dec += tdec;
    
    
    	alt_pnt  = 90. - dchain.TelElevation[0];
    	az_pnt    = dchain.TelAzimuth[0];
    	alt      = 90. - dchain.Ze;
    	az       = dchain.Az;
    
    	corex    = dchain.Xcore;
    	corey    = dchain.Ycore;
    //	core_err = dchain.;
    	xmax     = dchain.EmissionHeight;
    	xmax_err = dchain.EmissionHeightChi2;
    	energy   = dchain.Erec;
    	energy_err = dchain.dE;
    	mscw       = dchain.MSCW;
    //	mscw_err   = dchain.;
    	mscl       = dchain.MSCL;
    //	mscl_err   = dchain.;
    
    //        if( mscw > -1.2 && mscw < 0.35)
    //          if( mscl > -1.2 && mscl < 0.7)
    //	    if( energy > 0. ) rec2.write();
    
        }
    */
    //    rec2.write();
    //    rec3.write();
    
    //    rec2.finishWriting();
    //    rec3.finishWriting();
    
    mscwfile->Close();
    
    return 0;
    
}

int main( int argc, char* argv[] )
{

    if( argc < 2 )
    {
        cout << endl;
        cout << "./writeCTA_EVENTLIST <input file> <output name>" << endl;
        cout << endl;
        exit( -1 );
    }
    
    string iInputFile  = argv[1];
    string iOutputFile = argv[2];
    
    writeTELARRAY( iInputFile, iOutputFile );
    
    return 0;
    
}

