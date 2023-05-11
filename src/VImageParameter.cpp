/*! \class VImageParameter
    \brief class to store run and image parameter information (single telescope)

    variables called VARNAME_SC are in shower coordinates

*/

#include "VImageParameter.h"

VImageParameter::VImageParameter(
    unsigned int iShortTree,
    bool iWriteNImagePixels )
{
    fMC = false;
    tpars = 0;
    fShortTree = iShortTree;
    fWriteNImagePixels = iWriteNImagePixels;
    
    reset();
}


VImageParameter::~VImageParameter()
{
    delete tpars;
}


void VImageParameter::initTree( string iName, string iTitle, bool iMC, bool iLL, bool iMuon, bool iHough )
{

    tpars = new TTree( iName.c_str(), iTitle.c_str() );
    tpars->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    tpars->SetAutoSave( 150000000 );               // autosave when 150 Mbytes written
    
    if( fShortTree < 1 )
    {
        tpars->Branch( "telID", &fTelID, "telID/I" );
        tpars->Branch( "runNumber", &runNumber, "runNumber/I" );
        tpars->Branch( "MJD",  &MJD,  "MJD/I" );
        tpars->Branch( "Time",  &time,  "time/D" );
    }
    tpars->Branch( "eventType", &eventType, "eventType/s" );
    // tpars tree have same number of events as showerpars tree:
    // showerpars->AddFriend( "Tel_3/tpars", "3224.v3.root" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "eventNumber",  &eventNumber,  "eventNumber/I" );
    }
    tpars->Branch( "eventStatus", &eventStatus, "eventStatus/i" );
    
    if( fShortTree < 1 )
    {
        tpars->Branch( "fsumfirst", &fsumfirst, "fsumfirst/I" );
        tpars->Branch( "fsumwindow", &fsumwindow, "fsumwindow/I" );
        tpars->Branch( "fsumwindow_2", &fsumwindow_2, "fsumwindow_2/I" );
    }
    
    // MC parameters
    if( iMC && fShortTree < 1 )
    {
        tpars->Branch( "MCprim", &MCprimary, "MCprimary/s" );
        tpars->Branch( "MCe0", &MCenergy, "MCenergy/F" );
        tpars->Branch( "MCxcore", &MCxcore, "MCxcore/F" );
        tpars->Branch( "MCycore", &MCycore, "MCycore/F" );
        tpars->Branch( "MCxcos", &MCxcos, "MCxcos/F" );
        tpars->Branch( "MCycos", &MCycos, "MCycos/F" );
        tpars->Branch( "MCLocalTriggerTime", &MCLocalTriggerTime, "MCLocalTriggerTime/F" );
        tpars->Branch( "MCLocalDelayedTriggerTime", &MCLocalDelayedTriggerTime, "MCLocalDelayedTriggerTime/F" );
        tpars->Branch( "MCTel_Xoff", &MCTel_Xoff, "MCTel_Xoff/F" );
        tpars->Branch( "MCTel_Yoff", &MCTel_Yoff, "MCTel_Yoff/F" );
    }
    
    // image parameters
    if( fShortTree < 1 )
    {
        tpars->Branch( "meanPed_Image", &fmeanPed_Image, "meanPed_Image/F" );
    }
    tpars->Branch( "meanPedvar_Image", &fmeanPedvar_Image, "meanPedvar_Image/F" );
    tpars->Branch( "cen_x", &cen_x, "cen_x/F" );
    tpars->Branch( "cen_y", &cen_y, "cen_y/F" );
    tpars->Branch( "f_d", &f_d, "f_d/F" );
    tpars->Branch( "f_s", &f_s, "f_s/F" );
    tpars->Branch( "f_sdevxy", &f_sdevxy, "f_sdevxy/F" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "cen_x_trig", &cen_x_trig, "cen_x_trig/F" );
        tpars->Branch( "cen_y_trig", &cen_y_trig, "cen_y_trig/F" );
        tpars->Branch( "cen_x2_trig", &cen_x2_trig, "cen_x2_trig/F" );
        tpars->Branch( "cen_y2_trig", &cen_y2_trig, "cen_y2_trig/F" );
    }
    tpars->Branch( "length", &length, "length/F" );
    tpars->Branch( "width", &width, "width/F" );
    tpars->Branch( "size", &size, "size/F" );
    tpars->Branch( "size2", &size2, "size2/F" );
    tpars->Branch( "loss", &loss, "loss/F" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "sizeLL", &sizeLL, "sizeLL/F" );
        tpars->Branch( "size2LL", &size2LL, "size2LL/F" );
        tpars->Branch( "lossAndDead", &lossAndDead, "lossAndDead/F" );
    }
    tpars->Branch( "fui", &fui, "fui/F" );
    tpars->Branch( "fracLow", &fracLow, "fracLow/F" );
    tpars->Branch( "dist", &dist, "dist/F" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "azwidth", &azwidth, "azwidth/F" );
        tpars->Branch( "alpha", &alpha, "alpha/F" );
        tpars->Branch( "los", &los, "los/F" );
        tpars->Branch( "miss", &miss, "miss/F" );
        tpars->Branch( "phi", &phi, "phi/F" );
    }
    tpars->Branch( "cosphi", &cosphi, "cosphi/F" );
    tpars->Branch( "sinphi", &sinphi, "sinphi/F" );
    
    tpars->Branch( "ntubes", &ntubes, "ntubes/s" );
    // contains number of triggered pixels
    tpars->Branch( "trig_tubes", &trig_tubes, "trig_tubes/s" );
    tpars->Branch( "nzerosuppressed", &nzerosuppressed, "nzerosuppressed/s" );
    tpars->Branch( "nsat", &nsat, "nsat/s" );
    tpars->Branch( "nlowgain", &nlowgain, "nlowgain/s" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "ntubesBNI", &ntubesBrightNoImage, "ntubesBNI/s" );
        tpars->Branch( "max", &max, "max[3]/F" );
        tpars->Branch( "index_of_max", &index_of_max, "index_of_max[3]/s" );
    }
    tpars->Branch( "asymmetry", &asymmetry, "asymmetry/F" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "frac", &frac, "frac[3]/F" );
        tpars->Branch( "ntrig", &ntrig, "ntrig/s" );
        tpars->Branch( "ntrig_per_patch", &ntrig_per_patch, "ntrig_per_patch/s" );
    }
    tpars->Branch( "bad", &bad, "bad/s" );
    tpars->Branch( "badLow", &badLow, "badLow/s" );
    tpars->Branch( "tgrad_x", &tgrad_x, "tgrad_x/F" );
    if( fShortTree < 1 )
    {
        tpars->Branch( "tint_x", &tint_x, "tint_x/F" );
        tpars->Branch( "tgrad_dx", &tgrad_dx, "tgrad_dx/F" );
        tpars->Branch( "tint_dx", &tint_dx, "tint_dx/F" );
        tpars->Branch( "tchisq_x", &tchisq_x, "tchisq_x/F" );
    }
    
    // muon parameters (Iterative fit muon analysis)
    if( iMuon || iHough )
    {
        tpars->Branch( "muonX0", &muonX0, "muonX0/F" );
        tpars->Branch( "muonY0", &muonY0, "muonY0/F" );
        tpars->Branch( "muonXC", &muonXC, "muonXC/F" );
        tpars->Branch( "muonYC", &muonYC, "muonYC/F" );
        tpars->Branch( "muonRadius", &muonRadius, "muonRadius/F" );
        tpars->Branch( "muonRSigma", &muonRSigma, "muonRSigma/F" );
        tpars->Branch( "muonSize", &muonSize, "muonSize/F" );
        tpars->Branch( "muonIPCorrectedSize", &muonIPCorrectedSize, "muonIPCorrectedSize/F" );
        tpars->Branch( "muonValid", &muonValid, "muonValid/I" );
    }
    // muon parameters (Hough transform)
    if( iHough )
    {
    
        tpars->Branch( "houghAP", &houghAP, "houghAP/D" );
        tpars->Branch( "houghTD", &houghTD, "houghTD/D" );
        tpars->Branch( "houghNpix", &houghNpix, "houghNpix/I" );
        tpars->Branch( "houghCN", &houghCN, "houghCN/D" );
        tpars->Branch( "houghContained", &houghContained, "houghContained/D" );
        tpars->Branch( "houghMuonValid", &houghMuonValid, "houghMuonValid/I" );
        
    }
    
    // log likelihood fit errors/results
    tpars->Branch( "Fitstat", &Fitstat, "Fitstat/I" );
    if( iLL )
    {
        // number of dead tube values estimated by fit
        tpars->Branch( "ntfit", &ntfit, "ntfit/I" );
        tpars->Branch( "Fitmin", &Fitmin, "Fitmin/F" );
        tpars->Branch( "Fitedm", &Fitedm, "Fitedm/F" );
        // number of recovered dead channels
        tpars->Branch( "ntRec", &ntRec, "ntRec/I" );
        tpars->Branch( "ddist", &ddist, "ddist/F" );
        tpars->Branch( "dmiss", &dmiss, "dmiss/F" );
        tpars->Branch( "dalpha", &dalpha, "dalpha/F" );
        tpars->Branch( "dazwidth", &dazwidth, "dazwidth/F" );
        tpars->Branch( "rho", &rho, "rho/F" );
        tpars->Branch( "drho", &drho, "drho/F" );
        tpars->Branch( "sigmaX", &sigmaX, "sigmaX/F" );
        tpars->Branch( "dsigmaX", &dsigmaX, "dsigmaX/F" );
        tpars->Branch( "sigmaY", &sigmaY, "sigmaY/F" );
        tpars->Branch( "dsigmaY", &dsigmaY, "dsigmaY/F" );
        tpars->Branch( "signal", &signal, "signal/F" );
        tpars->Branch( "dsignal", &dsignal, "dsignal/F" );
    }
    tpars->Branch( "dcen_x", &dcen_x, "dcen_x/F" );
    tpars->Branch( "dcen_y", &dcen_y, "dcen_y/F" );
    tpars->Branch( "dlength", &dlength, "dlength/F" );
    tpars->Branch( "dwidth", &dwidth, "dwidth/F" );
    tpars->Branch( "dphi", &dphi, "dphi/F" );
    
    // image / border pixel list
    if( fWriteNImagePixels )
    {
        // image pixels
        tpars->Branch( "PixelListN", &PixelListN, "PixelListN/i" );
        tpars->Branch( "PixelID", PixelID, "PixelID[PixelListN]/i" );
        tpars->Branch( "PixelType", PixelType, "PixelType[PixelListN]/i" );
        tpars->Branch( "PixelIntensity", PixelIntensity, "PixelIntensity[PixelListN]/F" );
        tpars->Branch( "PixelTimingT0", PixelTimingT0, "PixelTimingT0[PixelListN]/F" );
        tpars->Branch( "PixelPE", PixelPE, "PixelPE[PixelListN]/F" );
    }
}


void VImageParameter::reset( unsigned int resetLevel )
{
    if( resetLevel == 0 )
    {
        fTelID = 0;
        runNumber = 0;
        eventNumber = 0;
        eventType = 0;
        time = 0.;
        MJD = 0;
        eventStatus = 0;
        
        // telescope parameters
        fsumfirst = 8;
        fsumwindow = 15;
        fsumwindow_2 = 5;
        fLocalTrigger = 0;
        
        // telescope positions in shower coordinates
        Tel_x_SC = 0.;
        Tel_y_SC = 0.;
        Tel_z_SC = 0.;
    }
    
    // MC parameters
    MCprimary = 0;
    MCenergy = 0.;
    MCxcore = 0.;
    MCycore = 0.;
    MCxcos = 0.;
    MCycos = 0.;
    MCLocalTriggerTime = 0.;
    MCLocalDelayedTriggerTime = 0.;
    MCTel_Xoff = 0.;
    MCTel_Yoff = 0.;
    //image parameters
    fTrig_type = 0;
    f_d = 0.;
    f_s = 0.;
    f_sdevxy = 0.;
    fmeanPed_Image = 0.;
    fmeanPedvar_Image = 0.;
    cen_x = 0.;
    cen_y = 0.;
    cen_x_trig = 0.;
    cen_y_trig = 0.;
    cen_x2_trig = 0.;
    cen_y2_trig = 0.;
    sigmaX = 0.;
    sigmaY = 0.;
    length = 0.;
    width = 0.;
    size = 0.;
    size2 = 0.;
    sizeLL = 0.;
    size2LL = 0.;
    loss = 0.;
    lossAndDead = 0.;
    fui = 0.;
    fracLow = 0.;
    dist = 0.;
    azwidth = 0.;
    alpha = 0.;
    los = 0.;
    miss = 0.;
    phi = 0.;
    sinphi = 0.;
    cosphi = 0.;
    asymmetry = 0.;
    
    nsamples = 0;
    ntubes = 0;
    ntubesBrightNoImage = 0;
    trig_tubes = 0;
    ntrig = 0;
    ntrig_per_patch = 0;
    nsat = 0;
    nlowgain = 0;
    nzerosuppressed = 0;
    bad = 0;
    badLow = 0;
    for( int i = 0; i < 3; i++ )
    {
        max[i] = 0.;
    }
    for( int i = 0; i < 3; i++ )
    {
        index_of_max[i] = 0;
    }
    for( int i = 0; i < 3; i++ )
    {
        frac[i] = 0.;
    }
    
    // timing parameters
    tint_x = 0.;
    tgrad_x = 0.;
    tint_dx = 0.;
    tgrad_dx = 0.;
    tchisq_x = 0.;
    tmin = 999.;
    tmax = 0.;
    tmean = 0.;
    
    // log likelihood fit parameters and errors
    ntfit = 0;
    Fitmin = 0.;
    Fitedm = 0.;
    Fitstat = -1;
    
    dcen_x = 0.;
    dcen_y = 0.;
    dsigmaX = 0.;
    dsigmaY = 0.;
    dlength = 0.;
    dwidth = 0.;
    ddist = 0.;
    dmiss = 0.;
    dphi = 0.;
    dalpha = 0.;
    dazwidth = 0.;
    rho = 0.;
    drho = 0.;
    signal = 0.;
    dsignal = 0.;
    
    // Iterative fit muon analysis
    muonX0 = 0.;
    muonY0 = 0.;
    muonXC = 0.;
    muonYC = 0.;
    muonRadius = 0.;
    muonRSigma = 0.;
    muonSize = 0.;
    muonIPCorrectedSize = 0.;
    muonValid = 0;
    
    // Hough transform muon parameters
    houghAP = 0.;
    houghTD = 0.;
    houghNpix = 0;
    houghCN = 0.;
    houghContained = 0.;
    houghMuonValid = 0;
    
    if( fWriteNImagePixels )
    {
        PixelListN = 0;
        memset( PixelID, 0, VDST_MAXCHANNELS * sizeof( unsigned int ) );
        memset( PixelType, 0, VDST_MAXCHANNELS * sizeof( unsigned int ) );
        memset( PixelIntensity, 0, VDST_MAXCHANNELS * sizeof( float ) );
        memset( PixelTimingT0, 0, VDST_MAXCHANNELS * sizeof( float ) );
        memset( PixelPE, 0, VDST_MAXCHANNELS * sizeof( float ) );
    }
}


void VImageParameter::printParameters()
{
    cout << "Image: " << fTelID << "\t" << cen_x << "\t" << cen_y << "\t" << dist << "\t" << length << "\t" << width << "\t" << alpha << "\t" << size << endl;
    cout << "Telescope: " << fTelID << "\t" << Tel_x_SC << "\t" << Tel_y_SC << "\t" << Tel_z_SC << endl;
}


bool VImageParameter::hasImage()
{
    if( ntubes > 0 )
    {
        return true;
    }
    
    return false;
}


void VImageParameter::fill()
{
    tpars->Fill();
}


void VImageParameter::setImageBorderPixelPosition( vector< float > iImageBorderPixelPosition_x, vector< float > iImageBorderPixelPosition_y )
{
    fImageBorderPixelPosition_x = iImageBorderPixelPosition_x;
    fImageBorderPixelPosition_y = iImageBorderPixelPosition_y;
}
