/*! \class VRadialAcceptance
 *  \brief fill radial acceptance histograms
 *
 *
 */

#include "VRadialAcceptance.h"

/*!
 *   use acceptance curve from simulation
 */
VRadialAcceptance::VRadialAcceptance()
{
    reset();
    
    fAcceptanceFunctionDefined = false;
}


/*!
 *  ******************************************************************
 *          reading acceptance curves
 *  ******************************************************************
 *  this constructor is called from anasum to get acceptance
 *
 */
VRadialAcceptance::VRadialAcceptance( string ifile, int irun )
{
    reset();

    // ignore acceptance files (getacceptance will always return 1)
    if( ifile == "IGNOREACCEPTANCE" || ifile == "simu" )
    {
        fAcceptanceFunctionDefined = false;
        return;
    }

    if( ifile.find( "/" ) == string::npos )
    {
        ifile = VUtilities::testFileLocation( ifile, "RadialAcceptances/", true );
    }
    fAccFile = new TFile( ifile.c_str() );
    if( fAccFile->IsZombie() )
    {
        cout << "VRadialAcceptance::VRadialAcceptance error reading acceptance file " << ifile << endl;
        exit( EXIT_FAILURE );
    }
    // get fit to 1D radial acceptance function
    char hname[200];
    if( irun < 0 )
    {
        sprintf( hname, "fRadialAcceptance" );
    }
    else
    {
        sprintf( hname, "fRadialAcceptance_perRun_%d", irun );
    }
    fRadialAcceptanceFit = ( TF1* )fAccFile->Get( hname );
    if( !fRadialAcceptanceFit )
    {
        cout << "VRadialAcceptance::VRadialAcceptance no radial acceptance fit function (";
        cout << hname << ") found in ";
        cout << ifile << endl;
        exit( EXIT_FAILURE );
    }
    fAcceptanceFunctionDefined = true;
    cout << "Reading radial acceptance function";
    if( irun > 0 )
    {
        cout << " for run " << irun;
    }
    cout << " from " << ifile;
    
    // count number of raw files used to calculate acceptances
    if( irun < 0 )
    {
        TKey* key;
        TIter nextkey( fAccFile->GetListOfKeys() );
        while( ( key = ( TKey* )nextkey() ) )
        {
            TObject* obj = key->ReadObj();
            string itemp = obj->GetName();
            if( itemp.find( "hRadialAcceptance_perRun_" ) < itemp.size()
             && itemp.find( "Fit" ) == string::npos  )
            {
                fNumberOfRawFiles++;
            }
        }
    }
    else
    {
        fNumberOfRawFiles = 1;
    }
    cout << ", calculated from " << fNumberOfRawFiles << " file(s)" << endl;
}


/*!

 *  ******************************************************************
 *          calculating (filling) acceptance curves
 *  ******************************************************************
 *   this constructor is called for determination of radial acceptance curves with makeAcceptance

 */
VRadialAcceptance::VRadialAcceptance( VGammaHadronCuts* icuts, VAnaSumRunParameter* irunpar, double iMaxDistanceAllowed )
{
    reset();
    
    fRunPar = irunpar;
    if( !fRunPar )
    {
        cout << "VRadialAcceptance error: no runparameter defined" << endl;
        cout << "exiting..";
        exit( EXIT_FAILURE );
    }
    
    fCuts = icuts;
    if( !fCuts )
    {
        cout << "VRadialAcceptance error: no gamma/hadron separation cuts defined" << endl;
        cout << "exiting..";
        exit( EXIT_FAILURE );
    }
    // maximum distance to camera center for which events are taken into account:
    if( iMaxDistanceAllowed > 0. )
    {
        fCut_CameraFiducialSize_max = iMaxDistanceAllowed;
        fCuts->fCut_CameraFiducialSize_max = iMaxDistanceAllowed;
        
    }
    else
    {
        fCut_CameraFiducialSize_max = fCuts->fCut_CameraFiducialSize_max;
    }
    
    
    // maximal offset [deg]
    double xymax = 5.0;
    int nxybin = 50;
    
    //////////////////////////////////////////////////////////////////
    // range used to normalise acceptance histograms (in bins)
    fAccZeFitMinBin = 2;
    fAccZeFitMaxBin = 5;
    
    //////////////////////////////////////////////////////////////////
    // list of histograms
    hList = new TList();
    // list of histograms which are written in production mode to file
    histListProductionIO = new TList();
    // list of histograms which will be normalized by area (1D)
    hListNormalizeHistograms = new TList();
    // list of histograms for which the 1D radial acceptance fit is applied
    hListFitHistograms = new TList();
    
    //////////////////////////////////////////////////////////////////
    // histogram used for the normalisation of the 1D acceptances
    hscale = new TH1F( "hscale", "", nxybin, 0., xymax );
    hscale->SetStats( 0 );
    hscale->SetXTitle( "distance to camera centre [deg]" );
    hscale->SetYTitle( "area / bin" );
    for( int i = 1; i < nxybin; i++ )
    {
        hscale->SetBinContent( i, TMath::Pi()*xymax * xymax / ( ( double )( nxybin * nxybin ) ) * ( 2 * i - 1 ) );
    }
    for( int i = 1; i < nxybin; i++ )
    {
        hscale->SetBinError( i, 0. );
    }
    hList->Add( hscale );
    
    //////////////////////////////////////////////////////
    // 1D radial acceptance
    // (main result)
    fRadialAcceptance = new TH1F( "hRadialAcceptance", "radial acceptance (average over all runs)", nxybin, 0., xymax );
    fRadialAcceptance->SetXTitle( "distance to camera center [deg]" );
    fRadialAcceptance->SetYTitle( "relative rate" );
    fRadialAcceptance->SetMarkerSize( 2 );
    fRadialAcceptance->SetLineWidth( 2 );
    fRadialAcceptance->Sumw2();
    hList->Add( fRadialAcceptance );
    hListNormalizeHistograms->Add( fRadialAcceptance );
    hListFitHistograms->Add( fRadialAcceptance );
    histListProductionIO->Add( fRadialAcceptance );
    
    char hname[200];
    char htitle[200];
    //////////////////////////////////////////////////////
    // run dependent scaling and acceptance curves histograms
    // (one histogram per run)
    for( unsigned int i = 0; i < fRunPar->fRunList.size(); i++ )
    {
        sprintf( hname, "hscaleRun_%d", fRunPar->fRunList[i].fRunOff );
        sprintf( htitle, "run %d", fRunPar->fRunList[i].fRunOff );
        hscaleRun.push_back( new TH1F( hname, htitle, nxybin, 0., xymax ) );
        hscaleRun.back()->SetXTitle( "distance to camera center [deg]" );
        hscaleRun.back()->SetYTitle( "number of counts" );
        hscaleRun.back()->SetMarkerSize( 2 );
        hscaleRun.back()->SetLineWidth( 2 );
        hscaleRun.back()->Sumw2();
        hList->Add( hscaleRun.back() );
        
        sprintf( hname, "hscaleRunRatio_%d", fRunPar->fRunList[i].fRunOff );
        sprintf( htitle, "run %d", fRunPar->fRunList[i].fRunOff );
        hscaleRunRatio.push_back( new TH1F( *hscale ) );
        hscaleRunRatio.back()->SetNameTitle( hname, htitle );
        hscaleRunRatio.back()->SetXTitle( "distance to camera center [deg]" );
        hscaleRunRatio.back()->SetYTitle( "ratio without/with exclusion regions" );
        hscaleRunRatio.back()->SetMarkerSize( 2 );
        hscaleRunRatio.back()->SetLineWidth( 2 );
        hList->Add( hscaleRunRatio.back() );
        
        // (large bin numbers required to achieve sufficent accuracy)
        sprintf( hname, "hAreaExcluded2D_%d", fRunPar->fRunList[i].fRunOff );
        sprintf( htitle, "run %d", fRunPar->fRunList[i].fRunOff );
        hAreaExcluded2D.push_back( new TH2F( hname, htitle, 1000, -xymax, xymax, 1000, -xymax, xymax ) );
        hAreaExcluded2D.back()->SetXTitle( "x_{off,derot} [deg]" );
        hAreaExcluded2D.back()->SetYTitle( "y_{off,derot} [deg]" );
        hAreaExcluded2D.back()->Sumw2();
        hList->Add( hAreaExcluded2D.back() );
        
        // radial acceptances
        sprintf( hname, "hRadialAcceptance_perRun_%d", fRunPar->fRunList[i].fRunOff );
        sprintf( htitle, "radial acceptance (run %d)", fRunPar->fRunList[i].fRunOff );
        fRadialAcceptance_perRun.push_back( new TH1F( hname, htitle, nxybin, 0., xymax ) );
        fRadialAcceptance_perRun.back()->SetXTitle( "distance to camera center [deg]" );
        fRadialAcceptance_perRun.back()->SetYTitle( "relative rate" );
        fRadialAcceptance_perRun.back()->SetMarkerSize( 2 );
        fRadialAcceptance_perRun.back()->SetLineWidth( 2 );
        fRadialAcceptance_perRun.back()->SetMarkerColor( 2 );
        fRadialAcceptance_perRun.back()->SetLineColor( 2 );
        fRadialAcceptance_perRun.back()->Sumw2();
        hList->Add( fRadialAcceptance_perRun.back() );
        hListNormalizeHistograms->Add( fRadialAcceptance_perRun.back() );
        histListProductionIO->Add( fRadialAcceptance_perRun.back() );
        hListFitHistograms->Add( fRadialAcceptance_perRun.back() );
        
        sprintf( hname, "hXYAccRun_%d", fRunPar->fRunList[i].fRunOff );
        sprintf( htitle, "run %d", fRunPar->fRunList[i].fRunOff );
        hXYAccRun.push_back( new TH2F( hname, htitle, 40, -2., 2., 40, -2., 2. ) );
        hXYAccRun.back()->SetXTitle( "x_{off} [deg]" );
        hXYAccRun.back()->SetYTitle( "y_{off} [deg]" );
        hList->Add( hXYAccRun.back() );
    }
    
    
    // phi angle (camera coordinates)
    hPhiDist = new TH1F( "hPhiDist", "", nxybin, -180., 180. );
    hPhiDist->SetXTitle( "azimuth (camera coordinates) [deg]" );
    hList->Add( hPhiDist );
    
    // 2D histogram (not normalized)
    
    double hr = 3.0 ;
    int hn = 40 ;
    hXYAccTotDeRot = new TH2F( "hXYAccTotDeRot", "",  hn, -hr, hr, hn, -hr, hr );
    hXYAccTotDeRot->SetXTitle( "x_{off,derot} [deg]" );
    hXYAccTotDeRot->SetYTitle( "y_{off,derot} [deg]" );
    hXYAccTotDeRot->Sumw2();
    hList->Add( hXYAccTotDeRot );
    histListProductionIO->Add( hXYAccTotDeRot );
    
    // azimuth dependent radial acceptance histograms
    fPhiMin.clear();
    fPhiMax.clear();
    /*	fPhiMin.push_back( 135.0 );
    	fPhiMax.push_back( -165.0 );
    	fPhiMin.push_back( 150.0 );
    	fPhiMax.push_back( -150.0 );
    	fPhiMin.push_back( -180. );
    	fPhiMax.push_back( -120. );
    	for( int i = 0; i < 13; i++ )
    	{
    		fPhiMin.push_back( fPhiMin.back() + 22.5 );
    		fPhiMax.push_back( fPhiMax.back() + 22.5 );
    	} */
    for( unsigned int i = 0; i < fPhiMin.size(); i++ )
    {
        // camera coordinates
        sprintf( hname, "hAccPhi_%d", i );
        sprintf( htitle, "%.0f < Phi < %.0f", fPhiMin[i], fPhiMax[i] );
        hAccPhi.push_back( new TH1F( hname, htitle, nxybin, 0., xymax ) );
        hAccPhi.back()->SetXTitle( "distance to camera center [deg]" );
        hAccPhi.back()->SetYTitle( "relative rate" );
        hAccPhi.back()->SetMarkerSize( 2 );
        hAccPhi.back()->SetLineWidth( 2 );
        hAccPhi.back()->Sumw2();
        hList->Add( hAccPhi.back() );
        hListNormalizeHistograms->Add( hAccPhi.back() );
        hListFitHistograms->Add( hAccPhi.back() );
        // derotated camera coordinates
        sprintf( hname, "hAccPhiDerot_%d", i );
        sprintf( htitle, "%.0f < Phi < %.0f (derot)", fPhiMin[i], fPhiMax[i] );
        hAccPhiDerot.push_back( new TH1F( hname, htitle, nxybin, 0., xymax ) );
        hAccPhiDerot.back()->SetXTitle( "distance to camera center [deg]" );
        hAccPhiDerot.back()->SetYTitle( "relative rate" );
        hAccPhiDerot.back()->SetMarkerSize( 2 );
        hAccPhiDerot.back()->SetLineWidth( 2 );
        hAccPhiDerot.back()->Sumw2();
        hList->Add( hAccPhiDerot.back() );
        hListNormalizeHistograms->Add( hAccPhiDerot.back() );
        hListFitHistograms->Add( hAccPhiDerot.back() );
    }
    
    
    // for PhiDependentSlice hists
    phi_minphi = 0.0 ;
    phi_maxphi = 2 * TMath::Pi() ;
    phi_nbins  = 18 ;
    phi_minradius = 1.25 ; // band will be from 1.25 deg to 1.75 deg
    phi_maxradius = 1.75 ;
    
    // for RadiusDependentSlice hists
    rad_minrad = 0.0 ;
    rad_maxrad = 2.0 ;
    rad_nbins  = 12  ;
    rad_phiwidth = TMath::Pi() / 6 ; // Slice000 will be from 330 deg to 30 deg, etc
    int centerdeg = 0;
    
    // global 1D Slice Histograms
    sprintf( hname, "hXYAccTotDeRotPhiDependentSlice" ) ;
    sprintf( htitle, "1D histogram from getXoff_derot() and getYoff_derot(), with All Events, PhiDependentSlice (doughnut)\n, %d Bins from %3.1f to %3.1f radians, from Radius %4.2f to %4.2f deg", phi_nbins, phi_minphi, phi_maxphi, phi_minradius, phi_maxradius ) ;
    hXYAccTotDeRotPhiDependentSlice = new TH1F( hname, htitle, phi_nbins, phi_minphi, phi_maxphi ) ;
    
    sprintf( hname, "hXYAccTotDeRotRadiusDependentSlice" ) ;
    sprintf( htitle, "1D histogram from getXoff_derot() and getYoff_derot(), with All Events, RadiusDependentSlice%03d (Pie Slice)\n, from Phi %3.1f to %3.1f radians", centerdeg, centerdeg - ( rad_phiwidth * 180 / 3.1415 ), centerdeg + ( rad_phiwidth * 180 / 3.1415 ) ) ;
    hXYAccTotDeRotRadiusDependentSlice000 = new TH1F( hname, htitle, rad_nbins, rad_minrad, rad_maxrad ) ;
    
    // 2D histograms, sorted by ImgSel
    // (max number of telescope == 4)
    TH2F* tmphist2 ;
    TH1F* tmphist1 ;
    for( int i = 0 ; i <= 15 ; i++ )
    {
        sprintf( hname, "hAccImgSel_%d", i ) ;
        sprintf( htitle, "2D histogram of getXoff_derot() and getYoff_derot(), with ImgSel=%d", i ) ;
        tmphist2 = new TH2F( hname, htitle, hn, -1.0 * hr, hr, hn, -1.0 * hr, hr ) ;
        hXYAccImgSel.push_back( tmphist2 ) ;
        
        sprintf( hname, "hAccImgSelPreDeRot_%d", i ) ;
        sprintf( htitle, "2D histogram of Xoff and Yoff, with ImgSel=%d", i ) ;
        tmphist2 = new TH2F( hname, htitle, hn, -1.0 * hr, hr, hn, -1.0 * hr, hr ) ;
        hXYAccImgSelPreDeRot.push_back( tmphist2 ) ;
        
        sprintf( hname, "hAccImgSelPhiDependentSlice_%d", i ) ;
        sprintf( htitle, "1D histogram from getXoff_derot() and getYoff_derot(), with ImgSel=%d, PhiDependentSlice\n, %d Bins from %3.1f to %3.1f radians, from Radius %4.2f to %4.2f deg", i, phi_nbins, phi_minphi, phi_maxphi, phi_minradius, phi_maxradius ) ;
        tmphist1 = new TH1F( hname, htitle, phi_nbins, phi_minphi, phi_maxphi ) ;
        hXYAccImgSelPhiDependentSlice.push_back( tmphist1 ) ;
        
        sprintf( hname, "hAccImgSelRadiusDependentSlice%03d_%d", centerdeg, i ) ;
        sprintf( htitle, "1D histogram from getXoff_derot() and getYoff_derot(), with ImgSel=%d, RadiusDependentSlice%03d (Pie Slice)\n, from Phi %3.1f to %3.1f radians", i, centerdeg, centerdeg - ( rad_phiwidth * 180 / 3.1415 ), centerdeg + ( rad_phiwidth * 180 / 3.1415 ) ) ;
        tmphist1 = new TH1F( hname, htitle, rad_nbins, rad_minrad, rad_maxrad ) ;
        hXYAccImgSelRadiusDependentSlice000.push_back( tmphist1 ) ;
        
    }
    
    
    // 2D histograms, sorted by NImages
    for( int i = 0 ; i <= 4 ; i++ )
    {
        sprintf( hname, "hAccNImages_%d", i ) ;
        sprintf( htitle, "2D histogram of getXoff_derot() and getYoff_derot(), with NImages=%d", i ) ;
        tmphist2 = new TH2F( hname, htitle, hn, -1.0 * hr, hr, hn, -1.0 * hr, hr ) ;
        hXYAccNImages.push_back( tmphist2 ) ;
        
        sprintf( hname, "hAccNImagesPreDeRot_%d", i ) ;
        sprintf( htitle, "2D histogram of Xoff and Yoff, with NImages=%d", i ) ;
        tmphist2 = new TH2F( hname, htitle, hn, -1.0 * hr, hr, hn, -1.0 * hr, hr ) ;
        hXYAccNImagesPreDeRot.push_back( tmphist2 ) ;
        
        sprintf( hname, "hAccNImagesPhiDependentSlice_%d", i ) ;
        sprintf( htitle, "1D histogram from getXoff_derot() and getYoff_derot(), with NImages=%d, PhiDependentSlice\n, %d Bins from %3.1f to %3.1f radians, from Radius %4.2f to %4.2f deg", i, phi_nbins, phi_minphi, phi_maxphi, phi_minradius, phi_maxradius ) ;
        tmphist1 = new TH1F( hname, htitle, phi_nbins, phi_minphi, phi_maxphi ) ;
        hXYAccNImagesPhiDependentSlice.push_back( tmphist1 ) ;
        
        sprintf( hname, "hAccNImagesRadiusDependentSlice%03d_%d", centerdeg, i ) ;
        sprintf( htitle, "2D histogram of getXoff_derot() and getYoff_derot(), with NImages=%d, RadiusDependentSlice%03d (Pi Slice)\n from %3.1f to %3.1f radians", i, centerdeg, centerdeg - ( rad_phiwidth * 180 / 3.1415 ), centerdeg + ( rad_phiwidth * 180 / 3.1415 ) ) ;
        tmphist1 = new TH1F( hname, htitle, rad_nbins, rad_minrad, rad_maxrad ) ;
        hXYAccNImagesRadiusDependentSlice000.push_back( tmphist1 ) ;
        
    }
    
}


VRadialAcceptance::~VRadialAcceptance()
{
    if( fAccFile )
    {
        delete fAccFile;
    }
}


void VRadialAcceptance::reset()
{
    setProductionIO();
    fNumberOfRawFiles = 0.;
    
    fAcceptanceFunctionDefined = false;
    fRadialAcceptanceFit = 0;
    
    fRunPar = 0;
    fCuts = 0;
    hList = 0;
    hListNormalizeHistograms = 0;
    hListFitHistograms = 0;
    
    fSourcePosition_X = 0.;
    fSourcePosition_Y = 0.;
    fSourcePosition_Radius = 0.;
    fMaxDistanceAllowed = 5.;
    fCut_CameraFiducialSize_max = fMaxDistanceAllowed;
    
    hscale = 0;
    hPhiDist = 0;
    hPhiDistDeRot = 0;
    hXYAccTotDeRot = 0;
    fAccFile = 0;
    
    f2DAcceptanceMode = 0 ;
    f2DBinNormalizationConstant = 0 ;
    
    eventcount = 0 ;
    
    hXYAccImgSel.clear();
    hXYAccNImages.clear();
    hXYAccImgSelPreDeRot.clear();
    hXYAccNImagesPreDeRot.clear();
    
    setEnergyReconstructionMethod();
    setAzCut();
    
    fExtraHistogramMode = 0 ;
    fExtraHistogramDir = "" ;
}


/*!

    get radial acceptance

    (ignore here any zenith angle acceptance)

    note: x,y are in derotated coordinates
 */
double VRadialAcceptance::getAcceptance( double x, double y )
{
    // use 1D radial acceptances
    if( f2DAcceptanceMode == 0 )
    {
        double iacc = 1.;
        
        if( fAcceptanceFunctionDefined && fRadialAcceptanceFit )
        {
            double idist = sqrt( x * x + y * y );
            if( idist > fRadialAcceptanceFit->GetXmax() )
            {
                iacc = 0.;
            }
            else
            {
                iacc = fRadialAcceptanceFit->Eval( idist );
            }
            if( iacc > 1. )
            {
                iacc = 1.;
            }
            if( iacc < 0. )
            {
                iacc = 0.;
            }
        }
        
        return iacc;
        
    }
    // use 2D acceptances
    // (use getXoff_derot() and getYoff_derot())
    else if( f2DAcceptanceMode == 1 )
    {
        double iacc = 1.0 ;
        int xbin = hXYAccTotDeRot->GetXaxis()->FindBin( x ) ;
        int ybin = hXYAccTotDeRot->GetYaxis()->FindBin( y ) ;
        iacc = hXYAccTotDeRot->GetBinContent( xbin, ybin ) / f2DBinNormalizationConstant ;
        if( iacc > 1. )
        {
            iacc = 1.;
        }
        if( iacc < 0. )
        {
            iacc = 0.;
        }
        return iacc;
    }
    else
    {
        cout << "ERROR: getAcceptance() not defined for f2DAcceptanceMode > 0 " << endl;
        exit( EXIT_FAILURE ) ;
    }
    
    return -1000.0 ;
}


/*!
 *    define here region in the sky which are excluded in the analysis
 *
 *     x,y are camera coordinates relative to camera centre
 *     (not wobble shifted)
 */
bool VRadialAcceptance::isExcluded( double x, double y )
{
    if( ( x * x + y * y ) > fMaxDistanceAllowed * fMaxDistanceAllowed )
    {
        return true;
    }
    
    return false;
}


/*!
 *    define here regions which are to be excluded from background analysis
 *
 *    - source region
 *    - out of camera
 *    - other excluded regions
 *
 *     x,y are de-rotated camera coordinates
 *     (not wobble shifted: relative to the camera center)
 */
bool VRadialAcceptance::isExcludedfromBackground( double x, double y )
{
    // event outside fiducial area
    if( isExcluded( x, y ) )
    {
        return true;
    }
    
    
    // source region (exclusion radius around source position)
    if( ( x - fSourcePosition_X ) * ( x - fSourcePosition_X )
            + ( y - fSourcePosition_Y ) * ( y - fSourcePosition_Y )
            < ( fSourcePosition_Radius ) * ( fSourcePosition_Radius ) )
    {
        return true;
    }
    
    // loop over regions to exclude from background
    // (read from e.g. runparameter file)
    for( unsigned int i = 0; i < fListOfExclusionRegions.size(); i++ )
    {
        if( fListOfExclusionRegions[i] )
        {
            if( fListOfExclusionRegions[i]->isInsideExclusionRegion( x, y ) )
            {
                return true;
            }
        }
    }
    
    return false;
}


/*
 *     x,y are camera coordinates (not wobble shifted)
 */
bool VRadialAcceptance::isExcludedfromSource( double x, double y )
{
    return isExcluded( x, y );
}

/*

   set the position of the potential gamma-ray source in camera coordinates

*/
void VRadialAcceptance::setSource( double x, double y, double r, double imaxdist )
{
    // coordinates of exclusion region around potential gamma-ray
    // source
    fSourcePosition_X = x;
    fSourcePosition_Y = y;
    fSourcePosition_Radius = r;
    // maximum allowed distance of an event
    // from the camera centre
    fMaxDistanceAllowed = imaxdist;
}

/*

     set list of exclusion regions

*/
void VRadialAcceptance::setRegionToExcludeAcceptance( vector< VListOfExclusionRegions* > iF )
{
    fListOfExclusionRegions = iF;
}


/*

    apply gamma/hadron cuts and fill radial acceptance histograms

*/
int VRadialAcceptance::fillAcceptanceFromData( CData* iData, int entry, double x_rotJ2000, double y_rotJ2000 )
{
    if( !iData )
    {
        cout << "VRadialAcceptance::fillAcceptanceFromData: no data tree defined" << endl;
        return -1;
    }
    
    double idist = 0;
    double i_Phi = 0.;
    bool bPassed = false;
    
    // apply some basic quality cuts
    if( fCuts->applyInsideFiducialAreaCut() && fCuts->applyStereoQualityCuts( fEnergyReconstructionMethod, false, entry, true ) )
    {
        // gamma/hadron cuts
        if( !fCuts->isGamma( entry, false ) )
        {
            return 0;
        }
        
        // az cut
        bool bFill = false;
        if( fAzCut_min < fAzCut_max )
        {
            if( iData->getAz() > fAzCut_min && iData->getAz() <= fAzCut_max )
            {
                bFill = true;
            }
        }
        else
        {
            if( iData->getAz() < fAzCut_max || iData->getAz() > fAzCut_min )
            {
                bFill = true;
            }
        }
        if( !bFill )
        {
            return 0;
        }
        /////////////////////////////////////////////////
        // no more cuts should be applied after this statement
        bPassed = true;
        /////////////////////////////////////////////////
        
        idist = sqrt( x_rotJ2000 * x_rotJ2000 + y_rotJ2000 * y_rotJ2000 );
        i_Phi = atan2( y_rotJ2000, x_rotJ2000 ); // radians
        if( i_Phi < 0.0 )
        {
            i_Phi += 2 * TMath::Pi() ;    // atan2 is from -pi to pi, we want 0 to 2pi
        }
        
        //////////////////////////////////////////
        // fill run dependent radial acceptances
        for( unsigned int j = 0; j < fRunPar->fRunList.size(); j++ )
        {
            if( iData->runNumber == fRunPar->fRunList[j].fRunOff )
            {
                if( idist > 0. && j < fRadialAcceptance_perRun.size() )
                {
                    fRadialAcceptance_perRun[j]->Fill( idist );
                }
                if( j < hXYAccRun.size() )
                {
                    hXYAccRun[j]->Fill( x_rotJ2000, iData->getYoff() );
                }
                break;
            }
        }
        
        // fill 2D distribution of events
        hXYAccTotDeRot->Fill( x_rotJ2000, y_rotJ2000 );
        
        hXYAccImgSel[iData->getImgSel()]->Fill( x_rotJ2000, y_rotJ2000 ) ;
        hXYAccImgSelPreDeRot[iData->getImgSel()]->Fill( x_rotJ2000, iData->getYoff() ) ;
        hXYAccNImages[iData->getNImages()]->Fill( x_rotJ2000, y_rotJ2000 ) ;
        hXYAccNImagesPreDeRot[iData->getNImages()]->Fill( x_rotJ2000, iData->getYoff() ) ;
        
        // PhiDependentSlice Fill
        if( idist > phi_minradius && idist < phi_maxradius )
        {
            hXYAccTotDeRotPhiDependentSlice->Fill( i_Phi ) ;
            hXYAccImgSelPhiDependentSlice[iData->getImgSel()]->Fill( i_Phi ) ;
            hXYAccNImagesPhiDependentSlice[iData->getNImages()]->Fill( i_Phi ) ;
        }
        
        // RadiusDependentSlice000 Fill
        if( i_Phi > 0.0 - rad_phiwidth && i_Phi < 0.0 + rad_phiwidth )
        {
            hXYAccTotDeRotRadiusDependentSlice000->Fill( idist ) ;
            hXYAccImgSelRadiusDependentSlice000[iData->getImgSel()]->Fill( idist ) ;
            hXYAccNImagesRadiusDependentSlice000[iData->getNImages()]->Fill( idist ) ;
        }
        
        // fill azimuth angle dependend histograms (camera coordinates)
        hPhiDist->Fill( i_Phi * TMath::RadToDeg() );
        
        for( unsigned int j = 0; j < fPhiMin.size(); j++ )
        {
            bFill = false;
            if( i_Phi > fPhiMin[j] && i_Phi < fPhiMax[j] )
            {
                bFill = true;
            }
            else
            {
                if( fPhiMin[j] > fPhiMax[j] )
                {
                    if( i_Phi < fPhiMin[j] && i_Phi > fPhiMax[j] )
                    {
                        bFill = false;
                    }
                    else
                    {
                        bFill = true;
                    }
                }
            }
            if( bFill && idist > 0. )
            {
                hAccPhi[j]->Fill( idist );
            }
        }
        // fill azimuth angle dependend histograms (derotated camera coordinates)
        i_Phi = atan2( y_rotJ2000, x_rotJ2000 ) * TMath::RadToDeg();
        for( unsigned int j = 0; j < fPhiMin.size(); j++ )
        {
            bool bFill = false;
            if( i_Phi > fPhiMin[j] && i_Phi < fPhiMax[j] )
            {
                bFill = true;
            }
            else
            {
                if( fPhiMin[j] > fPhiMax[j] )
                {
                    if( i_Phi < fPhiMin[j] && i_Phi > fPhiMax[j] )
                    {
                        bFill = false;
                    }
                    else
                    {
                        bFill = true;
                    }
                }
            }
            if( bFill &&  idist > 0. )
            {
                hAccPhiDerot[j]->Fill( idist );
            }
        }
        
    }
    
    if( bPassed )
    {
        return 1;
    }
    
    return 0;
}

/*
    normalize 1D acceptances, fit them and write them to disk

    called for making radial acceptances

*/
bool VRadialAcceptance::terminate( TDirectory* iDirectory )
{
    if( !iDirectory->cd() )
    {
        cout << "VRadialAcceptance::terminate() error accessing directory  " << iDirectory->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    /////////////////////////////////////
    // normalize radial acceptance histograms
    // scale everything to mean value of first three bins
    if( fAccZeFitMinBin == fAccZeFitMaxBin )
    {
        cout << "Error: normalisation range for acceptance curves not well defined: ";
        cout << fAccZeFitMinBin << "\t" << fAccZeFitMaxBin << endl;
    }
    else
    {
        double isc = 0.;
        double i_normBin = ( double )( fAccZeFitMaxBin - fAccZeFitMinBin );
        // scale all histograms in hListNormalizeHistograms
        unsigned int i = 0;
        TIter next( hListNormalizeHistograms );
        while( TH1F* h = ( TH1F* )next() )
        {
            isc = 0.;
            i_normBin = 0.;
            if( i == 0 && h->GetEntries() > 0 )
            {
                cout << "VRadialAcceptance::terminate: scaling histograms to bins ";
                cout << fAccZeFitMinBin << " to " << fAccZeFitMaxBin << endl;
            }
            for( unsigned int j = fAccZeFitMinBin; j < fAccZeFitMaxBin; j++ )
            {
                if( h->GetBinError( j ) > 0. )
                {
                    isc +=  h->GetBinContent( j ) / ( h->GetBinError( j ) * h->GetBinError( j ) );
                    i_normBin += 1. / h->GetBinError( j ) / h->GetBinError( j );
                }
            }
            if( i_normBin > 0. )
            {
                isc /= i_normBin;
            }
            if( isc > 0 )
            {
                h->Scale( 1. / isc );
            }
            i++;
        }
    }
 
    // try to copy with printouts from Minuit
    gPrintViaErrorHandler = kTRUE;
    
    /////////////////////////////////////
    // analyze and fit histograms
    string i_hname;
    TIter next( hListFitHistograms );
    while( TH1F* h = ( TH1F* )next() )
    {
        // require a few entries in the histogram for the fit
        if( h->GetEntries() < 5 )
        {
            continue;
        }
        // fit function
        i_hname = h->GetName();
        i_hname.replace( 0, 1, "f" );
        TF1* ffit = new TF1( i_hname.c_str(), VRadialAcceptance_fit_acceptance_function, 0., fCut_CameraFiducialSize_max, 5 );
        ffit->SetNpx( 1000 );
        ffit->SetParameter( 0, -0.3 );
        ffit->SetParameter( 1, -0.6 );
        ffit->SetParameter( 2, +0.6 );
        ffit->SetParameter( 3, -0.2 );
        ffit->SetParameter( 4, 0.0 );
        ffit->SetParLimits( 4, 0., 1.0 );
        hList->Add( ffit );
        histListProductionIO->Add( ffit );
        // fit histogram
        i_hname = h->GetName();
        i_hname += "Fit";
        TH1F* hfit = new TH1F( i_hname.c_str(), h->GetTitle(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
        hfit->SetXTitle( h->GetXaxis()->GetTitle() );
        hfit->SetYTitle( h->GetYaxis()->GetTitle() );
        hList->Add( hfit );
        
        // fill the fitting histogram
        for( int j = 1; j < hfit->GetNbinsX(); j++ )
        {
            hfit->SetBinContent( j, h->GetBinContent( j ) ) ;
        }
        
        // fit the data and fill histograms
        cout << "fitting acceptance curves (" << h->GetName() << ") ..." << endl << endl;
        double i_eval = 0.;
        gErrorIgnoreLevel = 5000;
        hfit->Fit( ffit, "0REMQ" );
        gErrorIgnoreLevel = 0;
        // hfit->Fit( ffit, "0REM" );
        //ffit->GetParameter(4);
        hfit->SetBins( 1000, 0., 5. );
        // replace bin content by values from the fit function (set max to 1 and min to 0)
        for( int j = 1; j < hfit->GetNbinsX(); j++ )
        {
            i_eval = ffit->Eval( hfit->GetBinCenter( j ) );
            if( hfit->GetBinCenter( j ) > ffit->GetXmax() )
            {
                i_eval = 0.;
            }
            else if( hfit->GetBinCenter( j ) < ffit->GetMaximumX() )
            {
                i_eval = 1.;
            }
            else if( i_eval < 0. )
            {
                i_eval = 0.;
            }
            else if( i_eval >= 1. )
            {
                i_eval = 1.;
            }
            
            hfit->SetBinContent( j, i_eval );
        }
    }
    
    // write total number of entries in each histogram (ze) to screen
    if( fRadialAcceptance->GetEntries() > 0. )
    {
        cout << "total number of events in radial acceptance histogram: ";
        cout << fRadialAcceptance->GetEntries() << endl;
        cout << endl << "writing acceptance curves to following directory: " << iDirectory->GetName() << endl;
    }
    
    if( fproduction_shortIO )
    {
        histListProductionIO->Write();
    }
    else
    {
        hList->Write();
    }
    
    // write cuts to disk
    if( fCuts )
    {
        fCuts->SetName( "GammaHadronCuts" );
        fCuts->Write();
    }
    
    return true;
}



/*
   fit function for acceptance curves

   pol5 with x < x0 == 1 and f(x0)' == 0

   p0 = a2
*/
Double_t VRadialAcceptance_fit_acceptance_function( Double_t* x, Double_t* par )
{
    double f = 0.;
    
    double a2 = par[0];
    double a3 = par[1];
    double a4 = par[2];
    double a5 = par[3];
    double x0 = par[4];
    
    double a1 = -1.*( x0 * ( 2.*a2 + x0 * ( 3.*a3 + x0 * ( 4.*a4 + x0 * 5. * a5 ) ) ) );
    double a0 = 1. - ( x0 * ( a1 + x0 * ( a2 + x0 * ( a3 + x0 * ( a4 + x0 * a5 ) ) ) ) );
    
    if( x[0] < x0 )
    {
        f = 1.;
    }
    else
    {
        f = a0 + x[0] * ( a1 + x[0] * ( a2 + x[0] * ( a3 + x[0] * ( a4 + x[0] * a5 ) ) ) );
    }
    return f;
}

/////////////////////////////////////
// get normalization factor for 2d histogram bins
// 2dacceptance = bincontent(Xoff,Yoff)/normfactor
// normfactor = avg bin content from all bins < 0.3 deg from center
// using 2d histogram of derotated Xoff and Yoff CData->getXoff_derot(), CData->getYoff_derot()
double VRadialAcceptance::calculate2DBinNormalizationConstant( double radius ) // radius in degrees
{
    double normconst = 1.0 ;
    int nbinsx = hXYAccTotDeRot->GetNbinsX() ;
    int nbinsy = hXYAccTotDeRot->GetNbinsY() ;
    int nbins = nbinsx * nbinsy ;
    int i_binx, i_biny, i_binz ;
    double bincentx, bincenty, iradius, binCont;
    int avgbincount = 0 ;
    double avgbintotal = 0.0 ;
    
    // loop over all bins in hXYAccTot
    for( int i_bin = 1 ; i_bin <= nbins ; i_bin++ )
    {
        // find bin center and content
        hXYAccTotDeRot->GetBinXYZ( i_bin, i_binx, i_biny, i_binz ) ;
        bincentx = hXYAccTotDeRot->GetXaxis()->GetBinCenter( i_binx ) ;
        bincenty = hXYAccTotDeRot->GetYaxis()->GetBinCenter( i_biny ) ;
        iradius  = sqrt( bincentx * bincentx + bincenty * bincenty ) ;
        binCont  = hXYAccTotDeRot->GetBinContent( i_binx, i_biny ) ;
        if( iradius < radius )  // bin center < 0.3 deg
        {
            avgbintotal += binCont ;
            avgbincount++ ;
        }
    }
    
    // Write extra histograms
    if( fExtraHistogramMode > 0 )
    {
        char buff[250] ;
        cout << "VRadialAcceptance::calculate2DBinNormalizationConstant - writing hists" << endl;
        // write 2d hist to text file
        if( hXYAccTotDeRot != 0 )
        {
            cout << "hXYAccTotDeRot go" << endl;
            sprintf( buff, "%s/2dAcceptanceHist", fExtraHistogramDir.c_str() ) ;
            string outname = buff ;
            Write2DHistToTextFile( hXYAccTotDeRot, outname ) ;
        }
        
        if( hXYAccTotDeRotRadiusDependentSlice000 != 0 )
        {
            sprintf( buff, "%s/2dAcceptanceHist.RadiusDependentSlice000", fExtraHistogramDir.c_str() ) ;
            string buff2 = buff ;
            Write1DHistToTextFile( hXYAccTotDeRotRadiusDependentSlice000, buff2, 1 ) ;
        }
        
        if( hXYAccTotDeRotPhiDependentSlice != 0 )
        {
            sprintf( buff, "%s/2dAcceptanceHist.PhiDependentSlice", fExtraHistogramDir.c_str() ) ;
            string buff2 = buff ;
            Write1DHistToTextFile( hXYAccTotDeRotPhiDependentSlice, buff2, 1 ) ;
        }
        
        int i ;
        // write ImgSel 2d hists to text files
        if( hXYAccImgSel.empty() == 0 )
        {
            cout << "hXYAccImgSel go" << endl;
            for( i = 0 ; i <= 15 ; i++ )
            {
                sprintf( buff, "%s/ImgSel%d", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write2DHistToTextFile( hXYAccImgSel[i], buff2 );
            }
        }
        
        // ImgSel RadiusDependentSlice
        if( hXYAccImgSelRadiusDependentSlice000.empty() == 0 )
        {
            cout << "hXYAccImgSelRadiusDependentSlice000 go" << endl;
            for( i = 0 ; i <= 15 ; i++ )
            {
                sprintf( buff, "%s/ImgSel%d.RadiusDependentSlice000", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write1DHistToTextFile( hXYAccImgSelRadiusDependentSlice000[i], buff2, 1 ) ;
            }
        }
        
        // ImgSel RadiusDependentSlice
        if( hXYAccNImagesRadiusDependentSlice000.empty() == 0 )
        {
            cout << "hXYAccNImagesRadiusDependentSlice000 go" << endl;
            for( i = 0 ; i <= 4 ; i++ )
            {
                sprintf( buff, "%s/NImages%d.RadiusDependentSlice000", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write1DHistToTextFile( hXYAccNImagesRadiusDependentSlice000[i], buff2, 1 ) ;
            }
        }
        
        // PhiDependentSlice
        if( hXYAccImgSelPhiDependentSlice.empty() == 0 )
        {
            cout << "hXYAccImgSelPhiDependentSlice go" << endl;
            for( i = 0 ; i <= 15 ; i++ )
            {
                sprintf( buff, "%s/ImgSel%d.PhiDependentSlice", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write1DHistToTextFile( hXYAccImgSelPhiDependentSlice[i], buff2 , 2 ) ;
            }
        }
        if( hXYAccNImagesPhiDependentSlice.empty() == 0 )
        {
            cout << "hXYAccNImagesPhiDependentSlice go" << endl;
            for( i = 0 ; i <= 4 ; i++ )
            {
                sprintf( buff, "%s/NImages%d.PhiDependentSlice", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write1DHistToTextFile( hXYAccNImagesPhiDependentSlice[i], buff2 , 2 ) ;
            }
        }
        
        // write ImgSel 2d hists to text files
        cout << "hXYAccImgSelPreDeRot prep" << endl;
        if( hXYAccImgSelPreDeRot.empty() == 0 )
        {
            for( i = 0 ; i <= 15 ; i++ )
            {
                sprintf( buff, "%s/ImgSelPreDeRot%d", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write2DHistToTextFile( hXYAccImgSelPreDeRot[i], buff2 );
            }
        }
        
        // write NImages 2d hists to text files
        cout << "hXYAccNImages prep" << endl;
        if( hXYAccNImages.empty() == 0 )
        {
            cout << "hXYAccNImages go" << endl;
            for( i = 0 ; i <= 4 ; i++ )
            {
                sprintf( buff, "%s/NImages%d", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write2DHistToTextFile( hXYAccNImages[i], buff2 );
            }
        }
        
        // write NImages 2d hists to text files
        cout << "hXYAccNImagesPreDeRot prep" << endl;
        if( hXYAccNImagesPreDeRot.empty() == 0 )
        {
            cout << "hXYAccNImagesPreDeRot go" << endl;
            for( i = 0 ; i <= 4 ; i++ )
            {
                sprintf( buff, "%s/NImagesPreDeRot%d", fExtraHistogramDir.c_str(), i ) ;
                string buff2 = buff ;
                Write2DHistToTextFile( hXYAccNImagesPreDeRot[i], buff2 );
            }
        }
    }
    
    if( avgbincount <= 0 )
    {
        cout << "Error: calculate2DBinNormalizationConstant(" << radius;
        cout << ") : no bins of VRadialAcceptance->hXYAccTotDeRot within radius " << radius << " of center." << endl;
        return 1.0 ;
    }
    else
    {
        normconst = avgbintotal / avgbincount ;
        f2DBinNormalizationConstant = normconst ;
        return normconst ;
    }
    return normconst ;
    
}

// write out 1d hist to a text file
// mostly so other programs (i.e. mathematica) can work with the data
void VRadialAcceptance::Write1DHistToTextFile( TH1F* hist, string& basename, int histtype )
{
    cout << "Write1DHistToTextFile:" << basename.c_str() << endl ;
    // setup
    ofstream datafile ;
    char filename[200] ;
    sprintf( filename, "%s.dat", basename.c_str() ) ;
    datafile.open( filename ) ;
    
    // loop over bins
    int binx ;
    double bincenx, bincont ;
    char dataline[100] ;
    int nxbins = hist->GetNbinsX() ;
    int i ;
    for( binx = 1 ; binx <= nxbins ; binx++ )
    {
        i = ( int )hist->GetBin( binx ) ;
        if( hist->IsBinOverflow( i ) || hist->IsBinUnderflow( i ) )
        {
            continue ;
        }
        bincenx = hist->GetXaxis()->GetBinCenter( binx ) ;
        bincont = hist->GetBinContent( binx ) ;
        sprintf( dataline, "%d %d %f %f", i, binx, bincenx, bincont ) ;
        datafile << dataline << endl ;
    }
    datafile.close() ;
    
    // write metafile, a text file for storing extra info about the histogram
    // (anything that is not bin coordinates or content)
    cout << "fRunPar->fRunList.size(): prep..." << endl ;
    //cout << "fRunPar->fRunList.size():" << fRunPar->fRunList.size() << endl ;
    //cout << ", calculated from " << fNumberOfRawFiles << " files)" << endl;
    //int nruns = fRunPar->fRunList.size() ;
    int nruns = ( int )fNumberOfRawFiles ;
    cout << "fRunPar->fRunList.size(): " << fNumberOfRawFiles << endl ;
    //sprintf( hname, "fRadialAcceptance_perRun_%d", fRunPar->fRunList[i].fRunOff );
    double nentries = hist->GetEntries() ;
    ofstream metafile ;
    sprintf( filename, "%s.meta", basename.c_str() ) ;
    metafile.open( filename ) ;
    metafile << "datafileheaders i binx bincentx bincont" << endl ;
    metafile << "nxbins " << nxbins << endl ;
    metafile << "nentries " << nentries << endl ;
    metafile << "nruns " << nruns << endl ;
    if( histtype == 1 )   // 1 = RadiusDependentSlice
    {
        metafile << "rad_minrad "   << rad_minrad << endl;
        metafile << "rad_maxrad "   << rad_maxrad << endl;
        metafile << "rad_phiwidth " << rad_phiwidth << endl;
        metafile << "notes Slice of Pie! Bin position depends on Radius" << endl;
    }
    else if( histtype == 2 )  // 2 = PhiDependentSlice
    {
        metafile << "phi_minphi " << phi_minphi << endl;
        metafile << "phi_maxphi " << phi_maxphi << endl;
        metafile << "phi_minrad " << phi_minradius << endl;
        metafile << "phi_maxrad " << phi_maxradius << endl;
        metafile << "notes Doughnut! Bin position depends on Phi" << endl;
    }
    metafile << "ImgSelMeta ImgSelNumber T1234" << endl ;
    metafile << "ImgSel 0 0000 T" << endl ;
    metafile << "ImgSel 1 1000 T1" << endl ;
    metafile << "ImgSel 2 0100 T2" << endl ;
    metafile << "ImgSel 3 1100 T12" << endl ;
    metafile << "ImgSel 4 0010 T3" << endl ;
    metafile << "ImgSel 5 1010 T13" << endl ;
    metafile << "ImgSel 6 0110 T23" << endl ;
    metafile << "ImgSel 7 1110 T123" << endl ;
    metafile << "ImgSel 8 0001 T4" << endl ;
    metafile << "ImgSel 9 1001 T14" << endl ;
    metafile << "ImgSel 10 0101 T24" << endl ;
    metafile << "ImgSel 11 1101 T124" << endl ;
    metafile << "ImgSel 12 0011 T34" << endl ;
    metafile << "ImgSel 13 1011 T134" << endl ;
    metafile << "ImgSel 14 0111 T234" << endl ;
    metafile << "ImgSel 15 1111 T1234" << endl ;
    metafile.close() ;
}

// write out 2d hist to a text file
// mostly so other programs (i.e. mathematica) can work with the data
void VRadialAcceptance::Write2DHistToTextFile( TH2F* hist, string& basename )
{
    cout << "Write2DHistToTextFile:" << basename.c_str() << endl ;
    // setup
    ofstream datafile ;
    char filename[200] ;
    sprintf( filename, "%s.dat", basename.c_str() ) ;
    datafile.open( filename ) ;
    
    // loop over bins
    int binx, biny ;
    double bincenx, binceny, bincont ;
    char dataline[100] ;
    int nxbins = hist->GetNbinsX() ;
    int nybins = hist->GetNbinsY() ;
    int i ;
    for( binx = 1 ; binx <= nxbins ; binx++ )
    {
        for( biny = 1 ; biny <= nybins ; biny++ )
        {
            i = ( int )hist->GetBin( binx, biny ) ;
            if( hist->IsBinOverflow( i ) || hist->IsBinUnderflow( i ) )
            {
                continue ;
            }
            bincenx = hist->GetXaxis()->GetBinCenter( binx ) ;
            binceny = hist->GetYaxis()->GetBinCenter( biny ) ;
            bincont = hist->GetBinContent( binx, biny ) ;
            sprintf( dataline, "%d %d %d %f %f %f", i, binx, biny, bincenx, binceny, bincont ) ;
            datafile << dataline << endl ;
        }
    }
    datafile.close() ;
    
    // write metafile, a text file for storing extra info about the histogram
    // (anything that is not bin coordinates or content)
    cout << "fRunPar->fRunList.size(): prep..." << endl ;
    //cout << "fRunPar->fRunList.size():" << fRunPar->fRunList.size() << endl ;
    //cout << ", calculated from " << fNumberOfRawFiles << " files)" << endl;
    //int nruns = fRunPar->fRunList.size() ;
    int nruns = ( int )fNumberOfRawFiles ;
    cout << "fRunPar->fRunList.size(): " << fNumberOfRawFiles << endl ;
    //sprintf( hname, "fRadialAcceptance_perRun_%d", fRunPar->fRunList[i].fRunOff );
    double nentries = hist->GetEntries() ;
    ofstream metafile ;
    sprintf( filename, "%s.meta", basename.c_str() ) ;
    metafile.open( filename ) ;
    metafile << "datafileheaders i binx biny bincentx bincenty bincont" << endl ;
    metafile << "nxbins " << nxbins << endl ;
    metafile << "nybins " << nybins << endl ;
    metafile << "nentries " << nentries << endl ;
    metafile << "nruns " << nruns << endl ;
    metafile << "ImgSelMeta ImgSelNumber T1234" << endl ;
    metafile << "ImgSel 0 0000 T" << endl ;
    metafile << "ImgSel 1 1000 T1" << endl ;
    metafile << "ImgSel 2 0100 T2" << endl ;
    metafile << "ImgSel 3 1100 T12" << endl ;
    metafile << "ImgSel 4 0010 T3" << endl ;
    metafile << "ImgSel 5 1010 T13" << endl ;
    metafile << "ImgSel 6 0110 T23" << endl ;
    metafile << "ImgSel 7 1110 T123" << endl ;
    metafile << "ImgSel 8 0001 T4" << endl ;
    metafile << "ImgSel 9 1001 T14" << endl ;
    metafile << "ImgSel 10 0101 T24" << endl ;
    metafile << "ImgSel 11 1101 T124" << endl ;
    metafile << "ImgSel 12 0011 T34" << endl ;
    metafile << "ImgSel 13 1011 T134" << endl ;
    metafile << "ImgSel 14 0111 T234" << endl ;
    metafile << "ImgSel 15 1111 T1234" << endl ;
    metafile.close() ;
    
}

/////////////////////////////////
// if set to >= 1, will use alternate 2D acceptance radial acceptances
int VRadialAcceptance::Set2DAcceptanceMode( int mode )
{
    if( mode >= 0 )
    {
        f2DAcceptanceMode = mode ;
    }
    
    if( f2DAcceptanceMode > 0 )
    {
    
        // load 2d hist hXYAccTotDeRot from file
        hXYAccTotDeRot = ( TH2F* )gDirectory->Get( "hXYAccTotDeRot" );
        if( ! hXYAccTotDeRot )
        {
            cout << "Error, Radial Acceptance File " << fAccFile->GetName() << " does not contain the TH2F histogram 'hXYAccTotDeRot', required for calculating 2D acceptance.  Suggest using a newer Acceptance File, or use 1D Radial Acceptance mode instead of 2D Acceptance mode." << endl;
            exit( -1 ) ;
        }
        
        // calculate normalization constant
        calculate2DBinNormalizationConstant() ;
    }
    
    return f2DAcceptanceMode ;
}

void VRadialAcceptance::SetExtraHistogramDirectory( string histdir )
{
    fExtraHistogramMode = 1 ;
    fExtraHistogramDir = histdir ;
}

int VRadialAcceptance::SetExtraHistogramMode( int ehm )
{
    fExtraHistogramMode = ehm ;
    return fExtraHistogramMode ;
}


/*

     get scaling histogram for 1D radial acceptances taking
     exclusion regions into account

*/
bool VRadialAcceptance::correctRadialAcceptancesForExclusionRegions( TDirectory* iDirectory,
        unsigned int iRunNumber )
{
    if( !iDirectory->cd() )
    {
        cout << "VRadialAcceptance::scale_histograms_run() error accessing directory  " << iDirectory->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    Double_t xeventssim = 0.;
    Double_t yeventssim = 0.;
    int ifocus = -1;
    char hname[200];
    sprintf( hname, "hscaleRun_%d", iRunNumber );
    for( unsigned int i = 0; i < hscaleRun.size(); i++ )
    {
        if( string( hscaleRun[i]->GetName() ) == string( hname ) )
        {
            ifocus = i;
            break;
        }
    }
    /////////////////////////////////////////////
    // fill a 2D map and the 1D area map taking
    // exclusion regions into account
    cout << "2D filling of acceptance ratios" << endl;
    double neventssim = ( double )( hAreaExcluded2D[ifocus]->GetNbinsX() * hAreaExcluded2D[ifocus]->GetNbinsY() );
    cout << "\t neventssim = " << neventssim << endl;
    cout << "\t nxbins = " << hAreaExcluded2D[ifocus]->GetNbinsX() << endl ;
    cout << "\t nybins = " << hAreaExcluded2D[ifocus]->GetNbinsY() << endl ;
    int nfilled = 0 ;
    for( int kx = 1; kx <= hAreaExcluded2D[ifocus]->GetNbinsX(); kx++ )
    {
        xeventssim = hAreaExcluded2D[ifocus]->GetXaxis()->GetBinCenter( kx );
        for( int ky = 1; ky <= hAreaExcluded2D[ifocus]->GetNbinsY(); ky++ )
        {
            yeventssim = hAreaExcluded2D[ifocus]->GetYaxis()->GetBinCenter( ky );
            if( !isExcludedfromBackground( xeventssim, yeventssim ) )
            {
                hAreaExcluded2D[ifocus]->Fill( xeventssim, yeventssim );
                hscaleRun[ifocus]->Fill( TMath::Sqrt( xeventssim * xeventssim + yeventssim * yeventssim ) );
                nfilled += 1 ;
            }
        }
    }
    cout << "\t nfilled = " << nfilled << endl;
    
    // multiply with total area and scale by number of simulated events
    if( neventssim > 0. )
    {
        hscaleRun[ifocus]->Scale( ( ( 2.*hscaleRun[ifocus]->GetXaxis()->GetXmax() )
                                    * ( 2.*hscaleRun[ifocus]->GetXaxis()->GetXmax() ) ) / neventssim );
    }
    // ratio histogram is not used, just
    // used to illustrate differences due to
    // exclusion regions
    hscaleRunRatio[ifocus]->Divide( hscaleRun[ifocus] );
    
    ///////////////////////////////////////
    // correct 1D acceptances for exclusion
    // regions
    TIter next( hListNormalizeHistograms );
    char numberstring[200];
    while( TH1F* h = ( TH1F* )next() )
    {
        for( unsigned int j = 0; j < fRunPar->fRunList.size(); j++ )
        {
            sprintf( numberstring, "_%d", iRunNumber );
            if( ( int )iRunNumber == ( int )fRunPar->fRunList[j].fRunOff && strstr( h->GetName(), numberstring ) != NULL )
            {
                h->Divide( hscaleRun[ifocus] );
            }
        }
    }
    
    return true;
}


/*

    calculate average radial acceptance curve from single runs

*/
int VRadialAcceptance::calculateAverageRadialAcceptanceCurveFromRuns( TDirectory* iDirectory )
{
    if( !iDirectory->cd() )
    {
        cout << "VRadialAcceptance::calculateAcceptanceCurveFromRuns()";
        cout << " error accessing directory  " << iDirectory->GetName() << endl;
        exit( EXIT_FAILURE );
    }
    
    TIter next( hListNormalizeHistograms );
    char numberstring[200];
    
    // add up all runwise radial acceptance curves
    while( TH1F* h = ( TH1F* )next() )
    {
    
        sprintf( numberstring, "%s_", "hRadialAcceptance_perRun" );
        if( strstr( h->GetName(), numberstring ) != NULL )
        {
            fRadialAcceptance->Add( h );
        }
    }
    
    return 0;
}
