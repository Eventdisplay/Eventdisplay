/*! class VPlotAnasumHistograms

    plotting class for anasum analysis results


*/

# include "VPlotAnasumHistograms.h"
ClassImp( VPlotAnasumHistograms )

/*!


*/
VPlotAnasumHistograms::VPlotAnasumHistograms()
{
    fDebug = false;
    fAnasumDataFile = "";
    fRunNumber = -1;
    fPlotMode = "colz";
    
    setPlottingCorrelatedHistograms();
    setPlottingUseHours();
    setPlottingPSF_radius_deg();
    
    fSkyMapCentreDecJ2000 = -9999.;
    fSkyMapCentreRAJ2000  = -9999.;
    fTargetShiftWest   = -9999.;
    fTargetShiftNorth  = -9999.;
    
    default_settings();
}


VPlotAnasumHistograms::VPlotAnasumHistograms( string ifile, int ion )
{

    fDebug = false;
    fAnasumDataFile = ifile;
    fRunNumber = ion;
    fPlotMode = "colz";
    
    setPlottingCorrelatedHistograms();
    setPlottingUseHours();
    setPlottingPSF_radius_deg();
    
    fSkyMapCentreDecJ2000 = -9999.;
    fSkyMapCentreRAJ2000  = -9999.;
    fTargetShiftWest   = -9999.;
    fTargetShiftNorth  = -9999.;
    
    default_settings();
    
    if( !openDataFile( ifile ) )
    {
        return;
    }
}

bool VPlotAnasumHistograms::openDataFile( string ifile )
{
    if( !openFile( ifile, fRunNumber, 1 ) )
    {
        return false;
    }
    
    fSkyMapCentreDecJ2000 = getSkyMapCentreDecJ2000();
    fSkyMapCentreRAJ2000  = getSkyMapCentreRAJ2000();
    fTargetShiftWest   = getTargetShiftWest();
    fTargetShiftNorth  = getTargetShiftNorth();
    
    readRunList();
    
    return true;
}

///////////////////////////////////////////////////////////////////////////////
/*
 *   Help Function
 *
 *
 */
void VPlotAnasumHistograms::help()
{

    cout << "VPlotAnasumHistograms, plotting tool for anasum results" << endl;
    cout << "=======================================================" << endl;
    cout << endl;
    cout << "plot_mscPlots(char *mscwfile=0, int irebin=2, double xmin=-2, double xmax=4)" << endl;
    cout << "   plot mean scaled with and length histograms" << endl << endl;
    cout << "plot_qualityHistograms(double iSt = 1., bool bUpper = true, int iMethod=0)" << endl;
    cout << "   plot quality histograms for mscw, mscl, and theta2" << endl << endl;
    cout << "plot_theta2(double t2min=0, double t2max=0.3, irbin = 5)" << endl;
    cout << "   plot theta2 histograms" << endl << endl;
    cout << "plot_skyPlots()" << endl;
    cout << "   plot on, background, and on-background count maps, plot normalisation factors" << endl << endl;
    cout << "plot_skyPlots_significance(char* filename = \"output.root\", int ion = -1, bool iCorrelated = false,";
    cout << "double rmax = -1., double zmin = -100., double zmax = -100., bool bPoster = false)" << endl;
    cout << "plot_UncorrelatedSkyPlots()" << endl;
    cout << "   plot on, background, and on-background uncorrelated count maps, plot normalisation factors" << endl << endl;
    cout << "plot_CorrelatedSkyPlots()" << endl;
    cout << "   plot on, background, and on-background correlated count maps, plot normalisation factors" << endl << endl;
    cout << "   plot on-background and significance maps" << endl << endl;
    cout << "plot_significanceDistributions(double rmax = 1.2, double rSource = 0.4, double xmin=-6.5, double xmax=10)" << endl;
    cout << "   plot 1D distribution of signifcances in sky maps" << endl << endl;
    cout << "plot_radec( int sPlot = 0, double rmax = -3., double zmin = -5., double zmax = -1000., double xcenter = 0.,";
    cout << "double ycenter = 0., bool bDrawSource = true, bool bSlices = false, double fSliceXmin = -0.1, double fSliceXmax = 0.1, bool bProjX = true )" << endl;
    cout << "   plot skymap" << endl << endl;
    cout << "plot_catalogue( TCanvas *c, string iCatalogue = \"Hipparcos_MAG8_1997.dat\", double iMaxBrightness = 6.5, string iBand = B,";
    cout << "int iColor = 1, int iLineStyle = 1, string hSkyMapName = hmap_stereo_sig_REFLECTED, double iTextAngle = 45., int iMarkerStyle = 5 )" << endl;
    cout << "   plot the stars from a given magnitude and catalogue on top of the skymap" << endl << endl;
    cout << "plot_RBM_ring(double r, double iA, double t2, double iN)" << endl;
    cout << "   plot the ring background" << endl << endl;
    cout << "plot_reflectedRegions(TCanvas *iC, int i, int j, int iColor=5)" << endl;
    cout <<  "  plot regions for background determination into sky plots (reflected region model only)" << endl << endl;
    cout << "plot_excludedRegions(TCanvas *c)" << endl;
    cout << "   plot excluded regions for background determination into sky plots" << endl << endl ;
    
    cout << "plot_deadTimes()" << endl;
    cout << "   plot delta t and dead time histograms" << endl << endl;
    cout << "plot_cumulativeSignificance()" << endl;
    cout << "   plot cumulative significance vs time" << endl << endl;
    cout << "printRunList()" << endl;
    cout << "print list of runs" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/*
 *  calulate RA/Dec J2000 for a given x_derot, y_derot in camera coordinates
 *
 *  x_deg on RA axis, y_deg in dec axis
 *
 *  all coordinates in degrees
 *
 *  input and output in [deg]
 *
 */

void VPlotAnasumHistograms::convert_derotated_RADECJ2000( double x_deg, double y_deg, double xerr_deg, double yerr_deg )
{
    cout << "Sky map centre: " << fSkyMapCentreRAJ2000 << " " << fSkyMapCentreDecJ2000 << endl;
    
    double ra      = 0.;
    double dec     = 0.;
    double ra_err_deg  = 0.;
    double dec_err_deg = 0.;
    
    // incorrect: see http://veritash.sao.arizona.edu:8081/Eventdisplay-WG/1912
    // x_deg *= -1.;
    
    VAstronometry::vlaDtp2s( x_deg * TMath::DegToRad(), y_deg * TMath::DegToRad(), fSkyMapCentreRAJ2000 * TMath::DegToRad(), fSkyMapCentreDecJ2000 * TMath::DegToRad(), &ra, &dec );
    VAstronometry::vlaDtp2s( ( x_deg + xerr_deg )*TMath::DegToRad(), ( y_deg + yerr_deg )*TMath::DegToRad(), fSkyMapCentreRAJ2000 * TMath::DegToRad(), fSkyMapCentreDecJ2000 * TMath::DegToRad(), &ra_err_deg, &dec_err_deg );
    
    ra          *= TMath::RadToDeg();
    dec         *= TMath::RadToDeg();
    ra_err_deg  *= TMath::RadToDeg();
    dec_err_deg *= TMath::RadToDeg();
    
    ra_err_deg  = TMath::Abs( ra - ra_err_deg );
    dec_err_deg = TMath::Abs( dec - dec_err_deg );
    
    cout << "(RA,Dec) (J2000) for (x,y)=( " << x_deg << "+/-" << xerr_deg << " , " << y_deg << "+-" << yerr_deg << " ): ";
    cout << "( " << ra << "+/-" << ra_err_deg << " , " << dec << "+/-" << dec_err_deg <<  " )" << endl;
    
    double hours = ( double )( ( int )( ra * 24. / 360. ) );
    double min   = ( double )( ( int )( 60.*( ra * 24. / 360. - hours ) ) );
    double sec   = ( ra - hours * 360. / 24. - min * 360. / 24. / 60. ) * 24. / 360.*60.*60.;
    double dec_d = ( double )( int )( dec );
    double dec_m = ( double )( int )( ( dec - dec_d ) * 60. );
    double dec_s = ( dec - dec_d - dec_m / 60. ) * 3600.;
    
    double hours_err = ( double )( ( int )( ra_err_deg * 24. / 360. ) );
    double min_err   = ( double )( ( int )( 60.*( ra_err_deg * 24. / 360. - hours_err ) ) );
    double sec_err   = ( ra_err_deg - hours_err * 360. / 24. - min_err * 360. / 24. / 60. ) * 24. / 360.*60.*60.;
    double dec_d_err = ( double )( int )( dec_err_deg );
    double dec_m_err = ( double )( int )( ( dec_err_deg - dec_d_err ) * 60. );
    double dec_s_err = ( dec_err_deg - dec_d_err - dec_m_err / 60. ) * 3600.;
    
    
    cout << "(RA,Dec) (J2000) for (x,y)=( " << x_deg << "+/-" << xerr_deg << " , " << y_deg << "+-" << yerr_deg << " ): ";
    cout << "( " << hours << "h " << min << "' " << sec << "''" << " +/- " << hours_err << "h " << min_err << "' " << sec_err << "''" ;
    cout << ", ";
    if( dec_d > 0 )
    {
        cout << "+";
    }
    cout << dec_d << " " << dec_m << " " << dec_s << " +/- " << dec_d_err << " " << dec_m_err << " " << dec_s_err << " ";
    cout << " )" << endl;
    
    // calculating and printing the offset of the position wrt camera center
    double offset = sqrt( x_deg * x_deg + y_deg * y_deg );
    cout << "Offset from camera center = " << offset << " deg" << endl;
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *
 *
 *
 */
void VPlotAnasumHistograms::plot_mscPlots( int irebin, double xmin, double xmax, string mscwfile )
{

    char hname[200];
    char htitle[200];
    
    hmscw_diff = ( TH1D* )getHistogram( "hmscw_diff", fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( hmscw_diff, 1, 2, 1, 1, irebin, 1001 );
    hmscw_on = ( TH1D* )getHistogram( "hmscw_on", fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( hmscw_on, 1, 2, 1, 1, irebin, 3001 );
    hmscw_off = ( TH1D* )getHistogram( "hmscw_off", fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( hmscw_off, 8, 2, 1, 1, irebin, 3003 );
    hmscl_diff = ( TH1D* )getHistogram( "hmscl_diff", fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( hmscl_diff, 1, 2, 1, 1, irebin, 1001 );
    hmscl_on = ( TH1D* )getHistogram( "hmscl_on", fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( hmscl_on, 1, 2, 1, 1, irebin, 3001 );
    hmscl_off = ( TH1D* )getHistogram( "hmscl_off", fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( hmscl_off, 8, 2, 1, 1, irebin, 3003 );
    
    
    if( hmscw_diff && hmscw_on && hmscw_off &&  hmscl_diff && hmscl_on && hmscl_off )
    {
    
        sprintf( hname, "c_mscw_%d", fRunNumber );
        sprintf( htitle, "mean scaled width (%d)", fRunNumber );
        TCanvas* c_mscw = new TCanvas( hname, htitle, 10, 10, 900, 380 );
        c_mscw->Divide( 2, 1 );
        
        c_mscw->cd( 1 );
        hmscw_on->SetAxisRange( xmin, xmax );
        hmscw_on->SetTitle( "" );
        
        hmscw_on->Draw( "hist e" );
        hmscw_off->Draw( "same hist e" );
        
        c_mscw->cd( 2 );
        hmscw_diff->SetAxisRange( xmin, xmax );
        hmscw_diff->SetTitle( "" );
        hmscw_diff->SetFillColor( 9 );
        
        hmscw_diff->Draw( "hist e " );
        cout << "Mean scaled width: " << hmscw_diff->GetMean() << " +- " << hmscw_diff->GetRMS() << endl;
        sprintf( hname, "%s_fit", hmscw_diff->GetName() );
        
        VDouble_gauss* fdouble_gauss = new VDouble_gauss();
        
        TF1* hmscw_diff_fit = new TF1( hname, fdouble_gauss, -1.5, 1.5, 4, "VDouble_gauss" );
        //TF1 *hmscw_diff_fit = new TF1( hname,double_gauss, -1.5, 1.5, 4);
        hmscw_diff_fit->SetParameter( 0, hmscw_diff->GetMaximum() );
        hmscw_diff_fit->SetParameter( 1, 0. );
        hmscw_diff_fit->SetParameter( 2, 0.3 );
        hmscw_diff_fit->SetParameter( 3, 0.3 );
        hmscw_diff->Fit( hmscw_diff_fit->GetName(), "0RME" );
        if( hmscw_diff_fit )
        {
            hmscw_diff_fit->SetLineColor( 2 );
            hmscw_diff_fit->SetLineStyle( 2 );
            hmscw_diff_fit->Draw( "same" );
            cout << "Mean scaled width fit: " << hmscw_diff_fit->GetParameter( 1 ) << " +- " << hmscw_diff_fit->GetParError( 1 ) << endl;
        }
        
        TLine* lmscw_diff = new TLine( 0., hmscw_diff->GetMinimum(), 0., hmscw_diff->GetMaximum() );
        lmscw_diff->SetLineStyle( 2 );
        lmscw_diff->Draw();
        
        // end of mscw
        cout << endl << "================================================================" << endl << endl;
        
        sprintf( hname, "c_mscl_%d", fRunNumber );
        sprintf( htitle, "mean scaled length(%d)", fRunNumber );
        TCanvas* c_mscl = new TCanvas( hname, htitle, 610, 10, 900, 380 );
        c_mscl->Divide( 2, 1 );
        
        c_mscl->cd( 1 );
        hmscl_on->SetAxisRange( -2., 10. );
        hmscl_on->SetTitle( "" );
        
        hmscl_on->Draw( "hist e" );
        hmscl_off->Draw( "same hist e" );
        
        c_mscl->cd( 2 );
        hmscl_diff->SetAxisRange( xmin, xmax );
        hmscl_diff->SetTitle( "" );
        hmscl_diff->SetFillColor( 9 );
        
        hmscl_diff->Draw( "hist e " );
        sprintf( hname, "%s_fit", hmscl_diff->GetName() );
        TF1* hmscl_diff_fit = new TF1( hname, fdouble_gauss, -1.5, 2.5, 4, "VDouble_gauss" );
        hmscl_diff_fit->SetParameter( 0, hmscl_diff->GetMaximum() );
        hmscl_diff_fit->SetParameter( 1, 0. );
        hmscl_diff_fit->SetParameter( 2, 0.3 );
        hmscl_diff_fit->SetParameter( 3, 0.3 );
        hmscl_diff->Fit( hmscl_diff_fit->GetName(), "0RME" );
        if( hmscl_diff_fit )
        {
            hmscl_diff_fit->SetLineColor( 2 );
            hmscl_diff_fit->SetLineStyle( 2 );
            hmscl_diff_fit->Draw( "same" );
            cout << "Mean scaled length fit: " << hmscl_diff_fit->GetParameter( 1 ) << " +- " << hmscl_diff_fit->GetParError( 1 ) << endl;
        }
        cout << "Mean scaled length: " << hmscl_diff->GetMean() << " +- " << hmscl_diff->GetRMS() << endl;
        
        TLine* lmscl_diff = new TLine( 0., hmscl_diff->GetMinimum(), 0., hmscl_diff->GetMaximum() );
        lmscl_diff->SetLineStyle( 2 );
        lmscl_diff->Draw();
        
        c_mscw->Update();
        c_mscl->Update();
        
        if( mscwfile.size() > 0 )
        {
            TFile* fmscw = new TFile( mscwfile.c_str() );
            if( !fmscw->IsZombie() )
            {
                TH1D* mw = ( TH1D* )fmscw->Get( "hMSCWsim" );
                setHistogramPlottingStyle( mw, 2, 2, 1, 20, 10 );
                TH1D* ml = ( TH1D* )fmscw->Get( "hMSCLsim" );
                setHistogramPlottingStyle( ml, 2, 2, 1, 20, 4 );
                
                if( mw && ml )
                {
                    c_mscw->cd( 1 );
                    mw->Scale( hmscw_off->GetMaximum() / mw->GetMaximum() );
                    mw->Scale( 1.10 );
                    mw->Draw( "same e" );
                    
                    c_mscl->cd( 1 );
                    ml->Scale( hmscl_off->GetMaximum() / ml->GetMaximum() );
                    ml->Scale( 1.0 );
                    ml->Draw( "same" );
                }
            }
        }
    }
    
}

///////////////////////////////////////////////////////////////////////////////////
/*
 *
 *
 */
void VPlotAnasumHistograms::plot_qualityHistograms( double iSourceStrength, bool bUpper, int iMethod )
{
    char hname[200];
    char htitle[200];
    
    hmscw_diff = ( TH1D* )getHistogram( "hmscw_diff", fRunNumber, "stereoParameterHistograms" );
    hmscw_on = ( TH1D* )getHistogram( "hmscw_on", fRunNumber, "stereoParameterHistograms" );
    hmscw_off = ( TH1D* )getHistogram( "hmscw_off", fRunNumber, "stereoParameterHistograms" );
    hmscl_diff = ( TH1D* )getHistogram( "hmscl_diff", fRunNumber, "stereoParameterHistograms" );
    hmscl_on = ( TH1D* )getHistogram( "hmscl_on", fRunNumber, "stereoParameterHistograms" );
    hmscl_off = ( TH1D* )getHistogram( "hmscl_off", fRunNumber, "stereoParameterHistograms" );
    htheta2_diff = ( TH1D* )getHistogram( "htheta2_diff", fRunNumber, "stereoParameterHistograms" );
    htheta2_on = ( TH1D* )getHistogram( "htheta2_on", fRunNumber, "stereoParameterHistograms" );
    htheta2_off = ( TH1D* )getHistogram( "htheta2_off", fRunNumber, "stereoParameterHistograms" );
    
    if( hmscw_on && hmscw_off && hmscw_diff )
    {
        TH1D* hq = doQfactors( hmscw_on, hmscw_off, hmscw_diff, bUpper, iMethod, iSourceStrength );
        setHistogramPlottingStyle( hq, 8, 2, 1, 1, 1, 3001 );
        
        sprintf( hname, "c_qqmscw_%d", fRunNumber );
        sprintf( htitle, "mscw quality plot, run %d", fRunNumber );
        TCanvas* c_qqmscw = new TCanvas( hname, htitle, 410, 10, 400, 400 );
        c_qqmscw->Draw();
        
        hq->SetAxisRange( -2., 2. );
        
        hq->Draw();
        cout << hq->GetName() << "\t" << hq->GetBinCenter( hq->GetMaximumBin() ) << endl;
    }
    if( hmscl_on && hmscl_off && hmscl_diff )
    {
        TH1D* hq = doQfactors( hmscl_on, hmscl_off, hmscl_diff, bUpper, iMethod, iSourceStrength );
        setHistogramPlottingStyle( hq, 8, 2, 1, 1, 1, 3001 );
        
        sprintf( hname, "c_qqmscl_%d", fRunNumber );
        sprintf( htitle, "mscl quality plot, run %d", fRunNumber );
        TCanvas* c_qqmscl = new TCanvas( hname, htitle, 810, 10, 400, 400 );
        c_qqmscl->Draw();
        
        hq->SetAxisRange( -2., 2. );
        
        hq->Draw();
        cout << hq->GetName() << "\t" << hq->GetBinCenter( hq->GetMaximumBin() ) << endl;
    }
    if( htheta2_on && htheta2_off && htheta2_diff )
    {
        TH1D* hq = doQfactors( htheta2_on, htheta2_off, htheta2_diff, bUpper, iMethod, iSourceStrength );
        setHistogramPlottingStyle( hq, 8, 2, 1, 1, 1, 3001 );
        
        sprintf( hname, "c_qqtheta2_%d", fRunNumber );
        sprintf( htitle, "theta2 quality plot, run %d", fRunNumber );
        TCanvas* c_qqtheta2 = new TCanvas( hname, htitle, 10, 10, 400, 400 );
        c_qqtheta2->Draw();
        gPad->SetGridx( 1 );
        
        hq->SetAxisRange( 0., 0.1 );
        hq->SetTitle( "" );
        
        hq->Draw();
        cout << hq->GetName() << "\t" << hq->GetBinCenter( hq->GetMaximumBin() ) << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////
/*
 * plot an overview of the most important sky maps
 * (on, off, excess, alpha, significance maps, 1D significance distribution
 *
 */
TCanvas* VPlotAnasumHistograms::plot_skyPlots( string iPlotMode, bool iSingleCanvases,
        double excess_min, double excess_max,
        double sig_min, double sig_max )
{
    char hname[200];
    char htitle[200];
    
    ///////////////////////////////////////
    // get histograms from file
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_on" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_on" );
    }
    TH2D* hmap_stereo_on = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_stereo_on, -9999. );
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_off" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_off" );
    }
    TH2D* hmap_stereo_off = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_stereo_off, -9999. );
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_diff" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_diff" );
    }
    TH2D* hmap_stereo_diff = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_stereo_diff, 1.5 );
    if( excess_min > -9999. )
    {
        hmap_stereo_diff->SetMinimum( excess_min );
    }
    if( excess_max > -9999. )
    {
        hmap_stereo_diff->SetMaximum( excess_max );
    }
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_alphaUC_on" );
    }
    else
    {
        sprintf( hname, "hmap_alpha_on" );
    }
    TH2D* hmap_alpha_on = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_alpha_on, -9999. );
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_alphaUC_off" );
    }
    else
    {
        sprintf( hname, "hmap_alpha_off" );
    }
    TH2D* hmap_alpha_off = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_alpha_off, -9999. );
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_alphaNormUC_off" );
    }
    else
    {
        sprintf( hname, "hmap_alphaNorm_off" );
    }
    TH2D* hmap_alpha_diff = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_alpha_diff, 1.5 );
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_sig" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_sig" );
    }
    TH2D* hmap_stereo_sig = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    setHistogramPlottingStyle( hmap_stereo_sig, 1.5 );
    if( sig_min > -999. )
    {
        hmap_stereo_sig->SetMinimum( sig_min );
    }
    if( sig_max > -999. )
    {
        hmap_stereo_sig->SetMaximum( sig_max );
    }
    
    //////////////////////////////////////////////
    // plotting canvas
    TCanvas* c_skyAll = 0;
    vector< TCanvas* > cSky;
    if( iSingleCanvases )
    {
        for( unsigned int i = 0; i < 8; i++ )
        {
            if( fPlotCorrelated )
            {
                sprintf( hname, "c_skyUC_%d_%d", fRunNumber, i );
                sprintf( htitle, "sky/camera map  run %d (uncorrelated plot %d )", fRunNumber, i );
            }
            else
            {
                sprintf( hname, "c_sky_%d_%d", fRunNumber, i );
                sprintf( htitle, "sky/camera map  run %d (correlated plots %d)", fRunNumber, i );
            }
            cSky.push_back( new TCanvas( hname, htitle, 10 + i * 50, 10 + i * 50, 400, 400 ) );
            setPadMargins( cSky.back(), 1, -9999, 0.12 );
        }
    }
    else
    {
        if( fPlotCorrelated )
        {
            sprintf( hname, "c_skyUC_%d", fRunNumber );
            sprintf( htitle, "sky/camera map  run %d (uncorrelated plots)", fRunNumber );
        }
        else
        {
            sprintf( hname, "c_sky_%d", fRunNumber );
            sprintf( htitle, "sky/camera map  run %d (correlated plots)", fRunNumber );
        }
        c_skyAll = new TCanvas( hname, htitle, 10, 10, 1200, 600 );
        c_skyAll->Divide( 4, 2 );
        setPadMargins( c_skyAll, 6, -9999, 0.12 );
        for( unsigned int i = 0; i < 8; i++ )
        {
            c_skyAll->cd( i + 1 );
            cSky.push_back( ( TCanvas* )gPad );
        }
    }
    if( cSky.size() < 8 )
    {
        return 0;
    }
    
    cSky[0]->cd();
    if( hmap_stereo_on )
    {
        hmap_stereo_on->Draw( iPlotMode.data() );
    }
    
    cSky[1]->cd();
    if( hmap_stereo_off )
    {
        hmap_stereo_off->Draw( iPlotMode.data() );
    }
    
    cSky[2]->cd();
    if( hmap_stereo_diff )
    {
        hmap_stereo_diff->Draw( iPlotMode.data() );
    }
    
    cSky[3]->cd();
    if( hmap_stereo_sig )
    {
        hmap_stereo_sig->Draw( iPlotMode.data() );
    }
    
    cSky[4]->cd();
    if( hmap_alpha_on )
    {
        hmap_alpha_on->Draw( iPlotMode.data() );
    }
    
    cSky[5]->cd();
    if( hmap_alpha_off )
    {
        hmap_alpha_off->Draw( iPlotMode.data() );
        plot_excludedRegions( cSky[5] );
    }
    
    cSky[6]->cd();
    if( hmap_alpha_diff )
    {
        hmap_alpha_diff->Draw( iPlotMode.data() );
    }
    
    cSky[7]->cd();
    plot_significanceDistributions( 2.0, 0.4, -6.5, 10., ( TCanvas* )gPad );
    
    if( !c_skyAll )
    {
        return cSky[3];
    }
    return c_skyAll;
}

//////////////////////////////////////////////////////////////////////////////////
/*
 * plot uncorrelated sky plots
 *
 */
void VPlotAnasumHistograms::plot_UncorrelatedSkyPlots()
{
    setPlottingCorrelatedHistograms( true );
    plot_skyPlots( "colz", false );
}


//////////////////////////////////////////////////////////////////////////////////
/*
 * plot correlated sky plots
 *
 */
void VPlotAnasumHistograms::plot_CorrelatedSkyPlots()
{
    setPlottingCorrelatedHistograms( false );
    plot_skyPlots( "colz", false );
}


//////////////////////////////////////////////////////////////////////////////////
/*
 *
 * Plot cumulative significance vs time
 */
TCanvas* VPlotAnasumHistograms::plot_cumulativeSignificance( bool doSqrtFit )
{

    // canvas that will contain the cumulative significance plot
    TCanvas* c_t2 = new TCanvas( "cCumSig", "Cumulative Significance", 10, 10, 600, 400 );
    
    TGraph* g = calcCumulativeSig( );
    
    //fit with sqrt function. Fit parameter p0 is "significance per sqrt(h)".
    if( doSqrtFit && g->GetN() > 1 )
    {
        TF1* f = new TF1( "f", "TMath::Sqrt(x/60)*[0]", 0, g->GetX()[g->GetN() - 1] );
        f->SetParName( 0, "Significance per sqrt(h)" );
        f->SetLineWidth( 2 );
        f->SetLineColor( kRed );
        g->Fit( f );
    }
    g->Draw( "ALP" );
    return c_t2;
}


//////////////////////////////////////////////////////////////////////////////////
/*
 *
 *
 */
TCanvas* VPlotAnasumHistograms::plot_theta2( double t2min, double t2max, int irbin )
{
    // int iPlotPSF = 0;
    double setYMax = -1.;
    
    char hname[200];
    char htitle[200];
    
    sprintf( hname, "htheta2_diff" );
    htheta2_diff = ( TH1D* )getHistogram( hname, fRunNumber, "stereoParameterHistograms" );
    sprintf( hname, "htheta2_on" );
    htheta2_on = ( TH1D* )getHistogram( hname, fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( htheta2_on, 1, 2, 1, 1, irbin, 0 );
    sprintf( hname, "htheta2_off" );
    htheta2_off = ( TH1D* )getHistogram( hname, fRunNumber, "stereoParameterHistograms" );
    setHistogramPlottingStyle( htheta2_off, 14, 2, 1, 1, irbin, 3001 );
    
    // canvas that will contain the theta2 and the theta2_diff plots
    TCanvas* c_t2 = 0;
    
    if( htheta2_diff && htheta2_on && htheta2_off )
    {
        sprintf( hname, "c_t2_%d", fRunNumber );
        sprintf( htitle, "theta2 (run %d)", fRunNumber );
        c_t2 = new TCanvas( hname, htitle, 10, 10, 900, 400 );
        c_t2->Divide( 2, 1 );
        //c_t2->SetLeftMargin( 0.12 );
        
        c_t2->cd( 1 );
        htheta2_on->SetXTitle( "#Theta^{2} [deg^{2}]" );
        htheta2_on->SetTitle( "" );
        htheta2_off->SetTitle( "#Theta^{2} Histogram" );
        htheta2_off->SetTitle( "" );
        if( setYMax < 0. )
        {
            htheta2_off->SetMaximum( htheta2_on->GetMaximum() * 1.1 );
        }
        else
        {
            htheta2_off->SetMaximum( setYMax );
        }
        sprintf( hname, "Number of events / %.3f deg^{2}", htheta2_off->GetXaxis()->GetBinWidth( 2 ) );
        htheta2_on->SetYTitle( hname );
        htheta2_off->SetYTitle( hname );
        htheta2_off->GetYaxis()->SetTitleOffset( 1.4 );
        htheta2_off->SetMinimum( 0. );
        
        htheta2_off->SetAxisRange( t2min, t2max );
        htheta2_off->Draw( "hist e" );
        htheta2_on->Draw( "hist e same" );
        
        // get 68% containment radius (up to theta2 = 0.05deg2)
        double nt2 = 0.;
        for( int i = 1; i < htheta2_diff->GetXaxis()->FindBin( 0.05 ); i++ )
        {
            nt2 += htheta2_diff->GetBinContent( i );
        }
        double nt2_68 = 0.;
        if( nt2 > 0. )
        {
            for( int i = 1; i < htheta2_diff->GetXaxis()->FindBin( 0.05 ); i++ )
            {
                nt2_68 += htheta2_diff->GetBinContent( i );
                if( nt2_68 / nt2 > 0.68 )
                {
                    cout << "Theta2 containment radius (68%, binning dependend): " << htheta2_diff->GetXaxis()->GetBinLowEdge( i ) << " deg2" << endl;
                    break;
                }
            }
        }
        
        c_t2->cd( 2 );
        htheta2_diff->SetFillColor( 8 );
        htheta2_diff->SetXTitle( "#Theta^{2} [deg^{2}]" );
        htheta2_diff->SetTitle( "" );
        sprintf( hname, "Number of events / %.3f deg^{2}", htheta2_diff->GetXaxis()->GetBinWidth( 2 ) );
        htheta2_diff->SetYTitle( hname );
        htheta2_diff->GetYaxis()->SetTitleOffset( 1.4 );
        setHistogramPlottingStyle( htheta2_diff, 1, 2, 1, 1, irbin, 1001 );
        htheta2_diff->SetAxisRange( t2min, t2max );
        htheta2_diff->Draw( "hist e" );
        
        
        if( fDebug )
        {
            sprintf( hname, "htheta2_debug%d", fRunNumber );
            TH1D* htheta2_debug = new TH1D( hname, "", 50, -50., 50. );
            htheta2_debug->SetXTitle( "on - off" );
            htheta2_debug->SetYTitle( "number of events" );
            htheta2_debug->SetLineWidth( 3 );
            htheta2_debug->Sumw2();
            
            int istart = htheta2_diff->GetXaxis()->FindBin( 0.15 );
            int istopp = htheta2_diff->GetXaxis()->FindBin( 1. );
            for( int i = istart; i < istopp; i++ )
            {
                htheta2_debug->Fill( htheta2_diff->GetBinContent( i ) );
            }
            sprintf( hname, "c_t2_debug_%d", fRunNumber );
            sprintf( htitle, "theta2 (on-off) (run %d)", fRunNumber );
            TCanvas* c_t2_debug = new TCanvas( hname, htitle, 810, 10, 400, 400 );
            c_t2_debug->Draw();
            htheta2_debug->Draw( "ehist" );
            htheta2_debug->Fit( "gaus" );
            htheta2_debug->Draw( "ehist same" );
        }
    }
    else
    {
        cout << "histograms not found" << endl;
    }
    
    return c_t2;
}

///////////////////////////////////////////////////////////////////////////////
/*
    plot 1D significance distribution for all bins, all but an user-given exclusion region and
    all but the exclusion regions used in the analysis

    rmax:        maximum distance to sky plot centre of bins to be taken into account
    rSource:     minimum distance to sky plot centre of bins to be taken into account (exclude the source region)
    xmin, xmax:  plotting range (significances)
    regioncode:  extra specifier for additional region cuts
		  - exclude x<0 (only use top    half of skymap)
		  - exclude x>0 (only use bottom half of skymap)
*/
TH1D* VPlotAnasumHistograms::plot_significanceDistributions( double rmax, double rSource,
        double xmin, double xmax,
        TCanvas* cCanvas, bool regioncodeflag )

{
    char hname[200];
    char htitle[200];
    
    /////////////////////////
    // get histograms
    TH2D* hmap_stereo_sig = 0;
    TH2D* hmap_stereo_on = 0;
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_sig" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_sig" );
    }
    hmap_stereo_sig = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    
    if( fPlotCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_on" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_on" );
    }
    hmap_stereo_on = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    
    /////////////////////////
    // get exclusion regions
    TTree* t = ( TTree* )getHistogram( "tExcludedRegions", -1, "" );
    if( !t )
    {
        cout << "tree with excluded regions not found" << endl;
        return 0;
    }
    float x = 0.;
    float y = 0.;
    //float r = 0.;
    float r1 = 0.;
    float r2 = -99.;
    float theta = 0.;
    bool  bOldStyleExclusionRegions = false;
    t->SetBranchAddress( "x", &x );
    t->SetBranchAddress( "y", &y );
    // backward compatability with circular exclusion region
    if( t->GetBranchStatus( "r" ) )
    {
        t->SetBranchAddress( "r", &r1 );
        bOldStyleExclusionRegions = true;
    }
    else
    {
        t->SetBranchAddress( "r1", &r1 );
        t->SetBranchAddress( "r2", &r2 );
        t->SetBranchAddress( "theta", &theta );
    }
    const int iN = t->GetEntries();
    float* v_x = new float[iN];
    float* v_y = new float[iN];
    float* v_r1 = new float[iN];
    float* v_r2 = new float[iN];
    float* v_theta = new float[iN];
    cout << "Found " << iN << " exclusion regions" << endl;
    
    for( int i = 0; i < iN; i++ )
    {
        t->GetEntry( i );
        v_x[i] = x;
        v_y[i] = y;
        v_r1[i] = r1;
        if( bOldStyleExclusionRegions )
        {
            v_r2[i] = r1;
            v_theta[i] = 0.;
        }
        else
        {
            v_r2[i] = r2;
            v_theta[i] = theta;
        }
    }
    
    /////////////////////////
    // get 1D significance distributions
    
    // all entries in sky map
    TH1D* hsig_1DAll  = get_Bin_Distribution( hmap_stereo_sig, fRunNumber, rmax, 0., false, hmap_stereo_on );
    setHistogramPlottingStyle( hsig_1DAll, 2, 2, 2, 1, 1, 0 );
    if( hsig_1DAll )
    {
        hsig_1DAll->SetStats( 0 );
    }
    cout << "Plot Legend:" << endl;
    cout << "  red  :   with source region" << endl;
    cout << "           (use entire significance skymap)" << endl;
    
    // all entries in sky map excluding the source (ON region)
    TH1D* hsig_1D  = get_Bin_Distribution( hmap_stereo_sig, fRunNumber, rmax, rSource, false, hmap_stereo_on );
    setHistogramPlottingStyle( hsig_1D, 4, 2, 1, 1, 1, 0 );
    if( hsig_1D )
    {
        hsig_1D->SetStats( 0 );
    }
    cout << "  blue :   without source region" << endl;
    cout << "           (use entire skymap, except the ON region)" << endl;
    
    // all entries in sky map excluding the source and the exclusion regions
    // (this is the sky map on which all the Gaussian fits are applied)
    TH1D* hsig_1DExcluded = get_Bin_Distribution( hmap_stereo_sig, fRunNumber, rmax, rSource, false, hmap_stereo_on,
                            iN, v_x, v_y, v_r1, v_r2, v_theta );
    setHistogramPlottingStyle( hsig_1DExcluded, 1, 2, 2, 1, 1, 0 );
    if( hsig_1DExcluded )
    {
        hsig_1DExcluded->SetStats( 1 );
    }
    cout << "  black:   without source region and exclusion regions" << endl;
    cout << "           (use entire skymap, except the ON region and the excluded regions)" << endl;
    cout << "  green-dot-dashed: Gauss(0,1)" << endl;
    
    ///////////////////////////
    // Top half of skymap ONLY
    TH1D* hsig_1DTopOnly = 0;
    if( regioncodeflag )
    {
        hsig_1DTopOnly = get_Bin_Distribution( hmap_stereo_sig, fRunNumber, rmax, rSource, false, hmap_stereo_on,
                                               iN, v_x, v_y, v_r1, v_r2, v_theta, "a" );
        setHistogramPlottingStyle( hsig_1DTopOnly, kMagenta, 2, 2, 1, 1, 0 );
        if( hsig_1DTopOnly )
        {
            hsig_1DTopOnly->SetStats( 1 );
        }
        cout << "light purple:  without source region and exclusion regions, Top Half Only" << endl;
        cout << "               (use top half of skymap, except the ON region and the excluded regions)" << endl;
    }
    
    ///////////////////////////
    // Bottom half of skymap ONLY
    TH1D* hsig_1DBottomOnly = 0;
    if( regioncodeflag )
    {
        hsig_1DBottomOnly = get_Bin_Distribution( hmap_stereo_sig, fRunNumber, rmax, rSource, false, hmap_stereo_on,
                            iN, v_x, v_y, v_r1, v_r2, v_theta, "b" );
        setHistogramPlottingStyle( hsig_1DBottomOnly, kGreen + 3, 2, 2, 1, 1, 0 );
        if( hsig_1DBottomOnly )
        {
            hsig_1DBottomOnly->SetStats( 1 );
        }
        cout << "dark green:  without source region and exclusion regions, Bottom Half Only" << endl;
        cout << "             (use bottom half of skymap, except the ON region and the excluded regions)" << endl;
    }
    
    delete[] v_x;
    delete[] v_y;
    delete[] v_r1;
    delete[] v_r2;
    delete[] v_theta;
    
    gStyle->SetOptStat( "mr" );
    gStyle->SetOptFit( 1111 );
    
    /////////////////////////
    // fit function is a normal distribution with mean 0 and width 1
    TF1* fND = new TF1( "fND", "gaus(0)", -5., 5. );
    fND->FixParameter( 1, 0. );
    fND->FixParameter( 2, 1. );
    fND->SetLineColor( 8 );
    fND->SetLineWidth( 2 );
    fND->SetLineStyle( 6 );
    
    /////////////////////////
    // plotting of significance distributions
    TCanvas* c_sig1D = 0;
    if( hsig_1D && hsig_1DAll && hsig_1DExcluded )
    {
        hsig_1DExcluded->SetXTitle( "significance #sigma" );
        hsig_1DExcluded->SetYTitle( "# of entries" );
        hsig_1DExcluded->GetXaxis()->SetRangeUser( xmin, xmax );
        if( !cCanvas )
        {
            if( fPlotCorrelated )
            {
                sprintf( hname, "c_sig1DUC_%d", fRunNumber );
                sprintf( htitle, "significance distribution (1D)  run %d (uncorrelated plots)", fRunNumber );
            }
            else
            {
                sprintf( hname, "c_sig1D_%d", fRunNumber );
                sprintf( htitle, "significance distribution (1D) run %d (correlated plots)", fRunNumber );
            }
            c_sig1D = new TCanvas( hname, htitle, 450, 10, 400, 400 );
            c_sig1D->Draw();
            if( hsig_1DExcluded->GetEntries() > 0 )
            {
                c_sig1D->SetLogy( 1 );
            }
        }
        else
        {
            if( hsig_1DExcluded->GetEntries() > 0 )
            {
                cCanvas->SetLogy( 1 );
            }
            cCanvas->cd();
        }
        
        // fit a Gaussian distribution
        hsig_1DExcluded->Fit( fND );
        hsig_1DExcluded->Draw( "e hist" );
        fND->Draw( "same" );
        
        hsig_1DExcluded->Draw( "e hist same" );
        
        hsig_1DAll->Draw( "e hist same" );
        hsig_1D->Draw( "e hist same" );
        
        if( regioncodeflag && hsig_1DTopOnly && hsig_1DBottomOnly )
        {
            hsig_1DTopOnly->Draw( "e hist same" );
            hsig_1DBottomOnly->Draw( "e hist same" );
        }
        
        plotHistogramTitle( hsig_1DExcluded );
        
        /////////////////////////
        // ratio plot
        if( !cCanvas )
        {
            if( fPlotCorrelated )
            {
                sprintf( hname, "c_diff1DUC_%d", fRunNumber );
                sprintf( htitle, "ratio distribution (1D)  run %d (uncorrelated plots)", fRunNumber );
            }
            else
            {
                sprintf( hname, "c_diff1D_%d", fRunNumber );
                sprintf( htitle, "ratio distribution (1D) run %d (correlated plots)", fRunNumber );
            }
            
            TCanvas* c_diff1D = new TCanvas( hname, htitle, 900, 10, 400, 400 );
            c_diff1D->Draw();
            
            sprintf( hname, "%s_cl", hsig_1DExcluded->GetName() );
            TH1D* hDiff = ( TH1D* )hsig_1DExcluded->Clone( hname );
            hDiff->SetStats( 0 );
            for( int i = 1; i < hsig_1DExcluded->GetNbinsX(); i++ )
            {
                if( hsig_1DExcluded->GetBinContent( i ) > 0. && fND->Eval( hsig_1DExcluded->GetXaxis()->GetBinCenter( i ) ) > 0. )
                {
                    hDiff->SetBinContent( i, hsig_1DExcluded->GetBinContent( i )
                                          / fND->Eval( hsig_1DExcluded->GetXaxis()->GetBinCenter( i ) ) );
                    hDiff->SetBinError( i, 0. );
                }
            }
            if( hDiff->GetEntries() > 0 )
            {
                c_diff1D->SetLogy( 1 );
            }
            hDiff->SetXTitle( "significance #sigma" );
            hDiff->SetYTitle( "ratio (measured to background expected)" );
            hDiff->GetXaxis()->SetRangeUser( xmin, xmax );
            hDiff->Draw( "e hist" );
            TLine* iLL = new TLine( xmin, 1., xmax, 1. );
            iLL->SetLineStyle( 2 );
            iLL->Draw();
        }
    }
    return hsig_1DExcluded;
}


///////////////////////////////////////////////////////////////////////////////////
/*
 *   sPlot = 0: plot significance
 *         = 1: plot excess
 *         = 2: on events
 *         = 3: off events
 *
 *   rmax > 0:  circulate cut for 2D histogram
 *   rmax < 0:  square cut for 2D histograms
*/
TCanvas* VPlotAnasumHistograms::plot_radec( int sPlot, double rmax, double zmin, double zmax, double xcenter, double ycenter,
        bool bSlices, double fSliceXmin, double fSliceXmax, bool bProjX )
{
    cout << "OBSERVE: right ascension axis might on local time settings; adjust offset hours with setPlottingUseHours( bool iB, int iZeroHours )" << endl;
    
    // different presentations
    fPlotMode = "A colz";
    
    //
    if( fabs( xcenter ) > 1.e-5 || fabs( ycenter ) > 1.e-5 )
    {
        cout << "WARNING: sky map shifting preliminary. Attention on the edge of the sky map" << endl;
    }
    
    char hname[200];
    char htitle[200];
    
    ycenter *= -1.;
    
    // get all histograms
    TH2D* hmap = 0;
    
    if( fPlotCorrelated )
    {
        if( sPlot == 0 )
        {
            sprintf( hname, "hmap_stereoUC_sig" );
        }
        else if( sPlot == 1 )
        {
            sprintf( hname, "hmap_stereoUC_diff" );
        }
        else if( sPlot == 2 )
        {
            sprintf( hname, "hmap_stereoUC_on" );
        }
        else if( sPlot == 3 )
        {
            sprintf( hname, "hmap_stereoUC_off" );
        }
    }
    else
    {
        if( sPlot == 0 )
        {
            sprintf( hname, "hmap_stereo_sig" );
        }
        else if( sPlot == 1 )
        {
            sprintf( hname, "hmap_stereo_diff" );
        }
        else if( sPlot == 2 )
        {
            sprintf( hname, "hmap_stereo_on" );
        }
        else if( sPlot == 3 )
        {
            sprintf( hname, "hmap_stereo_off" );
        }
    }
    
    hmap = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    if( rmax > 0. )
    {
        hmap = removeOuterRing( hmap, rmax, -9999. );
    }
    if( rmax < 0. )
    {
        rmax *= -1.;
    }
    setHistogramPlottingStyle( hmap, 1.5 );
    hmap = reflectXaxis( hmap );
    
    // significance
    TCanvas* c_skysig  = 0;
    if( hmap )
    {
        sprintf( hname, "c_skysig_%d_%d", fRunNumber, sPlot );
        sprintf( htitle, "sky map, run %d", fRunNumber );
        c_skysig  = new TCanvas( hname, htitle, 10, 10, 400, 400 );
        c_skysig->Draw();
        c_skysig->SetRightMargin( 0.14 );
        c_skysig->SetLeftMargin( 0.11 );
        
        hmap->SetTitle( "" );
        
        double iYRange = 0.;
        double x1 = -1.*rmax - xcenter;
        double x2 = rmax - xcenter;
        double y1 = -1.*rmax - ycenter;
        double y2 = rmax - ycenter;
        if( x1 < hmap->GetXaxis()->GetXmin() )
        {
            x1 = hmap->GetXaxis()->GetXmin();
        }
        if( x2 > hmap->GetXaxis()->GetXmax() )
        {
            x2 = hmap->GetXaxis()->GetXmax();
        }
        if( y1 < hmap->GetYaxis()->GetXmin() )
        {
            y1 = hmap->GetYaxis()->GetXmin();
        }
        if( y2 > hmap->GetYaxis()->GetXmax() )
        {
            y2 = hmap->GetYaxis()->GetXmax();
        }
        hmap->GetXaxis()->SetRangeUser( x1, x2 );
        hmap->GetYaxis()->SetRangeUser( y1, y2 );
        
        if( zmin > -1000. )
        {
            hmap->SetMinimum( zmin );
        }
        if( zmax > -1000. )
        {
            hmap->SetMaximum( zmax );
        }
        
        if( sPlot == 0 )
        {
            hmap->SetZTitle( "significance [#sigma]" );
        }
        else if( sPlot == 1 )
        {
            hmap->SetZTitle( "excess events" );
        }
        else if( sPlot == 2 )
        {
            hmap->SetZTitle( "events" );
        }
        else if( sPlot == 3 )
        {
            hmap->SetZTitle( "background events" );
        }
        
        hmap->Draw( fPlotMode.data() );
        
        if( fPlotDrawPSF > 0. )
        {
            drawPSF( c_skysig, "", hmap, fPlotDrawPSF );
        }
        
        // now do the axis
        double xmin, ymin, xmax, ymax, wmin, wmax;
        
        double dec = fSkyMapCentreDecJ2000;
        double ra  = fSkyMapCentreRAJ2000;
        cout << "(ra,dec)_J2000 = (" << ra << ", " << dec << ")" << endl;
        
        // dec axis
        xmin = hmap->GetXaxis()->GetBinLowEdge( hmap->GetXaxis()->FindBin( x1 ) );
        ymin = hmap->GetYaxis()->GetBinLowEdge( hmap->GetYaxis()->FindBin( y1 ) );
        xmax = hmap->GetXaxis()->GetBinLowEdge( hmap->GetXaxis()->FindBin( x1 ) );
        ymax = hmap->GetYaxis()->GetBinLowEdge( hmap->GetYaxis()->FindBin( y2 ) );
        iYRange = ymax - ymin;
        wmin = dec + ymin;
        wmax = dec + ymax;
        
        cout << "DECAXIS " << xmin << " " << xmax << " " << ymin << " " << ymax << " " << wmin << " " << wmax << " " << iYRange << endl;
        
        TGaxis* decAxis = new TGaxis( xmin, ymin, xmax, ymax, wmin, wmax, 505 );
        decAxis->SetTitleFont( hmap->GetYaxis()->GetTitleFont() );
        decAxis->SetTitleSize( hmap->GetYaxis()->GetTitleSize() );
        decAxis->SetTitleOffset( hmap->GetYaxis()->GetTitleOffset() );
        decAxis->SetLabelFont( hmap->GetYaxis()->GetLabelFont() );
        decAxis->SetLabelSize( hmap->GetYaxis()->GetLabelSize() );
        decAxis->SetLabelColor( hmap->GetYaxis()->GetLabelColor() );
        decAxis->SetTextColor( hmap->GetYaxis()->GetTitleColor() );
        decAxis->SetLineColor( hmap->GetYaxis()->GetAxisColor() );
        decAxis->Draw();
        decAxis->SetTitle( "declination_{J2000} [deg]" );
        decAxis->SetTitleOffset( 1.25 );
        
        //////////////////////////////////////////////////////////////////
        // ra axis (code copied from M.Beilicke)
        
        xmin = hmap->GetXaxis()->GetBinLowEdge( hmap->GetXaxis()->FindBin( x1 ) );
        ymin = hmap->GetYaxis()->GetBinLowEdge( hmap->GetYaxis()->FindBin( y1 ) );
        xmax = hmap->GetXaxis()->GetBinLowEdge( hmap->GetXaxis()->FindBin( x2 ) );
        ymax = hmap->GetYaxis()->GetBinLowEdge( hmap->GetYaxis()->FindBin( y1 ) );
        
        double Xmin = -1.*ra;
        double Xmax = -1.*ra;
        
        if( cos( ( dec + ymin ) * TMath::Pi() / 180. ) )
        {
            Xmin += xmin / cos( ( dec + ymin ) * TMath::Pi() / 180. );
            Xmax += xmax / cos( ( dec + ymin ) * TMath::Pi() / 180. );
        }
        cout << "RAXIS " << xmin << " " << xmax << " " << ymin << " " << ymax << " " << Xmin << " " << Xmax << endl;
        
        char TmpTimeFormat[200];
        if( fabs( rmax ) < 0.3 )
        {
            sprintf( TmpTimeFormat, "%%H^{h}%%M^{m}%%S^{s}" );
        }
        else
        {
            sprintf( TmpTimeFormat, "%%H^{h}%%M^{m}" );
        }
        
        TF1* IncValues = 0;
        int iRA_hrs = 0;
        int iRA_min = 0;
        int iRA_sec = 0;
        if( fPlotUseHours == true )
        {
            // convert angle to hours/min/seconds
            char iSign[10];
            int ihmsf[4];
            VAstronometry::vlaDr2tf( 4, -1.*Xmax * TMath::DegToRad(), iSign, ihmsf );
            iRA_hrs = ihmsf[0];
            iRA_min = ihmsf[1];
            iRA_sec = ihmsf[2];
            double iRA_frac = 0.;
            if( ihmsf[3] > 0 )
            {
                iRA_frac = 1. / ( double )ihmsf[3];
            }
            
            // convert to seconds
            Xmin *= 86400.0 / 360.0;
            Xmax *= 86400.0 / 360.0;
            IncValues = new TF1( "IncValues", "-x", iRA_frac, Xmax - Xmin + iRA_frac );
        }
        else
        {
            IncValues = new TF1( "IncValues", "-x", -Xmin, -Xmax );
        }
	IncValues->Print();
        
        TGaxis* raLowerAxis = new TGaxis( xmin, ymin, xmax, ymax, "IncValues", 4 );
        raLowerAxis->SetTitleFont( hmap->GetXaxis()->GetTitleFont() );
        raLowerAxis->SetTitleSize( hmap->GetXaxis()->GetTitleSize() );
        raLowerAxis->SetTitleOffset( hmap->GetXaxis()->GetTitleOffset() * 1.1 );
        raLowerAxis->SetLabelFont( hmap->GetXaxis()->GetLabelFont() );
        raLowerAxis->SetLabelSize( hmap->GetXaxis()->GetLabelSize() );
        raLowerAxis->SetLabelColor( hmap->GetXaxis()->GetLabelColor() );
        raLowerAxis->SetLabelOffset( 0.011 );
        if( fPlotCorrelated )
        {
            raLowerAxis->SetLabelOffset( 0.015 );
        }
        raLowerAxis->SetTextColor( hmap->GetXaxis()->GetTitleColor() );
        raLowerAxis->SetLineColor( hmap->GetXaxis()->GetAxisColor() );
        
        if( fPlotUseHours == true )
        {
            // set time offsets to midnight
            iRA_hrs += fPlotZeroHours;
            TTimeStamp dt( 2000, 01, 01, iRA_hrs, iRA_min, iRA_sec );
            raLowerAxis->SetTimeOffset( dt.GetSec(), "gmt" );
            raLowerAxis->SetTimeFormat( TmpTimeFormat );
            raLowerAxis->SetOption( "t" );
            raLowerAxis->SetNdivisions( 5 );
            raLowerAxis->SetTitle( "right ascension_{J2000} [hours]" );
        }
        else
        {
            raLowerAxis->SetTitle( "right ascension_{J2000} [deg]" );
            raLowerAxis->SetNdivisions( 510 );
        }
        raLowerAxis->Draw();
        
        ////////////////////////////////
        // plot slices
        ////////////////////////////////
        
        if( bSlices )
        {
        
            // get on or off histogram (in case of on or off plot)
            TH2D* hmap_opp = 0;
            TH2D* hmap_alpha = 0;
            
            if( fPlotCorrelated )
            {
                if( sPlot == 2 )
                {
                    sprintf( hname, "hmap_stereoUC_off" );
                }
                else if( sPlot == 3 )
                {
                    sprintf( hname, "hmap_stereoUC_on" );
                }
                sprintf( htitle, "hmap_alphaNormUC_off" );
            }
            else
            {
                if( sPlot == 2 )
                {
                    sprintf( hname, "hmap_stereo_off" );
                }
                else if( sPlot == 3 )
                {
                    sprintf( hname, "hmap_stereo_on" );
                }
                sprintf( htitle, "hmap_alphaNorm_off" );
            }
            hmap_opp = ( TH2D* )gDirectory->Get( hname );
            hmap_alpha = ( TH2D* )gDirectory->Get( htitle );
            hmap_alpha = reflectXaxis( hmap_alpha );
            
            sprintf( hname, "c_skySlice_%d_%d", fRunNumber, sPlot );
            sprintf( htitle, "sky slice, run %d", fRunNumber );
            TCanvas* c_skysigX  = new TCanvas( hname, htitle, 410, 10, 400, 400 );
            c_skysigX->Draw();
            c_skysigX->SetRightMargin( 0.14 );
            c_skysigX->SetLeftMargin( 0.11 );
            
            TH1D* h1 = 0;
            TH1D* h2 = 0;
            TH1D* h3 = 0;
            sprintf( hname, "%s_Slice_%d", hmap->GetName(), gRandom->Integer( 100 ) );
            sprintf( htitle, "%s_SliceO", hmap->GetName() );
            if( bProjX )
            {
                if( sPlot == 0 )
                {
                    h1 = hmap->ProfileY( hname, hmap->GetXaxis()->FindBin( fSliceXmin ), hmap->GetXaxis()->FindBin( fSliceXmax ) );
                    if( hmap_opp )
                    {
                        h2 = hmap_opp->ProfileY( htitle, hmap_opp->GetXaxis()->FindBin( fSliceXmin ), hmap_opp->GetXaxis()->FindBin( fSliceXmax ) );
                    }
                }
                else
                {
                    h1 = hmap->ProjectionY( hname, hmap->GetXaxis()->FindBin( fSliceXmin ), hmap->GetXaxis()->FindBin( fSliceXmax ) );
                    if( hmap_opp )
                    {
                        h2 = hmap_opp->ProjectionY( htitle, hmap->GetXaxis()->FindBin( fSliceXmin ), hmap->GetXaxis()->FindBin( fSliceXmax ) );
                    }
                }
                if( hmap_alpha )
                {
                    sprintf( hname, "%s_SliceAlpha", hmap_alpha->GetName() );
                    h3 = hmap_alpha->ProjectionY( hname, hmap_alpha->GetXaxis()->FindBin( fSliceXmin ), hmap_alpha->GetXaxis()->FindBin( fSliceXmax ) );
                    if( ( hmap_alpha->GetXaxis()->FindBin( fSliceXmax ) - hmap_alpha->GetXaxis()->FindBin( fSliceXmin ) ) > 0 )
                    {
                        h3->Scale( 1. / ( double )( 1. + hmap_alpha->GetXaxis()->FindBin( fSliceXmax ) - hmap_alpha->GetXaxis()->FindBin( fSliceXmin ) ) );
                    }
                }
            }
            else
            {
                if( sPlot == 0 )
                {
                    h1 = hmap->ProfileX( hname, hmap->GetYaxis()->FindBin( fSliceXmin ), hmap->GetYaxis()->FindBin( fSliceXmax ) );
                    if( hmap_opp )
                    {
                        h2 = hmap_opp->ProfileX( htitle, hmap_opp->GetYaxis()->FindBin( fSliceXmin ), hmap_opp->GetYaxis()->FindBin( fSliceXmax ) );
                    }
                }
                else
                {
                    h1 = hmap->ProjectionX( hname, hmap->GetYaxis()->FindBin( fSliceXmin ), hmap->GetYaxis()->FindBin( fSliceXmax ) );
                    if( hmap_opp )
                    {
                        h2 = hmap_opp->ProjectionX( htitle, hmap_opp->GetYaxis()->FindBin( fSliceXmin ), hmap_opp->GetYaxis()->FindBin( fSliceXmax ) );
                    }
                }
                if( hmap_alpha )
                {
                    sprintf( hname, "%s_SliceAlpha", hmap_alpha->GetName() );
                    h3 = hmap_alpha->ProjectionX( hname, hmap_alpha->GetYaxis()->FindBin( fSliceXmin ), hmap_alpha->GetYaxis()->FindBin( fSliceXmax ) );
                    if( ( hmap_alpha->GetYaxis()->FindBin( fSliceXmax ) - hmap_alpha->GetYaxis()->FindBin( fSliceXmin ) ) > 0 )
                    {
                        h3->Scale( 1. / ( double )( 1. + hmap_alpha->GetYaxis()->FindBin( fSliceXmax ) - hmap_alpha->GetYaxis()->FindBin( fSliceXmin ) ) );
                    }
                }
            }
            
            if( zmin > -1000. )
            {
                h1->SetMinimum( zmin );
            }
            if( zmax > -1000. )
            {
                h1->SetMaximum( zmax );
            }
            h1->SetStats( 0 );
            h1->SetLineWidth( 2 );
            
            if( h2 )
            {
                h2->SetStats( 0 );
                h2->SetLineWidth( 2 );
                h2->SetLineColor( 2 );
                if( sPlot == 2 || sPlot == 3 )
                {
                    if( h3 )
                    {
                        for( int i = 1; i <= h2->GetNbinsX(); i++ )
                        {
                            h2->SetBinContent( i, h2->GetBinContent( i ) * h3->GetBinContent( i ) );
                        }
                    }
                }
            }
            if( rmax > 0. )
            {
                h1->SetAxisRange( -1.*rmax, rmax );
            }
            else
            {
                h1->SetAxisRange( rmax, -1.*rmax );
            }
            
            h1->Draw( "e" );
            if( h2 )
            {
                h2->Draw( "same" );
            }
            c_skysigX->Update();
            
            // draw with right axis
            
            sprintf( hname, "c_skySliceRA_%d_%d", fRunNumber, sPlot );
            sprintf( htitle, "sky slice, run %d", fRunNumber );
            TCanvas* c_skysigY  = new TCanvas( hname, htitle, 810, 10, 400, 400 );
            c_skysigY->Draw();
            c_skysigY->SetRightMargin( 0.14 );
            c_skysigY->SetLeftMargin( 0.11 );
            
            h1->Draw( "AH e" );
            if( h2 )
            {
                h2->Draw( "same" );
            }
            c_skysigY->Update();
            
            iYRange  = h1->GetYaxis()->GetXmax() - h1->GetYaxis()->GetXmin();
            if( bProjX )
            {
                TGaxis* decAxis1D = ( TGaxis* )decAxis->Clone();
                decAxis1D->SetTitleOffset( 1.0 );
                
                decAxis1D->DrawAxis( c_skysigY->GetUxmin(), c_skysigY->GetUymin(), c_skysigY->GetUxmax(), c_skysigY->GetUymin(), dec - iYRange / 2., dec + iYRange / 2., 505 );
            }
            else
            {
                TGaxis* raLowerAxis1D = ( TGaxis* )raLowerAxis->Clone();
                
                raLowerAxis1D->Draw();
                
            }
            TGaxis* h1YAxis = new TGaxis( c_skysigY->GetUxmin(), c_skysigY->GetUymin(), c_skysigY->GetUxmin(), c_skysigY->GetUymax(), c_skysigY->GetUymin(), c_skysigY->GetUymax(), 505 );
            h1YAxis->SetNdivisions( decAxis->GetNdiv() );
            h1YAxis->SetTitleFont( decAxis->GetLabelFont() );
            h1YAxis->SetLabelFont( decAxis->GetLabelFont() );
            h1YAxis->SetTitleOffset( 1.3 );
            if( sPlot == 2 || sPlot == 3 )
            {
                h1YAxis->SetTitle( "number of events" );
            }
            else if( sPlot == 1 )
            {
                h1YAxis->SetTitle( "number of excess events" );
            }
            h1YAxis->Draw();
            
        }
        
    }
    
    return c_skysig;
}




//////////////////////////////////////////////////////////////////////////////////
/*
 * plot a circle into the sky map for each star of a given catalogue and magnitude
 *
 */

vector<sSource> VPlotAnasumHistograms::plot_catalogue( TCanvas* c, string iCatalogue, double iMaxBrightness, string iBand, double iStarRadius,
        int iColor, int iLineStyle, string hSkyMapName, double iTextAngle, int iMarkerStyle )
{
    // vector with coordinates of objects
    sSource i_sSource;
    vector< sSource > v_obj_XY;
    
    if( !c )
    {
        return v_obj_XY;
    }
    c->cd();
    
    if( iCatalogue.size() < 1 )
    {
        return v_obj_XY;
    }
    
    // define star catalogue
    // (MJD not important since catalogue and plotting coordinates are J2000)
    VStarCatalogue s;
    s.init( 55476., iCatalogue );
    
    double dec = fSkyMapCentreDecJ2000;
    double ra  = fSkyMapCentreRAJ2000;
    
    cout << "Sky map centre used in plot_catalogue: (ra,dec)_J2000 = (" << ra << ", " << dec << ")" << endl;
    
    // get sky map
    TH2D* h = ( TH2D* )c->GetListOfPrimitives()->FindObject( hSkyMapName.c_str() );
    if( !h )
    {
        return v_obj_XY;
    }
    
    // determine extension of sky map
    double iDecMin = fSkyMapCentreDecJ2000 + h->GetYaxis()->GetXmin();
    double iDecMax = fSkyMapCentreDecJ2000 + h->GetYaxis()->GetXmax();
    double iRAMin  = fSkyMapCentreRAJ2000;
    if( cos( dec * TMath::Pi() / 180. ) )
    {
        iRAMin  = fSkyMapCentreRAJ2000 + h->GetXaxis()->GetXmax() / cos( dec * TMath::Pi() / 180. );
    }
    double iRAMax  = fSkyMapCentreRAJ2000;
    if( cos( dec * TMath::Pi() / 180. ) )
    {
        iRAMax  = fSkyMapCentreRAJ2000 + h->GetXaxis()->GetXmin() / cos( dec * TMath::Pi() / 180. );
    }
    
    cout << "MINMAX  " << iDecMin << " " << iDecMax << " " << iRAMin << " " << iRAMax << endl;
    
    for( unsigned int i = 0; i < s.getNStar(); i++ )
    {
        if( s.getStarDec2000( i ) > iDecMin && s.getStarDec2000( i ) < iDecMax )
        {
            if( s.getStarRA2000( i ) < iRAMin && s.getStarRA2000( i ) > iRAMax )
            {
                double x = 0.;
                double y = 0.;
                int j = 0;
                VAstronometry::vlaDs2tp( s.getStarRA2000( i )*TMath::Pi() / 180., s.getStarDec2000( i )*TMath::Pi() / 180., ra * TMath::Pi() / 180., dec * TMath::Pi() / 180., &x, &y, &j );
                x *= -1. * 180. / TMath::Pi();
                y *= 180. / TMath::Pi();
                
                if( s.getStarBrightness( i, iBand ) > 9000 || s.getStarBrightness( i, iBand ) < iMaxBrightness )
                {
                
                    cout << "Object in FOV: ";
                    cout << "#" << i << " (ra,dec)_J2000 " << s.getStarRA2000( i ) << " " << s.getStarDec2000( i );
                    cout << " name=" << s.getStarName( i ) << " mag=" << s.getStarBrightness( i, iBand );
                    cout << ", significance " << h->GetBinContent( h->GetXaxis()->FindBin( x ), h->GetYaxis()->FindBin( y ) );
                    cout << ", (x,y) " << x << " " << y;
                    cout << endl;
                    // return vector with x,y coordinates
                    i_sSource.fX = x;
                    i_sSource.fY = y;
                    i_sSource.fStarName = s.getStarName( i );
                    v_obj_XY.push_back( i_sSource );
                    // draw markers
                    TMarker* m = new TMarker( x, y, TMath::Abs( iMarkerStyle ) );
                    m->SetMarkerColor( iColor );
                    if( s.getStarBrightness( i ) < 0.5 || s.getStarBrightness( i ) > 1000. )
                    {
                        m->Draw();
                    }
                    char hname[50];
                    s.getStarBrightness( i, iBand );
                    if( atoi( s.getStarName( i ).c_str() ) == ( int )s.getStarID( i ) && s.getStarBrightness( i, iBand ) < 9000 )
                    {
                        sprintf( hname, "%.1f", s.getStarBrightness( i, iBand ) );
                    }
                    else
                    {
                        sprintf( hname, "%s", s.getStarName( i ).c_str() );
                    }
                    if( iMarkerStyle > 0 )
                    {
                        TText*  t = new TText( x + 0.125, y + 0.125, hname );
                        t->SetTextAngle( iTextAngle );
                        t->SetTextSize( t->GetTextSize() * 0.5 );
                        t->SetTextColor( iColor );
                        t->Draw();
                        float iB = s.getStarBrightness( i, iBand );
                        if( iB > -50. && iB < 999 )
                        {
                            char hname1[100];
                            sprintf( hname1, "%sMag:%.1f", iBand.c_str(), iB );
                            TText*  t1 = new TText( x + 0.19, y + 0.06, hname1 );
                            t1->SetTextAngle( iTextAngle );
                            t1->SetTextSize( t1->GetTextSize() * 0.5 );
                            t1->SetTextColor( iColor );
                            t1->Draw();
                        }
                    }
                    // draw a circle for extended sources
                    if( s.getStarMajorDiameter( i ) > 0. )
                    {
                        TEllipse* e = new TEllipse( x, y, s.getStarMajorDiameter( i ) );
                        e->SetFillStyle( 0 );
                        e->SetLineColor( iColor );
                        e->SetLineStyle( iLineStyle );
                        e->Draw();
                    }
                    else if( s.getStarBrightness( i, iBand ) > 0. && s.getStarBrightness( i, iBand ) < 900. )
                    {
                        TEllipse* e = new TEllipse( x, y, iStarRadius );
                        e->SetFillStyle( 0 );
                        e->SetLineColor( iColor );
                        e->SetLineStyle( iLineStyle );
                        e->Draw();
                    }
                }
            }
        }
    }
    
    return v_obj_XY;
    
}


////////////////////////////////////////////////////////////////////////////////
/*!
 *   plot ring for ring background model
 *
 *   r:  ring radius
 *   iA: background to source area ratio
 *   t2: theta2
 *   iN: minimum distance of background events from source region
 */
void VPlotAnasumHistograms::plot_RBM_ring( double r, double iA, double t2, double iN )
{
    if( t2 == 0. )
    {
        return;
    }
    // get ring width
    double rw = iA / 4. / TMath::Pi() / r * t2 * 2.;
    rw /= 2.;
    cout << "theta " << TMath::Sqrt( t2 ) << endl;
    cout << "r " << r << endl;
    cout << "r_min " << r - rw << endl;
    cout << "r_max " << r + rw << endl;
    cout << "ring width " << 2.*rw << endl;
    
    // inner ring
    TEllipse* ringI = new TEllipse( 0., 0., r - rw );
    ringI->SetFillStyle( 0 );
    
    // outer ring
    TEllipse* ringO = new TEllipse( 0., 0., r + rw );
    ringO->SetFillStyle( 0 );
    
    // draw exclusion region
    TEllipse* ringE = new TEllipse( 0., 0., iN );
    ringE->SetFillStyle( 0 );
    ringE->SetLineStyle( 2 );
    ringE->SetLineColor( 2 );
    
    // draw theta2 ring
    TEllipse* ringT2 = new TEllipse( 0., 0., TMath::Sqrt( t2 ) );
    ringT2->SetLineStyle( 3 );
    ringT2->SetLineColor( 2 );
    ringT2->SetFillStyle( 0 );
    
    ringI->Draw();
    ringO->Draw();
    ringT2->Draw();
    //    ringE->Draw();
}



///////////////////////////////////////////////////////////////////////////
/*
 *   plot reflected regions from the Refl Reg Background estimation
 *
 *   Can be used both with reflected (ra_dec) and non-reflected (camera coordinates) histograms.
 */
void VPlotAnasumHistograms::plot_reflectedRegions( TCanvas* iC, int i, int j, int iColor )
{
    char itemp[200];
    if( fRunNumber >= 0 )
    {
        sprintf( itemp, "run_%d/stereo/debug/", fRunNumber );
    }
    else
    {
        cout << "plotting of reflected regions only possible for a particular run" << endl;
        cout << "specify with VPlotAnasumHistograms::setRunNumber( int iRun )" << endl;
        return;
    }
    
    TFile* f1 = 0;
    
    if( fAnasumDataFile.size() != 0 )
    {
        f1 = new TFile( fAnasumDataFile.c_str() );
        if( f1->IsZombie() )
        {
            return;
        }
    }
    if( !gDirectory->cd( itemp ) )
    {
        cout << "debug directory not found for run  " << fRunNumber << endl;
        return;
    }
    
    TH2D* hNull = 0;
    
    if( iC == 0 )
    {
        iC = new TCanvas( "iC", "reflected regions", 10, 10, 600, 600 );
        iC->Draw();
        
        hNull = new TH2D( "hNull", "", 100, -1.5, 1.5, 100, -1.5, 1.5 );
        hNull->SetStats( 0 );
        hNull->Draw();
    }
    else
    {
        iC->cd();
    }
    
    // get tree with reflected regions
    TTree* iT = ( TTree* )gDirectory->Get( "tRE" );
    if( !iT )
    {
        cout << "error: reflected region tree not found" << endl;
        return;
    }
    
    int n_r = 0;
    double x = 0.;
    double y = 0.;
    double x_n = 0.;
    double y_n = 0.;
    int x_bin = 0;
    int y_bin = 0;
    int x_bin_wobble = 0;
    int y_bin_wobble = 0;
    double r = 0.;
    double x_r[1000];
    double y_r[1000];
    double r_r[1000];
    int n_ex = 0;
    double x_ex[1000];
    double y_ex[1000];
    double r1_ex[1000];
    double r2_ex[1000];
    double theta_ex[1000];
    
    for( int iLoop = 0; iLoop < 1000; iLoop++ )
    {
        r2_ex[iLoop] = -1;
        theta_ex[iLoop] = 0;
    }
    
    iT->SetBranchAddress( "x_wobble", &x );
    iT->SetBranchAddress( "y_wobble", &y );
    iT->SetBranchAddress( "x", &x_n );
    iT->SetBranchAddress( "y", &y_n );
    iT->SetBranchAddress( "x_bin", &x_bin );
    iT->SetBranchAddress( "y_bin", &y_bin );
    iT->SetBranchAddress( "x_bin_wobble", &x_bin_wobble );
    iT->SetBranchAddress( "y_bin_wobble", &y_bin_wobble );
    iT->SetBranchAddress( "r", &r );
    iT->SetBranchAddress( "n_r", &n_r );
    iT->SetBranchAddress( "x_re", x_r );
    iT->SetBranchAddress( "y_re", y_r );
    iT->SetBranchAddress( "r_re", r_r );
    if( iT->GetBranch( "n_ex" ) )
    {
        iT->SetBranchAddress( "n_ex", &n_ex );
        iT->SetBranchAddress( "x_ex", x_ex );
        iT->SetBranchAddress( "y_ex", y_ex );
    }
    if( iT->GetBranch( "r1_ex" ) )
    {
        iT->SetBranchAddress( "r1_ex", r1_ex );
        iT->SetBranchAddress( "r2_ex", r2_ex );
        iT->SetBranchAddress( "ang_ex", &theta_ex );
    }
    if( iT->GetBranch( "r_ex" ) )
    {
        iT->SetBranchAddress( "r_ex", r1_ex );
    }
    
    
    double iSign = 1.;
    
    //figure out if we've plotted a reflected histogram (with a proper RA axis) or a histogram in derotated camera coordinates.
    //TList::FindObject( char * ) returns 0 if no object with that name is found. No wildcards.
    
    if( iC->GetListOfPrimitives()->FindObject( "hmap_stereo_sig_REFLECTED" ) ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereo_diff_REFLECTED" )
            ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereo_on_REFLECTED" ) ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereo_off_REFLECTED" )
            ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereoUC_sig_REFLECTED" ) ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereoUC_diff_REFLECTED" )
            ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereoUC_on_REFLECTED" ) ||  iC->GetListOfPrimitives()->FindObject( "hmap_stereoUC_off_REFLECTED" )
      )
    {
        iSign *= -1.;
    }
    
    
    
    
    bool bFound = false;
    for( int n = 0; n < iT->GetEntries(); n++ )
    {
        iT->GetEntry( n );
        
        if( i < 0 || ( i == x_bin_wobble && j == y_bin_wobble ) )
        {
            bFound = true;
            break;
        }
    }
    if( !bFound )
    {
        cout << "no reflected regions defined for this bin." << endl;
        return;
    }
    cout << "n_r \t r \t x \t y \t x_bin \t y_bin \t x_bin_wobble \t y_bin_wobble" << endl;
    cout << n_r << "\t" << r << "\t" << x << "\t" << y << "\t" << x_bin << "\t" << y_bin << "\t" << x_bin_wobble << "\t" << y_bin_wobble << endl;
    
    // source region
    TEllipse* iR = new TEllipse( iSign * x, y, r, r );
    iR->SetFillStyle( 3004 );
    iR->SetFillColor( 2 );
    iR->SetLineColor( 2 );
    iR->SetLineWidth( 2 );
    iR->Draw();
    TMarker* iRM = new TMarker( x, y, 5 );
    iRM->Draw();
    
    
    // region around camera center
    TEllipse* iLC = new TEllipse( iSign * ( x - x_n ), y - y_n, r, r );
    iLC->SetFillStyle( 0 );
    iLC->SetLineStyle( 2 );
    iLC->Draw();
    cout << "\t x_r \t y_r \t r_r" << endl;
    for( int n = 0; n < n_r; n++ )
    {
        cout << "\t" << n << "\t" << x_r[n] << "\t" << y_r[n] << "\t" << r_r[n] << endl;
        TEllipse* iL = new TEllipse( iSign * x_r[n], y_r[n], r_r[n], r_r[n] );
        iL->SetFillStyle( 0 );
        iL->SetLineWidth( 2 );
        iL->SetLineColor( iColor );
        iL->Draw();
    }
    // real source region
    TEllipse* iRR = new TEllipse( 0., 0., r, r );
    iRR->SetFillStyle( 0 );
    iRR->SetLineWidth( 2 );
    iRR->SetLineColor( 3 );
    iRR->Draw();
    
    // camera size
    /*   TEllipse *iCC = new TEllipse( x-x_n, y-y_n, 1.2, 1.2 );
       iCC->SetLineStyle( 3 );
       iCC->Draw(); */
    
    // excluded regions
    for( int e = 0; e < n_ex; e++ )
    {
        TEllipse* iEx = new TEllipse( iSign * x_ex[e], y_ex[e], r1_ex[e], r2_ex[e], 0, 360, iSign * theta_ex[e] );
        iEx->SetFillStyle( 0 );
        iEx->SetLineStyle( 3 );
        iEx->SetLineWidth( 2 );
        iEx->SetLineColor( 6 );
        iEx->Draw();
        cout << "\t exclusion regions: " << e << " " << x_ex[e] << " " << y_ex[e] << " " << r1_ex[e] << " " << r2_ex[e] << " " << theta_ex[e] << endl;
        
        if( e == 0 )
        {
            double rr = sqrt( x_ex[e] * x_ex[e] + y_ex[e] * y_ex[e] );
            TEllipse* iRRR = new TEllipse( 0., 0., rr, rr );
            iRRR->SetFillStyle( 0 );
            iRRR->Draw();
        }
    }
    
}


///////////////////////////////////////////////////////////////////////////
/*
 *   plot regions excluded from background calculations
 *
 *
 *
 */
void VPlotAnasumHistograms::plot_excludedRegions( TCanvas* c, int iLineColor )
{
    if( c )
    {
        c->cd();
    }
    char itemp[200];
    sprintf( itemp, "total_1/stereo/" );
    
    TFile* f1 = 0;
    if( fAnasumDataFile.size() != 0 )
    {
        f1 = new TFile( fAnasumDataFile.c_str() );
        if( f1->IsZombie() )
        {
            return;
        }
    }
    if( !gFile->cd( itemp ) )
    {
        cout << "directory " << itemp << " not found" << endl;
        return;
    }
    TTree* t = ( TTree* )gDirectory->Get( "tExcludedRegions" );
    if( !t )
    {
        cout << "tree with excluded regions not found" << endl;
        return;
    }
    float x = 0.;
    float y = 0.;
    float r1 = 0.;
    float r2 = 0.;
    float theta = 0.;
    float Vmag = 0.;
    float Bmag = 0.;
    t->SetBranchAddress( "x", &x );
    t->SetBranchAddress( "y", &y );
    // keep backwards compatibility to circular exclusion regions
    if( t->GetBranchStatus( "r" ) )
    {
        t->SetBranchAddress( "r", &r1 );
        r2 = -1;		//in the TEllipse constructor:  if (r2 <= 0) fR2 = fR1;
        theta = 0.;
    }
    else
    {
        t->SetBranchAddress( "r1", &r1 );
        t->SetBranchAddress( "r2", &r2 );
        t->SetBranchAddress( "theta", &theta );
    }
    t->SetBranchAddress( "Vmag", &Vmag );
    t->SetBranchAddress( "Bmag", &Bmag );
    
    double iSign = 1.;
    
    //figure out if we've plotted a reflected histogram (with a proper RA axis) or a histogram in derotated camera coordinates.
    //TList::FindObject( char * ) returns 0 if no object with that name is found. No wildcards.
    
    if( c->GetListOfPrimitives()->FindObject( "hmap_stereo_sig_REFLECTED" )
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereo_diff_REFLECTED" )
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereo_on_REFLECTED" ) 
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereo_off_REFLECTED" )
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereoUC_sig_REFLECTED" )
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereoUC_diff_REFLECTED" )
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereoUC_on_REFLECTED" ) 
            ||  c->GetListOfPrimitives()->FindObject( "hmap_stereoUC_off_REFLECTED" )
      )
    {
        iSign *= -1.;
    }
    
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        TEllipse* e = new TEllipse( iSign * x, y, r1, r2, 0, 360, iSign * theta );
        e->SetFillStyle( 0 );
        e->SetLineColor( iLineColor );
        e->Draw();
        cout << "#" << i << " Vmag " << Vmag << ", Bmag " << Bmag;
        cout << " ( x = " << x << ", y = " << y << ", r1 = " << r1;
        cout << ", r2 =  " << r2 << ", theta =  " << theta << ")" << endl;
    }
    if( f1 )
    {
        f1->Close();
    }
}


//////////////////////////////////////////////////////////////////////////////
/*!
 *
 * isourceStrength is a source scaling factor (0.1 = scale number of gamma-rays by 0.1)
 *
 * iMethod = 0:  calculate significance (using Li&Ma)
 * iMethod = 1:  calculate q-factor (epsilon_on/sqrt(epsilon_off))
 *
 */
TH1D* VPlotAnasumHistograms::doQfactors( TH1D* hon, TH1D* hoff, TH1D* hdiff, bool bUpper, int iMethod, double iSourceStrength )
{
    if( !hon || !hoff || !hdiff )
    {
        return 0;
    }
    
    hon->Add( hdiff, hoff, iSourceStrength, 1. );
    hdiff->Scale( iSourceStrength );
    
    double i_norm = 1.;
    double i_sumon = 0.;
    double i_sumoff = 0.;
    double i_sumdiff = 0.;
    int i_nbins;
    TH1D* hTemp = 0;
    
    // total number of gamma-ray like events
    double i_ngamma = 0.;
    i_nbins = hdiff->GetNbinsX();
    for( int i = 1; i <= i_nbins; i++ )
    {
        i_ngamma += hdiff->GetBinContent( i );
    }
    // total number of background events
    double i_noff = 0.;
    i_nbins = hoff->GetNbinsX();
    for( int i = 1; i <= i_nbins; i++ )
    {
        i_noff += hoff->GetBinContent( i );
    }
    
    // qfactor on low bound cut
    if( !bUpper )
    {
        hTemp = new TH1D( *( ( TH1D* )hon ) );
        hTemp->Reset();
        if( iMethod == 0 )
        {
            setTitles( hTemp, "qlo", " (Q-Factor for Low-Bound Cut)", "significance" );
        }
        else if( iMethod == 1 )
        {
            setTitles( hTemp, "qlo", " (Q-Factor for Low-Bound Cut)", "q-factor" );
        }
        
        i_sumon  = 0.;
        i_sumoff = 0.;
        i_sumdiff = 0.;
        i_nbins = hTemp->GetNbinsX();
        for( int i = i_nbins; i > 0; i-- )
        {
            i_sumon  += hon->GetBinContent( i );
            i_sumoff += hoff->GetBinContent( i );
            i_sumdiff += hdiff->GetBinContent( i );
            if( iMethod == 0 )
            {
                if( i_sumon + i_sumoff > 0 )
                {
                    hTemp->SetBinContent( i, VStatistics::calcSignificance( i_sumon, i_sumoff, i_norm, 9 ) );
                }
                else
                {
                    hTemp->SetBinContent( i, 0. );
                }
            }
            else if( iMethod == 1 )
            {
                if( i_ngamma > 0. && i_noff > 0. && i_sumoff > 0. )
                {
                    hTemp->SetBinContent( i, ( i_sumdiff / i_ngamma ) / sqrt( i_sumoff / i_noff ) );
                }
                else
                {
                    hTemp->SetBinContent( i, 0. );
                }
            }
            else
            {
                hTemp->SetBinContent( i, 0. );
            }
            hTemp->SetBinError( i, 0. );
        }
    }
    else
    {
        // qfactor on high bound cut
        hTemp = new TH1D( *( ( TH1D* )hon ) );
        hTemp->Reset();
        if( iMethod == 0 )
        {
            setTitles( hTemp, "qhi", " (Q-Factor for High-Bound Cut)", "significance" );
        }
        else if( iMethod == 1 )
        {
            setTitles( hTemp, "qhi", " (Q-Factor for High-Bound Cut)", "q-factor" );
        }
        
        i_sumon  = 0;
        i_sumoff = 0;
        i_sumdiff = 0.;
        i_nbins = hTemp->GetNbinsX();
        for( int i = 1; i <= i_nbins; i++ )
        {
            i_sumon  += hon->GetBinContent( i );
            i_sumoff += hoff->GetBinContent( i );
            i_sumdiff += hdiff->GetBinContent( i );
            if( iMethod == 0 )
            {
                if( i_sumon + i_sumoff > 0 )
                {
                    hTemp->SetBinContent( i, VStatistics::calcSignificance( i_sumon, i_sumoff, i_norm, 9 ) );
                }
                else
                {
                    hTemp->SetBinContent( i, 0. );
                }
            }
            else if( iMethod == 1 )
            {
                if( i_ngamma > 0. && i_noff > 0. && i_sumoff > 0. )
                {
                    hTemp->SetBinContent( i, ( i_sumdiff / i_ngamma ) / sqrt( i_sumoff / i_noff ) );
                }
                else
                {
                    hTemp->SetBinContent( i, 0. );
                }
            }
            else
            {
                hTemp->SetBinContent( i, 0. );
            }
            //	 cout << i << "\t" << i_sumon << "\t" << i_sumoff << "\t" << i_sumdiff << "\t" << i_norm;
            //	 cout << "\t" << calcSignificance(i_sumon,i_sumoff,i_norm,9);
            //	 cout << endl;
            hTemp->SetBinError( i, 0. );
        }
    }
    return hTemp;
}


//////////////////////////////////////////////////////////////////////
/*!
    return histogram with reflected x-axis (new name)
*/
TH2D* VPlotAnasumHistograms::reflectXaxis( TH2D* h, char* iNewName )
{
    if( !h )
    {
        return 0;
    }
    
    char hname[200];
    
    // temporary histogram
    if( iNewName )
    {
        sprintf( hname, "%s", iNewName );
    }
    else
    {
        sprintf( hname, "%s_REFLECTED", h->GetName() );
    }
    
    TH2D* hT = new TH2D( hname, "", h->GetNbinsX(), -1.*h->GetXaxis()->GetXmax(), -1.*h->GetXaxis()->GetXmin(), h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax() );
    hT->SetStats( 0 );
    hT->SetXTitle( h->GetXaxis()->GetTitle() );
    hT->SetYTitle( h->GetYaxis()->GetTitle() );
    
    for( int i = 0; i <= h->GetNbinsX() + 1; i++ )
    {
        for( int j = 0; j <= h->GetNbinsY() + 1; j++ )
        {
            hT->SetBinContent( h->GetNbinsX() + 1 - i, j, h->GetBinContent( i, j ) );
        }
    }
    return hT;
}

/*

    draw size of PSF into a sky map

    all units in [deg]
*/
void VPlotAnasumHistograms::drawPSF( TCanvas* c, string iFile, TH2D* h2, float rPSF_deg )
{
    if( !c )
    {
        return;
    }
    
    c->cd();
    TPad* bigPad = ( TPad* )gPad;
    TDirectory* currrentDir = gDirectory;
    
    // get range of big pad
    double b_xrange = 1. - bigPad->GetRightMargin() - bigPad->GetLeftMargin();
    double b_yrange = 1. - bigPad->GetTopMargin() - bigPad->GetBottomMargin();
    // get range of sky histogram
    double h2_xrange = h2->GetXaxis()->GetBinCenter( h2->GetXaxis()->GetLast() ) - h2->GetXaxis()->GetBinCenter( h2->GetXaxis()->GetFirst() );
    double h2_yrange = h2->GetYaxis()->GetBinCenter( h2->GetYaxis()->GetLast() ) - h2->GetYaxis()->GetBinCenter( h2->GetYaxis()->GetFirst() );
    if( iFile.size() == 0 )
    {
        if( b_xrange > 0. && b_yrange > 0. )
        {
            double x1 = -1.*h2_xrange / 2. + 0.85 * h2_xrange;
            double y1 = -1.*h2_yrange / 2. + 0.12 * h2_yrange;
            
            TEllipse* ePSF = new TEllipse( x1, y1, rPSF_deg, rPSF_deg );
            ePSF->SetLineColor( 1 );
            ePSF->SetLineWidth( 1 );
            ePSF->SetFillStyle( 0 );
            ePSF->Draw();
            
            TText* t = new TText( x1 + rPSF_deg * 1.25, y1 - rPSF_deg / 4., "PSF" );
            t->SetTextSize( 0.04 );
            t->SetTextColor( 1 );
            t->Draw();
        }
        return;
    }
    return;
    /////////////////////////////////////////////////////////////////////////
    // GM: don't know what happens below
    
    TFile f( iFile.c_str() );
    if( f.IsZombie() )
    {
        return;
    }
    
    f.cd( "total_1/stereo/skyHistograms" );
    
    char hname[400];
    sprintf( hname, "%s", h2->GetName() );
    TH2D* h = ( TH2D* )gDirectory->Get( hname );
    if( !h )
    {
        return;
    }
    
    // create a new pad on the lower left corner
    TPad* fPSF = new TPad( "fPSF", "", 0.7, 0.15, 0.85, 0.3 );
    fPSF->SetLeftMargin( 0. );
    fPSF->SetRightMargin( 0. );
    fPSF->SetBottomMargin( 0. );
    fPSF->SetTopMargin( 0. );
    fPSF->SetFrameLineColor( 0 );
    fPSF->SetFrameLineStyle( 1 );
    fPSF->SetFrameLineWidth( 5 );
    fPSF->Draw();
    fPSF->cd();
    
    // get scale for inset
    
    // i) scale maximum
    // (not necessary, since no axis drawn)
    //   if( h->GetMaximum() > 0. ) h->Scale( h2->GetMaximum() / h->GetMaximum() );
    
    // ii) scale x/y axis
    
    // get range of small pad
    double s_xrange = fPSF->GetWNDC();
    double s_yrange = fPSF->GetHNDC();
    // calculate range for small histogram
    double h_xrange = 0.;
    double h_yrange = 0.;
    if( b_xrange > 0. && b_yrange > 0. )
    {
        h_xrange = h2_xrange * s_xrange / b_xrange;
        h_yrange = h2_yrange * s_yrange / b_yrange;
        
        h->SetAxisRange( -1.*h_xrange / 2., h_xrange / 2., "X" );
        h->SetAxisRange( -1.*h_yrange / 2., h_yrange / 2., "Y" );
        h->SetTitle( "" );
        h->SetStats( 0 );
        h->DrawCopy( "A col" );
        gPad->Update();
        // draw a white box around the small pad
        TBox* b = new TBox( -1.*h_xrange / 2., -1.*h_yrange / 2., h_xrange / 2.,  h_yrange / 2. );
        b->SetFillStyle( 0 );
        b->SetLineColor( 0 );
        b->Draw();
    }
    bigPad->cd();
    currrentDir->cd();
}

/*!
    plot dead time histograms and show fit
*/
void VPlotAnasumHistograms::plot_deadTimes()
{
    if( fRunNumber < 0 )
    {
        cout << "VPlotAnasumHistograms::plot_deadTimes(): error no run number set" << endl;
        cout << endl;
        cout << "  this plotting routine works only for individual runs, use VPlotAnasumHistograms::setRunNumber( <run number> ); to select a run" <<  endl;
        return;
    }
    
    TH1D* hTimeDiff_on = ( TH1D* )getHistogram( "hTimeDiff_on", fRunNumber, "deadTimeHistograms" );
    TH1D* hTimeDiffLog_on = ( TH1D* )getHistogram( "hTimeDiffLog_on", fRunNumber, "deadTimeHistograms" );
    TF1* hFTimeDiff_on = ( TF1* )getHistogram( "hFTimeDiff_on", fRunNumber, "deadTimeHistograms" );
    TH1D* hTimeDiff_off = ( TH1D* )getHistogram( "hTimeDiff_off", fRunNumber, "deadTimeHistograms" );
    TH1D* hTimeDiffLog_off = ( TH1D* )getHistogram( "hTimeDiffLog_off", fRunNumber, "deadTimeHistograms" );
    TF1* hFTimeDiff_off = ( TF1* )getHistogram( "hFTimeDiff_off", fRunNumber, "deadTimeHistograms" );
    
    TGraphErrors* hgDeadTime_on = ( TGraphErrors* )getHistogram( "hgDeadTime_on", fRunNumber, "deadTimeHistograms" );
    TGraphErrors* hgDeadTime_off = ( TGraphErrors* )getHistogram( "hgDeadTime_off", fRunNumber, "deadTimeHistograms" );
    
    char hname[200];
    char htitle[200];
    
    if( hTimeDiff_on && hTimeDiffLog_on && hFTimeDiff_on )
    {
        sprintf( hname, "cDTime_on_%d", fRunNumber );
        sprintf( htitle, "dead time histograms, on (run %d)", fRunNumber );
        TCanvas* cDTime = new TCanvas( hname, htitle, 10, 10, 1200, 400 );
        cDTime->Divide( 3, 1 );
        
        cDTime->cd( 1 );
        //	hTimeDiff_on->Rebin( 10 );
        hTimeDiff_on->SetAxisRange( 0., 0.010 );
        if( hTimeDiff_on->GetEntries() > 0 )
        {
            gPad->SetLogy( 1 );
        }
        hTimeDiff_on->SetTitle( "" );
        hTimeDiff_on->GetYaxis()->SetTitleOffset( 1.3 );
        hTimeDiff_on->GetXaxis()->SetNdivisions( 505, true );
        hTimeDiff_on->Draw();
        hFTimeDiff_on->Draw( "same" );
        
        cDTime->cd( 2 );
        hTimeDiffLog_on->GetYaxis()->SetTitleOffset( 1.5 );
        hTimeDiffLog_on->SetStats( 0 );
        hTimeDiffLog_on->SetTitle( "" );
        hTimeDiffLog_on->Draw();
        
        cDTime->cd( 3 );
        TH1D* hTimeDiffFitDiff_on = ( TH1D* )get_ResidualHistogram_from_TF1( "hTimeDiffFitDiff_on", hTimeDiff_on, hFTimeDiff_on );
        if( hTimeDiffFitDiff_on )
        {
            hTimeDiffFitDiff_on->SetStats( 0 );
            hTimeDiffFitDiff_on->SetTitle( "" );
            hTimeDiffFitDiff_on->SetAxisRange( 0., 0.010 );
            hTimeDiffFitDiff_on->SetXTitle( hTimeDiff_on->GetXaxis()->GetTitle() );
            hTimeDiffFitDiff_on->SetYTitle( "fit residuals" );
            hTimeDiffFitDiff_on->SetMinimum( -1.5 );
            hTimeDiffFitDiff_on->SetMaximum( 1.5 );
            hTimeDiffFitDiff_on->Draw();
        }
    }
    if( hTimeDiff_off && hTimeDiffLog_off && hFTimeDiff_off )
    {
        sprintf( hname, "cDTime_off_%d", fRunNumber );
        sprintf( htitle, "dead time histograms, off (run %d)", fRunNumber );
        TCanvas* cDTime = new TCanvas( hname, htitle, 300, 50, 1200, 400 );
        cDTime->Divide( 3, 1 );
        
        cDTime->cd( 1 );
        hTimeDiff_off->SetTitle( "" );
        hTimeDiff_off->SetAxisRange( 0., 0.010 );
        if( hTimeDiff_off->GetEntries() > 0 )
        {
            gPad->SetLogy( 1 );
        }
        hTimeDiff_off->GetXaxis()->SetNdivisions( 505, true );
        hTimeDiff_off->Draw();
        hFTimeDiff_off->Draw( "same" );
        
        cDTime->cd( 2 );
        hTimeDiffLog_off->SetTitle( "" );
        hTimeDiffLog_off->Draw();
        
        cDTime->cd( 3 );
        TH1D* hTimeDiffFitDiff_off = ( TH1D* )get_ResidualHistogram_from_TF1( "hTimeDiffFitDiff_off", hTimeDiff_off, hFTimeDiff_off );
        if( hTimeDiffFitDiff_off )
        {
            hTimeDiffFitDiff_off->SetStats( 0 );
            hTimeDiffFitDiff_off->SetTitle( "" );
            hTimeDiffFitDiff_off->SetAxisRange( 0., 0.010 );
            hTimeDiffFitDiff_off->SetXTitle( hTimeDiff_off->GetXaxis()->GetTitle() );
            hTimeDiffFitDiff_off->SetYTitle( "fit residuals" );
            hTimeDiffFitDiff_off->SetMinimum( -1.5 );
            hTimeDiffFitDiff_off->SetMaximum( 1.5 );
            hTimeDiffFitDiff_off->Draw();
        }
    }
    if( hgDeadTime_on )
    {
        sprintf( hname, "cDTime_g_%d", fRunNumber );
        sprintf( htitle, "dead time vs run time (run %d)", fRunNumber );
        TCanvas* cDTimeG = new TCanvas( hname, htitle, 300, 200, 400, 400 );
        cDTimeG->Draw();
        
        hgDeadTime_on->Draw( "ap" );
        hgDeadTime_on->GetHistogram()->SetXTitle( "run time [min]" );
        hgDeadTime_on->GetHistogram()->SetYTitle( "dead time [%]" );
        if( hgDeadTime_off )
        {
            hgDeadTime_off->SetMarkerStyle( 24 );
            hgDeadTime_off->Draw( "p" );
        }
    }
}

bool VPlotAnasumHistograms::setRunNumber( int iRun )
{
    if( iRun < 0 )
    {
        fRunNumber = iRun;
        return true;
    }
    for( unsigned int i = 0; i < getRunList().size(); i++ )
    {
        if( getRunList()[i].runnumber == iRun )
        {
            fRunNumber = iRun;
            return true;
        }
    }
    cout << "VPlotAnasumHistograms::setRunNumber run number not found! " << endl;
    
    return false;
}

/*

    plot a sky map for each run in the anasum file
    and distributions of 1D width/mean for significances

*/
void VPlotAnasumHistograms::plot_skyPlots_perRun( string iHistoName, double rmax,
        double zmin, double zmax, double rSource,
        int nruns, unsigned int nstart )
{
    // no runs in runlist
    if( getRunList().size() == 0 )
    {
        return;
    }
    
    // useful
    char hname[400];
    char htitle[400];
    
    // total number of entries
    if( nruns < 0 )
    {
        nruns = ( int )getRunList().size();
    }
    if( ( unsigned int )( nruns + nstart ) > getRunList().size() && getRunList().size() > nstart )
    {
        nruns = ( int )( getRunList().size() - nstart );
    }
    
    // plotting mode for 2D plots
    fPlotMode = "colz";
    
    // canvas definition
    int nx = ( int )( sqrt( nruns ) );
    int ny = nruns / ( int )( sqrt( nruns ) );
    if( nx * ny < nruns )
    {
        ny += 1;
    }
    
    sprintf( hname, "cAllMaps_%s", iHistoName.c_str() );
    sprintf( htitle, "all maps (%s)", iHistoName.c_str() );
    TCanvas* allMaps = new TCanvas( hname, htitle, 10, 10, 800, 800 );
    allMaps->Divide( nx, ny );
    allMaps->Draw();
    
    sprintf( hname, "cAllDist_%s", iHistoName.c_str() );
    sprintf( htitle, "all distributions (%s)", iHistoName.c_str() );
    TCanvas* allDist = new TCanvas( hname, htitle, 810, 10, 800, 800 );
    allDist->Divide( nx, ny );
    allDist->Draw();
    gStyle->SetOptFit( 1111 );
    
    // histograms with distribution of fit parameters
    sprintf( hname, "hFit_mean_%s", iHistoName.c_str() );
    TH1D* hFit_mean = new TH1D( hname, "", 100, -0.5, 0.5 );
    hFit_mean->SetXTitle( "mean" );
    hFit_mean->SetYTitle( "# of runs" );
    setHistogramPlottingStyle( hFit_mean, 1, 2, 1 );
    
    sprintf( hname, "hFit_width_%s", iHistoName.c_str() );
    TH1D* hFit_width = new TH1D( hname, "", 100, 0., 2. );
    hFit_width->SetXTitle( "width" );
    hFit_width->SetYTitle( "# of runs" );
    setHistogramPlottingStyle( hFit_width, 1, 2, 1 );
    
    // Gaussian with mean = 0 and RMS = 1
    TF1* fG = new TF1( "fG", "gaus(0)", -5., 5. );
    fG->FixParameter( 1, 0. );
    fG->FixParameter( 2, 1. );
    fG->SetLineColor( 2 );
    fG->SetLineStyle( 3 );
    
    // loop over all runs
    unsigned int z = 0;
    for( unsigned int i = nstart; i < ( unsigned int )nruns + nstart; i++ )
    {
        cout << "now at run " << getRunList()[i].runnumber << endl;
        
        allMaps->cd( z + 1 );
        gPad->SetGridx( 0 );
        gPad->SetGridy( 0 );
        gPad->SetRightMargin( 0.03 );
        gPad->SetLeftMargin( 0.08 );
        gPad->SetTopMargin( 0.05 );
        gPad->SetBottomMargin( 0.05 );
        
        setRunNumber( getRunList()[i].runnumber );
        
        // 2D sky plot
        TH2D* h = ( TH2D* )getHistogram( iHistoName.c_str(), getRunList()[i].runnumber, "skyHistograms" );
        if( !h )
        {
            cout << "histogram not found! " << iHistoName.c_str() << endl;
            continue;
        }
        h->SetStats( 0 );
        h->SetTitle( "" );
        setHistogramPlottingStyle( h, 1.5 );
        if( zmax > -90. )
        {
            h->SetMaximum( zmax );
        }
        if( zmin > -90. )
        {
            h->SetMinimum( zmin );
        }
        h->Draw( fPlotMode.c_str() );
        
        TMarker* iMZ = new TMarker( 0., 0., 5 );
        iMZ->SetMarkerColor( 5 );
        iMZ->Draw();
        
        sprintf( hname, "%d", getRunList()[i].runnumber );
        TText* iT = new TText( 0.7, 0.8, hname );
        iT->SetTextSize( 0.10 );
        iT->SetNDC();
        iT->Draw();
        
        gPad->Update();
        
        //////////////////////////////////////////////////////////////////
        // significance distribution
        allDist->cd( z + 1 );
        gPad->SetGridx( 0 );
        gPad->SetGridy( 0 );
        gPad->SetRightMargin( 0.03 );
        gPad->SetLeftMargin( 0.08 );
        gPad->SetTopMargin( 0.05 );
        gPad->SetBottomMargin( 0.05 );
        
        TH1D* hsig_1D = plot_significanceDistributions( rmax, rSource, -5.5, 10., ( TCanvas* )gPad, false );
        
        if( hsig_1D && hsig_1D->GetEntries() > 0 )
        {
            hsig_1D->Fit( fG, "Q" );
            hsig_1D->Fit( "gaus", "Q" );
            cout << "RUN " << getRunList()[i].runnumber;
            if( hsig_1D->GetFunction( "gaus" ) )
            {
                hsig_1D->GetFunction( "gaus" )->SetLineColor( 8 );
                hsig_1D->GetFunction( "gaus" )->SetLineStyle( 2 );
                hsig_1D->GetFunction( "gaus" )->Draw( "same" );
                fG->Draw( "same" );
                cout << " ,fit results: mean " << setprecision( 3 ) << hsig_1D->GetFunction( "gaus" )->GetParameter( 1 );
                cout << " +- " << setprecision( 4 ) << hsig_1D->GetFunction( "gaus" )->GetParError( 1 );
                cout << ", RMS " << setprecision( 3 ) << hsig_1D->GetFunction( "gaus" )->GetParameter( 2 );
                cout << " +- " << setprecision( 4 ) <<  hsig_1D->GetFunction( "gaus" )->GetParError( 2 );
                cout << ", probability " << setprecision( 3 ) << hsig_1D->GetFunction( "gaus" )->GetProb();
                cout << ", Chi2/N " << setprecision( 3 ) << hsig_1D->GetFunction( "gaus" )->GetChisquare() / hsig_1D->GetFunction( "gaus" )->GetNDF();
                cout << endl;
            }
            else
            {
                cout << ", error retrieving fit function" << endl;
            }
            if( hsig_1D->GetFunction( "gaus" ) )
            {
                TLine* iM = new TLine( hsig_1D->GetFunction( "gaus" )->GetParameter( 1 ), hsig_1D->GetYaxis()->GetXmin(),
                                       hsig_1D->GetFunction( "gaus" )->GetParameter( 1 ), hsig_1D->GetYaxis()->GetXmax() );
                iM->SetLineStyle( 2 );
                iM->SetLineColor( 3 );
                iM->Draw();
            }
            TLine* iL = new TLine( 0., hsig_1D->GetYaxis()->GetXmin(), 0., hsig_1D->GetMaximum() );
            iL->SetLineStyle( 2 );
            iL->Draw();
        }
        
        // draw run number
        iT->Draw();
        
        gPad->Update();
        
        // fill histogram with fit parameters
        if( hsig_1D && hsig_1D->GetEntries() > 0 && hsig_1D->GetFunction( "gaus" ) )
        {
            hFit_mean->Fill( hsig_1D->GetFunction( "gaus" )->GetParameter( 1 ) );
            hFit_width->Fill( hsig_1D->GetFunction( "gaus" )->GetParameter( 2 ) );
        }
        // counter
        z++;
    }
    
    // plot 1D distributions for Gaus fit to all runs
    sprintf( hname, "cFit_mean_%s", iHistoName.c_str() );
    sprintf( htitle, "mean value of fit (%s)", iHistoName.c_str() );
    TCanvas* cFit_mean = new TCanvas( hname, htitle, 10, 10, 400, 400 );
    cFit_mean->SetGridx( 0 );
    cFit_mean->SetGridy( 0 );
    cFit_mean->Draw();
    
    hFit_mean->SetStats( 1 );
    hFit_mean->Draw();
    TLine* lFitMean = new TLine( 0., 0., 0., hFit_mean->GetMaximum() );
    lFitMean->SetLineStyle( 2 );
    lFitMean->Draw();
    
    sprintf( hname, "cFit_width_%s", iHistoName.c_str() );
    sprintf( htitle, "width value of fit (%s)", iHistoName.c_str() );
    TCanvas* cFit_width = new TCanvas( hname, htitle, 410, 10, 400, 400 );
    cFit_width->SetGridx( 0 );
    cFit_width->SetGridy( 0 );
    cFit_width->Draw();
    hFit_width->SetStats( 1 );
    hFit_width->Draw();
    TLine* lFitRMS = new TLine( 1., 0., 1., hFit_width->GetMaximum() );
    lFitRMS->SetLineStyle( 2 );
    lFitRMS->Draw();
    
    // reset run number
    setRunNumber( -1 );
}

///////////////////////////////////////////////////////////////////////////////////
/*
 *
 *
 */
void VPlotAnasumHistograms::plot_theta2Correction()
{
    if( fRunNumber < 0 )
    {
        cout << "VPlotAnasumHistograms::plot_theta2Correction(): error no run number set" << endl;
        cout << endl;
        cout << "  this plotting routine works only for individual runs, use VPlotAnasumHistograms::setRunNumber( <run number> ); to select a run" <<  endl;
        return;
    }
    char itemp[200];
    TH1D* hon = ( TH1D* )getHistogram( "hAux_theta2On", fRunNumber, "debug" );
    TH1D* hoff = ( TH1D* )getHistogram( "hAux_theta2Off", fRunNumber, "debug" );
    TH1D* hrat = ( TH1D* )getHistogram( "hAux_theta2Ratio", fRunNumber, "debug" );
    
    if( hon && hoff && hrat )
    {
        char ctitle[200];
        sprintf( itemp, "c_t2C_%d", fRunNumber );
        sprintf( ctitle, "theta2 correction (%s,%d)", fAnasumDataFile.c_str(), fRunNumber );
        TCanvas* c_t2C = new TCanvas( itemp, ctitle, 10, 10, 800, 400 );
        c_t2C->Divide( 2, 1 );
        
        setHistogramPlottingStyle( hon, 1, 2, 2, 20, 1 );
        setHistogramPlottingStyle( hoff, 2, 2, 2, 20, 1 );
        setHistogramPlottingStyle( hrat, 1, 2, 2, 20, 1 );
        cout << "On: black histogram" << endl;
        cout << "Off:  red histogram" << endl;
        
        c_t2C->cd( 1 );
        hon->SetMinimum( 1 );
        hon->Draw();
        hoff->Draw( "same" );
        
        c_t2C->cd( 2 );
        hrat->Draw();
        TLine* iL = new TLine( 0., 1., hrat->GetXaxis()->GetXmax(), 1. );
        iL->SetLineWidth( 2 );
        iL->SetLineStyle( 2 );
        iL->Draw();
    }
}


/*

    GM: obsolete, use VEnergySpectrum


void VPlotAnasumHistograms::fit_energy(double minE, double maxE)
{

    char hname[200];
    char htitle[200];

    TH1D *h = (TH1D*)getHistogram( "herec_diff", fRunNumber, "energyHistograms" );
    if( !h ) return;
    setHistogramPlottingStyle( h, 1, 3, 2, 20, 1 );
    h->SetMarkerStyle( 20 );
    h->SetTitle( "" );
    h->SetStats( 1 );

    sprintf( hname, "cg_efit_%d_%s", fRunNumber, fAnasumDataFile.c_str() );
    sprintf( htitle, "energy spectra (fit, %d, %s)", fRunNumber, fAnasumDataFile.c_str() );
    TCanvas *cg_efit = new TCanvas( hname, htitle, 10, 10, 400, 400 );
    cg_efit->Draw();
    gPad->SetLeftMargin( 0.12 );
    if( h->GetEntries() > 0 ) cg_efit->SetLogy( 1 );

    TF1 *e = new TF1( "e", "[0]*TMath::Power(TMath::Power(10,x),[1])", minE, maxE );
    e->SetParameter( 0, 1.e-7 );
    e->SetParameter( 1, -2. );

    h->SetAxisRange( minE, maxE );
    h->GetYaxis()->SetTitleOffset( 1.5 );

    h->Fit( "e", "RE" );

    cout << e->GetChisquare() << "\t" << e->GetNDF() << "\t";
    if( e->GetNDF() > 0 ) cout << e->GetChisquare()/ e->GetNDF();
    cout << endl;

    h->Draw();
} */

/////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Trigger pattern
 *
 */
TH1D* VPlotAnasumHistograms::plot_triggerpattern( int ntel, bool bPlot )
{

    TH1D* hTriggerPatternBeforeCuts_on = ( TH1D* )getHistogram( "hTriggerPatternBeforeCuts_on", fRunNumber, "rawRateHistograms" );
    TH1D* hTriggerPatternAfterCuts_on = ( TH1D* )getHistogram( "hTriggerPatternAfterCuts_on", fRunNumber, "rawRateHistograms" );
    //    TH1D *hImagePatternBeforeCuts_on = (TH1D*)getHistogram( "hImagePatternBeforeCuts_on", fRunNumber, "rawRateHistograms" );
    //    TH1D *hImagePatternAfterCuts_on = (TH1D*)getHistogram( "hImagePatternAfterCuts_on", fRunNumber, "rawRateHistograms" );
    
    char hname[200];
    char htitle[200];
    if( hTriggerPatternBeforeCuts_on && hTriggerPatternAfterCuts_on )
    {
    
        hTriggerPatternBeforeCuts_on->SetAxisRange( 1, TMath::Power( 2, ntel ) - 1 );
        hTriggerPatternBeforeCuts_on->SetStats( 0 );
        hTriggerPatternBeforeCuts_on->SetTitle( "" );
        hTriggerPatternBeforeCuts_on->SetLineWidth( 2 );
        hTriggerPatternAfterCuts_on->SetStats( 0 );
        hTriggerPatternAfterCuts_on->SetLineWidth( 2 );
        hTriggerPatternAfterCuts_on->SetLineColor( 2 );
        hTriggerPatternAfterCuts_on->SetTitle( "" );
        hTriggerPatternBeforeCuts_on->GetYaxis()->SetTitleOffset( 1.3 );
        
        if( bPlot )
        {
            sprintf( hname, "c_TPattern_%d", fRunNumber );
            sprintf( htitle, "trigger pattern ( run %d)", fRunNumber );
            TCanvas* c_TPattern_ = new TCanvas( hname, htitle, 10, 10, 600, 400 );
            c_TPattern_->Draw();
            
            if( hTriggerPatternBeforeCuts_on->GetEntries() > 0. )
            {
                hTriggerPatternBeforeCuts_on->Scale( 1. / hTriggerPatternBeforeCuts_on->GetEntries() );
            }
            if( hTriggerPatternAfterCuts_on->GetEntries() > 0. )
            {
                hTriggerPatternAfterCuts_on->Scale( 1. / hTriggerPatternAfterCuts_on->GetEntries() );
            }
            hTriggerPatternBeforeCuts_on->SetYTitle( "fraction of events" );
            hTriggerPatternBeforeCuts_on->SetMaximum( 1. );
            hTriggerPatternBeforeCuts_on->Draw( "histo" );
            hTriggerPatternAfterCuts_on->Draw( "histo same" );
            
            TLegend* iL = new TLegend( 0.13, 0.7, 0.4, 0.85 );
            iL->AddEntry( hTriggerPatternBeforeCuts_on, "before cuts", "l" );
            iL->AddEntry( hTriggerPatternAfterCuts_on, "after cuts", "l" );
            iL->Draw();
        }
        
    }
    
    return hTriggerPatternAfterCuts_on;
}


///////////////////////////////////////////////////////////////////////////////
/*
 *  rmax > 0: apply a circular cut around pointing direction
 *  rmax < 0: apply box cut around pointing direction
 */

TCanvas* VPlotAnasumHistograms::plot_skyPlots_significance( bool iCorrelated, double rmax, double zmin, double zmax, bool bPoster )
{
    bool bZetaTau = false;
    
    fPlotMode = "colz";
    
    char hname[600];
    char htitle[600];
    
    // get all histograms
    TH2D* hmap_stereo_sig = 0;
    TH2D* hmap_stereo_diff = 0;
    TH2D* hmap_alphaNorm_off = 0;
    
    if( iCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_diff" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_diff" );
    }
    hmap_stereo_diff = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    if( rmax > 0. )
    {
        hmap_stereo_diff = removeOuterRing( hmap_stereo_diff, rmax, -9999. );
    }
    setHistogramPlottingStyle( hmap_stereo_diff, 1.5 );
    
    if( iCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_sig" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_sig" );
    }
    hmap_stereo_sig = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    if( rmax > 0. )
    {
        hmap_stereo_sig = removeOuterRing( hmap_stereo_sig, rmax, -9999. );
    }
    setHistogramPlottingStyle( hmap_stereo_sig, 1.5 );
    
    if( iCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_sig" );
    }
    else
    {
        sprintf( hname, "hmap_alphaNorm_off" );
    }
    hmap_alphaNorm_off = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    if( rmax > 0. )
    {
        hmap_alphaNorm_off = removeOuterRing( hmap_alphaNorm_off, rmax, -9999. );
    }
    setHistogramPlottingStyle( hmap_alphaNorm_off, 1.5 );
    
    // ON-OFF
    if( hmap_stereo_diff )
    {
        if( iCorrelated )
        {
            sprintf( hname, "c_skyONFFUC_%s_%d", fAnasumDataFile.c_str(), fRunNumber );
            sprintf( htitle, "sky map (ON-OFF)  run %d (uncorrelated plots)", fRunNumber );
        }
        else
        {
            sprintf( hname, "c_skyONOFF_%s_%d", fAnasumDataFile.c_str(), fRunNumber );
            sprintf( htitle, "sky/camera map (ON-OFF) run %d (correlated plots)", fRunNumber );
        }
        TCanvas* c_skyOnOff = new TCanvas( hname, htitle, 10, 10, 400, 400 );
        c_skyOnOff->Draw();
        gPad->SetRightMargin( 0.12 );
        
        if( rmax != 0. )
        {
            hmap_stereo_diff->SetAxisRange( -1.*fabs( rmax ), fabs( rmax ), "X" );
            hmap_stereo_diff->SetAxisRange( -1.*fabs( rmax ), fabs( rmax ), "Y" );
        }
        
        hmap_stereo_diff->Draw( fPlotMode.c_str() );
    }
    
    // significance
    TCanvas* c_skySig = 0;
    if( hmap_stereo_sig )
    {
        if( iCorrelated )
        {
            sprintf( hname, "c_skySigUC_%s_%d", fAnasumDataFile.c_str(), fRunNumber );
            sprintf( htitle, "sky map (Significance)  run %d (uncorrelated plots)", fRunNumber );
        }
        else
        {
            sprintf( hname, "c_skySig_%s_%d", fAnasumDataFile.c_str(), fRunNumber );
            sprintf( htitle, "sky/camera map (Significance) run %d (correlated plots)", fRunNumber );
        }
        c_skySig = new TCanvas( hname, htitle, 430, 10, 600, 600 );
        c_skySig->Draw();
        gPad->SetRightMargin( 0.12 );
        
        if( zmax > -90. )
        {
            hmap_stereo_sig->SetMaximum( zmax );
        }
        if( zmin > -90. )
        {
            hmap_stereo_sig->SetMinimum( zmin );
        }
        if( rmax != 0. )
        {
            hmap_stereo_sig->SetAxisRange( -1.*fabs( rmax ), fabs( rmax ), "X" );
            hmap_stereo_sig->SetAxisRange( -1.*fabs( rmax ), fabs( rmax ), "Y" );
        }
        
        if( bPoster )
        {
            hmap_stereo_sig->SetTitle( "" );
            hmap_stereo_sig->SetZTitle( "significance" );
            hmap_stereo_sig->GetZaxis()->SetTitleOffset( 0.8 );
        }
        
        hmap_stereo_sig->Draw( fPlotMode.c_str() );
        if( bPoster )
        {
            drawPSF( 0, "", hmap_stereo_sig, 0.11 );
        }
        
        // plot ZetaTau?
        if( bZetaTau )
        {
            TGraph* gZetaTau = new TGraph( 1 );
            gZetaTau->SetPoint( 0, 0.7263087 , 0.8719 );
            gZetaTau->SetMarkerSize( 2 );
            gZetaTau->SetMarkerStyle( 24 );
            gZetaTau->SetMarkerColor( 5 );
            gZetaTau->Draw( "p" );
        }
    }
    TCanvas* c_skyAlphaNorm = 0;
    if( hmap_alphaNorm_off )
    {
        if( iCorrelated )
        {
            sprintf( hname, "c_skyAlphaNormUC_%s_%d", fAnasumDataFile.c_str(), fRunNumber );
            sprintf( htitle, "sky map (norm)  run %d (uncorrelated plots)", fRunNumber );
        }
        else
        {
            sprintf( hname, "c_skyAlphaNorm_%s_%d", fAnasumDataFile.c_str(), fRunNumber );
            sprintf( htitle, "sky/camera map (norm) run %d (correlated plots)", fRunNumber );
        }
        c_skyAlphaNorm = new TCanvas( hname, htitle, 860, 10, 400, 400 );
        c_skyAlphaNorm->Draw();
        gPad->SetRightMargin( 0.14 );
        
        hmap_alphaNorm_off->SetMinimum( 0. );
        hmap_alphaNorm_off->SetMaximum( hmap_alphaNorm_off->GetMaximum() * 1.1 );
        if( rmax != 0. )
        {
            hmap_alphaNorm_off->SetAxisRange( -1.*fabs( rmax ), fabs( rmax ), "X" );
            hmap_alphaNorm_off->SetAxisRange( -1.*fabs( rmax ), fabs( rmax ), "X" );
        }
        
        hmap_alphaNorm_off->Draw( fPlotMode.c_str() );
    }
    
    return c_skySig;
}


//////////////////////////////////////////////////////////////////////////////////
/*
 * plot_on
 *
 */
TCanvas* VPlotAnasumHistograms::plot_on( bool iCorrelated, double rmax )
{

    fPlotMode = "colz";
    
    char hname[200];
    char htitle[200];
    
    // get all histograms
    TH2D* hmap_stereo_on = 0;
    
    if( iCorrelated )
    {
        sprintf( hname, "hmap_stereoUC_on" );
    }
    else
    {
        sprintf( hname, "hmap_stereo_on" );
    }
    hmap_stereo_on = ( TH2D* )getHistogram( hname, fRunNumber, "skyHistograms" );
    if( rmax > 0. )
    {
        hmap_stereo_on = removeOuterRing( hmap_stereo_on, rmax, -9999. );
    }
    setHistogramPlottingStyle( hmap_stereo_on, 1.5 );
    
    TCanvas* c_skyOnOff = 0;
    
    if( hmap_stereo_on )
    {
        if( iCorrelated )
        {
            sprintf( hname, "c_skyONFFUC_%d", fRunNumber );
            sprintf( htitle, "sky map (ON-OFF)  run %d (uncorrelated plots)", fRunNumber );
        }
        else
        {
            sprintf( hname, "c_skyONOFF_%d", fRunNumber );
            sprintf( htitle, "sky/camera map (ON-OFF) run %d (correlated plots)", fRunNumber );
        }
        c_skyOnOff = new TCanvas( hname, htitle, 10, 10, 400, 400 );
        c_skyOnOff->Draw();
        gPad->SetRightMargin( 0.12 );
        
        hmap_stereo_on->SetAxisRange( -1.6, 1.6, "X" );
        hmap_stereo_on->SetAxisRange( -1.6, 1.6, "Y" );
        hmap_stereo_on->SetTitle( "" );
        
        hmap_stereo_on->Draw( fPlotMode.c_str() );
    }
    return c_skyOnOff;
}



