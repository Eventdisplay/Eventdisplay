/*! \class VSourceGeometryFitter
    \brief analyse source geometry (position and extension)


    TODO:

    plotter
    LL function
    arcsec calculation
    look at extension

*/

#include "VSourceGeometryFitter.h"

VSourceGeometryFitter::VSourceGeometryFitter()
{
    fDebug = false;
    
    fAnasumDataFile = "";
    fRunNumber      = -1;
    fHisSkyMap      = 0;
    fXStart         = 0;
    fYStart         = 0;
    fPSF            = 0.063;
    setFitterDefaultData();
    setFitter( "RadialAsymmetricSource_LL" );
    
}

VSourceGeometryFitter::VSourceGeometryFitter( string iAnaSumDataFile, int iRunNumber )
{
    fDebug = false;
    
    fAnasumDataFile = iAnaSumDataFile;
    fRunNumber      = iRunNumber;
    fHisSkyMap      = 0;
    fXStart         = 0;
    fYStart         = 0;
    fPSF            = 0.063;
    if( !openFile( fAnasumDataFile, fRunNumber, 1 ) )
    {
        return;
    }
    
    setFitterDefaultData();
    setFitter( "RadialAsymmetricSource_LL" );
    
}

void VSourceGeometryFitter::setFitterDefaultData()
{

    fDefaultFitterData.clear();
    
    
    ////////
    // SRC 1: radial symmetric 2D source
    ////////
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "RadialSymmetricSource_Chi2";
    fDefaultFitterData.back()->fFitterDescription = "Radial Symmetric 2D gaussian convolved with gaussian PSF, chisquare fit";
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaSRC" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "constant" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 1.e4 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 0. );
    
    
    
    
    ////////
    // SRC 2: radial asymmetric 2D source - does not work at the moment!!!
    ////////
    /* fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "RadialAsymmetricSource_Chi2";
    fDefaultFitterData.back()->fFitterDescription = "Radial Asymmetric 2D gaussian convolved with gaussian PSF, chisquare fit";
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back(  1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back(  1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaSRC_X" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaSRC_Y" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "theta" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.0 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1.e3 );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e3 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "constant" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 1.e5 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e5 );*/
    
    /*
    // 2D Normal distribution (LL)
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "2DAsymGauss_LL";
    fDefaultFitterData.back()->fFitterDescription = "Asymmetrical Gaussian";
    
    fDefaultFitterData.back()->fParameterName.push_back( "rho" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back(  0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back(  1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back(  1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaX" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back(  1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaY" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    */
    
    
    ////////
    // SRC 3: radial symmetric 2D source, LL
    ////////
    // 2D Normal distribution,convolved with PSF (LL)
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "RadialSymmetricSource_LL";
    fDefaultFitterData.back()->fFitterDescription = "Bi-variate normal function with rho=0 (correlation coefficient) convolved with PSF";
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigma" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    
    ////////
    // SRC 4: radial asymmetric 2D source
    ////////
    // 2D Normal distribution,convolved with PSF (LL)
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "RadialAsymmetricSource_LL";
    fDefaultFitterData.back()->fFitterDescription = "Bi-variate normal function convolved with PSF";
    
    //   fDefaultFitterData.back()->fParameterName.push_back( "rho" );
    fDefaultFitterData.back()->fParameterName.push_back( "angle" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaX" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaY" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    
    ////////
    // PSF 1: radial symmetric 2D gaussian with an offset (Chi2)
    ////////
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "2DGauss_Chi2";
    fDefaultFitterData.back()->fFitterDescription = "Radial symmetric 2D gaussian with an offset, chisquare fit";
    
    fDefaultFitterData.back()->fParameterName.push_back( "background" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 30. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1.e5 );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e5 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "constant" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 1.e5 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 0. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigma" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 1.e-3 );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.0 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    
    ////////
    // PSF 2: radial symmetric 2D gaussian (LL)
    ////////
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "2DGauss_LL";
    fDefaultFitterData.back()->fFitterDescription = "Radial Symmetric 2D gaussian, binned log-likelihood fit";
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigma" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    
    ////////
    // PSF 3: linear superposition of two gaussians (LL)
    ////////
    fDefaultFitterData.push_back( new VSourceGeometryFitterData() );
    fDefaultFitterData.back()->fFitterName = "LinearSuperposition2DGauss_LL";
    fDefaultFitterData.back()->fFitterDescription = "Superposition of two radial symmetric 2D gaussians (central spot + broad halo)";
    
    fDefaultFitterData.back()->fParameterName.push_back( "Xpos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "Ypos" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0. );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( -1. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigma" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.1 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "sigmaH" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.2 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1.e2 );
    
    fDefaultFitterData.back()->fParameterName.push_back( "alpha" );
    fDefaultFitterData.back()->fParameterInitValue.push_back( 0.9 );
    fDefaultFitterData.back()->fParameterStep.push_back( 1.e-7 );
    fDefaultFitterData.back()->fParameterLowerLimit.push_back( 0. );
    fDefaultFitterData.back()->fParameterUpperLimit.push_back( 1. );
    
    
}

void VSourceGeometryFitter::help()
{

    cout << "list of available fit functions (select with VSourceGeometryFitter::setFitter() ): " << endl;
    for( unsigned int i = 0; i < fDefaultFitterData.size(); i++ )
    {
        cout << "\t" << fDefaultFitterData[i]->fFitterName << "   " << fDefaultFitterData[i]->fFitterDescription << endl;
    }
}

TCanvas* VSourceGeometryFitter::plot( double rmax, double zmin, double zmax, string iPlotMode )
{
    if( !fHisSkyMap )
    {
        cout << "\t no sky map set" << endl;
        return 0;
    }
    // preliminary
    double xcenter = -1 * fXStart;
    double ycenter = -1 * fYStart;
    
    // style settings
    default_settings();
    setHistogramPlottingStyle( fHisSkyMap, 1.5 );
    
    // signi
    char hname[800];
    char htitle[800];
    sprintf( hname, "c_skymap_%d", fRunNumber );
    sprintf( htitle, "sky map, run %d", fRunNumber );
    TCanvas* c_sky = new TCanvas( hname, htitle, 10, 10, 400, 400 );
    c_sky->Draw();
    c_sky->SetRightMargin( 0.14 );
    c_sky->SetLeftMargin( 0.11 );
    c_sky->SetGridx( 1 );
    c_sky->SetGridy( 1 );
    
    fHisSkyMap->SetTitle( "" );
    
    double x1 = -1.*rmax - xcenter;
    double x2 = rmax - xcenter;
    double y1 = -1.*rmax - ycenter;
    double y2 = rmax - ycenter;
    if( x1 < fHisSkyMap->GetXaxis()->GetXmin() )
    {
        x1 = fHisSkyMap->GetXaxis()->GetXmin();
    }
    if( x2 > fHisSkyMap->GetXaxis()->GetXmax() )
    {
        x2 = fHisSkyMap->GetXaxis()->GetXmax();
    }
    if( y1 < fHisSkyMap->GetYaxis()->GetXmin() )
    {
        y1 = fHisSkyMap->GetYaxis()->GetXmin();
    }
    if( y2 > fHisSkyMap->GetYaxis()->GetXmax() )
    {
        y2 = fHisSkyMap->GetYaxis()->GetXmax();
    }
    fHisSkyMap->SetAxisRange( x1, x2, "X" );
    fHisSkyMap->SetAxisRange( y1, y2, "Y" );
    
    if( zmin > -1000. )
    {
        fHisSkyMap->SetMinimum( zmin );
    }
    if( zmax > -1000. )
    {
        fHisSkyMap->SetMaximum( zmax );
    }
    
    fHisSkyMap->Draw( iPlotMode.data() );
    
    // plot reconstructed source geometry
    plotSourceGeometry();
    
    // plot fit result
    plotFitResult();
    
    // return canvas
    return c_sky;
}

void VSourceGeometryFitter::plotFitResult()
{
    if( !fFitter || fFitter->fFitResult_Status < -10 )
    {
        return;
    }
    // check right number of parameters
    if( fFitter->fFitResult_Parameter.size() != fFitter->fParameterName.size() || fFitter->fFitResult_ParameterError.size() != fFitter->fParameterName.size() )
    {
        return;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // symmetric Gauss fit
    if( fFitter->fFitterName == "2DGauss_Chi2" )
    {
        if( fHisSkyMap )
        {
            char hname[800];
            char htitle[800];
            sprintf( hname, "c_skyFitResultmap_%d", fRunNumber );
            sprintf( htitle, "fit result, run %d", fRunNumber );
            TCanvas* c_skyFitResult = new TCanvas( hname, htitle, 610, 10, 400, 400 );
            c_skyFitResult->Draw();
            c_skyFitResult->SetRightMargin( 0.14 );
            c_skyFitResult->SetLeftMargin( 0.11 );
            c_skyFitResult->SetGridx( 1 );
            c_skyFitResult->SetGridy( 1 );
            
            sprintf( hname, "hFitResult_%d", fRunNumber );
            TH1D* hFitResult = new TH1D( hname, "", 4 * ( int )( sqrt( fFitter->fParameterUpperLimit[3]*fFitter->fParameterUpperLimit[3] + fFitter->fParameterUpperLimit[4]*fFitter->fParameterUpperLimit[4] ) / 0.025 / 0.025 ), 0., sqrt( fFitter->fParameterUpperLimit[3]*fFitter->fParameterUpperLimit[3] + fFitter->fParameterUpperLimit[4]*fFitter->fParameterUpperLimit[4] ) );
            hFitResult->SetXTitle( "#Theta^{2} [deg]" );
            hFitResult->SetYTitle( "dN/d#Theta^{2}" );
            setHistogramPlottingStyle( hFitResult, 1, 2. );
            
            // fill theta2 histogram
            
            double x = 0.;
            double y = 0.;
            double t2 = 0.;
            int nbinsX = fHisSkyMap->GetNbinsX();
            int nbinsY = fHisSkyMap->GetNbinsY();
            for( int i = 1; i <= nbinsX; i++ )
            {
                x = fHisSkyMap->GetXaxis()->GetBinCenter( i );
                for( int j = 1; j <= nbinsY; j++ )
                {
                    y = fHisSkyMap->GetYaxis()->GetBinCenter( j );
                    // calculate theta2
                    t2 = ( x - fFitter->fFitResult_Parameter[3] ) * ( x - fFitter->fFitResult_Parameter[3] ) + ( y - fFitter->fFitResult_Parameter[4] ) * ( y - fFitter->fFitResult_Parameter[4] );
                    // fill histogram
                    if( fHisSkyMap->GetBinContent( i, j ) > -90. )
                    {
                        hFitResult->Fill( t2, fHisSkyMap->GetBinContent( i, j ) );
                    }
                }
            }
            // plot histogram
            hFitResult->Draw();
            // plot fit function
            sprintf( hname, "fFitResult_%d", fRunNumber );
            sprintf( htitle, "%f + %f*TMath::Exp( -1.*x/ 2. / %f / %f)", fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[2] );
            TF1* fFitResult = new TF1( hname, htitle, hFitResult->GetXaxis()->GetXmin(), hFitResult->GetXaxis()->GetXmax() );
            fFitResult->SetLineColor( 2 );
            fFitResult->Draw( "same" );
        }
        
    }
}



TGraph* VSourceGeometryFitter::plotSourceGeometry( int iColor )
{
    if( !fFitter || fFitter->fFitResult_Status < -10 )
    {
        return 0;
    }
    // check right number of parameters
    if( fFitter->fFitResult_Parameter.size() != fFitter->fParameterName.size() || fFitter->fFitResult_ParameterError.size() != fFitter->fParameterName.size() )
    {
        return 0;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // PSF #1: radial symmetric Gauss fit + constant, Chi2
    if( fFitter->fFitterName == "2DGauss_Chi2" )
    {
        // define graph
        TGraphErrors* g = new TGraphErrors( 1 );
        setGraphPlottingStyle( g, iColor, 2., 7, 2. );
        
        g->SetPoint( 0, fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[4] );
        g->SetPointError( 0, fFitter->fFitResult_ParameterError[3], fFitter->fFitResult_ParameterError[4] );
        g->Draw( "p" );
        
        if( fFitter->fFitResult_Parameter[2] > 0. )
        {
            TEllipse* e = new TEllipse( fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[4], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[2] );
            e->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e->Draw();
        }
        
        return g;
    }
    /////////////////////////////////////////
    // PSF #2: radial symmetric Gauss fit, LL
    else if( fFitter->fFitterName == "2DGauss_LL" )
    {
    
        // define graph
        TGraphErrors* g = new TGraphErrors( 1 );
        setGraphPlottingStyle( g, iColor, 2., 7, 2. );
        
        g->SetPoint( 0, fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1] );
        g->SetPointError( 0, fFitter->fFitResult_ParameterError[0], fFitter->fFitResult_ParameterError[1] );
        g->Draw( "p" );
        
        
        if( fFitter->fFitResult_Parameter[2] > 0. )
        {
            TEllipse* e = new TEllipse( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[2] );
            e->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e->Draw();
        }
        
        return g;
    }
    ////////////////////////////////////////////
    // PSF #3: linear superposition of two radial symmetric Gaussians
    else if( fFitter->fFitterName == "LinearSuperposition2DGauss_LL" )
    {
    
        // define graph
        TGraphErrors* g = new TGraphErrors( 1 );
        setGraphPlottingStyle( g, iColor, 2., 7, 2. );
        
        g->SetPoint( 0, fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1] );
        g->SetPointError( 0, fFitter->fFitResult_ParameterError[0], fFitter->fFitResult_ParameterError[1] );
        g->Draw( "p" );
        
        
        if( fFitter->fFitResult_Parameter[2] > 0. )
        {
            TEllipse* e = new TEllipse( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[2] );
            e->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e->Draw();
            
            TEllipse* e2 = new TEllipse( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[3] );
            e2->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e2->Draw();
        }
        
        return g;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Radial Symmetric Source
    else if( fFitter->fFitterName == "RadialSymmetricSource_Chi2" || fFitter->fFitterName == "RadialSymmetricSource_LL" )
    {
        // define graph
        TGraphErrors* g = new TGraphErrors( 1 );
        setGraphPlottingStyle( g, iColor, 2., 7, 2. );
        
        g->SetPoint( 0, fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1] );
        g->SetPointError( 0, fFitter->fFitResult_ParameterError[0], fFitter->fFitResult_ParameterError[1] );
        g->Draw( "p" );
        
        if( fFitter->fFitResult_Parameter[2] > 0. )
        {
            TEllipse* e = new TEllipse( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[2] );
            e->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e->Draw();
        }
        
        return g;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Radial Asymmetric Source
    else if( fFitter->fFitterName == "RadialAsymmetricSource_Chi2" )
    {
        // define graph
        TGraphErrors* g = new TGraphErrors( 1 );
        setGraphPlottingStyle( g, iColor, 2., 7, 2. );
        
        g->SetPoint( 0, fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1] );
        g->SetPointError( 0, fFitter->fFitResult_ParameterError[0], fFitter->fFitResult_ParameterError[1] );
        g->Draw( "p" );
        
        if( fFitter->fFitResult_Parameter[2] > 0. )
        {
            TEllipse* e = new TEllipse( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[4] );
            e->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e->Draw();
        }
        
        return g;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Radial Asymmetric Source
    else if( fFitter->fFitterName == "RadialAsymmetricSource_LL" )
    {
        // define graph
        TGraphErrors* g = new TGraphErrors( 1 );
        setGraphPlottingStyle( g, iColor, 2., 7, 2. );
        
        g->SetPoint( 0, fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[3] );
        g->SetPointError( 0, fFitter->fFitResult_ParameterError[1], fFitter->fFitResult_ParameterError[3] );
        g->Draw( "p" );
        
        if( fFitter->fFitResult_Parameter[2] > 0. )
        {
            // double angle = 2*fFitter->fFitResult_Parameter[0]*fFitter->fFitResult_Parameter[2]*fFitter->fFitResult_Parameter[4];
            //angle /= (fFitter->fFitResult_Parameter[2]*fFitter->fFitResult_Parameter[2] - fFitter->fFitResult_Parameter[4]*fFitter->fFitResult_Parameter[4]);
            //angle = 1./2*atan(angle);
            
            double angle = fFitter->fFitResult_Parameter[0];
            
            TEllipse* e = new TEllipse( fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[4], 0, 360,  angle * 180. / acos( -1. ) );
            e->SetFillStyle( 0 );
            e->SetLineWidth( 2 );
            e->SetLineColor( kPink );
            e->Draw();
        }
        
        return g;
    }
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // 2D normal distribution
    /*  else if( fFitter->fFitterName == "2DAsymGauss_LL" )
    {
     // define graph
       TGraphErrors *g = new TGraphErrors( 1 );
       setGraphPlottingStyle( g, iColor, 2., 7, 2. );
    
       g->SetPoint( 0, fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[3] );
       g->SetPointError( 0, fFitter->fFitResult_ParameterError[1], fFitter->fFitResult_ParameterError[3] );
       g->Draw( "p" );
    
       if( fFitter->fFitResult_Parameter[2] > 0. )
       {
     TEllipse *e = new TEllipse( fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[2], fFitter->fFitResult_Parameter[4] );
     e->SetFillStyle( 0 );
     e->Draw();
       }
    
       return g;
    }
    */
    
    return 0;
}

bool VSourceGeometryFitter::setFitter( string iDescription )
{
    for( unsigned int i = 0; i < fDefaultFitterData.size(); i++ )
    {
        if( fDefaultFitterData[i]->fFitterName == iDescription )
        {
            fFitter = fDefaultFitterData[i];
            return true;
        }
    }
    fFitter = 0;
    
    return false;
}

void VSourceGeometryFitter::fitSource( string iHisName, double xStart, double yStart, double xyRange )
{

    fXStart = xStart;
    fYStart = yStart;
    fPSF = getPSF();
    
    if( !fFitter )
    {
        cout << "VSourceGeometryFitter::fitSource: undefined fitter" << endl;
        return;
    }
    
    // get sky map histogram
    fHisSkyMap = ( TH2D* )getHistogram( iHisName, fRunNumber, "skyHistograms" );
    
    if( !fHisSkyMap )
    {
        cout << "VSourceGeometryFitter::fitSource: histogram not found: " << iHisName << endl;
        return;
    }
    
    //////////////////////////////////////
    // define minuit
    //////////////////////////////////////
    ROOT::Minuit2::Minuit2Minimizer* fSourceGeometryFitter = new ROOT::Minuit2::Minuit2Minimizer();
    
    //////////////////////////////////////
    // set fit function
    //////////////////////////////////////
    
    // Source #1 radial symmetric source, Chi2
    VFun_SourceDescription_RadialSymmetricSource_Chi2 fcn_RadialSymmetricSource_Chi2( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange, fPSF );
    if( fFitter->fFitterName == "RadialSymmetricSource_Chi2" )
    {
        fSourceGeometryFitter->SetFunction( fcn_RadialSymmetricSource_Chi2 );
        
        // update parameters
        fFitter->fParameterInitValue[0]  = xStart;
        fFitter->fParameterLowerLimit[0] = xStart - xyRange;
        fFitter->fParameterUpperLimit[0] = xStart + xyRange;
        fFitter->fParameterInitValue[1]  = yStart;
        fFitter->fParameterLowerLimit[1] = yStart - xyRange;
        fFitter->fParameterUpperLimit[1] = yStart + xyRange;
    }
    
    // Source #2 radial asymmetric source, Chi2
    /*  VFun_SourceDescription_RadialAsymmetricSource_Chi2 fcn_RadialAsymmetricSource_Chi2( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange );
    if( fFitter->fFitterName == "RadialAsymmetricSource_Chi2" )
    {
      fSourceGeometryFitter_MINUIT->SetMinuitFCN( &fcn_RadialAsymmetricSource_Chi2 );
      // update parameters
      fFitter->fParameterInitValue[0]  = xStart;
      fFitter->fParameterLowerLimit[0] = xStart - xyRange;
      fFitter->fParameterUpperLimit[0] = xStart + xyRange;
      fFitter->fParameterInitValue[1]  = yStart;
      fFitter->fParameterLowerLimit[1] = yStart - xyRange;
      fFitter->fParameterUpperLimit[1] = yStart + xyRange;
      }*/
    
    
    // Source #3 radial symmetric source, LL
    VFun_SourceDescription_RadialSymmetricSource_LL fcn_RadialSymmetricSource_LL( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange, fPSF );
    if( fFitter->fFitterName == "RadialSymmetricSource_LL" )
    {
        fSourceGeometryFitter->SetFunction( fcn_RadialSymmetricSource_LL );
        
        // update parameters
        fFitter->fParameterInitValue[0]  = xStart;
        fFitter->fParameterLowerLimit[0] = xStart - xyRange;
        fFitter->fParameterUpperLimit[0] = xStart + xyRange;
        fFitter->fParameterInitValue[1]  = yStart;
        fFitter->fParameterLowerLimit[1] = yStart - xyRange;
        fFitter->fParameterUpperLimit[1] = yStart + xyRange;
    }
    
    // Source #4 radial asymmetric source, LL
    VFun_SourceDescription_RadialAsymmetricSource_LL fcn_RadialAsymmetricSource_LL( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange, fPSF );
    if( fFitter->fFitterName == "RadialAsymmetricSource_LL" )
    {
        fSourceGeometryFitter->SetFunction( fcn_RadialAsymmetricSource_LL );
        
        // update parameters
        fFitter->fParameterInitValue[1]  = xStart;
        fFitter->fParameterLowerLimit[1] = xStart - xyRange;
        fFitter->fParameterUpperLimit[1] = xStart + xyRange;
        fFitter->fParameterInitValue[3]  = yStart;
        fFitter->fParameterLowerLimit[3] = yStart - xyRange;
        fFitter->fParameterUpperLimit[3] = yStart + xyRange;
    }
    
    
    //// 0thers
    // Source description #3
    /*    VFun_SourceDescription_2DNormal_LL fcn_2DNormal_LL( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange );
        if( fFitter->fFitterName == "2DAsymGauss_LL" )
        {
           fSourceGeometryFitter_MINUIT->SetMinuitFCN( &fcn_2DNormal_LL );
    // update parameters
           fFitter->fParameterInitValue[1]  = xStart;
           fFitter->fParameterLowerLimit[1] = xStart - xyRange;
           fFitter->fParameterUpperLimit[1] = xStart + xyRange;
           fFitter->fParameterInitValue[3]  = yStart;
           fFitter->fParameterLowerLimit[3] = yStart - xyRange;
           fFitter->fParameterUpperLimit[3] = yStart + xyRange;
        }
    
    */
    
    // PSF description #1
    VFun_PSFDescription_2DGauss_Chi2 fcn_2DGauss_Chi2( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange );
    if( fFitter->fFitterName == "2DGauss_Chi2" )
    {
        fSourceGeometryFitter->SetFunction( fcn_2DGauss_Chi2 );
        
        // update parameters
        fFitter->fParameterInitValue[3]  = xStart;
        fFitter->fParameterLowerLimit[3] = xStart - xyRange;
        fFitter->fParameterUpperLimit[3] = xStart + xyRange;
        fFitter->fParameterInitValue[4]  = yStart;
        fFitter->fParameterLowerLimit[4] = yStart - xyRange;
        fFitter->fParameterUpperLimit[4] = yStart + xyRange;
    }
    
    // PSF description #2
    VFun_PSFDescription_2DGauss_LL fcn_2DGauss_LL( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange );
    if( fFitter->fFitterName == "2DGauss_LL" )
    {
        fSourceGeometryFitter->SetFunction( fcn_2DGauss_LL );
        
        // update parameters
        fFitter->fParameterInitValue[0]  = xStart;
        fFitter->fParameterLowerLimit[0] = xStart - xyRange;
        fFitter->fParameterUpperLimit[0] = xStart + xyRange;
        fFitter->fParameterInitValue[1]  = yStart;
        fFitter->fParameterLowerLimit[1] = yStart - xyRange;
        fFitter->fParameterUpperLimit[1] = yStart + xyRange;
    }
    
    // PSF description #3
    VFun_PSFDescription_LinearSuperposition2DGauss_LL fcn_LinearSuperposition2DGauss_LL( fHisSkyMap, xStart - xyRange, xStart + xyRange, yStart - xyRange, yStart + xyRange );
    if( fFitter->fFitterName == "LinearSuperposition2DGauss_LL" )
    {
        fSourceGeometryFitter->SetFunction( fcn_LinearSuperposition2DGauss_LL );
        
        // update parameters
        fFitter->fParameterInitValue[0]  = xStart;
        fFitter->fParameterLowerLimit[0] = xStart - xyRange;
        fFitter->fParameterUpperLimit[0] = xStart + xyRange;
        fFitter->fParameterInitValue[1]  = yStart;
        fFitter->fParameterLowerLimit[1] = yStart - xyRange;
        fFitter->fParameterUpperLimit[1] = yStart + xyRange;
    }
    
    // set parameters
    for( unsigned int i = 0; i < fFitter->fParameterName.size(); i++ )
    {
        fSourceGeometryFitter->SetLimitedVariable( i, fFitter->fParameterName[i].c_str(), fFitter->fParameterInitValue[i], fFitter->fParameterStep[i], fFitter->fParameterLowerLimit[i], fFitter->fParameterUpperLimit[i] );
    }
    
    
    // start minimizing
    // (default is kMigrad)
    if( fSourceGeometryFitter->Minimize() == true )	// convergence
    {
        fFitter->fFitResult_Status = 0;
    }
    else
    {
        fFitter->fFitResult_Status = fSourceGeometryFitter->Status();
    }
    
    cout << "Fit status " << fFitter->fFitResult_Status << endl;
    
    // retrieve parameters
    fFitter->fFitResult_Parameter.clear();
    fFitter->fFitResult_ParameterError.clear();
    for( unsigned int i = 0; i < fFitter->fParameterName.size(); i++ )
    {
        fFitter->fFitResult_Parameter.push_back( fSourceGeometryFitter->X()[i] );
        fFitter->fFitResult_ParameterError.push_back( fSourceGeometryFitter->Errors()[i] );
    }
    
    
    if( fFitter->fFitterName == "RadialAsymmetricSource_LL" )
    {
        //double rho = fFitter->fFitResult_Parameter[0];
        double angle = fFitter->fFitResult_Parameter[0];
        double sX = sqrt( fFitter->fFitResult_Parameter[2] * fFitter->fFitResult_Parameter[2] + fPSF * fPSF );
        double sY = sqrt( fFitter->fFitResult_Parameter[4] * fFitter->fFitResult_Parameter[4] + fPSF * fPSF );
        
        double angle_err = fFitter->fFitResult_ParameterError[0];
        double rho =  1. / 2. * tan( 2 * angle ) * ( sX * sX - sY * sY ) / sX / sY ;
        
        double p1 = sX * sX * sY * sY * ( 1 - rho * rho );
        p1 = p1 / ( sY * sY * pow( cos( angle ), 2 ) - 2 * rho * sX * sY * sin( angle ) * cos( angle ) + sX * sX * pow( sin( angle ), 2 ) );
        p1 = sqrt( p1 );
        
        double p2 = sX * sX * sY * sY * ( 1 - rho * rho );
        p2 = p2 / ( sY * sY * pow( sin( angle ), 2 ) + 2 * rho * sX * sY * sin( angle ) * cos( angle ) + sX * sX * pow( cos( angle ), 2 ) );
        p2 = sqrt( p2 );
        
        // approximate error, rho approximately 0
        //double errTanAngle =  pow(2*sX*sY/( sX * sX + sY * sY ) * rho_err, 2)  + pow(2*rho*sY/( sX * sX + sY * sY) + 2 * rho * sX * sY /(sX*sX+sY*sY)/(sX*sX+sY*sY)*2*sX ,2)*sX_err*sX_err ;
        //errTanAngle = sqrt(errTanAngle);
        //double errAngle = 1./2/(1+angle*angle)*errTanAngle;
        
        
        cout << "The semi-axis of the ellipse are p1 = " << p1 << " deg  and p2 = " << p2 << " deg " << endl;
        cout << "The correlation between X and Y is defined by the parameter rho = " << rho << endl;
        cout << "The angle a is the angle between the X axis and the semi-axis of the ellipse p1. Expressed in degrees the result of the fit gives an  angle a = " << angle * 180 / acos( -1 ) << " +/- " << angle_err * 180 / acos( -1 ) << " deg " << endl;
        
        
    }
    
    ///////////////////////////////////////////////////////////////
    // convert centroid from camera coordinates to sky coordinates
    ///////////////////////////////////////////////////////////////
    if( fFitter->fFitterName == "2DGauss_Chi2" )
    {
        VPlotAnasumHistograms myAna( fAnasumDataFile, -1 );
        myAna.convert_derotated_RADECJ2000( fFitter->fFitResult_Parameter[3], fFitter->fFitResult_Parameter[4], fFitter->fFitResult_ParameterError[3], fFitter->fFitResult_ParameterError[4] );
    }
    if( fFitter->fFitterName == "2DGauss_LL" || fFitter->fFitterName == "LinearSuperposition2DGauss_LL" )
    {
        VPlotAnasumHistograms myAna( fAnasumDataFile, -1 );
        myAna.convert_derotated_RADECJ2000( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_ParameterError[0], fFitter->fFitResult_ParameterError[1] );
    }
    if( fFitter->fFitterName == "RadialSymmetricSource_Chi2" || fFitter->fFitterName == "RadialSymmetricSource_LL" )
    {
        VPlotAnasumHistograms myAna( fAnasumDataFile, -1 );
        myAna.convert_derotated_RADECJ2000( fFitter->fFitResult_Parameter[0], fFitter->fFitResult_Parameter[1], fFitter->fFitResult_ParameterError[0], fFitter->fFitResult_ParameterError[1] );
    }
    if( fFitter->fFitterName == "RadialAsymmetricSource_LL" )
    {
        VPlotAnasumHistograms myAna( fAnasumDataFile, -1 );
        myAna.convert_derotated_RADECJ2000( fFitter->fFitResult_Parameter[1], fFitter->fFitResult_Parameter[3], fFitter->fFitResult_ParameterError[1], fFitter->fFitResult_ParameterError[3] );
    }
    
    
    // combined error calculation for 5 parameters (see Minuit manual Table 7.1)
    //    fMinuit.Command( "SET ERR 6.06" );
    
    delete fSourceGeometryFitter;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
     Functions for source position and extension fitting
*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



VFun_PSFDescription_2DGauss_Chi2::VFun_PSFDescription_2DGauss_Chi2( TH2D* iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax )
{
    hSkyMap = iSkymap;
    xmin = i_xmin;
    xmax = i_xmax;
    ymin = i_ymin;
    ymax = i_ymax;
}


VFun_PSFDescription_2DGauss_LL::VFun_PSFDescription_2DGauss_LL( TH2D* iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax )
{
    hSkyMap = iSkymap;
    xmin = i_xmin;
    xmax = i_xmax;
    ymin = i_ymin;
    ymax = i_ymax;
}


VFun_PSFDescription_LinearSuperposition2DGauss_LL::VFun_PSFDescription_LinearSuperposition2DGauss_LL( TH2D* iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax )
{
    hSkyMap = iSkymap;
    xmin = i_xmin;
    xmax = i_xmax;
    ymin = i_ymin;
    ymax = i_ymax;
}



/*
    radial symetrical source (Gaussian shape) (Chi2)

    fit parameters are:   source position, source extension, normalisation

    input:                TH2 histogram with on counts
    return value:         chi2, fit parameters

    units in [deg]
*/

VFun_SourceDescription_RadialSymmetricSource_Chi2::VFun_SourceDescription_RadialSymmetricSource_Chi2( TH2D* iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax, double i_psf )
{
    hSkyMap = iSkymap;
    xmin = i_xmin;
    xmax = i_xmax;
    ymin = i_ymin;
    ymax = i_ymax;
    sigmaPSF = i_psf;
}



/*
    radial asymetrical source (Gaussian shape) (Chi2)

    fit parameters are:   source position, source extension, normalisation

    input:                TH2 histogram with on counts
    return value:         chi2, fit parameters

    units in [deg]
*/

/*VFun_SourceDescription_RadialAsymmetricSource_Chi2::VFun_SourceDescription_RadialAsymmetricSource_Chi2( TH2D *iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax )
{
   hSkyMap = iSkymap;
   xmin = i_xmin;
   xmax = i_xmax;
   ymin = i_ymin;
   ymax = i_ymax;
   }*/



/*
    2D normal distribution (LL)

    fit parameters are:   source position, source extension, normalisation

    input:                TH2 histogram with on counts
    return value:         chi2, fit parameters

    units in [deg]
*/


VFun_SourceDescription_RadialSymmetricSource_LL::VFun_SourceDescription_RadialSymmetricSource_LL( TH2D* iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax, double i_psf )
{
    hSkyMap = iSkymap;
    xmin = i_xmin;
    xmax = i_xmax;
    ymin = i_ymin;
    ymax = i_ymax;
    sigmaPSF = i_psf;
}


VFun_SourceDescription_RadialAsymmetricSource_LL::VFun_SourceDescription_RadialAsymmetricSource_LL( TH2D* iSkymap, double i_xmin, double i_xmax, double i_ymin, double i_ymax, double i_psf )
{
    hSkyMap = iSkymap;
    xmin = i_xmin;
    xmax = i_xmax;
    ymin = i_ymin;
    ymax = i_ymax;
    sigmaPSF = i_psf;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* \class VSourceGeometryFitterData
   \brief containts fit defaults and results

*/

VSourceGeometryFitterData::VSourceGeometryFitterData()
{
    fFitResult_Status = -99;
}


