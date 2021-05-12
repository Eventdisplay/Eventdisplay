/*! \class VImageParameterFitter

    image fitter to calculate second moment paramters

    note that for efficiency reasons, only after calling calcLL 
    we check for zero pointers

*/

 #include <VImageParameterFitter.h>

VImageParameterFitter::VImageParameterFitter( VEvndispData* iData,
                                              bool iDebug,
                                              double iZeroTolerence )
{
    fDebug = iDebug;
    fLLDebug = false;
    if( fDebug )
    {
        cout << "VImageParameterFitter::VImageParameterFitter()" << endl;
    }
    iZeroTolerence = ZeroTolerence;

    // default data sources and result classes
    fData = iData;
    fParLL = 0;
    fParGeo = 0;
    fNormal2D = 0;

    // fit function: 
    // rotated normal distribution provides
    // with errors which are less correlated
    // bRotatedNormalDistributionFit = false;
    bRotatedNormalDistributionFit = true;
  
    resetFitParameters();
}

VImageParameterFitter::~VImageParameterFitter()
{
    if( fLLFitter )
    {
        delete fLLFitter;
    }
}


/*!
   \param iVmode Minuit print mode, 1 = quit, 2 - verbose

*/
void VImageParameterFitter::initMinuit( int iVmode )
{
    if( iVmode == 2 )
    {
        fLLDebug = true;
    }
    unsigned int iFitParameter = 6;
    if( fData->getRunParameter()->fMinimizeTimeGradient ) iFitParameter = 9;
    fLLFitter = new TMinuit( iFitParameter );
    // no minuit printouts
    if( iVmode == 1 )
    {
        fLLFitter->Command( "SET PRINT -1" );
    }
    // minuit printouts
    else if( iVmode == 2 )
    {
        fLLFitter->Command( "SET PRINT 1" );
    }
    fLLFitter->Command( "SET NOWA" );             // no warnings
    fLLFitter->Command( "SET ERR 0.5" );          // loglikelihood -> UP = 0.5  (solves 2*ll = chi2 != 1)
    fLLFitter->SetObjectFit( this );
    if( bRotatedNormalDistributionFit )
    {
        fLLFitter->SetFCN( get_LL_imageParameter_2DGaussRotated );
    }
    else
    {
        fLLFitter->SetFCN( get_LL_imageParameter_2DGauss );
    }
    // 2D function for integrated probability
    // - range is changed before calling
    // - number of parameters does not include charge normalisation
    // (optional and in testing only)
    fNormal2D = 0;
    //fNormal2D = new TF2( "normal2D", normal2DRotated, -20., 20., -20., 20., iFitParameter - 2 );
}


/*!
    The loglikelihood fit is fitting a 2D Gaussian to the image/border pixels

    FCN is get_LL_imageParameter_2DGauss()

    Only signals in image/border pixels are taken into account, all other are set to zero.

    Signals are estimated for dead pixels.
    If estimated signal is above image/border threshold, VEvndispData::getSums()[pixelID] is set to
    this value and the fLLEst[pixelID] is set true.
    A new value for the size is calculated.
    This is the old size from the geometrical calculation  plus the estimated signals from the
    dead pixels.

   \attention

    Check always in the analysis the status of the fit (FitStat).
    Values lower than 3 (approximate or no error matrix) should be excluded from the analysis.
    In this case there is a large probablitity that the parameters are wrong.

*/
vector<bool> VImageParameterFitter::calcLL( VImageParameter *iParGeo,
                                            VImageParameter *iParLL,
                                            bool iUseSums2, 
                                            bool i_reInitializeLL, 
                                            bool iEqualSummationWindows )
{
    fParGeo = iParGeo;
    fParLL = iParLL;
    if( !fData )
    {
        vector< bool > a;
        cout << "VImageParameterCalculation::calcLL error: data vector is zero" << endl;
        return a;
    }
    if( !fParGeo || !fParLL )
    {
        vector< bool > a;
        cout << "VImageParameterCalculation::calcLL error: no image parameters given " << endl;
        return a;
    }
    fLLDebug = false;
    if( fLLDebug )
    {
        cout << endl;
        cout << "=================================================================" << endl;
        cout << "Telescope " << fData->getTelID() + 1 << endl;
        cout << "=================================================================" << endl;
        cout << endl;
    }
    
    // reset fit variables
    resetFitParameters();

    // get pixel values, all nonimage/nonborder pixels have fData->getSums()=0.
    // if there are dead channels, size of these vectors is < fData->getSums().size() !!
    signal = fill_pixel_sums( iUseSums2 );
    
    // set start values for LL fit
    defineFitParameters();

    // now do the minimization
    fLLFitter->Command( "MIGRAD" );

    // fit statistics
    getFitStatistics();

    // call HESSE in case error calculation failed
    if( fParLL->Fitstat < 3 )
    {
        fLLFitter->Command( "HESSE" );
        getFitStatistics();
    }
    
    // fit results
    getFitResults();

    // calculate image parameters from fit results
    calculateImageParameters( iUseSums2, iEqualSummationWindows );

    return fLLEst;
}

/*
    calculate image brightness from 2D fit
*/
double VImageParameterFitter::calculatePixelBrightness( unsigned int iChannel, 
                                                        double rho, double phi,
                                                        double meanX, double sigmaX, 
                                                        double meanY, double sigmaY, 
                                                        double signal )
{
    double f = 0;
    // get channel coordinates
    double x = getDetectorGeometry()->getX()[iChannel];
    double y = getDetectorGeometry()->getY()[iChannel];
    // calculate 2D-gauss
    if( bRotatedNormalDistributionFit )
    {
        double x_p =     x*cos(phi) + y*sin(phi);
        double y_p = -1.*x*sin(phi) + y*cos(phi);
        double cx_p =     meanX*cos(phi) + meanY*sin(phi);
        double cy_p = -1.*meanX*sin(phi) + meanY*cos(phi);

        f  = (x_p - cx_p) * (x_p - cx_p) / sigmaX / sigmaX;
        f += (y_p - cy_p) * (y_p - cy_p) / sigmaY / sigmaY;
        f *= -1. / 2.;
        f  = exp( f );
        f *= 1. / 2. / TMath::Pi() / sigmaX / sigmaY;
    } 
    else
    {
        f  = ( x - meanX ) * ( x - meanX ) / sigmaX / sigmaX;
        f += ( y - meanY ) * ( y - meanY ) / sigmaY / sigmaY;
        f += -2. * rho * ( x - meanX ) / sigmaX * ( y - meanY ) / sigmaY;
        f *= -1. / 2. / ( 1. - rho * rho );
        f  = exp( f );
        f *= 1. / 2. / TMath::Pi() / sigmaX / sigmaY / sqrt( 1. - rho * rho );
    }
    f *= signal;
    
    return f;
}


/*!
   reduce angles to values in [0.,iMax]
*/
double VImageParameterFitter::redang( double iangle, double iMax )
{
    if( iMax == 0. )
    {
        return 0.;
    }
    if( iangle >= 0 )
    {
        iangle = iangle - int( iangle / iMax ) * iMax;
    }
    else
    {
        iangle = iMax + iangle + int( iangle / iMax ) * iMax;
    }
    
    return iangle;
}


/*! \fn get_LL_imageParameter_2DGauss( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
    \brief function for the loglikelihood calculation of image parameters in minuit (2dim-Gaus)

   gives same parameters as geometrical analysis, but with errors

   log likelihood: see Blobel p.195

   The image is described by a 2d-Gaussian. The signal in each tube is therefore:
   \f[
           S(x,y)=\frac{C}{2\pi\sigma_{x}\sigma_{y}\sqrt{1-\rho^{2}}}
              \exp\left\{-\frac{1}{2(1-\rho^{2})}
\left[\left(\frac{x-c_{x}}{\sigma_{x}}\right)^{2}
-2\rho\left(\frac{x-c_{x}}{\sigma_{x}}\right)\left(\frac{y-c_{y}}{\sigma_{y}}\right)
+\left(\frac{y-c_{y}}{\sigma_{y}}\right)^{2}\right]\right\}

\f]

Assume Poissonian (ignoring NSB) for measuring amplitude:

\f[
P(n_{i}; S ) = \frac{S^{n}}{n!}\exp^{-S}
\f]

This gives the extented log likelihood estimator:

\f[
LL = - \sum (n_{i}\ln S_{i} - S_{i} - n_{i}\ln n_{i} + n_{i})
\f]

Errors are calculated from this estimator by

\f[
2 LL \approx \chi^{2}
\f]

*/
void get_LL_imageParameter_2DGauss( Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag )
{
    double LL = 0.;
    double LL_t = 0.;
    double sum = 0.;
    
    npar = 6;
    gin = 0;
    iflag = 0;
    
    double x = 0.;
    double y = 0.;
    double n = 0.;
    double t = 0.;
    double tx = 0.;
    double phi = 0.;
    double t_sig = 4.;
    double t_sig_term = 0.;
    
    double rho_1 = -1. / 2. / ( 1. - par[0] * par[0] );
    double rho_s =  1. / 2. / M_PI / par[2] / par[4] / sqrt( 1. - par[0] * par[0] ) * par[5];
    
    VImageParameterFitter* iImageCalculation = ( VImageParameterFitter* )fLLFitter->GetObjectFit();

    if( iImageCalculation->minimize_time_gradient_for_this_event() )
    {
        t_sig = par[8];
        phi = atan2( 2.*par[0] * par[2] * par[4], 
                        par[2] * par[2] - par[4] * par[4] ) / 2.;
        t_sig_term = log( 1. / sqrt( 2. * M_PI * t_sig ) );
    }
    
    if( par[0] * par[0] < 1. && par[2] > 0. && par[4] > 0. )
    {
        unsigned int nSums = iImageCalculation->getLLSums().size();
        for( unsigned int i = 0; i < nSums; i++ )
        {
            n = iImageCalculation->getLLSums()[i];
            if( n > -999. )
            {
                x = iImageCalculation->getLLX()[i];
                y = iImageCalculation->getLLY()[i];
                sum  = ( x - par[1] ) * ( x - par[1] ) / par[2] / par[2];
                sum += ( y - par[3] ) * ( y - par[3] ) / par[4] / par[4];
                sum += -2. * par[0] * ( x - par[1] ) / par[2] * ( y - par[3] ) / par[4];
                sum  = rho_s * exp( sum * rho_1 );
                
                // assume Poisson fluctuations (neglecting background noise)
                if( n > 0. && sum > 0. )
                {
                    LL += n * log( sum ) - sum - n * log( n ) + n;
                }
                else
                {
                    LL += -1. * sum;
                }
                // Gauss with background noise
                //             iVar = sum+pedvars*pedvars;
                //             iVar = sum;
                //             LL += -0.5*log(2. * TMath::Pi()) - log(sqrt(iVar)) - (n-sum)*(n-sum)/2./iVar;
                // time gradient analysis
                // (line fit to time gradient)
                // - affects image orientation and centroid position
                if( iImageCalculation->minimize_time_gradient_for_this_event() )
                {
                    t = iImageCalculation->getLLT()[i];
                    if( t > 0. )
                    {
                        // calculate position along long axis
                        tx = (x-par[1]) * cos(phi) + (y-par[3]) * sin(phi);
                        LL_t += 
                            - 1./2./(t_sig*t_sig) 
                            * (t - par[6] - par[7] * tx )
                            * (t - par[6] - par[7] * tx );
                    }
                }

            }
        }
    }
    LL_t += t_sig_term;
    f = -1. * ( LL + LL_t );
}

/*
    decide if time gradient should be used

    not for
    - very flat images
    - small images (should be more pixel values than fit parameters

*/
bool VImageParameterFitter::minimize_time_gradient_for_this_event()
{
    if( minimize_time_gradient()
        && TMath::Abs( fParGeo->tgrad_x ) > fData->getRunParameter()->fMinimizeTimeGradient_minGradforFit
        && fParGeo->ntubes >= fData->getRunParameter()->fMinimizeTimeGradient_minNtubes )
    {
        return true;
    }
    return false;
}

/*
   LL optimisation: starting value for rho
*/
double VImageParameterFitter::getLL_startingvalue_rho()
{
    if( sigmaX > ZeroTolerence && sigmaY > ZeroTolerence )
    {
        rho = tan( 2. * fParGeo->phi );
        rho *= ( sigmaX * sigmaX - sigmaY * sigmaY );
        rho /= sigmaX / sigmaY / 2.;
    }
    else
    {
        rho = 0.;
    }
    return rho;
}

/*
  LL optimisation: parameter limits for rho
*/
double VImageParameterFitter::getLL_paramameterlimits_rho( double rho, double upper_sign )
{
    double rho_edge = upper_sign * 0.998;
    // TMP (this might lead to spurious results otherwise
    return rho_edge;
    // Case 1: image is well inside the camera and sufficiently large
    //         (don't expect large difference in phi)
    if( minimize_time_gradient_for_this_event() )
    {
        rho_edge = rho + upper_sign * 0.1;
    }
    if( upper_sign > 0. && rho_edge > 0.998 )  return 0.998;
    if( upper_sign < 0. && rho_edge < -0.998 ) return -0.998;

    return rho_edge;
}

/*
   LL optimisation: starting value for cen(_x or _y)
*/
double VImageParameterFitter::getLL_paramameterlimits_cen( double dist_limit,
                                                           double centroid, 
                                                           double fLL_StartingValue_sigma, 
                                                           double i_sign )
{
    // temporary settings for rotated normal distribution
/*    if( bRotatedNormalDistributionFit )
    {
        return i_sign * 5.;
    } */
    // Case 1: image is well inside the camera
    //         assumption is that geo centroids are good values
    if( fLL_StartingValue_sigma > 0. && fParGeo->loss < fData->getRunParameter()->fMinimizeTimeGradient_minLoss )
    {
        double camera_pixel_size = 0.1;
        // assume that zero pixel is a typical pixel
        if( fData->getDetectorGeometry()->getTubeRadius().size() > 0 )
        {
            camera_pixel_size = 5.*fData->getDetectorGeometry()->getTubeRadius()[0];
            if( camera_pixel_size < 0.1 ) camera_pixel_size = 5.*0.1;
        }
        return centroid + i_sign * camera_pixel_size;
    }
    // Case 2: image is partly at the edge of the FOV: determine box around image
    //         to search for centroid
    if( fLL_StartingValue_sigma > 0. )
    {
        dist_limit = centroid + i_sign* 2.*fLL_StartingValue_sigma;
        // make sure that this is inside the FOV (+10%)
        if( fData->getDetectorGeometry() &&  fData->getTelID() < fData->getDetectorGeometry()->getFieldofView().size() )
        {
            if( TMath::Abs( dist_limit ) > 1.1 * 0.5 * fData->getDetectorGeometry()->getFieldofView()[fData->getTelID()] )
            {
                if( dist_limit > 0 )
                {
                    dist_limit =  1.1 * 0.5 * fData->getDetectorGeometry()->getFieldofView()[fData->getTelID()];
                }
                else
                {
                    dist_limit = -1.1 * 0.5 * fData->getDetectorGeometry()->getFieldofView()[fData->getTelID()];
                }
            }
        }
    }
    else
    {
        // image centroid should not be outside of the FOV (by more than 10%)
        if( fData->getDetectorGeometry() &&  fData->getTelID() < fData->getDetectorGeometry()->getFieldofView().size() )
        {
            dist_limit = i_sign *1.1 * 0.5 * fData->getDetectorGeometry()->getFieldofView()[fData->getTelID()];
        }
        // should never land here
        else
        {
            dist_limit = i_sign * 5.;
        }
    }
    return dist_limit;

}
/*! \fn get_LL_imageParameter_2DGaussRotated( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
    \brief function for the loglikelihood calculation of image parameters in minuit (2dim-Gaus)

   optimised for disp reconstruction

   log likelihood: see Blobel p.195

   The image is described by a 2d-Gaussian. The signal in each tube is therefore:
   \f[
           S(x,y)=\frac{C}{2\pi\sigma_{x}\sigma_{y}}}}
              \exp\left\{-\frac{1}{2}
            \left[\left(\frac{x-c_{x}}{\sigma_{x}}\right)^{2}
            +\left(\frac{y-c_{y}}{\sigma_{y}}\right)^{2}\right]\right\}

\f]

Assume Poissonian (ignoring NSB) for measuring amplitude:

\f[
P(n_{i}; S ) = \frac{S^{n}}{n!}\exp^{-S}
\f]

This gives the extented log likelihood estimator:

\f[
LL = - \sum (n_{i}\ln S_{i} - S_{i} - n_{i}\ln n_{i} + n_{i})
\f]

Errors are calculated from this estimator by

\f[
2 LL \approx \chi^{2}
\f]

*/
void get_LL_imageParameter_2DGaussRotated( Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag )
{
    double LL = 0.;
    double LL_t = 0.;
    double sum = 0.;
    
    npar = 6;
    gin = 0;
    iflag = 0;
    
    double x = 0.;
    double y = 0.;
    double r = 0.;
    double x_p = 0.;
    double y_p = 0.;
    double cx_p =     par[1]*cos(par[0]) + par[3]*sin(par[0]);
    double cy_p = -1.*par[1]*sin(par[0]) + par[3]*cos(par[0]);
    double n = 0.;
    double t = 0.;
    double tx = 0.;
    double t_sig = 4.;
    double t_sig_term = 0.;

    // normalisation factor
    double rho_1 = -1. / 2.; 

    VImageParameterFitter* iImageCalculation = (VImageParameterFitter*)fLLFitter->GetObjectFit();

    if( iImageCalculation->minimize_time_gradient_for_this_event() )
    {
        t_sig = par[8];
        t_sig_term = log( 1. / sqrt( 2. * M_PI * t_sig ) );
    }

    // for probability densitivity integration
    // (optional)
    TF2 *fNormal2D = iImageCalculation->getNormal2D();
    if( fNormal2D )
    {
        fNormal2D->SetParameter( 0, cx_p );
        fNormal2D->SetParameter( 1, par[2] );
        fNormal2D->SetParameter( 2, cy_p );
        fNormal2D->SetParameter( 3, par[4] );
        fNormal2D->SetRange( par[1] - 2., par[3] - 2., par[1] + 2., par[3] + 2. );
    }


    if( par[2] > 0. && par[4] > 0. )
    {
        unsigned int nSums = iImageCalculation->getLLSums().size();
        for( unsigned int i = 0; i < nSums; i++ )
        {
            n = iImageCalculation->getLLSums()[i];
            if( n > -999. )
            {
                // work in rotated space
                x = iImageCalculation->getLLX()[i];
                y = iImageCalculation->getLLY()[i];

                // rotated coordinate system
                x_p = x*cos(par[0]) + y*sin(par[0]);
                y_p = -1.*x*sin(par[0]) + y*cos(par[0]);
                r = iImageCalculation->getLLR()[i] / 2.;

                // use integral of probability distribution
                // (integrate over pixels)
                // - test show that this method is inferior 
                //   to all others
                if( fNormal2D )
                {
                    sum = fNormal2D->Integral( x_p - r, x_p + r, 
                                               y_p - r, y_p + r );
                }
                // use probability densitiy at centre of 
                // pixel position
                else
                {
                    sum  = (x_p-cx_p)*(x_p-cx_p) / par[2] / par[2]
                         + (y_p-cy_p)*(y_p-cy_p) / par[4] / par[4];

                    sum  = 1. / 2. / M_PI / par[2] / par[4]
                         * exp( sum * rho_1 );
                }
                sum *= par[5];
                
                // assume Poisson fluctuations (neglecting background noise)
                // (Blobel p.196)
                if( n > 0. && sum > 0. )
                {
                    LL += n * log( sum ) - sum - n * log( n ) + n;
                }
                else
                {
                    LL += -1. * sum;
                }
                // Gauss with background noise
                //             iVar = sum+pedvars*pedvars;
                //             iVar = sum;
                //             LL += -0.5*log(2. * TMath::Pi()) - log(sqrt(iVar)) - (n-sum)*(n-sum)/2./iVar;
                // time gradient analysis
                // (line fit to time gradient)
                // - affects image orientation and centroid position
                if( iImageCalculation->minimize_time_gradient_for_this_event() )
                {
                    t = iImageCalculation->getLLT()[i];
                    if( t > 0. )
                    {
                        // calculate position along long axis
                        tx = (x-par[1]) * cos(par[0]) + (y-par[3]) * sin(par[0]);
                        LL_t += 
                            - 1./2./(t_sig*t_sig) 
                            * (t - par[6] - par[7] * tx )
                            * (t - par[6] - par[7] * tx );
                    }
                }

            }
        }
    }
    LL_t += t_sig_term;
    f = -1. * ( LL + LL_t );
}


/*
   fill pixel sums used during fitting period

   --> select only image, border, and image neighbour pixels
*/
double VImageParameterFitter::fill_pixel_sums( bool iUseSums2 )
{
    // pixel list
    fll_X.clear();
    fll_Y.clear();
    fll_R.clear();
    fll_Sums.clear();
    fll_T.clear();
    // will be true if sum in pixel is estimated by fit
    fLLEst.assign( fData->getSums().size(), false );
    // for limits of parameter fitting
    fdistXmin =  1.e99;
    fdistXmax = -1.e99;
    fdistYmin =  1.e99;
    fdistYmax = -1.e99;
    double i_sumMax = -1.e99;
    // loop over all pixels
    for( unsigned int j = 0; j < fData->getSums().size(); j++ )
    {
        // ignore dead channels
        if( fData->getDead( j, fData->getHiLo()[j] ) )
        {
            continue;
        }
        // only image/border pixels and their neighbours are used in the fit
        if( j < fData->getImageBorderNeighbour().size() && fData->getImageBorderNeighbour()[j] )
        {
            // pixel position in the camera
            double xi = getDetectorGeometry()->getX()[j];
            double yi = getDetectorGeometry()->getY()[j];
            double ri = fData->getDetectorGeo()->getTubeRadius()[j];
            fll_X.push_back( xi );
            fll_Y.push_back( yi );
            fll_R.push_back( ri );
            
            // get an estimate of the distance (limits to the fit)
            if( xi > fdistXmax ) fdistXmax = xi;
            if( xi < fdistXmin ) fdistXmin = xi;
            if( yi > fdistYmax ) fdistYmax = yi;
            if( yi < fdistYmin ) fdistYmin = yi;
            
            if( fData->getImage()[j] || fData->getBorder()[j] )
            {
                if( iUseSums2 )
                {
                    fll_Sums.push_back( fData->getSums2()[j] );
                }
                else
                {
                    fll_Sums.push_back( fData->getSums()[j] );
                }
                // LL fit using time gradient only for image
                // with GEO detected time gradient
                if( minimize_time_gradient_for_this_event() )
                {
                    fll_T.push_back( fData->getPulseTime()[j] );
                }
                else
                {
                    fll_T.push_back( -99. );
                }
            }
            else
            {
                fll_Sums.push_back( 0.1 );
                fll_T.push_back( -999. );
            }
            // image weighting with squared intensity
            // (non standard from traditional image calculation!)
            if( fData->getRunParameter() && fData->getRunParameter()->fSquaredImageCalculation )
            {
                fll_Sums.back() = fll_Sums.back() * fll_Sums.back();
            }
            if( fll_Sums.back() > i_sumMax )
            {
                i_sumMax = fll_Sums.back();
            }
        }
    }
    if( fLLDebug )
    {
        cout << "FLL FITTER limits:  xmax: " << fdistXmax << "  xmin: " << fdistXmin;
        cout << " ymax: " << fdistYmax << " ymin: " << fdistYmin << " #: " << fll_Sums.size() << endl;
    }
    return i_sumMax;
}

void VImageParameterFitter::resetFitParameters()
{
    rho = 0.;
    drho = 0.;
    cen_x = 0.;
    dcen_x = 0.;
    sigmaX = 0.;
    dsigmaX = 0.;
    cen_y = 0.;
    dcen_y = 0.;
    sigmaY = 0.;
    dsigmaY = 0.;
    dsignal = 0.;
    toffset = 0.;
    dtoffset = 0.;
    tgrad = 0.;
    dtgrad = 0.;
    tchi2 = 0.;
    dtchi2 = 0.;
    phi = 0.;
    dphi = 0.;
    if( fParGeo )
    {
        cen_x = fParGeo->cen_x;
        sigmaX = fParGeo->sigmaX;
        cen_y = fParGeo->cen_y;
        sigmaY = fParGeo->sigmaY;
        phi = fParGeo->phi;
        // fit will fail for width < this value
        if( fParGeo->width < 0.001 )
        {
            sigmaX = 0.1;
            sigmaY = 0.1;
        }
        // must be called last!
        rho = getLL_startingvalue_rho();
    }
}

void VImageParameterFitter::defineFitParameters()
{
    // don't know if this step parameter is the best
    double step = 1.e-4;
    fLLFitter->Release( 0 );
    fLLFitter->Release( 1 );
    fLLFitter->Release( 2 );
    fLLFitter->Release( 3 );
    fLLFitter->Release( 4 );
    if( bRotatedNormalDistributionFit )
    {
        fLLFitter->DefineParameter( 0, "phi",
                                    phi,
                                    step,
                                   -40.*M_PI,
                                    40.*M_PI );
    }
    else
    {
        fLLFitter->DefineParameter( 0, "rho",
                                    rho,
                                    step,
                                    getLL_paramameterlimits_rho( rho, -1. ),
                                    getLL_paramameterlimits_rho( rho, 1. ) );
    }
    // cen_x starting values
    fLLFitter->DefineParameter( 1, "meanX", 
                                cen_x,
                                step, 
                                getLL_paramameterlimits_cen( fdistXmin, cen_x, sigmaX, -1. ),
                                getLL_paramameterlimits_cen( fdistXmin, cen_x, sigmaX, 1. ) );
    // sigma x starting values
    double sigma_x_L = sigmaX / 4.;
    double sigma_x_U = 2.*sigmaX+ 1.;
    if( bRotatedNormalDistributionFit )
    {
        sigma_x_L = fParGeo->width * 0.5;
        sigma_x_U = fParGeo->length * 2.5;
    }
    fLLFitter->DefineParameter( 2, "sigmaX",
                                sigmaX,
                                step, 
                                sigma_x_L,
                                sigma_x_U );
    // cen_y starting values
    fLLFitter->DefineParameter( 3, "meanY", 
                                cen_y,
                                step, 
                                getLL_paramameterlimits_cen( fdistYmin, cen_y, sigmaY, -1. ),
                                getLL_paramameterlimits_cen( fdistYmin, cen_y, sigmaY, 1. ) );
    // sigma x starting values
    double sigma_y_L = sigmaY / 4.;
    double sigma_y_U = 2.*sigmaY + 1.;
    if( bRotatedNormalDistributionFit )
    {
        sigma_y_L = fParGeo->width * 0.5;
        sigma_y_U = fParGeo->length * 2.5;
    }
    fLLFitter->DefineParameter( 4, "sigmaY", 
                                sigmaY, 
                                step, 
                                sigma_y_L,
                                sigma_y_U );

    // signal
    fLLFitter->DefineParameter( 5, "signal", signal, step, 0., 1.e6 );
    // time gradient analysis
    if( fData->getRunParameter()->fMinimizeTimeGradient )
    {
        double toffset_start = 0.;
        if( fData->getXGraph( false ) )
        {
             toffset_start = fData->getXGraph( false )->Eval( 0. );
        }
        fLLFitter->DefineParameter( 6, "toffset", toffset_start, step, toffset_start-10., toffset_start+10. );
        fLLFitter->DefineParameter( 7, "tgrad", fParGeo->tgrad_x, step, -1.e3, 1.e3 );
        fLLFitter->DefineParameter( 8, "tchi2", sqrt( fParGeo->tchisq_x ), step, 0.1, 2.*sqrt( fParGeo->tchisq_x ) );
        if( minimize_time_gradient_for_this_event() )
        {
            fLLFitter->Release( 6 );
            fLLFitter->Release( 7 );
        }
        else
        {
            fLLFitter->FixParameter( 6 );
            fLLFitter->FixParameter( 7 );
        }
        fLLFitter->FixParameter( 8 );
    }
}
    
void VImageParameterFitter::getFitStatistics()
{
    double edm = 0.;
    double amin = 0.;
    double errdef = 0.;
    int nvpar = 0;
    int nparx = 0;
    int nstat = 0;
    fLLFitter->mnstat( amin, edm, errdef, nvpar, nparx, nstat );
    fParLL->Fitstat = nstat;
    fParLL->Fitmin = amin;
    fParLL->Fitedm = edm;
    
    if( fLLDebug )
    {
        cout << "FLLFITTER STAT " << nstat << endl;
    }
}

void VImageParameterFitter::getFitResults()
{
    if( bRotatedNormalDistributionFit )
    {
        fLLFitter->GetParameter( 0, phi, dphi );
    }
    else
    {
        fLLFitter->GetParameter( 0, rho, drho );
    }
    fLLFitter->GetParameter( 1, cen_x, dcen_x );
    fLLFitter->GetParameter( 2, sigmaX, dsigmaX );
    fLLFitter->GetParameter( 3, cen_y, dcen_y );
    fLLFitter->GetParameter( 4, sigmaY, dsigmaY );
    fLLFitter->GetParameter( 5, signal, dsignal );
    if( minimize_time_gradient_for_this_event() )
    {
        fLLFitter->GetParameter( 6, toffset, dtoffset );
        fLLFitter->GetParameter( 7, tgrad, dtgrad );
        fLLFitter->GetParameter( 8, tchi2, dtchi2 );
    }
    if( fLLDebug )
    {
        cout << "=======================================================================" << endl;
        cout << "Telescope " << fData->getTelID()+1 << " starting parameters" << endl;
        cout << "\t FITSTAT " << fParLL->Fitstat << endl;
        cout << "\t CENX " << cen_x << " +- " << dcen_x << endl;
        cout << "\t CENY " << cen_y << " +- " << dcen_y << endl;
        cout << "\t SIGX " << sigmaX << " +- " << dsigmaX << endl;
        cout << "\t SIGY " << sigmaY << " +- " << dsigmaY << endl;
        cout << "\t F signal " << signal << endl;
        if( minimize_time_gradient_for_this_event() )
        {
            cout << "\t TOFF " << toffset << " +- " << dtoffset << "\t" << endl;
            cout << "\t TGRAD " << tgrad << " +- " << dtgrad << "\t" << fParGeo->tgrad_x << endl;
            cout << "\t TCHI2 " << tchi2 << " +- " << dtchi2 << "\t" << fParGeo->tchisq_x << "\t" << sqrt( fParGeo->tchisq_x ) << endl;
        }
    }
    sigmaX = fabs( sigmaX );
    sigmaY = fabs( sigmaY );
    if( sigmaX < ZeroTolerence )
    {
        sigmaX = 0.;
    }
    if( sigmaY < ZeroTolerence )
    {
        sigmaY = 0.;
    }
    // copy parameters to output tree
    fParLL->rho = rho;
    fParLL->phi = phi;
    fParLL->cen_x = cen_x;
    fParLL->sigmaX = sigmaX;
    fParLL->cen_y = cen_y;
    fParLL->sigmaY = sigmaY;
    fParLL->signal = signal;
    fParLL->drho = drho;
    fParLL->dphi = dphi;
    fParLL->dcen_x = dcen_x;
    fParLL->dsigmaX = dsigmaX;
    fParLL->dsigmaY = dsigmaY;
    fParLL->dcen_y = dcen_y;
    fParLL->dsignal = dsignal;
}

void VImageParameterFitter::calculateImageParameters( bool iUseSums2,
                                                      bool iEqualSummationWindows )
{
    // image size
    calculate_image_size( iUseSums2, iEqualSummationWindows );

    if( bRotatedNormalDistributionFit )
    {
        if( sigmaX < sigmaY )
        {
            fParLL->width = sigmaX;
            fParLL->dwidth = dsigmaX;
            fParLL->length = sigmaY;
            fParLL->dlength = dsigmaY;
            phi = redang( fParLL->phi + M_PI/2.,  2. * M_PI );
            fParLL->phi = phi;
        }
        else
        {
            fParLL->width = sigmaY;
            fParLL->dwidth = dsigmaY;
            fParLL->length = sigmaX;
            fParLL->dlength = dsigmaX;
        }
        sigmaX *= cos(phi);
        sigmaY *= sin(phi);
        calculate_image_rho();
    }

    // covariance
    double sigmaXY = rho * sqrt( sigmaX * sigmaX * sigmaY * sigmaY );
    
    double d = sigmaY * sigmaY - sigmaX * sigmaX;
    double z = sqrt( d * d + 4. * sigmaXY * sigmaXY );

    double dsxxy2 = sigmaX * sigmaX * sigmaY * sigmaY;
    double dsigmaXY2 = dsxxy2 * drho * drho;
    dsigmaXY2 += rho * rho * sigmaX * sigmaX * sigmaY * sigmaY * sigmaY * sigmaY / dsxxy2 * dsigmaX * dsigmaX;
    dsigmaXY2 += rho * rho * sigmaX * sigmaX * sigmaX * sigmaX * sigmaY * sigmaY / dsxxy2 * dsigmaY * dsigmaY;
    double dd2 = 4.*sigmaY * sigmaY * dsigmaY * dsigmaY 
               + 4.*sigmaX * sigmaX * dsigmaX * dsigmaX;
    double dz2 = d * d / ( d * d + 4.*sigmaXY * sigmaXY ) * dd2
              + 16.*sigmaXY * sigmaXY / ( d * d + 4.*sigmaXY * sigmaXY ) * dsigmaXY2;

    if( !bRotatedNormalDistributionFit )
    {
        calculate_image_length( z, dz2 );
        calculate_image_width( z, dz2 );
        calculate_image_phi( dsxxy2 );
    }
    fParLL->f_d = d;
    fParLL->f_s = z;
    fParLL->f_sdevxy = sigmaXY;
    calculate_image_distance();

    //  fit was not successfull if width is close to zero -> take pargeo parameters
    if( fParLL->width < ZeroTolerence )
    {
        fParLL->width  = fParGeo->width;
        fParLL->length = fParGeo->length;
        fParLL->cen_x  = fParGeo->cen_x;
        fParLL->cen_y  = fParGeo->cen_y;
        fParLL->phi    = fParGeo->phi;
        fParLL->size = fParGeo->size;
        fParLL->size2 = fParGeo->size2;
    }

    // setters
    fParLL->los = 0.;
    if( fParLL->size != 0. )
    {
        fParLL->los = fParLL->length / fParLL->size;
    }
    if( fData->getRunParameter()->fMinimizeTimeGradient )
    {
        fParLL->tgrad_x = tgrad;
    }
    fParLL->azwidth = -1.;
    fParLL->dazwidth = 0.;
    fParLL->ntRec = 0;
    fParLL->ntubes = fParGeo->ntubes;
    fParLL->bad = fParGeo->bad;
}

/*    estimate image size
     integrate over all bins to get fitted size
*/
void VImageParameterFitter::calculate_image_size( bool iUseSums2,
                                                  bool iEqualSummationWindows )
{
    bool fWidthResetted = false;
    float  iSize = fParGeo->size;
    if( iUseSums2 )
    {
        iSize = fParGeo->size2;
    }
    // require a successfull fit
    if( fParLL->Fitstat > 0 )
    {
        // assume that all pixel are of the same geometrical size
        unsigned int iCentreTube = fData->getDetectorGeometry()->getCameraCentreTubeIndex();
        if( iCentreTube < 9999 )
        {
            // make sure that width is not close to zero (can happen when all pixels are on a line)
            if( fData->getDetectorGeometry()->getTubeRadius().size() > iCentreTube
             && fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] > 0. )
            {
                if( sigmaX < fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] * 0.1 )
                {
                    sigmaX = fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] * 0.1;
                    fWidthResetted = true;
                }
                if( sigmaY < fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] * 0.1 )
                {
                    sigmaY = fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] * 0.1;
                    fWidthResetted = true;
                }
            }
            // recentre fitted image to camera centre
            double cen_x_recentered = 0.;
            double cen_y_recentered = 0.;
            if( fData->getDetectorGeometry()->getTubeRadius().size() > iCentreTube
             && fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] > 0. )
            {
                cen_x_recentered = cen_x - fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] * 2.
                                   * ( int )( cen_x / fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] / 2. );
                cen_y_recentered = cen_y - fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] * 2.
                                   * ( int )( cen_y / fData->getDetectorGeometry()->getTubeRadius()[iCentreTube] / 2. );
            }
            if( fLLDebug )
            {
                cout << "FLLFITTER RECENTERED " << cen_x_recentered << "\t" << cen_y_recentered << endl;
            }
            ////////////////////////////////////////////////////////////////////
            // calculate new size from fit function (not from measured charges!)
            iSize = 0.;
            if( fLLDebug )
            {
                cout << "FLLFITTER SIZE CALCULATION ";
                cout << rho << "\t" << cen_x_recentered << "\t" << sigmaX << "\t" << cen_y_recentered << "\t";
                cout << sigmaY << "\t" << signal << endl;
            }
            for( unsigned int i = 0; i < fData->getSums().size(); i++ )
            {
                iSize +=  calculatePixelBrightness( i, rho, phi,
                                                       cen_x_recentered, sigmaX, 
                                                       cen_y_recentered, sigmaY, 
                                                       signal );
            }
        }
    }
    fParLL->size2LL = fParGeo->size2;
    fParLL->sizeLL =  fParGeo->size;
    if( iUseSums2 || iEqualSummationWindows )
    {
        if( fWidthResetted )
        {
            fParLL->size2 = fParGeo->size2;
        }
        else
        {
            fParLL->size2 = iSize;
        }
        fParLL->size = fParLL->size2;
    }
    // set size for !sums2
    else
    {
        if( fWidthResetted )
        {
            fParLL->size = fParGeo->size;
        }
        else
        {
            fParLL->size = iSize;
        }
        fParLL->size2 = fParGeo->size2;
    }
}      

void VImageParameterFitter::calculate_image_length( double z, double dz2 )
{
    double length = sqrt( 0.5 * ( sigmaX * sigmaX + sigmaY * sigmaY + z ) );

    double dlength = 0.;
    dlength  = 4.*sigmaX * sigmaX * dsigmaX * dsigmaX;
    dlength += 4.*sigmaY * sigmaY * dsigmaY * dsigmaY;
    dlength += 1. / 2. * dz2;
    dlength *= 1. / 8. / ( sigmaX * sigmaX + sigmaY * sigmaY + z );
    if( dlength > 0. )
    {
        dlength = sqrt( dlength );
    }

    fParLL->length = length;
    fParLL->dlength = dlength;
}

void VImageParameterFitter::calculate_image_width( double z, double dz2 )
{
    double width = sigmaX * sigmaX + sigmaY * sigmaY - z;
    if( width > ZeroTolerence )
    {
        width = sqrt( 0.5 * width );
    }
    double dwidth = 0.;
    dwidth  = 4.*sigmaX * sigmaX * dsigmaX * dsigmaX;
    dwidth += 4.*sigmaY * sigmaY * dsigmaY * dsigmaY;
    dwidth += dz2;
    dwidth *= 1. / 8. / ( sigmaX * sigmaX + sigmaY * sigmaY - z );
    if( dwidth > 0. )
    {
        dwidth  = sqrt( dwidth );
    }
    fParLL->width = width;
    fParLL->dwidth = dwidth;
}

void VImageParameterFitter::calculate_image_phi( double dsxxy2 )
{
    double dsx_y2 = sigmaX * sigmaX - sigmaY * sigmaY;

    phi = atan2( 2.*rho * sigmaX * sigmaY, dsx_y2 ) / 2.;

    if( dsx_y2 == 0. )
    {
        dphi = 1000.;
    }
    else
    {
        dphi  = 4.* dsxxy2 / dsx_y2 / dsx_y2 * drho * drho;
        dphi += ( 2.*rho * sigmaY * dsx_y2 - 4.*rho * sigmaY * sigmaX * sigmaX ) * ( 2.*rho * sigmaY * dsx_y2 - 4.*rho * sigmaY * sigmaX * sigmaX ) / dsx_y2 / dsx_y2 / dsx_y2 / dsx_y2 * dsigmaX * dsigmaX;
        dphi += ( 2.*rho * sigmaX * dsx_y2 + 4.*rho * sigmaY * sigmaY * sigmaX ) * ( 2.*rho * sigmaX * dsx_y2 + 4.*rho * sigmaY * sigmaY * sigmaX ) / dsx_y2 / dsx_y2 / dsx_y2 / dsx_y2 * dsigmaY * dsigmaY;
        dphi *= 1. / ( 1. + 2.*rho * sigmaX * sigmaY / dsx_y2 * 2.*rho * sigmaX * sigmaY / dsx_y2 );
        dphi *= 1. / ( 1. + 2.*rho * sigmaX * sigmaY / dsx_y2 * 2.*rho * sigmaX * sigmaY / dsx_y2 );
        dphi  = sqrt( dphi );
        dphi  = redang( dphi, M_PI / 2. );
   }
   fParLL->phi = phi;
   fParLL->cosphi = cos( fParLL->phi );
   fParLL->sinphi = sin( fParLL->phi );
   fParLL->dphi = dphi;
}

void VImageParameterFitter::calculate_image_distance()
{
    double dist = sqrt( cen_x * cen_x + cen_y * cen_y );
    double ddist = 0.;
    if( dist != 0. )
    {
        ddist = cen_x * cen_x / dist / dist * dcen_x * dcen_x + cen_y * cen_y / dist / dist * dcen_y * dcen_y;
    }
    ddist = sqrt( ddist );
    double miss = fabs( cen_y * fParLL->cosphi - cen_x * fParLL->sinphi );
    double dmiss = 0.;
    dmiss  = fParLL->cosphi * fParLL->cosphi * dcen_x * dcen_x
           + fParLL->sinphi * fParLL->sinphi * dcen_y * dcen_y;
    dmiss += fParLL->dphi * fParLL->dphi * ( fParLL->sinphi * fParLL->sinphi * cen_x * cen_x
                                           + fParLL->cosphi * fParLL->cosphi * cen_y * cen_y );
    fParLL->miss = miss;
    fParLL->dmiss = sqrt( dmiss );

    fParLL->alpha = 0.;
    fParLL->dalpha = 0.;
    if( dist != 0. )
    {
        fParLL->alpha = asin( miss / dist ) * TMath::RadToDeg();;
        double dalpha = 0.;
        dalpha  = dmiss * dmiss / dist / dist + miss * miss / dist / dist / dist / dist * ddist * dist;
        dalpha *= 1. / sqrt( 1. - miss * miss / dist / dist );
        dalpha = redang( sqrt( dalpha ), 2. * M_PI );
        fParLL->dalpha = dalpha * TMath::RadToDeg();
    }

    fParLL->dist = dist;
    fParLL->ddist = ddist;
}

void VImageParameterFitter::calculate_image_rho()
{
    rho = 0.;
    if( sigmaY > 0. && sigmaX > 0. )
    {
        rho = tan(2.*phi) * (sigmaX*sigmaX-sigmaY*sigmaY);
        rho = rho / (2. * sigmaY * sigmaX);
    }
    drho = 0.;

    fParLL->rho = rho;
    fParLL->drho = drho;
}


/*
    2D Normal distribution

    function used to integrate probability density
    over size of pixel

    optional, used only if fNormal2D != 0

*/
Double_t normal2DRotated( Double_t *x, Double_t *par )
{
    double f  = 0.;
    if( par[1] > 0. && par[3] > 0. )
    {
        f += (x[0] - par[0]) * (x[0] - par[0]) / par[1] / par[1];
        f += (x[1] - par[2]) * (x[1] - par[2]) / par[3] / par[3];
        f *= -1. / 2.;
        f  = exp( f );
        f *= 1. / 2. / TMath::Pi() / par[1] / par[3];
    }
    return f;
}
