/*! \class VImageParameterCalculation
    \brief calculation of the Hillas image parameters

    dead channels for low gains and LL

*/

#include <VImageParameterCalculation.h>

VImageParameterCalculation::VImageParameterCalculation( unsigned int iShortTree, VEvndispData* iData )
{
    fDebug = false;
    if( fDebug )
    {
        cout << "VImageParameterCalculation::VImageParameterCalculation()" << endl;
    }
    fData = iData;
    fParGeo = new VImageParameter( iShortTree, fData->getRunParameter()->fWriteImagePixelList );
    fParLL =  new VImageParameter( iShortTree, fData->getRunParameter()->fWriteImagePixelList );
    fboolCalcGeo = false;
    fboolCalcTiming = false;
    fDetectorGeometry = 0;
    fHoughTransform = 0;
    
    fImageFitter = new VImageParameterFitter( fData,
            fDebug );
}


VImageParameterCalculation::~VImageParameterCalculation()
{
    if( fDebug )
    {
        cout << "VImageParameterCalculation::~VImageParameterCalculation()" << endl;
    }
    delete fParGeo;
    delete fParLL;
    delete fImageFitter;
}


/*!
   \param iVmode Minuit print mode, 1 = quit, 2 - verbose

*/
void VImageParameterCalculation::initMinuit( int iVmode )
{
    fImageFitter->initMinuit( iVmode );
}

/*

   calculate timing parameters
   (time gradient along long axis of image)

*/
void VImageParameterCalculation::calcTimingParameters( bool iIsSecondPass )
{
    if( fDebug )
    {
        cout << "VImageParameterCalculation::calcTimingParameters" << endl;
    }
    if( !fboolCalcGeo )
    {
        return;
    }
    if( fData->getPulseTime().size() == 0 )
    {
        return;
    }
    
    vector< double > xpos;
    vector< double > t;
    vector< double > et;
    vector< bool > usePoint;
    double i_xmin = 1.e99;
    double i_xmax = -1.e99;
    // get average pulse times
    // depend on parameter set in getSumWindowStart_T_method():
    // 1: tzeros
    // 2: average pulse times
    // 3: trigger pixel time
    for( unsigned int i = 0; i < fData->getPulseTime().size(); i++ )
    {
        //  use only image and border and no low-gain channels (timing is different)
        // if( ( fData->getImage()[i] || fData->getBorder()[i] ) && fData->getPulseTime()[i] > 0. && !fData->getHiLo()[i] )
        if( ( fData->getImage()[i] || fData->getBorder()[i] ) && fData->getPulseTime()[i] > -998 && !fData->getHiLo()[i] )
        {
            double xi = getDetectorGeometry()->getX()[i];
            double yi = getDetectorGeometry()->getY()[i];
            // loop over image tubes
            double xpmt = xi - fParGeo->cen_x;
            double ypmt = yi - fParGeo->cen_y;
            // position along the major axis of the image (relative to centroid position)
            xpos.push_back( xpmt * fParGeo->cosphi + ypmt * fParGeo->sinphi );
            if( xpos.back() > i_xmax )
            {
                i_xmax =  xpos.back();
            }
            if( xpos.back() < i_xmin )
            {
                i_xmin =  xpos.back();
            }
            // timing parameter
            t.push_back( fData->getPulseTime()[i] );
            // timing error
            et.push_back( 0. );
            // (GM) (used before 20130115, not the default anymore)
            if( fData->getDoublePassErrorWeighting2005() )
            {
                //  timing resolution from variable laser pulse studies (run 751)
                et.back() = 13.0 * exp( -0.035 * ( fData->getSums()[i] + 30. ) ) + fData->getTOffsetvars()[i];
                // make that the timing resolution is not too small (important for MC)
                if( et.back() < 5.e-2 )
                {
                    et.back() = 0.3;
                }
            }
            // used now (note scaling in a later step)
            else
            {
                // use 1./sums as error for fitting (rescale errors later)
                if( fData->getSums()[i] > 0. )
                {
                    et.back() = 1. / fData->getSums()[i];
                }
                else
                {
                    et.back() = 0.3;
                }
            }
            // use this point
            usePoint.push_back( true );
            // min/max/mean times
            if( fData->getPulseTime()[i]  < fParGeo->tmin )
            {
                fParGeo->tmin = fData->getPulseTime()[i];
            }
            if( fData->getPulseTime()[i]  > fParGeo->tmax )
            {
                fParGeo->tmax = fData->getPulseTime()[i];
            }
            fParGeo->tmean += fData->getPulseTime()[i];
        }
    }
    
    ///////////////////////////////////////////////////
    // find and remove outliers
    if( xpos.size() > 9 && i_xmin < i_xmax )
    {
        TH1F h( "houtlier", "", 100, i_xmin, i_xmax );
        for( unsigned int i = 0; i < xpos.size(); i++ )
        {
            h.Fill( xpos[i] );
        }
        double i_a[] = { 0.25, 0.75 };
        double i_b[] = { 0.00, 0.00 };
        h.GetQuantiles( 2, i_b, i_a );
        double iqr = i_b[1] - i_b[0];
        double i_limin_min = i_b[0] - 1.5 * iqr;
        double i_limin_max = i_b[1] + 1.5 * iqr;
        for( unsigned int i = 0; i < xpos.size(); i++ )
        {
            if( xpos[i] < i_limin_min )
            {
                usePoint[i] = false;
            }
            else if( xpos[i] > i_limin_max )
            {
                usePoint[i] = false;
            }
        }
    }
    ///////////////////////////////////////////////////
    // rescale errors to max 2 samples
    // (do think that we can reconstruct the
    //  best integration window to 2 samples
    //  or better)
    double iEmax = 0.;
    for( unsigned int i = 0; i < xpos.size(); i++ )
    {
        if( usePoint[i] )
        {
            if( et[i] > iEmax )
            {
                iEmax = et[i];
            }
        }
    }
    if( iEmax > 0. )
    {
        for( unsigned int i = 0; i < xpos.size(); i++ )
        {
            if( usePoint[i] )
            {
                et[i] *= 2. / iEmax;
            }
        }
    }
    
    TGraphErrors* xgraph = fData->getXGraph( iIsSecondPass );
    if( xgraph )
    {
        int nclean = 0;
        for( unsigned int i = 0; i < usePoint.size(); i++ )
        {
            if( usePoint[i] )
            {
                nclean++;
            }
        }
        xgraph->Set( nclean );
        
        if( nclean > 2 )
        {
            // Fill the graphs for long (x) short(y) and radial (r) axis
            int z = 0;
            for( unsigned int i = 0; i < xpos.size(); i++ )
            {
                if( usePoint[i] )
                {
                    xgraph->SetPoint( z, xpos[i], t[i] );
                    xgraph->SetPointError( z, 0., et[i] );
                    z++;
                }
            }
            // ignore all the clutter from the tgraph fitting algorithm
            gErrorIgnoreLevel = 5000;
            // robust fitting for larger images
            // (only for first pass, which is relevant for the window placement)
            // robust fitting seems to be unstable for >1000 points
            if( z > 9 && z < 1000 && !iIsSecondPass )
            {
                xgraph->Fit( "pol1", "rob Q" );
            }
            else
            {
                xgraph->Fit( "pol1", "Q" );
            }
            
            TF1* xline = xgraph->GetFunction( "pol1" );
            
            gErrorIgnoreLevel = 0;
            
            // fill fit results
            fParGeo->tint_x   = xline->GetParameter( 0 );
            fParGeo->tgrad_x  = xline->GetParameter( 1 );
            fParGeo->tint_dx  = xline->GetParError( 0 );
            fParGeo->tgrad_dx = xline->GetParError( 1 );
            fParGeo->tchisq_x = xline->GetChisquare();
            if( nclean > 0. )
            {
                fParGeo->tmean = fParGeo->tmean / nclean;
            }
            else
            {
                fParGeo->tmean = -9999.;
            }
            fData->setXGraph( xgraph, iIsSecondPass );
        }
    }
    else
    {
        fParGeo->tint_x = -999.;
        fParGeo->tgrad_x = -999.;
        fParGeo->tint_dx = -999.;
        fParGeo->tgrad_dx = -999.;
        fParGeo->tchisq_x = -999.;
        fParGeo->tmin = -999.;
        fParGeo->tmax = -999.;
        fParGeo->tmean = -999.;
    }
    fboolCalcTiming = true;
}

/*****************************************************************************
muonRingFinder
input: pointer to pixels passing cleaning
output: return the likely location of x0, y0, and R with sigmaR
****************************************************************************
*/
void VImageParameterCalculation::muonRingFinder()
{
    if( fDebug )
    {
        cout << "VImageParameterCalculation::muonRingFinder" << endl;
    }
    if( fData->getSums().size() == 0 )
    {
        fParGeo->muonX0 = 0.0;
        fParGeo->muonY0 = 0.0;
        fParGeo->muonXC = 0.0;
        fParGeo->muonYC = 0.0;
        fParGeo->muonRadius = 0.0;
        fParGeo->muonRSigma = 0.0;
        return;
    }
    
    if( !getDetectorGeometry() )
    {
        cout << "VImageParameterCalculation::muonRingFinder error: detector geometry not defined" << endl;
        exit( EXIT_FAILURE );
    }
    
    int counter = 0, safty = 0, noChangeX = 0, noChangeY = 0;
    unsigned int i;
    double x0[2] =
    {
        0.0
    }
    , y0[2] =
    {
        0.0
    }
    , rBar[2] =
    {
        0.0
    }
    , rVariance[2] =
    {
        20.0
    };
    double rTotal = 0.0, rSquaredTotal = 0.0, tmp;
    double xi, yi;
    
    //calculate r to each point in the filtered binary image & thereby mean of r
    for( i = 0; i < fData->getSums().size(); i++ )
    {
        if( fData->getImage()[i] || fData->getBorder()[i] )
        {
            counter++;
            xi = getDetectorGeometry()->getX()[i];
            yi = getDetectorGeometry()->getY()[i];
            tmp = sqrt( pow( ( xi - x0[0] ), 2 ) + pow( ( yi - y0[0] ) , 2 ) );
            rTotal += tmp;
            rSquaredTotal += tmp * tmp;
        }
    }
    
    if( counter > 10 )
    {
        rVariance[0] = ( rSquaredTotal - rTotal * rTotal / counter ) / counter;
        rBar[0] = rTotal / counter;
        rBar[1] = 20.0;                           //so we have a starting value to enter the loop
    }
    else
    {
        fParGeo->muonX0 = 0.0;
        fParGeo->muonY0 = 0.0;
        fParGeo->muonXC = 0.0;
        fParGeo->muonYC = 0.0;
        fParGeo->muonRadius = 0.0;
        fParGeo->muonRSigma = 0.0;
        return;
    }
    ///****************************************************
    // LOOP through until good x0, y0 coordinates are found
    ///****************************************************
    while( ( ( noChangeX == 0 && noChangeY == 0 ) || safty < 20 ) && safty < 100 )
    {
        safty++;
        ///*********************************************
        // TAKE A STEP IN THE X-DRIRECTION ?
        ///*********************************************
        //try a different x0[1], see if sigma of r decreases
        x0[1] = x0[0] + 0.1 / ( 1 + pow( safty, .3 ) );
        y0[1] = y0[0];
        noChangeX = 0;                            //reset the end switch
        
        //calculate r to each point in the filtered binary image
        rTotal = 0;
        rSquaredTotal = 0;
        tmp = 0;
        for( i = 0; i < fData->getSums().size(); i++ )
        {
            if( fData->getImage()[i] || fData->getBorder()[i] )
            {
                xi = getDetectorGeometry()->getX()[i];
                yi = getDetectorGeometry()->getY()[i];
                tmp = sqrt( pow( ( xi - x0[1] ), 2 ) + pow( ( yi - y0[1] ) , 2 ) );
                rTotal += tmp;
                rSquaredTotal += tmp * tmp;
            }
        }
        rVariance[1] = ( rSquaredTotal - rTotal * rTotal / counter ) / counter;
        rBar[1] = rTotal / counter;
        
        //is rVariance[1] > rVariance[0] ? then try a step in the -x direction
        if( rVariance[1] > rVariance[0] )
        {
            //	x0[1] =  x0[0]-0.03/(1+pow(safty,.3));       //try a different x0[1], see if mean of r decreases
            x0[1] =  x0[0] - 0.05;                //try a different x0[1], see if mean of r decreases
            rTotal = 0;
            rSquaredTotal = 0;
            tmp = 0;
            for( i = 0; i < fData->getSums().size(); i++ )
            {
                if( fData->getImage()[i] || fData->getBorder()[i] )
                {
                    xi = getDetectorGeometry()->getX()[i];
                    yi = getDetectorGeometry()->getY()[i];
                    tmp = sqrt( pow( ( xi - x0[1] ), 2 ) + pow( ( yi - y0[1] ) , 2 ) );
                    rTotal += tmp;
                    rSquaredTotal += tmp * tmp;
                }
            }
            rVariance[1] = ( rSquaredTotal - rTotal * rTotal / counter ) / counter;
            rBar[1] = rTotal / counter;
            
            //is rVariance[1] > rVariance[0] ? then keep original coordinate
            if( rVariance[1] > rVariance[0] )
            {
                x0[1] = x0[0];
                rBar[1] = rBar[0];
                rVariance[1] = rVariance[0];
                noChangeX = 1;
            }
        }
        
        rBar[0] = rBar[1];
        rVariance[0] = rVariance[1];              // why was this commented out??
        x0[0] = x0[1];
        
        ///*********************************************
        // TAKE A STEP IN THE Y-DIRECTION ?
        ///*********************************************
        //try a different x0[1], see if sigma of r decreases
        x0[1] = x0[0];
        y0[1] = y0[0] + 0.1 / ( 1 + pow( safty, .3 ) );
        noChangeY = 0;                            //reset the end switch
        
        //calculate r to each point in the filtered binary image
        rTotal = 0;
        rSquaredTotal = 0;
        tmp = 0;
        for( i = 0; i < fData->getSums().size(); i++ )
        {
            if( fData->getImage()[i] || fData->getBorder()[i] )
            {
                xi = getDetectorGeometry()->getX()[i];
                yi = getDetectorGeometry()->getY()[i];
                tmp = sqrt( pow( ( xi - x0[1] ), 2 ) + pow( ( yi - y0[1] ) , 2 ) );
                rTotal += tmp;
                rSquaredTotal += tmp * tmp;
            }
        }
        rVariance[1] = ( rSquaredTotal - rTotal * rTotal / counter ) / counter;
        rBar[1] = rTotal / counter;
        
        //is rVariance[1] > rVariance[0] ? then try a step in the -x direction
        if( rVariance[1] > rVariance[0] )
        {
            //	y0[1] =  y0[0]-0.03/(1+pow(safty,.3));       //try a different x0[1], see if mean of r decreases
            y0[1] =  y0[0] - 0.05;                //try a different x0[1], see if mean of r decreases
            rTotal = 0;
            rSquaredTotal = 0;
            for( i = 0; i < fData->getSums().size(); i++ )
            {
                if( fData->getImage()[i] || fData->getBorder()[i] )
                {
                    xi = getDetectorGeometry()->getX()[i];
                    yi = getDetectorGeometry()->getY()[i];
                    tmp = sqrt( pow( ( xi - x0[1] ), 2 ) + pow( ( yi - y0[1] ) , 2 ) );
                    rTotal += tmp;
                    rSquaredTotal += tmp * tmp;
                }
            }
            rVariance[1] = ( rSquaredTotal - rTotal * rTotal / counter ) / counter;
            rBar[1] = rTotal / counter;
            
            //is rVariance[1] > rVariance[0] ? then keep original coordinate
            if( rVariance[1] > rVariance[0] )
            {
                y0[1] = y0[0];
                rBar[1] = rBar[0];
                rVariance[1] = rVariance[0];
                noChangeY = 1;
            }
        }
        
        rBar[0] = rBar[1];
        rVariance[0] = rVariance[1];
        y0[0] = y0[1];
        
    }
    
    fParGeo->muonX0 = x0[1];
    fParGeo->muonY0 = y0[1];
    fParGeo->muonRadius = rBar[1];
    fParGeo->muonRSigma = sqrt( rVariance[1] );
    
    //cout <<"x0 "<<fParGeo->muonX0<<",y0 "<<fParGeo->muonY0<<",R "<<fParGeo->muonRadius<<",RS "<<fParGeo->muonRSigma<<" steps: "<<safty<<endl;
    
}

/*****************************************************************************
sizeInMuonRing

input: pointer to  origion of masking ring, radius and sigmaR
output: sum of all pixel values in the ring +-2 sigmaR
notes: this needs to correct for the fraction of tubes off in the ring region
	   this is rather tricky because the intensity depends on the ring diamter and the blurring!
*****************************************************************************/
void VImageParameterCalculation::sizeInMuonRing()
{
    if( fDebug )
    {
        cout << "VImageParameterCalculation::sizeInMuonRing" << endl;
    }
    if( !fData )
    {
        return;
    }
    
    //if( fParGeo->muonValid == 0 && fParGeo->houghMuonValid == 0 )
    //{
    //	fParGeo->muonSize = 0.0;
    //	fParGeo->muonIPCorrectedSize = 0.0;
    //	return;
    //}
    
    if( fParGeo->muonValid == 0 && fParGeo->houghNpix == 0 )//Calculates size for muonValid or if HT was run.
    {
        fParGeo->muonSize = 0.0;
        fParGeo->muonIPCorrectedSize = 0.0;
        return;
    }
    
    if( fData->getSums().size() == 0 || fParGeo->muonRadius == 0.0 )
    {
        fParGeo->muonSize = 0.0;
        fParGeo->muonIPCorrectedSize = 0.0;
        return;
    }
    if( !getDetectorGeometry() )
    {
        cout << "VImageParameterCalculation::sizeInMuonRing error: detector geometry not defined" << endl;
        exit( 0 );
    }
    
    unsigned int i;
    double xi, yi, size = 0.0, x0, y0, radius, si;
    int totalPixels = 0;
    int offPixels = 0;
    double xc = 0.; //Centroid x coordinate
    double yc = 0.; //Centroid y coordinate
    x0 = fParGeo->muonX0;
    y0 = fParGeo->muonY0;
    radius = fParGeo->muonRadius;
    
    for( i = 0; i < fData->getSums().size(); i++ )
    {
        xi = getDetectorGeometry()->getX()[i];
        yi = getDetectorGeometry()->getY()[i];
        double rp = sqrt( pow( xi - x0 , 2 ) + pow( yi - y0, 2 ) );
        
        if( rp > radius - 0.15 && rp < radius + 0.15 )
        {
            if( fData->getBorder()[i] || fData->getImage()[i] )
            {
                si = ( double )fData->getSums()[i]; // charge (dc)
                size += si;
                
                //Calculating the centroid here
                xc += xi * si;
                yc += yi * si;
                
                totalPixels ++;
            }
            if( fData->getDead()[i] )
            {
                offPixels ++;    // list of pixels turned off
            }
        }
    }
    // Correct for fraction of tubes turned off:
    if( totalPixels > 0 )
    {
        size *= ( 1.*offPixels ) / ( 1.0 * totalPixels ) + 1. ;
    }
    
    //Finish calculating the centroid by dividing by size and save.
    fParGeo->muonXC = xc / size;
    fParGeo->muonYC = yc / size;
    
    //Fill size
    fParGeo->muonSize = size;
    
    //Fill impact parameter corrected size
    
    //Impact parameter corrected size
    float correctedSize = 0.0;
    
    //If GSL is installed, run the code otherwise set correctedSize to 0.0
#ifndef NOGSL
    correctedSize = size / this->correctSizeInMuonRing();
#else
    correctedSize = 0.0;
#endif
    
    //Fill the tree
    fParGeo->muonIPCorrectedSize = correctedSize;
    
}

/*****************************************************************************
calcPixelDistribution
input: pointer to binaryPicture, origin of masking ring, radius and sigmaR
output: return 1 if there are at least 1 pixels in each octant segment of the ring
        AND at least half the pixels fall within +-1 sigmaR,
        otherwise returns 0
*****************************************************************************/
void VImageParameterCalculation::muonPixelDistribution()
{
    if( fDebug )
    {
        cout << "VImageParameterCalculation::muonPixelDistribution" << endl;
    }
    if( !fData || fData->getSums().size() == 0 || fParGeo->muonRadius == 0.0 )
    {
        fParGeo->muonValid = 0;
        return;
    }
    if( !getDetectorGeometry() )
    {
        cout << "VImageParameterCalculation::muonPixelDistribution error: detector geometry not defined" << endl;
        exit( 0 );
    }
    
    int totalPixels = 0, inside = 0, pass = 0;
    unsigned int i;
    int q[8] = {0};
    double rp, phi, xi, yi, x0, y0, radius, rsigma;
    
    //Get previously determined muon parameters
    x0 = fParGeo->muonX0;
    y0 = fParGeo->muonY0;
    radius = fParGeo->muonRadius;
    rsigma = fParGeo->muonRSigma;
    
    for( i = 0; i < fData->getSums().size(); i++ )
        if( fData->getImage()[i] || fData->getBorder()[i] )
        {
            totalPixels++;
            xi = getDetectorGeometry()->getX()[i];
            yi = getDetectorGeometry()->getY()[i];
            rp = sqrt( pow( xi - x0 , 2 ) + pow( yi - y0, 2 ) );
            phi = 180.0 / TMath::Pi() * atan2( yi - y0 , xi - x0 ) + 180.;
            
            if( rp > radius - 1.5 * rsigma  && rp < radius + 1.5 * rsigma )
            {
                inside++;
                
                //	    cout << "i "<<i<<" phi "<<phi<<endl;
                
                if( phi > 0 && phi < 45 )
                {
                    q[0]++;
                }
                else if( phi > 45 && phi < 90 )
                {
                    q[1]++;
                }
                else if( phi > 90 && phi < 135 )
                {
                    q[2]++;
                }
                else if( phi > 135 && phi < 180 )
                {
                    q[3]++;
                }
                else if( phi > 180 && phi < 225 )
                {
                    q[4]++;
                }
                else if( phi > 225 && phi < 270 )
                {
                    q[5]++;
                }
                else if( phi > 270 && phi < 315 )
                {
                    q[6]++;
                }
                else if( phi > 315 && phi < 360 )
                {
                    q[7]++;
                }
                
            }
            
        }
        
    //cout <<""<<inside<<" total: "<<totalPixels<<" q0:"<<q[0]<<" q1:"<<q[1]<<" q2:"<<q[2]<<" q3:"<<q[3]<<" q4:"<<q[4]<<" q5:"<<q[5]<<" q6:"<<q[6]<<" q7:"<<q[7]<<endl;
    
    if( ( 1.0 * inside / totalPixels > 0.7 ) && q[0] > 1 && q[1] > 1 && q[2] > 1 && q[3] > 1 && q[4] > 1 && q[5] > 1 && q[6] > 1 && q[7] > 1 )
    {
        pass = 1;
    }
    
    fParGeo->muonValid = pass;
    
}

//Calculates impact parameter correction factor for size
float VImageParameterCalculation::correctSizeInMuonRing()
{

    if( fDebug )
    {
        cout << "VImageParameterCalculation::correctSizeInMuonRing" << endl;
    }
    
    //Calculate the distance from the center of the ring to the centroid of the ring.
    float dc = TMath::Sqrt( ( fParGeo->muonX0 - fParGeo->muonXC ) * ( fParGeo->muonX0 - fParGeo->muonXC )
                            + ( fParGeo->muonY0 - fParGeo->muonYC ) * ( fParGeo->muonY0 - fParGeo->muonYC ) );
                            
    float kRatio;
    kRatio = dc / fParGeo->muonRadius;
    
    float exi; //Imapct parameter correction factor
    exi = 0.0;
    
    float kTest = 0.0;
    kTest = 0.0;
    
    // Initialize the elliptic integral
    const int numSteps = 1000;
    float ngExi[numSteps];
    
    //Calculate the elliptic integral for various values of xi
    for( int i = 0; i < numSteps; i++ )
    {
    
        ngExi[i] = 0.0;
        
        
        //Calculate the elliptic integral. This requires GSL.
#ifndef NOGSL
        float xi_tmp = ( float )i / ( float )numSteps;
        ngExi[i] = ( 2.0 / TMath::Pi() ) * ( gsl_sf_ellint_E( TMath::Pi() / 2.0 , xi_tmp , 0 ) ); //Calculate elliptic integral
#endif
        
    }//End of elliptical integral for loop
    
    //Numerically solve the integral equation to calculate IP corrrection factor
    for( int i = 0; i < numSteps; i++ )
    {
    
        float xi_tmp = ( float )i / ( float )numSteps;
        float kTest_tmp = 2.0 * ngExi[i] * kRatio / xi_tmp;
        
        if( fabs( 1.0 - kTest_tmp ) < fabs( 1.0 - kTest ) )
        {
            exi = ngExi[i];
            kTest = kTest_tmp;
        } //End of if
        
    }//End of for loop
    
    return exi; //Return impact parameter correction factor
    
}

//Hough transform class instantiation
void VImageParameterCalculation::houghInitialization()
{

    fHoughTransform = new VHoughTransform( getDetectorGeometry() );
    
}



//Hough transform muon identification algorithm.
//Sets fParGeo->houghMuonValid to 1 if the Hough tansform based muon ID technique find a muon, sets fParGeo->houghMuonValid to 0 otherwise.
void VImageParameterCalculation::houghMuonPixelDistribution()
{

    //Show method if debug mode is run.
    if( fDebug )
    {
        cout << "VImageParameterCalculation::houghMuonPixelDistribution" << endl;
    }
    
    //Look for data
    if( fData->getSums().size() == 0 )
    {
    
        fParGeo->houghMuonValid = 0;
        return;
        
    }
    
    //Look for detector geometry
    if( !getDetectorGeometry() )
    {
    
        cout << "VImageParameterCalculation::houghMuonPixelDistribution error: detector geometry not defined" << endl;
        exit( 0 );
        
    }
    
    //Initial hit pixel cut.
    //Run the analysis if the number of image or border pixels exceeds 1.
    
    int iNpix = 0;
    
    // Loop over the pixels and calculate the hit pixel count.
    for( unsigned int i = 0; i < fData->getSums().size(); i++ )
    {
        //If the pixel is an image or border pixel.
        if( fData->getImage()[i] || fData->getBorder()[i] )
            //Increment hit pixel count
        {
            iNpix++;
        }
        
    }
    
    //Read pixel cut values from the Hough transform object
    int npixMaxVal = 0;
    int npixMinVal = 0;
    
    npixMaxVal = fHoughTransform->getNpixMax( fData->getTelID() );//Get the npixmax cut for the given telescope
    npixMinVal = fHoughTransform->getNpixMin( fData->getTelID() );//Get the npixmin cut for the given telescope
    
    //Initial pixel cut. Runthe Hough transform analysis if the event passes the pixel cut.
    if( iNpix >= npixMinVal && iNpix <= npixMaxVal )
    {
        fHoughTransform->analysis( fData, fParGeo );    //Pass pointers to the data and tree parameters to the Hough transform analysis method and perform analysis.
    }
    else
    {
    
        fParGeo->houghMuonValid = 0;
        fParGeo->houghAP = 0.;
        fParGeo->houghTD = 0.;
        fParGeo->houghNpix = 0;
        fParGeo->houghCN = 0.;
        fParGeo->houghContained = 0.;
        
    }
    
}


void VImageParameterCalculation::calcTriggerParameters( vector<bool> fTrigger )
{
    // Calculate the trigger-level centroids
    double sumx_trig = 0.;                        // MS
    double sumy_trig = 0.;                        // MS
    double sumx2_trig = 0.;                       // MS
    double sumy2_trig = 0.;                       // MS
    int trig_tubes = 0;
    
    // calculate trigger-level parameters
    // loop over all pixels
    for( unsigned int j = 0; j < fTrigger.size(); j++ )
    {
        // select trigger pixel
        if( fTrigger[j] )
        {
            trig_tubes += 1;
            
            double xi = getDetectorGeometry()->getX()[j];
            double yi = getDetectorGeometry()->getY()[j];
            
            sumx_trig += xi;
            sumy_trig += yi;
            sumx2_trig += xi * xi;
            sumy2_trig += yi * yi;
        }
    }
    fParGeo->trig_tubes = trig_tubes;
    
    if( fParGeo->trig_tubes != 0 )
    {
        const double xmean_trig = sumx_trig / trig_tubes;
        const double ymean_trig = sumy_trig / trig_tubes;
        const double x2mean_trig = sumx2_trig / trig_tubes;
        const double y2mean_trig = sumy2_trig / trig_tubes;
        
        fParGeo->cen_x_trig = xmean_trig;
        fParGeo->cen_y_trig = ymean_trig;
        fParGeo->cen_x2_trig = x2mean_trig;
        fParGeo->cen_y2_trig = y2mean_trig;
    }
    
}


/*
 *   calculate image parameters
 *   (classic style)
 *
 *   see Fegan, D.J. J. Phys. G: Nucl. Part. Phys. 23 (1997) 1013-1060
 *
*/
void VImageParameterCalculation::calcParameters()
{
    if( !fData )
    {
        return;
    }
    
    if( fDebug )
    {
        cout << "VImageParameterCalculation::calcParameters: ";
        cout << "# of pixels: " << fData->getImageBorderNeighbour().size();
        cout << endl;
    }
    const double ZeroTolerence = 1e-8;
    
    if( fData->getSums().size() == 0 )
    {
        if( fDebug )
        {
            cout << "\t VImageParameterCalculation::calcParameters: fData->getSums().size() " << fData->getSums().size() << endl;
        }
        return;
    }
    
    if( !getDetectorGeometry() )
    {
        cout << "VImageParameterCalculation::calcParameters error: detector geometry not defined" << endl;
        exit( EXIT_FAILURE );
    }
    
    double sumsig = 0;
    double sumsig_2 = 0.;
    double sumxsig = 0;
    double sumysig = 0;
    double sumx2sig = 0;
    double sumy2sig = 0;
    double sumxysig = 0;
    double sumx3sig = 0;
    double sumy3sig = 0;
    double sumx2ysig = 0;
    double sumxy2sig = 0;
    int pntubes = 0;
    int pntubesBrightNoImage = 0;
    double sumOuterRing = 0.;                     // sum signal of image in outer ring
    double sumDeadRing = 0.;                      // sum signal of image ring around dead pixel
    double sumLowGain = 0.;
    
    // calculate mean ped and pedvar
    fParGeo->fmeanPed_Image = 0.;
    fParGeo->fmeanPedvar_Image = 0.;
    
    if( fData->hasFADCData() )
    {
        double nPixPed = 0.;
        for( unsigned int j = 0; j < fData->getImageBorderNeighbour().size(); j++ )
        {
            if( fData->getImageBorderNeighbour()[j] )
            {
                // mean ped and pedvar over image
                if( fData )
                {
                    if( j < fData->getPeds().size() )
                    {
                        fParGeo->fmeanPed_Image += fData->getPeds()[j];
                        fParGeo->fmeanPedvar_Image += fData->getPedvars( fData->getSumWindow() )[j];
                    }
                    nPixPed++;
                }
            }
        }
        if( nPixPed > 0 )
        {
            fParGeo->fmeanPed_Image /= nPixPed;
            fParGeo->fmeanPedvar_Image /= nPixPed;
        }
        else
        {
            fParGeo->fmeanPed_Image = 0.;         // -> pead = 0 means that pedvar is meanpedvar over camera
            if( fData )
            {
                fParGeo->fmeanPedvar_Image = fData->getmeanPedvars( false, fData->getSumWindow() );
            }
            else
            {
                fParGeo->fmeanPedvar_Image = 0.;
            }
        }
    }
    
    /////////////////////////////////////////
    // calculate image parameters
    /////////////////////////////////////////
    
    setImageBorderPixelPosition( fParGeo );
    
    // loop over all pixels
    for( unsigned int j = 0; j < fData->getSums().size(); j++ )
    {
        if( fData->getBrightNonImage()[j] )
        {
            pntubesBrightNoImage++;
        }
        
        // loop over image tubes
        // select image or border pixel
        if( fData->getImage()[j] || fData->getBorder()[j] )
        {
            pntubes += 1;
            
            double xi = getDetectorGeometry()->getX()[j];
            double yi = getDetectorGeometry()->getY()[j];
            
            double si = ( double )fData->getSums()[j]; // charge (dc)
            double si2 = ( double )fData->getSums2()[j];
            // image weighting with squared intensity
            // (non standard from traditional image calculation!)
            if( fData->getRunParameter() && fData->getRunParameter()->fSquaredImageCalculation )
            {
                si *= si;
                si2 *= si2;
            }
            sumsig += si;
            sumsig_2 += si2;
            // sum in outer ring
            if( getDetectorGeometry()->isEdgePixel()[j] )
            {
                sumOuterRing += si;
            }
            // sum around dead pixels
            if( j < getDetectorGeometry()->getNeighbours().size() || getDetectorGeometry()->getNNeighbours()[j] < getDetectorGeometry()->getMaxNeighbour() )
            {
                bool iDead = false;
                for( unsigned int n = 0; n < getDetectorGeometry()->getNeighbours()[j].size(); n++ )
                {
                    unsigned int k = getDetectorGeometry()->getNeighbours()[j][n];
                    if( k < fData->getDead().size() && fData->getDead( k, fData->getHiLo()[k] ) )
                    {
                        sumDeadRing += si;
                        iDead = true;
                        break;              // each pixel should only be added once
                    }
                }
                if( !iDead && getDetectorGeometry()->getNNeighbours()[j] < getDetectorGeometry()->getMaxNeighbour() )
                {
                    sumDeadRing += si;
                }
            }
            if( fData->getHiLo()[j] )
            {
                sumLowGain += si2;
            }
            
            const double sixi = si * xi;
            const double siyi = si * yi;
            
            sumxsig += sixi;
            sumysig += siyi;
            
            const double sixi2 = sixi * xi;
            const double siyi2 = siyi * yi;
            const double sixiyi = sixi * yi;
            
            sumx2sig += sixi2;
            sumy2sig += siyi2;
            sumxysig += sixiyi;
            
            sumx3sig += sixi2 * xi;
            sumy3sig += siyi2 * yi;
            sumx2ysig += sixi2 * yi;
            sumxy2sig += siyi2 * xi;
        }
    }
    if( fDebug )
    {
        cout << "VImageParameterCalculation::calcParameters: ";
        cout << "# of bright pixels: " << pntubes << endl;
    }
    
    fParGeo->ntubes = pntubes;
    fParGeo->ntubesBrightNoImage = pntubesBrightNoImage;
    if( sumsig > 0. )
    {
        fParGeo->loss = sumOuterRing / sumsig;
        fParGeo->lossAndDead = sumDeadRing / sumsig;
    }
    else
    {
        fParGeo->loss = 0.;
        fParGeo->lossAndDead = 0.;
    }
    if( sumLowGain > 0. && sumsig_2 > 0. )
    {
        fParGeo->fracLow = sumLowGain / sumsig_2;
    }
    else
    {
        fParGeo->fracLow = 0.;
    }
    
    //! clumsy, but effective, algorithm for finding 3 largest sums and pixel indices
    double i_max[3];
    i_max[0] = -1000.;
    i_max[1] = -1000.;
    i_max[2] = -1000.;
    int i_index[3];
    i_index[0] = 0;
    i_index[1] = 0;
    i_index[2] = 0;
    
    for( unsigned int i = 0; i < fData->getSums().size(); i++ )
    {
        if( fData->getSums()[i] > i_max[0] )
        {
            i_max[2] = i_max[1];
            i_max[1] = i_max[0];
            i_max[0] = fData->getSums()[i];
            i_index[2] = i_index[1];
            i_index[1] = i_index[0];
            i_index[0] = i;
        }
        else if( fData->getSums()[i] > i_max[1] )
        {
            i_max[2] = i_max[1];
            i_max[1] = fData->getSums()[i];
            i_index[2] = i_index[1];
            i_index[1] = i;
        }
        else if( fData->getSums()[i] > i_max[2] )
        {
            i_max[2] = fData->getSums()[i];
            i_index[2] = i;
        }
    }
    
    for( unsigned int i = 0; i < 3; i++ )
    {
        fParGeo->max[i] = i_max[i];
        if( i_max[0] > 0.0 )
        {
            fParGeo->frac[i] = i_max[i] / i_max[0];
        }
        fParGeo->index_of_max[i] = i_index[i];
    }
    
    if( fParGeo->ntubes == 0 )
    {
        fParGeo->reset( 1 );   // reset level 1: don't reset run numbers, event numbers, etc.
    }
    else
    {
        double xmean = 0.;
        double ymean = 0.;
        double x2mean = 0.;
        double y2mean = 0.;
        double xymean = 0.;
        if( sumsig > 0. )
        {
            xmean = sumxsig / sumsig;
            ymean = sumysig / sumsig;
            x2mean = sumx2sig / sumsig;
            y2mean = sumy2sig / sumsig;
            xymean = sumxysig / sumsig;
        }
        
        const double xmean2 = xmean * xmean;
        const double ymean2 = ymean * ymean;
        const double meanxy = xmean * ymean;
        
        double sdevx2 = x2mean - xmean2;
        double sdevy2 = y2mean - ymean2;
        const double sdevxy = xymean - meanxy;
        
        const double dist    = sqrt( xmean2 + ymean2 );
        
        fParGeo->size = sumsig;
        fParGeo->size2 = sumsig_2;
        fParGeo->sizeLL = -1.;
        fParGeo->size2LL = -1.;
        fParGeo->cen_x = xmean;
        fParGeo->cen_y = ymean;
        fParGeo->dist = dist;
        if( sdevx2 < ZeroTolerence )
        {
            sdevx2 = 0.;
        }
        fParGeo->sigmaX = sqrt( sdevx2 );
        if( sdevy2 < ZeroTolerence )
        {
            sdevy2 = 0.;
        }
        fParGeo->sigmaY = sqrt( sdevy2 );
        
        ////////////////////////////////////////////////////////////////////////////
        /////////////// directional cosines of the semi-major axis /////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        const double d = sdevy2 - sdevx2;
        const double z = sqrt( d * d + 4.0 * sdevxy * sdevxy );
        
        fParGeo->f_d = d;
        fParGeo->f_s = z;
        fParGeo->f_sdevxy = sdevxy;
        
        double cosphi = 0;
        double sinphi = 0;
        
        if( fabs( sdevxy ) > ZeroTolerence )      // length != width, semi-major axis
        {
            // not along x or y
            const double ac = ( d + z ) * ymean + 2.0 * sdevxy * xmean;
            const double bc = 2.0 * sdevxy * ymean - ( d - z ) * xmean;
            const double cc = sqrt( ac * ac + bc * bc );
            cosphi = bc / cc;
            sinphi = ac / cc;
        }
        else if( z > ZeroTolerence )              // semi-major axis along x or y
        {
            // and length != width
            cosphi = ( sdevx2 > sdevy2 ) ? 1 : 0;
            sinphi = ( sdevx2 > sdevy2 ) ? 0 : 1;
        }
        else if( dist > 0 )                       // length = width so might as well have semi-major axis
        {
            // be consistant with miss = dist, ie alpha = 90
            cosphi = -ymean / dist;
            // There seems to be a strange FP problem with the code below..
            //      sinphi= xmean / dist;
            sinphi = sqrt( 1.0 - cosphi * cosphi );
        }
        else
        {
            cosphi = 1;
            sinphi = 0;
        }
        
        fParGeo->cosphi = cosphi;
        fParGeo->sinphi = sinphi;
        
        ////////////////////////////////////////////////////////////////////////////
        //////////////// length, width and miss - image parameters /////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        double length2 = ( sdevx2 + sdevy2 + z ) / 2.0;
        if( length2 < ZeroTolerence )
        {
            length2 = 0;
        }
        
        double width2  = ( sdevx2 + sdevy2 - z ) / 2.0;
        if( width2 < ZeroTolerence )
        {
            width2 = 0;
        }
        
        const double length  = sqrt( length2 );
        const double width   = sqrt( width2 );
        
#ifdef OLDMISS
        double miss2 = 0;
        if( z > ZeroTolerence )
        {
            const double u = 1 + d / z;
            const double v = 2 - u;
            miss2 = ( u * xmean2 + v * ymean2 ) / 2.0 - meanxy * ( 2.0 * sdevxy / z );
            
            if( miss2 < ZeroTolerence )
            {
                //if ( miss2 < -(ZeroTolerence) )
                //throw Error("Miss squared is less than -ZeroTolerence");
                
                miss2 = 0;
            }
        }
        else
        {
            miss2 = xmean2 + ymean2;              // ie = dist ^ 2
        }
        const double miss = sqrt( miss2 );
#else
        double miss    = fabs( -sinphi * xmean + cosphi * ymean );
        if( miss > dist )
        {
            miss = dist;    // Weird rounding error
        }
#endif
        
        fParGeo->length = length;
        fParGeo->width = width;
        fParGeo->miss = miss;
        
        fParGeo->los = length / sumsig;
        
        ///////////////////////////////////////////////////////////////////////////////////
        ///////// fraction of image/border pixels located under image ellipse /////////////
        ///////////////////////////////////////////////////////////////////////////////////
        
        fParGeo->fui = getFractionOfImageBorderPixelUnderImage( xmean, ymean, width, length, cosphi, sinphi );
        
        ////////////////////////////////////////////////////////////////////////////
        /////////////////////// orientation: sinalpha and alpha ////////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        const double sinalpha = ( dist > ZeroTolerence ) ? miss / dist : 0;
        
        const double alpha = fabs( TMath::RadToDeg() * asin( sinalpha ) );
        
        fParGeo->alpha = alpha;
        
        ////////////////////////////////////////////////////////////////////////////
        //////////////////////////////// Azwidth ///////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        double azwidth;
        if( width2 > ZeroTolerence )
        {
            const double s2a = sinalpha * sinalpha;
            const double c2a = 1.0 - s2a;
            const double azfactor =
                1.0 + ( ( sinalpha == 0 ) ? 0.0 : ( length2 - width2 ) / ( width2 + length2 * c2a / s2a ) );
            azwidth = width * sqrt( azfactor );
        }
        else
        {
            azwidth = length;
        }
        fParGeo->azwidth = azwidth;
        
        ////////////////////////////////////////////////////////////////////////////
        //////////////////////// asymmetry major and minor /////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        double asymmetry = 0;
        if( length2 > ZeroTolerence )
        {
            const double x3mean = sumx3sig / sumsig;
            const double y3mean = sumy3sig / sumsig;
            const double x2ymean = sumx2ysig / sumsig;
            const double xy2mean = sumxy2sig / sumsig;
            
            const double sdevx3 = x3mean - 3.0 * xmean * x2mean + 2.0 * xmean * xmean2;
            const double sdevy3 = y3mean - 3.0 * ymean * y2mean + 2.0 * ymean * ymean2;
            const double sdevx2y =
                x2ymean - 2.0 * xymean * xmean + 2.0 * xmean2 * ymean - x2mean * ymean;
            const double sdevxy2 =
                xy2mean - 2.0 * xymean * ymean + 2.0 * xmean * ymean2 - xmean * y2mean;
                
            const double cosphi2 = cosphi * cosphi;
            const double sinphi2 = sinphi * sinphi;
            
            const double asymmetry3length3 =
                sdevx3 * cosphi * cosphi2 + 3.0 * sdevx2y * sinphi * cosphi2 +
                3.0 * sdevxy2 * cosphi * sinphi2 + sdevy3 * sinphi * sinphi2;
                
            if( fabs( asymmetry3length3 ) > ZeroTolerence )
            {
                asymmetry = pow( fabs( asymmetry3length3 ), 0.333333333333 ) / length;
                if( asymmetry3length3 < 0 )
                {
                    asymmetry = -asymmetry;
                }
            }
        }
        
        fParGeo->asymmetry = asymmetry;
        
        fParGeo->phi = atan2( fParGeo->sinphi, fParGeo->cosphi );
    }
    
    fboolCalcGeo = true;
    
}


/*!
    The loglikelihood fit is fitting a 2D Gaussian to the image/border pixels

*/
vector<bool> VImageParameterCalculation::calcLL( bool iUseSums2, bool i_reInitializeLL, bool iEqualSummationWindows )
{
    // take geometrical values as start values (calculate if not already calculated)
    if( !fboolCalcGeo )
    {
        calcParameters();
    }
    setImageBorderPixelPosition( fParLL );
    
    // fit image
    vector< bool > imageList = fImageFitter->calcLL( fParGeo,
                               fParLL,
                               iUseSums2,
                               i_reInitializeLL,
                               iEqualSummationWindows );
                               
    // calculate same reamining parameters
    fParLL->fui = getFractionOfImageBorderPixelUnderImage( fParLL->cen_x, fParLL->cen_y,
                  fParLL->width, fParLL->length,
                  fParLL->cosphi, fParLL->sinphi );
                  
    return imageList;
}



/*!
    calculate the fraction of image/border pixels which are located under the image ellipse.

    For a perfect image this is 1, for many cosmic rays this is < 1
*/
double VImageParameterCalculation::getFractionOfImageBorderPixelUnderImage( double cen_x, double cen_y,
        double width, double length,
        double cosphi, double sinphi )
{
    float i_ImageCoverFactor = 2.;
    if( fData && fData->getRunParameter() )
    {
        i_ImageCoverFactor = fData->getRunParameter()->fImageAnalysisFUIFactor;
    }
    if( i_ImageCoverFactor <= 0. )
    {
        cout << "VImageParameterCalculation::getFractionOfImageBorderPixelUnderImage error: fui factor <= 0, using 2" << endl;
        i_ImageCoverFactor = 2.;
    }
    float i_ImageCoverNPixel = 0;
    float i_ImageCoverNPixelImageBorder = 0;
    
    if( length > 1.e-2 && width > 1.e-2 )
    {
        // loop over all image and check if they are close to the image ellipse
        for( unsigned int i = 0; i < fData->getImage().size(); i++ )
        {
            // pixel coordinates rotated into frame of image ellipse
            double xi =     cosphi * ( getDetectorGeometry()->getX()[i] - cen_x ) + sinphi * ( getDetectorGeometry()->getY()[i] - cen_y );
            double yi = -1.*sinphi * ( getDetectorGeometry()->getX()[i] - cen_x ) + cosphi * ( getDetectorGeometry()->getY()[i] - cen_y );
            
            // check if these pixels are inside the image ellipse
            if( xi * xi / length / length / i_ImageCoverFactor / i_ImageCoverFactor + yi * yi / width / width / i_ImageCoverFactor / i_ImageCoverFactor < 1. )
            {
                i_ImageCoverNPixel++;
                if( fData->getImage()[i] || fData->getBorder()[i] )
                {
                    i_ImageCoverNPixelImageBorder++;
                }
            }
        }
        if( fDebug )
        {
            cout << "Fraction of pixel covered: " << i_ImageCoverNPixel << "\t" << i_ImageCoverNPixelImageBorder << "\t";
            if( i_ImageCoverNPixel > 0. )
            {
                cout << i_ImageCoverNPixelImageBorder / i_ImageCoverNPixel;
            }
            else
            {
                cout << 0.;
            }
            cout << endl;
        }
    }
    else
    {
        return 1.;
    }
    
    if( i_ImageCoverNPixel > 0. )
    {
        return i_ImageCoverNPixelImageBorder / i_ImageCoverNPixel;
    }
    
    return 0.;
}


void VImageParameterCalculation::setImageBorderPixelPosition( VImageParameter* iPar )
{
    if( fData && iPar )
    {
        vector< float > i_x;
        vector< float > i_y;
        for( unsigned int i = 0; i < fData->getImage().size(); i++ )
        {
            if( fData->getImage()[i] || fData->getBorder()[i] )
            {
                if( fData->getDetectorGeometry() && i < fData->getDetectorGeometry()->getX().size() && i < fData->getDetectorGeometry()->getY().size() )
                {
                    i_x.push_back( fData->getDetectorGeometry()->getX()[i] );
                    i_y.push_back( fData->getDetectorGeometry()->getY()[i] );
                }
            }
        }
        iPar->setImageBorderPixelPosition( i_x, i_y );
    }
}


/*
    fill image/border pixel to image parameter tree
    (optional)

    PixelType == 0: Pe > 0 and not image and not border pixel
    PixelType == 1: image pixel
    PixelType == 2: border pixel
    PixelType == 3: neighbour pixel to image/border

*/
void VImageParameterCalculation::fillImageBorderPixelTree()
{
    if( !fParGeo )
    {
        return;
    }
    if( !fParGeo->isWriteNImagePixels() )
    {
        return;
    }
    fParGeo->PixelListN = 0;
    for( unsigned int i = 0; i < fData->getSums().size(); i++ )
    {
        // decide of pixel should be added to list written to tree
        fParGeo->PixelID[fParGeo->PixelListN] = i;
        fParGeo->PixelType[fParGeo->PixelListN] = 99;
        if( fData->getImage()[i] )
        {
            fParGeo->PixelType[fParGeo->PixelListN] = 1;
        }
        else if( fData->getBorder()[i] )
        {
            fParGeo->PixelType[fParGeo->PixelListN] = 2;
        }
        else if( fData->getImageBorderNeighbour()[i] )
        {
            fParGeo->PixelType[fParGeo->PixelListN] = 3;
        }
        else if( i < fData->getPE().size() && fData->getPE()[i] > 0 )
        {
            fParGeo->PixelType[fParGeo->PixelListN] = 0;
        }
        if( fParGeo->PixelType[fParGeo->PixelListN] < 99 )
        {
            fParGeo->PixelIntensity[fParGeo->PixelListN] = fData->getSums()[i];
            fParGeo->PixelTimingT0[fParGeo->PixelListN] = fData->getPulseTime()[i];
            if( i < fData->getPE().size() && fData->getPE()[i] > 0 )
            {
                fParGeo->PixelPE[fParGeo->PixelListN] = fData->getPE()[i];
            }
            // low-gain channel: add '10' to pixel type
            if( fData->getHiLo()[i] )
            {
                fParGeo->PixelType[fParGeo->PixelListN] += 10;
            }
            fParGeo->PixelListN++;
        }
    }
}
