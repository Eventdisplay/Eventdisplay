/*! \class VModel3D
     \brief VModel3D class for 3D-reconstruction of showers
            based on Lemoine-Goumard et al. 2006
            adapted by J. Grube and G. Gyuk
*/

#include "VModel3D.h"

VModel3D::VModel3D()
{
    fData3D = 0;
    fModel3DParameters = 0;
    fEmissionHeightCalculator = 0;
    fModel3DFn = 0;
    fModelLnL = 0;
    fInitialized3D = false;
}

VModel3D::~VModel3D()
{
    delete fData3D;
    delete fModel3DFn;
}

void VModel3D::initialize()
{
    fData3D = new VModel3DData();
    fModel3DParameters = new VModel3DParameters();
    fEmissionHeightCalculator = new VEmissionHeightCalculator();
    fModel3DFn = new VModel3DFn();
    fModelLnL = new VModelLnL( fData3D, fModel3DFn );
    fData3D->fexNoise3D = 0.4; // HARD-WIRED
    fData3D->fDebug3D = false; // set true for debug mode
}

void VModel3D::doModel3D()
{
    //// initialized only at first call in the analysis run ////
    if( !fInitialized3D )
    {
        initOutput();
        initModel3DTree(); // initialize the output tree
        vector<unsigned int> iNpix3D;
        iNpix3D.resize( fData->getNTel(), 0 );
        for( unsigned int i = 0; i < fData->getNTel(); i++ )
        {
            iNpix3D[i] = fData->getDetectorGeo()->getX( i ).size();
        }
        fData3D->initModel3D( fData->getNTel(), iNpix3D ); // initialize data
        readLnLTable();    // read LnL lookup table
        setGain();         // get dc/pe ratio
        getDetector();     // get detector configuration
        fInitialized3D = true;
    }
    fData3D->initEventModel3D();  /// initialize event
    
    if( fData3D->fDebug3D )
    {
        cout << "--- event " << fData->getShowerParameters()->eventNumber << " ---" << endl;
    }
    
    /////// use reconstruction method 0 as default ////////////
    /////// or set with MODEL3DSTARTID in runparameter file ///
    unsigned int iID = fRunPar->fIDstartDirectionModel3D;
    if( fRunPar->fIDstartDirectionModel3D > getEvndispReconstructionParameter()->getNReconstructionCuts() )
    {
        iID = 0;
    }
    /// reconstruction quality check ///
    if( fData->getShowerParameters()->fShowerNumImages[iID] < 2 )
    {
        if( fData3D->fDebug3D )
        {
            cout << "doModel3D: NImages < 2" << endl;
        }
        writeParameters3D();
        return;
    }
    /// setup Model3D ///
    calcPointing();         // setup the pointing
    calcCoord();            // convert coordinates
    calcPvector();          // calc pixel line of sights
    calcStartParameters();  // calc Model3D start parameters
    if( ! fData3D->fGoodEvent3D )
    {
        if( fData3D->fDebug3D )
        {
            cout << "doModel3D: bad start parameters" << endl;
        }
        writeParameters3D();
        return;
    }
    /// run Model3D analysis ///
    setSelectedPixels();    // select pixels used in Model3D fit
    startModel3D();         // model with start values, proceed with fit or not
    if( ! fData3D->fGoodEvent3D )
    {
        if( fData3D->fDebug3D )
        {
            cout << "doModel3D: bad initial GOF" << endl;
        }
        writeParameters3D();
        return;
    }
    doFit();
    
    /// below for testing only ///
    //gridModel3D();
    //if( fData->getRunParameter()->fUseDisplayModel3D )
    //{
    //  calcParameters( 0 );
    //}
    //////////////////////////////
    
    writeParameters3D();
    return;
}

/////////////////////////////////////////////////////////////////////

void VModel3D::initOutput()
{
    // check if root outputfile exist
    if( fOutputfile != 0 )
    {
        return;
    }
    
    if( fRunPar->foutputfileName != "-1" )
    {
        printf( "Model3D: attempt to create file\n" );
        char i_textTitle[300];
        sprintf( i_textTitle, "VERSION %d", getRunParameter()->getEVNDISP_TREE_VERSION() );
        fOutputfile = new TFile( fRunPar->foutputfileName.c_str(), "RECREATE", i_textTitle );
    }
}

bool VModel3D::initModel3DTree()
{
    char i_text[300];
    char i_textTitle[300];
    sprintf( i_text, "model3Dpars" );
    // tree versioning numbers used in mscw_energy
    sprintf( i_textTitle, "Model3D Parameters (VERSION %d)\n", getRunParameter()->getEVNDISP_TREE_VERSION() );
    // make sure that tree is created in top level directory of output file
    if( getOutputFile() && getOutputFile()->cd() )
    {
        fModel3DParameters->initTree( i_text, i_textTitle );
    }
    else
    {
        cout << "VModel3D::initModel3DTree(): error; output file directory not found" << endl;
        return false;
    }
    return true;
}

void VModel3D::terminate()
{
    getModel3DParameters()->getTree()->Write();
}

void VModel3D::createLnLTable()
{
    double pixel_var = 1.0; //need to create tables for many noise levels?
    double pixel_nsb = 0.;
    fModelLnL->createLnLTable( pixel_var, pixel_nsb );
}

void VModel3D::readLnLTable()
{
    string LnLTableFile = fData->getRunParameter()->fLnLTableFile;
    fModelLnL->readLnLTable( LnLTableFile );
}

void VModel3D::setGain()
{
    /// set gain for each telescope ///
    for( unsigned int iTel = 0; iTel < fData->getNTel(); iTel++ )
    {
        fData3D->fDCPE[iTel] = fRunPar->fEpochGain[iTel];
        cout << "telescope absolute gains: " << fRunPar->fInstrumentEpoch << " " << iTel + 1 << " " << fData3D->fDCPE[iTel] << endl;
    }
}

void VModel3D::getDetector()
{
    ///// HARD-WIRED to VERITAS /////////////
    /// get telescope locations on ground ///
    for( unsigned int iTel = 0; iTel < fData->getNTel(); iTel++ )
    {
    
        fData3D->flocTel3D[0][iTel] = fData->getDetectorGeo()->getTelXpos()[iTel];
        fData3D->flocTel3D[1][iTel] = fData->getDetectorGeo()->getTelYpos()[iTel];
        fData3D->flocTel3D[2][iTel] = fData->getDetectorGeo()->getTelZpos()[iTel];
        
        // fMarea3D[iTel] = getDetectorGeo()->getMirrorArea()[iTel];
        fData3D->fMarea3D[iTel] = 111.; // FIX: Hard-wired VERITAS mirror area (m^2)
        //// use full pixel diameter, or effective diameter?? ////
        //double pwid = getDetectorGeo()->getTubeRadius( iTel )[j];
        double pwid = 0.148 * ( TMath::Pi() / 180. ); // Hard-wired VERITAS PMT diameter
        
        for( unsigned int j = 0; j < fData3D->fomegapix3D[iTel].size(); j++ )
        {
            fData3D->fomegapix3D[iTel][j] = TMath::Pi() * ( pwid / 2. ) * ( pwid / 2. );
        }
    }
    
    //// for shower max estimation from Hillas reconstruction ///
    const unsigned int ntel = fData->getNTel();
    double TelX[ntel];
    double TelY[ntel];
    double TelZ[ntel];
    for( unsigned int iTel = 0; iTel < ntel; iTel++ )
    {
        TelX[iTel] = fData3D->flocTel3D[0][iTel];
        TelY[iTel] = fData3D->flocTel3D[1][iTel];
        TelZ[iTel] = fData3D->flocTel3D[2][iTel];
    }
    fEmissionHeightCalculator->setTelescopePositions( ntel, TelX, TelY, TelZ );
    
}

void VModel3D::calcPointing()
{
    /// elevation and azimuth of telescopes (in ground coordinates) ///
    /// get mean pointing from valid telescopes ////
    
    for( unsigned int i = 0; i < fData->getNTel(); i++ )
    {
        fData3D->fTze3D[i] = fData->getShowerParameters()->fTelElevation[i];
        fData3D->fTaz3D[i] = fData->getShowerParameters()->fTelAzimuth[i];
    }
    
    unsigned int Npoint = 0; // number of telescopes with non-zero pointing
    double mTze = 0;  // mean zenith of telescopes (non-zero pointing)
    double mTaz = 0;  // mean azimuth of telescopes (non-zero pointing)
    
    for( unsigned int i = 0; i < fData->getNTel(); i++ )
    {
        if( fData3D->fTze3D[i] != 0 && fData3D->fTaz3D[i] != 0 )
        {
            Npoint += 1;
            mTze = mTze + fData3D->fTze3D[i];
            mTaz = mTaz + fData3D->fTaz3D[i];
        }
    }
    
    if( Npoint > 0 )
    {
        mTze = mTze / ( double )Npoint;
        mTaz = mTaz / ( double )Npoint;
    }
    
    mTze = ( 90. - mTze ) * ( TMath::Pi() / 180. ); // zenith angle in radians
    mTaz = mTaz * ( TMath::Pi() / 180. ); // az in radians
    
    //// from spherical to cartesian (ground coordinates) ////
    fData3D->fT3D[0] = sin( mTze ) * sin( mTaz );
    fData3D->fT3D[1] = sin( mTze ) * cos( mTaz );
    fData3D->fT3D[2] = cos( mTze );
}

void VModel3D::calcCoord()
{
    /// telescope parameters ///
    
    /// sky unit base vectors in ground coordinate frame (base) ///
    /////// z (Same as T ) //////////////
    fData3D->fzsg3D[0] = fData3D->fT3D[0];
    fData3D->fzsg3D[1] = fData3D->fT3D[1];
    fData3D->fzsg3D[2] = fData3D->fT3D[2];
    //// should already be normalized ////
    fData3D->norm3D( fData3D->fzsg3D[0], fData3D->fzsg3D[1], fData3D->fzsg3D[2] );
    
    /////// x /////////////
    fData3D->cross3D( fData3D->fzg3D[0], fData3D->fzg3D[1], fData3D->fzg3D[2], fData3D->fT3D[0], fData3D->fT3D[1], fData3D->fT3D[2], fData3D->fxsg3D[0], fData3D->fxsg3D[1], fData3D->fxsg3D[2] );
    fData3D->norm3D( fData3D->fxsg3D[0], fData3D->fxsg3D[1], fData3D->fxsg3D[2] );
    
    /////// y /////////////
    fData3D->cross3D( fData3D->fzsg3D[0], fData3D->fzsg3D[1], fData3D->fzsg3D[2], fData3D->fxsg3D[0], fData3D->fxsg3D[1], fData3D->fxsg3D[2], fData3D->fysg3D[0], fData3D->fysg3D[1], fData3D->fysg3D[2] );
    fData3D->norm3D( fData3D->fysg3D[0], fData3D->fysg3D[1], fData3D->fysg3D[2] );
    ////// swap Y sign /////
    fData3D->fysg3D[0] = 0. - fData3D->fysg3D[0];
    fData3D->fysg3D[1] = 0. - fData3D->fysg3D[1];
    fData3D->fysg3D[2] = 0. - fData3D->fysg3D[2];
    
}

void VModel3D::calcPvector()
{
    //// get p-hat vector (line of sight) for each telescope and pixel /////
    for( unsigned int iTel = 0; iTel < fData3D->fNTel3D; iTel++ )
    {
        for( unsigned int iPix = 0; iPix < fData3D->fNpix3D[iTel]; iPix++ )
        {
        
            //// add pixel "p" vector to "T" vector in ground coordinates ////
            double tp[3]; // temp p vector
            //// convert from deg to rad ////
            tp[0] = fData->getDetectorGeo()->getX( iTel )[iPix] * ( TMath::Pi() / 180. );
            tp[1] = fData->getDetectorGeo()->getY( iTel )[iPix] * ( TMath::Pi() / 180. );
            tp[2] = 0;
            //// multiply by base vectors to be in ground system ////
            double tpx[3];
            double tpy[3];
            double tpz[3];
            tpx[0] = tp[0] * fData3D->fxsg3D[0];
            tpx[1] = tp[0] * fData3D->fxsg3D[1];
            tpx[2] = tp[0] * fData3D->fxsg3D[2];
            tpy[0] = tp[1] * fData3D->fysg3D[0];
            tpy[1] = tp[1] * fData3D->fysg3D[1];
            tpy[2] = tp[1] * fData3D->fysg3D[2];
            tpz[0] = tp[2] * fData3D->fzsg3D[0];
            tpz[1] = tp[2] * fData3D->fzsg3D[1];
            tpz[2] = tp[2] * fData3D->fzsg3D[2];
            
            //// add to the telescope vector ////
            fData3D->fpX3D[iTel][iPix] = fData3D->fT3D[0] + tpx[0] + tpy[0] + tpz[0];
            fData3D->fpY3D[iTel][iPix] = fData3D->fT3D[1] + tpx[1] + tpy[1] + tpz[1];
            fData3D->fpZ3D[iTel][iPix] = fData3D->fT3D[2] + tpx[2] + tpy[2] + tpz[2];
            
            fData3D->norm3D( fData3D->fpX3D[iTel][iPix], fData3D->fpY3D[iTel][iPix], fData3D->fpZ3D[iTel][iPix] );
            
            //// angle between pixel center and telescope axis (in ground frame) ////
            fData3D->fcosptheta3D[iTel][iPix] = fData3D->dot3D( fData3D->fpX3D[iTel][iPix], fData3D->fpY3D[iTel][iPix], fData3D->fpZ3D[iTel][iPix], fData3D->fT3D[0], fData3D->fT3D[1], fData3D->fT3D[2] );
            
        }
    }
}

void VModel3D::calcStartParameters()
{
    /////// Get model starting point ///////////
    /////// use reconstruction method 0 as default ////////////
    /////// or set with MODEL3DSTARTID in runparameter file ///
    unsigned int iID = fRunPar->fIDstartDirectionModel3D;
    if( fRunPar->fIDstartDirectionModel3D > getEvndispReconstructionParameter()->getNReconstructionCuts() )
    {
        iID = 0;
    }
    fData3D->fStartXcore3D = fData->getShowerParameters()->fShowerXcore[iID];
    fData3D->fStartYcore3D = fData->getShowerParameters()->fShowerYcore[iID];
    
    double corecut = 9999;
    if( fabs( fData3D->fStartXcore3D ) >= corecut || fabs( fData3D->fStartYcore3D ) >= corecut )
    {
        fData3D->fGoodEvent3D = false;
    }
    
    //// shower direction (sky to ground coordinates) ////
    //// multiply by base vectors to be in ground system ////
    double ts[3]; // temp s vector
    ts[0] = fData->getShowerParameters()->fShower_Xoffset[iID];
    ts[1] = fData->getShowerParameters()->fShower_Yoffset[iID];
    ts[2] = 0;
    if( ts[0] == 0 && ts[1] == 0 )
    {
        fData3D->fGoodEvent3D = false;
    }
    ts[1] = 0. - ts[1];   // swap Y sign
    ts[0] = ts[0] * ( TMath::Pi() / 180. ); // in radians
    ts[1] = ts[1] * ( TMath::Pi() / 180. ); // in radians
    ts[2] = ts[2] * ( TMath::Pi() / 180. ); // in radians
    double tsx[3];
    double tsy[3];
    double tsz[3];
    tsx[0] = ts[0] * fData3D->fxsg3D[0];
    tsx[1] = ts[0] * fData3D->fxsg3D[1];
    tsx[2] = ts[0] * fData3D->fxsg3D[2];
    tsy[0] = ts[1] * fData3D->fysg3D[0];
    tsy[1] = ts[1] * fData3D->fysg3D[1];
    tsy[2] = ts[1] * fData3D->fysg3D[2];
    tsz[0] = ts[2] * fData3D->fzsg3D[0];
    tsz[1] = ts[2] * fData3D->fzsg3D[1];
    tsz[2] = ts[2] * fData3D->fzsg3D[2];
    //// add to the telescope vector ////
    ts[0] = fData3D->fT3D[0] + tsx[0] + tsy[0] + tsz[0];
    ts[1] = fData3D->fT3D[1] + tsx[1] + tsy[1] + tsz[1];
    ts[2] = fData3D->fT3D[2] + tsx[2] + tsy[2] + tsz[2];
    //// normalize s ////
    fData3D->norm3D( ts[0], ts[1], ts[2] );
    //// get shower spherical parameters (ground coordinates) ////
    fData3D->fStartSaz3D = atan2( ts[1], ts[0] ) * ( 180. / TMath::Pi() ); // in deg
    fData3D->fStartSel3D = 90. - ( acos( ts[2] ) * ( 180. / TMath::Pi() ) ); // in deg
    
    //// calculated shower max ////
    const unsigned int ntel = fData->getNTel();
    double cen_x[ntel];
    double cen_y[ntel];
    double size[ntel];
    double ArrayPointing_Elevation = fData->getShowerParameters()->fArrayPointing_Elevation;
    double ArrayPointing_Azimuth   = fData->getShowerParameters()->fArrayPointing_Azimuth;
    double width[ntel];  // for 3D width estimate
    
    for( unsigned int iTel = 0; iTel < ntel; iTel++ )
    {
        initializeDataReader(); // inherited from VEvndispData
        setTelID( iTel );       // inherited from VEvndispData
        cen_x[iTel] = fData->getImageParameters()->cen_x;
        cen_y[iTel] = fData->getImageParameters()->cen_y;
        size[iTel] = fData->getImageParameters()->size;
        width[iTel] = fData->getImageParameters()->width;
    }
    fEmissionHeightCalculator->getEmissionHeight( cen_x, cen_y, size, ArrayPointing_Azimuth, ArrayPointing_Elevation );
    fData3D->fStartSmax3D = fEmissionHeightCalculator->getMeanEmissionHeight();
    
    fData3D->fStartsigmaL3D = 3.;  // 3DLength: use fixed starting value of 3 km ///
    
    /// 3DWidth estimate ///
    double iWidthWeight = 0;
    double iWidthWeightTemp = 0;
    double iWidth = 0;
    for( unsigned int i = 0; i < ntel; i++ )
    {
        if( size[i] > 0 )
        {
            iWidthWeightTemp = log10( size[i] );
            iWidthWeight += iWidthWeightTemp;
            iWidth += iWidthWeightTemp * fData3D->fStartSmax3D * tan( width[i] / ( 45. / atan( 1. ) ) );
        }
    }
    if( iWidthWeight > 0. )
    {
        fData3D->fStartsigmaT3D = iWidth / iWidthWeight;
    }
    else
    {
        fData3D->fStartsigmaT3D = 0;
    }
    fData3D->fStartsigmaT3D *= 1000.; // from km to m
    
    //// log(Nc) estimate: use rough scaling
    double TotSize = 0;
    unsigned int TotTel = 0;
    for( unsigned int i = 0; i < fData3D->fNTel3D; i++ )
    {
        if( size[i] > 0 )
        {
            TotSize += size[i];
            TotTel += 1;
        }
    }
    if( TotSize > 0 )
    {
        TotSize = ( ( double )fData3D->fNTel3D / ( double )TotTel ) * TotSize;
    }
    fData3D->fStartNc3D = 370 * TotSize * 0.97179; //Nc
    fData3D->fStartNc3D = log( fData3D->fStartNc3D ); //log(Nc)
    
    //// sanity cuts on starting values ////
    if( fData3D->fStartSmax3D <= 0 || fData3D->fStartSmax3D > 100. || fData3D->fStartsigmaL3D <= 0 || fData3D->fStartsigmaL3D > 100 || fData3D->fStartsigmaT3D <= 0 || fData3D->fStartsigmaT3D > 100 || fData3D->fStartNc3D <= 0 || fData3D->fStartNc3D > 100. )
    {
        fData3D->fGoodEvent3D = false;
    }
    
    if( fData3D->fDebug3D )
    {
        cout << "Model3D Start: " << fData3D->fStartSel3D << " " << fData3D->fStartSaz3D << " " << fData3D->fStartXcore3D << " " << fData3D->fStartYcore3D << " " << fData3D->fStartSmax3D << " " << fData3D->fStartsigmaL3D << " " << fData3D->fStartsigmaT3D << " " << fData3D->fStartNc3D << endl;
    }
    
}

void VModel3D::setSelectedPixels()
{
    unsigned int totPix = 0;
    unsigned int telPix = 0;
    
    for( unsigned int iTel = 0; iTel < fData3D->fNTel3D; iTel++ )
    {
        initializeDataReader(); // inherited from VEvndispData
        setTelID( iTel );       // inherited from VEvndispData
        
        for( unsigned int iPix = 0; iPix < fData3D->fNpix3D[iTel]; iPix++ )
        {
            /// add two rings of non-dead pixels around cleaned images
            if( getImage()[iPix] || getBorder()[iPix] || getImageBorderNeighbour()[iPix] || getBorderBorderNeighbour()[iPix] )
            {
                fData3D->fClean3D[iTel][iPix] = true;
                telPix += 1;
                fData3D->fMeasuredSum3D[iTel][iPix] = getData()->getSums()[iPix] / fData3D->fDCPE[iTel];
                fData3D->fPedvar3D[iTel][iPix] = getData()->getPedvars( fData->getCurrentSumWindow()[iPix] )[iPix] / fData3D->fDCPE[iTel];
            }
        }
        fData3D->fNDFTel3D[iTel] = telPix;
        totPix += telPix;
        telPix = 0;
    }
    fData3D->fNDF3D = totPix;
    
    if( fData3D->fDebug3D )
    {
        cout << "Model3D: total selected pixels: " << fData3D->fNDF3D << endl;
    }
}

void VModel3D::startModel3D()
{
    double lnl = 0;
    vector<double> beta;
    vector< vector< double > > alpha;
    vector<double> vp( 8 ); // Model3D uses 8 parameters
    vp[0] = fData3D->fStartSel3D;
    vp[1] = fData3D->fStartSaz3D;
    vp[2] = fData3D->fStartXcore3D;
    vp[3] = fData3D->fStartYcore3D;
    vp[4] = fData3D->fStartSmax3D;
    vp[5] = fData3D->fStartsigmaL3D;
    vp[6] = fData3D->fStartsigmaT3D;
    vp[7] = fData3D->fStartNc3D;
    
    fModelLnL->val( vp, lnl, beta, alpha );
    fData3D->fStartGOF3D = fData3D->fGOF3D;
    if( fData3D->fDebug3D )
    {
        cout << "Model3D Start GOF: " << fData3D->fStartGOF3D << endl;
    }
    double maxGOF = 1000000000; // (removed) quality cut for GOF
    if( fData3D->fStartGOF3D > maxGOF )
    {
        fData3D->fGoodEvent3D = false;
    }
}

void VModel3D::fillMu( const vector<double>& vp )
{
    /// for display ///
    fModel3DFn->setModel3DFn( fData3D->flocTel3D, vp );
    const unsigned int nTel = fData3D->fNTel3D;
    vector<double> dmuda( 8 );
    for( unsigned int iTel = 0; iTel < nTel; iTel++ )
    {
        for( unsigned int iPix = 0; iPix < fData3D->fNpix3D[iTel]; iPix++ )
        {
            double mu = 0;
            fModel3DFn->calcModel3DFn( iTel, iPix, fData3D->fpX3D[iTel][iPix], fData3D->fpY3D[iTel][iPix], fData3D->fpZ3D[iTel][iPix], mu, dmuda );
            double qterm1 = fData3D->fMarea3D[iTel] * fData3D->fomegapix3D[iTel][iPix] * fData3D->fcosptheta3D[iTel][iPix];
            mu = mu * qterm1;
            fData3D->fMu3D[iTel][iPix] = mu;
            fData3D->fMuTel3D[iTel] += mu;
        }
    }
}

void VModel3D::fillMuDisplay()
{
    /// for display ///
    for( unsigned int iTel = 0; iTel < fData3D->fNTel3D; iTel++ )
    {
        initializeDataReader(); // inherited from VEvndispData
        setTelID( iTel ); // inherited from VEvndispData
        valarray<double> a0( fData3D->fNpix3D[iTel] );
        vector<bool> clean( fData3D->fNpix3D[iTel] );
        for( unsigned int iPix = 0; iPix < fData3D->fNpix3D[iTel]; iPix++ )
        {
            a0[iPix] = fData3D->fMu3D[iTel][iPix] * fData3D->fDCPE[iTel];
            clean[iPix] = fData3D->fClean3D[iTel][iPix];
        }
        fData->setModel3DMu( a0 );
        fData->setModel3DClean( clean );
    }
}

void VModel3D::doFit()
{
    //// likelihood minimization ////
    VMinimizerFactory* mf = VMinimizerFactory::getInstance();
    
    vector<double> beta;
    vector< vector< double > > alpha;
    
    vector<double> vp( 8 ); // Model3D uses 8 parameters
    vp[0] = fData3D->fStartSel3D;
    vp[1] = fData3D->fStartSaz3D;
    vp[2] = fData3D->fStartXcore3D;
    vp[3] = fData3D->fStartYcore3D;
    vp[4] = fData3D->fStartSmax3D;
    vp[5] = fData3D->fStartsigmaL3D;
    vp[6] = fData3D->fStartsigmaT3D;
    vp[7] = fData3D->fStartNc3D;
    
    VMinimizer* minimizer = mf->getMinimizer( *fModelLnL, vp );
    
    minimizer->setLoBound( 0, vp[0] - 10. );
    minimizer->setHiBound( 0, vp[0] + 10. );
    minimizer->setLoBound( 1, vp[1] - 10. );
    minimizer->setHiBound( 1, vp[1] + 10. );
    minimizer->setLoBound( 2, vp[2] - 100. );
    minimizer->setHiBound( 2, vp[2] + 100. );
    minimizer->setLoBound( 3, vp[3] - 100. );
    minimizer->setHiBound( 3, vp[3] + 100. );
    minimizer->setLoBound( 4, 0.5 ) ;
    minimizer->setHiBound( 4, 50. );
    minimizer->setLoBound( 5, 0.5 );
    minimizer->setHiBound( 5, 50. );
    minimizer->setLoBound( 6, 1. );
    minimizer->setHiBound( 6, 100. );
    minimizer->setLoBound( 7, 1. );
    minimizer->setHiBound( 7, 50. );
    
    ///test freezing parameter (not used) ///
    //minimizer->freeze( 0 );
    
    minimizer->minimize(); // minimize!!
    
    vector<double> param = minimizer->getParam(); // get best-fit values
    
    /// azimuth from -180 to 180 ///
    if( param[1] < -180. )
    {
        param[1] += 360.;
    }
    if( param[1] >  180. )
    {
        param[1] -= 360.;
    }
    
    vector<double> err = minimizer->getError();
    fData3D->fConverged3D = minimizer->getConverged();
    
    //// for final GOF ////////////
    double lnl = 0;
    fModelLnL->val( param, lnl, beta, alpha );
    
    fData3D->fSel3D    = param[0];
    fData3D->fSaz3D    = param[1];
    fData3D->fXcore3D  = param[2];
    fData3D->fYcore3D  = param[3];
    fData3D->fSmax3D   = param[4];
    fData3D->fsigmaL3D = param[5];
    fData3D->fsigmaT3D = param[6];
    fData3D->fNc3D     = param[7];
    
    fData3D->fErrorSel3D    = err[0];
    fData3D->fErrorSaz3D    = err[1];
    fData3D->fErrorXcore3D  = err[2];
    fData3D->fErrorYcore3D  = err[3];
    fData3D->fErrorSmax3D   = err[4];
    fData3D->fErrorsigmaL3D = err[5];
    fData3D->fErrorsigmaT3D = err[6];
    fData3D->fErrorNc3D     = err[7];
    
    //////////////////////
    /// slant depth
    /// Fundamentals of Neutrino Physics and Astrophysics
    /// by Giunti and Kim (page 396)
    double h0 = 6.4;       // scale height of atmosphere (km)
    double X0 = 1300; // atomsperic depth at sea level
    double Ds;         // slant depth of shower maximum
    double Szen, hmax;
    //// reduced width: atmospheric density vs altitude
    double Po = 101325; // sea level standard pressure (Pa)
    double To = 288.15; // sea level standard temperature (K)
    double g = 9.80665; // gravitational constant (m s^-2)
    double L = 6.5;     // temperature lapse rate (deg K/km)
    double R = 8.31447; // gas constant (J / mol*deg K)
    double M = 28.9644; // molecular weight of dry air (gm/mol)
    double UP = ( g * M ) / ( R * L );
    double T, P, rho;
    double Dwidth, rWidth, ErWidth;
    /// slant depth ///
    Szen = ( 90. - fData3D->fSel3D ) * ( TMath::Pi() / 180. );
    hmax = fData3D->fSmax3D * cos( Szen );
    Ds = X0 * exp( -hmax / h0 );
    /// reduced 3D Width ///
    T = To - ( L * hmax );
    P = Po * pow( 1. - ( L * hmax / To ), UP );
    rho = ( P * M ) / ( R * T * 1000 ); // kg m^-3
    Dwidth = fData3D->fsigmaT3D * 100. * rho * 0.001;
    ErWidth = fData3D->fErrorsigmaT3D * 100. * rho * 0.001;
    if( Ds > 0 )
    {
        rWidth = Dwidth / Ds;
        ErWidth = ErWidth / Ds;
    }
    else
    {
        rWidth = 0;
        ErWidth = 0;
    }
    fData3D->fDepth3D = Ds;
    fData3D->fRWidth3D = rWidth * 1000.;   // e-3
    fData3D->fErrRWidth3D = ErWidth * 1000.; // e-3
    
    ///////////////////////
    if( fData->getRunParameter()->fUseDisplayModel3D )   //for display
    {
        fillMu( param );
        fillMuDisplay();
    }
    if( fData3D->fDebug3D )
    {
        cout << "Model3D Fit GOF:   " << fData3D->fGOF3D << " PAR: " << fData3D->fSel3D << " " << fData3D->fSaz3D << " " << fData3D->fXcore3D << " " << fData3D->fYcore3D << " " << fData3D->fSmax3D << " " << fData3D->fsigmaL3D << " " << fData3D->fsigmaT3D << " " << fData3D->fNc3D << endl;
    }
    
    return;
}

void VModel3D::fillInit3D()
{
    if( !fInitialized3D )
    {
        vector<unsigned int> iNpix3D;
        iNpix3D.resize( fData->getNTel(), 0 );
        for( unsigned int i = 0; i < fData->getNTel(); i++ )
        {
            iNpix3D[i] = fData->getDetectorGeo()->getX( i ).size();
        }
        fData3D->initModel3D( fData->getNTel(), iNpix3D ); // initialize data
    }
    fData3D->initEventModel3D();
    writeParameters3D();
    //for display
    if( fData->getRunParameter()->fUseDisplayModel3D )
    {
        fillMuDisplay();
    }
}

void VModel3D::writeParameters3D()
{
    fData->getModel3DParameters()->eventNumber = fData->getShowerParameters()->eventNumber;
    /// write start parameters to file ///
    fData->getModel3DParameters()->fStartSel3D = fData3D->fStartSel3D;
    fData->getModel3DParameters()->fStartSaz3D = fData3D->fStartSaz3D;
    fData->getModel3DParameters()->fStartXcore3D = fData3D->fStartXcore3D;
    fData->getModel3DParameters()->fStartYcore3D = fData3D->fStartYcore3D;
    fData->getModel3DParameters()->fStartSmax3D = fData3D->fStartSmax3D;
    fData->getModel3DParameters()->fStartsigmaL3D = fData3D->fStartsigmaL3D;
    fData->getModel3DParameters()->fStartsigmaT3D = fData3D->fStartsigmaT3D;
    fData->getModel3DParameters()->fStartNc3D = fData3D->fStartNc3D;
    /// write best-fit parameters to file ///
    fData->getModel3DParameters()->fSel3D = fData3D->fSel3D;
    fData->getModel3DParameters()->fSaz3D = fData3D->fSaz3D;
    fData->getModel3DParameters()->fXcore3D = fData3D->fXcore3D;
    fData->getModel3DParameters()->fYcore3D = fData3D->fYcore3D;
    fData->getModel3DParameters()->fSmax3D = fData3D->fSmax3D;
    fData->getModel3DParameters()->fsigmaL3D = fData3D->fsigmaL3D;
    fData->getModel3DParameters()->fsigmaT3D = fData3D->fsigmaT3D;
    fData->getModel3DParameters()->fNc3D = fData3D->fNc3D;
    /// fit quality ///
    fData->getModel3DParameters()->fStartGoodness3D = fData3D->fStartGOF3D;
    fData->getModel3DParameters()->fGoodness3D = fData3D->fGOF3D;
    fData->getModel3DParameters()->fConverged3D = fData3D->fConverged3D;
    /// slant depth and reduced width ///
    fData->getModel3DParameters()->fDepth3D = fData3D->fDepth3D;
    fData->getModel3DParameters()->fRWidth3D = fData3D->fRWidth3D;
    fData->getModel3DParameters()->fErrRWidth3D = fData3D->fErrRWidth3D;
    /// write errors in fit parameters to file ///
    fData->getModel3DParameters()->fErrorSel3D = fData3D->fErrorSel3D;
    fData->getModel3DParameters()->fErrorSaz3D = fData3D->fErrorSaz3D;
    fData->getModel3DParameters()->fErrorXcore3D = fData3D->fErrorXcore3D;
    fData->getModel3DParameters()->fErrorYcore3D = fData3D->fErrorYcore3D;
    fData->getModel3DParameters()->fErrorSmax3D = fData3D->fErrorSmax3D;
    fData->getModel3DParameters()->fErrorsigmaL3D = fData3D->fErrorsigmaL3D;
    fData->getModel3DParameters()->fErrorsigmaT3D = fData3D->fErrorsigmaT3D;
    fData->getModel3DParameters()->fErrorNc3D = fData3D->fErrorNc3D;
    /// model direction ///
    calcModelDirection();
    calcOmega();
    ///////////////////////
    fData->getModel3DParameters()->fOmega3D = fData3D->fOmega3D;
    fData->getModel3DParameters()->fXoffModel3D = fData3D->fXoffModel3D;
    // test fYoffModel3D for NaN
    if( TMath::IsNaN( fData3D->fYoffModel3D ) )
    {
        fData->getModel3DParameters()->fConverged3D = 0;
        fData->getModel3DParameters()->fYoffModel3D = -9999.;
    }
    else
    {
        fData->getModel3DParameters()->fYoffModel3D = fData3D->fYoffModel3D;
    }
    fData->getModel3DParameters()->fXoffDeRot3D = fData3D->fXoffDeRot3D;
    fData->getModel3DParameters()->fYoffDeRot3D = fData3D->fYoffDeRot3D;
    /// write to file ///
    fData->getModel3DParameters()->getTree()->Fill();
}

void VModel3D::calcModelDirection()
{
    //// check for valid Model3D shower azimuth and elevation //
    if( fData3D->fSel3D < 1 && fData3D->fSaz3D < 1 )
    {
        return;
    }
    //// convert Model3D shower azimuth and elevation to sky direction //
    double tSzen = 0; // temp S zenith
    double tSaz = 0; // temp S azimuth
    tSzen = ( 90. - fData3D->fSel3D ) * ( TMath::Pi() / 180. ); // in radians
    tSaz = fData3D->fSaz3D * ( TMath::Pi() / 180. ); // in radians
    
    double tsi[3]; // temp input s vector
    
    //// from spherical to cartesian (ground coordinates) ////
    tsi[0] = sin( tSzen ) * cos( tSaz );
    tsi[1] = sin( tSzen ) * sin( tSaz );
    tsi[2] = cos( tSzen );
    
    double ts[3]; // temp output s vector
    
    ts[0] = ( -fData3D->fysg3D[1] * ( tsi[0] - fData3D->fT3D[0] ) + fData3D->fysg3D[0] * ( tsi[1] - fData3D->fT3D[1] ) ) / ( -( fData3D->fysg3D[1] * fData3D->fxsg3D[0] ) + ( fData3D->fysg3D[0] * fData3D->fxsg3D[1] ) );
    
    ts[1] = ( ( tsi[0] - fData3D->fT3D[0] ) - ( ts[0] * fData3D->fxsg3D[0] ) ) / fData3D->fysg3D[0];
    ts[1] = 0. - ts[1];   ////// (swap Y sign) ///////
    ts[2] = 0;
    
    ts[0] = ts[0] / ( TMath::Pi() / 180. ); // from radians to deg
    ts[1] = ts[1] / ( TMath::Pi() / 180. ); // from radians to deg
    ts[2] = ts[2] / ( TMath::Pi() / 180. ); // from radians to deg
    
    fData3D->fXoffModel3D = ts[0];
    fData3D->fYoffModel3D = ts[1];
    
    /// derotate ///
    if( getShowerParameters()->isMC() )
    {
        fData3D->fXoffDeRot3D = fData3D->fXoffModel3D;
        fData3D->fYoffDeRot3D = fData3D->fYoffModel3D;
    }
    else
    {
        double iUTC = 0.;
        double xrot = 0.;
        double yrot = 0.;
        if( getArrayPointing() )
        {
            iUTC = VSkyCoordinatesUtilities::getUTC( getShowerParameters()->MJD, getShowerParameters()->time );
            getArrayPointing()->derotateCoords( iUTC, fData3D->fXoffModel3D, -1.*fData3D->fYoffModel3D, xrot, yrot );
        }
        fData3D->fXoffDeRot3D = xrot;
        fData3D->fYoffDeRot3D = -1.*yrot;
    }
    return;
}

void VModel3D::calcOmega()
{
    //// check for valid Model3D shower azimuth and elevation //
    if( fData3D->fXoffDeRot3D < -999 || fData3D->fYoffDeRot3D < -999 )
    {
        return;
    }
    ///// get angular distance between Hillas and Model3D ////
    double x3D = fData3D->fXoffDeRot3D;
    double y3D = fData3D->fYoffDeRot3D;
    double z3D = 0;
    x3D = x3D * ( TMath::Pi() / 180. ); // in radians
    y3D = y3D * ( TMath::Pi() / 180. ); // in radians
    ///cout<<"3D: (unnorm X,Y) "<< x3D <<", " << y3D <<endl;
    
    fData3D->norm3D( x3D, y3D, z3D );
    unsigned int iID = fRunPar->fIDstartDirectionModel3D;
    if( fRunPar->fIDstartDirectionModel3D > getEvndispReconstructionParameter()->getNReconstructionCuts() )
    {
        iID = 0;
    }
    double xH = fData->getShowerParameters()->fShower_XoffsetDeRot[iID];
    double yH = fData->getShowerParameters()->fShower_YoffsetDeRot[iID];
    double zH = 0;
    xH = xH * ( TMath::Pi() / 180. ); // in radians
    yH = yH * ( TMath::Pi() / 180. ); // in radians
    ///cout<<"Hi: (unnorm X,Y) "<< xH <<", " << yH <<endl;
    
    fData3D->norm3D( xH, yH, zH );
    double omega_a = 0;
    omega_a = fData3D->dot3D( x3D, y3D, z3D, xH, yH, zH );
    omega_a = acos( omega_a );
    double Omega3D = 0;
    if( omega_a > 1e-10 )
    {
        Omega3D = log10( omega_a );
    }
    else
    {
        Omega3D = -9999;
    }
    fData3D->fOmega3D = Omega3D;
}

//// TEST: step around grid at minimum
void VModel3D::gridModel3D()
{
    double lnl = 0;
    vector<double> beta;
    vector< vector< double > > alpha;
    vector<double> vp( 8 ); // Model3D uses 8 parameters
    vp[0] = fData3D->fSel3D;
    vp[1] = fData3D->fSaz3D;
    vp[2] = fData3D->fXcore3D;
    vp[3] = fData3D->fYcore3D;
    vp[4] = fData3D->fSmax3D;
    vp[5] = fData3D->fsigmaL3D;
    vp[6] = fData3D->fsigmaT3D;
    vp[7] = fData3D->fNc3D;
    
    cout << "===== Model3D: GRID SEARCH =====" << endl;
    /// keep Sel,Saz,Xcore,Ycore fixed at minimum
    unsigned int Smax3D_n = 9;
    double Smax3D_step = 0.025;
    double Smax3D_start = vp[4] - Smax3D_step * ( ( double )Smax3D_n - 1. ) / 2.;
    ////
    unsigned int sigmaL3D_n = 9;
    double sigmaL3D_step = 0.025;
    double sigmaL3D_start = vp[5] - sigmaL3D_step * ( ( double )sigmaL3D_n - 1. ) / 2.;
    ////
    unsigned int sigmaT3D_n = 9;
    double sigmaT3D_step = 0.025;
    double sigmaT3D_start = vp[6] - sigmaT3D_step * ( ( double )sigmaT3D_n - 1. ) / 2.;
    ////
    unsigned int Nc3D_n = 9;
    double Nc3D_step = 0.025;
    double Nc3D_start = vp[7] - Nc3D_step * ( ( double )Nc3D_n - 1. ) / 2.;
    //// size of parameter space
    const unsigned int nTot = Smax3D_n * sigmaL3D_n * sigmaT3D_n * Nc3D_n;
    cout << "total parameter space: " << nTot << endl;
    double bestpar[8];
    for( unsigned int iPar = 0; iPar < vp.size(); iPar++ )
    {
        bestpar[iPar] = 0.;
    }
    double startgood = fData3D->fGOF3D;
    double bettergood = startgood;
    bool printall = false;   //set to print goodness for all steps
    ////
    vector< double > vSmax3D;
    vector< double > vsigmaL3D;
    vector< double > vsigmaT3D;
    vector< double > vNc3D;
    vSmax3D.reserve( nTot );
    vsigmaL3D.reserve( nTot );
    vsigmaT3D.reserve( nTot );
    vNc3D.reserve( nTot );
    
    /// epic loop over parameter range
    for( unsigned int j = 0; j < Smax3D_n; j++ )
    {
        double iSmax3D = Smax3D_start + Smax3D_step * j;
        ///
        for( unsigned int k = 0; k < sigmaL3D_n; k++ )
        {
            double isigmaL3D = sigmaL3D_start + sigmaL3D_step * k;
            ///
            for( unsigned int l = 0; l < sigmaT3D_n; l++ )
            {
                double isigmaT3D = sigmaT3D_start + sigmaT3D_step * l;
                ///
                for( unsigned int m = 0; m < Nc3D_n; m++ )
                {
                    double iNc3D = Nc3D_start + Nc3D_step * m;
                    /////
                    vSmax3D.push_back( iSmax3D );
                    vsigmaL3D.push_back( isigmaL3D );
                    vsigmaT3D.push_back( isigmaT3D );
                    vNc3D.push_back( iNc3D );
                }
            }
        }
    }
    /// epic loop over parameter range
    for( unsigned int iTot = 0; iTot < nTot ; iTot++ )
    {
        vp[4] = vSmax3D[iTot];
        vp[5] = vsigmaL3D[iTot];
        vp[6] = vsigmaT3D[iTot];
        vp[7] = vNc3D[iTot];
        
        fModelLnL->val( vp, lnl, beta, alpha );
        if( printall )
        {
            cout << fData3D->fGOF3D << " ";
        }
        if( fData3D->fGOF3D < bettergood )
        {
            bettergood = fData3D->fGOF3D;
            for( unsigned int iPar3 = 0; iPar3 < vp.size(); iPar3++ )
            {
                bestpar[iPar3] = vp[iPar3];
            }
        }
    }
    /// epic
    if( !printall )
    {
        cout << "BEST GOOD " << bettergood << " PAR ";
        for( unsigned int iPar3 = 0; iPar3 < vp.size(); iPar3++ )
        {
            cout << bestpar[iPar3] << " ";
        }
        cout << endl;
    }
    
    if( printall )
    {
        cout << endl;
    }
    cout << "===== Model3D: DONE WITH GRID SEARCH =====" << endl;
}

/// TEST: Hillas parameterization of the best-fit Model3D images ////
void VModel3D::calcParameters( int iTel )
{
    const double ZeroTolerence = 1e-8;
    initializeDataReader(); // inherited from VEvndispData
    setTelID( iTel );       // inherited from VEvndispData
    /// FIX for Model3D
    if( fData3D->fMu3D[iTel].size() == 0 )
    {
        if( fData3D->fDebug3D )
        {
            cout << "Model3D Tel" << iTel + 1 << " fMu3D: " << fData->getSums().size() << endl;
        }
        return;
    }
    double sumsig = 0;
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
    
    /////////////////////////////////////////
    // calculate image parameters
    /////////////////////////////////////////
    
    // loop over all pixels
    for( unsigned int j = 0; j < fData3D->fMu3D[iTel].size(); j++ )
    {
        if( fData->getBrightNonImage()[j] )
        {
            pntubesBrightNoImage++;
        }
        
        // select image or border pixel
        if( fData->getImage()[j] || fData->getBorder()[j] )
        {
            pntubes += 1;
            
            // loop over image tubes
            double xi = fData->getDetectorGeo()->getX()[j];
            double yi = fData->getDetectorGeo()->getY()[j];
            
            const double si = ( double )fData3D->fMu3D[iTel][j] * fData3D->fDCPE[iTel]; // charge (dc)
            sumsig += si;
            // sum in outer ring
            if( fData->getDetectorGeo()->getNNeighbours()[j] < fData->getDetectorGeo()->getMaxNeighbour() )
            {
                sumOuterRing += si;
            }
            // sum around dead pixels
            if( j < fData->getDetectorGeo()->getNeighbours().size() || fData->getDetectorGeo()->getNNeighbours()[j] < fData->getDetectorGeo()->getMaxNeighbour() )
            {
                bool iDead = false;
                for( unsigned int n = 0; n < fData->getDetectorGeo()->getNeighbours()[j].size(); n++ )
                {
                    unsigned int k = fData->getDetectorGeo()->getNeighbours()[j][n];
                    if( k < fData->getDead().size() && fData->getDead( k, fData->getHiLo()[k] ) )
                    {
                        sumDeadRing += si;
                        iDead = true;
                        break;              // each pixel should only be added once
                    }
                }
                if( !iDead && getDetectorGeo()->getNNeighbours()[j] < getDetectorGeo()->getMaxNeighbour() )
                {
                    sumDeadRing += si;
                }
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
    
    //JG//fParGeo->ntubes = pntubes;
    //JG//fParGeo->ntubesBrightNoImage = pntubesBrightNoImage;
    //if( sumsig > 0. )
    //{
    //  //JG//fParGeo->loss = sumOuterRing / sumsig;
    //  //JG//fParGeo->lossAndDead = sumDeadRing / sumsig;
    //}
    //else
    //{
    //  //JG//fParGeo->loss = 0.;
    //  //JG//fParGeo->lossAndDead = 0.;
    //}
    //else
    //{
    //  //JG//fParGeo->fracLow = 0.;
    //}
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
    
    //for( unsigned int i = 0; i < 3; i++ )
    //{
    //  //JG//fParGeo->max[i] = i_max[i];
    //  if( i_max[0] > 0.0 )
    //  {
    //    //JG//fParGeo->frac[i] = i_max[i] / i_max[0];
    //  }
    //  //JG//fParGeo->index_of_max[i] = i_index[i];
    //}
    
    if( pntubes != 0 )
    {
        double xmean = 0.;
        if( sumsig > 0. )
        {
            xmean = sumxsig / sumsig;
        }
        double ymean = 0.;
        if( sumsig > 0. )
        {
            ymean = sumysig / sumsig;
        }
        double x2mean = 0.;
        if( sumsig > 0. )
        {
            x2mean = sumx2sig / sumsig;
        }
        double y2mean = 0.;
        if( sumsig > 0. )
        {
            y2mean = sumy2sig / sumsig;
        }
        double xymean = 0.;
        if( sumsig > 0. )
        {
            xymean = sumxysig / sumsig;
        }
        const double xmean2 = xmean * xmean;
        const double ymean2 = ymean * ymean;
        const double meanxy = xmean * ymean;
        
        double sdevx2 = x2mean - xmean2;
        double sdevy2 = y2mean - ymean2;
        const double sdevxy = xymean - meanxy;
        
        //fParGeo->size = sumsig;
        //fParGeo->size2 = sumsig_2;
        //fParGeo->sizeLL = -1.;
        //fParGeo->size2LL = -1.;
        //fParGeo->cen_x = xmean;
        //fParGeo->cen_y = ymean;
        //fParGeo->dist = dist;
        //if( sdevx2 < ZeroTolerence )
        //  {
        //	sdevx2 = 0.;
        //  }
        //fParGeo->sigmaX = sqrt( sdevx2 );
        //if( sdevy2 < ZeroTolerence )
        //  {
        //	sdevy2 = 0.;
        //  }
        //fParGeo->sigmaY = sqrt( sdevy2 );
        
        ////////////////////////////////////////////////////////////////////////////
        /////////////// directional cosines of the semi-major axis /////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        const double d = sdevy2 - sdevx2;
        const double z = sqrt( d * d + 4.0 * sdevxy * sdevxy );
        
        ////////////////////////////////////////////////////////////////////////////
        //////////////// length, width and miss - image parameters /////////////////
        ////////////////////////////////////////////////////////////////////////////
        
        double length2 = ( sdevx2 + sdevy2 + z ) / 2.0;
        if( length2 < ZeroTolerence )
        {
            //if ( length2 < -(ZeroTolerence) )
            //	throw Error("Length squared is less than -ZeroTolerence");
            length2 = 0;
        }
        
        double width2  = ( sdevx2 + sdevy2 - z ) / 2.0;
        if( width2 < ZeroTolerence )
        {
            //if ( width2 < -(ZeroTolerence) )
            //throw Error("Width squared is less than -ZeroTolerence");
            width2 = 0;
        }
        
        const double length  = sqrt( length2 );
        const double width   = sqrt( width2 );
        
        cout << "TEST: length = " << length << ", width = " << width << endl;
        
        //fParGeo->length = length;
        //fParGeo->width = width;
        //fParGeo->miss = miss;
        //fParGeo->los = length / sumsig;
        
    }
}
