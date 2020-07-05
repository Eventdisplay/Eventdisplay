/*********************************************
 *
 * EVENDISPLAY - FROGS INTERFACE Ver. 0.0
 *
 *********************************************/

/*! /class VFrogs

  \brief  Called from VEventLoop. Opens templates and does minimization.


*/

#include "VFrogs.h"
#include "frogs.h"

#define FROGSDEBUG 0

VFrogs::VFrogs()
{

    fFrogsParameters = new VFrogsParameters();
    fFrogsParameters->fNTel = getData()->fNTel;
    reset();
    
    frogsRecID		  = getRunParameter()->ffrogsRecID;
    templatelistname = getRunParameter()->ffrogstemplatelist ;
    fparamfile		  = getRunParameter()->ffrogsparameterfile ;
    processParamFile() ;
    frogs_seed_gsl_rng( ffrogsRandomSeed );
}

void VFrogs::processParamFile()
{
    cout << "VFrogs::processParamFile() !! " ;
    //istream
    // open file 'fparamfile'
    // loop over each line
    // 	if line begins with star, and
    // 	if line keyword matches *
    // 		save value to variable
    
    string iFROGS_PARAMETER = getRunParameter()->getDirectory_EVNDISPAnaData() + "/Frogs/" + fparamfile;
    ifstream is;
    is.open( iFROGS_PARAMETER.c_str(), ifstream::in );
    if( !is )
    {
        cerr << "\nError, could not open frogs parameter file:  " << iFROGS_PARAMETER << endl;
        cout << "\nexiting...." << endl;
        exit( EXIT_FAILURE ) ;
    }
    cout << iFROGS_PARAMETER << endl;
    
    string is_line ;
    string temp  ;
    string temp2 ;
    string tmpEpoch ;
    double tmpLowerThresh ;
    double tmpFirstParam  ;
    double tmpSecondParam ;
    double tmpDCtoPE ;
    double tmpPMTNoise ;
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            // print runparameter to stdout
            if( !(is_stream>>std::ws).eof() )
            {
                is_stream >> temp;
                if( temp == "DCTOPE" )
                {
                    is_stream >> tmpEpoch ;
                    is_stream >> tmpDCtoPE ;
                    if( tmpEpoch == "V4" )
                    {
                        frogsDCtoPE[4] = tmpDCtoPE ;
                    }
                    else if( tmpEpoch == "V5" )
                    {
                        frogsDCtoPE[5] = tmpDCtoPE ;
                    }
                    else if( tmpEpoch == "V6" )
                    {
                        frogsDCtoPE[6] = tmpDCtoPE ;
                    }
                    else
                    {
                        cerr << "Error, in file '" << fparamfile
                             << "'', keyword DCTOPE has a bad epoch (2nd column) value '"
                             << tmpEpoch << "', which VFrogs doesn't know how to handle (only supports V4, V5, and V6), so we're exiting!" << endl;
                    }
                }
                else if( temp == "MUCORRECT" )
                {
                    is_stream >> tmpEpoch ;
                    is_stream >> tmpLowerThresh ;
                    is_stream >> tmpFirstParam ; // TODO: needs better name!
                    is_stream >> tmpSecondParam ;
                    if( tmpEpoch == "V4" )
                    {
                        frogsLowerThresh[4] = tmpLowerThresh  ;
                        frogsFirstParam[4]  = tmpFirstParam   ;
                        frogsSecondParam[4] = tmpSecondParam  ;
                    }
                    else if( tmpEpoch == "V5" )
                    {
                        frogsLowerThresh[5] = tmpLowerThresh ;
                        frogsFirstParam[5]  = tmpFirstParam  ;
                        frogsSecondParam[5] = tmpSecondParam ;
                    }
                    else if( tmpEpoch == "V6" )
                    {
                        frogsLowerThresh[6] = tmpLowerThresh ;
                        frogsFirstParam[6]  = tmpFirstParam  ;
                        frogsSecondParam[6] = tmpSecondParam ;
                    }
                    else
                    {
                        cerr << "Error, in file '" << fparamfile << "'', keyword MUCORRECT has a bad epoch (2nd column) value '"
                             << tmpEpoch << "', which VFrogs doesn't know how to handle (only supports V4, V5, and V6), so we're exiting!" << endl;
                    }
                }
                else if( temp == "PMTNOISE" )
                {
                    is_stream >> tmpEpoch;
                    is_stream >> tmpPMTNoise;
                    if( tmpEpoch == "V4" )
                    {
                        frogsPMTNoise[4] = tmpPMTNoise;
                    }
                    else if( tmpEpoch == "V5" )
                    {
                        frogsPMTNoise[5] = tmpPMTNoise;
                    }
                    else if( tmpEpoch == "V6" )
                    {
                        frogsPMTNoise[6] = tmpPMTNoise;
                    }
                    else
                    {
                        cerr << "Error, in file '" << fparamfile << "'', keyword PMTNOISE has a bad epoch (2nd column) value '"
                             << tmpEpoch << "', which VFrogs doesn't know how to handle (only supports V4, V5, and V6), so we're exiting!" << endl;
                    }
                }
                else if( temp == "MINIMIZATION" )
                {
                    is_stream >> temp2;
                    if( temp2 == "ON" )
                    {
                        frogsMinimization = true;
                        is_stream >> frogsDeltaXS;
                        is_stream >> frogsDeltaYS;
                        is_stream >> frogsDeltaXP;
                        is_stream >> frogsDeltaYP;
                        is_stream >> frogsDeltaLog10e;
                        is_stream >> frogsDeltaLambda;
                    }
                    else if( temp2 == "OFF" )
                    {
                        frogsMinimization = false;
                        frogsDeltaXS      = 1E-15;
                        frogsDeltaYS      = 1E-15;
                        frogsDeltaXP      = 1E-15;
                        frogsDeltaYP      = 1E-15;
                        frogsDeltaLog10e  = 1E-15;
                        frogsDeltaLambda  = 1E-15;
                    }
                }
                else if( temp == "CHEATMODE" )
                {
                    frogsCheating = true;
                }
                else if( temp == "EXPORTDATA" )
                {
                    is_stream >> frogsNBEventCalib;
                }
                else if( temp == "RANDOMSEED" )
                {
                    is_stream >> ffrogsRandomSeed;
                }
                else if( temp == "INTERPOLATION" )
                {
                    is_stream >> frogsInterpOrder;
                }
            }
        }
    }
    // print out to check sanity
    for( int i = 0; i < VFROGSNEPOCH; i++ )
    {
        if( TMath::Abs( frogsLowerThresh[i] ) > 1.e-5
                && TMath::Abs( frogsFirstParam[i] ) > 1.e-5
                && TMath::Abs( frogsSecondParam[i] ) > 1.e-5
                && TMath::Abs( frogsDCtoPE[i] ) > 1.e-5 )
        {
            cout << "frogs params: epoch V" << i << " frogsLowerThresh[" << i << "]=" << frogsLowerThresh[i] << endl;
            cout << "frogs params: epoch V" << i << " frogsFirstParam[" << i << "] =" << frogsFirstParam[i] << endl;
            cout << "frogs params: epoch V" << i << " frogsSecondParam[" << i << "]=" << frogsSecondParam[i] << endl;
            cout << "frogs params: epoch V" << i << " frogsDCtoPE[" << i << "]=" << frogsDCtoPE[i] << endl;
        }
    }
    cout << "frogs params: minimization mode ( 1 = ON, 0 = OFF )    : " << frogsMinimization << endl;
    cout << "              stepsize (xs, ys, xp, yp, log10e, lambda): ( " << frogsDeltaXS << ", " << frogsDeltaYS << ", " << frogsDeltaXP << ", "
         << frogsDeltaYP << ", " << frogsDeltaLog10e << ", " << frogsDeltaLambda << " )" << endl;
    cout << "frogs params: interpolation mode ( 0 = no interpolation, 1 = linear, 2 = quadratic ): " << frogsInterpOrder << endl;
    cout << "frogs params: random seed for differential evolution" ;
    if( ffrogsRandomSeed > 0 )
    {
        cout << ": " << ffrogsRandomSeed << endl;
    }
    else
    {
        cout << " will be set to system time." << endl;
    }
}

//================================================================
//================================================================

void VFrogs::reset()
{
    frogsRecID       = 0;
    templatelistname = "";
    frogsEventID     = 0;
    frogsGSLConStat  = 0;
    frogsNB_iter     = 0;
    frogsNImages     = 0;
    frogsSelectedImages = 0;
    frogsXS          = 0.;
    frogsXSerr       = 0.;
    frogsYS          = 0.;
    frogsYSerr       = 0.;
    frogsXP          = 0.;
    frogsXPerr       = 0.;
    frogsYP          = 0.;
    frogsYPerr       = 0.;
    frogsXPGC        = 0.;
    frogsYPGC        = 0.;
    frogsEnergy      = 0.;
    frogsEnergyerr   = 0.;
    frogsLambda      = 0.;
    frogsLambdaerr   = 0.;
    frogsGoodnessImg = 0.;
    frogsNpixImg     = 0;
    frogsGoodnessBkg = 0.;
    frogsNpixBkg     = 0;
    frogsXPStart     = 0.;
    frogsYPStart     = 0.;
    frogsXPED        = 0.;
    frogsYPED        = 0.;
    frogsXSStart     = 0.;
    frogsYSStart     = 0.;
    frogsXS_derot = 0;
    frogsYS_derot = 0;
    frogsZe = 0;
    frogsAz = 0;
    fInitialized     = false;
    fStartEnergyLoop = 0;
    for( int i = 0; i < 500; i++ )
    {
        frogsTemplateMu0[i] = 0.;
        frogsTemplateMu1[i] = 0.;
        frogsTemplateMu2[i] = 0.;
        frogsTemplateMu3[i] = 0.;
    }
    for( int i = 0; i < 4; i++ )
    {
        frogsTelGoodnessImg[i] = 0.;
        frogsTelGoodnessBkg[i] = 0.;
    }
    for( int i = 0; i < VFROGSNEPOCH; i++ )
    {
        frogsLowerThresh[i] = 0.;
        frogsFirstParam[i]  = 0.;
        frogsSecondParam[i] = 0.;
        frogsDCtoPE[i]		  = 0.;
        frogsPMTNoise[i]	  = 0.;
    }
    frogsMinimization	= true;
    frogsDeltaXS		= 0.02;
    frogsDeltaYS		= 0.02;
    frogsDeltaXP		= 5.0;
    frogsDeltaYP		= 5.0;
    frogsDeltaLog10e	= 0.03;
    frogsDeltaLambda	= 0.2;
    frogsInterpOrder	= 2;
    frogsCheating		= false;
    frogsNBEventCalib = 0;
    ffrogsRandomSeed = 0;

    for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
    {
         frogsR[i] = 0.;
    }
}

//================================================================
//================================================================

VFrogs::~VFrogs()
{

}

//================================================================
//================================================================
void VFrogs::doFrogsStuff( int eventNumber, string fArrayEpoch )
{

    int i = 0;
    int j = 0;
    
    // only at first call in the analysis run: initialize data class, set trees
    if( !fInitialized )
    {
        initAnalysis();
        fInitialized = true;
        readTableFrogs();
        fStartEnergyLoop = 0;
    }
    
    // get energy from mscw analysis for the given event number
    double inEnergy = getFrogsStartEnergy( eventNumber );
    
    // Store data from the Eventdisplay analysis to a FROGS structure
    // In this example the function frogs_convert_from_ed takes
    // arguments corresponding the the Eventdisplay analysis. An equivalent
    // function should be written for any other analysis package to
    // make the FROGS analysis usable.
    if( inEnergy != FROGS_BAD_NUMBER )
    {
    
        struct frogs_imgtmplt_in d;
        d = frogs_convert_from_ed( eventNumber, inEnergy, fArrayEpoch );
        
        //Print out the data contained in the FROGS structure frogs_imgtmplt_in
        //This is useful when developing a frogs_convert_from_XXXX function
        
        if( FROGSDEBUG )
        {
            frogs_print_raw_event( d );
        }
        
        //Call the FROGS analysis
        struct frogs_imgtmplt_out output;
        //output = frogs_img_tmplt( &d );
        char templatelistnamecstr[FROGS_FILE_NAME_MAX_LENGTH] ;
        int maxchar = FROGS_FILE_NAME_MAX_LENGTH - 1 ;
        
        // 'formatbuff' is so we only put the first "FROGS_FILE_NAME_MAX_LENGTH-1" characters of the templatelistname string into the char array 'templatelistnamecstr'
        char formatbuff[20] ;
        sprintf( formatbuff, "%%.%ds", maxchar ) ;
        sprintf( templatelistnamecstr, formatbuff, templatelistname.c_str() ) ;
        
        output = frogs_img_tmplt( &d, templatelistnamecstr );
        
        frogsEventID     = output.event_id;
        frogsGSLConStat  = output.gsl_convergence_status;
        frogsNB_iter     = output.nb_iter;
        frogsNImages     = output.nb_images;
        frogsSelectedImages = output.selected_images;
        frogsXS          = output.cvrgpt.xs;
        frogsXSerr       = output.cvrgpterr.xs;
        frogsYS          = output.cvrgpt.ys;
        frogsYSerr       = output.cvrgpterr.ys;
        frogsXP          = output.cvrgpt.xp;
        frogsXPerr       = output.cvrgpterr.xp;
        frogsYP          = output.cvrgpt.yp;
        frogsYPerr       = output.cvrgpterr.yp;
        frogsEnergy      = output.cvrgpt.log10e;
        frogsEnergyerr   = output.cvrgpterr.log10e;
        frogsLambda      = output.cvrgpt.lambda;
        frogsLambdaerr   = output.cvrgpterr.lambda;
        frogsGoodnessImg = output.goodness_img;
        frogsNpixImg     = output.npix_img;
        frogsGoodnessBkg = output.goodness_bkg;
        frogsNpixBkg     = output.npix_bkg;
        
        //calculate core pos in ground coordinates, core distance, reconstructed Ze/Az, derotated coordinates.
        transformResults();
        
        for( j = 0; j < 4; j++ )
        {
            frogsTelGoodnessImg[j] = output.tel_goodnessImg[j];
            frogsTelGoodnessBkg[j] = output.tel_goodnessBkg[j];
        }

        // dc to pe
        // (TMPTMP)
        double i_temp_dcpe = 5.3;
        
        for( j = 0; j < 4; j++ )
        {
            for( i = 0; i < 499; i++ )
            {
                if( j == 0 )
                {
                    frogsTemplateMu0[i] = output.tmplt_tubes[j][i] * i_temp_dcpe;
                }
                if( j == 1 )
                {
                    frogsTemplateMu1[i] = output.tmplt_tubes[j][i] * i_temp_dcpe;
                }
                if( j == 2 )
                {
                    frogsTemplateMu2[i] = output.tmplt_tubes[j][i] * i_temp_dcpe;
                }
                if( j == 3 )
                {
                    frogsTemplateMu3[i] = output.tmplt_tubes[j][i] * i_temp_dcpe;
                }
            }
        }
        
        if( frogsRecID < 0 || frogsRecID >= ( int )getShowerParameters()->fNMethods )
        {
            cout << "VFrogs: error: invalid frogsRecID (should be in the range [0," << getShowerParameters()->fNMethods << "]" << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        frogsXPStart   = getShowerParameters()->fShowerXcore_SC[frogsRecID];
        frogsYPStart   = getShowerParameters()->fShowerYcore_SC[frogsRecID];
        frogsXPED      = getShowerParameters()->fShowerXcore[frogsRecID];
        frogsYPED      = getShowerParameters()->fShowerYcore[frogsRecID];
        frogsXSStart   = fData->getShowerParameters()->fShower_Xoffset[frogsRecID]; //TEMP GH
        //frogsXSStart = getShowerParameters()->fShower_Xoffset[frogsRecID];
        frogsYSStart   = -1.0 * fData->getShowerParameters()->fShower_Yoffset[frogsRecID];
        
        getFrogsParameters()->frogsEventID		= getFrogsEventID();
        getFrogsParameters()->frogsGSLConStat	= getFrogsGSLConStat();
        getFrogsParameters()->frogsNB_iter		= getFrogsNB_iter();
        getFrogsParameters()->frogsNImages		= getFrogsNImages();
        getFrogsParameters()->frogsSelectedImages	= getFrogsSelectedImages();
        getFrogsParameters()->frogsXS				= getFrogsXS();
        getFrogsParameters()->frogsXSerr			= getFrogsXSerr();
        getFrogsParameters()->frogsYS				= getFrogsYS();
        getFrogsParameters()->frogsYSerr			= getFrogsYSerr();
        getFrogsParameters()->frogsXP				= getFrogsXP();
        getFrogsParameters()->frogsXPerr			= getFrogsXPerr();
        getFrogsParameters()->frogsYP				= getFrogsYP();
        getFrogsParameters()->frogsXPGC			= getFrogsXPGC();
        getFrogsParameters()->frogsYPGC			= getFrogsYPGC();
        getFrogsParameters()->frogsYPerr			= getFrogsYPerr();
        getFrogsParameters()->frogsEnergy		= getFrogsEnergy();
        getFrogsParameters()->frogsEnergyerr	= getFrogsEnergyerr();
        getFrogsParameters()->frogsLambda		= getFrogsLambda();
        getFrogsParameters()->frogsLambdaerr	= getFrogsLambdaerr();
        getFrogsParameters()->frogsGoodnessImg = getFrogsGoodnessImg();
        getFrogsParameters()->frogsNpixImg		= getFrogsNpixImg();
        getFrogsParameters()->frogsGoodnessBkg = getFrogsGoodnessBkg();
        getFrogsParameters()->frogsNpixBkg		= getFrogsNpixBkg();
        
        getFrogsParameters()->frogsXPStart = getFrogsXPStart();
        getFrogsParameters()->frogsYPStart = getFrogsYPStart();
        getFrogsParameters()->frogsXPED	  = getFrogsXPED();
        getFrogsParameters()->frogsYPED	  = getFrogsYPED();
        getFrogsParameters()->frogsXSStart = getFrogsXSStart();
        getFrogsParameters()->frogsYSStart = getFrogsYSStart();
        
        getFrogsParameters()->frogsZe = frogsZe;
        getFrogsParameters()->frogsAz = frogsAz;
        getFrogsParameters()->frogsXS_derot = frogsXS_derot;
        getFrogsParameters()->frogsYS_derot = frogsYS_derot;
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            getFrogsParameters()->frogsR[i] = frogsR[i];
        }
        getFrogsParameters()->frogsTelGoodnessImg0 = getFrogsTelGoodnessImg( 0 );
        getFrogsParameters()->frogsTelGoodnessImg1 = getFrogsTelGoodnessImg( 1 );
        getFrogsParameters()->frogsTelGoodnessImg2 = getFrogsTelGoodnessImg( 2 );
        getFrogsParameters()->frogsTelGoodnessImg3 = getFrogsTelGoodnessImg( 3 );
        getFrogsParameters()->frogsTelGoodnessBkg0 = getFrogsTelGoodnessBkg( 0 );
        getFrogsParameters()->frogsTelGoodnessBkg1 = getFrogsTelGoodnessBkg( 1 );
        getFrogsParameters()->frogsTelGoodnessBkg2 = getFrogsTelGoodnessBkg( 2 );
        getFrogsParameters()->frogsTelGoodnessBkg3 = getFrogsTelGoodnessBkg( 3 );
        
        getFrogsParameters()->getTree()->Fill();
        

        /*
         *  templates - used for dipslay only
         *
         */
        valarray<double> a0( frogsTemplateMu0, 499 );
        valarray<double> a1( frogsTemplateMu1, 499 );
        valarray<double> a2( frogsTemplateMu2, 499 );
        valarray<double> a3( frogsTemplateMu3, 499 );
        for( j = 0; j < 4; j++ )
        {
            setTelID( j );
            if( j == 0 )
            {
                fData->setTemplateMu( a0 );
            }
            if( j == 1 )
            {
                fData->setTemplateMu( a1 );
            }
            if( j == 2 )
            {
                fData->setTemplateMu( a2 );
            }
            if( j == 3 )
            {
                fData->setTemplateMu( a3 );
            }
        }
        
        //Print out the results of the image template analysis
        //The values returned by frogs_img_tmplt really are to be stored in
        //structures used in the analysis package in which FROGS is being sed.
        if( FROGSDEBUG )
        {
            frogs_print_param_spc_point( output );
        }
        
    }
    //return FROGS_OK;
    return;
}
//================================================================
//================================================================
void VFrogs::initAnalysis()
{
    initOutput();
    initFrogsTree();
}
//================================================================
//================================================================
void VFrogs::initOutput()
{
    if( fDebug )
    {
        cout << "void VFrogs::initOuput()" << endl;
    }
    // check if root outputfile exist
    if( fOutputfile != 0 )
    {
        printf( "FROGSPUT: Frogs will write results to existing output file: %s \n", getOutputFile()->GetName() );
        return;
    }
    // otherwise create it
    if( fRunPar->foutputfileName != "-1" )
    {
        printf( "FROGSPUT: Frogs foutputfileName = -1 : Attempt to create file\n" );
        char i_textTitle[300];
        sprintf( i_textTitle, "VERSION %d", getRunParameter()->getEVNDISP_TREE_VERSION() );
        if( getRunParameter()->fShortTree )
        {
            sprintf( i_textTitle, "%s (short tree)", i_textTitle );
        }
        fOutputfile = new TFile( fRunPar->foutputfileName.c_str(), "RECREATE", i_textTitle );
    }
    
    if( FROGSDEBUG )
    {
        printf( "FROGSPUT Are we initOutput??\n" );
    }
    
}
//================================================================
//================================================================
void VFrogs::initFrogsTree()
{

    if( FROGSDEBUG )
    {
        printf( "FROGSPUT Check in initTree\n" );
    }
    
    char i_text[300];
    ostringstream iTreeTitle;
    sprintf( i_text, "frogspars" );
    // tree versioning numbers used in mscw_energy
    iTreeTitle << "FROGSPUT: Frogs Parameters (VERSION " << getRunParameter()->getEVNDISP_TREE_VERSION() << ")";
    if( getRunParameter()->fShortTree )
    {
        iTreeTitle << " (short tree)";
    }
    fFrogsParameters->initTree( i_text, iTreeTitle.str().c_str() );
    
}
//================================================================
//================================================================
// read energies from mscw file
// fill two vectors with event numbers and energies for each events
// (note that only these events are used in the analysis)
void VFrogs::readTableFrogs()
{

    int eventNumber = 0;
    double ErecS = 0.;
    
    string fmscwFrogsFile = getRunParameter()->ffrogsmscwfile;
    cout << "FROGS readTableFrogs for starting values from: " << getRunParameter()->ffrogsmscwfile << endl;
    
    TFile* mscwFrogsFile = new TFile( fmscwFrogsFile.c_str() , "READ" );
    
    if( mscwFrogsFile->IsZombie() )
    {
        cout << "VFrogs::readTableFrogs error: File " << fmscwFrogsFile << " does not exist!" << endl;
        exit( EXIT_FAILURE );
    }
    
    TTree* mscwTreeFrogs = ( TTree* )mscwFrogsFile->Get( "data" );
    mscwTreeFrogs->SetBranchAddress( "eventNumber", &eventNumber );
    mscwTreeFrogs->SetBranchAddress( "ErecS", &ErecS );
    
    int nentries = mscwTreeFrogs->GetEntries();
    fTableEventNumber.reserve( nentries );
    fTableEnergy.reserve( nentries );
    
    for( int i = 0 ; i < nentries; i++ )
    {
        mscwTreeFrogs->GetEntry( i );
        fTableEventNumber.push_back( eventNumber );
        fTableEnergy.push_back( ErecS );
    }
    
    mscwFrogsFile->Close();
    
    cout << "Finished Reading: " << fTableEventNumber.size() << "/" << fTableEnergy.size() << " events from mscw_energy file " << endl;
    
}
//================================================================
//================================================================
double VFrogs::getFrogsStartEnergy( int eventNumber )
{

    int mscwTotal = fTableEventNumber.size();

    for( int i = fStartEnergyLoop ; i < mscwTotal ; i++ )
    {
        if( eventNumber == fTableEventNumber[i] )
        {
            if( fStartEnergyLoop > 0 )
            {
                fStartEnergyLoop = i - 1;
            }
            return fTableEnergy[i];
        }
        else if( eventNumber < fTableEventNumber[i] )
        {
            return FROGS_BAD_NUMBER;
        }
    }
    return FROGS_BAD_NUMBER;
    
}
//================================================================
//================================================================
int VFrogs::getFrogsAnasumNumber( int eventNumber, int runNumber )
{

    int AnasumTotal = fAnasumRunNumber.size();
    
    for( int i = 0 ; i < AnasumTotal ; i++ )
    {
        if( runNumber == fAnasumRunNumber[i] && eventNumber == fAnasumEventNumber[i] )
        {
            return eventNumber;
        }
    }
    
    return FROGS_BAD_NUMBER;
    
}
//================================================================
//================================================================
float VFrogs::transformTelescopePosition( int iTel, float i_ze, float i_az, int axis )
{
    // transform telescope positions from ground into shower coordinates
    float i_xrot, i_yrot, i_zrot;
    float i_xcos = 0.;
    float i_ycos = 0.;
    
    // calculate direction cosine
    i_xcos = sin( i_ze / TMath::RadToDeg() ) * sin( ( i_az - 180. ) / TMath::RadToDeg() );
    i_ycos = sin( i_ze / TMath::RadToDeg() ) * cos( ( i_az - 180. ) / TMath::RadToDeg() );
    
    setTelID( iTel );
    // call to GrIsu routine
    tel_impact( i_xcos, i_ycos, getDetectorGeo()->getTelXpos()[iTel], getDetectorGeo()->getTelYpos()[iTel], getDetectorGeo()->getTelZpos()[iTel], &i_xrot, &i_yrot, &i_zrot, false );
    
    if( axis == 0 )
    {
        return i_xrot;
    }
    else if( axis == 1 )
    {
        return i_yrot;
    }
    else
    {
        return FROGS_BAD_NUMBER;
    }
}
//================================================================
//================================================================
float VFrogs::transformShowerPosition( float i_ze, float i_az, float xcore, float ycore, float zcore, int axis )
{
    // transform coordinates from shower coord system to ground coord
    // see also void VArrayAnalyzer::transformTelescopePosition in VArrayAnalyzer.cpp
    // and void VGrIsuAnalyzer::tel_impact in VGrIsuAnalyzer.cpp
    float i_xrot, i_yrot, i_zrot;
    float i_xcos = 0.;
    float i_ycos = 0.;
    
    // calculate direction cosine
    i_xcos = sin( i_ze / TMath::RadToDeg() ) * sin( ( i_az - 180. ) / TMath::RadToDeg() );
    i_ycos = sin( i_ze / TMath::RadToDeg() ) * cos( ( i_az - 180. ) / TMath::RadToDeg() );
    
    // call to GrIsu routine
    // the parameter true/false sets the transformation matrix
    // false: ground coordinate system -> shower coord system
    // true: the inverse matrix is ON. shower coord -> ground coordinate
    tel_impact( i_xcos, i_ycos, xcore, ycore, zcore, &i_xrot, &i_yrot, &i_zrot, true );
    
    if( axis == 0 )
    {
        return i_xrot;
    }
    else if( axis == 1 )
    {
        return i_yrot;
    }
    else
    {
        return FROGS_BAD_NUMBER;
    }
}
//================================================================
//================================================================
float VFrogs::transformPosition( float i_ze, float i_az, float x, float y, float z, int axis, bool bInv )
{
    // transform coordinates from shower coord system to ground coord
    // see also void VArrayAnalyzer::transformTelescopePosition in VArrayAnalyzer.cpp
    // and void VGrIsuAnalyzer::tel_impact in VGrIsuAnalyzer.cpp
    float i_xrot, i_yrot, i_zrot;
    float i_xcos = 0.;
    float i_ycos = 0.;
    
    // calculate direction cosine
    i_xcos = sin( i_ze / TMath::RadToDeg() ) * sin( ( i_az - 180. ) / TMath::RadToDeg() );
    i_ycos = sin( i_ze / TMath::RadToDeg() ) * cos( ( i_az - 180. ) / TMath::RadToDeg() );
    
    // call to GrIsu routine
    // the parameter true/false sets the transformation matrix
    // false: ground coordinate system -> shower coord system
    // true: the inverse matrix is ON. shower coord -> ground coordinate
    tel_impact( i_xcos, i_ycos, x, y, z, &i_xrot, &i_yrot, &i_zrot, bInv );
    
    if( axis == 0 )
    {
        return i_xrot;
    }
    else if( axis == 1 )
    {
        return i_yrot;
    }
    else if( axis == 2 )
    {
        return i_zrot;
    }
    else
    {
        return FROGS_BAD_NUMBER;
    }
}
//================================================================
//================================================================
void VFrogs::terminate()
{

    getFrogsParameters()->getTree()->Write();
    
}

void VFrogs::finishFrogs( TFile* f )
{

    //free random number generator
    frogs_free_gsl_rng();
    
    // Open outfile again to copy frogs tree to mscw file.
    
    // reopen mscw file
    string fmscwFrogsFile = getRunParameter()->ffrogsmscwfile;
    TFile* mscwFrogsFile = new TFile( fmscwFrogsFile.c_str() , "UPDATE" );
    if( mscwFrogsFile->IsZombie() )
    {
        cout << "Finish Frogs:" << endl;
        cout << "VFrogs::readTableFrogs error: File " << fmscwFrogsFile.c_str() << " does not exist!" << endl;
        exit( EXIT_FAILURE );
    }
    
    // Clone tree to mscw file checking it opened
    if( f->IsZombie() )
    {
        cout << "error: finish Frogs f file problem: "  << endl;
    }
    else
    {
        ( ( TTree* )f->Get( "frogspars" ) )->CloneTree()->Write();
    }
    
    // Close the files
    mscwFrogsFile->Close();
    
    return;
    
}
//================================================================
//================================================================
// convert starting parameters from eventdisplay analysis to
// frogs frame
struct frogs_imgtmplt_in VFrogs::frogs_convert_from_ed( int eventNumber, double inEnergy, string fArrayEpoch )
{
    /* The frogs_convert_from_grisu function is called with
       arguments containing all the data necessary to the image template
       analysis as per the structures used in grisu analysis. It can serve
       as an example for the interfacing of the template analysis with other
       analysis packages. It returns the data necessary to the template
       analysis in a structure that is appropriate.  */
    
    struct frogs_imgtmplt_in rtn;
    
    //Tracked elevation from telescope 0
    //rtn.elevation = fData->getShowerParameters()->fTelElevation[0];
    //rtn.azimuth   = fData->getShowerParameters()->fTelAzimuth[0];
    //rtn.event_id  = eventNumber;
    
    if( getArrayPointing() )
    {
        rtn.elevation = fData->getArrayPointing()->getTelElevation();
        rtn.azimuth   = fData->getArrayPointing()->getTelAzimuth();
    }
    else
    {
        rtn.elevation = 0.;
        rtn.azimuth   = 0.;
    }
    rtn.event_id				  = eventNumber;
    rtn.epoch_id				  = fArrayEpoch.c_str();
    string		 isub			  = fArrayEpoch.substr( 1, 1 );
    int			 epoch_as_int = atoi( isub.c_str() );
    rtn.lowerthresh			  = frogsLowerThresh[ epoch_as_int ];
    rtn.firstparam			  = frogsFirstParam [ epoch_as_int ];
    rtn.secondparam			  = frogsSecondParam[ epoch_as_int ];
    
    rtn.delta_xs	  = frogsDeltaXS;
    rtn.delta_ys	  = frogsDeltaYS;
    rtn.delta_xp	  = frogsDeltaXP;
    rtn.delta_yp	  = frogsDeltaXP;
    rtn.delta_log10e = frogsDeltaLog10e;
    rtn.delta_lambda = frogsDeltaLambda;
    
    rtn.interporder	  = frogsInterpOrder;
    rtn.nb_events_calib = frogsNBEventCalib;
    
    //Telescopes
    rtn.ntel = fData->getNTel(); //Number of telescopes
    
    rtn.scope				 = new struct frogs_telescope [rtn.ntel];
    rtn.nb_live_pix_total = 0;	//Total number or pixels in use
    for( int tel = 0; tel < rtn.ntel; tel++ )
    {
        initializeDataReader();
        setTelID( tel );
        
        //Telescope position in the shower coordinate coordinate system used in the reconstruction
        //rtn.scope[tel].xfield =
        //transformTelescopePosition( tel, 90. - fData->getShowerParameters()->fTelElevation[0], fData->getShowerParameters()->fTelAzimuth[0], 0 );
        //rtn.scope[tel].yfield =
        //transformTelescopePosition( tel, 90. - fData->getShowerParameters()->fTelElevation[0], fData->getShowerParameters()->fTelAzimuth[0], 1 );
        rtn.scope[tel].xfield =
            transformPosition( 90. - fData->getShowerParameters()->fArrayPointing_Elevation, fData->getShowerParameters()->fArrayPointing_Azimuth, getDetectorGeo()->getTelXpos()[tel], getDetectorGeo()->getTelYpos()[tel], getDetectorGeo()->getTelZpos()[tel], 0, false );
        rtn.scope[tel].yfield =
            transformPosition( 90. - fData->getShowerParameters()->fArrayPointing_Elevation, fData->getShowerParameters()->fArrayPointing_Azimuth, getDetectorGeo()->getTelXpos()[tel], getDetectorGeo()->getTelYpos()[tel], getDetectorGeo()->getTelZpos()[tel], 1, false );
        rtn.scope[tel].zfield =
            transformPosition( 90. - fData->getShowerParameters()->fArrayPointing_Elevation, fData->getShowerParameters()->fArrayPointing_Azimuth, getDetectorGeo()->getTelXpos()[tel], getDetectorGeo()->getTelYpos()[tel], getDetectorGeo()->getTelZpos()[tel], 2, false );
            
            
        if( FROGSDEBUG )
            printf( "TelSC %d | %.2f %.2f | %.2f %.2f | %.2f %.2f | %.2f %.2f | %.2f %.2f \n", tel,
                    fData->getShowerParameters()->fShowerZe[frogsRecID],
                    fData->getShowerParameters()->fShowerAz[frogsRecID],
                    getDetectorGeo()->getTelXpos()[tel] - fData->getShowerParameters()->fShowerXcore[frogsRecID],
                    getDetectorGeo()->getTelYpos()[tel] - fData->getShowerParameters()->fShowerYcore[frogsRecID],
                    getDetectorGeo()->getTelXpos()[tel],
                    getDetectorGeo()->getTelYpos()[tel],
                    rtn.scope[tel].xfield,
                    rtn.scope[tel].yfield,
                    ( 180.0 / 3.14159265 )*atan2( rtn.scope[tel].yfield, rtn.scope[tel].xfield ),
                    ( 180.0 / 3.14159265 )*atan2( getDetectorGeo()->getTelYpos()[tel], getDetectorGeo()->getTelXpos()[tel] ) );
                    
        //Telescope effective collection area
        //Number of pixels
        rtn.scope[tel].npix = fData->getDetectorGeo()->getNChannels()[tel];
        
        //Set the dimension of the pixel parameter arrays
        rtn.scope[tel].xcam		  = new float [rtn.scope[tel].npix];
        rtn.scope[tel].ycam		  = new float [rtn.scope[tel].npix];
        rtn.scope[tel].q		  = new float [rtn.scope[tel].npix];
        rtn.scope[tel].ped		  = new float [rtn.scope[tel].npix];
        rtn.scope[tel].exnoise	  = new float [rtn.scope[tel].npix];
        rtn.scope[tel].pixinuse	  = new int [rtn.scope[tel].npix];
        rtn.scope[tel].telpixarea = new float [rtn.scope[tel].npix];
        rtn.scope[tel].pixradius  = new float [rtn.scope[tel].npix];
        //		float		 foclen			  = 1000.0 * fData->getDetectorGeo()->getFocalLength()[tel]; //Focal length in mm
        
        //Initialize the number of live pixel in the telescope
        rtn.scope[tel].nb_live_pix = 0;
        
        //Loop on the pixels
        for( int pix = 0; pix < rtn.scope[tel].npix; pix++ )
        {
            //Pixel coordinates
            rtn.scope[tel].xcam[pix] = fData->getDetectorGeo()->getX()[pix];
            rtn.scope[tel].ycam[pix] = fData->getDetectorGeo()->getY()[pix];
            
            //Excess noise
            rtn.scope[tel].exnoise[pix] = frogsPMTNoise[epoch_as_int];
            //Pixel dead or alive
            rtn.scope[tel].pixinuse[pix] = FROGS_OK;
            //exclude channels that are dead or where the integration window falls outside the readout window or where the pedvar is 0
            if( fData->getDead( fData->getData()->getHiLo()[pix] )[pix] != 0
               || fData->getCurrentSumWindow_2()[pix] == 0
               || fData->getData()->getPedvars( fData->getData()->getHiLo()[pix], fData->getCurrentSumWindow_2()[pix] )[pix] == 0 )
               {
                rtn.scope[tel].pixinuse[pix] = FROGS_NOTOK;
            }
            //Increment the number of live pixels
            if( rtn.scope[tel].pixinuse[pix] == FROGS_OK )
            {
                rtn.scope[tel].nb_live_pix++;
            }
            //Pixel effective collecting area in square degrees
            //float tmppixarea = fData->getDetectorGeo()->getTubeRadius_MM( tel )[pix] * FROGS_DEG_PER_RAD / foclen; //still being used?
            //tmppixarea = FROGS_PI * tmppixarea * tmppixarea; //still being used?
            //rtn.scope[tel].telpixarea[pix] = telarea * tmppixarea * cone_eff; //still being used?
            //Pixel radius in degree
            //rtn.scope[tel].pixradius[pix] = fData->getDetectorGeo()->getTubeRadius_MM( tel )[pix] * FROGS_DEG_PER_RAD / foclen; //(SV)
            
            //Initialize the pixel signal and pedestal width to zero
            rtn.scope[tel].q[pix]	= 0;
            rtn.scope[tel].ped[pix] = 0;
            //Set them to their values in p.e. if the d.c./p.e. factor is non zero
            if( frogsDCtoPE[epoch_as_int] != 0 )
            {
                rtn.scope[tel].q[pix] = fData->getData()->getSums2()[pix] / frogsDCtoPE[epoch_as_int];
                //rtn.scope[tel].q[pix]=fData->getData()->getSums()[pix] / frogsDCtoPE[epoch_as_int];
                if( fData->getData()->getHiLo()[pix] == 1 )
                {
                    //rtn.scope[tel].q[pix]=fData->getData()->getSums()[pix]*fData->getData()->getLowGainMultiplier()[pix]/dc2pe;
                    //rtn.scope[tel].q[pix]=fData->getData()->getSums()[pix]/dc2pe; //(SV): getLowGainMultiplier removed
                    //rtn.scope[tel].ped[pix]=fData->getData()->getPedvars(true,18)[pix]*frogs_pedwidth_correction/dc2pe;
                    rtn.scope[tel].ped[pix] =
                        fData->getData()->getPedvars( fData->getData()->getHiLo()[pix], fData->getCurrentSumWindow_2()[pix] )[pix] / frogsDCtoPE[epoch_as_int];
                }
                else
                {
                    //rtn.scope[tel].q[pix]=fData->getData()->getSums()[pix]/dc2pe;
                    //rtn.scope[tel].ped[pix]=fData->getData()->getPedvars(false,18)[pix]*frogs_pedwidth_correction/dc2pe;
                    //rtn.scope[tel].ped[pix]=fData->getData()->getPedvars()[pix]*frogs_pedwidth_correction/dc2pe;
                    rtn.scope[tel].ped[pix] =
                        fData->getData()->getPedvars( fData->getData()->getHiLo()[pix], fData->getCurrentSumWindow_2()[pix] )[pix] / frogsDCtoPE[epoch_as_int];
                }
            }
        }
        //Total number of live pixels in the array
        rtn.nb_live_pix_total = rtn.nb_live_pix_total + rtn.scope[tel].nb_live_pix;
    }
    
    //Optimization starting point todo y -> -y ??
    rtn.startpt.xs =  1.0 * fData->getShowerParameters()->fShower_Xoffset[frogsRecID]; //(SV) starting points set to ED parameters
    rtn.startpt.ys = -1.0 * fData->getShowerParameters()->fShower_Yoffset[frogsRecID]; //(SV) starting points set to ED parameters
    
    rtn.startpt.xp = fData->getShowerParameters()->fShowerXcore_SC[frogsRecID]; //(SV) starting points set to ED parameters
    rtn.startpt.yp = fData->getShowerParameters()->fShowerYcore_SC[frogsRecID]; //(SV) starting points set to ED parameters
    
    rtn.startpt.log10e = inEnergy; //(SV) inEnergy from ED
    
    if( frogsCheating == true )
    {
        //the starting parameter values for the minimization are set to MC values
        rtn.startpt.xs		 = 1.0 * fData->getShowerParameters()->MCTel_Xoff;
        rtn.startpt.ys		 = -1.0 * fData->getShowerParameters()->MCTel_Yoff;
        rtn.startpt.xp		 = fData->getShowerParameters()->MCxcore_SC;
        rtn.startpt.yp		 = fData->getShowerParameters()->MCycore_SC;
        //rtn.startpt.xp	 =	fData->getShowerParameters()->MCxcore; //(Ground Coord) useful?
        //rtn.startpt.yp	 =	fData->getShowerParameters()->MCycore; //(Ground Coord) useful?
        rtn.startpt.log10e =  fData->getShowerParameters()->MCenergy;
    }
    
    if( FROGSDEBUG )
    {
        cout << " eventNumber  " << eventNumber << endl;
        cout << " ED: ShowerSC " << fData->getShowerParameters()->fShowerXcore_SC[frogsRecID] << " " << fData->getShowerParameters()->fShowerYcore_SC[frogsRecID] << endl;
        cout << " ED: Shower   " << fData->getShowerParameters()->fShowerXcore[frogsRecID]    << " " << fData->getShowerParameters()->fShowerYcore[frogsRecID] << endl;
        cout << " MC: ShowerSC " << fData->getShowerParameters()->MCxcore_SC                  << " " << fData->getShowerParameters()->MCycore_SC << " " << fData->getShowerParameters()->MCzcore_SC << endl;
        cout << " MC: Shower   " << fData->getShowerParameters()->MCxcore                     << " " << fData->getShowerParameters()->MCycore    << " " << fData->getShowerParameters()->MCzcore << endl;
    }
    
    rtn.startpt.lambda = 1.0; //(SV) We use a fixed value by lack of information.
    
    if( rtn.startpt.log10e > 0.0 )
    {
        rtn.startpt.log10e = log10( rtn.startpt.log10e );
    }
    else
    {
        rtn.startpt.log10e = FROGS_BAD_NUMBER;
    }
    
    //Decides if the event is worth analysing.
    rtn.worthy_event = FROGS_OK;
    //Log(0.06)=-1.2; Log(0.1)=-1; Log(0.15)=-0.824; Log(0.2)=-0.699; Log(0.25)=-0.602
    //Log(0.3)=-0.523; Log(0.35)=-0.456; Log(0.4)=-0.398; Log(30)=1.477
    //Energy large enough?
    if( rtn.startpt.log10e < -1.20 )
    {
        rtn.worthy_event = FROGS_NOTOK;
    }
    //Energy small enough?
    //if(rtn.startpt.log10e>1.0)   rtn.worthy_event=FROGS_NOTOK;
    if( rtn.startpt.log10e > 1.470 )
    {
        rtn.worthy_event = FROGS_NOTOK;
    }
    //Distance of the impact point small enough?
    if( sqrt( rtn.startpt.xp * rtn.startpt.xp + rtn.startpt.yp * rtn.startpt.yp ) > 450.0 )
    {
        rtn.worthy_event = FROGS_NOTOK;
    }
    //Count the number of telescopes with more than 100dc in their image
    int ngoodimages = 0;
    for( int tel = 0; tel < rtn.ntel; tel++ )
    {
        setTelID( tel );
        if( fData->getImageParameters()->size > 100.0 )
        {
            ngoodimages = ngoodimages + 1;
        }
    }
    //Require the number of telescopes with more than 100dc to be at least 3
    if( ngoodimages < 2 )
    {
        rtn.worthy_event = FROGS_NOTOK;
    }
    //Require the image to be fully contained into the camera
    //if (fData->getImageParameters()->loss>0.) rtn.worthy_event=FROGS_NOTOK;
    
    if( rtn.worthy_event == FROGS_OK )
    {
        return rtn;
    }
    
#ifdef DIFF_EVOLUTION
    //=======================================
    //    diff. evolution algorithm
    //=======================================
    //Decides if the event is worth analysing
    rtn.worthy_event = FROGS_NOTOK;
    if( fabs( rtn.startpt.log10e-FROGS_BAD_NUMBER ) < 1E-8 )
    {
    
        //Distance of the impact point small enough?
        double dummy = sqrt( rtn.startpt.xp * rtn.startpt.xp + rtn.startpt.yp * rtn.startpt.yp );
        
        //Count the number of telescopes with more than 200dc in their image
        ngoodimages = 0;
        for( int tel = 0; tel < rtn.ntel; tel++ )
        {
            setTelID( tel );
            if( fData->getImageParameters()->size > 200.0 )
            {
                ngoodimages = ngoodimages + 1;
            }
        }
        
        //Require the number of telescopes with more than 100dc to be at least 2
        if( ngoodimages > 1 && dummy < 200.0 )
        {
            rtn.worthy_event = FROGS_OK;
        }
    }
    //=======================================
    //    diff. evolution algorithm
    //=======================================
#endif //DIFF_EVOLUTION
    
    return rtn;
    
}

void VFrogs::transformResults()
{

    if( fabs( frogsXS ) > 50 )
    {
        frogsXPGC = FROGS_BAD_NUMBER;
        frogsYPGC = FROGS_BAD_NUMBER;
        frogsZe = FROGS_BAD_NUMBER;
        frogsAz = FROGS_BAD_NUMBER;
        frogsXS_derot = FROGS_BAD_NUMBER;
        frogsYS_derot = FROGS_BAD_NUMBER;
        for( unsigned int tel = 0; tel < fData->fNTel; tel++ )
        {
            frogsR[tel] = FROGS_BAD_NUMBER;
        }
        return;
    }
    
    else
    {
    
        double pointingEl = fData->getShowerParameters()->fArrayPointing_Elevation;
        double pointingAz = fData->getShowerParameters()->fArrayPointing_Azimuth;
        
        
        //Core in Ground Coordinates
        //z coordinate in shower system s.th. z=0 in ground coordinates. See VArrayAnalyzer::fillShowerCore
        double zp = frogsYP / tan( ( getArrayPointing()->getTelElevation() ) / TMath::RadToDeg() );
        frogsXPGC = transformShowerPosition( 90. - pointingEl, pointingAz, frogsXP, frogsYP, zp, 0 );
        frogsYPGC = transformShowerPosition( 90. - pointingEl, pointingAz, frogsXP, frogsYP, zp, 1 );
        
        
        
        //Event Zenith and Azimuth
        //code copied from VArrayAnalyzer.cpp & adjusted
        double ze = 0.;
        double az = 0.;
        if( fData->getArrayPointing() )
        {
            VAstronometry::vlaDtp2s( -1.*frogsXS * TMath::DegToRad(),
                                         frogsYS * TMath::DegToRad(),
                                         getArrayPointing()->getTelAzimuth() * TMath::DegToRad(),
                                         getArrayPointing()->getTelElevation() * TMath::DegToRad(),
                                         &az, &ze );
            az *= TMath::RadToDeg();
            ze = 90. - ze*TMath::RadToDeg();
        }
        if( TMath::IsNaN( ze ) ||  TMath::IsNaN( az ) )
        {
            ze = FROGS_BAD_NUMBER;
            az = FROGS_BAD_NUMBER;
        }
        else
        {
            az = VAstronometry::vlaDranrm( az * TMath::DegToRad() ) * TMath::RadToDeg();
        }
        frogsZe = ze;
        frogsAz = az;
        
        
        // calculate derotated shower directions. Note that frogs uses a different coordinate system than the rest of ED.
        if( !fReader->isMC() )
        {
            double iUTC = 0.;
            double xrot = 0.;
            double yrot = 0.;
            if( fData->getArrayPointing() )
            {
                iUTC = VSkyCoordinatesUtilities::getUTC( getShowerParameters()->MJD, getShowerParameters()->time );
                fData->getArrayPointing()->derotateCoords( iUTC, frogsXS, 1.*frogsYS, xrot, yrot );
            }
            
            frogsXS_derot = xrot;
            frogsYS_derot = 1.*yrot;
        }
        else
        {
            frogsXS_derot = frogsXS;
            frogsYS_derot = frogsYS;
        }
        
        
        //code copied from VTableLookupDataHandler::calcDistances and adjusted
        for( unsigned int tel = 0; tel < fData->fNTel; tel++ )
        {
            if( fabs( frogsXP ) < 999 )
            {
                // frogsR[tel] = VUtilities::line_point_distance( frogsYPGC, -1.*frogsXPGC, 0., frogsZe, frogsAz, getDetectorGeo()->getTelYpos()[tel], -1.*getDetectorGeo()->getTelXpos()[tel], getDetectorGeo()->getTelZpos()[tel] );
                //note: line_point_distance expects x-axis pointing north etc.
                //	In our coodinate system, the y-axis points north, the x-axis points east. So we have to rotate by 90 degrees.
                
                double xfield = transformPosition( 90. - pointingEl, pointingAz, getDetectorGeo()->getTelXpos()[tel], getDetectorGeo()->getTelYpos()[tel], getDetectorGeo()->getTelZpos()[tel], 0, false );
                double yfield =	transformPosition( 90. - pointingEl, pointingAz, getDetectorGeo()->getTelXpos()[tel], getDetectorGeo()->getTelYpos()[tel], getDetectorGeo()->getTelZpos()[tel], 1, false );
                
                frogsR[tel] = TMath::Sqrt( TMath::Power( xfield - frogsXP, 2 ) +  TMath::Power( yfield - frogsYP, 2 ) );
                
            }
            else
            {
                frogsR[tel] = FROGS_BAD_NUMBER;
            }
        }
    }
}

