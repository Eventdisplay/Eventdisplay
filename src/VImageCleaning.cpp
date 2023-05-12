/*! \class VImageCleaning

     collection of different image cleaning methods

 */

#include "VImageCleaning.h"

VImageCleaning::VImageCleaning( VEvndispData* iData )
{
    fData = iData;
    
    fProb4nnCurves = 0;
    fProb3nnrelCurves = 0;
    fProb2plus1Curves = 0;
    fProb2nnCurves = 0;
    fProbBoundCurves = 0;
    fIPRgraphs = 0;
    
    // initialize NN cleaning (TIMENEXTNEIGHBOUR)
    kInitNNImageCleaning = false;
    if( fData && fData->getImageCleaningParameter()->getImageCleaningMethod() == "TIMENEXTNEIGHBOUR" )
    {
        kInitNNImageCleaning = InitNNImageCleaning();
        fIPRgraphs = new TObjArray( VDST_MAXTELTYPES );
        fIPRgraphs_xmax.assign( VDST_MAXTELTYPES, 0. );
    }
    
    fWriteGraphToFileRecreate = true;
    fFakeImageProb = 0.;
    vector< bool > i_tempB;
    for( unsigned int i = 0; i < VDST_MAXTELTYPES; i++ )
    {
        kInitNNImgClnPerTelType[i] = false;
        fMinRate[i] = 0.;
        fifActiveNN.push_back( i_tempB );
    }
    fIPR_save_mincharge = -99.;
    fIPR_save_dT_from_probCurve = -99.;
    fIPR_save_telid = 99;
    fIPR_save_ProbCurve_par1 = -99.;
    fIPR_save_ProbCurve_par2 = -99.;
    
}

/*
 * simple error messager
 *
*/
void VImageCleaning::printDataError( string iFunctionName )
{
    cout << iFunctionName;
    cout << " error: no pointer to data class set" << endl;
    exit( EXIT_FAILURE );
}

/*!
  tailcut cleaning with fixed thresholds
   \par hithresh image threshold
   \par lothresh border threshold
   \par brightthresh bright pixel threshold
*/

void VImageCleaning::cleanImageFixed( VImageCleaningRunParameter* iImageCleaningParameters )
{
    ///////////////////////////////////////////////
    // setting of image cleaning run parameters
    if( !iImageCleaningParameters )
    {
        return;
    }
    double hithresh     = iImageCleaningParameters->fimagethresh;
    double lothresh     = iImageCleaningParameters->fborderthresh;
    double brightthresh = iImageCleaningParameters->fbrightnonimagetresh;
    
    cleanImageFixed( hithresh, lothresh, brightthresh );
}

void VImageCleaning::cleanImageFixed( double hithresh, double lothresh, double brightthresh )
{

    if( !fData )
    {
        printDataError( "VImageCleaning::cleanImageFixed" );
    }
    
    fData->setImage( false );
    fData->setBorder( false );
    fData->setBrightNonImage( false );
    fData->setImageBorderNeighbour( false );
    unsigned int i_nchannel = fData->getNChannels();
    
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        if( fData->getSums()[i] > hithresh )
        {
            if( fData->getDetectorGeo()->getAnaPixel()[i] > 0 && !fData->getDead( i, fData->getHiLo()[i] ) )
            {
                fData->setImage( i, true );
                for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[i].size(); j++ )
                {
                    unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                    if( k < fData->getImage().size() && fData->getSums()[k] > lothresh && !fData->getImage()[k] )
                    {
                        fData->setBorder( k, true );
                    }
                }
            }
        }
        if( fData->getSums()[i] > brightthresh )
        {
            if( fData->getDetectorGeo()->getAnaPixel()[i] > 0 && !fData->getDead( i, fData->getHiLo()[i] ) )
            {
                fData->setBrightNonImage( i, true );
            }
        }
    }
    
    // (preli) set the trigger vector in MC case (preli)
    // trigger vector are image/border tubes
    if( fData->getReader()->getDataFormatNum() == 1 || fData->getReader()->getDataFormatNum() == 4
            || fData->getReader()->getDataFormatNum() == 6 )
    {
        fData->getReader()->setTrigger( fData->getImage(), fData->getBorder() );
    }
    // (end of preli)
    if( fData->getRunParameter()->frecoverImagePixelNearDeadPixel )
    {
        recoverImagePixelNearDeadPixel();
    }
    if( fData->getRunParameter()->fFillImageBorderNeighbours )
    {
        fillImageBorderNeighbours();
    }
}


/*!

  signal-to-noise tailcut cleaning

   \par hithresh image threshold
   \par lothresh border threshold
   \par brightthresh bright pixel threshold

*/
void VImageCleaning::cleanImagePedvars( VImageCleaningRunParameter* iImageCleaningParameters )
{
    ///////////////////////////////////////////////
    // setting of image cleaning run parameters
    if( !iImageCleaningParameters )
    {
        return;
    }
    double hithresh     = iImageCleaningParameters->fimagethresh;
    double lothresh     = iImageCleaningParameters->fborderthresh;
    double brightthresh = iImageCleaningParameters->fbrightnonimagetresh;
    
    if( !fData )
    {
        printDataError( "VImageCleaning::cleanImagePedvars" );
    }
    
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::cleanImagePedvars " << fData->getTelID() << endl;
    }
    
    fData->setImage( false );
    fData->setBorder( false );
    fData->setBrightNonImage( false );
    fData->setImageBorderNeighbour( false );
    unsigned int i_nchannel = fData->getNChannels();
    double i_pedvars_i = 0.;
    double i_pedvars_k = 0.;
    unsigned int k = 0;
    
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        if( fData->getDetectorGeo()->getAnaPixel()[i] < 1 || fData->getDead( i, fData->getHiLo()[i] ) )
        {
            continue;
        }
        i_pedvars_i = fData->getPedvars( fData->getCurrentSumWindow()[i], fData->getHiLo()[i] )[i];
        
        if( fData->getSums()[i] > hithresh * i_pedvars_i )
        {
            fData->setImage( i, true );
            fData->setBorder( i, false );
            for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[i]; j++ )
            {
                k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( k < i_nchannel )
                {
                    i_pedvars_k = fData->getPedvars( fData->getCurrentSumWindow()[k], fData->getHiLo()[k] )[k];
                    if( !fData->getImage()[k] && fData->getSums()[k] > lothresh * i_pedvars_k )
                    {
                        fData->setBorder( k, true );
                    }
                }
            }
        }
        if( fData->getSums()[i] > brightthresh  * i_pedvars_i )
        {
            fData->setBrightNonImage( i, true );
        }
    }
    
    // (preli) set the trigger vector in MC case (preli)
    // trigger vector are image/border tubes
    if( fData->getReader() )
    {
        if( fData->getReader()->getDataFormatNum() == 1 || fData->getReader()->getDataFormatNum() == 4
                || fData->getReader()->getDataFormatNum() == 6 )
        {
            fData->getReader()->setTrigger( fData->getImage(), fData->getBorder() );
        }
    }
    // (end of preli)
    
    recoverImagePixelNearDeadPixel();
    fillImageBorderNeighbours();
}

/*!
   simple time cleaning

   Image cleaning method keyword in runparameter file: "TIMETWOLEVEL"

   \par hithresh     image threshold
   \par lothresh     border threshold
   \par brightthresh bright pixel threshold (for which time difference is calculated)
   \par timediff time constraint between next neighbor pixels

*/
void VImageCleaning::cleanImagePedvarsTimeDiff( VImageCleaningRunParameter* iImageCleaningParameters )
{
    ///////////////////////////////////////////////
    // setting of image cleaning run parameters
    if( !iImageCleaningParameters )
    {
        return;
    }
    double hithresh     = iImageCleaningParameters->fimagethresh;
    double lothresh     = iImageCleaningParameters->fborderthresh;
    double brightthresh =  iImageCleaningParameters->fbrightnonimagetresh;
    double timediff      = iImageCleaningParameters->ftimediff;
    
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::cleanImagePedvarsTimeDiff " << fData->getTelID() << endl;
    }
    
    fData->setImage( false );
    fData->setBorder( false );
    fData->setBrightNonImage( false );
    fData->setImageBorderNeighbour( false );
    unsigned int i_nchannel = fData->getNChannels();
    double i_pedvars_i = 0.;
    double i_pedvars_k = 0.;
    double i_pedvars_l = 0.;
    unsigned int k = 0;
    unsigned int l = 0;
    
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        // check if pixel is valid
        if( fData->getDetectorGeo()->getAnaPixel()[i] < 1 || fData->getDead( fData->getHiLo()[i] )[i] )
        {
            continue;
        }
        i_pedvars_i = fData->getPedvars( fData->getCurrentSumWindow()[i], fData->getHiLo()[i] )[i];
        
        //////////////////
        // image pixel
        if( fData->getSums()[i] > hithresh * i_pedvars_i )
        {
            // loop over all neighbours
            for( unsigned int z = 0; z < fData->getDetectorGeo()->getNNeighbours()[i]; z++ )
            {
                l = fData->getDetectorGeo()->getNeighbours()[i][z];
                if( l < i_nchannel )
                {
                    i_pedvars_l = fData->getPedvars( fData->getCurrentSumWindow()[l], fData->getHiLo()[l] )[l];
                    // image pixel has:
                    //   - one neighbour pixels above the border threshold AND
                    //     with a time difference smaller than timediff
                    if( fData->getSums()[l] > lothresh * i_pedvars_l && fabs( fData->getPulseTime()[i] - fData->getPulseTime()[l] ) < timediff )
                    {
                        fData->setImage( i, true );
                    }
                }
                // border pixel
                fData->setBorder( i, false );
                for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[i]; j++ )
                {
                    k = fData->getDetectorGeo()->getNeighbours()[i][j];
                    if( k < i_nchannel )
                    {
                        i_pedvars_k = fData->getPedvars( fData->getCurrentSumWindow()[k], fData->getHiLo()[k] )[k];
                        // border pixel has:
                        //   - one neighbour pixels above the border threshold AND
                        //     with a time difference smaller than timediff
                        if( !fData->getImage()[k] && fData->getSums()[k] > lothresh * i_pedvars_k
                                &&  fabs( fData->getPulseTime()[i] - fData->getPulseTime()[k] ) < timediff )
                        {
                            fData->setBorder( k, true );
                        }
                    }
                }
            }
        }
        // bright pixel threshold
        if( fData->getSums()[i] > brightthresh  * i_pedvars_i )
        {
            fData->setBrightNonImage( i, true );
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    // (preli) set the trigger vector in MC case (preli)
    // trigger vector are image/border tubes
    if( fData->getReader() )
    {
        if( fData->getReader()->getDataFormatNum() == 1 || fData->getReader()->getDataFormatNum() == 4
                || fData->getReader()->getDataFormatNum() == 6 )
        {
            fData->getReader()->setTrigger( fData->getImage(), fData->getBorder() );
        }
    }
    // (end of preli)
    
    recoverImagePixelNearDeadPixel();
    fillImageBorderNeighbours();
}

//*****************************************************************************************************
// NN image cleaning (Maxim)
//
//  Searches N-fold connected groups (4nn close-packed group) and 2+1 sparse groups
//  with dynamic cut in  GroupThresh-GroupCoincTime phase space
//
//  NOTE: IPR curves, used for cuts,  should be obtained with the SAME signal extraction method
//        as the image under study
// preliminary

/*
 * calculate maximum of a float array
 */
void VImageCleaning::LocMax( int n, float* ptr, float& max )
{
    if( n <= 0 )
    {
        return;
    }
    max = ptr[0];
    for( Int_t i = 1; i < n; i++ )
    {
        if( max < ptr[i] )
        {
            max = ptr[i];
        }
    }
}

/*
 * calculate minimum of a float array
 */
void VImageCleaning::LocMin( int n, float* ptr, float& min ) //ptr[i]>0
{
    if( n <= 0 )
    {
        return;
    }
    min = ptr[0];
    for( Int_t i = 1; i < n; i++ )
    {
        if( min > ptr[i] )
        {
            min = ptr[i];
        }
    }
}

//*****************************************************************************************************
// start NN image cleaning

/*
 *
 * init TIMENEXTNEIGHBOUR image cleaning method
 *
 */
bool VImageCleaning::InitNNImageCleaning()
{
    const int types = VDST_MAXTELTYPES;
    // arrays with rate contour curves
    // (length of number of telescope types)
    fProb4nnCurves    = new TObjArray( types );
    fProb3nnrelCurves = new TObjArray( types );
    fProb2plus1Curves = new TObjArray( types );
    fProb2nnCurves    = new TObjArray( types );
    fProbBoundCurves  = new TObjArray( types );
    
    // init IPR arrays
    // IPR[teltypes][IPRdim]
    IPR = new float* [VDST_MAXTELTYPES];
    for( int t = 0; t < VDST_MAXTELTYPES; t++ )
    {
        IPR[t] = new float[fIPRdim];
    }
    for( int t = 0; t < VDST_MAXTELTYPES; t++ )
    {
        for( unsigned int th = 0; th < fIPRdim; th++ )
        {
            IPR[t][th] = 0.;
        }
    }
    
    return true;
}


/*
    get IPR graphs from different sources (files), calculate probability contours

*/
bool VImageCleaning::InitNNImgClnPerTelType( unsigned int teltype )
{
    // reading IPR graphs
    TGraphErrors* IPRgraph = NULL;
    // get ready filled IPR from DST or pedestal file
    if( fData->getRunParameter()->ifReadIPRfromDSTFile )
    {
        cout << "VImageCleaning::InitNNImgClnPerTelType( int type ): getting IPR graph for (simtel) teltype: ";
        cout << teltype << " (type " << ( int )fData->getTelType( fData->getTelID() );
        cout << ", summation window " << fData->getSumWindow() << ") from DST file" << endl;
        IPRgraph = fData->getIPRGraph();
        if( IPRgraph )
        {
            cout << "\t found graph: " << IPRgraph->GetName() << endl;
        }
    }
    else if( !fData->getRunParameter()->ifReadIPRfromDatabase )
    {
        cout << "VImageCleaning::InitNNImgClnPerTelType( int type ): getting IPR graph for teltype: ";
        cout << teltype << " from external file (ped or DST file)" << endl;
        
        IPRgraph = fData->getIPRGraph();
        if( IPRgraph )
        {
            IPRgraph->SetName( Form( "IPRchargeTelType%d", teltype ) );
        }
    }
    // get IPR from database
    else
    {
        cout << "VImageCleaning::InitNNImgClnPerTelType( int type ): getting IPR graph for teltype:" << teltype << " from external IPR database" << endl;
        TFile* fgraphs = new TFile( fData->getRunParameter()->fIPRdatabaseFile, "READ" );
        TString  gname = Form( "IPRchargeTelType%d_TelID%d", ( int )fData->getTelType( fData->getTelID() ), ( int )fData->getTelID() );
        IPRgraph = ( TGraphErrors* )fgraphs->Get( gname.Data() );
        fgraphs->Close();
        fgraphs->Delete();
    }
    if( !IPRgraph )
    {
        cout << "VImageCleaning::InitNNImgClnPerTelType( int type ) ERROR: IPR graph is NULL for teltype " << teltype << " !!!" << endl;
        printDataError( "" );
    }
    fIPRgraphs->AddAt( IPRgraph, teltype );
    if( teltype < fIPRgraphs_xmax.size() )
    {
        fIPRgraphs_xmax[teltype] = IPRgraph->GetXaxis()->GetXmax();
    }
    else
    {
        cout << "Error filling IPR graphs required for NN cleaning ";
        cout << "(" << teltype << ")" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    cout << "IPR for TelType (TrigSim, TIMENEXTNEIGHBOUR): " << teltype << " read: " << IPRgraph->GetName();
    cout << endl;
    
    //initializing probability curves
    unsigned int NgroupTypes = 0;
    if( teltype < fifActiveNN.size() )
    {
        for( unsigned int i = 0; i < fifActiveNN[teltype].size(); i++ )
        {
            if( fifActiveNN[teltype][i] )
            {
                NgroupTypes++;
            }
        }
    }
    if( NgroupTypes == 0 )
    {
        cout << "VImageCleaning::InitNNImgClnPerTelType( int type ) ERROR:";
        cout << "no NN groups type selected for search (check image cleaning input card) ...exit" << endl;
        printDataError( "VImageCleaning::InitNNImgClnPerTelType( int type )" );
    }
    float upSampling = 1.;
    if( fData->getDigitalFilterMethod() > 0 && fData->getDigitalFilterUpSample() > 0 )
    {
        upSampling = ( float )fData->getDigitalFilterUpSample();
        // don't allow down sampling
        if( upSampling < 1. )
        {
            upSampling = 1.;
        }
    }
    float SimTime = fData->getDetectorGeo()->getLengthOfSampleTimeSlice( fData->getTelID() ) / upSampling;
    if( fData->getSearchWindowLast() < fData->getNSamples() )
    {
        SimTime *= ( fData->getSearchWindowLast() - fData->getSumFirst() );
    }
    else
    {
        SimTime *= ( fData->getNSamples() * upSampling - fData->getSumFirst() );
    }
    if( setExplicitSampleTimeSlice )
    {
        // simple peak sensing: sim_telarray uses the maximum bin only
        // (searched for in the window given by SimTime)
        // The values for sampleTimeSlice and nBinsADC are set in the cleaning parameter file
        // For example, for the currect (Apr 17) ASTRI simulation, it is set in sim_telarray as
        // fadc_mhz = 500 % MHz ==> sampleTimeSlice = 2 ns
        // fadc_sum_bins = nBinsADC = 25 % Number of ADC time intervals actually summed up.
        SimTime = sampleTimeSlice * nBinsADC;
    }
    float fMinRate = fFakeImageProb / ( SimTime * 1E-9 * float( NgroupTypes ) ); //ns
    
    // number of combinations
    // (this number is scaled to the size of the camera later)
    float CombFactor[5] = {60000.*20., 950000., 130000., 12000., 2.}; //for 2400 pixels
    float ChargeMax = fIPRgraphs_xmax[teltype];
    
    // define rate contour functions
    fProb4nnCurves->AddAt( defineRateContourFunction( teltype, "ProbCurve4nn", fMinRate, 4, CombFactor[0], 0, ChargeMax ), ( int )teltype );
    fProb3nnrelCurves->AddAt( defineRateContourFunction( teltype, "ProbCurve3nnrel", fMinRate, 3, CombFactor[2], 0, ChargeMax ), ( int )teltype );
    fProb2plus1Curves->AddAt( defineRateContourFunction( teltype, "ProbCurve2plus1", fMinRate, 3, CombFactor[1], 0, ChargeMax ), ( int )teltype );
    fProb2nnCurves->AddAt( defineRateContourFunction( teltype, "ProbCurve2nn", fMinRate, 2, CombFactor[3], 0, ChargeMax ), ( int )teltype );
    fProbBoundCurves->AddAt( defineRateContourBoundFunction( teltype, "ProbCurveBound", fMinRate, 4.0, CombFactor[4], 0, ChargeMax ), ( int )teltype );
    
    cout << "Fake image probability: " << fFakeImageProb << " NgroupTypes: " << NgroupTypes;
    cout << " teltype " << teltype << " ChargeMax:" << ChargeMax << " Min rate: " << fMinRate;
    cout << std::endl;
    
    /////////////////////////////////////////////////////////////////////////
    // writing file with graphs and prob curves if required in input card
    if( fData->getRunParameter()->ifWriteGraphsToFile )
    {
        string iFileN = "UPDATE";
        if( fWriteGraphToFileRecreate )
        {
            iFileN = "RECREATE";
        }
        TFile* fgraphs = new TFile( fData->getRunParameter()->fNNGraphsFile, iFileN.c_str() );
        if( fgraphs->IsZombie() )
        {
            printDataError( "VImageCleaning::InitNNImgClnPerTelType( int type ), error opening NNgraphs file" );
        }
        // IPR graphs
        IPRgraph->Write();
        if( fWriteGraphToFileRecreate )
        {
            // probability curves (note that the x-axis is not charge!)
            fProb4nnCurves->Write();
            fProb3nnrelCurves->Write();
            fProb2plus1Curves->Write();
            fProb2nnCurves->Write();
            fProbBoundCurves->Write();
        }
        fWriteGraphToFileRecreate = false;
        
        // probablity curves
        writeProbabilityCurve( ( TGraph* )fIPRgraphs->At( teltype ), ( TF1* )fProb4nnCurves->At( teltype ), fMinRate );
        writeProbabilityCurve( ( TGraph* )fIPRgraphs->At( teltype ), ( TF1* )fProb3nnrelCurves->At( teltype ), fMinRate );
        writeProbabilityCurve( ( TGraph* )fIPRgraphs->At( teltype ), ( TF1* )fProb2plus1Curves->At( teltype ), fMinRate );
        writeProbabilityCurve( ( TGraph* )fIPRgraphs->At( teltype ), ( TF1* )fProb2nnCurves->At( teltype ), fMinRate );
        
        // close output file
        fgraphs->Close();
        std::cout << "VImageCleaning::InitNNImgClnPerTelType( type ): writing graph root file to: " << fgraphs->GetName() << std::endl;
    }
    
    return true;
}

/*
 * fill and write contour plot of delta T vs minimal charge
 *
 * (Equation 2 and Fig 2 in ICRC 2013 proceedings)
 *
 */
void VImageCleaning::writeProbabilityCurve( TGraph* iIPR, TF1* iProb, double iRate )
{
    if( !iIPR || !iProb )
    {
        return;
    }
    TGraph* iGPP = new TGraph( iIPR->GetN() );
    double x = 0.;
    double y = 0.;
    for( int i = 0; i < iIPR->GetN(); i++ )
    {
        iIPR->GetPoint( i, x, y );
        float valIPR = iIPR->Eval( x, 0, "" );
        if( valIPR < 100. )
        {
            valIPR = 100.;
        }
        iGPP->SetPoint( i, x, iProb->Eval( valIPR ) );
    }
    char hname[400];
    sprintf( hname, "graph%s", iProb->GetName() );
    iGPP->SetName( hname );
    sprintf( hname, "Rate contour for %1.2e Hz (fake image prob:%1.2e)", iRate, fFakeImageProb );
    iGPP->SetTitle( hname );
    iGPP->GetXaxis()->SetTitle( "Threshold,  [FADC counts]" );
    iGPP->GetYaxis()->SetTitle( "#Delta t [ns]" );
    iGPP->Write();
}

/*
 * scale combinatorial factors by number of channels
 *
 *
*/
void VImageCleaning::ScaleCombFactors( unsigned int type, float scale )
{
    float CombFactor[5] = {60000.*20., 950000., 130000., 6000., 2.}; //for 2400 pixels
    TF1* f4nn   = ( TF1* )fProb4nnCurves->At( type );
    TF1* f3nnrel = ( TF1* )fProb3nnrelCurves->At( type );
    TF1* f2plus1 = ( TF1* )fProb2plus1Curves->At( type );
    TF1* f2nn   = ( TF1* )fProb2nnCurves->At( type );
    f4nn->SetParameter( 2, CombFactor[0]*scale );
    f3nnrel->SetParameter( 2, CombFactor[2]*scale );
    f2plus1->SetParameter( 2, CombFactor[1]*scale );
    f2nn->SetParameter( 2, CombFactor[3]*scale );
}

/*
 * reset combinatorial factors to default values
 *
*/
void VImageCleaning::ResetCombFactors( unsigned int type )
{
    ScaleCombFactors( type, fData->getNChannels() / 2400. ); // for this amount of pixels
}

bool VImageCleaning::BoundarySearch( unsigned int teltype, float thresh, TF1* fProbCurve, float refdT, int refvalidity, int idx )
{
    //idx - should be next neigbour of core!!!
    //skip core pix
    if( ( VALIDITYBUF[idx] > 1.9 && VALIDITYBUF[idx] < 6.1 ) )
    {
        return false;
    }
    
    // check for valid pixel number
    if( idx >= ( int )fData->getDetectorGeo()->getNeighbours().size() )
    {
        return false;
    }
    
    // get IPR graph
    TGraph* iIPR = ( TGraph* )fIPRgraphs->At( teltype );
    if( !iIPR )
    {
        cout << "VImageCleaning::BoundarySearch error, no IPR graph for telescope type " << teltype << endl;
        return 0;
    }
    float iIPR_max = fIPRgraphs_xmax[teltype];
    
    //    float TimeForReSearch = 0.;
    bool iffound = false;
    Int_t n = 0;
    float time = 0.;
    
    // reftime from core pixels
    for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[idx].size(); j++ )
    {
        const Int_t idx2 = fData->getDetectorGeo()->getNeighbours()[idx][j];
        Float_t t = TIMES[idx2];
        if( t > 0. && VALIDITYBUF[idx2] > 1.9 && VALIDITYBUF[idx2] < 5.1 )
        {
            time += t;
            n++;
        }
    }
    // check THIS pix and ring of neighbours
    if( n > 0.5 )
    {
        float TimeForReSearch = time / float( n );
        float dT = fabs( TIMES[idx] - TimeForReSearch );
        float times[2] = {refdT, dT};
        float maxtime = 1E6;
        LocMax( 2, times, maxtime );
        float charges[2] = {thresh, INTENSITY[idx]};
        float mincharge = 0;
        LocMin( 2, charges, mincharge );
        
        if( NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, maxtime, CoincWinLimit, true )
                && VALIDITY[idx] > 0.5 )
        {
            if( VALIDITYBOUND[idx] != refvalidity )
            {
                VALIDITYBOUND[idx] = refvalidity;
            }
            iffound = true;
        }
        
        for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[idx].size(); j++ )
        {
            const Int_t idx2 = fData->getDetectorGeo()->getNeighbours()[idx][j];
            if( ( TIMES[idx2] > 0. && VALIDITYBUF[idx2] > 1.9 && VALIDITYBUF[idx2] < 5.1 ) || VALIDITYBOUND[idx2] == refvalidity )
            {
                continue;
            }
            
            float dT2 = fabs( TIMES[idx2] - TimeForReSearch );
            times[1] = dT2;
            LocMax( 2, times, maxtime );
            charges[1] = INTENSITY[idx2];
            LocMin( 2, charges, mincharge );
            
            if( NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, maxtime, CoincWinLimit, true )
                    && VALIDITY[idx2] > 0.5 )
            {
                VALIDITYBOUND[idx2] = refvalidity;
                iffound = true;
            }
        }
    }
    return iffound;
}

/*
 * NN search
 *
 * if Nfold = 3 it will search for 2nn+1, including sparse groups (with the empty pix in between)
 *
*/
unsigned int VImageCleaning::NNGroupSearchProbCurve( unsigned int type, TF1* fProbCurve, float PreCut )
{
    if( !fProbCurve )
    {
        return 0;
    }
    
    // Nfold (e.g. 2 or 3)
    int NN = ( int )fProbCurve->GetParameter( 1 );
    
    TGraph* iIPR = ( TGraph* )fIPRgraphs->At( type );
    if( !iIPR )
    {
        cout << "VImageCleaning::NNGroupSearchProbCurve error, no IPR graph for telescope type " << type << endl;
        return 0;
    }
    float iIPR_max = fIPRgraphs_xmax[type];
    
    // (GM) unclear why this is hardwired here
    int NSBpix = 5;
    
    //////////////////////////////
    // loop over all pixels
    int numpix = fData->getDetectorGeo()->getNumChannels();
    for( int PixNum = 0; PixNum < numpix; PixNum++ )
    {
        // check validity of a pixel and apply pre cut on charge
        // if not: skip
        if( VALIDITY[PixNum] < 0.5 || INTENSITY[PixNum] < PreCut )
        {
            continue;
        }
        
        // access neighbour list and loop over all neighbours
        if( PixNum >= ( int )fData->getDetectorGeo()->getNeighbours().size() )
        {
            continue;
        }
        for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[PixNum].size(); j++ )
        {
            const Int_t PixNum2 = fData->getDetectorGeo()->getNeighbours()[PixNum][j];
            if( PixNum2 < 0 )
            {
                continue;
            }
            // apply validity and pre-cut to neighbour pixel
            if( VALIDITY[PixNum2] < 0.5 || INTENSITY[PixNum2] < PreCut )
            {
                continue;
            }
            
            // time difference between pixel and its neighbour
            Double_t dT = fabs( TIMES[PixNum] - TIMES[PixNum2] );
            
            // get minimum between pixel and neighbour charge
            float charges[2] = {INTENSITY[PixNum], INTENSITY[PixNum2]};
            float mincharge = 0;
            LocMin( 2, charges, mincharge );
            
            // apply charge and time cut
            if( !NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, dT, CoincWinLimit ) )
            {
                continue;
            }
            
            // validity of pixel and neighbour to 2
            if( VALIDITYBUF[PixNum] < 2.9 )
            {
                VALIDITYBUF[PixNum] = 2;
            }
            if( VALIDITYBUF[PixNum2] < 2.9 )
            {
                VALIDITYBUF[PixNum2] = 2;
            }
            // for NN 2 multiplicity: this is the end
            if( NN == 2 )
            {
                break;
            }
            
            ///////////////////////////////////////////////////
            // Boundary search (sparse group 2+1)
            // *************************
            if( VALIDITYBUF[PixNum2] == 2 && NN == 3 )
            {
                bool iffound = false;
                for( unsigned int k = 0; k < fData->getDetectorGeo()->getNeighbours()[PixNum].size(); k++ )
                {
                    if( BoundarySearch( type, mincharge, fProbCurve, dT, 3, fData->getDetectorGeo()->getNeighbours()[PixNum][k] ) )
                    {
                        iffound = true;
                    }
                }
                if( PixNum2 >= ( int )fData->getDetectorGeo()->getNeighbours().size() )
                {
                    continue;
                }
                for( unsigned int k = 0; k < fData->getDetectorGeo()->getNeighbours()[PixNum2].size(); k++ )
                {
                    if( BoundarySearch( type, mincharge, fProbCurve, dT, 3, fData->getDetectorGeo()->getNeighbours()[PixNum2][k] ) )
                    {
                        iffound = true;
                    }
                }
                if( iffound )
                {
                    if( VALIDITYBUF[PixNum] < 3.9 )
                    {
                        VALIDITYBUF[PixNum] = 3;
                    }
                    if( VALIDITYBUF[PixNum2] < 3.9 )
                    {
                        VALIDITYBUF[PixNum2] = 3;
                    }
                }
            }
            if( NN == 3 )
            {
                break;
            }
            
            ////////////////////////////////////////
            // *************************
            // NN = 4 and beyond
            if( NSBpix > 2 )
            {
                Double_t x = fData->getDetectorGeo()->getX()[PixNum2];
                Double_t y = fData->getDetectorGeo()->getY()[PixNum2];
                
                Int_t idxm = -1;
                Int_t idxp = -1;
                Int_t nn = 0;
                for( unsigned int kk = 0; kk < fData->getDetectorGeo()->getNeighbours()[PixNum].size(); kk++ )
                {
                    const Int_t k = fData->getDetectorGeo()->getNeighbours()[PixNum][kk];
                    if( k < 0 )
                    {
                        continue;
                    }
                    Double_t xx = x - fData->getDetectorGeo()->getX()[k];
                    Double_t yy = y - fData->getDetectorGeo()->getY()[k];
                    
                    Double_t dist = sqrt( xx * xx + yy * yy );
                    // assume that all pixel have the same tube radius
                    Double_t diam = 2.*fData->getDetectorGeo()->getTubeRadius()[1];
                    if( dist > 0.01 * diam && dist < 1.1 * diam )
                    {
                        if( nn )
                        {
                            idxp = k;
                        }
                        else
                        {
                            idxm = k;
                        }
                        nn++;
                    }
                }
                if( idxp < 0 )
                {
                    continue;
                }
                if( VALIDITY[idxp] < 0.5 || INTENSITY[idxp] < PreCut )
                {
                    continue;
                }
                
                dT = fabs( TIMES[PixNum] - TIMES[idxp] );
                double dt2 = fabs( TIMES[PixNum2] - TIMES[idxp] );
                float charges2[2] = {mincharge, INTENSITY[idxp]};
                LocMin( 2, charges2, mincharge );
                float times2[2] = {( float )dT, ( float )dt2 };
                float maxtime = 1E6;
                LocMax( 2, times2, maxtime );
                // apply charge and time cut
                if( !NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, maxtime, CoincWinLimit ) )
                {
                    continue;
                }
                
                if( VALIDITYBUF[PixNum] < 3.9 )
                {
                    VALIDITYBUF[PixNum] = 3;
                }
                if( VALIDITYBUF[PixNum2] < 3.9 )
                {
                    VALIDITYBUF[PixNum2] = 3;
                }
                if( VALIDITYBUF[idxp] < 3.9 )
                {
                    VALIDITYBUF[idxp] = 3;
                }
                if( NN == 3 )
                {
                    break;
                }
                if( NSBpix > 3 )
                {
                    if( idxm < 0 )
                    {
                        continue;
                    }
                    Double_t q3 = VALIDITY[idxm];
                    if( q3 < 0.5 || INTENSITY[idxm] < PreCut )
                    {
                        continue;
                    }
                    dT = fabs( TIMES[PixNum] - TIMES[idxm] );
                    dt2 = fabs( TIMES[PixNum2] - TIMES[idxm] );
                    double dt3 = fabs( TIMES[idxp] - TIMES[idxm] );
                    float charges3[2] = {mincharge, INTENSITY[idxm]};
                    LocMin( 2, charges3, mincharge );
                    float times3[4] = { ( float )maxtime, ( float )dT, ( float )dt2, ( float )dt3 };
                    LocMax( 4, times3, maxtime );
                    
                    // apply charge and time cut
                    if( !NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, maxtime, CoincWinLimit ) )
                    {
                        continue;
                    }
                    
                    VALIDITYBUF[PixNum] = 4;
                    VALIDITYBUF[PixNum2] = 4;
                    VALIDITYBUF[idxp] = 4;
                    VALIDITYBUF[idxm] = 4;
                    
                }
            }
        }
        
    }//end of for() loop
    
    /////////////////////////////////////////////
    // calculate number of NN groups
    
    Float_t nn2 = 0.;
    Float_t nn3 = 0.;
    Float_t nn4 = 0.;
    for( Int_t i = 0; i < numpix; i++ )
    {
        float q = VALIDITYBUF[i];
        if( q > 1.5 && q < 2.1 )
        {
            nn2 += 1. / 2.;
        }
        if( q > 2.5 && q < 3.1 )
        {
            nn3 += 1. / 3.;
        }
        if( q > 3.5 && q < 4.1 )
        {
            nn4 += 1. / 4.;
        }
        if( VALIDITYBUF[i] < 1.9 )
        {
            VALIDITYBUF[i] = -1;
        }
    }
    
    unsigned int ngroups = 0;
    if( NN == 2 )
    {
        ngroups = int( nn2 + 0.5 );
    }
    else if( NN == 3 )
    {
        ngroups = int( nn3 + 0.5 );
    }
    else if( NN == 4 )
    {
        ngroups = int( nn4 + 0.5 );
    }
    return ngroups;
}

/*
 * apply cut in charge and time space
 *
 */
bool VImageCleaning::NNChargeAndTimeCut(
    TGraph* iIPR, float iIPR_max, TF1* fProbCurve,
    float mincharge, float dT,
    float iCoincWinLimit,
    bool bInvert )
{
    if( !iIPR || !fProbCurve )
    {
        return false;
    }
    float valDT = 0.;
    // use previous result if charge and all other
    // parameter are the same
    if( fIPR_save_mincharge > -90.
            && fIPR_save_telid == fData->getTelID()
            && TMath::Abs( fIPR_save_ProbCurve_par1 - fProbCurve->GetParameter( 1 ) ) < 1.e-3
            && TMath::Abs( fIPR_save_ProbCurve_par2 - fProbCurve->GetParameter( 2 ) ) < 1.e-3
            && TMath::Abs( fIPR_save_mincharge - mincharge ) < 1.e-3 )
    {
        valDT = fIPR_save_dT_from_probCurve;
    }
    else
    {
        // get expected NSB frequency for this charge
        float valIPR = iIPR->Eval( mincharge, 0, "" );
        if( valIPR < 100. || mincharge > iIPR_max )
        {
            valIPR = 100.;   // Hz
        }
        valDT = fProbCurve->Eval( valIPR );
        fIPR_save_dT_from_probCurve = valDT;
        fIPR_save_mincharge = mincharge;
        fIPR_save_telid = fData->getTelID();
        fIPR_save_ProbCurve_par1 = fProbCurve->GetParameter( 1 );
        fIPR_save_ProbCurve_par2 = fProbCurve->GetParameter( 2 );
    }
    
    // apply cut in deltaT
    // (note cut on maximum coincidence limit (given in the cleaning parameter file))
    if( !bInvert )
    {
        if( dT > iCoincWinLimit || dT > valDT )
        {
            return false;
        }
    }
    else
    {
        if( dT < iCoincWinLimit && dT < valDT )
        {
            return true;
        }
    }
    
    return !bInvert;
}

/*
 *
 * relaxed NN group search
 *
 * used for NN 3 and 4
 */
unsigned int VImageCleaning::NNGroupSearchProbCurveRelaxed( unsigned int teltype, TF1* fProbCurve, float PreCut )
{
    if( !fProbCurve )
    {
        return 0;
    }
    
    // Nfold (e.g. 2 or 3)
    int NN = ( int )fProbCurve->GetParameter( 1 );
    
    TGraph* iIPR = ( TGraph* )fIPRgraphs->At( teltype );
    if( !iIPR )
    {
        cout << "VImageCleaning::NNGroupSearchProbCurveRelaxed error, no IPR graph for telescope type " << teltype << endl;
        return 0;
    }
    float iIPR_max = fIPRgraphs_xmax[teltype];
    fIPR_save_mincharge = -99.;
    fIPR_save_dT_from_probCurve = -99.;
    
    int NNcnt = 1;
    float dT = 0.;
    float dt2 = 0.;
    float dt3 = 0.;
    
    ///////////////////////////////////
    // loop over all pixels
    int numpix = fData->getDetectorGeo()->getNumChannels();
    for( int PixNum = 0; PixNum < numpix; PixNum++ )
    {
        int nng3[3];
        
        if( VALIDITY[PixNum] < 0.5 || INTENSITY[PixNum] < PreCut )
        {
            continue;
        }
        NNcnt = 1;
        int pix1 = 0, pix2 = 0, pix3 = 0, pix4 = 0;
        if( PixNum >= ( int )fData->getDetectorGeo()->getNeighbours().size() )
        {
            continue;
        }
        for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[PixNum].size(); j++ )
        {
            Int_t PixNum2 = fData->getDetectorGeo()->getNeighbours()[PixNum][j];
            if( PixNum2 < 0 )
            {
                continue;
            }
            if( VALIDITY[PixNum2] < 0.5 || INTENSITY[PixNum2] < PreCut )
            {
                continue;
            }
            dT = fabs( TIMES[PixNum] - TIMES[PixNum2] );
            float charges[2] = {INTENSITY[PixNum], INTENSITY[PixNum2]};
            float mincharge = 0;
            LocMin( 2, charges, mincharge );
            
            // apply charge and time cut
            if( !NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, dT, CoincWinLimit ) )
            {
                continue;
            }
            
            //////////////////////////////////////////
            float maxtime = 1E6;
            if( NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, dT, 1.e6, true )
                    && VALIDITY[PixNum2] > 0.5 && INTENSITY[PixNum2] > PreCut )
            {
                pix1 = PixNum;
                if( PixNum2 != pix2 )
                {
                    pix3 = PixNum2;
                }
                if( pix2 == 0 )
                {
                    pix2 = PixNum2;
                }
                dt2 = fabs( TIMES[pix2] - TIMES[pix3] );
                dt3 = fabs( TIMES[pix1] - TIMES[pix3] );
                float charges2[3] = {mincharge, INTENSITY[pix3], INTENSITY[pix2]};
                LocMin( 3, charges2, mincharge );
                float times2[3] = { dT, dt2, dt3};
                LocMax( 3, times2, maxtime );
                
                if( NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, mincharge, maxtime, 1.e6, true ) )
                {
                    NNcnt++;
                }
            }
            
            if( NNcnt > 2 )
            {
                if( VALIDITYBUF[pix1] < 4.9 )
                {
                    VALIDITYBUF[pix1] = 5;
                }
                nng3[0] = pix1;
                if( VALIDITYBUF[pix2] < 4.9 )
                {
                    VALIDITYBUF[pix2] = 5;
                }
                nng3[1] = pix2;
                if( VALIDITYBUF[pix3] < 4.9 )
                {
                    VALIDITYBUF[pix3] = 5;
                }
                nng3[2] = pix3;
                short VALIDITYLOCAL[numpix];
                memset( VALIDITYLOCAL, 0, sizeof( VALIDITYLOCAL ) );
                for( int n = 0; n < 3; n++ )
                {
                    VALIDITYLOCAL[nng3[n]] = 10;
                }
                if( NN == 3 )
                {
                    break;
                }
                //4 connected pixels
                for( int n = 0; n < 3; n++ )
                {
                    if( nng3[n] >= ( int )fData->getDetectorGeo()->getNeighbours().size() )
                    {
                        continue;
                    }
                    for( unsigned int jj = 0; jj < fData->getDetectorGeo()->getNeighbours()[nng3[n]].size(); jj++ )
                    {
                        const Int_t testpixnum = fData->getDetectorGeo()->getNeighbours()[nng3[n]][jj];
                        if( testpixnum < 0 || VALIDITYLOCAL[testpixnum] == 10 )
                        {
                            continue;
                        }
                        if( INTENSITY[testpixnum] < PreCut )
                        {
                            continue;
                        }
                        
                        float minchargeloc = 0.;
                        float maxtimeloc = 0.;
                        float charges3[2] = {mincharge, INTENSITY[testpixnum]};
                        LocMin( 2, charges3, minchargeloc );
                        
                        float times3[4] = { maxtime, fabs( TIMES[testpixnum] - TIMES[nng3[0]] ),
                                            fabs( TIMES[testpixnum] - TIMES[nng3[1]] ),
                                            fabs( TIMES[testpixnum] - TIMES[nng3[2]] )
                                          };
                        LocMax( 4, times3, maxtimeloc );
                        
                        if( NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurve, minchargeloc, maxtimeloc, CoincWinLimit, true ) )
                        {
                            NNcnt++;
                            pix4 = testpixnum;
                            float groupsize = INTENSITY[pix1] + INTENSITY[pix2] + INTENSITY[pix3] + INTENSITY[pix4];
                            if( groupsize > PreCut * 4.*1.5 )
                            {
                                if( VALIDITYBUF[pix1] < 5.9 )
                                {
                                    VALIDITYBUF[pix1] = 6;
                                }
                                if( VALIDITYBUF[pix2] < 5.9 )
                                {
                                    VALIDITYBUF[pix2] = 6;
                                }
                                if( VALIDITYBUF[pix3] < 5.9 )
                                {
                                    VALIDITYBUF[pix3] = 6;
                                }
                                if( VALIDITYBUF[pix4] < 5.9 )
                                {
                                    VALIDITYBUF[pix4] = 6;
                                }
                            }
                        }
                    }
                }
            }// end if NNcnt>2
        }
        
    }//end of for() loop
    
    //////////////////////
    // counting
    float NN3 = 0.;
    float NN4 = 0.;
    for( Int_t i = 0; i < numpix; i++ )
    {
        if( VALIDITYBUF[i] > 4.5 && VALIDITYBUF[i] < 5.1 )
        {
            NN3 += 1. / 3.;
        }
        else if( VALIDITYBUF[i] > 5.5 && VALIDITYBUF[i] < 6.1 )
        {
            NN4 += 1. / 4.;
        }
        if( VALIDITYBUF[i] < 1.9 )
        {
            VALIDITYBUF[i] = -1;
        }
    }
    if( NN3 > 2 && NN == 3 )
    {
        return ( unsigned int )( NN3 + 0.5 );
    }
    if( NN4 > 3 && NN == 4 )
    {
        return ( unsigned int )( NN4 + 0.5 );
    }
    return 0;
}

/*
 * discard any pixels with no neighbours
 *
 */
void VImageCleaning::DiscardIsolatedPixels()
{
    unsigned int numpix = fData->getDetectorGeo()->getNumChannels();
    unsigned int NumOfNeighbor = 0;
    int PixNum2 = 0;
    
    for( unsigned int PixNum = 0; PixNum < numpix; PixNum++ )
    {
        if( VALIDITY[PixNum] < 1.9 )
        {
            continue;
        }
        NumOfNeighbor = 0;
        if( PixNum >= fData->getDetectorGeo()->getNeighbours().size() )
        {
            continue;
        }
        for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[PixNum].size(); j++ )
        {
            PixNum2 = fData->getDetectorGeo()->getNeighbours()[PixNum][j];
            if( PixNum2 >= 0 && VALIDITY[PixNum2] > 1.9 )
            {
                NumOfNeighbor++;
            }
        }
        if( NumOfNeighbor < 0.5 )
        {
            VALIDITY[PixNum] = 0;
        }
    }
}

void VImageCleaning::DiscardLocalTimeOutlayers( float NNthresh[6] )
{
    unsigned int numpix = fData->getDetectorGeo()->getNumChannels();
    unsigned int nimagepix = 0;
    DiscardIsolatedPixels();
    for( unsigned int pixnum = 0; pixnum < numpix; pixnum++ )
    {
        if( VALIDITY[pixnum] < 1.9 )
        {
            continue;
        }
        nimagepix++;
    }
    // don't do anything for small images
    if( nimagepix <= 4 )
    {
        return;
    }
    // assume all pixels have same radius!
    Double_t diam = 2.*fData->getDetectorGeo()->getTubeRadius()[1];
    if( diam <= 0. )
    {
        return;
    }
    Double_t x = 0.;
    Double_t y = 0.;
    Double_t xx = 0.;
    Double_t yy = 0.;
    unsigned int pixcnt = 0;
    unsigned int pixzerocnt = 0;
    //******************************************************************
    // discard groups with no neighbouring group in the vicinity of 6pixels
    for( unsigned int pixnum = 0; pixnum < numpix; pixnum++ )
    {
        if( VALIDITY[pixnum] < 1.9 )
        {
            continue;
        }
        
        x = fData->getDetectorGeo()->getX()[pixnum] / diam; // coord in pixels units
        y = fData->getDetectorGeo()->getY()[pixnum] / diam; // coord in pixels units
        pixcnt = 0;
        pixzerocnt = 0;
        // loop over vicinity of 2 rings around pixnum
        for( unsigned int pp = 0; pp < numpix; pp++ )
        {
            if( VALIDITY[pp] < 1.9 || pp == pixnum )
            {
                continue;
            }
            xx = x - fData->getDetectorGeo()->getX()[pp] / diam; // coord in pixels units
            yy = y - fData->getDetectorGeo()->getY()[pp] / diam; // coord in pixels units
            Double_t dist = sqrt( xx * xx + yy * yy );
            if( dist > 6.1 )
            {
                continue;
            }
            if( dist < 2.1 )
            {
                pixcnt++;
            }
            else
            {
                pixzerocnt++;
            }
        }
        if( nimagepix > 7 && pixcnt < 7 && pixzerocnt == 0 )
        {
            VALIDITY[pixnum] = 0;
        }
    }
    //******************************************************************
    // loop over accepted pixels
    float SNRlimit = 5.0;
    float radicand = 0.;
    float sigmalimit = 0.11; // sigmaT>sigmalimit:  time clustering due to finite sampling rate
    unsigned int Tcnt = 0;
    float sigmaT = 0.;
    float meanT = 0.;
    for( unsigned int pixnum = 0; pixnum < numpix; pixnum++ )
    {
        if( VALIDITY[pixnum] < 1.9 )
        {
            continue;
        }
        
        Tcnt = 0;
        sigmaT = 0.;
        meanT = 0.;
        x = fData->getDetectorGeo()->getX()[pixnum] / diam; // coord in pixels units
        y = fData->getDetectorGeo()->getY()[pixnum] / diam; // coord in pixels units
        
        // loop over vicinity of 2 rings around pixnum
        for( unsigned int pp = 0; pp < numpix; pp++ )
        {
            if( VALIDITY[pp] < 1.9 || pp == pixnum )
            {
                continue;
            }
            xx = x - fData->getDetectorGeo()->getX()[pp] / diam; // coord in pixels units
            yy = y - fData->getDetectorGeo()->getY()[pp] / diam; // coord in pixels units
            if( xx * xx + yy * yy > 2.1 * 2.1 )
            {
                continue;
            }
            meanT  += TIMES[pp];
            sigmaT += TIMES[pp] * TIMES[pp];
            Tcnt++;
        }
        if( Tcnt > 1 && nimagepix > 4 )
        {
            meanT /= ( float )Tcnt;
            radicand = ( sigmaT - Tcnt * meanT * meanT ) / ( float( Tcnt ) - 1. );
            if( radicand > 0. )
            {
                sigmaT = sqrt( radicand );
            }
            else
            {
                sigmaT = 1E6;
            }
            
            SNRlimit = 5.0;
            if( Tcnt == 2 )
            {
                SNRlimit = 9.;
            }
            if( sigmaT > 0.
                    && fabs( TIMES[pixnum] - meanT ) / sigmaT > SNRlimit && sigmaT > sigmalimit
                    && INTENSITY[pixnum] < 1.5 * NNthresh[3] )
            {
                VALIDITY[pixnum] = 0;
            }
        }
    }
    //******************************************************************
}


/*
 * set ring of pixels around core pixels
 *
 * this defines the search region for boundary pixels
 *
 */
void  VImageCleaning::SetNeighborRings( unsigned short* VALIDITYBOUNDBUF, float* TIMESReSearch, float* REFTHRESH )
{
    unsigned int nfirstringpix = 0;
    unsigned int numpix = fData->getDetectorGeo()->getNumChannels();
    
    //Define search region, driven by found core pixels
    for( unsigned int p = 0; p < numpix; p++ )
    {
        TIMESReSearch[p] = 0.;
        REFTHRESH[p] = 0.;
        VALIDITYBOUNDBUF[p] = 0;
        // set core pixel
        if( VALIDITY[p] > 1.9 )
        {
            VALIDITYBOUNDBUF[p] = 2;
        }
    }
    // nRing is read from the cleaning parameter file
    for( unsigned int iRing = 0; iRing < nRings; iRing++ )
    {
        for( unsigned int idx = 0; idx < numpix; idx++ )
        {
            // skip core pixel
            if( VALIDITYBOUNDBUF[idx] == 2 )
            {
                continue;
            }
            if( ( iRing > 0 ) && ( VALIDITYBOUNDBUF[idx] < iRing + 7 ) && ( VALIDITYBOUNDBUF[idx] > 1.9 ) )
            {
                continue;
            }
            float time = 0.;
            float refthresh = 0.;
            int n = 0;
            if( idx >= fData->getDetectorGeo()->getNeighbours().size() )
            {
                continue;
            }
            for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[idx].size(); j++ )
            {
                int idx2 = fData->getDetectorGeo()->getNeighbours()[idx][j];
                if( idx2 < 0 || VALIDITYBOUNDBUF[idx2] < 1.9 )
                {
                    continue;
                }
                
                if( iRing == 0 )
                {
                    if( VALIDITYBOUNDBUF[idx2] < 1.9 || VALIDITYBOUNDBUF[idx2] == iRing + 7 )
                    {
                        continue;
                    }
                    if( TIMESReSearch[idx2] > 0.01 )
                    {
                        time += TIMESReSearch[idx2];
                        refthresh += INTENSITY[idx2];
                        n++;
                    }
                    else if( TIMES[idx2] > -500. )
                    {
                        time += TIMES[idx2];
                        refthresh += INTENSITY[idx2];
                        n++;
                    }
                }
                else
                {
                    if( VALIDITYBOUNDBUF[idx2] == iRing + 6 )
                    {
                        time += TIMESReSearch[idx2];
                        refthresh += REFTHRESH[idx2];
                        n++;
                    }
                }
            }
            if( iRing > 0 && n > 0.5 )
            {
                TIMESReSearch[idx] = time / float( n );
                REFTHRESH[idx] = refthresh / float( n );
            }
            if( n > 0.5 )
            {
                VALIDITYBOUNDBUF[idx] = iRing + 7;
                if( iRing == 0 )
                {
                    nfirstringpix++;
                }
            }
        }
    } // loop over rings
}

/*
 * optimized time-next-neighbour cleaning (main function)
 *
 *
 *  return number of groups ngroups
 */
int VImageCleaning::ImageCleaningCharge( unsigned int teltype )
{
    unsigned int numpix = fData->getDetectorGeo()->getNumChannels();
    
    // return value: number of groups
    int ngroups = 0;
    
    //////////////////////////////////////////////////////////////////
    // fill pre thresholds
    // (length must match VDST_MAXNNGROUPTYPES)
    //                 [p.e.]
    // (NOTE: replaced by FillPreThresholds() in the next line
    //  (unit then changed to d.c.)
    float PreThresh[6] = { 2.0,   // 4nn
                           3.0,   // 2+1
                           2.8,   // 3nn
                           5.2,   // 2nn
                           1.8,   // Bound.
                           4.0
                         }; // Bound RefCharge
                         
    FillPreThresholds( ( TGraph* )fIPRgraphs->At( teltype ), PreThresh );
    
    memset( VALIDITYBOUND, 0, sizeof( VALIDITYBOUND ) );
    memset( VALIDITY, 0, sizeof( VALIDITY ) );
    memset( VALIDITYBUF, 0, sizeof( VALIDITYBUF ) );
    ResetCombFactors( teltype );
    
    /////////////////////////////////////////////
    // apply pre-thresholds
    for( unsigned int p = 0; p < numpix; p++ )
    {
        if( INTENSITY[p] > PreThresh[4] )
        {
            VALIDITY[p] = 1;
            VALIDITYBUF[p] = 1;
        }
        else
        {
            VALIDITY[p] = -1;
            VALIDITYBUF[p] = -1;
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*************** Image Cleaning ***************************************************************************
    // first search for NN-groups  (CORE pixels, image cleaning trigger) if NN-group is activated in input card
    
    // 2NN
    if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][3] )
    {
        ngroups = NNGroupSearchProbCurve( teltype, ( TF1* )fProb2nnCurves->At( teltype ), PreThresh[3] );
        for( unsigned int p = 0; p < numpix; p++ )
        {
            if( VALIDITYBUF[p] == 2 )
            {
                VALIDITY[p] = 2;
            }
        }
    }
    
    // 2NNplus1
    if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][1] )
    {
        ngroups += NNGroupSearchProbCurve( teltype, ( TF1* )fProb2plus1Curves->At( teltype ), PreThresh[1] );
        for( unsigned int p = 0; p < numpix; p++ )
        {
            if( VALIDITYBUF[p] == 3 || VALIDITYBOUND[p] == 3 )
            {
                VALIDITY[p] = 3;
            }
        }
    }
    
    // 3NN (note: relaxed search)
    if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][2] )
    {
        ngroups += NNGroupSearchProbCurveRelaxed( teltype, ( TF1* )fProb3nnrelCurves->At( teltype ), PreThresh[2] );
        for( unsigned int p = 0; p < numpix; p++ )
        {
            if( VALIDITYBUF[p] == 5 )
            {
                VALIDITY[p] = 5;
            }
        }
    }
    
    // 4NN (note: relaxed search)
    if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][0] )
    {
        ngroups += NNGroupSearchProbCurveRelaxed( teltype, ( TF1* )fProb4nnCurves->At( teltype ), PreThresh[0] );
        for( unsigned int p = 0; p < numpix; p++ )
        {
            if( VALIDITYBUF[p] == 6 )
            {
                VALIDITY[p] = 6;
            }
        }
    }
    // (end of core pixel determination)
    
    ///////////////////////////////////////////////////////////////////////////////
    //*********************************************************************
    // Boundary pixels
    unsigned int ncorepix = 0;
    unsigned int ncore4nnpix = 0;
    unsigned int nboundsearchpix = 0;
    unsigned int nfirstringpix = 0;
    unsigned int nboundary = 0;
    float TIMESReSearch[numpix];
    unsigned short VALIDITYBOUNDBUF[numpix];
    float REFTHRESH[numpix];
    
    //*****************************************************************************************************************
    // boundary search (same as core search above but reduced search area (vicinity of core pixels) )
    // Define search region for boundary
    SetNeighborRings( &VALIDITYBOUNDBUF[0], &TIMESReSearch[0], &REFTHRESH[0] );
    
    // Reset validity buffers (Important)
    for( unsigned int p = 0; p < numpix; p++ )
    {
        if( VALIDITYBOUNDBUF[p] > 1.9 )
        {
            nboundsearchpix++;
            VALIDITYBUF[p] = 1;
            VALIDITYBOUND[p] = 0;
        }
        else
        {
            VALIDITYBUF[p] = 0;
            VALIDITY[p] = 0;
        }
    }
    
    // all found pixels are set also to CORE pixels!!!
    if( ngroups > 0 )
    {
        ScaleCombFactors( teltype, float( nboundsearchpix ) / ( numpix * 1.5 ) );
        if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][3] )
        {
            NNGroupSearchProbCurve( teltype, ( TF1* )fProb2nnCurves->At( teltype ), 0.8 * PreThresh[3] );
            for( unsigned int p = 0; p < numpix; p++ )
            {
                if( VALIDITY[p] > 1.9 )
                {
                    continue;
                }
                if( VALIDITYBUF[p] == 2 )
                {
                    VALIDITY[p] = 2;
                }
            }
        }
        
        if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][1] )
        {
            NNGroupSearchProbCurve( teltype, ( TF1* )fProb2plus1Curves->At( teltype ), 0.8 * PreThresh[1] );
            for( unsigned int p = 0; p < numpix; p++ )
            {
                if( VALIDITY[p] > 1.9 )
                {
                    continue;
                }
                if( VALIDITYBUF[p] == 3 || VALIDITYBOUND[p] == 3 )
                {
                    VALIDITY[p] = 3;
                }
            }
        }
        
        if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][2] )
        {
            NNGroupSearchProbCurveRelaxed( teltype, ( TF1* )fProb3nnrelCurves->At( teltype ), 0.8 * PreThresh[2] );
            for( unsigned int p = 0; p < numpix; p++ )
            {
                if( VALIDITY[p] > 1.9 )
                {
                    continue;
                }
                if( VALIDITYBUF[p] == 5 )
                {
                    VALIDITY[p] = 5;
                }
            }
        }
        
        if( teltype < fifActiveNN.size() &&  fifActiveNN[teltype][0] )
        {
            NNGroupSearchProbCurveRelaxed( teltype, ( TF1* )fProb4nnCurves->At( teltype ), 0.9 * PreThresh[0] );
            for( unsigned int p = 0; p < numpix; p++ )
            {
                if( VALIDITY[p] > 1.9 )
                {
                    continue;
                }
                if( VALIDITYBUF[p] == 6 )
                {
                    VALIDITY[p] = 6;
                }
            }
        }
        
        ResetCombFactors( teltype );
    }
    // update number of core pixels
    for( unsigned int p = 0; p < numpix; p++ )
    {
        if( VALIDITY[p] > 1.9 && VALIDITY[p] < 6.1 )
        {
            ncorepix++;
        }
    }
    
    // set rings of boundaries for newly found core pixels
    SetNeighborRings( &VALIDITYBOUNDBUF[0], &TIMESReSearch[0], &REFTHRESH[0] );
    for( unsigned int p = 0; p < numpix; p++ )
    {
        if( VALIDITYBOUNDBUF[p] == 7 )
        {
            nfirstringpix++;
        }
    }
    
    
    // BOUNDARY pixel search (usually very few pixels are found)
    // only first ring
    TF1* fProbCurveBound = ( TF1* )fProbBoundCurves->At( teltype );
    TGraph* iIPR = ( TGraph* )fIPRgraphs->At( teltype );
    if( !iIPR )
    {
        cout << "VImageCleaning::ImageCleaningCharge error, no IPR graph for telescope type " << teltype << endl;
        return 0;
    }
    float iIPR_max = fIPRgraphs_xmax[teltype];
    
    for( Int_t iRing = 0; iRing < 1; iRing++ )
    {
        for( UInt_t idx = 0; idx < numpix; idx++ )
        {
            if( ( VALIDITY[idx] < 5.1 ) && ( VALIDITY[idx] > 1.9 ) )
            {
                continue;
            }
            if( ( iRing > 0 ) && ( VALIDITY[idx] < iRing + 7 ) && ( VALIDITY[idx] > 1.9 ) )
            {
                continue;
            }
            int n = 0;
            float time = 0.;
            float charge = 0.;
            
            if( idx >= fData->getDetectorGeo()->getNeighbours().size() )
            {
                continue;
            }
            for( unsigned int j = 0; j < fData->getDetectorGeo()->getNeighbours()[idx].size(); j++ )
            {
                const Int_t idx2 = fData->getDetectorGeo()->getNeighbours()[idx][j];
                if( idx2 < 0 || VALIDITYBOUNDBUF[idx2] < 1.9 )
                {
                    continue;
                }
                
                if( iRing == 0 )
                {
                    if( VALIDITY[idx2] < 1.9 || VALIDITY[idx2] == iRing + 7 )
                    {
                        continue;
                    }
                    // allow for negative time differences in image cleaning
                    // ignore here as it leads to different (sub-optimal)
                    // cleaning results
                    // if( TIMESReSearch[idx2] > -500. )
                    if( TIMESReSearch[idx2] > 0. )
                    {
                        time += TIMESReSearch[idx2];
                        charge += REFTHRESH[idx2];
                        n++;
                    }
                    else
                    {
                        Float_t t = TIMES[idx2];
                        if( t > 0. )
                        {
                            time += t;
                            charge += INTENSITY[idx2];
                            n++;
                        }
                    }
                }
                if( iRing > 0 )
                {
                    if( VALIDITYBOUNDBUF[idx2] == iRing + 6 )
                    {
                        if( TIMESReSearch[idx2] > -500. )
                        {
                            time += TIMESReSearch[idx2];
                            charge += REFTHRESH[idx2];
                            n++;
                        }
                        else
                        {
                            Float_t t = TIMES[idx2];
                            if( t > 0. )
                            {
                                time += t;
                                charge += INTENSITY[idx2];
                                n++;
                            }
                        }
                    }
                }
            }
            if( n > 0.5 )
            {
                if( INTENSITY[idx] < PreThresh[4] )
                {
                    continue;
                }
                TIMESReSearch[idx] = time / float( n );
                charge /= float( n );
                if( charge > 2.*PreThresh[3] )
                {
                    charge = 2 * PreThresh[3];
                }
                float dT = fabs( TIMES[idx] - TIMESReSearch[idx] );
                float charges[2] = {INTENSITY[idx], ( float )charge };
                float refth = 0.;
                LocMin( 2, charges, refth );
                fProbCurveBound->SetParameter( 2, 2.*nfirstringpix );
                
                Double_t valIPRref = iIPR->Eval( charge, 0, "" );
                if( valIPRref < 100. || charge >= iIPR_max )
                {
                    valIPRref = 100.;
                }
                fProbCurveBound->SetParameter( 1, valIPRref );
                
                if( NNChargeAndTimeCut( iIPR, iIPR_max, fProbCurveBound, refth, dT, 0.6 * CoincWinLimit, true ) )
                {
                    VALIDITY[idx] = iRing + 7;
                }
            }
        }
    } // end of loop over rings
    
    
    ///////////////////////////////////////
    // search for outliers
    if( ncorepix > 4 )
    {
        // order is important
        DiscardLocalTimeOutlayers( PreThresh );
        DiscardIsolatedPixels();
    }
    
    for( unsigned int p = 0; p < numpix; p++ )
    {
        if( VALIDITY[p] > 6.1 )
        {
            nboundary++;
        }
        if( VALIDITY[p] == 6 )
        {
            ncore4nnpix++;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////
    // Discard small images and 4nn with no boundaries:
    // Important protection from fully fake images which can worsen
    // shower direction reconstuction
    // prob of any of fake group is <0.8%
    // after this discarding prob of fake image <0.05% !!!
    
    if( ( ( ncorepix + nboundary ) < 4 ) || ( ncore4nnpix == 4 && nboundary == 0 && ncorepix == 4 ) )
    {
        for( unsigned int p = 0; p < numpix; p++ )
        {
            if( VALIDITY[p] > 1.9 )
            {
                VALIDITY[p] = 1;
            }
        }
        ngroups = 0;
    }
    
    return ngroups;
}

/*********************************************************************************************
 * time-next-neighbour cleaning
 * (Maxim Shayduk)
 *
 * see http://arxiv.org/pdf/1307.4939v1.pdf for the NN image cleaning methodology
 */
void VImageCleaning::cleanNNImageFixed( VImageCleaningRunParameter* iImageCleaningParameters )
{
    if( !iImageCleaningParameters )
    {
        return;
    }
    if( !fData )
    {
        printDataError( "VImageCleaning::cleanNNImageFixed" );
    }
    if( !kInitNNImageCleaning )
    {
        printDataError( "VImageCleaning::cleanNNImageFixed image cleaning not initialized" );
    }
    /////////////////////////////////////////////////////////////
    // reset all valarrays regarding the image cleaning
    fData->setImage( false );
    fData->setBorder( false );
    fData->setBrightNonImage( false );
    fData->setImageBorderNeighbour( false );
    
    /////////////////////////////////////////////////////////////
    // NN cleaning is working with different telescope types than rest of evndisp
    
    // get telescope type
    unsigned int teltype = fData->getTelType_Counter( fData->getTelType( fData->getTelID() ) );
    if( teltype >= VDST_MAXTELTYPES )
    {
        char hname[200];
        sprintf( hname, "VImageCleaning::cleanNNImageFixed, invalid telescope type: %u", teltype );
        printDataError( hname );
    }
    
    /////////////////////////////////////////////////////////////
    // initialize IPRs and Prob. curves once for every telescope type
    
    // fake image probability (a typical value is 0.5%)
    fFakeImageProb = iImageCleaningParameters->fNNOpt_FakeImageProb;
    // Init options for NN cleaning
    nRings = iImageCleaningParameters->fNNOpt_nRings;
    CoincWinLimit = iImageCleaningParameters->fNNOpt_CoincWinLimit;
    NNoptNoTimeing = iImageCleaningParameters->fNNOpt_ifNNoptNoTimeing;
    setExplicitSampleTimeSlice = iImageCleaningParameters->fNNOpt_ifExplicitSampleTimeSlice;
    sampleTimeSlice = iImageCleaningParameters->fNNOpt_sampleTimeSlice;
    nBinsADC = iImageCleaningParameters->fNNOpt_nBinsADC;
    // init IPR graphs and rate countours
    if( !kInitNNImgClnPerTelType[teltype] )
    {
        // active multiplicities;
        if( teltype < fifActiveNN.size() )
        {
            fifActiveNN[teltype] = iImageCleaningParameters->fNNOpt_ActiveNN;
        }
        // initiate probability contours
        kInitNNImgClnPerTelType[teltype] = InitNNImgClnPerTelType( teltype );
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // timing parameters for image cleaning
    float FADCslice = fData->getDetectorGeo()->getLengthOfSampleTimeSlice( fData->getTelID() );
    // correct for upsampling in case digital filter is used
    if( fData->getDigitalFilterMethod() > 0 && fData->getDigitalFilterUpSample() > 0 )
    {
        FADCslice /= fData->getDigitalFilterUpSample();
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // fill intensities and times per pixel
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        INTENSITY[i] = 0.;
        TIMES[i]     = -500.;
        VALIDITY[i]  = 0;
        if( fData->getDetectorGeo()->getAnaPixel()[i] > 0
                && !fData->getDead( fData->getHiLo()[i] )[i] )
        {
            INTENSITY[i] = fData->getFADCtoPhe()[i] * fData->getSums()[i];
            //(first option is for tests of impact of timing)
            if( NNoptNoTimeing )
            {
                TIMES[i] = 0.;
            }
            else
            {
                // trace times are corrected for flasher determined time offsets
                TIMES[i] = fData->getPulseTime( true )[i] * FADCslice;
            }
            // make sure that pixels without timing are ignored
            if( TIMES[i] <= -500 && !NNoptNoTimeing )
            {
                INTENSITY[i] = 0;
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // optimized NN image cleaning
    int ngroups = ImageCleaningCharge( teltype );
    
    ///////////////////////////////////////////////////////////////////////////////////
    // set image flags for valid pixels
    // note: no border pixel flags are set
    if( ngroups > 0 )
    {
        // set image/border pixels
        for( unsigned int i = 0; i < fData->getNChannels(); i++ )
        {
            fData->setImage( i, false );
            fData->setBorder( i, false );
            if( fData->getDetectorGeo()->getAnaPixel()[i] > 0
                    && !fData->getDead( fData->getHiLo()[i] )[i] )
            {
                if( VALIDITY[i] > 1.9 )
                {
                    fData->setImage( i, true );
                }
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////
    // (end of optimized NN cleaning)
    
    // recover image pixels besides dead pixels
    recoverImagePixelNearDeadPixel();
    
    // fill list of pixels which are beside image/border pixels
    fillImageBorderNeighbours();
}

/*
 * define rate contour function
 *
 * equation 2 in ICRC 2013 proceedings
 *
 * parameter plane formed by minimal charge and maximal coincidence time delta T
 *
 */
TF1* VImageCleaning::defineRateContourFunction( unsigned int teltype, TString funcname, float iRate, float iNfold, float iCombFactor,
        float xlow, float xup )
{
    // (factor 1.e9: work in [ns])
    TF1* f1 = new TF1( funcname, "1.0E9 * TMath::Exp( 1. / ( [1] - 1 ) * TMath::Log( [0] / ( [2] * pow( x, [1] ) ) ) )", xlow, 3.*xup );
    f1->SetParameters( iRate, iNfold, iCombFactor );
    f1->SetParNames( "Rate", "Nfold", "CombFactor" );
    f1->SetName( ( TString )f1->GetName() + Form( "Type%d", teltype ) );
    
    return f1;
}

TF1* VImageCleaning::defineRateContourBoundFunction( unsigned int teltype, TString funcname, float iRate, float refThresh, float iCombFactor,
        float xlow, float xup )
{
    TF1* f1 = new TF1( funcname, "1.0E9 * TMath::Exp( 1. / ( 2. - 1. ) * TMath::Log( [0] / ( [2] * x * [1] ) ) )", xlow, 3.*xup );
    f1->SetParameters( iRate, refThresh, iCombFactor, xup );
    f1->SetParNames( "Rate", "RefThresh", "CombFactor", "ChargeMax" );
    f1->SetName( ( TString )f1->GetName() + Form( "Type%d", teltype ) );
    
    return f1;
}

/*
 * fill pre thresholds
 *
 * depends on NN configuration
 *
*/
void VImageCleaning::FillPreThresholds( TGraph* gipr, float NNthresh[6] )
{
    if( !gipr )
    {
        cout << "VImageCleaning::FillPreThresholds() error filling pre-thresholds" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    // assuming monotonous ipr curve (as it should be)
    //                    4nn   2+1   3nn   2nn   bound  ref for bound
    float ThreshFreq[6] = {8.5E6, 2.4E6, 4.2E6, 1.5E5, 1.1E7, 4.0E5};
    float gXres = ( gipr->GetXaxis()->GetXmax() - gipr->GetXaxis()->GetXmin() ) / gipr->GetN();
    for( int i = 0; i < 6; i++ )
    {
        for( int t = 0; t < gipr->GetN(); t++ )
        {
            float x = gipr->GetXaxis()->GetXmin() + gXres * ( float( t ) );
            float val = gipr->Eval( x );
            if( val <= ThreshFreq[i] )
            {
                NNthresh[i] = x;
                break;
            }
        }
    }
    if( NNthresh[4] < 0.001 )
    {
        NNthresh[4] = NNthresh[3] / 3.4;   //pre-threshold for boundary, if not set by IPR curve scanning due to IPR curve behaviour for small thresholds
    }
    if( NNthresh[0] < 0.001 )
    {
        NNthresh[0] = NNthresh[3] / 2.75;   //pre-threshold for boundary, if not set by IPR curve scanning
    }
}

void VImageCleaning::FillIPR( unsigned int teltype ) //tel type
{
    float  gIPRUp = 1500.; //charge in FADC counts
    float  gIPRLo = 0.;
    float  gIPRStep = ( gIPRUp - gIPRLo ) / float( fIPRdim );
    float RATE[fIPRdim];
    for( unsigned int thbin = 0; thbin < fIPRdim; thbin++ )
    {
        RATE[thbin] = 0;
    }
    
    //loop over pixels
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getDetectorGeo()->getAnaPixel()[i] > 0 && !fData->getDead( fData->getHiLo()[i] )[i] )
        {
            if( VALIDITY[i] < 1.9 )
            {
                // event counter
                IPR[teltype][0] += 1;
                for( unsigned int thbin = 1; thbin < fIPRdim; thbin++ )
                {
                    float val = gIPRLo + float( thbin ) * gIPRStep;
                    //if(INTENSITY[i]>val) RATE[thbin]+=1.;
                    IPR[0][thbin] = val;
                    if( fData->getSums()[i] > val )
                    {
                        RATE[thbin] += 1.;
                    }
                    if( fData->getSums()[i] > val )
                    {
                        IPR[teltype][thbin] += 1.;
                    }
                }
            }
        }
    }
}

TGraphErrors* VImageCleaning::GetIPRGraph( unsigned int teltype, float ScanWindow )
{
    float RATES[fIPRdim],  RATESERR[fIPRdim];
    float RATESX[fIPRdim], RATESXERR[fIPRdim];
    
    for( unsigned int i = 0; i < fIPRdim; i++ )
    {
        RATES[i] = IPR[teltype][i];
        float ConvToHz = ( ScanWindow * IPR[teltype][0] ) / 1E9;
        if( ConvToHz < 0.9E-9 )
        {
            break;
        }
        RATES[i] /= ConvToHz;
        RATESX[i] = IPR[0][i];
        RATESERR[i] = sqrt( RATES[i] * ConvToHz ) / ConvToHz;
        //RATESXERR[i]=0.1*float(i)/float(res);
        RATESXERR[i] = 0.;
    }
    TGraphErrors* gRate = new TGraphErrors( fIPRdim, RATESX, RATES, RATESXERR, RATESERR );
    gRate->SetTitle( "IPRcharge" );
    gRate->SetName( "IPRcharge" );
    gRate->GetXaxis()->SetTitle( "Threshold, FADC" );
    gRate->GetYaxis()->SetTitle( "Rate, Hz" );
    gRate->SetMinimum( 1 );
    
    TFile* fgraph = new TFile( "$CTA_USER_DATA_DIR/TestIPRtrunk.root", "RECREATE" );
    gRate->Write();
    fgraph->Close();
    std::cout << "[TriggerAnalogueSummation::GetIPRGraph()]: graph root file written:" << fgraph->GetName() << std::endl;
    return gRate;
}


// end of NN image cleaning
//*****************************************************************************************************

/*
   simple cluster cleaning.
   get image pixels (above threshold)
   find clusters of touching pixels
   keep N biggest clusters; cut on cluster size and number of pixels
   add border pixels
   in part derived from time cluster cleaning
*/
void VImageCleaning::cleanImageWithClusters( VImageCleaningRunParameter* iImageCleaningParameters, bool isFixed )
{

    ///////////////////////////////////////////////
    // setting of image cleaning run parameters
    if( !iImageCleaningParameters )
    {
        return;
    }
    double hithresh = iImageCleaningParameters->fimagethresh;
    double lothresh = iImageCleaningParameters->fborderthresh;
    double brightthresh = iImageCleaningParameters->fbrightnonimagetresh;
    int minNumPixel = iImageCleaningParameters->fminpixelcluster;
    int maxNumCluster = iImageCleaningParameters->fnmaxcluster;
    double minSizeCluster = iImageCleaningParameters->fminsizecluster;
    
    
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::cleanImageWithClusters " << fData->getTelID() << "\t" << brightthresh << endl;
    }
    
    // STEP 1: Select all pixels with a signal > hithresh
    fData->setImage( false );
    fData->setBorder( false );
    fData->setBrightNonImage( false );
    fData->setImageBorderNeighbour( false );
    double i_pedvars_i = 0.;
    
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getDetectorGeo()->getAnaPixel()[i] < 1 || fData->getDead( i, fData->getHiLo()[i] ) )
        {
            continue;
        }
        i_pedvars_i = fData->getPedvars( fData->getCurrentSumWindow()[i], fData->getHiLo()[i] )[i];
        
        if( ( isFixed && fData->getSums()[i] > hithresh ) || ( !isFixed && fData->getSums()[i] > hithresh * i_pedvars_i ) )
        {
            fData->setImage( i, true );
            fData->setBorder( i, false );
        }
        if( ( isFixed && fData->getSums()[i] > brightthresh ) || ( !isFixed && fData->getSums()[i] > brightthresh * i_pedvars_i ) )
        {
            fData->setBrightNonImage( i, true );
        }
    }
    
    // STEP 2: Make clusters
    // - group touching core pixels to clusters
    
    //in case there is still something in the memory, reset all channels to cluster id 0
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        fData->setClusterID( i, 0 );
    }
    
    unsigned int fNCluster = 0;
    fSizeCluster.clear();
    fNpixCluster.clear();
    
    //cluster 0 is for non-image pixels & non-analyzed pixels.
    fSizeCluster.push_back( 0 );
    fNpixCluster.push_back( 0 );
    
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getImage()[i] && fData->getClusterID()[i] == 0 ) //new cluster
        {
            fNCluster++;
            fSizeCluster.push_back( 0 );
            fNpixCluster.push_back( 0 );
            //recursively add channel & its neighbors to the cluster.
            addToCluster( fNCluster, i );
        }
    }
    
    // find clusters passing the size cuts
    vector<int> good_clusters;
    for( unsigned int i = 1; i <= fNCluster; i++ )
    {
        if( fSizeCluster.at( i ) >= minSizeCluster && fNpixCluster.at( i ) >= minNumPixel )
        {
            good_clusters.push_back( i );
        }
        else
        {
            removeCluster( i );
        }
    }
    
    //sort clusters by size. Bubblesort algorithm from wikipedia.
    unsigned int n = good_clusters.size();
    do
    {
        int newn = 0;
        for( unsigned int i = 1; i < n; i++ )
        {
            if( fSizeCluster.at( good_clusters.at( i - 1 ) ) < fSizeCluster.at( good_clusters.at( i ) ) )
            {
                unsigned int temp = good_clusters.at( i );
                good_clusters.at( i ) = good_clusters.at( i - 1 );
                good_clusters.at( i - 1 ) = temp;
                newn = i;
            }
        }
        n = newn;
    }
    while( n > 0 );
    
    //only keep the nmax largest clusters
    if( maxNumCluster > 0 )
    {
        for( unsigned int i = maxNumCluster; i < good_clusters.size() ; i++ )
        {
            removeCluster( good_clusters.at( i ) );
            good_clusters.erase( good_clusters.begin() + i );
            i--;
        }
    }
    
    //add border pixels
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        int i_ID = fData->getClusterID()[i];
        if( fData->getImage()[i] && fData->getClusterID()[i] > 0 )
        {
            for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[i]; j++ )
            {
                unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( fData->getImage()[k] || fData->getBorder()[k] )
                {
                    continue;
                }
                if( ( isFixed && fData->getSums()[k] > lothresh ) || ( !isFixed && fData->getSums()[k] > lothresh * fData->getPedvars( fData->getCurrentSumWindow()[k], fData->getHiLo()[k] )[k] ) )
                {
                    fData->setBorder( k, true );
                    fData->setClusterID( k, i_ID );
                }
            }
        }
    }
    
    // (preli) set the trigger vector in MC case (preli)
    // trigger vector are image/border tubes
    if( fData->getReader() )
    {
        if( fData->getReader()->getDataFormatNum() == 1 || fData->getReader()->getDataFormatNum() == 4	|| fData->getReader()->getDataFormatNum() == 6 )
        {
            fData->getReader()->setTrigger( fData->getImage(), fData->getBorder() );
        }
    }
    // (end of preli)
    
    recoverImagePixelNearDeadPixel();
    fillImageBorderNeighbours();
    
    return;
    
}

//recursively add a pixel and its neigboring image pixels to a cluster.
void VImageCleaning::addToCluster( unsigned int cID, unsigned int iChan )
{
    fData->setClusterID( iChan, cID );
    if( fSizeCluster.size() > ( unsigned int ) cID )
    {
        fSizeCluster.at( cID ) += fData->getSums()[iChan];
    }
    else
    {
        cout << "VImageCleaning::addToCluster warning: ClusterID " << cID << " not available in fSizeCluster vector (size " <<  fSizeCluster.size() << ")" << endl;
    }
    if( fNpixCluster.size() > ( unsigned int ) cID )
    {
        fNpixCluster.at( cID )++;
    }
    else
    {
        cout << "VImageCleaning::addToCluster warning: ClusterID " << cID << " not available in fNpixCluster vector (size " <<  fNpixCluster.size() << ")" << endl;
    }
    
    for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[iChan]; j++ )
    {
        unsigned int k = fData->getDetectorGeo()->getNeighbours()[iChan][j];
        if( fData->getImage()[k] && fData->getClusterID()[k] == 0 )
        {
            addToCluster( cID, k ) ;
        }
    }
    return;
}

//remove a cluster (ie. set all its pixels to non-image status)
void VImageCleaning::removeCluster( unsigned int cID )
{
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getClusterID()[i] == ( int )cID )
        {
            fData->setImage( i, false );
        }
    }
    return;
}


/*!
  Image cleaning routine using pixel timing information
  based on Nepomuk's PhD thesis time-cluster cleaning algorithm
   - uses fixed time differences for discrimination of pixels/clusters
   - adjusts time difference according to the time gradient
   - handles single core pixel
   - BrightNonImages not completely implemented yet (needs checks - but who is really using them?)

   \par hithresh image threshold
   \par lothresh border threshold
   \par brightthresh bright pixel threshold
   \par timeCutPixel time difference between pixels
   \par timeCutCluster time difference between clusters
   \par minNumPixel minimum number of pixels in a cluster
   \par loop_max number of loops
*/
void VImageCleaning::cleanImageFixedWithTiming( VImageCleaningRunParameter* iImageCleaningParameters )
{
    cleanImageWithTiming( iImageCleaningParameters, true );
}

void VImageCleaning::cleanImagePedvarsWithTiming( VImageCleaningRunParameter* iImageCleaningParameters )
{
    cleanImageWithTiming( iImageCleaningParameters, false );
}


void VImageCleaning::cleanImageWithTiming( VImageCleaningRunParameter* iImageCleaningParameters, bool isFixed )
{
    ///////////////////////////////////////////////
    // setting of image cleaning run parameters
    if( !iImageCleaningParameters )
    {
        return;
    }
    double hithresh     = iImageCleaningParameters->fimagethresh;
    double lothresh     = iImageCleaningParameters->fborderthresh;
    double brightthresh = iImageCleaningParameters->fbrightnonimagetresh;
    double timeCutPixel = iImageCleaningParameters->ftimecutpixel;
    double timeCutCluster = iImageCleaningParameters->ftimecutcluster;
    int minNumPixel     = iImageCleaningParameters->fminpixelcluster;
    int loop_max        = iImageCleaningParameters->floops;
    
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::cleanImageWithTiming " << fData->getTelID() << "\t" << brightthresh << endl;
    }
    
    ///////////////////////////////////////////////
    // check if time gradient was already calculated
    if( fData->getImageParameters()->tgrad_x == 0 && fData->getImageParameters()->tint_x == 0 )
    {
        timeCutPixel = 5.0;
        timeCutCluster = fData->getNSamples();         // = number of readout samples
    }
    else
    {
        float tpix = fData->getImageParameters()->tgrad_x * 0.15; // VERITAS pixel size (HARDCODED)
        if( tpix > timeCutPixel )
        {
            timeCutPixel = tpix + 0.1;                     // added some uncertainties adhoc
        }
    }
    
    fData->setImage( false );
    fData->setBorder( false );
    fData->setBrightNonImage( false );
    fData->setImageBorderNeighbour( false );
    unsigned int i_nchannel = fData->getNChannels();
    
    //////////////////////////////////////////////////////
    // STEP 1: Select all pixels with a signal > hithresh
    //
    
    double i_pedvars_i = 0.;
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        // select for dead channels
        // note: ignore all channels with highgain setting to 'dead'
        if( fData->getDetectorGeo()->getAnaPixel()[i] < 1
                || fData->getDead( fData->getHiLo()[i] )[i]
                || fData->getDead( false )[i] )
        {
            continue;
        }
        
        if( isFixed )
        {
            if( fData->getSums()[i] > hithresh )
            {
                fData->setImage( i, true );
            }
        }
        else
        {
            i_pedvars_i = fData->getPedvars( fData->getCurrentSumWindow()[i], fData->getHiLo()[i] )[i];
            if( fData->getSums()[i] > hithresh * i_pedvars_i )
            {
                fData->setImage( i, true );
            }
        }
    }
    
    //////////////////////////////////////////////////////////
    // STEP 2: Make clusters
    //
    //       - group touching core pixels to clusters
    //         (only if the pixels are coincident in time!)
    
    
    // REALLY NEEDED: in case there is still something in the memory, reset all channels to cluster id 0
    //                needs to be checked for memory leaks!!!
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        fData->setClusterID( i, 0 );
    }
    
    int i_cluster = 0;
    int c_id = 0;
    
    for( unsigned int i = 0; i < i_nchannel; i++ )
    {
        if( fData->getImage()[i] )
        {
            if( fData->getClusterID()[i] != 0 )
            {
                c_id = fData->getClusterID()[i];
            }
            else if( i_cluster == 0 || fData->getClusterID()[i] == 0 )
            {
                i_cluster++;
                c_id = i_cluster;
            }
            else
                cout << "WARNING: Something looks wrong - this should not happen\n"
                     << "[" << i_cluster << "/" << c_id << "] " << i << endl;
                     
            fData->setClusterID( i, c_id );
            
            if( i >= fData->getDetectorGeo()->getNeighbours().size() )
            {
                continue;
            }
            for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[i]; j++ )
            {
                unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( fData->getImage()[k] && fData->getClusterID()[k] == 0 )
                {
                    if( fabs( fData->getPulseTime()[i] - fData->getPulseTime()[k] ) < timeCutPixel )
                    {
                        fData->setClusterID( k, c_id );
                    }
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////
    // STEP 3: Calculate for each cluster weighted mean time
    //
    //   - each cluster gets its own "clustertime" and "clustersize"
    //   - the cluster with the largest size is the "main" cluster
    //
    // new: for time difference adjustment calculate "cenx" and "ceny"
    
    int i_clusterNpix;      // number of pixels in cluster
    double i_clustersize;   // size of cluster
    double i_clustertime;   // weighted mean time of all pixels in a cluster
    
    double i_cenx, i_ceny;
    
    double i_mainclustersize = 0; // size of the "main cluster"
    
    int cluster = 0;
    while( cluster <= i_cluster )
    {
    
        i_clusterNpix = 0;
        i_clustersize = 0.;
        i_clustertime = 0.;
        double i_clustercenx = 0.;
        i_cenx = 0.;
        double i_clusterceny = 0.;
        i_ceny = 0.;
        
        for( unsigned int i = 0; i < fData->getNChannels(); i++ )
        {
            if( fData->getClusterID()[i] == ( unsigned int ) cluster && fData->getImage()[i] )
            {
                i_clusterNpix++;
                
                i_clustersize += fData->getSums()[i];
                i_clustertime += ( fData->getSums()[i] * fData->getPulseTime()[i] );
                
                double xi = fData->getDetectorGeo()->getX()[i];
                double yi = fData->getDetectorGeo()->getY()[i];
                
                i_cenx += ( fData->getSums()[i] * xi );
                i_ceny += ( fData->getSums()[i] * yi );
            }
        }
        if( i_clustersize != 0 )
        {
            i_clustertime = i_clustertime / i_clustersize;
            i_clustercenx = i_cenx / i_clustersize;
            i_clusterceny = i_ceny / i_clustersize;
        }
        else
        {
            i_clustertime = -99;
            i_clustersize = -99;
            i_clustercenx = -99;
            i_clusterceny = -99;
        }
        
        fData->setClusterNpix( cluster, i_clusterNpix );
        fData->setClusterSize( cluster, i_clustersize );
        fData->setClusterTime( cluster, i_clustertime );
        fData->setClusterCenx( cluster, i_clustercenx );
        fData->setClusterCeny( cluster, i_clusterceny );
        
        if( i_clustersize == -99 )
        {
            fData->setMainClusterID( 0 );
        }
        else if( i_clustersize >= i_mainclustersize )
        {
            i_mainclustersize = i_clustersize;
            fData->setMainClusterID( cluster );
        }
        cluster++;
    }
    //     cout << "MAIN CLUSTER " << getMainClusterID()
    // 	 << ": Npix=" << getClusterNpix()[ getMainClusterID() ]
    // 	 << " Size=" << getClusterSize()[ getMainClusterID() ]
    // 	 << " Time=" << getClusterTime()[ getMainClusterID() ] << endl;
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // STEP 4: eliminate all clusters with time differences > Tcluster to the main cluster
    //
    //         NEW: use tgrad_x to calculate Tcluster
    
    
    int i_ID = 0;
    
    float i_clusterX, i_clusterY;  // c.o.g. (cluster)
    float i_clusterXpos;           // position on major axis
    float i_clusterXtime;          // time on major axis
    
    if( fData->getImageParameters()->tgrad_x != 0 && fData->getImageParameters()->tint_x != 0 )
    {
        // c.o.g. (main cluster)
        float i_mainX = fData->getClusterCenx()[ fData->getMainClusterID() ];
        // c.o.g. (main cluster)
        float i_mainY = fData->getClusterCeny()[ fData->getMainClusterID() ];
        // position on major axis
        float i_mainXpos = i_mainX * fData->getImageParameters()->cosphi + i_mainY * fData->getImageParameters()->sinphi;
        // time on major axis
        float i_mainXtime = fData->getImageParameters()->tgrad_x * i_mainXpos;
        
        i_clusterX = 0.;
        i_clusterY = 0.;
        i_clusterXpos = 0.;
        i_clusterXtime = 0.;
        
        cluster = 1;
        while( cluster <= i_cluster )
        {
            i_clusterX = fData->getClusterCenx()[ cluster ];
            i_clusterY = fData->getClusterCeny()[ cluster ];
            i_clusterXpos = i_clusterX * fData->getImageParameters()->cosphi + i_clusterY * fData->getImageParameters()->sinphi;
            i_clusterXtime = fData->getImageParameters()->tgrad_x * i_clusterXpos;
            
            for( unsigned int i = 0; i < fData->getNChannels(); i++ )
            {
                i_ID = fData->getClusterID()[i];
                if( i_ID == 0 || i_ID == fData->getMainClusterID() )
                {
                    continue;
                }
                
                if( i_ID == cluster && fData->getImage()[i] )
                {
                    if( fabs( ( fData->getClusterTime()[i_ID] - fData->getClusterTime()[fData->getMainClusterID()] ) - ( i_clusterXtime - i_mainXtime ) ) > timeCutCluster )
                    {
                        fData->setImage( i, false );
                        fData->setClusterID( i, -99 );
                    }
                }
            }
            cluster++;
        }
    }
    else
    {
        for( unsigned int i = 0; i < fData->getNChannels(); i++ )
        {
            if( fData->getClusterID()[i] == 0 )
            {
                continue;
            }
            
            if( fData->getImage()[i] )
            {
                i_ID = fData->getClusterID()[i];
                i_pedvars_i = fData->getPedvars( fData->getCurrentSumWindow()[i], fData->getHiLo()[i] )[i];
                
                if( i_ID != fData->getMainClusterID() && fabs( fData->getClusterTime()[i_ID] - fData->getClusterTime()[fData->getMainClusterID()] ) > timeCutCluster )
                {
                    fData->setImage( i, false );
                    fData->setClusterID( i, -99 );
                }
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////
    // STEP 5: find boundary pixels (and add them to the cluster)
    //
    //  selection criteria:
    //  - signal above lothresh
    //  - loop through already cleaned direct neighbors with time difference < Tpixel
    //
    //  this step can be reiterated more than once (number of loops)
    //
    
    for( int loop = 0; loop < loop_max; loop++ )
    {
        int counter = 0;
        int tmp_border[fData->getNChannels()];
        int tmp_cluster[fData->getNChannels()];
        
        double i_pedvars_k = 0.;
        for( unsigned int i = 0; i < fData->getNChannels(); i++ )
        {
            i_ID = fData->getClusterID()[i];
            if( fData->getImage()[i] || fData->getBorder()[i] )
            {
                for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[i]; j++ )
                {
                    unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                    
                    if( isFixed )
                    {
                        if( !fData->getImage()[k] && !fData->getBorder()[k] && fData->getSums()[k] > lothresh
                                && fabs( fData->getPulseTime()[i] - fData->getPulseTime()[k] ) < timeCutPixel )
                        {
                            tmp_border[counter] = k;
                            tmp_cluster[counter] = i_ID;
                            counter++;
                        }
                    }
                    else
                    {
                        i_pedvars_k = fData->getPedvars( fData->getCurrentSumWindow()[k], fData->getHiLo()[k] )[k];
                        if( !fData->getImage()[k] && !fData->getBorder()[k]
                                && fData->getSums()[k] > lothresh * i_pedvars_k
                                && fabs( fData->getPulseTime()[i] - fData->getPulseTime()[k] ) < timeCutPixel )
                        {
                            tmp_border[counter] = k;
                            tmp_cluster[counter] = i_ID;
                            counter++;
                        }
                    }
                }
            }
        }
        if( fData->getMainClusterID() != 0 )
        {
            for( int pixel = 0; pixel < counter; pixel++ )
            {
                fData->setBorder( tmp_border[pixel], true );
                fData->setClusterID( tmp_border[pixel], tmp_cluster[pixel] );
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////
    // STEP 6: merge touching clusters and reset cluster parameters
    //
    
    mergeClusters();
    
    for( int x = 0; x <= i_cluster; x++ )
    {
        fData->setClusterNpix( x, 0 );
        fData->setClusterSize( x, 0 );
        fData->setClusterTime( x, 0 );
    }
    
    i_mainclustersize = 0;
    cluster = 1;
    
    while( cluster <= i_cluster )
    {
    
        i_clusterNpix = 0;
        i_clustersize = 0.;
        i_clustertime = 0.;
        
        for( unsigned int i = 0; i < fData->getNChannels(); i++ )
        {
            if( fData->getClusterID()[i] == ( unsigned int ) cluster )
            {
                if( fData->getImage()[i] || fData->getBorder()[i] )
                {
                    i_clusterNpix++;
                    i_clustersize += fData->getSums()[i];
                    i_clustertime += fData->getSums()[i] * fData->getPulseTime()[i];
                }
                else
                {
                    fData->setClusterID( i, -99 );
                }
            }
        }
        if( i_clustersize != 0 )
        {
            i_clustertime = i_clustertime / i_clustersize;
        }
        else
        {
            i_clustertime = -99.;
            i_clustersize = -99;
        }
        
        fData->setClusterNpix( cluster, i_clusterNpix );
        fData->setClusterSize( cluster, i_clustersize );
        fData->setClusterTime( cluster, i_clustertime );
        
        if( i_clustersize == -99 )
        {
            fData->setMainClusterID( 0 );
        }
        else if( i_clustersize >= i_mainclustersize )
        {
            i_mainclustersize = i_clustersize;
            fData->setMainClusterID( cluster );
        }
        cluster++;
    }
    
    //  count number of clusters before removing (useful for hadron rejection?)
    set< int > tmp_counter_uncleaned;
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getImage()[i] || fData->getBorder()[i] )
        {
            i_ID = fData->getClusterID()[i];
            if( i_ID != 0 && i_ID != -99 )
            {
                tmp_counter_uncleaned.insert( i_ID );
            }
        }
    }
    fData->setNcluster_uncleaned( tmp_counter_uncleaned.size() );
    
    /////////////////////////////////////////////////////////////
    // STEP 7: eliminate clusters with less then XX pixels
    //         & clusters where one core pixel has less then 2 direct border pixels
    
    removeSmallClusters( minNumPixel );
    
    
    /////////////////////////////////////////////////////////////
    // FINAL STEP:
    //
    //  - count number of clusters (useful for hadron rejection?)
    
    set< int > tmp_counter_cleaned;
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getImage()[i] || fData->getBorder()[i] )
        {
            i_ID = fData->getClusterID()[i];
            if( i_ID != 0 && i_ID != -99 )
            {
                tmp_counter_cleaned.insert( i_ID );
            }
        }
    }
    fData->setNcluster_cleaned( tmp_counter_cleaned.size() );
    
    if( fData->getReader()->getDataFormatNum() == 1 || fData->getReader()->getDataFormatNum() == 4 || fData->getReader()->getDataFormatNum() == 6 )
    {
        fData->getReader()->setTrigger( fData->getImage(), fData->getBorder() );
    }
    
    
    fillImageBorderNeighbours();
}


void VImageCleaning::mergeClusters()
{
    int i_clusterID;
    int k_clusterID;
    
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getImage()[i] || fData->getBorder()[i] )
        {
            i_clusterID = fData->getClusterID()[i];
            unsigned int i_neighbour_size = fData->getDetectorGeo()->getNNeighbours()[i];
            
            for( unsigned int j = 0; j < i_neighbour_size; j++ )
            {
                unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                k_clusterID = fData->getClusterID()[k];
                
                if( ( fData->getImage()[k] || fData->getBorder()[k] ) && k_clusterID != 0 && ( unsigned int ) k_clusterID != fData->getClusterID()[i] )
                {
                    for( unsigned int n = 0; n < fData->getNChannels(); n++ )
                    {
                        if( ( fData->getImage()[n] || fData->getBorder()[n] ) && fData->getClusterID()[n] == fData->getClusterID()[k] )
                        {
                            fData->setClusterID( n, i_clusterID );
                        }
                    }
                }
            }
        }
    }
}


void VImageCleaning::removeSmallClusters( int minPix )
{
    int i_cluster = 0;
    
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        i_cluster = fData->getClusterID()[i];
        if( i_cluster == 0 || i_cluster == -99 )
        {
            continue;
        }
        
        // remove clusters with less then minPix
        if( fData->getClusterNpix()[i_cluster] < minPix )
        {
            if( fData->getImage()[i] )
            {
                fData->setImage( i, false );
                fData->setClusterID( i, -99 );
            }
            else if( fData->getBorder()[i] )
            {
                fData->setBorder( i, false );
                fData->setClusterID( i, -99 );
            }
        }
        
        // remove single core pixels with less than two direct border pixels
        int c1 = 0;
        int c2 = 0;
        
        bool dont_remove = false;
        
        if( fData->getImage()[i] )
        {
            unsigned int i_neighbour_size = fData->getDetectorGeo()->getNNeighbours()[i];
            for( unsigned int j = 0; j < i_neighbour_size; j++ )
            {
                unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( fData->getImage()[k] )
                {
                    dont_remove = true;
                }
                else if( fData->getBorder()[k] )
                {
                    c1++;
                }
                else if( k < fData->getDead().size()
                         && ( fData->getDead( fData->getHiLo()[k] )[k] || fData->getDead( false )[k] ) )
                {
                    for( unsigned l = 0; l < i_neighbour_size; l++ )
                    {
                        unsigned int m = fData->getDetectorGeo()->getNeighbours()[i][l];
                        if( m != i && m < fData->getBorder().size() && ( fData->getBorder()[m] || fData->getImage()[m] ) )
                        {
                            c2++;
                        }
                    }
                }
            }
            if( dont_remove )
            {
                continue;
            }
        }
        
        if( c1 + c2 < 2 && fData->getImage()[i] )
        {
            fData->setImage( i, false );
            fData->setBorder( i, false );
            fData->setBrightNonImage( i, true );
            fData->setClusterID( i, -99 );
            
            // remove the rest of the single core cluster (if it exists)
            unsigned int i_neighbour_size = fData->getDetectorGeo()->getNNeighbours()[i];
            for( unsigned int j = 0; j < i_neighbour_size; j++ )
            {
                unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( fData->getBorder()[k] )
                {
                    fData->setBorder( k, false );
                    
                    for( unsigned l = 0; l < fData->getDetectorGeo()->getNNeighbours()[k]; l++ )
                    {
                        unsigned int m = fData->getDetectorGeo()->getNeighbours()[k][l];
                        if( fData->getBorder()[m] )
                        {
                            fData->setBorder( m, false );
                        }
                    }
                }
            }
        }
    }
    
}

/*
 * fill a list with pixels which are neighbouring image or border pixels
 *
 */
void VImageCleaning::fillImageBorderNeighbours()
{
    fData->setImageBorderNeighbour( false );
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getImage()[i] || fData->getBorder()[i] )
        {
            // a pixel is its own neighbour :-)
            fData->getImageBorderNeighbour()[i] = true;
            // loop over all neighbours
            unsigned int i_neighbour_size = fData->getDetectorGeo()->getNNeighbours()[i];
            for( unsigned int j = 0; j < i_neighbour_size; j++ )
            {
                unsigned int k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( k < fData->getImageBorderNeighbour().size() && !fData->getDead()[k] )
                {
                    fData->getImageBorderNeighbour()[k] = true;
                }
            }
        }
    }
}

/*
 * loop again to remove isolated image pixels
 * if neighbour is dead, check neighbours of this dead channel (see e.g. run 329 event 709)
 */
void VImageCleaning::recoverImagePixelNearDeadPixel()
{
    bool i_neigh = false;
    unsigned int i_neighbour_size = 0;
    unsigned int k = 0;
    
    for( unsigned int i = 0; i < fData->getNChannels(); i++ )
    {
        if( fData->getImage()[i] )
        {
            i_neigh = false;
            i_neighbour_size = fData->getDetectorGeo()->getNNeighbours()[i];
            for( unsigned int j = 0; j < i_neighbour_size; j++ )
            {
                k = fData->getDetectorGeo()->getNeighbours()[i][j];
                if( k < fData->getBorder().size() && ( fData->getBorder()[k] || fData->getImage()[k] ) )
                {
                    fData->setImage( i, true );
                    i_neigh = true;
                    break;
                }
                else if( k < fData->getDead().size() && fData->getDead( k, fData->getHiLo()[k] ) )
                {
                    for( unsigned l = 0; l < i_neighbour_size; l++ )
                    {
                        unsigned int m = fData->getDetectorGeo()->getNeighbours()[i][l];
                        if( m != i && m < fData->getBorder().size() && ( fData->getBorder()[m] || fData->getImage()[m] ) )
                        {
                            fData->setImage( i, true );
                            i_neigh = true;
                            break;
                        }
                    }
                }
                if( !i_neigh )
                {
                    fData->setImage( i, false );
                }
            }
        }
        if( fData->getBrightNonImage()[i] )
        {
            if( fData->getImage()[i] )
            {
                fData->setBrightNonImage( i, false );
            }
            else if( fData->getBorder()[i] )
            {
                fData->setBrightNonImage( i, false );
            }
        }
    }
}

void VImageCleaning::addImageChannel( unsigned int i_channel )
{
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::addImageChannel" << endl;
    }
    if( i_channel < fData->getImage().size() )
    {
        fData->setImageUser( i_channel, 1 );
    }
}


void VImageCleaning::removeImageChannel( unsigned int i_channel )
{
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::removeImageChannel" << endl;
    }
    if( i_channel < fData->getImage().size() )
    {
        fData->setImageUser( i_channel, -1 );
    }
}


void VImageCleaning::resetImageChannel( unsigned int i_channel )
{
    if( fData->getDebugFlag() )
    {
        cout << "VImageCleaning::resetImageChannel" << endl;
    }
    if( i_channel < fData->getImage().size() )
    {
        fData->setImageUser( i_channel, 0 );
    }
}


/*
   trace correlation cleaning

*/
double getTraceCorrelationValue( double Amean, double Bmean,
                                 double Avar, double Bvar,
                                 vector < double > vA, vector < double > vB )
{
    if( Avar == 0. || Bvar == 0. )
    {
        return 0.;
    }
    double N = 0;
    
    for( unsigned int i = 0; i < vA.size() && i < vB.size() ; i++ )
    {
        N = N + ( vA[i] - Amean ) * ( vB[i] - Bmean );
    }
    
    return N / TMath::Sqrt( Avar * Bvar );
}

double getTraceMean( vector < double > vA )
{
    if( vA.size() == 0 )
    {
        return 0.;
    }
    
    double N = 0;
    for( unsigned int i = 0; i < vA.size(); i++ )
    {
        N = N + vA[i];
    }
    
    return N / double( vA.size() );
}

double getTraceVar( vector < double > vA, double Am )
{
    double N = 0;
    for( unsigned int i = 0; i < vA.size(); i++ )
    {
        N = N + ( vA[i] - Am ) * ( vA[i] - Am );
    }
    
    return N;
}

/*
   trace correlation cleaning

*/
void VImageCleaning::cleanImageTraceCorrelate( VImageCleaningRunParameter* iImageCleaningParameters )
{
    ///////////////////////////////////////////////
    // setting of image cleaning run parameters
    if( !iImageCleaningParameters )
    {
        return;
    }
    double sigNoiseThresh = iImageCleaningParameters->fCorrelationCleanBoardThresh;
    double corrThresh = iImageCleaningParameters->fCorrelationCleanCorrelThresh;
    double pixThresh = iImageCleaningParameters->fCorrelationCleanNpixThresh;
    
    fData->setBorderCorrelationCoefficient( 0. );
    
    vector < vector < double > > vImageTraces( fData->getDetectorGeo()->getNChannels( fData->getTelID() ), vector<double>( ( int )fData->getNSamples(), 0 ) );
    
    unsigned int nhits = fData->getReader()->getNumChannelsHit();
    for( unsigned int i = 0; i < nhits; i++ )
    {
        unsigned int i_channelHitID = 0;
        try
        {
            i_channelHitID = fData->getReader()->getHitID( i );
            if( i_channelHitID < fData->getHiLo().size() && !fData->getDead( fData->getHiLo()[i_channelHitID] )[i_channelHitID] )
            {
            
                fData->getTraceHandler()->setDigitalFilterParameters( fData->getDigitalFilterMethod(),
                        fData->getDigitalFilterUpSample(),
                        fData->getDigitalFilterPoleZero() );
                fData->getTraceHandler()->setTrace( fData->getReader(), fData->getNSamples(),
                                                    fData->getPeds( fData->getHiLo()[i_channelHitID] )[i_channelHitID],
                                                    i_channelHitID, i,
                                                    fData->getLowGainMultiplier_Trace()*fData->getHiLo()[i_channelHitID] );
                                                    
                vImageTraces[i_channelHitID] = fData->getTraceHandler()->getTrace();
            }
        }
        catch( ... )
        {
            if( fData->getDebugLevel() == 0 )
            {
                cout << "VImageCleaning::cleanImageTraceCorrelate, index out of range (fReader->getHitID) ";
                cout << i << "(Telescope " << fData->getTelID() + 1 << ", event " << fData->getReader()->getEventNumber() << ")" << endl;
                fData->setDebugLevel( 0 );
            }
            continue;
        }
        
    }
    
    if( vImageTraces.size() > 0 )
    {
        vector < double > avepulse( fData->getNSamples(), 0 );
        vector < unsigned int > ImagePixelList;
        vector < unsigned int > NearbyPixelList;
        double AvePulseMean = 0;
        double AvePulseVar = 0;
        int nimage = 0;
        unsigned int k = 0;
        
        for( unsigned int i = 0; i < vImageTraces.size(); i++ )
        {
            if( fData->getImage()[i] || fData->getBorder()[i] )
            {
                nimage++;
                ImagePixelList.push_back( i );
                for( unsigned int j = 0; j < fData->getNSamples(); j++ )
                {
                    avepulse[j] = avepulse[j] + vImageTraces[i][j];
                }
            }
        }
        
        if( nimage > 1 )
        {
            for( unsigned int j = 0; j < fData->getNSamples(); j++ )
            {
                avepulse[j] = avepulse[j] / double( nimage );
                AvePulseMean = AvePulseMean + avepulse[j];
            }
            AvePulseMean = AvePulseMean / double( fData->getNSamples() );
            
            for( unsigned int j = 0; j < fData->getNSamples(); j++ )
            {
                AvePulseVar = AvePulseVar + ( avepulse[j] - AvePulseMean ) * ( avepulse[j] - AvePulseMean );
            }
        }
        
        if( nimage > 3 && nimage < pixThresh )
        {
            //Build a list of pixels near the image pixels
            //Loop over image/border pixels
            for( unsigned int o = 0; o < ImagePixelList.size(); o++ )
            {
                unsigned int i = ImagePixelList[o];
                //Get the neighbours of this image pixel
                for( unsigned int j = 0; j < fData->getDetectorGeo()->getNNeighbours()[i]; j++ )
                {
                    //Check if it is already included in the neighbour list
                    bool have = false;
                    k = fData->getDetectorGeo()->getNeighbours()[i][j];
                    for( unsigned int p = 0; p < NearbyPixelList.size(); p++ )
                    {
                        if( NearbyPixelList[p] == k )
                        {
                            have = true;
                            break;
                        }
                    }
                    //Check if it is already included in the image pixel list
                    for( unsigned int p = 0; p < ImagePixelList.size(); p++ )
                    {
                        if( ImagePixelList[p] == k )
                        {
                            have = true;
                            break;
                        }
                    }
                    //Add to neighbour list if necessary
                    if( !have )
                    {
                        NearbyPixelList.push_back( k );
                    }
                }
            }
            
            
            for( unsigned int i = 0; i < NearbyPixelList.size(); i++ )
            {
                k =  NearbyPixelList[i];
                
                if( fData->getPedvars( fData->getCurrentSumWindow()[k], fData->getHiLo()[k] )[k] > 0. )
                {
                    double tMean = getTraceMean( vImageTraces[k] );
                    double tVar = getTraceVar( vImageTraces[k], tMean );
                    
                    
                    double corv = getTraceCorrelationValue( AvePulseMean, tMean, AvePulseVar, tVar, avepulse, vImageTraces[k] );
                    double sn = fData->getSums()[k] / fData->getPedvars( fData->getCurrentSumWindow()[k], fData->getHiLo()[k] )[k];
                    
                    // require that correlation coefficient and signal/noise is above certain thresholds
                    if( corv > corrThresh && sn > sigNoiseThresh )
                    {
                        fData->setBorder( k, true );
                        fData->setBorderCorrelationCoefficient( k, corv );
                    }
                    
                }
            }
        }
        recoverImagePixelNearDeadPixel();
        fillImageBorderNeighbours();
    }
}
