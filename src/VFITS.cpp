
#include "VFITS.h"

VFITS::VFITS( string anasum_file, string fits_file, string object_name, bool iOneFile, bool iPrint )
{
    fFile_anasum = anasum_file;
    fFile_FITS   = fits_file;
    fWriteOneFile = iOneFile;
    
    fTarget_Name = object_name;
    fTarget_Exposure = 0.;
    fTarget_RAJ2000 = 0.;
    fTarget_DecJ2000 = 0.;
    ctRunSum = 0;
    if( !readAnasumFile( iPrint ) )
    {
        exit( EXIT_FAILURE );
    }
}


/*****************************************************/
/* Print out cfitsio error messages and exit program */
/*****************************************************/
bool VFITS::printerror( int status )
{
    if( status )
    {
        fits_report_error( stderr, status );
    }
    return false;
}

/*******************************/
/* Read Anasum-Output-File     */
/*******************************/

bool VFITS::readAnasumFile( bool iPrint )
{

    if( !openFile( fFile_anasum , -1, false, iPrint ) )
    {
        cout << "VFITS::readAnasumFile() Error reading anasum file: " << fFile_anasum << endl;
        return false;
    }
    
    ctRunSum = getRunSummaryTree( 1 );
    ctRunSum->GetEntry( 0 );
    fTarget_RAJ2000 = ( float )ctRunSum->SkyMapCentreRAJ2000;
    fTarget_DecJ2000 = ( float )ctRunSum->SkyMapCentreDecJ2000;
    if( TMath::Abs( fTarget_RAJ2000 ) < 1.e-5 && TMath::Abs( fTarget_DecJ2000 ) < 1.e-5 )
    {
        fTarget_RAJ2000 = ( float )ctRunSum->TargetRAJ2000;
        fTarget_DecJ2000 = ( float )ctRunSum->TargetDecJ2000;
    }
    
    ctRunSum->GetEntry( ( ctRunSum->fChain->GetEntries() ) - 1 );
    if( ctRunSum->runOn == -1 )
    {
        fTarget_Exposure = ( float )ctRunSum->tOn;
    }
    else
    {
        cout << "VFITS::readAnasumFile() Error, no info for all runs " << endl;
        return false;
    }
    
    if( iPrint )
    {
        //  cout << "Histogram found with " << hSkyMap->GetEntries() << " entries" << endl;
        cout << "Pointing direction  (ra,dec): " << fTarget_RAJ2000 << " " << fTarget_DecJ2000 << endl;
        cout << "Exposure [s]: " << fTarget_Exposure << endl;
    }
    
    if( fWriteOneFile == true )
    {
        // create new FITS file
        int status = 0;
        string fFileName = fFile_FITS;
        //fFileName.insert(fFileName.size()-5,"_"+DiagName);
        if( fits_create_file( &fptr, fFileName.c_str(), &status ) )
        {
            return printerror( status );
        }
        if( fits_create_img( fptr,  DOUBLE_IMG, 0, 0, &status ) )
        {
            return printerror( status );
        }
        writeFITSInfo( iPrint );
        cout << "Created FITS file: " << fFileName << endl;
        
    }
    
    return true;
}

/*********************************************************************/
/* Calculate cumulative Significance and store as table in FITS file */
/*********************************************************************/

bool VFITS::writeCumSignificance( bool iPrint )
{
    cout << " Write cumulative significance to fits file" << endl;
    TGraph* gCumSig = calcCumulativeSig( 1 );
    if( iPrint )
    {
        cout << " Calculated the cumulative significance" << endl;
    }
    
    writeTGraphFits( gCumSig, "cumulativeSignificance", "time", "cumSig", "min", "sigma", iPrint );
    if( iPrint )
    {
        cout << " transformed cumulative significance into FITS-table" << endl;
    }
    
    return true;
}

/*********************************************************************/
/* Calculate integral flux in monthly bins and store as table in FITS file */
/*********************************************************************/

bool VFITS::writeMonthlyFlux( bool iPrint, string outfile )
{
    cout << " Write nightly flux to fits file" << endl;
    VFluxCalculation flux2( fFile_anasum );
    flux2.setSpectralParameters( 0.35, 1.0, -2.2 );
    flux2.setSignificanceParameters( -9, -9 );
    flux2.calculateIntegralFlux( 0.35 );
    
    VLightCurve flux;
    flux.initializeTeVLightCurve( fFile_anasum, 30, 54500, 59200 );
    flux.setSpectralParameters( 0.35, 1.0, -2.2 );
    flux.setSignificanceParameters( -9, -9 );
    flux.fill( 0.35 );
    flux.plotLightCurve();
    TGraphAsymmErrors* gFlux = flux.getLightCurveGraph();
    
    
    if( iPrint )
    {
        cout << " Calculated the integral flux " << endl;
    }
    writeTGraphAsymmErrorsFits( gFlux, "MonthlyFlux", "Date", "F(>0.35 TeV)", "MJD", "1/[cm^2*s]", iPrint );
    if( iPrint )
    {
        cout << " transformed monthly flux into FITS-table" << endl;
    }
    
    if( outfile != "" )
    {
        ofstream out;
        out.open( outfile.c_str() );
        out << "MJD\tF(>350GeV) [/cm^2/s]\tFlux in C.U." << endl;
        TString line;
        double mjd, flx, flxE, flxCrab, flxCrabE;
        
        for( int i = 0; i < gFlux->GetN(); i++ )
        {
            mjd = gFlux->GetX()[i];
            flx = gFlux->GetY()[i];
            flxE = ( gFlux->GetErrorYhigh( i ) + gFlux->GetErrorYlow( i ) ) / 2.0;
            
            flxCrab = flux2.getFluxVsCrab( flx, 0.35, 2.2 );
            flxCrabE = flux2.getFluxVsCrab( flxE, 0.35, 2.2 );
            
            line.Form( "%d\t%.3e\t%.3e\t%.2f\t%.2f\n", ( int )mjd, flx, flxE, flxCrab, flxCrabE );
            out << line ;
            
        }
        flux2.getFlux( -1, flx, flxE, mjd );
        flxCrab = flux2.getFluxVsCrab( flx, 0.35, 2.2 );
        flxCrabE = flux2.getFluxVsCrab( flxE, 0.35, 2.2 );
        
        line.Form( "Total:\t%.3e\t%.3e\t%.2f\t%.2f\n", flx, flxE, flxCrab, flxCrabE );
        out << line << endl;
    }
    
    return true;
}

/*********************************************************************/
/* Calculate integral flux in nightly bins and store as table in FITS file */
/*********************************************************************/

bool VFITS::writeNightlyFlux( bool iPrint, string outfile )
{
    cout << " Write nightly flux to fits file" << endl;
    
    VFluxCalculation flux( fFile_anasum );
    flux.setSpectralParameters( 0.2, 1.0, -2.5 );
    flux.setSignificanceParameters( -9., -9. );
    flux.calculateIntegralFlux( 0.2 );
    vector< VFluxDataPoint > iFluxData = flux.getFluxDataVector();
    const VFluxDataPoint* iFluxDataCombined = flux.getFluxDataCombined();
    
    if( iPrint )
    {
        cout << " Calculated the integral flux " << endl;
    }
    writeLightCurveFITS( iFluxData, "NightlyFlux", "Date", "F(>0.2 TeV)", "MJD", "1/[cm^2*s]", iPrint );
    if( iPrint )
    {
        cout << " transformed nightly flux into FITS-table" << endl;
    }
    
    if( outfile != "" )
    {
        ofstream out;
        out.open( outfile.c_str() );
        out << "MJD\tF(>200GeV) [/cm^2/s]\tFlux in C.U." << endl;
        TString line;
        double flxCrab = 0.;
        double flxCrabE = 0.;
        
        for( unsigned int i = 0; i < iFluxData.size(); i++ )
        {
            flxCrab  = VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( 0.2, iFluxData[i].fFlux, 2.5 );
            flxCrabE = VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( 0.2, iFluxData[i].fFluxCI_1sigma, 2.5 );
            
            line.Form( "%d\t%.3e\t%.3e\t%.2f\t%.2f\n", ( int )iFluxData[i].fMJD, iFluxData[i].fFlux, iFluxData[i].fFluxCI_1sigma, flxCrab, flxCrabE );
            out << line ;
            
        }
        flxCrab  = VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( 0.2, iFluxDataCombined->fFlux, 2.5 );
        flxCrabE = VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( 0.2, iFluxDataCombined->fFluxCI_1sigma, 2.5 );
        
        line.Form( "Total:\t%.3e\t%.3e\t%.2f\t%.2f\n", iFluxDataCombined->fFlux, iFluxDataCombined->fFluxCI_1sigma, flxCrab, flxCrabE );
        out << line << endl;
    }
    
    return true;
}

/********************************************************************/
/* Calculate 1-dim Significance Distribution and store in FITS file */
/********************************************************************/

bool VFITS::writeSignificanceDistribution( bool iPrint )
{
    cout << " Write significance distribution to fits file" << endl;
    // get all histograms
    TH2D* hmap_stereo_sig = 0;
    TH2D* hmap_stereo_on = 0;
    
    hmap_stereo_sig = ( TH2D* )getHistogram( "hmap_stereo_sig", 1, "skyHistograms" );
    hmap_stereo_on = ( TH2D* )getHistogram( "hmap_stereo_on", 1, "skyHistograms" );
    if( iPrint )
    {
        cout << " Got necessary histograms from file" << endl;
    }
    
    // get 1D histograms
    
    vector<pair<TH1D*, string> > vhist;
    // with source
    TH1D* hsig_1Dall  = get_Bin_Distribution( hmap_stereo_sig, -1, 1.2, 0., false, hmap_stereo_on );
    vhist.push_back( make_pair( hsig_1Dall, "xwy" ) );
    //without source
    TH1D* hsig_1D  = get_Bin_Distribution( hmap_stereo_sig, -1, 1.2, 0.4, false, hmap_stereo_on );
    vhist.push_back( make_pair( hsig_1D, "y" ) );
    //without stars -- not calculated yet -> use histo "with source"
    vhist.push_back( make_pair( hsig_1Dall, "y" ) );
    //without all -- not calculated yet -> use histo "without source"
    vhist.push_back( make_pair( hsig_1D, "y" ) );
    if( iPrint )
    {
        cout << " Calculated 1-dim significance distribution" << endl;
    }
    
    char* colNames[10] = {( char* )"Sigma", ( char* )"SigmaWidth", ( char* )"Counts", ( char* )"CountsMinusSource", ( char* )"CountsMinusStars", ( char* )"CountsMinusAll"};
    char* colUnits[10] = {( char* )"sigma", ( char* )"Delta sigma", ( char* )"counts", ( char* )"counts", ( char* )"counts", ( char* )"counts"};
    char* colDataFormat[10] = {( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D"};
    writeVecTH1DFits( vhist, "SignificanceDistributions", colNames, colUnits, colDataFormat, iPrint );
    if( iPrint )
    {
        cout << " transformed significance distribution into FITS-table" << endl;
    }
    
    //     vector<int> hdunums;
    //     hdunums.push_back(WriteTH1DFits(hsig_1D,"significanceDist", "significance", "w/o Source", "sigma", "counts ", iPrint));
    //     hdunums.push_back(WriteTH1DFits(hsig_1Dall,"significanceDistAll", "significance", "w/ Source", "sigma", "counts ", iPrint));
    //     hdunums.push_back(WriteTH1DFits(hsig_1Dtest,"significanceDistTest", "significance", "w/o Source", "sigma", "counts ", iPrint));
    
    return true;
}
/********************************************/
/* Get LightCurve and store in FITS file    */
/********************************************/

bool VFITS::writeLightCurve( bool iPrint )
{
    cout << " Write light curve to fits file" << endl;
    //get Rates from file
    TGraphErrors* gRates = ( TGraphErrors* )getHistogram( "gRateRun", 1, "ratePlots" );
    if( iPrint )
    {
        cout << " Got the rate TGraphErrors" << endl;
    }
    
    writeTGraphErrorsFits( gRates, "Rates", "run", "rate", " ", "1/min", iPrint );
    if( iPrint )
    {
        cout << " transformed rates into FITS-table" << endl;
    }
    
    
    return true;
}

/*****************************************/
/* Get ThetaSq and store in FITS file    */
/*****************************************/

bool VFITS::writeThetaSquareDistribution( bool iPrint )
{
    cout << " Write ThetaSquare distribution to fits file" << endl;
    //get ThetaSquare from file
    vector<pair<TH1D*, string> > vhist;
    //ON data
    TH1D* hThetaON = ( TH1D* )getHistogram( "htheta2_on", 1, "stereoParameterHistograms" );
    hThetaON->Rebin( 5 );
    vhist.push_back( make_pair( hThetaON, "xwy" ) );
    //OFF data
    TH1D* hThetaOFF = ( TH1D* )getHistogram( "htheta2_off", 1, "stereoParameterHistograms" );
    hThetaOFF->Rebin( 5 );
    vhist.push_back( make_pair( hThetaOFF, "y" ) );
    //difference (ON minus OFF)
    TH1D* hThetaDiff = ( TH1D* )getHistogram( "htheta2_diff", 1, "stereoParameterHistograms" );
    hThetaDiff->Rebin( 5 );
    vhist.push_back( make_pair( hThetaDiff, "ye" ) );
    if( iPrint )
    {
        cout << " Got the theta2 histograms" << endl;
    }
    
    char* colNames[10] = {( char* )"ThetaSq", ( char* )"ThetaSqWidth", ( char* )"ON", ( char* )"OFF", ( char* )"Excess", ( char* )"ExcessError"};
    char* colUnits[10] = {( char* )"deg^2", ( char* )"Delta deg^2", ( char* )"counts", ( char* )"counts", ( char* )"counts", ( char* )"counts"};
    char* colDataFormat[10] = {( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D"};
    writeVecTH1DFits( vhist, "ThetaSquare", colNames, colUnits, colDataFormat, iPrint );
    if( iPrint )
    {
        cout << " transformed rates into FITS-table" << endl;
    }
    
    
    return true;
}

/************************************************/
/* Get EnergySpectrum and store in FITS file    */
/************************************************/
bool VFITS::writeEnergySpectrum( bool iPrint )
{
    cout << " Write EnergySpectrum to fits file" << endl;
    //double integral = 0.;
    //get data from SummaryTree
    VEnergySpectrum fEspec( fFile_anasum );
    if( iPrint )
    {
        cout << " define VEnergySpectrum variable" << endl;
    }
    //combine runs and calculate energy Spectrum
    fEspec.combineRuns();
    if( iPrint )
    {
        cout << " combine all runs" << endl;
    }
    
    TGraphErrors* gEspec = ( TGraphErrors* )fEspec.getEnergySpectrumGraph();
    if( iPrint )
    {
        cout << " Got the energy spectrum TGraphErrors" << endl;
    }
    
    if( gEspec == 0 )
    {
        return false ;
    }
    //integral = getFluxIntegral(gEspec, 0.2, iPrint);
    
    
    //Transform to Fits
    writeTGraphErrorsFits( gEspec, "EnergySpectrum", "Energy", "dN_over_dE", "log10[TeV]", "1/[cm^2*TeV*s]", iPrint );
    if( iPrint )
    {
        cout << " transformed EnergySpectrum into FITS-table" << endl;
    }
    return true;
}


/**********************************************************/
/* Get SignificanceSkyMap-IMAGE and store in FITS file    */
/**********************************************************/
bool VFITS::writeSignificanceSkyMap( bool iPrint )
{
    cout << " Write significance sky map to fits file" << endl;
    //get Rates from file
    TH2D* hSkyMap = ( TH2D* )getHistogram( "hmap_stereo_sig", 1, "skyHistograms" );
    if( iPrint )
    {
        cout << " Got the SignificanceSkyMap TH2D" << endl;
    }
    
    //Transform to Fits
    createImageFitsFile( hSkyMap , "SignificanceMap", iPrint );
    if( iPrint )
    {
        cout << " transformed SkyMap image into FITS-table" << endl;
    }
    return true;
}


/****************************************************/
/* Get ExcessSkyMap-IMAGE and store in FITS file    */
/****************************************************/
bool VFITS::writeExcessSkyMap( bool iPrint )
{
    cout << " Write excess sky map to fits file" << endl;
    //get Rates from file
    TH2D* hSkyMap = ( TH2D* )getHistogram( "hmap_stereo_diff", 1, "skyHistograms" );
    if( iPrint )
    {
        cout << " Got the ExcessSkyMap TH2D" << endl;
    }
    
    //Transform to Fits
    createImageFitsFile( hSkyMap , "ExcessSkyMap", iPrint );
    if( iPrint )
    {
        cout << " transformed SkyMap image into FITS-table" << endl;
    }
    return true;
}

/****************************************************/
/* Get AlphaSkyMap-IMAGE and store in FITS file    */
/****************************************************/
bool VFITS::writeAlphaSkyMap( bool iPrint )
{
    cout << " Write alpha sky map to fits file" << endl;
    //get Rates from file
    TH2D* hSkyMap = ( TH2D* )getHistogram( "hmap_alphaNorm_off", 1, "skyHistograms" );
    if( iPrint )
    {
        cout << " Got the AlphaSkyMap TH2D" << endl;
    }
    
    //Transform to Fits
    createImageFitsFile( hSkyMap , "AlphaSkyMap", iPrint );
    if( iPrint )
    {
        cout << " transformed SkyMap image into FITS-table" << endl;
    }
    return true;
}

/****************************************************/
/* Get RawCountsSkyMap-IMAGE and store in FITS file    */
/****************************************************/
bool VFITS::writeRawCountsSkyMap( bool iPrint )
{
    cout << " Write raw counts sky map to fits file" << endl;
    //get Rates from file
    TH2D* hSkyMap = ( TH2D* )getHistogram( "hxyoff_stereo_diff", 1, "skyHistograms" );
    if( iPrint )
    {
        cout << " Got the RawCountsSkyMap TH2D" << endl;
    }
    
    //Transform to Fits
    createImageFitsFile( hSkyMap , "RawCountsSkyMap", iPrint );
    if( iPrint )
    {
        cout << " transformed SkyMap image into FITS-table" << endl;
    }
    return true;
}

//*******************NECESSARY METHODES************************

double VFITS::getFluxIntegral( TGraphErrors* gEspec, double minE, bool iPrint )
{
    double x, y;
    double integral = 0;
    for( int i = 1; i < gEspec->GetN(); i++ )
    {
        gEspec->GetPoint( i - 1, x, y );
        if( pow( 10, x ) < minE )
        {
            continue;
        }
        double minX = pow( 10, x );
        double yVal1 = y;
        gEspec->GetPoint( i, x, y );
        double maxX = pow( 10, x );
        double yVal2 = y;
        integral = integral + ( yVal1 + yVal2 ) / 2.*( maxX - minX );
        if( iPrint )
        {
            cout << "10^" << x << " = " << pow( 10, x );
            printf( "  Integral = %.19f \n", integral );
        }
    }
    return integral;
    
}


/**************************************************************/
/* Transform TH1D, TGraph or TGraphError into  FITS BinTable  */
/**************************************************************/


int VFITS::writeVecTH1DFits( vector<pair<TH1D*, string> > vhist, string DiagName,  char* tType[], char* tUnit[], char* tForm[] ,  bool iPrint )
{
    /*To use WriteVecTH1DFits() you have to specify which histograms
     and which columns of that histograms should be used. The strings
     in vector<pair> should include x,w,y and/or e:
     x is x_axis(1.col.)
     w is BinWidth(2.col.)
     y is y-Value(3.col.)
     e is y-Error(4.col.)
    */
    
    if( iPrint )
    {
        cout << " Write Vector of TH1Ds to fits file ..." << endl;
    }
    
    //create a table
    int size =  0;
    int HduNum = -99;
    vector< double> runValues;
    vector< vector<double> > table;
    
    //get the histograms and the column information from the vector
    for( int i = 0; i < vhist[0].first->GetNbinsX(); i++ )
    {
        runValues.clear();
        for( unsigned int k = 0; k < vhist.size(); k++ )
        {
            if( vhist[0].first->GetNbinsX() != vhist[k].first->GetNbinsX() )
            {
                cout << "Histograms don't have the same Num. of bins" << endl;
                return -1;
            }
            if( vhist[k].second.find( "x", 0 ) < 10 )
            {
                runValues.push_back( vhist[k].first->GetBinCenter( i ) );    //1.row = X values
            }
            if( vhist[k].second.find( "w", 0 ) < 10 )
            {
                runValues.push_back( vhist[k].first->GetBinWidth( i ) );    //2.row = bin width
            }
            if( vhist[k].second.find( "y", 0 ) < 10 )
            {
                runValues.push_back( vhist[k].first->GetBinContent( i ) );    //3.row = Y values
            }
            if( vhist[k].second.find( "e", 0 ) < 10 )
            {
                runValues.push_back( vhist[k].first->GetBinError( i ) );    //4.row = Y Error
            }
            
        }
        if( i == 0 )
        {
            size = int( runValues.size() );
        }
        table.push_back( runValues );
    }
    
    //change column names if they are equal
    if( iPrint )
    {
        cout << "   Size of table columns = " << size << endl;
    }
    for( int g = 0; g < size ; g++ )
    {
        string type1( tType[g] );
        for( int h = g; h < size; h++ )
        {
            if( h == g )
            {
                continue;
            }
            string type2( tType[h] );
            if( type1 == type2 )
            {
                cout << "Some colums have the same name! Give them different names, please!" << endl;
                return -1;
            }
        }
    }
    HduNum = createTableFitsFile( table, tType, tUnit, tForm, DiagName, iPrint );
    if( iPrint )
    {
        cout << "   Wrote table into FITS-table" << endl;
    }
    
    return HduNum;
}


//**************************************************************


int VFITS::writeTH1DFits( TH1D* h, string DiagName, string x_name, string y_name, string x_unit, string y_unit, bool iPrint )
{

    //create a table
    int HduNum = -99;
    vector< double> runValues;
    vector< vector<double> > table;
    
    for( int i = 0; i < h->GetNbinsX(); i++ )
    {
        runValues.clear();
        runValues.push_back( h->GetBinCenter( i ) ); //1.row = X values
        runValues.push_back( h->GetBinContent( i ) ); //2.row = Y values
        runValues.push_back( h->GetBinError( i ) ); //3.row = Y Error
        table.push_back( runValues );
        
    }
    if( iPrint )
    {
        cout << "   Got X and Y values from histogram and stored them in a table" << endl;
    }
    // define Names, Units and DataForms for new FITS BinTable
    // x_name.insert(0,"X ");
    // y_name.insert(0,"Y ");
    char* tType[3] = {const_cast<char*>( x_name.c_str() ), const_cast<char*>( y_name.c_str() ), ( char* )" Error"};
    char* tUnit[3] = {const_cast<char*>( x_unit.c_str() ), const_cast<char*>( y_unit.c_str() ), ( char* )" "};
    char* tForm[3] = {( char* )"1D", ( char* )"1D", ( char* )"1D"};
    if( iPrint )
    {
        cout << "   Set names, units and dataformats for different colums " << endl;
    }
    
    HduNum = createTableFitsFile( table, tType, tUnit, tForm, DiagName, iPrint );
    if( iPrint )
    {
        cout << "   Wrote table into FITS-table" << endl;
    }
    
    return HduNum;
}

//***************************************************

int VFITS::writeTGraphFits( TGraph* g, string DiagName, string x_name, string y_name, string x_unit, string y_unit, bool iPrint )
{

    //create a table
    int HduNum = -99;
    vector< double> runValues;
    vector< vector<double> > table;
    double x1 = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        double x = 0.;
        double y = 0.;
        g->GetPoint( i, x, y );
        runValues.clear();
        runValues.push_back( x );  //1.row = X values
        // Note: this is wrong and should be fixed (delta_x)
        runValues.push_back( x - x1 ); //2.row = delta_X values
        runValues.push_back( y );  //3.row = Y values
        table.push_back( runValues );
        x1 = x;
        
    }
    if( iPrint )
    {
        cout << "   Got X and Y values from TGraph and stored them in a table" << endl;
    }
    
    // define Names, Units and DataForms for new FITS BinTable
    string DeltaX( "Delta_" + x_name );
    char* tType[3] = {const_cast<char*>( x_name.c_str() ), const_cast<char*>( DeltaX.c_str() ) , const_cast<char*>( y_name.c_str() )};
    char* tUnit[3] = {const_cast<char*>( x_unit.c_str() ), const_cast<char*>( x_unit.c_str() ) , const_cast<char*>( y_unit.c_str() )};
    char* tForm[3] = {( char* )"1D", ( char* )"1D", ( char* )"1D"};
    if( iPrint )
    {
        cout << "   Set names, units and dataformats for different colums " << endl;
    }
    
    HduNum = createTableFitsFile( table, tType, tUnit, tForm, DiagName, iPrint );
    if( iPrint )
    {
        cout << "   Wrote table into FITS-table" << endl;
    }
    
    return HduNum;
}

//***************************************************

int VFITS::writeTGraphErrorsFits( TGraphErrors* g, string DiagName, string x_name, string y_name, string x_unit, string y_unit, bool iPrint )
{

    //create a table
    int HduNum = -99;
    vector< double> runValues;
    vector< vector<double> > table;
    double x1 = 0.;
    for( int i = 0; i < g->GetN(); i++ )
    {
        double x = 0.;
        double y = 0.;
        g->GetPoint( i, x, y );
        runValues.clear();
        runValues.push_back( x ); //1.row = X values
        runValues.push_back( x - x1 ); //2.row = delta_X values
        runValues.push_back( y ); //3.row = Y values
        runValues.push_back( g->GetErrorY( i ) ); //4.row = Error of Y
        table.push_back( runValues );
        x1 = x;
    }
    if( iPrint )
    {
        cout << "   Got X value, Y values and Y_errors from TGraphErrors and stored them in a table" << endl;
    }
    
    // define Names, Units and DataForms for new FITS BinTable
    string DeltaX( "Delta_" + x_name );
    string ErrorY( y_name + "_Error" );
    char* tType[4] = {const_cast<char*>( x_name.c_str() ), const_cast<char*>( DeltaX.c_str() ) , const_cast<char*>( y_name.c_str() ), const_cast<char*>( ErrorY.c_str() )};
    char* tUnit[4] = {const_cast<char*>( x_unit.c_str() ), const_cast<char*>( x_unit.c_str() ), const_cast<char*>( y_unit.c_str() ), const_cast<char*>( y_unit.c_str() ),};
    char* tForm[4] = {( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D"};
    if( iPrint )
    {
        cout << "   Set names, units and dataformats for different colums " << endl;
    }
    
    HduNum = createTableFitsFile( table, tType, tUnit, tForm, DiagName, iPrint );
    if( iPrint )
    {
        cout << "   Wrote table into FITS-table" << endl;
    }
    
    return HduNum;
}

int VFITS::writeLightCurveFITS( vector< VFluxDataPoint > iFluxData, string DiagName, string x_name, string y_name, string x_unit, string y_unit, bool iPrint )
{

    //create a table
    int HduNum = -99;
    vector< double> runValues;
    vector< vector<double> > table;
    for( unsigned int i = 0; i < iFluxData.size(); i++ )
    {
        runValues.clear();
        runValues.push_back( iFluxData[i].fMJD ); //1.row = X values
        runValues.push_back( iFluxData[i].fMJD_Stop - iFluxData[i].fMJD_Start ); //2.row = delta_X values
        runValues.push_back( iFluxData[i].fFlux ); //3.row = Y values
        runValues.push_back( iFluxData[i].fFluxCI_1sigma ); //4.row = Error of Y
        table.push_back( runValues );
    }
    if( iPrint )
    {
        cout << "   Got X value, Y values and Y_errors from TGraphErrors and stored them in a table" << endl;
    }
    
    // define Names, Units and DataForms for new FITS BinTable
    string DeltaX( "Delta_" + x_name );
    string ErrorY( y_name + "_Error" );
    char* tType[4] = {const_cast<char*>( x_name.c_str() ), const_cast<char*>( DeltaX.c_str() ) , const_cast<char*>( y_name.c_str() ), const_cast<char*>( ErrorY.c_str() )};
    char* tUnit[4] = {const_cast<char*>( x_unit.c_str() ), const_cast<char*>( x_unit.c_str() ), const_cast<char*>( y_unit.c_str() ), const_cast<char*>( y_unit.c_str() ),};
    char* tForm[4] = {( char* )"1D", ( char* )"1D", ( char* )"1D", ( char* )"1D"};
    if( iPrint )
    {
        cout << "   Set names, units and dataformats for different colums " << endl;
    }
    
    HduNum = createTableFitsFile( table, tType, tUnit, tForm, DiagName, iPrint );
    if( iPrint )
    {
        cout << "   Wrote table into FITS-table" << endl;
    }
    
    return HduNum;
}


/*************************************************/
/* Create FITS file with BinTable extension      */
/*************************************************/

int VFITS::createTableFitsFile( vector< vector<double> > Table , char* ttype[] , char* tunit[], char* tform[], string DiagName, bool iPrint )
{
    int hdunum = -99;
    int nRows = int( Table.size() );
    int nCol = int( Table[0].size() );
    //    cout<<"NCol: "<<nCol<<endl;
    int status = 0;
    string fFileName = fFile_FITS;
    fFileName.insert( fFileName.size() - 5, "_" + DiagName );
    
    
    if( fWriteOneFile == false )
    {
        // create new FITS file
        fitsfile* fFitsPtr = NULL;            /* pointer to the FITS file, defined in fitsio.h */
        fptr = fFitsPtr;
        if( fits_create_file( &fptr, fFileName.c_str(), &status ) )
        {
            return printerror( status );
        }
        if( fits_create_img( fptr,  DOUBLE_IMG, 0, 0, &status ) )
        {
            return printerror( status );
        }
        writeFITSInfo( iPrint );
        if( iPrint )
        {
            cout << "       Created FITS file: " << fFileName << endl;
        }
        //for (int j = 0; j<nCol; j++) { cout<<"tForm["<<j<<"] = "<<tform[j]<<endl;}
    }
    
    if( fits_create_tbl( fptr, BINARY_TBL, nRows, nCol, ttype, tform, tunit, DiagName.c_str() , &status ) )
    {
        return printerror( status );
    }
    ffghdn( fptr, &hdunum );
    if( iPrint )
    {
        cout << "      Created BinaryTable ..." << endl;
    }
    
    for( int k = 0; k < nCol; k++ )
    {
        for( int i = 0; i < nRows; i++ )
        {
            if( fits_write_col( fptr, TDOUBLE, k + 1, i + 1, 1, 1, &Table[i][k] , &status ) )
            {
                return printerror( status );
            }
            if( iPrint )
            {
                cout << "        Write Table[" << i << "][" << k << "] = " << Table[i][k] << endl;
            }
        }
        if( iPrint )
        {
            cout << "      Write " << k << ".Column in BinaryTable ..." << endl;
        }
    }
    if( fWriteOneFile == false )
    {
        if( fits_close_file( fptr, &status ) )
        {
            return printerror( status ) ;
        }
        if( iPrint )
        {
            cout << "      Close FITS file: " << fFileName << endl;
        }
    }
    return hdunum;
}


/********************************************/
/* Create FITS file with IMAGE extension    */
/********************************************/

int VFITS::createImageFitsFile( TH2D* hSkyMap , string DiagName, bool iPrint )
{

    if( !hSkyMap )
    {
        cout << "Error: no histogram available. Cannot find histograms: " << hSkyMap  << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    
    int hdunum = -99;
    long  fpixel, nelements;
    double* array[hSkyMap->GetNbinsX()];
    int status = 0;
    string fFileName = fFile_FITS;
    fFileName.insert( fFileName.size() - 5, "_" + DiagName );
    
    // initialize FITS image parameters
    int bitpix   =  DOUBLE_IMG;                   /* 16-bit unsigned short pixel values   */
    const long naxis    =   2;                    /* 2-dimensional image                  */
    long naxes[naxis] =                           /* image is 400 pixels wide by 400 rows */
    {
        hSkyMap->GetNbinsX(), hSkyMap->GetNbinsY()
    };
    
    // allocate memory for the whole image
    array[0] = ( double* )malloc( naxes[0] * naxes[1] * sizeof( double ) );
    
    // initialize pointers to the start of each row of the image
    for( int ii = 1; ii < naxes[1]; ii++ )
    {
        array[ii] = array[ii - 1] + naxes[0];
    }
    
    if( fWriteOneFile == false )
    {
        // create new FITS file
        fitsfile* fFitsPtr = NULL;            /* pointer to the FITS file, defined in fitsio.h */
        fptr = fFitsPtr;
        if( fits_create_file( &fptr, fFileName.c_str(), &status ) )
        {
            return printerror( status );
        }
        if( fits_create_img( fptr,  DOUBLE_IMG, 0, 0, &status ) )
        {
            return printerror( status );
        }
        writeFITSInfo( iPrint );
        if( iPrint )
        {
            cout << "       Created FITS file: " << fFileName << endl;
        }
    }
    
    if( fits_create_img( fptr,  bitpix, naxis, naxes, &status ) )
    {
        return printerror( status );
    }
    ffghdn( fptr, &hdunum );
    if( iPrint )
    {
        cout << "      Created image ..." << endl;
    }
    
    // initialize the values in the image with a linear ramp function
    for( int a = 0; a < naxes[0]; a++ )
    {
        for( int b = 0; b < naxes[1]; b++ )
        {
            array[a][b] = hSkyMap->GetBinContent( naxes[1] - b, a + 1 );
            if( array[a][b] < -9000 )
            {
                array[a][b] = NAN;
            }
        }
    }
    fpixel = 1;                                   /* first pixel to write      */
    nelements = naxes[0] * naxes[1];              /* number of pixels to write */
    
    // write the array of unsigned integers to the FITS file
    if( fits_write_img( fptr, TDOUBLE, fpixel, nelements, array[0], &status ) )
    {
        return printerror( status );
    }
    writeFITSimageInfo( naxes, hSkyMap , DiagName );
    
    free( array[0] );                             /* free previously allocated memory */
    
    if( fWriteOneFile == false )
    {
        if( fits_close_file( fptr, &status ) )
        {
            return printerror( status ) ;
        }
        if( iPrint )
        {
            cout << "      Close FITS file: " << fFileName << endl;
        }
    }
    return hdunum;
}

/*******************************************/
/* Write useful IMAGE infos into FITS file */
/*******************************************/

bool VFITS::writeFITSimageInfo( long* naxes, TH2D* hSkyMap , string DiagName )
{
    int status = 0;
    
    char extname[100];
    sprintf( extname, "%s", DiagName.c_str() );
    if( fits_write_key_str( fptr, "EXTNAME", extname , "Extension Name" , &status ) )
    {
        return printerror( status );
    }
    
    char radecsys[] = "FK5";
    if( fits_update_key( fptr, TSTRING, "RADECSYS", radecsys, ( char* )"WCS for this file", &status ) )
    {
        return printerror( status );
    }
    
    float equinoxsys = 2000.;
    if( fits_update_key( fptr, TFLOAT, "EQUINOX", &equinoxsys , ( char* )"Epoch of coordinate system", &status ) )
    {
        return printerror( status );
    }
    
    
    char ctype1[] = "RA---TAN";
    if( fits_update_key( fptr, TSTRING, "CTYPE1", ctype1, ( char* )"Axis type for dim 1 (RA)", &status ) )
    {
        return printerror( status );
    }
    
    if( fits_update_key( fptr, TFLOAT, "CRVAL1", &fTarget_RAJ2000, ( char* )"Sky coord of 1st axis (deg)", &status ) )
    {
        return printerror( status );
    }
    
    float pix1_origin = hSkyMap->GetNbinsX() - hSkyMap->GetXaxis()->FindFixBin( 0. ) + 1;
    if( fits_update_key( fptr, TFLOAT, "CRPIX1", &pix1_origin, ( char* )"Reference point of pixel location axis 1", &status ) )
    {
        return printerror( status );
    }
    
    float xbinwidth = -1.*hSkyMap->GetXaxis()->GetBinWidth( 1 );
    if( fits_update_key( fptr, TFLOAT, "CDELT1", &xbinwidth, ( char* )"X degrees per pixel", &status ) )
    {
        return printerror( status );
    }
    
    
    char ctype2[] = "DEC--TAN";
    if( fits_update_key( fptr, TSTRING, "CTYPE2", ctype2, ( char* )"Axis type for dim 2 (DEC)", &status ) )
    {
        return printerror( status );
    }
    
    if( fits_update_key( fptr, TFLOAT, "CRVAL2", &fTarget_DecJ2000, ( char* )"Sky coord of 2nd axis (deg)", &status ) )
    {
        return printerror( status );
    }
    
    
    float pix2_origin = hSkyMap->GetYaxis()->FindFixBin( 0. );
    if( fits_update_key( fptr, TFLOAT, "CRPIX2", &pix2_origin, ( char* )"Reference point of pixel location axis 2", &status ) )
    {
        return printerror( status );
    }
    
    float ybinwidth = hSkyMap->GetYaxis()->GetBinWidth( 1 );
    if( fits_update_key( fptr, TFLOAT, "CDELT2", &ybinwidth, ( char* )"Y degrees per pixel", &status ) )
    {
        return printerror( status );
    }
    
    
    return true;
}


/*************************************/
/* Write useful infos into FITS file */
/*************************************/

bool VFITS::writeFITSInfo( bool iPrint )
{

    TDatime* temps = new TDatime();
    temps->Set();
    
    int YYYY = temps->GetYear();
    int MM = temps->GetMonth();
    int DD = temps->GetDay();
    int HH = temps->GetHour();
    int MN = temps->GetMinute();
    int SS = temps->GetSecond();
    
    delete temps;
    
    if( iPrint )
    {
        cout << "--> current date is: " << DD << MM << YYYY << "   time:" << HH << ":" << MN << ":" << SS << endl;
    }
    
    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    
    int status = 0;
    
    char author[] = "VERITAS_Collaboration";
    if( fits_update_key( fptr, TSTRING, ( char* )"AUTHOR", author, ( char* )"Author", &status ) )
    {
        return printerror( status );
    }
    
    char name[] = "VERITAS";
    if( fits_update_key( fptr, TSTRING, ( char* )"TELESCOP", name, ( char* )"Telescope name", &status ) )
    {
        return printerror( status );
    }
    char iTarget[600];
    sprintf( iTarget, "%s", fTarget_Name.c_str() );
    if( fits_update_key( fptr, TSTRING, ( char* )"OBJECT", iTarget, ( char* )"Name of the object", &status ) )
    {
        return printerror( status );
    }
    
    if( fits_update_key( fptr, TFLOAT, ( char* )"RA_OBJ", &fTarget_RAJ2000, ( char* )"RA (deg) of the object", &status ) )
    {
        return printerror( status );
    }
    
    if( fits_update_key( fptr, TFLOAT, ( char* )"DEC_OBJ", &fTarget_DecJ2000, ( char* )"Dec (deg) of the Object", &status ) )
    {
        return printerror( status );
    }
    
    float equinoxsys = 2000.;
    if( fits_update_key( fptr, TFLOAT, ( char* )"EQUINOX", &equinoxsys , ( char* )"Epoch of coordinate system", &status ) )
    {
        return printerror( status );
    }
    
    if( fits_update_key( fptr, TFLOAT, ( char* )"EXPOSURE", &fTarget_Exposure, ( char* )"Total exposure time in sec", &status ) )
    {
        return printerror( status );
    }
    
    
    string str_temp;
    stringstream dd;
    char* date;
    dd << YYYY << "-" << MM << "-" << DD;
    dd >> str_temp;
    date = new char[str_temp.size() + 1];
    strcpy( date, str_temp.c_str() );
    
    if( fits_update_key( fptr, TSTRING, ( char* )"DATE", date, ( char* )"FITS file creation date (yyyy-mm-dd)", &status ) )
    {
        return printerror( status );
    }
    delete[] date;
    
    string str_temp2;
    stringstream tt;
    char* time;
    tt << HH << ":" << MN << ":" << SS;
    tt >> str_temp2;
    time = new char[str_temp2.size() + 1];
    strcpy( time, str_temp2.c_str() );
    
    if( fits_update_key( fptr, TSTRING, ( char* )"TIME", time, ( char* )"FITS file creation time (hh:mm:ss)", &status ) )
    {
        return printerror( status );
    }
    delete[] time;
    
    char iCreator[100];
    sprintf( iCreator, "EVENTDISPLAY_%s", fEVNDISPVersion.c_str() );
    if( fits_update_key( fptr, TSTRING, ( char* )"CREATOR", iCreator, ( char* )"Software package and version creating file", &status ) )
    {
        return printerror( status );
    }
    
    char iVersion[] = "0.2";
    if( fits_update_key( fptr, TSTRING, ( char* )"VERSION", iVersion, ( char* )"VERITAS FITS standard version", &status ) )
    {
        return printerror( status );
    }
    
    return true;
}



/*************************************/
/* Create FITS file at the End       */
/*************************************/

bool VFITS::writeFITSFile( bool iPrint )
{
    if( fWriteOneFile == true )
    {
        int status = 0;
        if( fits_close_file( fptr, &status ) )
        {
            return printerror( status ) ;
        }
        if( iPrint )
        {
            cout << "Close FITS file " << endl;
        }
    }
    cout << "END of VFITS" << endl;
    
    return true;
}

//*****************METHODES UNDER DEVELOPEMENT*****************

/***************************************/
/* Merge columns of FITS BinTables     */
/***************************************/
bool VFITS::mergeColumns( fitsfile* fPtr, vector<int> hdunums, vector<vector <int> > columns, int nRows, bool iPrint )
{
    if( !fWriteOneFile )
    {
        cout << "!!!!Merging of tables only possible if they are all in the same FITS file! Use 'fWriteOneFile=true'" << endl;
        return false;
    }
    if( hdunums.size() != columns.size() )
    {
        cout << "!!!Number of HDUs is not the same as the Number of column_configuration_string ('columns' Parameter) " << endl;
        return false;
    }
    
    char* tType[10];
    char* tUnit[10];
    char* tForm[10];
    
    vector< double> rows;
    vector< vector<double> > table;
    int hdutype = 0;
    int status = 0;
    int anynul = 0;
    int t = 0;
    int p = 0;
    int d = 0;
    //string value;
    char value[3000];
    char comment[10];
    char tname[10];
    for( unsigned int i = 0; i < hdunums.size(); i++ )
    {
        if( fits_movabs_hdu( fPtr, hdunums[i] - t, &hdutype, &status ) )
        {
            return printerror( status ) ;
        }
        if( iPrint )
        {
            cout << "  At hdunums[" << i << "] - t = " << hdunums[i] - t << endl;
        }
        for( unsigned int j = 0; j < columns[i].size(); j++ )
        {
            if( iPrint )
            {
                cout << "   At columns[" << i << "][" << j << "] = " << columns[i][j] << endl;
            }
            sprintf( tname, "TUNIT%d", columns[i][j] );
            if( fits_read_key( fPtr, TSTRING, tname, value + d, comment, &status ) )
            {
                return printerror( status ) ;
            }
            tUnit[p] = value + d;
            d = d + 30;
            sprintf( tname, "TFORM%d", columns[i][j] );
            if( fits_read_key( fPtr, TSTRING, tname, value + d, comment, &status ) )
            {
                return printerror( status ) ;
            }
            tForm[p] = value + d;
            d = d + 30;
            sprintf( tname, "TTYPE%d", columns[i][j] );
            if( fits_read_key( fPtr, TSTRING, tname, value + d, comment, &status ) )
            {
                return printerror( status ) ;
            }
            tType[p] = value + d;
            d = d + 30;
            rows.clear();
            for( int k = 0; k < nRows; k++ )
            {
                double test = 0.;
                if( fits_read_col( fPtr, TDOUBLE, columns[i][j] , k + 1, 1, 1, 0, &test, &anynul, &status ) )
                {
                    return printerror( status );
                }
                if( iPrint )
                {
                    cout << "     OutPut of table at Row " << k << " and column " << j << " = " << test << endl;
                }
                rows.push_back( test );
            }
            table.push_back( rows );
            p++;
        }
        if( fits_delete_hdu( fPtr, &hdutype, &status ) )
        {
            return printerror( status );
        }
        t++;
    }
    //  cout<<"nRows = "<<nRows<<endl;
    //  cout<<"NROWS = "<<table[1].size()<<endl;
    // cout<<"NCOLS = "<<table.size()<<endl;
    int f = 2;
    for( int g = 0; g < p; g++ )
    {
        string type1( tType[g] );
        for( int h = g; h < p; h++ )
        {
            if( h == g )
            {
                continue;
            }
            string type2( tType[h] );
            if( type1 == type2 )
            {
                sprintf( tType[h], "%s_%d", tType[h], f++ );
            }
        }
    }
    
    
    vector< double> cols;
    vector< vector<double> > table2;
    for( unsigned int k = 0; k < table[1].size(); k++ )
    {
        cols.clear();
        for( unsigned int j = 0; j < table.size(); j++ )
        {
            if( iPrint )
            {
                cout << "     OutPut of Table2 at Row " << k << " and column " << j << " = " << table[j][k] << endl;
            }
            cols.push_back( table[j][k] );
        }
        table2.push_back( cols );
    }
    
    createTableFitsFile( table2, tType, tUnit, tForm, "test", iPrint );
    
    return true;
}


