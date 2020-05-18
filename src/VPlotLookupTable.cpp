/*

   VPlotLookupTable

   missing compared to v35x: plot slices, interpolation of lookup tables

*/

#include "VPlotLookupTable.h"

VPlotLookupTable::VPlotLookupTable()
{

    // list of valid lookup table names
    fListOfTableNames.push_back( "mscw" );
    fListOfTableNames.push_back( "mscl" );
    fListOfTableNames.push_back( "energyER" );
    fListOfTableNames.push_back( "energySR" );
    
    setPlottingLogEnergyAxis();
    setPlottingLogSizeAxis();
    setPlottingDistanceAxis();
    setCanvasSize( 450, 400 );
}

void VPlotLookupTable::printLookupTables()
{
    for( unsigned int i = 0; i < fLookupTableData.size(); i++ )
    {
        cout << fLookupTableData[i]->fLookupTable << endl;
        cout << "\t file: " << fLookupTableData[i]->fLookupTableFileName << endl;
        cout << "\t ze: " << fLookupTableData[i]->fZe;
        cout << ", az: " << fLookupTableData[i]->fAz;
        cout << ", woff: " << fLookupTableData[i]->fWobbleOffset;
        cout << ", telid: " << fLookupTableData[i]->fTelID;
        cout << ", pedvar: " << fLookupTableData[i]->fNoise << endl;
        cout << endl;
    }
}

void VPlotLookupTable::plot2DHistogram( TH2F* h, unsigned int iSetID, string ititle, int iCanvasX, double i_min, double i_max, bool iZLog )
{
    if( !h )
    {
        cout << "VPlotLookupTable::plot2DHistogram error: histogram not found " << endl;
        return;
    }
    
    char hname[600];
    char htitle[600];
    sprintf( hname, "c%s_%d_%s_%d_%d_%d_%d", h->GetName(), iSetID, fLookupTableData[iSetID]->fLookupTable.c_str(), fLookupTableData[iSetID]->fZe,
             fLookupTableData[iSetID]->fAz, fLookupTableData[iSetID]->fNoise, fLookupTableData[iSetID]->fWobbleOffset );
    sprintf( htitle, "lookup table (%s, %s, %d deg, az %d, noise %d, woff %.2f", ititle.c_str(),
             fLookupTableData[iSetID]->fLookupTable.c_str(), fLookupTableData[iSetID]->fZe,
             fLookupTableData[iSetID]->fAz, fLookupTableData[iSetID]->fNoise, 1. / 100.*( double )fLookupTableData[iSetID]->fWobbleOffset );
             
    TCanvas* cE = new TCanvas( hname, htitle, iCanvasX, 10, fPlottingCanvasX, fPlottingCanvasY );
    cE->SetGridx( 0 );
    cE->SetGridy( 0 );
    cE->SetLeftMargin( 0.13 );
    cE->SetRightMargin( 0.16 );
    if( iZLog && ( i_min < -998. || i_min > 0. ) && h->GetEntries() > 0. )
    {
        cE->SetLogz( 1 );
    }
    cE->Draw();
    
    h->SetStats( 0 );
    h->SetTitle( "" );
    h->GetYaxis()->SetTitleOffset( 1.4 );
    h->GetZaxis()->SetTitleOffset( 1.3 );
    if( fLookupTableData[iSetID]->fLookupTable == "energyER" )
    {
        h->SetAxisRange( fLogEnergyAxis_min, fLogEnergyAxis_max, "X" );
    }
    else
    {
        h->SetAxisRange( fLogSizeAxis_min, fLogSizeAxis_max, "X" );
    }
    h->SetAxisRange( fDistanceAxis_min, fDistanceAxis_max, "Y" );
    
    if( i_min > -998. )
    {
        h->SetMinimum( i_min );
    }
    if( i_max > -998. )
    {
        h->SetMaximum( i_max );
    }
    
    h->Draw( "colz" );
}

TH2F* VPlotLookupTable::divide2DHistograms( TH2F* h1, TH2F* h2, char* hname )
{
    if( !h1 || !h2 )
    {
        return 0;
    }
    
    TH2F* hMM = new TH2F( hname, "", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(),
                          h1->GetNbinsY(), h1->GetYaxis()->GetXmin(), h1->GetYaxis()->GetXmax() );
    hMM->SetXTitle( h1->GetXaxis()->GetTitle() );
    hMM->SetYTitle( h1->GetYaxis()->GetTitle() );
    hMM->SetZTitle( h1->GetZaxis()->GetTitle() );
    
    int ii = 0;
    int jj = 0;
    
    for( int i = 1; i <= hMM->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= hMM->GetNbinsY(); j++ )
        {
            ii = h2->GetXaxis()->FindBin( h1->GetXaxis()->GetBinCenter( i ) );
            jj = h2->GetYaxis()->FindBin( h1->GetYaxis()->GetBinCenter( j ) );
            if( h2->GetBinContent( ii, jj ) > 0. )
            {
                hMM->SetBinContent( i, j, h1->GetBinContent( i, j ) / h2->GetBinContent( ii, jj ) );
            }
            else
            {
                hMM->SetBinContent( i, j, 0. );
            }
        }
    }
    
    return hMM;
}


void VPlotLookupTable::plotLookupTables( unsigned int iSetID, double i_ymin, bool iMedianOnly )
{
    if( iSetID >= fLookupTableData.size() )
    {
        cout << "VPlotLookupTable::plotLookupTables error: setID to large (should be <";
        cout << fLookupTableData.size() << ")" << endl;
        return;
    }
    plot2DHistogram( fLookupTableData[iSetID]->hmedian, iSetID, "median", 10, i_ymin, -999., ( fLookupTableData[iSetID]->fLookupTable == "energySR" ) );
    if( iMedianOnly )
    {
        return;
    }
    plot2DHistogram( fLookupTableData[iSetID]->hmean,   iSetID, "mean", 100, i_ymin, -999., ( fLookupTableData[iSetID]->fLookupTable == "energySR" ) );
    plot2DHistogram( fLookupTableData[iSetID]->hmpv,    iSetID, "mpv", 200, i_ymin, -999., ( fLookupTableData[iSetID]->fLookupTable == "energySR" ) );
    plot2DHistogram( fLookupTableData[iSetID]->hsigma,  iSetID, "sigma", 300, i_ymin, -999., false );
    plot2DHistogram( fLookupTableData[iSetID]->hnevents, iSetID, "number of events", 400, i_ymin, -999., true );
    
    // divide mean by median
    if( fLookupTableData[iSetID]->hmedian && fLookupTableData[iSetID]->hmean )
    {
        char hname[200];
        sprintf( hname, "hMM_%d", iSetID );
        
        TH2F* hMM = divide2DHistograms( fLookupTableData[iSetID]->hmean, fLookupTableData[iSetID]->hmedian, hname );
        if( hMM )
        {
            hMM->SetZTitle( "mean / median" );
        }
        
        plot2DHistogram( hMM, iSetID, "mean / median", 200, 0.95, 1.05 );
    }
    // divide MPV by median
    if( fLookupTableData[iSetID]->hmedian && fLookupTableData[iSetID]->hmpv )
    {
        char hname[200];
        sprintf( hname, "hMP_%d", iSetID );
        
        TH2F* hMP = divide2DHistograms( fLookupTableData[iSetID]->hmpv, fLookupTableData[iSetID]->hmedian, hname );
        if( hMP )
        {
            hMP->SetZTitle( "mpv/ median" );
        }
        
        plot2DHistogram( hMP, iSetID, "mpv / median", 200, 0.95, 1.05 );
    }
    
    
}

void VPlotLookupTable::plotRelativeTables( unsigned int iSetID1, unsigned int iSetID2, double iMin, double iMax )
{
    if( iSetID1 >= fLookupTableData.size() || iSetID2 >= fLookupTableData.size() )
    {
        cout << "VPlotLookupTable::plotLookupTables error: setID to large (should be <" << fLookupTableData.size() << endl;
        return;
    }
    if( fLookupTableData[iSetID1]->fLookupTable != fLookupTableData[iSetID2]->fLookupTable )
    {
        cout << "VPlotLookupTable::plotRelativeTables: error, trying to divide different histograms" << endl;
        cout << "\t" << fLookupTableData[iSetID1]->fLookupTable << "\t" << fLookupTableData[iSetID2]->fLookupTable << endl;
        return;
    }
    
    // median
    char hname[200];
    sprintf( hname, "hMM_DIV_median_%d_%d", iSetID1, iSetID2 );
    
    TH2F* hMM_median = divide2DHistograms( fLookupTableData[iSetID1]->hmedian, fLookupTableData[iSetID2]->hmedian, hname );
    sprintf( hname, "%s (%d/%d)", hMM_median->GetZaxis()->GetTitle(), iSetID1, iSetID2 );
    if( hMM_median )
    {
        hMM_median->SetZTitle( hname );
    }
    
    sprintf( hname, "median (%d/%d)", iSetID1, iSetID2 );
    plot2DHistogram( hMM_median, iSetID1, hname, 500, iMin, iMax );
    
    // sigma
    sprintf( hname, "hMM_DIV_sigma_%d_%d", iSetID1, iSetID2 );
    
    TH2F* hMM_sigma = divide2DHistograms( fLookupTableData[iSetID1]->hsigma, fLookupTableData[iSetID2]->hsigma, hname );
    if( hMM_sigma )
    {
        sprintf( hname, "%s (%d/%d)", hMM_sigma->GetZaxis()->GetTitle(), iSetID1, iSetID2 );
        hMM_sigma->SetZTitle( hname );
        sprintf( hname, "sigma (%d/%d)", iSetID1, iSetID2 );
        plot2DHistogram( hMM_sigma, iSetID1, hname, 550, iMin, iMax );
    }
    
    
    // mean
    sprintf( hname, "hMM_DIV_mean_%d_%d", iSetID1, iSetID2 );
    
    TH2F* hMM_mean = divide2DHistograms( fLookupTableData[iSetID1]->hmean, fLookupTableData[iSetID2]->hmean, hname );
    if( hMM_mean )
    {
        sprintf( hname, "%s (%d/%d)", hMM_mean->GetZaxis()->GetTitle(), iSetID1, iSetID2 );
        hMM_mean->SetZTitle( hname );
        sprintf( hname, "mean (%d/%d)", iSetID1, iSetID2 );
        plot2DHistogram( hMM_mean, iSetID1, hname, 600, iMin, iMax );
    }
    
    
}

/*
 * add a lookup table set to the list of lookup tables
 *
 */
bool VPlotLookupTable::addLookupTable( string iLookupTableFile, string iTable, int ze, int az, int telID, int noise, int woff )
{
    // add a new data set
    if( !checkTableName( iTable ) )
    {
        return false;
    }
    
    TFile* fI = new TFile( iLookupTableFile.c_str() );
    if( fI->IsZombie() )
    {
        cout << "VPlotLookupTable::addLookupTable error: file not found: " << iLookupTableFile << endl;
        return false;
    }
    char hname[600];
    // create full directory name
    sprintf( hname, "tel_%d/NOISE_%05d/ze_%03d/woff_%04d/az_%d/", telID, noise, ze * 10, woff, az );
    if( !fI->cd( hname ) )
    {
        cout << "VPlotLookupTable::addLookupTable error: directory  " << hname << " not found" << endl;
        return false;
    }
    cout << "directory: " << hname << endl;
    //  variable directory
    sprintf( hname, "%s", iTable.c_str() );
    if( !gDirectory->cd( hname ) )
    {
        cout << "VPlotLookupTable::addLookupTable error: directory for variable " << iTable << " not found" << endl;
        return false;
    }
    
    // now everything seems to be fine, fill new data entry
    fLookupTableData.push_back( new VPlotLookupTableData() );
    fLookupTableData.back()->fLookupTable = iTable;
    fLookupTableData.back()->fLookupTableFileName = iLookupTableFile;
    fLookupTableData.back()->fLookupTableFile = fI;
    fLookupTableData.back()->fZe = ze;
    fLookupTableData.back()->fAz = az;
    fLookupTableData.back()->fWobbleOffset = woff;
    fLookupTableData.back()->fTelID = telID;
    
    if( iTable == "mscw" )
    {
        fLookupTableData.back()->hmedian = ( TH2F* )gDirectory->Get( "width_median_tb" );
        fLookupTableData.back()->hmean = ( TH2F* )gDirectory->Get( "width_mean_tb" );
        fLookupTableData.back()->hmpv = ( TH2F* )gDirectory->Get( "width_mpv_tb" );
        fLookupTableData.back()->hsigma = ( TH2F* )gDirectory->Get( "width_sigma_tb" );
        fLookupTableData.back()->hnevents = ( TH2F* )gDirectory->Get( "width_nevents_tb" );
    }
    else if( iTable == "mscl" )
    {
        fLookupTableData.back()->hmedian = ( TH2F* )gDirectory->Get( "length_median_tb" );
        fLookupTableData.back()->hmean = ( TH2F* )gDirectory->Get( "length_mean_tb" );
        fLookupTableData.back()->hmpv = ( TH2F* )gDirectory->Get( "length_mpv_tb" );
        fLookupTableData.back()->hsigma = ( TH2F* )gDirectory->Get( "length_sigma_tb" );
        fLookupTableData.back()->hnevents = ( TH2F* )gDirectory->Get( "length_nevents_tb" );
    }
    else if( iTable == "energyER" )
    {
        fLookupTableData.back()->hmedian = ( TH2F* )gDirectory->Get( "hMedian_energy_tb" );
        fLookupTableData.back()->hsigma = ( TH2F* )gDirectory->Get( "hSigma_energy_tb" );
        fLookupTableData.back()->hmpv = 0;
        fLookupTableData.back()->hnevents = ( TH2F* )gDirectory->Get( "hNevents_energy_tb" );
        fLookupTableData.back()->hmean = ( TH2F* )gDirectory->Get( "hMean_energy_tb" );
    }
    else if( iTable == "energySR" )
    {
        fLookupTableData.back()->hmedian = ( TH2F* )gDirectory->Get( "energySR_median_tb" );
        fLookupTableData.back()->hmean = ( TH2F* )gDirectory->Get( "energySR_mean_tb" );
        fLookupTableData.back()->hmpv = ( TH2F* )gDirectory->Get( "energySR_mpv_tb" );
        fLookupTableData.back()->hsigma = ( TH2F* )gDirectory->Get( "energySR_sigma_tb" );
        fLookupTableData.back()->hnevents = ( TH2F* )gDirectory->Get( "energySR_nevents_tb" );
    }
    if( fLookupTableData.back()->hmedian && fLookupTableData.back()->hmean
            && fLookupTableData.back()->hsigma && fLookupTableData.back()->hnevents )
    {
        cout << "all histograms found..." << endl;
    }
    
    return true;
}

bool VPlotLookupTable::checkTableName( string iTableName )
{
    for( unsigned int i = 0; i < fListOfTableNames.size(); i++ )
    {
        if( iTableName == fListOfTableNames[i] )
        {
            return true;
        }
    }
    
    cout << "VPlotLookupTable::checkTableName unkown table name: " << iTableName << endl;
    cout << "\t allowed are: ";
    for( unsigned int i = 0; i < fListOfTableNames.size(); i++ )
    {
        cout << "\t" << fListOfTableNames[i];
    }
    cout << endl;
    cout << endl;
    return false;
    
}

/*
 * plot slices of the given lokup tables
 *
 */

void VPlotLookupTable::plotLookupTableSlice( double iLogSize, double iR )
{
    // prepare canvas
    char hname[600];
    char htitle[600];
    sprintf( hname, "cSlice_%d_%d", ( int )iLogSize * 100, ( int )iR );
    sprintf( htitle, "lookup table slice (log_{10} %.1f, R %.1f)", iLogSize, iR );
    TCanvas* cS = new TCanvas( hname, htitle, 500, 10, fPlottingCanvasX, fPlottingCanvasY );
    cS->SetGridx( 0 );
    cS->SetGridy( 0 );
    cS->Draw();
    
    // loop over all slices and plot them
    //
    for( unsigned int i = 0; i < fLookupTableData.size(); i++ )
    {
        if( fLookupTableData[i] && fLookupTableData[i]->hmedian )
        {
            sprintf( hname, "%s_%d_%d_%u", fLookupTableData[i]->hmedian->GetName(), ( int )iLogSize * 100, ( int )iR, i );
            TH1D* h = 0;
            if( iLogSize > 0.01 )
            {
                h = fLookupTableData[i]->hmedian->ProjectionY( hname,
                        fLookupTableData[i]->hmedian->GetXaxis()->FindBin( iLogSize ),
                        fLookupTableData[i]->hmedian->GetXaxis()->FindBin( iLogSize ) );
            }
            else
            {
                h = fLookupTableData[i]->hmedian->ProjectionX( hname,
                        fLookupTableData[i]->hmedian->GetYaxis()->FindBin( iR ),
                        fLookupTableData[i]->hmedian->GetYaxis()->FindBin( iR ) );
            }
            if( h )
            {
                h->SetLineColor( i + 1 );
                h->SetStats( 0 );
                h->GetYaxis()->SetTitle( fLookupTableData[i]->hmedian->GetZaxis()->GetTitle() );
                if( i == 0 )
                {
                    h->Draw();
                }
                else
                {
                    h->Draw( "same" );
                }
            }
        }
    }
    
}

/*
 * interpolate or smooth lookup tables
 */
bool VPlotLookupTable::smoothLookupTables( unsigned int iSetID, string iMethod, int iMinEvents )
{
    if( iSetID >= fLookupTableData.size() )
    {
        cout << "VPlotLookupTable::smoothLookupTables error: setID to large (should be <";
        cout << fLookupTableData.size() << ")" << endl;
        return false;
    }
    VInterpolate2DHistos i_inter;
    
    if( iMethod == "interpolate" )
    {
        fLookupTableData[iSetID]->hmedian =
            i_inter.doSimpleInterpolation( fLookupTableData[iSetID]->hmedian,
                                           iMethod,
                                           2, 2, false,
                                           fLookupTableData[iSetID]->hnevents,
                                           iMinEvents );
    }
    else if( iMethod.find( "fit" ) != string::npos )
    {
        fLookupTableData[iSetID]->hmedian =
            i_inter.doLogLinearExtrapolation( fLookupTableData[iSetID]->hmedian,
                                              iMethod,
                                              fLookupTableData[iSetID]->hnevents,
                                              iMinEvents );
                                              
    }
    else
    {
        return false;
    }
    
    return true;
    
}

// ====================================================================================

VPlotLookupTableData::VPlotLookupTableData()
{
    fLookupTable = "";
    fLookupTableFileName = "";
    fLookupTableFile = 0;
    fZe = 0;
    fAz = 0;
    fTelID = 0;
    fNoise = 0;
    fWobbleOffset = 0;
    
    hmedian = 0;
    hmean = 0;
    hmpv = 0;
    hsigma = 0;
    hnevents = 0;
}

