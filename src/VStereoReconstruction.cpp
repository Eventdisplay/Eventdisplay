/*! \file VStereoReconstruction.cpp
    \brief calculate and plot array reconstruction quality


*/

#include "VStereoReconstruction.h"

VStereoReconstructionData::VStereoReconstructionData()
{
    fTitle = "notitle";
    fColor = 1;
    fMarkerStyle = 20;
    fMarkerSize = 1.;
    fLineStyle = 1;
    fLineWidth = 1;
    fFillStyle = 1001;
    fPlotStyle = "cp";
    
    gData = new TGraphAsymmErrors( 1 );
}

void VStereoReconstructionData::draw()
{
    if( !gData )
    {
        return;
    }
    
    gData->SetLineColor( fColor );
    gData->SetMarkerColor( fColor );
    gData->SetFillColor( fColor );
    gData->SetLineWidth( fLineWidth );
    gData->SetLineStyle( fLineStyle );
    gData->SetMarkerStyle( fMarkerStyle );
    gData->SetMarkerSize( fMarkerSize );
    gData->SetFillStyle( fFillStyle );
    
    gData->Draw( fPlotStyle.c_str() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


VStereoReconstruction::VStereoReconstruction()
{
    bDebug = false;
    
    setName();
    setPlottingVariable( "direction" );
    setPlottingEnergyRangeLinear();
    setPlottingLogEnergyAxis();
    setPlottingYaxis();
    setPlottingCanvasSize();
}

/*
   todo: check if this is valid
*/
bool VStereoReconstruction::setPlottingVariable( string iVar )
{
    fPlottingVariable = iVar;
    
    return true;
}

TCanvas* VStereoReconstruction::plot( TCanvas* c )
{
    /*   if( fData.size() == 0 )
       {
          cout << "VStereoReconstruction::plot error: empty data vector" << endl;
          return 0;
       } */
    
    char hname[800];
    char htitle[800];
    
    TH1D* hNull = 0;
    
    if( c == 0 )
    {
        sprintf( hname, "c_%s_%s",  fName.c_str(), fPlottingVariable.c_str() );
        sprintf( htitle, "%s reconstruction (%s)", fPlottingVariable.c_str(), fName.c_str() );
        c = new TCanvas( hname, htitle, 10, 10, fPlottingCanvasSizeX, fPlottingCanvasSizeY );
        c->SetGridx( 0 );
        c->SetGridy( 0 );
        
        c->SetLeftMargin( 0.13 );
        
        sprintf( hname, "hnull_%s", fName.c_str() );
        hNull = new TH1D( hname, "", 100, log10( fPlottingMinEnergy ), log10( fPlottingMaxEnergy ) );
        hNull->SetMinimum( fPlottingYaxisMin );
        hNull->SetMaximum( fPlottingYaxisMax );
        hNull->SetStats( 0 );
        hNull->SetXTitle( "log_{10} energy [TeV]" );
        hNull->GetYaxis()->SetTitleOffset( 1.5 );
        if( fPlottingVariable == "direction" )
        {
            hNull->SetYTitle( "angular resolution [deg] (68%)" );
        }
        plot_nullHistogram( c, hNull, fPlottingLogEnergyAxis, false, hNull->GetYaxis()->GetTitleOffset(), fPlottingMinEnergy, fPlottingMaxEnergy );
    }
    fPlottingCanvas = c;
    
    fPlottingCanvas->cd();
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        fData[i]->draw();
    }
    
    return c;
    
    
}

bool VStereoReconstruction::setPlottingAtt( unsigned int iDataSet, string iPlotStyle, int iColor, int iMarkerStyle, float iMarkerSize, int iLineStyle, float iLineWidth, int iFillStyle )
{
    if( iDataSet >= fData.size() )
    {
        return false;
    }
    
    fData[iDataSet]->fColor = iColor;
    fData[iDataSet]->fMarkerStyle = iMarkerStyle;
    fData[iDataSet]->fMarkerSize = iMarkerSize;
    fData[iDataSet]->fLineStyle = iLineStyle;
    fData[iDataSet]->fLineWidth = ( Width_t )iLineWidth;
    fData[iDataSet]->fFillStyle = iFillStyle;
    fData[iDataSet]->fPlotStyle = iPlotStyle;
    
    return true;
}

/*
   expect root file: .root
   expect ascii file: .txt or .dat
*/
bool VStereoReconstruction::addDataSet( string iTitle, string iFile, double emin_lin, double emax_lin, bool bEnergyAxis_linear_GeV, bool bResolutionAxis_arcmin )
{
    unsigned int i_oldDataSize = fData.size();
    
    /////////////////////////////////////////////////////////////////////////
    // add a root file
    if( iFile.size() > 4 && iFile.substr( iFile.size() - 4, 4 ) == "root" )
    {
        TFile f( iFile.c_str() );
        if( f.IsZombie() )
        {
            cout << "VStereoReconstruction::addDataSet error reading root file: " << iFile << endl;
            return false;
        }
        if( fPlottingVariable == "direction" )
        {
            TGraphErrors* gE0 = 0;
            TH1F* AngRes = 0;
            if( f.Get( "direction" ) )
            {
                f.cd( "direction" );
                gE0 = ( TGraphErrors* )gDirectory->Get( "gAngE0_All" ); // from eventdisplay stereo_reconstruction
            }
            else
            {
                AngRes = ( TH1F* )f.Get( "AngRes" );
            }
            
            // add a eventdisplay style root file
            if( gE0 && gE0->GetN() > 0 )
            {
                cout << "Adding a eventdisplay style data set: " << iFile << endl;
                fData.push_back( new VStereoReconstructionData() );
                double x, y;
                int z = 0;
                for( int i = 0; i < gE0->GetN(); i++ )
                {
                    gE0->GetPoint( i, x, y );
                    if( emin_lin > 0. && log10( emin_lin ) > x )
                    {
                        continue;
                    }
                    if( emax_lin > 0. && log10( emax_lin ) < x )
                    {
                        continue;
                    }
                    // accept only 4 sigma points (preliminary)
                    if( gE0->GetErrorY( i ) > 0. && y / gE0->GetErrorY( i ) < 4. )
                    {
                        continue;
                    }
                    
                    fData.back()->gData->SetPoint( z, x, y );
                    fData.back()->gData->SetPointEYhigh( z, gE0->GetErrorY( i ) );
                    fData.back()->gData->SetPointEYlow( z, gE0->GetErrorY( i ) );
                    z++;
                }
                if( bDebug )
                {
                    fData.back()->gData->Print();
                }
            }
            // add a CTA production style root file
            else if( AngRes )
            {
                cout << "Adding a CTA production style data set " << iFile << endl;
                fData.push_back( new VStereoReconstructionData() );
                for( int i = 1; i <= AngRes->GetNbinsX(); i++ )
                {
                    if( emin_lin > 0. && log10( emin_lin ) > AngRes->GetXaxis()->GetBinCenter( i ) )
                    {
                        continue;
                    }
                    if( emax_lin > 0. && log10( emax_lin ) < AngRes->GetXaxis()->GetBinCenter( i ) )
                    {
                        continue;
                    }
                    fData.back()->gData->SetPoint( i - 1, AngRes->GetXaxis()->GetBinCenter( i ), AngRes->GetBinContent( i ) );
                    //		fData.back()->gData->SetPointEYhigh( i, AngRes->GetBinError( i ) );
                    //		fData.back()->gData->SetPointEYlow( i, AngRes->GetBinError( i ) );
                    fData.back()->gData->SetPointEYhigh( i - 1, 0. );
                    fData.back()->gData->SetPointEYlow( i - 1, 0. );
                }
                if( bDebug )
                {
                    fData.back()->gData->Print();
                }
            }
            else
            {
                cout << "VStereoReconstruction::addDataSet: no data found in " << iFile << endl;
                return false;
            }
        }
    }
    // assume that everything else is a text file
    else
    {
        TGraph g( iFile.c_str() );
        if( g.GetN() > 0 )
        {
            cout << "Adding a data set from ASCII: " << iFile << endl;
            fData.push_back( new VStereoReconstructionData() );
            double x, y;
            for( int i = 0; i < g.GetN(); i++ )
            {
                g.GetPoint( i, x, y );
                if( bEnergyAxis_linear_GeV )
                {
                    x = log10( x ) - 3.;
                }
                if( bResolutionAxis_arcmin )
                {
                    y /= 60.;
                }
                
                if( emin_lin > 0. && log10( emin_lin ) > x )
                {
                    continue;
                }
                if( emax_lin > 0. && log10( emax_lin ) < x )
                {
                    continue;
                }
                
                fData.back()->gData->SetPoint( i, x, y );
                fData.back()->gData->SetPointEYhigh( i, 0. );
                fData.back()->gData->SetPointEYlow( i, 0. );
            }
            if( bDebug )
            {
                fData.back()->gData->Print();
            }
        }
    }
    // do some automatic coloring/title, etc.
    if( i_oldDataSize < fData.size() )
    {
        fData.back()->fTitle = iTitle;
        fData.back()->fColor = fData.size();
        fData.back()->fMarkerStyle = 20 + fData.size() - 1;
    }
    
    return true;
}

bool VStereoReconstruction::removeDataSet( unsigned int iDataSet )
{
    if( fData.size() < iDataSet )
    {
        fData.erase( fData.begin() + iDataSet );
    }
    return false;
}

void VStereoReconstruction::plotLegend()
{
    if( fData.size() < 1 )
    {
        return;
    }
    
    TLegend* iLAng = 0;
    if( fData.size() < 3 )
    {
        iLAng = new TLegend( 0.5, 0.7, 0.94, 0.94 );
    }
    else
    {
        iLAng = new TLegend( 0.6, 0.6, 0.85, 0.85 );
    }
    
    for( unsigned int i = 0; i < fData.size(); i++ )
    {
        iLAng->AddEntry( fData[i]->gData, fData[i]->fTitle.c_str(), "pl" );
    }
    iLAng->Draw();
}

bool VStereoReconstruction::readDataSetsfromTextFile( string ifile, unsigned int iSet, bool bClearExistingDataSet )
{
    if( bClearExistingDataSet )
    {
        fData.clear();
    }
    
    // temporary data set
    string iFile;
    string iName;
    int     iColor;
    int     iMarkerStyle;
    float   iMarkerSize;
    int     iLineStyle;
    Width_t iLineWidth;
    int     iFillStyle = 1001;
    string  iPlotStyle;
    double emin = -1.;
    double emax = -1.;
    bool bEnergyAxis_linear_GeV = false;
    bool bResolutionAxis_arcmin = false;
    
    // open text file
    ifstream is;
    is.open( ifile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error opening input file list: " << ifile << endl;
        return false;
    }
    cout << "reading input file list for set " << iSet << " from " << ifile << endl;
    string is_line;
    string temp;
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
            is_stream >> temp;
            // check set number
            if( ( unsigned int )( atoi( temp.c_str() ) ) != iSet )
            {
                continue;
            }
            
            // read file name
            is_stream >> iFile;
            is_stream >> temp;
            emin = atof( temp.c_str() );
            is_stream >> temp;
            emax = atoi( temp.c_str() );
            is_stream >> temp;
            bEnergyAxis_linear_GeV = atoi( temp.c_str() );
            is_stream >> temp;
            bResolutionAxis_arcmin = atoi( temp.c_str() );
            is_stream >> iPlotStyle;
            is_stream >> temp;
            iColor = atoi( temp.c_str() );
            is_stream >> temp;
            iMarkerStyle = atoi( temp.c_str() );
            is_stream >> temp;
            iMarkerSize = atof( temp.c_str() );
            is_stream >> temp;
            iLineStyle = atoi( temp.c_str() );
            is_stream >> temp;
            iLineWidth = ( Width_t )( atof( temp.c_str() ) );
            iName = is_stream.str().substr( is_stream.tellg(), is_stream.str().size() );
            
            if( addDataSet( iName, iFile, emin, emax, bEnergyAxis_linear_GeV, bResolutionAxis_arcmin ) )
            {
                setPlottingAtt( fData.size() - 1, iPlotStyle, iColor, iMarkerStyle, iMarkerSize, iLineStyle, iLineWidth, iFillStyle );
            }
        }
    }
    
    return true;
}
