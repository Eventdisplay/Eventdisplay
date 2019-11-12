/*! \class VFitEffectiveArea
    \brief effective area fitting


*/

#include "VFitEffectiveArea.h"

VFitEffectiveArea::VFitEffectiveArea( string ifile )
{
    fInputFile = ifile;
    
    bZombie = false;
    
    gEff = 0;
    gEffLog = 0;
    cEffLog = 0;
    cEffLin = 0;
    bAMC = false;
    
    fCanvasMaximum = 1.e5;
    
    // his list for output tree
    hList = new TList();
    
    // fit function (on log Aeff scale)
    fEff = new TF1( "fEffLog", fun_eff, -1.5, 2.0, 7 );
    fEff->SetLineColor( 2 );
    // fit function (on lin Aeff scale)
    fEffLin = new TF1( "fEffLin", fun_effLin, -1.5, 2.0, 7 );
    fEffLin->SetLineColor( 2 );
    hList->Add( fEffLin );
    resetFitParameters();
    
    fZe = 0.;
    fAz = 0;
    fAzMin = 0.;
    fAzMax = 0.;
    fIndex = 0.;
    fFitMin = 0.;
    fFitMax = 0.;
    fTree = new TTree( "EffFit", "fit results for effective areas" );
    fTree->Branch( "Ze", &fZe, "Ze/D" );
    fTree->Branch( "Az", &fAz, "Az/I" );
    fTree->Branch( "AzMin", &fAzMin, "AzMin/D" );
    fTree->Branch( "AzMax", &fAzMax, "AzMax/D" );
    fTree->Branch( "Index", &fIndex, "Index/D" );
    fTree->Branch( "AMC", &fAMC, "AMC/I" );
    fTree->Branch( hList, 64000, 1 );
    fTree->Branch( "Fitxmin", &fFitMin, "Fitxmin/D" );
    fTree->Branch( "Fitxmax", &fFitMax, "Fitxmax/D" );
}


void VFitEffectiveArea::fill()
{
    if( fTree )
    {
        fTree->Fill();
    }
}


void VFitEffectiveArea::fit( double xmin, double xmax )
{
    if( !fEff || !gEff || !gEffLog )
    {
        return;
    }
    
    fFitMin = xmin;
    fFitMax = xmax;
    
    makeCanvases();
    
    if( cEffLog )
    {
        cEffLog->cd();
    }
    else
    {
        return;
    }
    
    gEffLog->Draw( "p" );
    
    gEffLog->Fit( fEff, "em", "", xmin, xmax );
    
    cEffLog->Update();
    
    if( cEffLin )
    {
        cEffLin->cd();
    }
    else
    {
        return;
    }
    setLinearFitParameters();
    
    gEff->Draw( "p" );
    fEffLin->DrawCopy( "same" );
}


void VFitEffectiveArea::makeCanvases()
{
    char hname[1000];
    
    gStyle->SetStatX( 0.8 );
    gStyle->SetStatY( 0.4 );
    
    cEffLog = new TCanvas( fName.c_str(), fTitle.c_str(), 10, 10, 600, 400 );
    cEffLog->SetGridx( 0 );
    cEffLog->SetGridy( 0 );
    cEffLog->Draw();
    
    sprintf( hname, "hnullLog_%s", fName.c_str() );
    hnullLog = new TH1D( hname, "", 100, -1.5, 2.0 );
    hnullLog->SetStats( 0 );
    if( bAMC )
    {
        hnullLog->SetXTitle( "log_{10} energy_{MC} [TeV]" );
    }
    else
    {
        hnullLog->SetXTitle( "log_{10} energy_{Rec} [TeV]" );
    }
    hnullLog->SetYTitle( "log_{10} effective area [m^{2}]" );
    hnullLog->SetMinimum( 0. );
    hnullLog->SetMaximum( 7. );
    hnullLog->Draw();
    
    sprintf( hname, "cLin_%s", fName.c_str() );
    cEffLin = new TCanvas( hname, fTitle.c_str(), 650, 10, 600, 400 );
    cEffLin->SetGridx( 0 );
    cEffLin->SetGridy( 0 );
    cEffLin->Draw();
    
    sprintf( hname, "hnullin_%s", fName.c_str() );
    hnullLin = new TH1D( hname, "", 100, -1.5, 2.0 );
    hnullLin->SetStats( 0 );
    if( bAMC )
    {
        hnullLin->SetXTitle( "log_{10} energy_{MC} [TeV]" );
    }
    else
    {
        hnullLin->SetXTitle( "log_{10} energy_{Rec} [TeV]" );
    }
    hnullLin->SetYTitle( "effective area [m^{2}]" );
    hnullLin->SetMinimum( 1. );
    hnullLin->SetMaximum( fCanvasMaximum );
    hnullLin->Draw();
    
}


void VFitEffectiveArea::resetFitParameters()
{
    if( !fEff )
    {
        return;
    }
    
    // set start parameters for fit
    fEff->SetParameter( 0, 3.5e4 );
    fEff->SetParameter( 1, 25. );
    // for two telescopes and Arec
    //       fEff->SetParLimits( 1, 15.0, 50.5 );
    //       fEff->SetParLimits( 1, 8.0, 50.5 );
    fEff->SetParLimits( 1, 5.0, 50.5 );
    fEff->SetParameter( 2, 1.2 );
    fEff->SetParLimits( 2, 0.01, 2.0 );
    fEff->SetParameter( 3, -0.5 );
    fEff->SetParameter( 4, 1.361931e+05 );
    fEff->SetParameter( 5, -6.917310e+04 );
    fEff->SetParameter( 6, 0.4 );
    // often the fit reaches that limit, but making it smaller worsens the fit
    //       fEff->SetParLimits( 6, -0.2, 0.7 );
    // for two telescopes and Arec
    //       fEff->SetParLimits( 6, -0.7, 0.6+i*0.1 );
    // for three telescopes and Arec
    fEff->SetParLimits( 6, -0.7, 0. );
    
    setLinearFitParameters();
    
}


void VFitEffectiveArea::setLinearFitParameters()
{
    if( !fEff || !fEffLin )
    {
        return;
    }
    
    // set parameters for linear Aeff axis
    fEffLin->SetParameter( 0, fEff->GetParameter( 0 ) );
    fEffLin->SetParameter( 1, fEff->GetParameter( 1 ) );
    fEffLin->SetParameter( 2, fEff->GetParameter( 2 ) );
    fEffLin->SetParameter( 3, fEff->GetParameter( 3 ) );
    fEffLin->SetParameter( 4, fEff->GetParameter( 4 ) );
    fEffLin->SetParameter( 5, fEff->GetParameter( 5 ) );
    fEffLin->SetParameter( 6, fEff->GetParameter( 6 ) );
}


/*!
    read effective areas for a certain zenith angle from file

    bAMC = true: use effective area vs true MC energy
    bAMC = false: use effective area vs reconstructed energy
*/
void VFitEffectiveArea::getEffectiveArea( double iZe, int az, double index, bool iAMC )
{
    bAMC = iAMC;
    fAMC = int( bAMC );
    TDirectory* iDir = gDirectory;
    TFile fIs( fInputFile.c_str() );
    if( fIs.IsZombie() )
    {
        cout << "error opening file" << endl;
        bZombie = true;
        return;
    }
    
    TTree* ts = ( TTree* )gDirectory->Get( "fEffArea" );
    if( !ts )
    {
        cout << "error finding effective area tree" << endl;
        return;
    }
    CEffArea c( ts );
    
    char hname[1000];
    sprintf( hname, "c_%d_%d_%d_%d", ( int )iZe, az, ( int )( index * 100 ), ( int )bAMC );
    fName = hname;
    sprintf( hname, "%d_%d_%d_%d", ( int )iZe, az, ( int )( index * 100 ), ( int )bAMC );
    fTitle = hname;
    
    int nentries = c.fChain->GetEntries();
    
    for( int i = 0; i < nentries; i++ )
    {
        c.GetEntry( i );
        
        if( fabs( iZe - c.ze ) < 0.1 && az == c.az && fabs( index - c.index ) < 0.05 )
        {
            if( bAMC )
            {
                gEff = ( TGraphAsymmErrors* )c.gEffAreaMC->Clone();
            }
            else
            {
                gEff = ( TGraphAsymmErrors* )c.gEffAreaRec->Clone();
            }
            iDir->Append( gEff );
            cout << "found effective area: " << gEff->GetName() << endl;
            setMarkers( gEff, 20, 2, 1 );
            hList->Add( gEff );
            getLogGraph();
            setMarkers( gEffLog, 20, 2, 1 );
            hList->Add( gEffLog );
            
            fZe = iZe;
            fAz = az;
            fIndex = index;
            fAzMin = c.azMin;
            fAzMax = c.azMax;
        }
    }
    fIs.Close();
    return;
}


void VFitEffectiveArea::setMarkers( TGraph* g, int iMarker, double iSize, int iColor )
{
    if( !g )
    {
        return;
    }
    
    g->SetMarkerSize( iSize );
    g->SetMarkerColor( iColor );
    g->SetLineColor( iColor );
    g->SetMarkerStyle( iMarker );
}


/*!
    fit function for effective areas

    following APP 15 (2001) 335 with additional pol3 for turnover at high energies

    transition between two function is fit parameter p[6]

*/
double fun_eff( double* e, double* p )
{
    double x = e[0];
    
    double f = 0.;
    
    // combined power law at low energies
    if( x < p[6] )
    {
        f = p[0] * TMath::Power( 10, ( x - p[3] ) * p[1] ) / ( 1 + TMath::Power( 10, ( x - p[3] ) * ( p[1] - p[2] ) ) );
    }
    else
    {
        // polinominal at high energies
        
        // match derivatives at x = p[6]
        double d  = ( 1 + TMath::Power( 10, ( p[6] - p[3] ) * ( p[1] - p[2] ) ) );
        d *= d;
        if( d != 0. )
        {
            d = 1. / d;
        }
        else
        {
            d = 1.;
        }
        
        d = d * p[0] * log( 10. ) * TMath::Power( 10., ( p[6] - p[3] ) * p[1] );
        d *= ( p[1] + p[2] * TMath::Power( 10., ( p[6] - p[3] ) * ( p[1] - p[2] ) ) );
        
        double f2 = 2.*p[4] * p[6] + 3.*p[5] * p[6] * p[6];
        if( f2 != 0. )
        {
            d -= f2;
        }
        else
        {
            d = 1.;
        }
        
        // match function values at x = p[6]
        double c = p[0] * TMath::Power( 10, ( p[6] - p[3] ) * p[1] ) / ( 1 + TMath::Power( 10, ( p[6] - p[3] ) * ( p[1] - p[2] ) ) );
        f2 = p[6] * d + p[6] * p[6] * p[4] + p[6] * p[6] * p[6] * p[5];
        if( f2 != 0. )
        {
            c -= f2;
        }
        else
        {
            c = 1.;
        }
        
        // calculate polynominal
        f = c + x * d + x * x * p[4] + x * x * x * p[5];
    }
    if( f > 0. )
    {
        return log10( f );
    }
    else
    {
        return 0.;
    }
    
    // should never arrive here
    return f;
}


double fun_effLin( double* e, double* p )
{
    double f = fun_eff( e, p );
    
    return TMath::Power( 10., f );
}


TGraphAsymmErrors* VFitEffectiveArea::getLogGraph()
{
    if( !gEff )
    {
        return 0;
    }
    
    if( gEffLog )
    {
        gEffLog->Set( gEff->GetN() );
    }
    else
    {
        gEffLog = new TGraphAsymmErrors( gEff->GetN() );
    }
    
    double x, y;
    for( int i = 0; i < gEff->GetN(); i++ )
    {
        gEff->GetPoint( i, x, y );
        if( y > 0. )
        {
            gEffLog->SetPoint( i, x, log10( y ) );
        }
        else
        {
            gEffLog->SetPoint( i, x, log10( 1. ) );
        }
        x = gEff->GetErrorYhigh( i );
        if( y > 0. )
        {
            gEffLog->SetPointEYhigh( i, x / y );
        }
        else
        {
            gEffLog->SetPointEYhigh( i, 0. );
        }
        x = gEff->GetErrorYlow( i );
        if( y > 0. )
        {
            gEffLog->SetPointEYlow( i, x / y );
        }
        else
        {
            gEffLog->SetPointEYlow( i, 0. );
        }
    }
    
    return gEffLog;
}


void VFitEffectiveArea::write( string ofile )
{
    TFile f( ofile.c_str(), "RECREATE" );
    if( !f.IsZombie() && fTree )
    {
        fTree->Write();
    }
    f.Close();
}


void VFitEffectiveArea::smoothGraph( double iOut, int iIter )
{
    if( !gEff )
    {
        return;
    }
    
    double x1, x2, x3;
    double y1, y2, y3;
    
    for( int t = 0; t < iIter; t++ )
    {
        for( int i = 1; i < gEff->GetN() - 1; i++ )
        {
            gEff->GetPoint( i - 1, x1, y1 );
            gEff->GetPoint( i, x2, y2 );
            gEff->GetPoint( i + 1, x3, y3 );
            
            if( y1 <= 0. || y2 <= 0. || y3 == 0. )
            {
                continue;
            }
            
            if( y2 / y1 > ( 1. + iOut ) && y2 / y3 > ( 1. + iOut ) )
            {
                y2 = 0.5 * ( y1 + y3 );
                cout << "\t" << y2 << endl;
            }
            if( y2 / y1 < ( 1. - iOut ) && y2 / y3 < ( 1. - iOut ) )
            {
                y2 = 0.5 * ( y1 + y3 );
            }
            
            gEff->SetPoint( i, x2, y2 );
        }
    }
    getLogGraph();
}
