/*! \class VPlotArrayReconstruction
    \brief plot core and direction reconstruction in sky or camera view

    #########################################################################################################################################

    This class was developed to visualize CTA/AGIS simulations

    (works currently for vertical showers only)

    Input:  eventdisplay output file from simulations

    Plots:

VPlotArrayReconstruction::groundPlot( int iEventNumber = 0, double xmin = -550., double xmax = 550., double iScale = 2., bool bImageAxes = false );

VPlotArrayReconstruction::cameraPlot( int iEventNumber = 0, double iFOV = 8., double iPixelSize = 0.1, bool blotImageAxes = false );

A negative event number will start a search for the next event with reconstructed parameters

#########################################################################################################################################
#########################################################################################################################################
*/

#include "VPlotArrayReconstruction.h"

VPlotArrayReconstruction::VPlotArrayReconstruction( string iname, string ifile )
{
    cout << "###############################################################" << endl;
    cout << " WORKS FOR VERTICAL SHOWERS ONLY" << endl;
    cout << "###############################################################" << endl;
    
    bZombie = false;
    
    fName = iname;
    
    ftcors = 0;
    fshowerpars = 0;
    
    fEventNumber = 0;
    
    hGround = 0;
    cGround = 0;
    
    hCamera = 0;
    cCamera = 0;
    
    setSizeCut();
    setPlotImageAxis();
    setPlotGroundCoordinates();
    setPlotFOV();
    setPlotPixelSize();
    setPlotTelescopeScale();
    setPlotLocalTrigger();
    
    fFile = new TFile( ifile.c_str() );
    if( fFile->IsZombie() )
    {
        cout << "VPlotArrayReconstruction: error opening file " << ifile << endl;
        bZombie = true;
        return;
    }
    
    TTree* t = ( TTree* )fFile->Get( "tcors" );
    if( t )
    {
        read_corsikaIOreader_file();
    }
    else
    {
        read_eventdisplay_file();
    }
}


void VPlotArrayReconstruction::read_corsikaIOreader_file()
{
    int n;
    double x[1000];
    double y[1000];
    double r[1000];
    ftcors = ( TTree* )fFile->Get( "tcors" );
    if( !ftcors )
    {
        cout << "VPlotArrayReconstruction: cannot find tcors tree" << endl;
        bZombie = true;
        return;
    }
    ftcors->SetBranchAddress( "telNumber", &n );
    ftcors->SetBranchAddress( "telXpos", x );
    ftcors->SetBranchAddress( "telYpos", y );
    ftcors->SetBranchAddress( "telR", r );
    if( ftcors->GetEntries() > 0 )
    {
        ftcors->GetEntry( 0 );
        for( int t = 0; t < n; t++ )
        {
            telPos_x.push_back( -1.*y[t] );
            telPos_y.push_back( x[t] );
            telPos_r.push_back( r[t] );
        }
    }
    cout << "total number of telescopes: " << telPos_x.size() << endl;
}


void VPlotArrayReconstruction::read_eventdisplay_file()
{
    float x = 0.;
    float y = 0.;
    TTree* t = ( TTree* )fFile->Get( "telconfig" );
    if( !t )
    {
        cout << "VPlotArrayReconstruction: cannot find telescope position tree" << endl;
        bZombie = true;
        return;
    }
    t->SetBranchAddress( "TelX", &x );
    t->SetBranchAddress( "TelY", &y );
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );
        telPos_x.push_back( x );
        telPos_y.push_back( y );
        telPos_r.push_back( 6. );
    }
    cout << "total number of telescopes: " << telPos_x.size() << endl;
    
    // showerpars tree
    fshowerpars = ( TTree* )fFile->Get( "showerpars" );
    if( !fshowerpars )
    {
        cout << "VPlotArrayReconstruction: cannot find showerpars tree" << endl;
        bZombie = true;
        return;
    }
    // tpars trees
    char hname[600];
    for( unsigned int i = 0; i < telPos_x.size(); i++ )
    {
        sprintf( hname, "Tel_%d", i + 1 );
        if( fFile->cd( hname ) )
        {
            ftpars.push_back( ( TTree* )gDirectory->Get( "tpars" ) );
            if( !ftpars.back() )
            {
                cout << "VPlotArrayReconstruction: cannot find tpars tree for telescope " << i + 1 << endl;
                cout << "(" << hname << ")" << endl;
                bZombie = true;
            }
        }
    }
}


void VPlotArrayReconstruction::printEvents()
{
    int icounter = 0;
    while( icounter >= 0 )
    {
        icounter = getEvent( -1 * ( icounter + 1 ) );
    }
}


void VPlotArrayReconstruction::printEvent( int iEventNumber )
{
    if( getEvent( iEventNumber ) < 0 )
    {
        return;
    }
    
    if( fsize.size() != telPos_x.size() )
    {
        return;
    }
    if( fsize.size() != telPos_y.size() )
    {
        return;
    }
    if( fsize.size() != fltrig.size() )
    {
        return;
    }
    if( fsize.size() != fwidth.size() )
    {
        return;
    }
    if( fsize.size() != flength.size() )
    {
        return;
    }
    if( fsize.size() != fdist.size() )
    {
        return;
    }
    if( fsize.size() != floss.size() )
    {
        return;
    }
    if( fsize.size() != fntubes.size() )
    {
        return;
    }
    
    cout << "\t\t ltrig \t ntubes \t size \t dist \t loss" << endl;
    for( unsigned int i = 0; i < fsize.size(); i++ )
    {
        cout << "Tel " << i + 1 << "\t" << telPos_x[i] << "\t" << telPos_y[i] << "\t" << fltrig[i] << "\t" << fntubes[i] << "\t" << fsize[i] << "\t" << fdist[i] << "\t" << floss[i] <<  endl;
    }
    
}


int VPlotArrayReconstruction::getEvent( int iEventNumber )
{
    if( bZombie )
    {
        return -1;
    }
    
    fImgSel.clear();
    fsize.clear();
    fltrig.clear();
    fcen_x.clear();
    fcen_y.clear();
    fphi.clear();
    flength.clear();
    fwidth.clear();
    floss.clear();
    fntubes.clear();
    fdist.clear();
    
    if( ftcors )
    {
        return get_corsikaIOreader_event( iEventNumber );
    }
    
    return get_eventdisplay_event( iEventNumber );
}


int VPlotArrayReconstruction::get_corsikaIOreader_event( int iEventNumber )
{
    if( !ftcors )
    {
        return -1;
    }
    
    fphotons.clear();
    
    double e0 = 0.;
    double az = 0.;
    double ze = 0.;
    double xc = 0.;
    double yc = 0.;
    int n = 0;
    double np[1000];
    
    ftcors->SetBranchAddress( "energy", &e0 );
    ftcors->SetBranchAddress( "az", &az );
    ftcors->SetBranchAddress( "ze", &ze );
    ftcors->SetBranchAddress( "telNumber", &n );
    ftcors->SetBranchAddress( "xCore", &xc );
    ftcors->SetBranchAddress( "yCore", &yc );
    ftcors->SetBranchAddress( "NCp", np );
    
    iEventNumber = TMath::Abs( iEventNumber );
    
    if( iEventNumber < ftcors->GetEntries() )
    {
        ftcors->GetEntry( iEventNumber );
        
        fMCEnergy = e0 / 1.e3;
        fMCxcore = -1.*yc;
        fMCycore = xc;
        fMCaz = az;
        fMCze = ze;
        for( int t = 0; t < n; t++ )
        {
            fphotons.push_back( np[t] );
        }
        
        cout << "found event " << iEventNumber << endl;
        cout << "\tMC energy [TeV]: " << fMCEnergy << endl;
        cout << "\tMC ze, az: " << fMCze << " " << fMCaz << endl;
        cout << "\tMC core: " << fMCxcore << " " << fMCycore << " " << sqrt( fMCxcore * fMCxcore + fMCycore * fMCycore ) << endl;
        
        fEventNumber = iEventNumber;
    }
    else
    {
        fEventNumber = -1 * iEventNumber;
    }
    
    return fEventNumber;
}


int VPlotArrayReconstruction::get_eventdisplay_event( int iEventNumber )
{
    if( fshowerpars )
    {
        fImgSel.clear();
        unsigned long int iimgsel[1000];
        fshowerpars->SetBranchAddress( "MCe0", &fMCEnergy );
        fshowerpars->SetBranchAddress( "NImages", fnimages );
        fshowerpars->SetBranchAddress( "Chi2", fchi2 );
        fshowerpars->SetBranchAddress( "Ze", fze );
        fshowerpars->SetBranchAddress( "Az", faz );
        fshowerpars->SetBranchAddress( "Xcore", fxcore );
        fshowerpars->SetBranchAddress( "Ycore", fycore );
        fshowerpars->SetBranchAddress( "Xoff", fxoff );
        fshowerpars->SetBranchAddress( "Yoff", fyoff );
        fshowerpars->SetBranchAddress( "MCaz", &fMCaz );
        fshowerpars->SetBranchAddress( "MCze", &fMCze );
        fshowerpars->SetBranchAddress( "MCxcore", &fMCxcore );
        fshowerpars->SetBranchAddress( "MCycore", &fMCycore );
        fshowerpars->SetBranchAddress( "MCxoff", &fMCxoff );
        fshowerpars->SetBranchAddress( "MCyoff", &fMCyoff );
        fshowerpars->SetBranchAddress( "ImgSel", &iimgsel );
        
        if( iEventNumber < 0 )
        {
            for( int i = -1 * iEventNumber; i < fshowerpars->GetEntries(); i++ )
            {
                fshowerpars->GetEntry( i );
                if( fnimages[0] > 1 )
                {
                    iEventNumber = i;
                    break;
                }
            }
        }
        else if( iEventNumber < fshowerpars->GetEntries() )
        {
            fshowerpars->GetEntry( iEventNumber );
        }
        else
        {
            fEventNumber = -1 * iEventNumber;
            return fEventNumber;
        }
        
        if( iEventNumber < 0 )
        {
            fEventNumber = iEventNumber;
            return fEventNumber;
        }
        bitset<8 * sizeof( unsigned long )> i_sel( iimgsel[0] );
        for( unsigned int t = 0; t < i_sel.size(); t++ )
        {
            fImgSel.push_back( i_sel.test( t ) );
        }
    }
    
    cout << "found event " << iEventNumber << endl;
    cout << "\tMC energy [TeV]: " << fMCEnergy << endl;
    cout << "\tNumber of images: " << fnimages[0] << endl;
    cout << "\tMC ze, az: " << fMCze << " " << fMCaz << endl;
    cout << "\tMC off: " << fMCxoff << " " << fMCyoff << " " << sqrt( fMCxoff * fMCxoff + fMCyoff * fMCyoff ) << endl;
    cout << "\tMC core: " << fMCxcore << " " << fMCycore << " " << sqrt( fMCxcore * fMCxcore + fMCycore * fMCycore ) << endl;
    if( fchi2[0] >= 0 )
    {
        cout << "\tRE ze, az: " << fze[0] << " " << faz[0] << endl;
        cout << "\tRE off: " << fxoff[0] << " " << fyoff[0] << " " << sqrt( fxoff[0]*fxoff[0] + fyoff[0]*fyoff[0] ) << endl;
        cout << "\terror direction: " << sqrt( ( fMCxoff - fxoff[0] ) * ( fMCxoff - fxoff[0] ) + ( fMCyoff - fyoff[0] ) * ( fMCyoff - fyoff[0] ) ) << endl;
        cout << "\tRE core: " << fxcore[0] << " " << fycore[0] << " " << sqrt( fxcore[0]*fxcore[0] + fycore[0]*fycore[0] ) << endl;
        cout << "\terror core: " << sqrt( ( fMCxcore - fxcore[0] ) * ( fMCxcore - fxcore[0] ) + ( fMCycore - fycore[0] ) * ( fMCycore - fycore[0] ) ) << endl;
        
    }
    
    float isize;
    short iltrig;
    float icenx;
    float iceny;
    float iphi;
    float ilength;
    float iwidth;
    float idist;
    float iloss;
    UShort_t intubes;
    
    for( unsigned int i = 0; i < ftpars.size(); i++ )
    {
        if( ftpars[i] && iEventNumber < ftpars[i]->GetEntries() )
        {
            ftpars[i]->SetBranchAddress( "size", &isize );
            ftpars[i]->SetBranchAddress( "cen_x", &icenx );
            ftpars[i]->SetBranchAddress( "cen_y", &iceny );
            ftpars[i]->SetBranchAddress( "phi", &iphi );
            ftpars[i]->SetBranchAddress( "length", &ilength );
            ftpars[i]->SetBranchAddress( "width", &iwidth );
            ftpars[i]->SetBranchAddress( "dist", &idist );
            ftpars[i]->SetBranchAddress( "ntubes", &intubes );
            ftpars[i]->SetBranchAddress( "loss", &iloss );
            ftpars[i]->SetBranchAddress( "LocalTrigger", &iltrig );
            
            ftpars[i]->GetEntry( iEventNumber );
            
            fsize.push_back( isize );
            fltrig.push_back( iltrig == 1 );
            fcen_x.push_back( icenx );
            fcen_y.push_back( iceny );
            fphi.push_back( iphi );
            flength.push_back( ilength );
            fwidth.push_back( iwidth );
            fdist.push_back( idist );
            floss.push_back( iloss );
            fntubes.push_back( intubes );
        }
    }
    for( unsigned int i = 0; i < fntubes.size(); i++ )
    {
        if( fntubes[i] > 0 )
        {
            cout << "\t Telescope " << i << "\t" << fntubes[i] << "\t" << fsize[i] << "\t" << fdist[i] << "\t" << fltrig[i] << endl;
        }
    }
    
    fEventNumber = iEventNumber;
    return fEventNumber;
}


VPlotArrayReconstruction::~VPlotArrayReconstruction()
{

}


TCanvas* VPlotArrayReconstruction::groundPlot( int iEventNumber )
{
    iEventNumber = getEvent( iEventNumber );
    if( iEventNumber < 0 )
    {
        return 0;
    }
    if( bZombie )
    {
        return 0;
    }
    
    cout << "plotting event " << iEventNumber << endl;
    
    char hname[800];
    char htitle[800];
    sprintf( hname, "hGround_%s_%d", fName.c_str(), iEventNumber );
    hGround = new TH2D( hname, "", 10, fPlotxmin, fPlotxmax, 10, fPlotymin, fPlotymax );
    hGround->SetStats( 0 );
    hGround->SetXTitle( "EAST [m]" );
    hGround->SetYTitle( "NORTH [m]" );
    hGround->SetZTitle( "log_{10} size" );
    hGround->GetZaxis()->SetTitleOffset( 1.2 );
    hGround->GetYaxis()->SetTitleOffset( 1.5 );
    // plot everything
    if( !cGround )
    {
        sprintf( hname, "cGround_%s", fName.c_str() );
        sprintf( htitle, "ground array (%s, event %d)", fName.c_str(), fEventNumber );
        cGround = new TCanvas( hname, htitle, 10, 10, 640, 600 );
        cGround->SetGridx( 0 );
        cGround->SetGridy( 0 );
        cGround->SetLeftMargin( 0.13 );
    }
    else
    {
        sprintf( htitle, "ground array (%s, event %d)", fName.c_str(), fEventNumber );
        cGround->SetTitle( htitle );
        cGround->cd();
    }
    
    hGround->Draw( );
    
    // plot telescope positions
    for( unsigned int i = 0; i < telPos_x.size(); i++ )
    {
        TEllipse* iE = new TEllipse( telPos_x[i], telPos_y[i], telPos_r[i], telPos_r[i] );
        iE->Draw();
    }
    
    // get min/max size or nphotons
    double iMin = 1.e20;
    double iMax = 0.;
    for( unsigned int i = 0; i < fsize.size(); i++ )
    {
        if( fsize[i] > 0 && log10( fsize[i] ) < iMin )
        {
            iMin = log10( fsize[i] );
        }
        if( fsize[i] > 0 && log10( fsize[i] ) > iMax )
        {
            iMax = log10( fsize[i] );
        }
    }
    for( unsigned int i = 0; i < fphotons.size(); i++ )
    {
        if( fphotons[i] > 0 && log10( fphotons[i] ) < iMin )
        {
            iMin = log10( fphotons[i] );
        }
        if( fphotons[i] > 0 && log10( fphotons[i] ) > iMax )
        {
            iMax = log10( fphotons[i] );
        }
    }
    setColorAxisDataVector_minmax( iMin * 0.98, iMax * 1.02 );
    
    // fill histogram
    for( unsigned int i = 0; i < telPos_x.size(); i++ )
    {
        if( i < fsize.size() && fsize[i] > 0. )
        {
            // plot telescopes
            TEllipse* iT = new TEllipse( telPos_x[i], telPos_y[i], telPos_r[i]*fPlotTelescopeScale, telPos_r[i]*fPlotTelescopeScale );
            iT->SetFillStyle( 1001 );
            iT->SetFillColor( getColorAxisColor( log10( fsize[i] ) ) );
            iT->SetLineColor( getColorAxisColor( log10( fsize[i] ) ) );
            iT->Draw();
        }
        else if( i < fphotons.size() && fphotons[i] > 0 )
        {
            TEllipse* iT = new TEllipse( telPos_x[i], telPos_y[i], telPos_r[i]*fPlotTelescopeScale, telPos_r[i]*fPlotTelescopeScale );
            iT->SetFillStyle( 1001 );
            iT->SetFillColor( getColorAxisColor( log10( fphotons[i] ) ) );
            iT->SetLineColor( getColorAxisColor( log10( fphotons[i] ) ) );
            iT->Draw();
        }
        if( i < fltrig.size() && fltrig[i] && bPlotLocalTrigger )
        {
            TMarker* iLT = new TMarker( telPos_x[i], telPos_y[i], 2 );
            iLT->Draw();
        }
    }
    // plot shower axis
    TGaxis* g = getColorAxisAxis( fPlotxmin + 0.05 * ( fPlotxmax - fPlotxmin ), fPlotxmin + 0.09 * ( fPlotxmax - fPlotxmin ), fPlotymin + 0.09 * ( fPlotymax - fPlotymin ) / 2., fPlotymin + 0.4 * ( fPlotymax - fPlotymin ), "log_{10} size", 10 );
    if( g && iMax > 0. )
    {
        g->Draw();
    }
    
    // plot core positions
    TMarker* icoreMC = new TMarker( fMCxcore, fMCycore, 5 );
    icoreMC->Draw();
    if( fchi2[0] >= 0 )
    {
        TMarker* icore   = new TMarker( fxcore[0], fycore[0], 24 );
        icore->Draw();
    }
    
    if( fshowerpars && bPlotImageAxes )
    {
        // plot image lines
        for( unsigned int i = 0; i < fsize.size(); i++ )
        {
            if( fsize[i] > fSizeCut )
            {
                double fPlotTelescopeScale = ( fPlotxmax - fPlotxmin ) / 2.;
                // distance between image centre and reconstruction direction
                double iMC_dist = sqrt( ( telPos_x[i] - fMCxcore ) * ( telPos_x[i] - fMCxcore ) + ( telPos_y[i] - fMCycore ) * ( telPos_y[i] - fMCycore ) );
                double iRE_dist = sqrt( ( telPos_x[i] - fxcore[0] ) * ( telPos_x[i] - fxcore[0] ) + ( telPos_y[i] - fycore[0] ) * ( telPos_y[i] - fycore[0] ) );
                if( fnimages[0] > 0 )
                {
                    fPlotTelescopeScale = max( iMC_dist, iRE_dist ) * 1.1;
                }
                else
                {
                    fPlotTelescopeScale = iMC_dist * 1.1;
                }
                
                TLine* iL = new TLine();
                float i_y = -1.*fyoff[0];
                float i_x = -1.*fxoff[0];
                float i_cen_x = fcen_x[i] + i_x;
                float i_cen_y = fcen_y[i] + i_y;
                float fMCSign = -1.;
                
                float i_X1 = telPos_x[i] - fPlotTelescopeScale * cos( atan2( fMCSign * i_cen_y, i_cen_x ) );
                float i_X2 = telPos_x[i] + ( fPlotxmax - fPlotxmin ) / 50. * cos( atan2( fMCSign * i_cen_y, i_cen_x ) );
                float i_Y1 = telPos_y[i] - fPlotTelescopeScale * sin( atan2( fMCSign * i_cen_y, i_cen_x ) );
                float i_Y2 = telPos_y[i] + ( fPlotymax - fPlotymin ) / 50. * sin( atan2( fMCSign * i_cen_y, i_cen_x ) );
                
                iL->SetX1( i_X1 );
                iL->SetX2( i_X2 );
                iL->SetY1( i_Y1 );
                iL->SetY2( i_Y2 );
                
                iL->SetLineWidth( 2 );
                iL->SetLineStyle( 2 );
                iL->SetLineColor( getColorAxisColor( log10( fsize[i] ) ) );
                
                iL->Draw();
            }
        }
    }
    return cGround;
}


void VPlotArrayReconstruction::printTelescopes()
{
    for( unsigned int i = 0; i < telPos_x.size(); i++ )
    {
        cout << "Tel " << i + 1 << "\t" << telPos_x[i] << "\t" << telPos_y[i] << "\t" << telPos_r[i] << endl;
    }
}


void VPlotArrayReconstruction::nextEvent()
{
    groundPlot( -1 * fEventNumber - 1 );
    cameraPlot( fEventNumber );
}


TCanvas* VPlotArrayReconstruction::cameraPlot( int iEventNumber )
{
    iEventNumber = getEvent( iEventNumber );
    if( iEventNumber < 0 )
    {
        return 0;
    }
    if( bZombie )
    {
        return 0;
    }
    
    default_settings();
    
    cout << "plotting event " << iEventNumber << endl;
    
    char hname[800];
    char htitle[800];
    if( !cCamera )
    {
        sprintf( hname, "cCamera_%s", fName.c_str() );
        sprintf( htitle, "camera (%s, event %d)", fName.c_str(), fEventNumber );
        cCamera = new TCanvas( hname, htitle, 710, 10, 640, 600 );
        cCamera->SetGridx( 0 );
        cCamera->SetGridy( 0 );
        cCamera->SetRightMargin( 0.15 );
    }
    else
    {
        sprintf( htitle, "camera (%s, event %d)", fName.c_str(), fEventNumber );
        cCamera->SetTitle( htitle );
        cCamera->cd();
    }
    
    sprintf( hname, "hCamera_%s_%d", fName.c_str(), iEventNumber );
    hCamera = new TH2D( hname, "", 10, -1.*fPlotFOV / 2.*1.1, fPlotFOV / 2.*1.1, 10, -1.*fPlotFOV / 2.*1.1, fPlotFOV / 2.*1.1 );
    hCamera->SetStats( 0 );
    hCamera->SetXTitle( "x [deg]" );
    hCamera->SetYTitle( "y [deg]" );
    hCamera->Draw();
    
    // plot camera size
    TEllipse* iC = new TEllipse( 0., 0., fPlotFOV / 2., fPlotFOV / 2. );
    iC->SetLineWidth( 1 );
    iC->SetLineStyle( 3 );
    iC->SetFillStyle( 1001 );
    iC->SetFillColor( 18 );
    iC->Draw();
    // plot central pixel
    TEllipse* iZ = new TEllipse( 0., 0., fPlotPixelSize / 2., fPlotPixelSize / 2. );
    iZ->SetFillStyle( 1001 );
    iZ->SetFillColor( 1 );
    iZ->Draw();
    // plot reconstructed core position
    TMarker* iM = new TMarker( fxoff[0], -1.*fyoff[0], 5 );
    iM->SetMarkerSize( 2 );
    iM->SetMarkerColor( 2 );
    iM->Draw();
    
    // draw analysis results
    if( fshowerpars )
    {
    
        // get min/max size
        double iMin = 1.e20;
        double iMax = 0.;
        for( unsigned int i = 0; i < fsize.size(); i++ )
        {
            if( fsize[i] > 0 && log10( fsize[i] ) < iMin )
            {
                iMin = log10( fsize[i] );
            }
            if( fsize[i] > 0 && log10( fsize[i] ) > iMax )
            {
                iMax = log10( fsize[i] );
            }
        }
        setColorAxisDataVector_minmax( iMin * 0.98, iMax * 1.02 );
        
        for( unsigned int i = 0; i < fsize.size(); i++ )
        {
            if( i < fsize.size() && fsize[i] )
            {
                // draw hillas ellipse
                TEllipse* iE = new TEllipse( fcen_x[i], fcen_y[i], flength[i], fwidth[i], 0., 360., fphi[i] * 180. / TMath::Pi() );
                iE->SetFillColor( getColorAxisColor( log10( fsize[i] ) ) );
                iE->SetFillStyle( 1001 );
                iE->SetLineColor( getColorAxisColor( log10( fsize[i] ) ) );
                iE->Draw();
                
                // draw length axis of images
                if( bPlotImageAxes )
                {
                    //		if( fImgSel[i] )
                    if( fsize[i] > fSizeCut )
                    {
                        double iScale = fPlotFOV / 2.;
                        // distance between image centre and reconstruction direction
                        double iMC_dist = sqrt( ( fcen_x[i] - fMCxoff ) * ( fcen_x[i] - fMCxoff ) + ( fcen_y[i] - fMCyoff ) * ( fcen_y[i] - fMCyoff ) );
                        double iRE_dist = sqrt( ( fcen_x[i] - fxoff[0] ) * ( fcen_x[i] - fxoff[0] ) + ( fcen_y[i] - fyoff[0] ) * ( fcen_y[i] - fyoff[0] ) );
                        iScale = max( iMC_dist, iRE_dist ) * 1.1;
                        double i_x1 = fcen_x[i] + fPlotFOV / 50. * cos( fphi[i] );
                        double i_x2 = fcen_x[i] - iScale * cos( fphi[i] );
                        double i_y1 = fcen_y[i] + fPlotFOV / 50.* sin( fphi[i] );
                        double i_y2 = fcen_y[i] - iScale * sin( fphi[i] );
                        
                        TLine* iL = new TLine( i_x1, i_y1, i_x2, i_y2 );
                        iL->SetLineColor( getColorAxisColor( log10( fsize[i] ) ) );
                        iL->SetLineStyle( 2 );
                        iL->Draw();
                    }
                }
            }
            // plot shower axis
            TGaxis* g = getColorAxisAxis( 0.85 * fPlotFOV / 2., 0.91 * fPlotFOV / 2., 0.6 * fPlotFOV / 2., 0.98 * fPlotFOV / 2., "log_{10} size", 10 );
            if( g && iMax > 0. )
            {
                g->Draw();
            }
        }
    }
    
    return cCamera;
}


void VPlotArrayReconstruction::setPlotGroundCoordinates( double x1, double x2, double y1, double y2 )
{
    fPlotxmin = x1;
    fPlotxmax = x2;
    fPlotymin = y1;
    fPlotymax = y2;
}


void VPlotArrayReconstruction::printToFile( string iFileName, string iFileType )
{
    if( cGround )
    {
        char hname[2000];
        sprintf( hname, "%s-%d-ground%s", iFileName.c_str(), fEventNumber, iFileType.c_str() );
        cGround->Print( hname );
    }
    if( cCamera )
    {
        char hname[2000];
        sprintf( hname, "%s-%d-camera%s", iFileName.c_str(), fEventNumber, iFileType.c_str() );
        cCamera->Print( hname );
    }
}


void VPlotArrayReconstruction::plotEvent( int iEventNumber )
{
    groundPlot( iEventNumber );
    cameraPlot( iEventNumber );
}
