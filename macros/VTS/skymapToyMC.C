/*   \file skymapToyMC
     \brief small toy MC to understand dependency of
            significance distributions on bin width,
            integration radius, number of events, etc.

     many hard coded values, please scan through the code
     before using it

     several sections directly copied from VStereoMaps

     e.g. example for events: 80 / bin or 450 / bin for moderate/soft cuts

     Question:

     - how does the significance distribution change as function of integration radius?
     - how large are the fluctuations?
     - how non-Gaussian are the distributions? (introduce check in fit Chi2/N)
     - how big is the impact from small disturbances from e.g. stars? (introduce stars)

     Missing:

     - ring-background model


*/

#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TPad.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TText.h"

#include <iostream>
#include <string>
#include <vector>

#include "../../inc/VStatistics.h"

using namespace std;


/*

    data structure for reflected regions

*/
struct sRE_REGIONS
{
    int noff;                                     //!< number of source regions
    vector< double > xoff;                        //!< x-position of center of off source region
    vector< double > yoff;                        //!< y-position of center of off source region
    vector< double > roff;                        //!< radius of off source region
};


class VToyMap
{
    private:
    
        double f_RE_AreaNorm;
        vector< vector< sRE_REGIONS > > fRE_off;
        
        void fillOn( double x, double y );
        void fillReflectedRegionOff( double x, double y );
        void fillExcessAndSignificanceMaps( bool iPrint );
        
    public:
        TList* hMapList;
        TH2D*  hMapOn;
        TH2D*  hMapOffRE;
        TH2D*  hMapAlpha;
        TH2D*  hMapExcessRE;
        TH2D*  hMapSignificanceRE;
        TH1D*  hSignificance1DRE;
        vector< TH1D* > hSignificance1DRE_dist;
        // Gaussian fit
        double fGausFit_Chi2;
        double fGausFit_Mean;
        double fGausFit_MeanError;
        double fGausFit_Sigma;
        double fGausFit_SigmaError;
        // Gaussian N(0,1) Chi2 test
        double fGausN01_Chi2;
        double fGausN01_Prob;
        // Extreme values
        double fSignificanceMax;
        double fSignificanceMin;
        double fSignificanceMax_dist;
        
        bool fCorrelatedMap;
        double fTheta2;
        double fMaxDistanceCut;
        
        int  fCentreBin_x;
        int  fCentreBin_y;
        
        VToyMap();
        ~VToyMap() {}
        void initializeHistograms( string iName, double iSkyMapSize_deg = 4, double iMapBinSize_deg = 0.01 );
        void initializeReflectedRegionModel( double iTheta2, unsigned int fRE_nMaxoffsource = 6, bool iCorrelatedMap = true );
        int  fill( double x, double y );
        void plot( int i_xbin = -1, int i_ybin = -1, int iPrint = 0, TH2D* hMapAcceptance = 0 );
        void reset();
        void setMaxDistanceCut( double iMaxDistanceCut_deg = 2. )
        {
            fMaxDistanceCut = iMaxDistanceCut_deg;
        }
};

VToyMap::VToyMap()
{
    hMapOn = 0;
    hMapAlpha = 0;
    hMapOffRE = 0;
    hMapExcessRE = 0;
    hMapSignificanceRE = 0;
    hSignificance1DRE = 0;
    f_RE_AreaNorm = 1.;
    fSignificanceMax = -99.;
    fSignificanceMin = 999.;
    fSignificanceMax_dist = -99.;
    
    fCorrelatedMap = true;
    fTheta2 = -1.;
    
    fCentreBin_x = 0;
    fCentreBin_y = 0;
    
    setMaxDistanceCut();
    
    // 1D significances Fit results
    fGausFit_Chi2 = 0.;
    fGausFit_Mean = 0.;
    fGausFit_MeanError = 0.;
    fGausFit_Sigma = 0.;
    fGausFit_SigmaError = 0.;
    fGausN01_Prob = 0.;
    fGausN01_Chi2 = 0.;
}

/*

   plot sky maps

*/
void VToyMap::plot( int i_xbin, int i_ybin, int iPrint, TH2D* hMapAcceptance )
{
    char hname[200];
    sprintf( hname, "cC%d_%d", ( int )( fTheta2 * 1000 ), ( int )fCorrelatedMap );
    char htitle[200];
    sprintf( htitle, "sky maps for theta2 = %.3f (correlated: %d)", fTheta2, ( int )fCorrelatedMap );
    
    TCanvas* c = new TCanvas( hname, htitle, 10, 10, 1400, 680 );
    c->Divide( 4, 2 );
    
    TPad* g = ( TPad* )c->cd( 1 );
    g->SetRightMargin( 0.15 );
    if( hMapOn )
    {
        hMapOn->Draw( "colz" );
        sprintf( hname, "On events (#Theta^{2} = %.3f deg^{2})", fTheta2 );
        TLatex* iTT = new TLatex( -1.6, 1.6, hname );
        iTT->Draw();
    }
    g = ( TPad* )c->cd( 2 );
    g->SetRightMargin( 0.15 );
    if( hMapOffRE )
    {
        hMapOffRE->Draw( "colz" );
        TText* iTT = new TText( -1.6, 1.6, "Off events" );
        iTT->Draw();
    }
    g = ( TPad* )c->cd( 3 );
    g->SetRightMargin( 0.15 );
    if( hMapAcceptance )
    {
        hMapAcceptance->Draw( "colz" );
        TText* iTT = new TText( -1.6, 1.6, "acceptance" );
        iTT->Draw();
    }
    g = ( TPad* )c->cd( 4 );
    g->SetRightMargin( 0.15 );
    if( hMapAlpha )
    {
        hMapAlpha->Draw( "colz" );
        TText* iTT = new TText( -1.6, 1.6, "# of background regions" );
        iTT->Draw();
    }
    g = ( TPad* )c->cd( 5 );
    g->SetRightMargin( 0.15 );
    if( hMapExcessRE )
    {
        hMapExcessRE->Draw( "colz" );
        if( i_xbin < 0 && i_ybin < 0 )
        {
            i_xbin = hMapExcessRE->GetXaxis()->FindBin( -0.5 );
            i_ybin = hMapExcessRE->GetYaxis()->FindBin( -0.5 );
        }
        // plot reflected regions
        if( i_xbin > 0 && i_ybin > 0
                && i_xbin < ( int )fRE_off.size()
                && i_ybin < ( int )fRE_off[i_xbin].size() )
        {
            cout << fRE_off[i_xbin][i_ybin].noff << endl;
            // source region
            TEllipse* iS = new TEllipse( hMapExcessRE->GetXaxis()->GetBinCenter( i_xbin ),
                                         hMapExcessRE->GetYaxis()->GetBinCenter( i_ybin ),
                                         sqrt( fTheta2 ) );
            iS->SetFillStyle( 0 );
            iS->SetLineColor( 2 );
            iS->Draw();
            // off regions
            for( int i = 0; i < fRE_off[i_xbin][i_ybin].noff; i++ )
            {
                TEllipse* iO = new TEllipse( fRE_off[i_xbin][i_ybin].xoff[i],
                                             fRE_off[i_xbin][i_ybin].yoff[i],
                                             fRE_off[i_xbin][i_ybin].roff[i] );
                iO->SetFillStyle( 0 );
                iO->SetLineColor( 5 );
                iO->Draw();
            }
        }
        TText* iTT = new TText( -1.6, 1.6, "excess events" );
        iTT->Draw();
    }
    g = ( TPad* )c->cd( 6 );
    g->SetRightMargin( 0.15 );
    if( hMapSignificanceRE )
    {
        hMapSignificanceRE->SetMaximum( 4. );
        hMapSignificanceRE->SetMinimum( -4. );
        hMapSignificanceRE->Draw( "colz" );
        TText* iTT = new TText( -1.6, 1.6, "significance" );
        iTT->Draw();
    }
    
    g = ( TPad* )c->cd( 7 );
    g->SetLogy( 1 );
    g->SetGridx( 0 );
    g->SetGridy( 0 );
    gStyle->SetOptStat( "eMR" );
    gStyle->SetOptFit( 111 );
    hSignificance1DRE->SetLineWidth( 2 );
    hSignificance1DRE->SetFillColor( 41 );
    hSignificance1DRE->SetFillStyle( 1001 );
    hSignificance1DRE->Draw( "histoe" );
    hSignificance1DRE->Draw( "e same" );
    hSignificance1DRE->Fit( "gaus" );
    /*     for( unsigned int i = 0; i < hSignificance1DRE_dist.size(); i++ )
         {
               hSignificance1DRE_dist[i]->Draw( "same" );
         } */
    // draw 0/1 Gaus
    TF1* f = hSignificance1DRE->GetFunction( "gaus" );
    if( f )
    {
        TF1* i_GG = new TF1( "g", "gaus", -5., 5. );
        i_GG->SetLineColor( 3 );
        i_GG->SetLineWidth( ( Width_t )0.5 );
        i_GG->SetParameters( f->GetParameter( 0 ), 0., 1. );
        i_GG->Draw( "same" );
        
        // draw ratio plot
        
        g = ( TPad* )c->cd( 8 );
        g->SetGridx( 0 );
        g->SetGridy( 0 );
        TH1D* h = new TH1D( "h", "", hSignificance1DRE->GetNbinsX(),
                            hSignificance1DRE->GetXaxis()->GetXmin(),
                            hSignificance1DRE->GetXaxis()->GetXmax() );
        h->SetStats( 0 );
        h->SetMinimum( 0.5 );
        h->SetMaximum( 2. );
        h->SetYTitle( "ratio to N(0,1)" );
        h->SetXTitle( "significance" );
        for( int j = 1; j <= hSignificance1DRE->GetNbinsX(); j++ )
        {
            if( hSignificance1DRE->GetBinContent( j ) < 1 )
            {
                continue;
            }
            double x = hSignificance1DRE->GetBinCenter( j );
            double y = i_GG->Eval( x );
            if( y > 0. )
            {
                h->SetBinContent( j, hSignificance1DRE->GetBinContent( j ) / y );
                h->SetBinError( j, hSignificance1DRE->GetBinError( j ) / y );
            }
        }
        h->Draw();
        TLine* iLR = new TLine( h->GetXaxis()->GetXmin(), 1., h->GetXaxis()->GetXmax(), 1. );
        iLR->Draw();
        sprintf( hname, "Chi2/NDF for N(0,1): %.1f", fGausN01_Chi2 );
        TText* iTT = new TText( -3., 1.8, hname );
        iTT->Draw();
    }
    c->Update();
    if( iPrint )
    {
        sprintf( hname, "cC%d_%d_%d.png", ( int )( fTheta2 * 1000 ), ( int )fCorrelatedMap, iPrint );
        c->Print( hname );
    }
}

/*

    use Li & Ma to calculate significances

*/

void VToyMap::fillExcessAndSignificanceMaps( bool iPrint )
{
    float non = 0.;
    float noff = 0.;
    float alpha = 1.;
    float ndiff = 0.;
    float sig = 0.;
    
    float ndiff_min = 1.e4;
    float nsig_min = 1.e4;
    
    fSignificanceMax = -99.;
    fSignificanceMin = 999.;
    fSignificanceMax_dist = -99.;
    
    double i_dist = 0.;
    
    for( int i = 1; i <= hMapOn->GetNbinsX(); i++ )
    {
        for( int j = 1; j <= hMapOn->GetNbinsX(); j++ )
        {
            hMapExcessRE->SetBinContent( i, j, -9999. );
            hMapSignificanceRE->SetBinContent( i, j, -9999. );
            // reflected region model
            non  = hMapOn->GetBinContent( i, j );
            noff = hMapOffRE->GetBinContent( i, j );
            alpha = ( float )fRE_off[i][j].noff * f_RE_AreaNorm;
            if( alpha > 0. && ( non > 0. || noff > 0. ) )
            {
                ndiff = non - 1. / alpha * noff;
                sig = VStatistics::calcSignificance( non, noff, 1. / alpha );
                hMapExcessRE->SetBinContent( i, j, ndiff );
                if( ndiff < ndiff_min )
                {
                    ndiff_min = ndiff;
                }
                if( sig < nsig_min )
                {
                    nsig_min = sig;
                }
                hMapSignificanceRE->SetBinContent( i, j, sig );
                hSignificance1DRE->Fill( sig );
                i_dist = sqrt( hMapOn->GetXaxis()->GetBinCenter( i ) * hMapOn->GetXaxis()->GetBinCenter( i )
                               + hMapOn->GetYaxis()->GetBinCenter( j ) * hMapOn->GetYaxis()->GetBinCenter( j ) );
                if( sig > fSignificanceMax )
                {
                    fSignificanceMax = sig;
                    fSignificanceMax_dist = i_dist;
                }
                if( sig < fSignificanceMin )
                {
                    fSignificanceMin = sig;
                }
                // fill significance distributions in distance rings
                for( unsigned int d = 0; d < hSignificance1DRE_dist.size(); d++ )
                {
                    if( i_dist > d * 0.5 && i_dist < ( d + 1 ) * 0.5 )
                    {
                        hSignificance1DRE_dist[d]->Fill( sig );
                    }
                }
            }
        }
    }
    hMapExcessRE->SetMinimum( ndiff_min );
    hMapSignificanceRE->SetMinimum( nsig_min );
    
    // fit with a Gaussian function
    TF1* iG = new TF1( "fG", "gaus", -6., 6. );
    hSignificance1DRE->Fit( "fG", "QNS" );
    fGausFit_Chi2 = iG->GetChisquare();
    if( iG->GetNDF() > 0. )
    {
        fGausFit_Chi2 /= iG->GetNDF();
    }
    fGausFit_Mean = iG->GetParameter( 1 );
    fGausFit_MeanError = iG->GetParError( 1 );
    fGausFit_Sigma = iG->GetParameter( 2 );
    fGausFit_SigmaError = iG->GetParError( 2 );
    
    // Chi2 test compared to a N(0,1) Gaussian
    TF1 i_GG( "Fg01", "gaus", -5., 5. );
    i_GG.SetParameters( iG->GetParameter( 0 ), 0., 1. );
    fGausN01_Chi2 = 0.;
    for( int i = 1; i < hSignificance1DRE->GetNbinsX(); i++ )
    {
        double iF = ( hSignificance1DRE->GetBinContent( i )
                      - i_GG.Eval( hSignificance1DRE->GetXaxis()->GetBinCenter( i ) ) );
        iF *= iF;
        iF /= i_GG.Eval( hSignificance1DRE->GetXaxis()->GetBinCenter( i ) );
        
        fGausN01_Chi2 += iF;
    }
    if( hSignificance1DRE->GetNbinsX() > 1 )
    {
        fGausN01_Prob = TMath::Prob( fGausN01_Chi2, hSignificance1DRE->GetNbinsX() - 1 );
        fGausN01_Chi2 /= ( hSignificance1DRE->GetNbinsX() - 1. );
    }
    
    if( iPrint )
    {
        cout << "Theta2: " << fTheta2;
        cout << " Mean: " << fGausFit_Mean << " +- " << fGausFit_MeanError;
        cout << " Sigma: " << fGausFit_Sigma << " +- " << fGausFit_SigmaError;
        cout << " (Fit Chi2: " << fGausFit_Chi2 << ")";
        cout << endl;
        cout << "\t Chi2 test for N(0,1): Chi2/NDF = " << fGausN01_Chi2;
        cout << ", probability " << fGausN01_Prob << endl;
        cout << "\t Maximum significance (nbins=" << hSignificance1DRE->GetNbinsX() << "):";
        cout << fSignificanceMax << " sigma (at " << fSignificanceMax << " deg distance)" << endl;
        cout << "\t Minimum significance (nbins=" << hSignificance1DRE->GetNbinsX() << "):";
        cout << fSignificanceMin << " sigma" << endl;
    }
}

void VToyMap::initializeHistograms( string iName, double iSkyMapSize_deg, double iMapBinSize_deg )
{
    if( iMapBinSize_deg == 0. )
    {
        iMapBinSize_deg = 0.01;
    }
    hMapList = new TList();
    
    char hname[500];
    
    int nbins = ( int )( iSkyMapSize_deg / iMapBinSize_deg );
    
    sprintf( hname, "hMapOn_%s", iName.c_str() );
    hMapOn = new TH2D( hname, "", nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg,
                       nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg );
    hMapOn->SetZTitle( "# of on events" );
    hMapList->Add( hMapOn );
    fCentreBin_x = hMapOn->GetXaxis()->FindBin( 0. );
    fCentreBin_y = hMapOn->GetXaxis()->FindBin( 0. );
    
    sprintf( hname, "hMapAlpha_%s", iName.c_str() );
    hMapAlpha = new TH2D( hname, "", nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg,
                          nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg );
    hMapAlpha->SetZTitle( "1/(normalization alpha)" );
    hMapAlpha->SetStats( 0 );
    
    sprintf( hname, "hMapOffRE_%s", iName.c_str() );
    hMapOffRE = new TH2D( hname, "", nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg,
                          nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg );
    hMapOffRE->SetZTitle( "# of off events" );
    hMapList->Add( hMapOffRE );
    
    sprintf( hname, "hMapExcessRE_%s", iName.c_str() );
    hMapExcessRE = new TH2D( hname, "", nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg,
                             nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg );
    hMapExcessRE->SetZTitle( "# of excess events" );
    hMapList->Add( hMapExcessRE );
    
    sprintf( hname, "hMapSignificanceRE_%s", iName.c_str() );
    hMapSignificanceRE = new TH2D( hname, "", nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg,
                                   nbins, -0.5 * iSkyMapSize_deg, 0.5 * iSkyMapSize_deg );
    hMapSignificanceRE->SetZTitle( "significance" );
    hMapList->Add( hMapSignificanceRE );
    
    sprintf( hname, "hSignificance1DRE_%s", iName.c_str() );
    hSignificance1DRE = new TH1D( hname, "", 60., -5.9, 6.1 );
    hSignificance1DRE->SetXTitle( "significance" );
    hSignificance1DRE->Sumw2( true );
    
    for( unsigned int i = 0; i < 6; i++ )
    {
        sprintf( hname, "hSignificance1DRE_%s_%d", iName.c_str(), i );
        hSignificance1DRE_dist.push_back( new TH1D( hname, "", 60., -6., 6. ) );
        hSignificance1DRE_dist.back()->SetXTitle( "significance" );
        hSignificance1DRE_dist.back()->Sumw2( true );
        hSignificance1DRE_dist.back()->SetLineColor( i + 2 );
        hSignificance1DRE_dist.back()->SetMarkerColor( i + 2 );
    }
    
    
    TIter next( hMapList );
    TH2D* h = 0;
    while( ( h = ( TH2D* )next() ) )
    {
        if( h )
        {
            h->SetStats( 0 );
        }
    }
}

int VToyMap::fill( double x, double y )
{
    /////////////////////////////////////////////////////////////
    // fill on map
    fillOn( x, y );
    /////////////////////////////////////////////////////////////
    // fill reflected region (off)
    fillReflectedRegionOff( x, y );
    
    return hMapOn->GetBinContent( fCentreBin_x, fCentreBin_y );
}

void VToyMap::fillOn( double x, double y )
{
    if( x * x + y * y > fMaxDistanceCut * fMaxDistanceCut )
    {
        return;
    }
    
    /////////////////////
    // uncorrelated on map
    if( !fCorrelatedMap )
    {
        hMapOn->Fill( x, y );
        return;
    }
    /////////////////////
    // correlated on maps
    
    // ** copied from VStereoMaps **
    
    // Constructs a 2D skymap of reconstructed source location on the camera plane
    // smoothed by a box (i.e accept all events within a given radius)
    // thetaCutMax is the radius of the smoothing box
    // (and at the same time sqrt(theta2Max)
    
    // bin numbers for current event
    int i_x = hMapOn->GetXaxis()->FindBin( x );
    int i_y = hMapOn->GetYaxis()->FindBin( y );
    
    // get upper and lower limit for loop (add 2 bins to be on the save side)
    int fn_r0X = int( sqrt( fTheta2 ) / hMapOn->GetXaxis()->GetBinWidth( 2 ) ) + 2;
    int fn_r0Y = int( sqrt( fTheta2 ) / hMapOn->GetYaxis()->GetBinWidth( 2 ) ) + 2;
    
    int ix_start = i_x - fn_r0X;
    if( ix_start < 0 )
    {
        ix_start = 0;
    }
    int ix_stopp = i_x + fn_r0X;
    if( ix_stopp > hMapOn->GetNbinsX() )
    {
        ix_stopp = hMapOn->GetNbinsX();
    }
    
    int iy_start = i_y - fn_r0Y;
    if( iy_start < 0 )
    {
        iy_start = 0;
    }
    int iy_stopp = i_y + fn_r0Y;
    if( iy_stopp > hMapOn->GetNbinsY() )
    {
        iy_stopp = hMapOn->GetNbinsY();
    }
    
    double i_xbin = 0.;
    double i_ybin = 0.;
    double i_r = 0.;
    
    for( int i = ix_start; i < ix_stopp; i++ )
    {
        for( int j = iy_start; j < iy_stopp; j++ )
        {
            i_xbin = hMapOn->GetXaxis()->GetBinCenter( i + 1 );
            i_ybin = hMapOn->GetYaxis()->GetBinCenter( j + 1 );
            // test if this position is inside maximum accepted distance from camera center
            if( sqrt( i_xbin * i_xbin + i_ybin * i_ybin ) > fMaxDistanceCut )
            {
                continue;
            }
            
            // theta2 cut
            i_r = ( ( x - i_xbin ) * ( x - i_xbin )
                    + ( y - i_ybin ) * ( y - i_ybin ) );
            if( i_r <= fTheta2 )
            {
                hMapOn->Fill( i_xbin, i_ybin );
            }
        }
    }
}

/*

    fill background maps for reflected region model

    copied from VStereoMaps

*/
void VToyMap::fillReflectedRegionOff( double x, double y )
{

    ///////////////////////////
    // filling
    
    // first check if (x,y) is inside the fiducal area in the camera
    // (fiducal area is defined as distance to center + ringradius < cameraradius)
    //
    double i_evDist = sqrt( x * x + y * y );
    if( i_evDist > fMaxDistanceCut )
    {
        return;
    }
    
    // now loop over the whole map to check if this event is in one of the off regions
    double i_cx = 0.;
    double i_cy = 0.;
    double i_binDist = 0.;
    unsigned int i_nr = 0;
    
    for( int i = 1; i <= hMapOffRE->GetNbinsX(); i++ )
    {
        i_cx =  hMapOffRE->GetXaxis()->GetBinCenter( i );
        
        for( int j = 1; j <= hMapOffRE->GetNbinsY(); j++ )
        {
            i_cy =  hMapOffRE->GetYaxis()->GetBinCenter( j );
            
            // check if event is in the same ring as this bin (all off regions are in a ring around the camera center)
            i_binDist = sqrt( i_cx * i_cx + i_cy * i_cy );
            if( i_evDist > i_binDist + sqrt( fTheta2 ) )
            {
                continue;
            }
            if( i_evDist < i_binDist - sqrt( fTheta2 ) )
            {
                continue;
            }
            
            for( unsigned int p = 0; p < fRE_off[i][j].xoff.size(); p++ )
            {
                // apply theta2 cut in background region
                double theta2 = ( x - fRE_off[i][j].xoff[p] ) * ( x - fRE_off[i][j].xoff[p] )
                                + ( y - fRE_off[i][j].yoff[p] ) * ( y - fRE_off[i][j].yoff[p] );
                                
                if( theta2 < fRE_off[i][j].roff[p]*fRE_off[i][j].roff[p] )
                {
                    hMapOffRE->Fill( i_cx, i_cy );
                }
            }
        }
    }
}

/*

    empty all sky maps

*/
void VToyMap::reset()
{
    TIter next( hMapList );
    TH2D* h = 0;
    while( ( h = ( TH2D* )next() ) )
    {
        if( h )
        {
            h->Reset();
        }
    }
    hSignificance1DRE->Reset();
    for( unsigned int i = 0; i < hSignificance1DRE_dist.size(); i++ )
    {
        hSignificance1DRE_dist[i]->Reset();
    }
}

/*

     calculate position of reflected regions

*/
void VToyMap::initializeReflectedRegionModel( double iTheta2, unsigned int fRE_nMaxoffsource, bool iCorrelatedMap )
{
    fTheta2 = iTheta2;
    fCorrelatedMap = iCorrelatedMap;
    
    cout << "initialize REFLECTED REGION MODEL for theta2 " << iTheta2;
    cout << ", " << fRE_nMaxoffsource << " off source regions (max)";
    cout << " (correlated: " << iCorrelatedMap << ")";
    cout << endl;
    
    // reflected region variables
    double x = 0.;
    double y = 0.;
    double r = 0.;                                // size of source region
    int n_r = 0;
    
    // off region parameters
    sRE_REGIONS i_off;
    i_off.noff = 0;
    // empty vector
    vector< double > i_Dempty;
    
    // source extension (equal to radius of off regions)
    double fRE_roffTemp = sqrt( fTheta2 );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // calculate off regions (following Zufelde, 2005, p.28)
    //
    //  for all bins in stereo maps
    //
    //    number of off regions depends on distance of bin to camera center (how many areas fit in?)
    //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // set number of bins to be covered
    int i_nbinsX = hMapOn->GetNbinsX();
    int i_nbinsY = hMapOn->GetNbinsY();
    
    // set up 2D vector of off source parameters ([n_x][n_y])
    fRE_off.clear();
    vector< sRE_REGIONS > i_x_off;
    for( int j = 0; j <= i_nbinsY; j++ )
    {
        i_x_off.push_back( i_off );
        i_x_off.back().noff = 0;
    }
    for( int i = 0; i <= i_nbinsX; i++ )
    {
        fRE_off.push_back( i_x_off );
    }
    
    // distance of bin to camera center
    double ids = 0.;
    
    // vectors with off source positions
    vector< double > r_off;
    vector< double > x_off;
    vector< double > y_off;
    // vectors with off source positions (temporary)
    vector< double > r_offTemp;
    vector< double > x_offTemp;
    vector< double > y_offTemp;
    
    // calculate maximum number of off source regions possible
    // (maximum number of off source regions at the edge of the camera)
    int n_max_RE = ( int )( TMath::Pi() / ( asin( fRE_roffTemp / fMaxDistanceCut ) ) );
    
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // loop over all bins in the stereo maps
    //////////////////////////////////////////////////////////////////////////////////////
    for( int i = 0; i <= i_nbinsX; i++ )
    {
        x = hMapOn->GetXaxis()->GetBinCenter( i );
        if( TMath::Abs( x ) < 1.e-5 )
        {
            x = 0.;
        }
        
        for( int j = 0; j <= i_nbinsY; j++ )
        {
            y = hMapOn->GetYaxis()->GetBinCenter( j );
            if( TMath::Abs( y ) < 1.e-5 )
            {
                y = 0.;
            }
            
            // distance of this bin from camera center
            ids = sqrt( x * x + y * y );
            
            r_off.clear();
            x_off.clear();
            y_off.clear();
            
            // initialize reflected region vector for this bin
            n_r = 0;
            fRE_off[i][j].roff = r_off;
            fRE_off[i][j].xoff = x_off;
            fRE_off[i][j].yoff = y_off;
            fRE_off[i][j].noff = 0;
            
            // bin is inside confidence region (distance of bin + off source radius)
            // and bin is not too close to center of camera
            if( ids < fMaxDistanceCut && ids > fRE_roffTemp && fRE_roffTemp != 0. )
            {
                // angular size of the on region seen from the observation position
                double w = asin( fRE_roffTemp / ids );
                // number of off source regions
                if( w > 0 )
                {
                    n_r = ( int )( TMath::Pi() / w );
                }
                else
                {
                    n_r = 0;
                }
                // set maximum number of off source region possible for this particular distance to camera center
                n_max_RE = n_r;
                
                // test if there are enough off source regions
                if( n_r >= 2 )
                {
                    // phi angle of on position
                    double phi_0 = atan2( y, x );
                    
                    // off positions
                    double phi_i = 0.;
                    double x_t = 0.;
                    double y_t = 0.;
                    
                    // try to fit at least fRunList.fRE_nMinoffsource into the available space, work by trial and error
                    // (rotate source start of source regions in XX degrees steps)
                    // (step size is a trade off between accuracy and speed)
                    for( double t = 0.; t < 360.; t += 1.0 )
                    {
                        // reset all previously filled vector
                        r_offTemp.clear();
                        x_offTemp.clear();
                        y_offTemp.clear();
                        // loop over all possible off source positions
                        for( int p = 0; p < n_r; p++ )
                        {
                            // get off-source positions
                            phi_i = phi_0 + TMath::Pi() + ( 2 * p + 1 - n_r ) * w;
                            phi_i += ( ( double )( t ) ) / TMath::RadToDeg();
                            x_t = ids * cos( phi_i );
                            y_t = ids * sin( phi_i );
                            
                            // check if off source region is not included in
                            // this off position (require to be at least 2.0*times theta2 circle away)
                            if( ( x_t - x ) * ( x_t - x )
                                    + ( y_t - y ) * ( y_t - y )
                                    > 4. * fTheta2 )
                            {
                                // fill an off region
                                r_offTemp.push_back( fRE_roffTemp );
                                x_offTemp.push_back( x_t );
                                y_offTemp.push_back( y_t );
                            }
                        }
                        // test if this configuration has more off source regions than a previous one
                        if( x_offTemp.size() > x_off.size() )
                        {
                            r_off = r_offTemp;
                            x_off = x_offTemp;
                            y_off = y_offTemp;
                        }
                        // test if there are enough off source regions
                        if( x_off.size() >= fRE_nMaxoffsource || ( x_off.size() > 10 && ( int )x_off.size() >= n_max_RE - 5 ) )
                        {
                            break;
                        }
                    }
                    n_r = ( int )x_off.size();
                    
                    // remove those regions which are furthest away
                    // (default)
                    double i_dist = 0.;
                    // TMPTMP
                    // fRE_nMaxoffsource = 2 + gRandom->Integer( 4 );
                    // TMPTMP
                    while( r_off.size() > fRE_nMaxoffsource )
                    {
                        double i_max = 0.;
                        unsigned int i_maxtt = 99999;
                        unsigned int tt = 0;
                        for( tt = 0; tt < r_off.size(); tt++ )
                        {
                            i_dist = sqrt( ( x_off[tt] - x ) * ( x_off[tt] - x ) + ( y_off[tt] - y ) * ( y_off[tt] - y ) );
                            if( i_dist > i_max )
                            {
                                i_max = i_dist;
                                i_maxtt = tt;
                            }
                        }
                        if( i_maxtt != 99999 && i_maxtt < r_off.size() )
                        {
                            r_off.erase( r_off.begin() + i_maxtt );
                            x_off.erase( x_off.begin() + i_maxtt );
                            y_off.erase( y_off.begin() + i_maxtt );
                        }
                    }
                    
                    // number of off source regions
                    n_r = ( int )x_off.size();
                }
            }
            // now fill off source regions
            if( n_r < 2 )
            {
                // set n_r to 0 (otherwise problems filling debug tree)
                n_r = 0;
                fRE_off[i][j].roff = i_Dempty;
                fRE_off[i][j].xoff = i_Dempty;
                fRE_off[i][j].yoff = i_Dempty;
                fRE_off[i][j].noff = 0;
            }
            else
            {
                fRE_off[i][j].roff = r_off;
                fRE_off[i][j].xoff = x_off;
                fRE_off[i][j].yoff = y_off;
                fRE_off[i][j].noff = n_r;
            }
            if( n_r > 0. && hMapAlpha )
            {
                hMapAlpha->SetBinContent( i, j, n_r );
            }
            
        }
    }
    // normalisation for uncorrelated plots
    // (off regions are still of source region, this is not the bin size)
    // (TTTH2)
    // add here energy dependent norm for theta2 cut
    f_RE_AreaNorm = 1.;
    if( !fCorrelatedMap && fTheta2 > 0. )
    {
        f_RE_AreaNorm = hMapOffRE->GetXaxis()->GetBinWidth( 2 )
                        * hMapOffRE->GetYaxis()->GetBinWidth( 2 )
                        / TMath::Pi() / fTheta2;
        if( f_RE_AreaNorm > 0. )
        {
            f_RE_AreaNorm = 1. / f_RE_AreaNorm;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    cout << "\t\t ....reflected regions initialized" << endl;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////

class VToyMapMaker
{
    private:
    
        vector< double > fTheta2;
        
        double fSkyMapSize_deg;
        double fSkyMapBinSize_deg;
        
        VToyMap* fMapUC;
        vector< VToyMap* > fMapCC;
        vector< TGraphErrors* > fSignificanceDistFit;
        
        // results
        vector< vector< double > > fExtreme_sig_max;
        vector< vector< double > > fExtreme_sig_min;
        vector< vector< double > > fGausN01_Chi2;
        
        // radial acceptance
        TF1*  fRadialAcceptance;
        TH2F* fRadialAcceptance2D;
        TH2D* hMapAcceptance;
        double fMaxDistanceCut_deg;
        
        // star
        TF1* fStar;
        double fStar_Magnitude;
        double fStar_Size;
        double fStar_x;
        double fStar_y;
        
        bool fUniformAcceptance;
        
        bool hole_due_to_star( double x, double y );
        
    public:
    
        VToyMapMaker( double iSkyMapSize_deg = 4, double iMapBinSize_deg = 0.01 );
        ~VToyMapMaker() {}
        void initialize( string iAcceptanceFile, bool iUse2DRadialAcceptance = false );
        void normalize2DAcceptance();
        void runMonteCarlo( unsigned int nloop = 10, int iCentralBinEvents = 100, bool iPrint = false );
        void plotExtremeValues( unsigned int ntoyMC = 1000, bool iPrint = false );
        void plotResults( bool iPrint = false );
        void setMaxDistanceCut( double iMaxDistanceCut_deg = 2. );
        void setStar( double iMagnitude = 0., double iSize = 0.1, double ix = -0.5, double iy = -0.5 )
        {
            fStar_Magnitude = iMagnitude;
            fStar_Size = iSize;
            fStar_x = ix;
            fStar_y = iy;
        }
        void setUniformAcceptance( bool iSetUniformAcceptance = false )
        {
            fUniformAcceptance = iSetUniformAcceptance;
        }
        void test( unsigned int nloop = 1, int iCentralBinEvents = 100, string iCut = "soft",
                   bool iUse2DRadialAcceptance = false, bool iPrint = false );
};

VToyMapMaker::VToyMapMaker( double iSkyMapSize_deg, double iMapBinSize_deg )
{
    setUniformAcceptance();
    
    fRadialAcceptance = 0;
    fRadialAcceptance2D = 0;
    hMapAcceptance = 0;
    
    // sky map sizes
    fSkyMapSize_deg = iSkyMapSize_deg;
    fSkyMapBinSize_deg = iMapBinSize_deg;
    if( fSkyMapBinSize_deg < 1.e-5 )
    {
        fSkyMapBinSize_deg = 0.01;
    }
    
    // star
    fStar = 0;
    fStar_Magnitude = 0.;
    fStar_Size = 0.1;
    fStar_x = -0.5;
    fStar_y = -0.5;
    
    fMapUC = new VToyMap();
    fMapUC->initializeHistograms( "UC", fSkyMapSize_deg, fSkyMapBinSize_deg );
    fMapUC->initializeReflectedRegionModel( 0.008, 6, false );
    
    // theta vector
    fTheta2.push_back( 0.004 );
    fTheta2.push_back( 0.008 );
    fTheta2.push_back( 0.050 );
    fTheta2.push_back( 0.090 );
    
    char hname[200];
    vector< double > iTemp;
    for( unsigned int i = 0; i < fTheta2.size(); i++ )
    {
        fMapCC.push_back( new VToyMap() );
        sprintf( hname, "CorrelatedMap_theta2_%d", ( int )( fTheta2[i] * 1000 ) );
        fMapCC.back()->initializeHistograms( hname, fSkyMapSize_deg, fSkyMapBinSize_deg );
        fMapCC.back()->initializeReflectedRegionModel( fTheta2[i], 6, true );
        
        fSignificanceDistFit.push_back( new TGraphErrors( 1 ) );
        fSignificanceDistFit.back()->SetName( "" );
        
        fExtreme_sig_max.push_back( iTemp );
        fExtreme_sig_min.push_back( iTemp );
        fGausN01_Chi2.push_back( iTemp );
    }
    setMaxDistanceCut();
}

/*
     set the radial acceptance function
     and fill the 2D acceptance histogram
     (possibility to add inhomogenities)
*/
void VToyMapMaker::initialize( string iAcceptanceFile, bool iUse2DRadialAcceptance )
{
    //////////////////////////
    // star definition
    fStar = new TF1( "fStar", "gaus", 0., 5. );
    fStar->SetParameter( 0, fStar_Magnitude );
    fStar->SetParameter( 1, 0. );
    fStar->SetParameter( 2, fStar_Size );
    cout << "definition of star: " << endl;
    fStar->Print();
    
    hMapAcceptance = new TH2D( "hMapAcceptance", "", 1000, -0.5 * fSkyMapSize_deg, 0.5 * fSkyMapSize_deg,
                               1000, -0.5 * fSkyMapSize_deg, 0.5 * fSkyMapSize_deg );
    hMapAcceptance->SetZTitle( "relative acceptance" );
    hMapAcceptance->SetStats( 0 );
    
    ////////////////////////////////////////
    // read radial acceptance from disk
    TFile iF( iAcceptanceFile.c_str() );
    if( iF.IsZombie() )
    {
        cout << "error reading acceptance file: " << iAcceptanceFile << endl;
        return;
    }
    // 1D acceptance
    fRadialAcceptance = ( TF1* )iF.Get( "fAccZe_0" );
    if( !fRadialAcceptance )
    {
        cout << "error reading acceptance" << endl;
        return;
    }
    // 2D acceptance
    fRadialAcceptance2D = ( TH2F* )iF.Get( "hXYAccTotDeRot" );
    if( !fRadialAcceptance2D )
    {
        cout << "error reading 2D acceptance" << endl;
        return;
    }
    normalize2DAcceptance();
    cout << "using radial acceptance from " << iAcceptanceFile << endl;
    
    // 2D radial acceptance map
    cout << "filling 2D acceptances" << endl;
    
    if( fRadialAcceptance && fRadialAcceptance2D && hMapAcceptance )
    {
        // inefficiency function
        TF1 i_iFeff( "i_iFeff", "1.-0.15-0.15*cos( x )", -1.*TMath::Pi(), TMath::Pi() );
        double i_dist = 0.;
        double x = 0.;
        double y = 0.;
        double p = 0.;
        double r = 0.;
        double norm = 1.;
        for( int i = 1; i <= hMapAcceptance->GetNbinsX(); i++ )
        {
            x = hMapAcceptance->GetXaxis()->GetBinCenter( i );
            for( int j = 1; j <= hMapAcceptance->GetNbinsY(); j++ )
            {
                y = hMapAcceptance->GetYaxis()->GetBinCenter( j );
                i_dist = sqrt( x * x + y * y );
                if( i_dist > fMaxDistanceCut_deg )
                {
                    hMapAcceptance->SetBinContent( i, j, 0. );
                }
                else
                {
                    // add a wave-like inefficiency
                    //                        p = TMath::ATan2( y, x );
                    //                        hMapAcceptance->SetBinContent( i, j, fRadialAcceptance->Eval( i_dist ) * i_iFeff.Eval( p ) );
                    // star
                    r = sqrt( ( x - fStar_x ) * ( x - fStar_x ) + ( y - fStar_y ) * ( y - fStar_y ) );
                    if( fStar )
                    {
                        norm = 1. - fStar->Eval( r );
                    }
                    // TMPTMPTMP
                    norm = 1.;
                    if( !iUse2DRadialAcceptance )
                    {
                        hMapAcceptance->SetBinContent( i, j, fRadialAcceptance->Eval( i_dist ) * norm );
                    }
                    else if( fRadialAcceptance2D )
                    {
                        hMapAcceptance->SetBinContent( i, j, fRadialAcceptance2D->Interpolate( x, y ) );
                    }
                }
            }
        }
    }
    iF.Close();
    
}

void VToyMapMaker::normalize2DAcceptance()
{
    if( !fRadialAcceptance2D )
    {
        return;
    }
    float iNormalizationRadius = 0.35;
    
    float i_nBins = 0;
    float i_nCounts = 0;
    for( int i = 1; i <= fRadialAcceptance2D->GetNbinsX(); i++ )
    {
        float x = fRadialAcceptance2D->GetXaxis()->GetBinCenter( i );
        for( int j = 1; j <= fRadialAcceptance2D->GetNbinsY(); j++ )
        {
            float y = fRadialAcceptance2D->GetYaxis()->GetBinCenter( j );
            
            if( x * x + y * y < iNormalizationRadius * iNormalizationRadius )
            {
                i_nBins++;
                i_nCounts += fRadialAcceptance2D->GetBinContent( i, j );
            }
        }
    }
    if( i_nCounts > 0. )
    {
        fRadialAcceptance2D->Scale( i_nBins / i_nCounts );
        
        cout << "Normalization: " << i_nBins / i_nCounts << endl;
    }
    cout << "smoothing radial acceptance map " << endl;
    fRadialAcceptance2D->Smooth( 1 );
}


bool VToyMapMaker::hole_due_to_star( double x, double y )
{
    // calculate distance to star
    double r = sqrt( ( x - fStar_x ) * ( x - fStar_x ) + ( y - fStar_y ) * ( y - fStar_y ) );
    
    if( fStar && fStar_Magnitude > 0. )
    {
        double u = 1. - fStar->Eval( r );
        double ir = gRandom->Uniform( 1. );
        if( gRandom->Uniform( 1. ) > u )
        {
            return true;
        }
    }
    
    return false;
}

/*
    soft cuts: 500 events in central bin
    moderate cuts: 80 events in central bin
*/
void VToyMapMaker::runMonteCarlo( unsigned int nloop, int iCentralBinEvents, bool iPrint )
{
    double x = 0.;
    double y = 0.;
    double s = 1.;
    double r = 0.;
    double p = 0.;
    double a = 0.;
    
    // reset extreme value vector
    for( unsigned int i = 0; i < fExtreme_sig_max.size(); i++ )
    {
        fExtreme_sig_max[i].clear();
        fExtreme_sig_min[i].clear();
        fGausN01_Chi2[i].clear();
    }
    
    ///////////////////////////////////////////////
    // main loop: produce nloop sky maps
    //
    for( unsigned int i = 0; i < nloop; i++ )
    {
        // reset all sky maps
        fMapUC->reset();
        for( unsigned int j = 0; j < fMapCC.size(); j++ )
        {
            fMapCC[j]->reset();
        }
        int iN = 0;
        while( iN < iCentralBinEvents )
        {
            // get a new random event
            if( fUniformAcceptance )
            {
                x = gRandom->Uniform( 2. );
                y = gRandom->Uniform( 2. );
                s = gRandom->Uniform();
                if( s > 0.25 && s <= 0.75 )
                {
                    x *= -1.;
                }
                if( s > 0.5 )
                {
                    y *= -1.;
                }
            }
            else
            {
                r = 2. * sqrt( gRandom->Uniform( 0, 1. ) );
                p = 2.*TMath::Pi() * gRandom->Uniform( 1. );
                x = r * sin( p );
                y = r * cos( p );
                // acceptance
                a = gRandom->Uniform( 0, 1. );
                //                          if( a > fRadialAcceptance->Eval( r ) )
                if( hMapAcceptance && a > hMapAcceptance->Interpolate( x, y ) )
                {
                    continue;
                }
            }
            // fill sky maps
            fMapUC->fill( x, y );
            for( unsigned int j = 0; j < fMapCC.size(); j++ )
            {
                iN = fMapCC[j]->fill( x, y );
            }
        }
        
        ///////////////////////////////////////////////////////
        // calculate excess and significances
        ///////////////////////////////////////////////////////
        cout << "Results for loop " << i << ":" << endl;
        fMapUC->fillExcessAndSignificanceMaps( true );
        for( unsigned int j = 0; j < fMapCC.size(); j++ )
        {
            fMapCC[j]->fillExcessAndSignificanceMaps( true );
            fSignificanceDistFit[j]->SetPoint( i, fMapCC[j]->fGausFit_Mean,
                                               fMapCC[j]->fGausFit_Sigma );
            fSignificanceDistFit[j]->SetPointError( i, fMapCC[j]->fGausFit_MeanError,
                                                    fMapCC[j]->fGausFit_SigmaError );
                                                    
            fExtreme_sig_max[j].push_back( fMapCC[j]->fSignificanceMax );
            fExtreme_sig_min[j].push_back( fMapCC[j]->fSignificanceMin );
            fGausN01_Chi2[j].push_back( fMapCC[j]->fGausN01_Chi2 );
        }
        
        // print results
        if( iPrint )
        {
            //                fMapUC->plot( -1., -1., i+1 );
            for( unsigned int j = 0; j < fMapCC.size(); j++ )
            {
                fMapCC[j]->plot( -1., -1., i + 1, hMapAcceptance );
            }
        }
    }
}


/*
    plot results for different loops

    sky maps are plotted for the last loop

*/
void VToyMapMaker::plotResults( bool iPrint )
{
    // (don't print uncorrelated map)
    fMapUC->plot( -1, -1, false, hMapAcceptance );
    for( unsigned int i = 0; i < fMapCC.size(); i++ )
    {
        fMapCC[i]->plot( -1, -1, iPrint, hMapAcceptance );
    }
    
    //////////////////////////////////////
    TCanvas* c = new TCanvas( "c", "fit results", 10, 100, 600 * fMapCC.size(), 425 );
    c->Divide( fMapCC.size(), 1 );
    
    for( unsigned int i = 0; i < fSignificanceDistFit.size(); i++ )
    {
        // first row: width vs mean distributions
        TPad* g = ( TPad* )c->cd( i + 1 );
        g->SetGridx( 0 );
        g->SetGridy( 0 );
        char hname[200];
        sprintf( hname, "h%d", i );
        TH1D* h = new TH1D( "h", "", 10, -0.2, 0.2 );
        h->SetStats( 0 );
        h->SetMaximum( 1.5 );
        h->SetMinimum( 0.5 );
        h->SetXTitle( "Mean from Gaussian fit" );
        h->SetYTitle( "Width from Gaussian fit" );
        h->Draw();
        TLine* iL1 = new TLine( 0., h->GetMinimum(), 0., h->GetMaximum() );
        iL1->SetLineStyle( 3 );
        iL1->Draw();
        TLine* iL2 = new TLine( h->GetXaxis()->GetXmin(), 1., h->GetXaxis()->GetXmax(), 1. );
        iL2->SetLineStyle( 3 );
        iL2->Draw();
        
        sprintf( hname, "#Theta^{2} = %.3f deg^{2}", fMapCC[i]->fTheta2 );
        TLatex* iTT = new TLatex( -0.1, 1.4, hname );
        iTT->Draw();
        
        fSignificanceDistFit[i]->Draw( "p" );
    }
    char hname[200];
    if( iPrint )
    {
        c->Print( "SkyToyMCFitResults.png" );
    }
    
    ////////////////////////
    // chi2 distributions
    TCanvas* c2 = new TCanvas( "c2", "chi2", 700, 10, 675, 675 );
    c2->SetGridx( 0 );
    c2->SetGridy( 0 );
    
    double chi2_max = 0.;
    for( unsigned int i = 0; i < fGausN01_Chi2.size(); i++ )
    {
        for( unsigned int j = 0; j < fGausN01_Chi2[i].size(); j++ )
        {
            if( fGausN01_Chi2[i][j] > chi2_max )
            {
                chi2_max = fGausN01_Chi2[i][j];
            }
        }
    }
    
    
    for( unsigned int i = 0; i < fGausN01_Chi2.size(); i++ )
    {
        sprintf( hname, "hChi2_%d", i );
        TH1D* h = new TH1D( hname, "", 100., -2., 3. );
        h->SetStats( 0 );
        h->SetLineWidth( 2 );
        h->SetLineColor( i + 2 );
        h->SetXTitle( "Chi2/NDF for N(0,1)" );
        for( unsigned int j = 0; j < fGausN01_Chi2[i].size(); j++ )
        {
            if( fGausN01_Chi2[i][j] > 0. )
            {
                h->Fill( log10( fGausN01_Chi2[i][j] ) );
            }
        }
        if( i == 0 )
        {
            h->Draw();
        }
        else
        {
            h->Draw( "same" );
        }
    }
    if( iPrint )
    {
        c2->Print( "SkyToyMCChi2Test.png" );
    }
}

/*

    plot extreme significances and compare
    to toy MC from a N(0,1) distribution

*/
void VToyMapMaker::plotExtremeValues( unsigned int ntoyMC, bool iPrint )
{
    /////////////////////////////////////////////////////////
    // determine binning for toy MC of N(0,1)
    // (assume that max distance cut is constant)
    int nBins = 1;
    if( fMapCC.size() > 0 && fMapCC[0] )
    {
        nBins = ( int )( TMath::Pi() * fMapCC[0]->fMaxDistanceCut
                         * fMapCC[0]-> fMaxDistanceCut
                         / fSkyMapBinSize_deg
                         / fSkyMapBinSize_deg );
    }
    
    cout << "total number of bins assumed in simulations: " << nBins << endl;
    
    //////////////////////
    // Gauss N(0,1)
    TH1D* hGausExtrema = new TH1D( "hGausExtrema", "", 200, 2., 7. );
    TF1 hG( "hG", "gaus", -10., 10. );
    hG.SetParameter( 0, 1. );
    hG.SetParameter( 1, 0. );
    hG.SetParameter( 2, 1. );
    
    for( unsigned int i = 0; i < ntoyMC; i++ )
    {
        double sigmax = -99.;
        double sig = 0.;
        
        for( int j = 0; j < nBins; j++ )
        {
            sig = hG.GetRandom();
            if( sig > sigmax )
            {
                sigmax = sig;
            }
        }
        hGausExtrema->Fill( sigmax );
    }
    
    // cumulative distribution
    TH1D* hGausExtremaCum = new TH1D( "hGausExtremaCum", "", hGausExtrema->GetNbinsX(),
                                      hGausExtrema->GetXaxis()->GetXmin(),
                                      hGausExtrema->GetXaxis()->GetXmax() );
    double i_sigCumMax = 1.;
    for( int b = 2; b < hGausExtrema->GetNbinsX(); b++ )
    {
        hGausExtremaCum->SetBinContent( b,  hGausExtremaCum->GetBinContent( b - 1 )
                                        + hGausExtrema->GetBinContent( b ) );
        i_sigCumMax = hGausExtremaCum->GetBinContent( b );
    }
    if( i_sigCumMax > 0. )
    {
        hGausExtremaCum->Scale( 1. / i_sigCumMax );
    }
    
    ///////////////////////////////////////////////////////////////////////
    // plot everything
    TCanvas* cEx = new TCanvas( "cEx", "extreme values of distribution", 300, 10, 1200, 600 );
    
    cEx->Divide( 2, 1 );
    
    TPad* g = ( TPad* )cEx->cd( 1 );
    g->SetGridx( 0 );
    g->SetGridy( 0 );
    g->SetLogy( 1 );
    hGausExtrema->SetStats( 0 );
    hGausExtrema->SetFillStyle( 1001 );
    hGausExtrema->SetFillColor( kGray );
    hGausExtrema->SetLineWidth( 2 );
    hGausExtrema->SetXTitle( "#sigma_{max}" );
    hGausExtrema->Draw();
    
    g = ( TPad* )cEx->cd( 2 );
    g->SetGridx( 0 );
    g->SetGridy( 0 );
    hGausExtremaCum->SetStats( 0 );
    hGausExtremaCum->SetXTitle( "#sigma_{max}" );
    hGausExtremaCum->SetYTitle( "cumulative distribution" );
    hGausExtremaCum->SetMaximum( 1.05 );
    hGausExtremaCum->SetFillStyle( 1001 );
    hGausExtremaCum->SetFillColor( kGray );
    hGausExtremaCum->SetLineWidth( 2 );
    hGausExtremaCum->Draw();
    TLine* iL1 = new TLine( hGausExtrema->GetXaxis()->GetXmin(), 1.,
                            hGausExtrema->GetXaxis()->GetXmax(), 1. );
    iL1->Draw();
    
    ////////////////////////////////////////////
    // now plot extreme sig from toyMap MC
    
    cout << "Extrema plotting: " << fExtreme_sig_max.size() << endl;
    
    char hname[200];
    for( unsigned int i = 0; i < fExtreme_sig_max.size(); i++ )
    {
        sprintf( hname, "hGausToyExtrema_%d", i );
        TH1D* h = new TH1D( hname, "", hGausExtrema->GetNbinsX(),
                            hGausExtrema->GetXaxis()->GetXmin(),
                            hGausExtrema->GetXaxis()->GetXmax() );
        h->SetStats( 0 );
        h->SetLineWidth( 2 );
        h->SetLineColor( i + 2 );
        for( unsigned int j = 0; j < fExtreme_sig_max[i].size(); j++ )
        {
            h->Fill( fExtreme_sig_max[i][j] );
        }
        cEx->cd( 1 );
        h->Draw( "same" );
        
        // cummulative histograms
        sprintf( hname, "hGausToyExtremaCum_%d", i );
        TH1D* hC = new TH1D( hname, "", hGausExtrema->GetNbinsX(),
                             hGausExtrema->GetXaxis()->GetXmin(),
                             hGausExtrema->GetXaxis()->GetXmax() );
        hC->SetStats( 0 );
        hC->SetLineWidth( 2 );
        hC->SetLineColor( i + 2 );
        for( int b = 2; b < h->GetNbinsX(); b++ )
        {
            hC->SetBinContent( b,  hC->GetBinContent( b - 1 )
                               + h->GetBinContent( b ) );
            i_sigCumMax = hC->GetBinContent( b );
        }
        if( i_sigCumMax > 0. )
        {
            hC->Scale( 1. / i_sigCumMax );
        }
        cEx->cd( 2 );
        hC->Draw( "same" );
    }
    
    if( iPrint )
    {
        cEx->Print( "SkyToyMCExtremaTest.png" );
    }
    
}

void VToyMapMaker::setMaxDistanceCut( double iMaxDistanceCut_deg )
{
    fMaxDistanceCut_deg = iMaxDistanceCut_deg;
    if( fMapUC )
    {
        fMapUC->setMaxDistanceCut( fMaxDistanceCut_deg );
    }
    for( unsigned int i = 0; i < fMapCC.size(); i++ )
    {
        if( fMapCC[i] )
        {
            fMapCC[i]->setMaxDistanceCut( fMaxDistanceCut_deg );
        }
    }
}


/*

   simplify the start of a test

*/
void VToyMapMaker::test( unsigned int nloop, int iCentralBinEvents, string iCut, bool iUse2DRadialAcceptance, bool iPrint )
{
    //    initialize( "/Users/maierg/Experiments/VERITAS/Data/EVNDISP_AnalysisFiles/trunk/VTS/RadialAcceptances/radialAcceptance-v470-auxv01-Cut-NTel2-Moderate-GEO-V6-T1234.root",
    if( iCut == "soft" )
    {
        initialize( "/Users/maierg/Experiments/VERITAS/Data/EVNDISP_AnalysisFiles/trunk/VTS/RadialAcceptances/radAccSoft.root",
                    iUse2DRadialAcceptance );
    }
    else if( iCut == "moderate" )
    {
        initialize( "/Users/maierg/Experiments/VERITAS/Data/EVNDISP_AnalysisFiles/trunk/VTS/RadialAcceptances/radAccModerate2Tel.root",
                    iUse2DRadialAcceptance );
    }
    cout << "running MC " << nloop << " times with " << iCentralBinEvents << " events " << endl;
    runMonteCarlo( nloop, iCentralBinEvents, iPrint );
    plotResults( iPrint );
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/*

    test at what point N(0,1) approximation for
    Li & Ma significances start to fail

    Todo:
    - check what happens for 0/0 in Li & Ma

*/

void LiMa( unsigned int nToyMC = 100 )
{
    double alpha = 1. / 6.;
    double non = 0.;
    double noff = 0.;
    double sig = 0.;
    double sig_max = 0.;
    unsigned int ndist = 125000;
    
    double non_max = 50.;
    int non_maxI = ( int )non_max;
    
    // result histogram
    TProfile* h1DGausWidth = new TProfile( "h1DGausWidth", "", 50, 0., non_max, 0, 5. );
    h1DGausWidth->SetStats( 0 );
    h1DGausWidth->SetXTitle( "mean number of on events" );
    h1DGausWidth->SetYTitle( "width" );
    TH2D* h1DGausWidth2D = new TH2D( "h1DGausWidth2D", "", 50, 0., non_max, 50, 0, 5. );
    h1DGausWidth2D->SetStats( 0 );
    h1DGausWidth2D->SetXTitle( "mean number of on events" );
    h1DGausWidth2D->SetYTitle( "Width2D" );
    
    TProfile* h1DGausChi2 = new TProfile( "h1DGausChi2", "", 50, 0., non_max, 0, 100. );
    h1DGausChi2->SetStats( 0 );
    h1DGausChi2->SetXTitle( "mean number of on events" );
    h1DGausChi2->SetYTitle( "Chi2/NDF" );
    TH2D* h1DGausChi22D = new TH2D( "h1DGausChi22D", "", 50, 0., non_max, 50, 0, 100. );
    h1DGausChi22D->SetStats( 0 );
    h1DGausChi22D->SetXTitle( "mean number of on events" );
    h1DGausChi22D->SetYTitle( "Chi2/NDF" );
    
    // 1D significance distributions
    // (filled in each cycle of the toy MC)
    TH1D* hDist = new TH1D( "hDist", "", 60, -5.90, 6.10 );
    hDist->SetStats( 0 );
    hDist->SetXTitle( "significance" );
    TF1* iG = new TF1( "fG", "gaus", -6., 6. );
    iG->SetParameter( 1, 0. );
    iG->SetParameter( 2, 0. );
    iG->SetLineWidth( 1 );
    iG->SetLineColor( 3 );
    iG->SetLineStyle( 2 );
    
    // extreme significances
    TH1D* hGausSigMax = new TH1D( "hGausSigMax", "", 100, 2., 7. );
    hGausSigMax->SetXTitle( "#sigma_{max}" );
    hGausSigMax->SetFillColor( kGray );
    hGausSigMax->SetFillStyle( 1001 );
    hGausSigMax->SetStats( 0 );
    
    //////////////////////
    // Gauss N(0,1)
    cout << "Calculate extrema for N(0,1) distributions: " << endl;
    TH1D* hGausExtrema = new TH1D( "hGausExtrema", "", 200, 2., 7. );
    hGausExtrema->SetXTitle( "#sigma_{max}" );
    hGausExtrema->SetYTitle( "cumulative distribution" );
    hGausExtrema->SetFillColor( kGray );
    hGausExtrema->SetStats( 0 );
    hGausExtrema->SetFillStyle( 1001 );
    TF1 hG( "hG", "gaus", -10., 10. );
    hG.SetParameter( 0, 1. );
    hG.SetParameter( 1, 0. );
    hG.SetParameter( 2, 1. );
    
    for( unsigned int i = 0; i < 1000; i++ )
    {
        double sigmax = -99.;
        double sig = 0.;
        
        for( unsigned int j = 0; j < ndist; j++ )
        {
            sig = hG.GetRandom();
            if( sig > sigmax )
            {
                sigmax = sig;
            }
        }
        hGausExtrema->Fill( sigmax );
    }
    
    char hname[200];
    vector< TH1D* > hLiMaSigMax;
    for( int i = 1; i < non_maxI; i++ )
    {
        sprintf( hname, "hLiMaSigMax_%d", i );
        hLiMaSigMax.push_back( new TH1D( hname, "", 100, 2., 7. ) );
        hLiMaSigMax.back()->SetStats( 0 );
        hLiMaSigMax.back()->SetLineWidth( 2 );
        hLiMaSigMax.back()->SetLineColor( 2 );
        hLiMaSigMax.back()->SetXTitle( "#sigma_{max}" );
    }
    
    // plot last example of distributions
    TCanvas* c1D = new TCanvas( "c1D", "1D distributions", 10, 10, 900, 900 );
    c1D->Divide( 4, 4 );
    
    // i corresponds to average Non
    for( int i = 1; i < non_maxI; i++ )
    {
        // mean number of on and off events
        double non_mean = ( double )i;
        double noff_mean = 1. / alpha * non_mean;
        
        cout << "Performing toy MC for non = " << non_mean;
        cout << ", noff = " << noff_mean << endl;
        
        // perform toy MC
        for( unsigned int m = 0; m < nToyMC; m++ )
        {
            hDist->Reset();
            sig_max = 0.;
            
            for( unsigned int t = 0; t < ndist; t++ )
            {
                non  = gRandom->PoissonD( non_mean );
                noff = gRandom->PoissonD( noff_mean );
                sig =  VStatistics::calcSignificance( non, noff, alpha );
                hDist->Fill( sig );
                if( sig > sig_max )
                {
                    sig_max = sig;
                }
            }
            
            hDist->Fit( "fG", "QNS" );
            h1DGausWidth->Fill( non_mean, iG->GetParameter( 2 ) );
            h1DGausWidth2D->Fill( non_mean, iG->GetParameter( 2 ) );
            
            if( iG->GetNDF() > 0. )
            {
                h1DGausChi2->Fill( non_mean, iG->GetChisquare() / iG->GetNDF() );
                h1DGausChi22D->Fill( non_mean, iG->GetChisquare() / iG->GetNDF() );
            }
            
            if( i - 1 < hLiMaSigMax.size() && hLiMaSigMax[i - 1] )
            {
                hLiMaSigMax[i - 1]->Fill( sig_max );
            }
            
            if( m == nToyMC - 1 && i <= 16 )
            {
                TPad* p = ( TPad* )c1D->cd( i );
                p->SetLogy( 1 );
                p->SetGridx( 0 );
                p->SetGridy( 0 );
                hDist->DrawCopy();
                iG->DrawCopy( "same" );
                sprintf( hname, "Non = %d", i );
                TText* iT = new TText( -4., 2000., hname );
            }
        }
    }
    
    
    // plot results
    TCanvas* c = new TCanvas( "cC", "Li & Ma test", 10, 10, 1200, 600 );
    c->Divide( 2, 1 );
    c->Draw();
    
    TPad* c1 = ( TPad* )c->cd( 1 );
    c1->SetGridx( 0 );
    c1->SetGridy( 0 );
    h1DGausWidth->SetMinimum( 0.90 );
    h1DGausWidth->SetMaximum( 1.10 );
    //     h1DGausWidth2D->Draw( "colz" );
    h1DGausWidth->SetMarkerStyle( 21 );
    h1DGausWidth->Draw();
    
    c1 = ( TPad* )c->cd( 2 );
    c1->SetGridx( 0 );
    c1->SetGridy( 0 );
    c1->SetLogy( 1 );
    //     h1DGausChi22D->Draw( "colz" );
    h1DGausChi2->SetMarkerStyle( 21 );
    h1DGausChi2->Draw();
    
    // cumulative distribution
    for( unsigned int i = 0; i < hLiMaSigMax.size(); i++ )
    {
        if( !hLiMaSigMax[i] )
        {
            continue;
        }
        
        double iNorm = 1.;
        for( int b = 2; b <= hLiMaSigMax[i]->GetNbinsX(); b++ )
        {
            hLiMaSigMax[i]->SetBinContent( b, hLiMaSigMax[i]->GetBinContent( b )
                                           + hLiMaSigMax[i]->GetBinContent( b - 1 ) );
            iNorm = hLiMaSigMax[i]->GetBinContent( b );
        }
        if( iNorm > 1. )
        {
            hLiMaSigMax[i]->Scale( 1. / iNorm );
        }
    }
    double iNorm = 1.;
    for( int b = 2; b <= hGausExtrema->GetNbinsX(); b++ )
    {
        hGausExtrema->SetBinContent( b, hGausExtrema->GetBinContent( b )
                                     + hGausExtrema->GetBinContent( b - 1 ) );
        iNorm = hGausExtrema->GetBinContent( b );
    }
    if( iNorm > 1. )
    {
        hGausExtrema->Scale( 1. / iNorm );
    }
    
    
    TCanvas* cCumu = new TCanvas( "cCumu", "cumulative distributions", 10, 10, 900, 900 );
    cCumu->Divide( 4, 4 );
    
    for( unsigned int i = 0; i < 16; i++ )
    {
        TPad* cU = ( TPad* )cCumu->cd( i + 1 );
        cU->SetGridx( 0 );
        cU->SetGridy( 0 );
        
        if( i < hLiMaSigMax.size() && hLiMaSigMax[i] )
        {
            hGausExtrema->SetMaximum( 1.07 );
            hGausExtrema->Draw();
            TLine* iL = new TLine( hGausExtrema->GetXaxis()->GetXmin(), 1.,
                                   hGausExtrema->GetXaxis()->GetXmax(), 1. );
            iL->Draw();
            sprintf( hname, "<Non>=%d", i + 1 );
            TText* iTT = new TText( 2.5, 0.8, hname );
            iTT->Draw();
            hLiMaSigMax[i]->Draw( "same" );
        }
    }
}
