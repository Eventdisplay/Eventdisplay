/* plot IPR charge distributions for data and sims
 *
 * (uses IPRcontours files produced by NN cleaning step in evndisp)
 *
 * usage:
 * VPlotIPRCharge a;
 * a.loadIPRGraphs( "58405.IPRcontours.root", false );
 * a.plot( 5.e2 );
 */

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TPad.h"

class VPlotIPRCharge
{
    private:
    
        vector< TGraphErrors* > fIPR_MC;
        vector< TGraphErrors* > fIPR_Data;
        
        
    public:
    
        VPlotIPRCharge();
        ~VPlotIPRCharge() {}
        void plot( double iMin = -99. );
        void loadIPRGraphs( string iFile, bool iMC );
        
};

VPlotIPRCharge::VPlotIPRCharge()
{

}

void VPlotIPRCharge::loadIPRGraphs( string iFile, bool iMC )
{
    TFile* iR = new TFile( iFile.c_str() );
    if( iR->IsZombie() )
    {
        return;
    }
    
    char hname[200];
    for( unsigned int i = 0; i < 4; i++ )
    {
        sprintf( hname, "IPRchargeTelType%d", i );
        if( iMC )
        {
            fIPR_MC.push_back( ( TGraphErrors* )iR->Get( hname ) );
        }
        else
        {
            fIPR_Data.push_back( ( TGraphErrors* )iR->Get( hname ) );
        }
    }
}

void VPlotIPRCharge::plot( double iMin )
{

    TCanvas* c = new TCanvas( "cIPR", "IPR charges", 10, 10, 1200, 700 );
    c->Divide( 2, 2 );
    
    char hname[200];
    for( unsigned int i = 0; i < 4; i++ )
    {
        TPad* p = ( TPad* )c->cd( i + 1 );
        p->SetLogy( 1 );
        
        double xmax = 100.;
        
        if( i < fIPR_Data.size() && fIPR_Data[i] )
        {
            if( iMin > 0. )
            {
                fIPR_Data[i]->SetMinimum( iMin );
            }
            fIPR_Data[i]->SetLineWidth( 2 );
            fIPR_Data[i]->Draw( "apl" );
            if( fIPR_Data[i]->GetHistogram() )
            {
                xmax = fIPR_Data[i]->GetHistogram()->GetXaxis()->GetXmax();
            }
        }
        
        if( i < fIPR_MC.size() && fIPR_MC[i] )
        {
            if( iMin > 0. )
            {
                fIPR_MC[i]->SetMinimum( iMin );
            }
            fIPR_MC[i]->SetLineWidth( 2 );
            fIPR_MC[i]->SetLineColor( 2 );
            fIPR_MC[i]->SetMarkerColor( 2 );
            fIPR_MC[i]->Draw( "pl" );
        }
        sprintf( hname, "Tel %d", i + 1 );
        TText* iT = new TText( 0.7 * xmax, 1.e6, hname );
        iT->SetTextSize( iT->GetTextSize() * 1.5 );
        iT->Draw();
    }
    
}

