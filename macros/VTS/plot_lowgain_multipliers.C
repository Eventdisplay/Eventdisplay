#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <vector>
#include <stdlib.h>
#include "TSystem.h"
#include "TF1.h"
#include "TMath.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>


ofstream logfile;

bool is_good_run( int run, int tel )
{
    /* 7083132 T4: completely wrong window.
    # 5769394, 5811314, 5880607 T2: Pulse very late -> not reliable for lpeds
    # 5423738, 5423941, 5508384, 5625859, 5625860 -> 64 sample readout
    # 7116263, 7217071, 7563739 7639596 7679293 T4 flasher to low
    # 7348488 T2,T3 flasher too low
    # 7422021 event loss
    # 7468182 T1, T2 flasher too low
    */
    
    if( ( run == 69475 || run == 69476 ) && ( tel == 1 || tel == 3 || tel == 4 ) )
    {
        return false;    //timing issue (start of trace not recorded)
    }
    if( run == 57384 )
    {
        return false;    //most events have Q=0??
    }
    if( run == 64688 )
    {
        return false;    // sw 12 not the right file
    }
    if( run < 67000 )
    {
        return false;
    }
    if( run == 70831 || run == 70832 || run == 71162 || run == 71163 || run == 72170 || run == 72171  || run == 75637 || run == 75639 || run == 76395 || run == 76396 || run == 76792 || run == 76793 || run == 73484 || run == 73488 || run == 74220 || run == 74221 || run == 74681 || run == 74682 || run == 75536 || run == 75537 )
    {
        return false;
    }
    return true;
    
    
}


TGraph* graph_channel( TTree* slopes, int tel, int channel, TString drawstring )
{
    TString name;
    name.Form( "status[0]==0 && low.status[1]==0 && channel==%d", channel );
    slopes->Draw( drawstring.Data(), name.Data() , "goff" );
    if( slopes->GetSelectedRows() < 2 )
    {
        return 0;
    }
    TGraph* g = new TGraph( slopes->GetSelectedRows(), slopes->GetV2(), slopes->GetV1() );
    name.Form( "tel%d_c%3d", tel, channel );
    for( int i = 0; i < g->GetN(); i++ )
    {
        if( !is_good_run( g->GetX()[i], tel ) )
        {
            g->RemovePoint( i );
            i--;
        }
    }
    g->SetName( name.Data() );
    g->SetTitle( name.Data() );
    
    return g;
    
}

TGraphErrors* graph_board( TTree* slopes, int tel, int min, int n, TString drawstring, double ymin = -99, double ymax = 999 )
{

    TString name;
    
    vector<double> run;
    vector<double> mult;
    vector<double> multE;
    
    //get possible run numbers
    
    for( int i = min; i < n + min; i++ )
    {
        TGraph* temp = graph_channel( slopes, tel, i, drawstring );
        if( !temp )
        {
            /*cout << "Help! " << i << endl; */continue;
        }
        for( int j = 0; j < temp->GetN() ; j++ )
        {
            bool found = false;
            for( unsigned int k = 0; k < run.size() ; k++ )
            {
                if( run.at( k ) == temp->GetX()[j] )
                {
                    found = true;
                    break;
                }
            }
            if( !found && is_good_run( temp->GetX()[j], tel ) )
            {
                run.push_back( temp->GetX()[j] );
            }
        }
    }
    TString bla = drawstring;
    bla.ReplaceAll( ":run", "" );
    for( unsigned int r = 0; r < run.size(); r++ )
    {
        name.Form( "status[0]==0 && low.status[1]==0 && channel>%d && channel<%d && run==%.0f && m[0]/low.m[1]<20  && m[0]/low.m[1]>0 && low.m[1]>0 && %s>%f && %s<%f", min, min + n, run.at( r ), bla.Data(), ymin, bla.Data(), ymax );
        slopes->Draw( drawstring.Data(), name.Data() , "goff" );
        
        double m = TMath::Mean( slopes->GetSelectedRows(), slopes->GetV1() );
        //double m = TMath::Median ( slopes->GetSelectedRows(), slopes->GetV1() );
        double mE = TMath::RMS( slopes->GetSelectedRows(), slopes->GetV1() );
        double N = slopes->GetSelectedRows();
        
        //cout << run.at(r) << " " << m << " " << mE << " " << N << endl;
        
        if( mE / sqrt( N ) > 0.5 * m )
        {
            cout << N << " " << mE << " " << m << " " << run.at( r ) << endl;
        }
        if( mE / sqrt( N ) > 0.5 * m )
        {
            for( int jk = 0; jk < slopes->GetSelectedRows(); jk++ )
            {
                cout << jk << " " << slopes->GetV1()[jk] << endl;
            }
        }
        
        if( N < 0.4 * n )
        {
            run.at( r ) = -1 ;
        }
        if( n < 400 && N < 0.8 * N )
        {
            run.at( r ) = -1 ;
        }
        mult.push_back( m );
        multE.push_back( mE / sqrt( N ) );
    }
    
    TGraphErrors* e = new TGraphErrors( run.size(), &( run[0] ), &( mult[0] ), 0, &( multE[0] ) );
    e->Sort();
    for( int r = 0; r < e->GetN() ; r++ )
    {
        if( e->GetX()[r] < 0 )
        {
            e->RemovePoint( r );
            r--;
        }
    }
    
    name.Form( "tel%d_c%3d-c%3d", tel, min, min + n - 1 );
    e->SetName( name.Data() );
    e->SetTitle( name.Data() );
    return e;
    
}


TH1D* hist_board( TTree* slopes, int tel, int min, int n, int minrun, int maxrun, TString drawstring, double ymin = -99, double ymax = 999 )
{
    TString cut;
    cut.Form( "status[0]==0 && low.status[1]==0 && channel>=%d && channel<%d && run>=%d && run<=%d && m[0]/low.m[1]<20  && m[0]>0 && low.m[1]>0 ", min, min + n, minrun, maxrun );
    TString name;
    name.Form( "hist_tel%d_c%3d-c%3d_r%d-r%d", tel, min, min + n - 1, minrun, maxrun );
    
    TString title;
    title.Form( "Tel %d;Low gain multiplier;# Pixels", tel );
    
    TH1D* hist = new TH1D( name.Data(), title.Data(), 100, ymin, ymax );
    TString temp;
    temp.Form( "%s>>%s", drawstring.Data(), name.Data() );
    slopes->Draw( temp.Data(), cut.Data() , "" );
    
    return hist;
    
}



void plot_by_group( int tel, int sw_high, int sw_low, int min, int n, int ngroup, TString infile, TString outdir, TFile* outfile, double minY, double maxY, TString drawstring = "m[0]/low.m[1]:run", TString prefix = "" )
{
    TString name;
    name.Form( infile.Data(), sw_high, tel );
    TFile* in = new TFile( name.Data(), "read" );
    if( !in )
    {
        exit( 1 );
    }
    name.Form( "slopes_%d_%d", tel, sw_high );
    TTree* slopes = ( TTree* ) in->Get( name.Data() );
    
    name.Form( infile.Data(), sw_low, tel );
    TFile* inlow = new TFile( name.Data(), "read" );
    if( !inlow )
    {
        exit( 1 );
    }
    name.Form( "slopes_%d_%d", tel, sw_low );
    TTree* slopes_low = ( TTree* ) inlow->Get( name.Data() );
    
    slopes->AddFriend( slopes_low, "low" );
    
    TMultiGraph* mg = new TMultiGraph;
    /*	TGraphErrors * t = graph_board( slopes, tel, 0, 249, drawstring);
    	if(t)
    	{
    		t->SetLineColor( 15 );
    		t->SetLineWidth(4);
    		mg->Add(t, "L");
    	}
    	t = graph_board( slopes, tel, 250, 498, drawstring);
    	if(t)
    	{
    		t->SetLineColor( 15 );
    		t->SetLineWidth(4);
    		mg->Add(t, "L");
    	}
    */	for( int i = 0; i < ngroup; i++ )
    {
        int minc = min + i * n;
        TGraphErrors* t = graph_board( slopes, tel, minc, n, drawstring, minY, maxY );
        if( !t )
        {
            continue;
        }
        
        int color = i + 1;
        if( i > 8 )
        {
            color = i + 2;
        }
        t->SetLineColor( color );
        t->SetLineWidth( 2 );
        mg->Add( t, "EL*" );
    }
    
    
    
    TCanvas* c2 = new TCanvas( "c", "c", 1200, 800 );
    c2->cd();
    mg->Draw( "alp" );
    //	mg->GetXaxis()->SetLimits(54000,90000);
    mg->GetXaxis()->SetLimits( 67000, 80000 );
    
    mg->SetMinimum( minY );
    mg->SetMaximum( maxY );
    
    TLegend* l = c2->BuildLegend( 0.55, 0.70, 0.95, 0.95 );
    l->SetNColumns( 2 );
    name.Form( "Sumwindow %d/%d", sw_high, sw_low );
    l->SetHeader( name.Data() );
    l->SetFillColor( kWhite );
    
    TLine* line = new TLine();
    line->SetLineStyle( 7 );
    line->SetLineWidth( 2 );
    line->SetLineColor( kGray );
    line->DrawLine( 63300, minY, 63300, maxY );
    line->DrawLine( 67500, minY, 67500, maxY );
    line->DrawLine( 69500, minY, 69500, maxY ); //T1 flasher
    line->DrawLine( 72190, minY, 72190, maxY ); //T3 flasher
    line->DrawLine( 72260, minY, 72260, maxY ); //T4 flasher
    
    c2->SetGrid();
    c2->Modified();
    c2->Update();
    
    name.Form( "%s/%sboard_lmult_sw%d_%d_tel%d_c%d-%d.pdf", outdir.Data(), prefix.Data(), sw_high, sw_low, tel, min, min + n * ngroup - 1 );
    c2->SaveAs( name.Data() );
    c2->SetName( name.Data() );
    outfile->cd();
    c2->Write();
    c2->Close();
    
    in->Close();
    inlow->Close();
}


void plot_by_channel( int tel, int sw_high, int sw_low , int min, int n, TString infile, TString outdir, TFile* outfile, double minY, double maxY, TString drawstring = "m[0]/low.m[1]:run", TString prefix = "" )
{
    TString name;
    name.Form( infile.Data(), sw_high, tel );
    TFile* in = new TFile( name.Data(), "read" );
    if( !in )
    {
        exit( 1 );
    }
    name.Form( "slopes_%d_%d", tel, sw_high );
    TTree* slopes = ( TTree* ) in->Get( name.Data() );
    
    name.Form( infile.Data(), sw_low, tel );
    TFile* inlow = new TFile( name.Data(), "read" );
    if( !inlow )
    {
        exit( 1 );
    }
    name.Form( "slopes_%d_%d", tel, sw_low );
    TTree* slopes_low = ( TTree* ) inlow->Get( name.Data() );
    
    slopes->AddFriend( slopes_low, "low" );
    
    
    TMultiGraph* mg = new TMultiGraph;
    /*	TGraphErrors * t = graph_board( slopes, tel, 0, 249, drawstring);
    	if(t)
    	{
    		t->SetLineColor( 15 );
    		t->SetLineWidth(4);
    		mg->Add(t, "L");
    	}
    	t = graph_board( slopes, tel, 250, 498, drawstring);
    	if(t)
    	{
    		t->SetLineColor( 15 );
    		t->SetLineWidth(4);
    		mg->Add(t, "L");
    	}
    */
    for( int i = 0; i < n; i++ )
    {
        int start = min;
        TGraph* t = graph_channel( slopes, tel, i + start, drawstring );
        if( !t )
        {
            continue;
        }
        int color = i + 1;
        if( i > 8 )
        {
            color = i + 2;
        }
        t->SetLineColor( color );
        t->SetLineWidth( 2 );
        mg->Add( t, "LP*" );
    }
    
    TCanvas* c2 = new TCanvas( "c", "c", 1200, 800 );
    c2->cd();
    mg->Draw( "alp" );
    mg->GetXaxis()->SetLimits( 54000, 90000 );
    mg->GetXaxis()->SetLimits( 67000, 80000 );
    mg->SetMinimum( minY );
    mg->SetMaximum( maxY );
    
    
    TLegend* l = c2->BuildLegend( 0.55, 0.70, 0.95, 0.95 );
    l->SetNColumns( 2 );
    name.Form( "Sumwindow %d/%d", sw_high, sw_low );
    l->SetHeader( name.Data() );
    l->SetFillColor( kWhite );
    
    TLine* line = new TLine();
    line->SetLineStyle( 7 );
    line->SetLineWidth( 2 );
    line->SetLineColor( kGray );
    line->DrawLine( 63300, minY, 63300, maxY );
    line->DrawLine( 67500, minY, 67500, maxY );
    line->DrawLine( 69500, minY, 69500, maxY ); //T1 flasher
    line->DrawLine( 72190, minY, 72190, maxY ); //T3 flasher
    line->DrawLine( 72260, minY, 72260, maxY ); //T4 flasher
    c2->SetGrid();
    
    c2->Modified();
    c2->Update();
    
    name.Form( "%s/%schannel_lmult_sw%d_%d_tel%d_c%d-%d.pdf", outdir.Data(), prefix.Data(), sw_high, sw_low, tel, min, min + n - 1 );
    c2->SaveAs( name.Data() );
    c2->SetName( name.Data() );
    outfile->cd();
    c2->Write();
    //c2->Close();
    
    in->Close();
    inlow->Close();
}

void plot_all( int sw_high, int sw_low , TString infile, TString outdir, TFile* outfile, double minY, double maxY, TString drawstring = "m[0]/low.m[1]:run", TString prefix = "" )
{

    TMultiGraph* mg = new TMultiGraph;
    
    TString name;
    
    for( int tel = 1; tel < 5; tel++ )
    {
    
        name.Form( infile.Data(), sw_high, tel );
        TFile* in = new TFile( name.Data(), "read" );
        if( !in )
        {
            exit( 1 );
        }
        name.Form( "slopes_%d_%d", tel, sw_high );
        TTree* slopes = ( TTree* ) in->Get( name.Data() );
        
        name.Form( infile.Data(), sw_low, tel );
        TFile* inlow = new TFile( name.Data(), "read" );
        if( !inlow )
        {
            exit( 1 );
        }
        name.Form( "slopes_%d_%d", tel, sw_low );
        TTree* slopes_low = ( TTree* ) inlow->Get( name.Data() );
        
        slopes->AddFriend( slopes_low, "low" );
        
        TGraphErrors* err = graph_board( slopes, tel, 0, 250 , drawstring, minY, maxY ) ;
        err->SetLineWidth( 2 );
        err->SetLineColor( 2 * tel - 1 );
        mg->Add( err );
        
        err = graph_board( slopes, tel, 250, 249 , drawstring, minY, maxY ) ;
        err->SetLineWidth( 2 );
        err->SetLineColor( 2 * tel - 0 );
        
        mg->Add( err );
        
        in->Close();
        inlow->Close();
    }
    
    TCanvas* c2 = new TCanvas( "c", "c", 1200, 800 );
    c2->cd();
    mg->Draw( "alp" );
    mg->GetXaxis()->SetLimits( 54000, 90000 );
    mg->GetXaxis()->SetLimits( 67000, 80000 );
    mg->SetMinimum( minY );
    mg->SetMaximum( maxY );
    
    TLegend* l = c2->BuildLegend( 0.55, 0.70, 0.95, 0.95 );
    l->SetNColumns( 2 );
    name.Form( "Sumwindow %d/%d", sw_high, sw_low );
    l->SetHeader( name.Data() );
    l->SetFillColor( kWhite );
    
    TLine* line = new TLine();
    line->SetLineStyle( 7 );
    line->SetLineWidth( 2 );
    line->SetLineColor( kGray );
    line->DrawLine( 63300, minY, 63300, maxY );
    line->DrawLine( 67500, minY, 67500, maxY );
    line->DrawLine( 69500, minY, 69500, maxY ); //T1 flasher
    line->DrawLine( 72190, minY, 72190, maxY ); //T3 flasher
    line->DrawLine( 72260, minY, 72260, maxY ); //T4 flasher
    c2->SetGrid();
    
    c2->Modified();
    c2->Update();
    
    
    name.Form( "%s/%sall_lmult_sw%d_%d.pdf", outdir.Data(), prefix.Data(), sw_high, sw_low );
    c2->SaveAs( name.Data() );
    c2->SetName( name.Data() );
    outfile->cd();
    c2->Write();
    c2->Close();
    
}


void plot_hist( int sw_high, int sw_low , TString infile, TString outdir, TFile* outfile, double min, double max, TString drawstring = "m[0]/low.m[1]", int minrun = 0, int maxrun = 0, TString prefix = "" )
{

    TString name;
    name.Form( "Runs %d to %d", minrun, maxrun );
    THStack hs( "hs", name.Data() );
    TPaveText pt( 0.6, 0.65, 0.95, 0.95, "NDC NB" );
    pt.SetFillColor( kWhite );
    
    for( int tel = 1; tel < 5; tel++ )
    {
    
        name.Form( infile.Data(), sw_high, tel );
        TFile* in = new TFile( name.Data(), "read" );
        if( !in )
        {
            exit( 1 );
        }
        name.Form( "slopes_%d_%d", tel, sw_high );
        TTree* slopes = ( TTree* ) in->Get( name.Data() );
        
        name.Form( infile.Data(), sw_low, tel );
        TFile* inlow = new TFile( name.Data(), "read" );
        if( !inlow )
        {
            exit( 1 );
        }
        name.Form( "slopes_%d_%d", tel, sw_low );
        TTree* slopes_low = ( TTree* ) inlow->Get( name.Data() );
        
        slopes->AddFriend( slopes_low, "low" );
        
        TH1D* hist = hist_board( slopes, tel, 0, 499,  minrun, maxrun, drawstring, min, max );
        hist->SetLineWidth( 2 );
        hist->SetLineColor( 2 * tel );
        
        hs.Add( hist, "E" );
        name.Form( "#color[%d]{Tel %d: mean = %.3f, #sigma = %.3f, N = %.0f}", 2 * tel, tel, hist->GetMean(), hist->GetRMS(), hist->Integral( 1, hist->GetNbinsX() ) );
        pt.AddText( name.Data() );
    }
    
    TCanvas* c2 = new TCanvas( "c", "c", 1200, 800 );
    c2->cd();
    hs.Draw( "nostack" );
    
    //	TLegend * l = c2->BuildLegend(0.8, 0.65, 0.95, 0.95);
    //	l->SetFillColor(kWhite);
    pt.Draw();
    
    c2->SetGrid();
    
    c2->Modified();
    c2->Update();
    
    
    name.Form( "%s/%shist_lmult_run%d%d_sw%d_%d.pdf", outdir.Data(), prefix.Data(), minrun, maxrun, sw_high, sw_low );
    c2->SaveAs( name.Data() );
    c2->SetName( name.Data() );
    outfile->cd();
    c2->Write();
    c2->Close();
    
}

void plot_tel( int sw_high, int sw_low , TString infile, TString outdir, TFile* outfile, double minY, double maxY, TString drawstring = "m[0]/low.m[1]:run", TString prefix = "" )
{

    TMultiGraph* mg = new TMultiGraph;
    
    TString name;
    
    TCanvas* c2 = new TCanvas( "c", "c", 1200, 800 );
    c2->cd();
    TPaveText* pt = new TPaveText( 0.1, 0.9, 0.4, 0.6, "NDC" );
    pt->SetFillStyle( 0 );
    pt->AddText( "V5" );
    TPaveText* pt2 = new TPaveText( 0.5, 0.9, 0.8, 0.6, "NDC" );
    pt2->AddText( "V6" );
    
    pt2->SetFillStyle( 0 );
    for( int tel = 1; tel < 5; tel++ )
    {
    
        int color = 2 * tel;
        
        name.Form( infile.Data(), sw_high, tel );
        TFile* in = new TFile( name.Data(), "read" );
        if( !in )
        {
            exit( 3 );
        }
        name.Form( "slopes_%d_%d", tel, sw_high );
        TTree* slopes = ( TTree* ) in->Get( name.Data() );
        if( !slopes )
        {
            exit( 5 );
        }
        name.Form( infile.Data(), sw_low, tel );
        TFile* inlow = new TFile( name.Data(), "read" );
        if( !inlow )
        {
            exit( 4 );
        }
        name.Form( "slopes_%d_%d", tel, sw_low );
        TTree* slopes_low = ( TTree* ) inlow->Get( name.Data() );
        if( !slopes_low )
        {
            exit( 6 );
        }
        slopes->AddFriend( slopes_low, "low" );
        
        TGraphErrors* err = graph_board( slopes, tel, 0, 499 , drawstring, minY, maxY ) ;
        if( !err )
        {
            exit( 8 );
        }
        err->SetLineWidth( 2 );
        err->SetLineColor( color );
        mg->Add( err );
        cout << "starting fit" << endl;
        cout << err->GetN() << endl;
        if( err->GetN() > 0 )
        
        {
            //fit old pmts (before run 63300)
            TF1* f = new TF1( "f", "[0]", err->GetX()[0], 63300 ) ;
            f->SetLineColor( color );
            f->SetParameter( 0, err->GetY()[0] );
            err->Fit( f, "Q" , "", err->GetX()[0], 63300 );
            cout << "done fit" << endl;
            TString fit = TString::Format( "#color[%d]{Tel %d: %1.3f#pm%1.3f, chi2/ndf %.1f}", color, tel, f->GetParameter( 0 ), f->GetParError( 0 ), ( ( f->GetNDF() > 0 ) ? f->GetChisquare() / f->GetNDF() : f->GetChisquare() ) );
            pt->AddText( fit.Data() );
            
            logfile << "Tel " << tel << "V5 swhigh " << sw_high << " swlow " << sw_low << " lmult " << f->GetParameter( 0 ) << endl;
            
            //fit old pmts (after run 67500)
            TF1* g = new TF1( "g", "[0]", 67500, err->GetX()[err->GetN() - 1 ] + 100 ) ;
            g->SetLineColor( color );
            g->SetParameter( 0, err->GetY()[ err->GetN() - 1 ] );
            
            err->Fit( g, "Q+" , "", 67500,  err->GetX()[err->GetN() - 1 ] + 100 );
            
            fit = TString::Format( "#color[%d]{Tel %d: %1.3f#pm%1.3f, chi2/ndf %.1f}", color, tel, g->GetParameter( 0 ), g->GetParError( 0 ), ( ( f->GetNDF() > 0 ) ? f->GetChisquare() / f->GetNDF() : f->GetChisquare() ) );
            pt2->AddText( fit.Data() );
            
            logfile << "Tel " << tel << "V6 swhigh " << sw_high << " swlow " << sw_low << " lmult " << g->GetParameter( 0 ) << endl;
            in->Close();
            inlow->Close();
        }
        else
        {
            cout << "Error in plot_tel( sw_high = " << sw_high << ", sw_low = " << sw_low << " ), no good runs. Check run quality and axis range!" << endl;
            in->Close();
            inlow->Close();
            return;
        }
    }
    
    mg->Draw( "ap" );
    mg->GetXaxis()->SetLimits( 54000, 90000 );
    mg->GetXaxis()->SetLimits( 67000, 80000 );
    mg->SetMinimum( minY );
    mg->SetMaximum( maxY );
    
    TLegend* l = c2->BuildLegend( 0.55, 0.70, 0.95, 0.95 );
    l->SetNColumns( 2 );
    name.Form( "Sumwindow %d/%d", sw_high, sw_low );
    l->SetHeader( name.Data() );
    l->SetFillColor( kWhite );
    
    TLine* line = new TLine();
    line->SetLineStyle( 7 );
    line->SetLineWidth( 2 );
    line->SetLineColor( kGray );
    line->DrawLine( 63300, minY, 63300, maxY );
    line->DrawLine( 67500, minY, 67500, maxY );
    line->DrawLine( 69500, minY, 69500, maxY ); //T1 flasher
    line->DrawLine( 72190, minY, 72190, maxY ); //T3 flasher
    line->DrawLine( 72260, minY, 72260, maxY ); //T4 flasher
    
    pt->Draw();
    pt2->Draw();
    c2->SetGrid();
    
    c2->Modified();
    c2->Update();
    
    
    name.Form( "%s/%stel_lmult_sw%d_%d.pdf", outdir.Data(), prefix.Data(), sw_high, sw_low );
    c2->SaveAs( name.Data() );
    c2->SetName( name.Data() );
    outfile->cd();
    c2->Write();
    c2->Close();
    
}

void plot_lowgain_multipliers( TString infileT = "/lustre/fs9/group/cta/users/fleish/LMULT/%d/Tel_%d/T%d.lmult.root", TString outdirT = "hilographs-%d-%d", TString log = "hilo.log" )
{

    gROOT->SetBatch( true );
    
    TString plotstring = "m[0]/low.m[1]:run";
    
    logfile.open( log.Data() );
    
    int high[2] = { 6, 16 };
    int low[15] = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
    
    
    double min[2][15] = {  { 5,	6,	5,	4,	4,	3.5,	3.5,	3.5,	3,	3,	3,	3,	3,	3,	3  } ,
        { 9,	8,	6,	6,	5,	4,	4,	4,	4,	4,	4.2,	4.2,	4.6,	4,	4.1 }
    };
    
    double max[2][15] = { { 13,	10,	7,	7,	6,	6,	6,	6,	5,	5,	5,	5,	5,	5,	5  },
        { 20,	18,	12,	12,	9,	9,	9,	9,	9,	9,	7.8,	7.8,	5.8,	7,	6.6 }
    } ;
    
    
    for( int ihigh = 1; ihigh < 2; ihigh++ )
    {
        for( int ilow = 0; ilow < 15; ilow++ )
        {
        
        
            double minY = min[ihigh][ilow];
            double maxY = max[ihigh][ilow];
            
            int sw_high = high[ihigh];
            int sw_low = low[ilow];
            
            cout << "sumwindows " << sw_high << " " << sw_low << endl;
            
            
            TString outdir = TString::Format( outdirT.Data(), sw_high, sw_low );
            
            gSystem->mkdir( outdir.Data() );
            int N = 10;
            
            TString outfilename = TString::Format( "%s/canvasses.root", outdir.Data() );
            TFile* outfile = new TFile( outfilename.Data(), "recreate" );
            
            //optional: histogram of low gain multiplier distribution
            //plot_hist( sw_high, sw_low, infileT, outdir, outfile, minY, maxY, "m[0]/low.m[1]" ,78348, 78349  );
            
            
            //low gain multiplier (mean) vs run number, by telescopes/camera halves.
            plot_all( sw_high, sw_low, infileT, outdir, outfile, minY, maxY, plotstring );
            plot_tel( sw_high, sw_low, infileT, outdir, outfile, minY, maxY, plotstring );
            
            //comment out to get tons of plots.
            continue;
            for( int tel = 1; tel < 5; tel++ )
            {
            
                //optional: low gain multiplier by channel vs run number
                for( int i = 0; i < 50; i++ )
                {
                    plot_by_channel( tel, sw_high, sw_low , i * N, N, infileT, outdir, outfile, minY, maxY, plotstring );
                    
                }
                
                //optional: low gain multiplier by FADC board (mean) vs run number
                for( int i = 0; i < 5; i++ )
                {
                    plot_by_group( tel, sw_high, sw_low , i * N * N, N, N, infileT, outdir, outfile, minY, maxY, plotstring );
                }
                
            }
            outfile->Close();
        }
    }
    
    logfile.close();
}
