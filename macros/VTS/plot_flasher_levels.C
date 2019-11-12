/*
Macro to plot monitor levels used for relative gain calculation.
Need to have analysed flasher run with -writeextracalibtree option enabled.

Call macro like this:

root -l -q -b "$EVNDISPSYS/macros/VTS/plot_flasher_levels.C( <runnumber>, \"$EVNDISPDATA/Calibration\", \"<outputdir>\" )"

(Assuming your .gain.root files are in the usual place. If not, the second argument will be the calibration directory.)

*/

void plot_flasher_levels( int run, TString calibdir, TString outdir )
{
    int nTel = 4;
    int nEvents = 0;
    TCanvas* c = new TCanvas( "c", "c", 1000, 1000 );
    c->Divide( 2, 2 );
    for( int iTel = 0; iTel < nTel; iTel++ )
    {
        c->cd( iTel + 1 );
        TString filename = TString::Format( "%s/Tel_%d/%d.gain.root", calibdir.Data(), iTel + 1, run );
        TFile* file = new TFile( filename.Data() , "read" );
        if( !file )
        {
            continue;
        }
        TString treename = TString::Format( "charges_%d", iTel + 1 );
        TTree* tree = ( TTree* )file->Get( treename.Data() );
        if( !tree )
        {
            continue;
        }
        TString titlestring = TString::Format( "Run %d, T%d;eventNumber;monitor charge", run, iTel + 1 );
        TString namestring = TString::Format( "hist_%d", iTel );
        TString plotstring = TString::Format( "QMon:eventNumber>>hist_%d", iTel );
        TH2D* bla = new TH2D( namestring.Data(), titlestring.Data(), 1000, 0, 36500, 900, 50, 950 );
        bla->SetBit( TH1::kCanRebin );
        bla->GetYaxis()->SetTitleOffset( 1.3 );
        TString cutstring = TString::Format( "eventNumber != 99999999 && %d==%d", iTel + 1, iTel + 1 );
        nEvents += tree->Draw( plotstring.Data(), cutstring.Data() );
        
    }
    TString savefile = TString::Format( "%s/%d.lightlevels.png",  outdir.Data(), run );
    if( nEvents > 0 )
    {
        c->SaveAs( savefile.Data() );
    }
    else
    {
        cout << "plot_flasher_levels warning: No flasher events." << endl;
    }
    return;
}
