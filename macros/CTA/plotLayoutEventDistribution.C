/* plot event distribution and direction reconstruction quality accross the array layout
 *
 * requires as input an effective area file produced with the full data tree feature
 *
 * - hard wired telescope types for South FG (prod3b)
 */

void plot( double iMaxCoreDistance = 1000., double iMinEnergy = 1.,
           string iEffFileFullDataTree = "gamma_onSource.S.3HB9-P1-A1L-FG_ID4.eff-0.root",
           string iTelConfigTreeFile = "../../../S.3HB9-TS-BB-FG/gamma_onSource/99HD_287_180deg.root" )
{
    // telconfig tree with telescope positions
     TFile *fTelConfig = 0;
     TTree *telconfig = 0;
     if( iTelConfigTreeFile.size() > 0 )
     {
        fTelConfig = new TFile( iTelConfigTreeFile.c_str() );
        if( !fTelConfig->IsZombie() )
        {
            telconfig = (TTree*)fTelConfig->Get( "telconfig" );
        }
     }
          


     TFile *fG = new TFile( iEffFileFullDataTree.c_str() );
     if( !fG )
     {
         return;
     }
     TTree *data = (TTree*)gDirectory->Get( "data" );
     if( !data )
     {
         return;
     }
     data->AddFriend( "fEventTreeCuts" );


     // plot core positions
     TH2D *hCore = new TH2D( "hCore", "", 50., -1.*iMaxCoreDistance, iMaxCoreDistance , 50., -1*iMaxCoreDistance, iMaxCoreDistance);
     hCore->SetTitle( "" );
     hCore->SetStats( 0 );
     hCore->SetXTitle( "core position x (m)" );
     hCore->SetYTitle( "core position y (m)" );
     hCore->GetYaxis()->SetTitleOffset( 1.5 );
     hCore->GetXaxis()->SetNdivisions( 505 );
     hCore->GetYaxis()->SetNdivisions( 505 );

     char hname[300];
     sprintf( hname, "(fEventTreeCuts.Class==5||fEventTreeCuts.Class==0)&&MCe0>%f", iMinEnergy );
     data->Project( "hCore", "Ycore:Xcore", hname );

     TCanvas *cCore = new TCanvas( "cCore", "core position", 100, 100, 1300, 1200 );
     cCore->SetRightMargin( 0.15 );
     cCore->SetLeftMargin( 0.12 );
     cCore->Draw();

     hCore->Draw( "colz" );

     if( telconfig )
     {
          telconfig->SetMarkerStyle( 25 );
          telconfig->SetMarkerSize( 1.3 );
          telconfig->Draw( "TelY:TelX", "TelType==10408618", "same" );
          telconfig->SetMarkerStyle( 24 );
          telconfig->Draw( "TelY:TelX", "TelType==201309316", "same" );
     }

     /////////////////////////
     // average direction error
     //
     //
     TCanvas *cAngRes = new TCanvas( "cAngRes", "angular resolution ", 500, 100, 1300, 1200 );
     cAngRes->SetRightMargin( 0.15 );
     cAngRes->SetLeftMargin( 0.12 );
     cAngRes->Draw();

     TProfile2D *hAngRes = new TProfile2D( "hAngRes", "", 50., -1.*iMaxCoreDistance, iMaxCoreDistance , 50., -1*iMaxCoreDistance, iMaxCoreDistance, 0., 2. );
     hAngRes->SetStats( 0 );
     hAngRes->SetXTitle( "core position x (m)" );
     hAngRes->SetYTitle( "core position y (m)" );
     hAngRes->GetYaxis()->SetTitleOffset( 1.5 );
     hAngRes->GetXaxis()->SetNdivisions( 505 );
     hAngRes->GetYaxis()->SetNdivisions( 505 );

     data->Project( "hAngRes", "sqrt(Xoff*Xoff+Yoff*Yoff):Ycore:Xcore", hname );

     hAngRes->SetMaximum( 0.1 );
     hAngRes->Draw( "colz" );

     if( telconfig )
     {
          telconfig->SetMarkerStyle( 25 );
          telconfig->SetMarkerSize( 1.3 );
          telconfig->Draw( "TelY:TelX", "TelType==10408618", "same" );
          telconfig->SetMarkerStyle( 24 );
          telconfig->Draw( "TelY:TelX", "TelType==201309316", "same" );
     }


}


