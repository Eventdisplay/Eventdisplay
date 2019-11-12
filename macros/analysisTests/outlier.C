/*
 *  test outlier detection for energy reconstruction in
 *  VTableCalculator
 *
 *  Goal is to remove telescopes which gives completely
 *  different energy estimates than other telescopes
 *
 *  Try different tests:
 *    - Iglewicz and Hoaglin Z-Score
 *    - Dixon's Q-score
 *
 *   Summary:
 *   - very litte advance compared to the current methods
 *     using either the median energy, or (as the lookup
 *     tables do) use the average weighted by the expected
 *     uncertainty read from a lookup table
 */

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TPad.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

/*
 * return limits for Dixon's Q-Test
 * from https://en.wikipedia.org/wiki/Dixon%27s_Q_test
 *
 * confidenceValue as integer
 * (addeded here only 90% values, see paper or wikipedia for 
 * 95 ad 99% CL values.
 */
double DixonsQTestLimits( unsigned int n_values, int confidenceValue )
{
       if( confidenceValue == 90 )
       {
           switch( n_values )
           {
               case 3:
                  return 0.941;
                  break;
               case 4:
                  return 0.765;
                  break;
               case 5:
                  return 0.642;
                  break;
               case 6:
                  return 0.560;
                  break;
               case 7:
                  return 0.507;
                  break;
               case 8:
                  return 0.468;
                  break;
               case 9:
                  return 0.437;
                  break;
               case 10:
                  return 0.412;
                  break;
               default:
                  return 0.;
           }
      }
      return 0.;
}


/*
 *   input file: mscw_energy file (chaining allowed)
 *
 */
void outlier( string iFile, int iN_Cut = 5 )
{
      TChain *data = new TChain( "data" );
      int iN = data->Add( iFile.c_str() );
      if( iN == 0 ) 
      {
          cout << "Error: no files in chain" << endl;
          return;
      }

      // define tree variables
      Double_t MSCW = 0.;
      Int_t    NImages = 0;
      Double_t ErecS = 0.;
      Double_t MCe0 = 0.;
      Double_t E[100];
      UInt_t   NImages_Ttype[100];

      data->SetBranchAddress( "NImages", &NImages );
      data->SetBranchAddress( "MSCW", &MSCW );
      data->SetBranchAddress( "ErecS", &ErecS );
      data->SetBranchAddress( "MCe0", &MCe0 );
      data->SetBranchAddress( "ES", &E );
      data->SetBranchAddress( "NImages_Ttype", NImages_Ttype );

      /////////////////////////////////////////
      // variables used in energy calculation
      
      // energy per telescope
      vector< double > energy;
      // absolute deviation from media
      vector< double > energyMA;
      // median energy
      double median = 0.;
      // median absolute deviation
      double medianE = 0.;

      // histogram Erec/MCe0 for current mscw_energy reconstruction
      // of energies
      TH1D *hErec_0 = new TH1D( "hErec_0", "", 100., 0., 4. );
      hErec_0->SetFillColor( 9 );
      // histogram Erec/median for reconstruction without any cuts
      TH1D *hErec_1 = new TH1D( "hErec_1", "", 100., 0., 4. );
      hErec_1->SetLineStyle( 2 );
      TH1D *hErec_1_mult = new TH1D( "hErec_1_mult", "", 10., 0., 10. );
      hErec_1_mult->SetLineStyle( 2 );

      // Iglewicz and Hoaglin Z-Score
      // (see http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm)
      TH1D *hZscore = new TH1D( "hZscore", "", 200., -10., 10. );
      hZscore->SetLineWidth( 2 );

      // Dixon's Q-Score vs reconstruction quality
      TH2D *hQscore_l = new TH2D( "hQscore_l", "", 100, 0., 4., 100, 0., 1. );
      hQscore_l->SetXTitle( "E_{rec}/E_{MC}" );
      hQscore_l->SetYTitle( "Dixon's Q-Score (low)" );
      TH2D *hQscore_h = new TH2D( "hQscore_h", "", 100, 0., 4., 100, 0., 1. );
      hQscore_h->SetXTitle( "E_{rec}/E_{MC}" );
      hQscore_h->SetYTitle( "Dixon's Q-Score (high)" );
      TH2D *hQscore_h_ES = new TH2D( "hQscore_h_ES", "", 100, 0., 4., 100, 0., 1. );
      hQscore_h_ES->SetXTitle( "E_{rec}/E_{MC}" );
      hQscore_h_ES->SetYTitle( "Dixon's Q-Score (per energy, removed)" );

      // histograms for different cuts on Iglewicz and Hoaglin Z-Score
      char hname[200];
      vector< TH1D* > hErec_2;
      vector< TH1D* > hErec_2_mult;
      vector< vector< double > > ErecS_2;
      vector< double > fZ;
      fZ.push_back( 5. );
      fZ.push_back( 3. );
      fZ.push_back( 1. );
      fZ.push_back( 0.5 );
      fZ.push_back( 0.2 );

      vector< double > itemp;
      for( unsigned int i = 0; i < fZ.size(); i++ )
      {
          sprintf( hname, "hErec_2_%d", i );
          hErec_2.push_back( new TH1D( hname, "", 100., 0., 4. ) );
          hErec_2.back()->SetLineColor( i+1 );
          hErec_2.back()->SetLineWidth( 2 );

          sprintf( hname, "hErec_2_mult_%d", i );
          hErec_2_mult.push_back( new TH1D( hname, "", 10., 0., 10 ) );
          hErec_2_mult.back()->SetLineColor( i+1 );
          hErec_2_mult.back()->SetLineWidth( 2 );

          ErecS_2.push_back( itemp );
      }

      /////////////////////////////////////////////////
      // loop over all entries in data tree
      for( unsigned int i = 0; i < data->GetEntries(); i++ )
      {
           data->GetEntry( i );

           // some basic cuts
           if( NImages >= 2
            && NImages_Ttype[0] == 1 && NImages_Ttype[1] == 4   // looking for outliers by MSTs
            && MSCW > -2.&& MSCW < 1. )
//            && ErecS > 0. && MCe0 < 0.1 )
            {
                  energy.clear();
                  energyMA.clear();

                  // fill energy vector
                  for( int n = 0; n < NImages; n++ )
                  {
                       if( E[n] > 0. )
                       {
                           energy.push_back( E[n] );
                       }
                  }
                  if( energy.size() > 1 )
                  {
                        // calculate median and median absolute deviation
                        median = TMath::Median( energy.size(), &energy[0] );
                        ///////////////////////////////////////
                        // sort variables for Dixon's Q test
                        sort( energy.begin(), energy.end() );
                        // calculate range
                        double range = energy.back() - energy[0];
                        double q_low = 0.;
                        double q_hig = 0.;
                        if( range > 0. )
                        {
                            q_low = (energy[1] - energy[0])/range;
                            if( energy.size() > 1 )
                            {
                                q_hig = (energy[energy.size()-1] - energy[energy.size()-2])/range;
                            }
                            hQscore_l->Fill( ErecS/MCe0, q_low );
                            hQscore_h->Fill( ErecS/MCe0, q_hig );
                        }
                        // apply Dixon's test
                        if( q_hig > DixonsQTestLimits( energy.size(), 90 ) )
                        {
                             hQscore_h_ES->Fill( energy.back()/median, q_hig );
                             cout << "REMOVED " << MCe0 << " median " << median << " e removed " << energy.back();
                             energy.pop_back();
                             // calculate median and median absolute deviation
                             median = TMath::Median( energy.size(), &energy[0] );
                             cout << " new median " << median << " ERECS " << ErecS << endl;
                        }
                        
//                        median = MCe0;
                        for( unsigned int e = 0; e < energy.size(); e++ )
                        {
                            energyMA.push_back( TMath::Abs( energy[e] - median ) );
                        }    
                        medianE = TMath::Median( energyMA.size(), &energyMA[0] );

                        if( medianE <= 0. )
                        {
                            continue;
                        }

                             

                        // reset variables
                        for( unsigned int p = 0; p < fZ.size(); p++ )
                        {
                              ErecS_2[p].clear();
                        }

                        /////////////////////////////////////////////////////////
                        // calculate energies, applying different outlier cuts
                        for( unsigned int e = 0; e < energy.size(); e++ )
                        {
                             if( energy.size() <= 2 ) continue;

                             // note: use absolute value here
                             double zscore = 0.6745 * (energy[e] - median) / medianE;
                             hZscore->Fill( zscore );

                             zscore = TMath::Abs( zscore );

                             // check different cuts on zcore
                             for( unsigned int p = 0; p < fZ.size(); p++ )
                             {
                                 if(  zscore < fZ[p] )
                                 {
                                     ErecS_2[p].push_back( energy[e] );
                                 }
                             }
                         }
                  }
                  if( ErecS > 0. && NImages > 2 )
                  {
                       hErec_0->Fill( ErecS / MCe0 );
                  }
                  if( ErecS > 0. && NImages > 2 && median > 0. )
                  {
                       hErec_1->Fill( ErecS / median );
                       hErec_1_mult->Fill( (double)energy.size() );
                  }
                  for( unsigned int p = 0; p < fZ.size(); p++ )
                  {
                       if( ErecS_2[p].size() > 0 && energy.size() > 2 )
                       {
                           double iE_median = TMath::Median( ErecS_2[p].size(), &ErecS_2[p][0] );
                           hErec_2[p]->Fill( iE_median / MCe0 );
                       }
                       hErec_2_mult[p]->Fill( (double)ErecS_2[p].size() );
                  }
            }
       }

       ////////////////////////////////////////
       // plot everything

       TCanvas *c = new TCanvas( "cZScore", "outlier detection", 10, 10, 1400, 800 );
       c->Divide( 3, 2 );

       TPad *c1 = (TPad*)c->cd( 1 );
       c1->SetLogy( 1 );
       c1->SetGridx( 0 );
       c1->SetGridy( 0 );
       hZscore->Draw();

       /* plot Erec/MCe0 for the different
        * cuts on zscore
        */
       TPad *c2 = (TPad*)c->cd( 2 );
       c2->SetLogy( 1 );
       c2->SetGridx( 0 );
       c2->SetGridy( 0 );

       hErec_0->Draw();
       hErec_1->Draw( "same" );
       cout << "Erec, mscw_energy (" << hErec_0->GetEntries() << ") " << hErec_0->GetMean() << " +- " << hErec_0->GetRMS() << endl;
       for( unsigned int p = 0; p < fZ.size(); p++ )
       {
           hErec_2[p]->Draw( "sames" );
           cout << "Erec2 ( " << fZ[p] << ", " << hErec_2[p]->GetEntries() << "): " << hErec_2[p]->GetMean() << " +- " << hErec_2[p]->GetRMS() << endl;
       }
      
       /* plot the number of telescopes contributing 
        * to the energy estimation
        */
       TPad *c3 = (TPad*)c->cd( 3 );
       c3->SetLogy( 0 );
       c3->SetGridx( 0 );
       c3->SetGridy( 0 );

       hErec_1_mult->Draw();
       for( unsigned int p = 0; p < fZ.size(); p++ )
       {
           hErec_2_mult[p]->Draw( "sames" );
       }

       TPad *c4 = (TPad*)c->cd( 4 );
       c4->SetLogy( 0 );
       c4->SetGridx( 0 );
       c4->SetGridy( 0 );

       hQscore_l->Draw( "colz" );

       TPad *c5 = (TPad*)c->cd( 5 );
       c5->SetLogy( 0 );
       c5->SetGridx( 0 );
       c5->SetGridy( 0 );

       hQscore_h->Draw( "colz" );

       TPad *c6 = (TPad*)c->cd( 6 );
       c6->SetLogy( 0 );
       c6->SetGridx( 0 );
       c6->SetGridy( 0 );

       hQscore_h_ES->Draw( "colz" );

}
