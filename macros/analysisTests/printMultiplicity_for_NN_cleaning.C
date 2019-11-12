/*
 * print probability for trigger multiplicity for optimised next-neighbour cleaning
 *
 * - during evndisp run only one nn-group should be activated
 * - expect that only one telescope type is analysed
 *
 *   input files are eventdisplay root files
*/

#include <sstream>
#include <string>

void printMultiplicity( string iFile, float iMultiplicity = 2. )
{
     // evndisp file
     TFile *fData = new TFile( iFile.c_str() );
     if( fData->IsZombie() )
     {
         return;
     }

     TTree *s = (TTree*)fData->Get( "showerpars" );
     if( !s )
     {
        cout << "Error: showerpars tree not found" << endl;
        return;
     }
     UInt_t NTel = 0;
     s->SetBranchAddress( "NTel", &NTel );
     s->GetEntry( 0 );

     cout << "Total number of telescopes: " << NTel << endl;

     // counters
     float ntot = 0.;
     vector< float > nm;
     for( unsigned int i = 0; i < 10; i++ )
     {
          nm.push_back( 0. );
     }

     // loop over all telescopes
     for( unsigned int i = 0; i < NTel; i++ )
     {
         ostringstream iAs;
         iAs << "Tel_" << i+1 << "/tpars";

         // get tpars tree
         TTree *t = (TTree*)fData->Get( iAs.str().c_str() );
         if( !t )
         {
             cout << "Error: showerpars tree " << iAs.str() << " not found" << endl;
             continue;
         }
         UShort_t ntubes;
         t->SetBranchAddress( "ntubes", &ntubes );

         cout << "Reading tree " << iAs.str() << " (" << t->GetEntries() << " entries) " << endl;

         // loop over all events in tpars tree
         for( int j = 0; j < t->GetEntries(); j++ )
         {
             t->GetEntry( j );

             ntot++;
             for( unsigned int m = 0; m < nm.size(); m++ )
             {
                  if( ntubes == m )
                  {
                       nm[m]++;
                  }
             }
         }


     }
     // summary
     cout << "Total number of events: " << ntot << endl;
     if( ntot > 0 )
     {
         for( unsigned int m = 0; m < nm.size(); m++ )
         {
              cout << "n" << m << ": " << nm[m];
              cout << ", fraction: " << nm[m]/ntot*100. << "%" << endl;
         }
         if( iMultiplicity > 0. )
         {
             cout << "Total fraction for multiplicity " << iMultiplicity << ": ";
             float iFR = 0.;
             for( unsigned int m = 0; m < nm.size(); m++ )
             {
                  iFR += nm[m]/ntot*100. * (float)m / iMultiplicity;
             }
             cout << iFR << "%" << endl;
         } 
     }


}
