/*! \file testEvndispOutput.cpp
 *
 *  test eventdisplay DL1 file for completeness:
 *  - telconfig tree: test telescope types and consistency with
 *    expectation for the given CTA 
 *  - showerpars tree
 *  - all expected tpars trees
 *
 * 
 */

#include <iostream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TError.h"
#include "TTree.h"

using namespace std;

void printHelp()
{
     cout << endl;
     cout << "./testEvndispOutput : test an evndisp file for all required trees and contents" << endl;
     cout << endl;
     cout << "to write a log file into a root file:" << endl;
     cout << "\t ./testEvndispOutput <evndisp root file (DL1)> [CTA production]" << endl;
     cout << endl;
     cout << "\t hardwired CTA productions: prod5-South, prod5-North, prod6-South, prod6-North" << endl;
     cout << endl;
     cout << "Important: checks for minimal number of required telescopes for a productions; ignores any additional telescopes" << endl;
     cout << endl;
}

bool checkTelType( ULong64_t teltype,
                   ULong64_t teltypematch,
                   string rootfile,
                   int telid )
{
   if( teltype == teltypematch )
   {
       return true;
   }

   cout << "Error: wrong telescope type for ";
   cout << teltype << " (tel " << telid << "): ";
   cout << rootfile << endl;
   exit( EXIT_FAILURE );

   return false;
}


int main( int argc, char* argv[] )
{
     // root file to be written to / read from
     string fRootFile;
     // CTA production
     string fProduction = "prod5-South";

     if( argc != 2 && argc != 3 )
     {
         printHelp();
         exit( EXIT_SUCCESS );
     } 

     fRootFile = argv[1];
     if( argc == 3 )
     {
         fProduction = argv[2];
     }

     // ignore all warnings
     // "Warning in <TClass::Init>: no dictionary for class .. available..
     gErrorIgnoreLevel = kError;

     ////////////////////////////////////////////////
     TFile fF( fRootFile.c_str(), "READ" );
     if( fF.IsZombie() )
     {
        cout << "Error: root file not found: " << fRootFile << endl;
        exit( EXIT_FAILURE );
     }
     // requirement 1: telconfig tree exist
     TTree *telconfig = (TTree*)fF.Get( "telconfig" );
     if( !telconfig )
     {
        cout << "Error: telconfig tree not found: " << fRootFile << endl;
        exit( EXIT_FAILURE );
     }
     int TelID = 0;
     ULong64_t TelType = 0;
     telconfig->SetBranchAddress( "TelID", &TelID );
     telconfig->SetBranchAddress( "TelType", &TelType );
     if( fProduction == "prod5-South" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );

             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 28 ) checkTelType( TelType, 10408618, fRootFile, TelID ); 
             else if( TelID >= 29 && TelID <= 98 ) checkTelType( TelType, 201409917, fRootFile, TelID );
             else if( TelID >= 99 && TelID <= 123 ) checkTelType( TelType, 10608418, fRootFile, TelID );
             else if( TelID >= 124 && TelID <= 126 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 127 && TelID <= 129 ) checkTelType( TelType, 10608418, fRootFile, TelID );
             else if( TelID >= 130 && TelID <= 170 ) checkTelType( TelType, 201409917, fRootFile, TelID );
         }
     }
     else if( fProduction == "prod5-South-BL-MSTF" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );

             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 28 ) checkTelType( TelType, 10408618, fRootFile, TelID ); 
             else if( TelID >= 29 && TelID <= 98 ) checkTelType( TelType, 201409917, fRootFile, TelID );
          }
     }
     else if( fProduction == "prod5-South-BL-MSTN" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );

             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 73 ) checkTelType( TelType, 201409917, fRootFile, TelID );
             else if( TelID >= 74 && TelID <= 98 ) checkTelType( TelType, 10608418, fRootFile, TelID ); 
          }
     }
     else if( fProduction == "prod5-North" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );

             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 18 ) checkTelType( TelType, 10608418, fRootFile, TelID ); 
             else if( TelID >= 19 && TelID <= 33 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 34 && TelID <= 46 ) checkTelType( TelType, 10608418, fRootFile, TelID );
             else if( TelID >= 47 && TelID <= 59 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 60 && TelID <= 62 ) checkTelType( TelType, 10608418, fRootFile, TelID );
             else if( TelID >= 63 && TelID <= 65 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 66 && TelID <= 74 ) checkTelType( TelType, 10608418, fRootFile, TelID );
             else if( TelID >= 75 && TelID <= 83 ) checkTelType( TelType, 10408618, fRootFile, TelID );
         }
     }
     else if( fProduction == "prod3b-paranal-SCT156Tel" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );

             if( TelID >= 0 && TelID <= 40 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 41 && TelID <= 114 ) checkTelType( TelType, 201309316, fRootFile, TelID );
             else if( TelID >= 114 && TelID <= 154 ) checkTelType( TelType, 207308707, fRootFile, TelID );
         }
     }
     else if( fProduction == "prod6-North" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );
             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 39 ) checkTelType( TelType, 10608418, fRootFile, TelID );
             else
             {
                 cout << "Error: unexpected telescope of type " << TelType;
                 cout << ", telID: " << TelID;
                 cout << " (" << fRootFile << ")" << endl;
                 exit( EXIT_FAILURE );
             }
         }
     }
     else if( fProduction == "prod6-South" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );
             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 36 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 37 && TelID <= 86 ) checkTelType( TelType,201109916, fRootFile, TelID );
             else
             {
                 cout << "Error: unexpected telescope of type " << TelType;
                 cout << ", telID: " << TelID;
                 cout << " (" << fRootFile << ")" << endl;
                 exit( EXIT_FAILURE );
             }
         }
     }
     else if( fProduction == "prod6-SouthSCT" )
     {
         for( unsigned int i = 0; i < telconfig->GetEntries(); i++ )
         {
             telconfig->GetEntry( i );
             if( TelID >= 0 && TelID <= 3 ) checkTelType( TelType, 138704810, fRootFile, TelID );
             else if( TelID >= 4 && TelID <= 36 ) checkTelType( TelType, 10408618, fRootFile, TelID );
             else if( TelID >= 37 && TelID <= 86 ) checkTelType( TelType,201109916, fRootFile, TelID );
             else if( TelID >= 87 && TelID <=105 ) checkTelType( TelType, 205008707, fRootFile, TelID );
             else
             {
                 cout << "Error: unexpected telescope of type " << TelType;
                 cout << ", telID: " << TelID;
                 cout << " (" << fRootFile << ")" << endl;
                 exit( EXIT_FAILURE );
             }
         }
     }
             
     // requirement 2: showerpars tree exist
     TTree *showerpars = (TTree*)fF.Get( "showerpars" );
     if( !showerpars )
     {
        cout << "Error: showerpars tree not found: " << fRootFile << endl;
        exit( EXIT_FAILURE );
     }
     Long64_t nentries = showerpars->GetEntries();

     // requirement 3: check if all tpars tree exist
     int ntel = (int)telconfig->GetEntries();
     for( int i = 0; i < ntel; i++ )
     {
         stringstream str_tel;
         str_tel << "Tel_" << i+1 << "/tpars";
         TTree *tpars = (TTree*)fF.Get( str_tel.str().c_str() );
         if( !tpars )
         {
            cout << "Error: tpars tree " << str_tel.str() << " not found: " << fRootFile << endl;
            exit( EXIT_FAILURE );
         }
         if( nentries != tpars->GetEntries() )
         {
            cout << "Error: tpars entries different to showerpars (Expected: ";
            cout << nentries << ", found " << tpars->GetEntries() << ": ";
            cout << fRootFile << endl;
         }
     }
     cout << "FILE OK: " << fRootFile << endl;
}
