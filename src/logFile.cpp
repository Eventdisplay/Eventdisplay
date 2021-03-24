/*! \file logFile.cpp
 *
 * write and read a log files into a root file
 *
 * can be used to lower the number of files produced
 * in large scale simulation productions
 *
 */

#include <iostream>
#include <string>

#include "TFile.h"
#include "TError.h"
#include "TMacro.h"

using namespace std;

void printHelp()
{
     cout << endl;
     cout << "./logFile: write or read a log file into / from a root file" << endl;
     cout << endl;
     cout << "to write a log file into a root file:" << endl;
     cout << "\t ./logFile <log name> <root file> <log file>" << endl;
     cout << endl;
     cout << "to read a log file from a root file: " << endl;
     cout << "\t ./logFile <log name> <root file>" << endl;
     cout << endl;
     cout << "\t examples are: " << endl;
     cout << "\t\t convLog, evndispLog" << endl;
     cout << "\t\t makeTableLog, makeTableFileList" << endl;
     cout << "\t\t smoothTableLog, mscwTableLog, mscwTableList" << endl;
     cout << "\t\t tmvaLog, tmvaRunparameter" << endl;
     cout << "\t\t effAreaLog, effAreaCuts, effAreaParameters" << endl;
     cout << "\t\t IRFLog" << endl;
     cout << endl;
}

int main( int argc, char* argv[] )
{
     // root file to be written to / read from
     string fRootFile;
     // log file to be written into root file
     string fLogFile;
     // log file name
     string fLogFileName;

     if( argc != 3 && argc != 4 )
     {
         printHelp();
         exit( EXIT_SUCCESS );
     } 

     // list of prefered logs to print
     vector< string > logObjectNames;
     logObjectNames.push_back( "evndispLog" );
     logObjectNames.push_back( "makeTableLog" );
     logObjectNames.push_back( "mscwTableLog" );
     logObjectNames.push_back( "tmvaLog" );
     logObjectNames.push_back( "effAreaLog" );
     logObjectNames.push_back( "IRFLog" );

     fLogFileName = argv[1];
     fRootFile = argv[2];

     // ignore all warnings
     // "Warning in <TClass::Init>: no dictionary for class .. available..
     gErrorIgnoreLevel = kError;

     ////////////////////////////////////////////////
     if( argc == 3 )
     {
           TFile fF( fRootFile.c_str(), "READ" );
           if( fF.IsZombie() )
           {
                cout << "Error: root file not found" << endl;
                exit( EXIT_FAILURE );
           }
           TMacro *iM = (TMacro*)fF.Get( fLogFileName.c_str() );
           // xml requires dedicated return if not found
           if( !iM && fLogFileName.find( "XML" ) != string::npos )
           {
               cout << "NOXML" << endl;
               exit( EXIT_SUCCESS );
           }
           if( !iM )
           {
               for( unsigned int i = 0; i < logObjectNames.size(); i++ )
               {
                    iM = (TMacro*)fF.Get( logObjectNames[i].c_str() );
                    if( iM )
                    {
                        break;
                    }
               }
           }
           if( iM )
           {
                iM->Print();
           }
           else
           {
                cout << "Error: log file object with name " << fLogFileName << " not found" << endl;
                exit( EXIT_FAILURE );
           }
     }
     else if( argc == 4 )
     {
           fLogFile = argv[3];
           TFile fF( fRootFile.c_str(), "update" );
           if( fF.IsZombie() )
           {
                cout << "Error: root file not found: " << fRootFile << endl;
                exit( EXIT_FAILURE );
           }
           TMacro *iM = new TMacro( fLogFile.c_str(), fLogFileName.c_str() );
           if( iM )
           {
               if( iM->GetListOfLines() && iM->GetListOfLines()->GetSize() > 0 )
               {
                   iM->Write( fLogFileName.c_str() );
               }
               else
               {
                   cout << "Error: log file not found: " << fLogFile << endl;
                   exit( EXIT_FAILURE );
               }
           }
           fF.Close();
     }
}
