/*! \file evndisp.cpp
    \brief eventdisplay

    event analysis and display

*/

#include <TApplication.h>
#include <TGClient.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include <string>
#include <vector>

#include "VDisplay.h"
#include "VEventLoop.h"
#include "VReadRunParameter.h"

using namespace std;

// fitter for log likelihood, has to be global
TMinuit* fLLFitter;

int main( int argc, char* argv[] )
{
    // some timing
    TStopwatch fStopWatch;
    fStopWatch.Start();
    
    // print version only
    if( argc == 2 )
    {
        string fCommandLine = argv[1];
        if( fCommandLine == "-v" || fCommandLine == "--version" )
        {
            VGlobalRunParameter fRunPara;
            cout << fRunPara.getEVNDISP_VERSION() << endl;
            exit( 0 );
        }
    }
    
    // read the command line parameters
    VReadRunParameter* fReadRunParameter = new VReadRunParameter();
    if( !fReadRunParameter->readCommandline( argc, argv ) )
    {
        exit( EXIT_FAILURE );
    }
    fReadRunParameter->getRunParameter()->print();
    
    // initialize main loop
    VEventLoop mainEventLoop( fReadRunParameter->getRunParameter() );
    if( !mainEventLoop.initEventLoop() )
    {
        exit( EXIT_FAILURE );
    }
    
    // no display, command line mode
    if( !fReadRunParameter->getRunParameter()->fdisplaymode )
    {
        mainEventLoop.loop( fReadRunParameter->getRunParameter()->fnevents );
        fStopWatch.Stop();
        fStopWatch.Print();
        mainEventLoop.shutdown();
    }
    // display mode
    else
    {
        // number of options set to one, otherwise TApplication prints sometimes help text (don't know how to switch that of in ROOT)
        Int_t targv = 1;
        TApplication app( "app", &targv, argv );
        VDisplay display( gClient->GetRoot(), fReadRunParameter->getRunParameter()->fw, fReadRunParameter->getRunParameter()->fh, &mainEventLoop );
        display.Draw();
        app.Run();
        
    }
    delete fReadRunParameter;
    
    return 0;
}
