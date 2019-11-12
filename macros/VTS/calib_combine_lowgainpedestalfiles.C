#include <fstream>
/*
Macro to combine low gain pedestal files from 2 calib runs to one common file per telescope, since each low gain calib run only has half the camera in low gain mode.

Author: Henrike Fleischhack

To be used with CINT.

*/
void calib_combine_lowgainpedestalfiles( TString file = "", TString calibdir = "" )
{
    gSystem->Load( "$EVNDISPSYS/lib/libVAnaSum.so" );
    if( file == "" )
    {
        cout << "Usage: root calib_combine_lowgainpedestalfiles\( TString list [, TString calibdir ] \)" << endl ;
        cout << "list should be the filename of a textfile containing two low gain calib runs per line " << endl;
        cout << "calibdir is the calibration directory, by default: $VERITAS_EVNDISP_AUX_DIR/Calibration" << endl;
        exit( 1 );
    }
    
    if( calibdir == "" )
    {
        calibdir = gSystem->Getenv( "VERITAS_EVNDISP_AUX_DIR" ) + "/Calibration/";
    }
    
    int ntel = 4;
    ifstream list( file.Data() );
    TString run1, run2;
    VPedestalLowGain a;
    while( list >> run1 >> run2 )
    {
        TString name = run1;
        name.Append( run2( 3, 4 ).Data() );
        cout << name << endl;
        a.combineLowGainPedestalFileForAllTelescopes( ntel , calibdir.Data(), run1.Data() , run2.Data() , name.Data() ) ;
    }
}

