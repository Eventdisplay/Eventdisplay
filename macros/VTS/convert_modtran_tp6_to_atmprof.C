/* read in tp6 file; write out atmprof file.
   uses VAtmosphereSoundings class (from VAnasum.so).


Author: Henrike Fleischhack
*/


void convert_modtran_tp6_to_atmprof( string name, int number )
{

    gSystem->Load( "$EVNDISPSYS/lib/libVAnaSum.so" );
    
    gROOT->SetBatch( true );
    
    VAtmosphereSoundings v;
    if( v.read_MODTRAN_Atmosphere( name ) < 0 )
    {
        cout << "Error, could not read modtran atmosphere " << name << endl;
        exit( 99 );
    }
    v.write_CORSIKA_UserProfile( 0, number, name );
    
    exit( 0 );
    
}
