/* \file convertSensitivityFilesToFITS.cpp
 *
 * convert a IRF ROOT file into a FITS file
 *
 */

#include "VFITSIRFs.h"

using namespace std;

int main( int argc, char* argv[] )
{
    string fFitsFileName = "test.fits";
    string fRootFile = "DESY.g20210610.V3.ID0NIM3LST3MST3SST4SCMST3.prod5-Paranal-20deg-sq10-LL.S-M6C5ax-14MSTs37SSTs-MSTF.180000s.root";
    TFile *fData = new TFile( fRootFile.c_str() );
    if( fData->IsZombie() )
    {
        cout << "file not found: " << fData->GetName() << endl;
        exit( EXIT_FAILURE );
    }

    VFITSIRFs a;
    a.open_fits_file( fFitsFileName );

    a.write_fits_header();

    a.write_histo2D( (TH2F*)fData->Get( "EffectiveAreaEtrueNoTheta2cut_offaxis" ),
                     "EFFECTIVE AREA",
                     (char *)"EFFAREA",
                     (char *)"m**2" );

    a.write_psf_gauss( (TH2F*)fData->Get( "AngResEtrue_offaxis" ) );

    a.write_fits_file();

    return 0;
}
