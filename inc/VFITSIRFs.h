/* convert IRFS from ROOT to FITS
 * (includes sensitivity IRFS)
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <fitsio.h>

#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"

using namespace std;

class VFITSIRFs
{
    private:
    fitsfile* fptr;

    vector< vector< float > > get_baseline_axes( TH2F *h );
    bool printerror( int status );
    bool write_table( vector< vector< float > > table );


    public:
    VFITSIRFs();
   ~VFITSIRFs() {}

    bool open_fits_file( string fits_file_name );
    bool write_fits_header();
    bool write_histo2D( TH2F *h, 
                        string name,
                        char* col_name,
                        char* col_unit );
    bool write_psf_gauss( TH2F *h );
    bool write_fits_file();
};
