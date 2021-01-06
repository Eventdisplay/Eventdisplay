/* \file writeCTAWPPhysSensitivityFiles write CTA WP Phys sensitivity files

   write a simple root file with histograms for sensitivities, effective areas,
   angular and energy resolution, migration matrix, etc.

   follow here CTA WP-Phys standards (with minor modifications)


*/

#include "VWPPhysSensitivityFile.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main( int argc, char* argv[] )
{
    /////////////////////
    // input parameters
    if( argc != 8 && argc != 7 && argc != 9 )
    {
        cout << endl;
        cout << "./writeCTAWPPhysSensitivityFiles <sub array> <observing time> <data directory>";
        cout << " <outputfile> <observatory (CTA/V5/V6)> <calculate offset sensitivites=0/1>";
        cout << "[recid (default=0)] [off-axis fine binning (default=0)]" << endl;
        cout << endl;
        cout << "\t observing time: give unit without space, e.g. 50h, 10m, 2s" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string fSubArray = argv[1];
    // observing time (translate from whatever unit into h)
    double fObservingTime_h = 0.;
    string iObstime = argv[2];
    if( iObstime.find( "h" ) != string::npos )
    {
        fObservingTime_h = atof( iObstime.substr( 0, iObstime.find( "h" ) ).c_str() );
    }
    else if( iObstime.find( "m" ) != string::npos )
    {
        fObservingTime_h = atof( iObstime.substr( 0, iObstime.find( "m" ) ).c_str() ) / 60.;
    }
    else if( iObstime.find( "s" ) != string::npos )
    {
        fObservingTime_h = atof( iObstime.substr( 0, iObstime.find( "s" ) ).c_str() ) / 3600.;
    }
    else
    {
        fObservingTime_h = atof( iObstime.c_str() );
    }
    string fDataDirectory = argv[3];
    string fOutputFile = argv[4];
    string fObservatory = argv[5];
    bool   fWriteOffsetFiles = atoi( argv[6] );
    int    fReconstructionID = 0;
    if( argc >= 8 )
    {
        fReconstructionID = atoi( argv[7] );
    }
    bool bOffAxisFineBinning = false;
    if( argc == 9 )
    {
        bOffAxisFineBinning = ( bool )( atoi( argv[8] ) );
    }
    
    /////////////////////
    // initialization
    VWPPhysSensitivityFile* iData = new VWPPhysSensitivityFile();
    iData->setDebug( true );
    iData->setObservatory( fObservatory );
    // Note: offset binning hardwired
    vector< double > iWobbleMin;
    vector< double > iWobbleMax;
    if( fWriteOffsetFiles )
    {
        // prod2 offsets
        /* iWobbleMin.push_back( 0.0 );
        iWobbleMax.push_back( 1.0 );
        iWobbleMin.push_back( 1.0 );
        iWobbleMax.push_back( 2.0 );
        iWobbleMin.push_back( 2.0 );
        iWobbleMax.push_back( 3.0 );
        iWobbleMin.push_back( 3.0 );
        iWobbleMax.push_back( 3.5 );
        iWobbleMin.push_back( 3.5 );
        iWobbleMax.push_back( 4.0 );
        iWobbleMin.push_back( 4.0 );
        iWobbleMax.push_back( 4.5 );
        iWobbleMin.push_back( 4.5 ); */
        // prod3 offsets
        if( !bOffAxisFineBinning )
        {
            iWobbleMin.push_back( 0.0 );
            iWobbleMax.push_back( 1.0 );
            iWobbleMin.push_back( 1.0 );
            iWobbleMax.push_back( 2.0 );
            iWobbleMin.push_back( 2.0 );
            iWobbleMax.push_back( 3.0 );
            iWobbleMin.push_back( 3.0 );
            iWobbleMax.push_back( 4.0 );
            iWobbleMin.push_back( 4.0 );
            iWobbleMax.push_back( 5.0 );
            iWobbleMin.push_back( 5.0 );
            iWobbleMax.push_back( 6.0 );
        }
        // prod3 offsets (fine binning)
        // (fine binning is achieved in analysis by overlapping
        //  bins in off-axis angle. The bins are given also below
        //  (commented out).
        else
        {
            iWobbleMin.push_back( 0.0 );
            iWobbleMax.push_back( 1.0 );
            iWobbleMin.push_back( 1.0 );
            iWobbleMax.push_back( 2.0 );
            iWobbleMin.push_back( 2.0 );
            iWobbleMax.push_back( 2.4 );
            iWobbleMin.push_back( 2.4 );
            iWobbleMax.push_back( 2.65 );
            iWobbleMin.push_back( 2.65 );
            iWobbleMax.push_back( 2.9 );
            iWobbleMin.push_back( 2.9 );
            iWobbleMax.push_back( 3.15 );
            iWobbleMin.push_back( 3.15 );
            iWobbleMax.push_back( 3.4 );
            iWobbleMin.push_back( 3.4 );
            iWobbleMax.push_back( 3.65 );
            iWobbleMin.push_back( 3.65 );
            iWobbleMax.push_back( 3.9 );
            iWobbleMin.push_back( 3.9 );
            iWobbleMax.push_back( 4.15 );
            iWobbleMin.push_back( 4.15 );
            iWobbleMax.push_back( 4.4 );
            iWobbleMin.push_back( 4.4 );
            iWobbleMax.push_back( 4.65 );
            iWobbleMin.push_back( 4.65 );
            iWobbleMax.push_back( 4.9 );
            iWobbleMin.push_back( 4.9 );
            iWobbleMax.push_back( 5.15 );
            iWobbleMin.push_back( 5.15 );
            iWobbleMax.push_back( 5.4 );
            iWobbleMin.push_back( 5.4 );
            iWobbleMax.push_back( 6.0 );
        }
    }
    // sub array
    iData->setDataFiles( fSubArray, fReconstructionID );
    // observing time
    iData->setObservationTime( fObservingTime_h );
    // output file
    if( !iData->initializeOutputFile( fOutputFile ) )
    {
        exit( EXIT_FAILURE );
    }
    // Crab spectrum from HEGRA (CTA default)
    iData->setCrabSpectrum( "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat", 5 );
    // CR spectra (protons + electrons)
    iData->setCosmicRaySpectrum( "$CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CR.dat", 0, 8 );
    
    ////////////////////////////////////
    // energy ranges and binning definitions
    double fEnergy_min = -1.9;
    double fEnergy_max =  2.3; // 200 TeV
    // standard binning with 5 logarithmic bins per decade
    int fNbins_diffBining = ( int )( ( fEnergy_max - fEnergy_min ) / 0.2 + 0.5 );
    int fNbins_fineBinning = (fNbins_diffBining-1) * 25;   // used e.g. for migration matrix
    
    cout << "IRF binning:" << endl;
    cout << "\t log_10 E [" << fEnergy_min << ", " << fEnergy_max << "]" << endl;
    cout << "\t binning: " << fNbins_diffBining << ", fine binning: " << fNbins_fineBinning << endl;
    
    /////////////////////////////////////
    // on source histograms
    // initialize histogram with the standard binning used in the CTA WP Phys group
    iData->initializeHistograms( fNbins_diffBining, fEnergy_min, fEnergy_max,
                                 fNbins_fineBinning, 9999 );
    if( !iData->fillHistograms1D( fDataDirectory, true ) )
    {
        cout << "error filling on source histograms" << endl;
        exit( EXIT_FAILURE );
    }
    
    /////////////////////////////////////
    // off source histograms
    for( unsigned int i = 0; i < iWobbleMin.size(); i++ )
    {
        cout << "Camera offset: " << i << "\t" << iWobbleMin[i] << "\t" << iWobbleMax[i] << endl;
        // initialize histogram with the standard binning used in the CTA WP Phys group
        iData->initializeHistograms( fNbins_diffBining, fEnergy_min, fEnergy_max,
                                     fNbins_fineBinning, i );
        if( !iData->fillHistograms1D( fDataDirectory, false ) )
        {
            // allow last bin to fail
            // ( might not have enough statistics for the
            //   very large offset)
            if( iWobbleMin.size() > 0 && i == iWobbleMin.size() - 1 )
            {
                cout << "Warning: wobble bin " << i << " not filled" << endl;
            }
            // fail for all other bins
            else
            {
                cout << "error filling off source histogram " << i << endl;
                exit( EXIT_FAILURE );
            }
        }
    }
    // copy all 1D off-axis histograms over to 2D/3D histograms
    if( iWobbleMin.size() > 0 )
    {
        iData->fillHistograms2D( iWobbleMin, iWobbleMax );
    }
    
    iData->terminate();
    cout << "end of calculating sensitivities" << endl;
}
