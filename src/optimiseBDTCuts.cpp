/*  optimiseBDTCuts
 *
 *  optimise BDT cuts
 *
 *  (used for VERITAS analysis)
 *
 *  this is a command line tool equivalent to macros/VTS/optimizeBDTcuts.C
 */


#include <iostream>
#include <string>
#include <vector>

#include "VTMVAEvaluator.h"

using namespace std;

void optimizeBDTcuts( string particleraterootfile, 
                      string weightFileDir, string weightFileName = "mva",
                      string MVAName = "BDT", unsigned int MVACounter = 0,
                      double observing_time_h = 20., 
                      int weightFileIndex_Emin = 0, int weightFileIndex_Emax = 4, 
                      int weightFileIndex_Zmin = 0, int weightFileIndex_Zmax = 2,
                      double significance = 3., int min_events = 5,
                      bool iPlotEfficiencyPlots = false, bool iPlotOptimisationResults = true,
                      string iWriteTMACuts = "" )
{
    // fixed parameters
    double min_sourcestrength_CU = 0.00001;
    double timeparticlerate = 3600.;
    double energyStepSize = 0.2;
    double min_backgroundrateratio = 1. / 5.;
    double min_backgroundevents = 0.;
    double signalefficiency = 0.90;
    // use energy bins defined for training
    // energyStepSize = -1.;

    
    VTMVAEvaluator a;
    a.setTMVAMethod( MVAName, MVACounter );
    a.setPrintPlotting( true );
    a.setPlotEfficiencyPlotsPerBin( iPlotEfficiencyPlots );
    a.setParticleNumberFile( particleraterootfile, timeparticlerate );

    a.setSensitivityOptimizationParameters( significance, min_events, observing_time_h, min_backgroundrateratio, min_backgroundevents );
    a.setSensitivityOptimizationFixedSignalEfficiency( signalefficiency );
    a.setSensitivityOptimizationSourceStrength( min_sourcestrength_CU );


    ostringstream iFullWeightFileName;
    iFullWeightFileName << weightFileDir << "/" << weightFileName;
    a.initializeWeightFiles( iFullWeightFileName.str().c_str(), 
                             weightFileIndex_Emin, weightFileIndex_Emax, 
                             weightFileIndex_Zmin, weightFileIndex_Zmax, 
                             // energyStepSize, "VTS", "UseInterpolatedCounts" );
                             energyStepSize, "VTS", "UseAveragedCounts" );

    if( iPlotOptimisationResults )
    {
        a.plotSignalAndBackgroundEfficiencies( false, 1.e-2, -0.2, 0.4 );
    }
    
    a.printOptimizedMVACutValues();
    if( iWriteTMACuts.size() > 0 )
    {
        a.writeOptimizedMVACutValues( iWriteTMACuts );
    }

}


int main( int argc, char* argv[] )
{
    cout << endl;
    cout << "optimiseBDTCuts (version ";
    cout << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "=================================================" << endl;
    cout << endl ;
    if( argc < 5 )
    {
        cout << "writeParticleRateFilesForTMVA" << endl;
        cout << "\t\t <particle rate file>" << endl;
        cout << "\t\t <BDT weight file directory> " << endl;
        cout << "\t\t <BDT weight file name> " << endl;
        cout << "\t\t <output file with optimised cuts> " << endl;
        cout << "\t\t [observing time (h)] " << endl;
        exit( EXIT_SUCCESS );
    }

    unsigned int MVACounter = 0;
    string MVAName = "BDT";
    int weightFileIndex_Emin = 0;
    int weightFileIndex_Emax = 4;
    int weightFileIndex_Zmin = 0;
    int weightFileIndex_Zmax = 2;
    double observing_time_h = 20.;
    if( argc == 6 )
    {
         observing_time_h = atof( argv[5] );
    } 
    double significance = 3.;
    int min_events = 5;

    optimizeBDTcuts( argv[1], argv[2], argv[3],
            MVAName, MVACounter,
            observing_time_h,
            weightFileIndex_Emin, weightFileIndex_Emax,
            weightFileIndex_Zmin, weightFileIndex_Zmax,
            significance, min_events,
            false, false, argv[4] );

}          
