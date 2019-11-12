/*! \file run next day analysis and produce text and FITS output

*/

#include "VAnalysisUtilities.h"
#include "VGlobalRunParameter.h"
#include "VFITS.h"
#include "VFluxCalculation.h"
#include "VRunList.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void help()
{
    cout << endl;
    cout << "VTS.next_day (version " << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "===========================" << endl;
    cout << endl;
    cout << "VTS.next_day <anasum input file> <output file (.dat and .fits)> <target name>";
    cout << "[merge information into one FITS-file=0/1 (default=1)] [debug=0/1 (default=0)]" << endl;
    cout << endl;
    exit( 0 );
}


int main( int argc, char* argv[] )
{
    bool fDebug = false;
    bool fMergeFITSFiles = true;
    
    /////////////////////////////////////////
    // read command line parameters
    if( argc == 1 )
    {
        help();
    }
    
    string fDataFile = argv[1];
    string fOUTFile = argv[2];
    string fTargetName = argv[3];
    if( argc >= 5 && atof( argv[4] ) == 0 )
    {
        fMergeFITSFiles = false;
    }
    if( argc == 6 && atof( argv[5] ) == 1 )
    {
        fDebug = true;
    }
    
    /////////////////////////////////////////
    // calculate total fluxes and upper limits
    double fMinEnergy = 0.2;
    double fGamma = 2.5;
    char ifile[800];
    sprintf( ifile, "%s", fDataFile.c_str() );
    
    //calculate fluxes for all runs
    VFluxCalculation* flux = new VFluxCalculation( ifile );
    if( flux->IsZombie() )
    {
        return 1;
    }
    flux->setDebug( fDebug );
    flux->calculateIntegralFlux( fMinEnergy );
    if( fDebug )
    {
        flux->printResults();
    }
    
    // read run list
    VAnalysisUtilities a;
    a.openFile( fDataFile, -1, true, fDebug );
    if( !a.IsZombie() )
    {
        CRunSummary* c = a.getRunSummaryTree( -1 );
        
        // open output stream
        ofstream fResults;
        fResults.open( ( fOUTFile + ".dat" ).c_str() );
        if( fResults && c )
        {
            double tOn_total = 0;
            for( int i = 0; i < c->fChain->GetEntries(); i++ )
            {
                c->GetEntry( i );
                if( i == 0 )
                {
                
                    fResults << setw( 12 ) << left << "MJD" ;
                    fResults << setw( 8 ) << left << "RA" ;
                    fResults << setw( 8 ) << left << "Dec" ;
                    fResults << setw( 9 ) << left << "runOn" ;
                    fResults << setw( 12 ) << left << "NOn" ;
                    fResults << setw( 12 ) << left << "NOff" ;
                    fResults << setw( 18 ) << left << "tOn w/o deadtime";
                    fResults << setw( 12 ) << left << "Signi" ;
                    fResults << setw( 8 ) << left << "OffNorm" ;
                    fResults << setw( 15 ) << left << "Rate" ;
                    fResults << setw( 15 ) << left << "RateE" ;
                    fResults << setw( 12 ) << left << "Flux" ;
                    fResults << setw( 12 ) << left << "FluxE";
                    fResults << setw( 14 ) << left << "FluxVsCrab";
                    fResults << setw( 12 ) << left << "FluxUL";
                    fResults << setw( 14 ) << left << "FluxULVsCrab";
                    fResults << endl;
                    fResults << setw( 12 ) << left << " " ;
                    fResults << setw( 8 ) << left << "[deg]" ;
                    fResults << setw( 8 ) << left << "[deg]" ;
                    fResults << setw( 9 ) << left << " " ;
                    fResults << setw( 12 ) << left << " " ;
                    fResults << setw( 12 ) << left << " " ;
                    fResults << setw( 18 ) << left << "[min]";
                    fResults << setw( 12 ) << left << "[sigma]" ;
                    fResults << setw( 8 ) << left << " " ;
                    fResults << setw( 15 ) << left << "[gammas/min]" ;
                    fResults << setw( 15 ) << left << "[gammas/min]" ;
                    fResults << setw( 12 ) << left << "[cm^-2 s^-1]" ;
                    fResults << setw( 12 ) << left << "[cm^-2 s^-1]";
                    fResults << setw( 14 ) << left << "[percent]";
                    fResults << setw( 12 ) << left << "[cm^-2 s^-1]";
                    fResults << setw( 14 ) << left << "[percent]";
                    fResults << endl;
                }
                const VFluxDataPoint* iFluxDataPoint = flux->getFluxDataPerRun( c->runOn );
                fResults.precision( 9 );
                if( c->runOn > 0 )
                {
                    fResults << setw( 12 ) << left << c->MJDOn;
                }
                else
                {
                    fResults << setw( 12 ) << left << "Total:";
                }
                fResults.precision( 4 );
                fResults << setw( 8 ) << left << c->TargetRAJ2000 ;
                fResults << setw( 8 ) << left << c->TargetDecJ2000;
                fResults << setw( 9 ) << left << c->runOn ;
                fResults << setw( 12 ) << left << c->NOn ;
                fResults << setw( 12 ) << left << c->NOff ;
                if( c->runOn > 0 )
                {
                    fResults << setw( 18 ) << left << c->tOn / 60.*( 1. - c->DeadTimeFracOn ) ;
                }
                else
                {
                    fResults << setw( 18 ) << left << tOn_total ;
                }
                fResults << setw( 12 ) << left << c->Signi;
                fResults << setw( 8 ) << left << c->OffNorm ;
                fResults << setw( 15 ) << left << c->Rate ;
                fResults << setw( 15 ) << left << c->RateE ;
                fResults << setw( 12 ) << left << iFluxDataPoint->fFlux;
                fResults << setw( 12 ) << left << iFluxDataPoint->fFluxCI_1sigma;
                fResults << setw( 14 ) << left << VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy, iFluxDataPoint->fFlux, fGamma );
                fResults << setw( 12 ) << left << iFluxDataPoint->fFluxUL;
                fResults << setw( 14 ) << left << VFluxAndLightCurveUtilities::convertPhotonFlux_to_CrabUnits( fMinEnergy, iFluxDataPoint->fFluxUL, fGamma );
                fResults << endl;
                tOn_total += c->tOn / 60.0 * ( 1.0 - c->DeadTimeFracOn ) ;
            }
        }
        
        /////////////////////////////////////////
        // convert to fits
        
        cout << "---------Start production of FITS ouput file ------" << endl;
        VFITS f( fDataFile, fOUTFile + ".fits" , fTargetName, fMergeFITSFiles, fDebug );
        f.writeCumSignificance( fDebug );
        f.writeSignificanceDistribution( fDebug );
        f.writeLightCurve( fDebug );
        f.writeThetaSquareDistribution( fDebug );
        f.writeSignificanceSkyMap( fDebug );
        f.writeExcessSkyMap( fDebug );
        f.writeNightlyFlux( fDebug, ( fOUTFile + ".flux" ).c_str() );
        f.writeAlphaSkyMap( fDebug );
        f.writeRawCountsSkyMap( fDebug );
        //f.writeEnergySpectrum( fDebug ); //Uncomment this to include the energy spectrum. Obs that this function can fail quietly, in which case the energy spectrum won't get written.
        f.writeFITSFile( fDebug );
        
    }
}
