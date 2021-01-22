/* writeParticleRateFilesFromEffectiveAreas

   write files with events rates for on (gamma) and off (protons+electrons) simulations

   files are needed e.g. for setting the optimal cut value for TMVA cuts

   Input:
       * effective areas (gammas/protons/electrons/ec)
       * cosmic ray spectra (read from files)

   Output:
       * root file with signal and background rates
       * [optional] angular resolution histograms

*/

#include "VHistogramUtilities.h"
#include "VInstrumentResponseFunctionReader.h"
#include "VWPPhysSensitivityFile.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

/*

     calculate and write angular angular resolution vs containment vs energy

     this allows to easily read e.g. the containment radius for 68% or 80% from the histogram

*/
void writeAngResHistogram( char* iMC_Gamma = 0, string iParticleNumberFile = "particleNumbers.tmp.root" )
{
    if( !iMC_Gamma )
    {
        return;
    }
    
    VInstrumentResponseFunctionReader iR;
    iR.fillData( iMC_Gamma );
    for( unsigned int i = 0; i < iR.fIRF_TreeNames.size(); i++ )
    {
        // get angular resolution tree (does not matter of 68 or 80% tree)
        if( iR.fIRF_TreeNames[i] == "t_angular_resolution_080p" )
        {
            TH2D* h = 0;
            if( i < iR.fIRF_Data.size() && iR.fIRF_Data[i] )
            {
                h = iR.fIRF_Data[i]->f2DHisto[VInstrumentResponseFunctionData::E_DIFF];
            }
            TH2D* hRes = VHistogramUtilities::calculateContainmentDistance( h, "AngResCumulative" );
            if( hRes )
            {
                TFile hh( iParticleNumberFile.c_str(), "update" );
                cout << "writing angular resolution histogram to " << hh.GetName() << endl;
                hRes->Write();
                hh.Close();
            }
            break;
        }
    }
}


/*
 *
 * calculate and write events rates for different particle types
 *
 * one should use the same CR spectra as in the other analysis
 *
 */
void writeParticleNumberFile( char* iMC_Gamma = 0, char* iMC_Proton = 0, char* iMC_Electron = 0,
                              unsigned int iCrabSpec_ID = 5, string iParticleNumberFile = "particleNumbers.tmp.root",
                              string iObservatory = "CTA" )
{
    string iESpecDataFile_CrabNebula = "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat";
    string iESpecDataFile_CosmicRays = "$" + iObservatory + "_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CR.dat";
    
    if( iMC_Gamma && iMC_Proton )
    {
        // use sensitivity calculator for differential flux calculation
        VSensitivityCalculator b;
        b.setDebug( false );             // creates lots of output (if set to true)
        // set Crab Nebula spectrum
        b.setEnergySpectrumfromLiterature( iESpecDataFile_CrabNebula, iCrabSpec_ID );
        // set output result file name
        b.setWriteParticleNumberFile( iParticleNumberFile );
        // CTA
        int i_Azbin_gamma = 0;
        double i_index_gamma = 2.5;
        int i_noise_gamma = 250;
        double i_woff_gamma = 0.;
        
        int i_Azbin_proton = 0;
        double i_index_proton = 2.5;
        int i_noise_proton = 250;
        double i_woff_proton = 0.;
        cout << "SETTING EFFECTIVE AREA SEARCH VALUES TO " << iObservatory << endl;
        //////////////////////////////////////////////////////////////////////////
        
        // gammas
        b.setMonteCarloParameters( 1, iESpecDataFile_CrabNebula, iCrabSpec_ID, iMC_Gamma, 20.,
                                   i_Azbin_gamma, i_woff_gamma, i_noise_gamma, i_index_gamma,
                                   -1.9 );
        // protons
        b.setMonteCarloParameters( 14, iESpecDataFile_CosmicRays, 0, iMC_Proton, 20.,
                                   i_Azbin_proton, i_woff_proton, i_noise_proton, i_index_proton,
                                   -1.9 );
        // electrons (spectral index?)
        if( iMC_Electron )
        {
            b.setMonteCarloParameters( 2, iESpecDataFile_CosmicRays, 8, iMC_Electron, 20., 0, 0.0, 250, 3.0, -1.9 );
        }
        // calculate differential fluxes for 5 bins per decade (0.2)
        b.calculateParticleNumberGraphs_MC( 0.2 );
    }
}

/*******************************************************************************************************

    write files with particle number spectra for on (gamma) and off (protons+electrons) counts

    files are needed e.g. for setting the optimal cut value for TMVA cuts

    (for all sub arrays)

*/
int main( int argc, char* argv[] )
{
    if( argc != 6 && argc != 7 )
    {
        cout << endl;
        cout << "writeParticleRateFilesFromEffectiveAreas <sub array> <onSource/cone>";
        cout <<  "<reconstruction ID> <directory with effective areas> ";
        cout << "<output directory> ";
        cout << "[directory with angular resolution files]" << endl;
        cout << argc << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    
    cout << endl;
    cout << "writeParticleRateFilesFromEffectiveAreas" << endl;
    cout << "========================================" << endl;
    cout << endl;
    
    string SubArray = argv[1];
    string iOnSource = argv[2];
    int    iRecID = atoi( argv[3] );
    string iDataDir = argv[4];
    string iODataDir = argv[5];
    string iAngResDir = "";
    if( argc == 7 )
    {
        iAngResDir = argv[6];
    }
    
    // hardwired total number of off source bins
    int iOffSetCounter = -1;
    if( iOnSource == "onSource" )
    {
        iOffSetCounter = 0;
    }
    else if( iOnSource == "cone" )
    {
        iOffSetCounter = 6;
    }
    // fine binning in off axis
    else if( iOnSource == "coneFB" )
    {
        iOffSetCounter = 17;
    }
    else
    {
        cout << iOnSource << " invalid off source descriptor; should be: onSource or cone" << endl;
        exit( EXIT_FAILURE );
    }
    
    // effective area file names
    string iMC_Gamma_onSource = "gamma_onSource";
    string iMC_Gamma_cone = "gamma_cone";
    string iMC_Proton = "proton";
    string iMC_Proton_onSource = "proton_onSource";
    string iMC_Electron = "electron";
    string iMC_Electron_onSource = "electron_onSource";
    
    char iGamma[800];
    char iProton[800];
    char iElectron[800];
    char iParticleNumberFile[800];
    
    cout << "STARTING SUBARRAY " << SubArray << endl;
    
    // on-axis rates
    if( iMC_Gamma_onSource.size() > 0 )
    {
        sprintf( iParticleNumberFile, "%s/ParticleNumbers.%s.00.root", iODataDir.c_str(), SubArray.c_str() );
        sprintf( iGamma, "%s/%s.%s_ID%d.eff-%d.root", iDataDir.c_str(), iMC_Gamma_onSource.c_str(), SubArray.c_str(), iRecID, 0 );
        sprintf( iProton, "%s/%s.%s_ID%d.eff-%d.root", iDataDir.c_str(), iMC_Proton_onSource.c_str(), SubArray.c_str(), iRecID, 0 );
        if( iMC_Electron_onSource.size() > 0 )
        {
            sprintf( iElectron, "%s/%s.%s_ID%d.eff-%d.root", iDataDir.c_str(), iMC_Electron_onSource.c_str(), SubArray.c_str(), iRecID, 0 );
            writeParticleNumberFile( iGamma, iProton, iElectron, 5, iParticleNumberFile );
        }
        else
        {
            writeParticleNumberFile( iGamma, iProton, 0, 5, iParticleNumberFile );
        }
        // angular resolution histogram
        if( iAngResDir.size() > 0 )
        {
            sprintf( iGamma, "%s/%s.%s_ID%d.eff-%d.root", iAngResDir.c_str(), iMC_Gamma_onSource.c_str(), SubArray.c_str(), iRecID, 0 );
            writeAngResHistogram( iGamma, iParticleNumberFile );
        }
    }
    
    // off-axis rates
    for( int j = 0; j < iOffSetCounter; j++ ) // use first bin on source particle file
    {
    
        sprintf( iParticleNumberFile, "%s/ParticleNumbers.%s.%d.root", iODataDir.c_str(), SubArray.c_str(), j );
        sprintf( iGamma, "%s/%s.%s_ID%d.eff-%d.root", iDataDir.c_str(), iMC_Gamma_cone.c_str(), SubArray.c_str(), iRecID, j );
        sprintf( iProton, "%s/%s.%s_ID%d.eff-%d.root", iDataDir.c_str(), iMC_Proton.c_str(), SubArray.c_str(), iRecID, j );
        
        if( iMC_Electron.size() > 0 )
        {
            sprintf( iElectron, "%s/%s.%s_ID%d.eff-%d.root", iDataDir.c_str(), iMC_Electron.c_str(), SubArray.c_str(), iRecID, j );
            writeParticleNumberFile( iGamma, iProton, iElectron, 5, iParticleNumberFile );
        }
        else
        {
            writeParticleNumberFile( iGamma, iProton, 0, 5, iParticleNumberFile );
        }
        // angular resolution histogram
        if( iAngResDir.size() > 0 )
        {
            sprintf( iGamma, "%s/%s.%s_ID%d.eff-%d.root", iAngResDir.c_str(), iMC_Gamma_cone.c_str(), SubArray.c_str(), iRecID, j );
            writeAngResHistogram( iGamma, iParticleNumberFile );
        }
    }
}

