/*! \file combineEffectiveAreas.cpp
    \brief combine effective areas calculated for a one combination of ze,az,woff... into a single large file


*/

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

#include <iostream>
#include <stdlib.h>
#include <string>

#include <VGammaHadronCuts.h>
#include <VGlobalRunParameter.h>
#include <VInstrumentResponseFunctionRunParameter.h>
#include <VEnergyThreshold.h>

using namespace std;

/*
 *
 *  combine several effective area files into one
 *
 */
void merge( string ifile, char* outputfile, bool bFull = false , bool bMergeLogs = true )
{
    char hname[2000];
    
    // chain with merged effective area values
    TChain f( "fEffArea" );
    if( ifile.find( ".root" ) != string::npos )
    {
        sprintf( hname, "%s", ifile.c_str() );
    }
    else
    {
        sprintf( hname, "%s.root", ifile.c_str() );
    }
    int i_nMerged = f.Add( hname );
    if( i_nMerged == 0 )
    {
        cout << "error: no files found to merge: " << endl;
        cout << "\t" << hname << endl;
        cout << "exiting.." << endl;
        exit( EXIT_FAILURE );
    }
    sprintf( hname, "%s.root", outputfile );
    cout << "merging " << i_nMerged << " files to " << hname << endl;
    
    // activate branches to be included in merged files
    // these are the branches needed for the anasum analysis
    if( !bFull )
    {
        f.SetBranchStatus( "*", 0 );
        f.SetBranchStatus( "ze", 1 );
        f.SetBranchStatus( "az", 1 );
        f.SetBranchStatus( "azMin", 1 );
        f.SetBranchStatus( "azMax", 1 );
        f.SetBranchStatus( "Woff", 1 );
        f.SetBranchStatus( "noise", 1 );
        f.SetBranchStatus( "pedvar", 1 );
        f.SetBranchStatus( "index", 1 );
        f.SetBranchStatus( "nbins", 1 );
        f.SetBranchStatus( "e0", 1 );
        f.SetBranchStatus( "eff", 1 );
		f.SetBranchStatus( "effNoTh2", 1 );
        f.SetBranchStatus( "esys_rel", 1 );
        f.SetBranchStatus( "Rec_eff", 1 );
        f.SetBranchStatus( "Rec_eff_error", 1 );
        f.SetBranchStatus( "Rec_angRes_p68", 1 );
        f.SetBranchStatus( "Rec_angRes_p80", 1 );
        f.SetBranchStatus( "Rec_angRes_kingSigma", 1 );
        f.SetBranchStatus( "Rec_angRes_kingGamma", 1 );
        // full histograms - mostly needed for DL3 step
        // f.SetBranchStatus( "hResponseMatrixFine", 1 );
        // f.SetBranchStatus( "hResponseMatrixFineNoDirectionCuts", 1 );
        f.SetBranchStatus( "hEsysMCRelative", 1 );
        f.SetBranchStatus( "hEsysMCRelative2D", 1 );
        f.SetBranchStatus( "hEsysMCRelative2DNoDirectionCut", 1 );
        f.SetBranchStatus( "hAngularLogDiffEmc_2D", 1 );
    }
    f.Merge( hname );
    cout << "done.." << endl;
    
    // get one example of hEmc
    // (this is needed later to get the binning right)
    TFile* fO = new TFile( hname, "UPDATE" );
    if( fO->IsZombie() )
    {
        cout << "error writing hEmc to output file" << endl;
        return;
    }
    TH1D* hEmc = 0;
    if( f.GetFile() )
    {
        cout << "reading hEmc from " << f.GetFile()->GetName() << endl;
        hEmc = ( TH1D* )f.GetFile()->Get( "hEmc" );
    }
    if( !hEmc && f.GetBranchStatus( "hEmc" ) )
    {
        f.SetBranchAddress( "hEmc", &hEmc );
        f.GetEntry( 0 );
    }
    fO->cd();
    if( hEmc )
    {
        hEmc->Reset();
        hEmc->Write();
    }
    else
    {
        cout << "Error: can not find required MC histogram in ";
        cout << fO->GetName();
        cout << endl;
    }
    // get one example of IRF-runparameters for later checks in the analysis
    // (this assumes they are the same in all merged files!)
    TFile* ifirst = f.GetFile();
    if( !ifirst )
    {
        cout << "error finding pointer to first file in chain" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    VInstrumentResponseFunctionRunParameter* iRunPara = ( VInstrumentResponseFunctionRunParameter* )ifirst->Get( "makeEffectiveArea_runparameter" );
    if( !iRunPara )
    {
        cout << "error copying VInstrumentResponseFunctionRunParameter to output file" << endl;
        cout << "could not find them in file: " << ifirst->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    iRunPara->Write();
    // get one example of the gamma-hadron cuts
    // (this assume they are the same in all merged files!)
    VGammaHadronCuts* iCuts = ( VGammaHadronCuts* )ifirst->Get( "GammaHadronCuts" );
    if( iCuts )
    {
        cout << "copying gamma/hadron cuts from first file (" << ifirst->GetName() << ") into the output file" << endl;
        iCuts->Write();
    }
    else
    {
        cout << "error copying gamma/hadron cuts into output file" << endl;
        cout << "could not find them in file: " << ifirst->GetName() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fO->Close();
    
    // merge all log files
    if( bMergeLogs )
    {
        if( ifile.find( ".root" ) != string::npos )
        {
            sprintf( hname, "cat %s*.log > %s.combine.log", ifile.substr( 0, ifile.size() - 5 ).c_str(), outputfile );
        }
        else
        {
            sprintf( hname, "cat %s*.log > %s.combine.log", ifile.c_str(), outputfile );
        }
        cout << "merge log files into " << hname << endl;
        if( system( hname ) != 0 )
	      {
	         cout << "error merging log files" << endl;
        }
    }
    else
    {
        cout << "due to command line argument, we are NOT merging log files..." << endl;
    }
    
    cout << "done..";
}


int main( int argc, char* argv[] )
{

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
    if( argc < 4 )
    {
        cout << endl;
        cout << "combineEffectiveAreas <effective area files> <combined file> <write all histograms (default value is=false)> <merge log files (default value is=true>" << endl;
        cout << endl;
        cout << "   <effective area files>    without .root suffix (e.g. effArea*. Note need of \"...\")" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    cout << endl;
    cout << "combineEffectiveAreas (" << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "------------------------------------" << endl;
    cout << endl;
    
    bool mergelogs = true ;
    if( argc == 5 )
    {
        if( strcmp( argv[4], "true" ) == 0 )
        {
            mergelogs = true  ;
        }
        if( strcmp( argv[4], "false" ) == 0 )
        {
            mergelogs = false ;
        }
    }
    
    merge( argv[1], argv[2], ( bool )atoi( argv[3] ), mergelogs );
    
    cout << endl << endl;
    cout << "new combined effective area file: " << endl;
    cout << argv[2] << endl;

    // calculate energy thresholds
    // - will be written to effective area file
    string iEnergyTresholdFile = string(argv[2]) + ".root";
    VEnergyThreshold iEnergyThreshold( iEnergyTresholdFile, "UPDATE" );
    iEnergyThreshold.openEffectiveAreaFile( string(argv[2]) + ".root" );
    iEnergyThreshold.calculateEnergyThreshold();
    iEnergyThreshold.writeResults();
    
    // (all done)
    cout << endl << "end combineEffectiveAreas" << endl;
    
}

