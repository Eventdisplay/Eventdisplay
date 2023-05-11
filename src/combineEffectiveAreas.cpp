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

using namespace std;

/*
 *
 *  combine several effective area files into one
 *
 */
void merge( vector< string > file_list,
            string outputfile,
            bool bFull = false ,
            bool bMergeLogs = true )
{
    if( file_list.size() == 0 )
    {
        cout << "error: no files found to merge" << endl;
        cout << "exiting.." << endl;
        exit( EXIT_FAILURE );
    }
    TChain f( "fEffArea" );
    for( unsigned int i = 0; i < file_list.size(); i++ )
    {
        f.Add( file_list[i].c_str() );
    }
    if( outputfile.find( ".root" ) == string::npos )
    {
        outputfile += ".root";
    }
    cout << "merging " << file_list.size() << " files to " << outputfile << endl;
    
    // set branches to be included in merged files
    if( !bFull )
    {
        f.SetBranchStatus( "*", 0 );
        f.SetBranchStatus( "ze", 1 );
        f.SetBranchStatus( "az", 1 );
        f.SetBranchStatus( "azMin", 1 );
        f.SetBranchStatus( "azMax", 1 );
        f.SetBranchStatus( "Xoff", 1 );
        f.SetBranchStatus( "Yoff", 1 );
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
    f.Merge( outputfile.c_str() );
    cout << "done.." << endl;
    
    // get one example of hEmc
    // (this is needed later to get the binning right)
    TFile* fO = new TFile( outputfile.c_str(), "UPDATE" );
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
    VInstrumentResponseFunctionRunParameter* iRunPara =
        ( VInstrumentResponseFunctionRunParameter* )ifirst->Get( "makeEffectiveArea_runparameter" );
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
        cout << "copying gamma/hadron cuts from first file (";
        cout << ifirst->GetName() << ") into the output file" << endl;
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
}

void write_log_files( vector< string > file_list, string outputfile )
{
    // merge all log files
    ostringstream i_sys;
    for( unsigned int i = 0; i < file_list.size(); i++ )
    {
        if( file_list[i].find( ".root" ) != string::npos )
        {
            i_sys << "cat " << file_list[i].substr( 0, file_list[i].size() - 5 ).c_str() << ".log > ";
        }
        else
        {
            i_sys << "cat " << file_list[i] << ".log > ";
        }
    }
    
    i_sys << outputfile << ".combine.log";
    cout << "merge log files into " << i_sys.str() << endl;
    system( i_sys.str().c_str() );
    cout << "done.." << endl;
}

/*
 * return list of effective area files
 * to be merged
 *
 */
vector< string > readListOfFiles( string iFile )
{
    vector< string > iList;
    
    ifstream is;
    is.open( iFile.c_str() );
    if( !is )
    {
        cout << "error while reading file list " << iFile << endl;
        cout << "exiting...." << endl;
        exit( EXIT_FAILURE );
    }
    string is_line;
    
    while( getline( is, is_line ) )
    {
        iList.push_back( is_line );
    }
    
    is.close();
    
    return iList;
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
            exit( EXIT_FAILURE );
        }
    }
    if( argc < 4 )
    {
        cout << endl;
        cout << "combineEffectiveAreas <effective area file list> <combined file> <tree type>" << endl;
        cout << endl;
        cout << "  <effective area file list>  list of effective files to be merged" << endl;
        cout << "  <tree type>  effective area tree type (defines size of combined tree)" << endl;
        cout << "                - DL3 (default): entries required for DL3 analyis (large)" << endl;
        cout << "                - all          : all entries of original trees (largest)" << endl;
        cout << "                - anasum       : entries required for anasum analysis only (smallest)" << endl;
        cout << "                - DL3reduced   : histograms are written as regular arrays for DL3 analysis" << endl;
        cout << endl;
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
    vector< string > file_list = readListOfFiles( argv[1] );
    
    merge( file_list, argv[2], ( bool )atoi( argv[3] ), mergelogs );
    // write_log_files( file_list, argv[2] );
    
    cout << endl << "end combineEffectiveAreas" << endl;
}
