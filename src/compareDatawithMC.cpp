/*! \file compareDatawithMC.cpp
 *  \brief compare MC gamma-ray distributions with on-off distributions from data
 *         (e.g. Crab Nebula or Mrk 421 flare data)
 *
*/

#include "VGlobalRunParameter.h"
#include "VDataMCComparision.h"

#include "TFile.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct sInputData
{
    string fType;     // 0 = sims, 1 = on, 2 = off
    string fFileName;
    int    fNTelescopes;
    double fWobbleNorth;
    double fWobbleEast;
    bool   fWobbleFromDataTree;
    vector< double > fTelX;
    vector< double > fTelY;
    vector< double > fTelZ;
    double fAz_deg_min;
    double fAz_deg_max;
    double fZe_deg_min;
    double fZe_deg_max;
};

vector< sInputData > fInputData;

/*
 * read telescope position from telconfig tree
 *
 * (very inefficient, as all files are chained - but positions should be the same)
 *
 */
double getTelescopePositions( string iF, vector< double >& iX, vector< double >& iY, vector< double >& iZ, int iNTel )
{
    cout << "reading telescope positions for " << iF << endl;
    double r_max = 0.;
    
    TChain* c = new TChain( "telconfig" );
    if( !c->Add( iF.c_str() ) )
    {
        cout << "error: no tree with telescope positions (telconfig) found" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    if( c->GetEntries() < iNTel )
    {
        cout << "error: invalid number of telescopes: expected " << iNTel << ", found " << c->GetEntries() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    float x = 0.;
    float y = 0.;
    float z = 0.;
    c->SetBranchAddress( "TelX", &x );
    c->SetBranchAddress( "TelY", &y );
    c->SetBranchAddress( "TelZ", &z );
    for( int i = 0; i < iNTel; i++ )
    {
        c->GetEntry( i );
        
        iX.push_back( ( double )x );
        iY.push_back( ( double )y );
        iZ.push_back( ( double )z );
        cout << "\t telescope " << i + 1 << "\t" << x << "\t" << y << "\t" << z << endl;
        if( sqrt( x * x + y * y ) > r_max )
        {
            r_max = sqrt( x * x + y * y );
        }
    }
    cout << endl;
    
    return r_max;
}

/*
 * read input parameters, MC and data files from a configuration file
 *
 */
void readInputfile( string fInputFile )
{
    ifstream is;
    is.open( fInputFile.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error: input file list not found" << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }
    cout << "reading input file list " << fInputFile << endl;
    cout << endl;
    string is_line;
    string temp;
    
    sInputData a;
    
    while( getline( is, is_line ) )
    {
        if( is_line.size() > 0 )
        {
            istringstream is_stream( is_line );
            is_stream >> temp;
            if( temp != "*" )
            {
                continue;
            }
            
            // check that there are enough parameters in this line
            istringstream is_check( is_line );
            int z = 0;
            while( !( is_check >> std::ws ).eof() )
            {
                is_check >> temp;
                z++;
            }
            if( z != 10 )
            {
                cout << "error reading input file, not enough parameters in this line: " << endl << is_line << endl;
                cout << "require 10, found " << z << endl;
                cout << "...exiting" << endl;
                exit( EXIT_FAILURE );
            }
            
            is_stream >> a.fType;
            is_stream >> a.fFileName;
            is_stream >> temp;
            a.fNTelescopes = atoi( temp.c_str() );
            is_stream >> temp;
            a.fWobbleNorth = atof( temp.c_str() );
            is_stream >> temp;
            a.fWobbleEast = atof( temp.c_str() );
            if( a.fWobbleNorth < -98. || a.fWobbleEast < -98. )
            {
                a.fWobbleFromDataTree = true;
            }
            else
            {
                a.fWobbleFromDataTree = false;
            }
            
            // az range
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> a.fAz_deg_min;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> a.fAz_deg_max;
            }
            // ze range
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> a.fZe_deg_min;
            }
            if( !( is_stream >> std::ws ).eof() )
            {
                is_stream >> a.fZe_deg_max;
            }
            
            getTelescopePositions( a.fFileName, a.fTelX, a.fTelY, a.fTelZ, a.fNTelescopes );
            
            fInputData.push_back( a );
        }
    }
    is.close();
}


/**********************************************************************
 **********************************************************************
 **********************************************************************
*/
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
            exit( EXIT_SUCCESS );
        }
    }
    cout << endl;
    cout << "compareDatawithMC (" << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "==========================" << endl << endl;
    
    if( argc < 4 )
    {
        cout << "compare MC simulations with excess events from data runs " << endl;
        cout << "(e.g. from Crab Nebula or Mrk 421 observations)" << endl;
        cout << endl;
        cout << endl;
        cout << "compareDatawithMC <input file list> <cut> <outputfile> [BDT gamma/hadron cuts] [shower max zenith angle (default=20deg)]" << endl;
        cout << endl;
        cout << "\t input file list: see example file COMPAREMC.runparameter in the parameter files directory" << endl;
        cout << "\t cuts: " << endl;
        cout << "\t\t cut=-3:        theta2 cut only (RECOMMENDED CUT)" << endl;
        cout << "\t\t in most cases, the following cuts should not be used: " << endl;
        cout << "\t\t cut=-2:        no cuts" << endl;
        cout << "\t\t cut=-1:        stereo cuts (MSCW, etc.)" << endl;
        cout << "\t\t cut=1,2,...:   single telescope cuts on telescope 1,2" << endl;
        cout << endl;
        cout << "\t output file:     results file" << endl;
        cout << "\t                  (use VPlotCompareDataWithMC (in shared library)for plotting)" << endl;
        cout << endl;
        cout << "\t use BDT cuts for gamma/hadron separation: 0 = no (default), 1 = yes" << endl;
        cout << "\t cut file needs to be indicated within VDataMCComparision::initialGammaHadronCuts()" << endl;
        cout << endl;
        cout << "Note: most cuts are hardwired in VDataMCComparision::fillHistograms()" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string fInputFile = argv[1];
    
    // read input parameters and files from parameter file
    readInputfile( fInputFile );
    
    int fSingleTelescopeCuts = atoi( argv[2] );
    
    string fOutputfile = argv[3];
    
    bool fCalculateMVACut = false;
    
    if( argc > 4 && atoi( argv[4] ) == 1 )
    {
        fCalculateMVACut = true;
    }
    double fShowerMaxZe_deg = 20.;
    if( argc > 5 )
    {
        fShowerMaxZe_deg = atof( argv[5] );
    }
    
    // test number of telescopes
    int iNT = 0;
    for( unsigned int i = 0; i < fInputData.size(); i++ )
    {
        if( i == 0 )
        {
            iNT = fInputData[i].fNTelescopes;
        }
        else
        {
            if( fInputData[i].fNTelescopes != iNT )
            {
                cout << "error: number of telescopes differ, comparision not possible" << endl;
                cout << "...exiting" << endl;
                exit( EXIT_FAILURE );
            }
        }
    }
    
    // -------- end of reading input parameters
    
    TH1D* hAzOn = 0;
    // output file
    TFile* fout = new TFile( fOutputfile.c_str(), "RECREATE" );
    
    ////////////////////////////////////////////////
    // get AZ weighted histogram
    for( unsigned int i = 0; i < fInputData.size(); i++ )
    {
        if( fInputData[i].fType == "ON" )
        {
            VDataMCComparision iTemp( fInputData[i].fType, fInputData[i].fNTelescopes, fCalculateMVACut );
            iTemp.setAzRange( fInputData[i].fAz_deg_min, fInputData[i].fAz_deg_max );
            iTemp.setZeRange( fInputData[i].fZe_deg_min, fInputData[i].fZe_deg_max );
            iTemp.setShowerMaximZe_deg( fShowerMaxZe_deg );
            hAzOn = iTemp.getAzimuthWeightingHistogram( fInputData[i].fFileName );
            if( hAzOn )
            {
                hAzOn->SetDirectory( 0 );
                hAzOn->Write();
            }
        }
    }
    fout->Close();
    
    if( hAzOn )
    {
        cout  << "Number of entries / mean of az weighting histogram: ";
        cout <<  hAzOn->GetEntries() << ", " << hAzOn->GetMean() << endl;
    }
    
    
    ////////////////////////////////////////////////
    // now analyse the data
    vector< VDataMCComparision* > fStereoCompare;
    VDataMCComparision* fStereoCompareOn = 0;
    VDataMCComparision* fStereoCompareOff = 0;
    
    for( unsigned int i = 0; i < fInputData.size(); i++ )
    {
        cout << fInputData[i].fType << endl;
        cout << "----" << endl;
        fStereoCompare.push_back( new VDataMCComparision( fInputData[i].fType, fInputData[i].fNTelescopes, fCalculateMVACut ) );
        fStereoCompare.back()->setAzRange( fInputData[i].fAz_deg_min, fInputData[i].fAz_deg_max );
        fStereoCompare.back()->setZeRange( fInputData[i].fZe_deg_min, fInputData[i].fZe_deg_max );
        // get telescope coordinates
        fStereoCompare.back()->resetTelescopeCoordinates();
        for( int t = 0; t < fInputData[i].fNTelescopes; t++ )
        {
            if( !fStereoCompare.back()->setTelescopeCoordinates( fInputData[i].fTelX[t], fInputData[i].fTelY[t], fInputData[i].fTelZ[t] ) )
            {
                exit( EXIT_FAILURE );
            }
        }
        if( fInputData[i].fWobbleFromDataTree )
        {
            fStereoCompare.back()->setWobbleFromDataTree();
        }
        // sims: histogram for AZ weighting
        if( fInputData[i].fType == "SIMS" )
        {
            fStereoCompare[i]->setAzimuthWeightingHistogram( hAzOn );
        }
        // fill histograms
        fStereoCompare.back()->fillHistograms( fInputData[i].fFileName, fSingleTelescopeCuts );
        fStereoCompare.back()->writeHistograms( fOutputfile );
        
        if( fInputData[i].fType == "ON" )
        {
            fStereoCompareOn = fStereoCompare.back();
        }
        else if( fInputData[i].fType == "OFF" )
        {
            fStereoCompareOff = fStereoCompare.back();
        }
        cout << endl;
    }
    
    ////////////////////////////////////////
    // calculate difference histograms
    cout << "DIFF" << endl;
    cout << "----" << endl;
    VDataMCComparision* fDiff = new VDataMCComparision( "DIFF", iNT, fCalculateMVACut );
    // assume 5 background regions
    fDiff->setOnOffHistograms( fStereoCompareOn, fStereoCompareOff, 1. / 5. );
    fDiff->writeHistograms( fOutputfile );
}
