/* \file writeCTAWPPhysSensitivityTree summary sensitivity tree for different sub arrays

*/

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

//////////////////////////////////////
// telescope data
class VTelescopeData
{
    public:
    
        ULong64_t fTelType;
        string    fTelTypeName;
        float     fTel_x;
        float     fTel_y;
        
        VTelescopeData();
        ~VTelescopeData() {}
};

VTelescopeData::VTelescopeData()
{
    fTelType = 0;
    fTelTypeName = "";
    fTel_x = 0.;
    fTel_y = 0.;
};


//////////////////////////////////////
// sensitivity tree with all results
class VSensitivityTree
{
    private:
    
        TFile* fOutputFile;
        TTree* fDataTree;
        
        //////////////////////////////////////
        // variables defining the array
        Char_t fArray[300];
        int fNTel;
        // hardwired number of telescopes
        // prod2: 0 = LST, 1 = MST, 2 = 7m SST, 3 = 4m SST, 4 = SC-MST
        // static const int fNumTelTypes = 5;
        // prod3: 0 = LST, 1 = F-MST, 2 = N-MST, 3 = G-SST, 4 = D-SST, 5 = A-SST, 6 = SC-MST
        // prod5: 0 = LST, 1 = F-MST, 2 = N-MST, 3 = SST, 4 = SC-MST
        static const int fNumTelTypes = 5;
        int fNTelType[fNumTelTypes];                    // types of telescopes
        float fAverageTelescopeDistance[fNumTelTypes];  // average telescope distance to closest telescope
        float fMedianTelescopeDistance[fNumTelTypes];  // median telescope distance to closest telescope
        int fObservingTime_s;
        float fOffset_deg;
        unsigned int fTelMult;
        unsigned int fTelMultLST;
        unsigned int fTelMultMST;
        unsigned int fTelMultSST;
        unsigned int fTelMultSCMST;
        
        string fPointingDirection;
        int    fAnalysisID;
        
        ///////////////////////////////////////
        // variables defining the performance
        
        // 5 bins per decade, from -2 to 2 = 21
        static const int fNEnergyBins = 21;
        float fEnergy_logTeV[fNEnergyBins];      // mid point of log energy bin
        
        vector< string > fVarName;
        vector< float* > fVar;
        vector< float* > fVarError;

        static const int fN2DBins = 10000000;
        vector< unsigned int > fVar2DNbins;
        vector< unsigned int > fVar2DNbinsx;
        vector< unsigned int > fVar2DNbinsy;
        vector< float* > fVar2DXaxis;
        vector< float* > fVar2DYaxis;
        vector< string > fVar2DName;
        vector< float* > fVar2D;
        
        void reset();
        
    public:
    
        VSensitivityTree();
        ~VSensitivityTree() {}
        bool fillEvent( string iSite, 
                        string iSubArrayName, 
                        int iObservingTime_s,
                        int iTelMultiplicity,
                        int iTelMultiplicityLST,
                        int iTelMultiplicityMST,
                        int iTelMultiplicitySST,
                        int iTelMultiplicitySCST,
                        string iWPPhysDirectory,
                        string iMSTDateID );
        bool initialize( string iOutputFileName, string iName );
        void setFileNameParameters( string iPointingDirection, int iAnalysisID )
        {
            fPointingDirection = iPointingDirection;
            fAnalysisID = iAnalysisID;
        }
        bool terminate();
};

VSensitivityTree::VSensitivityTree()
{
    fOutputFile = 0;
    fDataTree = 0;
    fObservingTime_s = 0;
    fOffset_deg = 0.;
    fTelMult = 0;
    fTelMultLST = 0;
    fTelMultMST = 0;
    fTelMultSST = 0;
    fTelMultSCMST = 0;
    
    fAnalysisID = 0;
    fPointingDirection = "_0deg";
    
    // variable name
    // --> correspond to file names in 
    //     sensitivity files
    // --> correspond to branch names in
    //     output tree
    fVarName.push_back( "DiffSens" );
    fVarName.push_back( "DiffSensCU" );
    fVarName.push_back( "IntSens" );
    fVarName.push_back( "IntSensCU" );
    fVarName.push_back( "AngRes" );
    fVarName.push_back( "AngRes80" );
    fVarName.push_back( "AngRes95" );
    fVarName.push_back( "ERes" );
    fVarName.push_back( "Ebias" );
    fVarName.push_back( "BGRate" );
    fVarName.push_back( "ProtRate" );
    fVarName.push_back( "ElecRate" );
    fVarName.push_back( "BGRatePerSqDeg" );
    fVarName.push_back( "ProtRateSqDeg" );
    fVarName.push_back( "ElecRateSqDeg" );
    fVarName.push_back( "EffectiveArea" );
    fVarName.push_back( "EffectiveAreaEtrue" );
    fVarName.push_back( "EffectiveArea80" );
    fVarName.push_back( "EffectiveAreaEtrueNoTheta2cut" );
    
    for( unsigned int i = 0; i < fVarName.size(); i++ )
    {
        fVar.push_back( new float[fNEnergyBins] );
        fVarError.push_back( new float[fNEnergyBins] );
    }

    // TH3F histgrams
    // (tricky, as variable binning: assuming
    //  a huge value)
    fVar2DName.push_back( "AngularPSF2DEtrue" );
    fVar2DName.push_back( "MigMatrixNoTheta2cut" );
    for( unsigned int i = 0; i < fVar2DName.size(); i++ )
    {
        fVar2DNbins.push_back( 0 );
        fVar2DNbinsx.push_back( 0 );
        fVar2DNbinsy.push_back( 0 );
        fVar2D.push_back( new float[fN2DBins] );
        fVar2DXaxis.push_back( new float[fN2DBins] );
        fVar2DYaxis.push_back( new float[fN2DBins] );
    }
    
    reset();
}

void VSensitivityTree::reset()
{
    sprintf( fArray, "noArrayName" );
    fNTel = 0;
    fObservingTime_s = 0;
    fOffset_deg = 0.;
    fTelMult = 0;
    fTelMultLST = 0;
    fTelMultMST = 0;
    fTelMultSST = 0;
    fTelMultSCMST = 0;
    for( int i = 0; i < fNumTelTypes; i++ )
    {
        fNTelType[i] = 0;
        fAverageTelescopeDistance[i] = 0.;
        fMedianTelescopeDistance[i] = 0.;
    }
    
    for( int i = 0; i < fNEnergyBins; i++ )
    {
        fEnergy_logTeV[i] = 0.;
    }
    for( unsigned int v = 0; v < fVar.size(); v++ )
    {
        for( int i = 0; i < fNEnergyBins; i++ )
        {
            fVar[v][i] = 0.;
            fVarError[v][i] = 0.;
        }
    }
    for( unsigned int v = 0; v < fVar2D.size(); v++ )
    {
        fVar2DNbins[v] = 0;
        fVar2DNbinsx[v] = 0;
        fVar2DNbinsy[v] = 0;
        for( int i = 0; i < fN2DBins; i++ )
        {
            fVar2D[v][i] = 0.;
            fVar2DXaxis[v][i] = 0.;
            fVar2DYaxis[v][i] = 0.;
        }
    }
}

/*
 *
 * read array layouts and sensitivity parameters from mscw and phys files
 *
 */
bool VSensitivityTree::fillEvent( string iSite,
                                  string iSubArrayName, 
                                  int iObservingTime_s,
                                  int iTelMultiplicity,
                                  int iTelMultiplicityLST,
                                  int iTelMultiplicityMST,
                                  int iTelMultiplicitySST,
                                  int iTelMultiplicitySCMST,
                                  string iWPPhysDirectory,
                                  string iMSTDateID )
{
    // directories and file names
    //string iDataDirectory = "$CTA_USER_DATA_DIR/analysis/AnalysisData/" + iSite;
    string iDataDirectory = "/lustre/fs22/group/cta/users/maierg/analysis/AnalysisData/" + iSite;

    ////////////////////////////////////////////////////////////////////
    // open WP Phys file - return if not available
    //
    // WP Physfile
    ostringstream iWPPhysFile_st;
    iWPPhysFile_st << "DESY." << iMSTDateID << ".V3.ID" << fAnalysisID;
    iWPPhysFile_st << fPointingDirection;
    iWPPhysFile_st << "NIM" << iTelMultiplicity;
    iWPPhysFile_st << "LST" << iTelMultiplicityLST;
    iWPPhysFile_st << "MST" << iTelMultiplicityMST;
    iWPPhysFile_st << "SST" << iTelMultiplicitySST;
    iWPPhysFile_st << "SCMST" << iTelMultiplicitySCMST;
    iWPPhysFile_st << "." << iSite << ".";
    string iWPPhysFile = iWPPhysFile_st.str();
    ostringstream iObservingTime;
    iObservingTime << iObservingTime_s << "s";

    ostringstream iWFile;
    // iWFile << iWPPhysDirectory << "-NIM" << iTelMultiplicity << "/";
    iWFile << iWPPhysDirectory << "/";
    iWFile << iWPPhysFile;
    iWFile << iSubArrayName << "." << iObservingTime.str();
    iWFile << ".root";

    cout << "reading IRFs from " << iWFile.str() << endl;

    // test if root file exists
    // (prevents printout of root errors / warnings)
    FILE *test_file = fopen(iWFile.str().c_str(), "r");
    if( !test_file )
    {
         cout << "File not found: " << iWFile.str() << endl;
         return false;
    }
    fclose(test_file);

    // open IRF file
    TFile iP( iWFile.str().c_str() );
    if( iP.IsZombie() )
    {
        cout << "File oppening error: " << iWFile.str() << endl;
        return false;
    }

    /////////////////////////////////////////////////////
    // get telconfig
    
    
    // dates in file names
    string iEffDate;
    string iMSCWDate;
    
    // directory with MSCW files
    ostringstream iMSCWDirectory_st;
    iMSCWDirectory_st << iDataDirectory << "/";
    iMSCWDirectory_st << iSubArrayName << "/Analysis-ID" << fAnalysisID << "-";
    iMSCWDirectory_st << iMSTDateID << "/";
    // MSCW file name
    // (only used for telconfig tree, so pointing directory is not relevant)
    ostringstream iMSCWFile_st;
    iMSCWFile_st << "gamma_onSource." << iSubArrayName << "_ID" << fAnalysisID;
    iMSCWFile_st << "_180deg-" << iSite << "-1.mscw.root";
    
    cout << iMSCWDirectory_st.str() << endl;
    cout << iMSCWFile_st.str() << endl;
    cout << iWPPhysFile << endl;

    reset();
    
    // array name & ID
    sprintf( fArray, "%s", iSubArrayName.c_str() );
    
    ////////////////////////////////////////////////
    // read telconfig tree from MSCW files
    
    string iMFile = iMSCWDirectory_st.str() + "/" + iMSCWFile_st.str();
    cout << "reading telconfig from " << iMFile << endl;
    
    TFile iT( iMFile.c_str() );
    if( iT.IsZombie() )
    {
        cout << "Error reading file " << iMFile << endl;
        return false;
    }
    TTree* t = ( TTree* )iT.Get( "telconfig" );
    if( !t )
    {
        cout << "Error telconfig tree from " << iMFile << endl;
        return false;
    }
    cout << "Telconfig tree found with " << t->GetEntries() << " telescopes" << endl;
    
    ULong64_t iTelType = 0;
    float iTelX = 0.;
    float iTelY = 0.;
    
    t->SetBranchAddress( "TelType", &iTelType );
    t->SetBranchAddress( "TelX", &iTelX );
    t->SetBranchAddress( "TelY", &iTelY );
    
    fNTel = t->GetEntries();
    
    // telescope data
    vector< vector< VTelescopeData* > > fTelescopeData;
    vector< VTelescopeData* > i_TempTD;
    for( int i = 0; i < fNumTelTypes; i++ )
    {
        fTelescopeData.push_back( i_TempTD );
    }
    // loop over all entries to count different telescope types
    for( int i = 0; i < t->GetEntries(); i++ )
    {
        t->GetEntry( i );

        // LSTs (Type 0: 23m-LSTs)
        if( iTelType == 138704810 )
        {
            fNTelType[0]++;
            if( fTelescopeData.size() > 0 )
            {
                fTelescopeData[0].push_back( new VTelescopeData() );
                fTelescopeData[0].back()->fTelTypeName = "LST";
                fTelescopeData[0].back()->fTel_x = iTelX;
                fTelescopeData[0].back()->fTel_y = iTelY;
            }
        }
        // F-MST
        else if( iTelType == 10408618 )
        {
            fNTelType[1]++;
            if( fTelescopeData.size() > 1 )
            {
                fTelescopeData[1].push_back( new VTelescopeData() );
                fTelescopeData[1].back()->fTelTypeName = "FlashCam-MST";
                fTelescopeData[1].back()->fTel_x = iTelX;
                fTelescopeData[1].back()->fTel_y = iTelY;
            }
        }
        // N-MST
        else if( iTelType == 10408418 || iTelType == 10608418 )
        {
            fNTelType[2]++;
            if( fTelescopeData.size() > 1 )
            {
                fTelescopeData[2].push_back( new VTelescopeData() );
                fTelescopeData[2].back()->fTelTypeName = "NectarCam-MST";
                fTelescopeData[2].back()->fTel_x = iTelX;
                fTelescopeData[2].back()->fTel_y = iTelY;
            }
        }
        // SST
        else if( iTelType == 201409917 || iTelType == 201109916 )
        {
            fNTelType[3]++;
            if( fTelescopeData.size() > 2 )
            {
                fTelescopeData[3].push_back( new VTelescopeData() );
                fTelescopeData[3].back()->fTelTypeName = "SST";
                fTelescopeData[3].back()->fTel_x = iTelX;
                fTelescopeData[3].back()->fTel_y = iTelY;
            }
        }
        // SC-MSTs (type 4)
        else if( iTelType == 207308707 || iTelType == 205008707 )
        {
            fNTelType[4]++;
            if( fTelescopeData.size() > 5 )
            {
                fTelescopeData[4].push_back( new VTelescopeData() );
                fTelescopeData[4].back()->fTelTypeName = "MST-SCT";
                fTelescopeData[4].back()->fTel_x = iTelX;
                fTelescopeData[4].back()->fTel_y = iTelY;
            }
        }
        else
        {
            cout << "unknown telescope type: " << iTelType << endl;
        }
    }
    iT.Close();
    
    ////////////////////////
    // calculate average distance to closest telescope
    float i_dist = 0.;
    
    // loop over all telescope types
    for( unsigned int i = 0; i < fTelescopeData.size(); i++ )
    {
        float dist_av = 0.;
        float dist_N = 0.;
        // distance to all other telescopes of this type
        vector< double > iTelescopeDistances;
        for( unsigned int j = 0; j < fTelescopeData[i].size(); j++ )
        {
            float i_min = 1.e10;
            for( unsigned int k = 0; k < fTelescopeData[i].size(); k++ )
            {
                if( j == k )
                {
                    continue;
                }
                i_dist  = ( fTelescopeData[i][j]->fTel_x - fTelescopeData[i][k]->fTel_x ) 
                        * ( fTelescopeData[i][j]->fTel_x - fTelescopeData[i][k]->fTel_x );
                i_dist += ( fTelescopeData[i][j]->fTel_y - fTelescopeData[i][k]->fTel_y )
                        * ( fTelescopeData[i][j]->fTel_y - fTelescopeData[i][k]->fTel_y );
                i_dist  = sqrt( i_dist );
                if( i_dist < i_min )
                {
                    i_min = i_dist;
                }
            }
            if( i_min < 1.e9 && i_min > 0. )
            {
                dist_av += i_min;
                dist_N++;
                iTelescopeDistances.push_back( i_min );
            }
        }
        if( dist_N > 0. )
        {
            fAverageTelescopeDistance[i] = dist_av / dist_N;
        }
        if( iTelescopeDistances.size() > 0 )
        {
            sort( iTelescopeDistances.begin(), iTelescopeDistances.end() );
            if( iTelescopeDistances.size() % 2 == 0 )
            {
                fMedianTelescopeDistance[i] = 0.5 * ( iTelescopeDistances[iTelescopeDistances.size() / 2 - 1 ]
                                                    + iTelescopeDistances[iTelescopeDistances.size() / 2] );
            }
            else
            {
                fMedianTelescopeDistance[i] = iTelescopeDistances[iTelescopeDistances.size() / 2];
            }
        }
        else
        {
            fMedianTelescopeDistance[i] = 0.;
        }
    }
    
    ////////////////////////////////////////////////
    // read sensitivities

    fObservingTime_s = iObservingTime_s;
    fTelMult = iTelMultiplicity;
    fTelMultLST = iTelMultiplicityLST;
    fTelMultMST = iTelMultiplicityMST;
    fTelMultSST = iTelMultiplicitySST;
    fTelMultSCMST = iTelMultiplicitySCMST;
    
    // check consistency between energy bins
    // (energy binning is fixed)
    // set energies
    TH1F* h = ( TH1F* )iP.Get( "DiffSens" );
    if( h )
    {
        if( h->GetNbinsX() != fNEnergyBins )
        {
            cout << "error reading IRFs, different number of energy bins: ";
            cout << fNEnergyBins << "\t" << h->GetNbinsX() << endl;
            return false;
        }
        for( int i = 0; i < h->GetNbinsX(); i++ )
        {
            fEnergy_logTeV[i] = h->GetBinCenter( i );
        }
    }
    //////////////////////////////////////////////
    // get on-axis performance values
    fOffset_deg = 0.;
    for( unsigned int i = 0; i < fVarName.size(); i++ )
    {
        h = ( TH1F* )iP.Get( fVarName[i].c_str() );
        if( !h )
        {
            cout << "error finding histogram " << fVarName[i] << endl;
            continue;
        }
        for( int e = 0; e < fNEnergyBins; e++ )
        {
            int i_bin = h->FindBin( fEnergy_logTeV[e] );
            fVar[i][e] = h->GetBinContent( i_bin );
            fVarError[i][e] = h->GetBinError( i_bin );
        }
    }
    for( unsigned int i = 0; i < fVar2DName.size(); i++ )
    {
        TH2F *h2 = ( TH2F* )iP.Get( fVar2DName[i].c_str() );
        if( !h2 )
        {
            cout << "error finding histogram " << fVar2DName[i] << endl;
            continue;
        }
        // fill 2D histogram into 1D array:
        // e1 + e2 * fxbins
        fVar2DNbinsx[i] = h2->GetXaxis()->GetNbins();
        fVar2DNbinsy[i] = h2->GetYaxis()->GetNbins();
        fVar2DNbins[i] = fVar2DNbinsx[i] * fVar2DNbinsy[i];
        for( unsigned int e1 = 0; e1 < fVar2DNbinsx[i]; e1++ )
        {
            fVar2DXaxis[i][e1] = h2->GetXaxis()->GetBinCenter( e1+1 );
            for( unsigned int e2 = 0; e2 < fVar2DNbinsy[i]; e2++ )
            {
                fVar2DYaxis[i][e2] = h2->GetYaxis()->GetBinCenter( e2+1 );

                fVar2D[i][e1+e2*fVar2DNbinsx[i]] = h2->GetBinContent( e1+1, e2+1 );
            }
        }
    }

    // fill results
    fDataTree->Fill();
    
    //////////////////////////////////////////////
    // get off-axis performances
    TH2F* h2 = ( TH2F* )iP.Get( "DiffSens_offaxis" );
    if( h2 )
    {
        for( int w = 1; w <= h2->GetNbinsY(); w++ )
        {
            fOffset_deg = h2->GetYaxis()->GetBinCenter( w );
            for( unsigned int i = 0; i < fVarName.size(); i++ )
            {
                string iHName = fVarName[i] + "_offaxis";
                h = ( TH1F* )iP.Get( iHName.c_str() );
                if( !h )
                {
                    cout << "error finding histogram " << iHName << endl;
                    continue;
                }
                for( int e = 0; e < fNEnergyBins; e++ )
                {
                    int i_bin = h->GetXaxis()->FindBin( fEnergy_logTeV[e] );
                    fVar[i][e] = h->GetBinContent( i_bin, w );
                    fVarError[i][e] = h->GetBinError( i_bin, w );
                }
            }
            for( unsigned int i = 0; i < fVar2DName.size(); i++ )
            {
                string iHName = fVar2DName[i] + "_offaxis";
                TH3F *h2 = ( TH3F* )iP.Get( iHName.c_str() );
                if( !h2 )
                {
                    cout << "error finding histogram " << fVar2DName[i] << endl;
                    continue;
                }
                // fill 2D histogram into 1D array:
                // e1 + e2 * fxbins
                fVar2DNbinsx[i] = h2->GetXaxis()->GetNbins();
                fVar2DNbinsy[i] = h2->GetYaxis()->GetNbins();
                fVar2DNbins[i] = fVar2DNbinsx[i] * fVar2DNbinsy[i];
                for( unsigned int e1 = 0; e1 < fVar2DNbinsx[i]; e1++ )
                {
                    fVar2DXaxis[i][e1] = h2->GetXaxis()->GetBinCenter( e1+1 );
                    for( unsigned int e2 = 0; e2 < fVar2DNbinsy[i]; e2++ )
                    {
                        fVar2DYaxis[i][e2] = h2->GetYaxis()->GetBinCenter( e2+1 );

                        fVar2D[i][e1+e2*fVar2DNbinsx[i]] = h2->GetBinContent( e1+1, e2+1, w );
                    }
                }
            }
            // fill results
            fDataTree->Fill();
        }
    }
    iP.Close();
    
    ////////////////////////////////////////////////
    
    return true;
}

/*
 *   open output file and initialize tree structure
 *
 */
bool VSensitivityTree::initialize( string iOutputFileName, string iName )
{
    fOutputFile = new TFile( iOutputFileName.c_str(), "RECREATE" );
    if( fOutputFile->IsZombie() )
    {
        cout << "VSensitivityTree::initialize error opening output file" << endl;
        return false;
    }
    
    // define data tree
    char hname[200];
    char htitle[200];
    
    fDataTree = new TTree( "IRFData", iName.c_str() );
    
    fDataTree->Branch( "Array", &fArray, "Array/C" );

    // telescope positions
    fDataTree->Branch( "NTel", &fNTel, "NTel/I" );
    sprintf( hname, "NTelType[%d]", fNumTelTypes );
    sprintf( htitle, "NTelType[%d]/I", fNumTelTypes );
    fDataTree->Branch( hname, fNTelType, htitle );
    
    sprintf( hname, "AverageTelescopeDistance[%d]", fNumTelTypes );
    sprintf( htitle, "AverageTelescopeDistance[%d]/F", fNumTelTypes );
    fDataTree->Branch( hname, fAverageTelescopeDistance, htitle );
    
    sprintf( hname, "MedianTelescopeDistance[%d]", fNumTelTypes );
    sprintf( htitle, "MedianTelescopeDistance[%d]/F", fNumTelTypes );
    fDataTree->Branch( hname, fMedianTelescopeDistance, htitle );
    // multiplicity cut parameters
    fDataTree->Branch( "TelMult", &fTelMult, "TelMult/i" );
    fDataTree->Branch( "TelMultLST", &fTelMultLST, "TelMultLST/i" );
    fDataTree->Branch( "TelMultMST", &fTelMultMST, "TelMultMST/i" );
    fDataTree->Branch( "TelMultSST", &fTelMultSST, "TelMultSST/i" );
    fDataTree->Branch( "TelMultSCMST", &fTelMultSCMST, "TelMultSCMST/i" );
    // observational details
    fDataTree->Branch( "ObsTime_s", &fObservingTime_s, "ObsTime_s/I" );
    fDataTree->Branch( "Offset_deg", &fOffset_deg, "Offset_deg/F" );
    // IRFs 
    sprintf( hname, "Energy_logTeV[%d]", fNEnergyBins );
    sprintf( htitle, "Energy_logTeV[%d]/F", fNEnergyBins );
    fDataTree->Branch( hname, fEnergy_logTeV, htitle );
    
    for( unsigned int i = 0; i < fVar.size(); i++ )
    {
        sprintf( hname, "%s[%d]", fVarName[i].c_str(), fNEnergyBins );
        sprintf( htitle, "%s[%d]/F", fVarName[i].c_str(), fNEnergyBins );
        fDataTree->Branch( hname, fVar[i], htitle );
        
        sprintf( hname, "%sError[%d]", fVarName[i].c_str(), fNEnergyBins );
        sprintf( htitle, "%sError[%d]/F", fVarName[i].c_str(), fNEnergyBins );
        fDataTree->Branch( hname, fVarError[i], htitle );
    }
    for( unsigned int i = 0; i < fVar2D.size(); i++ )
    {
        sprintf( hname, "Nbins_%s", fVar2DName[i].c_str() );
        sprintf( htitle, "Nbins_%s/i", fVar2DName[i].c_str() );
        fDataTree->Branch( hname, &fVar2DNbins[i], htitle );

        sprintf( hname, "Nbins_%s_x", fVar2DName[i].c_str() );
        sprintf( htitle, "Nbins_%s_x/i", fVar2DName[i].c_str() );
        fDataTree->Branch( hname, &fVar2DNbinsx[i], htitle );

        sprintf( hname, "Axis_%s_x[Nbins_%s_x]", fVar2DName[i].c_str(), fVar2DName[i].c_str() );
        sprintf( htitle, "Axis_%s_x[Nbins_%s_x]/F", fVar2DName[i].c_str(), fVar2DName[i].c_str() );
        fDataTree->Branch( hname, fVar2DXaxis[i], htitle );

        sprintf( hname, "Nbins_%s_y", fVar2DName[i].c_str() );
        sprintf( htitle, "Nbins_%s_y/i", fVar2DName[i].c_str() );
        fDataTree->Branch( hname, &fVar2DNbinsy[i], htitle );

        sprintf( hname, "Axis_%s_y[Nbins_%s_y]", fVar2DName[i].c_str(), fVar2DName[i].c_str() );
        sprintf( htitle, "Axis_%s_y[Nbins_%s_y]/F", fVar2DName[i].c_str(), fVar2DName[i].c_str() );
        fDataTree->Branch( hname, fVar2DYaxis[i], htitle );

        sprintf( hname, "%s[Nbins_%s]", fVar2DName[i].c_str(), fVar2DName[i].c_str() );
        sprintf( htitle, "%s[Nbins_%s]/F", fVar2DName[i].c_str(), fVar2DName[i].c_str() );
        fDataTree->Branch( hname, fVar2D[i], htitle );
    }
    
    return true;
}

bool VSensitivityTree::terminate()
{
    if( fOutputFile && fOutputFile->cd() )
    {
        if( fDataTree )
        {
            fDataTree->Write();
        }
        
        cout << "closing output file " << fOutputFile->GetName() << endl;
        fOutputFile->Close();
    }
    
    return true;
}


/////////////////////////////////////////////////////
// main
/////////////////////////////////////////////////////


int main( int argc, char* argv[] )
{
    /////////////////////
    // input parameters
    if( argc != 8 )
    {
        cout << endl;
        cout << "./writeCTAWPPhysSensitivityTree <list of sub array> <data set> <output file> <pointing direction> <ID> <PHYS file directory> <MSCW data identifier>" << endl;
        cout << endl;
        cout << "  combine a list of PHYS root files into one single tree" << endl;
        cout << endl;
        cout << endl;
        cout << "<pointing direction> (e.g. _180deg, _0deg, Average )" << endl;
        cout << "<ID>                 (e.g. 0, 1, 2, ...)" << endl;
        cout << "<PHYS file directory> relative to analysis directory" << endl;
        cout << "<MSCW data identifier> MSCW data used (e.g., g20210921 for Analysis-ID0-g20210921" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    cout << endl;
    cout << "writeCTAWPPhysSensitivityTree"  << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    string fSubArrayList = argv[1];
    string fSite = argv[2];
    string fOutputfile = argv[3];
    string fPointingDirection = argv[4];
    if( fPointingDirection == "Average" )
    {
        fPointingDirection = "";
    }
    else if( fPointingDirection == "0deg" )
    {
        fPointingDirection = "_0deg";
    }
    else if( fPointingDirection == "180deg" )
    {
        fPointingDirection = "_180deg";
    }
    int    fAnalysisID = atoi( argv[5] );
    string fInputDirectory = argv[6];
    string fMSTDateID = argv[7];
    
    // loop over all sub arrays and load files
    vector< string > fSubArray;
    
    ifstream is;
    is.open( fSubArrayList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "error reading list of arrays from: " << endl;
        cout << fSubArrayList << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "reading list of arrays from " << endl;
    cout << fSubArrayList << endl;
    string iLine;
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            fSubArray.push_back( iLine );
        }
    }
    is.close();
    
    cout << "total number of subarrays: " << fSubArray.size() << endl;

    /////////////////////////////////////////
    // initialize and fill data tree
    
    VSensitivityTree* fData = new VSensitivityTree();
    if( !fData->initialize( fOutputfile, "test data" ) )
    {
        exit( EXIT_FAILURE );
    }
    fData->setFileNameParameters( fPointingDirection, fAnalysisID );
    
    // hardwired observing time
    vector< int > iObservingTimeVector;
    iObservingTimeVector.push_back( 180000 );
    iObservingTimeVector.push_back( 1800 );

    // hardwired telescope multiplicity
    //
    // (LST multiplicity is overwritten below)
    vector< int > iTelMultiplicityLST;
    iTelMultiplicityLST.push_back( 3 );
    vector< int > iTelMultiplicityMST;
    iTelMultiplicityMST.push_back( 2 );
    iTelMultiplicityMST.push_back( 3 );
    iTelMultiplicityMST.push_back( 4 );
    iTelMultiplicityMST.push_back( 5 );
    iTelMultiplicityMST.push_back( 6 );
    vector< int > iTelMultiplicitySST;
    iTelMultiplicitySST.push_back( 2 );
    iTelMultiplicitySST.push_back( 3 );
    iTelMultiplicitySST.push_back( 4 );
    iTelMultiplicitySST.push_back( 5 );
    iTelMultiplicitySST.push_back( 6 );
    // (SCMST multiplicity is overwritten below)
    vector< int > iTelMultiplicitySCMST;
    iTelMultiplicitySCMST.push_back( 3 );

    /////////////////////////////////////////////////////////
    // fill events for complete parameter space
    // array loop
    for( unsigned int i = 0; i < fSubArray.size(); i++ )
    {
        // observing time loop
        for( unsigned int t = 0; t < iObservingTimeVector.size(); t++ )
        {
            // telescope multiplicity (per telescope type)
            for( unsigned l = 0; l < iTelMultiplicityLST.size(); l++ )
            {
                for( unsigned m = 0; m < iTelMultiplicityMST.size(); m++ )
                {
                    for( unsigned s = 0; s < iTelMultiplicitySST.size(); s++ )
                    {
                        for( unsigned c = 0; c < iTelMultiplicitySCMST.size(); c++ )
                        {
                            // TMPTMPTMP
                            // (requires NMSTs == NSSTs)
                            // if( iTelMultiplicityMST[m] != iTelMultiplicitySST[s] ) continue;

                            // TMPTMPTMP
                            iTelMultiplicityLST[l] = std::min( iTelMultiplicityMST[m], iTelMultiplicitySST[s] );
                            iTelMultiplicitySCMST[c] = std::min( iTelMultiplicityMST[m], iTelMultiplicitySST[s] );
                            // TMPTMPTMP

                            unsigned int iTelMultiplicity = 
                                        std::min( iTelMultiplicityLST[l], 
                                                std::min( iTelMultiplicityMST[m],
                                                std::min( iTelMultiplicitySST[s], iTelMultiplicitySCMST[c] ) ) );

                            cout << endl << endl;
                            cout << "now filling array layout " << fSubArray[i];
                            cout << " (obs time " << iObservingTimeVector[t] << " s), ";
                            cout << "telescope multiplicity " << iTelMultiplicity;
                            cout << ", LST: " << iTelMultiplicityLST[l];
                            cout << ", MST: " << iTelMultiplicityMST[m];
                            cout << ", SST: " << iTelMultiplicitySST[s];
                            cout << ", SCMST: " << iTelMultiplicitySCMST[c] << ")" << endl;
                            cout << "=========================================" << endl;
                            cout << endl;
                            
                            fData->fillEvent( fSite, fSubArray[i], iObservingTimeVector[t], 
                                   iTelMultiplicity,
                                   iTelMultiplicityLST[l],
                                   iTelMultiplicityMST[m],
                                   iTelMultiplicitySST[s],
                                   iTelMultiplicitySCMST[c],
                                   fInputDirectory, fMSTDateID );
                        }
                     }
                }
            }
        }
    }
        
    fData->terminate();
    
    return 0;
}
