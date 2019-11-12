/* \file writeCTAWPPhysSensitivityTree summary sensitivity tree for different sub arrays

*/

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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
        unsigned int fArrayID;
        unsigned int fScaling;
        int fNTel;
        // hardwired number of telescopes
        // prod2: 0 = LST, 1 = MST, 2 = 7m SST, 3 = 4m SST, 4 = SC-MST
        // static const int fNumTelTypes = 5;
        // prod3: 0 = LST, 1 = F-MST, 2 = N-MST, 3 = G-SST, 4 = D-SST, 5 = A-SST, 6 = SC-MST
        static const int fNumTelTypes = 7;
        int fNTelType[fNumTelTypes];                    // types of telescopes
        float fAverageTelescopeDistance[fNumTelTypes];  // average telescope distance to closest telescope
        float fMedianTelescopeDistance[fNumTelTypes];  // median telescope distance to closest telescope
        int fObservingTime_s;
        float fOffset_deg;
        
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
        
        void reset();
        
    public:
    
        VSensitivityTree();
        ~VSensitivityTree() {}
        bool fillEvent( string iSite, string iSubArrayName, int iObservingTime_s, unsigned int iSubArrayID, unsigned int iScaling );
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
    
    fAnalysisID = 0;
    fPointingDirection = "_0deg";
    
    // variable name
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
    
    for( unsigned int i = 0; i < fVarName.size(); i++ )
    {
        fVar.push_back( new float[fNEnergyBins] );
        fVarError.push_back( new float[fNEnergyBins] );
    }
    
    reset();
}

void VSensitivityTree::reset()
{
    sprintf( fArray, "noArrayName" );
    fArrayID = 0;
    fNTel = 0;
    fObservingTime_s = 0;
    fOffset_deg = 0.;
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
}

/*
 *
 * read array layouts and sensitivity parameters from mscw and phys files
 *
 */
bool VSensitivityTree::fillEvent( string iSite, string iSubArrayName, int iObservingTime_s, unsigned int iSubArrayID, unsigned int iScaling )
{
    ////////////////////////////////////
    // HARDWIRED VALUES - START //
    // PROD2
    /*    string iMSCWDirectory = "/lustre/fs9/group/cta/users/maierg/CTA/analysis/AnalysisData/prod2-Aar2014-NS/ARRAYLAYOUT/Analysis-ID0-d20150322-onAxis/";
        string iMSCWFile = "gamma_onSource.ARRAYLAYOUT_ID0POINTINGDIRECTION-prod2-Aar2014-NS-1.mscw.root";
        string iWPPhysDirectory = "/lustre/fs9/group/cta/users/maierg/CTA//analysis/WPPhys.20150901/";
    //    string iWPPhysFile = "DESY.d20150824.Erec1.L1.ID0POINTINGDIRECTIONNIM2.prod2-Aar2014-NN.";
        string iWPPhysFile = "DESY.d20150824.Erec1.L1.ID0POINTINGDIRECTIONNIM2.prod2-Aar-40deg-NN.";
        string iTelescopePointingDirection = "_0deg";  */
    // PROD3
    char hname[2000];
    
    string iDataDirectory = "/lustre/fs9/group/cta/users/maierg/CTA/analysis/AnalysisData/";
    string iWPPhysDirectory = "/lustre/fs16/group/cta/users/maierg/CTA/analysis/WPPhys/";
    
    // dates in file names
    string iEffDate;
    string iMSCWDate;
    if( iSite == "prod3-paranalp05-NN" )
    {
        iEffDate = "d20160304";
        iMSCWDate = "d20160304";
    }
    else if( iSite == "prod3-paranalp05-40deg-NN" )
    {
        iEffDate = "d20160304";
        iMSCWDate = "d20160224";
    }
    else if( iSite == "prod3-LaPalmap05-NN" )
    {
        iEffDate = "d20160304";
        iMSCWDate = "d20160304";
    }
    else
    {
        cout << "Site not found" << endl;
        exit( EXIT_FAILURE );
    }
    
    
    // directory with MSCW files
    ostringstream iMSCWDirectory_st;
    iMSCWDirectory_st << iDataDirectory << "/" << iSite << "/";
    iMSCWDirectory_st << "ARRAYLAYOUT/Analysis-ID" << fAnalysisID << "-";
    iMSCWDirectory_st << iMSCWDate  << "/";
    string iMSCWDirectory = iMSCWDirectory_st.str();
    // MSCW file name
    ostringstream iMSCWFile_st;
    iMSCWFile_st << "gamma_onSource.ARRAYLAYOUT_ID" << fAnalysisID;
    iMSCWFile_st << "POINTINGDIRECTION-" << iSite << "-1.mscw.root";
    string iMSCWFile = iMSCWFile_st.str();
    // WP Physfile
    ostringstream iWPPhysFile_st;
    iWPPhysFile_st << "DESY." << iEffDate << ".V5.ID" << fAnalysisID;
    iWPPhysFile_st << "POINTINGDIRECTIONNIM2." << iSite << ".";
    string iWPPhysFile = iWPPhysFile_st.str();
    
    cout << iMSCWDirectory << endl;
    cout << iMSCWFile << endl;
    cout << iWPPhysFile << endl;
    
    ////////////////////////////////////
    // HARDWIRED VALUES - STOPP //
    
    reset();
    
    // array name & ID
    sprintf( fArray, "%s", iSubArrayName.c_str() );
    fArrayID = iSubArrayID;
    fScaling = iScaling;
    
    ////////////////////////////////////////////////
    // read telconfig tree from MSCW files
    
    string iMFile = iMSCWDirectory.replace( iMSCWDirectory.find( "ARRAYLAYOUT" ), std::string( "ARRAYLAYOUT" ).length(), iSubArrayName );
    iMFile += iMSCWFile.replace( iMSCWFile.find( "ARRAYLAYOUT" ), std::string( "ARRAYLAYOUT" ).length(), iSubArrayName );
    if( fPointingDirection.size() > 0 )
    {
        iMFile  = iMFile.replace( iMFile.find( "POINTINGDIRECTION" ), std::string( "POINTINGDIRECTION" ).length(), fPointingDirection );
    }
    else
    {
        iMFile  = iMFile.replace( iMFile.find( "POINTINGDIRECTION" ), std::string( "POINTINGDIRECTION" ).length(), "_180deg" );
    }
    
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
        
        ///////////////////////////////////////////
        // prod2 (5 different telescope types)
        if( fNumTelTypes == 5 )
        {
            // LSTs (Type 0: 23m-LSTs)
            if( iTelType == 138704810 || iTelType == 141305009 || iTelType == 141305109 )
            {
                fNTelType[0]++;
                if( fTelescopeData.size() > 0 )
                {
                    fTelescopeData[0].push_back( new VTelescopeData() );
                    fTelescopeData[0].back()->fTelTypeName = "23m-LST";
                    fTelescopeData[0].back()->fTel_x = iTelX;
                    fTelescopeData[0].back()->fTel_y = iTelY;
                }
            }
            // standard MSTs (type 1: 12m-MSTs)
            else if( iTelType == 10007818 || iTelType == 10408418 || iTelType == 10008118 )
            {
                fNTelType[1]++;
                if( fTelescopeData.size() > 1 )
                {
                    fTelescopeData[1].push_back( new VTelescopeData() );
                    fTelescopeData[1].back()->fTelTypeName = "12m-MST";
                    fTelescopeData[1].back()->fTel_x = iTelX;
                    fTelescopeData[1].back()->fTel_y = iTelY;
                }
            }
            // standard SSTs (type 2: 4m SSTs)
            else if( iTelType == 201509515 )
            {
                fNTelType[2]++;
                if( fTelescopeData.size() > 2 )
                {
                    fTelescopeData[2].push_back( new VTelescopeData() );
                    fTelescopeData[2].back()->fTelTypeName = "4m-SC-SST";
                    fTelescopeData[2].back()->fTel_x = iTelX;
                    fTelescopeData[2].back()->fTel_y = iTelY;
                }
            }
            // large SSTs (type 3: 7m SSTs)
            else if( iTelType == 3709725 || iTelType == 3709425 || iTelType == 3710125 )
            {
                fNTelType[3]++;
                if( fTelescopeData.size() > 3 )
                {
                    fTelescopeData[3].push_back( new VTelescopeData() );
                    fTelescopeData[3].back()->fTelTypeName = "7m-DC-SST";
                    fTelescopeData[3].back()->fTel_x = iTelX;
                    fTelescopeData[3].back()->fTel_y = iTelY;
                }
            }
            // SC-MSTs (type 4)
            else if( iTelType == 201509515 )
            {
                fNTelType[4]++;
                if( fTelescopeData.size() > 4 )
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
        // prod3 (7 telescope types)
        else
        {
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
            else if( iTelType == 10408418 )
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
            // G-SST
            else if( iTelType == 201309415 || iTelType == 201309316 )
            {
                fNTelType[3]++;
                if( fTelescopeData.size() > 2 )
                {
                    fTelescopeData[3].push_back( new VTelescopeData() );
                    fTelescopeData[3].back()->fTelTypeName = "G-SST";
                    fTelescopeData[3].back()->fTel_x = iTelX;
                    fTelescopeData[3].back()->fTel_y = iTelY;
                }
            }
            // D-SST
            else if( iTelType == 908924 || iTelType == 909924 )
            {
                fNTelType[4]++;
                if( fTelescopeData.size() > 3 )
                {
                    fTelescopeData[4].push_back( new VTelescopeData() );
                    fTelescopeData[4].back()->fTelTypeName = "D-SST";
                    fTelescopeData[4].back()->fTel_x = iTelX;
                    fTelescopeData[4].back()->fTel_y = iTelY;
                }
            }
            // A-SST
            else if( iTelType == 201510718 || iTelType == 201511619 )
            {
                fNTelType[5]++;
                if( fTelescopeData.size() > 4 )
                {
                    fTelescopeData[5].push_back( new VTelescopeData() );
                    fTelescopeData[5].back()->fTelTypeName = "A-SST";
                    fTelescopeData[5].back()->fTel_x = iTelX;
                    fTelescopeData[5].back()->fTel_y = iTelY;
                }
            }
            // SC-MSTs (type 4)
            else if( iTelType == 207308707 )
            {
                fNTelType[6]++;
                if( fTelescopeData.size() > 5 )
                {
                    fTelescopeData[6].push_back( new VTelescopeData() );
                    fTelescopeData[6].back()->fTelTypeName = "MST-SCT";
                    fTelescopeData[6].back()->fTel_x = iTelX;
                    fTelescopeData[6].back()->fTel_y = iTelY;
                }
            }
            else
            {
                cout << "unknown telescope type: " << iTelType << endl;
            }
        }
    }
    iT.Close();
    
    ////////////////////////
    // calculate average distance to closest telescope 3 telescopes
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
                i_dist  = ( fTelescopeData[i][j]->fTel_x - fTelescopeData[i][k]->fTel_x ) * ( fTelescopeData[i][j]->fTel_x - fTelescopeData[i][k]->fTel_x );
                i_dist += ( fTelescopeData[i][j]->fTel_y - fTelescopeData[i][k]->fTel_y ) * ( fTelescopeData[i][j]->fTel_y - fTelescopeData[i][k]->fTel_y );
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
    sprintf( hname, "%ds", iObservingTime_s );
    string iObservingTime = hname;
    
    string iWFile  = iWPPhysDirectory;
    iWFile        += iWPPhysFile.replace( iWPPhysFile.find( "POINTINGDIRECTION" ), std::string( "POINTINGDIRECTION" ).length(), fPointingDirection );
    iWFile        += iSubArrayName + "." + iObservingTime + ".root";
    
    cout << "reading IRFs from " << iWFile << endl;
    
    // open IRF file
    TFile iP( iWFile.c_str() );
    if( iP.IsZombie() )
    {
        cout << "Error reading file " << iWFile << endl;
        return false;
    }
    
    
    // check energy bins
    TH1F* h = ( TH1F* )iP.Get( "DiffSens" );
    if( h )
    {
        if( h->GetNbinsX() != fNEnergyBins )
        {
            cout << "error reading IRFs, different number of energy bins: " << fNEnergyBins << "\t" << h->GetNbinsX() << endl;
            return false;
        }
        for( int i = 0; i < h->GetNbinsX(); i++ )
        {
            fEnergy_logTeV[i] = h->GetBinCenter( i );
        }
    }
    //////////////////////////////////////////////
    // get on-axis sensitivies
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
    // fill results
    fDataTree->Fill();
    
    //////////////////////////////////////////////
    // get off-axis sensitivies
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
    fDataTree->Branch( "ArrayID", &fArrayID, "ArrayID/i" );
    fDataTree->Branch( "Scaling", &fScaling, "Scaling/i" );
    fDataTree->Branch( "NTel", &fNTel, "NTel/I" );
    sprintf( hname, "NTelType[%d]", fNumTelTypes );
    sprintf( htitle, "NTelType[%d]/I", fNumTelTypes );
    fDataTree->Branch( hname, fNTelType, htitle );
    fDataTree->Branch( "ObsTime_s", &fObservingTime_s, "ObsTime_s/I" );
    fDataTree->Branch( "Offset_deg", &fOffset_deg, "Offset_deg/F" );
    
    sprintf( hname, "AverageTelescopeDistance[%d]", fNumTelTypes );
    sprintf( htitle, "AverageTelescopeDistance[%d]/F", fNumTelTypes );
    fDataTree->Branch( hname, fAverageTelescopeDistance, htitle );
    
    sprintf( hname, "MedianTelescopeDistance[%d]", fNumTelTypes );
    sprintf( htitle, "MedianTelescopeDistance[%d]/F", fNumTelTypes );
    fDataTree->Branch( hname, fMedianTelescopeDistance, htitle );
    
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
    if( argc != 7 )
    {
        cout << endl;
        cout << "./writeCTAWPPhysSensitivityTree <list of sub array> <site> <nscalings> <output file> <pointing direction> <ID>" << endl;
        cout << endl;
        cout << "nscalings = 5 (prod3 South) and = 1 (prod3 North)" << endl;
        cout << "NOTE SEVERAL HARDWIRED PARAMETERS" << endl;
        cout << "search for all files in hard wired PHYS directory" << endl;
        cout << endl;
        cout << "<pointing direction> (e.g. _180deg, _0deg, \"\" )" << endl;
        cout << "<ID>                 (e.g. 0, 1, 2, ...)" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    cout << endl;
    cout << "writeCTAWPPhysSensitivityTree"  << endl;
    cout << "-----------------------------" << endl;
    cout << endl;
    string fSubArrayList = argv[1];
    string fSite = argv[2];
    int    fNN_Min = 1;
    int    fNN_Max = atoi( argv[3] );
    string fOutputfile = argv[4];
    string fPointingDirection = argv[5];
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
    int    fAnalysisID = atoi( argv[6] );
    
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
    iObservingTimeVector.push_back( 18000 );
    iObservingTimeVector.push_back( 1800 );
    /*    iObservingTimeVector.push_back( 1800000 );
        iObservingTimeVector.push_back( 360000 );
        iObservingTimeVector.push_back( 180000 );
        iObservingTimeVector.push_back( 72000 );
        iObservingTimeVector.push_back( 36000 );
        iObservingTimeVector.push_back( 18000 );
        iObservingTimeVector.push_back( 7200 );
        iObservingTimeVector.push_back( 600 );
        iObservingTimeVector.push_back( 300 );
        iObservingTimeVector.push_back( 60 ); */
    //        iObservingTimeVector.push_back( 20 );
    
    
    /////////////////////////////////////////////////////////
    // fill events for complete parameter space
    
    // array loop
    for( unsigned int i = 0; i < fSubArray.size(); i++ )
    {
        // array scaling loop
        for( int a = fNN_Min; a <= fNN_Max; a++ )
        {
            stringstream iSubArray;
            iSubArray << fSubArray[i];
            if( fSite.find( "LaPalma" ) == string::npos )
            {
                iSubArray << "-" << a;
            }
            
            // observing time loop
            for( unsigned int t = 0; t < iObservingTimeVector.size(); t++ )
            {
                cout << endl << endl;
                cout << "now filling array layout " << iSubArray.str();
                cout << " (obs time " << iObservingTimeVector[t] << " s)" << endl;
                cout << "=========================================" << endl;
                cout << endl;
                fData->fillEvent( fSite, iSubArray.str(), iObservingTimeVector[t], i, a );
            }
        }
    }
    
    fData->terminate();
    
    return 0;
}

