/*! \file  trainTMVAforGammaHadronSeparation.cpp
    \brief  use TMVA methods for gamma/hadron separation

*/

#include "TChain.h"
#include "TCut.h"
#include "TEventList.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "VEvndispRunParameter.h"
#include "VTMVARunData.h"

using namespace std;

/*
 * check settings for number of training events;
 * if required reset event numbers
 *
 * note! assume name number of training and testing events
 */
string resetNumberOfTrainingEvents( string a, Long64_t n, bool iSignal, unsigned int iRequiredEvents )
{
    Long64_t n_fix = ( Long64_t )( n * 0.5 - 1 );
    stringstream r_a;
    if( n_fix > iRequiredEvents )
    {
        n_fix = iRequiredEvents;
    }
    if( iSignal )
    {
        r_a << a << ":nTrain_Signal=" << n_fix << ":nTest_Signal=" << n_fix;
    }
    else
    {
        r_a << a << ":nTrain_Background=" << n_fix << ":nTest_Background=" << n_fix;
    }
    return r_a.str();
}

/*
 * For cases when DispNImages < Images, return indexes for which ImgSel_list matches DispTelList_T
 */
vector< unsigned int> get_matching_indexes(unsigned int DispNImages, int NImages, UInt_t* DispTelList_T, UInt_t* ImgSel_list)
{
    vector< unsigned int> matching_indexes;

    for( unsigned int d = 0; d < DispNImages; d++ )
    {
        for( int i = 0; i < NImages; i++ )
        {
            if( DispTelList_T[d] == ImgSel_list[i] )
            {
                matching_indexes.push_back(i);
            }
        }
    }
    return matching_indexes;
}


bool get_largest_weight_disp_entries(
        unsigned int d_n, unsigned int DispNImages, int NImages, UInt_t* DispTelList_T, UInt_t *ImgSel_list,
        Float_t* DispWoff_T, Float_t* Disp_T, Float_t* cross, Float_t* loss, Float_t* size, Float_t *width, Float_t *length,
        Double_t *R, Double_t *ES, Double_t ErecS,
        Float_t *d_size, Float_t *d_width, Float_t *d_length, Float_t *d_cross,
        Float_t *d_disp, Float_t *d_disp_weigth, Float_t *d_loss, Float_t *d_erec_mc, Float_t *d_R )
{

    // Create a vector of pairs: {DispWoff_T value, original index}
    vector<pair<Float_t, unsigned int>> indexed_disp_woff(DispNImages);
    for (unsigned int i = 0; i < DispNImages; ++i)
    {
        indexed_disp_woff[i] = {DispWoff_T[i], i};
    }

    // Sort the vector in ascending order based on DispWoff_T value (smallest values first)
    sort(indexed_disp_woff.begin(), indexed_disp_woff.end(),
         [](const pair<Float_t, unsigned int>& a,
            const pair<Float_t, unsigned int>& b) {
            return a.first < b.first; // Sort ascending
           });

    // Rare events have DispNImages < Images
    // careful: some arrays are of [Images], others of [DispNImages]
    if( DispNImages < (unsigned int)NImages )
    {
        vector< unsigned int > matching_indexes = get_matching_indexes(DispNImages, NImages, DispTelList_T, ImgSel_list);
        for( unsigned int i = 0; i < matching_indexes.size(); i++ )
        {
            size[i] = size[matching_indexes[i]];
            width[i] = width[matching_indexes[i]];
            length[i] = length[matching_indexes[i]];
            loss[i] = loss[matching_indexes[i]];
            cross[i] = cross[matching_indexes[i]];
            R[i] = R[matching_indexes[i]];
            ES[i] = ES[matching_indexes[i]];
        }
    }
    // this case should never happen
    else if( DispNImages > (unsigned int)NImages )
    {
        return false;
    }

    // --- Step 2: Populate d_ arrays based on sorted indices ---

    unsigned int entries_to_copy = min(d_n, DispNImages);

    // Copy the largest 'entries_to_copy' values
    for (unsigned int i = 0; i < entries_to_copy; ++i) {
        unsigned int original_index = indexed_disp_woff[i].second;

        d_size[i] = log10(size[original_index]);
        d_width[i] = width[original_index];
        d_length[i] = length[original_index];
        d_cross[i] = cross[original_index];
        d_disp[i] = Disp_T[original_index];
        d_disp_weigth[i] = DispWoff_T[original_index];
        d_loss[i] = loss[original_index];
        d_R[i] = R[original_index];
        d_erec_mc[i] = (ES[original_index] - ErecS) / ES[original_index];
        if( d_erec_mc[i] > 20. ) d_erec_mc[i] = 20.;
        if( d_erec_mc[i] < -20. ) d_erec_mc[i] = -20.;
    }

    // --- Step 3: Handle DispNImages < d_n (duplicate values) ---
    // If DispNImages is less than d_n, duplicate the existing entries
    for (unsigned int i = entries_to_copy; i < d_n; ++i) {
        unsigned int source_index = (entries_to_copy > 0) ? (i % entries_to_copy) : 0;

        d_size[i] = d_size[source_index];
        d_width[i] = d_width[source_index];
        d_length[i] = d_length[source_index];
        d_cross[i] = d_cross[source_index];
        d_disp[i] = d_disp[source_index];
        d_disp_weigth[i] = d_disp_weigth[source_index];
        d_loss[i] = d_loss[source_index];
        d_R[i] = d_R[source_index];
        d_erec_mc[i] = d_erec_mc[source_index];
    }
    return true;
}

/*
 * prepare training / testing trees with reduced number of events
 *
 *   - apply pre-cuts here
 *   - copy only variables which are needed for TMVA into new tree
 *   - delete full trees (IMPORTANT)
 *
 */
TTree* prepareSelectedEventsTree( VTMVARunData* iRun, TCut iCut,
                                  bool iSignal, Long64_t iResetEventNumbers )
{
    if(!iRun )
    {
        return 0;
    }
    vector< TChain* > iTreeVector;
    string iDataTree_reducedName;
    if( iSignal )
    {
        cout << "Preparing reduced signal trees" << endl;
        iTreeVector = iRun->fSignalTree;
        iDataTree_reducedName = "data_signal";
    }
    else
    {
        cout << "Preparing reduced background trees" << endl;
        iTreeVector = iRun->fBackgroundTree;
        iDataTree_reducedName = "data_background";
    }
    // reduced tree (and name)
    TTree* iDataTree_reduced = 0;
    // list of variables copied.
    // must include at least the variables used for the training
    Double_t MSCW = 0.;
    Double_t MSCL = 0.;
    Double_t ErecS = 0.;
    Double_t EChi2S = 0.;
    Double_t Xcore = 0.;
    Double_t Ycore = 0.;
    Double_t Xoff_derot = 0.;
    Double_t Yoff_derot = 0.;
    Float_t Xoff_intersect = 0.;
    Float_t Yoff_intersect = 0.;
    Int_t NImages = 0;
    UInt_t DispNImages = 0;
    // fixed max number of telescope types
    UInt_t NImages_Ttype[20];
    // fixed max number of telescopes
    UInt_t ImgSel_list[200];
    UInt_t DispTelList_T[200];
    for( unsigned int i = 0; i < 20; i++ )
    {
        NImages_Ttype[i] = 0;
    }
    Float_t EmissionHeight = 0.;
    Float_t EmissionHeightChi2 = 0.;
    Double_t SizeSecondMax = 0.;
    Double_t DispDiff = 0.;
    Float_t DispAbsSumWeigth = 0.;
    Double_t MCe0 = 0.;
    Double_t MCxoff = 0.;
    Double_t MCyoff = 0.;
    Float_t disp_mc_error = 0.;
    Float_t intersect_mc_error = 0.;
    Float_t size[200];
    Float_t width[200];
    Float_t length[200];
    Float_t cross[200];
    Float_t disp[200];
    Float_t disp_weight[200];
    Float_t loss[200];
    Double_t R[200];
    Double_t ES[200];

    // add d_n most important variables to reduced tree
    const int d_n = 4;
    Float_t d_size[d_n];
    Float_t d_width[d_n];
    Float_t d_length[d_n];
    Float_t d_cross[d_n];
    Float_t d_disp[d_n];
    Float_t d_disp_weigth[d_n];
    Float_t d_loss[d_n];
    Float_t d_R[d_n];
    Float_t d_erec_es[d_n];

    char htemp[200];

    iDataTree_reduced = new TTree( iDataTree_reducedName.c_str(), iDataTree_reducedName.c_str() );
    iDataTree_reduced->Branch( "MSCW", &MSCW, "MSCW/D" );
    iDataTree_reduced->Branch( "MSCL", &MSCL, "MSCL/D" );
    iDataTree_reduced->Branch( "ErecS", &ErecS, "ErecS/D" );
    iDataTree_reduced->Branch( "EChi2S", &EChi2S, "EChi2S/D" );
    iDataTree_reduced->Branch( "Xcore", &Xcore, "Xcore/D" );
    iDataTree_reduced->Branch( "Ycore", &Ycore, "Ycore/D" );
    iDataTree_reduced->Branch( "NImages", &NImages, "NImages/I" );
    iDataTree_reduced->Branch( "NImages_Ttype", NImages_Ttype, "NImages_Ttype[20]/i" );
    iDataTree_reduced->Branch( "EmissionHeight", &EmissionHeight, "EmissionHeight/F" );
    iDataTree_reduced->Branch( "EmissionHeightChi2", &EmissionHeightChi2, "EmissionHeightChi2/F" );
    iDataTree_reduced->Branch( "SizeSecondMax", &SizeSecondMax, "SizeSecondMax/D" );
    iDataTree_reduced->Branch( "DispDiff", &DispDiff, "DispDiff/D" );
    iDataTree_reduced->Branch( "DispAbsSumWeigth", &DispAbsSumWeigth, "DispAbsSumWeigth/F" );
    iDataTree_reduced->Branch( "MCe0", &MCe0, "MCe0/D" );
    iDataTree_reduced->Branch( "disp_mc_error", &disp_mc_error, "disp_mc_error/F" );
    iDataTree_reduced->Branch( "intersect_mc_error", &intersect_mc_error, "intersect_mc_error/F" );
    sprintf( htemp, "d_size[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_size", d_size, htemp );
    sprintf( htemp, "d_width[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_width", d_width, htemp );
    sprintf( htemp, "d_length[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_length", d_length, htemp );
    sprintf( htemp, "d_cross[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_cross", d_cross, htemp );
    sprintf( htemp, "d_disp[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_disp", d_disp, htemp );
    sprintf( htemp, "d_disp_weigth[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_disp_weigth", d_disp_weigth, htemp );
    sprintf( htemp, "d_loss[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_loss", d_loss, htemp );
    sprintf( htemp, "d_R[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_R", d_R, htemp );
    sprintf( htemp, "d_erec_es[%d]/F", d_n );
    iDataTree_reduced->Branch( "d_erec_es", d_erec_es, htemp );

    Long64_t n = 0;

    for( unsigned  int i = 0; i < iTreeVector.size(); i++ )
    {
        if( iTreeVector[i] )
        {
            iTreeVector[i]->SetBranchAddress( "MSCW", &MSCW );
            iTreeVector[i]->SetBranchAddress( "MSCL", &MSCL );
            iTreeVector[i]->SetBranchAddress( "ErecS", &ErecS );
            iTreeVector[i]->SetBranchAddress( "EChi2S", &EChi2S );
            iTreeVector[i]->SetBranchAddress( "Xcore", &Xcore );
            iTreeVector[i]->SetBranchAddress( "Ycore", &Ycore );
            iTreeVector[i]->SetBranchAddress( "Xoff_derot", &Xoff_derot );
            iTreeVector[i]->SetBranchAddress( "Yoff_derot", &Yoff_derot );
            iTreeVector[i]->SetBranchAddress( "Xoff_intersect", &Xoff_intersect );
            iTreeVector[i]->SetBranchAddress( "Yoff_intersect", &Yoff_intersect );
            iTreeVector[i]->SetBranchAddress( "NImages", &NImages );
            iTreeVector[i]->SetBranchAddress( "ImgSel_list", ImgSel_list );
            iTreeVector[i]->SetBranchAddress( "NImages_Ttype", NImages_Ttype );
            iTreeVector[i]->SetBranchAddress( "EmissionHeight", &EmissionHeight );
            iTreeVector[i]->SetBranchAddress( "EmissionHeightChi2", &EmissionHeightChi2 );
            iTreeVector[i]->SetBranchAddress( "SizeSecondMax", &SizeSecondMax );
            iTreeVector[i]->SetBranchAddress( "DispDiff", &DispDiff );
            iTreeVector[i]->SetBranchAddress( "DispAbsSumWeigth", &DispAbsSumWeigth );
            iTreeVector[i]->SetBranchAddress( "DispNImages", &DispNImages );
            iTreeVector[i]->SetBranchAddress( "size", size );
            iTreeVector[i]->SetBranchAddress( "width", width );
            iTreeVector[i]->SetBranchAddress( "length", length );
            iTreeVector[i]->SetBranchAddress( "cross", cross );
            iTreeVector[i]->SetBranchAddress( "Disp_T", disp );
            iTreeVector[i]->SetBranchAddress( "DispWoff_T", disp_weight );
            iTreeVector[i]->SetBranchAddress( "DispTelList_T", DispTelList_T );
            iTreeVector[i]->SetBranchAddress( "loss", loss );
            iTreeVector[i]->SetBranchAddress( "R", R );
            iTreeVector[i]->SetBranchAddress( "ES", ES );
            if( iTreeVector[i]->GetBranchStatus( "MCe0" ) )
            {
                iTreeVector[i]->SetBranchAddress( "MCe0", &MCe0 );
                iTreeVector[i]->SetBranchAddress( "MCxoff", &MCxoff );
                iTreeVector[i]->SetBranchAddress( "MCyoff", &MCyoff );
            }
            if(!iDataTree_reduced )
            {
                cout << "Error preparing reduced tree" << endl;
                cout << "exiting..." << endl;
                exit( EXIT_FAILURE );
            }
            iTreeVector[i]->Draw( ">>elist", iCut, "entrylist" );
            TEntryList* elist = ( TEntryList* )gDirectory->Get( "elist" );
            if( elist )
            {
                for( Long64_t el = 0; el < elist->GetN(); el++ )
                {
                    Long64_t treeEntry = elist->GetEntry( el );
                    iTreeVector[i]->GetEntry( treeEntry );

                    // error in angular reconstruction
                    disp_mc_error = sqrt( (Xoff_derot-MCxoff)*(Xoff_derot-MCxoff) + (Yoff_derot-MCyoff)*(Yoff_derot-MCyoff) );
                    intersect_mc_error = sqrt( (Xoff_intersect-MCxoff)*(Xoff_intersect-MCxoff) + (Yoff_intersect-Yoff_intersect)*(Yoff_intersect-Yoff_intersect) );

                    // d_n largest weight entries
                    bool good_event = get_largest_weight_disp_entries(
                           d_n, DispNImages, NImages, DispTelList_T, ImgSel_list,
                           disp_weight, disp, cross, loss, size, width, length, R, ES, ErecS,
                           d_size, d_width, d_length, d_cross, d_disp, d_disp_weigth, d_loss, d_erec_es, d_R
                    );
                    if( !good_event )
                    {
                        continue;
                    }

                    iDataTree_reduced->Fill();
                    n++;
                }
            }
            // remove this tree
            if( iSignal )
            {
                iRun->fSignalTree[i]->Delete();
                iRun->fSignalTree[i] = 0;
            }
            else
            {
                iRun->fBackgroundTree[i]->Delete();
                iRun->fBackgroundTree[i] = 0;
            }
            // factor of 2: here for training and testing events
            if( iResetEventNumbers > 0 && n > iResetEventNumbers * 2 )
            {
                cout << "\t reached required ";
                if( iSignal )
                {
                    cout << "signal";
                }
                else
                {
                    cout << "background";
                }
                cout << " event numbers ";
                cout << "(" << iResetEventNumbers << ")";
                cout << " after " << i + 1 << " tree(s)" << endl;
                cout << "\t applied cut: " << iCut << endl;
                break;
            }
        }
    }
    if( iSignal && iDataTree_reduced )
    {
        cout << "\t Reduced signal tree entries: " << iDataTree_reduced->GetEntries() << endl;
    }
    else if( iDataTree_reduced )
    {
        cout << "\t Reduced background tree entries: " << iDataTree_reduced->GetEntries() << endl;
    }
    else
    {
        cout << "Error in reducing data trees (missing tree)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    // cleanup all remaining trees
    for( unsigned int i = 0; i < iTreeVector.size(); i++ )
    {
        if( iSignal && iRun->fSignalTree[i] )
        {
            iRun->fSignalTree[i]->Delete();
            iRun->fSignalTree[i] = 0;
        }
        else if(!iSignal && iRun->fBackgroundTree[i] )
        {
            iRun->fBackgroundTree[i]->Delete();
            iRun->fBackgroundTree[i] = 0;
        }
    }
    return iDataTree_reduced;
}


bool train( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin, string iRunMode )
{
    // sanity checks
    if(!iRun )
    {
        return false;
    }
    if( iRun->fEnergyCutData.size() <= iEnergyBin || iRun->fOutputFile.size() <= iEnergyBin )
    {
        cout << "error during training: energy bin out of range " << iEnergyBin << endl;
        return false;
    }
    if( iRun->fZenithCutData.size() < iZenithBin || iRun->fOutputFile[0].size() < iZenithBin )
    {
        cout << "error during training: zenith bin out of range " << iZenithBin << endl;
        return false;
    }
    // quality cuts before training
    TCut iCutSignal = iRun->fQualityCuts && iRun->fQualityCutsSignal
                      && iRun->fMCxyoffCut && iRun->fAzimuthCut &&
                      iRun->fMultiplicityCuts &&
                      iRun->fEnergyCutData[iEnergyBin]->fEnergyCut
                      && iRun->fZenithCutData[iZenithBin]->fZenithCut;

    TCut iCutBck = iRun->fQualityCuts && iRun->fQualityCutsBkg
                   && iRun->fAzimuthCut
                   && iRun->fMultiplicityCuts
                   && iRun->fEnergyCutData[iEnergyBin]->fEnergyCut
                   && iRun->fZenithCutData[iZenithBin]->fZenithCut;

    if(!iRun->fMCxyoffCutSignalOnly )
    {
        iCutBck = iCutBck && iRun->fMCxyoffCut;
    }

    // adding training variables
    if( iRun->fTrainingVariable.size() != iRun->fTrainingVariableType.size() )
    {
        cout << "train: error: training-variable vectors have different size" << endl;
        return false;
    }

    // prepare trees for training and testing with selected events only
    // this step is necessary to minimise the memory impact for the BDT
    // training
    if( iRunMode == "WriteTrainingEvents" )
    {
        prepareSelectedEventsTree( iRun,
                              iCutSignal, true,
                              iRun->fResetNumberOfTrainingEvents );
        prepareSelectedEventsTree( iRun,
                                  iCutBck, false,
                                  iRun->fResetNumberOfTrainingEvents );

        if( iRun->getTLRunParameter() )
        {
            iRun->getTLRunParameter()->Write();
        }
        cout << "Writing reduced event lists for training: ";
        cout << gDirectory->GetName() << endl;
        exit( EXIT_SUCCESS );
    }

    ////////////////////////////////////////////////////////////////
    // Prepare TMVA instances
    TMVA::Tools::Instance();
    gSystem->mkdir( iRun->fOutputDirectoryName.c_str() );
    TString iOutputDirectory( iRun->fOutputDirectoryName.c_str() );
    gSystem->ExpandPathName( iOutputDirectory );
    ( TMVA::gConfig().GetIONames() ).fWeightFileDir = iOutputDirectory;

    //////////////////////////////////////////
    // defining training class
    string mva_options = "V:!DrawProgressBar";
    if( iRunMode == "TrainReconstructionQuality" )
    {
        mva_options = "V:!DrawProgressBar:!Color:!Silent:AnalysisType=Regression:VerboseLevel=Debug:Correlations=True";
    }
    TMVA::Factory* factory = new TMVA::Factory( iRun->fOutputFile[iEnergyBin][iZenithBin]->GetTitle(),
        iRun->fOutputFile[iEnergyBin][iZenithBin],
        mva_options.c_str() );
    TMVA::DataLoader* dataloader = new TMVA::DataLoader( "" );

    // training preparation
    cout << "Reading training / testing trees from ";
    cout << iRun->fSelectedEventTreeName << endl;
    TFile* iF = new TFile( iRun->fSelectedEventTreeName.c_str() );
    if( iF->IsZombie() )
    {
        cout << "Error open file with pre-selected events: ";
        cout << iRun->fSelectedEventTreeName << endl;
        exit( EXIT_FAILURE );
    }
    TTree *iSignalTree_reduced = ( TTree* )iF->Get( "data_signal" );
    TTree *iBackgroundTree_reduced = ( TTree* )iF->Get( "data_background" );
    if(!iSignalTree_reduced || !iBackgroundTree_reduced )
    {
        cout << "Error: failed preparing traing / testing trees" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_SUCCESS );
    }
    // check for number of training and background events
    cout << "Reading number of events before / after cuts: " << endl;
    Long64_t nEventsSignal = iSignalTree_reduced->GetEntries();
    cout << "\t total number of signal events before cuts: ";
    cout << nEventsSignal << endl;
    iSignalTree_reduced->Draw( ">>elist", iRun->fEnergyCutData[iEnergyBin]->fEnergyCut && iRun->fMultiplicityCuts, "entrylist" );
    TEntryList* elist = ( TEntryList* )gDirectory->Get( "elist" );
    if( elist )
    {
        nEventsSignal = elist->GetN();
        cout << "\t total number of signal events after cuts:  ";
        cout << elist->GetN();
        cout << " (required are " << iRun->fMinSignalEvents << ")" << endl;
    }
    Long64_t nEventsBck = iBackgroundTree_reduced->GetEntries();
    cout << "\t total number of background events before cuts: ";
    cout << nEventsBck << endl;
    iBackgroundTree_reduced->Draw( ">>elist", iRun->fEnergyCutData[iEnergyBin]->fEnergyCut && iRun->fMultiplicityCuts, "entrylist" );
    elist = ( TEntryList* )gDirectory->Get( "elist" );
    if( elist )
    {
        nEventsBck = elist->GetN();
        cout << "\t total number of background events after cuts:  ";
        cout << elist->GetN();
        cout << " (required are " << iRun->fMinSignalEvents << ")" << endl;
        cout << " (required are " << iRun->fMinBackgroundEvents << ")" << endl;
    }
    if( nEventsSignal < iRun->fMinSignalEvents || nEventsBck < iRun->fMinBackgroundEvents )
    {
        cout << "Error: not enough training events" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_SUCCESS );
    }
    // check if this number if consistent with requested number of training events
    // required at least 10 events
    if( iRun->fResetNumberOfTrainingEvents > 0 )
    {
        cout << "Checking number of training events: " << endl;
        iRun->fPrepareTrainingOptions = resetNumberOfTrainingEvents( iRun->fPrepareTrainingOptions,
                                        nEventsSignal, true, iRun->fResetNumberOfTrainingEvents );
        iRun->fPrepareTrainingOptions = resetNumberOfTrainingEvents( iRun->fPrepareTrainingOptions,
                                        nEventsBck, false, iRun->fResetNumberOfTrainingEvents );
        cout << "\t Updated training options: " <<  iRun->fPrepareTrainingOptions << endl;
    }

    ////////////////////////////
    // train gamma/hadron separation
    if( iRunMode == "TrainGammaHadronSeparation" )
    {
        dataloader->AddSignalTree( iSignalTree_reduced, iRun->fSignalWeight );
        dataloader->AddBackgroundTree( iBackgroundTree_reduced, iRun->fBackgroundWeight );
    }
    ////////////////////////////
    // train for angular reconstruction method
    else if( iRunMode == "TrainAngularReconstructionMethod" )
    {
        TCut signalCut = "intersect_mc_error < 1.e5 && (disp_mc_error > intersect_mc_error)";
        TCut backgrCut = "disp_mc_error < intersect_mc_error";

        iRun->fOutputFile[iEnergyBin][iZenithBin]->cd();

        TTree* sigTree = iSignalTree_reduced->CopyTree(signalCut);
        TTree* bkgTree = iSignalTree_reduced->CopyTree(backgrCut);
        sigTree->SetName("data_signal");
        bkgTree->SetName("data_background");

        dataloader->AddSignalTree(sigTree, iRun->fSignalWeight);
        dataloader->AddBackgroundTree(bkgTree, iRun->fBackgroundWeight);
    }
    ////////////////////////////
    // train reconstruction quality
    else if( iRunMode == "TrainReconstructionQuality" )
    {
        dataloader->AddRegressionTree( iSignalTree_reduced, iRun->fSignalWeight );
        dataloader->AddTarget( iRun->fReconstructionQualityTarget.c_str(), 'F' );
    }

    // loop over all training variables and add them to TMVA
    for( unsigned int i = 0; i < iRun->fTrainingVariable.size(); i++ )
    {
        if( iRun->fTrainingVariable[i].rfind("d_", 0 ) == 0 )
        {
            dataloader->AddVariablesArray(iRun->fTrainingVariable[i].c_str(), 4 );
        }
        else
        {
            dataloader->AddVariable( iRun->fTrainingVariable[i].c_str(), iRun->fTrainingVariableType[i] );
        }
    }
    // adding spectator variables
    for( unsigned int i = 0; i < iRun->fSpectatorVariable.size(); i++ )
    {
        dataloader->AddSpectator( iRun->fSpectatorVariable[i].c_str() );
    }

    //////////////////////////////////////////
    // prepare training events
    // (cuts are already applied at an earlier stage)
    cout << "Preparing training and test tree" << endl;
    // cuts after pre-selection
    TCut iCutSignal_post = iRun->fEnergyCutData[iEnergyBin]->fEnergyCut
                           && iRun->fMultiplicityCuts;
    TCut iCutBck_post = iRun->fEnergyCutData[iEnergyBin]->fEnergyCut
                        && iRun->fMultiplicityCuts;

    if( iRunMode == "TrainReconstructionQuality" )
    {
        dataloader->PrepareTrainingAndTestTree( iCutSignal_post, iRun->fPrepareTrainingOptions );
    }
    else
    {
        dataloader->PrepareTrainingAndTestTree( iCutSignal_post,
                                                iCutBck_post,
                                                iRun->fPrepareTrainingOptions );
    }

    //////////////////////////////////////////
    // book all methods
    char htitle[6000];

    for( unsigned int i = 0; i < iRun->fMVAMethod.size(); i++ )
    {
        TMVA::Types::EMVA i_tmva_type = TMVA::Types::kBDT;
        if( iRun->fMVAMethod[i] == "BDT" )
        {
            i_tmva_type = TMVA::Types::kBDT;
        }
        else if( iRun->fMVAMethod[i] == "MLP" )
        {
            i_tmva_type = TMVA::Types::kMLP;
        }

        //////////////////////////
        if( iRunMode == "TrainGammaHadronSeparation" )
        {
            sprintf( htitle, "%s_%u", iRun->fMVAMethod[i].c_str(), i );
        }
        else if( iRunMode == "TrainAngularReconstructionMethod" )
        {
            sprintf( htitle, "%s_RecMethod_%u", iRun->fMVAMethod[i].c_str(), i );
        }
        else if( iRunMode == "TrainReconstructionQuality" )
        {
            sprintf( htitle, "%s_RecQuality_%u", iRun->fMVAMethod[i].c_str(), i );
        }
        if( i < iRun->fMVAMethod_Options.size() )
        {
            cout << "Booking method " << htitle << endl;
            factory->BookMethod( dataloader, i_tmva_type, htitle, iRun->fMVAMethod_Options[i].c_str() );
        }
        else
        {
            factory->BookMethod( dataloader, i_tmva_type, htitle );
        }
    }


    //////////////////////////////////////////
    // start training

    factory->TrainAllMethods();

    //////////////////////////////////////////
    // evaluate results

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

    factory->Delete();

    return true;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

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
    cout << "trainTMVAforGammaHadronSeparation " << VGlobalRunParameter::getEVNDISP_VERSION() << endl;
    cout << "----------------------------------------" << endl;
    if( argc != 2 && argc != 3 )
    {
        cout << endl;
        cout << "./trainTMVAforGammaHadronSeparation <configuration file> [WRITETRAININGEVENTS]" << endl;
        cout << endl;
        cout << "  (an example for a configuration file can be found in " << endl;
        cout << "   $CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDT.runparameter )" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    cout << endl;

    //////////////////////////////////////
    // data object
    VTMVARunData* fData = new VTMVARunData();
    fData->fName = "OO";

    //////////////////////////////////////
    // read run parameters from configuration file
    if(!fData->readConfigurationFile( argv[1] ) )
    {
        cout << "error opening or reading run parameter file (";
        cout << argv[1];
        cout << ")" << endl;
        exit( EXIT_FAILURE );
    }
    // randomize list of input files
    fData->shuffleFileVectors();
    fData->print();

    //////////////////////////////////////
    // read and prepare data files
    if(!fData->openDataFiles( false ) )
    {
        cout << "error opening data files" << endl;
        exit( EXIT_FAILURE );
    }

    //////////////////////////////////////
    // train MVA
    // (one training per energy and zenith bin)
    cout << "Run mode " << fData->fRunMode << endl;
    cout << "Number of energy bins: " << fData->fEnergyCutData.size();
    cout << ", number of zenith bins: " << fData->fZenithCutData.size();
    cout << endl;
    cout << "================================" << endl << endl;
    for( unsigned int i = 0; i < fData->fEnergyCutData.size(); i++ )
    {
        for( unsigned int j = 0; j < fData->fZenithCutData.size(); j++ )
        {
            if( fData->fEnergyCutData[i]->fEnergyCut && fData->fZenithCutData[j]->fZenithCut )
            {
                cout << "Training energy bin " << fData->fEnergyCutData[i]->fEnergyCut;
                cout << " zenith bin " << fData->fZenithCutData[j]->fZenithCut << endl;
                cout << "===================================================================================" << endl;
                cout << endl;
            }
            // training
            if( !train( fData, i, j, fData->fRunMode) )
            {
                cout << "Error during training...exiting" << endl;
                exit( EXIT_FAILURE );
            }
            stringstream iTempS;
            stringstream iTempS2;
            if( fData->fEnergyCutData.size() > 1 && fData->fZenithCutData.size() > 1 )
            {
                iTempS << fData->fOutputDirectoryName << "/" << fData->fOutputFileName << "_" << i << "_" << j << ".bin.root";
                iTempS2 << "/" << fData->fOutputFileName << "_" << i << "_" << j << ".root";
            }
            else if( fData->fEnergyCutData.size() > 1 && fData->fZenithCutData.size() <= 1 )
            {
                iTempS << fData->fOutputDirectoryName << "/" << fData->fOutputFileName << "_" << i << ".bin.root";
                iTempS2 << "/" << fData->fOutputFileName << "_" << i << ".root";
            }
            else if( fData->fZenithCutData.size() > 1 &&  fData->fEnergyCutData.size() <= 1 )
            {
                iTempS << fData->fOutputDirectoryName << "/" << fData->fOutputFileName << "_0_" << j << ".bin.root";
                iTempS2 << "/" << fData->fOutputFileName << "_0_" << j << ".root";
            }
            else
            {
                iTempS << fData->fOutputDirectoryName << "/" << fData->fOutputFileName << ".bin.root";
                iTempS2 << fData->fOutputFileName << ".root";
            }

            // prepare a short root file with the necessary values only
            // write energy & zenith cuts, plus signal and background efficiencies
            TFile* root_file = fData->fOutputFile[i][j];
            if(!root_file )
            {
                cout << "Error finding tvma root file " << endl;
                continue;
            }
            TFile* short_root_file = TFile::Open( iTempS.str().c_str(), "RECREATE" );
            if(!short_root_file->IsZombie() )
            {
                VTMVARunDataEnergyCut* fDataEnergyCut = ( VTMVARunDataEnergyCut* )root_file->Get( "fDataEnergyCut" );
                VTMVARunDataZenithCut* fDataZenithCut = ( VTMVARunDataZenithCut* )root_file->Get( "fDataZenithCut" );
                TH1D* MVA_effS = 0;
                TH1D* MVA_effB = 0;

                char hname[200];
                for( unsigned int d = 0; d < fData->fMVAMethod.size(); d++ )
                {
                    // naming of directories is different for different TMVA versions
                    sprintf( hname, "Method_%s_%u/%s_%u/MVA_%s_%u_effS",
                             fData->fMVAMethod[d].c_str(), d,
                             fData->fMVAMethod[d].c_str(), d,
                             fData->fMVAMethod[d].c_str(), d );
                    if(( TH1D* )root_file->Get( hname ) )
                    {
                        MVA_effS = ( TH1D* )root_file->Get( hname );
                        sprintf( hname, "Method_%s_%u/%s_%u/MVA_%s_%u_effB",
                                 fData->fMVAMethod[d].c_str(), d,
                                 fData->fMVAMethod[d].c_str(), d,
                                 fData->fMVAMethod[d].c_str(), d );
                        MVA_effB = ( TH1D* )root_file->Get( hname );
                    }
                    else
                    {
                        sprintf( hname, "Method_%s/%s_%u/MVA_%s_%u_effS",
                                 fData->fMVAMethod[d].c_str(),
                                 fData->fMVAMethod[d].c_str(), d,
                                 fData->fMVAMethod[d].c_str(), d );
                        MVA_effS = ( TH1D* )root_file->Get( hname );
                        sprintf( hname, "Method_%s/%s_%u/MVA_%s_%u_effB",
                                 fData->fMVAMethod[d].c_str(),
                                 fData->fMVAMethod[d].c_str(), d,
                                 fData->fMVAMethod[d].c_str(), d );
                        MVA_effB = ( TH1D* )root_file->Get( hname );
                    }

                    if( fDataEnergyCut )
                    {
                        fDataEnergyCut->Write();
                    }
                    if( fDataZenithCut )
                    {
                        fDataZenithCut->Write();
                    }
                    sprintf( hname, "Method_%s_%u", fData->fMVAMethod[d].c_str(), d );
                    TDirectory* Method_MVA = short_root_file->mkdir( hname );
                    Method_MVA->cd();
                    sprintf( hname, "%s_%u", fData->fMVAMethod[d].c_str(), d );
                    TDirectory* MVA = Method_MVA->mkdir( hname );
                    MVA->cd();
                    if( MVA_effS )
                    {
                        MVA_effS->Write();
                    }
                    if( MVA_effB )
                    {
                        MVA_effB->Write();
                    }
                    short_root_file->GetList();
                    short_root_file->Write();
                    short_root_file->cd();
                }
                short_root_file->Close();
            }
            else
            {
                cout << "Error: could not create file with energy cuts " << iTempS.str().c_str() << endl;
            }
            // copy complete TMVA output root-file to another directory
            string iOutputFileName( fData->fOutputDirectoryName + "/" + iTempS2.str() );
            string iOutputFileNameCompleteSubDir( "complete_BDTroot" );
            string iOutputFileNameCompleteDir( fData->fOutputDirectoryName + "/" + iOutputFileNameCompleteSubDir + "/" );
            gSystem->mkdir( iOutputFileNameCompleteDir.c_str() );
            string iOutputFileNameComplete( iOutputFileNameCompleteDir + iTempS2.str() );
            rename( iOutputFileName.c_str(), iOutputFileNameComplete.c_str() );
            cout << "Complete TMVA output root-file moved to: " << iOutputFileNameComplete << endl;

            // rename .bin.root file to .root-file
            string iFinalRootFileName( iTempS.str() );
            string iBinRootString( ".bin.root" );
            iFinalRootFileName.replace( iFinalRootFileName.find( iBinRootString ), iBinRootString.length(), ".root" );
            rename( iTempS.str().c_str(), iFinalRootFileName.c_str() );
        }
    }
    return 0;
}
