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

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "VEvndispRunParameter.h"
#include "VTMVARunData.h"

using namespace std;

bool train( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin, bool iGammaHadronSeparation, string fRunOption );
bool trainGammaHadronSeparation( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin, string fRunOption );
bool trainReconstructionQuality( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin, string fRunOption );

/*
 * check settings for number of training events;
 * if required reset event numbers
 *
 * note! assume name number of training and testing events
 */
string resetNumberOfTrainingEvents( string a, Long64_t n, bool iSignal, unsigned int iRequiredEvents )
{
   Long64_t n_fix = (Long64_t)(n*0.5-1);
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
    if( !iRun )
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
    TTree *iDataTree_reduced = 0;
    // list of variables copied.
    // must be at least the variables used for the training
    Double_t MSCW = 0.;
    Double_t MSCL = 0.;
    Double_t ErecS = 0.;
    Double_t EChi2S = 0.;
    Double_t dES = 0.;
    // fixed max number of telescope types
    UInt_t NImages_Ttype[20];
    for( unsigned int i = 0; i < 20; i++ ) NImages_Ttype[i] = 0;
    Float_t EmissionHeight = 0.;
    Float_t EmissionHeightChi2 = 0.;
    Double_t SizeSecondMax = 0.;
    Double_t DispDiff = 0.;
    Double_t MCe0 = 0.;
    iDataTree_reduced = new TTree( iDataTree_reducedName.c_str(), iDataTree_reducedName.c_str() );
    iDataTree_reduced->Branch( "MSCW", &MSCW, "MSCW/D" );
    iDataTree_reduced->Branch( "MSCL", &MSCL, "MSCL/D" );
    iDataTree_reduced->Branch( "ErecS", &ErecS, "ErecS/D" );
    iDataTree_reduced->Branch( "EChi2S", &EChi2S, "EChi2S/D" );
    iDataTree_reduced->Branch( "dES", &dES, "dES/D" );
    iDataTree_reduced->Branch( "NImages_Ttype", NImages_Ttype, "NImages_Ttype[20]/i" );
    iDataTree_reduced->Branch( "EmissionHeight", &EmissionHeight, "EmissionHeight/F" );
    iDataTree_reduced->Branch( "EmissionHeightChi2", &EmissionHeightChi2, "EmissionHeightChi2/F" );
    iDataTree_reduced->Branch( "SizeSecondMax", &SizeSecondMax, "SizeSecondMax/D" );
    iDataTree_reduced->Branch( "DispDiff", &DispDiff, "DispDiff/D" );
    iDataTree_reduced->Branch( "MCe0", &MCe0, "MCe0/D" );

    Long64_t n = 0;
    
    for( unsigned  int i = 0; i < iTreeVector.size(); i++ )
    {
        if( iTreeVector[i] )
        {
             iTreeVector[i]->SetBranchAddress( "MSCW", &MSCW );
             iTreeVector[i]->SetBranchAddress( "MSCL", &MSCL );
             iTreeVector[i]->SetBranchAddress( "ErecS", &ErecS );
             iTreeVector[i]->SetBranchAddress( "EChi2S", &EChi2S );
             iTreeVector[i]->SetBranchAddress( "dES", &dES );
             iTreeVector[i]->SetBranchAddress( "NImages_Ttype", NImages_Ttype );
             iTreeVector[i]->SetBranchAddress( "EmissionHeight", &EmissionHeight );
             iTreeVector[i]->SetBranchAddress( "EmissionHeightChi2", &EmissionHeightChi2 );
             iTreeVector[i]->SetBranchAddress( "SizeSecondMax", &SizeSecondMax );
             iTreeVector[i]->SetBranchAddress( "DispDiff", &DispDiff );
             iTreeVector[i]->SetBranchAddress( "MCe0", &MCe0 );
             if( !iDataTree_reduced )
             {
                 cout << "Error preparing reduced tree" << endl;
                 cout << "exiting..." << endl;
                 exit( EXIT_FAILURE );
             }
             iTreeVector[i]->Draw( ">>elist", iCut, "entrylist" );
             TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
             if( elist )
             {
                  for( Long64_t el = 0; el < elist->GetN(); el++)
                  {
                       Long64_t treeEntry = elist->GetEntry( el );
                       iTreeVector[i]->GetEntry( treeEntry );
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
                 cout << " after " << i+1 << " tree(s)" << endl;
                 cout << "\t applied cut: " << iCut << endl;
                break;
             }
        }
    }
    if( iSignal && iDataTree_reduced )
    {
        cout << "\t Reduced signal tree entries: " << iDataTree_reduced->GetEntries() << endl;
        cout << "Removing signal tree(s)" << endl;
    }
    else if( iDataTree_reduced )
    {
        cout << "\t Reduced background tree entries: " << iDataTree_reduced->GetEntries() << endl;
        cout << "Removing background tree(s)" << endl;
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
        else if( !iSignal && iRun->fBackgroundTree[i] )
        {
            iRun->fBackgroundTree[i]->Delete();
            iRun->fBackgroundTree[i] = 0;
        }
    }

    return iDataTree_reduced;
}

/*
   check if a training variable is constant

   (constant variables are removed from set of training variables)

   return values:

   -1:  value is variable
   0-N: value is constant (array identifier)

*/
double checkIfVariableIsConstant( VTMVARunData* iRun, TCut iCut, string iVariable, bool iSignal, bool iSplitBlock )
{
    char hname[2000];
    TH1D* h = 0;
    TH1I* hI = 0;
    cout << "initializing TMVA variables: checking";
    if( iSignal )
    {
        cout << " signal";
    }
    else
    {
        cout << " background";
    }
    cout << " variable " << iVariable << " for consistency " << endl;
    vector< TChain* > iTreeVector;
    if( iSignal )
    {
        iTreeVector = iRun->fSignalTree;
    }
    else
    {
        iTreeVector = iRun->fBackgroundTree;
    }
    
    for( unsigned  int i = 0; i < iTreeVector.size(); i++ )
    {
        h = 0;
        hI = 0;
        if( iTreeVector[i] )
        {
            Long64_t iNEntriesBlock = 0;
            if( iSplitBlock )
            {
                iNEntriesBlock = iTreeVector[i]->GetEntries() / 2;
            }
            else
            {
                iNEntriesBlock = iTreeVector[i]->GetEntries();
            }
            // fill a histogram with the variable to be checked
            sprintf( hname, "hXX_%u", i );
            if( iVariable.find( "NImages_Ttype" ) != string::npos )
            {
                hI = new TH1I( hname, "", 500, 0., 500. );
                iTreeVector[i]->Project( hI->GetName(), iVariable.c_str(), iCut, "", iNEntriesBlock );
                if( hI->GetRMS() > 1.e-5 )
                {
                    cout << "\t variable " << iVariable << " ok, RMS: " << hI->GetRMS() << ", tree: " << i;
                    cout << ", entries " << hI->GetEntries();
                    cout << endl;
                    hI->Delete();
                    return -9999.;
                }
            }
            else
            {
                h = new TH1D( hname, "", 100, -1.e5, 1.e5 );
                iTreeVector[i]->Project( h->GetName(), iVariable.c_str(), iCut, "", iNEntriesBlock );
                if( h->GetRMS() > 1.e-5 )
                {
                    cout << "\t variable " << iVariable << " ok, RMS: " << h->GetRMS() << ", tree: " << i;
                    cout << ", entries " << h->GetEntries();
                    cout << endl;
                    h->Delete();
                    return -9999.;
                }
            }
        }
        if( i < iTreeVector.size() - 1 )
        {
            if( h )
            {
                h->Delete();
            }
            if( hI )
            {
                hI->Delete();
            }
        }
    }
    // means: variable is in all trees constant
    cout << "\t warning: constant variable  " << iVariable << " in ";
    if( iSignal )
    {
        cout << " signal tree";
    }
    else
    {
        cout << " background tree";
    }
    if( h )
    {
        cout << " (mean " << h->GetMean() << ", RMS " << h->GetRMS() << ", entries " << h->GetEntries() << ")";
    }
    else if( hI )
    {
        cout << " (mean " << hI->GetMean() << ", RMS " << hI->GetRMS() << ", entries " << hI->GetEntries() << ")";
    }
    cout << ", checked " << iTreeVector.size() << " trees";
    cout << endl;
    double i_mean = -9999.;
    if( h )
    {
        i_mean = h->GetMean();
        h->Delete();
    }
    else if( hI )
    {
        i_mean = hI->GetMean();
        hI->Delete();
    }
    
    return i_mean;
}

/*!

     train the MVA

*/

bool trainGammaHadronSeparation( VTMVARunData* iRun, 
                                 unsigned int iEnergyBin, unsigned int iZenithBin, 
                                 string fRunOption )
{
    return train( iRun, iEnergyBin, iZenithBin, true, fRunOption );
}

bool trainReconstructionQuality( VTMVARunData* iRun, 
                                 unsigned int iEnergyBin, unsigned int iZenithBin, 
                                 string fRunOption )
{
    return train( iRun, iEnergyBin, iZenithBin, false, fRunOption );
}


bool train( VTMVARunData* iRun, 
            unsigned int iEnergyBin, unsigned int iZenithBin, 
            bool iTrainGammaHadronSeparation,
            string fRunOption )
{
    // sanity checks
    if( !iRun )
    {
        return false;
    }
    if( iRun->fEnergyCutData.size() <= iEnergyBin || iRun->fOutputFile.size() <= iEnergyBin )
    {
        cout << "error in train: energy bin out of range " << iEnergyBin;
        return false;
    }
    if( iRun->fZenithCutData.size() < iZenithBin || iRun->fOutputFile[0].size() < iZenithBin )
    {
        cout << "error in train: zenith bin out of range " << iZenithBin;
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

    if( !iRun->fMCxyoffCutSignalOnly )
    {
        iCutBck = iCutBck && iRun->fMCxyoffCut;
    }
    
    // adding training variables
    if( iRun->fTrainingVariable.size() != iRun->fTrainingVariableType.size() )
    {
        cout << "train: error: training-variable vectors have different size" << endl;
        return false;
    }

    // check split mode
    bool iSplitBlock = false;
    if( iRun->fPrepareTrainingOptions.find( "SplitMode=Block" ) != string::npos )
    {
        cout << "train: use option SplitMode=Block" << endl;
        iSplitBlock = true;
    }
     // prepare trees for training and testing with selected events only
     // this step is necessary to minimise the memory impact for the BDT
     // training
     TTree *iSignalTree_reduced = 0;
     TTree *iBackgroundTree_reduced = 0;
     if( fRunOption == "WRITETRAININGEVENTS" )
     {
         iSignalTree_reduced = prepareSelectedEventsTree( iRun,
                           iCutSignal, true,
                           iRun->fResetNumberOfTrainingEvents );
         iBackgroundTree_reduced = prepareSelectedEventsTree( iRun, 
                           iCutBck, false,
                           iRun->fResetNumberOfTrainingEvents );
     
         if( iSignalTree_reduced ) iSignalTree_reduced->Write();
         if( iBackgroundTree_reduced ) iBackgroundTree_reduced->Write();
         if( iRun->getTLRunParameter() ) iRun->getTLRunParameter()->Write();
         cout << "Writing reduced event lists for training: ";
         cout << gDirectory->GetName() << endl;
         exit( EXIT_SUCCESS );
     }
     else
     {
         cout << "Reading training / testing trees from ";
         cout << iRun->fSelectedEventTreeName << endl;
         TFile *iF = new TFile( iRun->fSelectedEventTreeName.c_str() );
         if( iF->IsZombie() )
         {
             cout << "Error open file with pre-selected events: ";
             cout << iRun->fSelectedEventTreeName << endl;
             exit( EXIT_FAILURE );
         }
         iSignalTree_reduced = (TTree*)iF->Get( "data_signal" );
         iBackgroundTree_reduced = (TTree*)iF->Get( "data_background" );
     }
     if( !iSignalTree_reduced || !iBackgroundTree_reduced )
     {
         cout << "Error: failed preparing traing / testing trees" << endl;
         cout << "exiting..." << endl;
         exit( EXIT_SUCCESS );
     }
    // check for number of training and background events
     Long64_t nEventsSignal = iSignalTree_reduced->GetEntries();
     Long64_t nEventsBck = iBackgroundTree_reduced->GetEntries();
     cout << "Reading number of events after cuts: " << endl;
     cout << "\t total number of signal events after cuts: ";
     cout << nEventsSignal;
     cout << " (required are " << iRun->fMinSignalEvents << ")" << endl;
     cout << "\t total number of background events after cuts: ";
     cout << nEventsBck;
     cout << " (required are " << iRun->fMinBackgroundEvents << ")" << endl;
     if( nEventsSignal < iRun->fMinSignalEvents || nEventsBck < iRun->fMinBackgroundEvents  )
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
             
    TMVA::Tools::Instance();
    
    // set output directory
    gSystem->mkdir( iRun->fOutputDirectoryName.c_str() );
    TString iOutputDirectory( iRun->fOutputDirectoryName.c_str() );
    gSystem->ExpandPathName( iOutputDirectory );
    ( TMVA::gConfig().GetIONames() ).fWeightFileDir = iOutputDirectory;
    
    //////////////////////////////////////////
    // defining training class
    TMVA::Factory* factory = new TMVA::Factory( iRun->fOutputFile[iEnergyBin][iZenithBin]->GetTitle(),
            iRun->fOutputFile[iEnergyBin][iZenithBin],
            "V:!DrawProgressBar" );
    TMVA::DataLoader* dataloader = new TMVA::DataLoader( "" );
    ////////////////////////////
    // train gamma/hadron separation
    if( iTrainGammaHadronSeparation )
    {
        // adding signal and background tree
        dataloader->AddSignalTree( iSignalTree_reduced, iRun->fSignalWeight );
        dataloader->AddBackgroundTree( iBackgroundTree_reduced, iRun->fBackgroundWeight );
    }
    ////////////////////////////
    // train reconstruction quality
    else
    {
        dataloader->AddSignalTree( iSignalTree_reduced, iRun->fSignalWeight );
        dataloader->AddRegressionTarget( iRun->fReconstructionQualityTarget.c_str(), iRun->fReconstructionQualityTargetName.c_str() );
    }
    
    // loop over all trainingvariables and add them to TMVA
    // (test first if variable is constant, TMVA will stop when a variable
    //  is constant)
    for( unsigned int i = 0; i < iRun->fTrainingVariable.size(); i++ )
    {
        if( iRun->fTrainingVariable[i].find( "NImages_Ttype" ) != string::npos )
        {
            for( int j = 0; j < iRun->fNTtype; j++ )
            {
                ostringstream iTemp;
                iTemp << iRun->fTrainingVariable[i] << "[" << j << "]";
                ostringstream iTempCut;
                // require at least 2 image per telescope type
                iTempCut << iTemp.str() << ">1";
                char *cstr = new char [iTempCut.str().length()+1];
                std::strcpy (cstr, iTempCut.str().c_str() );
                TCut iCutCC = cstr;
                
                double iSignalMean = 1.;
                double iBckMean    = -1.;
                if( iRun->fCheckValidityOfInputVariables )
                {
                    iSignalMean = checkIfVariableIsConstant( iRun, iCutSignal && iCutCC, iTemp.str(), true, iSplitBlock );
                    iBckMean    = checkIfVariableIsConstant( iRun, iCutBck && iCutCC, iTemp.str(), false, iSplitBlock );
                    cout << "\t mean values(1), signal: " << iSignalMean << " background: " << iBckMean << endl;
                }
                if( ( TMath::Abs( iSignalMean - iBckMean ) > 1.e-6
                        || TMath::Abs( iSignalMean + 9999. ) < 1.e-2 || TMath::Abs( iBckMean + 9999. ) < 1.e-2 )
                        && iSignalMean != 0 && iBckMean != 0 )
                {
                    dataloader->AddVariable( iTemp.str().c_str(), iRun->fTrainingVariableType[i] );
                }
                else
                {
                    cout << "warning: removed constant variable " << iTemp.str() << " from training (added to spectators)" << endl;
                    dataloader->AddSpectator( iTemp.str().c_str() );
                }
                delete[] cstr;
            }
        }
        else
        {
            // check if the training variable is constant
            double iSignalMean = 1.;
            double iBckMean    = -1.;
            if( iRun->fCheckValidityOfInputVariables )
            {
                iSignalMean = checkIfVariableIsConstant( iRun, iCutSignal, iRun->fTrainingVariable[i].c_str(), true, iSplitBlock );
                iBckMean    = checkIfVariableIsConstant( iRun, iCutBck, iRun->fTrainingVariable[i].c_str(), false, iSplitBlock );
                cout << "\t mean values (2), signal: " << iSignalMean << " background: " << iBckMean << endl;
            }
            
            if( TMath::Abs( iSignalMean - iBckMean ) > 1.e-6
                    || TMath::Abs( iSignalMean + 9999. ) < 1.e-2 || TMath::Abs( iBckMean + 9999. ) < 1.e-2 )
            {
                dataloader->AddVariable( iRun->fTrainingVariable[i].c_str(), iRun->fTrainingVariableType[i] );
            }
            else
            {
                cout << "warning: removed constant variable " << iRun->fTrainingVariable[i] << " from training (added to spectators)" << endl;
                dataloader->AddSpectator( iRun->fTrainingVariable[i].c_str() );
            }
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
    if( iTrainGammaHadronSeparation )
    {
        cout << "Preparing training and test tree" << endl;
        // cuts after pre-selection
        TCut iCutSignal_post = iRun->fEnergyCutData[iEnergyBin]->fEnergyCut
                            && iRun->fMultiplicityCuts;
        TCut iCutBck_post = iRun->fEnergyCutData[iEnergyBin]->fEnergyCut
                            && iRun->fMultiplicityCuts;

        dataloader->PrepareTrainingAndTestTree( iCutSignal_post,
                                                iCutBck_post, 
                                                iRun->fPrepareTrainingOptions );
    }
    else
    {
        dataloader->PrepareTrainingAndTestTree( "", iRun->fPrepareTrainingOptions );
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
        // BOOSTED DECISION TREES
        if( iRun->fMVAMethod[i] != "BOXCUTS" )
        {
            if( iTrainGammaHadronSeparation )
            {
                sprintf( htitle, "%s_%u", iRun->fMVAMethod[i].c_str(), i );
            }
            else
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
        //////////////////////////
        // BOX CUTS
        // (note: box cuts needs additional checking, as the code might be outdated)
        else if( iRun->fMVAMethod[i] == "BOXCUTS" )
        {
            stringstream i_opt;
            i_opt << iRun->fMVAMethod_Options[i].c_str();
            for( unsigned int i = 0; i < iRun->fTrainingVariable_CutRangeMin.size(); i++ )
            {
                i_opt << ":CutRangeMin[" << i << "]=" << iRun->fTrainingVariable_CutRangeMin[i];
            }
            for( unsigned int i = 0; i < iRun->fTrainingVariable_CutRangeMax.size(); i++ )
            {
                i_opt << ":CutRangeMax[" << i << "]=" << iRun->fTrainingVariable_CutRangeMax[i];
            }
            for( unsigned int i = 0; i < iRun->fTrainingVariable_VarProp.size(); i++ )
            {
                i_opt << ":VarProp[" << i << "]=" << iRun->fTrainingVariable_VarProp[i];
            }
            sprintf( htitle, "BOXCUTS_%u_%u", iEnergyBin, iZenithBin );
            factory->BookMethod( dataloader, TMVA::Types::kCuts, htitle, i_opt.str().c_str() );
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
    if( !fData->readConfigurationFile( argv[1] ) )
    {
        cout << "error opening or reading run parameter file (";
        cout << argv[1];
        cout << ")" << endl;
        exit( EXIT_FAILURE );
    }
    string fRunOption = "TRAIN";
    if( argc == 3 )
    {
       fRunOption = "WRITETRAININGEVENTS";
    }
    // randomize list of input files
    fData->shuffleFileVectors();
    fData->print();
    
    //////////////////////////////////////
    // read and prepare data files
    if( !fData->openDataFiles( false ) )
    {
        cout << "error opening data files" << endl;
        exit( EXIT_FAILURE );
    }
    
    //////////////////////////////////////
    // train MVA
    // (one training per energy bin)
    cout << "Total number of energy bins: " << fData->fEnergyCutData.size() << endl;
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
            if( fData->fTrainGammaHadronSeparation )
            {
                trainGammaHadronSeparation( fData, i, j, fRunOption );
            }
            if( fData->fTrainReconstructionQuality )
            {
                trainReconstructionQuality( fData, i, j, fRunOption );
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
            if( !root_file )
            {
                cout << "Error finding tvma root file " << endl;
                continue;
            }
            TFile* short_root_file = TFile::Open( iTempS.str().c_str() , "RECREATE" );
            if( !short_root_file->IsZombie() )
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
                    if( ( TH1D* )root_file->Get( hname ) )
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
            rename( iTempS.str().c_str() , iFinalRootFileName.c_str() );
        }
    }
    return 0;
}

