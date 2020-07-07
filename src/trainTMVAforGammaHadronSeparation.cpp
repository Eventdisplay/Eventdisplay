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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "VEvndispRunParameter.h"
#include "VTMVARunData.h"

using namespace std;

bool train( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin, bool iGammaHadronSeparation );
bool trainGammaHadronSeparation( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin );
bool trainReconstructionQuality( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin );

/*
 * get number of requested training events
 *
 */
Long64_t getNumberOfRequestedTrainingEvents( string a, bool iSignal )
{
   string b = "nTrain_Signal=";
   if( !iSignal )
   {
      b = "nTrain_Background=";
   }
   // extensive string gymnastics
   Long64_t nS = (Long64_t)atoi( a.substr( a.find( b ) + b.size(), a.find( ":", a.find( b ) - a.find( b ) - b.size() ) ).c_str() );

   return nS;
}

/*
 * check settings for number of training events;
 * if required reset event numbers
 *
 * note! assume name number of training and testing events
 */
string resetNumberOfTrainingEvents( string a, Long64_t n, bool iSignal )
{
   Long64_t nS = getNumberOfRequestedTrainingEvents( a, iSignal );
   if( nS == 0 )
   {
       return a;
   }
   cout << "\t requested are " << nS;
   cout << ", found are " << n << endl;
   if( n > nS )
   {
      cout << "\t (resetting to zero)" << endl;;
      stringstream r;
      r << nS;
      size_t pos = a.find( r.str() );
      while( pos != string::npos)
      {
           a.replace( pos, r.str().size(), "0" );
           pos = a.find( r.str() );
      }
   }
   return a;
}
   

/*
 * get number of signal or background events after cuts
 *
 */
Long64_t getNumberOfEventsAfterCuts( VTMVARunData* iRun, TCut iCut, bool iSignal, bool iResetEventNumbers )
{
    if( !iRun )
    {
        return 0;
    }
    Long64_t nS = getNumberOfRequestedTrainingEvents( iRun->fPrepareTrainingOptions, iSignal );
    vector< TChain* > iTreeVector;
    Long64_t n = 0;
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
        if( iTreeVector[i] )
        {
             TTree *t = (TTree*)iTreeVector[i]->Clone("a");
             t->Draw( ">>elist", iCut );
             TEventList *elist = (TEventList*)gDirectory->Get("elist");
             if( elist )
             {
                 n += elist->GetN();
                 delete elist;
             }
             delete t;
             if( !iResetEventNumbers && n > nS )
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
                 cout << "(" << nS << ")";
                 cout << " after " << i+1 << " tree(s)" << endl;
                 cout << "\t applied cut: " << iCut << endl;
                 return n;
             }
        }
    }
    return n;
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

bool trainGammaHadronSeparation( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin )
{
    return train( iRun, iEnergyBin, iZenithBin, true );
}

bool trainReconstructionQuality( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin )
{
    return train( iRun, iEnergyBin, iZenithBin, false );
}


bool train( VTMVARunData* iRun, unsigned int iEnergyBin, unsigned int iZenithBin, bool iTrainGammaHadronSeparation )
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
        // adding signal and background trees
        for( unsigned int i = 0; i < iRun->fSignalTree.size(); i++ )
        {
            dataloader->AddSignalTree( iRun->fSignalTree[i], iRun->fSignalWeight );
        }
        for( unsigned int i = 0; i < iRun->fBackgroundTree.size(); i++ )
        {
            dataloader->AddBackgroundTree( iRun->fBackgroundTree[i], iRun->fBackgroundWeight );
        }
    }
    ////////////////////////////
    // train reconstruction quality
    else
    {
        for( unsigned int i = 0; i < iRun->fSignalTree.size(); i++ )
        {
            dataloader->AddRegressionTree( iRun->fSignalTree[i], iRun->fSignalWeight );
        }
        dataloader->AddRegressionTarget( iRun->fReconstructionQualityTarget.c_str(), iRun->fReconstructionQualityTargetName.c_str() );
    }
    
    // quality cuts before training
    TCut iCutSignal = iRun->fQualityCuts && iRun->fQualityCutsSignal 
                   && iRun->fMCxyoffCut && iRun->fAzimuthCut &&
                   iRun->fEnergyCutData[iEnergyBin]->fEnergyCut 
                   && iRun->fZenithCutData[iZenithBin]->fZenithCut;

    TCut iCutBck = iRun->fQualityCuts && iRun->fQualityCutsBkg 
                && iRun->fAzimuthCut 
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
    // check for number of training and background events
     cout << "Reading number of events after cuts: " << endl;
     Long64_t nEventsSignal = getNumberOfEventsAfterCuts( iRun, 
                           iCutSignal, true,
                           iRun->fResetNumberOfTrainingEvents );
     cout << "\t total number of signal events after cuts: " << nEventsSignal << endl;
     Long64_t nEventsBck = getNumberOfEventsAfterCuts( iRun, 
                           iCutBck, false,
                           iRun->fResetNumberOfTrainingEvents );
     cout << "\t total number of background events after cuts: " << nEventsBck << endl;
     if( nEventsSignal < 10 || nEventsBck < 10 )
     {
         cout << "Error: not enough training events" << endl;
         cout << "exiting..." << endl;
         exit( EXIT_SUCCESS );
     }
     // check if this number if consistent with requested number of training events
     // required at least 10 events
     if( iRun->fResetNumberOfTrainingEvents )
     {
          cout << "Checking number of training events: " << endl;
          iRun->fPrepareTrainingOptions = resetNumberOfTrainingEvents( iRun->fPrepareTrainingOptions, nEventsSignal, true );
          iRun->fPrepareTrainingOptions = resetNumberOfTrainingEvents( iRun->fPrepareTrainingOptions, nEventsBck, false );
          cout << "Updated training options: " <<  iRun->fPrepareTrainingOptions << endl;
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
                TCut iCutCC = iTempCut.str().c_str();
                
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
    if( iTrainGammaHadronSeparation )
    {
        cout << "Preparing training and test tree" << endl;
        dataloader->PrepareTrainingAndTestTree( iCutSignal, iCutBck, iRun->fPrepareTrainingOptions );
    }
    else
    {
        dataloader->PrepareTrainingAndTestTree( iCutSignal, iRun->fPrepareTrainingOptions );
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
    if( argc != 2 )
    {
        cout << endl;
        cout << "./trainTMVAforGammaHadronSeparation <configuration file>" << endl;
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
    fData->print();
    
    //////////////////////////////////////
    // read and prepare data files
    if( !fData->openDataFiles() )
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
                trainGammaHadronSeparation( fData, i, j );
            }
            if( fData->fTrainReconstructionQuality )
            {
                trainReconstructionQuality( fData, i, j );
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

