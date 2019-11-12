//! VTMVARunData run parameter data class for TMVA optimization

#ifndef VTMVARunData_H
#define VTMVARunData_H

#include "TChain.h"
#include "TCut.h"
#include "TFile.h"
#include "TEntryList.h"
#include "TMath.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "VTMVARunDataEnergyCut.h"
#include "VTMVARunDataZenithCut.h"
#include "VUtilities.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// data class with all signal/background files names and run parameters

class VTMVARunData : public TNamed
{
    private:
    
        bool              fDebug;
        
    public:
    
        string            fName;
        
        // run type
        bool fTrainGammaHadronSeparation;
        bool fTrainReconstructionQuality;
        
        // output file
        string            fOutputFileName;
        string            fOutputDirectoryName;
        vector< vector< TFile* > >  fOutputFile;
        
        // training options
        bool              fCheckValidityOfInputVariables;
        
        // training data
        double            fSignalWeight;
        vector< string >  fSignalFileName;
        vector< TChain* > fSignalTree;
        double            fBackgroundWeight;
        vector< string >  fBackgroundFileName;
        vector< TChain* > fBackgroundTree;
        
        // list of training variables
        vector< string >  fTrainingVariable;
        vector< char >    fTrainingVariableType;
        vector< float >   fTrainingVariable_CutRangeMin;
        vector< float >   fTrainingVariable_CutRangeMax;
        vector< string >  fTrainingVariable_VarProp;
        
        // spectator variables
        vector< string > fSpectatorVariable;
        
        // quality and energy and zenith cuts
        unsigned int      fMinSignalEvents;
        unsigned int      fMinBackgroundEvents;
        TCut              fQualityCuts;
        TCut              fQualityCutsBkg;
        TCut              fQualityCutsSignal;
        TCut              fMCxyoffCut;
        bool              fMCxyoffCutSignalOnly;
        TCut              fAzimuthCut;
        string            fPrepareTrainingOptions;
        vector< VTMVARunDataEnergyCut* > fEnergyCutData;
        vector< VTMVARunDataZenithCut* > fZenithCutData;
        
        // analysis variables
        int               fNTtype;
        
        // MVA methods
        vector< string >  fMVAMethod;
        vector< string >  fMVAMethod_Options;
        
        // reconstruction quality target
        string            fReconstructionQualityTarget;
        string            fReconstructionQualityTargetName;
        
        VTMVARunData();
        ~VTMVARunData() {}
        void print();
        bool readConfigurationFile( char* );
        bool openDataFiles();
        void setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        void setName( string iN )
        {
            fName = iN;
        }
        
        ClassDef( VTMVARunData, 9 );
};

#endif
