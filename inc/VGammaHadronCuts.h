//! VGammaHadronCuts: class for parameter cut definitions

#ifndef VGammaHadronCuts_H
#define VGammaHadronCuts_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "CData.h"
#include "VAnalysisUtilities.h"
#include "VGammaHadronCutsStatistics.h"
#include "VInstrumentResponseFunctionData.h"
#include "VTMVAEvaluator.h"
#include "VUtilities.h"

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#define VANACUTS_PROBSELECTIONCUTS_MAX 1000

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// class for telescope type dependent multiplicity  cut
////////////////////////////////////////////////////////////////////////////////
class VNTelTypeCut : public TNamed
{
    public:
    
        vector< unsigned int > fTelType_counter;
        vector< ULong64_t >    fTelTypeID;
        unsigned int           fNTelType_min;
        
        VNTelTypeCut();
        ~VNTelTypeCut() {}
        void print();
        void purgeTelTypeIDs();
        bool test( CData* );
        
        ClassDef( VNTelTypeCut, 5 );
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
   dummy class (for compatibility reasons)

*/

class VGammaHadronCutsStats : public TNamed
{
    private:
    
        vector< string > fName;
        
    public:
    
        vector< unsigned int > fN;
        
        VGammaHadronCutsStats() {}
        ~VGammaHadronCutsStats() {}
        
        void printCutStatistics() {}
        void reset() {}
        
        ClassDef( VGammaHadronCutsStats, 3 );
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class VGammaHadronCuts : public VAnalysisUtilities
{
    private:
    
        bool   fDebug;                               // lots of debug output
        
        CData* fData;                                       //! transient
        string fDataDirectory;
        string fInstrumentEpoch;
        
        // cut selector
        int fGammaHadronCutSelector;                            // see description at beginning of VGammaHadronCuts.cpp
        int fDirectionCutSelector;
        
        // array characteristics (number of telescopes, centre of array)
        unsigned int fNTel;
        double       fArrayCentre_X;
        double       fArrayCentre_Y;
        
        // number of possible telescope combinations
        unsigned int fNLTrigs;
        
        // telescope used in analysis (optional)
        vector< unsigned int > fTelToAnalyze;
        
        // values calculated from shower/image parameter
        // (used also in VStereoAnalysis)
        double fMeanImageDistance;
        double fMeanImageLength;
        double fMeanImageWidth;
        
        // event by event cuts (read in from an additional friend tree, used by random forest analysis, pulsar analysis, etc)
        TFile* fProbabilityCut_File;                  //!
        TTree* fProbabilityCut_Tree;                  //!
        int fProbabilityCut_QualityFlag;              // quality flag for probability cut 0: cut estimation failed; >0 successful probability estimation
        unsigned int fProbabilityCut_NSelectors;      // number of elements in fProbabilityCut_SelectionCut[]
        unsigned int fProbabilityCut_ProbID;          // array element to be used from fProbabilityCut_SelectionCut[]
        double fProbabilityCut_SelectionCut[VANACUTS_PROBSELECTIONCUTS_MAX];    // selection cut
        
        //////////////////////////
        // TMVA evaluator
        VTMVAEvaluator* fTMVAEvaluator;                             //!
        string          fTMVA_MVAMethod;
        unsigned int    fTMVA_MVAMethodCounter;
        string          fTMVAWeightFile;
        unsigned int    fTMVAWeightFileIndex_Emin;
        unsigned int    fTMVAWeightFileIndex_Emax;
        unsigned int    fTMVAWeightFileIndex_Zmin;
        unsigned int    fTMVAWeightFileIndex_Zmax;
        double          fTMVAEnergyStepSize;
        map< unsigned int, double > fTMVASignalEfficiency;
        map< unsigned int, double > fTMVA_MVACut;
        string          fTMVA_MVACutGraphFileName;
        double          fTMVA_MVACutGraphSmoothing;
        double          fTMVA_MVACutGraphSmoothingMax;
        double          fTMVA_MVACutGraphConstantCutEnergy_TeV;
        double          fTMVAProbabilityThreshold;
        string          fTMVAOptimizeSignalEfficiencyParticleNumberFile;
        double          fTMVAParticleNumberFile_Conversion_Rate_to_seconds;
        double          fTMVAOptimizeSignalEfficiencySignificance_Min;
        double          fTMVAOptimizeSignalEfficiencySignalEvents_Min;
        double          fTMVAOptimizeSignalEfficiencyObservationTime_h;
        double          fTMVAFixedSignalEfficiencyMax;
        double          fTMVAMinSourceStrength;
        double          fTMVAFixedThetaCutMin;
        double          fTMVA_EvaluationResult;
        VTMVAEvaluatorResults* fTMVAEvaluatorResults;
        vector< TGraphAsymmErrors* > fMVACutGraphs;
        
        // parameters for energy dependent theta2 cuts
        // (implemented for MC only)
        string fFileNameAngRes;
        TFile* fFileAngRes;                                         //!
        string fF1AngResName;
        TF1*   fF1AngRes;
        double       fAngRes_ScalingFactor;
        double       fAngRes_AbsoluteMinimum;
        double       fAngRes_AbsoluteMaximum;
        unsigned int fAngResContainmentProbability;
        
        //////////////////////////
        // energy dependent cuts
        map< string, TGraph* > fEnergyDependentCut;
        
        // cut statistics
        VGammaHadronCutsStatistics* fStats;                       //!
        
        bool   applyProbabilityCut( int i, bool fIsOn );
        bool   applyFrogsCut();
        bool   applyDeepLearnerCut();
        double getEnergyDependentCut( double energy_TeV, TGraph* iG, bool bUseEvalue = true, bool bMaxCut = true );
        TGraph* getEnergyDependentCut( string iCutName );
        bool   getEnergyDependentCutFromFile( string iFileName, string iVariable );
        double getMeanGoodness( double, double, double, double, int );
        bool   initAngularResolutionFile();
        bool   initProbabilityCuts( int irun );
        bool   initProbabilityCuts( string iDir );
        bool   initTMVAEvaluator( string iTMVAFile, unsigned int iTMVAWeightFileIndex_Emin, unsigned int iTMVAWeightFileIndex_Emax, unsigned int iTMVAWeightFileIndex_Zmin, unsigned int iTMVAWeightFileIndex_Zmax, double iTMVAEnergy_StepSize );
        string getTelToAnalyzeString();
        
        
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        
    public:
        E_ReconstructionType fReconstructionType;
        
        // stereo cuts
        double fCut_MeanImageDistance_min;
        double fCut_MeanImageDistance_max;
        double fCut_MeanImageLength_min;
        double fCut_MeanImageLength_max;
        double fCut_MeanImageWidth_min;
        double fCut_MeanImageWidth_max;
        
        double fCut_Theta2_min;
        double fCut_Theta2_max;
        double fCut_Chi2_min;
        double fCut_Chi2_max;
        double fCut_Dir_max;
        double fCut_Size_min;
        double fCut_Size_max;
        double fCut_MSCW_min;
        double fCut_MSCW_max;
        double fCut_MSCL_min;
        double fCut_MSCL_max;
        double fCut_MSW_min;
        double fCut_MSW_max;
        double fCut_MSL_min;
        double fCut_MSL_max;
        vector< int > fCut_ImgSelect;
        double fCut_CameraFiducialSize_min;
        double fCut_CameraFiducialSize_max;
        bool   bMCCuts;
        double fCut_CameraFiducialSize_MC_min;
        double fCut_CameraFiducialSize_MC_max;
        double fCut_AverageCoreDistanceToTelescopes_min;
        double fCut_AverageCoreDistanceToTelescopes_max;
        int    fCut_AverageCoreDistanceToTelescopes_nimages;
        double fCut_MinimumCoreDistanceToTelescopes_max;
        double fCut_dE_min;
        double fCut_dE_max;
        double fCut_EmultiplicityRatio_min;
        double fCut_EChi2_min;
        double fCut_EChi2_max;
        double fCut_Erec_min;
        double fCut_Erec_max;
        double fCut_Emmission_min;
        double fCut_Emmission_max;
        int    fCut_NImages_min;
        int    fCut_NImages_max;
        unsigned int    fCut_DispNImages_min;
        unsigned int    fCut_DispNImages_max;
        double fCut_SizeSecondMax_min;
        double fCut_SizeSecondMax_max;
        double fProbabilityCut;
        vector <double> fProbabilityCutRangeLower;
        vector <double> fProbabilityCutRangeUpper;
        
        vector< VNTelTypeCut* > fNTelTypeCut;
        
        VGammaHadronCuts();
        ~VGammaHadronCuts();
        
        bool   applyDirectionCuts( bool bCount = false, double x0 = -99999., double y0 = -99999. );
        bool   applyEnergyReconstructionQualityCuts( unsigned int iEnergyReconstructionMethod = 0, bool bCount = false );
        bool   applyInsideFiducialAreaCut( bool bCount = false );
        bool   applyInsideFiducialAreaCut( float Xoff, float Yoff, bool bCount = false );
        bool   applyMCXYoffCut( double x, double y, bool bCount = false );
        bool   applyMeanReducedScaledStereoShapeCuts();
        bool   applyMeanStereoShapeCuts();
        bool   applyMeanScaledStereoShapeCuts();
        bool   applyStereoQualityCuts( unsigned int iEnergyReconstructionMethod = 0, bool bCount = false, int iEntry = 0, bool fIsOn = false );
        bool   applyStereoQualityCuts_Chi2( bool bCount = false );
        bool   applyStereoQualityCuts_TableVariables( bool bCount = false, unsigned int iEnergyReconstructionMethod = 0 );
        bool   applyStereoQualityCuts_NImages( bool bCount = false );
        bool   applyStereoQualityCuts_ImageSelection( bool bCount = false );
        bool   applyStereoQualityCuts_TelescopeCoreDistance( bool bCount = false );
        bool   applyStereoQualityCuts_SecondMaxCut( bool bCount = false );
        bool   applyStereoShapeCuts();
        bool   applyTMVACut( int i );
        bool   applyTelTypeTest( bool bCount = false );
        
        TF1*   getAngularResolutionFunction()
        {
            return fF1AngRes;
        }
        double getAngularResolutionAbsoluteMinimum()
        {
            return fAngRes_AbsoluteMinimum;
        }
        double getAngularResolutionAbsoluteMaximum()
        {
            return fAngRes_AbsoluteMaximum;
        }
        double getAngularResolutionScaleFactor()
        {
            return fAngRes_ScalingFactor;
        }
        double getArrayCentre_X()
        {
            return fArrayCentre_X;
        }
        double getArrayCentre_Y()
        {
            return fArrayCentre_Y;
        }
        int    getDirectionCutSelector()
        {
            return fDirectionCutSelector;
        }
        int    getGammaHadronCutSelector()
        {
            return fGammaHadronCutSelector;
        }
        double getMeanImageDistance()
        {
            return fMeanImageDistance;
        }
        double getMeanImageLength()
        {
            return fMeanImageLength;
        }
        double getMeanImageWidth()
        {
            return fMeanImageWidth;
        }
        unsigned int getAngularResolutionContainmentRadius()
        {
            return fAngResContainmentProbability;
        }
        double getProbabilityCut_Selector( unsigned int iID = 0 )
        {
            if( iID < fProbabilityCut_NSelectors )
            {
                return fProbabilityCut_SelectionCut[iID];
            }
            else
            {
                return -1;
            }
        }
        double getProbabilityCutAlpha( bool fIsOn );
        double getTheta2Cut_min( double e = 0.1 )
        {
            if( e > 0. )
            {
                return fCut_Theta2_min;
            }
            else
            {
                return 0.;
            }
        }
        double getTheta2Cut_max()
        {
            return fCut_Theta2_max;
        }
        double getTheta2Cut_max( double e );                           // get theta2 max cut (might be energy dependent)    [TeV] energy (linear)
        TGraph* getTheta2Cut_TMVA_max()
        {
            return getEnergyDependentCut( "TMVABoxCut_Theta2_max" );
        }
        TGraph* getTheta2Cut_IRF_Max()
        {
            return getEnergyDependentCut( "IRFAngRes" );
        }
        double getTMVA_EvaluationResult()
        {
            return fTMVA_EvaluationResult;
        }
        VTMVAEvaluatorResults* getTMVAEvaluatorResults()
        {
            return fTMVAEvaluatorResults;
        }
        vector< TGraphAsymmErrors* > getMVACutGraphs()
        {
            return fMVACutGraphs;
        }
        void   initialize();
        bool   isGamma( int i = 0, bool bCount = false, bool fIsOn = true );
        bool   isMCCuts()
        {
            return bMCCuts;
        }
        void   newEvent( bool iFillStats = true );
        void   printCutSummary();
        void   printCutStatistics()
        {
            if( fStats )
            {
                fStats->printCutStatistics();
            }
        }
        void   printDirectionCuts();
        void   printEnergyDependentCuts();
        void   printSignalEfficiency();
        void   printTMVA_MVACut();
        bool   readCuts( string i_cutfilename, int iPrint = 1 );
        void   resetCutValues();
        void   resetCutStatistics();
        void   initializeCuts( int irun = -1, string iDir = "" );
        void   setDataDirectory( string id )
        {
            fDataDirectory = id;
        }
        bool   setDataTree( CData* idata );
        void   setDebug( bool iB = false )
        {
            fDebug = iB;
        }
        void   setEnergyCuts( double imin, double imax )
        {
            fCut_Erec_min = imin;
            fCut_Erec_max = imax;
        }
        void   setInstrumentEpoch( string iEpoch )
        {
            fInstrumentEpoch = iEpoch;
        }
        bool   setIRFGraph( TGraphErrors* g );
        void   setNTel( unsigned int itel,  double iX = 0., double iY = 0. )
        {
            fNTel = itel;
            fArrayCentre_X = iX;
            fArrayCentre_Y = iY;
        }
        void   setTelToAnalyze( vector< unsigned int > iTelToAnalyze )
        {
            fTelToAnalyze = iTelToAnalyze;
        }
        void   setTheta2Cut( double it2 )
        {
            fCut_Theta2_max = it2;
        }
        void   terminate( bool iShort = false );
        bool   useDeepLearnerCuts()
        {
            return ( fReconstructionType == DEEPLEARNER );
        }
        bool   useFrogsCuts()
        {
            return ( fReconstructionType == FROGS );
        }
        bool   useTMVACuts()
        {
            return ( fGammaHadronCutSelector / 10 == 4 );
        }
        void setReconstructionType( E_ReconstructionType type )
        {
            fReconstructionType = type;
        }
        ClassDef( VGammaHadronCuts, 67 );
};
#endif
