//! VTableLookupRunParameter parameter storage class

//   ClassDef VERSION NUMBER HAS TO BE INCREASED BY ONE AFTER ANY CHANGES HERE

#ifndef VTABLELOOKUPRUNPARAMTER_H
#define VTABLELOOKUPRUNPARAMTER_H

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TNamed.h>
#include <TSystem.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Ctelconfig.h"
#include "VEvndispRunParameter.h"
#include "VEvndispReconstructionParameter.h"
#include "VGlobalRunParameter.h"
#include "VUtilities.h"

using namespace std;

class VTableLookupTelToAnalyze
{
    public:
        unsigned int fTelToAnalyze;
        unsigned int fTelID;
        unsigned int fTelID_hyperArray;
        ULong64_t    fTelType;
        double       fWeight;
        
        VTableLookupTelToAnalyze();
        virtual ~VTableLookupTelToAnalyze() {}
        void print();
        
        ClassDef( VTableLookupTelToAnalyze, 1 ); //for any changes to this file: increase this number
};


class VTableLookupRunParameter : public TNamed, public VGlobalRunParameter
{
    private:
    
        bool fillInputFile_fromList( string iList );
        void fillTelescopeTypeDependentWeights();
        void printCTA_MC_offaxisBins();
        bool readTelescopeToAnalyze( string iEvndispRootFile );
        bool readTelescopeToAnalyze( string iTelescopeList_sim_telarray_Counting,
                                     string iEvndispRootFile,
                                     bool iSimTelArrayCounting );
        bool readRunParameters( string iFile );
        bool readTelTypeDepdendentWeights( string iFile );
        void setCTA_MC_offaxisBins();
        
    public:
    
        // debug levels 0 = off, 1 = default debug level, 2 = detailed
        unsigned int fDebug;
        // debug printing
        string printpara;
        
        // true for table filling, false for table reading
        bool         fWriteTables;
        
        // list of evndisp input files
        vector< string > inputfile;
        // name of lookup table file
        string tablefile;
        
        // reconstructed method (quality cut) read from
        // evndisp file
        int  rec_method;
        // input size type is 'pe' (not [dc])
        bool fPE;
        // use median energy (instead of mean)
        int fUseMedianEnergy;
        // use image selection applied already in evndisp for the table reading/writing
        bool  fUseSelectedImagesOnly;
        // telescopes used in analysis (old style)
        vector< unsigned int > fTelToAnalyse;
        // telescopes to be used in analysis (new style)
        // should be always of size ntel
        vector< VTableLookupTelToAnalyze* > fTelToAnalyzeData;
        // maximum number of events read
        Long64_t fNentries;
        // event selection cuts for table filling and reading (only used for energy reconstruction)
        double fEventSelectionCut_lossCutMax;
        // event selection cuts for table filling and reading (only used for energy reconstruction)
        double fEventSelectionCut_distanceCutMax;
        
        // run parameter file
        string fRunParameterFile;
        // list of telescopes (subarrays) to be active in analysis
        string fTelescopeList_sim_telarray_Counting;
        // file with telescope type dependent weights
        string fTelescopeType_weightFile;
        // list with telescope type dependent weights
        map< ULong64_t, double > fTelescopeType_weight;
        // quality cut level
        unsigned int fQualityCutLevel;
        
        // use lookup tables for time gradient (optional)
        bool fUsetimeGradientLookupTables;
        // use lookup tables for FROGS analysis (optional)
        bool fUsefrogsGoodnessTables;
        
        //////////////////////////////////////////
        // parameters for table filling (writing) only
        // minimum number of showers required per table bin
        float fMinRequiredShowerPerBin;
        // write (for debugging) all 1D distribution to
        // lookup table file
        bool fWrite1DHistograms;
        // spectral index used to re-weight events while filling the
        // lookup tables
        double fSpectralIndex;
        // use core position and direction from model3D analysis
        bool bUseModel3DStereoParameters;
        // parameters set to fix lookup table file directory structure
        // zenith angle
        double ze;
        // wobble offset
        int fWobbleOffset;
        // NSB (pedvars) level
        int fNoiseLevel;
        // cuts applied for table filling
        unsigned int fTableFillingCut_NImages_min;
        double       fTableFillingCut_WobbleCut_max;
        double       fMC_distance_to_cameracenter_min;
        double       fMC_distance_to_cameracenter_max;
        double       fmaxdist;   // note: same for all telescope types
        double       fmaxloss;   // note: same for all telescope types
        double       fminfui ;   // note: same for all telescope types
        double       fminsize;   // note: same for all telescope times
        // seed for random selection of showers before table filling
        double       fSelectRandom;
        int          fSelectRandomSeed;
        // definition of offaxis bins (CTA only)
        vector< double > fCTA_MC_offaxisBin_min;
        vector< double > fCTA_MC_offaxisBin_max;
        
        //////////////////////////////////////////
        // parameters for table reading only
        // name of mscw output file
        string outputfile;
        // this is a MC file
        bool isMC;
        // overwrite or update the current output file
        string writeoption;
        // write triggered but not reconstructed events
        bool bNoNoTrigger;
        // write reconstructed events only
        int  bWriteReconstructedEventsOnly;
        // short version of the output tree
        bool bShortTree;
        // copy MC tree over to output file (from evndisp files)
        bool bWriteMCPars;
        // maximum time (in s) used of this run
        double fMaxRunTime;
        // parameters to be used in anasum
        double meanpedvars;                       // mean pedvar
        vector< double > pedvars;                 // mean pedvar per telescope
        // use MC histogram with this spectral index
        vector< double > fAddMC_spectral_index;
        
        // rerun stereo reconstruction (MC only)
        bool  fRerunStereoReconstruction;
        double fRerunStereoReconstruction_minAngle;
        string fRerunStereoReconstruction_BDTFileName;
        unsigned int fRerunStereoReconstruction_BDTNImages_max;
        string fEnergyReconstruction_BDTFileName;
        string fCoreReconstruction_BDTFileName;
        string fDispError_BDTFileName;
        float  fDispError_BDTWeight;
        
        // functions...
        VTableLookupRunParameter();
        ~VTableLookupRunParameter() {}
        
        bool fillParameters( int argc, char* argv[] );
        void print( int iB = 0 );
        void printHelp();
        
        ClassDef( VTableLookupRunParameter, 46 ); //for any changes to this file: increase this number
};
#endif
