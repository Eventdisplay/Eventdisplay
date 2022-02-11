//!< VDL2Writer calculate effective areas and energy spectra

#ifndef VDL2Writer_H
#define VDL2Writer_H

#include "CData.h"
#include "VDL2Tree.h"
#include "VTMVAEvaluator.h"

#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

#define VDST_MAXTELESCOPES  180

class VTMVA_eval_dist
{
   public:
    VTMVAEvaluator *fTMVAEvaluator;
    string fInstrumentEpoch;
    string fTMVA_MVAMethod;
    unsigned int fTMVA_MVAMethodCounter;
    unsigned int fTMVAWeightFileIndex_Emin;
    unsigned int fTMVAWeightFileIndex_Emax;
    unsigned int fTMVAWeightFileIndex_Zmin;
    unsigned int fTMVAWeightFileIndex_Zmax;
    double fTMVAEnergyStepSize;
    string fTMVAWeightFile;

    bool fIsZombie;

    VTMVA_eval_dist( string i_fTMVA_MVAMethod,
                     unsigned int i_fTMVA_MVAMethodCounter,
                     string i_fTMVAWeightFile,
                     unsigned int i_fTMVAWeightFileIndex_Emin, unsigned int i_fTMVAWeightFileIndex_Emax,
                     unsigned int i_fTMVAWeightFileIndex_Zmin, unsigned int i_fTMVAWeightFileIndex_Zmax,
                     double i_fTMVAEnergyStepSize, string i_fInstrumentEpoch,
                     CData *d );
   ~VTMVA_eval_dist();
    double evaluate();
    
};

class VDL2Writer
{
    private:

        vector< string > fdatafile;

       // distance (xoff bins)
        vector< string > dist_mean;
        vector< float > dist_min;
        vector< float > dist_max;
        bool fCopyDataTree;

        // TMVA
        vector< VTMVA_eval_dist* > fTMVA;
        string fInstrumentEpoch;
        string fTMVA_MVAMethod;
        unsigned int fTMVA_MVAMethodCounter;
        unsigned int fTMVAWeightFileIndex_Emin;
        unsigned int fTMVAWeightFileIndex_Emax;
        unsigned int fTMVAWeightFileIndex_Zmin;
        unsigned int fTMVAWeightFileIndex_Zmax;
        double fTMVAEnergyStepSize;
        vector< string > fTMVAWeightFile;
        
        // event data
        VDL2Tree *fDL2DataTree;

        unsigned int get_tmva_distance_bin( float distance );
        bool initializeTMVAEvaluators( CData *d );
        bool readConfigFile( string iConfigFile );
        
    public:
    
        VDL2Writer( string iConfigFile );
       ~VDL2Writer();
        bool  copyDataTree() { return fCopyDataTree; }
        bool  fill( CData* d );
        vector< string > getDataFile() { return fdatafile; }
        TTree* getEventDataTree()
        {
             if( fDL2DataTree ) return fDL2DataTree->getDL2Tree();
             return 0;
        }
};
#endif
