//!< VDL2Writer calculate effective areas and energy spectra

#ifndef VDL2Writer_H
#define VDL2Writer_H

#include "CData.h"
#include "VTMVAEvaluator.h"

#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

class VDL2Writer
{
    private:

        string fdatafile;

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
        
        // event data
        TTree *fEventTreeCuts;
        float fCut_MVA;

        void fillEventDataTree( float iMVA );
        bool initializeTMVAEvaluator( CData *d );
        bool readConfigFile( string iConfigFile );
        
    public:
    
        VDL2Writer( string iConfigFile );
       ~VDL2Writer();
        bool  fill( CData* d );
        string getDataFile() { return fdatafile; }
        TTree* getEventCutDataTree()
        {
             return fEventTreeCuts;
        }
};
#endif
