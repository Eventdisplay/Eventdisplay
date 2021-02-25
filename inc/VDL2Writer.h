//!< VDL2Writer calculate effective areas and energy spectra

#ifndef VDL2Writer_H
#define VDL2Writer_H

#include "CData.h"
#include "VGammaHadronCuts.h"
#include "VInstrumentResponseFunctionRunParameter.h"

#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"

#include <iomanip>
#include <string>

using namespace std;

class VDL2Writer
{
    private:
        VInstrumentResponseFunctionRunParameter* fRunPara;
        VGammaHadronCuts* fCuts;
        bool fIgnoreEnergyReconstruction;
        bool fIsotropicArrivalDirections;
        bool fTelescopeTypeCutsSet;
        
        // event data
        TTree *fEventTreeCuts;
        int fCut_Class;
        float fCut_MVA;

        void               fillEventDataTree( int iCutClass, float iMVA );
        
    public:
    
        VDL2Writer( VInstrumentResponseFunctionRunParameter*, VGammaHadronCuts* );
       ~VDL2Writer();
        bool  fill( CData* d, unsigned int iMethod );
        
        TTree* getEventCutDataTree()
        {
             return fEventTreeCuts;
        }
        void setIgnoreEnergyReconstructionCuts( bool iB = false )
        {
            fIgnoreEnergyReconstruction = iB;
        }
        void setIsotropicArrivalDirections( bool iB = false )
        {
            fIsotropicArrivalDirections = iB;
        }
        void setTelescopeTypeCuts( bool iB = true )
        {
            fTelescopeTypeCutsSet = iB;
        }
};
#endif
