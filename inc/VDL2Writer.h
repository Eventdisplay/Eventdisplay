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

#define VDST_MAXTELESCOPES  180

class VTMVA_eval_dist
{
    public:
        VTMVAEvaluator* fTMVAEvaluator;
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
                         CData* d );
        ~VTMVA_eval_dist();
        double evaluate();
        
};

class VDL2Writer
{
    private:
    
        string fdatafile;
        
        // distance (xoff bins)
        vector< string > dist_mean;
        vector< float > dist_min;
        vector< float > dist_max;
        
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
        TTree* fDL2DataTree;
        
        Int_t runNumber;
        Int_t eventNumber;
        Float_t         ArrayPointing_Elevation;
        Float_t         ArrayPointing_Azimuth;
        Double_t        MCe0;
        Double_t        MCxcore;
        Double_t        MCycore;
        Double_t        MCaz;
        Double_t        MCze;
        Double_t        MCxoff;
        Double_t        MCyoff;
        Int_t           NImages;
        ULong64_t       ImgSel;
        Double_t        Ze;
        Double_t        Az;
        Double_t        Xoff;
        Double_t        Yoff;
        Double_t        Xoff_derot;
        Double_t        Yoff_derot;
        Float_t         Xoff_intersect;
        Float_t         Yoff_intersect;
        Double_t        Xcore;
        Double_t        Ycore;
        Double_t        MSCW;
        Double_t        MSCL;
        Double_t        Chi2;
        Double_t        ErecS;
        Double_t        EChi2S;
        Double_t        dES;
        Double_t        ES[VDST_MAXTELESCOPES];
        Double_t        SizeSecondMax;
        Float_t         EmissionHeight;
        Float_t         EmissionHeightChi2;
        UInt_t          DispNImages;
        Int_t           NTtype;
        ULong64_t       TtypeID[VDST_MAXTELESCOPES];
        UInt_t          NImages_Ttype[VDST_MAXTELESCOPES];
        Double_t        R[VDST_MAXTELESCOPES];
        float fCut_MVA;
        
        void fillEventDataTree( float iMVA );
        bool initializeTMVAEvaluators( CData* d );
        bool readConfigFile( string iConfigFile );
        
    public:
    
        VDL2Writer( string iConfigFile );
        ~VDL2Writer();
        bool  fill( CData* d );
        string getDataFile()
        {
            return fdatafile;
        }
        TTree* getEventDataTree()
        {
            return fDL2DataTree;
        }
};
#endif
