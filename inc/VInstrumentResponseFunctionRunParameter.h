//! VInstrumentResponseFunctionRunParameter run parameters for response function calculator (effective areas)

#ifndef VInstrumentResponseFunctionRunParameter_H
#define VInstrumentResponseFunctionRunParameter_H

#include "Ctelconfig.h"
#include "VEvndispRunParameter.h"
#include "VMonteCarloRunHeader.h"
#include "VTableLookupRunParameter.h"
#include "CData.h"
#include "VEnergySpectrumfromLiterature.h"

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TChain.h"
#include "TF1.h"
#include "TNamed.h"

using namespace std;

class VInstrumentResponseFunctionRunParameter : public TNamed
{
    private:
    
        bool            readRunParameters( string ifilename );
        bool            readCRSpectralParameters();
        
    public:
    
        unsigned int    fFillingMode;              // filling mode
        bool            fEffArea_short_writing;    // short/long tree writing
        
        string          fCutFileName;
        string          fInstrumentEpoch;
        vector< unsigned int > fTelToAnalyse;             // telescopes used in analysis (optional, not always filled)
        int             fGammaHadronCutSelector;
        int             fDirectionCutSelector;
        
        E_ReconstructionType fReconstructionType;
        unsigned int    fEnergyReconstructionMethod;
        unsigned int    fEnergyAxisBins_log10;
        bool            fIgnoreEnergyReconstructionQuality;
        unsigned int    fNSpectralIndex;
        double          fSpectralIndexMin;
        double          fSpectralIndexStep;
        vector< double > fSpectralIndex;
        double          fMCEnergy_min;
        double          fMCEnergy_max;
        double          fMCEnergy_index;
        bool            fFillMCHistograms;
        bool            fgetXoff_Yoff_afterCut;
                bool            fWriteEventdatatrees;
        
        string          fCoreScatterMode;
        double          fCoreScatterRadius;
        
        double          fViewcone_min;
        double          fViewcone_max;
        
        bool            fAzimuthBins;
        bool            fIsotropicArrivalDirections;
        float           fIgnoreFractionOfEvents;
        
        bool            fTelescopeTypeCuts;
        
        string          fdatafile;
        string          fMCdatafile_tree;
        string          fMCdatafile_histo;
        string          fGammaHadronProbabilityFile;
        
        double          fze;
        int             fnoise;
        double          fpedvar;
        double          fXoff;
        double          fYoff;
        vector< double > fAzMin;
        vector< double > fAzMax;
        
        double          fWobbleIsotropic;
        
        unsigned int    telconfig_ntel;
        double          telconfig_arraycentre_X;
        double          telconfig_arraycentre_Y;
        double          telconfig_arraymax;
        
        string          fCREnergySpectrumFile;
        unsigned int    fCREnergySpectrumID;
        TF1*            fCREnergySpectrum;
        
        
        VInstrumentResponseFunctionRunParameter();
        ~VInstrumentResponseFunctionRunParameter() {}
        
        void                  print();
        VMonteCarloRunHeader* readMCRunHeader();
        bool                  readRunParameterFromTextFile( string iFile );
        bool                  testRunparameters();
        
        ClassDef( VInstrumentResponseFunctionRunParameter, 17 );
};

#endif
