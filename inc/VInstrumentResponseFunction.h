//! VInstrumentResponseFunction calculate response function (e.g. angular resolution)

#ifndef VInstrumentResponseFunction_H
#define VInstrumentResponseFunction_H

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <vector>

#include "CData.h"
#include "VGammaHadronCuts.h"
#include "VInstrumentResponseFunctionData.h"
#include "VInstrumentResponseFunctionRunParameter.h"
#include "VSpectralWeight.h"

using namespace std;

class VInstrumentResponseFunction
{
    private:
    
        bool  fDebug;
        
        string fName;
        string fType;
        VInstrumentResponseFunctionRunParameter* fRunPara;
        
        // output file
        TFile* fOutputFile;
        
        // data tree
        CData*   fData;
        
        // return data tree
        TTree*    fDataProduct;
        VInstrumentResponseFunctionData* fIRFData_Tree;
        
        // histograms are not re-filled but duplicated
        unsigned int fDuplicationID;
        
        unsigned int fEnergyReconstructionMethod;
        
        // histograms and data
        vector< vector< VInstrumentResponseFunctionData* > > fIRFData;
        
        // cuts
        vector< VGammaHadronCuts* > fAnaCuts;
        bool              fTelescopeTypeCutsSet;
        
        // effective area calculation
        vector< double > fVMinAz;
        vector< double > fVMaxAz;
        // spectral weighting
        vector< double > fVSpectralIndex;
        VSpectralWeight* fSpectralWeight;
        
        // containment probabilities
        double  fContainmentProbability;
        double  fContainmentProbabilityError;
        
        bool    defineHistograms();
        bool    fillEventData();
        VGammaHadronCuts* getGammaHadronCuts( CData* );
        
    public:
    
        VInstrumentResponseFunction();
        ~VInstrumentResponseFunction();
        bool   doNotDuplicateIRFs()
        {
            if( fDuplicationID == 9999 )
            {
                return true;
            }
            return false;
        }
        bool   fill();
        bool   fillResolutionGraphs( vector< vector< VInstrumentResponseFunctionData* > > iIRFData );
        double getContainmentProbability()
        {
            return fContainmentProbability;
        }
        TTree* getDataProduct()
        {
            return fDataProduct;
        }
        vector< TH2D* > getAngularResolution2D( unsigned int iAzBin, unsigned int iSpectralIndexBin );
        TGraphErrors* getAngularResolutionGraph( unsigned int iAzBin, unsigned int iSpectralIndexBin );
        unsigned int getDuplicationID()
        {
            return fDuplicationID;
        }
        vector< vector< VInstrumentResponseFunctionData* > > getIRFData()
        {
            return fIRFData;
        }
        string getName()
        {
            return fName;
        }
        string getResolutionType()
        {
            return fType;
        }
        bool   initialize( string iName, string iType, unsigned int iNTel, double iMCMaxCoreRadius,
                           double iZe, int iNoise, double iPedvars, double iXoff, double iYoff );
        void   printCutStatistics();
        void   setDuplicationID( unsigned int iDuplicationID = 9999 );
        void   setEnergyReconstructionMethod( unsigned int iMethod );
        void   setCuts( vector< VGammaHadronCuts* > iCuts );
        void   setContainmentProbability( double iP = 0.68, double iPError = 0.95 )
        {
            fContainmentProbability = iP;
            fContainmentProbabilityError = iPError;
        }
        void   setDataTree( CData* iData );
        void   setMonteCarloEnergyRange( double iMin, double iMax, double iMCIndex = 2. );
        void   setOutputFile( TFile* iOutputFile )
        {
            fOutputFile = iOutputFile;
        }
        void   setRunParameter( VInstrumentResponseFunctionRunParameter* iRunPar = 0 );
        void   setTelescopeTypeCuts( bool iB = true )
        {
            fTelescopeTypeCutsSet = iB;
        }
};

#endif
