//! VDataMCComparision

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"

#include "CData.h"
#include "VEvndispRunParameter.h"
#include "VGammaHadronCuts.h"
#include "VMonteCarloRunHeader.h"
#include "VSpectralWeight.h"
#include "VUtilities.h"

#include <bitset>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VDataMCComparisionHistogramData
{
    public:
    
        string fVarName;
        string fHistogramType;
        unsigned int fTelescopeID;                    // 0 = array variable
        
        TH1D*  fHis1D;
        TH2D*  fHis2D_Erec;
        TH2D*  fHis2D_ntubes;
        TH2D*  fHis2D_size;
        TH2D*  fHis2D_sizeHG;
        TH2D*  fHis2D_sizeLG;
        
        VDataMCComparisionHistogramData( string iVarName = "", string iHistogramType = "", unsigned int iTelescopeID = 0 );
        ~VDataMCComparisionHistogramData() {}
        bool   initHistogram( string iXTitle, int iNbins, double ix_min, double ix_max );
        TH2D*  newHistogram( string iName, string iYTitle, string iXTitle, int iNbins, double ix_min, double ix_max, double iy_min, double iy_max );
        void   fill( double iV, double iWeight = 1., double iLogEnergy_TeV = -99., int i_ntubes = -99,
                     double i_sizeLog10 = -99., double i_sizeHGLog10 = -99. );
};

class VDataMCComparision
{
    private:
    
        enum E_varname { ELENGTH, EWIDTH, EDIST, EALPHA, ENTUBES, ENLOWGAIN, ESIZE, ESIZEHG,
                         ESIZELG, EFRACLOW, EMAX1, EMAX2, EMAX3, ELOSS, ELOS, EASYM,
                         ECENX, ECENY, ETGRADX, EMSCWT, EMSCLT, EMWRT, EMLTT,
                         ETELDIST, ETHETA2, ELTHETA2,
                         EMSCW, EMSCL, EMWR, EMLR, EXCORE, EYCORE, EEREC, ENIMAGES, ERECRAT, EPEDVAR, EPEDVART,
                         EAEL, EAAZ,
                         EIMGSEL, EEMISSIONHEIGHT, EMVA,
                         ESIGMAT3D, ENC3D, ESMAX3D, EERRSIGMAT3D, EOMEGA3D, EDEPTH3D, ERWIDTH3D, EERRRWIDTH3D
                       };
                       
        string fName;
        int fNTel;
        
        vector< double > fTel_x;
        vector< double > fTel_y;
        vector< double > fTel_z;
        
        // wobble north offset
        double fWobbleNorth;
        double fWobbleEast;
        bool   fWobbleFromDataTree;
        
        double fAzMin;
        double fAzMax;
        bool fAzRange;
        double fZeMin;
        double fZeMax;
        
        // spectral weighting
        VSpectralWeight* fSpectralWeight;
        
        // data tree
        CData* fData;
        
        // cuts
        VGammaHadronCuts* fCuts;
        bool fCalculateMVAValues;
        
        // lists with all histograms
        TList* hisList;
        vector<TH1D* > hTel;
        vector<TH2D* > hTel2D;
        
        // histogram classes
        map< E_varname, vector< VDataMCComparisionHistogramData* > > fHistoSingleTel;
        map< E_varname, VDataMCComparisionHistogramData* > fHistoArray;
        
        // stereo histograms
        TH2D* hXYcore;
        TH2D* hAzYcore;
        TH2D* hYt2;
        vector<TH2D* > hcen_xy;
        vector< TH2D* > hdistR;
        
        // histogram for azimuth weighting
        TH1D* hAzWeight;
        
        // angle for shower max correction
        double fShowerMaxZe_deg;
        
        void setEntries( TH1D* );
        void setEntries( TH2D* );
        
        double getCorrectedEmissionHeight( double iEM, double iZe );
        void initialGammaHadronCuts();
        
    public:
    
        VDataMCComparision( string, int, bool );
        ~VDataMCComparision() {}
        void defineHistograms();
        bool fillHistograms( string ifile, int iSingleTelescopeCuts );
        TH1D* getAzimuthWeightingHistogram( string ifile );
        void resetTelescopeCoordinates();
        void scaleHistograms( string );
        void setAzimuthWeightingHistogram( TH1D* hAz )
        {
            hAzWeight = hAz;
        }
        void setAzRange( double iAzMin, double iAzMax );
        void setZeRange( double iZeMin, double iZeMax );
        bool setOnOffHistograms( VDataMCComparision*, VDataMCComparision*, double norm );
        void setShowerMaximZe_deg( double iZe = 20. )
        {
            fShowerMaxZe_deg = iZe;
        }
        bool setTelescopeCoordinates( double x, double y, double z = 0. );
        void setWobbleFromDataTree()
        {
            fWobbleFromDataTree = true;
        }
        bool writeHistograms( string iOutFile );
};
