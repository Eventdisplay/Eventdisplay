//!< VEffectiveAreaCalculator calculate effective areas and energy spectra

#ifndef VEffectiveAreaCalculator_H
#define VEffectiveAreaCalculator_H

#include "CData.h"
#include "VGammaHadronCuts.h"
#include "VAnaSumRunParameter.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "TEfficiency.h"
#include "VHistogramUtilities.h"
#include "VInstrumentResponseFunctionRunParameter.h"
#include "VStatistics.h"
#include "VSpectralWeight.h"
#include "VUtilities.h"

#include "TChain.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VEffectiveAreaCalculator
{
    private:
    
        vector< vector< double > > fEffArea_time;           // effective area vs time [time bin][energy bin]
	vector< vector< double > > fEffAreaMC_time;           // effective area vs time [time bin][energy bin] (MC Energy)
        vector< double > timebins;                          // centre of time bin [time bin]
        
        float fMC_ScatterArea;                                  // scatter area in MC (read from CORSIKA run header)
        
        bool bNOFILE;
        TDirectory* fGDirectory;
        
        vector< double > fZe;
        vector< double > fMCZe;
        vector< vector< double > > fEff_WobbleOffsets;
        vector< vector< vector< double > > > fEff_Noise;
        vector< vector< vector< vector< double > > > > fEff_SpectralIndex;
        
        // effective areas (reading of effective areas)
        vector< float >                      fEff_E0;
        map< unsigned int, vector< float > > fEffArea_map;
	map< unsigned int, vector< float > > fEffAreaMC_map;
        map< unsigned int, vector< float > > fEff_EsysMCRelative;
	map< unsigned int, vector< float > > fe_MC_Res_map;
	map< unsigned int, vector< float > > fe_Rec_Res_map;
	map< unsigned int, vector< float > > fe_Rec_Res_Err_map;
	map< unsigned int, unsigned int > fEntry_map;
        
        // mean and mean time binned effective areas
        unsigned int     fNTimeBinnedMeanEffectiveArea;
        vector< double > fVTimeBinnedMeanEffectiveArea;
	unsigned int     fNTimeBinnedMeanEffectiveAreaMC;
	vector< double > fVTimeBinnedMeanEffectiveAreaMC;
        TGraphAsymmErrors* gMeanEffectiveArea;
        TGraph2DErrors*    gTimeBinnedMeanEffectiveArea;
        
        TGraphErrors* gMeanSystematicErrorGraph;
        
	// For Binned Likelihood
	TGraphAsymmErrors* gMeanEffectiveAreaMC;
	TH2D*			   hMeanResponseMatrix;
        // unique event counting
        map< unsigned int, unsigned short int> fUniqueEventCounter;
        
        vector< double > fAreaRadius;
        vector< string > fScatterMode;
        vector< double > fXWobble;                //!< wobble offset in camera coordinates (grisudet)
        vector< double > fYWobble;                //!< wobble offset in camera coordinates (grisudet)
        vector< int >    fNoise;
        vector< double > fPedVar;
        
        double fEnergyAxis_minimum_defaultValue;
        double fEnergyAxis_maximum_defaultValue;
        
        VInstrumentResponseFunctionRunParameter* fRunPara;
        
        // gamma hadron cuts
        VGammaHadronCuts* fCuts;
        bool fIgnoreEnergyReconstruction;
        bool fIsotropicArrivalDirections;
        bool fTelescopeTypeCutsSet;
        
        // effective area calculation
        vector< double > fVMinAz;
        vector< double > fVMaxAz;
        // spectral weighting
        vector< double > fVSpectralIndex;
        VSpectralWeight* fSpectralWeight;
        
        // list of histograms
        // (yes, it is a mess!)
        
        TList* hisVList;
        vector< vector< TH1D* > > hVEmc;
        vector< vector< TH1D* > > hVEcut;
        vector< vector< TH1D* > > hVEcutRec;
        vector< vector< TH1D* > > hVEcutUW;
        vector< vector< TH1D* > > hVEcutRecUW;
        vector< vector< TH1D* > > hVEcutNoTh2;
        vector< vector< TH1D* > > hVEcutRecNoTh2;
        vector< vector< TProfile* > > hVEmcSWeight;
        vector< vector< TH1D* > > hVEcut500;
        vector< vector< TProfile* > > hVEsysRec;
        vector< vector< TProfile* > > hVEsysMC;
        vector< vector< TProfile* > > hVEsysMCRelative;
        vector< vector< TH2D* > > hVEsysMCRelativeRMS;
        vector< vector< TH2D* > > hVEsysMCRelative2D;
        vector< vector< TH2D* > > hVEsysMCRelative2DNoDirectionCut;
        vector< vector< TH2D* > > hVEsys2D;
        vector< vector< TH2D* > > hVResponseMatrix;
        vector< vector< TH2D* > > hVResponseMatrixFine;
        vector< vector< TH2D* > > hVResponseMatrixQC;
        vector< vector< TH2D* > > hVResponseMatrixFineQC;
        vector< vector< TH2D* > > hVResponseMatrixNoDirectionCut;
        vector< vector< TH2D* > > hVResponseMatrixFineNoDirectionCut;
        vector< vector< TH1D* > > hVWeightedRate;
        vector< vector< TH1D* > > hVWeightedRate005;
        vector< vector< vector< TH1D* > > > hVEcutSub; // [index][ecut_index][az]
        vector< vector< TH2D* > > hVAngErec2D;            // direction reconstruction
        vector< vector< TH2D* > > hVAngMC2D;            // direction reconstruction
        
        // angular resolution graphs (vector in az)
        vector< TGraphErrors* > fGraph_AngularResolution68p;
        vector< TGraphErrors* > fGraph_AngularResolution80p;
        vector< TGraphErrors* > fGraph_AngularResolutionKingSigma;
        vector< TGraphErrors* > fGraph_AngularResolutionKingGamma;
        vector< TH2D* >         hVAngularDiff_2D;
        vector< TH2D* >         hVAngularDiffEmc_2D;
        vector< TH2D* >         hVAngularLogDiff_2D;
        vector< TH2D* >         hVAngularLogDiffEmc_2D;
        
        // the following histograms are written to the output file
        // (into the fEffArea tree)
        TList* hisTreeList;
        TList* hisTreeListofHistograms;
        TH1D* hEmc;
        TH1D* hEcut;
        TH1D* hEcut500;            // fine binned effective areas for energy threshold determination
        TH1D* hEcutRec;
        TH1D* hEcutUW;
        TH1D* hEcutRecUW;
        TH1D* hEcutNoTh2;
        TH1D* hEcutRecNoTh2;
        TGraphAsymmErrors* gEffAreaMC;
        TGraphAsymmErrors* gEffAreaRec;
        TGraphAsymmErrors* gEffAreaNoTh2MC;
        TGraphAsymmErrors* gEffAreaNoTh2Rec;
        TProfile* hEmcSWeight;
        TProfile* hEsysRec;
        TProfile* hEsysMC;
        TProfile* hEsysMCRelative;
        TH2D* hEsysMCRelativeRMS;
        TH2D* hEsysMCRelative2D;
        TH2D* hEsysMCRelative2DNoDirectionCut;
        TH2D* hEsys2D;
        TH2D* hResponseMatrix;
        TH2D* hResponseMatrixFine;
        TH2D* hResponseMatrixQC;
        TH2D* hResponseMatrixFineQC;
        TH2D* hResponseMatrixNoDirectionCut;
        TH2D* hResponseMatrixFineNoDirectionCut;
        TH1D* hWeightedRate;
        TH1D* hWeightedRate005;
        vector< TH1D* > hEcutSub;                //! events after individual cuts
        TH2D *hAngularDiff_2D;
        TH2D *hAngularDiffEmc_2D;
        TH2D *hAngularLogDiff_2D;
        TH2D *hAngularLogDiffEmc_2D;
        
        int fEffectiveAreaVsEnergyMC;            // 0 = vs MC energy, 1 = vs rec energy (approx. method), 2 = vs rec energy (default)
        TTree* fEffArea;
        float ze;
        int fAzBin;                               //!< az bin: definitions see getEffectiveArea
        float fMinAz;
        float fMaxAz;
        float fXoff;
        float fYoff;
        float fWoff;
        float fSpectralIndex;
        int   fTNoise;
        float fTNoisePE;
        float fTPedvar;
        int   nbins;
        float e0[1000];
        float eff[1000];
	int nbins_MC;
	float e0_MC[1000];
	float eff_MC[1000];
        float seff_L[1000];
        float seff_U[1000];
        float eff_error[1000];
        float esys_rel[1000];
        int   Rec_nbins;
        float Rec_e0[1000];
        float Rec_eff[1000];
        float Rec_seff_L[1000];
        float Rec_seff_U[1000];
        float Rec_eff_error[1000];
        float Rec_angRes_p68[1000];
        float Rec_angRes_p80[1000];
        float Rec_angRes_kingSigma[1000];
        float Rec_angRes_kingGamma[1000];
        int nbins_MC_Res;
        float e_MC_Res[1000];
        float e_Rec_Res[1000];
        float e_Rec_Res_Err[1000];
        
        TTree* fAcceptance_AfterCuts_tree;       //Information for all the events after cuts to construct the background map
        double fXoff_aC;
        double fYoff_aC;
        double fXoff_derot_aC;
        double fYoff_derot_aC;
        double fErec;
        double fEMC;
        double fCRweight;                         // #/s/sr (the right unit for the ctools acceptance map) This normalise the map to the CR spectrum
        // Needs option ESPECTRUM_FOR_WEIGHTING to be turned on, which only make sense for CR
        bool fsolid_angle_norm_done;
        double fsolid_angle_norm;                   // solid angle normalisation needed for the CRweight filled in
        // fAcceptance_AfterCuts_tree (for the histogram it is done later in VSensitivityCalculator)
        void Calculate_Bck_solid_angle_norm();
        
        /////////////////////////////
        // event data
        TTree *fEventTreeCuts;
        int fCut_Class;
        float fCut_MVA;

        // effective area smoothing
        int fSmoothIter;
        double fSmoothThreshold;
        
        // mean values from getEffectiveAreas
        double fEffectiveAreas_meanZe;
        double fEffectiveAreas_meanWoff;
        double fEffectiveAreas_meanPedVar;
        double fEffectiveAreas_meanIndex;
        double fEffectiveAreas_meanN;
        
        // Gaussian function for approximating the response matrix
        TF1 *fGauss;
        // Bool to handle if likelihood analysis is required
        bool bLikelihoodAnalysis;
        bool bIsOn;

        TGraphAsymmErrors* applyResponseMatrix( TH2* h, TGraphAsymmErrors* g );
        bool               binomialDivide( TGraphAsymmErrors* g, TH1D* hrec, TH1D* hmc );
        void               copyProfileHistograms( TProfile*,  TProfile* );
        void               copyHistograms( TH1*,  TH1*, bool );
        void               fillAngularResolution( unsigned int i_az, bool iContaintment_80p );
        void               fillEventDataTree( int iCutClass, float iMVA );
        void               fillEcutSub( double iE, unsigned int iIndex );
        float              getAzMean( float azmin, float azmax );
        double             getCRWeight( double iEMC_TeV_lin, TH1* h , bool for_back_map = false, TH1* hF = 0 );
        double             getEffectiveAreasFromHistograms( double erec, double ze, double woff, double iPedVar,
                double iSpectralIndex, bool bAddtoMeanEffectiveArea = true );
        bool               getMonteCarloSpectra( VEffectiveAreaCalculatorMCHistograms* );
        double             getMCSolidAngleNormalization();
        vector< unsigned int > getUpperLowBins( vector< double > i_values, double d );
        bool   initializeEffectiveAreasFromHistograms( TTree*, TH1D*, double azmin, double azmax, double ispectralindex, double ipedvar );
        vector< float >    interpolate_effectiveArea( double iV, double iVLower, double iVupper,
               vector< float > iEL, vector< float > iEU, bool iCos = true );
        void               multiplyByScatterArea( TGraphAsymmErrors* g );
        void               reset();
        void               resetTimeBin();
        void               resetHistograms( unsigned int iZe );
        void               resetHistogramsVectors();
        void               smoothEffectiveAreas( map< unsigned int, vector< float > > );
        bool               testAzimuthInterval( CData* d, double iZe, double iMinAz, double iMaxAz );
        
    public:
    
        VEffectiveAreaCalculator( string ieffFile, double azmin, double azmax, double iPedVar, double iIndex,
                                  vector< double> fMCZe, int iSmoothIter = -1, double iSmoothThreshold = 1.,
				  int iEffectiveAreaVsEnergyMC = 2, bool iLikelihoodAnalysis = false, bool iIsOn = false );          // constructor for reading
        VEffectiveAreaCalculator( VInstrumentResponseFunctionRunParameter*, VGammaHadronCuts* );       // constructor for filling
        ~VEffectiveAreaCalculator();
        
        void               cleanup();
        bool               fill( CData* d, VEffectiveAreaCalculatorMCHistograms* iMC_histo, unsigned int iMethod );
        TH1D*              getHistogramhEmc();
        TGraphErrors*      getMeanSystematicErrorHistogram();
        TTree*             getEffectiveAreaTree()
        {
            return fEffArea;
        }
                TTree*             getEventCutDataTree()
                {
                        return fEventTreeCuts;
                }
        TTree*             getAcceptance_AfterCuts()
        {
            return fAcceptance_AfterCuts_tree;
        }
        double             getEnergyAxis_minimum_defaultValue()
        {
            return fEnergyAxis_minimum_defaultValue;
        }
        double             getEnergyAxis_maximum_defaultValue()
        {
            return fEnergyAxis_maximum_defaultValue;
        }
        double             getEffectiveArea( double erec, double ze, double iWoff, double iPedvar, double iSpectralIndex = -2.5,
                                             bool bAddtoMeanEffectiveArea = true );
        TGraphAsymmErrors*  getMeanEffectiveArea();
        TGraph2DErrors*     getTimeBinnedMeanEffectiveArea();
        TGraphAsymmErrors*  getMeanEffectiveAreaMC();

        TH1D*               getMCHistogram()
        {
            return hEmc;
        }
        vector< TH1D* >     initializeHistogramsVectorH1D( TH1D* h, string iName, unsigned int i );
        vector< TH2D* >     initializeHistogramsVectorH2D( TH2D* h, string iName, unsigned int i );
        vector< TProfile* > initializeHistogramsVectorHProfile( TProfile* h, string iName, unsigned int i );
        void                initializeHistograms( vector< double > iAzMin, vector< double > iAzMax, vector< double > iSpectralIndex );
        void                setTimeBinnedMeanEffectiveArea( double i_time );
        void                setTimeBinnedMeanEffectiveAreaMC( double i_time );
        void                setAngularResolution2D( unsigned int i_az, vector< TH2D* > );
        void                setAngularResolutionGraph( unsigned int i_az, TGraphErrors* g, bool iAngContainment_80p );
        void                setAngularResolutionKingSigmaGraph( unsigned int i_az, TGraphErrors* g );
        void                setAngularResolutionKingGammaGraph( unsigned int i_az, TGraphErrors* g );
        void                setAzimuthCut( int iAzBin, double iAzMin, double iAzMax );
        void                setEffectiveArea( int iMC )
        {
            fEffectiveAreaVsEnergyMC = iMC;
        }
        void                setIgnoreEnergyReconstructionCuts( bool iB = false )
        {
            fIgnoreEnergyReconstruction = iB;
        }
        void                setIsotropicArrivalDirections( bool iB = false )
        {
            fIsotropicArrivalDirections = iB;
        }
        bool                setMonteCarloEnergyRange( double iMin, double iMax, double iMCIndex = 2. );
        void                setNoiseLevel( int iN, double iP );
        void                setTelescopeTypeCuts( bool iB = true )
        {
            fTelescopeTypeCutsSet = iB;
        }
        void                setWobbleOffset( double x, double y );
       
        void 				addMeanResponseMatrix( vector <float> i_emc, vector <float> i_erec , vector <float> i_erec_err );
        TH2D* 				getMeanResponseMatrix()
        {
            VHistogramUtilities::normalizeTH2D_y(hMeanResponseMatrix);
            return (TH2D*)hMeanResponseMatrix->Clone();
        }
};
#endif
