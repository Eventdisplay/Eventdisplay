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

#define VMAXBINS 1000

using namespace std;

class VEffectiveAreaCalculator
{
    private:
    
        vector< vector< double > > fEffArea_time;
     	vector< vector< double > > fEffAreaMC_time;
        vector< double > timebins;
        
        float fMC_ScatterArea;
        
        bool bNOFILE;
        TDirectory* fGDirectory;
        
        vector< double > fZe;
        vector< double > fMCZe;
        vector< vector< double > > fEff_WobbleOffsets;
        vector< vector< vector< double > > > fEff_Noise;
        vector< vector< vector< vector< double > > > > fEff_SpectralIndex;
        
        // effective areas (reading of effective areas)
		unsigned int fNBins;                    // bins in the true energy of MC (fEff_E0)
        unsigned int fBiasBin;                  // bins in the energy bias
        unsigned int fhistoNEbins;              // energy bins for histograms only
        unsigned int fResponseMatricesEbinning; // fine bins for response matrices. Likelihood analysis.
        unsigned int fLogAngularBin;            // bins for the log10(angular diff R,MC [deg])

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

        double fLogAngular_minimum_defaultValue;
        double fLogAngular_maximum_defaultValue;
        
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
        enum E_HISTYPE { E_1D, E_1P, E_2D };
        enum E_HIS1D { E_Emc, E_Ecut, E_EcutUW, E_EcutNoTh2,
                E_Ecut500, E_EcutRec, E_EcutRecUW, E_EcutRecNoTh2,
                E_WeightedRate, E_WeightedRate005,
                E_EcutTrigger, E_EcutFiducialArea, 
                E_EcutStereoQuality, E_EcutTelType,
                E_EcutDirection, E_EcutEnergyReconstruction,
                E_EcutGammaHadron, E_EmcUW };
        enum E_HIS1P { E_EmcSWeight, E_EsysMCRelative };
        enum E_HIS2D {E_EsysMCRelativeRMS, E_EsysMCRelative2D, 
               E_EsysMCRelative2DNoDirectionCut, E_Esys2D,
               E_ResponseMatrix, E_ResponseMatrixFine,
               E_ResponseMatrixQC, E_ResponseMatrixFineQC,
               E_ResponseMatrixNoDirectionCut,
               E_ResponseMatrixFineNoDirectionCut};
        
        TList* hisVList;
        // dimension: [spectral index][az]
        map< int, vector< vector< TH1D* > > > hV_HIS1D;
        map< int, vector< vector< TProfile* > > > hV_HIS1P;
        map< int, vector< vector< TH2D* > > > hV_HIS2D;

        // angular resolution graphs (vector in az)
        vector< TGraphErrors* > fGraph_AngularResolution68p;
        vector< TGraphErrors* > fGraph_AngularResolution80p;
        vector< TGraphErrors* > fGraph_AngularResolutionKingSigma;
        vector< TGraphErrors* > fGraph_AngularResolutionKingGamma;

        vector< TH2D* >         hVAngularLogDiffEmc_2D;
        TH2D *hAngularLogDiffEmc_2D;
        
        // the following histograms are written to the output file
        // (into the fEffArea tree)
        TList* hisTreeList;
        TList* hisTreeListofHistograms;
        map< int, TH1D* > h_HIS1D;
        map< int, TProfile* > h_HIS1P;
        map< int, TH2D* > h_HIS2D;

        TGraphAsymmErrors* gEffAreaMC;
        TGraphAsymmErrors* gEffAreaRec;
        TGraphAsymmErrors* gEffAreaNoTh2MC;
        TGraphAsymmErrors* gEffAreaNoTh2Rec;

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
        float e0[VMAXBINS];
        float eff[VMAXBINS];
        int nbins_MC;
        float e0_MC[VMAXBINS];
        float eff_MC[VMAXBINS];
        float eff_error[VMAXBINS];
        float effNoTh2[VMAXBINS];
        float effNoTh2_error[VMAXBINS];
        float esys_rel[VMAXBINS];
        float Rec_eff[VMAXBINS];
        float Rec_eff_error[VMAXBINS];
        float Rec_effNoTh2[VMAXBINS];
        float Rec_effNoTh2_error[VMAXBINS];
        float Rec_angRes_p68[VMAXBINS];
        float Rec_angRes_p80[VMAXBINS];
        float Rec_angRes_kingSigma[VMAXBINS];
        float Rec_angRes_kingGamma[VMAXBINS];
        int nbins_MC_Res;
        float e_MC_Res[VMAXBINS];
        float e_Rec_Res[VMAXBINS];
        float e_Rec_Res_Err[VMAXBINS];
        
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
        bool               binomialDivide( TGraphAsymmErrors* g, TH1D* hrec, TH1D* hmc,
                               float *eff = 0, float* eff_error = 0 );
        void               copyProfileHistograms( TProfile*,  TProfile* );
        void               copyHistograms( TH1*,  TH1*, bool );
        void               fillAngularResolution( unsigned int i_az, bool iContaintment_80p );
        void               fillEventDataTree( int iCutClass, float iMVA );
        void               fillEcutSub( double iE, enum E_HIS1D iCutIndex );
        void               fillHistogram( int iHisType, int iHisN, 
                                          unsigned int s, unsigned i_az, 
                                          double i_x, double i_w, double i_w2D = 1. );
        float              getAzMean( float azmin, float azmax );
        double             getCRWeight( double iEMC_TeV_lin, TH1* h , bool for_back_map = false, TH1* hF = 0 );
        double             getEffectiveAreasFromHistograms( double erec, double ze, double woff, double iPedVar,
                double iSpectralIndex, bool bAddtoMeanEffectiveArea = true );
        string             getEffectiveAreaNamefromEnumInt( int i, string iType );
        bool               getMonteCarloSpectra( VEffectiveAreaCalculatorMCHistograms* );
        double             getMCSolidAngleNormalization();
        vector< unsigned int > getUpperLowBins( vector< double > i_values, double d );
        bool   initializeEffectiveAreasFromHistograms( TTree*, TH1D*, double azmin, double azmax, double ispectralindex, double ipedvar );
        vector< float >    interpolate_effectiveArea( double iV, double iVLower, double iVupper,
               vector< float > iEL, vector< float > iEU, bool iCos = true );
        bool               newEffectiveAreaHistogram( string iType, int iHisN, string iHisTitle, 
                                                      string iTitleX, string iTitleY, 
                                                      int i_nbins, double i_xmin, double i_xmax,
                                                      int i_nbins_y = -1,
                                                      double i_ymin = -1., double i_ymax = 1.,
                                                      string iPOpt = "" );
        void               reset();
        void               resetEffAreaArray( float *v );
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
        TGraphErrors*      getMeanSystematicErrorHistogram();
        TTree*             getEffectiveAreaTree()
        {
            return fEffArea;
        }
        int getEnergyAxis_nbins_defaultValue()
        {
                return fhistoNEbins;
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
            if( h_HIS1D.find( E_Emc ) != h_HIS1D.end() )
            {
                return h_HIS1D[E_Emc];
            }
            return 0;
        }
        TH1D*               getMCHistogramUnWeighted()
        {
            if( h_HIS1D.find( E_EmcUW ) != h_HIS1D.end() )
            {
                return h_HIS1D[E_EmcUW];
            }
            return 0;
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
