//! VImageCleaning  image cleaning routines

#ifndef VIMAGECLEANING_H
#define VIMAGECLEANING_H

#include <VEvndispData.h>
#include <TGraphErrors.h>

#include "VImageCleaningRunParameter.h"

using namespace std;

class VImageCleaning
{
    private:
    
        VEvndispData* fData;
        bool fWriteGraphToFileRecreate;
        
        void cleanImageWithTiming( VImageCleaningRunParameter* iImageCleaningParameters, bool isFixed );
        void fillImageBorderNeighbours();
        void fillBorderBorderNeighbours();
        void recoverImagePixelNearDeadPixel();
        void printDataError( string iFunctionName );
        
        // cluster cleaning
        void mergeClusters();
        void removeSmallClusters( int );
        void addToCluster( unsigned int cID, unsigned int iChan );
        void removeCluster( unsigned int cID ) ;
        vector<int> fNpixCluster;
        vector<double> fSizeCluster;
        
        // NN image cleaning
        bool  kInitNNImageCleaning;
        bool  kInitNNImgClnPerTelType[VDST_MAXTELTYPES];
        const static unsigned int fIPRdim = 200;
        bool  NNoptNoTimeing;
        TObjArray* fProb4nnCurves;
        TObjArray* fProb3nnrelCurves;
        TObjArray* fProb2plus1Curves;
        TObjArray* fProb2nnCurves;
        TObjArray* fProbBoundCurves;
        TObjArray* fIPRgraphs;
        vector< vector< bool > > fifActiveNN;                      // [nteltypes][nngroups]
        bool ifActiveNN[VDST_MAXNNGROUPTYPES][VDST_MAXTELTYPES];   // if  NN groups is searched in NN-image cleaning procedure
        int   VALIDITY[VDST_MAXCHANNELS];      //   Flags for pixels, accepted by nn-image cleaning. VALIDITY[i]=2-6 : core pixels, VALIDITY[i]>6 :boundary pixels
        int   VALIDITYBOUND[VDST_MAXCHANNELS]; //
        int   VALIDITYBUF[VDST_MAXCHANNELS];   //
        unsigned int nRings;                   //
        float CoincWinLimit;                   // ns (maximum coicidence for the NN group)
        float fMinRate[VDST_MAXTELTYPES];      //
        double fFakeImageProb;
        bool   setExplicitSampleTimeSlice;      // Set the sample time slice and number of ADC bins to read explicitly
        float  sampleTimeSlice;                // Size of time slice in ns (usually 1 or 2 ns)
        unsigned int  nBinsADC;                // Number of ADC bins summed up, each bin the size of sampleTimeSlice
        
        float INTENSITY[VDST_MAXCHANNELS];     //
        float TIMES[VDST_MAXCHANNELS];         //
        float** IPR;                           // IPR[TelType][ScanDim] scan. Not used TelType==0 is filled with DT values
        
        int   LocMin( int n, float* ptr, float& min );
        int   LocMax( int n, float* ptr, float& max );
        
        // main functions
        bool  BoundarySearch( unsigned int TrigSimTelType, float thresh, TF1* fProbCurve, float refdT, int refvalidity, int idx );
        unsigned int   NNGroupSearchProbCurve( unsigned int TrigSimTelType, TF1* fProbCurve, float PreCut );
        unsigned int   NNGroupSearchProbCurveRelaxed( unsigned int TrigSimTelType, TF1* fProbCurve, float PreCut );
        bool  NNChargeAndTimeCut( TGraph* IPR, TF1* fProbCurve, float charge, float dT,
                                  float iCoincWinLimit, bool bInvert = false );
        void  ScaleCombFactors( unsigned int TrigSimTelType, float scale );
        void  ResetCombFactors( unsigned int TrigSimTelType );
        int   ImageCleaningCharge( unsigned int TrigSimTelType );
        bool  InitNNImageCleaning();
        bool  InitNNImgClnPerTelType( unsigned int TrigSimTelType );
        void  DiscardTimeOutlayers( unsigned int TrigSimTelType );
        void  DiscardLocalTimeOutlayers( float NNthresh[6] ); // use this function
        void  DiscardIsolatedPixels();
        void  FillIPR( unsigned int TrigSimTelType );
        void  FillPreThresholds( TGraph* gipr, float NNthresh[6] ); // defines pre-search thresholds for nn-groups (below this threshold group is not searched)
        TGraphErrors* GetIPRGraph( unsigned int TrigSimTelType, float ScanWidow );
        void  SetNeighborRings( unsigned short* VALIDITYBOUNDBUF, float* TIMESReSearch, float* REFTHRESH );
        
        
        TF1* defineRateContourFunction( unsigned int type, TString funcname, float iRate, float iNfold, float iCombFactor,
                                        float xlow, float xup );
        TF1* defineRateContourBoundFunction( unsigned int type, TString funcname, float iRate, float iNfold, float iCombFactor,
                                             float xlow, float xup );
        void writeProbabilityCurve( TGraph* iIPR, TF1* iProb, double iRate );
        
    public:
    
        VImageCleaning( VEvndispData* iData = 0 );
        ~VImageCleaning() {}
        
        // tailcut cleaning
        void cleanImageFixed( VImageCleaningRunParameter* iImageCleaningParameters );
        void cleanImageFixed( double hithresh, double lothresh, double brightthresh = -999. );
        void cleanImagePedvars( VImageCleaningRunParameter* iImageCleaningParameters );
        
        // time tailcut cleaning
        void cleanImagePedvarsTimeDiff( VImageCleaningRunParameter* iImageCleaningParameters );
        
        // time cluster cleaning
        void cleanImageFixedWithTiming( VImageCleaningRunParameter* iImageCleaningParameters );
        void cleanImagePedvarsWithTiming( VImageCleaningRunParameter* iImageCleaningParameters );
        
        // cluster cleaning
        void cleanImageWithClusters( VImageCleaningRunParameter* iImageCleaningParameters, bool isFixed );
        
        // trace correlation cleaning
        void cleanImageTraceCorrelate( VImageCleaningRunParameter* iImageCleaningParameters );
        
        // Optimized NN image cleaning
        void  cleanNNImageFixed( VImageCleaningRunParameter* iImageCleaningParameters );
        
        
        void addImageChannel( unsigned int iChannel );                      // add this pixel to image
        void removeImageChannel( unsigned int iChannel );                   // remove this pixel from image
        void resetImageChannel( unsigned int iChannel );                    // reset this pixel to standard value
        
};
#endif
