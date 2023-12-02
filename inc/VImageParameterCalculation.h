//! VImageParameterCalculation   calculation of Hillas Parameters

#ifndef VImageParameterCalculation_H
#define VImageParameterCalculation_H

#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>

#include <VDetectorGeometry.h>
#include <VEvndispData.h>
#include <VHoughTransform.h>
#include <VImageParameter.h>
#include <VImageParameterFitter.h>

#include "TError.h"
#include "TMath.h"
#include "TObject.h"
#include <TH2.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>

//Include the GSL elliptical integral header file if GSL is installed
//This is for the impact parameter corrected muon size calculation
//#ifndef NOGSL
//#include "gsl/gsl_sf_ellint.h"
//#endif

using namespace std;

class VImageParameterCalculation
{
    private:
    
        bool fDebug;
        
        VDetectorGeometry* fDetectorGeometry;     //< the detector geometry
        VImageParameter* fParGeo;                 //!< image parameters (geo.)
        bool fboolCalcGeo;                        //!< switch to show if there was a geometrical analysis
        bool fboolCalcTiming;                     //!< switch to show if there was a timing analysis
        VImageParameter* fParLL;                  //!< image parameters (log lik.)
        VImageParameterFitter* fImageFitter;
        
        VEvndispData* fData;
        
        double getFractionOfImageBorderPixelUnderImage( double, double, double, double, double, double );
        void   setImageBorderPixelPosition( VImageParameter* iPar );
        
        // Hough transform
        VHoughTransform* fHoughTransform;
        
    public:
    
        VImageParameterCalculation( unsigned int iShortTree = 0, VEvndispData* iData = 0 );
        ~VImageParameterCalculation();
        vector<bool> calcLL( bool iUseSums2 = false,
                             bool i_reInitializeLL = false,
                             bool iEqualSummationWindows = false );          //!< calculate image parameters (log like)
        void muonRingFinder();                     //!< fit a single ring to the image to look for muons
        void sizeInMuonRing();                     //! calculate the brightness of the muon ring
        void muonPixelDistribution();              //!< determine the distribution of pixels in the image
        
        //Impact parameter correction factor for size
        float correctSizeInMuonRing();
        
        //Hough transform
        void houghInitialization();
        void houghMuonPixelDistribution();          //!< determine the distribution of pixels in the image
        
        void calcTriggerParameters( vector<bool> fTrigger );                                   //!< calculate trigger-level image parameters
        void calcParameters();                                                                 //!< calculate image parameters (geo.)
        void calcTimingParameters( bool iIsSecondPass );
        VImageParameter* getParameters()              //!< get image parameters
        {
            return fParGeo;
        }
        VImageParameter* getLLParameters()            //!< get image parameters from loglikelihood
        {
            return fParLL;
        }
        bool getboolCalcGeo()                     //!< get image parameters calculated flag
        {
            return fboolCalcGeo;
        }
        bool getboolCalcTiming()                  //!< get timing parameters calculated flag
        {
            return fboolCalcTiming;
        }
        VDetectorGeometry* getDetectorGeo()       //!< get the detector geometry
        {
            return fDetectorGeometry;
        }
        VDetectorGeometry* getDetectorGeometry()  //!< get the detector geometry
        {
            return fDetectorGeometry;
        }
        //!< fill image/border list (optional)
        void fillImageBorderPixelTree();
        //!< initialize minuit
        void initMinuit( int );
        //!< set the detector geometry
        void setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        void setDetectorGeometry( VDetectorGeometry* iDet )
        {
            fDetectorGeometry = iDet;
        }
        void setParameters( VImageParameter* iP )     //!< set pointer to parameter class
        {
            fParGeo = iP;
        }
        void setParametersLogL( VImageParameter* iP ) //!< set pointer to parameter class
        {
            fParLL = iP;
        }
};
#endif
