//! VImageParameterFitter

#ifndef VImageParameterFitter_H
#define VImageParameterFitter_H

#include <iostream>
#include <vector>

#include <VEvndispData.h>

#include "TF2.h"
#include "TMinuit.h"

using namespace std;

// global functions and pointers for fitting
extern void get_LL_imageParameter_2DGauss( Int_t&, Double_t*, Double_t&, Double_t*, Int_t );
extern void get_LL_imageParameter_2DGaussRotated( Int_t&, Double_t*, Double_t&, Double_t*, Int_t );
Double_t normal2DRotated( Double_t *x, Double_t *par );
extern TMinuit* fLLFitter;

class VImageParameterFitter : public TObject
{
   private:
      // data class
      VEvndispData* fData;
      VImageParameter* fParLL;
      VImageParameter* fParGeo;
      TF2 *fNormal2D;

      bool fDebug;
      bool fLLDebug;

      // tolerance value - usually not to be modified
      double ZeroTolerence;

      // image fitting function
      bool bRotatedNormalDistributionFit;

      //  starting values for fit
      double fdistXmin;
      double fdistXmax;
      double fdistYmin;
      double fdistYmax;

      vector<double> fll_X;                     //!< data vector for minuit function (x-coordinate of pmt)
      vector<double> fll_Y;                     //!< data vector for minuit function (y-coordinate of pmt)
      vector<double> fll_R;                     //!< tube radius of pmt
      vector<double> fll_Sums;                  //!< data vector for minuit function
      vector<bool> fLLEst;                      //!< true if channel has an estimated sum from the LL fit
      vector<double> fll_T;                     //!< data vector for minuit function (time)

      // fit parameter
      double rho;
      double drho;
      double cen_x;
      double dcen_x;
      double sigmaX;
      double dsigmaX;
      double cen_y;
      double dcen_y;
      double sigmaY;
      double dsigmaY;
      double signal;
      double dsignal;
      double toffset;
      double dtoffset;
      double tgrad;
      double dtgrad;
      double tchi2;
      double dtchi2;
      double phi;
      double dphi;
        
      //!< return value of 2d-gaus at channel iChannel
      double fill_pixel_sums( bool iUseSums2 );
      double getLL_startingvalue_rho();
      double getLL_paramameterlimits_rho( double, double );
      double getLL_paramameterlimits_cen( double, double, double, double );
      bool minimize_time_gradient()
      {
          if( fData && fData->getRunParameter() ) return fData->getRunParameter()->fMinimizeTimeGradient;
          return false;
      }
      double redang( double angle, double maxI );  //!< reduce angle to interval [0.,maxI]

      // fitting
      void defineFitParameters();
      void getFitStatistics();
      void getFitResults();
      void resetFitParameters();

      // image parameters
      void calculateImageParameters( bool iUseSums2, bool iEqualSummationWindows );
      void calculate_image_distance();
      void calculate_image_size( bool iUseSums2, bool iEqualSummationWindows );
      void calculate_image_length( double z, double dz2 );
      void calculate_image_phi( double dsxxy2 );
      void calculate_image_rho();
      void calculate_image_width( double z, double dz2 );
      double calculatePixelBrightness( unsigned int iChannel, double, double, double, double, double, double, double );

    public:

     VImageParameterFitter( VEvndispData* iData,
                            bool iDebug = false,
                            double ZeroTolerence = 1.e-8 );
    ~VImageParameterFitter();

     vector<bool> calcLL( VImageParameter *iParGeo,
                          VImageParameter *iParLL,
                          bool iUseSums2 = false, 
                          bool i_reInitializeLL = false, 
                          bool iEqualSummationWindows = false );          //!< calculate image parameters (log like)
     TF2* getNormal2D() { return fNormal2D; }
     void initMinuit( int );
     
     VDetectorGeometry* getDetectorGeometry()
     {
         if( fData ) return fData->getDetectorGeometry();
         return 0;
     }

     // public functions for image fitters
     vector<double>& getLLSums()                 //!< return data vector for minuit function
     {
        return fll_Sums;
     }
     vector<double>& getLLX()                    //!< return data vector for minuit function
     {
        return fll_X;
     }
     vector<double>& getLLY()                    //!< return data vector for minuit function
     {
        return fll_Y;
     }
     vector<double>& getLLR()
     {
        return fll_R;
     }
     vector<double>& getLLT()
     {
        return fll_T;
     }
     bool minimize_time_gradient_for_this_event();

};



#endif
