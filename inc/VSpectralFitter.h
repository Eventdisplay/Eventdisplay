//! VSpectralFitter fitter class for energy spectra (fit functions are predefined)

#ifndef VSpectralFitter_H
#define VSpectralFitter_H

#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"

#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

class VSpectralFitter : public TObject
{
    private:

        TF1*   fFitFunction;                             // fit functions (log energy axis)
        TF1*   fFitFunction_lin;                         // function for flux integration (lin energy axis)
        double* fFitFunction_CovarianceMatrix;         // covariance matrix from fit
        TFitResult* fFitResult;
        string fFitName;

        int    fSpectralFitFunction;
        double fSpectralFitFluxNormalisationEnergy;      // [TeV] linear axis

        double fSpectralFitEnergy_min;                   // [TeV] linear axis
        double fSpectralFitEnergy_max;                   // [TeV] linear axis

        // plotting variables
        int    fPlottingEnergySpectrumLineColor;
        int    fPlottingEnergySpectrumLineStyle;
        float  fPlottingEnergySpectrumLineWidth;

        bool   defineFitFunction();
        void   updateFitFunction_lin();

    public:

        VSpectralFitter( string fitname = "fit" );
        ~VSpectralFitter() {}

        TF1*   fit( TGraph* g, string fitname = "" );
        double getIntegralFlux( double iMinEnergy_TeV, double iMaxEnergy_TeV = 1.e6 );
        double getIntegralFluxError( double iMinEnergy_TeV, double iMaxEnergy_TeV = 1.e6 );
        TF1*   getSpectralFitFunction()
        {
            return fFitFunction;
        }
        double getSpectralFitNormalisationEnergy()
        {
            return fSpectralFitFluxNormalisationEnergy;
        }
        void   print();
        void   setSpectralFitFunction( int iD  = 0 )
        {
            fSpectralFitFunction = iD;
        }
        void   setSpectralFitFluxNormalisationEnergy( double iE_TeV = 1. )
        {
            fSpectralFitFluxNormalisationEnergy = iE_TeV;
        }
        void   setSpectralFitRangeLin( double xmin = 0.1, double xmax = 10. )
        {
            fSpectralFitEnergy_min = xmin;
            fSpectralFitEnergy_max = xmax;
        }
        void   setPlottingStyle( int iColor = 1, int iStyle = 1, float iWidth = 2. )
        {
            fPlottingEnergySpectrumLineColor = iColor;
            fPlottingEnergySpectrumLineStyle = iStyle;
            fPlottingEnergySpectrumLineWidth = iWidth;
        }

        ClassDef( VSpectralFitter, 1 );
};
#endif
