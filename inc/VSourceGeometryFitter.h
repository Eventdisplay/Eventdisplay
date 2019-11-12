//! VSourceGeometryFitter analyse source geometry (position and extension)

#ifndef VSourceGeometryFitter_H
#define VSourceGeometryFitter_H

#include "VAnalysisUtilities.h"
#include "VPlotUtilities.h"
#include "VPlotAnasumHistograms.h"

#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "Minuit2/FCNBase.h"
#include "TMinuit.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <iostream>
#include <string>

using namespace std;


///////////////////////////////////////////////////////////////////////////////
class VSourceGeometryFitterData
{
    public:
    
        string           fFitterName;
        string           fFitterDescription;
        vector< string > fParameterName;
        vector< double > fParameterInitValue;
        vector< double > fParameterStep;
        vector< double > fParameterLowerLimit;
        vector< double > fParameterUpperLimit;
        int              fFitResult_Status;
        vector< double > fFitResult_Parameter;
        vector< double > fFitResult_ParameterError;
        
        VSourceGeometryFitterData();
        ~VSourceGeometryFitterData() {}
};

class VSourceGeometryFitter : public VAnalysisUtilities, public VPlotUtilities
{
    private:
    
        bool   fDebug;
        
        string fAnasumDataFile;
        int    fRunNumber;
        double fXStart;
        double fYStart;
        double fPSF;
        
        // sky map to be fitted
        TH2D* fHisSkyMap;
        
        // default fitter data
        vector< VSourceGeometryFitterData* > fDefaultFitterData;
        
        // fitter used
        VSourceGeometryFitterData*            fFitter;
        
        void setFitterDefaultData();
        
    public:
    
        VSourceGeometryFitter();
        VSourceGeometryFitter( string iAnaSumDataFile, int irun = -1 );
        ~VSourceGeometryFitter() {}
        
        void     fitSource( string iHisName = "hmap_stereoUC_diff", double xStart = 0., double yStart = 0., double xyRange = 0.15 );
        TH2D*    getSkyMap()
        {
            return fHisSkyMap;
        }
        void     help();
        TCanvas* plot( double rmax = 0.2, double zmin = -1000., double zmax = -1000., string iPlotMode = "colz" );
        void     plotFitResult();
        TGraph*  plotSourceGeometry( int iColor = 1 );
        void     setDebug( bool iB = true )
        {
            fDebug = iB;
        }
        bool     setFitter( string iFitter );
        void     setPSF( double psf )
        {
            fPSF = psf;
        }
        double   getPSF()
        {
            return fPSF;
        }
        
        ClassDef( VSourceGeometryFitter, 1 );
};




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Functions for the Point Spread Function
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  PSF description (1) - radial symmetric PSF with an offset
//
///////////////////////////////////////////////////////////////////////////////
class VFun_PSFDescription_2DGauss_Chi2 : public ROOT::Math::IBaseFunctionMultiDim
{
    private:
        TH2D* hSkyMap;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        
        /////////////////////////////////
        // function to be minimized
        double DoEval( const double* par ) const
        {
            // set variables
            double x = 0.;
            double y = 0.;
            
            double sum = 0.;
            double fT = 0.;
            double fH = 0.;
            double fHErr = 0.;
            
            double t2 = 0.;
            double sigmaSource2 = par[2] * par[2];
            
            // loop over sky map
            if( hSkyMap )
            {
                int nbinsX = hSkyMap->GetNbinsX();
                int nbinsY = hSkyMap->GetNbinsY();
                for( int i = 1; i <= nbinsX; i++ )
                {
                    x = hSkyMap->GetXaxis()->GetBinCenter( i );
                    // check x-range
                    if( x > xmax )
                    {
                        continue;
                    }
                    if( x < xmin )
                    {
                        continue;
                    }
                    for( int j = 1; j <= nbinsY; j++ )
                    {
                        y = hSkyMap->GetYaxis()->GetBinCenter( j );
                        // check y-range
                        if( y > ymax )
                        {
                            continue;
                        }
                        if( y < ymin )
                        {
                            continue;
                        }
                        
                        // // skip empty bins
                        // if( hSkyMap->GetBinContent( i, j ) <= 0. )
                        // {
                        // 	continue;
                        // }
                        
                        // calculate theta2
                        t2 = ( x - par[3] ) * ( x - par[3] ) + ( y - par[4] ) * ( y - par[4] );
                        
                        // calculate expectation from model function
                        fT = par[0] + par[1] * TMath::Exp( -1.*t2 / 2. / sigmaSource2 );
                        if( TMath::IsNaN( fT ) )
                        {
                            continue;
                        }
                        
                        // get value and error in histogram
                        fH = hSkyMap->GetBinContent( i, j );
                        fHErr = hSkyMap->GetBinError( i, j );
                        
                        // calculate chi2
                        if( fHErr > 0. && fH > -90. )
                        {
                            sum += ( fT - fH ) * ( fT - fH ) / fHErr / fHErr;
                        }
                    }
                }
            }
            return sum;
        }
        
    public:
        VFun_PSFDescription_2DGauss_Chi2( TH2D* iSkymap = 0, double i_xmin = -1., double i_xmax = 1.,
                                          double i_ymin = -1., double i_ymax = 1. );
                                          
        unsigned int NDim() const
        {
            return 5;
        }
        
        ROOT::Math::IBaseFunctionMultiDim* Clone() const
        {
            ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_PSFDescription_2DGauss_Chi2( hSkyMap, xmin, xmax, ymin, ymax );
            return clone;
        }
};



///////////////////////////////////////////////////////////////////////////////
// PSF description (2) - radial symmetric gaussian
//
//
///////////////////////////////////////////////////////////////////////////////
class VFun_PSFDescription_2DGauss_LL : public ROOT::Math::IBaseFunctionMultiDim
{
    private:
    
        TH2D* hSkyMap;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        
        /////////////////////////////////
        // function to be minimized
        double DoEval( const double* par ) const
        {
            // initialize variables
            double  LL = 0.;
            double  sum = 0.;
            double  meanX  = par[0];
            double  meanY  = par[1];
            double  sigma = par[2];
            
            double x = 0.;
            double y = 0.;
            double n = 0.;                                // measured sum in channel i
            
            if( sigma >= 0. ) // && sigmaH >= 0. )
            {
                int nbinsX = hSkyMap->GetNbinsX();
                int nbinsY = hSkyMap->GetNbinsY();
                for( int i = 1; i <= nbinsX; i++ )
                {
                    x = hSkyMap->GetXaxis()->GetBinCenter( i );
                    // check x-range
                    if( x > xmax )
                    {
                        continue;
                    }
                    if( x < xmin )
                    {
                        continue;
                    }
                    for( int j = 1; j <= nbinsY; j++ )
                    {
                        y = hSkyMap->GetYaxis()->GetBinCenter( j );
                        // check y-range
                        if( y > ymax )
                        {
                            continue;
                        }
                        if( y < ymin )
                        {
                            continue;
                        }
                        
                        n = hSkyMap->GetBinContent( i, j );
                        
                        // check for valid entries
                        if( n > -999. )
                        {
                            // calculate log-likelihood
                            sum  = ( x - meanX ) * ( x - meanX ) / sigma / sigma ;
                            sum += ( y - meanY ) * ( y - meanY ) / sigma / sigma ;
                            sum *= -1. / 2. ;
                            sum  = 1. / 2. / M_PI / sigma / sigma * exp( sum );
                            
                            // assume Poisson fluctuations (neglecting background noise)
                            if( n > 0. && sum > 0. )
                            {
                                LL += n * log( sum ) - sum - n * log( n ) + n;
                            }
                            else
                            {
                                LL += -1. * sum;
                            }
                        }
                    }
                }
            }
            
            return -1. * LL;
        }
        
    public:
    
        VFun_PSFDescription_2DGauss_LL( TH2D* iSkymap = 0, double i_xmin = -1., double i_xmax = 1., double i_ymin = -1., double i_ymax = 1. );
        
        unsigned int NDim() const
        {
            return 3;
        }
        
        ROOT::Math::IBaseFunctionMultiDim* Clone() const
        {
            ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_PSFDescription_2DGauss_LL( hSkyMap, xmin, xmax, ymin, ymax );
            return clone;
        }
};



///////////////////////////////////////////////////////////////////////////////
// PSF description (3) - central spot + broad halo
//  radial symmetry is assumed, superposition of two radial symmetric gaussians
//
//
///////////////////////////////////////////////////////////////////////////////
class VFun_PSFDescription_LinearSuperposition2DGauss_LL: public ROOT::Math::IBaseFunctionMultiDim
{
    private:
    
        TH2D* hSkyMap;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        
        /////////////////////////////////
        // function to be minimized
        double DoEval( const double* par ) const
        {
            // initialize variables
            double  LL = 0.;
            double  sum = 0.;
            double  sum1 = 0.; // central spot
            double  sum2 = 0.; // broad halo
            double  meanX  = par[0];
            double  meanY  = par[1];
            double  sigma1 = par[2];
            double  sigma2 = par[3];
            double  alpha  = par[4]; // relative importance of each component , alpha = 1 only central spot matters
            
            double x = 0.;
            double y = 0.;
            double n = 0.;                                // measured sum in channel i
            
            if( sigma1 >= 0. && sigma2 >= 0. )
            {
                int nbinsX = hSkyMap->GetNbinsX();
                int nbinsY = hSkyMap->GetNbinsY();
                for( int i = 1; i <= nbinsX; i++ )
                {
                    x = hSkyMap->GetXaxis()->GetBinCenter( i );
                    // check x-range
                    if( x > xmax )
                    {
                        continue;
                    }
                    if( x < xmin )
                    {
                        continue;
                    }
                    for( int j = 1; j <= nbinsY; j++ )
                    {
                        y = hSkyMap->GetYaxis()->GetBinCenter( j );
                        // check y-range
                        if( y > ymax )
                        {
                            continue;
                        }
                        if( y < ymin )
                        {
                            continue;
                        }
                        
                        n = hSkyMap->GetBinContent( i, j );
                        
                        // check for valid entries
                        if( n > -999. )
                        {
                            // calculate log-likelihood
                            sum1  = ( x - meanX ) * ( x - meanX ) / sigma1 / sigma1 ;
                            sum1 += ( y - meanY ) * ( y - meanY ) / sigma1 / sigma1 ;
                            sum1 *= -1. / 2. ;
                            
                            sum2  = ( x - meanX ) * ( x - meanX ) / sigma2 / sigma2 ;
                            sum2 += ( y - meanY ) * ( y - meanY ) / sigma2 / sigma2 ;
                            sum2 *= -1. / 2. ;
                            
                            sum  = alpha * sqrt( 1. / 2. / M_PI / sigma1 / sigma1 ) * exp( sum1 );
                            sum += ( 1 - alpha ) * sqrt( 1. / 2. / M_PI / sigma2 / sigma2 ) * exp( sum2 );
                            
                            
                            // assume Poisson fluctuations (neglecting background noise)
                            if( n > 0. && sum > 0. )
                            {
                                LL += n * log( sum ) - sum - n * log( n ) + n;
                            }
                            else
                            {
                                LL += -1. * sum;
                            }
                        }
                    }
                }
            }
            return -1. * LL;
        }
        
    public:
    
        VFun_PSFDescription_LinearSuperposition2DGauss_LL( TH2D* iSkymap = 0, double i_xmin = -1., double i_xmax = 1., double i_ymin = -1., double i_ymax = 1. );
        
        unsigned int NDim() const
        {
            return 5;
        }
        
        ROOT::Math::IBaseFunctionMultiDim* Clone() const
        {
            ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_PSFDescription_LinearSuperposition2DGauss_LL( hSkyMap, xmin, xmax, ymin, ymax );
            return clone;
        }
};



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Functions for source position and extension fitting
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Source Description (1); radial symmetric source, Chi2
///////////////////////////////////////////////////////////////////////////////
class VFun_SourceDescription_RadialSymmetricSource_Chi2 : public ROOT::Math::IBaseFunctionMultiDim
{
    private:
    
        TH2D* hSkyMap;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double sigmaPSF;
        
        /////////////////////////////////
        // function to be minimized
        double DoEval( const double* par ) const
        {
            // set variables
            double x = 0.;
            double y = 0.;
            
            double sum = 0.;
            double fT = 0.;
            double fH = 0.;
            double fHErr = 0;
            
            double t2 = 0.;
            double sigmaSRC = par[2];
            
            // loop over sky map
            if( hSkyMap )
            {
                int nbinsX = hSkyMap->GetNbinsX();
                int nbinsY = hSkyMap->GetNbinsY();
                for( int i = 1; i <= nbinsX; i++ )
                {
                    x = hSkyMap->GetXaxis()->GetBinCenter( i );
                    // check x-range
                    if( x > xmax )
                    {
                        continue;
                    }
                    if( x < xmin )
                    {
                        continue;
                    }
                    for( int j = 1; j <= nbinsY; j++ )
                    {
                        y = hSkyMap->GetYaxis()->GetBinCenter( j );
                        // check y-range
                        if( y > ymax )
                        {
                            continue;
                        }
                        if( y < ymin )
                        {
                            continue;
                        }
                        
                        // // skip empty bins
                        // if( hSkyMap->GetBinContent( i, j ) <= 0. )
                        // {
                        // 	continue;
                        // }
                        
                        // calculate theta2
                        t2 = ( x - par[0] ) * ( x - par[0] ) + ( y - par[1] ) * ( y - par[1] );
                        
                        // calculate expectation from model function
                        fT = par[3] * TMath::Exp( -1.*t2 / 2. / ( sigmaSRC * sigmaSRC + sigmaPSF * sigmaPSF ) );
                        if( TMath::IsNaN( fT ) )
                        {
                            continue;
                        }
                        
                        // get value and error in histogram
                        fH = hSkyMap->GetBinContent( i, j );
                        fHErr = hSkyMap->GetBinError( i, j );
                        
                        // calculate chi2
                        if( fHErr > 0. && fH > -90. )
                        {
                            sum += ( fT - fH ) * ( fT - fH ) / fHErr / fHErr;
                        }
                    }
                }
            }
            return sum;
        }
        
    public:
    
        VFun_SourceDescription_RadialSymmetricSource_Chi2( TH2D* iSkymap = 0, double i_xmin = -1., double i_xmax = 1., double i_ymin = -1., double i_ymax = 1., double i_psf = 0.063 );
        
        unsigned int NDim() const
        {
            return 4;
        }
        
        ROOT::Math::IBaseFunctionMultiDim* Clone() const
        {
            ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_SourceDescription_RadialSymmetricSource_Chi2( hSkyMap, xmin, xmax, ymin, ymax, sigmaPSF );
            return clone;
        }
};



///////////////////////////////////////////////////////////////////////////////
// Source Description (2); radial asymmetric source, Chi2  NOT WORKING YET
///////////////////////////////////////////////////////////////////////////////
/*class VFun_SourceDescription_RadialAsymmetricSource_Chi2 : public ROOT::Math::IBaseFunctionMultiDim
{
   private:

   TH2D *hSkyMap;
   double xmin;
   double xmax;
   double ymin;
   double ymax;

/////////////////////////////////
// function to be minimized
   double DoEval( const double* par ) const
   {
       // set variables
       double sigmaPSF=0.063;

       double x = 0.;
       double y = 0.;

       double sum = 0.;
       double fT = 0.;
       double fH = 0.;
       double meanX      = par[0];
       double meanY      = par[1];
       //  double sigmaSRC_X = par[2];
       // double sigmaSRC_Y = par[3];
       double theta = par[4]*acos(-1)/180.;
       //double sigmaX = sqrt( sigmaPSF * sigmaPSF + sigmaSRC_X * sigmaSRC_X );
       //double sigmaY = sqrt( sigmaPSF * sigmaPSF + sigmaSRC_Y * sigmaSRC_Y );
       double sigmaX = par[2];
       double sigmaY = par[3];


       double a = cos(theta) * cos(theta) / 2. / sigmaX / sigmaX  +  sin(theta) * sin(theta) / 2. / sigmaY / sigmaY;
       double b = (-1) * sin(2*theta) / 4. / sigmaX / sigmaX + sin(2*theta) / 4. / sigmaY / sigmaY;
       double c = sin(theta) * sin(theta) / 2./ sigmaX / sigmaX + cos(theta) * cos(theta)/ 2./ sigmaY / sigmaY;



// loop over sky map
       if( hSkyMap )
       {
	  int nbinsX = hSkyMap->GetNbinsX();
	  int nbinsY = hSkyMap->GetNbinsY();
	  for( int i = 1; i <= nbinsX; i++ )
	  {
	      x = hSkyMap->GetXaxis()->GetBinCenter( i );
// check x-range
	      if( x > xmax ) continue;
      	      if( x < xmin ) continue;
	      for( int j = 1; j <= nbinsY; j++ )
	      {
		  y = hSkyMap->GetYaxis()->GetBinCenter( j );
// check y-range
		  if( y > ymax ) continue;
		  if( y < ymin ) continue;

// skip empty bins
		  if( hSkyMap->GetBinContent( i, j ) <= 0. ) continue;


// calculate expectation from model function
		  fT = a * (x - meanX) * (x - meanX) + 2 * b *(x - meanX)*(y - meanY) +  c * (y - meanY) * (y - meanY);
		  fT = par[5]*TMath::Exp(-1 * fT);

		  if( TMath::IsNaN( fT ) ) continue;

// get value and error in histogram
		  fH = hSkyMap->GetBinContent( i, j );

// calculate chi2
		  if( fH != 0. && fH > -90. ) sum += ( fT - fH ) * ( fT - fH ) / fH;
	       }
	   }
       }
       return sum;
   }

   public:

   VFun_SourceDescription_RadialAsymmetricSource_Chi2( TH2D *iSkymap = 0, double i_xmin = -1., double i_xmax = 1., double i_ymin = -1., double i_ymax = 1. );

		unsigned int NDim() const
		{
			return 6;
		}

		ROOT::Math::IBaseFunctionMultiDim * Clone() const
		{
			ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_SourceDescription_RadialAsymmetricSource_Chi2( hSkyMap, xmin, xmax, ymin, ymax );
			return clone;
		}
};
*/



///////////////////////////////////////////////////////////////////////////////
// Source Description (3): Radial Symmetric Sources, LL
// TODO: needs more work to take zero and negative bins into account
///////////////////////////////////////////////////////////////////////////////
class VFun_SourceDescription_RadialSymmetricSource_LL: public ROOT::Math::IBaseFunctionMultiDim
{
    private:
    
        TH2D* hSkyMap;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double sigmaPSF;
        
        /////////////////////////////////
        // function to be minimized
        double DoEval( const double* par ) const
        {
            // initialize variables
            double  LL = 0.;
            double  sum = 0.;
            double  meanX = par[0];
            double  meanY = par[1];
            double  sigmaSRC = par[2];
            
            double x = 0.;
            double y = 0.;
            double n = 0.;                                // measured sum in channel i
            
            if( sigmaSRC > 0. )
            {
                int nbinsX = hSkyMap->GetNbinsX();
                int nbinsY = hSkyMap->GetNbinsY();
                for( int i = 1; i <= nbinsX; i++ )
                {
                    x = hSkyMap->GetXaxis()->GetBinCenter( i );
                    // check x-range
                    if( x > xmax )
                    {
                        continue;
                    }
                    if( x < xmin )
                    {
                        continue;
                    }
                    for( int j = 1; j <= nbinsY; j++ )
                    {
                        y = hSkyMap->GetYaxis()->GetBinCenter( j );
                        // check y-range
                        if( y > ymax )
                        {
                            continue;
                        }
                        if( y < ymin )
                        {
                            continue;
                        }
                        
                        n = hSkyMap->GetBinContent( i, j );
                        
                        // check for valid entries
                        if( n > -999. )
                        {
                            // calculate log-likelihood
                            sum  = ( x - meanX ) * ( x - meanX ) / ( sigmaSRC * sigmaSRC + sigmaPSF * sigmaPSF );
                            sum += ( y - meanY ) * ( y - meanY ) / ( sigmaSRC * sigmaSRC + sigmaPSF * sigmaPSF );
                            sum *= -1. / 2.;
                            sum  = exp( sum );
                            sum *= 1. / 2. / M_PI / ( sigmaSRC * sigmaSRC + sigmaPSF * sigmaPSF );
                            
                            // assume Poisson fluctuations (neglecting background noise)
                            if( n > 0. && sum > 0. )
                            {
                                LL += n * log( sum ) - sum - n * log( n ) + n;
                            }
                            else
                            {
                                LL += -1. * sum;
                            }
                        }
                    }
                }
            }
            return -1. * LL;
        }
        
    public:
    
        VFun_SourceDescription_RadialSymmetricSource_LL( TH2D* iSkymap = 0, double i_xmin = -1., double i_xmax = 1., double i_ymin = -1., double i_ymax = 1., double i_psf = 0.063 );
        
        unsigned int NDim() const
        {
            return 3;
        }
        
        ROOT::Math::IBaseFunctionMultiDim* Clone() const
        {
            ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_SourceDescription_RadialSymmetricSource_LL( hSkyMap, xmin, xmax, ymin, ymax, sigmaPSF );
            return clone;
        }
};



///////////////////////////////////////////////////////////////////////////////
// Source Description (4): Radial asymmetric gaussian, convolved with simple PSF
// TODO: needs more work to take zero and negative bins into account
///////////////////////////////////////////////////////////////////////////////
class VFun_SourceDescription_RadialAsymmetricSource_LL: public ROOT::Math::IBaseFunctionMultiDim
{
    private:
    
        TH2D* hSkyMap;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double sigmaPSF;
        
        /////////////////////////////////
        // function to be minimized
        double DoEval( const double* par ) const
        {
            // initialize variables
            double  LL = 0.;
            double  sum = 0.;
            //    double  rho =  par[0];
            double  meanX = par[1];
            double  sigmaX = par[2];
            double  meanY = par[3];
            double  sigmaY = par[4];
            double  angle = par[0];
            double  rho = 1. / 2. * tan( 2 * angle ) * ( sigmaX * sigmaX - sigmaY * sigmaY ) / sqrt( sigmaX * sigmaX + sigmaPSF * sigmaPSF ) / sqrt( sigmaY * sigmaY + sigmaPSF * sigmaPSF );
            
            
            double x = 0.;
            double y = 0.;
            double n = 0.;                                // measured sum in channel i
            
            if( rho * rho < 1. && sigmaX > 0. && sigmaY > 0. )
            {
                int nbinsX = hSkyMap->GetNbinsX();
                int nbinsY = hSkyMap->GetNbinsY();
                for( int i = 1; i <= nbinsX; i++ )
                {
                    x = hSkyMap->GetXaxis()->GetBinCenter( i );
                    // check x-range
                    if( x > xmax )
                    {
                        continue;
                    }
                    if( x < xmin )
                    {
                        continue;
                    }
                    for( int j = 1; j <= nbinsY; j++ )
                    {
                        y = hSkyMap->GetYaxis()->GetBinCenter( j );
                        // check y-range
                        if( y > ymax )
                        {
                            continue;
                        }
                        if( y < ymin )
                        {
                            continue;
                        }
                        
                        n = hSkyMap->GetBinContent( i, j );
                        
                        // check for valid entries
                        if( n > -999. )
                        {
                            // calculate log-likelihood
                            sum  = ( x - meanX ) * ( x - meanX ) / ( sigmaX * sigmaX + sigmaPSF * sigmaPSF );
                            sum += ( y - meanY ) * ( y - meanY ) / ( sigmaY * sigmaY + sigmaPSF * sigmaPSF );
                            sum += -2. * rho * ( x - meanX ) / sqrt( sigmaX * sigmaX + sigmaPSF * sigmaPSF ) * ( y - meanY ) / sqrt( sigmaY * sigmaY + sigmaPSF * sigmaPSF );
                            sum *= -1. / 2. / ( 1. - rho * rho );
                            sum  = exp( sum );
                            sum *= 1. / 2. / M_PI / sqrt( sigmaX * sigmaX + sigmaPSF * sigmaPSF ) /  sqrt( sigmaY * sigmaY + sigmaPSF * sigmaPSF ) / sqrt( 1. - rho * rho );
                            
                            // assume Poisson fluctuations (neglecting background noise)
                            if( n > 0. && sum > 0. )
                            {
                                LL += n * log( sum ) - sum - n * log( n ) + n;
                            }
                            else
                            {
                                LL += -1. * sum;
                            }
                        }
                    }
                }
            }
            return -1. * LL;
        }
        
    public:
    
        VFun_SourceDescription_RadialAsymmetricSource_LL( TH2D* iSkymap = 0, double i_xmin = -1., double i_xmax = 1., double i_ymin = -1., double i_ymax = 1., double i_psf = 0.063 );
        
        unsigned int NDim() const
        {
            return 5;
        }
        
        ROOT::Math::IBaseFunctionMultiDim* Clone() const
        {
            ROOT::Math::IBaseFunctionMultiDim* clone = new VFun_SourceDescription_RadialAsymmetricSource_LL( hSkyMap, xmin, xmax, ymin, ymax, sigmaPSF );
            return clone;
        }
};


#endif  // end of header

