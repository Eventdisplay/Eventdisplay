/*!  \brief collection of useful analysis functions for significance and upper limit calculations
    \author
      Jamie Holder
      Gernot Maier

*/

#ifndef LIANDMA_C
#define LIANDMA_C

#include "TF1.h"
#include "TFeldmanCousins.h"
#include "TMath.h"
#include "TRolke.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TLegend.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

namespace VStatistics
{

    inline void liandma( double Non, double Noff, double alpha, double& Nsig, double& Sig5, double& Sig9, double& Sig17 )
    {
        if( alpha == 0. )
        {
            Sig17 = 0.;
            return;
        }
        
        double alphasq = alpha * alpha;
        double oneplusalpha = 1.0 + alpha;
        double oneplusalphaoveralpha = oneplusalpha / alpha;
        
        Nsig   = Non - alpha * Noff;
        double Ntot   = Non + Noff;
        
        if( Non + alphasq * Noff > 0. )
        {
            Sig5 = Nsig / sqrt( Non + alphasq * Noff );
        }
        else
        {
            Sig5 = 0.;
        }
        if( alpha * Ntot > 0. )
        {
            Sig9 = Nsig / sqrt( alpha * Ntot );
        }
        else
        {
            Sig9 = 0.;
        }
        if( Ntot == 0. )
        {
            Sig17 = 0.;
        }
        else if( Non == 0 && Noff != 0. )
        {
            Sig17 = sqrt( 2.*( Noff * log( oneplusalpha * ( Noff / Ntot ) ) ) );
        }
        else if( Non != 0 && Noff == 0. )
        {
            Sig17 = sqrt( 2.*( Non * log( oneplusalphaoveralpha * ( Non / Ntot ) ) ) );
        }
        else
        {
            Sig17 = 2.*( Non * log( oneplusalphaoveralpha * ( Non / Ntot ) ) + Noff * log( oneplusalpha * ( Noff / Ntot ) ) );
            // value in brackets can be a small negative number
            if( TMath::Abs( Sig17 ) < 1.e-15 )
            {
                Sig17 = 0.;
            }
            else
            {
                Sig17 = sqrt( Sig17 );
            }
        }
        if( Nsig < 0 )
        {
            Sig17 = -Sig17;
        }
        
        return;
    }
    
    inline double calcSignificance( double nOn, double nOff, double norm, int iLiMaForm = 17 )
    {
        double nSig = 0.;
        double limaSig5 = 0.;
        double limaSig9 = 0.;
        double limaSig17 = 0.;
        
        if( fabs( nOn ) < 1.e-5 && fabs( nOff ) < 1.e-5 )
        {
            return 0.;
        }
        
        liandma( nOn, nOff, norm, nSig, limaSig5, limaSig9, limaSig17 );
        //  if( !isnormal( limaSig17 ) ) return 0.;
        
        if( iLiMaForm == 5 )
        {
            return limaSig5;
        }
        if( iLiMaForm == 9 )
        {
            return limaSig9;
        }
        
        return limaSig17;
    }
    
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    //
    // upper limit calculations
    //
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
    /*
        Gaussian approximation, equ 10 and following paragraph of Helene (1983)
    */
    inline double Helene( double nOn, double nOff, double alpha, double CL )
    {
        double nDiff = nOn - alpha * nOff;
        double sigma = sqrt( nOn + alpha * alpha * nOff );
        
        double ulim = 0.;
        
        if( sigma > 0. )
        {
            ulim = TMath::ErfcInverse( ( 1. - CL ) * TMath::Erfc( -1. * nDiff / sigma / sqrt( 2. ) ) ) * sigma * sqrt( 2. ) + nDiff;
        }
        
        return ulim;
    }
    
    inline Double_t funcg( Double_t* x, Double_t* par )
    {
        double a = x[0];                          //!< upper limit
        double c = par[0];                        //!< on source counts
        double b = par[1] * par[2];               //!< off source counts * ratio
        
        return TMath::Poisson( c, a + b );
    }
    
    inline Double_t funcf( Double_t* x, Double_t* par )
    {
        double n1 = 1.;                           //!< Normalisation
        double p = 0.;                            //!< Probability
        double eps = 1E-5;                        //!< Accuracy
        double ul = x[0];                         //!< Upper limit
        double c = par[0];                        //<! Counts in the Peak Region
        double b = par[1] * par[2];               //<! Background Events * ratio
        double alpha = par[3];                    //<! Confidence Level
        
        double big = 10.*sqrt( b + c ) + fabs( c - b ) + 10.; //! 10 sigma noise
        TF1 g( "myfuncg", funcg, 0.0, big, 5 );
        g.SetParameters( par );
#ifndef ROOT6
        n1 = g.Integral( 0., big, par, eps );
#else
        n1 = g.Integral( 0., big, eps );
#endif
        if( n1 > 0. )
        {
#ifndef ROOT6
            p = g.Integral( ul, big, par, eps ) / n1;
#else
            p = g.Integral( ul, big, eps ) / n1;
#endif
        }
        else
        {
            p = 0.;
        }
        
        return alpha - p;
    }
    
    /*
        calculate upper limit using Helene, Feldman & Cousins, Rolke et al
    
        iMethod == 0: Helene
        iMethod == 1: Helene eq 2
        iMethod == 2: Helene eq 6
        iMethod == 3: Feldman & Cousins
        iMethod == 4: Rolke Model 3
        iMethod == 5: Rolke Model 4
    
        ratio = ratio between on and off exposure
    
        CL    = confidence limit (e.g. 0.99 for 99% probability)
    */
    inline double calcUpperLimit( double nOn, double nOff, double ratio, double CL, int iMethod = 0, int iBoundedLimits = true )
    {
        // Helene taking ratio into account
        if( iMethod == 0 || iMethod == 2 )
        {
            return Helene( nOn, nOff, ratio, CL );
        }
        // Helene, equ (2)
        else if( iMethod == 1 )
        {
            double mypars[5];
            mypars[0] = nOn;
            mypars[1] = nOff;
            mypars[2] = ratio;
            mypars[3] = 1. - CL;
            mypars[4] = ( double )iMethod;
            
            //! ICP and JB agree this is correct:
            if( nOn < nOff )
            {
                nOn = nOff;
            }
            
            TF1 f( "myfuncf", funcf, 0.0, 100000, 5 );
            f.SetParameters( mypars );
            f.SetNpx( 10000 );
            return f.GetX( 0.0 , 0.0, 100000, 1.e-10, 1000 );
        }
        // Feldman & Cousins
        else if( iMethod == 3 )
        {
            TFeldmanCousins i_FeldmanCousins( CL );
            if( nOn > 20 || ratio * nOff > 20 )
            {
                i_FeldmanCousins.SetMuMax( 100. );
            }
            if( nOn > 60 || ratio * nOff > 60 )
            {
                i_FeldmanCousins.SetMuMax( 1000. );
            }
            if( nOn > 1000 )
            {
                cout << "VStatistics::calcUpperLimit() warning: Feldman Cousins maximum values set to 1000" << endl;
                cout << "  --> results are not valid" << endl;
            }
            if( nOn > 40 )
            {
                i_FeldmanCousins.SetMuStep( 0.5 );
            }
            return i_FeldmanCousins.CalculateUpperLimit( nOn, ratio * nOff );
        }
        // Rolke Model 3 Background - Gaussian, Efficiency - Gaussian
        //
        // method4_em    efficiency
        // method4_sdem  standard deviation of efficiency
        else if( iMethod == 4 )
        {
            TRolke i_Rolke;
            i_Rolke.SetCL( CL );
            i_Rolke.SetBounding( iBoundedLimits );
            
            double sdb = ratio * sqrt( nOff );
            
            i_Rolke.SetGaussBkgGaussEff( ( int )nOn, ratio * nOff, 1., 0.3, sdb );
            
            return i_Rolke.GetUpperLimit();
        }
        // Rolke Model 4 Background - Poisson, Efficiency - known
        // method5_e    efficiency
        else if( iMethod == 5 && ratio > 0. )
        {
            TRolke i_Rolke;
            i_Rolke.SetCL( CL );
            i_Rolke.SetBounding( iBoundedLimits );
            i_Rolke.SetPoissonBkgKnownEff( ( int )nOn, ( int )nOff, 1. / ratio, 1. );
            
            return i_Rolke.GetUpperLimit();
        }
        else
        {
            cout << "unknown upper limit method: " << iMethod << endl;
            return 0;
        }
        return 0.;
    }
    
    /*!
        interpolate between two zenith angles (iCos = true )
        interpolate between two values (iCos = false )
    
        ze [deg]
    
        weighted by cos ze
    */
    inline double interpolate( double w1, double ze1, double w2, double ze2, double ze, bool iCos = false,
                               double iLimitforInterpolation = 0.5, double iMinValidValue = -90. )
    {
        // don't interpolate if both values are not valid
        if( w1 < iMinValidValue && w2 < iMinValidValue )
        {
            return -99.;
        }
        
        // same x-value, don't interpolate
        if( fabs( ze1 - ze2 ) < 1.e-3 )
        {
            if( w1 < iMinValidValue )
            {
                return w2;
            }
            else if( w2 < iMinValidValue )
            {
                return w1;
            }
            return ( w1 + w2 ) / 2.;
        }
        
        // interpolate
        double id = 0.;
        double f1 = 0.;
        double f2 = 0.;
        if( iCos )
        {
            id = cos( ze1 * TMath::DegToRad() ) - cos( ze2 * TMath::DegToRad() );
            f1 = 1. - ( cos( ze1 * TMath::DegToRad() ) - cos( ze * TMath::DegToRad() ) ) / id;
            f2 = 1. - ( cos( ze * TMath::DegToRad() ) - cos( ze2 * TMath::DegToRad() ) ) / id;
        }
        else
        {
            id = ze1 - ze2;
            f1 = 1. - ( ze1 - ze ) / id;
            f2 = 1. - ( ze  - ze2 ) / id;
        }
        
        // one of the values is not valid:
        // return valid value only when f1 or f2 > iLimitforInterPolation
        if( w1 > iMinValidValue && w2 < iMinValidValue )
        {
            if( f1 > iLimitforInterpolation )
            {
                return w1;
            }
            else
            {
                return -99.;
            }
        }
        else if( w1 < iMinValidValue && w2 > iMinValidValue )
        {
            if( f2 > iLimitforInterpolation )
            {
                return w2;
            }
            else
            {
                return -99.;
            }
        }
        
        return ( w1 * f1 + w2 * f2 );
    }
    
    /*
    
       mean
    
    */
    inline double getMean( vector< double >& x )
    {
        double sum = 0.;
        for( unsigned int i = 0; i < x.size(); i++ )
        {
            sum += x[i];
        }
        if( x.size() > 0. )
        {
            return sum / ( ( double )x.size() );
        }
        
        return -999.;
    }
    
    /*
    
      weighted mean
    
    */
    inline double getWeightedMean( vector< double >& x, vector< double >& w )
    {
        if( x.size() != w.size() )
        {
            return -999.;
        }
        
        double sum = 0.;
        double weight = 0.;
        
        for( unsigned int i = 0; i < x.size(); i++ )
        {
            sum += x[i] * w[i];
            weight += w[i];
        }
        if( weight > 0. )
        {
            return sum / ( weight );
        }
        
        return -999.;
    }
    
    /*
    
      median
    
    */
    inline double getMedian( vector< double > x )
    {
        return TMath::Median( x.size(), &x[0] );
    }
    
    /*
    
       median absolute deviation
    
       (note for a normal distribution holds sigma = 1.4826 * medianAbsoluteDeviation)
    
    */
    inline double getMedianAbsoluteDeviation( vector< double > x )
    {
        double iMedian = getMedian( x );
        vector< double > MAD( x.size(), 0. );
        for( unsigned int i = 0; i < MAD.size(); i++ )
        {
            MAD[i] = TMath::Abs( x[i] - iMedian );
        }
        
        return getMedian( MAD );
    }
    
    /*
    
       mean absolute error
    
    */
    inline double getMeanAbsoluteError( vector< double >& x )
    {
        double mean = getMean( x );
        double iAbs = 0.;
        for( unsigned int i = 0; i < x.size(); i++ )
        {
            iAbs += TMath::Abs( x[i] - mean );
        }
        if( x.size() > 0. )
        {
            return iAbs / ( ( double )x.size() );
        }
        
        return 0.;
    }
    /*
     *
     * median absolute error
     *
     */
    inline double getMedianAbsoluteError( vector< double >& x, double median )
    {
        double iAbs = 0.;
        for( unsigned int i = 0; i < x.size(); i++ )
        {
            iAbs += TMath::Abs( x[i] - median );
        }
        if( x.size() > 0. )
        {
            return iAbs / ( ( double )x.size() );
        }
        
        return 0.;
    }
    /*
    
       trimmed mean absolute error
    
           trimMLow = trim N elements on lower side
           trimMLow = trim N elements on upper side
    
           (Note: this is nowhere used that this point and therefore untested)
    
    */
    inline double getMeanAbsoluteError( vector< double > x, unsigned int trimMLow, unsigned int trimMUp )
    {
        // sort vector
        std::sort( x.begin(), x.end() );
        
        // trimmed mean
        if( x.size() > trimMLow + trimMUp + 1 )
        {
            vector< double > y( x.begin() + trimMLow, x.begin() + x.size() - trimMUp );
            double mean = getMean( y );
            double iAbs = 0.;
            for( unsigned int i = 0; i < y.size(); i++ )
            {
                iAbs += TMath::Abs( y[i] - mean );
            }
            if( y.size() > 0 )
            {
                return iAbs / ( ( double )y.size() );
            }
        }
        
        return 0.;
    }
    
}
#endif
