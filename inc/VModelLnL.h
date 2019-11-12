//! VModelLnL (adapted from Matthew Wood's SLAC analysis)
//! Poisson likelihood function for a model analysis.

#ifndef VMODELLNL_H
#define VMODELLNL_H

#include "TFile.h" // for lookup table
#include "TH2D.h"  // for lookup table
#include "TH1D.h"  // for lookup table

#include "VMinimizer.h"
#include "VModel3DData.h"
#include "VModel3DFn.h"

using namespace std;

class VModelLnL : public VMinimizerFn
{
    private:
        VModel3DData* fData3D;    // pointer to data3D class
        VModel3DFn* fModel3DFn;   // pointer to Model3D function class
        
        void scopeLnL( unsigned int iscope, const vector<double>& a, double& lnl, vector<double>& beta, vector< vector<double> >& alpha, double& gof ) const; //! sum of pixel log-likelihoods for the given telescope
        
        void pixelLnLGauss( double s, double pixel_var, double pixel_nsb, double mu, double& lnl, double& dlnl, double& d2lnl ) const; //! Compute the pixel likelihood using a Gaussian distribution for a given signal s and model amplitude mu
        
        static void pixelLnL( double s, double pixel_var, double pixel_nsb, double singlepe_var, double mu, double& lnl, double& dlnl, double& d2lnl ); //! Compute the pixel likelihood for a given signal s and model amplitude mu
        
        void pixelLnL_lookup( double s, double pixel_var, double pixel_nsb, double mu, double& lnl, double& dlnl, double& d2lnl ) const;
        
        bool fDebug;             // print debug info
        unsigned int  fNParam;   // number of parameters in Model3D
        double fSinglePEvar;     // probably not needed for VERITAS
        
        /// lookup table ///
        TFile* fLnLFile;
        TH2D* m_lnl_lmu;
        TH2D* m_dlnl_lmu;
        TH2D* m_d2lnl_lmu;
        TH2D* m_lnl_smu;
        TH2D* m_dlnl_smu;
        TH2D* m_d2lnl_smu;
        TH1D* m_lnl_exp;
        
    public:
    
        VModelLnL( VModel3DData* iVModel3DData, VModel3DFn* iVModel3DFn );
        ~VModelLnL();
        
        void readLnLTable( string LnLTableFile );
        void createLnLTable( double pixel_var, double pixel_nsb );
        
        void val( const vector<double>& a, double& lnl, vector<double>& beta, vector< vector< double > >& alpha ) const;
        
};

#endif // VMODELLNL_H
