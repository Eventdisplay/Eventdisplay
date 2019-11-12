//! VMinimizer (adapted from Matthew Wood's SLAC analysis)

#ifndef VMINIMIZER_H
#define VMINIMIZER_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <memory>
#include <limits> //(numeric_limits)

using namespace std;

class VMinimizerFn
{
    private:
        unsigned int fNParam;
        
    public:
        VMinimizerFn( unsigned int NParam ): fNParam( NParam ) { }
        virtual ~VMinimizerFn() { }
        unsigned int getNParam() const
        {
            return fNParam;
        }
        virtual void val( const vector<double>& a, double& lnl, vector<double>& beta, vector< vector< double > >& alpha ) const { }
};

/////////////////////////////////////////////////////////////////

class VMinimizer
{
    private:
        VMinimizer( const VMinimizer& o );
        VMinimizer& operator=( const VMinimizer& o );
        const VMinimizerFn&  fFn;
        
        vector<bool>   fHasLoBound;
        vector<double> fLoBound;
        vector<bool>   fHasHiBound;
        vector<double> fHiBound;
        vector<double> fParam;
        vector<double> fError;
        vector< vector< double > > fCov;
        double         fFnVal;
        bool           fConverged;
        
    protected:
        unsigned       fNFit;   //!< Number of free parameters
        vector<bool>   fFit;    //!< Vector designating free parameters
        
    public:
        VMinimizer( const VMinimizerFn& fn, const vector<double>& param, const vector<double>& err );
        virtual ~VMinimizer() { }
        //// Get Methods
        virtual unsigned getNParam() const
        {
            return fParam.size();    //! number of parameters.
        }
        virtual unsigned getNFit() const
        {
            return fNFit;    //! number of fitted parameters.
        }
        virtual const vector<double>& getParam() const
        {
            return fParam;    //! vector of parameters.
        }
        virtual double getParam( unsigned ip ) const
        {
            return fParam[ip];    //! parameter value.
        }
        virtual const vector<double>& getError() const
        {
            return fError;    //! vector of parameter errors.
        }
        virtual double getError( unsigned ip ) const
        {
            return fError[ip];    //! parameter error.
        }
        bool getIsFree( unsigned ip ) const
        {
            return fFit[ip];    //! free parameters.
        }
        virtual const vector< vector< double > >& getCov() const
        {
            return fCov;    //! vectors of parameter covariance
        }
        virtual double getCov( unsigned ip1, unsigned ip2 ) const
        {
            return fCov[ip1][ip2];    //! covariance of two parameters.
        }
        virtual double getRho( unsigned ip1, unsigned ip2 ) const
        {
            return getCov( ip1, ip2 ) / ( getError( ip1 ) * getError( ip2 ) );    //! correlation coefficient of two parameters.
        }
        virtual double getFnVal() const
        {
            return fFnVal;    //! value of object function.
        }
        const VMinimizerFn& getFn() const
        {
            return fFn;    //! object of function.
        }
        bool getConverged() const
        {
            return fConverged;    //! fit converged or not
        }
        double getHiBound( unsigned ip ) const
        {
            return fHiBound[ip];    //! high bound
        }
        bool getHasHiBound( unsigned ip ) const
        {
            return fHasHiBound[ip];    //! has high bound
        }
        double getLoBound( unsigned ip ) const
        {
            return fLoBound[ip];    //! low bound
        }
        bool getHasLoBound( unsigned ip ) const
        {
            return fHasLoBound[ip];    //! has low bound
        }
        //// Set Methods
        void setParam( const vector<double>& param )
        {
            fParam = param;
        }
        void setParam( unsigned ip, double val )
        {
            fParam[ip] = val;
        }
        void setError( unsigned ip, double err )
        {
            fError[ip] = err;
        }
        void setCovariance( const vector< vector< double > >& cov )
        {
            fCov = cov;
            for( unsigned int iparm = 0; iparm < getNParam(); iparm++ )
            {
                setError( iparm, sqrt( fCov[iparm][iparm] ) );
            }
        }
        virtual void setLoBound( unsigned ip, double val )
        {
            fHasLoBound[ip] = true;    //! Set a lower bound on a parameter.
            fLoBound[ip] = val;
        }
        virtual void setHiBound( unsigned ip, double val )
        {
            fHasHiBound[ip] = true;    //! Set an upper bound on a parameter.
            fHiBound[ip] = val;
        }
        void setFnVal( double val )
        {
            fFnVal = val;
        }
        void setConverged( bool converged )
        {
            fConverged = converged;
        }
        //// other functions
        void freeze( unsigned ip )
        {
            if( ip < fFit.size() )
            {
                fFit[ip] = false;
            }
            fNFit = 0;
            for( unsigned ifit = 0; ifit < fFit.size(); ifit++ )
            {
                if( fFit[ifit] )
                {
                    fNFit++;
                }
            }
        }
        virtual void minimize() = 0;
};

/////////////////////////////////////////////////////////////////

class VLevenbergMarquardt : public VMinimizer
{
    private:
        void initialize();
        void iterate();
        bool mrqcof( const vector<double>& param, vector< vector< double > >& alpha, vector<double>& beta, double& fnval );
        vector<double>              fDa;     //!< Step in each of the parameters.
        vector< vector< double > >  fAlpha;  //!< Approximation to Hessian matrix
        vector<double>              fBeta;   //!< Vector of first derivatives
        vector< vector< double > >  fCov;    //!< Covariance matrix
        unsigned int                fNiter;
        double                      fLambda; //!< Scaling parameter for diagonal of Hessian
        double                      fFnVal;
        double                      fTol;    //!< Tolerance for function minimization
        unsigned int                fMaxIterations;
        bool                        fVerbose;
        
    public:
        VLevenbergMarquardt( const VMinimizerFn& fn, const vector<double>& param, const vector<double>& err, double tol = 0.001, unsigned max_iterations = 1000, bool verbose = false );
        virtual void minimize();
};

////////// VSVD: Singular Value Decomposition ////////////////////

class VSVD
{
    public:
        VSVD( const vector< vector< double > >& a );  //! Construct SVD of matrix a.
        void solveSVD( const vector<double>& b, vector<double>& x, double thresh = -1.0 );
        vector< vector< double > > inverse() const; //! Return the pseudoinverse of the input matrix.
        
    private:
        void decompose();
        void reorder();
        double pythag( const double a, const double b );
        void setThresh( double thresh );
        double SIGN( const double a, const double b )
        {
            return ( b >= 0.0 ) ? fabs( a ) : -fabs( a );
        }
        int                        fNRow;
        int                        fNCol;
        vector< vector< double > > fU;
        vector< vector< double > > fV;
        vector<double>             fW; //!< Vector of singular values.
        double                     fEps;
        double                     fThresh;
};

/////////////////////////////////////////////////////////////////

class VMinimizerFactory
{
    private:
        static unique_ptr<VMinimizerFactory> fInstance;
        static string       fDefaultMinimizer;
        static double       fDefaultTolerance;
        static unsigned int fDefaultMaxIterations;
        static bool         fDefaultVerbose;
        
    protected:
        VMinimizerFactory() { }
        
    public:
        virtual ~VMinimizerFactory();
        VMinimizer* getMinimizer( const VMinimizerFn& fn, const vector<double>& param, const vector<double>& err, unsigned int max_iterations = fDefaultMaxIterations, double tolerance = fDefaultTolerance ) const; //! Minimizer factory method to generate default minimizer.
        VMinimizer* getMinimizer( const VMinimizerFn& fn, const vector<double>& param, unsigned int max_iterations = fDefaultMaxIterations, double tolerance = fDefaultTolerance ) const;   //! Minimizer factory method to generate default minimizer.
        static VMinimizerFactory* getInstance();
        static void configure();
};

#endif // VMINIMIZER_H
