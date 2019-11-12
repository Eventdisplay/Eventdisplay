//! VMinimizer (adapted from Matthew Wood's SLAC analysis)

#include "VMinimizer.h"

using namespace std;

// ----------------------------------------------------------------------------
// VMinimizer

VMinimizer::VMinimizer( const VMinimizerFn& fn, const vector<double>& param, const vector<double>& error ): fFn( fn )
{
    fHasLoBound.resize( param.size() );
    fLoBound.resize( param.size() );
    fHasHiBound.resize( param.size() );
    fHiBound.resize( param.size() );
    fParam = param;
    fError = error;
    fCov.resize( param.size() );
    for( unsigned int i = 0; i < param.size(); i++ )
    {
        fCov[i].resize( param.size() );
    }
    //fFnVal()
    //fConverged()
    fNFit = fn.getNParam();
    fFit.resize( fn.getNParam(), true );
}

// ----------------------------------------------------------------------------
// VLevenbergMarquardt

VLevenbergMarquardt::VLevenbergMarquardt( const VMinimizerFn& fn, const vector<double>& param, const vector<double>& err, double tol, unsigned max_iterations, bool verbose ): VMinimizer( fn, param, err )
{
    fCov.resize( fn.getNParam() );
    for( unsigned int i = 0; i < fn.getNParam(); i++ )
    {
        fCov[i].resize( fn.getNParam() );
    }
    fTol = tol;
    fMaxIterations = max_iterations;
    fVerbose = verbose;
}

void VLevenbergMarquardt::initialize()
{
    fVerbose = false; // set true for debug printing
    fLambda = 1E-3;
    fFnVal = 0;
    fNiter = 0;
    setConverged( false );
    
    const unsigned nf = getNFit();
    
    //JG: vectors
    fAlpha.resize( nf );
    fBeta.resize( nf, 0.0 );
    fDa.resize( nf, 0.0 );
    
    for( unsigned int i = 0; i < nf; i++ )
    {
        fAlpha[i].resize( nf );
        for( unsigned int j = 0; j < nf; j++ )
        {
            fAlpha[i][j] = 0.0;
        }
    }
    
    mrqcof( VMinimizer::getParam(), fAlpha , fBeta, fFnVal ); //JG:
    
}

void VLevenbergMarquardt::minimize()
{
    initialize();
    
    double dfnval = 0;
    double fnval = fFnVal;
    
    while( 1 )
    {
    
        if( fNiter >= fMaxIterations )
        {
            break;
        }
        fnval = fFnVal;
        
        iterate();
        
        dfnval = ( fFnVal - fnval );
        if( fVerbose )
        {
            cout << "VMinimizer: iter " << fNiter << " min " << setprecision( 8 ) << fnval << " difmin " << dfnval << " lamba " << fLambda << endl;
        }
        
        if( fNiter >= 2 && fabs( dfnval ) < fTol && dfnval < 0 )
        {
            setConverged( true );
            break;
        }
        
        fNiter++;
    }
    fLambda = 0;
    
    iterate();
    
    setFnVal( fFnVal );
    if( fVerbose && !getConverged() )
    {
        cout << "VLevenbergMarquardt::minimize(): warning not converged" << endl;
    }
}

void VLevenbergMarquardt::iterate()
{
    const unsigned nf = getNFit();
    
    vector< vector< double > > alpha;
    alpha.resize( nf );
    for( unsigned int i = 0; i < nf; i++ )
    {
        alpha[i].resize( nf );
        for( unsigned int j = 0; j < nf; j++ )
        {
            alpha[i][j] = 0.0;
        }
    }
    const unsigned getNParam = VMinimizer::getNParam();
    
    for( unsigned ifit = 0; ifit < nf; ifit++ )
    {
        for( unsigned jfit = 0; jfit < nf; jfit++ )
        {
            alpha[ifit][jfit] = fAlpha[ifit][jfit];
        }
        alpha[ifit][ifit] = fAlpha[ifit][ifit] * ( 1 + fLambda );
    }
    
    vector<double> beta = fBeta;
    
    
    // ----------------------------------------------------------------------
    // Solve for inverse of alpha and solution to normal equations
    VSVD svd( alpha );
    svd.solveSVD( beta, fDa );
    
    // Evaluate covariance matrix -------------------------------------------
    if( fLambda == 0 )
    {
        VSVD svd( alpha );
        vector< vector< double > > cov = svd.inverse();
        
        for( unsigned int ifit = 0; ifit < nf; ifit++ )
            for( unsigned int jfit = 0; jfit < nf; jfit++ )
            {
                fCov[ifit][jfit] = cov[ifit][jfit];
            }
        for( unsigned int iparm = nf; iparm < getNParam; iparm++ )
            for( unsigned int jparm = 0; jparm < iparm + 1; jparm++ )
            {
                fCov[iparm][jparm] = fCov[jparm][iparm] = 0.0;
            }
        unsigned int kparm = nf - 1;
        for( int jparm = getNParam - 1; jparm >= 0; jparm-- )
            if( fFit[jparm] )
            {
                for( unsigned iparm = 0; iparm < getNParam; iparm++ )
                {
                    swap( fCov[iparm][kparm], fCov[iparm][jparm] );
                }
                for( unsigned iparm = 0; iparm < getNParam; iparm++ )
                {
                    swap( fCov[kparm][iparm], fCov[jparm][iparm] );
                }
                kparm--;
            }
            
        setCovariance( fCov );
        
        return;
    }
    
    double fnval = 0;
    vector<double> atry( getNParam );
    for( unsigned int iparm = 0, jparm = 0; iparm < getNParam; iparm++ )
    {
        if( fFit[iparm] )
        {
            atry[iparm] = VMinimizer::getParam( iparm ) + fDa[jparm++];
            if( getHasLoBound( iparm ) && atry[iparm] <= getLoBound( iparm ) )
            {
                atry[iparm] = getLoBound( iparm );
            }
            if( getHasHiBound( iparm ) && atry[iparm] >= getHiBound( iparm ) )
            {
                atry[iparm] = getHiBound( iparm );
            }
        }
        else
        {
            atry[iparm] = VMinimizer::getParam( iparm );
        }
    }
    
    mrqcof( atry, alpha, beta, fnval );
    
    if( fnval < fFnVal )
    {
        fLambda *= 0.1;
        fFnVal = fnval;
        
        fBeta = beta;
        fAlpha = alpha;
        setParam( atry );
    }
    else
    {
        fLambda *= 10.;
    }
}

bool VLevenbergMarquardt::mrqcof( const vector<double>& param, vector< vector< double > >& alpha, vector<double>& beta, double& fnval )
{
    for( unsigned int i = 0; i < alpha.size(); i++ )
    {
        for( unsigned int j = 0; j < alpha[i].size(); j++ )
        {
            alpha[i][j] = 0.0;
        }
    }
    for( unsigned int i = 0; i < beta.size(); i++ )
    {
        beta[i] = 0.0;
    }
    vector< vector< double > > alpha2;
    vector<double> beta2;
    
    VMinimizer::getFn().val( param, fnval, beta2, alpha2 );
    const unsigned NParam = VMinimizer::getNParam();
    
    for( unsigned iparm = 0, ifit = 0; iparm < NParam; iparm++ )
    {
        if( fFit[iparm] )
        {
            for( unsigned jparm = 0, jfit = 0; jparm < NParam; jparm++ )
                if( fFit[jparm] )
                {
                    alpha[ifit][jfit++] = alpha2[iparm][jparm];
                }
            beta[ifit++] = beta2[iparm];
        }
    }
    
    const unsigned nf = getNFit();
    
    // Fill in symmetric side
    for( unsigned ifit = 1; ifit < nf; ifit++ )
        for( unsigned jfit = 0; jfit < ifit; jfit++ )
        {
            alpha[jfit][ifit] = alpha[ifit][jfit];
        }
        
    return true;
}


// ----------------------------------------------------------------------------
// VSVD (from LinearAlgebra)

VSVD::VSVD( const vector< vector< double > >& a ): fNRow( a.size() ), fNCol( a[0].size() ), fU( a ), fV(), fW( fNCol ), fEps( numeric_limits<double>::epsilon() ), fThresh()
{
    fV.resize( fNCol );
    for( int i = 0; i < fNCol; i++ )
    {
        fV[i].resize( fNCol );
    }
    
    decompose();
    reorder();
    setThresh( -1.0 );
}

void VSVD::solveSVD( const vector<double>& b, vector<double>& x, double thresh )
{
    if( b.size() != unsigned( fNRow ) || x.size() != unsigned( fNCol ) )
    {
        return;
    }
    
    vector<double> temp;
    temp.resize( fNCol );
    setThresh( thresh );
    
    for( int jcol = 0; jcol < fNCol; jcol++ )
    {
        double s = 0.0;
        if( fW[jcol] > fThresh )
        {
            for( int irow = 0; irow < fNRow; irow++ )
            {
                s += fU[irow][jcol] * b[irow];
            }
            s /= fW[jcol];
        }
        temp[jcol] = s;
    }
    
    for( int jcol = 0; jcol < fNCol; jcol++ )
    {
        double s = 0.0;
        for( int jjcol = 0; jjcol < fNCol; jjcol++ )
        {
            s += fV[jcol][jjcol] * temp[jjcol];
        }
        x[jcol] = s;
    }
    
}

vector< vector< double > > VSVD::inverse() const
{
    double thresh = fThresh;
    
    vector< vector< double > > v;
    v.resize( fV.size() );
    for( unsigned int i = 0; i < fV.size(); i++ )
    {
        v[i].resize( fV[i].size() );
        for( unsigned int j = 0; j < fV[i].size(); j++ )
        {
            v[i][j] = fV[i][j];
        }
    }
    vector< vector< double > > ut;
    ut.resize( fU.size() );
    for( unsigned int i = 0; i < fU.size(); i++ )
    {
        ut[i].resize( fU[i].size() );
    }
    for( unsigned int i = 0; i < fU.size(); i++ )
    {
        for( unsigned int j = 0; j < fU[i].size(); j++ )
        {
            ut[j][i] = fU[i][j]; // transpose
        }
    }
    vector< vector< double > > w;
    w.resize( fNCol );
    for( int i = 0; i < fNCol; i++ )
    {
        w[i].resize( fNCol );
    }
    
    for( int icol = 0; icol < fNCol; icol++ )
        if( fabs( fW[icol] ) > thresh )
        {
            w[icol][icol] = 1 / fW[icol];
        }
        else
        {
            w[icol][icol] = 0;
        }
        
    for( unsigned int i = 0; i < fV.size(); i++ )
    {
        for( unsigned int j = 0; j < fV[i].size(); j++ )
        {
            v[i][j] *= w[i][j];
            v[i][j] *= ut[i][j];
        }
    }
    
    return v;
}

void VSVD::decompose()
{
    bool flag;
    int i, its, j, jj, k, l = 0, nm = 0;
    double anorm, c, f, g, h, s, scale, x, y, z;
    
    vector<double> rv1;
    rv1.resize( fNCol );
    
    g = scale = anorm = 0.0;
    
    for( i = 0; i < fNCol; i++ )
    {
        l = i + 2;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        
        if( i < fNRow )
        {
            for( k = i; k < fNRow; k++ )
            {
                scale += fabs( fU[k][i] );
            }
            if( scale != 0.0 )
            {
                for( k = i; k < fNRow; k++ )
                {
                    fU[k][i] /= scale;
                    s += fU[k][i] * fU[k][i];
                }
                f = fU[i][i];
                g = -SIGN( sqrt( s ), f );
                
                h = f * g - s;
                fU[i][i] = f - g;
                
                for( j = l - 1; j < fNCol; j++ )
                {
                    for( s = 0.0, k = i; k < fNRow; k++ )
                    {
                        s += fU[k][i] * fU[k][j];
                    }
                    f = s / h;
                    for( k = i; k < fNRow; k++ )
                    {
                        fU[k][j] += f * fU[k][i];
                    }
                }
                for( k = i; k < fNRow; k++ )
                {
                    fU[k][i] *= scale;
                }
            }
        }
        fW[i] = scale * g;
        g = s = scale = 0.0;
        if( i + 1 <= fNRow && i + 1 != fNCol )
        {
            for( k = l - 1; k < fNCol; k++ )
            {
                scale += fabs( fU[i][k] );
            }
            if( scale != 0.0 )
            {
                for( k = l - 1; k < fNCol; k++ )
                {
                    fU[i][k] /= scale;
                    s += fU[i][k] * fU[i][k];
                }
                f = fU[i][l - 1];
                g = -SIGN( sqrt( s ), f );
                h = f * g - s;
                fU[i][l - 1] = f - g;
                for( k = l - 1; k < fNCol; k++ )
                {
                    rv1[k] = fU[i][k] / h;
                }
                for( j = l - 1; j < fNRow; j++ )
                {
                    for( s = 0.0, k = l - 1; k < fNCol; k++ )
                    {
                        s += fU[j][k] * fU[i][k];
                    }
                    for( k = l - 1; k < fNCol; k++ )
                    {
                        fU[j][k] += s * rv1[k];
                    }
                }
                for( k = l - 1; k < fNCol; k++ )
                {
                    fU[i][k] *= scale;
                }
            }
        }
        anorm = max( anorm, ( fabs( fW[i] ) + fabs( rv1[i] ) ) );
    }
    
    for( i = fNCol - 1; i >= 0; i-- )
    {
        if( i < fNCol - 1 )
        {
            if( g != 0.0 )
            {
                for( j = l; j < fNCol; j++ )
                {
                    fV[j][i] = ( fU[i][j] / fU[i][l] ) / g;
                }
                for( j = l; j < fNCol; j++ )
                {
                    for( s = 0.0, k = l; k < fNCol; k++ )
                    {
                        s += fU[i][k] * fV[k][j];
                    }
                    for( k = l; k < fNCol; k++ )
                    {
                        fV[k][j] += s * fV[k][i];
                    }
                }
            }
            for( j = l; j < fNCol; j++ )
            {
                fV[i][j] = fV[j][i] = 0.0;
            }
        }
        fV[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    
    for( i = min( fNRow, fNCol ) - 1; i >= 0; i-- )
    {
        l = i + 1;
        g = fW[i];
        for( j = l; j < fNCol; j++ )
        {
            fU[i][j] = 0.0;
        }
        if( g != 0.0 )
        {
            g = 1.0 / g;
            for( j = l; j < fNCol; j++ )
            {
                for( s = 0.0, k = l; k < fNRow; k++ )
                {
                    s += fU[k][i] * fU[k][j];
                }
                f = ( s / fU[i][i] ) * g;
                for( k = i; k < fNRow; k++ )
                {
                    fU[k][j] += f * fU[k][i];
                }
            }
            for( j = i; j < fNRow; j++ )
            {
                fU[j][i] *= g;
            }
        }
        else for( j = i; j < fNRow; j++ )
            {
                fU[j][i] = 0.0;
            }
        ++fU[i][i];
    }
    
    for( k = fNCol - 1; k >= 0; k-- )
    {
        for( its = 0; its < 30; its++ )
        {
            flag = true;
            for( l = k; l >= 0; l-- )
            {
                nm = l - 1;
                if( l == 0 || fabs( rv1[l] ) <= fEps * anorm )
                {
                    flag = false;
                    break;
                }
                if( fabs( fW[nm] ) <= fEps * anorm )
                {
                    break;
                }
            }
            if( flag )
            {
                c = 0.0;
                s = 1.0;
                for( i = l; i < k + 1; i++ )
                {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if( fabs( f ) <= fEps * anorm )
                    {
                        break;
                    }
                    g = fW[i];
                    h = pythag( f, g );
                    fW[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for( j = 0; j < fNRow; j++ )
                    {
                        y = fU[j][nm];
                        z = fU[j][i];
                        fU[j][nm] = y * c + z * s;
                        fU[j][i] = z * c - y * s;
                    }
                }
            }
            z = fW[k];
            if( l == k )
            {
                if( z < 0.0 )
                {
                    fW[k] = -z;
                    for( j = 0; j < fNCol; j++ )
                    {
                        fV[j][k] = -fV[j][k];
                    }
                }
                break;
            }
            if( its == 29 )
            {
                cout << "VSVD::decompose: warning no convergence after 30 iterations" << endl;
            }
            
            x = fW[l];
            nm = k - 1;
            y = fW[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y );
            g = pythag( f, 1.0 );
            f = ( ( x - z ) * ( x + z ) + h * ( ( y / ( f + SIGN( g, f ) ) ) - h ) ) / x;
            c = s = 1.0;
            for( j = l; j <= nm; j++ )
            {
                i = j + 1;
                g = rv1[i];
                y = fW[i];
                h = s * g;
                g = c * g;
                z = pythag( f, h );
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for( jj = 0; jj < fNCol; jj++ )
                {
                    x = fV[jj][j];
                    z = fV[jj][i];
                    fV[jj][j] = x * c + z * s;
                    fV[jj][i] = z * c - x * s;
                }
                z = pythag( f, h );
                fW[j] = z;
                if( z )
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for( jj = 0; jj < fNRow; jj++ )
                {
                    y = fU[jj][j];
                    z = fU[jj][i];
                    fU[jj][j] = y * c + z * s;
                    fU[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            fW[k] = x;
        }
    }
    
}

void VSVD::reorder()
{
    int i, j, k, s, inc = 1;
    double sw;
    vector<double> su;
    vector<double> sv;
    su.resize( fNRow );
    sv.resize( fNCol );
    
    do
    {
        inc *= 3;
        inc++;
    }
    while( inc <= fNCol );
    do
    {
        inc /= 3;
        for( i = inc; i < fNCol; i++ )
        {
            sw = fW[i];
            for( k = 0; k < fNRow; k++ )
            {
                su[k] = fU[k][i];
            }
            for( k = 0; k < fNCol; k++ )
            {
                sv[k] = fV[k][i];
            }
            j = i;
            while( fW[j - inc] < sw )
            {
                fW[j] = fW[j - inc];
                for( k = 0; k < fNRow; k++ )
                {
                    fU[k][j] = fU[k][j - inc];
                }
                for( k = 0; k < fNCol; k++ )
                {
                    fV[k][j] = fV[k][j - inc];
                }
                j -= inc;
                if( j < inc )
                {
                    break;
                }
            }
            fW[j] = sw;
            for( k = 0; k < fNRow; k++ )
            {
                fU[k][j] = su[k];
            }
            for( k = 0; k < fNCol; k++ )
            {
                fV[k][j] = sv[k];
            }
        }
    }
    while( inc > 1 );
    for( k = 0; k < fNCol; k++ )
    {
        s = 0;
        for( i = 0; i < fNRow; i++ ) if( fU[i][k] < 0. )
            {
                s++;
            }
        for( j = 0; j < fNCol; j++ ) if( fV[j][k] < 0. )
            {
                s++;
            }
        if( s > ( fNRow + fNCol ) / 2 )
        {
            for( i = 0; i < fNRow; i++ )
            {
                fU[i][k] = -fU[i][k];
            }
            for( j = 0; j < fNCol; j++ )
            {
                fV[j][k] = -fV[j][k];
            }
        }
    }
}

double VSVD::pythag( const double a, const double b )
{
    double fabsa = fabs( a ), fabsb = fabs( b );
    return ( fabsa > fabsb ? fabsa * sqrt( 1.0 + std::pow( fabsb / fabsa, 2 ) ) : ( fabsb == 0.0 ? 0.0 : fabsb * sqrt( 1.0 + std::pow( fabsa / fabsb, 2 ) ) ) );
}

void VSVD::setThresh( double thresh )
{
    if( thresh >= 0.0 )
    {
        fThresh = thresh;
    }
    else
    {
        fThresh = 0.5 * sqrt( double( fNRow + fNCol ) + 1.0 ) * fW[0] * fEps;
    }
}

// ----------------------------------------------------------------------------
// VMinimizerFactory

auto_ptr<VMinimizerFactory> VMinimizerFactory::fInstance;

string VMinimizerFactory::fDefaultMinimizer = "lbmq";
double VMinimizerFactory::fDefaultTolerance = 0.02;
unsigned VMinimizerFactory::fDefaultMaxIterations = 100;
bool VMinimizerFactory::fDefaultVerbose = false;

VMinimizerFactory::~VMinimizerFactory()
{
}

VMinimizerFactory* VMinimizerFactory::getInstance()
{
    if( !fInstance.get() )
    {
        fInstance.reset( new VMinimizerFactory );
    }
    return fInstance.get();
}

VMinimizer* VMinimizerFactory::getMinimizer( const VMinimizerFn& fn, const vector<double>& param, unsigned int max_iterations, double tolerance ) const
{
    vector<double> err;
    err.resize( param.size() );
    for( unsigned int i = 0; i < param.size(); i++ )
    {
        err[i] = 0.01 * param[i];
    }
    
    return getMinimizer( fn, param, err, max_iterations, tolerance );
}

VMinimizer* VMinimizerFactory::getMinimizer( const VMinimizerFn& fn, const vector<double>& param, const vector<double>& err, unsigned int max_iterations, double tolerance ) const
{
    return new VLevenbergMarquardt( fn, param, err, tolerance, max_iterations, fDefaultVerbose );
}

void VMinimizerFactory::configure()
{
}
