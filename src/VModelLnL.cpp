//! VModelLnL (adapted from Matthew Wood's SLAC analysis)

#include "VModelLnL.h"

using namespace std;

// ----------------------------------------------------------------------------
// VModelLnL

VModelLnL::VModelLnL( VModel3DData* iVModel3DData, VModel3DFn* iVModel3DFn ): VMinimizerFn( 8 ), fData3D( iVModel3DData ), fModel3DFn( iVModel3DFn )
{
    fSinglePEvar = 0.4 * 0.4; // Hard-coded for VERITAS
    fNParam = 8;  //Model3D number of parameters
    fDebug = false;
}

VModelLnL::~VModelLnL()
{
    if( fLnLFile )
    {
        fLnLFile->Close();
    }
}

void VModelLnL::val( const vector<double>& a, double& lnl, vector<double>& beta, vector< vector< double > >& alpha ) const
{
    fModel3DFn->setModel3DFn( fData3D->flocTel3D, a ); // run once per iteration on model parameters (not needed for each pixel calculation)
    
    lnl = 0;
    beta.resize( fNParam, 0.0 );
    alpha.resize( fNParam );
    for( unsigned int i = 0; i < fNParam; i++ )
    {
        alpha[i].resize( fNParam );
        for( unsigned int j = 0; j < fNParam; j++ )
        {
            alpha[i][j] = 0.0;
        }
    }
    vector<double> scope_beta;
    vector< vector< double > > scope_alpha;
    
    double totGOF = 0; //JG: 3D model
    
    for( unsigned int iscope = 0; iscope < fData3D->fNTel3D; iscope++ )
    {
    
        if( fData3D->fNDFTel3D[iscope] == 0 )
        {
            continue;    //skip telescopes with no selected pixels
        }
        
        double scope_lnl = 0;
        double scopeGOF = 0; //JG: 3D model
        
        scope_beta.resize( fNParam, 0.0 );
        scope_alpha.resize( fNParam );
        for( unsigned int i = 0; i < fNParam; i++ )
        {
            scope_alpha[i].resize( fNParam );
            for( unsigned int j = 0; j < fNParam; j++ )
            {
                scope_alpha[i][j] = 0.0;
            }
        }
        
        scopeLnL( iscope, a, scope_lnl, scope_beta, scope_alpha, scopeGOF );
        
        lnl += scope_lnl;
        totGOF += scopeGOF;
        
        for( unsigned int i = 0; i < fNParam; i++ )
        {
            beta[i] += scope_beta[i];
            for( unsigned int j = 0; j < fNParam; j++ )
            {
                alpha[i][j] += scope_alpha[i][j];
            }
        }
        
    }
    
    // Flip sign of log likelihood
    lnl *= -1;
    for( unsigned int i = 0; i < fNParam; i++ )
    {
        for( unsigned int j = 0; j < fNParam; j++ )
        {
            alpha[i][j] *= -1;
        }
    }
    
    if( fData3D->fNDF3D > 0 )
    {
        fData3D->fGOF3D = totGOF / ( double )fData3D->fNDF3D;    // chi^2
    }
    /// subtract number of parameters and get GOF ///
    /// double totNDF = (double)fData3D->fNDF3D - (double)fNParam;
    /// if( totNDF > 0 ) fData3D->fGOF3D = totGOF / sqrt(2.* totNDF);  // de Naurois
    else
    {
        fData3D->fGOF3D = -99.;
    }
}

void VModelLnL::scopeLnL( unsigned int iscope, const vector<double>& a, double& lnl, vector<double>& beta, vector< vector<double> >& alpha, double& gof ) const
{
    double pixel_nsb = 0;
    double pixel_var = 0;  //use pedvar for Gauss LnL
    
    double scope_chi2 = 0; //JG: chi^2
    //double scope_gof = 0;  //JG: deNaurois
    
    const unsigned int nchan = fData3D->fNpix3D[iscope];
    
    for( unsigned ichan = 0; ichan < nchan; ichan++ )
    {
        /// only include pixels selected through image cleaning ///
        if( fData3D->fClean3D[iscope][ichan] )
        {
            double s = fData3D->fMeasuredSum3D[iscope][ichan];
            double mu = 0;
            pixel_var = fData3D->fPedvar3D[iscope][ichan];
            vector<double> dmuda;
            dmuda.resize( fNParam );
            
            fModel3DFn->calcModel3DFn( iscope, ichan, fData3D->fpX3D[iscope][ichan], fData3D->fpY3D[iscope][ichan], fData3D->fpZ3D[iscope][ichan], mu, dmuda );
            double qterm1 = fData3D->fMarea3D[iscope] * fData3D->fomegapix3D[iscope][ichan] * fData3D->fcosptheta3D[iscope][ichan];
            mu *= qterm1;
            for( unsigned int par = 0; par < fNParam; par++ )
            {
                dmuda[par] *= qterm1;
            }
            //// TEST !!!!
            //if( mu < 0.01 || fabs( mu - s ) > 900 ) //ORIG
            if( mu < 0.1 || fabs( mu - s ) > 900 ) //TEST
            {
                mu = 0.;    /// Poisson set small mu to 0
            }
            
            double pix_lnl = 0;
            double dlnl = 0;
            double d2lnl = 0;
            // double pixGOF = 0; //JG Gaussian
            
            ///pixelLnLGauss(s,pixel_var,pixel_nsb,mu,pix_lnl,dlnl,d2lnl ); //JG: Gauss
            pixelLnL_lookup( s, pixel_var, pixel_nsb, mu, pix_lnl, dlnl, d2lnl ); //JG: Poisson
            /////pixelLnL(s,pixel_var,pixel_nsb,fSinglePEvar,mu,pix_lnl,dlnl,d2lnl); //JG: Poisson
            
            if( TMath::IsNaN( mu ) )
            {
                cout << "Mu NaN: T" << iscope + 1 << " ch" << ichan << " mu " << mu << " s " << s << " var " << pixel_var << " lnl " << pix_lnl << " dlnl " << dlnl << " d2lnl " << d2lnl << endl;
                pix_lnl = 0;
                continue;
            }
            
            if( !isfinite( pix_lnl ) )
            {
                cout << "LnL notfinite: T" << iscope + 1 << " ch" << ichan << " mu " << mu << " s " << s << " var " << pixel_var << " lnl " << pix_lnl << " dlnl " << dlnl << " d2lnl " << d2lnl << endl;
                pix_lnl = 0;
                continue;
            }
            lnl += pix_lnl;
            
            //// GOF ////
            scope_chi2 += pow( s - mu, 2 ) / ( mu + pixel_nsb + pixel_var ); //JG: chi^2
            // double lnl_exp = 0;
            // /// if( mu <= 0 ) logmu = m_lnl_exp->GetXaxis()->GetXmin();
            // if( mu > 0 ) lnl_exp = m_lnl_exp->Interpolate( log10(mu) );
            // else lnl_exp = 0.;
            // scope_gof += (-2.*pix_lnl) - (-lnl_exp); //JG: deNaurois
            // double exnoise = 0.4; //IMPLEMENT THIS !!!!
            // double munoise = pixel_var*pixel_var + mu*(1 + exnoise*exnoise);
            // double Lexpect = 1. + log(2.*TMath::Pi()) + log(munoise);
            // scope_gof += (-2.*pix_lnl) - Lexpect;
            ///////////////
            
            for( unsigned int i = 0; i < fNParam; i++ )
            {
                beta[i] += dmuda[i] * dlnl;
            }
            for( unsigned ip = 0; ip < fNParam; ip++ )
            {
                for( unsigned jp = 0; jp < fNParam; jp++ )
                {
                    alpha[ip][jp] += dmuda[ip] * dmuda[jp] * d2lnl;
                }
            }
            
        } //image cleaning
    }
    
    gof = scope_chi2; //JG chi_2
    //gof = scope_gof;  //JG: deNaurois
    ///  m_scope_lnl[iscope] = lnl;
    return;
}

void VModelLnL::pixelLnLGauss( double s, double pixel_var, double pixel_nsb, double mu, double& lnl, double& dlnl, double& d2lnl ) const
{
    if( mu <= 0 )
    {
        lnl = -0.5 * log( 2.*TMath::Pi() * pixel_var ) - pow( s, 2 ) / ( 2 * pixel_var );
        dlnl = 0;
        d2lnl = 0;
        ///double Lexpect = 1. + log(2.*TMath::Pi()) + log(pixel_var*pixel_var);
        ///pixGOF = (-2.*lnl) - Lexpect;
        return;
    }
    
    //// TEST: Gaussian test (good for high signal) /////////
    double exnoise = 0.4; //TEST: 40% resolution from HESS
    double ped = pixel_var; // use pedvar
    
    double smu = s - mu;
    double munoise = ped * ped + mu * ( 1 + exnoise * exnoise );
    double pdf = exp( -( smu * smu ) / ( 2 * munoise ) ) / sqrt( 2 * TMath::Pi() * munoise );
    
    ///// NOTE: to match Wood don't multiply by 2 and keep positive ////
    lnl = log( pdf );
    
    //// GOF ////
    //double Lexpect = 1. + log(2.*TMath::Pi()) + log(munoise);
    //pixGOF = (-2.*lnl) - Lexpect;
    ///cout<<"pixLnL "<< lnl <<" pixLex "<< Lexpect <<" pixGOF "<< pixGOF <<endl; //TEST
    
    //// first and second derivatives lnl(mu), positive not times 2 ////
    dlnl = exp( ( smu * smu ) / ( 2 * munoise ) ) * sqrt( 2 * TMath::Pi() * munoise ) * ( -( ( exp( -( smu * smu ) / ( 2 * munoise ) ) * ( 1 + exnoise * exnoise ) ) / ( 2 * pow( munoise, 1.5 ) * sqrt( 2 * TMath::Pi() ) ) ) + ( ( exp( -( smu * smu ) / ( 2 * munoise ) ) * ( ( smu / munoise ) + ( ( 1 + exnoise * exnoise ) * smu * smu ) / ( 2 * munoise * munoise ) ) ) / sqrt( 2 * TMath::Pi() * munoise ) ) );
    
    d2lnl = ( ( exp( pow( smu, 2 ) / ( 2.*munoise ) ) * ( 1 + pow( exnoise, 2 ) ) * sqrt( TMath::Pi() / 2. ) * ( -( 1 + pow( exnoise, 2 ) ) / ( 2.*exp( pow( smu, 2 ) / ( 2.*munoise ) ) * pow( munoise, 1.5 ) * sqrt( 2 * TMath::Pi() ) ) + ( ( smu / munoise ) + ( ( 1 + pow( exnoise, 2 ) ) * pow( smu, 2 ) ) / ( 2.*pow( munoise, 2 ) ) ) / ( exp( pow( smu, 2 ) / ( 2.*munoise ) ) * sqrt( munoise ) * sqrt( 2 * TMath::Pi() ) ) ) ) / sqrt( munoise ) ) + exp( pow( smu, 2 ) / ( 2.*munoise ) )   * sqrt( munoise ) * sqrt( 2 * TMath::Pi() ) * ( -( ( smu / munoise ) ) - ( ( 1 + pow( exnoise, 2 ) ) * pow( smu, 2 ) ) / ( 2.*pow( munoise, 2 ) ) ) * ( -( 1 + pow( exnoise, 2 ) ) / ( 2.*exp( pow( smu, 2 ) / ( 2.*munoise ) ) * pow( munoise, 1.5 ) * sqrt( 2 * TMath::Pi() ) ) + ( ( smu / munoise ) + ( ( 1 + pow( exnoise, 2 ) ) * pow( smu, 2 ) ) / ( 2.*pow( munoise, 2 ) ) ) / ( exp( pow( smu, 2 ) / ( 2.*munoise ) ) * sqrt( munoise ) * sqrt( 2 * TMath::Pi() ) ) ) + exp( pow( smu, 2 ) / ( 2.*munoise ) ) * sqrt( munoise ) * sqrt( 2 * TMath::Pi() ) * ( ( 3 * pow( 1 + pow( exnoise, 2 ), 2 ) ) / ( 4.*exp( pow( smu, 2 ) / ( 2.*munoise ) ) * pow( munoise, 2.5 ) * sqrt( 2 * TMath::Pi() ) ) + ( -( 1 / munoise ) - ( 2 * ( 1 + pow( exnoise, 2 ) ) * smu ) / pow( munoise, 2 ) - ( pow( 1 + pow( exnoise, 2 ), 2 ) * pow( smu, 2 ) ) / pow( munoise, 3 ) ) / ( exp( pow( smu, 2 ) / ( 2.*munoise ) ) * sqrt( munoise ) * sqrt( 2 * TMath::Pi() ) ) - ( ( 1 + pow( exnoise, 2 ) ) * ( ( smu / munoise ) + ( ( 1 + pow( exnoise, 2 ) ) * pow( smu, 2 ) ) / ( 2.*pow( munoise, 2 ) ) ) ) / ( exp( pow( smu, 2 ) / ( 2.*( munoise ) ) ) * pow( munoise, 1.5 ) * sqrt( 2 * TMath::Pi() ) ) + pow( ( smu / munoise ) + ( ( 1 + pow( exnoise, 2 ) ) * pow( smu, 2 ) ) / ( 2.*pow( munoise, 2 ) ), 2 ) / ( exp( pow( smu, 2 ) / ( 2.*( munoise ) ) ) * sqrt( munoise ) * sqrt( 2 * TMath::Pi() ) ) );
    
}

void VModelLnL::pixelLnL( double s, double pixel_var, double pixel_nsb, double singlepe_var, double mu, double& lnl, double& dlnl, double& d2lnl )
{

    mu += pixel_nsb; //where does NSB come from (setting to zero)?
    s += pixel_nsb;  //where does NSB come from (setting to zero)?
    
    if( mu <= 0 )
    {
        double var = pixel_var;
        lnl = -0.5 * log( 2. * TMath::Pi() * var ) - pow( s, 2 ) / ( 2 * var );
        dlnl = 0;
        d2lnl = 0;
        return;
    }
    
    int ncenter = min( ( int )lround( mu ), ( int )lround( s ) );
    
    while( 1 )
    {
        int j = ncenter;
        double v2 = pixel_var + j * singlepe_var;
        double dldn = - 0.5 * singlepe_var * pow( v2, -1. ) + singlepe_var * pow( s - j, 2 ) / ( 2 * v2 * v2 ) + ( s - j ) / v2 + log( mu ) - log( j ) - 1. / ( 2 * j );
        
        if( dldn < 0 )
        {
            break;
        }
        
        ncenter++;
    }
    
    double var0 = pixel_var + ncenter * singlepe_var;
    double nvar = max( var0, ( double )ncenter );
    
    int nmin = max( 0, ( int )lround( ncenter - 5 * sqrt( nvar ) ) );
    int nmax = max( 5, ( int )lround( ncenter + 5 * sqrt( nvar ) ) ) + 1;
    
    double lsum = ncenter * log( mu ) - lgamma( ncenter + 1 ) - mu - pow( s - ( double )ncenter, 2 ) / ( 2 * var0 );
    
    const double lgam = lgamma( ncenter + 1 );
    
    double lmu = 0;
    double dldmu = 0;
    double d2ldmu = 0;
    
    for( int n = nmin; n < nmax; n++ )
    {
        const double var = pixel_var + n * singlepe_var;
        const double exp0 = -pow( s - ( double )n, 2 ) / ( 2 * var ) + pow( s - ( double )ncenter, 2 ) / ( 2 * var0 );
        const double exp1 = -lgamma( n + 1 ) + lgam;
        const double exp2 = ( n - ncenter ) * log( mu );
        const double gfact = 1. / sqrt( 2.*TMath::Pi() * var );
        double w = gfact * exp( exp0 + exp1 + exp2 );
        
        lmu += w;
        dldmu += ( n / mu - 1.0 ) * w;
        d2ldmu += ( 1 - 2 * n / mu + ( n * n - n ) / ( mu * mu ) ) * w;
    }
    
    lnl = lsum + log( lmu );
    dlnl = ( dldmu / lmu );
    d2lnl = ( d2ldmu / lmu - dldmu * dldmu / ( lmu * lmu ) );
    
}

void VModelLnL::pixelLnL_lookup( double s, double pixel_var, double pixel_nsb, double mu, double& lnl, double& dlnl, double& d2lnl ) const
{
    if( mu <= 0 )
    {
        pixelLnL( s, pixel_var, pixel_nsb, fSinglePEvar, mu, lnl, dlnl, d2lnl );
        return;
    }
    else if( s < m_lnl_smu->GetXaxis()->GetXmax() && log10( mu ) < m_lnl_smu->GetYaxis()->GetXmax() )
    {
        if( s <= -5 )
        {
            s = -4.9;
        }
        lnl = m_lnl_smu->Interpolate( s, log10( mu ) );
        dlnl = m_dlnl_smu->Interpolate( s, log10( mu ) );
        d2lnl = m_d2lnl_smu->Interpolate( s, log10( mu ) );
        return;
    }
    else if( s >= m_lnl_smu->GetXaxis()->GetXmax() && log10( mu ) < m_lnl_lmu->GetYaxis()->GetXmax() )
    {
        if( s <= -3 )
        {
            s = -2.9;
        }
        lnl = m_lnl_lmu->Interpolate( log10( s ), log10( mu ) );
        dlnl = m_dlnl_lmu->Interpolate( log10( s ), log10( mu ) );
        d2lnl = m_d2lnl_lmu->Interpolate( log10( s ), log10( mu ) );
        return;
    }
    else
    {
        if( fDebug )
        {
            cout << "VModelLnL: else s " << s << " mu " << mu << endl;
        }
        pixelLnL( s, pixel_var, pixel_nsb, fSinglePEvar, mu, lnl, dlnl, d2lnl );
        return;
    }
}

void VModelLnL::readLnLTable( string LnLTableFile )
{
    //// read in Lookup table and get histograms ////
    cout << "LnL Table: " << LnLTableFile.c_str() << endl;
    
    fLnLFile = new TFile( LnLTableFile.c_str() , "READ" );
    
    if( fLnLFile->IsZombie() )
    {
        cout << "VModelLnL::readLnLTable error: File " << LnLTableFile.c_str() << " does not exist!" << endl;
        exit( -1 );
    }
    
    m_lnl_lmu = ( TH2D* )fLnLFile->Get( "m_lnl_lmu" );
    m_dlnl_lmu = ( TH2D* )fLnLFile->Get( "m_dlnl_lmu" );
    m_d2lnl_lmu = ( TH2D* )fLnLFile->Get( "m_d2lnl_lmu" );
    m_lnl_smu = ( TH2D* )fLnLFile->Get( "m_lnl_smu" );
    m_dlnl_smu = ( TH2D* )fLnLFile->Get( "m_dlnl_smu" );
    m_d2lnl_smu = ( TH2D* )fLnLFile->Get( "m_d2lnl_smu" );
    
    m_lnl_exp = ( TH1D* )fLnLFile->Get( "m_lnl_exp" );
    
    //NOTE: fLnLFile is closed in the destructor
}


void VModelLnL::createLnLTable( double pixel_var, double pixel_nsb )
{
    ////// Create and fill Lookup Table ///////
    TFile* fTableFile = new TFile( "table_LnL.root", "RECREATE" ); // create output file
    cout << "creating table_LnL.root" << endl;
    if( fTableFile->IsZombie() )
    {
        cout << "VModelLnL::createLnLTable error creating table, exiting" << endl;
        exit( 0 );
    }
    
    const int nbinMuL = 1200; // number of bins in log(mu) in large range
    const int nbinSiL = 400;  // number of bins in log(signal)
    double logMuLoL = -3; // lowest log(mu) in large range
    double logMuHiL = 4;  // highest log(mu) in large range
    double logSiLo = 1;  // lowest log(signal) in large range
    double logSiHi = 4;  // highest log(signal) in large range
    
    const int nbinMuS = 1000; // number of bins in log(mu) in small range
    const int nbinSiS = 400;  // number of bins in signal in small range
    double logMuLoS = -3; // lowest log(mu) in small range
    double logMuHiS = 3;  // highest log(mu) in small range
    double SiLo = -5;  // lowest signal in small range
    double SiHi = 10;  // highest signal in small range
    
    m_lnl_lmu = new TH2D( "m_lnl_lmu", "m_lnl_lmu", nbinSiL, logSiLo, logSiHi, nbinMuL, logMuLoL, logMuHiL );
    m_dlnl_lmu = new TH2D( "m_dlnl_lmu", "m_dlnl_lmu", nbinSiL, logSiLo, logSiHi, nbinMuL, logMuLoL, logMuHiL );
    m_d2lnl_lmu = new TH2D( "m_d2lnl_lmu", "m_d2lnl_lmu", nbinSiL, logSiLo, logSiHi, nbinMuL, logMuLoL, logMuHiL );
    
    m_lnl_smu = new TH2D( "m_lnl_smu", "m_lnl_smu", nbinSiS, SiLo, SiHi, nbinMuS, logMuLoS, logMuHiS );
    m_dlnl_smu = new TH2D( "m_dlnl_smu", "m_dlnl_smu", nbinSiS, SiLo, SiHi, nbinMuS, logMuLoS, logMuHiS );
    m_d2lnl_smu = new TH2D( "m_d2lnl_smu", "m_d2lnl_smu", nbinSiS, SiLo, SiHi, nbinMuS, logMuLoS, logMuHiS );
    
    m_lnl_exp = new TH1D( "m_lnl_exp", "m_lnl_exp", nbinMuL, logMuLoL, logMuHiL );
    
    double lnl, dlnl, d2lnl = 0;
    double s, mu = 0;
    
    cout << "VModelLnL Computing LnL Lookup Table (Large Range)" << endl;
    
    for( int iSiL = 0; iSiL < nbinSiL; iSiL++ )
    {
        for( int iMuL = 0; iMuL < nbinMuL; iMuL++ )
        {
            s = m_lnl_lmu->GetXaxis()->GetBinCenter( iSiL + 1 );
            mu = m_lnl_lmu->GetYaxis()->GetBinCenter( iMuL  + 1 );
            s = pow( 10, s );
            mu = pow( 10, mu );
            
            pixelLnL( s, pixel_var, pixel_nsb, fSinglePEvar, mu, lnl, dlnl, d2lnl );
            
            m_lnl_lmu->SetBinContent( iSiL + 1, iMuL + 1, lnl );
            m_dlnl_lmu->SetBinContent( iSiL + 1, iMuL + 1, dlnl );
            m_d2lnl_lmu->SetBinContent( iSiL + 1, iMuL + 1, d2lnl );
        }
    }
    
    cout << "VModelLnL Computing LnL Lookup Table (Small Range)" << endl;
    
    for( int iSiS = 0; iSiS < nbinSiS; iSiS++ )
    {
        for( int iMuS = 0; iMuS < nbinMuS; iMuS++ )
        {
            s = m_lnl_smu->GetXaxis()->GetBinCenter( iSiS + 1 );
            mu = m_lnl_smu->GetYaxis()->GetBinCenter( iMuS  + 1 );
            mu = pow( 10, mu );
            
            pixelLnL( s, pixel_var, pixel_nsb, fSinglePEvar, mu, lnl, dlnl, d2lnl );
            
            m_lnl_smu->SetBinContent( iSiS + 1, iMuS + 1, lnl );
            m_dlnl_smu->SetBinContent( iSiS + 1, iMuS + 1, dlnl );
            m_d2lnl_smu->SetBinContent( iSiS + 1, iMuS + 1, d2lnl );
        }
    }
    
    cout << "VModelLnL Computing LnL Expectation Table" << endl;
    
    // Compute expecation value of likelihood m_lnl_exp
    
    for( int iMuL = 0; iMuL < nbinMuL; iMuL++ )
    {
        mu = m_lnl_exp->GetBinCenter( iMuL + 1 );
        mu = pow( 10, mu );
        
        double step = max( sqrt( mu ) / 10., sqrt( pixel_var ) );
        double smin = mu - 100 * step;
        double smax = mu + 100 * step;
        
        double lnl_exp = 0;
        
        for( double s = smin; s < smax; s += step )
        {
        
            pixelLnL( s, pixel_var, pixel_nsb, fSinglePEvar, mu, lnl, dlnl, d2lnl );
            
            lnl_exp += lnl * exp( lnl );
        }
        
        lnl_exp *= step;
        
        m_lnl_exp->SetBinContent( iMuL + 1, lnl_exp );
    }
    
    if( fTableFile )
    {
        cout << endl << "writing data to " << fTableFile->GetName() << endl;
        fTableFile->cd();
        m_lnl_lmu->Write();
        m_dlnl_lmu->Write();
        m_d2lnl_lmu->Write();
        m_lnl_smu->Write();
        m_dlnl_smu->Write();
        m_d2lnl_smu->Write();
        m_lnl_exp->Write();
    }
    if( fTableFile )
    {
        fTableFile->Close();
    }
    
}
