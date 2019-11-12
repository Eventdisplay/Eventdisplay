/*! \class VModel3DFn
     based on Lemoine-Goumard et al. 2006
     adapted by J. Grube and G. Gyuk
*/

#include "VModel3DFn.h"

VModel3DFn::VModel3DFn()
{
}

VModel3DFn::~VModel3DFn()
{
}

void VModel3DFn::setModel3DFn( const vector< vector< double > >& locTel3D, const vector<double>& param )
{
    flocTel3D = locTel3D;
    
    /// read in model parameters ///
    if( param.size() == 8 )
    {
        fSel     = param[0];
        fSaz     = param[1];
        fXcore   = param[2];
        fYcore   = param[3];
        fSmax    = param[4];
        fsigmaL  = param[5];
        fsigmaT  = param[6];
        fNc      = param[7];
    }
    else
    {
        cout << "MODEL PARAMETERS DON'T MATCH" << endl;
        return;
    }
    /// azimuth from -180 to 180 ///
    if( fSaz < -180. )
    {
        fSaz += 360.;
    }
    if( fSaz >  180. )
    {
        fSaz -= 360.;
    }
    
    /// scale fSmax and fsigmaL by 1e3 ///
    fSmax *= 1000.;
    fsigmaL *= 1000.;
    
    ////////////////////////////////////////////////
    fNdim3D = locTel3D.size();    // x,y,z
    fNTel3D = locTel3D[0].size(); // same number of telescopes for x,y,z
    /// geometry for model fit ///
    fs.resize( fNdim3D, 0 );
    fxBo.resize( fNdim3D, 0 );
    fxB.resize( fNdim3D );
    for( unsigned int i = 0; i < fNdim3D; i++ )
    {
        fxB[i].resize( fNTel3D );
        for( unsigned int j = 0; j < fNTel3D; j++ )
        {
            fxB[i][j] = 0;
        }
    }
    
    /// get "s" vector ///
    double tSzen = 0; // temp S zenith
    double tSaz = 0; // temp S azimuth
    tSzen = ( 90. - fSel ) * ( TMath::Pi() / 180. ); // in radians
    tSaz = fSaz * ( TMath::Pi() / 180. ); // in radians
    
    /// from spherical to cartesian (ground coordinates) ///
    fs[0] = sin( tSzen ) * cos( tSaz );
    fs[1] = sin( tSzen ) * sin( tSaz );
    fs[2] = cos( tSzen );
    
    /// calculate xB using height of shower max and core location ///
    /// apply magnitude (shower max) to get shower max point ///
    fxBo[0] = fSmax * fs[0];
    fxBo[1] = fSmax * fs[1];
    fxBo[2] = fSmax * fs[2];
    
    /// move shower to land at core (on ground) relative to (0,0,0) ///
    fxBo[0] = fXcore + fxBo[0];
    fxBo[1] = fYcore + fxBo[1];
    fxBo[2] = 0. + fxBo[2];
    
    /// telescope xB ///
    for( unsigned int i = 0; i < fNdim3D; i++ )
    {
        for( unsigned int j = 0; j < fNTel3D; j++ )
        {
            fxB[i][j] = fxBo[i] - locTel3D[i][j];
        }
    }
    
}


void VModel3DFn::calcModel3DFn( unsigned int iTel, unsigned int iPix, const double ipX, const double ipY, const double ipZ, double& val, vector<double>& grad ) const
{
    double tSzen = ( 90. - fSel ) * ( TMath::Pi() / 180. ); // in radians
    double tSaz = fSaz * ( TMath::Pi() / 180. ); // in radians
    
    ////////////////////////////////////////////////////////////////
    //// Lemoine-Goumard 2006, equation 5 //////////////////////////
    
    //// get Bs, Bp, u, epsilon, deltaB, sigmaU, sigmaD ///////
    
    double Bs = dot( fxB[0][iTel], fxB[1][iTel], fxB[2][iTel], fs[0], fs[1], fs[2] );
    double Bp = dot( fxB[0][iTel], fxB[1][iTel], fxB[2][iTel], ipX, ipY, ipZ ); // CHECK fp
    
    double u = dot( fs[0], fs[1], fs[2], ipX, ipY, ipZ );
    
    double epsilon = acos( u ); // in radians
    
    ///////fepsilon[iTel][iPix] = acos( u ); // in radians
    ///Note: deltaB gives nan error (negative in sqrt) ///
    // use sqare of deltaB //
    
    double deltaBsq = ( ( fxB[0][iTel] * fxB[0][iTel] ) + ( fxB[1][iTel] * fxB[1][iTel] ) + ( fxB[2][iTel] * fxB[2][iTel] ) ) - ( Bp * Bp );
    
    double sigmaU = sqrt( ( fsigmaT * fsigmaT * u * u ) + ( fsigmaL * fsigmaL * ( 1 - ( u * u ) ) ) );
    
    double sigmaD = sqrt( ( fsigmaL * fsigmaL ) - ( fsigmaT * fsigmaT ) );
    
    //// Eq.5 (first term) ///////
    double Nc = exp( fNc ); // get Nc from ln(Nc)
    
    ///// Don't need freq term /////
    ///double Cx = -( ( (fsigmaL*fsigmaL *Bp) - (sigmaD*sigmaD *u *Bs) ) / (sigmaU * fsigmaT * fsigmaL ) );
    ///double C = 1. - freq( Cx ); // very close to 1.0
    ///double Eq5t1 = (Nc * C) / (2.*TMath::Pi() * sigmaU * fsigmaT);
    
    double Eq5t1 = ( Nc ) / ( 2.*TMath::Pi() * sigmaU * fsigmaT );
    
    //// Eq.5 (second term) ///////
    double Eq5t2p1 = deltaBsq / ( fsigmaT * fsigmaT );
    double Eq5t2p2 = ( sigmaD * sigmaD ) / ( fsigmaT * fsigmaT * sigmaU * sigmaU );
    double Eq5t2p3 = ( ( u * Bp ) - Bs ) * ( ( u * Bp ) - Bs );
    
    double Eq5t2 = Eq5t2p1 - ( Eq5t2p2 * Eq5t2p3 );
    
    ////  Eq.5 ///////
    double Eq5 = Eq5t1 * exp( -0.5 * Eq5t2 );
    
    //////////////////////////////////////////////////////////////////////////
    //// I(epsilon) probability of Cherenkov photon per unit solid angle /////
    /// Lemoine ///
    double eta = 0.015 * sqrt( cos( tSzen ) );
    
    double Ie = 0; // I(epsilon)
    double K = 1. / ( 9 * TMath::Pi() * ( eta * eta ) ); // normalization
    
    if( epsilon <= eta )
    {
        Ie = K;
    }
    else if( epsilon > eta )
    {
        Ie = K * ( eta / epsilon ) * exp( -( ( epsilon - eta ) / ( 4 * eta ) ) );
    }
    
    ////////////////////////////////////////////
    
    val = Ie * Eq5;
    
    /////// get partial derivatives //////////////////////////////////
    
    double dEq5dSel = -( exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *   pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +  pow( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +  pow( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -  pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +    ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) +    ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) * ( -2 * pow( fsigmaL, 2 ) * ( -( ( ipX * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) +   2 * pow( fsigmaT, 2 ) * ( -( ( ipX * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) ) ) / ( 4.*fsigmaT * TMath::Pi() * pow( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 1.5 ) ) - ( 0.07957747154594767 * exp( fNc -  0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *  pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) * ( ( ( 2 * fSmax * TMath::Pi() * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) *  sin( tSzen ) ) / 180. - ( 2 * fSmax * TMath::Pi() * cos( tSzen ) * cos( tSaz ) * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) ) / 180. - ( 2 * fSmax * TMath::Pi() * cos( tSzen ) * sin( tSaz ) * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ) / 180. - 2 * ( -( ( ipX * fSmax * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) /   180. ) + ( ipZ * fSmax * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * fSmax * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) / pow( fsigmaT, 2 )  - ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -( ( ( flocTel3D[0][iTel] - fXcore ) * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) /  180. ) + ( flocTel3D[2][iTel] * TMath::Pi() * sin( tSzen ) ) / 180. - ( flocTel3D[1][iTel] * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. + ( fYcore * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. + ( -( ( ipX * fSmax * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) / 180. ) + ( ipZ * fSmax * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * fSmax * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ) ) + ( -( ( ipX * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) * ( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +  flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -  fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *  pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) + ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -2 * pow( fsigmaL, 2 ) * ( -( ( ipX * TMath::Pi() * cos( tSzen ) *  cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ) ) +  2 * pow( fsigmaT, 2 ) * ( -( ( ipX * TMath::Pi() * cos( tSzen ) *  cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ) ) ) *    pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +  flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -  fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * pow( pow( fsigmaT, 2 ) *  pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 2 ) ) ) ) / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dSaz =  -( exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *  pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *     pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +     pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore +      fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore +      fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) * ( -2 * pow( fsigmaL, 2 ) * ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) +    2 * pow( fsigmaT, 2 ) * ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) ) ) / ( 4.*fsigmaT * TMath::Pi() * pow( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +      ipX * cos( tSaz ) * sin( tSzen ) +      ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 1.5 ) ) - ( 0.07957747154594767 * exp( fNc -   0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +    flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -    fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *    pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +    pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore +     fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore +     fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) /      pow( fsigmaT, 2 ) ) ) * ( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -2 * pow( fsigmaL, 2 ) * ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) + 2 * pow( fsigmaT, 2 ) * ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) ) *     pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore +     fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore +     fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * pow( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 2 ) ) - ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore +     fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore +     fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) * ( ( flocTel3D[1][iTel] * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( fYcore * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ( flocTel3D[0][iTel] - fXcore ) * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. + ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ( ipY * fSmax * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * fSmax * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) + ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore +     fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore +     fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) + ( ( -2 * fSmax * TMath::Pi() * sin( tSzen ) * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) * sin( tSaz ) ) / 180. + ( 2 * fSmax * TMath::Pi() * cos( tSaz ) * sin( tSzen ) * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ) /      180. - 2 * ( ( ipY * fSmax * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * fSmax * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) / pow( fsigmaT, 2 ) ) )      / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +      ipX * cos( tSaz ) * sin( tSzen ) +      ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dXcore = ( -0.07957747154594767 * exp( fNc -  0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *   pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +   pow( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +   pow( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -   pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) * ( ( 2 * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) - 2 * ipX * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) / pow( fsigmaT, 2 ) - ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -( cos( tSaz ) * sin( tSzen ) ) +  ipX * ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) ) * ( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +  flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -  fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *  pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) ) / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dYcore = ( -0.07957747154594767 * exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *   pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *   pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +  pow( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +  pow( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -  pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +    ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) +    ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) /    pow( fsigmaT, 2 ) ) ) * ( ( 2 * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) -    2 * ipY * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) / pow( fsigmaT, 2 ) - ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -( sin( tSzen ) * sin( tSaz ) ) + ipY * ( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ) ) ) * ( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +    ipX * ( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ) ) +    ipY * ( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) ) / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dSmax = ( -0.07957747154594767 * exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *   pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +  flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -  fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +     ipX * cos( tSaz ) * sin( tSzen ) +     ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +     ipX * ( -flocTel3D[0][iTel] + fXcore +        fSmax * cos( tSaz ) * sin( tSzen ) ) +     ipY * ( -flocTel3D[1][iTel] + fYcore +        fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *  pow( ipZ * cos( tSzen ) +     ipX * cos( tSaz ) * sin( tSzen ) +     ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 -     pow( ipZ * cos( tSzen ) +       ipX * cos( tSaz ) * sin( tSzen ) +       ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +  pow( -flocTel3D[0][iTel] + fXcore +    fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +  pow( -flocTel3D[1][iTel] + fYcore +    fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -  pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +    ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) +    ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) /            pow( fsigmaT, 2 ) ) ) * ( ( 2 * cos( tSzen ) * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +            2 * cos( tSaz ) * sin( tSzen ) * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +            2 * sin( tSzen ) * sin( tSaz ) * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) -            2 * ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) / pow( fsigmaT, 2 ) - ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( -1 + pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) * ( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +    ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) +    ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ) ) ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) ) / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +            ipX * cos( tSaz ) * sin( tSzen ) +            ipY * sin( tSzen ) * sin( tSaz ), 2 ) +         pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dsigmaL =  -( exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *  pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +    ipX * ( -flocTel3D[0][iTel] + fXcore +       fSmax * cos( tSaz ) * sin( tSzen ) ) +    ipY * ( -flocTel3D[1][iTel] + fYcore +       fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 -    pow( ipZ * cos( tSzen ) +      ipX * cos( tSaz ) * sin( tSzen ) +      ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) )  + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) /  pow( fsigmaT, 2 ) ) ) * fsigmaL * ( 1 -   pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) / ( 2.*fsigmaT * TMath::Pi() * pow( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 1.5 ) ) - ( 0.07957747154594767 * exp( fNc -  0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +    flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -    fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +   ipX * ( -flocTel3D[0][iTel] + fXcore +      fSmax * cos( tSaz ) * sin( tSzen ) ) +   ipY * ( -flocTel3D[1][iTel] + fYcore +      fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *    pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) +    pow( fsigmaL, 2 ) * ( 1 -   pow( ipZ * cos( tSzen ) +     ipX * cos( tSaz ) * sin( tSzen ) +     ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +    pow( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +    pow( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -    pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) * ( ( 2 * fsigmaL * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) *    pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * pow( pow( fsigmaT, 2 ) *   pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +    ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 2 ) ) - ( 2 * fsigmaL * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *   pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +    ipX * cos( tSaz ) * sin( tSzen ) +     ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) ) / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dsigmaT = -( exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +  flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -  fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) *  pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) )  + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) /  pow( fsigmaT, 2 ) ) ) * pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) / ( 2.*TMath::Pi() * pow( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 1.5 ) ) - exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) *   pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) +  ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) +  ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 -  pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +  pow( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +  pow( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -  pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) / ( 2.*pow( fsigmaT, 2 ) * TMath::Pi() * sqrt( pow( fsigmaT, 2 ) *  pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) - ( 0.07957747154594767 * exp( fNc -  0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +   fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +   fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) * ( ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +   ipY * sin( tSzen ) * sin( tSaz ), 2 ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( fsigmaT * pow( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ), 2 ) ) + ( 2 * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( fsigmaT * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) + ( 2 * ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) +   flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) -   fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 3 ) * ( pow( fsigmaT, 2 ) *   pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +   pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) - ( 2 * ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) +   pow( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ), 2 ) +   pow( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ), 2 ) -   pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) ) /  pow( fsigmaT, 3 ) ) ) / ( fsigmaT * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) +  pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    double dEq5dNc = exp( fNc - 0.5 * ( -( ( ( pow( fsigmaL, 2 ) - pow( fsigmaT, 2 ) ) * pow( -fSmax + flocTel3D[2][iTel] * cos( tSzen ) + ( flocTel3D[0][iTel] - fXcore ) * cos( tSaz ) * sin( tSzen ) + flocTel3D[1][iTel] * sin( tSzen ) * sin( tSaz ) - fYcore * sin( tSzen ) * sin( tSaz ) + ( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * ( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ) ) ), 2 ) ) / ( pow( fsigmaT, 2 ) * ( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) ) ) + ( pow( -flocTel3D[2][iTel] + fSmax * cos( tSzen ), 2 ) + pow( -flocTel3D[0][iTel] + fXcore + fSmax * cos( tSaz ) * sin( tSzen ), 2 ) + pow( -flocTel3D[1][iTel] + fYcore + fSmax * sin( tSzen ) * sin( tSaz ), 2 ) - pow( ipZ * ( -flocTel3D[2][iTel] + fSmax * cos( tSzen ) ) + ipX * ( -flocTel3D[0][iTel] + fXcore +  fSmax * cos( tSaz ) * sin( tSzen ) ) + ipY * ( -flocTel3D[1][iTel] + fYcore +  fSmax * sin( tSzen ) * sin( tSaz ) ), 2 ) ) / pow( fsigmaT, 2 ) ) ) / ( 2.*fsigmaT * TMath::Pi() * sqrt( pow( fsigmaT, 2 ) * pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) + pow( fsigmaL, 2 ) * ( 1 - pow( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) );
    
    
    //// include dIedSel and dIedSaz ////
    
    double dIedSel = ( exp( ( eta - acos( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) ) / ( 4.*eta ) ) * eta * K * ( -( ( ipX * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) ) / ( pow( acos( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ), 2 ) * sqrt( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) + ( exp( ( eta - acos( ipZ * cos( tSzen ) +  ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) ) / ( 4.*eta ) ) * K * ( -( ( ipX * TMath::Pi() * cos( tSzen ) * cos( tSaz ) ) / 180. ) + ( ipZ * TMath::Pi() * sin( tSzen ) ) / 180. - ( ipY * TMath::Pi() * cos( tSzen ) * sin( tSaz ) ) / 180. ) ) / ( 4.*acos( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) * sqrt( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) );
    
    double dIedSaz = ( exp( ( eta - acos( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) ) / ( 4.*eta ) ) * eta * K * ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) ) / ( pow( acos( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ), 2 ) * sqrt( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) ) + ( exp( ( eta - acos( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ) ) ) / ( 4.*eta ) ) * K * ( ( ipY * TMath::Pi() * cos( tSaz ) * sin( tSzen ) ) / 180. - ( ipX * TMath::Pi() * sin( tSzen ) * sin( tSaz ) ) / 180. ) ) / ( 4.*acos( ipZ * cos( tSzen ) +   ipX * cos( tSaz ) * sin( tSzen ) +  ipY * sin( tSzen ) * sin( tSaz ) ) *  sqrt( 1 - pow( ipZ * cos( tSzen ) + ipX * cos( tSaz ) * sin( tSzen ) + ipY * sin( tSzen ) * sin( tSaz ), 2 ) ) );
    
    
    ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    
    if( epsilon <= eta )
    {
        grad[0] = Ie * dEq5dSel;
        grad[1] = Ie * dEq5dSaz;
    }
    else if( epsilon > eta )
    {
        //grad[0] = Ie * dEq5dSel;
        //grad[1] = Ie * dEq5dSaz;
        /// product rule for dSel and dSaz ///
        grad[0] = ( Ie * dEq5dSel ) + ( Eq5 * dIedSel );
        grad[1] = ( Ie * dEq5dSaz ) + ( Eq5 * dIedSaz );
    }
    grad[2] = Ie * dEq5dXcore;
    grad[3] = Ie * dEq5dYcore;
    grad[4] = Ie * dEq5dSmax;
    grad[5] = Ie * dEq5dsigmaL;
    grad[6] = Ie * dEq5dsigmaT;
    grad[7] = Ie * dEq5dNc;
    
    grad[4] *= 1000.; // rescale fSmax
    grad[5] *= 1000.; // rescale fsigmaL
    
    ////////// check for nan ////////////////
    // if(!isfinite(val)) {
    //   val = 0;
    // }
    // for (unsigned int iPar = 0; iPar < grad.size(); iPar++) grad[iPar] = 0;
    
}

//// standard normal cumulative distribution function (CDF) ////

double VModel3DFn::freq( double x ) const
{
    // Handbook of Mathematical Functions by Abramowitz and Stegun //
    
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    // Save the sign of x
    int sign = 1;
    if( x < 0 )
    {
        sign = -1;
    }
    
    x = fabs( x ) / sqrt( 2.0 );
    
    // A&S formula 7.1.26
    double t = 1.0 / ( 1.0 + p * x );
    double y = 1.0 - ( ( ( ( ( a5 * t + a4 ) * t ) + a3 ) * t + a2 ) * t + a1 ) * t * exp( -x * x );
    
    return 0.5 * ( 1.0 + sign * y );
}
