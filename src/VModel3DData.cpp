/*! \class VModel3DData
     VModel3DData class for data in VFitModel3D
     based on Lemoine-Goumard et al. 2006
     adapted by J. Grube and G. Gyuk
*/

#include "VModel3DData.h"

VModel3DData::VModel3DData()
{
    fNdim3D = 3; // 3 dimensional space (x,y,z)
}

VModel3DData::~VModel3DData()
{
}

//// normalize the vector ////
void VModel3DData::norm3D( double& ax, double& ay, double& az )
{
    double mag = sqrt( ( ax * ax ) + ( ay * ay ) + ( az * az ) );
    if( mag > 0 )
    {
        ax = ax / mag;
        ay = ay / mag;
        az = az / mag;
    }
}

//// compute the cross product ////
void VModel3DData::cross3D( double ax, double ay, double az, double bx, double by, double bz, double& cx, double& cy, double& cz )
{
    cx = ( ay * bz ) - ( az * by );
    cy = ( az * bx ) - ( ax * bz );
    cz = ( ax * by ) - ( ay * bx );
}

//// compute the dot product ////
double VModel3DData::dot3D( double ax, double ay, double az, double bx, double by, double bz )
{
    return ( ax * bx ) + ( ay * by ) + ( az * bz );
}

//// initialize data: called once ////
void VModel3DData::initModel3D( unsigned int iNTel3D, vector<unsigned int>& iNpix3D )
{
    fGoodEvent3D = true;
    fConverged3D = false;
    
    fNTel3D = iNTel3D;
    fNpix3D = iNpix3D;
    
    // general coordinate frame conversion
    fxg3D.resize( fNdim3D );
    fyg3D.resize( fNdim3D );
    fzg3D.resize( fNdim3D );
    fxg3D[0] = 1.;
    fxg3D[1] = 0.;
    fxg3D[2] = 0.;
    fyg3D[0] = 0.;
    fyg3D[1] = 1.;
    fyg3D[2] = 0.;
    fzg3D[0] = 0.;
    fzg3D[1] = 0.;
    fzg3D[2] = 1.;
    
    fxsg3D.resize( fNdim3D, 0 );
    fysg3D.resize( fNdim3D, 0 );
    fzsg3D.resize( fNdim3D, 0 );
    
    // telescope configuration parameters
    fDCPE.resize( fNTel3D, 0 );
    
    fMarea3D.resize( fNTel3D, 0 );
    
    fomegapix3D.resize( fNTel3D );
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        fomegapix3D[i].resize( fNpix3D[i] );
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fomegapix3D[i][j] = 0;
        }
    }
    flocTel3D.resize( fNdim3D );
    for( unsigned int i = 0; i < fNdim3D; i++ )
    {
        flocTel3D[i].resize( fNTel3D );
        for( unsigned int j = 0; j < fNTel3D; j++ )
        {
            flocTel3D[i][j] = 0;
        }
    }
    
    
    // pointing for each telescope
    fTze3D.resize( fNTel3D, 0 );
    fTaz3D.resize( fNTel3D, 0 );
    fT3D.resize( fNdim3D, 0 );
    
    fpX3D.resize( fNTel3D );
    fpY3D.resize( fNTel3D );
    fpZ3D.resize( fNTel3D );
    fcosptheta3D.resize( fNTel3D );
    
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        fpX3D[i].resize( fNpix3D[i] );
        fpY3D[i].resize( fNpix3D[i] );
        fpZ3D[i].resize( fNpix3D[i] );
        fcosptheta3D[i].resize( fNpix3D[i] );
        
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fpX3D[i][j] = 0;
            fpY3D[i][j] = 0;
            fpZ3D[i][j] = 0;
            fcosptheta3D[i][j] = 0;
        }
    }
    
    // Model3D calculated signals
    fMu3D.resize( fNTel3D );
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        fMu3D[i].resize( fNpix3D[i] );
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fMu3D[i][j] = 0.0;
        }
    }
    fMuTel3D.resize( fNTel3D, 0.0 );
    
    // measured signal and pedvar
    fClean3D.resize( fNTel3D );
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        fClean3D[i].resize( fNpix3D[i] );
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fClean3D[i][j] = false;
        }
    }
    fMeasuredSum3D.resize( fNTel3D );
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        fMeasuredSum3D[i].resize( fNpix3D[i] );
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fMeasuredSum3D[i][j] = 0;
        }
    }
    fPedvar3D.resize( fNTel3D );
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        fPedvar3D[i].resize( fNpix3D[i] );
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fPedvar3D[i][j] = 0;
        }
    }
    
    // model start parameters
    fStartSel3D = 0;
    fStartSaz3D = 0;
    fStartXcore3D = 0;
    fStartYcore3D = 0;
    fStartSmax3D = 0;
    fStartsigmaL3D = 0;
    fStartsigmaT3D = 0;
    fStartNc3D = 0;
    
    // model fit parameters
    fSel3D = 0;
    fSaz3D = 0;
    fXcore3D = 0;
    fYcore3D = 0;
    fSmax3D = 0;
    fsigmaL3D = 0;
    fsigmaT3D = 0;
    fNc3D = 0;
    
    // slant depth and reduced width
    fDepth3D = 0;
    fRWidth3D = 0;
    fErrRWidth3D = 0;
    
    // model and hillas: direction difference
    fOmega3D = -9999;
    
    // error in model start parameters
    fErrorSel3D = 0;
    fErrorSaz3D = 0;
    fErrorXcore3D = 0;
    fErrorYcore3D = 0;
    fErrorSmax3D = 0;
    fErrorsigmaL3D = 0;
    fErrorsigmaT3D = 0;
    fErrorNc3D = 0;
    
    // shower direction
    fXoffModel3D = -9999;
    fYoffModel3D = -9999;
    fXoffDeRot3D = -9999;
    fYoffDeRot3D = -9999;
    
    // goodness of fit
    fNDF3D = 0;
    fNDFTel3D.resize( fNTel3D, 0 );
    fStartGOF3D = 0;
    fGOF3D = 0;
    
}

//// initialize data: called for every event ////

void VModel3DData::initEventModel3D()
{
    fGoodEvent3D = true;
    fConverged3D = false;
    
    // general coordinate frame conversion
    for( unsigned int i = 0; i < fNdim3D; i++ )
    {
        fxsg3D[i] = 0;
        fysg3D[i] = 0;
        fzsg3D[i] = 0;
    }
    
    // pointing for each telescope
    for( unsigned int j = 0; j < fNTel3D; j++ )
    {
        fTze3D[j] = 0;
        fTaz3D[j] = 0;
    }
    for( unsigned int i = 0; i < fNdim3D; i++ )
    {
        fT3D[i] = 0;
    }
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fcosptheta3D[i][j] = 0;
        }
    }
    
    // pixel pointing
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fpX3D[i][j] = 0;
            fpY3D[i][j] = 0;
            fpZ3D[i][j] = 0;
        }
    }
    
    // measured signal and pedvar
    for( unsigned int i = 0; i < fNTel3D; i++ )
    {
        for( unsigned int j = 0; j < fNpix3D[i]; j++ )
        {
            fClean3D[i][j] = false;
            fMu3D[i][j] = 0.0;
            fMeasuredSum3D[i][j] = 0;
            fPedvar3D[i][j] = 0;
        }
    }
    
    // model start parameters
    fStartSel3D = 0;
    fStartSaz3D = 0;
    fStartXcore3D = 0;
    fStartYcore3D = 0;
    fStartSmax3D = 0;
    fStartsigmaL3D = 0;
    fStartsigmaT3D = 0;
    fStartNc3D = 0;
    
    // model fit parameters
    fSel3D = 0;
    fSaz3D = 0;
    fXcore3D = 0;
    fYcore3D = 0;
    fSmax3D = 0;
    fsigmaL3D = 0;
    fsigmaT3D = 0;
    fNc3D = 0;
    
    // slant depth and reduced width
    fDepth3D = 0;
    fRWidth3D = 0;
    fErrRWidth3D = 0;
    
    // model and hillas: direction difference
    fOmega3D = -9999;
    
    // error in model start parameters
    fErrorSel3D = 0;
    fErrorSaz3D = 0;
    fErrorXcore3D = 0;
    fErrorYcore3D = 0;
    fErrorSmax3D = 0;
    fErrorsigmaL3D = 0;
    fErrorsigmaT3D = 0;
    fErrorNc3D = 0;
    
    // shower direction
    fXoffModel3D = -9999;
    fYoffModel3D = -9999;
    fXoffDeRot3D = -9999;
    fYoffDeRot3D = -9999;
    
    // goodness of fit
    fNDF3D = 0;
    fNDFTel3D.resize( fNTel3D, 0 );
    fStartGOF3D = 0;
    fGOF3D = 0;
    
}
