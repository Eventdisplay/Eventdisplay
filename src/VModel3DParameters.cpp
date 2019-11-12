/*! \class VModel3DParameters
    \brief  VModel3Parameters storage class
*/

#include "VModel3DParameters.h"

VModel3DParameters::VModel3DParameters()
{
    fTreeModel3D = 0;
    reset();
}

void VModel3DParameters::initTree( string iName, string iTitle )
{
    cout << "Model3D: creating tree " << iName.c_str() << "  " << iTitle.c_str() << endl;
    
    fTreeModel3D = new TTree( iName.c_str(), iTitle.c_str() );
    fTreeModel3D->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    fTreeModel3D->SetAutoSave( 1000 );
    /// branches ///
    fTreeModel3D->Branch( "eventNumber",  &eventNumber,  "eventNumber/I" );
    fTreeModel3D->Branch( "Sel3D", &fSel3D, "Sel3D/F" );
    fTreeModel3D->Branch( "Saz3D", &fSaz3D, "Saz3D/F" );
    fTreeModel3D->Branch( "Xcore3D", &fXcore3D, "Xcore3D/F" );
    fTreeModel3D->Branch( "Ycore3D", &fYcore3D, "Ycore3D/F" );
    fTreeModel3D->Branch( "Smax3D", &fSmax3D, "Smax3D/F" );
    fTreeModel3D->Branch( "sigmaL3D", &fsigmaL3D, "sigmaL3D/F" );
    fTreeModel3D->Branch( "sigmaT3D", &fsigmaT3D, "sigmaT3D/F" );
    fTreeModel3D->Branch( "Nc3D", &fNc3D, "Nc3D/F" );
    fTreeModel3D->Branch( "Xoff3D", &fXoffModel3D, "Xoff3D/F" );
    fTreeModel3D->Branch( "Yoff3D", &fYoffModel3D, "Yoff3D/F" );
    fTreeModel3D->Branch( "XoffDeRot3D", &fXoffDeRot3D, "XoffDeRot3D/F" );
    fTreeModel3D->Branch( "YoffDeRot3D", &fYoffDeRot3D, "YoffDeRot3D/F" );
    fTreeModel3D->Branch( "Goodness3D", &fGoodness3D, "Goodness3D/F" );
    fTreeModel3D->Branch( "Omega3D", &fOmega3D, "Omega3D/F" );
    fTreeModel3D->Branch( "Depth3D", &fDepth3D, "Depth3D/F" );
    fTreeModel3D->Branch( "RWidth3D", &fRWidth3D, "RWidth3D/F" );
    fTreeModel3D->Branch( "ErrRWidth3D", &fErrRWidth3D, "ErrRWidth3D/F" );
    fTreeModel3D->Branch( "Converged3D", &fConverged3D, "Converged3D/O" );
    /// debug below ///
    fTreeModel3D->Branch( "StartGoodness3D", &fStartGoodness3D, "StartGoodness3D/F" );
    fTreeModel3D->Branch( "StartSel3D",   &fStartSel3D,   "StartSel3D/F" );
    fTreeModel3D->Branch( "StartSaz3D",   &fStartSaz3D,   "StartSaz3D/F" );
    fTreeModel3D->Branch( "StartXcore3D", &fStartXcore3D, "StartXcore3D/F" );
    fTreeModel3D->Branch( "StartYcore3D", &fStartYcore3D, "StartYcore3D/F" );
    fTreeModel3D->Branch( "StartSmax3D",  &fStartSmax3D,  "StartSmax3D/F" );
    fTreeModel3D->Branch( "StartsigmaL3D", &fStartsigmaL3D, "StartsigmaL3D/F" );
    fTreeModel3D->Branch( "StartsigmaT3D", &fStartsigmaT3D, "StartsigmaT3D/F" );
    fTreeModel3D->Branch( "StartNc3D",    &fStartNc3D,    "StartNc3D/F" );
    fTreeModel3D->Branch( "ErrorSel3D",   &fErrorSel3D,   "ErrorSel3D/F" );
    fTreeModel3D->Branch( "ErrorSaz3D",   &fErrorSaz3D,   "ErrorSaz3D/F" );
    fTreeModel3D->Branch( "ErrorXcore3D", &fErrorXcore3D, "ErrorXcore3D/F" );
    fTreeModel3D->Branch( "ErrorYcore3D", &fErrorYcore3D, "ErrorYcore3D/F" );
    fTreeModel3D->Branch( "ErrorSmax3D",  &fErrorSmax3D,  "ErrorSmax3D/F" );
    fTreeModel3D->Branch( "ErrorsigmaL3D", &fErrorsigmaL3D, "ErrorsigmaL3D/F" );
    fTreeModel3D->Branch( "ErrorsigmaT3D", &fErrorsigmaT3D, "ErrorsigmaT3D/F" );
    fTreeModel3D->Branch( "ErrorNc3D",    &fErrorNc3D,    "ErrorNc3D/F" );
    
}

void VModel3DParameters::reset()
{
    eventNumber = 0;
    fSel3D = 0.;
    fSaz3D = 0.;
    fXcore3D = 0.;
    fYcore3D = 0.;
    fSmax3D = 0.;
    fsigmaL3D = 0.;
    fsigmaT3D = 0.;
    fNc3D = 0.;
    fXoffModel3D = 0.;
    fYoffModel3D = 0.;
    fXoffDeRot3D = 0.;
    fYoffDeRot3D = 0.;
    fGoodness3D = 0.;
    fOmega3D = 0.;
    fDepth3D = 0.;
    fRWidth3D = 0.;
    fErrRWidth3D = 0.;
    fConverged3D = false;
    
    fStartGoodness3D = 0.;
    fStartSel3D = 0.;
    fStartSaz3D = 0.;
    fStartXcore3D = 0.;
    fStartYcore3D = 0.;
    fStartSmax3D = 0.;
    fStartsigmaL3D = 0.;
    fStartsigmaT3D = 0.;
    fStartNc3D = 0.;
    fErrorSel3D = 0.;
    fErrorSaz3D = 0.;
    fErrorXcore3D = 0.;
    fErrorYcore3D = 0.;
    fErrorSmax3D = 0.;
    fErrorsigmaL3D = 0.;
    fErrorsigmaT3D = 0.;
    fErrorNc3D = 0.;
}
