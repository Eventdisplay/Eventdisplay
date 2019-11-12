/*! \class VPETree
    \brief reader class for the pe trees


*/

#include "VPETree.h"

VPETree::VPETree( TTree* t )
{
    fTree = t;
    
    fReadFullTree = true;
}


bool VPETree::initialize( bool iF )
{
    fReadFullTree = iF;
    
    fEventNumber = 0;
    fPrimaryEnergy = 0.;
    fXcore = 0.;
    fYcore = 0.;
    fXcos = 0.;
    fYcos = 0.;
    fXsource = 0.;
    fYsource = 0.;
    v_f_ID = new std::vector< unsigned int >();
    v_f_time = new std::vector< float >();
    b_v_f_ID = 0;
    b_v_f_time = 0;
    
    if( !fTree )
    {
        return true;
    }
    
    if( fReadFullTree )
    {
        fTree->SetBranchAddress( "eventNumber", &fEventNumber );
        fTree->SetBranchAddress( "primaryEnergy", &fPrimaryEnergy );
        fTree->SetBranchAddress( "Xcore", &fXcore );
        fTree->SetBranchAddress( "Ycore", &fYcore );
        fTree->SetBranchAddress( "Xcos", &fXcos );
        fTree->SetBranchAddress( "Ycos", &fYcos );
        fTree->SetBranchAddress( "Xsource", &fXsource );
        fTree->SetBranchAddress( "Ysource", &fYsource );
    }
    fTree->SetBranchAddress( "pixelID", &v_f_ID, &b_v_f_ID );
    fTree->SetBranchAddress( "time", &v_f_time, &b_v_f_time );
    
    return true;
}


int VPETree::getEntry( int iEntry )
{
    if( !fTree )
    {
        return 0;
    }
    
    return fTree->GetEntry( iEntry );
}


int VPETree::getEntries()
{
    if( !fTree )
    {
        return 0;
    }
    
    return fTree->GetEntries();
}
