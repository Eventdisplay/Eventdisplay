/*! \class VMCParameters
    \brief MC data storage class

*/

#include "VMCParameters.h"

VMCParameters::VMCParameters( bool iDebug )
{
    fDebug = iDebug;
    
    fTree = 0;
    
    reset();
}


void VMCParameters::initTree()
{
    fTree = new TTree( "MCpars", "Monte Carlo Parameters" );
    fTree->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );
    fTree->SetAutoSave( 100000000 );              // autosave when 100 Mbytes written
    
    fTree->Branch( "runNumber", &runNumber, "runNumber/i" );
    fTree->Branch( "eventNumber",  &eventNumber,  "eventNumber/i" );
    fTree->Branch( "MCprim", &MCprimary, "MCprimary/s" );
    fTree->Branch( "MCe0", &MCenergy, "MCenergy/F" );
    fTree->Branch( "MCxcore", &MCxcore, "MCxcore/F" );
    fTree->Branch( "MCycore", &MCycore, "MCycore/F" );
    fTree->Branch( "MCze", &MCze, "MCze/F" );
    fTree->Branch( "MCaz", &MCaz, "MCaz/F" );
    fTree->Branch( "MCxoff", &MCTel_Xoff, "MCxoff/F" );
    fTree->Branch( "MCyoff", &MCTel_Yoff, "MCyoff/F" );
    fTree->Branch( "MCCorsikaShowerID", &MCCorsikaShowerID, "MCCorsikaShowerID/I" );
    fTree->Branch( "MCCorsikaRunID", &MCCorsikaRunID, "MCCorsikaRunID/I" );
    fTree->Branch( "MCFirstInteractionHeight", &MCFirstInteractionHeight, "MCFirstInteractionHeight/F" );
    fTree->Branch( "MCFirstInteractionDepth", &MCFirstInteractionDepth, "MCFirstInteractionDepth/F" );
    fTree->Branch( "ArrayTrigger", &ArrayTrigger, "ArrayTrigger/s" );
    
}


void VMCParameters::reset()
{
    runNumber = 0;
    eventNumber = 0;
    MCprimary = 0;
    MCenergy = 0.;
    MCxcore = 0.;
    MCycore = 0.;
    MCze = 0.;
    MCaz = 0.;
    MCTel_Xoff = 0.;
    MCTel_Yoff = 0.;
    MCCorsikaShowerID = 0;
    MCFirstInteractionHeight = 0;
    MCFirstInteractionDepth = 0;
    MCCorsikaRunID = 0;
    ArrayTrigger = 0;
}
