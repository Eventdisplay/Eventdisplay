//////////////////////////////////////////////////////////
// from TTree model3Dpars
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Adjusted to mscw_energy
//   DO NOT OVERWRITE BY DOING SIMPLY model3Dpars->MakeClass !!!!
////////////////////////////////////////////

#ifndef Cmodel3Dpars_h
#define Cmodel3Dpars_h

#include "VGlobalRunParameter.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Cmodel3Dpars
{
    public :
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leave types
        Int_t           eventNumber;
        // Model3D (JG)
        bool            bModel3D;
        Float_t         Smax3D;
        Float_t         sigmaL3D;
        Float_t         sigmaT3D;
        Float_t         Nc3D;
        Float_t         Xcore3D;
        Float_t         Ycore3D;
        Float_t         Xoff3D;
        Float_t         Yoff3D;
        Float_t         XoffDeRot3D;
        Float_t         YoffDeRot3D;
        Float_t         Goodness3D;
        Float_t         Omega3D;
        Float_t         Depth3D;
        Float_t         RWidth3D;
        Float_t         ErrRWidth3D;
        Float_t         ErrorsigmaT3D;
        bool            Converged3D;
        
        // List of branches
        TBranch*        b_eventNumber;            //!
        // Model3D (JG)
        TBranch*        b_Smax3D;
        TBranch*        b_sigmaL3D;
        TBranch*        b_sigmaT3D;
        TBranch*        b_Nc3D;
        TBranch*        b_Xcore3D;
        TBranch*        b_Ycore3D;
        TBranch*        b_Xoff3D;
        TBranch*        b_Yoff3D;
        TBranch*        b_XoffDeRot3D;
        TBranch*        b_YoffDeRot3D;
        TBranch*        b_Goodness3D;
        TBranch*        b_Omega3D;
        TBranch*        b_Depth3D;
        TBranch*        b_RWidth3D;
        TBranch*        b_ErrRWidth3D;
        TBranch*        b_ErrorsigmaT3D;
        TBranch*        b_Converged3D;
        
        
        // JG added Model3D
        Cmodel3Dpars( TTree* tree = 0 );
        virtual ~Cmodel3Dpars();
        virtual Int_t    GetEntry( Long64_t entry );
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
};
#endif

#ifdef Cmodel3Dpars_cxx
//JG: added Model3D
Cmodel3Dpars::Cmodel3Dpars( TTree* tree )
{
    if( !tree )
    {
        return;
    }
    Init( tree );
}


Cmodel3Dpars::~Cmodel3Dpars()
{
    if( !fChain )
    {
        return;
    }
    delete fChain->GetCurrentFile();
}


Int_t Cmodel3Dpars::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    
    int a = fChain->GetEntry( entry );
    
    return a;
}


Long64_t Cmodel3Dpars::LoadTree( Long64_t entry )
{
    // Set the environment to read one entry
    if( !fChain )
    {
        return -5;
    }
    Long64_t centry = fChain->LoadTree( entry );
    if( centry < 0 )
    {
        return centry;
    }
    if( fChain->IsA() != TChain::Class() )
    {
        return centry;
    }
    TChain* chain = ( TChain* )fChain;
    if( chain->GetTreeNumber() != fCurrent )
    {
        fCurrent = chain->GetTreeNumber();
        Notify();
    }
    return centry;
}


void Cmodel3Dpars::Init( TTree* tree )
{

    // Set branch addresses
    if( tree == 0 )
    {
        return;
    }
    fChain = tree;
    fCurrent = -1;
    //	fChain->SetMakeClass( 1 );
    
    fChain->SetBranchAddress( "eventNumber", &eventNumber );
    // Model3D parameters (JG)
    fChain->SetBranchAddress( "Smax3D", &Smax3D );
    fChain->SetBranchAddress( "sigmaL3D", &sigmaL3D );
    fChain->SetBranchAddress( "sigmaT3D", &sigmaT3D );
    fChain->SetBranchAddress( "Nc3D", &Nc3D );
    fChain->SetBranchAddress( "Xcore3D", &Xcore3D );
    fChain->SetBranchAddress( "Ycore3D", &Ycore3D );
    fChain->SetBranchAddress( "Xoff3D", &Xoff3D );
    fChain->SetBranchAddress( "Yoff3D", &Yoff3D );
    fChain->SetBranchAddress( "XoffDeRot3D", &XoffDeRot3D );
    fChain->SetBranchAddress( "YoffDeRot3D", &YoffDeRot3D );
    fChain->SetBranchAddress( "Goodness3D", &Goodness3D );
    fChain->SetBranchAddress( "Omega3D", &Omega3D );
    fChain->SetBranchAddress( "Depth3D", &Depth3D );
    fChain->SetBranchAddress( "RWidth3D", &RWidth3D );
    fChain->SetBranchAddress( "ErrRWidth3D", &ErrRWidth3D );
    fChain->SetBranchAddress( "ErrorsigmaT3D", &ErrorsigmaT3D );
    fChain->SetBranchAddress( "Converged3D", &Converged3D );
    
    Notify();
}

Bool_t Cmodel3Dpars::Notify()
{
    b_eventNumber = fChain->GetBranch( "eventNumber" );
    // Model3D (JG)
    b_Smax3D = fChain->GetBranch( "Smax3D" );
    b_sigmaL3D = fChain->GetBranch( "sigmaL3D" );
    b_sigmaT3D = fChain->GetBranch( "sigmaT3D" );
    b_Nc3D = fChain->GetBranch( "Nc3D" );
    b_Xcore3D = fChain->GetBranch( "Xcore3D" );
    b_Ycore3D = fChain->GetBranch( "Ycore3D" );
    b_Xoff3D = fChain->GetBranch( "Xoff3D" );
    b_Yoff3D = fChain->GetBranch( "Yoff3D" );
    b_XoffDeRot3D = fChain->GetBranch( "XoffDeRot3D" );
    b_YoffDeRot3D = fChain->GetBranch( "YoffDeRot3D" );
    b_Goodness3D = fChain->GetBranch( "Goodness3D" );
    b_Omega3D  = fChain->GetBranch( "Omega3D" );
    b_Depth3D  = fChain->GetBranch( "Depth3D" );
    b_RWidth3D = fChain->GetBranch( "RWidth3D" );
    b_ErrRWidth3D = fChain->GetBranch( "ErrRWidth3D" );
    b_Converged3D = fChain->GetBranch( "Converged3D" );
    
    return kTRUE;
}


void Cmodel3Dpars::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}


#endif                                            // #ifdef Cmodel3Dpars_cxx
