//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May  2 10:23:50 2007 by ROOT version 5.13/06
// from TTree EffFit/fit results for effective areas
// found on file: EffectiveAreas/effectiveArea_w0.5_ID08_ana12Fit.root
//////////////////////////////////////////////////////////

#ifndef cEffFit_h
#define cEffFit_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class cEffFit
{
    public :
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leave types
        Double_t        Ze;
        Int_t           AMC;
        TF1*             fEff;
        TGraphAsymmErrors* gEffAreaLog;
        TGraphAsymmErrors* gEffArea;
        Double_t        Fitxmin;
        Double_t        Fitxmax;
        
        // List of branches
        TBranch*        b_Ze;                     //!
        TBranch*        b_AMC;                    //!
        TBranch*        b_fEff;                   //!
        TBranch*        b_gEffAreaLog;            //!
        TBranch*        b_gEffArea;               //!
        TBranch*        b_Fitxmin;                //!
        TBranch*        b_Fitxmax;                //!
        
        cEffFit( TTree* tree = 0 );
        virtual ~cEffFit();
        virtual Int_t    Cut( Long64_t entry );
        virtual Int_t    GetEntry( Long64_t entry );
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
};
#endif

#ifdef cEffFit_cxx
cEffFit::cEffFit( TTree* tree )
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if( tree == 0 )
    {
        TFile* f = ( TFile* )gROOT->GetListOfFiles()->FindObject( "EffectiveAreas/effectiveArea_w0.5_ID08_ana12Fit.root" );
        if( !f )
        {
            f = new TFile( "EffectiveAreas/effectiveArea_w0.5_ID08_ana12Fit.root" );
        }
        tree = ( TTree* )gDirectory->Get( "EffFit" );
        
    }
    Init( tree );
}


cEffFit::~cEffFit()
{
    if( !fChain )
    {
        return;
    }
    delete fChain->GetCurrentFile();
}


Int_t cEffFit::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    return fChain->GetEntry( entry );
}


Long64_t cEffFit::LoadTree( Long64_t entry )
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
    if( !fChain->InheritsFrom( TChain::Class() ) )
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


void cEffFit::Init( TTree* tree )
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normaly not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set object pointer
    fEff = 0;
    gEffAreaLog = 0;
    gEffArea = 0;
    // Set branch addresses and branch pointers
    if( !tree )
    {
        return;
    }
    fChain = tree;
    fCurrent = -1;
    //	fChain->SetMakeClass( 1 );
    
    fChain->SetBranchAddress( "Ze", &Ze, &b_Ze );
    fChain->SetBranchAddress( "AMC", &AMC, &b_AMC );
    fChain->SetBranchAddress( "fEff", &fEff, &b_fEff );
    fChain->SetBranchAddress( "gEffAreaLog", &gEffAreaLog, &b_gEffAreaLog );
    fChain->SetBranchAddress( "gEffArea", &gEffArea, &b_gEffArea );
    fChain->SetBranchAddress( "Fitxmin", &Fitxmin, &b_Fitxmin );
    fChain->SetBranchAddress( "Fitxmax", &Fitxmax, &b_Fitxmax );
    Notify();
}


Bool_t cEffFit::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normaly not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}


void cEffFit::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}


Int_t cEffFit::Cut( Long64_t entry )
{
    entry = 1;
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif                                            // #ifdef cEffFit_cxx
