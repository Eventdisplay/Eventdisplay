//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 28 10:52:52 2007 by ROOT version 5.04/00
// from TTree tRunSummary/anasum results
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef CRunSummary_h
#define CRunSummary_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class CRunSummary : public TObject
{
    public :
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leave types
        Int_t           runOn;
        Int_t           runOff;
        Double_t        MJDOn;
        Double_t        MJDOn_runStart;
        Double_t        MJDOn_runStopp;
        Double_t        RunDurationOn;
        Double_t        MJDOff;
        Double_t        MJDOff_runStart;
        Double_t        MJDOff_runStopp;
        Double_t        RunDurationOff;
        Double_t        TargetRA;
        Double_t        TargetDec;
        Double_t        TargetRAJ2000;
        Double_t        TargetDecJ2000;
        Double_t        SkyMapCentreRAJ2000;
        Double_t        SkyMapCentreDecJ2000;
        Double_t        TargetShiftRAJ2000;
        Double_t        TargetShiftDecJ2000;
        Double_t        TargetShiftWest;
        Double_t        TargetShiftNorth;
        Double_t        WobbleNorth;
        Double_t        WobbleWest;
        UInt_t          NTel;
   Char_t          TelList[300];
        Double_t        tOn;
        Double_t        tOff;
        Double_t        elevationOn;
        Double_t        azimuthOn;
        Double_t        elevationOff;
        Double_t        azimuthOff;
   Double_t        Theta2Max;
        Double_t        RawRateOn;
        Double_t        RawRateOff;
        Double_t        pedvarsOn;
        Double_t        pedvarsOff;
        Double_t        NOn;
        Double_t        NOff;
        Double_t        NOffNorm;
        Double_t        OffNorm;
        Double_t        Signi;
        Double_t        Rate;
        Double_t        RateE;
        Double_t        RateOff;
        Double_t        RateOffE;
        Double_t        DeadTimeFracOn;
        Double_t        DeadTimeFracOff;
        Double_t        MaxSigni;
        Double_t        MaxSigniX;
        Double_t        MaxSigniY;
        
        // List of branches
        TBranch*        b_runOn;                  //!
        TBranch*        b_runOff;                 //!
        TBranch*        b_MJDOn;                  //!
        TBranch*        b_MJDOn_runStart;         //!
        TBranch*        b_MJDOn_runStopp;         //!
        TBranch*        b_RunDurationOn;          //!
        TBranch*        b_MJDOff;                 //!
        TBranch*        b_MJDOff_runStart;         //!
        TBranch*        b_MJDOff_runStopp;         //!
        TBranch*        b_RunDurationOff;          //!
        TBranch*        b_TargetRA;               //!
        TBranch*        b_TargetDec;              //!
        TBranch*        b_TargetRAJ2000;          //!
        TBranch*        b_TargetDecJ2000;         //!
        TBranch*        b_SkyMapCentreRAJ2000;    //!
        TBranch*        b_SkyMapCentreDecJ2000;   //!
        TBranch*        b_TargetShiftRAJ2000;     //!
        TBranch*        b_TargetShiftDecJ2000;    //!
        TBranch*        b_TargetShiftWest;        //!
        TBranch*        b_TargetShiftNorth;       //!
        TBranch*        b_WobbleNorth;            //!
        TBranch*        b_WobbleWest;             //!
        TBranch*        b_NTel;                   //!
   TBranch        *b_TelList;   //!
        TBranch*        b_tOn;                    //!
        TBranch*        b_tOff;                   //!
        TBranch*        b_elevationOn;            //!
        TBranch*        b_azimuthOn;              //!
        TBranch*        b_elevationOff;           //!
        TBranch*        b_azimuthOff;             //!
   TBranch        *b_Theta2Max;   //!
        TBranch*        b_RawRateOn;              //!
        TBranch*        b_RawRateOff;             //!
        TBranch*        b_pedvarsOn;              //!
        TBranch*        b_pedvarsOff;             //!
        TBranch*        b_NOn;                    //!
        TBranch*        b_NOff;                   //!
        TBranch*        b_NOffNorm;               //!
        TBranch*        b_OffNorm;                //!
        TBranch*        b_Signi;                  //!
        TBranch*        b_Rate;                   //!
        TBranch*        b_RateE;                  //!
        TBranch*        b_RateOff;                //!
        TBranch*        b_RateOffE;               //!
        TBranch*        b_DeadTimeFracOn;         //!
        TBranch*        b_DeadTimeFracOff;        //!
        TBranch*        b_MaxSigni;               //!
        TBranch*        b_MaxSigniX;              //!
        TBranch*        b_MaxSigniY;              //!
        
        CRunSummary( TTree* tree = 0 );
        virtual ~CRunSummary();
//        virtual Int_t    Cut( Long64_t entry );
        virtual Int_t    GetEntry( Long64_t entry );
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
        
        ClassDef( CRunSummary, 1 );
};
#endif

#ifdef CRunSummary_cxx

CRunSummary::CRunSummary( TTree* tree )
{
    fChain = 0;
    Init( tree );
}


CRunSummary::~CRunSummary()
{
    if( fChain && fChain->GetCurrentFile() )
    {
        delete fChain->GetCurrentFile();
    }
}


Int_t CRunSummary::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    Int_t a = fChain->GetEntry( entry );
    
    // ensure kind of backwards compatibility
    if( !b_MJDOn_runStart )
    {
        MJDOn_runStart = MJDOn - ( tOn / 2. ) / 86400.;
        MJDOn_runStopp = MJDOn + ( tOn / 2. ) / 86400.;
        RunDurationOn  = tOn;
        MJDOff_runStart = MJDOff - ( tOff / 2. ) / 86400.;
        MJDOff_runStopp = MJDOff + ( tOff / 2. ) / 86400.;
        RunDurationOff  = tOff;
    }
    
    return a;
}


Long64_t CRunSummary::LoadTree( Long64_t entry )
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


void CRunSummary::Init( TTree* tree )
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses of the tree
    // will be set. It is normaly not necessary to make changes to the
    // generated code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running with PROOF.
    
    // Set branch addresses
    if( tree == 0 )
    {
        return;
    }
    fChain = tree;
    fCurrent = -1;
    
    fChain->SetBranchAddress( "runOn", &runOn );
    fChain->SetBranchAddress( "runOff", &runOff );
    fChain->SetBranchAddress( "MJDOn", &MJDOn );
    fChain->SetBranchAddress( "MJDOn_runStart", &MJDOn_runStart );
    fChain->SetBranchAddress( "MJDOn_runStopp", &MJDOn_runStopp );
    fChain->SetBranchAddress( "RunDurationOn", &RunDurationOn );
    fChain->SetBranchAddress( "MJDOff", &MJDOff );
    fChain->SetBranchAddress( "MJDOff_runStart", &MJDOff_runStart );
    fChain->SetBranchAddress( "MJDOff_runStopp", &MJDOff_runStopp );
    fChain->SetBranchAddress( "RunDurationOff", &RunDurationOff );
    fChain->SetBranchAddress( "TargetRA", &TargetRA );
    fChain->SetBranchAddress( "TargetDec", &TargetDec );
    fChain->SetBranchAddress( "TargetRAJ2000", &TargetRAJ2000 );
    fChain->SetBranchAddress( "TargetDecJ2000", &TargetDecJ2000 );
    fChain->SetBranchAddress( "SkyMapCentreRAJ2000", &SkyMapCentreRAJ2000 );
    fChain->SetBranchAddress( "SkyMapCentreDecJ2000", &SkyMapCentreDecJ2000 );
    fChain->SetBranchAddress( "TargetShiftRAJ2000", &TargetShiftRAJ2000 );
    fChain->SetBranchAddress( "TargetShiftDecJ2000", &TargetShiftDecJ2000 );
    fChain->SetBranchAddress( "TargetShiftWest", &TargetShiftWest );
    fChain->SetBranchAddress( "TargetShiftNorth", &TargetShiftNorth );
    fChain->SetBranchAddress( "WobbleNorth", &WobbleNorth );
    fChain->SetBranchAddress( "WobbleWest", &WobbleWest );
    fChain->SetBranchAddress( "NTel", &NTel );
   fChain->SetBranchAddress("TelList", TelList, &b_TelList);
    fChain->SetBranchAddress( "tOn", &tOn );
    fChain->SetBranchAddress( "tOff", &tOff );
    fChain->SetBranchAddress( "elevationOn", &elevationOn );
    fChain->SetBranchAddress( "azimuthOn", &azimuthOn );
    fChain->SetBranchAddress( "elevationOff", &elevationOff );
    fChain->SetBranchAddress( "azimuthOff", &azimuthOff );
   fChain->SetBranchAddress("Theta2Max", &Theta2Max, &b_Theta2Max);
    fChain->SetBranchAddress( "RawRateOn", &RawRateOn );
    fChain->SetBranchAddress( "RawRateOff", &RawRateOff );
    if( fChain->GetBranchStatus( "pedvarsOn" ) )
    {
        fChain->SetBranchAddress( "pedvarsOn", &pedvarsOn );
        fChain->SetBranchAddress( "pedvarsOff", &pedvarsOff );
    }
    // no pedvars given, assume galactic source
    else
    {
        pedvarsOn = 8.1;
        pedvarsOff = 8.1;
    }
    fChain->SetBranchAddress( "NOn", &NOn );
    fChain->SetBranchAddress( "NOff", &NOff );
    fChain->SetBranchAddress( "NOffNorm", &NOffNorm );
    fChain->SetBranchAddress( "OffNorm", &OffNorm );
    fChain->SetBranchAddress( "Signi", &Signi );
    fChain->SetBranchAddress( "Rate", &Rate );
    fChain->SetBranchAddress( "RateE", &RateE );
    fChain->SetBranchAddress( "RateOff", &RateOff );
    fChain->SetBranchAddress( "RateOffE", &RateOffE );
    fChain->SetBranchAddress( "DeadTimeFracOn", &DeadTimeFracOn );
    fChain->SetBranchAddress( "DeadTimeFracOff", &DeadTimeFracOff );
    fChain->SetBranchAddress( "MaxSigni", &MaxSigni );
    fChain->SetBranchAddress( "MaxSigniX", &MaxSigniX );
    fChain->SetBranchAddress( "MaxSigniY", &MaxSigniY );
    Notify();
}


Bool_t CRunSummary::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. Typically here the branch pointers
    // will be retrieved. It is normaly not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed.
    
    // Get branch pointers
    b_runOn = fChain->GetBranch( "runOn" );
    b_runOff = fChain->GetBranch( "runOff" );
    b_MJDOn = fChain->GetBranch( "MJDOn" );
    b_MJDOn_runStart = fChain->GetBranch( "MJDOn_runStart" );
    b_MJDOn_runStopp = fChain->GetBranch( "MJDOn_runStopp" );
    b_RunDurationOn = fChain->GetBranch( "RunDurationOn" );
    b_MJDOff = fChain->GetBranch( "MJDOff" );
    b_MJDOff_runStart = fChain->GetBranch( "MJDOff_runStart" );
    b_MJDOff_runStopp = fChain->GetBranch( "MJDOff_runStopp" );
    b_RunDurationOff = fChain->GetBranch( "RunDurationOff" );
    b_TargetRA = fChain->GetBranch( "TargetRA" );
    b_TargetDec = fChain->GetBranch( "TargetDec" );
    b_TargetRAJ2000 = fChain->GetBranch( "TargetRAJ2000" );
    b_TargetDecJ2000 = fChain->GetBranch( "TargetDecJ2000" );
    b_SkyMapCentreRAJ2000 = fChain->GetBranch( "SkyMapCentreRAJ2000" );
    b_SkyMapCentreDecJ2000 = fChain->GetBranch( "SkyMapCentreDecJ2000" );
    b_TargetShiftRAJ2000 = fChain->GetBranch( "TargetShiftRAJ2000" );
    b_TargetShiftDecJ2000 = fChain->GetBranch( "TargetShiftDecJ2000" );
    b_TargetShiftWest = fChain->GetBranch( "TargetShiftWest" );
    b_TargetShiftNorth = fChain->GetBranch( "TargetShiftNorth" );
    b_WobbleNorth = fChain->GetBranch( "WobbleNorth" );
    b_WobbleWest = fChain->GetBranch( "WobbleWest" );
    b_NTel = fChain->GetBranch( "NTel" );
    b_tOn = fChain->GetBranch( "tOn" );
    b_tOff = fChain->GetBranch( "tOff" );
    if( fChain->GetBranchStatus( "pedvarsOn" ) )
    {
        b_pedvarsOn = fChain->GetBranch( "pedvarsOn" );
        b_pedvarsOff = fChain->GetBranch( "pedvarsOff" );
    }
    else
    {
        b_pedvarsOn = 0;
        b_pedvarsOff = 0;
    }
    b_elevationOn = fChain->GetBranch( "elevationOn" );
    b_azimuthOn = fChain->GetBranch( "azimuthOn" );
    b_elevationOff = fChain->GetBranch( "elevationOff" );
    b_azimuthOff = fChain->GetBranch( "azimuthOff" );
    b_RawRateOn = fChain->GetBranch( "RawRateOn" );
    b_RawRateOff = fChain->GetBranch( "RawRateOff" );
    b_NOn = fChain->GetBranch( "NOn" );
    b_NOff = fChain->GetBranch( "NOff" );
    b_NOffNorm = fChain->GetBranch( "NOffNorm" );
    b_OffNorm = fChain->GetBranch( "OffNorm" );
    b_Signi = fChain->GetBranch( "Signi" );
    b_Rate = fChain->GetBranch( "Rate" );
    b_RateE = fChain->GetBranch( "RateE" );
    b_RateOff = fChain->GetBranch( "RateOff" );
    b_RateOffE = fChain->GetBranch( "RateOffE" );
    b_DeadTimeFracOn = fChain->GetBranch( "DeadTimeFracOn" );
    b_DeadTimeFracOff = fChain->GetBranch( "DeadTimeFracOff" );
    b_MaxSigni = fChain->GetBranch( "MaxSigni" );
    b_MaxSigniX = fChain->GetBranch( "MaxSigniX" );
    b_MaxSigniY = fChain->GetBranch( "MaxSigniY" );
    
    return kTRUE;
}


void CRunSummary::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}

#endif                                            // #ifdef CRunSummary_cxx
