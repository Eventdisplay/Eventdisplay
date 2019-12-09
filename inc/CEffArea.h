//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 10 17:00:24 2008 by ROOT version 5.18/00
// from TTree fEffArea/effective area values
// found on file: effectiveArea.root
//////////////////////////////////////////////////////////

#ifndef CEffArea_h
#define CEffArea_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>

class CEffArea : public TObject
{
    public :
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leaf types
        Float_t        ze;
        Int_t           az;
        Float_t        azMin;
        Float_t        azMax;
        Float_t        Xoff;
        Float_t        Yoff;
        Float_t        Woff;
        Int_t           noise;
        Float_t        noisePE;
        Float_t        pedvar;
        Float_t        index;
        Int_t           nbins;
        Float_t        e0[1000];                 //[nbins]
        Float_t        eff[1000];                //[nbins]
        Float_t        seff_L[1000];             //[nbins]
        Float_t        seff_U[1000];             //[nbins]
        Float_t        eff_error[1000];                //[nbins]
        Float_t        esys_rel[1000];                //[nbins]
        Int_t           Rec_nbins;
        Float_t         Rec_e0[1000];             //[Rec_nbins]
        Float_t         Rec_eff[1000];            //[Rec_nbins]
        Float_t         Rec_seff_L[1000];         //[Rec_nbins]
        Float_t         Rec_seff_U[1000];         //[Rec_nbins]
        Float_t         Rec_eff_error[1000];            //[nbins]
        Float_t         Rec_angRes_p68[1000];     //[Rec_nbins]
        Float_t         Rec_angRes_p80[1000];     //[Rec_nbins]
        Float_t         Rec_angRes_kingSigma[1000]; //[Rec_nbins]
        Float_t         Rec_angRes_kingGamma[1000]; //[Rec_nbins]
        TH1D*            hEmc;
        TH1D*            hEcut;
        TH1D*            hEcutUW;
        TH1D*            hEcut500;
        TH1D*            hEcutRec;
        TH1D*            hEcutRecUW;
        TGraphAsymmErrors* gEffAreaMC;
        TGraphAsymmErrors* gEffAreaRec;
        TGraphAsymmErrors* gEffAreaNoTh2MC;
        TGraphAsymmErrors* gEffAreaNoTh2Rec;
        TProfile*        hEmcSWeight;
        TProfile*        hEsysRec;
        TProfile*        hEsysMC;
        TProfile*        hEsysMCRelative;
        TH2D*            hEsysMCRelativeRMS;
        TH2D*            hEsysMCRelative2D;
        TH2D*            hEsysMCRelative2DNoDirectionCut;
        TH2D*            hEsys2D;
        TH2D*            hResponseMatrix;
        TH2D*            hResponseMatrixFine;
        TH2D*            hResponseMatrixFineQC;
        TH2D*            hResponseMatrixQC;
        TH2D*            hResponseMatrixNoDirectionCuts;
        TH2D*            hResponseMatrixFineNoDirectionCuts;
        TH2D*            hAngularDiff_2D;
        TH2D*            hAngularDiffEmc_2D;
        TH2D*            hAngularLogDiff_2D;
        TH2D*            hAngularLogDiffEmc_2D;
        TH1D*            hhEcutTrigger;
        TH1D*            hhEcutFiducialArea;
        TH1D*            hhEcutStereoQuality;
        TH1D*            hhEcutTelType;
        TH1D*            hhEcutDirection;
        TH1D*            hhEcutGammaHadron;
        TH1D*            hhEcutEnergyReconstruction;
        TH1D*            hWeightedRate;
        TH1D*            hWeightedRate005;
        
        // List of branches
        TBranch*        b_ze;                     //!
        TBranch*        b_az;                     //!
        TBranch*        b_azMin;                  //!
        TBranch*        b_azMax;                  //!
        TBranch*        b_Xoff;                   //!
        TBranch*        b_Yoff;                   //!
        TBranch*        b_Woff;                   //!
        TBranch*        b_noise;                  //!
        TBranch*        b_noisePE;                //!
        TBranch*        b_pedvar;                 //!
        TBranch*        b_index;                  //!
        TBranch*        b_nbins;                  //!
        TBranch*        b_e0;                     //!
        TBranch*        b_eff;                    //!
        TBranch*        b_eff_error;                    //!
        TBranch*        b_esys_rel;                    //!
        TBranch*        b_seff_L;                 //!
        TBranch*        b_seff_U;                 //!
        TBranch*        b_Rec_nbins;              //!
        TBranch*        b_Rec_e0;                 //!
        TBranch*        b_Rec_angRes_p68;         //!
        TBranch*        b_Rec_angRes_p80;         //!
        TBranch*        b_Rec_angRes_kingSigma;   //!
        TBranch*        b_Rec_angRes_kingGamma;   //!
        TBranch*        b_Rec_eff;                //!
        TBranch*        b_Rec_eff_error;                //!
        TBranch*        b_Rec_seff_L;             //!
        TBranch*        b_Rec_seff_U;             //!
        TBranch*        b_hEmc;                   //!
        TBranch*        b_hEcut;                  //!
        TBranch*        b_hEcutUW;                  //!
        TBranch*        b_hEcut500;                  //!
        TBranch*        b_hEcutRec;               //!
        TBranch*        b_hEcutRecUW;               //!
        TBranch*        b_gEffAreaMC;             //!
        TBranch*        b_gEffAreaRec;            //!
        TBranch*        b_gEffAreaNoTh2MC;             //!
        TBranch*        b_gEffAreaNoTh2Rec;            //!
        TBranch*        b_hEmcSWeight;            //!
        TBranch*        b_hEsysRec;               //!
        TBranch*        b_hEsysMC;                //!
        TBranch*        b_hEsysMCRelative;        //!
        TBranch*        b_hEsysMCRelativeRMS;        //!
        TBranch*        b_hEsysMCRelative2D;        //!
        TBranch*        b_hEsysMCRelative2DNoDirectionCut;        //!
        TBranch*        b_hEsys2D;                //!
        TBranch*        b_hResponseMatrix;                //!
        TBranch*        b_hResponseMatrixFine;                //!
        TBranch*        b_hResponseMatrixQC;                //!
        TBranch*        b_hResponseMatrixFineQC;                //!
        TBranch*        b_hResponseMatrixNoDirectionCuts;                //!
        TBranch*        b_hResponseMatrixFineNoDirectionCuts;                //!
        TBranch*        b_hAngularDiff_2D;                //!
        TBranch*        b_hAngularDiff_Emc2D;                //!
        TBranch*        b_hAngularLogDiff_2D;                //!
        TBranch*        b_hAngularLogDiff_Emc2D;                //!
        TBranch*        b_hhEcutTrigger;   //!
        TBranch*        b_hhEcutFiducialArea;   //!
        TBranch*        b_hhEcutStereoQuality;   //!
        TBranch*        b_hhEcutTelType;   //!
        TBranch*        b_hhEcutDirection;   //!
        TBranch*        b_hhEcutGammaHadron;   //!
        TBranch*        b_hhEcutEnergyReconstruction;   //!
        TBranch*        b_hWeightedRate; //!
        TBranch*        b_hWeightedRate005; //!
        
        CEffArea( TTree* tree = 0 );
        virtual ~CEffArea();
        virtual Int_t    Cut( Long64_t entry );
        virtual Int_t    GetEntry( Long64_t entry );
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
        ClassDef( CEffArea, 4 );
};
#endif

#ifdef CEffArea_cxx

ClassImp( CEffArea )

CEffArea::CEffArea( TTree* tree )
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if( tree == 0 )
    {
        TFile* f = ( TFile* )gROOT->GetListOfFiles()->FindObject( "effectiveArea.root" );
        if( !f )
        {
            f = new TFile( "effectiveArea.root" );
        }
        tree = ( TTree* )gDirectory->Get( "fEffArea" );
        
    }
    Init( tree );
}


CEffArea::~CEffArea()
{
    if( !fChain )
    {
        return;
    }
    delete fChain->GetCurrentFile();
}


Int_t CEffArea::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    return fChain->GetEntry( entry );
}


Long64_t CEffArea::LoadTree( Long64_t entry )
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


void CEffArea::Init( TTree* tree )
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set object pointer
    hEmc = 0;
    hEcut = 0;
    hEcutUW = 0;
    hEcut500 = 0;
    hEcutRec = 0;
    hEcutRecUW = 0;
    gEffAreaMC = 0;
    gEffAreaRec = 0;
    gEffAreaNoTh2MC = 0;
    gEffAreaNoTh2Rec = 0;
    hEmcSWeight = 0;
    hEsysRec = 0;
    hEsysMC = 0;
    hEsysMCRelative = 0;
    hEsysMCRelativeRMS = 0;
    hEsysMCRelative2D = 0;
    hEsysMCRelative2DNoDirectionCut = 0;
    hEsys2D = 0;
    hResponseMatrix = 0;
    hResponseMatrixFine = 0;
    hResponseMatrixQC = 0;
    hResponseMatrixFineQC = 0;
    hResponseMatrixNoDirectionCuts = 0;
    hResponseMatrixFineNoDirectionCuts = 0;
    hAngularDiff_2D = 0;
    hAngularDiffEmc_2D = 0;
    hAngularLogDiff_2D = 0;
    hAngularLogDiffEmc_2D = 0;
    hhEcutTrigger = 0;
    hhEcutFiducialArea = 0;
    hhEcutStereoQuality = 0;
    hhEcutTelType = 0;
    hhEcutDirection = 0;
    hhEcutGammaHadron = 0;
    hhEcutEnergyReconstruction = 0;
    hWeightedRate = 0;
    hWeightedRate005 = 0;
    // Set branch addresses and branch pointers
    if( !tree )
    {
        return;
    }
    fChain = tree;
    fCurrent = -1;
    //	fChain->SetMakeClass( 1 );
    
    fChain->SetBranchAddress( "ze", &ze, &b_ze );
    fChain->SetBranchAddress( "az", &az, &b_az );
    if( fChain->GetBranchStatus( "azMin" ) )
    {
        fChain->SetBranchAddress( "azMin", &azMin, &b_azMin );
        fChain->SetBranchAddress( "azMax", &azMax, &b_azMax );
    }
    else
    {
        azMin = 0;
        azMax = 0;
    }
    if( fChain->GetBranchStatus( "Xoff" ) )
    {
        fChain->SetBranchAddress( "Xoff", &Xoff, &b_Xoff );
        fChain->SetBranchAddress( "Yoff", &Yoff, &b_Yoff );
    }
    else
    {
        Xoff  = 0.;
        Yoff  = 0.;
    }
    fChain->SetBranchAddress( "Woff", &Woff, &b_Woff );
    fChain->SetBranchAddress( "noise", &noise, &b_noise );
    if( fChain->GetBranchStatus( "noisePE" ) )
    {
        fChain->SetBranchAddress( "noisePE", &noisePE, &b_noisePE );
    }
    else
    {
        noisePE = 0.;
    }
    fChain->SetBranchAddress( "pedvar", &pedvar, &b_pedvar );
    fChain->SetBranchAddress( "index", &index, &b_index );
    fChain->SetBranchAddress( "nbins", &nbins, &b_nbins );
    fChain->SetBranchAddress( "e0", e0, &b_e0 );
    fChain->SetBranchAddress( "eff", eff, &b_eff );
    if( fChain->GetBranchStatus( "eff_error" ) )
    {
        fChain->SetBranchAddress( "eff_error", eff_error, &b_eff_error );
        fChain->SetBranchAddress( "esys_rel", esys_rel, &b_esys_rel );
    }
    else
    {
        for( int i = 0; i < 1000; i++ )
        {
            eff_error[i] = 0.;
            esys_rel[i] = 0.;
            Rec_eff_error[i] = 0;
        }
    }
    if( fChain->GetBranchStatus( "Rec_eff_error" ) )
    {
        fChain->SetBranchAddress( "Rec_eff_error", Rec_eff_error, &b_Rec_eff_error );
    }
    else
    {
        for( int i = 0; i < 1000; i++ )
        {
            Rec_eff_error[i] = 0;
        }
    }
    
    if( fChain->GetBranchStatus( "seff_L" ) )
    {
        fChain->SetBranchAddress( "seff_L", seff_L, &b_seff_L );
        fChain->SetBranchAddress( "seff_U", seff_U, &b_seff_U );
        fChain->SetBranchAddress( "Rec_seff_L", Rec_seff_L, &b_Rec_seff_L );
        fChain->SetBranchAddress( "Rec_seff_U", Rec_seff_U, &b_Rec_seff_U );
    }
    else
    {
        for( int i = 0; i < 1000; i++ )
        {
            seff_L[i] = 0.;
            seff_U[i] = 0.;
            Rec_seff_L[i] = 0.;
            Rec_seff_U[i] = 0.;
        }
    }
    fChain->SetBranchAddress( "Rec_nbins", &Rec_nbins, &b_Rec_nbins );
    fChain->SetBranchAddress( "Rec_e0", Rec_e0, &b_Rec_e0 );
    fChain->SetBranchAddress( "Rec_eff", Rec_eff, &b_Rec_eff );
    
    if( fChain->GetBranchStatus( "Rec_angRes_p68" ) )
    {
        fChain->SetBranchAddress( "Rec_angRes_p68", Rec_angRes_p68, &b_Rec_angRes_p68 ) ;
    }
    if( fChain->GetBranchStatus( "Rec_angRes_p80" ) )
    {
        fChain->SetBranchAddress( "Rec_angRes_p80", Rec_angRes_p80, &b_Rec_angRes_p80 ) ;
    }
    if( fChain->GetBranchStatus( "Rec_angRes_kingSigma" ) )
    {
        fChain->SetBranchAddress( "Rec_angRes_kingSigma", Rec_angRes_kingSigma, &b_Rec_angRes_kingSigma ) ;
    }
    if( fChain->GetBranchStatus( "Rec_angRes_kingGamma" ) )
    {
        fChain->SetBranchAddress( "Rec_angRes_kingGamma", Rec_angRes_kingGamma, &b_Rec_angRes_kingGamma ) ;
    }
    if( fChain->GetBranchStatus( "hResponseMatrixFine" ) )
    {
        fChain->SetBranchAddress( "hResponseMatrixFine", &hResponseMatrixFine, &b_hResponseMatrixFine );
    }
    else
    {
        hResponseMatrixFine = 0;
    }
    if( fChain->GetBranchStatus( "hEmc" ) )
    {
        fChain->SetBranchAddress( "hEmc", &hEmc, &b_hEmc );
    }
    else
    {
        hEmc = 0;
    }
    if( fChain->GetBranchStatus( "hEcut" ) )
    {
        fChain->SetBranchAddress( "hEcut", &hEcut, &b_hEcut );
        if( fChain->GetBranchStatus( "hEcutUW" ) )
        {
            fChain->SetBranchAddress( "hEcutUW", &hEcutUW, &b_hEcutUW );
        }
        else
        {
            hEcutUW = 0;
        }
        fChain->SetBranchAddress( "hEcut500", &hEcut500, &b_hEcut500 );
        fChain->SetBranchAddress( "hEcutRec", &hEcutRec, &b_hEcutRec );
        if( fChain->GetBranchStatus( "hEcutRecUW" ) )
        {
            fChain->SetBranchAddress( "hEcutRecUW", &hEcutRecUW, &b_hEcutRecUW );
        }
        else
        {
            hEcutRecUW = 0;
        }
        fChain->SetBranchAddress( "gEffAreaMC", &gEffAreaMC, &b_gEffAreaMC );
        
        fChain->SetBranchAddress( "hEmcSWeight", &hEmcSWeight, &b_hEmcSWeight );
        fChain->SetBranchAddress( "hEsysRec", &hEsysRec, &b_hEsysRec );
        fChain->SetBranchAddress( "hEsysMC", &hEsysMC, &b_hEsysMC );
        fChain->SetBranchAddress( "hEsys2D", &hEsys2D, &b_hEsys2D );
        fChain->SetBranchAddress( "hResponseMatrix", &hResponseMatrix, &b_hResponseMatrix );
        if( fChain->GetBranchStatus( "hEmcCutCTA" ) )
        {
            fChain->SetBranchAddress( "hEmcCutCTA", &hResponseMatrixFine, &b_hResponseMatrixFine );
        }
        else if( fChain->GetBranchStatus( "hResponseMatrixFine" ) )
        {
            fChain->SetBranchAddress( "hResponseMatrixFine", &hResponseMatrixFine, &b_hResponseMatrixFine );
        }
        else
        {
            hResponseMatrixFine = 0;
        }
        
        if( fChain->GetBranchStatus( "hResponseMatrixQC" ) )
        {
            fChain->SetBranchAddress( "hResponseMatrixQC", &hResponseMatrixQC, &b_hResponseMatrixQC );
        }
        else
        {
            hResponseMatrixQC = 0;
        }
        if( fChain->GetBranchStatus( "hResponseMatrixFineQC" ) )
        {
            fChain->SetBranchAddress( "hResponseMatrixFineQC", &hResponseMatrixFineQC, &b_hResponseMatrixFineQC );
        }
        else
        {
            hResponseMatrixFineQC = 0;
        }
        if( fChain->GetBranchStatus( "hResponseMatrixNoDirectionCut" ) )
        {
            fChain->SetBranchAddress( "hResponseMatrixNoDirectionCut", &hResponseMatrixNoDirectionCuts, &b_hResponseMatrixNoDirectionCuts );
        }
        else
        {
            hResponseMatrixNoDirectionCuts = 0;
        }
        if( fChain->GetBranchStatus( "hResponseMatrixFineNoDirectionCut" ) )
        {
            fChain->SetBranchAddress( "hResponseMatrixFineNoDirectionCut", &hResponseMatrixFineNoDirectionCuts, &b_hResponseMatrixFineNoDirectionCuts );
        }
        else
        {
            hResponseMatrixFineNoDirectionCuts = 0;
        }
    }
    else
    {
        hEmc = 0;
        hEcut = 0;
        hEcutUW = 0;
        hEcut500 = 0;
        hEcutRec = 0;
        hEcutRecUW = 0;
        gEffAreaMC = 0;
        gEffAreaRec = 0;
        gEffAreaNoTh2MC = 0;
        gEffAreaNoTh2Rec = 0;
        hEmcSWeight = 0;
        hEsysRec = 0;
        hEsysMC = 0;
        hEsysMCRelative = 0;
        hEsysMCRelativeRMS = 0;
        hEsysMCRelative2D = 0;
        hEsysMCRelative2DNoDirectionCut = 0;
        hEsys2D = 0;
        hResponseMatrix = 0;
        hResponseMatrixFine = 0;
        hResponseMatrixQC = 0;
        hResponseMatrixFineQC = 0;
        hResponseMatrixNoDirectionCuts = 0;
        hResponseMatrixFineNoDirectionCuts = 0;
    }
    if( fChain->GetBranchStatus( "hEsysMCRelative" ) )
    {
        fChain->SetBranchAddress( "hEsysMCRelative", &hEsysMCRelative, &b_hEsysMCRelative );
    }
    else
    {
        hEsysMCRelative = 0;
    }
    if( fChain->GetBranchStatus( "hEsysMCRelativeRMS" ) )
    {
        fChain->SetBranchAddress( "hEsysMCRelativeRMS", &hEsysMCRelativeRMS, &b_hEsysMCRelativeRMS );
    }
    else
    {
        hEsysMCRelativeRMS = 0;
    }
    if( fChain->GetBranchStatus( "hEsysMCRelative2D" ) )
    {
        fChain->SetBranchAddress( "hEsysMCRelative2D", &hEsysMCRelative2D, &b_hEsysMCRelative2D );
    }
    else
    {
        hEsysMCRelative2D = 0;
    }
    if( fChain->GetBranchStatus( "hEsysMCRelative2DNoDirectionCut" ) )
    {
        fChain->SetBranchAddress( "hEsysMCRelative2DNoDirectionCut", &hEsysMCRelative2DNoDirectionCut,
                                  &b_hEsysMCRelative2DNoDirectionCut );
    }
    else
    {
        hEsysMCRelative2DNoDirectionCut = 0;
    }
    if( fChain->GetBranchStatus( "hAngularDiff_2D" ) )
    {
        fChain->SetBranchAddress( "hAngularDiff_2D", &hAngularDiff_2D, &b_hAngularDiff_2D );
    }
    else
    {
        hAngularDiff_2D = 0;
    }
    if( fChain->GetBranchStatus( "hAngularDiffEmc_2D" ) )
    {
        fChain->SetBranchAddress( "hAngularDiffEmc_2D", &hAngularDiffEmc_2D, &b_hAngularDiff_Emc2D );
    }
    else
    {
        hAngularDiffEmc_2D = 0;
    }
    if( fChain->GetBranchStatus( "hAngularLogDiff_2D" ) )
    {
        fChain->SetBranchAddress( "hAngularLogDiff_2D", &hAngularLogDiff_2D, &b_hAngularLogDiff_2D );
    }
    else
    {
        hAngularLogDiff_2D = 0;
    }
    if( fChain->GetBranchStatus( "hAngularLogDiffEmc_2D" ) )
    {
        fChain->SetBranchAddress( "hAngularLogDiffEmc_2D", &hAngularLogDiffEmc_2D, &b_hAngularLogDiff_Emc2D );
    }
    else
    {
        hAngularLogDiffEmc_2D = 0;
    }
    if( fChain->GetBranchStatus( "gEffAreaRec" ) )
    {
        fChain->SetBranchAddress( "gEffAreaRec", &gEffAreaRec, &b_gEffAreaRec );
    }
    else
    {
        gEffAreaRec = 0;
    }
    if( fChain->GetBranchStatus( "gEffAreaNoTh2MC" ) )
    {
        fChain->SetBranchAddress( "gEffAreaNoTh2MC", &gEffAreaNoTh2MC, &b_gEffAreaNoTh2MC );
        fChain->SetBranchAddress( "gEffAreaNoTh2Rec", &gEffAreaNoTh2Rec, &b_gEffAreaNoTh2Rec );
    }
    else
    {
        gEffAreaNoTh2MC = 0;
        gEffAreaNoTh2Rec = 0;
    }
    if( fChain->GetBranchStatus( "hhEcutTrigger" ) )
    {
        fChain->SetBranchAddress( "hhEcutTrigger", &hhEcutTrigger, &b_hhEcutTrigger );
        fChain->SetBranchAddress( "hhEcutFiducialArea", &hhEcutFiducialArea, &b_hhEcutFiducialArea );
        fChain->SetBranchAddress( "hhEcutStereoQuality", &hhEcutStereoQuality, &b_hhEcutStereoQuality );
        fChain->SetBranchAddress( "hhEcutTelType", &hhEcutTelType, &b_hhEcutTelType );
        fChain->SetBranchAddress( "hhEcutDirection", &hhEcutDirection, &b_hhEcutDirection );
        fChain->SetBranchAddress( "hhEcutEnergyReconstruction", &hhEcutEnergyReconstruction, &b_hhEcutEnergyReconstruction );
        fChain->SetBranchAddress( "hhEcutGammaHadron", &hhEcutGammaHadron, &b_hhEcutGammaHadron );
    }
    else
    {
        hhEcutTrigger = 0;
        hhEcutFiducialArea = 0;
        hhEcutStereoQuality = 0;
        hhEcutTelType = 0;
        hhEcutDirection = 0;
        hhEcutGammaHadron = 0;
        hhEcutEnergyReconstruction = 0;
    }
    if( fChain->GetBranchStatus( "hWeightedRate" ) )
    {
        fChain->SetBranchAddress( "hWeightedRate", &hWeightedRate, &b_hWeightedRate );
    }
    else
    {
        hWeightedRate = 0;
    }
    if( fChain->GetBranchStatus( "hWeightedRate005" ) )
    {
        fChain->SetBranchAddress( "hWeightedRate005", &hWeightedRate005, &b_hWeightedRate005 );
    }
    else
    {
        hWeightedRate005 = 0;
    }
    
    Notify();
}


Bool_t CEffArea::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}


void CEffArea::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}


Int_t CEffArea::Cut( Long64_t entry )
{
    entry = 1;
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif                                            // #ifdef CEffArea_cxx
