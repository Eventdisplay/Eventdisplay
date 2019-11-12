//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  3 15:41:50 2010 by ROOT version 5.22/00a
// from TTree telconfig/detector configuration
// found on file: 10020.root
//////////////////////////////////////////////////////////

#ifndef Ctelconfig_h
#define Ctelconfig_h

#include "VGlobalRunParameter.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

#include <vector>

class Ctelconfig
{
    public :
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leaf types
        UInt_t          NTel;
        Int_t           TelID;
        ULong64_t       TelType;
        UInt_t          TelID_hyperArray;
        Float_t         TelX;
        Float_t         TelY;
        Float_t         TelZ;
        Int_t           NMirrors;
        Float_t         MirrorArea;
        Float_t         FOV;
        Float_t         FocalLength;
        Float_t         CameraScaleFactor;
        Float_t         CameraCentreOffset;
        Float_t         CameraRotation;
        UInt_t          NPixel;
        UInt_t          NSamples;
        Float_t         Sample_time_slice;
        UInt_t          NGains;
        Float_t         HiLoScale;
        Int_t           HiLoThreshold;
        Float_t         HiLoOffset;
        Float_t         XTubeMM[VDST_MAXCHANNELS];//[NPixel]
        Float_t         YTubeMM[VDST_MAXCHANNELS];//[NPixel]
        Float_t         RTubeMM[VDST_MAXCHANNELS];//[NPixel]
        //[NPixel]
        Float_t         XTubeDeg[VDST_MAXCHANNELS];
        //[NPixel]
        Float_t         YTubeDeg[VDST_MAXCHANNELS];
        //[NPixel]
        Float_t         RTubeDeg[VDST_MAXCHANNELS];
        
        // List of branches
        TBranch*        b_NTel;                   //!
        TBranch*        b_TelID;                  //!
        TBranch*        b_TelType;                //!
        TBranch*        b_TelID_hyperArray;       //!
        TBranch*        b_TelX;                   //!
        TBranch*        b_TelY;                   //!
        TBranch*        b_TelZ;                   //!
        TBranch*        b_NMirrors;   //!
        TBranch*        b_MirrorArea;   //!
        TBranch*        b_FOV;   //!
        TBranch*        b_FocalLength;            //!
        TBranch*        b_CameraScaleFactor;      //!
        TBranch*        b_CameraCentreOffset;     //!
        TBranch*        b_CameraRotation;         //!
        TBranch*        b_NPixel;                 //!
        TBranch*        b_NSamples;               //!
        TBranch*        b_Sample_time_slice;      //!
        TBranch*        b_NGains;   //!
        TBranch*        b_HiLoScale;   //!
        TBranch*        b_HiLoThreshold;   //!
        TBranch*        b_HiLoOffset;   //!
        TBranch*        b_XTubeMM;                //!
        TBranch*        b_YTubeMM;                //!
        TBranch*        b_RTubeMM;                //!
        TBranch*        b_XTubeDeg;               //!
        TBranch*        b_YTubeDeg;               //!
        TBranch*        b_RTubeDeg;               //!
        
        Ctelconfig( TTree* tree = 0 );
        virtual ~Ctelconfig();
        virtual double   getArrayCentreX();
        virtual double   getArrayCentreY();
        virtual double   getArrayMaxSize();
        virtual Int_t    GetEntry( Long64_t entry );
        virtual unsigned int getNTel();
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual bool     IsZombie();
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
};
#endif

#ifdef Ctelconfig_cxx
Ctelconfig::Ctelconfig( TTree* tree )
{
    Init( tree );
}


Ctelconfig::~Ctelconfig()
{
    if( !fChain )
    {
        return;
    }
    delete fChain->GetCurrentFile();
}

bool Ctelconfig::IsZombie()
{
    if( fChain )
    {
        return false;
    }
    
    return true;
}


Int_t Ctelconfig::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    return fChain->GetEntry( entry );
}


Long64_t Ctelconfig::LoadTree( Long64_t entry )
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


void Ctelconfig::Init( TTree* tree )
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set branch addresses and branch pointers
    if( !tree )
    {
        return;
    }
    fChain = tree;
    fCurrent = -1;
    //	fChain->SetMakeClass( 1 );
    
    NTel = 0;
    fChain->SetBranchAddress( "NTel", &NTel, &b_NTel );
    if( fChain->GetBranchStatus( "TelID" ) )
    {
        fChain->SetBranchAddress( "TelID", &TelID, &b_TelID );
    }
    else
    {
        TelID = 0;
    }
    if( fChain->GetBranchStatus( "TelType" ) )
    {
        fChain->SetBranchAddress( "TelType", &TelType, &b_TelType );
    }
    else
    {
        TelType = 1;
    }
    if( fChain->GetBranchStatus( "TelID_hyperArray" ) )
    {
        fChain->SetBranchAddress( "TelID_hyperArray", &TelID_hyperArray, &b_TelID_hyperArray );
    }
    else
    {
        TelID_hyperArray = 0;
    }
    fChain->SetBranchAddress( "TelX", &TelX, &b_TelX );
    fChain->SetBranchAddress( "TelY", &TelY, &b_TelY );
    fChain->SetBranchAddress( "TelZ", &TelZ, &b_TelZ );
    if( fChain->GetBranchStatus( "NMirrors" ) )
    {
        fChain->SetBranchAddress( "NMirrors", &NMirrors, &b_NMirrors );
    }
    else
    {
        NMirrors = 0;
    }
    if( fChain->GetBranchStatus( "MirrorArea" ) )
    {
        fChain->SetBranchAddress( "MirrorArea", &MirrorArea, &b_MirrorArea );
    }
    else
    {
        MirrorArea = 0;
    }
    if( fChain->GetBranchStatus( "FOV" ) )
    {
        fChain->SetBranchAddress( "FOV", &FOV, &b_FOV );
    }
    else
    {
        FOV = 0;
    }
    fChain->SetBranchAddress( "FocalLength", &FocalLength, &b_FocalLength );
    fChain->SetBranchAddress( "CameraScaleFactor", &CameraScaleFactor, &b_CameraScaleFactor );
    fChain->SetBranchAddress( "CameraCentreOffset", &CameraCentreOffset, &b_CameraCentreOffset );
    fChain->SetBranchAddress( "CameraRotation", &CameraRotation, &b_CameraRotation );
    fChain->SetBranchAddress( "NPixel", &NPixel, &b_NPixel );
    fChain->SetBranchAddress( "NSamples", &NSamples, &b_NSamples );
    if( fChain->GetBranchStatus( "Sample_time_slice" ) )
    {
        fChain->SetBranchAddress( "Sample_time_slice", &Sample_time_slice, &b_Sample_time_slice );
    }
    else
    {
        Sample_time_slice = 0.;
    }
    if( fChain->GetBranchStatus( "NGains" ) )
    {
        fChain->SetBranchAddress( "NGains", &NGains, &b_NGains );
    }
    else
    {
        NGains = 0;
    }
    if( fChain->GetBranchStatus( "HiLoScale" ) )
    {
        fChain->SetBranchAddress( "HiLoScale", &HiLoScale, &b_HiLoScale );
    }
    else
    {
        HiLoScale = 0.;
    }
    if( fChain->GetBranchStatus( "HiLoThreshold" ) )
    {
        fChain->SetBranchAddress( "HiLoThreshold", &HiLoThreshold, &b_HiLoThreshold );
    }
    else
    {
        HiLoThreshold = 0;
    }
    if( fChain->GetBranchStatus( "HiLoOffset" ) )
    {
        fChain->SetBranchAddress( "HiLoOffset", &HiLoOffset, &b_HiLoOffset );
    }
    else
    {
        HiLoOffset = 0.;
    }
    fChain->SetBranchAddress( "XTubeMM", XTubeMM, &b_XTubeMM );
    fChain->SetBranchAddress( "YTubeMM", YTubeMM, &b_YTubeMM );
    fChain->SetBranchAddress( "RTubeMM", RTubeMM, &b_RTubeMM );
    fChain->SetBranchAddress( "XTubeDeg", XTubeDeg, &b_XTubeDeg );
    fChain->SetBranchAddress( "YTubeDeg", YTubeDeg, &b_YTubeDeg );
    fChain->SetBranchAddress( "RTubeDeg", RTubeDeg, &b_RTubeDeg );
    Notify();
}


Bool_t Ctelconfig::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}


void Ctelconfig::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}

/*
 * get centre of gravity of
 * telecope positions
 */
double Ctelconfig::getArrayCentreX()
{
    if( fChain )
    {
        double iX = 0.;
        double iN = 0.;
        for( Long64_t i = 0; i < fChain->GetEntries(); i++ )
        {
            fChain->GetEntry( i );
            
            iX += TelX;
            iN++;
        }
        if( iN > 0. )
        {
            return iX / iN;
        }
    }
    return -1.e99;
}

/*
 * get centre of gravity of
 * telecope positions
 */
double Ctelconfig::getArrayCentreY()
{
    if( fChain )
    {
        double iY = 0.;
        double iN = 0.;
        for( Long64_t i = 0; i < fChain->GetEntries(); i++ )
        {
            fChain->GetEntry( i );
            
            iY += TelY;
            iN++;
        }
        if( iN > 0. )
        {
            return iY / iN;
        }
    }
    return -1.e99;
}

/*
 * get maximum distance of a telescope
 * from the COG of the arra
 *
 */
double Ctelconfig::getArrayMaxSize()
{
    if( fChain )
    {
        double iXc = getArrayCentreX();
        double iYc = getArrayCentreY();
        double iMax = 0.;
        for( Long64_t i = 0; i < fChain->GetEntries(); i++ )
        {
            fChain->GetEntry( i );
            
            if( TMath::Sqrt( ( TelX - iXc ) * ( TelX - iXc ) + ( TelY - iYc ) * ( TelY - iYc ) ) > iMax )
            {
                iMax = TMath::Sqrt( ( TelX - iXc ) * ( TelX - iXc ) + ( TelY - iYc ) * ( TelY - iYc ) );
            }
        }
        
        return iMax;
    }
    
    return 0.;
}

unsigned int Ctelconfig::getNTel()
{
    if( NTel > 0 )
    {
        return NTel;
    }
    
    if( fChain && fChain->GetEntries() > 0 )
    {
        fChain->GetEntry( 0 );
        return NTel;
    }
    
    return 0;
}

#endif                                            // #ifdef Ctelconfig_cxx
