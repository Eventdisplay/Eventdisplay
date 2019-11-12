//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  8 14:35:23 2007 by ROOT version 5.10/00
// from TTree showerpars/Shower Parameters
// found on file: output/32855.root
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Adjusted to mscw_energy
//
//   DO NOT OVERWRITE BY DOING SIMPLY A showerpars->MakeClass !!!!
//
//
//
//   (GM)
//
////////////////////////////////////////////

#ifndef Cshowerpars_h
#define Cshowerpars_h

#include "VGlobalRunParameter.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Cshowerpars
{
    public :
        bool            bMC;
        bool            bShort;
        bool            bUseArrayPointing;
        
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leave types
        Int_t           runNumber;
        Int_t           eventNumber;
        UInt_t          eventStatus;
        Int_t           MJD;
        Double_t        Time;
        Int_t           dataFormat;
        UInt_t          NTel;
        Int_t           traceFit;
        //[NTel]
        Float_t         TelElevation[VDST_MAXTELESCOPES];
        //[NTel]
        Float_t         TelAzimuth[VDST_MAXTELESCOPES];
        //[NTel]
        Float_t         TelElevationVBF[VDST_MAXTELESCOPES];
        //[NTel]
        Float_t         TelAzimuthVBF[VDST_MAXTELESCOPES];
        //[NTel]
        Float_t         TelPointingMismatch[VDST_MAXTELESCOPES];
        //[NTel]
        Float_t         TelDec[VDST_MAXTELESCOPES];
        Float_t         TelRA[VDST_MAXTELESCOPES];//[NTel]
        //[NTel]
        Float_t         PointingErrorX[VDST_MAXTELESCOPES];
        //[NTel]
        Float_t         PointingErrorY[VDST_MAXTELESCOPES];
        Float_t         Tel_x_SC[VDST_MAXTELESCOPES];
        Float_t         Tel_y_SC[VDST_MAXTELESCOPES];
        Float_t         Tel_z_SC[VDST_MAXTELESCOPES];
        Float_t         ArrayPointing_Elevation;
        Float_t         ArrayPointing_Azimuth;
        Float_t         WobbleN;
        Float_t         WobbleE;
        UInt_t          NTrig;
        ULong64_t       LTrig;
        //[NTrig]
        UShort_t        Trig_list[VDST_MAXTELESCOPES];
        UShort_t        Trig_type[VDST_MAXTELESCOPES];
        UInt_t          NMethods;
        //[NMethods]
        UShort_t        MethodID[VDST_MAXRECMETHODS];
        //[NMethods]
        UShort_t        NImages[VDST_MAXRECMETHODS];
        //[NMethods]
        ULong64_t       ImgSel[VDST_MAXRECMETHODS];
        //[NMethods]
        UChar_t         ImgSel_list[VDST_MAXRECMETHODS][VDST_MAXTELESCOPES];
        //[NMethods]
        Float_t         img2_ang[VDST_MAXRECMETHODS];
        Float_t         Ze[VDST_MAXRECMETHODS];   //[NMethods]
        Float_t         Az[VDST_MAXRECMETHODS];   //[NMethods]
        Float_t         Xoff[VDST_MAXRECMETHODS]; //[NMethods]
        Float_t         Yoff[VDST_MAXRECMETHODS]; //[NMethods]
        //[NMethods]
        Float_t         XoffDeRot[VDST_MAXRECMETHODS];
        //[NMethods]
        Float_t         YoffDeRot[VDST_MAXRECMETHODS];
        Float_t         stds[VDST_MAXRECMETHODS]; //[NMethods]
        Float_t         dec[VDST_MAXRECMETHODS];  //[NMethods]
        Float_t         ra[VDST_MAXRECMETHODS];   //[NMethods]
        Float_t         Xcore[VDST_MAXRECMETHODS];//[NMethods]
        Float_t         Ycore[VDST_MAXRECMETHODS];//[NMethods]
        //[NMethods]
        Float_t         Xcore_SC[VDST_MAXRECMETHODS];
        //[NMethods]
        Float_t         Ycore_SC[VDST_MAXRECMETHODS];
        Float_t         stdp[VDST_MAXRECMETHODS]; //[NMethods]
        Float_t         Chi2[VDST_MAXRECMETHODS]; //[NMethods]
        Float_t         DispDiff[VDST_MAXRECMETHODS];
        Int_t           MCprim;
        Float_t        MCe0;
        Float_t        MCxcore;
        Float_t        MCycore;
        Float_t        MCxcos;
        Float_t        MCycos;
        Float_t        MCze;
        Float_t        MCaz;
        Float_t        MCxoff;
        Float_t        MCyoff;
        Float_t        MCxcore_SC;
        Float_t        MCycore_SC;
        Float_t        MCzcore_SC;
        Int_t	       MCCorsikaRunID;
        Int_t	       MCCorsikaShowerID;
        Float_t        MCFirstInteractionHeight;
        Float_t	       MCFirstInteractionDepth;
        
        
        // List of branches
        TBranch*        b_runNumber;              //!
        TBranch*        b_eventNumber;            //!
        TBranch*        b_MJD;                    //!
        TBranch*        b_Time;                   //!
        TBranch*        b_dataFormat;             //!
        TBranch*        b_NTel;                   //!
        TBranch*        b_traceFit;               //!
        TBranch*        b_TelElevation;           //!
        TBranch*        b_TelAzimuth;             //!
        TBranch*        b_TelElevationVBF;        //!
        TBranch*        b_TelAzimuthVBF;          //!
        TBranch*        b_TelPointingMismatch;    //!
        TBranch*        b_TelDec;                 //!
        TBranch*        b_TelRA;                  //!
        TBranch*        b_Tel_x_SC;               //!
        TBranch*        b_Tel_y_SC;               //!
        TBranch*        b_Tel_z_SC;               //!
        TBranch*        b_ArrayPointing_Elevation; //!
        TBranch*        b_ArrayPointing_Azimuth; //!
        TBranch*        b_WobbleN;                //!
        TBranch*        b_WobbleE;                //!
        TBranch*        b_NTrig;                  //!
        TBranch*        b_LTrig;                  //!
        TBranch*        b_Trig_list;              //!
        TBranch*        b_Trig_type;              //!
        TBranch*        b_NMethods;               //!
        TBranch*        b_MethodID;               //!
        TBranch*        b_NImages;                //!
        TBranch*        b_ImgSel;                 //!
        TBranch*        b_ImgSel_list;            //!
        TBranch*        b_img2_ang;               //!
        TBranch*        b_Ze;                     //!
        TBranch*        b_Az;                     //!
        TBranch*        b_Xoff;                   //!
        TBranch*        b_Yoff;                   //!
        TBranch*        b_XoffDeRot;                   //!
        TBranch*        b_YoffDeRot;                   //!
        TBranch*        b_stds;                   //!
        TBranch*        b_dec;                    //!
        TBranch*        b_ra;                     //!
        TBranch*        b_Xcore;                  //!
        TBranch*        b_Ycore;                  //!
        TBranch*        b_Xcore_SC;               //!
        TBranch*        b_Ycore_SC;               //!
        TBranch*        b_stdp;                   //!
        TBranch*        b_Chi2;                   //!
        TBranch*        b_DispDiff;               //!
        TBranch*        b_MCprim;                 //!
        TBranch*        b_MCe0;                   //!
        TBranch*        b_MCxcore;                //!
        TBranch*        b_MCycore;                //!
        TBranch*        b_MCxcos;                 //!
        TBranch*        b_MCycos;                 //!
        TBranch*        b_MCze;                   //!
        TBranch*        b_MCaz;                   //!
        TBranch*        b_MCxoff;                 //!
        TBranch*        b_MCyoff;                 //!
        TBranch*        b_MCxcore_SC;             //!
        TBranch*        b_MCycore_SC;             //!
        TBranch*        b_MCzcore_SC;             //!
        TBranch*        b_MCCorsikaRunID;         //!
        TBranch*        b_MCCorsikaShowerID;      //!
        TBranch*        b_MCFirstInteractionHeight;    //!
        TBranch*        b_MCFirstInteractionDepth;     //!
        
        Cshowerpars( TTree* tree = 0, bool iMC = false, bool iShort = false );
        virtual ~Cshowerpars();
        virtual Int_t    GetEntry( Long64_t entry );
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
        bool             isMC()
        {
            return bMC;
        }
        bool             isShort()
        {
            return bShort;
        }
};
#endif

#ifdef Cshowerpars_cxx

Cshowerpars::Cshowerpars( TTree* tree, bool iMC, bool iShort )
{
    if( !tree )
    {
        return;
    }
    bUseArrayPointing = false;
    
    bMC = iMC;
    bShort = iShort;
    
    Init( tree );
}


Cshowerpars::~Cshowerpars()
{
    if( !fChain )
    {
        return;
    }
    delete fChain->GetCurrentFile();
}


Int_t Cshowerpars::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    
    int a = fChain->GetEntry( entry );
    
    if( bUseArrayPointing )
    {
        for( UInt_t i = 0; i < NTel; i++ )
        {
            TelElevation[i] = ArrayPointing_Elevation;
            TelAzimuth[i] = ArrayPointing_Azimuth;
        }
    }
    return a;
}


Long64_t Cshowerpars::LoadTree( Long64_t entry )
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


void Cshowerpars::Init( TTree* tree )
{

    // Set branch addresses
    if( tree == 0 )
    {
        return;
    }
    fChain = tree;
    fCurrent = -1;
    //	fChain->SetMakeClass( 1 );
    
    if( tree->GetBranchStatus( "MCe0" ) )
    {
        bMC = true;
    }
    
    fChain->SetBranchAddress( "runNumber", &runNumber );
    fChain->SetBranchAddress( "eventNumber", &eventNumber );
    fChain->SetBranchAddress( "MJD", &MJD );
    fChain->SetBranchAddress( "Time", &Time );
    if( fChain->GetBranchStatus( "dataFormat" ) )
    {
        fChain->SetBranchAddress( "dataFormat", &dataFormat );
    }
    else
    {
        dataFormat = 0;
    }
    fChain->SetBranchAddress( "NTel", &NTel );
    if( !bShort )
    {
        fChain->SetBranchAddress( "traceFit", &traceFit );
    }
    else
    {
        traceFit = 0;
    }
    if( fChain->GetBranchStatus( "TelElevation" ) )
    {
        fChain->SetBranchAddress( "TelElevation", TelElevation );
        fChain->SetBranchAddress( "TelAzimuth", TelAzimuth );
    }
    else
    {
        bUseArrayPointing = true;
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            TelElevation[i] = 0;
            TelAzimuth[i] = 0;
        }
    }
    
    
    if( fChain->GetBranchStatus( "TelDec" ) )
    {
        fChain->SetBranchAddress( "TelDec", TelDec );
        fChain->SetBranchAddress( "TelRA", TelRA );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            TelDec[i] = 0.;
            TelRA[i] = 0;
        }
    }
    
    if( !bMC && !bShort )
    {
        fChain->SetBranchAddress( "TelElevationVBF", TelElevationVBF );
        fChain->SetBranchAddress( "TelAzimuthVBF", TelAzimuthVBF );
        fChain->SetBranchAddress( "TelPointingMismatch", TelPointingMismatch );
        fChain->SetBranchAddress( "Tel_x_SC", Tel_x_SC );
        fChain->SetBranchAddress( "Tel_y_SC", Tel_y_SC );
        fChain->SetBranchAddress( "Tel_z_SC", Tel_z_SC );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            TelElevationVBF[i] = 0.;
            TelAzimuthVBF[i] = 0.;
            TelPointingMismatch[i] = 0.;
            Tel_x_SC[i] = 0.;
            Tel_y_SC[i] = 0.;
            Tel_z_SC[i] = 0.;
        }
    }
    if( fChain->GetBranchStatus( "ArrayPointingElevation" ) )
    {
        fChain->SetBranchAddress( "ArrayPointingElevation", &ArrayPointing_Elevation );
        fChain->SetBranchAddress( "ArrayPointingAzimuth", &ArrayPointing_Azimuth );
    }
    else
    {
        ArrayPointing_Elevation = 0.;
        ArrayPointing_Azimuth = 0.;
    }
    
    fChain->SetBranchAddress( "WobbleN", &WobbleN );
    fChain->SetBranchAddress( "WobbleE", &WobbleE );
    fChain->SetBranchAddress( "NTrig", &NTrig );
    fChain->SetBranchAddress( "LTrig", &LTrig );
    fChain->SetBranchAddress( "Trig_list", Trig_list, &b_Trig_list );
    if( fChain->GetBranchStatus( "Trig_type" ) )
    {
        fChain->SetBranchAddress( "Trig_type", Trig_type, &b_Trig_type );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            Trig_type[i] = 0;
        }
    }
    fChain->SetBranchAddress( "NMethods", &NMethods );
    if( !bShort )
    {
        fChain->SetBranchAddress( "MethodID", MethodID );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXRECMETHODS; i++ )
        {
            MethodID[i] = 0;
        }
    }
    
    fChain->SetBranchAddress( "NImages", NImages );
    fChain->SetBranchAddress( "img2_ang", img2_ang );
    fChain->SetBranchAddress( "ImgSel", ImgSel );
    fChain->SetBranchAddress( "ImgSel_list", ImgSel_list, &b_ImgSel_list );
    fChain->SetBranchAddress( "Ze", Ze );
    fChain->SetBranchAddress( "Az", Az );
    fChain->SetBranchAddress( "Xoff", Xoff );
    fChain->SetBranchAddress( "Yoff", Yoff );
    if( tree->GetBranchStatus( "XoffDeRot" ) )
    {
        fChain->SetBranchAddress( "XoffDeRot", XoffDeRot );
        fChain->SetBranchAddress( "YoffDeRot", YoffDeRot );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXRECMETHODS; i++ )
        {
            XoffDeRot[i] = 0.;
            YoffDeRot[i] = 0.;
        }
    }
    if( !bShort )
    {
        fChain->SetBranchAddress( "stds", stds );
        if( fChain->GetBranchStatus( "dec" ) )
        {
            fChain->SetBranchAddress( "dec", dec );
            fChain->SetBranchAddress( "ra", ra );
        }
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXRECMETHODS; i++ )
        {
            stds[i] = 0.;
            dec[i] = 0.;
            ra[i] = 0.;
        }
    }
    fChain->SetBranchAddress( "Xcore", Xcore );
    fChain->SetBranchAddress( "Ycore", Ycore );
    if( !bShort )
    {
        fChain->SetBranchAddress( "Xcore_SC", Xcore_SC );
        fChain->SetBranchAddress( "Ycore_SC", Ycore_SC );
        fChain->SetBranchAddress( "stdp", stdp );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXRECMETHODS; i++ )
        {
            Xcore_SC[i] = 0.;
            Ycore_SC[i] = 0.;
            stdp[i] = 0.;
        }
    }
    fChain->SetBranchAddress( "Chi2", Chi2 );
    if( fChain->GetBranchStatus( "DispDiff" ) )
    {
        fChain->SetBranchAddress( "DispDiff", DispDiff );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXRECMETHODS; i++ )
        {
            DispDiff[i] = 0.;
        }
    }
    if( bMC )
    {
        if( fChain->GetBranchStatus( "MCprim" ) )
        {
            fChain->SetBranchAddress( "MCprim", &MCprim );
        }
        else
        {
            MCprim = 0;
        }
        fChain->SetBranchAddress( "MCe0", &MCe0 );
        fChain->SetBranchAddress( "MCxcore", &MCxcore );
        fChain->SetBranchAddress( "MCycore", &MCycore );
        if( !bShort )
        {
            fChain->SetBranchAddress( "MCxcos", &MCxcos );
            fChain->SetBranchAddress( "MCycos", &MCycos );
        }
        else
        {
            MCxcos = 0.;
            MCycos = 0.;
        }
        fChain->SetBranchAddress( "MCze", &MCze );
        fChain->SetBranchAddress( "MCaz", &MCaz );
        fChain->SetBranchAddress( "MCxoff", &MCxoff );
        fChain->SetBranchAddress( "MCyoff", &MCyoff );
        if( !bShort )
        {
            fChain->SetBranchAddress( "MCxcore_SC", &MCxcore_SC );
            fChain->SetBranchAddress( "MCycore_SC", &MCycore_SC );
            fChain->SetBranchAddress( "MCzcore_SC", &MCzcore_SC );
        }
        else
        {
            MCxcore_SC = 0.;
            MCycore_SC = 0.;
            MCzcore_SC = 0.;
        }
        if( fChain->GetBranch( "MCCorsikaRunID" ) )
        {
            fChain->SetBranchAddress( "MCCorsikaRunID", &MCCorsikaRunID );
            fChain->SetBranchAddress( "MCCorsikaShowerID", &MCCorsikaShowerID );
            fChain->SetBranchAddress( "MCFirstInteractionHeight", &MCFirstInteractionHeight );
            fChain->SetBranchAddress( "MCFirstInteractionDepth", &MCFirstInteractionDepth );
        }
        else
        {
            MCCorsikaRunID = -1;
            MCCorsikaShowerID = -1;
            MCFirstInteractionDepth = -1;
            MCFirstInteractionHeight = -1;
            
        }
    }
    Notify();
}


Bool_t Cshowerpars::Notify()
{
    b_runNumber = fChain->GetBranch( "runNumber" );
    b_eventNumber = fChain->GetBranch( "eventNumber" );
    b_MJD = fChain->GetBranch( "MJD" );
    b_Time = fChain->GetBranch( "Time" );
    if( !bShort )
    {
        b_dataFormat = fChain->GetBranch( "dataFormat" );
    }
    else
    {
        b_dataFormat = 0;
    }
    b_NTel = fChain->GetBranch( "NTel" );
    b_traceFit = fChain->GetBranch( "traceFit" );
    b_TelElevation = fChain->GetBranch( "TelElevation" );
    b_TelAzimuth = fChain->GetBranch( "TelAzimuth" );
    if( !bMC )
    {
        b_TelDec = fChain->GetBranch( "TelDec" );
        b_TelRA = fChain->GetBranch( "TelRA" );
    }
    else
    {
        b_TelDec = 0;
        b_TelRA = 0;
    }
    if( !bShort )
    {
        b_TelElevationVBF = fChain->GetBranch( "TelElevationVBF" );
        b_TelAzimuthVBF = fChain->GetBranch( "TelAzimuthVBF" );
        b_TelPointingMismatch = fChain->GetBranch( "TelPointingMismatch" );
        b_Tel_x_SC = fChain->GetBranch( "Tel_x_SC" );
        b_Tel_y_SC = fChain->GetBranch( "Tel_y_SC" );
        b_Tel_z_SC = fChain->GetBranch( "Tel_z_SC" );
    }
    else
    {
        b_TelElevationVBF = 0;
        b_TelAzimuthVBF = 0;
        b_TelPointingMismatch = 0;
        b_Tel_x_SC = 0;
        b_Tel_y_SC = 0;
        b_Tel_z_SC = 0;
    }
    b_WobbleN = fChain->GetBranch( "WobbleN" );
    b_WobbleE = fChain->GetBranch( "WobbleE" );
    b_ArrayPointing_Elevation = fChain->GetBranch( "ArrayPointing_Elevation" );
    b_ArrayPointing_Azimuth = fChain->GetBranch( "ArrayPointing_Azimuth" );
    b_NTrig = fChain->GetBranch( "NTrig" );
    b_LTrig = fChain->GetBranch( "LTrig" );
    b_NMethods = fChain->GetBranch( "NMethods" );
    if( !bShort )
    {
        b_MethodID = fChain->GetBranch( "MethodID" );
    }
    else
    {
        b_MethodID = 0;
    }
    b_NImages = fChain->GetBranch( "NImages" );
    b_ImgSel = fChain->GetBranch( "ImgSel" );
    b_ImgSel_list = fChain->GetBranch( "ImgSel_list" );
    if( !bShort )
    {
        b_img2_ang = fChain->GetBranch( "img2_ang" );
    }
    else
    {
        b_img2_ang = 0;
    }
    b_Ze = fChain->GetBranch( "Ze" );
    b_Az = fChain->GetBranch( "Az" );
    b_Xoff = fChain->GetBranch( "Xoff" );
    b_Yoff = fChain->GetBranch( "Yoff" );
    b_XoffDeRot = fChain->GetBranch( "XoffDeRot" );
    b_YoffDeRot = fChain->GetBranch( "YoffDeRot" );
    b_stds = fChain->GetBranch( "stds" );
    b_dec = fChain->GetBranch( "dec" );
    b_ra = fChain->GetBranch( "ra" );
    b_Xcore = fChain->GetBranch( "Xcore" );
    b_Ycore = fChain->GetBranch( "Ycore" );
    b_Xcore_SC = fChain->GetBranch( "Xcore_SC" );
    b_Ycore_SC = fChain->GetBranch( "Ycore_SC" );
    b_stdp = fChain->GetBranch( "stdp" );
    b_Chi2 = fChain->GetBranch( "Chi2" );
    b_DispDiff = fChain->GetBranch( "DispDiff" );
    if( bMC )
    {
        if( !bShort )
        {
            b_MCprim = fChain->GetBranch( "MCprim" );
        }
        else
        {
            b_MCprim = 0;
        }
        b_MCe0 = fChain->GetBranch( "MCe0" );
        b_MCxcore = fChain->GetBranch( "MCxcore" );
        b_MCycore = fChain->GetBranch( "MCycore" );
        if( !bShort )
        {
            b_MCxcos = fChain->GetBranch( "MCxcos" );
            b_MCycos = fChain->GetBranch( "MCycos" );
        }
        else
        {
            b_MCxcos = 0;
            b_MCycos = 0;
        }
        b_MCze = fChain->GetBranch( "MCze" );
        b_MCaz = fChain->GetBranch( "MCaz" );
        b_MCxoff = fChain->GetBranch( "MCxoff" );
        b_MCyoff = fChain->GetBranch( "MCyoff" );
        if( !bShort )
        {
            b_MCxcore_SC = fChain->GetBranch( "MCycore_SC" );
            b_MCycore_SC = fChain->GetBranch( "MCycore_SC" );
            b_MCycore_SC = fChain->GetBranch( "MCycore_SC" );
        }
        else
        {
            b_MCxcore_SC = 0;
            b_MCycore_SC = 0;
            b_MCycore_SC = 0;
        }
        b_MCCorsikaRunID = fChain->GetBranch( "MCCorsikaRunID" );
        b_MCCorsikaShowerID = fChain->GetBranch( "MCCorsikaShowerID" );
        b_MCFirstInteractionHeight = fChain->GetBranch( "MCFirstInteractionHeight" );
        b_MCFirstInteractionDepth = fChain->GetBranch( "MCFirstInteractionDepth" );
        
    }
    
    return kTRUE;
}


void Cshowerpars::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}


#endif                                            // #ifdef Cshowerpars_cxx
