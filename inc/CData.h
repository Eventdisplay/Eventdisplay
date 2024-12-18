//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 22 23:32:42 2006 by ROOT version 5.06/00
// from TTree data/MSWC and energy lookup results
// found on file: 31391.mscw.root
//
//  heavily modified by hand!!
//////////////////////////////////////////////////////////

#ifndef CData_h
#define CData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "VGlobalRunParameter.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <assert.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// reconstruction types
// note: reconstruction types determine which values are written from the mscw root
//       trees
//       e.g. for energy, should it be the lookup table energy
//       e.g. for direction, should it be the disp result, or the classical result
//       (note that not all reconstruction types are still available today)
////////////////////////////////////////////////////////////////////////////////
enum E_ReconstructionType { NOT_SET = -1, GEO = 0, FROGSDIR = 1, FROGS = 2, MODEL3D = 3, ENERGY_ER = 4, NN = 5, TL = 6, DEEPLEARNER = 7 };


class CData
{

    public :
    
        E_ReconstructionType  fReconstructionType ;
        bool            fMC;
        bool            fDeepLearner;
        
        bool            fShort;
        TTree*          fChain;                   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent;                 //!current Tree number in a TChain
        
        // Declaration of leave types
        Int_t           runNumber;
        Int_t           eventNumber;
        Int_t           MJD;
        Double_t        Time;
        Double_t        TelElevation[VDST_MAXTELESCOPES];
        Double_t        TelAzimuth[VDST_MAXTELESCOPES];
        Double_t        TelDec[VDST_MAXTELESCOPES];
        Double_t        TelRA[VDST_MAXTELESCOPES];
        Float_t         ArrayPointing_Elevation;
        Float_t         ArrayPointing_Azimuth;
        // MC parameters
        Int_t           MCprimary;
        Double_t        MCe0;
        Double_t        MCxcore;
        Double_t        MCycore;
        Double_t        MCxcore_SC;
        Double_t        MCycore_SC;
        Double_t        MCxcos;
        Double_t        MCycos;
        Double_t        MCaz;
        Double_t        MCze;
        Double_t        MCxoff;
        Double_t        MCyoff;
        Int_t	       MCCorsikaRunID;
        Int_t	       MCCorsikaShowerID;
        Float_t        MCFirstInteractionHeight;
        Float_t	       MCFirstInteractionDepth;
        
        ULong64_t       LTrig;
        UInt_t          NTrig;
        Int_t           NImages;
        ULong64_t       ImgSel;
        UInt_t          ImgSel_list[VDST_MAXTELESCOPES];
        Int_t           NTtype;
        ULong64_t       TtypeID[VDST_MAXTELESCOPES];
        UInt_t		NImages_Ttype[VDST_MAXTELESCOPES];
        Double_t        img2_ang;
        Double_t        Ze;
        Double_t        Az;
        Double_t        ra;
        Double_t        dec;
        Double_t        Xoff;
        Double_t        Yoff;
        Double_t        Xoff_derot;
        Double_t        Yoff_derot;
        Float_t        Xoff_intersect;
        Float_t        Yoff_intersect;
        Double_t        stdS;
        Double_t        theta2;
        Double_t        Xcore;
        Double_t        Ycore;
        Double_t        Xcore_SC;
        Double_t        Ycore_SC;
        Double_t        stdP;
        Double_t        Chi2;
        Float_t         meanPedvar_Image;
        Float_t         meanPedvar_ImageT[VDST_MAXTELESCOPES];
        Float_t         dist[VDST_MAXTELESCOPES];
        Float_t         size[VDST_MAXTELESCOPES];
        Double_t        size2[VDST_MAXTELESCOPES];
        Double_t        fraclow[VDST_MAXTELESCOPES];
        Float_t         loss[VDST_MAXTELESCOPES];
        Float_t         fui[VDST_MAXTELESCOPES];
        Float_t         cross[VDST_MAXTELESCOPES];
        Double_t        max1[VDST_MAXTELESCOPES];
        Double_t        max2[VDST_MAXTELESCOPES];
        Double_t        max3[VDST_MAXTELESCOPES];
        Int_t           maxindex1[VDST_MAXTELESCOPES];
        Int_t           maxindex2[VDST_MAXTELESCOPES];
        Int_t           maxindex3[VDST_MAXTELESCOPES];
        Float_t         width[VDST_MAXTELESCOPES];
        Float_t         length[VDST_MAXTELESCOPES];
        Int_t           ntubes[VDST_MAXTELESCOPES];
        Int_t           ntubesBNI[VDST_MAXTELESCOPES];
        UShort_t        nsat[VDST_MAXTELESCOPES];
        UShort_t        nlowgain[VDST_MAXTELESCOPES];
        Double_t        alpha[VDST_MAXTELESCOPES];
        Double_t        los[VDST_MAXTELESCOPES];
        Float_t         asym[VDST_MAXTELESCOPES];
        Float_t         cen_x[VDST_MAXTELESCOPES];
        Float_t         cen_y[VDST_MAXTELESCOPES];
        Float_t         cosphi[VDST_MAXTELESCOPES];
        Float_t         sinphi[VDST_MAXTELESCOPES];
        Float_t         tgrad_x[VDST_MAXTELESCOPES];
        Double_t        tchisq_x[VDST_MAXTELESCOPES];
        Int_t           Fitstat[VDST_MAXTELESCOPES];
        Double_t        R[VDST_MAXTELESCOPES];
        Double_t        MSCWT[VDST_MAXTELESCOPES];
        Double_t        MSCLT[VDST_MAXTELESCOPES];
        Double_t        ES[VDST_MAXTELESCOPES];
        Int_t           NMSCW;
        Double_t        MSCW;
        Double_t        MSCL;
        Float_t         MWR;
        Float_t         MLR;
        Double_t        Erec;
        Double_t        EChi2;
        Double_t        dE;        // Error on Erec
        Double_t        ErecS;
        Double_t        EChi2S;
        Double_t        dES;       // Error on ErecS
        Double_t        SizeSecondMax;
        Double_t        theta2_All[25];
        Float_t         EmissionHeight;
        Float_t         EmissionHeightChi2;
        UInt_t          NTelPairs;
        //[NTelPairs]
        Float_t         EmissionHeightT[VDST_MAXTELESCOPES* VDST_MAXTELESCOPES];
        Double_t        DispDiff;  // from disp method
        Float_t         DispAbsSumWeigth;
        UInt_t          DispNImages;
        // Deep Learner Parameters
        Double_t         dl_gammaness;
        Bool_t         dl_isGamma;
        
        // List of branches
        TBranch*        b_runNumber;              //!
        TBranch*        b_eventNumber;            //!
        TBranch*        b_MJD;                    //!
        TBranch*        b_Time;                   //!
        TBranch*        b_TelElevation;           //!
        TBranch*        b_TelAzimuth;             //!
        TBranch*        b_TelDec;                 //!
        TBranch*        b_TelRA;                  //!
        TBranch*        b_ArrayPointing_Elevation;  //!
        TBranch*        b_ArrayPointing_Azimuth;  //!
        // MC parameter
        TBranch*        b_MCprimary;
        TBranch*        b_MCe0;                   //!
        TBranch*        b_MCxcore;                //!
        TBranch*        b_MCycore;                //!
        TBranch*        b_MCxcore_SC;             //!
        TBranch*        b_MCycore_SC;             //!
        TBranch*        b_MCxcos;                 //!
        TBranch*        b_MCycos;                 //!
        TBranch*        b_MCaz;                   //!
        TBranch*        b_MCze;                   //!
        TBranch*        b_MCxoff;                 //!
        TBranch*        b_MCyoff;                 //!
        TBranch*        b_MCCorsikaRunID;         //!
        TBranch*        b_MCCorsikaShowerID;      //!
        TBranch*        b_MCFirstInteractionHeight;    //!
        TBranch*        b_MCFirstInteractionDepth;     //!
        
        TBranch*        b_LTrig;                  //!
        TBranch*        b_NTrig;                  //!
        TBranch*        b_NImages;                //!
        TBranch*        b_ImgSel;                 //!
        TBranch*        b_img2_ang;               //!
        TBranch*        b_Ze;                     //!
        TBranch*        b_Az;                     //!
        TBranch*        b_ra;                     //!
        TBranch*        b_dec;                    //!
        TBranch*        b_Xoff;                   //!
        TBranch*        b_Yoff;                   //!
        TBranch*        b_Xoff_derot;             //!
        TBranch*        b_Yoff_derot;             //!
        TBranch*        b_Xoff_intersect;             //!
        TBranch*        b_Yoff_intersect;             //!
        TBranch*        b_stdS;                   //!
        TBranch*        b_theta2;                 //!
        TBranch*        b_Xcore;                  //!
        TBranch*        b_Ycore;                  //!
        TBranch*        b_Xcore_SC;               //!
        TBranch*        b_Ycore_SC;               //!
        TBranch*        b_stdP;                   //!
        TBranch*        b_Chi2;                   //!
        TBranch*        b_meanPedvar_Image;       //!
        TBranch*        b_meanPedvar_ImageT;      //!
        TBranch*        b_dist;                   //!
        TBranch*        b_size;                   //!
        TBranch*        b_size2;                  //!
        TBranch*        b_fraclow;                //!
        TBranch*        b_loss;                   //!
        TBranch*        b_fui;                   //!
        TBranch*        b_cross;                   //!
        TBranch*        b_max1;                   //!
        TBranch*        b_max2;                   //!
        TBranch*        b_max3;                   //!
        TBranch*        b_maxindex1;              //!
        TBranch*        b_maxindex2;              //!
        TBranch*        b_maxindex3;              //!
        TBranch*        b_width;                  //!
        TBranch*        b_length;                 //!
        TBranch*        b_ntubes;                 //!
        TBranch*        b_ntubesBNI;              //!
        TBranch*        b_nsat;                   //!
        TBranch*        b_nlowgain;               //!
        TBranch*        b_alpha;                  //!
        TBranch*        b_los;                    //!
        TBranch*        b_asym;                   //!
        TBranch*        b_cen_x;                  //!
        TBranch*        b_cen_y;                  //!
        TBranch*        b_cosphi;                 //!
        TBranch*        b_sinphi;                 //!
        TBranch*        b_tgrad_x;                //!
        TBranch*        b_Fitstat;                //!
        TBranch*        b_tchisq_x;               //!
        TBranch*        b_R;                      //!
        TBranch*        b_MSCWT;                  //!
        TBranch*        b_MSCLT;                  //!
        TBranch*        b_ES;                     //!
        TBranch*        b_NMSCW;                  //!
        TBranch*        b_MSCW;                   //!
        TBranch*        b_MSCL;                   //!
        TBranch*        b_MWR;                    //!
        TBranch*        b_MLR;                    //!
        TBranch*        b_Erec;                   //!
        TBranch*        b_EChi2;                  //!
        TBranch*        b_ErecS;                  //!
        TBranch*        b_EChi2S;                 //!
        TBranch*        b_SizeSecondMax;          //!
        TBranch*        b_theta2_All;             //!
        TBranch*        b_EmissionHeight;         //!
        TBranch*        b_EmissionHeightChi2;     //!
        TBranch*        b_NTelPairs;              //!
        TBranch*        b_EmissionHeightT;        //!
        TBranch*        b_DispDiff; //disp
        TBranch*        b_DispAbsSumWeigth;
        TBranch*        b_DispNImages; //disp
        // deep learner parameters
        TBranch*        b_dl_gammaness;             //!
        TBranch*        b_dl_isGamma;             //!
        
        CData( TTree* tree = 0, bool bMC = false, bool bShort = false );
        virtual ~CData();
        virtual Int_t    Cut( Long64_t entry );
        virtual Int_t    GetEntry( Long64_t entry );
        virtual Long64_t LoadTree( Long64_t entry );
        virtual void     Init( TTree* tree );
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show( Long64_t entry = -1 );
        bool             isDeepLearner()
        {
            return fDeepLearner;
        }
        bool             isMC()
        {
            return fMC;
        }
        void setReconstructionType( E_ReconstructionType type )
        {
            fReconstructionType  = type;
            
            // general check if analysis types
            if( fReconstructionType != GEO
                    && fReconstructionType != NN
                    && fReconstructionType != TL
                    && fReconstructionType != ENERGY_ER
                    && fReconstructionType != DEEPLEARNER )
            {
                cout << "CData::setReconstructionType Error: unknown analysis type: ";
                cout << fReconstructionType << endl;
                exit( EXIT_FAILURE );
            }
        }
        double getEnergy_TeV()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Erec;
            }
            if( fReconstructionType  == GEO )
            {
                return ErecS;
            }
            else
            {
                return ErecS;
            }
        }
        double getEnergy_Log10()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                if( Erec <= 0 )
                {
                    return -99;
                }
                else
                {
                    return log10( Erec );
                }
            }
            if( ErecS <= 0 )
            {
                return -99;
            }
            else
            {
                return log10( ErecS );
            }
        }
        double getXcore_M()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Xcore;
            }
            return Xcore;
        }
        double getYcore_M()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Ycore;
            }
            return Ycore;
        }
        double getXoff()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Xoff;
            }
            return Xoff;
        }
        double getYoff()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Yoff;
            }
            return Yoff;
        }
        double getXoff_derot()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Xoff_derot;
            }
            return Xoff_derot;
        }
        double getYoff_derot()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return Yoff_derot;
            }
            return Yoff_derot;
        }
        double getEnergyChi2()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return EChi2;
            }
            return EChi2S;
        }
        
        double getEnergyDelta()
        {
            if( fReconstructionType  == ENERGY_ER )
            {
                return dE;
            }
            return dES;
        }
        double getChi2()
        {
            return Chi2;
        }
        float getDirectionReconstructionDifference()
        {
            if( Xoff > -998. && Yoff > -998. && Xoff_intersect > -998. && Yoff_intersect > -998. )
            {
                if( ( Xoff - Xoff_intersect ) * ( Xoff - Xoff_intersect ) + ( Yoff - Yoff_intersect ) * ( Yoff - Yoff_intersect ) > 0. )
                {
                    return log10( sqrt( ( Xoff - Xoff_intersect ) * ( Xoff - Xoff_intersect ) + ( Yoff - Yoff_intersect ) * ( Yoff - Yoff_intersect ) ) );
                }
            }
            return -999.;
        }
        int getNImages()
        {
            return NImages ;
        }
        ULong64_t getImgSel()
        {
            return ImgSel;
        }
        
        UInt_t* getImgSel_list()
        {
            return ImgSel_list; //preliminary
        }
        
        double getZe()
        {
            return Ze;
        }
        double getAz()
        {
            return Az;
        }
        double* getR()
        {
            return R;
        }
        double* getEnergy_per_telescope()
        {
            return ES;
        }
        unsigned getNImages_used_in_EnergyReconstruction()
        {
            unsigned int n = 0;
            for( int i = 0; i < NImages; i++ )
            {
                if( ES[i] > 0. )
                {
                    n++;
                }
            }
            return n;
        }
        
        
};
#endif

#ifdef CData_cxx

CData::CData( TTree* tree, bool bMC, bool bShort )
{
    fMC = bMC;
    fShort = bShort;
    fDeepLearner = false;
    fReconstructionType = GEO;
    Init( tree );
}


CData::~CData()
{
    if( !fChain )
    {
        return;
    }
    delete fChain->GetCurrentFile();
}


Int_t CData::GetEntry( Long64_t entry )
{
    // Read contents of entry.
    if( !fChain )
    {
        return 0;
    }
    
    int a = fChain->GetEntry( entry );
    
    return a;
}


Long64_t CData::LoadTree( Long64_t entry )
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


void CData::Init( TTree* tree )
{

    // Set branch addresses
    if( tree == 0 )
    {
        return;
    }
    
    // test if this is a MC file
    if( tree->GetBranchStatus( "MCe0" ) )
    {
        fMC = true;
    }
    // test if deep learner variables exists
    if( tree->GetBranchStatus( "dl_gammaness" ) )
    {
        fDeepLearner = true;
    }
    
    fChain = tree;
    fCurrent = -1;
    
    fChain->SetBranchAddress( "runNumber", &runNumber );
    fChain->SetBranchAddress( "eventNumber", &eventNumber );
    if( !fShort )
    {
        fChain->SetBranchAddress( "MJD", &MJD );
        fChain->SetBranchAddress( "Time", &Time );
    }
    else
    {
        MJD = 0;
        Time = 0;
    }
    if( fChain->GetBranchStatus( "TelElevation" ) && fChain->GetBranchStatus( "TelAzimuth" ) )
    {
        fChain->SetBranchAddress( "TelElevation", TelElevation );
        fChain->SetBranchAddress( "TelAzimuth", TelAzimuth );
    }
    else
    {
        for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            TelElevation[i] = 0.;
            TelAzimuth[0] = 0.;
        }
    }
    
    if( fChain->GetBranchStatus( "ArrayPointing_Azimuth" ) )
    {
        fChain->SetBranchAddress( "ArrayPointing_Azimuth", &ArrayPointing_Azimuth );
    }
    else
    {
        ArrayPointing_Azimuth = 0.;
    }
    if( fChain->GetBranchStatus( "ArrayPointing_Elevation" ) )
    {
        fChain->SetBranchAddress( "ArrayPointing_Elevation", &ArrayPointing_Elevation );
    }
    else
    {
        ArrayPointing_Elevation = 0.;
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
            TelRA[i] = 0.;
        }
    }
    
    // MC tree
    if( fMC )
    {
        if( fChain->GetBranch( "MCprimary" ) )
        {
            fChain->SetBranchAddress( "MCprimary", &MCprimary );
        }
        else
        {
            MCprimary = 0;
        }
        fChain->SetBranchAddress( "MCe0", &MCe0 );
        fChain->SetBranchAddress( "MCxcore", &MCxcore );
        fChain->SetBranchAddress( "MCycore", &MCycore );
        if( !fShort )
        {
            fChain->SetBranchAddress( "MCxcore_SC", &MCxcore_SC );
            fChain->SetBranchAddress( "MCycore_SC", &MCycore_SC );
            fChain->SetBranchAddress( "MCxcos", &MCxcos );
            fChain->SetBranchAddress( "MCycos", &MCycos );
        }
        else
        {
            MCxcore_SC = MCycore_SC = MCxcos = MCycos = 0.;
        }
        fChain->SetBranchAddress( "MCaz", &MCaz );
        fChain->SetBranchAddress( "MCze", &MCze );
        fChain->SetBranchAddress( "MCxoff", &MCxoff );
        fChain->SetBranchAddress( "MCyoff", &MCyoff );
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
    else
    {
        MCprimary = -99;
        MCe0 = 0.;
        MCxcore = 0.;
        MCycore = 0.;
        MCxcore_SC = 0.;
        MCycore_SC = 0.;
        MCxcos = 0.;
        MCycos = 0.;
        MCaz = 0.;
        MCze = 0.;
        MCxoff = 0.;
        MCyoff = 0.;
        MCCorsikaShowerID = 0;
        MCFirstInteractionHeight = 0;
        MCFirstInteractionDepth = 0;
        MCCorsikaRunID = 0;
        
    }
    
    
    if( fChain->GetBranchStatus( "LTrig" ) )
    {
        fChain->SetBranchAddress( "LTrig", &LTrig );
        fChain->SetBranchAddress( "NTrig", &NTrig );
    }
    else
    {
        LTrig = 0;
        NTrig = 0;
    }
    fChain->SetBranchAddress( "NImages", &NImages );
    fChain->SetBranchAddress( "ImgSel", &ImgSel );
    if( fChain->GetBranchStatus( "img2_ang" ) )
    {
        fChain->SetBranchAddress( "img2_ang", &img2_ang );
    }
    else
    {
        img2_ang = 0.;
    }
    if( fChain->GetBranchStatus( "Ze" ) &&  fChain->GetBranchStatus( "Az" ) )
    {
        fChain->SetBranchAddress( "Ze", &Ze );
        fChain->SetBranchAddress( "Az", &Az );
    }
    else
    {
        Ze = 0.;
        Az = 0.;
    }
    if( !fShort )
    {
        fChain->SetBranchAddress( "ra", &ra );
        fChain->SetBranchAddress( "dec", &dec );
    }
    else
    {
        ra = dec = 0.;
    }
    fChain->SetBranchAddress( "Xoff", &Xoff );
    fChain->SetBranchAddress( "Yoff", &Yoff );
    if( fChain->GetBranchStatus( "Xoff_derot" ) )
    {
        fChain->SetBranchAddress( "Xoff_derot", &Xoff_derot );
    }
    if( fChain->GetBranchStatus( "Yoff_derot" ) )
    {
        fChain->SetBranchAddress( "Yoff_derot", &Yoff_derot );
    }
    if( fChain->GetBranchStatus( "Xoff_intersect" ) )
    {
        fChain->SetBranchAddress( "Xoff_intersect", &Xoff_intersect );
    }
    else
    {
        Xoff_intersect = 0.;
    }
    if( fChain->GetBranchStatus( "Yoff_intersect" ) )
    {
        fChain->SetBranchAddress( "Yoff_intersect", &Yoff_intersect );
    }
    else
    {
        Yoff_intersect = 0.;
    }
    
    if( fChain->GetBranchStatus( "stdS" ) )
    {
        fChain->SetBranchAddress( "stdS", &stdS );
    }
    else
    {
        stdS = 0.;
    }
    if( !fShort )
    {
        fChain->SetBranchAddress( "theta2", &theta2 );
    }
    else
    {
        theta2 = 0.;
    }
    fChain->SetBranchAddress( "Xcore", &Xcore );
    fChain->SetBranchAddress( "Ycore", &Ycore );
    if( !fShort )
    {
        fChain->SetBranchAddress( "Xcore_SC", &Xcore_SC );
        fChain->SetBranchAddress( "Ycore_SC", &Ycore_SC );
    }
    else
    {
        Xcore_SC = Ycore_SC = 0.;
    }
    if( fChain->GetBranchStatus( "stdP" ) )
    {
        fChain->SetBranchAddress( "stdP", &stdP );
    }
    else
    {
        stdP = 0.;
    }
    
    fChain->SetBranchAddress( "Chi2", &Chi2 );
    if( fChain->GetBranchStatus( "meanPedvar_Image" ) )
    {
        fChain->SetBranchAddress( "meanPedvar_Image", &meanPedvar_Image );
    }
    else
    {
        meanPedvar_Image = 0.;
    }
    if( !fShort )
    {
        fChain->SetBranchAddress( "meanPedvar_ImageT", meanPedvar_ImageT );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            meanPedvar_ImageT[i] = 0.;
        }
    }
    
    fChain->SetBranchAddress( "SizeSecondMax", &SizeSecondMax );
    
    if( fChain->GetBranchStatus( "theta2_All" ) )
    {
        fChain->SetBranchAddress( "theta2_All", &theta2_All );
    }
    else
    {
        for( unsigned int dex = 0; dex < 25; dex++ )
        {
            theta2_All[dex] = 99.0;
        }
    }
    
    if( fChain->GetBranchStatus( "NTtype" ) )
    {
        fChain->SetBranchAddress( "ImgSel_list", ImgSel_list );
        fChain->SetBranchAddress( "NTtype", &NTtype );
        fChain->SetBranchAddress( "NImages_Ttype", NImages_Ttype );
    }
    else
    {
        NTtype = 0;
        for( unsigned int tt = 0; tt < VDST_MAXTELESCOPES; tt++ )
        {
            NImages_Ttype[tt] = 0;
        }
    }
    if( fChain->GetBranchStatus( "TtypeID" ) )
    {
        fChain->SetBranchAddress( "TtypeID", TtypeID );
    }
    else
    {
        for( unsigned int tt = 0; tt < VDST_MAXTELESCOPES; tt++ )
        {
            TtypeID[tt] = 0;
        }
    }
    fChain->SetBranchAddress( "dist", dist );
    fChain->SetBranchAddress( "size", size );
    fChain->SetBranchAddress( "loss", loss );
    fChain->SetBranchAddress( "asym", asym );
    fChain->SetBranchAddress( "tgrad_x", tgrad_x );
    if( fChain->GetBranchStatus( "fui" ) )
    {
        fChain->SetBranchAddress( "fui", fui );
    }
    else
    {
        for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            fui[i] = 0.;
        }
    }
    if( fChain->GetBranchStatus( "cross" ) )
    {
        fChain->SetBranchAddress( "cross", cross );
    }
    else
    {
        for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            cross[i] = 0.;
        }
    }
    
    
    if( !fShort )
    {
        if( fChain->GetBranchStatus( "size2" ) )
        {
            fChain->SetBranchAddress( "size2", size2 );
        }
        else
        {
            for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
            {
                size2[i] = 0.;
            }
        }
        if( fChain->GetBranchStatus( "fracLow" ) )
        {
            fChain->SetBranchAddress( "fracLow", fraclow );
        }
        else
        {
            for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
            {
                fraclow[i] = 0.;
            }
        }
        fChain->SetBranchAddress( "max1", max1 );
        fChain->SetBranchAddress( "max2", max2 );
        fChain->SetBranchAddress( "max3", max3 );
        fChain->SetBranchAddress( "maxindex1", maxindex1 );
        fChain->SetBranchAddress( "maxindex2", maxindex2 );
        fChain->SetBranchAddress( "maxindex3", maxindex3 );
        fChain->SetBranchAddress( "width", width );
        fChain->SetBranchAddress( "length", length );
        fChain->SetBranchAddress( "ntubes", ntubes );
        fChain->SetBranchAddress( "nsat", nsat );
        if( fChain->GetBranchStatus( "nlowgain" ) )
        {
            fChain->SetBranchAddress( "nlowgain", nlowgain );
        }
        else
        {
            for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
            {
                nlowgain[i] = 0;
            }
        }
        if( fChain->GetBranchStatus( "ntubesBNI" ) )
        {
            fChain->SetBranchAddress( "ntubesBNI", ntubesBNI );
        }
        else
        {
            for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
            {
                ntubesBNI[i] = 0;
            }
        }
        fChain->SetBranchAddress( "alpha", alpha );
        fChain->SetBranchAddress( "los", los );
        fChain->SetBranchAddress( "cen_x", cen_x );
        fChain->SetBranchAddress( "cen_y", cen_y );
        fChain->SetBranchAddress( "cosphi", cosphi );
        fChain->SetBranchAddress( "sinphi", sinphi );
        fChain->SetBranchAddress( "Fitstat", Fitstat );
        fChain->SetBranchAddress( "tchisq_x", tchisq_x );
    }
    else
    {
        for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            size2[i] = 0.;
            fraclow[i] = 0.;
            max1[i] = 0.;
            max2[i] = 0.;
            max3[i] = 0.;
            maxindex1[i] = 0;
            maxindex2[i] = 0;
            maxindex3[i] = 0;
            width[i] = 0.;
            length[i] = 0.;
            ntubes[i] = 0;
            nsat[i] = 0;
            nlowgain[i] = 0;
            ntubesBNI[i] = 0;
            alpha[i] = 0.;
            los[i] = 0.;
            cen_x[i] = 0.;
            cen_y[i] = 0.;
            cosphi[i] = 0.;
            sinphi[i] = 0.;
            Fitstat[i] = 0;
            tchisq_x[i] = 0.;
        }
    }
    fChain->SetBranchAddress( "R", R );
    if( !fShort )
    {
        fChain->SetBranchAddress( "MSCWT", MSCWT );
        fChain->SetBranchAddress( "MSCLT", MSCLT );
    }
    else
    {
        for( int i = 0; i < VDST_MAXTELESCOPES; i++ )
        {
            MSCWT[i] = 0.;
            MSCLT[i] = 0.;
        }
    }
    fChain->SetBranchAddress( "ES", ES );
    if( !fShort )
    {
        fChain->SetBranchAddress( "NMSCW", &NMSCW );
    }
    else
    {
        NMSCW = 0;
    }
    fChain->SetBranchAddress( "MSCW", &MSCW );
    fChain->SetBranchAddress( "MSCL", &MSCL );
    fChain->SetBranchAddress( "MWR", &MWR );
    fChain->SetBranchAddress( "MLR", &MLR );
    fChain->SetBranchAddress( "ErecS", &ErecS );
    fChain->SetBranchAddress( "EChi2S", &EChi2S );
    fChain->SetBranchAddress( "dES", &dES );
    if( fChain->GetBranchStatus( "Erec" ) )
    {
        fChain->SetBranchAddress( "Erec", &Erec );
        fChain->SetBranchAddress( "EChi2", &EChi2 );
        fChain->SetBranchAddress( "dES", &dES );
    }
    else
    {
        Erec = -99.;
        EChi2 = -99.;
        dES = -99.;
    }
    EmissionHeight = -99.;
    fChain->SetBranchAddress( "EmissionHeight", &EmissionHeight );
    fChain->SetBranchAddress( "EmissionHeightChi2", &EmissionHeightChi2 );
    if( fChain->GetBranchStatus( "NTelPairs" ) )
    {
        fChain->SetBranchAddress( "NTelPairs", &NTelPairs );
    }
    else
    {
        NTelPairs = 0;
    }
    if( !fShort )
    {
        fChain->SetBranchAddress( "EmissionHeightT", EmissionHeightT );
    }
    else
    {
        for( unsigned int i = 0; i < VDST_MAXTELESCOPES * VDST_MAXTELESCOPES; i++ )
        {
            EmissionHeightT[i] = 0.;
        }
    }
    if( fChain->GetBranchStatus( "DispDiff" ) )
    {
        fChain->SetBranchAddress( "DispDiff", &DispDiff );
    }
    else
    {
        DispDiff = 0.;
    }
    if( fChain->GetBranchStatus( "DispAbsSumWeigth" ) )
    {
        fChain->SetBranchAddress( "DispAbsSumWeigth", &DispAbsSumWeigth );
    }
    else
    {
        DispAbsSumWeigth = 0.;
    }
    if( fChain->GetBranchStatus( "DispNImages" ) )
    {
        fChain->SetBranchAddress( "DispNImages", &DispNImages );
    }
    else
    {
        DispNImages = 0;
    }
    if( fDeepLearner )
    {
        fChain->SetBranchAddress( "dl_gammaness", &dl_gammaness );
        fChain->SetBranchAddress( "dl_isGamma", &dl_isGamma );
    }
    else
    {
        dl_gammaness = 0.;
        dl_isGamma = 0;
    }
    
    Notify();
}


Bool_t CData::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. Typically here the branch pointers
    // will be retrieved. It is normaly not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed.
    
    // Get branch pointers
    b_runNumber = fChain->GetBranch( "runNumber" );
    b_eventNumber = fChain->GetBranch( "eventNumber" );
    b_MJD = fChain->GetBranch( "MJD" );
    b_Time = fChain->GetBranch( "Time" );
    b_TelElevation = fChain->GetBranch( "TelElevation" );
    b_ArrayPointing_Elevation = fChain->GetBranch( "ArrayPointing_Elevation" );
    b_ArrayPointing_Azimuth = fChain->GetBranch( "ArrayPointing_Azimuth" );
    b_TelAzimuth = fChain->GetBranch( "TelAzimuth" );
    b_TelDec = fChain->GetBranch( "TelDec" );
    b_TelRA = fChain->GetBranch( "TelRA" );
    
    if( fMC )
    {
        b_MCprimary = fChain->GetBranch( "MCprimary" );
        b_MCe0 = fChain->GetBranch( "MCe0" );
        b_MCxcore  = fChain->GetBranch( "MCxcore" );
        b_MCycore  = fChain->GetBranch( "MCycore" );
        b_MCxcore_SC  = fChain->GetBranch( "MCxcore_SC" );
        b_MCycore_SC  = fChain->GetBranch( "MCycore_SC" );
        b_MCxcos = fChain->GetBranch( "MCxcos" );
        b_MCycos = fChain->GetBranch( "MCycos" );
        b_MCaz = fChain->GetBranch( "MCaz" );
        b_MCze = fChain->GetBranch( "MCze" );
        b_MCxoff = fChain->GetBranch( "MCxoff" );
        b_MCyoff = fChain->GetBranch( "MCyoff" );
        b_MCCorsikaRunID = fChain->GetBranch( "MCCorsikaRunID" );
        b_MCCorsikaShowerID = fChain->GetBranch( "MCCorsikaShowerID" );
        b_MCFirstInteractionHeight = fChain->GetBranch( "MCFirstInteractionHeight" );
        b_MCFirstInteractionDepth = fChain->GetBranch( "MCFirstInteractionDepth" );
    }
    
    b_LTrig = fChain->GetBranch( "LTrig" );
    b_NTrig = fChain->GetBranch( "NTrig" );
    b_NImages = fChain->GetBranch( "NImages" );
    b_ImgSel = fChain->GetBranch( "ImgSel" );
    b_img2_ang = fChain->GetBranch( "img2_ang" );
    b_Ze = fChain->GetBranch( "Ze" );
    b_Az = fChain->GetBranch( "Az" );
    b_ra = fChain->GetBranch( "ra" );
    b_dec = fChain->GetBranch( "dec" );
    b_Xoff_derot = fChain->GetBranch( "Xoff_derot" );
    b_Yoff_derot = fChain->GetBranch( "Yoff_derot" );
    b_Xoff_intersect = fChain->GetBranch( "Xoff_intersect" );
    b_Yoff_intersect = fChain->GetBranch( "Yoff_intersect" );
    b_Xoff = fChain->GetBranch( "Xoff" );
    b_Yoff = fChain->GetBranch( "Yoff" );
    b_stdS = fChain->GetBranch( "stdS" );
    b_theta2 = fChain->GetBranch( "theta2" );
    b_Xcore = fChain->GetBranch( "Xcore" );
    b_Ycore = fChain->GetBranch( "Ycore" );
    b_Xcore_SC = fChain->GetBranch( "Xcore_SC" );
    b_Ycore_SC = fChain->GetBranch( "Ycore_SC" );
    b_stdP = fChain->GetBranch( "stdP" );
    b_Chi2 = fChain->GetBranch( "Chi2" );
    b_meanPedvar_Image = fChain->GetBranch( "meanPedvar_Image" );
    b_meanPedvar_ImageT = fChain->GetBranch( "meanPedvar_ImageT" );
    
    b_SizeSecondMax = fChain->GetBranch( "SizeSecondMax" );
    
    b_theta2_All = fChain->GetBranch( "theta2_All" );
    b_dist = fChain->GetBranch( "dist" );
    b_size = fChain->GetBranch( "size" );
    b_size2 = fChain->GetBranch( "size2" );
    b_fraclow = fChain->GetBranch( "fraclow" );
    b_max1 = fChain->GetBranch( "max1" );
    b_max2 = fChain->GetBranch( "max2" );
    b_max3 = fChain->GetBranch( "max3" );
    b_maxindex1 = fChain->GetBranch( "maxindex1" );
    b_maxindex2 = fChain->GetBranch( "maxindex2" );
    b_maxindex3 = fChain->GetBranch( "maxindex3" );
    b_width = fChain->GetBranch( "width" );
    b_length = fChain->GetBranch( "length" );
    b_ntubes = fChain->GetBranch( "ntubes" );
    b_ntubesBNI = fChain->GetBranch( "ntubesBNI" );
    b_alpha = fChain->GetBranch( "alpha" );
    b_los = fChain->GetBranch( "los" );
    b_asym = fChain->GetBranch( "asym" );
    b_cen_x = fChain->GetBranch( "cen_x" );
    b_cen_y = fChain->GetBranch( "cen_y" );
    b_cosphi = fChain->GetBranch( "cosphi" );
    b_sinphi = fChain->GetBranch( "sinphi" );
    b_tgrad_x = fChain->GetBranch( "tgrad_x" );
    b_Fitstat = fChain->GetBranch( "Fitstat" );
    b_tchisq_x = fChain->GetBranch( "tchisq_x" );
    b_R = fChain->GetBranch( "R" );
    b_MSCWT = fChain->GetBranch( "MSCWT" );
    b_MSCLT = fChain->GetBranch( "MSCLT" );
    b_ES = fChain->GetBranch( "ES" );
    b_NMSCW = fChain->GetBranch( "NMSCW" );
    b_MSCW = fChain->GetBranch( "MSCW" );
    b_MSCL = fChain->GetBranch( "MSCL" );
    b_MWR = fChain->GetBranch( "MWR" );
    b_MLR = fChain->GetBranch( "MLR" );
    b_Erec = fChain->GetBranch( "Erec" );
    b_EChi2 = fChain->GetBranch( "EChi2" );
    b_ErecS = fChain->GetBranch( "ErecS" );
    b_EChi2S = fChain->GetBranch( "EChi2S" );
    b_EmissionHeight = fChain->GetBranch( "EmissionHeight" );
    b_EmissionHeightChi2 = fChain->GetBranch( "EmissionHeightChi2" );
    b_NTelPairs = fChain->GetBranch( "NTelPairs" );
    b_EmissionHeightT = fChain->GetBranch( "EmissionHeightT" );
    if( fChain->GetBranchStatus( "DispDiff" ) )
    {
        b_DispDiff = fChain->GetBranch( "DispDiff" );
    }
    else
    {
        b_DispDiff = 0;
    }
    if( fChain->GetBranchStatus( "DispAbsSumWeigth" ) )
    {
        b_DispAbsSumWeigth = fChain->GetBranch( "DispAbsSumWeigth" );
    }
    else
    {
        b_DispAbsSumWeigth = 0;
    }
    if( fChain->GetBranchStatus( "DispNImages" ) )
    {
        b_DispNImages = fChain->GetBranch( "DispNImages" );
    }
    else
    {
        b_DispNImages = 0;
    }
    if( fDeepLearner )
    {
        b_dl_gammaness = fChain->GetBranch( "dl_gammaness" );
        b_dl_isGamma = fChain->GetBranch( "dl_isGamma" );
    }
    else
    {
        b_dl_gammaness = 0;
        b_dl_isGamma = 0;
    }
    
    return kTRUE;
}


void CData::Show( Long64_t entry )
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if( !fChain )
    {
        return;
    }
    fChain->Show( entry );
}


Int_t CData::Cut( Long64_t entry )
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    entry = 0;
    
    return 1;
}
#endif                                            // #ifdef CData_cxx
