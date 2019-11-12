//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  9 09:01:50 2015 by ROOT version 5.34/32
// from TTree IRFData/test data
// found on file: DESY.d20151114.V5.ID0_0degNIM2.prod3-paranalp05-NN.sensitivityTree.v2.root
//////////////////////////////////////////////////////////

#ifndef cIRFData_h
#define cIRFData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class cIRFData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Char_t          Array[14];
   UInt_t          ArrayID;
   UInt_t          Scaling;
   Int_t           NTel;
   Int_t           NTelType[7];
   Int_t           ObsTime_s;
   Float_t         Offset_deg;
   Float_t         AverageTelescopeDistance[7];
   Float_t         MedianTelescopeDistance[7];
   Float_t         Energy_logTeV[21];
   Float_t         DiffSens[21];
   Float_t         DiffSensError[21];
   Float_t         DiffSensCU[21];
   Float_t         DiffSensCUError[21];
   Float_t         IntSens[21];
   Float_t         IntSensError[21];
   Float_t         IntSensCU[21];
   Float_t         IntSensCUError[21];
   Float_t         AngRes[21];
   Float_t         AngResError[21];
   Float_t         AngRes80[21];
   Float_t         AngRes80Error[21];
   Float_t         ERes[21];
   Float_t         EResError[21];
   Float_t         Ebias[21];
   Float_t         EbiasError[21];
   Float_t         BGRate[21];
   Float_t         BGRateError[21];
   Float_t         ProtRate[21];
   Float_t         ProtRateError[21];
   Float_t         ElecRate[21];
   Float_t         ElecRateError[21];
   Float_t         BGRatePerSqDeg[21];
   Float_t         BGRatePerSqDegError[21];
   Float_t         ProtRateSqDeg[21];
   Float_t         ProtRateSqDegError[21];
   Float_t         ElecRateSqDeg[21];
   Float_t         ElecRateSqDegError[21];
   Float_t         EffectiveArea[21];
   Float_t         EffectiveAreaError[21];
   Float_t         EffectiveAreaEtrue[21];
   Float_t         EffectiveAreaEtrueError[21];
   Float_t         EffectiveArea80[21];
   Float_t         EffectiveArea80Error[21];

   // List of branches
   TBranch        *b_Array;   //!
   TBranch        *b_ArrayID;   //!
   TBranch        *b_Scaling;   //!
   TBranch        *b_NTel;   //!
   TBranch        *b_NTelType;   //!
   TBranch        *b_ObsTime_s;   //!
   TBranch        *b_Offset_deg;   //!
   TBranch        *b_AverageTelescopeDistance;   //!
   TBranch        *b_MedianTelescopeDistance;   //!
   TBranch        *b_Energy_logTeV;   //!
   TBranch        *b_DiffSens;   //!
   TBranch        *b_DiffSensError;   //!
   TBranch        *b_DiffSensCU;   //!
   TBranch        *b_DiffSensCUError;   //!
   TBranch        *b_IntSens;   //!
   TBranch        *b_IntSensError;   //!
   TBranch        *b_IntSensCU;   //!
   TBranch        *b_IntSensCUError;   //!
   TBranch        *b_AngRes;   //!
   TBranch        *b_AngResError;   //!
   TBranch        *b_AngRes80;   //!
   TBranch        *b_AngRes80Error;   //!
   TBranch        *b_ERes;   //!
   TBranch        *b_EResError;   //!
   TBranch        *b_Ebias;   //!
   TBranch        *b_EbiasError;   //!
   TBranch        *b_BGRate;   //!
   TBranch        *b_BGRateError;   //!
   TBranch        *b_ProtRate;   //!
   TBranch        *b_ProtRateError;   //!
   TBranch        *b_ElecRate;   //!
   TBranch        *b_ElecRateError;   //!
   TBranch        *b_BGRatePerSqDeg;   //!
   TBranch        *b_BGRatePerSqDegError;   //!
   TBranch        *b_ProtRateSqDeg;   //!
   TBranch        *b_ProtRateSqDegError;   //!
   TBranch        *b_ElecRateSqDeg;   //!
   TBranch        *b_ElecRateSqDegError;   //!
   TBranch        *b_EffectiveArea;   //!
   TBranch        *b_EffectiveAreaError;   //!
   TBranch        *b_EffectiveAreaEtrue;   //!
   TBranch        *b_EffectiveAreaEtrueError;   //!
   TBranch        *b_EffectiveArea80;   //!
   TBranch        *b_EffectiveArea80Error;   //!

   cIRFData(TTree *tree=0);
   virtual ~cIRFData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef cIRFData_cxx
cIRFData::cIRFData(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DESY.d20151114.V5.ID0_0degNIM2.prod3-paranalp05-NN.sensitivityTree.v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DESY.d20151114.V5.ID0_0degNIM2.prod3-paranalp05-NN.sensitivityTree.v2.root");
      }
      f->GetObject("IRFData",tree);

   }
   Init(tree);
}

cIRFData::~cIRFData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t cIRFData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t cIRFData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void cIRFData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Array", Array, &b_Array);
   fChain->SetBranchAddress("ArrayID", &ArrayID, &b_ArrayID);
   fChain->SetBranchAddress("Scaling", &Scaling, &b_Scaling);
   fChain->SetBranchAddress("NTel", &NTel, &b_NTel);
   fChain->SetBranchAddress("NTelType[7]", NTelType, &b_NTelType);
   fChain->SetBranchAddress("ObsTime_s", &ObsTime_s, &b_ObsTime_s);
   fChain->SetBranchAddress("Offset_deg", &Offset_deg, &b_Offset_deg);
   fChain->SetBranchAddress("AverageTelescopeDistance[7]", AverageTelescopeDistance, &b_AverageTelescopeDistance);
   fChain->SetBranchAddress("MedianTelescopeDistance[7]", MedianTelescopeDistance, &b_MedianTelescopeDistance);
   fChain->SetBranchAddress("Energy_logTeV[21]", Energy_logTeV, &b_Energy_logTeV);
   fChain->SetBranchAddress("DiffSens[21]", DiffSens, &b_DiffSens);
   fChain->SetBranchAddress("DiffSensError[21]", DiffSensError, &b_DiffSensError);
   fChain->SetBranchAddress("DiffSensCU[21]", DiffSensCU, &b_DiffSensCU);
   fChain->SetBranchAddress("DiffSensCUError[21]", DiffSensCUError, &b_DiffSensCUError);
   fChain->SetBranchAddress("IntSens[21]", IntSens, &b_IntSens);
   fChain->SetBranchAddress("IntSensError[21]", IntSensError, &b_IntSensError);
   fChain->SetBranchAddress("IntSensCU[21]", IntSensCU, &b_IntSensCU);
   fChain->SetBranchAddress("IntSensCUError[21]", IntSensCUError, &b_IntSensCUError);
   fChain->SetBranchAddress("AngRes[21]", AngRes, &b_AngRes);
   fChain->SetBranchAddress("AngResError[21]", AngResError, &b_AngResError);
   fChain->SetBranchAddress("AngRes80[21]", AngRes80, &b_AngRes80);
   fChain->SetBranchAddress("AngRes80Error[21]", AngRes80Error, &b_AngRes80Error);
   fChain->SetBranchAddress("ERes[21]", ERes, &b_ERes);
   fChain->SetBranchAddress("EResError[21]", EResError, &b_EResError);
   fChain->SetBranchAddress("Ebias[21]", Ebias, &b_Ebias);
   fChain->SetBranchAddress("EbiasError[21]", EbiasError, &b_EbiasError);
   fChain->SetBranchAddress("BGRate[21]", BGRate, &b_BGRate);
   fChain->SetBranchAddress("BGRateError[21]", BGRateError, &b_BGRateError);
   fChain->SetBranchAddress("ProtRate[21]", ProtRate, &b_ProtRate);
   fChain->SetBranchAddress("ProtRateError[21]", ProtRateError, &b_ProtRateError);
   fChain->SetBranchAddress("ElecRate[21]", ElecRate, &b_ElecRate);
   fChain->SetBranchAddress("ElecRateError[21]", ElecRateError, &b_ElecRateError);
   fChain->SetBranchAddress("BGRatePerSqDeg[21]", BGRatePerSqDeg, &b_BGRatePerSqDeg);
   fChain->SetBranchAddress("BGRatePerSqDegError[21]", BGRatePerSqDegError, &b_BGRatePerSqDegError);
   fChain->SetBranchAddress("ProtRateSqDeg[21]", ProtRateSqDeg, &b_ProtRateSqDeg);
   fChain->SetBranchAddress("ProtRateSqDegError[21]", ProtRateSqDegError, &b_ProtRateSqDegError);
   fChain->SetBranchAddress("ElecRateSqDeg[21]", ElecRateSqDeg, &b_ElecRateSqDeg);
   fChain->SetBranchAddress("ElecRateSqDegError[21]", ElecRateSqDegError, &b_ElecRateSqDegError);
   fChain->SetBranchAddress("EffectiveArea[21]", EffectiveArea, &b_EffectiveArea);
   fChain->SetBranchAddress("EffectiveAreaError[21]", EffectiveAreaError, &b_EffectiveAreaError);
   fChain->SetBranchAddress("EffectiveAreaEtrue[21]", EffectiveAreaEtrue, &b_EffectiveAreaEtrue);
   fChain->SetBranchAddress("EffectiveAreaEtrueError[21]", EffectiveAreaEtrueError, &b_EffectiveAreaEtrueError);
   fChain->SetBranchAddress("EffectiveArea80[21]", EffectiveArea80, &b_EffectiveArea80);
   fChain->SetBranchAddress("EffectiveArea80Error[21]", EffectiveArea80Error, &b_EffectiveArea80Error);
   Notify();
}

Bool_t cIRFData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void cIRFData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t cIRFData::Cut(Long64_t entry)
{
    entry = 0;
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef cIRFData_cxx
