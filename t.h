//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 29 18:35:28 2015 by ROOT version 6.02/05
// from TTree t/t
// found on file: ntuple_output.root
//////////////////////////////////////////////////////////

#ifndef t_h
#define t_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class t {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           NJets;
   Float_t         JetE[25];
   Float_t         JetPx[25];
   Float_t         JetPy[25];
   Float_t         JetPz[25];
   Float_t         JetPt[25];
   Float_t         JetEt[25];
   Float_t         JetEta[25];
   Float_t         JetPhi[25];
   Float_t         EleE[25];
   Float_t         ElePx[25];
   Float_t         ElePy[25];
   Float_t         ElePz[25];
   Float_t         ElePt[25];
   Float_t         EleEt[25];
   Float_t         EleEta[25];
   Float_t         ElePhi[25];
   Float_t         PhE[25];
   Float_t         PhPx[25];
   Float_t         PhPy[25];
   Float_t         PhPz[25];
   Float_t         PhPt[25];
   Float_t         PhEt[25];
   Float_t         PhEta[25];
   Float_t         PhPhi[25];
   Float_t         MuE[25];
   Float_t         MuPx[25];
   Float_t         MuPy[25];
   Float_t         MuPz[25];
   Float_t         MuPt[25];
   Float_t         MuEt[25];
   Float_t         MuEta[25];
   Float_t         MuPhi[25];
   Float_t         ST;
   Float_t         mBH;
   Float_t         Met;
   Float_t         MetPx;
   Float_t         MetPy;
   Float_t         MetPhi;
   Float_t         Sphericity;
   Float_t         Jet[4];
   Float_t         Ele[4];
   Float_t         Mu[4];
   Float_t         Ph[4];
   Int_t           NPV;
   Int_t           NTracks;
   Int_t           NElectrons;
   Int_t           NPhotons;
   Int_t           NMuons;
   Char_t          NoScrap;
   Int_t           Multiplicity;
   Char_t          isLeptonPhoton;
   Char_t          isEleChannel;
   Char_t          isMuChannel;
   Char_t          isPhChannel;
   Float_t         ResJetEta;
   Float_t         ResJetPhi;
   Float_t         ResJetM;
   Float_t         ResJetPt;
   Float_t         ResJetEt;
   Float_t         ResEleEta;
   Float_t         ResElePhi;
   Float_t         ResEleM;
   Float_t         ResElePt;
   Float_t         ResPhEta;
   Float_t         ResPhPhi;
   Float_t         ResPhM;
   Float_t         ResPhPt;
   Float_t         ResMuEta;
   Float_t         ResMuPhi;
   Float_t         ResMuM;
   Float_t         ResMuPt;
   Float_t         ResLepEta;
   Float_t         ResLepPhi;
   Float_t         ResLepM;
   Float_t         ResLepPt;
   Float_t         ResObjEta;
   Float_t         ResObjPhi;
   Float_t         ResObjM;
   Float_t         ResObjPt;
   Int_t           runno;
   Int_t           evtno;
   Int_t           lumiblock;
   Int_t           isRealData;
   Float_t         muon_d0;
   Bool_t          firedHLT_PFHT475_v2;
   Bool_t          firedHLT_PFHT800_v2;
   Bool_t          passed_CSCTightHaloFilter;
   Bool_t          passed_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          passed_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          passed_goodVertices;
   Bool_t          passed_eeBadScFilter;
   Bool_t          passed_METFilters;

   // List of branches
   TBranch        *b_NJets;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_JetPx;   //!
   TBranch        *b_JetPy;   //!
   TBranch        *b_JetPz;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_EleE;   //!
   TBranch        *b_ElePx;   //!
   TBranch        *b_ElePy;   //!
   TBranch        *b_ElePz;   //!
   TBranch        *b_ElePt;   //!
   TBranch        *b_EleEt;   //!
   TBranch        *b_EleEta;   //!
   TBranch        *b_ElePhi;   //!
   TBranch        *b_PhE;   //!
   TBranch        *b_PhPx;   //!
   TBranch        *b_PhPy;   //!
   TBranch        *b_PhPz;   //!
   TBranch        *b_PhPt;   //!
   TBranch        *b_PhEt;   //!
   TBranch        *b_PhEta;   //!
   TBranch        *b_PhPhi;   //!
   TBranch        *b_MuE;   //!
   TBranch        *b_MuPx;   //!
   TBranch        *b_MuPy;   //!
   TBranch        *b_MuPz;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_ST;   //!
   TBranch        *b_mBH;   //!
   TBranch        *b_Met;   //!
   TBranch        *b_MetPx;   //!
   TBranch        *b_MetPy;   //!
   TBranch        *b_MetPhi;   //!
   TBranch        *b_Sphericity;   //!
   TBranch        *b_JetArr;   //!
   TBranch        *b_EleArr;   //!
   TBranch        *b_MuArr;   //!
   TBranch        *b_PhArr;   //!
   TBranch        *b_NPV;   //!
   TBranch        *b_NTracks;   //!
   TBranch        *b_NElectrons;   //!
   TBranch        *b_NPhotons;   //!
   TBranch        *b_NMuons;   //!
   TBranch        *b_NoScrap;   //!
   TBranch        *b_Multiplicity;   //!
   TBranch        *b_isLeptonPhoton;   //!
   TBranch        *b_isEleChannel;   //!
   TBranch        *b_isMuChannel;   //!
   TBranch        *b_isPhChannel;   //!
   TBranch        *b_ResJetEta;   //!
   TBranch        *b_ResJetPhi;   //!
   TBranch        *b_ResJetM;   //!
   TBranch        *b_ResJetPt;   //!
   TBranch        *b_ResJetEt;   //!
   TBranch        *b_ResEleEta;   //!
   TBranch        *b_ResElePhi;   //!
   TBranch        *b_ResEleM;   //!
   TBranch        *b_ResElePt;   //!
   TBranch        *b_ResPhEta;   //!
   TBranch        *b_ResPhPhi;   //!
   TBranch        *b_ResPhM;   //!
   TBranch        *b_ResPhPt;   //!
   TBranch        *b_ResMuEta;   //!
   TBranch        *b_ResMuPhi;   //!
   TBranch        *b_ResMuM;   //!
   TBranch        *b_ResMuPt;   //!
   TBranch        *b_ResLepEta;   //!
   TBranch        *b_ResLepPhi;   //!
   TBranch        *b_ResLepM;   //!
   TBranch        *b_ResLepPt;   //!
   TBranch        *b_ResObjEta;   //!
   TBranch        *b_ResObjPhi;   //!
   TBranch        *b_ResObjM;   //!
   TBranch        *b_ResObjPt;   //!
   TBranch        *b_runno;   //!
   TBranch        *b_evtno;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_isRealData;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_firedHLT_PFHT475_v2;   //!
   TBranch        *b_firedHLT_PFHT800_v2;   //!
   TBranch        *b_passed_CSCTightHaloFilter;   //!
   TBranch        *b_passed_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_passed_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_passed_goodVertices;   //!
   TBranch        *b_passed_eeBadScFilter;   //!
   TBranch        *b_passed_METFilters;   //!

   t(TTree *tree=0);
   virtual ~t();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef t_cxx
t::t(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuple_output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ntuple_output.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ntuple_output.root:/bhana");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

t::~t()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t t::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t t::LoadTree(Long64_t entry)
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

void t::Init(TTree *tree)
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

   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("JetPx", JetPx, &b_JetPx);
   fChain->SetBranchAddress("JetPy", JetPy, &b_JetPy);
   fChain->SetBranchAddress("JetPz", JetPz, &b_JetPz);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEt", JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("EleE", EleE, &b_EleE);
   fChain->SetBranchAddress("ElePx", ElePx, &b_ElePx);
   fChain->SetBranchAddress("ElePy", ElePy, &b_ElePy);
   fChain->SetBranchAddress("ElePz", ElePz, &b_ElePz);
   fChain->SetBranchAddress("ElePt", ElePt, &b_ElePt);
   fChain->SetBranchAddress("EleEt", EleEt, &b_EleEt);
   fChain->SetBranchAddress("EleEta", EleEta, &b_EleEta);
   fChain->SetBranchAddress("ElePhi", ElePhi, &b_ElePhi);
   fChain->SetBranchAddress("PhE", PhE, &b_PhE);
   fChain->SetBranchAddress("PhPx", PhPx, &b_PhPx);
   fChain->SetBranchAddress("PhPy", PhPy, &b_PhPy);
   fChain->SetBranchAddress("PhPz", PhPz, &b_PhPz);
   fChain->SetBranchAddress("PhPt", PhPt, &b_PhPt);
   fChain->SetBranchAddress("PhEt", PhEt, &b_PhEt);
   fChain->SetBranchAddress("PhEta", PhEta, &b_PhEta);
   fChain->SetBranchAddress("PhPhi", PhPhi, &b_PhPhi);
   fChain->SetBranchAddress("MuE", MuE, &b_MuE);
   fChain->SetBranchAddress("MuPx", MuPx, &b_MuPx);
   fChain->SetBranchAddress("MuPy", MuPy, &b_MuPy);
   fChain->SetBranchAddress("MuPz", MuPz, &b_MuPz);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEt", MuEt, &b_MuEt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("ST", &ST, &b_ST);
   fChain->SetBranchAddress("mBH", &mBH, &b_mBH);
   fChain->SetBranchAddress("Met", &Met, &b_Met);
   fChain->SetBranchAddress("MetPx", &MetPx, &b_MetPx);
   fChain->SetBranchAddress("MetPy", &MetPy, &b_MetPy);
   fChain->SetBranchAddress("MetPhi", &MetPhi, &b_MetPhi);
   fChain->SetBranchAddress("Sphericity", &Sphericity, &b_Sphericity);
   fChain->SetBranchAddress("Jet", Jet, &b_JetArr);
   fChain->SetBranchAddress("Ele", Ele, &b_EleArr);
   fChain->SetBranchAddress("Mu", Mu, &b_MuArr);
   fChain->SetBranchAddress("Ph", Ph, &b_PhArr);
   fChain->SetBranchAddress("NPV", &NPV, &b_NPV);
   fChain->SetBranchAddress("NTracks", &NTracks, &b_NTracks);
   fChain->SetBranchAddress("NElectrons", &NElectrons, &b_NElectrons);
   fChain->SetBranchAddress("NPhotons", &NPhotons, &b_NPhotons);
   fChain->SetBranchAddress("NMuons", &NMuons, &b_NMuons);
   fChain->SetBranchAddress("NoScrap", &NoScrap, &b_NoScrap);
   fChain->SetBranchAddress("Multiplicity", &Multiplicity, &b_Multiplicity);
   fChain->SetBranchAddress("isLeptonPhoton", &isLeptonPhoton, &b_isLeptonPhoton);
   fChain->SetBranchAddress("isEleChannel", &isEleChannel, &b_isEleChannel);
   fChain->SetBranchAddress("isMuChannel", &isMuChannel, &b_isMuChannel);
   fChain->SetBranchAddress("isPhChannel", &isPhChannel, &b_isPhChannel);
   fChain->SetBranchAddress("ResJetEta", &ResJetEta, &b_ResJetEta);
   fChain->SetBranchAddress("ResJetPhi", &ResJetPhi, &b_ResJetPhi);
   fChain->SetBranchAddress("ResJetM", &ResJetM, &b_ResJetM);
   fChain->SetBranchAddress("ResJetPt", &ResJetPt, &b_ResJetPt);
   fChain->SetBranchAddress("ResJetEt", &ResJetEt, &b_ResJetEt);
   fChain->SetBranchAddress("ResEleEta", &ResEleEta, &b_ResEleEta);
   fChain->SetBranchAddress("ResElePhi", &ResElePhi, &b_ResElePhi);
   fChain->SetBranchAddress("ResEleM", &ResEleM, &b_ResEleM);
   fChain->SetBranchAddress("ResElePt", &ResElePt, &b_ResElePt);
   fChain->SetBranchAddress("ResPhEta", &ResPhEta, &b_ResPhEta);
   fChain->SetBranchAddress("ResPhPhi", &ResPhPhi, &b_ResPhPhi);
   fChain->SetBranchAddress("ResPhM", &ResPhM, &b_ResPhM);
   fChain->SetBranchAddress("ResPhPt", &ResPhPt, &b_ResPhPt);
   fChain->SetBranchAddress("ResMuEta", &ResMuEta, &b_ResMuEta);
   fChain->SetBranchAddress("ResMuPhi", &ResMuPhi, &b_ResMuPhi);
   fChain->SetBranchAddress("ResMuM", &ResMuM, &b_ResMuM);
   fChain->SetBranchAddress("ResMuPt", &ResMuPt, &b_ResMuPt);
   fChain->SetBranchAddress("ResLepEta", &ResLepEta, &b_ResLepEta);
   fChain->SetBranchAddress("ResLepPhi", &ResLepPhi, &b_ResLepPhi);
   fChain->SetBranchAddress("ResLepM", &ResLepM, &b_ResLepM);
   fChain->SetBranchAddress("ResLepPt", &ResLepPt, &b_ResLepPt);
   fChain->SetBranchAddress("ResObjEta", &ResObjEta, &b_ResObjEta);
   fChain->SetBranchAddress("ResObjPhi", &ResObjPhi, &b_ResObjPhi);
   fChain->SetBranchAddress("ResObjM", &ResObjM, &b_ResObjM);
   fChain->SetBranchAddress("ResObjPt", &ResObjPt, &b_ResObjPt);
   fChain->SetBranchAddress("runno", &runno, &b_runno);
   fChain->SetBranchAddress("evtno", &evtno, &b_evtno);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("muon_d0", &muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("firedHLT_PFHT475_v2", &firedHLT_PFHT475_v2, &b_firedHLT_PFHT475_v2);
   fChain->SetBranchAddress("firedHLT_PFHT800_v2", &firedHLT_PFHT800_v2, &b_firedHLT_PFHT800_v2);
   fChain->SetBranchAddress("passed_CSCTightHaloFilter", &passed_CSCTightHaloFilter, &b_passed_CSCTightHaloFilter);
   fChain->SetBranchAddress("passed_EcalDeadCellTriggerPrimitiveFilter", &passed_EcalDeadCellTriggerPrimitiveFilter, &b_passed_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("passed_EcalDeadCellBoundaryEnergyFilter", &passed_EcalDeadCellBoundaryEnergyFilter, &b_passed_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("passed_goodVertices", &passed_goodVertices, &b_passed_goodVertices);
   fChain->SetBranchAddress("passed_eeBadScFilter", &passed_eeBadScFilter, &b_passed_eeBadScFilter);
   fChain->SetBranchAddress("passed_METFilters", &passed_METFilters, &b_passed_METFilters);
   Notify();
}

Bool_t t::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void t::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t t::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef t_cxx
