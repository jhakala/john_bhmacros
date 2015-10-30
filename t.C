#define t_cxx
#include "t.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void t::Loop()
{
//   In a ROOT session, you can do:
//      root> .L t.C
//      root> t t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   bool checkTriggers = false;
   bool checkJets = false;
   bool checkElectrons = false;
   bool checkPhotons = false; 
   bool checkMuons = false;
   bool checkMETs = true;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (checkTriggers) cout << "Event number " << ientry << endl;
      if (checkTriggers) cout << "   firedHLT_PFHT475_v2 is " << (firedHLT_PFHT475_v2 ? "true" : "false") << endl;
      if (checkTriggers) cout << "   firedHLT_PFHT800_v2 is " << (firedHLT_PFHT800_v2 ? "true" : "false") << endl;
      if (checkJets && NJets>=1) cout << "Event number " << ientry << endl << "    NJets  = " << NJets << endl;
      if (checkElectrons && NElectrons>=1) cout << "Event number " << ientry << endl << "    NElectrons  = " << NElectrons << endl;
      if (checkPhotons && NPhotons>=1) cout << "Event number " << ientry << endl << "    NPhotons  = " << NPhotons << endl;
      if (checkMuons && NMuons>=1) cout << "Event number " << ientry << endl << "    NMuons  = " << NMuons << endl;

      if (checkJets && NJets>=1) {
        for (int i = 0; i<sizeof(JetPt)/sizeof(JetPt[0]); ++i) {
          cout << "    JetE[" << i << "] = "<< JetE[i] << endl;
          cout << "    JetEt[" << i << "] = "<< JetEt[i] << endl;
          cout << "    JetPt[" << i << "] = "<< JetPt[i] << endl;
          cout << "    JetPx[" << i << "] = "<< JetPx[i] << endl;
          cout << "    JetPy[" << i << "] = "<< JetPy[i] << endl;
          cout << "    JetPz[" << i << "] = "<< JetPz[i] << endl;
          cout << "    JetEta[" << i << "] = "<< JetEta[i] << endl;
          cout << "    JetPhi[" << i << "] = "<< JetPhi[i] << endl;
        }
      }
      if (checkElectrons && NElectrons>=1) {
        for (int i = 0; i<sizeof(ElePt)/sizeof(ElePt[0]); ++i) {
          cout << "    EleE[" << i << "] = "<< EleE[i] << endl;
          cout << "    EleEt[" << i << "] = "<< EleEt[i] << endl;
          cout << "    ElePt[" << i << "] = "<< ElePt[i] << endl;
          cout << "    ElePx[" << i << "] = "<< ElePx[i] << endl;
          cout << "    ElePy[" << i << "] = "<< ElePy[i] << endl;
          cout << "    ElePz[" << i << "] = "<< ElePz[i] << endl;
          cout << "    EleEta[" << i << "] = "<< EleEta[i] << endl;
          cout << "    ElePhi[" << i << "] = "<< ElePhi[i] << endl;
        }
      }
      if (checkPhotons && NPhotons>=1) {
        for (int i = 0; i<sizeof(PhPt)/sizeof(PhPt[0]); ++i) {
          cout << "    PhE[" << i << "] = "<< PhE[i] << endl;
          cout << "    PhEt[" << i << "] = "<< PhEt[i] << endl;
          cout << "    PhPt[" << i << "] = "<< PhPt[i] << endl;
          cout << "    PhPx[" << i << "] = "<< PhPx[i] << endl;
          cout << "    PhPy[" << i << "] = "<< PhPy[i] << endl;
          cout << "    PhPz[" << i << "] = "<< PhPz[i] << endl;
          cout << "    PhEta[" << i << "] = "<< PhEta[i] << endl;
          cout << "    PhPhi[" << i << "] = "<< PhPhi[i] << endl;
        }
      }
      if (checkMuons && NMuons>=1) {
        for (int i = 0; i<sizeof(MuPt)/sizeof(MuPt[0]); ++i) {
          cout << "    MuE[" << i << "] = "<< MuE[i] << endl;
          cout << "    MuEt[" << i << "] = "<< MuEt[i] << endl;
          cout << "    MuPt[" << i << "] = "<< MuPt[i] << endl;
          cout << "    MuPx[" << i << "] = "<< MuPx[i] << endl;
          cout << "    MuPy[" << i << "] = "<< MuPy[i] << endl;
          cout << "    MuPz[" << i << "] = "<< MuPz[i] << endl;
          cout << "    MuEta[" << i << "] = "<< MuEta[i] << endl;
          cout << "    MuPhi[" << i << "] = "<< MuPhi[i] << endl;
        }
      }
      if (checkMETs) {
        cout << "    Met = "    << Met << endl;
        cout << "    MetPx = "  << MetPx << endl;
        cout << "    MetPy = "  << MetPy << endl;
        cout << "    MetPhi = " << MetPhi << endl;
      }
   }
}
