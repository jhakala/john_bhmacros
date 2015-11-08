#include <stdio.h>
#include <iostream>
#include "Riostream.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include <fstream>
#include <string>

#include <TROOT.h>

void SThist_data(std::string inFilename, std::string outFilename);
float dR(float eta1, float phi1, float eta2, float phi2);

void SThist_data(std::string inFilename, std::string outFilename) {
  bool debugFlag=false;
  int debugMaxEvents = 10;
  bool runOnEosFiles=false;

  // define output file and output histogram
  TFile *outfile = new TFile(outFilename.c_str(),"RECREATE");
  TH1F stHist = TH1F("stHist", "ST", 100, 700, 9700);
  int multMax = 12;

  // loop to create ST histograms for inclusive and exclusive multiplicities from 2 up to multMax
  TH1F *stIncHist[multMax-2];
  TH1F *stExcHist[multMax-2];
  char *histTitle = new char[11];
  int mult=2;
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    sprintf(histTitle, "stInc%02dHist", mult);
    stIncHist[iHist] = new TH1F(histTitle, "Inclusive ST", 100, 700, 9700);
    sprintf(histTitle, "stExc%02dHist", mult);
    stExcHist[iHist] = new TH1F(histTitle, "Exclusive ST", 100, 700, 9700);
    ++mult;
  }

  // variables calculated in the loop over events
  float OurMet      = 0.   ;
  float Px          = 0.   ;
  float Py          = 0.   ;
  float ST          = 0.   ;
  int multiplicity  = 0    ;
  bool passIso      = true ;
  int nPassedEvents = 0    ;

  // variables accessed from the tree
  Bool_t   firedHLT_PFHT800_v2       ;
  Bool_t   passed_CSCTightHaloFilter ;
  Bool_t   passed_goodVertices       ;
  Bool_t   passed_eeBadScFilter      ;
  int      runno                     ;
  int      evtno                     ;
  int      lumiblock                 ;
  float    JetEt [25]                ;
  float    JetPx [25]                ;
  float    JetPy [25]                ;
  Float_t  JetEta[25]                ;
  Float_t  JetPhi[25]                ;
  float    EleEt[25]                 ;
  float    ElePx[25]                 ;
  float    ElePy[25]                 ;
  Float_t  EleEta[25]                ;
  Float_t  ElePhi[25]                ;
  float    PhEt[25]                  ;
  float    PhPx[25]                  ;
  float    PhPy[25]                  ;
  Float_t  PhEta[25]                 ;
  Float_t  PhPhi[25]                 ;
  float    MuEt[25]                  ;
  float    MuPx[25]                  ;
  float    MuPy[25]                  ;
  Float_t  MuEta[25]                 ;
  Float_t  MuPhi[25]                 ;
  Float_t  Met                       ; 

  // tree branches
  TBranch  *b_firedHLT_PFHT800_v2       ;
  TBranch  *b_passed_CSCTightHaloFilter ;
  TBranch  *b_passed_goodVertices       ;
  TBranch  *b_passed_eeBadScFilter      ;
  TBranch  *b_runno                     ;
  TBranch  *b_evtno                     ;
  TBranch  *b_lumiblock                 ;
  TBranch  *b_JetEt                     ;
  TBranch  *b_JetPx                     ;
  TBranch  *b_JetPy                     ;
  TBranch  *b_JetEta                    ;
  TBranch  *b_JetPhi                    ;
  TBranch  *b_EleEt                     ;
  TBranch  *b_ElePx                     ;
  TBranch  *b_ElePy                     ;
  TBranch  *b_EleEta                    ;
  TBranch  *b_ElePhi                    ;
  TBranch  *b_PhEt                      ;
  TBranch  *b_PhPx                      ;
  TBranch  *b_PhPy                      ;
  TBranch  *b_PhEta                     ;
  TBranch  *b_PhPhi                     ;
  TBranch  *b_MuEt                      ;
  TBranch  *b_MuPx                      ;
  TBranch  *b_MuPy                      ;
  TBranch  *b_MuEta                     ;
  TBranch  *b_MuPhi                     ;
  TBranch  *b_Met                       ;

  //create a chain by looping over the input filename
  TChain chain("bhana/t");
  if (runOnEosFiles) {
    ifstream infile;
    infile.open(inFilename.c_str()); 
    std::string buffer;
    const char *eosURL = "root://eoscms.cern.ch/";
    //chain.SetMakeClass(1);
    while (std::getline(infile, buffer)) {
      std::string ntupleURL = eosURL + buffer; 
      chain.Add(ntupleURL.c_str());
    }
  }
  else chain.Add(inFilename.c_str());

  cout << "Opened chain: " << chain.GetName() << endl;

  // set all branch addresses
  chain.SetBranchAddress( "firedHLT_PFHT800_v2"       , &firedHLT_PFHT800_v2       , &b_firedHLT_PFHT800_v2       );
  chain.SetBranchAddress( "passed_CSCTightHaloFilter" , &passed_CSCTightHaloFilter , &b_passed_CSCTightHaloFilter );
  chain.SetBranchAddress( "passed_goodVertices"       , &passed_goodVertices       , &b_passed_goodVertices       );
  chain.SetBranchAddress( "passed_eeBadScFilter"      , &passed_eeBadScFilter      , &b_passed_eeBadScFilter      );
  chain.SetBranchAddress( "runno"                     , &runno                     , &b_runno                     );
  chain.SetBranchAddress( "evtno"                     , &evtno                     , &b_evtno                     );
  chain.SetBranchAddress( "lumiblock"                 , &lumiblock                 , &b_lumiblock                 );
  chain.SetBranchAddress( "JetEt"                     , JetEt                      , &b_JetEt                     );
  chain.SetBranchAddress( "JetPx"                     , JetPx                      , &b_JetPx                     );
  chain.SetBranchAddress( "JetPy"                     , JetPy                      , &b_JetPy                     );
  chain.SetBranchAddress( "JetEta"                    , JetEta                     , &b_JetEta                    );
  chain.SetBranchAddress( "JetPhi"                    , JetPhi                     , &b_JetPhi                    );
  chain.SetBranchAddress( "EleEt"                     , EleEt                      , &b_EleEt                     );
  chain.SetBranchAddress( "ElePx"                     , ElePx                      , &b_ElePx                     );
  chain.SetBranchAddress( "ElePy"                     , ElePy                      , &b_ElePy                     );
  chain.SetBranchAddress( "EleEta"                    , EleEta                     , &b_EleEta                    );
  chain.SetBranchAddress( "ElePhi"                    , ElePhi                     , &b_ElePhi                    );
  chain.SetBranchAddress( "PhEt"                      , PhEt                       , &b_PhEt                      );
  chain.SetBranchAddress( "PhPx"                      , PhPx                       , &b_PhPx                      );
  chain.SetBranchAddress( "PhPy"                      , PhPy                       , &b_PhPy                      );
  chain.SetBranchAddress( "PhEta"                     , PhEta                      , &b_PhEta                     );
  chain.SetBranchAddress( "PhPhi"                     , PhPhi                      , &b_PhPhi                     );
  chain.SetBranchAddress( "MuEt"                      , MuEt                       , &b_MuEt                      );
  chain.SetBranchAddress( "MuPx"                      , MuPx                       , &b_MuPx                      );
  chain.SetBranchAddress( "MuPy"                      , MuPy                       , &b_MuPy                      );
  chain.SetBranchAddress( "MuEta"                     , MuEta                      , &b_MuEta                     );
  chain.SetBranchAddress( "MuPhi"                     , MuPhi                      , &b_MuPhi                     );
  chain.SetBranchAddress( "Met"                       , &Met                       , &b_Met                       );

  const int nEvents = chain.GetEntries();
  cout << "Number of events in chain is: " << nEvents << endl;

  // loop over all events
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // reset variables
    if (iEvent%50000==0) cout << "Processed " << iEvent << " events." << endl;
    OurMet       = 0.   ;
    Px           = 0.   ;
    Py           = 0.   ;
    ST           = 0.   ;
    multiplicity = 0    ;
    passIso      = true ;

    chain.GetEntry(iEvent);
    // apply trigger and filter requirements
    if (    !firedHLT_PFHT800_v2 || !passed_CSCTightHaloFilter 
        || !passed_goodVertices || !passed_eeBadScFilter      ) continue;
    // apply isolation requirement and calculate ST and MET.
    if (debugFlag) cout << "For event number " << iEvent << endl;
    for (int iJet = 0; iJet < 25; ++iJet) {
      if (JetEt[iJet]>50.) {

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>50 && dR(JetEta[iJet],JetPhi[iJet], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    JetEt for jet number " << iJet << " is: " << JetEt[iJet] << endl;
        if (JetEt[iJet] > 50) {
          ST += JetEt[iJet];
          multiplicity+=1;
        }
        Px += JetPx[iJet];
        Py += JetPy[iJet];
      }
      else break;
    }
    for (int iElectron = 0; iElectron < 25; ++iElectron) {
      if (EleEt[iElectron]>50.) {
        for (int iJet = 0; iJet < 25; ++iJet ) {
          if (JetEt[iJet]>50 && dR(EleEta[iElectron],ElePhi[iElectron], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (PhEt[iPhoton]>50 && dR(EleEta[iElectron],ElePhi[iElectron], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>50 && dR(EleEta[iElectron],ElePhi[iElectron], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    EleEt for electron number " << iElectron << " is: " << EleEt[iElectron] << endl;
        if (EleEt[iElectron] > 50) {
          ST += EleEt[iElectron];
          multiplicity+=1;
        }
        Px += ElePx[iElectron];
        Py += ElePy[iElectron];
      }
      else break;
    }
    for (int iPhoton = 0; iPhoton < 25; ++iPhoton) {
      if (PhEt[iPhoton]>50.) {
        for (int iJet = 0; iJet < 25; ++iJet ) {
          if (JetEt[iJet]>50 && dR(PhEta[iPhoton],PhPhi[iPhoton], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;


        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>50 && dR(PhEta[iPhoton], PhPhi[iPhoton], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    PhEt for photon number " << iPhoton << " is: " << PhEt[iPhoton] << endl;
        if (PhEt[iPhoton] > 50) {
          ST += PhEt[iPhoton];
          multiplicity+=1;
        }
        Px += PhPx[iPhoton];
        Py += PhPy[iPhoton];
      }
      else break;
    }
    for (int iMuon = 0; iMuon < 25; ++iMuon) {
      if (MuEt[iMuon]>50.) {
        if (debugFlag) cout << "    MuEt for muon number " << iMuon << " is: " << MuEt[iMuon] << endl;
        if (MuEt[iMuon] > 50) {
          ST += MuEt[iMuon];
          multiplicity+=1;
        }
        Px += MuPx[iMuon];
        Py += MuPy[iMuon];
      }
      else break;
    }
    if (debugFlag) cout << "    Met from PAT collection is: " << Met << endl;
    OurMet = std::sqrt(Px*Px + Py*Py);
    if (debugFlag) cout << "    Met calculated according to my recipe is: " << OurMet << endl;
    ST += OurMet;
    stHist.Fill(ST);
    for (int iHist = 0; iHist<multMax-2; ++iHist) {
      if (multiplicity == iHist+2) stExcHist[iHist]->Fill(ST);
      if (multiplicity >= iHist+2) stIncHist[iHist]->Fill(ST);
    }
    nPassedEvents+=1;
    if (debugFlag && debugMaxEvents==nPassedEvents) break;
  }
  // write the histogram and the output file
  outfile->cd();
  stHist.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHist[iHist]->Write();
    stIncHist[iHist]->Write();
  }
}

// function to calculate dR between two objects
float dR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt( ( eta1 - eta2 )*( eta1 - eta2 ) + ( phi1 - phi2 )*( phi1 - phi2 ) );
}
