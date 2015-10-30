#include <stdio.h>
#include <iostream>
#include "Riostream.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include <fstream>
#include <string>

#include <TROOT.h>

void lookForBigSTevts(std::string inFilename, std::string outFilename);
float dR(float eta1, float phi1, float eta2, float phi2);

void lookForBigSTevts(std::string inFilename, std::string outFilename) {
  bool debugFlag = true;

  // define output textfile
  ofstream outFile;
  outFile.open(outFilename.c_str());
  // variables calculated in the loop
  float ST          = 0.   ;
  bool passIso      = true ;
  char *messageBuffer = new char[100];

  // variables accessed from the tree
  Char_t firedHLT_PFHT800_v2       ;
  Char_t passed_CSCTightHaloFilter ;
  Char_t passed_goodVertices       ;
  Char_t passed_eeBadScFilter      ;
  int    runno                     ;
  int    evtno                     ;
  int    lumiblock                 ;
  float  JetPt [15]                ;
  float  JetEta[15]                ;
  float  JetPhi[15]                ;
  float  ElePt[25]                 ;
  float  EleEta[25]                ;
  float  ElePhi[25]                ;
  float  PhPt[25]                  ;
  float  PhEta[25]                 ;
  float  PhPhi[25]                 ;
  float  MuPt[25]                  ;
  float  MuEta[25]                 ;
  float  MuPhi[25]                 ;
  float  MetPt                     ; 

  // tree branches
  TBranch  *b_firedHLT_PFHT800_v2       ;
  TBranch  *b_passed_CSCTightHaloFilter ;
  TBranch  *b_passed_goodVertices       ;
  TBranch  *b_passed_eeBadScFilter      ;
  TBranch  *b_JetPt    ;
  TBranch  *b_JetEta   ;
  TBranch  *b_JetPhi   ;
  TBranch  *b_ElePt    ;
  TBranch  *b_EleEta   ;
  TBranch  *b_ElePhi   ;
  TBranch  *b_PhPt     ;
  TBranch  *b_PhEta    ;
  TBranch  *b_PhPhi    ;
  TBranch  *b_MuPt     ;
  TBranch  *b_MuEta    ;
  TBranch  *b_MuPhi    ;
  TBranch  *b_MetPt    ;
  TBranch  *b_runno    ;
  TBranch  *b_evtno    ;
  TBranch  *b_lumiblock;

  //create a chain by looping over the input filename
  TChain chain("bhana/t");
  ifstream infile;
  infile.open(inFilename.c_str()); 
  std::string buffer;
  const char *eosURL = "root://eoscms.cern.ch/";
  chain.SetMakeClass(1);
  while (std::getline(infile, buffer)) {
    std::string ntupleURL = eosURL + buffer; 
    chain.Add(ntupleURL.c_str());
  }

  cout << "Opened chain: " << chain.GetName() << endl;

  // set all branch addresses
  chain.SetBranchAddress( "firedHLT_PFHT800_v2"       , &firedHLT_PFHT800_v2       , &b_firedHLT_PFHT800_v2       );
  chain.SetBranchAddress("passed_CSCTightHaloFilter"  , &passed_CSCTightHaloFilter , &b_passed_CSCTightHaloFilter );
  chain.SetBranchAddress("passed_goodVertices"        , &passed_goodVertices       , &b_passed_goodVertices       );
  chain.SetBranchAddress("passed_eeBadScFilter"       , &passed_eeBadScFilter      , &b_passed_eeBadScFilter      );
  chain.SetBranchAddress( "runno",      &runno,     &b_runno  );
  chain.SetBranchAddress( "evtno",      &evtno,     &b_evtno  );
  chain.SetBranchAddress( "lumiblock",  &lumiblock, &b_lumiblock  );
  chain.SetBranchAddress( "JetPt",      JetPt,      &b_JetPt  );
  chain.SetBranchAddress( "JetEta",     JetEta,     &b_JetEta );
  chain.SetBranchAddress( "JetPhi",     JetPhi,     &b_JetPhi );
  chain.SetBranchAddress( "ElePt",      ElePt,      &b_ElePt  );
  chain.SetBranchAddress( "EleEta",     EleEta,     &b_EleEta );
  chain.SetBranchAddress( "ElePhi",     ElePhi,     &b_ElePhi );
  chain.SetBranchAddress( "PhPt",       PhPt,       &b_PhPt   );
  chain.SetBranchAddress( "PhEta",      PhEta,      &b_PhEta  );
  chain.SetBranchAddress( "PhPhi",      PhPhi,      &b_PhPhi  );
  chain.SetBranchAddress( "MuPt",       MuPt,       &b_MuPt   );
  chain.SetBranchAddress( "MuEta",      MuEta,      &b_MuEta  );
  chain.SetBranchAddress( "MuPhi",      MuPhi,      &b_MuPhi  );
  chain.SetBranchAddress( "MetPt",      &MetPt,     &b_MetPt  );

  const int nEvents = chain.GetEntries();
  cout << "Number of events in chain is: " << nEvents << endl;

  // loop over all events
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // reset variables
    if (iEvent%50000==0) cout << "Scanned " << iEvent << " events." << endl;
    ST      = 0.   ;
    passIso = true ;

    chain.GetEntry(iEvent);
    if (debugFlag) cout << "firedHLT_PFHT800_v2 is: " << firedHLT_PFHT800_v2 << endl;
    // apply trigger and filter requirements
    if (    !firedHLT_PFHT800_v2 || !passed_CSCTightHaloFilter 
        || !passed_goodVertices || !passed_eeBadScFilter      ) continue;

    // apply isolation requirement and calculate ST.
    //cout << "For event number " << iEvent << endl;
    for (int iJet = 0; iJet < 15; ++iJet) {
      if (JetPt[iJet]>20.) {
        for (int iElectron = 0; iElectron < 25; ++iElectron ) {
          if (dR(JetEta[iJet],JetPhi[iJet], EleEta[iElectron], ElePhi[iElectron]) < 0.3) {
            passIso = false;
            break;
          }
        }

        if (!passIso) continue;
        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (dR(JetEta[iJet],JetPhi[iJet], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (dR(JetEta[iJet],JetPhi[iJet], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    JetPt for jet number " << iJet << " is: " << JetPt[iJet] << endl;
        ST += JetPt[iJet];
      }
      else break;
    }
    for (int iElectron = 0; iElectron < 25; ++iElectron) {
      if (ElePt[iElectron]>20.) {
        for (int iJet = 0; iJet < 15; ++iJet ) {
          if (dR(EleEta[iElectron],ElePhi[iElectron], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (dR(EleEta[iElectron],ElePhi[iElectron], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (dR(EleEta[iElectron],ElePhi[iElectron], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    ElePt for electron number " << iElectron << " is: " << ElePt[iElectron] << endl;
        ST += ElePt[iElectron];
      }
      else break;
    }
    for (int iPhoton = 0; iPhoton < 25; ++iPhoton) {
      if (PhPt[iPhoton]>20.) {
        for (int iJet = 0; iJet < 15; ++iJet ) {
          if (dR(PhEta[iPhoton],PhPhi[iPhoton], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iElectron = 0; iElectron < 25; ++iElectron ) {
          if (dR(PhEta[iPhoton], PhPhi[iPhoton], EleEta[iElectron],ElePhi[iElectron]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (dR(PhEta[iPhoton], PhPhi[iPhoton], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    PhPt for photon number " << iPhoton << " is: " << PhPt[iPhoton] << endl;
        ST += PhPt[iPhoton];
      }
      else break;
    }
    for (int iMuon = 0; iMuon < 25; ++iMuon) {
      if (MuPt[iMuon]>20.) {
        for (int iJet = 0; iJet < 15; ++iJet ) {
          if (dR(MuEta[iMuon],MuPhi[iMuon], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iElectron = 0; iElectron < 25; ++iElectron ) {
          if (dR(MuEta[iMuon], MuPhi[iMuon], EleEta[iElectron],ElePhi[iElectron]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (dR( MuEta[iMuon], MuPhi[iMuon], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    MuPt for muon number " << iMuon << " is: " << MuPt[iMuon] << endl;
        ST += MuPt[iMuon];
      }
      else break;
    }
    //
    if (debugFlag) cout << "    MetPt is: " << MetPt << endl;
    ST += MetPt;
    if (debugFlag) cout << "In run number " << runno << " lumi section " << lumiblock << " event number " << evtno << " ST is:" << ST << endl;
    if (ST>5000) {
      sprintf(messageBuffer, "In run number %d lumi section %d event number %d ST is %f\n", runno, lumiblock, evtno, ST);
      outFile << messageBuffer;
    }
    //cout << "    mult is: " << mult << endl;
    //if (ST>9000) {
    //  cout << "Quit working since we found a monster event!" << endl;
    //  break;
    //}
  }
  // write output file
  outFile.close();
}

// function to calculate dR between two objects
float dR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt( ( eta1 - eta2 )*( eta1 - eta2 ) + ( phi1 - phi2 )*( phi1 - phi2 ) );
}
