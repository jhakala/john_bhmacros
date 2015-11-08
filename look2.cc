#include <stdio.h>
#include <iostream>
#include "Riostream.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include <fstream>
#include <string>

#include <TROOT.h>

void look2(std::string inFilename, std::string outFilename);
float dR(float eta1, float phi1, float eta2, float phi2);

void look2(std::string inFilename, std::string outFilename) {
  bool debugFlag = false;

  // define output textfile
  ofstream outFile;
  outFile.open(outFilename.c_str());
  // variables calculated in the loop
  float OurMet      = 0.   ;
  float Px          = 0.   ;
  float Py          = 0.   ;
  float ST          = 0.   ;
  int multiplicity  = 0    ;
  bool passIso      = true ;
  char *messageBuffer = new char[100];

  // variables accessed from the tree
  Bool_t firedHLT_PFHT800_v2       ;
  Bool_t passed_CSCTightHaloFilter ;
  Bool_t passed_goodVertices       ;
  Bool_t passed_eeBadScFilter      ;
  int    runno                     ;
  int    evtno                     ;
  int    lumiblock                 ;
  float  JetEt [25]                ;
  float  JetPx [25]                ;
  float  JetPy [25]                ;
  float  JetEta[25]                ;
  float  JetPhi[25]                ;
  float  EleEt[25]                 ;
  float  ElePx[25]                 ;
  float  ElePy[25]                 ;
  float  EleEta[25]                ;
  float  ElePhi[25]                ;
  float  PhEt[25]                  ;
  float  PhPx[25]                  ;
  float  PhPy[25]                  ;
  float  PhEta[25]                 ;
  float  PhPhi[25]                 ;
  float  MuEt[25]                  ;
  float  MuPx[25]                  ;
  float  MuPy[25]                  ;
  float  MuEta[25]                 ;
  float  MuPhi[25]                 ;
  float  Met                       ; 

  // tree branches
  TBranch  *b_firedHLT_PFHT800_v2       ;
  TBranch  *b_passed_CSCTightHaloFilter ;
  TBranch  *b_passed_goodVertices       ;
  TBranch  *b_passed_eeBadScFilter      ;
  TBranch  *b_JetEt    ;
  TBranch  *b_JetPx    ;
  TBranch  *b_JetPy    ;
  TBranch  *b_JetEta   ;
  TBranch  *b_JetPhi   ;
  TBranch  *b_EleEt    ;
  TBranch  *b_ElePx    ;
  TBranch  *b_ElePy    ;
  TBranch  *b_EleEta   ;
  TBranch  *b_ElePhi   ;
  TBranch  *b_PhEt     ;
  TBranch  *b_PhPx    ;
  TBranch  *b_PhPy    ;
  TBranch  *b_PhEta    ;
  TBranch  *b_PhPhi    ;
  TBranch  *b_MuEt     ;
  TBranch  *b_MuPx    ;
  TBranch  *b_MuPy    ;
  TBranch  *b_MuEta    ;
  TBranch  *b_MuPhi    ;
  TBranch  *b_Met    ;
  TBranch  *b_runno    ;
  TBranch  *b_evtno    ;
  TBranch  *b_lumiblock;

  //create a chain by looping over the input filename
  TChain chain("bhana/t");
  //ifstream infile;
  //infile.open(inFilename.c_str()); 
  //std::string buffer;
  //const char *eosURL = "root://eoscms.cern.ch/";
  //chain.SetMakeClass(1);
  //while (std::getline(infile, buffer)) {
  //  std::string ntupleURL = eosURL + buffer; 
    chain.Add(inFilename.c_str());
  //}

  cout << "Opened chain: " << chain.GetName() << endl;

  // set all branch addresses
  chain.SetBranchAddress( "firedHLT_PFHT800_v2"       , &firedHLT_PFHT800_v2       , &b_firedHLT_PFHT800_v2       );
  chain.SetBranchAddress("passed_CSCTightHaloFilter"  , &passed_CSCTightHaloFilter , &b_passed_CSCTightHaloFilter );
  chain.SetBranchAddress("passed_goodVertices"        , &passed_goodVertices       , &b_passed_goodVertices       );
  chain.SetBranchAddress("passed_eeBadScFilter"       , &passed_eeBadScFilter      , &b_passed_eeBadScFilter      );
  chain.SetBranchAddress( "runno",      &runno,     &b_runno  );
  chain.SetBranchAddress( "evtno",      &evtno,     &b_evtno  );
  chain.SetBranchAddress( "lumiblock",  &lumiblock, &b_lumiblock  );
  chain.SetBranchAddress( "JetEt",      JetEt,      &b_JetEt  );
  chain.SetBranchAddress( "JetPx",      JetPx,      &b_JetPx  );
  chain.SetBranchAddress( "JetPy",      JetPy,      &b_JetPy  );
  chain.SetBranchAddress( "JetEta",     JetEta,     &b_JetEta );
  chain.SetBranchAddress( "JetPhi",     JetPhi,     &b_JetPhi );
  chain.SetBranchAddress( "EleEt",      EleEt,      &b_EleEt  );
  chain.SetBranchAddress( "ElePx",      ElePx,      &b_ElePx  );
  chain.SetBranchAddress( "ElePy",      ElePy,      &b_ElePy  );
  chain.SetBranchAddress( "EleEta",     EleEta,     &b_EleEta );
  chain.SetBranchAddress( "ElePhi",     ElePhi,     &b_ElePhi );
  chain.SetBranchAddress( "PhEt",       PhEt,       &b_PhEt   );
  chain.SetBranchAddress( "PhPx",       PhPx,       &b_PhPx  );
  chain.SetBranchAddress( "PhPy",       PhPy,       &b_PhPy  );
  chain.SetBranchAddress( "PhEta",      PhEta,      &b_PhEta  );
  chain.SetBranchAddress( "PhPhi",      PhPhi,      &b_PhPhi  );
  chain.SetBranchAddress( "MuEt",       MuEt,       &b_MuEt   );
  chain.SetBranchAddress( "MuPx",       MuPx,       &b_MuPx  );
  chain.SetBranchAddress( "MuPy",       MuPy,       &b_MuPy  );
  chain.SetBranchAddress( "MuEta",      MuEta,      &b_MuEta  );
  chain.SetBranchAddress( "MuPhi",      MuPhi,      &b_MuPhi  );
  chain.SetBranchAddress( "Met",        &Met,       &b_Met  );

  const int nEvents = chain.GetEntries();
  cout << "Number of events in chain is: " << nEvents << endl;

  // loop over all events
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // reset variables
    if (iEvent%50000==0) cout << "Scanned " << iEvent << " events." << endl;
    OurMet  = 0.   ;
    Px      = 0.   ;
    Py      = 0.   ;
    ST      = 0.   ;
    multiplicity = 0;
    passIso = true ;

    chain.GetEntry(iEvent);
    //  if (runno != 257613 || lumiblock != 614 || evtno != 962690898) continue;
    // apply trigger and filter requirements
    if (    !firedHLT_PFHT800_v2 || !passed_CSCTightHaloFilter 
        || !passed_goodVertices || !passed_eeBadScFilter      ) continue;
    // apply isolation requirement and calculate ST and MET.
    //cout << "For event number " << iEvent << endl;
    for (int iJet = 0; iJet < 25; ++iJet) {
      if (JetEt[iJet]>20.) {
        //for (int iElectron = 0; iElectron < 25; ++iElectron ) {
        //  if (EleEt[iElectron]>20 && dR(JetEta[iJet],JetPhi[iJet], EleEta[iElectron], ElePhi[iElectron]) < 0.3) {
        //    passIso = false;
        //    break;
        //  }
        //}

        //if (!passIso) continue;
        //for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
        //  if (PhEt[iPhoton]>20 && dR(JetEta[iJet],JetPhi[iJet], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
        //    passIso = false;
        //    break;
        //  }
        //}
        //if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>20 && dR(JetEta[iJet],JetPhi[iJet], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
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
      if (EleEt[iElectron]>20.) {
        for (int iJet = 0; iJet < 25; ++iJet ) {
          if (JetEt[iJet]>20 && dR(EleEta[iElectron],ElePhi[iElectron], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (PhEt[iPhoton]>20 && dR(EleEta[iElectron],ElePhi[iElectron], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>20 && dR(EleEta[iElectron],ElePhi[iElectron], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
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
      if (PhEt[iPhoton]>20.) {
        for (int iJet = 0; iJet < 25; ++iJet ) {
          if (JetEt[iJet]>20 && dR(PhEta[iPhoton],PhPhi[iPhoton], JetEta[iJet], JetPhi[iJet]) < 0.3) {
            passIso = false;
            break;
          }
        }
        if (!passIso) continue;

        //for (int iElectron = 0; iElectron < 25; ++iElectron ) {
        //  if (EleEt[iElectron]>20 && dR(PhEta[iPhoton], PhPhi[iPhoton], EleEta[iElectron],ElePhi[iElectron]) < 0.3) {
        //    passIso = false;
        //    break;
        //  }
        //}
        //if (!passIso) continue;

        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>20 && dR(PhEta[iPhoton], PhPhi[iPhoton], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
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
      if (MuEt[iMuon]>20.) {
        //for (int iJet = 0; iJet < 25; ++iJet ) {
        //  if (JetEt[iJet]>20 && dR(MuEta[iMuon],MuPhi[iMuon], JetEta[iJet], JetPhi[iJet]) < 0.3) {
        //    passIso = false;
        //    break;
        //  }
        //}
        //if (!passIso) continue;

        //for (int iElectron = 0; iElectron < 25; ++iElectron ) {
        //  if (EleEt[iElectron]>20 && dR(MuEta[iMuon], MuPhi[iMuon], EleEta[iElectron],ElePhi[iElectron]) < 0.3) {
        //    passIso = false;
        //    break;
        //  }
        //}
        //if (!passIso) continue;

        //for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
        //  if (PhEt[iPhoton]>20 && dR( MuEta[iMuon], MuPhi[iMuon], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
        //    passIso = false;
        //    break;
        //  }
        //}
        //if (!passIso) continue;

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
    //
    if (debugFlag) cout << "    Met from PAT collection is: " << Met << endl;
    OurMet = std::sqrt(Px*Px + Py*Py);
    if (debugFlag) cout << "    Met calculated according to my recipe is: " << OurMet << endl;
    ST += OurMet;
    if (ST>3000 && multiplicity>10) {
      cout << "In run number " << runno << " lumi section " << lumiblock << " event number " << evtno << " ST is:" << ST << endl;
        sprintf(messageBuffer, "In run number %d lumi section %d event number %d ST is %f and multiplicity is %d\n", runno, lumiblock, evtno, ST, multiplicity);
        outFile << messageBuffer;
        if (debugFlag) cout << messageBuffer;
      for (int j=0; j<25; ++j) {
        if (JetEt[j]>50) {
          sprintf(messageBuffer, "    Jet %d has Et=%f, Eta=%f, Phi=%f\n", j, JetEt[j], JetEta[j], JetPhi[j]);
          outFile << messageBuffer;
          if (debugFlag) cout << messageBuffer;
        }
        if (EleEt[j]>50) {
          outFile << messageBuffer;
          if (debugFlag) cout  << messageBuffer;
        } 
        if (PhEt[j]>50) {
          sprintf(messageBuffer, "    Photon %d has Et=%f, Eta=%f, Phi=%f\n", j, PhEt[j], PhEta[j], PhPhi[j]);
          outFile << messageBuffer;
          if (debugFlag) cout  << messageBuffer;
        }
        if (MuEt[j]>50) {
          sprintf(messageBuffer, "    Muon %d has Et=%f, Eta=%f, Phi=%f\n", j, MuEt[j], MuEta[j], MuPhi[j]);
          outFile << messageBuffer;
          if (debugFlag) cout  << messageBuffer;
        } 
      }
      sprintf(messageBuffer, "    our MET is=%f\n", OurMet);
      outFile << messageBuffer;
      if (debugFlag) cout  << messageBuffer;
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
