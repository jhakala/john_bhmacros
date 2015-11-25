#include <stdio.h>
#include <iostream>
#include "Riostream.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include <string>
#include <TROOT.h>
#include <TMath.h>

void STnewIso(std::string inFilename, std::string outFilename);
float dR(float eta1, float phi1, float eta2, float phi2);

void STnewIso(std::string inFilename, std::string outFilename) {
  bool debugFlag     = true ;
  int  eventsToDump  = 100    ;  // if debugFlag is true, then stop once the number of dumped events reaches eventsToDump
  int  nDumpedEvents = 0     ;

  // define output root file
  TFile* outFile = new TFile(outFilename.c_str(), "RECREATE");

  // define output histograms
  // loop to create ST histograms for inclusive and exclusive multiplicities from 2 up to multMax
  TH1F stHist = TH1F("stHist", "ST", 100, 500, 10500);
  int mult=2;
  int multMax = 12;
  TH1F *stIncHist[multMax-2];
  TH1F *stExcHist[multMax-2];
  TH1F stHistMHT = TH1F("stHistMHT", "ST using MHT", 100, 500, 10500);
  TH1F *stIncHistMHT[multMax-2];
  TH1F *stExcHistMHT[multMax-2];
  char *histTitle = new char[20];
  // These use pat::slimmedMETs
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    sprintf(histTitle, "stInc%02dHist", mult);
    stIncHist[iHist] = new TH1F(histTitle, "Inclusive ST", 100, 500, 10500);
    sprintf(histTitle, "stExc%02dHist", mult);
    stExcHist[iHist] = new TH1F(histTitle, "Exclusive ST", 100, 500, 10500);
    ++mult;
  }
  mult=2;
  // These use MHT
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    sprintf(histTitle, "stInc%02dHistMHT", mult);
    stIncHistMHT[iHist] = new TH1F(histTitle, "Inclusive ST using MHT", 100, 500, 10500);
    sprintf(histTitle, "stExc%02dHistMHT", mult);
    stExcHistMHT[iHist] = new TH1F(histTitle, "Exclusive ST using MHT", 100, 500, 10500);
    ++mult;
  }

  // variables calculated in the loop
  float MHT             = 0.            ;
  float Px              = 0.            ;
  float Py              = 0.            ;
  float ST              = 0.            ;
  float STMHTnoMET      = 0.            ;
  int multiplicity      = 0             ;
  bool passIso          = true          ;
  char *messageBuffer   = new char[400] ;
  bool eventHasMuon     = false         ;
  bool eventHasPhoton   = false         ;
  bool eventHasElectron = false         ;
  float JetMuonEt       = 0.            ;
  float JetElectronEt   = 0.            ;
  float JetPhotonEt     = 0.            ;

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
  TBranch  *b_runno                     ;
  TBranch  *b_evtno                     ;
  TBranch  *b_lumiblock                 ;

  //create a chain using the input filename
  TChain chain("bhana/t");
  chain.Add(inFilename.c_str());

  cout << "Opened chain: " << chain.GetName() << endl;

  // set all branch addresses
  chain.SetBranchAddress( "firedHLT_PFHT800_v2"       , &firedHLT_PFHT800_v2       , &b_firedHLT_PFHT800_v2       );
  chain.SetBranchAddress("passed_CSCTightHaloFilter"  , &passed_CSCTightHaloFilter , &b_passed_CSCTightHaloFilter );
  chain.SetBranchAddress("passed_goodVertices"        , &passed_goodVertices       , &b_passed_goodVertices       );
  chain.SetBranchAddress("passed_eeBadScFilter"       , &passed_eeBadScFilter      , &b_passed_eeBadScFilter      );

  chain.SetBranchAddress( "runno",      &runno,     &b_runno      );
  chain.SetBranchAddress( "lumiblock",  &lumiblock, &b_lumiblock  );
  chain.SetBranchAddress( "evtno",      &evtno,     &b_evtno      );
  chain.SetBranchAddress( "JetEt",      JetEt,      &b_JetEt      );
  chain.SetBranchAddress( "JetPx",      JetPx,      &b_JetPx      );
  chain.SetBranchAddress( "JetPy",      JetPy,      &b_JetPy      );
  chain.SetBranchAddress( "JetEta",     JetEta,     &b_JetEta     );
  chain.SetBranchAddress( "JetPhi",     JetPhi,     &b_JetPhi     );
  chain.SetBranchAddress( "EleEt",      EleEt,      &b_EleEt      );
  chain.SetBranchAddress( "ElePx",      ElePx,      &b_ElePx      );
  chain.SetBranchAddress( "ElePy",      ElePy,      &b_ElePy      );
  chain.SetBranchAddress( "EleEta",     EleEta,     &b_EleEta     );
  chain.SetBranchAddress( "ElePhi",     ElePhi,     &b_ElePhi     );
  chain.SetBranchAddress( "PhEt",       PhEt,       &b_PhEt       );
  chain.SetBranchAddress( "PhPx",       PhPx,       &b_PhPx       );
  chain.SetBranchAddress( "PhPy",       PhPy,       &b_PhPy       );
  chain.SetBranchAddress( "PhEta",      PhEta,      &b_PhEta      );
  chain.SetBranchAddress( "PhPhi",      PhPhi,      &b_PhPhi      );
  chain.SetBranchAddress( "MuEt",       MuEt,       &b_MuEt       );
  chain.SetBranchAddress( "MuPx",       MuPx,       &b_MuPx       );
  chain.SetBranchAddress( "MuPy",       MuPy,       &b_MuPy       );
  chain.SetBranchAddress( "MuEta",      MuEta,      &b_MuEta      );
  chain.SetBranchAddress( "MuPhi",      MuPhi,      &b_MuPhi      );
  chain.SetBranchAddress( "Met",        &Met,       &b_Met        );

  const int nEvents = chain.GetEntries();
  cout << "Number of events in chain is: " << nEvents << endl;

  // loop over all events
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (iEvent%100==0) cout << "Scanned " << iEvent << " events." << endl;

    // reset variables
    ST               = 0.    ;
    multiplicity     = 0     ;
    Px               = 0.    ;
    Py               = 0.    ;
    MHT              = 0.    ;
    eventHasMuon     = false ;
    eventHasPhoton   = false ;
    eventHasElectron = false ;

    chain.GetEntry(iEvent);
    // apply trigger and filter requirements
    if (    !firedHLT_PFHT800_v2 || !passed_CSCTightHaloFilter
        || !passed_goodVertices || !passed_eeBadScFilter      ) continue;

    if (debugFlag) {
      sprintf(messageBuffer, "\n\n--------------------------------------------------------------\n");
      sprintf(messageBuffer, "\n\nEntry number %d -- this passed the trigger/METFilter criteria.\n", iEvent);
      cout << messageBuffer;
    }
    // apply isolation requirement and calculate ST and MHT.
    //Jets
    for (int iJet = 0; iJet < 25; ++iJet) {
      passIso=true;
      JetMuonEt     =0;
      JetElectronEt =0;
      JetPhotonEt   =0;
      if (JetEt[iJet]>50.) {
        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>50) {
            eventHasMuon = true;
            if (JetEt[iJet] && dR(JetEta[iJet],JetPhi[iJet], MuEta[iMuon], MuPhi[iMuon]) < 0.05) {
              JetMuonEt+=MuEt[iMuon];
              if (JetMuonEt>0.8*JetEt[iJet]) {
                passIso = false;
                if (debugFlag) {
                  sprintf(messageBuffer, "Jet number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %d\n", iJet, iMuon, runno, lumiblock, evtno);
                  cout << messageBuffer;
                }
                break;
              }
            }
          }
        }
        for (int iElectron = 0; iElectron < 25; ++iElectron ) {
          if (EleEt[iElectron]>50) {
            eventHasElectron = true;
            if (dR(JetEta[iJet],JetPhi[iJet], EleEta[iElectron], ElePhi[iElectron]) < 0.05) {
              JetElectronEt+=EleEt[iElectron];
              if (JetElectronEt > 0.7*JetEt[iJet] ) {
                passIso = false;
                if (debugFlag) {
                  sprintf(messageBuffer, "Jet number %d failed isolation with Electron number %d  in run number %d lumi section %d event number %d\n", iJet, iElectron, runno, lumiblock, evtno);
                  cout << messageBuffer;
                }
                break;
              }
            }
          }
        }
        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (PhEt[iPhoton]>50) {
            eventHasPhoton = true;
            if (dR(JetEta[iJet],JetPhi[iJet], PhEta[iPhoton], PhPhi[iPhoton]) < 0.05) {
              JetPhotonEt+=PhEt[iPhoton];
              if (JetPhotonEt>0.5*JetEt[iJet] ) {
                passIso = false;
                if (debugFlag) {
                  sprintf(messageBuffer, "Jet number %d failed isolation with Photon number %d  in run number %d lumi section %d event number %d\n", iJet, iPhoton, runno, lumiblock, evtno);
                  cout << messageBuffer;
                }
                break;
              }
            }
          }
        }
        if (!passIso) continue;

        if (debugFlag) cout << "    JetEt for jet number " << iJet << " is: " << JetEt[iJet] << endl;
        ST += JetEt[iJet];
        multiplicity+=1;
        if (debugFlag) {
          sprintf(messageBuffer, "Jet number %d passed isolation in run number %d lumi section %d event number %d.\n       It had Px=%f and Py=%f\n", iJet, runno, lumiblock, evtno, JetPx[iJet], JetPy[iJet]);
          cout << messageBuffer;
        }
        Px += JetPx[iJet];
        Py += JetPy[iJet];
        if (debugFlag) {
          sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
          cout << messageBuffer;
        }
      }
      else break;
    }

    //Electrons
    if (eventHasElectron) {
      for (int iElectron = 0; iElectron < 25; ++iElectron) {
        passIso=true;
        if (EleEt[iElectron]>50.) {
          for (int iJet = 0; iJet < 25; ++iJet ) {
            if (JetEt[iJet]>50 && dR(EleEta[iElectron],ElePhi[iElectron], JetEta[iJet], JetPhi[iJet]) < 0.05) {
              if (EleEt[iElectron]<0.7*JetEt[iJet]) {
                passIso = false;
                if (debugFlag) {
                  sprintf(messageBuffer, "Electron number %d failed isolation with Jet number %d  in run number %d lumi section %d event number %d\n", iElectron, iJet, runno, lumiblock, evtno);
                  cout << messageBuffer;
                }
                break;
              }
            }
          }
          if (!passIso) continue;
         
          // Throw away electron if there's an electron/muon overlap.
          for (int iMuon = 0; iMuon < 25; ++iMuon ) {
            if (MuEt[iMuon]>50 && dR(EleEta[iElectron],ElePhi[iElectron], MuEta[iMuon], MuPhi[iMuon]) < 0.05) {
              passIso = false;
              if (debugFlag) {
                sprintf(messageBuffer, "Electron number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %d\n", iElectron, iMuon, runno, lumiblock, evtno);
                cout << messageBuffer;
              }
              break;
            }
          }
          if (!passIso) continue;

          if (debugFlag) cout << "    EleEt for electron number " << iElectron << " is: " << EleEt[iElectron] << endl;
          ST += EleEt[iElectron];
          multiplicity+=1;
          if (debugFlag) {
            sprintf(messageBuffer, "Ele number %d passed isolation in run number %d lumi section %d event number %d.      \n It had Px=%f and Py=%f\n", iElectron, runno, lumiblock, evtno, ElePx[iElectron], ElePy[iElectron]);
            cout << messageBuffer;
          }
          Px += ElePx[iElectron];
          Py += ElePy[iElectron];
          if (debugFlag) {
            sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
            cout << messageBuffer;
          }
        }
        else break;
      }
    }

    //Photons
    if (eventHasPhoton) {
      for (int iPhoton = 0; iPhoton < 25; ++iPhoton) {
        passIso=true;
        if (PhEt[iPhoton]>50.) {
          for (int iJet = 0; iJet < 25; ++iJet ) {
            if (JetEt[iJet]>50 && dR(PhEta[iPhoton],PhPhi[iPhoton], JetEta[iJet], JetPhi[iJet]) < 0.05) {
              if (PhEt[iPhoton]<0.5*JetEt[iJet]) {
                passIso = false;
                break;
              }
            }
          }
          if (!passIso) continue;

          // Throw out photon if there's a photon/muon overlap
          for (int iMuon = 0; iMuon < 25; ++iMuon ) {
            if (MuEt[iMuon]>50 && dR(PhEta[iPhoton], PhPhi[iPhoton], MuEta[iMuon], MuPhi[iMuon]) < 0.05) {
              if (debugFlag) {
                sprintf(messageBuffer, "Photon number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %d\n", iPhoton, iMuon, runno, lumiblock, evtno);
                cout << messageBuffer;
              }
              passIso = false;
              break;
            }
          }
          if (!passIso) continue;
          
          // Throw out photon if there's a photon/electron overlap
          for (int iElectron = 0; iElectron < 25; ++iElectron ) {
            if (EleEt[iElectron]>50 && dR(PhEta[iPhoton], PhPhi[iPhoton], EleEta[iElectron], ElePhi[iElectron]) < 0.05) {
              if (debugFlag) {
                sprintf(messageBuffer, "Photon number %d failed isolation with Electron number %d  in run number %d lumi section %d event number %d\n", iPhoton, iElectron, runno, lumiblock, evtno);
                cout << messageBuffer;
              }
              passIso = false;
              break;
            }
          }
          if (!passIso) continue;

          if (debugFlag) cout << "    PhEt for photon number " << iPhoton << " is: " << PhEt[iPhoton] << endl;
          ST += PhEt[iPhoton];
          multiplicity+=1;
          if (debugFlag) {
            sprintf(messageBuffer, "Photon number %d passed isolation in run number %d lumi section %d event number %d.\n      It had Px=%f and Py=%f\n", iPhoton, runno, lumiblock, evtno, PhPx[iPhoton], PhPy[iPhoton]);
            cout << messageBuffer;
          }
          Px += PhPx[iPhoton];
          Py += PhPy[iPhoton];
          if (debugFlag) {
            sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
            cout << messageBuffer;
          }
        }
        else break;
      }
    }

    //Muons
    if (eventHasMuon) {
      for (int iMuon = 0; iMuon < 25; ++iMuon) {
        passIso=true;
        if (MuEt[iMuon]>50.) {
          if (debugFlag) cout << "    MuEt for muon number " << iMuon << " is: " << MuEt[iMuon] << endl;
          ST += MuEt[iMuon];
          multiplicity+=1;
          if (debugFlag) {
            sprintf(messageBuffer, "Muon number %d passed isolation in run number %d lumi section %d event number %d.\n       It had Px=%f and Py=%f\n", iMuon, runno, lumiblock, evtno, MuPx[iMuon], MuPy[iMuon]);
          }
          Px += MuPx[iMuon];
          Py += MuPy[iMuon];
        }
        else break;
      }
    }


    // MET / MHT
    if (debugFlag) cout << "    Met from PAT collection is: " << Met << endl;
    MHT = std::sqrt(Px*Px + Py*Py);
    if (debugFlag) cout << "    MHT is: " << MHT << endl;
    STMHTnoMET = ST + MHT;
    ST += Met;
    if (debugFlag) cout << "    ST is: " << ST << endl;
    if (debugFlag) cout << "    ST using MHT is: " << STMHTnoMET << endl;
    if (debugFlag) cout << "    multiplicity is: " << multiplicity << endl;

    // Fill all the histograms
    stHist.Fill(ST);
    stHistMHT.Fill(STMHTnoMET);
    for (int iHist = 0; iHist<multMax-2; ++iHist) {
      if (multiplicity == iHist+2) stExcHist[iHist]->Fill(ST);
      if (multiplicity >= iHist+2) stIncHist[iHist]->Fill(ST);
    }
    for (int iHist = 0; iHist<multMax-2; ++iHist) {
      if (multiplicity == iHist+2) stExcHistMHT[iHist]->Fill(STMHTnoMET);
      if (multiplicity >= iHist+2) stIncHistMHT[iHist]->Fill(STMHTnoMET);
    }
    ++nDumpedEvents;
    sprintf(messageBuffer, "\n\n--------------------------------------------------------------\n");
    cout << messageBuffer;
    if (debugFlag && nDumpedEvents==eventsToDump) break;
  }
  // write output root file
  outFile->cd();
  stHist.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHist[iHist]->Write();
    stIncHist[iHist]->Write();
  }
  stHistMHT.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHistMHT[iHist]->Write();
    stIncHistMHT[iHist]->Write();
  }
  outFile->Close();
}

// function to calculate dR between two objects
float dR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt( ( eta1 - eta2 )*( eta1 - eta2 ) + std::pow(TMath::ATan2(TMath::Sin( phi1 - phi2), TMath::Cos(phi1-phi2)),2) );
}
