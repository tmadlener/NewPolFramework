#define prepareEvents_cxx
#include "prepareEvents.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// lifetimecuts, hardcoded here, toggle application with boolean

bool ctauCut(const double ctau, const double ctauErr)
{
  const double ctauMax = 1000.0;
  const double ctauMin = -1000.0;
  const double ctauSigMax = 2.0;
  const double ctauSigMin = -1000.0;

  const double ctauSig = std::abs(ctau / ctauErr);

  if (ctau < ctauMin) return false;
  if (ctau > ctauMax) return false;
  if (ctauSig < ctauSigMin) return false;
  if (ctauSig > ctauSigMax) return false;

  return true;
}

bool fiducialCuts(const double pt, const double absEta)
{
  if (absEta < 1.2 && pt < 4.5) return false;
  if (absEta > 1.2 && absEta < 1.4 && pt < 3.5) return false;
  if (absEta > 1.4 && absEta < 1.6 && pt < 3.0) return false;

  return true;
}

TH1F *Reco_StatEv;
TTree *treeOut;
TLorentzVector *lepP, *lepN;
TH2F *Reco_Onia_rap_pT;
TH1F *Reco_Onia_mass;

void prepareEvents::Loop(bool RequestTrigger, bool rejectCowboys, bool applyCtauCut, bool fidCuts)
{
  if (fChain == 0) return;

  Double_t Nch;
  treeOut->Branch("Nch", &Nch, "Nch/D");

  double ctauOut;
  treeOut->Branch("ctau", &ctauOut);
  double ctauErrOut;
  treeOut->Branch("ctauErr", &ctauErrOut);

  double mupPt;
  double munPt;
  double mupEta;
  double munEta;

  treeOut->Branch("mupPt", &mupPt);
  treeOut->Branch("munPt", &munPt);
  treeOut->Branch("mupEta", &mupEta);
  treeOut->Branch("munEta", &munEta);


  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0, countRecEvent = 0;

  // nentries=25000000;
  for (Long64_t jentry=33000000; jentry<nentries;jentry++) {
    if(jentry % 1000000 == 0) printf("event %d out of %d\n", (Int_t) jentry, (Int_t) nentries);


    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;



    if(onia->Pt() > 990.)
      continue;

    if(JpsiVprob < 0.01)
      continue;

    Reco_StatEv->Fill(0.5);

    if(RequestTrigger) {
      if (!(
            //HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
            //HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
            //HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
            // HLT_Dimuon5_Upsilon_Barrel_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
            // HLT_Dimuon5_Upsilon_Barrel_v2 == 1 || //  1e33: 166346
            // HLT_Dimuon5_Upsilon_Barrel_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 2)
            // HLT_Dimuon5_Upsilon_Barrel_v5 == 1 || //2E33 (no cowboys)
            // HLT_Dimuon7_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
            HLT_Dimuon9_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ)
            // HLT_Dimuon7_Upsilon_Barrel_v4 == 1 || //5E33 (becomes inactive for Linst >= 5E33)
            HLT_Dimuon9_Upsilon_Barrel_v4 == 1 //5E33
            )) {
        continue;
      }
    }

    Double_t onia_pt = onia->Pt();
    if (!(
          // HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
          // HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
          // HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
          // (HLT_Dimuon5_Upsilon_Barrel_v1 == 1 && onia->Pt() > 5.5) || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
          // (HLT_Dimuon5_Upsilon_Barrel_v2 == 1 && onia->Pt() > 5.5)  || //  1e33: 166346
          // (HLT_Dimuon5_Upsilon_Barrel_v3 == 1 && onia->Pt() > 5.5)  ||  //1.4e33: 167078 - 167913 (prescale of 2)
          // (HLT_Dimuon5_Upsilon_Barrel_v5 == 1 && onia->Pt() > 5.5)   || //2E33 (no cowboys)
          // (HLT_Dimuon7_Upsilon_Barrel_v1 == 1 && onia->Pt() > 7.5)  || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
          (HLT_Dimuon9_Upsilon_Barrel_v1 == 1 && onia_pt > 9.5)  || //3E33 (L1_DoubleMu0_HighQ)
          // (HLT_Dimuon7_Upsilon_Barrel_v4 == 1 && onia->Pt() > 7.5)  || //5E33 (becomes inactive for Linst >= 5E33)
          (HLT_Dimuon9_Upsilon_Barrel_v4 == 1 && onia_pt > 9.5) //5E33
          )) {
      continue;
    }

    Reco_StatEv->Fill(1.5); //count all events

    Double_t onia_mass = onia->M();
    //    Double_t onia_P = onia->P();
    //    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();

    //    if(TMath::Abs(onia_rap) > onia::rapYPS)
    if(TMath::Abs(onia_rap) > 1.2)
      continue;
    Reco_StatEv->Fill(2.5);

    if (onia_pt < 10.0) continue;

    Reco_StatEv->Fill(3.5);

    Double_t deltaPhi = muNeg->Phi() - muPos->Phi();
    if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

    if(rejectCowboys)
      if(deltaPhi < 0.) continue;


    Reco_StatEv->Fill(4.5);

    if(onia_mass < 8.4 || onia_mass > 11.6) //all Upsilons triggered
      continue;
    Reco_StatEv->Fill(5.5);


    if (std::isnan(Jpsict)) continue; // sanity check that sometimes is necessary
    if (applyCtauCut && !ctauCut(Jpsict, JpsictErr)) continue;
    Reco_StatEv->Fill(6.5);

    if (fidCuts) {
      if (!fiducialCuts(muPos->Pt(), std::abs(muPos->PseudoRapidity())) ||
          !fiducialCuts(muNeg->Pt(), std::abs(muNeg->PseudoRapidity()))) {
        continue;
      }
    }
    Reco_StatEv->Fill(7.5);

    Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    Reco_Onia_mass->Fill(onia_mass);

    countRecEvent++;
    Nch=vertexWeight;
    lepP = muPos;
    lepN = muNeg;

    ctauOut = Jpsict;
    ctauErrOut = JpsictErr;

    mupPt = lepP->Pt();
    munPt = lepN->Pt();
    mupEta = lepP->Eta();
    munEta = lepN->Eta();

    treeOut->Fill();

  }

  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}
