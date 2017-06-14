#define prepareEvents_cxx
#include "prepareEvents.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TH1F *Reco_StatEv;
TTree *treeOut;
TLorentzVector *lepP, *lepN;
Double_t Nch;
TH2F *Reco_Onia_rap_pT;
TH1F *Reco_Onia_mass;

void prepareEvents::Loop(bool RequestTrigger, bool rejectCowboys)
{

   if (fChain == 0) return;
   
   Double_t Nch;
   treeOut->Branch("Nch", &Nch, "Nch/D");   

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0, countRecEvent = 0;
   
   
//   nentries=2500000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
      if(jentry % 1000000 == 0) printf("event %d out of %d\n", (Int_t) jentry, (Int_t) nentries);
      
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      

      
      if(onia->Pt() > 990.)
      continue;

      if(JpsiVprob < 0.01)
        continue;
      
      Reco_StatEv->Fill(0.5);
      
      Int_t trigDecision = -99;
      Int_t trigPtDecision = -99;
      
      //Upsilon trigger paths:
   if(//HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
      //HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
      //HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
      HLT_Dimuon5_Upsilon_Barrel_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
      HLT_Dimuon5_Upsilon_Barrel_v2 == 1 || //  1e33: 166346
      HLT_Dimuon5_Upsilon_Barrel_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 2)
      HLT_Dimuon5_Upsilon_Barrel_v5 == 1 || //2E33 (no cowboys)
      HLT_Dimuon7_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ)
      HLT_Dimuon7_Upsilon_Barrel_v4 == 1 || //5E33 (becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v4 == 1) //5E33
      trigDecision = 1;

   if(trigDecision != 1 && RequestTrigger)
     continue;


   if(//HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
      //HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
      //HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
      HLT_Dimuon5_Upsilon_Barrel_v1 == 1 && onia->Pt() < 5.5 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
      HLT_Dimuon5_Upsilon_Barrel_v2 == 1 && onia->Pt() < 5.5  || //  1e33: 166346
      HLT_Dimuon5_Upsilon_Barrel_v3 == 1 && onia->Pt() < 5.5  ||  //1.4e33: 167078 - 167913 (prescale of 2)
      HLT_Dimuon5_Upsilon_Barrel_v5 == 1 && onia->Pt() < 5.5  || //2E33 (no cowboys)
      HLT_Dimuon7_Upsilon_Barrel_v1 == 1 && onia->Pt() < 7.5  || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v1 == 1 && onia->Pt() < 9.5  || //3E33 (L1_DoubleMu0_HighQ)
      HLT_Dimuon7_Upsilon_Barrel_v4 == 1 && onia->Pt() < 7.5  || //5E33 (becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v4 == 1 && onia->Pt() < 9.5 ) //5E33
      trigPtDecision = 1;

   if(trigPtDecision == 1)
     continue;
     
    Reco_StatEv->Fill(1.5); //count all events
      
	Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
//    Double_t onia_P = onia->P();
//    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_cpm = vertexWeight;
    
//    if(TMath::Abs(onia_rap) > onia::rapYPS)
	if(TMath::Abs(onia_rap) > 1.2)
      continue;
    Reco_StatEv->Fill(2.5);
    
    Double_t deltaPhi = muNeg->Phi() - muPos->Phi();
	if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
	else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

//	if(rejectCowboys)
//		  if(deltaPhi < 0.) continue;
		  
		  
	Reco_StatEv->Fill(3.5);
	
	if(onia_mass < 8.4 || onia_mass > 11.6) //all Upsilons triggered
      continue;
            
    Reco_StatEv->Fill(4.5);
    
    Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    Reco_Onia_mass->Fill(onia_mass);
    
    countRecEvent++;
	Nch=vertexWeight;
    lepP = muPos;
    lepN = muNeg;
    treeOut->Fill();
    
   }
   
printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}
