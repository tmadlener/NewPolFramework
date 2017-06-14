#include "rootIncludes.inc"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include <vector>
//#include <random>
#include <algorithm>
#include <iostream>


inline int getRandom(const int N)
{
  return gRandom->Integer(N);
}


void ReshuffleNch(){

Char_t fileNameIn[100];
sprintf(fileNameIn, "tmpFiles/data.root");
TFile *fIn = new TFile(fileNameIn);
TTree *treeIn = (TTree *) gDirectory->Get("selectedData");
Char_t fileNameOut[100];
sprintf(fileNameOut, "tmpFiles/data_ReshuffledNch.root");
TFile *fOut = new TFile(fileNameOut, "RECREATE");
TTree *treeOut =  new TTree ("reshuffledData", "reshuffled events");

Double_t pT_, mass_, rap_, costh_CS_, costh_HX_, phi_CS_, phi_HX_, Nch_, w_Y1S_, w_Y2S_, w_Y3S_;
Double_t pT, mass, rap, costh_CS, costh_HX, phi_CS, phi_HX, Nch, w_Y1S, w_Y2S, w_Y3S;
treeIn->SetBranchAddress("pT", &pT);
treeIn->SetBranchAddress("mass", &mass);
treeIn->SetBranchAddress("rap", &rap);
treeIn->SetBranchAddress("costh_HX", &costh_HX);
treeIn->SetBranchAddress("costh_CS", &costh_CS);
treeIn->SetBranchAddress("phi_CS", &phi_CS);
treeIn->SetBranchAddress("phi_HX", &phi_HX);
treeIn->SetBranchAddress("Nch", &Nch);
treeIn->SetBranchAddress("w_Y1S", &w_Y1S);
treeIn->SetBranchAddress("w_Y2S", &w_Y2S);
treeIn->SetBranchAddress("w_Y3S", &w_Y3S);

treeOut->Branch("pT", &pT_, "pT/D");
treeOut->Branch("mass", &mass_, "mass/D");
treeOut->Branch("rap", &rap_, "rap/D");
treeOut->Branch("costh_HX", &costh_HX_, "costh_HX/D");
treeOut->Branch("phi_HX", &phi_HX_, "phi_HX/D");
treeOut->Branch("costh_CS", &costh_CS_, "costh_CS/D");
treeOut->Branch("phi_CS", &phi_CS_, "phi_CS/D");
treeOut->Branch("Nch", &Nch_, "Nch/D");
treeOut->Branch("w_Y1S", &w_Y1S_, "w_Y1S/D");
treeOut->Branch("w_Y2S", &w_Y2S_, "w_Y2S/D");
treeOut->Branch("w_Y3S", &w_Y3S_, "w_Y3S/D");

Long64_t nentries=treeIn->GetEntries();
std::vector<double> vecNch;
vecNch.reserve(nentries);
treeIn->LoadBaskets();

for(int iEn = 0; iEn < nentries; iEn++){ 
  
  Long64_t iEntry = treeIn->LoadTree(iEn);
  
  treeIn->GetEntry(iEntry);
  
  vecNch.push_back(Nch);
  
}

  if (gRandom) delete gRandom;
  gRandom = new TRandom3();
  std::random_shuffle(vecNch.begin(), vecNch.end(), getRandom);
cout<<"Nch vec filled and randomized"<<endl;

cout<<"Filling trees with kinematic/angular info with reshuffled Nch"<<endl;

for(int iEn = 0; iEn < nentries; iEn++){ 
  
  Long64_t iEntry = treeIn->LoadTree(iEn);
  
    treeIn->GetEntry(iEntry);
    if(iEn % 100000 == 0)
	  cout << "entry " << iEntry << " out of " << nentries << endl;

	pT_=pT;
	mass_=mass;
	rap_=rap;
	costh_CS_=costh_CS;
	costh_HX_=costh_HX;
	phi_CS_=phi_CS;
	phi_HX_=phi_HX;
	w_Y1S_=w_Y1S;
	w_Y2S_=w_Y2S;
	w_Y3S_=w_Y3S;
	Nch_=vecNch[iEn];
	
	treeOut->Fill();
}
  

  fOut->cd();
  treeOut->Write();
  fOut->Close();
  fIn->Close();

}
