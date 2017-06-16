#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TChain.h"
#include "rootIncludes.inc"
#include "prepareEvents.C"

void BookHistosReco();
void WriteHistosReco(Char_t *fNameOut);
//==========================================
int main(int argc, char** argv){

  Char_t *fNameOut = (char*)"tmpFiles/selEvents_data_Ups.root";
  bool rejectCowboys = false;
  bool RequestTrigger=true;

  Char_t const *inputTree1 = "Default";

  int inputTrees=0;


  for( int i=0;i < argc; ++i ) {
    if(std::string(argv[i]).find("rejectCowboys=false") != std::string::npos) {rejectCowboys=kFALSE;}
    if(std::string(argv[i]).find("RequestTrigger=1") != std::string::npos) {RequestTrigger=true;}
    if(std::string(argv[i]).find("inputTree1") != std::string::npos) {inputTrees++; char* inputTree1char = argv[i]; char* inputTree1char2 = strtok (inputTree1char, "="); inputTree1 = inputTree1char2; cout<<"inputTree1 = "<<inputTree1<<endl;}
    //          if(std::string(argv[i]).find("inputTree2") != std::string::npos) {inputTrees++; char* inputTree2char = argv[i]; char* inputTree2char2 = strtok (inputTree2char, "="); inputTree2 = inputTree2char2; cout<<"inputTree2 = "<<inputTree2<<endl;}
    //          if(std::string(argv[i]).find("inputTree3") != std::string::npos) {inputTrees++; char* inputTree3char = argv[i]; char* inputTree3char2 = strtok (inputTree3char, "="); inputTree3 = inputTree3char2; cout<<"inputTree3 = "<<inputTree3<<endl;}
    //          if(std::string(argv[i]).find("inputTree4") != std::string::npos) {inputTrees++; char* inputTree4char = argv[i]; char* inputTree4char2 = strtok (inputTree4char, "="); inputTree4 = inputTree4char2; cout<<"inputTree4 = "<<inputTree4<<endl;}
  }

  cout<<"Number of Input Trees = "<<inputTrees<<endl;

  TChain *chain = new TChain("data");
  chain->Add(inputTree1);

  TTree *tree = chain;

  TFile *fOut = fOut = new TFile(fNameOut, "RECREATE");

  prepareEvents treeReco(tree);
  BookHistosReco();
  treeReco.Loop(RequestTrigger, rejectCowboys);
  printf("writing out the histograms\n");
  WriteHistosReco(fNameOut);

  fOut->Close();

  return 0;

}
//==========================================
void BookHistosReco(){

  //mass
  Int_t nBinsMass = 320;
  Double_t massMin = 8.4, massMax = 11.6;
  //pt
  Int_t nBinsPt = 1000;
  Double_t pTMin = 0., pTMaxOnia = 100.;
  //rap
  Int_t nBinsRap = 100;
  Double_t rapMin = -2.5, rapMax = 2.5;

  Char_t name[100], title[300];
  //statistics
  Reco_StatEv = new TH1F("Reco_StatEv", "", 12, 0., 12.);

  //reconstruction variables for the Onia

  sprintf(name, "Reco_Onia_mass");
  sprintf(title, ";M(#mu#mu) [GeV];");
  Reco_Onia_mass = new TH1F(name, title, nBinsMass,massMin,massMax);
  Reco_Onia_mass->Sumw2();

  sprintf(name, "Reco_Onia_rap_pt");
  sprintf(title, ";y(#mu#mu);p_{T}^{#mu#mu} [GeV]");
  Reco_Onia_rap_pT = new TH2F(name, title, nBinsRap,rapMin,rapMax, nBinsPt,pTMin,pTMaxOnia);
  Reco_Onia_rap_pT->Sumw2();


  //prepare the branches for the output tree
  treeOut = new TTree ("selectedData", "selected events");
  treeOut->SetAutoSave(0);

  lepP = new TLorentzVector();
  lepN = new TLorentzVector();

  treeOut->Branch("lepP", "TLorentzVector", &lepP);
  treeOut->Branch("lepN", "TLorentzVector", &lepN);


}

//==========================================
void WriteHistosReco(Char_t *fNameOut){

  treeOut->Write();

  Reco_StatEv->Write();
  Reco_Onia_rap_pT->Write();
  Reco_Onia_mass->Write();


}
