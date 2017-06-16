#include "rootIncludes.inc"
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TROOT.h"

#include "upsilon_2StepFit.C"



//====================================
int main(){

  Double_t nSigma = 3.; //sigma_mass cut for preparation of figures
  Char_t *fileNameIn = (char*)"tmpFiles/selEvents_data_Ups.root";


  gROOT->ProcessLine(".L upsilon_2StepFit.C+");

  upsilon_2StepFit(nSigma, fileNameIn);

  return 0;
}
