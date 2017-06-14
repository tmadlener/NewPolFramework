#include "rootIncludes.inc"
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TROOT.h"

#include "ReshuffleNch.C"



//====================================
int main(){

	  gROOT->ProcessLine(".L ReshuffleNch.C+");

        ReshuffleNch();	

  return 0;
}
