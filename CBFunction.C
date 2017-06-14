#include "TMath.h"

Double_t CBFunction(Double_t *x, Double_t *par){

  Double_t CB = 0.;
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  Double_t n = par[4];

  if(((x[0] - mean)/sigma) > -alpha){
    CB = TMath::Exp(-(pow(x[0]-mean,2)/(2.*sigma*sigma)));
  }
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB = A*pow(B - (x[0]-mean)/sigma, -n);
  } 
  return par[0]*CB; 
}
