#include "rootIncludes.inc"
#include "CBFunction.C"
#include "TMath.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFrame.h"

Double_t peak_min = 8.9, peak_max = 10.6; //veto this range when fitting the continuum

const Double_t massPDG1S = 9.460;
const Double_t massPDG2S = 10.023;
const Double_t massPDG3S = 10.355;

Int_t const kNbSpecies = 3;
Char_t const *specName[kNbSpecies] = {"1S", "2S", "3S"};
enum {UPS1S, UPS2S, UPS3S, BG};
Int_t colour[kNbSpecies+1] = {kRed-9,kGreen-3,kBlue-8,kGray};
TH1F *hMass;
Double_t binWidth;
TF1* fRECO;
Double_t fitParBG[3];
Double_t fitParBGerr[3];
Double_t intCB[kNbSpecies]; //integral values of a CB with N=1 and fixed alpha, n, sigma, width
bool PlotSimplistic=false;

Double_t massMin[kNbSpecies], massMax[kNbSpecies];
TF1 *fUps1S, *fUps2S, *fUps3S, *fBG;
Double_t fracBG[kNbSpecies];
Double_t nY[kNbSpecies];

double sigma1S_save;

void GetHisto(Char_t *fileNameIn);
void FitSignalBG(Double_t nSigma);
Double_t fitPolyCrystal3(Double_t *x, Double_t *par);
Double_t fitContinuum(Double_t *x, Double_t *par);
Double_t DrawContinuum(Double_t *x, Double_t *par);
//void DrawFit(Double_t nSigma);
void SaveCBParameters(Double_t alpha, Double_t n, Double_t alphaErr, Double_t nErr);
void SaveFitPars();
//==============================
void upsilon_2StepFit(
                      Double_t nSigma = 2.,
                      Char_t *fileNameIn = (char*) "RootFiles/selEvents_data_Ups_2Aug2011.root"){

  GetHisto(fileNameIn);

  if(hMass->GetEntries() < 1.){
    printf("\n\n\nskip processing this bin, because the number of entries is smaller than 1\n\n\n");
    return;
  }

  FitSignalBG(nSigma);
  //  DrawFit(nSigma);
}

//===============================
#if 0
void DrawFit(Double_t nSigma){

  gStyle->SetFillColor(0);

  Char_t name[100];
  gStyle->SetFrameBorderMode(0);

  //prepare the drawing of the individual components:
  fBG->SetFillColor(colour[BG]);
  fBG->SetLineColor(colour[BG]);
  fBG->SetFillStyle(1001);
  fBG->SetNpx(1000);
  TH1 *hBG = fBG->GetHistogram();

  fUps1S->SetNpx(1000);
  fUps1S->SetFillColor(colour[UPS1S]);
  fUps1S->SetLineColor(colour[UPS1S]);
  fUps1S->SetFillStyle(1001);

  fUps2S->SetNpx(1000);
  fUps2S->SetFillColor(colour[UPS2S]);
  fUps2S->SetLineColor(colour[UPS2S]);
  fUps2S->SetFillStyle(1001);

  fUps3S->SetNpx(1000);
  fUps3S->SetFillColor(colour[UPS3S]);
  fUps3S->SetLineColor(colour[UPS3S]);
  fUps3S->SetFillStyle(1001);

  TH1 *hUps1S = fUps1S->GetHistogram();
  TH1 *hUps2S = fUps2S->GetHistogram();
  TH1 *hUps3S = fUps3S->GetHistogram();

  THStack *hStack = new THStack("hMass_Stack", "");
  hStack->Add(hBG);
  hStack->Add(hUps3S);
  hStack->Add(hUps2S);
  hStack->Add(hUps1S);
  if(!PlotSimplistic) hStack->Draw("same");

  hMass->Draw("same");
  if(!PlotSimplistic) fRECO->Draw("same");

  TLine *line[3];
  Double_t max[3] = {1., 0.5, 0.3};
  for(int iL = 0; iL < 3; iL++){
    line[iL]= new TLine(massMin[iL], 0.1, massMin[iL], max[iL]*hUps1S->GetMaximum());
    line[iL]->SetLineStyle(2); line[iL]->SetLineColor(colour[iL]);
    line[iL]->SetLineWidth(2);
    if(!PlotSimplistic) line[iL]->Draw();
    if(!PlotSimplistic) line[iL]->DrawLine(massMax[iL], 0.1, massMax[iL], max[iL]*hUps1S->GetMaximum());
  }


  /*  double MassScan[13]={8.6,8.95,9.3,9.45,9.6,9.85,10.0125,10.175,10.3425,10.51,10.8,11.1,11.4};
      TLine *line[13];
      for(int iL = 0; iL < 13; iL++){
      line[iL]= new TLine(MassScan[iL], 0.1, MassScan[iL], 1.1*hUps1S->GetMaximum());
      line[iL]->SetLineStyle(2); line[iL]->SetLineColor(kWhite);
      line[iL]->SetLineWidth(2);
      line[iL]->Draw();
      }

      TLatex *texMassScan[13];
      char MassScanName[200];
      for(int iL = 0; iL < 13; iL++){
      sprintf(MassScanName,"%d",iL+1);
      texMassScan[iL] = new TLatex((MassScan[iL]+MassScan[iL+1])/2., 0.015*hStack->GetMaximum(), MassScanName);
      texMassScan[iL]->SetTextSize(0.03);
      texMassScan[iL]->SetTextColor(kWhite);
      texMassScan[iL]->Draw();
      }
  */
  if(iRapBin == 0) sprintf(name, "|y| < %1.1f", onia::rapYPS);
  else if(iRapBin == 1) sprintf(name, "|y| < %1.1f", onia::rapForPTRange[iRapBin]);
  else if(iRapBin > 1)  sprintf(name, "%1.1f < |y| < %1.1f", onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin]);
  double xText=10.25;
  TLatex *tex = new TLatex(xText, hStack->GetMaximum(), name);
  tex->SetTextSize(0.04);
  if(!PlotSimplistic) tex->Draw();

  if(iPTBin == 0) sprintf(name, "all p_{T}");
  else if(iPTBin > 0) sprintf(name, "%1.1f < p_{T} < %1.1f", onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin]);
  if(!PlotSimplistic)  tex->DrawLatex(xText, 0.94*hStack->GetMaximum(), name);

  if(iCPMBin == 0) sprintf(name, "all N_{ch}");
  else if(iCPMBin > 0) sprintf(name, "%1.1f < N_{ch} < %1.1f", onia::cpmRange[iCPMBin-1], onia::cpmRange[iCPMBin]);
  if(!PlotSimplistic)  tex->DrawLatex(xText, 0.86*hStack->GetMaximum(), name);


  sprintf(name, "frac(BG) in #pm %1.1f#sigma:", nSigma);
  if(!PlotSimplistic)  tex->DrawLatex(xText, 0.80*hStack->GetMaximum(), name);
  sprintf(name, "%1.2f, %1.2f, %1.2f", fracBG[0], fracBG[1], fracBG[2]);
  if(!PlotSimplistic) tex->DrawLatex(xText, 0.74*hStack->GetMaximum(), name);

  sprintf(name, "Figures/massFit_rap%d_pT%d_cpm%d.pdf", iRapBin, iPTBin, iCPMBin);
  //  gPad->SetLogy(kTRUE);
  if(iRapBin == 0 && iPTBin == 0) gPad->Print(name);
  else if(iRapBin > 0 && iPTBin > 0) gPad->Print(name);




  if(iPTBin > -1 && iRapBin > -10 && iCPMBin >- 1){


    /// produce pedagogical plot

    cout<<"Plot pedagogical"<<endl;

    double mean1S_draw = fUps1S->GetParameter(1);
    double sigma1S_draw = fUps1S->GetParameter(2);

    double nSigmaMin=0;
    double nSigmaMax=2.5;
    int nIntegrals=100;

    double nSigmaCenter[nIntegrals];
    double lSig[nIntegrals];
    double lBkg[nIntegrals];
    double lSig_[nIntegrals];
    double lBkg_[nIntegrals];
    double lSigOVERBkg[nIntegrals];

    double maxSig=2.5;


    int whichBinIsAtOne=1/(nSigmaMax-nSigmaMin)*nIntegrals;
    cout<<"whichBinIsAtOne "<<whichBinIsAtOne<<endl;

    for(int nIter=0;nIter<nIntegrals;nIter++){
      nSigmaCenter[nIter]= (nSigmaMax-nSigmaMin)/double(nIntegrals)*double(nIter+1);
      lSig_[nIter]= fUps1S->Integral(mean1S_draw-nSigmaCenter[nIter]*sigma1S_draw, mean1S_draw+nSigmaCenter[nIter]*sigma1S_draw);
      lBkg_[nIter]= fBG->Integral(mean1S_draw-nSigmaCenter[nIter]*sigma1S_draw, mean1S_draw+nSigmaCenter[nIter]*sigma1S_draw);

    }

    for(int nIter=0;nIter<nIntegrals;nIter++){
      nSigmaCenter[nIter]= (nSigmaMax-nSigmaMin)/double(nIntegrals)*double(nIter+1);
      lSig[nIter]= fUps1S->Integral(mean1S_draw-nSigmaCenter[nIter]*sigma1S_draw, mean1S_draw+nSigmaCenter[nIter]*sigma1S_draw)/lSig_[whichBinIsAtOne-1];
      lBkg[nIter]= fBG->Integral(mean1S_draw-nSigmaCenter[nIter]*sigma1S_draw, mean1S_draw+nSigmaCenter[nIter]*sigma1S_draw)/lBkg_[whichBinIsAtOne-1];
    }

    for(int nIter=0;nIter<nIntegrals;nIter++){
      lSigOVERBkg[nIter]= lSig[nIter]/lBkg[nIter];
      cout<<"lSigOVERBkg[nIter] "<<lSigOVERBkg[nIter]<<endl;

    }

    TGraphErrors *nSigma_Sig = new TGraphErrors(nIntegrals,nSigmaCenter,lSig,0,0);
    TGraphErrors *nSigma_Bkg = new TGraphErrors(nIntegrals,nSigmaCenter,lBkg,0,0);
    TGraphErrors *nSigma_SigOVERBkg = new TGraphErrors(nIntegrals,nSigmaCenter,lSigOVERBkg,0,0);


    TCanvas *SystCanvas = new TCanvas("SystCanvas","SystCanvas",1000,800);
    //                gStyle->SetPalette(1);
    //                gPad->SetFillColor(kWhite);
    //            gPad->SetLeftMargin(0.15);

    SystCanvas->SetFillColor(kWhite);
    SystCanvas->SetGrid();
    SystCanvas->GetFrame()->SetFillColor(kWhite);
    SystCanvas->GetFrame()->SetBorderSize(0);
    SystCanvas->SetRightMargin(0.05) ;


    TLegend* plotLegend=new TLegend(0.1,0.7,0.45,0.9);
    plotLegend->SetFillColor(kWhite);
    plotLegend->SetTextFont(72);
    plotLegend->SetTextSize(0.03);
    plotLegend->SetBorderSize(1);
    char legendentry[200];

    double lineWidth=3;
    char drawGraphStyle[200];
    sprintf(drawGraphStyle,"LX");


    TH1F *SystHisto = new TH1F;
    SystHisto = SystCanvas->DrawFrame(nSigmaMin,0,nSigmaMax,maxSig);
    SystHisto->SetXTitle("n_{#sigma}");
    SystHisto->GetYaxis()->SetTitleOffset(1.5);

    nSigma_Sig->SetLineColor(kGreen-2);
    nSigma_Sig->SetLineWidth(lineWidth);
    nSigma_Sig->Draw(drawGraphStyle);
    sprintf(legendentry,"Signal yield");
    plotLegend->AddEntry(nSigma_Sig,legendentry,"l");

    nSigma_Bkg->SetLineColor(kBlue);
    nSigma_Bkg->SetLineWidth(lineWidth);
    nSigma_Bkg->Draw(drawGraphStyle);
    sprintf(legendentry,"Background yield");
    plotLegend->AddEntry(nSigma_Bkg,legendentry,"l");

    nSigma_SigOVERBkg->SetLineColor(kRed);
    nSigma_SigOVERBkg->SetLineWidth(lineWidth);
    nSigma_SigOVERBkg->Draw(drawGraphStyle);
    sprintf(legendentry,"S/B ratio");
    plotLegend->AddEntry(nSigma_SigOVERBkg,legendentry,"l");

    plotLegend->Draw();

    sprintf(name,"Figures/pedagogical_rap%d_pT%d_cpm%d.pdf",iRapBin,iPTBin, iCPMBin);
    if(iRapBin == 0 && iPTBin == 0) SystCanvas->SaveAs(name);
    else if(iRapBin > 0 && iPTBin > 0) SystCanvas->SaveAs(name);

  }


}
#endif

//===============================
void FitSignalBG(Double_t nSigma){

  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat(kFALSE);

  Char_t name[100];
  //1.) perform the fit to the continuum, using the sidebands
  sprintf(name, "c1");
  TCanvas *c1 = new TCanvas(name);
  sprintf(name, "Counts per %1.0f MeV", 1000.*binWidth);
  c1->SetFillColor(0);
  c1->SetLeftMargin(0.15);
  hMass->SetYTitle(name);
  hMass->SetXTitle("M_{#mu#mu} [GeV]");
  hMass->SetTitleOffset(1.5, "y");
  hMass->SetStats(0);
  hMass->Draw();

  //starting values for fit:
  Double_t a = -2.0e6, b = 5.e5, c = -2.4e4;
  Double_t range_min = 8.6, range_max = 11.4;
  Double_t normY1S = hMass->GetMaximum() / binWidth;
  Double_t normY2S = 0.3*normY1S;
  Double_t normY3S = 0.15*normY1S;
  Double_t sigma1S = 0.1, mean1S = 9.45;
  Double_t alpha = 1.33, n = 6.6; //CB-tail parameters

  printf("will be fitting the continuum between %1.2f < M < %1.2f\n", range_min, range_max);
  sprintf(name, "fCont");
  TF1* fCONT = new TF1(name, fitContinuum, range_min, range_max, 3);
  fCONT->SetParameters(a, b, c);
  hMass->Fit(fCONT, "0", "", range_min, range_max);
  fCONT = hMass->GetFunction(name);

  Double_t chisqrd = fCONT->GetChisquare();
  Int_t NDF = fCONT->GetNDF();
  Double_t redChi2 = chisqrd/NDF;
  printf("\nChisqrd = %1.3f, NDF = %d, Chisqrd/NDF = %1.3f, Prob = %1.3f\n\n", chisqrd, NDF, redChi2, TMath::Prob(chisqrd,NDF));

  fCONT->GetParameters(fitParBG);
  for(int iPar = 0; iPar < 3; iPar++){
    fitParBGerr[iPar]=fCONT->GetParError(iPar);
  }

  //2.) fit the peaks on the top of a fixed continuum
  sprintf(name,"fPeaks");
  Int_t const npar = 10;
  fRECO = new TF1(name, fitPolyCrystal3, peak_min, peak_max, npar);
  fRECO->SetParNames("normY1S", "mass_Ups1S", "sigma_Ups1S", "normY2S", "normY3S", "n", "alpha", "a", "b", "c");
  fRECO->SetParameters(normY1S, mean1S, sigma1S, normY2S, normY3S, n, alpha, a, b, c);
  fRECO->FixParameter(7,fitParBG[0]);
  fRECO->FixParameter(8,fitParBG[1]);
  fRECO->FixParameter(9,fitParBG[2]);
  //fix alpha and n from the fit to all bins
  /*
    Char_t fileName[100];
    sprintf(fileName, "tmpFiles/CBParameters.root");
    TFile *fIn = new TFile(fileName);
    TTree *treeIn = (TTree *) gDirectory->Get("CBPars");
    Double_t alphaAll, nAll;
    TBranch *b_alphaAll, *b_nAll;
    treeIn->SetBranchAddress("alphaAll", &alphaAll, &b_alphaAll);
    treeIn->SetBranchAddress("nAll", &nAll, &b_nAll);
    Long64_t iEntry = treeIn->LoadTree(0);
    treeIn->GetEntry(0);
    printf("alpha and n from File: %1.3f, %1.3f\n", alphaAll, nAll);
    fRECO->FixParameter(5, nAll);
    fRECO->FixParameter(6, alphaAll);
  */
  hMass->Fit(fRECO, "0", "", peak_min, peak_max);

  fRECO = hMass->GetFunction(name);
  fRECO->SetLineWidth(1);

  Double_t fitParTot[npar];
  fRECO->GetParameters(fitParTot);
  normY1S = fitParTot[0];   printf("normY1S = %1.3e\n", normY1S);
  mean1S = fitParTot[1];
  sigma1S = fitParTot[2];
  normY2S = fitParTot[3];
  normY3S = fitParTot[4];
  n = fitParTot[5];  printf("n = %1.3f\n", n);
  alpha = fitParTot[6];  printf("alpha = %1.3f\n", alpha);
  sigma1S_save=sigma1S;

  chisqrd = fRECO->GetChisquare();
  NDF = fRECO->GetNDF();
  redChi2 = chisqrd/NDF;
  printf("\nChisqrd = %1.3f, NDF = %d, Chisqrd/NDF = %1.3f, Prob = %1.3e\n", chisqrd, NDF, redChi2, TMath::Prob(chisqrd,NDF));

  //save alpha and n parameters if fit is
  //for integrated bins in y and pT:
  SaveCBParameters(alpha, n, fRECO->GetParError(6), fRECO->GetParError(5));


  Double_t mean2S = mean1S*(massPDG2S/massPDG1S);
  Double_t mean3S = mean1S*(massPDG3S/massPDG1S);
  Double_t sigma2S = sigma1S*(massPDG2S/massPDG1S);
  Double_t sigma3S = sigma1S*(massPDG3S/massPDG1S);

  printf("=========================================\n");
  printf("Calculate the number of Y's in the sample\n");
  printf("=========================================\n");

  TF1 *CB[kNbSpecies];
  //  Double_t intCBFit[kNbSpecies]; //integral values of a CB with N=Nfit and fixed alpha, n, sigma, width

  for(int iUps = 0; iUps < kNbSpecies; iUps++){
    sprintf(name, "CB_%d", iUps);
    CB[iUps] = new TF1(name, CBFunction, range_min, range_max, 5);
    CB[iUps]->SetParameter(0, 1.);
    if(iUps == 0){
      CB[iUps]->FixParameter(1, mean1S);
      CB[iUps]->FixParameter(2, sigma1S);
    }
    else if(iUps == 1){
      CB[iUps]->FixParameter(1, mean2S);
      CB[iUps]->FixParameter(2, sigma2S);
    }
    else if(iUps == 2){
      CB[iUps]->FixParameter(1, mean3S);
      CB[iUps]->FixParameter(2, sigma3S);
    }
    CB[iUps]->FixParameter(3, alpha);
    CB[iUps]->FixParameter(4, n);
    intCB[iUps] = CB[iUps]->Integral(range_min, range_max);
  }

  nY[UPS1S] = normY1S * intCB[UPS1S];
  nY[UPS2S] = normY2S * intCB[UPS2S];
  nY[UPS3S] = normY3S * intCB[UPS3S];

  printf("the integral of the fitted CB for the %s is: %1.3e --> #%s = %1.3e\n", specName[UPS1S], intCB[UPS1S], specName[UPS1S], nY[UPS1S]);
  printf("the integral of the fitted CB for the %s is: %1.3e --> #%s = %1.3e\n", specName[UPS2S], intCB[UPS2S], specName[UPS2S], nY[UPS2S]);
  printf("the integral of the fitted CB for the %s is: %1.3e --> #%s = %1.3e\n", specName[UPS3S], intCB[UPS3S], specName[UPS3S], nY[UPS3S]);

  //calculate the fraction of BG in a given mass interval
  massMin[UPS1S] = mean1S - nSigma*sigma1S;
  massMin[UPS2S] = mean2S - nSigma*sigma2S;
  massMin[UPS3S] = mean3S - nSigma*sigma3S;
  massMax[UPS1S] = mean1S + nSigma*sigma1S;
  massMax[UPS2S] = mean2S + nSigma*sigma2S;
  massMax[UPS3S] = mean3S + nSigma*sigma3S;
  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++){
    printf("integrating histos between %1.3f and %1.3f GeV (+- %1.1f sigma window)\n",
           massMin[iSpecies], massMax[iSpecies], nSigma);
  }

  fUps1S =  new TF1("fUps1S", CBFunction, range_min, range_max, 5);
  fUps1S->FixParameter(0, normY1S * binWidth);
  fUps1S->FixParameter(1, mean1S);
  fUps1S->FixParameter(2, sigma1S);
  fUps1S->FixParameter(3, alpha);
  fUps1S->FixParameter(4, n);

  fUps2S = new TF1("fUps2S", CBFunction, range_min, range_max, 5);
  fUps2S->FixParameter(0, normY2S * binWidth);
  fUps2S->FixParameter(1, mean2S);
  fUps2S->FixParameter(2, sigma2S);
  fUps2S->FixParameter(3, alpha);
  fUps2S->FixParameter(4, n);

  fUps3S = new TF1("fUps3S", CBFunction, range_min, range_max, 5);
  fUps3S->FixParameter(0, normY3S * binWidth);
  fUps3S->FixParameter(1, mean3S);
  fUps3S->FixParameter(2, sigma3S);
  fUps3S->FixParameter(3, alpha);
  fUps3S->FixParameter(4, n);


  //  Double_t nUps[kNbSpecies];
  //  nUps[UPS1S] = fUps1S->Integral(massMin[UPS1S], massMax[UPS1S]);
  //  nUps[UPS2S] = fUps2S->Integral(massMin[UPS2S], massMax[UPS2S]);
  //  nUps[UPS3S] = fUps3S->Integral(massMin[UPS3S], massMax[UPS3S]);

  fBG = new TF1("fBG", DrawContinuum, range_min, range_max, 3);
  for(int iPar = 0; iPar < 3; iPar++){
    fBG->FixParameter(iPar, fitParBG[iPar]);
    fBG->SetParError(iPar, fitParBGerr[iPar]);
    //    printf("fBG par %f +- %f\n", fBG->GetParameter(iPar), fBG->GetParError(iPar));

  }
  //  Double_t nBG[kNbSpecies];
  //  nBG[UPS1S] = fBG->Integral(massMin[UPS1S], massMax[UPS1S]);
  //  nBG[UPS2S] = fBG->Integral(massMin[UPS2S], massMax[UPS2S]);
  //  nBG[UPS3S] = fBG->Integral(massMin[UPS3S], massMax[UPS3S]);


  printf("1S: mass = %1.3f, sigma = %1.3f\n", mean1S, sigma1S);
  printf("1S: mass = %1.3f, sigma = %1.3f\n", mean2S, sigma2S);
  printf("3S: mass = %1.3f, sigma = %1.3f\n", mean3S, sigma3S);
  Double_t massMinSB[2], massMaxSB[2];
  massMinSB[0] = 8.6;
  massMaxSB[0] = mean1S - 3*sigma1S;
  massMinSB[1] = mean3S + 3*sigma3S;
  massMaxSB[1] = 11.4;
  printf("--> L mass window: %1.3f < M < %1.3f GeV\n", massMinSB[0], massMaxSB[0]);
  printf("--> R mass window: %1.3f < M < %1.3f GeV\n", massMinSB[1], massMaxSB[1]);




  SaveFitPars();

}

//==============================
void GetHisto(Char_t *fileNameIn){
  TFile *fin = new TFile(fileNameIn);
  Char_t name[100];
  sprintf(name, "Reco_Onia_mass");
  hMass = (TH1F*) fin->Get(name);

  hMass->Rebin(2);
  binWidth = hMass->GetBinWidth(1); //valid only for an equal bin histogram!
  printf("binwidth = %1.2e\n", binWidth);
}

//==============================
void SaveCBParameters(Double_t alpha, Double_t n, Double_t alphaErr, Double_t nErr){

  Char_t name[100];
  sprintf(name, "tmpFiles/CBpars.root");
  TFile *fOut = new TFile(name, "RECREATE");
  TTree *treeOut = new TTree("CBPars", "");
  Double_t alphaAll = alpha, nAll = n;
  Double_t alphaAllErr = alphaErr, nAllErr = nErr;
  treeOut->Branch("alphaAll", &alphaAll, "alphaAll/D");
  treeOut->Branch("nAll", &nAll, "nAll/D");
  treeOut->Branch("alphaAllErr", &alphaAllErr, "alphaAllErr/D");
  treeOut->Branch("nAllErr", &nAllErr, "nAllErr/D");
  treeOut->Fill();
  treeOut->Write();
  fOut->Close();
}

//=========================
void SaveFitPars(){

  Char_t name[100];
  sprintf(name, "tmpFiles/massFitParameters_Ups.root");
  TFile *fOut = new TFile(name, "RECREATE");
  TTree *treeOut = new TTree("massFitParameters", "");
  Int_t bufsize = 32000; //default = 32000
  Int_t splitlevel = 0; //recommended by R. Brun
  treeOut->Branch("fUps1S", "TF1", &fUps1S, bufsize, splitlevel);
  treeOut->Branch("fUps2S", "TF1", &fUps2S, bufsize, splitlevel);
  treeOut->Branch("fUps3S", "TF1", &fUps3S, bufsize, splitlevel);
  treeOut->Branch("fBG", "TF1", &fBG, bufsize, splitlevel);
  treeOut->Fill();

  treeOut->Write();
  fOut->Close();
}

//=========================
Double_t fitPolyCrystal3(Double_t *x, Double_t *par){

  Double_t normY1S = par[0];
  Double_t mean1S = par[1];
  Double_t sigma1S = par[2];
  Double_t normY2S = par[3];
  Double_t normY3S = par[4];
  Double_t n = par[5];
  Double_t alpha = par[6];
  Double_t a = par[7];
  Double_t b = par[8];
  Double_t c = par[9];

  Double_t mean2S = mean1S*(massPDG2S/massPDG1S);
  Double_t mean3S = mean1S*(massPDG3S/massPDG1S);
  Double_t sigma2S = sigma1S*(massPDG2S/massPDG1S);
  Double_t sigma3S = sigma1S*(massPDG3S/massPDG1S);

  Double_t poly2 = a + b*x[0] + c*x[0]*x[0];

  Double_t CB1 = 0.;
  if(((x[0] - mean1S)/sigma1S) > -alpha)
    CB1 = TMath::Exp(-(pow(x[0]-mean1S,2)/(2.*sigma1S*sigma1S)));
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB1 = A*pow(B - (x[0]-mean1S)/sigma1S, -n);
  }

  Double_t CB2 = 0.;
  if(((x[0] - mean2S)/sigma2S) > -alpha)
    CB2 = TMath::Exp(-(pow(x[0]-mean2S,2)/(2.*sigma2S*sigma2S)));
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB2= A*pow(B - (x[0]-mean2S)/sigma2S, -n);
  }
  Double_t CB3 = 0.;
  if(((x[0] - mean3S)/sigma3S) > -alpha)
    CB3 = TMath::Exp(-(pow(x[0]-mean3S,2)/(2.*sigma3S*sigma3S)));
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB3 = A*pow(B - (x[0]-mean3S)/sigma3S, -n);
  }

  Double_t result = poly2  + normY1S*CB1 + normY2S *CB2 + normY3S * CB3;

  result *= binWidth; //correct for the bin width
  return result;
}

//==================================
Double_t fitContinuum(Double_t *x, Double_t *par){

  if (x[0] > peak_min && x[0] < peak_max) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t result = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  result *= binWidth; //correct for the bin width
  return result;
}

//==================================
Double_t DrawContinuum(Double_t *x, Double_t *par){

  Double_t result = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  result *= binWidth; //correct for the bin width
  return result;
}
