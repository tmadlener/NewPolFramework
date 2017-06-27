#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgList.h"
// #include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooCBShape.h"

#include <string>
// #include <vector>
#include <iostream>
#include <array>
#include <sstream>
#include <regex>

void buildModel(RooWorkspace* ws)
{
  RooRealVar* mass = ws->var("mass");

  mass->setRange("fitRange", 8.6, 11.4);

  RooRealVar a0("a0", "a0", 0.5, -1, 1); a0.setConstant(false);
  RooRealVar a1("a1", "a1", 0, -1, 1); a1.setConstant(false);
  RooRealVar a2("a2", "a2", 0, -1, 1); a2.setConstant(false);
  RooChebychev bkgPoly("bkgPoly", "polynomial background", *mass, RooArgList(a0, a1, a2));

  constexpr double mPDG1S = 9.460;
  constexpr double mPDG2S = 10.023;
  constexpr double mPDG3S = 10.355;

  RooRealVar r2S1S("r2S1S", "r2S1S", mPDG2S / mPDG1S);
  RooRealVar r3S1S("r3S1S", "r3S1S", mPDG3S / mPDG1S);

  RooRealVar mean1S("mean1S", "mean1S", mPDG1S, 8.6, 11.4);
  RooRealVar sigma1S("sigma1S", "sigma1S", 0.1, 0, 2.5);
  RooRealVar alpha("alpha", "alpha", 1.33, 0, 2.5);
  RooRealVar n("n", "n", 6.6, 0, 10);

  RooFormulaVar mean2S("mean2S", "mean2S", "mean1S * r2S1S", RooArgList(mean1S, r2S1S));
  RooFormulaVar sigma2S("sigma2S", "sigma2S", "sigma1S * r2S1S", RooArgList(sigma1S, r2S1S));

  RooFormulaVar mean3S("mean3S", "mean3S", "mean1S * r3S1S", RooArgList(mean1S, r3S1S));
  RooFormulaVar sigma3S("sigma3S", "sigma3S", "sigma1S * r3S1S", RooArgList(sigma1S, r3S1S));

  RooCBShape sigCB1S("sigCB1S", "sigCB1S", *mass, mean1S, sigma1S, alpha, n);
  RooCBShape sigCB2S("sigCB2S", "sigCB2S", *mass, mean2S, sigma2S, alpha, n);
  RooCBShape sigCB3S("sigCB3S", "sigCB3S", *mass, mean3S, sigma3S, alpha, n);

  ws->import(RooArgList(bkgPoly, sigCB1S, sigCB2S, sigCB3S));
  ws->factory("SUM:fullModel(fBkg[0.5,0,1] * bkgPoly, f1S[0.2,0,1]*sigCB1S, f2S[0.15,0,1]*sigCB2S, sigCB3S)");
}

void plotModel(RooWorkspace* ws, const std::string& snapshot)
{
  using namespace RooFit;

  auto* mass = ws->var("mass");
  auto* frame = mass->frame(Range("fitRange"));
  auto* data = ws->data("fullData");
  auto* fullModel = ws->pdf("fullModel");

  ws->loadSnapshot("snap_fullData");

  data->plotOn(frame);
  fullModel->plotOn(frame);

  auto* can = new TCanvas("c", "c", 1000, 1000);
  can->cd();
  frame->Draw();

  can->SaveAs("fitResults.pdf");
}

template<typename V>
std::string getBinExpr(const V& binning, const size_t bin, const std::string& var)
{
  std::stringstream sstr;
  sstr << "(" << var << " > " << binning[bin - 1] << " && " << var << " < " << binning[bin] << ")";
  return sstr.str();
}

template<typename V>
std::string getBinName(const V& binning, const size_t bin, const std::string& var)
{
  std::stringstream sstr;
  sstr << var << "_" << binning[bin - 1] << "to" << binning[bin];

  return std::regex_replace(sstr.str(), std::regex("([0-9]+).([0-9]+)"), "$1p$2");
}

void costhBinFits(RooWorkspace* ws, const std::string&& fullDataName)
{
  using namespace RooFit;

  constexpr std::array<double, 10> absCosThEdges = {0.0, 0.1, 0.2, 0.3, 0.4,
                                                    0.5, 0.6, 0.7, 0.8, 1.0};

  auto* fullData = ws->data(fullDataName.c_str());
  auto* model = ws->pdf("fullModel");
  auto* params = (RooArgSet*) model->getParameters(*(ws->var("mass")));

  std::cout << fullData << " " << model << " " << params << "\n";

  for (size_t i = 1; i < absCosThEdges.size(); ++i) {
    const auto cutString = getBinExpr(absCosThEdges, i, "TMath::Abs(costh_HX)");
    const auto binName = getBinName(absCosThEdges, i, "absCosth");

    auto* binData = fullData->reduce(cutString.c_str());
    binData->SetName(("data_" + binName).c_str());
    ws->import(*binData);

    auto* rlt = model->fitTo(*binData, Minos(false), NumCPU(4), Range("fitRange"),
                             Save(true));

    rlt->SetName(("fitResults_" + binName).c_str());
    ws->import(*rlt);

    ws->saveSnapshot(("snap_" + binName).c_str(), *params, true);

    // delete binData;
  }

  // delete fullData;
}

void massFits_costh(std::string fn)
{
  using namespace RooFit;

  TFile* f = TFile::Open(fn.c_str());
  TTree* t = static_cast<TTree*>(f->Get("selectedData"));

  RooRealVar pT("pT", "p_{T}", 10, 70);
  RooRealVar mass("mass", "m_{B}", 8.4, 11.6);
  RooRealVar Nch("Nch", "Nch", 0, 180);
  RooRealVar costh("costh_HX", "cos#theta^{HX}", -1, 1);
  RooRealVar phi("phi_HX", "phi^{HX}", -180, 180);
  RooRealVar ctau("ctau", "c#tau", -40, 40);
  RooRealVar ctauErr("ctauErr", "#sigma_{c#tau}", 0, 5);

  RooDataSet fullData("fullData", "dataset without cuts", t,
                      RooArgList(pT, mass, Nch, costh, phi, ctau, ctauErr));

  RooWorkspace* ws = new RooWorkspace("workspace", "workspace");
  ws->import(fullData);

  buildModel(ws);

  auto* model = ws->pdf("fullModel");
  auto* params = (RooArgSet*) model->getParameters(mass);


  auto* fitData = dynamic_cast<RooDataSet*>(fullData.reduce("(TMath::Abs(ctau) / ctauErr) < 2.0 && pT > 20.0"));
  fitData->SetName("fitData");
  ws->import(*fitData);

  RooFitResult* rlt = model->fitTo(*fitData, Minos(false), NumCPU(4), Range("fitRange"),
                                   Save(true));

  ws->saveSnapshot("snap_fullData", *params, true);
  // std::cout << "First Fit done ==================================================\n";
  // RooFitResult* rlt = ws->pdf("fullModel")->fitTo(fullData, NumCPU(4), Save(true),
  //                                                 Minos(true), Range("fitRange"));

  std::cout << rlt->status() << " " << rlt->covQual() << "\n";
  ws->import(*rlt);

  costhBinFits(ws, "fitData");
  ws->writeToFile("ws_fit_result_ctau2p0_pT20.root");

  // plotModel(ws, "snap_fullData");
}

#if !(defined(__CINT__) || defined(__CLING__))
int main(int argc, char *argv[])
{
  std::string filename = argv[1];

  massFits_costh(filename);

  return 0;
}
#endif
