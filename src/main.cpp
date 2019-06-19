#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <math.h>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include "RadSolver.hpp"
#include "isrutils.hpp"

namespace po = boost::program_options;

typedef struct {
  double thsd;
  std::string gname;
  std::string fcnname;
  std::string ifname;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help",
                      "A simple tool, designed to find numerical"
                      "solution of the Kuraev-Fadin equation.")(
      "thsd", po::value<double>(&(opts->thsd)), "Threshold (GeV).")(
      "gname", po::value<std::string>(&(opts->gname)),
      "Name of the measured cross section graph.")(
      "fcnname", po::value<std::string>(&(opts->fcnname)),
      "Name of the Born cross section function.")(
      "ifname", po::value<std::string>(&(opts->ifname)), "Path to input file.")(
      "ofname", po::value<std::string>(&(opts->ofname)),
      "Path to output file.");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options:");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  if (vmap.count("thsd") && vmap.count("gname") && vmap.count("ifname") &&
      vmap.count("ofname")) {
    TMatrixT<double> invErrM;
    TMatrixT<double> integralOperatorMatrix;
    auto fl0 = TFile::Open(opts.ifname.c_str(), "read");
    auto gr = reinterpret_cast<TGraphErrors*>(fl0->Get(opts.gname.c_str()));
    RadSolver rs(gr, opts.thsd * opts.thsd);
    fl0->Close();
    auto gborn = rs.getBornCS(invErrM, integralOperatorMatrix);
    auto fl1 = TFile::Open(opts.ofname.c_str(), "recreate");
    fl1->cd();
    rs.visible_cs->Write("visible");
    gborn->Write("cborn");
    invErrM.Write("invBornErrMatrix");
    integralOperatorMatrix.Write("integralOperatorMatrix");
    fl1->Close();
  } else if (vmap.count("fcnname") && vmap.count("ifname") &&
             vmap.count("ofname")) {
    TGraph* rad;
    TFile* fl = TFile::Open(opts.ifname.c_str(), "read");
    TF1* bfcn =
        dynamic_cast<TF1*>(fl->Get(opts.fcnname.c_str())->Clone("bfcn"));
    fl->Close();
    if (bfcn) {
      int n = 1000;
      double a = 1.e-3 * bfcn->GetXaxis()->GetXmin();
      double b = 1.e-3 * bfcn->GetXaxis()->GetXmax();
      double h = (b - a) / n;
      double min_s = (a + h) * (a + h);
      TFile* fout = TFile::Open(opts.ofname.c_str(), "recreate");
      fout->cd();
      rad = new TGraph();
      rad->SetName("radcorr");
      int c = 0;
      double en;
      double s;
      double bcs;
      std::function<double(double)> bornfcn = [bfcn, min_s](double s) {
        if (s <= min_s) {
          return 0.;
        }
        return bfcn->Eval(1.e+3 * sqrt(s));
      };
      for (int i = 1; i < n; ++i) {
        en = a + h * i;
        s = en * en;
        bcs = bornfcn(s);
        if (bcs > 0) {
          rad->SetPoint(c, en,
                        integral_KuraevFadin(s, bornfcn, 1. - min_s / s) / bcs);
          c++;
        }
      }
      rad->Write();
      fout->Close();
      delete rad;
    } else {
      std::cout << "The fcnname flag doesn't correspond to any TF1 function."
                << std::endl;
      return 1;
    }

  } else {
    help(desc);
  }
  return 0;
}
