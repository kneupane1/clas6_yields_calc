/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "branches.hpp"
#include "color.hpp"
#include "constants.hpp"
#include "reaction.hpp"

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using THnSparse_ptr = std::shared_ptr<THnSparse>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram {
 protected:
  std::shared_ptr<TFile> RootOutputFile;
  std::shared_ptr<TCanvas> def;

  // W and Q^2
  int bins = 500;
  float w_min = 0;
  float w_max = 3.25;
  float q2_min = 0;
  float q2_max = 5;

  TH2D_ptr WvsQ2_hist = std::make_shared<TH2D>("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max, bins, q2_min, q2_max);
  TH1D_ptr W_hist = std::make_shared<TH1D>("W", "W", bins, w_min, w_max);
  std::vector<TH1D_ptr> W_hist_sec;
  TH1D_ptr Q2_hist = std::make_shared<TH1D>("Q2", "Q2", bins, q2_min, w_max);

  static const short W_BIN_CHECK_NUM = 11;

  std::string W_BIN_CHECK_NAME[W_BIN_CHECK_NUM] = {" All_W_range ", " <1.4W<1.425 ", " 1.425<W<1.45 ", " 1.45<W<1.475 ",
                                                   " 1.475<W<1.5 ", " 1.5<W<1.525 ", " 1.525<W<1.55 ", " 1.55<W<1.575 ",
                                                   " 1.575<W<1.6 ", " 1.6<W<1.625 ", " 1.625<W<1.8 "};

  TH1D_ptr mm2_mPim_hist_check[11];
  TH1D_ptr W_hist_check[11];
  TH2D_ptr W_VS_Q2_hist_check[11];

  static const short q2_bin = 6;
  //   float q2_low_values[5] = {2.0, 2.40, 3.0, 3.5, 4.2};
  //   float q2_up_values[5] = {2.40, 3.0, 3.5, 4.2, 5.0};

  // int q2_bining(float q2) {
  //   // if (q2 < 2.4)
  //     return 0;
  //   // else if (q2 < 3.0)
  //   //   return 1;
  //   // else if (q2 < 3.5)
  //   //   return 2;
  //   // else if (q2 < 4.2)
  //   //   return 3;
  //   // else if (q2 < 5.0)
  //   //   return 4;
    
  // }

  static const short w_bin = 29;
  THnSparse *sevenDHist_pim[q2_bin][w_bin];
  THnSparse *sevenDHist_pip[q2_bin][w_bin];
  THnSparse *sevenDHist_prot[q2_bin][w_bin];


 public:
  Histogram(const std::string& output_file);
  ~Histogram();

  // W and Q^2
  void WvsQ2_make_hists();
  void WvsQ2_Fill(float W, float Q2, int sec);
  void WvsQ2_Write();

  void makeHistsWBinCheck();
  void Fill_W_bin_check(const std::shared_ptr<Reaction>& _e);
  void Write_WBinCheck();


  void Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e);
  void Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e);
  void Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e);

  void writeHists7D_pim();
  void writeHists7D_prot();
  void writeHists7D_pip();

  void Write();
};

#endif
