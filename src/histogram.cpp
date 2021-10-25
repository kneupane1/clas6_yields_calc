/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = std::make_shared<TCanvas>("def");
  WvsQ2_make_hists();
  makeHistsWBinCheck();
  const Int_t ndims_5D = 5;
  Int_t bins_5D[ndims_5D] = {14, 14, 10, 10, 10};
  // Mlower = mh1 + mh2
  Double_t xmin_5D[ndims_5D] = {1.1, 0.3, 0., 0., 0.};
  Double_t xmax_5D[ndims_5D] = {2.0, 1.1, 180, 360, 360};

  for (short q2 = 0; q2 < 6; q2++) {
    for (short w = 0; w < 30; w++) {
      // Mupper(W) = W − mh3

      auto name = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (2.0 + q2 * 0.5), (2.5 + q2 * 0.5),
                       (1.4 + 0.025 * w), (1.4 + 0.025 * w + 0.025));

      sevenDHist_prot[q2][w] = new THnSparseD(name, name, ndims_5D, bins_5D, xmin_5D, xmax_5D);
      sevenDHist_prot[q2][w]->Sumw2();

      sevenDHist_pip[q2][w] = new THnSparseD(name, name, ndims_5D, bins_5D, xmin_5D, xmax_5D);

      sevenDHist_pip[q2][w]->Sumw2();

      sevenDHist_pim[q2][w] = new THnSparseD(name, name, ndims_5D, bins_5D, xmin_5D, xmax_5D);

      sevenDHist_pim[q2][w]->Sumw2();

      // const Int_t ndims_5D = 7;
      // Int_t bins_5D[ndims_5D] = {29,6,14, 14, 10, 10, 10};
      // // Mlower = mh1 + mh2
      // Double_t xmin_5D[ndims_5D] = {1.4,2.0, 1.1, 0.3, 0., 0., 0.};

      //   Double_t xmax_5D[ndims_5D] = {2.125,5.0,2.0,1.1 , 180, 360, 360};

      //   auto name =
      //       Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", 2.0, 5.0, 1.4 , 2.125);

      //   sevenDHist_prot = new THnSparseD(name, name, ndims_5D, bins_5D, xmin_5D, xmax_5D);
      //   sevenDHist_prot->Sumw2();
      //   sevenDHist_pip = new THnSparseD(name, name, ndims_5D, bins_5D, xmin_5D, xmax_5D);
      //   sevenDHist_pip->Sumw2();
      //   sevenDHist_pim = new THnSparseD(name, name, ndims_5D, bins_5D, xmin_5D, xmax_5D);
      //   sevenDHist_pim->Sumw2();
    }
  }
}
Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;

  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  TDirectory* WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  WvsQ2_Write();

  std::cerr << BOLDBLUE << "WBinCheck()" << DEF << std::endl;
  TDirectory* WBinCheck_folder = RootOutputFile->mkdir("WBinCheck");
  WBinCheck_folder->cd();
  Write_WBinCheck();

  // TDirectory *THnSparse_7D_prot_folder =
  //     RootOutputFile->mkdir("THnSparse_7D_prot");
  // THnSparse_7D_prot_folder->cd();
  // writeHists7D_prot();

  // TDirectory *THnSparse_7D_pip_folder =
  //     RootOutputFile->mkdir("THnSparse_7D_pip");
  // THnSparse_7D_pip_folder->cd();
  // writeHists7D_pip();

  TDirectory* THnSparse_7D_pim_folder = RootOutputFile->mkdir("THnSparse_7D_pim");
  THnSparse_7D_pim_folder->cd();
  writeHists7D_pim();

  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::WvsQ2_make_hists() {
  W_hist_sec.reserve(6);
  for (int i = 0; i < 6; i++)
    W_hist_sec[i] = std::make_shared<TH1D>(Form("W_%d", i + 1), Form("W_%d", i + 1), bins, w_min, w_max);
}

void Histogram::WvsQ2_Fill(float W, float Q2, int sec) {
  WvsQ2_hist->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);
  if (sec > 0 && sec <= 6) {
    W_hist_sec[sec - 1]->Fill(W);
  }
}

void Histogram::WvsQ2_Write() {
  WvsQ2_hist->SetXTitle("W (GeV)");
  WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_hist->SetOption("COLZ");
  WvsQ2_hist->Write();

  W_hist->SetXTitle("W (GeV)");
  W_hist->Write();

  for (int i = 0; i < 6; i++) {
    W_hist_sec[i]->SetXTitle("W (GeV)");
    W_hist_sec[i]->Write();
  }

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist->Write();
}

void Histogram::makeHistsWBinCheck() {
  for (short i = 0; i < 11; i++) {
    mm2_mPim_hist_check[i] =
        std::make_shared<TH1D>(Form("mm2_mPim_%6.12s_", W_BIN_CHECK_NAME[i].c_str()),
                               Form("mm2_mPim_%6.12s_", W_BIN_CHECK_NAME[i].c_str()), bins, -0.1, 0.1);

    W_VS_Q2_hist_check[i] = std::make_shared<TH2D>(Form("W_VS_Q2_hist_check_%6.12s_", W_BIN_CHECK_NAME[i].c_str()),
                                                   Form("W_VS_Q2_hist_check_%6.12s_", W_BIN_CHECK_NAME[i].c_str()),
                                                   bins, 1.25, 1.85, bins, 2.9, 3.6);

    W_hist_check[i] =
        std::make_shared<TH1D>(Form("W_hist_check_%6.12s_", W_BIN_CHECK_NAME[i].c_str()),
                               Form("W_hist_check_%6.12s_", W_BIN_CHECK_NAME[i].c_str()), bins, 1.25, 1.85);
  }
}

void Histogram::Fill_W_bin_check(const std::shared_ptr<Reaction>& _e) {
  if (_e->MM_cut()) {
    if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
      W_hist_check[0]->Fill(_e->W(), _e->weight());
      mm2_mPim_hist_check[0]->Fill(_e->MM2(), _e->weight());
      W_VS_Q2_hist_check[0]->Fill(_e->W(), _e->Q2(), _e->weight());

      if (_e->W() >= 1.4 && _e->W() <= 1.425) {
        W_hist_check[1]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[1]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[1]->Fill(_e->W(), _e->Q2(), _e->weight());
      }

      else if (_e->W() >= 1.425 && _e->W() <= 1.450) {
        W_hist_check[2]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[2]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[2]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.45 && _e->W() <= 1.475) {
        W_hist_check[3]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[3]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[3]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.475 && _e->W() <= 1.50) {
        W_hist_check[4]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[4]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[4]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.5 && _e->W() <= 1.525) {
        W_hist_check[5]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[5]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[5]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.525 && _e->W() <= 1.55) {
        W_hist_check[6]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[6]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[6]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.55 && _e->W() <= 1.575) {
        W_hist_check[7]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[7]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[7]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.575 && _e->W() <= 1.60) {
        W_hist_check[8]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[8]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[8]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.6 && _e->W() <= 1.625) {
        W_hist_check[9]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[9]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[9]->Fill(_e->W(), _e->Q2(), _e->weight());
      } else if (_e->W() >= 1.625 && _e->W() <= 2.125) {
        W_hist_check[10]->Fill(_e->W(), _e->weight());
        mm2_mPim_hist_check[10]->Fill(_e->MM2(), _e->weight());
        W_VS_Q2_hist_check[10]->Fill(_e->W(), _e->Q2(), _e->weight());
      }
    }
  }
}
void Histogram::Write_WBinCheck() {
  for (short i = 0; i < 11; i++) {
    W_hist_check[i]->SetXTitle("W (GeV)");
    W_hist_check[i]->Write();

    mm2_mPim_hist_check[i]->SetXTitle("MM^{2} (GeV^{2})");
    mm2_mPim_hist_check[i]->Write();

    // W_VS_Q2_hist_check[i]->SetYTitle("Q^{2} (GeV^{2})");
    // W_VS_Q2_hist_check[i]->SetXTitle("W (GeV)");
    // W_VS_Q2_hist_check[i]->SetOption("COLZ");
    // W_VS_Q2_hist_check[i]->Write();
  }
}

void Histogram::Fill_histSevenD_prot(const std::shared_ptr<Reaction>& _e) {
  // fill it
  const Int_t ndims = 5;
  Double_t x[ndims];
  // x[0] = _e->W();
  // x[1] = _e->Q2();
  x[0] = _e->inv_Ppip();
  x[1] = _e->inv_pip_pim();
  x[2] = _e->prot_theta();
  x[3] = _e->prot_Phi();
  x[4] = _e->alpha_pippim_pipf();
  if (_e->MM_cut()) {  // abs(mmsq<0.03)
    if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
      {
        TThread::Lock();
        sevenDHist_prot[int((_e->Q2() - 2.0) / 0.5)][int((_e->W() - 1.4) / 0.025)]->Fill(x, 1);
        TThread::UnLock();
        sevenDHist_prot[int((_e->Q2() - 2.0) / 0.5)][int((_e->W() - 1.4) / 0.025)]->GetNbins();
      }
    }
  }
}

void Histogram::writeHists7D_prot() {
  for (short q2 = 0; q2 < 6; q2++) {
    for (short w = 0; w < 30; w++) {
      sevenDHist_prot[q2][w]->Write();
    }
  }
 }

    void Histogram::Fill_histSevenD_pip(const std::shared_ptr<Reaction>& _e) {
      // fill it
      const Int_t ndims = 5;
      Double_t x[ndims];
      // x[0] = _e->W();
      // x[1] = _e->Q2();
      x[0] = _e->inv_Ppim();
      x[1] = _e->inv_pip_pim();
      x[2] = _e->pip_theta();
      x[3] = _e->pip_Phi();
      x[4] = _e->alpha_ppim_pipip();
      if (_e->MM_cut()) {  // abs(mmsq<0.03)
        if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
          TThread::Lock();
          sevenDHist_pip[int((_e->Q2() - 2.0) / 0.5)][int((_e->W() - 1.4) / 0.025)]->Fill(x, 1);
          TThread::UnLock();
          sevenDHist_pip[int((_e->Q2() - 2.0) / 0.5)][int((_e->W() - 1.4) / 0.025)]->GetNbins();
        }
      }
    }
    void Histogram::writeHists7D_pip() {
      for (short q2 = 0; q2 < 6; q2++) {
        for (short w = 0; w < 30; w++) {
          sevenDHist_pip[q2][w]->Write();
        }
      }
}
void Histogram::Fill_histSevenD_pim(const std::shared_ptr<Reaction>& _e) {
  // fill it
  const Int_t ndims = 5;
  Double_t x[ndims];
  // x[0] = _e->W();
  // x[1] = _e->Q2();
  x[0] = _e->inv_Ppip();
  x[1] = _e->inv_pip_pim();
  x[2] = _e->pim_theta();
  x[3] = _e->pim_Phi();
  x[4] = _e->alpha_ppip_pipim();
  if (_e->MM_cut()) {  // abs(mmsq<0.03)
    if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
      TThread::Lock();
      sevenDHist_pim[int((_e->Q2() - 2.0) / 0.5)][int((_e->W() - 1.4) / 0.025)]->Fill(x, 1.0);
      TThread::UnLock();
      sevenDHist_pim[int((_e->Q2() - 2.0) / 0.5)][int((_e->W() - 1.4) / 0.025)]->GetNbins();
    }
  }
}
void Histogram::writeHists7D_pim() {
  for (short q2 = 0; q2 < 6; q2++) {
    for (short w = 0; w < 30; w++) {
      sevenDHist_pim[q2][w]->Write();
    }
  }
}
// void Histogram::Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e) {
//   // fill it
//   const Int_t ndims = 7;
//   Double_t x[ndims];
//   x[0] = _e->W();
//    x[1] = _e->Q2();
//   x[2] = _e->inv_Ppip();
//   x[3] = _e->inv_pip_pim();
//   x[4] = _e->prot_theta();
//   x[5] = _e->prot_Phi();
//   x[6] = _e->alpha_pippim_pipf();
//   if (_e->MM_cut()) {  // abs(mmsq<0.03)
//     if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
//       {
//         TThread::Lock();
//         sevenDHist_prot->Fill(x, _e->weight());
//         TThread::UnLock();
//         sevenDHist_prot->GetNbins();
//       }
//     }
//   }
// }

// void Histogram::writeHists7D_prot() {
//       sevenDHist_prot->Write();
// }

// void Histogram::Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e) {
//   // fill it
//   const Int_t ndims = 7;
//   Double_t x[ndims];
//   x[0] = _e->W();
//    x[1] = _e->Q2();
//   x[2] = _e->inv_Ppim();
//   x[3] = _e->inv_pip_pim();
//   x[4] = _e->pip_theta();
//   x[5] = _e->pip_Phi();
//   x[6] = _e->alpha_ppim_pipip();
//   if (_e->MM_cut()) {  // abs(mmsq<0.03)
//     if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
//       TThread::Lock();
//       sevenDHist_pip->Fill(x, _e->weight());
//       TThread::UnLock();
//       sevenDHist_pip->GetNbins();
//     }
//   }
// }
// void Histogram::writeHists7D_pip() {

//       sevenDHist_pip->Write();

// }

// void Histogram::Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e) {
//   // fill it
//   const Int_t ndims = 7;
//   Double_t x[ndims];
//   x[0] = _e->W();
//    x[1] = _e->Q2();
//   x[2] = _e->inv_Ppip();
//   x[3] = _e->inv_pip_pim();
//   x[4] = _e->pim_theta();
//   x[5] = _e->pim_Phi();
//   x[6] = _e->alpha_ppip_pipim();
//   if (_e->MM_cut()) {  // abs(mmsq<0.03)
//     if (_e->W() <= 2.125 && _e->W() >= 1.4 && _e->Q2() >= 2.0 && _e->Q2() <= 5.0) {
//       TThread::Lock();
//       sevenDHist_pim->Fill(x, _e->weight());
//        TThread::UnLock();
//        sevenDHist_pim->GetNbins();
//     }
//   }
// }
// void Histogram::writeHists7D_pim() {

//       sevenDHist_pim->Write();

// }
