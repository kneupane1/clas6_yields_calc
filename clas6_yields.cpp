#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "constants.hpp"
#include <mutex>
#include "TMath.h"

static const int W_bins_no = 64;

float q2_low_values[11] = {1.0, 2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.2, 7.4, 8.6, 9.8};
float q2_up_values[11] = {2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.2, 7.4, 8.6, 9.8, 12.0};
float q2_binning_size[11] = {1.0, 0.4, 0.6, 0.5, 0.7, 0.8, 1.2, 1.2, 1.2, 1.2, 2.2};

THnSparseD *h_exp_prot[W_bins_no];
THnSparseD *h_exp_pip[W_bins_no];
THnSparseD *h_exp_pim[W_bins_no];

float w_bin_size = 0.05;

const float E_beam =5.754;
int n = 64;
// TCanvas *can1 =
//     new TCanvas("Nine_differential_Yields", "Nine_diff_Yields", 1900, 1500);


//The subroutine calculates the virtual photon flux
float flux(float w, float q2_value)
{
    // float w = 1.0 + w_bin * 0.025;
    // float q2 = 1.0 + q2_bin * 1.0;
    // float q2_value = q2_low_values[q2_bin];
    // std::cout << "q2 low in flux function : " << q2_value << std::endl;
    // std::cout << "w low in flux function : " << w << std::endl;

    float omega = (w * w + q2_value - MASS_P * MASS_P) / (2. * MASS_P);
    // std::cout << "omega =  " << omega << std::endl;

    float en_elp = E_beam - omega;
    float th_elp = 2 * asin(sqrt(q2_value / 4. / E_beam / en_elp));
    // std::cout << "theta =  " << th_elp << std::endl;

    float epsilon = 1 / (1. + 2. * (1. + omega * omega / q2_value) * (tan(th_elp / 2.)) * (tan(th_elp / 2.)));

    std::cout << "epsilon = : " << epsilon << std::endl;

    float flux_calc = (omega - q2_value / 2. / MASS_P) / 137.;

    flux_calc = flux_calc / 2. / (PI) / E_beam / q2_value / (1 - epsilon);
    flux_calc = flux_calc * w / E_beam / MASS_P;
    // std::cout << "flux  =  " << flux_calc << std::endl;
    // std::cout << "pi  =  " << PI << std::endl;

    return flux_calc;
}

int main(int argc, char **argv)
{
  std::string inFileName = "/Users/krishnaneupane/Desktop/2021/oct_2021/e16_yield_tests.root";

  // sept_2021/sim_gemc_4_4_all_data_weight_1_09_22_2021.root"; // without weight

  gStyle->SetHistMinimumZero();
  TFile *root_data = new TFile(inFileName.c_str());

  // define and set style
  gStyle->SetTitleSize(0.07, "t");
  gStyle->SetOptStat(0);
  // gStyle->SetStatColor(0);
  // gStyle->SetPaperSize(18, 24);
  gStyle->SetLabelSize(0.05, "XY");
  gStyle->SetTitleOffset(0.70, "X");
  gStyle->SetTitleOffset(0.9, "Y");
  gStyle->SetTitleSize(0.06, "XY");
  // gStyle->SetStatFontSize(0.07);
  // gStyle->SetTitleFont(62, "XY");
  gStyle->SetTitleFontSize(0.06);
  // gStyle->SetLabelFont(60, "XY");
  // gStyle->SetTextFont(72);
  gStyle->SetLegendFont(72);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.34);
  // gStyle->SetMarkerColor(4);
  gStyle->SetErrorX(0);

  // TFile *fiveD = new TFile("fiveD.root","RECREATE");
  for (short q2 = 3; q2 < 4; q2++) {
    float q2_lower_lim = q2_low_values[q2];
    float q2_upper_lim = q2_up_values[q2];
    float q2_bin_size = q2_binning_size[q2];
    float q2_mid_value = (q2_low_values[q2] + q2_up_values[q2]) / 2;

    std::cout << " q2_bin_size " << q2_bin_size << "\n";
    std::cout << " q2 values for flux " << q2_mid_value << "\n";

    for (short w = 10; w < 11; w++) {
      // if (w == 19 || w == 20)
      // continue;
      float w_bin_for_flux = 1.0 + w * 0.05 + 0.025;
      std::cout << " w low " << (1.0 + 0.05 * w) << "\n";
      std::cout << " w bin size " << w_bin_size << "\n";
      std::cout << " w_bin for flux " << w_bin_for_flux << "\n";

      TCanvas *can1 = new TCanvas("Nine_differential_Yields", "Nine_diff_Yields", 1900, 1500);
      can1->Divide(3, 3);

      auto name = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim), (1.0 + 0.05 * w),
                       (1.0 + 0.05 * w + 0.05));

      auto output_name = Form("nine_1D_cs_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim),
                              (1.0 + 0.05 * w), (1.0 + 0.05 * w + 0.05));
      // auto output_name = Form("nine-1D_cs_%.1f<=Q2<=%.1f GeV2_%.3f<=W<=%.3f GeV", (q2_lower_lim), (q2_upper_lim),
      // (1.0 + 0.025 * w), (1.0 + 0.025 * w + 0.025));

      int no_w_bins = w;

      // From THnsparse Prot
      h_exp_prot[no_w_bins] = (THnSparseD *)root_data->Get(Form("THnSparse_7D_prot/%s", name));

      // From THnsparse pip

      h_exp_pip[no_w_bins] = (THnSparseD *)root_data->Get(Form("THnSparse_7D_pip/%s", name));

      // From THnsparse pim

      h_exp_pim[no_w_bins] = (THnSparseD *)root_data->Get(Form("THnSparse_7D_pim/%s", name));

      {
        // 1. faraday cup charge = 0.038 ..C
        h_exp_prot[no_w_bins]->Scale(1. / (0.02129/3));

        // 3, 4 Normalization for luminosity, flux
        h_exp_prot[no_w_bins]->Scale(7.55314965e-11);
        h_exp_prot[no_w_bins]->Scale(1. / flux(w_bin_for_flux, q2_mid_value));

        std::cout << " flux : " << flux(w_bin_for_flux, q2_mid_value) << "\n";

        // 5 Projection and Normalization : 3 factors
        TH1D *yield_inv_Ppip = h_exp_prot[no_w_bins]->Projection(0);
        TH1D *yield_inv_pip_pim = h_exp_prot[no_w_bins]->Projection(1);
        TH1D *yield_theta_prot = h_exp_prot[no_w_bins]->Projection(2);
        TH1D *yield_alpha_prot = h_exp_prot[no_w_bins]->Projection(4);

        Float_t factor_inv_mass = yield_inv_Ppip->GetBinWidth(5);
        std::cout << "scale factor inv mass Ppip: " << factor_inv_mass << std::endl;
        Float_t factor_alpha_angle = PI * yield_alpha_prot->GetBinWidth(5) / 180;
        std::cout << "scale factor_alpha_ angle prot: " << factor_alpha_angle << std::endl;

        yield_inv_Ppip->Scale(1. / factor_inv_mass / (w_bin_size * q2_bin_size));
        yield_inv_pip_pim->Scale(1. / factor_inv_mass / (w_bin_size * q2_bin_size));
        yield_alpha_prot->Scale(1. / factor_alpha_angle / (w_bin_size * q2_bin_size));

        TH1D *hnew = new TH1D("hnew", "#Theta_{p'}", 10, 0, 180.);
        float cosine_value;
        for (int i = 1; i <= yield_theta_prot->GetXaxis()->GetNbins(); i++)
        {
            cosine_value = -1 / cos(yield_theta_prot->GetBinContent(i));
            hnew->SetBinContent(i, abs(cosine_value * yield_theta_prot->GetBinContent(i)));
        }
        TH1D *h_cos_th_prot;
        Double_t temp_prot;
        //Int_t n_theta_bins;
        h_cos_th_prot = new TH1D("h_cos_th_prot", "h_cos_th_prot", yield_theta_prot->GetXaxis()->GetNbins(), 0.,
        180.); for (Int_t j = 1; j <= yield_theta_prot->GetXaxis()->GetNbins(); j++)
        {
            temp_prot = cos((yield_theta_prot->GetBinLowEdge(j)) * PI / 180.) - cos(PI / 180. *
            (yield_theta_prot->GetBinLowEdge(j) + yield_theta_prot->GetBinWidth(j))); h_cos_th_prot->SetBinContent(j,
            temp_prot); h_cos_th_prot->SetBinError(j, 0.);
        }
        yield_theta_prot->Divide(h_cos_th_prot);
        yield_theta_prot->Scale(1. / (w_bin_size * q2_bin_size));

        can1->cd(1);
        yield_inv_Ppip->SetMarkerColor(7);
        yield_inv_Ppip->SetLineColor(7);
        yield_inv_Ppip->SetTitle("M(#pi^{+}p')");
        yield_inv_Ppip->SetMinimum(0.);
        yield_inv_Ppip->Draw("Z ");

        can1->cd(3);

        yield_inv_pip_pim->SetMarkerColor(7);
        yield_inv_pip_pim->SetLineColor(7);

        yield_inv_pip_pim->SetTitle("M(#pi^{+}#pi^{-})");
        yield_inv_pip_pim->SetXTitle("M#pi^{+}#pi^{-} (GeV)");
        yield_inv_pip_pim->SetMinimum(0.);
        yield_inv_pip_pim->Draw("Z");
        yield_inv_pip_pim->SetYTitle("Yield");

        can1->cd(6);
        yield_theta_prot->SetMarkerColor(7);
        yield_theta_prot->SetLineColor(7);
        yield_theta_prot->SetTitle("#Theta_{p'}");
        yield_theta_prot->SetXTitle("#Theta_{p'}(deg)");

        yield_theta_prot->SetMinimum(0.);
        yield_theta_prot->Draw("Z");
        yield_theta_prot->SetYTitle("Yield/d(-cos#theta)");  //("Yield/d(-cos#theta)");

        can1->cd(9);
        yield_alpha_prot->SetMarkerColor(4);
        yield_alpha_prot->SetMarkerColor(7);
        yield_alpha_prot->SetLineColor(7);
        yield_alpha_prot->SetTitle("#alpha_{p'}");
        yield_alpha_prot->SetXTitle("#alpha_{p'}(deg)");
        yield_alpha_prot->Draw("Z");
        yield_alpha_prot->SetYTitle("Yield");
        yield_alpha_prot->Draw("Z same");
      }
      {
        h_exp_pip[no_w_bins]->Scale(1. / (0.02129/3));

        // here i am trying to scale for luminosity
        h_exp_pip[no_w_bins]->Scale(7.55314965e-11);
        h_exp_pip[no_w_bins]->Scale(1. / flux(w_bin_for_flux, q2_mid_value));
        TH1D *yield_inv_Ppim = h_exp_pip[no_w_bins]->Projection(0);
        TH1D *yield_theta_pip = h_exp_pip[no_w_bins]->Projection(2);
        TH1D *yield_alpha_pip = h_exp_pip[no_w_bins]->Projection(4);

        Float_t factor_inv_mass = yield_inv_Ppim->GetBinWidth(5);
        std::cout << "scale factor inv mass: " << factor_inv_mass << std::endl;
        Float_t factor_alpha_angle = PI * yield_alpha_pip->GetBinWidth(5) / 180;
        std::cout << "scale factor_alpha_ angle pip : " << factor_alpha_angle << std::endl;

        TH1D *h_cos_th_pip;
        Double_t temp_pip;
        // Int_t n_theta_bins;
        h_cos_th_pip = new TH1D("h_cos_th_pip", "h_cos_th_pip",
                                yield_theta_pip->GetXaxis()->GetNbins(), 0., 180.);
        for (Int_t j = 1; j <= yield_theta_pip->GetXaxis()->GetNbins(); j++) {
          temp_pip = cos((yield_theta_pip->GetBinLowEdge(j)) * PI / 180.) -
                     cos(PI / 180. *
                         (yield_theta_pip->GetBinLowEdge(j) +
                          yield_theta_pip->GetBinWidth(j)));
          h_cos_th_pip->SetBinContent(j, temp_pip);
          h_cos_th_pip->SetBinError(j, 0.);
        }
        yield_inv_Ppim->Scale(1. / (w_bin_size * q2_bin_size));
        yield_inv_Ppim->Scale(1. / factor_inv_mass);

        yield_theta_pip->Divide(h_cos_th_pip);
        yield_theta_pip->Scale(1. / (w_bin_size * q2_bin_size));

        yield_alpha_pip->Scale(1. / factor_alpha_angle / (w_bin_size * q2_bin_size));

        can1->cd(2);
        yield_inv_Ppim->SetMarkerColor(7);
        yield_inv_Ppim->SetLineColor(7);

        yield_inv_Ppim->SetTitle("M(#pi^{-}p')");
        yield_inv_Ppim->SetStats(0);
        yield_inv_Ppim->SetXTitle("M#pi^{-}p'(GeV)");
        yield_inv_Ppim->Draw("Z");
        yield_inv_Ppim->SetYTitle("Yield");

        can1->cd(5);
        yield_theta_pip->SetMarkerColor(7);
        yield_theta_pip->SetLineColor(7);

        yield_theta_pip->SetTitle("#Theta_{#pi^{+}}");
        yield_theta_pip->SetStats(0);
        yield_theta_pip->SetXTitle("#Theta_{#pi^{+}}(deg)");
        yield_theta_pip->Draw("Z");
        yield_theta_pip->SetYTitle("Yield/d(-cos#theta)");

        can1->cd(8);
        yield_alpha_pip->SetMarkerColor(7);
        yield_alpha_pip->SetLineColor(7);
        yield_alpha_pip->SetTitle("#alpha_{#pi^{+}}");
        yield_alpha_pip->SetStats(0);
        yield_alpha_pip->SetXTitle(" #alpha_{#pi^{+}}(deg)");
        yield_alpha_pip->Draw("Z");
        yield_alpha_pip->SetYTitle("Yield");
      }
      {
        // faraday cup charge = 0.021 ..C

        // L = 28.18
        // 1/L = 0.035486 * 10 ^39 = 3.5486 * 10 ^-11 microbarn

        // or L = 1.3419 without including charge
        // 1/L = 0.755 
        h_exp_pim[no_w_bins]->Scale(1. / (0.02129/3));

        // here i am trying to scale for luminosity
        h_exp_pim[no_w_bins]->Scale(7.55314965e-11);
        h_exp_pim[no_w_bins]->Scale(1. / flux(w_bin_for_flux, q2_mid_value));

        TH1D *yield_theta_pim = h_exp_pim[no_w_bins]->Projection(2);
        TH1D *yield_alpha_pim = h_exp_pim[no_w_bins]->Projection(4);

        Float_t factor_alpha_angle = PI * yield_alpha_pim->GetBinWidth(5) / 180;
        std::cout << "scale factor_alpha_ angle pim : " << factor_alpha_angle << std::endl;

        TH1D *h_cos_th_pim;
        Double_t temp_pim;
        // Int_t n_theta_bins;
        h_cos_th_pim = new TH1D("h_cos_th_pim", "h_cos_th_pim",
                                yield_theta_pim->GetXaxis()->GetNbins(), 0., 180.);
        for (Int_t j = 1; j <= yield_theta_pim->GetXaxis()->GetNbins(); j++) {
          temp_pim = cos((yield_theta_pim->GetBinLowEdge(j)) * PI / 180.) -
                     cos(PI / 180. *
                         (yield_theta_pim->GetBinLowEdge(j) +
                          yield_theta_pim->GetBinWidth(j)));
          h_cos_th_pim->SetBinContent(j, temp_pim);
          h_cos_th_pim->SetBinError(j, 0.);
        }
        yield_theta_pim->Divide(h_cos_th_pim);
        yield_theta_pim->Scale(1. / (w_bin_size * q2_bin_size));

        yield_alpha_pim->Scale(1. / factor_alpha_angle / (w_bin_size * q2_bin_size));

        can1->cd(4);
        yield_theta_pim->SetMarkerColor(7);
        yield_theta_pim->SetLineColor(7);
        yield_theta_pim->SetTitle("#Theta_{#pi^{-}}");
        yield_theta_pim->SetStats(0);
        yield_theta_pim->SetXTitle("#Theta_{#pi^{-}}(deg)");
        yield_theta_pim->SetYTitle(" Yield");
        yield_theta_pim->Draw("Z");
        yield_theta_pim->SetYTitle("Yield/d(-cos#theta)");

        can1->cd(7);
        yield_alpha_pim->SetMarkerColor(7);
        yield_alpha_pim->SetLineColor(7);
        yield_alpha_pim->SetTitle("#alpha_{#pi^{-}}");
        yield_alpha_pim->SetStats(0);
        yield_alpha_pim->SetXTitle(" #alpha_{#pi^{-}}(deg)");
        yield_alpha_pim->Draw("Z");
        yield_alpha_pim->SetYTitle("Yield");
      }
      // Int[i] = (Int_1[i] + Int_2[i] + Int_3[i]) / 3.;
      // // h_w_int->Fill(1.0 + w * 0.05 + 0.025, Int[w]);

      // // TDirectory *output_png = RootOutputFile->mkdir("output_png");
      // // output_png->cd();
      can1->SaveAs(Form("mPim_clas6_yield_%s.png", output_name));
    }
    std::cout << "In_   " << q2_mid_value << "\n";

    // for (short w = 6; w < 64; w++)
    // {
    //         // std::cout << std::setprecision(3);
    //         // std::cout << estimate1[w] << "\n";

    //         Int[w] = (Int_1[w] + Int_2[w] + Int_3[w]) / 3.;
    //         // // // w_for_int[w] = 1.0 + w * 0.05 + 0.025;
    //         std::cout << Int[w] << "\n";
    //         // // std::cout << " w val " << w_for_int[w] << "\n";

    //         // // h_w_int->Fill(1.0 + w * 0.05 + 0.025, Int[w]);

    //         // // h_w_int->Draw("e1pX same");
    // }
    // std::cout << " w val "
    //           << "\n";

    // for (short w = 6; w < 64; w++)
    // {

    //         w_for_int[w] = 1.0 + w * 0.05 + 0.025;
    //         std::cout << w_for_int[w] << "\n";

    //         // h_w_int->Fill(1.0 + w * 0.05 + 0.025, Int[w]);

    //         // h_w_int->Draw("e1pX same");
    // }
    std::cout << " q2 values for flux " << q2_mid_value << "\n";

    // TGraph *gr1 = new TGraph(n, w_for_int, Int);
    // // gr1->SetName("gr1");
    // gr1->SetTitle("Integrated Yields");
    // gr1->SetMarkerColor(1);
    // gr1->SetMarkerStyle(20);
    // gr1->SetDrawOption("AP");
    // gr1->SetFillStyle(0);
    // TAxis *axis = gr1->GetXaxis();
    // axis->SetLimits(2, 3.0); // along X
    // // gr1->GetYaxis()->SetRangeUser(0.2, 300.7);
    // // gr1->SetMarkerSize(2.5);
    // // gr1->SetTitle("1/R vs W");
    // gr1->GetXaxis()->SetTitle("W (GeV)");
    // gr1->GetYaxis()->SetTitle(" Int Yield");
    // gr1->Draw("AP");
    // mg->Draw("P");
    // can1->Print("multigraphleg2.png");

    // can1->SaveAs(Form("%s.png", name_int));
    }
    // std::string outfilename;
    //
    // std::shared_ptr<TFile> RootOutputFile;
    // RootOutputFile = std::make_shared<TFile>(outfilename, "RECREATE");

    // can1->SaveAs(Form("%s.png",can1->GetName()));

    // return 0;
}

