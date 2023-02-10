
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include "dlUtility.h"
#include "Utility/TreeSetting.h"
#include "TH1.h"
#include "Fit_Funcs.h"

using namespace std;
double GetFWHM(TF1 *f){
  double start = 0.0;
  double end = 12000.0;
  double x1, x2;
  double halfmax = f->GetMaximum()/2.;
  double max_x = f->GetMaximumX();
  if (f->Eval(start) < halfmax){
    x1 = f->GetX(halfmax, start, max_x);
  }
  else x1 = f->GetMinimum(start, max_x);

  if (f->Eval(end) < halfmax){
    x2 = f->GetX(halfmax, max_x, end);
  }
  else x2 = end;

  double fwhm = x2 - x1;
  return fwhm;
}


double GetFWHM(TH1D *h){
  double fwhm = 0.0;
  int f = 0;
  int nb = h->GetNbinsX();
  double halfmax = h->GetMaximum()/2.;
  if (halfmax < 10) return 0.0;
  for (int i = 1; i < nb+1; i++){
    if (f == 0){
      if (h->GetBinContent(i) > halfmax) {
        fwhm = h->GetBinLowEdge(i);
        f++;
      }
    }
    else if (f == 1){
      if (h->GetBinContent(i) < halfmax){
	fwhm = h->GetBinLowEdge(i) - fwhm;
        f++;
      }
    }
    else break;

  }
  return fwhm;
}
void sepdmip_draw(TFile *fin, int sector)
{
  gStyle->SetOptStat(0);

  char sector_string[5];
  if (sector < 10) sprintf(sector_string,"s0%d",sector);
  else sprintf(sector_string,"s%d",sector);
  
  char outdir[100];
  sprintf(outdir, "sepd_mip_plots/cosmics_%s",sector_string);

  TH1D *hist_test_coincAndIsoFullCut[nCh];

  TF1 *f_landau[nCh];
  TF1 *f_langaus[nCh];
  TF1 *f_gaus[nCh];
  TF1 *f_landau_fixed[nCh];

  
  float_t channel, tile, nsector;
  float_t mean_gaus, mean_gaus_err, width_gaus, width_gaus_err;
  float_t mpv_lan, mpv_lan_err, width_lan, width_lan_err;
  float_t mpv_langaus, mpv_langaus_err, width_langaus, width_langaus_err;
  float_t fwhm_func, fwhm_dist;
  float_t mpv_lan_fixed, mpv_lan_fixed_err, width_lan_fixed, width_lan_fixed_err;

  vector<float> vtile;
  vector<float> landau_width_mpv;
  vector<float> langaus_width_mpv;
  vector<float> langaus_FWHM_mean;
  vector<float> gaus_sigma_mean;
  vector<float> gaus_FWHM_sigma;
  vector<float> landau_fixed_width_mpv;

  TCanvas *c[nCh];
  TText *t = new TText(.5,.5,Form("%s Cosmics Testing", sector_string));
  t->SetTextAlign(22);
  t->SetTextColor(kBlack);
  t->SetTextFont(43);
  t->SetTextSize(40);
  t->Draw();
  
  for (int i = 0; i < nCh; i++){
    c[i]= new TCanvas(Form("c_%d",i),"", 800, 800);
    c[i]->SetGrid(1,1);
    c[i]->cd();
    hist_test_coincAndIsoFullCut[i] = (TH1D*) fin->Get(Form("hist_test_coincAndIsoFullCut_ch%d",i*2+1));
    f_landau[i] = (TF1*) fin->Get(Form("f_landau_ch%d",i*2+1));
    f_langaus[i] = (TF1*) fin->Get(Form("f_langaus_ch%d",i*2+1));
    f_gaus[i] = (TF1*) fin->Get(Form("f_gaus_ch%d",i*2+1));
    f_landau_fixed[i] = (TF1*) fin->Get(Form("f_landau_fixed_ch%d",i*2+1));
    
    SetyjPadStyle();
    cout<<"Tile: "<<i<<endl;

    if (!(f_gaus[i])) {
      cout<<" No function for gaus " <<endl;
      continue;
    }
    if (!(f_landau[i])) {
      cout<<" No function for landau"<<endl;
      continue;
    }
    if (!(f_landau_fixed[i])) {
      cout<<" No function for landaufixed"<<endl;
      continue;
    }
    if (!(f_langaus[i])) {
      cout<<" No function for langaus"<<endl;
      continue;
    }

    mean_gaus = f_gaus[i]->GetParameter(1);
    width_gaus = f_gaus[i]->GetParameter(2);
    mpv_lan = f_landau[i]->GetParameter(1);
    
    width_lan = f_landau[i]->GetParameter(2);
    mpv_lan_err = f_landau[i]->GetParError(1);

    width_lan_err = f_landau[i]->GetParError(2);
    mpv_langaus_err = f_langaus[i]->GetParError(1);
    width_langaus_err = f_langaus[i]->GetParError(0);
    mean_gaus_err = f_gaus[i]->GetParError(1);
    width_gaus_err = f_gaus[i]->GetParError(2);

    mpv_langaus = f_langaus[i]->GetParameter(1);
    width_langaus = f_langaus[i]->GetParameter(0);
    width_lan_fixed = f_landau_fixed[i]->GetParameter(2);
    width_lan_fixed_err = f_landau_fixed[i]->GetParError(2);
    mpv_lan_fixed = f_landau_fixed[i]->GetParameter(1);
    mpv_lan_fixed_err = f_langaus[i]->GetParError(1);
    fwhm_func = GetFWHM(f_langaus[i]);
    fwhm_dist = GetFWHM(hist_test_coincAndIsoFullCut[i]);

    hist_test_coincAndIsoFullCut[i]->SetTitle(";ADC Peak - Pedestal; counts");
  
    SetLineAtt(hist_test_coincAndIsoFullCut[i], kBlack, 2, 1);
    hist_test_coincAndIsoFullCut[i]->SetMaximum(f_gaus[i]->GetMaximum()*2.);
    hist_test_coincAndIsoFullCut[i]->SetMinimum(1);
    hist_test_coincAndIsoFullCut[i]->GetFunction(f_gaus[i]->GetName())->SetBit(TF1::kNotDraw);
    hist_test_coincAndIsoFullCut[i]->Draw();
    SetLineAtt(f_langaus[i], kRed, 3, 1);
    SetLineAtt(f_gaus[i], kViolet, 3, 1);
    SetLineAtt(f_landau[i], kBlue+2, 3, 1);
    f_landau[i]->Draw("same");
    f_langaus[i]->Draw("same");
    f_gaus[i]->Draw("same");


    vtile.push_back(i*2+1);
    landau_width_mpv.push_back(width_lan/mpv_lan);
    langaus_width_mpv.push_back(width_langaus/mpv_langaus);
    langaus_FWHM_mean.push_back((fwhm_func/2.355)/mean_gaus);
    gaus_sigma_mean.push_back(width_gaus/mean_gaus);
    gaus_FWHM_sigma.push_back((fwhm_func/2.355)/width_gaus);
    landau_fixed_width_mpv.push_back(width_lan_fixed/mpv_langaus);

    std::vector<std::string> lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};

    lines.push_back(Form("%s tile %d",sector_string, i*2+1));
    
    lines.push_back(Form("(FWHM/2.355)/Mean_{G}: %.2f", (fwhm_func/2.355)/(mean_gaus)));
    lines.push_back(Form("Width/MPV_{L}: %.2f", (width_lan)/(mpv_lan)));

    MakeTextPrint(lines, 0.53, 0.87, 0.04);
    

    TLegend *tl = new TLegend(0.15, 0.6, 0.5, 0.9);

    tl->AddEntry(hist_test_coincAndIsoFullCut[i], "Coincidence w/ Iso");
    tl->AddEntry(f_gaus[i], Form("Gaussian Mean: %.0f #pm %.0f", mean_gaus, mean_gaus_err));
    tl->AddEntry(f_langaus[i], Form("Landau-Gaus FWHM: %.0f", fwhm_func));    
    tl->AddEntry(f_landau[i], Form("Landau MPV: %.0f #pm %.0f", mpv_lan, mpv_lan_err));


    tl->Draw();

    c[i]->SaveAs(Form("%s/sepd_mip_ch%d.pdf",outdir, i*2+1));
    c[i]->SaveAs(Form("%s/sepd_mip_ch%d.png",outdir, i*2+1));
  }
  
 

  const int bt = vtile.size();
  float  a_tile[bt];
  float  a_landau_width_mpv[bt];
  float  a_langaus_width_mpv[bt];
  float  a_langaus_FWHM_mean[bt];
  float  a_gaus_sigma_mean[bt];
  float  a_gaus_FWHM_sigma[bt];  
  float  a_landau_fixed_width_mpv[bt];

  for (int i = 0; i < bt; i++){
    a_tile[i] = vtile[i];
    a_landau_width_mpv[i] = landau_width_mpv[i];
    a_langaus_width_mpv[i] = langaus_width_mpv[i];
    a_langaus_FWHM_mean[i] = langaus_FWHM_mean[i];
    a_gaus_sigma_mean[i] = gaus_sigma_mean[i];
    a_gaus_FWHM_sigma[i] = gaus_FWHM_sigma[i];
    a_landau_fixed_width_mpv[i] = landau_fixed_width_mpv[i];
  }
  TGraph *g_gaus_sigma_mean = new TGraph(bt, a_tile, a_gaus_sigma_mean);
  TGraph *g_gaus_fwhm_sigma =new TGraph(bt, a_tile, a_gaus_FWHM_sigma);
  TGraph *g_landau_width_mpv =new TGraph(bt, a_tile, a_landau_width_mpv);
  TGraph *g_langaus_width_mpv =new TGraph(bt, a_tile, a_langaus_width_mpv);
  TGraph *g_langaus_FWHM_mean =new TGraph(bt, a_tile, a_langaus_FWHM_mean);
  TGraph *g_landau_fixed_width_mpv = new TGraph(bt, a_tile, a_landau_fixed_width_mpv);

  std::vector<std::string> lines;
  
  TCanvas *c6 = new TCanvas("c6", "c6");
  SetyjPadStyle();
  SetMarkerAtt(g_langaus_FWHM_mean, kBlue, 2, 8);
  g_langaus_FWHM_mean->SetTitle(";Tile Num.; (FWHM/2.355) / Mean_{G}");
  g_langaus_FWHM_mean->Draw("AP");
  g_langaus_FWHM_mean->GetHistogram()->SetMinimum(.1);
  g_langaus_FWHM_mean->GetHistogram()->SetMaximum(.4);
  lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};
  lines.push_back(Form("%s - (FWHM/2.355)/Mean_{Gaus}",sector_string));
  MakeTextPrint(lines, 0.6, 0.82, 0.04);
  
  c6->SaveAs(Form("%s/langaus_fwhm_mean_Sector%d.pdf",outdir, sector));
  c6->SaveAs(Form("%s/langaus_fwhm_mean_Sector%d.png",outdir, sector));
  
  return;
}


void sEPD_Cosmics_Draw_modified(int testsec=24, int refsec=5){
  cout << "Start sEPD_Cosmics_Draw_modified.C ... " << endl;
  TFile *fin = new TFile(Form("root/fit_out_Ref%d_Test%d.root", refsec, testsec),"read");
  sepdmip_draw(fin, testsec);
  fin->Close();

  cout << "End sEPD_Cosmics_Draw_modieid.C ... " << endl;
  return;
}
