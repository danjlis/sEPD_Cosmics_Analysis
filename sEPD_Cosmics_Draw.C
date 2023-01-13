
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include "dlUtility.h"
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
void sepdmip_draw(
		  TFile *fin, int sector, int sector_oe, int minadc
	     )
{
  gStyle->SetOptStat(0);

  char sector_string[5];
  if (sector < 10) sprintf(sector_string,"s0%d",sector);
  else sprintf(sector_string,"s%d",sector);
  
  char outdir[100];
  sprintf(outdir, "sepd_mip_plots/cosmics_%s/",sector_string);

  TH1D *h_sepdmip[31];
  TH1D *h_sepdmip_coincidence[31];
  TH1D *h_sepdmip_coincidence_w_iso[31];
  TF1 *f_landau[31];
  TF1 *f_langaus[31];
  TF1 *f_gaus[31];
  TF1 *f_landau_fixed[31];
  int ntile[31];
  int useData[31];
  
  for (int j = 0; j < 31; j++) ntile[j] = j+1;
  
  
  for (int i = 0; i < 31; i++){
    if (!(fin->Get(Form("h_sepdmip_coincidence_w_iso_%s_%d",sector_string,ntile[i])))){
      useData[i] = 0;
      continue;
    }
    if (!(fin->Get(Form("f_gaus_%s_%d",sector_string,ntile[i])))){
      useData[i] = 0;
      continue;
    }
    useData[i] = 1;
    h_sepdmip[i] = (TH1D*) fin->Get(Form("h_sepdmip_%s_%d",sector_string,ntile[i]));
    h_sepdmip_coincidence[i] = (TH1D*) fin->Get(Form("h_sepdmip_coincidence_%s_%d",sector_string,ntile[i]));
    h_sepdmip_coincidence_w_iso[i] = (TH1D*) fin->Get(Form("h_sepdmip_coincidence_w_iso_%s_%d",sector_string,ntile[i]));
    f_landau[i] = (TF1*) fin->Get(Form("f_landau_%s_%d",sector_string,ntile[i]));
    f_langaus[i] = (TF1*) fin->Get(Form("f_langaus_%s_%d",sector_string,ntile[i]));
    f_gaus[i] = (TF1*) fin->Get(Form("f_gaus_%s_%d",sector_string,ntile[i]));
    f_landau_fixed[i] = (TF1*) fin->Get(Form("f_landau_fixed_%s_%d",sector_string,ntile[i]));

  }

  
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

  TCanvas *c1 = new TCanvas("c1","c1", 800, 800);
  TText *t = new TText(.5,.5,Form("%s Cosmics Testing", sector_string));
  t->SetTextAlign(22);
  t->SetTextColor(kBlack);
  t->SetTextFont(43);
  t->SetTextSize(40);
  t->Draw();
  c1->Print(Form("cosmics_%s.pdf(",sector_string),"Title");

  SetyjPadStyle();
  


  
  for (int i = 0; i < 32; i++){
    SetyjPadStyle();
    c1->SetGrid(1,1);
    cout<<"Tile: "<<i<<endl;
    if (useData[i] == 0) continue;

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
    fwhm_dist = GetFWHM(h_sepdmip_coincidence_w_iso[i]);

    h_sepdmip_coincidence_w_iso[i]->SetTitle(";ADC Peak - Pedestal; counts");
  
    SetLineAtt(h_sepdmip_coincidence_w_iso[i], kBlack, 2, 1);
    h_sepdmip_coincidence_w_iso[i]->SetMaximum(f_gaus[i]->GetMaximum()*2.);
    h_sepdmip_coincidence_w_iso[i]->SetMinimum(1);
    h_sepdmip_coincidence_w_iso[i]->GetFunction(f_gaus[i]->GetName())->SetBit(TF1::kNotDraw);
    h_sepdmip_coincidence_w_iso[i]->Draw();
    SetLineAtt(f_langaus[i], kRed, 3, 1);
    SetLineAtt(f_gaus[i], kViolet, 3, 1);
    SetLineAtt(f_landau[i], kBlue+2, 3, 1);
    f_landau[i]->Draw("same");
    f_langaus[i]->Draw("same");
    f_gaus[i]->Draw("same");


    vtile.push_back(ntile[i]);
    landau_width_mpv.push_back(width_lan/mpv_lan);
    langaus_width_mpv.push_back(width_langaus/mpv_langaus);
    langaus_FWHM_mean.push_back((fwhm_func/2.355)/mean_gaus);
    gaus_sigma_mean.push_back(width_gaus/mean_gaus);
    gaus_FWHM_sigma.push_back((fwhm_func/2.355)/width_gaus);
    landau_fixed_width_mpv.push_back(width_lan_fixed/mpv_langaus);

    std::vector<std::string> lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};

    lines.push_back(Form("%s tile %d",sector_string, ntile[i]));
    
    lines.push_back(Form("(FWHM/2.355)/Mean_{G}: %.2f", (fwhm_func/2.355)/(mean_gaus)));
    lines.push_back(Form("Width/MPV_{L}: %.2f", (width_lan)/(mpv_lan)));

    MakeTextPrint(lines, 0.53, 0.87, 0.04);
    

    TLegend *tl = new TLegend(0.15, 0.6, 0.5, 0.9);

    tl->AddEntry(h_sepdmip_coincidence_w_iso[i], "Coincidence w/ Iso");
    tl->AddEntry(f_gaus[i], Form("Gaussian Mean: %.0f #pm %.0f", mean_gaus, mean_gaus_err));
    tl->AddEntry(f_langaus[i], Form("Landau-Gaus FWHM: %.0f", fwhm_func));    
    tl->AddEntry(f_landau[i], Form("Landau MPV: %.0f #pm %.0f", mpv_lan, mpv_lan_err));


    tl->Draw();

    c1->Print(Form("cosmics_%s.pdf",sector_string),Form("Tile %d", ntile[i]));
    c1->SaveAs(Form("%s/sepd_mip_%s_%d.pdf",outdir, sector_string, ntile[i]));
    c1->SaveAs(Form("%s/sepd_mip_%s_%d.png",outdir, sector_string, ntile[i]));
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
  /*
  TCanvas *c2 = new TCanvas("c2", "c2");
  SetyjPadStyle();
  SetMarkerAtt(g_gaus_sigma_mean, kBlue, 2, 8);
  g_gaus_sigma_mean->SetTitle(";Tile Num.; #sigma / Mean");
  g_gaus_sigma_mean->Draw("AP");
 
  lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};
  lines.push_back(Form("%s - #sigma/Mean",sector_string));
  MakeTextPrint(lines, 0.6, 0.82, 0.04);
  c2->Print(Form("cosmics_%s.pdf",sector_string),"sigma over mean");

  TCanvas *c3 = new TCanvas("c3", "c3");
  SetyjPadStyle();
  SetMarkerAtt(g_gaus_fwhm_sigma, kBlue, 2, 8);
  g_gaus_fwhm_sigma->SetTitle(";Tile Num.; FWHM / #sigma");
  g_gaus_fwhm_sigma->Draw("AP");
  lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};
  lines.push_back(Form("%s - FWHM/#sigma",sector_string));
  MakeTextPrint(lines, 0.6, 0.82, 0.04);
  c3->Print(Form("cosmics_%s.pdf",sector_string),"FWHM over sigma");

  TCanvas *c4 = new TCanvas("c4", "c4");
  SetyjPadStyle();
  SetMarkerAtt(g_landau_width_mpv, kBlue, 2, 8);
  g_landau_width_mpv->SetTitle(";Tile Num.; Width_{landau} / MPV_{landau}");
  g_landau_width_mpv->Draw("AP");
  lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};
  lines.push_back(Form("%s - Width_{landau} / MPV_{landau}",sector_string));
  MakeTextPrint(lines, 0.6, 0.82, 0.04);
  c4->Print(Form("cosmics_%s.pdf",sector_string),"width over mpv landau");

  TCanvas *c5 = new TCanvas("c5", "c5");
  SetyjPadStyle();
  SetMarkerAtt(g_langaus_width_mpv, kBlue, 2, 8);
  g_langaus_width_mpv->SetTitle(";Tile Num.; Width_{langaus} / MPV_{langaus}");
  g_langaus_width_mpv->Draw("AP");
  lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};
  lines.push_back(Form("%s - Width_{langaus} / MPV_{langaus}",sector_string));
  MakeTextPrint(lines, 0.6, 0.82, 0.04);
  c5->Print(Form("cosmics_%s.pdf",sector_string),"width over mpv langaus");
  */
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
  /*
  TCanvas *c7 = new TCanvas("c7", "c7");
  SetyjPadStyle();
  SetMarkerAtt(g_landau_fixed_width_mpv, kBlue, 2, 8);
  g_landau_fixed_width_mpv->SetTitle(";Tile Num.; Width_{Landau, fixed} / MPV_{langaus}");
  g_landau_fixed_width_mpv->Draw("AP");
  lines = {"#bf{sPHENIX} sEPD Cosmics","Two sector stack"};
  lines.push_back(Form("%s - Width_{L, fixed}/MPV_{LG}",sector_string));
  MakeTextPrint(lines, 0.6, 0.82, 0.04);
  */

  /*
  c2->SaveAs(Form("%s/gaus_sigma_mean_%s.pdf",outdir, sector_string));
  c2->SaveAs(Form("%s/gaus_sigma_mean_%s.png",outdir, sector_string));
  c3->SaveAs(Form("%s/gaus_fwhm_sigma_%s.pdf",outdir, sector_string));
  c3->SaveAs(Form("%s/gaus_fwhm_sigma_%s.png",outdir, sector_string));
  c4->SaveAs(Form("%s/landau_width_mpv_%s.pdf",outdir, sector_string));
  c4->SaveAs(Form("%s/landau_width_mpv_%s.png",outdir, sector_string));
  c5->SaveAs(Form("%s/langaus_width_mpv_%s.pdf",outdir, sector_string));
  c5->SaveAs(Form("%s/langaus_width_mpv_%s.png",outdir, sector_string));
  */ 
  c6->SaveAs(Form("%s/langaus_fwhm_mean_%s.pdf",outdir, sector_string));
  c6->SaveAs(Form("%s/langaus_fwhm_mean_%s.png",outdir, sector_string));
  // c7->SaveAs(Form("%s/landau_fixed_width_mpv_%s.pdf",outdir, sector_string));
  //c7->SaveAs(Form("%s/landau_fixed_width_mpv_%s.png",outdir, sector_string));

  c6->Print(Form("cosmics_%s.pdf)",sector_string),"fwhm over mpv");
  //c7->Print(Form("cosmics_%s.pdf",sector_string),"fixed width over mpv");
  return;
}


void sEPD_Cosmics_Draw(const string config_file = "cosmics_config.config"){
  TEnv *config_p = new TEnv(config_file.c_str());
  if (!config_p) {
    cout<<"No configuration file..."<<endl;
    return;
  }
  int top_sector = config_p->GetValue("TOPSECTOR", 8 );
  int bottom_sector = config_p->GetValue("BOTTOMSECTOR", 5 );
  int top_oe = config_p->GetValue("TOPODDOREVEN", 1 );
  int bottom_oe = config_p->GetValue("BOTTOMODDOREVEN", 1 );

  int minadc = config_p->GetValue("MINADC", 1000);
  const string file_stack = config_p->GetValue("FILESTACK","");
  const string file_shift = config_p->GetValue("FILESHIFT","");
  bool debug = config_p->GetValue("DODEBUG", false);
  int bad_sipm = config_p->GetValue("BADSIPM", 1);
  const string outdir = config_p->GetValue("OUTDIR",".");

  TFile *fin = new TFile(Form("%s/hists_and_fits_%d_%d.root", outdir.c_str(), top_sector, bottom_sector), "r");
  sepdmip_draw(fin, top_sector, top_oe, minadc);
  fin->Close();

  return;
}
