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
float findLow(TH1D *h, int buff_max = 2){
  int nb = h->GetNbinsX();
  TH1D *hc = (TH1D*) h->Clone("hc");
  int x = 0;
  int i = 1;

  while (i <= nb){
    int j = i;
    hc->SetBinContent(i, 0);
    while (j < i+5 && j <= nb) {
      hc->Fill(hc->GetBinCenter(i), h->GetBinContent(j));
      j++;
    }

    if (i > 2 && hc->GetBinContent(i) >= hc->GetBinContent(i - 1)) x++; 
    else x = 0;
    if (x == buff_max)  return hc->GetBinCenter(i - buff_max);
    i++;
  }
  return 0;
}
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

void sepdmip_run(
		 TFile *fin, TFile *fout, int sector, int minadc, bool debug = true
	     )
{
  
  gStyle->SetOptStat(0);
  
  char sector_string[5];
  if (sector < 10) sprintf(sector_string,"s0%d",sector);
  else sprintf(sector_string,"s%d",sector);
  
  int xmax = 12000;
  
  //Histogram declaration

  int useData[31];

  TH1D *h_sepdmip[31];
  TH1D *h_sepdmip_coincidence[31];
  TH1D *h_sepdmip_coincidence_w_iso[31];
 
  int ntile[31];

  for (int j = 0; j < 31; j++) ntile[j] = j+1;
  
  for (int i = 0; i < 31; i++){
    if (!(fin->Get(Form("h_sepdmip_coincidence_w_iso_%s_%d",sector_string,ntile[i])))){
      useData[i] = 0;
      cout<<"Not using "<<ntile[i]<<endl;
      continue;
    }
    useData[i] = 1;
    h_sepdmip[i] = (TH1D*) fin->Get(Form("h_sepdmip_%s_%d",sector_string,ntile[i]));
    h_sepdmip_coincidence[i] = (TH1D*) fin->Get(Form("h_sepdmip_coincidence_%s_%d",sector_string,ntile[i]));
    h_sepdmip_coincidence_w_iso[i] = (TH1D*) fin->Get(Form("h_sepdmip_coincidence_w_iso_%s_%d",sector_string,ntile[i]));
  }
  
  if (debug) cout<< "Made all of my histograms..."<<endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();
  
  float fwhm_dist, fwhm_func;
  float mpv_lan, width_lan, mpv_lan_err, width_lan_err;
  float mean_gaus, width_gaus, mean_gaus_err, width_gaus_err;
  float mpv_langaus, width_langaus, mpv_langaus_err, width_langaus_err;
  float width_lan_fixed;
  float xv, w, minv, maxv, pp;
  for (int i = 0; i < 31; i++){
    fout->cd();
    if (useData[i] == 0) continue;
    printf("----- Tile %d ----- \n", ntile[i]);
    char lan_name[20];
    sprintf(lan_name, "f_landau_%s_%d", sector_string, ntile[i]);
    char langaus_name[20];
    sprintf(langaus_name, "f_langaus_%s_%d", sector_string, ntile[i]);
    char gaus_name[20];
    sprintf(gaus_name, "f_gaus_%s_%d", sector_string, ntile[i]);
    char lan_fix_name[20];
    sprintf(lan_fix_name, "f_landau_fixed_%s_%d", sector_string, ntile[i]);

    float xlow = findLow(h_sepdmip_coincidence_w_iso[i]);
    if (xlow < 100.0) {
      printf("Cannot find low point ...\n");
      xlow = 100.;
      continue;
    }
    if (xlow > 5000) xlow = 5000;
    printf("Using low point... %4.1f Now fitting... \n", xlow);
    fwhm_dist = GetFWHM(h_sepdmip_coincidence_w_iso[i]);
    TF1 *mylandau = new TF1(lan_name,"landau",xlow,xmax);
    mylandau->SetParLimits(1, xlow, xmax);
    h_sepdmip_coincidence_w_iso[i]->Fit(lan_name,"RH LL");
    mylandau->Write();

    mpv_lan = mylandau->GetParameter(1);
    mpv_lan_err = mylandau->GetParError(1);
    width_lan = mylandau->GetParameter(0);
    width_lan_err= mylandau->GetParError(0);

    // Setting fit range and start values
    double fr[2];
    double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=xlow;
    fr[1]=xmax;

    pllo[0]=10; pllo[1]=xlow; pllo[2]=1.0; pllo[3]=0.0;
    plhi[0]=3000; plhi[1]=9000; plhi[2]=10000000.0; plhi[3]=5000.0;
    sv[0]=width_lan; sv[1]=mpv_lan; sv[2]=50000.0; sv[3]=300.0;

    double chisqr;
    int    ndf;
    TF1 *fitsnr;
    double SNRPeak, SNRFWHM;


    fitsnr = langaufit(h_sepdmip_coincidence_w_iso[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

    fwhm_func = GetFWHM(fitsnr);

    mpv_langaus = fitsnr->GetParameter(1);
    mpv_langaus_err = fitsnr->GetParError(1);
    width_langaus = fitsnr->GetParameter(2);
    width_langaus_err = fitsnr->GetParError(2);

    langaupro(fp,SNRPeak,SNRFWHM);
    fitsnr->SetName(langaus_name);
    fitsnr->Write();

    mylandau->FixParameter(1, mpv_langaus);
    h_sepdmip_coincidence_w_iso[i]->Fit(lan_fix_name,"RH LL");
    mylandau->SetName(lan_fix_name);
    mylandau->Write();


    if (fitsnr->GetMaximumX() < 1000) {
      xv = 1500;
      w = 500;
    }
    else {
      xv = fitsnr->GetMaximumX();
      w = fwhm_func/2.355;
    }
    if (xv - w < xlow) minv = xlow;
    else minv = xv - w;
    if (xv + w > xmax) maxv = xmax;
    else maxv = xv + w;

    TF1 *mygaus = new TF1(gaus_name,"gaus",minv,maxv);
    mygaus->SetParLimits(1, minv, maxv);
    h_sepdmip_coincidence_w_iso[i]->Fit(gaus_name,"RH LL");

    h_sepdmip_coincidence_w_iso[i]->Draw();

    SetLineAtt(mygaus, kBlue, 2, 1);
    SetLineAtt(mylandau, kRed, 2, 1);
    SetLineAtt(fitsnr, kViolet, 2, 1);

    mygaus->Draw("same");
    fitsnr->Draw("same");
    mylandau->Draw("same");
    mygaus->Write();
    h_sepdmip[i]->Write();
    h_sepdmip_coincidence[i]->Write();
    h_sepdmip_coincidence_w_iso[i]->Write();
    mean_gaus = mygaus->GetParameter(1);
    mean_gaus_err = mygaus->GetParError(1);
    width_gaus = mygaus->GetParameter(2);
    width_gaus_err = mygaus->GetParError(2);
    printf("| tile | mean_gaus | width_gaus | mpv_lan | width_lan | mpv_langaus | fwhm_func | fwhm_dist |  \n");
    printf("| %d | %4.1f | %4.1f | %4.1f | %4.1f |%4.1f | %4.1f | %4.1f |\n",ntile[i], mean_gaus, width_gaus, mpv_lan, width_lan, mpv_langaus, fwhm_func, fwhm_dist);

  }
  

  return;
}


void sEPD_Cosmics_Fitting(const string config_file = "cosmics_config.config"){
  cout << "Start sEPD_Cosmics_Fitting.C ... " << endl;
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
  int maxadc = config_p->GetValue("MAXADC", 1000);
  const string file_stack = config_p->GetValue("FILESTACK","");
  const string file_shift = config_p->GetValue("FILESHIFT","");
  bool debug = config_p->GetValue("DODEBUG", false);
  int bad_sipm = config_p->GetValue("BADSIPM", 1);
  const string outdir = config_p->GetValue("OUTDIR",".");
  const string caption = config_p->GetValue("CAPTION","");
  
  //TFile *fin = new TFile(Form("%s/outfile_%d_%d.root", outdir.c_str(), top_sector, bottom_sector), "read");
  TFile *fin = new TFile(Form("%s/outfile_%d_%d%s.root", outdir.c_str(), top_sector, bottom_sector, (caption=="" ? caption.data() : ("_"+caption).data()) ), "read");
  TFile *fout = new TFile(Form("%s/hists_and_fits_%d_%d%s.root", outdir.c_str(), top_sector, bottom_sector, (caption=="" ? caption.data() : ("_"+caption).data())), "recreate");

  sepdmip_run(fin, fout, top_sector, minadc, debug);

  config_p->Write("config", TObject::kOverwrite);
  fout->Close();

  cout << "End sEPD_Cosmics_Fitting.C ... " << endl;
  return;
}
