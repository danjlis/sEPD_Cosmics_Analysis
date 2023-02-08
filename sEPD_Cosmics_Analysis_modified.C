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

void sepdmip_run(TFile *fout, const string fname, int top_sector, int bottom_sector, int top_oe, int bottom_oe, int minadc = 0, int maxadc = 0, int IsoCutADCMaxSwitch=100, bool shift = false, int bad_sipm = 1 ,bool debug = true)
{

  gStyle->SetOptStat(0);


  const int nADC=384;
  const int nADC_=31;


  char top_sector_string[5];
  char bottom_sector_string[5];

  if (top_sector < 10) sprintf(top_sector_string,"s0%d",top_sector);
  else sprintf(top_sector_string,"s%d",top_sector);
  if (bottom_sector < 10) sprintf(bottom_sector_string,"s0%d",bottom_sector);
  else sprintf(bottom_sector_string,"s%d",bottom_sector);

  int xmax = 12000;
  
  vector<double> pmpRef;
  vector<double> pmpTest;
  vector<bool> isCoincOld;
  vector<bool> isIsoOld;

  TTree* wt = new TTree("tree","");
  wt->Branch("pmpRef",&pmpRef);
  wt->Branch("pmpTest",&pmpTest);
  wt->Branch("isCoincOld",&isCoincOld);
  wt->Branch("isIsoOld",&isIsoOld);

  TFile *fin = new TFile(fname.c_str(), "r");

  if (!fin){
    cout<<"Could not find file "<<fname<<"..."<<endl;
    return;
  }

  TTree *w = (TTree*) fin->Get("W");

  if (!w){
    cout<<"TTree W not found in file..."<<endl;
  }

  int adc[nADC][nADC_];

  float peak[nADC];
  float pedestal[nADC];

  w->SetBranchAddress("adc", adc);
  w->SetBranchAddress("peak", peak);
  w->SetBranchAddress("pedestal", pedestal);

  // First and second sipm index start
  const int first = 0;
  const int second = 16;
  const int nHist = 32;
  const int nComp=2;

  // nentries
  int nentries = w->GetEntries();
  cout << "nentries : " << nentries << endl;

  //Histogram declaration

  TH1D *h_sepdmip[nHist];
  TH1D *h_sepdmip_coincidence[nHist];
  TH1D *h_sepdmip_coincidence_w_iso[nHist];


  int nHistBins=100;
  int bin_low_adc = 0.;
  int bin_high_adc = 12000.;
  float tmp_max[nComp];
  float pmp[nComp][second];
  int tile[nHist];
  for (int i = 0; i < nComp; i++){
    tile[i*second] = 1;
    for (int j = 1; j < second; j++) tile[i*second+j] = (i?bottom_oe:top_oe) + j*nComp;
  }

  for (int i = 0; i < nHist; i++){

    h_sepdmip[i] = new TH1D(Form("h_sepdmip_%s_%d",(i < second ? bottom_sector_string: top_sector_string),tile[i]),"", nHistBins, bin_low_adc, bin_high_adc);
    h_sepdmip_coincidence[i] = new TH1D(Form("h_sepdmip_coincidence_%s_%d",(i < second ? bottom_sector_string: top_sector_string), tile[i]),"", nHistBins, bin_low_adc, bin_high_adc);
    h_sepdmip_coincidence_w_iso[i] = new TH1D(Form("h_sepdmip_coincidence_w_iso_%s_%d",(i < second ? bottom_sector_string: top_sector_string), tile[i]),"", nHistBins, bin_low_adc, bin_high_adc);
  }

  if (debug) cout<< "Made all of my histograms..."<<endl;


  int count = 0;
  int iup = 0;



  //for (int e = 0; e < nentries; e++){
  for (auto e : ROOT::TSeqI(nentries)){
    count = 0;
    if (e%10000 == 0) std::cout<<"Progress: "<< 100.*static_cast<float>(e)/static_cast<float>(nentries)<<" % ..."<<std::endl;

    pmpRef.clear();
    pmpTest.clear();
    isCoincOld.clear();
    isIsoOld.clear();

    w->GetEntry(e);

    for (int i = 0 ; i < second; i++){
      for (int j = 0; j < nComp; j++){
        pmp[j][i] = peak[i+(j<<4)] - pedestal[i+(j<<4)];
        j==0 ? pmpRef.push_back(pmp[j][i]) : pmpTest.push_back(pmp[j][i]);
        if(pmp[j][i] > minadc) h_sepdmip[i + (j<<4)]->Fill(pmp[j][i]);
      }
    }
    if (shift) {
      if (debug) {
        cout<<"Event "<<e<<" pmp: "<<endl;
        for (int i = 0; i < 2; i++){
          if(i == 0) cout<<top_sector_string<<" |";
          else cout<<bottom_sector_string<<" |";

          for (int j = 0; j < 16; j++){
            printf("| %5.0f |",pmp[1-i][j]);
          }
          cout<<"|"<<endl;
        }
      }

      // Coincidence cut
      if (pmp[0][bad_sipm+1] > minadc && pmp[1][bad_sipm] > minadc) {
        h_sepdmip_coincidence[bad_sipm + 1 + (0<<4)]->Fill(pmp[0][bad_sipm+1]);
        h_sepdmip_coincidence[bad_sipm + (1<<4)]->Fill(pmp[1][bad_sipm]);

        tmp_max[0]= (maxadc > 100? maxadc : (static_cast<float>(maxadc)/100.)*pmp[0][bad_sipm+1]);
        tmp_max[1] = (maxadc > 100? maxadc : (static_cast<float>(maxadc)/100.)*pmp[1][bad_sipm]);

        // isolation cut  
        if (pmp[0][bad_sipm+2] < tmp_max[0] && pmp[1][bad_sipm+1] < tmp_max[1] && pmp[0][bad_sipm] < tmp_max[0] && pmp[1][bad_sipm-1] < tmp_max[1]){
          h_sepdmip_coincidence_w_iso[bad_sipm+1 + (0<<4)]->Fill(pmp[0][bad_sipm+1]);
          h_sepdmip_coincidence_w_iso[bad_sipm + (1<<4)]->Fill(pmp[1][bad_sipm]);
        }
      }
    }
    else{
      if (debug) {
        cout<<"Event "<<e<<" pmp: "<<endl;
        for (int i = 0; i < 2; i++){
          if(i == 0) cout<<top_sector_string<<" |";
          else cout<<bottom_sector_string<<" |";

          for (int j = 0; j < 16; j++){
            printf("| %5.0f |",pmp[1-i][j]); 
          }
          cout<<"|"<<endl;
        }
      }

      // Coincidence cut
      for (int i = 0 ; i < 16; i++){
        if (i == 0 && top_oe == 0) continue;
        isCoincOld.push_back(false);
        isIsoOld.push_back(false);

        if (i == bad_sipm){
          if (pmp[0][i+1] > minadc && pmp[1][i] > minadc) {
            isCoincOld.push_back(true);
            h_sepdmip_coincidence[i + 1 + (0<<4)]->Fill(pmp[0][i+1]);
            h_sepdmip_coincidence[i + (1<<4)]->Fill(pmp[1][i]);

            tmp_max[0] = (maxadc > IsoCutADCMaxSwitch? maxadc : (static_cast<float>(maxadc)/100.)*pmp[0][i+1]);
            tmp_max[1] = (maxadc > IsoCutADCMaxSwitch? maxadc : (static_cast<float>(maxadc)/100.)*pmp[1][i]);

            // isolation cut                                                                               
            if ( (pmp[0][i+2] < tmp_max[0] &&  pmp[0][i] < tmp_max[0]) && (pmp[1][i+1] < tmp_max[1] && pmp[1][i-1] < tmp_max[1]) ){
              h_sepdmip_coincidence_w_iso[i + (0<<4) + 1] -> Fill(pmp[0][i+1]);
              h_sepdmip_coincidence_w_iso[i + (1<<4)]     -> Fill(pmp[1][i]);
              isIsoOld.push_back(true);
            }
          }
        }
        else if (pmp[0][i] > minadc && pmp[1][i] > minadc) {
          isCoincOld.push_back(true);
          h_sepdmip_coincidence[i + (0<<4)]->Fill(pmp[0][i]);
          h_sepdmip_coincidence[i + (1<<4)]->Fill(pmp[1][i]);
          tmp_max[0] = (maxadc > 100? maxadc : (static_cast<float>(maxadc)/100.)*pmp[0][i]);
          tmp_max[1] = (maxadc > 100? maxadc : (static_cast<float>(maxadc)/100.)*pmp[1][i]);

          if (top_oe == 1){
            if ((i != 0? pmp[0][i-1] < tmp_max[0] && pmp[1][i-1] < tmp_max[1]: true) && (i != 15? pmp[0][i+1] < tmp_max[0] && pmp[1][i+1] < tmp_max[1]: true)){
              isIsoOld.push_back(true);
              h_sepdmip_coincidence_w_iso[i + (0<<4)]->Fill(pmp[0][i]);
              h_sepdmip_coincidence_w_iso[i + (1<<4)]->Fill(pmp[1][i]);
            }
          }
          else {
            if ((i != 1? pmp[0][i-1] < tmp_max[0] && pmp[1][i-1] < tmp_max[1]: true) && (i != 15? pmp[0][i+1] < tmp_max[0] && pmp[1][i+1] < tmp_max[1] : true)){
              isIsoOld.push_back(true);
              h_sepdmip_coincidence_w_iso[i + (0<<4)]->Fill(pmp[0][i]);
              h_sepdmip_coincidence_w_iso[i + (1<<4)]->Fill(pmp[1][i]);
            }
          }
        }
      }
      /*
      // isolation cut
      for (int i = 0; i < 16; i++){
      if(i==bad_sipm) continue;
      tmp_max[0]= (maxadc > 100? maxadc : (static_cast<float>(maxadc)/100.)*pmp[0][i]);
      tmp_max[1] = (maxadc > 100? maxadc : (static_cast<float>(maxadc)/100.)*pmp[1][i]);

      if (pmp[0][i] < minadc && pmp[1][i] < minadc){
      h_sepdmip_coincidence_w_iso[i + (0<<4)]->Fill(pmp[0][i]);
      h_sepdmip_coincidence_w_iso[i + (1<<4)]->Fill(pmp[1][i]);
      }


      if(i==0){
      if (pmp[0][0] > minadc && pmp[1][0] > minadc && pmp[0][1] < minadc&&pmp[1][1] < minadc){
      h_sepdmip_coincidence_w_iso[i + (0<<4)]->Fill(pmp[0][i]);
      h_sepdmip_coincidence_w_iso[i + (1<<4)]->Fill(pmp[1][i]);
      }
      }
      else if (i==15){
      if (pmp[0][15] > minadc && pmp[1][15] > minadc && pmp[0][14] < minadc&&pmp[1][14] < minadc){
      h_sepdmip_coincidence_w_iso[i + (0<<4)]->Fill(pmp[0][i]);
      h_sepdmip_coincidence_w_iso[i + (1<<4)]->Fill(pmp[1][i]);
      }
      }
      else {
      if (pmp[0][i] > minadc && pmp[1][i] > minadc && pmp[0][i+1] < minadc&&pmp[1][i+1] < minadc&&pmp[0][i-1] < minadc&&pmp[1][i-1] < minadc){
      h_sepdmip_coincidence_w_iso[i + (0<<4)]->Fill(pmp[0][i]);
      h_sepdmip_coincidence_w_iso[i + (1<<4)]->Fill(pmp[1][i]);
      }
      }
      }
      */
    }
    wt->Fill();
  }

  for (int i = 0; i < 32; i++){
    // Check if the file goes through a shifted test for a shifted coincidence
    //if(!shift && (i==bad_sipm||i==(bad_sipm + 16))) continue;
    if(shift && (i!=bad_sipm + 16)) continue;
    //if(!top_oe && (i%16 == 0)) continue;

    fout->cd();
    h_sepdmip[i]->Write();
    h_sepdmip_coincidence[i]->Write();
    h_sepdmip_coincidence_w_iso[i]->Write();
  }
  wt->Write();

  return;
}


void sEPD_Cosmics_Analysis_modified(const string config_file = "cosmics_config.config"){
  cout << "Start sEPD_Cosmics_Analysis.C ... " << endl;
  TEnv *config_p = new TEnv(config_file.c_str());
  if (!config_p) {
    cout<<"No configuration file..."<<endl;
    return;
  }
  int top_sector = config_p->GetValue("TOPSECTOR", 8 );
  int bottom_sector = config_p->GetValue("BOTTOMSECTOR", 5 );
  int oe = config_p->GetValue("ODDOREVEN", 1 );
  int even_run = config_p->GetValue("EVENRUN", 0);
  int odd_run = config_p->GetValue("ODDRUN", 0);
  int minadc = config_p->GetValue("MINADC", 1000);
  int maxadc = config_p->GetValue("MAXADC", 1000);
  int IsoCutADCMaxSwitch = 100;
  
  const string file_stack_odd = config_p->GetValue("FILESTACKODD","");
  const string file_shift_odd = config_p->GetValue("FILESHIFTODD","");
  const string file_stack_even = config_p->GetValue("FILESTACKEVEN","");
  const string file_shift_even = config_p->GetValue("FILESHIFTEVEN","");
 
  bool debug = config_p->GetValue("DODEBUG", false);
  int bad_sipm = config_p->GetValue("BADSIPM", 1);
  const string outdir = config_p->GetValue("OUTDIR",".");
  const string caption = config_p->GetValue("CAPTION","");

  void * dirf = gSystem->OpenDirectory(outdir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(outdir.c_str(), kTRUE);}

  TFile *fout = new TFile(Form("%s/outfile_%d_%d%s.root", outdir.c_str(), top_sector, bottom_sector, (caption=="" ? caption.data() : ("_"+caption).data()) ), "recreate");

  if (odd_run){
    if (!file_stack_odd.empty()) sepdmip_run(fout, file_stack_odd, top_sector, bottom_sector, 1, 1, minadc, maxadc, IsoCutADCMaxSwitch, false, bad_sipm, debug);
    if (!file_shift_odd.empty()) sepdmip_run(fout, file_shift_odd, top_sector, bottom_sector, 1, 1, minadc, maxadc, IsoCutADCMaxSwitch, true, bad_sipm, debug);
  }

  if (even_run){
    if (!file_stack_even.empty()) sepdmip_run(fout, file_stack_even, top_sector, bottom_sector, 0, 0, minadc, maxadc, IsoCutADCMaxSwitch, false, bad_sipm, debug);
    if (!file_stack_odd.empty()) sepdmip_run(fout, file_shift_even, top_sector, bottom_sector, 0, 0, minadc,  maxadc, IsoCutADCMaxSwitch, true, bad_sipm, debug);
  }

  config_p->Write("config", TObject::kOverwrite);
  fout->Close();

  cout << "End sEPD_Cosmics_Analysis.C ... " << endl;
  return;
}
