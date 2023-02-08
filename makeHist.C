#include <iostream>
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TreeSetting.h"

using namespace std;

void makeHist(int testsec=24, int refsec=5)
{
  TFile *f = new TFile(Form("root/outfile_%d_%d.root",testsec,refsec),"read");

  TTree* t = (TTree*) f->Get("tree");

  SetTree settree_;
  settree_.TreeSetting(t);

  const int nHist=16;
  map<TString, TH1D*> href;
  map<TString, TH1D*> htest;
  for(int i=0; i<nHist; i++){
    href[i] = new TH1D(Form("hist_raw_ref_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),100,0,12000);
    htest[i] = new TH1D(Form("hist_raw_test_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),100,0,12000);
//    hrefcoinc[i] = new TH1D(Form("hist_coinc_test_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),100,0,12000);
//    htestcoinciso[i] = new TH1D(Form("hist_coinciso_test_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),100,0,12000);
  }

  for (auto i : ROOT::TSeqI(t->GetEntries())){
    t->GetEntry(i);
  //  cout << "nHist - pmpRef.size() : " << nHist << ", " <<  pmpRef->size() << endl; 
    for(int j=0; j<nHist; j++){
 //     cout << " j : " << j << " pmpRef->at(j) : " << pmpRef->at(j) << " - pmpTest->at(j) :" << pmpTest->at(j) << endl;
      href[j]->Fill(pmpRef->at(j));
      htest[j]->Fill(pmpTest->at(j));
  //    if(isCoincOld.at(j)==true){
  //      hrefcoinc[i]->Fill(pmpTest.at(j));
  //      if(isIsoOld.at(j)==true) htestcoinciso[i]->Fill(pmpTest.at(j));
  //    }
    }
  }

  TCanvas* cref = new TCanvas("cref","",700,700);
  TCanvas* ctest = new TCanvas("ctest","",700,700);
  //TCanvas* ccoinc = new TCanvas("ccoinc","",700,700);
  //TCanvas* ciso = new TCanvas("ciso","",700,700);

  cref->Divide(4,4);
  ctest->Divide(4,4);
  for(int i=0; i<nHist; i++){
    cref->cd(i+1); gPad->SetLogy(); href[i]->Draw();
    ctest->cd(i+1); gPad->SetLogy(); htest[i]->Draw();
  }
  cref->SaveAs(Form("sepd_mip_plots/cosmics_s%d/RawDist_refSec%d.pdf",testsec,refsec));
  ctest->SaveAs(Form("sepd_mip_plots/cosmics_s%d/RawDist_testSec%d.pdf",testsec,testsec));

}

