#include <iostream>
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "Utility/TreeSetting.h"
#include "Utility/Style_jaebeom.h"
#include "Utility/commonUtility.h"

using namespace std;

void makeHist(int testsec=24, int refsec=5)
{
  gStyle->SetOptStat(0);

  string outdir = Form("sepd_mip_plots/cosmics_s%d",testsec);
  void * dirf = gSystem->OpenDirectory(outdir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(outdir.c_str(), kTRUE);}

  TFile *f = new TFile(Form("root/outfile_%d_%d.root",testsec,refsec),"read");

  TTree* t = (TTree*) f->Get("tree");

  SetTree settree_;
  settree_.TreeSetting(t);

  ifstream openFile(Form("Utility/adcpedcuttest_sec%d.txt",testsec));
  string line;
  int ii=0;
  while (getline(openFile, line))
  {
    istringstream iss(line);
    vector<double> vec(4);
    iss >> vec[0] >> vec[1] >> vec[2] >> vec[3]; 
    adcpedcuttest[ii*4+0] = vec[0];
    adcpedcuttest[ii*4+1] = vec[1];
    adcpedcuttest[ii*4+2] = vec[2];
    adcpedcuttest[ii*4+3] = vec[3];
    ii++;
  };


  TFile* wf = new TFile(Form("root/hist_out_Ref%d_Test%d.root",refsec,testsec),"recreate");

  const int nBins=100;
  double xmin=0; double xmax=12000;
  map<TString, TH1D*> href;
  map<TString, TH1D*> htest;
  map<TString, TH1D*> htest_refCut;
  map<TString, TH1D*> htest_coincCut;
  map<TString, TH1D*> htest_coincCut_Upper;
  map<TString, TH1D*> htest_coincCut_Lower;
  map<TString, TH1D*> htest_coincAndIsoRefCut;
  map<TString, TH1D*> htest_coincAndIsoFullCut;
  for(int i=0; i<nCh; i++){
    href[i] = new TH1D(Form("hist_raw_ref_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest[i] = new TH1D(Form("hist_raw_test_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest_refCut[i] = new TH1D(Form("hist_test_refCut_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest_coincCut[i] = new TH1D(Form("hist_test_coincCut_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest_coincCut_Upper[i] = new TH1D(Form("hist_test_coincCut_Upper_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest_coincCut_Lower[i] = new TH1D(Form("hist_test_coincCut_Lower_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest_coincAndIsoRefCut[i] = new TH1D(Form("hist_test_coincAndIsoRefCut_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);
    htest_coincAndIsoFullCut[i] = new TH1D(Form("hist_test_coincAndIsoFullCut_ch%d",i*2+1),Form(";ADC;Counts (Ch:%d)",i*2+1),nBins,xmin,xmax);

    SetHistStyleSmall(href[i], 2, 0);
  }

  for (auto i : ROOT::TSeqI(t->GetEntries())){
    t->GetEntry(i);
    if(i%100000==0) cout << ">>>>> EVENT " << i << " / " << t->GetEntries() <<  " ("<<(int)(100.*i/t->GetEntries()) << "%)" << endl;
    for(int j=0; j<nCh; j++){
      href[j]->Fill(pmpRef->at(j));
      htest[j]->Fill(pmpTest->at(j));
      
      //Pedestal cut on reference sector
      if(!settree_.ApplyPedCutRef(j)) continue;
      htest_refCut[j]->Fill(pmpTest->at(j));
      
      //Pedestal cut on test sector
      if(!settree_.ApplyPedCutTest(j)) continue;
      htest_coincCut[j]->Fill(pmpTest->at(j));
      
      j<nCh-1 ? htest_coincCut_Upper[j]->Fill(pmpTest->at(j+1)) : htest_coincCut_Upper[j]->Fill(-10);
      j>0 ? htest_coincCut_Lower[j]->Fill(pmpTest->at(j-1)) :  htest_coincCut_Lower[j]->Fill(-10);

      //Isolation cut on reference sector
      if(!settree_.IsoRef(j)) continue; 
      htest_coincAndIsoRefCut[j]->Fill(pmpTest->at(j));
      
      //Isolation cut on test sector
      if(!settree_.IsoTest(j)) continue;
      htest_coincAndIsoFullCut[j]->Fill(pmpTest->at(j));


    }
  }

  TCanvas* cref = new TCanvas("cref","",850,850);
  TCanvas* cref_wl = new TCanvas("cref_wl","",850,850);
  TCanvas* ctest = new TCanvas("ctest","",850,850);
  TCanvas* ctest_wl = new TCanvas("ctest_wl","",850,850);
  TCanvas* ctest_refCut = new TCanvas("ctest_refCut","",850,850);
  TCanvas* ctest_CoincCut = new TCanvas("ctest_CoincCut","",850,850);
  TCanvas* ctest_CoincCut_IsoRefCut = new TCanvas("ctest_CoincCut_IsoRefCut","",850,850);
  TCanvas* ctest_CoincCut_IsoFullCut = new TCanvas("ctest_CoincCut_IsoFullCut","",850,850);

  cref->Divide(4,4);
  cref_wl->Divide(4,4);
  ctest->Divide(4,4);
  ctest_refCut->Divide(4,4);
  ctest_CoincCut->Divide(4,4);
  ctest_CoincCut_IsoRefCut->Divide(4,4);
  ctest_CoincCut_IsoFullCut->Divide(4,4);

  TCanvas* ctest_CoincCut_IsoScan[nCh];

  double xpos=0.49; double ypos=0.84;
  double yposdiff= 0.06;
  int textsize=14;
  double xpostitle = 0.12; double ypostitle = 0.9793; int titlesize = 22;
  
  wf->cd();
  for(int i=0; i<nCh; i++)
  {
    cref->cd(i+1); gPad->SetLogy(); href[i]->Draw(); 
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    cref_wl->cd(i+1); gPad->SetLogy(); href[i]->Draw(); solidLine(adcpedcutref[i],0,adcpedcutref[i],href[i]->GetMaximum(),kRed,1);
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    ctest->cd(i+1); gPad->SetLogy(); htest[i]->Draw(); 
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    ctest_wl->cd(i+1); gPad->SetLogy(); htest[i]->Draw(); solidLine(adcpedcuttest[i],0,adcpedcuttest[i],htest[i]->GetMaximum(),kRed,1);
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    ctest_refCut->cd(i+1); gPad->SetLogy(); htest_refCut[i]->Draw();
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    ctest_CoincCut->cd(i+1); gPad->SetLogy(); htest_coincCut[i]->Draw();
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    ctest_CoincCut_IsoRefCut->cd(i+1); gPad->SetLogy(); htest_coincAndIsoRefCut[i]->Draw();
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    ctest_CoincCut_IsoFullCut->cd(i+1); gPad->SetLogy(); htest_coincAndIsoFullCut[i]->Draw();
    drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize);
    
    //Isolation Scan 
    ctest_CoincCut_IsoScan[i] = new TCanvas(Form("ctest_CoincCut_IsoScan_ch#%d",i*2+1),"",850,450);
    ctest_CoincCut_IsoScan[i]->Divide(3,1);
    ctest_CoincCut_IsoScan[i]->cd(2);  htest_coincCut[i]->Draw(); drawText(Form("channel #%d", i*2+1),xpos, ypos, 1, textsize); drawText("Coincidence dist",xpos, ypos-yposdiff, 1, textsize);
    ctest_CoincCut_IsoScan[i]->cd(1); gPad->SetLogy();htest_coincCut_Upper[i]->Draw(); 
    if(i>0){drawText(Form("channel #%d", (i-1)*2+1),xpos, ypos, 1, textsize); drawText("Previous channel",xpos, ypos-yposdiff, 1, textsize); solidLine(adcpedcuttest[i-1],0,adcpedcuttest[i-1],htest_coincCut_Upper[i-1]->GetMaximum(),kRed,1);}
    ctest_CoincCut_IsoScan[i]->cd(3); gPad->SetLogy(); htest_coincCut_Lower[i]->Draw(); 
    if(i<nCh-1){drawText(Form("channel #%d", (i+1)*2+1),xpos, ypos, 1, textsize); drawText("Next channel",xpos, ypos-yposdiff, 1, textsize); solidLine(adcpedcuttest[i+1],0,adcpedcuttest[i+1],htest_coincCut_Lower[i+1]->GetMaximum(),kRed,1);}

    ctest_CoincCut_IsoScan[i]->cd();  drawText(Form("Isolation scan test sector #%d - Ref. #%d & Test pedestal cut",testsec, refsec), xpostitle, ypostitle-yposdiff, kBlack, titlesize);
    ctest_CoincCut_IsoScan[i]->SaveAs(Form("sepd_mip_plots/cosmics_s%d/IsoScan_Ch%d_CoincidencCut.pdf",testsec,i*2+1));
    ctest_CoincCut_IsoScan[i]->Write();
    
    href[i]->Write();
    htest[i]->Write();
    htest_refCut[i]->Write();
    htest_coincCut[i]->Write();
    htest_coincAndIsoRefCut[i]->Write();
    htest_coincAndIsoFullCut[i]->Write();

  }

  cref->cd(); drawText(Form("Reference Sector #%d raw",refsec), xpostitle, ypostitle, kBlack, titlesize);
  cref_wl->cd(); drawText(Form("Reference Sector #%d raw",refsec), xpostitle, ypostitle, kBlack, titlesize);
  ctest->cd(); drawText(Form("Test Sector #%d (Ref. #%d) raw",testsec, refsec), xpostitle, ypostitle, kBlack, titlesize);
  ctest_wl->cd(); drawText(Form("Test Sector #%d (Ref. #%d) raw",testsec, refsec), xpostitle, ypostitle, kBlack, titlesize);
  ctest_refCut->cd(); drawText(Form("Test Sector #%d - Ref. #%d pedestal cut",testsec, refsec), xpostitle, ypostitle, kBlack, titlesize);
  ctest_CoincCut->cd(); drawText(Form("Test Sector #%d - Ref. #%d & Test pedestal cut",testsec, refsec), xpostitle, ypostitle, kBlack, titlesize);
  ctest_CoincCut_IsoRefCut->cd(); drawText(Form("Test Sector #%d - (Ref. #%d & Test pedestal cut) + (Ref. #%d Iso cut)",testsec, refsec,refsec), xpostitle, ypostitle, kBlack, titlesize);
  ctest_CoincCut_IsoFullCut->cd(); drawText(Form("Test Sector #%d - (Ref. #%d & Test pedestal cut) + (Ref. #%d & Test Iso cut)",testsec, refsec,refsec), xpostitle, ypostitle, kBlack, titlesize);

  cref->SaveAs(Form("sepd_mip_plots/cosmics_s%d/RawDist_refSec%d.pdf",testsec,refsec));
  cref_wl->SaveAs(Form("sepd_mip_plots/cosmics_s%d/RawDist_refSec%d_wline.pdf",testsec,refsec));
  ctest->SaveAs(Form("sepd_mip_plots/cosmics_s%d/RawDist_testSec%d.pdf",testsec,testsec));
  ctest_wl->SaveAs(Form("sepd_mip_plots/cosmics_s%d/RawDist_testSec%d_wline.pdf",testsec,testsec));
  ctest_refCut->SaveAs(Form("sepd_mip_plots/cosmics_s%d/TestDist_refPedestalCut.pdf",testsec));
  ctest_CoincCut->SaveAs(Form("sepd_mip_plots/cosmics_s%d/TestDist_CoincidenceCut.pdf",testsec));
  ctest_CoincCut_IsoRefCut->SaveAs(Form("sepd_mip_plots/cosmics_s%d/TestDist_CoincidenceCut_IsoRefCut.pdf",testsec));
  ctest_CoincCut_IsoFullCut->SaveAs(Form("sepd_mip_plots/cosmics_s%d/TestDist_CoincidenceCut_IsoFullCut.pdf",testsec));

  cref->Write();
  cref_wl->Write();
  ctest->Write();
  ctest_wl->Write();
  ctest_refCut->Write();
  ctest_CoincCut->Write();
  ctest_CoincCut_IsoRefCut->Write();
  ctest_CoincCut_IsoFullCut->Write();

  wf->Close();
}
