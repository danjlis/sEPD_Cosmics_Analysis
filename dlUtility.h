#ifndef UTILITY_Daniel_H
#define UTILITY_Daniel_H

//TREE,HIST,GRAPH,VECTOR ...
    // std::clock()
//etc.
//#include "tdrstyle.C"   // std::clock()
using namespace std;
void SetLineAtt(TLine *h, Color_t color, int width, int style){
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
}

void SetLineAtt(TF1 *h, Color_t color, int width, int style){
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
}


void SetLineAtt(TGraph *h, Color_t color, int width, int style){
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
}

void SetMarkerAtt(TGraph *h, Color_t color, int size, int style){
  h->SetMarkerColor(color);
  h->SetMarkerSize(size);
  h->SetMarkerStyle(style);
}

void SetLineAtt(TH1D *h, Color_t color, int width, int style){
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
}
void SetLineAtt(TProfile *h, Color_t color, int width, int style){
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
}

void SetMarkerAtt(TProfile *h, Color_t color, int size, int style){
  h->SetMarkerColor(color);
  h->SetMarkerSize(size);
  h->SetMarkerStyle(style);
}

void SetMarkerAtt(TH1D *h, Color_t color, int size, int style){
  h->SetMarkerColor(color);
  h->SetMarkerSize(size);
  h->SetMarkerStyle(style);
}

void SetLineAtt(TEfficiency *h, Color_t color, int width, int style){
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
}

void SetMarkerAtt(TEfficiency *h, Color_t color, int size, int style){
  h->SetMarkerColor(color);
  h->SetMarkerSize(size);
  h->SetMarkerStyle(style);
}


void SetyjPadStyle(){
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

//void SetHistTitleStyle(double titlesize=0.06, double titleoffset=0.04, double labelsize = 0.05, double labeloffset=0.01){
//    gStyle->SetTitleSize( titlesize, "X" ); gStyle->SetTitleOffset(titleoffset, "X");
//   gStyle->SetTitleSize( titlesize, "Y" ); gStyle->SetTitleOffset(titleoffset, "Y");
//  gStyle->SetLabelSize( labelsize, "X" ); gStyle->SetLabelOffset(labeloffset, "X");
// gStyle->SetLabelSize( labelsize, "Y" ); gStyle->SetLabelOffset(labeloffset, "Y");
//}

void thisPadStyle(){
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.08);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
}
void SetPadStyle(){
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
}

void drawText(const char *text, float xp, float yp, bool isRightAlign=0, int textColor=kBlack, double textSize=0.04, int textFont = 42, bool isNDC=true){
  // when textfont 42, textSize=0.04
  // when textfont 43, textSize=18
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  //   if(bold)tex->SetTextFont(43);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  if(isNDC) tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
}
void jumSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
  TLine* t1 = new TLine(x1,y1,x2,y2);
  t1->SetLineWidth(width);
  t1->SetLineStyle(7);
  t1->SetLineColor(color);
  t1->Draw();
}

void onSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
  TLine* t1 = new TLine(x1,y1,x2,y2);
  t1->SetLineWidth(width);
  t1->SetLineStyle(1);
  t1->SetLineColor(color);
  t1->Draw();
}
double findCross(TH1* h1, TH1* h2, double& frac, double& effi, double& fracErr, double& effiErr){
  Int_t nBins = h1->GetNbinsX();
  double crossVal =0;
  int binAt0 = h1->FindBin(0);
  for(Int_t ix=binAt0; ix<=nBins ;ix++){
    float yy1 = h1->GetBinContent(ix);
    float yy2 = h2->GetBinContent(ix);
    if(yy2>yy1) {
      crossVal= h1->GetBinLowEdge(ix);
      break;
    }
  }
  int crossBin = h1->FindBin(crossVal);
  frac = 1 - (h2->Integral(1,crossBin) / h1->Integral(1,crossBin) );
  effi = ( h1->Integral(1,crossBin) / h1->Integral() );
  fracErr = frac * TMath::Sqrt( (1./h2->Integral(1,crossVal)) + (1./h1->Integral(1,crossVal)) );
  effiErr = ( TMath::Sqrt(h1->Integral(1,crossVal)) / h1->Integral() ) * TMath::Sqrt(1 - (h1->Integral(1,crossVal)/h1->Integral()) );

  return crossVal;
}

void ratioPanelCanvas(TCanvas*& canv,
    const Float_t divRatio=0.4,
    //const Float_t leftOffset=0.,
    //const Float_t bottomOffset=0.,
    const Float_t leftMargin=0.17,
    const Float_t bottomMargin=0.3,
    const Float_t edge=0.05) {
  if (canv==0) {
    return;
  }
  canv->Clear();


  TPad* pad1 = new TPad("pad1","",0.0,divRatio,1.0,1.0);
  canv->cd();
  pad1->SetLeftMargin(leftMargin);
  pad1->SetRightMargin(edge);
  pad1->SetTopMargin(edge*2);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  pad1->SetNumber(1);

  TPad* pad2 = new TPad("pad2","",0.0,0.0,1.0,divRatio);
  canv->cd();
  pad2->SetLeftMargin(leftMargin);
  pad2->SetRightMargin(edge);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(bottomMargin);
  pad2->Draw();
  pad2->cd();
  pad2->SetNumber(2);
}

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
    const Int_t rows,
    const Float_t leftOffset=0.,
    const Float_t bottomOffset=0.,
    const Float_t leftMargin=0.2,
    const Float_t bottomMargin=0.2,
    const Float_t edge=0.05,
    const Float_t edge_bottom=0.,
    const Float_t edge_top=0.) {
  if (canv==0) {
    //Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
        (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
        (1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  // Xlow[0] = 0;
  Xup[0] = leftOffset + PadWidth/(1.0-leftOffset);
  //Xup[0] = leftOffset + PadWidth/(1.0-leftOffset);
  //Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
          Xlow[i],Ylow[j],Xup[i],Yup[j]);
      // if(i==0) pad[i][j]->SetLeftMargin(-leftMargin+edge);
      if(i==0) pad[i][j]->SetLeftMargin(edge);
      else pad[i][j]->SetLeftMargin(edge);
      //else pad[i][j]->SetLeftMargin(0);
      //else pad[i][j]->SetLeftMargin(0.03);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);
      //else pad[i][j]->SetRightMargin(edge);
      // else pad[i][j]->SetRightMargin(0);
      //else pad[i][j]->SetRightMargin(0.03);

      // if(j==0) pad[i][j]->SetTopMargin(edge_top);
      if(j==0) pad[i][j]->SetTopMargin(bottomMargin+edge_top);
      else pad[i][j]->SetTopMargin(edge_top);
      //else pad[i][j]->SetTopMargin(edge_top);
      //else pad[i][j]->SetTopMargin(0.03);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin+edge_bottom);
      else pad[i][j]->SetBottomMargin(edge_bottom);
      //else pad[i][j]->SetBottomMargin(0);

      //pad[i][j]->SetFrameFillStyle(4000);
      pad[i][j]->SetFillColor(0);
      pad[i][j]->SetFillStyle(0);
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);

    }
  }
}



Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  //3.141592653589
  if ( dphi > TMath::Pi() )
    dphi = dphi - 2. * TMath::Pi();
  if ( dphi <= -TMath::Pi() )
    dphi = dphi + 2. * TMath::Pi();

  if ( TMath::Abs(dphi) > TMath::Pi() ) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
    std::cout << " " << phi1 << ", " << phi2 << ", " << dphi << std::endl;
  }

  return dphi;
}

Double_t cleverRange(TH1* h,TH1* h2, Float_t fac=1.2)
{
  Float_t maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  Float_t maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());

  Float_t minY1 =  (2.0-fac) * h->GetBinContent(h->GetMinimumBin());
  Float_t minY2 =  (2.0-fac) * h2->GetBinContent(h2->GetMinimumBin());
  //cout <<" range will be set as " << minY1 << " ~ " << minY2 << endl;
  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
  h->SetAxisRange(min(minY1,minY2),max(maxY1,maxY2),"Y");
  h2->SetAxisRange(min(minY1,minY2),max(maxY1,maxY2),"Y");
  return min(minY1,minY2);
  //return max(maxY1,maxY2);
}


Double_t cleverRange_Zero(TH1* h,TH1* h2, TH1* h3, Float_t fac=1.2)
{
  Float_t maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  Float_t maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());
  Float_t maxY3 =  fac * h3->GetBinContent(h3->GetMaximumBin());

  //Float_t minY1 =  (2.0-fac) * h->GetBinContent(h->GetMinimumBin());
  //Float_t minY2 =  (2.0-fac) * h2->GetBinContent(h2->GetMinimumBin());
  //Float_t minY3 =  (2.0-fac) * h3->GetBinContent(h3->GetMinimumBin());
  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;                             
  //Float_t firstmin = min(minY1,minY2);
  Float_t firstmax = max(maxY1,maxY2);
  Float_t finalmin = 0; //min(firstmin,minY3);
  Float_t finalmax = max(firstmax,maxY3);
  h->SetAxisRange(finalmin,finalmax,"Y");
  h2->SetAxisRange(finalmin,finalmax,"Y");
  h3->SetAxisRange(finalmin,finalmax,"Y");
  return finalmin;
}

Double_t cleverRange(TH1* h,TH1* h2, TH1* h3, Float_t fac=1.2)
{
  Float_t maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  Float_t maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());
  Float_t maxY3 =  fac * h3->GetBinContent(h3->GetMaximumBin());

  Float_t minY1 =  (2.0-fac) * h->GetBinContent(h->GetMinimumBin());
  Float_t minY2 =  (2.0-fac) * h2->GetBinContent(h2->GetMinimumBin());
  Float_t minY3 =  (2.0-fac) * h3->GetBinContent(h3->GetMinimumBin());
  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
  Float_t firstmin = min(minY1,minY2);
  Float_t firstmax = max(maxY1,maxY2);
  Float_t finalmin = min(firstmin,minY3);
  Float_t finalmax = max(firstmax,maxY3);
  h->SetAxisRange(finalmin,finalmax,"Y");
  h2->SetAxisRange(finalmin,finalmax,"Y");
  h3->SetAxisRange(finalmin,finalmax,"Y");
  return finalmin;
}


void SetHistColor(TH1* h, Int_t color=1)
{
  h->SetMarkerColor(color);
  h->SetLineColor(color);
}

float mean(float data[], int n)
{
  float mean=0.0;
  int i;
  for(i=0; i<n;++i)
  {
    mean+=data[i];
  }
  mean=mean/n;
  return mean;
}


void saveHistogramsToPicture(TH1* h, const char* fileType="pdf", const char* caption="", const char* directoryToBeSavedIn="figures", const char* text = "", int styleIndex=1, int rebin =1){
  TCanvas* c1=new TCanvas();
  if(rebin!=1)
  {
    h->Rebin(rebin);
  }

  if(styleIndex==1)
  {
    h->Draw("E");
  }
  else
  {
    h->Draw();
    if(h->InheritsFrom("TH2"))
    {
      h->Draw("COLZ TEXT");    // default plot style for TH2 histograms
    }
  }
  drawText(text,0.7,0.7);
  if(strcmp(directoryToBeSavedIn, "") == 0)   // save in the current directory if no directory is specified
  {
    c1->SaveAs(Form("%s_%s.%s" ,h->GetName(),caption, fileType));  // name of the file is the name of the histogram
  }
  else
  {
    c1->SaveAs(Form("%s/%s_%s.%s", directoryToBeSavedIn ,h->GetName(),caption, fileType));
  }
  c1->Close();
}


TGraphAsymmErrors* scale_graph(TGraphAsymmErrors* gr=0, Float_t s=0){
    int np = gr->GetN();
    TGraphAsymmErrors* new_gr = (TGraphAsymmErrors*) gr->Clone(Form("%s_scaled",gr->GetName()));
    //TGraphAsymmErrors* new_gr = new TGraphAsymmErrors(np);
    new_gr->SetName(Form("%s_scaled",gr->GetName()));
    for (int p=0; p<np; ++p) {
        double Xmean = 0;
        double Ymean = 0;
        gr->GetPoint(p,Xmean,Ymean);
        new_gr->SetPoint(p,Xmean,Ymean*s);
        double Xerr_l= gr->GetErrorXlow(p);
        double Xerr_h= gr->GetErrorXhigh(p);
        double Yerr_l= gr->GetErrorYlow(p);
        double Yerr_h= gr->GetErrorYhigh(p);
        new_gr->SetPointError(p,Xerr_l,Xerr_h,Yerr_l*s,Yerr_h*s);
    }
    return new_gr;
}

TGraphAsymmErrors* divide_graph_by_hist(TGraphAsymmErrors* gr=0, TH1* h1=0){
    int np = gr->GetN();
    TGraphAsymmErrors* new_gr = (TGraphAsymmErrors*) gr->Clone(Form("%s_dividedBy_%s",gr->GetName(),h1->GetName()));
    //TGraphAsymmErrors* new_gr = new TGraphAsymmErrors(np);
    new_gr->SetName(Form("%s_dividedBy_%s",gr->GetName(),h1->GetName()));
    for (int p=0; p<np; ++p) {
        double Xmean = 0;
        double Ymean = 0;
        double scaleF = h1->GetBinContent(p+1);
        gr->GetPoint(p,Xmean,Ymean);
        new_gr->SetPoint(p,Xmean,Ymean/scaleF);
        double Xerr_l= gr->GetErrorXlow(p);
        double Xerr_h= gr->GetErrorXhigh(p);
        double Yerr_l= gr->GetErrorYlow(p);
        double Yerr_h= gr->GetErrorYhigh(p);
        new_gr->SetPointError(p,Xerr_l,Xerr_h,Yerr_l/scaleF,Yerr_h/scaleF);
    }
    return new_gr;
}

TGraphAsymmErrors* divide_graph_by_graph(TGraph* gr_num=0, TGraph* gr_den=0){
    int np = gr_den->GetN();
    TGraphAsymmErrors* new_gr = (TGraphAsymmErrors*) gr_den->Clone(Form("%s_dividedBy_%s",gr_num->GetName(),gr_den->GetName()));
    for (int p=0; p<np; ++p) {
        double Xmean_den = 0;
        double Ymean_den = 0;
        double Xmean_num = 0;
        double Ymean_num = 0;
        //double scaleF = h1->GetBinContent(p+1);
        gr_den->GetPoint(p,Xmean_den,Ymean_den);
        gr_num->GetPoint(p,Xmean_num,Ymean_num);
        new_gr->SetPoint(p,Xmean_den,Ymean_num/Ymean_den);
        //double Xerr_l= gr->GetErrorXlow(p);
        //double Xerr_h= gr->GetErrorXhigh(p);
        //double Yerr_l= gr->GetErrorYlow(p);
        //double Yerr_h= gr->GetErrorYhigh(p);
        //new_gr->SetPointError(p,Xerr_l,Xerr_h,Yerr_l/scaleF,Yerr_h/scaleF);
    }
    return new_gr;
}

void hist_to_graph(TGraphAsymmErrors* gr=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, bool doSelfScale=1){
    int np = gr->GetN();
    for (int p=0; p<np; ++p) {
        double Xmean = h1->GetBinCenter(p+1);
        double bin_width = h1->GetBinLowEdge(p+2) - h1->GetBinLowEdge(p+1);
        double nominal = h1->GetBinContent(p+1);
        double sysDown = h2->GetBinContent(p+1);
        double sysUp = h3->GetBinContent(p+1);
        double dy = abs(nominal-sysDown);
        if(dy<abs(sysUp-nominal)) dy = abs(sysUp-nominal);
        //cout << " nominal = " << nominal << ", sysDown = " << sysDown << ", sysUp = " << sysUp << endl;
        double dy_up = abs(nominal-sysUp);
        double dy_down = abs(nominal-sysDown);

        if(doSelfScale){
            gr->SetPoint(p,Xmean,nominal/nominal);
            gr->SetPointError(p,(bin_width)/2.,(bin_width)/2.,dy_up/nominal,dy_down/nominal);
        } else{
            gr->SetPoint(p,Xmean,nominal);
            gr->SetPointError(p,(bin_width)/2.,(bin_width)/2.,dy_up,dy_down);
        }
    }
}

void MakeTextPrint(vector<string> lines, double x, double y, double dy){

  int s = lines.size();
  for (int i = 0; i < s; i++){
    drawText(lines.at(i).c_str(), x, y - i*dy, 0, kBlack, 0.03, 42);
  }
}
#endif
