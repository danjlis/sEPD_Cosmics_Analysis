#ifndef TreeSetting_h
#define TreeSetting_h

#include <iostream>
#include <sstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TChain.h"

using namespace std;

const long int maxBranchSize = 100000;


// import the tree to the RooDataSet


vector<double>* pmpRef;
vector<double>* pmpTest;
vector<bool>* isCoincOld;
vector<bool>* isIsoOld;
const int nCh=16;
double adcpedcutref[nCh] = {800, 900, 1050, 1070, 
  1000, 1000, 950, 1230, 
  1210, 1130, 1110, 1150, 
  1100, 1050, 1100, 1100};

double adcpedcuttest[nCh] = {800, 900, 1050, 1070, 
  1000, 1000, 950, 1230, 
  1210, 1130, 1110, 1150, 
  1100, 1050, 1100, 1100
};

double xmaxPerc=1.0;

class SetTree
{ 
  public:
    SetTree(){};

    virtual ~SetTree();
    virtual void TreeSetting(TTree* tree);
    Bool_t ApplyPedCutRef(int i);
    Bool_t ApplyPedCutTest(int i);
    Bool_t IsoRef(int i);
    Bool_t IsoTest(int i);
};

SetTree::~SetTree()
{
}


void SetTree::TreeSetting(TTree* tree)
{
  tree->SetBranchAddress("pmpRef", &pmpRef);
  tree->SetBranchAddress("pmpTest", &pmpTest);
  tree->SetBranchAddress("isCoincOld", &isCoincOld);
  tree->SetBranchAddress("isIsoOld", &isIsoOld);
};

Bool_t SetTree::ApplyPedCutRef(int i)
{
  return (pmpRef->at(i) > adcpedcutref[i]);
};

Bool_t SetTree::ApplyPedCutTest(int i)
{
    return (pmpTest->at(i) > adcpedcuttest[i]);
};

Bool_t SetTree::IsoRef(int i)
{
  if(i==0) return ( (pmpRef->at(i+1) < adcpedcutref[i+1]*xmaxPerc) ); 
  else if(i>0 && i<nCh-1) return ( (pmpRef->at(i+1) < adcpedcutref[i+1]*xmaxPerc) && (pmpRef->at(i-1) < adcpedcutref[i-1]*xmaxPerc) );
  else if(i==nCh-1) return ( (pmpRef->at(i-1) < adcpedcutref[i-1]*xmaxPerc) );
};

Bool_t SetTree::IsoTest(int i)
{
  if(i==0) return ( (pmpTest->at(i+1) < adcpedcuttest[i+1]*xmaxPerc) ); 
  else if(i>0 && i<nCh-1) return ( (pmpTest->at(i+1) < adcpedcuttest[i+1]*xmaxPerc) && (pmpTest->at(i-1) < adcpedcuttest[i-1]*xmaxPerc) );
  else if(i==nCh-1) return ( (pmpTest->at(i-1) < adcpedcuttest[i-1]*xmaxPerc) );
};

#endif
