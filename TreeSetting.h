#ifndef TreeSetting_h
#define TreeSetting_h

#include <iostream>
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

class SetTree
{ 
  public:
    SetTree(){};

    virtual ~SetTree();
    virtual void TreeSetting(TTree* tree);
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

#endif
