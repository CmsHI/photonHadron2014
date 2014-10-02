///////////////////////////////////////////////////////////////////                                
// makeTrialForestFiles.C                                        //
// Creator : Yongsun Kim (KU), kimy@cern.ch                      //
// Function : Extract a small portion of a Forest tree           //
///////////////////////////////////////////////////////////////////         

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <TMath.h>
#include "../../HiForestAnalysisPostHP/hiForest.h"
#include "../CutAndBinCollection2012.h"
#include <time.h>


using namespace std;

static const long MAXTREESIZE = 10000000000;





void makeTrialForestFile(int nMax= 500){ 
  //collisionType
  collisionType cMode = cPbPb;
  TString inputFile_="forestFiles/0.root";
  TString outname = "0_1000files_pbpb_minbiasData.root";
  
  
  // start from here
  // path length histogram
  
  HiForest * t;
  t = new HiForest(inputFile_.Data(),"",cMode);
   
  t->SetOutputFile(outname.Data());
  // LOOP!!
  t->InitTree();

  for (Long64_t jentry = 0 ; jentry < nMax ; jentry++) {
    if (jentry% 1000 == 0)  {
      cout <<jentry<<" / "<<nMax<<" "<<setprecision(2)<<(double)jentry/nMax*100<<endl;
    }
    t->GetEntry(jentry);
    t->FillOutput();
  }
  delete t;
}

