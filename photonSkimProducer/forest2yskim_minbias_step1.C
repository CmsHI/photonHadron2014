///////////////////////////////////////////////////////////////////                                
// forest2yskim.C                                                //                                                 
// Creator : Yongsun Kim (MIT), jazzitup@mit.edu                 //                                                 
// Function : Transform hiForest files into yskim file           //
// yskims for MinBias1, Minbias2 and photon jet skims            //
///////////////////////////////////////////////////////////////////         
//d
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
#include <TRandom3.h>

using namespace std;

static const long MAXTREESIZE = 10000000000;





void forest2yskim_minbias_tracks(TString inputFile_="forestFiles/HiForest4/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21_0.root",
				   sampleType colli=kHIDATA,
				   int maxEvent = -1
				   )
{ 
  
  bool isMC=true;
  if ((colli==kPPDATA)||(colli==kPADATA)||(colli==kHIDATA))
    isMC=false;
  
   int seconds = time(NULL);  
  cout << " time = " <<seconds%10000<< endl;
  TRandom3 rand(seconds%10000);

  TString sampleString = "kPPDATA";
  if (colli==kPADATA) sampleString = "kPADATA";
  if (colli==kHIDATA) sampleString = "kHIDATA";
  if (colli==kPPMC) sampleString = "kPPMC";
  if (colli==kPAMC) sampleString = "kPAMC";
  if (colli==kHIMC) sampleString = "kHIMC";


  TString outname =  inputFile_(0, inputFile_.Last('/')+1) +  "trackSkim_collId_"+ sampleString + "_"+inputFile_(inputFile_.Last('/')+1,200);
  HiForest *c;
  if  ((colli==kPADATA)||(colli==kPAMC))      c = new HiForest(inputFile_.Data(), "forest", cPPb, isMC);
  else if  ((colli==kPPDATA)||(colli==kPPMC)) c = new HiForest(inputFile_.Data(), "forest", cPP, isMC);
  else if  ((colli==kHIDATA)||(colli==kHIMC)) c = new HiForest(inputFile_.Data(), "forest", cPbPb, isMC);
  else {
    cout << " Error!  No such collision type" << endl;
    return;
  }
  c->InitTree();
  
  // Second forest files for track mixing 
  // output file
  TFile* newfile_data = new TFile(outname,"recreate");
   
  // Track tree retrieved on Feb 11 2014
  int nTrk;
  static const int MAXTRK  = 5000;   // This must be  enough.
  float trkPt[MAXTRK];
  float trkEta[MAXTRK];
  float trkPhi[MAXTRK];
  int   trkPurity[MAXTRK];
  int   trkAlgo[MAXTRK];

  EvtSel evt;
  TTree* newtreeTrkJet[200][nVtxBin+1];
  
  int nCentBins =  nCentBinSkim;
  if ((colli==kPADATA)||(colli==kPAMC)) {
    nCentBins = nCentBinSkimPA;
  }
  
  
  for( int icent = 0 ; icent< nCentBins ; icent++) { 
    for( int ivz = 1 ; ivz<=nVtxBin ; ivz++) {
      newtreeTrkJet[icent][ivz] = new TTree(Form("trkAndJets_first_icent%d_ivz%d",icent,ivz),"track and jets");
      newtreeTrkJet[icent][ivz]->SetMaxTreeSize(MAXTREESIZE);

      newtreeTrkJet[icent][ivz]->Branch("evt",&evt.run,evtLeaves.Data());
      newtreeTrkJet[icent][ivz]->Branch("nTrk",&nTrk,"nTrk/I");
      newtreeTrkJet[icent][ivz]->Branch("trkPt",trkPt,"trkPt[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkEta",trkEta,"trkEta[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkPhi",trkPhi,"trkPhi[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkPurity",trkPurity,"trkPurity[nTrk]/I");
      newtreeTrkJet[icent][ivz]->Branch("trkAlgo",trkAlgo,"trkAlgo[nTrk]/I");

    } 
  }
  
  // vertex histogram 
  float vzCut = vtxCutPhotonAna;
  TH1F* hvz = new TH1F("hvz","",nVtxBin,-vzCut,vzCut);
  // event plane hitogram
  TH1F* hEvtPlnBin = new TH1F("hEvtPlnBin", "", nPlnBin, -PI/2., PI/2.);

  
  /// LOOP!!
  int nentries = c->GetEntries();
  if ( maxEvent > 0 ) 
    nentries = maxEvent;
  cout << "number of entries = " << nentries << endl;
  
  for (Long64_t jentry=0 ; jentry<nentries;jentry++) {
    if (jentry% 100 == 0)  {
      cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
    }
    c->GetEntry(jentry);
    evt.clear();
    evt.run   = c->evt.run;
    evt.evt = c->evt.evt;
    evt.hf4Pos = c->evt.hiHFplusEta4;
    evt.hf4Neg = c->evt.hiHFminusEta4;
    evt.hf4Sum = evt.hf4Pos + evt.hf4Neg;
    evt.cBin = -99;
    evt.pBin   = -99 ;
    if ((colli==kHIDATA)||(colli==kHIMC)) {
      evt.cBin = getCbinFrom200(c->evt.hiBin);
      evt.pBin   = hEvtPlnBin->FindBin( c->evt.hiEvtPlanes[theEvtPlNumber] ) ;
    }
    else if ((colli==kPADATA)||(colli==kPAMC))   {
      evt.cBin =  getHfBin(evt.hf4Sum);
      if (  ((evt.cBin) < 0) || (evt.cBin) >= nCentBinSkimPA )  
	cout << " Check the pA centrality..  cbin = " << evt.cBin << endl;
    }
    
    evt.vtxCentWeight = 1;
    evt.vz = c->evt.vz;
    
    // Event selection

    if ( ( (colli==kHIDATA)||(colli==kHIMC)||(colli==kPADATA)||(colli==kPAMC) || (colli==kPPMC) ) && ( c->selectEvent() == 0 ))
      continue;
    //    if ( ( (colli==kPADATA)||(colli==kPPDATA) ) && ( c->skim.pVertexFilterCutGplus ==0 ) ) // No Pile up events But it's not on the forest file anymore (Jul 4 2014)
    //    continue;
       
    // vertex bin and cut!! 
    int vzBin = hvz->FindBin(evt.vz)  ;
    hvz->Fill(evt.vz)  ;
    if ( (vzBin<1) || ( vzBin > nVtxBin) ) 
      continue;
    
    ///////////////////////////// Tracks //////////////////////////////
    nTrk = 0; 
    for (int it=0; it < c->track.nTrk; it++ ) { 
      if ( c->track.trkPt[it] < cuttrkPtSkim )   continue;
      if (  fabs(c->track.trkEta[it]) > cuttrkEtaSkim ) continue;
      if ( c->selectTrack(it)  == false ) continue;
      trkPt[nTrk]  = c->track.trkPt[it];
      trkEta[nTrk] = c->track.trkEta[it];
      trkPhi[nTrk] = c->track.trkPhi[it]; 
      trkPurity[nTrk] = c->track.highPurity[it];
      trkAlgo[nTrk] = c->track.trkAlgo[it];
      
      nTrk++;
    }
    
    newtreeTrkJet[evt.cBin][vzBin]->Fill();
  }
  newfile_data->Write();
  cout << " Done! "<< endl;
}



