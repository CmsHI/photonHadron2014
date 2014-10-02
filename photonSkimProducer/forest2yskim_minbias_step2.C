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





void forest2yskim_minbias_forestV3(TString inputFile_="forestFiles/HiForest4/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21_0.root",
				   TString inputFileTrack="forestFiles/HiForest4/trackSkim_collId_kHIDATA_HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21_0.root",
				   sampleType colli=kHIDATA,
				   int maxEvent = 1000,
				   bool useGenJetColl = 0,
				   TString jetAlgo="akPu3PF"
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


  TString outname =  inputFile_(0, inputFile_.Last('/')+1) +  "skim_collId_"+ sampleString + "_jetAlgo_"+jetAlgo+"_"+inputFile_(inputFile_.Last('/')+1,200);


  
  
  HiForest *c;
  if  ((colli==kPADATA)||(colli==kPAMC))      c = new HiForest(inputFile_.Data(), "forest", cPPb, isMC);
  else if  ((colli==kPPDATA)||(colli==kPPMC)) c = new HiForest(inputFile_.Data(), "forest", cPP, isMC);
  else if  ((colli==kHIDATA)||(colli==kHIMC)) c = new HiForest(inputFile_.Data(), "forest", cPbPb, isMC);
  else {
    cout << " Error!  No such collision type" << endl;
    return;
  }
  c->InitTree();
  
  
  // output file
  TFile* newfile_data = new TFile(outname,"recreate");
   
  // Track tree retrieved on Feb 11 2014
  int nTrk;
  static const int MAXTRK  = 10000;   // This must be  enough.
  float trkPt[MAXTRK];
  float trkEta[MAXTRK];
  float trkPhi[MAXTRK];
  int   trkPurity[MAXTRK];
  int    trkAlgo[MAXTRK];
  float trkAsJetPt[MAXTRK];  // associated Jet pT 
  float trkAsJetEta[MAXTRK];  // associated Jet pT 
  float trkAsJetPhi[MAXTRK];  // associated Jet pT 
  float trkAsJetDR[MAXTRK];  // associated Jet pT 

  int nMTrk;
  float mTrkPt[MAXTRK];
  float mTrkEta[MAXTRK];
  float mTrkPhi[MAXTRK];
  int  mTrkPurity[MAXTRK];
  int   mTrkAlgo[MAXTRK];
  float mTrkAsJetPt[MAXTRK];  // associated Jet pT 
  float mTrkAsJetEta[MAXTRK];  // associated Jet pT 
  float mTrkAsJetPhi[MAXTRK];  // associated Jet pT 
  float mTrkAsJetDR[MAXTRK];  // associated Jet pT 

  // Jet tree
  int nJet;
  static const int MAXJET  = 200;   // This must be  enough.
  float jetPt[MAXJET];
  float jetEta[MAXJET];
  float jetPhi[MAXJET];
  int jetSubid[MAXJET];
  float jetRefPt[MAXJET];
  float jetRefEta[MAXJET];
  float jetRefPhi[MAXJET];
  float jetRefPartonPt[MAXJET];
  int  jetRefPartonFlv[MAXJET];

  // genparticle tree
  
  const int MAXCh = 10000;
  int nCh;
  int chChg[MAXCh];
  int chPdg[MAXCh];
  float chPt[MAXCh];
  float chEta[MAXCh];
  float chPhi[MAXCh];


  EvtSel evt;
  TTree* newtreeTrkJet[200][nVtxBin+1];
  
  int nCentBins =  nCentBinSkim;
  if ((colli==kPADATA)||(colli==kPAMC)) {
    nCentBins = nCentBinSkimPA;
  }
  
    
  for( int icent = 0 ; icent< nCentBins ; icent++) { 
    for( int ivz = 1 ; ivz<=nVtxBin ; ivz++) {
      newtreeTrkJet[icent][ivz] = new TTree(Form("trkAndJets_second_icent%d_ivz%d",icent,ivz),"track and jets");
      newtreeTrkJet[icent][ivz]->SetMaxTreeSize(MAXTREESIZE);
      newtreeTrkJet[icent][ivz]->Branch("evt",&evt.run,evtLeaves.Data());
      
      // jets
      newtreeTrkJet[icent][ivz]->Branch("nJet",&nJet,"nJet/I");
      newtreeTrkJet[icent][ivz]->Branch("jetPt",jetPt,"jetPt[nJet]/F");
      newtreeTrkJet[icent][ivz]->Branch("jetEta",jetEta,"jetEta[nJet]/F");
      newtreeTrkJet[icent][ivz]->Branch("jetPhi",jetPhi,"jetPhi[nJet]/F");
      if ( isMC )  {
	newtreeTrkJet[icent][ivz]->Branch("subid",jetSubid,"subid[nJet]/I");
	newtreeTrkJet[icent][ivz]->Branch("refPt",jetRefPt,"refPt[nJet]/F");
	newtreeTrkJet[icent][ivz]->Branch("refEta",jetRefEta,"refEta[nJet]/F");
	newtreeTrkJet[icent][ivz]->Branch("refPhi",jetRefPhi,"refPhi[nJet]/F");
	newtreeTrkJet[icent][ivz]->Branch("refPartonPt",jetRefPartonPt,"refPartonPt[nJet]/F");
	newtreeTrkJet[icent][ivz]->Branch("refPartonFlv",jetRefPartonFlv,"refPartonFlv[nJet]/I");
      }
      
      // charged particles
      if ( isMC) { 
	newtreeTrkJet[icent][ivz]->Branch("nCh",&nCh,"nCh/I");
	newtreeTrkJet[icent][ivz]->Branch("chPt",chPt,"chPt[nCh]/F");
	newtreeTrkJet[icent][ivz]->Branch("chEta",chEta,"chEta[nCh]/F");
	newtreeTrkJet[icent][ivz]->Branch("chPhi",chPhi,"chPhi[nCh]/F");
	newtreeTrkJet[icent][ivz]->Branch("chPdg",chPdg,"chPdg[nCh]/I");
	newtreeTrkJet[icent][ivz]->Branch("chChg",chChg,"chChg[nCh]/I");
      }
      
      // tracks
      newtreeTrkJet[icent][ivz]->Branch("nTrk",&nTrk,"nTrk/I");
      newtreeTrkJet[icent][ivz]->Branch("trkPt",trkPt,"trkPt[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkEta",trkEta,"trkEta[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkPhi",trkPhi,"trkPhi[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkPurity",trkPurity,"trkPurity[nTrk]/I");
      newtreeTrkJet[icent][ivz]->Branch("trkAlgo",trkAlgo,"trkAlgo[nTrk]/I");
      newtreeTrkJet[icent][ivz]->Branch("trkAsJetPt",trkAsJetPt,"trkAsJetPt[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkAsJetEta",trkAsJetEta,"trkAsJetEta[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkAsJetPhi",trkAsJetPhi,"trkAsJetPhi[nTrk]/F");
      newtreeTrkJet[icent][ivz]->Branch("trkAsJetDR",trkAsJetDR,"trkAsJetDR[nTrk]/F");
      
      if ((colli==kHIDATA)||(colli==kHIMC))    { 
	newtreeTrkJet[icent][ivz]->Branch("nMTrk",&nMTrk,"nMTrk/I");
	newtreeTrkJet[icent][ivz]->Branch("mTrkPt",mTrkPt,"mTrkPt[nMTrk]/F");
	newtreeTrkJet[icent][ivz]->Branch("mTrkEta",mTrkEta,"mTrkEta[nMTrk]/F");
	newtreeTrkJet[icent][ivz]->Branch("mTrkPhi",mTrkPhi,"mTrkPhi[nMTrk]/F");
	newtreeTrkJet[icent][ivz]->Branch("mTrkPurity",mTrkPurity,"mTrkPurity[nMTrk]/I");
	newtreeTrkJet[icent][ivz]->Branch("mTrkAlgo",mTrkAlgo,"mTrkAlgo[nMTrk]/I");
	newtreeTrkJet[icent][ivz]->Branch("mTrkAsJetPt",mTrkAsJetPt,"mTrkAsJetPt[nMTrk]/F");
	newtreeTrkJet[icent][ivz]->Branch("mTrkAsJetEta",mTrkAsJetEta,"mTrkAsJetEta[nMTrk]/F");
	newtreeTrkJet[icent][ivz]->Branch("mTrkAsJetPhi",mTrkAsJetPhi,"mTrkAsJetPhi[nMTrk]/F");
	newtreeTrkJet[icent][ivz]->Branch("mTrkAsJetDR",mTrkAsJetDR,"mTrkAsJetDR[nMTrk]/F");
      }
    }
  }

  // tracks for the seoncd minbias track skim.
  EvtSel evtImb;
  TBranch        *b_evt;
  Int_t           nTrkImb;
  Float_t         trkPtImb[MAXTRK];
  Float_t         trkEtaImb[MAXTRK];
  Float_t         trkPhiImb[MAXTRK];
  int          trkPurityImb[MAXTRK];
  int          trkAlgoImb[MAXTRK];
  TBranch        *b_nTrkImb;
  TBranch        *b_trkPtImb;
  TBranch        *b_trkEtaImb;
  TBranch        *b_trkPhiImb;
  TBranch        *b_trkPurityImb;
  TBranch        *b_trkAlgoImb;
  
  TChain   *tjmbImb[200][nVtxBin+1];
  int nMB[200][nVtxBin+1];
  int mbItr[200][nVtxBin+1];
  
  if ( inputFileTrack != "" ) { 
    for( int icent = 0 ; icent< nCentBins ; icent++) {
      for( int ivz = 1 ; ivz<=nVtxBin ; ivz++) {
	tjmbImb[icent][ivz] = new TChain(Form("trkAndJets_first_icent%d_ivz%d",icent,ivz));
	tjmbImb[icent][ivz]->Add(inputFileTrack.Data());
	tjmbImb[icent][ivz]->SetBranchAddress("evt", &evtImb,&b_evt);
	
	tjmbImb[icent][ivz]->SetBranchAddress("nTrk",   &nTrkImb,   &b_nTrkImb);
	tjmbImb[icent][ivz]->SetBranchAddress("trkPt",  &trkPtImb,  &b_trkPtImb);
	tjmbImb[icent][ivz]->SetBranchAddress("trkEta", &trkEtaImb, &b_trkEtaImb);
	tjmbImb[icent][ivz]->SetBranchAddress("trkPhi", &trkPhiImb, &b_trkPhiImb);
	tjmbImb[icent][ivz]->SetBranchAddress("trkPurity", &trkPurityImb, &b_trkPurityImb);
	tjmbImb[icent][ivz]->SetBranchAddress("trkAlgo", &trkAlgoImb, &b_trkAlgoImb);

	nMB[icent][ivz] = tjmbImb[icent][ivz]->GetEntries();
	cout << "number of evetns in (icent = " << icent << ", ivtxZ = "<< ivz << ")  = " << nMB[icent][ivz] << endl;
        int primeSeed = rand.Integer(37359); // 37357 is an arbitrary prime number set by YS March 17th 2014
        mbItr[icent][ivz] = primeSeed%(nMB[icent][ivz]);
        cout <<" initial itr = " << mbItr[icent][ivz] << endl;
	tjmbImb[icent][ivz]->GetEntry(mbItr[icent][ivz]);
      }
    }
  }
  // vertex histogram 
  float vzCut = vtxCutPhotonAna;
  TH1F* hvz = new TH1F("hvz","",nVtxBin,-vzCut,vzCut);
  // event plane hitogram
  TH1F* hEvtPlnBin = new TH1F("hEvtPlnBin", "", nPlnBin, -PI/2., PI/2.);
  // jet algos  
             
  Jets* theJet;
  if ( jetAlgo == "akPu3PF")  {
    theJet = &(c->akPu3PF) ;   cout << "Using akPu3PF Jet Algo" << endl<<endl;
  }
  else if ( jetAlgo == "akVs3PF") {
    theJet = &(c->akVs3PF) ;   cout << "Using ak3PF Jet Algo, Voronoi Subtraction method" << endl<<endl;
  } 
  
  
  /// LOOP!!
  int nentries = c->GetEntries();
  //  int nentries = 300;
  if ( maxEvent > 0 ) 
    nentries = maxEvent;
  cout << "number of entries = " << nentries << endl;
  
  for (Long64_t jentry=0 ; jentry<nentries;jentry++) {
    if (jentry% 1000 == 0)  {
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
    if ((colli==kHIDATA) || (colli==kHIMC)) {
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
    
    int cBin = evt.cBin; 
    // vertex bin and cut!! 
    int vzBin = hvz->FindBin(evt.vz)  ;
    hvz->Fill(evt.vz)  ;
    
    if ( ( (colli==kHIDATA)||(colli==kHIMC)||(colli==kPADATA)||(colli==kPAMC) || (colli==kPPMC) ) && ( c->selectEvent() == 0 ))
      continue;
    if ( (vzBin<1) || ( vzBin > nVtxBin) ) 
      continue;
    
    ///////////// Collection of jets  /////////////////////////////////////////////
    nJet = 0 ;
    int jetEntries = 0;
    if (useGenJetColl )    jetEntries = theJet->ngen;
    else                   jetEntries = theJet->nref;

    for (int ij=0; ij< jetEntries ; ij++) {
      if (  useGenJetColl )   {   
	jetPt[nJet] = theJet->genpt[ij];
	jetEta[nJet] = theJet->geneta[ij];
	jetPhi[nJet] = theJet->genphi[ij];
      }
      else  {
	jetPt[nJet] = theJet->jtpt[ij];
	jetEta[nJet] = theJet->jteta[ij];
	jetPhi[nJet] = theJet->jtphi[ij];
      }


      if ( jetPt[nJet] < cutjetPtSkim)
	continue;
      if ( fabs( jetEta[nJet] ) > cutjetEtaSkim )
        continue;
      
      if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
        jetEta[nJet] = -jetEta[nJet] ; 
      }
      if (  useGenJetColl )   {
	jetSubid[nJet] = -9999;
	jetRefPt[nJet] = -9999;
	jetRefEta[nJet] = -9999;
	jetRefPhi[nJet] = -9999;
	jetRefPartonPt[nJet] = -9999;
	jetRefPartonFlv[nJet] = -9999;
      }
      else { 
	jetSubid[nJet] = theJet->subid[ij];
	jetRefPt[nJet] = theJet->refpt[ij];
	jetRefEta[nJet] = theJet->refeta[ij];
	jetRefPhi[nJet] = theJet->refphi[ij];
	jetRefPartonPt[nJet] = theJet->refparton_pt[ij];
	jetRefPartonFlv[nJet] = theJet->refparton_flavor[ij];
      }
      
      nJet++ ;
    }
    
    // charged particles 
    nCh = 0;
    if ( isMC) {
      for (int it=0; it< c->genparticle.mult ; it++) {
        if ( c->genparticle.pt[it] < cuttrkPtSkim )   continue;
        if (  fabs(c->genparticle.eta[it]) > cuttrkEtaSkim ) continue;
        if ( c->genparticle.chg[it] == 0 ) continue;
        //      if ( c->genparticle.sube[it] != 0 ) continue;
	chPdg[nCh]  = c->genparticle.pdg[it];
        chChg[nCh]  = c->genparticle.chg[it];
        chPt[nCh]  = c->genparticle.pt[it];
        chEta[nCh]  = c->genparticle.eta[it];
        chPhi[nCh]  = c->genparticle.phi[it];
        nCh++;}
    }
    
    ///////////////////////////// Tracks //////////////////////////////
    nTrk = 0; 
    for (int it=0; it < c->track.nTrk; it++ ) { 
      if ( c->track.trkPt[it] < cuttrkPtSkim )   continue;
      if (  fabs(c->track.trkEta[it]) > cuttrkEtaSkim ) continue;
      if ( c->selectTrack(it)  == false) continue;
      trkPt[nTrk]  = c->track.trkPt[it];
      trkEta[nTrk] = c->track.trkEta[it];
      trkPhi[nTrk] = c->track.trkPhi[it]; 
      trkPurity[nTrk] = c->track.highPurity[it]; 
      trkAlgo[nTrk] = c->track.trkAlgo[it]; 
      //  trkWeight[nTrk] = c->getTrackCorrection(it);
      int assocJetId = matchedJetFinder( theJet, trkEta[nTrk], trkPhi[nTrk]);
      if ( assocJetId < 0 )  {
	trkAsJetPt[nTrk] = -1; 
	trkAsJetEta[nTrk] = -1; 
	trkAsJetPhi[nTrk] = -1; 
	trkAsJetDR[nTrk] = 100; 
      } 
      else { 
	trkAsJetPt[nTrk] = theJet->jtpt[assocJetId];
	trkAsJetEta[nTrk] = theJet->jteta[assocJetId];
	trkAsJetPhi[nTrk] = theJet->jtphi[assocJetId];
	trkAsJetDR[nTrk] =getDR( trkEta[nTrk], trkPhi[nTrk], theJet->jteta[assocJetId], theJet->jtphi[assocJetId]) ;
      }
      
      nTrk++;
    }
    
    ///////////////////////////// background tracks //////////////////////////////
    // Step 1.  Find the centrality matched events 
    if ( inputFileTrack != "" )  {   
      bool findFlag = false;
      for ( int ie = 1 ; ie<=nMB[cBin][vzBin] ; ie++) { 
	mbItr[cBin][vzBin] = mbItr[cBin][vzBin]++;
	if ( mbItr[cBin][vzBin] == nMB[cBin][vzBin] ) 
	  mbItr[cBin][vzBin] = mbItr[cBin][vzBin] -  nMB[cBin][vzBin] ; 
	tjmbImb[cBin][vzBin]->GetEntry(mbItr[cBin][vzBin]);
	if ( c->evt.evt == evtImb.evt ) // Identical event
	  continue;
	// OK we found the matched event
	findFlag = true;
	break;
      }
      if (!findFlag) { 
	cout << " WARNING!!! Could not find the centrlaity matched events" << endl;
	return;
      }
      nMTrk = 0; 
      for (int it=0; it < nTrkImb ; it++ ) { 
	if ( trkPtImb[it] < cuttrkPtSkim )   continue;
	if (  fabs(trkEtaImb[it]) > cuttrkEtaSkim ) continue;
	mTrkPt[nMTrk]  = trkPtImb[it];
	mTrkEta[nMTrk] = trkEtaImb[it];
	mTrkPhi[nMTrk] = trkPhiImb[it];
	mTrkPurity[nMTrk] = trkPurityImb[it];
	mTrkAlgo[nMTrk] = trkAlgoImb[it];
	int assocJetId = matchedJetFinder( theJet, mTrkEta[nMTrk], mTrkPhi[nMTrk]);
	if ( assocJetId < 0 )  {
	  mTrkAsJetPt[nMTrk] = -1; 
	  mTrkAsJetEta[nMTrk] = -1; 
	  mTrkAsJetPhi[nMTrk] = -1; 
	  mTrkAsJetDR[nMTrk] = 100; 
	} 
	else { 
	  mTrkAsJetPt[nMTrk] = theJet->jtpt[assocJetId];
	  mTrkAsJetEta[nMTrk] = theJet->jteta[assocJetId];
	  mTrkAsJetPhi[nMTrk] = theJet->jtphi[assocJetId];
	  mTrkAsJetDR[nMTrk] =getDR( mTrkEta[nMTrk], mTrkPhi[nMTrk], theJet->jteta[assocJetId], theJet->jtphi[assocJetId]) ;
	}
	nMTrk++;
      }
    }
    
    newtreeTrkJet[cBin][vzBin]->Fill();
  }
  
  newfile_data->Write();
  cout << " Done! "<< endl;  
}

