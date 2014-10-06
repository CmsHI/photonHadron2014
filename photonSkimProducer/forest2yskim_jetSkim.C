///////////////////////////////////////////////////////////////////
// forest2yskim.C                                                //
// Creator : Yongsun Kim (MIT), jazzitup@mit.edu                 //
// Function : Transform hiForest files into yskim file           //
// yskims for MinBias1, Minbias2 and photon jet skims            //
///////////////////////////////////////////////////////////////////
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
#include <TRandom3.h>

using namespace std;

static const int MAXTRK  = 10000;   // This is very enough.
static const int MAXMTRK = 80000;   // Again this is very enough for 10 mixing
static const int MAXMJET = 2000;
static const long MAXTREESIZE = 50000000000;



vector<jetKinem> nullVec;

double getPhiSmear( int smearingCentBin, float jetPt);
double getPtSmear(  int smearingCentBin, float jetPt);
float  getHiResCorr( EvtSel evt, float jetPt);
float  getL2L3Corr(EvtSel evt,  sampleType colli, float jetPt, float jetEta, TF1*  fgaus);

float getPpResCorr( float jetPt=0);
float normalAngle(float phi=0 );


// vertex and centrality vtxCentWeighting
TFile* fWeight = new TFile("../vertexReweightingHistogram_pthatweighted.root");

TH1D* hWeight_vtx_data_pp = (TH1D*)fWeight->Get("vertexHistoData_pp");
TH1D* hWeight_vtx_mc_pp = (TH1D*)fWeight->Get("vertexHistoMC_pp");

TH1D* hWeight_vtx_data_ppb = (TH1D*)fWeight->Get("vertexHistoData_ppb");
TH1D* hWeight_vtx_mc_ppb = (TH1D*)fWeight->Get("vertexHistoMC_ppb");

TH1D* hWeight_vtx_data_pbpb = (TH1D*)fWeight->Get("vertexHistoData_pbpb");
TH1D* hWeight_vtx_mc_pbpb = (TH1D*)fWeight->Get("vertexHistoMC_pbpb");
TH1D* hWeight_cent_data_pbpb = (TH1D*)fWeight->Get("centBinHistoData_pbpb");
TH1D* hWeight_cent_mc_pbpb = (TH1D*)fWeight->Get("centBinHistoMC_pbpb");

// L2L3 correction
TFile* fL2L3pp = new TFile("../corrL2L3/Casym_pp_double_hcalbins_algo_ak3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
TH1D * c_etapp=(TH1D*)fL2L3pp->Get("C_asym");
TF1* fptpp = new TF1("fptpp","1-[0]/pow(x,[1])",20,300);

TFile* fL2L3pA = new TFile("../corrL2L3/Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
TH1D * c_etapA=(TH1D*)fL2L3pA->Get("C_asym");
TF1* fptpA = new TF1("fptpA","1-[0]/pow(x,[1])",20,300);

TFile* fL2L3Ap = new TFile("../corrL2L3/Casym_Pbp_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root");
TH1D * c_etaAp=(TH1D*)fL2L3Ap->Get("C_asym");
TF1* fptAp = new TF1("fptAp","1-[0]/pow(x,[1])",20,300);

// pA MC
TF1 * fgaus=new TF1("fgaus_pA","gaus(0)",-20,20);
TF1 * fsmear_pA = new TF1("fsmear_pA","[0]/pow(x,[1])",50,300);

void drawL2L3Correciton(); 


void forest2yskim_jetSkim(TString inputFile_="forestFiles/HiForest4/hiForest_Photon40_GR_R_53_LV6_25Feb2014_1530CET_Track8_Jet15.root",
				   std::string MinbiasFname = "forestFiles/HiForest4/skim_collId_kHIDATA_jetAlgo_akPu3PF_HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21_0.root",
				   float cutphotonPt  = 35,  // default value dropped to 35GeV  for later photon energy smearing/scaling
				   sampleType colli=kHIDATA,
				   TString jetAlgo="akPu3PF",
				   bool doMix = true,
				   TString triggerSelection="",
				   bool doJetResCorrection = 0,  // = L2L3 * MC nonclosure correction  jet energy correction is done by default from Oct 19th (YS)
				   int smearingCentBin = -1, //0=0-10%, 1=10-30%, 2=30-50%, 3=50-100%, 4=0-30%, 5=30-100%  : Jet pT and phi smearing!
				   bool useGenJetColl = 0
				   )

{
  bool fillSecondMtrk = false;


  TString sampleString = "kPPDATA";
  if (colli==kPADATA) sampleString = "kPADATA";
  if (colli==kHIDATA) sampleString = "kHIDATA";
  if (colli==kPPMC) sampleString = "kPPMC";
  if (colli==kPAMC) sampleString = "kPAMC";
  if (colli==kHIMC) sampleString = "kHIMC";
  
  fgaus->SetParameters(1,0,1);
  fsmear_pA->SetParameters(1.052,0.5261);
  fptpp->SetParameters(0.06971,0.8167);
  fptpA->SetParameters(0.3015, 0.8913);
  fptAp->SetParameters(0.3015, 0.8913);


  bool isMC=true;
  if ((colli==kPPDATA)||(colli==kPADATA)||(colli==kHIDATA))
    isMC=false;

  int seconds = time(NULL);            //   cout << " time = " <<seconds%10000<< endl;
  TRandom3 rand(seconds%10000);
  TString datafname  = "";
  float cutphotonEta = 1.44;
  float preCutPhotonEt = 15;

  const int nMaxPho = 100;



  HiForest *c;
  if((colli==kPADATA)||(colli==kPAMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPPb, isMC );
  }
  else if  ((colli==kPPDATA)||(colli==kPPMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPP, isMC );
  }
  else if  ((colli==kHIDATA)||(colli==kHIMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPbPb, isMC );
    c->GetEnergyScaleTable("../photonEnergyScaleTable_lowPt_v6.root");
  }
  else {
    cout << " Error!  No such collision type" << endl;
    return;
  }
  
  c->InitTree();
  
  drawL2L3Correciton();

  // Create a new root file
  
  TString outname =  inputFile_(0, inputFile_.Last('/')+1) +  "skim_collId_" + sampleString + "_jetAlgo_" + jetAlgo + "_20mixed_" + inputFile_(inputFile_.Last('/')+1,200)  ;
  TFile* newfile_data = new TFile(outname.Data(),"recreate");

  TTree* newtreePhoton;
  float newPt[nMaxPho];  // <<= temporary space
  int order[nMaxPho];
  float corrPt[nMaxPho];
  newtreePhoton = c->photonTree->CloneTree(0);
  newtreePhoton->SetName("yPhotonTree");
  newtreePhoton->SetMaxTreeSize(MAXTREESIZE);
  newtreePhoton->Branch("order",  order, "order[nPhotons]/I");
  newtreePhoton->Branch("corrPt", corrPt,"corrPt[nPhotons]/F");
  
  // 1.1 jet tree!
  int nJet;
  const int MAXJET = 200; // to accomodate 100 smeared jets, need to be careful with ram
  float jetPt[MAXJET];
  float jetEta[MAXJET];
  float jetPhi[MAXJET];
  float jetDphi[MAXJET];
  int jetSubid[MAXJET];
  float jetRefPt[MAXJET];
  float jetRefEta[MAXJET];
  float jetRefPhi[MAXJET];
  float jetRefDphi[MAXJET];
  float jetRefPartonPt[MAXJET];
  int  jetRefPartonFlv[MAXJET];

  TTree *newtreeJet = new TTree("yJet","jets");
  newtreeJet->SetMaxTreeSize(MAXTREESIZE);
  newtreeJet->Branch("nJet",&nJet,"nJet/I");
  newtreeJet->Branch("pt",jetPt,"pt[nJet]/F");
  newtreeJet->Branch("eta",jetEta,"eta[nJet]/F");
  newtreeJet->Branch("phi",jetPhi,"phi[nJet]/F");
  newtreeJet->Branch("dphi",jetDphi,"dphi[nJet]/F");
  if ( isMC )  {
    newtreeJet->Branch("subid",jetSubid,"subid[nJet]/I");
    newtreeJet->Branch("refPt",jetRefPt,"refPt[nJet]/F");
    newtreeJet->Branch("refEta",jetRefEta,"refEta[nJet]/F");
    newtreeJet->Branch("refPhi",jetRefPhi,"refPhi[nJet]/F");
    newtreeJet->Branch("refDphi",jetRefDphi,"refDphi[nJet]/F");
    newtreeJet->Branch("refPartonPt",jetRefPartonPt,"refPartonPt[nJet]/F");
    newtreeJet->Branch("refPartonFlv",jetRefPartonFlv,"refPartonFlv[nJet]/I");
  }

  // 1.2 track tree!
  int nTrk;
  static const int MAXTRK  = 5000;   // This must be  enough.
  float trkPt[MAXTRK];
  float trkEta[MAXTRK];
  float trkPhi[MAXTRK];
  float trkDphi[MAXTRK];
  int   trkPurity[MAXTRK];
  int   trkAlgo[MAXTRK];
  float trkAsJetPt[MAXTRK];  // associated Jet pT
  float trkAsJetEta[MAXTRK];  // associated Jet pT
  float trkAsJetPhi[MAXTRK];  // associated Jet pT
  float trkAsJetDR[MAXTRK];  // associated Jet pT

  TTree *newtreeTrk = new TTree("yTrk","ytrks");
  newtreeTrk->SetMaxTreeSize(MAXTREESIZE);
  newtreeTrk->Branch("nTrk",&nTrk,"nTrk/I");
  newtreeTrk->Branch("pt",trkPt,"pt[nTrk]/F");
  newtreeTrk->Branch("eta",trkEta,"eta[nTrk]/F");
  newtreeTrk->Branch("phi",trkPhi,"phi[nTrk]/F");
  newtreeTrk->Branch("dphi",trkDphi,"dphi[nTrk]/F");
  newtreeTrk->Branch("purity",trkPurity,"purity[nTrk]/I");
  newtreeTrk->Branch("algo",trkAlgo,"algo[nTrk]/I");
  newtreeTrk->Branch("ajPt",trkAsJetPt,"ajPt[nTrk]/F"); // associated jet pt 
  newtreeTrk->Branch("ajEta",trkAsJetEta,"ajEta[nTrk]/F");
  newtreeTrk->Branch("ajPhi",trkAsJetPhi,"ajPhi[nTrk]/F");
  newtreeTrk->Branch("ajDR",trkAsJetDR,"ajDR[nTrk]/F");

  
  const int MAXCh = 10000;
  int nCh;
  int chSube[MAXCh];
  int chChg[MAXCh];
  int chPdg[MAXCh];
  float chPt[MAXCh];
  float chEta[MAXCh];
  float chPhi[MAXCh];
  float chDphi[MAXCh];
  
  TTree *treeChg = new TTree("ch","charged particles stable");
  treeChg->SetMaxTreeSize(MAXTREESIZE);
  treeChg->Branch("nCh",&nCh,"nCh/I");
  treeChg->Branch("pt",chPt,"pt[nCh]/F");
  treeChg->Branch("eta",chEta,"eta[nCh]/F");
  treeChg->Branch("phi",chPhi,"phi[nCh]/F");
  treeChg->Branch("dphi",chDphi,"dphi[nCh]/F");
  treeChg->Branch("pdg",chPdg,"pdg[nCh]/I");
  treeChg->Branch("sube",chSube,"sube[nCh]/I");
  treeChg->Branch("chg",chChg,"chg[nCh]/I");


  
  // 2.0. Background gen charged particle tree
  int nMCh;
  int  mChPdg[MAXMTRK];  // associaated Jet pT
  int  mChChg[MAXMTRK];  // associated Jet pT
  float mChPt[MAXMTRK];
  float mChEta[MAXMTRK];
  float mChPhi[MAXMTRK];
  float mChDphi[MAXMTRK];

  TTree * tmixChg = new TTree("mCh","charged gen particles from minbias events");
  if ( isMC ) { 
    tmixChg->SetMaxTreeSize(MAXTREESIZE);
    tmixChg->Branch("nCh",&nMCh,"nCh/I");
    tmixChg->Branch("pdg",mChPdg,"pdg[nCh]/I");
    tmixChg->Branch("chg",mChChg,"chg[nCh]/I");
    tmixChg->Branch("pt",mChPt,"pt[nCh]/F");
    tmixChg->Branch("eta",mChEta,"eta[nCh]/F");
    tmixChg->Branch("phi",mChPhi,"phi[nCh]/F");
    tmixChg->Branch("dphi",mChDphi,"dphi[nCh]/F");
  }

  // 2.1 Background jet tree 
  int nMJet;
  float mJetPt[MAXMJET];
  float mJetEta[MAXMJET];
  float mJetPhi[MAXMJET];
  float mJetDphi[MAXMJET];
  int   mJetMixingItr[MAXMJET];
  int             mJetSubid[MAXMJET];
  float           mJetRefPt[MAXMJET];
  float           mJetRefEta[MAXMJET];
  float           mJetRefPhi[MAXMJET];
  float           mJetRefPartonPt[MAXMJET];
  int             mJetRefPartonFlv[MAXMJET];

  TTree * tmixJet = new TTree("mJet","Jet from minbias events");
  tmixJet->SetMaxTreeSize(MAXTREESIZE);
  tmixJet->Branch("nJet",&nMJet,"nJet/I");
  tmixJet->Branch("pt",mJetPt,"pt[nJet]/F");
  tmixJet->Branch("eta",mJetEta,"eta[nJet]/F");
  tmixJet->Branch("phi",mJetPhi,"phi[nJet]/F");
  tmixJet->Branch("dphi", mJetDphi, "dphi[nJet]/F");
  tmixJet->Branch("mItr", mJetMixingItr, "mItr[nJet]/I");
  if (isMC) {
    tmixJet->Branch("subid",mJetSubid,"subid[nJet]/F");
    tmixJet->Branch("refPt",mJetRefPt,"refPt[nJet]/F");
    tmixJet->Branch("refEta",mJetRefEta,"refEta[nJet]/F");
    tmixJet->Branch("refPhi",mJetRefPhi,"refPhi[nJet]/F");
    tmixJet->Branch("refPartonPt",mJetRefPartonPt,"refPartonPt[nJet]/F");
    tmixJet->Branch("refPartonFlv",mJetRefPartonFlv,"refPartonFlv[nJet]/F");
  }
  // 2.2 Background track tree
  // 2.2.1. Background tracks from the MinBias events 
  int nMTrk;
  float mTrkPt[MAXMTRK];
  float mTrkEta[MAXMTRK];
  float mTrkPhi[MAXMTRK];
  int   mTrkPurity[MAXMTRK];
  int   mTrkAlgo[MAXMTRK];
  float mTrkDphi[MAXMTRK];
  float mTrkAsJetPt[MAXMTRK];  // associated Jet pT
  float mTrkAsJetEta[MAXMTRK];  // associated Jet pT
  float mTrkAsJetPhi[MAXMTRK];  // associated Jet pT
  float mTrkAsJetDR[MAXMTRK];  // associated Jet pT

  TTree * tmixTrk = new TTree("mTrk","Trk from minbias events");
  tmixTrk->SetMaxTreeSize(MAXTREESIZE);
  tmixTrk->Branch("nTrk",&nMTrk,"nTrk/I");
  tmixTrk->Branch("pt",mTrkPt,"pt[nTrk]/F");
  tmixTrk->Branch("eta",mTrkEta,"eta[nTrk]/F");
  tmixTrk->Branch("phi",mTrkPhi,"phi[nTrk]/F");
  tmixTrk->Branch("purity",mTrkPurity,"purity[nTrk]/I");
  tmixTrk->Branch("algo",mTrkAlgo,"algo[nTrk]/I");
  tmixTrk->Branch("dphi", mTrkDphi, "dphi[nTrk]/F");
  tmixTrk->Branch("ajPt",mTrkAsJetPt,"ajPt[nTrk]/F"); // associated jet pt 
  tmixTrk->Branch("ajEta",mTrkAsJetEta,"ajEta[nTrk]/F");
  tmixTrk->Branch("ajPhi",mTrkAsJetPhi,"ajPhi[nTrk]/F");
  tmixTrk->Branch("ajDR",mTrkAsJetDR,"ajDR[nTrk]/F");

  // 2.2.2. 2nd order Background tracks from the MinBias events 
  int nMmTrk;
  float mmTrkPt[MAXMTRK];
  float mmTrkEta[MAXMTRK];
  float mmTrkPhi[MAXMTRK];
  int   mmTrkPurity[MAXMTRK];
  int   mmTrkAlgo[MAXMTRK];
  float mmTrkDphi[MAXMTRK];
  float mmTrkAsJetPt[MAXMTRK];  // associated Jet pT
  float mmTrkAsJetEta[MAXMTRK];  // associated Jet pT
  float mmTrkAsJetPhi[MAXMTRK];  // associated Jet pT
  float mmTrkAsJetDR[MAXMTRK];  // associated Jet pT

  TTree * tmmixTrk = new TTree("mmTrk","Trk from 2nd order minbias events");
  tmmixTrk->SetMaxTreeSize(MAXTREESIZE);
  tmmixTrk->Branch("nTrk",&nMmTrk,"nTrk/I");
  tmmixTrk->Branch("pt",mmTrkPt,"pt[nTrk]/F");
  tmmixTrk->Branch("eta",mmTrkEta,"eta[nTrk]/F");
  tmmixTrk->Branch("phi",mmTrkPhi,"phi[nTrk]/F");
  tmmixTrk->Branch("purity",mmTrkPurity,"purity[nTrk]/I");
  tmmixTrk->Branch("algo",mmTrkAlgo,"algo[nTrk]/I");
  tmmixTrk->Branch("dphi", mmTrkDphi, "dphi[nTrk]/F");
  tmmixTrk->Branch("ajPt",mmTrkAsJetPt,"ajPt[nTrk]/F"); // associated jet pt 
  tmmixTrk->Branch("ajEta",mmTrkAsJetEta,"ajEta[nTrk]/F");
  tmmixTrk->Branch("ajPhi",mmTrkAsJetPhi,"ajPhi[nTrk]/F");
  tmmixTrk->Branch("ajDR",mmTrkAsJetDR,"ajDR[nTrk]/F");


  // 3. Background jets/tracks from Input Minbias Skim
  // "Imb" means Input MinBias events
  EvtSel evtImb;
  TBranch        *b_evt;

  // 3.0. gen particles
  const int nMaxCh=3000;
  Int_t           nChImb;
  Float_t         chPtImb[nMaxCh];   //Imb[nCh]
  Float_t         chEtaImb[nMaxCh];   //Imb[nCh]
  Float_t         chPhiImb[nMaxCh];   //Imb[nCh]
  Int_t           chPdgImb[nMaxCh];   //Imb[nCh]
  Int_t           chChgImb[nMaxCh];   //Imb[nCh]

  TBranch        *b_nChImb;   //!
  TBranch        *b_chPtImb;   //!
  TBranch        *b_chEtaImb;   //!
  TBranch        *b_chPhiImb;   //!
  TBranch        *b_chPdgImb;   //!
  TBranch        *b_chChgImb;   //!

  // 3.1. Jets
  Int_t           nJetImb;
  Float_t         jetPtImb[MAXJET];
  Float_t         jetEtaImb[MAXJET];
  Float_t         jetPhiImb[MAXJET];
  int             jetSubidImb[MAXJET];
  float           jetRefPtImb[MAXJET];
  float           jetRefEtaImb[MAXJET];
  float           jetRefPhiImb[MAXJET];
  float           jetRefPartonPtImb[MAXJET];
  int             jetRefPartonFlvImb[MAXJET];

  TBranch        *b_nJetImb;
  TBranch        *b_jetPtImb;
  TBranch        *b_jetEtaImb;
  TBranch        *b_jetPhiImb;
  TBranch        *b_jetSubidImb;
  TBranch        *b_jetRefPtImb;
  TBranch        *b_jetRefEtaImb;
  TBranch        *b_jetRefPhiImb;
  TBranch        *b_jetRefPartonPtImb;
  TBranch        *b_jetRefPartonFlvImb;


  
  // 3.2. Tracks  
  // 3.2.1. Tracks in the MB background
  Int_t           nTrkImb;
  Float_t         trkPtImb[MAXTRK];
  Float_t         trkEtaImb[MAXTRK];
  Float_t         trkPhiImb[MAXTRK];
  int             trkPurityImb[MAXTRK];
  int             trkAlgoImb[MAXTRK];
  TBranch        *b_nTrkImb;
  TBranch        *b_trkPtImb;
  TBranch        *b_trkEtaImb;
  TBranch        *b_trkPhiImb;
  TBranch        *b_trkPurityImb;
  TBranch        *b_trkAlgoImb;

  // 3.2.2. Tracks underlyng the jets in MB events
  Int_t           nMTrkImb;
  Float_t         mTrkPtImb[MAXTRK];
  Float_t         mTrkEtaImb[MAXTRK];
  Float_t         mTrkPhiImb[MAXTRK];
  int             mTrkPurityImb[MAXTRK];
  int             mTrkAlgoImb[MAXTRK];
  Float_t         mTrkAsJetPtImb[MAXTRK];
  Float_t         mTrkAsJetEtaImb[MAXTRK];
  Float_t         mTrkAsJetPhiImb[MAXTRK];
  Float_t         mTrkAsJetDRImb[MAXTRK];

  TBranch        *b_nMTrkImb;
  TBranch        *b_mTrkPtImb;
  TBranch        *b_mTrkEtaImb;
  TBranch        *b_mTrkPhiImb;
  TBranch        *b_mTrkPurityImb;
  TBranch        *b_mTrkAlgoImb;
  TBranch        *b_mTrkAsJetPtImb;
  TBranch        *b_mTrkAsJetEtaImb;
  TBranch        *b_mTrkAsJetPhiImb;
  TBranch        *b_mTrkAsJetDRImb;

  
  
  
  int nCentBins =  nCentBinSkim;
  if ((colli==kPADATA)||(colli==kPAMC)) {
    nCentBins = nCentBinSkimPA;
  }

  TChain   *tjmb[200][nVtxBin+1];
  int nMB[200][nVtxBin+1];
  int mbItr[200][nVtxBin+1];
  if ( doMix ) {
    cout <<"  Tree initialization for MinBias mixing" << endl;
    for( int icent = 0 ; icent< nCentBins ; icent++) {
      for( int ivz = 1 ; ivz<=nVtxBin ; ivz++) {

	//	tjmb[icent][ivz] = new TChain(Form("trkAndJets_first_icent%d_ivz%d",icent,ivz));
	tjmb[icent][ivz] = new TChain(Form("trkAndJets_second_icent%d_ivz%d",icent,ivz));
	tjmb[icent][ivz]->Add(MinbiasFname.data());
	tjmb[icent][ivz]->SetBranchAddress("evt", &evtImb,&b_evt);

	//  3.0. charged gen particles
	if ( isMC) {
	  tjmb[icent][ivz]->SetBranchAddress("nCh",   &nChImb,   &b_nChImb);
	  tjmb[icent][ivz]->SetBranchAddress("chPt",  &chPtImb,  &b_chPtImb);
	  tjmb[icent][ivz]->SetBranchAddress("chEta",  &chEtaImb,  &b_chEtaImb);
	  tjmb[icent][ivz]->SetBranchAddress("chPhi",  &chPhiImb,  &b_chPhiImb);
	  tjmb[icent][ivz]->SetBranchAddress("chPdg",  &chPdgImb,  &b_chPdgImb);
	  tjmb[icent][ivz]->SetBranchAddress("chChg",  &chChgImb,  &b_chChgImb);
	}
	//  3.1. jets
	tjmb[icent][ivz]->SetBranchAddress("nJet",   &nJetImb,   &b_nJetImb);
	tjmb[icent][ivz]->SetBranchAddress("jetPt",  &jetPtImb,  &b_jetPtImb);
	tjmb[icent][ivz]->SetBranchAddress("jetEta", &jetEtaImb, &b_jetEtaImb);
	tjmb[icent][ivz]->SetBranchAddress("jetPhi", &jetPhiImb, &b_jetPhiImb);
	if ( isMC )  {
	  tjmb[icent][ivz]->SetBranchAddress("jetSubid", &jetSubidImb, &b_jetSubidImb);
	  tjmb[icent][ivz]->SetBranchAddress("jetRefPt", &jetRefPtImb, &b_jetRefPtImb);
	  tjmb[icent][ivz]->SetBranchAddress("jetRefEta", &jetRefEtaImb, &b_jetRefEtaImb);
	  tjmb[icent][ivz]->SetBranchAddress("jetRefPhi", &jetRefPhiImb, &b_jetRefPhiImb);
	  tjmb[icent][ivz]->SetBranchAddress("jetRefPartonPt", &jetRefPartonPtImb, &b_jetRefPartonPtImb);
	  tjmb[icent][ivz]->SetBranchAddress("jetRefPartonFlv",  &jetRefPartonFlvImb, &b_jetRefPartonFlvImb);
	}
	
	// 3.2.1 Tracks  - 1st kind tracks
	tjmb[icent][ivz]->SetBranchAddress("nTrk",      &nTrkImb,      &b_nTrkImb);
	tjmb[icent][ivz]->SetBranchAddress("trkPt",     trkPtImb,     &b_trkPtImb);
	tjmb[icent][ivz]->SetBranchAddress("trkEta",    trkEtaImb,    &b_trkEtaImb);
	tjmb[icent][ivz]->SetBranchAddress("trkPhi",    trkPhiImb,    &b_trkPhiImb);
	tjmb[icent][ivz]->SetBranchAddress("trkPurity", trkPurityImb, &b_trkPurityImb);
	tjmb[icent][ivz]->SetBranchAddress("trkAlgo",   trkAlgoImb,   &b_trkAlgoImb);
	// 3.2.2 Tracks  - 2nd kind tracks
	tjmb[icent][ivz]->SetBranchAddress("nMTrk",     &nMTrkImb,      &b_nMTrkImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkPt",     mTrkPtImb,     &b_mTrkPtImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkEta",    mTrkEtaImb,    &b_mTrkEtaImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkPhi",    mTrkPhiImb,    &b_mTrkPhiImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkPurity", mTrkPurityImb, &b_mTrkPurityImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkAlgo",   mTrkAlgoImb,   &b_mTrkAlgoImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkAsJetPt",mTrkAsJetPtImb,   &b_mTrkAsJetPtImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkAsJetEta",mTrkAsJetEtaImb,   &b_mTrkAsJetEtaImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkAsJetPhi",mTrkAsJetPhiImb,   &b_mTrkAsJetPhiImb);
	tjmb[icent][ivz]->SetBranchAddress("mTrkAsJetDR",mTrkAsJetDRImb,   &b_mTrkAsJetDRImb);

	nMB[icent][ivz] = tjmb[icent][ivz]->GetEntries();
	cout << "number of evetns in (icent = " << icent << ", ivtxZ = "<< ivz << ")  = " << nMB[icent][ivz] << endl;
	int primeSeed = rand.Integer(37357); // 37357 is an arbitrary prime number set by YS March 17th 2014
	mbItr[icent][ivz] = primeSeed%(nMB[icent][ivz]);
	cout <<" initial itr = " << mbItr[icent][ivz] << endl;
      }
    }
  }
  else
    cout << endl << endl << "  Mixing process is skipped" << endl << endl << endl ;
  
  
  // reweighting factor should go here
  int eTot(0), eSel(0);
  EvtSel evt;
  GammaJet gj;
  Isolation isol;
  TTree *tgj;
  tgj = new TTree("tgj","gamma jet tree");
  tgj->SetMaxTreeSize(MAXTREESIZE);
  tgj->Branch("evt",&evt.run,evtLeaves.Data());
  tgj->Branch("lpho",&gj.photonEt,gammaJetLeaves.Data());
  tgj->Branch("isolation",&isol.cc1, isolLeaves.Data());

  
  // Vertex
  TH1F* hvz = new TH1F("hvz","",nVtxBin,-vtxCutPhotonAna,vtxCutPhotonAna);
  // Event Plane
  TH1F* hEvtPlnBin = new TH1F("hEvtPlnBin", "", nPlnBin, -PI/2., PI/2.);
  // jet algos
  Jets* theJet;
  Jets* genJetTree;

  if ( jetAlgo == "akPu3PF")  {
    theJet = &(c->akPu3PF) ;   cout << "Using akPu3PF Jet Algo" << endl<<endl;
  }
  else if ( jetAlgo == "akVs3PF") {
    theJet = &(c->akVs3PF) ;   cout << "Using ak3PF Jet Algo, Voronoi Subtraction method" << endl<<endl;
  }
  else if ( jetAlgo == "ak3PF")  {
    theJet = &(c->ak3PF) ;   cout << "Using ak3PF Jet Algo" << endl<<endl;
  }
  else if ( jetAlgo == "akVs3Calo") {
    theJet = &(c->akVs3Calo) ;   cout << "Using akVs3Calo Jet Algo, Voronoi Subtraction method" << endl<<endl;
  }
  else {
    cout <<" Jet algo  " << jetAlgo << "is not available" << endl;
    return;
  }
  genJetTree = &(c->akPu3PF);
  
 
  // Trigger bit
  Int_t triggerFlag; 
  TBranch        *b_triggerFlag;
  TChain   *HltChain = new TChain("hltanalysis/HltTree");
  if ( triggerSelection != "" )   {
    HltChain->Add(inputFile_.Data());
    HltChain->SetBranchAddress( triggerSelection.Data(), &triggerFlag, &b_triggerFlag);   // HLT_PAPhoton20_NoCaloIdVL_v1
  }
  
  // Ready to go into the loop!! 
  int nentries = c->GetEntries();
  //  nentries = 500;
  cout << "number of entries = " << nentries << endl;
  for (Long64_t jentry = 0 ; jentry < nentries; jentry++) {
    eTot++;
    if (jentry% 2 == 0)  {
      cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
    }
    
    c->GetEntry(jentry);

    if ( triggerSelection == "" )            triggerFlag = true; 
    else                            HltChain->GetEntry(jentry);   // trigger tree
    
    if ( triggerFlag == false)  {  //      cout << " Trigger is not fired" << endl;
      continue;}
    // Select events with a generated photon in mid-rapidity
    evt.clear();
    evt.run   = c->evt.run;
    evt.evt = c->evt.evt;
    evt.hf4Pos = c->evt.hiHFplusEta4;
    evt.hf4Neg = c->evt.hiHFminusEta4;
    evt.hf4Sum = evt.hf4Pos + evt.hf4Neg;
    evt.cBin = -99;
    evt.pBin   = -99 ;
    if ( (colli==kHIDATA) || (colli==kHIMC) ) {
      evt.cBin = c->evt.hiBin;
      evt.pBin   = hEvtPlnBin->FindBin( c->evt.hiEvtPlanes[theEvtPlNumber] ) ;
    }
    else if ((colli==kPADATA)||(colli==kPAMC))   {
      evt.cBin =  getHfBin(evt.hf4Sum);
      if (  ((evt.cBin) < 0) || (evt.cBin) > 18 )
        cout << " Check the pA centrality..  cbin = " << evt.cBin << endl;
    }

    evt.vz = c->evt.vz;

    int cBinSkim = evt.cBin;
    if ((colli==kHIDATA)||(colli==kHIMC))  
       cBinSkim = getCbinFrom200(evt.cBin);
    
    int vzBin = hvz->FindBin(evt.vz)  ;
    hvz->Fill(evt.vz) ;

    if ( ( (colli==kHIDATA)||(colli==kHIMC)||(colli==kPADATA)||(colli==kPAMC) || (colli==kPPMC) ) && ( c->selectEvent() == 0 ))
      continue;
    //    if ( ( (colli==kPADATA)||(colli==kPPDATA) ) && ( c->skim.pVertexFilterCutGplus ==0 ) ) // No Pile up events But it's not on the forest file anymore (Jul 4 2014)
    //    continue;
    if ( (vzBin<1) || ( vzBin > nVtxBin) )
      continue;


    eSel++;      // OK.  This event is a collisional and no-noise event.
      
    // Reweight for vertex and centrality of MC 
    evt.vtxCentWeight = 1;
    double wVtx=1;
    double wCent=1;
    if (colli ==kHIMC) {
      int vBin =  hWeight_vtx_data_pbpb->FindBin(evt.vz);
      wVtx =   hWeight_vtx_data_pbpb->GetBinContent(vBin) / hWeight_vtx_mc_pbpb->GetBinContent(vBin) ;
      wCent =  hWeight_cent_data_pbpb->GetBinContent(c->evt.hiBin/4 +1) / hWeight_cent_mc_pbpb->GetBinContent(c->evt.hiBin/4 +1) ;
    }
    else if ( colli ==kPPMC) {
      int vBin =  hWeight_vtx_data_pp->FindBin(evt.vz);
      wVtx =   hWeight_vtx_data_pp->GetBinContent(vBin) / hWeight_vtx_mc_pp->GetBinContent(vBin) ;
    }
    else if ( colli ==kPAMC) {
      int vBin =  hWeight_vtx_data_ppb->FindBin(evt.vz);
      wVtx =   hWeight_vtx_data_ppb->GetBinContent(vBin) / hWeight_vtx_mc_ppb->GetBinContent(vBin) ;
    }
    evt.vtxCentWeight = wVtx * wCent;
      
    evt.ptHat = -1;
    evt.ptHatWeight = 1;
    evt.ptHat = c->photon.ptHat; 
    if (colli ==kHIMC) {
      if ( evt.ptHat < 50  )       evt.ptHatWeight = 9008/16237. ; 
      else if ( evt.ptHat < 80  )       evt.ptHatWeight = 3750/85438. ; 
      else   evt.ptHatWeight =  1191/140432. ; 
    }
    else if ( colli == kPPMC) { // pp has only 4 pthat samples 
      if ( evt.ptHat < 50  )       evt.ptHatWeight = 9008/9008. ; 
      else if ( evt.ptHat < 80  )       evt.ptHatWeight = 3750/40109. ; 
      else   evt.ptHatWeight = 1191/66934. ; 
    }
    
    for (int j=0;j< c->photon.nPhotons;j++) {
      
      if (  ( c->photon.pt[j] > preCutPhotonEt ) && ( fabs( c->photon.eta[j] ) < cutphotonEta ) ) {
	newPt[j] =c->photon.pt[j] ;   //  c->getCorrEt(j); Until the new energy scale table is ready (Oct 6th 2014)
      }
      else
	newPt[j] = c->photon.pt[j] - 10000;

      if ( (c->isSpike(j)) || (c->photon.hadronicOverEm[j]>0.2) ) //||  (c->photon.isEle[j]) )
	newPt[j] = newPt[j] - 20000;
      if (c->photon.seedTime[j] ==0 )   // clustering bug
	newPt[j] = newPt[j] - 30000;

      corrPt[j] = newPt[j];
    }

    TMath::Sort(c->photon.nPhotons, newPt, order);

  
    // Select the leading photon
    gj.clear();
    int leadingIndex=-1;

    for (int j=0;j<c->photon.nPhotons;j++) {
      if ( c->photon.pt[j]  < preCutPhotonEt ) continue;
      if ( fabs(c->photon.eta[j]) > cutphotonEta ) continue;
      if (c->isSpike(j)) continue;
      //   if (!(c->isLoosePhoton(j))) continue;
      if (c->photon.hadronicOverEm[j]>0.1) continue;
      if ((c->photon.rawEnergy[j]/c->photon.energy[j])<0.5) continue;


      
      // sort using corrected photon pt
      float theCorrPt= corrPt[j];
      if ( theCorrPt > gj.photonEt) {
	gj.photonEt = theCorrPt;
	leadingIndex = j;
      }
    }
    //    if ( (gj.photonEt < cutphotonPt) )    <== This cut ruins the ptHat weighting factor
    //  continue;


    /// Save leading photons
    if (leadingIndex!=-1) {
      gj.photonRawEt=c->photon.pt[leadingIndex];
      gj.photonEta=c->photon.eta[leadingIndex];
      gj.photonPhi=c->photon.phi[leadingIndex];
      gj.hovere=c->photon.hadronicOverEm[leadingIndex];
      gj.r9=c->photon.r9[leadingIndex];
      gj.sigmaIetaIeta=c->photon.sigmaIetaIeta[leadingIndex];
      gj.sumIsol = (c->photon.cr4[leadingIndex]+c->photon.cc4[leadingIndex]+c->photon.ct4PtCut20[leadingIndex]) / 0.9;
      gj.genIso   = c->photon.genCalIsoDR04[leadingIndex];
      gj.genPhotonEt = c->photon.genMatchedPt[leadingIndex];
      gj.genMomId = c->photon.genMomId[leadingIndex];
      isol.Set(c,leadingIndex);
    }
    else {
      gj.clear();
    }
    if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
      gj.photonEta = - gj.photonEta;
    }

    
    ///////////////////// 1.0. Gen Particle Tree ///////////////////////////////////
    nCh = 0;
    if ( isMC) { 
      for (int it=0; it< c->genparticle.mult ; it++) {
	if ( c->genparticle.pt[it] < cuttrkPtSkim )   continue;
	if (  fabs(c->genparticle.eta[it]) > cuttrkEtaSkim ) continue;
	if ( c->genparticle.chg[it] == 0 ) continue;
	//	if ( c->genparticle.sube[it] != 0 ) continue;
	
	chPdg[nCh]  = c->genparticle.pdg[it];
	chChg[nCh]  = c->genparticle.chg[it];
	chSube[nCh]  = c->genparticle.sube[it];
	chPt[nCh]  = c->genparticle.pt[it];
	chEta[nCh]  = c->genparticle.eta[it];
	chPhi[nCh]  = c->genparticle.phi[it];
	chDphi[nCh] = getAbsDphi( chPhi[nCh], gj.photonPhi) ;
	nCh++;}
    }  
    
    ///////////////////// Jet tree ///////////////////////////////////
    nJet = 0 ;
    int jetEntries = 0;
    if (useGenJetColl )    jetEntries = theJet->ngen;
    else                   jetEntries = theJet->nref;
    
    int nSmear = 1;
    if(smearingCentBin != -1)
      nSmear = 100;
    
    for (int ij=0; ij< jetEntries ; ij++) {
      if ( gj.photonEt < 0 )    continue ;    // If there is no photon in this event
      for(int iSmear =0; iSmear < nSmear; iSmear++){
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
      	
	// Smear phi
	Double_t newPhi = jetPhi[nJet] ;
	if( smearingCentBin != -1 )
	{
	  Double_t phiSmear  = getPhiSmear(smearingCentBin, jetPt[nJet]); 				    
	  newPhi  =  jetPhi[nJet] +  rand.Gaus(0, phiSmear);
	}
	jetPhi[nJet] = normalAngle(newPhi);
	
	
	// smear the jet pT
	Double_t smeared = jetPt[nJet];
	if( smearingCentBin != -1 )     {
	  Double_t smearSigma = getPtSmear( smearingCentBin,  jetPt[nJet]) ; 
	  smeared = jetPt[nJet] * rand.Gaus(1, smearSigma);
	}
	// then multiply residual correction factor
	float resCorrection =1. ;
	float l2l3Corr =1. ;
	if  (doJetResCorrection)  {
	  
	  if   ((colli==kHIDATA)||(colli==kHIMC)) 
	    resCorrection = getHiResCorr( evt,  jetPt[nJet]); 
	  else if ((colli==kPPDATA)||(colli==kPPMC))
	    resCorrection = getPpResCorr( jetPt[nJet]); 
	  
	  // L2L3 correction!
	  l2l3Corr = getL2L3Corr(evt, colli,  jetPt[nJet], jetEta[nJet], fgaus);
	}
	jetPt[nJet] = smeared * l2l3Corr /resCorrection;


	// reflect eta! 
	if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
	  jetEta[nJet] = -jetEta[nJet];
	}
      
	// Jet kinematic cuts //////////////////////////////////////////////
	if ( jetPt[nJet] < cutjetPtSkim) 
	  continue;
	if ( fabs( jetEta[nJet] ) > cutjetEtaSkim )
	  continue;
	if ( getDR( jetEta[nJet], jetPhi[nJet], gj.photonEta, gj.photonPhi) < 0.5 )
	  continue;
	///////////////////////////////////////////////////////////////////////////////

	if (jetPt[nJet] >0)
	  jetDphi[nJet] = getAbsDphi( jetPhi[nJet], gj.photonPhi) ;
	else
	  jetDphi[nJet] = -1;


	if (  useGenJetColl )   {
	  jetSubid[nJet] = -9999;
	  jetRefPt[nJet] = -9999;
	  jetRefEta[nJet] = -9999;
	  jetRefPhi[nJet] = -9999;
	  jetRefPt[nJet] =  -9999;
	  jetRefPartonPt[nJet] = -9999;
	  jetRefPartonFlv[nJet] = -9999;
	}
	else {
	  jetRefPt[nJet] = theJet->refpt[ij];
	  jetRefEta[nJet] = theJet->refeta[ij];
	  jetRefPhi[nJet] = theJet->refphi[ij];
	  if (jetRefPt[nJet] >0)
	    jetRefDphi[nJet] = getAbsDphi( jetRefPhi[nJet] , gj.photonPhi) ;
	  else
	    jetRefDphi[nJet] = -1;
	  jetRefPartonPt[nJet] = theJet->refparton_pt[ij];
	  jetRefPartonFlv[nJet] = theJet->refparton_flavor[ij];

	  jetSubid[nJet] = -9999; 
	  if (jetRefPt[nJet] >0)  {   // Find the collisional subid of this gen jet!!
	    for ( int igen=0; igen < genJetTree->ngen ; igen++) { 
	      if ( jetRefPt[nJet] == genJetTree->genpt[igen] )
		jetSubid[nJet] = genJetTree->gensubid[igen] ;
	    }
	    
	    if ( jetSubid[nJet] == -9999 ) // No genJet matched! 
	      cout << " WARNING!  This reco jet was not matched to anyone in the gen jet collection!!! " << endl;
	    
	  }
	}
	
      }
      
      nJet++ ;
    }
  
    ///////////////////// 1.2. Track tree ///////////////////////////////////
    nTrk = 0 ;
    for (int it=0; it< c->track.nTrk ; it++) {
      if ( c->track.trkPt[it] < cuttrkPtSkim )   continue;
      if (  fabs(c->track.trkEta[it]) > cuttrkEtaSkim ) continue;
      if ( c->selectTrack(it)  == false) continue;
      trkPt[nTrk]  = c->track.trkPt[it];
      trkEta[nTrk] = c->track.trkEta[it];
      trkPhi[nTrk] = c->track.trkPhi[it];
      trkPurity[nTrk] = c->track.highPurity[it];
      trkAlgo[nTrk] = c->track.trkAlgo[it];
      trkDphi[nTrk] = getAbsDphi( trkPhi[nTrk], gj.photonPhi) ;
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
      
      
      nTrk++;}
    
    
    nMJet = 0;
    nMTrk = 0;
    nMmTrk = 0;
    nMCh   =0 ;
    bool noSuchEvent = false;
    int iMix=0;
    int loopCounter=0;
    if ( (!doMix) || ( gj.photonEt < 0) )
      iMix = nMixing1+1;   // Mixing step will
    
    
    while (iMix<nMixing1)  {
      loopCounter++;
      if ( loopCounter > nMB[cBinSkim][vzBin]+1) { 
	// If the photon-jet event are not matched with any of
	// MinBias Skim events (Centrality/vertex/event plane etc.etc.)
	iMix = 999999 ;
	noSuchEvent = true;
	cout << " no such event!!";
	cout << "cBinSkim = " << cBinSkim << ",  vzBin = " << vzBin;
	cout << ",  pBin = " << evt.pBin << endl;
	continue;
      }
    
      mbItr[cBinSkim][vzBin] = mbItr[cBinSkim][vzBin] + 1;
      if ( mbItr[cBinSkim][vzBin] == nMB[cBinSkim][vzBin] )
	mbItr[cBinSkim][vzBin] =  mbItr[cBinSkim][vzBin] - nMB[cBinSkim][vzBin];
      
      tjmb[cBinSkim][vzBin]->GetEntry(mbItr[cBinSkim][vzBin]);
      // ok found the event!! ///////////
      loopCounter =0;  // Re-initiate loopCounter
      
      // Jet mixing
      for (int it = 0 ; it < nJetImb ; it++) {
	// Smear phi
	Double_t newPhi = jetPhiImb[it];
	if( smearingCentBin != -1 )     {
	  Double_t phiSmear  = getPhiSmear(smearingCentBin, jetPtImb[nJet]);
          newPhi  =  jetPhiImb[it] +   rand.Gaus(0, phiSmear);
	}
	jetPhiImb[nJet] = normalAngle(newPhi);
	
	// jet pT residual correction and smearing 
	Double_t smeared = jetPtImb[it];
	if( smearingCentBin != -1 )     {
          Double_t smearSigma = getPtSmear( smearingCentBin,  jetPtImb[it]);
          smeared = jetPt[nJet] * rand.Gaus(1, smearSigma);
        }
	float resCorrection =1. ;	float l2l3Corr =1 ;
	if  (doJetResCorrection)   {
	  // Correction from MC closure
	  if   ((colli==kHIDATA)||(colli==kHIMC))
            resCorrection = getHiResCorr( evt,  jetPtImb[it]);
          else if ((colli==kPPDATA)||(colli==kPPMC))
            resCorrection = getPpResCorr( jetPtImb[it]);
	  // L2L3 correction for pp 
          l2l3Corr = getL2L3Corr( evt, colli,  jetPtImb[it], jetEtaImb[it], fgaus);
	}
	float smearedCorrected  = smeared *l2l3Corr / resCorrection; 
	/////////////////////////////////////////////////////////////
	
	
	if ( smearedCorrected < cutjetPtSkim )  // double cutjetPtSkim = 15; Oct 19th
	  continue;
	if ( fabs( jetEtaImb[it] ) > cutjetEtaSkim )   // double cutjetEtaSkim = 3.0; Oct 19th
	  continue;
	if ( getDR( jetEtaImb[it], jetPhiImb[it], gj.photonEta, gj.photonPhi) < 0.5 )  // This cut added for consistency ; Oct 19th
	  continue;

	mJetPt[nMJet]    = smearedCorrected;
	mJetEta[nMJet]   = jetEtaImb[it];
	mJetPhi[nMJet]   = jetPhiImb[it];
	mJetMixingItr[nMJet]   = iMix;
	if  ( mJetPt[nMJet]>0 )  mJetDphi[nMJet] =  getAbsDphi(mJetPhi[nMJet], gj.photonPhi) ;
	else            	 mJetDphi[nMJet] = -1;
	if ( isMC) {
	  mJetSubid[nMJet] = jetSubidImb[it];
	  mJetRefPt[nMJet] = jetRefPtImb[it];
	  mJetRefEta[nMJet] = jetRefEtaImb[it];
	  mJetRefPhi[nMJet] = jetRefPhiImb[it];
	  mJetRefPartonPt[nMJet] = jetRefPartonPtImb[it];
	  mJetRefPartonFlv[nMJet] = jetRefPartonFlvImb[it];
	}
	nMJet++; // < == Important!
      }

          
      // 3.2.0 Gen Particles
      if ( isMC) {
	for (int it = 0 ; it < nChImb ; it++) {
	  if (  chPtImb[it] < cuttrkPtSkim )   continue;
	  if (  fabs(chEtaImb[it]) > cuttrkEtaSkim ) continue;
	  if (  chChgImb[it] == 0 ) continue;
	  mChPdg[nMCh]  = chPdgImb[it];
	  mChChg[nMCh]  = chChgImb[it];
	  mChPt[nMCh]  = chPtImb[it];
	  mChEta[nMCh]  = chEtaImb[it];
	  mChPhi[nMCh]  = chPhiImb[it];
	  mChDphi[nMCh] = getAbsDphi( chPhiImb[it], gj.photonPhi) ;
	  nMCh++;}
      }
          
      // 3.2.1 Tracks  - 1st kind tracks
      for (int it = 0 ; it < nTrkImb ; it++) {
	if ( trkPtImb[it] < cuttrkPtSkim )   continue;
	if (  fabs(trkEtaImb[it]) > cuttrkEtaSkim ) continue;
	mTrkPt[nMTrk]  = trkPtImb[it];
	mTrkEta[nMTrk] = trkEtaImb[it];
	mTrkPhi[nMTrk] = trkPhiImb[it];
	mTrkDphi[nMTrk] =  getAbsDphi(mTrkPhi[nMTrk], gj.photonPhi) ;
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
	
	nMTrk++;}

      
      // 3.2.1 Tracks  - 2nd kind tracks   
      if ( fillSecondMtrk ) 
	nMTrkImb = 0;

      for (int it = 0 ; it < nMTrkImb ; it++) {
	if ( mTrkPtImb[it] < cuttrkPtSkim )   continue;
	  if (  fabs(mTrkEtaImb[it]) > cuttrkEtaSkim ) continue;
	  mmTrkPt[nMmTrk]  = mTrkPtImb[it];
	  mmTrkEta[nMmTrk] = mTrkEtaImb[it];
	  mmTrkPhi[nMmTrk] = mTrkPhiImb[it];
	  mmTrkDphi[nMmTrk] =  getAbsDphi(mmTrkPhi[nMmTrk], gj.photonPhi) ;
	  mmTrkPurity[nMmTrk] = mTrkPurityImb[it];
	  mmTrkAlgo[nMmTrk] = mTrkAlgoImb[it];
	  int assocJetId = matchedJetFinder( theJet, mmTrkEta[nMmTrk], mmTrkPhi[nMmTrk]);
	  if ( assocJetId < 0 )  {
	    mmTrkAsJetPt[nMmTrk] = -1;
	    mmTrkAsJetEta[nMmTrk] = -1;
	    mmTrkAsJetPhi[nMmTrk] = -1;
	    mmTrkAsJetDR[nMmTrk] = 100;
	  }
	  else {
	    mmTrkAsJetPt[nMmTrk] = theJet->jtpt[assocJetId];
	    mmTrkAsJetEta[nMmTrk] = theJet->jteta[assocJetId];
	    mmTrkAsJetPhi[nMmTrk] = theJet->jtphi[assocJetId];
	    mmTrkAsJetDR[nMmTrk] =getDR( mmTrkEta[nMmTrk], mmTrkPhi[nMmTrk], theJet->jteta[assocJetId], theJet->jtphi[assocJetId]) ;
	  }
	  nMmTrk++;}
      
      iMix++;}
    
    
    tgj->Fill();
    newtreeJet->Fill();
    newtreeTrk->Fill();
    tmixJet->Fill();
    tmixTrk->Fill();
    if ( fillSecondMtrk ) 
      tmmixTrk->Fill();

    newtreePhoton->Fill();
    if ( isMC )   {
      treeChg->Fill();
      tmixChg->Fill();
    }
  }  
  
  
  newfile_data->Write();
  //   newfile_data->Close();   // <<=== If there is close() function. writing stucks in the middle of looping.. I don't know why!!
  cout << " Done! "<< endl;
  cout << "    " << eSel<<" out of total "<<eTot<<" events were analyzed."<<endl;
}











double getPhiSmear( int smearingCentBin, float jetPt){
  if ( smearingCentBin < 0 ) {
    cout << " forest2yskim_jetSkim_forest::getPhiSmear   WARNING!!!!  negative centrality bin" << endl;
    return -1;
  }
  
  return TMath::Sqrt((cphi_pbpb[smearingCentBin]*cphi_pbpb[smearingCentBin] - cphi_pp*cphi_pp)
		     + (sphi_pbpb[smearingCentBin]*sphi_pbpb[smearingCentBin] - sphi_pp*sphi_pp)/jetPt
		     + (nphi_pbpb[smearingCentBin]*nphi_pbpb[smearingCentBin] - nphi_pp*nphi_pp)/(jetPt*jetPt));
}

double getPtSmear( int smearingCentBin, float jetPt) {
  TMath::Sqrt((c_pbpb[smearingCentBin]*c_pbpb[smearingCentBin] - c_pp*c_pp)
	      + (s_pbpb[smearingCentBin]*s_pbpb[smearingCentBin] - s_pp*s_pp)/jetPt
	      + (n_pbpb[smearingCentBin]*n_pbpb[smearingCentBin] - n_pp*n_pp)/(jetPt*jetPt));
}



float getHiResCorr( EvtSel evt, float jetPt) {
  float resCorrection=1;
  if ( evt.cBin  < 60 )   // central 30%
    resCorrection  =  1.04503 -1.6122  /(sqrt(jetPt)) + 9.27212 / (jetPt);  //1.04503    -1.6122    9.27212
  else                  // peripheral
    resCorrection  =  1.00596 -0.653191/(sqrt(jetPt)) + 4.35373 / (jetPt);  //1.00596     -0.653191  4.35373
  return resCorrection;
}

float getPpResCorr( float jetPt) {
  return 0.993609 + 0.158418/(sqrt(jetPt)) + 0.335479/jetPt ;
}

float  getL2L3Corr(EvtSel evt,  sampleType colli, float jetPt, float jetEta, TF1* fgaus) {
  float l2l3Corr = 1;
  
  if ( colli == kPPDATA)   {
    l2l3Corr = c_etapp->GetBinContent(c_etapp->FindBin(jetEta)) * fptpp->Eval( jetPt);
  }
  else if ( colli == kPADATA)   {
    if ( evt.run > 211256 )
      l2l3Corr = c_etapA->GetBinContent(c_etapA->FindBin(jetEta))  * fptpA->Eval( jetPt);
    else
      l2l3Corr = c_etaAp->GetBinContent(c_etaAp->FindBin(jetEta))  * fptAp->Eval( jetPt);
  }
  else if ( colli == kPAMC)
    l2l3Corr = 1 + (fsmear_pA->Eval( jetPt )) * fgaus->GetRandom()  ;
  else {
    cout << " Error in forest2yskim_jetSkim_forestV3::getL2L3Corr!!  l2l3 correction is not defined in this collisions" << endl;
    l2l3Corr = 1;
  }
  
  return l2l3Corr;
}

void drawL2L3Correciton() {
  TCanvas* c11 = new TCanvas("c11","",1200,400);   // valiation of smearing factors
  c11->Divide(3,1);
  c11->cd(1);
  c_etapp->Draw();
  c11->cd(2);
  c_etaAp->Draw();
  c11->cd(3);
  c_etapA->Draw();
  c11->SaveAs("f1.gif");
}

