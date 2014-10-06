#include "../CutAndBinCollection2012.h"



//combinations
const int kPhoCand=1;
const int kPhoDecay=2;
const int kPhoSig = 0;

const int kTrkRaw=1;
const int kTrkBkg=2;
const int kTrkSig=0;

const int kJetRaw=1;
const int kJetBkg=2;
const int kJetSig=0;

//kPure is shared with photon




class GjSpectra
{
 public:
  
  void init(TString hName);
  
  ~GjSpectra(){};
  
  TH1D *hPtPhoCand;
  TH1D *hPtPhoDecay;
  TH1D *hPtPhoSig;
  
};


void GjSpectra::init(TString hName)  { 
  hPtPhoCand = new TH1D(Form("hPhoPtCand_%s",hName.Data()),";p_{T} (GeV);Entries", 60,0,300);
  hPtPhoDecay = (TH1D*)hPtPhoCand->Clone(Form("hPhoPtDecay_%s",hName.Data()));
  hPtPhoSig     =  (TH1D*)hPtPhoCand->Clone(Form("hPhoPtSig_%s",hName.Data()));
}


class corrFunctionTrk
{
 public:
  void init(GjSpectra* gjSpec_, int fsampleType, double fPurity, TH1D* hTemplate);
  void setNorm();
  void calCorr();
  
  ~corrFunctionTrk(){};
  
  //      int GetBinNumber(double c);                                                                                   
  TH1D* Func[3][3];
  
 private:
  float nPhoCand;
  float nPhoDecay;
  float nPhoSig;
  int sampleType;
  float photonPurity;
   
};

void corrFunctionTrk::init(GjSpectra* gjSpec_, int fsampleType, double fPurity, TH1D* hTemplate) {
  TH1::SetDefaultSumw2();
  
  sampleType = fsampleType;
  photonPurity = fPurity;
  nPhoCand =0;
  nPhoDecay =0;
  nPhoSig =0;
  
  if (  ( sampleType == kHIMC )   ||( sampleType == kPPMC )  ||( sampleType == kPAMC )   ) 
    photonPurity = 1 ;
    
  TString hName1;
  if ( sampleType == kHIMC ) 
    hName1 = "hi_mc";
  else if ( sampleType == kHIDATA) {
    hName1 = "hi_data";
  }
  else if ( sampleType == kPPDATA) {
    hName1 = "pp_data";
  }
  else if ( sampleType == kPPMC) {
    hName1 = "pp_mc";
  }
  else if ( sampleType == kPADATA) {
    hName1 = "pa_data";
  }
  else if ( sampleType == kPAMC) {
    hName1 = "pa_mc";
  }
  
  // count the number of candidates 
  nPhoCand  = gjSpec_->hPtPhoCand->Integral(1, gjSpec_->hPtPhoCand->GetNbinsX());
  nPhoDecay = gjSpec_->hPtPhoDecay->Integral(1, gjSpec_->hPtPhoDecay->GetNbinsX());
  setNorm();
  
  hTemplate->Reset();
  Func[kPhoCand][kTrkRaw] = (TH1D*)hTemplate->Clone(Form("%s_phoCand_rawTrk",hTemplate->GetName()));
  Func[kPhoCand][kTrkBkg] = (TH1D*)hTemplate->Clone(Form("%s_phoCand_bkgTrk",hTemplate->GetName()));
  Func[kPhoCand][kTrkSig] = (TH1D*)hTemplate->Clone(Form("%s_phoCand_sigTrk",hTemplate->GetName()));
  
  Func[kPhoDecay][kTrkRaw] = (TH1D*)hTemplate->Clone(Form("%s_phoDecay_rawTrk",hTemplate->GetName()));
  Func[kPhoDecay][kTrkBkg] = (TH1D*)hTemplate->Clone(Form("%s_phoDecay_bkgTrk",hTemplate->GetName()));
  Func[kPhoDecay][kTrkSig] = (TH1D*)hTemplate->Clone(Form("%s_phoDecay_sigTrk",hTemplate->GetName()));
  
  Func[kPhoSig][kTrkRaw] = (TH1D*)hTemplate->Clone(Form("%s_phoSig_rawTrk",hTemplate->GetName()));
  Func[kPhoSig][kTrkBkg] = (TH1D*)hTemplate->Clone(Form("%s_phoSig_bkgTrk",hTemplate->GetName()));
  Func[kPhoSig][kTrkSig] = (TH1D*)hTemplate->Clone(Form("%s_final",hTemplate->GetName()));

}



void corrFunctionTrk::setNorm() {

  cout << " ============= Photon counting ============================ " << endl;
    
  if (nPhoCand==0)  cout << "ERROR :  strange!!! n_photon = 0 " << endl;
  
  if  (sampleType == kHIMC){
    nPhoSig = nPhoCand;
    cout << " --------PbPb MC------------------------" << endl;
  }
  if  (sampleType == kPAMC){
    nPhoSig = nPhoCand;
    cout << " -------- pPb MC------------------------" << endl;
  }
  if  (sampleType == kPPMC){
    nPhoSig = nPhoCand;
    cout << " --------  pp MC------------------------" << endl;
  }
  else if ( sampleType == kHIDATA) {
    nPhoSig = nPhoCand * photonPurity;
    cout << " --------PbPb DATA------------------------" << endl;
  }
  else if ( sampleType == kPPDATA) {
    nPhoSig = nPhoCand * photonPurity;
    cout << " --------  pp DATA------------------------" << endl;
  }
  else if ( sampleType == kPADATA) {
    nPhoSig = nPhoCand * photonPurity;
    cout << " -------- pA DATA------------------------" << endl;
  }

  cout << "purity = " << photonPurity << endl;
  cout << "Pure photons  = " << nPhoSig  << "  (out of " << nPhoCand << " candidates)" <<endl;


  cout << " Number of Sideband Photons = " << nPhoDecay << endl;
  cout << " ========================================================== " << endl;
  
}


void corrFunctionTrk::calCorr() {

  // First Divide by repetition number of data mixing                                                
  float nMixing = nMixing1;
  Func[kPhoCand][kTrkBkg]->Scale(1./nMixing);
  Func[kPhoDecay][kTrkBkg]->Scale(1./nMixing);
  
  Func[kPhoCand][kTrkRaw]->Scale(1./nPhoCand);
  Func[kPhoCand][kTrkBkg]->Scale(1./nPhoCand);
  
  if ( (sampleType == kPPDATA) || (sampleType == kPADATA) || (sampleType == kHIDATA) ) {
    Func[kPhoDecay][kTrkRaw]->Scale(1./nPhoDecay);
    Func[kPhoDecay][kTrkBkg]->Scale(1./nPhoDecay);
  }
  else { 
    Func[kPhoDecay][kTrkRaw]->Reset();
    Func[kPhoDecay][kTrkBkg]->Reset();
    cout << " it's MC.. so decay photon histograms are reset to zero" << endl;
  }
  // If it is pp 
  if ( (sampleType == kPPDATA) ||  (sampleType == kPPMC) ) {
    cout << " it's pp data/mc.  background jet histogram will be set as zero.  No background subtraction for jets" << endl;
    Func[kPhoDecay][kTrkBkg]->Reset();
    Func[kPhoCand][kTrkBkg] ->Reset();
  }
  
  //  
  // pure tracks                                                                                                           
  Func[kPhoCand][kTrkSig]->Reset();
  Func[kPhoCand][kTrkSig]->Add(Func[kPhoCand][kTrkRaw],1);
  Func[kPhoCand][kTrkSig]->Add(Func[kPhoCand][kTrkBkg],-1);
  
  Func[kPhoDecay][kTrkSig]->Reset();
  Func[kPhoDecay][kTrkSig]->Add(Func[kPhoDecay][kTrkRaw],1);
  Func[kPhoDecay][kTrkSig]->Add(Func[kPhoDecay][kTrkBkg],-1);

  Func[kPhoSig][kTrkSig]->Reset();
  Func[kPhoSig][kTrkSig]->Add(Func[kPhoCand][kTrkSig],1);
  Func[kPhoSig][kTrkSig]->Add(Func[kPhoDecay][kTrkSig], -1 * (1-photonPurity) );
  Func[kPhoSig][kTrkSig]->Scale( 1./photonPurity);

  // other tracks
  Func[kPhoSig][kTrkRaw]->Reset();
  Func[kPhoSig][kTrkRaw]->Add(Func[kPhoCand][kTrkRaw],1);
  Func[kPhoSig][kTrkRaw]->Add(Func[kPhoDecay][kTrkRaw], -1 * (1-photonPurity) );
  Func[kPhoSig][kTrkRaw]->Scale( 1./photonPurity);
  
  Func[kPhoSig][kTrkBkg]->Reset();
  Func[kPhoSig][kTrkBkg]->Add(Func[kPhoCand][kTrkBkg],1);
  Func[kPhoSig][kTrkBkg]->Add(Func[kPhoDecay][kTrkBkg], -1 * (1-photonPurity) );
  Func[kPhoSig][kTrkBkg]->Scale( 1./photonPurity);
  
  
}

