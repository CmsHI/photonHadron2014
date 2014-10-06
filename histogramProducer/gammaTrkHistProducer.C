///////////////////////////////////////////////////////////////////                                
// gammaJetHistProducer.C                                        //
// Creator : Yongsun Kim (Korea Univ.), kimy@cern.ch             //
// Function : Make the correlation function between Photon and   //
//            tracks.  Signal and Background are calculated      //
//            separately, and then background subtraction is     //
//            done.  Output histogrmas are saved in output file  //
///////////////////////////////////////////////////////////////////         
#include <TStyle.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TCut.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include "../CutAndBinCollection2012.h"
#include "corrFunctionJetTrk.h"


GjSpectra* nullGj;

void gammaTrkSingle(     GjSpectra* gjSpec_=nullGj,
			 multiTreeUtil* tObj[3]=0,
			 corrFunctionTrk* corr_=0,
			 float purity=0,
			 sampleType collision_=kHIDATA,
			 TString var_ = "",
			 TCut cut_ ="",
			 TString theWeight="",
			 TCut phoCandCut_="",
			 TCut phoDecayCut_="",
			 TH1D* hTemplate_ =0,
			 TString outName="");

fitResult getPurity(TString fname="", sampleType collision=kHIDATA, TCut candEvtCut = "", TCut sbEvtCut="", TString ccanvasName="",float photonPtThr=60, float photonPtThrUp=9999);


void gammaTrkHistProducer(sampleType collision = kPPDATA, float photonPtThr=20, float photonPtThrUp=9999, int icent =7){
  TH1::SetDefaultSumw2();
  
  TString stringSampleType = getSampleName(collision); "";
  
  TDatime* date = new TDatime();
  TString  outName=  Form("photonTrackCorr_%s_output_photonPtThr%d_to_%d_%d.root",stringSampleType.Data(),(int)photonPtThr, (int)photonPtThrUp,  date->GetDate());
  delete date;
  
 
  int lowerCent(0),  upperCent(0); 
  TCut centCut  = "";
  if ( (collision ==kHIDATA) || (collision ==kHIMC) ) {
     lowerCent = ((icent/100)%100) *2 ;
     upperCent =  (icent%100)*2 -1 ; 
     centCut = Form("cBin >= %d && cBin<= %d",lowerCent,upperCent);
  }
  else if (  (collision ==kPPDATA) || (collision==kPPMC)  ){  // if it's pp 
    centCut = "(1==1)";
    //    icent = 7;   // for pp, centrality is set as the smearing 
  }
  else { // pPb
    centCut = Form( "hf4Sum > %f && hf4Sum <= %f", (float)centBinPa[icent-1], (float)centBinPa[icent]);
  }
  
  cout <<" centrality : " << centCut.GetTitle() << endl;
      
  ///////// Photon cut //////////////////////////////////////////////////////////////////////////////
  
  cout <<" photon pt >" << photonPtThr << " GeV" << endl;
  TCut ptPhoCut  = Form("photonEt>%.1f && photonEt<%.1f", (float)photonPtThr, (float)photonPtThrUp  );
  TCut caloIso;
  
  if ( (collision==kPPMC) || (collision==kPPDATA) ) 
    caloIso = "(ecalIso < 4.2  &&  hcalIso < 2.2  &&  trackIso < 2) && hovere<0.1";
  else if ( (collision==kHIMC) || (collision==kHIDATA) )
    caloIso = "(sumIso<5) && hovere<0.1";
  else {
    caloIso = "ecalIso < 4.2  &&  hcalIso < 2.2  &&  trackIso < 2 && hovere<0.1";
  }
  
  TCut sbIso   = "(sumIso>5) && (sumIso<20) && hovere<0.1";
  //  if ( (collision==kPPMC) || (collision==kPPDATA) || (collision==kPAMC) || (collision==kPADATA)  )
  //  sbIso   = "ecalIso < 4.2  &&  hcalIso < 2.2 && trackIso > 2 && trackIso < 5 && hovere<0.1";

  TCut basicPhoCut = centCut && ptPhoCut && caloIso ;
  TCut sbPhoCut    = centCut && ptPhoCut && sbIso   ;
  TCut evtSeltCut = basicPhoCut;
  TCut sbSeltCut  = sbPhoCut;

  TCut phoCandCut   = "sigmaIetaIeta<0.010";
  TCut phoDecayCut  = "(sigmaIetaIeta>0.011) && (sigmaIetaIeta<0.017)";
  if  ( ( collision == kHIMC ) || (collision == kPPMC) || (collision == kPAMC))  
    phoCandCut = phoCandCut && "genIso<5 && abs(genMomId)<=22";
  
  
  TString fname = "";
  if ( collision == kHIDATA)      fname = fnameHIDATA; //
  else if ( collision == kPADATA) fname = fnamePADATA;
  else if ( collision == kPPDATA) {
    if  ( icent == 7 ) fname = fnamePPDATA;
    else if ( icent == 10010 ) fname = fnamePPDATA0010;
    else if ( icent == 11030 ) fname = fnamePPDATA1030;
    else if ( icent == 13050 ) fname = fnamePPDATA3050;
    else if ( icent == 15099 ) fname = fnamePPDATA5099;
    else if ( icent == 10030 ) fname = fnamePPDATA0030;
    else if ( icent == 13099 ) fname = fnamePPDATA30100;
  }  
  else fname = "";
  
  multiTreeUtil* tgj = new multiTreeUtil();
  tgj->addFile(fname,  "tgj",  evtSeltCut,  1);
  tgj->AddFriend("yTrk");
 
  float purity(0);
  
  TString canvasName = Form("gifs/purity_%s_output_icent%d_photonPtThr%d-%d", stringSampleType.Data(),  (int)icent, (int)photonPtThr, (int)photonPtThrUp);
  
  if ( collision == kPPDATA) {  
    purity = 0.85;
  }
  else {
    fitResult fitr = getPurity(fname, collision, evtSeltCut, sbSeltCut, canvasName, photonPtThr, photonPtThrUp);
    purity = fitr.purity010;
  }
  GjSpectra* gSpec = new GjSpectra();
  gSpec->init(Form("icent%d",(int)icent) );
  tgj->Draw2(gSpec->hPtPhoCand,  "photonEt", phoCandCut, "");
  tgj->Draw2(gSpec->hPtPhoDecay, "photonEt", phoDecayCut, "");
  
  // Obtain background subtracted spectra
  
  float candInt = gSpec->hPtPhoCand->Integral();
  float decayInt = gSpec->hPtPhoDecay->Integral();
  gSpec->hPtPhoSig->Reset();
  gSpec->hPtPhoSig->Add(gSpec->hPtPhoCand);
  gSpec->hPtPhoSig->Add(gSpec->hPtPhoDecay, -(1. - purity) * candInt / decayInt);
  gSpec->hPtPhoSig->Scale(1./purity ) ;
  
  TFile outf = TFile(Form("ffFiles/%s",outName.Data()),"update");
  gSpec->hPtPhoCand->Write();
  gSpec->hPtPhoDecay->Write();
  gSpec->hPtPhoSig->Write();
  outf.Close();
  
  
  // Objects
  multiTreeUtil* tObj[3];
  tObj[kTrkRaw] = new multiTreeUtil();
  tObj[kTrkBkg] = new multiTreeUtil();
  
  tObj[kTrkRaw]->addFile(fname,  "yTrk",  evtSeltCut,  1);
  tObj[kTrkBkg]->addFile(fname,  "mTrk",  evtSeltCut,  1);
  
  tObj[kTrkRaw]->AddFriend("tgj");
  tObj[kTrkBkg]->AddFriend("tgj");
  TCut trkCut     =  Form("abs(eta)<%f && pt>%f", (float)cuttrkEta, (float)cuttrkPt );
  
  TString varTrkDphi         = Form("dphi");

  corrFunctionTrk* cTrkDphi = new corrFunctionTrk();
  TH1D* hTrkDphi = new TH1D(Form("dphi_icent%d",icent),";#Delta#phi_{Jet,#gamma} ;dN/d#Delta#phi",20,0,3.141592);
  gammaTrkSingle( gSpec,  tObj, cTrkDphi,  purity, 
		  collision, varTrkDphi, trkCut, "pt",
		  phoCandCut, phoDecayCut,  hTrkDphi, outName);

  TH1D* hTrkDphi1to2GeV = new TH1D(Form("dphi_icent%d_pt1to2GeV",icent),";#Delta#phi_{Jet,#gamma} ;dN/d#Delta#phi",20,0,3.141592);
  corrFunctionTrk* cTrkDphi1to2GeV = new corrFunctionTrk();
  gammaTrkSingle( gSpec,  tObj, cTrkDphi1to2GeV,  purity, 
		  collision, varTrkDphi, trkCut && "pt > 1 && pt <=2", "pt",
		  phoCandCut, phoDecayCut,  hTrkDphi, outName);

  TH1D* hTrkDphi2to4GeV = new TH1D(Form("dphi_icent%d_pt2to4GeV",icent),";#Delta#phi_{Jet,#gamma} ;dN/d#Delta#phi",20,0,3.141592);
  corrFunctionTrk* cTrkDphi2to4GeV = new corrFunctionTrk();
  gammaTrkSingle( gSpec,  tObj, cTrkDphi2to4GeV,  purity, 
		  collision, varTrkDphi, trkCut && "pt > 2 && pt <=4", "pt",
		  phoCandCut, phoDecayCut,  hTrkDphi, outName);

  TH1D* hTrkDphi4to8GeV = new TH1D(Form("dphi_icent%d_pt4to8GeV",icent),";#Delta#phi_{Jet,#gamma} ;dN/d#Delta#phi",20,0,3.141592);
  corrFunctionTrk* cTrkDphi4to8GeV = new corrFunctionTrk();
  gammaTrkSingle( gSpec,  tObj, cTrkDphi4to8GeV,  purity, 
		  collision, varTrkDphi, trkCut && "pt > 4 && pt <=8", "pt",
		  phoCandCut, phoDecayCut,  hTrkDphi, outName);

  TH1D* hTrkDphi8andHighGeV = new TH1D(Form("dphi_icent%d_pt8andHighGeV",icent),";#Delta#phi_{Jet,#gamma} ;dN/d#Delta#phi",20,0,3.141592);
  corrFunctionTrk* cTrkDphi8andHighGeV = new corrFunctionTrk();
  gammaTrkSingle( gSpec,  tObj, cTrkDphi8andHighGeV,  purity, 
		  collision, varTrkDphi, trkCut && "pt > 8", "pt",
		  phoCandCut, phoDecayCut,  hTrkDphi, outName);


}


void gammaTrkSingle(     GjSpectra* gSpec,  multiTreeUtil* tObj[3],   corrFunctionTrk* corr, 
			 float purity,       sampleType collision, 	     TString var,     
			 TCut cut, 		TString theWeight,	 TCut phoCandCut,   TCut phoDecayCut,
			 TH1D* hTemplate,  TString outfName)
{
  TH1::SetDefaultSumw2();
  
  TCanvas* c1 = new TCanvas(Form("canvas_%s",hTemplate->GetName()),"",800,650);
  makeMultiPanelCanvas(c1,2,2,0.0,0.0,0.2,0.15,0.02);
  c1->cd(1);
  corr->init(gSpec, collision, purity, hTemplate);
  cout << "Filling raw jets" << endl;
  tObj[kTrkRaw]->Draw2(corr->Func[kPhoCand][kTrkRaw],  var,  phoCandCut  && cut, theWeight);
  if ( tObj[kTrkBkg]!=0)   {
  cout << "Filling mixed jets" << endl;
    tObj[kTrkBkg]->Draw2(corr->Func[kPhoCand][kTrkBkg],  var,  phoCandCut  && cut, theWeight);
  }
  
  if ( (collision==kHIDATA) || (collision==kPPDATA) || (collision==kPADATA) ) {
    cout << "Filling decay photon - jet correlation" << endl;
    tObj[kTrkRaw]->Draw2(corr->Func[kPhoDecay][kTrkRaw], var, phoDecayCut && cut,  theWeight);
    if (tObj[kTrkBkg]!=0)   {
      cout << "Filling decay photon - mixed jet correlation" << endl;
      tObj[kTrkBkg]->Draw2(corr->Func[kPhoDecay][kTrkBkg], var, phoDecayCut && cut,  theWeight);
    }
  }
  TH1ScaleByWidth( corr->Func[kPhoCand][kTrkRaw]); 
  TH1ScaleByWidth( corr->Func[kPhoCand][kTrkBkg]);
  TH1ScaleByWidth( corr->Func[kPhoDecay][kTrkRaw]);
  TH1ScaleByWidth( corr->Func[kPhoDecay][kTrkBkg]);
  
  corr->calCorr();
  c1->cd(1); 
  handsomeTH1(corr->Func[kPhoCand][kTrkRaw],1);
  handsomeTH1(corr->Func[kPhoCand][kTrkBkg],1);
  handsomeTH1(corr->Func[kPhoCand][kTrkSig],2);
  corr->Func[kPhoCand][kTrkRaw]->Draw();
  corr->Func[kPhoCand][kTrkBkg]->Draw("same hist");
  corr->Func[kPhoCand][kTrkSig]->Draw("same");
  gPad->SetLogy();

  c1->cd(2); 
  handsomeTH1(corr->Func[kPhoDecay][kTrkRaw],1);
  handsomeTH1(corr->Func[kPhoDecay][kTrkBkg],1);
  handsomeTH1(corr->Func[kPhoDecay][kTrkSig],4);
  corr->Func[kPhoDecay][kTrkRaw]->Draw();
  corr->Func[kPhoDecay][kTrkBkg]->Draw("same hist");
  corr->Func[kPhoDecay][kTrkSig]->Draw("same");
  gPad->SetLogy();
  
  c1->cd(3);
  corr->Func[kPhoCand][kTrkSig]->Draw();
  corr->Func[kPhoDecay][kTrkSig]->Draw("same");

  c1->cd(4);
  handsomeTH1(corr->Func[kPhoSig][kTrkSig],1);
  corr->Func[kPhoSig][kTrkSig]->Draw();
  c1->SaveAs(Form("gifs/%s_%s.gif",outfName.Data(),c1->GetName()) );
  
  TFile outf = TFile(Form("ffFiles/%s",outfName.Data()),"update");
  corr->Func[kPhoCand][kTrkRaw]->Write();
  corr->Func[kPhoCand][kTrkBkg]->Write();
  corr->Func[kPhoCand][kTrkSig]->Write();
  corr->Func[kPhoDecay][kTrkRaw]->Write();
  corr->Func[kPhoDecay][kTrkBkg]->Write();
  corr->Func[kPhoDecay][kTrkSig]->Write();
  corr->Func[kPhoSig][kTrkRaw]->Write();
  corr->Func[kPhoSig][kTrkBkg]->Write();
  corr->Func[kPhoSig][kTrkSig]->Write();
  
  outf.Close();  
}



fitResult getPurity(TString fname, sampleType collision, TCut evtSeltCut, TCut sbEvtCut, TString canvasName, float photonPtThr, float photonPtThrUp) {
  double purity(0);
  
  multiTreeUtil* tgj = new multiTreeUtil();
  multiTreeUtil* tgjMC = new multiTreeUtil();
  cout << " Calculating Purity....." << endl;           
  tgj->addFile(fname,  "tgj",  "",  1);         //  tgj->AddFriend("yTrk");
  if   (collision==kPPDATA) { 
    tgjMC->addFile(fnamePPMC_AllQcdPho30to50,    "tgj", "", wPPMC_AllQcdPho30to50 ); 
    tgjMC->addFile(fnamePPMC_AllQcdPho50to80,    "tgj", "", wPPMC_AllQcdPho50to80 ); 
    tgjMC->addFile(fnamePPMC_AllQcdPho80to120,   "tgj", "", wPPMC_AllQcdPho80to120); 
    tgjMC->addFile(fnamePPMC_AllQcdPho120to9999, "tgj", "", wPPMC_AllQcdPho120to9999); 
  }
  else if (collision==kPADATA)     {
    tgjMC->addFile(fnamePAMC_AllQcdPho30to50,     "tgj","", wPAMC_AllQcdPho30to50);
    tgjMC->addFile(fnamePAMC_AllQcdPho50to80,     "tgj","", wPAMC_AllQcdPho50to80);
    tgjMC->addFile(fnamePAMC_AllQcdPho80to120,    "tgj","", wPAMC_AllQcdPho80to120);
    tgjMC->addFile(fnamePAMC_AllQcdPho120to9999,  "tgj","", wPAMC_AllQcdPho120to9999);
  }
  else if (collision==kHIDATA)  {
     /*    tgjMC->addFile(fnameHIMC_AllQcdPho30to50,     "tgj", "",wHIMC_AllQcdPho30to50);
	   tgjMC->addFile(fnameHIMC_AllQcdPho50to80,     "tgj", "",wHIMC_AllQcdPho50to80);
	   tgjMC->addFile(fnameHIMC_AllQcdPho80to9999,   "tgj", "",wHIMC_AllQcdPho80to9999);       */
     tgjMC->addFile(fnameHIMC_AllQcdPho30to9999,   "tgj", "",wHIMC_AllQcdPho80to9999);
  }   
  else { 
    cout << " Error: getPurity.  check the type of the collision!  " << endl;
    fitResult fitr0;
    return fitr0;
  }
  
  TH1D* hCand = new TH1D("cand","",25,0,0.025);
  TH1D* hBkg = (TH1D*)hCand->Clone("bkg");  TH1D* hSig = (TH1D*)hCand->Clone("sig");
  
  tgj->Draw2(   hCand, "sigmaIetaIeta", evtSeltCut , "");
  tgj->Draw2(   hBkg, "sigmaIetaIeta", sbEvtCut , "");
  tgjMC->Draw2( hSig, "sigmaIetaIeta", evtSeltCut && "genIso<5 && abs(genMomId)<=22", "");

  handsomeTH1(hCand,1);
  handsomeTH1(hSig,2);
  handsomeTH1(hBkg,4);
  
  hCand->Draw();
  hSig->Draw("same");
  hBkg->Draw("same hist");

  TCanvas* cPurity = new TCanvas("cpurity","",500,500);
  fitResult  fitr = doFit ( hSig, hBkg, hCand, 0.005, 0.025);
  drawText(Form("Purity : %.2f", (float)fitr.purity010), 0.5680963,0.429118);
  drawText(Form("p_{T}^{#gamma}: %d-%d GeV", (int)photonPtThr, (int)photonPtThrUp), 0.568,0.529118);
  cPurity->SaveAs( Form("%s.pdf",canvasName.Data() ) );
  //  gPad->SetLogy();   cPurity->SaveAs( Form("%s_logScale.pdf",canvasName.Data() ) );
   
  TCanvas* c1 = new TCanvas("c1","",100,100);
  
  
  delete tgj;
  delete tgjMC;
  delete hSig;
  delete hBkg;
  delete hCand;
  
  return fitr;
  
  
}

