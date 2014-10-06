#include "CutAndBinCollection2012.h"
#include <TRandom3.h>
#include <time.h>

void drawPAFigures(int xNorm = 1) { 
  
  const int nCentBinPa = 3;
  double centBinPa[nCentBinPa+1] = {0,20,30,100};
  double AvePtBin[nCentBinPa+1] = { 10, 25, 65, 9999};
  
  TH1D* hxjg[7][1000]; // [species][centrality] 
  TH1D* hxjgNorm[7][1000];
  TH1D* hdphi[7][1000];
  
  valPair meanX[7][1000];
  valPair rjg[7][1000];
  TGraphAsymmErrors* gr[7];
  TGraphAsymmErrors* gx[7];
  for (int i=0;i<7;i++) {  
    gr[i] = new TGraphAsymmErrors();
    gr[i]->SetName(Form("gr_%d",i));
    gx[i] = new TGraphAsymmErrors();
    gx[i]->SetName(Form("gx_%d",i));
  }
  TFile* fff[7][1000];
  
  for ( int icent = 1 ; icent<=nCentBinPa  ; icent++) { 
    int iSpecies =kPADATA ; 
    
    int iglb = icent + iSpecies*50.;
    
    
    fff[iSpecies][iglb] = new TFile("ffFiles/photonTrackCorr_ppbMC_output_photonPtThr60_to_9999_jetPtThr30_20130813.root");
    hdphi[iSpecies][iglb] = (TH1D*)fff[iSpecies][iglb]->Get(Form("jetDphi_icent%d_final",icent));
    
    hxjg[iSpecies][iglb] = (TH1D*)fff[iSpecies][iglb]->Get(Form("jetXjg_icent%d_final",icent));
    hxjgNorm[iSpecies][iglb] = (TH1D*)hxjg[iSpecies][iglb]->Clone(Form("norm_%s",hxjg[iSpecies][iglb]->GetName()) );
    
    meanX[iSpecies][iglb].val = hxjg[iSpecies][iglb]->GetMean();
    meanX[iSpecies][iglb].err = hxjg[iSpecies][iglb]->GetMeanError();
    gx[iSpecies]->SetPoint     (icent, AvePtBin[icent-1], meanX[iSpecies][iglb].val);
    gx[iSpecies]->SetPointError(icent, AvePtBin[icent-1] - centBinPa[icent-1], centBinPa[icent] - AvePtBin[icent-1], meanX[iSpecies][iglb].err, meanX[iSpecies][iglb].err);
    
    double rError;
    rjg[iSpecies][iglb].val = hxjg[iSpecies][iglb]->IntegralAndError(1, hxjg[iSpecies][iglb]->GetNbinsX(), rError, "width");
    rjg[iSpecies][iglb].err = rError;
    gr[iSpecies]->SetPoint     (icent, AvePtBin[icent-1], rjg[iSpecies][iglb].val);
    gr[iSpecies]->SetPointError(icent,  AvePtBin[icent-1] - centBinPa[icent-1], centBinPa[icent] - AvePtBin[icent-1], rjg[iSpecies][iglb].err, rjg[iSpecies][iglb].err);
    
    
  }
  
  TCanvas* c1 = new TCanvas("c1","",900,700);
  makeMultiPanelCanvas(c1,3,2,0.0,0.0,0.2,0.15,0.02);
  for ( int icent = 1 ; icent<=nCentBinPa  ; icent++) {
    
    
    c1->cd(icent);
    
    int iSpecies = kPADATA;
    int iglb = icent + iSpecies*50.;
    handsomeTH1(hxjgNorm[iSpecies][iglb],1);
    hxjgNorm[iSpecies][iglb]->Scale(0.1);
    if (xNorm==1) scaleInt(hxjgNorm[iSpecies][iglb]);
    hxjgNorm[iSpecies][iglb]->SetXTitle("x_{J,#gamma}");
    hxjgNorm[iSpecies][iglb]->SetYTitle("#frac{dN}{dx} #frac{1}{N}");
    if (xNorm==1)     hxjgNorm[iSpecies][iglb]->SetAxisRange(-0.015,0.3,"Y");
    else              hxjgNorm[iSpecies][iglb]->SetAxisRange(-0.06,1.5,"Y");
    hxjgNorm[iSpecies][iglb]->SetMarkerStyle(20);
    hxjgNorm[iSpecies][iglb]->Draw();

    
    float dx1;    
    ((icent==1)||(icent==4))? dx1=0.15 : dx1=0 ;
    
    if ( icent == nCentBinPa ) 
      drawText(Form("E_{T}^{HF|#eta|>4} > %dGeV, ", (int)centBinPa[icent]), 0.12+dx1+0.25,0.889118,1,15);
    else 
      drawText(Form("%dGeV < E_{T}^{HF|#eta|>4} < %dGeV, ", (int)centBinPa[icent-1], (int)centBinPa[icent]), 0.12+dx1,0.889118,1,15);

      
    TLegend *l1 = new TLegend(0.6365615,0.6445304,0.9577623,0.846736,NULL,"brNDC");
    easyLeg(l1,"5.02TeV");
    l1->AddEntry(hxjgNorm[kPADATA][icent + kPADATA*50],"pPb ","p");
    if ( icent==1 )   l1->Draw();
    onSun(0,0,2,0);
    jumSun(1,0,1,1);
  }

  TCanvas* c11 = new TCanvas("c11","",500,500);
  for ( int icent = 1 ; icent<=nCentBinPa  ; icent++) {
    int iSpecies = kPADATA;
    int iglb = icent + iSpecies*50.;
    handsomeTH1(hxjgNorm[iSpecies][iglb],ycolor[icent]);

    if ( icent == 1 ) hxjgNorm[iSpecies][iglb]->Draw();
    else     hxjgNorm[iSpecies][iglb]->Draw("same");
    onSun(0,0,2,0);
  }


  TCanvas* cDphi = new TCanvas("cDphi","",900,700);
  makeMultiPanelCanvas(cDphi,3,2,0.0,0.0,0.2,0.15,0.02);
  for ( int icent = 1 ; icent<=nCentBinPa  ; icent++) {
    cDphi->cd(icent);
    
    int iSpecies = kPADATA;
    int iglb = icent + iSpecies*50.;
    handsomeTH1(hdphi[iSpecies][iglb],1);
    hdphi[iSpecies][iglb]->SetXTitle("#Delta#phi_{J,#gamma}");
    hdphi[iSpecies][iglb]->Scale(1./hdphi[iSpecies][iglb]->Integral("width"));
    hdphi[iSpecies][iglb]->SetAxisRange(0.01,30,"Y");
    //    hdphi[iSpecies][iglb]->SetMarkerStyle(24);
    handsomeTH1(hdphi[iSpecies][iglb],ycolor[icent]);
    hdphi[iSpecies][iglb]->Draw();

    gPad->SetLogy();
    float dx1;    
    ((icent==1)||(icent==4))? dx1=0.15 : dx1=0 ;
    

    if ( icent == nCentBinPa )
      drawText(Form("E_{T}^{HF|#eta|>4} > %dGeV", (int)centBinPa[icent]), 0.32+dx1+0.25,0.889118,1,15);
    else
      drawText(Form("%dGeV < E_{T}^{HF|#eta|>4} < %dGeV", (int)centBinPa[icent-1], (int)centBinPa[icent]), 0.32+dx1,0.889118,1,15);
    

  }

  
  TCanvas* c2 = new TCanvas("c2","",800,400);
  c2->Divide(2,1);
  c2->cd(2);
  TH1D* hr = new TH1D("hr",";E_{T}^{HF|#eta|>4} (GeV);r_{J,#gamma}",200,0,100);
  handsomeTH1(hr,1);
  hr->SetAxisRange(0.0,1.1,"Y");
  //  hr->SetNdivisions(5);
  hr->Draw();
  handsomeTGraph(gr[kPADATA],1);
  gr[kPADATA]->SetMarkerStyle(20);
  gr[kPADATA]->Draw("same p");

  jumSun(0,1,100,1);


  
  c2->cd(1);
  TH1D* hx = (TH1D*)hr->Clone("hx");
  handsomeTH1(hx,1);
  hx->SetYTitle("<x_{J,#gamma}>");
  hx->SetAxisRange(0.6,1.05,"Y");
  //  hx->SetNdivisions(5);
  hx->Draw();
  handsomeTGraph(gx[kPADATA],1);
  gx[kPADATA]->SetMarkerStyle(20);
  gx[kPADATA]->Draw("same p");
  jumSun(0,1,100,1);

  TLegend *l1=   new TLegend(0.5226071,0.6300415,0.9597724,0.8332584,NULL,"brNDC");
  easyLeg(l1,"5.02TeV");
  l1->AddEntry(gr[kPADATA],"pPb ","p");
  l1->Draw();


}

