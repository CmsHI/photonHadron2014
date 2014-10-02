#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>

class histFunction2
{ 
public:
  histFunction2(TH1D *h, TH1D *h2);
  ~histFunction2(){ 
    delete histBck;
    delete histSig;
  };
  
  //      int GetBinNumber(double c);
  Double_t evaluate(Double_t *x, Double_t *par);
  TH1D *histSig;
      TH1D *histBck;
      Double_t lowEdge;
      Double_t highEdge;
      Double_t nbin;
};

histFunction2::histFunction2(TH1D *h,TH1D *h2)
{
   histSig=(TH1D*)h->Clone();
   histBck=(TH1D*)h2->Clone();
   
   nbin=h->GetNbinsX();
   lowEdge=h->GetBinLowEdge(1);
   highEdge=h->GetBinLowEdge(nbin+1);
   
   histSig->SetName("hSig");
   histSig->Scale(1./histSig->Integral(1,histSig->GetNbinsX()+1));
   histBck->SetName("hBck");
   histBck->Scale(1./histBck->Integral(1,histBck->GetNbinsX()+1));
}

Double_t histFunction2::evaluate(Double_t *x, Double_t *par) {
     
   Double_t xx = x[0];
   int binNum=histSig->FindBin(xx);  //
   
   return par[0]*(histSig->GetBinContent(binNum)*par[1]+histBck->GetBinContent(binNum)*(1-par[1]));
}

/* example: 
   histFunction2 *myfun = new histFunction2(hist);
   TF1 *test = new TF1("landau",myfun,&histFunction2::evaluate,0,0.7,3);
   
*/



/*                                                                                                        
  int histFunction2::GetBinNumber(double c)                                                               
  {                                                                                                       
                                                                                                          
  double delta = c - lowEdge;                                                                             
  double n = delta/(highEdge-lowEdge)*nbin;                                                               
  return (int)(n+1);                                                                                      
  }                                                                                                       
*/
