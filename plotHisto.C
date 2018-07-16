void plotHisto(const char * filename1,const char * filename2){
  
  
  TFile*  f1= TFile::Open(filename1);
  TFile*  f2= TFile::Open(filename2);
  
  TH1F* hr1=(TH1F*)f1->Get("hr_amp");
  TH1F* hl1=(TH1F*)f1->Get("hl_amp");
  
  TH1F* hr2=(TH1F*)f2->Get("hr_amp");
  TH1F* hl2=(TH1F*)f2->Get("hl_amp");

  TCanvas* c1= new TCanvas("mycanvas","",1000,500);
  c1->Divide(2,1);
  
  c1->cd(1)->SetLogy();
  
  hl1->SetLineColor(1);
  hl1->DrawNormalized();

  hl2->SetLineColor(2);
  hl2->DrawNormalized("same");
  
  
  c1->cd(2)->SetLogy();
  hr1->SetLineColor(3);
  hr1->DrawNormalized();
  hr2->SetLineColor(4);
  hr2->DrawNormalized("same");
  
}
