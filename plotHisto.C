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
  TLegend* l1=new TLegend(0.1,0.7,0.48,0.9);
  l1->SetHeader("Configuration","C");
  l1->AddEntry(hl1,"Conf1.3");
  l1->AddEntry(hl2,"Conf4.1");
  hl1->SetLineColor(4);
  hl1->DrawNormalized();
  hl2->SetLineColor(2);
  hl2->DrawNormalized("same");
  l1->Draw();
  
  c1->cd(2)->SetLogy();
  TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
  l2->SetHeader("Configuration","C");
  l2->AddEntry(hr1,"Conf1.3");
  l2->AddEntry(hr2,"Conf4.1");
  hr1->SetLineColor(4);
  hr1->DrawNormalized();
  hr2->SetLineColor(2);
  hr2->DrawNormalized("same");
  l2->Draw();
}
