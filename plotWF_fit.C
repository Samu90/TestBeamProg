
void plotWF_fit(const char * filename){
  
  
  TFile *  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[54];
  int k;
  Double_t max=0;
  
 
 

  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",500,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",500,0.0,1);
  TF1 *fit_r = new TF1("f_r","landau",0.13,2);
  TF1 *fit_l = new TF1("f_l","landau",0.13,2);
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
  }//chiudo for k
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%3000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }//chiudo for k
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Clear();
  cout<< max << endl;

  gStyle->SetOptFit();
  wf_c->Divide(2,1);
  wf_c->cd(1)->SetLogy();
  hr_amp->Fit("f_r","RV"); 
  hr_amp->Draw("");
  wf_c->cd(2)->SetLogy();
  hl_amp->Fit("f_l","RV"); 
  hl_amp->Draw("");

  /*
  TString histoname = "";
  histoname.Append("histo");
  histoname.Append(filename);

  
 TFile* f2 = new TFile(histoname.Data(),"RECREATE");
  f2->cd();
  hr_amp->Write();
  hl_amp->Write();
  f2->Close();
  // wf_c->cd(2)->SetLogy();
  //hr_amp->Draw();
  //hl_amp->Draw("same");
  */
}

