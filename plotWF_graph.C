
using namespace std;
void plotWF_graph(const char * filename){
  
  
  TFile *  file= TFile::Open(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[54];
  int k;
  Double_t max=0;
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Divide(2,1);
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",500,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",500,0.0,1);

  digiTree->SetBranchAddress("amp_max",&amp_max);
  
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
  }//chiudo for k
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%100==0) cout<<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }//chiudo for k   
  
  hl_amp->SetLineColor(kRed);
  
  cout<< max << endl;
  wf_c->cd(1)->SetLogy();
  hr_amp->DrawNormalized();
  wf_c->cd(2)->SetLogy();
  hl_amp->DrawNormalized();
  
  
  // wf_c->cd(2)->SetLogy();
  //hr_amp->Draw();
  //hl_amp->Draw("same");
  
}

