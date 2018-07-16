//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3




TH1F cut(const TH1F&  h, Double_t m){

  int nbins = h.GetNbinsX();
  float xmin = h.GetXaxis()->GetXmin();
  float xmax = h.GeXaxis()->GetXmax();
  TH1F* h1_new = new TH1F("h_cut", "", nbins, xmin, xmax);
  
  int i;
  for (i=0;i<500; i++){
    if (i < m* 0.8*500 || i> m*3*500){
      h1_new->SetBinContent(i,0);
    } else {
      h1_new->SetBinContent(i,h->GetBinContent(i));
    }
  }  return h1_new;
}

void plotWF_cut(const char * filename){
  
  
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
  
  // TH1F *hl_cut =new TH1F("hr_cut","histos_cut",500,0.0,1);
  
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
  Double_t mip = fit_r ->GetParameter(1);
  
  TH1F hr_cut = cut(hr_amp,mip);
  hr_cut->Draw("same");
  
  
  hr_amp->Draw("");
  wf_c->cd(2)->SetLogy();
  hl_amp->Fit("f_l","RV"); 
  hl_amp->Draw("");


}


