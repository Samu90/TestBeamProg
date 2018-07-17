//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_tamp(const char * filename){
  
  
  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");


  
  Float_t amp_max[54], time[54];
  int k;
  bool debug=false;
  Double_t max=0;
  
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",500,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",500,0.0,1);
  TF1 *fit_r = new TF1("f_r","landau",0.10,2);
  TF1 *fit_l = new TF1("f_l","landau",0.10,2);
  
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  // digiTree->SetBranchAddress("LED30",&LED30);
  //digiTree->SetBranchAddress("LED50",&LED50);
  
  for(k=0; k<digiTree->GetEntries(); k++){
    digiTree->GetEntry(k);
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
  }/*chiudo for k*/
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%3000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for k*/
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  cout<< max << endl;
  
  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");
  
  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",300,0,0.5,500,7.5,15); 
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",300,0,0.5,500,7.5,15); 
 
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*fit_l->GetParameter(1) < amp_max[3]/max && amp_max[3]/max < 3*fit_l->GetParameter(1))
      { 
	h2_l->Fill(amp_max[3]/max,time[3]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[4]-time[0]);
	if(debug)cout<< amp_max[4]/max << "   "<< time[4]-time[0]<<endl;
      }
    
  }//chiudo for k
  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Divide(2,1);
  wf_c->cd(1);
  h2_l->Draw("COLZ");
  wf_c->cd(2);
  h2_r->Draw("COLZ");
  
  
  
}


