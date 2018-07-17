//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_tamp(const char * filename){
  
  
  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");


  
  Float_t amp_max[54], time[54];
  int k,nbin,maxbin_l,maxbin_r;
  Float_t rxmin,rxmax,rymin,rymax;
  bool debug=true;
  Double_t max=0;
  rxmin=0;
  rxmax=0.5;
  rymin=7.5;
  rymax=18;
  nbin=800;

  Float_t x_r[nbin],y_r[nbin], x_l[nbin],y_l[nbin];

  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbin,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbin,0.0,1);
  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);
  
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  // digiTree->SetBranchAddress("LED30",&LED30);
  //digiTree->SetBranchAddress("LED50",&LED50);
  
  for(k=0; k<digiTree->GetEntries(); k++){
    digiTree->GetEntry(k);
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
  }/*chiudo for k*/
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for k*/
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  cout<< max << endl;
  
  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");
  
  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbin,rxmin,rxmax,nbin,rymin,rymax); 
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbin,rxmin,rxmax,nbin,rymin,rymax); 
 
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      { 
	h2_l->Fill(amp_max[3]/max,time[3]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[4]-time[0]);
	
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4]-time[0] <<endl;
      }
    
  }//chiudo for k
  
  TH1D* histotemp_l;
  TH1D* histotemp_r;

  for(k=0;k<nbin;k++){
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    
    maxbin_l=histotemp_l->GetMaximumBin();
    maxbin_r=histotemp_r->GetMaximumBin();
    
    x_r[k]=(rxmax-rxmin)/nbin*k;
    y_r[k]=rymin+(rymax-rymin)/nbin*maxbin_r;
    
    x_l[k]=(rxmax-rxmin)/nbin*k;
    y_l[k]=rymin+(rymax-rymin)/nbin*maxbin_l;
    
    delete histotemp_l;
    delete histotemp_r;
    
    if(k%20==0) cout << k << " / " << nbin << endl;
  }//chiudo for k

  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  TGraphErrors* graph_r=new TGraphErrors(nbin-1,x_r,y_r,0,0);
  TGraphErrors* graph_l=new TGraphErrors(nbin-1,x_l,y_l,0,0);
  
  wf_c->Divide(2,1);
  
  wf_c->cd(1);
  h2_l->Draw("COLZ");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("P");
  
  
  wf_c->cd(2);
  h2_r->Draw("COLZ");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);  
  graph_r->Draw("P");
  
  
  
  
  
}


