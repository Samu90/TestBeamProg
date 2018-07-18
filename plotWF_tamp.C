//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_tamp(const char * filename){
  
  
  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");


  
  Float_t amp_max[54], time[54];
  int k,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax;
  bool debug=false;
  Double_t max=0,tmax=0;
  rxmin=0;
  rxmax=0.5;
  rymin_l=5;
  rymax_l=21;
  rymin_r=7;
  rymax_r=17;
  tymin=6;
  tymax=17;
  const Int_t  nbin=100;

  Float_t x_r[nbin],y_r[nbin], x_l[nbin],y_l[nbin],rmsy_l[nbin],rmsy_r[nbin];
  Float_t xt[nbin],yt[nbin],rmsyt[nbin];
  

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
    if(time[3]-time[4]>tmax && time[3]-time[4]<10) {tmax = time[3]-time[4];}

  }/*chiudo for */

  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for */
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  cout << tmax <<endl;
  cout<< max << endl;
  
  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");
  
  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbin,rxmin,rxmax,nbin,rymin_l,rymax_l); 
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbin,rxmin,rxmax,nbin,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbin,-0.4,0.8,nbin,tymin,tymax); 
 
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)) && (time[3]-time[4])<7 && time[3]-time[4]>0)
      { 
	h2_l->Fill(amp_max[3]/max,time[3]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[4]-time[0]);
	h2_t->Fill((time[3]-time[4])/tmax,(time[3]+time[4])/2-time[0]);

	
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4]-time[0] <<endl;
      }
    
  }//chiudo for k
  
 
  for(k=0;k<nbin;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;
  
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);
    
    
    maxbin_l=histotemp_l->GetMaximumBin();
    maxbin_r=histotemp_r->GetMaximumBin();
    maxbin_t=histotemp_t->GetMaximumBin();
    
    xt[k]=-0.4+(Float_t)(0.8-(-0.4))/nbin*k;
    yt[k]=tymin+(Float_t)(tymax-tymin)/nbin*maxbin_t;
    rmsyt[k]=histotemp_t->GetRMS();
       
    x_l[k]=(rxmax-rxmin)/nbin*k;
    y_l[k]=rymin_l+(rymax_l-rymin_l)/nbin*maxbin_l;
    rmsy_l[k]=histotemp_l->GetRMS();

    x_r[k]=(rxmax-rxmin)/nbin*k;
    y_r[k]=rymin_r+(rymax_r-rymin_r)/nbin*maxbin_r;
    rmsy_r[k]=histotemp_r->GetRMS();

    
    
    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    
    if(k%20==0) cout << k << " / " << nbin << endl;
  }//chiudo for k

  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,550);
  TGraphErrors* graph_r=new TGraphErrors(nbin-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbin-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbin-1,xt,yt,0,rmsyt);
  TF1* hyp_r = new TF1("hyp_r","[0] + [1]/(x**[3]-[2])",rxmin,rxmax);
  TF1* hyp_l = new TF1("hyp_l","[0] + [1]/(x**[3]-[2])",rxmin,rxmax);
  TF1* hyp_t = new TF1("hyp_t","[0] + [1]*x",-0.5,1);
  hyp_r->SetParameter(0,7);
  hyp_l->SetParameter(0,7);
  hyp_r->SetParameter(1,10);
  hyp_l->SetParameter(0,10);
  
  wf_c->Divide(3,1);
  
  wf_c->cd(1);
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_r","R");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("P");
  
  
  wf_c->cd(2);
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_l","R");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);  
  graph_r->Draw("P");

  wf_c->cd(3);
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","R");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);  
  graph_t->Draw("P");
  
  
  
  
  
}


