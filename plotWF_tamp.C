//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_tamp(const char * filename){
  
  
  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");


  
  Float_t amp_max[54], time[54];
  int k;
  Double_t  meanbin_l,meanbin_r,meanbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax;
  bool debug=false;
  Double_t max=0,tmax=0;
  rxmin=0;
  rxmax=0.5;

 
  const Int_t  nbinx=100, nbiny=500;

  rymin_l=12;
  rymax_l=25;
  rymin_r=12;
  rymax_r=20;
  tymin=13;
  tymax=20;
 


  Float_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Float_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  

  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

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
  

  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l); 
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx,-0.4,0.8,nbinx,tymin,tymax); 
 
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
  
 
  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;
  
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);
    
    
    meanbin_l=histotemp_l->GetMean();
    meanbin_r=histotemp_r->GetMean();
    meanbin_t=histotemp_t->GetMean();
    


    xt[k]=-0.4+(Float_t)(0.8-(-0.4))/nbinx*k;
    yt[k]=meanbin_t;
    rmsyt[k]=histotemp_t->GetRMS();
       
    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=meanbin_l;
    rmsy_l[k]=histotemp_l->GetRMS();

    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=meanbin_r;
    rmsy_r[k]=histotemp_r->GetRMS();

    
    
    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    
    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k

  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,550);

  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);

  TF1* hyp_r = new TF1("hyp_r","[0] - [1]/(x**[3]-[2])",0,0.4);
  TF1* hyp_l = new TF1("hyp_l","[0] - [1]/(x**[3]-[2])",0.0,0.45);
  TF1* hyp_t = new TF1("hyp_t","[0] + [1]*x",-0.1,1);
  //  hyp_r->SetParameter(0,10);
  // hyp_l->SetParameter(0,10);

  
  hyp_r->SetParameter(0,1e-3);
  hyp_r->SetParameter(1,-6);
  hyp_r->SetParameter(2,-3e-1);
  hyp_r->SetParameter(3,1.3);

  hyp_l->SetParameter(0, -6e1);
  hyp_l->SetParameter(1, -2e1);
  hyp_l->SetParameter(2, 7e-1);
  hyp_l->SetParameter(3, 2e-2);
  
  
  
  


  // hyp_r->SetParLimits(0,1,8);
  gStyle->SetOptStat("");
  
 

  wf_c->Divide(3,1);
  
  wf_c->cd(1);
  h2_l->SetXTitle("max_amp left");
  h2_l->SetYTitle("t_l - t_MCP (ns)");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","R");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("P");
  
  
  wf_c->cd(2);
  h2_r->SetXTitle("max_amp right");
  h2_r->SetYTitle("t_r - t_MCP (ns)");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","R");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);  
  graph_r->Draw("P");

  wf_c->cd(3);
  h2_t->SetXTitle("t_left-t_right (ns)");
  h2_t->SetYTitle("t_l - t_MCP (ns)");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","R");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);  
  graph_t->Draw("P");
  
  
  wf_c->SaveAs("plot.pdf");
  
  
}


