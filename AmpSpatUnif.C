//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void AmpSpatUnif(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0;
  rxmax=0.5;

  Int_t nentries=digiTree->GetEntries(), counter1=0, counter2=0, counter3=0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];

  const Int_t  nbinx=100,nbiny=120;


  
  txmin=-0.3;
  txmax=0.8;


  Double_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.08,1);
  TF1 *fit_l = new TF1("f_l","landau",0.08,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED300",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);

  digiTree->GetEntry(3);
  LEDi=LED300;
  
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    //cout<<"HERE"<<endl;
    if(time[1+LEDi]-time[0]<15 && time[1+LEDi]-time[0]>0) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0];
    }
    //  else{Times1[k]=Times2[k-1];}
    if(time[2+LEDi]-time[0]<15 && time[2+LEDi]-time[0]>0) {
      counter2++;
      Times2[counter2]=time[2+LEDi]-time[0];}
    // else{Times2[k]=Times2[k-1];}  
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<15 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>-10) {
      counter3++;
      Times3[counter3]=(time[1+LEDi]+time[2+LEDi])/2-time[0];}
    // else{Times3[k]=Times3[k-1];}
  }
  
  Double_t mean1=TMath::Mean(counter1,Times1);
  Double_t rms1=TMath::RMS(counter1,Times1);
  cout<<mean1<<"________"<<rms1<<endl;
  Double_t mean2=TMath::Mean(counter2,Times2);
  Double_t rms2=TMath::RMS(counter2,Times2);
  cout<<mean2<<"________"<<rms2<<endl;
  Double_t mean3=TMath::Mean(counter3,Times3);
  Double_t rms3=TMath::RMS(counter3,Times3);
  cout<<mean3<<"________"<<rms3<<endl;
  
  rymin_l=mean1-1.1*rms1;
  rymax_l=mean1+0.8*rms1;
  rymin_r=mean2-1.1*rms2;
  rymax_r=mean2+0.8*rms2;
    
  

  tymin=mean3-1.1*rms3;
  tymax=mean3+0.8*rms3;
  
  


  txmin=-0.3;
  txmax=0.8;

  
  max=4096;


  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
    mcp_amp->Fill(amp_max[0]/max);

  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);


    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > mcp_amp->GetMean()-1*mcp_amp->GetRMS() && amp_max[0]/max < mcp_amp->GetMean()+1*mcp_amp->GetRMS())
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	if (amp_max[4]/max < 0.35)h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if

  }//chiudo for k


  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);
   
    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();
    RMS[2][k]= histotemp_t->GetRMS();

    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();
    RMS[0][k]= histotemp_l->GetRMS();
    
    
    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();
    RMS[1][k]= histotemp_r->GetRMS();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
    
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,1100);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);

  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x+[1])",0.135,0.35);
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x+[1])",0.11,0.35);

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");


  /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  /* hyp_l->SetParameter(3, -2.43e-2);
  */
  hyp_r->SetParameter(0, 7);
  hyp_r->SetParameter(1, -9e-2);
  hyp_r->SetParameter(2, -1e-1);


  wf_c->Divide(3,2);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
  h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0RL");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("P");
  hyp_l->DrawF1(0,1,"same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","R0L");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("P");
  hyp_r->DrawF1(0,1,"same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","RL0");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("P");
  hyp_t->DrawF1(txmin,txmax,"same");

  rymin_lc=rymin_l-hyp_l->Eval(0.25)+hyp_l->GetParameter(0);
  rymax_lc=rymax_l-hyp_l->Eval(0.25)+hyp_l->GetParameter(0);
  rymin_rc=rymin_r-hyp_r->Eval(0.25)+hyp_r->GetParameter(0);
  rymax_rc=rymax_r-hyp_r->Eval(0.25)+hyp_r->GetParameter(0);
  tymin_c=tymin;
  tymax_c=tymax;

  
  TH2F* hc_l= new TH2F("hc_l", "histo hc_l",nbinx,-1,1,nbiny,0,1);
  TH2F* hc_r= new TH2F("hc_r", "histo hc_r",nbinx,-1,1,nbiny,0,1);
  


  
   for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
     {
       hc_l->Fill(time[1+LEDi]/*-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)*/-time[2+LEDi]/*-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))*/, amp_max[3]/max);
       hc_r->Fill(time[1+LEDi]/*-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)*/-time[2+LEDi]/*-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))*/, amp_max[4]/max);
    }
    
   }//chiudo for k
   cout<<"_______________________________________"<<hc_r->GetEntries()<<endl;
   wf_c->cd(4);
   hc_l->Draw("COLZ");
   wf_c->cd(5);
   hc_r->Draw("COLZ");
   
   Float_t xmin,xmax;
   xmin=-1;
   xmax=1;
   
   TH1D* histotemp_l;
   TH1D* histotemp_r;
   
   int newbin=nbinx;

   for(k=0;k<newbin;k++){
     
     histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
     histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
     
     x_l[k]=xmin+(xmax-xmin)/newbin*k;
     y_l[k]=histotemp_l->GetMean();
     RMS[0][k]= histotemp_l->GetMeanError();
     
     x_r[k]=xmin+(xmax-xmin)/newbin*k; 
     y_r[k]=histotemp_r->GetMean();
     RMS[1][k]= histotemp_r->GetMeanError();
     
     
     delete histotemp_l;
     delete histotemp_r;
     
     
     
   }//chiudo for k
   
   
  TGraphErrors* graph_rt=new TGraphErrors(newbin-1,x_r,y_r,0,RMS[0]);
  TGraphErrors* graph_lt=new TGraphErrors(newbin-1,x_l,y_l,0,RMS[1]);  
  
  //TCanvas* imma = new TCanvas("mycanv","title",1000,600);
  TLegend* l1=new TLegend(0.2,0.3,0.2,0.3);
  l1->SetHeader("SiPM","C");
  l1->AddEntry(graph_rt,"SiPM right");
  l1->AddEntry(graph_lt,"SiPM left");
  
  graph_rt->GetXaxis()->SetTitle("t_{left}-t_{right} [ns]");
  graph_rt->GetYaxis()->SetTitle("amp_{max} [mV]");
  
  graph_rt->SetMarkerStyle(8);
  graph_lt->SetMarkerStyle(8);
  
  graph_rt->SetMarkerSize(.8);
  graph_lt->SetMarkerSize(.8);
  
  graph_rt->SetMarkerColor(kRed);
  graph_lt->SetMarkerColor(kBlue);
  wf_c->cd(5);
  graph_rt->Draw("SAMEP");
  wf_c->cd(4);
  graph_lt->Draw("SAMEP");

  TCanvas* newcanv = new TCanvas();
  graph_rt->Draw("AP");
  graph_rt->GetXaxis()->SetLimits(-0.3,0.7);
  graph_rt->GetYaxis()->SetLimits(0.12,0.24);
  graph_lt->Draw("SAMEP");

  l1->Draw();
   
  }
