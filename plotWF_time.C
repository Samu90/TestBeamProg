//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_time(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0;
  rxmax=0.5;


  const Int_t  nbinx=200,nbiny=500;

  rymin_l=7.6;
  rymax_l=8.8;
  rymin_r=7.6;
  rymax_r=8.8;
 
  
  rymin_lc=-5;
  rymax_lc=5;
  rymin_rc=-5;
  rymax_rc=5;

  tymin=7;
  tymax=9;

  txmin=-0.6;
  txmax=0.6;


  Float_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Float_t xt[nbinx],yt[nbinx],rmsyt[nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED300",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);

  digiTree->GetEntry(3);
  LEDi=LED300;
  
  max=4096;

  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_tt",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);
      }//chiudo if
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      

  }//chiudo for k


  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("ht_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("ht_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("ht_projY",k,k);


    //    maxbin_l=histotemp_l->GetMean();
    //  maxbin_r=histotemp_r->GetMaximumBin();
    // maxbin_t=histotemp_t->GetMaximumBin();



    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();

    x_l[k]=rxmin+(rxmax-rxmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();

    x_r[k]=rxmin+(rxmax-rxmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();



    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k



  // TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,550);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
 
  
  TF1* hyp_r = new TF1("hyp_r","[0]-[1]*log(x+[2])",0.125,0.35);
  TF1* hyp_l = new TF1("hyp_l","[0]-[1]*log(x+[2])",0.11,0.35);
  TF1* hyp_t = new TF1("hyp_t","[0]*x+[1]",0.7,txmax);
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 1.2);
  hyp_l->SetParameter(2, 5);
  /* hyp_l->SetParameter(3, -2.43e-2);
  */
  hyp_r->SetParameter(0, 8.51);

  // gStyle->SetOptStat("");



  /*  wf_c->Divide(3,1);
   wf_c->cd(1);
   h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
   h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
   h2_l->Draw("COLZ");*/
   graph_l->Fit("hyp_l","0RL");
   /*  graph_l->SetMarkerStyle(8);
   graph_l->SetMarkerSize(.5);
   graph_l->Draw("P");
   hyp_l->Draw("same");
  
   

   wf_c->cd(2);
   h2_r->Draw("COLZ");*/
   graph_r->Fit("hyp_r","0RL");
   /* graph_r->SetMarkerStyle(8);
   graph_r->SetMarkerSize(.5);
   graph_r->Draw("P");
   hyp_r->Draw("same");
   wf_c->cd(3);
   h2_t->Draw("COLZ");
   // graph_t->Fit("hyp_t","RL");

   graph_t->SetMarkerStyle(8);
   graph_t->SetMarkerSize(.5);
   // graph_t->Draw("P");
   */
   delete graph_l;
   delete graph_r;
   delete graph_t;

   delete h2_l;
   delete h2_r;
   delete h2_t;

   
  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_tt",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      {
	h2_l->Fill(time[1+LEDi]-time[2+LEDi],time[1+LEDi]-time[0]);
	h2_r->Fill(time[1+LEDi]-time[2+LEDi],time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);
      }//chiudo if
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      

  }//chiudo for k


  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("ht_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("ht_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("ht_projY",k,k);


    //    maxbin_l=histotemp_l->GetMean();
    //  maxbin_r=histotemp_r->GetMaximumBin();
    // maxbin_t=histotemp_t->GetMaximumBin();



    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();

    x_l[k]=txmin+(Float_t)(txmax-txmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();

    x_r[k]=txmin+(Float_t)(txmax-txmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();



    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  

  TCanvas* time_p =new TCanvas("time","Plot time",600,550);
  TLegend* l1= new TLegend(0.1,0.7,0.48,0.9);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);

   
   
  l1->AddEntry(graph_l,"t_{left}-t_{MCP}","P");
  l1->AddEntry(graph_r,"t_{right}-t_{MCP}","P");
  l1->AddEntry(graph_t,"t_{ave}-t_{MCP}","P");
  graph_l->GetYaxis()->SetRangeUser(7.5,8.6);
  graph_l->GetXaxis()->SetRangeUser(-0.08,0.6);
  graph_l->GetYaxis()->SetTitle("t-t_{MCP}(ns)");
  graph_l->GetXaxis()->SetTitle("t_{left}-t_{right}(ns)");
  graph_l->SetMarkerColor(kRed);
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("AP");
  
  graph_r->SetMarkerColor(kBlue);
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
    
    
    
  graph_t->SetMarkerColor(kBlack);
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");

  l1->Draw();
  
   TH2F* ht_l= new TH2F("hc_l", "histo hc_l",nbinx,txmin,txmax,nbiny,tymin,tymax);
   TH2F* ht_r= new TH2F("hc_r", "histo hc_r",nbinx,txmin,txmax,nbiny,tymin,tymax);
   TH2F* ht= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin,tymax);


  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      {
	ht_l->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
	
        ht_r->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));
	ht->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[3]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[4]/max)-hyp_l->GetParameter(0))/2);
	
	

	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      }

  }//chiudo for k
   
    for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=ht_l->ProjectionY("ht_lprojY",k,k);
    histotemp_r=ht_r->ProjectionY("ht_rprojY",k,k);
    histotemp_t=ht->ProjectionY("ht_projY",k,k);


    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();

    x_l[k]=txmin +(Float_t)(txmax-txmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();

    x_r[k]=txmin+(Float_t)(txmax-txmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();



    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k



  TCanvas* time_pc =new TCanvas("timecorr","Plot timecorr",600,550);
  TLegend* l2= new TLegend(0.1,0.7,0.48,0.9);
  TGraphErrors* g_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* g_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* g_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);

   
   
  l2->AddEntry(graph_l,"t_{left}-t_{MCP}","P");
  l2->AddEntry(graph_r,"t_{right}-t_{MCP}","P");
  l2->AddEntry(graph_t,"t_{ave}-t_{MCP}","P"); 
  g_l->GetYaxis()->SetRangeUser(7.3,8.6);
  g_l->GetXaxis()->SetRangeUser(-0.08,0.6);
  g_l->GetYaxis()->SetTitle("t-t_{MCP}(amp.walk corr)(ns)");
  g_l->GetXaxis()->SetTitle("t_{left}-t_{right}(amp.walk corr)(ns)");
  g_l->SetMarkerColor(kRed);
  g_l->SetMarkerStyle(8);
  g_l->SetMarkerSize(.5);
  g_l->Draw("AP");
  
  g_r->SetMarkerColor(kBlue);
  g_r->SetMarkerStyle(8);
  g_r->SetMarkerSize(.5);
  g_r->Draw("SAMEP");
    
    
    
  g_t->SetMarkerColor(kBlack);
  g_t->SetMarkerStyle(8);
  g_t->SetMarkerSize(.5);
  g_t->Draw("SAMEP");
  l2->Draw();
    /*  TH1D* histo_cl;
   TH1D* histo_cr;
   TH1D* histo_ct;
   TF1* gaus_cl = new TF1("gaus_cl","gaus",-2.5,-0.7);
   TF1* gaus_cr = new TF1("gaus_cr","gaus",-2.5,-0.7);
   TF1* gaus_ct = new TF1("gaus_ct","gaus",-2.5,-0.7);
   histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
   histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);
   histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);


   histo_ct->SetLineColor(kBlack);
   histo_cl->SetLineColor(kBlue);
   histo_cr->SetLineColor(kRed);

   TCanvas * timeres = new TCanvas("timeres","plot_timeres",600,550);
   gStyle->SetOptStat("");
   histo_ct->Draw();
   gaus_ct->SetParameter(0,500);
   histo_cl->Fit("gaus_cl");
   histo_cr->Fit("gaus_cr");
   histo_ct->Fit("gaus_ct");

   histo_cr->Draw("same");
  
   histo_cl->Draw("same");
  */
}

