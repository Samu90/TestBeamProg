//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void GraficiTdiff(const char * filename){


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



  const Int_t  nbinx=200,nbiny=150;

 int i;
 Double_t sigma[50],erry[50],cut[50],errx[50];
  rymin_l=7.6;
  rymax_l=8.8;
  rymin_r=7.6;
  rymax_r=8.8;
  
  rymin_lc=7.4;
  rymax_lc=8;
  rymin_rc=7.4;
  rymax_rc=8;

  
  tymin=7.7;
  tymax=8.4;
  
  tymin_c=7.4;

  tymax_c=8;

    txmin=-0.3;
    txmax=0.8;
  
  

  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.14,1);
  TF1 *fit_l = new TF1("f_l","landau",0.14,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED300",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);

  digiTree->GetEntry(3);
  LEDi=LED300;
  
  max=4096;
  /*for(k=0; k<digiTree->GetEntries(); k++){
    digiTree->GetEntry(k);
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
  //if(time[3]-time[4]>tmax && time[3]-time[4]<10) {tmax = time[3]-time[4];}

  } chiudo for */

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


    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	if (amp_max[4]/max < 0.35)h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      

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
  /*
  for(k=0;k<nbinx;k++){
    for(j=0;j<nbiny;j++){
      //      if (k>20 && k<70) cout <<"  "<< rymin_l+(rymax_l-rymin_l)/nbiny*j << "<" << y_l[k]-3*RMS[0][k] <<"     "<< rymin_l+(rymax_l-rymin_l)/nbiny*j << ">" << y_l[k]+3*RMS[0][k] <<endl; 
      if (rymin_l+(rymax_l-rymin_l)/nbiny*j < y_l[k]-3*RMS[0][k] || rymin_l+(rymax_l-rymin_l)/nbiny*j > y_l[k]+3*RMS[0][k] )h2_l->SetBinContent(k,j,0);
      if (rymin_r+(rymax_r- rymin_r)/nbiny*j < y_r[k]-3*RMS[1][k] || rymin_r+(rymax_r-rymin_r)/nbiny*j > y_r[k]+3*RMS[1][k] ) h2_r->SetBinContent(k,j,0);
      // if (tymin+(tymax-tymin)/nbiny*j < yt[k]-3*RMS[2][k] || tymin+(tymax-tymin)/nbiny*j > yt[k]+3*RMS[2][k] ) h2_t->SetBinContent(k,j,0);
    }
    }*/
  
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
  hyp_r->SetParameter(0, 8.51);


  /* hyp_r->SetParameter(1, -1.54e1);
  hyp_r->SetParameter(2, 4.28e-2);
  hyp_r->SetParameter(3, -2.43e-2);
  */

 
 wf_c->Divide(3,2);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");

   h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0RL");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("SAMEP");
  hyp_l->Draw("same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","R0L");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
  hyp_r->Draw("same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","RL0M");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");
  hyp_t->Draw("same");



  
  TH2F* hc_l= new TH2F("hc_l", "histo hc_l",nbinx,txmin,txmax,nbiny,rymin_lc,rymax_lc);
  TH2F* hc_r= new TH2F("hc_r", "histo hc_r",nbinx,txmin,txmax,nbiny,rymin_rc,rymax_rc);
  TH2F* hc_t= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2F* hc_tdiff= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  
  
  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
        hc_r->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));
	hc_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0))/2);

	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      }

  }//chiudo for k

    for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
    histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
    histotemp_t=hc_t->ProjectionY("hc_tprojY",k,k);
   
    xt[k]=txmin+(txmax-txmin)/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();
    
    x_l[k]=txmin+(txmax-txmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    RMS[0][k]= histotemp_l->GetMeanError();
    
    x_r[k]=txmin+(txmax-txmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    RMS[1][k]= histotemp_r->GetMeanError();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    
  }//chiudo for k


   TGraphErrors* graph_lc = new TGraphErrors(nbinx-1,x_l,y_l,0,RMS[0]);
   TGraphErrors* graph_rc = new TGraphErrors(nbinx-1,x_r,y_r,0,RMS[1]);
   TGraphErrors* graph_tc = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x",-0.1,0.65);
   TF1* fit_lc = new TF1("fit_lc","[0]+[1]*x",-0.1,0.65);
   TF1* fit_rc = new TF1("fit_rc","[0]+[1]*x",-0.1,0.65);
  
   graph_tc->Fit("fit_tdiff","R");
   graph_rc->Fit("fit_rc","R");
   graph_lc->Fit("fit_lc","R");


   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_tdiff->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi])+fit_tdiff->GetParameter(0)-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0))/2);
	
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      }

  }//chiudo for k
    for(k=0;k<nbinx;k++){
   
   
    TH1D* histotemp_t;

   
   
    histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
   
   
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();

    delete histotemp_t;

    
  }//chiudo for k

    TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
    wf_c->cd(4);
    hc_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
    hc_l->GetXaxis()->SetTitle("t_left-t_right [mV]");
    hc_l->Draw("COLZ");
    graph_lc->SetMarkerStyle(8);
    graph_lc->SetMarkerSize(.5);
    graph_lc->Draw("P");

  // graph_l->Fit("hyp_l","R");
   //   graph_l->SetMarkerStyle(8);
   // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

   wf_c->cd(5);

   hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
   hc_r->GetXaxis()->SetTitle("t_left-t_right [mV]");
   hc_r->Draw("COLZ");
   graph_rc->SetMarkerStyle(8);
   graph_rc->SetMarkerSize(.5);
   graph_rc->Draw("P");

  // graph_l->Fit("hyp_l","R");
   // graph_l->SetMarkerStyle(8);
   // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

   wf_c->cd(6);

   hc_tdiff->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
   hc_tdiff->GetXaxis()->SetTitle("t_left-t_right [ns]");
   hc_tdiff->Draw("COLZ");
   graph_tcdiff->SetMarkerStyle(8);
   graph_tcdiff->SetMarkerSize(.5);
   graph_tcdiff->Draw("P");

  // graph_l->Fit("hyp_l","R");
  // graph_l->SetMarkerStyle(8);
  // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

   TCanvas* plottini = new TCanvas("altroCanvas","",800,600);
   plottini->Divide(2,1);
   hc_r->Reset();
   hc_l->Reset();
   
   for(k=0;k<digiTree->GetEntries();k++){
     digiTree->GetEntry(k);
     
     if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
       {
	 hc_l->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-fit_lc->Eval((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))))+fit_lc->GetParameter(0));
	 
	 hc_r->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0)-fit_rc->Eval((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))))+fit_rc->GetParameter(0));
	 
	 if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
       }
     
   }//chiudo for k
   
   plottini->cd(1);
   hc_l->Draw("COLZ");
   plottini->cd(2);
   hc_r->Draw("COLZ");
   
   TH1D* histo_cl;
   TH1D* histo_cr;
   TH1D* histo_ct;
   TH1D* histo_ctdiff;
   TF1* gaus_cl = new TF1("gaus_cl","gaus",-2.5,-0.7);
   TF1* gaus_cr = new TF1("gaus_cr","gaus",-2.5,-0.7);
   TF1* gaus_ct = new TF1("gaus_ct","gaus",-2.5,-0.7);
   TF1* gaus_ctdiff = new TF1("gaus_ctdiff","gaus",6,8);
   histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
   histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);
   histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);
   histo_ctdiff = hc_tdiff->ProjectionY("histo_ctdiff",0,nbinx);


   histo_ct->SetLineColor(kBlack);
   histo_cl->SetLineColor(kBlue);
   histo_cr->SetLineColor(kRed);
   gaus_ct->SetLineColor(kBlack);

   TCanvas * timeres = new TCanvas("timeres","plot_timeres",600,550);

   TLegend* l1=new TLegend(0.1,0.7,0.48,0.9);
   l1->SetHeader("time stamps");
   l1->AddEntry(histo_cl,"t_left-t_MCP");
   l1->AddEntry(histo_cr,"t_right-t_MCP");
   l1->AddEntry(histo_cr,"t_ave-t_MCP");
   
   gStyle->SetOptStat("");
   histo_ct->Draw();
   gaus_ct->SetParameter(0,500);
   histo_ct->GetYaxis()->SetTitle("counts");
   histo_ct->GetXaxis()->SetTitle("t_{ave}-t_{MCP}(ns)");
   histo_cl->SetLineColor(kBlue);

   histo_cl->Fit("gaus_cl");
   histo_cr->Fit("gaus_cr");
   if (blind==true) histo_ct->Fit("gaus_ct");

   histo_cr->Draw("same");
  
   histo_cl->Draw("same");

   l1->Draw();

   TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
   TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
   l2->SetHeader("time stamps");
   l2->AddEntry(histo_ct,"t_ave-t_MCP");
   l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
   gaus_ctdiff->SetLineColor(kGreen);
   histo_ctdiff->Fit("gaus_ctdiff");
   histo_ctdiff->SetLineColor(kGreen);
   histo_ct->Draw();
   histo_ctdiff->Draw("same");
   l2->Draw();
 
   
   cut[0]=0.4;
   cut[1]=0.2;
   cut[2]=0.1;
   cut[3]=0.05;
   cut[4]=0.02;
   cut[5]=0.01;
   
   TCanvas* canvino = new TCanvas("canvino","",1200,600);
   TF1* fittino[3][6];
   canvino->Clear();
   canvino->Divide(3,2);
   
   TH1D* istogrammi[3][6];
   
   for (i=0;i<6;i++){
     fittino[0][i] = new TF1("argaerg"+i ,"gaus",6,8);
     fittino[1][i] = new TF1("xzcbxcn"+i ,"gaus",6,8);
     fittino[2][i] = new TF1("jklnjbm"+i ,"gaus",6,8);
     
     istogrammi[0][i]=hc_tdiff->ProjectionY("ghijklxz"+i,0 /*hc_tdiff->GetXaxis()->FindBin(0.0)*/, hc_tdiff->GetXaxis()->FindBin(cut[i]));
     istogrammi[1][i]=hc_l->ProjectionY("abcdefuv"+i, 0/*hc_tdiff->GetXaxis()->FindBin(0.0)*/, hc_tdiff->GetXaxis()->FindBin(cut[i]));
     istogrammi[2][i]=hc_r->ProjectionY("mnopqrst"+i, 0/*hc_tdiff->GetXaxis()->FindBin(0.0)*/, hc_tdiff->GetXaxis()->FindBin(cut[i]));

     istogrammi[0][i]->SetLineColor(kBlack);
     istogrammi[1][i]->SetLineColor(kBlue);
     istogrammi[2][i]->SetLineColor(kRed);
     
     cout<<"__________________"<< hc_tdiff->GetXaxis()->FindBin(0.0)<<"     "<<hc_tdiff->GetXaxis()->FindBin(cut[i])<<endl;

     //canvino->cd(i+1);
     istogrammi[0][i]->Fit("argaerg"+i,"R0");
     istogrammi[1][i]->Fit("xzcbxcn"+i,"R0");
     istogrammi[2][i]->Fit("jklnjbm"+i,"R0");
     
     //istogrammi[0][i]->Draw();
     //istogrammi[1][i]->Draw("SAME");
     //istogrammi[2][i]->Draw("SAME");
     
     cout<< fittino[0][0]->GetParameter(2)<<"   "<< istogrammi[0][i]->GetName()<<  " "<<istogrammi[0][i]->GetEntries()<<endl;

     
   }

   canvino->cd(1);
   istogrammi[0][0]->Draw();
   istogrammi[1][0]->Draw("SAME");
   istogrammi[2][0]->Draw("SAME");
   fittino[0][0]->Draw("SAME");
   fittino[1][0]->Draw("SAME");
   fittino[2][0]->Draw("SAME");
   
   canvino->cd(2);
   istogrammi[0][1]->Draw();
   istogrammi[1][1]->Draw("SAME");
   istogrammi[2][1]->Draw("SAME");
   fittino[0][1]->Draw("SAME");
   fittino[1][1]->Draw("SAME");
   fittino[2][1]->Draw("SAME");
   
   canvino->cd(3);
   istogrammi[0][2]->Draw();
   istogrammi[1][2]->Draw("SAME");
   istogrammi[2][2]->Draw("SAME");
   fittino[0][2]->Draw("SAME");
   fittino[1][2]->Draw("SAME");
   fittino[2][2]->Draw("SAME");
   
   canvino->cd(4);
   istogrammi[1][3]->Draw();
   istogrammi[2][3]->Draw("SAME");
   istogrammi[0][3]->Draw("SAME");
   fittino[0][3]->Draw("SAME");
   fittino[1][3]->Draw("SAME");
   fittino[2][3]->Draw("SAME");
   
   canvino->cd(5);
   istogrammi[1][4]->Draw();
   istogrammi[2][4]->Draw("SAME");
   istogrammi[0][4]->Draw("SAME");
   fittino[0][4]->Draw("SAME");
   fittino[1][4]->Draw("SAME");
   fittino[2][4]->Draw("SAME");
   
   canvino->cd(6);
   istogrammi[1][5]->Draw();
   istogrammi[2][5]->Draw("SAME");
   istogrammi[0][5]->Draw("SAME");
   fittino[0][5]->Draw("SAME");
   fittino[1][5]->Draw("SAME");
   fittino[2][5]->Draw("SAME");
}
