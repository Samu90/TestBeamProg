//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void TResAmp(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,j,maxbin_l,maxbin_r,maxbin_t;
  Double_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0.8;
  rxmax=1.5;



  const Int_t  nbinx=200,nbiny=150;

 int i;
 Double_t sigma[50],erry[50],cut[50],errx[50];
  rymin_l=7.6;
  rymax_l=8.8;
  rymin_r=7.6;
  rymax_r=8.8;
  
  rymin_lc=-1;
  rymax_lc=2;
  rymin_rc=-1;
  rymax_rc=2;

  tymin=7.7;
  tymax=8.4;
  
  tymin_c=6.8;

  tymax_c=8.9;


  txmin=-0.3;
  txmax=0.8;


  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];


  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1D *mcp_amp =new TH1D("mcp_amp","histomcp_ampl",nbinx,0.0,1);


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


  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* h2_t= new TH2D("h2_t", "histo h2_tot",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);


    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/(max*fit_l->GetParameter(1)),time[1+LEDi]-time[0]);
	h2_r->Fill(amp_max[4]/(max*fit_r->GetParameter(1)),time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
	

  }//chiudo for k
cout<<"HERE"<<endl;
  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    //    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);
   
    /*   xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();
    RMS[2][k]= histotemp_t->GetRMS();
    */
    x_l[k]=rxmin+(rxmax-rxmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();
    RMS[0][k]= histotemp_l->GetRMS();
    
    
    x_r[k]=rxmin+(rxmax-rxmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();
    RMS[1][k]= histotemp_r->GetRMS();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  cout<<"HERE"<<endl;
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

  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x+[1])",0.8,3);
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x+[1])",0.8,3);

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");


  /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  // hyp_l->SetParameter(1, 5);
  // hyp_l->SetParameter(2, 1.2);
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



  
  TH2F* hc_l= new TH2F("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2F* hc_r= new TH2F("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2F* hc_t= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2F* hc_tdiff= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tl= new TH2D("hc_tl", "histo hc_tl",nbinx,0.8,3,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_tr= new TH2D("hc_tr", "histo hc_tr",nbinx,0.8,3,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_tot= new TH2D("hc_tr", "histo hc_tr",nbinx,0.8,3,nbiny,rymin_lc,rymax_lc);
    
  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill(amp_max[3]/(max*fit_l->GetParameter(1)),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))));
	hc_r->Fill(amp_max[4]/(max*fit_r->GetParameter(1)),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))));
	hc_t->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1)))-time[2+LEDi]+hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1)))+hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))))/2);

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
   
    xt[k]=txmin +(Float_t)(txmax-txmin)*k/nbinx;
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();

    x_l[k]=rxmin+(Float_t)(rxmax-rxmin)*k/nbinx;
    y_l[k]=histotemp_l->GetMean();
    RMS[0][k]= histotemp_l->GetMeanError();

    x_r[k]=rxmin+(Float_t)(rxmax-rxmin)*k/nbinx;
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
   // fit_tdiff->SetParameter(1,8);
   graph_tc->Fit("fit_tdiff");

   wf_c->cd(4);
    hc_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
    hc_l->GetXaxis()->SetTitle("max.amplitude [mV]");
    hc_l->Draw("COLZ");
    graph_lc->SetMarkerStyle(8);
    graph_lc->SetMarkerSize(.5);
    graph_lc->Draw("P");
 
    for(k=0;k<nbinx;k++){
   
   
    TH1D* histotemp_t;

   
   
    histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
   
   
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();

    delete histotemp_t;

    
  }//chiudo for k

    TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   


 
  // graph_l->Fit("hyp_l","R");
   //   graph_l->SetMarkerStyle(8);
   // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

  wf_c->cd(5);

   hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
   hc_r->GetXaxis()->SetTitle("max.amplitude [mV]");
   hc_r->Draw("COLZ");
   graph_rc->SetMarkerStyle(8);
   graph_rc->SetMarkerSize(.5);
   graph_rc->Draw("P");

     for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	
	hc_tdiff->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1)))-time[2+LEDi]+hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))),(time[1+LEDi]+time[2+LEDi])/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi])+fit_tdiff->GetParameter(0)-time[0]-(hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1)))+hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))))/2);
	hc_tl->Fill(amp_max[3]/(max*fit_l->GetParameter(1)),(time[1+LEDi]+time[2+LEDi])/2/*-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi])*/-time[0]-(hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1)))+hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))))/2);
	hc_tr->Fill(amp_max[4]/(max*fit_r->GetParameter(1)),(time[1+LEDi]+time[2+LEDi])/2/*-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi])*/-time[0]-(hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1)))+hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))))/2);

	
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      }

  }//chiudo for k
    hc_tot->Add(hc_tl,hc_tr);
    
   for (i=0;i<nbinx/10;i++){
     cut[i] =0.8+(Float_t)(1.5-0.8)*(i*10)/nbinx;  
   }
   for (i=0;i<nbinx/10;i++){
     TF1* fit = new TF1("fit","gaus",-1,2);
     TH1D* histotemp_t;
     histotemp_t=hc_tot->ProjectionY("hc_totprojY",hc_tot->GetXaxis()->FindBin(cut[i]),hc_tot->GetXaxis()->FindBin(cut[i+1]));
     cout << "____________" << hc_tot->GetXaxis()->FindBin(cut[i]) << hc_tot->GetXaxis()->FindBin(cut[i+1]) << endl;
     
    histotemp_t->Fit("fit","0");
    sigma[i]=sqrt((fit->GetParameter(2))*(fit->GetParameter(2))-0.015*0.015);
    erry[i]=fit->GetParError(2);
    errx[i]= (1.5-0.8)*10/(2*nbinx);
    
    delete histotemp_t;
    delete fit;
    // delete fit;
   }
				     
   TCanvas* c_TresAmp = new TCanvas("c_rest","c_rest_plot",600,550);
   TGraphErrors* TResAmp = new TGraphErrors(nbinx/10,cut,sigma,errx,erry);
   TF1* fitramp = new TF1("fitramp","[0]+[1]/sqrt(x)",rxmin,rxmax-0.06);
   gStyle->SetOptFit(1111);
   fitramp->SetParameter(0,0.029);
   fitramp->SetParameter(1,1e-2);
   TResAmp->Fit("fitramp","L0R");
   TResAmp ->GetXaxis()->SetTitle("amp/mip peak");
   TResAmp ->GetYaxis()->SetTitle("sigma_{t_{ave}}(ns)");
   TResAmp ->SetMarkerStyle(8);
   TResAmp ->SetMarkerSize(.8);
   
  
   TResAmp ->Draw("AP");
   fitramp->DrawF1(rxmin,3,"same");

   // graph_l->Fit("hyp_l","R");
   // graph_l->SetMarkerStyle(8);
   // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

  
  // graph_l->Fit("hyp_l","R");
  // graph_l->SetMarkerStyle(8);
  // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

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
   
}
