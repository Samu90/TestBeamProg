//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_corr13(const char * filename){


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


  const Int_t  nbinx=100,nbiny=600;

  rymin_l=2;
  rymax_l=7;
  rymin_r=0;
  rymax_r=7;
  
  rymin_lc=-5;
  rymax_lc=5;
  rymin_rc=-5;
  rymax_rc=5;

  tymin=1.5;
  tymax=7;

  txmin=0;
  txmax=3;


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
  LEDi=LED30;
  
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

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      {
	if(amp_max[3]/max<0.35) h2_l->Fill(amp_max[3]/max,time[3+LEDi]-time[0]);
	if(amp_max[4]/max<0.35) h2_r->Fill(amp_max[4]/max,time[4+LEDi]-time[0]);
	h2_t->Fill((time[3+LEDi]-time[4+LEDi]),(time[3+LEDi]+time[4+LEDi])/2-time[0]);

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


    maxbin_l=histotemp_l->GetMaximumBin();
    maxbin_r=histotemp_r->GetMaximumBin();
    maxbin_t=histotemp_t->GetMaximumBin();



    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=tymin+(Float_t)(tymax-tymin)/nbiny*maxbin_t;
    rmsyt[k]=histotemp_t->GetRMS();

    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=rymin_l+(rymax_l-rymin_l)/nbiny*maxbin_l;
    rmsy_l[k]=histotemp_l->GetRMS();

    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=rymin_r+(rymax_r-rymin_r)/nbiny*maxbin_r;

    rmsy_r[k]=histotemp_r->GetRMS();



    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k



  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,1100);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
  TF1* hyp_r = new TF1("hyp_r","[0]-[1]/(x**[2]+[3])",0.125,0.35);
  TF1* hyp_l = new TF1("hyp_l","[0]-[1]/(x**[2]+[3])",0.135,0.35);
  TF1* hyp_t = new TF1("hyp_t","[0]*x+[1]",0.7,txmax);
  
  gStyle->SetOptStat("");



  wf_c->Divide(3,2);

  wf_c->cd(1);
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","RL");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("P");


  wf_c->cd(2);
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","RL");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("P");

  wf_c->cd(3);
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","RL");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("P");

  TH2F* hc_l= new TH2F("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2F* hc_r= new TH2F("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2F* hc_t= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,-3,4.5);


  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      {
	if(amp_max[3]/max<0.35) hc_l->Fill(amp_max[3]/max,time[3+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
	if(amp_max[4]/max<0.35) hc_r->Fill(amp_max[4]/max,time[4+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));
	hc_t->Fill((time[3+LEDi]-time[4+LEDi]),(time[3+LEDi]+time[4+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[3]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[4]/max)-hyp_r->GetParameter(0))/2);


	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      }

  }//chiudo for k


    wf_c->cd(4);
    hc_l->Draw("COLZ");
  // graph_l->Fit("hyp_l","R");
   //   graph_l->SetMarkerStyle(8);
   // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

   wf_c->cd(5);
   hc_r->Draw("COLZ");
  // graph_l->Fit("hyp_l","R");
   // graph_l->SetMarkerStyle(8);
   // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

   wf_c->cd(6);
   hc_t->Draw("COLZ");
  // graph_l->Fit("hyp_l","R");
  // graph_l->SetMarkerStyle(8);
  // graph_l->SetMarkerSize(.5);
  // graph_l->Draw("P");

   TH1D* histo_cl;
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
   
  }

