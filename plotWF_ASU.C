//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_ASU(const char * filename){


  TFile*  file= TFile::Open(filename);
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=-1;
  rxmax=1;



  const Int_t  nbinx=60,nbiny=100;


  rymin_l=0;
  rymax_l=0.3;
  rymin_r=0;
  rymax_r=0.3;
 
 
  Double_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[2][nbinx];


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
	h2_l->Fill(time[1+LEDi]-time[2+LEDi],amp_max[4]/max);
	h2_r->Fill(time[1+LEDi]-time[2+LEDi],amp_max[3]/max);
	
      }//chiudo if
  }//chiudo for k
  


  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_t;
  
  for(k=0;k<nbinx;k++){
        
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
          
    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();
    

 
    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();
    


    delete histotemp_l;
    delete histotemp_r;


    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k

  

  TGraphErrors* graph_rt=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_lt=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);  
  
  TCanvas* imma = new TCanvas("mycanv","title",1000,600);
  TLegend* l1=new TLegend(0.7,.8,0.7,.8);
  l1->SetHeader("time stamps","C");
  l1->AddEntry(graph_rt,"SiPM right");
  l1->AddEntry(graph_lt,"SiPM left");
  
  graph_rt->GetXaxis()->SetTitle("t_left-t_right [ns]");
  graph_rt->GetYaxis()->SetTitle("amp_max [mV]");
  
  graph_rt->SetMarkerStyle(8);
  graph_lt->SetMarkerStyle(8);
  
  graph_rt->SetMarkerSize(.8);
  graph_lt->SetMarkerSize(.8);
  
  graph_rt->SetMarkerColor(kRed);
  graph_lt->SetMarkerColor(kBlue);
  
  graph_rt->Draw("AP");
  graph_lt->Draw("SAMEP");
  
  l1->Draw();
}
