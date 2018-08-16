#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TSystem.h>


void Proiezione(TH2F* hl,TH2F* hr,TH2F* ht,Int_t nxbin,Float_t DCR){
  TH1D* proj[3];
  proj[0]=hl->ProjectionY("uncorrR",0,nxbin);
  proj[1]=hr->ProjectionY("uncorrL",0,nxbin);
  proj[2]=ht->ProjectionY("uncorrT",0,nxbin);
  
  TCanvas* pro = new TCanvas("leproj","",500,400);
  pro->Divide(3,1);
  pro->cd(1);
  proj[0]->GetXaxis()->SetTitle("t_{left}-t_{MCP}");
  proj[0]->GetYaxis()->SetTitle("counts");
  proj[0]->SetTitle("SiPM_left");
  proj[0]->Draw();
  
  pro->cd(2);
  proj[1]->GetXaxis()->SetTitle("t_{right}-t_{MCP}");
  proj[1]->GetYaxis()->SetTitle("counts");
  proj[1]->SetTitle("SiPM_right");
  proj[1]->Draw("SAME");

  pro->cd(3);
  proj[2]->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
  proj[2]->GetYaxis()->SetTitle("counts");
  proj[2]->SetTitle("t_{ave}");
  proj[2]->Draw("SAME");
  
  pro->SaveAs(("HDCRPlot/ProiezioniLR"+to_string((int)DCR)+".pdf").c_str());
  cout<<"PRIMA DEI DELETE"<<endl;
  
  pro->Close();
  delete pro;
   
}


void Ris(TFile* file,Float_t* Resol,Float_t* errResol,Float_t index){


  //TFile*  file= TFile::Open(filename);
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



  const Int_t  nbinx=250,nbiny=200;

  int i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-0.9;
  txmax=1.4;
  
  
  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];
  
  
  Int_t nentries=digiTree->GetEntries(), counter1=0,counter2=0, counter3=0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];
  
  for(k=0;k<digiTree->GetEntries();k++){
    Times1[k]=0;
    Times2[k]=0;
    Times3[k]=0;
  }  
  
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.04,1);
  TF1 *fit_l = new TF1("f_l","landau",0.04,1);


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
    
    if(time[1+LEDi]-time[0]<15 && time[1+LEDi]-time[0]>0) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0];
    }
    
    if(time[2+LEDi]-time[0]<15 && time[2+LEDi]-time[0]>0) {
      counter2++;
      Times2[counter2]=time[2+LEDi]-time[0];
    }
    
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<15 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>-1) {
      counter3++;
      Times3[counter3]=(time[1+LEDi]+time[2+LEDi])/2-time[0];
    }
    
  }
  
  Double_t mean1=TMath::Mean(counter1,Times1);
  Double_t rms1=TMath::RMS(counter1,Times1);
  cout<<mean1<<"_________"<<rms1<<endl;
  Double_t mean2=TMath::Mean(counter2,Times2);
  Double_t rms2=TMath::RMS(counter2,Times2);
  cout<<mean2<<"_________"<<rms2<<endl;
  Double_t mean3=TMath::Mean(counter3,Times3);
  Double_t rms3=TMath::RMS(counter3,Times3);
  cout<<mean3<<"_________"<<rms3<<endl;

  rymin_l=mean1-1.2*rms1;
  rymax_l=mean1+0.8*rms1;
  rymin_r=mean2-1.2*rms2;
  rymax_r=mean2+0.8*rms2;
    
  

  tymin=mean3-1.2*rms3;
  tymax=mean3+1.2*rms3;
 


  
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

  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");

  bool DebugLand=true;
  
  if(DebugLand){
    TCanvas* LandCanv = new TCanvas("mycanvas","",700,500);
    LandCanv->Divide(2,1);
    
    LandCanv->cd(1)->SetLogy();
    hl_amp->Draw();
    fit_l->Draw("SAME");
    
    LandCanv->cd(2)->SetLogy();
    hr_amp->Draw();
    fit_r->Draw("SAME");
    
    LandCanv->SaveAs(("HDCRPlot/Landau"+to_string((int)index)+".pdf").c_str());
    LandCanv->Close();
  }

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

  Proiezione(h2_l,h2_r,h2_t,nbinx,index);

  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x+[1])",0.8*fit_r->GetParameter(1),1.5*fit_r->GetParameter(1));
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x+[1])",0.5*fit_l->GetParameter(1),1.2*fit_l->GetParameter(1));

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");


    /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  /* hyp_l->SetParameter(3, -2.43e-2);
  */
  hyp_r->SetParameter(0, 8.51);
  hyp_r->SetParameter(1, 5);
  hyp_r->SetParameter(2, 1.2);

 
  wf_c->Divide(3,2);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");

   h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("SAMEP");
  hyp_l->DrawF1(0,0.5,"same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
  hyp_r->DrawF1(0,0.5,"same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","0");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");
  hyp_t->DrawF1(txmin,txmax,"same");

  rymin_lc=rymin_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+hyp_l->GetParameter(0);
  rymax_lc=rymax_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+hyp_l->GetParameter(0);
  rymin_rc=rymin_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+hyp_r->GetParameter(0);
  rymax_rc=rymax_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+hyp_r->GetParameter(0);
  tymin_c=tymin-(hyp_l->Eval(0.25)-hyp_l->GetParameter(0)+hyp_r->Eval(0.25)-hyp_r->GetParameter(0))/2;
  tymax_c=tymax-(hyp_l->Eval(0.25)-hyp_l->GetParameter(0)+hyp_r->Eval(0.25)-hyp_r->GetParameter(0))/2;

  
  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  

  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
        hc_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));
	hc_t->Fill(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2);
      }

  }//chiudo for k

    for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
    histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
    histotemp_t=hc_t->ProjectionY("hc_tprojY",k,k);
   
   
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();

    y_l[k]=histotemp_l->GetMean();
    RMS[0][k]= histotemp_l->GetMeanError();
    
    y_r[k]=histotemp_r->GetMean();
    RMS[1][k]= histotemp_r->GetMeanError();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    
  }//chiudo for k


   TGraphErrors* graph_lc = new TGraphErrors(nbinx-1,x_l,y_l,0,RMS[0]);
   TGraphErrors* graph_rc = new TGraphErrors(nbinx-1,x_r,y_r,0,RMS[1]);
   TGraphErrors* graph_tc = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x",-0.2,0.4);
  
   graph_tc->Fit("fit_tdiff","R0L");
   
   for(k=0;k<digiTree->GetEntries();k++){
     
     
     digiTree->GetEntry(k);
     if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
       {
	 hc_tdiff->Fill(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0))+fit_tdiff->GetParameter(0));
       }
     
   }//chiudo for k
   
  

 
   
   TF1* retta = new TF1("retta","[0]+[1]*x",txmin+0.3,txmax-0.3);
   

   for(k=0;k<nbinx;k++){
     
     
     TH1D* histotemp_t;
     histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
     
   
     yt[k]=histotemp_t->GetMean();
     RMS[2][k]= histotemp_t->GetMeanError();
     
     delete histotemp_t;
     
     
   }//chiudo for k
    
   TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   graph_tcdiff->Fit("retta","0R");
   
   wf_c->cd(4);
   hc_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
   hc_l->GetXaxis()->SetTitle("max.amplitude [mV]");
   hc_l->Draw("COLZ");
   graph_lc->SetMarkerStyle(8);
   graph_lc->SetMarkerSize(.5);
   graph_lc->Draw("P");

   
   wf_c->cd(5);
   
   hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
   hc_r->GetXaxis()->SetTitle("max.amplitude [mV]");
   hc_r->Draw("COLZ");
   graph_rc->SetMarkerStyle(8);
   graph_rc->SetMarkerSize(.5);
   graph_rc->Draw("P");

  
   wf_c->cd(6);

   hc_tdiff->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
   hc_tdiff->GetXaxis()->SetTitle("t_left-t_right [ns]");
   hc_tdiff->Draw("COLZ");
   graph_tcdiff->SetMarkerStyle(8);
   graph_tcdiff->SetMarkerSize(.5);
   graph_tcdiff->Draw("P");
   
   wf_c->SaveAs(("HDCRPlot/Controllo"+to_string((int)index)+".pdf").c_str());
   wf_c->Close();
   /*#################################################################################################################################*/ 
   
   TH1D* histo_cl;
   TH1D* histo_cr;
   TH1D* histo_ct;
   TH1D* histo_ctdiff;
   
   histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
   histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);
   
   histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);
   histo_ctdiff = hc_tdiff->ProjectionY("histo_ctdiff",0,nbinx);

   TF1* gaus_cl = new TF1("gaus_cl","gaus");
   TF1* gaus_cr = new TF1("gaus_cr","gaus");
   TF1* gaus_ct = new TF1("gaus_ct","gaus", rymin_rc, rymax_rc);
   TF1* gaus_ctdiff = new TF1("gaus_ctdiff","gaus", rymin_rc, rymax_rc);

   histo_ct->SetLineColor(kBlack);
   histo_cl->SetLineColor(kBlue);
   histo_cr->SetLineColor(kRed);
   gaus_ct->SetLineColor(kBlack);
      
   TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
   TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
   
   gaus_ct->SetParameter(0,80);
   gaus_ct->SetParameter(1, (rymin_rc+rymax_rc)/2);
   gaus_ct->SetParameter(2,0.03);
   histo_ct->Fit("gaus_ct");
   gaus_ctdiff->SetLineColor(kGreen);
   
   gaus_ctdiff->SetParameter(0,gaus_ct->GetParameter(0));
   gaus_ctdiff->SetParameter(1,gaus_ct->GetParameter(1));
   gaus_ctdiff->SetParameter(2,gaus_ct->GetParameter(2));
   
   histo_ctdiff->Fit("gaus_ctdiff");
   histo_ctdiff->SetLineColor(kGreen);
   histo_ctdiff->Draw("SAME");
   histo_ct->Draw("SAME");
   
   *Resol=gaus_ctdiff->GetParameter(2);
   *errResol=gaus_ctdiff->GetParError(2);

   l2->SetHeader("t_{ave}-t_{MCP} distrib");
   l2->AddEntry(histo_ct,"t_ave-t_MCP");
   l2->AddEntry(gaus_ct,("#sigma="+to_string(gaus_ct->GetParameter(2))).c_str());
   l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
   l2->AddEntry(gaus_ctdiff,("#sigma="+to_string(gaus_ctdiff->GetParameter(2))).c_str());
   l2->Draw();
   tdiff->SaveAs(("HDCRPlot/Gaussiane"+to_string((int)index)+".pdf").c_str());
   tdiff->Close();
     
   bool ConfrontoTdiff=true;
   
   if(ConfrontoTdiff){
     TCanvas* ConfCanv = new TCanvas("confronto","",1200,800);
     gStyle->SetOptFit();
     ConfCanv->Divide(3,1);
     
     ConfCanv->cd(1);
     h2_t->Draw("COLZ");
     graph_t->Draw("P");

     ConfCanv->cd(2);
     hc_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
     hc_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
     hc_t->Draw("COLZ");
     graph_tc->Draw("P");
     fit_tdiff->Draw("SAME");

     ConfCanv->cd(3);
     hc_tdiff->Draw("COLZ");
     graph_tcdiff->Draw("P");
     graph_tcdiff->Fit("retta");
     
     ConfCanv->SaveAs(("HDCRPlot/Confronto"+to_string((int)index)+".pdf").c_str());
     ConfCanv->Close();
   }

   delete h2_r;
   delete h2_l;
   delete h2_t;
   delete hc_r;
   delete hc_l;
   delete hc_t;
   delete hc_tdiff;
   delete hr_amp;
   delete hl_amp;
   delete mcp_amp;
   
}


void RisVsDcr(){
  //Int_t nfiles=19;
  Int_t nfiles=5;
  
  TFile* myfile[nfiles];
  Float_t sigma[nfiles],errsigma[nfiles];
  Int_t i;
  Float_t bias[nfiles],NINOthr[nfiles],DCR[nfiles];
  Float_t biasPl[nfiles],NINOthrPl[nfiles],DCRPl[nfiles];
  
  gSystem->Exec("rm -r -f HDCRPlot");
  gSystem->Exec("mkdir HDCRPlot");

  for(i=0;i<nfiles;i++){
    myfile[i]=TFile::Open(("Pd/DCR10072/"+to_string(i)+".root").c_str());    
  }
  
  TTree* info[nfiles];
  
  for(i=0;i<nfiles;i++){
    info[i]=(TTree*)myfile[i]->Get("info");
    info[i]->SetBranchAddress("SiPMCurrent_bar",&DCR[i]);
    info[i]->SetBranchAddress("NINOthr_bar",&NINOthr[i]);
    info[i]->SetBranchAddress("Vbias_bar",&bias[i]);
  }
gROOT->SetBatch(kTRUE);
  for(i=0;i<nfiles;i++){
    info[i]->GetEntry(1);
    biasPl[i]=bias[i];
    NINOthrPl[i]= NINOthr[i];
    DCRPl[i]=DCR[i];
    Ris(myfile[i],&sigma[i],&errsigma[i],DCRPl[i]);
    
}
gROOT->SetBatch(kFALSE);


  TCanvas* defcanv = new TCanvas("RisvsDCR","",600,400);
  defcanv->SetLogy();
  TGraphErrors* graph = new TGraphErrors(nfiles,DCRPl,sigma,0,errsigma);
  graph->GetXaxis()->SetTitle("DCR [#muA]");
  graph->GetYaxis()->SetTitle("t_{res} [ns]");
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(.8);

   graph->Draw("AP");
  
   defcanv->SaveAs("HDCRPlot/plot.pdf");

  TLatex* tex[nfiles];
  cout<<"HERE"<<endl;
  for(i=0;i<nfiles;i++){
    cout<<"HERE"<<endl;
    tex[i]= new TLatex(DCRPl[i],sigma[i],( (to_string((int)NINOthrPl[i]))+"   "+to_string((int)biasPl[i])).c_str());
    tex[i]->SetTextSize(0.035);
    tex[i]->SetLineWidth(2);
    tex[i]->SetLineColor(i);
    
    
    tex[i]->Draw();
  }
  
  
 
  
 

}
