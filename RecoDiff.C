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
#include <TROOT.h> 

void GaussianAmpWalk(TFile* file,Float_t DCRval,Float_t* NAWRes,Float_t* errNAWRes){
TTree* digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[54], time[54];
  int k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,txmin,txmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  Double_t max;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0;
  rxmax=0.5;

  const Int_t  nbinx=290,nbiny=150;

  int i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-1.1;
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
  
  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1D *mcp_amp =new TH1D("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.04,1);
  TF1 *fit_l = new TF1("f_l","landau",0.04,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED300",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);

     
  max=4096;

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
    

  tymin=mean3-1.3*rms3;
  tymax=mean3+1.3*rms3;

  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
    mcp_amp->Fill(amp_max[0]/max);

  }/*chiudo for */

  LEDi=LED300;
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");

  bool DebugLand=true;
  
  if(DebugLand){
    TCanvas* LandCanv = new TCanvas("mycanvas","",1200,700);
    LandCanv->Divide(2,1);
    
    LandCanv->cd(1)->SetLogy();
    hl_amp->Draw();
    fit_l->Draw("SAME");
    
    LandCanv->cd(2)->SetLogy();
    hr_amp->Draw();
    fit_r->Draw("SAME");
    
    LandCanv->SaveAs(("ControlliRecoDiff/Landaubis"+to_string((int)DCRval)+".pdf").c_str());
    LandCanv->Close();
  }

  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* h2_t= new TH2D("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	if (amp_max[4]/max < 0.35) h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if

  }//chiudo for k

  TH1D* histohampl;
  TH1D* histohampr;
  TF1* gaush_r;
  TF1* gaush_l;

  TCanvas* ProiezioniTemp;
  Int_t npoint=50;
  
  gSystem->Exec(("mkdir ControlliRecoDiff/Gauss"+to_string((int)DCRval)).c_str());
  
  for(i=0;i<npoint;i++){
    x_r[i]=0;
    x_l[i]=0;
    y_l[i]=0;
    y_r[i]=0;
    rmsy_r[i]=0;
    rmsy_l[i]=0;
  }
    
  for(i=0;i<npoint;i++){
    
    histohampl=h2_l->ProjectionY("FettaAmpL",h2_l->FindBin(rxmin+(rxmax-rxmin)/npoint*i),h2_l->FindBin(rxmin+(rxmax-rxmin)/npoint*(i+1)));
    histohampr=h2_r->ProjectionY("FettaAmpR",h2_r->FindBin(rxmin+(rxmax-rxmin)/npoint*i),h2_r->FindBin(rxmin+(rxmax-rxmin)/npoint*(i+1)));
  
    x_r[i]=rxmin+(rxmax-rxmin)/npoint*i+((rxmax-rxmin)/npoint)/2;
    x_l[i]=rxmin+(rxmax-rxmin)/npoint*i+((rxmax-rxmin)/npoint)/2;
    
    gaush_r = new TF1("gaushr","gaus",7.7,9);
    gaush_l = new TF1("gaushl","gaus",7.7,9);
    ProiezioniTemp = new TCanvas("mycanv","",1200,800);
    ProiezioniTemp->Divide(2,1);
    
    ProiezioniTemp->cd(1);
    
    histohampl->Draw();
    histohampl->Fit("gaushl","R");
    y_l[i]=gaush_l->GetParameter(1);
    rmsy_l[i]=gaush_l->GetParError(1);

    ProiezioniTemp->cd(2);
    
    histohampr->Draw();
    histohampr->Fit("gaushr","R");
    y_r[i]=gaush_r->GetParameter(1);
    rmsy_r[i]=gaush_r->GetParError(1);
   
    
    ProiezioniTemp->SaveAs(("ControlliRecoDiff/Gauss"+to_string((int)DCRval)+"/GaussianeAmp"+to_string(i)+".pdf").c_str());
    ProiezioniTemp->Close();
    delete ProiezioniTemp;
    delete gaush_l;
    delete gaush_r;
    //delete mcp_amp;
  
  }//chiudo for
  
  TGraphErrors* graphGAmp_r=new TGraphErrors(npoint,x_r,y_r,0,rmsy_r);
  TGraphErrors* graphGAmp_l=new TGraphErrors(npoint,x_l,y_l,0,rmsy_l);

  TF1* FitLog_r = new TF1("FitLog_r","[0]+[2]*log(x+[1])",0.8*fit_r->GetParameter(1),1.5*fit_r->GetParameter(1));
  TF1* FitLog_l = new TF1("FitLog_l","[0]+[2]*log(x+[1])",0.5*fit_l->GetParameter(1),1.2*fit_l->GetParameter(1));
  
  gStyle->SetOptFit();
  
     /* SetParameters*/
  FitLog_l->SetParameter(0, 8.51);
  FitLog_l->SetParameter(1, 5);
  FitLog_l->SetParameter(2, 1.2);
  
  FitLog_r->SetParameter(0, 8.51);
  FitLog_r->SetParameter(1, 5);
  FitLog_r->SetParameter(2, 1.2);
  

  bool controlFit=false;
  
  if(controlFit)gROOT->SetBatch(kFALSE);
 
  TCanvas* CanvAmpGauss = new TCanvas("canv amplitude gauss","",1200,800);
  CanvAmpGauss->Divide(3,2);
  
  CanvAmpGauss->cd(1);
  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
  h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graphGAmp_l->SetMarkerStyle(8);
  graphGAmp_l->SetMarkerSize(.8);
  graphGAmp_l->Draw("SAMEP");
  graphGAmp_l->Fit("FitLog_l");

  CanvAmpGauss->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graphGAmp_r->SetMarkerStyle(8);
  graphGAmp_r->SetMarkerSize(.8);
  graphGAmp_r->Draw("SAMEP");
  graphGAmp_r->Fit("FitLog_r");
  
  CanvAmpGauss->cd(3);
  h2_t->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_t->Draw("COLZ");
  


  rymin_lc=rymin_l-FitLog_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+FitLog_l->GetParameter(0);
  rymax_lc=rymax_l-FitLog_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+FitLog_l->GetParameter(0);
  rymin_rc=rymin_r-FitLog_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+FitLog_r->GetParameter(0);
  rymax_rc=rymax_r-FitLog_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+FitLog_r->GetParameter(0);
  tymin_c=tymin-(FitLog_l->Eval(0.25)-FitLog_l->GetParameter(0)+FitLog_r->Eval(0.25)-FitLog_r->GetParameter(0))/2;
  tymax_c=tymax-(FitLog_l->Eval(0.25)-FitLog_l->GetParameter(0)+FitLog_r->Eval(0.25)-FitLog_r->GetParameter(0))/2;
  


  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);

  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]-FitLog_l->Eval(amp_max[3]/max)+FitLog_l->GetParameter(0));
        hc_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]-FitLog_r->Eval(amp_max[4]/max)+FitLog_r->GetParameter(0));
	hc_t->Fill(time[1+LEDi]-time[2+LEDi],(time[1+LEDi]+time[2+LEDi])/2-time[0]-(FitLog_r->Eval(amp_max[4]/max)-FitLog_r->GetParameter(0)+FitLog_l->Eval(amp_max[3]/max)-FitLog_l->GetParameter(0))/2);
      }
    
  }//chiudo for k
  gSystem->Exec("mkdir ControlliRecoDiff/NAW");

  gStyle->SetOptFit();
  
  CanvAmpGauss->cd(4);
  hc_l->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  hc_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  hc_l->Draw("COLZ");

  CanvAmpGauss->cd(5);
  hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  hc_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  hc_r->Draw("COLZ");

  CanvAmpGauss->cd(6);
  hc_t->GetYaxis()->SetTitle("t_{ave}-t_MCP [ns]");
  hc_t->GetXaxis()->SetTitle("t_{left}-t_{right} [ns]");
  hc_t->Draw("COLZ");
  
  if(controlFit) gPad->WaitPrimitive();
  CanvAmpGauss->SaveAs(("ControlliRecoDiff/GaussianeAmp"+to_string((int)DCRval)+".pdf").c_str());
  CanvAmpGauss->Close();
  if(controlFit) gROOT->SetBatch(kTRUE);

  
  TH1D* tProjection = hc_t->ProjectionY("ProjectionTave",0,nbinx);
  TCanvas* RisDef = new TCanvas("RisDef","",1200,700);

  TF1* fitgaussiano = new TF1("fitgaussiano","gaus");
  gStyle->SetOptFit(0111);
  fitgaussiano->SetRange(tProjection->GetBinCenter(tProjection->GetMaximumBin())-1, tProjection->GetBinCenter(tProjection->GetMaximumBin())+1);
  
  tProjection->Fit("fitgaussiano","R");
  tProjection->Draw("SAME");
  cout<< "E IL FIT??_________________________________________________________________________________________________________"<<endl;
  RisDef->SaveAs(("ControlliRecoDiff/NewAWTimeRes"+to_string((int)DCRval)+".pdf").c_str());
  RisDef->Close();
  
  *NAWRes=fitgaussiano->GetParameter(2);
  *errNAWRes=fitgaussiano->GetParError(2);

  delete hr_amp;
  delete hl_amp;
  delete h2_l;
  delete h2_r;
  delete h2_t;
  delete hc_r;
  delete hc_l;
  delete hc_t;
  delete RisDef;
  delete fitgaussiano;

}//CHIUDO FUNZIONE




void RecoDiff(){
  
  Int_t nfiles;
  TFile* FNoCorr[nfiles]
  TFile* FWindowCorr[nfiles];
  TFile* FSmartAnalysis[nfiles];
    
  TTree* infoNoCorr[3][nfiles];
    
  Double_t Ris[3][nfiles], errRis[3][nfiles];
  Float_t DCR[3][nfilens];
  Double_t DCRPl[3][nfilens];

  for(i=0;i<nfiles;i++){
    FNoCorr[i]=TFile::Open(("Pd/DCR10072/"+to_string(i)+".root").c_str());
    FWindowCorr[i]=TFile::Open(("Pd/NewRecoDcr/"+to_string(i)+".root").c_str());
    FSmartAnalysis[i]==TFile::Open(("Pd/SmartAnalysisDcr/"+to_string(i)+".root").c_str());;
  }
  

  
  for(i=0;i<nfiles;i++){
    info[0][i]=(TTree*)FNoCorr[i]->Get("info");
    info[1][i]=(TTree*)FWindowCorr[i]->Get("info");
    info[2][i]=(TTree*)FSmartAnalysis[i]->Get("info");
    
    info[0][i]->SetBranchAddress("SiPMCurrent_bar",&DCR[0][i]);
    info[1][i]->SetBranchAddress("SiPMCurrent_bar",&DCR[1][i]);
    info[2][i]->SetBranchAddress("SiPMCurrent_bar",&DCR[2][i]);
    
    DCRPl[0][i]=DCR[0][i];
    DCRPl[1][i]=DCR[1][i];
    DCRPl[2][i]=DCR[2][i];
    
  }
  //(TFile* file,Float_t DCRval,Float_t* NAWRes,Float_t* errNAWRes)
  
  for(i=0;i<nfiles;i++){
    
    GaussianAmpWalk(FNoCorr[i],DCRPl[0][i],&Ris[0][i],&errRis[0][i]);
    GaussianAmpWalk(FWindowCorr[i],DCRPl[1][i],&Ris[1][i],&errRis[1][i]);
    GaussianAmpWalk(FSmartAnalysis[i],DCRPl[2][i],&Ris[2][i],&errRis[2][i]);
  
  }

  
  
  
  
  



}//chiudo main
