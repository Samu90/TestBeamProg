//pioooo
//pippo
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
#include <TLine.h>

void plotWF(const char * filename,int i){
  
  TFile *  file= new TFile(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");
 
  Float_t tot_WF[6144],wf_time[6144],digi_time[60];
  Float_t ch[6][1024], time[1024],startime[2];
  int wf_sample,step[6144];
  int j,k,LED300;

  
  WFTree->SetBranchAddress("WF_time",&wf_time);
  WFTree->SetBranchAddress("WF_val",&tot_WF);

  WFTree->SetBranchAddress("WF_samples",&wf_sample);
  digiTree->SetBranchAddress("time",&digi_time);
  digiTree->SetBranchAddress("LED300",&LED300);
  int p;
  p=i;

  
  
  for(k=0;k<WFTree->GetEntries();k++){

    if (i==0) p=k;
    cout<<p<<endl;
    WFTree->GetEntry(p);
    digiTree->GetEntry(p);
    
    for(j=0;j<wf_sample;j++){
      step[j]= j;
      ch[(int)j/1024][j%1024]=tot_WF[j];
      if(j<1024) time[j]= wf_time[j];
    }//chiudo for j

    startime[0]=digi_time[1+LED300];
    startime[1]=digi_time[2+LED300];
    
    TCanvas* wf_c =new TCanvas("wf","Plot wf",1000,650);
    TGraph* graph[6];
    TLine* start[2];
    wf_c->Divide(3,2);

    TLegend* leg[2];
    
    for (j=0;j<6;j++){  
      graph[j] = new TGraph(wf_sample/6,wf_time,ch[j]);
      wf_c->cd(j+1);
      //graph[j]->GetXaxis()->SetTitleSize(10);
      graph[j]->GetXaxis()->SetTitle("time [ns]");
      //graph[j]->GetYaxis()->SetTitleSize(10);
      graph[j]->GetYaxis()->SetTitle("amplitude [mV]");
      graph[j]->Draw();
      if(j==1){
	start[0]=new TLine(startime[0],0,startime[0],500);
	start[0]->SetLineColor(kRed);
	start[0]->Draw("same");
	leg[0]=new TLegend();
	leg[0]->AddEntry(start[0],(to_string(startime[0])).c_str());
	leg[0]->Draw("SAME");
	cout<< startime[0]<<endl;
      }
      if(j==2){
	start[1]=new TLine(startime[1],0,startime[1],500);
	start[1]->SetLineColor(kRed);
	start[1]->Draw("same");
	leg[1]=new TLegend();
	leg[1]->AddEntry(start[1],(to_string(startime[1])).c_str());
	leg[1]->Draw("SAME");
	cout<< startime[1]<<endl;
      }
    }//chiudo for j
    
    
    wf_c->Update();  
    if(i!=0) break;
    gPad->WaitPrimitive();
   
  }//chiudo for k
}

