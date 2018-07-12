#include <iostream>
#include <TH1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

using namespace std;
void plotWF_graph(const char * filename){
  
  
  TFile *  file= TFile::Open(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[54];
  int k;
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1000,650);
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",100,0.0,1000);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",100,0.0,1000);

  digiTree->SetBranchAddress("amp_max",&amp_max);
  


  
  for(k=0;k<digiTree->GetEntries();k++){
   
    if (k%100==0) cout<<k<<endl;
    digiTree->GetEntry(k);
    
    
    hr_amp->Fill(amp_max[3]);
    hl_amp->Fill(amp_max[4]);

  }//chiudo for k   
  
    
  cout << "here!"<<endl;
  
  hl_amp->SetLineColor(kRed);
  hr_amp->DrawNormalized();
  hl_amp->DrawNormalized("same");
  cout << "here!"<<endl;
}

