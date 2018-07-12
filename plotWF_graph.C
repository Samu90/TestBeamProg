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
  Float_t amp_max[54],time[54];
  // Float_t a_r[13551],a_l[13551],deltat[13551];
  int n_timetypes;
  int j,k;
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1000,650);
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",100,0.0,400);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",100,0.0,400);
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);


  
  for(k=0;k<digiTree->GetEntries();k++){
    amp_max[3]=0;
    amp_max[4]=0;
   
    // if (k%100==0) cout<<k<<endl;
    digiTree->GetEntry(k);
    cout << k <<"           " << amp_max[3]<< "         " << amp_max[4]<< endl;
    
    hr_amp->Fill(amp_max[3]);
    hl_amp->Fill(amp_max[4]);
    //    a_r[k]=amp_max[3];
    // a_l[k]=amp_max[4];
    //  deltat[k]= time[3]-time[4];
  }//chiudo for k   
  
    
  cout << "here!"<<endl;
  
  hl_amp->SetLineColor(kRed);
  hr_amp->DrawNormalized();
  hl_amp->DrawNormalized("same");
    // TGraph* graph_r = new TGraph(digiTree->GetEntries(),delta_t,a_r);
    // TGraph* graph_l = new TGraph(digiTree->GetEntries(),delta_t,a_l);
    //    graph_r->Draw();
    // graph_l->Draw("same");
  cout << "here!"<<endl;
}

