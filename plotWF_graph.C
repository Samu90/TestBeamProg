void plotWF_graph(const char * filename){
  
  TFile *  file= new TFile(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");
  Float_t amp_max[5],time[5];
  Float_t a_r[13551],a_l[13551],delta_t[13551];
  int n_timetypes;
  int j,k;

  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);

  for(int i=0;i<digiTree->GetEntries();i++){
    
    a_r[i]=0;
    a_l[i]=0;
    delta_t[i]=0;
  }
  
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%100==0) cout<<k<<endl;
    digiTree->GetEntry(k);
    
    a_r[k]=amp_max[3];
    a_l[k]=amp_max[4];
    //  deltat[k]= time[3]-time[4];
   
    TCanvas* wf_c =new TCanvas("wf","Plot wf",1000,650);
    TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",100,0.0,400);
    hr_amp->Fill(a_r[k]);
    TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",100,0.0,400);
    hl_amp->Fill(a_l[k]);
    hr_amp->DrawNormalized();
    hl_amp->DrawNormalized("same");
    // TGraph* graph_r = new TGraph(digiTree->GetEntries(),delta_t,a_r);
    // TGraph* graph_l = new TGraph(digiTree->GetEntries(),delta_t,a_l);
    //    graph_r->Draw();
    // graph_l->Draw("same");
  }//chiudo for k

}

