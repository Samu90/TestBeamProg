
void plotWF(const char * filename,int i){
  
  TFile *  file= new TFile(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
 
  Float_t tot_WF[6144],wf_time[6144];
  Float_t ch[6][1024], time[1024];
  int wf_sample,step[6144];
  int j,k;

  
  WFTree->SetBranchAddress("WF_time",&wf_time);
  WFTree->SetBranchAddress("WF_val",&tot_WF);
  // WFTree->SetBranchAddress("Instance",&step);
  WFTree->SetBranchAddress("WF_samples",&wf_sample);
  int p;
  p=i;
  
  for(k=0;k<WFTree->GetEntries();k++){

    if (i==0) p=k;

    WFTree->GetEntry(p);
    
    for(j=0;j<wf_sample;j++){
      step[j]= j;
      ch[(int)j/1024][j%1024]=tot_WF[j];
      if(j<1024) time[j]= wf_time[j];
    }//chiudo for j

    TCanvas* wf_c =new TCanvas("wf","Plot wf",1000,650);
    TGraph* graph[6];
    wf_c->Divide(3,2);
    
    for (j=0;j<6;j++){  
      graph[j] = new TGraph(wf_sample/6,wf_time,ch[j]);
      wf_c->cd(j+1);
      graph[j]->Draw();
    }//chiudo for j
    
    
    wf_c->Update();
    if(i!=0) break;
    gPad->WaitPrimitive();
  }//chiudo for k

    
      string c = NULL;
      cin>> c;
      cout<< "-"<< c <<endl;
      if(c=="s") {break;}
      
    
  }//chiudo for k
}

