

void plotWF_distrib(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");
  Int_t ndat=digiTree->GetEntries();
  Float_t amp_max[54], time[54], x[ndat], y[ndat];
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi,i;
  

  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED300",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);


  
  TH1D* histo1 = new TH1D("myhisto1","Histo_1",200,0,50);
  TH1D* histo2 = new TH1D("myhisto2","Histo_2",500,-120, -60);
  
  
  
  for(i=0;i<digiTree->GetEntries();i++){
    
    digiTree->GetEntry(i);
    
    histo1->Fill(time[2+LED300]-time[0]);

    histo2->Fill(time[2+LED300]-time[5]-time[0]);

    if(time[0]<1000 && time[0]>-1000){x[i]=time[0];}
    else{x[i]=0;}
    
    if(time[0]<1000 && time[0]>-1000) {y[i]=amp_max[0];}
    else{y[i]=0;}
    
    if(i % 10000==0)cout << i << " / " << digiTree->GetEntries()<< endl; 
  }

  TGraph* grafico = new TGraph(ndat,x,y);
  
  TCanvas* CanvHisto = new TCanvas("Mycanvas","",1000,600);
  CanvHisto->Divide(3,1);

  CanvHisto->cd(1);
  histo1->Draw("HISTO");

  CanvHisto->cd(2);
  histo2->Draw("HISTO");

  CanvHisto->cd(3);
  grafico->SetMarkerStyle(8);
  grafico->SetMarkerSize(.2);
  grafico->Draw("AP");
  
}
