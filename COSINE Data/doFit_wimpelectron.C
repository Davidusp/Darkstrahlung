//#include "./src"

void doFit_wimpelectron(int mdm)
{                                                                                
  double cmass[5] = {9.15,9.16,18.0,12.5,12.5};

  TH1F * hdata[5];
  TH1F * hdatafit[5];
  TH1F * hhypo[5];
  double data=0.0;
  double err=0.0;
  TH1F * hbkg[5][33];
  TH1F * hbkgfit[5][33];
  double bkg=0.0;
  TH1F * hunc[5];
  double unc[5][33];
  float rate[5];
  float energy;
  Int_t crystals[5] = {2,3,4,6,7};
  TFile* f_SET2 = new TFile("/home/LuisFranca/wimp_electron/fit/bkgdmodel_SET2.root"); 

  for(int k=0;k<5;k++){
    FILE *s1 = fopen(Form("/home/LuisFranca/wimp_electron/correctedspectra/set3res/energy/COSINE100_med_m%i_heavy_vector.txt",mdm),"r"); 
    hdata[k] = (TH1F*)f_SET2->Get(Form("hData_CRY%d_CH1", crystals[k]));
    hdatafit[k]=new TH1F(Form("hdatafit%d",k+1),"",16,1.0,5.0);
    hhypo[k]=new TH1F("hhypo","",16,1.0,5.0);
    for(int bin=0;bin<16;bin++){
      fscanf(s1,"%f %f %f %f %f %f",&energy,&rate[0],&rate[1],&rate[2],&rate[3],&rate[4]);
      data=hdata[k]->GetBinContent(bin+5);
      err=hdata[k]->GetBinError(bin+5);
      hdatafit[k]->SetBinContent(bin+1,data);
      hdatafit[k]->SetBinError(bin+1,err);
      if (k==0) hhypo[k]->SetBinContent(bin+1,rate[k]*1.5337*365.25*cmass[k]*0.25);
      if (k==1) hhypo[k]->SetBinContent(bin+1,rate[k]*1.5921*365.25*cmass[k]*0.25);
      if (k==2) hhypo[k]->SetBinContent(bin+1,rate[k]*1.5921*365.25*cmass[k]*0.25);
      if (k==3) hhypo[k]->SetBinContent(bin+1,rate[k]*1.5921*365.25*cmass[k]*0.25);
      if (k==4) hhypo[k]->SetBinContent(bin+1,rate[k]*1.3427*365.25*cmass[k]*0.25);
      cout<<k<<":"<<energy<<":"<<rate[k]<<endl;
    }
  }

  for(int i=0;i<5;i++){
    int counter=0;
    for(int j=0; j<34; j++){
      if(j==12 or j==16) continue;
      hbkg[i][counter] = (TH1F*)f_SET2->Get(Form("hBkgd_CRY%d_CH1_COMPONENT%d", crystals[i],j+1));
      hbkgfit[i][counter]=new TH1F(Form("hbkgfit%d%d",i+1,j+1),"",16,1.0,5.0);
      for(int bin=4;bin<20;bin++){
        bkg=hbkg[i][counter]->GetBinContent(bin+1);
        hbkgfit[i][counter]->SetBinContent(bin-3,bkg);
      }
      counter++;
    }
    hbkg[i][32] = (TH1F*)f_SET2->Get(Form("hBkgd_CRY%d_CH1_COMPONENT%d", crystals[i],45));
    hbkgfit[i][32]=new TH1F(Form("hbkgfit%d%d",i+1,45),"",16,1.0,5.0);
    for(int bin=4;bin<20;bin++){
      bkg=hbkg[i][32]->GetBinContent(bin+1);
      hbkgfit[i][32]->SetBinContent(bin-3,bkg);
    }  
  }

  for(int i = 0;i<5;i++){
    int counter=0;
    hunc[i] = (TH1F*)f_SET2->Get(Form("hUncertainty_CRY%d", crystals[i]));
    for(int j=0; j<34; j++){
      if(j==12 or j==16) continue;
      unc[i][counter]=hunc[i]->GetBinContent(j+1);
      counter++;
    }
    unc[i][32]=hunc[i]->GetBinContent(45);
    cout<<counter<<":"<<unc[0][30]<<endl;
  }

  TFile fout(Form("/data/COSINE/WORK/LuisFranca/wimp_electron/fit/output/out_wimpelectron_SET2_%i.root",mdm), "recreate");
  
  CAT * fitter = new CAT();
  for(int i = 0; i < 5; i++){
    fitter->SetData(hdatafit[i], i);
    fitter->SetHypo(hhypo[i], i, 1.0, 0.001);
    for(int j=0;j<33;j++){
      fitter->SetBkgd(hbkgfit[i][j],i,j,1.0,unc[i][j],1);
    }
  }

  fitter->SetNGeneration(5000000);
  //fitter->SetStepSize(0.01);
  fitter->SetStepSizeVariationOn();
  fitter->SetMarginalOn();
  
  cout << "Doing Fit" << endl;
  fitter->DoFit();

  TTree * tree = (TTree*)fitter->GetOutput();
  tree->Write();
  fout.Close();
}