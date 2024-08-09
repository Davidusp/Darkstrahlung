//#include "./src"
#include "Rtypes.h"
#include "TMathBase.h"
#include "TError.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <float.h>

void doFit_DPbremss_SET3(long int mdm, long int kappa)
{                                                                                
  double cmass[5] = {9.15,9.16,18.0,12.5,12.5};

  TH1F * hdatafit[5];
  TH1F * hhypo[5];
  float data[5];
  TH1F * hbkgfit[5][10];
  TH1F * htotbkg[5];
  double bkg[10];
  float unc[5][10];
  double rate[5];
  double xpe;
  Int_t crystals[5] = {2,3,4,6,7};
  FILE *datafile = fopen("/data/COSINE/WORK/dfreitas/single_hit/SET3_data_100keV.txt","r"); 
  FILE *s1 = fopen(Form("/data/COSINE/WORK/dfreitas/single_hit/FinalSpec100keV/Specorr100keV_m=%ld_k=%ld.txt",mdm,kappa),"r");
  printf("Original kappa: %ld\n", kappa); // Adjust .10lf for desired precision

  hdatafit[0]=new TH1F("hdatafit1","",739,8.0,1486.0);
  hdatafit[1]=new TH1F("hdatafit2","",739,8.0,1486.0);
  hdatafit[2]=new TH1F("hdatafit3","",739,8.0,1486.0);
  hdatafit[3]=new TH1F("hdatafit4","",739,8.0,1486.0);
  hdatafit[4]=new TH1F("hdatafit5","",739,8.0,1486.0);
  hhypo[0]=new TH1F("hhypo1","",739,8.0,1486.0);
  hhypo[1]=new TH1F("hhypo2","",739,8.0,1486.0);
  hhypo[2]=new TH1F("hhypo3","",739,8.0,1486.0);
  hhypo[3]=new TH1F("hhypo4","",739,8.0,1486.0);
  hhypo[4]=new TH1F("hhypo5","",739,8.0,1486.0);
  for(int bin=0;bin<739;bin++){
    fscanf(datafile,"%f %f %f %f %f",&data[0],&data[1],&data[2],&data[3],&data[4]);
    //cout<<s1<<endl;
    hdatafit[0]->SetBinContent(bin+1,data[0]);
    hdatafit[1]->SetBinContent(bin+1,data[1]);
    hdatafit[2]->SetBinContent(bin+1,data[2]);
    hdatafit[3]->SetBinContent(bin+1,data[3]);
    hdatafit[4]->SetBinContent(bin+1,data[4]);
    //printf("%.8lf",mdm);
    fscanf(s1,"%lf %lf %lf %lf %lf %lf",&xpe,&rate[0],&rate[1],&rate[2],&rate[3],&rate[4]);
    cout<<xpe<<endl;
    hhypo[0]->SetBinContent(bin+1,rate[0]*1033.91*cmass[0]*1E+10);
    hhypo[1]->SetBinContent(bin+1,rate[1]*1055.91*cmass[1]*1E+10);
    hhypo[2]->SetBinContent(bin+1,rate[2]*1055.91*cmass[2]*1E+10);
    hhypo[3]->SetBinContent(bin+1,rate[3]*1055.91*cmass[3]*1E+10);
    hhypo[4]->SetBinContent(bin+1,rate[4]*982.04*cmass[4]*1E+10);
  }

  for(int i=0;i<5;i++){
    FILE *bkgfile = fopen(Form("/data/COSINE/WORK/dfreitas/single_hit/Bkg_data/SET3_bkg_C%d_8xpe_poisson_nonpr.txt",crystals[i]),"r");
    for(int j=0; j<10; j++){
      hbkgfit[i][j]=new TH1F(Form("hbkgfit%d%d",(i+1),(j+1)),"",739,8.0,1486.0);
    }
    for(int bin=0;bin<739;bin++){
      fscanf(bkgfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&bkg[0],&bkg[1],&bkg[2],&bkg[3],&bkg[4],&bkg[5],&bkg[6],&bkg[7],&bkg[8],&bkg[9]);       
      for(int j=0; j<10; j++){
        hbkgfit[i][j]->SetBinContent(bin+1,bkg[j]);
      }
    }
  }

 double errC2[10]={0.00759902,0.00333234,2.0e-01,2.0,2.434e-01,8.9055e-01,0.04397709,0.06929293,0.06971072,0.36164221};
 double errC3[10]={0.01244558,0.00835733,1.6e-01,1.198e-01,2.1428e-01,3.1356e-01,0.03213859,0.02322728,0.08347595,0.35584793};
 double errC4[10]={0.0089812,0.00484939,1.539e-01,7.597e-02,1.62398e-01,4.0116e-01,0.05122261,0.02955826,0.06465133,0.04540802};
 double errC6[10]={0.01804092,0.00296248,3.637e-01,4.6737e-01,2.15395e-01,6.3348e-02,0.02364468,0.03704874,0.06268287,0.01492749};
 double errC7[10]={0.01764252,0.00329033,4.444e-01,3.3263e-01,2.64362e-01,8.4337e-02,0.0222003,0.03490571,0.13615742,0.05917625};

 for(int j=0; j<10; j++){
   unc[0][j]=errC2[j];
   unc[1][j]=errC3[j];
   unc[2][j]=errC4[j];
   unc[3][j]=errC6[j];
   unc[4][j]=errC7[j];    
 }

  TFile fout(Form("/data/COSINE/WORK/dfreitas/single_hit/DataFit/out_DP_SET3_bremss_m=%ld_k=%ld.root",mdm,kappa), "recreate");
  
  CAT * fitter = new CAT();
  for(int i=0;i<5;i++){
    fitter->SetData(hdatafit[i], i);
    fitter->SetHypo(hhypo[i], i, 1.0, 1e-1);
    for(int j=0;j<10;j++){
      fitter->SetBkgd(hbkgfit[i][j],i,j,1.0,unc[i][j],1);
    }
  }

  fitter->SetNGeneration(1000000);
  //fitter->SetStepSize(0.01);
  fitter->SetStepSizeVariationOn();
  fitter->SetMarginalOn();
  
  cout << "Doing Fit" << endl;
  fitter->DoFit();

  TTree * tree = (TTree*)fitter->GetOutput();
  tree->Write();
  //hdatafit[0]->Write();
  //hhypo[0]->Write();
  fout.Close();
}