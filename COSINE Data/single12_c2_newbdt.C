#include "TCut.h"
#include "TCanvas.h"
#include "TTree.h"
#include <fstream>
#include <iostream>

void single12_c2_newbdt(int bin, int cnum){

double loli [10] = {0,1,2,3,4,5,6,7,8,9};
double hili [10] = {0,2,3,4,5,6,7,8,9,10};

TString indir = Form("/data/COSINE/WORK/LuisFranca/mod/c%d/",cnum);
TChain *ch = new TChain("smrgd");
int nFile;

if (bin>0){nFile = ch->Add(indir+"trim_c2*");}

double q0[1+8] = {0,0,7778.244,8021.927,7421.372,0,7357.201,7485.456,0};

float p0new[1+8] = {0,0,1.52e-06,1.50e-07,1.50e-09,0,6.01e-08,6.01e-08,0};
float p1new[1+8] = {0,0,-107.840,-120.000,-110.000,0,-120.177,-120.177,0};
float p2new[1+8] = {0,0,   -1.94,    0.0,   0.3,0,   0.8,   -0.24,0};
float p3new[1+8] = {0,0,  -29.288,   -20.0,  -10.00,0, -21.615,  -31.615,0};

float line1[8], line2[8], line3[8];
float deff1[8], deff2[8];
float li1[76], li2[76], li3[76], li4[76];

FILE *litime = fopen("/home/LuisFranca/c3/timelimit.txt","r");
FILE *myline = fopen("/home/LuisFranca/c3/linecut.txt","r");

for (int i=0; i<76; i++){
  fscanf(litime,"%f %f %f %f",&li1[i],&li2[i],&li3[i],&li4[i]);
  }

for (int i=0; i<8; i++){
  fscanf(myline,"%f %f %f",&line1[i],&line2[i],&line3[i]);
  }

int count, countcheck = 0.;
int norm = 1451606400.;
int day = 86400.;
int binsize = 432;
float min = 0.0001;

char crate1[2000], crate2[2000], crate3[2000], crate4[2000], crate5[2000], crate6[2000], crate7[2000], crate8[2000], crate9[2000];
Double_t rate1[2000], rate2[2000], rate3[2000], rate4[2000], rate5[2000], rate6[2000], rate7[2000], rate8[2000], rate9[2000];

long eventtime, iEvtSec, fEvtSec, subDuration;
ch->SetBranchAddress("eventtime", &eventtime);
ch->SetBranchAddress("iEvtSec", &iEvtSec);
ch->SetBranchAddress("fEvtSec", &fEvtSec);
ch->SetBranchAddress("subDuration", &subDuration);

TH1F *hTime[2000];

int total = ch->GetEntries();
ch->GetEntry(0);
long iTime = (iEvtSec-norm)/86400.;
ch->GetEntry(total-1);
long fTime = (fEvtSec-norm)/86400.;
const int nBin = (int)(fTime-iTime);
Int_t timeBin = (fTime-iTime)*2.;
//cout<<"==========================="<<endl;
Double_t iTimeBin[2000], fTimeBin[2000];

FILE *fp = fopen(Form("/home/LuisFranca/c3/log_axion/c2/ls_80/test_single_c%d_newbdt_axion.txt",cnum),"a");

for (int i=li1[bin]; i<li1[bin]+1; i++){
  iTimeBin[i] = li3[bin];
  fTimeBin[i] = li4[bin];
  hTime[i] = new TH1F(Form("hTime_%1.f",li1[bin]),"",binsize,li2[i],li2[i]+18);

  //(New BDT)
  sprintf(crate1,"((eventtime-%d)/%d)>%f && ((eventtime-%d)/%d)<%f && (cEnergy_%d>=1 && cEnergy_%d<=2) && pmtnc1_%d>1 && pmtnc2_%d>1 && cxrqcn_%d>-1 && ((%f*TMath::Exp(%f*bdt_%d)+%f-(%f*bdt_%d)) < (cEnergy_%d)) && MuontotaldeltaT0/1e6>30 && cxt1_%d>2.2 && (1-(cxnx2_%d-cxnx1_%d))/2>%f && ((1-(cxnx2_%d-cxnx1_%d))/2-%f*bdt_%d < %f) && LScharge/143.8<=80 && LSCoin==1",
norm,day,iTimeBin[i],norm,day,fTimeBin[i],cnum,cnum,cnum,cnum,cnum,p0new[cnum],p1new[cnum],cnum,p2new[cnum],p3new[cnum],cnum,cnum,cnum,cnum,cnum,line3[cnum-1],cnum,cnum,line1[cnum-1],cnum,line2[cnum-1]); 
  
  sprintf(crate2,"((eventtime-%d)/%d)>%f && ((eventtime-%d)/%d)<%f && (cEnergy_%d>=2 && cEnergy_%d<=3) && pmtnc1_%d>1 && pmtnc2_%d>1 && cxrqcn_%d>-1 && ((%f*TMath::Exp(%f*bdt_%d)+%f-(%f*bdt_%d)) < (cEnergy_%d)) && MuontotaldeltaT0/1e6>30 && cxt1_%d>2.2 && (1-(cxnx2_%d-cxnx1_%d))/2>%f && ((1-(cxnx2_%d-cxnx1_%d))/2-%f*bdt_%d < %f) && LScharge/143.8<=80 && LSCoin==1",
norm,day,iTimeBin[i],norm,day,fTimeBin[i],cnum,cnum,cnum,cnum,cnum,p0new[cnum],p1new[cnum],cnum,p2new[cnum],p3new[cnum],cnum,cnum,cnum,cnum,cnum,line3[cnum-1],cnum,cnum,line1[cnum-1],cnum,line2[cnum-1]);

  sprintf(crate3,"((eventtime-%d)/%d)>%f && ((eventtime-%d)/%d)<%f && (cEnergy_%d>=3 && cEnergy_%d<=4) && pmtnc1_%d>1 && pmtnc2_%d>1 && cxrqcn_%d>-1 && ((%f*TMath::Exp(%f*bdt_%d)+%f-(%f*bdt_%d)) < (cEnergy_%d)) && MuontotaldeltaT0/1e6>30 && cxt1_%d>2.2 && (1-(cxnx2_%d-cxnx1_%d))/2>%f && ((1-(cxnx2_%d-cxnx1_%d))/2-%f*bdt_%d < %f) && LScharge/143.8<=80 && LSCoin==1",
norm,day,iTimeBin[i],norm,day,fTimeBin[i],cnum,cnum,cnum,cnum,cnum,p0new[cnum],p1new[cnum],cnum,p2new[cnum],p3new[cnum],cnum,cnum,cnum,cnum,cnum,line3[cnum-1],cnum,cnum,line1[cnum-1],cnum,line2[cnum-1]);

  sprintf(crate4,"((eventtime-%d)/%d)>%f && ((eventtime-%d)/%d)<%f && (cEnergy_%d>=4 && cEnergy_%d<=5) && pmtnc1_%d>1 && pmtnc2_%d>1 && cxrqcn_%d>-1 && ((%f*TMath::Exp(%f*bdt_%d)+%f-(%f*bdt_%d)) < (cEnergy_%d)) && MuontotaldeltaT0/1e6>30 && cxt1_%d>2.2 && (1-(cxnx2_%d-cxnx1_%d))/2>%f && ((1-(cxnx2_%d-cxnx1_%d))/2-%f*bdt_%d < %f) && LScharge/143.8<=80 && LSCoin==1",
norm,day,iTimeBin[i],norm,day,fTimeBin[i],cnum,cnum,cnum,cnum,cnum,p0new[cnum],p1new[cnum],cnum,p2new[cnum],p3new[cnum],cnum,cnum,cnum,cnum,cnum,line3[cnum-1],cnum,cnum,line1[cnum-1],cnum,line2[cnum-1]);

  sprintf(crate5,"((eventtime-%d)/%d)>%f && ((eventtime-%d)/%d)<%f && (cEnergy_%d>=5 && cEnergy_%d<=6) && pmtnc1_%d>1 && pmtnc2_%d>1 && cxrqcn_%d>-1 && ((%f*TMath::Exp(%f*bdt_%d)+%f-(%f*bdt_%d)) < (cEnergy_%d)) && MuontotaldeltaT0/1e6>30 && cxt1_%d>2.2 && (1-(cxnx2_%d-cxnx1_%d))/2>%f && ((1-(cxnx2_%d-cxnx1_%d))/2-%f*bdt_%d < %f) && LScharge/143.8<=80 && LSCoin==1",
norm,day,iTimeBin[i],norm,day,fTimeBin[i],cnum,cnum,cnum,cnum,cnum,p0new[cnum],p1new[cnum],cnum,p2new[cnum],p3new[cnum],cnum,cnum,cnum,cnum,cnum,line3[cnum-1],cnum,cnum,line1[cnum-1],cnum,line2[cnum-1]);

  rate1[i] = ch->GetEntries(crate1);
  rate2[i] = ch->GetEntries(crate2);
  rate3[i] = ch->GetEntries(crate3);
  rate4[i] = ch->GetEntries(crate4);
  rate5[i] = ch->GetEntries(crate5);
  cout<<Form("%1.f %.3f %1.f %1.f %1.f %1.f %1.f %1.f %1.f %1.f %1.f\n",li1[bin],li3[bin]+(li4[bin]-li3[bin])/2,rate1[i],rate2[i],rate3[i],rate4[i],rate5[i],rate6[i],rate7[i],rate8[i],rate9[i])<<endl;
  fprintf(fp,"%1.f %.3f %1.f %1.f %1.f %1.f %1.f %1.f %1.f %1.f %1.f\n",li1[bin],li3[bin]+(li4[bin]-li3[bin])/2,rate1[i],rate2[i],rate3[i],rate4[i],rate5[i],rate6[i],rate7[i],rate8[i],rate9[i]);
  }

cout<<"***** Finished *****"<<endl;
fclose(fp);
}
