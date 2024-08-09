#include "smalltree.h"
#include <climits>
#include "TH1D.h"
#include <string>
#include <iostream>
#include <fstream>

void crystalpreselection(int run, int sub){
TChain *intree = new TChain("ntp");
//intree->Add(Form("/data/COSINE/WORK/hafizh/COSINE/MERGED/phys/V00-04-04/mrgd_M%06d.root.%03d",run,sub));

// ***** DAMA STUDY VERSION 00.04.15
//TString indir  = "/data/COSINE/WORK/hafizh/COSINE/NTP/phys/V00-04-15/";
//TString indir  = "/data/COSINE/ONLINE/MRGD/V00-04-15/";

TString dataname   = Form("mrgd_M00%04d.root.%03d",run,sub);
TString outputname = Form("sntp_%04d.root.%03d",run,sub);

TString indir  = Form("/data/COSINE/MRGD/V00-04-20_p04/%04d/",run);
//TString indir  = Form("/data/COSINE/WORK/hafizh/COSINE/MERGED/phys/V00-04-%02d/",version);	// hafizh directory

// ***** OUTPUT DATA
TString output = "/data/COSINE/WORK/dfreitas/Cdata_less4MeV/";
//TString output = Form("/data/COSINE/WORK/hafizh/COSINE/MERGED/trim_phys/V00-04-%02d/c%d/",version,i+1);
cout<<dataname<<endl;
intree->Add( (indir+dataname).Data() );
TFile fout( (output+outputname).Data(),"recreate");

//=====================================================
if (intree) cout<<"***** ntuple opened *****"<<endl;

Int_t nevent = intree->GetEntries();
cout<<"Total events : "<<nevent<<endl;

float q0[8], eff0[8], eff1[8];
int i;
//float p0 = 1.79006e+06;
//float p1 = -123.731;
//float p2 = 2.61407;
//float p3 = -11.9022;

Int_t eventnumber;
Long64_t eventtime;
Long64_t fadctrgtime;
Long64_t MuondeltaT0;
Long64_t MuondeltaTrgTime0;
Long64_t MuontotaldeltaT0;
Long64_t MuondiffDeltaC[1+8];
Long64_t MuondiffTrgTimeC[1+8];
Long64_t MuondiffTotalC[1+8];
Long64_t LStime;

Int_t pmtnc1[1+8], pmtnc2[1+8];
Float_t pmtqc51[1+8], pmtqc52[1+8];
Float_t pmtmt1[1+8], pmtmt2[1+8];
Float_t pmt1t0[1+8], pmt2t0[1+8];
Float_t trgbit1[1+8], trgbit2[1+8];

Int_t cxnc[1+8];
Float_t cxqc[1+8];
Float_t cxt0[1+8], cxt1[1+8], cxx1[1+8], cxx2[1+8];
Float_t cxnx1[1+8], cxnx2[1+8], cxqc5[1+8];
Float_t cxmt[1+8], cEnergy[1+8], cEnergyD[1+8];
Float_t cxnq600[1+8], cxrqcn[1+8], cxtzero[1+8];

Float_t bdtA[1+8];
Float_t bdt[1+8];
Float_t LScharge;
Int_t LSCoin;

Float_t pmt11rqtD1_5, pmt21rqtD1_5, pmt31rqtD1_5, pmt41rqtD1_5, pmt51rqtD1_5, pmt61rqtD1_5, pmt71rqtD1_5,pmt81rqtD1_5;
Float_t pmt12rqtD1_5, pmt22rqtD1_5, pmt32rqtD1_5, pmt42rqtD1_5, pmt52rqtD1_5, pmt62rqtD1_5, pmt72rqtD1_5,pmt82rqtD1_5;
TTree *outtree = new TTree("smrgd","smrgd");

intree->GetLeaf("bdt")->                        SetAddress(bdt);
intree->GetLeaf("bdtA")->                       SetAddress(bdtA);
intree->GetLeaf("BLSVeto","Charge")->           SetAddress(&LScharge);
intree->GetLeaf("BLSVeto","isCoincident")->     SetAddress(&LSCoin);
intree->GetLeaf("BLSVeto","Time")->             SetAddress(&LStime);
intree->GetLeaf("eventsec")->                   SetAddress(&eventtime);
intree->GetLeaf("eventNumber")->                SetAddress(&eventnumber);
intree->GetLeaf("trgtime")->                    SetAddress(&fadctrgtime);
intree->GetLeaf("BMuon","deltaT0")->            SetAddress(&MuondeltaT0);
intree->GetLeaf("BMuon","deltaTrgTime0")->      SetAddress(&MuondeltaTrgTime0);
intree->GetLeaf("BMuon","totalDeltaT0")->       SetAddress(&MuontotaldeltaT0);

for (int i=0; i<8; i++){
intree->GetLeaf(Form("crystal%d",i+1),"nc")->		SetAddress(&cxnc[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"energy")->       SetAddress(&cEnergy[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"energyD")->      SetAddress(&cEnergyD[i+1]);
//intree->GetLeaf(Form("crystal%d",i+1),"nq600")->        SetAddress(&cxnq600[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"nmt500")->       SetAddress(&cxmt[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"rqcn")->         SetAddress(&cxrqcn[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"x1")->           SetAddress(&cxx1[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"x2")->           SetAddress(&cxx2[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"nx1")->          SetAddress(&cxnx1[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"nx2")->          SetAddress(&cxnx2[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"t0")->           SetAddress(&cxt0[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"t1")->           SetAddress(&cxt1[i+1]);
intree->GetLeaf(Form("crystal%d",i+1),"tzero")->        SetAddress(&cxtzero[i+1]);
intree->GetLeaf(Form("pmt%d1",i+1),"nc")->              SetAddress(&pmtnc1[i+1]);
intree->GetLeaf(Form("pmt%d2",i+1),"nc")->              SetAddress(&pmtnc2[i+1]);
intree->GetLeaf(Form("pmt%d1",i+1),"trgbit")->          SetAddress(&trgbit1[i+1]);
intree->GetLeaf(Form("pmt%d2",i+1),"trgbit")->          SetAddress(&trgbit2[i+1]);
//intree->GetLeaf("BMuonDiff",Form("deltaC%d",i+1))->     	SetAddress(&MuondiffDeltaC[i+1]);
//intree->GetLeaf("BMuonDiff",Form("deltaTrgTimeC%d",i+1))->	SetAddress(&MuondiffTrgTimeC[i+1]);
//intree->GetLeaf("BMuonDiff",Form("totalDeltaC%d",i+1))->	SetAddress(&MuondiffTotalC[i+1]);

outtree->Branch(Form("cEnergy_%d",i+1),           &cEnergy[i+1]);
outtree->Branch(Form("cEnergyD_%d",i+1),          &cEnergyD[i+1]);
outtree->Branch(Form("cxnq600_%d",i+1),             &cxnq600[i+1]);
outtree->Branch(Form("cxmt_%d",i+1),                &cxmt[i+1]);
outtree->Branch(Form("cxrqcn_%d",i+1),              &cxrqcn[i+1]);
//outtree->Branch(Form("cxx1_%d",i+1),                &cxx1[i+1l]);
//outtree->Branch(Form("cxx2_%d",i+1),                &cxx2[i+1]);
outtree->Branch(Form("cxnx1_%d",i+1),               &cxnx1[i+1]);
outtree->Branch(Form("cxnx2_%d",i+1),               &cxnx2[i+1]);
outtree->Branch(Form("cxt0_%d",i+1),              &cxt0[i+1]);
outtree->Branch(Form("cxt1_%d",i+1),              &cxt1[i+1]);
outtree->Branch(Form("cxtzero_%d",i+1),           &cxtzero[i+1]);
outtree->Branch(Form("bdt_%d",i+1),               &bdt[i+1]);
outtree->Branch(Form("bdtA_%d",i+1),              &bdtA[i+1]);
outtree->Branch(Form("trgbit1_%d",i+1),              &trgbit1[i+1]);
outtree->Branch(Form("trgbit2_%d",i+1),              &trgbit2[i+1]);
outtree->Branch(Form("pmtnc1_%d",i+1),              &pmtnc1[i+1]);
outtree->Branch(Form("pmtnc2_%d",i+1),              &pmtnc2[i+1]);
//outtree->Branch(Form("MuondiffDeltaC_%d",i+1),      &MuondiffDeltaC[i+1]);
//outtree->Branch(Form("MuondiffTrgTimeC_%d",i+1),    &MuondiffTrgTimeC[i+1]);
outtree->Branch(Form("MuondiffTotalC_%d",i+1),      &MuondiffTotalC[i+1]);
outtree->Branch("LScharge", &LScharge);
outtree->Branch("LSCoin", &LSCoin);
//outtree->Branch("LStime", &LStime);
outtree->Branch("eventtime", &eventtime);
outtree->Branch("eventnumber", &eventnumber);
outtree->Branch("fadctrgtime", &fadctrgtime);
//outtree->Branch("MuondeltaT0", &MuondeltaT0);
//outtree->Branch("MuondeltaTrgTime0", &MuondeltaTrgTime0);
outtree->Branch("MuontotaldeltaT0", &MuontotaldeltaT0);
outtree->Branch(Form("cxnc_%d",i+1),            &cxnc[i+1]);
//outtree->Branch(Form("cxqc_%d",i+1),            &cxqc[i+1]);
//outtree->Branch(Form("cEnergy_%d",i+1),         &cEnergy[i+1]);
//outtree->Branch(Form("cEnergyD_%d",i+1),        &cEnergyD[i+1]);
//outtree->Branch(Form("cxnq600_%d",i+1),         &cxnq600[i+1]);
//outtree->Branch(Form("cxrqcn_%d",i+1),          &cxrqcn[i+1]);
//outtree->Branch(Form("cxmt_%d",i+1),            &cxmt[i+1]);
//outtree->Branch(Form("trgbit1_%d",i+1),         &trgbit1[i+1]);
//outtree->Branch(Form("trgbit2_%d",i+1),         &trgbit2[i+1]);
//outtree->Branch(Form("pmtnc1_%d",i+1),          &pmtnc1[i+1]);
//outtree->Branch(Form("pmtnc2_%d",i+1),          &pmtnc2[i+1]);
//outtree->Branch(Form("pmtqc51_%d",i+1),         &pmtqc51[i+1]);
//outtree->Branch(Form("pmtqc52_%d",i+1),         &pmtqc52[i+1]);
//outtree->Branch(Form("pmt1t0_%d",i+1),          &pmt1t0[i+1]);
//outtree->Branch(Form("pmt2t0_%d",i+1),          &pmt2t0[i+1]);
//outtree->Branch(Form("pmtmt1_%d",i+1),          &pmtmt1[i+1]);
//outtree->Branch(Form("pmtmt2_%d",i+1),          &pmtmt2[i+1]);
//outtree->Branch(Form("cxnx1_%d",i+1),           &cxnx1[i+1]);
//outtree->Branch(Form("cxx1_%d",i+1),            &cxx1[i+1]);
//outtree->Branch(Form("cxnx2_%d",i+1),           &cxnx2[i+1]);
//outtree->Branch(Form("cxx2_%d",i+1),            &cxx2[i+1]);
//outtree->Branch(Form("cxqc5_%d",i+1),           &cxqc5[i+1]);
//outtree->Branch(Form("cxt0_%d",i+1),            &cxt0[i+1]);
//outtree->Branch(Form("cxt1_%d",i+1),            &cxt1[i+1]);
//outtree->Branch(Form("cxtzero_%d",i+1),         &cxtzero[i+1]);
//outtree->Branch(Form("bdt_%d",i+1),             &bdt[i+1]);
//outtree->Branch(Form("bdtA_%d",i+1),            &bdtA[i+1]);
}

Long64_t iEvtSec,fEvtSec, subDuration;

outtree->Branch("iEvtSec", &iEvtSec, "iEvtSec/L");
outtree->Branch("fEvtSec", &fEvtSec, "fEvtSec/L");
outtree->Branch("subDuration", &subDuration, "subDuration/L");

intree->GetEntry(0);
iEvtSec = eventtime;
intree->GetEntry(nevent-1);
fEvtSec = eventtime;
subDuration = fEvtSec - iEvtSec;
int count=0;
//DAMAES[i] = Form("TMath::Exp(%f + (%f*(1-(cxnx2_%d-cxnx1_%d))/2)) > (cxnq600_%d-0.0001)/%f",deff1[cnum-1],deff2[cnum-1],cnum,cnum,cnum,q0[cnum-1]);

for (int i=0; i<nevent; i++){
 intree->GetEntry(i);
 if (i%20000 ==0) cout<<i<<"th events passed"<<endl;
 //if (trgbit1[i+1]>0 && trgbit2[i+1]>0){
 if ((cEnergy[1]<4000 && cEnergy[1]>1.0) || (cEnergy[2]<4000 && cEnergy[2]>1.0) || (cEnergy[3]<4000 && cEnergy[3]>1.0) 
 || (cEnergy[4]<4000 && cEnergy[4]>1.0) || (cEnergy[5]<4000 && cEnergy[5]>1.0) || (cEnergy[6]<4000 && cEnergy[6]>1.0)
 || (cEnergy[7]<4000 && cEnergy[7]>1.0) || (cEnergy[8]<4000 && cEnergy[8]>1.0)){
 count++; 
 outtree->Fill();
//}
}
}
cout<<"Selected events : "<<count<<endl;
cout<<"***** ntuple closed *****"<<endl;
delete intree;
outtree->Write();
//}

}