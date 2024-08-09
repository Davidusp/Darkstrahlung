//#include "./src"

void extractfit()
{     

  int mDP[43] = {10,20,30,40,50,70,100,130,170,220,290,390,520,690,910,1210,1600,2120,2810,3730,4940,6550,8690,11510,
                15260,11514,20240,26830,35560,47150,62510,82860,109850,145630,193070,255950,339320,449840,596360,790600,
                1048110,1389500,1842070};   // DP mass

  long int kappa[49] = {10,15,24,37,56,87,133,205,316,487,750,1155,1778,2738,4217,6494,10000,15399,23714,36517,56234,86596,
                133352,205353,316228,486968,749894,1154782,1778279,2738420,4216965,6493816,10000000,15399265,23713737,
                36517413,56234133,86596432,133352143,205352503,316227766,486967525,749894209,1154781985,1778279410,
                2738419634,4216965034,6493816316,10000000000};  // kappa

  TH1D * hsig[43][49];
  TH1D * hbkg[5][10];

  // Nested loops to iterate through arrays
  for (int i=0;i<43;i++){
    for (int j=0;j<49;j++){      
    TChain * chain = new TChain("t");
    chain->Add(Form("/data/COSINE/WORK/dfreitas/single_hit/DataFit/out_DP_SET3_bremss_m=%d_k=%ld.root",mDP[i],kappa[j]));

    gStyle->SetOptStat(000);

    double maxsig = chain->GetMaximum("sigma");
    cout<<maxsig<<endl;
    hsig[i][j] = new TH1D(Form("hsig%d%d",i,j), "", 1000, 0, maxsig*1.01);
    chain->Draw(Form("sigma>>hsig%d%d",i,j));
    delete chain;
    }
  }
  TFile fout("/data/COSINE/WORK/dfreitas/single_hit/sigma/sigma_DP_SET3_bremss.root", "recreate");
  for (int i=0;i<43;i++){
    for (int j=0;j<49;j++){  
      hsig[i][j]->Write();
    }
  }
  fout.Close();
}