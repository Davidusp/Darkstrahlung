void doFit(long int mdm, long int kappa){
   gROOT->ProcessLine(".L /data/COSINE/WORK/dfreitas/single_hit/src/CAT.hh");
   gROOT->ProcessLine(".L /data/COSINE/WORK/dfreitas/single_hit/src/CAT.cc");
   gROOT->ProcessLine(Form(".x /data/COSINE/WORK/dfreitas/single_hit/doFit_DPbremss_SET3.C(%ld, %ld)",mdm,kappa));
}