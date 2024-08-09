#!/bin/bash/

dir=/data/COSINE/WORK/dfreitas/single_hit

for mdm in {10,20,30,40,50,70,100,130,170,220,290,390,520,690,910,1210,1600,2120,2810,3730,4940,6550,8690,11510,15260,11514,20240,26830,35560,47150,62510,82860,109850,145630,193070,255950,339320,449840,596360,790600,1048110,1389500,1842070}

do
	for kappa in	{10,15,24,37,56,87,133,205,316,487,750,1155,1778,2738,4217,6494,10000,15399,23714,36517,56234,86596,133352,205353,316228,486968,749894,1154782,1778279,2738420,4216965,6493816,10000000,15399265,23713737,36517413,56234133,86596432,133352143,205352503,316227766,486967525,749894209,1154781985,1778279410,2738419634,4216965034,6493816316,10000000000}

	do
		jobfile=$dir/jobs/doFit_${mdm}_${kappa}.pbs
		errfile=$dir/logs/doFit_${mdm}_${kappa}.err
		outfile=$dir/logs/doFit_${mdm}_${kappa}.txt

		echo "#PBS -l nodes=1:ppn=1" > $jobfile
		echo "#PBS -l walltime=02:00:00" >> $jobfile
		echo "#PBS -V" >> $jobfile
		echo "#PBS -q very_short" >> $jobfile
		echo "#PBS -e $errfile" >> $jobfile
		echo "#PBS -o $outfile" >> $jobfile

		echo "root -l -q -b \"$dir/doFit.C($mdm,$kappa)\"" >> $jobfile
		qsub $jobfile

		#root -l -q -b doFit.C\($mdm\)
	done
done
