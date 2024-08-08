#!/bin/bash/

dir=/home/dfreitas

for m in {0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,0.000039,0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000}

do
	jobfile=$dir/jobs/CompNumbeta_${m}.pbs
	errfile=$dir/logs/CompNumbeta_${m}.err
	outfile=$dir/logs/CompNumbeta_${m}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=12:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "python $dir/DPnum_Comp/totalDPs/CompNumbeta.py ${m}" >> $jobfile
	qsub $jobfile
done
