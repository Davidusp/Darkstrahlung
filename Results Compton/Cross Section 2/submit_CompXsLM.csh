#!/bin/bash/

dir=/home/dfreitas

for m in {0.0000000001000,0.000000001668,0.000000002783,0.000000004642,0.000000007743,0.000000012915,0.000000021544,0.000000035938,0.000000059948,0.000000100000}

do
	jobfile=$dir/jobs/DPCompton-xs_${m}.pbs
	errfile=$dir/logs/DPCompton-xs_${m}.err
	outfile=$dir/logs/DPCompton-xs_${m}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=12:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "python $dir/DPnum_Comp/DPCompton-xs.py ${m}" >> $jobfile
	qsub $jobfile
done
