#!/bin/bash/

dir=/home/dfreitas

for m in {1.2,1.57079914,2.05617494,2.69153152,3.52321282,4.61188305,6.03695159,7.90236529,10.34419048,13.54053789,17.72455436,23.20142891,30.37065375,39.75516391,52.03948096,68.11964314,89.16856387,116.72158596,152.78847205,200.}

do
	jobfile=$dir/jobs/DPCompton-xs_${m}.pbs
	errfile=$dir/logs/DPCompton-xs_${m}.err
	outfile=$dir/logs/DPCompton-xs_${m}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=48:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q medium" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "python $dir/DPnum_Comp/DPCompton-xs.py ${m}" >> $jobfile
	qsub $jobfile
done
