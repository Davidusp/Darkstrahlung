#!/bin/bash/

dir=/home/dfreitas

for subrun in {0..481}

do

	jobfile=$dir/jobs/Run1544_${subrun}.pbs
	errfile=$dir/logs/Run1544_${subrun}.err
	outfile=$dir/logs/Run1544_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1544,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..129}

do

	jobfile=$dir/jobs/Run1546_${subrun}.pbs
	errfile=$dir/logs/Run1546_${subrun}.err
	outfile=$dir/logs/Run1546_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1546,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..649}

do

	jobfile=$dir/jobs/Run1616_${subrun}.pbs
	errfile=$dir/logs/Run1616_${subrun}.err
	outfile=$dir/logs/Run1616_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1616,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..226}

do

	jobfile=$dir/jobs/Run1617_${subrun}.pbs
	errfile=$dir/logs/Run1617_${subrun}.err
	outfile=$dir/logs/Run1617_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1617,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..6}

do

	jobfile=$dir/jobs/Run1625_${subrun}.pbs
	errfile=$dir/logs/Run1625_${subrun}.err
	outfile=$dir/logs/Run1625_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1625,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..8}

do

	jobfile=$dir/jobs/Run1626_${subrun}.pbs
	errfile=$dir/logs/Run1626_${subrun}.err
	outfile=$dir/logs/Run1626_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1626,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..104}

do

	jobfile=$dir/jobs/Run1634_${subrun}.pbs
	errfile=$dir/logs/Run1634_${subrun}.err
	outfile=$dir/logs/Run1634_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1634,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..802}

do

	jobfile=$dir/jobs/Run1652_${subrun}.pbs
	errfile=$dir/logs/Run1652_${subrun}.err
	outfile=$dir/logs/Run1652_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1652,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..216}

do

	jobfile=$dir/jobs/Run1654_${subrun}.pbs
	errfile=$dir/logs/Run1654_${subrun}.err
	outfile=$dir/logs/Run1654_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1654,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..297}

do

	jobfile=$dir/jobs/Run1666_${subrun}.pbs
	errfile=$dir/logs/Run1666_${subrun}.err
	outfile=$dir/logs/Run1666_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1666,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..301}

do

	jobfile=$dir/jobs/Run1671_${subrun}.pbs
	errfile=$dir/logs/Run1671_${subrun}.err
	outfile=$dir/logs/Run1671_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1671,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..955}

do

	jobfile=$dir/jobs/Run1672_${subrun}.pbs
	errfile=$dir/logs/Run1672_${subrun}.err
	outfile=$dir/logs/Run1672_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1672,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..707}

do

	jobfile=$dir/jobs/Run1678_${subrun}.pbs
	errfile=$dir/logs/Run1678_${subrun}.err
	outfile=$dir/logs/Run1678_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1678,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..74}

do

	jobfile=$dir/jobs/Run1683_${subrun}.pbs
	errfile=$dir/logs/Run1683_${subrun}.err
	outfile=$dir/logs/Run1683_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1683,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..5}

do

	jobfile=$dir/jobs/Run1690_${subrun}.pbs
	errfile=$dir/logs/Run1690_${subrun}.err
	outfile=$dir/logs/Run1690_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1690,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..929}

do

	jobfile=$dir/jobs/Run1718_${subrun}.pbs
	errfile=$dir/logs/Run1718_${subrun}.err
	outfile=$dir/logs/Run1718_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1718,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..992}

do

	jobfile=$dir/jobs/Run1719_${subrun}.pbs
	errfile=$dir/logs/Run1719_${subrun}.err
	outfile=$dir/logs/Run1719_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1719,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..97}

do

	jobfile=$dir/jobs/Run1720_${subrun}.pbs
	errfile=$dir/logs/Run1720_${subrun}.err
	outfile=$dir/logs/Run1720_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1720,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..479}

do

	jobfile=$dir/jobs/Run1771_${subrun}.pbs
	errfile=$dir/logs/Run1771_${subrun}.err
	outfile=$dir/logs/Run1771_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1771,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..228}

do

	jobfile=$dir/jobs/Run1777_${subrun}.pbs
	errfile=$dir/logs/Run1777_${subrun}.err
	outfile=$dir/logs/Run1777_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1777,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..994}

do

	jobfile=$dir/jobs/Run1858_${subrun}.pbs
	errfile=$dir/logs/Run1858_${subrun}.err
	outfile=$dir/logs/Run1858_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1858,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..991}

do

	jobfile=$dir/jobs/Run1859_${subrun}.pbs
	errfile=$dir/logs/Run1859_${subrun}.err
	outfile=$dir/logs/Run1859_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1859,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..998}

do

	jobfile=$dir/jobs/Run1860_${subrun}.pbs
	errfile=$dir/logs/Run1860_${subrun}.err
	outfile=$dir/logs/Run1860_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1860,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..986}

do

	jobfile=$dir/jobs/Run1861_${subrun}.pbs
	errfile=$dir/logs/Run1861_${subrun}.err
	outfile=$dir/logs/Run1861_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1861,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..993}

do

	jobfile=$dir/jobs/Run1862_${subrun}.pbs
	errfile=$dir/logs/Run1862_${subrun}.err
	outfile=$dir/logs/Run1862_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1862,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..315}

do

	jobfile=$dir/jobs/Run1863_${subrun}.pbs
	errfile=$dir/logs/Run1863_${subrun}.err
	outfile=$dir/logs/Run1863_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1863,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..341}

do

	jobfile=$dir/jobs/Run1865_${subrun}.pbs
	errfile=$dir/logs/Run1865_${subrun}.err
	outfile=$dir/logs/Run1865_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1865,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..57}

do

	jobfile=$dir/jobs/Run1866_${subrun}.pbs
	errfile=$dir/logs/Run1866_${subrun}.err
	outfile=$dir/logs/Run1866_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1866,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..1000}

do

	jobfile=$dir/jobs/Run1868_${subrun}.pbs
	errfile=$dir/logs/Run1868_${subrun}.err
	outfile=$dir/logs/Run1868_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1868,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..963}

do

	jobfile=$dir/jobs/Run1873_${subrun}.pbs
	errfile=$dir/logs/Run1873_${subrun}.err
	outfile=$dir/logs/Run1873_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1873,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..129}

do

	jobfile=$dir/jobs/Run1883_${subrun}.pbs
	errfile=$dir/logs/Run1883_${subrun}.err
	outfile=$dir/logs/Run1883_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1883,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..677}

do

	jobfile=$dir/jobs/Run1919_${subrun}.pbs
	errfile=$dir/logs/Run1919_${subrun}.err
	outfile=$dir/logs/Run1919_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1919,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..6}

do

	jobfile=$dir/jobs/Run1932_${subrun}.pbs
	errfile=$dir/logs/Run1932_${subrun}.err
	outfile=$dir/logs/Run1932_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1932,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..990}

do

	jobfile=$dir/jobs/Run1939_${subrun}.pbs
	errfile=$dir/logs/Run1939_${subrun}.err
	outfile=$dir/logs/Run1939_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1939,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..57}

do

	jobfile=$dir/jobs/Run1942_${subrun}.pbs
	errfile=$dir/logs/Run1942_${subrun}.err
	outfile=$dir/logs/Run1942_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1942,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..159}

do

	jobfile=$dir/jobs/Run1951_${subrun}.pbs
	errfile=$dir/logs/Run1951_${subrun}.err
	outfile=$dir/logs/Run1951_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1951,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..24}

do

	jobfile=$dir/jobs/Run1954_${subrun}.pbs
	errfile=$dir/logs/Run1954_${subrun}.err
	outfile=$dir/logs/Run1954_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1954,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..84}

do

	jobfile=$dir/jobs/Run1968_${subrun}.pbs
	errfile=$dir/logs/Run1968_${subrun}.err
	outfile=$dir/logs/Run1968_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1968,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..90}

do

	jobfile=$dir/jobs/Run1970_${subrun}.pbs
	errfile=$dir/logs/Run1970_${subrun}.err
	outfile=$dir/logs/Run1970_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1970,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..102}

do

	jobfile=$dir/jobs/Run1972_${subrun}.pbs
	errfile=$dir/logs/Run1972_${subrun}.err
	outfile=$dir/logs/Run1972_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1972,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..25}

do

	jobfile=$dir/jobs/Run1993_${subrun}.pbs
	errfile=$dir/logs/Run1993_${subrun}.err
	outfile=$dir/logs/Run1993_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1993,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..370}

do

	jobfile=$dir/jobs/Run1994_${subrun}.pbs
	errfile=$dir/logs/Run1994_${subrun}.err
	outfile=$dir/logs/Run1994_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1994,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..590}

do

	jobfile=$dir/jobs/Run1996_${subrun}.pbs
	errfile=$dir/logs/Run1996_${subrun}.err
	outfile=$dir/logs/Run1996_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(1996,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..116}

do

	jobfile=$dir/jobs/Run2012_${subrun}.pbs
	errfile=$dir/logs/Run2012_${subrun}.err
	outfile=$dir/logs/Run2012_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2012,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..575}

do

	jobfile=$dir/jobs/Run2025_${subrun}.pbs
	errfile=$dir/logs/Run2025_${subrun}.err
	outfile=$dir/logs/Run2025_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2025,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..64}

do

	jobfile=$dir/jobs/Run2030_${subrun}.pbs
	errfile=$dir/logs/Run2030_${subrun}.err
	outfile=$dir/logs/Run2030_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2030,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..518}

do

	jobfile=$dir/jobs/Run2039_${subrun}.pbs
	errfile=$dir/logs/Run2039_${subrun}.err
	outfile=$dir/logs/Run2039_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2039,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..13}

do

	jobfile=$dir/jobs/Run2045_${subrun}.pbs
	errfile=$dir/logs/Run2045_${subrun}.err
	outfile=$dir/logs/Run2045_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2045,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..579}

do

	jobfile=$dir/jobs/Run2046_${subrun}.pbs
	errfile=$dir/logs/Run2046_${subrun}.err
	outfile=$dir/logs/Run2046_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2046,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..568}

do

	jobfile=$dir/jobs/Run2048_${subrun}.pbs
	errfile=$dir/logs/Run2048_${subrun}.err
	outfile=$dir/logs/Run2048_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2048,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..482}

do

	jobfile=$dir/jobs/Run2050_${subrun}.pbs
	errfile=$dir/logs/Run2050_${subrun}.err
	outfile=$dir/logs/Run2050_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2050,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..108}

do

	jobfile=$dir/jobs/Run2054_${subrun}.pbs
	errfile=$dir/logs/Run2054_${subrun}.err
	outfile=$dir/logs/Run2054_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2054,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..296}

do

	jobfile=$dir/jobs/Run2079_${subrun}.pbs
	errfile=$dir/logs/Run2079_${subrun}.err
	outfile=$dir/logs/Run2079_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2079,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..291}

do

	jobfile=$dir/jobs/Run2080_${subrun}.pbs
	errfile=$dir/logs/Run2080_${subrun}.err
	outfile=$dir/logs/Run2080_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2080,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..252}

do

	jobfile=$dir/jobs/Run2083_${subrun}.pbs
	errfile=$dir/logs/Run2083_${subrun}.err
	outfile=$dir/logs/Run2083_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2083,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..96}

do

	jobfile=$dir/jobs/Run2103_${subrun}.pbs
	errfile=$dir/logs/Run2103_${subrun}.err
	outfile=$dir/logs/Run2103_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2103,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..68}

do

	jobfile=$dir/jobs/Run2118_${subrun}.pbs
	errfile=$dir/logs/Run2118_${subrun}.err
	outfile=$dir/logs/Run2118_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2118,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..122}

do

	jobfile=$dir/jobs/Run2119_${subrun}.pbs
	errfile=$dir/logs/Run2119_${subrun}.err
	outfile=$dir/logs/Run2119_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2119,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..516}

do

	jobfile=$dir/jobs/Run2126_${subrun}.pbs
	errfile=$dir/logs/Run2126_${subrun}.err
	outfile=$dir/logs/Run2126_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2126,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..75}

do

	jobfile=$dir/jobs/Run2128_${subrun}.pbs
	errfile=$dir/logs/Run2128_${subrun}.err
	outfile=$dir/logs/Run2128_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2128,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..124}

do

	jobfile=$dir/jobs/Run2134_${subrun}.pbs
	errfile=$dir/logs/Run2134_${subrun}.err
	outfile=$dir/logs/Run2134_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2134,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..168}

do

	jobfile=$dir/jobs/Run2135_${subrun}.pbs
	errfile=$dir/logs/Run2135_${subrun}.err
	outfile=$dir/logs/Run2135_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2135,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..295}

do

	jobfile=$dir/jobs/Run2142_${subrun}.pbs
	errfile=$dir/logs/Run2142_${subrun}.err
	outfile=$dir/logs/Run2142_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2142,${subrun}\)" >> $jobfile
	qsub $jobfile
done

for subrun in {0..348}

do

	jobfile=$dir/jobs/Run2148_${subrun}.pbs
	errfile=$dir/logs/Run2148_${subrun}.err
	outfile=$dir/logs/Run2148_${subrun}.txt

	echo "#PBS -l nodes=1:ppn=1" > $jobfile
	echo "#PBS -l walltime=02:00:00" >> $jobfile
	echo "#PBS -V" >> $jobfile
	echo "#PBS -q very_short" >> $jobfile
	echo "#PBS -e $errfile" >> $jobfile
	echo "#PBS -o $outfile" >> $jobfile

	echo "root -l -q -b $dir/Rdata/crystalpreselection.C\(2148,${subrun}\)" >> $jobfile
	qsub $jobfile
done
