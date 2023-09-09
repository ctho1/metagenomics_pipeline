#!/bin/bash

in_files=$(find ./fastq -type f -name "*R1_001.fastq.gz" -print|sort)

for R1 in $in_files; do
	R2=`echo $R1 | sed 's/R1_001/R2_001/'`
	R1_base=`basename $R1 .fastq.gz`
	echo "submitting job ${R1_base}"
	sbatch ./slurm_script.sh $R1 $R2
done
