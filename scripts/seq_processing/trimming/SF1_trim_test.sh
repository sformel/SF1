#!/bin/bash

#SBATCH --job-name=SF1_trim_test
#SBATCH --output=job_reports/SF1_trim_test.output
#SBATCH --error=job_reports/SF1_trim_test.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#test trimmomatic settings

DATAPATH='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data'
OUTPATH='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_test'

java -jar /lustre/project/svanbael/steve/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 $DATAPATH/SF1_01_1.fq.gz $DATAPATH/SF1_01_2.fq.gz $OUTPATH/output_forward_paired.fq.gz  $OUTPATH/output_forward_unpaired.fq.gz $OUTPATH/output_reverse_paired.fq.gz $OUTPATH/output_reverse_unpaired.fq.gz ILLUMINACLIP:/lustre/project/svanbael/steve/software/Trimmomatic-0.38/adapters/SF1_novogene-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
