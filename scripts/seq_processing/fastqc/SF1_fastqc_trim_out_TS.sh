#!/bin/bash

#SBATCH --job-name=SF1_fastqc_trim_out_TS
#SBATCH --output=job_reports/SF1_fastqc_trim_out_TS.output
#SBATCH --error=job_reports/SF1_fastqc_trim_out_TS.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#run multiprocessed version of fastqc on output from trimmomatic after using TruSeq adapter seqs.

ls /lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out_TS/*.fq.gz | xargs -P 20 -n 1 fastqc
