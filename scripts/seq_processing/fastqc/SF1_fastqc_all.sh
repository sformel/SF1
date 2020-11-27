#!/bin/bash

#SBATCH --job-name=SF1_fastqc_all
#SBATCH --output=job_reports/SF1_fastqc_all.output
#SBATCH --error=job_reports/SF1_fastqc_all.error
#SBATCH --qos=normal
#SBATCH --time=0-24:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#run multiprocessed version of fastqc

ls /lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/*.fq.gz | xargs -P 20 -n 1 fastqc
