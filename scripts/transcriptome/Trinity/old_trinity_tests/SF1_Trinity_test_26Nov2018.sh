#!/bin/bash

#SBATCH --qos=long
#SBATCH --error job_reports/SF1_Trinity_26Nov2018.error          
#SBATCH --output job_reports/SF1_Trinity_26Nov2018.output  
#SBATCH --job-name SF1_Trinity_26Nov2018
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load trinity/2.4.0
module load bowtie2/2.3.3 

F='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_01_1_paired.fq.gz'
R='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_01_2_paired.fq.gz'
O='Trinity_test_out_SF1_26Nov2018'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --left $F --right $R --output $O
