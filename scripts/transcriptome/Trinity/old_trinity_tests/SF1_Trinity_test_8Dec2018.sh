#!/bin/bash

#SBATCH --qos=long
#SBATCH --error job_reports/SF1_Trinity_8Dec2018.error          
#SBATCH --output job_reports/SF1_Trinity_8Dec2018.output  
#SBATCH --job-name SF1_Trinity_test_8Dec2018
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load bowtie2/2.3.3 
module load samtools/1.5
module load anaconda3/5.1.0

F='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1/SF1_01_1_paired.fq.gz'
R='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R2/SF1_01_2_paired.fq.gz'
O='Trinity_test_out_SF1_7Dec2018'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --left $F --right $R --output $O --monitoring
