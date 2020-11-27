#!/bin/bash

#SBATCH --qos=long
#SBATCH --job-name SF1_Trinity_test_25Jan2019
#SBATCH --error job_reports/SF1_Trinity_test_25Jan2019.error          
#SBATCH --output job_reports/SF1_Trinity_test_25Jan2019.output  
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load bowtie2/2.3.3 
module load samtools/1.5
module load anaconda3/5.1.0

F='/lustre/project/svanbael/steve/SF1/test_Trinity_Assembly/reads.left.fq.gz'
R='/lustre/project/svanbael/steve/SF1/test_Trinity_Assembly/reads.right.fq.gz'
O='Trinity_test_out_SF1_25Jan2019'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --left $F --right $R --output $O --monitoring --monitor_sec 10 --verbose

# Exit the Script

exit 0

