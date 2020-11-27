#!/bin/bash

#SBATCH --qos=long
#SBATCH --job-name SF1_Trinity_test_3Feb2019_test
#SBATCH --error job_reports/SF1_Trinity_test_3Feb2019_test.error          
#SBATCH --output job_reports/SF1_Trinity_test_3Feb2019_test.output  
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load bowtie2/2.3.3 
module load samtools/1.5
module load anaconda3/5.1.0
module load trinity/2.4.0

#forward reads
F='/lustre/project/svanbael/steve/SF1/test_Trinity_Assembly/reads.left.fq.gz'

#reverse reads
R='/lustre/project/svanbael/steve/SF1/test_Trinity_Assembly/reads.right.fq.gz'


#output folder
O='Trinity_test_out_SF1_3Feb2019'

#sample list
S='/lustre/project/svanbael/steve/SF1/SF1_Trinity_samples_list.txt'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --left $F --right $R  --output $TMPDIR/SF_Trinity_3Feb2019_test --full_cleanup --grid_exec "/lustre/project/svanbael/steve/software/HPC_grid_runner/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl --grid_conf /lustre/project/svanbael/steve/software/HPC_grid_runner/HpcGridRunner-1.0.2/hpc_conf/SLURM.Cypress.SF1.conf -c"

mkdir ./$O

cp -r /local/tmp/sformel.*/SF_Trinity_3Feb2019_test.Trinity.fasta ./$O/


