#!/bin/bash

#SBATCH --qos=long
#SBATCH --job-name SF1_Trinity_test_10Feb2019_test
#SBATCH --error job_reports/SF1_Trinity_test_10Feb2019_test.error          
#SBATCH --output job_reports/SF1_Trinity_test_10Feb2019_test.output  
#SBATCH --time=1-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load trinity/2.8.4

#forward reads
#F='/lustre/project/svanbael/steve/SF1/test_Trinity_Assembly/reads.left.fq.gz'

#reverse reads
#R='/lustre/project/svanbael/steve/SF1/test_Trinity_Assembly/reads.right.fq.gz'


#output folder
O='Trinity_test_out_SF1_10Feb2019_test'

#sample list
S='/lustre/project/svanbael/steve/SF1/sample_lists/SF1_Trinity_samples_list_test1.txt'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --samples_file $S --output $O 

#--full_cleanup --grid_exec "hpc_cmds_GridRunner.pl --grid_conf /lustre/project/svanbael/steve/software/HPC_grid_runner/HpcGridRunner-1.0.2/hpc_conf/SLURM.Cypress.SF1.conf -c"


