#!/bin/bash

#SBATCH --qos=long
#SBATCH --job-name SF1_Trinity_test_1Feb2019_3
#SBATCH --error job_reports/SF1_Trinity_test_1Feb2019_3.error          
#SBATCH --output job_reports/SF1_Trinity_test_1Feb2019_3.output  
#SBATCH --time=2-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

module load bowtie2/2.3.3 
module load samtools/1.5
module load anaconda3/5.1.0
module load trinity/2.4.0

#forward reads
#F='/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/R1'

#reverse reads
#R='/lustre/project/svanbael/steve/Data/C101HW18101466/raw_data/trim_out/paired/R2'

#output folder
O='Trinity_test_out_SF1_1Feb2019_3'

#sample list
S='/lustre/project/svanbael/steve/SF1/SF1_Trinity_samples_list.txt'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --samples_file $S  --output $O --full_cleanup --grid_exec "/lustre/project/svanbael/steve/software/HPC_grid_runner/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl --grid_conf /lustre/project/svanbael/steve/software/HPC_grid_runner/HpcGridRunner-1.0.2/hpc_conf/SLURM.Cypress.SF1.conf -c"


