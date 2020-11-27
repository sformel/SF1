#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name rm_partitions
#SBATCH --error job_reports/rm_partitions.error          
#SBATCH --output job_reports/rm_partitions.output  
#SBATCH --time=4-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu


rm -r /lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test/read_partitions
