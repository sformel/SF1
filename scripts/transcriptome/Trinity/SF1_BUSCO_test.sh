#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_BUSCO_test
#SBATCH --error job_reports/SF1_BUSCO_test.error          
#SBATCH --output job_reports/SF1_BUSCO_test.output  
#SBATCH --time=4-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu
	
#trinity output

TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

#BUSCO db
BDB='/lustre/project/svanbael/steve/software/BUSCO_git/'

module load BUSCO/3
	
run_BUSCO.py -i $TO/Trinity.fasta -o SF1_BUSCO -l $BDB/embryophyta_odb9 -m tran --cpu 20
	

