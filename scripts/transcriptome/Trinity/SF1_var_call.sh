#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name SF1_var_call
#SBATCH --error job_reports/SF1_var_call.error
#SBATCH --output job_reports/SF1_var_call.output
#SBATCH --time=4-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

	
module load trinity/2.8.4
module load star/2.5.2a

#trinity output
TO='/lustre/project/svanbael/steve/SF1/Trinity_test_out_SF1_20Feb2019_test'

${TRINITY_HOME}/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa trinity_genes.fasta \
--st_gtf trinity_genes.gtf \
-S sample_lists/SF1_Trinity_samples_list.txt \
-t 20 \
-o $TO/variant_calling \
-m 250000000000 \
--STAR_genomeGenerate_opts '--genomeChrBinNbits 10 --limitSjdbInsertNsj 1527404 --limitOutSJcollapsed 5000000 --limitIObufferSize 300000000'
