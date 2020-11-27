#!/bin/bash

#SBATCH --qos=long
#SBATCH --error job_reports/SF1_Trinity_18Jan2019.error          
#SBATCH --output job_reports/SF1_Trinity_18Jan2019.output  
#SBATCH --job-name SF1_Trinity_test_18Jan2019
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
O='Trinity_test_out_SF1_18Jan2019'

Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 20 --bflyCalculateCPU --left $F --right $R --output $O --monitoring --full_cleanup

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa trinity_out_dir/Trinity.fasta 90

    ./test_FL.sh --query trinity_out_dir/Trinity.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse
fi
