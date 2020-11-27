#!/bin/bash
	
#SBATCH --qos=long
#SBATCH --job-name bbduk_bac_contam
#SBATCH --error job_reports/bbduk_bac_contam.error
#SBATCH --output job_reports/bbduk_bac_contam.output
#SBATCH --time=3-23:00:00
#SBATCH --nodes=1
#SBATCH --mem=128000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu
	
module load BBtools/38.41
	
WORKDIR=/lustre/project/svanbael/steve/SF1/Data/C101HW18101466/raw_data/trim_out/paired/
	
ANNODIR=/lustre/project/svanbael/steve/SF1/BLAST_db/bac_contam
DATADIR=/lustre/project/svanbael/steve/SF1/
	
SAMPLIST='SF1_01 SF1_04 SF1_10 SF1_13 SF1_17 SF1_19 SF1_26 SF1_30 SF1_34 SF1_37 SF1_40 SF1_44'
PPN=20
	
contaminants=$ANNODIR/all_NCBI_bac_genomes.fa
	
for sample in $SAMPLIST
do
bbduk.sh ref=${contaminants} \
threads=${PPN} \
ordered=t \
k=31 \
in=${DATADIR}/R1/${sample}_1_paired.fq.gz \
in2=${DATADIR}/R2/${sample}_2_paired.fq.gz \
out=${DATADIR}/bbduk/${sample}.R1.fq \
out2=${DATADIR}/bbduk/${sample}.R2.fq \
outm=${DATADIR}/bbduk/${sample}.bad.R1.fq \
outm2=${DATADIR}/bbduk/${sample}.bad.R2.fq
done 
