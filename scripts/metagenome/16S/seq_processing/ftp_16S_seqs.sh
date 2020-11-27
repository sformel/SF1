#!/bin/bash

#SBATCH --qos=normal
#SBATCH --job-name  ftp_16S_seqs
#SBATCH --error  job_reports/ftp_16S_seqs.error
#SBATCH --output  job_reports/ftp_16S_seqs.output
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sformel@tulane.edu

#downloading 16S seqs from Novogene ftp
#downloaded on Sept 6, 2019

wget -r -l 20 --user="X202SC19080294-Z01_09_05_19_xBOA" --password="vpKtRt2v" ftp://usftp1.novogene.com

#move to VB lab box.com

module load Rclone/1.47

SF1=/lustre/project/svanbael/steve/SF1/soil_16S
VB_box=VBlab_box:Stephen_Formel/SF1/soil_16S
	
#check to make sure VV box path exists
rclone mkdir $VB_box
	
rclone copy $SF1 $VB_box  --tpslimit 10
	
#check that rclone worked
rclone check $SF1 $VB_box --log-file $SF1/ftp_16S_seqs.log
