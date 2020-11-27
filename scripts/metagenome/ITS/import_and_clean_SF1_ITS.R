#Script for SF1
#Purpose: Import ITS phyloseq object and replace/clean sample data

#Last updated by Steve FOrmel
#Jan 29, 2020

#load libraries---
library(readxl)
library(phyloseq)
library(plyr)

#Import phyloseq objects-----
#made in script called "ITS_seqs_BI.R"

SF1 <- readRDS("../../../Grad Students/Steve/S5_SF1_ITS_seqs_BI/SF1_ITS_ps.rds")

#Replace sample data-----

samdf <- as.data.frame(read_excel("data/soil_chem/SF1_soil_chem_all_data.xls", 
                                  sheet = "data_cleaned"))

#do sampleID match?
any(sample_data(SF1)$sampleID!=samdf$sampleID)

#match rownames
rownames(samdf) <- rownames(sample_data(SF1))

#double check that rownames matches sampleID
any(rownames(samdf)!=samdf$sampleID)


#cleaning and fixing----

#reorder factors
samdf$time <- factor(samdf$time, levels = c("pre-exp","post-exp"))

#rename factors  

samdf$soil <- revalue(samdf$soil, c("L"="Live Soil", "S"="Autoclaved Soil"))
samdf$oil <- revalue(samdf$oil, c("Y"="Oiled", "N"="No Oil"))
samdf$time <- revalue(samdf$time, c("pre-exp"="Pre-Experiment", "post-exp"="Post-Experiment"))

#rename chem titles
colnames(samdf) <- c("sampleID", "soil", "oil", "time", "Calcium", "Copper", "Magnesium", "pH", "Phosphorus", "Potassium",  "Sodium", "Sulfur", "Zinc", "Percent Carbon", "Percent Nitrogen")

samdf$Treatment <- interaction(samdf$soil, samdf$oil, sep = " , ")
samdf$Treatment <- factor(samdf$Treatment, levels = c("Live Soil , No Oil","Live Soil , Oiled", "Autoclaved Soil , No Oil", "Autoclaved Soil , Oiled"))

#put new data into pseq object
sample_data(SF1) <- samdf

rm(samdf)
