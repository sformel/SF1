#Script for SF1
#Purpose: Import 16S phyloseq object and replace/clean sample data

#Last updated by Steve FOrmel
#Jan 21, 2020

#make phyloseq object-----
library(phyloseq)
library(stringr)

#read in data
seqtab.nochim <- readRDS("data/metagenome/16S/dada2/16S_soil_otu.rds")
samdf <- readRDS("data/metagenome/16s/dada2/16S_soil_env.rds")
taxa <- readRDS("data/metagenome/16s/dada2/16S_soil_taxa.rds")

SF1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#add sampleID to match ITS
sample_data(SF1)$sampleID <- paste0("SF1_", rownames(sample_data(SF1)))
sample_data(SF1)$sampleID <- str_replace_all(sample_data(SF1)$sampleID,
                                             "LN",
                                             "SLN")
sample_data(SF1)$sampleID <- str_replace_all(sample_data(SF1)$sampleID,
                                             "LY",
                                             "SLY")
sample_data(SF1)$sampleID <- str_replace_all(sample_data(SF1)$sampleID,
                                             "SN",
                                             "SSN")
sample_data(SF1)$sampleID <- str_replace_all(sample_data(SF1)$sampleID,
                                             "SY",
                                             "SSY")

#Replace sample data-----
library(readxl)

samdf <- as.data.frame(read_excel("data/soil_chem/SF1_soil_chem_all_data.xls", 
                                  sheet = "data_cleaned"))

#reorder to match
samdf <- samdf[match(sample_data(SF1)$sampleID, samdf$sampleID),]

#do any sampleID not match?
any(sample_data(SF1)$sampleID!=samdf$sampleID)

#replace rownames on both
rownames(samdf) <- samdf$sampleID
sample_names(SF1) <- sample_data(SF1)$sampleID

#double check that rownames matches sampleID
any(rownames(samdf)!=samdf$sampleID)

#cleaning and fixing----

#reorder factors
samdf$time <- factor(samdf$time, levels = c("pre-exp","post-exp"))

#rename factors  
library(plyr)

samdf$soil <- revalue(samdf$soil, c("L"="Live Soil", "S"="Autoclaved Soil"))
samdf$oil <- revalue(samdf$oil, c("Y"="Oiled", "N"="No Oil"))
samdf$time <- revalue(samdf$time, c("pre-exp"="Pre-Experiment", "post-exp"="Post-Experiment"))

#rename chem titles
colnames(samdf) <- c("sampleID", "soil", "oil", "time", "Calcium", "Copper", "Magnesium", "pH", "Phosphorus", "Potassium",  "Sodium", "Sulfur", "Zinc", "Percent Carbon", "Percent Nitrogen")

samdf$Treatment <- interaction(samdf$soil, samdf$oil, sep = " , ")
samdf$Treatment <- factor(samdf$Treatment, levels = c("Live Soil , No Oil","Live Soil , Oiled", "Autoclaved Soil , No Oil", "Autoclaved Soil , Oiled"))

#put new data into pseq object
sample_data(SF1) <- samdf


#clean up taxa names to something readable instead of sequence
taxa_names(SF1) <- paste0("Seq", seq(ntaxa(SF1)))

rm(samdf)
rm(seqtab.nochim)
rm(taxa)