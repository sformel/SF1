#Make heat maps
#script received from Sarah Khalil on Jan 22, 2020

library(optparse)
library(DESeq2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(data.table)
library(reshape)
library(ggrepel)
library(apeglm)
library(pheatmap)
library(GeneNet)

#heatmap examples#

test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)
annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)), 
  Time = 1:5
)
rownames(annotation_col) = paste("Test", 1:10, sep = "")

annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")

pheatmap(test, annotation_col = annotation_col)


#MY DATA for HEATMAPS
#import metadatafile

sample_id <- dir(file.path("scripts/transcriptome/R/example_heatmap_from_SK/"))
kal_dirs <- file.path("scripts/transcriptome/R/example_heatmap_from_SK/", sample_id)
directory <- "scripts/transcriptome/R/example_heatmap_from_SK"
s2c <- read.table("scripts/transcriptome/R/example_heatmap_from_SK/sample_table_RBFW_RNAseq_modSF.csv", header = TRUE, stringsAsFactors = TRUE, sep = ",", fileEncoding="UTF-8-BOM")

s2c<-dplyr::filter(s2c,sampleName!="F84") #bad sample

s2c <- s2c %>% mutate(path = paste("scripts/transcriptome/R/example_heatmap_from_SK/",fileName, sep = ""))

hm_hts<-DESeqDataSetFromHTSeqCount(sampleTable=s2c, directory=directory, design=~PhenotypePart)
hm_dds<-DESeq(hm_hts)
hm.vsd <- vst(hm_dds, blind=TRUE)
head(assay(hm.vsd), 20)
counts<-(assay(hm.vsd)) #get count data per gene - this does not have any negatives, so i will z-transform to get that
counts_scale<-scale(counts) #z-tranform count data

#make column annotations for pheatmap
column_annotation<-s2c[,-1]
rownames(column_annotation)<-s2c[,1]
column_annotation_samplephenotype<- column_annotation %>% select(samplePhenotype)

###subset to only include genes different (made up by SF) to subset counts dataframe

liver_dm_bm <- as.vector(as.matrix(read.delim("scripts/transcriptome/R/example_heatmap_from_SK/SF_DE_gene_list", header=T)))

counts_liver_dm_bm_sig_up_genes<-counts[liver_dm_bm] #subset for just liver 
pheatmap(counts_liver_dm_bm_sig_up_genes, annotation_col = column_annotation_samplephenotype)

counts_scale_liver_dm_bm_sig_up_genes<-counts_scale[liver_dm_bm_sig_up_genes_list,26:37] #subset for just liver 
pheatmap(counts_scale_liver_dm_bm_sig_up_genes, annotation_col = column_annotation_samplephenotype)

#################################################################################
###subset to only include genes different in liver of dull and Timplanted males###
#################################################################################
liver_dm_tdm<-read.delim("output/deseq_DE_output/liver_timplantdm_vs_dm_DE_genes.txt", header=T)
liver_dm_tdm$Lpadj<--log10(liver_dm_tdm$padj) #create -log10 (FDR)
liver_dm_tdm<-na.omit(liver_dm_tdm)
liver_dm_tdm$padj.t<-ifelse(liver_dm_tdm$log2FoldChange>0,liver_dm_tdm$Lpadj*-1,liver_dm_tdm$Lpadj)

#filter genes that are signifcantly upregulated
liver_dm_tdm_sig_up<-dplyr::filter(liver_dm_tdm,padj.t<(-1)) 
liver_dm_tdm_sig_up_genes<-liver_dm_tdm_sig_up[1]
liver_dm_tdm_sig_up_genes<-unique(liver_dm_tdm_sig_up_genes)
#make a list of UPREGULATED genes between dm and bm in liver, to subset counts dataframe
liver_dm_tdm_sig_up_genes_list<-as.list(liver_dm_tdm_sig_up_genes[,1]) 
liver_dm_tdm_sig_up_genes_list<-as.vector(as.matrix(liver_dm_tdm_sig_up_genes))

counts_scale_liver_dm_tdm_sig_up_genes<-counts_scale[liver_dm_tdm_sig_up_genes_list,26:37] #subset for just liver 
pheatmap(counts_scale_liver_dm_tdm_sig_up_genes, annotation_col = column_annotation_samplephenotype)

#filter genes that are signifcantly DE
liver_dm_tdm_sig<-dplyr::filter(liver_dm_tdm,Lpadj>1) 
liver_dm_tdm_sig_genes<-liver_dm_tdm_sig[1]
liver_dm_tdm_sig_genes<-unique(liver_dm_tdm_sig_genes)
#make a list of differentially expressed genes between DM and TDM in liver, to subset counts dataframe
liver_dm_tdm_sig_genes_list<-as.list(liver_dm_tdm_sig_genes[,1]) 
liver_dm_tdm_sig_genes_list<-as.vector(as.matrix(liver_dm_tdm_sig_genes))
#subset for just liver (not feathers) 
counts_scale_liver_dm_tdm_sig_genes<-counts_scale[liver_dm_tdm_sig_genes_list,26:37]

#replace row names with gene names
df.counts_scale_liver_dm_tdm_sig_genes<-as.data.frame(counts_scale_liver_dm_tdm_sig_genes) #make into dataframe
df.counts_scale_liver_dm_tdm_sig_genes<-setDT(df.counts_scale_liver_dm_tdm_sig_genes, keep.rownames = TRUE)[] #turn row name to colum
df.merge.counts_scale_liver_dm_tdm_sig_genes<-merge(df.counts_scale_liver_dm_tdm_sig_genes,annie_rowname_unique,by.x="rn",by.y="gene_name") #merge with gene names
df.merge.counts_scale_liver_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_liver_dm_tdm_sig_genes,rn!="RBFW011550")

#make gene name row names
#but first replace CYP2J2 with CYP2J19
df.merge.counts_scale_liver_dm_tdm_sig_genes$blast[df.merge.counts_scale_liver_dm_tdm_sig_genes$blast=="CYP2J2"]<-"CYP2J19"
rownames(df.merge.counts_scale_liver_dm_tdm_sig_genes)<-df.merge.counts_scale_liver_dm_tdm_sig_genes$blast
counts_scale_liver_dm_tdm_sig_genes_rename<-df.merge.counts_scale_liver_dm_tdm_sig_genes[,-1]
counts_scale_liver_dm_tdm_sig_genes_rename<-counts_scale_liver_dm_tdm_sig_genes_rename[,-13]


pheatmap(counts_scale_liver_dm_tdm_sig_genes, annotation_col = column_annotation_samplephenotype)
pheatmap(counts_scale_liver_dm_tdm_sig_genes_rename, annotation_col = column_annotation_samplephenotype)


##use just count data without z tranform##
counts_liver_dm_tdm_sig_genes<-counts[liver_dm_tdm_sig_genes_list,26:37]
#replace row names with gene names
df.counts_liver_dm_tdm_sig_genes<-as.data.frame(counts_liver_dm_tdm_sig_genes) #make into dataframe
df.counts_liver_dm_tdm_sig_genes<-setDT(df.counts_liver_dm_tdm_sig_genes, keep.rownames = TRUE)[] #turn row name to colum
df.merge.counts_liver_dm_tdm_sig_genes<-merge(df.counts_liver_dm_tdm_sig_genes,annie_rowname_unique,by.x="rn",by.y="gene_name") #merge with gene names
df.merge.counts_liver_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_liver_dm_tdm_sig_genes,rn!="RBFW011550")

#make gene name row names
#but first replace CYP2J2 with CYP2J19
df.merge.counts_liver_dm_tdm_sig_genes$blast[df.merge.counts_liver_dm_tdm_sig_genes$blast=="CYP2J2"]<-"CYP2J19"
rownames(df.merge.counts_liver_dm_tdm_sig_genes)<-df.merge.counts_liver_dm_tdm_sig_genes$blast
counts_liver_dm_tdm_sig_genes_rename<-df.merge.counts_liver_dm_tdm_sig_genes[,-1]
counts_liver_dm_tdm_sig_genes_rename<-counts_liver_dm_tdm_sig_genes_rename[,-13]

pheatmap(counts_liver_dm_tdm_sig_genes_rename, annotation_col = column_annotation_samplephenotype)
#qualitatively the same as z-transformed data

#################################################################################
###subset to only include genes different in liver of dull and bright males###
#################################################################################
liver_dm_bm<-read.delim("output/deseq_DE_output/liver_dm_vs_bm_DE_genes.txt", header=T)
liver_dm_bm$Lpadj<--log10(liver_dm_bm$padj) #create -log10 (FDR)
liver_dm_bm<-na.omit(liver_dm_bm)
liver_dm_bm$padj.t<-ifelse(liver_dm_bm$log2FoldChange>0,liver_dm_bm$Lpadj*-1,liver_dm_bm$Lpadj)

#filter genes that are signifcantly upregulated
liver_dm_bm_sig_up<-dplyr::filter(liver_dm_bm,padj.t<(-1)) 
liver_dm_bm_sig_up_genes<-liver_dm_bm_sig_up[1]
liver_dm_bm_sig_up_genes<-unique(liver_dm_bm_sig_up_genes)
#make a list of UPREGULATED genes between dm and bm in liver, to subset counts dataframe
liver_dm_bm_sig_up_genes_list<-as.list(liver_dm_bm_sig_up_genes[,1]) 
liver_dm_bm_sig_up_genes_list<-as.vector(as.matrix(liver_dm_bm_sig_up_genes))

counts_scale_liver_dm_bm_sig_up_genes<-counts_scale[liver_dm_bm_sig_up_genes_list,26:37] #subset for just liver 
pheatmap(counts_scale_liver_dm_bm_sig_up_genes, annotation_col = column_annotation_samplephenotype)

#filter genes that are signifcantly DE
liver_dm_bm_sig<-dplyr::filter(liver_dm_bm,Lpadj>1) 
liver_dm_bm_sig_genes<-liver_dm_bm_sig[1]
liver_dm_bm_sig_genes<-unique(liver_dm_bm_sig_genes)
#make a list of differentially expressed genes between DM and TDM in liver, to subset counts dataframe
liver_dm_bm_sig_genes_list<-as.list(liver_dm_bm_sig_genes[,1]) 
liver_dm_bm_sig_genes_list<-as.vector(as.matrix(liver_dm_bm_sig_genes))
counts_scale_liver_dm_bm_sig_genes<-counts_scale[liver_dm_bm_sig_genes_list,26:37] #subset for just liver 

#replace row names with gene names
df.counts_scale_liver_dm_bm_sig_genes<-as.data.frame(counts_scale_liver_dm_bm_sig_genes) #make into dataframe
df.counts_scale_liver_dm_bm_sig_genes<-setDT(df.counts_scale_liver_dm_bm_sig_genes, keep.rownames = TRUE)[] #turn row name to colum
df.merge.counts_scale_liver_dm_bm_sig_genes<-merge(df.counts_scale_liver_dm_bm_sig_genes,annie_rowname_unique,by.x="rn",by.y="gene_name") #merge with gene names
df.merge.counts_scale_liver_dm_bm_sig_genes<-dplyr::filter(df.merge.counts_scale_liver_dm_bm_sig_genes,rn!="RBFW011550")
#make gene name row names
#but first replace CYP2J2 with CYP2J19
df.merge.counts_scale_liver_dm_bm_sig_genes$blast[df.merge.counts_scale_liver_dm_bm_sig_genes$blast=="CYP2J2"]<-"CYP2J19"
#duplicated(df.merge.counts_scale_liver_dm_bm_sig_genes$blast) #one of them is duplicated
##remove the row with 2 of the same gene names
df.merge.counts_scale_liver_dm_bm_sig_genes<-df.merge.counts_scale_liver_dm_bm_sig_genes[-c(104),]
rownames(df.merge.counts_scale_liver_dm_bm_sig_genes)<-df.merge.counts_scale_liver_dm_bm_sig_genes$blast
counts_scale_liver_dm_bm_sig_genes_rename<-df.merge.counts_scale_liver_dm_bm_sig_genes[,-1]
counts_scale_liver_dm_bm_sig_genes_rename<-counts_scale_liver_dm_bm_sig_genes_rename[,-13]


pheatmap(counts_scale_liver_dm_bm_sig_genes, annotation_col = column_annotation_samplephenotype)
pheatmap(counts_scale_liver_dm_bm_sig_genes_rename, annotation_col = column_annotation_samplephenotype)
#this is a hot mess, don't use

###for replacing row names###
annie_rowname<-annie
annie_rowname_unique<-unique(annie_rowname)
###

#################################################################################
###subset to only include genes different in mantle feathers of dull and Timplanted males###
#################################################################################
mantle_dm_tdm<-read.delim("output/deseq_DE_output/mantle_red_vs_brown_Timplant_DE_genes.txt", header=T)
mantle_dm_tdm$Lpadj<--log10(mantle_dm_tdm$padj) #create -log10 (FDR)
mantle_dm_tdm<-na.omit(mantle_dm_tdm)
mantle_dm_tdm$padj.t<-ifelse(mantle_dm_tdm$log2FoldChange>0,mantle_dm_tdm$Lpadj*-1,mantle_dm_tdm$Lpadj)

#filter genes that are signifcantly DE
mantle_dm_tdm_sig<-dplyr::filter(mantle_dm_tdm,Lpadj>1) 
mantle_dm_tdm_sig_genes<-mantle_dm_tdm_sig[1]
mantle_dm_tdm_sig_genes<-unique(mantle_dm_tdm_sig_genes)
#make a list of differentially expressed genes between DM and TDM in mantle, to subset counts dataframe
mantle_dm_tdm_sig_genes_list<-as.list(mantle_dm_tdm_sig_genes[,1]) 
mantle_dm_tdm_sig_genes_list<-as.vector(as.matrix(mantle_dm_tdm_sig_genes))
mantle_feathers<-c("F109","F91","F38","F41","F55","F82","F101","F87","F95","F34","F47","F53")
counts_scale_mantle_dm_tdm_sig_genes<-counts_scale[mantle_dm_tdm_sig_genes_list,mantle_feathers] #subset for just mantle 

df.counts_scale_mantle_dm_tdm_sig_genes<-as.data.frame(counts_scale_mantle_dm_tdm_sig_genes) #make into dataframe
df.counts_scale_mantle_dm_tdm_sig_genes<-setDT(df.counts_scale_mantle_dm_tdm_sig_genes, keep.rownames = TRUE)[] #turn row name to colum
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-merge(df.counts_scale_mantle_dm_tdm_sig_genes,annie_rowname_unique,by.x="rn",by.y="gene_name") #merge with gene names
# 2 rows were added during merging, because of the same gene number being associated with different gene names
# remove the duplicate manually
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_mantle_dm_tdm_sig_genes,blast!="Prps2")
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_mantle_dm_tdm_sig_genes,rn!="RBFW001191")
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_mantle_dm_tdm_sig_genes,blast!="Gpr158")
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_mantle_dm_tdm_sig_genes,blast!="Ranbp3l")
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_mantle_dm_tdm_sig_genes,blast!="Espn")
df.merge.counts_scale_mantle_dm_tdm_sig_genes<-dplyr::filter(df.merge.counts_scale_mantle_dm_tdm_sig_genes,blast!="elavl2")
#make gene name row names
rownames(df.merge.counts_scale_mantle_dm_tdm_sig_genes)<-df.merge.counts_scale_mantle_dm_tdm_sig_genes$blast
counts_scale_mantle_dm_tdm_sig_genes_rename<-df.merge.counts_scale_mantle_dm_tdm_sig_genes[,-1]
counts_scale_mantle_dm_tdm_sig_genes_rename<-counts_scale_mantle_dm_tdm_sig_genes_rename[,-13]



pheatmap(counts_scale_mantle_dm_tdm_sig_genes, annotation_col = column_annotation_samplephenotype)
pheatmap(counts_scale_mantle_dm_tdm_sig_genes_rename, annotation_col = column_annotation_samplephenotype)


