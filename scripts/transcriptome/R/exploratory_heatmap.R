#SF1 script
#last updated by Steve Formel
#Jan 22, 2020

#Based on script received from Sarah Khalil on Jan 22, 2020

#make metadata for figures
substitle = "Exploratory Analysis ; exploratory_heatmap.R"

#load packages-----

library(DESeq2)
library(tidyverse)
library(data.table)
library(pheatmap)

#Examples of pheatmap------

# test = matrix(rnorm(200), 20, 10)
# test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
# test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
# test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
# colnames(test) = paste("Test", 1:10, sep = "")
# rownames(test) = paste("Gene", 1:20, sep = "")
# 
# pheatmap(test)
# pheatmap(test, kmeans_k = 2)
# pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
# pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# pheatmap(test, cluster_row = FALSE)
# pheatmap(test, legend = FALSE)
# annotation_col = data.frame(
#   CellType = factor(rep(c("CT1", "CT2"), 5)), 
#   Time = 1:5
# )
# rownames(annotation_col) = paste("Test", 1:10, sep = "")
# 
# annotation_row = data.frame(
#   GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
# )
# rownames(annotation_row) = paste("Gene", 1:20, sep = "")
# 
# pheatmap(test, annotation_col = annotation_col)

#import data----

#it turns out that salmon estimates countso my numbers have decimals.  Therefore I have to import differently than Sarah because I have to account for this using tximport example on: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#Also see:https://groups.google.com/forum/#!topic/trinityrnaseq-users/8yPuoZh1E5I actually scaled TPM

#I'm having problems understanding the tx2gene file format.  Note that there are two outputs from salmon: isoform and gene level quant files, I can use gene if I'm having problems

library("tximport")
library("readr")

dir <- file.path("data/transcriptome/SF1_Trinity_24Mar2019_output/salmon_output")
samdf <- read.csv("data/transcriptome/DEseq_sample_list_for_R.csv", row.names = 1)

files <- file.path(dir, samdf$name, "quant.sf")
names(files) <- samdf$name

#make DEseq data set (only need to do once)----

#import gene map
#tx2gene <- read_tsv(file.path(dir, "../Trinity.fasta.gene_trans_map.cleaned"), col_names = FALSE)
#colnames(tx2gene) <- c("GeneID", "isoformID")

#flipflop isoform and gene- important that whatever is being analyzed is first
#tx2gene <- tx2gene %>%
#  select(isoformID, GeneID)

#txi <- tximport(files, 
#                type="salmon", 
#                tx2gene=tx2gene,
#                txOut = TRUE)  #a two-column data.frame linking transcript id (column 1) to gene id (column 2). the column names are not relevant, but this column order must be used. this argument is required for gene-level summarization for methods that provides transcript-level estimates only (kallisto, Salmon, Sailfish), look here for options: https://rdrr.io/bioc/tximport/man/tximport.html

#import into DEseq format - took about 35 min to run at SF desktop
#df_dds <- ddsTxi <- DESeqDataSetFromTximport(txi,
#                                             colData = samdf,
#                                             design = ~ oil*soil)

#start <- Sys.time()

#Run DEseq
#dds<-DESeq(df_dds)

#end <- Sys.time()

#save as RDS

#saveRDS(dds, file = "data/transcriptome/DEseq_df_isoform.rds")

#load DEseq data RDS-----
dds <- readRDS("data/transcriptome/DEseq_df_isoform.rds")

#Variance Stabilizing Transformation https://rdrr.io/bioc/DESeq2/man/varianceStabilizingTransformation.html
dds.vsd <- vst(dds, blind=TRUE)

counts<-(assay(dds.vsd)) #get count data per gene - this does not have any negatives, so i will z-transform to get that
counts_scale<-scale(counts) #z-tranform count data


#import DEseq analysis done in Trinity-----
LNLY <- read.delim(file.path(dir, "../DEseq_output/salmon.isoform.counts.matrix.soil_L_oil_N_vs_soil_L_oil_Y.DESeq2.DE_results"))
LNLY$isoformID <- rownames(LNLY)
rownames(LNLY) <- c()

LNSN <- read.delim(file.path(dir, "../DEseq_output/salmon.isoform.counts.matrix.soil_L_oil_N_vs_soil_S_oil_N.DESeq2.DE_results"))
LNSN$isoformID <- rownames(LNSN)
rownames(LNSN) <- c()

LNSY <- read.delim(file.path(dir, "../DEseq_output/salmon.isoform.counts.matrix.soil_L_oil_N_vs_soil_S_oil_Y.DESeq2.DE_results"))
LNSY$isoformID <- rownames(LNSY)
rownames(LNSY) <- c()

LYSN <- read.delim(file.path(dir, "../DEseq_output/salmon.isoform.counts.matrix.soil_L_oil_Y_vs_soil_S_oil_N.DESeq2.DE_results"))
LYSN$isoformID <- rownames(LYSN)
rownames(LYSN) <- c()

LYSY <- read.delim(file.path(dir, "../DEseq_output/salmon.isoform.counts.matrix.soil_L_oil_Y_vs_soil_S_oil_Y.DESeq2.DE_results"))
LYSY$isoformID <- rownames(LYSY)
rownames(LYSY) <- c()

SNSY <- read.delim(file.path(dir, "../DEseq_output/salmon.isoform.counts.matrix.soil_S_oil_N_vs_soil_S_oil_Y.DESeq2.DE_results"))
SNSY$isoformID <- rownames(SNSY)
rownames(SNSY) <- c()


#Oil addition DGE----

oil.live <- LNLY
oil.autoclave <- SNSY

#transform p value for human-readable numbers
oil.live$Lpadj<- -log10(oil.live$padj) #create -log10 (FDR)  anything with padj = 0.05 = 1.3 ; 0.1 = 1 (DEseq standard)
oil.live <- na.omit(oil.live)

oil.autoclave$Lpadj<- -log10(oil.autoclave$padj) #create -log10 (FDR)
oil.autoclave <- na.omit(oil.autoclave)

oil.live.DE <- subset(oil.live, Lpadj>1.3) 
oil.autoclave.DE <- subset(oil.autoclave, Lpadj>1.3) 

#Subset to isoforms that are in both live and autoclaved
oil.DE <- intersect(oil.live.DE$isoformID, oil.autoclave.DE$isoformID)

#subset count data to oil gene isoforms
counts_scale.genes_that_go_up_with_oil <- subset(counts_scale, rownames(counts) %in% oil.DE)

#make color scheme for treatments
annotation_colors = list(
  oil = c(Oiled="black", Not_Oiled="gray"),
  soil = c(Autoclaved="burlywood4", Live="darkgreen"))

pheatmap(counts_scale.genes_that_go_up_with_oil,
         show_rownames=F,
         annotation_col= samdf[,2:3],
         annotation_colors = annotation_colors,
         main = "SF1 DEgenes when oil is present \n exploratory_heatmap.R",
         filename = "figs/transcriptome/exploratory/heatmap_oil_added_genes.png", width = 9, height = 7)

#autoclave DGE----

autoclave.no.oil <- LNSN
autoclave.oil.added <- LYSY

#transform p value for human-readable numbers
autoclave.no.oil$Lpadj<- -log10(autoclave.no.oil$padj) #create -log10 (FDR)  anything with padj = 0.05 = 1.3 ; 0.1 = 1 (DEseq standard)
autoclave.no.oil <- na.omit(autoclave.no.oil)

autoclave.oil.added$Lpadj<- -log10(autoclave.oil.added$padj) #create -log10 (FDR)
autoclave.oil.added <- na.omit(autoclave.oil.added)

autoclave.no.oil.DE <- subset(autoclave.no.oil, Lpadj>1.3) 
autoclave.oil.added.DE <- subset(autoclave.oil.added, Lpadj>1.3) 

#Subset to isoforms that are in both live and autoclaved
autoclave.DE <- intersect(autoclave.no.oil.DE$isoformID, autoclave.oil.added.DE$isoformID)

#subset count data to autoclave gene isoforms
counts_scale.genes_that_go_up_with_autoclave <- subset(counts_scale, rownames(counts) %in% autoclave.DE)

#make color scheme for treatments
annotation_colors = list(
  oil = c(Oiled="black", Not_Oiled="gray"),
  soil = c(Autoclaved="burlywood4", Live="darkgreen"))

pheatmap(counts_scale.genes_that_go_up_with_autoclave,
         show_rownames=F,
         annotation_col= samdf[,2:3],
         annotation_colors = annotation_colors,
         main = "SF1 DEgenes when autoclave is present \n exploratory_heatmap.R",
         filename = "figs/transcriptome/exploratory/heatmap_autoclave_added_genes.png", width = 9, height = 7)

#Both treatments DGE----

#Subset to isoforms that are in both in oil addition and autoclaving lists, above
both.DE <- intersect(autoclave.DE, oil.DE)

#subset count data to autoclave gene isoforms
counts_scale.genes_that_go_up_with_autoclave <- subset(counts_scale, rownames(counts) %in% both.DE)

#make color scheme for treatments
annotation_colors = list(
  oil = c(Oiled="black", Not_Oiled="gray"),
  soil = c(Autoclaved="burlywood4", Live="darkgreen"))

pheatmap(counts_scale.genes_that_go_up_with_autoclave,
         show_rownames=F,
         annotation_col= samdf[,2:3],
         annotation_colors = annotation_colors,
         main = "SF1 DEgenes when both treatments are considered \n exploratory_heatmap.R",
         filename = "figs/transcriptome/exploratory/heatmap_both_trt_genes.png", width = 9, height = 7)
