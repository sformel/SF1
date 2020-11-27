#Make plot of DE genes from oil vs DE genes from soil for SF1 project


#by Steve Formel
#last updated: Oct 25, 2019
library(data.table)
library(tidyverse)


#import and combine all DE genes
filenames.oil <- list.files(path = "Data/DEseq_output_oil",full.names = TRUE, pattern = "*results$")

DE_oil <- do.call("rbind", lapply(filenames.oil, fread, data.table = FALSE))

filenames.soil <- list.files(path = "Data/DEseq_output_soil",full.names = TRUE, pattern = "*results$")

DE_soil <- do.call("rbind", lapply(filenames.soil, fread, data.table = FALSE))

length(DE_oil$V1[DE_oil$padj < 0.05])
length(DE_soil$V1[DE_soil$padj < 0.05])

int_gene <- intersect(DE_oil$V1[DE_oil$padj < 0.05], DE_soil$V1[DE_soil$padj < 0.05])

#make Venn DIagram of gene overlap
library(VennDiagram)

Oil <- DE_oil$V1[DE_oil$padj < 0.05]
Soil <- DE_soil$V1[DE_soil$padj < 0.05]

VD <- calculate.overlap(list(Oil, Soil))

grid.newpage();
venn.plot <- draw.pairwise.venn(area1 = 1023, 
                   area2 = 657,
                   cross.area = 21, 
                   category = c("Oil Treatment", "Autoclave Treatment"), 
                   fill = c("black", "red"),
                   alpha = c(0.5,0.5), 
                   cex = c(2,2,2), 
                   cat.cex = c(2,2), 
                   cat.pos = c(0,0), cat.dist = c(0.05,0.05))

png(file="DE_main_gene_venn.png")
grid.draw(venn.plot);
dev.off()



#Group by main effects of treatments------

DE_genes <- rbind(DE_oil[DE_oil$V1 %in% int_gene,] , DE_soil[DE_soil$V1 %in% int_gene,])

#clean things up----
colnames(DE_genes)[1] <- "Gene"
DE_genes$Gene <- str_sub(DE_genes$Gene, 9)

DE_genes$log2FoldChange <- round(DE_genes$log2FoldChange,2)

#Make plot----

#view by gene
DE_genes$COMP <- as.factor(paste(DE_genes$sampleA, DE_genes$sampleB, sep = "_"))
#levels(DE_genes$COMP)[1] <- "Unoiled_to_Oiled"
#levels(DE_genes$COMP)[2] <- "Live_to_Sterile"

#This is silly but rught now this is the only way my brain can see to do this
plot.data <- DE_genes %>%
  select(Gene,COMP, log2FoldChange) %>%
  spread(key = Gene, value = log2FoldChange) %>%
  t()

colnames(plot.data) <- plot.data[1,]
plot.data <- as.data.frame(plot.data[-1,])
plot.data$Gene <- row.names(plot.data)

#bar plots

ggplot(data = DE_genes, aes(x = Gene, y = -(log2FoldChange))) +
  geom_bar(aes(fill = COMP), stat = "identity", position = position_dodge(), alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("black", "red"), labels = c("Oil Treatment", "Autoclave Treatment")) +
  ylab("log2 Fold Change")

ggsave("DE_main_barplot.png", width = 5, height = 4, units = "in")

#point plots
plot.data$oil_N_oil_Y <- as.numeric(as.character(plot.data$oil_N_oil_Y))
plot.data$soil_L_soil_S <- as.numeric(as.character(plot.data$soil_L_soil_S))

ggplot(data = plot.data, aes(x = -(oil_N_oil_Y), y = -(soil_L_soil_S), color = Gene)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Oil Treatment") +
  ylab("Autoclave Treatment") +
  annotate("text", x = 10, y = 15, label = "Oil and Autoclave \nincrease expression") +
  annotate("text", x = 10, y = -15, label = "Oil increases and \nAutoclave decreases \nexpression") +
  annotate("text", x = -10, y = 15, label = "Oil decreases and \nAutoclave increases \nexpression") +
  annotate("text", x = -10, y = -15, label = "Oil and Autoclave \ndecrease expression") +
  theme(legend.position = "none")

ggsave("Genes_by_main_trt.png", width = 5, height = 4, units = "in")

#barplot diagram to compare relative change in treatments----

#import and combine all DE genes
filenames <- list.files(path = "Data/SF1_Trinity_24Mar2019_output/DEseq_output",full.names = TRUE, pattern = "*2.DE_results$")

DE <- do.call("rbind", lapply(filenames, fread, data.table = FALSE))

#trim gene name
DE$Gene <- str_sub(DE$V1, 9)
#Group by main effects of treatments------

#DE_all <- DE[DE$Gene %in% plot.data$Gene,]

#clean things up----
DE_all$log2FoldChange <- round(DE_all$log2FoldChange,2)

#Make plot----

#view by gene
DE_all$COMP <- as.factor(paste(DE_all$sampleA, DE_all$sampleB, sep = "|"))
levels(DE_all$COMP)[1] <- "LNLY"
levels(DE_all$COMP)[2] <- "LNSN"
levels(DE_all$COMP)[3] <- "LNSY"
levels(DE_all$COMP)[4] <- "LYSN"
levels(DE_all$COMP)[5] <- "LYSY"
levels(DE_all$COMP)[6] <- "SNSY"

#order COMP
#DE_all <- subset(DE_all, COMP =="LNLY" | COMP =="LNSN" | COMP =="LNSY")

ggplot(data = DE_all, aes(x = Gene, y = log2FoldChange, fill = COMP)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45))

#Make PCA of pairwise comparisons----
#I think if I make a PCA of pairwise comparisons based on DE genes, I can find which pairwise comparisons are the closest and most different

DE$COMP <- as.factor(paste(DE$sampleA, DE$sampleB, sep = "|"))
levels(DE$COMP)[1] <- "LNLY"
levels(DE$COMP)[2] <- "LNSN"
levels(DE$COMP)[3] <- "LNSY"
levels(DE$COMP)[4] <- "LYSN"
levels(DE$COMP)[5] <- "LYSY"
levels(DE$COMP)[6] <- "SNSY"

FX <- DE %>%
  select(Gene, log2FoldChange, COMP)

#This is silly but rught now this is the only way my brain can see to do this
FX.data <- DE %>%
  select(Gene,COMP, log2FoldChange) %>%
  spread(key = Gene, value = log2FoldChange) %>%
  t()

colnames(FX.data) <- FX.data[1,]
FX.data <- FX.data[-1,]
FX.data <- na.omit(FX.data)
FX.data <- t(FX.data)
class(FX.data) <- "numeric"
FX.data <- as.data.frame(FX.data)

#negate some signs to see what happens
FX.data.neg <- -(FX.data)
rownames(FX.data.neg) <- paste0("N_",rownames(FX.data.neg))
FX.data.both <- rbind(FX.data, FX.data.neg)

fx.pca <- prcomp(FX.data.both)

df_out <- as.data.frame(fx.pca$x)
df_out$COMP <- row.names(df_out)

ggplot(data = df_out, aes(x = PC1, y = PC2, color = COMP)) +
  geom_text(aes(label = COMP)) +
  annotate("text", x = -20, y = 0, label = "PC1 explains the effect of the oil on the live trt") +
  annotate("text", x = -20, y = 100, label = "PC2 explains the effect of the oil on the sterile trt") +
  annotate("text", x = -20, y = -150, label = "A combination of the two to the downward right \n explains the difference in soil \n combined with the difference in oil")

ggsave("PCA_pairwise_comp_DE_NEG.png")

#Can I break down the principal components to find the names of these genes?  Essentially it looks like oil did not effect gene expression in the sterile treatments nearly as extremely as it affected gene expression in the live treatments.  This is assumed because PC1 is stronger (supposedly) than PC2, so the distance counts for less.  But I need to check this out.  But after thinking about it, the direction of the pairwise comparison is really important, so I may need to flip some signs.

