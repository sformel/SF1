#Analysis of soil chemistry results for SF1
#Basic chemistry was analyzed by LSU AgCenter.  This is to test whether there was a difference among the treatments at the beginning of the experiment and how conditions changed over the course of the experiment.

#last updated Oct 28, 2020 by Steve Formel

subtitle <- "Exploratory Analysis, soil_chem_analysis.R"

#import data and load libraries-----
library(ggfortify)
library(readxl)
library(plyr)
library(vegan)
library(tidyverse)
library(compositions)
library(ggplot2)

df <- read_excel("Data/soil_chem/SF1_soil_chem_all_data.xls", 
                                    sheet = "data_cleaned")

#clean data----

#rename chem titles
colnames(df) <- c("sampleID", "soil", "oil", "time", "Calcium", "Copper", "Magnesium", "pH", "Phosphorus", "Potassium",  "Sodium", "Sulfur", "Zinc", "Percent Carbon", "Percent Nitrogen")

df$soil <- plyr:::revalue(df$soil, c("L"="Live", "S"="Autoclaved"))
df$oil <- plyr:::revalue(df$oil, c("N"="No Oil", "Y"="Oil Added"))

df$time <- factor(df$time, levels = c("pre-exp", "post-exp"))
df$Treatment <- interaction(df$soil, df$oil)

SF1_aes <- aes(shape = soil,
      fill = oil)

my.scale_aes <- list(scale_shape_manual(values = c(21,24)),
    scale_fill_manual(values = c("white", "black")),
    theme_bw(),
    guides(fill = guide_legend(override.aes=list(shape=21, size = 4)),
           shape = guide_legend(override.aes=list(size = 4))),
    labs(shape = "Soil",
         fill = "Oil"))

#Pre-experiment-----

df.pre <- subset(df, df$time=="pre-exp")

#PCA and plot pre-experiment only as biplot
#from https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff

df.PCA <- prcomp(df.pre[c(5:15)], center = TRUE, scale = TRUE)

p <- autoplot(df.PCA, 
         data = df.pre,
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         loadings.label = TRUE,
         loadings.label.vjust = -1.5,
         loadings.label.hjust = -0.05,
         loadings.label.size = 3)

p$layers <- p$layers[-1]

pre <- p+ geom_point(mapping = SF1_aes) +
  my.scale_aes +
  geom_abline(intercept = -0.1, slope = -0.7, color="gray", 
              linetype="dashed", size=0.5)

ggsave(filename = "figs/manuscript/soil_chem/SF1_biplot_pre_exp_V1.png", width = 8, height = 6, units = "in")

#Post-Experiment----

df.post <- subset(df, df$time=="post-exp")

#PCA and plot post-experiment only as biplot
#from https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff

df.PCA <- prcomp(df.post[c(5:15)], center = TRUE, scale = TRUE)


p <- autoplot(df.PCA, 
         data = df.post,
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         loadings.label = TRUE,
         loadings.label.vjust = -1.5,
         loadings.label.hjust = -0.05,
         loadings.label.size = 3)
  
p$layers <- p$layers[-1] 
  
post <- p + geom_point(mapping = SF1_aes) +
  my.scale_aes +
  geom_abline(intercept = -0.1, slope = -0.6, color="gray", 
              linetype="dashed", size=0.5) +
  geom_abline(intercept = 0.1, slope = 1.5, color="gray", 
              linetype="dashed", size=0.5)

#save
ggsave(filename = "figs/manuscript/soil_chem/SF1_biplot_post_exp_V1.png", width = 8, height = 6, units = "in")


#Plot pre and post together-----

library(cowplot)

top <- plot_grid(pre + theme(legend.position = "none"), post + theme(legend.position = "none"), ncol = 2, labels = "AUTO")

plot_grid(top, get_legend(post + theme(legend.position = "bottom")), ncol = 1, rel_heights = c(1,0.2))

ggsave("figs/manuscript/soil_chem/SF1_pre_post_chem_V1.png", width = 11, height = 6, units = "in")

#Aitchison distance of chemistry-----


#CLR Transformations

#This transforms each number into it's relative abundance in that sample
comps <- acomp(df[,c(5:7,9:13)])

RA.df <- cbind(comps, df)

RA.df.gathered <- RA.df %>%
  gather(key = "Element", value = "Proportion", c(1:8))

RA.df.gathered %>%
  ggplot(aes(x = Element,
             y = Proportion,
             color = Element)) +
  geom_boxplot() +
  facet_wrap(~ time*soil*oil)
  
#But as I understand it, it's inappropriate to do the stats on this because it's closed.  Need to translate to CLR.

# 7 Multivariate Methods
# The central idea of the package ??? following the coordinate approach of [Pawlowsky-Glahn(2003)]
# and [Pawlowsky-Glahn and Mateu-Figueras(2005)] ??? is to transform the data by one
# of transforms into a classical multivariate dataset, to apply classical multivariate
# statistics and to back transform or interpreted the results afterwards in the original
# space

#According to https://www.rdocumentation.org/packages/rgr/versions/1.1.15/topics/clr all units should be the same, so the percents should be converted to ppm (or whatever is used).  I'm not sure about pH, but this paper (10.1007/s12665-019-8248-6) just used clr transformed pH , but its not clear what they used as the geometric mean.  So maybe the thing to do is CLR transform the other chemicals and then include the pH in the PCA.

#Convert %C and N to ppm
df$`Percent Carbon` <- df$`Percent Carbon`*10000
df$`Percent Nitrogen` <- df$`Percent Nitrogen`*10000

clr.df <- clr(df[,c(5:7,9:15)])

clr.df <- as.data.frame(clr.df)
clr.df$pH <- df$pH

pc <- princomp(x = clr.df)

pc$Loadings # The loadings as compositional vector
pc$loadings # The loadings in clr-space
df.pca <- pc$scores

aitchison.pca.plot <- cbind(na.omit(df), df.pca) %>%
  ggplot(aes(x = Comp.1,
             y = Comp.2)) +
  geom_point(mapping = SF1_aes) +
    facet_wrap(~ time) +
    my.scale_aes

p <- autoplot(pc, 
         data = df,
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         loadings.label = TRUE,
         loadings.label.vjust = -1.5,
         loadings.label.hjust = -0.05,
         loadings.label.size = 3)
  
p$layers <- p$layers[-1] 

#set aspect ratio for ploitting labels

p + geom_point(mapping = SF1_aes) +
  my.scale_aes +
  geom_abline(intercept = -0.075, slope = -0.6, color="gray", 
              linetype="dashed", size=0.5) +
  geom_abline(intercept = 0.1, slope = 1.5, color="gray", 
              linetype="dashed", size=0.5) +
  coord_fixed() +
  geom_text(data=data.frame(x = 0.15,y = 0.36), aes(x, y), label="Oiled", hjust = 0.5, vjust=0.5, angle = 58) +
  geom_text(data=data.frame(x = 0.19,y = 0.35), aes(x, y), label="Not Oiled", hjust = 0.5, vjust=0.5, angle = 58) +
  geom_text(data=data.frame(x = 0.32,y = -0.295), aes(x, y), label="Pre-Experiment", hjust = 0.5, vjust=0.5, angle = -30) +
  geom_text(data=data.frame(x = 0.35,y = -0.265), aes(x, y), label="Post-Experiment", hjust = 0.5, vjust=0.5, angle = -30)

ggsave("figs/manuscript/soil_chem/SF1_biplot_V2.png", width = 8, height = 6, units = "in")

#PERMANOVA----
df.clr <- cbind(df[,c(1:4)], clr.df)

set.seed(1)
adonis(formula = clr.df ~ soil*oil*time, data = df.clr, permutations = 9999, method = "euclidean")

#Boxplot-----

elements <- df.clr %>%
  gather(key = "Element", value = "clr_val", Calcium:'Percent Nitrogen') %>%
  ggplot(aes(x = fct_reorder(Element, clr_val, .fun='median'),
         y = clr_val)) +
  geom_point(mapping = SF1_aes,
              size = 3,
             position = position_dodge(width = 0.5)) +
  my.scale_aes +
  theme_bw() +
  labs(y = "Centered Log-Ratio",
       x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ time)

ggsave("figs/manuscript/soil_chem/SF1_clr_stripplot_V1.png", width = 8, height = 6, units = "in")

pH <- df.clr %>%
  ggplot(aes(x = time,
         y = pH)) +
  geom_point(mapping = SF1_aes,
              size = 3,
             position = position_dodge(width = 0.5)) +
  my.scale_aes +
  theme_bw() +
  labs(y = "pH",
       x = NULL)

ggsave("figs/manuscript/soil_chem/SF1_pH_stripplot_V1.png", width = 8, height = 6, units = "in")

#plot Fig S1

library(cowplot)

top <- plot_grid(elements + theme(legend.position = "None"), pH + theme(legend.position = "None"), rel_widths = c(1,0.5), labels = "AUTO")

plot_grid(top, get_legend(pH + theme(legend.position = "bottom")), nrow = 2, rel_heights = c(1,0.1))

ggsave("figs/manuscript/soil_chem/SF1_soil_chem_strip_V1.png", width = 8, height = 6, units = "in")


#OLD BELOW------

#Look into LDA, Factor analysis-----

#https://www.researchgate.net/post/What_is_the_difference_between_PCA_FA_and_LDA
#http://www.statpower.net/Content/312/Handout/Fabrigar1999.pdf

#After some more reading about EFA and LDA, it seems like they are just alternatives to PCA for characterizing the soil.  But it's pretty obvious that the soil breaks out into the treatments, and can be characterized by describing the C and N content, the pH, and the metal content
         