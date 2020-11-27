#Analysis of soil chemistry results for SF1
#Basic chemistry was analyzed by LSU AgCenter.  This is to test whether there was a difference among the treatments at the beginning of the experiment and how conditions changed over the course of the experiment.

#last updated Jan. 21, 2020 by Steve Formel

subtitle <- "Exploratory Analysis, soil_chem_analysis.R"

#import data-----

library(readxl)
df <- read_excel("Data/soil_chem/SF1_soil_chem_all_data.xls", 
                                    sheet = "data_cleaned")

#clean data----
library(plyr)

#rename chem titles
colnames(df) <- c("sampleID", "soil", "oil", "time", "Calcium", "Copper", "Magnesium", "pH", "Phosphorus", "Potassium",  "Sodium", "Sulfur", "Zinc", "Percent Carbon", "Percent Nitrogen")

df$soil <- revalue(df$soil, c("L"="Live", "S"="Autoclaved"))
df$oil <- revalue(df$oil, c("N"="No Oil", "Y"="Oil Added"))

df$time <- factor(df$time, levels = c("pre-exp", "post-exp"))
df$Treatment <- interaction(df$soil, df$oil)

#Post-Experiment----

df.post <- subset(df, df$time=="post-exp")

#PCA and plot post-experiment only as biplot
#from https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff

df.PCA <- prcomp(df.post[c(5:15)], center = TRUE, scale = TRUE)

library(ggfortify)

autoplot(df.PCA, 
         data = df.post,
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         loadings.label = TRUE,
         loadings.label.vjust = -1.5,
         loadings.label.hjust = -0.05,
         loadings.label.size = 3) +
  geom_point(color = "black",
             shape = 21,
             size = 7, 
             aes(fill = Treatment)) +
  scale_fill_manual(values = c("#000000", "#009E73", "#e79f00", "#9ad0f3")) +
  geom_abline(intercept = -0.1, slope = -0.6, color="gray", 
              linetype="dashed", size=0.5) +
  geom_abline(intercept = 0.1, slope = 1.5, color="gray", 
              linetype="dashed", size=0.5) +
  theme_bw() +
  labs(title = "SF1 Post Experiment Soil Chem", subtitle = subtitle)

#rearrange to put lines under labels
#p$layers <- p$layers[c(1,2,5,6,3,4)]

#save
ggsave(filename = "figs/soil_chem/exploratory/biplot_post_exp.png", width = 9, height = 7, units = "in")

#Stat Test for how the groups differed----

Calcium <- lm(Calcium ~ soil*oil, data = df.post) #sig for main effect of Oil
summary(Calcium)
hist(resid(Calcium))
shapiro.test(df.post$Calcium)

Copper<- lm(Copper ~ soil*oil, data = df.post) #Big time oil main effect
summary(Copper)
hist(resid(Copper))
shapiro.test(df.post$Copper) #failed

Magnesium <- lm(Magnesium ~ soil*oil, data = df.post) #sig for main effect of Soil
summary(Magnesium)
hist(resid(Magnesium))
shapiro.test(df.post$Magnesium)

pH <- lm(pH ~ soil*oil, data = df.post) #sig for main effect of Oil
summary(pH)
hist(resid(pH))
shapiro.test(df.post$pH)

Phosphorus <- lm(Phosphorus ~ soil*oil, data = df.post) #sig for main effect of Oil
summary(Phosphorus)
hist(resid(Phosphorus))
shapiro.test(df.post$Phosphorus)

Potassium <- lm(Potassium ~ soil*oil, data = df.post) #sig for main effect of Soil
summary(Potassium)
hist(resid(Potassium))
shapiro.test(df.post$Potassium)

Sodium <- lm(Sodium ~ soil*oil, data = df.post) #sig for main effect of Soil
summary(Sodium)
hist(resid(Sodium))
shapiro.test(df.post$Sodium)

Sulfur <- lm(Sulfur ~ soil*oil, data = df.post) #sig for main effect of Soil and Oil
summary(Sulfur)
hist(resid(Sulfur))
shapiro.test(df.post$Sulfur)

Zinc <- lm(Zinc ~ soil*oil, data = df.post) #sig for main effect of Oil
summary(Zinc)
hist(resid(Zinc))
shapiro.test(df.post$Zinc)

#Permanova-----
library(vegan)

dist.mat <- vegdist(df.post[,5:13], method = "euclidean")

adonis(dist.mat ~ soil*oil, data = df.post, method = "euclidean") #sig for main effects, interaction not quite sig


#Pre-experiment-----

df.pre <- subset(df, df$time=="pre-exp")

#PCA and plot pre-experiment only as biplot
#from https://towardsdatascience.com/principal-component-analysis-pca-101-using-r-361f4c53a9ff

df.PCA <- prcomp(df.pre[c(5:15)], center = TRUE, scale = TRUE)

#clean data
library(plyr)
df.pre$soil <- revalue(df.pre$soil, c("L"="Live", "S"="Autoclaved"))
df.pre$oil <- revalue(df.pre$oil, c("N"="No Oil", "Y"="Oil Added"))

df.pre$time <- factor(df.pre$time, levels = c("pre-exp", "post-exp"))
df.pre$Treatment <- interaction(df.pre$soil, df.pre$oil)



library(ggfortify)

autoplot(df.PCA, 
         data = df.pre,
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         loadings.label = TRUE,
         loadings.label.vjust = -1.5,
         loadings.label.hjust = -0.05,
         loadings.label.size = 3) +
  geom_point(color = "black",
             shape = 21,
             size = 7, 
             aes(fill = Treatment)) +
  scale_fill_manual(values = c("#000000", "#009E73", "#e79f00", "#9ad0f3")) +
  geom_abline(intercept = -0.1, slope = -0.6, color="gray", 
              linetype="dashed", size=0.5) +
  geom_abline(intercept = 0.1, slope = 1.5, color="gray", 
              linetype="dashed", size=0.5) +
  theme_bw() +
  labs(title = "SF1 Pre-experiment Soil Chem", subtitle = subtitle)

#save
ggsave(filename = "figs/soil_chem/exploratory/biplot_pre_exp.png", width = 9, height = 7, units = "in")


#Permanova-----
library(vegan)

dist.mat <- vegdist(df.pre[,5:13], method = "euclidean")

adonis(dist.mat ~ soil*oil, data = df.pre, method = "euclidean") #Not significant

#Pre to Post Exp Plot----

df.PCA <- prcomp(df[c(5:15)], center = TRUE, scale = TRUE)

#clean data
library(plyr)
df$soil <- revalue(df$soil, c("L"="Live", "S"="Autoclaved"))
df$oil <- revalue(df$oil, c("N"="No Oil", "Y"="Oil Added"))

df$time <- factor(df$time, levels = c("pre-exp", "post-exp"))
df$Treatment <- interaction(df$soil, df$oil)

#get group centroids
df$TT <- interaction(df$Treatment, df$time)
df.PCA.x <- as.data.frame(df.PCA$x)
df.PCA.x$groups <- df$TT
pca.centroids <- aggregate(df.PCA.x[,1:9], list(Treatment = df.PCA.x$groups), mean)

pca.centroids$time <- c(rep("pre-exp", 4), rep("post-exp", 4))
pca.centroids$time <- factor(pca.centroids$time, levels = c("pre-exp", "post-exp"))

pca.centroids$Treatment <- gsub(pattern = "\\.p.*[exp]",replacement = "",x = pca.centroids$Treatment)

library(ggfortify)

p <- autoplot(df.PCA, 
         data = df,
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label.colour = 'black',
         loadings.label = TRUE,
         loadings.label.vjust = 0.4,
         loadings.label.hjust = -0.3,
         loadings.label.size = 3)

#remove orig points
p$layers <- p$layers[c(2:3)]
  
#plot centroids
p1 <- p + 
  geom_point(data = pca.centroids, 
             color = "black",
             size = 7, 
             aes(x = PC1/5, 
                 y = PC2/5, 
                 shape = time, 
                 fill = Treatment)) +
  scale_fill_manual(values = c("#000000", "#009E73", "#e79f00", "#9ad0f3")) +
  scale_shape_manual(values = c(21,24)) +
  geom_abline(intercept = 0, slope = 2.3, color="gray", 
              linetype="dashed", size=0.5) +
  geom_abline(intercept = 0, slope = -0.5, color="gray", 
              linetype="dashed", size=0.5) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  labs(title = "SF1 Pre to Post Experiment Soil Chem", subtitle = subtitle)

#rearrange to put lines under labels
p1$layers <- p1$layers[c(3,1,4,5,2)]

p1 

#save
ggsave(filename = "figs/soil_chem/exploratory/biplot_pre_to_post_exp.png", width = 9, height = 7, units = "in")


#Make barplot-----
library(tidyverse)

df.traits <- df %>%
  gather(key = "Trait", value = "value", Calcium:'Percent Nitrogen')

library(Rmisc)

df.trait.means <- summarySE(data = df.traits, measurevar = "value", groupvars = c("Trait", "time", "soil","oil"))

ggplot(data = df.trait.means,
       aes(x= interaction(time, soil, oil, sep = " ; "),
           y = value,
           fill = interaction(soil, oil, sep = " ; "))) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.2) +
  geom_bar(stat = "identity", alpha = 0.5) +
  facet_wrap(~ Trait,
             scales = "free",
             ncol = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 50, b = 0, l = 0))) +
  labs(title = "SF1 Soil Chem",
       subtitle = subtitle,
       caption = "n for pre-exp = 2, n for post-exp = 3 ; for all elements unit is ppm",
       x = "Treatments",
       y = "",
       fill = "Treatments")

ggsave("figs/soil_chem/exploratory/barplot_all_traits.png", width = 14, height = 10, units = "in")

#ANOVA on chem traits-----

df.traits.post <- df[df$time=="post-exp",]

#replace spaces in colnames
colnames(df.traits.post) <- str_replace_all(colnames(df.traits.post), " ", ".")

A <- lm(Calcium ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Copper ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Magnesium ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Percent.Carbon ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Percent.Nitrogen ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(pH ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Phosphorus ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Potassium ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Sodium ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Sulfur ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)

A <- lm(Zinc ~ soil*oil, data = df.traits.post)
summary(A)
shapiro.test(A$residuals)
#Look into LDA, Factor analysis-----

#https://www.researchgate.net/post/What_is_the_difference_between_PCA_FA_and_LDA
#http://www.statpower.net/Content/312/Handout/Fabrigar1999.pdf

#After some more reading about EFA and LDA, it seems like they are just alternatives to PCA for characterizing the soil.  But it's pretty obvious that the soil breaks out into the treatments, and can be characterized by describing the C and N content, the pH, and the metal content