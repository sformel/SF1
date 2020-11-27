library(Rmisc)
library(readxl)
SF1_GC_data <- read_excel("Data/SF1_GC_data.xlsx", 
                          sheet = "biomass", col_types = c("text", 
                                                           "text", "text", "text", "numeric"), na = "NA")

SF1_GC_data <- na.omit(SF1_GC_data)
summary(aov(formula = dry_mass ~ soil*oil, data = subset(SF1_GC_data,Tissue %in% "AG")))
summary(aov(formula = dry_mass ~ soil*oil, data = subset(SF1_GC_data,Tissue %in% "IR")))
summary(aov(formula = dry_mass ~ soil*oil, data = subset(SF1_GC_data,Tissue %in% "OT")))

total_BG_dry_mass <- aggregate(dry_mass ~ plantID + soil + oil, data = SF1_GC_data[SF1_GC_data$Tissue=="OT" | SF1_GC_data$Tissue=="OT",], sum)
total_dry_mass <- aggregate(dry_mass ~ plantID + soil + oil, SF1_GC_data, sum)

summary(aov(formula = dry_mass ~ soil*oil, data = total_BG_dry_mass))
summary(aov(formula = dry_mass ~ soil*oil, data = total_dry_mass))

summarySE(data = total_dry_mass, measurevar = "dry_mass", groupvars = c("oil", "soil"))

#reorder treatments
SF1_GC_data$Treatment <- interaction(SF1_GC_data$soil, SF1_GC_data$oil)

#rename
library(plyr)
SF1_GC_data$Treatment <- revalue(SF1_GC_data$Treatment, c("L.N"="Live, No Oil", "S.N"="Sterile, No Oil", "L.Y" = "Live, Oiled", "S.Y" = "Sterile, Oiled"))

#reorder levels
SF1_GC_data$Treatment <- factor(SF1_GC_data$Treatment, levels = c("Live, No Oil", "Live, Oiled", "Sterile, No Oil", "Sterile, Oiled"))

#a little unexpected, soil and oil are significant, but no interactoins are.

library(Rmisc)
summarySE(data = SF1_GC_data, measurevar = "dry_mass", groupvars = c("Tissue","oil", "soil"))



#plot it---------
library(ggplot2)

ggplot(data = SF1_GC_data, aes(x = oil, y = dry_mass, fill = Treatment)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(~ Tissue) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black" ))

#ggsave(filename = mass_by_tissue.png, width = 8, height = 6, units = "in")

#sum all the mass
ggplot(data = total_dry_mass, aes(x = oil, y = dry_mass, fill = interaction(oil, soil))) +
  geom_boxplot(alpha = 0.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#F0E442","black" ))

#ggsave(filename = total_mass.png, width = 8, height = 6, units = "in")

total_dry_mass$soil <- factor(total_dry_mass$soil)
total_dry_mass$oil <- factor(total_dry_mass$oil)

total_dry_mass$interaction <- interaction(total_dry_mass$soil, total_dry_mass$oil)

t.test(total_dry_mass$dry_mass ~ total_dry_mass$soil, var.equal = FALSE)
t.test(total_dry_mass$dry_mass ~ total_dry_mass$oil, var.equal = FALSE)

#A <- aov(formula = dry_mass ~ interaction, data = total_dry_mass)

#to do type iii, need the car package
library(car)

Anova(lm(dry_mass ~ interaction, data = total_dry_mass, contrasts=list(soil=contr.sum, oil=contr.sum)), type=3)

#TukeyHSD(x = A, which = "interaction", conf.level = 0.95)

