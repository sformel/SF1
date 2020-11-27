#SF1 by Steve Formel
#Last updated: Oct 24, 2019

#Biomass
library(readxl)
library(tidyverse)
library(ggplot2)


biomass <- read_excel("Data/SF1_biomass.xlsx", 
                          sheet = "data")

biomass$plant_mass_g <- as.numeric(biomass$plant_mass_g)

hist(biomass$plant_mass_g)
shapiro.test(biomass$plant_mass_g)

biomass$log_mass <- log(biomass$plant_mass_g)

hist(biomass$log_mass)
shapiro.test(biomass$log_mass)

#totally not normal

ggplot(data = na.omit(biomass), aes(x = oil, y = plant_mass_g, fill = type)) +
  geom_boxplot() +
  facet_grid(~ ifelse(soil=="L", "Live Soil", "Sterile Soil")) +
  xlab("Oil Added")

#ANOVA
summary(aov(plant_mass_g ~ oil*soil, data = subset(biomass, type=="AG")))
summary(aov(plant_mass_g ~ oil*soil, data = subset(biomass, type=="IR")))        
summary(aov(plant_mass_g ~ oil*soil, data = subset(biomass, type=="OT")))        
summary(aov(plant_mass_g ~ oil*soil, data = biomass))        

#KW
kruskal.test(plant_mass_g ~ oil, data = subset(biomass, type=="AG")) #p = 0.03
kruskal.test(plant_mass_g ~ soil, data = subset(biomass, type=="AG")) #p = 0.0002

kruskal.test(plant_mass_g ~ oil, data = subset(biomass, type=="IR")) #p = 0.024
kruskal.test(plant_mass_g ~ soil, data = subset(biomass, type=="IR")) #p = 0.33

kruskal.test(plant_mass_g ~ oil, data = subset(biomass, type=="OT")) #p = 0.07
kruskal.test(plant_mass_g ~ soil, data = subset(biomass, type=="OT")) #p = 0.08

kruskal.test(plant_mass_g ~ oil, data = biomass) #p = 0.035
kruskal.test(plant_mass_g ~ soil, data = biomass) #p = 0.032

#oil only
ggplot(data = na.omit(biomass), aes(x = oil, y = plant_mass_g, fill = soil)) +
  geom_boxplot() +
  xlab("Oil Added")

