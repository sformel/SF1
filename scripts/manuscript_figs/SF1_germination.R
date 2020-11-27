#SF1 by Steve Formel
#Last updated: Oct 30, 2020

# Post-experiment plant Biomass to traits relationship for SF1 project

#load libraries----
library(readxl)
library(tidyverse)
library(ggplot2)

#load and clean data----
df <- read_excel("data/plant/SF1_GC_data.xlsx", 
                          sheet = "germination")

df %>%
  ggplot(aes(x = day_number,
           y = 12-not_germinated,
           color = interaction(soil,oil))) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", formula=y~x-1 ) +
  facet_wrap(~soil*oil)

#filter to final week
df %>%
  filter(day_number==28) %>%
  ggplot(aes(x = soil,
           y = 12-not_germinated,
           color = oil)) +
  geom_boxplot()

#Poisson Regressoin

df.lastday <- df %>%
  filter(day_number==28) 

#However, underdispersed
ggplot(df.lastday, aes(12-not_germinated, fill = interaction(soil,oil))) +
  geom_histogram(binwidth=.5, position="dodge")

summary(m1 <- glm(12-not_germinated ~ soil*oil, family="poisson", data=df.lastday))

AER:::dispersiontest(object = m1, trafo = NULL, alternative = c("less"))
