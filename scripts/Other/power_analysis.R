#power analysis for Steve Formel's 3rd chapter (Growth Chamber)

#April 16, 2018

#https://www.statmethods.net/stats/power.html

install.packages("pwr2")
library(pwr2)

#Cohen suggests (for one-way ANOVA) that f values of 0.1, 0.25, and 0.4 represent small, medium, and large effect sizes respectively.

#power test for balanced two-way ANOVA
pwr.2way(a=2,b=2,alpha=0.05, size.A = 12, size.B = 12, f.A = 0.4, f.B = 0.4)

#sample size calculation for balanced two-way ANOVA
ss.2way(a=2, b=2, alpha = 0.05, beta = 0.4, f.A = 0.25, f.B = 0.25, B = 100)


#so sample size for small effect size (0.1) and high power (0.9) = 101
#medium effect size (0.25) = 43
#large effect (0.4) = 17 

#for power = 0.8 and large effect, n = 13