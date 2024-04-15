## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## small, species-specific effects on pollen chemistry and plant growth across
  ## flowering plant species

# script to plot and analyze greenhouse CO2 logger data

library("tidyverse")

# Experiment 1
C1<-read.csv("CO2_log_exp1.csv")

r1 <- subset(C1, datenum > 6000 & datenum < 8000)
r2 <- subset(C1, datenum > 9500 & datenum < 11000)

par(mfcol = c(1,2))
boxplot(as.numeric(r1$co2_1), as.numeric(r1$co2_2), 
        outline = FALSE, col= c('#becbca', '#69a1ff'), 
        axes=FALSE, xlab="", ylab = "CO2 (ppm)", ylim = c(350,800))
axis(1, labels = c('ambient', 'elevated'), at = c(1,2))
axis(2)
title('Round 1')

boxplot(as.numeric(r2$co2_2), as.numeric(r2$co2_1), 
        outline = FALSE, col= c('#becbca', '#69a1ff'), 
        axes=FALSE, ann = FALSE, ylim = c(350,800))
axis(1, labels = c('ambient', 'elevated'), at = c(1,2))
title('Round 2')

dev.off()

#Experiment 2
#read in data from experiment 2
C2<-read.csv("CO2_log_exp2.csv")

# separate data by round
r1 <- C2 %>% filter(!is.na(Round.1)) # round 1 goes from Nov 17 2021 - Mar 18 2022
r2 <- C2 %>% filter(!is.na(Round.2)) # round 2 goes from Jan 20 2022 until end. 

# remove outliers as determined by +/- 1.5 IQR 
C2 <- C2 %>% filter(is.na(outliers))

# plot CO2 levels 
par(mfcol = c(1,2))
boxplot(as.numeric(r1$aco2), as.numeric(r1$eco2), 
        outline = FALSE, col= c('#becbca', '#69a1ff'), 
        axes=FALSE, xlab="", ylab = "CO2 (ppm)", 
        ylim = c(300,800))
axis(1, labels = c('ambient', 'elevated'), at = c(1,2))
axis(2)
title('Round 1')

boxplot(as.numeric(r2$aco2), as.numeric(r2$eco2), 
        outline = FALSE, col= c('#becbca', '#69a1ff'), 
        axes=FALSE, ann = FALSE, 
        ylim = c(300,800))
axis(1, labels = c('ambient', 'elevated'), at = c(1,2))
title('Round 2')

dev.off()

