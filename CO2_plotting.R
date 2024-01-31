## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## inconsistent effects on pollen chemistry and plant growth across
  ## flowering plant species

# script to plot and analyze greenhouse CO2 logger data

library("tidyverse")
library("ggplot2")
library("chron")
library("jsonlite")
library("lubridate")

# Experiment 1
#Create pre-allocated data frame
interval <- 50 #Skipping interval for files: i.e., '10' loads and saves every tenth file within a folder
nrows <- round(220*(1440/interval))

#Make empty data frame
comb_data <- as.data.frame(matrix(ncol = 5, nrow = nrows))
colnames(comb_data) <- c('time', 'co2_a', 'co2_b', 'temp_a', 'temp_b')

#Set iterator
itr <- 1

#home directory
hd <- '/Volumes/Harvard Greenhouse/greenhouse/sensor_logs'
setwd(hd)
dirs <- dir()

for(j in dirs){
  setwd(file.path(hd, j))
  files <- list.files(pattern=".json")
  
  i <- 10
  
  for(i in seq(1,length(files), interval)){
    
    
    fname <-file.path(files[i])
    try({
      data<-fromJSON(fname)
      if('co2b' %in% names(data)){
        comb_data$time[itr] <- data$timestamp
        comb_data$co2_a[itr] <- data$co2a$co2
        comb_data$temp_a[itr] <- data$tempa
        comb_data$co2_b[itr] <- data$co2b$co2
        comb_data$temp_b[itr] <- data$tempb
        itr <- itr + 1
      }
    }, silent = TRUE)
    
    
    
    
    
    
  }
  print(j)
}


### Add timestamp

start.date <- parse_date_time('2020-01-01 00:00:00', "%Y-%m-%d %H:%M:%S") #Create a new variable for quantitative time, counted in days since October 1, 2022
comb_data$datenum <- as.numeric(difftime(parse_date_time(comb_data$time, 'ymd_HMS'), start.date))

r1 <- subset(comb_data, datenum > 6000 & datenum < 8000)
r2 <- subset(comb_data, datenum > 9500 & datenum < 11000)

par(mfcol = c(1,2))
boxplot(as.numeric(r1$co2_a), as.numeric(r1$co2_b), 
        outline = FALSE, col= c('#becbca', '#69a1ff'), 
        axes=FALSE, xlab="", ylab = "CO2 (ppm)", ylim = c(350,700))
axis(1, labels = c('ambient', 'elevated'), at = c(1,2))
axis(2)
title('Round 1')

boxplot(as.numeric(r2$co2_a), as.numeric(r2$co2_b), 
        outline = FALSE, col= c('#69a1ff', '#becbca'), 
        axes=FALSE, ann = FALSE, ylim = c(350,700))
axis(1, labels = c('elevated', 'ambient'), at = c(1,2))
title('Round 2')

dev.off()

#Experiment 2
#read in data from experiment 2
C2<-read.csv("CO2_log_exp2.csv")

# separate data by round
r1 <- C2 %>% filter(!is.na(Round.1)) # round 1 goes from Nov 17 2021 - Mar 18 2022
r2 <- C2 %>% filter(!is.na(Round.2)) # round 2 goes from Jan 20 2022 until end. 

# determine IQR to calculate outliers
# round 1
quantile(r1$aco2, na.rm = TRUE)
#       0%      25%      50%      75%     100% 
#   0.0000  427.6050  460.2765  518.1548  873.3140 
# IQR = 135.8247
# low outliers are below 324.4518
# high outliers are above 596.1012

quantile(r1$eco2, na.rm = TRUE)
# 0%  25%    50%   75%   100% 
# 5   620   645    678   1023  
# IQR = 87
# low outliers are below 558
# high outliers are above 732

# round 2
quantile(r2$aco2, na.rm = TRUE)
#       0%      25%      50%      75%     100% 
#   0.0000  417.6485  439.3290  491.1130 1034.7500  
# IQR = 110.197
# low outliers are below 329.132
# high outliers are above 549.526

quantile(r2$eco2, na.rm = TRUE)
# 0%    25%     50%   75%   100% 
# 0     607     640   675   1023  
# IQR = 102
# low outliers are below 538
# high outliers are above 742

# now that outliers are identified, add column in excel in data sheet to filter out outliers
# remove outliers
C2 <- C2 %>% filter(is.na(outliers))

# plot CO2 levels 
par(mfcol = c(1,2))
boxplot(as.numeric(r1$aco2), as.numeric(r1$CO2_LiCor), 
        outline = FALSE, col= c('#becbca', '#69a1ff'), 
        axes=FALSE, xlab="", ylab = "CO2 (ppm)", 
        ylim = c(300,800))
axis(1, labels = c('ambient', 'elevated'), at = c(1,2))
axis(2)
title('Round 1')

boxplot(as.numeric(r2$eco2), as.numeric(r2$aco2), 
        outline = FALSE, col= c('#69a1ff', '#becbca'), 
        axes=FALSE, ann = FALSE, 
        ylim = c(300,800))
axis(1, labels = c('elevated', 'ambient'), at = c(1,2))
title('Round 2')

dev.off()

