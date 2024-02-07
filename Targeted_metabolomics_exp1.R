## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## small, species-specific effects on pollen chemistry and plant growth across
  ## flowering plant species

# script to plot and analyze targeted metabolomics from experiment 1

library(ggplot2)
library(tidyverse)

#Load and clean data
data <- read.csv('Pollen_metabolomics_exp1.csv')
data[data == "N/F"] <- NaN
data.1 <- data %>% filter(Round == "1")
data.2 <- data %>% filter(Round == "2")

## Round 1 
# remove Squash data, sample sizes too small
data.1 <- data.1 %>% filter(species != "squash")
# also remove tomato for small sample sizes
data.1 <- data.1 %>% filter(species != "tomato")

#Get list of columns with reasonable data representation across species
cls <- c(6, 12, 13, 14, 15, 16, 19, 21, 22, 23)

#Convert columns to numeric
for(i in 1:length(cls)){
  
  data.1[,cls[i]] <- as.numeric(data.1[,cls[i]])
}

#Run principal components analysis
pc.data.1 <- data.1[,cls]
ind <- complete.cases(pc.data.1)
pc.data.1 <- pc.data.1[ind,]
pca.1 <- prcomp(pc.data.1, scale= TRUE)

#Write PCs back to original data frame
data.1$PC1[ind] <- pca.1$x[,1]
data.1$PC2[ind] <- pca.1$x[,2]

#Generate figure and write to current working directory
data$CO2.treatment <- as.factor(data$CO2.treatment)
data$species <- as.factor(data$species)
dev.off()
ggplot(data.1, aes(x = PC1, y = PC2, colour = CO2.treatment, shape = species))+
  geom_point()+
  stat_ellipse()+
  theme_classic()+
  scale_color_manual(values = c("#becbca","#69a1ff"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_shape_manual(values = c(15,16,17,18,0),
                     labels = c("Buckwheat", 
                                "Poppy", "Sunflower"),
                     name = "Plant species")

#Test whether CO2 and/or plant species were significant predictors of secondary metabolomics
model <- manova(cbind(PC1, PC2)~CO2.treatment*species, data = data.1)
summary(model)
# both CO2 and species were significant predictors of secondary chemistry but no interaction between the two

## Round 2 
#Get list of columns with reasonable data representation across species
cls <- c(6, 12, 13, 14, 15, 16, 19, 21, 22, 23)

#Convert columns to numeric
for(i in 1:length(cls)){
  
  data.2[,cls[i]] <- as.numeric(data.2[,cls[i]])
}

#Run principal components analysis
pc.data.2 <- data.2[,cls]
ind <- complete.cases(pc.data.2)
pc.data.2 <- pc.data.2[ind,]
pca.2 <- prcomp(pc.data.2, scale= TRUE)

#Write PCs back to original data frame
data.2$PC1[ind] <- pca.2$x[,1]
data.2$PC2[ind] <- pca.2$x[,2]

#Generate figure and write to current working directory
data.2$CO2.treatment <- as.factor(data.2$CO2.treatment)
data.2$species <- as.factor(data.2$species)
dev.off()
ggplot(data.2, aes(x = PC1, y = PC2, colour = CO2.treatment, shape = species))+
  geom_point()+
  stat_ellipse()+
  theme_classic()+
  scale_color_manual(values = c("#becbca","#69a1ff"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_shape_manual(values = c(15,1,17),
                     labels = c("Buckwheat", 
                                "Squash", "Sunflower"),
                     name = "Plant species")


#Test whether CO2 and/or plant species were significant predictors of secondary metabolomics
model <- manova(cbind(PC1, PC2)~CO2.treatment*species, data = data.2)
summary(model)
# both CO2 and species were significant predictors of secondary chemistry with an interaction between the two

# Compare each species individually
# Buckwheat
bw <- data.2 %>% filter(species == "buckwheat")
model <- manova(cbind(PC1, PC2)~CO2.treatment, data = bw)
summary(model)
# CO2 not significant for BW

# Sunflower
sf <- data.2 %>% filter(species == "sunflower")
model <- manova(cbind(PC1, PC2)~CO2.treatment, data = sf)
summary(model)
# CO2 not significant for SF

#Squash
s <- data.2 %>% filter(species == "squash")
model <- manova(cbind(PC1, PC2)~CO2.treatment, data = s)
summary(model)
# CO2 significant for S


