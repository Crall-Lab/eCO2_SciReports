## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## small, species-specific effects on pollen chemistry and plant growth across
  ## flowering plant species

# script to plot and analyze targeted metabolomics from experiment 1

library(ggplot2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(DHARMa)

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

#Separate out amino acids 
cls <- c(6, 12, 13, 15, 16, 19, 21, 22, 23)

#Convert columns to numeric
for(i in 1:length(cls)){
  
  data.1[,cls[i]] <- as.numeric(data.1[,cls[i]])
}

#Run principal components analysis
pc.data.1 <- data.1[,cls]
ind <- complete.cases(pc.data.1)
pc.data.1 <- pc.data.1[ind,]
pca.1 <- prcomp(pc.data.1, scale= TRUE)
summary(pca.1)

#Write PCs back to original data frame
data.1$PC1[ind] <- pca.1$x[,1]
data.1$PC2[ind] <- pca.1$x[,2]

#Generate figure and write to current working directory
data$CO2.treatment <- as.factor(data$CO2.treatment)
data$species <- as.factor(data$species)
dev.off()
pdf(file="round1.pdf", 
    width = 4.75, 
    height = 4)
ggplot(data.1, aes(x = PC1, y = PC2, colour = CO2.treatment, shape = species))+
  geom_point()+
  stat_ellipse()+
  theme_classic()+
  scale_color_manual(values = c("#becbca","#69a1ff"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_shape_manual(values = c(15,16,17, 18, 0),
                     labels = c("Buckwheat", 
                                "Poppy", "Sunflower"),
                     name = "Plant species")+
  labs(y = "PC2 (20.5%)", x = "PC1 (53.8%)")

#Test whether CO2 and/or plant species were significant predictors of secondary metabolomics
model <- manova(cbind(PC1, PC2)~CO2.treatment*species, data = data.1)
summary(model)
# both CO2 and species were significant predictors of secondary chemistry but no interaction between the two

## Round 2 
#Get list of columns with reasonable data representation across species
cls <- c(6, 12, 13, 15, 16, 19, 21, 22, 23)

#Convert columns to numeric
for(i in 1:length(cls)){
  
  data.2[,cls[i]] <- as.numeric(data.2[,cls[i]])
}

#Run principal components analysis
pc.data.2 <- data.2[,cls]
ind <- complete.cases(pc.data.2)
pc.data.2 <- pc.data.2[ind,]
pca.2 <- prcomp(pc.data.2, scale= TRUE)
summary(pca.2)

#Write PCs back to original data frame
data.2$PC1[ind] <- pca.2$x[,1]
data.2$PC2[ind] <- pca.2$x[,2]

#Generate figure and write to current working directory
data.2$CO2.treatment <- as.factor(data.2$CO2.treatment)
data.2$species <- as.factor(data.2$species)
dev.off()
pdf(file="round2.pdf", 
    width = 4.75, 
    height = 4)
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
                     name = "Plant species")+
  labs(y = "PC2 (29.4%)", x = "PC1 (58.5%)")

#Test whether CO2 and/or plant species were significant predictors of secondary metabolomics
model <- manova(cbind(PC1, PC2)~CO2.treatment*species, data = data.2)
summary(model)
# only species was a significant predictors of amino acids in pollen

# Compare each remaining compound individually
data$Round <- as.factor(data$Round)
# caffeine
# filter out caffeine samples, and then only keep samples with 3+ per plant, co2 treatment, and round
caf <- data %>% filter(Caffeine != "NaN")
caf <- caf %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
caf$Caffeine <- as.numeric(caf$Caffeine)
# analyze data and test fit 
m.caf <- lm(log(Caffeine)~CO2.treatment*species + Round, data = caf)
plot(simulationOutput <- simulateResiduals(fittedModel = m.caf, plot = F))
anova(m.caf) # round and species significant
  
# chlorogenic acid
# filter out chlorogenic acid samples, and then only keep samples with 3+ per plant, co2 treatment, and round
cha <- data %>% filter(Chlorogenic.acid != "NaN")
cha <- cha %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
cha$Chlorogenic.acid <- as.numeric(cha$Chlorogenic.acid)
# analyze data and test fit 
m.cha <- lm(log(Chlorogenic.acid)~CO2.treatment*species + Round, data = cha)
plot(simulationOutput <- simulateResiduals(fittedModel = m.cha, plot = F))
anova(m.cha) # round and species, co2 x species significant

# investigate CO2 x species independently
# buckwheat
cha.bw <- cha %>% filter(species == "buckwheat")
m.cha.bw <- lm(log(Chlorogenic.acid) ~ CO2.treatment + Round, data = cha.bw)
plot(simulationOutput <- simulateResiduals(fittedModel = m.cha.bw, plot = F))
anova(m.cha.bw) # co2 significant

# sunflower
cha.sf <- cha %>% filter(species == "sunflower")
m.cha.sf <- lm(log(Chlorogenic.acid) ~ CO2.treatment + Round, data = cha.sf)
plot(simulationOutput <- simulateResiduals(fittedModel = m.cha.sf, plot = F))
anova(m.cha.sf) # round significant

# tomato
cha.t <- cha %>% filter(species == "tomato") # only from round 1
t.test(Chlorogenic.acid ~ CO2.treatment, data = cha.t) # not significant. 

# cinnamic acid
# filter out cinnamic acid samples, and then only keep samples with 3+ per plant, co2 treatment, and round
cma <- data %>% filter(Cinnamic.acid != "NaN")
cma <- cma %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
cma$Cinnamic.acid <- as.numeric(cma$Cinnamic.acid)
# analyze data and test fit 
m.cma <- lm(log(Cinnamic.acid)~CO2.treatment*species + Round, data = cma)
plot(simulationOutput <- simulateResiduals(fittedModel = m.cma, plot = F))
anova(m.cma) # round and species, co2 x species significant

# eugenol
# filter out eugenol samples, and then only keep samples with 3+ per plant, co2 treatment, and round
eug <- data %>% filter(Eugenol != "NaN")
eug <- eug %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
eug$Eugenol <- as.numeric(eug$Eugenol)
# analyze data and test fit 
m.eug <- lm(Eugenol~CO2.treatment*species, data = eug) # not enough data from both rounds to include in model
plot(simulationOutput <- simulateResiduals(fittedModel = m.eug, plot = F))
anova(m.eug) # nothing significant

# gallic acid
# filter out gallic acid samples, and then only keep samples with 3+ per plant, co2 treatment, and round
gal <- data %>% filter(Gallic.acid != "NaN")
gal <- gal %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
gal$Gallic.acid <- as.numeric(gal$Gallic.acid)
# analyze data and test fit 
m.gal <- lm(Gallic.acid~CO2.treatment*species + Round, data = gal) 
plot(simulationOutput <- simulateResiduals(fittedModel = m.gal, plot = F))
anova(m.gal) # round and species significant

# kaempferol
# filter out kaempferol samples, and then only keep samples with 3+ per plant, co2 treatment, and round
kae <- data %>% filter(Kaempferol != "NaN")
kae <- kae %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
kae$Kaempferol <- as.numeric(kae$Kaempferol)
# analyze data and test fit 
m.kae <- lm(log(Kaempferol)~CO2.treatment*species + Round, data = kae) 
plot(simulationOutput <- simulateResiduals(fittedModel = m.kae, plot = F))
anova(m.kae) # round and species, co2 x species significant

# nicotine
# filter out nicotine samples, and then only keep samples with 3+ per plant, co2 treatment, and round
nic <- data %>% filter(Nicotine != "NaN")
nic <- nic %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
nic$Nicotine <- as.numeric(nic$Nicotine)
# buckwheat in round 1 only samples w/nicotine
t.test(Nicotine~CO2.treatment, data = nic) # no significant difference. 

# p-coumaric acid
# filter out p-coumaric acid samples, and then only keep samples with 3+ per plant, co2 treatment, and round
pc.a <- data %>% filter(P.coumaric.acid != "NaN")
pc.a <- pc.a %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
pc.a$P.coumaric.acid <- as.numeric(pc.a$P.coumaric.acid)
# analyze data and test fit 
m.pc.a <- lm(log(P.coumaric.acid)~CO2.treatment*species + Round, data = pc.a) 
plot(simulationOutput <- simulateResiduals(fittedModel = m.pc.a, plot = F))
anova(m.pc.a) # round and species

# quercetin
# filter out quercetin samples, and then only keep samples with 3+ per plant, co2 treatment, and round
que <- data %>% filter(Quercitin != "NaN")
que <- que %>% group_by(species, Round, CO2.treatment) %>% filter(n()>2) %>% ungroup()
que$Quercitin <- as.numeric(que$Quercitin)
# analyze data and test fit 
m.que <- lm(log(Quercitin)~CO2.treatment*species + Round, data = que) 
plot(simulationOutput <- simulateResiduals(fittedModel = m.que, plot = F))
anova(m.que) # round and species

# compare each species separately
# buckwheat
que.bw <- que %>% filter(species == "buckwheat")
m.que.bw <- lm(Quercitin ~ CO2.treatment + Round, data = que.bw)
plot(simulationOutput <- simulateResiduals(fittedModel = m.que.bw, plot = F))
anova(m.que.bw) # co2 significant

# poppy
que.p <- que %>% filter(species == "poppy") # only round 1 data
t.test(Quercitin ~ CO2.treatment, data = que.p) # significant effect of CO2 - eCO2 increased Quercitin

# squash
que.sq <- que %>% filter(species == "squash") # only round 1 data
t.test(Quercitin ~ CO2.treatment, data = que.sq) # not significant
 
# sunflower
que.sf <- que %>% filter(species == "sunflower") # only round 1 data
t.test(Quercitin ~ CO2.treatment, data = que.sf) # not significant


