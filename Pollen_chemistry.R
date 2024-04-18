## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## small, species-specific effects on pollen chemistry and plant growth across
  ## flowering plant species

library("tidyverse")
library("ggplot2")
library("ggpubr")
library("lme4")
library("DHARMa")
library("lmerTest")
library("emmeans")
library("multcomp")

##############################
## read in the data frames
W.pollen<-read.csv(file = "Pollen_chemistry_exp2.csv", header = T, na.strings = c("", "NULLL"))
H.pollen <- read.csv(file = "Pollen_chemistry_exp1.csv", header = T, na.string = c("", "NULL"))

# clean data frames
# Pollen nutrition exp 1
H.pollen <- H.pollen %>% filter(is.na(Omit))
H.pollen$CO2 <- as.factor(H.pollen$CO2)
H.pollen$Round <- as.factor(H.pollen$Round)
H.pollen$Plant <- as.factor(H.pollen$Plant)
H.pollen$Chamber <- as.factor(H.pollen$Chamber)

# Pollen nutrition exp 2
W.pollen$Plant_SP <- factor(W.pollen$Plant_SP, levels = c("B", "BW", "C", "D", "LP", "N", "PP", "SA", "SF"))
W.pollen$Chamber <- factor(W.pollen$Chamber, levels = c("60", "63", "62", "61"))
# omit samples that were too low and duplicates
W.pollen <- W.pollen %>% filter(is.na(OMIT))
# add C:N ratio column
W.pollen$ratio <- W.pollen$C/W.pollen$N

############################
# summarize and plot data 
# remove partridge pea, too few samples
W.po.short <- W.pollen %>% filter(Plant_SP != "PP")

# remove plants, rounds, and chambers w/fewer than 3 samples
W.po.short <- W.po.short %>% group_by(Plant_SP, Chamber, Round) %>% filter(n()>2) %>% ungroup()
# ok good this didn't change anything, small sample sizes already omitted

# plot it
plants.2 <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia", "Nasturtium", "Sweet Alyssum", "Sunflower")
names(plants.2) <- c("B", "BW", "C", "LP", "N", "SA", "SF")

W.po.short$CO2 <- as.factor(W.po.short$CO2)

N<-ggplot(W.po.short, aes(x= Plant_SP, y = log(N), fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.3, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants.2)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))

C<- ggplot(W.po.short, aes(x= Plant_SP, y = C, fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.5, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants.2)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))

CN<-ggplot(W.po.short, aes(x= Plant_SP, y = log(ratio), fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.5, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants.2)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, size = .2),
        axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))

# add all to one plot
ggarrange(N, C, CN, nrow=3, ncol = 1)

### analysis
# N model
W.po.short1 <- W.po.short %>% filter(!is.na(N))
m.N <- lm(log(N) ~ CO2*Plant_SP + Round + Chamber, data = W.po.short)
plot(simulationOutput <- simulateResiduals(fittedModel = m.N, plot = F))
summary(m.N)
anova(m.N) # plant sp, round significant

emm_model1 <- emmeans(m.N, pairwise ~ CO2|Plant_SP)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1) 


## C:N ratio
m.CN <- lm(log(ratio) ~ CO2*Plant_SP + Round + Chamber, data = W.po.short)
plot(simulationOutput <- simulateResiduals(fittedModel = m.CN, plot = F))
summary(m.CN)
anova(m.CN) # plant species, co2 x plant species significant 

emm_model1 <- emmeans(m.CN, pairwise ~ CO2|Plant_SP)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

####
# %C
m.C <- lm(C ~ CO2*Plant_SP + Round + Chamber, data = W.po.short)
plot(simulationOutput <- simulateResiduals(fittedModel = m.C, plot = F))
summary(m.C)
anova(m.C) # significant effect of Round and Plant species

emm_model1 <- emmeans(m.C, pairwise ~ CO2|Plant_SP)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

###################
## Experiment 1 

# remove plants, rounds, and chambers w/fewer than 3 samples
H.short <- H.pollen %>% group_by(Plant, CO2, Round) %>% filter(n()>2) %>% ungroup()
H.short <- H.short %>% filter(is.na(Omit))
H.short$ratio <- H.short$C/H.short$N

plants.1 <- c("Buckwheat","Melon", "Partridge pea", "Poppy", "Squash", "Sunflower", "Tomatillo", "Tomato")
names(plants.1) <- c("buckwheat", "melon", "partridge pea", "poppy", "squash", "sunflower", "tomatillo", "tomato")

N<-ggplot(H.short, aes(x= Plant, y = log(N), fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.5, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants.1)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, size = .2),
        axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))

C<- ggplot(H.short, aes(x= Plant, y = C, fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.5, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants.1)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, size = .2),
        axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))

CN<-ggplot(H.short, aes(x= Plant, y = ratio, fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.5, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants.1)+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, size = .2),
        axis.title.x=element_blank())

ggarrange(N, C, CN, nrow=3, ncol =1, align = "h")

# Compare all plant species together 
# N
m.N <- lm(N ~ CO2*Plant + Round + Chamber, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N, plot =F))
summary(m.N)
anova(m.N)

emm_model1 <- emmeans(m.N, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

# C:N ratio
m.CN <- lm(ratio ~ CO2*Plant + Round + Chamber, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.CN, plot =F))
summary(m.CN)
anova(m.CN)

emm_model1 <- emmeans(m.CN, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

# %C
m.C <- lm(C ~ CO2*Plant + Round + Chamber, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C, plot =F))
summary(m.C)
anova(m.C)

emm_model1 <- emmeans(m.C, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)
