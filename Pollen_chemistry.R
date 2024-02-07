## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## small, species-specific effects on pollen chemistry and plant growth across
  ## flowering plant species

library("tidyverse")
library("ggplot2")
library("ggpubr")
library("lme4")
library("DHARMa")
library("lmerTest")

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
m.N <- lmer(N ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.N, plot = F))
summary(m.N)
anova(m.N) # plant sp, round significant

## C:N ratio
m.CN <- lmer(ratio ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.CN, plot = F))
summary(m.CN)
anova(m.CN) # plant species, co2 x plant species significant 

####
# %C
m.C <- lmer(C ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.C, plot = F))
summary(m.C)
anova(m.C) # significant effect of Round and Plant species

# compare each plant species independently

## Borage
# %N (proxy for pollen protein)
B.po <- W.pollen %>% filter(Plant_SP == "B")

m.B1 <- lmer(N ~ CO2 + Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))
m.B2 <- lmer(N ~ CO2 * Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B2, plot =F))
summary(m.B1)
summary(m.B2)
anova(m.B1, m.B2) # no sig difference, lower AIC for m.B1, so use that one 
anova(m.B1)
anova(m.B2)

# C:N ratio
m.B3 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B3, plot =F))
m.B4 <- lmer(ratio ~ CO2 * Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B4, plot =F))
summary(m.B3)
summary(m.B4)
anova(m.B4, m.B3) # no difference, but m.B2 has lower AIC 
anova(m.B4) # Again, only Round is significant

# C 
m.B5 <- lmer(C ~ CO2 + Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B5, plot = F))
summary(m.B5)
anova(m.B5) # only round is significant

## Buckwheat
# %N (proxy for pollen protein)
BW.po <- W.pollen %>% filter(Plant_SP == "BW")

m.BW1 <- lmer(N ~ CO2 + Round + (1|Chamber), data = BW.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1, plot =F))
m.BW2 <- lmer(N ~ CO2 * Round + (1|Chamber), data = BW.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW2, plot =F))
summary(m.BW1)
anova(m.BW1)# significant effect of CO2, CO2 increases %N
anova(m.BW2)
anova(m.BW1, m.BW2) # virtually no difference between models

# C 
m.BW3 <- lmer(C ~ CO2 + Round + (1|Chamber), data = BW.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.BW3, plot = F))
summary(m.BW3)
anova(m.BW3) # nothing significant

# CN ratio
m.BW4 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = BW.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.BW4, plot = F))
summary(m.BW4)
anova(m.BW4) # nothing significant

## Red Clover
# %N (proxy for pollen protein)
C.po <- W.pollen %>% filter(Plant_SP == "C")

m.C1 <- lmer(N ~ CO2+Round + (1|Chamber), data = C.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1, plot =F))
m.C2 <- lmer(N ~ CO2*Round + (1|Chamber), data = C.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C2, plot =F))
summary(m.C1)
anova(m.C1) # No significant difference by round or CO2 level
anova(m.C1, m.C2) # no difference, slightly lower AIC in m.C1

# C:N ratio
m.C3 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = C.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.C3, plot = F))
summary(m.C3)
anova(m.C3) # No significant difference by round or CO2 level

# C 
m.C4 <- lmer(C ~ CO2 + Round + (1|Chamber), data = C.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.C4, plot = F))
summary(m.C4)
anova(m.C4)

## Lacy Phacelia
# %N (proxy for pollen protein)
LP.po <- W.pollen %>% filter(Plant_SP == "LP")

m.LP1 <- lmer(N ~ CO2+Round + (1|Chamber), data = LP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
m.LP2 <- lmer(N ~ CO2*Round + (1|Chamber), data = LP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP2, plot =F))
summary(m.LP1)
anova(m.LP1) # significant difference by round 
anova(m.LP1, m.LP2) # no difference, slightly lower AIC in m.LP1

# C:N ratio
m.LP3 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = LP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP3, plot =F))
summary(m.LP3)
anova(m.LP3) # No significant difference by CO2 level, but effect of round

# C
m.LP4 <- lmer(C ~ CO2 + Round + (1|Chamber), data = LP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP4, plot =F))
summary(m.LP4)
anova(m.LP4) # effect of round

## Nasturtium
# %N (proxy for pollen protein)
N.po <- W.pollen %>% filter(Plant_SP == "N")

m.N1 <- lmer(N ~ CO2 + Round + (1|Chamber), data = N.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
m.N2 <- lmer(N ~ CO2 * Round + (1|Chamber), data = N.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N2, plot =F))
summary(m.N1)
anova(m.N1) # No significant difference by round or CO2 level, but effect of round
anova(m.N1, m.N2) # no significant difference, slightly lower AIC in m.N1

# C:N ratio
m.N3 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = N.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N3, plot =F))
summary(m.N3)
anova(m.N3) # No significant difference by round or CO2 level, but effect of round

# C 
m.N4 <- lmer(C ~ CO2 + Round + (1|Chamber), data = N.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N4, plot =F))
summary(m.N4)
anova(m.N4) # nothing is significant

## Partridge pea
# %N (proxy for pollen protein)
PP.po <- W.pollen %>% filter(Plant_SP == "PP")

m.PP1 <- lmer(N ~ CO2 + (1|Chamber), data = PP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1, plot =F))
summary(m.PP1)

m.PP2 <- lmer(C ~ CO2 + (1|Chamber), data = PP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP2, plot =F))
summary(m.PP2)

m.PP3 <- lmer(ratio ~ CO2 + (1|Chamber), data = PP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP3, plot =F))
summary(m.PP3)

## Sweet Alyssum 
# %N (proxy for pollen protein)
SA.po <- W.pollen %>% filter(Plant_SP == "SA")

m.SA1 <- lmer(N ~ CO2 + Round + (1|Chamber), data = SA.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
m.SA2 <- lmer(N ~ CO2 * Round + (1|Chamber), data = SA.po, REML = F)
anova(m.SA1, m.SA2) # No significant difference by round or CO2 level, but effect of round
anova(m.SA1)

# C:N ratio
m.SA3 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = SA.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA3, plot =F))
summary(m.SA3)
anova(m.SA3) # No significant difference by round or CO2 level, but effect of round

# C 
m.SA4 <- lmer(C ~ CO2 + Round + (1|Chamber), data = SA.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA4, plot =F))
summary(m.SA4)
anova(m.SA4) # nothing is significant

## Sunflower 
# %N (proxy for pollen protein)
SF.po <- W.pollen %>% filter(Plant_SP == "SF")

m.SF1 <- lmer(N ~ CO2 * Round + (1|Chamber), data = SF.po, REML=F) # better fit
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
m.SF2 <- lmer(N ~ CO2 + Round + (1|Chamber), data = SF.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SFs, plot =F))
anova(m.SF1, m.SF2)
summary(m.SF1)
anova(m.SF1) # CO2, Round, and the interaction are significant
anova(m.SF2)

# C:N ratio
m.SF3 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = SF.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF3, plot =F))
summary(m.SF3)
anova(m.SF3) # CO2, Round, and the interaction are significant

# % C
m.SF4 <- lmer(C ~ CO2 + Round + (1|Chamber), data = SF.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF4, plot =F))
summary(m.SF4)
anova(m.SF4) # nothing is significant

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
  scale_x_discrete(labels = plants)+
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
  scale_x_discrete(labels = plants)+
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
m.N <- lm(N ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N, plot =F))
summary(m.N)
anova(m.N)

# C:N ratio
m.CN <- lm(ratio ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.CN, plot =F))
summary(m.CN)
anova(m.CN)

# %C
m.C <- lm(C ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C, plot =F))
summary(m.C)
anova(m.C)

# compare each plant species separately 
## Melon - only R1 
M.po <- H.short %>% filter(Plant == "melon")
# N
t.test(N~CO2, M.po)
# C:N ratio
t.test(ratio~CO2, M.po)
# C
t.test(C~CO2, M.po)

## Poppy - only R1 
P.po <- H.short %>% filter(Plant == "poppy")
P.po <- P.po %>% filter(Round =="1")
# N
t.test(N~CO2, P.po)
# C:N ratio
t.test(ratio~CO2, P.po)
# C
t.test(C~CO2, P.po)

## Squash - both rounds 
S.po <- H.short %>% filter(Plant == "squash")
m.S.N <- lm(N ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.S.N, plot =F))
summary(m.S.N)
anova(m.S.N)
# C
m.S.C <- lm(C ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.S.C, plot =F))
summary(m.S.C)
anova(m.S.C)
# C:N
m.S.CN <- lm(ratio ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.S.CN, plot =F))
summary(m.S.CN)
anova(m.S.CN)

## Sunflower 
SF.po <- H.short %>% filter(Plant == "sunflower")
SF.po <- SF.po %>% filter(Round == "1")

# %N
t.test(N~CO2, SF.po)
# C:N ratio
t.test(ratio~CO2, SF.po)
# C
t.test(C~CO2, SF.po)

# Tomatillo 
TT.po <- H.short %>% filter(Plant == "tomatillo")
# %N
t.test(N~CO2, TT.po)
# C:N ratio
t.test(ratio~CO2, TT.po)
# C
t.test(C~CO2, TT.po)

# Tomato 
T.po <- H.short %>% filter(Plant == "tomato")
# %N
t.test(N~CO2, T.po)
# C:N ratio
t.test(ratio~CO2, T.po)
# C
t.test(C~CO2, T.po)

