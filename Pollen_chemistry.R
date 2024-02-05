## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
## inconsistent effects on pollen chemistry and plant growth across
## flowering plant species

library("tidyverse")
library("ggplot2")
library("lme4")
library("lmerTest")
library("multilevelTools")
library("JWileymisc")
library("nlme")
library("scales")
library("beeswarm")
library("see")
library("DHARMa")
library("GLMMadaptive")
library("ggpubr")

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
W.po.short <- W.pollen %>% filter(Plant_SP != "PP")

# remove plants, rounds, and chambers w/fewer than 3 samples
W.po.short <- W.po.short %>% group_by(Plant_SP, Chamber, Round) %>% filter(n()>2) %>% ungroup()
# ok good this didn't change anything, small sample sizes already omitted

# plot w/just CO2 treatments
po.sum <- W.po.short %>% group_by(Plant_SP, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(log(N), na.rm = T),
    sd = sd(log(N), na.rm = T)
  )
po.sum$se <- po.sum$sd/sqrt(po.sum$count)

po.sum$CO2 <-factor(po.sum$CO2, levels = c(0, 1))
po.sum$chemistry <- "N"

po.sum1 <- W.po.short %>% group_by(Plant_SP, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(log(C), na.rm = T),
    sd = sd(log(C), na.rm = T)
  )
po.sum1$se <- po.sum1$sd/sqrt(po.sum1$count)

po.sum1$CO2 <-factor(po.sum1$CO2, levels = c(0, 1))
po.sum1$chemistry <- "C"

# log transfer to bring buckwheat back with everyone else
W.po.short$log_ratio <- log(W.po.short$ratio)

po.sum2 <- W.po.short %>% group_by(Plant_SP, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(ratio, na.rm = T),
    sd = sd(ratio, na.rm = T)
  )
po.sum2$se <- po.sum2$sd/sqrt(po.sum2$count)

po.sum2$CO2 <-factor(po.sum2$CO2, levels = c(0, 1))
po.sum2$chemistry <- "log(C:N)"


# put them all together
WI.pc <- rbind(po.sum, po.sum1, po.sum2)
WI.pc$chemistry <- factor(WI.pc$chemistry, levels = c("N", "C", "log(C:N)"))

# plot it
plants <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia", "Nasturtium", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "LP", "N", "SA", "SF")

ggplot(WI.pc, aes(x=Plant_SP, y = mean, color = CO2))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))+
  labs(y="Pollen chemistry ± se")+
  scale_color_manual(values = c("grey","navyblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_x_discrete(labels = plants)+
  facet_wrap(vars(chemistry), 
             ncol = 1, 
             scales = "free")

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
  scale_x_discrete(labels = plants)+
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
  scale_x_discrete(labels = plants)+
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

#ggplot(W.po.short, aes(x= Plant_SP, y = C, fill = CO2))+
#  geom_boxplot(aes(fill = CO2))+
#  theme_classic()+
#  scale_fill_manual(values = c("grey", "cornflowerblue"))+
#  scale_x_discrete(labels = plants)+
#  theme(legend.position = "none",
#        axis.text.x = element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.title.x=element_blank())

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


#ggplot(W.po.short, aes(x= Plant_SP, y = log(ratio), fill = CO2))+
#geom_boxplot(aes(fill = CO2))+
#theme_classic()+
#scale_fill_manual(values = c("grey", "cornflowerblue"))+
#theme(legend.position = "none",
#      axis.text.x = element_text(angle = 90, hjust=1, size = .2),
#      axis.title.x=element_blank())

ggarrange(N, C, CN, nrow=3, ncol = 1)


# same but log transform the data? Maybe, come back to this later if needed

po.sum <- W.po.short %>% group_by(Plant_SP, Chamber, Round) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(N, na.rm = T),
    sd = sd(N, na.rm = T)
  )
po.sum$se <- po.sum$sd/sqrt(po.sum$count)

# plot it
plants <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia", "Nasturtium", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "LP", "N", "SA", "SF")

po.sum$Round <-as.factor(po.sum$Round)

ggplot(po.sum, aes(x=Round, y = mean, shape = Chamber))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = c(.9,.2))+
  labs(y="%N +/- se", x = "Experimental Round")+
  scale_shape_manual(values = c(0,15,1,16), labels = c("60 - eCO2", "63 - eCO2", "62 - aCO2", "61 - aCO2"))+
  facet_wrap(vars(Plant_SP), scales = "free",
             labeller = labeller(Plant_SP = plants),
             ncol = 4)

po.sum1 <- W.po.short %>% group_by(Plant_SP, Chamber, Round) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(ratio, na.rm = T),
    sd = sd(ratio, na.rm = T)
  )
po.sum1$se <- po.sum1$sd/sqrt(po.sum1$count)

po.sum1$Round <-as.factor(po.sum1$Round)

ggplot(po.sum1, aes(x=Round, y = mean, shape = Chamber))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = c(.9,.2))+
  labs(y="C:N ratio +/- se", x = "Experimental Round")+
  scale_shape_manual(values = c(0,15,1,16), labels = c("60 - eCO2", "63 - eCO2", "62 - aCO2", "61 - aCO2"))+
  facet_wrap(vars(Plant_SP), scales = "free",
             labeller = labeller(Plant_SP = plants),
             ncol = 4)

# need to omit plant+chamber+round situations where there were less than three observations (pollen collections)
# all PP

# ok let's think about how to analyze this statistically
# We need to include chamber as a random effect
# Round as a fixed effect
# Plant species 

# N content as outcome variable
# CO2 level as predictor variable
# Plant species indicate which N observations belong to which plant

# how many unique plants are there?
length(unique(W.po.short$Plant)) # 8 - not pollen for dandelions, this is right

# briefly explore each variable (really only would want to explore N)
summary(W.po.short$N)

# distribution of the variables visually
tmp <- meanDecompose(N ~ Plant_SP, data = W.po.short)
str(tmp)

plot(testDistribution(tmp[["N by Plant_SP"]]$X,
                      extremevalues = "theoretical", ev.perc = .001),
     varlab = "Between Plant N%")

plot(testDistribution(tmp[["N by residual"]]$X,
                      extremevalues = "theoretical", ev.perc = .001),
     varlab = "Within plant N%")

###
# N model
# full model
W.po.short1 <- W.po.short %>% filter(!is.na(N))
m <- lmer(N ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m, plot = F))
plotResiduals(simulationOutput, form = W.po.short$Plant_SP)
testDispersion(simulationOutput, alternative = "less")
summary(simulationOutput)
summary(m)
anova(m)
AIC(m) # 762.359
# only plant species is significant and Round x Plant species, best fitting model

# test w/o CO2 and see how it performs 
m1 <- lmer(N ~ Round*Plant_SP + (1|Chamber), data = W.po.short, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m1, plot = F))
summary(m1)
anova(m1)
AIC(m1) # 746.531


fitme(N~CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, family = COMPoission())
m.binom <- glmer(N ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, family = COMPoission(nu=))
summary(m.binom) # AIC = 1187.6summary(m.binom) # AIC = 1187.6summary(m.binom) # AIC = 1187.6
plot(simulationOutput <- simulateResiduals(fittedModel = m.binom, plot = F))


anova(m, m1) # no significant difference

# test w/formula James suggested
m.1 <-lmer(N ~ CO2*Plant_SP + Round + (1|Chamber),  data = W.po.short, REML=F)
summary(m.1)
anova(m.1) # Plant SP and Round significant 
AIC(m.1) # 827.3099

# compare between the two:
anova(m, m.1) # m is a significantly better fit, so we will stick with this for now 

## C:N ratio
m <- lmer(ratio ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
testDispersion(m)
simulationOutput <- simulateResiduals(fittedModel = m, plot = F)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testDispersion(simulationOutput, alternative = "less")
summary(m)
anova(m) # only plant species significant 
AIC(m) # 1488.357

# remove CO2
m1 <- lmer(ratio ~ Round*Plant_SP + (1|Chamber), data = W.po.short, REML=F)
summary(m1)
anova(m1) # only plant species significant 
AIC(m1) # 1491.561

anova(m, m1) # full model significantly better, ok so keep CO2 in

# remove Round
m2 <- lmer(ratio ~ CO2*Plant_SP + (1|Chamber), data = W.po.short, REML=F)
summary(m2)
anova(m2) #  plant species and CO2xPlant species significant - so something is happening w/CO2 at the plant species level 
AIC(m2) # 1478.852

anova(m, m2) # no difference between this and full model
anova(m, m1, m2) # no difference between the three

# forumla james suggested
m3 <- lmer(ratio ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
summary(m3)
anova(m3) #  plant species and CO2xPlant species significant - so something is happening w/CO2 at the plant species level 
AIC(m3) # 1480.653

anova(m,m3) # no difference from full model

####
# %C
qqPlot(W.po.short$C)
m <- lmer(C ~ CO2*Round*Plant_SP + (1|Chamber), data = W.po.short, REML=F)
testDispersion(m)
simulationOutput <- simulateResiduals(fittedModel = m, plot = F)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
testDispersion(simulationOutput, alternative = "greater")
plot(simulationOutput)
summary(m)
anova(m) # significant effect of Round and Plant species
AIC(m) # 1190.176

# remove CO2
m1 <- lmer(C ~ Round*Plant_SP + (1|Chamber), data = W.po.short, REML=F)
summary(m1)
anova(m1) # significant effect of Round and Plant species and Round x Plant sp
AIC(m1) # 1167.712

anova(m, m1) # no difference between models

# formula james suggested
m2 <- lmer(C ~ CO2*Plant_SP + Round + (1|Chamber), data = W.po.short, REML=F)
summary(m2)
anova(m2) # significant effect of Round and Plant species and Round x Plant sp
AIC(m2) # 1190.213

anova(m, m2) # full model significantly better predictor than this. 

# ok what if we compare each plant separately?

########
## Borage
########
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

# ok new model shows effect of CO2 and of Round, compare %N
B.po.1 <- B.po %>% filter(Round == "1")
t.test(N~CO2, data = B.po.1)
##t(21.597) = -1.4346, p = 0.1657 # round 1 not significant 
B.po.2 <- B.po %>% filter(Round == "2")
t.test(N~CO2, data = B.po.2)
##t(20.653) = -2.1287, p = 0.04549 # round 2 significant 

sum <- B.po %>% group_by(CO2, Round) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(N, na.rm = T),
    sd = sd(N, na.rm = T)
  )
sum$se <- sum$sd/sqrt(sum$count)

# C:N ratio
m.B2 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B2, plot =F))
m.B3 <- lmer(ratio ~ CO2 * Round + (1|Chamber), data = B.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B3, plot =F))
summary(m.B2)
summary(m.B3)
anova(m.B2, m.B3) # no difference, but m.B2 has lower AIC 
anova(m.B2) # Again, only Round is significant
AIC(m.B2) # -5.510952

# ok new model shows effect of CO2 and of Round, compare %N
B.po.1 <- B.po %>% filter(Round == "1")
t.test(ratio~CO2, data = B.po.1)
##t(21.597) = -1.4346, p = 0.1657 # round 1 not significant 
B.po.2 <- B.po %>% filter(Round == "2")
t.test(ratio~CO2, data = B.po.2)
##t(20.653) = -2.1287, p = 0.04549 # round 2 significant 

sum <- B.po %>% group_by(CO2, Round) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(ratio, na.rm = T),
    sd = sd(ratio, na.rm = T)
  )
sum$se <- sum$sd/sqrt(sum$count)


# C 
m.B1 <- lmer(C ~ CO2 + Round + (1|Chamber), data = B.po, REML=F)
testDispersion(m.B1)
simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
summary(m.B1)
anova(m.B1) # only round is significant
AIC(m.B1) # 132.6573


#######
## Buckwheat
#########
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
AIC(m.BW1) # 45.4868

t.test(N~CO2, data = BW.po)
# significant difference of N between CO2 treatments
# t(35.032) = -2.1971, p = 0.03472

BW.p.sum <- BW.po %>% group_by(CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(N, na.rm = T),
    sd = sd(N, na.rm = T)
  )
BW.p.sum$se <- BW.p.sum$sd/sqrt(BW.p.sum$count)

# C 
m.B1 <- lmer(C ~ CO2 + Round + (1|Chamber), data = BW.po, REML=F)
testDispersion(m.B1)
simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
summary(m.B1)
anova(m.B1) # nothing significant
AIC(m.B1) # 219.5807

# CN ratio
m.B1 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = BW.po, REML=F)
testDispersion(m.B1)
simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
summary(m.B1)
anova(m.B1) # nothing significant
AIC(m.B1) # 322.9815

######
## Red Clover
######
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
m.C2 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = C.po, REML=F)
testDispersion(m.C2)
simulationOutput <- simulateResiduals(fittedModel = m.C2, plot = F)
plot(simulationOutput)
summary(m.C2)
anova(m.C2) # No significant difference by round or CO2 level
AIC(m.C2) # 69.4536

# C 
m.C2 <- lmer(C ~ CO2 + Round + (1|Chamber), data = C.po, REML=F)
testDispersion(m.C2)
simulationOutput <- simulateResiduals(fittedModel = m.C2, plot = F)
plot(simulationOutput)
summary(m.C2)
summary(m.C2)
anova(m.C2)
AIC(m.C2) # 81.84848

#######
## Lacy Phacelia
########
# %N (proxy for pollen protein)
LP.po <- W.pollen %>% filter(Plant_SP == "LP")

m.LP1 <- lmer(N ~ CO2+Round + (1|Chamber), data = LP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
m.LP2 <- lmer(N ~ CO2*Round + (1|Chamber), data = LP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP2, plot =F))
summary(m.LP1)
anova(m.LP1) # significant difference by round 
anova(m.LP1, m.LP2) # no difference, slightly lower AIC in m.LP1

# try without CO2
m.LP2 <- lmer(N ~ Round + (1|Chamber), data = LP.po, REML=F)
summary(m.LP2)
anova(m.LP2) # round still significant
AIC(m.LP2) # 230.8186

anova(m.LP1, m.LP2) # no difference

# C:N ratio
m.LP2 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = LP.po, REML=F)
testDispersion(m.LP2)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP2, plot =F))
summary(m.LP2)
anova(m.LP2) # No significant difference by CO2 level, but effect of round
AIC(m.LP2) # 257.6186

# remove CO2
m.LP3 <- lmer(ratio ~ Round + (1|Chamber), data = LP.po, REML=F)
summary(m.LP3)
anova(m.LP3) # sig effect of round
AIC(m.LP3) # 254.9079

anova(m.LP2, m.LP3) # no difference

# C
m.LP2 <- lmer(C ~ CO2 + Round + (1|Chamber), data = LP.po, REML=F)
testDispersion(m.LP2)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP2, plot =F))
summary(m.LP2)
anova(m.LP2) # effect of round
AIC(m.LP2) # 243.7541

# remove CO2
m.LP3 <- lmer(C ~ Round + (1|Chamber), data = LP.po, REML=F)
summary(m.LP3)
anova(m.LP3) # effect of round
AIC(m.LP3) # 240.0096

anova(m.LP2, m.LP3) # no difference

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
m.N2 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = N.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N2, plot =F))
summary(m.N2)
anova(m.N2) # No significant difference by round or CO2 level, but effect of round

# C 
m.N2 <- lmer(C ~ CO2 + Round + (1|Chamber), data = N.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N2, plot =F))
summary(m.N2)
anova(m.N2) # nothing is significant

## Partridge pea
# %N (proxy for pollen protein)
PP.po <- W.pollen %>% filter(Plant_SP == "PP")

m.N1 <- lmer(N ~ CO2 + (1|Chamber), data = PP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
summary(m.N1)

m.N1 <- lmer(C ~ CO2 + (1|Chamber), data = PP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
summary(m.N1)

m.N1 <- lmer(ratio ~ CO2 + (1|Chamber), data = PP.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
summary(m.N1)

## Sweet Alyssum 
# %N (proxy for pollen protein)
SA.po <- W.pollen %>% filter(Plant_SP == "SA")

m.SA1 <- lmer(N ~ CO2 + Round + (1|Chamber), data = SA.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
m.SA2 <- lmer(N ~ CO2 * Round + (1|Chamber), data = SA.po, REML = F)
anova(m.SA1, m.SA2) # No significant difference by round or CO2 level, but effect of round

# C:N ratio
m.SA2 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = SA.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA2, plot =F))
summary(m.SA2)
anova(m.SA2) # No significant difference by round or CO2 level, but effect of round

# C 
m.SA2 <- lmer(C ~ CO2 + Round + (1|Chamber), data = SA.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA2, plot =F))
summary(m.SA2)
anova(m.SA2) # nothing is significant

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
m.SF2 <- lmer(ratio ~ CO2 + Round + (1|Chamber), data = SF.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF2, plot =F))
summary(m.SF2)
anova(m.SF2) # CO2, Round, and the interaction are significant

# % C
m.SF1 <- lmer(C ~ CO2 + Round + (1|Chamber), data = SF.po, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
summary(m.SF1)
anova(m.SF1) # nothing is significant

###################
## Experiment 1 

# remove plants, rounds, and chambers w/fewer than 3 samples
H.short <- H.pollen %>% group_by(Plant, CO2, Round) %>% filter(n()>2) %>% ungroup()
H.short <- H.short %>% filter(is.na(Omit))

plants.1 <- c("Buckwheat","Melon", "Partridge pea", "Poppy", "Squash", "Sunflower", "Tomatillo", "Tomato")
names(plants.1) <- c("buckwheat", "melon", "partridge pea", "poppy", "squash", "sunflower", "tomatillo", "tomato")

H.short$ratio <- H.short$C/H.short$N

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
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)

# C:N ratio
m.B <- lm(ratio ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)

# %C
m.B <- lm(C ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
AIC(m.B) # 254.8034
# sig effect of C, Plant species, so compare each plant species separately

# each plant species separately 
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
m.B <- lm(N ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
# C
m.B <- lm(C ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
# C:N
m.B <- lm(ratio ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)

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

