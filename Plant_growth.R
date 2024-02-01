## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
## inconsistent effects on pollen chemistry and plant growth across
## flowering plant species

# script to plot and analyze plant growth data from exp 1

library("tidyverse")
library("ggplot2")
library("lubridate")

library("chron")
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

##############################
## read in all of the data frames you need
W.plant<-read.csv(file="Plant_growth_exp2.csv",header=T, na.strings = c("","NULL"))
W.area<-read.csv(file="Flower_size_exp2.csv", header=T, na.strings=c("","NULL"))
W.FF <- read.csv(file = "Flowering_initiation_exp2.csv", header = T, na.strings = c("", "NULL"))
biomass <- read.csv(file = "Biomass_exp2.csv", header = T, na.string = c("", "NULL"))

############################
## clean up data frames

## growth data (height, leaf #, flower #)
W.plant$Date <- as.Date(W.plant$Date, "%m/%d/%y")
W.plant$Chamber <- factor(W.plant$Chamber, levels = c(60, 63, 62, 61))
W.plant$Plant <- as.factor(W.plant$Plant)

# set start planting date
W.plant0 <- W.plant %>% group_by(Plant, Chamber, Round) %>% filter(n()>2) %>% ungroup()
W.plant1 <- W.plant0 %>% filter(Round == "1")
start.date <- parse_date_time('2021-11-17 00:00:00', "%Y-%m-%d %H:%M:%S")
W.plant1$datenum <- as.numeric(difftime(parse_date_time(W.plant1$Date, 'ymd'),start.date))

W.plant2 <- W.plant0 %>% filter(Round == "2")
start.date <- parse_date_time('2022-01-20 00:00:00', "%Y-%m-%d %H:%M:%S")
W.plant2$datenum <- as.numeric(difftime(parse_date_time(W.plant2$Date, 'ymd'),start.date))

W.plant <- rbind(W.plant1, W.plant2)

#remove missing values
height <- W.plant %>% filter(!is.na(Height_cm))
leaves <- W.plant %>% filter(!is.na(Leaf_no))
flowers <- W.plant%>% filter(!is.na(Flower_no))

## first flower
W.FF$Date <- as.Date(W.FF$Date, "%m/%d/%y")
W.FF <- W.FF %>% filter(!is.na(Date))
W.FF$Chamber <- factor(W.FF$Chamber, levels = c(60, 63, 62, 61))


## flower area
W.area$Chamber <- factor(W.area$Chamber, levels = c(60, 63, 62, 61))
W.area$Plant <- as.factor(W.area$Plant)
W.area$CO2 <- as.factor(W.area$CO2)

## biomass
biomass$CO2 <- as.factor(biomass$CO2)
biomass$Plant <- as.factor(biomass$Plant)

############################

## height

# remove plants and dates with fewer than 3 observations, remove plants with zero height (seeds that had not germinated yet)
height <- height %>% group_by(Plant, Chamber, Round, Date) %>% filter(n()>2) %>% ungroup()
height <- height %>% filter (Height_cm != "0")

# set list of plant names
plants.height <- c("Borage", "Buckwheat","Red Clover", "Dandelion", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants.height) <- c("B", "BW", "C", "D", "LP", "N", "PP", "SA", "SF")

height$CO2 <- as.factor(height$CO2)

# plot height by week by CO2 treatments
ggplot(height, aes(x=datenum, y = Height_cm, color = CO2))+
  geom_point(aes(color = CO2), alpha = 0.4, size = 0.5)+
  geom_smooth(method = "loess", aes(fill = CO2))+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Height (cm) +/- se", x = "Weeks since planting")+
  scale_fill_manual(values = c("black", "#6c6cff"))+
  scale_color_manual(values = c("grey","cornflowerblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_x_continuous(breaks = c(14, 21, 28, 35, 42, 49, 56, 63, 70), 
                     labels = c("2", "3", "4", "5", "6", "7", "8", "9", "10"),
                     limits = c(14,70))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants.height),
             ncol = 3)

# convert date since planting to experimental week
height$week <- height$datenum/7

# full model
m.1 <- lmer(Height_cm ~ CO2*Plant + Round + week + (1|Chamber), data = height, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.1, plot =F))
summary(m.1)
anova(m.1)
# try ln-transforming data to improve normality
m.2 <- lmer(log(Height_cm) ~ CO2*Plant + Round + week + (1|Chamber), data = height, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.2, plot =F))
anova(m.1, m.2) # equivalent model fit, lower AIC in ln-transformed model - go with that one. 
anova(m.2) 

# compare each individual plant species separately
## Borage
h.B <- height %>% filter(Plant == "B")
m.B1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))

# try ln-transforming height data
m.B1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log, plot =F))

anova(m.B1, m.B1.log) # no difference, but much lower AIC for m.B1.log
summary(m.B1.log) # No significant difference by CO2 

## Buckwheat
h.BW <- height %>% filter(Plant == "BW")
m.BW1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1, plot =F))

# try ln-transforming height data
m.BW1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1.log, plot =F))
anova(m.BW1, m.BW1.log) # no difference but AIC is lower for ln-transformed data
summary(m.BW1.log) # CO2 is significant

## Clover
h.C <- height %>% filter(Plant == "C")
m.C1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1, plot =F))
m.C1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1.log, plot =F))
anova(m.C1, m.C1.log) # = explanation, lower AIC w/log-transformed
summary(m.C1.log) # no significant effect of round or of CO2 

# what if we remove round from the model
m.C1.log1 <- lmer(log(Height_cm) ~ CO2 + week + (1|Chamber), data = h.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1.log1, plot =F))
anova(m.C1.log, m.C1.log1) # = explanation, lower AIC w/o round
summary(m.C1.log1) # CO2 not significant 

## Dandelion
h.D <- height %>% filter(Plant == "D")
h.D <- h.D %>% filter(Round == "1") # Remove round 2 dandelion starts, flowering took took too long to continue in round 2

m.D1 <- lmer(Height_cm ~ CO2 + week + (1|Chamber), data = h.D, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.D1, plot =F))
m.D1.log <- lmer(log(Height_cm) ~ CO2 + week + (1|Chamber), data = h.D, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.D1.log, plot =F))
anova(m.D1, m.D1.log) # no difference in model fit, but AIC lower for ln-transformed data
summary(m.D1.log) # No significant difference by CO2 level, only week

#######
## Lacy Phacelia
######
h.LP <- height %>% filter(Plant == "LP")

m.LP1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
m.LP1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1.log, plot =F))
summary(m.LP1)
summary(m.LP1.log)
anova(m.LP1,m.LP1.log) # same as all others
AIC(m.LP1) # 6406.778
# no significant effect of CO2 but it's close. 

#########
## Nasturtium
##########
h.N <- height %>% filter(Plant == "N")

m.N1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
m.N1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1.log, plot =F))
summary(m.N1)
summary(m.N1.log)
anova(m.N1, m.N1.log)
AIC(m.N1) # 3245.996
# No significant difference by CO2 level

######
## PP
######
h.PP <- height %>% filter(Plant == "PP")
h.PP <- h.PP %>% filter(Round == "1")

m.PP1 <- lmer(Height_cm ~ CO2+week + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1, plot =F))
m.PP1.log <- lmer(log(Height_cm) ~ CO2 + week + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1.log, plot =F))
summary(m.PP1)
summary(m.PP1.log)
anova(m.PP1, m.PP1.log)
AIC(m.PP1) # 1754.496
# everything is significant

########
## SA
#######
h.SA <- height %>% filter(Plant == "SA")

m.SA1 <- lmer(Height_cm ~ CO2 + Round+ week + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
m.SA1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1.log, plot =F))
summary(m.SA1)
summary(m.SA1.log)
anova(m.SA1, m.SA1.log)
AIC(m.SA1) # 4477.474
# No significant difference by CO2 level, only round and date

#########
## SF
#########
h.SF <- height %>% filter(Plant == "SF")

m.SF1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
m.SF1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1.log, plot =F))
summary(m.SF1)
summary(m.SF1.log)
anova(m.SF1, m.SF1.log)
AIC(m.SF1) # 5843.401
# No significant difference by CO2 level, only round and date

# is it better w/o co2?
m.SF2 <- lmer(Height_cm ~ Round*Date + (1|Chamber), data = h.SF, REML=F)
summary(m.SF2)
anova(m.SF2)
AIC(m.SF2) # 5845.847

anova(m.SF1, m.SF2) # full model is better. 

##############
## leaves
##############

# remove any plant, date, round, chambers w/o 3+ observations
leaves <- leaves %>% group_by(Plant, Chamber, Round, Date) %>% filter(n()>2) %>% ungroup()

## summarize and plot
l.sum <- leaves %>% group_by(Plant, CO2, datenum) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Leaf_no, na.rm = T),
    sd = sd(Leaf_no, na.rm = T)
  )
l.sum$se <- l.sum$sd/sqrt(l.sum$count)

l.sum$CO2 <- factor(l.sum$CO2, levels = c("0", "1"))
l.sum <- l.sum %>% filter(datenum < 77)

ggplot(l.sum, aes(x=datenum, y = mean, color = CO2))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), width=0.2)+
  geom_line()+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Leaves +/- se", x = "Weeks since planting")+
  scale_color_manual(values = c("grey","navyblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_x_continuous(breaks = c(14, 21, 28, 35, 42, 49, 56, 63, 70), 
                     labels = c("2", "3", "4", "5", "6", "7", "8", "9", "10"))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants),
             ncol = 3)

leaves$CO2 <- as.factor(leaves$CO2)

ggplot(leaves, aes(x=datenum, y = Leaf_no, color = CO2))+
  geom_point(aes(color = CO2), alpha = 0.4, size = 0.5)+
  geom_smooth(method = "loess", aes(fill = CO2))+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Leaf # +/- se", x = "Weeks since planting")+
  scale_fill_manual(values = c("black", "#6c6cff"))+
  scale_color_manual(values = c("grey","cornflowerblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_x_continuous(breaks = c(14, 21, 28, 35, 42, 49, 56, 63, 70), 
                     labels = c("2", "3", "4", "5", "6", "7", "8", "9", "10"),
                     limits = c(14,70))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants),
             ncol = 3)


## summarize and plot
l.sum <- leaves %>% group_by(Plant, Chamber, Round, Date) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Leaf_no, na.rm = T),
    sd = sd(Leaf_no, na.rm = T)
  )
l.sum$se <- l.sum$sd/sqrt(l.sum$count)

plants <- c("Borage", "Buckwheat","Red Clover", "Dandelion", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "D", "LP", "N", "PP", "SA", "SF")

l.sum$Round <-as.factor(l.sum$Round)
l.sum1 <- l.sum %>% filter(Round == "1")

ggplot(l.sum1, aes(x=Date, y = mean, shape = Chamber))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), width=0.2)+
  geom_line()+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Leaf # +/- se", x = "Date")+
  scale_shape_manual(values = c(0,15,1,16), labels = c("60 - eCO2", "63 - eCO2", "62 - aCO2", "61 - aCO2"))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants),
             ncol = 3)

l.sum2 <- l.sum %>% filter(Round == "2")

ggplot(l.sum2, aes(x=Date, y = mean, shape = Chamber))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), width=0.2)+
  geom_line()+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Height (cm) +/- se", x = "Date")+
  scale_shape_manual(values = c(0,15,1,16), labels = c("60 - eCO2", "63 - eCO2", "62 - aCO2", "61 - aCO2"))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants),
             ncol = 3)

leaves$week <- leaves$datenum/7

m.B1 <- lmer(Leaf_no ~ CO2*Plant + Round + week + (1|Chamber), data = leaves, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))
m.B1.log <- lmer(log(Leaf_no) ~ CO2*Plant + Round + week + (1|Chamber), data = leaves, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log, plot =F))
qqPlot(log(leaves$Leaf_no+1))
anova(m.B1, m.B1.log)
summary(m.B1)
anova(m.B1)
anova(m.B1.log)

########
## Borage
########
h.B <- leaves %>% filter(Plant == "B")

m.B1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = h.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))
m.B1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = h.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log, plot =F))
summary(m.B1)
summary(m.B1.log)
anova(m.B1, m.B1.log)
m.B1.log2 <- lmer(log(Leaf_no) ~CO2 + week + (1|Chamber), data = h.B, REML = F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log2, plot =F))
anova(m.B1.log, m.B1.log2)
summary(m.B1.log2)
AIC(m.B1) # 3407.825
# No significant difference by CO2 level, only round and week

######
## Buckwheat
######
h.BW <- leaves %>% filter(Plant == "BW")

m.BW1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = h.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1, plot =F))
m.BW1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = h.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1.log, plot =F))
summary(m.BW1.log)
summary(m.BW1)
anova(m.BW1, m.BW1.log)
anova(m.BW1.log)
AIC(m.BW1) # 4312.423
# only week is significant - no effect of CO2

#########
## Clover
#########
h.C <- leaves %>% filter(Plant == "C")

m.C1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = h.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1, plot =F))
m.C1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = h.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1.log, plot =F))
summary(m.C1)
summary(m.C1.log)
anova(m.C1, m.C1.log)
AIC(m.C1) # 6435.151
# everything is significant

#######
## Dandelion
######
h.D <- leaves %>% filter(Plant == "D")
h.D <- h.D %>% filter(Round == "1")

m.D1 <- lmer(Leaf_no ~ CO2 + week + (1|Chamber), data = h.D, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.D1, plot =F))
m.D1.log <- lmer(log(Leaf_no) ~ CO2 + week + (1|Chamber), data = h.D, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.D1.log, plot =F))
summary(m.D1.log)
summary(m.D1)
anova(m.D1, m.D1.log)
AIC(m.D1) # 2110.201
# week significant

#######
## Lacy Phacelia
#######
h.LP <- leaves %>% filter(Plant == "LP")

m.LP1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
m.LP1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1.log, plot =F))
summary(m.LP1)
summary(m.LP1.log)
anova(m.LP1, m.LP1.log)
m.LP1.log2 <- lmer(log(Leaf_no) ~ CO2 + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1.log2, plot =F))
anova(m.LP1.log2, m.LP1.log) #AIC is 0.1 less without Round. Going to leave it as is for now. 
AIC(m.LP1) # 5969.193
# everything is significant

########
## Nasturtium
#######
h.N <- leaves %>% filter(Plant == "N")

m.N1 <- lmer(Leaf_no ~ CO2+Round+week + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
m.N1.log <- lmer(log(Leaf_no) ~ CO2+Round+week + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1.log, plot =F))
summary(m.N1)
summary(m.N1.log)
anova(m.N1, m.N1.log)
AIC(m.N1) # 4120.195
# round, week, roundxweek significant

########
## PP
########
h.PP <- leaves %>% filter(Plant == "PP")
h.PP <- h.PP %>% filter(Round == "1")

m.PP1 <- lmer(Leaf_no ~ CO2 + week + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1, plot =F))
m.PP1.log <- lmer(log(Leaf_no) ~ CO2 + week + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1.log, plot =F))
summary(m.PP1)
summary(m.PP1.log)
anova(m.PP1, m.PP1.log)
AIC(m.PP1) # 2194.852
# week is significant

########
## SA
########
h.SA <- leaves %>% filter(Plant == "SA")

m.SA1 <- lmer(Leaf_no ~ CO2+Round+week + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
m.SA1.log <- lmer(log(Leaf_no) ~ CO2+Round+week + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1.log, plot =F))
summary(m.SA1)
summary(m.SA1.log)
anova(m.SA1, m.SA1.log)
AIC(m.SA1) # 6336.782
# No significant difference by CO2 level, only round and date and not of co2 x round, but everything else 

#########
## SF
########
h.SF <- leaves %>% filter(Plant == "SF")

m.SF1 <- lmer(Leaf_no ~ CO2+Round+week + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
m.SF1.log <- lmer(log(Leaf_no) ~ CO2+Round+week + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1.log, plot =F))
summary(m.SF1)
summary(m.SF1.log)
anova(m.SF1, m.SF1.log)
AIC(m.SF1) #3883.972
# CO2, Co2xRound, Co2xweek not significant

##########################
## flowers
##############
# remove any plant, date, round, chambers w/o 3+ observations
flowers <- flowers %>% group_by(Plant, Chamber, Round, Date) %>% filter(n()>2) %>% ungroup()
flowers$week <- flowers$datenum/7
flowers <- flowers %>% filter(Plant != "D")
# omit dandelions as they never flowered
hist(log(flowers$Flower_no), breaks = 50)
# ok lots of zeros, zero-inflated count data

## summarize and plot
f.sum <- leaves %>% group_by(Plant, CO2, datenum) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Flower_no, na.rm = T),
    sd = sd(Flower_no, na.rm = T)
  )
f.sum$se <- f.sum$sd/sqrt(f.sum$count)
f.sum <- f.sum %>% filter(Plant != "D")

f.sum$CO2 <- factor(f.sum$CO2, levels = c("0", "1"))
#f.sum <- f.sum %>% filter(datenum < 77)

plants1 <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants1) <- c("B", "BW", "C", "LP", "N", "PP", "SA", "SF")

# remove weeks with no flowers
write.csv(f.sum, "flower.summary.csv")

f.sum.short<-read.csv("flower.summary.csv", header = T)
f.sum.short$CO2 <- as.factor(f.sum.short$CO2)

flowers$CO2 <- as.factor(flowers$CO2)

# fig for flowers - used this one 
ggplot(flowers, aes(x=week, y = Flower_no, color = CO2))+
  geom_smooth(method = "loess", aes(fill = CO2))+
  theme_classic()+
  geom_point(aes(color = CO2), alpha = 0.4, size = 0.5)+
  theme(legend.position = "bottom")+
  labs(y="Flowers +/- se", x = "Weeks since planting")+
  scale_fill_manual(values = c("black", "#6c6cff"))+
  scale_color_manual(values = c("grey", "cornflowerblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants1),
             ncol = 3)



#scale_x_continuous(breaks = c(28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105), 
#labels = c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

# sunflower should really be shown as a proportion of flowers blooming than like this, make second 
# figure and merge in AI
# omit weeks where mean is zero
f.sum.short <- f.sum.short %>% filter(mean !=0)

ggplot(f.sum.short, aes(x=datenum, y = mean, color = CO2))+
  geom_smooth(method = "loess", aes(fill = CO2))+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Flowers +/- se", x = "Weeks since planting")+
  scale_fill_manual(values = c("black", "#6c6cff"))+
  scale_color_manual(values = c("grey", "navyblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_x_continuous(breaks = c(28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105), 
                     labels = c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants1),
             ncol = 3)

## summarize and plot
f.sum <- flowers %>% group_by(Plant, Chamber, Round, Date) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Flower_no, na.rm = T),
    sd = sd(Flower_no, na.rm = T)
  )
f.sum$se <- f.sum$sd/sqrt(f.sum$count)

plants <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "LP", "N", "PP", "SA", "SF")

f.sum$Round <-as.factor(f.sum$Round)
f.sum1 <- f.sum %>% filter(Round == "1")

ggplot(f.sum1, aes(x=Date, y = mean, shape = Chamber))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), width=0.2, color = "gray")+
  geom_line()+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Flower # +/- se", x = "Date")+
  scale_shape_manual(values = c(0,15,1,16), labels = c("60 - eCO2", "63 - eCO2", "62 - aCO2", "61 - aCO2"))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants),
             ncol = 3)

f.sum2 <- f.sum %>% filter(Round == "2")
#omit PP
f.sum2 <- f.sum2 %>% filter(Plant != "PP")
# omit the zero points way out on March 17 of SA
f.sum2 <- f.sum2 %>% filter(Date != "2022-03-17")

ggplot(f.sum2, aes(x=Date, y = mean, shape = Chamber))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), width=0.2, color = "gray")+
  geom_line()+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs(y="Flower (cm) +/- se", x = "Date")+
  scale_shape_manual(values = c(0,15,1,16), labels = c("60 - eCO2", "63 - eCO2", "62 - aCO2", "61 - aCO2"))+
  facet_wrap(vars(Plant), scales = "free",
             labeller = labeller(Plant = plants),
             ncol = 3)

m.B1 <- lmer(Flower_no+1 ~ CO2*Plant + Round + datenum + (1|Chamber), data = flowers, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))
summary(m.B1)
anova(m.B1)

#remove data where there are no flowers
flowers1 <- flowers %>% filter(Flower_no != 0)
hist(flowers1$Flower_no)
qqPlot(log(flowers1$Flower_no))

m.short <- glmer.nb(Flower_no ~ CO2*Plant + Round + datenum + (1|Chamber),
                    data = flowers)
plot(simulationOutput <- simulateResiduals(fittedModel=m.short, plot =F))
car::Anova(m.short)

# where zeros are included
m.B2 <- glmer.nb(Flower_no ~ CO2*Plant + Round + datenum + (1|Chamber),
                 data = flowers)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B2, plot =F))
#check for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m.B2) # no over dispersion
car::Anova(m.B2) # report this

######
## Borage
#######
h.B <- flowers %>% filter(Plant == "B")

hist(h.B$Flower_no) # zero-inflated
m.B2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                 data = h.B)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B2, plot =F))
#check for overdispersion
overdisp_fun(m.B2) # yes overdispersion - wait actual underdispersion
testDispersion(m.B2, alternative="less")
car::Anova(m.B2) # report this
AIC(m.B2) # lower AIC than lmer, so go with this. 

#######
## Buckwheat
#########
h.BW <- flowers %>% filter(Plant == "BW")

m.BW1 <- lmer(Flower_no ~ CO2+Round+datenum + (1|Chamber), data = h.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1, plot =F))
summary(m.BW1) #AIC 9860.9
anova(m.BW1)
# significant difference by CO2 level
hist(h.BW$Flower_no) # zero-inflated
m.BW2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = h.BW)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW2, plot =F))
#check for overdispersion
testDispersion(m.BW2, alternative="greater") # no over or under dispersion according to this test
car::Anova(m.BW2) # report this
AIC(m.BW2) # 7291.22 lower AIC than lmer, so go with this. 
summary(m.BW2)


########
## Clover
########
h.C <- flowers %>% filter(Plant == "C")

m.C1 <- lmer(Flower_no ~ CO2+Round+datenum + (1|Chamber), data = h.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1, plot =F))
summary(m.C1) #AIC 8705.2
anova(m.C1)
# sig effect of CO2 - fewer flowers in eco2
hist(h.C$Flower_no) # zero-inflated
m.C2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                 data = h.C)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C2, plot =F))
#check for overdispersion
testDispersion(m.C2, alternative="less") # underdispersion according to this test
car::Anova(m.C2) # report this
AIC(m.C2) # 6239.437 lower AIC than lmer, so go with this. 
summary(m.C2)

########
## Dandelion - never flowered
#########

#######
## Lacy Phacelia
#######
h.LP <- flowers %>% filter(Plant == "LP")

m.LP1 <- lmer(Flower_no ~ CO2+Round+datenum + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
summary(m.LP1) #AIC 17846.7
anova(m.LP1)
# sig effect of CO2 - fewer flowers in eco2
hist(h.LP$Flower_no) # zero-inflated
m.LP2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = h.LP)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP2, plot =F))
#check for overdispersion
testDispersion(m.LP2, alternative="greater") # no over or underdispersion according to this test
car::Anova(m.LP2) # report this
AIC(m.LP2) # 13485.57 lower AIC than lmer, so go with this. 
summary(m.LP2)
# No significant difference by CO2 level, only round and date

#########
## Nasturtium
#########
h.N <- flowers %>% filter(Plant == "N")

m.N1 <- lmer(Flower_no ~ CO2 + Round+ datenum + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
summary(m.N1) #AIC 6432.1
anova(m.N1)
# sig effect of CO2 - fewer flowers in eco2
hist(h.N$Flower_no) # zero-inflated
m.N2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                 data = h.N)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N2, plot =F))
#check for overdispersion
testDispersion(m.N2, alternative="less") # no over or underdispersion according to this test
car::Anova(m.N2) # report this
AIC(m.N2) # 5037.598 lower AIC than lmer, so go with this. 
summary(m.N2)

########
## PP
########
h.PP <- flowers %>% filter(Plant == "PP")
h.PP <- h.PP %>% filter(Round == "1")

m.PP1 <- lmer(Flower_no ~ CO2 + datenum + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1, plot =F))
summary(m.PP1) #AIC 851.6
anova(m.PP1)
# sig effect of CO2 - fewer flowers in eco2
hist(h.PP$Flower_no) # zero-inflated
m.PP2 <- glmer.nb(Flower_no ~ CO2 + datenum + (1|Chamber),
                  data = h.PP)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP2, plot =F))
#check for overdispersion
testDispersion(m.PP2, alternative="greater") # no over or underdispersion according to this test
car::Anova(m.PP2) # report this
AIC(m.PP2) # 130.2388 lower AIC than lmer, so go with this. 
summary(m.PP2)

#########
## SA
#########
h.SA <- flowers %>% filter(Plant == "SA")

m.SA1 <- lmer(Flower_no ~ CO2+ Round+datenum + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
summary(m.SA1) #AIC 13594.0
anova(m.SA1)
# sig effect of CO2 - fewer flowers in eco2
hist(h.SA$Flower_no) # zero-inflated
m.SA2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = h.SA)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA2, plot =F))
#check for overdispersion
testDispersion(m.SA2, alternative="greater") # no over or underdispersion according to this test
car::Anova(m.SA2) # report this
AIC(m.SA2) # 8349.376 lower AIC than lmer, so go with this. 
summary(m.SA2)

#######
## SF
#######
h.SF <- flowers %>% filter(Plant == "SF")

m.SF1 <- lmer(Flower_no ~ CO2+Round+datenum + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
summary(m.SF1) #AIC 461.4
anova(m.SF1)
# sig effect of CO2 - fewer flowers in eco2
hist(h.SF$Flower_no) # zero-inflated
m.SF2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = h.SF)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF2, plot =F))
#check for overdispersion
testDispersion(m.SF2, alternative="less") # no over or underdispersion according to this test
car::Anova(m.SF2) # report this
AIC(m.SA2) # 8349.376 lower AIC than lmer, so go with this. 
summary(m.SA2)
summary(m.SF1)
anova(m.SF1)
# No significant difference by CO2 level, only round and date

###############
## days to first flower
##############
# remove any plant, date, round, chambers w/o 3+ observations

W.FF <- W.FF %>% group_by(Plant, Chamber, Round) %>% filter(n()>2) %>% ungroup()
ff.1 <- W.FF %>% filter(Round == "1")
start.date <- parse_date_time('2021-11-17 00:00:00', "%Y-%m-%d %H:%M:%S")
ff.1$datenum <- as.numeric(difftime(parse_date_time(ff.1$Date, 'ymd'),start.date))

ff.2 <- W.FF %>% filter(Round == "2")
start.date <- parse_date_time('2022-01-20 00:00:00', "%Y-%m-%d %H:%M:%S")
ff.2$datenum <- as.numeric(difftime(parse_date_time(ff.2$Date, 'ymd'),start.date))

ff<- rbind(ff.1, ff.2)

plants <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia", "Nasturtium", "Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "LP", "N", "PP", "SA", "SF")

## summarize and plot
ff.sum <- ff %>% group_by(Plant, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(datenum, na.rm = T),
    sd = sd(datenum, na.rm = T)
  )
ff.sum$se <- ff.sum$sd/sqrt(ff.sum$count)

ff.sum$CO2 <- factor(ff.sum$CO2, levels = c(0, 1))

ggplot(ff.sum, aes(x=Plant, y = mean, color = CO2))+
  geom_point(position=position_dodge(0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), width=0.2, position=position_dodge(0.75))+
  theme_classic()+
  theme(legend.position = c(0.1, 0.9),
        axis.text.x = element_text(angle = 90))+
  labs(y="Days since planting ± se")+
  scale_color_manual(values = c("grey","navyblue"), 
                     labels = c("aCO2", "eCO2"),
                     name = "Treatment")+
  scale_x_discrete(labels = c(plants))

W.FF$CO2 <- as.factor(W.FF$CO2)
W.FF$Plant <- as.factor(W.FF$Plant)

ggplot(ff, aes(x=Plant, y = log10(datenum), fill = CO2))+
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
               position=position_dodge(0.9))+
  labs(y = "log(Days since planting)")

## Compare all 
hist(log(ff$datenum))
m.B1 <- lmer(datenum ~ CO2*Plant+Round + (1|Chamber), data = ff, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2*Plant+Round + (1|Chamber), data = ff, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but log(datenum) with a substantially lower AIC. 

######
## Borage
#######
ff.B <- ff %>% filter(Plant == "B")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
anova(m.B2)
summary(m.B2)
anova(m.B1, m.B2) # no difference, but AIC lower when log-transformed
# no significance 

######
## Buckwheat
######
ff.B <- ff %>% filter(Plant == "BW")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference between model coverage, but lower AIC when log-transformed

#######
## Clover
########
ff.B <- ff %>% filter(Plant == "C")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, AIC lower for log transformed

#######
## Lacy Phacelia
########
ff.B <- ff %>% filter(Plant == "LP")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no sig dif, but lower AIC when log-transformed

#######
## Nasturtium
#######
ff.B <- ff %>% filter(Plant == "N")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no sig dif between models, but log transformed w/lower AIC

#######
## Partridge Pea
#####
ff.B <- ff %>% filter(Plant == "PP")
ff.B <- ff.B %>% filter(Round == "1")

m.B1 <- lmer(datenum ~ CO2 + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no sig dif, but log-transformed w/lower AIC

########
## Sweet alyssum
########
ff.B <- ff %>% filter(Plant == "SA")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference but lower AIC when log transformed

########
## Sunflower
######
ff.B <- ff %>% filter(Plant == "SF")

m.B1 <- lmer(datenum ~ CO2+Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference by log transformed model w/lower AIC

######################
## flower diameter
#####################
# remove NAs'
W.area <- W.area %>% filter(!is.na(Plant))

# remove any plant, date, round, chambers w/o 3+ observations
area <- W.area %>% group_by(Plant, CO2) %>% filter(n()>2) %>% ungroup()

## summarize and plot
a.sum <- area %>% group_by(Plant, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(log(Mean), na.rm = T),
    sd = sd(log(Mean), na.rm = T)
  )
a.sum$se <- a.sum$sd/sqrt(a.sum$count)
a.sum$CO2 <- factor(a.sum$CO2, levels = c(0, 1))

plants <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "LP", "N", "PP", "SA", "SF")

ggplot(a.sum, aes(x=Plant, y = mean, color = CO2))+
  geom_point(position=position_dodge(0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), 
                width=0.3, 
                position=position_dodge(0.75))+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))+
  labs(y="ln(Flower diameter) ± se", x = "Plant species")+
  scale_color_manual(values = c("grey","navyblue"), 
                     labels = c("aCO2", "eCO2"))+
  scale_x_discrete(labels = plants)

ggplot(area, aes(x=Plant, y = log(Mean), fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.3, alpha = 0.3, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1, size = .5),
        axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))+
  labs(y = "ln(Flower diameter)")


## Compare all 
qqPlot(log(area$Mean))
hist(log(area$Mean))
m.B1 <- lmer(Mean ~ CO2*Plant + Round + (1|Chamber), data = area, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2*Plant + Round + (1|Chamber), data = area, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but log-transformed has lower AIC

#######
## Borage
########
ff.B <- area %>% filter(Plant == "B")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

########
## Buckwheat
######
ff.B <- area %>% filter(Plant == "BW")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

##########
## Clover
########
ff.B <- area %>% filter(Plant == "C")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

######
## Lacy Phacelia
######
ff.B <- area %>% filter(Plant == "LP")

m.B1 <- lmer(Mean ~ CO2+Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

#####
## Nasturtium
######
ff.B <- area %>% filter(Plant == "N")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

#########
## Partridge Pea
#########
ff.B <- area %>% filter(Plant == "PP")
ff.B <- ff.B %>% filter(Round == "1")

m.B1 <- lmer(Mean ~ CO2 + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

#########
## Sweet alyssum
##########
ff.B <- area %>% filter(Plant == "SA")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

##########
## Sunflower
###########
ff.B <- area %>% filter(Plant == "SF")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
summary(m.B1)
anova(m.B1)
summary(m.B2)
anova(m.B2)
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed

# Round 1
ff.B1 <- ff.B %>% filter(Round == "1")
t.test(Mean~CO2, data = ff.B1)
# t(21.998) = -2.0989, p = 0.04753 -  significant
# mean aCO2 = 97.39394, eCO2 = 110.58974, eCO2 flowers are bigger than aCO2
# compare by chamber
res.aov<-aov(Mean~Chamber, data = ff.B1)
summary(res.aov) # no significant difference by chamber F(3)=1.938, p = 0.156
# Round 2
ff.B1 <- ff.B %>% filter(Round == "2")
t.test(Mean~CO2, data = ff.B1)
# t(48.377) = -0.72369, p = 0.4727 -  not significant 
res.aov<-aov(Mean~Chamber, data = ff.B1)
summary(res.aov) # no significant difference by chamber F(3)=0.978, p = 0.415


##############
## biomass
#############
# only from round 1
BM <- biomass %>% group_by(Plant, CO2) %>% filter(n()>2) %>% ungroup()
# ok all are good to go 

BM.sum <- BM %>% group_by(Plant, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Plant.biomass, na.rm = T),
    sd = sd(Plant.biomass, na.rm = T)
  )
BM.sum$se <- BM.sum$sd/sqrt(BM.sum$count)

# plot it
plants <- c("Borage", "Buckwheat","Red Clover", "Dandelion", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants) <- c("B", "BW", "C", "D", "LP", "N", "PP", "SA", "SF")

ggplot(BM, aes(x=Plant, y = log(Plant.biomass), fill = CO2))+
  geom_violin()+
  geom_point(position=position_jitterdodge(), size = 0.5, alpha = 0.5, aes(color=CO2))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_color_manual(values = c("black", "navy"),
                     labels = c("aCO2", "eCO2"),
                     name = "CO2 Treatment")+
  scale_x_discrete(labels = plants)#+
theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust=1, size = .2),
      axis.title.x=element_blank())+
  stat_summary(fun = "mean",
               fun.args = list(mult = 1),
               geom="crossbar", 
               color = "black",
               width=.75,
               position=position_dodge(0.9))+
  labs(y="ln(Above-ground dry biomass)")

# compare all 
hist(log(BM$Plant.biomass))
str(BM)
BM$Chamber <- as.factor(BM$Chamber)
BM <- BM %>% filter(Plant.biomass != 0)
m.B <- lmer(log(Plant.biomass) ~ CO2*Plant+(1|Chamber), data = BM, REML = F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B, plot = F))
summary(m.B)
anova(m.B)
AIC(m.B) # 823.3703 
# only sig. effect of plant, can quickly compare w/t-test

# Borage
b <- BM %>% filter(Plant == "B")
t.test(Plant.biomass~CO2, data = b)
# not significant 

# Buckwheat
bw <- BM %>% filter(Plant == "BW")
t.test(Plant.biomass~CO2, data = bw)
# not significant 

# Clover
c <- BM %>% filter(Plant == "C")
t.test(Plant.biomass~CO2, data = c)
# not significant 

# Dandelion
d <- BM %>% filter(Plant == "D")
t.test(Plant.biomass~CO2, data = d)
# not significant 

# Lacy Phacelia
l <- BM %>% filter(Plant == "LP")
t.test(Plant.biomass~CO2, data = l)
# not significant 

# Nasturtium
n <- BM %>% filter(Plant == "N")
t.test(Plant.biomass~CO2, data = n)
# not significant 

# Partridge pea
p <- BM %>% filter(Plant == "PP")
t.test(Plant.biomass~CO2, data = p)
# significant

# Sweet Alyssum
sa <- BM %>% filter(Plant == "SA")
t.test(Plant.biomass~CO2, data = sa)
# not significant 

# Sunflower
sf <- BM %>% filter(Plant == "SF")
t.test(Plant.biomass~CO2, data = sf)
# not significant 

###################
## Harvard data 
##################

##################
## pollen - %N & C:N ratio
#################
# summarize and plot data 

# remove plants, rounds, and chambers w/fewer than 3 samples
H.short <- H.pollen %>% group_by(Plant, CO2, Round) %>% filter(n()>2) %>% ungroup()
# Takes us down to 110 observations

H.sum <- H.short %>% group_by(Plant, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(N, na.rm = T),
    sd = sd(N, na.rm = T)
  )
H.sum$se <- H.sum$sd/sqrt(H.sum$count)
H.sum$Variable <- "N"

# plot it
plants <- c("Buckwheat","Melon", "Partridge pea", "Poppy", "Squash", "Sunflower", "Tomatillo", "Tomato")
names(plants) <- c("buckwheat", "melon", "partridge pea", "poppy", "squash", "sunflower", "tomatillo", "tomato")

ggplot(H.sum, aes(x=Plant, y = mean, shape = CO2))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = c(.1,.8))+
  labs(y="%N +/- se", x = "Plant species")+
  scale_shape_manual(values = c(1,16), labels = c("aCO2", "eCO2"), name = "Treatment")

H.short$ratio <- H.short$C/H.short$N

H.sum1 <- H.short %>% group_by(Plant, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(ratio, na.rm = T),
    sd = sd(ratio, na.rm = T)
  )
H.sum1$se <- H.sum1$sd/sqrt(H.sum1$count)
H.sum1$Variable <- "C:N"

ggplot(H.sum1, aes(x=Plant, y = mean, shape = CO2))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = c(.9,.8))+
  labs(y="C:N ratio +/- se", x = "Plant species")+
  scale_shape_manual(values = c(1,16), labels = c("aCO2", "eCO2"), name = "Treatment")

H.sum2 <- H.short %>% group_by(Plant, CO2) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(C, na.rm = T),
    sd = sd(C, na.rm = T)
  )
H.sum2$se <- H.sum2$sd/sqrt(H.sum2$count)
H.sum2$Variable <- "C"

ggplot(H.sum2, aes(x=Plant, y = mean, shape = CO2))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = c(.9,.2))+
  labs(y="%C +/- se", x = "Plant species")+
  scale_shape_manual(values = c(1,16), labels = c("aCO2", "eCO2"), name = "Treatment")

H.p <- rbind(H.sum, H.sum1, H.sum2)

H.p$Variable <- factor(H.p$Variable, levels = c("N", "C", "C:N"))

ggplot(H.p, aes(x=Plant, y = mean, color = CO2))+
  geom_point(position = position_dodge(w = 0.75), size = 2)+
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(w = 0.75), width=0.2)+
  theme_classic()+
  theme(legend.position = "bottom")+
  labs( y = "Mean +/- s.e.")+
  scale_color_manual(values = c("grey","navyblue"), 
                     labels = c("aCO2", "eCO2"))+
  scale_x_discrete(labels = c(plants))+
  facet_wrap(vars(Variable), scales = "free", 
             ncol = 1)


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

ggarrange(N, C, CN, nrow=3, ncol = 1, align = "h")

ggplot(H.short, aes(x= Plant, y = N, fill = CO2))+
  geom_violin()+
  geom_boxplot(width=0.1, aes(fill = CO2), position=position_dodge(.9))+
  theme_classic()+
  scale_fill_manual(values = c("grey", "cornflowerblue"),
                    labels = c("aCO2", "eCO2"),
                    name = "CO2 Treatment")+
  scale_x_discrete(labels = plants)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())


# N
# modeling to compare values:
m.B <- lm(N ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
AIC(m.B) #124.2365
# significant effect of Round, Plant, RoundxPlant, so compare each plant separately!
# AIC

m.B2 <- lmer(N~CO2*Plant + (1|Round), data = H.short)
summary(m.B2)
anova(m.B2)
AIC(m.B2) #154.5247 - ok lm better fit than lmer
# only plant significant


# C:N ratio
m.B <- lm(ratio ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
AIC(m.B) #194.1173
# sig effect of Round, Plant, round x plant, co2 x plant

# %C
m.B <- lm(C ~ CO2*Plant + Round, data = H.short)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
AIC(m.B) # 254.8034
# sig effect of C, Plant species, so compare each plant species separately

# is it a better model fit without Round?
m.B1 <- lm(C ~ CO2*Plant, data = H.short)
summary(m.B1)
anova(m.B1)
AIC(m.B1) # 271.9053 - no

anova(m.B, m.B1) # not significantly better. 

# Melon - only R1 
# %N (proxy for pollen protein)
M.po <- H.short %>% filter(Plant == "melon")
# only one factor being compared here, do not need to use a model, instead compare using t-test
t.test(N~CO2, M.po)
# not significant
# t(2.6825) = 1.3438, p = 0.2814

# C:N ratio
t.test(ratio~CO2, M.po)
# not significant
# t(2.5884) = -1.4126, p = 0.2661

# C
t.test(C~CO2, M.po)
# not significant
# t(2.9693) = 1.3939, p = 0.2585


# Poppy - only R1 
# %N (proxy for pollen protein)
P.po <- H.short %>% filter(Plant == "poppy")
P.po <- P.po %>% filter(Round =="1")
# only one factor being compared here, do not need to use a model, instead compare using t-test
t.test(N~CO2, P.po)
# not significant
# t(13.603) = -0.96607, p = 0.3509

# C:N ratio
t.test(ratio~CO2, P.po)
# not significant
# t(15.01) = 1.0785, p = 0.2978

# C
t.test(C~CO2, P.po)
# not significant
# t(15.01) = 1.0785, p = 0.2978

# Squash - both rounds 
# %N (proxy for pollen protein)
S.po <- H.short %>% filter(Plant == "squash")
m.B <- lm(N ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
# nothing significant
m.B <- lm(ratio ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
# nothing significant
m.B <- lm(C ~ CO2 + Round, data = S.po)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B, plot =F))
summary(m.B)
anova(m.B)
# Round significant

# compare between rounds
# round 1
S.po1 <- S.po %>% filter(Round == 1)
t.test(C~CO2, S.po1)
# not significant
# t(4.3331) = -0.85743, p = 0.4361
# round 2
S.po2 <- S.po %>% filter(Round == 2)
t.test(C~CO2, S.po2)
# not significant
# t(8.0162) = 0.47031, p = 0.6507


# Sunflower 
# %N (proxy for pollen protein)
SF.po <- H.short %>% filter(Plant == "sunflower")
SF.po <- SF.po %>% filter(Round == "1")

# %N
t.test(N~CO2, SF.po)
# not significant, but barely 
# t(10.989) = 2.1374, p = 0.05589

# C:N ratio
t.test(ratio~CO2, SF.po)
# not significant
# t(10.671) = -2.0801, p = 0.06246

# C
t.test(C~CO2, SF.po)
# not significant
# t(10.933) = 0.38058, p = 0.7108



# Tomatillo 
# %N (proxy for pollen protein)
TT.po <- H.short %>% filter(Plant == "tomatillo")

# %N
t.test(N~CO2, TT.po)
# not significant, but barely 
# t(11.396) = -1.0187, p = 0.3295

# C:N ratio
t.test(ratio~CO2, TT.po)
# not significant
# t(12.141) = 1.5465, p = 0.1476

# C
t.test(C~CO2, TT.po)
# not significant
# t(11.777) = 1.2489, p = 0.2359



# Tomato 
# %N (proxy for pollen protein)
T.po <- H.short %>% filter(Plant == "tomato")

# %N
t.test(N~CO2, T.po)
# not significant, but barely 
# t(9.9433) = -0.18633, p = 0.8559

# C:N ratio
t.test(ratio~CO2, T.po)
# not significant
# t(9.9864) = 0.15093, p = 0.883

# C
t.test(C~CO2, T.po)
# not significant
# t(7.9272) = -0.19124, p = 0.8531

