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
W.FF$CO2 <- as.factor(W.FF$CO2)
W.FF$Plant <- as.factor(W.FF$Plant)


## flower area
W.area$Chamber <- factor(W.area$Chamber, levels = c(60, 63, 62, 61))
W.area$Plant <- as.factor(W.area$Plant)
W.area$CO2 <- as.factor(W.area$CO2)

## biomass
biomass$CO2 <- as.factor(biomass$CO2)
biomass$Plant <- as.factor(biomass$Plant)

###################

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

## Lacy Phacelia
h.LP <- height %>% filter(Plant == "LP")

m.LP1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
m.LP1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1.log, plot =F))
anova(m.LP1,m.LP1.log) # = explanation, lower AIC when ln-transformed
summary(m.LP1.log) # no significant effect of CO2 

## Nasturtium
h.N <- height %>% filter(Plant == "N")

m.N1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
m.N1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1.log, plot =F))
anova(m.N1, m.N1.log) # = explanation, lower AIC when ln-transformed
summary(m.N1.log)# No significant difference by CO2 level

## Partridge pea
h.PP <- height %>% filter(Plant == "PP")
h.PP <- h.PP %>% filter(Round == "1") # high mortality in round 2, did not continue w/this plant in round 2

m.PP1 <- lmer(Height_cm ~ CO2+week + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1, plot =F))
m.PP1.log <- lmer(log(Height_cm) ~ CO2 + week + (1|Chamber), data = h.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1.log, plot =F))
anova(m.PP1, m.PP1.log) # = explanation, lower AIC when ln-transformed
summary(m.PP1.log) # only week is significant

## Sweet alyssum
h.SA <- height %>% filter(Plant == "SA")

m.SA1 <- lmer(Height_cm ~ CO2 + Round+ week + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
m.SA1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1.log, plot =F))
anova(m.SA1, m.SA1.log) # = explanation, lower AIC when ln-transformed
summary(m.SA1.log) # No significant difference by CO2 level, only round and date

## Sunflower 
h.SF <- height %>% filter(Plant == "SF")

m.SF1 <- lmer(Height_cm ~ CO2 + Round + week + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
m.SF1.log <- lmer(log(Height_cm) ~ CO2 + Round + week + (1|Chamber), data = h.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1.log, plot =F))
anova(m.SF1, m.SF1.log)
summary(m.SF1.log) # No significant difference by CO2 level, only round and date

##############
## leaves

# remove any plant, date, round, chambers w/o 3+ observations
leaves <- leaves %>% group_by(Plant, Chamber, Round, Date) %>% filter(n()>2) %>% ungroup()
leaves$CO2 <- as.factor(leaves$CO2)
plants.leaves <- c("Borage", "Buckwheat","Red Clover", "Dandelion", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants.leaves) <- c("B", "BW", "C", "D", "LP", "N", "PP", "SA", "SF")

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
             labeller = labeller(Plant = plants.leaves),
             ncol = 3)

leaves$week <- leaves$datenum/7

# full model looking at leaf number, week, plant species, and CO2 treament 
m.B1 <- lmer(Leaf_no ~ CO2*Plant + Round + week + (1|Chamber), data = leaves, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))
# see if ln-transformation is a better fit
m.B1.log <- lmer(log(Leaf_no) ~ CO2*Plant + Round + week + (1|Chamber), data = leaves, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log, plot =F))
anova(m.B1, m.B1.log) # yes, AIC lower after ln-transformation, = explanation of variation in data
anova(m.B1.log) # CO2 not significant, but Plant species and Co2 x Plant species, so compare each species alone

## Borage
l.B <- leaves %>% filter(Plant == "B")

m.B1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = l.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1, plot =F))
m.B1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = l.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log, plot =F))
anova(m.B1, m.B1.log) # equal explanation, lower AIC when ln-transformed
summary(m.B1.log) # round not significant, compare w/o round

m.B1.log2 <- lmer(log(Leaf_no) ~CO2 + week + (1|Chamber), data = l.B, REML = F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log2, plot =F))
anova(m.B1.log, m.B1.log2) # model slightly lower AIC w/o Round, go with that. 
summary(m.B1.log2) # No significant difference by CO2 level, only week

## Buckwheat
l.BW <- leaves %>% filter(Plant == "BW")

m.BW1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = l.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1, plot =F))
m.BW1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = l.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW1.log, plot =F))
anova(m.BW1, m.BW1.log) # ln-transformed w/lower AIC
summary(m.BW1.log) # no effect of CO2

## Clover
l.C <- leaves %>% filter(Plant == "C")

m.C1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = l.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1, plot =F))
m.C1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = l.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C1.log, plot =F))
anova(m.C1, m.C1.log) # ln-transformed w/lower AIC
summary(m.C1.log) # CO2 not significant 

## Dandelion
l.D <- leaves %>% filter(Plant == "D")
l.D <- l.D %>% filter(Round == "1") # only had plants in round 1

m.D1 <- lmer(Leaf_no ~ CO2 + week + (1|Chamber), data = l.D, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.D1, plot =F))
m.D1.log <- lmer(log(Leaf_no) ~ CO2 + week + (1|Chamber), data = l.D, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.D1.log, plot =F))
anova(m.D1, m.D1.log) #ln-transformed lower AIC
summary(m.D1.log) # week significant

## Lacy Phacelia
l.LP <- leaves %>% filter(Plant == "LP")

m.LP1 <- lmer(Leaf_no ~ CO2 + Round + week + (1|Chamber), data = l.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1, plot =F))
m.LP1.log <- lmer(log(Leaf_no) ~ CO2 + Round + week + (1|Chamber), data = l.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1.log, plot =F))
anova(m.LP1, m.LP1.log) #ln-transformed lower AIC
summary(m.LP1.log) # round not significant compare models w/o round

m.LP1.log2 <- lmer(log(Leaf_no) ~ CO2 + week + (1|Chamber), data = l.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP1.log2, plot =F))
anova(m.LP1.log2, m.LP1.log) #AIC is 0.1 less without Round. Going to leave it as is for now. 

## Nasturtium
l.N <- leaves %>% filter(Plant == "N")

m.N1 <- lmer(Leaf_no ~ CO2+Round+week + (1|Chamber), data = l.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1, plot =F))
m.N1.log <- lmer(log(Leaf_no) ~ CO2+Round+week + (1|Chamber), data = l.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1.log, plot =F))
anova(m.N1, m.N1.log) #ln-transformed lower AIC
summary(m.N1.log) # Co2, week significant

# try without round
m.N1.log2 <- lmer(log(Leaf_no) ~ CO2+week + (1|Chamber), data = l.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N1.log2, plot =F))
anova(m.N1.log2, m.N1.log) # lower AIC w/o round
summary(m.N1.log2)

## Partridge Pea
l.PP <- leaves %>% filter(Plant == "PP")
l.PP <- l.PP %>% filter(Round == "1") # only had plants in round 1

m.PP1 <- lmer(Leaf_no ~ CO2 + week + (1|Chamber), data = l.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1, plot =F))
m.PP1.log <- lmer(log(Leaf_no) ~ CO2 + week + (1|Chamber), data = l.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP1.log, plot =F))
anova(m.PP1, m.PP1.log) # lower AIC w/ln-transformation
summary(m.PP1.log) # co2 not significant 

## Sweet alyssum 
l.SA <- leaves %>% filter(Plant == "SA")

m.SA1 <- lmer(Leaf_no ~ CO2+Round+week + (1|Chamber), data = l.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1, plot =F))
m.SA1.log <- lmer(log(Leaf_no) ~ CO2+Round+week + (1|Chamber), data = l.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA1.log, plot =F))
anova(m.SA1, m.SA1.log) # lower AIC w/ln-transformation
summary(m.SA1.log)n# No significant difference by CO2 level

## Sunflower
l.SF <- leaves %>% filter(Plant == "SF")

m.SF1 <- lmer(Leaf_no ~ CO2+Round+week + (1|Chamber), data = l.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1, plot =F))
m.SF1.log <- lmer(log(Leaf_no) ~ CO2+Round+week + (1|Chamber), data = l.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF1.log, plot =F))
anova(m.SF1, m.SF1.log) # lower AIC w/ln-transformation
summary(m.SF1.log) # CO2 not significant

###################
## flowers

# remove any plant, date, round, chambers w/o 3+ observations
flowers <- flowers %>% group_by(Plant, Chamber, Round, Date) %>% filter(n()>2) %>% ungroup()
flowers$week <- flowers$datenum/7
flowers <- flowers %>% filter(Plant != "D") # omit dandelions as they never flowered

plants.flowers <- c("Borage", "Buckwheat","Red Clover", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants.flowers) <- c("B", "BW", "C", "LP", "N", "PP", "SA", "SF")

flowers$CO2 <- as.factor(flowers$CO2)

# plot flowers
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
             labeller = labeller(Plant = plants.flowers),
             ncol = 3)

# full model w/all plant species
m.3 <- glmer.nb(Flower_no ~ CO2*Plant+Round+datenum+ (1|Chamber),
                 data = flowers)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B2, plot =F))
car::Anova(m.3) # report this

# compare each plant species separately 
## Borage
f.B <- flowers %>% filter(Plant == "B")

m.B2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                 data = f.B)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B2, plot =F))
car::Anova(m.B2) # report this

## Buckwheat
f.BW <- flowers %>% filter(Plant == "BW")

m.BW2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = f.BW)
plot(simulationOutput <- simulateResiduals(fittedModel=m.BW2, plot =F))
car::Anova(m.BW2) # report this

## Clover
f.C <- flowers %>% filter(Plant == "C")

m.C2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                 data = f.C)
plot(simulationOutput <- simulateResiduals(fittedModel=m.C2, plot =F))
car::Anova(m.C2) # report this

## Lacy Phacelia
f.LP <- flowers %>% filter(Plant == "LP")

m.LP2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = f.LP)
plot(simulationOutput <- simulateResiduals(fittedModel=m.LP2, plot =F))
car::Anova(m.LP2) # report this

## Nasturtium
f.N <- flowers %>% filter(Plant == "N")

m.N2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                 data = f.N)
plot(simulationOutput <- simulateResiduals(fittedModel=m.N2, plot =F))
car::Anova(m.N2) # report this

## Partridge pea
f.PP <- flowers %>% filter(Plant == "PP")
f.PP <- f.PP %>% filter(Round == "1")

m.PP2 <- glmer.nb(Flower_no ~ CO2 + datenum + (1|Chamber),
                  data = h.PP)
plot(simulationOutput <- simulateResiduals(fittedModel=m.PP2, plot =F))
car::Anova(m.PP2) # report this

## Sweet alyssum 
f.SA <- flowers %>% filter(Plant == "SA")

m.SA2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = f.SA)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SA2, plot =F))
car::Anova(m.SA2) # report this

## Sunflower
f.SF <- flowers %>% filter(Plant == "SF")

m.SF2 <- glmer.nb(Flower_no ~ CO2 + Round + datenum + (1|Chamber),
                  data = f.SF)
plot(simulationOutput <- simulateResiduals(fittedModel=m.SF2, plot =F))
car::Anova(m.SF2) # report this

############
## days to first flower

# remove any plant, date, round, chambers w/o 3+ observations
W.FF <- W.FF %>% group_by(Plant, Chamber, Round) %>% filter(n()>2) %>% ungroup()
ff.1 <- W.FF %>% filter(Round == "1")
start.date <- parse_date_time('2021-11-17 00:00:00', "%Y-%m-%d %H:%M:%S")
ff.1$datenum <- as.numeric(difftime(parse_date_time(ff.1$Date, 'ymd'),start.date))

ff.2 <- W.FF %>% filter(Round == "2")
start.date <- parse_date_time('2022-01-20 00:00:00', "%Y-%m-%d %H:%M:%S")
ff.2$datenum <- as.numeric(difftime(parse_date_time(ff.2$Date, 'ymd'),start.date))

ff<- rbind(ff.1, ff.2)

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
  scale_x_discrete(labels = plants.flowers)+
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

## Compare all plant species together
m.4 <- lmer(datenum ~ CO2*Plant+Round + (1|Chamber), data = ff, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.4, plot = F))
m.4.log <- lmer(log(datenum) ~ CO2*Plant+Round + (1|Chamber), data = ff, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.4.log, plot = F))
anova(m.4, m.4.log) # no difference, but log(datenum) with a substantially lower AIC.
summary(m.4.log)
anova(m.4.log)
 
# compare each plant species independently
## Borage
ff.B <- ff %>% filter(Plant == "B")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but AIC lower when log-transformed
summary(m.B2)

## Buckwheat
ff.BW <- ff %>% filter(Plant == "BW")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference between model coverage, but lower AIC when log-transformed
summary(m.B2)

## Clover
ff.C <- ff %>% filter(Plant == "C")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, AIC lower for log transformed
summary(m.B2)

## Lacy Phacelia
ff.LP <- ff %>% filter(Plant == "LP")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no sig dif, but lower AIC when log-transformed
summary(m.B2)

## Nasturtium
ff.N <- ff %>% filter(Plant == "N")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no sig dif between models, but log transformed w/lower AIC
summary(m.B2)

## Partridge Pea
#####
ff.PP <- ff %>% filter(Plant == "PP")
ff.PP <- ff.PP %>% filter(Round == "1")

m.B1 <- lmer(datenum ~ CO2 + (1|Chamber), data = ff.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + (1|Chamber), data = ff.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no sig dif, but log-transformed w/lower AIC
summary(m.B2)

## Sweet alyssum
ff.SA <- ff %>% filter(Plant == "SA")

m.B1 <- lmer(datenum ~ CO2 + Round + (1|Chamber), data = ff.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference but lower AIC when log transformed
summary(m.B2)

## Sunflower
ff.SF <- ff %>% filter(Plant == "SF")

m.B1 <- lmer(datenum ~ CO2+Round + (1|Chamber), data = ff.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(datenum) ~ CO2 + Round + (1|Chamber), data = ff.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference by log transformed model w/lower AIC
summary(m.B2)

######################
## flower diameter

# remove NAs'
W.area <- W.area %>% filter(!is.na(Plant))
# remove any plant, date, round, chambers w/o 3+ observations
area <- W.area %>% group_by(Plant, CO2) %>% filter(n()>2) %>% ungroup()

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
  scale_x_discrete(labels = plants.flowers)+
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


## Compare all plant species together
m.B1 <- lmer(Mean ~ CO2*Plant + Round + (1|Chamber), data = area, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2*Plant + Round + (1|Chamber), data = area, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but log-transformed has lower AIC
summary(m.B2)

# compare each plant species separately
## Borage
fa.B <- area %>% filter(Plant == "B")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = fa.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.B, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Buckwheat
######
fa.BW <- area %>% filter(Plant == "BW")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = fa.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.BW, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Clover
fa.C <- area %>% filter(Plant == "C")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = fa.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.C, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Lacy Phacelia
fa.LP <- area %>% filter(Plant == "LP")

m.B1 <- lmer(Mean ~ CO2+Round + (1|Chamber), data = fa.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.LP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Nasturtium
fa.N <- area %>% filter(Plant == "N")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = fa.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.N, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Partridge Pea
fa.PP <- area %>% filter(Plant == "PP")
fa.PP <- fa.PP %>% filter(Round == "1")

m.B1 <- lmer(Mean ~ CO2 + (1|Chamber), data = fa.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + (1|Chamber), data = fa.PP, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Sweet alyssum
fa.SA <- area %>% filter(Plant == "SA")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = fa.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.SA, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

## Sunflower
fa.SF <- area %>% filter(Plant == "SF")

m.B1 <- lmer(Mean ~ CO2 + Round + (1|Chamber), data = fa.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B1, plot = F))
m.B2 <- lmer(log(Mean) ~ CO2 + Round + (1|Chamber), data = fa.SF, REML=F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B1, m.B2) # no difference, but lower AIC when log-transformed
summary(m.B2)

##############
## biomass
# only from round 1
BM <- biomass %>% group_by(Plant, CO2) %>% filter(n()>2) %>% ungroup()
# ok all are good to go 

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
  scale_x_discrete(labels = plants.height)+
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

# compare all plant species together
BM$Chamber <- as.factor(BM$Chamber)
BM <- BM %>% filter(Plant.biomass != 0)
m.B <- lmer(log(Plant.biomass) ~ CO2*Plant+(1|Chamber), data = BM, REML = F)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B, plot = F))
summary(m.B)
anova(m.B)

# compare each plant species independently with t-tests
# Borage
b <- BM %>% filter(Plant == "B")
t.test(Plant.biomass~CO2, data = b) # not significant 

# Buckwheat
bw <- BM %>% filter(Plant == "BW")
t.test(Plant.biomass~CO2, data = bw) # not significant 

# Clover
c <- BM %>% filter(Plant == "C")
t.test(Plant.biomass~CO2, data = c) # not significant 

# Dandelion
d <- BM %>% filter(Plant == "D")
t.test(Plant.biomass~CO2, data = d) # not significant 

# Lacy Phacelia
l <- BM %>% filter(Plant == "LP")
t.test(Plant.biomass~CO2, data = l) # not significant 

# Nasturtium
n <- BM %>% filter(Plant == "N")
t.test(Plant.biomass~CO2, data = n) # not significant 

# Partridge pea
p <- BM %>% filter(Plant == "PP")
t.test(Plant.biomass~CO2, data = p) # significant

# Sweet Alyssum
sa <- BM %>% filter(Plant == "SA")
t.test(Plant.biomass~CO2, data = sa) # not significant 

# Sunflower
sf <- BM %>% filter(Plant == "SF")
t.test(Plant.biomass~CO2, data = sf) # not significant 
