## Bernauer et al. (in prep for SciReports) Elevated atmospheric CO2 has 
  ## small, species-specific effects on pollen chemistry and plant growth across
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
library("see")
library("DHARMa")
library("GLMMadaptive")
library("emmeans")
library("multcomp")

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
  geom_smooth(method = "lm", aes(fill = CO2))+
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
m.1 <- lm(log(Height_cm) ~ CO2*Plant + Round + week + Chamber, data = height)
plot(simulationOutput <- simulateResiduals(fittedModel=m.1, plot =F))
anova(m.1) 

emm_model1 <- emmeans(m.1, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

# predict changes in height at week 10 for plants with significant height response to co2 treatment
# Make predictions using fixed effect only and then compare differences at week 10
pred.h <- height %>% 
  mutate(pred_fixef = predict(m.1, newdata = ., re.form = NA),
         pred_ranef = predict(m.1, newdata = ., re.form = NA))

pred.h$pred_fixef_cm <- exp(pred.h$pred_fixef)

pred.h.sum <- pred.h %>% group_by(Plant, CO2, week) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(pred_fixef_cm, na.rm = T),
    sd = sd(pred_fixef_cm, na.rm = T)
  )
pred.h.sum$se <- pred.h.sum$sd/sqrt(pred.h.sum$count)

##############
## leaves

# remove any plant, date, round, chambers w/o 3+ observations
leaves <- leaves %>% group_by(Plant, Chamber, Round, Date) %>% filter(n()>2) %>% ungroup()
leaves$CO2 <- as.factor(leaves$CO2)
plants.leaves <- c("Borage", "Buckwheat","Red Clover", "Dandelion", "Lacy Phacelia","Nasturtium","Partridge Pea", "Sweet Alyssum", "Sunflower")
names(plants.leaves) <- c("B", "BW", "C", "D", "LP", "N", "PP", "SA", "SF")

ggplot(leaves, aes(x=datenum, y = Leaf_no, color = CO2))+
  geom_point(aes(color = CO2), alpha = 0.4, size = 0.5)+
  geom_smooth(method = "lm", aes(fill = CO2))+
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
m.B1.log <- lm(log(Leaf_no) ~ CO2*Plant + Round + week + Chamber, data = leaves)
plot(simulationOutput <- simulateResiduals(fittedModel=m.B1.log, plot =F))
anova(m.B1.log) # CO2 not significant, but Plant species and Co2 x Plant species, so compare each species alone

emm_model1 <- emmeans(m.B1.log, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

## predict leaves at week 10
pred.l <- leaves %>% 
  mutate(pred_fixef = predict(m.B1.log, newdata = ., re.form = NA),
         pred_ranef = predict(m.B1.log, newdata = ., re.form = NA))

pred.l$pred_fixef_cm <- exp(pred.l$pred_fixef)

pred.l.sum <- pred.l %>% group_by(Plant, CO2, week) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(pred_fixef_cm, na.rm = T),
    sd = sd(pred_fixef_cm, na.rm = T)
  )
pred.l.sum$se <- pred.l.sum$sd/sqrt(pred.l.sum$count)

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
  geom_smooth(method = "lm", aes(fill = CO2))+
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

# round date to nearest week
flowers.1 <- flowers %>% filter(week %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17))

# full model w/all plant species
m.3 <- glmer.nb(Flower_no ~ CO2*Plant+Round+Chamber + (1|week),
                 data = flowers.1)
plot(simulationOutput <- simulateResiduals(fittedModel=m.3, plot =F))
car::Anova(m.3) # report this

emm_model1 <- emmeans(m.3, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

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
m.4.log <- lm(log(datenum) ~ CO2*Plant+Round + Chamber, data = ff)
plot(simulationOutput <- simulateResiduals(fittedModel = m.4.log, plot = F))
summary(m.4.log)
anova(m.4.log)

emm_model1 <- emmeans(m.4, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

 
######################
## flower diameter

# remove NAs'
W.area <- W.area %>% filter(!is.na(Plant))
# remove any plant, date, round, chambers w/o 3+ observations
area <- W.area %>% group_by(Plant, CO2) %>% filter(n()>2) %>% ungroup()
# convert to mean per plant
area.means <- area %>% group_by(PlantID, CO2, Round, Chamber, Plant) %>%
  summarise(mean = mean(Mean))

ggplot(area.means, aes(x=Plant, y = log(mean), fill = CO2))+
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
m.B2 <- lm(log(mean) ~ CO2*Plant + Round + Chamber, data = area.means)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B2, plot = F))
anova(m.B2) # no difference, but log-transformed has lower AIC
summary(m.B2)

emm_model1 <- emmeans(m.B2, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)

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
m.B <- lm(log(Plant.biomass) ~ CO2*Plant+Chamber, data = BM)
plot(simulationOutput <- simulateResiduals(fittedModel = m.B, plot = F))
summary(m.B)
anova(m.B)

emm_model1 <- emmeans(m.B, pairwise ~ CO2|Plant)
groups_emm_model1 <-cld(emm_model1, level = 0.05)
summary(groups_emm_model1)
pairs(emm_model1)