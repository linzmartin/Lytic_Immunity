#lytic immunity - ZOI analysis

#last updated: 06/30/23 LEM
################################################
# clear existing workspace
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")

#set wd to your project folder
getwd() #check working directory

#sessionInfo()
###########################

#####################
#load libraries needed:
library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(rstatix)
library(car) 
library(ggpubr)
#library(lme4)
#library(Rcpp)
#library(lmerTest)
#library(effects)
#library(broom)
#library(fitdistrplus)
library(emmeans)


#######################################################
#import the data and clean it up:


#import the data:
ZOI_data <- read_xlsx("ZOI_assay_sample_log_092723_LEM.xlsx")

str(ZOI_data)
head(ZOI_data)

#variables of interest:


#format the data:
ZOI_data$Age <- as.numeric(ZOI_data$Age)
ZOI_data$Temperature <- as.numeric(ZOI_data$Temperature)
ZOI_data$Treatment <- as.factor(ZOI_data$Treatment)
ZOI_data$Diameter_avg <- as.numeric(ZOI_data$Diameter_avg)
ZOI_data$ZOI_area <- as.numeric(ZOI_data$ZOI_area)
str(ZOI_data)

ZOI_data$Sample_ID <- as.factor(ZOI_data$Sample_ID)

levels(ZOI_data$Treatment)

#relevel treatment to make Naive come first:
ZOI_data <- ZOI_data %>%
  mutate(Treatment = fct_relevel(Treatment, 
                                 "Naïve","LB","E_coli","M_luteus","Water_control","NA"))

#omit NA treatments (no sample)
ZOI_data <- subset(ZOI_data, Treatment != "NA")

#subset data to exclude samples that did not have enough:
#remove NA values for temp and age
#omit NA treatments (no sample)
ZOI_data <- subset(ZOI_data, Enough!= 0)

#omit water group:
ZOI_wateronly <- subset(ZOI_data,Treatment == "Water_control")
ZOI_data <- subset(ZOI_data, Treatment != "Water_control")

str(ZOI_data)

##################

#center temp and aged data so that mean approximates zero intercept in model later:

#now center (take average and then subtract individ value):
ZOI_data$Temperature_centered <- scale(ZOI_data$Temperature, center=TRUE,scale=TRUE)
ZOI_data$Age_centered <- scale(ZOI_data$Age, center=TRUE,scale=TRUE)

str(ZOI_data)

#change character type to just numeric instead of scaled variable
ZOI_data$Temperature_centered <- as.numeric(ZOI_data$Temperature_centered)
ZOI_data$Age_centered <- as.numeric(ZOI_data$Age_centered)

str(ZOI_data)

#round digits for each:
ZOI_data$Temperature_centered <- round(ZOI_data$Temperature_centered, 2)
str(ZOI_data$Temperature_centered)
ZOI_data$Age_centered <-round(ZOI_data$Age_centered, 2)
str(ZOI_data$Age_centered)

###
#how centering works
unique(ZOI_data$Temperature_centered)
#check:
mean(ZOI_data$Temperature_centered) #approximates zero
mean(ZOI_data$Temperature)
(sum(27+30+32))/3

unique(ZOI_data$Age_centered)
mean(ZOI_data$Age_centered) #approximates zero
mean(ZOI_data$Age)
(sum(1+5+10+15))/4

################################
#create labels:
templabs <- c("27˚C","30˚C","32˚C")
#names(templabs)<- c("-2.52","2.48","0.48") #c("27","30","32") 
names(templabs)<- c("-1.15","1.13","0.22") #c("27","30","32") 


templabs_original <- c("27˚C","30˚C","32˚C")
names(templabs_original)<- c("27","30","32") 

agelabs <- c("1 day","5 days", "10 days", "15 days","Control")
#names(agelabs)<-c("-6","-2","3","8","NA") #c("1","5","10","15")
names(agelabs)<-c("-1.16","-0.39","0.58","1.55","NA")

agelabs_original <- c("1 day","5 days", "10 days", "15 days","Control")
names(agelabs_original)<-c("1","5","10","15","NA")

treatlabs <-c("Naïve","LB","E. coli","M. luteus","Control")
names(treatlabs)<-c("Naïve","LB","E_coli","M_luteus","Water_control")

#relevel treatment to make Naive come first:

#ZOI_data <- ZOI_data %>%
 # mutate(Treatment = fct_relevel(Treatment, 
  #                               "Naïve","LB","E_coli","M_luteus","Water_control","NA"))
#ZOI_data$Treatment <- relevel(ZOI_data$Treatment, ref="Naïve")
levels(ZOI_data$Treatment)
levels(ZOI_data$Age_centered)

#ZOI_data <- ZOI_data %>%
# mutate(Age = fct_relevel(Age, "1","5","10","15","NA"))



#summary stats:

ZOI_summary <- ZOI_data %>%
  group_by(Treatment,Age,Temperature,Sample_ID) %>%
  reframe(#Temperature = Temperature, 
            #Age = Age,
            #Treatment=Treatment,
            #Sample_ID = Sample_ID,
            Technical_Rep = Technical_Rep, 
            #E_coli_infection = E_coli_infection,
            #M_luteus_infection = M_luteus_infection,
            mean_diameter = mean(Diameter_avg),
            mean_area = mean(ZOI_area),
            sd_diameter = sd(Diameter_avg),
            n_diameter = n(),
            SE_diameter = sd(Diameter_avg)/sqrt(n()),
            SE_area = sd(ZOI_area)/sqrt(n()))

ZOI_summary <- as.data.frame(ZOI_summary)

#check output and format:
str(ZOI_summary)

levels(ZOI_summary$Treatment) 

#then count number of replicates

count_replicates <- ZOI_summary %>%
  group_by(Temperature,Treatment,Age,Technical_Rep)%>%
  summarise(n=n())

count_replicates<- as.data.frame(count_replicates)
head(count_replicates)

write_xlsx(count_replicates, "number_reps_per_treat_ZOI_assays.xlsx")

count_replicates_biological <- subset(count_replicates, Technical_Rep == 1)

require(EnvStats)
#plot number of replicates per group:
count_replicates_biological %>%
  ggplot(aes(x=Treatment,y=n,color=Treatment))+
  geom_bar(aes(fill=Treatment),color="black",
           stat="identity",
           position=position_dodge())+
  facet_grid(Age~Temperature,labeller = labeller(Temperature=templabs_original,Age=agelabs_original,Treatment=treatlabs))+
  #labs(title="Number of Replicates per Group",
  #caption="Error bars represent mean +/- standard error")+
  scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  scale_y_continuous(breaks = c(0,2,4,6))+
  ylab(expression("Number of Replicates"))+ theme_bw()+
  theme(legend.position = "none")+ #stat_n_text(size=3,fontface = "bold",y.pos=1.5,y.expand.factor = 0.5)+
  theme(axis.text.x = element_text(angle=45, hjust=1))


#now take average of technical reps:

ZOI_summary_bioreps <- ZOI_summary %>%
  group_by(Treatment,Age,Temperature) %>%
  reframe(#Temperature = Temperature, 
           # Age = Age,
            #Treatment=Treatment,
            #Sample_ID = Sample_ID,
            #Technical_Rep = Technical_Rep, 
            #E_coli_infection = E_coli_infection,
            #M_luteus_infection = M_luteus_infection,
            mean_diameter_all = mean(mean_diameter),
            mean_area_all = mean(mean_area),
            sd_diameter_all = sd(mean_diameter),
            n_diameter_all = n(),
            SE_diameter_all = sd(mean_diameter)/sqrt(n()),
            SE_area_all = sd(mean_area)/sqrt(n()))

ZOI_summary_bioreps <- as.data.frame(ZOI_summary_bioreps)

###
ZOI_summary <- ZOI_data %>%
  group_by(Treatment,Age_centered,Temperature_centered,Sample_ID) %>%
  reframe(#Temperature = Temperature, 
    #Age = Age,
    #Treatment=Treatment,
    #Sample_ID = Sample_ID,
    #Technical_Rep = Technical_Rep, 
    #E_coli_infection = E_coli_infection,
    #M_luteus_infection = M_luteus_infection,
    mean_diameter = mean(Diameter_avg),
    mean_area = mean(ZOI_area),
    sd_diameter = sd(Diameter_avg),
    n_diameter = n(),
    SE_diameter = sd(Diameter_avg)/sqrt(n()),
    SE_area = sd(ZOI_area)/sqrt(n()))

ZOI_summary <- as.data.frame(ZOI_summary)
####
head(ZOI_data)
ZOI_data %>% ggplot(aes(x=Temperature_centered,y=ZOI_area))+
  geom_point(color=ZOI_data$Technical_Rep)+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age_centered, labeller = labeller(Age_centered=agelabs, Temperature_centered=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Area of the Zone of Inhibition (mm2)"))+
  theme_bw()

ZOI_data %>% ggplot(aes(x=Temperature_centered,y=ZOI_area,group=Temperature_centered))+
  geom_point(color=ZOI_data$Technical_Rep)+ geom_boxplot()+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age_centered, labeller = labeller(Age_centered=agelabs, Temperature_centered=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Area of the Zone of Inhibition (mm2)"))+
  theme_bw()

ZOI_data %>% ggplot(aes(x=Temperature_centered,y=Diameter_avg,group=Temperature_centered))+
  geom_point(color=ZOI_data$Technical_Rep)+ geom_boxplot()+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age_centered, labeller = labeller(Age_centered=agelabs, Temperature_centered=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Diameter of Zone of Inhibition (mm)"))+
  theme_bw()


#########################################################################

#check data for normality:

normtestresults <- ZOI_data %>%
  # mutate(mean_area = Delta_OD) %>%
  group_by(Temperature_centered, Treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ZOI_area)))
head(normtestresults)

library(broom)
normtest.glance <- normtestresults %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
normtest.glance #values are not normal! p < 0.05.

boxplot(ZOI_area ~ Temperature_centered+Age_centered+Treatment,
        col=c("white","lightgray"),
        data = ZOI_data)

ZOI_data %>%
  ggplot(aes(x=Temperature_centered,y=ZOI_area,color=Age,group=Temperature_centered))+
  facet_grid(Treatment~Age_centered, labeller = labeller(Age_centered=agelabs, Temperature_centered=templabs))+
  labs(x="Temperature (˚C)")+
  #ylab(expression(~italic("E. coli")~" Norm Mel Area"))+
  geom_boxplot(notch = FALSE)+
  theme_bw()+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme(legend.position = "none")

library(rcompanion)
plotNormalHistogram(ZOI_data$Diameter_avg)
ggplot(ZOI_data, aes((Diameter_avg))) + geom_histogram()

plotNormalHistogram(ZOI_data$ZOI_area)
ggplot(ZOI_data, aes((ZOI_area))) + geom_histogram()


#log-transform values

#ZOI_data$log_area <- log(ZOI_data$ZOI_area)

ggplot(ZOI_data, aes(log(ZOI_area))) + geom_histogram()
ggplot(ZOI_data, aes(log(ZOI_area+0.1))) + geom_histogram()
ggplot(ZOI_data, aes(log(ZOI_area+1))) + geom_histogram()


#model + anova
library(lme4)
lm_full <- lmer(log(ZOI_area) ~ Temperature+Age+Treatment + 
                  poly(Temperature,2)+poly(Age,2)+
                  Temperature*Age+
                  Temperature*Treatment+
                  Age*Treatment+
                  Temperature*Age*Treatment+
                (1|Sample_ID)+
                (1|Plate_ID_Total)+
                (1|Batch_Number),
            data=ZOI_data,
            REML=FALSE)

summary(lm_full) #can possibly drop a random effect
Anova(lm_full,type=2,test.statistic = "Chisq") #can't seem to drop any main or interactive effects
plot(lm_full) #mostly random

drop1(lm_full,test="Chisq") #signif - cannot remove 3way or any below

library(lmerTest)
rand(lm_full)#can get rid of sample ID and plate ID

##

lm_mod1 <- lmer(log(ZOI_area) ~ Temperature+Age+Treatment + 
                   Temperature*Age+
                   Temperature*Treatment+
                   Age*Treatment+
                  # Temperature*Age*Treatment+
                  # (1|Sample_ID)+
                  (1|Plate_ID_Total)+
                  (1|Batch_Number),
                 data=ZOI_data,
                 REML=FALSE)

summary(lm_mod1)
plot(lm_mod1)
Anova(lm_mod1,type=2,test.statistic = "Chisq") 

drop1(lm_mod1)
rand(lm_mod1)

lm_mod2 <- lmer(log(ZOI_area) ~ Temperature+Age+Treatment + 
                  #Temperature*Age+
                  Temperature*Treatment+
                  #Age*Treatment+
                  # Temperature*Age*Treatment+
                  # (1|Sample_ID)+
                  (1|Plate_ID_Total)+
                  (1|Batch_Number),
                data=ZOI_data,
                REML=FALSE)

summary(lm_mod2)
plot(lm_mod2)
Anova(lm_mod2,type=2,test.statistic = "Chisq") 

#interaction of temperature and treatment seems strong- what if we separate treatments to get at the effects of temp and age on each treatment?

##################
#confidence intervals for model estimates:
confint(lm_mod1)
exp(confint(lm_mod1))

library("merDeriv")
library("parameters")
modelvalues <- parameters::model_parameters(lm_mod1, exponentiate = TRUE)
modelvalues

modelvalues<- as.data.frame(modelvalues)

#####################
#what does continuous data look like?
str(ZOI_data)
#ZOIcontinuous <- ZOI_data

#ZOIcontinuous$Temperature <- as.numeric(ZOIcontinuous$Temperature)
#ZOIcontinuous$Age <- as.numeric(ZOIcontinuous$Age)
#str(ZOIcontinuous)

#ggplot(ZOIcontinuous, aes(log(ZOI_area))) + geom_histogram()

lm_full_continuous <- lmer(log(ZOI_area) ~ Temperature+Age+Treatment + 
                  Temperature*Age+
                  Temperature*Treatment+
                  Age*Treatment+
                  Temperature*Age*Treatment+
                  (1|Sample_ID)+
                   (1|Plate_ID_Total)+
                   (1|Batch_Number),
                data=ZOI_data,
                REML=FALSE)

plot(lm_full_continuous)

summary(lm_full_continuous)

Anova(lm_full_continuous,type=2,test.statistic = "Chisq") 

#basically the same thing

lm_full_continuous2 <- lmer(log(ZOI_area) ~ poly(Temperature,2)+
                              poly(Age,2)+
                              Treatment +
                             Temperature*Age+
                             Temperature*Treatment+
                             Age*Treatment+
                             #Temperature*Age*Treatment+
                             # (1|Sample_ID)+
                             (1|Plate_ID_Total)+
                             (1|Batch_Number),
                           data=ZOI_data,
                           REML=FALSE)
plot(lm_full_continuous2)

summary(lm_full_continuous2)

Anova(lm_full_continuous2,type=2,test.statistic = "Chisq") 
drop1(lm_full_continuous2)
rand(lm_full_continuous2)

lm_full_continuous3 <- lmer(log(ZOI_area) ~ poly(Temperature,2)+
                              poly(Age,2)+
                              Treatment +
                              Temperature*Age+
                              Temperature*Treatment+
                              #Age*Treatment+
                              #Temperature*Age*Treatment+
                              # (1|Sample_ID)+
                              (1|Plate_ID_Total)+
                              (1|Batch_Number),
                            data=ZOI_data,
                            REML=FALSE)
plot(lm_full_continuous3)

summary(lm_full_continuous3)


Anova(lm_full_continuous3,type=2,test.statistic = "Chisq") 



######################
#treatment has a major effect.

#stratify categorical by treatment:
#ZOI_data <- subset(ZOI_data,Technical_Rep == 1)

ZOI_data$Temperature_centered <- ZOI_data$Temperature
ZOI_data$Age_centered <- ZOI_data$Age

Naive_ZOI <- subset(ZOI_data,Treatment == "Naïve")
LB_ZOI <- subset(ZOI_data,Treatment == "LB")
Ecoli_ZOI <- subset(ZOI_data,Treatment == "E_coli")
Mlut_ZOI <- subset(ZOI_data,Treatment == "M_luteus")
Naive_and_LB <- subset(ZOI_data,Treatment == "LB" |Treatment == "Naïve")
Naive_and_Ecoli <- subset(ZOI_data,Treatment == "Naïve"|Treatment == "E_coli")
Naive_and_Mluteus <- subset(ZOI_data,Treatment == "Naïve"|Treatment == "M_luteus")


#try analysis by treatment:

lmer_Naive_1 <- lmer(log(ZOI_area) ~ Temperature_centered+
                       Age_centered + 
                       Temperature_centered*Age_centered+
                        (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Batch_Number),
                     data=Naive_ZOI,
                     REML=FALSE)

summary(lmer_Naive_1)

rand(lmer_Naive_1)
drop1(lmer_Naive_1)

lm_Naive_2 <- lmer(log(ZOI_area) ~ Temperature_centered+
                     Age_centered + 
                     Temperature_centered*Age_centered+
                     #(1|Sample_ID)+
                     #(1|Plate_ID_Total)+
                     (1|Batch_Number),
                   data=Naive_ZOI,
                   REML=FALSE)
summary(lm_Naive_2)
plot(lm_Naive_2)
Anova(lm_Naive_2,type=2)

drop1(lm_Naive_2)

lm_Naive_3 <-  lmer(log(ZOI_area) ~ Temperature_centered+
                      Age_centered + 
                      #Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      #(1|Plate_ID_Total)+
                      (1|Batch_Number),
                    data=Naive_ZOI,
                    REML=FALSE)
anova(lm_Naive_2,lm_Naive_3)
AIC(lm_Naive_2,lm_Naive_3)

summary(lm_Naive_3)

plot(lm_Naive_3)
Anova(lm_Naive_3,type=2)

#save outputs:
sink("Naive_summary.txt")
print(summary(lm_Naive_3))
sink()  # returns output to the console

sink("Naive_ANOVAoutput.txt")
print(Anova(lm_Naive_3,type="2",test.statistic = "Chisq"))
sink()


naive_anova <- parameters::model_parameters(Anova(lm_Naive_3),type=2)
naive_anova<-as.data.frame(naive_anova)
#naive_anova <- as.data.frame(Anova(lm_Naive_3,type="2",test.statistic = "Chisq"))
head(naive_anova)
write_xlsx(naive_anova, "Naive_ANOVA.xlsx")

modelvalues_naive <- parameters::model_parameters(lm_Naive_3, exponentiate = TRUE)
head(modelvalues_naive)

write_xlsx(modelvalues_naive, "Naive_model_expcoefandCIs.xlsx")

library(effectsize)

Naive_effectsize_partial <- eta_squared(lm_Naive_3,partial=TRUE)
Naive_effectsize_partial
Naive_effectsize_partial_dataframe <- as.data.frame(Naive_effectsize_partial)

write_xlsx(Naive_effectsize_partial_dataframe,"Naive_partialeffectsizes.xlsx")

#takeaway: temp affects naive, but not age or interaction

###############################
lmer_LB_1 <- lmer(log(ZOI_area) ~ Temperature_centered +
                    Age_centered + 
                    Temperature_centered*Age_centered+
                    (1|Sample_ID)+
                    (1|Plate_ID_Total)+
                    (1|Seeded_plates_OD),
                     data=LB_ZOI,
                     REML=FALSE)

summary(lmer_LB_1)
rand(lmer_LB_1) #can get rid of all random effects except seeded plate OD
drop1(lmer_LB_1) #remove two way int.

lm_LB_2 <- lmer(log(ZOI_area) ~ Temperature_centered +
                  Age_centered+
                 # Temperature_centered*Age_centered+
                  #(1|Sample_ID)+
                  #(1|Plate_ID_Total)+
                  (1|Seeded_plates_OD),
                data=LB_ZOI)
summary(lm_LB_2)
plot(lm_LB_2)
rand(lm_LB_2)
Anova(lm_LB_2,type=2)
drop1(lm_LB_2)

#lm_LB_3 <- lm(log(ZOI_area) ~ Temperature_centered +
 #               Age_centered,
  #            data=LB_ZOI)

#summary(lm_LB_3)
#plot(lm_LB_3)

#anova(lm_LB_2,lm_LB_3)
#AIC(lm_LB_2,lm_LB_3)

#no effects on injury group


lm_LB_4 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                poly(Age_centered,2)+
                  (1|Seeded_plates_OD),
              data=LB_ZOI)

plot(lm_LB_4)
Anova(lm_LB_2,type=2)
AIC(lm_LB_2,lm_LB_4) #2 is better

#try glm:
lm_LB5 <- glm(ZOI_area ~ Temperature_centered + Age_centered,
              data=LB_ZOI, family=gaussian(link="log"))

summary(lm_LB5)
plot(lm_LB5) #not great
AIC(lm_LB_2,lm_LB5) 
logLik(lm_LB5)
logLik(lm_LB_2)
#2 is best

##
#save outputs:
sink("LB_summary.txt")
print(summary(lm_LB_2))
sink()  # returns output to the console

sink("LB_ANOVAoutput.txt")
print(Anova(lm_LB_2,type="2",test.statistic = "Chisq"))
sink()

LB_anova <- parameters::model_parameters(Anova(lm_LB_2),type=2)
LB_anova<-as.data.frame(LB_anova)

write_xlsx(LB_anova,"LB_ANOVA.xlsx")

modelvalues_LB <- parameters::model_parameters(lm_LB_2, exponentiate = TRUE)
head(modelvalues_LB)

write_xlsx(modelvalues_LB, "LB_model_expcoefandCIs.xlsx")

LB_effectsize_partial <- eta_squared(lm_LB_2)
LB_effectsize_partial_dataframe <- as.data.frame(LB_effectsize_partial)
LB_effectsize_partial_dataframe
write_xlsx(LB_effectsize_partial_dataframe,"LB_partialeffectsizes.xlsx")



###########################


#infection:
lmer_Ecoli_1 <- lmer(log(ZOI_area) ~  Temperature_centered +
                       Age_centered + 
                       Temperature_centered*Age_centered+
                       (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Seeded_plates_OD),
                     data=Ecoli_ZOI,
                     REML=FALSE)

summary(lmer_Ecoli_1)

rand(lmer_Ecoli_1) #get rid of plate ID

lmer_Ecoli_2 <- lm(log(ZOI_area) ~  Temperature_centered +
                       Age_centered + 
                       Temperature_centered*Age_centered,
                   data=Ecoli_ZOI)

lmer_Ecoli_2 <- lmer(log(ZOI_area) ~  Temperature_centered +
                     Age_centered + 
                     Temperature_centered*Age_centered+
                       (1|Sample_ID)+
                       #(1|Plate_ID),
                       (1|Batch_Number),
                     data=Ecoli_ZOI)

summary(lmer_Ecoli_2)
Anova(lmer_Ecoli_2,type=2)
plot(lmer_Ecoli_2)
rand(lmer_Ecoli_2)

drop1(lmer_Ecoli_2)



AIC(lmer_Ecoli_2)

lmer_Ecoli_3 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       Age_centered + 
                       Temperature_centered*Age_centered+
                       (1|Sample_ID)+
                     #(1|Plate_ID),
                     (1|Batch_Number),
                     data=Ecoli_ZOI)
summary(lmer_Ecoli_3)
Anova(lmer_Ecoli_3)
plot(lmer_Ecoli_3)


lmer_Ecoli_4 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       poly(Age_centered,2) + 
                       Temperature_centered*Age_centered+
                       (1|Sample_ID) + (1|Batch_Number),
                     #(1|Plate_ID)
                     data=Ecoli_ZOI)
#temp has impact on E. coli zoi but not age or int of temp and age

summary(lmer_Ecoli_4)
AIC(lmer_Ecoli_4)
plot(lmer_Ecoli_4)

lmer_Ecoli_5 <- lmer(log(ZOI_area) ~  Temperature_centered +
                       poly(Age_centered,2) + 
                       Temperature_centered*Age_centered+
                       (1|Sample_ID) + (1|Batch_Number),
                     #(1|Plate_ID),
                     #(1|Batch_Number),
                     data=Ecoli_ZOI)

summary(lmer_Ecoli_5)
plot(lmer_Ecoli_5)
Anova(lmer_Ecoli_5,type=2)

AIC(lmer_Ecoli_5)


lmer_Ecoli_6 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                     poly(Age_centered,2)+
                   poly(Temperature_centered*Age_centered,2)+
                     #(1|Sample_ID) + 
                     (1|Batch_Number),
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)

summary(lmer_Ecoli_6)
plot(lmer_Ecoli_6)
rand(lmer_Ecoli_6)
Anova(lmer_Ecoli_6,type=2)


AIC(lmer_Ecoli_3,lmer_Ecoli_4,lmer_Ecoli_5,lmer_Ecoli_6)
#go with 6

lmer_Ecoli_8 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                   poly(Age_centered,3)+
                     (1|Sample_ID) + (1|Batch_Number),
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)
plot(lmer_Ecoli_8)

AIC(lmer_Ecoli_3,lmer_Ecoli_4,lmer_Ecoli_5,lmer_Ecoli_6,lmer_Ecoli_8)

library(modelr)
grid <- Ecoli_ZOI %>% 
  data_grid(Temperature_centered, .model = lmer_Ecoli_6) %>% 
  add_predictions(lmer_Ecoli_6)

ggplot(grid, aes(Temperature_centered, pred)) + 
  geom_point()
####################################################################

ggplot(Ecoli_ZOI, aes(Temperature_centered, ZOI_area)) + 
  geom_point() + geom_boxplot(aes(group=Temperature_centered))
ggplot(Ecoli_ZOI, aes(Temperature_centered, log(ZOI_area))) + 
  geom_point() + geom_boxplot(aes(group=Temperature_centered))
ggplot(Ecoli_ZOI, aes(Age_centered, ZOI_area)) + 
  geom_point()+ geom_boxplot(aes(group=Age_centered))
ggplot(Ecoli_ZOI, aes(Age_centered, log(ZOI_area))) + 
  geom_point() + geom_boxplot(aes(group=Age_centered))


grid <- Ecoli_ZOI %>% 
  data_grid(Temperature_centered, .model = lmer_Ecoli_6) %>% 
  add_predictions(lmer_Ecoli_6)

ggplot(Ecoli_ZOI, aes(x=Temperature_centered, 
                      y=log(ZOI_area),group=Temperature_centered)) + 
 geom_boxplot()+ 
  geom_point(aes(x=Temperature_centered,y=pred),data=grid,colour="red",size=2)


grid_age <- Ecoli_ZOI %>% 
  data_grid(Age_centered, .model = lmer_Ecoli_6) %>% 
  add_predictions(lmer_Ecoli_3)

ggplot(Ecoli_ZOI, aes(x=Age_centered, 
                      y=log(ZOI_area),group=Age_centered)) + 
  geom_boxplot()+ 
  geom_point(aes(x=Age_centered,y=pred),data=grid_age,colour="red",size=2)

#####
Ecoli_ZOI <- Ecoli_ZOI %>%
  add_residuals(lmer_Ecoli_3)

Ecoli_ZOI %>%
  ggplot(aes(Temperature_centered, resid,color=as.factor(Age_centered)))+
  geom_ref_line(h=0)+geom_point()
Ecoli_ZOI %>%
  ggplot(aes(Age_centered, resid, color=as.factor(Temperature_centered)))+
  geom_ref_line(h=0)+geom_point()

grid_tandage <- ZOI_data %>% 
  data_grid(Age_centered, Temperature_centered) %>% 
  add_predictions(lmer_Ecoli_3, "log_ZOI_area")



ggplot(Ecoli_ZOI, aes(x=Age_centered, 
                      y=log(ZOI_area),group=Age_centered)) + 
  geom_boxplot()+
  geom_point(aes(x=Age_centered,y=log_ZOI_area),data=grid_tandage,colour="red",size=2)+
  facet_wrap(~Temperature_centered)

###
library(MASS)
mod3 <- MASS::rlm(ZOI_area ~ poly(Temperature_centered,2)*Age_centered,
                 data=Ecoli_ZOI)

ZOI_data%>% 
  add_residuals(mod3, "resid") %>% 
  ggplot(aes(Temperature_centered, resid, color=as.factor(Temperature_centered))) + 
  geom_hline(yintercept = 0, size = 2, colour = "white") +
  facet_wrap(~Age_centered)+
  geom_point()



#
lmer_Mlut_1 <- lmer(log(ZOI_area) ~ Temperature+
                       Age + 
                       Temperature*Age+
                       (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Batch_Number),
                     data=Mlut_ZOI,
                     REML=FALSE)

summary(lmer_Mlut_1)

rand(lmer_Mlut_1) #get rid of sample ID and batch number

lmer_Mlut_2 <- lmer(log(ZOI_area) ~ Temperature+
                      Age + 
                      Temperature*Age+
                      #(1|Sample_ID)+
                      (1|Plate_ID_Total),
                      #(1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=FALSE)

summary(lmer_Mlut_2)

plot(lmer_Mlut_2)

Anova(lmer_Mlut_2)
drop1(lmer_Mlut_2)

lmer_Mlut_3 <- lmer(log(ZOI_area) ~ Temperature+
                                     Age + poly(Age,2)+
                                     #Temperature*Age+
                                     #(1|Sample_ID)+
                                     (1|Plate_ID_Total),
                                   #(1|Batch_Number),
                                   data=Mlut_ZOI,
                                   REML=FALSE)

summary(lmer_Mlut_3)

plot(lmer_Mlut_3)

Anova(lmer_Mlut_3)

lmer_Mlut_4 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                      poly(Age_centered,2)+
                      Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      (1|Plate_ID_Total),
                    #(1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=FALSE)

summary(lmer_Mlut_4)
Anova(lmer_Mlut_4)
drop1(lmer_Mlut_4)
#temperature and age both significantly shape m.lut response but do not interact together
########################
#Naive vs. others:

#infection:
lmer_NaivevEcoli_1 <- lmer(log(ZOI_area) ~  Temperature_centered +
                       Age_centered + Treatment+
                       Temperature_centered*Age_centered +
                     Temperature_centered*Treatment+
                     Age_centered*Treatment+
                       Temperature_centered*Age_centered*Treatment+
                       (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Seeded_plates_OD),
                     data=Naive_and_Ecoli,
                     REML=FALSE)

plot(lmer_NaivevEcoli_1)
summary(lmer_NaivevEcoli_1)
drop1(lmer_NaivevEcoli_1)
rand(lmer_NaivevEcoli_1)

lmer_NaivevEcoli_2 <- lmer(log(ZOI_area) ~  Temperature_centered +
                             Age_centered + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             Age_centered*Treatment+
                             #Temperature_centered*Age_centered*Treatment+
                             (1|Sample_ID)+
                             #(1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Ecoli,
                           REML=FALSE)

plot(lmer_NaivevEcoli_2)
summary(lmer_NaivevEcoli_2)
drop1(lmer_NaivevEcoli_2)
rand(lmer_NaivevEcoli_2)

lmer_NaivevEcoli_3 <- lmer(log(ZOI_area) ~  Temperature_centered +
                             Age_centered + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             #Age_centered*Treatment+
                             #Temperature_centered*Age_centered*Treatment+
                            # (1|Sample_ID)+
                             #(1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Ecoli,
                           REML=FALSE)

plot(lmer_NaivevEcoli_3)
summary(lmer_NaivevEcoli_3)
drop1(lmer_NaivevEcoli_3)
rand(lmer_NaivevEcoli_3)

Anova(lmer_NaivevEcoli_3)

lmer_NaivevEcoli_4 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                             Age_centered + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             #Age_centered*Treatment+
                             #Temperature_centered*Age_centered*Treatment+
                             #(1|Sample_ID)+
                             #(1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Ecoli,
                           REML=FALSE)


plot(lmer_NaivevEcoli_4)
summary(lmer_NaivevEcoli_4)

qqnorm(resid(lmer_NaivevEcoli_4))
qqline(resid(lmer_NaivevEcoli_4))
drop1(lmer_NaivevEcoli_4)
rand(lmer_NaivevEcoli_4)

lmer_NaivevEcoli_5 <- lmer(log(ZOI_area) ~  poly(Temperature,2) +
                             poly(Age,2) + Treatment+
                             Temperature*Age +
                             Temperature*Treatment+
                             #Age_centered*Treatment+
                             #Temperature_centered*Age_centered*Treatment+
                             (1|Sample_ID)+
                             #(1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Ecoli,
                           REML=FALSE)
plot(lmer_NaivevEcoli_5)
summary(lmer_NaivevEcoli_5)
rand(lmer_NaivevEcoli_5)

lmer_NaivevEcoli_5 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                             poly(Age_centered,2) + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             #Age_centered*Treatment+
                             #Temperature_centered*Age_centered*Treatment+
                             (1|Sample_ID)+
                             #(1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Ecoli,
                           REML=FALSE)
plot(lmer_NaivevEcoli_5)
summary(lmer_NaivevEcoli_5)

lmer_NaivevEcoli_6 <- lmer(log(ZOI_area) ~  Temperature_centered +
                             poly(Age_centered,2) + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             #Age_centered*Treatment+
                             #Temperature_centered*Age_centered*Treatment+
                             (1|Sample_ID)+
                             #(1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Ecoli,
                           REML=FALSE)


AIC(lmer_NaivevEcoli_1,lmer_NaivevEcoli_2,lmer_NaivevEcoli_3,lmer_NaivevEcoli_4,lmer_NaivevEcoli_5,lmer_NaivevEcoli_6)

#go with lmer_NaivevEcoli_5 to have low AIC but also include parabolic effects of age and temp, include sample

Anova(lmer_NaivevEcoli_5, type=2)
#####

#Mlut infection:
lmer_NaivevMlut_1 <- lmer(log(ZOI_area) ~  Temperature_centered +
                             Age_centered + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             Age_centered*Treatment+
                             Temperature_centered*Age_centered*Treatment+
                             (1|Sample_ID)+
                             (1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_Mluteus,
                           REML=FALSE)

summary(lmer_NaivevMlut_1)
drop1(lmer_NaivevMlut_1)
rand(lmer_NaivevMlut_1)
plot(lmer_NaivevMlut_1)

lmer_NaivevMlut_2 <- lmer(log(ZOI_area) ~  Temperature_centered +
                            Age_centered + Treatment+
                            Temperature_centered*Age_centered +
                            Temperature_centered*Treatment+
                            Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            (1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_2)
drop1(lmer_NaivevMlut_2)
rand(lmer_NaivevMlut_2)
plot(lmer_NaivevMlut_2)

lmer_NaivevMlut_3 <- lmer(log(ZOI_area) ~  Temperature_centered +
                            Age_centered + Treatment+
                            Temperature_centered*Age_centered +
                            Temperature_centered*Treatment+
                            #Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            #(1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_3)
drop1(lmer_NaivevMlut_3)
rand(lmer_NaivevMlut_3)
plot(lmer_NaivevMlut_3)

lmer_NaivevMlut_4 <- lmer(log(ZOI_area) ~  Temperature_centered +
                            Age_centered + Treatment+
                            Temperature_centered*Age_centered +
                            Temperature_centered*Treatment+
                            #Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            #(1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_4)
drop1(lmer_NaivevMlut_4)
rand(lmer_NaivevMlut_4)
plot(lmer_NaivevMlut_4)

AIC(lmer_NaivevMlut_3,lmer_NaivevMlut_4)

lmer_NaivevMlut_5 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                            Age_centered + Treatment+
                            Temperature_centered*Age_centered +
                            Temperature_centered*Treatment+
                            #Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            #(1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_5)
drop1(lmer_NaivevMlut_5)
rand(lmer_NaivevMlut_5)
plot(lmer_NaivevMlut_5)
Anova(lmer_NaivevMlut_5,type=2)

lmer_NaivevMlut_6 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                            poly(Age_centered,2) + Treatment+
                            Temperature_centered*Age_centered +
                            Temperature_centered*Treatment+
                            #Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            #(1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_6)
drop1(lmer_NaivevMlut_6)
rand(lmer_NaivevMlut_6)
plot(lmer_NaivevMlut_6)
Anova(lmer_NaivevMlut_6)
qqnorm(resid(lmer_NaivevMlut_6))
qqline(resid(lmer_NaivevMlut_6))
AIC(lmer_NaivevMlut_5,lmer_NaivevMlut_6)

lmer_NaivevMlut_7 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                            poly(Age_centered,3) + Treatment+
                            Temperature_centered*Age_centered +
                            Temperature_centered*Treatment+
                            #Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            #(1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_7)
drop1(lmer_NaivevMlut_7)
rand(lmer_NaivevMlut_7)
plot(lmer_NaivevMlut_7)
Anova(lmer_NaivevMlut_7)

lmer_NaivevMlut_8 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                            poly(Age_centered,2) + Treatment+
                            poly(Temperature_centered,2)*poly(Age_centered,2)+
                            Temperature_centered*Treatment+
                            #Age_centered*Treatment+
                            #Temperature_centered*Age_centered*Treatment+
                            #(1|Sample_ID)+
                            #(1|Plate_ID_Total)+
                            (1|Seeded_plates_OD),
                          data=Naive_and_Mluteus,
                          REML=FALSE)

summary(lmer_NaivevMlut_8)
drop1(lmer_NaivevMlut_8)
rand(lmer_NaivevMlut_8)
plot(lmer_NaivevMlut_8)
Anova(lmer_NaivevMlut_8)
qqnorm(resid(lmer_NaivevMlut_8))
qqline(resid(lmer_NaivevMlut_8))

AIC(lmer_NaivevMlut_4,lmer_NaivevMlut_5,lmer_NaivevMlut_6,lmer_NaivevMlut_7,lmer_NaivevMlut_8)
BIC(lmer_NaivevMlut_4,lmer_NaivevMlut_5,lmer_NaivevMlut_6,lmer_NaivevMlut_7)

#go with lmer_NaivevMlut_6
###
#Naive vs. Injury:
lmer_NaivevLB_1 <- lmer(log(ZOI_area) ~  Temperature_centered +
                             Age_centered + Treatment+
                             Temperature_centered*Age_centered +
                             Temperature_centered*Treatment+
                             Age_centered*Treatment+
                             Temperature_centered*Age_centered*Treatment+
                             (1|Sample_ID)+
                             (1|Plate_ID_Total)+
                             (1|Seeded_plates_OD),
                           data=Naive_and_LB,
                           REML=FALSE)

plot(lmer_NaivevLB_1)
summary(lmer_NaivevLB_1)
drop1(lmer_NaivevLB_1)
rand(lmer_NaivevLB_1)

lmer_NaivevLB_2 <- lmer(log(ZOI_area) ~  Temperature_centered +
                          Age_centered + Treatment+
                          Temperature_centered*Age_centered +
                          Temperature_centered*Treatment+
                          Age_centered*Treatment+
                          Temperature_centered*Age_centered*Treatment+
                          #(1|Sample_ID)+
                          #(1|Plate_ID_Total)+
                          (1|Seeded_plates_OD),
                        data=Naive_and_LB,
                        REML=FALSE)


plot(lmer_NaivevLB_2)
summary(lmer_NaivevLB_2)
drop1(lmer_NaivevLB_2)
rand(lmer_NaivevLB_2)

lmer_NaivevLB_3 <- lmer(log(ZOI_area) ~  Temperature_centered +
                          Age_centered + Treatment+
                          Temperature_centered*Age_centered +
                          Temperature_centered*Treatment+
                          Age_centered*Treatment+
                          #Temperature_centered*Age_centered*Treatment+
                          #(1|Sample_ID)+
                          #(1|Plate_ID_Total)+
                          (1|Seeded_plates_OD),
                        data=Naive_and_LB,
                        REML=FALSE)

plot(lmer_NaivevLB_3)
summary(lmer_NaivevLB_3)
drop1(lmer_NaivevLB_3)
rand(lmer_NaivevLB_3)

lmer_NaivevLB_4 <- lmer(log(ZOI_area) ~  Temperature_centered +
                          Age_centered + Treatment+
                          Temperature_centered*Age_centered +
                          Temperature_centered*Treatment+
                          #Age_centered*Treatment+
                          #Temperature_centered*Age_centered*Treatment+
                          #(1|Sample_ID)+
                          #(1|Plate_ID_Total)+
                          (1|Seeded_plates_OD),
                        data=Naive_and_LB,
                        REML=FALSE)

plot(lmer_NaivevLB_4)
summary(lmer_NaivevLB_4)
drop1(lmer_NaivevLB_4)
rand(lmer_NaivevLB_4)

lmer_NaivevLB_5 <- lmer(log(ZOI_area) ~  Temperature_centered +
                          Age_centered + Treatment+
                          Temperature_centered*Age_centered +
                          #Temperature_centered*Treatment+
                          #Age_centered*Treatment+
                          #Temperature_centered*Age_centered*Treatment+
                          #(1|Sample_ID)+
                          #(1|Plate_ID_Total)+
                          (1|Seeded_plates_OD),
                        data=Naive_and_LB,
                        REML=FALSE)

plot(lmer_NaivevLB_5)
summary(lmer_NaivevLB_5)
drop1(lmer_NaivevLB_5)
rand(lmer_NaivevLB_5)

lmer_NaivevLB_6 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                          Age_centered + Treatment+
                          Temperature_centered*Age_centered +
                          #Temperature_centered*Treatment+
                          #Age_centered*Treatment+
                          #Temperature_centered*Age_centered*Treatment+
                          #(1|Sample_ID)+
                          #(1|Plate_ID_Total)+
                          (1|Seeded_plates_OD),
                        data=Naive_and_LB,
                        REML=FALSE)

plot(lmer_NaivevLB_6)
summary(lmer_NaivevLB_6)
drop1(lmer_NaivevLB_6)
rand(lmer_NaivevLB_6)

AIC(lmer_NaivevLB_1,lmer_NaivevLB_2,lmer_NaivevLB_3,lmer_NaivevLB_4,lmer_NaivevLB_5,lmer_NaivevLB_6)
BIC(lmer_NaivevLB_1,lmer_NaivevLB_2,lmer_NaivevLB_3,lmer_NaivevLB_4,lmer_NaivevLB_5,lmer_NaivevLB_6)

#go with lmer_NaivevLB_5 for simplicity

############################
#now check model fits, plot emmeans, get model coefs and make comparisons:
#final models:
summary(lmer_NaivevLB_5)
summary(lmer_NaivevEcoli_4)
summary(lmer_NaivevMlut_6)


#find exponentiated model coefficients and confidence intervals
#save outputs:
modelvalues_NvLB <- parameters::model_parameters(lmer_NaivevLB_5, exponentiate = TRUE)
write_xlsx(modelvalues_NvLB, "Naive_v_LB_expcoefandCIs.xlsx")

modelvalues_NvLB_nonexp <- parameters::model_parameters(lmer_NaivevLB_5, exponentiate = FALSE)
write_xlsx(modelvalues_NvLB_nonexp, "Naive_v_LB_coefandCIs_nonexp.xlsx")


sink("Naive_v_LB_model_summary.txt")
print(summary(lmer_NaivevLB_5))
sink()  # returns output to the console

sink("Naive_v_LB_ANOVAoutput.txt")
print(Anova(lmer_NaivevLB_5,type="2",test.statistic = "Chisq"))
sink()

Naive_v_LB_ANOVA_param <- parameters::model_parameters(Anova(lmer_NaivevLB_5),type=2)
write_xlsx(Naive_v_LB_ANOVA_param,"Naive_v_LB_ANOVAexcel.xlsx")

####
#Ecoli:
modelvalues_NvEcoli <- parameters::model_parameters(lmer_NaivevEcoli_5, exponentiate = TRUE)
write_xlsx(modelvalues_NvEcoli, "Naive_v_Ecoli_expcoefandCIs.xlsx")

modelvalues_NvEcoli_nonexp <- parameters::model_parameters(lmer_NaivevEcoli_5, exponentiate = FALSE)
write_xlsx(modelvalues_NvEcoli_nonexp, "Naive_v_Ecoli_coefandCIs_nonexp.xlsx")


sink("Naive_v_Ecoli_model_summary.txt")
print(summary(lmer_NaivevEcoli_5))
sink()  # returns output to the console

sink("Naive_v_Ecoli_ANOVAoutput.txt")
print(Anova(lmer_NaivevEcoli_5,type="2",test.statistic = "Chisq"))
sink()

Naive_v_Ecoli_ANOVA_param <- parameters::model_parameters(Anova(lmer_NaivevEcoli_5),type=2)
write_xlsx(Naive_v_Ecoli_ANOVA_param,"Naive_v_Ecoli_ANOVAexcel.xlsx")

#####
#M luteus:
modelvalues_NvMlut <- parameters::model_parameters(lmer_NaivevMlut_6, exponentiate = TRUE)
write_xlsx(modelvalues_NvMlut, "Naive_v_Mlut_expcoefandCIs.xlsx")

modelvalues_NvMlut_nonexp <- parameters::model_parameters(lmer_NaivevMlut_6, exponentiate = FALSE)
write_xlsx(modelvalues_NvMlut_nonexp, "Naive_v_Mlut_coefandCIs_nonexp.xlsx")


sink("Naive_v_Mlut_model_summary.txt")
print(summary(lmer_NaivevMlut_6))
sink()  # returns output to the console

sink("Naive_v_Mlut_ANOVAoutput.txt")
print(Anova(lmer_NaivevMlut_6,type="2",test.statistic = "Chisq"))
sink()

Naive_v_Mlut_ANOVA_param <- parameters::model_parameters(Anova(lmer_NaivevMlut_6),type=2)
write_xlsx(Naive_v_Mlut_ANOVA_param,"Naive_v_Mlut_ANOVAexcel.xlsx")

###################################################
#https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/lmm.html
library(emmeans)

#LB emmens:
ref_grid(lmer_NaivevLB_5)
Naive_v_LB_emm <- emmeans(lmer_NaivevLB_5,specs=c("Temperature_centered", "Age_centered","Treatment"),
                          type="response", at = list(Temperature_centered = c(27,30,32),
                                                     Age_centered = c(1,5,10,15)))

#Naive_v_LB_emm <- emmeans(lmer_NaivevLB_5,specs=c("Temperature_centered", "Age_centered","Treatment"),
 #                         type="response", at = list(Temperature_centered = c(-1.18,0.24,1.18),
  #                                                   Age_centered = c(-1.24,-0.46,0.52,1.50)))
Naive_v_LB_emm
Naive_v_LB_emm_dataframe <- as.data.frame(Naive_v_LB_emm)
Naive_v_LB_emm_dataframe

str(Naive_v_LB_emm_dataframe)
head(Naive_v_LB_emm_dataframe)
write_xlsx(Naive_v_LB_emm_dataframe,"Naive_v_LB_emmeans.xlsx")


Naive_v_LB_TempandTreatment <- Naive_v_LB_emm_dataframe %>%
  group_by(Treatment,Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_LB_TempandTreatment

Naive_v_LB_TempAlone <- Naive_v_LB_emm_dataframe %>%
  group_by(Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_LB_TempAlone

Naive_v_LB_AgeandTreatment <- Naive_v_LB_emm_dataframe %>%
  group_by(Treatment,Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_LB_AgeandTreatment

Naive_v_LB_AgeAlone <- Naive_v_LB_emm_dataframe %>%
  group_by(Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_LB_AgeAlone


Naive_v_LB_TreatmentAlone <- Naive_v_LB_emm_dataframe %>%
  group_by(Treatment)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_LB_TreatmentAlone

write_xlsx(Naive_v_LB_TempAlone,"Naive_v_LB_TempAlone.xlsx")
write_xlsx(Naive_v_LB_TempandTreatment,"Naive_v_LB_TempandTreatment.xlsx")
write_xlsx(Naive_v_LB_AgeAlone,"Naive_v_LB_AgeAlone_emmeans.xlsx")
write_xlsx(Naive_v_LB_AgeandTreatment,"Naive_v_LB_AgeandTreatment_emmeans.xlsx")
write_xlsx(Naive_v_LB_TreatmentAlone,"Naive_v_LB_Treatmentalone.xlsx")

'''
#use this below if want to make planned contrasts
Naive_emm_27_1 <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_30_1 <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_32_1 <- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_27_5 <- c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_30_5 <- c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_32_5 <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_27_10 <- c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_30_10 <- c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_32_10 <- c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_27_15 <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_30_15 <- c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
Naive_emm_32_15 <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
LB_emm_27_1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
LB_emm_30_1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)
LB_emm_32_1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
LB_emm_27_5 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
LB_emm_30_5 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
LB_emm_32_5 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
LB_emm_27_10 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
LB_emm_30_10 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0)
LB_emm_32_10 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
LB_emm_27_15 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
LB_emm_30_15 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
LB_emm_32_15 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)


# exp5c_emm # print in console to get row numbers
# set the mean as the row number from the emmeans table
scr_basal <- c(1,0,0,0)
siNR4A3_basal <- c(0,1,0,0)
scr_eps <- c(0,0,1,0)
siNR4A3_eps <- c(0,0,0,1)

exp5c_m1_planned <- contrast(exp5c_m1_emm,
                             method = list(
                               "(Scr EPS) - (Scr Basal)" = c(scr_eps - scr_basal),
                               "(siNR4A3 EPS) - (siNR4A3 Basal)" = c(siNR4A3_eps - siNR4A3_basal),
                               "Interaction" = c(siNR4A3_eps - siNR4A3_basal) -
                                 c(scr_eps - scr_basal)
                               
                             ),
                             adjust = "none"
) %>%
  summary(infer = TRUE)

#https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/lmm#plotting-models-fit-to-batched-data
'''
##


#create labels:
templabs <- c("27˚C","30˚C","32˚C")
#names(templabs)<- c("-1.18","0.24","1.18")
names(templabs)<- c("27","30","32")

agelabs <- c("1 day","5 days", "10 days", "15 days")
names(agelabs)<-c("1","5","10","15")
#names(agelabs)<-c("-1.24","-0.46","0.52","1.5")

NLB_treatlabs <-c("Naïve","Injury")
names(NLB_treatlabs)<-c("Naïve","LB")

#make pairwise contrasts:
pairwisecontrasts.NaivevLB <- contrast(Naive_v_LB_emm, 
                                        type="response",
                                        method = "pairwise", 
                                        adjust = "sidak") %>%
  summary(infer=TRUE)

pairwisecontrasts.NaivevLB 
as.data.frame(pairwisecontrasts.NaivevLB)
head(pairwisecontrasts.NaivevLB)

pairwisecontrasts.NaivevLB_reverse <- contrast(Naive_v_LB_emm, 
                                                type="response",
                                                method = "revpairwise", 
                                                adjust = "sidak") %>%
  summary(infer=TRUE)


write_xlsx(pairwisecontrasts.NaivevLB, "pairwisecontrasts.NaivevLB.xlsx")
write_xlsx(pairwisecontrasts.NaivevLB_reverse, "pairwisecontrasts.NaivevLB_reverse.xlsx")


#plot emmeans:
str(Naive_v_LB_emm_dataframe)
Naive_v_LB_emm_dataframe$Age_centered <- as.factor(Naive_v_LB_emm_dataframe$Age_centered)
Naive_v_LB_emm_dataframe$Temperature_centered <- as.factor(Naive_v_LB_emm_dataframe$Temperature_centered)
 

png(filename = "Naive_v_LB_emmeans_tempattop.png", width = 6, height = 5, units = "in", res = 300)
Naive_v_LB_emm_dataframe %>% 
ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(color=Treatment,fill=Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(Treatment~Temperature_centered, labeller = labeller(Treatment=NLB_treatlabs,Temperature_centered=templabs))+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(breaks = c(-1.24,-0.46,0.52,1.50),
                  labels = c("1","5","10","15"),
                 name="Adult Age (days)")+
  scale_color_manual(name="Treatment",
                     labels = c("Naïve","Injury"),
                     values=c("#F8766D","#7CAE00"))+
  scale_fill_manual(name="Treatment",
                     labels = c("Naïve","Injury"),
                     values=c("#F8766D","#7CAE00"))+
    ylab(expression("Estimated Area of Zone of Inhibition"~(mm^2)~""))+ 
  scale_y_continuous(limits=c(0,40))+
  theme_pubr()+ 
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

png(filename = "Naive_v_LB_emmeans_ageattop.png", width = 6, height = 5, units = "in", res = 300)
Naive_v_LB_emm_dataframe %>% 
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(Treatment~Age_centered, labeller = labeller(Treatment=NLB_treatlabs,Age_centered=agelabs,Temperature_centered=templabs))+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(breaks = c(-1.18,0.24,1.18),
                   labels = c("27","30","32"),
                   name="Temperature (˚C)")+
  scale_color_manual(name="Treatment",
                     labels = c("Naïve","Injury"),
                     values=c("#F8766D","#7CAE00"))+
  scale_fill_manual(name="Treatment",
                    labels = c("Naïve","Injury"),
                    values=c("#F8766D","#7CAE00"))+
  ylab(expression("Estimated Area of Zone of Inhibition"~(mm^2)~""))+ 
  scale_y_continuous(limits=c(0,40))+
  theme_pubr()+ 
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

########################################################
#Naive v Ecoli:
ref_grid(lmer_NaivevEcoli_5)

Naive_v_Ecoli_emm <- emmeans(lmer_NaivevEcoli_5,specs=c("Temperature_centered", "Age_centered","Treatment"),
                          type="response", at = list(Temperature_centered = c(27,30,32),
                                                     Age_centered = c(1,5,10,15)))
#Naive_v_Ecoli_emm <- emmeans(lmer_NaivevEcoli_5,specs=c("Temperature_centered", "Age_centered","Treatment"),
 #                            type="response", at = list(Temperature_centered = c(-1.18,0.24,1.18),
  #                                                      Age_centered = c(-1.24,-0.46,0.52,1.50)))

Naive_v_Ecoli_emm
Naive_v_Ecoli_emm_dataframe <- as.data.frame(Naive_v_Ecoli_emm)
Naive_v_Ecoli_emm_dataframe

str(Naive_v_Ecoli_emm_dataframe)
head(Naive_v_Ecoli_emm_dataframe)
write_xlsx(Naive_v_Ecoli_emm_dataframe,"Naive_v_Ecoli_emmeans.xlsx")
######
Naive_v_Ecoli_TempandTreatment <- Naive_v_Ecoli_emm_dataframe %>%
  group_by(Treatment,Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Ecoli_TempandTreatment

Naive_v_Ecoli_TempAlone <- Naive_v_Ecoli_emm_dataframe %>%
  group_by(Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Ecoli_TempAlone

Naive_v_Ecoli_AgeandTreatment <- Naive_v_Ecoli_emm_dataframe %>%
  group_by(Treatment,Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Ecoli_AgeandTreatment

Naive_v_Ecoli_AgeAlone <- Naive_v_Ecoli_emm_dataframe %>%
  group_by(Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Ecoli_AgeAlone


Naive_v_Ecoli_TreatmentAlone <- Naive_v_Ecoli_emm_dataframe %>%
  group_by(Treatment)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Ecoli_TreatmentAlone

write_xlsx(Naive_v_Ecoli_TempAlone,"Naive_v_Ecoli_TempAlone_emmeans.xlsx")
write_xlsx(Naive_v_Ecoli_TempandTreatment,"Naive_v_Ecoli_TempandTreatment_Emmeans.xlsx")
write_xlsx(Naive_v_Ecoli_AgeAlone,"Naive_v_Ecoli_AgeAlone_emmeans.xlsx")
write_xlsx(Naive_v_Ecoli_AgeandTreatment,"Naive_v_Ecoli_AgeandTreatment_emmeans.xlsx")
write_xlsx(Naive_v_Ecoli_TreatmentAlone,"Naive_v_Ecoli_Treatmentalone_emmeans.xlsx")


#make pairwise contrasts:
pairwisecontrasts.NaivevEcoli <- contrast(Naive_v_Ecoli_emm, 
                                       type="response",
                                       method = "pairwise", 
                                       adjust = "sidak") %>%
  summary(infer=TRUE)

pairwisecontrasts.NaivevEcoli
as.data.frame(pairwisecontrasts.NaivevEcoli)
head(pairwisecontrasts.NaivevEcoli)

pairwisecontrasts.NaivevEcoli_reverse <- contrast(Naive_v_Ecoli_emm, 
                                               type="response",
                                               method = "revpairwise", 
                                               adjust = "sidak") %>%
  summary(infer=TRUE)


write_xlsx(pairwisecontrasts.NaivevEcoli, "pairwisecontrasts.NaivevEcoli.xlsx")
write_xlsx(pairwisecontrasts.NaivevEcoli_reverse, "pairwisecontrasts.NaivevEcoli_reverse.xlsx")

#graph:
Naive_v_Ecoli_emm_dataframe$Age_centered <- as.factor(Naive_v_Ecoli_emm_dataframe$Age_centered)
Naive_v_Ecoli_emm_dataframe$Temperature_centered <- as.factor(Naive_v_Ecoli_emm_dataframe$Temperature_centered)

#create labels:
Naive_v_Ecoli_emm_dataframe_italics <- Naive_v_Ecoli_emm_dataframe

Naive_v_Ecoli_emm_dataframe_italics$Treatment <- factor(Naive_v_Ecoli_emm_dataframe_italics$Treatment,    # Change factor labels
                                        labels = c("Naïve","italic(`E. coli`)"))

Naive_v_Ecoli_emm_dataframe_italics$Temperature_centered <- factor(Naive_v_Ecoli_emm_dataframe_italics$Temperature_centered,
                                          labels = c("`27˚C`","`30˚C`","`32˚C`"))
Naive_v_Ecoli_emm_dataframe_italics$Age_centered <- factor(Naive_v_Ecoli_emm_dataframe_italics$Age_centered,
                                                           labels = c("1","5","10","15"))
png(filename = "Naive_v_Ecoli_emmeans_tempattop.png", width = 6, height = 5, units = "in", res = 300)
Naive_v_Ecoli_emm_dataframe_italics %>% 
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(color=Treatment,fill=Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(Treatment~Temperature_centered, labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(name="Adult Age (days)")+
  scale_color_manual(name="Treatment",
                     values=c("#F8766D","#00BFC4"))+
  scale_fill_manual(name="Treatment",
                    values=c("#F8766D","#00BFC4"))+
  ylab(expression("Estimated Area of Zone of Inhibition"~(mm^2)~""))+ 
  theme_pubr()+ 
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,60))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()


Naive_v_Ecoli_emm_dataframe_italics$Age_centered <- factor(Naive_v_Ecoli_emm_dataframe_italics$Age_centered,
                                                           labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
Naive_v_Ecoli_emm_dataframe_italics$Temperature_centered <- factor(Naive_v_Ecoli_emm_dataframe_italics$Temperature_centered,
                                                                   labels = c("27","30","32"))

png(filename = "Naive_v_Ecoli_emmeans_ageattop.png", width = 6, height = 5, units = "in", res = 300)
Naive_v_Ecoli_emm_dataframe_italics %>% 
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(color=Treatment,fill=Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(Treatment~Age_centered, labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(name="Temperature (˚C)")+
  scale_color_manual(name="Treatment",
                     values=c("#F8766D","#00BFC4"))+
  scale_fill_manual(name="Treatment",
                    values=c("#F8766D","#00BFC4"))+
  ylab(expression("Estimated Area of Zone of Inhibition"~(mm^2)~""))+ 
  theme_pubr()+ 
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,50))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

######################
#M luteus:
#Naive v Mluteus:
ref_grid(lmer_NaivevMlut_6)

Naive_v_Mluteus_emm <- emmeans(lmer_NaivevMlut_6,specs=c("Temperature_centered", "Age_centered","Treatment"),
                             type="response", at = list(Temperature_centered = c(-1.18,0.24,1.18),
                                                        Age_centered = c(-1.24,-0.46,0.52,1.50)))
Naive_v_Mluteus_emm 
Naive_v_Mluteus_emm_dataframe <- as.data.frame(Naive_v_Mluteus_emm)


str(Naive_v_Mluteus_emm_dataframe)
head(Naive_v_Mluteus_emm_dataframe)

Naive_v_Mluteus_TempandTreatment <- Naive_v_Mluteus_emm_dataframe %>%
  group_by(Treatment,Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Mluteus_TempandTreatment

Naive_v_Mluteus_TempAlone <- Naive_v_Mluteus_emm_dataframe %>%
  group_by(Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Mluteus_TempAlone

Naive_v_Mluteus_AgeandTreatment <- Naive_v_Mluteus_emm_dataframe %>%
  group_by(Treatment,Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Mluteus_AgeandTreatment

Naive_v_Mluteus_AgeAlone <- Naive_v_Mluteus_emm_dataframe %>%
  group_by(Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Mluteus_AgeAlone


Naive_v_Mluteus_TreatmentAlone <- Naive_v_Mluteus_emm_dataframe %>%
  group_by(Treatment)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

Naive_v_Mluteus_TreatmentAlone

write_xlsx(Naive_v_Mluteus_TempAlone,"Naive_v_Mluteus_TempAlone_emmeans.xlsx")
write_xlsx(Naive_v_Mluteus_TempandTreatment,"Naive_v_Mluteus_TempandTreatment_Emmeans.xlsx")
write_xlsx(Naive_v_Mluteus_AgeAlone,"Naive_v_Mluteus_AgeAlone_emmeans.xlsx")
write_xlsx(Naive_v_Mluteus_AgeandTreatment,"Naive_v_Mluteus_AgeandTreatment_emmeans.xlsx")
write_xlsx(Naive_v_Mluteus_TreatmentAlone,"Naive_v_Mluteus_Treatmentalone_emmeans.xlsx")



write_xlsx(Naive_v_Mluteus_emm_dataframe,"Naive_v_Mluteus_emmeans.xlsx")

#make pairwise contrasts:
pairwisecontrasts.NaivevMluteus <- contrast(Naive_v_Mluteus_emm, 
                                          type="response",
                                          method = "pairwise", 
                                          adjust = "sidak") %>%
  summary(infer=TRUE)

pairwisecontrasts.NaivevMluteus
as.data.frame(pairwisecontrasts.NaivevMluteus)
head(pairwisecontrasts.NaivevMluteus)

pairwisecontrasts.NaivevMluteus_reverse <- contrast(Naive_v_Mluteus_emm, 
                                                  type="response",
                                                  method = "revpairwise", 
                                                  adjust = "sidak") %>%
  summary(infer=TRUE)


write_xlsx(pairwisecontrasts.NaivevMluteus, "pairwisecontrasts.NaivevMluteus.xlsx")
write_xlsx(pairwisecontrasts.NaivevMluteus_reverse, "pairwisecontrasts.NaivevMluteus_reverse.xlsx")

#graph:
Naive_v_Mluteus_emm_dataframe$Age_centered <- as.factor(Naive_v_Mluteus_emm_dataframe$Age_centered)
Naive_v_Mluteus_emm_dataframe$Temperature_centered <- as.factor(Naive_v_Mluteus_emm_dataframe$Temperature_centered)

#labels:
Naive_v_Mluteus_emm_dataframe_italics <- Naive_v_Mluteus_emm_dataframe

Naive_v_Mluteus_emm_dataframe_italics$Treatment <- factor(Naive_v_Mluteus_emm_dataframe_italics$Treatment,    # Change factor labels
                                                        labels = c("Naïve","italic(`M. luteus`)"))

Naive_v_Mluteus_emm_dataframe_italics$Temperature_centered <- factor(Naive_v_Mluteus_emm_dataframe_italics$Temperature_centered,
                                                                   labels = c("`27˚C`","`30˚C`","`32˚C`"))
Naive_v_Mluteus_emm_dataframe_italics$Age_centered <- factor(Naive_v_Mluteus_emm_dataframe_italics$Age_centered,
                                                           labels = c("1","5","10","15"))
png(filename = "Naive_v_Mluteus_emmeans_tempattop.png", width = 6, height = 5, units = "in", res = 300)
Naive_v_Mluteus_emm_dataframe_italics %>% 
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(color=Treatment,fill=Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(Treatment~Temperature_centered, labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(name="Adult Age (days)")+
  scale_color_manual(name="Treatment",
                     values=c("#F8766D","#C77CFF"))+
  scale_fill_manual(name="Treatment",
                    values=c("#F8766D","#C77CFF"))+
  ylab(expression("Estimated Area of Zone of Inhibition"~(mm^2)~""))+ 
  scale_y_continuous(limits=c(0,40))+
  theme_pubr()+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()


Naive_v_Mluteus_emm_dataframe_italics$Age_centered <- factor(Naive_v_Mluteus_emm_dataframe_italics$Age_centered,
                                                           labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
Naive_v_Mluteus_emm_dataframe_italics$Temperature_centered <- factor(Naive_v_Mluteus_emm_dataframe_italics$Temperature_centered,
                                                                   labels = c("27","30","32"))

png(filename = "Naive_v_Mluteus_emmeans_ageattop.png", width = 6, height = 5, units = "in", res = 300)
Naive_v_Mluteus_emm_dataframe_italics %>% 
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(color=Treatment,fill=Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(Treatment~Age_centered, labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(name="Temperature (˚C)")+
  scale_color_manual(name="Treatment",
                     values=c("#F8766D","#C77CFF"))+
  scale_fill_manual(name="Treatment",
                    values=c("#F8766D","#C77CFF"))+
  ylab(expression("Estimated Area of Zone of Inhibition"~(mm^2)~""))+ 
  theme_pubr()+ 
  scale_y_continuous(limits=c(0,40))+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()


#################################

#effect sizes:
library(effectsize)

LB_effectsize_partial <- eta_squared(lmer_NaivevLB_5,partial=TRUE)
LB_effectsize_partial_df <- as.data.frame(LB_effectsize_partial)
LB_effectsize_partial_df
#save partial effect size:
write_xlsx(LB_effectsize_partial_df,"LB_effectsize_partial.xlsx")



#effect sizes for mixed models:
effects<-parameters::model_parameters(lmer_NaivevEcoli_5, effects = "fixed", ci_method = "satterthwaite")
effects <-as.data.frame(effects)
effects
Naive_v_Ecoli_ANOVA_param
t_to_eta2(2.5495358, df_error = 2)


#######################################
#combine all and graph emmeans:

allem_agetreat <- rbind(Naive_v_LB_emm_dataframe,Naive_v_Ecoli_emm_dataframe,Naive_v_Mluteus_emm_dataframe) %>%
  group_by(Treatment, Age_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))
allem_agetreat

allem_temptreat <- rbind(Naive_v_LB_emm_dataframe,Naive_v_Ecoli_emm_dataframe,Naive_v_Mluteus_emm_dataframe) %>%
  group_by(Treatment, Temperature_centered)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))

allem_treatalone <- rbind(Naive_v_LB_emm_dataframe,Naive_v_Ecoli_emm_dataframe,Naive_v_Mluteus_emm_dataframe) %>%
  group_by(Treatment)%>%
  reframe(mean_area = mean(response),
          SE_area = sd(response)/sqrt(n()))
allem_treatalone

allem_allfactors <- rbind(Naive_v_LB_emm_dataframe,Naive_v_Ecoli_emm_dataframe,Naive_v_Mluteus_emm_dataframe)

#######
allem_agetreat$Treatment <- factor(allem_agetreat$Treatment,    # Change factor labels
                                                          labels = c("Naïve","Injury","italic(`E. coli`)","italic(`M. luteus`)"))
allem_agetreat$Age_centered <- factor(allem_agetreat$Age_centered,
                                                             labels = c("1","5","10","15"))


allem_temptreat$Temperature_centered <- factor(allem_temptreat$Temperature_centered,
                                                                     labels = c("27","30","32"))
allem_temptreat$Treatment <- factor(allem_temptreat$Treatment,    # Change factor labels
                                       labels = c("Naïve","Injury","italic(`E. coli`)","italic(`M. luteus`)"))


allem_treatalone$Treatment<- factor(allem_treatalone$Treatment,    # Change factor labels
                          labels = c("Naïve","Injury","italic(`E. coli`)","italic(`M. luteus`)"))

allem_treatalone

allem_temptreat

png(filename = "ZOI_emmeans_TEMP_effect.png",width = 4, height = 4, units = "in", res = 300)
emp1 <- allem_temptreat %>% 
  ggplot(aes(x=Temperature_centered,y=mean_area))+
  geom_bar(aes(color=Treatment,fill=Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  facet_grid(~Treatment, labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
                width=0.6,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(name="Temperature (˚C)")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  theme_pubr()+ 
  scale_y_continuous(limits=c(0,40))+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
emp1
dev.off()

allem_agetreat
png(filename = "ZOI_emmeans_Age_effect.png",width = 4, height = 4, units = "in", res = 300)
emp2<-allem_agetreat%>%
  group_by(Treatment)%>%
  ggplot(aes(x=Age_centered,y=mean_area))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(~Treatment,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
                width=0.4,position=position_dodge(1),
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Adult Age (days)") +
  #geom_jitter(data=ZOI_Temperature_italics, aes(x=Temperature,y=ZOI_area),#color=ZOI_italic$Technical_Rep,
  #           position = position_dodge(0.5),size=1)+
  theme_pubr()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y = element_text(size=10))+
  theme(legend.position = "none")
emp2
dev.off()

allem_treatalone
png(filename = "ZOI_emmeans_Treatment_effect.png",width = 4, height = 4, units = "in", res = 300)
emp3<-allem_treatalone %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Treatment,y=mean_area))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  #facet_grid(~Treatment,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
                width=0.4,position=position_dodge(1),
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Immune Treatment") +
  scale_x_discrete(limits = c("Naïve","Injury","italic(`E. coli`)","italic(`M. luteus`)"),labels=c("Naïve","Injury",bquote(italic("E. coli")), bquote(italic("M. luteus"))))+
  #geom_jitter(data=ZOI_Temperature_italics, aes(x=Temperature,y=ZOI_area),#color=ZOI_italic$Technical_Rep,
  #           position = position_dodge(0.5),size=1)+
  theme_pubr()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y = element_text(size=10))+
  theme(legend.position = "none")
emp3
dev.off()

p123 <- ggarrange(emp1,emp2,emp3,
                  ncol = 3, nrow = 1,
                  heights = c(12,12,12),
                  labels = c("A", "B", "C"))
p123

png(filename = "ZOI_all_effects_emmeans.png",width = 12, height =5, units = "in", res = 300)
p123
dev.off()
