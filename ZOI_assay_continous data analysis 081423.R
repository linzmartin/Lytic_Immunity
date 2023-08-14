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
library(ggpubr)

#library(emmeans)


#######################################################
#import the data and clean it up:


#import the data:
ZOI_data <- read_xlsx("ZOI_assay_sample_log_081123_LEM.xlsx")

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

unique(ZOI_data$Age_centered)
mean(ZOI_data$Age_centered) #approximates zero
mean(ZOI_data$Age)


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
normtest.glance #values are not normal!

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

ZOI_data$log_area <- log(ZOI_data$ZOI_area)

ggplot(ZOI_data, aes(log(ZOI_area))) + geom_histogram()

#ggplot(ZOI_data, aes(log(ZOI_area+0.1))) + geom_histogram()


#model + anova
library(lme4)
lm_full <- lmer(log(ZOI_area) ~ Temperature+Age+Treatment + 
                  Temperature*Age+
                  Temperature*Treatment+
                  Age*Treatment+
                  Temperature*Age*Treatment+
              (1|Sample_ID)+(1|Plate_ID_Total)+(1|Batch_Number),
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

#ZOIcontinuous <- ZOI_data

#ZOIcontinuous$Temperature <- as.numeric(ZOIcontinuous$Temperature)
#ZOIcontinuous$Age <- as.numeric(ZOIcontinuous$Age)
#str(ZOIcontinuous)

#ggplot(ZOIcontinuous, aes(log(ZOI_area))) + geom_histogram()

lm_full_continuous <- lmer(log(ZOI_area) ~ Temperature_centered+Age+Treatment + 
                  Temperature*Age+
                  Temperature*Treatment+
                  Age*Treatment+
                  Temperature*Age*Treatment+
                 # (1|Sample_ID)+
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




######################

#stratify categorical by treatment:

Naive_ZOI <- subset(ZOI_data,Treatment == "Naïve")
LB_ZOI <- subset(ZOI_data,Treatment == "LB")
Ecoli_ZOI <- subset(ZOI_data,Treatment == "E_coli")
Mlut_ZOI <- subset(ZOI_data,Treatment == "M_luteus")

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
                      #Temperature*Age+
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

modelvalues_naive <- parameters::model_parameters(lm_Naive_3, exponentiate = TRUE)
head(modelvalues_naive)

write_xlsx(modelvalues_naive, "Naive_model_expcoefandCIs.xlsx")

library(effectsize)

Naive_effectsize_partial <- eta_squared(lm_Naive_3,partial=TRUE)
Naive_effectsize_partial_dataframe <- as.data.frame(Naive_effectsize_partial)

write_xlsx(Naive_effectsize_partial_dataframe,"Naive_partialeffectsizes.xlsx")

#takeaway: temp affects naive, but not age or interaction

###############################
lmer_LB_1 <- lmer(log(ZOI_area) ~ Temperature_centered +
                    Age_centered + 
                    Temperature_centered*Age_centered+
                    (1|Sample_ID)+
                    (1|Plate_ID_Total)+
                    (1|Batch_Number),
                     data=LB_ZOI,
                     REML=FALSE)

summary(lmer_LB_1)
rand(lmer_LB_1) #can get rid of all random effects 

lm_LB_2 <- lm(log(ZOI_area) ~  Temperature_centered +
                Age_centered + 
                Temperature_centered*Age_centered,
                 data=LB_ZOI)

summary(lm_LB_2)
plot(lm_LB_2)

Anova(lm_LB_2,type=2)
drop1(lm_LB_2)

lm_LB_3 <- lm(log(ZOI_area) ~ Temperature_centered +
                Age_centered,
              data=LB_ZOI)

summary(lm_LB_3)
plot(lm_LB_3)

anova(lm_LB_2,lm_LB_3)
AIC(lm_LB_2,lm_LB_3)

#no effects on injury group


lm_LB_4 <- lm(log(ZOI_area) ~ poly(Temperature_centered,2) +
                poly(Age_centered,2),
              data=LB_ZOI)

plot(lm_LB_4)
AIC(lm_LB_3,lm_LB_4)

lm_LB5 <- glm(ZOI_area ~ Temperature_centered + Age_centered,
              data=LB_ZOI, family=gaussian(link="log"))

summary(lm_LB5)
plot(lm_LB5)
AIC(lm_LB_3,lm_LB5)

##
#save outputs:
sink("LB_summary.txt")
print(summary(lm_LB_3))
sink()  # returns output to the console

sink("LB_ANOVAoutput.txt")
print(Anova(lm_LB_3,type="2",test.statistic = "Chisq"))
sink()

modelvalues_LB <- parameters::model_parameters(lm_LB_3, exponentiate = TRUE)
head(modelvalues_LB)

write_xlsx(modelvalues_LB, "LB_model_expcoefandCIs.xlsx")

LB_effectsize_partial <- eta_squared(lm_LB_3,partial=TRUE)
LB_effectsize_partial_dataframe <- as.data.frame(LB_effectsize_partial)

write_xlsx(LB_effectsize_partial_dataframe,"LB_partialeffectsizes.xlsx")



###########################


#infection:
lmer_Ecoli_1 <- lmer(log(ZOI_area) ~  Temperature_centered +
                       Age_centered + 
                       Temperature_centered*Age_centered+
                       (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Batch_Number),
                     data=Ecoli_ZOI,
                     REML=FALSE)

summary(lmer_Ecoli_1)

rand(lmer_Ecoli_1) #get rid of all random effects


lmer_Ecoli_2 <- lm(log(ZOI_area) ~  Temperature_centered +
                     Age_centered + 
                     Temperature_centered*Age_centered,
                       #(1|Sample_ID)+
                       #(1|Plate_ID),
                       #(1|Batch_Number),
                     data=Ecoli_ZOI)

summary(lmer_Ecoli_2)
Anova(lmer_Ecoli_2,type=2)
plot(lmer_Ecoli_2)

drop1(lmer_Ecoli_2)

lmer_Ecoli_3 <- lm(log(ZOI_area) ~  Temperature_centered +
                     Age_centered,  
                    # Temperature_centered*Age_centered,
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)


summary(lmer_Ecoli_3)
Anova(lmer_Ecoli_3,type=2)
plot(lmer_Ecoli_3)

AIC(lmer_Ecoli_2, lmer_Ecoli_3)


lmer_Ecoli_4 <- lm(log(ZOI_area) ~  poly(Temperature_centered,2) +
                     Age_centered,  
                   # Temperature_centered*Age_centered,
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)
#temp has impact on E. coli zoi but not age or int of temp and age

summary(lmer_Ecoli_4)

lmer_Ecoli_5 <- lm(log(ZOI_area) ~  Temperature_centered +
                     poly(Age_centered,2),
                   # Temperature_centered*Age_centered,
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)

summary(lmer_Ecoli_5)
plot(lmer_Ecoli_5)
Anova(lmer_Ecoli_5,type=2)

AIC(lmer_Ecoli_5)

lmer_Ecoli_6 <- lm(log(ZOI_area) ~  poly(Temperature_centered,2) +
                     poly(Age_centered,2),
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)

summary(lmer_Ecoli_6)

lmer_Ecoli_7 <- lm(log(ZOI_area) ~  poly(Temperature_centered,2) +
                     poly(Age_centered,2)+
                   poly(Temperature_centered*Age_centered,2),
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)

summary(lmer_Ecoli_7)

AIC(lmer_Ecoli_3,lmer_Ecoli_4,lmer_Ecoli_5,lmer_Ecoli_6,lmer_Ecoli_7)


lmer_Ecoli_8 <- lm(log(ZOI_area) ~  poly(Temperature_centered,2) +
                   poly(Age_centered,3),
                   #(1|Sample_ID)+
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)

AIC(lmer_Ecoli_3,lmer_Ecoli_4,lmer_Ecoli_5,lmer_Ecoli_6,lmer_Ecoli_7,lmer_Ecoli_8)

library(modelr)
grid <- Ecoli_ZOI %>% 
  data_grid(Temperature_centered, .model = lmer_Ecoli_3) %>% 
  add_predictions(lmer_Ecoli_3)

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
  data_grid(Temperature_centered, .model = lmer_Ecoli_3) %>% 
  add_predictions(lmer_Ecoli_3)

ggplot(Ecoli_ZOI, aes(x=Temperature_centered, 
                      y=log(ZOI_area),group=Temperature_centered)) + 
 geom_boxplot()+ 
  geom_point(aes(x=Temperature_centered,y=pred),data=grid,colour="red",size=2)


grid_age <- Ecoli_ZOI %>% 
  data_grid(Age_centered, .model = lmer_Ecoli_3) %>% 
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

ggplot(daily, aes(wday, n)) +
  geom_boxplot() + 
  geom_point(data = grid, colour = "red") + 
  facet_wrap(~ term)

ggplot(Ecoli_ZOI, aes(x=Age_centered, 
                      y=log(ZOI_area),group=Age_centered)) + 
  geom_boxplot()+
  geom_point(aes(x=Age_centered,y=log_ZOI_area),data=grid_tandage,colour="red",size=2)+
  facet_wrap(~Temperature_centered)

###
library(MASS)
mod3 <- MASS::rlm(n ~ wday * term, data = daily)

mod3 <- MASS::rlm(ZOI_area ~ Temperature_centered*Age_centered,
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
                                     Age + I(Age^2)+
                                     #Temperature*Age+
                                     #(1|Sample_ID)+
                                     (1|Plate_ID_Total),
                                   #(1|Batch_Number),
                                   data=Mlut_ZOI,
                                   REML=FALSE)

summary(lmer_Mlut_3)

plot(lmer_Mlut_3)

Anova(lmer_Mlut_3)

#temperature and age both significantly shape m.lut response but do not interact together
########################


#################################
#################################

#import the data:
gene_expression_data <- read_xlsx("lytic_expression.xlsx",sheet="Cecropin_ref")

str(gene_expression_data)

gene_expression_data$Age <- as.factor(gene_expression_data$Age)
gene_expression_data$Temperature <- as.factor(gene_expression_data$Temperature)
gene_expression_data$Treatment <- as.factor(gene_expression_data$Treatment)
gene_expression_data$Gene <- as.factor(gene_expression_data$Gene)

str(gene_expression_data)
gene_expression_data <- as.data.frame(gene_expression_data)



###########
#make Naive and 1day the reference groups:

#relevel treatment to make Naive come first:
gene_expression_data <- gene_expression_data %>%
  mutate(Treatment = fct_relevel(Treatment, 
                                 "Naïve","Live_Ecoli"))

gene_expression_data <- gene_expression_data %>%
  mutate(Age = fct_relevel(Age, "1","10"))

#find average of two trials:

gene_avg_summary <- gene_expression_data %>%
  group_by(Treatment,Age,Temperature,Gene)%>%
  reframe(meanexp = mean(Relative_Expression),
    sd_exp = sd(Relative_Expression),
    n_exp = n(),
    SE_exp = sd(Relative_Expression)/sqrt(n()))

gene_avg_summary <- as.data.frame(gene_avg_summary)
str(gene_avg_summary)
    

#italicize:
##
gene_expression_data_italic <- gene_expression_data

gene_expression_data_italic$Treatment <- factor(gene_expression_data_italic$Treatment,    # Change factor labels
                               labels = c("Naïve",
                                          "italic(`E. coli`)"))
gene_expression_data_italic$Age <- factor(gene_expression_data_italic$Age,
                         labels = c("1","10"))

gene_expression_data_italic$Temperature <- factor(gene_expression_data_italic$Temperature,
                                 labels = c("`27˚C`","`30˚C`"))


gene_avg_summary_italic <- gene_avg_summary
gene_avg_summary_italic$Treatment <- factor(gene_avg_summary_italic$Treatment,    # Change factor labels
                                                labels = c("Naïve",
                                                           "italic(`E. coli`)"))
gene_avg_summary_italic$Age <- factor(gene_avg_summary_italic$Age,
                                          labels = c("1","10"))

gene_avg_summary_italic$Temperature <- factor(gene_avg_summary_italic$Temperature,
                                                  labels = c("`27˚C`","`30˚C`"))



####
#plots with italic labels (use these!):

treatlabs <-c("E. coli","Naïve")
names(treatlabs)<-c("E_coli","Naïve")

#labels:
gene_avg_summary_italic <- gene_avg_summary
gene_avg_summary_italic$Treatment <- factor(gene_avg_summary_italic$Treatment,    # Change factor labels
                                            labels = c("Naïve",
                                                       "italic(`E. coli`)"))
gene_avg_summary_italic$Age <- factor(gene_avg_summary_italic$Age,
                                      labels = c("1","10"))

gene_avg_summary_italic$Temperature <- factor(gene_avg_summary_italic$Temperature,
                                              labels = c("`27˚C`","`30˚C`"))


gene_avg_summary_italic$Gene <- factor(gene_avg_summary_italic$Gene,
                                       labels=c("italic(`CECA`)",
                                       "italic(`RPS17`)"))

png(filename = "gene_expression_cecropin.png", width = 5, height = 4, units = "in", res = 300)
gene_avg_summary_italic %>%
  group_by(Gene)%>%
  ggplot(aes(x=Age,y=meanexp))+
  geom_bar(aes(fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_fill_manual(name="Treatment",
                    labels=c("Naïve", expression(italic("E. coli"))), 
                    values=c("tomato2", "darkturquoise"))+
  facet_grid(Gene~Temperature,labeller = label_parsed)+
  geom_errorbar(aes(ymin=meanexp - SE_exp,
                    ymax=meanexp + SE_exp,
                    group=Treatment),
                width=0.5,position=position_dodge(1),
                color="black")+
  ylab("Relative Gene Expression")+ 
  xlab("Age (days)") +
  theme_pubr()+
  theme(legend.position = "right")
dev.off()
