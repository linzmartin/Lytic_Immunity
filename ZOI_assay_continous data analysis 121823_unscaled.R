#lytic immunity - ZOI analysis
#individual: naive, injured, infected
#last updated: 01/03/23 LEM
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
library(EnvStats)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(rstatix)
library(car) 
library(ggpubr)
library(lme4)
library(broom)
library(rcompanion)
library(modelr)
library(emmeans)
library(lmerTest)
library(effectsize)
library(ggsci)
library(multcompView)
library(emmeans)
library(multcomp)
#######################################################
#import the data and clean it up:


#import the data:
ZOI_data <- read_xlsx("ZOI_assay_sample_log_092723_LEM.xlsx")

str(ZOI_data)
head(ZOI_data)

#variables of interest:


#format the data: 
#data is numeric / continuous
ZOI_data$Age <- as.numeric(ZOI_data$Age)
ZOI_data$Temperature <- as.numeric(ZOI_data$Temperature)
ZOI_data$Treatment <- as.factor(ZOI_data$Treatment)
ZOI_data$Diameter_avg <- as.numeric(ZOI_data$Diameter_avg)
ZOI_data$ZOI_area <- as.numeric(ZOI_data$ZOI_area)
str(ZOI_data)

ZOI_data$Sample_ID <- as.factor(ZOI_data$Sample_ID)

levels(ZOI_data$Treatment)
str(ZOI_data)

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

levels(ZOI_data$Treatment)
levels(ZOI_data$Age_centered)

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
  group_by(Treatment,Age_centered,Temperature_centered) %>%
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
str(ZOI_summary)
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
  group_by(Temperature_centered,Treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ZOI_area)))
head(normtestresults)

library(broom)
normtest.glance <- normtestresults %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
normtest.glance 
#values are not normal! p < 0.05.

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

ZOI_data$log_area <- log((ZOI_data$ZOI_area)+1) #go with this

ggplot(ZOI_data, aes((ZOI_area))) + geom_histogram()

ggplot(ZOI_data, aes(log(ZOI_area))) + geom_histogram()
ggplot(ZOI_data,aes(log_area))+geom_histogram()

ggplot(ZOI_data, aes(log(ZOI_area+0.1))) + geom_histogram()
ggplot(ZOI_data, aes(log(ZOI_area+1))) + geom_histogram()

#does log transformation achieve normality?

normtestresults <- ZOI_data %>%
  # mutate(mean_area = Delta_OD) %>%
  group_by(Temperature_centered, Treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$log_area)))
head(normtestresults)

library(broom)
normtest.glance <- normtestresults %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
normtest.glance
#mostly normal now

normtestresults <- ZOI_data %>%
  # mutate(mean_area = Delta_OD) %>%
  group_by(Age_centered,Treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$log_area)))
head(normtestresults)

library(broom)
normtest.glance <- normtestresults %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
normtest.glance
#also mostly normal

#model + anova
library(lme4)
lm_full <- lmer(log(ZOI_area) ~ 
                  Treatment + 
                  poly(Temperature_centered,2)+
                  poly(Age_centered,2)+
                  Temperature_centered*Treatment+
                  Age_centered*Treatment+
                  poly(Temperature_centered,2)*poly(Age_centered,2)+
                  #Temperature_centered*Age_centered*Treatment+
                  poly(Temperature_centered,2)*poly(Age_centered,2)*Treatment+
                  (1|Sample_ID)+
                  (1|Plate_ID_Total)+
                  (1|Batch_Number),
                data=ZOI_data,
                REML=FALSE)
coef(lm_full)
summary(lm_full) #can possibly drop a random effect
drop1(lm_full,test="Chisq") #drop 3 way int
rand(lm_full) #drop sample ID (technical rep)
Anova(lm_full,type=2)


lm_2 <- lmer(log(ZOI_area) ~ 
                  Treatment + 
                  poly(Temperature_centered,2)+
                  poly(Age_centered,2)+
                  Temperature_centered*Treatment+
                  Age_centered*Treatment+
                  poly(Temperature_centered,2)*poly(Age_centered,2)+
                  Temperature_centered*Age_centered*Treatment+
                  (1|Sample_ID)+
                  (1|Plate_ID_Total)+
                  (1|Batch_Number),
                data=ZOI_data,
                REML=FALSE)

summary(lm_2) #can possibly drop a random effect
drop1(lm_2,test="Chisq") #drop 3 way
rand(lm_2)
anova(lm_full,lm_2)


lm_3 <- lmer(log(ZOI_area) ~ 
               Treatment + 
               poly(Temperature_centered,2)+
               poly(Age_centered,2)+
               Temperature_centered*Treatment+
               Age_centered*Treatment+
               poly(Temperature_centered,2)*poly(Age_centered,2)+
               #Temperature_centered*Age_centered*Treatment+
              # (1|Sample_ID)+
               (1|Plate_ID_Total)+
               (1|Batch_Number),
             data=ZOI_data,
             REML=FALSE)

summary(lm_3) #can possibly drop a random effect
drop1(lm_3,test="Chisq") #drop treat*age
rand(lm_3) #leave remaining random effects in model
anova(lm_2,lm_3,lm_full)


lm_4 <- lmer(log(ZOI_area) ~ 
               Treatment + 
               poly(Temperature_centered,2)+
               poly(Age_centered,2)+
               Temperature_centered*Treatment+
               #Age_centered*Treatment+
               poly(Temperature_centered,2)*poly(Age_centered,2)+
               #Temperature_centered*Age_centered*Treatment+
               # (1|Sample_ID)+
               (1|Plate_ID_Total)+
               (1|Batch_Number),
             data=ZOI_data,
             REML=FALSE)

lm_5 <- lmer(log(ZOI_area) ~ 
               #Treatment + 
               poly(Temperature_centered,2)+
               poly(Age_centered,2)+
               #Temperature_centered*Treatment+
               #Age_centered*Treatment+
               poly(Temperature_centered,2)*poly(Age_centered,2)+
               #Temperature_centered*Age_centered*Treatment+
               # (1|Sample_ID)+
               (1|Plate_ID_Total)+
               (1|Batch_Number),
             data=ZOI_data,
             REML=FALSE)

summary(lm_4) #can possibly drop a random effect
drop1(lm_4,test="Chisq") #cannot simplify more
rand(lm_4) #leave remaining random effects in model
anova(lm_2,lm_3,lm_4,lm_full)

anova(lm_4,lm_5)

plot(lm_4)
#go with lm_4

#re-fit lm4 with REML=TRUE for more accurate estimates of variance:
lm_4 <- lmer(log(ZOI_area) ~ 
               Treatment + 
               poly(Temperature_centered,2)+
               poly(Age_centered,2)+
               Temperature_centered*Treatment+
               #Age_centered*Treatment+
               poly(Temperature_centered,2)*poly(Age_centered,2)+
               #Temperature_centered*Age_centered*Treatment+
               # (1|Sample_ID)+
               (1|Plate_ID_Total)+
               (1|Batch_Number),
             data=ZOI_data,
             REML=TRUE)

library(lmerTest)
#Anova(lm_4,type=2,test.statistic = "Chisq") 


#interaction of temperature and treatment seems strong-
#what if we separate treatments to get at the effects of temp and age on each treatment?


##################
#save outputs:
sink("fullmodel_summary.txt")
print(summary(lm_4))
sink()  # returns output to the console

#Anova(lm_4,type=2,test.statistic = "F")

library(lmerTest)

sink("fullmodel_ANOVAoutput.txt")
print(anova(lm_4,type=2,ddf ="Kenward-Roger"))
sink()

full_anova <- parameters::model_parameters(anova(lm_4,type=2,ddf ="Kenward-Roger"))
full_anova<-as.data.frame(full_anova)
full_anova

write_xlsx(full_anova,"fullmodel_ANOVA.xlsx")

modelvalues_fullmodel <- parameters::model_parameters(lm_4, exponentiate = TRUE)
head(modelvalues_fullmodel)

write_xlsx(modelvalues_fullmodel, "fullmodel_expcoefandCIs.xlsx")

#library(effectsize)
#fullmodel_effectsize_partial <- eta_squared(lm_4,partial=TRUE)

#eta_squared(lm_4)
#epsilon_squared(lm_4)

#eta_squared(car::Anova(lm_4, type = 2), partial = TRUE)

#fullmodel_effectsize_partial_dataframe <- as.data.frame(fullmodel_effectsize_partial)
#fullmodel_effectsize_partial_dataframe
#write_xlsx(fullmodel_effectsize_partial_dataframe,"fullmodel_partialeffectsizes.xlsx")

# Effect size based on chi-squared value
#install.packages("esc")
#library(esc)
#esc_chisq(chisq = 136.2, totaln = 303)


#save outputs:

summary(lm_4)
ref_grid(lm_4)

emmip(lm_4,Age_centered~Temperature_centered | Treatment, cov.reduce = range)


full_emm <- emmeans(lm_4,specs=c("Temperature_centered", "Age_centered","Treatment"),
                  type="response",
                  at = list(Temperature_centered = c(-1.2,1.19,0.23),
                            Age_centered = c(-1.25,-0.47,0.5,1.48)))
full_emm

pairs(full_emm, adjust="tukey")

full_dv = cld(full_emm,
            alpha   = 0.05,
            Letters = letters,         ###  Use lowercase letters for .group
            adjust  = "sidak")         ###  sidak-adjusted comparisons
full_dv
str(full_dv)

full_posthoc_comparisons_dv <- as.data.frame(full_dv)
write_xlsx(full_posthoc_comparisons_dv,"posthoc_comparisons_fullmodel.xlsx")

str(full_posthoc_comparisons_dv)
#make as factor to graph:

full_posthoc_comparisons_dv$Temperature_centered <- as.factor(full_posthoc_comparisons_dv$Temperature_centered)
full_posthoc_comparisons_dv$Age_centered <- as.factor(full_posthoc_comparisons_dv$Age_centered)

head(full_posthoc_comparisons_dv)
levels(full_posthoc_comparisons_dv$Age_centered)

#make labels:
full_posthoc_comparisons_dv$Age_centered <- factor(full_posthoc_comparisons_dv$Age_centered,
                                                 labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

full_posthoc_comparisons_dv$Temperature_centered <- factor(full_posthoc_comparisons_dv$Temperature_centered,
                                                         labels = c("27","30","32"))


full_posthoc_comparisons_dv$Treatment <- factor(full_posthoc_comparisons_dv$Treatment,    # Change factor labels
                                             labels = c("Naïve","Injury",
                                                        "italic(`E. coli`)",
                                                        "italic(`M. luteus`)"))
full_emmeans1<- full_posthoc_comparisons_dv %>%
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(fill=Treatment),
           stat="identity",
           width=0.8)+
  facet_grid(Treatment~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
full_emmeans1

ggsave("full_emmeans1_png.png",plot=full_emmeans1,width = 6, height = 5, units = "in", dpi = 600)
ggsave("full_emmeans1_eps.eps",plot=full_emmeans1,width = 6, height = 5, units = "in", dpi = 600)
ggsave("full_emmeans1_pdf.pdf",plot=full_emmeans1,width = 6, height = 5, units = "in", dpi = 600)


full_posthoc_comparisons_dv$Age_centered <- factor(full_posthoc_comparisons_dv$Age_centered,
                                                 labels = c("1","5","10","15"))

full_posthoc_comparisons_dv$Temperature_centered <- factor(full_posthoc_comparisons_dv$Temperature_centered,
                                                         labels = c("`27˚C`","`30˚C`","`32˚C`"))

full_posthoc_comparisons_dv$Treatment <- factor(full_posthoc_comparisons_dv$Treatment,    # Change factor labels
                                                labels = c("Naïve","Injury",
                                                           "italic(`E. coli`)",
                                                           "italic(`M. luteus`)"))

full_emmeans2 <- full_posthoc_comparisons_dv %>%
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(fill=Treatment),
           stat="identity",
           width=0.8)+
  facet_grid(Treatment~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+ 
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
full_emmeans2


ggsave("full_emmeans2_png.png",plot=full_emmeans2,width = 6, height = 5, units = "in", dpi = 600)
ggsave("full_emmeans2_eps.eps",plot=full_emmeans2,width = 6, height = 5, units = "in", dpi = 600)
ggsave("full_emmeans2_pdf.pdf",plot=full_emmeans2,width = 6, height = 5, units = "in", dpi = 600)
ggsave("full_emmeans2_jpg.jpg",plot=full_emmeans2,width = 6, height = 5, units = "in", dpi = 600)

#mean by temp or age or treatment only


full_Ecemm <- as.data.frame(full_emm)
head(full_Ecemm)

full_Ecemm_TEMPonly <- full_Ecemm %>%
  group_by(Temperature_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(full_Ecemm_TEMPonly) #gives average ZOI for each temp across all ages

full_Ecemm_AGEonly <- full_Ecemm %>%
  group_by(Age_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(full_Ecemm_AGEonly) #gives average ZOI for each age across all temps

full_Ecemm_TREATMENTonly <- full_Ecemm %>%
  group_by(Treatment) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(full_Ecemm_TREATMENTonly) 

full_Ecemm_TEMPonly <- as.data.frame(full_Ecemm_TEMPonly)
full_Ecemm_AGEonly <- as.data.frame(full_Ecemm_AGEonly)
full_Ecemm_TREATMENTonly <- as.data.frame(full_Ecemm_TREATMENTonly)

write_xlsx(full_Ecemm_TEMPonly, "full_Ecemm_TEMPonly.xlsx")
write_xlsx(full_Ecemm_AGEonly, "full_Ecemm_AGEonly.xlsx")
write_xlsx(full_Ecemm_TREATMENTonly, "full_Ecemm_TREATMENTonly.xlsx")


full_Ecemm_TEMPonly$Temperature_centered <- as.factor(full_Ecemm_TEMPonly$Temperature_centered)
full_Ecemm_AGEonly$Age_centered <- as.factor(full_Ecemm_AGEonly$Age_centered)

full_Ecemm_AGEonly$Age_centered <- factor(full_Ecemm_AGEonly$Age_centered,
                                        labels = c("1","5","10","15"))

full_Ecemm_TEMPonly$Temperature_centered <- factor(full_Ecemm_TEMPonly$Temperature_centered,
                                                 labels = c("27","30","32"))

full_Temponly <- full_Ecemm_TEMPonly %>%
  ggplot(aes(x=Temperature_centered,y=mean))+
  geom_bar(aes(fill=Temperature_centered),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
full_Temponly

ggsave("full_Temponly_png.png",plot=full_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("full_Temponly_eps.eps",plot=full_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("full_Temponly_pdf.pdf",plot=full_Temponly,width = 4, height = 4, units = "in", dpi = 600)

full_Ageonly <- full_Ecemm_AGEonly %>%
  ggplot(aes(x=Age_centered,y=mean))+
  geom_bar(aes(fill= Age_centered),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+   scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
full_Ageonly

ggsave("full_Ageonly_png.png",plot=full_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("full_Ageonly_eps.eps",plot=full_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("full_Ageonly_pdf.pdf",plot=full_Ageonly,width = 4, height = 4, units = "in", dpi = 600)


full_TREATMENTonly <- full_Ecemm_TREATMENTonly %>%
  ggplot(aes(x=Treatment,y=mean))+
  geom_bar(aes(fill= Treatment),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Immune Treatment") +
  scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","Injury",bquote(italic("E. coli")), bquote(italic("M. luteus"))))+
  theme_pubr()+ 
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
full_TREATMENTonly

ggsave("full_TREATMENTonly_png.png",plot=full_TREATMENTonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("full_TREATMENTonly_eps.eps",plot=full_TREATMENTonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("full_TREATMENTonly_pdf.pdf",plot=full_TREATMENTonly,width = 4, height = 4, units = "in", dpi = 600)

####
#compare naive/injury vs. ecoli/mlut for full model:
treat_comparison<-pairs(full_emm,simple=c("Treatment"),combine=TRUE,adjust="Tukey")
treat_comparison <- as.data.frame(treat_comparison)
treat_comparison
write_xlsx(treat_comparison,"Treatmenteffects_fullmodelemm.xlsx")

#contrast(full_emm,"eff",by="Treatment",enhance.levels="Temperature_centered")

#contrast(warp.emm,"eff",by= NULL, enhance.levels=c("wool", "tension"))

full_emm_treatment <- emmeans(lm_4,specs=c("Treatment"),
                    type="response")
                    #at = list(Temperature_centered = c(-1.2,1.19,0.23),
                     #         Age_centered = c(-1.25,-0.47,0.5,1.48)))
full_emm_treatment

#emmeans(lm_4,specs=pairwise ~"Treatment",type="response")

full_dv_treat = cld(full_emm_treatment,
              alpha   = 0.05,
              Letters = letters,         ###  Use lowercase letters for .group
              adjust  = "sidak")         ###  sidak-adjusted comparisons
full_dv_treat
str(full_dv_treat)

#full_posthoc_comparisons_dv <- as.data.frame(full_dv)
#write_xlsx(full_posthoc_comparisons_dv,"posthoc_comparisons_fullmodel.xlsx")


######################
#treatment has a major effect.

#stratify categorical by treatment:
#ZOI_data <- subset(ZOI_data,Technical_Rep == 1)

#ZOI_data$Temperature_centered <- ZOI_data$Temperature
#ZOI_data$Age_centered <- ZOI_data$Age

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


#refit w/ REML = TRUE
lm_Naive_3 <-  lmer(log(ZOI_area) ~ Temperature_centered+
                      Age_centered + 
                      #Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      #(1|Plate_ID_Total)+
                      (1|Batch_Number),
                    data=Naive_ZOI,
                    REML=TRUE)
lm_Naive_3_factor <-  lmer(log(ZOI_area) ~ factor(Temperature_centered)+
                      factor(Age_centered) + 
                      #Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      #(1|Plate_ID_Total)+
                      (1|Batch_Number),
                    data=Naive_ZOI,
                    REML=TRUE)

library(lmerTest)

#save outputs:
sink("Naive_summary.txt")
print(summary(lm_Naive_3))
sink()  # returns output to the console

sink("Naive_ANOVAoutput.txt")
print(anova(lm_Naive_3,type="2",ddf ="Kenward-Roger"))
sink()

naive_anova <- parameters::model_parameters(anova(lm_Naive_3,type=2,ddf ="Kenward-Roger"))
naive_anova<-as.data.frame(naive_anova)
naive_anova
#naive_anova <- as.data.frame(Anova(lm_Naive_3,type="2",test.statistic = "Chisq"))
head(naive_anova)
write_xlsx(naive_anova, "Naive_ANOVA.xlsx")

modelvalues_naive <- parameters::model_parameters(lm_Naive_3, exponentiate = TRUE)
head(modelvalues_naive)

write_xlsx(modelvalues_naive, "Naive_model_expcoefandCIs.xlsx")

library(effectsize)

Naive_effectsize_partial <- eta_squared(lm_Naive_3,partial=TRUE)

#eta_squared(car::Anova(lm_Naive_3, type = 2), partial = FALSE)
Naive_effectsize_partial
Naive_effectsize_partial_dataframe <- as.data.frame(Naive_effectsize_partial)

write_xlsx(Naive_effectsize_partial_dataframe,"Naive_partialeffectsizes.xlsx")

#takeaway: temp affects naive, but not age or interaction


#save outputs:
lm_Naive_3$Temperature_centered <- factor(lm_Naive_3$Temperature_centered,
                                                        labels = c("27","30","32"))

plot.data <- emmip(lm_Naive_3,Age_centered~Temperature_centered,style="factor",col="black",
                   linearg = list(), dotarg = list(size = 3),
                   cov.reduce=FALSE,type="response",CIs = FALSE, plotit = FALSE)

str(lm_Naive_3)
summary(lm_Naive_3)
require("ggplot2")


Naive_emmip <- emmip(lm_Naive_3,Age_centered~Temperature_centered,style="factor",col="black",
      linearg = list(), dotarg = list(size = 3),
      cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(1,5,10,15))+
  scale_linetype(labels=c(1,5,10,15))+
  guides(shape = guide_legend(title = "Adult Age (days)"),
         linetype = guide_legend(title = "Adult Age (days)"))+
  scale_x_discrete(breaks = c(-1.2,0.23,1.19),labels=c(27,30,32))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +  theme(legend.position = "right")+
  #scale_y_continuous(limits=c(5,12))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))
Naive_emmip

ggsave("Naive_interaction.png",plot=Naive_emmip,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Naive_interaction.eps",plot=Naive_emmip,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Naive_interaction.pdf",plot=Naive_emmip,width = 5.5, height = 4, units = "in", dpi = 600)


Naive_emm <- emmeans(lm_Naive_3,specs=c("Temperature_centered","Age_centered"),
                  type="response",
                  at = list(Temperature_centered = c(-1.2,1.19,0.23),
                            Age_centered = c(-1.25,-0.47,0.5,1.48)))

Naive_emm <- emmeans(lm_Naive_3,~Temperature_centered|Age_centered,
                     type="response",
                     at = list(Temperature_centered = c(-1.2,1.19,0.23),
                               Age_centered = c(-1.25,-0.47,0.5,1.48)))
Naive_emm
plot(Naive_emm)
head(Naive_emm)

Naive_emmeans_data <- as.data.frame(Naive_emm)
Naive_emmeans_data
write_xlsx(Naive_emmeans_data,"Naive_emmeans.xlsx")

Naive_emmeans_pwc<-pairs(Naive_emm,type="response", adjust="sidak")
Naive_emmeans_pwc <- as.data.frame(Naive_emmeans_pwc)
head(Naive_emmeans_pwc)
write_xlsx(Naive_emmeans_pwc, "posthoc_comparisons_Naive_pairwise.xlsx")

Naive_dv = cld(Naive_emm,
            alpha   = 0.05,
            Letters = letters,         ###  Use lowercase letters for .group
            adjust  = "sidak")         ###  sidak-adjusted comparisons
Naive_dv
str(Naive_dv)

Naive_posthoc_comparisons_dv <- as.data.frame(Naive_dv)
write_xlsx(Naive_posthoc_comparisons_dv,"posthoc_comparisons_Naive.xlsx")

#make as factor to graph:
Naive_posthoc_comparisons_dv$Temperature_centered <- as.factor(Naive_posthoc_comparisons_dv$Temperature_centered)
Naive_posthoc_comparisons_dv$Age_centered <- as.factor(Naive_posthoc_comparisons_dv$Age_centered)

head(Naive_posthoc_comparisons_dv)
levels(Naive_posthoc_comparisons_dv$Age_centered)

#make labels:
Naive_posthoc_comparisons_dv$Age_centered <- factor(Naive_posthoc_comparisons_dv$Age_centered,
                                                 labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

Naive_posthoc_comparisons_dv$Temperature_centered <- factor(Naive_posthoc_comparisons_dv$Temperature_centered,
                                                         labels = c("27","30","32"))

Naive_emmeans1<- Naive_posthoc_comparisons_dv %>%
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(fill="#F8766D"),
           stat="identity",
           width=0.8)+
  facet_grid(~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,45))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Naive_emmeans1

ggsave("Naive_emmeans1_png.png",plot=Naive_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_emmeans1_eps.eps",plot=Naive_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_emmeans1_pdf.pdf",plot=Naive_emmeans1,width = 6, height = 4, units = "in", dpi = 600)


Naive_posthoc_comparisons_dv$Age_centered <- factor(Naive_posthoc_comparisons_dv$Age_centered,
                                                 labels = c("1","5","10","15"))

Naive_posthoc_comparisons_dv$Temperature_centered <- factor(Naive_posthoc_comparisons_dv$Temperature_centered,
                                                         labels = c("`27˚C`","`30˚C`","`32˚C`"))

Naive_emmeans2 <- Naive_posthoc_comparisons_dv %>%
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(fill="#F8766D"),
           stat="identity",
           width=0.8)+
  facet_grid(~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+ 
  #scale_fill_npg()+
  scale_fill_manual(values="#F8766D")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,45))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Naive_emmeans2


ggsave("Naive_emmeans2_png.png",plot=Naive_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_emmeans2_eps.eps",plot=Naive_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_emmeans2_pdf.pdf",plot=Naive_emmeans2,width = 6, height = 4, units = "in", dpi = 600)

#mean by temp or age only

ref_grid(lm_Naive_3)
Naive_emm_temp <- emmeans(lm_Naive_3,specs=c("Temperature_centered"),
                     type="response",
                     at = list(Temperature_centered = c(-1.2,1.19,0.23)))


pairs(Naive_emm_temp, adjust="sidak")

Naive_dv_temp = cld(Naive_emm_temp,
               alpha   = 0.05,
               Letters = letters,         ###  Use lowercase letters for .group
               adjust  = "sidak")         ###  sidak-adjusted comparisons
Naive_dv_temp

Naive_posthoc_comparisons_dv_temp <- as.data.frame(Naive_dv_temp)
write_xlsx(Naive_posthoc_comparisons_dv_temp,"posthoc_comparisons_Naive_temp.xlsx")

Naive_emm_age <- emmeans(lm_Naive_3,specs=c("Age_centered"),
                     type="response",
                     at = list(Age_centered = c(-1.25,-0.47,0.5,1.48)))


pairs(Naive_emm_age, adjust="sidak")

Naive_dv_age = cld(Naive_emm_age,
                    alpha   = 0.05,
                    Letters = letters,         ###  Use lowercase letters for .group
                    adjust  = "sidak")         ###  sidak-adjusted comparisons
Naive_dv_age

Naive_posthoc_comparisons_dv_age <- as.data.frame(Naive_dv_age)
Naive_posthoc_comparisons_dv_age
write_xlsx(Naive_posthoc_comparisons_dv_age,"posthoc_comparisons_Naive_age.xlsx")

#emtrends(lm_Naive_3,pairwise~as.factor(Temperature_centered),var="Age_centered")
contrast(Naive_emm, "pairwise", 
         at = list(Temperature_centered = c(-1.2,1.19,0.23), Age_centered = c(-1.25,-0.47,0.5,1.48)),
         type="response",
         simple = "each", combine = TRUE, adjust = "mvt")
#mvcontrast(Naive_emm, "pairwise", mult.name = c("Temperature_centered", "Age_centered"),type="response")
#contrast(noise.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

emtrends(lm_Naive_3, pairwise~factor(Temperature_centered),
         var = "Age_centered",
         at=list(Age_centered = c(-1.25,-0.47,0.5,1.48)),
         adjust = "mvt")
emtrends(
  contcont,
  pairwise ~ effort,
  var = "hours",
  at = mylist,
  adjust = "none"
)


Naive_emm[1]

mvcontrast(Naive_emm, "pairwise", type="response", adjust="sidak",mult.name = c("Age_centered"),show.ests = TRUE)
mvcontrast(Naive_emm, "pairwise", type="response", adjust="sidak",mult.name = c("Temperature_centered"),show.ests = TRUE)


contrast(emm_s.t[[1]], "poly")   ## 'by = "type"' already in previous result 

########
#plot:
Naive_Ecemm <- as.data.frame(Naive_emm)
head(Naive_Ecemm)

Naive_Ecemm_TEMPonly <- Naive_Ecemm %>%
  group_by(Temperature_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(Naive_Ecemm_TEMPonly) #gives average ZOI for each temp across all ages

Naive_Ecemm_AGEonly <- Naive_Ecemm %>%
  group_by(Age_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(Naive_Ecemm_AGEonly) #gives average ZOI for each age across all temps

Naive_Ecemm_TEMPonly <- as.data.frame(Naive_Ecemm_TEMPonly)
Naive_Ecemm_AGEonly <- as.data.frame(Naive_Ecemm_AGEonly)

write_xlsx(Naive_Ecemm_TEMPonly, "Naive_Ecemm_TEMPonly.xlsx")
write_xlsx(Naive_Ecemm_AGEonly, "Naive_Ecemm_AGEonly.xlsx")

Naive_Ecemm_TEMPonly$Temperature_centered <- as.factor(Naive_Ecemm_TEMPonly$Temperature_centered)
Naive_Ecemm_AGEonly$Age_centered <- as.factor(Naive_Ecemm_AGEonly$Age_centered)

Naive_Ecemm_AGEonly$Age_centered <- factor(Naive_Ecemm_AGEonly$Age_centered,
                                        labels = c("1","5","10","15"))

Naive_Ecemm_TEMPonly$Temperature_centered <- factor(Naive_Ecemm_TEMPonly$Temperature_centered,
                                                 labels = c("27","30","32"))


Naive_Temponly <- Naive_Ecemm_TEMPonly %>%
  ggplot(aes(x=Temperature_centered,y=mean))+
  geom_bar(aes(fill=c("#4D6FAE","#6F9F51", "#CC763B")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  scale_fill_manual(values=#"#F8766D",
                    c("#4D6FAE","#6F9F51", "#CC763B"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,45))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Naive_Temponly

ggsave("Naive_Temponly_png.png",plot=Naive_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Naive_Temponly_eps.eps",plot=Naive_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Naive_Temponly_pdf.pdf",plot=Naive_Temponly,width = 4, height = 4, units = "in", dpi = 600)

Naive_Ageonly <- Naive_Ecemm_AGEonly %>%
  ggplot(aes(x=Age_centered,y=mean))+
  geom_bar(aes(fill=c("#DCD1E9","#BAA4D3","#9776BE","#7549A8")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+  
  scale_fill_manual(values= c("#DCD1E9","#BAA4D3","#9776BE","#7549A8"))+
                      #"#F8766D"
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,45))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Naive_Ageonly

ggsave("Naive_Ageonly_png.png",plot=Naive_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Naive_Ageonly_eps.eps",plot=Naive_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Naive_Ageonly_pdf.pdf",plot=Naive_Ageonly,width = 4, height = 4, units = "in", dpi = 600)



######################################################################
######################################################################
###############################
#LB:
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
                data=LB_ZOI,REML=FALSE)
summary(lm_LB_2)
plot(lm_LB_2)
rand(lm_LB_2)
Anova(lm_LB_2,type=2)
drop1(lm_LB_2)

nullLB <-lm(log(ZOI_area)~1,data=LB_ZOI)
anova(lm_LB_2,nullLB)
AIC(lm_LB_2)

#no effects on injury group


lm_LB_4 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                poly(Age_centered,2)+
                  (1|Seeded_plates_OD),REML=FALSE,
              data=LB_ZOI)

plot(lm_LB_4)
Anova(lm_LB_2,type=2)
AIC(lm_LB_2,lm_LB_4,nullLB) #2 is better

#try glm:
lm_LB5 <- glm(ZOI_area ~ Temperature_centered + Age_centered,
              data=LB_ZOI, family=gaussian(link="log"))

summary(lm_LB5)
plot(lm_LB5) #not great
AIC(lm_LB_2,lm_LB5) 
logLik(lm_LB5)
logLik(lm_LB_2)
#2 is best

#re-run as REML=TRUE
lm_LB_2 <- lmer(log(ZOI_area) ~ Temperature_centered +
                  Age_centered+
                  # Temperature_centered*Age_centered+
                  #(1|Sample_ID)+
                  #(1|Plate_ID_Total)+
                  (1|Seeded_plates_OD),
                data=LB_ZOI,REML=TRUE)

##
library(lmerTest)
#save outputs:
sink("LB_summary.txt")
print(summary(lm_LB_2))
sink()  # returns output to the console

sink("LB_ANOVAoutput.txt")
print(anova(lm_LB_2,type="2",ddf="Kenward-Roger"))
sink()

LB_anova <- parameters::model_parameters(anova(lm_LB_2,type="2",ddf="Kenward-Roger"))
LB_anova<-as.data.frame(LB_anova)

write_xlsx(LB_anova,"LB_ANOVA.xlsx")

modelvalues_LB <- parameters::model_parameters(lm_LB_2, exponentiate = TRUE)
head(modelvalues_LB)

write_xlsx(modelvalues_LB, "LB_model_expcoefandCIs.xlsx")

LB_effectsize_partial <- eta_squared(lm_LB_2)
LB_effectsize_partial_dataframe <- as.data.frame(LB_effectsize_partial)
LB_effectsize_partial_dataframe
write_xlsx(LB_effectsize_partial_dataframe,"LB_partialeffectsizes.xlsx")


#save outputs:


ref_grid(lm_LB_2)

LB_emmip<-emmip(lm_LB_2,Age_centered~Temperature_centered,style="factor",col="black",
      linearg = list(), dotarg = list(size = 3),
      cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(1,5,10,15))+
  scale_linetype(labels=c(1,5,10,15))+
  guides(shape = guide_legend(title = "Adult Age (days)"),
         linetype = guide_legend(title = "Adult Age (days)"))+
  scale_x_discrete(breaks = c(-1.2,0.23,1.19),labels=c(27,30,32))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +  theme(legend.position = "right")+
  scale_y_continuous(limits=c(5,12))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))
LB_emmip

ggsave("LB_interaction.png",plot=LB_emmip,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("LB_interaction.eps",plot=LB_emmip,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("LB_interaction.pdf",plot=LB_emmip,width = 5.5, height = 4, units = "in", dpi = 600)


LB_emm <- emmeans(lm_LB_2,specs=c("Temperature_centered", "Age_centered"),
                  type="response",
                  at = list(Temperature_centered = c(-1.2,1.19,0.23),
                            Age_centered = c(-1.25,-0.47,0.5,1.48)))
LB_emm

LB_emm_df <- as.data.frame(LB_emm)
write_xlsx(LB_emm_df,"LB_emmeans.xlsx")

LB_emmeans_pwc <- pairs(LB_emm, adjust="sidak")
LB_emmeans_pwc <- as.data.frame(LB_emmeans_pwc)
write_xlsx(LB_emmeans_pwc, "posthoc_comparisons_LB_pairwise.xlsx")

LB_dv = cld(LB_emm,
             alpha   = 0.05,
             Letters = letters,         ###  Use lowercase letters for .group
             adjust  = "sidak")         ###  sidak-adjusted comparisons
LB_dv
str(LB_dv)

LB_posthoc_comparisons_dv <- as.data.frame(LB_dv)
write_xlsx(LB_posthoc_comparisons_dv,"posthoc_comparisons_LB.xlsx")

#make as factor to graph:
LB_posthoc_comparisons_dv$Temperature_centered <- as.factor(LB_posthoc_comparisons_dv$Temperature_centered)
LB_posthoc_comparisons_dv$Age_centered <- as.factor(LB_posthoc_comparisons_dv$Age_centered)

head(LB_posthoc_comparisons_dv)
levels(LB_posthoc_comparisons_dv$Age_centered)

#make labels:
LB_posthoc_comparisons_dv$Age_centered <- factor(LB_posthoc_comparisons_dv$Age_centered,
                                                    labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

LB_posthoc_comparisons_dv$Temperature_centered <- factor(LB_posthoc_comparisons_dv$Temperature_centered,
                                                            labels = c("27","30","32"))

LB_emmeans1<- LB_posthoc_comparisons_dv %>%
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(fill="#7CAE00"),
           stat="identity",
           width=0.8)+
  facet_grid(~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  scale_fill_manual(values="#7CAE00")+
  theme_pubr()+ 
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
LB_emmeans1

ggsave("LB_emmeans1_png.png",plot=LB_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("LB_emmeans1_eps.eps",plot=LB_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("LB_emmeans1_pdf.pdf",plot=LB_emmeans1,width = 6, height = 4, units = "in", dpi = 600)


LB_posthoc_comparisons_dv$Age_centered <- factor(LB_posthoc_comparisons_dv$Age_centered,
                                                    labels = c("1","5","10","15"))

LB_posthoc_comparisons_dv$Temperature_centered <- factor(LB_posthoc_comparisons_dv$Temperature_centered,
                                                            labels = c("`27˚C`","`30˚C`","`32˚C`"))

LB_emmeans2 <- LB_posthoc_comparisons_dv %>%
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(fill="#7CAE00"),
           stat="identity",
           width=0.8)+
  facet_grid(~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+ 
  scale_fill_manual(values="#7CAE00")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
LB_emmeans2


ggsave("LB_emmeans2_png.png",plot=LB_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("LB_emmeans2_eps.eps",plot=LB_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("LB_emmeans2_pdf.pdf",plot=LB_emmeans2,width = 6, height = 4, units = "in", dpi = 600)

#mean by temp or age only


LB_Ecemm <- as.data.frame(LB_emm)
head(LB_Ecemm)

LB_Ecemm_TEMPonly <- LB_Ecemm %>%
  group_by(Temperature_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(LB_Ecemm_TEMPonly) #gives average ZOI for each temp across all ages

LB_Ecemm_AGEonly <- LB_Ecemm %>%
  group_by(Age_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(LB_Ecemm_AGEonly) #gives average ZOI for each age across all temps

LB_Ecemm_TEMPonly <- as.data.frame(LB_Ecemm_TEMPonly)
LB_Ecemm_AGEonly <- as.data.frame(LB_Ecemm_AGEonly)

write_xlsx(LB_Ecemm_TEMPonly, "LB_Ecemm_TEMPonly.xlsx")
write_xlsx(LB_Ecemm_AGEonly, "LB_Ecemm_AGEonly.xlsx")

LB_Ecemm_TEMPonly$Temperature_centered <- as.factor(LB_Ecemm_TEMPonly$Temperature_centered)
LB_Ecemm_AGEonly$Age_centered <- as.factor(LB_Ecemm_AGEonly$Age_centered)

LB_Ecemm_AGEonly$Age_centered <- factor(LB_Ecemm_AGEonly$Age_centered,
                                           labels = c("1","5","10","15"))

LB_Ecemm_TEMPonly$Temperature_centered <- factor(LB_Ecemm_TEMPonly$Temperature_centered,
                                                    labels = c("27","30","32"))


LB_Temponly <- LB_Ecemm_TEMPonly %>%
  ggplot(aes(x=Temperature_centered,y=mean))+
  geom_bar(aes(fill= c("#4D6FAE","#6F9F51", "#CC763B")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  #scale_fill_npg()+
  scale_fill_manual(values=c("#4D6FAE","#6F9F51", "#CC763B"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
LB_Temponly

ggsave("LB_Temponly_png.png",plot=LB_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("LB_Temponly_eps.eps",plot=LB_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("LB_Temponly_pdf.pdf",plot=LB_Temponly,width = 4, height = 4, units = "in", dpi = 600)

LB_Ageonly <- LB_Ecemm_AGEonly %>%
  ggplot(aes(x=Age_centered,y=mean))+
  geom_bar(aes(fill= c("#DCD1E9","#BAA4D3","#9776BE","#7549A8")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+   
  #scale_fill_npg()+
  scale_fill_manual(values=c("#DCD1E9","#BAA4D3","#9776BE","#7549A8"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
LB_Ageonly

ggsave("LB_Ageonly_png.png",plot=LB_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("LB_Ageonly_eps.eps",plot=LB_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("LB_Ageonly_pdf.pdf",plot=LB_Ageonly,width = 4, height = 4, units = "in", dpi = 600)




ref_grid(lm_LB_2)
LB_emm_temp <- emmeans(lm_LB_2,specs=c("Temperature_centered"),
                          type="response",
                          at = list(Temperature_centered = c(-1.2,1.19,0.23)))


pairs(LB_emm_temp, adjust="sidak")

LB_emm_temp = cld(LB_emm_temp,
                    alpha   = 0.05,
                    Letters = letters,         ###  Use lowercase letters for .group
                    adjust  = "sidak")         ###  sidak-adjusted comparisons
LB_emm_temp

LB_posthoc_comparisons_temp <- as.data.frame(LB_emm_temp)
write_xlsx(LB_posthoc_comparisons_temp,"posthoc_comparisons_LB_temp.xlsx")

LB_emm_age <- emmeans(lm_LB_2,specs=c("Age_centered"),
                         type="response",
                         at = list(Age_centered = c(-1.25,-0.47,0.5,1.48)))


pairs(LB_emm_age, adjust="sidak")

LB_dv_age = cld(LB_emm_age,
                   alpha   = 0.05,
                   Letters = letters,         ###  Use lowercase letters for .group
                   adjust  = "sidak")         ###  sidak-adjusted comparisons
LB_dv_age

LB_posthoc_comparisons_dv_age <- as.data.frame(LB_dv_age)
LB_posthoc_comparisons_dv_age
write_xlsx(LB_posthoc_comparisons_dv_age,"posthoc_comparisons_LB_age.xlsx")


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
plot(lmer_Ecoli_1) #bad residuals
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
#need to account for nonlinear relationship:
lmer_Ecoli_3 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       Age_centered + 
                       Temperature_centered*Age_centered+
                       (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                     (1|Batch_Number),
                     data=Ecoli_ZOI)
summary(lmer_Ecoli_3)
drop1(lmer_Ecoli_3)
rand(lmer_Ecoli_3)

Anova(lmer_Ecoli_3)
plot(lmer_Ecoli_3)


lmer_Ecoli_4 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       poly(Age_centered,2) + 
                       poly(Temperature_centered,2)*poly(Age_centered,2)+
                       (1|Sample_ID) + 
                       (1|Batch_Number)+
                       (1|Plate_ID_Total),
                     data=Ecoli_ZOI,
                     REML=FALSE)

summary(lmer_Ecoli_4)
AIC(lmer_Ecoli_4)
plot(lmer_Ecoli_4)
rand(lmer_Ecoli_4) #can get rid of plate id total
drop1(lmer_Ecoli_4) #interaction is significant
Anova(lmer_Ecoli_4)

#lmer_Ecoli_g4 <- glmer((ZOI_area) ~  poly(Temperature_centered,2) +
 #                      poly(Age_centered,2) + 
  #                     poly(Temperature_centered,2)*poly(Age_centered,2)+
   #                    (1|Sample_ID) + 
    #                   (1|Batch_Number),
     #                  family=gaussian(link="log"), 
      #                                  control=glmerControl(optimizer="bobyqa",
       #                                optCtrl=list(maxfun=1e5)),
        #             data=Ecoli_ZOI)

lmer_Ecoli_4b <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       poly(Age_centered,2) + 
                       Temperature_centered*Age_centered+
                       poly(Temperature_centered,2)*poly(Age_centered,2)+
                       (1|Sample_ID) + 
                       (1|Batch_Number)+
                       (1|Plate_ID_Total),
                     data=Ecoli_ZOI,
                     REML=FALSE)
summary(lmer_Ecoli_4b)
AIC(lmer_Ecoli_4b)
plot(lmer_Ecoli_4b)


lmer_Ecoli_4c <- lm(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       poly(Age_centered,2) + 
                       poly(Temperature_centered,2)*poly(Age_centered,2),
                     data=Ecoli_ZOI)

summary(lmer_Ecoli_4c)
AIC(lmer_Ecoli_4c)
plot(lmer_Ecoli_4c)

AIC(lmer_Ecoli_4,lmer_Ecoli_4b,lmer_Ecoli_4c)

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
                     (1|Sample_ID) + 
                     (1|Batch_Number),
                   #(1|Plate_ID_Total),
                   data=Ecoli_ZOI)

summary(lmer_Ecoli_6)
plot(lmer_Ecoli_6)
rand(lmer_Ecoli_6)
Anova(lmer_Ecoli_6,type=2)

lmer_Ecoli_7 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       poly(Age_centered,2)+
                      Temperature_centered*Age_centered+
                       (1|Sample_ID) + 
                       (1|Batch_Number),
                     #(1|Plate_ID_Total),
                     data=Ecoli_ZOI)
summary(lmer_Ecoli_7)
plot(lmer_Ecoli_7)
rand(lmer_Ecoli_7)
Anova(lmer_Ecoli_7,type=2)

lmer_Ecoli_8 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                   poly(Age_centered,3)+
                     (1|Sample_ID) + (1|Batch_Number),
                   #(1|Plate_ID),
                   #(1|Batch_Number),
                   data=Ecoli_ZOI)
plot(lmer_Ecoli_8)

AIC(lmer_Ecoli_3,lmer_Ecoli_4,lmer_Ecoli_5,lmer_Ecoli_6,lmer_Ecoli_4b,lmer_Ecoli_4c,lmer_Ecoli_7,lmer_Ecoli_8)
#go with lmer_Ecoli_6

#re-run as REML=TRUE
lmer_Ecoli_6 <- lmer(log(ZOI_area) ~  poly(Temperature_centered,2) +
                       poly(Age_centered,2)+
                       poly(Temperature_centered*Age_centered,2)+
                       (1|Sample_ID) + 
                       (1|Batch_Number),
                     #(1|Plate_ID_Total),
                     data=Ecoli_ZOI,
                     REML=TRUE)
##

library(modelr)
grid <- Ecoli_ZOI %>% 
  data_grid(Temperature_centered, .model = lmer_Ecoli_6) %>% 
  add_predictions(lmer_Ecoli_6)

ggplot(grid, aes(Temperature_centered, pred)) + 
  geom_point()


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
  add_predictions(lmer_Ecoli_6)

ggplot(Ecoli_ZOI, aes(x=Age_centered, 
                      y=log(ZOI_area),group=Age_centered)) + 
  geom_boxplot()+ 
  geom_point(aes(x=Age_centered,y=pred),data=grid_age,colour="red",size=2)

#####
Ecoli_ZOI <- Ecoli_ZOI %>%
  add_residuals(lmer_Ecoli_6)

Ecoli_ZOI %>%
  ggplot(aes(Temperature_centered, resid,color=as.factor(Age_centered)))+
  facet_grid(~Age_centered)+
  geom_ref_line(h=0)+geom_point()
Ecoli_ZOI %>%
  ggplot(aes(Age_centered, resid, color=as.factor(Temperature_centered)))+
  facet_grid(~Temperature_centered)+
  geom_ref_line(h=0)+geom_point()

############################
#save outputs:
library(lmerTest)

#save outputs:
sink("Ecoli_summary.txt")
print(summary(lmer_Ecoli_6))
sink()  # returns output to the console

sink("Ecoli_ANOVAoutput.txt")
#print(Anova(lmer_Ecoli_6,type="2",test.statistic = "Chisq"))
print(anova(lmer_Ecoli_6,type="2",ddf="Kenward-Roger"))
sink()

Ecoli_anova <- parameters::model_parameters(anova(lmer_Ecoli_6,type="2",ddf="Kenward-Roger"))
Ecoli_anova<-as.data.frame(Ecoli_anova)
Ecoli_anova
write_xlsx(Ecoli_anova,"Ecoli_ANOVA.xlsx")

modelvalues_Ec <- parameters::model_parameters(lmer_Ecoli_6, exponentiate = TRUE)
head(modelvalues_Ec)

write_xlsx(modelvalues_Ec, "Ec_model_expcoefandCIs.xlsx")

Ec_effectsize_partial <- eta_squared(lmer_Ecoli_6)
Ec_effectsize_partial_dataframe <- as.data.frame(Ec_effectsize_partial)
Ec_effectsize_partial_dataframe
write_xlsx(Ec_effectsize_partial_dataframe,"Ec_partialeffectsizes.xlsx")

ref_grid(lmer_Ecoli_6)

emmip(lmer_Ecoli_6,poly(Age_centered,2)~poly(Temperature_centered,2),style="factor",
      cov.reduce=FALSE,type="response",CIs = FALSE)+ theme_pubr()+
  ylab(expression("Predicted ZOI Area"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +  #theme(legend.position = "none")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))


Ecoli_emp<- emmip(lmer_Ecoli_6,poly(Age_centered,2)~poly(Temperature_centered,2),
      style="factor",col="black",
      linearg = list(), dotarg = list(size = 3),
      cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(1,5,10,15))+
  scale_linetype(labels=c(1,5,10,15))+
  guides(shape = guide_legend(title = "Adult Age (days)"),
         linetype = guide_legend(title = "Adult Age (days)"))+
  scale_y_continuous(limits=c(0,35),breaks = c(0,10,20,30))+
  scale_x_discrete(breaks = c(-1.2,0.23,1.19),labels=c(27,30,32))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +  theme(legend.position = "right")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))

ggsave("Ecoli_interaction.png",plot=Ecoli_emp,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_interaction.eps",plot=Ecoli_emp,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_interaction.pdf",plot=Ecoli_emp,width = 5.5, height = 4, units = "in", dpi = 600)

Ecoli_emp2<- emmip(lmer_Ecoli_6,poly(Temperature_centered,2)~poly(Age_centered,2),
           style="factor",col="black",
           linearg = list(), dotarg = list(size = 3),
           cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(27,30,32))+
  scale_linetype(labels=c(27,30,32))+
  guides(shape = guide_legend(title = "Temperature (˚C)"),
         linetype = guide_legend(title = "Temperature (˚C)"))+
  scale_y_continuous(limits=c(0,35),breaks = c(0,10,20,30))+
  scale_x_discrete(breaks = c(-1.25,-0.47,0.5,1.48),labels=c(1,5,10,15))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Adult Age (days)") +  theme(legend.position = "right")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))

ggsave("Ecoli_interaction2.png",plot=Ecoli_emp2,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_interaction2.eps",plot=Ecoli_emp2,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_interaction2.pdf",plot=Ecoli_emp2,width = 5.5, height = 4, units = "in", dpi = 600)


str(Ecoli_ZOI$Temperature_centered)
str(Ecoli_ZOI$Age_centered)
Ec_emm <- emmeans(lmer_Ecoli_6,specs=c("Temperature_centered", "Age_centered"),
                          type="response",
                  at = list(Temperature_centered = c(-1.2,1.19,0.23),
                           Age_centered = c(-1.25,-0.47,0.5,1.48)))
Ec_emm

Ec_emm_df <- as.data.frame(Ec_emm)
write_xlsx(Ec_emm_df,"EC_emmeans.xlsx")

#Ec_emm_temp <- emmeans(lmer_Ecoli_6,specs=c("Temperature_centered","Age_centered"),
 #                 type="response",
  #                at = list(Temperature_centered = c(-1.2,1.19,0.23)))

#head(Ec_emm_temp) #gives zoi at each temp for mean centered age

#mean(Ecoli_ZOI$Temperature_centered)
#find emmeans for post hoc comparisons:
library(multcompView)
library(emmeans)
library(multcomp)


Ec_emmeans_pwc <- pairs(Ec_emm, adjust="sidak")
Ec_emmeans_pwc <- as.data.frame(Ec_emmeans_pwc)
write_xlsx(Ec_emmeans_pwc, "posthoc_comparisons_EC_pairwise.xlsx")
head(Ec_emmeans_pwc)

Sum_dv = cld(Ec_emm,
             alpha   = 0.05,
             Letters = letters,         ###  Use lowercase letters for .group
             adjust  = "sidak")         ###  sidak-adjusted comparisons
Sum_dv
#str(Sum_dv)

Ecoli_posthoc_comparisons_dv <- as.data.frame(Sum_dv)
write_xlsx(Ecoli_posthoc_comparisons_dv,"posthoc_comparisons_Ecoli.xlsx")

str(Ecoli_posthoc_comparisons_dv)

#make as factor to graph:
Ecoli_posthoc_comparisons_dv$Temperature_centered <- as.factor(Ecoli_posthoc_comparisons_dv$Temperature_centered)
Ecoli_posthoc_comparisons_dv$Age_centered <- as.factor(Ecoli_posthoc_comparisons_dv$Age_centered)

head(Ecoli_posthoc_comparisons_dv)
levels(Ecoli_posthoc_comparisons_dv$Age_centered)

#make labels:
Ecoli_posthoc_comparisons_dv$Age_centered <- factor(Ecoli_posthoc_comparisons_dv$Age_centered,
                         labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

Ecoli_posthoc_comparisons_dv$Temperature_centered <- factor(Ecoli_posthoc_comparisons_dv$Temperature_centered,
                                 labels = c("27","30","32"))


Ecoli_emmeans1<- Ecoli_posthoc_comparisons_dv %>%
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(fill="#00BFC4"),
    stat="identity",
    width=0.8)+
  facet_grid(~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  scale_fill_manual(values="#00BFC4")+
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Ecoli_emmeans1

ggsave("Ecoli_emmeans1_png.png",plot=Ecoli_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_emmeans1_eps.eps",plot=Ecoli_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_emmeans1_pdf.pdf",plot=Ecoli_emmeans1,width = 6, height = 4, units = "in", dpi = 600)


Ecoli_posthoc_comparisons_dv$Age_centered <- factor(Ecoli_posthoc_comparisons_dv$Age_centered,
                                    labels = c("1","5","10","15"))

Ecoli_posthoc_comparisons_dv$Temperature_centered <- factor(Ecoli_posthoc_comparisons_dv$Temperature_centered,
                                 labels = c("`27˚C`","`30˚C`","`32˚C`"))

Ecoli_emmeans2 <- Ecoli_posthoc_comparisons_dv %>%
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(fill="#00BFC4"),
           stat="identity",
           width=0.8)+
  facet_grid(~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+ 
  scale_fill_manual(values="#00BFC4")+
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Ecoli_emmeans2


ggsave("Ecoli_emmeans2_png.png",plot=Ecoli_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_emmeans2_eps.eps",plot=Ecoli_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Ecoli_emmeans2_pdf.pdf",plot=Ecoli_emmeans2,width = 6, height = 4, units = "in", dpi = 600)

#mean by temp or age only

Ecoli_Ecemm <- as.data.frame(Ec_emm)
head(Ecoli_Ecemm)

Ecoli_Ecemm_TEMPonly <- Ecoli_Ecemm %>%
  group_by(Temperature_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(Ecoli_Ecemm_TEMPonly) #gives average ZOI for each temp across all ages

Ecoli_Ecemm_AGEonly <- Ecoli_Ecemm %>%
  group_by(Age_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(Ecoli_Ecemm_AGEonly) #gives average ZOI for each age across all temps

Ecoli_Ecemm_TEMPonly <- as.data.frame(Ecoli_Ecemm_TEMPonly)
Ecoli_Ecemm_AGEonly <- as.data.frame(Ecoli_Ecemm_AGEonly)

write_xlsx(Ecoli_Ecemm_TEMPonly, "Ecoli_Ecemm_TEMPonly.xlsx")
write_xlsx(Ecoli_Ecemm_AGEonly, "Ecoli_Ecemm_AGEonly.xlsx")

Ecoli_Ecemm_TEMPonly$Temperature_centered <- as.factor(Ecoli_Ecemm_TEMPonly$Temperature_centered)
Ecoli_Ecemm_AGEonly$Age_centered <- as.factor(Ecoli_Ecemm_AGEonly$Age_centered)

Ecoli_Ecemm_AGEonly$Age_centered <- factor(Ecoli_Ecemm_AGEonly$Age_centered,
                                                    labels = c("1","5","10","15"))

Ecoli_Ecemm_TEMPonly$Temperature_centered <- factor(Ecoli_Ecemm_TEMPonly$Temperature_centered,
                                                            labels = c("27","30","32"))


Ec_Temponly <- Ecoli_Ecemm_TEMPonly %>%
  ggplot(aes(x=Temperature_centered,y=mean))+
  geom_bar(aes(fill= c("#4D6FAE","#6F9F51", "#CC763B")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  scale_fill_manual(values=c("#4D6FAE","#6F9F51", "#CC763B"))+
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Ec_Temponly

ggsave("Ec_Temponly_png.png",plot=Ec_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Ec_Temponly_eps.eps",plot=Ec_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Ec_Temponly_pdf.pdf",plot=Ec_Temponly,width = 4, height = 4, units = "in", dpi = 600)


Ec_Ageonly <- Ecoli_Ecemm_AGEonly %>%
  ggplot(aes(x=Age_centered,y=mean))+
  geom_bar(aes(fill= c("#DCD1E9","#BAA4D3","#9776BE","#7549A8")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+   
  scale_fill_manual(values=c("#DCD1E9","#BAA4D3","#9776BE","#7549A8"))+
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Ec_Ageonly

ggsave("Ec_Ageonly_png.png",plot=Ec_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Ec_Ageonly_eps.eps",plot=Ec_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Ec_Ageonly_pdf.pdf",plot=Ec_Ageonly,width = 4, height = 4, units = "in", dpi = 600)

ref_grid(lmer_Ecoli_6)
EC_emm_temp <- emmeans(lmer_Ecoli_6,poly(Temperature_centered,2)|poly(Age_centered,2),
                       type="response",
                       at = list(poly(Temperature_centered,2) = c(-1.2,1.19,0.23)))
EC_emm_temp

EC_emm_temp <- emmeans(lmer_Ecoli_6,~Temperature_centered*Age_centered,
                  type="response",
                  at = list(Temperature_centered = c(-1.2,1.19,0.23)))

pairs(EC_emm_temp, adjust="sidak")

EC_emm_temp = cld(EC_emm_temp,
                  alpha   = 0.05,
                  Letters = letters,         ###  Use lowercase letters for .group
                  adjust  = "sidak")         ###  sidak-adjusted comparisons
EC_emm_temp

EC_posthoc_comparisons_temp <- as.data.frame(EC_emm_temp)
EC_posthoc_comparisons_temp
write_xlsx(EC_posthoc_comparisons_temp,"posthoc_comparisons_EC_temp.xlsx")

EC_emm_age <- emmeans(lmer_Ecoli_6,~Temperature_centered*Age_centered,
                      type="response",
                      at = list(Age_centered = c(-1.25,-0.47,0.5,1.48)))


pairs(EC_emm_age, adjust="sidak")

EC_dv_age = cld(EC_emm_age,
                alpha   = 0.05,
                Letters = letters,         ###  Use lowercase letters for .group
                adjust  = "sidak")         ###  sidak-adjusted comparisons
EC_dv_age

EC_posthoc_comparisons_dv_age <- as.data.frame(EC_dv_age)
EC_posthoc_comparisons_dv_age
write_xlsx(EC_posthoc_comparisons_dv_age,"posthoc_comparisons_EC_age.xlsx")


###########################

######################################################################
######################################################################
###
#M luteus:
#full model
lmer_Mlut_1 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                      poly(Age_centered,2)+
                      poly(Temperature_centered*Age_centered,2)+
                      (1|Sample_ID)+
                      (1|Plate_ID_Total)+
                    (1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=FALSE)

summary(lmer_Mlut_1)
plot(lmer_Mlut_1)
rand(lmer_Mlut_1)
Anova(lmer_Mlut_1)
drop1(lmer_Mlut_1)


lmer_Mlut_2 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                      poly(Age_centered,2)+
                      poly(Temperature_centered*Age_centered,2)+
                      #(1|Sample_ID)+
                      (1|Plate_ID_Total),
                      #(1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=FALSE)

summary(lmer_Mlut_2)
plot(lmer_Mlut_2)
rand(lmer_Mlut_2)
Anova(lmer_Mlut_2)
drop1(lmer_Mlut_2)

lmer_Mlut_3 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                      poly(Age_centered,2)+
                      Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      (1|Plate_ID_Total),
                    #(1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=FALSE)

summary(lmer_Mlut_3)
plot(lmer_Mlut_3)
rand(lmer_Mlut_3)
Anova(lmer_Mlut_3)
drop1(lmer_Mlut_3)

lmer_Mlut_4 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                      poly(Age_centered,2)+
                      #Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      (1|Plate_ID_Total),
                    #(1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=FALSE)
summary(lmer_Mlut_4)
plot(lmer_Mlut_4)
rand(lmer_Mlut_4)
Anova(lmer_Mlut_4)
drop1(lmer_Mlut_4)

anova(lmer_Mlut_3,lmer_Mlut_4)

BIC(lmer_Mlut_1,lmer_Mlut_2,lmer_Mlut_3,lmer_Mlut_4)
AIC(lmer_Mlut_1,lmer_Mlut_2,lmer_Mlut_3,lmer_Mlut_4)
anova(lmer_Mlut_1,lmer_Mlut_2,lmer_Mlut_3,lmer_Mlut_4)
#4 is best
#no interaction btwn temp and age for mlut

#temperature and age both significantly shape m.lut response but do not interact together

#go with lmer_Mlut_4

#rerun as REML= TRUE
lmer_Mlut_4 <- lmer(log(ZOI_area) ~ poly(Temperature_centered,2) +
                      poly(Age_centered,2)+
                      #Temperature_centered*Age_centered+
                      #(1|Sample_ID)+
                      (1|Plate_ID_Total),
                    #(1|Batch_Number),
                    data=Mlut_ZOI,
                    REML=TRUE)


#save outputs:
sink("Mlut_summary.txt")
print(summary(lmer_Mlut_4))
sink()  # returns output to the console

sink("Mlut_ANOVAoutput.txt")
print(anova(lmer_Mlut_4,type="2",ddf="Kenward-Roger"))
sink()

Mlut_anova <- parameters::model_parameters(anova(lmer_Mlut_4,type="2",ddf="Kenward-Roger"))
Mlut_anova<-as.data.frame(Mlut_anova)
Mlut_anova
write_xlsx(Mlut_anova,"Mlut_ANOVA.xlsx")

modelvalues_Mlut <- parameters::model_parameters(lmer_Mlut_4, exponentiate = TRUE)
head(modelvalues_Mlut)

write_xlsx(modelvalues_Mlut, "Mlut_model_expcoefandCIs.xlsx")

Mlut_effectsize_partial <- eta_squared(lmer_Mlut_4)
Mlut_effectsize_partial_dataframe <- as.data.frame(Mlut_effectsize_partial)
Mlut_effectsize_partial_dataframe
write_xlsx(Mlut_effectsize_partial_dataframe,"Mlut_partialeffectsizes.xlsx")

ref_grid(lmer_Mlut_4)

#emmip(lmer_Mlut_4,poly(Age_centered,2)~Temperature_centered, cov.reduce = range)
#emmip(lmer_Mlut_4,poly(Age_centered,1)~Temperature_centered, cov.reduce = range)

Mlut_emp1 <- emmip(lmer_Mlut_4,poly(Age_centered,2)~poly(Temperature_centered,2),style="factor",col="black",
      linearg = list(), dotarg = list(size = 3),
      cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(1,5,10,15))+
  scale_linetype(labels=c(1,5,10,15))+
  guides(shape = guide_legend(title = "Adult Age (days)"),
         linetype = guide_legend(title = "Adult Age (days)"))+
  scale_y_continuous(limits=c(0,35),breaks = c(0,10,20,30))+
  scale_x_discrete(breaks = c(-1.2,0.23,1.19),labels=c(27,30,32))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +  theme(legend.position = "right")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))
Mlut_emp1

Mlut_emp2 <- emmip(lmer_Mlut_4,poly(Temperature_centered,2)~poly(Age_centered,2),style="factor",col="black",
                   linearg = list(), dotarg = list(size = 3),
                   cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(1,5,10,15))+
  scale_linetype(labels=c(1,5,10,15))+
  guides(shape = guide_legend(title = "Adult Age (days)"),
         linetype = guide_legend(title = "Adult Age (days)"))+
  scale_y_continuous(limits=c(0,35),breaks = c(0,10,20,30))+
  scale_x_discrete(breaks = c(-1.2,0.23,1.19),labels=c(27,30,32))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +  theme(legend.position = "right")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))
Mlut_emp1

Mlut_emp2<- emmip(lmer_Mlut_4,poly(Temperature_centered,2)~poly(Age_centered,2),
                   style="factor",col="black",
                   linearg = list(), dotarg = list(size = 3),
                   cov.reduce=FALSE,type="response",CIs = FALSE) +
  theme_pubr()+ 
  scale_shape(labels=c(27,30,32))+
  scale_linetype(labels=c(27,30,32))+
  guides(shape = guide_legend(title = "Temperature (˚C)"),
         linetype = guide_legend(title = "Temperature (˚C)"))+
  scale_y_continuous(limits=c(0,35),breaks = c(0,10,20,30))+
  scale_x_discrete(breaks = c(-1.25,-0.47,0.5,1.48),labels=c(1,5,10,15))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Adult Age (days)") +  theme(legend.position = "right")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))

Mlut_emp2

ggsave("Mlut_interaction1.png",plot=Mlut_emp1,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Mlut_interaction1.eps",plot=Mlut_emp1,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Mlut_interaction1.pdf",plot=Mlut_emp1,width = 5.5, height = 4, units = "in", dpi = 600)

ggsave("Mlut_interaction2.png",plot=Mlut_emp2,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Mlut_interaction2.eps",plot=Mlut_emp2,width = 5.5, height = 4, units = "in", dpi = 600)
ggsave("Mlut_interaction2.pdf",plot=Mlut_emp2,width = 5.5, height = 4, units = "in", dpi = 600)


###

Mlut_emm <- emmeans(lmer_Mlut_4,specs=c("Temperature_centered", "Age_centered"),
                  type="response",
                  at = list(Temperature_centered = c(-1.2,1.19,0.23),
                            Age_centered = c(-1.25,-0.47,0.5,1.48)))
Mlut_emm

#find emmeans for post hoc comparisons:
library(multcompView)
library(emmeans)
library(multcomp)

#marginal
Mlut_pwc<-pairs(Mlut_emm,
      adjust="sidak")

Mlut_pwc<- as.data.frame(Mlut_pwc)
write_xlsx(Mlut_pwc,"Mluteus_emmeans_pariwisecomparisons.xlsx")


Mlut_dv = cld(Mlut_emm,
             alpha   = 0.05,
             Letters = letters,         ###  Use lowercase letters for .group
             adjust  = "sidak")         ###  sidak-adjusted comparisons
Mlut_dv
str(Mlut_dv)

Mlut_posthoc_comparisons_dv <- as.data.frame(Mlut_dv)
write_xlsx(Mlut_posthoc_comparisons_dv,"posthoc_comparisons_Mlut.xlsx")

#make as factor to graph:
Mlut_posthoc_comparisons_dv$Temperature_centered <- as.factor(Mlut_posthoc_comparisons_dv$Temperature_centered)
Mlut_posthoc_comparisons_dv$Age_centered <- as.factor(Mlut_posthoc_comparisons_dv$Age_centered)

head(Mlut_posthoc_comparisons_dv)
levels(Mlut_posthoc_comparisons_dv$Age_centered)

#make labels:
Mlut_posthoc_comparisons_dv$Age_centered <- factor(Mlut_posthoc_comparisons_dv$Age_centered,
                                                    labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

Mlut_posthoc_comparisons_dv$Temperature_centered <- factor(Mlut_posthoc_comparisons_dv$Temperature_centered,
                                                            labels = c("27","30","32"))


Mlut_emmeans1<- Mlut_posthoc_comparisons_dv %>%
  ggplot(aes(x=Temperature_centered,y=response))+
  geom_bar(aes(fill="#C77CFF"),
           stat="identity",
           width=0.8)+
  facet_grid(~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  scale_fill_manual(values="#C77CFF")+
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Mlut_emmeans1

ggsave("Mlut_emmeans1_png.png",plot=Mlut_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Mlut_emmeans1_eps.eps",plot=Mlut_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Mlut_emmeans1_pdf.pdf",plot=Mlut_emmeans1,width = 6, height = 4, units = "in", dpi = 600)


Mlut_posthoc_comparisons_dv$Age_centered <- factor(Mlut_posthoc_comparisons_dv$Age_centered,
                                                    labels = c("1","5","10","15"))

Mlut_posthoc_comparisons_dv$Temperature_centered <- factor(Mlut_posthoc_comparisons_dv$Temperature_centered,
                                                            labels = c("`27˚C`","`30˚C`","`32˚C`"))

Mlut_emmeans2 <- Mlut_posthoc_comparisons_dv %>%
  ggplot(aes(x=Age_centered,y=response))+
  geom_bar(aes(fill="#C77CFF"),
           stat="identity",
           width=0.8)+
  facet_grid(~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=response - SE,
                    ymax=response + SE),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+ 
  scale_fill_manual(values="#C77CFF")+
 # scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,40))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Mlut_emmeans2


ggsave("Mlut_emmeans2_png.png",plot=Mlut_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Mlut_emmeans2_eps.eps",plot=Mlut_emmeans2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Mlut_emmeans2_pdf.pdf",plot=Mlut_emmeans2,width = 6, height = 4, units = "in", dpi = 600)

########################################################################
#Mlut_emcont <- emmeans(lmer_Mlut_4, "Temperature_centered")
#ref_grid(Mlut_emcont)
#ply <- contrast(Mlut_emcont,"poly")
#ply
lmer_Mlut_4
library(multcomp)
summary(glht(lm_Naive_3, linfct = mcp(factor(Temperature_centered) = "Tukey")),test=adjusted("holm"))
library(emmeans)
#pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#emm = emmeans(lmer_Mlut_4, ~ Temperature_centered + Age_centered,)
pwpp(Mlut_emm)
str(Mlut_emm)
pwpp(Mlut_emm, 
     method = XXXX,
     adjust="sidak",
     type = "response", side = ">",aes=my.aes)

# custom aesthetics:
my.aes <- list(point = list(shape = "square"), 
               segment = list(linetype = "dashed", color = "red"),
               label = list(family = "serif", fontface = "italic"))
my.pal <- c("darkgreen", "blue", "magenta", "orange")
pwpp(Mlut_emm, aes = my.aes) + ggplot2::scale_color_manual(values = my.pal)


###################################################
#https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html
contrast(Mlut_emm, simple=list(c("Temperature_centered","Age_centered")),combine = TRUE, adjust = "sidak")
#contrast(noise.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")
Mlut_emm2 <- emmeans(lmer_Mlut_4,specs=c(factor("Temperature_centered"), factor("Age_centered")),
                    type="response")
#Mlut_emm2

mvcontrast(Naive_emm, "pairwise", type="response", adjust="sidak",mult.name = c("Age_centered"),show.ests = TRUE)
mvcontrast(Naive_emm, "pairwise", type="response", adjust="sidak",mult.name = c("Temperature_centered"),show.ests = TRUE)


mvcontrast(Ec_emm, "pairwise", type="response", adjust="sidak",mult.name = c("Age_centered"),show.ests = TRUE)
mvcontrast(Ec_emm, "pairwise",type="response", adjust="sidak", mult.name = c("Temperature_centered"),show.ests = TRUE)
mvcontrast(Mlut_emm, "pairwise", type="response", adjust="sidak",mult.name = c("Age_centered"),show.ests = TRUE)
mvcontrast(Mlut_emm, "pairwise",type="response", adjust="sidak", mult.name = c("Temperature_centered"),show.ests = TRUE)




str(lm_Naive_3_factor)
emtrends(lm_Naive_3_factor,pairwise~Temperature_centered,var="Age_centered")
emtrends(fiber.lm, pairwise ~ machine, var = "diameter")
emtrends(org.int, pairwise ~ variety, var = "price1", mult.name = "variety")

pigs.poly <- lm(conc ~ poly(percent, degree = 3), data = pigs)
emt <- emtrends(pigs.poly, ~ degree | percent, "percent", max.degree = 3,
                at = list(percent = c(9, 13.5, 18)))
# note: 'degree' is an extra factor created by 'emtrends'

summary(emt, infer = c(TRUE, TRUE))

######################################################
Mlut_emm_$Temperature_centered <- factor(Mlut_emm$Temperature_centered,
                                                    labels = c("27","30","32"))

#plot pairwise p value comparisons:
pwpp(Mlut_emm, by = "Temperature_centered", type = "response")+
  # ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(y="Adult Age (days)")+ 
  theme_pubr()+
  facet_grid(~Temperature_centered,labeller=c(27,30,32))
#  scale_y_continuous(limits=c(0,35),breaks = c(0,10,20,30))+
  scale_y_discrete(breaks = c(-1.25,-0.47,0.5,1.48),labels=c(1,5,10,15))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.2, "lines"))+
  theme(legend.position = "none")


#mean by temp or age only

Mlut_Ecemm <- as.data.frame(Mlut_emm)
head(Mlut_Ecemm)

Mlut_Ecemm_TEMPonly <- Mlut_Ecemm %>%
  group_by(Temperature_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(Mlut_Ecemm_TEMPonly) #gives average ZOI for each temp across all ages

Mlut_Ecemm_AGEonly <- Mlut_Ecemm %>%
  group_by(Age_centered) %>%
  reframe(
    mean = mean(response),
    n = n(),
    SE = sd(response)/sqrt(n()))
head(Mlut_Ecemm_AGEonly) #gives average ZOI for each age across all temps

Mlut_Ecemm_TEMPonly <- as.data.frame(Mlut_Ecemm_TEMPonly)
Mlut_Ecemm_AGEonly <- as.data.frame(Mlut_Ecemm_AGEonly)

write_xlsx(Mlut_Ecemm_TEMPonly, "Mlut_Ecemm_TEMPonly.xlsx")
write_xlsx(Mlut_Ecemm_AGEonly, "Mlut_Ecemm_AGEonly.xlsx")

Mlut_Ecemm_TEMPonly$Temperature_centered <- as.factor(Mlut_Ecemm_TEMPonly$Temperature_centered)
Mlut_Ecemm_AGEonly$Age_centered <- as.factor(Mlut_Ecemm_AGEonly$Age_centered)

Mlut_Ecemm_AGEonly$Age_centered <- factor(Mlut_Ecemm_AGEonly$Age_centered,
                                           labels = c("1","5","10","15"))

Mlut_Ecemm_TEMPonly$Temperature_centered <- factor(Mlut_Ecemm_TEMPonly$Temperature_centered,
                                                    labels = c("27","30","32"))


Mlut_Temponly <- Mlut_Ecemm_TEMPonly %>%
  ggplot(aes(x=Temperature_centered,y=mean))+
  geom_bar(aes(fill= c("#4D6FAE","#6F9F51", "#CC763B")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  #scale_fill_npg()+
  scale_fill_manual(values=c("#4D6FAE","#6F9F51", "#CC763B"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Mlut_Temponly

ggsave("Mlut_Temponly_png.png",plot=Mlut_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Mlut_Temponly_eps.eps",plot=Mlut_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Mlut_Temponly_pdf.pdf",plot=Mlut_Temponly,width = 4, height = 4, units = "in", dpi = 600)


Mlut_Ageonly <- Mlut_Ecemm_AGEonly %>%
  ggplot(aes(x=Age_centered,y=mean))+
  geom_bar(aes(fill= c("#DCD1E9","#BAA4D3","#9776BE","#7549A8")),
           stat = "identity", 
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  geom_errorbar(aes(ymin=(mean - SE),
                    ymax=(mean + SE)),
                width=0.4,
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+   
  scale_fill_manual(values=c("#DCD1E9","#BAA4D3","#9776BE","#7549A8"))+
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
Mlut_Ageonly

ggsave("Mlut_Ageonly_png.png",plot=Mlut_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Mlut_Ageonly_eps.eps",plot=Mlut_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("Mlut_Ageonly_pdf.pdf",plot=Mlut_Ageonly,width = 4, height = 4, units = "in", dpi = 600)

##################################

ref_grid(lmer_Mlut_4)

Mlut_emm_temp <- emmeans(lmer_Mlut_4,specs=c("Temperature_centered"),
                       type="response",
                       at = list(Temperature_centered = c(-1.2,1.19,0.23)))
Mlut_emm_temp

pairs(Mlut_emm_temp, adjust="sidak")

Mlut_emm_temp = cld(Mlut_emm_temp,
                  alpha   = 0.05,
                  Letters = letters,         ###  Use lowercase letters for .group
                  adjust  = "sidak")         ###  sidak-adjusted comparisons
Mlut_emm_temp

Mlut_posthoc_comparisons_temp <- as.data.frame(Mlut_emm_temp)
Mlut_posthoc_comparisons_temp
write_xlsx(Mlut_posthoc_comparisons_temp,"posthoc_comparisons_Mlut_temp.xlsx")

Mlut_emm_age <- emmeans(lmer_Mlut_4,specs=c("Age_centered"),
                      type="response",
                      at = list(Age_centered = c(-1.25,-0.47,0.5,1.48)))


pairs(Mlut_emm_age, adjust="sidak")

Mlut_dv_age = cld(Mlut_emm_age,
                alpha   = 0.05,
                Letters = letters,         ###  Use lowercase letters for .group
                adjust  = "sidak")         ###  sidak-adjusted comparisons
Mlut_dv_age

Mlut_posthoc_comparisons_dv_age <- as.data.frame(Mlut_dv_age)
Mlut_posthoc_comparisons_dv_age
write_xlsx(Mlut_posthoc_comparisons_dv_age,"posthoc_comparisons_Mlut_age.xlsx")






######################################################################
######################################################################
###





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
