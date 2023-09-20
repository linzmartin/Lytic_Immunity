#water control comparison

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
ZOI_data <- read_xlsx("ZOI_assay_sample_log_082923_LEM_controlsduplicated.xlsx")

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
ZOI_waterandnaive <- subset(ZOI_data,Treatment == "Water_control" | Treatment=="Naïve")

ZOI_data_allgroups <- ZOI_data
ZOI_data <- subset(ZOI_data, Treatment != "Water_control")

str(ZOI_data)

#############
ZOI_wateronly 


ZOI_summary_water <- ZOI_wateronly %>%
  #group_by(Treatment,Age,Temperature) %>%
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

ZOI_sum_allgroups <- ZOI_data_allgroups %>%
  group_by(Treatment,Age,Temperature) %>%
  reframe(
    mean_area_all = mean(ZOI_area),
    n_all = n(),
    SE_area_all = sd(ZOI_area)/sqrt(n()))

str(ZOI_sum_allgroups)

ZOI_sum_allgroups$Age <- as.factor(ZOI_sum_allgroups$Age)
ZOI_sum_allgroups$Temperature <- as.factor(ZOI_sum_allgroups$Temperature)

#labels:
#italicize:
##
ZOI_italic <- ZOI_data_allgroups

ZOI_italic$Treatment <- factor(ZOI_italic$Treatment,    # Change factor labels
                               labels = c("Naïve","Injury",
                                          "italic(`E. coli`)",
                                          "italic(`M. luteus`)",
                                          "`Water Control`"))
#ZOI_italic$Age <- factor(ZOI_italic$Age,
#                                    labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
ZOI_italic$Age <- factor(ZOI_italic$Age,
                         labels = c("1","5","10","15"))

ZOI_italic$Temperature <- factor(ZOI_italic$Temperature,
                                 labels = c("`27˚C`","`30˚C`","`32˚C`"))
head(ZOI_italic)
levels(ZOI_italic$Treatment)


ZOI_sum_allgroups_italic <- ZOI_sum_allgroups
ZOI_sum_allgroups_italic$Treatment <- factor(ZOI_sum_allgroups_italic$Treatment,    # Change factor labels
                                               labels = c("Naïve","Injury",
                                                          "italic(`E. coli`)",
                                                          "italic(`M. luteus`)",
                                                          "`Water Control`"))
#ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
#                        labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
ZOI_sum_allgroups_italic$Age <- factor(ZOI_sum_allgroups_italic$Age,
                                         labels = c("1","5","10","15"))

ZOI_sum_allgroups_italic$Temperature <- factor(ZOI_sum_allgroups_italic$Temperature,
                                                 labels = c("`27˚C`","`30˚C`","`32˚C`"))

####
#plots with italic labels (use these!):

#mean area with jitter for points:
png(filename = "ZOI_area_jitter_with_watercontrol.png", width = 5, height = 6, units = "in", res = 300)
ZOI_sum_allgroups_italic %>%
  group_by(Treatment) %>%
  ggplot(aes(x=Age,y=mean_area_all))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Temperature,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area_all - SE_area_all,
                    ymax=mean_area_all + SE_area_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Age (days)") +
  theme_pubr()+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Age,y=ZOI_area),#color=ZOI_italic$Technical_Rep,
              position = position_dodge(0.5),size=1)+
  scale_color_manual(name="Treatment",
                     values=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#999999"))+
  scale_fill_manual(name="Treatment",
                     values=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#999999"))+
    theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

#change to put age at the top:
#need to adjust labels first:
#italicize:
##
ZOI_italic$Age <- factor(ZOI_italic$Age,
                         labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
                               
ZOI_italic$Temperature <- factor(ZOI_italic$Temperature,
                                 labels = c("27","30","32"))

ZOI_sum_allgroups_italic$Age <- factor(ZOI_sum_allgroups_italic$Age,
                                         labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

ZOI_sum_allgroups_italic$Temperature <- factor(ZOI_sum_allgroups_italic$Temperature,
                                                 labels = c("27","30","32"))

png(filename = "ZOI_area_jitter_ageattop_withwatercontrol.png", width = 6, height = 7, units = "in", res = 300)
ZOI_sum_allgroups_italic %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Temperature,y=mean_area_all))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area_all - SE_area_all,
                    ymax=mean_area_all + SE_area_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+
  theme(legend.position = "none")+
  scale_color_manual(name="Treatment",
                     values=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#999999"))+
  scale_fill_manual(name="Treatment",
                    values=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#999999"))+
    geom_jitter(data=ZOI_italic, aes(x=Temperature,y=ZOI_area),#color=ZOI_italic$Technical_Rep,
              position = position_dodge(0.5),size=1)+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()
