#lytic immunity - ZOI analysis

#last updated: 04/03/23 LEM
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
ZOI_data <- read_xlsx("ZOI_assay_sample_log_040323.xlsx")

str(ZOI_data)
head(ZOI_data)

#variables of interest:


#format the data:
ZOI_data$Age <- as.factor(ZOI_data$Age)
ZOI_data$Temperature <- as.factor(ZOI_data$Temperature)
ZOI_data$Treatment <- as.factor(ZOI_data$Treatment)
str(ZOI_data)

levels(ZOI_data$Treatment)

#omit water group:
ZOI_data <- subset(ZOI_data, Treatment != "Water")

#create labels:
templabs <- c("27˚C","30˚C","32˚C")
names(templabs)<- c("27","30","32")

agelabs <- c("1 day","5 days", "10 days", "15 days")
names(agelabs)<-c("1","5","10","15")

treatlabs <-c("E. coli","LB","M. luteus","Naïve")
names(treatlabs)<-c("E_coli","LB","M_luteus","Naïve")

#relevel treatment to make Naive come first:

ZOI_data$Treatment <- relevel(ZOI_data$Treatment, ref="Naïve")
levels(ZOI_data$Treatment)
levels(ZOI_data$Age)

ZOI_data <- ZOI_data %>%
  mutate(Age = fct_relevel(Age, "1","5","10","15","NA"))


####
head(ZOI_data)
ZOI_data %>% ggplot(aes(x=Temperature,y=ZOI_area))+
  geom_point()+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Area of the Zone of Inhibition (mm2)"))+
  theme_bw()

ZOI_data %>% ggplot(aes(x=Temperature,y=ZOI_area))+
  geom_point()+ geom_boxplot()+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Area of the Zone of Inhibition (mm2)"))+
  theme_bw()
