#lytic immunity - ZOI analysis

#last updated: 02/24/23 LEM
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
Imaging_Data <- read_xlsx("LEM_melanization_experimental_image_analysis_03_09_23.xlsx",
                          sheet = "Binary in ROI")

str(Imaging_Data)
Imaging_Data <- subset(Imaging_Data, select = -c(Item,Source,FieldID,RoiID,BinaryID,Threshold_Value,MaxIntensity,SumIntensity))
str(Imaging_Data)


#variables of interest: