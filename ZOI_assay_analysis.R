#lytic immunity - ZOI analysis

#last updated: 06/13/23 LEM
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
ZOI_data <- read_xlsx("ZOI_assay_sample_log_061323_LEM.xlsx")

str(ZOI_data)
head(ZOI_data)

#variables of interest:


#format the data:
ZOI_data$Age <- as.factor(ZOI_data$Age)
ZOI_data$Temperature <- as.factor(ZOI_data$Temperature)
ZOI_data$Treatment <- as.factor(ZOI_data$Treatment)
str(ZOI_data)

levels(ZOI_data$Treatment)

#relevel treatment to make Naive come first:
ZOI_data <- ZOI_data %>%
  mutate(Treatment = fct_relevel(Treatment, 
                                 "Naïve","LB","E_coli","M_luteus","Water_control","NA"))

#omit NA treatments (no sample)
ZOI_data <- subset(ZOI_data, Treatment != "NA")

#create labels:
templabs <- c("27˚C","30˚C","32˚C")
names(templabs)<- c("27","30","32")

agelabs <- c("1 day","5 days", "10 days", "15 days","Control")
names(agelabs)<-c("1","5","10","15","NA")

treatlabs <-c("E. coli","LB","M. luteus","Naïve","Control")
names(treatlabs)<-c("E_coli","LB","M_luteus","Naïve","Water_control")

#relevel treatment to make Naive come first:

#ZOI_data <- ZOI_data %>%
 # mutate(Treatment = fct_relevel(Treatment, 
  #                               "Naïve","LB","E_coli","M_luteus","Water_control","NA"))
#ZOI_data$Treatment <- relevel(ZOI_data$Treatment, ref="Naïve")
levels(ZOI_data$Treatment)
levels(ZOI_data$Age)

ZOI_data <- ZOI_data %>%
  mutate(Age = fct_relevel(Age, "1","5","10","15","NA"))

#subset data to exclude samples that did not have enough:
#remove NA values for temp and age
#omit NA treatments (no sample)
ZOI_data <- subset(ZOI_data, Enough!= 0)

#omit water group:
ZOI_wateronly <- subset(ZOI_data,Treatment == "Water_control")
ZOI_data <- subset(ZOI_data, Treatment != "Water_control")



#summary stats:

ZOI_summary <- ZOI_data %>%
  group_by(Treatment,Age,Temperature,Sample_ID) %>%
  reframe(#Temperature = Temperature, 
            #Age = Age,
            #Treatment=Treatment,
            #Sample_ID = Sample_ID,
           # Technical_Rep = Technical_Rep, 
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
  facet_grid(Age~Temperature,labeller = labeller(Temperature=templabs,Age=agelabs,Treatment=treatlabs))+
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


####
head(ZOI_data)
ZOI_data %>% ggplot(aes(x=Temperature,y=ZOI_area))+
  geom_point(color=ZOI_data$Technical_Rep)+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Area of the Zone of Inhibition (mm2)"))+
  theme_bw()

ZOI_data %>% ggplot(aes(x=Temperature,y=ZOI_area))+
  geom_point(color=ZOI_data$Technical_Rep)+ geom_boxplot()+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Area of the Zone of Inhibition (mm2)"))+
  theme_bw()




##############################################################

png(filename = "ZOI_diameter_summary.png", width = 12, height = 11, units = "in", res = 300)
ZOI_summary_bioreps %>%
  ggplot(aes(x=Treatment,y=mean_diameter_all,color=Treatment))+
  geom_bar(aes(fill=Treatment),color="black",
           stat="identity",
           position=position_dodge())+
  facet_grid(Age~Temperature,labeller = labeller(Age=agelabs,Treatment=treatlabs,Temperature=templabs))+
  geom_errorbar(aes(ymin=mean_diameter_all - SE_diameter_all,
                    ymax=mean_diameter_all + SE_diameter_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean diameter of ZOI"))+ 
  theme_bw()+
  theme(text = element_text(size = 26), axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")
dev.off()

png(filename = "ZOI_area_summary.png", width = 12, height = 11, units = "in", res = 300)
ZOI_summary_bioreps %>%
  ggplot(aes(x=Treatment,y=mean_area_all,color=Treatment))+
  geom_bar(aes(fill=Treatment),color="black",
           stat="identity",
           position=position_dodge())+
  facet_grid(Age~Temperature,labeller = labeller(Age=agelabs,Treatment=treatlabs,Temperature=templabs))+
  geom_errorbar(aes(ymin=mean_area_all - SE_area_all,
                    ymax=mean_area_all + SE_area_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean area of ZOI (mm2)"))+ 
  theme_bw()+
  theme(text = element_text(size = 26), axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")
dev.off()


#################

ZOI_data_with_dots_diameter_bytemp <- ZOI_data %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Temperature,y=Diameter_avg))+
  geom_point(color=ZOI_data$Technical_Rep)+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs,Treatment=treatlabs))+
  labs(x="Temperature (˚C)")+ ylab(expression("Zone of Inhibition Diameter (mm)"))+
  theme_bw()+
  guides(color = "none")+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=rel(1)))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="blue", width=0.2)+
  stat_summary(fun=mean, geom="point", color="blue",size=1)

ZOI_data_with_dots_diameter_bytemp

ZOI_summary_bioreps %>%
  ggplot(aes(x=Treatment,y=mean_area_all,color=Treatment))+
  geom_bar(aes(fill=Treatment),
           stat="identity",
           position=position_dodge())+
  facet_grid(Age~Temperature,labeller = labeller(Age=agelabs,Treatment=treatlabs,Temperature=templabs))+
  geom_errorbar(aes(ymin=mean_area_all - SE_area_all,
                    ymax=mean_area_all + SE_area_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean area of ZOI (mm2)"))+ 
  theme_bw()+
  #theme(text = element_text(size = 26))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_data, aes(x=Treatment,y=ZOI_area),position = position_dodge(0.5),color="black")


#mean area with jitter for points:
ZOI_summary_bioreps %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Age,y=mean_area_all))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Temperature,labeller = labeller(Age=agelabs,Treatment=treatlabs,Temperature=templabs))+
  geom_errorbar(aes(ymin=mean_area_all - SE_area_all,
                    ymax=mean_area_all + SE_area_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  #scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean area of ZOI"~(mm^2)~""))+ 
  theme_bw()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_data, aes(x=Age,y=ZOI_area),
              position = position_dodge(0.5),color=ZOI_data$Technical_Rep)

#mean diameter with jitter for points:
ZOI_summary_bioreps %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Age,y=mean_diameter_all))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Temperature,labeller = labeller(Age=agelabs,Treatment=treatlabs,Temperature=templabs))+
  geom_errorbar(aes(ymin=mean_diameter_all - SE_diameter_all,
                    ymax=mean_diameter_all + SE_diameter_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  #scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean diameter of ZOI"~(mm)~""))+ 
  theme_bw()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_data, aes(x=Age,y=Diameter_avg),
              position = position_dodge(0.5),color=ZOI_data$Technical_Rep)



###########################################
#italicize:
##
ZOI_italic <- ZOI_data

ZOI_italic$Treatment <- factor(ZOI_italic$Treatment,    # Change factor labels
                                           labels = c("Naïve","Injury",
                                                      "italic(`E. coli`)",
                                                      "italic(`M. luteus`)"))
ZOI_italic$Age <- factor(ZOI_italic$Age,
                                     labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
ZOI_italic$Age <- factor(ZOI_italic$Age,
                                     labels = c("1","5","10","15"))

ZOI_italic$Temperature <- factor(ZOI_italic$Temperature,
                                             labels = c("`27˚C`","`30˚C`","`32˚C`"))
head(ZOI_italic)
levels(ZOI_italic$Treatment)


ZOI_summary_bioreps_italic <- ZOI_summary_bioreps
ZOI_summary_bioreps_italic$Treatment <- factor(ZOI_summary_bioreps_italic$Treatment,    # Change factor labels
                               labels = c("Naïve","Injury",
                                          "italic(`E. coli`)",
                                          "italic(`M. luteus`)"))
ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
                         labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
                         labels = c("1","5","10","15"))

ZOI_summary_bioreps_italic$Temperature <- factor(ZOI_summary_bioreps_italic$Temperature,
                                 labels = c("`27˚C`","`30˚C`","`32˚C`"))

####
#plots with italic labels (use these!):

#mean diameter with jitter for points:
ZOI_summary_bioreps_italic %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Age,y=mean_diameter_all))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Temperature,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_diameter_all - SE_diameter_all,
                    ymax=mean_diameter_all + SE_diameter_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  #scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean diameter of ZOI"~(mm)~""))+ 
  theme_bw()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Age,y=Diameter_avg),
              position = position_dodge(0.5),color=ZOI_italic$Technical_Rep)


#mean area with jitter for points:
ZOI_summary_bioreps_italic %>%
  group_by(Treatment)%>%
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
  #scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Mean area of ZOI"~(mm^2)~""))+ 
  theme_bw()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Age,y=ZOI_area),
              position = position_dodge(0.5),color=ZOI_italic$Technical_Rep)


#########################################################################

#check data for normality:

#log-transform values

#model + anova

#posthoc