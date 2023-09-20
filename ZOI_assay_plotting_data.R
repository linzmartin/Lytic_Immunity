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
ZOI_data <- read_xlsx("ZOI_assay_sample_log_082923_LEM.xlsx")

str(ZOI_data)
head(ZOI_data)

#variables of interest:


#format the data:
ZOI_data$Age <- as.factor(ZOI_data$Age)
ZOI_data$Temperature <- as.factor(ZOI_data$Temperature)
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
ZOI_data <- subset(ZOI_data, Treatment != "NA")

#omit water group:
ZOI_wateronly <- subset(ZOI_data,Treatment == "Water_control")
ZOI_data <- subset(ZOI_data, Treatment != "Water_control")

str(ZOI_data)




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

###
ZOI_summary <- ZOI_data %>%
  group_by(Treatment,Age,Temperature) %>%
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

ZOI_data %>% ggplot(aes(x=Temperature,y=Diameter_avg))+
  geom_point(color=ZOI_data$Technical_Rep)+ geom_boxplot()+
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs, Treatment=treatlabs))+
  labs(title="Zone of Inhibition of Hemolymph",
       x="Temperature")+
  ylab(expression("Diameter of Zone of Inhibition (mm)"))+
  theme_bw()




##############################################################

png(filename = "ZOI_diameter_summary_081123.png", width = 12, height = 11, units = "in", res = 300)
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

png(filename = "ZOI_area_summary_081123.png", width = 12, height = 11, units = "in", res = 300)
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

png(filename = "ZOI_area_summary_081123_ageattop.png", width = 12, height = 11, units = "in", res = 300)
ZOI_summary_bioreps %>%
  ggplot(aes(x=Treatment,y=mean_area_all,color=Treatment))+
  geom_bar(aes(fill=Treatment),color="black",
           stat="identity",
           position=position_dodge())+
  facet_grid(Temperature~Age,labeller = labeller(Age=agelabs,Treatment=treatlabs,Temperature=templabs))+
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
#ZOI_italic$Age <- factor(ZOI_italic$Age,
 #                                    labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
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
#ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
 #                        labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
                         labels = c("1","5","10","15"))

ZOI_summary_bioreps_italic$Temperature <- factor(ZOI_summary_bioreps_italic$Temperature,
                                 labels = c("`27˚C`","`30˚C`","`32˚C`"))

####
#plots with italic labels (use these!):


#mean diameter with jitter for points:
png(filename = "ZOI_diameter_jitter_083023.png", width = 5, height = 4, units = "in", res = 300)
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
  ylab(expression("Diameter of Zone of Inhibition"~(mm)~""))+ xlab("Age (days)") +
  theme_pubr()+
  #theme_bw()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Age,y=Diameter_avg),#color=ZOI_italic$Technical_Rep,
              position = position_dodge(0.5),size=1)+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

#mean area with jitter for points:
png(filename = "ZOI_area_jitter_083023.png", width = 5, height = 4, units = "in", res = 300)
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
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Age (days)") +
  theme_pubr()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Age,y=ZOI_area),#color=ZOI_italic$Technical_Rep,
              position = position_dodge(0.5),size=1)+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

#change to put age at the top:
#need to adjust labels first:
#italicize:
##
ZOI_italic <- ZOI_data

ZOI_italic$Treatment <- factor(ZOI_italic$Treatment,    # Change factor labels
                               labels = c("Naïve","Injury",
                                          "italic(`E. coli`)",
                                          "italic(`M. luteus`)"))
ZOI_italic$Age <- factor(ZOI_italic$Age,
                                    labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
#ZOI_italic$Age <- factor(ZOI_italic$Age,
 #                        labels = c("1","5","10","15"))

#ZOI_italic$Temperature <- factor(ZOI_italic$Temperature,
 #                                labels = c("bquote(`27˚C`)","(`30˚C`)","(`32˚C`)"))
ZOI_italic$Temperature <- factor(ZOI_italic$Temperature,
                                 labels = c("27","30","32"))
head(ZOI_italic)
levels(ZOI_italic$Treatment)


ZOI_summary_bioreps_italic <- ZOI_summary_bioreps
ZOI_summary_bioreps_italic$Treatment <- factor(ZOI_summary_bioreps_italic$Treatment,    # Change factor labels
                                               labels = c("Naïve","Injury",
                                                          "italic(`E. coli`)",
                                                          "italic(`M. luteus`)"))
ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
                        labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))
#ZOI_summary_bioreps_italic$Age <- factor(ZOI_summary_bioreps_italic$Age,
 #                                        labels = c("1","5","10","15"))

ZOI_summary_bioreps_italic$Temperature <- factor(ZOI_summary_bioreps_italic$Temperature,
                                                 labels = c("27","30","32"))

####

png(filename = "ZOI_diameter_jitter_083023_ageattop.png", width = 5, height = 4, units = "in", res = 300)
ZOI_summary_bioreps_italic %>%
  group_by(Treatment)%>%
  ggplot(aes(x=Temperature,y=mean_diameter_all))+
  geom_bar(aes(color = Treatment, 
               fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_shape_identity(guide="legend")+
  facet_grid(Treatment~Age,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_diameter_all - SE_diameter_all,
                    ymax=mean_diameter_all + SE_diameter_all),
                width=0.8,position=position_dodge(0.9),
                color="black")+
  #scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Diameter of Zone of Inhibition"~(mm)~""))+ xlab("Temperature (˚C)") +
  theme_pubr()+
  #theme_bw()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Temperature,y=Diameter_avg),#color=ZOI_italic$Technical_Rep,
              position = position_dodge(0.5),size=1)+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()

#mean area with jitter for points:
png(filename = "ZOI_area_jitter_083023_ageattop.png", width = 5, height = 4, units = "in", res = 300)
ZOI_summary_bioreps_italic %>%
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
  #scale_x_discrete(limits = c("Naïve","LB","E_coli","M_luteus"),labels=c("Naïve","LB","E. coli","M. luteus"))+
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+
  #theme(text = element_text(size = 26))+
  #theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme(legend.position = "none")+
  geom_jitter(data=ZOI_italic, aes(x=Temperature,y=ZOI_area),#color=ZOI_italic$Technical_Rep,
              position = position_dodge(0.5),size=1)+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.5, "lines"))
dev.off()
#########################################################################

#check data for normality:

normtestresults <- ZOI_data %>%
  # mutate(mean_area = Delta_OD) %>%
  group_by(Temperature, Treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~ shapiro.test(.x$ZOI_area)))
head(normtestresults)

library(broom)
normtest.glance <- normtestresults %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
normtest.glance #values are not normal!

boxplot(ZOI_area ~ Temperature+Age+Treatment,
        col=c("white","lightgray"),
        data = ZOI_data)

ZOI_data %>%
  ggplot(aes(x=Temperature,y=ZOI_area,color=Age))+
  facet_grid(Treatment~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
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

ZOIcontinuous <- ZOI_data

ZOIcontinuous$Temperature <- as.numeric(ZOIcontinuous$Temperature)
ZOIcontinuous$Age <- as.numeric(ZOIcontinuous$Age)
str(ZOIcontinuous)

ggplot(ZOIcontinuous, aes(log(ZOI_area))) + geom_histogram()

lm_full_continuous <- lmer(log(ZOI_area) ~ Temperature+Age+Treatment + 
                  Temperature*Age+
                  Temperature*Treatment+
                  Age*Treatment+
                  #Temperature*Age*Treatment+
                 # (1|Sample_ID)+
                   (1|Plate_ID_Total)+
                   (1|Batch_Number),
                data=ZOI_data,
                REML=FALSE)

plot(lm_full_continuous)

summary(lm_full_continuous)

Anova(lm_full_continuous,type=2,test.statistic = "Chisq") 

#basically the same thing

lm_full_continuous2 <- lmer(log(ZOI_area) ~ poly(Temperature,2)+poly(Age,2)+Treatment +
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

lmer_Naive_1 <- lmer(log(ZOI_area) ~ Temperature+
                       Age + 
                       Temperature*Age+
                        (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Batch_Number),
                     data=Naive_ZOI,
                     REML=FALSE)

summary(lmer_Naive_1)

rand(lmer_Naive_1)
drop1(lmer_Naive_1)

lm_Naive_2 <- lmer(log(ZOI_area) ~ Temperature+
                     Age + 
                     Temperature*Age+
                     #(1|Sample_ID)+
                     #(1|Plate_ID_Total)+
                     (1|Batch_Number),
                   data=Naive_ZOI,
                   REML=FALSE)
summary(lm_Naive_2)
plot(lm_Naive_2)
Anova(lm_Naive_2,type=2)

drop1(lm_Naive_2)

lm_Naive_3 <-  lmer(log(ZOI_area) ~ Temperature+
                      Age + 
                      #Temperature*Age+
                      #(1|Sample_ID)+
                      #(1|Plate_ID_Total)+
                      (1|Batch_Number),
                    data=Naive_ZOI,
                    REML=FALSE)
anova(lm_Naive_2,lm_Naive_3)
summary(lm_Naive_3)

plot(lm_Naive_3)
Anova(lm_Naive_3,type=2)
#takeaway: temp affects naive, but not age or interaction

#
lmer_LB_1 <- lmer(log(ZOI_area) ~ Temperature +
                    Age + 
                    Temperature*Age+
                    (1|Sample_ID)+
                    (1|Plate_ID_Total)+
                    (1|Batch_Number),
                     data=LB_ZOI,
                     REML=FALSE)

summary(lmer_LB_1)
rand(lmer_LB_1) #can get rid of all random effects 

lm_LB_2 <- lm(log(ZOI_area) ~ Temperature+
                   Age + 
                   Temperature*Age,
                 data=LB_ZOI)

summary(lm_LB_2)
plot(lm_LB_2)

Anova(lm_LB_2,type=2)
drop1(lm_LB_2)

lm_LB_3 <- lm(log(ZOI_area) ~ Temperature+
                Age,
                #Temperature*Age,
              data=LB_ZOI)

summary(lm_LB_3)
plot(lm_LB_3)

anova(lm_LB_2,lm_LB_3)
AIC(lm_LB_2,lm_LB_3)

#Age affects injury but not temp or interaction


#infection:
lmer_Ecoli_1 <- lmer(log(ZOI_area) ~ Temperature +
                       Age + 
                       Temperature*Age+
                       (1|Sample_ID)+
                       (1|Plate_ID_Total)+
                       (1|Batch_Number),
                     data=Ecoli_ZOI,
                     REML=FALSE)

summary(lmer_Ecoli_1)

rand(lmer_Ecoli_1) #get rid of all random effects


lmer_Ecoli_2 <- lm(log(ZOI_area)~ Temperature+
                       Age + 
                       Temperature*Age,
                       #(1|Sample_ID)+
                       #(1|Plate_ID),
                       #(1|Batch_Number),
                     data=Ecoli_ZOI)

summary(lmer_Ecoli_2)
Anova(lmer_Ecoli_2,type=2)
plot(lmer_Ecoli_2)

#temp has impact on E. coli zoi but not age or int of temp and age


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
