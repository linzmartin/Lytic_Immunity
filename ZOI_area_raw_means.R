#antimicrobial activity

#raw data means (not estimated marginal means) without water control

#last updated: 05/07/24 LEM
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
ZOI_data <- read_xlsx("Antimicrobial-activity-rawdata_LEM042524.xlsx")

str(ZOI_data)
head(ZOI_data)

#format the data: 
#data is numeric / continuous
ZOI_data$Age <- as.numeric(ZOI_data$Age)
ZOI_data$Temperature <- as.numeric(ZOI_data$Temperature)
ZOI_data$Treatment <- as.factor(ZOI_data$Treatment)
ZOI_data$Diameter_avg <- as.numeric(ZOI_data$Diameter_avg)
ZOI_data$ZOI_area <- as.numeric(ZOI_data$ZOI_area)
str(ZOI_data)

ZOI_data$Injected_CFU_calculated_dose <- as.numeric(ZOI_data$Injected_CFU_calculated_dose)
ZOI_data$Infectious_Dose_for_Plates <- as.numeric(ZOI_data$Infectious_Dose_for_Plates)

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

###############################
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

ZOI_summary$Age_centered <- factor(ZOI_summary$Age_centered,
                                                   labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

ZOI_summary$Temperature_centered <- factor(ZOI_summary$Temperature_centered,
                                                           labels = c("27","30","32"))


ZOI_summary$Treatment <- factor(ZOI_summary$Treatment,    # Change factor labels
                                                labels = c("Naïve","Injury",
                                                           "italic(`E. coli`)",
                                                           "italic(`M. luteus`)"))
raw_means1<- ZOI_summary %>%
  ggplot(aes(x=Temperature_centered,y=mean_area))+
  geom_bar(aes(fill=Treatment),
           stat="identity",
           width=0.8)+
  facet_grid(Treatment~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
                width=0.4,stat="identity",
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  xlab("Temperature (˚C)") +
  theme_pubr()+ 
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  #scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")
raw_means1


ggsave("raw_means1_png.png",plot=raw_means1,width = 6, height = 5, units = "in", dpi = 600)
ggsave("raw_means1_eps.eps",plot=raw_means1,width = 6, height = 5, units = "in", dpi = 600)
ggsave("raw_means1_pdf.pdf",plot=raw_means1,width = 6, height = 5, units = "in", dpi = 600)


ZOI_summary$Age_centered <- factor(ZOI_summary$Age_centered,
                                                   labels = c("1","5","10","15"))

ZOI_summary$Temperature_centered <- factor(ZOI_summary$Temperature_centered,
                                                           labels = c("`27˚C`","`30˚C`","`32˚C`"))

ZOI_summary$Treatment <- factor(ZOI_summary$Treatment,    # Change factor labels
                                                labels = c("Naïve","Injury",
                                                           "italic(`E. coli`)",
                                                           "italic(`M. luteus`)"))

raw_means2 <- ZOI_summary %>%
  ggplot(aes(x=Age_centered,y=mean_area))+
  geom_bar(aes(fill=Treatment),
           stat="identity",
           width=0.8)+
  facet_grid(Treatment~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
                width=0.4,
                color="black")+
  #geom_text(aes(label = .group, y = response + SE), size=3.5,vjust = -0.5) +
  ylab(expression("Area of Zone of Inhibition"~(mm^2)~""))+ 
  labs(x="Adult Age (days)")+ 
  theme_pubr()+ 
  #scale_fill_npg()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  #scale_y_continuous(limits=c(0,50))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(10)),
        axis.text.y=element_text(size=10))+
  theme(legend.position = "none")

raw_means2

ggsave("raw_means2_png.png",plot=raw_means2,width = 6, height = 5, units = "in", dpi = 600)
ggsave("raw_means2_eps.eps",plot=raw_means2,width = 6, height = 5, units = "in", dpi = 600)
ggsave("raw_means2_pdf.pdf",plot=raw_means2,width = 6, height = 5, units = "in", dpi = 600)





Naive_ZOI <- subset(ZOI_data,Treatment == "Naïve")
LB_ZOI <- subset(ZOI_data,Treatment == "LB")
Ecoli_ZOI <- subset(ZOI_data,Treatment == "E_coli")
Mlut_ZOI <- subset(ZOI_data,Treatment == "M_luteus")



Naive_summary <- Naive_ZOI %>%
  group_by(Age_centered,Temperature_centered) %>%
  reframe(mean_diameter = mean(Diameter_avg),    mean_area = mean(ZOI_area),
    sd_diameter = sd(Diameter_avg),
    n_diameter = n(),
    SE_diameter = sd(Diameter_avg)/sqrt(n()),
    SE_area = sd(ZOI_area)/sqrt(n()))

Naive_summary <- as.data.frame(Naive_summary)
write_xlsx(Naive_summary,"Naive_raw_means.xlsx")


Naive_summary$Age_centered <- factor(Naive_summary$Age_centered,
                                   labels = c("`1 day`","`5 days`","`10 days`","`15 days`"))

Naive_summary$Temperature_centered <- factor(Naive_summary$Temperature_centered,
                                           labels = c("27","30","32"))

Naive_raw_means1<- Naive_summary%>%
  ggplot(aes(x=Temperature_centered,y=mean_area))+
  geom_bar(aes(fill="#F8766D"),
           stat="identity",
           width=0.8)+
  facet_grid(~Age_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
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
Naive_raw_means1

ggsave("Naive_emmeans1_png.png",plot=Naive_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_emmeans1_eps.eps",plot=Naive_emmeans1,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_emmeans1_pdf.pdf",plot=Naive_emmeans1,width = 6, height = 4, units = "in", dpi = 600)




Naive_summary$Age_centered <- factor(Naive_summary$Age_centered,
                                                    labels = c("1","5","10","15"))

Naive_summary$Temperature_centered <- factor(Naive_summary$Temperature_centered,
                                                            labels = c("`27˚C`","`30˚C`","`32˚C`"))

Naive_raw_means2 <- Naive_summary %>%
  ggplot(aes(x=Age_centered,y=mean_area))+
  geom_bar(aes(fill="#F8766D"),
           stat="identity",
           width=0.8)+
  facet_grid(~Temperature_centered,labeller = label_parsed)+
  geom_errorbar(aes(ymin=mean_area - SE_area,
                    ymax=mean_area + SE_area),
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
Naive_raw_means2


ggsave("Naive_raw_means2_png.png",plot=Naive_raw_means2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_raw_means2_eps.eps",plot=Naive_raw_means2,width = 6, height = 4, units = "in", dpi = 600)
ggsave("Naive_raw_means2_pdf.pdf",plot=Naive_raw_means2,width = 6, height = 4, units = "in", dpi = 600)

###########

Naive_ZOI <- as.data.frame(Naive_ZOI)
head(Naive_ZOI)

Naive_ZOI_TEMPonly <- Naive_ZOI %>%
  group_by(Temperature_centered) %>%
  reframe(
    mean = mean(ZOI_area),
    n = n(),
    SE = sd(ZOI_area)/sqrt(n()))
head(Naive_ZOI_TEMPonly) #gives average ZOI for each temp across all ages

Naive_ZOI_AGEonly <- Naive_ZOI %>%
  group_by(Age_centered) %>%
  reframe(
    mean = mean(ZOI_area),
    n = n(),
    SE = sd(ZOI_area)/sqrt(n()))
head(Naive_ZOI_AGEonly) #gives average ZOI for each age across all temps

Naive_ZOI_TEMPonly <- as.data.frame(Naive_ZOI_TEMPonly)
Naive_ZOI_AGEonly <- as.data.frame(Naive_ZOI_AGEonly)

write_xlsx(Naive_ZOI_TEMPonly, "Naive_ZOI_TEMPonly.xlsx")
write_xlsx(Naive_ZOI_AGEonly, "Naive_ZOI_AGEonly.xlsx")

Naive_ZOI_TEMPonly$Temperature_centered <- as.factor(Naive_ZOI_TEMPonly$Temperature_centered)
Naive_ZOI_AGEonly$Age_centered <- as.factor(Naive_ZOI_AGEonly$Age_centered)

Naive_ZOI_AGEonly$Age_centered <- factor(Naive_ZOI_AGEonly$Age_centered,
                                           labels = c("1","5","10","15"))

Naive_ZOI_TEMPonly$Temperature_centered <- factor(Naive_ZOI_TEMPonly$Temperature_centered,
                                                    labels = c("27","30","32"))


Naive_Temponly <- Naive_ZOI_TEMPonly %>%
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

ggsave("RAW-Naive_Temponly_png.png",plot=Naive_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("RAW-Naive_Temponly_eps.eps",plot=Naive_Temponly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("RAW-Naive_Temponly_pdf.pdf",plot=Naive_Temponly,width = 4, height = 4, units = "in", dpi = 600)

Naive_Ageonly <- Naive_ZOI_AGEonly %>%
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

ggsave("RAW-Naive_Ageonly_png.png",plot=Naive_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("RAW-Naive_Ageonly_eps.eps",plot=Naive_Ageonly,width = 4, height = 4, units = "in", dpi = 600)
ggsave("RAW-Naive_Ageonly_pdf.pdf",plot=Naive_Ageonly,width = 4, height = 4, units = "in", dpi = 600)



######################################################################
#LB

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



