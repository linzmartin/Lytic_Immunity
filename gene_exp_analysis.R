# clear existing workspace
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")

#set wd to your project folder
getwd() #check working directory

#sessionInfo()


#import the data:
gene_expression_data <- read_xlsx("lytic_expression.xlsx",sheet="all_genes")

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
                                 "Naïve","Live_Ecoli","Heat_killed_Ecoli"))

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
                                                           "italic(`E. coli`)",
                                                           "italic(`Heat killed E. coli`)"))
gene_expression_data_italic$Age <- factor(gene_expression_data_italic$Age,
                                          labels = c("1","10"))

gene_expression_data_italic$Temperature <- factor(gene_expression_data_italic$Temperature,
                                                  labels = c("`27˚C`","`30˚C`"))


gene_avg_summary_italic <- gene_avg_summary
gene_avg_summary_italic$Treatment <- factor(gene_avg_summary_italic$Treatment,    # Change factor labels
                                            labels = c("Naïve",
                                                       "italic(`E. coli`)",
                                                       "italic(`Heat killed E. coli`)"))
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
                                                       "italic(`E. coli`)",
                                                       "italic(`Heat killed E. coli`)"))
gene_avg_summary_italic$Age <- factor(gene_avg_summary_italic$Age,
                                      labels = c("1","10"))

gene_avg_summary_italic$Temperature <- factor(gene_avg_summary_italic$Temperature,
                                              labels = c("`27˚C`","`30˚C`"))

levels(gene_avg_summary$Gene)
gene_avg_summary_italic$Gene <- factor(gene_avg_summary_italic$Gene,
                                       labels=c("italic(`CECA`)",
                                                "italic(`DEF1`)",
                                                "italic(`LYSC1`)",
                                                "italic(`RPS17`)"))

png(filename = "gene_expression_all.png", width = 6, height = 6, units = "in", res = 300)
gene_avg_summary_italic %>%
  group_by(Gene)%>%
  ggplot(aes(x=Age,y=meanexp))+
  geom_bar(aes(fill = Treatment),
           stat = "identity", 
           position = position_dodge(1),
           width = 0.8) +
  scale_fill_manual(name="Treatment",
                    labels=c("Naïve", "Live Ec","Heat killed Ec"), 
                    values=c("skyblue","blue","darkblue"))+
  facet_grid(Gene~Temperature,labeller = label_parsed)+
  geom_errorbar(aes(ymin=meanexp - SE_exp,
                    ymax=meanexp + SE_exp,
                    group=Treatment),
                width=0.5,position=position_dodge(1),
                color="black")+
  ylab("Relative gene expression")+ 
  xlab("Adult age (days)") +
  theme_pubr()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(panel.spacing = unit(0.6, "lines"))+
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=(12)),
        axis.text.y=element_text(size=12))+
  theme(legend.position = "bottom")
dev.off()


################################
#Manova
library(car)

#import the data:
gene_expression_data_man <- read_xlsx("lytic_expression.xlsx",sheet="manova")

str(gene_expression_data_man)

gene_expression_data_man$Age <- as.factor(gene_expression_data_man$Age)
gene_expression_data_man$Temperature <- as.factor(gene_expression_data_man$Temperature)
gene_expression_data_man$Treatment <- as.factor(gene_expression_data_man$Treatment)

gene_expression_data_man <- as.data.frame(gene_expression_data_man)

#do a MANOVA:
model_1 <- lm(cbind(LysC1,CecA,Def1,rps17)~Temperature*Age*Treatment,data=gene_expression_data_man)

manova<-Manova(model_1, test.statistic = "Pillai")

summary(manova)

manova

model_1
