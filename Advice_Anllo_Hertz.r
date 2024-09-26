

##SCRIPT FOR 
#Prereg analyses for 
# "Experience and advice consequences shape information sharing strategies"
# by subcomarc


############### ################################
#Preloading settings, directories and libraries#
############### ################################

#housekeeping
rm(list = ls())
# setwd() #set your working directory here

#make sure you have the following files in the directory
#RawData_Adv_ANLLO_HERTZ.csv #the raw data
#RW_RL.stan #the stan RL model
#RWpe.r #the R script for the PE and choice probability calculation

#make sure you load the following libraries
library(lme4)
library(emmeans)
library("ggplot2")
library("viridis")
library("dplyr")
library("ggeffects")
library("car")
library("ggthemes")
library("RColorBrewer")
library("MetBrewer")
library("tidyverse")
library("rstan")
library("brms")
library("optimx")
library("dfoptim")
library("sdamr")
library("DHARMa")
# library("rethinking")
library("httpgd") #not to call if running in RSTUDIO
library("multimode")
library("mclust")

#set hgd to true to produce better plots #NOT TO RUN IF IN RSTUDIO
hgd()
hgd_browse()

### LOAD RAW DATA
#####

data <- read.csv(file='RawData_Adv_ANLLO_HERTZ.csv') # load the dataframe not to do this again

##### 
####BEGIN CLEANUP
##### 

#practice exclusions as described on the preprint
dataClean <- data
# Remove individual trials with choice/advice response time of above 15 seconds.
dataClean <- dataClean %>% filter(rt_choice<15000, rt_yesnochoice<15000)

##### 
####RL: invoke STAN
##### 

dataClean$blocktype<-as.numeric(as.character(dataClean$blocktype))
dataClean$id<-as.numeric(as.character(dataClean$id))

dat <- as.list(dataClean[,c("id", "choice","reward", "blocktype","trialnum")])
K=1
for (sub_id in unique(dat$id)){
  dat$id_2[dat$id==sub_id]=K
  K=K+1
}
dat$id=dat$id_2
dat$reward=dat$reward/10
dat$choice=dat$choice+1
dat$N <- nrow(dataClean)
dat$N_id <- length(unique(dataClean$id))

m <- stan(file = "RW_RL.stan", data = dat, iter = 8000, cores = 7, chains = 4, refresh = 10, control = list(adapt_delta = 0.9, max_treedepth = 14))

# saveRDS(m, "RL_RW_DATA.rds") #if you want to save the model results

LR_summary <- summary(m, pars = c("phi"), probs = c(0.1, 0.9))$summary #alpha, learning rate
Temp_summary <- summary(m, pars = c("lambda"), probs = c(0.1, 0.9))$summary #temperature
library(rethinking)

WAIC(m)

RF_Param<-data.frame(LR=LR_summary[,1],Temp=Temp_summary[,1],id=seq(1,K-1))
#write.csv(RF_Param,'StanRLparam.csv') #save the parameters

#Done with STAN

#Now let's integrate the parameters, 
#and use them to calculate PE and Choice_Prob

source("RWpe.r") #load the PE and choice probability calculation script
bcpdat<-dataClean
dataClean <- merge(dataClean, RF_Param, by="id")
dataClean$trialnum <- dataClean$trialnum+1

all_result <- data.frame("id"="", "trialnum"= "", "blocktype"="", "Q_choice"= "", "Q_not_choice"= "","p_choice"= "","PE"="")

for (sub_id in unique(dataClean$id)){
  
  #filter only the choices of the current subject, use %>% and filter function
  subject_choices <- dataClean %>% filter(id == sub_id)
  
  #fit models per each subject's choices, use optim, with initial values beta=1, alpha(s)=0.5
  #lower and upper set the boundaries for the free params
  
  subject_result<-RWmodelPE(subject_choices, c(subject_choices$Temp[1],subject_choices$LR[1]))
  colnames(subject_result) <- c("Q_choice", "Q_not_choice","p_choice","PE")
  subject_result <- data.frame("id"=subject_choices$id, "trialnum"=subject_choices$trialnum, "blocktype"= subject_choices$blocktype, subject_result)
  all_result<-rbind(all_result, subject_result)
}

# DF containing for each choice:"Q_choice", "Q_not_choice","p_choice","PE"
all_result$id <- as.numeric(all_result$id)
dataClean<-merge(dataClean,all_result, by=c("id","trialnum","blocktype"))
dataClean$p_choice<- as.numeric(dataClean$p_choice)
dataClean$PE<- as.numeric(dataClean$PE)
dataClean$exptype <- factor(dataClean$exptype)
dataClean$blocktype <- factor(dataClean$blocktype, labels=c("Hard","Medium","Easy"))

#write.csv(dataClean,'CleanDataWithRL.csv') #save this entire dataframe
# dataClean <- read.csv(file='CleanDataWithRL.csv') #load this entire dataframe

#Final manipulation: shift prediction error 
#so that we may see the effect of previous trial 
#PE on current trial advice-giving

PlayDat <- dataClean
PlayDat$shiftedPE <- 100
PlayDat <- sort_df(PlayDat, vars=c("blocktype","id","trialnum"))
ShiftedPEDF <- c()
for(l in 1:length(unique(PlayDat$id))){
 ThisPart <- PlayDat[PlayDat$id %in% unique(PlayDat$id)[l],]
  for(j in 1:length(unique(ThisPart$blocktype))){
       for(i in 1:length(unique(ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j],]$trialnum))){
         if(ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j],]$trialnum[i] == min(ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j],]$trialnum)){
           ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j] & ThisPart$trialnum==unique(ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j],]$trialnum)[i],]$shiftedPE = 0
         }else{
           ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j] & ThisPart$trialnum==unique(ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j],]$trialnum)[i],]$shiftedPE = ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j] & ThisPart$trialnum==unique(ThisPart[ThisPart$blocktype==unique(ThisPart$blocktype)[j],]$trialnum)[i-1],]$PE
         }
           }
         }
 ShiftedPEDF <- rbind(ShiftedPEDF, ThisPart)
}
dataClean <- ShiftedPEDF

#write.csv(dataClean,'CleanDataWithRL_SHIFTED.csv') #save this entire dataframe
# dataClean <- read.csv(file='CleanDataWithRL_SHIFTED.csv') #load this entire dataframe

#Robustness check: simulate data from the model and compare it to actual performance
#This is done by comparing the actual data to the model's predictions
for (sub_id in unique(dataClean$id)){
  #filter only the choices of the current subject, use %>% and filter function
  subject_choices <- dataClean %>% filter(id == sub_id)
  #fit models per each subject's choices, use optim, with initial values beta=1, alpha(s)=0.5
  #lower and upper set the boundaries for the free params
  subject_result<-RWmodelPE(subject_choices, c(subject_choices$Temp[1],subject_choices$LR[1]))
  colnames(subject_result) <- c("Q_choice", "Q_not_choice","p_choice","PE")
  subject_result <- data.frame("id"=subject_choices$id, "trialnum"=subject_choices$trialnum, "blocktype"= subject_choices$blocktype, subject_result)
  all_result<-rbind(all_result, subject_result)
}


##### 
#IDENTIFY OVERADVICE AND SET SEPARATE DFs FOR STUDY 1 AND 2
#####
#(and correct some typos)

dataExp1 <- dataClean %>% filter(exptype %in% "No cost") #load this entire dataframe
dataExp1$exptype <- factor(dataExp1$exptype)
dataExp1$blocktype <- factor(dataExp1$blocktype)

#
OverAdv <- dataExp1
DF <- dataExp1 %>% group_by(id) %>% summarise(id=id, MeanAdv = mean(yesnochoicenum)) %>% distinct()
DF <- DF %>% mutate(OA = ifelse(MeanAdv==1, 1,0))
dataExp1 <- merge(dataExp1, DF[,c(1,3)], by="id") 
dataExp1$OA <- factor(dataExp1$OA)

dataExp2 <- dataClean %>% filter(exptype %in% c("No cost (replication)", "Money","Reputation","Responsibilty")) #load this entire dataframe
dataExp2[dataExp2$exptype %in% "No cost (replication)",]$exptype <- "No Cost (replication)"
dataExp2[dataExp2$exptype %in% "Responsibilty",]$exptype <- "Responsibility"
dataExp2$exptype <- factor(dataExp2$exptype)
dataExp2$blocktype <- factor(dataExp2$blocktype)

#
OverAdv <- dataExp2
DF <- dataExp2 %>% group_by(id) %>% reframe(id=id, MeanAdv = mean(yesnochoicenum)) %>% distinct()
DF <- DF %>% mutate(OA = ifelse(MeanAdv==1, 1,0))
dataExp2 <- merge(dataExp2, DF[,c(1,3)], by="id") 
dataExp2$OA <- factor(dataExp2$OA)

#dataExp1 corresponds to Study 1
#dataExp2 corresponds to Study 2

##### 
#FIGURES AND ANALYSES FOR MAIN TEXT
#####

#STUDY 1
#####

#FIGURE 2: PANEL A

#The corresponding stats
ModPreg1 <- glmer(yesnochoicenum ~ blocktype*trialnum+scaledPTM+scaledSPIN+(1+blocktype+trialnum | id),family =  binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",  optCtrl = list(maxfun = 1000000)),  data=dataExp1[dataExp1$OA %in% c("0"),])
Anova(ModPreg1)

#prepare data for plot
ModA1forPlot <- dataExp1[dataExp1$OA %in% c("0"),] %>% group_by(blocktype,trialnum) %>% summarise(Accuracy=mean(choice), Advising = mean(yesnochoicenum),
                                                                        AccuracySD=sd(replicate(1000, mean(sample(choice, replace=T)))),
                                                                        AdvisingSD=sd(replicate(1000, mean(sample(yesnochoicenum, replace=T)))),
                                                                        PTM_all=first(PTM_all), SPIN_fear=first(SPIN_fear)) %>% distinct()

ModA1forPlot$blocktype <- fct_relevel(ModA1forPlot$blocktype, "Easy","Medium","Hard") 

#plot

ggplot() +
           geom_errorbar(data=ModA1forPlot, aes(x=trialnum, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD, color=blocktype), width=0.3) +
           geom_line(data=ModA1forPlot, aes(x=trialnum, y=Accuracy, color=blocktype), linewidth=3, size=2) +
           geom_point(data=ModA1forPlot, aes(x=trialnum, y=Accuracy, fill=blocktype), color="black", size=3, stroke=3, shape=21) +
           geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
           geom_errorbar(data=ModA1forPlot, aes(x=trialnum, y=Advising, ymax=Advising + AdvisingSD, ymin=Advising - AdvisingSD, color=blocktype), width=0.3, alpha=0.2) +
           geom_line(data=ModA1forPlot, aes(x=trialnum, y=Advising, color=blocktype), linewidth=3, size=2, alpha=0.5) +
           geom_point(data=ModA1forPlot, aes(x=trialnum, y=Advising, fill=blocktype), color="black", size=3, stroke=3, shape=21, alpha=0.5) +
          scale_fill_manual(aesthetics = c("color", "fill"), values = canva_palettes$`Primary colors with a vibrant twist`[1:3]) +
  facet_wrap(vars(blocktype)) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
labs(face="plain", family="arial",
     title = "Figure 2 panel A",
     x = "Trial", 
     y = "P(correct)")

#FIGURE 2: PANEL B

#The corresponding stats

ModPreg3 <- glmer(yesnochoicenum ~ p_choice*(scaledSPIN+scaledPTM)+shiftedPE+(1 | id),family =  binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",  optCtrl = list(maxfun = 1000000)),  data=dataExp1[dataExp1$OA %in% c("0"),])
Anova(ModPreg3)

#prepare data for plot
#this includes plotting binned behavioral data (just for illustrative purposes)
#and the model predictions/marginals (to understand the stat effects)

tertilesSPIN <- quantile(dataExp1[dataExp1$OA %in% "0",]$scaledSPIN, probs = c(1/3, 2/3, 3/3))
data3int <- dataExp1[dataExp1$OA %in% "0",]
data3int$pchoicebin <- 1 #again: illustration purposes
data3int[data3int$p_choice >= 0.8,]$pchoicebin <- 1
data3int[data3int$p_choice <= 0.8 & data3int$p_choice >= 0.6,]$pchoicebin <- 0.75
data3int[data3int$p_choice < 0.6 & data3int$p_choice >= 0.4,]$pchoicebin <- 0.50
data3int[data3int$p_choice <0.4,]$pchoicebin <- 0.25
data3int$SPINclass <- 1
data3int[data3int$scaledSPIN <= tertilesSPIN[[1]],]$SPINclass <- "Low"
data3int[data3int$scaledSPIN > tertilesSPIN[[1]] & data3int$scaledSPIN <= tertilesSPIN[[2]],]$SPINclass <- "Medium"
data3int[data3int$scaledSPIN > tertilesSPIN[[2]],]$SPINclass <- "High"
Plotchoicebin <- data3int %>% group_by(pchoicebin, SPINclass) %>% summarise(MeanAdvbin = mean(yesnochoicenum), 
                                                          SEPchoicebin =  sd(replicate(1000, mean(sample(yesnochoicenum, replace=T))))) %>% distinct()

DFforModel1 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy1", scaledSPIN=-0.5118994, scaledPTM=-0.02, SPINclass="Low")
DFforModel2 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy2", scaledSPIN=0.4782218, scaledPTM=-0.02, SPINclass="Medium")
DFforModel3 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy3", scaledSPIN=2.034127, scaledPTM=-0.02, SPINclass="High")

DFforModel<- rbind(DFforModel1,DFforModel2,DFforModel3)


#get marginals and ribbons (i.e. confidence interval for glmer fit)
#as prescribed by Ben Bolker in GLMM FAQ
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

ModA3bis <- glm(yesnochoicenum ~ p_choice*scaledSPIN,family =  binomial(link = "logit"),  data=dataExp1[dataExp1$OA %in% "0",]) 
DFforModel$modeledY <- predict(ModA3bis, type = 'response', newdata = DFforModel,  allow.new.levels = TRUE )
DFforModel$yesnochoicenum <- DFforModel$modeledY
des = model.matrix(formula(ModA3bis), DFforModel) #get the model matrix                                             
predvar = diag( des %*% vcov(ModA3bis) %*% t(des) ) # use matrix multiplication on the model matrix and variance-covariance matrix. Get diagonal of the resulting matrix
DFforModel$ConflowerRib = with(DFforModel, modeledY - 2*sqrt(predvar) ) #2 is fine, since SEM comes from 1.96
DFforModel$ConfupperRib = with(DFforModel, modeledY + 2*sqrt(predvar) )

#order factors for nicer plot
Plotchoicebin$SPINclass <- fct_relevel(Plotchoicebin$SPINclass, c("Low", "Medium","High")) #SPIN
DFforModel$SPINclass <- fct_relevel(DFforModel$SPINclass, c("Low", "Medium","High"))#SPIN

#and plot
ggplot() +
  geom_line(data=DFforModel, aes(x=p_choice, y=modeledY, color=SPINclass), linewidth=3, size=2, alpha=1) + #PTM plot
  geom_ribbon(data=DFforModel, aes(x=p_choice, ymin=ifelse(ConflowerRib<0, 0, ConflowerRib), ymax=ifelse(ConfupperRib>1,1,ConfupperRib), fill=SPINclass), linewidth=3, size=2, alpha=0.1) +
  geom_errorbar(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, ymax=MeanAdvbin + SEPchoicebin, ymin=MeanAdvbin - SEPchoicebin, color=SPINclass), linewidth=1, width=0.02, alpha=1) +
  geom_point(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, fill=SPINclass), color="black", size=3, stroke=2, shape=21, alpha=1) +
  geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
  scale_fill_manual(aesthetics = c("color","fill"), values = c("#2dbd20","#197c10","#0d3309")) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
  labs(face="plain", family="arial",
       title = "Figure 2 Panel B",
       x = "P(choice)", 
       y = "P(advice)")

#FIGURE 2: PANEL C

DtAvsI <- dataExp1[dataExp1$OA %in% c("0"),] %>% group_by(trialnum,yesnochoice) %>% summarise(Accuracy=mean(choice),
                                                                        AccuracySD=sd(replicate(1000, mean(sample(choice, replace=T)))),
                                                                        PTM_all=first(PTM_all), SPIN_fear=first(SPIN_fear)) %>% distinct()
ggplot() +
            geom_errorbar(data=DtAvsI[DtAvsI$yesnochoice %in% "0",], aes(x=trialnum, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD), color="black", width=0.3) +
            geom_line(data=DtAvsI[DtAvsI$yesnochoice %in% "0",], aes(x=trialnum, y=Accuracy), color="black", linewidth=3, size=2) +
            geom_point(data=DtAvsI[DtAvsI$yesnochoice %in% "0",], aes(x=trialnum, y=Accuracy), fill="black", color="black", size=3, stroke=3, shape=21) +
            geom_errorbar(data=DtAvsI[DtAvsI$yesnochoice %in% "1",], aes(x=trialnum, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD), color="black", alpha=0.3, width=0.3) +
            geom_line(data=DtAvsI[DtAvsI$yesnochoice %in% "1",], aes(x=trialnum, y=Accuracy), color="black", linewidth=3, alpha=0.3, size=2) +
            geom_point(data=DtAvsI[DtAvsI$yesnochoice %in% "1",], aes(x=trialnum, y=Accuracy), fill="black", color="black", size=3, alpha=0.3, stroke=3, shape=21) +
           geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
          scale_fill_manual(aesthetics = c("color", "fill"), values = canva_palettes$`Primary colors with a vibrant twist`[1:3]) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
labs(face="plain", family="arial",
     title = "Figure 2 panel C",
     x = "Trial", 
     y = "P(correct)")


# FIGURE 3

#The corresponding stats
ModPreg4 <- glmer(choice ~ exptype*(scaledPTM+scaledSPIN+trialnum)+(1 | id),family =  binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",  optCtrl = list(maxfun = 1000000)),  data=dataExp2[dataExp2$OA %in% c("0"),])
Anova(ModPreg4)

ModPreg7 <- glmer(yesnochoicenum ~ p_choice*exptype*(scaledPTM+scaledSPIN)+(1 | id),family =  binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",  optCtrl = list(maxfun = 1000000)),  data=dataExp2[dataExp2$OA %in% c("0"),])
Anova(ModPreg7) #these are also good for Figures 5

#FIGURE 3: Panel A

#prepare data for plot
ModA1forPlot <- dataExp2[dataExp2$OA %in% c("0"),] %>% group_by(trialnum, exptype) %>% summarise(Accuracy=mean(choice), Advising = mean(yesnochoicenum),
                                                                        AccuracySD=sd(replicate(1000, mean(sample(choice, replace=T)))),
                                                                        AdvisingSD=sd(replicate(1000, mean(sample(yesnochoicenum, replace=T)))),
                                                                        PTM_all=first(PTM_all), SPIN_fear=first(SPIN_fear)) %>% distinct()
                                                                        
ggplot() +
           geom_errorbar(data=ModA1forPlot, aes(x=trialnum, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD, color=exptype), width=0.3) +
           geom_line(data=ModA1forPlot, aes(x=trialnum, y=Accuracy, color=exptype), linewidth=3, size=2) +
           geom_point(data=ModA1forPlot, aes(x=trialnum, y=Accuracy, fill=exptype), color="black", size=3, stroke=3, shape=21) +
           geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
           geom_errorbar(data=ModA1forPlot, aes(x=trialnum, y=Advising, ymax=Advising + AdvisingSD, ymin=Advising - AdvisingSD, color=exptype), width=0.3, alpha=0.2) +
           geom_line(data=ModA1forPlot, aes(x=trialnum, y=Advising, color=exptype), linewidth=3, size=2, alpha=0.5) +
           geom_point(data=ModA1forPlot, aes(x=trialnum, y=Advising, fill=exptype), color="black", size=3, stroke=3, shape=21, alpha=0.5) +
          scale_fill_manual(aesthetics = c("color", "fill"), values = c("#14616B", "#ff2768", "#DD8016", "#b94b95")) +
  facet_wrap(vars(exptype), nrow = 4, ncol = 3) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        # legend.position = c(0.9,0.9),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
labs(face="plain", family="arial",
     title = "Figure 3 Panel A",
     x = "Trial", 
     y = "Ratio")


#FIGURE 3: Panel B
ModA1forPlotIND$interaCAT <- interaction(ModA1forPlotIND$exptype, ModA1forPlotIND$OA)
ggplot() +
           geom_bar(data=ModA1forPlotIND,  aes(x=exptype, fill=interaCAT), position="fill", alpha=0.5, width = 0.07) +
            scale_fill_manual(aesthetics = c("fill"), values = c("#6eb6c0", "#f080a2", "#d6ac7b", "#c990b6",
                                                                "#14616B", "#ff2768", "#DD8016", "#b94b95")) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
labs(face="plain", family="arial",
    title = "Figure 3 Panel B",
     x = "", 
     y = "P(advice)")

#FIGURE 3 PANEL C

#Prepare dataset

PlotEDF <- dataExp2[dataExp2$OA %in% c("0"),] #for illustration purposes(!!!), a segmented bin of behavioral data
PlotEDF$pchoicebin <- 1 
PlotEDF[PlotEDF$p_choice >= 0.8,]$pchoicebin <- 1
PlotEDF[PlotEDF$p_choice <= 0.8 & PlotEDF$p_choice >= 0.6,]$pchoicebin <- 0.75
PlotEDF[PlotEDF$p_choice < 0.6 & PlotEDF$p_choice >= 0.4,]$pchoicebin <- 0.50
PlotEDF[PlotEDF$p_choice <0.4,]$pchoicebin <- 0.25
PlotEDF$pchoicebin <- factor(PlotEDF$pchoicebin) 

Plotchoicebin <- PlotEDF %>% group_by(exptype, pchoicebin) %>% summarise(MeanAdvbin = mean(yesnochoicenum), SEPchoicebin =  sd(replicate(1000, mean(sample(yesnochoicenum, replace=T))))) %>% distinct() 
Plotchoicebin$pchoicebin <- as.numeric(as.character(Plotchoicebin$pchoicebin))

#get fit and calculate ribbons for confidence interval for the fit
#same approach as above, Bolker GLMM FAQ
ModA3bis <- glm(yesnochoicenum ~ exptype*p_choice,family =  binomial(link = "logit"),  data=PlotEDF) 
DFforModel <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy")
DFforModel$exptype <- fct_relevel(DFforModel$exptype, "Money","No Cost (replication)","Reputation","Responsibility") 
DFforModel$modeledY <- predict(ModA3bis, type = 'response', newdata = DFforModel,  allow.new.levels = TRUE )
DFforModel$yesnochoicenum <- DFforModel$modeledY
des = model.matrix(formula(ModA3bis), DFforModel) #get ribbons                                             
predvar = diag( des %*% vcov(ModA3bis) %*% t(des) )
DFforModel$ConflowerRib = with(DFforModel, modeledY - 2*sqrt(predvar) )
DFforModel$ConfupperRib = with(DFforModel, modeledY + 2*sqrt(predvar) )

ggplot() +  
  geom_line(data=DFforModel, aes(x=p_choice, y=modeledY, color=exptype), linewidth=3, size=2, alpha=1) +
  geom_ribbon(data=DFforModel, aes(x=p_choice, ymin=ConflowerRib, ymax=ConfupperRib, fill=exptype), linewidth=3, size=2, alpha=0.1) +
  geom_errorbar(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, ymax=MeanAdvbin + SEPchoicebin, ymin=MeanAdvbin - SEPchoicebin, color=exptype), linewidth=1, width=0.02, alpha=1) +
  geom_point(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, fill=exptype), color="black", size=3, stroke=2, shape=21, alpha=1) +
  scale_fill_manual(aesthetics = c("color","fill"), values = c("#14616B", "#ff2768", "#DD8016", "#b94b95")) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        # legend.title = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
  labs(face="plain", family="arial",
       title = "Figure 3 Panel C",
       x = "P(choice)", 
       y = "P(advice)") 

#FIGURE 4

DtAvs2 <- dataExp2[dataExp2$OA %in% c("0"),] %>% group_by(trialnum,yesnochoice, exptype) %>% summarise(Accuracy=mean(choice),
                                                                        AccuracySD=sd(replicate(1000, mean(sample(choice, replace=T)))),
                                                                        PTM_all=first(PTM_all), SPIN_fear=first(SPIN_fear)) %>% distinct()
ggplot() +
           geom_errorbar(data=DtAvs2[DtAvs2$yesnochoice %in% "0",], aes(x=trialnum, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD, color=exptype), width=0.3) +
           geom_line(data=DtAvs2[DtAvs2$yesnochoice %in% "0",], aes(x=trialnum, y=Accuracy, color=exptype), linewidth=3, size=2) +
           geom_point(data=DtAvs2[DtAvs2$yesnochoice %in% "0",], aes(x=trialnum, y=Accuracy, fill=exptype), color="black", size=3, stroke=3, shape=21) +
           geom_errorbar(data=DtAvs2[DtAvs2$yesnochoice %in% "1",], aes(x=trialnum, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD, color=exptype), width=0.3, alpha=0.3) +
           geom_line(data=DtAvs2[DtAvs2$yesnochoice %in% "1",], aes(x=trialnum, y=Accuracy, color=exptype), linewidth=3, size=2, alpha=0.3) +
           geom_point(data=DtAvs2[DtAvs2$yesnochoice %in% "1",], aes(x=trialnum, y=Accuracy, fill=exptype), color="black", size=3, stroke=3, shape=21, alpha=0.3) +
           geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
          scale_fill_manual(aesthetics = c("color", "fill"), values = c("#14616B", "#ff2768", "#DD8016", "#b94b95")) +
  facet_wrap(vars(exptype)) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
labs(face="plain", family="arial",
     title = "Figure 4",
     x = "Trial", 
     y = "P(correct)")

#FIGURE 5

#FIGURE 5 Panel A

tertilesSPIN <- quantile(dataExp2[dataExp2$OA %in% "0",]$scaledSPIN, probs = c(1/3, 2/3, 3/3))
data3int <- dataExp2[dataExp2$OA %in% "0",]
data3int$pchoicebin <- 1 #again: illustration purposes
data3int[data3int$p_choice >= 0.8,]$pchoicebin <- 1
data3int[data3int$p_choice <= 0.8 & data3int$p_choice >= 0.6,]$pchoicebin <- 0.75
data3int[data3int$p_choice < 0.6 & data3int$p_choice >= 0.4,]$pchoicebin <- 0.50
data3int[data3int$p_choice <0.4,]$pchoicebin <- 0.25
data3int$SPINclass <- 1
data3int[data3int$scaledSPIN <= tertilesSPIN[[1]],]$SPINclass <- "Low"
data3int[data3int$scaledSPIN > tertilesSPIN[[1]] & data3int$scaledSPIN <= tertilesSPIN[[2]],]$SPINclass <- "Medium"
data3int[data3int$scaledSPIN > tertilesSPIN[[2]],]$SPINclass <- "High"
Plotchoicebin <- data3int %>% group_by(pchoicebin, SPINclass) %>% summarise(MeanAdvbin = mean(yesnochoicenum), 
                                                          SEPchoicebin =  sd(replicate(1000, mean(sample(yesnochoicenum, replace=T))))) %>% distinct()

DFforModel1 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy1", scaledSPIN=-0.5118994, scaledPTM=-0.02, SPINclass="Low")
DFforModel2 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy2", scaledSPIN=0.4782218, scaledPTM=-0.02, SPINclass="Medium")
DFforModel3 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy3", scaledSPIN=2.034127, scaledPTM=-0.02, SPINclass="High")


DFforModel<- rbind(DFforModel1,DFforModel2,DFforModel3)

#get marginals and ribbons (i.e. confidence interval for glmer fit)
#as prescribed by Ben Bolker in GLMM FAQ
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

ModA3bis <- glm(yesnochoicenum ~ p_choice*scaledSPIN,family =  binomial(link = "logit"),  data=dataExp2[dataExp2$OA %in% "0",]) 
DFforModel$modeledY <- predict(ModA3bis, type = 'response', newdata = DFforModel,  allow.new.levels = TRUE )
DFforModel$yesnochoicenum <- DFforModel$modeledY
des = model.matrix(formula(ModA3bis), DFforModel) #get the model matrix                                             
predvar = diag( des %*% vcov(ModA3bis) %*% t(des) ) # use matrix multiplication on the model matrix and variance-covariance matrix. Get diagonal of the resulting matrix
DFforModel$ConflowerRib = with(DFforModel, modeledY - 2*sqrt(predvar) ) #2 is fine, since SEM comes from 1.96
DFforModel$ConfupperRib = with(DFforModel, modeledY + 2*sqrt(predvar) )

#order factors for nicer plot
Plotchoicebin$SPINclass <- fct_relevel(Plotchoicebin$SPINclass, c("Low", "Medium","High")) #SPIN
DFforModel$SPINclass <- fct_relevel(DFforModel$SPINclass, c("Low", "Medium","High"))#SPIN

#and plot
ggplot() +
  geom_line(data=DFforModel, aes(x=p_choice, y=modeledY, color=SPINclass), linewidth=3, size=2, alpha=1) + #PTM plot
  geom_ribbon(data=DFforModel, aes(x=p_choice, ymin=ifelse(ConflowerRib<0, 0, ConflowerRib), ymax=ifelse(ConfupperRib>1,1,ConfupperRib), fill=SPINclass), linewidth=3, size=2, alpha=0.1) +
  geom_errorbar(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, ymax=MeanAdvbin + SEPchoicebin, ymin=MeanAdvbin - SEPchoicebin, color=SPINclass), linewidth=1, width=0.02, alpha=1) +
  geom_point(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, fill=SPINclass), color="black", size=3, stroke=2, shape=21, alpha=1) +
  geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
  scale_fill_manual(aesthetics = c("color","fill"), values = c("#2dbd20","#197c10","#0d3309")) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
  labs(face="plain", family="arial",
       title = "Figure 5 Panel A",
       x = "P(choice)", 
       y = "P(advice)")

#FIGURE 5 Panel B
tertilesSPIN <- quantile(dataExp2[dataExp2$OA %in% "0",]$scaledSPIN, probs = c(1/3, 2/3, 3/3))
data3int <- dataExp2[dataExp2$OA %in% "0",]

data3int$pchoicebin <- 1 
data3int[data3int$p_choice >= 0.8,]$pchoicebin <- 1
data3int[data3int$p_choice <= 0.8 & data3int$p_choice >= 0.6,]$pchoicebin <- 0.75
data3int[data3int$p_choice < 0.6 & data3int$p_choice >= 0.4,]$pchoicebin <- 0.50
data3int[data3int$p_choice <0.4,]$pchoicebin <- 0.25

data3int$SPINclass <- 1
data3int[data3int$scaledSPIN <= tertilesSPIN[[1]],]$SPINclass <- "Low"
data3int[data3int$scaledSPIN > tertilesSPIN[[1]] & data3int$scaledSPIN <= tertilesSPIN[[2]],]$SPINclass <- "Medium"
data3int[data3int$scaledSPIN > tertilesSPIN[[2]],]$SPINclass <- "High"

Plotchoicebin <- data3int %>% group_by(pchoicebin, SPINclass, exptype) %>% summarise(MeanAdvbin = mean(yesnochoicenum), 
                                                          SEPchoicebin =  sd(replicate(1000, mean(sample(yesnochoicenum, replace=T))))) %>% distinct()                                                           



DFforModel1 <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy1", scaledSPIN=-0.5118994, scaledPTM=-0.02, SPINclass="Low")
DFforModel2 <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy2", scaledSPIN=0.4782218, scaledPTM=-0.02, SPINclass="Medium")
DFforModel3 <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy3", scaledSPIN=2.034127, scaledPTM=-0.02, SPINclass="High")


DFforModel<- rbind(DFforModel1,DFforModel2,DFforModel3)
DFforModel$exptype <- factor(DFforModel$exptype)
DFforModel$exptype <- fct_relevel(DFforModel$exptype, "Money","No Cost (replication)","Reputation","Responsibility") 
ModA3bis <- glm(yesnochoicenum ~ p_choice*exptype*scaledSPIN,family =  binomial(link = "logit"),  data=dataExp2[dataExp2$OA %in% "0",])  #
DFforModel$modeledY <- predict(ModA3bis, type = 'response', newdata = DFforModel,  allow.new.levels = TRUE )
DFforModel$yesnochoicenum <- DFforModel$modeledY


des = model.matrix(formula(ModA3bis), DFforModel) #get ribbons                                             
predvar = diag( des %*% vcov(ModA3bis) %*% t(des) )
DFforModel$ConflowerRib = with(DFforModel, modeledY - 2*sqrt(predvar) )
DFforModel$ConfupperRib = with(DFforModel, modeledY + 2*sqrt(predvar) )
Plotchoicebin$SPINclass <- fct_relevel(Plotchoicebin$SPINclass, c("Low", "Medium","High")) #SPIN
DFforModel$SPINclass <- fct_relevel(DFforModel$SPINclass, c("Low", "Medium","High"))#SPIN

ggplot() +
  geom_line(data=DFforModel, aes(x=p_choice, y=modeledY, color=exptype), linewidth=3, size=2, alpha=1) +
  geom_ribbon(data=DFforModel, aes(x=p_choice, ymin=ifelse(ConflowerRib<0, 0, ConflowerRib), ymax=ifelse(ConfupperRib>1,1,ConfupperRib), fill=exptype), linewidth=3, size=2, alpha=0.1) +
  geom_errorbar(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, ymax=MeanAdvbin + SEPchoicebin, ymin=MeanAdvbin - SEPchoicebin, color=exptype), linewidth=1, width=0.02, alpha=1) +
  geom_point(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, fill=exptype), color="black", size=3, stroke=2, shape=21, alpha=1) +
  geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
  scale_fill_manual(aesthetics = c("color","fill"), values = c("#14616B", "#ff2768", "#DD8016", "#b94b95")) +
  facet_wrap(vars(factor(SPINclass)))+ #for the PTM plot
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
  labs(face="plain", family="arial",
       title = "Figure 5 Panel B",
       x = "P(choice)", 
       y = "P(advice)")

#FIGURE 5 Panel C

tertilesPTM <- quantile(dataExp2[dataExp2$OA %in% "0",]$scaledPTM, probs = c(1/3, 2/3, 3/3))
data3int <- dataExp2[dataExp2$OA %in% "0",]
data3int$pchoicebin <- 1 #again: illustration purposes
data3int[data3int$p_choice >= 0.8,]$pchoicebin <- 1
data3int[data3int$p_choice <= 0.8 & data3int$p_choice >= 0.6,]$pchoicebin <- 0.75
data3int[data3int$p_choice < 0.6 & data3int$p_choice >= 0.4,]$pchoicebin <- 0.50
data3int[data3int$p_choice <0.4,]$pchoicebin <- 0.25
data3int$PTMclass <- 1
data3int[data3int$scaledPTM <= tertilesPTM[[1]],]$PTMclass <- "Low"
data3int[data3int$scaledPTM > tertilesPTM[[1]] & data3int$scaledPTM < tertilesPTM[[2]],]$PTMclass <- "Medium"
data3int[data3int$scaledPTM >= tertilesPTM[[2]],]$PTMclass <- "High"
Plotchoicebin <- data3int %>% group_by(pchoicebin, PTMclass) %>% summarise(MeanAdvbin = mean(yesnochoicenum), 
                                                          SEPchoicebin =  sd(replicate(1000, mean(sample(yesnochoicenum, replace=T))))) %>% distinct()

DFforModel1 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy1", scaledSPIN=0.02746083, scaledPTM=-0.6962764, PTMclass="Low")
DFforModel2 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy2", scaledSPIN=0.02746083, scaledPTM=0.3075088, PTMclass="Medium")
DFforModel3 <- data.frame("p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy3", scaledSPIN=0.02746083, scaledPTM=2.315079, PTMclass="High")
DFforModel<- rbind(DFforModel1,DFforModel2,DFforModel3)

#get marginals and ribbons (i.e. confidence interval for glmer fit)
#as prescribed by Ben Bolker in GLMM FAQ
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

ModA3bis <- glm(yesnochoicenum ~ p_choice*scaledPTM,family =  binomial(link = "logit"),  data=dataExp2[dataExp2$OA %in% "0",]) 
DFforModel$modeledY <- predict(ModA3bis, type = 'response', newdata = DFforModel,  allow.new.levels = TRUE )
DFforModel$yesnochoicenum <- DFforModel$modeledY
des = model.matrix(formula(ModA3bis), DFforModel) #get the model matrix                                             
predvar = diag( des %*% vcov(ModA3bis) %*% t(des) ) # use matrix multiplication on the model matrix and variance-covariance matrix. Get diagonal of the resulting matrix
DFforModel$ConflowerRib = with(DFforModel, modeledY - 2*sqrt(predvar) ) #2 is fine, since SEM comes from 1.96
DFforModel$ConfupperRib = with(DFforModel, modeledY + 2*sqrt(predvar) )

#order factors for nicer plot
Plotchoicebin$PTMclass <- fct_relevel(Plotchoicebin$PTMclass, c("Low", "Medium","High")) #PTM
DFforModel$PTMclass <- fct_relevel(DFforModel$PTMclass, c("Low", "Medium","High"))#PTM

#and plot
ggplot() +
  geom_line(data=DFforModel, aes(x=p_choice, y=modeledY, color=PTMclass), linewidth=3, size=2, alpha=1) + #PTM plot
  geom_ribbon(data=DFforModel, aes(x=p_choice, ymin=ifelse(ConflowerRib<0, 0, ConflowerRib), ymax=ifelse(ConfupperRib>1,1,ConfupperRib), fill=PTMclass), linewidth=3, size=2, alpha=0.1) +
  geom_errorbar(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, ymax=MeanAdvbin + SEPchoicebin, ymin=MeanAdvbin - SEPchoicebin, color=PTMclass), linewidth=1, width=0.02, alpha=1) +
  geom_point(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, fill=PTMclass), color="black", size=3, stroke=2, shape=21, alpha=1) +
  geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
  scale_fill_manual(aesthetics = c("color","fill"), values = c("#2dbd20","#197c10","#0d3309")) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
  labs(face="plain", family="arial",
       title = "Figure 5 Panel C",
       x = "P(choice)", 
       y = "P(advice)")
#FIGURE 5 Panel D

tertilesPTM <- quantile(dataExp2[dataExp2$OA %in% "0",]$scaledPTM, probs = c(1/3, 2/3, 3/3))
data3int <- dataExp2[dataExp2$OA %in% "0",]

data3int$pchoicebin <- 1 
data3int[data3int$p_choice >= 0.8,]$pchoicebin <- 1
data3int[data3int$p_choice <= 0.8 & data3int$p_choice >= 0.6,]$pchoicebin <- 0.75
data3int[data3int$p_choice < 0.6 & data3int$p_choice >= 0.4,]$pchoicebin <- 0.50
data3int[data3int$p_choice <0.4,]$pchoicebin <- 0.25

data3int$PTMclass <- 1
data3int[data3int$scaledPTM <= tertilesPTM[[1]],]$PTMclass <- "Low"
data3int[data3int$scaledPTM > tertilesPTM[[1]] & data3int$scaledPTM <= tertilesPTM[[2]],]$PTMclass <- "Medium"
data3int[data3int$scaledPTM > tertilesPTM[[2]],]$PTMclass <- "High"

Plotchoicebin <- data3int %>% group_by(pchoicebin, PTMclass, exptype) %>% summarise(MeanAdvbin = mean(yesnochoicenum), 
                                                          SEPchoicebin =  sd(replicate(1000, mean(sample(yesnochoicenum, replace=T))))) %>% distinct()                                                           



DFforModel1 <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy1", scaledSPIN=0.02746083, scaledPTM=-0.6962764, PTMclass="Low")
DFforModel2 <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy2", scaledSPIN=0.02746083, scaledPTM=0.3075088, PTMclass="Medium")
DFforModel3 <- data.frame("exptype"=c(rep("No Cost (replication)",60),rep("Money",60), rep("Reputation",60), rep("Responsibility", 60)), "p_choice"=rep(seq(0.05,1,by=0.05),3),
                         "trialnum"=rep(1:20,3), "id"="dummy3", scaledSPIN=0.02746083, scaledPTM=2.315079, PTMclass="High")


DFforModel<- rbind(DFforModel1,DFforModel2,DFforModel3)
DFforModel$exptype <- factor(DFforModel$exptype)
DFforModel$exptype <- fct_relevel(DFforModel$exptype, "Money","No Cost (replication)","Reputation","Responsibility") 
ModA3bis <- glm(yesnochoicenum ~ p_choice*exptype*scaledPTM,family =  binomial(link = "logit"),  data=dataExp2[dataExp2$OA %in% "0",])  #
DFforModel$modeledY <- predict(ModA3bis, type = 'response', newdata = DFforModel,  allow.new.levels = TRUE )
DFforModel$yesnochoicenum <- DFforModel$modeledY


des = model.matrix(formula(ModA3bis), DFforModel) #get ribbons                                             
predvar = diag( des %*% vcov(ModA3bis) %*% t(des) )
DFforModel$ConflowerRib = with(DFforModel, modeledY - 2*sqrt(predvar) )
DFforModel$ConfupperRib = with(DFforModel, modeledY + 2*sqrt(predvar) )
Plotchoicebin$PTMclass <- fct_relevel(Plotchoicebin$PTMclass, c("Low", "Medium","High")) #PTM
DFforModel$PTMclass <- fct_relevel(DFforModel$PTMclass, c("Low", "Medium","High"))#PTM

ggplot() +
  geom_line(data=DFforModel, aes(x=p_choice, y=modeledY, color=exptype), linewidth=3, size=2, alpha=1) +
  geom_ribbon(data=DFforModel, aes(x=p_choice, ymin=ifelse(ConflowerRib<0, 0, ConflowerRib), ymax=ifelse(ConfupperRib>1,1,ConfupperRib), fill=exptype), linewidth=3, size=2, alpha=0.1) +
  geom_errorbar(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, ymax=MeanAdvbin + SEPchoicebin, ymin=MeanAdvbin - SEPchoicebin, color=exptype), linewidth=1, width=0.02, alpha=1) +
  geom_point(data=Plotchoicebin, aes(x=pchoicebin, y=MeanAdvbin, fill=exptype), color="black", size=3, stroke=2, shape=21, alpha=1) +
  geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
  scale_fill_manual(aesthetics = c("color","fill"), values = c("#14616B", "#ff2768", "#DD8016", "#b94b95")) +
  facet_wrap(vars(factor(PTMclass)))+ #for the PTM plot
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
  labs(face="plain", family="arial",
       title = "Figure 5 Panel D",
       x = "P(choice)", 
       y = "P(advice)")

##A convenient function to calculate effect sizes in glmer (i.e. R2 difference between full and null) based on models
# Calculate effect size (R² difference) for both main effects and interaction terms
# Main effects are removed along with their interaction terms, but interaction terms are only removed without affecting main effects
# as it makes little sense to make a model with an interaction term and not the main effects

library(MuMIn)

effect_size_R2 <- function(model_full) {
  
  # Get the fixed effect terms from the model
  terms <- attr(terms(model_full), "term.labels")
  
  # Initialize a list to store effect sizes
  effect_sizes <- list()
  
  # Calculate marginal R² for the full model
  r2_full <- r.squaredGLMM(model_full)[1]  # Marginal R2
  
  # Iterate over each term and calculate R2 difference
  for (term in terms) {
    
    if (grepl(":", term)) {
      # Interaction term: remove only the interaction term, but keep the main effects
      formula_reduced <- as.formula(paste(". ~ . -", term))
      
    } else {
      # Main effect: remove the main effect and any related interaction terms
      interaction_terms <- grep(term, terms, value = TRUE)
      formula_reduced <- as.formula(paste(". ~ . -", paste(interaction_terms, collapse = " - ")))
    }
    
    # Fit the reduced model
    model_reduced <- update(model_full, formula_reduced)
    
    # Calculate marginal R² for the reduced model
    r2_reduced <- r.squaredGLMM(model_reduced)[1]
    
    # Calculate the effect size (R2 difference)
    effect_size <- r2_full - r2_reduced
    
    # Store the result in the list
    effect_sizes[[term]] <- effect_size
  }
  
  # Convert the result to a data frame for easier viewing
  effect_sizes_df <- as.data.frame(do.call(rbind, effect_sizes))
  colnames(effect_sizes_df) <- "Effect Size (R² Difference)"
  
  return(effect_sizes_df)
}

# Call the function to calculate effect sizes for both main effects and interaction terms
effect_size_R2(ModPreg7)
effect_size_R2(ModPreg3)

#Get confidence intervals for the coefficients for the models
cl <- makeCluster(30) #or whatever number you want depending on your hardware
conf_intervals <- confint(ModPreg7, method = "boot", parallel = "snow", ncpus = 30, cl = cl) #ncpus = your cores for the task
