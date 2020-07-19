### This code requires use of the PECARN public use dataset for the study
### "Identifying Children at Very Low Risk of Clinically Important Blunt Abdominal Injuries"
### available at http://pecarn.org/studyDatasets/StudyDetails?studyID=8 after reviewing
### and accepting terms of agreement


#Load libraries; install pacakges as needed
library(gmodels)
library(Hmisc)
library(missForest)
library(splitstackshape)
library(SuperLearner)
library(ranger)
library(xgboost)
library(MASS)
library(arm)
library(kernlab)
library(rpart)
library(earth)
library(glmnet)
library(devtools)
library(cutpointr)
library(epiR)
library(pROC)
library(ROCR)

setwd("~/Research/IAI/IAIP/Datasets/CSV/")  ##set working directory

temp = list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))),
         read.csv), envir = .GlobalEnv)

IAI <- merge(x = demographics, y = form1, by = c("SubjectID"), all.x = TRUE)
##19 predictors are used in the study by Pennell
# Age                                         ageinyrs
# Abdominal distension                        AbdDistention
# Sex                                         sex
# Absent bowel sounds                         BowelSounds
# Heart rate                                  InitHeartRate
# Abdominal tenderness                        AbdomenTender
# Respiratory rate                            InitRespRange
# Peritoneal signs                            PeritonIrrit
# Glasgow Coma Scale                          GCSScore
# Visible thoracic trauma                     ThoracicTrauma
# Dyspnea                                     ShortBreath
# Flank pain                                  FlankTender
# Emesis                                      VomitWretch
# Pelvic pain                                 PelvicTender
# Visible abdominal trauma                    AbdTrauma
# Unstable pelvis                             PevisUnstable
# Seatbelt sign                               SeatBeltSign
# Occult rectal blood                         RectalBlood
# Abdominal pain                              AbdomenPain

##perform some initial cleanup of predictors. Convert sex to a number, GCS to a dichotomous variable,
##convert all categorical variables to factors and all continuous variables to numeric
IAI$GCSScore <- ifelse(IAI$GCSScore<15,1,0)
IAI$sex <- ifelse(IAI$SEX=="F",2,1)
IAI$AbdDistention[IAI$AbdDistention==4] <- NA
IAI$AbdDistention[IAI$BowelSounds==4] <- NA
IAI$PeritonIrrit[IAI$PeritonIrrit==4] <- NA
IAI$ThoracicTrauma[IAI$ThoracicTrauma==3] <- NA
IAI$ShortBreath[IAI$ShortBreath==4] <- NA
IAI$FlankTender[IAI$FlankTender==4] <- NA
IAI$VomitWretch[IAI$VomitWretch==4] <- NA
IAI$PelvicTender[IAI$PelvicTender==4] <- NA
IAI$AbdTrauma[IAI$AbdTrauma==4] <- NA
IAI$PelvisUnstable[IAI$PelvisUnstable==4] <- NA
IAI$RectalBlood[IAI$RectalBlood==4] <- NA
IAI$AbdomenPain[IAI$AbdomenPain==4] <- NA
IAI$AbdomenTender[IAI$AbdomenTender==4] <- NA


IAI$ageinyrs <-as.numeric(IAI$ageinyrs)
IAI$sex <-as.factor(IAI$sex)
IAI$BowelSounds <-as.factor(IAI$BowelSounds)
IAI$InitHeartRate <-as.numeric(IAI$InitHeartRate)
IAI$InitRespRange <-as.numeric(IAI$InitRespRange)
IAI$AbdDistention <-as.factor(IAI$AbdDistention)
IAI$PeritonIrrit <-as.factor(IAI$PeritonIrrit)
IAI$GCSScore <-as.factor(IAI$GCSScore)
IAI$ThoracicTrauma <-as.factor(IAI$ThoracicTrauma)
IAI$ShortBreath <-as.factor(IAI$ShortBreath)
IAI$FlankTender <-as.factor(IAI$FlankTender)
IAI$VomitWretch <-as.factor(IAI$VomitWretch)
IAI$PelvicTender <-as.factor(IAI$PelvicTender)
IAI$AbdTrauma <-as.factor(IAI$AbdTrauma)
IAI$PelvisUnstable <-as.factor(IAI$PelvisUnstable)
IAI$SeatBeltSign <-as.factor(IAI$SeatBeltSign)
IAI$RectalBlood <-as.factor(IAI$RectalBlood)
IAI$AbdomenPain <-as.factor(IAI$AbdomenPain)
IAI$AbdomenTender <-as.factor(IAI$AbdomenTender)


## Next let's work on the outcomes.
# IAI                                                                form6a$IAIinED1
# AND any of the following:

# 1) death,                                                          form6a, DeathCause =1 or =2
# 2) required therapeutic angiography                                form4bother_abdangio$AbdAngioVessel and form4bother_abdangio$PelAngioVessel                           
# 3) laparotomy                                                      form6c$IntervenDurLap
# 4) blood transfusion,                                              BldTransfusion$Bldtransfusion
# 5) admission to the hospital for two or more nights for IVF        form6b$IVfluids

## We will create an additional dataframe called "Outcomes" to temporarily identify the outcome.

form4bother_abdangio <- subset(form4bother_abdangio,AbdAngioVessel==1,select=c(subjectid,AbdAngioVessel))
form4bother_abdangio <- form4bother_abdangio[!duplicated(form4bother_abdangio$subjectid), ]

form4bother_pelangio <- subset(form4bother_pelangio,PelAngioVessel==1,select=c(subjectid,PelAngioVessel))
form4bother_pelangio <- form4bother_pelangio[!duplicated(form4bother_pelangio$subjectid), ]
##there are no patients with therapeutic pelvic angiography, so this is ignored in further steps

form6c$laparotomy[form6c$IntervenDurLap==1] <- 1
form6c <- subset(form6c,laparotomy==1,select=c(subjectid,laparotomy))
form6c <- form6c[!duplicated(form6c$subjectid), ]  

BldTransfusion <- subset(form6b,BldTransfusion==1,select=c(SubjectID,BldTransfusion))
BldTransfusion <- BldTransfusion[!duplicated(BldTransfusion$SubjectID), ]  

IVFluids <- subset(form6b,IVFluids==1,select=c(SubjectID,IVFluids))
IVFluids <- IVFluids[!duplicated(IVFluids$SubjectID), ]  
form6a$AbdDeath[form6a$DeathCause==1 | form6a$DeathCause==2] <- 1
Outcome <- merge(x = form6a, y = form4bother_abdangio, by.x = c("subjectid"), by.y = c("subjectid"), all.x = TRUE)  ##angio
Outcome <- merge(x = Outcome, y = form6c, by.x = c("subjectid"), by.y = c("subjectid"), all.x = TRUE)  ##laparotomy
Outcome <- merge(x = Outcome, y = BldTransfusion, by.x = c("subjectid"), by.y = c("SubjectID"), all.x = TRUE)  ##transfusion
Outcome <- merge(x = Outcome, y = IVFluids, by.x = c("subjectid"), by.y = c("SubjectID"), all.x = TRUE)  ##adission for ivf
Outcome <- subset(Outcome,select=c(subjectid,IAIinED1,AbdDeath,AbdAngioVessel,laparotomy,BldTransfusion,IVFluids))
Outcome[is.na(Outcome)] <- 0


Outcome$seq <- Outcome$AbdDeath + Outcome$AbdAngioVessel + Outcome$laparotomy + Outcome$BldTransfusion + Outcome$IVFluids
Outcome$outcome <- ifelse(Outcome$IAIinED1==1 & Outcome$seq>0,1,0)
Outcome$outcome <-as.factor(Outcome$outcome)

describe(as.factor(Outcome$IAIinED1)) #761 with IAI (same as parent study)
CrossTable(Outcome$outcome,Outcome$IAIinED1) #203 patients with IAI-I (same as parent study)

IAI <- merge(x = IAI, y = Outcome, by.x = c("SubjectID"), by.y = c("subjectid"), all.x = TRUE)  ##merge predictor to outcome

IAI_learner <- subset(IAI,select=c(ageinyrs,AbdDistention,sex,BowelSounds,AbdomenTender,InitHeartRate,InitRespRange,PeritonIrrit,GCSScore,ThoracicTrauma,
                                   ShortBreath,FlankTender,VomitWretch,PelvicTender,AbdTrauma,PelvisUnstable,SeatBeltSign,RectalBlood,AbdomenPain,outcome
))


summary(IAI_learner)  ##investigate the variables

set.seed(108)
missF <- missForest(IAI_learner)
imputedData <- missF$ximp
imputedData <- imputedData[complete.cases(imputedData), ]


## take a stratified outcome to ensure equal proportions of outcome in each group.
set.seed(108)
out <- stratified(imputedData, c("outcome"), 0.75,bothSets = T,keep.rownames=F)
Model_Train <- as.data.frame(out$SAMP1)
Model_Test <- as.data.frame(out$SAMP2)

#compare the two groups by chi-squared tests or wilcoxon rank-sum tests
a <- Model_Test
a$group <- "Test"
b<- Model_Train
b$group <- "Train"
compare <- rbind(a,b)
CrossTable(compare$AbdDistention,compare$group, chisq = T)
CrossTable(compare$BowelSounds,compare$group, chisq = T)
CrossTable(compare$AbdomenTender,compare$group, chisq = T)
CrossTable(compare$PeritonIrrit,compare$group, chisq = T)
CrossTable(compare$AbdomenTender,compare$group, chisq = T)
CrossTable(compare$GCSScore,compare$group, chisq = T)
CrossTable(compare$ThoracicTrauma,compare$group, chisq = T)
CrossTable(compare$ShortBreath,compare$group, chisq = T)
CrossTable(compare$FlankTender,compare$group, chisq = T)
CrossTable(compare$VomitWretch,compare$group, chisq = T)
CrossTable(compare$PelvicTender,compare$group, chisq = T)
CrossTable(compare$AbdTrauma,compare$group, chisq = T)
CrossTable(compare$PelvisUnstable,compare$group, chisq = T)
CrossTable(compare$SeatBeltSign,compare$group, chisq = T)
CrossTable(compare$RectalBlood,compare$group, chisq = T)
CrossTable(compare$AbdomenPain,compare$group, chisq = T)
wilcox.test(ageinyrs ~ group, data = compare, exact = FALSE)
wilcox.test(InitHeartRate ~ group, data = compare, exact = FALSE)
wilcox.test(InitRespRange ~ group, data = compare, exact = FALSE)

## prepare the training dataset for SuperLearner
base_outcomes <- as.numeric(Model_Train$outcome)-1
predictors <- Model_Train[,1:19]

et.seed(108)
superlearner <- SuperLearner(Y = base_outcomes, X = predictors , family = binomial(), verbose=TRUE,method="method.NNloglik",
                                   SL.library = c("SL.xgboost","SL.ranger","SL.stepAIC","SL.bayesglm","SL.ksvm","SL.rpart","SL.earth","SL.glmnet","SL.glm")) 


Model_Test$predictions_base <- predict.SuperLearner(superlearner, newdata=subset(Model_Test,select=c(1:19)))$pred
Model_Train$predictions_base <- predict.SuperLearner(superlearner, newdata=subset(Model_Train,select=c(1:19)))$pred

Cutpoints <- cutpointr(Model_Test, predictions_base, outcome, 
                       direction = ">=", pos_class = 1,
                       neg_class = 0,method=minimize_metric,
                       metric = misclassification_cost, cost_fp=1, cost_fn=250)

Z<-epi.tests(as.table(matrix(c(length(which(Model_Test$outcome==1 & Model_Test$predictions_base>=as.numeric(Cutpoints[2]))),
                               length(which(Model_Test$outcome==0 & Model_Test$predictions_base>=as.numeric(Cutpoints[2]))),
                               length(which(Model_Test$outcome==1 & Model_Test$predictions_base<as.numeric(Cutpoints[2]))),
                               length(which(Model_Test$outcome==0 & Model_Test$predictions_base<as.numeric(Cutpoints[2])))),
                             nrow = 2, byrow = TRUE)))

paste("Sensitivity",round(Z$elements$se*100,1),"(",round(Z$elements$se.low*100,1),round(Z$elements$se.up*100,1),")")
paste("Specificity",round(Z$elements$sp*100,1),"(",round(Z$elements$sp.low*100,1),round(Z$elements$sp.up*100,1),")")
paste("PPV",round(Z$elements$ppv*100,1),"(",round(Z$elements$ppv.low*100,1),round(Z$elements$ppv.up*100,1),")")
paste("NPV",round(Z$elements$npv*100,1),"(",round(Z$elements$npv.low*100,1),round(Z$elements$npv.up*100,1),")")
roc(Model_Test$outcome, Model_Test$predictions_base, percent=TRUE, ci=TRUE)

#all predictions
allmodels <- data.frame(predict.SuperLearner(superlearner, newdata=subset(Model_Test,select=c(1:19)))$library.predict)
allmodels_train <- data.frame(predict.SuperLearner(superlearner, newdata=subset(Model_Train,select=c(1:19)))$library.predict)




par(pty="s",mgp=c(2,1,0))
pred <- prediction( Model_Test$predictions_base, Model_Test$outcome )
perf <- performance( pred, "tpr", "fpr" )
pred2 <- prediction( allmodels$SL.xgboost_All, Model_Test$outcome )
perf2 <- performance( pred2, "tpr", "fpr" )
pred3 <- prediction( allmodels$SL.ranger_All, Model_Test$outcome )
perf3 <- performance( pred3, "tpr", "fpr" )
pred4 <- prediction( allmodels$SL.stepAIC_All, Model_Test$outcome )
perf4 <- performance( pred4, "tpr", "fpr" )
pred5 <- prediction( allmodels$SL.bayesglm_All, Model_Test$outcome )
perf5 <- performance( pred5, "tpr", "fpr" )
pred6 <- prediction( allmodels$SL.ksvm_All, Model_Test$outcome )
perf6 <- performance( pred6, "tpr", "fpr" )
pred7 <- prediction( allmodels$SL.rpart_All, Model_Test$outcome )
perf7 <- performance( pred7, "tpr", "fpr" )
pred8 <- prediction( allmodels$SL.earth_All, Model_Test$outcome )
perf8 <- performance( pred8, "tpr", "fpr" )
pred9 <- prediction( allmodels$SL.glmnet_All, Model_Test$outcome )
perf9 <- performance( pred9, "tpr", "fpr" )
pred10 <- prediction( allmodels$SL.glm_All, Model_Test$outcome )
perf10 <- performance( pred10, "tpr", "fpr" )

discrete<- "black"
xgboost<- "#cb9b44"
ranger<- "#b459c1" 
stepAIC<- "#74b54a"
bayesglm<- "#d03f68"
ksvm<- "#4bb194"
rpart<- "#c8603f"
earth<- "#6e7ecc"
glmnet<- "#737e39"
glm<- "#c16b94"
plot( perf, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=discrete,xlab="1 - Specificity", ylab="Sensitivity")
plot( perf2, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=xgboost,add=T,lty=3)
plot( perf3, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=ranger,add=T,lty=4)
plot( perf4, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=stepAIC,add=T,lty=5)
plot( perf5, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=bayesglm,add=T,lty=6)
plot( perf6, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=ksvm,add=T,lty=7)
plot( perf7, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=rpart,add=T,lty=8)
plot( perf8, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=earth,add=T,lty=9)
plot( perf9, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=glmnet,add=T,lty=10)
plot( perf10, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=glm,add=T,lty=11)
abline(coef = c(0,1),col="#696969",lty=2)




par(pty="s",mgp=c(2,1,0))
pred <- prediction( Model_Train$predictions_base, Model_Train$outcome )
perf <- performance( pred, "tpr", "fpr" )
pred2 <- prediction( allmodels_train$SL.xgboost_All, Model_Train$outcome )
perf2 <- performance( pred2, "tpr", "fpr" )
pred3 <- prediction( allmodels_train$SL.ranger_All, Model_Train$outcome )
perf3 <- performance( pred3, "tpr", "fpr" )
pred4 <- prediction( allmodels_train$SL.stepAIC_All, Model_Train$outcome )
perf4 <- performance( pred4, "tpr", "fpr" )
pred5 <- prediction( allmodels_train$SL.bayesglm_All, Model_Train$outcome )
perf5 <- performance( pred5, "tpr", "fpr" )
pred6 <- prediction( allmodels_train$SL.ksvm_All, Model_Train$outcome )
perf6 <- performance( pred6, "tpr", "fpr" )
pred7 <- prediction( allmodels_train$SL.rpart_All, Model_Train$outcome )
perf7 <- performance( pred7, "tpr", "fpr" )
pred8 <- prediction( allmodels_train$SL.earth_All, Model_Train$outcome )
perf8 <- performance( pred8, "tpr", "fpr" )
pred9 <- prediction( allmodels_train$SL.glmnet_All, Model_Train$outcome )
perf9 <- performance( pred9, "tpr", "fpr" )
pred10 <- prediction( allmodels_train$SL.glm_All, Model_Train$outcome )
perf10 <- performance( pred10, "tpr", "fpr" )

discrete<- "black"
xgboost<- "#cb9b44"
ranger<- "#b459c1" 
stepAIC<- "#74b54a"
bayesglm<- "#d03f68"
ksvm<- "#4bb194"
rpart<- "#c8603f"
earth<- "#6e7ecc"
glmnet<- "#737e39"
glm<- "#c16b94"

plot( perf, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=discrete,xlab="1 - Specificity", ylab="Sensitivity")
plot( perf2, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=xgboost,add=T,lty=3)
plot( perf3, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=ranger,add=T,lty=4)
plot( perf4, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=stepAIC,add=T,lty=5)
plot( perf5, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=bayesglm,add=T,lty=6)
plot( perf6, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=ksvm,add=T,lty=7)
plot( perf7, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=rpart,add=T,lty=8)
plot( perf8, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=earth,add=T,lty=9)
plot( perf9, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=glmnet,add=T,lty=10)
plot( perf10, colorize = F,xaxs="i",yaxs="i",lwd=1.5,col=glm,add=T,lty=11)
abline(coef = c(0,1),col="#696969",lty=2)
