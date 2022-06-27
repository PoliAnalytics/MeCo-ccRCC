# install.packages("pROC")
# install.packages("MLmetrics")

library(MASS)
library(pROC)
library(MLmetrics)


# set up a dataset with the variables of interest:
surv_data <- data.frame(pheno_data$submitter_id.samples) 
#age
surv_data$age <- pheno_data$age_at_initial_pathologic_diagnosis
#time to event
surv_data$time <- NA
time_to_death <- pheno_data$days_to_death.demographic
time_to_last_FU <- pheno_data$days_to_last_follow_up.diagnoses

for (i in 1:length(sequenced_patients)){
  if (is.na(time_to_death[i]))
    surv_data$time[i] <- time_to_last_FU[i]
  else
    surv_data$time[i] <- time_to_death[i]
}

#status
surv_data$status <- as.numeric(as.factor(pheno_data$vital_status.demographic)) # 1 = patient alive, 2 = patient dead
#metastasis
surv_data$metastasis <- NA
meta <- pheno_data$additional_surgery_metastatic_procedure

for (i in 1:length(sequenced_patients)){
  if (meta[i]=="")
    surv_data$metastasis[i] <- 'NO'
  else
    surv_data$metastasis[i] <- 'YES'
}
surv_data$metastasis <- as.factor(surv_data$metastasis)
#sex
surv_data$gender <- as.factor(pheno_data$gender.demographic)
#tumor stage
surv_data$stage <- as.factor(pheno_data$tumor_stage.diagnoses)
# MeCo scores
surv_data$MeCo <- MecoScore$MeSc
surv_data$MeCo_Ant <- MecoScoreAnt$MeSc
surv_data$MeCo_Ch <- MecoScoreCh$MeSc
surv_data$MeCo_ECM <- MecoScoreECM$MeSc
surv_data$MeCo_Inf <- MecoScoreInf$MeSc
surv_data$MeCo_Pro <- MecoScorePro$MeSc
# laterality
surv_data$laterality <- as.factor(pheno_data$laterality)

# since in logistic model 0<=y<=1:
surv_data$status[which(surv_data$status==1)] <- 0 # alive
surv_data$status[which(surv_data$status==2)] <- 1 # dead

# fit the logistic model by gender
logistic_gender<- glm(status~ gender, family=binomial(link=logit), data = surv_data)
summary(logistic_gender) 
# being male has no effect on the probability of death 

# fit the logistic model by laterality
logistic_lat<- glm(status~ laterality, family=binomial(link=logit), data = surv_data)
summary(logistic_lat) 
# laterality has no effect on the probability of death

# fit the logistic model by stage
logistic_stage<- glm(status~ stage, family=binomial(link=logit), data = surv_data)
summary(logistic_stage) 
exp(logistic_stage$coefficients)
# being of stage 2 has no effect.
# being of stage 3 is a significant risk factor (p-value = 8.85e-07) that increases the probability of death
# of 3.37 times with respect to being of stage 1.
# being of stage 4 is a significant risk factor (p-value = <2e-16) that increases the probability of death
# of 21 times with respect to being of stage 1.

# fit the logistic model by metastasis
logistic_metastasis <- glm(status~ metastasis, family=binomial(link=logit), data = surv_data)
summary(logistic_metastasis) 
exp(logistic_metastasis$coefficients)
# Having metastasis is a significant (p-value = 3.76e-09) risk factor that increases the probability fo death of
# 3.27 times with respect to not having them

# fit the logistic model by age
logistic_age <- glm(status~ age, family=binomial(link=logit), data = surv_data)
summary(logistic_age) 
exp(logistic_age$coefficients)
# age is a significant risk factor (p val 3.42e-06): unitary increase of age increases the probability of death
# of 3.9%

# fit the logistic model by age and sex
log_age_sex <- glm(status~ age+gender, family=binomial(link=logit), data = surv_data)
summary(log_age_sex) 
exp(log_age_sex$coefficients)
# being male remains with no effect and age still has the same effect as considering it alone


# fit the logistic model by MeCo (non refined)
logistic_MeCo <- glm(status~ MeCo, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCo) 
# Meco has no effect: p-va 0.07

# fit the logistic model by MeCo and age
log_MeCo_age <- glm(status~ MeCo+age, family=binomial(link=logit), data = surv_data)
summary(log_MeCo_age) 
exp(log_MeCo_age$coefficients)
# meco has no effect (0.09) while age has the same as alone

# fit the logistic model by MeCo_Ant 
logistic_MeCoAnt <- glm(status~ MeCo_Ant, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoAnt) 
# MecoAnt has no effect: p val 0.42

# fit the logistic model by MeCo_Ch
logistic_MeCoCh <- glm(status~ MeCo_Ch, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoCh) 
# MecoCh has no effect: p val 0.61

# fit the logistic model by MeCo_ECM
logistic_MeCoECM <- glm(status~ MeCo_ECM, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoECM) 
# MecoECM has no effect: p val 0.14

# fit the logistic model by MeCo_Inf
logistic_MeCoInf <- glm(status~ MeCo_Inf, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoInf) 
# mecoInf has no effect on the probability of survival (p val 0.97)

# fit the logistic model by MeCo_Pro
logistic_MeCoPro <- glm(status~ MeCo_Pro, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoPro) 
# MecoPro has no effect: p val 0.08

# fit the logistic model by age, stage, metastasis
log_age_st_meta <- glm(status~ age+stage+metastasis, family=binomial(link=logit), data = surv_data)
summary(log_age_st_meta) 
exp(log_age_st_meta$coefficients)
# being  of stage 2 wrt to stage 1 has no effect.
# age is a significant risk factor (p val 8.94e-07): its unitary increase increases the rish of death od 4.8%
# being of stage 3 is a significant risk factor (p-val 0.000101) that increases the probability of death of 2.74 times;
# being of stage 4 is a significant risk factor (p-val < 2e-16) that increases the probability of death of 19.6 times;
# having metastasis is a significant risk factor (p val 0.007165) that increases the probability of death of 90%.

anova(log_age_st_meta, logistic_age,test="LRT")
# low p-val -> full model is more informative

# fit the logistic model by age, stage, metastasis and MeCO (non refined)
log_meco <- glm(status~ age+stage+metastasis+MeCo, family=binomial(link=logit), data = surv_data)
summary(log_meco) 
# age, stage 3, stage 4, meta remain meaningful while meco not 

# fit the logistic model by age, stage, metastasis and all MeCos
log_meco2 <- glm(status~ age+stage+metastasis+MeCo_Ant+MeCo_Ch+MeCo_Inf+MeCo_Pro+MeCo_ECM, family=binomial(link=logit), data = surv_data)
summary(log_meco2) 
# none meco is significant 

# fit the logistic model by age, stage, metastasis and MeCo_Pro
log_mecoPro <- glm(status~ age+stage+metastasis+MeCo_Pro, family=binomial(link=logit), data = surv_data)
summary(log_mecoPro) 
# mecoPro p-value 0.22 (peggio che considerato da solo)

# fit the logistic model only with refined mecos
log_meco3 <- glm(status~ MeCo_Ant+MeCo_Ch+MeCo_Inf+MeCo_Pro+MeCo_ECM, family=binomial(link=logit), data = surv_data)
summary(log_meco3) 
# none 

# fit the logistic model only with meco ant and meco pro
log_meco4 <- glm(status~ MeCo_Ant+MeCo_Pro, family=binomial(link=logit), data = surv_data)
summary(log_meco4) 
exp(log_meco4$coefficients)
# unitary increase of meco pro is s risk factor (p val 0.048) that increases the prob of death of 93%

# fit the logistic model to predict metastasis
log_meta <- glm(metastasis~ age+stage+MeCo_Pro, family=binomial(link=logit), data = surv_data)
summary(log_meta) 
# only stage 3 and 4 are significant risk factors

# fit the logistic model to predict metastasis
log_meta2 <- glm(metastasis~ age+stage+MeCo_Pro+MeCo_Ant, family=binomial(link=logit), data = surv_data)
summary(log_meta2) 
# only stage 3 and 4 are significant risk factors

p_threshold = 0.5

Y.hat <- ifelse(log_age_st_meta$fitted.values<p_threshold, 0, 1) 
Y.hat

# Confusion Matrix
table(Predicted = Y.hat, Observed = surv_data$status)
N <- nrow(surv_data)

# Compute the misclassification rate
errors <- (Y.hat != surv_data$status)
MIS_Rate  <- sum(errors)/N
MIS_Rate 
#21.5% of the patients are wrongly collocated

Specificity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 92.7% specificity

Sensitivity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 49% sensitivity

# use the empirical threshold: p=sum(observedy)/N
p_threshold = sum(surv_data$status)/N
# 0.32

Y.hat <- ifelse(log_age_st_meta$fitted.values<p_threshold, 0, 1) 
Y.hat

table(Predicted = Y.hat, Observed = surv_data$status)

errors <- (Y.hat != surv_data$status)
MIS_Rate  <- sum(errors)/N
MIS_Rate 
# 24% misclassified

Specificity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 79.1% specificity

Sensitivity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 69.2% sensitivity

ROC_curve <- roc(response = surv_data$status, predictor = log_age_st_meta$fitted.values,
                 levels = c('0','1'),
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")
auc(ROC_curve)
# 0.8106

coords(ROC_curve, x="best", transpose = TRUE)
# threshold specificity sensitivity 
# 0.3685864   0.8384401   0.6627907 

p_threshold <- 0.3685864

#  computation of sensitivity, specificity and AUC performances using the empirical threshold and 10-fold cross-validation.
K = 10
folds <- cut(seq(1,N), breaks=K ,labels=FALSE)#Create K equally size folds (if possible)
set.seed(1234)
folds <- sample(folds)#Randomly shuffle the observations
table(folds)

sensitivity<-NULL
specificity<-NULL
AUC<-NULL
for(k in 1:10){
  train.data <- surv_data[which(folds!=k),]
  test.data <- surv_data[which(folds==k),]

  log_age_st_meta.k <- glm(status~ age+stage+metastasis , family=binomial(link=logit), data = train.data)
  p.hat.k <- predict( log_age_st_meta.k, newdata = data.frame(test.data), type='response' )
  Y.hat.k <- Y.hat <- ifelse(p.hat.k <p_threshold, 0, 1) 
  
  sensitivity <- c(sensitivity,
                   Sensitivity(y_true =  test.data$status, y_pred = Y.hat.k, positive = 1)
  )
  specificity <- c(specificity,
                   Specificity(y_true =  test.data$status, y_pred = Y.hat.k, positive = 1)
  )  
  
  AUC <- c(AUC,
           roc(response =  test.data$status, predictor = Y.hat.k,
               levels = c('0','1'),
               smooth=FALSE, plot=F, print.auc=F)$auc
  )
  
}

mean(sensitivity)
# 0.64

mean(specificity)
# 0.83

sd(specificity)
# 0.098

mean(AUC)
# 0.736





# build a logistic regression model for a multiclass problem

library(nnet)

# predict the stage from status, metastasis, age and meco_pro
mult_logistic_stage <- multinom(stage ~ status+metastasis+age+MeCo_Pro, data = surv_data)
z <- summary(mult_logistic_stage)$coefficients/summary(mult_logistic_stage)$standard.errors
z

p <- (1 - pnorm(abs(z), 0, 1)) * 2
p
exp(coef(mult_logistic_stage))
# for stage 2 the only feature with an effect is metastasis (p val: 6.400487e-02): having the metastasis increases the probability
# of being of stage 2, wrt stage 1, of 88%
# for stage 3 both status and metastasis have an effect: being dead increases the probability
# of being of stage 3, wrt stage 1, of 2.7 times. having the metastasis increases the probability
# of being of stage 3, wrt stage 1, of 2.5 times.
# for stage 4 both status and metastasis have an effect: being dead increases the probability
# of being of stage 4, wrt stage 1, of 20 times. having the metastasis increases the probability
# of being of stage 4, wrt stage 1, of 5.2 times.

# predict the stage from status, metastasis, age and meco_pro
mult_logistic_stage2 <- multinom(stage ~ MeCo_Pro, data = surv_data)
z <- summary(mult_logistic_stage2)$coefficients/summary(mult_logistic_stage2)$standard.errors
z

p <- (1 - pnorm(abs(z), 0, 1)) * 2
p
# meco pro not significant
