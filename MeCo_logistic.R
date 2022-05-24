# performing Logistic Regression
# having found that the MeCo_Pro score has some predictive power for survival, test also whether it is capable of predicting:
# status(dead or alive), metastasis(presence or absence), stage(i, ii, iii or iv)

library(MASS)
library(pROC)
library(MLmetrics)

# test for status
# dead patients as 1, alive as 0
surv_data$status2 <- NA
for (i in 1:length(surv_data$status)){
  if (surv_data$status[i]==1)
    surv_data$status2[i]=0
  else 
    surv_data$status2[i]=1
} 

mod.status <- glm(status2 ~ MeCo_Pro+MeCo_Ant+MeCo_Inf+MeCo_ECM+MeCo_Ch, family=binomial(link=logit), data = surv_data)    
summary(mod.status)      
p_threshold = 0.3
Y.hat <- ifelse(mod.status$fitted.values<p_threshold, 0, 1) 
confusion.matrix <- table(Predicted=Y.hat, Observed=surv_data$status2)
confusion.matrix
Sensitivity(y_true = surv_data$status2, y_pred = Y.hat, positive = 1)
Specificity(y_true = surv_data$status2, y_pred = Y.hat, positive = 1)
ROC_curve <- roc(response = surv_data$status2, predictor = mod.status$fitted.values,
                 levels = c('0','1'),
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")
coords(ROC_curve, x="best", transpose = TRUE)
# the model including all the refined mecos is not performing much better than a random classifier for the status

# test for metastasis
surv_data$metastasis2 <- surv_data$metastasis
levels(surv_data$metastasis2) <- c(0,1)

mod.meta <- glm(metastasis2 ~ MeCo_Pro+MeCo_Ant+MeCo_Inf+MeCo_ECM+MeCo_Ch, family=binomial(link=logit), data = surv_data)    
summary(mod.meta)   
p_threshold = 0.3
Y.hat <- ifelse(mod.meta$fitted.values<p_threshold, 0, 1) 
confusion.matrix <- table(Predicted=Y.hat, Observed=surv_data$metastasis2)
confusion.matrix
Sensitivity(y_true = surv_data$metastasis2, y_pred = Y.hat, positive = 1)
Specificity(y_true = surv_data$metastasis2, y_pred = Y.hat, positive = 1)
ROC_curve <- roc(response = surv_data$metastasis2, predictor = mod.meta$fitted.values,
                 levels = c('0','1'),
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")
coords(ROC_curve, x="best", transpose = TRUE)
# the model including all the refined mecos is not performing much better than a random classifier for the metastasis
# NOTE: in this case as what we are predicting is no longer the survival of the patients, the pvalues have changed
# in the overall model the most significant pvalue is associated with MeCo_ECM 

# test for tumor stages
library(nnet)
mod.stage <- multinom(stage ~ MeCo_Pro+MeCo_Ant+MeCo_Inf+MeCo_ECM+MeCo_Ch, data = surv_data)
summary(mod.stage)

# test tumor stage as either i-ii or iii-iv
surv_data$stage2 <- NA
for (i in 1:length(surv_data$stage)){
  if (surv_data$stage[i]=='stage i' | surv_data$stage[i]=='stage ii')
    surv_data$stage2[i]=0
  else 
    surv_data$stage2[i]=1
}

mod.stagecat <- glm(stage2 ~ MeCo_Pro+MeCo_Ant+MeCo_Inf+MeCo_ECM+MeCo_Ch, family=binomial(link=logit), data = surv_data)    
summary(mod.stagecat)  
# NOTE: in this case as what we are predicting is no longer the survival of the patients, the pvalues have changed
# in the overall model the most significant pvalue is associated with MeCo_Inf

