# Logistic Regression

# install.packages("MASS")
# install.packages("pROC")
# install.packages("MLmetrics")

library(MASS)
library(pROC)
library(MLmetrics)

# since in logistic model 0<=y<=1:
surv_data$status[which(surv_data$status==1)] <- 0 # alive
surv_data$status[which(surv_data$status==2)] <- 1 # dead

# fit the logistic model by stage
logistic_stage<- glm(status~ stage, family=binomial(link=logit), data = surv_data)
summary(logistic_stage) 
exp(logistic_stage$coefficients)
# being of stage 2 has no effect on the probability of death.
# being of stage 3 is a significant risk factor that increases the probability of death
# of 3.5 times with respect to being of stage 1.
# being of stage 4 is a significant risk factor that increases the probability of death
# of 21.8 times with respect to being of stage 1.

# fit the logistic model by metastasis
logistic_metastasis <- glm(status~ metastasis, family=binomial(link=logit), data = surv_data)
summary(logistic_metastasis) 
exp(logistic_metastasis$coefficients)
# Having metastasis is a significant risk factor that increases the probability of death of
# 3.24 times with respect to not having them.

# fit the logistic model by age
logistic_age <- glm(status~ age, family=binomial(link=logit), data = surv_data)
summary(logistic_age) 
exp(logistic_age$coefficients)
# age is a significant risk factor: its unitary increase increases the probability of death
# of 3.7%

# Define the vector x of limits for the AGE classes
x   <- seq(min(surv_data$age),max(surv_data$age),1)

# Central point in the intervals of ages
mid <- c((x[-1]+x[-7])/2)

# Divide the data in classes
GRAGE <- cut(surv_data$age, breaks=x)    
tab <- table(GRAGE)
tab

y <- tapply(surv_data$status, GRAGE, mean)
y

plot(surv_data$age, surv_data$status, pch=19, col=surv_data$status+1, ylim=c(0,1), xlab = 'Age',ylab = 'Status', main='Status by age')
abline(v=x, col='grey', lty=2)
points(mid, y, col="blue", pch=3, lwd=2)

library(tidyverse)

d <- data.frame(age =surv_data$age, fitted=logistic_age$fitted)
s <- arrange(d, age)

lines(s$age, s$fitted, col='blue')

# OR for an age increment of 10 years
exp(10*coef(logistic_age)[2])   
# the probability of death after 10 years is 44% higher

# confidence intervals for the coefficients as:
cis <- confint.default(logistic_age)
cis

# CI for the OR of age for an increment of 10 years
exp(10*cis[2,])
# the probability of death after 10 years increases between the 22.8% and the 70%

# confidence intervals for the prediction of either logit(p|AGE) or p|AGE:
x.new=60
pred <- predict(logistic_age, data.frame(age=x.new), se=TRUE) 

# The logit(p|AGE=60) is:
pred$fit 

# its standard error
pred$se.fit

# Wald type confidence interval for logit(p) (level 95%)
alpha <- .05
c('Inf'=pred$fit-qnorm(1-alpha/2)*pred$se, ## predicted value - quantile*standard_deviation
  'Center'=pred$fit,
  'Up'=pred$fit+qnorm(1-alpha/2)*pred$se)
# CI: [-0.9833312, -0.6034349]

par(mfrow=c(1,2))
# Representation in the space (AGE,logit(p))
plot(surv_data$age, coef(logistic_age)[1]+coef(logistic_age)[2]*surv_data$age, type='l',ylab='logit(p)',xlab='Age',
     main="CI in the space (AGE,logit(p))")
points(x.new,pred$fit,pch=19)
segments(x.new,pred$fit-qnorm(1-alpha/2)*pred$se,
         x.new,pred$fit+qnorm(1-alpha/2)*pred$se, col='red')
points(x.new,pred$fit-qnorm(1-alpha/2)*pred$se,pch='-',col='red')
points(x.new,pred$fit+qnorm(1-alpha/2)*pred$se,pch='-',col='red')

# Representation in the space (AGE,p)
gl <- binomial(link=logit)    # link function

d <- data.frame(age =surv_data$age, lpred=logistic_age$linear.predictors)
s <- arrange(d, age)

plot(s$age, gl$linkinv(s$lpred), type='l',ylab='p',xlab='Age',
     main="CI in the space (AGE,p)")
points(x.new,gl$linkinv(pred$fit),pch=19)
segments(x.new,gl$linkinv(pred$fit-qnorm(1-alpha/2)*pred$se),
         x.new,gl$linkinv(pred$fit+qnorm(1-alpha/2)*pred$se), col='red')
points(x.new,gl$linkinv(pred$fit-qnorm(1-alpha/2)*pred$se),pch='-',col='red')
points(x.new,gl$linkinv(pred$fit+qnorm(1-alpha/2)*pred$se),pch='-',col='red')


# fit the logistic model by MeCo (non refined)
logistic_MeCo <- glm(status~ MeCo, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCo) 

# since the MeCo score varies between -0.3449 and 0.4179 talking about a unitary increase doesn't
# make sense but we have to look at the effect of an increase of 0.01 of the score since it
# is much more plausible 

exp(0.01*logistic_MeCo$coefficients)
# MeCo has a significant effect on the probability of death: its 0.01 increase
# decreases the probability of death of 2.9%

# Define the vector x of limits for the MeCo classes
x <- seq(min(surv_data$MeCo),max(surv_data$MeCo),0.01)

# Central point in the intervals of MeCo
mid <- c((x[-1]+x[-77])/2)

# Divide the data in classes
GRAGE <- cut(surv_data$MeCo, breaks=x)    
tab <- table(GRAGE)
tab

y <- tapply(surv_data$status, GRAGE, mean)
y

plot(surv_data$MeCo, surv_data$status, pch=19, col=surv_data$status+1, ylim=c(0,1), xlab = 'MeCo',ylab = 'Status', main='Status by MeCo')
abline(v=x, col='grey', lty=2)
points(mid, y, col="blue", pch=3, lwd=2)

d <- data.frame(MeCo =surv_data$MeCo, fitted=logistic_MeCo$fitted)
s <- arrange(d, MeCo)

lines(s$MeCo, s$fitted, col='blue')


# confidence intervals for the coefficients as:
cis <- confint.default(logistic_MeCo)
cis

# CI for the OR of 0.01 of MeCo score
exp(0.01*cis[2,])
# the probability of death in case of 0.01 increment decreases between the 0.01% and the 5%

# confidence intervals for the prediction of either logit(p|MeCo) or p|MeCo:
x.new=0.05
pred <- predict(logistic_MeCo, data.frame(MeCo=x.new), se=TRUE) 

# The logit(p|MeCo=0.05) is:
pred$fit 

# its standard error
pred$se.fit

# Wald type confidence interval for logit(p) (level 95%)
alpha <- .05
c('Inf'=pred$fit-qnorm(1-alpha/2)*pred$se, ## predicted value - quantile*standard_deviation
  'Center'=pred$fit,
  'Up'=pred$fit+qnorm(1-alpha/2)*pred$se)
# CI: [-0.6588410, 0.2059039]

par(mfrow=c(1,2))
# Representation in the space (MeCo,logit(p))
plot(surv_data$MeCo, coef(logistic_MeCo)[1]+coef(logistic_MeCo)[2]*surv_data$MeCo, type='l',ylab='logit(p)',xlab='MeCo',
     main="CI in the space (MeCo,logit(p))")
points(x.new,pred$fit,pch=19)
segments(x.new,pred$fit-qnorm(1-alpha/2)*pred$se,
         x.new,pred$fit+qnorm(1-alpha/2)*pred$se, col='red')
points(x.new,pred$fit-qnorm(1-alpha/2)*pred$se,pch='-',col='red')
points(x.new,pred$fit+qnorm(1-alpha/2)*pred$se,pch='-',col='red')

# Representation in the space (MeCo,p)
gl <- binomial(link=logit)    # link function

d <- data.frame(MeCo =surv_data$MeCo, lpred=gl$linkinv(logistic_MeCo$linear.predictors))
s <- arrange(d, MeCo)

plot(s$MeCo, s$lpred, type='l',ylab='p',xlab='MeCo',
     main="CI in the space (MeCo,p)")
points(x.new,gl$linkinv(pred$fit),pch=19)
segments(x.new,gl$linkinv(pred$fit-qnorm(1-alpha/2)*pred$se),
         x.new,gl$linkinv(pred$fit+qnorm(1-alpha/2)*pred$se), col='red')
points(x.new,gl$linkinv(pred$fit-qnorm(1-alpha/2)*pred$se),pch='-',col='red')
points(x.new,gl$linkinv(pred$fit+qnorm(1-alpha/2)*pred$se),pch='-',col='red')


# fit the logistic model by MeCo_reg
logistic_MeCoReg <- glm(status~ MeCo_reg, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoReg) 

# since the MeCo regulation score varies between -0.62646 and 0.45316 talking about a unitary increase doesn't
# make sense but we have to look at the effect of an increase of 0.01 of the score since it
# is much more plausible 

exp(0.01*logistic_MeCoReg$coefficients)
# Meco regulation has a significant effect on the probability of death: its 0.01 increase 
# decreases the probability of death of 3.2%

# Define the vector x of limits for the MeCo_reg classes
x <- seq(min(surv_data$MeCo_reg),max(surv_data$MeCo_reg),0.01)

# Central point in the intervals of MeCo_reg
mid <- c((x[-1]+x[-108])/2)

# Divide the data in classes
GRAGE <- cut(surv_data$MeCo_reg, breaks=x)    
tab <- table(GRAGE)
tab

y <- tapply(surv_data$status, GRAGE, mean)
y

plot(surv_data$MeCo_reg, surv_data$status, pch=19, col=surv_data$status+1, ylim=c(0,1), xlab = 'MeCo_reg',ylab = 'Status', main='Status by MeCo regulation')
abline(v=x, col='grey', lty=2)
points(mid, y, col="blue", pch=3, lwd=2)

d <- data.frame(MeCo_reg =surv_data$MeCo_reg, fitted=logistic_MeCoReg$fitted)
s <- arrange(d, MeCo_reg)

lines(s$MeCo_reg, s$fitted, col='blue')


# confidence intervals for the coefficients as:
cis <- confint.default(logistic_MeCoReg)
cis

# CI for the OR of 0.01 of MeCo regulation score
exp(0.01*cis[2,])
# the probability of death in case of 0.01 increment decreases between the 0.1% and the 5%

# confidence intervals for the prediction of either logit(p|MeCo_reg) or p|MeCo_reg:
x.new=-0.05
pred <- predict(logistic_MeCoReg, data.frame(MeCo_reg=x.new), se=TRUE) 

# The logit(p|MeCo_reg=-0.05) is:
pred$fit 

# its standard error
pred$se.fit

# Wald type confidence interval for logit(p) (level 95%)
alpha <- .05
c('Inf'=pred$fit-qnorm(1-alpha/2)*pred$se, ## predicted value - quantile*standard_deviation
  'Center'=pred$fit,
  'Up'=pred$fit+qnorm(1-alpha/2)*pred$se)
# CI: [-0.61753718, -0.02623533]

par(mfrow=c(1,2))
# Representation in the space (MeCo_reg,logit(p))
plot(surv_data$MeCo_reg, coef(logistic_MeCoReg)[1]+coef(logistic_MeCoReg)[2]*surv_data$MeCo_reg, type='l',ylab='logit(p)',xlab='MeCo_reg',
     main="CI in the space (MeCo_reg,logit(p))")
points(x.new,pred$fit,pch=19)
segments(x.new,pred$fit-qnorm(1-alpha/2)*pred$se,
         x.new,pred$fit+qnorm(1-alpha/2)*pred$se, col='red')
points(x.new,pred$fit-qnorm(1-alpha/2)*pred$se,pch='-',col='red')
points(x.new,pred$fit+qnorm(1-alpha/2)*pred$se,pch='-',col='red')

# Representation in the space (MeCo_reg,p)
gl <- binomial(link=logit)    # link function

d <- data.frame(MeCo_reg =surv_data$MeCo_reg, lpred=gl$linkinv(logistic_MeCoReg$linear.predictors))
s <- arrange(d, MeCo_reg)

plot(s$MeCo_reg, s$lpred, type='l',ylab='p',xlab='MeCo_reg',
     main="CI in the space (MeCo_reg,p)")
points(x.new,gl$linkinv(pred$fit),pch=19)
segments(x.new,gl$linkinv(pred$fit-qnorm(1-alpha/2)*pred$se),
         x.new,gl$linkinv(pred$fit+qnorm(1-alpha/2)*pred$se), col='red')
points(x.new,gl$linkinv(pred$fit-qnorm(1-alpha/2)*pred$se),pch='-',col='red')
points(x.new,gl$linkinv(pred$fit+qnorm(1-alpha/2)*pred$se),pch='-',col='red')


# fit the logistic model by age, stage, MeCo
log1 <- glm(status~ age+stage+MeCo, family=binomial(link=logit), data = surv_data)
summary(log1) 

exp(log1$coefficients)
# being  of stage 2 with respect to stage 1 has no effect on the probability of death.
# age is a significant risk factor: its unitary increase increases the probability of death of 4.9%
# being of stage 3 is a significant risk factor that increases the probability of death of 2.9 times
# with respect to being of stage 1;
# being of stage 4 is a significant risk factor that increases the probability of death of 24 times
# with respect to being of stage 1.

exp(0.01*-2.851267)
# MeCo is a significant protective factor: its 0.01 increase decreases the probability of death 
# of 2.9%

anova(log1,logistic_MeCo, test='LRT')
# Since the p-value of LRT test is low, the full model- with age, stage and MeCo as predictors - is more
# informative

# fit the logistic model by age, stage, meco_reg
log2 <- glm(status~ age+stage+MeCo_reg, family=binomial(link=logit), data = surv_data)
summary(log2) 

exp(log2$coefficients)
# being  of stage 2 with respect to stage 1 has no effect on the probability of death.
# age is a significant risk factor: its unitary increase increases the probability of death of 4.9%
# being of stage 3 is a significant risk factor that increases the probability of death of 2.9 times
# with respect to being of stage 1;
# being of stage 4 is a significant risk factor that increases the probability of death of 22.7 times
# with respect to being of stage 1.

exp(0.01*-2.244392)
# MeCo regulation is a significant protective factor: its 0.01 increase decreases the probability of death 
# of 2.3%

anova(log2,logistic_MeCoReg,test='LRT')
# Since the p-value of LRT test is low, the full model- with age, stage and MeCo regulation 
# as predictors - is more informative

# Goodness Of Fit (GOF) of the model with age, stage and MeCo predictors
p_threshold = 0.5

Y.hat <- ifelse(log1$fitted.values<p_threshold, 0, 1) 
Y.hat

# Confusion Matrix
table(Predicted = Y.hat, Observed = surv_data$status)

N <- nrow(surv_data)

# Compute the misclassification rate
errors <- (Y.hat != surv_data$status)
MIS_Rate  <- sum(errors)/N
MIS_Rate 
#22.8% of the patients are wrongly collocated

Specificity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 92% specificity

Sensitivity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 45.9% sensitivity

# use the empirical threshold to improve the sensitivity: p=sum(observed)/N
p_threshold = sum(surv_data$status)/N
# 0.32

Y.hat <- ifelse(log1$fitted.values<p_threshold, 0, 1) 
Y.hat

table(Predicted = Y.hat, Observed = surv_data$status)

errors <- (Y.hat != surv_data$status)
MIS_Rate  <- sum(errors)/N
MIS_Rate 
# 24.3% of the patients are wrongly collocated

Specificity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 80% specificity

Sensitivity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 66.5% sensitivity

ROC_curve <- roc(response = surv_data$status, predictor = log1$fitted.values,
                 levels = c('0','1'),
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")
auc(ROC_curve)
# 0.8028

coords(ROC_curve, x="best", transpose = TRUE)
# threshold specificity sensitivity 
# 0.3818019   0.8543417   0.6352941 


# Eventually we perform the computation of sensitivity, specificity and AUC performances 
# using the empirical threshold and 10-fold cross-validation.
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

  log.k <- glm(status~ age+stage+MeCo, family=binomial(link=logit), data = train.data)
  p.hat.k <- predict( log.k, newdata = data.frame(test.data), type='response' )
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
# 0.66

mean(specificity)
# 0.80

sd(specificity)
# 0.048

mean(AUC)
# 0.73


# Goodness Of Fit (GOF) of the model with age, stage and MeCo regulation predictors
p_threshold = 0.5

Y.hat <- ifelse(log2$fitted.values<p_threshold, 0, 1) 
Y.hat

# Confusion Matrix
table(Predicted = Y.hat, Observed = surv_data$status)

N <- nrow(surv_data)

# Compute the misclassification rate
errors <- (Y.hat != surv_data$status)
MIS_Rate  <- sum(errors)/N
MIS_Rate 
#21.4% of the patients are wrongly collocated

Specificity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 92.4% specificity

Sensitivity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 49.4% sensitivity

# use the empirical threshold to improve the sensitivity: p=sum(observed)/N
p_threshold = sum(surv_data$status)/N
# 0.32

Y.hat <- ifelse(log2$fitted.values<p_threshold, 0, 1) 
Y.hat

table(Predicted = Y.hat, Observed = surv_data$status)

errors <- (Y.hat != surv_data$status)
MIS_Rate  <- sum(errors)/N
MIS_Rate 
# 24.9% of the patients are wrongly collocated

Specificity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 79.2% specificity

Sensitivity(y_true = surv_data$status, y_pred = Y.hat, positive = 1)
# 66.5% sensitivity

ROC_curve <- roc(response = surv_data$status, predictor = log2$fitted.values,
                 levels = c('0','1'),
                 smooth=FALSE, plot=TRUE, print.auc=TRUE, auc.polygon=TRUE,
                 main="ROC Curve")
auc(ROC_curve)
# 0.8026

coords(ROC_curve, x="best", transpose = TRUE)
# threshold specificity sensitivity 
# 0.3651942   0.8487395   0.6411765 


# Eventually we perform the computation of sensitivity, specificity and AUC performances 
# using the empirical threshold and 10-fold cross-validation.
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
  
  log.k <- glm(status~ age+stage+MeCo_reg, family=binomial(link=logit), data = train.data)
  p.hat.k <- predict( log.k, newdata = data.frame(test.data), type='response' )
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
# 0.68

mean(specificity)
# 0.79

sd(specificity)
# 0.061

mean(AUC)
# 0.732





# build a logistic regression model for a multiclass problem

library(nnet)

# predict the stage from status, metastasis, age and meco_pro
mult_logistic_stage <- multinom(stage ~ status+age+MeCo_reg, data = surv_data)
z <- summary(mult_logistic_stage)$coefficients/summary(mult_logistic_stage)$standard.errors
z

p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

exp(coef(mult_logistic_stage))
# the only feature that affects being of stage 2 is status: being dead increases the probability
# of being of stage 2, with respect to stage 1, of 53%
exp(0.01*coef(mult_logistic_stage))
# for stage 3 both status, age and MeCo regulation have an effect: being dead increases the probability
# of being of stage 3, with respect to stage 1, of 2.9 times. Unitary increases of age increases the probability
# of being of stage 3, with respect to 1, of 2.5%, 0.01 increases of MeCo regulation decreases the probability
# of being of stage 3, with respect to 1, of 4.3%
# for stage 4 both status and MeCo regulation have an effect: being dead increases the probability
# of being of stage 4, with respect to stage 1, of 22 times. 0.01 increases of MeCo regulation
# decreases the probability of being of stage 4, with respect to to 1, of 4.3%
