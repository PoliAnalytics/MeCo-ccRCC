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
x   <- seq(min(surv_data$age),max(surv_data$age),10)

# Central point in the intervals of ages
mid <- c((x[-1]+x[-7])/2)

# Divide the data in classes
GRAGE <- cut(surv_data$age, breaks=x)   # cut is a very convenient command to 
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
exp(logistic_MeCo$coefficients)
# MeCo has a significant effect on the probability of death: its unitary increase
# decreases the probability of death of 94%

# OR for an age increment of 0.1 of MeCo score
exp(0.1*coef(logistic_MeCo)[2])   
# the probability of death in case of 0.1 increment of MeCo score is 25% lower

# confidence intervals for the coefficients as:
cis <- confint.default(logistic_MeCo)
cis

# CI for the OR of 0.1 of MeCo score
exp(0.1*cis[2,])
# the probability of death in case of 0.1 increment decreases between the 6% and the 40%

# confidence intervals for the prediction of either logit(p|AGE) or p|AGE:
x.new=0.2
pred <- predict(logistic_MeCo, data.frame(MeCo=x.new), se=TRUE) 

# The logit(p|MeCo=0.2) is:
pred$fit 

# its standard error
pred$se.fit

# Wald type confidence interval for logit(p) (level 95%)
alpha <- .05
c('Inf'=pred$fit-qnorm(1-alpha/2)*pred$se, ## predicted value - quantile*standard_deviation
  'Center'=pred$fit,
  'Up'=pred$fit+qnorm(1-alpha/2)*pred$se)
# CI: [-0.8529719, -0.4675687]


# fit the logistic model by MeCo_reg
logistic_MeCoReg <- glm(status~ MeCo_reg, family=binomial(link=logit), data = surv_data)
summary(logistic_MeCoReg) 
exp(logistic_MeCoReg$coefficients)
# Meco regulation has a significant effect on the probability of death: its unitary increase 
# decreases the probability of death of 96%

# OR for an age increment of 0.1 of MeCo regulation score
exp(0.1*coef(logistic_MeCoReg)[2])   
# the probability of death in case of 0.1 increment of MeCo regulation score is 27% lower

# confidence intervals for the coefficients as:
cis <- confint.default(logistic_MeCoReg)
cis

# CI for the OR of 0.1 of MeCo regulation score
exp(0.1*cis[2,])
# the probability of death in case of 0.1 increment decreases between the 13% and the 40%

# confidence intervals for the prediction of either logit(p|AGE) or p|AGE:
x.new=0.2
pred <- predict(logistic_MeCohReg, data.frame(MeCo_reg=x.new), se=TRUE) 

# The logit(p|MeCo_reg=0.2) is:
pred$fit 

# its standard error
pred$se.fit

# Wald type confidence interval for logit(p) (level 95%)
alpha <- .05
c('Inf'=pred$fit-qnorm(1-alpha/2)*pred$se, ## predicted value - quantile*standard_deviation
  'Center'=pred$fit,
  'Up'=pred$fit+qnorm(1-alpha/2)*pred$se)
# CI: [-1.4098515, -0.8314904]


# fit the logistic model by age, stage, MeCo
log1 <- glm(status~ age+stage+MeCo, family=binomial(link=logit), data = surv_data)
summary(log1) 
exp(log1$coefficients)
# being  of stage 2 with respect to stage 1 has no effect on the probability of death.
# age is a significant risk factor: its unitary increase increases the probability of death of 4.9%
# being of stage 3 is a significant risk factor that increases the probability of death of 2.9 times
# with respect to being of stage 1;
# being of stage 4 is a significant risk factor that increases the probability of death of 24 times
# with respect to being of stage 1;
# MeCo is a significant protective factor: its unitary increase decreases the probability of death 
# of 94%

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
# with respect to being of stage 1;
# MeCo regulation is a significant protective factor: its unitary increase decreases the probability of death 
# of 90%

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
# for stage 3 both status, age and MeCo regulation have an effect: being dead increases the probability
# of being of stage 3, with respect to stage 1, of 2.9 times. Unitary increases of age increases the probability
# of being of stage 3, with respect to 1, of 2.5%, unitary increases of MeCo regulation decreases the probability
# of being of stage 3, with respect to 1, of 99%
# for stage 4 both status and MeCo regulation have an effect: being dead increases the probability
# of being of stage 4, with respect to stage 1, of 22 times. Unitary increases of MeCo regulation
# decreases the probability of being of stage 4, with respect to to 1, of 99%
