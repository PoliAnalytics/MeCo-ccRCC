# Survival Analysis
BiocManager::install('survival')
BiocManager::install('survminer')
library(survival)
library(survminer)
# set up a dataset with the variables of interest:
surv_data <- data.frame(sequenced_patients) 
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

# Perform survival analysis
# fit the Kaplan-Meier estimator and plot the curve considering no stratification
fit <- survfit(Surv(time, status==2) ~ 1, data = surv_data)
summary(fit)
print(fit)
plot(fit, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for ccRCC Cancer Survival based on sex")
grid()

# investigate survival by gender
fit.sex <- survfit(Surv(time, status) ~ gender, data = surv_data)
summary(fit.sex)
print(fit.sex)
plot(fit.sex, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on gender",col=c("orchid2","dodgerblue2"))
legend('topright', legend=c('Male','Female'), lty=c(1,1), col=c("dodgerblue2", "orchid2"))
grid()
# perform log rank test
survdiff(Surv(time, status) ~ gender, data = surv_data) 
# p=0.6 there is strong statistical evidence against a difference between the survival of the two sexes

# investigate survival based on laterality of the tumor
fit.lat <- survfit(Surv(time, status==2) ~ laterality, data = surv_data)
summary(fit.lat)
print(fit.lat)
plot(fit.lat, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col=c('red','green','blue'),
     main="Kaplan-Meier Curve for ccRCC Cancer Survival based on laterality")
legend('topright', legend=c('Bilateral','Left','Right'), lty=c(1,1), col=c('red','green','blue'))
grid()
# perform log rank test
survdiff(Surv(time, status) ~ laterality, data = surv_data) 
# p=0.06 statistical evidence is not sufficient to state a difference in the survival of patients based on laterality of the tumor

# investigate survival based on tumor stage
fit.stage <- survfit(Surv(time, status==2) ~ stage, data = surv_data)
summary(fit.stage)
print(fit.stage)
plot(fit.stage, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col=c('red','green','blue','orange'),
     main="Kaplan-Meier Curve for ccRCC Cancer Survival based on tumor stage")
legend('topright', legend=c('I','II','III','IV'), lty=c(1,1), col=c('red','green','blue','orange'))
grid()
# perform log rank test
survdiff(Surv(time, status) ~ stage, data = surv_data) 
# p<2e-16 there is strong statistical evidence for a difference in the survival of patients at different tumor stages

# investigate survival by presence or absence of metastasis
fit.meta <- survfit(Surv(time, status) ~ metastasis, data = surv_data)
summary(fit.meta)
print(fit.meta)
plot(fit.meta, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on metastasis",col=c("green","red"))
legend('topright', legend=c('NO','YES'), lty=c(1,1), col=c("green", "red"))
grid()
# perform log rank test
survdiff(Surv(time, status) ~ metastasis, data = surv_data) 
# p=6e-06 there is strong statistical evidence for a difference in the survival of patients with and without metastasis

# investigate survival by age
# in order to do so, visualize the distribution of the variable to establish a cut-off point so to categorize the data
hist(surv_data$age, xlab='Age', main='Histogram of Age in ccRCC data')
summary(surv_data$age)
# consider the median as a cut-off
surv_data$agecat <- cut(surv_data$age, breaks=c(0,61,Inf), label=c('young','old'))
fit.agecat <- survfit(Surv(time,status)~agecat, data=surv_data)
summary(fit.agecat)
print(fit.agecat)
plot(fit.agecat, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on age",col=c("blue","green"))
legend('topright', legend=c('Young','Old'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survdiff(Surv(time,status)~agecat, data=surv_data)
# p=2e-04 there is strong statistical evidence for a difference in the survival of younger and older patients 

# investigate survival by age and sex 
fit.agecat.sex <- survfit(Surv(time,status)~agecat+gender, data=surv_data)
print(fit.agecat.sex)
plot(fit.agecat.sex, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on age and sex",col=c("pink","blue","red","green"))
legend('topright', legend=c('Young Female','Young Male','Old Female','Old Male'), lty=c(1,1), col=c("pink","blue","red","green"))
grid()
# perform Log-Rank test between sexes 
survdiff(Surv(time,status)~gender , data=surv_data[surv_data$agecat=='young',])
survdiff(Surv(time,status)~gender , data=surv_data[surv_data$agecat=='old',])
# within the same age group there seems to be no difference in survival based on sex
# consider the case of localized tumor and spread tumor
survdiff(Surv(time,status)~gender, data=surv_data[surv_data$agecat=='young'& surv_data$metastasis=='NO',])
survdiff(Surv(time,status)~gender, data=surv_data[surv_data$agecat=='young'& surv_data$metastasis=='YES',])
survdiff(Surv(time,status)~gender, data=surv_data[surv_data$agecat=='old'& surv_data$metastasis=='NO',])
survdiff(Surv(time,status)~gender, data=surv_data[surv_data$agecat=='old'& surv_data$metastasis=='YES',])
# within the same age group and considering the same tumor condition there seems to be no difference in survival based on sex 

# investigate MeCo (non refined)
hist(surv_data$MeCo)
summary(surv_data$MeCo)
# cut-off=median(MeCo)=0.19318
surv_data$mecocat <- cut(surv_data$MeCo, breaks=c(-Inf,0.19318,Inf), labels=c('Low MeCo','High MeCo'))
fit.mecocat <- survfit(Surv(time,status)~mecocat, data=surv_data)
summary(fit.mecocat)
print(fit.mecocat)
plot(fit.mecocat, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo",col=c("blue","green"))
legend('topright', legend=c('Low MeCo','High MeCo'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survdiff(Surv(time,status)~mecocat, data=surv_data)
# p=0.9 there is no statistical evidence of a difference in the survival of patients with higher and lower meco scores

# investigate survival by MeCo and sex
fit.mecocat.sex <- survfit(Surv(time, status)~mecocat+gender, data=surv_data)
print(fit.mecocat.sex)
plot(fit.agecat.sex, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo and sex",col=c("pink","blue","red","green"))
legend('topright', legend=c('Low MeCo Female','Low MeCo Male','High MeCo Female','High MeCo Male'), lty=c(1,1), col=c("pink","blue","red","green"))
grid()
# perform Log-Rank test between MeCos
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$gender=='male',])
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$gender=='female',])
# within the same sex there seems to be no difference in survival based on MeCo

# investigate survival by MeCo and age
fit.mecocat.agecat <- survfit(Surv(time, status)~mecocat+agecat, data=surv_data)
print(fit.mecocat.agecat)
plot(fit.mecocat.agecat, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo and Age",col=c("pink","blue","red","green"))
legend('topright', legend=c('Low MeCo Young','Low MeCo Old','High MeCo Young','High MeCo Old'), lty=c(1,1), col=c("pink","blue","red","green"))
grid()
# perform Log-Rank test between MeCos
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$agecat=='young',])
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$agecat=='old',])
# within the same age group there seems to be no difference in survival based on MeCo

# In order to better evaluate the relationship between covariates and the role of quantitative predictors, such as MeCo scores
# now build a multivariate Cox Proportional Hazard model 

# firstly, consider age, gender, stage, metastasis
cox.mecoless<- coxph(Surv(time, status) ~ age+gender+stage+metastasis, data = surv_data)
summary(cox.mecoless)
print(cox.mecoless)
ggforest(cox.mecoless, data=surv_data)
# among the covariates only age and stage are statistically significant, yet as the CI of age, stage I and II include 1
# they only make a smaller difference when compared to stage III and stage IV

# include MeCo (not refined)
cox.meco<- coxph(Surv(time, status) ~ age+gender+stage+metastasis+MeCo, data = surv_data)
summary(cox.meco)
print(cox.meco)
ggforest(cox.meco, data=surv_data)
# MeCo score appears to be more useful than gender and metastasis in predicting survival, yet it is not particularly significant
# moreover, the confidence interval for the hazard ratio is the widest

# investigate only refined MeCos
cox.mecoonly <- coxph(Surv(time, status)~MeCo_Pro+MeCo_Ant+MeCo_Ch+MeCo_ECM+MeCo_Inf, data=surv_data)
summary(cox.mecoonly)
print(cox.mecoonly)
ggforest(cox.mecoonly, data=surv_data)
# also the refined MeCo scores bear no usefulness in predicting the survival of patients 
# again note the wide confidence intervals

# visualize MeCo scores  
par(mfrow=c(3,2))
hist(surv_data$MeCo)
hist(surv_data$MeCo_Pro)
hist(surv_data$MeCo_ECM)
hist(surv_data$MeCo_Ant)
hist(surv_data$MeCo_Ch)
hist(surv_data$MeCo_Inf)
# standardize values
standard <- function(x){
  (x-mean(x))/sd(x)
}

cox.mecoonly2 <- coxph(Surv(time, status)~standard(MeCo_Pro)+standard(MeCo_Ant)+standard(MeCo_Ch)+standard(MeCo_ECM)+standard(MeCo_Inf), data=surv_data)
summary(cox.mecoonly2)
print(cox.mecoonly2)
ggforest(cox.mecoonly2, data=surv_data)
# standardization has narrowed confidence intervals yet no score is showing improved pvalues.
