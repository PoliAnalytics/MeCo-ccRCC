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

# remove time=0
which(surv_data$time==0)
surv_data <- surv_data[-c(220,302,429,440),]

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
# p=0.7 there is strong statistical evidence against a difference between the survival of the two sexes

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
# p=0.08 statistical evidence is not sufficient to state a difference in the survival of patients based on laterality of the tumor

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
# p=7e-06 there is strong statistical evidence for a difference in the survival of patients with and without metastasis

# investigate survival by age
# in order to do so, visualize the distribution of the variable to establish a cut-off point so to categorize the data
hist(surv_data$age, xlab='Age', main='Histogram of Age in ccRCC data')
summary(surv_data$age)
# consider the median as a cut-off
surv_data$agecat <- cut(surv_data$age, breaks=c(0,60,Inf), label=c('young','old'))
fit.agecat <- survfit(Surv(time,status)~agecat, data=surv_data)
summary(fit.agecat)
print(fit.agecat)
plot(fit.agecat, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on age",col=c("blue","green"))
legend('topright', legend=c('Young','Old'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survdiff(Surv(time,status)~agecat, data=surv_data)
# p=4e-04 there is strong statistical evidence for a difference in the survival of younger and older patients 

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
surv_data$mecocat <- cut(surv_data$MeCo, breaks=c(-Inf,0.19343,Inf), labels=c('Low MeCo','High MeCo'))
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

# investigate MeCo_Ant
hist(surv_data$MeCo_Ant)
summary(surv_data$MeCo_Ant)
# cut-off=median(MeCo_Ant)=-0.04365
surv_data$mecocat_ant <- cut(surv_data$MeCo_Ant, breaks=c(-Inf,-0.04365,Inf), labels=c('Low MeCo_Ant','High MeCo_Ant'))
fit.mecocat_ant <- survival::survfit(survival::Surv(time,status)~mecocat_ant, data=surv_data)
summary(fit.mecocat_ant)
print(fit.mecocat_ant)
plot(fit.mecocat_ant, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo_Ant",col=c("blue","green"))
legend('topright', legend=c('Low MeCo','High MeCo'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survival::survdiff(survival::Surv(time,status)~mecocat_ant, data=surv_data)
# p=0.6 there is no statistical evidence of a difference in the survival of patients with higher and lower meco scores

# investigate MeCo_Ch
hist(surv_data$MeCo_Ch)
summary(surv_data$MeCo_Ch)
# cut-off=median(MeCo_Ch)=-0.4313
surv_data$mecocat_Ch <- cut(surv_data$MeCo_Ch, breaks=c(-Inf,-0.9699,Inf), labels=c('Low MeCo_Ch','High MeCo_Ch'))
fit.mecocat_Ch <- survival::survfit(survival::Surv(time,status)~mecocat_Ch, data=surv_data)
summary(fit.mecocat_Ch)
print(fit.mecocat_Ch)
plot(fit.mecocat_Ch, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo_Ch",col=c("blue","green"))
legend('topright', legend=c('Low MeCo','High MeCo'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survival::survdiff(survival::Surv(time,status)~mecocat_Ch, data=surv_data)
# p=0.6 there is no statistical evidence of a difference in the survival of patients with higher and lower meco scores

# investigate MeCo_ECM
hist(surv_data$MeCo_ECM)
summary(surv_data$MeCo_ECM)
# cut-off=median(MeCo_Ch)=-0.2611
surv_data$mecocat_ECM <- cut(surv_data$MeCo_ECM, breaks=c(-Inf,-0.2602,Inf), labels=c('Low MeCo_ECM','High MeCo_ECM'))
fit.mecocat_ECM <- survival::survfit(survival::Surv(time,status)~mecocat_ECM, data=surv_data)
summary(fit.mecocat_ECM)
print(fit.mecocat_ECM)
plot(fit.mecocat_ECM, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo_ECM",col=c("blue","green"))
legend('topright', legend=c('Low MeCo','High MeCo'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survival::survdiff(survival::Surv(time,status)~mecocat_ECM, data=surv_data)
# p=0.5 there is no statistical evidence of a difference in the survival of patients with higher and lower meco scores

# investigate MeCo_Inf
hist(surv_data$MeCo_Inf)
summary(surv_data$MeCo_Inf)
# cut-off=median(MeCo_Ch)=0.6201
surv_data$mecocat_Inf <- cut(surv_data$MeCo_Inf, breaks=c(-Inf,0.6227,Inf), labels=c('Low MeCo_Inf','High MeCo_Inf'))
fit.mecocat_Inf <- survival::survfit(survival::Surv(time,status)~mecocat_Inf, data=surv_data)
summary(fit.mecocat_Inf)
print(fit.mecocat_Inf)
plot(fit.mecocat_Inf, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo_Inf",col=c("blue","green"))
legend('topright', legend=c('Low MeCo','High MeCo'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survival::survdiff(survival::Surv(time,status)~mecocat_Inf, data=surv_data)
# p=0.7 there is no statistical evidence of a difference in the survival of patients with higher and lower meco scores

# investigate MeCo_Pro
hist(surv_data$MeCo_Pro)
summary(surv_data$MeCo_Pro)
# cut-off=median(MeCo_Ch)=0.5060
surv_data$mecocat_Pro <- cut(surv_data$MeCo_Pro, breaks=c(-Inf,0.5066,Inf), labels=c('Low MeCo_Pro','High MeCo_Pro'))
fit.mecocat_Pro <- survival::survfit(survival::Surv(time,status)~mecocat_Pro, data=surv_data)
summary(fit.mecocat_Pro)
print(fit.mecocat_Pro)
plot(fit.mecocat_Pro, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', main="Kaplan-Meier Curve for ccRCC Cancer Survival based on MeCo_Pro",col=c("blue","green"))
legend('topright', legend=c('Low MeCo','High MeCo'), lty=c(1,1), col=c("blue", "green"))
grid()
# perform Log-Rank test
survival::survdiff(survival::Surv(time,status)~mecocat_Pro, data=surv_data)
# p=0.07 there is no statistical evidence of a difference in the survival of patients with higher and lower meco scores

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

cox.mecoonly2 <- survival::coxph(survival::Surv(time, status)~standard(MeCo_Pro)+standard(MeCo_Ant)+standard(MeCo_Ch)+standard(MeCo_ECM)+standard(MeCo_Inf), data=surv_data)
summary(cox.mecoonly2)
print(cox.mecoonly2)
ggforest(cox.mecoonly2, data=surv_data)
# standardization has narrowed confidence intervals yet no score is showing improved pvalues.

# investigate MeCo_Pro as it is the most promising MeCo score for prediction 
# based on the previous analysis by KM
cox.meco_Pro <- survival::coxph(survival::Surv(time, status)~MeCo_Pro, data=surv_data)
summary(cox.meco_Pro)
print(cox.meco_Pro)
ggforest(cox.meco_Pro, data=surv_data)
# MeCo_Pro appers to be significant at a level of 0.0521

# now verify the assumptions of the cox proportional hazard model 

# firstly assess goodness of fit:
# plot the martingale residuals and verify that they have mean=0 
# residuals vs linear predictions
ggcoxdiagnostics(cox.meco_Pro, type = "martingale",linear.predictions = T)
# residuals vs observation id
ggcoxdiagnostics(cox.meco_Pro, type = "martingale",linear.predictions = F)
# the condition seems to hold

# deviance residuals
ggcoxdiagnostics(cox.meco_Pro, type = "deviance",linear.predictions = T)
# mean remains close to zero, yet no simmetry is observed

# now assess proportional hazards:
# Schoenefeld residuals
diag.ph <- cox.zph(cox.meco_Pro)
ggcoxzph(diag.ph)
# residuals are centered around zero and they do not show any particular pattern along time, thus the PH assumption seems to hold

# now plot $log(-log(KM(t)))$ vs. $t$ or $log(t)$ and look for parallelism:
plot(fit.mecocat_Pro, fun='cloglog', 
     col=c("deeppink2","dodgerblue2"), lwd=2, lty=1,
     ylab="log(-log(Survival Probability))")
grid()
legend('topleft', c("Meco Low", "Meco High"),
       lty=c(1,1), lwd=c(2,2), col=c("deeppink2","dodgerblue2"))
# curves seem to be parallel: PH assumption holds

# finally, return pvalue of the scaled Shoenfeld test:
# +   $H_0$: Hazards are proportional
# +   $H_1$: Hazards are NOT proportional
diag.ph
# p=0.95, so there is strong evidence of proportionality of hazards
