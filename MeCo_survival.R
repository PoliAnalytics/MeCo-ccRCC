# Survival Analysis

BiocManager::install('survival')
BiocManager::install('survminer')
library(survival)
library(survminer)
library(ggplot2)

# fit the Kaplan-Meier estimator and plot the curve considering no stratification
fit <- survfit(Surv(time, status==2) ~ 1, data = surv_data)
summary(fit)
print(fit)
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           conf.int = T,
           title="Kaplan-Meier Curve")

median_surv<-fit$time[fit$surv<=0.5][1]
median_surv
# 3554

# cumulative incidence
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           fun='event', # function event -> we tell ggplot to plot the cumulative incidence and not the prob of survival
           conf.int = T,
           title="Cumulative Incidence Curve")

# cumulative hazard
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           fun='cumhaz',
           conf.int = T,
           title="Cumulative Hazard Curve")


# investigate survival based on tumor stage
fit.stage <- survfit(Surv(time, status==2) ~ stage, data = surv_data)
summary(fit.stage)
print(fit.stage)
ggsurvplot(fit.stage,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           conf.int = T,
           title="Kaplan-Meier Curve by stage")
# perform log rank test
survdiff(Surv(time, status) ~ stage, data = surv_data) 
# p<2e-16 -> there is strong statistical evidence for a difference in the survival of 
# patients at different tumor stages

# investigate survival by presence or absence of metastasis
fit.meta <- survfit(Surv(time, status) ~ metastasis, data = surv_data)
summary(fit.meta)
print(fit.meta)
# perform log rank test
survdiff(Surv(time, status) ~ metastasis, data = surv_data) 
# p=7e-06 -> there is strong statistical evidence for a difference in the survival of 
# patients with and without metastasis

ggsurvplot(fit.meta,risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=250,
           conf.int = T,
           title="Kaplan-Meier Curve for Metastasis")


# investigate survival by age
# in order to do so, visualize the distribution of the variable to establish a cut-off point so to categorize the data
hist(surv_data$age, xlab='Age', main='Histogram of Age in ccRCC data')
summary(surv_data$age)
# consider the median as a cut-off
surv_data$agecat <- cut(surv_data$age, breaks=c(0,60,Inf), label=c('young','old'))
fit.agecat <- survfit(Surv(time,status)~agecat, data=surv_data)
summary(fit.agecat)
print(fit.agecat)

# perform Log-Rank test
survdiff(Surv(time,status)~agecat, data=surv_data)
# p=4e-04 -> there is strong statistical evidence for a difference in the survival of younger 
# and older patients 

ggsurvplot(fit.agecat,risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=250,
           conf.int = T,
           title="Kaplan-Meier Curve for Age Groups")

# perform Log-Rank test by gender
survdiff(Surv(time,status)~gender, data=surv_data)
# 0.7 -> no significant difference in survival in female patients vs male ones

# perform Log-Rank test by tumor laterality
survdiff(Surv(time,status)~laterality, data=surv_data)
# 0.08 -> no significant difference in survival depending on the laterality of the tumor

# investigate survival by age and sex 
fit.agecat.sex <- survfit(Surv(time,status)~agecat+gender, data=surv_data)
print(fit.agecat.sex)

# perform Log-Rank test between sexes 
survdiff(Surv(time,status)~gender , data=surv_data[surv_data$agecat=='young',])
survdiff(Surv(time,status)~gender , data=surv_data[surv_data$agecat=='old',])
# within the same age group there seems to be no difference in survival based on sex

survdiff(Surv(time,status)~agecat , data=surv_data[surv_data$gender=='male',])
survdiff(Surv(time,status)~agecat , data=surv_data[surv_data$gender=='female',])
# within patients of the same sex there is difference depending on age

# investigate MeCo (non refined)
hist(surv_data$MeCo)
summary(surv_data$MeCo)
# cut-off=median(MeCo)=0.2390
surv_data$mecocat <- cut(surv_data$MeCo, breaks=c(-0.3449, 0.1878, 0.2390, 0.2825, 0.4179))
fit.meco <- survfit(Surv(time,status)~mecocat, data=surv_data)
summary(fit.mecocat)

# perform Log-Rank test
survdiff(Surv(time,status)~mecocat, data=surv_data)
# p=0.007 -> there is statistical evidence of a difference in the survival of patients with 
# different range of MeCo score

ggsurvplot(fit.meco,risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=250,
           conf.int = T,
           title="Kaplan-Meier Curve for MeCo")

# investigate survival by MeCo and age
fit.mecocat.sex <- survfit(Surv(time, status)~mecocat+age, data=surv_data)
print(fit.mecocat.sex)

# perform Log-Rank test between MeCos
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$agecat=='young',])
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$agecat=='old',])
# within patients of the same age there seems to be no difference in survival based on MeCo

# investigate survival by MeCo and stage
fit.mecocat.stage <- survfit(Surv(time, status)~mecocat+stage, data=surv_data)
print(fit.mecocat.agecat)

# perform Log-Rank test between MeCos
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$stage== 'stage i',])
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$stage== 'stage ii',])
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$stage== 'stage iii',])
survdiff(Surv(time, status)~mecocat, data=surv_data[surv_data$stage== 'stage iv',])
# within patients of the same stage there seems to be no difference in survival based on MeCo

# investigate MeCo regulation
hist(surv_data$MeCo_reg)
summary(surv_data$MeCo_reg)
# cut-off=median(MeCo_reg)=0.08867
surv_data$mecocat_reg <- cut(surv_data$MeCo_reg, breaks=c(-0.62646, 0.02177, 0.08867, 0.15343, 0.45316))
fit.mecocat_reg <- survival::survfit(survival::Surv(time,status)~mecocat_reg, data=surv_data)
summary(fit.mecocat_reg)
print(fit.mecocat_reg)

# perform Log-Rank test
survival::survdiff(survival::Surv(time,status)~mecocat_reg, data=surv_data)
# p=0.02 -> there is statistical evidence of a difference in the survival of patients with higher 
# and lower MeCo regulation scores

ggsurvplot(fit.mecocat_reg,risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=250,
           conf.int = T,
           title="Kaplan-Meier Curve for MeCo regulation")


# In order to better evaluate the relationship between covariates and the role of the different predictors, such as MeCo scores,
# now build a multivariate Cox Proportional Hazard model 

# firstly, consider age, gender, stage
cox.mecoless<- coxph(Surv(time, status) ~ age+stage+gender, data = surv_data)
summary(cox.mecoless)
print(cox.mecoless)
ggforest(cox.mecoless, data=surv_data)
# among the covariates only age and stage have a statistically significant effect on the probability of survival.
# unitary increase of age increases the risk of death of 3.36%, 
# being of stage 3 increases the risk of death of 2.63 times with respect to being if stage 1
# being of stage 4 increases the risk of death of 7.17 times with respect to being if stage 1


# include MeCo (not refined)
cox.meco <- coxph(Surv(time, status) ~ age+stage+MeCo, data = surv_data)
summary(cox.meco)
print(cox.meco)
ggforest(cox.meco, data=surv_data)
# meco, age, stage 3 and stage 4 have a statistically significant effect on the probability of survival.
# unitary increase of age increases the risk of death of 3.6%
# unitary increase of the MeCo score decreases the risk of death of 91% 
# being of stage 3 increases the risk of death of 2.5 times with respect to being if stage 1
# being of stage 4 increases the risk of death of 6.73 times with respect to being if stage 1

# since the MeCo score varies between -0.3449 and 0.4179 talking about a unitary increase doesn't
# make sense but we have to look at the effect of an increase of 0.01 of the score since it
# is much more plausible 

exp(-2.371970*0.01)
# an increase of the MeCo score of 0.01 decreases the risk of death of 2.4%

# cox model with only MeCo score as covariate
cox.meco_only <- coxph(Surv(time, status) ~ MeCo, data = surv_data)
summary(cox.meco_only)
print(cox.meco_only)
ggforest(cox.meco_only, data=surv_data)
exp(0.01*cox.meco_only$coefficients)
# MeCo is a significant protective factor: its increase of 0.01 decreases the risk of death of 3%

# cox model with only MeCo regulation score as covariate
cox.mecoreg_only <- coxph(Surv(time, status) ~ MeCo_reg, data = surv_data)
summary(cox.mecoreg_only)
print(cox.mecoreg_only)
ggforest(cox.mecoreg_only, data=surv_data)
exp(0.01*cox.mecoreg_only$coefficients)
# MeCo regulation is a significant protective factor: its increase of 0.01 decreases the risk of 
# death of 2.8%

# include Meco regulation with only significant covariates according to KM
cox.meco2<- coxph(Surv(time, status) ~ age+stage+MeCo_reg, data = surv_data)
summary(cox.meco2)
print(cox.meco2)
ggforest(cox.meco2, data=surv_data)
# age, stage 3 and 4 and MeCo_reg have a statistically significant effect on the probability of survival.
# unitary increase of age increases the risk of death of 3.6%
# unitary increase of the MeCo_reg score decreases the risk of death of 86% 
# being of stage 3 increases the risk of death of 2.42 times with respect to being if stage 1
# being of stage 4 increases the risk of death of 6.44 times with respect to being if stage 1

# since the MeCo regulation score varies between -0.62646 and 0.45316 talking about a unitary increase doesn't
# make sense but we have to look at the effect of an increase of 0.01 of the score since it
# is much more plausible 

exp(-2.011810*0.01)
# an increase of the MeCo regulation score of 0.01 decreases the risk of death of 2%

# Predicted adjusted survival probability from a Cox model - stage 3 and MeCo general
df_stage3 <- with(surv_data,
                  data.frame(age = c(60,60,60),
                             stage = as.factor(rep('stage iii',3)),
                             MeCo = c(0.1,0.24,0.35))
)

fit.df_stage3 <- survfit(cox.meco, newdata = df_stage3)
fit.df_stage3

plot(fit.df_stage3, conf.int=T,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('bottomleft', c("MeCo = 0.1", "MeCo = 0.24", "MeCo = 0.35"),
       lty=c(1,1,1), lwd=c(2,2,2), cex=1,col=c("dodgerblue2","navy","darkmagenta"))

# Predicted adjusted survival probability from a Cox model - stage 4 ang MeCo general
df_stage4 <- with(surv_data,
                  data.frame(age = c(60,60,60),
                             stage = as.factor(rep('stage iv',3)),
                             MeCo = c(0.1,0.24,0.35))
)

fit.df_stage4 <- survfit(cox.meco, newdata = df_stage4)
fit.df_stage4

plot(fit.df_stage4, conf.int=T,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('topright', c("MeCo = 0.1", "MeCo = 0.24", "MeCo = 0.35"),
       lty=c(1,1,1), lwd=c(2,2,2), cex = 1, col=c("dodgerblue2","navy","darkmagenta"))


# Predicted adjusted survival probability from a Cox model - stage 3 and MeCo regulation
df_stage3 <- with(surv_data,
                  data.frame(age = c(60,60,60),
                             stage = as.factor(rep('stage iii',3)),
                             MeCo_reg = c(0.05,0.08,0.13))
)

fit.df_stage3 <- survfit(cox.meco2, newdata = df_stage3)
fit.df_stage3

plot(fit.df_stage3, conf.int=T,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('bottomleft', c("MeCo_reg = 0.05", "MeCo_reg = 0.08", "MeCo_reg = 0.13"),
       lty=c(1,1,1), lwd=c(2,2,2), cex=1,col=c("dodgerblue2","navy","darkmagenta"))

# Predicted adjusted survival probability from a Cox model - stage 4 and MeCo regulation
df_stage4 <- with(surv_data,
                  data.frame(age = c(60,60,60),
                             stage = as.factor(rep('stage iv',3)),
                             MeCo_reg = c(0.05,0.08,0.13))
)

fit.df_stage4 <- survfit(cox.meco2, newdata = df_stage4)
fit.df_stage4

plot(fit.df_stage4, conf.int=T,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('topright', c("MeCo_reg = 0.05", "MeCo_reg = 0.08", "MeCo_reg = 0.13"),
       lty=c(1,1,1), lwd=c(2,2,2), cex = 1, col=c("dodgerblue2","navy","darkmagenta"))



# now verify the assumptions of the cox proportional hazard model 

# firstly assess goodness of fit:
# plot the martingale residuals and verify that they have mean=0 
# residuals vs linear predictions
ggcoxdiagnostics(cox.meco, type = "martingale",linear.predictions = T, main='Martingale residuals')
# the condition seems to hold

ggcoxdiagnostics(cox.meco2, type = "martingale",linear.predictions = T, main='Martingale residuals')
# the condition seems to hold

# plot the deviance residuals and verify that they are roughly symmetrically distributed 
# about zero with a standard deviation of 1.
ggcoxdiagnostics(cox.meco, type = "deviance",linear.predictions = T, main='Deviance residuals')
# Even if the residuals are distributed around 0 there are many positive values - individuals who died 
# too soon on respect to the survival model expectancy - and many negative values - individual that “lived too long”.
# Also there are some large or small values: they are outliers, which are poorly predicted by the model.

ggcoxdiagnostics(cox.meco2, type = "deviance",linear.predictions = T, main='Deviance residuals')
# the result is very similar to the one of cox.best1

# now assess proportional hazards:
# Schoenefeld residuals are expected to be centered about zero
diag.ph <- cox.zph(cox.meco)
ggcoxzph(diag.ph)
# residuals are centered around zero and they do not show any particular pattern along time, 
# thus the PH assumption seems to hold

# Schoenefeld residuals 
diag.ph2 <- cox.zph(cox.meco2)
ggcoxzph(diag.ph2)
# also for this second model residuals are centered around zero and they do not show any particular 
# pattern along time, thus the PH assumption seems to hold

# Test for PH using scaled Schoenfeld test for PH
# H0: Hazards are proportional
# H1: Hazards are NOT proportional
# cox.zph() return tests for each covariate and for the global model.
diag.ph
# in the cox model with predictors age, stage and MeCo, the PH assumption is met by all the 
# covariates but stage 

diag.ph2
# in the cox model with predictors age, stage and MeCo_reg, the PH assumption is met by all the 
# covariates but stage 


# stratified cox model by stage using MeCo: build a different baseline hazard function for each single 
# level of stage 
mod.cox.strata <- coxph(Surv(time, status) ~ age + MeCo + strata(stage) , data =  surv_data)
summary(mod.cox.strata)
exp(0.01*-2.31770)
# both age and MeCo result statistically significant
# unitary increase of age increases the risk of death of 3.6%
# 0.01 increase of MeCo reduces the risk of death of 2.3%

test.ph.strata <- cox.zph(mod.cox.strata)
test.ph.strata
# now PH is met for all the covariates of the model

fit<-survfit(mod.cox.strata, data=surv_data)

plot(fit, conf.int = F,
     col=1:6, lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probabilities')
grid()
legend('topright', legend=unique(summary(fit)$strata), col=1:3, lwd=2, lty=1, title='Strata k')

# Predicted adjusted survival probability from a stratified Cox model - MeCo 
df_stage <- with(surv_data,
                  data.frame(age = c(60,60,60),
                             MeCo = c(0.05,0.08,0.13))
)

fit.df_stage <- survfit(mod.cox.strata, newdata = df_stage)
fit.df_stage

plot(fit.df_stage, 
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('topright', c("MeCo = 0.05", "MeCo = 0.08", "MeCo = 0.13"),
       lty=c(1,1,1), lwd=c(2,2,2), cex=1,col=c("dodgerblue2","navy","darkmagenta"))

# stratified cox model by stage using MeCo regulation: build a different baseline hazard function for 
# each single level of stage 
mod.cox.strata2 <- coxph(Surv(time, status) ~ age + MeCo_reg + strata(stage) , data =  surv_data)
summary(mod.cox.strata2)
exp(0.01*-1.976300)
# both age and MeCo regulation result statistically significant
# unitary increase of age increases the risk of death of 3.6%
# 0.01 increase of MeCo regulation reduces the risk of death of 2%

test.ph.strata2 <- cox.zph(mod.cox.strata2)
test.ph.strata2
# now PH is met for all the covariates of the model

fit2<-survfit(mod.cox.strata2, data=surv_data)

plot(fit2, conf.int = F,
     col=1:6, lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probabilities')
grid()
legend('topright', legend=unique(summary(fit2)$strata), col=1:4, lwd=2, lty=1, title='Strata k')

# Predicted adjusted survival probability from a stratified Cox model - MeCo regulation
df_stage_reg <- with(surv_data,
                 data.frame(age = c(60,60,60),
                            MeCo_reg = c(0.05,0.08,0.13))
)

fit.df_stage_reg <- survfit(mod.cox.strata2, newdata = df_stage_reg)
fit.df_stage_reg

plot(fit.df_stage_reg, 
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('bottomleft', c("MeCo_reg = 0.05", "MeCo_reg = 0.08", "MeCo_reg = 0.13"),
       lty=c(1,1,1), lwd=c(2,2,2), cex=1,col=c("dodgerblue2","navy","darkmagenta"))


# fit a Exponential AFT model with MeCo:
fit_exp <- survreg(Surv(time,status) ~  age+MeCo+stage, data = surv_data,
                   dist = "exponential")
summary(fit_exp)
# age, MeCo, stage 3 and stage 4 result significant in their affect on the survival time

exp(coef(fit_exp))
# unitary increase of age contracts the survival time of 0.97 times
# being of stage 3 wrt stage 1 contracts thesurvival time of 0.4 times
# being of stage 4 wrt stage 1 contracts the survival time of 0.15 times

exp(0.01*2.32783)
# 0.01 increase of MeCo extends the survival time of 1.02 times

# PH interpretation
exp(-1 * coef(fit_exp))

# In order to check whether the assumed distribution of survival times fits the observed data 
# sufficiently well, we can use the residuals of the model. However, because the residuals are 
# calculated based on the observed event times, they will also be censored. Hence, to 
# investigate whether they follow a particular distribution, we will need to employ a 
# graphical procedure accounting for censoring. One way to achieve this is by using the 
# Kaplan-Meier estimator of the residuals.

fitted_values <- fit_exp$linear.predictors
residuals <- log(fit_exp$y[, 1]) - fitted_values # scale=1

# We then compute the Kaplan-Meier estimate of these residuals, and we plot it
# in the code below we superimpose in the graph the assumed exponential distribution.
resKM <- survfit(Surv(residuals, status) ~ 1, data = surv_data)
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability", main = 'Check the exponential distribution of survival times')
xx <- seq(min(residuals), max(residuals), length.out = 35)
yy <- exp(- exp(xx))
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")


# fit a Exponential AFT model with MeCo regulation:
fit_exp <- survreg(Surv(time,status) ~  age+MeCo_reg+stage, data = surv_data,
                   dist = "exponential")
summary(fit_exp)
# age, MeCo_reg, stage 3 and stage 4 result significant in their affect on the survival time

exp(coef(fit_exp))
# unitary increase of age contracts the survival time of 0.96 times
# being of stage 3 wrt stage 1 contracts the survival time of 0.4 times
# being of stage 4 wrt stage 1 contracts the survival time of 0.15 times

exp(0.01*2.02869)
# 0.01 increase of MeCo regulation extends the survival time 1.02 times

# PH interpretation
exp(-1 * coef(fit_exp))

fitted_values <- fit_exp$linear.predictors
residuals <- log(fit_exp$y[, 1]) - fitted_values # scale=1

# We then compute the Kaplan-Meier estimate of these residuals, and we plot it
# in the code below we superimpose in the graph the assumed exponential distribution.
resKM <- survfit(Surv(residuals, status) ~ 1, data = surv_data)
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability", main = 'Check the exponential distribution of survival times')
xx <- seq(min(residuals), max(residuals), length.out = 35)
yy <- exp(- exp(xx))
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")


# We now fit AFT model using Weibull distribution with MeCo:
fit_weibull <- survreg(Surv(time,status) ~  age+MeCo+stage, data = surv_data)
summary(fit_weibull)
# age, MeCo, stage 3 and stage 4 result significant in their affect on the survival time

exp(coef(fit_weibull))
# unitary increase of age contracts the survival time of 0.96 times
# being of stage 3 wrt stage 1 contracts the survival time of 0.4 times
# being of stage 4 wrt stage 1 contracts the survival time of 0.14 times

exp(0.01*2.35734)
# 0.01 increase of MeCo extends the survival time of 1.02 times

# PH interpretation
shapeParameter <- 1 / fit_weibull$scale
shapeParameter
exp(-1 * shapeParameter * coef(fit_weibull))

# As before we plot the Kaplan-Meier estimator of the residuals.
fitted_values <- fit_weibull$linear.predictors
residuals <- (log(fit_weibull$y[, 1]) - fitted_values) / fit_weibull$scale

# We then compute the Kaplan-Meier estimate of these residuals, and we plot it. 
# In the code below we superimpose in the graph the assumed Weibull distribution.
resKM <- survfit(Surv(residuals, status) ~ 1, data = surv_data)
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability", main = 'Check the Weibull distribution of survival times')
xx <- seq(min(residuals), max(residuals), length.out = 35)
yy <- exp(- exp(xx))
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")

# Eventually we fit AFT model using Weibull distribution with MeCo regulation:
fit_weibull <- survreg(Surv(time,status) ~  age+MeCo_reg+stage, data = surv_data)
summary(fit_weibull)
# age, MeCo_reg, stage 3 and stage 4 result significant in their affect on the survival time

exp(coef(fit_weibull))
# unitary increase of age contracts the survival time of 0.96 times
# being of stage 3 wrt stage 1 contracts the survival time of 0.4 times
# being of stage 4 wrt stage 1 contracts the survival time of 0.15 times

exp(0.01*2.05573)
# 0.01 increase of MeCo regulation extends the survival time of 1.02 times

# PH interpretation
shapeParameter <- 1 / fit_weibull$scale
shapeParameter
exp(-1 * shapeParameter * coef(fit_weibull))

# As before we plot the Kaplan-Meier estimator of the residuals.
fitted_values <- fit_weibull$linear.predictors
residuals <- (log(fit_weibull$y[, 1]) - fitted_values) / fit_weibull$scale

# We then compute the Kaplan-Meier estimate of these residuals, and we plot it. 
# In the code below we superimpose in the graph the assumed Weibull distribution.
resKM <- survfit(Surv(residuals, status) ~ 1, data = surv_data)
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability", main = 'Check the Weibull distribution of survival times')
xx <- seq(min(residuals), max(residuals), length.out = 35)
yy <- exp(- exp(xx))
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")
