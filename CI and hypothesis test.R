# set up a dataset with the variables of interest:
surv_data <- data.frame(pheno_data2$sequenced_patients) 
#age
surv_data$age <- pheno_data2$age_at_initial_pathologic_diagnosis
#time to event
surv_data$time <- NA
time_to_death <- pheno_data2$days_to_death.demographic
time_to_last_FU <- pheno_data2$days_to_last_follow_up.diagnoses

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

attach(surv_data)



# test , for each type of meco score, if it increases with age

# general meco
df <- data.frame(age=surv_data$age,MeCo)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo[which(df$age>61)]
meco_young <- df$MeCo[which(df$age<=61)]
boxplot(meco_old,meco_young,names = c('old','young'))
plot(age_old,meco_old)
plot(age_young,meco_young)
mean_mecoold <- mean(meco_old)
mean_mecoy <- mean(meco_young)

# number of observations is high -> x1 and x2 are normal 
var.test(meco_old,meco_young)
# equal variance 

# they are normal with same unknown variance -> t test using the pooled variance

# H0 -> mu_old = mu_young
# H1 -> mu_old > mu_young

s2.pooled=((length(meco_old)-1)*sd(meco_old)^2+(length(meco_young)-1)*sd(meco_young)^2)/(length(meco_old)+length(meco_young)-2)
s2.pooled

t.test(meco_old,meco_young,mu=0,paired=FALSE,var.equal=TRUE,alternative='two.sided')
# significantly different mean 
t.test(meco_old,meco_young,mu=0,paired=FALSE,var.equal=TRUE,alternative='greater')
# old people have higher meco

plot(age_old,meco_old, xlim=c(0,100))
points(age_young,meco_young, col = "green")

# CI of the difference

t <- qt(0.025,length(meco_old)+length(meco_young)-2, lower.tail = F)
CI_MeCo_age <- c((mean_mecoold-mean_mecoy)-t*sqrt(s2.pooled) , (mean_mecoold-mean_mecoy)+t*sqrt(s2.pooled))
CI_MeCo_age


# meco ant
df <- data.frame(age= surv_data$age,MeCo_Ant)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeSc_Ant[which(df$age>61)]
meco_young <- df$MeSc_Ant[which(df$age<=61)]
boxplot(meco_old,meco_young,names = c('old','young'))
plot(age_old,meco_old)
plot(age_young,meco_young)
mean_mecoold <- mean(meco_old)
mean_mecoy <- mean(meco_young)

# number of observations is high -> x1 and x2 are normal 
var.test(meco_old,meco_young)
# different variance 

# they are normal with different unknown variance -> z test using the sample variances

# H0 -> mu_old = mu_young
# H1 -> mu_old > mu_young

library(BSDA)
z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='two.sided')
# significantly different mean 
z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='greater')
# old people have higher meco

plot(age_old,meco_old, xlim=c(0,100))
points(age_young,meco_young, col = "green")

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCoAnt_age <- c((mean_mecoold-mean_mecoy)-z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))),(mean_mecoold-mean_mecoy)+z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))))
CI_MeCoAnt_age


# meco ch
df <- data.frame(age = surv_data$age,MeCo_Ch)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_Ch[which(df$age>61)]
meco_young <- df$MeCo_Ch[which(df$age<=61)]
boxplot(meco_old,meco_young,names = c('old','young'))
plot(age_old,meco_old)
plot(age_young,meco_young)
mean_mecoold <- mean(meco_old)
mean_mecoy <- mean(meco_young)

# number of observations is high -> x1 and x2 are normal 
var.test(meco_old,meco_young)
# different variance 

# they are normal with different unknown variance -> z test using the sample variances

# H0 -> mu_old = mu_young
# H1 -> mu_old > mu_young

z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='two.sided')
# significantly different mean 
z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='greater')
# old people have higher meco

plot(age_old,meco_old, xlim=c(0,100))
points(age_young,meco_young, col = "green")

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoCh_age <- c((mean_mecoold-mean_mecoy)-z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))),(mean_mecoold-mean_mecoy)+z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))))
CI_MeCoCh_age


# meco ECM
df <- data.frame(age = surv_data$age,MeCo_ECM)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_ECM[which(df$age>61)]
meco_young <- df$MeCo_ECM[which(df$age<=61)]
boxplot(meco_old,meco_young,names = c('old','young'))
plot(age_old,meco_old)
plot(age_young,meco_young)
mean_mecoold <- mean(meco_old)
mean_mecoy <- mean(meco_young)

# number of observations is high -> x1 and x2 are normal 
var.test(meco_old,meco_young)
# different variance 

# they are normal with different unknown variance -> z test using the sample variances

# H0 -> mu_old = mu_young
# H1 -> mu_old > mu_young

z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='two.sided')
# p.value 0.13 -> not significantly different mean 
z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='greater')
# p-value 0.06 -> old people don't have higher meco

plot(age_old,meco_old, xlim=c(0,100))
points(age_young,meco_young, col = "green")

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoECM_age <- c((mean_mecoold-mean_mecoy)-z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))),(mean_mecoold-mean_mecoy)+z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))))
CI_MeCoECM_age


# meco Inf
df <- data.frame(age = surv_data$age,MeCo_Inf)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_Inf[which(df$age>61)]
meco_young <- df$MeCo_Inf[which(df$age<=61)]
boxplot(meco_old,meco_young,names = c('old','young'))
plot(age_old,meco_old)
plot(age_young,meco_young)
mean_mecoold <- mean(meco_old)
mean_mecoy <- mean(meco_young)

# number of observations is high -> x1 and x2 are normal 
var.test(meco_old,meco_young)
# different variance 

# they are normal with different unknown variance -> z test using the sample variances

# H0 -> mu_old = mu_young
# H1 -> mu_old > mu_young

z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='two.sided')
# significantly different mean 
z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='greater')
# old people have higher meco

plot(age_old,meco_old, xlim=c(0,100))
points(age_young,meco_young, col = "green")

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoInf_age <- c((mean_mecoold-mean_mecoy)-z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))),(mean_mecoold-mean_mecoy)+z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))))
CI_MeCoInf_age


# meco Pro
df <- data.frame(age = surv_data$age,MeCo_Pro)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_Pro[which(df$age>61)]
meco_young <- df$MeCo_Pro[which(df$age<=61)]
boxplot(meco_old,meco_young,names = c('old','young'))
plot(age_old,meco_old)
plot(age_young,meco_young)
mean_mecoold <- mean(meco_old)
mean_mecoy <- mean(meco_young)

# number of observations is high -> x1 and x2 are normal 
var.test(meco_old,meco_young)
# different variance 

# they are normal with different unknown variance -> z test using the sample variances

# H0 -> mu_old = mu_young
# H1 -> mu_old > mu_young

z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='two.sided')
# significantly different mean 
z.test(meco_old,meco_young,mu=0,sigma.x = sd(meco_old),sigma.y = sd(meco_young),alternative='greater')
# old people have higher meco

plot(age_old,meco_old, xlim=c(0,100))
points(age_young,meco_young, col = "green")

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_age <- c((mean_mecoold-mean_mecoy)-z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))),(mean_mecoold-mean_mecoy)+z*sqrt((sd(meco_old)**2/length(meco_old))+(sd(meco_young)**2/length(meco_young))))
CI_MeCoPro_age



# meco decreases as stage increases

# meco general 

df2 <- data.frame(MeCo,stage)
meco_1 <- df2$MeCo[which(df2$stage =='stage i')]
meco_2 <- df2$MeCo[which(df2$stage =='stage ii')]
meco_3 <- df2$MeCo[which(df2$stage =='stage iii')]
meco_4 <- df2$MeCo[which(df2$stage =='stage iv')]
boxplot(meco_1*100,meco_2*100,meco_3*100,meco_4*100, names = c('stage1','stage2','stage3','stage4'))

var.test(meco_1,meco_2)
# different variance -> z test

z.test(meco_1,meco_2,mu=0,sigma.x = sd(meco_1),sigma.y = sd(meco_2),alternative='two.sided')
# significantly different mean 
z.test(meco_1,meco_2,mu=0,sigma.x = sd(meco_1),sigma.y = sd(meco_2),alternative='less')
# people of stage 1 have a lower general meco score than people of stage 2

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCo_stage12 <- c((mean(meco_1)-mean(meco_2))-z*sqrt((sd(meco_1)**2/length(meco_1))+(sd(meco_2)**2/length(meco_2))),(mean(meco_1)-mean(meco_2))+z*sqrt((sd(meco_1)**2/length(meco_1))+(sd(meco_2)**2/length(meco_2))))
CI_MeCo_stage12


var.test(meco_3,meco_2)
# different variance -> z test

z.test(meco_2,meco_3,mu=0,sigma.x = sd(meco_2),sigma.y = sd(meco_3),alternative='two.sided')
# not significantly different mean 
z.test(meco_2,meco_3,mu=0,sigma.x = sd(meco_2),sigma.y = sd(meco_3),alternative='less')
# people of stage 2 don't have a lower general meco score than people of stage 3

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCo_stage23 <- c((mean(meco_2)-mean(meco_3))-z*sqrt((sd(meco_2)**2/length(meco_2))+(sd(meco_3)**2/length(meco_3))),(mean(meco_2)-mean(meco_3))+z*sqrt((sd(meco_2)**2/length(meco_2))+(sd(meco_3)**2/length(meco_3))))
CI_MeCo_stage23



var.test(meco_3,meco_4)
# different variance -> z test

z.test(meco_3,meco_4,mu=0,sigma.x = sd(meco_3),sigma.y = sd(meco_4),alternative='two.sided')
# significantly different mean 
z.test(meco_3,meco_4,mu=0,sigma.x = sd(meco_3),sigma.y = sd(meco_4),alternative='less')
# people of stage 3 have a lower general meco score than people of stage 4

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCo_stage34 <- c((mean(meco_3)-mean(meco_4))-z*sqrt((sd(meco_3)**2/length(meco_3))+(sd(meco_4)**2/length(meco_4))),(mean(meco_3)-mean(meco_4))+z*sqrt((sd(meco_3)**2/length(meco_3))+(sd(meco_4)**2/length(meco_4))))
CI_MeCo_stage34



# meco pro 

df2 <- data.frame(MeCo_Pro,stage)
meco_1 <- df2$MeCo_Pro[which(df2$stage =='stage i')]
meco_2 <- df2$MeCo_Pro[which(df2$stage =='stage ii')]
meco_3 <- df2$MeCo_Pro[which(df2$stage =='stage iii')]
meco_4 <- df2$MeCo_Pro[which(df2$stage =='stage iv')]
boxplot(meco_1*100,meco_2*100,meco_3*100,meco_4*100, names = c('stage1','stage2','stage3','stage4'))

var.test(meco_1,meco_2)
# different variance -> z test

z.test(meco_1,meco_2,mu=0,sigma.x = sd(meco_1),sigma.y = sd(meco_2),alternative='two.sided')
# significantly different mean 
z.test(meco_1,meco_2,mu=0,sigma.x = sd(meco_1),sigma.y = sd(meco_2),alternative='less')
# people of stage 1 have a lower meco_pro score than people of stage 2

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_stage12 <- c((mean(meco_1)-mean(meco_2))-z*sqrt((sd(meco_1)**2/length(meco_1))+(sd(meco_2)**2/length(meco_2))),(mean(meco_1)-mean(meco_2))+z*sqrt((sd(meco_1)**2/length(meco_1))+(sd(meco_2)**2/length(meco_2))))
CI_MeCoPro_stage12



var.test(meco_3,meco_2)
# different variance -> z test

z.test(meco_2,meco_3,mu=0,sigma.x = sd(meco_2),sigma.y = sd(meco_3),alternative='two.sided')
#  significantly different mean 
z.test(meco_2,meco_3,mu=0,sigma.x = sd(meco_2),sigma.y = sd(meco_3),alternative='greater')
# people of stage 2 have a higher meco_pro score than people of stage 3


# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_stage23 <- c((mean(meco_2)-mean(meco_3))-z*sqrt((sd(meco_2)**2/length(meco_2))+(sd(meco_3)**2/length(meco_3))),(mean(meco_2)-mean(meco_3))+z*sqrt((sd(meco_2)**2/length(meco_2))+(sd(meco_3)**2/length(meco_3))))
CI_MeCoPro_stage23


var.test(meco_3,meco_4)
# different variance -> z test

z.test(meco_3,meco_4,mu=0,sigma.x = sd(meco_3),sigma.y = sd(meco_4),alternative='two.sided')
# not significantly different mean 
z.test(meco_3,meco_4,mu=0,sigma.x = sd(meco_3),sigma.y = sd(meco_4),alternative='less')
# people of stage 3 don't have a lower meco_pro score than people of stage 4

# CI difference 
z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_stage34 <- c((mean(meco_3)-mean(meco_4))-z*sqrt((sd(meco_3)**2/length(meco_3))+(sd(meco_4)**2/length(meco_4))),(mean(meco_3)-mean(meco_4))+z*sqrt((sd(meco_3)**2/length(meco_3))+(sd(meco_4)**2/length(meco_4))))
CI_MeCoPro_stage34


# meco score in alive people is lower than in dead

# general meco

df3 <- data.frame(MeCo,status)
meco_alive <- df3$MeCo[which(df3$status == 1)]
meco_dead <- df3$MeCo[which(df3$status == 2)]
boxplot(meco_alive, meco_dead, names = c('alive','dead'))

var.test(meco_alive,meco_dead)
# different variance

z.test(meco_dead, meco_alive,mu=0,sigma.x = sd(meco_dead),sigma.y = sd(meco_alive),alternative='two.sided')
# significantly different mean 
z.test(meco_dead, meco_alive,mu=0,sigma.x = sd(meco_dead),sigma.y = sd(meco_alive),alternative='greater')
# dead people have a greater general meco score than alive people

# CI difference
 
z <- qnorm(0.025, lower.tail = F)
CI_MeCo_status <- c((mean(meco_dead)-mean(meco_alive))-z*sqrt((sd(meco_dead)**2/length(meco_dead))+(sd(meco_alive)**2/length(meco_alive))),(mean(meco_dead)-mean(meco_alive))+z*sqrt((sd(meco_dead)**2/length(meco_dead))+(sd(meco_alive)**2/length(meco_alive))))
CI_MeCo_status



# meco pro

df3 <- data.frame(MeCo_Pro,status)
meco_alive <- df3$MeCo_Pro[which(df3$status == 1)]
meco_dead <- df3$MeCo_Pro[which(df3$status == 2)]
boxplot(meco_alive, meco_dead, names = c('alive','dead'))

var.test(meco_alive,meco_dead)
# different variance

z.test(meco_dead, meco_alive,mu=0,sigma.x = sd(meco_dead),sigma.y = sd(meco_alive),alternative='two.sided')
# significantly different mean 
z.test(meco_dead, meco_alive,mu=0,sigma.x = sd(meco_dead),sigma.y = sd(meco_alive),alternative='greater')
# dead people have a greater meco_pro score than alive people

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_status <- c((mean(meco_dead)-mean(meco_alive))-z*sqrt((sd(meco_dead)**2/length(meco_dead))+(sd(meco_alive)**2/length(meco_alive))),(mean(meco_dead)-mean(meco_alive))+z*sqrt((sd(meco_dead)**2/length(meco_dead))+(sd(meco_alive)**2/length(meco_alive))))
CI_MeCoPro_status


# difference of meco in the two sexs

# general meco

df4 <- data.frame(MeCo,gender)
meco_female <- df4$MeCo[which(df4$gender =='female')]
meco_male <- df4$MeCo[which(df4$gender =='male')]
boxplot(meco_female, meco_male, names = c('female','male'))

var.test(meco_female,meco_male)
# different variance

z.test(meco_female, meco_male,mu=0,sigma.x = sd(meco_female),sigma.y = sd(meco_male),alternative='two.sided')
# significantly different mean 
z.test(meco_female, meco_male,mu=0,sigma.x = sd(meco_female),sigma.y = sd(meco_male),alternative='less')
# female people have a lower general meco score than male

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCo_gender <- c((mean(meco_male)-mean(meco_female))-z*sqrt((sd(meco_male)**2/length(meco_male))+(sd(meco_female)**2/length(meco_female))),(mean(meco_male)-mean(meco_female))+z*sqrt((sd(meco_male)**2/length(meco_male))+(sd(meco_female)**2/length(meco_female))))
CI_MeCo_gender

# meco Pro

df4 <- data.frame(MeCo_Pro,gender)
meco_female <- df4$MeCo_Pro[which(df4$gender =='female')]
meco_male <- df4$MeCo_Pro[which(df4$gender =='male')]
boxplot(meco_female, meco_male, names = c('female','male'))

var.test(meco_female,meco_male)
# different variance

z.test(meco_female, meco_male,mu=0,sigma.x = sd(meco_female),sigma.y = sd(meco_male),alternative='two.sided')
# significantly different mean 
z.test(meco_female, meco_male,mu=0,sigma.x = sd(meco_female),sigma.y = sd(meco_male),alternative='less')
# female people have a lower general meco score than male

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_gender <- c((mean(meco_male)-mean(meco_female))-z*sqrt((sd(meco_male)**2/length(meco_male))+(sd(meco_female)**2/length(meco_female))),(mean(meco_male)-mean(meco_female))+z*sqrt((sd(meco_male)**2/length(meco_male))+(sd(meco_female)**2/length(meco_female))))
CI_MeCoPro_gender


# meco and laterality of the tumor

# general meco
df5 <- data.frame(MeCo,laterality)
meco_bilateral <- df5$MeCo[which(df5$laterality =='Bilateral')]
meco_left <- df5$MeCo[which(df5$laterality =='Left')]
meco_right <- df5$MeCo[which(df5$laterality =='Right')]
boxplot(meco_female, meco_male, names = c('female','male'))

var.test(meco_bilateral,meco_left)
# different variance

z.test(meco_left, meco_bilateral,mu=0,sigma.x = sd(meco_left),sigma.y = sd(meco_bilateral),alternative='two.sided')
# significantly different mean 
z.test(meco_left, meco_bilateral,mu=0,sigma.x = sd(meco_left),sigma.y = sd(meco_bilateral),alternative='less')
# people with left tumor have a lower general meco score than people with bilateral

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCo_bil_left <- c((mean(meco_left)-mean(meco_bilateral))-z*sqrt((sd(meco_left)**2/length(meco_left))+(sd(meco_bilateral)**2/length(meco_bilateral))),(mean(meco_left)-mean(meco_bilateral))+z*sqrt((sd(meco_left)**2/length(meco_left))+(sd(meco_bilateral)**2/length(meco_bilateral))))
CI_MeCo_bil_left 


var.test(meco_bilateral,meco_right)
# different variance

z.test(meco_right, meco_bilateral,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_bilateral),alternative='two.sided')
# significantly different mean 
z.test(meco_right, meco_bilateral,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_bilateral),alternative='less')
# people with right tumor have a lower general meco score than people with bilateral

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCo_bil_right <- c((mean(meco_right)-mean(meco_bilateral))-z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_bilateral)**2/length(meco_bilateral))),(mean(meco_right)-mean(meco_bilateral))+z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_bilateral)**2/length(meco_bilateral))))
CI_MeCo_bil_right



var.test(meco_left,meco_right)
# different variance

z.test(meco_right, meco_left,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_left),alternative='two.sided')
# significantly different mean 
z.test(meco_right, meco_left,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_left),alternative='greater')
# people with right tumor have a greater general meco score than people with left

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCo_left_right <- c((mean(meco_right)-mean(meco_left))-z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_left)**2/length(meco_left))),(mean(meco_right)-mean(meco_left))+z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_left)**2/length(meco_left))))
CI_MeCo_left_right


# meco Pro

df5 <- data.frame(MeCo_Pro,laterality)
meco_bilateral <- df5$MeCo_Pro[which(df5$laterality =='Bilateral')]
meco_left <- df5$MeCo_Pro[which(df5$laterality =='Left')]
meco_right <- df5$MeCo_Pro[which(df5$laterality =='Right')]
boxplot(meco_female, meco_male, names = c('female','male'))

var.test(meco_bilateral,meco_left)
# different variance

z.test(meco_left, meco_bilateral,mu=0,sigma.x = sd(meco_left),sigma.y = sd(meco_bilateral),alternative='two.sided')
# significantly different mean 
z.test(meco_left, meco_bilateral,mu=0,sigma.x = sd(meco_left),sigma.y = sd(meco_bilateral),alternative='less')
# people with left tumor have a lower general meco score than people with bilateral

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_bil_left <- c((mean(meco_left)-mean(meco_bilateral))-z*sqrt((sd(meco_left)**2/length(meco_left))+(sd(meco_bilateral)**2/length(meco_bilateral))),(mean(meco_left)-mean(meco_bilateral))+z*sqrt((sd(meco_left)**2/length(meco_left))+(sd(meco_bilateral)**2/length(meco_bilateral))))
CI_MeCoPro_bil_left 


var.test(meco_bilateral,meco_right)
# different variance

z.test(meco_right, meco_bilateral,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_bilateral),alternative='two.sided')
# significantly different mean 
z.test(meco_right, meco_bilateral,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_bilateral),alternative='less')
# people with right tumor have a lower general meco score than people with bilateral

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_bil_right <- c((mean(meco_right)-mean(meco_bilateral))-z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_bilateral)**2/length(meco_bilateral))),(mean(meco_right)-mean(meco_bilateral))+z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_bilateral)**2/length(meco_bilateral))))
CI_MeCoPro_bil_right



var.test(meco_left,meco_right)
# different variance

z.test(meco_right, meco_left,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_left),alternative='two.sided')
# significantly different mean 
z.test(meco_right, meco_left,mu=0,sigma.x = sd(meco_right),sigma.y = sd(meco_left),alternative='greater')
# people with right tumor have a greater general meco score than people with left

# CI difference

z <- qnorm(0.025, lower.tail = F)
CI_MeCoPro_left_right <- c((mean(meco_right)-mean(meco_left))-z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_left)**2/length(meco_left))),(mean(meco_right)-mean(meco_left))+z*sqrt((sd(meco_right)**2/length(meco_right))+(sd(meco_left)**2/length(meco_left))))
CI_MeCoPro_left_right


# meco is higher in those with metastasis?

# general meco
df6 <- data.frame(MeCo,metastasis)
meco_Yes <- df6$MeCo[which(df6$metastasis =='YES')]
meco_No <- df6$MeCo[which(df6$metastasis =='NO')]
boxplot(meco_Yes, meco_No, names = c('metastasis','no metastasis'))

var.test(meco_Yes,meco_No)
# equal variance -> t test

s2.pooled=((length(meco_Yes)-1)*sd(meco_Yes)^2+(length(meco_No)-1)*sd(meco_No)^2)/(length(meco_Yes)+length(meco_No)-2)
s2.pooled

t.test(meco_Yes,meco_No,mu=0,paired=FALSE,var.equal=TRUE,alternative='two.sided')
# not significantly different mean 
t.test(meco_Yes,meco_No,mu=0,paired=FALSE,var.equal=TRUE,alternative='greater')
# people with metastasis don't have higher meco

# CI of the difference

t <- qt(0.025,length(meco_Yes)+length(meco_No)-2, lower.tail = F)
CI_MeCo_meta <- c((mean(meco_Yes)-mean(meco_No))-t*sqrt(s2.pooled) , (mean(meco_Yes)-mean(meco_No))+t*sqrt(s2.pooled) )
CI_MeCo_meta


# meco Pro 

df6 <- data.frame(MeCo_Pro,metastasis)
meco_Yes <- df6$MeCo_Pro[which(df6$metastasis =='YES')]
meco_No <- df6$MeCo_Pro[which(df6$metastasis =='NO')]
boxplot(meco_Yes, meco_No, names = c('metastasis','no metastasis'))

var.test(meco_Yes,meco_No)
# equal variance -> t test

s2.pooled=((length(meco_Yes)-1)*sd(meco_Yes)^2+(length(meco_No)-1)*sd(meco_No)^2)/(length(meco_Yes)+length(meco_No)-2)
s2.pooled

t.test(meco_Yes,meco_No,mu=0,paired=FALSE,var.equal=TRUE,alternative='two.sided')
# not significantly different mean 
t.test(meco_Yes,meco_No,mu=0,paired=FALSE,var.equal=TRUE,alternative='greater')
# people with metastasis don't have higher meco

# CI of the difference

t <- qt(0.025,length(meco_Yes)+length(meco_No)-2, lower.tail = F)
CI_MeCoPro_meta <- c((mean(meco_Yes)-mean(meco_No))-t*sqrt(s2.pooled) , (mean(meco_Yes)-mean(meco_No))+t*sqrt(s2.pooled) )
CI_MeCoPro_meta

