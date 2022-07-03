attach(surv_data)

# FIRST HYPOTHESIS TEST: test, for each type of MeCo score, if it increases with age

# general MeCo score
df <- data.frame(age=surv_data$age,MeCo)
summary(age)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo[which(df$age>61)]
meco_young <- df$MeCo[which(df$age<=61)]
shapiro.test(meco_old)
shapiro.test(meco_young)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_old,meco_young,names = c('old','young'), main = ('MeCo distribution in young and old people'))
plot(age_old,meco_old)
plot(age_young,meco_young)

# H0: No tendency of MeCo score depending on the age
# H1: Tendency of some type of MeCo score depending on the age
wilcox.test(meco_old, meco_young, paired=F, alternative='two.sided', conf.int = T)
# 0.0003 -> significant difference
# CI: 0.01092387 0.03632707
wilcox.test(meco_old, meco_young, paired=F, alternative='greater', conf.int = T)
# 0.0001 ->  old people have higher MeCo
# 0.01302886        Inf


# MeCo regulation
df <- data.frame(age= surv_data$age,MeCo_reg)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_reg[which(df$age>61)]
meco_young <- df$MeCo_reg[which(df$age<=61)]
shapiro.test(meco_old)
shapiro.test(meco_young)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_old,meco_young,names = c('old','young'), main = ('MeCo regulation distribution in young and old people'))
plot(age_old,meco_old)
plot(age_young,meco_young)

# H0: No tendency of MeCo regulation score depending on the age
# H1: Tendency of some type of MeCo regulation score depending on the age
wilcox.test(meco_old, meco_young, paired=F, alternative='two.sided', conf.int = T)
# 0.013 -> significant difference
# CI: 0.004560885 0.038791022
wilcox.test(meco_old, meco_young, paired=F, alternative='greater', conf.int = T)
# 0.007 ->  old people have higher MeCo regulation
# 0.007327603        Inf


# MeCo stimulus
df <- data.frame(age = surv_data$age,MeCo_st)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_st[which(df$age>61)]
meco_young <- df$MeCo_st[which(df$age<=61)]
shapiro.test(meco_old)
shapiro.test(meco_young)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_old,meco_young,names = c('old','young'),main='MeCo stimulus distribution in young and old people')
plot(age_old,meco_old)
plot(age_young,meco_young)

# H0: No tendency of MeCo stimulus score depending on the age
# H1: Tendency of some type of MeCo stimulus score depending on the age
wilcox.test(meco_old, meco_young, paired=F, alternative='two.sided', conf.int = T)
# 0.0001 -> significant difference
# CI:  0.01649375 0.04945930
wilcox.test(meco_old, meco_young, paired=F, alternative='greater', conf.int = T)
# 0.00007 ->  old people have higher MeCo stimulus
# 0.01925956         Inf


# MeCo development
df <- data.frame(age = surv_data$age,MeCo_dev)
age_old <- df$age[which(df$age>61)]
age_young <- df$age[which(df$age<=61)]
meco_old <- df$MeCo_dev[which(df$age>61)]
meco_young <- df$MeCo_dev[which(df$age<=61)]
shapiro.test(meco_old)
shapiro.test(meco_young)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_old,meco_young,names = c('old','young'), main='MeCo development distribution in young and old people')
plot(age_old,meco_old)
plot(age_young,meco_young)

# H0: No tendency of MeCo development score depending on the age
# H1: Tendency of some type of MeCo development score depending on the age
wilcox.test(meco_old, meco_young, paired=F, alternative='two.sided', conf.int = T)
# 0.0004 -> significant difference
# CI:  0.01303680 0.04421459
wilcox.test(meco_old, meco_young, paired=F, alternative='greater', conf.int = T)
# 0.0002 ->  old people have higher MeCo development
# 0.01551131         Inf



# SECOND HYPOTHESIS TEST: test whether each MeCo type decreases as stage increases

# MeCo general 
df2 <- data.frame(MeCo=surv_data$MeCo,stage=surv_data$stage)
meco_1 <- df2$MeCo[which(df2$stage =='stage i')]
meco_2 <- df2$MeCo[which(df2$stage =='stage ii')]
meco_3 <- df2$MeCo[which(df2$stage =='stage iii')]
meco_4 <- df2$MeCo[which(df2$stage =='stage iv')]
shapiro.test(meco_1)
# p-value is high so this variable is normally distributed
shapiro.test(meco_2)
shapiro.test(meco_3)
shapiro.test(meco_4)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_1*100,meco_2*100,meco_3*100,meco_4*100, names = c('stage1','stage2','stage3','stage4'), main='MeCo distribution in the different stages')

# H0: No tendency of MeCo score depending on the stage
# H1: Tendency of some type of MeCo score depending on the stage
wilcox.test(meco_1, meco_2, paired=F, alternative='two.sided', conf.int = T)
# 0.4 -> no significant difference
wilcox.test(meco_1, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.2 -> no significant difference
wilcox.test(meco_1, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.18 -> no significant difference

wilcox.test(meco_2, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.1 -> no significant difference
wilcox.test(meco_2, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.14 -> no significant difference

wilcox.test(meco_3, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.8 -> no significant difference


# MeCo regulation
df2 <- data.frame(MeCo_reg=surv_data$MeCo_reg,stage=surv_data$stage)
meco_1 <- df2$MeCo_reg[which(df2$stage =='stage i')]
meco_2 <- df2$MeCo_reg[which(df2$stage =='stage ii')]
meco_3 <- df2$MeCo_reg[which(df2$stage =='stage iii')]
meco_4 <- df2$MeCo_reg[which(df2$stage =='stage iv')]
shapiro.test(meco_1)
shapiro.test(meco_2)
# p-value are high so these variables are normally distributed
shapiro.test(meco_3)
shapiro.test(meco_4)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_1*100,meco_2*100,meco_3*100,meco_4*100, names = c('stage1','stage2','stage3','stage4'),  main='MeCo regulation distribution in the different stages')

var.test(meco_1,meco_2)
# different variance -> z test

# HO: there is no difference between the mean MeCo regulation score between people of stage 1 and 2
# HO: there is difference between the mean MeCo regulation score between people of stage 1 and 2
z.test(meco_1,meco_2,mu=0,sigma.x = sd(meco_1),sigma.y = sd(meco_2),alternative='two.sided')
# 0.9 -> not significantly different mean 

# H0: No tendency of MeCo regulation score depending on the stage
# H1: Tendency of some type of MeCo regulation score depending on the stage
wilcox.test(meco_1, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.0007 -> significant difference
# ci: 0.01511899 0.05561281
wilcox.test(meco_1, meco_3, paired=F, alternative='greater', conf.int = T)
# 0.0004 -> MeCo regulation is higher in stage 1 than stage 3
#  0.01852165        Inf
wilcox.test(meco_1, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 9.527e-05 -> significant difference
# ci: 0.02585632 0.07489848
wilcox.test(meco_1, meco_4, paired=F, alternative='greater', conf.int = T)
# 4.764e-05 -> MeCo regulation is higher in stage 1 than stage 4
# ci: 0.02958501        Inf

wilcox.test(meco_2, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.03 ->  significant difference
# ci: 0.001765502 0.066255799
wilcox.test(meco_2, meco_3, paired=F, alternative='greater', conf.int = T)
# 0.019 -> MeCo regulation is higher in stage 2 than stage 3
# ci: 0.007487142         Inf
wilcox.test(meco_2, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.009 -> significant difference
#  0.01265027 0.08624226
wilcox.test(meco_2, meco_4, paired=F, alternative='greater', conf.int = T)
# 0.004 -> MeCo regulation is higher in stage 2 than stage 4
# 0.01891998        Inf
wilcox.test(meco_3, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.27 -> no significant difference


# MeCo stimulus
df2 <- data.frame(MeCo_st=surv_data$MeCo_st,stage=surv_data$stage)
meco_1 <- df2$MeCo_st[which(df2$stage =='stage i')]
meco_2 <- df2$MeCo_st[which(df2$stage =='stage ii')]
meco_3 <- df2$MeCo_st[which(df2$stage =='stage iii')]
meco_4 <- df2$MeCo_st[which(df2$stage =='stage iv')]
shapiro.test(meco_1)
shapiro.test(meco_2)
# p-value are high so these variables are normally distributed
shapiro.test(meco_3)
shapiro.test(meco_4)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_1*100,meco_2*100,meco_3*100,meco_4*100, names = c('stage1','stage2','stage3','stage4'),  main='MeCo stimulus distribution in the different stages')

var.test(meco_1,meco_2)
# same variance -> t-test

# HO: there is no difference between the mean MeCo stimulus score between people of stage 1 and 2
# HO: there is difference between the mean MeCo stimulus score between people of stage 1 and 2
s2.pooled=((length(meco_1)-1)*sd(meco_1)^2+(length(meco_2)-1)*sd(meco_2)^2)/(length(meco_1)+length(meco_2)-2)
s2.pooled

t.test(meco_1,meco_2,mu=0,paired=FALSE,var.equal=TRUE,alternative='two.sided')
# 0.9 -> not significantly different mean 

# H0: No tendency of MeCo stimulus score depending on the stage
# H1: Tendency of some type of MeCo stimulus score depending on the stage
wilcox.test(meco_1, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.006 -> significant difference
# ci: 0.008434417 0.049883093
wilcox.test(meco_1, meco_3, paired=F, alternative='greater', conf.int = T)
# 0.003 -> MeCo stimulus score is greater in stage 1 than stage 3
#  0.01177697         Inf
wilcox.test(meco_1, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.1 -> no significant difference

wilcox.test(meco_2, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.046 ->  significant difference
# ci: 0.0003593156 0.0611257491
wilcox.test(meco_2, meco_3, paired=F, alternative='greater', conf.int = T)
# 0.023 -> MeCo stimulus score is greater in stage 2 than stage 3
# ci: 0.005343695         Inf
wilcox.test(meco_2, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.2 -> no significant difference

wilcox.test(meco_3, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.6 -> no significant difference


# MeCo development
df2 <- data.frame(MeCo_dev=surv_data$MeCo_dev,stage=surv_data$stage)
meco_1 <- df2$MeCo_dev[which(df2$stage =='stage i')]
meco_2 <- df2$MeCo_dev[which(df2$stage =='stage ii')]
meco_3 <- df2$MeCo_dev[which(df2$stage =='stage iii')]
meco_4 <- df2$MeCo_dev[which(df2$stage =='stage iv')]
shapiro.test(meco_1)
shapiro.test(meco_2)
# p-value are high so these variables are normally distributed
shapiro.test(meco_3)
shapiro.test(meco_4)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_1*100,meco_2*100,meco_3*100,meco_4*100, names = c('stage1','stage2','stage3','stage4'),  main='MeCo development distribution in the different stages')

var.test(meco_1,meco_2)
# same variance -> t-test

# HO: there is no difference between the mean MeCo development score between people of stage 1 and 2
# HO: there is difference between the mean MeCo development score between people of stage 1 and 2
s2.pooled=((length(meco_1)-1)*sd(meco_1)^2+(length(meco_2)-1)*sd(meco_2)^2)/(length(meco_1)+length(meco_2)-2)
s2.pooled

t.test(meco_1,meco_2,mu=0,paired=FALSE,var.equal=TRUE,alternative='two.sided')
# 0.87 -> not significantly different mean 

# H0: No tendency of MeCo development score depending on the stage
# H1: Tendency of some type of MeCo development score depending on the stage
wilcox.test(meco_1, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.089 -> no significant difference
wilcox.test(meco_1, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.5 -> no significant difference

wilcox.test(meco_2, meco_3, paired=F, alternative='two.sided', conf.int = T)
# 0.3 ->  no significant difference
wilcox.test(meco_2, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.6 -> no significant difference

wilcox.test(meco_3, meco_4, paired=F, alternative='two.sided', conf.int = T)
# 0.5 -> no significant difference



# THIRD HYPOTHESIS TEST: Test whether each type of MeCo score is higher in alive people than in dead ones

# general MeCo
df3 <- data.frame(MeCo,status)
meco_alive <- df3$MeCo[which(df3$status == 0)]
meco_dead <- df3$MeCo[which(df3$status == 1)]
shapiro.test(meco_alive)
shapiro.test(meco_dead)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_alive, meco_dead, names = c('alive','dead'), main='MeCo distribution in alive and dead people')

# H0: No tendency of MeCo score depending on the status
# H1: Tendency of some type of MeCo score depending on the status
wilcox.test(meco_alive, meco_dead, paired=F, alternative='two.sided', conf.int = T)
# 0.16 -> no significant difference


# MeCo regulation
df3 <- data.frame(MeCo_reg,status)
meco_alive <- df3$MeCo_reg[which(df3$status == 0)]
meco_dead <- df3$MeCo_reg[which(df3$status == 1)]
shapiro.test(meco_alive)
# this variable is normally distributed
shapiro.test(meco_dead)
# this variable is not normally distributed 

boxplot(meco_alive, meco_dead, names = c('alive','dead'), main='MeCo regulation distribution in alive and dead people')

# H0: No tendency of MeCo regulation score depending on the status
# H1: Tendency of some type of MeCo regulation score depending on the status
wilcox.test(meco_alive, meco_dead, paired=F, alternative='two.sided', conf.int = T)
# 0.0035 -> significant difference
# ci:  0.008864337 0.045091981
wilcox.test(meco_alive, meco_dead, paired=F, alternative='greater', conf.int = T)
# 0.002 -> MeCo regulation is higher in alive people than in dead ones
#ci: 0.01179347        Inf

# MeCo stimulus
df4 <- data.frame(MeCo_st,status)
meco_alive <- df4$MeCo_st[which(df4$status == 0)]
meco_dead <- df4$MeCo_st[which(df4$status == 1)]
shapiro.test(meco_alive)
# this variable is normally distributed
shapiro.test(meco_dead)
# this variable is not normally distributed 

boxplot(meco_alive, meco_dead, names = c('alive','dead'), main='MeCo stimulus distribution in alive and dead people')

# H0: No tendency of MeCo stimulus score depending on the status
# H1: Tendency of some type of MeCo stimulus score depending on the status
wilcox.test(meco_alive, meco_dead, paired=F, alternative='two.sided', conf.int = T)
# 0.12  -> no significant difference


# MeCo development
df3 <- data.frame(MeCo_dev,status)
meco_alive <- df3$MeCo_dev[which(df3$status == 0)]
meco_dead <- df3$MeCo_dev[which(df3$status == 1)]
shapiro.test(meco_alive)
# this variable is normally distributed
shapiro.test(meco_dead)
# this variable is not normally distributed

boxplot(meco_alive, meco_dead, names = c('alive','dead'), main='MeCo development distribution in alive and dead people')

# H0: No tendency of MeCo development score depending on the status
# H1: Tendency of some type of MeCo development score depending on the status
wilcox.test(meco_alive, meco_dead, paired=F, alternative='two.sided', conf.int = T)
# 0.5 -> no significant difference



# FOURTH HYPOTHESIS TEST: Test for each type of MeCo score whether it has a different level in the two sexes

# general MeCo
df4 <- data.frame(MeCo,gender)
meco_female <- df4$MeCo[which(df4$gender =='female')]
meco_male <- df4$MeCo[which(df4$gender =='male')]
shapiro.test(meco_female)
shapiro.test(meco_male)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_female, meco_male, names = c('female','male'), main='MeCo distribution by gender')

# H0: No tendency of MeCo score depending on the gender
# H1: Tendency of some type of MeCo score depending on the gender
wilcox.test(meco_female, meco_male, paired=F, alternative='two.sided', conf.int = T)
# 0.03 -> statistical difference
# ci:  0.001195015 0.027812322
wilcox.test(meco_female, meco_male, paired=F, alternative='greater', conf.int = T)
# 0.017 -> MeCo is higher in female
# ci: 0.003274702         Inf


# MeCo regulation
df4 <- data.frame(MeCo_reg,gender)
meco_female <- df4$MeCo_reg[which(df4$gender =='female')]
meco_male <- df4$MeCo_reg[which(df4$gender =='male')]
shapiro.test(meco_female)
shapiro.test(meco_male)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_female, meco_male, names = c('female','male'), main='MeCo regulation distribution by gender')

# H0: No tendency of MeCo regulation score depending on the gender
# H1: Tendency of some type of MeCo regulation score depending on the gender
wilcox.test(meco_female, meco_male, paired=F, alternative='two.sided', conf.int = T)
# 1.23e-05 -> significant difference
# ci:  0.02239838 0.05730498
wilcox.test(meco_female, meco_male, paired=F, alternative='greater', conf.int = T)
# 6.152e-06 -> MeCo regulation is higher in female
# ci: 0.0252336         Inf

# MeCo stimulus
df4 <- data.frame(MeCo_st,gender)
meco_female <- df4$MeCo_st[which(df4$gender =='female')]
meco_male <- df4$MeCo_st[which(df4$gender =='male')]
shapiro.test(meco_female)
shapiro.test(meco_male)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_female, meco_male, names = c('female','male'), main='MeCo stimulus distribution by gender')

# H0: No tendency of MeCo stimulus score depending on the gender
# H1: Tendency of some type of MeCo regulation score depending on the gender
wilcox.test(meco_female, meco_male, paired=F, alternative='two.sided', conf.int = T)
# 0.00017 -> significant difference
# ci:  0.01647529 0.05265989
wilcox.test(meco_female, meco_male, paired=F, alternative='greater', conf.int = T)
# 8.44e-05 -> MeCo stimulus is higher in female
# ci: 0.01939692          Inf

# MeCo development
df4 <- data.frame(MeCo_dev,gender)
meco_female <- df4$MeCo_dev[which(df4$gender =='female')]
meco_male <- df4$MeCo_dev[which(df4$gender =='male')]
shapiro.test(meco_female)
shapiro.test(meco_male)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_female, meco_male, names = c('female','male'), main='MeCo development distribution by gender')

# H0: No tendency of MeCo development score depending on the gender
# H1: Tendency of some type of MeCo development score depending on the gender
wilcox.test(meco_female, meco_male, paired=F, alternative='two.sided', conf.int = T)
# 0.007 -> significant difference
# ci:  0.006276724 0.038185415
wilcox.test(meco_female, meco_male, paired=F, alternative='greater', conf.int = T)
# 0.003338 -> MeCo development is higher in female
# ci: 0.008859827          Inf



# FIFTH HYPOTHESIS TEST: Test, for each type of MeCO score, whether it changes depending on the laterality of the tumor

# general MeCo
df5 <- data.frame(MeCo,laterality)
meco_left <- df5$MeCo[which(df5$laterality =='Left')]
meco_right <- df5$MeCo[which(df5$laterality =='Right')]
shapiro.test(meco_left)
shapiro.test(meco_right)
# the p-values of the shapiro test are low so these two variables can't be assumed to 
# be normally distributed 

boxplot(meco_left, meco_right, names = c('Left','Right'), main='MeCo distribution by laterality')

# H0: No tendency of MeCo score depending on the laterality
# H1: Tendency of some type of MeCo score depending on the laterality
wilcox.test(meco_right, meco_left, paired=F, alternative='two.sided', conf.int = T)
# 0.2 -> no significant difference

# MeCo regulation
df5 <- data.frame(MeCo_reg,laterality)
meco_left <- df5$MeCo_reg[which(df5$laterality =='Left')]
meco_right <- df5$MeCo_reg[which(df5$laterality =='Right')]
shapiro.test(meco_left)
# this variable is not normally distributed
shapiro.test(meco_right)
# this variable is normally distributed

boxplot(meco_left, meco_right, names = c('Left','Right'), main='MeCo regultion distribution by laterality')

# H0: No tendency of MeCo regulation score depending on the laterality
# H1: Tendency of some type of MeCo regulation score depending on the laterality
wilcox.test(meco_right, meco_left, paired=F, alternative='two.sided', conf.int = T)
# 0.8 -> no significant difference


# MeCo stimulus
df5 <- data.frame(MeCo_st,laterality)
meco_left <- df5$MeCo_st[which(df5$laterality =='Left')]
meco_right <- df5$MeCo_st[which(df5$laterality =='Right')]
shapiro.test(meco_left)
# this variable is not normally distributed
shapiro.test(meco_right)
# this variable is normally distributed

boxplot(meco_left, meco_right, names = c('Left','Right'), main='MeCo stimulus distribution by laterality')

# H0: No tendency of MeCo stimulus score depending on the laterality
# H1: Tendency of some type of MeCo stimulus score depending on the laterality
wilcox.test(meco_right, meco_left, paired=F, alternative='two.sided', conf.int = T)
# 0.9 -> no significant difference


# MeCo development
df5 <- data.frame(MeCo_dev,laterality)
meco_left <- df5$MeCo_dev[which(df5$laterality =='Left')]
meco_right <- df5$MeCo_dev[which(df5$laterality =='Right')]
shapiro.test(meco_left)
# this variable is not normally distributed
shapiro.test(meco_right)
# this variable is normally distributed

boxplot(meco_left, meco_right, names = c('Left','Right'), main='MeCo development distribution by laterality')

# H0: No tendency of MeCo development score depending on the laterality
# H1: Tendency of some type of MeCo development score depending on the laterality
wilcox.test(meco_right, meco_left, paired=F, alternative='two.sided', conf.int = T)
# 0.1 -> no significant difference
