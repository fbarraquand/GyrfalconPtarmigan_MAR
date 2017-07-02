################################################################################################################
################# Analysis Gyrfalcon-Partmigan data - time series of abundance plots and MAR modelling #########
### FBarraquand 06/07/2015, with O. Nielsen 
### Starting with very basic plots (reproducing Nielsen 2011 report) and fitting MAR models
### Using abundance data only -> no fledging prod., no complex production of floaters or functional responses###
### Modified 15/11/2016, added some Granger causality testing. 
### Modified 14/02/2017, add weather data and simplify the MAR modelling
### Modified 07/03/2017, add weather data averaged over the whole area (3 stations for temp and 6 for rainfall)
### Title changed -- analyses cleaned up. 
### 17/03/2017 Using MARSS package to have better model selection criteria and take some values to zero?
### 23/03/2017 Add winter weather to ptarmigan growth analyses. 
################################################################################################################

### Initializing
rm(list=ls())
graphics.off()

### Reading data on gyr-ptarmigan
DGP<-read.csv("/home/frederic/Documents/MAR_modelling/Gyr/Gyrfalcon_Data.csv")
par(mfrow=c(2,1))
plot(DGP$Year,DGP$Occupied,type="b")
plot(DGP$Year,DGP$Ptarmigan,type="b")
## Need to correct the gyr portion for observational effort
DGP$Occupancy=DGP$Occupied/DGP$N

par(mfrow=c(2,1))
plot(DGP$Year,DGP$Occupancy,type="b",main="Percentage territories occupied Gyrfalcon")
plot(DGP$Year,DGP$Ptarmigan,type="b",main="Mean density Ptarmigan")

#Plot both standardized
png(file="GyrPtarDensities.png",width=10,height=8,res=300,units="in")
DGP$OccStd=(DGP$Occupancy-mean(DGP$Occupancy))/sd(DGP$Occupancy)
DGP$PtarStd=(DGP$Ptarmigan-mean(DGP$Ptarmigan))/sd(DGP$Ptarmigan)
par(mfrow=c(1,1),cex=1.5)
plot(DGP$Year,DGP$OccStd,type="b",col="red",ylim=c(-3,3),ylab = "Stdized population density",xlab="Year",lwd=3,pch=20)
lines(DGP$Year,DGP$PtarStd,type="b",lwd=3,pch=20)
dev.off()

##############################################################################################
### Easier to interpret the coefficients once the data is logged and standardized. Do that now. 
################# Standardized data analysis ##############################################
xbis=matrix(0,nrow=2,ncol=nrow(DGP))
xbis[1,]<-log(DGP$Ptarmigan)#1 is prey
xbis[2,]<-log(DGP$Occupancy)#2 is predator
m1=mean(xbis[1,])
m2=mean(xbis[2,])
s1=sd(xbis[1,])
s2=sd(xbis[2,])
xbis[1,]<-(xbis[1,]-m1)/s1
xbis[2,]<-(xbis[2,]-m2)/s2
growth_rate = xbis[,2:nrow(DGP)]-xbis[,1:(nrow(DGP)-1)]
x=xbis[,1:(nrow(DGP)-1)]

############################################################################################
###### Simple analyses using lm() and Granger tests - without climate data for now
############################################################################################

############ Loading packages to investigate Granger causality 
library(lmtest)
library(vars)

###########  Causal approach with several lags. 
varpp<-VAR(y=data.frame(t(xbis)), type="none",lag.max=5)
## Considering 5 maximum lags and using VAR(p) estimation with "vars" package (MAR(p) in ecology)
varpp #Yields model order=3
### Note this is somewhat consistent with previous results, at least concerning AIC (see below for BIC)
## http://raunvisindastofnun.hi.is/sites/raunvisindastofnun.hi.is/files/rh-18-2003.pdf
# I also tried a lag.max = 20, gives lag 7
# clearly this is overparameterized and would need some more work on model selection

###### Looking at several model selection criteria
var_order_select=VARselect(y=data.frame(t(xbis)), type="none",lag.max=5)
var_order_select

### Now compute AICc and BIC for completeness and coherence with ecological standards
AIC_values=var_order_select$criteria[1,] ## Are these row AIC values or Delta AIC?
### Looks like they are divided by the sample size perhaps?
# https://www.rdocumentation.org/packages/vars/versions/1.5-2/topics/VARselect
n=c(33,32,31,30,29) ## different data length for the different models MAR(1), MAR(2), ...
# Should I take it into account? 
34*AIC_values ### still not there... 
# May look like the AIC values also lack a 2 in front of the log-likelihood, 
# but actually this might just be because we look at the error here and not the LL. 
# Conclusion: there is a division by TS length here and there are also terms from the LL not shown
# (so this is not true AIC... and not comparable across packages)
# Compute AICc and BIC
T=34
k_param=c(4,2*4,3*4,4*4,5*4)
BIC_values = T*AIC_values + log(T)*k_param-2*k_param
BIC_values
BIC_values/T ### this is the SC(n) 

# We select model 1 there, if I am not mistaken. 
AICc_values = T*AIC_values + 2*(k_param+1)*(k_param+2)/(T-k_param-2)
### Model 2 is consistently chosen once AICc and BIC are considered. 

####################### More analyses ################################################
########### Checking results with lm() // similar but assumes diagonal error matrix
lm1=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,])# pb
lm2=lm(growth_rate[2,] ~ 0 + x[1,]+x[2,])
##########
lm1
# Call:
#   lm(formula = growth_rate[1, ] ~ 0 + x[1, ] + x[2, ])
# 
# Coefficients:
#   x[1, ]   x[2, ]  
# -0.2290  -0.2333  
lm2
# Call:
#   lm(formula = growth_rate[2, ] ~ 0 + x[1, ] + x[2, ])
# 
# Coefficients:
#   x[1, ]   x[2, ]  
# 0.2150  -0.3399  
########## Similar to MARSS and other analyses later on ##################

#### Testing Granger causality at several lags, using yet another algorithm
grangertest(x[2,],x[1,],order = 1) # Effect of predator on prey
grangertest(x[1,],x[2,],order = 1) # Effect of prey on predator

grangertest(x[2,],x[1,],order = 2) # "Pb" here at order 2 -- no more effect of predator on prey
grangertest(x[1,],x[2,],order = 2) # Still an effect of prey on predator. 

grangertest(x[2,],x[1,],order = 3) ### Pb here at order 3 as well. 
grangertest(x[1,],x[2,],order = 3)

### Not 100% sure if the causal effect of predator on prey is weak or the MAR(1) model is just the right representation. 
### This can be tested by fit and / or simulation of the dynamics. Could be / will be done in JAGS too if need be. 
### Testing the simulated dynamics of the MAR(2) model will probably suffice here, given higher-order lags are likely overparameterized. 
### The MAR(2) model is already estimating 4*4 autoregressive coefficients on a 34 steps time series. 

######## MAR(2) model #################
varpp2<-VAR(y=data.frame(t(xbis)), p=2, type="none")
causality(varpp2,cause="X2") ## check this. 

######## MAr(1) model ################# 
varpp1<-VAR(y=data.frame(t(xbis)), p=1, type="none") # Looks better but...

### Extracting coefficients of MAR(2) model for simulations 
### (I later put to zero non-significant coefficients, but not for now)
coef(varpp2)
coef(varpp2)$X1
coef(varpp2)$X1[1,]
a_11_l1=coef(varpp2)$X1[1,1]
a_12_l1=coef(varpp2)$X1[2,1]
a_11_l2=coef(varpp2)$X1[3,1]
a_12_l2=coef(varpp2)$X1[4,1]

a_21_l1=coef(varpp2)$X2[1,1]
a_22_l1=coef(varpp2)$X2[2,1]
a_21_l2=coef(varpp2)$X2[3,1]
a_22_l2=coef(varpp2)$X2[4,1]

### Diagnostics of MAR(2) model
summary(varpp2)
residuals(varpp2)
acf(residuals(varpp2))
varpp2$varresult
Sigma=summary(varpp2)$covres ### Extracting variance-covariance matrix

################ Simulations of our fitted models ##################################################
#### Full MAR(2) model
z=matrix(0,nrow=2,ncol=nrow(DGP))
z[,1]=runif(2,0,1)
z[,2]=runif(2,0,1)
z
n=nrow(DGP)
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + a_12_l1*z[2,t] + a_11_l2*z[1,t-1] + a_12_l2*z[2,t-1]+eps[1]
  z[2,t+1] = a_21_l1*z[1,t] + a_22_l1*z[2,t] + a_21_l2*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}
matplot(t(z))
matlines(t(z))

### Not that bad...
### Often species 2 follows species 1, as in the data. 

### Now let's set to zero these non-significant coeff to understand what's going on
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + 0*z[2,t] + a_11_l2*z[1,t-1] + 0*z[2,t-1]+eps[1]
  z[2,t+1] = a_21_l1*z[1,t] + a_22_l1*z[2,t] + 0*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}
matplot(t(z))
matlines(t(z))
### A little less realistic, clearly

### Two separate AR(2) models
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + 0*z[2,t] + a_11_l2*z[1,t-1] + 0*z[2,t-1]+eps[1]
  z[2,t+1] = 0*z[1,t] + a_22_l1*z[2,t] + 0*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}
matplot(t(z))
matlines(t(z))

### Keep the delayed effect of the predator on the prey even though non-significant
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + 0*z[2,t] + a_11_l2*z[1,t-1] + a_12_l2*z[2,t-1]+eps[1]
  z[2,t+1] = 0*z[1,t] + a_22_l1*z[2,t] + 0*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}
matplot(t(z))
matlines(t(z)) # not enough

### Keep the delayed effect of the prey on the predator even though non-significant
### With the prey that is internally driven as an AR(2) model. 
for (t in 2:(n-1)){
  eps=mvrnorm(mu=c(0,0),Sigma=Sigma)
  z[1,t+1] = a_11_l1*z[1,t] + 0*z[2,t] + a_11_l2*z[1,t-1] + 0*z[2,t-1]+eps[1]
  z[2,t+1] = 0*z[1,t] + a_22_l1*z[2,t] + a_21_l2*z[1,t-1] + a_22_l2*z[2,t-1]+eps[2]
}
matplot(t(z))
matlines(t(z))
### That seems to be trick, i.e. the element in the model that yields realistic-looking simulations. 
### Growth rate of predator depends on the prey a year before. 
### This could be an artefact, or it could also be a genuine mechanisms with age structure in the predator population.
### This MAR(2) model with cyclic prey and bottom-up driven predator growth is overparameterized, as not selected by AICc/BIC
### But this model should be compared at some point to the MAR(1) with a strong penalty for parameter numbers. 
### Hereafter I will refer to this model as the "bottom-up MAR(2) model". 

###################################################################################
#### Uses MARSS package to check the results produced by "vars"
###################################################################################
# ### Estimation with MARSS package
# library('MARSS')
# ### See e.g http://cran.r-project.org/web/packages/MARSS/vignettes/Quick_Start.pdf
# ### Setting matrices
# B1=matrix(list("b11","b21","b12","b22"),2,2)
#### WARNING! #####################################################################
### In the vignette B1=matrix(list("b11","b12","b21","b22"),2,2)
### There was an error in the vignette here because R fills matrices by columns....
###################################################################################
# U1=matrix(0,2,1)
# Q1=matrix(c("q11",0,0,"q22"),2,2) ##assume no correlated noise
# Z1=matrix(c(1,0,0,1),2,2)
# A1=matrix(list(0,0),2,1)
# R1=matrix(list(0,0,0,0),2,2)
# pi1=matrix(0,2,1); V1=diag(1,2)
# # Estimation
# data<-xbis
# model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
# fit=MARSS(data, model=model.list)
# MARSSparamCIs(fit)
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -65.92927 
# AIC: 145.8585   AICc: 147.7252   
# 
# ML.Est Std.Err   low.CI   up.CI
# B.b11  0.7710  0.1115  0.55249  0.9894
# B.b21  0.2150  0.1091  0.00113  0.4289
# B.b12 -0.2333  0.1116 -0.45207 -0.0146
# B.b22  0.6601  0.1092  0.44596  0.8742
# Q.q11  0.3961  0.0776  0.23297  0.6306
# Q.0    0.0436  0.0542 -0.06778  0.2178
# Q.q22  0.3844  0.0752  0.22474  0.6029
# Previously the matrix was c(0.77, 0.22,-0.24,0.66) which is similar. The Q matrix looks different although these are variances...

### Try correlated noise 
# B1=matrix(list("b11","b21","b12","b22"),2,2)
#### WARNING! #####################################################################
### In the vignette B1=matrix(list("b11","b12","b21","b22"),2,2)
### There was an error in the vignette here because R fills matrices by columns....
###################################################################################
# Q1=matrix(c("q11","q21","q12","q22"),2,2) ##assume correlated noise
# model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)#did not work
# model.list=list(B=B1,U=U1,Q="unconstrained",Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
# fit2=MARSS(data, model=model.list)
# MARSSparamCIs(fit2)
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -65.92927 
# AIC: 145.8585   AICc: 147.7252   
# 
# ML.Est Std.Err   low.CI   up.CI
# B.b11    0.7710  0.1115  0.55249  0.9894
# B.b21    0.2150  0.1091  0.00113  0.4289
# B.b12   -0.2333  0.1116 -0.45207 -0.0146
# B.b22    0.6601  0.1092  0.44596  0.8742
# Q.(1,1)  0.3961  0.0776  0.23297  0.6306
# Q.(2,1)  0.0436  0.0542 -0.06778  0.2178
# Q.(2,2)  0.3844  0.0752  0.22474  0.6029
### seems to be the exact same thing, with a covariance given previously as well (which is weak)...

##############################################################################
## Model with no interactions
# B1ter=matrix(list("b11",0,0,"b22"),2,2)
# B1ter #OK
# model.list=list(B=B1ter,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
# fit3=MARSS(data, model=model.list)
# MARSSparamCIs(fit3)
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -69.71731 
# AIC: 149.4346   AICc: 150.4024   
# 
# ML.Est Std.Err  low.CI up.CI
# B.b11 0.6991  0.1222  0.4597 0.939
# B.b22 0.7301  0.1190  0.4969 0.963
# Q.q11 0.4497  0.0829  0.2640 0.724
# Q.0   0.0657  0.0638 -0.0762 0.285
# Q.q22 0.4302  0.0788  0.2498 0.675
# 
# CIs calculated at alpha = 0.05 via method=hessian 
### So the residual variances are higher without the interactions again and AIC/AIC lower. 

#### Compute and compare the BICs ############################
AIC1=145.8585
AIC2=145.8585 ### equal because only free parameters considered?
AIC3=149.4346 
k1=4 # Because it is a 2x2 matrix
k2=4
k3=2
BIC1=AIC1+2*log(ncol(xbis))*k1-2*k1
BIC2=AIC2+2*log(ncol(xbis))*k2-2*k2
BIC3=AIC3+2*log(ncol(xbis))*k3-2*k3
BIC1
BIC2
BIC3 # Prefers the third model without interaction. 
### Check a few things with Log-likelihood. 
LL1=-65.92927 
LL2=-65.92927 
LL3=-69.71731 
(-2)*LL1+2*k1
(-2)*LL2+2*k2
(-2)*LL3+2*k3
#WTF? Looks like we don't have the right number of parameters
(-2)*LL1+2*(k1+3)
(-2)*LL2+2*(k2+3)
(-2)*LL3+2*(k3+3)
### That's because they count the var-covar parameters and also the initial values. 
### NB the VAR formulation of AIC in vars looks like what Ives et al. 2003 had
### https://www.rdocumentation.org/packages/vars/versions/1.5-2/topics/VARselect
### https://courses.maths.ox.ac.uk/node/view_material/924
### The AIC in vars is therefore divided by sample size. 
### See http://stats.stackexchange.com/questions/191531/different-aic-definitions

################################################################################################
########### Let's now assume the MAR(1) model is OK ############################################
########### It has been selected by both AICc and BIC ##########################################
########### More model fitting in MARSS is required though #####################################

### Let's stay in the MAR(1) / ARX framework and add in the weather 
### We will later contrast the MAR(1) framework to the "bottom-up" MAR(2) 

### Adding in the new weather data, averaged for log(rainfall) and temperature.   
# Previously using 
#DB=read.table("weather_iceland/Stod_495_Grimsstadir.ManMedal.txt",header=T,encoding = "latin1")
#DB2=read.table("weather_iceland/Stod_422_Akureyri.ManMedal.txt",header=T,encoding = "latin1")
#DB3=read.table("weather_iceland/Stod_479_Manarbakki.ManMedal.txt",header=T,encoding = "latin1")
#names(DB)[1:3]=c("site","year","month")
#names(DB2)[1:3]=c("site","year","month")
#names(DB3)[1:3]=c("site","year","month")

### Loading in 
DB=read.csv("weather_iceland/average_weatherNEiceland.csv",header=T)
names(DB)[6]="r" #renaming logRainfall
### Years considered for weather
DGP$Year[1:(nrow(DGP)-1)] # 1981 to 2013, since the last growth rate is arriving in 2014

### Weather data for prey
year_minus_1=DGP$Year[1:(nrow(DGP)-1)]-1
# Average temperature 1 year earlier in May 
tempMay_year_minus1=DB$temp[(DB$year %in% year_minus_1)&(DB$month==5)]
# Log-rainfall 1 year earlier in May
rainMay_year_minus1=log(DB$r[(DB$year %in% year_minus_1)&(DB$month==5)])

### Weather data for predator
year_minus_4=DGP$Year[1:(nrow(DGP)-1)]-4
# Average temperature 4 years earlier in April
tempApril_year_minus4=DB$temp[(DB$year %in% year_minus_4)&(DB$month==4)]
# Log-rainfall 4 years earlier in April
rainApril_year_minus4=log(DB$r[(DB$year %in% year_minus_4)&(DB$month==4)])

### Standardize all those variables to be able to compare something
tempMay_year_minus1=(tempMay_year_minus1-mean(tempMay_year_minus1))/sd(tempMay_year_minus1)
tempApril_year_minus4=(tempApril_year_minus4-mean(tempApril_year_minus4))/sd(tempApril_year_minus4)
rainMay_year_minus1=scale(rainMay_year_minus1,scale=TRUE) #tired of writing
rainApril_year_minus4=scale(rainApril_year_minus4,scale=TRUE)

lm1=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,])# 
lm2=lm(growth_rate[2,] ~ 0 + x[1,]+x[2,])

lm1_null=lm(growth_rate[1,] ~ 0 + x[1,])# 
lm2_null=lm(growth_rate[2,] ~ 0 + x[2,])

lm1_temp=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempMay_year_minus1)# 
lm1_temp # no effect of temp in May would be logical but... 
lm2_temp=lm(growth_rate[2,] ~ 0 + x[1,]+x[2,]+tempApril_year_minus4)
lm2_temp # positive effect of temp in April on pred growth with a 4-year timelag

lm1_temp_only=lm(growth_rate[1,] ~ 0 + x[1,]+tempMay_year_minus1)# 
lm1_temp # no effect of temp in May would be logical but... 
lm2_temp_only=lm(growth_rate[2,] ~ 0 + x[2,]+tempApril_year_minus4)
lm2_temp # positive effect of temp in April on pred growth with a 4-year timelag

lm1_rain=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+rainMay_year_minus1)# 
lm1_rain # at least the effect is smallish and of logical sign
lm2_rain=lm(growth_rate[2,] ~ 0 + x[1,]+x[2,]+rainApril_year_minus4)
lm2_rain

lm1_rain_only=lm(growth_rate[1,] ~ 0 + x[1,]+rainMay_year_minus1)# 
lm1_rain # at least the effect is smallish and of logical sign
lm2_rain_only=lm(growth_rate[2,] ~ 0 + x[2,]+rainApril_year_minus4)
lm2_rain

lm1_2cov=lm(growth_rate[1,] ~ 0 + x[1,]+tempMay_year_minus1+rainMay_year_minus1)# 
lm1_2cov # no effect of temp in May logical
lm2_2cov=lm(growth_rate[2,] ~ 0 + x[2,]+tempApril_year_minus4+rainApril_year_minus4)
lm2_2cov # there could be interactions but unclear if we should try. 

lm1_full=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempMay_year_minus1+rainMay_year_minus1)# 
lm1_full # no effect of temp in May logical
lm2_full=lm(growth_rate[2,] ~ 0 + x[1,]+x[2,]+tempApril_year_minus4+rainApril_year_minus4)
lm2_full # there could be interactions but unclear if we should try. 

AIC(lm1_null,lm1,lm1_temp,lm1_rain,lm1_full) #lm1 winner
AIC(lm2_null,lm2,lm2_temp,lm2_rain,lm2_full) #lowest AIC for lm2_temp
# But all of these are very very close - we will probably need different criteria to contrast those. 

####### Checking coefficients associated to the various models
summary(lm1_rain)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# x[1, ]              -0.22340    0.11769  -1.898   0.0673 .
# x[2, ]              -0.22990    0.11708  -1.964   0.0589 .
# rainMay_year_minus1  0.04011    0.11840   0.339   0.7371  
summary(lm1_temp)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# x[1, ]               -0.2554     0.1191  -2.145   0.0402 *
#   x[2, ]               -0.2028     0.1205  -1.683   0.1027  
# tempMay_year_minus1  -0.1081     0.1227  -0.881   0.3855  
### There might a minor though NS effect of temperature in May on ptarmigan -- however it has the wrong sign... 

### What if we consider that interactions might not be necessary?
AIC(lm1_null,lm1,lm1_temp,lm1_temp_only,lm1_rain,lm1_rain_only,lm1_full,lm1_2cov) #lm1 winner, clearly lm1_2cov the worst
AIC(lm2_null,lm2,lm2_temp,lm2_temp_only,lm2_rain,lm2_rain_only,lm2_full,lm2_2cov) #lowest AIC for lm2_temp
summary(lm2_temp)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# x[1, ]                  0.1830     0.1081   1.693  0.10085   
# x[2, ]                 -0.3177     0.1077  -2.950  0.00611 **
# tempApril_year_minus4   0.2323     0.1084   2.143  0.04032 * 
#### Positive and non-negligible effect of lagged 4 years ago temperature. 

### Check the view of BIC
BIC(lm1_null,lm1,lm1_temp,lm1_temp_only,lm1_rain,lm1_rain_only,lm1_full,lm1_2cov) #lm1 winner, clearly lm1_2cov and lm1_full the worst
BIC(lm2_null,lm2,lm2_temp,lm2_temp_only,lm2_rain,lm2_rain_only,lm2_full,lm2_2cov) #lowest BIC for lm2_temp_only but not far from lm2_temp
summary(lm2_temp_only)# won't yield correct dynamics. 

library(AICcmodavg) ## for AICc
### AICc ##################
AICc(lm1) # best
AICc(lm1_null) # not so bad
AICc(lm1_temp)
AICc(lm1_temp_only) #worst
AICc(lm1_2cov)

### AICc ##################
AICc(lm2) # OK
AICc(lm2_null) # 
AICc(lm2_temp) # best
AICc(lm2_temp_only) #a little off only - no good model in other respects though. 
AICc(lm2_2cov)
AICc(lm2_full)

####### Clearly there should be reciprocal feedback + some effect of delayed temp on gyr pop growth rate. 
####### No effect of temperature of the year before in May on ptarmigan dynamics. 

########### More analyses regarding ptarmigan that seems poorly characterized so far. 
### Check the effect of wind speed on ptarmigan -- not now because average wind speed not useful. 
### We would need to get daily data and e.g. reconstruct the effect of snowstorms on chicks. 

### But perhaps no timelag, i.e. weather of the year is informative as well? 
year_now=DGP$Year[1:(nrow(DGP)-1)]
# Average temperature  in May 
tempMay_year=DB$temp[(DB$year %in% year_now)&(DB$month==5)]
# Log-rainfall in May
rainMay_year=log(DB$r[(DB$year %in% year_now)&(DB$month==5)])
### windMay_year=log(DB$f[(DB$year %in% year_now)&(DB$month==5)]) #Not considered anymore... 

lm1_temp0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempMay_year)# 
lm1_rain0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+rainMay_year)# 
#lm1_wind0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+windMay_year)# 
lm1_full0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempMay_year+rainMay_year)# 

AIC(lm1,lm1_null,lm1_temp0,lm1_rain0,lm1_full0)
# All of these have higher AIC than lm1. But the null is the worst. 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# x[1, ]       -0.228751   0.116780  -1.959   0.0595 .
# x[2, ]       -0.234391   0.117740  -1.991   0.0557 .
# tempMay_year  0.001859   0.025102   0.074   0.9415
summary(lm1_temp0)
# Now what does happen with AICc and BIC here? 

BIC(lm1,lm1_null,lm1_temp0,lm1_rain0,lm1_full0)
### BIC also favors lm1 - not by much compared to the null though

AICc(lm1) ## still number 1
AICc(lm1_null)
AICc(lm1_temp0)
AICc(lm1_rain0)
AICc(lm1_full0)

################### Now implement the suggestion of Olafur - looking at July's data #######################################
### The growth rate from year t to t+1 (in spring) is considered. Hence if chick survival is a key factor, it has to be 
### July from the year t (or t-1, but let's start with t, assuming the chicks contribute immediately to pop growth). 

### Creating new weather variables
tempJuly_year=DB$temp[(DB$year %in% year_now)&(DB$month==7)]
# Log-rainfall in July
rainJuly_year=log(DB$r[(DB$year %in% year_now)&(DB$month==7)])

### 
lm1_temp0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempJuly_year)# 
lm1_rain0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+rainJuly_year)# 
#lm1_wind0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+windMay_year)# 
lm1_full0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempJuly_year+rainJuly_year)# 

AIC(lm1,lm1_temp0,lm1_rain0,lm1_full0)
summary(lm1_temp0) # very weak effect of July temp and not significant 
### I tend to conclude there is no effect of July temp on ptarmigan pop growth. 

### June weather

### Creating new weather variables
tempJune_year=DB$temp[(DB$year %in% year_now)&(DB$month==6)]
# Log-rainfall in June
rainJune_year=log(DB$r[(DB$year %in% year_now)&(DB$month==6)])
### Fuck I did not scale my variables... do that later. 

### 
lm1_temp0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempJune_year)# 
lm1_rain0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+rainJune_year)# 
#lm1_wind0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+windMay_year)# 
lm1_full0=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+tempJune_year+rainJune_year)# 

AIC(lm1,lm1_temp0,lm1_rain0,lm1_full0)
summary(lm1_temp0) # very weak effect of June temp and not significant 
### I tend to conclude there is no effect of June temp on ptarmigan pop growth. 

### New idea (cf. discussion Olafur 21/03/2017), to look at winter weather (we talked about this before, then dismissed it)
### Oli pointed out that the age ratios might hint at low winter survival for the chicks
### Need some kind of proxy of winter severity. 

### Creating new weather variables
December_temp=DB$temp[(DB$year %in% year_minus_1)&(DB$month==12)]
January_temp = DB$temp[(DB$year %in% year_now)&(DB$month==1)] 
February_temp = DB$temp[(DB$year %in% year_now)&(DB$month==2)] 
March_temp = DB$temp[(DB$year %in% year_now)&(DB$month==3)] 
winter_temp=cbind(December_temp,January_temp,February_temp,March_temp)
## Use monthly average and min month temp
avMonthly_tempWinter=rowMeans(winter_temp)
avMonthly_tempWinter=scale(avMonthly_tempWinter)
minOfMonths_tempWinter=pmin(December_temp,January_temp,February_temp,March_temp)
minOfMonths_tempWinter=scale(minOfMonths_tempWinter)

matplot(winter_temp)
matlines(winter_temp)


## Rainfall can be cumulated over all winter (and then logged)
December_rain=DB$r[(DB$year %in% year_minus_1)&(DB$month==12)]
January_rain = DB$r[(DB$year %in% year_now)&(DB$month==1)] 
February_rain = DB$r[(DB$year %in% year_now)&(DB$month==2)] 
March_rain = DB$r[(DB$year %in% year_now)&(DB$month==3)] 
winter_rain=cbind(December_rain,January_rain,February_rain,March_rain)
matplot(winter_rain)
matlines(winter_rain)
log_winter_rain = log(rowSums(winter_rain)) ## of course some of that rain is snow
log_winter_rain = scale(log_winter_rain)
plot(log_winter_rain,type="o")

lm1_temp_mean=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+avMonthly_tempWinter)# 
lm1_temp_min=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+minOfMonths_tempWinter)# 
lm1_total_rain=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+log_winter_rain)# 
lm1_full=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+minOfMonths_tempWinter+log_winter_rain)# 

AIC(lm1_temp_mean,lm1_temp_min,lm1_total_rain,lm1_full,lm1,lm1_null)
summary(lm1_temp_mean) # very weak effect of temp and not significant
summary(lm1_temp_min) # very weak effect of temp and not significant
summary(lm1_total_rain) # very weak effect of rain and not significant
# just to check
lm1_interaction=lm(growth_rate[1,] ~ 0 + x[1,]+x[2,]+minOfMonths_tempWinter:log_winter_rain)# 
summary(lm1_interaction) #nothing there as well. 


############## 02/05/ 2017 - Add analyses using varpp ##################################################

### Loading in 
DB=read.csv("weather_iceland/average_weatherNEiceland.csv",header=T)
names(DB)[6]="r" #renaming logRainfall
### Years considered for weather
DGP$Year[1:(nrow(DGP))] # 1981 to 2014

### Weather data for prey
year_minus_1=DGP$Year[1:(nrow(DGP))]-1
# Average temperature 1 year earlier in May 
tempMay_year_minus1=DB$temp[(DB$year %in% year_minus_1)&(DB$month==5)]
# Log-rainfall 1 year earlier in May
rainMay_year_minus1=log(DB$r[(DB$year %in% year_minus_1)&(DB$month==5)])

### Weather data for predator
year_minus_4=DGP$Year[1:(nrow(DGP))]-4
# Average temperature 4 years earlier in April
tempApril_year_minus4=DB$temp[(DB$year %in% year_minus_4)&(DB$month==4)]
# Log-rainfall 4 years earlier in April
rainApril_year_minus4=log(DB$r[(DB$year %in% year_minus_4)&(DB$month==4)])

exogen_vec<-cbind(tempMay_year_minus1,tempApril_year_minus4)
varpp2w<-VAR(y=data.frame(t(xbis)), p=2, type="none",exogen=exogen_vec)
varpp2w
coef(varpp2w)
summary(varpp2w) ### NB the weather for the predator is only vaguely significant. 

### check this is not a likely difference of years, since previously the growth rate was dependent upon weather for years ago
### Weather data for predator
year_minus_5=DGP$Year[1:(nrow(DGP))]-5
# Average temperature 4 years earlier in April
tempApril_year_minus5=DB$temp[(DB$year %in% year_minus_5)&(DB$month==4)]
# Log-rainfall 4 years earlier in April
rainApril_year_minus5=log(DB$r[(DB$year %in% year_minus_5)&(DB$month==4)])

exogen_vec<-cbind(tempMay_year_minus1,tempApril_year_minus5)
varpp2w<-VAR(y=data.frame(t(xbis)), p=2, type="none",exogen=exogen_vec)
varpp2w
coef(varpp2w)
summary(varpp2w) ### it is not. 

causality(varpp2w,cause="X1")
causality(varpp2w,cause="X2")

###########################################################################################################################
########## More complex model fitting analyses in JAGS (Bayesian)                                           ################
########## They show almost identical results to lm() and VAR for the MAR(1) model. 
########## MAR(2) models are more complicated in MARSS (but MARSS provide AIC, BIC, AICc more easily...)
########## We could compare MAR(1) and bottom-up MAR(2) in JAGS later on though. 
########## Predictive abilities might not be easier to evaluate though. Bayesian P-values have their limits. 
###########################################################################################################################

jags.data<-list(t_max=ncol(xbis),x=xbis)

library("R2jags")      # Load R2jags package
{
sink("MAR_GyrPtar.txt")
cat("
    model 
{
    ### Definition of priors
    b11 ~ dunif(-1,1) # regulation prey (add possibility for overcompensation here)
    b21 ~ dunif(0,1) # effect prey on predator
    b22 ~ dunif(0,1) # regulation predator (no overcompensation assumed)
    b12 ~ dunif(-1,0) # effect predator on prey. Assumed negative. 
    
    ### Priors for the matrix noise parameters - multiple possible choices here. 
    #tau[1] ~ dgamma(.01,.01)
    #tau[2] ~ dgamma(.01,.01)
    for (i in 1:2){
    sigma[i] ~ dunif(0,1)
    tau[i]<-pow(sigma[i],-2)
    }
    
    ### B matrix definition
    B[1,1]<-b11
    B[1,2]<-b12
    B[2,1]<-b21
    B[2,2]<-b22
    
    ### Priors for initial values of x. No latent states -> no need to initialize them?
    x10 ~ dnorm(0,.01) #Those are on a log-scale, remember...
    x20 ~ dnorm(0,.01)
    
    x[1,1] ~ dnorm(x10,tau[1]) ## put the same precision as for the growth rate for the initial log-abundance
    x[2,1] ~ dnorm(x20,tau[2])
    
    #for (i in 1:2){r0[i] ~ dnorm(1,0.01)} ## prior on intercept of growth rate
    
    ### Loop over the years
    
    for (t in 1:(t_max-1)){
    for (i in 1:2){
    x[i,t+1] ~ dnorm(mu[i,t],tau[i]) ### 
    mu[i,t]<- B[i,1]*x[1,t]+B[i,2]*x[2,t] ### Line-column product
    x.rep[i,t+1] ~ dnorm(mu[i,t],tau[i])
    #xsim[i,t+1] <- x.rep[i,t+1] #can't store it for some reason
    }
    }
    # Compute Bayesian p-value
    for (i in 1:2){
    d_rep[i]<-sum(pow(x.rep[i,2:t_max]-mu[i,],2))/abs(mean(mu[i,]))
    d_data[i]<-sum(pow(mu[i,]-x[i,2:t_max],2))/abs(mean(mu[i,]))
    d_error[i]<-sum(pow(mu[i,]-x[i,2:t_max],2))/pow(sigma[i],2) 
    ## Latter distance from Ives et al. 2003 Likelihood for indep. error cov matrix - whether this can sensibly made into a R2 remains to be demonstrated. 
    rss[i]<-sum(pow(mu[i,]-x[i,2:t_max],2)) #Should be directly related to sigma through averaging...
    total_var[i]<-sum(pow(mean(x[i,2:t_max])-x[i,2:t_max],2))
    r2[i]<-rss[i]/total_var[i]
    # percentage of variance explained per species, still informative
    }
    D_rep<-sum(d_rep[])
    D_data<-sum(d_data[])
    
}
    
    ",fill=TRUE)
sink()
}

inits <- function(){list(b11=runif(1,0.5,1),b22=runif(1,0.5,1),b12=runif(1,-1,0),b21=runif(1,0,0.5),sigma=runif(2,0,1),
                         x10=rnorm(1,0,1),x20=rnorm(1,0,1))}

parameters <- c("b11","b22","b21","b12","sigma","x10","x20","D_rep","D_data","rss","r2","mu")#"xsim" does not work
# above line to get mu
#parameters <- c("b11","b22","b21","b12","sigma","x10","x20","x40","D_rep","D_data","rss","r2","r0")

# MCMC settings
nc <- 3 #number of chains
nb <- 10000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-20000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "MAR_GyrPtar.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())

print(out2, dig = 2)
#sigma[1]     0.67
#sigma[2]     0.66

jags.sum<-out2$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_GyrPtarData_stdized.txt")

## Bayesian P-value
par(mfrow=c(1,1))
mean(out2$BUGSoutput$sims.list$D_rep>out2$BUGSoutput$sims.list$D_data)
plot(out2$BUGSoutput$sims.list$D_rep,out2$BUGSoutput$sims.list$D_data)
abline(0,1)
##dev.off()

### Check predicted dynamics
mu=out2$BUGSoutput$mean$mu[,]
#xsim=out1$BUGSoutput$mean$xsim[,]
matplot(t(xbis),type="b")
matlines(t(mu),type="b") #very close to the current point... A good plot puts the prediction one steps ahead, see below

#Better plots including predictions one step ahead
par(mfrow=c(2,1))
plot(DGP$Year,xbis[1,],type="o",pch = 16,bg="black",xlab="Year",ylab="Ptarmigan density")
lines(DGP$Year[2:34],mu[1,],type="p")
arrows(DGP$Year[1:33],xbis[1,1:33],DGP$Year[2:34],mu[1,],length = 0.05)
plot(DGP$Year,xbis[2,],type="o",col="red",pch = 16,bg = "red",xlab="Year",ylab="Gyrfalcon log occupancy")
lines(DGP$Year[2:34],mu[2,],type="p",col="red")
arrows(DGP$Year[1:33],xbis[2,1:33],DGP$Year[2:34],mu[2,],col="red",length = 0.05)

#Better plots including predictions - scaled back to the normal scale
par(mfrow=c(2,1))
plot(DGP$Year,exp(xbis[1,]*s1+m1),type="o",pch = 16,bg="black",xlab="Year",ylab="Ptarmigan density")
lines(DGP$Year[2:34],exp(mu[1,]*s1+m1),type="p")
arrows(DGP$Year[1:33],exp(xbis[1,1:33]*s1+m1),DGP$Year[2:34],exp(mu[1,]*s1+m1),length = 0.05)
plot(DGP$Year,exp(xbis[2,]*s2+m2),type="o",col="red",pch = 16,bg = "red",xlab="Year",ylab="Gyrfalcon occupancy")
lines(DGP$Year[2:34],exp(mu[2,]*s2+m2),type="p",col="red")
arrows(DGP$Year[1:33],exp(xbis[2,1:33]*s2+m2),DGP$Year[2:34],exp(mu[2,]*s2+m2),col="red",length = 0.05)

### Now simulate dynamics
require(MASS)
library(Matrix)
t_max=34 # like in the data
x=matrix(0,nrow=2,ncol=t_max)
x[,1]=c(0,0) # we assume reasonable starting values based on the observed
#A=c(0.12,-0.28) #stdized now!!
mu=c(0,0)

### B matrix, point estimates
B=matrix(c(0.77, 0.22,-0.24,0.66), nrow=2,ncol=2)
lam_values=eigen(B)$values
max(abs(lam_values)) ### Spectral radius (not Re(lambda) here)
### Stable

## Noise 
Sigma = Diagonal(2, x = c(0.67,0.66)) ## Less variability on 
pdf("Simulated_MAR_dynamics_stdized.pdf",height=28,width=14)
par(mfrow=c(4,2))#,cex=1.5
plot(1:t_max,xbis[1,],type="o",col="black",ylim=c(-3,3),ylab="Real data")#ylim=c(-3,3)
lines(1:t_max,xbis[2,],type="o",col="red")
ccf(x[1,],x[2,],ylab = "cross-correlation")
for(nrep in 1:7)
{
  ### Initial values
  x[,1]=c(-0.8,-1.8)
  for (t in 1:(t_max-1))
  {
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B %*% x[,t] + epsilon
  }
  
  #plot(1:t_max,x[1,],type="b",col="black",ylim=c(-3,3),ylab="Unstandardized log densities")
  #lines(1:t_max,x[2,],type="b",col="red")
  # 
  
  plot(1:t_max,x[1,],type="o",col="black",ylab="(log(N)-mean)/SD",ylim=c(-3,3),main=paste("Simulation",nrep))
  lines(1:t_max,x[2,],type="o",col="red")
  ccf(x[1,],x[2,],ylab = "cross-correlation")
}
dev.off()

###

### Simulated phase-planes
Sigma = Diagonal(2, x = c(0.67,0.66)) ## Less variability on 
pdf("Simulated_PhasePlane_stdized.pdf",height=28,width=14)
par(mfrow=c(4,2))#,cex=1.5
plot(xbis[1,],xbis[2,],type="o",col="black",ylab="Predator",xlab="Prey",main="Real data")
arrows(xbis[1,1:33],xbis[2,1:33],xbis[1,2:34],xbis[2,2:34],length = 0.2)
#lines(1:t_max,xbis[2,],type="o",col="red")
ccf(x[1,],x[2,],ylab = "cross-correlation")
for(nrep in 1:7)
{
  ### Initial values
  x[,1]=c(-0.8,-1.8)
  for (t in 1:(t_max-1))
  {
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B %*% x[,t] + epsilon
  }

  plot(x[1,],x[2,],type="o",col="black",ylab="Predator",xlab="Prey",main=paste("Simulation",nrep))
  arrows(x[1,1:(t_max-1)],x[2,1:(t_max-1)],x[1,2:t_max],x[2,2:t_max],length = 0.2)
  ccf(x[1,],x[2,],ylab = "cross-correlation")
}
dev.off()
### Need for a correlated matrix?

{
sink("MAR_GyrPtar_Null.txt") ## Null model
cat("
    model 
{
    ### Definition of priors
    b11 ~ dunif(-1,1) # regulation prey (add possibility for overcompensation here)
    #b21 ~ dunif(0,1) # effect prey on predator
    b22 ~ dunif(0,1) # regulation predator (no overcompensation assumed)
    #b12 ~ dunif(-1,0) # effect predator on prey. Assumed negative. 
    
    ### Priors for the matrix noise parameters - multiple possible choices here. 
    #tau[1] ~ dgamma(.01,.01)
    #tau[2] ~ dgamma(.01,.01)
    for (i in 1:2){
    sigma[i] ~ dunif(0,1)
    tau[i]<-pow(sigma[i],-2)
    }
    
    ### B matrix definition
    B[1,1]<-b11
    B[1,2]<-0
    B[2,1]<-0
    B[2,2]<-b22
    
    ### Priors for initial values of x. No latent states -> no need to initialize them?
    x10 ~ dnorm(0,.01) #Those are on a log-scale, remember...
    x20 ~ dnorm(0,.01)
    
    x[1,1] ~ dnorm(x10,tau[1]) ## put the same precision as for the growth rate for the initial log-abundance
    x[2,1] ~ dnorm(x20,tau[2])
    
    #for (i in 1:2){r0[i] ~ dnorm(1,0.01)} ## prior on intercept of growth rate
    
    ### Loop over the years
    
    for (t in 1:(t_max-1)){
    for (i in 1:2){
    x[i,t+1] ~ dnorm(mu[i,t],tau[i]) ### 
    mu[i,t]<- B[i,1]*x[1,t]+B[i,2]*x[2,t] ### Line-column product
    x.rep[i,t+1] ~ dnorm(mu[i,t],tau[i])
    #xsim[i,t+1] <- x.rep[i,t+1] #can't store it for some reason
    }
    }
    # Compute Bayesian p-value
    for (i in 1:2){
    d_rep[i]<-sum(pow(x.rep[i,2:t_max]-mu[i,],2))/abs(mean(mu[i,]))
    d_data[i]<-sum(pow(mu[i,]-x[i,2:t_max],2))/abs(mean(mu[i,]))
    d_error[i]<-sum(pow(mu[i,]-x[i,2:t_max],2))/pow(sigma[i],2) 
    ## Latter distance from Ives et al. 2003 Likelihood for indep. error cov matrix - whether this can sensibly made into a R2 remains to be demonstrated. 
    rss[i]<-sum(pow(mu[i,]-x[i,2:t_max],2)) #Should be directly related to sigma through averaging...
    total_var[i]<-sum(pow(mean(x[i,2:t_max])-x[i,2:t_max],2))
    r2[i]<-rss[i]/total_var[i]
    # percentage of variance explained per species, still informative
    }
    D_rep<-sum(d_rep[])
    D_data<-sum(d_data[])
    
}
    
    ",fill=TRUE)
sink()
}

inits <- function(){list(b11=runif(1,0.5,1),b22=runif(1,0.5,1),sigma=runif(2,0,1),
                         x10=rnorm(1,0,1),x20=rnorm(1,0,1))}

parameters <- c("b11","b22","sigma","x10","x20","D_rep","D_data","rss","r2","mu")#"xsim" does not work
# above line to get mu
#parameters <- c("b11","b22","b21","b12","sigma","x10","x20","x40","D_rep","D_data","rss","r2","r0")

# MCMC settings
nc <- 3 #number of chains
nb <- 10000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-20000
nt <- 10 # “thinning”

# run model
out3 <- jags(jags.data, inits, parameters, "MAR_GyrPtar_Null.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())

print(out3, dig = 2)
## DIC = 147.0
## Previously DIC = 142.6, which was better. 

#### estimated sd(residuals) 
# sigma[1]     0.71     0.09     0.55     0.64     0.70     0.76     0.90    1  2800
# sigma[2]     0.69     0.09     0.54     0.63     0.68     0.74     0.88    1  1900
#### Both are above the previous values, not much above because CIs overlap, but still above. 
#### Likely that predator granger cause prey and prey granger cause predator. 

### Check predicted dynamics
mu=out3$BUGSoutput$mean$mu[,]

#Better plots including predictions one step ahead
par(mfrow=c(2,1))
plot(DGP$Year,xbis[1,],type="o",pch = 16,bg="black",xlab="Year",ylab="Ptarmigan density")
lines(DGP$Year[2:34],mu[1,],type="p")
arrows(DGP$Year[1:33],xbis[1,1:33],DGP$Year[2:34],mu[1,],length = 0.05)
plot(DGP$Year,xbis[2,],type="o",col="red",pch = 16,bg = "red",xlab="Year",ylab="Gyrfalcon log occupancy")
lines(DGP$Year[2:34],mu[2,],type="p",col="red")
arrows(DGP$Year[1:33],xbis[2,1:33],DGP$Year[2:34],mu[2,],col="red",length = 0.05)
### Seems less well. 

jags.sum<-out3$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_GyrPtarData_stdized_nullModel.txt")

## Bayesian P-value
par(mfrow=c(1,1))
mean(out3$BUGSoutput$sims.list$D_rep>out3$BUGSoutput$sims.list$D_data)
plot(out3$BUGSoutput$sims.list$D_rep,out3$BUGSoutput$sims.list$D_data)
abline(0,1)
#################################################################################################

#################################################################################################
######### Make the correlated error model Bayesian as well - unfinished code for info. 
#################################################################################################

# sink("MAR_GyrPtar_matrix.txt")
# cat("
# data{
#   PresMatrix[1,1]<-0.1
#   PresMatrix[1,2]<-1
#   PresMatrix[2,1]<-1
#   PresMatrix[2,2]<-0.1
#     }
#     model 
# {
#     ### Definition of priors
#     b11 ~ dunif(-1,1) # regulation prey (add possibility for overcompensation here)
#     b21 ~ dunif(0,1) # effect prey on predator
#     b22 ~ dunif(0,1) # regulation predator (no overcompensation assumed)
#     b12 ~ dunif(-1,0) # effect predator on prey. Assumed negative. 
#     
#     ### Priors for the matrix noise parameters - multiple possible choices here. 
#    
#    tau ~ dwish(PresMatrix,3)
#    Msigma<-inverse(tau)
#     
#     ### B matrix definition
#     B[1,1]<-b11
#     B[1,2]<-b12
#     B[2,1]<-b21
#     B[2,2]<-b22
#     
#     ### Priors for initial values of x. No latent states -> no need to initialize them?
#     x0[1] ~ dnorm(0,.01) #Those are on a log-scale, remember...
#     x0[2] ~ dnorm(0,.01)
#     
#     x[1:2,1] ~ dmnorm(x0[1:2],tau[,])
#     
#     ### Loop over the years
#     
#     for (t in 1:(t_max-1)){
#     
#     x[1:2,t+1] ~ dmnorm(mu[1:2,t],tau[,]) ### tau[,]?
#    
#     mu[1:2,t]<- B %*% x[1:2,t] ### Line-column product, use B[,]?
#     x.rep[1:2,t+1] ~ dmnorm(mu[1:2,t],tau[,]) ### tau[,]?
# 
#     }
# 
# }
#     
#     ",fill=TRUE)
# sink()


### Another try inverting dimensions to see what's wrong here
# 
# sink("MAR_GyrPtar_matrix.txt")
# cat("
#     model 
# {
#     ### Definition of priors
#     b11 ~ dunif(-1,1) # regulation prey (add possibility for overcompensation here)
#     b21 ~ dunif(0,1) # effect prey on predator
#     b22 ~ dunif(0,1) # regulation predator (no overcompensation assumed)
#     b12 ~ dunif(-1,0) # effect predator on prey. Assumed negative. 
#     
#     ### Priors for the matrix noise parameters - multiple possible choices here. 
#     for (i in 1:2){
#     for (j in 1:2){
#     PresMatrix[i,j]<-0.01
#     }
#     }
#     tau ~ dwish(PresMatrix,2)
#     Msigma<-inverse(tau)
#     
#     ### B matrix definition
#     B[1,1]<-b11
#     B[1,2]<-b12
#     B[2,1]<-b21
#     B[2,2]<-b22
#     
#     ### Priors for initial values of x. No latent states -> no need to initialize them?
#     x0[1] ~ dnorm(0,.01) #Those are on a log-scale, remember...
#     x0[2] ~ dnorm(0,.01)
#     
#     x[1,1:2] ~ dmnorm(x0[1:2],tau[,])
#     
#     ### Loop over the years
#     
#     for (t in 1:(t_max-1)){
#     
#     x[t+1,1:2] ~ dmnorm(mu[t,1:2],tau)
#     mu[t,1:2]<- x[t,1:2] %*% t(B) ###  use B[,]? Line-column product before, need to refine dimensions again?
#     x.rep[t+1,1:2] ~ dmnorm(mu[t,1:2],tau) ### tau[,]?
#     
#     }
#     
# }
#     
#     ",fill=TRUE)
# sink()
# 
# jags.data<-list(t_max=ncol(xbis),x=t(xbis))
# 
# jags.data<-list(t_max=ncol(xbis),x=xbis)
# 
# inits <- function(){list(b11=runif(1,0.5,1),b22=runif(1,0.5,1),b12=runif(1,-1,0),b21=runif(1,0,0.5),tau=matrix(runif(4,0,0.1),2,2),x0=rnorm(2,0,1))}
# 
# #parameters <- c("b11","b22","b21","b12","Msigma","x10","x20","mu")#"xsim" does not work
# # above line to get mu
# parameters <- c("b11","b22","b21","b12","Msigma","x10","x20")
# 
# # MCMC settings
# nc <- 3 #number of chains
# nb <- 10000 # “burn in”
# #ni <- 14000# “number of iterations” # that's for a symmetric distrib...
# ni<-20000
# nt <- 10 # “thinning”
# 
# # run model
# out <- jags(jags.data, inits, parameters, "MAR_GyrPtar_matrix.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
# 
# print(out, dig = 2)

### Again invalid parent values, whatever I put as symmetric matrix for the precision values
# Error in jags.model(model.file, data = data, inits = init.values, n.chains = n.chains,  : 
#                       Error in node x[1:2,1]
#                     Invalid parent values
### I need to investigate by doing smaller models with the mvnorm distribution. 
