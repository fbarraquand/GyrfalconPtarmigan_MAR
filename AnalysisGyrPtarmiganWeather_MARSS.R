################################################################################################################
################# Analysis Gyrfalcon-Partmigan-Weather data - MAR modelling with MARSS package #################
### FBarraquand 26/03/2017, analyses started 06/07/2015, with O. Nielsen #######################################
### Updated version producing all the figures and results of the paper in a reproducible workflow 25/07/2017 ###
### Cleaned up commments and removed some redundant pieces of code 29/08/2018 ##################################
### FB 04/03/2020 Adding more bottom-up models and correcting previous error  ##################################
################################################################################################################

### Initializing
rm(list=ls())
graphics.off()

### Setting the seed to keep the same simulations
set.seed(42) # What else?

options(digits=4)

library(MASS)

### Reading data on gyr-ptarmigan
DGP<-read.csv("Gyrfalcon_Data.csv")
par(mfrow=c(2,1))
plot(DGP$Year,DGP$Occupied,type="b")
plot(DGP$Year,DGP$Ptarmigan,type="b")
## Need to correct the gyr portion for observational effort // number of territories surveyed 
DGP$Occupancy=DGP$Occupied/DGP$N

par(mfrow=c(2,1))
plot(DGP$Year,DGP$Occupancy,type="b",main="Percentage territories occupied Gyrfalcon")
plot(DGP$Year,DGP$Ptarmigan,type="b",main="Mean density Ptarmigan")

#Plot both densities standardized
png(file="GyrPtarDensities.png",width=10,height=8,res=300,units="in")
DGP$OccStd=(DGP$Occupancy-mean(DGP$Occupancy))/sd(DGP$Occupancy)
DGP$PtarStd=(DGP$Ptarmigan-mean(DGP$Ptarmigan))/sd(DGP$Ptarmigan)
par(mfrow=c(1,1),cex=1.5)
plot(DGP$Year,DGP$OccStd,type="b",col="red",ylim=c(-3,3),ylab = "Stdized population density",xlab="Year",lwd=3,pch=20)
lines(DGP$Year,DGP$PtarStd,type="b",lwd=3,pch=20)
dev.off()

#####################################################################################################
### Easier to interpret the coefficients once the data is both logged and standardized. Do that now. 
################# Standardized data analysis ########################################################

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
x=xbis[,1:(nrow(DGP)-1)] ### restricted dataset (not for MARSS, for custom analyses)

################# Weather variables #################################################################

### Loading in the weather data 
DB=read.csv("weather_iceland/average_weatherNEiceland.csv",header=T)
names(DB)[6]="r" #renaming "r" the logRainfall variable
### Years considered for weather
DGP$Year[1:(nrow(DGP))] # 1981 to 2014, - last growth rate is 2013-2014

##################### Weather data for prey ########################

### Defining a lagged year
year_minus_1=DGP$Year[1:(nrow(DGP))]-1

# Average temperature 1 year earlier in May 
tempMay_year_minus1=DB$temp[(DB$year %in% year_minus_1)&(DB$month==5)]
### This is in fact weather of the spring t-1 affecting growth from spring t-1 to spring t 

# Average temperature of the year in May 
# This makes sense because chicks of the previous year survival depend on the year's temperature in the spring. 
# Using the previous year encapsulate effects of temperature on reproduction and early chick survival
DGP$Year
# This is weather of the spring t affecting growth from spring t-1 to spring t (later on in the MAR model)
tempMay_year=DB$temp[(DB$year %in% DGP$Year)&(DB$month==5)]

# Log-rainfall of the year in May
rainMay_year=log(DB$r[(DB$year %in% DGP$Year)&(DB$month==5)])
# Log-rainfall 1 year earlier in May
rainMay_year_minus1=log(DB$r[(DB$year %in% year_minus_1)&(DB$month==5)])

################## Weather data for predator ###################################
year_minus_4=DGP$Year[1:(nrow(DGP))]-5 ### if we want the weather four years back affecting gyr growth, 
## this is what we have to write because MARSS includes the effects of covariate at the same time as "end" abundance
## Hence the growth of the predator between 2013 and 2014 depends on weather in 2009 (2009 + 4 = 2013)
### NB: if we don't correct this, there's no effect and this is also something we might want to discuss.
## we don't do this correction for the grouse because the life-history differs

# Average temperature 4 years earlier in April
tempApril_year_minus4=DB$temp[(DB$year %in% year_minus_4)&(DB$month==4)]
# Log-rainfall 4 years earlier in April
rainApril_year_minus4=log(DB$r[(DB$year %in% year_minus_4)&(DB$month==4)])

### --- Standardize all those variables to be able to compare something --- ###
tempMay_year_minus1=(tempMay_year_minus1-mean(tempMay_year_minus1))/sd(tempMay_year_minus1)
tempMay_year = (tempMay_year-mean(tempMay_year))/sd(tempMay_year)
tempApril_year_minus4=(tempApril_year_minus4-mean(tempApril_year_minus4))/sd(tempApril_year_minus4)

rainMay_year_minus1=scale(rainMay_year_minus1,scale=TRUE) #tired of writing
rainMay_year=scale(rainMay_year,scale=TRUE) #tired of writing
rainApril_year_minus4=scale(rainApril_year_minus4,scale=TRUE)

################### Now implement the suggestion of Olafur - looking at July's data #######################################
### If chick survival to adulthood is a key factor, it has to be July from the year t-1
### ----  Note the ptarmigans are counted in spring, that is, in May and early June. 
### Makes therefore no sense to look at how ptarmigan growth from May_{t-1} to May_t is affected by June_t or July_t

### Creating new weather variables
tempJuly_year=DB$temp[(DB$year %in% year_minus_1)&(DB$month==7)]
# Log-rainfall in July
rainJuly_year=log(DB$r[(DB$year %in% year_minus_1)&(DB$month==7)])
### June weather
tempJune_year=DB$temp[(DB$year %in% year_minus_1)&(DB$month==6)]
# Log-rainfall in June
rainJune_year=log(DB$r[(DB$year %in% year_minus_1)&(DB$month==6)])
### Scale these
tempJune_year=scale(tempJune_year)
tempJuly_year=scale(tempJuly_year)

## rainfall
rainJune_year=scale(rainJune_year)
rainJuly_year=scale(rainJuly_year)


### Winter weather now ################################

### New idea (cf. discussion Olafur 21/03/2017), to look at winter weather (we talked about this before, then dismissed it)
### Oli pointed out that the age ratios might hint at low winter survival for the chicks
### Need some kind of proxy of winter severity. 

year_now=DGP$Year[1:(nrow(DGP))]
### Now the weather in the year t might affect growth from t-1 to t if we consider the end of winter.  

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


################# Estimation of MAR models #############################################

### Estimation with MARSS package
library('MARSS')
### See e.g http://cran.r-project.org/web/packages/MARSS/vignettes/Quick_Start.pdf


### State-space, observation part - never changes 
Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
### Initial values
pi1=matrix(0,2,1); #Inial values
V1=diag(1,2)

### Process model part
### Setting matrices
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
#Q1=matrix(c("q11",0,0,"q22"),2,2) ### Assumes no correlated noise (- or does it?)
Q1="diagonal and unequal"

# Estimation
data<-xbis
model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
mar1.full=MARSS(data, model=model.list)
CIs.mar1.full=MARSSparamCIs(mar1.full)

### Produce Fig 1 -- time series with model predictions 
B=matrix(mar1.full$par$B,nrow=2)#value=CIs.mar1.full$par$B
Q=diag(as.vector(mar1.full$par$Q))
t_max=34
xsimrepeats=array(data=0,dim=c(2,100,t_max-1))
mu=c(0,0)
for (t in 1:(t_max-1))
{
  for(nrep in 1:100)
  {
    x=xbis
    epsilon=mvrnorm(n = 1, mu, Q)
    xnew = B %*% xbis[,t] + epsilon
    xsimrepeats[,nrep,t] = xnew
  }
}

### Plotting the whole thing
pdf(file = "PredictedLogAbundances.pdf",width=10,height=6)
par(mfrow=c(1,1),cex=1.5)
plot(DGP$Year,xbis[1,],type="b",col="black",ylim=c(-3,3),ylab = "Stdized log(population density)",xlab="Year",lwd=3,pch=20)
lines(DGP$Year,xbis[2,],type="b",col="red",lwd=3,pch=20)
for (k in 1:50){
  points(DGP$Year[2:t_max],xsimrepeats[1,k,],col="black",pch=".")
  points(DGP$Year[2:t_max],xsimrepeats[2,k,],col="red",pch=".")
}
dev.off()
############## End of plotting for Fig. 1 ####################

### Now include the same full model but with a correlated noise matrix
Q1=matrix(c("q11","q21","q12","q22"),2,2) ##assume correlated noise
model.list=list(B=B1,U=U1,Q="unconstrained",Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
mar1.full2=MARSS(data, model=model.list)
MARSSparamCIs(mar1.full2)

### Now we assume a non-correlated matrix is OK. 
Q1="diagonal and unequal"
### Now include a null model, diagonal MAR(1) model without interactions
B1=matrix(list("b11",0,0,"b22"),2,2,byrow = T) ### Interaction matrix
U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
mar1.null=MARSS(data, model=model.list)
CIs.mar1.null=MARSSparamCIs(mar1.null)

### Fit quality
mar1.null$logLik
mar1.null$AIC
mar1.null$AICc

### Now include a model with temperature covariates
covar=t(as.matrix(cbind(tempMay_year,tempApril_year_minus4)))
## C1=matrix(c("tempMay_year",0,0,"tempApril_year_minus4"),2,2,byrow=T) ## creates problems
C1=matrix(list("tempMay_year",0,0,"tempApril_year_minus4"),2,2,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.temp=MARSS(data, model=model.list)
CIs.mar1.temp=MARSSparamCIs(mar1.temp)

####### Storage of interaction matrices B and environmental matrices C ################

## function for storage 
mar_results <- function(model_name,CIs_model_name){
  
  value=model_name$par$B#value=CIs_model_name$par$B
  SE=CIs_model_name$par.se$B
  lower=CIs_model_name$par.lowCI$B
  upper=CIs_model_name$par.upCI$B
  mar1.results=data.frame(value,SE,lower,upper)
  
  value=model_name$par$U ### C is U
  SE=CIs_model_name$par.se$U
  lower=CIs_model_name$par.lowCI$U
  upper=CIs_model_name$par.upCI$U
  mar1.results=rbind(mar1.results,data.frame(value,SE,lower,upper))
  
  value=model_name$par$Q ### C is U
  SE=CIs_model_name$par.se$Q
  lower=CIs_model_name$par.lowCI$Q
  upper=CIs_model_name$par.upCI$Q
  mar1.results=rbind(mar1.results,data.frame(value,SE,lower,upper))
  return(mar1.results)
} 

# Storage full model
mar1.data=mar_results(mar1.full,CIs.mar1.full)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.full.csv")

mar1.null$par$B
#Storage null model
mar1.data=mar_results(mar1.null,CIs.mar1.null)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.null.csv")

mar1.temp$par$B
#Storage full model with effect of temperature
mar1.data = mar_results(mar1.temp,CIs.mar1.temp)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.temp.csv")

covar=t(as.matrix(cbind(tempMay_year,tempApril_year_minus4)))
C1=matrix(list("tempMay_year",0,0,"tempApril_year_minus4"),2,2,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11",0,0,"b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.temp.only=MARSS(data, model=model.list)
CIs.mar1.temp.only=MARSSparamCIs(mar1.temp.only)

mar1.data=mar_results(mar1.temp.only,CIs.mar1.temp.only)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.temp.only.csv")

### Delay one more year the temperature for ptarmigan. 
covar=t(as.matrix(cbind(tempMay_year_minus1,tempApril_year_minus4)))
C1=matrix(list("tempMay_year",0,0,"tempApril_year_minus4"),2,2,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.temp1=MARSS(data, model=model.list)
CIs.mar1.temp1=MARSSparamCIs(mar1.temp1)

#Storage full model with effect of temperature one year delayed
mar1.data=mar_results(mar1.temp1,CIs.mar1.temp1)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.temp1.csv")

### Now a model with rainfall
covar=t(as.matrix(cbind(rainMay_year,rainApril_year_minus4)))
C1=matrix(list("rainMay_year",0,0,"rainApril_year_minus4"),2,2,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.rain=MARSS(data, model=model.list)
CIs.mar1.rain = MARSSparamCIs(mar1.rain)
# No effect. 

#Storage full model with effect of (cumulated) rainfall
mar1.data=mar_results(mar1.rain,CIs.mar1.rain)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.rain.csv")

### With temperature + rainfall
covar=t(as.matrix(cbind(tempMay_year,rainMay_year,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("tempMay_year","rainMay_year",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.both=MARSS(data, model=model.list)
CIs.mar1.both=MARSSparamCIs(mar1.both)
mar1.both

mar1.data=mar_results(mar1.both,CIs.mar1.both)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.both.csv")

### First comparison of AICc/BIC for all models 
mar.bic <- function(my.mar)
{
  my.mar$BIC=0  
  my.mar$BIC=-2*my.mar$logLik + my.mar$num.params*log(my.mar$samp.size/2)
}

## Initialize
mar1.null$BIC=mar.bic(mar1.null)
aic.table=data.frame(mar1.null$logLik,mar1.null$AIC,mar1.null$AICc,mar1.null$BIC)
names(aic.table)=c("logLik","AIC","AICc","BIC")

mar1.full$BIC=mar.bic(mar1.full)
aic.table.temp=data.frame(mar1.full$logLik,mar1.full$AIC,mar1.full$AICc,mar1.full$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

mar1.temp$BIC=mar.bic(mar1.temp)
aic.table.temp=data.frame(mar1.temp$logLik,mar1.temp$AIC,mar1.temp$AICc,mar1.temp$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

mar1.temp.only$BIC=mar.bic(mar1.temp.only)
aic.table.temp=data.frame(mar1.temp.only$logLik,mar1.temp.only$AIC,mar1.temp.only$AICc,mar1.temp.only$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

mar1.temp1$BIC=mar.bic(mar1.temp1)
aic.table.temp=data.frame(mar1.temp1$logLik,mar1.temp1$AIC,mar1.temp1$AICc,mar1.temp1$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

mar1.rain$BIC=mar.bic(mar1.rain)
aic.table.temp=data.frame(mar1.rain$logLik,mar1.rain$AIC,mar1.rain$AICc,mar1.rain$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

mar1.both$BIC=mar.bic(mar1.both)
aic.table.temp=data.frame(mar1.both$logLik,mar1.both$AIC,mar1.both$AICc,mar1.both$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

rownames(aic.table)=c("mar1.null","mar1.full","mar1.temp","mar1.temp.only","mar1.temp1","mar1.rain","mar1.both")

### AIC table shows the full model is better -- did we find that earlier? Not by much though
aic.table
write.csv(format(aic.table,digits=4),file="mar1/aic.table.csv")
### we keep in the set of nested models mar1.full, mar1.temp and mar1.temp1

#############################################################
### Alternate models with June and July weather to check
##############################################################
covar=t(as.matrix(cbind(tempJune_year,rainJune_year,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("tempJune_year","rainJune_year",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.both.june=MARSS(data, model=model.list)
CIs.mar1.both.june=MARSSparamCIs(mar1.both.june)

### Store model with June temperature
mar1.data=mar_results(mar1.both.june,CIs.mar1.both.june)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.both.june.csv")

covar=t(as.matrix(cbind(tempJuly_year,rainJuly_year,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("tempJuly_year","rainJuly_year",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.both.july=MARSS(data, model=model.list)
CIs.mar1.both.july=MARSSparamCIs(mar1.both.july)

mar1.data=mar_results(mar1.both.july,CIs.mar1.both.july)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.both.july.csv")

nrow(aic.table) ## How many models have we now

mar1.both.july$BIC=mar.bic(mar1.both.july)
aic.table.temp=data.frame(mar1.both.july$logLik,mar1.both.july$AIC,mar1.both.july$AICc,mar1.both.july$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

mar1.both.june$BIC=mar.bic(mar1.both.june)
aic.table.temp=data.frame(mar1.both.june$logLik,mar1.both.june$AIC,mar1.both.june$AICc,mar1.both.june$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

aic.table[8:9,]
rownames(aic.table)[8:9]=c("mar1.both.july","mar1.both.june")
write.csv(format(aic.table,digits=4),file="mar1/aic.table.csv")

#############################################################
### Alternate models with winter weather to check
##############################################################

### Variables are (the first two ones are correlated)
#avMonthly_tempWinter
#minOfMonths_tempWinter
#log_winter_rain

covar=t(as.matrix(cbind(avMonthly_tempWinter,log_winter_rain,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("avMonthly_tempWinter","log_winter_rain",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
### Stronger controls to check the warnings
iter_min=100 # a minimum of 100 iterations
cntl.list=list(conv.test.slope.tol=0.001,minit=iter_min,maxit=500,abstol=0.001,trace=1,safe=T)
#,MCInit=TRUE,silent=2)#,numInits=iter_estimate,numInitSteps=10) 
#Took values from Griffiths for minit,maxit,abstol. 
# conv.test.slope is recommended in MARSS User Guide itself, 
# numInits and numInitStep are taken as rules of thumbs 
#(quick search on the Internet, numInits=500 elsewhere but it seems really big to me
mar1.both.winter=MARSS(data, model=model.list,control=cntl.list)
MARSSparamCIs(mar1.both.winter)
mar1.both.winter$AICc
CIs.mar1.both.winter=MARSSparamCIs(mar1.both.winter)

covar=t(as.matrix(cbind(avMonthly_tempWinter,log_winter_rain,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("minOfMonths_tempWinter","log_winter_rain",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.both.winter2=MARSS(data, model=model.list)
MARSSparamCIs(mar1.both.winter2)
CIs.mar1.both.winter2=MARSSparamCIs(mar1.both.winter2)

#### Check again those models
mar1.data=mar_results(mar1.both.winter,CIs.mar1.both.winter)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.both.winter.csv")

### Same for model 2
mar1.data=mar_results(mar1.both.winter2,CIs.mar1.both.winter2)
write.csv(format(mar1.data,digits=4),file="mar1/mar1.both.winter2.csv")

### AIC Table winter

mar1.both.winter$BIC=mar.bic(mar1.both.winter)
aic.table=data.frame(mar1.both.winter$logLik,mar1.both.winter$AIC,mar1.both.winter$AICc,mar1.both.winter$BIC)
names(aic.table)=c("logLik","AIC","AICc","BIC")

mar1.both.winter2$BIC=mar.bic(mar1.both.winter2)
aic.table.temp=data.frame(mar1.both.winter2$logLik,mar1.both.winter2$AIC,mar1.both.winter2$AICc,mar1.both.winter2$BIC)
names(aic.table.temp)=c("logLik","AIC","AICc","BIC")
aic.table=rbind(aic.table,aic.table.temp)

rownames(aic.table)=c("mar1.both.winter","mar1.both.winter2")
aic.table
write.csv(format(aic.table,digits=4),file="mar1/aic.table.winter.csv")

#################################################################################################################
# Now considering MAR(2) models -- what if prey has an impact on the predator but not the other way around ? 
# What if both population are in fact independent (unlikely based on simulations but we'll try that)
##################################################################################################################

### Format generally the model
### See Chap 16 here - https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
#### Use MARSS example with two AR(2) models ########
TT=50
true.2=c(r=0,b1=-1.5,b2=-0.75,q=1)
temp1=arima.sim(n=TT,list(ar=true.2[c("b1","b2")]),sd=sqrt(true.2["q"]))
temp2=arima.sim(n=TT,list(ar=true.2[c("b1","b2")]),sd=sqrt(true.2["q"]))
sim.mar2=rbind(temp1,temp2)

Z=matrix(c(1,0,0,1,0,0,0,0),2,4)
B1=matrix(list(0),2,2); diag(B1)="b1"
B2=matrix(list(0),2,2); diag(B2)="b2"
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
U=matrix(0,4,1)
Q=matrix(list(0),4,4)
Q[1,1]="q"; Q[2,2]="q"
A=matrix(0,2,1)
R=matrix(0,2,2)
pi=matrix(c(sim.mar2[,2],sim.mar2[,1]),4,1)
V=matrix(0,4,4)
model.list.2m=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2=MARSS(sim.mar2[,2:TT],model=model.list.2m)

##################### Now on our data #############################################

### State-space, observation part - never changes 
Z=matrix(c(1,0,0,1,0,0,0,0),2,4) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
A=matrix(0,2,1)### Intercept state space = 0
R=matrix(0,2,2) ### Error matrix state-space = 0

#Z=matrix(c(1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),4,4,byrow=T) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
#A=matrix(0,4,1)### Intercept state space = 0
#R=matrix(0,4,4) ### Error matrix state-space = 0

### Initial values
V=matrix(0,4,4)
pi=matrix(c(xbis[,2],xbis[,1]),4,1)
pi

### Interaction matrix
B1=matrix(list("b11_1","b12_1","b21_1","b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2","b21_2","b22_2"),2,2,byrow = T)
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

### Process error matrix
U=matrix(0,4,1)
Q=matrix(list(0),4,4)
Q[1,1]="q11"; Q[2,2]="q22"
Q
#stacked_data=rbind(xbis[,2:nrow(DGP)],xbis[,1:(nrow(DGP)-1)])

### Model call
model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
### 
mar2.full=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
CIs.mar2.full=MARSSparamCIs(mar2.full)
# we don't use stacked_data here, just the abundance data that is "observed at t"

### Store data
mar1.data=mar_results(mar2.full,CIs.mar2.full)
write.csv(format(mar1.data,digits=4),file="mar2/mar2.full.csv")

#### Bottom-up model that we highlighted before -- with an effect of prey on predator growth the next year. 
### Interaction matrix
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T) 
B2=matrix(list("b11_2",0,"b21_2","b22_2"),2,2,byrow = T) ## error was here, now corrected
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.bottom.up=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
CIs.mar2.bottom.up=MARSSparamCIs(mar2.bottom.up)

### Store data
mar1.data=mar_results(mar2.bottom.up,CIs.mar2.bottom.up)
write.csv(format(mar1.data,digits=4),file="mar2/mar2.bottom.up.csv")


#### Bottom-up model -- more direct effect of prey on predator growth 
### Interaction matrix
B1=matrix(list("b11_1",0,"b21_1","b22_1"),2,2,byrow = T) 
B2=matrix(list("b11_2",0,"b21_2","b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.BUv1=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
CIs.mar2.BUv1=MARSSparamCIs(mar2.BUv1)
### coefficients make moderate sense (B.b21_1 negative)

### Store data
mar1.data=mar_results(mar2.BUv1,CIs.mar2.BUv1)
write.csv(format(mar1.data,digits=4),file="mar2/mar2.BUv1.csv")

### Effect of prey on predator growth of that year and no delayed predator dd
B1=matrix(list("b11_1",0,"b21_1","b22_1"),2,2,byrow = T) 
B2=matrix(list("b11_2",0,0,0),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.BUv2=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
CIs.mar2.BUv2=MARSSparamCIs(mar2.BUv2)

### Store data
mar1.data=mar_results(mar2.BUv2,CIs.mar2.BUv2)
write.csv(format(mar1.data,digits=4),file="mar2/mar2.BUv2.csv")

### Effect of prey on predator growth of that year and delayed predator dd
B1=matrix(list("b11_1",0,"b21_1","b22_1"),2,2,byrow = T) 
B2=matrix(list("b11_2",0,0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.BUv3=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
CIs.mar2.BUv3=MARSSparamCIs(mar2.BUv3)

### Store data
mar1.data=mar_results(mar2.BUv3,CIs.mar2.BUv3)
write.csv(format(mar1.data,digits=4),file="mar2/mar2.BUv3.csv")

### Independent AR(2) models
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2",0,0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.indep=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
MARSSparamCIs(mar2.indep)
CIs.mar2.indep=MARSSparamCIs(mar2.indep)

### Store data
mar1.data=mar_results(mar2.indep,CIs.mar2.indep)
write.csv(format(mar1.data,digits=4),file="mar2/mar2.indep.csv")

################### Full MAR(1) with same number of data points // TS length affects AIC ##################
### State-space, observation part - never changes 
Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
### Initial values
pi1=matrix(0,2,1); #Inial values
V1=diag(1,2)

### Process model part
### Setting matrices
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
Q1="diagonal and unequal"

# Estimation
data<-xbis
model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
mar1.full.bis=MARSS(xbis[,2:nrow(DGP)], model=model.list)
MARSSparamCIs(mar1.full.bis) ### AIC 140, AICc 142, less of a difference... 
### But still 10 AIC points. I need to check what the BIC might be saying. 

### Keep on going there. 
### Now include a null model, diagonal MAR(1) model without interactions
B1=matrix(list("b11",0,0,"b22"),2,2,byrow = T) ### Interaction matrix
U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
mar1.null.bis=MARSS(xbis[,2:nrow(DGP)], model=model.list)
MARSSparamCIs(mar1.null.bis)


#########################################################################################
######################### MAR(2) models with weather covariates #########################
#########################################################################################
### Let's add some weather on those 
covar=t(as.matrix(cbind(tempMay_year[2:34],tempApril_year_minus4[2:34])))
C1=matrix(list("tempMay_year",0,0,"tempApril_year_minus4",0,0,0,0),4,2,byrow=T)
C1
### Initial values
V=matrix(0,4,4)
pi=matrix(c(xbis[,2],xbis[,1]),4,1)
pi
model.list.2lags=list(Z=Z,B=B,U=U,C=C1,c=covar,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.indep.temp=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags) #,MCInit=TRUE

######## Try some other model fitting
mar2.indep.temp=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags,method="BFGS") 

MARSSparamCIs(mar2.indep.temp)
# MARSSaic(mar2.indep.temp, output = "AICbp")
# AICbp calculation in progress...
# |2%      |20%      |40%      |60%      |80%      |100%
#   Progress: ||||||||||||||||||||||||||||||||||||||||||||||||||
#   
#   MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Estimation converged in 17 iterations. 
# Log-likelihood: -56.36994 
# AIC: 130.7399   AICc: 133.9542   AICbp(param): 139.654   
# 
# Estimate
# B.b11_1                   1.0833
# B.b22_1                   0.6966
# B.b11_2                  -0.4530
# B.b22_2                  -0.0709
# Q.q11                     0.3299
# Q.q22                     0.3522
# C.tempMay_year           -0.0371
# C.0                       0.1061
# C.tempApril_year_minus4   0.1937
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.
# 
# MARSSkem warnings. Type MARSSinfo() for help.
# MARSSkem: The soln became unstable and logLik DROPPED.
#Use control$trace=1 to generate a more detailed error report.

##################################################################################################
###### Compare fit of models
##################################################################################################
mar2.full$AICc
mar2.bottom.up$AICc
mar2.BUv1$AICc
mar2.BUv2$AICc
mar2.BUv3$AICc
mar2.indep$AICc 
mar1.full.bis$AICc
mar1.null.bis$AICc #even higher, how is it possible when compared to VAR? 

mar2.full$AIC
mar2.bottom.up$AIC
mar2.BUv1$AIC
mar2.BUv2$AIC
mar2.BUv3$AIC
mar2.indep$AIC 
mar1.full.bis$AIC
mar1.null.bis$AIC
### They are quite close to each other, but clearly below MAR(1)
### T=33
### BIC=AIC-2*kparam+kparam*log(T)

#### Creates table with AIC and BIC values
## Initialize
mar1.null.bis$BIC=mar.bic(mar1.null.bis)
aic.table2=data.frame(mar1.null.bis$logLik,mar1.null.bis$AIC,mar1.null.bis$AICc,mar1.null.bis$BIC)
names(aic.table2)=c("logLik","AIC","AICc","BIC")

mar1.full.bis$BIC=mar.bic(mar1.full.bis)
aic.table2.temp=data.frame(mar1.full.bis$logLik,mar1.full.bis$AIC,mar1.full.bis$AICc,mar1.full.bis$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.full$BIC=mar.bic(mar2.full)
aic.table2.temp=data.frame(mar2.full$logLik,mar2.full$AIC,mar2.full$AICc,mar2.full$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.bottom.up$BIC=mar.bic(mar2.bottom.up)
aic.table2.temp=data.frame(mar2.bottom.up$logLik,mar2.bottom.up$AIC,mar2.bottom.up$AICc,mar2.bottom.up$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.BUv1$BIC=mar.bic(mar2.BUv1)
aic.table2.temp=data.frame(mar2.BUv1$logLik,mar2.BUv1$AIC,mar2.BUv1$AICc,mar2.BUv1$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.BUv2$BIC=mar.bic(mar2.BUv2)
aic.table2.temp=data.frame(mar2.BUv2$logLik,mar2.BUv2$AIC,mar2.BUv2$AICc,mar2.BUv2$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.BUv3$BIC=mar.bic(mar2.BUv3)
aic.table2.temp=data.frame(mar2.BUv3$logLik,mar2.BUv3$AIC,mar2.BUv3$AICc,mar2.BUv3$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.indep$BIC=mar.bic(mar2.indep)
aic.table2.temp=data.frame(mar2.indep$logLik,mar2.indep$AIC,mar2.indep$AICc,mar2.indep$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.indep.temp$BIC=mar.bic(mar2.indep.temp)
aic.table2.temp=data.frame(mar2.indep.temp$logLik,mar2.indep.temp$AIC,mar2.indep.temp$AICc,mar2.indep.temp$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

rownames(aic.table2)=c("mar1.null.bis","mar1.full.bis","mar2.full","mar2.bottom.up","mar2.BUv1","mar2.BUv2","mar2.BUv3","mar2.indep","mar2.indep.temp")

aic.table2
write.csv(format(aic.table2,digits=4),file="mar2/aic.table2.csv")

############ Previous estimates using AICb for memory ################################################
#MARSSaic(mar2.full, output = "AICbp") ### AIC: 129.5646   AICc: 133.5646   AICbp(param): 143.6595 
#MARSSaic(mar2.bottom.up, output = "AICbp") ### AIC: 129.8294   AICc: 131.7604   AICbp(param): 137.406   
#MARSSaic(mar2.indep,output = "AICbp") ### AIC: 129.5638   AICc: 130.9876   AICbp(param): 133.7121   
#MARSSaic(mar1.full.bis,output = "AICbp") ### AIC: 140.8294   AICc: 142.2531   AICbp(param): 148.7717   
#MARSSaic(mar1.null.bis,output = "AICbp") ### AIC: 143.1355   AICc: 143.7913   AICbp(param): 147.0655   
#MARSSaic(mar2.indep.temp,output="AICbp") ### AIC: 129.5638   AICc: 130.9876   AICbp(param): 134.6404   
######################################################################################################

######################################################################################################
################# Using VAR to check the model order #################################################
######################################################################################################

library(vars)
varpp<-VAR(y=data.frame(t(xbis)), type="none",lag.max=5)
## Considering 5 maximum lags and using VAR(p) estimation with "vars" package (MAR(p) in ecology)
varpp #Yields model order=3
### Note this is somewhat consistent with previous results, at least concerning AIC (see below for BIC)
## http://raunvisindastofnun.hi.is/sites/raunvisindastofnun.hi.is/files/rh-18-2003.pdf
# I also tried a lag.max = 20, gives lag 7 -- clearly this is overparameterized (lag.max = 20 is nonsensical here)

###### Looking at several model selection criteria
var_order_select=VARselect(y=data.frame(t(xbis)), type="none",lag.max=5)
var_order_select

### The AIC selects 3 lags and the BIC and HQ, that are generally more conservative, two lags. 

### Do these models see causality in the 2x2 model?

######## MAR(2) model #################
varpp2<-VAR(y=data.frame(t(xbis)), p=2, type="none")
causality(varpp2,cause="X2") ## No effect of X2 on X1
causality(varpp2,cause="X1") ## close to reject the hypothesis of non-GC

######## MAR(1) model ################# 
varpp1<-VAR(y=data.frame(t(xbis)), p=1, type="none") 
causality(varpp1,cause="X2") ## reject
causality(varpp1,cause="X1") ## reject
# Looks better for causality but that's not the model that's selected... 
#######################################################################################################

#####################################################################################
##################### Simulation-based model selection ##############################
####################################################################################

###### 1. Can we get the right cross-correlation patterns with the best-fitting models? 

# Output the cross-correlation pattern for the data

# Simulations of the fitted models -- only for MAR(1) null, MAR(1) full, MAR(2) full, MAR(2) bottom-up, MAR(2) indep

# -- Extract data from MAR(1) null model ---
data1=read.csv(file="mar1/mar1.null.csv")
data1
B = diag(c(data1[1,]$value,data1[2,]$value))
Sigma = diag(c(data1[3,]$value,data1[4,]$value)) ## Less variability on 
# Simulation of model and computation of cross-correlation
t_max=34
# pdf("Simulated_MAR_dynamics_stdized_MAR1_null.pdf",height=28,width=14)
# par(mfrow=c(4,2))#,cex=1.5
# plot(1:t_max,xbis[1,],type="o",col="black",ylim=c(-3,3),ylab="Real data")#ylim=c(-3,3)
# lines(1:t_max,xbis[2,],type="o",col="red")
# ccf(xbis[1,],xbis[2,],ylab = "cross-correlation")
# for(nrep in 1:7)
# {
#   ### Initial values
#   x[,1]=c(-0.8,-1.8)
#   for (t in 1:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B %*% x[,t] + epsilon
#   }
#   
#   plot(1:t_max,x[1,],type="o",col="black",ylab="(log(N)-mean)/SD",ylim=c(-3,3),main=paste("Simulation",nrep))
#   lines(1:t_max,x[2,],type="o",col="red")
#   ccf(x[1,],x[2,],ylab = "cross-correlation")
# }
# dev.off()

### -- Better representation with real cross-corr overlayed onto simulated cross-corr--- ##
# Cross-correlation for the real data 
c_hat = ccf(xbis[1,],xbis[2,],ylab = "cross-correlation",plot=FALSE)
str(c_hat)
cc_hat = c_hat$acf
n_cc = length(cc_hat)
y_cc_mar1null=matrix(0,100,n_cc)
### Use 100 values to get a simulation enveloppe 
x=matrix(0,2,t_max)
for(nrep in 1:100)
{
  ### Initial values
  x[,1]=xbis[,1]
  for (t in 1:(t_max-1))
  {
    mu=c(0,0)
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B %*% x[,t] + epsilon
  }
  
  y_cc_mar1null[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation")$acf)
}

# Start plotting --
pdf("crossCorrelations_simulatedModels.pdf",height=16,width=16)
par(mfrow=c(3,2),cex=1.4)#,cex=1.5
# --- First panel -----
matplot(c_hat$lag,t(y_cc_mar1null), pch=".",ylab="Cross-correlation",xlab="Time lag",main= "MAR(1) - no interactions")
matlines(c_hat$lag,t(y_cc_mar1null))
lines(c_hat$lag,cc_hat,col="black",lwd=4)
add_label_legend <- function(pos = "topleft", label, ...) {
  legend(pos, inset=c(-0.2,-0.4), label, bty = "n", xpd=TRUE, ...)
}
add_label_legend(label="A",cex=1.8)
#add_label_legend(c(0.9,0.7),label="A")
#add_label_legend(c(0.8,0.6),label="A")
# --- end of plot on MAR(1) null model --- # 

# --- Extract data from MAR(1) full model --- #
data.mar1.full=read.csv(file="mar1/mar1.full.csv")
data.mar1.full
B = matrix(data.mar1.full$value[1:4],2,2)
Sigma = diag(data.mar1.full[5:6,]$value) ## Less variability on 
# Simulation of model and computation of cross-correlation
t_max=34
# pdf("Simulated_MAR_dynamics_stdized_MAR1_full.pdf",height=28,width=14)
# par(mfrow=c(4,2))#,cex=1.5
# plot(1:t_max,xbis[1,],type="o",col="black",ylim=c(-3,3),ylab="Real data")#ylim=c(-3,3)
# lines(1:t_max,xbis[2,],type="o",col="red")
# ccf(xbis[1,],xbis[2,],ylab = "cross-correlation")
# for(nrep in 1:7)
# {
#   ### Initial values
#   x[,1]=c(-0.8,-1.8)
#   for (t in 1:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B %*% x[,t] + epsilon
#   }
#   
#   plot(1:t_max,x[1,],type="o",col="black",ylab="(log(N)-mean)/SD",ylim=c(-3,3),main=paste("Simulation",nrep))
#   lines(1:t_max,x[2,],type="o",col="red")
#   ccf(x[1,],x[2,],ylab = "cross-correlation")
# }
# dev.off()

### -- Better representation with cross-corr overlayed --- ##
# Cross-correlation for the real data 
c_hat = ccf(xbis[1,],xbis[2,],ylab = "cross-correlation",plot=FALSE)
str(c_hat)
cc_hat = c_hat$acf
n_cc = length(cc_hat)
y_cc_mar1full=matrix(0,100,n_cc)
### Use 100 values to get a simulation enveloppe 
x=matrix(0,2,t_max)
for(nrep in 1:100)
{
  ### Initial values
  x[,1]=xbis[,1]
  for (t in 1:(t_max-1))
  {
    mu=c(0,0)
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B %*% x[,t] + epsilon
  }
  
  y_cc_mar1full[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
}

matplot(c_hat$lag,t(y_cc_mar1full), pch=".",ylab="Cross-correlation",xlab="Time lag",main = "MAR(1) - full interactions")
matlines(c_hat$lag,t(y_cc_mar1full))
lines(c_hat$lag,cc_hat,col="black",lwd=4)
add_label_legend(label="B",cex=1.8)

### ---- Simulation MAR(2) indep model ### 
data.mar2.indep=read.csv(file="mar2/mar2.indep.csv")
data.mar2.indep
B1 = diag(data.mar2.indep$value[1:2])
B2 = diag(data.mar2.indep$value[3:4],2,2)
Sigma = diag(data.mar2.indep[5:6,]$value) ## Less variability on 
y_cc_mar2indep=matrix(0,100,n_cc)
x=matrix(0,2,t_max)
for(nrep in 1:100)
{
  ### Initial values
  x[,1]=xbis[,1]
  x[,2]=xbis[,2]
  for (t in 2:(t_max-1))
  {
    mu=c(0,0)
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
  }
  
  y_cc_mar2indep[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
}
matplot(c_hat$lag,t(y_cc_mar2indep), pch=".",ylab="Cross-correlation",xlab="Time lag",main = "MAR(2) - no interactions")
matlines(c_hat$lag,t(y_cc_mar2indep))
lines(c_hat$lag,cc_hat,col="black",lwd=4)
add_label_legend(label="C",cex=1.8)

### ---- Simulation MAR(2) bottom-up model ### 
data.mar2.bottom.up=read.csv(file="mar2/mar2.bottom.up.csv")
data.mar2.bottom.up
B1 = diag(data.mar2.bottom.up$value[1:2],2,2)
B2 = matrix(c(data.mar2.bottom.up$value[3],0,data.mar2.bottom.up$value[4:5]),2,2,byrow=TRUE)
Sigma = diag(data.mar2.bottom.up[6:7,]$value) ## Less variability on 
y_cc_mar2bottomup=matrix(0,100,n_cc)
x=matrix(0,2,t_max)
for(nrep in 1:100)
{
  ### Initial values
  x[,1]=xbis[,1]
  x[,2]=xbis[,2]
  for (t in 2:(t_max-1))
  {
    mu=c(0,0)
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
  }
  
  y_cc_mar2bottomup[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
}
matplot(c_hat$lag,t(y_cc_mar2bottomup), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - cycles & bottom-up predator dyn.")
matlines(c_hat$lag,t(y_cc_mar2bottomup))
lines(c_hat$lag,cc_hat,col="black",lwd=4)
add_label_legend(label="D",cex=1.8)

### ---- Simulation MAR(2) bottom-up variant model ### 
data.mar2.BUv3=read.csv(file="mar2/mar2.BUv3.csv")
data.mar2.BUv3
B1 = matrix(c(data.mar2.BUv2$value[1],0,data.mar2.BUv2$value[2:3]),2,2,byrow = TRUE)
B2 = matrix(c(data.mar2.BUv2$value[4],0,0,data.mar2.BUv2$value[5]),2,2,byrow=TRUE)
Sigma = diag(data.mar2.bottom.up[6:7,]$value) ## Less variability on 
y_cc_mar2bottomup=matrix(0,100,n_cc)
x=matrix(0,2,t_max)
for(nrep in 1:100)
{
  ### Initial values
  x[,1]=xbis[,1]
  x[,2]=xbis[,2]
  for (t in 2:(t_max-1))
  {
    mu=c(0,0)
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
  }
  
  y_cc_mar2bottomup[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
}
matplot(c_hat$lag,t(y_cc_mar2bottomup), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - Bottom-up variant")
matlines(c_hat$lag,t(y_cc_mar2bottomup))
lines(c_hat$lag,cc_hat,col="black",lwd=4)
add_label_legend(label="E",cex=1.8)


### ---- Simulation MAR(2) full model ### 
data.mar2.full=read.csv(file="mar2/mar2.full.csv")
B1 = matrix(data.mar2.full$value[1:4],2,2)
B2 = matrix(data.mar2.full$value[5:8],2,2)
Sigma = diag(data.mar2.full[9:10,]$value) 
y_cc_mar2full=matrix(0,100,n_cc)
x=matrix(0,2,t_max)
for(nrep in 1:100)
{
  ### Initial values
  x[,1]=xbis[,1]
  x[,2]=xbis[,2]
  for (t in 2:(t_max-1))
  {
    mu=c(0,0)
    epsilon=mvrnorm(n = 1, mu, Sigma)
    x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
  }
  
  y_cc_mar2full[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
}
matplot(c_hat$lag,t(y_cc_mar2full), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - full model")
matlines(c_hat$lag,t(y_cc_mar2full))
lines(c_hat$lag,cc_hat,col="black",lwd=4)
add_label_legend(label="F",cex=1.8)

dev.off()
# 
# # Start plotting --
# pdf("crossCorrelations_simulatedModels_moreBUmodels.pdf",height=20,width=16)
# par(mfrow=c(3,2),cex=1.4)#,cex=1.5
# # --- First panel -----
# 
# ### ---- Simulation MAR(2) bottom-up model ### 
# data.mar2.bottom.up=read.csv(file="mar2/mar2.bottom.up.csv")
# data.mar2.bottom.up
# B1 = diag(data.mar2.bottom.up$value[1:2],2,2)
# B2 = matrix(c(data.mar2.bottom.up$value[3],0,data.mar2.bottom.up$value[4:5]),2,2,byrow=TRUE)
# Sigma = diag(data.mar2.bottom.up[6:7,]$value) ## Less variability on 
# y_cc_mar2bottomup=matrix(0,100,n_cc)
# x=matrix(0,2,t_max)
# for(nrep in 1:100)
# {
#   ### Initial values
#   x[,1]=xbis[,1]
#   x[,2]=xbis[,2]
#   for (t in 2:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
#   }
#   
#   y_cc_mar2bottomup[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
# }
# matplot(c_hat$lag,t(y_cc_mar2bottomup), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - cycles & bottom-up predator dyn.")
# matlines(c_hat$lag,t(y_cc_mar2bottomup))
# lines(c_hat$lag,cc_hat,col="black",lwd=4)
# add_label_legend(label="A",cex=1.8)
# 
# ### ---- Simulation MAR(2) bottom-up variant model ### 
# B1 = diag(data.mar2.bottom.up.variant1$value[1:2],2,2)
# B2 = matrix(c(data.mar2.bottom.up$value[3],data.mar2.bottom.up$value[4],0,0),2,2)
# Sigma = diag(data.mar2.bottom.up[6:7,]$value) ## Less variability on 
# y_cc_mar2bottomup_variant=matrix(0,100,n_cc)
# x=matrix(0,2,t_max)
# for(nrep in 1:100)
# {
#   ### Initial values
#   x[,1]=xbis[,1]
#   x[,2]=xbis[,2]
#   for (t in 2:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
#   }
#   
#   y_cc_mar2bottomup_variant[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
# }
# matplot(c_hat$lag,t(y_cc_mar2bottomup_variant), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - BU variant 2, b_22^2:=0")
# matlines(c_hat$lag,t(y_cc_mar2bottomup_variant))
# lines(c_hat$lag,cc_hat,col="black",lwd=4)
# add_label_legend(label="B",cex=1.8)
# 
# ##### BU variant 1
# 
# data.mar2.BUv1=read.csv(file="mar2/mar2.BUv1.csv")
# data.mar2.BUv1
# B1 = matrix(c(data.mar2.BUv2$value[1],0,data.mar2.BUv2$value[2:3]),2,2,byrow = T)
# B2 = matrix(c(data.mar2.BUv2$value[4],0,data.mar2.BUv2$value[5:6]),2,2,byrow = T)
# Sigma = diag(data.mar2.bottom.up[7:8,]$value) ## Less variability on 
# y_cc_mar2bottomup=matrix(0,100,n_cc)
# x=matrix(0,2,t_max)
# for(nrep in 1:100)
# {
#   ### Initial values
#   x[,1]=xbis[,1]
#   x[,2]=xbis[,2]
#   for (t in 2:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
#   }
#   
#   y_cc_mar2bottomup[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
# }
# matplot(c_hat$lag,t(y_cc_mar2bottomup), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - BU v1")
# matlines(c_hat$lag,t(y_cc_mar2bottomup))
# lines(c_hat$lag,cc_hat,col="black",lwd=4)
# add_label_legend(label="C",cex=1.8)
# 
# 
# ##### BU variant 2
# 
# data.mar2.BUv2=read.csv(file="mar2/mar2.BUv2.csv")
# data.mar2.BUv2
# B1 = matrix(c(data.mar2.BUv2$value[1],0,data.mar2.BUv2$value[2:3]),2,2,byrow = T)
# B2 = matrix(c(data.mar2.BUv2$value[4],0,0,0),2,2,byrow=TRUE)
# Sigma = diag(data.mar2.bottom.up[5:6,]$value) ## Less variability on 
# y_cc_mar2bottomup=matrix(0,100,n_cc)
# x=matrix(0,2,t_max)
# for(nrep in 1:100)
# {
#   ### Initial values
#   x[,1]=xbis[,1]
#   x[,2]=xbis[,2]
#   for (t in 2:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
#   }
#   
#   y_cc_mar2bottomup[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
# }
# matplot(c_hat$lag,t(y_cc_mar2bottomup), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - BU v2")
# matlines(c_hat$lag,t(y_cc_mar2bottomup))
# lines(c_hat$lag,cc_hat,col="black",lwd=4)
# add_label_legend(label="D",cex=1.8)
# 
# ##### BU variant 3
# 
# data.mar2.BUv3=read.csv(file="mar2/mar2.BUv3.csv")
# data.mar2.BUv3
# B1 = matrix(c(data.mar2.BUv2$value[1],0,data.mar2.BUv2$value[2:3]),2,2,byrow = T)
# B2 = matrix(c(data.mar2.BUv2$value[4],0,0,data.mar2.BUv2$value[5]),2,2,byrow=TRUE)
# Sigma = diag(data.mar2.bottom.up[6:7,]$value) ## Less variability on 
# y_cc_mar2bottomup=matrix(0,100,n_cc)
# x=matrix(0,2,t_max)
# for(nrep in 1:100)
# {
#   ### Initial values
#   x[,1]=xbis[,1]
#   x[,2]=xbis[,2]
#   for (t in 2:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
#   }
#   
#   y_cc_mar2bottomup[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
# }
# matplot(c_hat$lag,t(y_cc_mar2bottomup), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - BU v3")
# matlines(c_hat$lag,t(y_cc_mar2bottomup))
# lines(c_hat$lag,cc_hat,col="black",lwd=4)
# add_label_legend(label="E",cex=1.8)
# 
# ### ---- Simulation MAR(2) full model ### 
# data.mar2.full=read.csv(file="mar2/mar2.full.csv")
# B1 = matrix(data.mar2.full$value[1:4],2,2)
# B2 = matrix(data.mar2.full$value[5:8],2,2)
# Sigma = diag(data.mar2.full[9:10,]$value) 
# y_cc_mar2full=matrix(0,100,n_cc)
# x=matrix(0,2,t_max)
# for(nrep in 1:100)
# {
#   ### Initial values
#   x[,1]=xbis[,1]
#   x[,2]=xbis[,2]
#   for (t in 2:(t_max-1))
#   {
#     mu=c(0,0)
#     epsilon=mvrnorm(n = 1, mu, Sigma)
#     x[,t+1]=B1 %*% x[,t] +B2 %*% x[,t-1] + epsilon
#   }
#   
#   y_cc_mar2full[nrep,]=as.vector(ccf(x[1,],x[2,],ylab = "cross-correlation",plot=F)$acf)
# }
# matplot(c_hat$lag,t(y_cc_mar2full), pch=".",ylab="Cross-correlation",xlab="Time lag",main="MAR(2) - full model")
# matlines(c_hat$lag,t(y_cc_mar2full))
# lines(c_hat$lag,cc_hat,col="black",lwd=4)
# add_label_legend(label="F",cex=1.8)
# 
# dev.off()
#End of cross-correlations for simulations under the fitted models 

#RQ: Should I do the same thing for models with covariates? 

###################################################################################
###### 2. Can the models be correctly identified - given the time series length? 
###################################################################################

######################### Explo analysis on single time series ############################
### Simulate the data according to a MAR(1) and see which model fits best
sim.data=MARSSsimulate(mar1.full, nsim=1, tSteps=100)$sim.data
# We simulate 100 years as a first try. 

### NB Other thing
residuals(mar2.indep.temp)$model.residuals ## why the f** are these 0? 

plot(1:100,sim.data[1,,],type="o")
lines(1:100,sim.data[2,,],type="o",col="red")

xsim=matrix(c(sim.data[1,,1],sim.data[2,,1]),nrow=2,byrow=T) ## Put that into a matrix

### Now try to MAR(2) and MAR(1) analyze this data

#MAR(1)
Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
### Initial values
pi1=matrix(0,2,1); #Initial values
V1=diag(1,2)

### Process model part
### Setting matrices
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
Q1="diagonal and unequal"

# Estimation
model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
mar1.full.sim=MARSS(xsim[,2:ncol(xsim)], model=model.list)
CIs.mar1.full.sim=MARSSparamCIs(mar1.full.sim) 

### Now MAR(2)
### State-space, observation part - never changes 
Z=matrix(c(1,0,0,1,0,0,0,0),2,4) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
A=matrix(0,2,1)### Intercept state space = 0
R=matrix(0,2,2) ### Error matrix state-space = 0

### Initial values
V=matrix(0,4,4)
pi=matrix(c(xsim[,2],xsim[,1]),4,1)
pi

### Interaction matrix
B1=matrix(list("b11_1","b12_1","b21_1","b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2","b21_2","b22_2"),2,2,byrow = T)
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

### Process error matrix
U=matrix(0,4,1)
Q=matrix(list(0),4,4)
Q[1,1]="q11"; Q[2,2]="q22"
Q
#stacked_data=rbind(xbis[,2:nrow(DGP)],xbis[,1:(nrow(DGP)-1)])

### Model call
model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.full.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
CIs.mar2.full.sim=MARSSparamCIs(mar2.full.sim)

### Interaction matrix
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2",0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B
model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.bottom.up.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
CIs.mar2.bottom.up.sim=MARSSparamCIs(mar2.bottom.up.sim) 

############### Old results -- different simulation ########################################################"
### Warning: takes ages --
#MARSSaic(mar1.full.sim, output = "AICbp") ### AIC: 410.4315   AICc: 410.8713   AICbp(param): 412.8782 
#MARSSaic(mar2.full.sim, output = "AICbp") ### AIC: 410.7696   AICc: 411.9461   AICbp(param): 415.0599   
#MARSSaic(mar2.bottom.up.sim, output = "AICbp") ### AIC: 427.4286   AICc: 428.0181   AICbp(param): 429.6933   
### So here we select quite clearly the right model (here, the MAR(1) model) with AICbp. 
#############################################################################################################


############################## New check: whether the range of values in simulated models are OK ###

### The data on occupancy was first logged and then standardized
x2new=(xsim[2,])*s2+m2
new_occupancy=exp(x2new)

plot(DGP$Year,DGP$Occupancy,type="b",main="Percentage territories occupied Gyrfalcon")
plot(new_occupancy,type="o") ## a bit too high but nothing wrong with the order of magnitude
plot(new_occupancy*80,type="o") ## between 60 and 100 birds, nothing crazy. 

sim.data2=MARSSsimulate(mar2.bottom.up, nsim=1, tSteps=100)$sim.data
xsim2=matrix(c(sim.data2[1,,1],sim.data2[2,,1]),nrow=2,byrow=T) ## Put that into a matrix
x2new=(xsim2[2,])*s2+m2
new_occupancy2=exp(x2new)
### Same thing for the bottom-up model
plot(DGP$Year,DGP$Occupancy,type="b",main="Percentage territories occupied Gyrfalcon")
plot(new_occupancy2,type="o") ## a bit too high but nothing wrong with the order of magnitude
plot(new_occupancy2*80,type="o") ## between 60 and 100 birds, nothing crazy. 
#########################################################################################################


#############################################################################################
### More formal test: simulate 100 times the data according to a MAR(1) and MAR(2) bottom-up
### Then evaluate which model fits best most of the time. 
### Is a scenario (top-down or bottom-up) easier or more difficult to evaluate?
### Can one be identified with 35 years and another needs 100 - or more? 
#############################################################################################


######## Now simulate for only 35 years #####################################################

#### Make a figure showing the simulated predator-prey and the simulated bottom-up
sim.data=MARSSsimulate(mar1.full, nsim=1000, tSteps=35)$sim.data ### MAR(1) full
sim.data2=MARSSsimulate(mar2.bottom.up, nsim=1000, tSteps=35)$sim.data ### MAR(2) bottom-up

### Simulated 1000 times because the estimates for 100 repeats were unstable. 

### Plot both models for two repeats - difficult to tell. 
pdf("SimulatedModels35ts.pdf",width=6,height=8)
par(pch=19,cex=1.5,lwd=3,mfrow=c(2,1))
plot(1:35,sim.data[1,,1],type="o",xlab="Years",ylab="Stdized pop. density")
lines(1:35,sim.data[2,,1],type="o",col="red")
plot(1:35,sim.data2[1,,2],type="o",xlab="Years",ylab="Stdized pop. density")
lines(1:35,sim.data2[2,,2],type="o",col="blue")
dev.off()

### Now try fit MAR(1) and MAR(2) to analyze this data

# Initializing IC criteria
AIC_mar1_sim1=AIC_mar2_sim1=AIC_mar1_sim2=AIC_mar2_sim2=NA
AICc_mar1_sim1=AICc_mar2_sim1=AICc_mar1_sim2=AICc_mar2_sim2=NA
BIC_mar1_sim1=BIC_mar2_sim1=BIC_mar1_sim2=BIC_mar2_sim2=NA

IC_simData_MAR=data.frame(AIC_mar1_sim1,AIC_mar2_sim1,AIC_mar1_sim2,AIC_mar2_sim2,AICc_mar1_sim1,AICc_mar2_sim1,AICc_mar1_sim2,AICc_mar2_sim2,BIC_mar1_sim1,BIC_mar2_sim1,BIC_mar1_sim2,BIC_mar2_sim2)

for (krep in 1:1000){ # for all the repeats
  
  #Temporary data structure
  AIC_mar1_sim1=AIC_mar2_sim1=AIC_mar1_sim2=AIC_mar2_sim2=NA
  AICc_mar1_sim1=AICc_mar2_sim1=AICc_mar1_sim2=AICc_mar2_sim2=NA
  BIC_mar1_sim1=BIC_mar2_sim1=BIC_mar1_sim2=BIC_mar2_sim2=NA
  #Store all those in a dataframe
  IC_simData_MAR_temp=data.frame(AIC_mar1_sim1,AIC_mar2_sim1,AIC_mar1_sim2,AIC_mar2_sim2,AICc_mar1_sim1,AICc_mar2_sim1,AICc_mar1_sim2,AICc_mar2_sim2,BIC_mar1_sim1,BIC_mar2_sim1,BIC_mar1_sim2,BIC_mar2_sim2)
  
  ######### 1. Fit MAR(1) model to MAR(1) simulation
  # MAR(1) definition 
  Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
  A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
  R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
  ### Initial values
  pi1=matrix(0,2,1); #Initial values
  V1=diag(1,2)
  ### Process model part
  ### Setting matrices
  B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
  U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
  Q1="diagonal and unequal"
  
  # Estimation
  model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
  xsim=matrix(c(sim.data[1,,krep],sim.data[2,,krep]),nrow=2,byrow=T) ## already defined above
  mar1.full.sim=MARSS(xsim[,2:ncol(xsim)], model=model.list) #MARSSparamCIs(mar1.full.sim) 
  IC_simData_MAR_temp$AIC_mar1_sim1=mar1.full.sim$AIC
  IC_simData_MAR_temp$AICc_mar1_sim1=mar1.full.sim$AICc
  IC_simData_MAR_temp$BIC_mar1_sim1=mar.bic(mar1.full.sim)
  
  ######## 2. Fit MAR(2) bottom-up sim to MAR(1) simulation
  ### State-space, observation part - never changes 
  Z=matrix(c(1,0,0,1,0,0,0,0),2,4) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
  A=matrix(0,2,1)### Intercept state space = 0
  R=matrix(0,2,2) ### Error matrix state-space = 0
  
  ### Initial values
  V=matrix(0,4,4)
  pi=matrix(c(xsim[,2],xsim[,1]),4,1)
  pi
  
  B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
  B2=matrix(list("b11_2",0,"b21_2","b22_2"),2,2,byrow = T) 
  B=matrix(list(0),4,4)
  B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
  B
  
  ### Process error matrix
  U=matrix(0,4,1)
  Q=matrix(list(0),4,4)
  Q[1,1]="q11"; Q[2,2]="q22"
  Q
  
  ### Model call
  model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
  #xsim=matrix(c(sim.data[1,,1],sim.data[2,,1]),nrow=2,byrow=T) ## already defined above
  mar2.bottom.up.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
  ###MARSSparamCIs(mar2.bottom.up.sim) 
  IC_simData_MAR_temp$AIC_mar2_sim1=mar2.bottom.up.sim$AIC
  IC_simData_MAR_temp$AICc_mar2_sim1=mar2.bottom.up.sim$AICc
  IC_simData_MAR_temp$BIC_mar2_sim1=mar.bic(mar2.bottom.up.sim)
  
  ######## 3. Fit MAR(1) model to MAR(2) bottom-up sim
  # MAR(1) definition 
  Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
  A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
  R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
  ### Initial values
  pi1=matrix(0,2,1); #Initial values
  V1=diag(1,2)
  ### Process model part
  ### Setting matrices
  B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
  U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
  Q1="diagonal and unequal"
  
  # Estimation
  model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
  xsim2=matrix(c(sim.data2[1,,krep],sim.data2[2,,krep]),nrow=2,byrow=T) ## Put that into a matrix
  mar1.full.sim2=MARSS(xsim2[,2:ncol(xsim2)], model=model.list) #MARSSparamCIs(mar1.full.sim) 
  IC_simData_MAR_temp$AIC_mar1_sim2=mar1.full.sim2$AIC
  IC_simData_MAR_temp$AICc_mar1_sim2=mar1.full.sim2$AICc
  IC_simData_MAR_temp$BIC_mar1_sim2=mar.bic(mar1.full.sim2)
  
  ######## 4. Fit MAR(2) bottom-up sim to MAR(2) sim
  ### State-space, observation part - never changes 
  Z=matrix(c(1,0,0,1,0,0,0,0),2,4) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
  A=matrix(0,2,1)### Intercept state space = 0
  R=matrix(0,2,2) ### Error matrix state-space = 0
  
  ### Initial values
  V=matrix(0,4,4)
  pi=matrix(c(xsim2[,2],xsim2[,1]),4,1)
  pi
  
  B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
  B2=matrix(list("b11_2",0,"b21_2","b22_2"),2,2,byrow = T) 
  B=matrix(list(0),4,4)
  B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
  B
  
  ### Process error matrix
  U=matrix(0,4,1)
  Q=matrix(list(0),4,4)
  Q[1,1]="q11"; Q[2,2]="q22"
  Q
  
  ### Model call
  model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
  #xsim2=matrix(c(sim.data2[1,,1],sim.data2[2,,1]),nrow=2,byrow=T) ## already defined. 
  mar2.bottom.up.sim2=MARSS(xsim2[,2:ncol(xsim2)],model=model.list.2lags)
  ###MARSSparamCIs(mar2.bottom.up.sim) 
  IC_simData_MAR_temp$AIC_mar2_sim2=mar2.bottom.up.sim2$AIC
  IC_simData_MAR_temp$AICc_mar2_sim2=mar2.bottom.up.sim2$AICc
  IC_simData_MAR_temp$BIC_mar2_sim2=mar.bic(mar2.bottom.up.sim2)
  
  ### Store data
  if (krep==1){IC_simData_MAR=IC_simData_MAR_temp} else {IC_simData_MAR=rbind(IC_simData_MAR,IC_simData_MAR_temp)}
} # end of loop on krep

write.csv(format(IC_simData_MAR,digits=4),file="aic.table.simulated.t35.csv")
### Choice by AIC
# Proportion of correct assignation of model to sim MAR(1)
sum(IC_simData_MAR$AIC_mar2_sim1>IC_simData_MAR$AIC_mar1_sim1)/1000
#0.52 # with error
#0.468     # corrected
# Proportion of correct assignation of model to sim MAR(2)
sum(IC_simData_MAR$AIC_mar1_sim2>IC_simData_MAR$AIC_mar2_sim2)/1000
#0.984 # with error
# 0.994      # corrected

### Choice by AICc
# Proportion of correct assignation of model to sim MAR(1)
sum(IC_simData_MAR$AICc_mar2_sim1>IC_simData_MAR$AICc_mar1_sim1)/1000
#0.565 # with error
#0.511     # corrected

# Proportion of correct assignation of model to sim MAR(2)
sum(IC_simData_MAR$AICc_mar1_sim2>IC_simData_MAR$AICc_mar2_sim2)/1000
#0.981 # with error
#0.993 # corrected

### Choice by BIC
# Proportion of correct assignation of model to sim MAR(1)
sum(IC_simData_MAR$BIC_mar2_sim1>IC_simData_MAR$BIC_mar1_sim1)/1000
#0.64 # with error
#0.586 # corrected

# Proportion of correct assignation of model to sim MAR(2)
sum(IC_simData_MAR$BIC_mar1_sim2>IC_simData_MAR$BIC_mar2_sim2)/1000
#0.97 # with error
#0.988 # corrected

hist(IC_simData_MAR$BIC_mar1_sim2) # just a check

####################################################################################
### Exact same analyses with tSim = 100 ############################################
####################################################################################

#### Make a figure showing the simulated predator-prey and the simulated bottom-up
sim.data=MARSSsimulate(mar1.full, nsim=1000, tSteps=100)$sim.data ### MAR(1) full
sim.data2=MARSSsimulate(mar2.bottom.up, nsim=1000, tSteps=100)$sim.data ### MAR(2) bottom-up

# Initializing IC criteria
AIC_mar1_sim1=AIC_mar2_sim1=AIC_mar1_sim2=AIC_mar2_sim2=NA
AICc_mar1_sim1=AICc_mar2_sim1=AICc_mar1_sim2=AICc_mar2_sim2=NA
BIC_mar1_sim1=BIC_mar2_sim1=BIC_mar1_sim2=BIC_mar2_sim2=NA

IC_simData_MAR=data.frame(AIC_mar1_sim1,AIC_mar2_sim1,AIC_mar1_sim2,AIC_mar2_sim2,AICc_mar1_sim1,AICc_mar2_sim1,AICc_mar1_sim2,AICc_mar2_sim2,BIC_mar1_sim1,BIC_mar2_sim1,BIC_mar1_sim2,BIC_mar2_sim2)

for (krep in 1:1000){ # for all the repeats
  
  #Temporary data structure
  AIC_mar1_sim1=AIC_mar2_sim1=AIC_mar1_sim2=AIC_mar2_sim2=NA
  AICc_mar1_sim1=AICc_mar2_sim1=AICc_mar1_sim2=AICc_mar2_sim2=NA
  BIC_mar1_sim1=BIC_mar2_sim1=BIC_mar1_sim2=BIC_mar2_sim2=NA
  #Store all those in a dataframe
  IC_simData_MAR_temp=data.frame(AIC_mar1_sim1,AIC_mar2_sim1,AIC_mar1_sim2,AIC_mar2_sim2,AICc_mar1_sim1,AICc_mar2_sim1,AICc_mar1_sim2,AICc_mar2_sim2,BIC_mar1_sim1,BIC_mar2_sim1,BIC_mar1_sim2,BIC_mar2_sim2)
  
  ######### 1. Fit MAR(1) model to MAR(1) simulation
  # MAR(1) definition 
  Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
  A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
  R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
  ### Initial values
  pi1=matrix(0,2,1); #Initial values
  V1=diag(1,2)
  ### Process model part
  ### Setting matrices
  B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
  U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
  Q1="diagonal and unequal"
  
  # Estimation
  model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
  xsim=matrix(c(sim.data[1,,krep],sim.data[2,,krep]),nrow=2,byrow=T) ## already defined above
  mar1.full.sim=MARSS(xsim[,2:ncol(xsim)], model=model.list) #MARSSparamCIs(mar1.full.sim) 
  IC_simData_MAR_temp$AIC_mar1_sim1=mar1.full.sim$AIC
  IC_simData_MAR_temp$AICc_mar1_sim1=mar1.full.sim$AICc
  IC_simData_MAR_temp$BIC_mar1_sim1=mar.bic(mar1.full.sim)
  
  ######## 2. Fit MAR(2) bottom-up sim to MAR(1) simulation
  ### State-space, observation part - never changes 
  Z=matrix(c(1,0,0,1,0,0,0,0),2,4) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
  A=matrix(0,2,1)### Intercept state space = 0
  R=matrix(0,2,2) ### Error matrix state-space = 0
  
  ### Initial values
  V=matrix(0,4,4)
  pi=matrix(c(xsim[,2],xsim[,1]),4,1)
  pi
  
  B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
  B2=matrix(list("b11_2",0,"b21_2","b22_2"),2,2,byrow = T) 
  B=matrix(list(0),4,4)
  B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
  B
  
  ### Process error matrix
  U=matrix(0,4,1)
  Q=matrix(list(0),4,4)
  Q[1,1]="q11"; Q[2,2]="q22"
  Q
  
  ### Model call
  model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
  #xsim=matrix(c(sim.data[1,,1],sim.data[2,,1]),nrow=2,byrow=T) ## already defined above
  mar2.bottom.up.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
  ###MARSSparamCIs(mar2.bottom.up.sim) 
  IC_simData_MAR_temp$AIC_mar2_sim1=mar2.bottom.up.sim$AIC
  IC_simData_MAR_temp$AICc_mar2_sim1=mar2.bottom.up.sim$AICc
  IC_simData_MAR_temp$BIC_mar2_sim1=mar.bic(mar2.bottom.up.sim)
  
  ######## 3. Fit MAR(1) model to MAR(2) bottom-up sim
  # MAR(1) definition 
  Z1=diag(1,2)                 ### Diagonal matrix from intrinsic to observed variables
  A1=matrix(list(0,0),2,1)     ### Intercept state space = 0
  R1=matrix(list(0,0,0,0),2,2) ### Error matrix state-space = 0
  ### Initial values
  pi1=matrix(0,2,1); #Initial values
  V1=diag(1,2)
  ### Process model part
  ### Setting matrices
  B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
  U1=matrix(0,2,1)                  ### Intercept is zero because data is centered. 
  Q1="diagonal and unequal"
  
  # Estimation
  model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)
  xsim2=matrix(c(sim.data2[1,,krep],sim.data2[2,,krep]),nrow=2,byrow=T) ## Put that into a matrix
  mar1.full.sim2=MARSS(xsim2[,2:ncol(xsim2)], model=model.list) #MARSSparamCIs(mar1.full.sim) 
  IC_simData_MAR_temp$AIC_mar1_sim2=mar1.full.sim2$AIC
  IC_simData_MAR_temp$AICc_mar1_sim2=mar1.full.sim2$AICc
  IC_simData_MAR_temp$BIC_mar1_sim2=mar.bic(mar1.full.sim2)
  
  ######## 4. Fit MAR(2) bottom-up sim to MAR(2) sim
  ### State-space, observation part - never changes 
  Z=matrix(c(1,0,0,1,0,0,0,0),2,4) ### Diagonal matrix from intrinsic to observed variables, 0 for delayed variables
  A=matrix(0,2,1)### Intercept state space = 0
  R=matrix(0,2,2) ### Error matrix state-space = 0
  
  ### Initial values
  V=matrix(0,4,4)
  pi=matrix(c(xsim2[,2],xsim2[,1]),4,1)
  pi
  
  B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
  B2=matrix(list("b11_2",0,"b21_2","b22_2"),2,2,byrow = T) 
  B=matrix(list(0),4,4)
  B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
  B
  
  ### Process error matrix
  U=matrix(0,4,1)
  Q=matrix(list(0),4,4)
  Q[1,1]="q11"; Q[2,2]="q22"
  Q
  
  ### Model call
  model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
  #xsim2=matrix(c(sim.data2[1,,1],sim.data2[2,,1]),nrow=2,byrow=T) ## already defined. 
  mar2.bottom.up.sim2=MARSS(xsim2[,2:ncol(xsim2)],model=model.list.2lags)
  ###MARSSparamCIs(mar2.bottom.up.sim) 
  IC_simData_MAR_temp$AIC_mar2_sim2=mar2.bottom.up.sim2$AIC
  IC_simData_MAR_temp$AICc_mar2_sim2=mar2.bottom.up.sim2$AICc
  IC_simData_MAR_temp$BIC_mar2_sim2=mar.bic(mar2.bottom.up.sim2)
  
  ### Store data
  if (krep==1){IC_simData_MAR=IC_simData_MAR_temp} else {IC_simData_MAR=rbind(IC_simData_MAR,IC_simData_MAR_temp)}
} # end of loop on krep

write.csv(format(IC_simData_MAR,digits=4),file="aic.table.simulated.t100.csv")
### Choice by AIC
# Proportion of correct assignation of model to sim MAR(1)
sum(IC_simData_MAR$AIC_mar2_sim1>IC_simData_MAR$AIC_mar1_sim1)/1000
#0.91 # with error
#0.866 # corrected

# Proportion of correct assignation of model to sim MAR(2)
sum(IC_simData_MAR$AIC_mar1_sim2>IC_simData_MAR$AIC_mar2_sim2)/1000
#1  # with error
#1 # corrected

### Choice by AICc
# Proportion of correct assignation of model to sim MAR(1)
sum(IC_simData_MAR$AICc_mar2_sim1>IC_simData_MAR$AICc_mar1_sim1)/1000
#0.92 # with error
#0.871 # corrected

# Proportion of correct assignation of model to sim MAR(2)
sum(IC_simData_MAR$AICc_mar1_sim2>IC_simData_MAR$AICc_mar2_sim2)/1000
#1 # with error
#1 # corrected

### Choice by BIC
# Proportion of correct assignation of model to sim MAR(1)
sum(IC_simData_MAR$BIC_mar2_sim1>IC_simData_MAR$BIC_mar1_sim1)/1000
#0.95 # with error
#0.937    # corrected

# Proportion of correct assignation of model to sim MAR(2)
sum(IC_simData_MAR$BIC_mar1_sim2>IC_simData_MAR$BIC_mar2_sim2)/1000
#1 # with error
#0.998 # corrected

### analyses stop here ### 

