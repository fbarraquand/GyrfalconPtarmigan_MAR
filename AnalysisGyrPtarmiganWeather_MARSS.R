################################################################################################################
################# Analysis Gyrfalcon-Partmigan-Weather data - MAR modelling with MARSS package #################
### FBarraquand 26/03/2017, analyses started 06/07/2015, with O. Nielsen #######################################
### Updated version producing all the figures and results of the paper in a reproducible workflow 25/07/2017 ###
################################################################################################################

### Initializing
rm(list=ls())
graphics.off()
### Setting the seed to keep the same simulations
set.seed(42) # What else?

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

### Now include the same full model but with a correlated noise matrix
Q1=matrix(c("q11","q21","q12","q22"),2,2) ##assume correlated noise
#model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1)#did not work at some point. 
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

### Diagnostics
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

# Storage full model
value=mar1.full$par$B#value=CIs.mar1.full$par$B
SE=CIs.mar1.full$par.se$B
lower=CIs.mar1.full$par.lowCI$B
upper=CIs.mar1.full$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.full$par$U ### C is U
SE=CIs.mar1.full$par.se$U
lower=CIs.mar1.full$par.lowCI$U
upper=CIs.mar1.full$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.full$par$Q ### C is U
SE=CIs.mar1.full$par.se$Q
lower=CIs.mar1.full$par.lowCI$Q
upper=CIs.mar1.full$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.full.csv")

mar1.null$par$B
#Storage null model
value=mar1.null$par$B#value=CIs.mar1.null$par$B
SE=CIs.mar1.null$par.se$B
lower=CIs.mar1.null$par.lowCI$B
upper=CIs.mar1.null$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.null$par$U
SE=CIs.mar1.null$par.se$U
lower=CIs.mar1.null$par.lowCI$U
upper=CIs.mar1.null$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.null$par$Q
SE=CIs.mar1.null$par.se$Q
lower=CIs.mar1.null$par.lowCI$Q
upper=CIs.mar1.null$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.null.csv")

mar1.temp$par$B
#Storage full model with effect of temperature
value=mar1.temp$par$B#value=CIs.mar1.temp$par$B
SE=CIs.mar1.temp$par.se$B
lower=CIs.mar1.temp$par.lowCI$B
upper=CIs.mar1.temp$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.temp$par$U
SE=CIs.mar1.temp$par.se$U
lower=CIs.mar1.temp$par.lowCI$U
upper=CIs.mar1.temp$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.temp$par$Q
SE=CIs.mar1.temp$par.se$Q
lower=CIs.mar1.temp$par.lowCI$Q
upper=CIs.mar1.temp$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.temp.csv")

covar=t(as.matrix(cbind(tempMay_year,tempApril_year_minus4)))
C1=matrix(list("tempMay_year",0,0,"tempApril_year_minus4"),2,2,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11",0,0,"b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.temp.only=MARSS(data, model=model.list)
CIs.mar1.temp.only=MARSSparamCIs(mar1.temp.only)

value=mar1.temp.only$par$B#value=CIs.mar1.temp.only$par$B
SE=CIs.mar1.temp.only$par.se$B
lower=CIs.mar1.temp.only$par.lowCI$B
upper=CIs.mar1.temp.only$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.temp.only$par$U
SE=CIs.mar1.temp.only$par.se$U
lower=CIs.mar1.temp.only$par.lowCI$U
upper=CIs.mar1.temp.only$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.temp.only$par$Q
SE=CIs.mar1.temp.only$par.se$Q
lower=CIs.mar1.temp.only$par.lowCI$Q
upper=CIs.mar1.temp.only$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.temp.only.csv")


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
value=mar1.temp1$par$B#value=CIs.mar1.temp1$par$B
SE=CIs.mar1.temp1$par.se$B
lower=CIs.mar1.temp1$par.lowCI$B
upper=CIs.mar1.temp1$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.temp1$par$U
SE=CIs.mar1.temp1$par.se$U
lower=CIs.mar1.temp1$par.lowCI$U
upper=CIs.mar1.temp1$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.temp1$par$Q
SE=CIs.mar1.temp1$par.se$Q
lower=CIs.mar1.temp1$par.lowCI$Q
upper=CIs.mar1.temp1$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.temp1.csv")

### Now a model with rainfall
covar=t(as.matrix(cbind(rainMay_year,rainApril_year_minus4)))
C1=matrix(list("rainMay_year",0,0,"rainApril_year_minus4"),2,2,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.rain=MARSS(data, model=model.list)
MARSSparamCIs(mar1.rain)
# No effect. 

#Storage full model with effect of (cumulated) rainfall
value=mar1.rain$par$B#value=CIs.mar1.rain$par$B
SE=CIs.mar1.rain$par.se$B
lower=CIs.mar1.rain$par.lowCI$B
upper=CIs.mar1.rain$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.rain$par$U
SE=CIs.mar1.rain$par.se$U
lower=CIs.mar1.rain$par.lowCI$U
upper=CIs.mar1.rain$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.rain$par$Q
SE=CIs.mar1.rain$par.se$Q
lower=CIs.mar1.rain$par.lowCI$Q
upper=CIs.mar1.rain$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.rain.csv")

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

value=mar1.both$par$B#value=CIs.mar1.both$par$B
SE=CIs.mar1.both$par.se$B
lower=CIs.mar1.both$par.lowCI$B
upper=CIs.mar1.both$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.both$par$U
SE=CIs.mar1.both$par.se$U
lower=CIs.mar1.both$par.lowCI$U
upper=CIs.mar1.both$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.both$par$Q
SE=CIs.mar1.both$par.se$Q
lower=CIs.mar1.both$par.lowCI$Q
upper=CIs.mar1.both$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.both.csv")

### First comparison of AICc/BIC for all models 
mar.bic <- function(my.mar)
{
  my.mar$BIC=0  
  my.mar$BIC=-2*my.mar$logLik + my.mar$num.params*log(my.mar$samp.size/2)
}

mar.list=list(mar1.null=mar1.null,mar1.full=mar1.full,mar1.temp=mar1.temp,mar1.temp.only=mar1.temp.only,mar1.temp1=mar1.temp1,mar1.rain=mar1.rain,mar1.both=mar1.both)
names(mar.list[1]) #check

for (k in 1:length(mar.list)){
  mar.model=mar.list[k]$AIC
  mar.model$BIC=mar.bic(mar.model)
  aic.table.temp=data.frame(mar.model$logLik,mar.model$AIC,mar.model$AICc,mar.model$BIC)
  if (k==1){
    aic.table=aic.table.temp
  } else {
    aic.table=rbind(aic.table,aic.table.temp)
  }
}
### Can't make that f*** loop work -- mar.model$AIC undefined

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
write.csv(aic.table,file="mar1/aic.table.csv")
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
value=mar1.both.june$par$B
SE=CIs.mar1.both.june$par.se$B
lower=CIs.mar1.both.june$par.lowCI$B
upper=CIs.mar1.both.june$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.both.june$par$U
SE=CIs.mar1.both.june$par.se$U
lower=CIs.mar1.both.june$par.lowCI$U
upper=CIs.mar1.both.june$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.both.june$par$Q
SE=CIs.mar1.both.june$par.se$Q
lower=CIs.mar1.both.june$par.lowCI$Q
upper=CIs.mar1.both.june$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.both.june.csv")

covar=t(as.matrix(cbind(tempJuly_year,rainJuly_year,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("tempJuly_year","rainJuly_year",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.both.july=MARSS(data, model=model.list)
CIs.mar1.both.july=MARSSparamCIs(mar1.both.july)

value=mar1.both.july$par$B
SE=CIs.mar1.both.july$par.se$B
lower=CIs.mar1.both.july$par.lowCI$B
upper=CIs.mar1.both.july$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar1.both.july$par$U
SE=CIs.mar1.both.july$par.se$U
lower=CIs.mar1.both.july$par.lowCI$U
upper=CIs.mar1.both.july$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar1.both.july$par$Q
SE=CIs.mar1.both.july$par.se$Q
lower=CIs.mar1.both.july$par.lowCI$Q
upper=CIs.mar1.both.july$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar1/mar1.both.july.csv")

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
write.csv(aic.table,file="mar1/aic.table.csv")

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

covar=t(as.matrix(cbind(avMonthly_tempWinter,log_winter_rain,tempApril_year_minus4,rainApril_year_minus4)))
C1=matrix(list("minOfMonths_tempWinter","log_winter_rain",0,0,0,0,"tempApril_year_minus4","rainApril_year_minus4"),2,4,byrow=T)
C1
Q1="diagonal and unequal"
B1=matrix(list("b11","b12","b21","b22"),2,2,byrow = T) ### Interaction matrix
model.list=list(B=B1,U=U1,C=C1,c=covar,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)
mar1.both.winter2=MARSS(data, model=model.list)
MARSSparamCIs(mar1.both.winter2)

#### Check again those models

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
value=mar2.full$par$B
SE=CIs.mar2.full$par.se$B
lower=CIs.mar2.full$par.lowCI$B
upper=CIs.mar2.full$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar2.full$par$U
SE=CIs.mar2.full$par.se$U
lower=CIs.mar2.full$par.lowCI$U
upper=CIs.mar2.full$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar2.full$par$Q
SE=CIs.mar2.full$par.se$Q
lower=CIs.mar2.full$par.lowCI$Q
upper=CIs.mar2.full$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar2/mar2.full.csv")


#### Bottom-up model that we highlighted before -- with an effect of prey on predator growth the next year. 
### Interaction matrix
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2",0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B

model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.bottom.up=MARSS(xbis[,2:nrow(DGP)],model=model.list.2lags)
CIs.mar2.bottom.up=MARSSparamCIs(mar2.bottom.up)

### Store data
value=mar2.bottom.up$par$B
SE=CIs.mar2.bottom.up$par.se$B
lower=CIs.mar2.bottom.up$par.lowCI$B
upper=CIs.mar2.bottom.up$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar2.bottom.up$par$U
SE=CIs.mar2.bottom.up$par.se$U
lower=CIs.mar2.bottom.up$par.lowCI$U
upper=CIs.mar2.bottom.up$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar2.bottom.up$par$Q
SE=CIs.mar2.bottom.up$par.se$Q
lower=CIs.mar2.bottom.up$par.lowCI$Q
upper=CIs.mar2.bottom.up$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar2/mar2.bottom.up.csv")

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
value=mar2.indep$par$B
SE=CIs.mar2.indep$par.se$B
lower=CIs.mar2.indep$par.lowCI$B
upper=CIs.mar2.indep$par.upCI$B
mar1.data=data.frame(value,SE,lower,upper)
value=mar2.indep$par$U
SE=CIs.mar2.indep$par.se$U
lower=CIs.mar2.indep$par.lowCI$U
upper=CIs.mar2.indep$par.upCI$U
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
value=mar2.indep$par$Q
SE=CIs.mar2.indep$par.se$Q
lower=CIs.mar2.indep$par.lowCI$Q
upper=CIs.mar2.indep$par.upCI$Q
mar1.data=rbind(mar1.data,data.frame(value,SE,lower,upper))
write.csv(mar1.data,file="mar2/mar2.indep.csv")

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

## MARSSkem warnings. Type MARSSinfo() for help.
## MARSSkem: The soln became unstable and logLik DROPPED.
### where the hell should I mention control$trace=1 ???

########### Check that --- 

##################################################################################################
###### Compare fit of models
##################################################################################################
mar2.full$AICc
mar2.bottom.up$AICc
mar2.indep$AICc 
mar1.full.bis$AICc
mar1.null.bis$AICc #even higher, how is it possible when compared to VAR? 

mar2.full$AIC
mar2.bottom.up$AIC
mar2.indep$AIC 
mar1.full.bis$AIC
mar1.null.bis$AIC
### They are quite close to each other, but clearly below  MAR(1)
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

mar2.indep$BIC=mar.bic(mar2.indep)
aic.table2.temp=data.frame(mar2.indep$logLik,mar2.indep$AIC,mar2.indep$AICc,mar2.indep$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)

mar2.indep.temp$BIC=mar.bic(mar2.indep.temp)
aic.table2.temp=data.frame(mar2.indep.temp$logLik,mar2.indep.temp$AIC,mar2.indep.temp$AICc,mar2.indep.temp$BIC)
names(aic.table2.temp)=c("logLik","AIC","AICc","BIC")
aic.table2=rbind(aic.table2,aic.table2.temp)


rownames(aic.table2)=c("mar1.null.bis","mar1.full.bis","mar2.full","mar2.bottom.up","mar2.indep","mar2.indep.temp")

aic.table2
write.csv(aic.table2,file="mar2/aic.table2.csv")

############ Previous estimates #####################################################################
#MARSSaic(mar2.full, output = "AICbp") ### AIC: 129.5646   AICc: 133.5646   AICbp(param): 143.6595 
#MARSSaic(mar2.bottom.up, output = "AICbp") ### AIC: 129.8294   AICc: 131.7604   AICbp(param): 137.406   
#MARSSaic(mar2.indep,output = "AICbp") ### AIC: 129.5638   AICc: 130.9876   AICbp(param): 133.7121   
#MARSSaic(mar1.full.bis,output = "AICbp") ### AIC: 140.8294   AICc: 142.2531   AICbp(param): 148.7717   
#MARSSaic(mar1.null.bis,output = "AICbp") ### AIC: 143.1355   AICc: 143.7913   AICbp(param): 147.0655   
#MARSSaic(mar2.indep.temp,output="AICbp") ### AIC: 129.5638   AICc: 130.9876   AICbp(param): 134.6404   
######################################################################################################

########################## Note and other analyses ###########"""""""###############################
#### I'll be tempted to do some Bayesian comparison for the MAR(1), MAR(2) indep and bottom-up, with classical diagonal covariance structure. 
# - And some fit of models with lm() too

growth_rate_bis = growth_rate[,2:ncol(x)]
x_current = x[,2:ncol(x)]
x_delayed =  x[,1:(ncol(x)-1)]

lm1=lm(growth_rate_bis[1,] ~ 0 + x_current[1,]+x_current[2,])# pb
lm2=lm(growth_rate_bis[2,] ~ 0 + x_current[1,]+x_current[2,])

lm1.delayed=lm(growth_rate_bis[1,] ~ 0 + x_current[1,]+x_current[2,]+x_delayed[1,]+x_delayed[2,])# pb
lm2.delayed=lm(growth_rate_bis[2,] ~ 0 + x_current[1,]+x_current[2,]+x_delayed[1,]+x_delayed[2,])
### These have few significant coefficients, but that's perhaps just the TS length. 

### Would it be that the delayed component of lm2 has to be removed?  

AIC(lm1,lm1.delayed)
BIC(lm1,lm1.delayed)
AIC(lm2,lm2.delayed)
BIC(lm2,lm2.delayed)

### BIC favors mostly the MAR(1) here, taking the other species as a covariate. 

################# We use now VAR to check the model order ############################################
library(vars)
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

### The AIC selects 3 lags and the BIC and HQ, that are more conservative, two lags. 

### Does this fit sees causality in the 2x2 model?

######## MAR(2) model #################
varpp2<-VAR(y=data.frame(t(xbis)), p=2, type="none")
causality(varpp2,cause="X2") ## No effect of X2 on X1
causality(varpp2,cause="X1") ## close to reject the hypothesis of non-GC

######## MAR(1) model ################# 
varpp1<-VAR(y=data.frame(t(xbis)), p=1, type="none") 
causality(varpp1,cause="X2") ## No effect of X2 on X1
causality(varpp1,cause="X1")
# Looks better for causality but that's not the model that's selected... 
#######################################################################################################

#####################################################################################
##################### Simulation-based model selection ##############################
####################################################################################

###### 1. Can we get the right cross-correlation patterns with the best-fitting models? 

# Output the cross-correlation pattern for the data

# Simulations of the fitted models 

# --- old code to improve --- #

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


# ----------------------------#


# Cross-correlations for simulations under the fitted models 

####### ---- Code to UPDATE --- #####
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

###### 2. Can the models be correctly identified - given the time series length? 

### Simulate the data according to a MAR(1) and see which model fits best
sim.data=MARSSsimulate(mar1.full, nsim=1, tSteps=100)$sim.data

### Other thing
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
MARSSparamCIs(mar1.full.sim) 

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

### Interaction matrix
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2",0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B
model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.bottom.up.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
MARSSparamCIs(mar2.bottom.up.sim) 
### MARSSparamCIs(mar2.bottom.up.sim)
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -206.7143 
# AIC: 427.4286   AICc: 428.0181   
# 
# ML.Est Std.Err low.CI    up.CI
# B.b11_1  0.87476  0.1027  0.673  1.07603
# B.b22_1  0.65843  0.0868  0.488  0.82853
# B.b11_2 -0.08725  0.1044 -0.292  0.11740
# B.b12_2 -0.18226  0.0879 -0.355 -0.00994
# B.b22_2 -0.00345  0.0633 -0.128  0.12064
# Q.q11    0.51215  0.0511  0.379  0.66559
# Q.q22    0.45477  0.0482  0.336  0.59101
# 
# CIs calculated at alpha = 0.05 via method=hessian 


############### Old results -- different simulation ########################################################"
MARSSaic(mar1.full.sim, output = "AICbp") ### AIC: 410.4315   AICc: 410.8713   AICbp(param): 412.8782 
MARSSaic(mar2.full.sim, output = "AICbp") ### AIC: 410.7696   AICc: 411.9461   AICbp(param): 415.0599   
MARSSaic(mar2.bottom.up.sim, output = "AICbp") ### AIC: 427.4286   AICc: 428.0181   AICbp(param): 429.6933   
### So here we select quite clearly the right model (here, the MAR(1) model) with AICbp. 
#############################################################################################################

### This tends to suggest that the bottom-up model is appropriate on the real data - or that we don't know. 
### I need to check the coefficients of these models in quite some details. 

############################## check whether the range of values in simulated models are OK ###

### The data on occupancy was first logged and then standardized
x2new=(xsim[2,]+m2)*s2
new_occupancy=exp(x2new)

plot(DGP$Year,DGP$Occupancy,type="b",main="Percentage territories occupied Gyrfalcon")
plot(new_occupancy,type="o") ## a bit too high but nothing wrong with the order of magnitude
plot(new_occupancy*80,type="o") ## between 60 and 100 birds, nothing crazy. 

sim.data2=MARSSsimulate(mar2.bottom.up, nsim=1, tSteps=100)$sim.data
xsim2=matrix(c(sim.data2[1,,1],sim.data2[2,,1]),nrow=2,byrow=T) ## Put that into a matrix
x2new=(xsim2[2,]+m2)*s2
new_occupancy2=exp(x2new)
### Same thing for the bottom-up model
plot(DGP$Year,DGP$Occupancy,type="b",main="Percentage territories occupied Gyrfalcon")
plot(new_occupancy2,type="o") ## a bit too high but nothing wrong with the order of magnitude
plot(new_occupancy2*80,type="o") ## between 60 and 100 birds, nothing crazy. 

######### Very difficult to differentiate based on those two models. 
# Ballpark estimate? We have ~ 100 predators, 100 000 prey at best on the study area. 
# But they decrease at best by 20000 in a single year. So about 200 have to be eaten by predators
# (we don't count reproduction here, but we do not count death by other causes either)
# Is this plausible?
# Would it make more sense to directly formulate a mechanistic model where there is a phenom component for the gyr
# and a mechanistic model - with predation included - for the ptarmigan. Perhaps based on Erla's models? 
# Knowing that there is also a lot of hunting, and other predators of ptarmigan. 

######## Now simulate for only 33 years ##################
### Other stuff, e.g. simulate the data according to a MAR(1) and see what happens

### Plot those
plot(1:35,sim.data[1,,],type="o")
lines(1:35,sim.data[2,,],type="o",col="red")
xsim=matrix(c(sim.data[1,,1],sim.data[2,,1]),nrow=2,byrow=T) ## Put that into a matrix

#### Make a figure showing the simulated predator-prey and the simulated bottom-up
sim.data=MARSSsimulate(mar1.full, nsim=1, tSteps=35)$sim.data
sim.data2=MARSSsimulate(mar2.bottom.up, nsim=1, tSteps=35)$sim.data
### Plot both
pdf("SimulatedModels35ts.pdf",width=6,height=8)
par(pch=19,cex=1.5,lwd=3,mfrow=c(2,1))
plot(1:35,sim.data[1,,],type="o",xlab="Years",ylab="Stdized pop. density")
lines(1:35,sim.data[2,,],type="o",col="red")
plot(1:35,sim.data2[1,,],type="o",xlab="Years",ylab="Stdized pop. density")
lines(1:35,sim.data2[2,,],type="o",col="red")
dev.off()

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
MARSSparamCIs(mar1.full.sim) 

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

### Interaction matrix
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2",0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B
model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.bottom.up.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
MARSSparamCIs(mar2.bottom.up.sim) 

MARSSaic(mar1.full.sim, output = "AICbp") ###AIC: 111.6506   AICc: 113.0277   AICbp(param): 116.3878   
# Coeffs
# Estimate
# B.b11           0.6389
# B.b21           0.0842
# B.b12          -0.5065
# B.b22           0.2907
# Q.(X.Y1,X.Y1)   0.2711
# Q.(X.Y2,X.Y2)   0.2236
MARSSaic(mar2.full.sim, output = "AICbp") ###AIC: 104.5234   AICc: 108.3831   AICbp(param): 115.959   
# Coeffs
# Estimate
# B.b11_1  0.96058
# B.b21_1  0.00368
# B.b12_1 -0.46595
# B.b22_1  0.29607
# B.b11_2 -0.27712
# B.b21_2  0.05380
# B.b12_2  0.53499
# B.b22_2 -0.17746
# Q.q11    0.20422
# Q.q22    0.21743
MARSSaic(mar2.bottom.up.sim, output = "AICbp") ### AIC: 104.5126   AICc: 106.3793   AICbp(param): 108.3304   
# Coeffs
# Estimate
# B.b11_1    1.103
# B.b22_1    0.316
# B.b11_2   -0.445
# B.b12_2    0.464
# B.b22_2   -0.203
# Q.q11      0.243
# Q.q22      0.219
mar1_sim_35steps_bic=mar.bic(mar1.full.sim)
mar2_sim_35steps_bic=mar.bic(mar2.full.sim)
mar2_bottomUp_35steps_bic=mar.bic(mar2.bottom.up.sim)
# > mar1_sim_35steps_bic
# [1] 120.8088
# > mar2_sim_35steps_bic
# [1] 119.7871
# > mar2_bottomUp_35steps_bic
# [1] 115.1971

###########################################################################################################################
### Argh -- on this very short time series length, MAR(1) vs MAR(2) models are not identifiable. 
### And this is not even something we can use generally to discuss results of Vik because they had 100 years of data... 
### Should I check whether we get the same thing with different simulations? YES

###### New simulation

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
MARSSparamCIs(mar1.full.sim) 

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

### Interaction matrix
B1=matrix(list("b11_1",0,0,"b22_1"),2,2,byrow = T)
B2=matrix(list("b11_2","b12_2",0,"b22_2"),2,2,byrow = T) 
B=matrix(list(0),4,4)
B[1:2,1:2]=B1; B[1:2,3:4]=B2; B[3:4,1:2]=diag(1,2)
B
model.list.2lags=list(Z=Z,B=B,U=U,Q=Q,A=A,R=R,x0=pi,V0=V,tinitx=1)
mar2.bottom.up.sim=MARSS(xsim[,2:ncol(xsim)],model=model.list.2lags)
MARSSparamCIs(mar2.bottom.up.sim) 

MARSSaic(mar1.full.sim, output = "AICbp") ### AIC: 111.6506   AICc: 113.0277   AICbp(param): 115.2834   
MARSSaic(mar2.full.sim, output = "AICbp") ### AIC: 104.5234   AICc: 108.3831   AICbp(param): 114.7151   
MARSSaic(mar2.bottom.up.sim, output = "AICbp") ### AIC: 104.5126   AICc: 106.3793   AICbp(param): 108.0672   

mar1_sim_35steps_bic=mar.bic(mar1.full.sim)
mar2_sim_35steps_bic=mar.bic(mar2.full.sim)
mar2_bottomUp_35steps_bic=mar.bic(mar2.bottom.up.sim)
# > mar1_sim_35steps_bic
# [1] 120.8088
# > mar2_sim_35steps_bic
# [1] 119.7871
# > mar2_bottomUp_35steps_bic
# [1] 115.1971

###### Conclusion (preliminary?) - likely that MAR(2) bottom-up model is best whenever we simulate with MAR(1) full. 
### NB I could try to analyze these data with lm() - would be interesting (but more of a modelling project?)

### rewriting the data structure in the exact same fashion as for the true data. 
growth_rate = xsim[,2: ncol(xsim)]-xsim[,1:(ncol(xsim)-1)]
xs=xsim[,1:(ncol(xsim)-1)]
growth_rate_bis = growth_rate[,2:ncol(xs)]
x_current = xsim[,2:ncol(xs)]
x_delayed =  xsim[,1:(ncol(xs)-1)]

lm1=lm(growth_rate_bis[1,] ~ 0 + x_current[1,]+x_current[2,])# pb
lm2=lm(growth_rate_bis[2,] ~ 0 + x_current[1,]+x_current[2,])

lm1.delayed=lm(growth_rate_bis[1,] ~ 0 + x_current[1,]+x_current[2,]+x_delayed[1,]+x_delayed[2,])# pb
lm2.delayed=lm(growth_rate_bis[2,] ~ 0 + x_current[1,]+x_current[2,]+x_delayed[1,]+x_delayed[2,])
### These have few significant coefficients, but that's perhaps just the TS length. 

### Would it be that the delayed component of lm2 has to be removed?  
AIC(lm1,lm1.delayed) #delayed best
BIC(lm1,lm1.delayed) #delayed best
AIC(lm2,lm2.delayed) #nondelayed best
BIC(lm2,lm2.delayed) #nondelayed best
### Inconsistent model selection but quite clear that AIC and BIC can select the delayed model even though the model is not delayed. 

########################################################################################################################################
########## Conclusion: A simulated MAR(1) model on 35 timesteps can be mistaken for a MAR(2) [need even stronger penalties for params?]
################## Stopped there // FB 14/04/2017 ######################################################################################



