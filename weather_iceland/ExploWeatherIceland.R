### Exploratory analyses of monthly data weather from Iceland - FBarraquand 03/12/2016

################## Biological significance of the variables ###########################
### March-April for Gyrfalcon (April in Olafur's report) and June-July for ptarmigan. 
### Positive effect of high temp in June on Ptarmigan from Watson and Moss Ecology 2000
### Would appear with a two-year timelag. 
### Looks like March-April weather is influential for Gyr from Olafur's analyses
### Lag of 3 to 4 years are expected, given production is affected. 
### Positive effect for temp, negative for precipitation. 
#######################################################################################

rm(list=ls())
graphics.off()

#Stod_422_Akureyri.ManMedal.txt
#Stod_479_Manarbakki.ManMedal.txt

DB=read.table("Stod_495_Grimsstadir.ManMedal.txt",header=T,encoding = "latin1")
head(DB)
ndb=names(DB)
names(DB)[1:3]=c("site","year","month")
names(DB)

### March weather
min_t = min(DB$tnn[DB$month==3],na.rm=T)
max_t = max(DB$txx[DB$month==3],na.rm=T)
plot(DB$year[DB$month==3],DB$t[DB$month==3],type="o",ylim=c(min_t,max_t))
lines(DB$year[DB$month==3],DB$tx[DB$month==3],type="o",col="orange")
lines(DB$year[DB$month==3],DB$tn[DB$month==3],type="o",col="cyan")
lines(DB$year[DB$month==3],DB$txx[DB$month==3],type="o",col="red")
lines(DB$year[DB$month==3],DB$tnn[DB$month==3],type="o",col="blue")
### Patterns of average temp are consistent with extremes. 
### Lowest temp are more variable though. 

### Correlation with April
plot(DB$year[DB$month==3],DB$t[DB$month==3],type="o",ylim=c(min_t,max_t))
lines(DB$year[DB$month==4],DB$t[DB$month==4],type="o",col="green")
### Some relation but nonlinear...


plot(DB$year[DB$month==3],DB$r[DB$month==3],type="o")
plot(DB$t[DB$month==3],DB$r[DB$month==3])
# Precipitation... Slight trend towards colder years with more snow I guess. 

### June weather
min_t = min(DB$tnn[DB$month==6],na.rm=T)
max_t = max(DB$txx[DB$month==6],na.rm=T)
plot(DB$year[DB$month==6],DB$t[DB$month==6],type="o",ylim=c(min_t,max_t))
lines(DB$year[DB$month==6],DB$tx[DB$month==6],type="o",col="orange")
lines(DB$year[DB$month==6],DB$tn[DB$month==6],type="o",col="cyan")
lines(DB$year[DB$month==6],DB$txx[DB$month==6],type="o",col="red")
lines(DB$year[DB$month==6],DB$tnn[DB$month==6],type="o",col="blue")
### The highest temp are more variable. 
plot(DB$year[DB$month==6],DB$r[DB$month==6],type="o")

### Check, say June trends at the other stations
DB2=read.table("Stod_422_Akureyri.ManMedal.txt",header=T,encoding = "latin1")
DB3=read.table("Stod_479_Manarbakki.ManMedal.txt",header=T,encoding = "latin1")
names(DB2)[1:3]=c("site","year","month")
names(DB3)[1:3]=c("site","year","month")

### March
min_t = min(DB$tnn[DB$month==3],na.rm=T)
max_t = max(DB3$txx[DB$month==3],na.rm=T)
plot(DB$year[DB$month==3],DB$t[DB$month==3],type="o",ylim=c(min_t,max_t))
lines(DB2$year[DB2$month==3],DB2$t[DB2$month==3],type="o",col="green")
lines(DB3$year[DB3$month==3],DB3$t[DB3$month==3],type="o",col="blue")
# lines(DB$year[DB$month==3],DB$tx[DB$month==3],type="o",col="black",lty=2)
# lines(DB2$year[DB2$month==3],DB2$tx[DB2$month==3],type="o",col="green",lty=2)
# lines(DB3$year[DB3$month==3],DB3$tx[DB3$month==3],type="o",col="blue",lty=2)
# #
lines(DB$year[DB$month==3],DB$txx[DB$month==3],type="o",col="black",lty=3)
lines(DB2$year[DB2$month==3],DB2$txx[DB2$month==3],type="o",col="green",lty=3)
lines(DB3$year[DB3$month==3],DB3$txx[DB3$month==3],type="o",col="blue",lty=3)
#
# lines(DB$year[DB$month==3],DB$tn[DB$month==3],type="o",col="black",lty=2)
# lines(DB2$year[DB2$month==3],DB2$tn[DB2$month==3],type="o",col="green",lty=2)
# lines(DB3$year[DB3$month==3],DB3$tn[DB3$month==3],type="o",col="blue",lty=2)
#
lines(DB$year[DB$month==3],DB$tnn[DB$month==3],type="o",col="black",lty=3)
lines(DB2$year[DB2$month==3],DB2$tnn[DB2$month==3],type="o",col="green",lty=3)
lines(DB3$year[DB3$month==3],DB3$tnn[DB3$month==3],type="o",col="blue",lty=3)

### April
min_t = min(DB$tnn[DB$month==4],na.rm=T)
max_t = max(DB$txx[DB$month==4],na.rm=T)
plot(DB$year[DB$month==4],DB$t[DB$month==4],type="o",ylim=c(min_t,max_t))
lines(DB2$year[DB2$month==4],DB2$t[DB2$month==4],type="o",col="green")
lines(DB3$year[DB3$month==4],DB3$t[DB3$month==4],type="o",col="blue")
### The blue (Manarbakki, tip of a peninsula) is closer to green, more oceanic as well. 

### June
min_t = min(DB$tnn[DB$month==6],na.rm=T)
max_t = max(DB$txx[DB$month==6],na.rm=T)
plot(DB$year[DB$month==6],DB$t[DB$month==6],type="o",ylim=c(min_t,max_t))
lines(DB2$year[DB2$month==6],DB2$t[DB2$month==6],type="o",col="green")
lines(DB3$year[DB3$month==6],DB3$t[DB3$month==6],type="o",col="blue")
### This time Grim and Manabakki are more synchronous. 

### Need to compare precipitation as well. 
### April
### 
plot(DB$year[DB$month==4],DB$r[DB$month==4],type="o",ylim=c(0,100))
lines(DB2$year[DB2$month==4],DB2$r[DB2$month==4],type="o",col="green")
lines(DB3$year[DB3$month==4],DB3$r[DB3$month==4],type="o",col="blue")

lines(DB$year[DB$month==4],DB$rx[DB$month==4],type="o",col="black",lty=2)
lines(DB2$year[DB2$month==4],DB2$rx[DB2$month==4],type="o",col="green",lty=2)
lines(DB3$year[DB3$month==4],DB3$rx[DB3$month==4],type="o",col="blue",lty=2)

lines(DB$year[DB$month==4],DB$rh[DB$month==4],type="o",col="black",lty=3)
lines(DB2$year[DB2$month==4],DB2$rh[DB2$month==4],type="o",col="green",lty=3)
lines(DB3$year[DB3$month==4],DB3$rh[DB3$month==4],type="o",col="blue",lty=3)

### June
plot(DB$year[DB$month==6],DB$r[DB$month==6],type="o",ylim=c(0,100))
lines(DB2$year[DB2$month==6],DB2$r[DB2$month==6],type="o",col="green")
lines(DB3$year[DB3$month==6],DB3$r[DB3$month==6],type="o",col="blue")

lines(DB$year[DB$month==6],DB$rx[DB$month==6],type="o",col="black",lty=2)
lines(DB2$year[DB2$month==6],DB2$rx[DB2$month==6],type="o",col="green",lty=2)
lines(DB3$year[DB3$month==6],DB3$rx[DB3$month==6],type="o",col="blue",lty=2)

lines(DB$year[DB$month==6],DB$rh[DB$month==6],type="o",col="black",lty=3)
lines(DB2$year[DB2$month==6],DB2$rh[DB2$month==6],type="o",col="green",lty=3)
lines(DB3$year[DB3$month==6],DB3$rh[DB3$month==6],type="o",col="blue",lty=3)

### Will we need to use more integrated indices as well? (NAO and the likes?)

### Need to compare with all the months (use Date)
### to get a feel of within- vs between-year variability

DB$date=paste(DB$year,DB$month,15,sep="-")
DB2$date=paste(DB2$year,DB2$month,15,sep="-")
DB3$date=paste(DB3$year,DB3$month,15,sep="-")
as.Date(DB$date)
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o")
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green")
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue")

##### Plot only spring, from March to June
plot(as.Date(DB$date[DB$year>1980 & (DB$month>2&DB$month<7)]),DB$t[DB$year>1980& (DB$month>2&DB$month<7)],type="o")
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month>2&DB2$month<7)]),DB2$t[DB2$year>1980 & (DB2$month>2&DB2$month<7)],type="o",col="green")
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month>2&DB3$month<7)]),DB3$t[DB3$year>1980 & (DB3$month>2&DB3$month<7)],type="o",col="blue")

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),DB$t[DB$year>1980& (DB$month==4)],col="black",pch=16)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),DB2$t[DB2$year>1980 & (DB2$month==4)],col="green",pch=16)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),DB3$t[DB3$year>1980 & (DB3$month==4)],col="blue",pch=16)

### I could make a line appear connecting the April dots, overlayed on the general trend
pdf("Temp_NEIceland_April.pdf",width = 14,height=8)
### All the months
par(cex=1.3)
##### Plot all
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on April)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),DB$t[DB$year>1980& (DB$month==4)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),DB2$t[DB2$year>1980 & (DB2$month==4)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),DB3$t[DB3$year>1980 & (DB3$month==4)],col="blue",pch=16,lwd=3)
dev.off()

### Same with June
##### Plot all
pdf("Temp_NEIceland_June.pdf",width = 14,height=8)
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on June)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for June weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==6)]),DB$t[DB$year>1980& (DB$month==6)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==6)]),DB2$t[DB2$year>1980 & (DB2$month==6)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==6)]),DB3$t[DB3$year>1980 & (DB3$month==6)],col="blue",pch=16,lwd=3)
dev.off()

spectrum(DB$t,method = "ar")
## First peak at 0.083=1/12 OK
## Second peak at 0.17, about 1/6, what is it, echo. 
## Very little peak at 0.045 = 22 months
plot(spectrum(DB$t,method = "ar"),xlim=c(0,0.1))

spectrum(DB$t[DB$month==4],method = "ar")
spectrum(DB$t[DB$month==6],method = "ar") #Peak at 3.5, 4 years?
### Should I use wavelets there? Perhaps fishing for patterns - be careful
### Not much signal

pdf("Temp_NEIceland_March.pdf",width = 14,height=8)
### All the months
par(cex=1.3)
##### Plot all
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on March)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==3)]),DB$t[DB$year>1980& (DB$month==3)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==3)]),DB2$t[DB2$year>1980 & (DB2$month==3)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==3)]),DB3$t[DB3$year>1980 & (DB3$month==3)],col="blue",pch=16,lwd=3)
dev.off()



spectrum(DB$t, spans = 3) #robust result


### lty=2?

### Same for precipitation
##### Plot only spring, from March to June
plot(as.Date(DB$date[DB$year>1980 & (DB$month>2&DB$month<7)]),DB$r[DB$year>1980& (DB$month>2&DB$month<7)],type="o",ylim=c(0,120))
### Some problems there??? 
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month>2&DB2$month<7)]),DB2$r[DB2$year>1980 & (DB2$month>2&DB2$month<7)],type="o",col="green")
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month>2&DB3$month<7)]),DB3$r[DB3$year>1980 & (DB3$month>2&DB3$month<7)],type="o",col="blue")

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),DB$r[DB$year>1980& (DB$month==4)],col="black",pch=20)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),DB2$r[DB2$year>1980 & (DB2$month==4)],col="green",pch=20)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),DB3$r[DB3$year>1980 & (DB3$month==4)],col="blue",pch=20)

### All the months
plot(as.Date(DB$date[DB$year>1980]),DB$r[DB$year>1980],type="o",lty=2,ylim=c(0,120))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$r[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$r[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),DB$r[DB$year>1980& (DB$month==4)],col="black",pch=16)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),DB2$r[DB2$year>1980 & (DB2$month==4)],col="green",pch=16)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),DB3$r[DB3$year>1980 & (DB3$month==4)],col="blue",pch=16)


### Log-scale
pdf("LogRainfall_NEIceland_April.pdf",width = 14,height=8)
### All the months
par(cex=1.3)
plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on April)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),log(DB$r[DB$year>1980& (DB$month==4)]),col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),log(DB2$r[DB2$year>1980 & (DB2$month==4)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),log(DB3$r[DB3$year>1980 & (DB3$month==4)]),col="blue",pch=16,lwd=3)
dev.off()

pdf("LogRainfall_NEIceland_June.pdf",width = 14,height=8)
### Focus on June weather
par(cex=1.3)
plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on June)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

### Add points for June weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==6)]),log(DB$r[DB$year>1980& (DB$month==6)]),col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==6)]),log(DB2$r[DB2$year>1980 & (DB2$month==6)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==6)]),log(DB3$r[DB3$year>1980 & (DB3$month==6)]),col="blue",pch=16,lwd=3)
dev.off()

############### Compound plots ##################################
pdf("Temp_NEIceland.pdf",width = 14,height=8)
par(cex=1.3)#mar=c(4,1)

##### Plot all
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on March)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==3)]),DB$t[DB$year>1980& (DB$month==3)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==3)]),DB2$t[DB2$year>1980 & (DB2$month==3)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==3)]),DB3$t[DB3$year>1980 & (DB3$month==3)],col="blue",pch=16,lwd=3)

##### Plot all
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on April)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),DB$t[DB$year>1980& (DB$month==4)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),DB2$t[DB2$year>1980 & (DB2$month==4)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),DB3$t[DB3$year>1980 & (DB3$month==4)],col="blue",pch=16,lwd=3)

### Same with May
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on May)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for June weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),DB$t[DB$year>1980& (DB$month==5)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),DB2$t[DB2$year>1980 & (DB2$month==5)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),DB3$t[DB3$year>1980 & (DB3$month==5)],col="blue",pch=16,lwd=3)

### Same with June
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on June)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)

### Add points for June weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==6)]),DB$t[DB$year>1980& (DB$month==6)],col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==6)]),DB2$t[DB2$year>1980 & (DB2$month==6)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==6)]),DB3$t[DB3$year>1980 & (DB3$month==6)],col="blue",pch=16,lwd=3)
dev.off()

######### Rainfall
### Log-scale
pdf("LogRainfall_NEIceland.pdf",width = 14,height=8)
### All the months
par(cex=1.3)
# March
plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on March)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

### Add points for March weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==3)]),log(DB$r[DB$year>1980& (DB$month==3)]),col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==3)]),log(DB2$r[DB2$year>1980 & (DB2$month==3)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==3)]),log(DB3$r[DB3$year>1980 & (DB3$month==3)]),col="blue",pch=16,lwd=3)


plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on April)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

### Add points for April weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==4)]),log(DB$r[DB$year>1980& (DB$month==4)]),col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==4)]),log(DB2$r[DB2$year>1980 & (DB2$month==4)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==4)]),log(DB3$r[DB3$year>1980 & (DB3$month==4)]),col="blue",pch=16,lwd=3)

plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on May)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

### Add points for May weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),log(DB$r[DB$year>1980& (DB$month==5)]),col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),log(DB2$r[DB2$year>1980 & (DB2$month==5)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),log(DB3$r[DB3$year>1980 & (DB3$month==5)]),col="blue",pch=16,lwd=3)



### Focus on June weather
plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on June)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki"),pch=21,pt.bg=c("black","green","blue"))
lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

### Add points for June weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==6)]),log(DB$r[DB$year>1980& (DB$month==6)]),col="black",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==6)]),log(DB2$r[DB2$year>1980 & (DB2$month==6)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==6)]),log(DB3$r[DB3$year>1980 & (DB3$month==6)]),col="blue",pch=16,lwd=3)
dev.off()

######## More data on rainfall
DB4=read.table("Stod_468_Reykjahlid.ManMedal.txt",header=T,encoding = "latin1")
DB5=read.table("Stod_462_Myri.ManMedal.txt",header=T,encoding = "latin1")
DB6=read.table("Stod_473_Stadarholl.ManMedal.txt",header=T,encoding = "latin1")
DB7=read.table("Stod_502_Raufarhofn.ManMedal.txt",header=T,encoding = "latin1") # very coastal, a bit outside. 
DB8=read.table("Stod_502_Raufarhofn.ManMedal.txt",header=T,encoding = "latin1")
DB9=read.table("Stod_448_Lerkihlid.ManMedal.txt",header=T,encoding = "latin1")
### We remove Akureyri which is actually not in the study area. 

names(DB4)[1:3]=c("site","year","month")
names(DB5)[1:3]=c("site","year","month")
names(DB6)[1:3]=c("site","year","month")
names(DB7)[1:3]=c("site","year","month")
names(DB8)[1:3]=c("site","year","month")
names(DB9)[1:3]=c("site","year","month")
DB4$date=paste(DB4$year,DB4$month,15,sep="-")
DB5$date=paste(DB5$year,DB5$month,15,sep="-")
DB6$date=paste(DB6$year,DB6$month,15,sep="-")
DB7$date=paste(DB7$year,DB7$month,15,sep="-")
DB8$date=paste(DB8$year,DB8$month,15,sep="-")
DB9$date=paste(DB9$year,DB9$month,15,sep="-")

plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on May)")
legend("bottomright",c("Grimsstadir","Manarbakki","Reykjahlid","Myri","Stadarholl","Lerkihlid"),pch=21,pt.bg=c("black","blue","red","pink","orange","green"))
#lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)
lines(as.Date(DB4$date[DB4$year>1980]),log(DB4$r[DB5$year>1980]),type="o",col="red",lty=2)
lines(as.Date(DB5$date[DB5$year>1980]),log(DB5$r[DB5$year>1980]),type="o",col="pink",lty=2)
lines(as.Date(DB6$date[DB6$year>1980]),log(DB6$r[DB6$year>1980]),type="o",col="orange",lty=2)
lines(as.Date(DB9$date[DB9$year>1980]),log(DB9$r[DB9$year>1980]),type="o",col="green",lty=2)

### Add points for May weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),log(DB$r[DB$year>1980& (DB$month==5)]),col="black",pch=16,lwd=3)
#lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),log(DB2$r[DB2$year>1980 & (DB2$month==5)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),log(DB3$r[DB3$year>1980 & (DB3$month==5)]),col="blue",pch=16,lwd=3)
lines(as.Date(DB4$date[DB4$year>1980 & (DB4$month==5)]),log(DB4$r[DB4$year>1980 & (DB4$month==5)]),col="red",pch=16,lwd=3)
lines(as.Date(DB5$date[DB5$year>1980 & (DB5$month==5)]),log(DB5$r[DB5$year>1980 & (DB5$month==5)]),col="pink",pch=16,lwd=3)
lines(as.Date(DB6$date[DB6$year>1980 & (DB6$month==5)]),log(DB6$r[DB6$year>1980 & (DB6$month==5)]),col="orange",pch=16,lwd=3)
lines(as.Date(DB9$date[DB9$year>1980 & (DB9$month==5)]),log(DB9$r[DB9$year>1980 & (DB9$month==5)]),col="green",pch=16,lwd=3)

logRainfall=cbind(log(DB$r[DB$year>1980]),log(DB3$r[DB3$year>1980]),log(DB4$r[DB4$year>1980]),log(DB5$r[DB5$year>1980]),log(DB6$r[DB6$year>1980]),log(DB9$r[DB9$year>1980]))
#av_logRainfall=rowMeans(logRainfall,na.rm=T)
### OK, differing length there, we need to keep a data structure. 
#plot(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]))
#lines(as.Date(DB3$date[B3$year>1980]),av_logRainfall,type="b")

nrow(DB3[DB3$year>1980,])
nrow(DB4[DB4$year>1980,])#1 less
nrow(DB5[DB5$year>1980,])
nrow(DB6[DB6$year>1980,])#10 less
nrow(DB9[DB9$year>1980,])#1 more!

DB[nrow(DB),]
DB4[nrow(DB4),]

#Creating new data structure from 1980
DB_av=DB[DB$year>1980,]
logRainfall=data.frame(DB_av$date,DB_av$year,DB_av$month,DB_av$r)
logRainfall=cbind(logRainfall,DB3$r[DB$year>1980])
which(is.na(DB4$r))

cbind(DB3$date[DB3$year>1980],DB9$date[DB9$year>1980])
#what? 
DB3$year[DB3$year>1980]
DB3$date[DB3$year>1980]

cbind(DB3$year[DB3$year>1980],DB3$month[DB3$year>1980],DB3$date[DB3$year>1980])

plot(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)

c=rep(1,12)
for (i in 1:12){
c[i]=length(DB3$date[(DB3$month==i)&(DB3$year>1980)])
}
c
DB3$date[(DB3$month==i)&(DB3$year>1980)]# not yet updated. 

c=rep(1,12)
for (i in 1:12){
  c[i]=length(DB9$date[(DB9$month==i)&(DB9$year>1980)])
}
c
#WTF
DB9$date[(DB9$month==1)&(DB9$year>1980)] #with jan 2017

c=rep(1,12)
for (i in 1:12){
  c[i]=length(DB6$date[(DB6$month==i)&(DB6$year>1980)])
}
c
DB6$date[DB6$year>1980]
# e.g. 
DB6$date[(DB6$month==1)&(DB6$year>1980)] # no 2012
DB6$date[(DB6$month==6)&(DB6$year>1980)] # no 1989, no 2011


############# Allright, this is way too messy, we're going to loop through months and years to reconstruct the data ########

# Loop through years and then months
for (year in 1975:2016){
  for (month in 1:12){
  ### Take averages of temp and log Rainfall
    
  # Temp
  t1 = DB$t[(DB$year==year)& (DB$month==month)]
  if(length(t1)==0){t1=NA}
  t2 = DB2$t[(DB2$year==year)& (DB2$month==month)]
  if(length(t2)==0){t2=NA}
  t3 = DB3$t[(DB3$year==year)& (DB3$month==month)]
  if(length(t3)==0){t3=NA}
  t=c(t1,t2,t3)
  av.temp = mean(t,na.rm = T)
  
  # Rainfall 
  lr1 = log(DB$r[(DB$year==year)& (DB$month==month)])
  if(length(lr1)==0){lr1=NA}
  lr2 = log(DB3$r[(DB3$year==year)& (DB3$month==month)])
  if(length(lr2)==0){lr2=NA}
  lr3 = log(DB4$r[(DB4$year==year)& (DB4$month==month)])
  if(length(lr3)==0){lr3=NA}
  lr4 = log(DB5$r[(DB5$year==year)& (DB5$month==month)])
  if(length(lr4)==0){lr4=NA}
  lr5 = log(DB6$r[(DB6$year==year)& (DB6$month==month)])
  if(length(lr5)==0){lr5=NA}
  lr6 = log(DB9$r[(DB9$year==year)& (DB9$month==month)])
  if(length(lr6)==0){lr6=NA}
  # Removes numeric(0) and transform in NAs
  lr=c(lr1,lr2,lr3,lr4,lr5,lr6)
  av.logRainfall = mean(lr,na.rm = T)
  date=paste(year,month,15,sep="-")
  ### Construct data structure here
  if ((year== 1975)&(month==1)){
  av.weather=data.frame(date,year,month,av.temp,av.logRainfall)
  }
  else {# add to data structure
  new_row = data.frame(date,year,month,av.temp,av.logRainfall)
  av.weather=rbind(av.weather,new_row)  
    }
  }
}

#av.weather$av.temp[is.nan(av.weather$av.temp)]=NA
#av.weather$av.logRainfall[is.nan(av.weather$av.logRainfall)]=NA

# Verif
year=1987
month=4
lr1 = log(DB$r[(DB$year==year)& (DB$month==month)])
if(length(lr1)==0){lr1=NA}
lr2 = log(DB3$r[(DB3$year==year)& (DB3$month==month)])
if(length(lr2)==0){lr2=NA}
lr3 = log(DB4$r[(DB4$year==year)& (DB4$month==month)])
if(length(lr3)==0){lr3=NA}
lr4 = log(DB5$r[(DB5$year==year)& (DB5$month==month)])
if(length(lr4)==0){lr4=NA}
lr5 = log(DB6$r[(DB6$year==year)& (DB6$month==month)])
if(length(lr5)==0){lr5=NA}
lr6 = log(DB9$r[(DB9$year==year)& (DB9$month==month)])
if(length(lr6)==0){lr6=NA}
# Removes numeric(0) and transform in NAs
lr=c(lr1,lr2,lr3,lr4,lr5,lr6)
av.logRainfall = mean(lr,na.rm = T)
# seems OK

names(av.weather)[4:5]=c("temp","logRainfall")
head(av.weather)
write.csv(av.weather,"average_weatherNEIceland.csv")

### Plots
pdf("logRainfall_NEIceland_May_withAverage.pdf",width = 14,height=8)
plot(as.Date(DB$date[DB$year>1980]),log(DB$r[DB$year>1980]),type="o",lty=2,ylim=c(0,5.2),xlab="Time (years)",ylab="log(rainfall) (focus on May)")
legend("bottomright",c("Grimsstadir","Manarbakki","Reykjahlid","Myri","Stadarholl","Lerkihlid","Average"),pch=21,pt.bg=c("cyan","blue","red","pink","orange","green","black"))
#lines(as.Date(DB2$date[DB2$year>1980]),log(DB2$r[DB2$year>1980]),type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),log(DB3$r[DB3$year>1980]),type="o",col="blue",lty=2)
lines(as.Date(DB4$date[DB4$year>1980]),log(DB4$r[DB4$year>1980]),type="o",col="red",lty=2)
lines(as.Date(DB5$date[DB5$year>1980]),log(DB5$r[DB5$year>1980]),type="o",col="pink",lty=2)
lines(as.Date(DB6$date[DB6$year>1980]),log(DB6$r[DB6$year>1980]),type="o",col="orange",lty=2)
lines(as.Date(DB9$date[DB9$year>1980]),log(DB9$r[DB9$year>1980]),type="o",col="green",lty=2)
lines(as.Date(av.weather$date),av.weather$logRainfall,type="o",col="black",lty=2)

### Add points for May weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),log(DB$r[DB$year>1980& (DB$month==5)]),col="cyan",pch=16,lwd=3)
#lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),log(DB2$r[DB2$year>1980 & (DB2$month==5)]),col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),log(DB3$r[DB3$year>1980 & (DB3$month==5)]),col="blue",pch=16,lwd=3)
lines(as.Date(DB4$date[DB4$year>1980 & (DB4$month==5)]),log(DB4$r[DB4$year>1980 & (DB4$month==5)]),col="red",pch=16,lwd=3)
lines(as.Date(DB5$date[DB5$year>1980 & (DB5$month==5)]),log(DB5$r[DB5$year>1980 & (DB5$month==5)]),col="pink",pch=16,lwd=3)
lines(as.Date(DB6$date[DB6$year>1980 & (DB6$month==5)]),log(DB6$r[DB6$year>1980 & (DB6$month==5)]),col="orange",pch=16,lwd=3)
lines(as.Date(DB9$date[DB9$year>1980 & (DB9$month==5)]),log(DB9$r[DB9$year>1980 & (DB9$month==5)]),col="green",pch=16,lwd=3)
lines(as.Date(av.weather$date[av.weather$month==5]),av.weather$logRainfall[av.weather$month==5],type="o",col="black",lwd=3)
dev.off()


pdf("Temp_NEIceland_May_withAverage.pdf",width = 14,height=8)
### Same with temp
plot(as.Date(DB$date[DB$year>1980]),DB$t[DB$year>1980],type="o",lty=2,xlab="Time (years)",ylab="Average temperature (focus on May)")
legend("bottomright",c("Grimsstadir","Akureyri","Manarbakki","Average"),pch=21,pt.bg=c("cyan","green","blue","black"))
lines(as.Date(DB2$date[DB2$year>1980]),DB2$t[DB2$year>1980],type="o",col="green",lty=2)
lines(as.Date(DB3$date[DB3$year>1980]),DB3$t[DB3$year>1980],type="o",col="blue",lty=2)
lines(as.Date(av.weather$date),av.weather$logRainfall,type="o",col="black",lty=2)

### Add points for May weather
lines(as.Date(DB$date[DB$year>1980 & (DB$month==5)]),DB$t[DB$year>1980& (DB$month==5)],col="cyan",pch=16,lwd=3)
lines(as.Date(DB2$date[DB2$year>1980 & (DB2$month==5)]),DB2$t[DB2$year>1980 & (DB2$month==5)],col="green",pch=16,lwd=3)
lines(as.Date(DB3$date[DB3$year>1980 & (DB3$month==5)]),DB3$t[DB3$year>1980 & (DB3$month==5)],col="blue",pch=16,lwd=3)
lines(as.Date(av.weather$date[av.weather$month==5]),av.weather$temp[av.weather$month==5],type="o",col="black",lwd=3)
dev.off()
