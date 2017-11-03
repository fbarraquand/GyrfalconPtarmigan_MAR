#### New fig 1

# Old code for model in JAGS
mu=out3$BUGSoutput$mean$mu[,]
#plots including predictions one step ahead
par(mfrow=c(2,1))
plot(DGP$Year,xbis[1,],type="o",pch = 16,bg="black",xlab="Year",ylab="Ptarmigan density")
lines(DGP$Year[2:34],mu[1,],type="p")
arrows(DGP$Year[1:33],xbis[1,1:33],DGP$Year[2:34],mu[1,],length = 0.05)
plot(DGP$Year,xbis[2,],type="o",col="red",pch = 16,bg = "red",xlab="Year",ylab="Gyrfalcon log occupancy")
lines(DGP$Year[2:34],mu[2,],type="p",col="red")
arrows(DGP$Year[1:33],xbis[2,1:33],DGP$Year[2:34],mu[2,],col="red",length = 0.05)
#plots including predictions - scaled back to the normal scale
par(mfrow=c(2,1))
plot(DGP$Year,exp(xbis[1,]*s1+m1),type="o",pch = 16,bg="black",xlab="Year",ylab="Ptarmigan density")
lines(DGP$Year[2:34],exp(mu[1,]*s1+m1),type="p")
arrows(DGP$Year[1:33],exp(xbis[1,1:33]*s1+m1),DGP$Year[2:34],exp(mu[1,]*s1+m1),length = 0.05)
plot(DGP$Year,exp(xbis[2,]*s2+m2),type="o",col="red",pch = 16,bg = "red",xlab="Year",ylab="Gyrfalcon occupancy")
lines(DGP$Year[2:34],exp(mu[2,]*s2+m2),type="p",col="red")
arrows(DGP$Year[1:33],exp(xbis[2,1:33]*s2+m2),DGP$Year[2:34],exp(mu[2,]*s2+m2),col="red",length = 0.05)

### New code 

### Should I do it with the median (perhaps first) or attempt to derive a prediction interval (should be Gaussian)
# https://www.rdocumentation.org/packages/MARSS/versions/3.9/topics/predict.marssMLE

