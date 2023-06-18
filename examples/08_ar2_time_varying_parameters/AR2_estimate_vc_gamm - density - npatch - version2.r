library(cyclar)

# SIMULATION
npoints <- 8

coord.mat <- matrix(
c(
0,-0.2,
0.75, -0.2,
0.75, -0.8,
0, -0.8,
-1, -0.8,
-1, -0.2,
0.4, 0.3,
0,-0.5
),
ncol=2, nrow=npoints, byrow=TRUE, dimnames=list(1:npoints,c("x","y")))


route <- c(1,4,6,5,1)
route <- c(4,8,1)
route <- c(1,2)
route <- c(2,8,5)


route <- c(5,1,3) # Short-no-long
route <- c(3,1,5) # Long-no-short

route <- c(3,4,5) # shortening, no change in a2, high var
route <- c(2,1,6) # shortening, no change in a2, low var

nsteps <- length(route)-1
#yrs <- 100
yrs <- 25
yrs <- 10

vals <- matrix(c(
	rep(yrs,nsteps),
	rep(0,nsteps),
	rep(0.2,nsteps)),
nrow=nsteps, ncol=3, byrow=FALSE)

ps <- ar2.sim.tv.gen(npoints=8, coord.mat=coord.mat, route=route, vals=vals, gam=TRUE)

a0.0 <- 2 # CHANGE THIS ONE, AFFECTS GAM FIT A LOT!!! MUCH WORSE FIT IF THIS VALUE >0!?
a0 <- rep(a0.0, length(ps$a1))

################################################################################
# DATA GENERATION & PARAMETER ESTIMATION
n <- length(a0)
npatch <- 30
ar.var <- 0.2
z <- sapply(1:npatch, function(i) ar2.sim.tv(n=n, a0=a0, a1=ps$a1, a2=ps$a2, sd=sqrt(ar.var), plot.log=TRUE, vals=vals)$ts.x)
z <- ts(z, start=1)

# Plot
plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0,lty=2)
for (i in 1:ncol(z)) points(z[,i], type="l", col=i)

################################################################################
# Construct data for parameter estimation:
lgs <- 2

dat <- lapply(1:ncol(z), function(i) growth.rate(exp(z[,i]),lags=lgs))
datr <- data.frame(do.call("rbind", dat))

inds <- 1:(n+lgs)
datr$time <- rep(inds,npatch)
datr$patch <- factor(rep(1:npatch, each=(n+lgs)))
datr <- na.omit(datr)

################################################################################
# ADDED A SMOOTHER FOR TREND ALSO, without that one and a0.0 > 0, fit is totally WRONG!?

b.gamm <- gamm(Lag0 ~ s(time, bs="ps", m=c(2,2), by=Lag1) + s(time, bs="ps", m=c(2,2), by=Lag2) + s(time),
			data=datr, random=list(patch=~1))

b.gamm2 <- gamm(Lag0 ~ s(time, bs="ps", m=c(2,2), by=Lag1) + s(time, bs="ps", m=c(2,2), by=Lag2),
			data=datr, random=list(patch=~1))

par(mfrow=c(1,3))
plot(b.gamm$gam)
par(mfrow=c(1,1))

b.gamm.fit <- ar2.vc.gam.fit(b.gamm$gam, method="raw", model="gamm")
b.gamm2.fit <- ar2.vc.gam.fit(b.gamm2$gam, method="raw", model="gamm")

ar2.vc.gam.plot(x=b.gamm.fit,fit="x",inds=c(1,2))

ar2.plot(k=c(3,4,5,7),windows=TRUE,text.lab=FALSE,ylims=c(-1,0),triangle=FALSE,winw=11,winh=7)
ar2.arrows(ps$a1,ps$a2)
ar2.arrows(b.gamm.fit[[1]][,"fit"], b.gamm.fit[[2]][,"fit"], col=4, length=0.05)
ar2.arrows(b.gamm2.fit[[1]][,"fit"], b.gamm2.fit[[2]][,"fit"], col=3, length=0.05)
legend("topright",legend=c("simulated","estimated + time","estimated wo time"), col=c(2,4,3), lty=c(1,1,1))

