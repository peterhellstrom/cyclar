library(cyclar)

# Make parameters a function of time:
# Example function

a.tv <- function(x,a,b,c) a / (1 + exp(-(x - b) / -c))

n <- 50
a1.t <- a.tv(x=0:n, a=1.3, b=-0.0025, c=10) # Direct d-d
a2.t <- a.tv(x=0:n, a=-0.8, b=30, c=8) # Delayed d-d
a0 <- 2.1
a0.t <- rep(2.1,(n+1))

# Investigate properties of function:
bi <- seq(-10,40,5)
ci <- seq(1,20,2)


plot(x=c(0,50),y=c(-1,1),type="n",xlab="Time",ylab="Parameter value", font.lab=2, las=1, main="b")
for (i in 1:length(bi)) curve(a.tv(x=x, a=-0.8, b=bi[i], c=8), from=0, to=50, n=101, col=i, lwd=2, add=T)
legend("topleft", legend=bi, col=1:length(bi), lwd=2, bty="n", title="b =")


plot(x=c(0,50),y=c(-1,1),type="n",xlab="Time",ylab="Parameter value", font.lab=2, las=1, main="c")
for (i in 1:length(ci)) curve(a.tv(x=x, a=-0.8, b=30, c=ci[i]), from=0, to=50, n=101, col=i, lwd=2, add=T)
legend("topleft", legend=ci, col=1:length(ci), lwd=2, bty="n", title="c =")

# Plot changes in parameter values over time:
par(mfrow=c(1,2))
plot(0:n, a1.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a1", bty="l")
plot(0:n, a2.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a2", bty="l")
par(mfrow=c(1,1))

ar2.plot(xlims=c(-2,2),ylims=c(-1,0),triangle=FALSE)
ar2.arrows(a1=a1.t, a2=a2.t)

# Estimate ar-coefficients from time-varying simulations:
ar2.plot(ylims=c(-1,0),triangle=FALSE,winw=7,winh=5)
nrep <- 100

for (i in 1:nrep) {

	x <- ar2.sim.tv(n=50, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2))$ts.x
	#x <- growth.rate(exp(x), lags=1)[,"R"]
	mod <- ar.univar(x, log=FALSE, constant=0)
	points(mod$coef[1,1], mod$coef[1,2], col=4, pch=16, cex=1)
}

ar2.arrows(a1=a1.t, a2=a2.t)

######
npatch <- 6
z <- sapply(1:npatch, function(i) ar2.sim.tv(n=50, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2))$ts.x)

# Plot
plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0.t,lty=2)
for (i in 1:ncol(z)) points(z[,i], type="l", col=i)

plot(x=c(0,n),y=range(exp(z)), type="n", xlab="Time", ylab="Population density", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0.t,lty=2)
for (i in 1:ncol(z)) points(exp(z[,i]), type="l", col=i)

# Plot changes in variability
z.s <- sapply(1:npatch, function(i) swin.var(x=z[,i], swin=6))
matplot(z.s,type="l",xlab="Time (center of time-window)",ylab="S-index",main="Changes in crude variability over time", font.lab=2, las=1, bty="l")

# Construct data for parameter estimation in WinBUGS:

dat <- lapply(1:ncol(z), function(i) growth.rate(exp(z[,i]),lags=1))
dat <- lapply(1:ncol(z), function(i) na.omit(dat[[i]]))
sapply(1:ncol(z), function(i) nrow(dat[[i]])) # Check number of observations per patch

growth <- sapply(1:ncol(z), function(i) dat[[i]][,"R"])
lag1 <- sapply(1:ncol(z), function(i) dat[[i]][,"Lag0"])
lag2 <- sapply(1:ncol(z), function(i) dat[[i]][,"Lag1"])

ntime <- nrow(growth)
npatch <- ncol(growth)
