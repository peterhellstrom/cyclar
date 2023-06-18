library(cyclar)

n <- 50

a1.t <- 0.7 - 0.03*(0:n)

a2 <- -0.7
a2.t <- rep(a2,(n+1)) # Delayed d-d

a0 <- 2.1
a0.t <- rep(a0,(n+1))


# Plot changes in parameter values over time:
par(mfrow=c(1,2))
plot(0:n, a1.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a1", bty="l")
plot(0:n, a2.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a2", bty="l")
par(mfrow=c(1,1))


npatch <- 10
z <- lapply(1:npatch, function(i) ar2.sim.tv(n=50, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2)))

plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0.t,lty=2)
for (i in 1:length(z)) points(z[[i]]$ts.x, type="l", col=i)


# Plot changes in variability
z.s <- sapply(1:npatch, function(i) swin.var(x=z[[i]]$ts.x, swin=6))
matplot(z.s,type="l",xlab="Time (center of time-window)",ylab="S-index",main="Changes in crude variability over time", font.lab=2, las=1, bty="l")

# Create input files (batch)
inp.list <- lapply(1:npatch, function(i) vc.inp(x=z[[i]], method="x", file = paste("vc_ar2_",i,".csv", sep="")))

# RUN ESTIMATION MANUALLY!!!
# "C:/Program Files/Varying Coefficients Estimation/VC.exe"

# Check output (batch)
out.list <- lapply(1:npatch, function(i) vc.out(vc.dat=inp.list[[i]], method="x", file = paste("vc_ar2_",i,"-VC.csv", sep="")))


# Simulated
ar2.plot()
ar2.arrows(a1.t, a2.t)

# Estimated (collect all in one plot)
ar2.plot()
for (i in 1:npatch) ar2.arrows(a1=out.list[[i]][,"a_x1"], a2=out.list[[i]][,"a_x2"], col=i)


# Check some "weird" series manually?
#ar2.plot()
#ar2.arrows(a1=out.list[[10]][,"a_x1"], a2=out.list[[10]][,"a_x2"], col=2)
#
#plot(z[[10]]$x)

