################################################################################
library(cyclar)
################################################################################

# Simulate time-series
p <- ar2.parms(k=3.5, v=2)
x <- ar2.sim(n=500, mu=2.1, a1=p$a1, a2=p$a2, sd=sqrt(0.2))
ts.summary(x)

x.ar <- ar.univar(x, log=FALSE, constant=0)
x.ar$coef
p

# Convert numeric series to binary sequence
x.bin <- ts.binary(x)
ts.summary(x)
ts.summary(x.bin,plot=TRUE)
ts.spec(x)

# Check the acf's
par(mfrow=c(1,2))
	acf(x, main="Simulated series")
	acf(x.bin, main="Binary series")
par(mfrow=c(1,1))

# Check correlation between acf for simulated and binary converted series
# Repeat the simulation nrep number of times
nrep <- 500
rsq <- numeric(nrep) # Vector to hold output values (r-squared)

# Simulation parameters
mu <- 2.1
p <- ar2.parms(k=4, v=2)
a1 <- p$a1
a2 <- p$a2
ar.var <- 0.2

for (i in 1:nrep){
	x <- ar2.sim(n=50, mu=mu, a1=a1, a2=a2, sd=sqrt(ar.var))
	x.ar <- ar.univar(x, log=FALSE, constant=0)
	x.bin <- ts.binary(x)

	dat.bin <- data.frame(x = acf(x, plot=FALSE)$acf, y = acf(x.bin, plot=FALSE)$acf)
	rsq[i] <- summary(lm(y ~ x, data=dat.bin))$r.squared
	#plot(dat.bin[,"x"], dat.bin[,"y"], xlab = "Numeric series", ylab="Binary series", font.lab=2, las=1)
}

hist(sqrt(rsq),breaks=20, col="lightgrey",
	main="Correlation between ACF from numeric and binary ts", font.lab=2, xlab="Correlation coefficient")

mean(rsq)
mean(sqrt(rsq))

# Simulate ~3.5 year cycle, convert to binary series and check distribution of dominant lag among nrep replicate simulations
# Increasing the variance v in ar2.parms(k,v) results in higher proportions of "mis-classifications".
nrep <- 1000

# Simulation parameters
mu <- 0
p <- ar2.parms(k=3.5, v=3)
a1 <- p$a1
a2 <- p$a2
ar.var <- 0.2

out <- sapply(1:nrep, function(i) {
	x <- ar2.sim(n=50, mu=mu, a1=a1, a2=a2, sd=sqrt(ar.var))
	x.bin <- ts.binary(x)
	temp <- ts.summary(x.bin, base="raw")
	temp }
)

acf.bin.tab <- table(sapply(1:nrep, function(i) as.numeric(out[[i]][7])))/nrep
sum(acf.bin.tab)

plot(acf.bin.tab, xlab="Dominant lag", ylab="Proportion", font.lab=2, las=1)

