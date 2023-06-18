#########################################################################
library(cyclar)
#########################################################################

# LINEAR TREND, LOESS and GAM models

x <- ts(runif(30), start=1979, deltat=1)

ts.detr(x)
ts.detr(x, span=1, add.mean=TRUE, constant=1, log=TRUE, title="test")

# Replace some values with NA
x[c(5)] <- NA
#x[c(5,21)] <- NA
ts.detr(x)

# ARIMA-model
x <- arima.sim(n=100, model=list(ar=c(0,-0.8), ma=c(-0.5,1)))
plot(x)

ts.detr(x)

# Detrend gypsy moth data

Year <- 1911:1934
EggMass <- c(5214,5407,2635,3658,4400,3751,3702,2273,4032,2134,2387,402,110,50,61,127,303,722,407,66,42,112,344,181)
density <- c(8.559102594,8.595449689,7.876638461,8.204671829,8.38935982,8.22977775,8.216628493,7.728855824,8.30201781,7.665753432,7.777792626,5.996452089,4.700480366,3.912023005,4.110873864,4.844187086,5.713732806,6.582025139,6.008813185,4.189654742,3.737669618,4.718498871,5.840641657,5.198497031)

gypsy <- data.frame(Year=Year, EggMass=EggMass, density=density)

gypsy.eggs <- ts(EggMass, start=Year[1])

plot(gypsy$Year, gypsy$EggMass, xlab="Year", ylab="Egg mass", font.lab=2, main="Original time series, all years", bty="l", cex=1.3, pch=16)
lines(gypsy$Year, gypsy$EggMass)

# Only analyze 1911-1921
gypsy.eggs <- window(gypsy.eggs,1911,1921)

# Following Berryman. 1999. Principles of Population Dynamics. p. 108-109
ts.detr(x=gypsy.eggs, span=0.75, type="mts")
ts.detr(x=gypsy.eggs, span=1, type="mts")

ts.detr(x=gypsy.eggs, add.mean=FALSE, span=1)

test <- ts.detr(x=gypsy.eggs, add.mean=FALSE, span=5, plot=FALSE, type="mts")
apply(test,2,mean)

# STEP-DETREND

# Example 1
ngr <- 3 # Number of groups
nrep <- c(10,10,10) # Number of data points in each group
g <- rep(1:ngr, times=nrep) # Grouping variable
ms <- c(3,6,3) # Means for each segment
sds <- c(0.8,0.8,0.8) # Standard deviations for each segment
N <- rnorm(g,ms[g],sds[g])

N <- ts(N, start=1964)

# Transform
test <- ts.detr.step(x=N, g=g, plot=TRUE)
test

test2 <- ts.detr.step(x=N, g=g, fun="median", plot=TRUE)
test2

N[c(5,18)] <- NA
test <- ts.detr.step(x=N, g=g, plot=TRUE)
test


# Example 2
ngr <- 4 # Number of groups
nrep <- c(10,20,5,10) # Number of data points in each group
g <- rep(1:ngr, nrep) # Grouping variable
ms <- c(3,6,1,4) # Means for each segment
sds <- c(0.8,0.8,0.8,0.8) # Standard deviations for each segment
N <- rnorm(g,ms[g],sds[g])

# Transform
test <- ts.detr.step(x=N, g=g, plot=TRUE)
test

test <- ts.detr.step(x=N, g=g, k=c(3, 6, 1, 4), base=2, plot=TRUE)

# Leave base out
test <- ts.detr.step(x=N, g=g, k=c(3, 6, 1, 4), plot=TRUE)

#########################################################################

