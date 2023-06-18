library(cyclar)

p <- ar2.parms(k=4, v=4, output=list)
ar2.plot()
with(p, points(1+a1,a2,pch=16,col=2,cex=1.3))
with(p, ar2.period(a1=a1, a2=a2))

x <- with(p, ar2.sim(n=50, mu=0, a1=a1, a2=a2, sd=sqrt(0.2)))
ts.summary(x=x, plot=TRUE)
ts.spec(x=x)

################################################################################
# Wavelet analysis
library(sowas)

nvoice <- 125
noctave <- 2
nreal <- 100
nvoice*noctave+1 # scales

test <- wsp(x, s0=2, nvoice=nvoice, noctave=noctave, nreal=nreal, units="years", arealsiglevel=0)
test <- wsp(exp(x), s0=2, nvoice=nvoice, noctave=noctave, nreal=nreal, units="years", arealsiglevel=0)

data(lynx)
ts.summary(log(lynx))
test <- wsp(log(lynx), s0=4, nvoice=nvoice, noctave=2, nreal=nreal, units="years", arealsiglevel=0)
str(test)

#####
x <- dget("C:/WORK/ANALYSER/-= R Development =-/cyclar/Examples/10 - Wavelets/wavelet_example_data")
plot(x, main="Simulated series with two distinct regimes", bty="l", xlab="log(x)", font.lab=2, las=1)
abline(v=60, col=2, lty=2)

ts.spec(x) # Periodogram
spec.pgram(x, spans=c(3,3), log="no", plot=TRUE)

marks <- ar2.period(a1=0.2, a2=-.9, plot=TRUE)[,"period"]

test1 <- wsp(x, s0=2, nvoice=200, noctave=3, nreal=125, units="years", arealsiglevel=0, markt=60, marks=marks, labtext="nv = 200")
test1a <- wsp(x, s0=2, nvoice=200, noctave=3, nreal=125, units="years", arealsiglevel=0, markt=60, marks=marks, labtext="nv = 200")

test2 <- wsp(x, s0=2, nvoice=100, noctave=3, nreal=125, units="years", arealsiglevel=0.9, markt=60, marks=marks, labtext="nv = 100")


# Wavelet coherency
co <- dget("C:/WORK/ANALYSER/-= R Development =-/cyclar/Examples/10 - Wavelets/wavelet_coherency_example")
plot(co)
cor(co)

test5a <- wsp(co[,1], s0=2, nvoice=100, noctave=2, nreal=125, units="years", arealsiglevel=0)
test5b <- wsp(co[,2], s0=2, nvoice=100, noctave=2, nreal=125, units="years", arealsiglevel=0)

test5 <- wco(co[,1],co[,2], tw=1.5, sw=0.5)
