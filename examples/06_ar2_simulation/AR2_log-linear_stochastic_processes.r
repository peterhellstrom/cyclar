################################################################################
library(cyclar)
################################################################################

k <- 3.5
v <- 2
p <- ar2.parms(k=k, v=v)
mu <- 2
a1 <- p$a1
a2 <- p$a2

n <- 1000
nrepl <- 10000
z <- replicate(ar2.sim2(n=n, mu=mu, a1=a1, a2=a2),n=nrepl)
dim(z) # Time = rows, Series = columns

# Series means
dsz <- density(colMeans(z))
# Time means, expectation of the stochastic process at time t
dtz <- density(rowMeans(z))
xlim <- range(c(dsz$x,dtz$x))
ylim <- range(c(dsz$y,dtz$y))

plot(x=xlim,y=ylim,type="n",ylab="density",xlab="mean")
lines(dsz, col=1)
lines(dtz, col=2)

abline(v=mean(colMeans(z)), col=1, lty=2)
abline(v=mean(rowMeans(z)), col=2, lty=2)
abline(v=mu, col=2, lty=2, lwd=2)

legend("topleft", c("Series means", "Time means"), lty=1, col=1:2, bty="n")

summary(colMeans(z))
summary(rowMeans(z))
