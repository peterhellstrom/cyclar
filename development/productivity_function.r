pr.fn <- function(x,a,b,theta) a * (1 - (1/(1+(x/b)^theta)))
pr.fn2 <- function(x,a,b,theta) a * x^theta / (b^theta + x^theta)

b <- seq(1,12,2)
#theta <- seq(1,10,1)
theta <- c(1,5,10)

# Vary threshold, b
plot(x=c(0,15),y=c(0,1.8), type="n", xlab="Prey density", ylab="Number of young produced", font.lab=2, bty="l")
for(i in b) curve(pr.fn(x,a=1.75,b=i,theta=3), from=0, to=15, col=i, add=TRUE)
legend("bottomright", legend=b, col=b, lty=1, bty="n", title="theta")
for(i in b) curve(pr.fn2(x,a=1.75,b=i,theta=3), from=0, to=15, col=2, add=TRUE, lty=2)

# Vary theta, steepness
plot(x=c(0,15),y=c(0,1.8), type="n", xlab="Prey density", ylab="Number of young produced", font.lab=2, bty="l")
for(i in theta) curve(pr.fn(x,a=1.75,b=6,theta=i), from=0, to=15, col=i, add=TRUE)
legend("bottomright", legend=theta, col=theta, lty=1, bty="n", title="theta")
#for(i in theta) curve(pr.fn2(x,a=1.75,b=6,theta=i), from=0, to=15, col=2, add=TRUE, lty=2)

# Vary both b and theta, double for-loop
plot(x=c(0,15),y=c(0,1.8), type="n", xlab="Prey density", ylab="Number of young produced", font.lab=2, bty="l")
for(i in b) {
	for (j in theta) {
		curve(pr.fn(x,a=1.75,b=i,theta=j), from=0, to=15, col=j, add=TRUE)
	}
}

# Same as above, both use expand.grid instead of double for-loop
zz <- expand.grid(b=b, theta=theta)
plot(x=c(0,15), y=c(0,1.8), type="n", xlab="Prey density", ylab="Number of young produced", font.lab=2, bty="l")
for(i in 1:nrow(zz)) curve(pr.fn(x,a=1.75,b=zz$b[i],theta=zz$theta[i]), from=0, to=15, col=zz$theta[i], add=TRUE)
legend("bottomright", legend=unique(zz$theta), col=unique(zz$theta), lty=1, bty="n", title="theta")
	
colors()[1:20]
col2rgb(colors()[1:20])

palette()
palette()[1:3]

colors()