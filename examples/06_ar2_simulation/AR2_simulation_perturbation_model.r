library(cyclar)
library(zoo)
library(pastecs)

# Perturbation
perturbation <- function(spec = "unif", p = 0.1, a = 0.4, b = 0.6, n = 100, ...) {
	rnd <- runif(n)
	inds <- which(rnd <= p)
	u <- rep(1, n)
	u[inds] <- rtrunc(n = length(inds), spec = spec, a, b, ...)
	u
}

perturbation()


# Distribution of perturbation, u
n <- 100000
system.time(u <- perturbation(spec = "unif", p = 0.1, a = 0.4, b = 0.6, n = n))
hist(u[which(u != 1)], col = "steelblue", breaks = 30, freq = FALSE, xlab = "u", main = "Distribution of perturbations, u")
length(u[which(u != 1)]) / n # p


n <- 100000
mu <- 0.5
sigma <- 0.07
a <- mu - 0.15
b <- mu + 0.15

system.time(u <- perturbation(spec = "norm", mean = mu, sd = sigma, p = 0.2, a = a, b = b, n = n))
hist(u[which(u != 1)],  col = "steelblue",  breaks = 30, freq = FALSE, xlab = "u", main = "Distribution of perturbations, u")
curve(dtrunc(x, spec = "norm", mean = mu, sd = sigma, a = a, b = b), add = TRUE, n = 1001)
curve(dnorm(x, mean = mu, sd = sigma), add = TRUE, n = 1001, col = 2)
range(u[which(u != 1)])
a
b
length(u[which(u != 1)]) / n # p

# Distribution of probability of perturbation, p
n <- 1000
nrepl <- 10000
xpt <- sapply(1:nrepl, function(i) {
	mu.sim <- perturbation(n = n)
	length(which(mu.sim != 1)) / n # Should give approx p = 0.1
})

hist(xpt, col = "steelblue", breaks = 30, freq = FALSE)
lines(density(xpt))
abline(v = mean(xpt), col = 2, lwd = 2, lty = 2)

###
# Delayed Ricker model
r <- 0.9
a1 <- 0.005
a1 <- 0.015
# a1 <- -0.4
a2 <- -0.1

r <- 0.4
a1 <- 0.12
a2 <- -0.2

xbar <- -r / (a1 + a2) # Equilibrium population level
xbar

exp(r + (a1+a2)*xbar)

n <- 50
x <- numeric(n+1)

if (xbar < 0) x[1:2] <- runif(2, 1.25*xbar, 0.75*xbar)
if (xbar >= 0) x[1:2] <- runif(2, 0.75*xbar, 1.25*xbar)

# Without perturbation
for (t in 2:n) {
	x[t+1] <- x[t] * exp(r + a1*x[t] + a2*x[t-1])
}

x1 <- numeric(n+1)
x1[1:2] <- x[1:2]

for (t in 2:n) {
	x1[t+1] <- x1[t] * exp(r + a1*x1[t] + a2*x1[t-1])
}

plot(zoo(x))
ts.summary(x, plot = TRUE)
ts.spec(as.ts(x))

# With perturbation
u <- perturbation(n = n)

for (t in 2:n) {
	x[t+1] <- u[t] * x[t] * exp(r + a1*x[t] + a2*x[t-1])
}


par(mfrow = c(1,2))
	plot(x,type = "l")
	abline(h = xbar,col = 2,lty = 2)
	acf(x)
par(mfrow = c(1,1))


plot(scale(x),type = "l")
abline(h = 0, col = 2, lty = 2)

hist(x, col = "steelblue", breaks = 30)

####
# Nonlinear Royama model
a1 <- -0.0732
a2 <- -0.8819

n <- 50
x <- numeric(n + 1)
init <- runif(2)
x[1:2] <- exp(init)
for (t in 2:n) {
	x[t+1] <- x[t] * exp(1 - 1 / ( (x[t]^a1) * (x[t-1]^a2)))
}

x1 <- numeric(n + 1)
x1[1:2] <- init
for (t in 2:n) {
	x1[t+1] <- x1[t] + (1 - exp(-a1*x1[t]-a2*x1[t-1]))
}


x
exp(x1)



u <- perturbation(n=n,p=0.2)

for (t in 2:n) {
	x[t+1] <- u[t] * x[t] * exp(1 - 1 / ( (x[t]^a1) * (x[t-1]^a2)))
}

par(mfrow=c(1,2))
	plot(scale(x), type="l")
	acf(x)
par(mfrow=c(1,1))

plot(scale(x),type="l")
abline(h=0, col=2, lty=2)

hist(x, col="steelblue", breaks=30)

# Include time-varying p and mu
# Initial values?

# Stability analysis
r <- 0.9
a1 <- 0.01
a2 <- -0.1
xbar <- - r/(a1+a2)

mat <- matrix(c(a1*xbar + 1, a2*xbar, 1, 0), ncol=2, nrow=2, byrow=TRUE)
A <- eigen(mat)
a <- Re(A$values) # Real part
b <- Im(A$values) # Imaginary part
R <- Mod(A$values) # Modulus
phi <- acos(a[1]/R[1])
per <- 2*pi/phi # Period

unit.circle <- function() {
	#dev.new(width=6.375,height=7)
	plot(x=c(-1.1,1.1), y=c(-1.1,1.1), type="n", xlab="Re", ylab="Im", font.lab=2, las=1,main="Unit circle in complex plane")
	abline(v=0, lty=2)
	abline(h=0, lty=2)
	curve(sqrt(1-x^2), from=-1, to=1, add=T)
	curve(-sqrt(1-x^2), from=-1, to=1, add=T)
}

unit.circle()
segments(x0=0,y0=0,x1=a[1],y1=b[1], col=2, lty=2)
points(a,b,col=2,pch=16)

a1 <- 1 + -0.5
a2 <- -0.7


ar2.period(a1=a1-1,a2=a2,print.all=TRUE,plot=TRUE)
a <- a1/2
b <- (1/2) * sqrt(-a1^2-4*a2)
R <- sqrt(a^2 + b^2)
a;b;R

unit.circle()
