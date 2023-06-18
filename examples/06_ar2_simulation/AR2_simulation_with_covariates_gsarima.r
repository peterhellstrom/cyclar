# garsim function in gsarima:

library(gsarima)
n <- 100
ar.coef <- c(1.5060442, -0.67043)
beta.x <- 1.5
intercept <- 2
# Covariate
x <- rnorm(n)
# Design matrix
X <- matrix(c(rep(intercept, n+length(ar.coef)), rep(0, length(ar.coef)), x), ncol=2)
dim(X)

xv <- X %*% matrix(c(intercept, beta.x),ncol=1)

y.sim <- garsim(n=n+length(ar.coef), phi=ar.coef, beta=c(1, beta.x), sd=sqrt(1))
y.sim <- y.sim[(1+length(ar.coef)):(n+length(ar.coef))]

plot(ts(y.sim))