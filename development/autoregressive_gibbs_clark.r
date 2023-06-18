# p. 92 - 
# Clark, Lab Manual

n <- 50
x <- cbind(rep(1,n), runif(n,0,20))
b0 <- 0.5
b1 <- 2
sig <- 1
rho <- -0.5
beta <- matrix(c(b0,b1),2,1)

inv.rmat <- function(rho) {
	rinv <- diag((1+rho^2),n,n)
	rinv[1,1] <- 1
	rinv[n,n] <- 1
	rinv[row(rinv) == (col(rinv)-1)] <- -rho
	rinv[row(rinv) == (col(rinv)+1)] <- -rho
	return(rinv)
}

solve(inv.rmat(rho))
sinv <- inv.rmat(rho)/sig
smat <- sig*solve(inv.rmat(rho))

y <- matrix(0,n,1)
error <- rnorm(1,0,sqrt(sig))
y[1,1] <- t(x[1,]) %*% beta + error
for (t in 2:n) {
	error <- rnorm(1,rho*error,sqrt(sig)) # AR(1) process
	y[t,1] <- t(x[t,]) %*% beta + error
}

library(MASS)
?mvrnorm
y <- t(mvrnorm(1, (x %*% beta), smat))

bprior <- c(0,0) # priors for regression parameters
vinvert <- solve(diag(1000,2)) # inverse covariance for beta's
ppart <- vinvert %*% bprior
smu <- 2
s1 <- 2
s2 <- smu*(s1-1)
