# ar2.sim
# Function that simulates an AR2-process
# Mainly a wrapper for arima.sim, with extra plot options.
# This function is rather basic, and more options are available
# in arima.sim (such as better control of innovations and distribution of innovations).
# I have created another function called ar2.sim2 that is very similar to arima.sim.

# mu = long-term mean of time series
# AR(2)-parameters a1 & a2 refer to Royama's notation.
# a1 is direct density-dependence, a2 is delayed density-dependence
# arima.sim uses the following input form:
# arima.sim(model=list(ar=c(1+a1,a2)),n=n) + mu

# AR(2)-parameters can not be entered as vectors.
# Use lapply to "loop" over vectors with AR(2)-parameters instead.
# Attributes p (simulation input parameters) and a0 (intercept) are added
# to the output series.

# sd = used only to generate innovations? Effect on time-series? Test this!

# Output is a time series x, on log-scale.
#' @export
ar2.sim <- function(n=100, mu=0, a1, a2, sd=1, n.start=1000, plot=FALSE) {

  x <- arima.sim(model=list(ar=c((1+a1),a2)), n=n, sd=sd, n.start=n.start) + mu

	# PLOT
	if (plot==TRUE) {
		par(mfrow=c(2,2))

		# 1) Triangle + coefficients
		plot(c(-2,2), c(-1,1), type="n",
		  xlab=expression(paste("1+ ",beta,"1")), ylab=expression(paste(beta,"2")),
		  main="AR2 parameters",
		  font.lab=2, las=1)
		points(c(-2,0,2,-2), c(-1,1,-1,-1), type="l")
		lines(c(0,0), c(-1,1))
		curve(-0.25*x^2, from=-2, to=2, add=T)
		points((1+a1), a2, col="red", cex=1.4, pch=16)
		# 2) Simulated time series
		plot.ts(x, main=paste("Simulated time series"), xlab="Time", ylab="log(x)", font.lab=2, las=1)
		# 3) Auto-correlation
		acf.c <- acf(x, type="correlation", main="Autocorrelation", font.lab=2, las=1)
		# 4) Partial rate autocorrelation
		acf.p <- acf(x, type="partial", main="Partial autocorrelation", font.lab=2, las=1)

		par(mfrow=c(1,1))
	}

	# Print data:
	a0 <- mu*(-(a1+a2))
	p <- c(mu=as.numeric(mu), a0=as.numeric(a0), a1=as.numeric(a1), a2=as.numeric(a2), sd=sd)
	attr(x, "p") <- p
	x

}

# ar2.sim2
# Custom function for generation of AR(2)-data.
# Syntax & arguments is similar to the main arima.sim function.
# See ?arima,sim for info on arguments rand.gen, innov & n.start
# This function uses two methods. method=c("filter","recursive")
# method="loop" uses an recursive update function (for-loop), whereas
# method="filter" filters the innovations linearly with the parameters (1+a1) & a2.

# mu is the long-term mean of series x
# a1 is the strength of direct density-dependence
# a2 is the strength of delayed density-dependence

# Output is a time series x, on log-scale.
#' @export
ar2.sim2 <- function(mu=0, a1, a2, n=100, rand.gen=rnorm, innov=rand.gen(n, ...),
	n.start=NA, start.innov = rand.gen(n.start, ...), method=c("filter","recursive"), ...) {

	method <- match.arg(method)
	if (!method %in% c("filter","recursive")) stop("Method must be either 'filter' or 'recursive'")

	# If n.start is NA, calculate a reasonable number of start values
	# Code borrowed from arima.sim
	if (is.na(n.start)) {
		ar <- c(1+a1, a2)
		p <- length(ar)
		minroots <- min(Mod(polyroot(c(1, -ar))))
		if (minroots <= 1) stop("'ar' part of model is not stationary")
		n.start <- p + ifelse(p > 0, ceiling(6/log(minroots)), 0)
		if (n.start < p) stop("burn-in 'n.start' must be as long as 'ar'")
	}

	if (method=="recursive") {
		eps <- c(start.innov[1L:n.start], innov[1L:n])
		x <- numeric(length(eps))
		x[1] <- eps[1] + mu
		x[2] <- mu*(-a1) + (1+a1)*(x[1]) + eps[2]
		a0 <- mu*(-(a1+a2))
		for (t in 3:(n+n.start)) x[t] <- a0 + (1+a1)*x[t-1] + a2*x[t-2] + eps[t]
		if (n.start > 0) x <- x[-(1L:n.start)]
	}

	if (method=="filter") {
		x <- ts(c(start.innov[1L:n.start], innov[1L:n]), start = 1 - n.start)
		x <- filter(x, filter=c((1+a1),a2), method="recursive") + mu
		if (n.start > 0) x <- x[-(1L:n.start)]
	}

	x <- as.ts(x)
	a0 <- mu*(-(a1+a2))
	p <- c(mu=as.numeric(mu), a0=as.numeric(a0), a1=as.numeric(a1), a2=as.numeric(a2))
	attr(x, "p") <- p
	x
}

# use pars, with etc.
# n.start might not be working for all models, perhaps use only a numerical value input? Default perhaps n.start = 1000
#' @export
ar2.sim.model <- function(mu=0, a1, a2, n=100, rand.gen=rnorm, innov=rand.gen(n, ...),
	n.start=NA, start.innov = rand.gen(n.start, ...), model="linear", ...) {

	# If n.start is NA, calculate a reasonable number of start values
	# Code borrowed from arima.sim
	if (is.na(n.start)) {
		ar <- c(1+a1,a2)
		p <- length(ar)
		minroots <- min(Mod(polyroot(c(1, -ar))))
		if (minroots <= 1) stop("'ar' part of model is not stationary")
		n.start <- p + ifelse(p > 0, ceiling(6/log(minroots)), 0)
		if (n.start < p) stop("burn-in 'n.start' must be as long as 'ar'")
	}

	if (model=="linear") {
		eps <- c(start.innov[1L:n.start], innov[1L:n])
		x <- numeric(length(eps))
		x[1] <- eps[1] + mu
		x[2] <- mu*(-a1) + (1+a1)*x[1] + eps[2]
		a0 <- mu*(-(a1+a2))
		for (t in 3:(n+n.start)) x[t] <- a0 + (1+a1)*x[t-1] + a2*x[t-2] + eps[t]
		if (n.start > 0) x <- x[-(1L:n.start)]

		x <- as.ts(x)
		a0 <- mu*(-(a1+a2))
		p <- c(mu=as.numeric(mu), a0=as.numeric(a0), a1=as.numeric(a1), a2=as.numeric(a2))
		attr(x, "p") <- p
		x
	}

	if (model=="nonlinear") {


	}

	if (model=="competition") {
		# input
		# variables : x (focal species), y (competitor)
		# parameters: a1, a2, b1, b2
		# noise: eps

		# Deviations from their means?
		# Intercept term?
		# Add mu AFTER simulation

		eps <- c(start.innov[1L:n.start], innov[1L:n])
		x <- numeric(length(eps))

		x[1] <- eps[1] + mu
		x[2] <- mu*(-a1) + (1+a1)*x[1] + b1*y[1]+ eps[2]
		a0 <- mu*(-(a1+a2))

		for (t in 3:(n+n.start)) x[t] <- a0 + (1+a1)*x[t-1] + a2*x[t-2] + b1*y[t-1] + b2*y[t-2] + eps[t]
		if (n.start > 0) x <- x[-(1L:n.start)]

		x <- as.ts(x)
		a0 <- mu*(-(a1+a2))
		p <- c(mu=as.numeric(mu), a0=as.numeric(a0), a1=as.numeric(a1), a2=as.numeric(a2))
		attr(x, "p") <- p
		x
	}

}


# ar2.sim.npatch
# Simulate correlated time series (based on correlated noise, w[t] = Moran effect).
# n = number of data points
# npatch = number of patches
# rho = correlation coefficient (beween -1 & 1)
# sd = noise of the AR(2) - process.
# nstart = number of burn-in steps
# mu = mean of time series
# a1, a2 = AR(2)-parameters

# The arguments mu and sd can either be single values or vectors with the same
# length as the number of simulated patches (controlled by npatch).
# If mu and sd are vectors of length 1, all patches have the same sd and mean.
# For cases where habitats or patches are correlated and have the same density-dependent structure, but might have different
# sd and mean (for instance along productivity gradients), mu and sd can be entered as vectors.

# ToDo: COMPARE PERFORMANCE OF THIS FUNCTION compared to ar2.sim()

# NOTE: at present, only mu & sd are allowed to change between patches.
# It would however be easy to update the function (if needed) so that a1 & a2 can change as well.
#' @export
ar2.sim.npatch <- function(n, npatch, rho, mu, a1, a2, sd=1, n.start=100, plot=TRUE) {

	# Covariance matrix
	Sigma <- covmat(rho=rho, sigma=sd, n=npatch)

	# Generate the error part
	# w is a matrix with random white noise, err ~ N(mu,Sigma)
	# i.e. w is drawn from a multinormal distribution with mean 0 and
	# covariance matrix Sigma. The function mvrnorm in MASS is used here.
	w <- mvrnorm(n=n + n.start, mu=rep(0,npatch), Sigma=Sigma)
	# check: cor(w) (should be close to rho for non-diagonal elements)

	# Filter the white noise series
	if (length(mu) == 1) x <- filter(w, filter=c((1+a1),a2), method="recursive") + mu
	if (length(mu) > 1) {
		if (length(mu) != npatch) stop("Number of patches and length of 'mu' vector must be equal")
		x <- filter(w, filter=c((1+a1),a2), method="recursive")
		x <- t(t(x) + mu)
	}

	# Remove burn-in sequence
	x <- window(x,n.start+1,n.start+n)
	w <- window(w,n.start+1,n.start+n)
	colnames(w) <- colnames(x)

	if (plot == TRUE) {
		op <- par(no.readonly = TRUE)
		matplot(x, type="l", xlab="Time", ylab="x, AR(2) on log-scale", font.lab=2, las=1, bty="l", lty=1,
			main=paste("Simulated time series, rho =",rho))
		abline(h=mu, col=1:length(mu), lty=2)
		par(op)
	}

	out <- list(x = as.ts(x), w=ts(w), Sigma=Sigma,
		p=c(npatch=as.numeric(npatch), rho=as.numeric(rho), mu=as.numeric(mu), a1=as.numeric(a1), a2=as.numeric(a2), sd=as.numeric(sd))
	)

	out
}
