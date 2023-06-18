# UPDATE - NOT CORRECT AT THE MOMENT!

# Simulate time-series with time-varying parameters
# a0, a1, a2 are vectors of length n, with values for time-varying parameters
# ar.var = sigma^2
# n = desired length of sequence, n.start = length of burnin sequence

# Q: what happens if the noise (innovations) w are nonstationary or time-varying?
# Two different ways of losing cycles:
# time-varying AR(2)-coefficients - changes in the density dependent structure
# OR changes in the density-independent part (environmental "perturbations")

#' @export
ar2.sim.tv <- function(n, a0, a1, a2, sd, innov=TRUE, vals=NULL, n.start=500, plot=TRUE) {

	ns <- c(length(a0), length(a1), length(a2))
	if (length(unique(ns)) > 1) stop("Vectors of unequal length")
	if (unique(ns) != n) stop("Length of input vectors must equal n")

	op <- options(warn = (-1)) # Suppress warnings (zero length arrows)

	# Generate burn-in sequence
	# Iterate the system n.start steps with initial value of a0, a1 & a2.
	n.start <- n.start
	wb <- rnorm(n.start,0,sd) # Error
	xb <- filter(wb, filter=c((1+a1)[1], a2[1]), method="recursive") + a0[1]
	init.vals <- tail(xb,2) # Store final values as initial values for the recursive ar-filtering.

	if (is.logical(innov) == TRUE) {
			# Generate sequence to use for "real" simulation
			w <- rnorm(n, 0, sd)
	}

	if (is.numeric(innov) == TRUE) {
			w <- innov
			if (length(w) != n) stop("w must match n")
	}

	x <- numeric(length(n))
	mu.v <- numeric(length(n))

	# For.loop
	for (i in 1:n) {
			a0t <- a0[i]
			a1t <- a1[i]
			a2t <- a2[i]

			z <- sum(a1t,a2t)
			mu.v[i] <- mu <- a0t*(1-z)

			if (i == 1) x[i] <- mu + a1t*init.vals[2] + a2t*init.vals[1] + w[i]
			if (i == 2) x[i] <- mu + a1t*x[1] + a2t*init.vals[2] + w[i]
			if (i > 2) x[i] <- mu + a1t*x[i-1] + a2t*x[i-2] + w[i]
	}

	exp.x <- exp(x)

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

		ar2.arrows(a1, a2, col="red", length=0.05)

		# 2) Simulated time series
		plot.ts(x, main=paste("Simulated time series"), xlab="Time", ylab="log(x)", font.lab=2, las=1)
		if (is.null(vals) == FALSE) abline(v=cumsum(vals[,1]), lty=2, col=2)
		# 3) Auto-correlation
		acf.c <- acf(x, type="correlation", main="Autocorrelation", font.lab=2, las=1)
		# 4) Partial rate autocorrelation
		acf.p <- acf(x, type="partial", main="Partial autocorrelation", font.lab=2, las=1)

		par(mfrow=c(1,1))

	}

	# Print list with data:
	out <- list(
		"n" = n,
		"p" = cbind(a0=as.numeric(a0), a1=as.numeric(a1), a2=as.numeric(a2)),
		"ts.x" = ts(x),
		"exp.ts.x" = ts(exp.x))

	op <- options(warn = (0))
	out

}


# UPDATE THIS FUNCTION WITH THE FOLLOWING FEATURES:
# allow for differences in mean and sd
#' @export
ar2.sim.tv.gen <- function(npoints, coord.mat, route, vals, gam=TRUE) {

	if (npoints != nrow(coord.mat)) stop("npoints must equal number of rows in coord.mat")
	if ((length(route) - 1) != nrow(vals)) stop("Number of steps must equal number of rows in matrix vals")

	nsteps <- length(route) - 1
	nstops <- length(route)

	arrow.mat <- matrix(c(
		x0 = c(coord.mat[route[-nstops],1]),
		y0 = c(coord.mat[route[-nstops],2]),
		x1 = c(coord.mat[route[2:nstops],1]),
		y1 = c(coord.mat[route[2:nstops],2])
		), nrow=4, ncol=nsteps, byrow=TRUE)

	ar2.plot()
		for (i in 1:npoints) points(coord.mat[i,1], coord.mat[i,2], col=2, cex=1.3, pch=16)
		for (i in 1:npoints) text(coord.mat[i,1]+0.1, coord.mat[i,2], labels=rownames(coord.mat)[i], font=2)

		arrows(
			x0 = arrow.mat[1,], y0 = arrow.mat[2,], x1 = arrow.mat[3,], y1 = arrow.mat[4,],
				col=2, length=0.15, angle=20, lwd=2, code=2)

	val.mat <- matrix(nrow=nsteps,ncol=3,dimnames=list(
		paste(route[-length(route)],"-",route[-1],sep=""),
		c("duration","mean","sd")))


	val.mat[1:nsteps,] <- vals

	time <- 1:sum(val.mat[,"duration"])
	trans.mat <- coord.mat[route,]

	a1 <- lapply(1:nsteps, function(i) seq(trans.mat[i,1],trans.mat[i+1,1],length.out=val.mat[i]))
	a1 <- unlist(a1)

	a2 <- lapply(1:nsteps, function(i) seq(trans.mat[i,2],trans.mat[i+1,2],length.out=val.mat[i]))
	a2 <- unlist(a2)

	if (gam == TRUE) {
		gam.a1 <- gam(a1 ~ s(time))
		gam.a2 <- gam(a2 ~ s(time))

		par(mfrow=c(1,2))
			plot(time, a1, type="l", pch=16, font.lab=2, las=1, main="Direct density-dependence")
				points(time,predict(gam.a1),type="l",col=2)
			plot(time, a2, type="l",pch=16, font.lab=2, las=1, main="Delayed density-dependence")
				points(time, predict(gam.a2),type="l",col=2)
		par(mfrow=c(1,1))

		a1 <- predict(gam.a1)
		a2 <- predict(gam.a2)
		}

	if (gam == FALSE) {
		par(mfrow=c(1,2))
			plot(time, a1, type="l", pch=16, font.lab=2, las=1, main="Direct density-dependence")
			plot(time, a2, type="l",pch=16, font.lab=2, las=1, main="Delayed density-dependence")
		par(mfrow=c(1,1))
		}

	ntime <- sum(val.mat[,"duration"])

	list(a1 = ts(a1), a2 = ts(a2), ntime = ntime)
}
