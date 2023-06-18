# ts.summary
# Extract basic descriptive statistics of a time series
# Mean, sd, quantiles (25%, 50%, 75%), n (length of x), dominant acf & pacf

# x = time series
# base = c("exp","log","raw"). Default is "raw", when the entered series x is used.
# option "log" uses the log(x)-series, and "exp" takes exponent of x for analyses.
# quantile = vector with quantiles (default is c(.25, .50, .75)).
# plot = TRUE/FALSE (TRUE returns 2*2 plot with time series, smoothed periodogram, acf & pacf.
# main = optional vector for time series name to be plotted (defaults to Series: [name of x]).

# Dependencies: ts.diag.acf, {zoo}

# Known issues:
# ar2.sim & ar2.sim2 adds extra attributes to the time series.
# The function as.zoo() can not handle these attributes (positions 3 & 4), but it works with zoo().
#' @export
ts.summary <- function(x, quantile = c(0.25,0.5,0.75), log=FALSE, base=exp(1), constant=0,
			plot=FALSE, main=NULL) {

	if (any(is.na(x))) cat("Time series 'x' contains missing values","\n")

	n <- length(x)

	if (is.null(main)) main.str <- paste("Series: ", deparse(substitute(x)), sep="")

	if (!log && constant != 0) x <- x + constant

	if (log) {
		if (any(x + constant <= 0)) stop("Time series 'x' contains zeros or negative values")
		x <- log(x + constant, base=base)
	}

	m <- mean(x, na.rm=T)
	s <- sd(x, na.rm=T)

	out <- c(
		mean = m,
		sd = s,
		cv = s/m,
		min = min(x, na.rm=T),
		max = max(x, na.rm=T),
		quantile(x, quantile, na.rm=T),
		n = n,
		ts.diag.acf(x, f="acf", plot=FALSE, max.only=TRUE),
		#crit.acf = qnorm((1 + 0.95)/2) / sqrt(n),
		ts.diag.acf(x, f="pacf", plot=FALSE, max.only=TRUE)
    )

	if (plot) {
		par(mfrow=c(2,2))

		# 1) Time series
		plot(zoo(x), main=main.str, xlab="Time", ylab="x", font.lab=2, las=1)

		# 2) Spectrogram, TRUE raw periodogram
		spec.pgram(x, taper=0, fast=FALSE, detrend=FALSE, log="no", main=main)

		# 3) Auto-correlation
		acf.c <- acf(x, type="correlation", main="Autocorrelation", font.lab=2, las=1)

		# 4) Partial rate autocorrelation
		acf.p <- acf(x, type="partial", main="Partial autocorrelation", font.lab=2, las=1)

		par(mfrow=c(1,1))
	}

	round(out,4)
}

# ts.cobweb
# Crete a Ricker (a.k.a. cobweb) plot.
# Inputs are a (time) series x and an optional argument for K = carrying capacity
# if K is supplied as a numerical value, a horizontal value is added at K.

# Dependencies: {zoo}
#' @export
ts.cobweb <- function (x, K=NULL) {

	n <- length(x)

	par(mfrow=c(1,2))
	# Time series plot
	plot(zoo(x), lty=2, xlab="Time", ylab="Population size", las=1, font.lab=2, main="Time series", bty="l")
	points(zoo(x), pch=16)
	if (is.numeric(K)) abline(h=K, col=2, lty=2)
	# Ricker plot
	plot(x[-n], x[-1], xlab=expression(x[" t"]), ylab=expression(x[" t+1"]), font.lab=2,las=1,
		main="Ricker (cobweb) plot", pch=16, bty="l")
	abline(0,1)
	lines(
		x=rep(x[-n], each=2),
		y=c(0,rep(x[c(-1,-n)], each=2),x[n]), col=2, pch=16
	)
	par(mfrow=c(1,1))
}
