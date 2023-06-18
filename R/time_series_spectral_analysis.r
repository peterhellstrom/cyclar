# ts.spec
# Plots spectral analysis of a time series.
# Also return smoothed periodogram and spectral density of an AR-model.
# Output is invisible if not assigned.
# Essentially a wrapper for plot.ts, spec.pgram, spec.ar, and cp.gram.
# x = a time series
# spans = vector of odd integers giving the widths of modified Daniell smoothers to be used to smooth the periodogram
# see ?spec.pgram
# sig.spec = significance level for smoothed periodogram
# plot
#' @export
ts.spec <- function(x, spans=c(3,3), sig.spec=0.1,
				main=paste("Series: ", deparse(substitute(x)), sep=""), plot=TRUE) {

	if (class(x) != "ts") stop("Input must be a time series of class ts")

	if (plot == TRUE) {
		par(mfrow=c(2,2))
			# plot time series
			plot.ts(x); title(main)
			# nonparametric spectral estimate; also see spectrum()
			spec.pgram(x, spans=spans, log="no", plot=TRUE)
			abline(h=sig.spec, col=2, lty=2)
			# parametric spectral estimate
			spec.ar(x, log="no", plot=TRUE)
			# cumulative periodogram
			cpgram(x)
		par(mfrow=c(1,1))
	}

	# Create output
	out <- list(
		x = x,
		spec.pgram = spec.pgram(x, spans=spans, log="no", plot=FALSE),
		spec.ar = spec.ar(x, log="no", plot=FALSE)
	)

	invisible(out)

}

# ts.spec.peaks
# Finds the peak frequencies in the periodogram,
# x is an object fitted with spec or spec.pgram
# span = width of window (number of data points) where the peak is sought.

# Dependencies: peaks
#' @export
ts.spec.peaks <- function(x, span=9, plot=TRUE) {

	if (class(x) != "spec") stop("Object 'x' must be of class 'spec'")
	inds <- peaks(x$spec, span=span)
	z <- cbind(freq=x$freq[inds],period=1/x$freq[inds],spec=x$spec[inds])
	z <- z[order(z[,3],decreasing=TRUE), ]

	if (plot) {
		plot(x)
		abline(v=z[1,"freq"], col=2, lty=2)
		points(z[1,1],z[1,3],col=2,pch=16,cex=1.3)
	}
	z
}

# ts.spec.ar2
# Evaluates the spectral density of an AR(2)-process.
# Adapted from Shumway & Stoffer's function spec.arma (issued with the 3rd ed of their book).
# Their formula is more general, since it includes MA and also higher-order AR components.
# This function can only handle AR(2)-processes.
# The location of the peak (frequency) of the spectral density is calculated according
# to formula 13.5.9 in Cryer & Chan (2008). The original reference is Jenkins & Watts (1968), p. 229.
#' @export
ts.spec.ar2 <- function(a1, a2, sd=1, n.freq=1001, plot=FALSE, ...) {
	ar <- c(1+a1,a2)
	var.noise <- sd^2
	ar.order <- length(ar)
	freq <- seq.int(0, 0.5, length.out = n.freq)
            cs.ar <- outer(freq, 1:ar.order, function(x, y) cos(2 *
                pi * x * y)) %*% ar
            sn.ar <- outer(freq, 1:ar.order, function(x, y) sin(2 *
                pi * x * y)) %*% ar
    spec <- var.noise * (1 / ((1 - cs.ar)^2 + sn.ar^2))

	f <- acos(-((1+a1)*(1-a2))/(4*a2)) / (2*pi)
	period <- 1 / f
	out <- list(freq=freq,spec=spec,f=f,period=period)
	class(out) <- "spec"

	if (plot) {
		plot(out, ci=0, xlab="Frequency",ylab="Spectral density",font.lab=2,main="Spectral density of an AR(2)-process", ...)
		abline(v=f, col=2, lty=2)
	}

	invisible(out)
}

#' @export
ts.spec.ar22 <- function(phi1, phi2, sd=1, n.freq=1001, plot=FALSE, ...) {
	ar <- c(phi,phi2)
	var.noise <- sd^2
	ar.order <- length(ar)
	freq <- seq.int(0, 0.5, length.out = n.freq)
            cs.ar <- outer(freq, 1:ar.order, function(x, y) cos(2 *
                pi * x * y)) %*% ar
            sn.ar <- outer(freq, 1:ar.order, function(x, y) sin(2 *
                pi * x * y)) %*% ar
    spec <- var.noise * (1 / ((1 - cs.ar)^2 + sn.ar^2))

	f <- acos(-(phi1*(1-phi2))/(4*phi2)) / (2*pi)
	period <- 1 / f
	out <- list(freq=freq,spec=spec,f=f,period=period)
	class(out) <- "spec"

	if (plot) {
		plot(out, ci=0, xlab="Frequency",ylab="Spectral density",font.lab=2,main="Spectral density of an AR(2)-process", ...)
		abline(v=f, col=2, lty=2)
	}

	invisible(out)
}
