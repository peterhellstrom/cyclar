# ts.detr
# Detrend time series by fitting linear, polynomial, loess and GAM regression models.
# The residuals are used as the detrended time series.
# span, loess.degree & family are loess parameters (default settings are approximately equal to lowess defaults).
# span = the parameter (alpha) which controls the degree of smoothing, i.e. the proportion of the points
# that are used in each window.
# loess.degree = degree of local polynomial, usually 1 or 2. Default for ts.detr is 1.
# family = c("gaussian","symmetric"). See ?loess
# poly.degree = degree/order of polynomial
# add.mean = TRUE/FALSE, whether the mean of time series x should be added to the detrended series or not.
# constant = numeric value of a constant that is added to the time series (before detrending). Defaults to 0.
# log = TRUE/FALSE, option TRUE returns the natural logarithm of (x + constant)
# plot = TRUE/FALSE. If plot=TRUE, a four-panel plot is drawn. The original time series is plotted with trend lines in panel 1.
# Panels 2-4 contain the detrended series (one panel each for the linear, loess, and gam models).
# title = optional string with name of the time series, is only added to the plot.
# output = c("mts","list") output type can be either a multiple time series object or a list.

# Dependencies: {zoo}, {mgcv}

# for gam settings: see ?s
# The settings used here are:
# s(time, bs="cr", fx=FALSE, k=-1), which means a penalized cubic regression smooth term.
# The degree of smoothing is controlled automatically, unlike loess or lowess.
#' @export
ts.detr <- function(x, span=0.75, loess.degree=1, family="symmetric", poly.degree=3,
  add.mean=TRUE, constant=0, log=FALSE, plot=TRUE, title=NULL, output="mts") {

	if (class(x) != "ts") stop("Input series must be a ts-object")

	if (log==FALSE) x <- x + constant
	if (log==TRUE) x <- log(x + constant)

	t <- time(x)

	# Detrend
	fm1 <- lm(x ~ t)
	fm2 <- loess(x ~ t, span=span, degree=loess.degree, family=family)
	fm3 <- gam(x ~ s(t, bs="cr", fx=FALSE, k=-1))
	fm4 <- lm(x ~ poly(t, poly.degree, raw=TRUE))

	# "Intersect" original time series with detrended, necessary if ts contains missing values (NA)
	x.detr.lm <- add.NAs(x=x, obj=fm1)
	x.detr.loess <- add.NAs(x=x, obj=fm2)
	x.detr.gam <- add.NAs(x=x, obj=fm3)
	x.detr.poly <- add.NAs(x=x, obj=fm4)

	if (add.mean == TRUE) {
		mu <- mean(x, na.rm=TRUE)
		x.detr.lm <- x.detr.lm + mu
		x.detr.loess <- x.detr.loess + mu
		x.detr.gam <- x.detr.gam + mu
		x.detr.poly <- x.detr.poly + mu
	}

	if (output=="list") {
		out <- list(
			t = t,
			x = x,
			x.lm = x.detr.lm,
			x.loess = x.detr.loess,
			x.gam = x.detr.gam,
			x.poly = x.detr.poly
	)}

	if (output=="mts") {
		out <- ts.union(
			x = x,
			x.lm = x.detr.lm,
			x.loess = x.detr.loess,
			x.gam = x.detr.gam,
			x.poly = x.detr.poly
	)}

	# Plot
	if (plot == TRUE) {

		# Graph1, plot original series x together with trend lines
		if (is.null(title) == FALSE) main.lab <- paste("Time series, ", title, sep="")
		if (is.null(title) == TRUE) main.lab <- c("Time series")

		op <- par(bty="l")

		plot.detr(x=x, pch=16, p.col=1, p.cex=1, l.col=2, lwd=1, ylim=range(x, na.rm=TRUE), main=main.lab)
		abline(h=mean(x),col=1, lty=2)
		lines(x, lty=2)

		xv <- seq(from=start(x)[1], to=end(x)[1], length.out=2*length(x))
		lines(xv, predict(fm1, list(t=xv)), col=2, lwd=2, lty=1)
		lines(xv, predict(fm2, newdata=(t=xv)), col=4, lwd=2, lty=1)
		lines(xv, predict(fm3, list(t=xv)), col=3, lwd=2, lty=1)
		lines(xv, predict(fm4, list(t=xv)), col=5, lwd=2, lty=1)

		legend(0, max(x)*1.45, c("Linear", "LOESS", "GAM", "Polynomial"), col=c(2,4,3,5), lwd=c(2,2,2,2),
			bty="n", horiz=TRUE, xpd=TRUE, cex=0.9)

		par(op)

		# Graph 2, plot each detrended series in a separate panel
		ylims <- range(c(x.detr.lm, x.detr.loess, x.detr.gam), na.rm=TRUE)

		par(mfrow=c(2,2))

		plot.detr(x=x.detr.lm, pch=16, p.col=1, p.cex=1, l.col=2, lwd=1, ylim=ylims, main="Linear detrend")
		plot.detr(x=x.detr.loess, pch=16, p.col=1, p.cex=1, l.col=4, lwd=1, ylim=ylims, main="LOESS detrend")
		plot.detr(x=x.detr.gam, pch=16, p.col=1, p.cex=1, l.col=3, lwd=1, ylim=ylims, main="GAM detrend")
		plot.detr(x=x.detr.poly, pch=16, p.col=1, p.cex=1, l.col=5, lwd=1, ylim=ylims,
			main=paste("Polynomial detrend, order =", poly.degree))

		par(mfrow=c(1,1))
	}

	# Print output
	invisible(out)
}

### Intermediate functions ###
# "Intersect" original time series with detrended, necessary if ts contains missing values (NA)
#' @export
add.NAs <- function(x, obj) {

	if (any(is.na(x)) == TRUE) {
		x.detr <- obj$residuals
		NAs <- obj$na.action
		names(x.detr) <- (1:length(x))[-as.numeric(NAs)]
		NAs[1:length(NAs)] <- NA
		x.detr <- c(x.detr, NAs)
		x.detr <- x.detr[order(as.numeric(names(x.detr)))]
		x.detr <- ts(as.numeric(x.detr))
		tsp(x.detr) <- tsp(x)
	}

	if (all(!is.na(x)) == TRUE) {
		x.detr <- ts(obj$residuals)
		tsp(x.detr) <- tsp(x)
	}

	x.detr
}

# Plot (detrended) time series
#' @export
plot.detr <- function (x, pch=16, p.col=1, p.cex=1, l.col=1, lwd=1, lty=2, ylim=range(x, na.rm=TRUE), main="") {
      plot(zoo(x), type="n", xlab="Time", ylab="x", font.lab=2,
		main=main, bty="l", ylim=c(ylim[1], ylim[2]))
		lines(x, col=l.col, lwd=lwd, lty=lty)
		points(x, pch=pch, col=p.col, cex=p.cex)
		abline(h=mean(x, na.rm=T),col=1, lty=2)
}

# ts.detr.lo
# A function comparing the performance of loess and lowess.
# loess is a newer version of lowess, with different defaults.
# However, loess and lowess are not equivivalent.
# But almost identical results can be obtained by setting
# degree = 1 & family = "symmetric".

# x = a time series
# span = [from loess help] For a < 1, the neighbourhood includes proportion a of the points, and these have tricubic weighting
# (proportional to (1 - (dist/maxdist)^3)^3. For a > 1, all points are used, with the ‘maximum distance’
# assumed to be a^1/p times the actual maximum distance for p explanatory variables
# degree = the degree of the polynomials to be used, normally 1 or 2
# family = if "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M estimator is
# used with Tukey's biweight function

# See ?loess & & lowess for further information.

# This function is mainly for learning purposes, the main function for detrending
# is ts.detr.

# Dependencies: loess, lowess
#' @export
ts.detr.lo <- function(x, span=0.75, degree=2, family="symmetric") {

	if (class(x) != "ts") stop("Input object x must be a time series")

	fm1 <- loess(x ~ time(x), span=span, degree=degree, family=family)
	fm2 <- lowess(x ~ time(x), f=span)

	op <- par(list(mfrow=c(2,1),mar=c(4.5,4.5,2.25,1)))
	main.str <- paste("span = ", span, ", points used = ", length(x)*span, ", degree = ", degree, ", family = ", family, sep="")

	plot(x, xlab="time", ylab="x", font.lab=2, type="b", pch=16, lty=2, cex.lab=0.9)
	title(main = main.str, cex.main=0.9)
	points(fm1$x, predict(fm1), type="l", col=2, lwd=1)
	points(fm2$x, fm2$y, type="l", col=4, lwd=1)
	legend("topleft", c("loess","lowess"), col=c(2,4), lwd=c(2,2), bty="n", bg="white", cex=0.75)

	plot(residuals(fm1), type="l", col=2, ylim=range(c(residuals(fm1), x - fm2$y)), ylab="x'", cex.lab=0.9)
	title(main="Detrended series", cex.main=0.9)
	points(fm2$x, x - fm2$y, type="l", col=4)
	par(op)

	out <- ts.intersect(x, x.loess=residuals(fm1), x.lowess=x - fm2$y)
	invisible(out)
}

# ts.detr.step
# Step detrend a time series
# x = time series to be detrended
# g = a vector with group identities
# k = a vector with (known) levels of g. Length must match length(unique(g))
# base = if a vector k, has been supplied, base can be used to set the base level to be added.
# fun = name of function to use for each group identity, default ="mean"
# base refers to the position in vector k. So if k=c(100,200) and base=2, k[2] (= 200) is added to each count
# (after subtracting the group means).
# plot = TRUE returns a plot with the original and detrended series
#' @export
ts.detr.step <- function(x, g, k=NULL, base=NULL, fun="mean", plot=TRUE) {

  if (any(is.na(x))) cat("Time series contains missing values", "\n")
	if (class(x) != "ts") x <- ts(x)

  FUN <- match.fun(fun)

	if (!is.null(k)) {
		if (length(unique(g)) != length(k)) stop("Levels of group vector 'g' must match length of 'k'")
		x.g <- k[g]
		if (is.null(base)) x.tr <- (x - x.g)
		if (!is.null(base)) {
			if (base %in% 1:length(k) == FALSE) stop("The position of 'base' can not be found in vector 'k'")
			x.tr <- (x - x.g) + k[base]
		}
	}

	if (is.null(k)) {
		x.g <- tapply(x, g, FUN, na.rm=TRUE)[g]
		x.tr <- (x - x.g) + FUN(x,na.rm=TRUE)
	}

	# Calculate break-points for graph
	inds <- which(diff(g)==1)
	breaks <- time(x)[inds+1]
	xbr <- c(start(x)[1], rep(breaks, each=2), end(x)[1])

	if (is.null(k)) {means <- tapply(x,g,mean,na.rm=T); ybr <- rep(means, each=2)}
	if (!is.null(k)) ybr <- rep(k, each=2)

	# Plot
	if (plot == TRUE) {

		if (!is.null(base)) {
			ylims1 <- c(min(c(x,x.tr),na.rm=TRUE), max(c(x,x.tr),na.rm=TRUE)); ylims2 <- ylims1
		}
		if (is.null(base)) {
			if (!is.null(k)) {ylims1 <- range(x,na.rm=TRUE); ylims2 <- range(x.tr,na.rm=TRUE)}
			if (is.null(k)) {ylims1 <- c(min(c(x,x.tr),na.rm=TRUE), max(c(x,x.tr),na.rm=TRUE)); ylims2 <- ylims1}
		}

		par(mfrow=c(1,2))

		# Plot original series
		plot(x, type="n", xlab="Time", ylab="Population size", font.lab=2, bty="l",
			ylim=ylims1, main="Original time series")
		points(x, pch=16, cex=1.3)
		lines(x)
		lines(xbr,ybr,col=2, lty=2, lwd=2)

		# Plot detrended series
		plot(x.tr, type="n", pch=16, cex=1.3, xlab="Time", ylab="Population size", font.lab=2, bty="l",
			ylim=ylims2, main="Step-detrended time series")
		points(x.tr, pch=16, cex=1.3)
		lines(x.tr)
		abline(h=mean(x.tr,na.rm=T), lty=2, col=2, lwd=2)
		par(mfrow=c(1,1))
	}

	# Create output
	out <- ts.union(x, x.tr, x.g, g)
	invisible(out)
}
