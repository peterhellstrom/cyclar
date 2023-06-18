# ar2.lag
# calculates growth rate, r, and creates a data.frame with lagged variables, x[t], x[t-1], ..., x[t-n], where n equals the argument lag.max.
# x = a time series object.
# x must be a time series object, and should preferably bepopulation densities on log-scale.
# The growth rate calculated here is r =(log(lambda)), and is calculated with the diff function.
# Using diff is only valid of the population counts/densities are entered on log-scale, as
# r[t] = log(N[t+1]/N[t]) = log(N[t+1]) - log(N[t]) = x[t+1] - x[t]

# If the input variable is on natural scale, i.e. Nt = exp(x[t]), use the arguments log, constant & base to transform the time series.
# log = TRUE/FALSE, defaults to FALSE. if log = TRUE, the input variable x is log-transformed.
# The logarithmic base can be changed with the base argument. Default is base=exp(1), i.e. the natural logarithm.
# If the time series N[t] contains zeros, a constant can be added prior to log-transformation with the constant argument.
# constant = numeric value added to x if log = TRUE, so that the analyzed variable is x = log(x + constant, base)
# Default is constant=0. Note that arguments constant & base are only relevant if log=TRUE.
# lag.max = number of lags [note that the lags created by ar2.lag is equal to negative k in the commonly used lag-function, lag(x,k)].
# output: a data.frame with lag.max + 1 columns. The left lag.max-number of columns contain the lagged variables, named Lag0 to Lag[lag.max].
# The rightmost column is called r, and contains growth rates on log scale.

# To get the lags correct, the lagged series must be tied together.
# R has several functions for such operations: lag(), embed(), ts.intersect(), ts.union()
# In ar2.lag, a custom function (lagmatrix) is used to create the lagged variables.
# The lagmatrix-function calls embed() in order to create the lagged variables.
# embed does not output a time series (rather a matrix), so the lagged variable-matrix
# is converted to a time series object within ar2.lag.
# lagmatrix is also available for independent use.

# Dependencies: lagmatrix
#' @export
ar2.lag <- function(x, lag.max=4, log=FALSE, constant=0, base=exp(1)) {

	if (class(x) != "ts") stop("object 'x' must be a time series")
	if (log) x <- log(x + constant, base=base)

	# Create lagged variables
	x.lags <- lagmatrix(x=x, lag.max=lag.max)
	x.lags <- ts(x.lags,start=start(x)[1])

	# Calculate growth rate
	r <- ts(diff(x), start=start(x)[1])

	# Create output
	out <- cbind(x.lags, r)
	colnames(out) <- c(paste("Lag",0:lag.max,sep=""),"r")
	out
}

# ar2.lag.plot
# Plots x[t] against x[t-n] or r[t] against x[t-n]
# x is preferably densities/counts on log scale. Calculation of r (growth rate) assumes that x is on log-scale.
# Natural scale-data can also be supplied, but if so the data must be log-transformed with the arguments
# log, constant & base.

# x = a time series (must be a time series object)
# variable = c("x","r") x=population counts/densities, r=population growth rate
# If variable = "x", ar2.lag.plot produces of x[t] against x[t-1], x[t-2], ..., x[t-n]
# If variable = "r", ar2.lag.plot produces of r[t] against x[t-1], x[t-2], ..., x[t-n]
# lag.max = maximum number of lags.
# (ar2.lag.plot1 does not have a lag.max argument, but instead a lag argument with the same meaning).
# log = TRUE/FALSE. If TRUE, the series x is log transformed according to x <- log(x+constant,base)
# constant = constant added to x prior to log-transformation.
# base = logarithmic, default is exp(1), natural logarithm.
# do.lines = TRUE/FALSE. If do.lines=TRUE, a loess polynomial is fitted to each plot.
# The degree of smoothing is controlled with the arguments span & degree.
# span = proportion of points used in each neighbourhood
# degree = c(0,1,2). The degree of the polynomials to be used, normally 1 or 2.
# See ?loess for further info

# NOTE: Only tested (and so far only intended) for annual data!

# Two functions are used, the main function is ar2.lag.plot1, which is called repeatedly by ar2.lag.plot.
# The intended use is to call ar2.lag.plot, but ar2.lag.plot1 can be used independently.
#' @export
ar2.lag.plot1 <- function(x, variable="r", lag, log=FALSE, constant=0, base=exp(1),
                          do.lines=TRUE, span=0.75, degree=1, ...) {

	if (class(x) != "ts") stop("Object 'x' must be a time series")

	if (log == TRUE) x <- log(x + constant, base=base)

	if (variable =="r") {
		r.x <- diff(x)
		# Set time series attributes so that growth rate is correctly aligned, starting at start(x)
		tsp(r.x) <- c(start(x)[1], end(x)[1]-1, 1)
		# Use ts.intersect together with lag to construct variables for plotting. Note that k must be negative in call to lag!
		out <- ts.intersect(lag(x, k=-lag), r.x)
		# y-axis label
		y.str <- expression(r[t])
	}

	if (variable == "x") {
		# See comments under variable=="r"
		out <- ts.intersect(lag(x, k=-lag), x)
		y.str <- expression(x[t])
	}

	# Create plot/output data
	xx <- as.numeric(out[,1])
	yy <- as.numeric(out[,2])

	# x-axis labels
	if (lag == 0) x.str <- expression(x[t])
	if (lag != 0) x.str <- substitute(x[t-a], list(a=lag))

	# Create plot
	plot(xx,yy,
		xlab = x.str,
		ylab = y.str,
		...)
	if (variable == "r") abline(h=0, lty=2) # add line at h=0
	# Add loess fit
	if (do.lines) {
		fm1 <- loess(yy ~ xx, span=span, degree=degree)
		xv <- seq(from=min(xx,na.rm=T), to=max(xx,na.rm=T), length.out=2*length(x))
		yv <- predict(fm1, newdata=(x=xv))
		lines(xv,yv,col=2,lty=2,lwd=1)
	}

	invisible(list(x=xx,y=yy))
}

#' @export
ar2.lag.plot <- function(x, variable="r", lag.max, log=FALSE, constant=0, base=exp(1),
                         do.lines=TRUE, span=0.75, degree=1, ...) {

	if (variable == "r") {
		n <- lag.max + 1
		op <- par(mfrow = n2mfrow(n))
		sapply(1:n, function(i) ar2.lag.plot1(x=x, variable=variable, lag=i-1, log=log, constant=constant, base=base,
			do.lines=do.lines, span=span, degree=degree, ...))
		par(op)
	}

	if (variable == "x") {
		n <- lag.max
		op <- par(mfrow = n2mfrow(n))
		sapply(1:n, function(i) ar2.lag.plot1(x=x, variable=variable, lag=i, log=log, constant=constant, base=base,
			do.lines=do.lines, span=span, degree=degree, ...))
		par(op)
	}

}

# ar2.gr
# In contrast to ar2.lag & ar2.lag.plot, the input series x should be on natural scale
# (which motivates the inclusion of ar2.gr in cyclar).
# Calculate growth rate, input is a vector with population counts
# Choose between method="lambda" or method="r".

# This function is available in consEcol.demogR.r (but called gr in that "package")
#' @export
ar2.gr <- function(x, method="r", plot=TRUE, ...) {

	if (class(x) != "ts") stop("Object 'x' must be a time series object")

	na.inds <- which(is.na(x))
	if (length(na.inds) > 0) message("Object 'x' contains missing/NA values")

	g <- switch(method,
		lambda = x[-1] / x[-length(x)],
		r = diff(log(x))
	)

	# Time component of the variable g is wrong (it is associated with t instead of t-1)
	# and must therefore be changed
	g <- ts(g, start=start(x)[1], deltat=deltat(x))

	title.text <- switch(method,
		lambda = expression(bold(paste("Growth rate =", lambda))),
		r = expression(bold(paste("Growth rate = log(", lambda, ")")))
	)

	if (plot) {

		# LOESS for x
		loess.x <- loess(x ~ time(x))
		pred.loess.x <- predict(loess.x, newdata=time(x))
		# LOESS for g
		loess.g <- loess(g ~ time(g))
		pred.loess.g <- predict(loess.g, newdata=time(g))
		# LOESS for dd in g
		dd.x <- x[-length(x)]
		loess.dd <- loess(g ~ dd.x)
		pred.dd.x <- seq(min(x,na.rm=T),max(x,na.rm=T),length=101)
		pred.loess.dd <- predict(loess.dd, newdata=pred.dd.x)

		# Density dependence in growth rates
		lm.dd <- lm(g ~ x[-length(x)])
		coefs <- summary(lm.dd)$coefficients[,1]
		p <- summary(lm.dd)$coefficients[2,4]
		lm.text <- paste(round(coefs[1],6), " + ", round(coefs[2],6), "x, p =", round(p,6), sep="")

		par(mfrow=c(2,2))
		plot(x,las=1,type="n",main="Population counts",xlab="Time", ylab="Counts")
			lines(x,lty=2)
			points(x,pch=16,cex=1)
			lines(ts(pred.loess.x,start=start(x)[1]),col=2,...)
		plot(g,las=1,type="n",main=title.text,xlab="Time",ylab="Growth rate")
			lines(g,lty=2)
			points(g,pch=16,cex=1)
			lines(ts(pred.loess.g,start=start(x)[1]),col=2,...)
		plot(x[-length(x)],las=1,g,pch=16,cex=1,xlab="Population density",ylab="Growth rate",main="Density-dependence in growth rates")
			points(pred.dd.x,pred.loess.dd,col=2,type="l",...)
		acf(g, main="Auto-correlation in growth rates", na.action=na.pass)
		if (length(na.inds) > 0) title(sub="'x' contains NA values")
		par(mfrow=c(1,1))
	}
	# Print output
	invisible(list(TimeSeries=x, GrowthRate=g))
}
