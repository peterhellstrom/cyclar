# ts.moving
# Custom function for calculating window-based (moving) statistics.
# Moving averages can be calculated with several packages.
# For instance, the zoo-package has functions such as rollapply, rollmean, rollmedian, etc.
# ts.moving is yet another quick function that can be used to calculate moving/sliding window statistics.
# Input arguments:
# x = a (time) series
# width = number of data points in window (assumed to be centered along the midpoint of the window)
# fun = name of the function used to calculate the corresponding statistic in each window (e.g. mean, sd, var or also custom functions).
# output = c("vector","ts"). Output can either be as a vector or as a time series object.
# names = c("mid";"span"). If output = "vector", names can either be printed as the midpoint of the window, or as the span of the window.
#' @export
ts.moving <- function(x, width=3, fun="mean", output="vector", names="mid", ...) {

	FUN = match.fun(fun)

	N <- length(x)
	n <- width
	if (class(x) != "ts") {
		x <- ts(x)
		tsp(x) <- c(1, N, 1)
    }
	if (class(x) == "ts") i <- time(x)
	j <- N-n+1 # number of windows
	w <- cbind(i[1:j], i[1:j]+n-1) # bind start- & endpoints together

	out <- sapply(1:j, function(k) FUN(window(x, start=w[k,1], end=w[k,2]), ...))
	mds <- floor(sapply(1:j, function(k) median(w[k,1]:w[k,2]))) # midpoints of window
	sps <- sapply(1:j, function(k) paste(w[k,1],"-",w[k,2], sep="")) # range of each window

	if (output == "ts") {
		out <- ts(out, start=mds[1], end=mds[j])
	}

	if (output == "vector") {
		if (names=="mid") names(out) <- mds
		if (names=="span") names(out) <- sps
	}

	out
}


# ts.window
# Calculate a given statistic (fun) over a specified sliding window (specified by width).
# x is the time series
# This functions is basically a wrapper function built on rollapply in package zoo
# rollapply(x, width, FUN)
# This function is only intended for single time series.
# Use the companion function ts.window.n for multiple time series.
# Two plot options: plot.type="single" plots the original series in the same graph: x on the left y-axis and the
# sliding window statistic on the right y-axis. plot.type="multiple" plots the two series in separate panels.

# NOTE: ts.window & ts.window.n can also use custom-made functions as input,
# see example file for examples!

# Dependencies: {zoo}
#' @export
ts.window <- function(x, width, fun, plot=TRUE, plot.type="single",
              type=c("l","l"), col=c(1,2), pch=c(15,16), lty=c(1,1), cex=c(1,1), ...) {

	if (class(x)[1] == "mts") stop("Use the the ts.window.n - function for multiple series")
	if (!is.null(ncol(x))) stop("Use the the ts.window.n - function for multiple series")

	FUN <- match.fun(fun)
	x <- zoo(x)
	y <- rollapply(x, width=width, FUN=FUN, ...)

	if (plot == TRUE) {

		z <- merge(x,y)

		if (plot.type == "single") {
			op <- par(mar=c(5,5,4,5))
			plot(z$x, xlab="Time", ylab="x", main=paste("Window =",width),
				font.lab=2, type=type[1], lty=lty[1], col=col[1], pch=pch[1], cex=cex[1])
			par(new=TRUE)
			plot(z$y, ann=FALSE, xaxt="n", yaxt="n", type=type[2], lty=lty[2], col=col[2], pch=pch[2], cex=cex[2])
			axis(side=4)
			mtext(side=4, line=2.5, text=deparse(substitute(fun)), font=2)
			par(op)
		}
		if (plot.type == "multiple") {
			plot(z, xlab="Time", ylab=c("x",deparse(substitute(fun))), main=paste("Window =",width),
				plot.type="multiple", col=col, type=type, pch=pch, lty=lty, font.lab=2)
		}
	}

	out <- list(statistic=fun, width=width, x=x, window=y)
	invisible(out)
}

# ts.window.n
# Calculates a given statistic (entered as a function name, argument 'fun') over a sliding window (argument 'width').
# This function is only used for multiple time series.
# See ts.window for explanation of input arguments.
# plot.type="single" plots a 2-panel plot (original series x and sliding window series).
# plot.type="multiple" plots each series in a separate panel.

# Dependencies: {zoo}
#' @export
ts.window.n <- function(x, width, fun, plot=TRUE, plot.type="single",
              type=c("l","l"), col=c(1:dim(x)[2]), pch=c(1:dim(x)[2]), lty=c(1,1), cex=c(1,1), ...) {

	n <- dim(x)[2]
	FUN <- match.fun(fun)
	x <- zoo(x)
	y <- rollapply(x, width=width, FUN=FUN)

	if (plot == TRUE) {
		if (plot.type == "single") {
			par(mfrow=c(1,2))
			plot(x, plot.type=plot.type, type=type[1], col=col, pch=pch, lty=rep(lty[1], n), cex=cex[1],
				xlab="Time", ylab="x", main="Time series",font.lab=2)
			plot(y, xlim=range(time(x)), plot.type=plot.type, type=type[2], col=col, pch=pch, lty=rep(lty[2], n), cex=cex[2],
				xlab="Time", ylab=deparse(substitute(fun)), main=paste("Window =",width), font.lab=2)
		par(mfrow=c(1,1))
		}

		if (plot.type == "multiple") {
			z <- merge(x,y)
			plot(z, col=col, xlab="Time", ylab=c(rep("x",n),rep(deparse(substitute(fun)),n)), font.lab=2, main="Sliding window statistic")
		}
	}

	out <- list(statistic=fun, width=width, x=x, window=y)
	invisible(out)
}


# Calculate changes over time in crude variability
# CHECK WINDOW LENGTH AND CALCULATIONS!
# Notice that number of used points here is swin + 1
# which is not the case in the funciton swin.var (where number of points is swin).
#' @export
swin.acf <- function(x, swin, lag.max=16) {

	# x is a time series, swin is size of sliding window

	if (class(x) != "ts") stop("Input variable x must be a time series object")
	x <- na.omit(x)

	start <- start(x)[1]
	end <- end(x)[1]

	n <- length(x)
	n.win <- length(start:(end - swin))

	tinds <- start:end
	n.used <- swin + 1

	labs <- sapply(1:n.win, function(i) paste(tinds[i], "-", tinds[i + swin],sep=""))

	acf.win <- lapply(1:n.win, function(i) acf(window(x, start=tinds[i], end=tinds[i + swin]), lag.max=lag.max, plot=FALSE, na.action=na.omit))

	check <- sapply(1:n.win, function(i) length(acf.win[[i]]$acf))
	if (length(which(check < lag.max)) > 0) stop("Check start of calculations, number of data points < lag.max")

	acf.tab <- sapply(1:n.win, function(i) acf.win[[i]]$acf)

	rownames(acf.tab) <- 0:lag.max # Code breaks here (if acf.tab is a list and not a matrix or data.frame)
	colnames(acf.tab) <- labs

	npoints <- sapply(1:n.win, function(i) acf.win[[i]]$n.used)

	list(acf = acf.tab, npoints = npoints)

}

# Calculate changes over a specified sliding window:
# There's a slow for loop in here that should be re-coded!!!
#' @export
swin.var <- function(x, swin, stat=sd, plot=TRUE, col=2, lty=1,
				mtext=mtext(side=3, line=0, deparse(substitute(x)))) {

# x is a time series, swin is size of sliding window

	start <- start(x)[1]
	end <- end(x)[1]

	if ((class(x)[1]=="ts")==TRUE) n.series <-  1
	if ((class(x)[1]=="mts")==TRUE) n.series <- ncol(x)

	n.win <- length(start:(end - swin)) + 1
	tinds <- start:end
	labs <- sapply(1:n.win, function(i) paste(tinds[i], "-", tinds[i + (swin-1)],sep=""))

	if (n.series == 1) {
		res <- sapply(1:n.win, function(i) stat(window(x, start=tinds[i], end=tinds[i + (swin-1)])))
		names(res) <- labs
		}
	if (n.series > 1) {
		res <- matrix(nrow=n.win, ncol=ncol(x), dimnames=list(labs, colnames(x)))
		#for (j in 1:n.series) {
		#	res[,j] <- t(sapply(1:n.win, function(i) stat(window(x[,j], start=tinds[i], end=tinds[i + (swin-1)]))))
		#	}
		res <- sapply(1:n.series, function(j) {
					t(sapply(1:n.win, function(i) stat(window(x[,j], start=tinds[i], end=tinds[i + (swin-1)]))))
					})
		rownames(res) <- labs
		colnames(res) <- colnames(x)
		}

	if (plot==TRUE & n.series == 1) {
		plot(res,type="l",xaxt="n",col=col,lty=lty, bty="l",
		xlab="Center of time window", ylab=deparse(substitute(stat)), font.lab=2,
		main=paste("Sliding window (", start,"-",end,"): ",swin," steps, statistic: ", deparse(substitute(stat)), sep=""))
		axis(side=1,at=1:n.win,labels=names(res))
		}

	if (plot==TRUE & n.series > 1) {
		plot(1:n.win,type="n",xaxt="n",col=col,lty=lty, bty="l", ylim=range(res,na.rm=T),
		xlab="Center of time window", ylab=deparse(substitute(stat)), font.lab=2,
		main=paste("Sliding window (", start,"-",end,"): ",swin," steps, statistic: ", deparse(substitute(stat)), sep=""))
		axis(side=1,at=1:n.win,labels=labs)
			for (j in 1:ncol(x)) points(1:n.win, res[,j], type="l", col=j, lty=lty)
		if(n.series < 10) legend("topright", legend=colnames(x), col=1:ncol(x), lty=lty, bty="n", cex=0.75)
		}

	res
}
