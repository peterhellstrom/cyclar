# ts.diag.acf
# Extract the dominant period of an acf
# x is a time-series
# fun is the name of a function (within quotes), either "acf" or "pacf"
# ci = desired siginficance level, default is 0.95
# plot = TRUE returns a plot
# if max.only=TRUE, only maximum of acf or pacf is printed
# lag.max is the maximum number of lags
# out = desired construction of output, use list, rbind or cbind

# Method 2 is included ONLY for development purposes at the moment.
# Use method == 1

# Dependencies: pastecs
#' @export
ts.diag.acf <- function (x, fun="acf", ci=0.95, method=1, max.only=FALSE, lag.max=NULL, output=rbind,
    plot=TRUE, main=paste("Series: ", deparse(substitute(x)), sep="")) {

	if (class(x)[1] == "mts") stop("This function can not handle multiple time series")
	if (method == 2) cat("Warning: method == 2 does not return correct estimates")

	FUN <- match.fun(fun)
	x.acf <- FUN(x, lag.max=lag.max, type="correlation", plot=FALSE, na.action=na.pass)

	z <- cbind(Lag=x.acf$lag, ACF=x.acf$acf)
	# Turnpoins of ACF
	x.acf.tp <- turnpoints(z[,"ACF"])

	if (method == 1) {
		# Peaks
		pk <- which(x.acf.tp$peaks==TRUE)
		z.pk <- cbind(Lag=z[,1][pk], ACF=z[,2][pk])
		max.index <- which(z.pk[,2]==max(z.pk[,2]))
		max.acf <- cbind(Lag=z.pk[,1][max.index],ACF=z.pk[,2][max.index])
		names(max.acf) <- c("Lag","ACF")
		# Pits
		pt <- which(x.acf.tp$pits==TRUE)
		z.pt <- cbind(Lag=z[,1][pt], ACF=z[,2][pt])
		# Extract max & min
		min.index <- which(z.pt[,2]==min(z.pt[,2]))
		min.acf <- cbind(Lag=z.pt[,1][min.index],ACF=z.pt[,2][min.index])
		names(min.acf) <- c("Lag","ACF")
	}

	if (method == 2) {
		z.pk <- z[x.acf.tp$peaks,]
		z.pk <- z.pk[order(z.pk[,2],decreasing=TRUE), ]
		z.pt <- z[x.acf.tp$pits,]
		z.pt <- z.pt[order(z.pt[,2],decreasing=FALSE), ]
		# Extract max & min
		max.acf <- z.pk[1, ]
		min.acf <- z.pt[1, ]
	}

	# Significance level
	p.crit <- qnorm((1 + ci)/2) / sqrt(x.acf$n.used)

	# Plot
	if (plot == TRUE) {
		leg.str <- list(
			paste("max","[",round(max.acf[1],2),"] = ", round(max.acf[2],4),sep=""),
			paste("min","[",round(min.acf[1],2),"] = ", round(min.acf[2],4),sep="")
		)

		plot(x.acf, main=main)
		legend("topright", do.call("expression", leg.str), pch=c(16,16), col=c(2,4), bty="n",
			title=paste("Significance level: +/-", round(p.crit,4)))
		points(max.acf[1],max.acf[2],col=2,pch=16)
		points(min.acf[1],min.acf[2],col=4,pch=16)
	}

	if (max.only == TRUE) {
		if (fun == "acf") out <- c(max.acf, p.crit = p.crit)
		if (fun == "pacf") out <- c(min.acf, p.crit = -p.crit)
	}

	else # Output
	out <- output(
		max = c(max.acf, p.crit = p.crit),
		min = c(min.acf, p.crit = -p.crit)
	)

	out
}
