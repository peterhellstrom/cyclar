# ar.acf
# Calculate and display the theoretical autocorrelation function for an ARMA-process
# ar & ma are vectors with ar(ma)-parameters
# lag.max = maximum number of lags
# pacf = TRUE/FALSE. If TRUE, the pacf is calculated, if FALSE the acf is calculated.
# n = ?

# Dependencies: plotrix

# Note that the ar component must be given as c((1+a1),a2) for an AR(2)-process.
#' @export
ar.acf <- function(ar, ma=0, lag.max=16, pacf=FALSE, n=FALSE, plot=TRUE,
                   title=NULL, output=list) {

	theor.acf <- ARMAacf(ar = ar, ma = ma, lag.max = lag.max, pacf = pacf)

	if (n!=FALSE & is.numeric(n)==TRUE) {
		acf.ci <- qnorm((1 + 0.95)/2) / sqrt(n)
	}

	x.acf.tp <- turnpoints(theor.acf)
	pk <- which(x.acf.tp$peaks==TRUE)
	z.pk <- cbind(Lag=as.numeric(names(theor.acf[pk])), ACF=theor.acf[pk])
	max.index <- which(z.pk[,2]==max(z.pk[,2]))
	max.acf <- cbind(Lag=z.pk[,1][max.index],ACF=z.pk[,2][max.index])
	dom.per <- max.acf[1]

	if (plot == TRUE) {
		ncoefs <- length(ar)
		labs <- paste("a", 1:ncoefs, sep="")
		coef.labs <- paste(labs,"=", round(ar,3))

		if (is.null(title) == FALSE) main.lab <- paste("Theoretical ACF, ", title, sep="")
		if (is.null(title) == TRUE) main.lab <- c("Theoretical ACF")

		plot(names(theor.acf),theor.acf, type="h", xlab="Lag", ylab="Autocorrelation",
			main=main.lab, font.lab=2, las=1)
		points(max.acf[1],max.acf[2],col=2,pch=16)
		legend("topright",coef.labs,bty="n",cex=1.3)
		abline(h=0)

		if (n!=FALSE & is.numeric(n)==TRUE) {
			abline(h=acf.ci, lty=2, col=4)
			abline(h=-acf.ci, lty=2, col=4)
		}
	}

	out <- output(
	acf = theor.acf,
	'Dominant period' = dom.per)
	out
}

# acf.calc
# This function shows how the acf is calculated
# Included only for learning purposes, use {stats} function acf.

# See e.g. Royama (1992) p. 17, equations 1.15 a-d
#' @export
acf2 <- function(x, lag.max=16) {

	mu <- mean(x)
	x <- x - mu # Center the series x
	var <- sum(x^2)
	x.lag <- embed(x=c(rep(NA,lag.max),x), dimension=lag.max+1)
	#autocov <- colSums(sweep(x.lag, 1, x.lag[,1], "*"), na.rm=T) # alternative, using sweep
	autocov <- colSums(x.lag * rep(x.lag[,1],ncol(x.lag)),na.rm=T)
	autocor <- autocov / var # var == autocov[1]
	names(autocor) <- names(autocov) <- 0:lag.max
	list(mean=mu, var=var, autocov=autocov, autocor=autocor)
}

#z <- ar2.sim2(mu=2, a1=-1, a2=-0.75)
#acf2(z,lag.max=20)
#acf(z,plot=FALSE)

# ar2.acf
# Autocorrelation function for an AR(2)-process
# A more explicit alternative to ARMAacf(ar=c(a1,a2)) (or ar.acf in this package)
# Cryer & Chan (2008), p. 73, equation 4.3.17 omit the tan(theta) term in psi
# Fuller (1996), p. 56, equation 2.5.11 has the correct formula (see also Priestley 1981)
# The formula given by Fuller matches the results from R's ARMAacf function.

# a1 & a2 are AR(2) parameters
# lags is an optional vector with lags for which the acf should be calculated.
# lag.max is the maximum number of lags
# subplot = TRUE/FALSE, adds an subplot of the a1/a2-space if plot = TRUE
# vadj & hadj = vertical and horizontal adjustment of subplot
#' @export
ar2.acf <- function(lags=NULL, a1, a2, lag.max=40, plot=FALSE,
                    subplot=TRUE, vadj=1.25, hadj=-1.5) {

	if ((1+a1)^2 + 4*a2 > 0) stop("Roots are not complex")

	if (is.null(lags)) k <- 0:lag.max
	if (!is.null(lags)) k <- lags

	R <- sqrt(-a2) # damping factor
	theta <- acos((1+a1)/(2*sqrt(-a2))) # frequency
	# 2*pi/theta # quasiperiod
	psi <- atan(((1-a2)/(1+a2)) * tan(theta)) # phase
	out <- (R^k)*(sin(theta*k + psi) / sin(psi))

	if (plot==TRUE) {
		curve((R^x)*(sin(theta*x + psi) / sin(psi)), from=0, to=lag.max, n=1001, las=1,
			xlab="", ylab=expression(rho(italic(k))), main="ACF for a stationary AR(2) with complex roots")
		mtext(text="k", side=1, line=2.5, font=3)
		abline(h=0, lty=2)
		str <- bquote(paste((1+beta[1]) == .(1+a1), ", ", beta[2] == .(a2)))
		mtext(text=str, side=1, cex=0.8, line=4)
		title(main=paste("Quasi-period =", round(2*pi/theta,3)), line=0.5)
		if (subplot == TRUE) {
		subplot(ar2.plot.simple(a1=a1,a2=a2,cex.axis=0.75,cex.lab=0.75,las=1),
        x="top", vadj=vadj, hadj=hadj)
	}
	}

	names(out) <- k
	out
}

#' @export
ar2.acf2 <- function(lags=NULL, phi1, phi2, lag.max=40, plot=FALSE,
                     subplot=TRUE, vadj=1.25, hadj=-1.5) {

	if (phi1^2 + 4*phi2 > 0) stop("Roots are not complex")

	if (is.null(lags)) k <- 0:lag.max
	if (!is.null(lags)) k <- lags

	R <- sqrt(-phi2) # damping factor
	theta <- acos(phi1/(2*sqrt(-phi2))) # frequency
	# 2*pi/theta # quasiperiod
	psi <- atan(((1-phi2)/(1+phi2)) * tan(theta)) # phase
	out <- (R^k)*(sin(theta*k + psi) / sin(psi))

	if (plot==TRUE) {
		curve((R^x)*(sin(theta*x + psi) / sin(psi)), from=0, to=lag.max, n=1001, las=1,
			xlab="", ylab=expression(rho(italic(k))), main="ACF for a stationary AR(2) with complex roots")
		mtext(text="k", side=1, line=2.5, font=3)
		abline(h=0, lty=2)
		str <- bquote(paste(phi[1] == .(phi1), ", ", phi[2] == .(phi2)))
		mtext(text=str, side=1, cex=0.8, line=4)
		title(main=paste("Quasi-period =", round(2*pi/theta,3)), line=0.5)
		if (subplot == TRUE) {
		subplot(ar2.plot.simple(a1=phi1-1,a2=phi2,cex.axis=0.75,cex.lab=0.75,las=1,
      xlab=expression(phi[1]), ylab=expression(phi[2]), par.name="phi"),
        x="top", vadj=vadj, hadj=hadj)
	}
	}

	names(out) <- k
	out
}

# ar2.acf.test
# This function repeatedly simulates an AR(2)-process with parameters selected
# along a k-contour (period length), and calculates the dominant peak of the ACF-function for
# each simulation.

# k = periode length
# npoints = number of points to generate along k-contour k
# n = length of each arima-simulation
# nrepl = number of replicates per level
# k.per = the period lengths to extract and plot

# Dependencies: ar2.k.contour, ar2.period, ar2.sim, ts.diag.acf
#' @export
ar2.acf.test <- function(k=3, npoints=20, n=100, nrepl=100, k.per=c("3","6","7","9"), from=-0.95, to=-0.05, ...) {

	sim.parms <- ar2.k.gen(k=k, length.out=npoints, from=from,to=to, ...)
	xv <- ar2.period(a1=sim.parms$a1,a2=sim.parms$a2)

	ar2.plot(k=k)
	points(1+sim.parms$a1,sim.parms$a2,pch=16,col=2)

	out <- sapply(1:npoints, function(i) {
		sim.dat <- t(replicate({
			sim.ts <- ar2.sim(n=n, mu=0, a1=sim.parms$a1[i], a2=sim.parms$a2[i], sd=sqrt(0.2), plot=FALSE)
			period <- ts.diag.acf(fun="acf", x=sim.ts, plot=FALSE, max.only=TRUE)
			as.numeric(period[1])}, n=nrepl))
		sim.dat
	})

	out.tab <- apply(out,2,table)
	m <- seq(min(out),max(out),1)

	out.tab.prop <- t(sapply(1:npoints, function(i) {
		temp <- as.data.frame(out.tab[[i]])
		colnames(temp) <- c("m","Freq")
		temp <- merge(m, temp, by=1, all=TRUE)
		temp[,2]/nrepl
	}))
	colnames(out.tab.prop) <- m
	out.tab.prop

	# Plot estimated dominant acf-period and distribution
	# Create a new matrix, with only the lag specified by k.per
	find.mat <- out.tab.prop[,k.per]

	# Plot output
	main.lab <- k
	plot(x=range(xv[,2]),y=c(0,1),type="n",bty="l",xlab="Delayed density-dependence",ylab="Prob (k = dominant ACF)",
		main=paste("Expected cycle period =", main.lab), font.lab=2,las=1)
		for (i in 1:length(k.per)) points(xv[,2], find.mat[,i], pch=i, col=i)
	legend("topright", legend=k.per, pch=1:length(k.per), col=1:length(k.per), bty="n", cex=1.3, title="Cycle period (k)")

	out.list <- list(out, out.tab.prop, find.mat)
	invisible(out.list)
}
