# Estimate subsets of the maximal ar2-model
# x is an object created by the function growth.rate().
#' @export
ar2.lm <- function(x, method="raw", model="lm", p=2, q=0) {

	if (method=="growth rate") {
		y <- x[,"R"]
		x1 <- x[,"Lag0"]
		x2 <- x[,"Lag1"]
		}

	if (method=="raw") {
		y <- x[,"Lag0"]
		x1 <- x[,"Lag1"]
		x2 <- x[,"Lag2"]
		}

	if (model=="lm") {
		fm0 <- lm(y ~ 1)
		fm1 <- lm(y ~ x1)
		fm2 <- lm(y ~ x2)
		fm12 <- lm(y ~ x1 + x2)
		}
	if (model=="glm") {
		fm0 <- glm(y ~ 1)
		fm1 <- glm(y ~ x1)
		fm2 <- glm(y ~ x2)
		if (method=="growth rate") fm12 <- glm(y ~ x1 + x2)
		if (method=="raw") fm12 <- glm(y ~ offset(x1) + x1 + x2)
		}
	if (model=="gls") {
		fm0 <- gls(y ~ 1, corr=corARMA(p=p, q=q))
		fm1 <- gls(y ~ x1, corr=corARMA(p=p, q=q))
		fm2 <- gls(y ~ x2, corr=corARMA(p=p, q=q))
		fm12 <- gls(y ~ x1 + x2, corr=corARMA(p=p, q=q))
		}

	mod <- list(fm0,fm1,fm2,fm12)
	names(mod) <- paste("fm",c(0,1,2,12),sep="")

	coefs <- rbind(fm0=c(coef(fm0),NA,NA), fm1=c(coef(fm1),NA), fm2=c(coef(fm2),NA), fm12=coef(fm12))
	colnames(coefs) <- c("a0","a1","a2")

	nobs <- length(y)

	out <- list(
	data = x,
	method = method,
	model = model,
	ictab = ICtab(fm0,fm1,fm2,fm12, weights=TRUE, delta=TRUE, base=TRUE, sort=TRUE, type="AICc", nobs=nobs),
	models = mod,
	coefs = coefs)
out
}

# arima() reports the mean as the intercept. But the correct intercept
# can be obtained with this function.
# Note that ar.ols() reports the "correct" intercept.
#' @export
arima.coef <- function(x) {
	if (class(x) != "Arima") stop("Input data must be an Arima-class object")
	coefs <- coef(x)
	n <- length(coefs)
	z <- sum(coefs[1:(n-1)])
	b0 <- coefs["intercept"]*(1-z)
	out <- c(coefs[1:(n-1)], b0)
	out
}

# ar2.intercept
# arima() confusingly reports the mean as the intercept.
# But the correct intercept can be obtained with this function.
# Note that ar.ols() reports the "correct" intercept.
# Enter either mu OR a0, and a1 & a2.

# mu = mean of time series on log-scale
# a0 = intercept of AR(2)-model

# Source: Shumway & Stoffer's website has an example for an AR(1)-process.
# Add better description...

# Note: the intercept is used also in ar2.sim2.
# The calculations in that function should be based on this function.
#' @export
ar2.intercept <- function(mu=NULL,a0=NULL,a1,a2) {

  if (is.null(mu) & is.null(a0)) stop("Supply mu or a0")

	if (!is.null(mu) & is.null(a0)) {
		a0 <- mu * (-(a1+a2))
		out <- c(a0=a0)
    }
    else if (!is.null(a0) & is.null(mu)) {
		mu <- a0 / (-(a1+a2))
		out <- c(mu=mu)
    }

  out
}

#' @export
ar2.intercept2 <- function(mu=NULL,phi0=NULL,phi1,phi2) {

  if (is.null(mu) & is.null(phi0)) stop("Supply mu or phi0")

	if (!is.null(mu) & is.null(phi0)) {
		phi0 <- mu * (1 - (phi1 + phi2))
		out <- c(phi0=phi0)
    }
    else if (!is.null(phi0) & is.null(mu)) {
		mu <- phi0 / (1 - (phi1 + phi2))
		out <- c(mu=mu)
    }

  out
}

#Example: (NOT RUN)
#ar.coef <- c(-1.7679548, -0.5897545, 0.2)
#ar2.period(ar.coef[1], ar.coef[2])
#mu <- 2.1
#x <- arima.sim(list(order=c(2,0,0), ar=c(1+ar.coef[1], ar.coef[2])), sd=sqrt(ar.coef[3]), n=10000) + mu
#plot(x)
#mean(x)
#arima(x, order = c(2, 0, 0))  #Intercept is really the mean!
#ar.ols(x, order=2, demean=F, intercept=T)
#Intercept: mu * (-(a1+a2))
#mu * -(sum(ar.coef))


# Fit ar2-model to univariate time series
# UPDATE THIS FUNCTION TO INCLUDE A COVARIATE, xreg, in arima part
# ALSO check "include.mean" argument and if/how it works, and connection with demean.
#' @export
ar.univar <- function(x, log=TRUE, constant=0, aic=FALSE, order.max=2,
				method="arima", d=0, q=0, include.mean=TRUE, plot=TRUE) {

	if (log == TRUE) X <- log(x + constant)
	if (log == FALSE) X <- (x + constant)

	ar <- switch(method,
		"arima" = arima(X, order=c(order.max,d,q), include.mean=include.mean),
		"ar.yw" = ar.yw(X, aic=aic, order.max=order.max),
		"ar" = ar(X, aic=aic, order.max=order.max, method="yule-walker"),
		"ar.ols" = ar.ols(X, aic=aic, order.max=order.max, demean=FALSE, intercept=TRUE))

	# Note that Beta1 in output for coefficients is really (1 + Beta1)! Perhaps change this in output?
	if (method == "arima") {
		coefs <- matrix(ncol=length(ar$coef), nrow=2)
		rownames(coefs) <- c("Estimate","SE")
		colnames(coefs) <- c(paste("Beta",1:(ncol(coefs)-1),sep=""), "Intercept")
		coefs[1,] <- ar$coef
		coefs[2,] <- sqrt(diag(ar$var.coef)) # getAnywhere(print.Arima)
		sigma2 <- ar$sigma2
	}

	if (method == "ar.yw" | method == "ar" | method == "ar.ols") {
		coefs <- matrix(ncol=length(ar$ar)+1, nrow=2)
		rownames(coefs) <- c("Estimate","SE")
		colnames(coefs) <- c(paste("Beta",1:(ncol(coefs)-1),sep=""), "Intercept")

		if (method == "ar.ols") {
			coefs[1,] <- c(ar$ar, ar$x.intercept)
			coefs[2,] <- c(ar$asy.se$ar, ar$asy.se$x.mean)
		}
		if (method == "ar.yw" | method == "ar") {
			coefs[1,] <- c(ar$ar, ar$x.mean)
			coefs[2,] <- rep(NA, ncol(coefs)) # ar$asy.var.coef
		}

		sigma2 <- ar$var.pred
	}

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
		points(coefs[1,1], coefs[1,2], col="red", cex=1.4, pch=16)

		# 2) Time series
		plot.ts(X, main=paste("Time series"), xlab="Time", ylab="Population", font.lab=2, las=1)

		# 3) Auto-correlation
		acf.c <- acf(X, type="correlation", main="Autocorrelation", font.lab=2, las=1)

		# 4) Partial rate autocorrelation
		acf.p <- acf(X, type="partial", main="Partial autocorrelation", font.lab=2, las=1)

		par(mfrow=c(1,1))
	}

	out <- list(
		ar = ar,
		coef = coefs,
		sigma2 = as.numeric(sigma2),
		log = log,
		constant = constant,
		aic = aic,
		order.max = order.max,
		method = method,
		x=x,
		X=X,
		n=length(x))

out
}
