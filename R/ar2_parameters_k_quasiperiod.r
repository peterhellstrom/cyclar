# ar2.period
# Find period length for given ar2-parameters (a1 & a2)
# It is possible to enter vectors of a1 & a2-values (vectors must be of equal length).
# See Hamilton 1994, p. 10-18
# See also Hamilton 1994 p. 704 and onwards for trigonometry
# Two equal methods. method = c("eigen", "polyroot")
# method = "polyroot" is substantially faster for large datasets.
# Plot is only available for method="eigen", which return a plot of the eigenvalues in complex plane.

# Dependencies: ar2.ipv
#' @export
ar2.period <- function(a1, a2, method=c("eigen","polyroot"),
                       output=cbind, print.all=FALSE, plot=FALSE) {

  if (identical(length(a1),length(a2)) == FALSE) stop("Input vectors a1 and a2 must be of equal length")
  n <- length(a1)
  per <- numeric(n)
  a1 <- as.numeric(a1)
  a2 <- as.numeric(a2)

  method <- match.arg(method)

  if (method == "eigen") {
		for (i in 1:n) {
			mat <- matrix(c((1+a1[i]),a2[i],1,0),ncol=2,nrow=2,byrow=T)
			A <- eigen(mat)
			a <- Re(A$values) # Real part
			b <- Im(A$values) # Imaginary part
			R <- Mod(A$values) # Modulus
			phi <- acos(a[1]/R[1])
			# Period, general solution:
			# (2*pi)/(acos((1+a1)/(2*sqrt(-a2))))
			per[i] <- (2*pi)/phi

      if (!print.all) out <- output(a1=a1, a2=a2, period=per, ipv=ar2.ipv(a1,a2))
      if (print.all) out <- output(a1=a1, a2=a2, period=per, ipv=ar2.ipv(a1,a2), Re=a, Im=b, Mod=R, phi=phi)

			if (plot == TRUE) {
				plot(x=c(-1.1,1.1), y=c(-1.1,1.1), type="n", xlab="Re", ylab="Im", main="Eigenvalues in complex plane",
					font.lab=2, las=1)
				abline(v=0, lty=2)
				abline(h=0, lty=2)
				curve(sqrt(1-x^2), from=-1, to=1, add=T)
				curve(-sqrt(1-x^2), from=-1, to=1, add=T)
				segments(x0=0,y0=0,x1=a[1],y1=b[1], col=2, lty=2)
				points(a,b,col=2,pch=16)
			}
		}
	}

	if (method == "polyroot") {
		for (i in 1:n) {
		u <- polyroot(c(1,-(1+a1[i]),-a2[i]))
		per[i] <- 2*pi/Arg(u[1])
		}
		out <- output(a1=a1, a2=a2, period=per, ipv=ar2.ipv(a1,a2))
	}

	out
}

#' @export
ar2.period2 <- function(phi1, phi2, method="eigen", output=cbind, plot=FALSE) {

  if (identical(length(phi1),length(phi2)) == FALSE) stop("Input vectors phi1 and phi2 must be of equal length")
  n <- length(phi1)
  per <- numeric(n)
  phi1 <- as.numeric(phi1)
  phi2 <- as.numeric(phi2)

  if (method == "eigen") {
		for (i in 1:n) {
			mat <- matrix(c(phi1[i],phi2[i],1,0),ncol=2,nrow=2,byrow=T)
			A <- eigen(mat)
			a <- Re(A$values) # Real part
			b <- Im(A$values) # Imaginary part
			R <- Mod(A$values) # Modulus
			phi <- acos(a[1]/R[1])
			# Period, general solution:
			# (2*pi)/(acos(phi1/(2*sqrt(-phi2))))
			per[i] <- (2*pi)/phi

			if (plot == TRUE) {
				plot(x=c(-1.1,1.1), y=c(-1.1,1.1), type="n", xlab="Re", ylab="Im", font.lab=2, las=1)
				abline(v=0, lty=2)
				abline(h=0, lty=2)
				curve(sqrt(1-x^2), from=-1, to=1, add=T)
				curve(-sqrt(1-x^2), from=-1, to=1, add=T)
				segments(x0=0,y0=0,x1=a[1],y1=b[1], col=2, lty=2)
				points(a,b,col=2,pch=16)
			}
		}
	}

	if (method == "polyroot") {
		for (i in 1:n) {
		u <- polyroot(c(1,-phi1[i],-phi2[i]))
		per[i] <- 2*pi/Arg(u[1])
		}
	}

	out <- output(phi1=phi1, phi2=phi2, period=per, ipv=ar2.ipv2(phi1,phi2))
	out
}

# ar2.k
# Contour lines for periods (k-contours)

# The formula for the k-contours have been derived in two different ways:
# For method="royama", see p. 56-60 in Royama (1992). The solution is presented in eq. 2.17 (p. 60)
# For method="jenkins", see p. 338 in Cryer & Chan 2008, equation 13.5.9. Note that the two methods are not
# equivivalent, although they converge at high intrinsic process variance (i.e. as delayed density-dependence increases).
# Jenkins & Watts - contours show the contours for where the peak of the spectrum for an AR(2)-process is located.
# Royama's method is the correct method for obtaining parameters for e.g. simulations with arima.sim
# (since Jenkins & Watts' contours can not retrieve the input parameters, period and variance).

# Eq 2.17 in Royama, solved for a2:
# tan(w) = sqrt(abs((1+a1)^2 + 4*a2)) / (1+a1)
# tan(w) * (1+a1) = sqrt(abs((1+a1)^2 + 4*a2))
# tan(w)^2 * (1+a1)^2 = abs((1+a1)^2 + 4*a2)
# tan(w)^2 * (1+a1)^2 = -(1+a1)^2 - 4*a2
# a2 = (tan(w)^2 * (1+a1)^2 + (1+a1)^2) / -4
# a2 = -((1+a1)^2 / 4) * (tan(w)^2 + 1)
# set w = 2*pi/k, where k is the quasi-period
#' @export
ar2.k <- function(a1, k, kmethod="royama") {
	switch(kmethod,
		"royama" = (-0.25*(1+a1)^2) * (1 + tan(2*pi/k)^2),
		#"royama" = ar2.contour <- function(a1,k) -(1+a1)^2 / (4*cos(2*pi/k)^2) # Equivivalent!
		"jenkins" = -(1+a1) / (4*cos(2*pi/k) - (1+a1))
	)
}

#' @export
ar2.k2 <- function(phi1, k, kmethod="royama") {
	switch(kmethod,
		"royama" = (-0.25*phi1^2) * (1 + tan(2*pi/k)^2),
		#"royama" = ar2.contour <- function(a1,k) -phi1^2/(4*cos(2*pi/k)^2) # Equivivalent!
		"jenkins" = -phi1 / (4*cos(2*pi/k) - phi1)
	)
}

# ar2.k.inv
# see Example > 03 Plot > AR(2)_variance_contours_inverse_functions.r
#' @export
ar2.k.inv <- function(x, k) {
	# if k < 4, cos(2*pi/k) is negative
	# if k > 4, cos(2*pi/k) is positive
	# if k = 4, cos(2*pi/k is zero)
	#y <- -2*cos(2*pi/k)*sqrt(-x)
	y <- 2*cos(2*pi/k)*sqrt(-x) # The second solution
	y
}

# ar2.k.add
# Function that can be used to add extra k-contours after the ar2.plot has been created.
# Vectors can not be used, i.e. k = must be a single value of length 1.
# This is done by calling the curve-function.
# Can take optional arguments such as lwd, col, lty ...

# Dependencies: ar2.k
#' @export
ar2.k.add <- function(k, kmethod="royama", n=1001, output=list, ...) {

  if (k != 4) {
		x.max <- 2*cos(2*pi/k)
		z <- curve(ar2.k(a1=x-1,k=k,kmethod=kmethod),
			from=ifelse(k <= 4, x.max, 0), to=ifelse(k <=4, 0, x.max),
			ylim=c(-1,0), n=n, add=TRUE, ...)
	}

	if (k == 4) {
		lines(c(0,0), c(-1,0), ...)
		x <- rep(0,n)
		# ar2.k returns zero for all x=0 (a1=-1)
		z <- list(x=x, y=ar2.k(x-1,k,kmethod=kmethod))
	}

	invisible(output(a1=z$x - 1, a2=z$y, phi1=z$x, phi2=z$y))

}

# ar2.k.gen
# generate values of (1+a1) & a2 along a k-contour (k=quasi-period)
# method = "x" creates equally spaced values along the x (1+a1) - axis.
# method = "y" creates equally spaced values along the y (a2) - axis.
# from (defaults to -1) , to (defaults to 0):
# only valid for method="y", generate values only in the range specified by from,to arguments
# default option is the generate values in the entire range with solutions with complex roots,
# that is under the parabola in ar2.plot.

# ... = optional arguments must be given. Either length.out or n, which are
# used in an internal call to seq (used for generating sequences).
#' @export
ar2.k.gen <- function(k,method="y",from=NULL,to=NULL,plot=TRUE, ...) {

	if (!method %in% c("x","y")) stop("Method must be method='x' or method='y'")

	if (method == "y") {

		y.min <- -1
		y.max <- 0

		if (!is.null(from) && from > y.max) stop(paste("'from' must not be greater than", y.max))
		if (!is.null(to) && to < y.min) stop(paste("'to' must not be less than", y.min))

		if (is.null(from) & is.null(to)) y <- seq(y.min,y.max,...)
		if (!is.null(from) & !is.null(to)) {
			if (from < y.min) from <- y.min
			if (to > y.max) to <- y.max
			y <- seq(from,to,...)
		}
		if (!is.null(from) & is.null(to)) {
			if (from < y.min) from <- y.min
			y <- seq(from,y.max,...)
		}
		if (is.null(from) & !is.null(to)) {
			if (to > y.max) to <- y.max
			y <- seq(y.min,to,...)
		}

		x <- 2*cos(2*pi/k)*sqrt(-y)
	}

	if (method == "x") {

		if (k == 4) stop("x-method not useful for k = 4")

		if (k != 4) {
			if (k < 4) { x.min <- 2*cos(2*pi/k); x.max <- 0 }
			if (k > 4) { x.min <- 0; x.max <- 2*cos(2*pi/k) }

			if (!is.null(from) && from > x.max) stop(paste("'from' must not be greater than", x.max))
			if (!is.null(to) && to < x.min) stop(paste("'to' must not be less than", x.min))

			if (is.null(from) & is.null(to)) x <- seq(x.min,x.max,...)
			if (!is.null(from) & !is.null(to)) {
				if (from < x.min) from <- x.min
				if (to > x.max) to <- x.max
				x <- seq(from,to,...)
			}
			if (!is.null(from) & is.null(to)) {
				if (from < x.min) from <- x.min
				x <- seq(from,x.max,...)
			}
			if (is.null(from) & !is.null(to)) {
				if (to > x.max) to <- x.max
				x <- seq(x.min,to,...)
			}

			y <- -x^2/(4*cos(2*pi/k)^2)
		}
	}

	if (plot) {ar2.plot(k=k,v=NULL,kcol=2,klty=2); points(x,y,col=2,pch=16)}

	out <- list(a1=x-1,a2=y,phi1=x,phi2=y)
	invisible(out)
}

# ar2.k.rand

# Create a random sample of (1+a1) & a2-values for a given period length (k).
# n = number of pairs (a1 & a2) to generate

# The random sample is created based on a specified distribution for either x (1+a1) or y (a2).
# The corresponding x or y values are then calculated.
# method = c("x","y").
# If method="y" [default] sample from distribution (see below) along y-axis (a2)
# If method="x" [default] sample from distribution (see below) along x-axis (1+a1)

# ... spec = name of distribution + additional arguments to distribution function
# Examples: for normal distribution: spec="norm". Necessary extra arguments are mean & sd.
# Examples: for uniform distribution: spec="unif". Necessary extra arguments are min & max.

# a = lower truncation point
# b = upper truncation point

# Dependencies: trunc [rtrunc, dtrunc, qtrunc, ptrunc], ar2.plot.simple, ar2.k.add, subplot{TeachingDemos}
#' @export
ar2.k.rand <- function(k=5, n=50,
	method="y", a=NULL, b=NULL,
	output=list, plot=TRUE, xlims=c(-2,2), ylims=c(-1,0), ...) {

	if (method == "y") {
		# x = (1+a1)
		# y = a2
		if (is.null(a)) a <- -1
		if (is.null(b)) b <- 0
		y <- rtrunc(n=n, a=a, b=b, ...)
		x <- (2*sqrt(-y) * cos(2*pi/k))
	}

	if (method == "x") {
		if (k==4) stop("method = 'x' is not appropriate for k = 4, use method = 'y'")

		if (is.null(a)) {
			if (k > 4) a <- 0
			if (k < 4) a <- 2 * cos(2*pi/k)
		}
		if (is.null(b)) {
			if (k > 4) b <- 2 * cos(2*pi/k)
			if (k < 4) b <- 0
		}
		x <- rtrunc(n=n, a=a, b=b, ...)
		y <- -x^2/(4*cos(2*pi/k)^2)
	}

	if (plot) {
		ar2.plot.simple(main="k-contour")
		ar2.k.add(k=k,col=2)
		if (k == 4) points(rep(0,length(x)), y)
		if (k != 4) points(x, y)

		if (method=="y") x.str <- expression(beta[2])
		if (method=="x") x.str <- expression(1+beta[1])

		subplot(
			fun={
				curve(dtrunc(x=x, a=a, b=b, ...), from=a, to=b, xlim=c(a,b), n=1001,
					ylab="density", xlab=x.str)
				if (method=="x") lines(density(x),col=2)
				if (method=="y") lines(density(y),col=2)
				abline(v=c(a,b),col=2,lty=2)
			},
			x="top", hadj=-2, vadj=1.25)

	}

	if (k == 4) out <- output(a1=rep(-1,length(x)), a2=y, phi1=rep(0,length(x)), phi2=y, a=a, b=b)
	if (k != 4) out <- output(a1=x-1, a2=y, phi1=x, phi2=y, a=a, b=b)

	invisible(out)
}
