# ar2.parms
# Find parameters a1 & a2 that gives a certain period length k, given the intrinsic process variance v
# Uses uniroot to find the a1(x)-values. The function for which the root is found is specified
# in ar2.parms.root.
# It is possible to enter vectors (of equal length for k and v).

# Dependencies: ar2.parms.root, ar2.plot, matrix.to.list

# Equivivalent to older function ar.parms.old.
#' @export
ar2.parms <- function(k, v, plot = FALSE, list = TRUE, ...) {

	if (identical(length(k), length(v)) == FALSE) stop("Input vectors a1 and a2 must be of equal length")
	if (any(k < 2)) stop("k < 2 is not possible") # Why?
	n <- length(k)

	out <- t(sapply(1:n, function(i) ar2.parms.root(k[i], v[i], output = c)))

	if (plot == TRUE) {
		ar2.plot(k = unique(k), v = unique(v), ...)
		points(1 + out[,"a1"], out[,"a2"], pch = 16, col = 2)
	}

	if (list == TRUE) {
		out <- matrix.to.list(out, margin = 2)
	}

	out
}

# ar2.parms.old
# Find parameters a1 & a2 that gives a certain period length k, given the intrinsic process variance v
# Uses optimization to find the a1(x)-values.
# It is possible to enter vectors (of equal length for k and v).

# Dependencies: ar2.plot, ar2.ipv
#' @export
ar2.parms.old <- function(k,v,kmethod="royama",plot=FALSE,output=cbind, ...) {

	if (identical(length(k),length(v)) == FALSE) stop("Input vectors a1 and a2 must be of equal length")
  if (any(k < 2)) stop("k < 2 is not possible")

	n <- length(k)

	a1.out <- numeric(n)
	a2.out <- numeric(n)

	# Functions to use with optimize().
	# optimize() finds a max or min within the specified interval, in this case a minimum.
	# The optimization function is defined as:
	# 1) Set x = (1+a1)
	# 2) Given x, calculate a2.
	# 3) Use this set of (1+a1) & a2 to calculate the intrinsic process variance.
	# 4) Step 3 is subject to a constraint, so that optimize finds the set of (1+a1) & a2
	# that minimizes the squared difference between the desired intrinsic process variance v,
	# and the ipv given the current set of (1+a1) & a2.
	# Both the k-contours given by Royama and by Jenkins & Watts can be used.
	# However, if you plug the results from this function into ar2.period, you'll find that
	# only kmethod="royama" retrieves the input parameters (method="jenkins" & method="royama" converge
	# only at high ipv-values).

	if (kmethod == "royama") {
		f <- function(x,k,v) {
			a2 <- -x^2/(4*cos(2*pi/k)^2)
			(v - ar2.ipv(a1=x-1, a2=a2))^2
		}
	}

	if (kmethod == "jenkins") {
		f <- function(x,k,v) {
		a2 <- -x / (4*cos(2*pi*(1/k)) - x)
		(v - ar2.ipv(a1=x-1, a2=a2))^2
		}
	}

	for (i in 1:n) {

		k.temp <- k[i]
		v.temp <- v[i]

		if (k.temp != 4) {

			# Do not search for minima of function f below a2=-1
			# The (1+a1) value (given k) where a2=-1 is found by
			#(set the equation for the k-contour [a2 = -(1+a1)^2 / 4*cos(2*pi/k)^2 equal to zero and solve for x):
			x.lim <- 2 * cos(2*pi/k.temp)
			if (k.temp < 4) interval.lims <- c(x.lim,0) # Left side of parabola
			if (k.temp > 4) interval.lims <- c(0,x.lim) # Right side of parabola

			# Optimize finds (1+a1), so subtract 1 in order to get a1
			a1 <- optimize(f=f, interval=interval.lims, tol=1e-15, k=k.temp, v=v.temp)$minimum - 1
			# Calculate the corresponding a2 value
			if (kmethod == "royama") a2 <- -0.25*(1+a1)^2 / cos(2*pi/k.temp)^2
			if (kmethod == "jenkins") a2 <- -(1+a1) / (4*cos(2*pi*(1/k.temp)) - (1+a1))
			# Send to output vectors
			a1.out[i] <- a1
			a2.out[i] <- a2
		}
		# If quasi-period (k) = 4, it is not necessary to use optimize, since (1+a1)=-1 ==> a1=0
		# The a2-value at (1+a1=0 can be found by
		# setting x (1+a1) in the expression for ipv equal to zero and then solve for a2 (=y):
		# (1-y) / ((1-x-y)*(1-y+x)*(1+y)) = v
		# (1-y) / ((1-y)*(1-y)*(1+y)) = v
		# 1 / ((1-y)*(1+y)) = v
		# 1 / (1-y^2) = v
		# y^2 = (v-1)/v
		# y = +/- sqrt(1 - (1/v))
		else if (k.temp == 4) {
			a1 <- -1
			a2 <- -sqrt(1 - (1/v.temp))
			a1.out[i] <- a1
			a2.out[i] <- a2
		}
	}

	if (plot == TRUE) {
		ar2.plot(k=unique(k), v=unique(v), kmethod=kmethod, ...)
		points(1 + a1.out, a2.out, pch=16, col=2)
	}

	output(a1=a1.out, a2=a2.out, phi1=a1.out+1, phi2=a2.out)
}

# ar2.parms.root
# See AR2_find_parameters_uniroot.r in the Examples folder 04 Parameters.
# k = (quasi-)period
# v = intrinsic process variance
#' @export
ar2.parms.root <- function(k,v,output=list,tol=1e-15) {
	if (k == 4) {
		x <- -sqrt(1-(1/v))
		y <- 2*cos(2*pi/k)*sqrt(-x)
	}
	if (k != 4) {
		# Use inverse functions
		f.root <- function(x) x^2 - 2*x + 1 - ((1-x)/((1+x)*v)) + 4*cos(2*pi/k)^2*x
		x.min <- -sqrt(1-(1/v))
		x <- uniroot(f.root, interval=c(x.min,0), tol=tol)$root
		y <- 2*cos(2*pi/k)*sqrt(-x)
	}

	output(a1=y-1,a2=x,phi1=y,phi2=x)
}
