# ar2.ipv
# Intrinsic process variance
# Note that the first-order component is (1+a1) [as in Royama 1992],
# which is different from the references cited in ar2.ipv.method,
# where phi1 = 1+a1.
#' @export
ar2.ipv <- function(a1,a2,sd=1) as.numeric((sd^2) * (1-a2) / ((1-(1+a1)-a2)*(1-a2+(1+a1))*(1+a2)))

#' @export
ar2.ipv2 <- function(phi1,phi2,sd=1) as.numeric((sd^2) * (1-phi2) / ((1-phi1-phi2)*(1-phi2+phi1)*(1+phi2)))

# A compilation of various sources found in the literature.
# All methods return the same result.

#' @export
ar2.ipv.method <- function(a1,a2,sd=1,method="Kendall-Ord") {

	sigma2 <- sd^2

	as.numeric(switch(method,
		#
		'rho' = {rho1 <- (1+a1)/(1-a2); rho2 <- (1+a1)*rho1 + a2; as.numeric(sigma2 / (1 - (1+a1)*rho1 - a2*rho2))},
		# Kendall & Ord (1990), p. 59
		'Kendall-Ord' = sigma2 * (1-a2) / ((1-(1+a1)-a2)*(1-a2+(1+a1))*(1+a2)),
		# Hamilton (1994) Time Series Analysis: p. 58, equation 3.4.30
		'Hamilton' = sigma2 * (1-a2) / ((1+a2) * ((1-a2)^2 - (1+a1)^2)),
		# Cryer & Chan (2008) Time Series Analysis - with Applications in R (2nd ed): p. 75, equation 4.3.20
		'Cryer-Chan' = ((1-a2)/(1+a2)) * (1 / ((1-a2)^2 - (1+a1)^2)) * sigma2,
		# Royama (1992) Analytical Population Dynamics: p. 117, equation 3.34
		'Royama' = sigma2 * (1-a2) / ((1+a2)*(a1+a2)*(a2-a1-2)),
		# Priestley (1981) Priestley 1981 Spectral Analysis and Time Series - Vol I: p. 128, equation 3.5.34
		# Note that Priestley have defined the a1-a2 space so that a2 is positive.
		'Priestley' = {a2 <- -a2; sigma2 * (1+a2) / ((1-a2)*(1-(1+a1)+a2)*(1+(1+a1)+a2))}
	))
}

# ar2.v.gen
# generate values of (1+a1) & a2 along a v-contour (v = intrisic process variance, see ar2.ipv)
# method = "x" is not implemented.
# method = "y" creates equally spaced values along the y (a2) - axis.
# vregion c("parabola","triangle"). If "parabola", contours are drawn in the region a1^2 + 4*a2 < 0.
# If "triangle", the contours are drawn also above the parabola.
# side = c("pos","neg","both"). If "pos", only positive (1+a1) values are created,
# if "neg", only negative (1+a1) are created.
# if side = "both", [default] (1+a1) values are generated first for postive (1+a1) values, then
# duplicated at negative (1+a1) values.

# ... = optional arguments must be given. Either length.out or n, which are
# used in an internal call to seq (used for generating sequences).

#' @export
ar2.v.gen <- function(v,method="y",vregion="parabola", from=NULL, to=NULL, side="both", plot=TRUE, ...) {

	if (!method %in% c("y")) stop("Method must be method='y', method='x' is not implemented")

	if (method == "y") {
		y.min <- -sqrt(1-(1/v))
		if (vregion=="parabola") y.max <- (sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3)-1/(3*v*(sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3))-1
		if (vregion=="triangle") y.max <- -y.min

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

		x2 <- y^2 - 2*y + 1 - ((1 - y)/((1 + y)*v))
		x <- sqrt(abs(x2))
	}

	if (side=="pos") out <- list(a1=x-1,a2=y,phi1=x,phi2=y)
	if (side=="neg") out <- list(a1=-(x-1),a2=-y,phi1=-x,phi2=-y)
	if (side=="both") {
		x <- c(rev(-x),x)
		y <- c(rev(y),y)
		out <- list(a1=x-1,a2=y,phi1=x,phi2=y)
	}

	if (plot) {ar2.plot(k=NULL,v=v,vcol=4,vlty=2); points(x,y,col=2,pch=16)}

	invisible(out)
}

# ar2.v.rand

# Create a random sample of (1+a1) & a2-values for a given period length (v).
# v = vector (scalar) with the desired intrinsic process variance v.
# See function ar2.ipv for how v depends on parameters a1 & a2.

# n = number of pairs (a1 & a2) to generate

# vregion c("parabola","triangle"). If "parabola", random values are generated in the region a1^2 + 4*a2 < 0.
# If "triangle", random points are generated also above the parabola.

# The random sample is created based on a specified distribution y (a2).
# The corresponding x are then calculated.

# ... spec = name of distribution + additional arguments to distribution function
# Examples: for normal distribution: spec="norm". Necessary extra arguments are mean & sd.
# Examples: for uniform distribution: spec="unif". Necessary extra arguments are min & max.

# a = lower truncation point
# b = upper truncation point

# Dependencies: ar2.plot.simple, ar2.v.inv, ar2.v.add, subplot{TeachingDemos}, trunc [rtrunc, dtrunc, qtrunc, ptrunc].
#' @export
ar2.v.rand <- function(v, n=50, vregion="parabola", side="both",
	plot=TRUE, output=list, ...) {

	if (vregion == "triangle") {
		a <- -sqrt(1-(1/v))
		b <- -a
		y <- rtrunc(n=n, a=a, b=b, ...)
	}
	if (vregion=="parabola") {
		a <- -sqrt(1-(1/v))
		b <- (sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3)-1/(3*v*(sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3))-1
		y <- rtrunc(n=n, a=a, b=b, ...)
	}

	x <- ar2.v.inv(x=y, v=v, na.rm=TRUE) # only positive values

	if (side == "both") {
		inds.pos <- sample(1:n, n/2, replace=FALSE)
		x.pos <- x[1:n %in% inds.pos]
		x.neg <- x[!1:n %in% inds.pos]
		x <- c(x.pos, -x.neg)
		y <- c(y[1:n %in% inds.pos], y[!1:n %in% inds.pos])
		out <- output(a1=x-1, a2=y, phi1=x, phi2=y)
	}

	if (side == "pos") out <- output(a1=x-1, a2=y, phi1=x, phi2=y)
	if (side == "neg") out <- output(a1=-(x-1), a2=y, phi1=-x, phi2=y)

	if (plot) {
		ar2.plot.simple(main="v-contour")
		ar2.v.add(v=v, vregion=vregion, col=2)
		points(out$phi1,out$phi2)

		subplot(
			fun={curve(dtrunc(x=x, a=a, b=b, ...), from=a, to=b, xlim=c(a,b), n=1001,
				ylab="density", xlab=expression(1 + beta[2]))
				lines(density(y),col=2)
				abline(v=c(a,b),col=2,lty=2)
			},
	     x="top", hadj=-2, vadj=1.25)
    }

	invisible(out)
}

# ar2.v.add
# Add a variance-contour (specified by v) to an existing ar2.plot.
# Vectors can not be used, i.e. k = must be a single value of length 1.
# Can take optional arguments such as lwd, col, lty ...

# Dependencies: ar2.v.inv
#' @export
ar2.v.add <- function(v, vregion="parabola", n=1001, output=list, ...) {
	# y = a2, x = (1+a1)
	if (vregion == "triangle") {
		y.min <- -sqrt(1-(1/v))
		y.max <- -y.min
		yv <- seq(from=y.min, to=y.max, length.out=n)
	}
	if (vregion=="parabola") {
		y.min <- -sqrt(1-(1/v))
		y.max <- (sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3)-1/(3*v*(sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3))-1
		yv <- seq(from=y.min, to=y.max, length.out=n)
	}

	xv <- ar2.v.inv(x=yv,v=v,na.rm=TRUE)
	# xv only contains positive (1+a1)-values, but the curve is symmetric.
	points(xv,yv,type="l",...)
	points(-xv,yv,type="l",...) # plot negative (1+a1) values.
	invisible(output(a1=c(-xv-1,xv-1), a2=c(yv,yv), x=c(-xv,xv), y=c(yv,yv)))
}


# ar2.v.inv
# see Example > 03 Plot > AR(2)_variance_contours_inverse_functions.r
#' @export
ar2.v.inv <- function(x, v, na.rm=TRUE, out="y") {

	if (any(out %in% c("y","both")) == FALSE) stop("out must be one of = c('both','y')")

	#y2 <- (-1 + x) * (1 + v*(-1 + x^2)) / (v*(1+x))
	y2 <- x^2 - 2*x + 1 - ((1 - x)/((1 + x)*v))

	inds <- which(sign(y2)==-1)

	if (na.rm == FALSE) y <- sqrt(y2)
	else if (na.rm == TRUE) y <- sqrt(abs(y2))

	if (out == "both") list(y2=y2, y=y, NAs=inds, y2NAs=y2[inds])
	else if (out == "y") y
}
