# ar2.plot
# Plot AR(2)-space
# k = vector with contour lines for quasi-period length in the region a1^2 + 4*a2 < 0
# v = vector with contour lines for intrinsic process variance
# kcontours = TRUE/FALSE. If TRUE, the contours specified by argument k are drawn
# kcol, klty & klwd are col, lty & lwd parameters for k-contours.
# kcol.lab = color for k-labels
# klab = TRUE/FALSE. TRUE = k-labels are added to the graph
# kmethod = c("royama","jenkins"). Two different ways of drawing the k-contours.
# if kmethod = "royama", the caclulations follows p. 56-60 in Royama (1992).
# if kmethod = "jenkins", the contours are drawn according to Jenkins & Watts (1968, cited by Cryer & Chan [2008], p. 338).
# The latter option is based on the peak frequency of the spectral density of an AR(2)-process.
# vcontours = TRUE/FALSE. If TRUE, the contours specifided by argument v are drawn.
# vcol, vlty & vlwd are col, lty & lwd parameters for v-contours.
# vcol.lab = color for v-labels
# vregion = c("parabola","triangle"). If "parabola", v-contours are drawn only under the parabola.
# If "triangle" v-contours are drawn in the triangle delineating the statinary regions for an AR(2)-process.
# vlab = TRUE/FALSE. If labels for v should be added at position vlabels.pos
# vlab.pos = position of v contours (drawn on a period contour, defaults to k=2)
# vlab.both = TRUE/FALSE. If TRUE, v-contour labels are added at both the left and right side of the parabola
# (at intersection between parabola and v-contour).
# cex.lab = size of labels
# cex.axis = size of axes notation
# cex.contours = size of contour labels

# par.name = name of greek letter (within quotes) to use for x- and y-axes labels.
# defaults to par.name="beta", can be changed to any greek letter available for expressions [see ?plotmath].
# e.g. par.name="Omega" or par.name="alpha".
# text.lab = TRUE/FALSE. Two different options for x- and y-labels:
# If TRUE, xlab = 1 + beta[1] (Direct density dependence), ylab = beta[2] (Delayed density dependence).
# If FALSE: xlab = 1 + beta[1], ylab = beta[2]
# x-labels and y-labels can be changed to other expressions, use the arguments xlab and ylab.
# If you want to change the default x-label, you could for instance use xlab=expression(phi[1]).
# formula = TRUE/FALSE. If TRUE, the equation for an AR(2)-process is written on the plot at position (1,1)
# formula.pos = vector of length two specifying position of formula on graph. Default is c(x=1, y=1).
# triangle = TRUE/FALSE. If TRUE, the region where the AR(2) is stationary is outlined with a "boundary triangle".
# xlim = extent of x-axis
# ylim = extent of y-axis

# Dependencies: ar2.parms, ar2.ipv, ar2.v.add, ar2.parabola, ar2.k, ar2.k.add, ar2.v.inv
#' @export
ar2.plot <- function(
        k = c(3, 3.5, 4, 4.5, 5, 6, 7), v = c(2, 4, 10),
				kcontours = TRUE, kcol = 2, klty = 1, klwd = 1, klab = TRUE, kcol.lab = 1, kmethod = "royama",
        vcontours = TRUE, vregion = "parabola", vcol = 4, vlty = 2, vlwd = 1,
        vlab = TRUE, vlab.pos = 2, vlab.both = FALSE, vcol.lab = "darkgrey",
        cex.contours = 1.2, cex.lab = 1.2, cex.axis = 1.2,
        xlab = expression(paste(1 + beta[1])), ylab = expression(paste(beta[2])),
        par.name = "beta", text.lab = TRUE, formula = FALSE, formula.pos = c(1,1),
        triangle = TRUE, xlim = c(-2,2), ylim = c(-1,1),
				main = "") {

	if (!is.null(v)) if (min(v) <= 1) stop("Variance equal or less than unity not allowed!")

	if (text.lab) {
		xlab <- expression(paste(1 + beta[1], " (Direct density-dependence)"))
		ylab <- expression(paste(beta[2], " (Delayed density-dependence)"))
		if (any(gregexpr(par.name,xlab)[[1]] == -1)) xlabs <- parse(text = gsub("beta", par.name, xlab))
		if (any(gregexpr(par.name,ylab)[[1]] == -1)) ylabs <- parse(text = gsub("beta", par.name, ylab))
	} else {
		if (any(gregexpr(par.name,xlab)[[1]] == -1)) xlab <- parse(text = gsub("beta", par.name, xlab))
		if (any(gregexpr(par.name,ylab)[[1]] == -1)) ylab <- parse(text = gsub("beta", par.name, ylab))
	}

	# Basic plot
		op <- par(mar = c(4.5, 4.5, 1, 2))
		plot(x = xlim, y = ylim, type = "n", bty = "l",
		xlab = xlab,
		ylab = ylab,
		main = main,
		font.lab = 2, las = 1, cex.lab = cex.lab, cex.axis = cex.axis)

		if (triangle) {
			points(c(-2, 0, 2, -2), c(-1, 1, -1, -1), type = "l") # Outer triangle
			lines(c(0, 0), c(-1, 1), col = 1) # Vertical line at k=4
		} else {
			points(c(-2, 2), c(-1, -1), type = "l")
			lines(c(0, 0), c(-1, 0), col = 1) # Vertical line at k=4
		}

		curve(ar2.parabola(a1 = x - 1), from = -2, to = 2, add = TRUE)
		par(op)

	# Draw k- & v-contours & associated labels
	# k-contours (period)
	if (!is.null(k)) {

		k <- unique(k)

		if (kcontours) {
			for (i in 1:length(k)) ar2.k.add(k = k[i], kmethod = kmethod, col =  kcol,lty = klty, lwd = klwd)
		}

		if (klab) {
			x <- 2*sqrt(-min(ylim)) * cos(2*pi/k)
			y <- -x^2 / (4*cos(2*pi/k)^2)
			points(x = x, y = y, pch = 15, col = "white", cex = cex.contours*2)
			text(x = x, y = y, labels = round(k,2), font = 2, cex = cex.contours, col = kcol.lab)
		}
	}

	# v-contour (intrinsic process variance)
	if (!is.null(v)) {

		v <- unique(v)

		if (vcontours) {
			for (i in 1:length(v)) ar2.v.add(v = v[i], vregion = vregion, col = vcol, lty = vlty, lwd = vlwd)
		}

		if (vlab) {
			ipv.lab.pos <- ar2.parms(k = rep(vlab.pos,length(v)), v = v, plot = FALSE)
			points(x = 1 + ipv.lab.pos$a1, y = ipv.lab.pos$a2, pch = 15, col = "white", cex = cex.contours*2)
			text(x = 1 + ipv.lab.pos$a1, y = ipv.lab.pos$a2, labels = round(v,2), font = 2, cex = cex.contours, col = vcol.lab)
			if (vlab.both) {
				points(x = -ipv.lab.pos$a1 - 1, y = ipv.lab.pos$a2, pch = 15, col = "white", cex = cex.contours*2)
				text(x = -ipv.lab.pos$a1 - 1, y = ipv.lab.pos$a2, labels = round(v,2), font = 2, cex = cex.contours, col = vcol.lab)
			}
		}
	}

	# Add equation for AR2-model
	if (formula) {
		expr <- expression(x[t] == beta[0] + (1 + beta[1])*x[t-1] + beta[2]*x[t-2] + epsilon[t])
		#expr <- expression((x[t]-mu) == (1 + beta[1])*(x[t-1]-mu) + beta[2]*(x[t-2]-mu) + epsilon[t]), cex = cex.axis)
		if (any(gregexpr(par.name,expr)[[1]] == -1)) expr <- parse(text = gsub("beta", par.name, expr))
		text(formula.pos[1], formula.pos[2], expr, cex = cex.axis)
	}

}

# ar2.parabola

# Parabola = a1^2 + 4*a2 = 0, solved for a2:
#' @export
ar2.parabola <- function(a1) -0.25*(1+a1)^2

#' @export
ar2.parabola2 <- function(phi1) -0.25*phi^2

# ar2.plot.simple
# A basic version of the ar2.plot.
# A "skeleton" function used as a subplot in ar2.acf, but can be used independently as well.
# No contours are drawn with this function
# It is possible to enter points by supplying vectors with a1(x) and a2(y)-values.
# Other input arguments (where available) follows the main function ar2.plot.
#' @export
ar2.plot.simple <- function(xlims=c(-2,2),ylims=c(-1,1),a1=NULL,a2=NULL,
  xlab=expression(1 + beta[1]), ylab=expression(beta[2]), par.name="beta", ...) {

  if (any(gregexpr(par.name,xlab)[[1]]==-1)) xlab <- parse(text=gsub("beta", par.name, xlab))
	if (any(gregexpr(par.name,ylab)[[1]]==-1)) ylab <- parse(text=gsub("beta", par.name, ylab))

	plot(c(xlims[1],xlims[2]), c(ylims[1],ylims[2]), type="n", bty="l", xlab=xlab, ylab=ylab, ...)

	# Triangle
	points(c(-2,0,2,-2), c(-1,1,-1,-1), type="l") # Outer triangle
	lines(c(0,0), c(-1,1)) # Vertical line at k=4

	# Parabola
	curve(-0.25*x^2, from=-2, to=2, add=TRUE) # Parabola = (1 + a1)^2 + 4*a2 = 0, solve for a2
	points((1+a1),a2,pch=16,col=2)
}

#' @export
ar2.plot.simple2 <- function(xlims=c(-2,2),ylims=c(-1,1),phi1=NULL,phi2=NULL,
  xlab=expression(phi[1]), ylab=expression(phi[2]), par.name="phi", ...) {

  if (any(gregexpr(par.name,xlab)[[1]]==-1)) xlab <- parse(text=gsub("phi", par.name, xlab))
	if (any(gregexpr(par.name,ylab)[[1]]==-1)) ylab <- parse(text=gsub("phi", par.name, ylab))

	plot(c(xlims[1],xlims[2]), c(ylims[1],ylims[2]), type="n", bty="l", xlab=xlab, ylab=ylab, ...)

	# Triangle
	points(c(-2,0,2,-2), c(-1,1,-1,-1), type="l") # Outer triangle
	lines(c(0,0), c(-1,1)) # Vertical line at k=4

	# Parabola
	curve(-0.25*x^2, from=-2, to=2, add=TRUE) # Parabola = phi1^2 + 4*a2 = 0, solve for a2
	points(phi1,phi2,pch=16,col=2)
}

# ar2.arrows
# Add time-varying parameter estimates as arrows to an existing ar2.plot.
# (see Bierman et al 2006 Am Nat for an example).
# a1 and a2 are vectors with parameters in ar(2)-space.
# col, lwd, angle & length are arguments passed to the arrows-function.
# end.only = FALSE/TRUE. If TRUE, only the final a1/a2 pair is drawn with an arrowhead,
# the previous points are added as lines. If FALSE, an arrowhead is drawn at each
# a1/a2 pair.
#' @export
ar2.arrows <- function(a1, a2, col=2, lwd=1, angle=30, length=0.075, end.only=FALSE) {
	n1 <- length(a1); n2 <- length(a2)
	if (n1 != n2) stop("Vectors do no match, unequal lengths")
	n <- n1
	if (n < 2) stop("Number of data points must exceed 2")

	if (end.only == FALSE) {
		arrows((1+a1)[-n],a2[-n],(1+a1)[-1],a2[-1],col=col,lwd=lwd,angle=angle,length=length)
	}
	if (end.only == TRUE) {
		a1 <- as.numeric(a1)
		a2 <- as.numeric(a2)
		lines((1+a1),a2,col=col,lwd=lwd,lty=1)
		arrows((1+a1)[n-1], a2[n-1], (1+a1)[n], a2[n], col=col,lwd=lwd,angle=angle,length=length)
		points((1+a1)[1], a2[1],col=col,pch=16,cex=1)
	}

}

#' @export
ar2.arrows2 <- function(phi, phi2, col=2, lwd=1, angle=30, length=0.075, end.only=FALSE) {
	n1 <- length(phi1); n2 <- length(phi2)
	if (n1 != n2) stop("Vectors do no match, unequal lengths")
	n <- n1
	if (n < 2) stop("Number of data points must exceed 2")

	if (end.only == FALSE) {
		arrows(phi1[-n],phi2[-n],phi1[-1],phi2[-1],col=col,lwd=lwd,angle=angle,length=length)
	}
	if (end.only == TRUE) {
		phi1 <- as.numeric(phi1)
		phi2 <- as.numeric(phi2)
		lines(phi1,phi2,col=col,lwd=lwd,lty=1)
		arrows(phi1[n-1], phi2[n-1], phi1[n], phi2[n], col=col,lwd=lwd,angle=angle,length=length)
		points(phi1[1], phi2[1],col=col,pch=16,cex=1)
	}

}

# Pink Floyd meets autoregressive models...
#' @export
ar2.dsotm <- function(a1=NULL, a2=NULL, n=100, add.series=TRUE, col.series="darkgreen", main=NULL,
	scatter.smooth=FALSE,...) {

	# Open graph window, set background to black
	op <- par(list(mar=c(1,1,1,1),bg="black"))

	# Prism / triangle
	xmax <- 2.5 # Limits for plot

	if (!scatter.smooth)
		plot(x=c(-xmax,xmax),y=c(-xmax/2,xmax/2),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

	if (scatter.smooth) {
		x1 <- seq(-1,0, by=0.01)
		x2 <- seq(0,1, by=0.01)
		x3 <- seq(-1,1, by=0.01)
		x <- c(x1, x2, x3)

		est1 <- 1 + 2*x1
		est2 <- 1- 2*x2
		est3 <- rep(-1, length(x3))
		est <- c(est1, est2, est3)

		z <- cbind(x,est)
		Lab.palette <- colorRampPalette(c("black", "grey", "white"), space = "Lab")
		smoothScatter(z, colramp = Lab.palette, nrpoints=0, ...,
			xlim=c(-xmax,xmax), ylim=c(-xmax/2,xmax/2), xaxt="n",yaxt="n",bty="n")
	}

	if (is.null(main)) title(main=list("Why Pink Floyd?", col="white"))
	if (!is.null(main)) title(main=list(main, col="white"))

	# Add triangle
	# points(c(-1,0,1,-1), c(-1,1,-1,-1), type="l")
	curve(-x^2, from=-1,to=1,add=T,n=1001,lty=2,col="darkgrey")

	# Ray of light from left
	# abline(v=0, lty=3)
	# abline(h=0, lty=3)
	# abline(h=1, lty=3, col=4)

	a <- 0.2
	b <- 0.35
	# abline(a,b,col=2)
	# abline(a,-b,col=4)
	# abline(v=(a-1)/b, lty=3)
	# abline(v=-(a-1)/b, lty=3)
	# points(0,a,col=4,pch=15)
	segments(x0=-xmax, y0= a + b*-xmax, x1=(a-1)/(2-b), y1=1+2*(a-1)/(2-b), col="white", lwd=3)
	# abline(1,2,col=2,lty=2)
	# abline(1,-2,col=4,lty=2)
	# points((a-1)/(2-b), 1+2*((a-1)/(2-b)), pch=16, col=2)

	# "Rainbow" at right side: polygons

	# Find equation for lines:

	xv <- (a-1) / b # x-position where all lines interesect (negative)
	dy1 <- 0.35 # difference in y between first and seventh polygon
	z1 <- (1-a)/(b+2) # start x-position for polygons
	z2 <- xmax

	# Calculate slopes for top (7) and bottom line (1)
	b7 <- 2*z1 / (xv-z1)
	b1 <- (2*z1 + dy1) / (xv-z1)
	bv <- seq(b7,b1,length.out=7) # generate a sequence of slopes
	av <- 1 - bv*xv # calculate intercepts

	# Coordinates for polygons
	x.start <- (1-av) / (bv + 2)
	x.stop <- rep(xmax, length(x.start))
	y.start <- av + bv*x.start
	y.stop <- av + bv*xmax

	cols <- c("red","orange","yellow","green","cyan","purple")
	# Polygons
	for (i in 1:length(bv)) {
		polygon(
			x = c(x.start[i], xmax, xmax, x.start[i+1], x.start[i]),
			y = c(y.start[i], y.stop[i], y.stop[i+1], y.start[i+1], y.start[i]),
			col=cols[i],border=cols[i])
	}

	# Add triangle
	if (!scatter.smooth) points(c(-1,0,1,-1), c(-1,1,-1,-1), type="l",col="white")

	if (add.series) {
			# Add a time series (lynx) in the "rainbow" section:
			# Create a trend in a time series
			# 1) Negative slope
			# 2) Increasing variance with time (so it fits "rainbow section")
			# 2) Not implemented (yet?)

		if (is.null(a1) | is.null(a2)) series <- ts(log(lynx),start=1)
		if (!is.null(a1) && !is.null(a2)) series <- ar2.sim2(n=n, mu=0, a1=a1, a2=a2)

		x.coords <- seq(from=x.start[4], to=xmax, length.out=length(series))
		trend <- ts(av[4] + bv[4]*x.coords,start=1)
		# scale series
		scaled.series <- scale(series)
		# Check this section so that series exactly fits in the rainbow region
		scaled.series <- scaled.series / (diff(range(scaled.series)) * (1/(sqrt(2)*dy1)))
		# scaled.series <- scaled.series / (diff(range(scaled.series)) * (1/dy1))

		# mean(scaled.series); sd(scaled.series); diff(range(scaled.series))
		# add trend to series
		y <- scaled.series + trend
		#plot(y); #points(trend)
		# lm(y ~ x.coords)

		points(x.coords, y, type="l", lwd=2, col=col.series)

		if (!is.null(a1) & !is.null(a2)) {
			points((1+a1)/2, a2, pch=16, col="white")
		}

	}

	par(op)
}
