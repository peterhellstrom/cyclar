# Function for calculating measures of central tendency
#' @export
central <- function(x) {
	c1 <- mean(x)
	c2 <- exp(sum(log(x))/length(x))
	c3 <- 1/mean(1/x)
	c4 <- median(x)
	c5 <- range(x)

	cat("Arithmetic mean = ",round(c1,5),"\n")
	cat("Geomethric mean = ",round(c2,5),"\n")
	cat("Harmonic mean   = ",round(c3,5),"\n")
	cat("Median          = ",round(c4,5),"\n")
	cat("Range from",round(c5[1],5), "to", round(c5[2],5),"\n")
}

# cv
# x = a (time) series

cv <- function(x) sd(x) / abs(mean(x))

# covmat
# Generate a diagonal covariance matrix, i.e. independence between the columns of z,
# with non-zero off-diagonal terms
# Useful when simulating random correlated variables

# See intro to copulas from Mathworks
# http://www.mathworks.se/products/statistics/demos.html?file=/products/demos/shipping/stats/copulademo.html

# Dependencies: {MASS}
#' @export
covmat <- function(rho,sigma,n=NULL) {
	# require(MASS)
	if (length(sigma) == 1) var.x <- rep(sigma^2, n)
	if (length(sigma) > 1) var.x <- sigma^2
	Sigma <- rho * sqrt(var.x) %o% sqrt(var.x)
	diag(Sigma) <- var.x
	Sigma
}

#' @export
Square <- function(A){ # Test whether matrix A is square

    s <- sqrt(floor(length(A)))
    {if (!(length(A) == s^2)){
        print('Warning! Matrix is not square.',quote=FALSE)
        R <- FALSE}
     else{R <- TRUE}}
    R
}

# matrix.to.list
# x = a matrix
# margin = split matrix by rows (margin=1) or by columns (margin=2).
#' @export
matrix.to.list <- function(x, margin = 2) {
	if (!"matrix" %in% class(x)) stop("Input object x must be a matrix")
	if (margin == 1) {
		out <- lapply(seq_len(nrow(x)), function(i) x[i,])
		if (!is.null(rownames(x))) names(out) <- rownames(x)
	}
	if (margin == 2) {
		out <- lapply(seq_len(ncol(x)), function(i) x[,i])
		if (!is.null(colnames(x))) names(out) <- colnames(x)
	}
	out
}

# EXAMPLE
#n <- 1E5
#z <- matrix(runif(n), ncol = n / 1000, nrow = n / 100)
#dim(z)
#class(z)

#system.time(p0 <- matrix.to.list(z))
#system.time(p1 <- matrix.to.list(z, margin = 1))
#system.time(p2 <- matrix.to.list(z, margin = 2))
#' @export
matrix.to.vector <- function(x, margin=1) {
	if (margin == 1) x <- as.numeric(t(x))
	if (margin == 2) x <- as.numeric(x)
	x
}

#' @export
bugs.out <- function(coda, output="both") {

	est0 <- summary(coda) # outputs a list

	if (output=="both") out <- est0
	if (output=="statistics") out <- est0$statistics # Extract mean & sd's.
	if (output=="quantiles") out <- est0$quantiles # Extract quantiles.

	out
}

#' @export
niche.mat <- function(x) {
	x <- t(x)
	prop.x <- x / rowSums(x)
	out <- matrix(NA, nrow=nrow(prop.x), ncol=nrow(prop.x))
	for (i in 1:nrow(prop.x)) {
		for (j in 1:nrow(prop.x)) {
			zz <- prop.x[c(i,j),]
			out[i,j] <- sum(apply(zz,2,prod)) / sqrt(prod(rowSums(zz^2)))
		}
	}
	out
}
