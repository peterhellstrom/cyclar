# Tests for randomness

# 3 different tests
# turning.point.test (tests for positive/negative correlations at lag 1)
# rank.test (tests for linear trend)
# difference.sign.test (tests for trend. [But a periodic sequence can pass this test])

# See chapter 9 in Brockwell & Davis, p. 312-313

# x = a time series
# alpha = a number between 0 and 1 which indicate the level of the test.

# These function are available in {hacks}, and code is taken from that library.
# Earlier versions of the functions in cyclar was adopted from code by Kyrre Lekve.

# library(hacks)
# Code also found at:
# http://www.uta.fi/~al18853/aika1.html
# turning point test
# rank.test (slightly changed calculation of STATISTIC, included here):
#' @export
turning.point.test <- function(x) {
	DNAME <- deparse(substitute(x))
	n <- length(x)
	METHOD <- "Turning point test"
	X <- embed(x,3)
	STATISTIC <- sum((X[,2] > X[,1] & X[,2] > X[,3])|(X[,2] < X[,1] & X[,2] < X[,3]))
	mu <- 2*(n-2)/3
	sigma2 <- (16*n-29)/90
	PVAL <- 2 * (1-pnorm(abs(STATISTIC-mu)/sqrt(sigma2)))
	PARAMETER <- c(mu,sigma2)
	names(STATISTIC) <- "normal"
	names(PARAMETER) <- c("mu", "sigma2")
	structure(list(statistic = STATISTIC, parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME),
        class = "htest")
}

#' @export
rank.test <- function(x) {
    DNAME <- deparse(substitute(x))
    n <- length(x)
    METHOD <- "Rank Test (NULL:iid)"
    STATISTIC <- 0
    # for (i in 1:(n - 1)) for (j in i:n) if (x[j] > x[i]) STATISTIC <- STATISTIC + 1
	STATISTIC <- sum(outer(x,x,"<")[outer(1:n,1:n,"<")])
	mu <- n * (n - 1)/4
    sigma2 <- n * (n - 1) * (2 * n + 5)/72
    PVAL <- 2 * (1 - pnorm(abs(STATISTIC - mu)/sqrt(sigma2)))
    PARAMETER <- c(mu, sigma2)
    names(STATISTIC) <- "normal"
    names(PARAMETER) <- c("mu", "sigma2")
    structure(list(statistic = STATISTIC, parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME),
        class = "htest")
}

#' @export
difference.sign.test <- function(x) {
    DNAME <- deparse(substitute(x))
    n <- length(x)
    METHOD <- "Difference-Sign Test (NULL:iid)"
    X <- embed(x, 2)
    STATISTIC <- sum(X[, 2] < X[, 1]) # sign changed compared to {hacks}
	# STATISTIC is in this case equal to:
	# STATISTIC <- length(which(diff(x) > 0))
    mu <- (n - 1)/2
    sigma2 <- (n + 1)/12
    PVAL <- 2 * (1 - pnorm(abs(STATISTIC - mu)/sqrt(sigma2)))
    PARAMETER <- c(mu, sigma2)
    names(STATISTIC) <- "normal"
    names(PARAMETER) <- c("mu", "sigma2")
    structure(list(statistic = STATISTIC, parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME),
        class = "htest")
}

