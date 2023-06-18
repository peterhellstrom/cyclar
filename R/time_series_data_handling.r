# ts.stand
# Standardize time series x to zero mean and unit variance:
#' @export
ts.stand <- function(x, log=FALSE, constant=0) {

	x <- na.omit(x)

	if (log==TRUE) x.temp <- log(x + constant)
	if (log==FALSE) x.temp <- x + constant

	x.stand <- (x.temp - mean(x.temp)) / sd(x.temp)
	# scale not used, because if class(x)=="ts", scale returns a numeric vector.
	#x.stand <- scale(x.temp)
	x.stand
}

# ts.binary
# Convert a numeric (time) series to a binary series
# See e.g. Steen et al 1990 (Oikos)
# Values larger than or equal to the specified value gets 1, below gets 0.
# fun = name of function used to calculate the value val
# Defaults to using fun="mean", i.e. mean of time series as dividing point.
#' @export
ts.binary <- function(x, val=NULL, fun="mean") {
  FUN <- match.fun(fun)
	if (is.null(val)) val <- FUN(x, na.rm=TRUE)
	ifelse(x >= val, 1, 0)
}

# lagmatrix
# Creates a matrix with lagged variables.
# x = a time series
# The number of lags is controlled with the lag.max argument.
# For example, if lag.max=2 lagmatrix returns a matrix with lag.max+1 = 3 columns.
# The first column contains the original series, the following columns are
# lagged variables padded with NA in top row(s). Number of NAs added at the top
# of each column equals the lag.
# The number of rows in the output matrix is equal to the length of series x.
# Note that the output is of class matrix, and not a time series object.
# However, a 'mts' object (multiple time series) is easily created by for instance
# ts(lagmatrix(x, lag.max=2))

# This function is called by ar2.lag, but can also be used independently.
#' @export
lagmatrix <- function(x,lag.max) {
  embed(x=c(rep(NA,lag.max),x), dimension=lag.max+1)
}
