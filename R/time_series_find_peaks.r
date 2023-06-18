# Functions to find peaks in a (time) series
# I used the function turnpoints from {pastecs} in e.g. ts.diag.acf
# I haven't fully tested these function, but if the can perhaps be used in the
# Markov-processes functions in cyclar.

# see ?max.col

# http://finzi.psych.upenn.edu/R/Rhelp02a/archive/33097.html
#' @export
peaks.old <- function(x,span=3) {
	# span has to be odd number
	if (span %% 2 == 0) stop("span has to be odd number")
	z <- embed(x, span)
	s <- span %/% 2
	v <- max.col(z) == 1 + s
	out <- c(rep(FALSE,s), v)
	out <- out[1:(length(out)-s)]
	out
}

# http://tolstoy.newcastle.edu.au/R/e2/help/07/02/10145.html
#' @export
peaks <- function(x, span=3, ties.method = "first") {
	if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
	z <- embed(x, span)
	s <- span %/% 2
	v <- max.col(z, ties.method=ties.method) == 1 + s
	pad <- rep(FALSE, s)
	out <- c(pad, v, pad)
	out
}

#' @export
peaks.plot <- function(x,span=3,pch=16,col=2,...) {

	inds <- as.numeric(peaks(x,span=span))

	plot(x,...)
	points(x,
	pch = ifelse(inds == 1, pch, 0),
	col = ifelse(inds == 1, col, 0))
}
