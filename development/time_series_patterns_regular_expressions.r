# Finding patterns in time series using regular expressions
# http://www.r-bloggers.com/finding-patterns-in-time-series-using-regular-expressions/

x <- lynx
delta <- (sign(diff(x)) == 1) + 0
delta_str <- paste(delta, collapse = "")

matches <- gregexpr("00+", delta_str, perl = TRUE)[[1]]
recessions <- sapply(1:length(attr(matches, "match.length")), function(ind) matches[ind] + 0:(m.length[ind]))

hl <- function(inds) lines(time(x)[inds], x[inds], col = "red", lwd = 3)

plot(x, main = "x", ylab = "furs", xlab = "Year")
tmp <- sapply(recessions, hl)  # Used for side-effects only

