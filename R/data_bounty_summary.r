# Function that summarizes the bounty data from Norway
#' @export
bounty.summary <- function(x, species=NULL, plot=TRUE, scale=FALSE, plot.type="multiple",...) {

	if (!is.null(species)) x <- x[[species]]
	valid <- !is.na(x)
	#valid <- x > 0 & !(is.na(x))

	out <- rbind(
		sapply(1:ncol(x), function(i) length(which(valid[,1]==TRUE)) ),
		sapply(1:ncol(x), function(i) length(which(valid[,1]==FALSE)) ),

		sapply(1:ncol(valid), function(i) {
			t.temp <- time(x)[which(valid[,i]==TRUE)]
			range(t.temp)
		}),
		apply(x,2,sum,na.rm=T),
		apply(x,2,mean,na.rm=T),
		apply(x,2,sd,na.rm=T),
		apply(x,2,range,na.rm=T)
	)

	rownames(out) <- c("Valid","NA","Start","End","Sum","Mean","SD","min","max")

	if (plot) {

		if (scale) x <- scale(x)

		if (plot.type == "multiple") {
			plot(x[,1:9],nc=3,main=species,ylab=paste(species,"bounties"),...)
			plot(x[,10:ncol(x)],nc=3,main=species,ylab=paste(species,"bounties"),...)
		}

		if (plot.type == "single") {
			plot(x, plot.type="single",col=1:ncol(x),main=species,ylab=paste(species,"bounties"),...)
		}
	}

	out
}
