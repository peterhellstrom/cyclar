#' @export
ar2.vc.gam.fit <- function(x, n=100, method="xx", model="gam", se.mult=2) {

	nsmooth <- length(x$smooth)
	out.dat <- vector("list",nsmooth)

	for (i in 1:nsmooth) {

		if (model == "gam") {

			Raw <- x$model[x$smooth[[i]]$term]
			if (method == "xx") n <- n
			if (method == "raw") n <- nrow(Raw)
		}

		if (model=="gamm") {

			tab <- table(x$model[x$smooth[[i]]$term])
			n.tab <- length(tab)
			un.tab <- unique(tab)
			Raw <- as.numeric(names(tab))
			if (method == "xx") n <- n
			if (method == "raw") n <- length(as.numeric(names(tab)))
		}

		xx <- seq(min(Raw), max(Raw), length = n)

		by <- rep(1, n)
		if (method == "xx") dat <- data.frame(x = xx, by = by)
		if (method == "raw") dat <- data.frame(x = Raw, by = by)
		names(dat) <- c(x$smooth[[i]]$term, x$smooth[[i]]$by)

		X <- PredictMat(x$smooth[[i]], dat)
			first <- x$smooth[[i]]$first.para
			last <- x$smooth[[i]]$last.para
		p <- x$coefficients[first:last]
		fit <- X %*% p
		se.fit <- sqrt(rowSums((X %*% x$Vp[first:last, first:last]) * X))
		ll <- fit - se.mult*se.fit
		ul <- fit + se.mult*se.fit

		out <- data.frame(
				x = xx,
				fit = fit,
				se.fit = se.fit,
				ll = ll,
				ul = ul)

		out.dat[[i]] <- out
	}
	out.dat
}

#' @export
ar2.vc.gam.plot <- function(x, fit="x", inds=c(1,2)) {

	x <- lapply(inds, function(i) x[[i]])

	if (fit == "R") {
			x[[inds[1]]][,"fit"] <- 1 + x[[inds[1]]][,"fit"]
			x[[inds[1]]][,"ll"] <- 1 + x[[inds[1]]][,"ll"]
			x[[inds[1]]][,"ul"] <- 1 + x[[inds[1]]][,"ul"]
	}

	ylims <- range(sapply(1:2, function(i) t(x[[i]][,c("ll","ul")])))

	main.labs1 <- c(
		"Smoother - direct d-d",
		"Smoother - delayed d-d")
	main.labs2 <- c(
		"Parameter estimates and CI's",
		"Parameter estimates and CI's")

	op <- par(list(
			mar=c(5,4,3,2),
			mfrow=c(2,2)))


		for (i in 1:2) {
			plot(x[[i]]$x, x[[i]]$fit, type="l",lty=1,col=1,ylim=ylims,xlab="Time",
				ylab=ifelse(i == 1, expression(1 + beta[1]), expression(beta[2])),
				font.lab=2,las=1,bty="l",main=main.labs1[i])
			points(x[[i]]$x, x[[i]]$ll, type="l", lty=2)
			points(x[[i]]$x, x[[i]]$ul, type="l",lty=2)
		}

		for (i in 1:2) {
			plot(x[[i]]$x, x[[i]]$fit, type="n", ylim=ylims,xlab="Time",
				ylab=ifelse(i == 1, expression(1 + beta[1]), expression(beta[2])),
				font.lab=2,las=1,bty="l",main=main.labs2[i])
			abline(h=0, lty=2, col=2)
			points(x[[i]]$x, x[[i]]$fit,pch=16,cex=1)
			arrows(x0=x[[i]]$x, y0=x[[i]]$ll, x1=x[[i]]$x, y1=x[[i]]$ul, length=0, col=1, lwd=1)
		}
	par(op)
}
