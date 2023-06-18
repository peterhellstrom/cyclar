# THIS FUNCTION CAN NOT HANDLE DATA FRAME WITH NAs
# CHANGE THIS SOMETIME...
#' @export
ar2.pred.plot <- function(obj, sp=TRUE, col=c(1,3,4,5,2), lwd=c(2,1,1,1,1), lty=c(1,1,1,1,1)) {

	dat <- na.omit(obj$data)

	if (obj$method=="growth rate") {
		y <- dat[,"R"]
		x1 <- dat[,"Lag0"]
		x2 <- dat[,"Lag1"]
		}

	if (obj$method=="raw") {
		y <- dat[,"Lag0"]
		x1 <- dat[,"Lag1"]
		x2 <- dat[,"Lag2"]
		}

		mtime <- start(dat)[1]:end(dat)[1]

		aic <- c(NA,sapply(1:length(obj$models), function(i)AIC(obj$models[[i]])))
		mod.labs <- c("Observed","y ~ 1","y ~ y ~ x1","y ~ x2", "y ~ x1 + x2")
		mod <- obj$model

		if (mod=="lm") {
			r.squared <- c(NA,sapply(1:length(obj$models), function(i) summary(obj$models[[i]])$r.squared))
			leg.labs <- paste(mod.labs, round(r.squared,2), round(aic,2), sep=", ")
		}

		if (mod!="lm") {
			leg.labs <- paste(mod.labs, round(aic,2), sep=", ")
		}

		if (sp==FALSE) {
		plot(y,main=paste("Linear AR-models, model =",mod,sep=" "),
			bty="l",lty=lty[1], lwd=lwd[1], xlab="Time", ylab="Fitted value")

		points(mtime,predict(obj$models$fm0),col=col[2],lwd=lwd[2],type="l")
		points(mtime,predict(obj$models$fm1),col=col[3],lwd=lwd[3],type="l")
		points(mtime,predict(obj$models$fm2),col=col[4],lwd=lwd[4],type="l")
		points(mtime,predict(obj$models$fm12),col=col[5],lwd=lwd[5],type="l")

		legend("topleft",legend=, leg.labs,
			col=col, lty=lty, lwd=lwd, bty="n", cex=0.8)
		}

		if (sp==TRUE) {

		par(mfrow=c(2,2))

		for (i in 1:length(obj$models)) {

		fit.lab <- deparse(substitute(obj$models[[i]]))

		plot(y,main=leg.labs[i+1],
			bty="l",lty=lty[1], lwd=1, xlab="Time", ylab="Fitted value")

		points(mtime,predict(obj$models[[i]]),col=2,lwd=1,type="l")

		legend("topleft",legend=leg.labs[c(1,i+1)],
			col=c(1,2), lty=c(1,1), lwd=c(1,1), bty="n", cex=0.8)
		}
		}
}
