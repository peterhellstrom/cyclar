library(cyclar)

# function curve3d, ?curve3d

# Contours
# phi1 = (1+a1), phi2 = a2

# Period
period.surf <- function(phi1,phi2,method="all") {
	if (method == "complex") {
		if (phi1^2 + 4*phi2 > 0) NA
		else (2*pi)/(acos(phi1/(2*sqrt(-phi2))))
	}
	else if (method == "all") (2*pi)/(acos(phi1/(2*sqrt(-phi2))))

}

k <- c(2,3,3.5,4,4.5,5:10,20)
f.period <- curve3d(expr=period.surf(phi1=x,phi2=y, method="complex"), from=c(-2,-1), to=c(2,1), n=c(201,201), sys3d="contour", levels=k)
ar2.plot(k=k, kcol=1, klty=1, vlab=FALSE, vcontours=FALSE)
contour(f.period,add=TRUE,col=2,lty=3,levels=k, drawlabels=TRUE)

# Variance
variance.surf <- function(phi1,phi2,method="all") {
	if (method == "complex") {
		if (phi1^2 + 4*phi2 > 0) NA
		else (1-phi1)/((1-phi1-phi2)*(1-phi2+phi1)*(1+phi2))
	}
	else if (method == "all") (1-phi2)/((1-phi1-phi2)*(1-phi2+phi1)*(1+phi2))

}

v <- c(1.005,1.01,1.05,1.1,1.25,1.5,1.75,2,2.5,3,4,5,6,8,10)
f.variance <- curve3d(expr=variance.surf(phi1=x,phi2=y,method="all"), from=c(-2,-1), to=c(2,1), n=c(201,201), sys3d="contour", levels=v)
ar2.plot(v=v, klab=FALSE, kcontours=FALSE, vcontours=TRUE, vcol=2, vlty=1, cex.contours=0.7)
f.variance$z[which(f.variance$z > max(v))] <- NA
contour(f.variance,add=TRUE,col=4,lty=3,levels=v,drawlabels=TRUE)

# Royama's correlations

# Xt+1
royama.surfx <- function(a1,a2) sqrt((1+a2)*((a1+a2)*(a2-a1-2))) / sqrt(1-a2)

fx <- curve3d(expr=royama.surfx(a1=x,a2=y), from=c(-3,-1), to=c(1,1), n=c(201,201), sys3d="contour")
fx$x <- fx$x + 1
ar2.plot(klab=FALSE, vlab=FALSE, kcontours=FALSE, vcontours=FALSE)
contour(fx,add=TRUE,col=2,lty=2, levels=seq(0,1,0.1))

# Rt
royama.surfr <- function(a1,a2) sqrt((1+a2)*(2+a1-a2)/2)

fr <- curve3d(expr=royama.surfr(a1=x,a2=y), from=c(-3,-1), to=c(1,1), n=c(201,201), sys3d="contour", levels=seq(0,1,0.1))
fr$x <- fr$x + 1
ar2.plot(klab=FALSE, vlab=FALSE, kcontours=FALSE, vcontours=FALSE)
abline(1,-2, lwd=2)
contour(fr,add=TRUE,col=2,lty=2,levels=seq(0,1,0.1))
