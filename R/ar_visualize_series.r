# http://www.r-bloggers.com/visualizing-autoregressive-time-series/
#' @export
ar1.graph <- function(phi,n=500) {
	nf <- layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow=TRUE), respect=TRUE)
	e=rnorm(n)
	X=rep(0,n)
	for(t in 2:n) X[t]=phi*X[t-1]+e[t]
	plot(X[1:n],type="l",ylab="")
	abline(h=mean(X),lwd=2,col="red")
	abline(h=mean(X)+2*sd(X),lty=2,col="red")
	abline(h=mean(X)-2*sd(X),lty=2,col="red")
	u=seq(-1,1,by=.001)
	plot(0:1,0:1,col="white",xlab="",ylab="",axes=FALSE,ylim=c(-2,2),xlim=c(-2.5,2.5))
	polygon(c(u,rev(u)),c(sqrt(1-u^2),rev(-sqrt(1-u^2))),col="light yellow")
	abline(v=0,col="grey")
	abline(h=0,col="grey")
	points(1/phi,0,pch=19,col="red",cex=1.3)
	plot(0:1,0:1,col="white",xlab="",ylab="",axes=FALSE,ylim=c(-.2,.2),xlim=c(-1,1))
	axis(1)
	points(phi,0,pch=19,col="red",cex=1.3)
	acf(X,lwd=3,col="blue",main="",ylim=c(-1,1))
	pacf(X,lwd=3,col="blue",main="",ylim=c(-1,1),xlim=c(0,16))
}

#' @export
ar2.graph <- function(phi1, phi2, n=500) {
	nf <- layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow=TRUE), respect=TRUE)

	e=rnorm(n)
	X=rep(0,n)
	for(t in 3:n) X[t]=phi1*X[t-1]+phi2*X[t-2]+e[t]

	plot(X[1:n],type="l",ylab="",main="AR(2)")
	abline(h=mean(X),lwd=2,col="red")
	abline(h=mean(X)+2*sd(X),lty=2,col="red")
	abline(h=mean(X)-2*sd(X),lty=2,col="red")

	P=polyroot(c(1,-phi1,-phi2))

	u=seq(-1,1,by=.001)
	plot(0:1,0:1,col="white",xlab="",ylab="",axes=FALSE,ylim=c(-2,2),xlim=c(-2.5,2.5))
	polygon(c(u,rev(u)),c(sqrt(1-u^2),rev(-sqrt(1-u^2))),col="light yellow")
	abline(v=0,col="grey")
	abline(h=0,col="grey")
	points(P,pch=19,col="red",cex=1.3)

	plot(0:1,0:1,col="white",xlab="",ylab="",axes=FALSE,xlim=c(-2.1,2.1),ylim=c(-1.2,1.2))
	polygon(c(-2,0,2,-2),c(-1,1,-1,-1),col="light green")
	u=seq(-2,2,by=.001)
	lines(u,-u^2/4)
	#abline(v=seq(-2,2,by=.2),col="grey",lty=2)
	#abline(h=seq(-1,1,by=.2),col="grey",lty=2)
	segments(0,-1,0,1)
	axis(1)
	axis(2)
	points(phi1,phi2,pch=19,col="red",cex=1.3)

	acf(X,lwd=3,col="blue",main="ACF",ylim=c(-1,1))
	pacf(X,lwd=3,col="blue",main="PACF",ylim=c(-1,1),xlim=c(0,16))
}

#ar2.graph(1,-0.7,n=100)

#ar1.graph(0.8)
