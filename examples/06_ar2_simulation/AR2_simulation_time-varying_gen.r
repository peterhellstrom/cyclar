################################################################################
library(cyclar)
################################################################################

npoints <- 2
coord.mat <- matrix(
	c(
		1.2,-0.9,
		0.1,-1
	),
ncol=2, nrow=npoints, byrow=TRUE, dimnames=list(1:npoints,c("x","y")))
route <- c(1,1,2,2)
nsteps <- length(route)-1
vals <- matrix(c(
	c(60,0,40),
	c(1.5,1.5,4),
	c(0.15,0.15,0.05)),
nrow=nsteps, ncol=3, byrow=FALSE)

ps <- ar2.sim.tv.gen(npoints=npoints, coord.mat=coord.mat, route=route, vals=vals, gam=FALSE)
test.tv <- ar2.sim.tv(n=ps$ntime, a0=rep(c(1.5,4),times=c(60,40)), a1=ps$a1, a2=ps$a2, sd=sqrt(0.15), plot=TRUE)

####
# Repeat with different sd per "strata" and convert mu to a0
w.sd <- rep(c(sqrt(0.15),sqrt(0.05)),times=c(60,40))
w <- rnorm(100,0,w.sd)

a01 <- ar2.a0(mu=1.5,a1=1.2,a2=-.9)
a02 <- ar2.a0(mu=4,a1=0.1,a2=-1)

ar2.mu(a01,a1=1.2,a2=-.9)
ar2.mu(a02,a1=0.1,a2=-1)

test.tvw <- ar2.sim.tv(n=ps$ntime, a0=rep(c(a01,a02),times=c(60,40)), a1=ps$a1, a2=ps$a2, sd=sqrt(0.15), innov=w, plot=TRUE)
############

w <- rnorm(60,0,sqrt(0.15))
test.tvw2 <- ar2.sim.tv(n=60, a0=rep(1.5,60), a1=rep(1.2,60), a2=rep(-0.9,60), sd=sqrt(0.2), innov=w, plot=TRUE, plot.log=TRUE)

n1 <- 60
n2 <- 40

test2 <- ar2.sim(n=n1, a0=1.5, a1=1.2, a2=-0.9, sd=sqrt(0.15), plot=TRUE, plot.log=TRUE)
test3 <- ar2.sim(n=n2, a0=4, a1=0.1, a2=-.99, sd=sqrt(0.05), plot=TRUE, plot.log=TRUE)

test.tv2 <- ar2.sim.tv(n=n1, a0=rep(1.5,n1), a1=rep(1.2,n1), a2=rep(-0.9,n1), sd=sqrt(0.15), plot=TRUE, plot.log=TRUE)
test.tv3 <- ar2.sim.tv(n=n2, a0=rep(4,n2), a1=rep(0.1,n2), a2=rep(-.99,n2), sd=sqrt(0.05), plot=TRUE, plot.log=TRUE)

ts.stat(test2$ts.x)
ts.stat(test3$ts.x)
ts.stat(test.tv2$ts.x)
ts.stat(test.tv3$ts.x)
############

w <- c(rnorm(n1,0,0.2),rnorm(n2,0,0.2))
p <- ar2.parms(k=5,v=4)
a1 <- p[1]
a2 <- p[2]
a01 <- 2.2
a02 <- 0
test.tv4 <- ar2.sim.tv(n=n1+n2, a0=c(rep(a01,n1),rep(a02,n2)), a1=c(rep(a1,n1),rep(a1,n2)), a2=c(rep(a2,n1),rep(a2,n2)),
	sd=sqrt(0.2), innov=w, plot=TRUE, plot.log=TRUE, windows=TRUE)

plot(test.tv4$ts.x)

library(sowas)

nvoice <- 125
noctave <- 2
nreal <- 100
nvoice*noctave+1 # scales

x <- test.tv4$ts.x
test <- wsp(x, s0=2, nvoice=nvoice, noctave=noctave, nreal=nreal, units="years", arealsiglevel=0.9)

############
npoints <- 8
coord.mat <- matrix(
c(
0,-0.2,
0.75, -0.2,
0.75, -0.8,
0, -0.8,
-1, -0.8,
-1, -0.2,
0.4, 0.3,
0,-0.5
),
ncol=2, nrow=npoints, byrow=TRUE, dimnames=list(1:npoints,c("x","y")))


route <- c(1,4,6,5,1)

route <- c(2,1,6)
route <- c(5,8,2)

route <- c(4,1,3,2,4)

route <- c(4,1,7)

route <- c(4,4,6,5,5)

route <- c(2,8,5) # Good performance...
route <- c(2,5)

route <- c(6,8,3) # Rather poor performance, sometimes spot-on...

route <- c(4,8,1) # Rather poor performance

nsteps <- length(route)-1
yrs <- 25
vals <- matrix(c(
	rep(yrs,nsteps),
	rep(0,nsteps),
	rep(0.2,nsteps)),
nrow=nsteps, ncol=3, byrow=FALSE)

ps <- ar2.sim.tv.gen(npoints=8, coord.mat=coord.mat, route=route, vals=vals, gam=TRUE)
##
test <- ar2.sim.tv(n=ps$ntime, a0=rep(0,ps$ntime), a1=ps$a1, a2=ps$a2, sd=0.15, plot=TRUE)

ts.stat(test$ts.x, plot=TRUE)
spec.plot(test$ts.x)
################################################################################
