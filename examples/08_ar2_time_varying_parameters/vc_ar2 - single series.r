library(cyclar)

a.tv <- function(x,a,b,c) a / (1 + exp(-(x - b) / -c))

n <- 50
a1.t <- a.tv(x=0:n, a=1.3, b=-0.0025, c=10) # Direct d-d
a2.t <- a.tv(x=0:n, a=-0.8, b=30, c=8) # Delayed d-d
a0 <- 2.1
a0.t <- rep(a0,(n+1))


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

route <- c(1,4,6)

nsteps <- length(route)-1
yrs <- 25
vals <- matrix(c(
	rep(yrs,nsteps),
	rep(0,nsteps),
	rep(0.2,nsteps)),
nrow=nsteps, ncol=3, byrow=FALSE)

ps <- ar2.sim.tv.gen(npoints=8, coord.mat=coord.mat, route=route, vals=vals, gam=TRUE)
a1.t <- ps$a1
a2.t <- ps$a2
##
x <- ar2.sim.tv(n=ps$ntime, a0=rep(0,ps$ntime), a1=ps$a1, a2=ps$a2, sd=0.15, plot=FALSE)

plot(x$ts.x, font.lab=2, las=1, xlab="Time", ylab="x")
abline(v=cumsum(vals[,1]), col=2, lty=2)

#ar2.plot(ylims=c(-1,0),triangle=FALSE)
#ar2.arrows(a1.t, a2.t)

z <- vc.inp(x=x, method="x")
# z <- vc.inp(x=x, method="r")

# RUN ESTIMATION:
# "C:/Program Files/Varying Coefficients Estimation/VC.exe"

# Check output
z.out <- vc.out(vc.dat=z, method="x")
# vc.out(x=vc.dat, method="r")

plot(z.out[,"a_x1"], z.out[,"a_x2"], type="n")
arrows(z.out[,"a_x1"][-nrow(z.out)], z.out[,"a_x2"][-nrow(z.out)], z.out[,"a_x1"][-1], z.out[,"a_x2"][-1],code=2, length=0.1)
points(z.out[,"a_x1"], z.out[,"a_x2"], col=2, pch=16)

