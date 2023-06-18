################################################################################
library(cyclar)
################################################################################

# Find parameters for a given period (k) & intrinsic process variance (v).
# I have implemented this in cyclar with the common root finding technique.

# A test
k <- 3
v <- 2
p <- ar2.parms(k=k,v=v)
ar2.period(p$a1,p$a2)

# The solution to this problem was mostly outlined in the following file:
# AR(2)_variance_contours_inverse_functions.r (in Examples > 03 Plot).

# We know the equation for the k-contours:
# a2 <- -(1+a1)^2/(4*cos(2*pi/k)^2)

# We also know the intrinsic process variance v as a function of (1+a1) & a2:
# v <- (1-a2) / ((1-(1+a1)-a2)*(1-a2+(1+a1))*(1+a2))
# I couldn't solve the variance expression for a2, but for (1+a1):
# Setting x = (1+a1), y = a2
# x2 <- y^2 - 2*y + 1 - ((1 - y)/((1 + y)*v))
# The full solution is the sqrt of the above equation.
# Switch x's and y's (this is the function in ar2.v.inv):
# y2 <- x^2 - 2*x + 1 - ((1 - x)/((1 + x)*v))

# Do the same for the k-conturs.
# Set y = a2, x = (1+a1) and solve for x^2:
# x^2 <- -4*y*cos(2*pi/k)^2
# Switch x's and y's
# y^2 <- -4*x*cos(2*pi/k)^2
# y1 <- 2*cos(2*pi/k)*sqrt(-x)
# y2 <- -2*cos(2*pi/k)*sqrt(-x)

# We know have two expressions for y^2 that can be combined, see below!

ar2.k.inv(p$a2,k=k)
ifelse(k < 4, -ar2.v.inv(p$a2,v=v), ar2.v.inv(p$a2,v=v))

# Plot 'inverse. functions
plot(x=c(-1,0),y=c(-2,2), type="n", xlab=expression(beta[2]), ylab=expression(1+beta[1]), main="Find intersection between k & v")
curve(ar2.k.inv(x=x,k=k),n=1001,from=-1,to=0,col=2,add=TRUE)
curve(ar2.v.inv(x=x,v=v),n=1001,from=-sqrt(1-(1/v)),to=0,col=4,lty=2,add=TRUE)
curve(-ar2.v.inv(x=x,v=v),n=1001,from=-sqrt(1-(1/v)),to=0,col=4,lty=2,add=TRUE)
abline(h=0, col=1)
lines(x=c(p$a2,p$a2,-1),
	y=c(-2,ifelse(k < 4, -ar2.v.inv(p$a2,v=v), ar2.v.inv(p$a2,v=v)),ifelse(k < 4, -ar2.v.inv(p$a2,v=v), ar2.v.inv(p$a2,v=v))),
	col=3,lty=2)
points(p$a2,1 + p$a1, pch=16,col=2)
legend("bottomright",legend=c("k","v"),col=c(2,4),lty=1:2,bty="n",title="Contours")

# Use root-finding to get the parameters, cont.

# For v.inv
# y^2 = x^2 - 2*x + 1 - ((1 - x )/((1 + x)*v))

# For k.inv
# y = 2*cos(2*pi/k)*sqrt(-x)
# y^2 = -4*cos(2*pi/k)^2*x

# So we have
# x^2 - 2*x + 1 - ((1 - x )/((1 + x)*v)) = -4*cos(2*pi/k)^2*x
# x^2 - 2*x + 1 - ((1 - x )/((1 + x)*v)) + 4*cos(2*pi/k)^2*x = 0

# Solving for x is possible, but doing this is rather difficult...
# I tried to get a solution with Mathematica, but the result was long and hard to interpret.
# Another option is to find the root of the function
# f(x) = x^2 - 2*x + 1 - ((1 - x )/((1 + x)*v)) + 4*cos(2*pi/k)^2*x
# This can easily be solved in R with uniroot.
# Note that the root only needs to be sought in the interval c(-sqrt(1-(1/v)),0)
# Too see this, solve the equation [v.inv] x^2 - 2*x + 1 - ((1 - x )/((1 + x)*v)) = 0 for x.

# I have implemented this in function ar2.parms.root().
# The main function to call is ar2.parms(), which calls ar2.parms.root().
# ar2.parms can handle vectors with k & v as inputs and has a plot options (calls ar2.plot).

# Plot function f.root
f.root <- function(x,k,v) x^2 - 2*x + 1 - ((1-x)/((1+x)*v)) + 4*cos(2*pi/k)^2*x

# Change k
k <- 2:10
v <- 4
p <- ar2.parms(k=k,v=rep(v,length(k)),cex.contours=0.8,triangle=FALSE,ylim=c(-1,0),main="Change quasi-period")

plot(x=c(-1,0), y=c(-2,2), type="n", xlab=expression(beta[2]), ylab=expression(f(beta[2],k,v)), main="Find AR(2)-parameters")
abline(h=0, col=2, lty=2)
for (i in 1:length(k)) curve(f.root(x=x,k=k[i],v=v),col=k[i],from=-1,to=0,add=TRUE)
points(p$a2,rep(0,length(k)),col=k,pch=16)
legend("bottomright",legend=k,col=k,lty=1,bty="n",cex=0.8,title="k")

# Change v
v <- 2:8
k <- 3
p <- ar2.parms(k=rep(k,length(v)),v=v,cex.contours=0.8,triangle=FALSE,ylim=c(-1,0),main="Change variance")

plot(x=c(-1,0), y=c(-2,2), type="n", xlab=expression(beta[2]), ylab=expression(f(beta[2],k,v)), main="Find AR(2)-parameters")
abline(h=0, col=2, lty=2)
for (i in 1:length(v)) curve(f.root(x=x,k=k,v=v[i]),col=v[i],from=-1,to=0,add=TRUE)
points(p$a2,rep(0,length(v)),col=v,pch=16)
legend("bottomright",legend=v,col=v,lty=1,bty="n",cex=0.8,title="v")

################################################################################

# Check that input parameters can be retrieved
k <- 3:10
v <- rep(7,length(k))
p <- ar2.parms(k=k,v=v,formula=TRUE,par.name="Omega",cex.contours=0.7)
ar2.period(p$a1,p$a2) # check, should equal vectors k & v

################################################################################

# Compare with previous function ar2.parms.old
ar2.parms.old(k=5,v=3)
ar2.parms.old(k=4,v=8)

ar2.parms(k=5,v=3,list=FALSE)
ar2.parms(k=4,v=8,list=FALSE)

k <- 2:20
v <- rep(6,length(k))

p1 <- ar2.parms(k=k,v=v,plot=TRUE,list=FALSE)
p2 <- ar2.parms.old(k=k,v=v,plot=TRUE,output=cbind)
p1; p2
all.equal(p1,p2,tol=0)


v <- seq(1.5,10,0.25)
k <- runif(n=length(v),min=2,max=10)

system.time(p3 <- ar2.parms(k=k,v=v,plot=TRUE,list=FALSE,cex.contours=0.7))
system.time(p4 <- ar2.parms.old(k=k,v=v,plot=TRUE,output=cbind,cex.contours=0.7))
p3; p4
all.equal(p3,p4,tol=0)

################################################################################
