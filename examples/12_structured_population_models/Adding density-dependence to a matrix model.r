########################################################################
# TUTORIAL: Adding density dependence to a matrix model
########################################################################
library(cyclar)

# Density dependence (Ricker) - overcompensatory
dd.r <- function(rate,beta,x) {
  rate*exp(-beta*x)
}

# Density dependence (Beverton-Holt) - compensatory
dd.bh <- function(rate,beta,x) {
  rate / (1 + beta*x)
}

# Plot to show the shape of these two functions:
curve(dd.r(rate=0.9, beta=0.05, x=x), from=0, to=100, font.lab=2, las=1, lwd=2, cex.lab=1.3, cex.axis=1.3,
xlab="Density", ylab="Vital rate", main="Density-dependence in vital rates")
curve(dd.bh(rate=0.9, beta=0.05, x=x), add=T, col=2, lwd=2)
legend("topright",c("Ricker","Beverton-Holt"),col=c(1,2), lty=c(1,1), lwd=c(2,2), bty="n", cex=1.3)

# Next, the plot recruitment curves.
# This example assumes that survival is density dependent,
# The resulting graph shows the number of surving eggs in relation to inital egg numbers.
# Notice the difference between the Ricker and the Beverton-Holt functions!

x <- 0:100 # set up vector with number in inital stage
# Calculate the number surviving to the next stage
surv.r <- x*dd.r(rate=0.9,beta=0.05,x=x)
surv.bh <- x*dd.bh(rate=0.9,beta=0.05,x=x)
# Plot the calculations
plot(x=c(0,100), y=c(0,16), type="n", xlab="Initial number in stage", ylab="Number surviving", main="Density-dependent recruitment curves",
font.lab=2, las=1, cex.lab=1.3, cex.axis=1.3)
points(x,surv.r,type="l",col=1,lwd=2)
points(x,surv.bh,type="l",col=2,lwd=2)
#abline(0,1,col=2,lty=2)
legend("topleft", c("Ricker","Beverton-Holt"),col=c(1,2), lty=c(1,1), lwd=c(2,2), bty="n", cex=1.3)

# Now, three examples.
# Assume a pre-breeding census (juvenily mortality must be incorporated in top row of matrix A!),
# and three censused stages. Stages 2 and 3 are reproductive.

# 1) Without density-dependence
# 2) With density-dependent juvenile survival
# 3) With density dependent fecundities, d-d effect more pronounced for stage 2.

# It is of course possible to introduce d-d in all stages, but data is most often
# not available to support this.
# First, without density-dependence.

# Define the vital rates:

# m is fecundity for each stage class (total number of offspring, we consider just females here, so this number is multiplied by 0.5 when calculating matrix entries.
m1 <- 0
m2 <- 4
m3 <- 6
# p is survival for each stage class
p0 <- 0.45
p1 <- 0.5
p2 <- 0.6
p3 <- 0.65

# Calculate the matrix

A <- matrix(c(
0.5*m1*p0, 0.5*m2*p0, 0.5*m3*p0,
p1, 0 ,0,
0, p2, p3),
nrow=3,ncol=3,byrow=T)

# Matrix projection
nsim <- 40
N <- matrix(ncol=3,nrow=nsim+1)
n0 <- c(10,10,10)
n <- n0

N[1,] <- n

for (t in 1:nsim) {
	n <- A %*% n
	N[t+1,] <- n
}

N
plot(rowSums(N),type="b",xlab="Time",ylab="Population density",font.lab=2,las=1,pch=16, main="No density depence - geometric growth")

# Now, add density dependence in juvenile survival
# Density-dependence is of Beverton-Holt type

# Static rates
m1 <- 0
m2 <- 4
m3 <- 6
p1 <- 0.5
p2 <- 0.6
p3 <- 0.65

# Dynamic rates (p0 = juvenile survival); Beverton-Holt d-d
# *.ini is rate when x=0, beta.* is the parameter that controls the strength of density-dependence
p0.ini <- 0.45
beta.p0 <- 0.0005

# The shape of the density dependence is:
curve(dd.bh(rate=p0.ini, beta=beta.p0, x=x), from=0, to=1000, font.lab=2, las=1, lwd=2, cex.lab=1.3, cex.axis=1.3,
xlab="Density", ylab="Vital rate", main="Density-dependence in juvenile survival")

# And the recruitment curve is:
x <- 0:200
surv.bh <- x*dd.bh(rate=p0.ini,beta=beta.p0,x=x)

plot(x=c(0,200), y=c(0,100), type="n", xlab="Initial number in stage", ylab="Number surviving", main="Recruitment curve",
font.lab=2, las=1, cex.lab=1.3, cex.axis=1.3)
points(x,surv.bh,type="l",col=1,lwd=2)
abline(0,1,col=2,lty=2)

# Set up the simulation
nsim <- 40 # number of simulations
N.p <- matrix(ncol=3,nrow=nsim+1) # Output matrix to hold the data
n0 <- c(10,10,10) # Initial stage vector
n <- n0

N.p[1,] <- n # Send initial vector to first row in output matrix

for (t in 1:nsim) {

	# Calculate density (I assume that the total population size is the relevant density here!)
	# Why not use use just x <- sum(n) as the density?
	# This would ignore the juvenile (non-censused) part (class 0-1) of the population, meaning that the age class 0-1 year olds
	# are only competing with individuals in the other (older) classes and not with individuals within their age class.
	# Consequently, density is underestimated, and the subsequent density-dependent survival rate is too high.
	x <- 0.5*(m2*n[2] + m3*n[3]) + sum(n) # The first term is the number of female offspring produced just after the census

	p0 <- dd.bh(rate=p0.ini, beta=beta.p0, x=x) # Calculate the value of the dynamic rate (juvenile survival), dependent on density x
	# Construct the matrix A, notice that the value for p0 is dynamic and changes for each iteration
	A <- matrix(c(
			0.5*m1*p0, 0.5*m2*p0, 0.5*m3*p0,
			p1, 0 ,0,
			0, p2, p3),
			nrow=3,ncol=3,byrow=T)
	# Multiply the transition matrix A with the stage vector
	n <- A %*% n
	N.p[t+1,] <- n # send result to the output matrix
}

N.p # View the results
# Plot the total population size against time
plot(rowSums(N.p),type="b",xlab="Time",ylab="Population density",font.lab=2,las=1,pch=16,
main="Density-dependence in juvenile survival")

# Now, add density dependence in fecundity
# m's are fecundities for each stage class, p's are survival rates

# Static rates
m1 <- 0
p0 <- 0.45
p1 <- 0.5
p2 <- 0.6
p3 <- 0.65

# Dynamic rates; Beverton-Holt d-d
# Set the density-dependent effect on second class fecundity to be slightly stronger:
m2.ini <- 4
m3.ini <- 6
beta.m2 <- 0.002
beta.m3 <- 0.001

curve(dd.bh(rate=m3.ini, beta=beta.m3, x=x), from=0, to=1000, font.lab=2, las=1, lwd=2, cex.lab=1.3, cex.axis=1.3,
xlab="Density", ylab="Vital rate", main="Density-dependent fecundity", ylim=c(0,6))
curve(dd.bh(rate=m2.ini, beta=beta.m2, x=x), add=T, col=2, lwd=2)
legend("bottomleft",c("m3 = Third stage class fecundity","m2 = Second stage class fecundity"),col=c(1,2), lty=c(1,1), lwd=c(2,2), bty="n", cex=1.3)

# Set up simulation
nsim <- 40
N.m <- matrix(ncol=3,nrow=nsim+1)
n0 <- c(10,10,10)
n <- n0

N.m[1,] <- n

for (t in 1:nsim) {

	x <- sum(n[2:3]) # Calculate density (I assume that the only the reproductive classes are competing!)
	m2 <- dd.bh(rate=m2.ini, beta=beta.m2, x=x) # Calculate the dynamic rates, m2 & m3
	m3 <- dd.bh(rate=m3.ini, beta=beta.m3, x=x)
	# Construct the matrix A, notice that the values for m2 and m3 are dynamic and change for each iteration
	A <- matrix(c(
			0.5*m1*p0, 0.5*m2*p0, 0.5*m3*p0,
			p1, 0 ,0,
			0, p2, p3),
			nrow=3,ncol=3,byrow=T)
	# Multiply the transition matrix A with the stage vector
	n <- A %*% n
	N.m[t+1,] <- n  #send result to the output matrix
}

N.m
plot(rowSums(N.m),type="b",xlab="Time",ylab="Population density",font.lab=2,las=1,pch=16,
main="Density-dependence in reproductive output")

# Compare the three different scenarios (by plotting the total population size):
plot(x=c(0,nsim), y=c(0,1000), type="n", xlab="Time", ylab="Population density", font.lab=2, las=1, main="Compare different scenarios")
points(0:nsim,rowSums(N),col=1,type="b",pch=15)
points(0:nsim,rowSums(N.p),col=2,type="b",pch=16)
points(0:nsim,rowSums(N.m),col=3,type="b",pch=17)
legend("bottomright",c("No density-dependence","D-d in juvenile survival","D-d in fecundity"), pch=15:17, col=1:3, lty=1, bty="n", cex=1.2)

# Also see the salmon example (Salmon_density_dependence.r), that illustrates the difference between
# the Ricker and Beverton-Holt functions!

# More details on density-dependent models are given in the books by Caswell (chapter 16) and Morris & Doak (chapter 8).
