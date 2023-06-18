########################################################################
# Salmon example
########################################################################

# This code was adopted from Morris & Doak (2002) "Quantitative Conservation Biology", p. 322
# The translation from their MATLAB-code to R-code was done by Peter Hellström, Stockholm University
# The original data and parameter estimates are from Ratner et al. (1997) Conservation Biology 11:879-889.

# The model is based on a pre-breeding census, with 5 age classes (ages 1-5).

# Density dependence (Ricker) - overcompensatory
dd.r <- function(rate,beta,x) { rate*exp(-beta*x) }
# Density dependence (Beverton-Holt) - compensatory
dd.bh <- function(rate,beta,x) { rate / (1 + beta*x) }

# Set up parameters
n0 <- c(0.7667, 0.6163, 0.4945, 0.3524, 0.1309) # Initial population vector

s0.max <- 0.002267 # Maximum egg survival
beta.s0 <- 0.001 # Density-dependent egg survival parameter

sx <- c(0.8, 0.8, 0.8, 0.8) # Survival probabilities
bx <- c(0, 0, 0.112, 0.532, 1) # Probability of breeding at each age
fx <- c(0, 0, 3185, 3940, 4336) # Eggs per breeding female
fx.mult <- c(1,2,4,8,16) # Sequential doubling of fecundities

# Plot
plot(x=c(0,100),y=c(0,16),type="n",
	xlab="Initial number of eggs", ylab="Number surviving to age 1", font.lab=2, las=1,
	main="Density-dependent survival of salmon eggs")
curve(x*dd.bh(rate=0.9,beta=0.05,x=x), from=0, to=100, add=T, col=1, lwd=2)
curve(x*dd.r(rate=0.9,beta=0.05,x=x), from=0, to=100, add=T, col=2, lwd=2, lty=2)
legend("topleft",c("Beverton-Holt","Ricker"),col=1:2, lwd=2, lty=1:2,bty="n")

########################################################################
# Set up graph window
dev.new(width=12,height=6)
op <- par(list(
mfrow=c(2,3),
mar=c(5,5,2,2)))

for (i in 1:length(fx.mult)) { # Start of outer for-loop

# For each iteration i, multiply the original parameter estimates fx with a sequential doubling factor
fx <- c(0, 0, 3185, 3940, 4336) # Eggs per breeding female
fx <- fx*fx.mult[i]

# Fill the transition matrix A with survival values
# Since all spawners died after spawning, survival rates
# must be adjusted for probability of breeding.
A <- matrix(0,ncol=5,nrow=5)
A[2,1] <- sx[1]
A[3,2] <- sx[2]
A[4,3] <- sx[3]*(1-bx[3])
A[5,4] <- sx[4]*(1-bx[4])

A.r <- A # Matrix for Ricker scenario
A.bh <- A # Matrix for Beverton-Holt scenario

# Set up the simulation
tmax <- 100 # Time horizon

# Create output matrix for Ricker scenario
N.r <- matrix(ncol=5,nrow=tmax+1)
N.r[1,] <- n0

# Create output matrix for Beverton-Holt scenario
N.bh <- matrix(ncol=5,nrow=tmax+1)
N.bh[1,] <- n0

# Set inital vectors
n.r <- n0 # Start at n0
n.bh <- n0 # Start at n0

	for (t in 1:tmax) { # Start of inner for-loop
		# Compute total eggs produced (sum all eggs produced for females aged 3-5)
		eggs.r <- bx*fx*n.r
		eggs.bh <- bx*fx*n.bh
		# Update survival rates, given the number of eggs
		s0.r <- dd.r(rate=s0.max, beta=beta.s0, x=sum(eggs.r))
		s0.bh <- dd.bh(rate=s0.max, beta=beta.s0, x=sum(eggs.bh))
		# Calculate fecundities
		f.r <- bx*fx*s0.r
		f.bh <- bx*fx*s0.bh
		# Update matrix with fecundities
		A.r[1,] <- f.r
		A.bh[1,] <- f.bh
		# Multiply the A's with stage vectors:
		n.r <- A.r %*% n.r
		n.bh <- A.bh %*% n.bh
		# Send data to output matrices:
		N.r[t+1,] <- n.r
		N.bh[t+1,] <- n.bh
	} 	# End of inner for-loop

# Calculate number of spawners in each year
spawners.r <- rowSums(t(sapply(1:nrow(N.r), function(i) bx*N.r[i,])))
spawners.bh <- rowSums(t(sapply(1:nrow(N.bh), function(i) bx*N.bh[i,])))

# Plot number of spawners against time
plot(x=c(0,100),y=c(0,1.1), type="n", xlab="Time (years)", ylab="Spawners", font.lab=2, las=1, cex.axis=1.2, cex.lab=1.2,
main=paste("Fecundities multiplied by", fx.mult[i]),bty="l")
points(0:100, spawners.r, type="l", col=1, lwd=2)
points(0:100, spawners.bh, type="l", col=2, lwd=2)
if (i==1) legend("bottomleft",c("Ricker","Beverton-Holt"),col=c(1,2),lwd=c(2,2),bty="n",cex=1.2)

} # End of outer for-loop

# Plot juvenile survival rate against initial egg density:
curve(dd.r(rate=s0.max, beta=beta.s0, x=x), from=0, to=25000, font.lab=2, lwd=2, cex.lab=1.3, cex.axis=1.3, bty="l",
xlab="Initial egg density", ylab="Survival from egg to 1-year old", main="Density-dependent juvenile survival")
curve(dd.bh(rate=s0.max, beta=beta.s0, x=x), add=T, col=2, lwd=2)
abline(h=s0.max,lty=2)
legend("topright",c("Ricker","Beverton-Holt"),col=c(1,2), lty=c(1,1), lwd=c(2,2), bty="n", bg="white", cex=1.3)

# Restore graphical parameters
par(op)
