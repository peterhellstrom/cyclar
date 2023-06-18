######################################################################
# Peregrine falcon model
# Wooton & Bell (1992) Ecol Appl 2:307-321

######################################################################
# Basic model (bm): No spatial structure
# Pre-breeding census with two classes

# Variables
# N = number of paired, territorial females
# F = number of age 1 (nonterritorial females)

# Parameters (with standard deviations, *.sd)
b <- 0.69 # annual number of female fledglings per paired, territorial females
b.sd <- 0.11

y <- 0.36 # Probability of fledglings surviving to age 1
y.sd <- 0.05

r <- 0.72 # Probability of age 1 birds surviving and recruiting to territorial population
r.sd <- 0.05

s <- 0.77 # probability of survival by a territorial adult to the following year
s.sd <- 0.08

i <- c(0,50) # 0 or 50, number of captive-reared fledglings introduced


# Transition matrix
A <- matrix(c(
	0, y*b,
	r, s),
ncol=2,nrow=2,byrow=T)

nyears <- 100

# N.bm - model without hacking

N.bm <- array(dim=c(nyears+1,ncol(A),length(i)))

n0 <- c(55,90) # F = 55, N = 90
N.bm[1,,] <- n0
n <- n0

for (j in 1:length(i)) {
	n <- n0
	for (t in 1:nyears) {
		n <- (A %*% n) + matrix(c(y*0.5*i[j],0),ncol=1,nrow=nrow(A))
		N.bm[t+1,,j] <- n
		}
}


######################################################################
# Metapopulation model

# Northern population = n
# Southern population = s

# Variables
# Nn = number of paired, territorial females, north
# Fn = number of age 1 (nonterritorial females), north
# Ns = number of paired, territorial females, south
# Fs = number of age 1 (nonterritorial females), south

# Parameters (with standard deviations, *.sd)
bn <- 0.71 # annual number of female fledglings per paired, territorial females
bn.sd <- 0.11

bs <- 0.53
bs.sd <- 0.24

yn <- 0.36 # Probability of fledglings surviving to age 1
yn.sd <- 0.05

ys <- 0.36
ys.sd <- 0.05

rn <- 0.72 # Probability of age 1 birds surviving and recruiting to territorial population
rn.sd <- 0.05

rs <- 0.72
rs.sd <- 0.05

sn <- 0.77 # probability of survival by a territorial adult to the following year
sn.sd <- 0.08

ss <- 0.77
ss.sd <- 0.08


m <- 0.27 # probability of a nonterritorial female migrating into the other population
m.sd <- 0.05

i <- c(0,50) # 0 or 50, number of captive-reared fledglings introduced
f <- c(0, 0.5, 1) # fraction of captive-reared fledglings introduced into the northern subpopulation

h <- (1-m) # Probability of remaining in population (h = (1-m) )
mc <- 1 # Cost of moving between subpopulations (negligible)

# Transition matrix
As <- matrix(c(
	0, yn*bn, 0, 0,
	rn*h, sn, rn*m*mc, 0,
	0, 0, 0, ys*bs,
	rs*m*mc, 0, rs*h, ss),
ncol=4,nrow=4,byrow=T)

nyears <- 100

# N.mp - metapopulation model

N.mp <- array(dim=c(nyears+1, ncol(As), length(i)))

n0 <- c(19,60,20,29) # Fn = 19, Nn = 60, Fs = 20, Ns = 29)
N.mp[1,,] <- n0
n <- n0

# First simulation, only include a fraction of 0.5 captive-reared fledglings introduced in each subpopulation
for (j in 1:length(i)) {
	n <- n0
	for (t in 1:nyears) {
		n <- (As %*% n) + matrix(c(f[2]*yn*0.5*i[j], 0, (1-f[2])*ys*0.5*i[j], 0),ncol=1,nrow=nrow(As))
		N.mp[t+1,,j] <- n
		}
}

# Plot number of breeding pairs
# Figure 3 in original paper (p. 313).
dev.new()
op <- par(mar=c(5,5,4,2))
plot(x=c(0,nyears),y=c(0,175), type="n", xlab="Time (years)", ylab="", font.lab=2, las=1, cex.axis=1.3, cex.lab=1.3, main="California peregrine population model",bty="l",yaxt="n")

axis(2,at=(seq(0,175,25)),font=1,las=1,cex.axis=1.3)
mtext(side=2, line=4, text="No. of breeding pairs", font=2, cex=1.3)

points(0:nyears, N.bm[,,1][,2], pch=15, type="b")
points(0:nyears, N.bm[,,2][,2], type="b")

points(0:nyears, rowSums(N.mp[,,1][,c(2,4)]),pch=15, col=2, type="b")
points(0:nyears, rowSums(N.mp[,,2][,c(2,4)]), col=2, type="b")

legend("topleft",c("No spatial structure", "No spatial structure + 50 fledglings", "Spatial structure", "Spatial structure + 25 fledglings in each subpop"), col=c(1,1,2,2), pch=c(15,1,15,1), bty="n", bg="white")

# Use the same model to evaluate differences in the strategy of allocating 50 introduced fledglings to either or both of the two subpopulations
# See Figure 4, p. 314

# Transition matrix
As <- matrix(c(
	0, yn*bn, 0, 0,
	rn*h, sn, rn*m*mc, 0,
	0, 0, 0, ys*bs,
	rs*m*mc, 0, rs*h, ss),
ncol=4,nrow=4,byrow=T)

nyears <- 100

# N.mp - metapopulation model

N.mp <- array(dim=c(nyears+1, ncol(As), length(f)))

n0 <- c(19,60,20,29) # Fn = 19, Nn = 60, Fs = 20, Ns = 29)
N.mp[1,,] <- n0
n <- n0

# Second simulation, vary the fraction of 50 fledglings introduced in each subpopulation
# Let number of introduced fledglings be 50 in all cases
for (j in 1:length(f)) {
	n <- n0
	for (t in 1:nyears) {
		n <- (As %*% n) + matrix(c(f[j]*yn*0.5*i[2], 0, (1-f[j])*ys*0.5*i[2], 0),ncol=1,nrow=nrow(As))
		N.mp[t+1,,j] <- n
	}
}


dev.new()
op <- par(mar=c(5,5,4,2))
plot(x=c(0,nyears),y=c(85,105), type="n", xlab="Time (years)", ylab="No. of breeding pairs", font.lab=2, las=1, cex.axis=1.3, cex.lab=1.3, main="California peregrine population model",bty="l")

points(0:nyears, rowSums(N.mp[,,1][,c(2,4)]),pch=15, col=1, type="b")
points(0:nyears, rowSums(N.mp[,,2][,c(2,4)]), col=2, type="b")
points(0:nyears, rowSums(N.mp[,,3][,c(2,4)]), pch=16, col=4, type="b")

legend("topleft",c("0 N : 50 S", "25 N : 25 S", "50 N : 0 S"), col=c(1,2,4), pch=c(15,1,16), bty="n", bg="white")

# Metapopulation with limited number of territories
# Source-sink model (N.sos)
# Figs. 6 and 7

# Adjusted parameters (p. 316)
sn <- 0.91
sn.sd <- 0.07

ss <- 0.77 # They estimate the parameter ss to 0.8 +/- 0.32, but conclude that it is not different from 0.77. But using 0.77 or 0.8 produces quite different output.
ss.sd <- 0.08
# Additional parameter
dn <- sn # Annual survivorship of floaters (p. 310)

# Maximum number of pairs (breeding females) in northern population
Tn <- 100
# Maximum number of pairs (breeding females) in southern population
Ts <- 120

# Matrix for years when Tn < 100
As <- matrix(c(
	0, yn*bn, 0, 0,
	rn*h, sn, rn*m*mc, 0,
	0, 0, 0, ys*bs,
	rs*m*mc, 0, rs*h, ss),
ncol=4,nrow=4,byrow=T)

# If Tn >= 100, the model changes to (p. 310):
Ad <- matrix(c(
	dn*h, yn*bn, dn*m*mc, 0,
	0, sn, 0, 0,
	0, 0, 0, ys*bs,
	rs*m, 0, rs*h, ss),
ncol=4,nrow=4,byrow=T)

nyears <- 200

N.sos <- array(dim=c(nyears+1, ncol(As), 1))

n0 <- c(19,60,20,29) # Fn = 19, Nn = 60, Fs = 20, Ns = 29)
N.sos[1,,] <- n0
n <- n0

# Assume no captive breeding, and local extinction of the northern population at t=100:
for (t in 1:nyears) {

	temp <- As %*% n

	if (temp[2] >= Tn) {
		n[2] <- Tn
		n <- (Ad %*% n) + matrix(c(-(1-sn)*n[2], (1-sn)*n[2], 0, 0), ncol=1, nrow=nrow(Ad))
	} else if (temp[2] < Tn) {
		n <- (As %*% n)
	}

	if (t == 100) n[1:2] <- c(0,0)

	N.sos[t+1,,] <- n
}

# Plot figure 6
dev.new(width=10, height=6)
op <- par(mar=c(5,5,4,2))
plot(x=c(0,nyears),y=c(0,100), type="n", xlab="Time (years)", ylab="No. of breeding pairs", font.lab=2, las=1, cex.axis=1.3, cex.lab=1.3,
	main="California peregrine source-sink population model",bty="l")

points(0:nyears, N.sos[,,1][,2],pch=16, col=1, type="p")
points(0:nyears, N.sos[,,1][,4],pch=15, col=2, type="p")

legend("bottomleft", c("North","South"), col=c(1,2), pch=c(16,15), bty="n", cex=1.3)
abline(v=100, lty=2)
abline(h=0,lty=2)

# Figure 7
# Vary migration rate
m <- c(0.27,0) # Migration = 0.27 or migration = 0

nyears <- 100
N.sos <- array(dim=c(nyears+1, ncol(As), length(m)))

n0 <- c(19,60,20,29) # Fn = 19, Nn = 60, Fs = 20, Ns = 29)
N.sos[1,,] <- n0
n <- n0


for (j in 1:length(m)) {

	n <- n0

	As <- matrix(c(
			0, yn*bn, 0, 0,
			rn*h, sn, rn*m[j]*mc, 0,
			0, 0, 0, ys*bs,
			rs*m[j]*mc, 0, rs*h, ss),
			ncol=4,nrow=4,byrow=T)

	# If Tn >= 100, the model changes to (p. 310):
	Ad <- matrix(c(
			dn*h, yn*bn, dn*m[j]*mc, 0,
			0, sn, 0, 0,
			0, 0, 0, ys*bs,
			rs*m[j], 0, rs*h, ss),
			ncol=4,nrow=4,byrow=T)

	# inner loop (over t)
	for (t in 1:nyears) {

		temp <- As %*% n

		if (temp[2] >= Tn) {
			n[2] <- Tn
			n <- (Ad %*% n) + matrix(c(-(1-sn)*n[2], (1-sn)*n[2], 0, 0), ncol=1, nrow=nrow(Ad))
		} else if (temp[2] < Tn) { n <- (As %*% n) }

		N.sos[t+1,,j] <- n

	} # end loop over t
} # End loop over j


# Plot figure 7A, Source-Sink dynamics
dev.new()
op <- par(mar=c(5,5,4,2))
plot(x=c(0,nyears),y=c(0,100), type="n", xlab="Time (years)", ylab="No. of breeding pairs", font.lab=2, las=1, cex.axis=1.3, cex.lab=1.3, main="California peregrine source-sink population model",bty="l")

points(0:nyears, N.sos[,,1][,2],pch=1, col=1, type="p")
points(0:nyears, N.sos[,,1][,4],pch=1, col=2, type="p")

points(0:nyears, N.sos[,,2][,2],pch=15, col=1, type="p")
points(0:nyears, N.sos[,,2][,4],pch=15, col=2, type="p")

legend(x=60, y=40, c("N","S", "N - No migration", "S - No migration"), col=c(1,2,1,2), pch=c(1,1,15,15), bty="n", cex=1.3)

# Stage-class composition
# Plot figure 8A and 8B
dev.new(width=8, height=6)
op <- par(mar=c(5,5,4,2))
plot(x=c(0,nyears),y=c(0,200), type="n", xlab="Time (years)", ylab="No. of breeding pairs", font.lab=2, las=1, cex.axis=1.3, cex.lab=1.3, main="California peregrine source-sink population model",bty="l")

polygon(c(0,0:nyears,nyears), c(0,rowSums(N.sos[,,1][,1:2]),0), density=25, angle=45)
polygon(c(0,0:nyears,nyears), c(0,N.sos[,,1][,2],0), col="darkgrey")

legend("topleft", legend=c("Non-breeders","Breeders"), fill=c("black","darkgrey"), density=c(25,100), angle=c(45,0), bty="n")
text(nyears,200,"A",font=2,cex=1.5)


dev.new(width=8,height=6)
op <- par(mar=c(5,5,4,2))
plot(x=c(0,nyears),y=c(0,200), type="n", xlab="Time (years)", ylab="No. of breeding pairs", font.lab=2, las=1, cex.axis=1.3, cex.lab=1.3, main="California peregrine source-sink population model",bty="l")

polygon(c(0,0:nyears,nyears), c(0,rowSums(N.sos[,,1][,3:4]),0), density=25, angle=45)
polygon(c(0,0:nyears,nyears), c(0,N.sos[,,1][,4],0), col="darkgrey")

legend("topleft", legend=c("Non-breeders","Breeders"), fill=c("black","darkgrey"), density=c(25,100), angle=c(45,0), bty="n")
text(nyears,200,"B",font=2,cex=1.5)
