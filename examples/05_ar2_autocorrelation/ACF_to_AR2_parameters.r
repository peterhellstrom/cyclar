library(cyclar)

# Check if the "correct" input simulation parameters
# can be retrieved from the acf (by using the acf2ar-function):

# Set up simulation
nrep <- 1000 # Number of simulations
ntime <- 300 # Length of each simulation
out <- matrix(nrow=nrep, ncol=6) # Output matrix
colnames(out) <- c("a1","a2","acf.a1.N","acf.a2.N","acf.a1.R","acf.a2.R")

p <- ar2.parms(k=5,v=2,output=list) # Change values of k and v as desired
ar.var <- 0.2
ar2.period(p$a1,p$a2)

for (i in 1:nrep) {

	x <- ar2.sim(n=ntime, mu=0, a1=p$a1, a2=p$a2, sd=sqrt(ar.var))
	gr.x <- diff(x) # growth rate
	
	acf.x <- acf(na.omit(x),plot=FALSE)$acf
	acf.gr.x <- acf(na.omit(gr.x),plot=FALSE)$acf

	out[i,1:2] <- c(p$a1, p$a2)
	out[i,3:4] <- acf2AR(acf.x)[2,1:2]
	out[i,5:6] <- acf2AR(acf.gr.x)[2,1:2]
}

ar2.plot()
legend("topright",legend=c("Simulated","acf.x","acf.R"),pch=c(16,1,2),col=c(1,2,3),bty="n")

points(out[,3],out[,4],col=2,cex=0.75,pch=1) #x[t]
points(out[,5],out[,6],col=3,cex=0.75,pch=2) #r[t]
points(1+out[,1],out[,2],cex=1.5,pch=16,col=1) # input parameters
