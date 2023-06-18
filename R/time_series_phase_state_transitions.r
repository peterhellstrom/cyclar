# phase
# Define cycle phases (four states/phases):
# Two methods:
# qualitative states (following Haydon et al 2003) - arguments method="for" and method="logic"
# quantitative states (based on three quantile limits)

# Note on how to read the transition matrix for a Markov chain:
# Column is start, row is end. Columns sum to one.
# Note that e.g. Haydon et al uses a transposed matrix, where rows is start and columns is end.
# Both ways are of course equal, just remember which way you are going...
# I used the column-end transpose, since this is similar to Leslie matrices.

# x is a time series
# The available methods [method = c("for","logic","quantile")] are:
# for (uses a for-loop)
# logic
# quantile

# Dependencies: state.logic()
#' @export
phase <- function(x, method=c("for","logic","quantile","Henden"), qs=c(0.25, 0.5, 0.75), plot=TRUE) {

	if (class(x) != "ts") x <- ts(x)
	if (class(x) == "mts") stop("Multiple time series are not implemented")

	# NAs are removed - not a very good option if NAs not are at the ends of a time series
	# (but only method currently implemented)
	inds.na <- which(is.na(x))
	if (length(inds.na) > 0) {
		message(paste("The following observations are NA and will be removed:", inds.na))
		x <- ts(na.omit(x)) # Omit NAs
	}

	# Calculate changes in ts with diff, then take sign to get a series with -1,0 and 1.
	x.diff <- sign(diff(x))
	n <- length(x)

	method <- match.arg(method)
	# What about first and last position?
	# For last position: It should be possible to determine if increase/decrease, but not turning point
	if (method == "for") {
		z <- numeric(length(x))
		z[c(1,length(z))] <- NA
		# Classification of phases (states) follows Haydon et al 2003 (Proc R Soc B 270:435-445)
		for (i in 2:(length(x)-1)) {
			if(x[i-1] > x[i] & x[i] <= x[i+1]) z[i] <- 1 # Low (through)
			if(x[i-1] <= x[i] & x[i] <= x[i+1]) z[i] <- 2 # Increase (intermediate)
			if(x[i-1] <= x[i] & x[i] > x[i+1]) z[i] <- 3 # Peak
			if(x[i-1] > x[i] & x[i] > x[i+1]) z[i] <- 4 # Decrease (intermediate)
		}
	} else if (method=="logic") {
		n.diff <- length(x.diff)
		xv.diff <- cbind(x.diff[1:(n.diff-1)], x.diff[2:n.diff])
		z <- state.logic(xv.diff) # Call to function state.logic()
		z <- c(NA,as.numeric(z),NA)
	} else if (method=="quantile") {
		if (length(qs) != 3) stop("Give exactly three quantile ranges")
		qs <- quantile(x, qs)
		z <- numeric(length(x))
		z[which(x < qs[1])] <- 1 # Low
		z[which(x >= qs[1] & x < qs[2])] <- 2 # Intermediate 1
		z[which(x >= qs[2] & x < qs[3])] <- 3 # Intermediate 2
		z[which(x >= qs[3])] <- 4 # High
	} else if (method == "Henden") {
		z <- numeric(length(x))
		z[1:2] <- NA
		for (i in 3:(length(x))) {
			if(x[i] < x[i-1] & x[i-1] < x[i-2]) z[i] <- 1 # Low (through)
			if(x[i] > x[i-1] & x[i-1] < x[i-2]) z[i] <- 2 # Increase (intermediate)
			if(x[i] > x[i-1] & x[i-1] > x[i-2]) z[i] <- 3 # Peak
			if(x[i] < x[i-1] & x[i-1] > x[i-2]) z[i] <- 4 # Decrease (intermediate)
		}
	}

	# Estimate phat (probability of increase between two time points)
	phat <- table(x.diff)/(n-1)
	phat <- as.numeric(phat["1"])

	zres <- phase.tr(z)

	if (plot==TRUE) {

    # Time vector, x-axis
    xv <- seq(1:length(x))
    # Set up plot layout
    nf <- layout(matrix(c(3,4,3,1,3,2,3,5), nrow=4, ncol=2, byrow=TRUE),
			widths=c(1,10,10,10), heights=c(1.5,8,8,1.5))

    op <- par(mar=c(2,2,2,1))

		# Create first plot
    	if (method != "quantile") {
			plot(xv,x,type="l",xlab="",ylab="",font.lab=2,las=1,main="",cex.axis=1.4)
				for (i in 1:length(unique(z))) points(xv[z==i],x[z==i],col=i,pch=16,cex=1.4)
			par(xpd=NA)
			tmp.u <- par('usr')
			legend(tmp.u[1], tmp.u[4], xjust=-4, yjust=0,
				paste(1:4, c("low","increase","peak","decrease")),col=1:4,cex=1.4,pch=16,bty="n")
		}
		if (method == "quantile") {
			plot(xv,x,type="l",xlab="",ylab="",font.lab=2,las=1,main="",cex.axis=1.4)
				for (i in 1:4) points(xv[z==i],x[z==i],col=i,pch=16,cex=1.4)
			abline(h=qs[1],lty=2,col=1)
			abline(h=qs[2],lty=2,col=2)
			abline(h=qs[3],lty=2,col=4)
			par(xpd=NA)
			tmp.u <- par('usr')
			legend(tmp.u[1], tmp.u[4], xjust=-3, yjust=0,
				paste(1:4, c("low","intermediate 1","intermediate 2","high")),col=1:4,cex=1.4,pch=16,bty="n")
		}

		# Create second plot
		par(xpd=FALSE)
		plot(xv,x,type="l",xlab="",ylab="",font.lab=2,las=1,main="",cex.axis=1.4)
		points(xv[turnpoints(x)$tppos],x[turnpoints(x)$tppos],col=1,pch=16,cex=1.4)
		abline(h=quantile(x, c(0.25, 0.5, 0.75)),lty=3,col=2)
		plot.new()
			mtext("Population density",line=-1,side=2,font=2,cex=1.3)
		plot.new()
			mtext("Phase classification",line=-1.5,font=2,cex=1.3)
		plot.new()
			mtext("Time",font=2,cex=1.3)
		par(op)
	}

	out <- list(
			n = n,
			states = n-2,
			transitions = n-3,
			x = x, # Original time series
			z = ts(z), # Vector time series with states (1,2,3,4),
			T = zres$Tf, # Frequency table of state transitions
			Tf = zres$T, # Proportion of state transitions
			phat = phat # Proportion of increases between two time points
	)
	out
}

# phase.tr
# Function that tabulates state transitions, given a series z
# containing a vector with state transitions.
# Can be used independently, but intended use is to call it from
# the main function phase().
# Can also be used together with markov.sim, to estimate state transitions
# for a simulated Markov chain
# Inputs: z  = sequence of environment states
#' @export
phase.tr <- function(z) {

	# leaders
	z1 <- paste(z[-length(z)])
	# followers
	z2 <- paste(z[-1])

	z.list <- list(followers=z2,leaders=z1)
	Tf <- table(z.list)
	# prop.table converts a count table to proportions, either by rows
	# or by columns, depending on 2nd argument (margin) - see its help file
	# if margin=2, prop.table assumes that P matrix is column-wise!
	T <- prop.table(Tf,2)
	ll <- sum(Tf[Tf>0] * log(T[T>0]))

	list(
		Tf = Tf,
		T  = T,
		logLik = ll)
}

# state.cycle
# Function to calculate number of cycles in a sequence,
# as well as number of steps between a given phase (default 1, i.e. through),
# average cycle length and frequence table of cycle lengths:
# x is a time series fitted with phase!

# Previous name: mean.cycle
#' @export
state.cycle <- function(x, phase=1) {

	tempx <- x$x
	tempz <- x$z

	inds <- which(tempz == phase)
	n <- length(inds)
	period <- diff(inds)
	n.period <- length(period)
	mean.period <- mean(period)
	tab <- table(period)

	out <- list(
			x = tempx,
			phase = tempz,
			phase = phase,
			n.phase = n,
			period = period,
			n.period = n.period,
			mean.period = mean.period,
			sd.period = sd(period),
			table = tab,
			prop.table = tab/sum(tab)
			)
	out
}

# state.logic
# A function that logically evaluates the state of a time series,
# given a sign-differentiated input matrix with two columns. This function is called by the main function phase.
#' @export
state.logic <- function(x) {
	y <- numeric(nrow(x))
	y[which(x[,1]==-1 & x[,2]==1 | x[,1]==-1 & x[,2]==0)] <- 1 # 1 = Low
	y[which(x[,1]==1 & x[,2]==1 | x[,1]==0 & x[,2]==1 | x[,1]==1 & x[,2]==0 | x[,1]==0 & x[,2]==0)] <- 2 # 2 = Increase
	y[which(x[,1]==1 & x[,2]==-1 | x[,1]==0 & x[,2]==-1)] <- 3 # 3 = Peak
	y[which(x[,1]==-1 & x[,2]==-1)] <- 4 # 4 = Decrease
	y
}
