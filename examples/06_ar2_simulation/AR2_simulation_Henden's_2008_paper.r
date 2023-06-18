################################################################################
library(cyclar)
################################################################################
# Function that simulates an AR2-process is available in the cyclar-"package".
# ar2.sim

# Simulate rodent series, with varying period, mean, variance
# Parameters in John-André's first paper, J Appl Ecol 2008
p <- ar2.parms(k=3:5, v=rep(2,3), output=list) # get the correct parameters
ar2.plot(k=3:5, v=2) # plot in a1/a2 space
points(1+p$a1, p$a2, pch=16, col=2) # add points used in simulations

mu <- c(2.1, 2.3, 2.5) # Mean
a1 <- p$a1 # direct d-d, Period 3, 4, 5
a2 <- p$a2 # delayed d-d, Period 3, 4, 5
ar.var <- c(0.1,0.2,0.3)
period <- c(3,4,5)
n <- 10000

# Combine into a "long" table for later simulations:
# Should be 27 rows
design <- expand.grid(Period=period, Mean=mu, Variance=ar.var)

# Assign values of dd parameters
p <- ar2.parms(k=design$Period, v=rep(2,nrow(design)), list=FALSE)
design <- cbind(design, p[,1:2])
# Sort data frame. Feed this data frame to the function ar2
design <- design[order(design$Period, design$Mean),]
rownames(design) <- 1:nrow(design)

# simulate 27 time series:
system.time(
ts.sim <- lapply(1:nrow(design), function(i) ar2.sim(n=n, mu=design$Mean[i], a1=design$a1[i], a2=design$a2[i], sd=sqrt(design$Variance)[i])))

#str(ts.sim)

p.est <- t(sapply(1:length(ts.sim), function(i) coef(arima(ts.sim[[i]], order=c(2,0,0)))))
ar2.plot(k=3:5, v=2)
points(p.est[,"ar1"], p.est[,"ar2"])

# Create labels (i.e. variable names) for each time-series
# The labels are created as follow (using the paste function):
# Each label starts with ts.sim (=simulated time series)
# After that follows a string of four numbers.
# The first gives the period (3,4,5)
# The second number shows the Mean (1,3,5). Since all simulated series starts with 2, that information is not necessary to include
# The third is a separator (0)
# The second gives the variance (0.1,0.2,0.3). Only 1, 2 or 3 is shown, not the leading zero.
labs <- paste("ts.sim", design$Period, substr(design$Mean,3,3),  paste(0, substr(design$Variance,3,3), sep=""), sep="")
names(ts.sim) <- labs

# Extract statistics from each time series (median, and first and third quantiles)
# Could also use ts.summary()
extr.stat <- t(sapply(1:nrow(design), function(i) c(
	Median = median(ts.sim[[i]]),
	quantile(ts.sim[[i]], c(0.25,0.75))
)))

# Make a summary table of the simulations, by combining the tables design and extr.stat
ts.sim.tab <- data.frame(design, extr.stat)
# rownames(ts.sim.tab) <- labs

# The ts.sim.tab table can be re-ordered if you wish
ts.sim.tab[order(ts.sim.tab$Period, ts.sim.tab$Variance),]

################################################################################
# Write output to external files (txt-files and R-object for future use)
write.table(ts.sim.tab, "AR2_simulations_table.txt", row.names=F, col.names=T, sep="\t", quote=F)

# Output all the time series to a dat-file (can be opened in Notepad, Excel etc.). Or read back into R.

# Since the time series are stored in a list, convert them to a matrix (and they are no longer time series objects!)
# Additional attributes stored in the time series in list form are also lost in the process!
z1 <- sapply(ts.sim, cbind)
write.table(z1, "Simulated_AR(2)_data.dat", row.names=F, col.names=T, sep="\t", quote=F)
z1.back <- read.table("Simulated_AR(2)_data.dat", header=T) # check
#str(z1.back)

# Print output to an Excel-sheet
# Requires
library(xlsReadWrite)
write.xls(z1,"Simulated_AR(2)_data.xls")

# Store all simulated time-series as an R-object for future use.
dput(ts.sim,"ts_sim")
z2 <- dget("ts_sim") # check
# str(z2)

################################################################################
