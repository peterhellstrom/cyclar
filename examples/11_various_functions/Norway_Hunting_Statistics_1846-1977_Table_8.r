library(cyclar)

# 14 species
# 18 counties, 131 years
# Columns 1 & 2; county & year
# Should be 18*131 = 2358 rows

dat <- read.delim("Table 8.txt", header=T, na.strings=c("NA","-"))
str(dat)
dim(dat)

colnames(dat)
colnames(dat)[3:ncol(dat)]


dat.species <- lapply(3:ncol(dat), function(i) data.frame(County=dat$County, Year=dat$Year, Bounties=dat[,i]))
names(dat.species) <- colnames(dat)[3:ncol(dat)]
str(dat.species)


ts.species <- lapply(1:length(dat.species), function (i) {
		z <- dat.species[[i]]
		z.ts <- ts(tapply(z[,3],z[,c(2,1)],c),start=1846)
		z.ts
	}
)
names(ts.species) <- names(dat.species)

dput(ts.species,"Table 8.r")

# BASIC DESCRIPTIVE STATS

bounty.summary(x=bountyNorway, species="Bear")
bounty.summary(x=bountyNorway, species="Wolf") # Returns an error, because final "year" of payment is 1866-1870 (which can not be converted to a numeric value)
bounty.summary(x=bountyNorway, species="Wolverine")
bounty.summary(x=bountyNorway, species="Lynx")
bounty.summary(x=bountyNorway, species="Otter")
bounty.summary(x=bountyNorway, species="Fox")
bounty.summary(x=bountyNorway, species="Marten")
bounty.summary(x=bountyNorway, species="Mink")
bounty.summary(x=bountyNorway, species="Eagle")
bounty.summary(x=bountyNorway, species="Eagle.owl")
bounty.summary(x=bountyNorway, species="Falcon")
bounty.summary(x=bountyNorway, species="Goshawk")
bounty.summary(x=bountyNorway, species="Sparrowhawk")
bounty.summary(x=bountyNorway, species="Loon")

# Wolverine and mink are the only species with no payments in at least one county

# Alternate version of plot
bounty.summary(x=bountyNorway, species="Fox", plot.type="single", scale=TRUE)

# Plot Elton's lynx data, do not use the species argument
bounty.summary(x=lynx.elton, plot.type="single", scale=TRUE)

# Example on how to subset the fox data and plot only one time series:

# Fox example
x <- bountyNorway$Fox[,"Nordland"]
plot(window(x,1879,1976), xlab="",ylab="", main="Fox in Nordland 1879-1976",font.lab=2, type="l", lwd=1, lty=1, col=2, bty="l", xaxt="n")
axis(1,at=seq(1870,1980,10))

# Goshawk example
x <- bountyNorway$Goshawk[,"Troms"]
plot(window(x,1846,1932),xlab="",ylab="", main="Goshawk in Troms 1846-1932",font.lab=2, type="l", lwd=1, lty=1, col=2, bty="l", xaxt="n")
axis(1,at=seq(1840,1940,10))
