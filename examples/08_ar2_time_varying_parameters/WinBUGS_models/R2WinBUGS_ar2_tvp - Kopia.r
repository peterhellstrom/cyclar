# Load required libraries for WinBUGS
library(R2WinBUGS)
library(coda)
library(emdbook)
library(cyclar)
library(xlsReadWrite)

# Read Excel sheet
dat <- read.xls(normalizePath("FoxData1961-2008.xls"), colNames = TRUE, sheet=1, colClasses=rep("numeric",10))
names(dat)

# Convert to time series:
fox <- vector("list",ncol(dat)-1)
for (i in 2:ncol(dat)) fox[[i-1]] <- ts(dat[,i], start=dat$year[1], frequency=1)
names(fox) <- paste(c(rep("rr",6), rep("fr",3)),
	c(c("z","ac","bd","w","x","y"),c("z","ac","bd")), sep=".")


i <- 2
dat <- growth.rate(fox[[i]],lags=1,log=TRUE)
dat <- na.omit(dat)
ntime <- nrow(dat)

growth <- as.numeric(dat[,"R"])
lag1 <- as.numeric(dat[,"Lag0"])
lag2 <- as.numeric(dat[,"Lag1"])

inits <- list(
list(a1 = runif(ntime,-1,1), a2 = runif(ntime,-1,0), tau.data = 5, tau.smooth1 = 5),
list(a1 = runif(ntime,-1,1), a2 = runif(ntime,-1,0), tau.data = 4, tau.smooth1 = 4))

b.bugs <- bugs(
      data=list(growth=growth, lag1=lag1, lag2=lag2, ntime=ntime),
      inits,
      parameters.to.save = c("a1","a2","sigma.data","diffa1","diffa2","rep01","rep02"),
      n.thin = 1,
      n.burnin = 20000,
      model.file = "C:/WORK/ANALYSER/-= Population dynamics (time series) =-/Autoregressive models/Time-varying parameters/WinBUGS_models/ar2_tvp.bug",
      n.chains = length(inits),
      n.iter = 100000,
      codaPkg = TRUE,
      debug = TRUE
)

model0.coda <- read.bugs(b.bugs)
bugs.out(model0.coda)

# Some plotting...

plot(model0.coda)


densityplot(model0.coda)


densityplot(model0.coda)[2]


densityplot(model0.coda)[3]
