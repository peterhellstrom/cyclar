################################################################################
# Load required libraries for WinBUGS
library(R2WinBUGS)
library(coda)
library(emdbook)
################################################################################

inits <- list(
list(
a1 = a1.t[3:(ntime+2)], a2 = a2.t[3:(ntime+2)], 
#rep01 = rep(0.01,ntime-2), rep02 = rep(0.01,ntime-2),
tau.data = 5, tau.smooth1 = 5))

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


################################################################################

model0.coda <- read.bugs(b.bugs)
bugs.out(model0.coda)

# Some plotting...

plot(model0.coda)


densityplot(model0.coda)


densityplot(model0.coda)[2]


densityplot(model0.coda)[3]
