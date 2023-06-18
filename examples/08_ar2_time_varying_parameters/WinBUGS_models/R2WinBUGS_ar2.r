################################################################################
# Load required libraries for WinBUGS
library(R2WinBUGS)
library(coda)
library(emdbook)
################################################################################

inits <- list(
list(a1 = 0, a2 = -0.8, tau.data = 5),
list(a1 = 0.05, a2 = -0.7, tau.data = 4),
list(a1 = -0.05, a2 = -0.65, tau.data = 4)
)

b.bugs <- bugs(
      data=list(growth=growth, lag1=lag1, lag2=lag2, ntime=ntime),
      inits,
      parameters.to.save = c("a1","a2","sigma.data"),
      n.thin = 1,
      n.burnin = 2000,
      model.file = "C:/WORK/ANALYSER/-= Population dynamics (time series) =-/Autoregressive models/Time-varying parameters/WinBUGS_models/ar2.bug",
      n.chains = length(inits),
      n.iter = 20000,
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
