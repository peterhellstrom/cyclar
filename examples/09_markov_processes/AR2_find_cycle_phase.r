################################################################################
library(cyclar)
################################################################################
# Phase classification:
# 1 = Peak, 2 = Decline, 3 = Low, 4 = Increase

a0 <- 2.1 # mean
p <- ar2.parms(k=4, v=2)
ar.var <- 0.2 # variance
n <- 50 # number of years

x <- ar2.sim(n=n, mu=a0, a1=p$a1, a2=p$a2, sd=sqrt(ar.var))
ts.summary(x, plot=TRUE)
# x <- exp(x)

# write.table(x, "Markov_test_file_to_Excel.txt", sep="\t", quote=F)

# Check computation time:
system.time(out1 <- phase(x, method="for", plot=TRUE))
system.time(out2 <- phase(x, method="logic", plot=TRUE))
system.time(out3 <- phase(x, method="quantile", plot=TRUE))
system.time(out3 <- phase(x, method="Henden", plot=TRUE))

# Turning plots off (plot=FALSE) reduces computation time
system.time(out1 <- phase(x, method="for", plot=FALSE))
system.time(out2 <- phase(x, method="logic", plot=FALSE)) # Fastest
system.time(out3 <- phase(x, method="quantile", plot=FALSE))

# On my laptop, logic-method took 0.24 seconds for a time series of length 5000
# Repeat this 10 000 times ~ approx time 40 mins.
# Repeat shorter series (n=100) 10 000 times ~ 3.5 mins

# Check correct alignment
# ts(data.frame(x=out2$x, xd=out2$xd, z=out2$z, ztr=out2$ztr))

# Note on system time:
# Placing and declaring the function row.matches outside of find.phase speeds up,
# as does switching the plot option to FALSE. The option "match" is also considerably
# quicker than the for-loop )option "for").

table(out1$z==out2$z) # Methods for and logic should yield identical results!


# Check correlation between state and density:
plot(as.factor(out1$z), out1$x)
fm1 <- lm(out1$x ~ as.factor(out1$z)-1)
summary(fm1)
anova(fm1)

par(mfrow=c(2,2))
plot(fm1)
par(mfrow=c(1,1))

out1.narm <- na.omit(cbind(x=out1$x, z=as.factor(out1$z)))
str(out1.narm)

library(nlme) # gls, check if model with variance per phase fits better
fm1.gls <- gls(x ~ as.factor(z), data=out1.narm)
fm2.gls <- gls(x ~ as.factor(z), weights=varIdent(form = ~1|as.factor(z)), data=out1.narm)
anova(fm1.gls, fm2.gls)

# Check average cycle length
state.cycle(x=out2, phase=1)
barplot(table(state.cycle(x=out2, phase=1)$period))
sapply(1:4, function(i) state.cycle(x=out2, phase=i)$mean.period) # Differs for phases 1:4

# Fox litters Helags
litters <- ts(c(1,3,6,0,1,7,1,0,5,12,0,1,4,13,1,0,5,3,4,1,1,4,4,1,1,2,4,1,4,6,1,9,14,0), start=1976)
plot(litters,main="Arctic fox litters, Helags")

hel <- phase(litters)
phase.tr(litters)
state.cycle(hel, phase=3)

hel1 <- phase(litters, method="for", plot=TRUE)
hel2 <- phase(litters, method="logic", plot=TRUE)
table(hel1$z==hel2$z)

# Compare with a random-walk:
n <- 100
rts <- ts(runif(n,0,4))
mean(rts)
phase(rts)
state.cycle(phase(rts))
