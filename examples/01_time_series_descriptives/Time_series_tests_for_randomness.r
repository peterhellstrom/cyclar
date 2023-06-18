library(cyclar)

# Autocorrelated, AR(2)
z <- ar2.sim2(n=100, mu=2, a1=0, a2=-0.5)

# Random, uniform, series
# z <- runif(n=50)

plot(z)

turning.point.test(z)
difference.sign.test(z)
rank.test(z)

table(diff(z)>0)

library(pastecs)
ztp <- turnpoints(z)
ztp
summary(ztp)
extract.turnpoints(ztp)
plot(z)
lines(ztp,median=FALSE)
