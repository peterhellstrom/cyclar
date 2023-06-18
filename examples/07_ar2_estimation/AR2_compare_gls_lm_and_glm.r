library(cyclar)

p <- ar2.parms(k=3,v=4)
a1 <- p$a1
a2 <- p$a2

test <- ar2.sim(n=48, mu=2, a1=a1, a2=a2, sd=sqrt(0.3))
ts.summary(test)

ar.univar(test, log=FALSE)

# Compare gls with glm fits
fm.gls <- ar2.lm(lags(test,2), model="gls")
summary(fm.gls$models$fm12)

fm.lm <- ar2.lm(lags(test,2), model="lm")
summary(fm.lm$models$fm12)

fm.glm <- ar2.lm(lags(test,2), model="glm")
summary(fm.glm$models$fm12)

ar2.pred.plot(fm.gls, sp=TRUE)
ar2.pred.plot(fm.lm, sp=TRUE)
ar2.pred.plot(fm.glm, sp=TRUE)


fm.gls$coef
fm.lm$coef

confint(fm.gls$models$fm12)
confint(fm.lm$models$fm12)

ar2.plot(ylim=c(-1,0), triangle=FALSE)
points(1+a1,a2,pch=16,col=2)
points(fm.gls$coef[4,2],fm.gls$coef[4,3],col=3,pch=15)
points(fm.lm$coef[4,2],fm.glm$coef[4,3],col=4,pch=15)

fm.lm.R <- ar2.lm(R,model="lm",method="growth rate")
summary(fm.lm$models$fm12)
ar2.pred.plot(fm.lm.R, sp=TRUE) # Doesn't work yet...
