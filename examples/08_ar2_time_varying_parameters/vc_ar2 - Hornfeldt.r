library(cyclar)

cr.vindeln <- plot.birger(study.area="vindeln", study.species="Cr", output=TRUE)
cg.vindeln <- plot.birger(study.area="vindeln", study.species="Cg", output=TRUE)
ma.vindeln <- plot.birger(study.area="vindeln", study.species="Ma", output=TRUE)


gr.cr.spring <- growth.rate(cr.vindeln$Spring, constant=1, lags=2)
gr.cr.fall <- growth.rate(cr.vindeln$Fall, constant=1, lags=2)

gr.cg.spring <- growth.rate(cg.vindeln$Spring, constant=1, lags=2)
gr.cg.fall <- growth.rate(cg.vindeln$Fall, constant=1, lags=2)

gr.ma.spring <- growth.rate(ma.vindeln$Spring, constant=1, lags=2)
gr.ma.fall <- growth.rate(na.vindeln$Fall, constant=1, lags=2)

gr.vindeln <- list(
	gr.cr.spring = gr.cr.spring,
	gr.cr.fall = gr.cr.fall,
	gr.cg.spring = gr.cg.spring,
	gr.cg.fall = gr.cg.fall,
	gr.ma.spring = gr.ma.spring,
	gr.na.fall = gr.ma.fall
	)

npatch <- length(gr.vindeln)

z <- sapply(1:npatch, function(i) gr.vindeln[[i]][,"R"])

plot(x=c(1971,2009),y=range(z,na.rm=TRUE), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Observed time series")
for (i in 1:length(z)) points(z[[i]], type="l", col=i)


# Plot changes in variability
z.s <- sapply(1:npatch, function(i) swin.var(x=z[[i]], swin=6))
plot(x=c(0,39),y=c(-1,2), type="n",xlab="Time (center of time-window)",ylab="S-index",main="Changes in crude variability over time", font.lab=2, las=1, bty="l")
for (i in 1:npatch) points(z.s[[i]], type="l", col=i)

# Create input files (batch)
inp.list.x <- lapply(1:npatch, function(i) vc.inp(x=gr.vindeln[[i]], method="x", input="growth.rate", file = paste("vc_ar2_x_",i,".csv", sep="")))

inp.list.r <- lapply(1:npatch, function(i) vc.inp(x=gr.vindeln[[i]], method="r", input="growth.rate", file = paste("vc_ar2_r_",i,".csv", sep="")))

# RUN ESTIMATION:
# "C:/Program Files/Varying Coefficients Estimation/VC.exe"

# Check output
out.list.x <- lapply(1:npatch, function(i) vc.out(vc.dat=inp.list.r[[i]], method="x", input="data", file = paste("vc_ar2_x_",i,"-VC.csv", sep="")))

out.list.r <- lapply(1:npatch, function(i) vc.out(vc.dat=inp.list.r[[i]], method="r", input="data", file = paste("vc_ar2_r_",i,"-VC.csv", sep="")))


