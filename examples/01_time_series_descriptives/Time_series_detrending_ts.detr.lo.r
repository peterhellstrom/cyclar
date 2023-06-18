library(cyclar)

n <- 30
x <- arima.sim(n=n, model=list(ar=c(0,-.71)))

ts.detr.lo(x)

ts.detr.lo(x, degree=1)

fm1 <- ts.detr.lo(x, degree=1)
fm1[,"x.loess"] - fm1[,"x.lowess"]

ts.detr.lo(x, degree=1, family="gaussian")

ts.detr.lo(x, degree=1, span=0.2)

# How to obtain similar results for loess and lowess?
# 1) family = "symmetric"
# 2) degree = 1
