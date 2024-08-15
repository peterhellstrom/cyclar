library(cyclar)

# Two functions that generates equally spaced parameter values for (1+a1) & a2.
# ar2.k.gen = for k-contours (quasi-period)
# ar2.v.gen = for v-contours (intrinsic process variance)

# ar2.k.gen has two method options, method=c("x","y")
# ar2.v.gen has one method option, method=c("y")

# The character 'x' or 'y' specifies along which axis that the values should be generated.
# Use either length.out or by as extra input arguments to ar2.k.gen/ar2.v.gen
# these arguments are necessary for the internal generation of sequences with seq().

# Test ar2.k.gen

n <- 20

# For method="y", use length.out as argument in k.gen (if you want exactly the same spacing between y [a2] for different values of k)
# Method y
k <- 7
p <- ar2.k.gen(k=k,method="y",length.out=n)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
diff(p$phi2) # constant

k <- 5
p <- ar2.k.gen(k=k,method="y",length.out=n)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
diff(p$phi2) # constant

ar2.period(p$a1,p$a2)
plot(p$a2, ar2.period(p$a1,p$a2)[,4])


k <- 5
p <- ar2.k.gen(k=k,method="y",by=0.25)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
diff(p$phi2) # constant

k <- 7
p <- ar2.k.gen(k=k,method="y",by=0.25)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
diff(p$phi2) # constant

# For method="x", use by as argument in k.gen (if you want exactly the same spacing between x (1+a1) for different values of k).
# Method x
k <- 7
p <- ar2.k.gen(k=k,method="x",by=0.05)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
length(p$phi1)
diff(p$phi1) # constant
diff(p$phi2) # increases with x (1+a1)

k <- 5
p <- ar2.k.gen(k=k,method="x",by=0.05)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
length(p$phi1)
diff(p$phi1)

k <- 5
p <- ar2.k.gen(k=k,method="x",length.out=20)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
length(p$phi1)
diff(p$phi1)

k <- 5
p <- ar2.k.gen(k=k,method="x",length.out=20,from=0,to=0.5)
ar2.plot(k=k,v=NULL)
points(p$phi1,p$phi2,col=2,pch=16)
length(p$phi1)
diff(p$phi1)

# Test ar2.v.gen ----

n <- 20
v <- 4

p <- ar2.v.gen(v=4,length.out=n)
ar2.plot(v=v)
points(p$phi1,p$phi2,col=2,pch=16)
length(p$phi1)
diff(p$phi1)
diff(p$phi2) # constant

n <- 50
v <- 4
p <- ar2.v.gen(v=4,vregion="triangle",length.out=n)
ar2.plot(v=v,vregion="triangle")
points(p$phi1,p$phi2,col=2,pch=16)
length(p$phi1)
diff(p$phi1)
diff(p$phi2) # constant
