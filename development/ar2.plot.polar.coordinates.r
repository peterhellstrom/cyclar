library(cyclar)

# So far this only works for k > 4.
# If k < 4, the equation for the straight line must be re-calculated.

k <- 6
b <- 0.75 # Slope of line starting at (0,-1)
b <- ifelse(k < 4, -b, b)

ar2.plot(k=k, v=NULL)
abline(-1,b,col=4)

# Find intersection
z <- 4*cos(2*pi/k)^2
x.star <- ifelse(k < 4, -(sqrt(b^2 * z^2 + 4*z) + b*z)/2, (sqrt(b^2 * z^2 + 4*z) - b*z)/2)
y.star <- b*x.star - 1
points(x.star,y.star,col=2,pch=16)

# Convert from cartesian to polar coordinates
# tan(theta) = y/x ==> theta = arctan(y/x)
r <- sqrt(x.star^2 + y.star^2)
theta <- atan(y.star / x.star)
# Convert from radian to degress
theta * (180/pi)

r; theta
y.star; r*sin(theta) # negative sign i f k < 4
x.star; r*cos(theta) # negative sign i f k < 4
y.star/x.star; tan(theta)

# should equal b
(r*sin(theta) + 1) / (r*cos(theta)); b # if k > 4
(-r*sin(theta) + 1) / (-r*cos(theta)); b # if k < 4

# Variance
v <- ar2.ipv(x.star-1,y.star)
v

ar2.v.add(v=v,col=2)
(1-r*sin(theta)) / (1+r*sin(theta))

# Find tangent for k-contour at x.star
# http://en.wikipedia.org/wiki/Tangent

slope <- -x.star/(2*cos(2*pi/k)^2)
y.at.x <- ar2.k(a1=x.star-1, k=k)
a.tg <- y.at.x - slope*x.star
b.tg <- slope
abline(a.tg,b.tg,lty=2)

a.tg; b.tg

asin(y.star/r); theta
acos(x.star/r); theta

k <- 6
v <- 1.25
p <- ar2.parms(k=k,v=v)
p

r <- as.numeric(sqrt((1+p$a1)^2 + p$a2^2))
theta <- as.numeric(atan(p$a2/(1+p$a1)))

r; theta
theta*180/pi

b <- as.numeric((1 + p$a2) / (1 + p$a1))
b

ar2.plot(k=k, v=v)
points(1+p$a1, p$a2, pch=16, col=2)
abline(-1,b,col=4)
